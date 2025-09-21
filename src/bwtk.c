/*
 *   bwtk: bigWig Toolkit
 *   Copyright (C) 2025  Benjamin Jean-Marie Tremblay
 *
 *   This program is free software: you can redistribute it and/or modify
 *   it under the terms of the GNU General Public License as published by
 *   the Free Software Foundation, either version 3 of the License, or
 *   (at your option) any later version.
 *
 *   This program is distributed in the hope that it will be useful,
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *   GNU General Public License for more details.
 *
 *   You should have received a copy of the GNU General Public License
 *   along with this program.  If not, see <https://www.gnu.org/licenses/>.
 *
 */

#include <stdlib.h>
#include <stdint.h>
#include <stdio.h>
#include <string.h>
#include <stdbool.h>
#include <errno.h>
#include <math.h>
#include "libBigWig/bigWig.h"
#include "ksort.h"
#include "kseq.h"
#include "khash.h"

#define BWTK_VERSION "1.0"
#define BWTK_YEAR "2025"

// common ----------------------------------------------------------------------

KSTREAM_INIT(gzFile, gzread, 8192)
KHASH_MAP_INIT_STR(str_m, int64_t)

#include <sys/resource.h>
static long peak_mem(void) {
  struct rusage r_mem;
  getrusage(RUSAGE_SELF, &r_mem);
#ifdef __linux__
  return r_mem.ru_maxrss * 1024;
#else
  return r_mem.ru_maxrss;
#endif
}

static void *calloc_or_die(size_t size, const char *func_name) {
  void *result = calloc(1, size);
  if (result == NULL) {
    fprintf(stderr, "[E::%s] Out of memory (requested %lu B).\n",
      func_name, size);
    exit(EXIT_FAILURE);
  }
  return result;
}
#define alloc(size) calloc_or_die((size), __func__)

#ifndef kroundup32
#define kroundup32(x) (--(x), (x)|=(x)>>1, (x)|=(x)>>2, (x)|=(x)>>4, (x)|=(x)>>8, (x)|=(x)>>16, ++(x))
#endif

// Some BED handling code adapted from HTSlib/bedidx.c

#define LIDX_SHIFT 13
#define ALL 0
#define FILTERED 1

typedef struct {
    uint32_t beg, end;
    char strand;
    int64_t namei;
} bedItem_t;

static inline int lt_bedItem(bedItem_t a, bedItem_t b) {
    if (a.beg == b.beg) return a.end < b.end;
    return a.beg < b.beg;
}
KSORT_INIT_STATIC(bedItem_t, bedItem_t, lt_bedItem)

typedef struct {
    int n, m;
    bedItem_t *a;
    int *idx;
    int filter;
} bedList_t;
KHASH_MAP_INIT_STR(bedHash, bedList_t)

typedef struct {
    int64_t n;
    uint32_t size;
    int n_names, m_names;
    char **names;
    kh_bedHash_t *data;
} bed_t;

#define BED_NAME_SIZE 1024
#define BED_INIT_ROWS 1024

bed_t *bed_read(const char *fn, const bool constant_size, const uint32_t size, const bool resize_left, const bool resize_right) {
    bed_t *bed = alloc(sizeof(bed_t));
    bed->m_names = BED_INIT_ROWS;
    bed->names = alloc(sizeof(char *) * bed->m_names);
    kh_bedHash_t *h = kh_init(bedHash);
    gzFile fp;
    kstream_t *ks = NULL;
    int dret;
    unsigned int line = 0, save_errno;
    kstring_t str = { 0, 0, NULL };

    if (NULL == h) {
        fprintf(stderr, "[E::bed_read] Out of memory\n");
        return NULL;
    }
    fp = gzopen(fn, "r");
    if (fp == 0) {
        fprintf(stderr, "[E::bed_read] Unable to open BED file (-b)\n");
        return NULL;
    }
    ks = ks_init(fp);
    if (NULL == ks) {
        fprintf(stderr, "[E::bed_read] Error opening BED file (-b)\n");
        goto fail;
    }
    int ks_len;
    uint32_t prev_size = -1;
    int64_t name_index = 0;
    while ((ks_len = ks_getuntil(ks, KS_SEP_LINE, &str, &dret)) >= 0) {
        char *ref = str.s, *ref_end;
        char strand;
        char name[BED_NAME_SIZE];
        uint32_t beg = 0, end = 0;
        int num = 0;
        khint_t k;
        bedList_t *p;

        if (ks_len == 0)
            continue; // skip blank lines

        line++;
        while (*ref && isspace(*ref)) ref++;
        if ('\0' == *ref) continue;  // Skip blank lines
        if ('#'  == *ref) continue;  // Skip BED file comments
        ref_end = ref;   // look for the end of the reference name
        while (*ref_end && !isspace(*ref_end)) ref_end++;
        if ('\0' != *ref_end) {
            *ref_end = '\0';  // terminate ref and look for start, end, name, score, strand
            num = sscanf(ref_end + 1, "%"SCNu32" %"SCNu32" %s %*s %c",
                         &beg, &end, name, &strand);
        }
        if (1 == num) {  // VCF-style format
            end = beg--; // Counts from 1 instead of 0 for BED files
        }
        if (num < 1 || end < beg) {
            // These two are special lines that can occur in BED files.
            // Check for them here instead of earlier in case someone really
            // has called their reference "browser" or "track".
            if (0 == strcmp(ref, "browser")) continue;
            if (0 == strcmp(ref, "track")) continue;
            if (num < 1) {
                fprintf(stderr,
                        "[E::bed_read] Parse error reading \"%s\" at line %u\n",
                        fn, line);
            } else {
                fprintf(stderr,
                        "[E::bed_read] Parse error reading \"%s\" at line %u : "
                        "end (%"PRIu32") must not be less "
                        "than start (%"PRIu32")\n",
                        fn, line, end, beg);
            }
            errno = 0; // Prevent caller from printing misleading error messages
            goto fail;
        }
        if (num < 4) strand = '.';
        // Put reg in the hash table if not already there
        k = kh_get(bedHash, h, ref);
        if (k == kh_end(h)) { // absent from the hash table
            int ret;
            char *s = strdup(ref);
            if (NULL == s) {
                fprintf(stderr, "[E::bed_read] Out of memory\n");
                goto fail;
            }
            k = kh_put(bedHash, h, s, &ret);
            if (-1 == ret) {
                free(s);
                fprintf(stderr, "[E::bed_read] Out of memory\n");
                goto fail;
            }
            memset(&kh_val(h, k), 0, sizeof(bedList_t));
        }
        p = &kh_val(h, k);

        // Add begin,end to the list
        if (p->n == p->m) {
            p->m = p->m ? p->m<<1 : 4;
            bedItem_t *new_a = realloc(p->a, p->m * sizeof(p->a[0]));
            if (NULL == new_a) {
                fprintf(stderr, "[E::bed_read] Out of memory");
                goto fail;
            }
            p->a = new_a;
        }
        if (size) {
            uint32_t up = size / 2;
            uint32_t down = size / 2 + size % 2;
            if (!resize_left && !resize_right) {
                uint32_t range_size = end - beg;
                beg += range_size / 2;
                end -= range_size / 2;
                if (strand == '-') {
                    beg += range_size % 2;
                    beg = down > beg ? 0 : beg - down;
                    end += up;
                } else {
                    end -= range_size % 2;
                    beg = up > beg ? 0 : beg - up;
                    end += down;
                }
            } else if (resize_left) {
                if (strand == '-') {
                    beg = down > end ? 0 : end - down;
                    end += up;
                } else {
                    end = beg + down;
                    beg = up > beg ? 0 : beg - up;
                }
            } else if (resize_right) {
                if (strand == '-') {
                    end = beg + 1;
                    beg = down > end ? 0 : end - down;
                    end = end + up;
                } else {
                    beg = end - 1;
                    end = beg + down;
                    beg = up > beg ? 0 : beg - up;
                }
            }
        } else if (constant_size) {
            if (prev_size == -1) {
                prev_size = end - beg;
            } else {
                if (prev_size != (end - beg)) {
                    fprintf(stderr, "[E::bed_read] Ranges are not equal size, use -s\n");
                    goto fail;
                }
            }
        }
        if (num < 3 || strcmp(name, ".") == 0) {
            int sret = snprintf(name, BED_NAME_SIZE, "%s:%u-%u", ref, beg, end);
            if (sret < 0 || sret > BED_NAME_SIZE) {
                fprintf(stderr,
                        "[E::bed_read] Failed to store range name in \"%s\" at line %u\n",
                        fn, line);
                errno = 0;
                goto fail;
            }
        }
        if (bed->n_names == bed->m_names) {
            bed->m_names *= 2;
            char **tmp_names = realloc(bed->names, bed->m_names * sizeof(char *));
            if (tmp_names == NULL) {
                fprintf(stderr, "[E::bed_read] Out of memory\n");
                errno = 0;
                goto fail;
            }
            bed->names = tmp_names;
        }
        bed->names[bed->n_names++] = strdup(name);

        p->a[p->n].strand = strand;
        p->a[p->n].namei = name_index++;
        p->a[p->n].beg = beg;
        p->a[p->n++].end = end;
        bed->n++;
    }
    if (!bed->n) {
        fprintf(stderr, "[E::bed_read] Found 0 ranges in BED file\n");
        goto fail;
    }
    if (size) {
        bed->size = size;
    } else if (constant_size) {
        bed->size = prev_size;
    }

    if (gzclose(fp) != Z_OK) {
        fp = NULL;
        fprintf(stderr, "[E::bed_read] Error closing BED file (-b)\n");
        goto fail;
    }
    ks_destroy(ks);
    free(str.s);
    bed->data = h;
    return bed;
 fail:
    save_errno = errno;
    if (ks) ks_destroy(ks);
    if (fp) gzclose(fp);
    free(str.s);
    errno = save_errno;
    return NULL;
}

// subset ----------------------------------------------------------------------

static int *bed_index_core(int n, bedItem_t *a) {
    int i, j, l, *idx, *new_idx;
    l = 0; idx = 0;
    for (i = 0; i < n; ++i) {
        uint32_t beg, end;
        beg = a[i].beg >> LIDX_SHIFT; end = a[i].end >> LIDX_SHIFT;
        if (l < end + 1) {
            int old_l = l;
            l = end + 1;
            kroundup32(l);
            new_idx = realloc(idx, l * sizeof(*idx));
            if (!new_idx) {
                free(idx);
                return NULL;
            }
            idx = new_idx;

            for (j = old_l; j < l; ++j)
                idx[j] = -1;
        }

        for (j = beg; j < end+1; ++j)
            if (idx[j] < 0)
                idx[j] = i;
    }
    return idx;
}

static void bed_index(kh_bedHash_t *h) {
    khint_t k;
    for (k = 0; k < kh_end(h); ++k) {
        if (kh_exist(h, k)) {
            bedList_t *p = &kh_val(h, k);
            if (p->idx) free(p->idx);
            ks_introsort(bedItem_t, p->n, p->a);
            p->idx = bed_index_core(p->n, p->a);
        }
    }
}

void bed_unify(kh_bedHash_t *h) {

    int i, j, new_n;
    bedList_t *p;

    for (i = kh_begin(h); i < kh_end(h); i++) {
        if (!kh_exist(h,i) || !(p = &kh_val(h,i)) || !(p->n))
            continue;

        for (new_n = 0, j = 1; j < p->n; j++) {
            if (p->a[new_n].end < p->a[j].beg) {
                p->a[++new_n] = p->a[j];
            } else {
                if (p->a[new_n].end < p->a[j].end)
                    p->a[new_n].end = p->a[j].end;
            }
        }

        p->n = ++new_n;

        for (int k = 0; k < p->n; k++) {
          p->a[k].strand = '.';
          p->a[k].namei = -1;
        }
    }
}

// values ----------------------------------------------------------------------

static void help_values(void) {
    printf(
        "bwtk v%s  Copyright (C) %s  Benjamin Jean-Marie Tremblay\n"
        "bwtk values [options] -i <file.bw> -b <ranges.bed> -o <values.tsv>\n"
        "    -i    Input bigWig\n"
        "    -b    BED file with ranges to extract values from\n"
        "    -o    Output values (use '-' for stdout)\n"
        "    -s    Desired size of ranges, will be resized from the centre\n"
        "    -l    Resize ranges from the left with -s (based on strand)\n"
        "    -r    Resize ranges from the right with -s (based on strand)\n"
        "    -h    Print this message and exit\n"
        , BWTK_VERSION, BWTK_YEAR
    );
}

static int values(int argc, char **argv) {
    if (argc == 1) {
        fprintf(stderr, "[E::values] No arguments provided\n");
        help_values();
        return EXIT_FAILURE;
    }
    char *bedfn = NULL;
    bed_t *bed;
    bigWigFile_t *bw = NULL;
    int size = 0;
    bool left = false, right = false;
    FILE *fout = NULL;
    int opt;
    while ((opt = getopt(argc, argv, "i:b:o:s:lrh")) != -1) {
        switch (opt) {
            case 'i':
                if (!bwIsBigWig(optarg, NULL)) {
                    fprintf(stderr, "[E::values] Not a bigWig file (-i): '%s'\n", optarg);
                    return EXIT_FAILURE;
                }
                bw = bwOpen(optarg, NULL, "r");
                if (!bw) {
                    fprintf(stderr, "[E::values] Unable to open bigWig (-i): '%s'\n", optarg);
                    return EXIT_FAILURE;
                }
                break;
            case 'b':
                bedfn = optarg;
                break;
            case 'o':
                fout = strcmp(optarg, "-") ? fopen(optarg, "w") : fdopen(1, "w");
                if (fout == NULL) {
                    fprintf(stderr, "[E::values] Unable to open file (-o): %s\n", strerror(errno));
                    return EXIT_FAILURE;
                }
                break;
            case 's':
                size = atoi(optarg);
                if (size <= 0) {
                    fprintf(stderr, "[E::values] -s must be a positive integer\n");
                    return EXIT_FAILURE;
                }
                break;
            case 'l':
                left = true;
                break;
            case 'r':
                right = true;
                break;
            case 'h':
                help_values();
                return EXIT_SUCCESS;
            default:
                return EXIT_FAILURE;
        }
    }
    if (bedfn != NULL) {
        // If I don't index, I can keep the ranges in their original order
        // per chromosome. But after that, there is no way to return to the
        // original order! So might as well index. For consistency, when
        // printing the results, use the order of chromosomes listed in the
        // bigWig.
        bed = bed_read(bedfn, true, (uint32_t) size, left, right);
        if (bed == NULL) return EXIT_FAILURE;
        bed_index(bed->data);
    } else {
        fprintf(stderr, "[E::values] Missing -b\n");
        return EXIT_FAILURE;
    }
    if (bw == NULL) {
        fprintf(stderr, "[E::values] Missing -i\n");
        return EXIT_FAILURE;
    }
    if (fout == NULL) {
        fprintf(stderr, "[E::values] Missing -o\n");
        return EXIT_FAILURE;
    }

    fputs("Range", fout);
    for (int32_t i = 0; i < bed->size; i++) {
        fprintf(fout, "\t%u", i + 1);
    }
    fputc('\n', fout);

    int64_t nranges = 0;
    uint32_t chromLen;
    char *chromName;
    khint_t k;
    bedList_t *b;
    bwOverlappingIntervals_t *bwVals;
    double v;
    for (int64_t i = 0; i < bw->cl->nKeys; i++) {
        chromName = bw->cl->chrom[i];
        chromLen = bw->cl->len[i];
        k = kh_get(bedHash, bed->data, chromName);
        if (k != kh_end(bed->data)) {
            b = &kh_val(bed->data, k);
            // TODO: Lots of malloc()/free() in this libBigWig code, could be improved
            for (int j = 0; j < b->n; j++) {
                nranges++;
                bwVals = bwGetValues(bw, chromName, b->a[j].beg, b->a[j].end, 1);
                if (bwVals == NULL) {
                    fprintf(stderr, "[E::values] Failed to extract ranges from bigWig: %s:%u-%u\n", chromName, b->a[j].beg, b->a[j].end);
                    return EXIT_FAILURE;
                }
                fputs(bed->names[b->a[j].namei], fout);
                if (b->a[j].beg == 0 && (b->a[j].end - b->a[j].beg) < size) {
                    for (uint32_t h = 0; h < size - (b->a[j].end - b->a[j].beg); h++) {
                        fputs("\tnan", fout);
                    }
                }
                if (b->a[j].strand == '-') {
                    for (uint32_t h = bwVals->l; h > 0; h--) {
                        v = isnan(bwVals->value[h - 1]) ? 0.0 : (double) bwVals->value[h - 1];
                        fprintf(fout, "\t%f", v);
                    }
                } else {
                    for (uint32_t h = 0; h < bwVals->l; h++) {
                        v = isnan(bwVals->value[h]) ? 0.0 : (double) bwVals->value[h];
                        fprintf(fout, "\t%f", v);
                    }
                }
                if (b->a[j].end == chromLen && (b->a[j].end - b->a[j].beg) < size) {
                    for (uint32_t h = 0; h < size - (b->a[j].end - b->a[j].beg); h++) {
                        fputs("\tnan", fout);
                    }
                }
                fputc('\n', fout);
                bwDestroyOverlappingIntervals(bwVals);
            }
        }
    }
    if (!nranges) {
        fprintf(stderr, "[E::values] Found zero valid BED ranges in bigWig\n");
        return EXIT_FAILURE;
    }
    if (nranges < bed->n) {
        fprintf(stderr, "[W::values] Found only %lld/%lld ranges in bigWig\n", nranges, bed->n);
    }

    fclose(fout);
    bwClose(bw);
    return EXIT_SUCCESS;
}

// bw2bg -----------------------------------------------------------------------

static void help_bw2bg(void) {
    printf(
        "bwtk v%s  Copyright (C) %s  Benjamin Jean-Marie Tremblay\n"
        "bwtk bg2bw [options] -i <file.bw> -o <file.bedGraph.gz>\n"
        "    -i    Input bigWig\n"
        "    -o    Output gzipped bedGraph, use '-' for stdout (ungzipped)\n"
        "    -m    Multiply scores by this value [1]\n"
        "    -a    Add this value to scores [0]\n"
        "    -l    log10-transform scores\n"
        "    -h    Print this message and exit\n"
        , BWTK_VERSION, BWTK_YEAR
    );
}

static union fout_t {
    FILE   *f;
    gzFile gz;
} fout;

static int bw2bg(int argc, char **argv) {
    if (argc == 1) {
        fprintf(stderr, "[E::bw2bg] No arguments provided, use -h for usage\n");
        help_bw2bg();
        return EXIT_FAILURE;
    }
    int opt;
    bool tofile = true, do_log10 = false;
    float mult = 1.0f, add = 0.0f;
    fout.f = NULL;
    bigWigFile_t *bw = NULL;
    while ((opt = getopt(argc, argv, "i:o:m:a:lh")) != -1) {
        switch (opt) {
            case 'i':
                if (!bwIsBigWig(optarg, NULL)) {
                    fprintf(stderr, "[E::bw2bg] Not a bigWig file (-i): '%s'\n", optarg);
                    return EXIT_FAILURE;
                }
                bw = bwOpen(optarg, NULL, "r");
                if (!bw) {
                    fprintf(stderr, "[E::bw2bg] Unable to open bigWig (-i): '%s'\n", optarg);
                    return EXIT_FAILURE;
                }
                break;
            case 'o':
                if (strcmp(optarg, "-") == 0) {
                    tofile = false;
                    fout.f = fdopen(1, "w");
                    if (fout.f == NULL) {
                        fprintf(stderr, "[E::bw2bg] Failed to fdopen stdout: %s\n", strerror(errno));
                        return EXIT_FAILURE;
                    }
                } else {
                    fout.gz = gzopen(optarg, "wb");
                    if (fout.gz == NULL) {
                        int e;
                        fprintf(stderr, "[E::bw2bg] Failed to open file '%s': %s\n", optarg, gzerror(fout.gz, &e));
                        return EXIT_FAILURE;
                    }
                }
                break;
            case 'm':
                mult = atof(optarg);
                if (mult == 0 && errno == ERANGE) {
                    fprintf(stderr, "[E::bw2bg] Unable to parse '-m': %s\n", strerror(errno));
                    return EXIT_FAILURE;
                } else if (mult == 0) {
                    fprintf(stderr, "[E::bw2bg] '-m' must be nonzero\n");
                    return EXIT_FAILURE;
                }
                break;
            case 'a':
                add = atof(optarg);
                if (add == 0 && errno == ERANGE) {
                    fprintf(stderr, "[E::bw2bg] Unable to parse '-a': %s\n", strerror(errno));
                    return EXIT_FAILURE;
                }
                break;
            case 'l':
                do_log10 = true;
                break;
            case 'h':
                help_bw2bg();
                return EXIT_SUCCESS;
            default:
                return EXIT_FAILURE;
        }
    }
    if (fout.f == NULL) {
        fprintf(stderr, "[E::bw2bg] Missing -o\n");
        return EXIT_FAILURE;
    }
    if (bw == NULL) {
        fprintf(stderr, "[E::bw2bg] Missing -i\n");
        return EXIT_FAILURE;
    }

    // TODO: the actual conversion

    if (tofile && gzclose (fout.gz) != Z_OK) {
        int e;
        fprintf(stderr, "[E::bw2bg] Failed to close output file: %s\n", gzerror(fout.gz, &e));
        return EXIT_FAILURE;
    } else if (fclose(fout.f) != 0) {
        fprintf(stderr, "[E::bw2bg] Failed to fclose output stream: %s\n", strerror(errno));
        return EXIT_FAILURE;
    }
    bwClose(bw);
    return EXIT_SUCCESS;
}

// bg2bw -----------------------------------------------------------------------

#define CHROMSIZES_MEM   128
#define CHROMNAME_SIZE  1024
#define RANGES_MEM      8192  // how many bedGraph rows to keep in memory
#define PRESET_GENOMES  "tair10"

static void help_bg2bw(void) {
    printf(
        "bwtk v%s  Copyright (C) %s  Benjamin Jean-Marie Tremblay\n"
        "bwtk bg2bw [options] -g <chrom.sizes> -i <file.bedGraph[.gz]> -o <file.bw>\n"
        "    -i    Input bedGraph (can be gzipped), use '-' for stdin\n"
        "    -o    Output bigWig\n"
        "    -g    Genome chrom.sizes or genome.fa.fai file\n"
        "    -m    Multiply scores by this value [1]\n"
        "    -a    Add this value to scores [0]\n"
        "    -l    log10-transform scores\n"
        "    -p    Use a preset genome instead of -g [%s]\n"
        "    -u    When using -p, use UCSC-style names (default: Ensembl)\n"
        "    -h    Print this message and exit\n"
        , BWTK_VERSION, BWTK_YEAR, PRESET_GENOMES
    );
}

typedef struct {
    int m, n;
    char **chroms;
    uint32_t *sizes;
} chromSizes_t;

typedef struct {
    int m, n, times_added;
    char **chrom;
    uint32_t *start, *end;
    float *val;
} ranges_t;

static inline void addBwInterval(bigWigFile_t *bw, ranges_t *ranges) {
    if (ranges->times_added == 0) {
        bwAddIntervals(bw, (const char * const *) ranges->chrom, ranges->start, ranges->end, ranges->val, ranges->n);
        ranges->times_added = 1;
    } else {
        bwAppendIntervals(bw, ranges->start, ranges->end, ranges->val, ranges->n);
        ranges->times_added++;
    }
    ranges->n = 0;
}

static void initChromSizes(chromSizes_t *cs, const int m) {
    cs->m = m;
    cs->chroms = alloc(sizeof(char *) * cs->m);
    cs->sizes = alloc(sizeof(uint32_t) * cs->m);
}

static chromSizes_t *readChromSizes(gzFile cs) {
    chromSizes_t *chromSizes = alloc(sizeof(chromSizes_t));
    initChromSizes(chromSizes, CHROMSIZES_MEM);

    kstream_t *ks_cs = ks_init(cs);
    kstring_t str = { 0, 0, NULL };
    int dret, ks_len;
    int cs_lines = 0;
    while ((ks_len = ks_getuntil(ks_cs, KS_SEP_LINE, &str, &dret)) >= 0) {
        if (ks_len == 0) continue;
        char *ref = str.s;
        while (*ref && isspace(*ref)) ref++;
        if (*ref == 0 || *ref == '#') continue;
        char chr[CHROMNAME_SIZE];
        uint32_t size;
        int num = sscanf(ref, "%s %"SCNu32, chr, &size);
        if (num < 2) {
            fprintf(stderr, "[E::bg2bw] Expected at least two columns in chrom.sizes (-g)\n");
            return (chromSizes_t *) NULL;
        }
        if (chromSizes->n + 1 > chromSizes->m) {
            chromSizes->m *= 2;
            chromSizes->chroms = realloc(chromSizes->chroms, sizeof(char *) * chromSizes->m);
            chromSizes->sizes = realloc(chromSizes->sizes, sizeof(uint32_t) * chromSizes->m);
            if (chromSizes->chroms == NULL || chromSizes->sizes == NULL) {
                fprintf(stderr, "[E::readChromSizes] Out of memory\n");
                return (chromSizes_t *) NULL;
            }
        }
        chromSizes->chroms[chromSizes->n] = strdup(chr);
        chromSizes->sizes[chromSizes->n++] = size;
        cs_lines++;
    }
    if (!cs_lines) {
        fprintf(stderr, "[E::readChromSizes] Empty chrom.sizes file (-g)\n");
        return (chromSizes_t *) NULL;
    }
    ks_destroy(ks_cs);
    free(str.s);
    gzclose(cs);

    return chromSizes;
}

static void copyStrings(const char **source, char **dest, const int n) {
    for (int i = 0; i < n; i++) {
        dest[i] = strdup(source[i]);
    }
}

static void copyInts(const uint32_t *source, uint32_t *dest, const int n) {
    for (int i = 0; i < n; i++) {
        dest[i] = source[i];
    }
}

static chromSizes_t *chromSizesFromPreset(const char *preset, const bool ucsc) {
    chromSizes_t *chromSizes = alloc(sizeof(chromSizes_t));
    if (strcmp(preset, "tair10") == 0) {
#include "chromsizes/tair10.c"
    } else {
        fprintf(stderr, "[E::chromSizesFromPreset] Unknown preset genome '%s'\n", preset);
        return (chromSizes_t *) NULL;
    }
    return chromSizes;
}

static int bg2bw(int argc, char **argv) {
    if (argc == 1) {
        fprintf(stderr, "[E::bg2bw] No arguments provided, use -h for usage\n");
        help_bg2bw();
        return EXIT_FAILURE;
    }
    int opt;
    gzFile bg = NULL, cs = NULL;
    bigWigFile_t *bw = NULL;
    float mult = 1.0f, add = 0.0f;
    char *preset = NULL;
    bool do_log10 = false, ucsc = false;
    while ((opt = getopt(argc, argv, "i:o:g:m:a:lp:uh")) != -1) {
        switch (opt) {
            case 'i':
                bg = strcmp(optarg, "-") ? gzopen(optarg, "r") : gzdopen(fileno(stdin), "r");
                if (bg == 0) {
                    fprintf(stderr, "[E::bg2bw] Unable to open input (-i) '%s'\n", optarg);
                    return EXIT_FAILURE;
                }
                break;
            case 'o':
                bw = bwOpen(optarg, NULL, "w");
                if (!bw) {
                    fprintf(stderr, "[E::bg2bw] Unable to create output (-o) '%s'\n", optarg);
                    return EXIT_FAILURE;
                }
                break;
            case 'g':
                cs = gzopen(optarg, "r");
                if (cs == 0) {
                    fprintf(stderr, "[E::bg2bw] Unable to open chrom.sizes (-g) '%s'\n", optarg);
                    return EXIT_FAILURE;
                }
                break;
            case 'm':
                mult = atof(optarg);
                if (mult == 0 && errno == ERANGE) {
                    fprintf(stderr, "[E::bg2bw] Unable to parse '-m': %s\n", strerror(errno));
                    return EXIT_FAILURE;
                } else if (mult == 0) {
                    fprintf(stderr, "[E::bg2bw] '-m' must be nonzero\n");
                    return EXIT_FAILURE;
                }
                break;
            case 'a':
                add = atof(optarg);
                if (add == 0 && errno == ERANGE) {
                    fprintf(stderr, "[E::bg2bw] Unable to parse '-a': %s\n", strerror(errno));
                    return EXIT_FAILURE;
                }
                break;
            case 'l':
                do_log10 = true;
                break;
            case 'p':
                preset = optarg;
                break;
            case 'u':
                ucsc = true;
                break;
            case 'h':
                help_bg2bw();
                return EXIT_SUCCESS;
            default:
                return EXIT_FAILURE;
        }
    }
    if (bg == NULL) {
        fprintf(stderr, "[E::bg2bw] Missing input (-i)\n");
        return EXIT_FAILURE;
    }
    if (bw == NULL) {
        fprintf(stderr, "[E::bg2bw] Missing output (-o)\n");
        return EXIT_FAILURE;
    }
    if (cs == NULL && preset == NULL) {
        fprintf(stderr, "[E::bg2bw] Missing chrom.sizes (-g or -p)\n");
        return EXIT_FAILURE;
    }

    ranges_t *ranges = alloc(sizeof(ranges_t));
    ranges->m = RANGES_MEM;
    ranges->chrom = alloc(sizeof(char *) * ranges->m);
    ranges->start = alloc(sizeof(uint32_t) * ranges->m);
    ranges->end = alloc(sizeof(uint32_t) * ranges->m);
    ranges->val = alloc(sizeof(float) * ranges->m);

    kstring_t str = { 0, 0, NULL };
    int dret, ks_len;

    chromSizes_t *chromSizes;
    if (preset != NULL) {
        chromSizes = chromSizesFromPreset(preset, ucsc);
    } else {
        chromSizes = readChromSizes(cs);
    }
    if (chromSizes == NULL) {
        return EXIT_FAILURE;
    }

    if (bwCreateHdr(bw, 10)) {  // Lowering this value doesn't lower mem usage in bwClose()
        fprintf(stderr, "[E::bg2bw] Failed to init bigWig header\n");
        return EXIT_FAILURE;
    }
    bw->cl = bwCreateChromList((const char * const *) chromSizes->chroms, chromSizes->sizes, chromSizes->n);
    if (!bw->cl) {
        fprintf(stderr, "[E::bg2bw] Failed to create bigWig chrom list\n");
        return EXIT_FAILURE;
    }
    if (bwWriteHdr(bw)) {
        fprintf(stderr, "[E::bg2bw] Failed to write bigWig header\n");
        return EXIT_FAILURE;
    }

    khash_t(str_m) *bgChroms = kh_init(str_m);
    int absent;
    khint64_t k;
    for (int i = 0; i < chromSizes->n; i++) {
        k = kh_put(str_m, bgChroms, chromSizes->chroms[i], &absent);
        if (absent == -1) {
            fprintf(stderr, "[E::bg2bw] Failed to create bgChroms hash table\n");
            return EXIT_FAILURE;
        }
        kh_val(bgChroms, k) = -1;
    }

    kstream_t *ks_bg = ks_init(bg);
    char *chrom_last = NULL;
    int64_t frow = 0, erow = 0, nranges = 0, bwdumps = 0;
    uint32_t end_last = 0;
    while ((ks_len = ks_getuntil(ks_bg, KS_SEP_LINE, &str, &dret)) >= 0) {
        frow++;
        if (ks_len == 0) continue;
        char *ref = str.s;
        while (*ref && isspace(*ref)) ref++;
        if (*ref == 0 || *ref == '#') continue;
        char chr[CHROMNAME_SIZE];
        uint32_t start, end;
        float val;
        int num = sscanf(ref, "%s %"SCNu32" %"SCNu32" %f", chr, &start, &end, &val);
        if (num < 4) {
            fprintf(stderr, "[E::bg2bw] Error: Malformed bedGraph data on line %lld\n", frow);
            return EXIT_FAILURE;
        }
        val += add;
        val *= mult;
        if (do_log10) val = log10f(val);
        erow++;
        k = kh_get(str_m, bgChroms, chr);
        if (k == kh_end(bgChroms)) {
            fprintf(stderr, "[E::bg2bw] Found unknown chrom '%s' on line %lld\n", chr, frow);
            return EXIT_FAILURE;
        }
        if (kh_val(bgChroms, k) != -1 && kh_val(bgChroms, k) < (erow - 1)) {
            fprintf(stderr, "[E::bg2bw] bedGraph file must be sorted (found out of order chroms)\n");
            return EXIT_FAILURE;
        }
        kh_val(bgChroms, k) = erow;
        if (end <= start || start < 0) {
            fprintf(stderr, "[E::bg2bw] Bad bedGraph range on line %lld (start=%u end=%u)\n", frow, start, end);
            return EXIT_FAILURE;
        }
        // TODO: Make sure end doesn't go beyond chrom limit, libBigWig does no input checks
        if (chrom_last == NULL || strcmp(chrom_last, chr) == 0) {
            if (ranges->n == ranges->m) {
                addBwInterval(bw, ranges);
                bwdumps++;
            }
            if (end_last > 0 && start < end_last) {
                fprintf(stderr, "[E::bg2bw] bedGraph must be sorted and/or not contain overlapping ranges (lines %lld-%lld)\n", frow - 1, frow);
                return EXIT_FAILURE;
            }
            ranges->chrom[ranges->n] = (char *) kh_key(bgChroms, k);
            ranges->start[ranges->n] = start;
            ranges->end[ranges->n] = end;
            ranges->val[ranges->n++] = val;
        } else {
            if (ranges->n) {
                addBwInterval(bw, ranges);
                bwdumps++;
            }
            ranges->chrom[0] = (char *) kh_key(bgChroms, k);
            ranges->start[0] = start;
            ranges->end[0] = end;
            ranges->val[0] = val;
            ranges->n = 1;
            ranges->times_added = 0;
        }
        end_last = end;
        chrom_last = (char *) kh_key(bgChroms, k);
        nranges++;
    }
    free(str.s);
    ks_destroy(ks_bg);
    gzclose(bg);

    if (ranges->n) {
        addBwInterval(bw, ranges);
        bwdumps++;
    }
    if (!bwdumps) {
        fprintf(stderr, "[E::bg2bw] Failed to find any ranges in bedGraph\n");
        return EXIT_FAILURE;
    }

    fprintf(stderr, "Read %lld ranges, dumped to bw %lld times.\n", nranges, bwdumps);
    fprintf(stderr, "Peak mem before indexing: %'.2f MB\n", ((double) peak_mem() / 1024.0) / 1024.0);

    // TODO: Inside bwClose() is bwFinalize() -> constructZoomLevels(), which for a full
    // bedGraph increases the memory usage by 100x
    bwClose(bw);
    fprintf(stderr, "Peak mem after indexing: %'.2f MB\n", ((double) peak_mem() / 1024.0) / 1024.0);

    return EXIT_SUCCESS;
}

// main ------------------------------------------------------------------------

static void help(void) {
        /* "    merge      Average multiple bigWig files together\n" */
        /* "    adjust     Perform an operation on bigWig values\n" */
    printf(
        "bwtk v%s  Copyright (C) %s  Benjamin Jean-Marie Tremblay\n"
        "Usage:  bwtk <subcommand> [options]\n"
        "Available subcommands:\n"
        "    bw2bg      Convert a bigWig file to bedGraph\n"
        "    bg2bw      Convert a bedGraph file to bigWig\n"
        "    values     Return bigWig values from overlapping BED ranges\n"
        "    subset     Subset a bigWig using a BED file\n"
        "    help       Print this message and exit\n"
        "    version    Print the version number and exit\n"
        , BWTK_VERSION, BWTK_YEAR
    );
}

int main(int argc, char **argv) {

    if (argc < 2) {
        fprintf(stderr, "[E::main] Missing subcommand, 'help' for usage\n");
        return EXIT_FAILURE;
    }
    if (strcmp(argv[1], "version") == 0) {
        fprintf(stderr, "v%s\n", BWTK_VERSION);
        return EXIT_SUCCESS;
    } else if (strcmp(argv[1], "help") == 0) {
        help();
        return EXIT_SUCCESS;
    } else if (strcmp(argv[1], "bw2bg") == 0) {
        return bw2bg(argc - 1, argv + 1);
    } else if (strcmp(argv[1], "bg2bw") == 0) {
        return bg2bw(argc - 1, argv + 1);
    } else if (strcmp(argv[1], "values") == 0) {
        return values(argc - 1, argv + 1);
    } else if (strcmp(argv[1], "subset") == 0) {
    } else {
        fprintf(stderr, "[E::main] Unknown subcommand '%s'; 'help' for usage\n", argv[1]);
        return EXIT_FAILURE;
    }

    return EXIT_SUCCESS;
}

