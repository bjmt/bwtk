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

#define BWTK_VERSION "1.6.0"
#define BWTK_YEAR "2025"

// common ----------------------------------------------------------------------

static union fout_t {
    FILE   *f;
    gzFile gz;
} fout;

typedef struct {
    int64_t x, y;
} int64_t2;

KSTREAM_INIT(gzFile, gzread, 8192)
KHASH_MAP_INIT_STR(str_m, int64_t2)
KHASH_MAP_INIT_STR(len_m, int32_t)

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

#define min(a, b) ((a) > (b) ? (b) : (a))
#define max(a, b) ((a) < (b) ? (b) : (a))

// Increasing this value only increases the memory usage without actually speeding
// anything up.
#define BW_ITERATOR_CHUNKS 5

// Some BED handling code adapted from HTSlib/bedidx.c

#define LIDX_SHIFT 13
#define ALL 0
#define FILTERED 1

typedef struct {
    int m, n;
    char **chroms;
    uint32_t *sizes;
} chromSizes_t;

#define RANGES_MEM      8192

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

static bed_t *readBED(const char *fn, const bool constant_size, const uint32_t size, const bool resize_left, const bool resize_right) {
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
        fprintf(stderr, "[E::readBED] Out of memory\n");
        return NULL;
    }
    fp = strcmp(fn, "-")? gzopen(fn, "r") : gzdopen(fileno(stdin), "r");
    if (fp == 0) {
        fprintf(stderr, "[E::readBED] Unable to open BED file (-b)\n");
        return NULL;
    }
    ks = ks_init(fp);
    if (NULL == ks) {
        fprintf(stderr, "[E::readBED] Error opening BED file (-b)\n");
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
            continue;

        line++;
        while (*ref && isspace(*ref)) ref++;
        if ('\0' == *ref) continue;
        if ('#'  == *ref) continue;
        ref_end = ref;
        while (*ref_end && !isspace(*ref_end)) ref_end++;
        if ('\0' != *ref_end) {
            *ref_end = '\0';
            num = sscanf(ref_end + 1, "%"SCNu32" %"SCNu32" %s %*s %c",
                         &beg, &end, name, &strand);
        }
        if (1 == num) {
            end = beg--;
        }
        if (num < 1 || end < beg) {
            if (0 == strcmp(ref, "browser")) continue;
            if (0 == strcmp(ref, "track")) continue;
            if (num < 1) {
                fprintf(stderr,
                        "[E::readBED] Parse error reading \"%s\" at line %u\n",
                        fn, line);
            } else {
                fprintf(stderr,
                        "[E::readBED] Parse error reading \"%s\" at line %u : "
                        "end (%"PRIu32") must not be less "
                        "than start (%"PRIu32")\n",
                        fn, line, end, beg);
            }
            errno = 0;
            goto fail;
        }
        if (num < 4) strand = '.';
        k = kh_get(bedHash, h, ref);
        if (k == kh_end(h)) {
            int ret;
            char *s = strdup(ref);
            if (NULL == s) {
                fprintf(stderr, "[E::readBED] Out of memory\n");
                goto fail;
            }
            k = kh_put(bedHash, h, s, &ret);
            if (-1 == ret) {
                free(s);
                fprintf(stderr, "[E::readBED] Out of memory\n");
                goto fail;
            }
            memset(&kh_val(h, k), 0, sizeof(bedList_t));
        }
        p = &kh_val(h, k);

        if (p->n == p->m) {
            p->m = p->m ? p->m<<1 : 4;
            bedItem_t *new_a = realloc(p->a, p->m * sizeof(p->a[0]));
            if (NULL == new_a) {
                fprintf(stderr, "[E::readBED] Out of memory");
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
                    fprintf(stderr, "[E::readBED] Ranges are not equal size, use -s\n");
                    goto fail;
                }
            }
        }
        if (num < 3 || strcmp(name, ".") == 0) {
            int sret = snprintf(name, BED_NAME_SIZE, "%s:%u-%u", ref, beg, end);
            if (sret < 0 || sret > BED_NAME_SIZE) {
                fprintf(stderr,
                        "[E::readBED] Failed to store range name in \"%s\" at line %u\n",
                        fn, line);
                errno = 0;
                goto fail;
            }
        }
        if (bed->n_names == bed->m_names) {
            bed->m_names *= 2;
            char **tmp_names = realloc(bed->names, bed->m_names * sizeof(char *));
            if (tmp_names == NULL) {
                fprintf(stderr, "[E::readBED] Out of memory\n");
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
        fprintf(stderr, "[E::readBED] Found 0 ranges in BED file\n");
        goto fail;
    }
    if (size) {
        bed->size = size;
    } else if (constant_size) {
        bed->size = prev_size;
    }

    if (gzclose(fp) != Z_OK) {
        fp = NULL;
        fprintf(stderr, "[E::readBED] Error closing BED file (-b)\n");
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

// adjust ----------------------------------------------------------------------

static int bed_unify(kh_bedHash_t *h) {

    int64_t total_n = 0;
    int i, j, new_n;
    bedList_t *p;

    for (i = kh_begin(h); i < kh_end(h); i++) {
        if (!kh_exist(h,i) || !(p = &kh_val(h,i)) || !(p->n))
            continue;

        for (new_n = 0, j = 1; j < p->n; j++) {
            if (p->a[new_n].end < p->a[j].beg) {
                p->a[++new_n] = p->a[j];
            } else if (p->a[new_n].end < p->a[j].end) {
                p->a[new_n].end = p->a[j].end;
            }
        }

        p->n = ++new_n;

        for (int k = 0; k < p->n; k++) {
          p->a[k].strand = '.';
          p->a[k].namei = -1;
        }
        total_n += (int64_t) p->n;
    }
    return total_n;
}

static void help_adjust(void) {
    printf(
        "bwtk v%s  Copyright (C) %s  Benjamin Jean-Marie Tremblay\n"
        "bwtk adjust [options] -i <in.bw> -o <out.bw>\n"
        "    -i    Input bigWig file\n"
        "    -o    Output bigWig file\n"
        "    -B    Output as bedGraph.gz ('-o-' for ungzipped stdout)\n"
        "    -b    Subset to ranges in a BED file ('-' for stdin)\n"
        "    -r    Subset to a single range (chrName:X-Y)\n"
        "    -a    Add this value to scores [0]\n"
        "    -m    Multiply scores by this value [1]\n"
        "    -l    log10-transform scores\n"
        "    -t    Trim values above this max [Inf]\n"
        "    -s    Step size for binning [0]\n"
        "    -h    Print this message and exit\n"
        "Order of operations: a -> m -> l -> t -> s\n"
        , BWTK_VERSION, BWTK_YEAR
    );
}

static void addBgInterval(ranges_t *ranges, const bool tofile) {
    for (int64_t i = 0; i < ranges->n; i++) {
        if (tofile) {
            gzprintf(fout.gz, "%s\t%"PRIu32"\t%"PRIu32"\t%g\n", ranges->chrom[i], ranges->start[i], ranges->end[i], (double) ranges->val[i]);
        } else {
            fprintf(fout.f, "%s\t%"PRIu32"\t%"PRIu32"\t%g\n", ranges->chrom[i], ranges->start[i], ranges->end[i], (double) ranges->val[i]);
        }
    }
    ranges->n = 0;
}

// Some more code adapted from HTSlib:

static inline uint64_t pushDigit(uint64_t i, char c) {
    // ensure subtraction occurs first, avoiding overflow for >= MAX-48 or so
    int digit = c - '0';
    return 10 * i + digit;
}

uint32_t parseDecimal(const char *str, char **strend) {
    uint64_t n = 0;
    int digits = 0, decimals = 0, e = 0, lost = 0;
    char sign = '+', esign = '+';
    const char *s, *str_orig = str;

    while (isspace(*str)) str++;
    s = str;

    if (*s == '+' || *s == '-') sign = *s++;
    while (*s)
        if (isdigit(*s)) digits++, n = pushDigit(n, *s++);
        else if (*s == ',') s++;
        else break;

    if (*s == '.') {
        s++;
        while (isdigit(*s)) decimals++, digits++, n = pushDigit(n, *s++);
    }

    switch (*s) {
        case 'e': case 'E':
            s++;
            if (*s == '+' || *s == '-') esign = *s++;
            while (isdigit(*s)) e = pushDigit(e, *s++);
            if (esign == '-') e = -e;
            break;

        case 'k': case 'K': e += 3; s++; break;
        case 'm': case 'M': e += 6; s++; break;
        case 'g': case 'G': e += 9; s++; break;
    }

    e -= decimals;
    while (e > 0) n *= 10, e--;
    while (e < 0) lost += n % 10, n /= 10, e++;

    if (lost > 0) {
        fprintf(stderr, "[E::parseDecimal] Discarding fractional part of %.*s", (int)(s - str), str);
    }

    if (strend) {
        // Set to the original input str pointer if not valid number syntax
        *strend = (digits > 0)? (char *)s : (char *)str_orig;
    } else if (digits == 0) {
        fprintf(stderr, "[E::parseDecimal] Invalid numeric value %.8s[truncated]", str);
    } else if (*s) {
        fprintf(stderr, "[E::parseDecimal] Ignoring unknown characters after %.*s[%s]", (int)(s - str), str, s);
    }

    if (n >= 0xFFFFFFFF) {
        fprintf(stderr, "[E::parseDecimal] Coordinate value exceeds uint32_t limit: %lld\n", n);
    }
    return (sign == '+')? n : -n;
}

static uint32_t parseSingleRange(const char *s, const bigWigFile_t *bw, uint32_t *beg, uint32_t *end) {
    char *hyphen;
    char *colon = strrchr(s, ':');
    uint32_t tid = -1;
    if (colon == NULL) {
        tid = bwGetTid(bw, s);
        if (tid == -1) {
            fprintf(stderr, "[E::parseSingleRange] Could not find sequence name in bigWig: '%s'\n", s);
            return -1;
        }
        *beg = 0; *end = bw->cl->len[tid];
        return tid;
    }
    *colon = 0;
    tid = bwGetTid(bw, s);
    if (tid == -1) {
        fprintf(stderr, "[E::parseSingleRange] Could not find sequence name in bigWig: '%s'\n", s);
        return -1;
    }

    *beg = parseDecimal(colon+1, &hyphen) - 1;
    if (*beg < 0) *beg = 0;

    if (*hyphen == '\0') {
        *end = *beg + 1;
    } else if (*hyphen == '-') {
        *end = parseDecimal(hyphen+1, NULL);
    } else {
        fprintf(stderr, "Unknown separator: '%c'\n", *hyphen);
        return -1;
    }

    if (*beg >= *end) {
        fprintf(stderr, "[E::parseSingleRange] Range start cannot be equal to or greater than end");
        return -1;
    }
    return tid;
}

static int adjust(int argc, char **argv) {
    if (argc == 1) {
        fprintf(stderr, "[E::adjust] No arguments provided\n");
        help_adjust();
        return EXIT_FAILURE;
    }
    if (bwInit(1<<17) != 0) {
        fprintf(stderr, "[E::adjust] Unable to init curl buffer\n");
        return EXIT_FAILURE;
    }
    char *bedfn = NULL, *outfn = NULL, *rRange = NULL;
    bed_t *bed;
    bigWigFile_t *bw_in = NULL, *bw_out = NULL;
    bool do_log10 = false, use_trim = false;
    bool bg = false, tofile = true;
    float mult = 1.0f, add = 0.0f, step = 0.0f, trim;
    int opt;
    while ((opt = getopt(argc, argv, "i:b:r:o:m:a:lt:s:Bh")) != -1) {
        switch (opt) {
            case 'i':
                if (!bwIsBigWig(optarg, NULL)) {
                    fprintf(stderr, "[E::adjust] Not a bigWig file (-i): '%s'\n", optarg);
                    return EXIT_FAILURE;
                }
                bw_in = bwOpen(optarg, NULL, "r");
                if (!bw_in) {
                    fprintf(stderr, "[E::adjust] Unable to open bigWig (-i): '%s'\n", optarg);
                    return EXIT_FAILURE;
                }
                break;
            case 'b':
                bedfn = optarg;
                break;
            case 'r':
                rRange = optarg;
                break;
            case 'o':
                outfn = optarg;
                break;
            case 'm':
                mult = atof(optarg);
                if (mult == 0 && errno == ERANGE) {
                    fprintf(stderr, "[E::adjust] Unable to parse '-m': %s\n", strerror(errno));
                    return EXIT_FAILURE;
                } else if (mult == 0) {
                    fprintf(stderr, "[E::adjust] '-m' must be nonzero\n");
                    return EXIT_FAILURE;
                }
                break;
            case 'a':
                add = atof(optarg);
                if (add == 0 && errno == ERANGE) {
                    fprintf(stderr, "[E::adjust] Unable to parse '-a': %s\n", strerror(errno));
                    return EXIT_FAILURE;
                }
                break;
            case 'l':
                do_log10 = true;
                break;
            case 't':
                use_trim = true;
                trim = atof(optarg);
                if (trim == 0 && errno == ERANGE) {
                    fprintf(stderr, "[E::adjust] Unable to parse '-t': %s\n", strerror(errno));
                    return EXIT_FAILURE;
                }
                break;
            case 's':
                step = atof(optarg);
                if (step == 0 && errno == ERANGE) {
                    fprintf(stderr, "[E::adjust] Unable to parse '-s': %s\n", strerror(errno));
                    return EXIT_FAILURE;
                }
                break;
            case 'B':
                bg = true;
                break;
            case 'h':
                help_adjust();
                return EXIT_SUCCESS;
            default:
                return EXIT_FAILURE;
        }
    }
    if (outfn != NULL) {
        if (bg) {
            fout.f = NULL;
            if (strcmp(outfn, "-") == 0) {
                tofile = false;
                fout.f = fdopen(1, "w");
                if (fout.f == NULL) {
                    fprintf(stderr, "[E::adjust] Failed to fdopen stdout: %s\n", strerror(errno));
                    return EXIT_FAILURE;
                }
            } else {
                fout.gz = gzopen(outfn, "wb");
                if (fout.gz == NULL) {
                    int e;
                    fprintf(stderr, "[E::adjust] Failed to open file '%s': %s\n", outfn, gzerror(fout.gz, &e));
                    return EXIT_FAILURE;
                }
            }
        } else {
            bw_out = bwOpen(outfn, NULL, "w");
            if (!bw_out) {
                fprintf(stderr, "[E::adjust] Unable to create output (-o) '%s'\n", outfn);
                return EXIT_FAILURE;
            }
        }
    } else {
        fprintf(stderr, "[E::adjust] Missing -o\n");
        return EXIT_FAILURE;
    }
    if (bw_in == NULL) {
        fprintf(stderr, "[E::adjust] Missing -i\n");
        return EXIT_FAILURE;
    }
    if (bedfn != NULL) {
        bed = readBED(bedfn, false, 0, 0, 0);
        if (bed == NULL) return EXIT_FAILURE;
        bed_index(bed->data);
        bed->n = bed_unify(bed->data);
    } else {
        bed = alloc(sizeof(bed_t));
        kh_bedHash_t *h = kh_init(bedHash);
        khint64_t k;
        int ret;
        if (rRange == NULL) {
            bed->n_names = bw_in->cl->nKeys;
            bed->names = alloc(sizeof(char *) * bed->n_names);
            for (int64_t i = 0; i < bw_in->cl->nKeys; i++) {
                bed->names[i] = bw_in->cl->chrom[i];
                k = kh_put(bedHash, h, bw_in->cl->chrom[i], &ret);
                if (ret == -1) {
                    fprintf(stderr, "[E::adjust] Failed to create chrom hash table\n");
                    return EXIT_FAILURE;
                }
                bedList_t *p = &kh_val(h, k);
                p->a = alloc(sizeof(p->a[0]));
                p->a[0].namei = i;
                p->a[0].beg = 0;
                p->a[0].end = bw_in->cl->len[i];
                p->n = 1;
            }
        } else {
            uint32_t beg, end;
            uint32_t tid = parseSingleRange(rRange, bw_in, &beg, &end);
            if (tid == -1) {
                return EXIT_FAILURE;
            }
            bed->n_names = 1;
            bed->names = alloc(sizeof(char *) * bed->n_names);
            bed->names[0] = bw_in->cl->chrom[tid];
            k = kh_put(bedHash, h, bw_in->cl->chrom[tid], &ret);
            if (ret == -1) {
                fprintf(stderr, "[E::adjust] Failed to create chrom hash table\n");
                return EXIT_FAILURE;
            }
            bedList_t *p = &kh_val(h, k);
            p->a = alloc(sizeof(p->a[0]));
            p->a[0].namei = tid;
            p->a[0].beg = beg;
            p->a[0].end = end;
            p->n = 1;
        }
        bed->data = h;
    }

    if (!bg) {
        if (bwCreateHdr(bw_out, 10)) {
            fprintf(stderr, "[E::adjust] Failed to init output bigWig header\n");
            return EXIT_FAILURE;
        }
        bw_out->cl = bwCreateChromList((const char * const *)bw_in->cl->chrom, bw_in->cl->len, bw_in->cl->nKeys);
        if (!bw_out->cl) {
            fprintf(stderr, "[E::adjust] Failed to create output bigWig chrom list\n");
            return EXIT_FAILURE;
        }
        if (bwWriteHdr(bw_out)) {
            fprintf(stderr, "[E::adjust] Failed to write output bigWig header\n");
            return EXIT_FAILURE;
        }
    }

    ranges_t *ranges = alloc(sizeof(ranges_t));
    ranges->m = RANGES_MEM;
    ranges->chrom = alloc(sizeof(char *) * ranges->m);
    ranges->start = alloc(sizeof(uint32_t) * ranges->m);
    ranges->end = alloc(sizeof(uint32_t) * ranges->m);
    ranges->val = alloc(sizeof(float) * ranges->m);

    uint32_t start, end, bed_start, bed_end;
    char *chromName;
    khint_t k;
    bedList_t *b;
    float val;
    for (int64_t i = 0; i < bw_in->cl->nKeys; i++) {
        chromName = bw_in->cl->chrom[i];
        k = kh_get(bedHash, bed->data, chromName);
        if (k != kh_end(bed->data)) {
            b = &kh_val(bed->data, k);
            ranges->n = ranges->times_added = 0;
            for (int j = 0; j < b->n; j++) {
                bed_start = b->a[j].beg;
                bed_end = b->a[j].end;
                bwOverlapIterator_t *bwIt = bwOverlappingIntervalsIterator(bw_in, chromName, bed_start, bed_end, BW_ITERATOR_CHUNKS);
                if (bwIt == NULL) {
                    fprintf(stderr, "[E::adjust] Error traversing bigWig\n");
                    return EXIT_FAILURE;
                }
                while (bwIt->data) {
                    for (int64_t h = 0; h < bwIt->intervals->l; h++) {
                        start = bwIt->intervals->start[h];
                        end = bwIt->intervals->end[h];
                        if (start < bed_start) start = bed_start;
                        if (end > bed_end) end = bed_end;
                        val = bwIt->intervals->value[h];
                        val += add;
                        val *= mult;
                        if (do_log10) val = log10f(val);
                        if (use_trim) val = min(trim, val);
                        if (step != 0.0f) val = roundf(val / step) * step;
                        if (ranges->n && ranges->end[ranges->n - 1] == start && ranges->val[ranges->n - 1] == val) {
                            ranges->end[ranges->n - 1] = end;
                        } else {
                            if (ranges->n == ranges->m) {
                                if (!bg) {
                                    addBwInterval(bw_out, ranges);
                                } else {
                                    addBgInterval(ranges, tofile);
                                }
                            }
                            ranges->chrom[ranges->n] = chromName;
                            ranges->start[ranges->n] = start;
                            ranges->end[ranges->n] = end;
                            ranges->val[ranges->n++] = val;
                        }
                    }
                    bwIt = bwIteratorNext(bwIt);
                    if (bwIt == NULL) {
                        fprintf(stderr, "[E::adjust] Error traversing bigWig\n");
                        return EXIT_FAILURE;
                    }
                }
                bwIteratorDestroy(bwIt);
            }
            if (ranges->n) {
                if (!bg) {
                    addBwInterval(bw_out, ranges);
                } else {
                    addBgInterval(ranges, tofile);
                }
            }
        }
    }

    if (!bg) {
        bwClose(bw_out);
    } else {
        if (tofile && gzclose (fout.gz) != Z_OK) {
            int e;
            fprintf(stderr, "[E::adjust] Failed to close output file: %s\n", gzerror(fout.gz, &e));
            return EXIT_FAILURE;
        } else if (!tofile && fclose(fout.f) != 0) {
            fprintf(stderr, "[E::adjust] Failed to close output stream: %s\n", strerror(errno));
            return EXIT_FAILURE;
        }
    }
    bwClose(bw_in);
    bwCleanup();
    return EXIT_SUCCESS;
}

// score ----------------------------------------------------------------------

static void help_score(void) {
    printf(
        "bwtk v%s  Copyright (C) %s  Benjamin Jean-Marie Tremblay\n"
        "bwtk score [options] -i <file.bw> -o <scores.tsv>\n"
        "    -i    Input bigWig\n"
        "    -o    Output scores in TSV format (use '-' for stdout)\n"
        "    -b    BED file to score, otherwise scores chromosomes ('-' for stdin)\n"
        "    -B    Return a BED instead with this stat in the score column\n"
        "    -h    Print this message and exit\n"
        , BWTK_VERSION, BWTK_YEAR
    );
}

enum BWSTAT {
    NONE,
    SIZE,
    COVERED,
    SUM,
    MEAN0,
    MEAN,
    MIN,
    MAX,
};

static int score(int argc, char **argv) {
    if (argc == 1) {
        fprintf(stderr, "[E::score] No arguments provided\n");
        help_score();
        return EXIT_FAILURE;
    }
    if (bwInit(1<<17) != 0) {
        fprintf(stderr, "[E::score] Unable to init curl buffer\n");
        return EXIT_FAILURE;
    }
    char *bedfn = NULL, *bedStatString = NULL;
    enum BWSTAT bedStat = NONE;
    bed_t *bed;
    bigWigFile_t *bw = NULL;
    FILE *fout = NULL;
    int opt;
    while ((opt = getopt(argc, argv, "i:b:o:B:h")) != -1) {
        switch (opt) {
            case 'i':
                if (!bwIsBigWig(optarg, NULL)) {
                    fprintf(stderr, "[E::score] Not a bigWig file (-i): '%s'\n", optarg);
                    return EXIT_FAILURE;
                }
                bw = bwOpen(optarg, NULL, "r");
                if (!bw) {
                    fprintf(stderr, "[E::score] Unable to open bigWig (-i): '%s'\n", optarg);
                    return EXIT_FAILURE;
                }
                break;
            case 'b':
                bedfn = optarg;
                break;
            case 'o':
                fout = strcmp(optarg, "-") ? fopen(optarg, "w") : fdopen(1, "w");
                if (fout == NULL) {
                    fprintf(stderr, "[E::score] Unable to open file (-o): %s\n", strerror(errno));
                    return EXIT_FAILURE;
                }
                break;
            case 'B':
                bedStatString = optarg;
                break;
            case 'h':
                help_score();
                return EXIT_SUCCESS;
            default:
                return EXIT_FAILURE;
        }
    }
    if (bedStatString != NULL) {
        if (strcmp(bedStatString, "size") == 0) {
            bedStat = SIZE;
        } else if (strcmp(bedStatString, "covered") == 0) {
            bedStat = COVERED;
        } else if (strcmp(bedStatString, "sum") == 0) {
            bedStat = SUM;
        } else if (strcmp(bedStatString, "mean0") == 0) {
            bedStat = MEAN0;
        } else if (strcmp(bedStatString, "mean") == 0) {
            bedStat = MEAN;
        } else if (strcmp(bedStatString, "min") == 0) {
            bedStat = MIN;
        } else if (strcmp(bedStatString, "max") == 0) {
            bedStat = MAX;
        } else {
            fprintf(stderr, "[E::score] Unknown stat '%s': expected one of\n", bedStatString);
            fprintf(stderr, "[E::score]   covered, sum, mean0, mean, min, max\n");
            return EXIT_FAILURE;
        }
    }
    if (bw == NULL) {
        fprintf(stderr, "[E::score] Missing -i\n");
        return EXIT_FAILURE;
    }
    if (bedfn != NULL) {
        bed = readBED(bedfn, false, 0, 0, 0);
        if (bed == NULL) return EXIT_FAILURE;
        bed_index(bed->data);
    } else {
        bed = alloc(sizeof(bed_t));
        bed->n_names = bw->cl->nKeys;
        bed->names = alloc(sizeof(char *) * bed->n_names);
        kh_bedHash_t *h = kh_init(bedHash);
        khint64_t k;
        int ret;
        for (int64_t i = 0; i < bw->cl->nKeys; i++) {
            bed->names[i] = bw->cl->chrom[i];
            k = kh_put(bedHash, h, bw->cl->chrom[i], &ret);
            if (ret == -1) {
                fprintf(stderr, "[E::score] Failed to create chrom hash table\n");
                return EXIT_FAILURE;
            }
            bedList_t *p = &kh_val(h, k);
            p->a = alloc(sizeof(p->a[0]));
            p->a[0].namei = i;
            p->a[0].beg = 0;
            p->a[0].end = bw->cl->len[i];
            p->n = 1;
        }
        bed->data = h;
    }
    if (fout == NULL) {
        fprintf(stderr, "[E::score] Missing -o\n");
        return EXIT_FAILURE;
    }

    if (bedStat == NONE) {
        fputs("name\tsize\tcovered\tsum\tmean0\tmean\tmin\tmax\n", fout);
    }

    int64_t nranges = 0;
    uint32_t bed_start, bed_end, start, end;
    char *chromName;
    khint_t k;
    bedList_t *b;
    bwOverlapIterator_t *bwVals;
    uint32_t size, covered;
    double sum, mean0, mean, min_val, max_val;
    bool first = false;
    for (int64_t i = 0; i < bw->cl->nKeys; i++) {
        chromName = bw->cl->chrom[i];
        k = kh_get(bedHash, bed->data, chromName);
        if (k != kh_end(bed->data)) {
            b = &kh_val(bed->data, k);
            for (int j = 0; j < b->n; j++) {
                nranges++;
                bed_start = b->a[j].beg;
                bed_end = b->a[j].end;
                size = bed_end - bed_start;
                bwVals = bwOverlappingIntervalsIterator(bw, chromName, bed_start, bed_end, BW_ITERATOR_CHUNKS);
                if (bwVals == NULL) {
                    fprintf(stderr, "[E::score] Error traversing bigWig\n");
                    return EXIT_FAILURE;
                }
                first = true;
                mean = mean0 = covered = sum = min_val = max_val = 0;
                while (bwVals->data) {
                    for (int64_t h = 0; h < bwVals->intervals->l; h++) {
                        if (first) {
                            max_val = min_val = (double) bwVals->intervals->value[h];
                            first = false;
                        }
                        min_val = min(min_val, (double) bwVals->intervals->value[h]);
                        max_val = max(max_val, (double) bwVals->intervals->value[h]);
                        start = max(bed_start, bwVals->intervals->start[h]);
                        end = min(bed_end, bwVals->intervals->end[h]);
                        covered += end - start;
                        sum += (double) bwVals->intervals->value[h] * (double) (end - start);
                    }
                    bwVals = bwIteratorNext(bwVals);
                    if (bwVals == NULL) {
                        fprintf(stderr, "[E::score] Error traversing bigWig\n");
                        return EXIT_FAILURE;
                    }
                }
                if (covered) {
                    mean = sum / (double) covered;
                    mean0 = sum / (double) size;
                }
                switch (bedStat) {
                    case NONE:
                        fprintf(fout, "%s\t%u\t%u\t%g\t%g\t%g\t%g\t%g\n", bed->names[b->a[j].namei], size, covered, sum, mean0, mean, min_val, max_val);
                        break;
                    case SIZE:
                        fprintf(fout, "%s\t%u\t%u\t%s\t%u\t%c\n", chromName, bed_start, bed_end, bed->names[b->a[j].namei], size, b->a[j].strand);
                        break;
                    case COVERED:
                        fprintf(fout, "%s\t%u\t%u\t%s\t%u\t%c\n", chromName, bed_start, bed_end, bed->names[b->a[j].namei], covered, b->a[j].strand);
                        break;
                    case SUM:
                        fprintf(fout, "%s\t%u\t%u\t%s\t%g\t%c\n", chromName, bed_start, bed_end, bed->names[b->a[j].namei], sum, b->a[j].strand);
                        break;
                    case MEAN0:
                        fprintf(fout, "%s\t%u\t%u\t%s\t%g\t%c\n", chromName, bed_start, bed_end, bed->names[b->a[j].namei], mean0, b->a[j].strand);
                        break;
                    case MEAN:
                        fprintf(fout, "%s\t%u\t%u\t%s\t%g\t%c\n", chromName, bed_start, bed_end, bed->names[b->a[j].namei], mean, b->a[j].strand);
                        break;
                    case MIN:
                        fprintf(fout, "%s\t%u\t%u\t%s\t%g\t%c\n", chromName, bed_start, bed_end, bed->names[b->a[j].namei], min_val, b->a[j].strand);
                        break;
                    case MAX:
                        fprintf(fout, "%s\t%u\t%u\t%s\t%g\t%c\n", chromName, bed_start, bed_end, bed->names[b->a[j].namei], max_val, b->a[j].strand);
                        break;
                }
                bwIteratorDestroy(bwVals);
            }
        }
    }
    if (!nranges) {
        fprintf(stderr, "[E::score] Found zero valid BED ranges in bigWig\n");
        return EXIT_FAILURE;
    }
    if (nranges < bed->n) {
        fprintf(stderr, "[W::score] Found only %lld/%lld ranges in bigWig (missing chroms)\n", nranges, bed->n);
    }

    if (fclose(fout) != 0) {
        fprintf(stderr, "[E::score] Failed to close output stream: %s\n", strerror(errno));
        return EXIT_FAILURE;
    }
    bwClose(bw);
    bwCleanup();
    return EXIT_SUCCESS;
}

// values ----------------------------------------------------------------------

static void help_values(void) {
    printf(
        "bwtk v%s  Copyright (C) %s  Benjamin Jean-Marie Tremblay\n"
        "bwtk values [options] -i <file.bw> -b <ranges.bed> -o <values.tsv>\n"
        "    -i    Input bigWig\n"
        "    -b    BED file with ranges to extract values from ('-' for stdin)\n"
        "    -o    Output values in TSV format (use '-' for stdout)\n"
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
    if (bwInit(1<<17) != 0) {
        fprintf(stderr, "[E::values] Unable to init curl buffer\n");
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
        bed = readBED(bedfn, true, (uint32_t) size, left, right);
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
                        fprintf(fout, "\t%g", v);
                    }
                } else {
                    for (uint32_t h = 0; h < bwVals->l; h++) {
                        v = isnan(bwVals->value[h]) ? 0.0 : (double) bwVals->value[h];
                        fprintf(fout, "\t%g", v);
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
        fprintf(stderr, "[W::values] Found only %lld/%lld ranges in bigWig (missing chroms)\n", nranges, bed->n);
    }

    if (fclose(fout) != 0) {
        fprintf(stderr, "[E::values] Failed to close output stream: %s\n", strerror(errno));
        return EXIT_FAILURE;
    }
    bwClose(bw);
    bwCleanup();
    return EXIT_SUCCESS;
}

// chromsizes ------------------------------------------------------------------

// TODO: Either add an option to sort by name/length or do it by default
static void help_chromsizes(void) {
    printf(
        "bwtk v%s  Copyright (C) %s  Benjamin Jean-Marie Tremblay\n"
        "bwtk chroms [options] -i <file.bw> -o <chrom.sizes>\n"
        "    -i    Input bigWig\n"
        "    -o    Output chrom.sizes file (use '-' for stdout)\n"
        "    -h    Print this message and exit\n"
        , BWTK_VERSION, BWTK_YEAR
    );
}

static int chromsizes(int argc, char **argv) {
    if (argc == 1) {
        fprintf(stderr, "[E::chroms] No arguments provided\n");
        help_chromsizes();
        return EXIT_FAILURE;
    }
    if (bwInit(1<<17) != 0) {
        fprintf(stderr, "[E::chroms] Unable to init curl buffer\n");
        return EXIT_FAILURE;
    }
    int opt;
    FILE *fout = NULL;
    bigWigFile_t *bw = NULL;
    while ((opt = getopt(argc, argv, "i:o:h")) != -1) {
        switch (opt) {
            case 'i':
                if (!bwIsBigWig(optarg, NULL)) {
                    fprintf(stderr, "[E::chroms] Not a bigWig file (-i): '%s'\n", optarg);
                    return EXIT_FAILURE;
                }
                bw = bwOpen(optarg, NULL, "r");
                if (!bw) {
                    fprintf(stderr, "[E::chroms] Unable to open bigWig (-i): '%s'\n", optarg);
                    return EXIT_FAILURE;
                }
                break;
            case 'o':
                fout = strcmp(optarg, "-") ? fopen(optarg, "w") : fdopen(1, "w");
                if (fout == NULL) {
                    fprintf(stderr, "[E::chroms] Unable to open file (-o): %s\n", strerror(errno));
                    return EXIT_FAILURE;
                }
                break;
                break;
            case 'h':
                help_chromsizes();
                return EXIT_SUCCESS;
            default:
                return EXIT_FAILURE;
        }
    }
    if (bw == NULL) {
        fprintf(stderr, "[E::chroms] Missing -i\n");
        return EXIT_FAILURE;
    }
    if (fout == NULL) {
        fprintf(stderr, "[E::chroms] Missing -o\n");
        return EXIT_FAILURE;
    }
    if (!bw->cl->nKeys) {
        fprintf(stderr, "[E::chroms] Zero chromosomes in bigWig\n");
        return EXIT_FAILURE;
    }
    for (int64_t i = 0; i < bw->cl->nKeys; i++) {
        fprintf(fout, "%s\t%"PRIu32"\n", bw->cl->chrom[i], bw->cl->len[i]);
    }
    bwClose(bw);
    fclose(fout);
    bwCleanup();
    return EXIT_SUCCESS;
}

// merge -----------------------------------------------------------------------

static void help_merge(void) {
        /* "    -c    Merge chromosomes in chunks to save memory\n" */
    printf(
        "bwtk v%s  Copyright (C) %s  Benjamin Jean-Marie Tremblay\n"
        "bwtk merge [options] -o <out.bw> <file1.bw> <file2.bw> [...]\n"
        "    -o    Output bigWig\n"
        "    -S    Sum values instead of averaging\n"
        "    -M    Take the max value instead of averaging\n"
        "    -a    Add this value to scores [0]\n"
        "    -m    Multiply scores by this value [1]\n"
        "    -l    log10-transform scores\n"
        "    -t    Trim values above this max [Inf]\n"
        "    -s    Step size for binning [0]\n"
        "    -h    Print this message and exit\n"
        "Order of operations: avg|sum|max -> a -> m -> l -> t -> s\n"
        , BWTK_VERSION, BWTK_YEAR
    );
}

#define BW_MERGE_CHUNK_SIZE 1000000

static int addConsensusChromsToBw(bigWigFile_t **bw_in, const int n_bw, bigWigFile_t *bw_out) {
    khash_t(str_m) *conChroms = kh_init(str_m);
    khint64_t k;
    int absent;
    bool warn_about_size_diffs = false;
    for (int i = 0; i < n_bw; i++) {
        for (int j = 0; j < bw_in[i]->cl->nKeys; j++) {
            k = kh_get(str_m, conChroms, bw_in[i]->cl->chrom[j]);
            if (k == kh_end(conChroms)) {
                k = kh_put(str_m, conChroms, bw_in[i]->cl->chrom[j], &absent);
                if (absent == -1) {
                    fprintf(stderr, "[E::addConsensusChroms] Failed to create conChroms hash table\n");
                    return 1;
                }
                kh_val(conChroms, k).x = 0;
                kh_val(conChroms, k).y = (int64_t) bw_in[i]->cl->len[j];
            } else {
                if (kh_val(conChroms, k).y != (int64_t) bw_in[i]->cl->len[j]) {
                    warn_about_size_diffs = true;
                    kh_val(conChroms, k).y = min(kh_val(conChroms, k).y, (int64_t) bw_in[i]->cl->len[j]);
                }
            }
            kh_val(conChroms, k).x++;
        }
    }
    if (warn_about_size_diffs) {
        fprintf(stderr, "[W::addConsensusChroms] Found identical chromosomes across files with differing lengths\n");
    }
    int total_chroms = kh_size(conChroms);
    int final_chroms = 0;
    for (k = kh_begin(conChroms); k!= kh_end(conChroms); k++) {
        if (kh_exist(conChroms, k)) {
            int64_t n_chrom = kh_val(conChroms, k).x;
            if (n_chrom == (int64_t) n_bw) {
                final_chroms++;
            }
        }
    }
    if (!final_chroms) {
        fprintf(stderr, "[E::addConsensusChroms] Found %d possible chromosomes, but none are present in all bigWigs\n", total_chroms);
        return EXIT_FAILURE;
    }
    if (final_chroms < total_chroms) {
        fprintf(stderr, "[W::addConsensusChroms] Found %d possible chromosomes, but only %d are present in all bigWigs\n", total_chroms, final_chroms);
    }
    char **chromNames = alloc(sizeof(char *) * final_chroms);
    uint32_t *chromLens = alloc(sizeof(uint32_t) * final_chroms);
    int chrom_i = 0;
    for (k = kh_begin(conChroms); k!= kh_end(conChroms); k++) {
        if (kh_exist(conChroms, k) && kh_val(conChroms, k).x == (int64_t) n_bw) {
            chromNames[chrom_i] = strdup(kh_key(conChroms, k));
            chromLens[chrom_i++] = (uint32_t) kh_val(conChroms, k).y;
        }
    }
    bw_out->cl = bwCreateChromList((const char * const *)chromNames, chromLens, final_chroms);
    if (!bw_out->cl) {
        fprintf(stderr, "[E::merge] Failed to create output bigWig chrom list\n");
        return EXIT_FAILURE;
    }
    kh_destroy(str_m, conChroms);
    free(chromLens);
    for (int i = 0; i < final_chroms; i++) {
        free(chromNames[i]);
    }
    free(chromNames);
    return 0;
}

static int merge(int argc, char **argv) {
    if (argc == 1) {
        fprintf(stderr, "[E::merge] No arguments provided\n");
        help_merge();
        return EXIT_FAILURE;
    }
    if (bwInit(1<<17) != 0) {
        fprintf(stderr, "[E::merge] Unable to init curl buffer\n");
        return EXIT_FAILURE;
    }
    int opt, n_bw;
    bigWigFile_t *bw_out = NULL;
    bigWigFile_t **bw_in = NULL;
    float mult = 1.0f, add = 0.0f, step = 1.0f, trim;
    bool do_log10 = false, use_trim = false;
    bool sum_only = false, take_max = false;
    /* bool chunk = false; */
    while ((opt = getopt(argc, argv, "o:SMa:m:lt:s:h")) != -1) {
        switch (opt) {
            case 'o':
                bw_out = bwOpen(optarg, NULL, "w");
                if (!bw_out) {
                    fprintf(stderr, "[E::merge] Unable to create output (-o) '%s'\n", optarg);
                    return EXIT_FAILURE;
                }
                break;
            case 'S':
                sum_only = true;
                break;
            case 'M':
                take_max = true;
                break;
            case 'm':
                mult = atof(optarg);
                if (mult == 0 && errno == ERANGE) {
                    fprintf(stderr, "[E::merge] Unable to parse '-m': %s\n", strerror(errno));
                    return EXIT_FAILURE;
                } else if (mult == 0) {
                    fprintf(stderr, "[E::merge] '-m' must be nonzero\n");
                    return EXIT_FAILURE;
                }
                break;
            case 'a':
                add = atof(optarg);
                if (add == 0 && errno == ERANGE) {
                    fprintf(stderr, "[E::merge] Unable to parse '-a': %s\n", strerror(errno));
                    return EXIT_FAILURE;
                }
                break;
            case 'l':
                do_log10 = true;
                break;
            case 't':
                use_trim = true;
                trim = atof(optarg);
                if (trim == 0 && errno == ERANGE) {
                    fprintf(stderr, "[E::merge] Unable to parse '-t': %s\n", strerror(errno));
                    return EXIT_FAILURE;
                }
                break;
            case 's':
                step = atof(optarg);
                if (step == 0 && errno == ERANGE) {
                    fprintf(stderr, "[E::merge] Unable to parse '-s': %s\n", strerror(errno));
                    return EXIT_FAILURE;
                }
                break;
            case 'h':
                help_merge();
                return EXIT_SUCCESS;
            default:
                return EXIT_FAILURE;
        }
    }
    if (sum_only && take_max) {
        fprintf(stderr, "[E::merge] Cannot set both -s and -M\n");
        return EXIT_FAILURE;
    }
    if (bw_out == NULL) {
        fprintf(stderr, "[E::merge] Missing output (-o)\n");
        return EXIT_FAILURE;
    }
    if (optind == argc) {
        fprintf(stderr, "[E::merge] Missing input bigWigs\n");
        return EXIT_FAILURE;
    }
    n_bw = argc - optind;
    if (n_bw < 2) {
        fprintf(stderr, "[E::merge] Expected at least two input bigWigs\n");
        return EXIT_FAILURE;
    }
    bw_in = alloc(sizeof(bigWigFile_t) * n_bw);
    for (int i = 0; i < n_bw; i++) {
        if (!bwIsBigWig(argv[optind + i], NULL)) {
            fprintf(stderr, "[E::merge] File is not a bigWig: '%s'\n", argv[optind + i]);
            return EXIT_FAILURE;
        }
        bw_in[i] = bwOpen(argv[optind + i], NULL, "r");
        if (!bw_in[i]) {
            fprintf(stderr, "[E::merge] Unable to open bigWig: '%s'\n", argv[optind + i]);
            return EXIT_FAILURE;
        }
    }
    if (bwCreateHdr(bw_out, 10)) {
        fprintf(stderr, "[E::merge] Failed to init output bigWig header\n");
        return EXIT_FAILURE;
    }
    if (addConsensusChromsToBw(bw_in, n_bw, bw_out)) {
        return EXIT_FAILURE;
    }
    if (bwWriteHdr(bw_out)) {
        fprintf(stderr, "[E::merge] Failed to write output bigWig header\n");
        return EXIT_FAILURE;
    }

    ranges_t *ranges = alloc(sizeof(ranges_t));
    ranges->m = RANGES_MEM;
    ranges->chrom = alloc(sizeof(char *) * ranges->m);
    ranges->start = alloc(sizeof(uint32_t) * ranges->m);
    ranges->end = alloc(sizeof(uint32_t) * ranges->m);
    ranges->val = alloc(sizeof(float) * ranges->m);

    bwOverlapIterator_t *bwIt;
    uint32_t max_len = bw_out->cl->len[0];
    for (int i = 1; i < bw_out->cl->nKeys; i++) {
        max_len = max(max_len, bw_out->cl->len[i]);
    }
    float *values = alloc(sizeof(float) * max_len);
    for (int i = 0; i < bw_out->cl->nKeys; i++) {
        memset(values, 0, bw_out->cl->len[i] * sizeof(float));
        for (int j = 0; j < n_bw; j++) {
            bwIt = bwOverlappingIntervalsIterator(bw_in[j], bw_out->cl->chrom[i], 0, bw_out->cl->len[i], BW_ITERATOR_CHUNKS);
            if (bwIt == NULL) {
                fprintf(stderr, "[E::merge] Error traversing bigWig\n");
                return EXIT_FAILURE;
            }
            while (bwIt->data) {
                for (int64_t h = 0; h < bwIt->intervals->l; h++) {
                    for (uint32_t m = bwIt->intervals->start[h]; m < bwIt->intervals->end[h]; m++) {
                        if (take_max) {
                            values[m] = max(values[m], bwIt->intervals->value[h]);
                        } else if (sum_only) {
                            values[m] += bwIt->intervals->value[h];
                        } else {
                            values[m] += bwIt->intervals->value[h] / (float) n_bw;
                        }
                    }
                }
                bwIt = bwIteratorNext(bwIt);
                if (bwIt == NULL) {
                    fprintf(stderr, "[E::merge] Error traversing bigWig\n");
                    return EXIT_FAILURE;
                }
            }
            bwIteratorDestroy(bwIt);
        }
        ranges->times_added = 0;
        for (uint32_t j = 0; j < bw_out->cl->len[i]; j++) {
            values[j] += add;
            values[j] *= mult;
            if (do_log10) values[j] = log10f(values[j]);
            if (use_trim) values[j] = min(trim, values[j]);
            if (step != 0) values[j] = roundf(values[j] / step) * step;
        }
        for (uint32_t j = 0; j < bw_out->cl->len[i]; j++) {
            if (ranges->n == ranges->m) {
                addBwInterval(bw_out, ranges);
            }
            ranges->chrom[ranges->n] = bw_out->cl->chrom[i];
            ranges->start[ranges->n] = j;
            ranges->val[ranges->n] = values[j];
            while (j + 1 < bw_out->cl->len[i] && values[j + 1] == values[j]) j++;
            ranges->end[ranges->n++] = j + 1;
        }
        if (ranges->n) {
            addBwInterval(bw_out, ranges);
        }
    }
    free(values);
    /*  TODO: implement chunking (instead of alloc()'ing whole chromosome)
    } else {
        fprintf(stderr, "[E::merge] -c not yet implemented.\n");
        return EXIT_FAILURE;
    }
    */

    for (int i = 0; i < n_bw; i++) {
        bwClose(bw_in[i]);
    }
    bwClose(bw_out);
    bwCleanup();
    return EXIT_SUCCESS;

}

// bg2bw -----------------------------------------------------------------------

#define CHROMSIZES_MEM   128
#define CHROMNAME_SIZE  1024
#define PRESET_GENOMES  "tair10"

static void help_bg2bw(void) {
    printf(
        "bwtk v%s  Copyright (C) %s  Benjamin Jean-Marie Tremblay\n"
        "bwtk bg2bw [options] -g <chrom.sizes> -i <file.bedGraph[.gz]> -o <file.bw>\n"
        "    -i    Input bedGraph (can be gzipped), use '-' for stdin\n"
        "    -o    Output bigWig\n"
        "    -g    Genome chrom.sizes or genome.fa.fai file\n"
        "    -p    Use a preset genome instead of -g [%s]\n"
        "    -u    When using -p, use UCSC-style names (default: Ensembl)\n"
        "    -S    Ignore chromosomes found in bedGraph but not chrom.sizes\n"
        "    -a    Add this value to scores [0]\n"
        "    -m    Multiply scores by this value [1]\n"
        "    -l    log10-transform scores\n"
        "    -t    Trim values above this max [Inf]\n"
        "    -s    Step size for binning [0]\n"
        "    -h    Print this message and exit\n"
        "Order of operations: a -> m -> l -> t -> s\n"
        , BWTK_VERSION, BWTK_YEAR, PRESET_GENOMES
    );
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
        if (*ref == 0 || *ref == '#') continue; // TODO: Maybe some sequences can start with #?
        char chr[CHROMNAME_SIZE];
        uint32_t size;
        int num = sscanf(ref, "%s %"SCNu32, chr, &size);
        if (num < 2) {
            fprintf(stderr, "[E::readChromSizes] Expected at least two columns in chrom.sizes (-g)\n");
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
        fprintf(stderr, "[E::bg2bw] No arguments provided\n");
        help_bg2bw();
        return EXIT_FAILURE;
    }
    int opt;
    gzFile bg = NULL, cs = NULL;
    bigWigFile_t *bw = NULL;
    char *preset = NULL;
    bool ucsc = false;
    bool ignore_unknown_chr = false;
    bool do_log10 = false, use_trim = false;
    float mult = 1.0f, add = 0.0f, step = 0.0f, trim;
    while ((opt = getopt(argc, argv, "i:o:g:p:uSm:a:lt:s:h")) != -1) {
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
            case 'p':
                preset = optarg;
                break;
            case 'u':
                ucsc = true;
                break;
            case 'S':
                ignore_unknown_chr = true;
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
            case 't':
                use_trim = true;
                trim = atof(optarg);
                if (trim == 0 && errno == ERANGE) {
                    fprintf(stderr, "[E::bg2bw] Unable to parse '-t': %s\n", strerror(errno));
                    return EXIT_FAILURE;
                }
                break;
            case 's':
                step = atof(optarg);
                if (step == 0 && errno == ERANGE) {
                    fprintf(stderr, "[E::bg2bw] Unable to parse '-s': %s\n", strerror(errno));
                    return EXIT_FAILURE;
                }
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
        kh_val(bgChroms, k).x = -1;
        kh_val(bgChroms, k).y = (int64_t) chromSizes->sizes[i];
    }

    kstream_t *ks_bg = ks_init(bg);
    char *chrom_last = NULL;
#ifdef DEBUG
    int64_t nranges = 0;
#endif
    int64_t frow = 0, erow = 0, bwdumps = 0;
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
        k = kh_get(str_m, bgChroms, chr);
        if (k == kh_end(bgChroms)) {
            if (!ignore_unknown_chr) {
                fprintf(stderr, "[E::bg2bw] Found unknown chrom '%s' on line %lld\n", chr, frow);
                return EXIT_FAILURE;
            } else {
                continue;
            }
        }
        erow++;
        val += add;
        val *= mult;
        if (do_log10) val = log10f(val);
        if (use_trim) val = min(trim, val);
        if (step != 0.0f) val = roundf(val / step) * step;
        if (kh_val(bgChroms, k).x != -1 && kh_val(bgChroms, k).x < (erow - 1)) {
            fprintf(stderr, "[E::bg2bw] bedGraph file must be sorted (found out of order chroms)\n");
            return EXIT_FAILURE;
        }
        kh_val(bgChroms, k).x = erow;
        if (end <= start || start < 0) {
            fprintf(stderr, "[E::bg2bw] Bad bedGraph range on line %lld (start=%u end=%u)\n", frow, start, end);
            return EXIT_FAILURE;
        }
        if (end > kh_val(bgChroms, k).y) {
            fprintf(stderr, "[E::bg2bw] bedGraph range extends beyond chromosome limit on line %lld\n", frow);
            fprintf(stderr, "[E::bg2bw] Range: %s:%u-%u\n", chr, start, end);
            fprintf(stderr, "[E::bg2bw] Chromosome size: %lld\n", kh_val(bgChroms, k).y);
            return EXIT_FAILURE;
        }
        if (chrom_last == NULL || strcmp(chrom_last, chr) == 0) {
            if (end_last > 0 && start < end_last) {
                fprintf(stderr, "[E::bg2bw] bedGraph must be sorted and/or not contain overlapping ranges (lines %lld-%lld)\n", frow - 1, frow);
                return EXIT_FAILURE;
            }
            if (ranges->n && ranges->end[ranges->n - 1] == start && ranges->val[ranges->n - 1] == val) {
                ranges->end[ranges->n - 1] = end;
            } else {
                if (ranges->n == ranges->m) {
                    addBwInterval(bw, ranges);
                    bwdumps++;
                }
                ranges->chrom[ranges->n] = (char *) kh_key(bgChroms, k);
                ranges->start[ranges->n] = start;
                ranges->end[ranges->n] = end;
                ranges->val[ranges->n++] = val;
            }
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
#ifdef DEBUG
        nranges++;
#endif
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

#ifdef DEBUG
    fprintf(stderr, "Read %lld ranges, dumped to bw %lld times.\n", nranges, bwdumps);
    fprintf(stderr, "Peak mem before indexing: %'.2f MB\n", ((double) peak_mem() / 1024.0) / 1024.0);
#endif

    // TODO: Inside bwClose() is bwFinalize() -> constructZoomLevels(), which for a full
    // bedGraph increases the memory usage by 100x
    bwClose(bw);
#ifdef DEBUG
    fprintf(stderr, "Peak mem after indexing: %'.2f MB\n", ((double) peak_mem() / 1024.0) / 1024.0);
#endif

    return EXIT_SUCCESS;
}

// main ------------------------------------------------------------------------

static void help(void) {
    printf(
        "bwtk v%s  Copyright (C) %s  Benjamin Jean-Marie Tremblay\n"
        "Usage:  bwtk <subcommand> [options]\n"
        "Available subcommands:\n"
        "    bg2bw      Convert a bedGraph file to bigWig\n"
        "    adjust     Perform an operation on a bigWig\n"
        "    merge      Average multiple bigWig files together\n"
        "    values     Return bigWig values from overlapping BED ranges\n"
        "    score      Get summary scores of bigWig values from BED ranges\n"
        "    chroms     Print a chrom.sizes file from a bigWig header\n"
        "    help       Print this message and exit\n"
        "    version    Print the version number and exit\n"
        "For subcommand usage, try: bwtk <subcommand> -h\n"
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
    } else if (strcmp(argv[1], "bg2bw") == 0) {
        return bg2bw(argc - 1, argv + 1);
    } else if (strcmp(argv[1], "adjust") == 0) {
        return adjust(argc - 1, argv + 1);
    } else if (strcmp(argv[1], "merge") == 0) {
        return merge(argc - 1, argv + 1);
    } else if (strcmp(argv[1], "values") == 0) {
        return values(argc - 1, argv + 1);
    } else if (strcmp(argv[1], "chroms") == 0) {
        return chromsizes(argc - 1, argv + 1);
    } else if (strcmp(argv[1], "score") == 0) {
        return score(argc - 1, argv + 1);
    } else {
        fprintf(stderr, "[E::main] Unknown subcommand '%s'; try 'help' for usage\n", argv[1]);
        return EXIT_FAILURE;
    }

    return EXIT_SUCCESS;
}

