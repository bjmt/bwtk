#include <stdlib.h>
#include <stdint.h>
#include <stdio.h>
#include <string.h>
#include <stdbool.h>
#include <errno.h>
#include "libBigWig/bigWig.h"

#include "kseq.h"
KSTREAM_INIT(gzFile, gzread, 8192)

#include "khash.h"
KHASH_MAP_INIT_STR(str_m, int64_t)

#define BWTK_VERSION "1.0"
#define BWTK_YEAR "2025"

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

// Subcommands:
// - merge -> average some bw together 
// - bw2bg
// - bg2bw
// - values -> get 1-bp resolution values over ranges in a bed file (unsorted BED ok)
// - scale -> multiply values by a constant factor (also as an option for bg2bw)
// - version

static void help_bg2bw(void) {
    // remember to check that the bg is sorted and doesn't have overlapping ranges
        //"    -1    Create the bigWig in one pass (i.e. store it all in memory)\n"
    // Other options: addition, -log10, log2 (change scaling to mult?) [have an order of operations: add/mult->log]
    printf(
        "bwtk v%s  Copyright (C) %s  Benjamin Jean-Marie Tremblay\n"
        "bwtk bg2bw [options] -g <chrom.sizes> -i <file.bedGraph> -o <file.bw>\n"
        "    -i    Input bedGraph (can be gzipped), use '-' for stdin\n"
        "    -o    Output bigWig\n"
        "    -g    Genome chrom.sizes file\n"
        "    -s    Scaling factor [1]\n"
        "    -p    Use a preset genome\n"
        "    -h    Print this message and exit\n"
        , BWTK_VERSION, BWTK_YEAR
    );
}

static void help(void) {
    printf(
        "bwtk v%s  Copyright (C) %s  Benjamin Jean-Marie Tremblay\n"
        "Usage:  bwtk <subcommand> [options]\n"
        "Available subcommands:\n"
        "    merge      Combine multiple bigWig files\n"
        "    bw2bg      Convert a bigWig file to bedGraph\n"
        "    bg2bw      Convert a bedGraph file to bigWig\n"
        "    values     Return bigWig values from overlapping ranges\n"
        "    scale      Perform an operation on bigWig values\n"
        "    help       Print this message and exit\n"
        "    version    Print the version number and exit\n"
        , BWTK_VERSION, BWTK_YEAR
    );
}

void *calloc_or_die(size_t size, const char *func_name) {
  void *result = calloc(1, size);
  if (result == NULL) {
    fprintf(stderr, "[E::%s] Out of memory (requested %lu B).\n",
      func_name, size);
    exit(EXIT_FAILURE);
  }
  return result;
}
#define alloc(size) calloc_or_die((size), __func__)

/*
typedef struct {
    int64_t beg, end;
    int chr_i;
    char strand;
} bed_entry_t;

typedef struct {
    // add a hash table for chr names
    int64_t n, m;
    
    bed_entry_t **b;
} bed_t;

static inline void coord_resize(int64_t *s1, int64_t *s2, const int64_t up, const int64_t down, const char s) {
    int64_t len = *s2 - *s1;
    *s2 -= len / 2;
    *s1 += len / 2;
    if (s == '-') {
        *s1 += len % 2;
        *s1 = down > *s1 ? 0: *s1 - down;
        *s2 += up;
    } else {
        *s2 -= len % 2;
        *s1 = up > *s1 ? 0 : *s1 - up;
        *s2 += down;
    }
}

static void resize_bed(bed_t *bed, const int64_t size) {
    int64_t up = size / 2, down = size;
    down = down / 2 + down % 2 - 1;
    for (int64_t i = 0; i < bed->n; i++) {
        coord_resize(&bed->b[i]->beg, &bed->b[i]->end, up, down, bed->b[i]->strand);
    }
}

// Code modified from HTSlib
static bed_t *bed_read(const char *fn, const int64_t size)
{
    bed_t *bed = alloc(sizeof(bed_t));
    bed->b = alloc(sizeof(bed_entry_t) * 256);
    bed->m = bed->n = 256;
    gzFile fp;
    kstream_t *ks = NULL;
    int dret;
    unsigned int line = 0;
    kstring_t str = { 0, 0, NULL };

    if (NULL == bed) return NULL;
    // fp = strcmp(fn, "-")? gzopen(fn, "r") : gzdopen(fileno(stdin), "r");
    fp = gzopen(fn, "r");
    if (fp == 0) return NULL;
    ks = ks_init(fp);
    if (NULL == ks) goto fail;
    int ks_len;
    int64_t range_len;
    bool needs_resize = false;
    int64_t n = 0;
    while ((ks_len = ks_getuntil(ks, KS_SEP_LINE, &str, &dret)) >= 0) {
        char *ref = str.s, *ref_end;
        char strand;
        uint64_t beg = 0, end = 0;
        int num = 0;

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
            num = sscanf(ref_end + 1, "%"SCNu64" %"SCNu64" %*s %*s %c",
                         &beg, &end, &strand);
        }
        if (1 == num) {
            end = beg--;
        }
        if (num < 1 || end < beg) {
            if (0 == strcmp(ref, "browser")) continue;
            if (0 == strcmp(ref, "track")) continue;
            if (num < 1) {
                fprintf(stderr,
                        "[E::bed_read] Parse error reading \"%s\" at line %u\n",
                        fn, line);
            } else {
                fprintf(stderr,
                        "[E::bed_read] Parse error reading \"%s\" at line %u : "
                        "end (%"PRIu64") must not be less "
                        "than start (%"PRIu64")\n",
                        fn, line, end, beg);
            }
            goto fail;
        }
        if (num < 3) strand = '.';
        if (strand != '+' && strand != '-' && strand != '.') {
            fprintf(stderr, "[E::bed_read] Strand must be blank or one of +/-/.\n");
            goto fail;
        }

        if (n > bed->m) {
            bed->m *= 2;
            bed->b = realloc(bed->b, bed->m * sizeof(bed_entry_t));
            if (bed->b == NULL) {
                fprintf(stderr, "[E::bed_read] Out of memory.\n");
                goto fail;
            }
        }
        bed->n++;
        bed->b[n]->beg = beg;
        bed->b[n]->end = end;
        bed->b[n]->strand = strand;
        bed->b[n]->chr = strdup(ref);
        if (size) {
            if (end - beg != size) {
                needs_resize = true;
            }
        } else {
            if (!n) {
                range_len = end - beg;
            } else if (end - beg != range_len) {
                needs_resize = true;
                range_len = range_len < end - beg ? end - beg : range_len;
            }
        }
        n++;

    }

    if (gzclose(fp) != Z_OK) {
        fp = NULL;
        goto fail;
    }
    ks_destroy(ks);
    free(str.s);
    if (needs_resize) {
        resize_bed(bed, size ? size : range_len);
    }
    return bed;
 fail:
    if (ks) ks_destroy(ks);
    if (fp) gzclose(fp);
    free(str.s);
    return NULL;
}
*/

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

#define CHROMSIZES_MEM   128
#define CHROMNAME_SIZE  1024
#define RANGES_MEM      8192
static int bg2bw(int argc, char **argv) {
    int opt;
    gzFile bg = NULL, cs = NULL;
    bigWigFile_t *bw = NULL;
    float scale = 1.0f;
    while ((opt = getopt(argc, argv, "i:o:g:s:p:h")) != -1) {
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
            case 's':
                scale = atof(optarg);
                if (scale == 0 && errno == ERANGE) {
                    fprintf(stderr, "[E::bg2bw] Unable to parse '-s': %s\n", strerror(errno));
                    return EXIT_FAILURE;
                } else if (scale == 0) {
                    fprintf(stderr, "[E::bg2bw] '-s' must be nonzero\n");
                    return EXIT_FAILURE;
                }
                break;
            case 'p':
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
    if (cs == NULL) {
        fprintf(stderr, "[E::bg2bw] Missing chrom.sizes (-g)\n");
        return EXIT_FAILURE;
    }

    chromSizes_t *chromSizes = alloc(sizeof(chromSizes_t));
    chromSizes->m = CHROMSIZES_MEM;
    chromSizes->chroms = alloc(sizeof(char *) * chromSizes->m);
    chromSizes->sizes = alloc(sizeof(uint32_t) * chromSizes->m);

    ranges_t *ranges = alloc(sizeof(ranges_t));
    ranges->m = RANGES_MEM;
    ranges->chrom = alloc(sizeof(char *) * ranges->m);
    ranges->start = alloc(sizeof(uint32_t) * ranges->m);
    ranges->end = alloc(sizeof(uint32_t) * ranges->m);
    ranges->val = alloc(sizeof(float) * ranges->m);

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
            fprintf(stderr, "[E::bg2bw] Error: Expected at least two columns in chrom.sizes (-g)\n");
            return EXIT_FAILURE;
        }
        if (chromSizes->n + 1 > chromSizes->m) {
            chromSizes->m *= 2;
            chromSizes->chroms = realloc(chromSizes, sizeof(char *) * chromSizes->m);
            chromSizes->sizes = realloc(chromSizes, sizeof(uint32_t) * chromSizes->m);
            if (chromSizes->chroms == NULL || chromSizes->sizes == NULL) {
                fprintf(stderr, "[E::bg2bw] Out of memory\n");
            }
        }
        chromSizes->chroms[chromSizes->n] = strdup(chr);
        chromSizes->sizes[chromSizes->n++] = size;
        cs_lines++;
    }
    if (!cs_lines) {
        fprintf(stderr, "[E::bg2bw] Empty chrom.sizes file (-g)\n");
        return EXIT_FAILURE;
    }
    ks_destroy(ks_cs);

    if (bwCreateHdr(bw, 10)) {
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
        val *= scale;
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
        if (end <= start) {
            fprintf(stderr, "[E::bg2bw] Bad bedGraph range on line %lld (start=%u end=%u)\n", frow, start, end);
            return EXIT_FAILURE;
        }
        // Check start/end are within [0, size of chrom]?
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

    // TODO: Inside bwClose() is bwFinalize() -> writeIndex(), which for a full
    // bedGraph increases the memory usage by 100x
    bwClose(bw);
    gzclose(bg);
    gzclose(cs);

    fprintf(stderr, "Peak mem after indexing: %'.2f MB\n", ((double) peak_mem() / 1024.0) / 1024.0);

    return EXIT_SUCCESS;
}

int main(int argc, char **argv) {

    if (argc < 2) {
        fprintf(stderr, "Error: Missing subcommand.\n");
        help();
        return EXIT_FAILURE;
    }
    if (strcmp(argv[1], "version") == 0) {
        fprintf(stderr, "v%s\n", BWTK_VERSION);
        return EXIT_SUCCESS;
    } else if (strcmp(argv[1], "help") == 0) {
        help();
        return EXIT_SUCCESS;
    } else if (strcmp(argv[1], "merge") == 0) {
    } else if (strcmp(argv[1], "bw2bg") == 0) {
    } else if (strcmp(argv[1], "bg2bw") == 0) {
        return bg2bw(argc - 1, argv + 1);
    } else if (strcmp(argv[1], "values") == 0) {
    } else if (strcmp(argv[1], "scale") == 0) {
    } else {
        fprintf(stderr, "Error: Unknown subcommand '%s'. Run 'bwtk help' for usage.\n", argv[1]);
        return EXIT_FAILURE;
    }

    return EXIT_SUCCESS;
}

