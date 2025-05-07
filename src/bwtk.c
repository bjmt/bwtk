#include <stdlib.h>
#include <stdio.h>
#include "libBigWig/bigWig.h"

int main(int argc, char **argv) {

    char *bg_fn = argv[1];
    char *bw_fn = argv[2];

    bigWigFile_t *bw = NULL;

    if (bwInit(1<<17) != 0) {
        fprintf(stderr, "Error loading libBigWig\n");
        return EXIT_FAILURE;
    }

    bw = bwOpen("test.bw", NULL, "w");
    if (!bw) {
        fprintf(stderr, "Error creating bigWig\n");
        return EXIT_FAILURE;
    }

    if (bwCreateHdr(bw, 10)) {
        fprintf(stderr, "Error creating bigWig header\n");
        return EXIT_FAILURE;
    }

    const char* const chr[] = { "1", "2", "Pt" };
    const uint32_t l[] = { 500, 400, 100 };

    bw->cl = bwCreateChromList(chr, l, 3);
    if (!bw->cl) {
        fprintf(stderr, "Error creating chrom list\n");
        return EXIT_FAILURE;
    }

    if (bwWriteHdr(bw)) {
        fprintf(stderr, "Error writing bigWig header\n");
        return EXIT_FAILURE;
    }

    const uint32_t start[] = { 50, 200, 20 };
    const uint32_t end[] = { 100, 275, 80 };
    const float values[] = { 25.0f, 86.0f, 0.5f };
    bwAddIntervals(bw, chr, start, end, values, 3);

    bwClose(bw);
    bwCleanup();

    return EXIT_SUCCESS;

}
