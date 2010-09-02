#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "imgio.h"
#include "FFT2D.h"
#include "imgfilter.h"

#define FILENAME_LENGTH 100

int main(int argc, char **argv) {

    FILE *a = fopen(argv[1], "r");
    FILE *b = fopen(argv[2], "r");
    struct portImage *pa = readImage(a);
    struct portImage *pb = readImage(b);

    printf("Root-Mean-Squared Difference: %.2f\n", imageCompare(pa, pb));
    printf("Different Pixels:             %d\n", imageDifferentPixels(pa, pb));
    printf("50%% Low Pass Difference       %.2f\n", fftSimilarity(pa, pb, 0.5));

    freeImage(pa);
    freeImage(pb);
    fclose(a);
    fclose(b);
    return 0;
} // End Main
