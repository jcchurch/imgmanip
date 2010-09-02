#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "portable.h"

#define FILENAME_LENGTH 100

int main(int argc, char **argv) {

    int i;
    int window;
    int threshold;
    FILE *a;
    double noise[] = {0.0,0.05,0.10,0.15,0.20,0.25,0.30,0.35,0.40,0.45,0.50};
    double *scores;
    portImage *pa;
    a = fopen(argv[1], "r");
    sscanf(argv[2], "%d", &window);
    pa = readImage(a);

    scores = (double*)malloc( sizeof(double) * window * window );
    if (!scores) {
        fprintf(stderr, "Memory allocation failed.\n");
        return 1;
    }

    for (i = 0; i < 11; i++) { 
        printf("%s,%2.2f,", argv[1], noise[i]);
        calculateAverageDifference(pa, window, noise[i], scores);
        for (threshold = 0; threshold < window*window; threshold++)
            printf("%.3f,",scores[threshold]);
        printf("\n");
    }

    free(scores);
    freeImage(pa);
    fclose(a);
    return 0;
} // End Main
