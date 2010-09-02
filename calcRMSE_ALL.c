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
    portImage *na;
    a = fopen(argv[1], "r");
    sscanf(argv[2], "%d", &window);
    pa = readImage(a);

    scores = (double*)malloc( sizeof(double) * window * window );
    if (!scores) {
        fprintf(stderr, "Memory allocation failed.\n");
        return 1;
    }

    for (i = 0; i < 11; i++) { 
        na = copyImage(pa);
        addSaltPepperNoise(na, noise[i]);
        printf("%s,noise=%2.2f,", argv[1], noise[i]);
        calculatePotentialRMSE(pa, na, window, scores);

        for (threshold = 0; threshold < window*window; threshold++)
            printf("smf(%d)=%.3f,",threshold+1,scores[threshold]);

        calculatePotentialRMSE_vmf(pa, na, window, scores);
        for (threshold = 0; threshold < window*window; threshold++)
            printf("vmf(%d)=%.3f,",threshold+1,scores[threshold]);

        printf("cmf=%.3f,smf=%.3f,mean=%.3f", calculatePotentialRMSE_cmf(pa, na, window), calculatePotentialRMSE_smf(pa, na, window), calculatePotentialRMSE_mean(pa, na, window));
        printf("\n");
        freeImage(na);
    }

    free(scores);
    freeImage(pa);
    fclose(a);
    return 0;
} // End Main
