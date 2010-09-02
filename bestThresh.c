#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "portable.h"

#define FILENAME_LENGTH 100

int main(int argc, char **argv) {

    int i;
    int window;
    FILE *a;
    double noise[] = {0.0,0.05,0.10,0.15,0.20,0.25,0.30,0.35,0.40,0.45,0.50};
    portImage *pa;
    a = fopen(argv[1], "r");
    sscanf(argv[2], "%d", &window);
    pa = readImage(a);

    printf("%s,", argv[1]);
    for (i = 0; i < 11; i++)
        printf("%d,",findBestThreshold(pa, window, noise[i]));
    printf("\n");

    freeImage(pa);
    fclose(a);
    return 0;
} // End Main
