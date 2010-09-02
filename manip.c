#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "imgio.h"
#include "minkowski.h"
#include "FFT2D.h"
#include "portable.h"

#define FILENAME_LENGTH 100

int main(int argc, char **argv) {

    int commandArgs = 1;
    int i;
    char x;
    char performContrast = 0;
    char fft_filename[FILENAME_LENGTH];
    char cdf_filename[FILENAME_LENGTH];
    char histo_filename[FILENAME_LENGTH];
    float contrastLow = 0.0;
    float contrastHigh = 0.0;
    float highpasslevel;
    float lowpasslevel;
    float percentCorrupt;
    float sigma;
    float brightness;
    float sat;
    int m;
    int replaceWithM;
    int performHistogram = 0;
    int performCDF = 0;
    int performFFT = 0;
    int performVectorMedianFilter = 0;
    int performMedianFilter = 0;
    int performMeanFilter = 0;
    int performSpacialFilter = 0;
    int performLevelSlicing = 0;
    int performEqualize = 0;
    int performColorScale = 0;
    int performSpatialReduceWidth  = 0;
    int performSpatialReduceHeigth = 0;
    int performHighpass = 0;
    int performLowpass = 0;
    int performComponentMedianFilter = 0;
    int performVectorOrderStatistic = 0;
    int performVectorSpacialOrderStatistic = 0;
    int performVectorMedianOrderStatistic = 0;
    int performMinkowskiAddition = 0;
    int performMinkowskiSubtraction = 0;
    int performMinkowskiOpening = 0;
    int performMinkowskiClosing = 0;
    int performEdgeDetectOne = 0;
    int performEdgeDetectTwo = 0;
    int performEdgeDetectThree = 0;
    int performEdgeDetectFour = 0;
    int performAddGaussianNoise = 0;
    int performAddSaltPepperNoise = 0;
    int performSetBrightness = 0;
    int performSetSaturation = 0;
    int performBrightFilter = 0;
    int imageWriteNecessary = 0;
    int maskwidth = 0;
    int maskheight = 0;
    int maskrepeat = 0;
    FILE *in = stdin;
    FILE *out = stdout;
    struct portImage *pi;

    if (argc < 3) {
        printf("\n");
        printf("  Usage: %s [inputFile] ([outputFile]) [option] ([option] ...)\n", argv[0]);
        printf("\n");
        printf("  InputFile: Either a filename or '-' for stdin.\n");
        printf("  OutputFile: Either a filename or '-' for stdout. (Not needed if no output necessary.)\n");
        printf("\n");
        printf("  Options:\n");
        printf("  -ghisto FILENAME       Graph Histogram\n");
        printf("  -gcdf FILENAME         Graph Cumulative Distribution\n");
        printf("  -gfft FILENAME         Graph FFT plot\n");
        printf("  -color n               Reduce color scale to n\n");
        printf("  -spatial WIDTH-HEIGHT  Perform spacial reduction to Width and Height\n");
        printf("  -level n               Perform level slicing from graylevel n to graylevel n+10\n");
        printf("  -con LOW-HIGH          Scale image contrast from LOW graylevel percentage to HIGH graylevel percentage\n");
        printf("  -equ                   Histogram Eqaulization\n");
        printf("  -medianf n             Simple Median Filter of window size n*n\n");
        printf("  -meanf n               Simple Mean Filter of window size n*n\n");
        printf("  -cmf n                 Component Median Filter of window size n*n\n");
        printf("  -vmf n                 Vector Median Filter of window n*n\n");
        printf("  -sf                    Spacial Filter\n");
        printf("  -vos n v               Vector Order Stat of window size n*n and value v\n");
        printf("  -vmos n m [01]         Vector Median Order Stat of window size n*n, m threshold, and 0 or 1(True) replace with m\n");
        printf("  -vsos n m [01]         Vector Spacial Order Stat of window size n*n, m threshold, and 0 or 1(True) replace with m\n");
        printf("  -brightf n             Perform an Brightness filter using HSV colorspace on size n*n window.\n");
        printf("  -bright %%              Set brightness to %% percent\n");
        printf("  -sat %%                 Set saturation to %% percent\n");
        printf("  -hp %%                  Highpass filter of %% percent\n");
        printf("  -lp %%                  Lowpass filter of %% percent\n");
        printf("  -ma NxM R              Perform Minkowski Addition using NxM mask, repeated R times.\n");
        printf("  -ms NxM R              Perform Minkowski Subtraction using NxM mask, repeated R times.\n");
        printf("  -mo NxM R              Perform Minkowski Opening using NxM mask, repeated R times.\n");
        printf("  -mc NxM R              Perform Minkowski Closing using NxM mask, repeated R times.\n");
        printf("  -e1                    Perform Edge Detection using X/(X – B)\n");
        printf("  -e2                    Perform Edge Detection using (X + B)/X\n");
        printf("  -e3                    Perform Edge Detection using [(X+B)/(X-B)]-B\n");
        printf("  -e4 n                  Experimental Edge Detection on Color Images using n*n window.\n");
        printf("  -noiseG p s            Add Gaussian noise to p (0 to 1 floating) percent of image with s sigma noise.\n");
        printf("  -noiseSP p             Add Salt and Pepper noise to p (0 to 1 floating) percent of image.\n");
        printf("\n");
        return(1);
    }

    if (strcmp(argv[commandArgs], "-") != 0) {
        in = fopen(argv[1],"r");

        if (in == NULL) {
            fprintf(stderr, "File '%s' failed to open for reading.\n", argv[1]);
            exit(1);
        }
    }
    commandArgs++;

    if (strcmp(argv[commandArgs], "-") != 0 && argv[commandArgs][0] != '-') {
        commandArgs++;
        out = fopen(argv[2],"w");

        if (out == NULL) {
            fprintf(stderr, "File '%s' failed to open for writing.\n", argv[2]);
            exit(1);
        }
    }

    for (; commandArgs < argc; commandArgs++) {
        if (strcmp(argv[commandArgs], "-color") == 0) {
            imageWriteNecessary = 1;
            commandArgs++;
            sscanf(argv[commandArgs], "%d", &performColorScale);
        }

        if (strcmp(argv[commandArgs], "-spatial") == 0) {
            imageWriteNecessary = 1;
            commandArgs++;
            sscanf(argv[commandArgs], "%d%c%d", &performSpatialReduceWidth, &x, &performSpatialReduceHeigth);
        }

        if (strcmp(argv[commandArgs], "-level") == 0) {
            imageWriteNecessary = 1;
            commandArgs++;
            sscanf(argv[commandArgs], "%d", &performLevelSlicing);
        }

        if (strcmp(argv[commandArgs], "-con") == 0) {
            imageWriteNecessary = 1;
            commandArgs++;
            sscanf(argv[commandArgs], "%f%c%f", &contrastLow, &performContrast, &contrastHigh);
        }

        if (strcmp(argv[commandArgs], "-vos") == 0) {
            imageWriteNecessary = 1;
            commandArgs++;
            sscanf(argv[commandArgs], "%d", &performVectorOrderStatistic);
            commandArgs++;
            sscanf(argv[commandArgs], "%d", &m);
        }

        if (strcmp(argv[commandArgs], "-vmf") == 0) {
            imageWriteNecessary = 1;
            commandArgs++;
            sscanf(argv[commandArgs], "%d", &performVectorMedianFilter);
        }

        if (strcmp(argv[commandArgs], "-sf") == 0) {
            imageWriteNecessary = 1;
            performSpacialFilter = 1;
        }

        if (strcmp(argv[commandArgs], "-vmos") == 0) {
            imageWriteNecessary = 1;
            commandArgs++;
            sscanf(argv[commandArgs], "%d", &performVectorMedianOrderStatistic);
            commandArgs++;
            sscanf(argv[commandArgs], "%d", &m);
            commandArgs++;
            sscanf(argv[commandArgs], "%d", &replaceWithM);
        }

        if (strcmp(argv[commandArgs], "-vsos") == 0) {
            imageWriteNecessary = 1;
            commandArgs++;
            sscanf(argv[commandArgs], "%d", &performVectorSpacialOrderStatistic);
            commandArgs++;
            sscanf(argv[commandArgs], "%d", &m);
            commandArgs++;
            sscanf(argv[commandArgs], "%d", &replaceWithM);
        }

        if (strcmp(argv[commandArgs], "-cmf") == 0) {
            imageWriteNecessary = 1;
            commandArgs++;
            sscanf(argv[commandArgs], "%d", &performComponentMedianFilter);
        }

        if (strcmp(argv[commandArgs], "-medianf") == 0) {
            imageWriteNecessary = 1;
            commandArgs++;
            sscanf(argv[commandArgs], "%d", &performMedianFilter);
        }

        if (strcmp(argv[commandArgs], "-meanf") == 0) {
            imageWriteNecessary = 1;
            commandArgs++;
            sscanf(argv[commandArgs], "%d", &performMeanFilter);
        }

        if (strcmp(argv[commandArgs], "-equ") == 0) {
            imageWriteNecessary = 1;
            performEqualize = 1;
        }

        if (strcmp(argv[commandArgs], "-hp") == 0) {
            imageWriteNecessary = 1;
            commandArgs++;
            performHighpass = 1;
            sscanf(argv[commandArgs], "%f", &highpasslevel);
        }

        if (strcmp(argv[commandArgs], "-lp") == 0) {
            imageWriteNecessary = 1;
            commandArgs++;
            performLowpass = 1;
            sscanf(argv[commandArgs], "%f", &lowpasslevel);
        }

        if (strcmp(argv[commandArgs], "-brightf") == 0) {
            imageWriteNecessary = 1;
            commandArgs++;
            sscanf(argv[commandArgs], "%d", &performBrightFilter);
        }

       if (strcmp(argv[commandArgs], "-bright") == 0) {
            imageWriteNecessary = 1;
            commandArgs++;
            performSetBrightness = 1;
            sscanf(argv[commandArgs], "%f", &brightness);
        }

        if (strcmp(argv[commandArgs], "-sat") == 0) {
            imageWriteNecessary = 1;
            commandArgs++;
            performSetSaturation = 1;
            sscanf(argv[commandArgs], "%f", &sat);
        }

        if (strcmp(argv[commandArgs], "-ghisto") == 0) {
            commandArgs++;
            performHistogram = 1;
            strncpy(histo_filename, argv[commandArgs], FILENAME_LENGTH);
        }

        if (strcmp(argv[commandArgs], "-gcdf") == 0) {
            commandArgs++;
            performCDF = 1;
            strncpy(cdf_filename, argv[commandArgs], FILENAME_LENGTH);
        }

        if (strcmp(argv[commandArgs], "-gfft") == 0) {
            commandArgs++;
            performFFT = 1;
            strncpy(fft_filename, argv[commandArgs], FILENAME_LENGTH);
        }

        if (strcmp(argv[commandArgs], "-ma") == 0) {
            imageWriteNecessary = 1;
            performMinkowskiAddition = 1;
            commandArgs++;
            sscanf(argv[commandArgs], "%d%c%d", &maskwidth, &x, &maskheight);
            commandArgs++;
            sscanf(argv[commandArgs], "%d", &maskrepeat);
        }

        if (strcmp(argv[commandArgs], "-ms") == 0) {
            imageWriteNecessary = 1;
            performMinkowskiSubtraction = 1;
            commandArgs++;
            sscanf(argv[commandArgs], "%d%c%d", &maskwidth, &x, &maskheight);
            commandArgs++;
            sscanf(argv[commandArgs], "%d", &maskrepeat);
        }

        if (strcmp(argv[commandArgs], "-mo") == 0) {
            imageWriteNecessary = 1;
            performMinkowskiOpening = 1;
            commandArgs++;
            sscanf(argv[commandArgs], "%d%c%d", &maskwidth, &x, &maskheight);
            commandArgs++;
            sscanf(argv[commandArgs], "%d", &maskrepeat);
        }

        if (strcmp(argv[commandArgs], "-mc") == 0) {
            imageWriteNecessary = 1;
            performMinkowskiClosing = 1;
            commandArgs++;
            sscanf(argv[commandArgs], "%d%c%d", &maskwidth, &x, &maskheight);
            commandArgs++;
            sscanf(argv[commandArgs], "%d", &maskrepeat);
        }

        if (strcmp(argv[commandArgs], "-e1") == 0) {
            imageWriteNecessary = 1;
            performEdgeDetectOne = 1;
        }

        if (strcmp(argv[commandArgs], "-e2") == 0) {
            imageWriteNecessary = 1;
            performEdgeDetectTwo = 1;
        }

        if (strcmp(argv[commandArgs], "-e3") == 0) {
            imageWriteNecessary = 1;
            performEdgeDetectThree = 1;
        }

        if (strcmp(argv[commandArgs], "-e4") == 0) {
            imageWriteNecessary = 1;
            commandArgs++;
            sscanf(argv[commandArgs], "%d", &performEdgeDetectFour);
        }

        if (strcmp(argv[commandArgs], "-noiseG") == 0) {
            imageWriteNecessary = 1;
            performAddGaussianNoise = 1;
            commandArgs++;
            sscanf(argv[commandArgs], "%f", &percentCorrupt);
            commandArgs++;
            sscanf(argv[commandArgs], "%f", &sigma);
        }

        if (strcmp(argv[commandArgs], "-noiseSP") == 0) {
            imageWriteNecessary = 1;
            performAddSaltPepperNoise = 1;
            commandArgs++;
            sscanf(argv[commandArgs], "%f", &percentCorrupt);
        }
    }

    pi = readImage(in);

    if (performHighpass || performLowpass || performFFT) {
        FFT2D(pi);
        if (performHighpass)  highpass(pi, highpasslevel);
        if (performLowpass)   lowpass(pi, lowpasslevel);
        if (performFFT)       graph_fftlogplot(pi, fft_filename);
        IFFT2D(pi);
    }

    if (performEdgeDetectOne || performEdgeDetectTwo || performEdgeDetectThree ||
        performMinkowskiAddition || performMinkowskiSubtraction || performMinkowskiOpening || performMinkowskiClosing)
        thresholdImage(pi);

    if (performAddGaussianNoise)             addGaussianNoise(pi, percentCorrupt, sigma);
    if (performAddSaltPepperNoise)           addSaltPepperNoise(pi, percentCorrupt);
    if (performMedianFilter)                 simpleMedianFilter(pi, performMedianFilter);
    if (performMeanFilter)                   simpleMeanFilter(pi, performMeanFilter);
    if (performComponentMedianFilter)        componentMedianFilter(pi, performComponentMedianFilter);
    if (performVectorOrderStatistic)         vectorOrderStatistic(pi, performVectorOrderStatistic, m);
    if (performVectorSpacialOrderStatistic)  vectorSpacialOrderStatistic(pi, performVectorSpacialOrderStatistic, m, replaceWithM);
    if (performVectorMedianOrderStatistic)   vectorMedianOrderStatistic(pi, performVectorMedianOrderStatistic, m, replaceWithM);
    if (performVectorMedianFilter)           vectorMedianFilter(pi, performVectorMedianFilter);
    if (performSpacialFilter)                spacialFilter(pi);
    if (performBrightFilter)                 HSV_ValueFilter(pi, performBrightFilter);
    if (performColorScale)                   scale_reduce(pi,performColorScale);
    if (performSetBrightness)                setBrightness(pi,brightness);
    if (performSetSaturation)                setSaturation(pi,sat);
    if (performSpatialReduceWidth)           spacial_reduce(pi,performSpatialReduceWidth, performSpatialReduceHeigth);
    if (performContrast)                     contrast_stretching(pi, contrastLow, contrastHigh);
    if (performLevelSlicing)                 level_slice(pi, performLevelSlicing);
    if (performEqualize)                     equalize(pi);
    if (performHistogram)                    graph_histogram(pi, histo_filename);
    if (performCDF)                          graph_cdf(pi, cdf_filename);
    if (performEdgeDetectFour)               spacialEdgeDetect(pi, performEdgeDetectFour);

    if (performMinkowskiAddition)
        for (i = 0; i < maskrepeat; i++)
            minkowskiAddition(pi, maskwidth, maskheight);

    if (performMinkowskiSubtraction)
        for (i = 0; i < maskrepeat; i++)
            minkowskiSubtraction(pi, maskwidth, maskheight);

    if (performMinkowskiOpening)
        for (i = 0; i < maskrepeat; i++)
            minkowskiOpening(pi, maskwidth, maskheight);

    if (performMinkowskiClosing)
        for (i = 0; i < maskrepeat; i++)
            minkowskiClosing(pi, maskwidth, maskheight);

    if (performEdgeDetectOne) {
        struct portImage *pc = copyImage(pi);
        imageWriteNecessary = 1;

        minkowskiSubtraction(pc, 3, 3);
        minkowskiDivision(pi, pc);

        freeImage(pc);
    }

    if (performEdgeDetectTwo) {
        struct portImage *pc = copyImage(pi);
        imageWriteNecessary = 1;

        minkowskiAddition(pi, 3, 3);
        minkowskiDivision(pi, pc);

        freeImage(pc);
    }

    if (performEdgeDetectThree) {
        struct portImage *pd = copyImage(pi);
        maskrepeat = 3;
        imageWriteNecessary = 1;

        for (i = 0; i < maskrepeat; i++) {
            minkowskiAddition(pi, 3, 3);
            minkowskiSubtraction(pd, 3, 3);
        }

        minkowskiDivision(pi, pd);
        minkowskiSubtraction(pi, 3, 3);

        freeImage(pd);
    }

    if (imageWriteNecessary)
        writeImage(pi, out);

    freeImage(pi);
    if (in  != stdin)  fclose(in);
    if (out != stdout) fclose(out);
    return 0;
} /* End Main */
