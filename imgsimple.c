#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include "imgsimple.h"

void scale_reduce(struct portImage *pi, int newscale) {

    int i;
    int reduction = (pi->scale+1) / (newscale);

    for (i = 0; i < pi->size; i++)
        pi->f[i] /= reduction;

    pi->scale = newscale-1;
}

void spacial_reduce(struct portImage *pi, int resWidth, int resHeigth) {

    int reduce_h = pi->height / resWidth;
    int reduce_w = pi->width  / resHeigth;
    int w    = pi->width;
    if (pi->color == csRGB) {
        int i, j, x, y, average_red, average_green, average_blue;
        w *= 3;

        for (j = 0; j < pi->height; j += reduce_h) {
            for (i = 0; i < w; i += reduce_w*3) {

                average_red = 0;
                average_green = 0;
                average_blue = 0;
                for (y = 0; y < reduce_h; y++)
                    for (x = 0; x < reduce_w; x++) {
                        average_red   += pi->f[(j+y)*w + (i+x) + 0];
                        average_green += pi->f[(j+y)*w + (i+x) + 1];
                        average_blue  += pi->f[(j+y)*w + (i+x) + 2];
                    }

                average_red   /= reduce_h*reduce_w;
                average_green /= reduce_h*reduce_w;
                average_blue  /= reduce_h*reduce_w;

                for (y = 0; y < reduce_h; y++)
                    for (x = 0; x < reduce_w; x++) {
                        pi->f[(j+y)*w + (i+x) + 0] = average_red;
                        pi->f[(j+y)*w + (i+x) + 1] = average_green;
                        pi->f[(j+y)*w + (i+x) + 2] = average_blue;
                    }
            }
        }

    }
    else if (pi->color == csGREYSCALE) {
        int i, j, x, y, average = 0;
        for (j = 0; j < pi->height; j += reduce_h) {
            for (i = 0; i < w; i += reduce_w) {

                average = 0;
                for (y = 0; y < reduce_h; y++)
                    for (x = 0; x < reduce_w; x++)
                        average += pi->f[(j+y)*w + (i+x)];

                average /= reduce_h*reduce_w;

                for (y = 0; y < reduce_h; y++)
                    for (x = 0; x < reduce_w; x++)
                        pi->f[(j+y)*w + (i+x)] = average;
            }
        }
    }
}

void histogram(struct portImage *pi, int *hist) {
    int i;
    if (pi->color != csGREYSCALE)
        return;

    for (i = 0; i <= pi->scale; i++)
        hist[i] = 0;

    for (i = 0; i < (pi->height * pi->width); i++)
        hist[ (int)pi->f[i] ]++;
}

void cumulative_distribution(int scales, int *hist, int *cdf) {
    int i;
    int total = 0;
    for (i = 0; i < scales; i++) {
        total += hist[i];
        cdf[i] = total;
    }
}

void graph_histogram(struct portImage *pi, char *filename) {
    int i, j;
    int max;
    int width = pi->scale + 1;
    int height;
    int *hist;
    double scale;
    FILE *out;
    struct portImage histogram_out;

    histogram_out.height = width;
    histogram_out.width  = width;
    histogram_out.scale  = 1;
    histogram_out.color  = csGREYSCALE;
    histogram_out.type   = 2;
    histogram_out.c      = 0;

    histogram_out.f = (unsigned char*)malloc(sizeof(unsigned char) * width * width);
    hist = (int*)malloc(width * sizeof(int));

    if (!hist || !histogram_out.f) {
        fprintf(stderr, "ERROR: Memory for memory allocation of histogram failed.\n");
        return;
    }

    histogram(pi, hist);

    max = hist[0];
    for (i = 1; i < width; i++) /* For each bar. */
        if (hist[i] > max)
            max = hist[i];

    scale = (double)width / (double)max;

    for (i = 0; i < width; i++) { /* For each bar. */
        height = (int) ( scale * (double)hist[i] );

        for (j = width-height; j <= pi->scale; j++)
            histogram_out.f[j*width+i] = 1;
    }

    out = fopen(filename, "w");
    if (out == 0) {
        fprintf(stderr, "ERROR: Failed to open %s for writing histogram.\n", filename);
        exit(1);
    }

    writeImage(&histogram_out, out);
    fclose(out);
    free(hist);
    free(histogram_out.f);
}

void graph_cdf(struct portImage *pi, char *filename) {
    int i, j;
    int max;
    int width = pi->scale + 1;
    int height;
    int *hist;
    int *cdf;
    double scale;
    FILE *out;
    struct portImage cdf_out;

    cdf_out.height = width;
    cdf_out.width  = width;
    cdf_out.scale  = 1;
    cdf_out.color  = csGREYSCALE;
    cdf_out.type   = 2;
    cdf_out.c      = 0;

    cdf_out.f = (unsigned char*)malloc(sizeof(unsigned char) * width * width);
    hist = (int*)malloc(width * sizeof(int));
    cdf  = (int*)malloc(width * sizeof(int));

    if (!hist || !cdf || !cdf_out.f) {
        fprintf(stderr, "Memory allocation for CDF failed.\n");
        return;
    }

    histogram(pi, hist);
    cumulative_distribution(width, hist, cdf);

    max = cdf[width-1];
    scale = (double)width / (double)max;

    for (i = 0; i < width; i++) { /* For each bar. */
        height = (int) ( scale * (double)cdf[i] );

        for (j = width-height; j <= pi->scale; j++)
            cdf_out.f[j*width+i] = 1;
    }

    out = fopen(filename, "w");
    if (out == 0) {
        fprintf(stderr, "ERROR: Failed to open %s for writing histogram.\n", filename);
        exit(1);
    }

    writeImage(&cdf_out, out);
    fclose(out);
    free(cdf_out.f);
    free(hist);
    free(cdf);
}

void level_slice(struct portImage *pi, int level) {
    int i;
    for (i = 0; i < (pi->height * pi->width); i++)
        if  (level <= pi->f[i] && pi->f[i] <= level+10)
            pi->f[i] = pi->scale;
        else
            pi->f[i] = 0;
}

void equalize(struct portImage *pi) {
    int i;
    int *hist;
    int *cdf;
    int *nhist;

    hist  = (int*)malloc( (pi->scale + 1) * sizeof(int) );
    cdf   = (int*)malloc( (pi->scale + 1) * sizeof(int) );
    nhist = (int*)malloc( (pi->scale + 1) * sizeof(int) );

    if (!hist || !cdf || !nhist) {
        fprintf(stderr, "Memory allocation for image equalization failed.\n");
        return;
    }

    if (pi->color == csGREYSCALE) {
        histogram(pi, hist);
        cumulative_distribution(pi->scale+1, hist, cdf);

        for (i = 0; i <= pi->scale; i++) {
            double Crr = ( (double)cdf[i] / (double)(pi->height * pi->width) ) * (double)pi->scale;

            double Crr_dec = Crr - (double)(int)Crr;

            nhist[i] = (int)Crr;

            if (Crr_dec >= 0.5)
                nhist[i] += 1;
        }

        for (i = 0; i < (pi->height * pi->width); i++)
            pi->f[i] = nhist[ pi->f[i] ];
    }

    free(hist);
    free(cdf);
    free(nhist);
}

void contrast_stretching(struct portImage *pi, float low_bound, float high_bound) {
    int i;
    int *hist;
    int *cdf;
    int low_scale  = -1;
    int high_scale = -1;
    int low_pixels  = (int) ( (float)(pi->height * pi->width) * low_bound  );
    int high_pixels = (int) ( (float)(pi->height * pi->width) * high_bound );
    float factor;

    hist  = (int*)malloc( (pi->scale + 1) * sizeof(int) );
    cdf   = (int*)malloc( (pi->scale + 1) * sizeof(int) );

    if (!hist || !cdf) {
        fprintf(stderr, "Memory allocation for image equalization failed.\n");
        return;
    }

    histogram(pi, hist);
    cumulative_distribution(pi->scale+1, hist, cdf);

    for (i = 0; i <= pi->scale; i++) {
        if (low_scale  == -1 && cdf[i] > low_pixels)  low_scale  = i-1;
        if (high_scale == -1 && cdf[i] > high_pixels) high_scale = i;
    }

    factor = (float)(pi->scale) / (float)(high_scale - low_scale);

    for (i = 0; i < (pi->height * pi->width); i++) {
        /* Set all pixels in the bottom lowbound to black */
        if      (pi->f[i] < low_scale)
            pi->f[i] = 0;

        /* Set all of the pixels in the upper highbound to white */
        else if (pi->f[i] > high_scale)
            pi->f[i] = pi->scale;

        /* Scale whatever is left. */
        else
            pi->f[i] = (int)(factor * (float)(pi->f[i] - low_scale));
    }
    free(hist);
    free(cdf);
}
