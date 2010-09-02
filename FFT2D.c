#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include "imgtypes.h"
#include "FFT.h"
#include "FFT2D.h"

double logMag(struct Complex p) {
    return log10( sqrt( p.real * p.real + p.imag * p.imag ) );
}

void writeComplex(struct portImage *pi, FILE* out) {
    int i, j, k=0;
    int w = pi->width;
    double point;
    double max;
    double scale;
    double *buffer;

    if (pi->c == 0) {
        fprintf(stderr, "No complex values to print.\n");
        return;
    }

    if (pi->color)
        w *= 3;

    buffer = (double*)malloc( sizeof(double) * pi->size );

    if (!buffer) {
        fprintf(stderr, "Error: allocation of memory for Complex write failed.\n");
        return;
    }

    max = logMag(pi->c[0]);
    for (i = 1; i < pi->size; i++) { /* For each pixel. */
        point = logMag(pi->c[i]);
        if (point > max)
            max = point;
    }

    scale = (double)pi->scale / (double)max;

    for (i = pi->height / 2; i < pi->height; i++) {
        for (j = pi->width / 2; j < pi->width; j++)
            buffer[k++] = scale * logMag(pi->c[i*w+j]);
        for (j = 0; j < pi->width / 2; j++)
            buffer[k++] = scale * logMag(pi->c[i*w+j]);
    }

    for (i = 0; i < pi->height / 2; i++) {
        for (j = pi->width / 2; j < pi->width; j++)
            buffer[k++] = scale * logMag(pi->c[i*w+j]);
        for (j = 0; j < pi->width / 2; j++)
            buffer[k++] = scale * logMag(pi->c[i*w+j]);
    }

    /* Header: */
    fprintf(out, "P%d\n#\n%d %d\n%d\n", pi->type, pi->width, pi->height, pi->scale);

    /* ASCII Format */
    if (pi->type == 1 || pi->type == 2 || pi->type == 3)
        for (i = 0; i < pi->height; i++) {
            for (j = 0; j < w; j++)
                fprintf(out, "%3d ", (unsigned char)buffer[i*w+j]);

            fprintf(out, "\n");
        }

    /* Binary Format */
    else
        for (i = 0; i < pi->size; i++) {
            unsigned char x = (unsigned char) buffer[i];
            fwrite(&x, sizeof(unsigned char), 1, out);
        }

    free(buffer);
}

void FFT2D(struct portImage *pi) {
    int i, j;
    struct Complex *buffer;

    /* If already computed, don't recompute. */
    if (pi->c)
        return;

    /* Works only for greyscale images for now */
    if (pi->color != csGREYSCALE)
        return;

    pi->c  = (struct Complex*) malloc( sizeof(struct Complex) * pi->size   );
    buffer = (struct Complex*) malloc( sizeof(struct Complex) * pi->height );

    if (!(pi->c) || !buffer) {
        fprintf(stderr, "ERROR: memory allocation for Complex buffer failed.");
        return;
    }

    /* Transfer image to complex datatype */
    for (i = 0; i < pi->size; i++) {
        pi->c[i].real = (double) pi->f[i];
        pi->c[i].imag = 0.0;
    }

    /* Compute the FFT for each row */
    for (i = 0; i < pi->height; i++) {
        
        FFT(&pi->c[i * pi->width], pi->width); /*Uncomment */
    }

    /* Compute the FFT for each column */
    for (i = 0; i < pi->width; i++) {
        for (j = 0; j < pi->height; j++)
            buffer[j] = pi->c[j * pi->width + i];

        FFT(buffer, pi->width); /*Uncomment */

        for (j = 0; j < pi->height; j++)
            pi->c[j * pi->width + i] = buffer[j];
    }

    free(buffer);
}

void IFFT2D(struct portImage *pi) {
    int i, j;
    struct Complex *buffer;

    /* If original not computed, don't compute */
    if (pi->c == 0)
        return;

    /* Works only for greyscale images for now */
    if (pi->color != csGREYSCALE)
        return;

    buffer = (struct Complex*) malloc( sizeof(struct Complex) * pi->height );

    if (!buffer) {
        fprintf(stderr, "ERROR: memory allocation for Complex buffer failed.");
        return;
    }

    /* Compute the FFT for each row */
    for (i = 0; i < pi->height; i++)
        IFFT(&pi->c[i * pi->width], pi->width); /*Uncomment */

    /* Compute the FFT for each column */
    for (i = 0; i < pi->width; i++) {
        for (j = 0; j < pi->height; j++)
            buffer[j] = pi->c[j * pi->width + i];

        IFFT(buffer, pi->width); /*Uncomment */

        for (j = 0; j < pi->height; j++)
            pi->c[j * pi->width + i] = buffer[j];
    }

    for (i = 0; i < pi->size; i++)
        pi->f[i] = (unsigned char)pi->c[i].real;

    free(buffer);
}

void graph_fftlogplot(struct portImage *pi, char *filename) {
    FILE *out;

    struct portImage fftlogplot_out;

    if (pi->color != csGREYSCALE)
        return;

    out = fopen(filename, "w");
    if (out == 0) {
        fprintf(stderr, "ERROR: Failed to open %s for writing histogram.\n", filename);
        exit(1);
    }

    fftlogplot_out.height = pi->height;
    fftlogplot_out.width  = pi->width;
    fftlogplot_out.scale  = 1;
    fftlogplot_out.color  = csGREYSCALE;
    fftlogplot_out.type   = 2;
    fftlogplot_out.c      = 0;

    fftlogplot_out.f = (unsigned char*)malloc(sizeof(unsigned char) * pi->size);
    if (fftlogplot_out.f == 0) {
        fprintf(stderr, "ERROR: Memory for memory allocation of fftlogplot failed.\n");
        return;
    }

    FFT2D(pi);
    writeComplex(pi, out);

    if (pi->c == 0) {
        fprintf(stderr, "FFT Failed.\n");
        return;
    }

    fclose(out);
    free(fftlogplot_out.f);
}

void lowpass(struct portImage *pi, double percent) {
    int i, j;
    int pass_width  = (int)(percent * 0.5 * (double)pi->width);
    int pass_height = (int)(percent * 0.5 * (double)pi->height);

    if (pi->color != csGREYSCALE)
        return;

    for (i = 0; i < pi->height; i++) {
        for (j = 0; j < pi->width; j++) {
            if ( (i >= pass_height && (pi->height - i) >= pass_height) ||
                 (j >= pass_width  && (pi->width  - j) >= pass_width ) ) {
                 pi->c[i * pi->width + j].real = 0.0;
                 pi->c[i * pi->width + j].imag = 0.0;
            }
        }
    }
}

void highpass(struct portImage *pi, double percent) {
    int i, j;
    int pass_width  = (int)(percent * 0.5 * (double)pi->width);
    int pass_height = (int)(percent * 0.5 * (double)pi->height);
    int half_width  = pi->width >> 1;
    int half_height = pi->height >> 1;

    if (pi->color != csGREYSCALE)
        return;

    for (i = 0; i < pi->height; i++)
        for (j = 0; j < pi->width; j++)
            if ( (i <= (half_height-pass_height) ||  i >= (half_height+pass_height)) ||
                 (j <= (half_width -pass_width ) ||  j >= (half_width +pass_width)) ) {
                 pi->c[i * pi->width + j].real = 0.0;
                 pi->c[i * pi->width + j].imag = 0.0;
            }

}

double fftSimilarity(struct portImage *pa, struct portImage *pb, double percent) {

    FFT2D(pa);
    FFT2D(pb);

    highpass(pa, percent);
    highpass(pb, percent);
    /*
    lowpass(pa, percent);
    lowpass(pb, percent);
    */

    int i;
    double accumulator = 0.0;
    for (i = 0; i < pa->size; i++) {
        accumulator += (pa->c[i].real - pb->c[i].real)*(pa->c[i].real - pb->c[i].real);
    }

    return accumulator / (double) pa->size;
}
