#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include "imgtypes.h"
#include "minkowski.h"

void thresholdImage(struct portImage *pi) {
    int half = pi->scale >> 1;
    int i;
    for (i = 0; i < pi->size; i++)
        if (pi->f[i] > half)
            pi->f[i] = pi->scale;
        else
            pi->f[i] = 0;
}

void minkowskiDivision(struct portImage *pi_n, struct portImage *pi_d) {
    int i;

    if (pi_n->color == 1)
        return;

    for (i = 0; i < pi_n->size; i++)
        if (pi_n->f[i] != pi_d->f[i])
            pi_n->f[i] = pi_n->scale;
        else
            pi_n->f[i] = 0;
}

void minkowskiAddition(struct portImage *pi, int maskwidth, int maskheight) {

    int i,j,k,l;
    unsigned char *b;
    int maskheight_half = maskheight >> 1;
    int maskwidth_half  = maskwidth >> 1;

    if (pi->color == 1)
        return;

    b = (unsigned char*)malloc( sizeof(unsigned char) * pi->size );
    if (!b) {
        fprintf(stderr, "ERROR: failed to allocate memory for minkowski operation.\n");
        return;
    }

    for (j = 0; j < pi->height; j++)
        for (i = 0; i < pi->width; i++)
            b[j * pi->width + i] = 0;

    for (j = maskheight_half; j < pi->height - maskheight_half; j++) {
        for (i = maskwidth_half; i < pi->width - maskwidth_half; i++) {
            int any = 0;
            for (k = j - maskheight_half; k <= j + maskheight_half; k++) {
                for (l = i - maskwidth_half; l <= i + maskwidth_half; l++) {
                    if (pi->f[k * pi->width + l] == pi->scale) {
                        any = 1;
                        break;
                    }
                }
            }

            if (any == 1) {
                b[j * pi->width + i] = pi->scale;
            }
        }
    }

    for (i = 0; i < pi->size; i++)
        pi->f[i] = b[i];

    free(b);
}

void minkowskiSubtraction(struct portImage *pi, int maskwidth, int maskheight) {

    int i,j,k,l;
    unsigned char *b;
    int maskheight_half = maskheight >> 1;
    int maskwidth_half  = maskwidth >> 1;

    if (pi->color == 1)
        return;

    b = (unsigned char*)malloc( sizeof(unsigned char) * pi->size );
    if (!b) {
        fprintf(stderr, "ERROR: failed to allocate memory for minkowski operation.\n");
        return;
    }

    for (j = 0; j < pi->height; j++)
        for (i = 0; i < pi->width; i++)
            b[j * pi->width + i] = 0;

    for (j = maskheight_half; j < pi->height - maskheight_half; j++) {
        for (i = maskwidth_half; i < pi->width - maskwidth_half; i++) {
            int all = 1;
            for (k = j - maskheight_half; k <= j + maskheight_half; k++) {
                for (l = i - maskwidth_half; l <= i + maskwidth_half; l++) {
                    if (pi->f[k * pi->width + l] != pi->scale) {
                        all = 0;
                        break;
                    }
                }
            }

            if (all == 1) {
                b[j * pi->width + i] = pi->scale;
            }
        }
    }

    for (i = 0; i < pi->size; i++)
        pi->f[i] = b[i];

    free(b);
}

void minkowskiOpening(struct portImage *pi, int maskwidth, int maskheight) {
    minkowskiSubtraction(pi, maskwidth, maskheight);
    minkowskiAddition(pi, maskwidth, maskheight);
}

void minkowskiClosing(struct portImage *pi, int maskwidth, int maskheight) {
    minkowskiAddition(pi, maskwidth, maskheight);
    minkowskiSubtraction(pi, maskwidth, maskheight);
}
