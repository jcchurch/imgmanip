#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include "noise.h"

double randDouble() {
    return rand() / (double) RAND_MAX;
}

void addSaltPepperNoise(struct portImage *pi, double percentCorrupt) {
    int i;
    int x;
    int noise;

    for (i = 0; i < pi->width * pi->height; i++) {
        if (randDouble() < percentCorrupt) {
            x     = (int)(randDouble() * 3.0);                     /* Random number from 0 to 2 */
            noise = (int)(randDouble() * (double)(pi->scale + 1)); /* 0 to 255 */
            pi->f[i*3+x] = noise;
        }
    }
}

/* This is a bad noise model. Do not use it. */
void addGaussianNoise(struct portImage *pi, double percentCorrupt, double sigma) {
    int i;
    int pixelsToCorrupt = (int) ((double)pi->height * pi->width * percentCorrupt);
    double rand_num1, rand_num2;
    int new_value, pos;

    for (i = 0; i < pixelsToCorrupt; i++) {
        rand_num1 = randDouble();
        rand_num2 = randDouble();
        pos = rand() % pi->size;
        new_value = pi->f[ pos ] + (int) (sqrt( -2*log(rand_num1) * sin(TWO_PI*rand_num2) ) * sigma);
        if (new_value < 0)
            pi->f[pos] = 0;
        else if (new_value > pi->scale)
            pi->f[pos] = 255;
        else
            pi->f[pos] = new_value;
    }
}
