#include "FFT.c"
#include <math.h>
#include <stdio.h>

#define PERIOD  256.0
#define TWO_PI  6.28318531
#define SAMPLES 8

double f(double x) { return 1 + cos( (TWO_PI * x)/PERIOD ) + cos( (TWO_PI * x * 2.0)/PERIOD ) +  cos( (TWO_PI * x * 3.0)/PERIOD ); }

int main(int argc, char **argv) {

    Complex samples[SAMPLES];
    int i, j;
    double x = 0.0;

    for (i = 0; i < SAMPLES; i++) {
        samples[i].real = f(x);
        samples[i].imag = 0.0;
        x += PERIOD/ (double)SAMPLES;
    }

    FFT(samples, SAMPLES);

    for (i = 0; i < SAMPLES; i++) {
        samples[i].real /= (double)SAMPLES;
        samples[i].imag /= (double)SAMPLES;
    }

    for (i = 0; i < 256; i++) {
        double spatial = 0.0;
        int freq = -4;
        for (j = 4; j < 4+SAMPLES; j++) {
            double theta = (TWO_PI * (double)i * (double)freq) / (PERIOD);
            double real = cos(theta);
            double imag = sin(theta);
            spatial += samples[j%8].real*real - samples[j%8].imag*imag;
            freq++;
        }
        printf("%3.3f,\n", spatial);
    }

    return 0;
} // End Main
