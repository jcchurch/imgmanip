#ifndef __IMGMANIP_FFT2D
#define __IMGMANIP_FFT2D
#include "imgtypes.h"
#include "FFT.h"

double logMag(struct Complex p);
void writeComplex(struct portImage *pi, FILE* out);
void FFT2D(struct portImage *pi);
void IFFT2D(struct portImage *pi);
void graph_fftlogplot(struct portImage *pi, char *filename);
void lowpass(struct portImage *pi, double percent);
void highpass(struct portImage *pi, double percent);
double fftSimilarity(struct portImage* a, struct portImage* b, double percent);

#endif
