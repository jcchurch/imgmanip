#ifndef __IMGMANIP_NOISE
#define __IMGMANIP_NOISE
#include "imgtypes.h"

double randDouble();
void addSaltPepperNoise(struct portImage *pi, double percentCorrupt);
void addGaussianNoise(struct portImage *pi, double percentCorrupt, double sigma);

#endif
