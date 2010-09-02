#ifndef __JAMES_NOISE
#define __JAMES_NOISE
#include "imgtypes.h"

double randDouble();
void addSaltPepperNoise(struct portImage *pi, double percentCorrupt);
void addGaussianNoise(struct portImage *pi, double percentCorrupt, double sigma);

#endif
