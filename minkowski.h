#ifndef __IMGMANIP_MINKOWSKI
#define __IMGMANIP_MINKOWSKI
#include "imgtypes.h"

void thresholdImage(struct portImage *pi);
void minkowskiDivision(struct portImage *pi_n, struct portImage *pi_d);
void minkowskiAddition(struct portImage *pi, int maskwidth, int maskheight);
void minkowskiSubtraction(struct portImage *pi, int maskwidth, int maskheight);
void minkowskiOpening(struct portImage *pi, int maskwidth, int maskheight);
void minkowskiClosing(struct portImage *pi, int maskwidth, int maskheight);

#endif
