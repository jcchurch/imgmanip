#ifndef __JAMES_COLORSPACE
#define __JAMES_COLORSPACE
#include "imgtypes.h"
#include "colorspace.h"

struct hsv RGBtoHSV(struct rgb p);
struct rgb HSVtoRGB(struct hsv x);
void HSV(struct portImage *pi);
void RGB(struct portImage *pi);
void setBrightness(struct portImage *pi, float brightness);
void setSaturation(struct portImage *pi, float sat);
void HSV_ValueFilter(struct portImage *pi, int window);

#endif
