#ifndef __IMGMANIP_IMAGE_SIMPLE
#define __IMGMANIP_IMAGE_SIMPLE
#include "imgtypes.h"

void scale_reduce(struct portImage *pi, int newscale);
void spacial_reduce(struct portImage *pi, int resWidth, int resHeigth);
void histogram(struct portImage *pi, int *hist);
void cumulative_distribution(int scales, int *hist, int *cdf);
void graph_histogram(struct portImage *pi, char *filename);
void graph_cdf(struct portImage *pi, char *filename);
void level_slice(struct portImage *pi, int level);
void equalize(struct portImage *pi);
void contrast_stretching(struct portImage *pi, float low_bound, float high_bound);

#endif
