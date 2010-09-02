#ifndef __JAMES_PORTABLE_IMAGE
#define __JAMES_PORTABLE_IMAGE
#include "imgtypes.h"
#include "FFT.h"

int compare (const void * a, const void * b);
void thresholdImage(struct portImage *pi);
double vectorMag(struct vector *v);
struct vector vectorInit(struct portImage *pi, int point);
void insertionSort(struct pair *a, int length, double value, int p, int q);
struct hsv RGBtoHSV(struct rgb p);
struct rgb HSVtoRGB(struct hsv x);
void HSV(struct portImage *pi);
void RGB(struct portImage *pi);
void setBrightness(struct portImage *pi, float brightness);
void setSaturation(struct portImage *pi, float sat);
void simpleMeanFilter(struct portImage *pi, int window);
void vectorOrderStatistic(struct portImage *pi, int window, int value);
void simpleMedianFilter(struct portImage *pi, int window);
int imageDifferentPixels(struct portImage *pa, struct portImage *pb);
double imageAverageError(struct portImage *pa, struct portImage *pb);
double imageCompare(struct portImage *pa, struct portImage *pb);
void componentOrderStatistic(struct portImage *pi, int window, int value);
void componentMedianFilter(struct portImage *pi, int window);
int vectorMedianSinglePixel(struct portImage *pi, int i, int j, int window, struct pair *neighbors);
void vectorMedianOrderStatistic(struct portImage *pi, int window, int m, int replaceWithM);
void vectorMedianFilter(struct portImage *pi, int window);
double randDouble();
void addSaltPepperNoise(struct portImage *pi, double percentCorrupt);
void addGaussianNoise(struct portImage *pi, double percentCorrupt, double sigma);
int vectorSpacialSinglePixel(struct portImage *pi, int i, int j, int window, struct pair *neighbors);
void vectorSpacialOrderStatistic(struct portImage *pi, int window, int m, int replaceWithM);
void spacialFilter(struct portImage *pi);
void HSV_ValueFilter(struct portImage *pi, int window);
void spacialEdgeDetect(struct portImage *pi, int window);
int findBestThreshold(struct portImage *pi, int window, double noiseCorruption);
double calculatePotentialRMSE_mean(struct portImage *pi, struct portImage *noise, int window);
double calculatePotentialRMSE_smf(struct portImage *pi, struct portImage *noise, int window);
double calculatePotentialRMSE_cmf(struct portImage *pi, struct portImage *noise, int window);
void calculatePotentialRMSE_vmf(struct portImage *pi, struct portImage *noise, int window, double *scores);
void calculatePotentialRMSE(struct portImage *pi, struct portImage *noise, int window, double *scores);
void calculateAverageDifference(struct portImage *pi, int window, double noiseCorruption, double *scores);
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
