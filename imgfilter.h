#ifndef __JAMES_IMAGE_FILTERS
#define __JAMES_IMAGE_FILTERS
#include "imgtypes.h"

int compare (const void * a, const void * b);
double vectorMag(struct vector *v);
struct vector vectorInit(struct portImage *pi, int point);
void insertionSort(struct pair *a, int length, double value, int p, int q);
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
int vectorSpacialSinglePixel(struct portImage *pi, int i, int j, int window, struct pair *neighbors);
void vectorSpacialOrderStatistic(struct portImage *pi, int window, int m, int replaceWithM);
void spacialFilter(struct portImage *pi);
int findBestThreshold(struct portImage *pi, int window, double noiseCorruption);
double calculatePotentialRMSE_mean(struct portImage *pi, struct portImage *noise, int window);
double calculatePotentialRMSE_smf(struct portImage *pi, struct portImage *noise, int window);
double calculatePotentialRMSE_cmf(struct portImage *pi, struct portImage *noise, int window);
void calculatePotentialRMSE_vmf(struct portImage *pi, struct portImage *noise, int window, double *scores);
void calculatePotentialRMSE(struct portImage *pi, struct portImage *noise, int window, double *scores);
void calculateAverageDifference(struct portImage *pi, int window, double noiseCorruption, double *scores);

#endif
