#ifndef __JAMES_PORTABLE_IMAGE
#define __JAMES_PORTABLE_IMAGE
#include "FFT.h"

#define MAX_LINE    100
#define TWO_PI      6.28318531

struct portImage {
    int height;       // Image Height
    int width;        // Image Width
    int size;         // Width * Height * Color Scales (so I don't have to keep recalculating it.)
    int scale;        // Max Color Scale value (usually 255 for full color or full grey scale images)
    int color;        // Boolean. 1 = COLOR IMAGE, 0 = GREY SCALE IMAGE
    int type;         // P-type of the image. More on this here: http://en.wikipedia.org/wiki/Portable_pixmap
    unsigned char *f; // Stores the raw RED, GREEN, BLUE scalers for each pixel in image. Very traditional.
    Complex       *c; // Stores the complex values after an FFT of the image is compute. Only used when FFT operations are needed.
    struct hsv    *p; // Perceptual Color Space values for each pixel. HSV stands for Hue, Saturation, and Value.
};

struct hsv {
    double h; /* Hue */
    double s; /* Saturation */
    double v; /* Value */
};

struct rgb {
    int r;
    int g;
    int b;
};

struct vector {
    double r;   /* Red */
    double g;   /* Green */
    double b;   /* Blue */
    double ur;  /* Unit Red */
    double ug;  /* Unit Green */
    double ub;  /* Unit Blue */
    double mag; /* Vector Magnitude */
};

struct pair {
    double v; /* Vector Length */
    int i; /* Horizontal Direction */
    int j; /* Vertical Direction */
};

void freeImage(struct portImage *pi);
int compare (const void * a, const void * b);
double logMag(Complex p);
struct portImage* copyImage(struct portImage *pi);
void thresholdImage(struct portImage *pi);
void minkowskiDivision(struct portImage *pi_n, struct portImage *pi_d);
void minkowskiAddition(struct portImage *pi, int maskwidth, int maskheight);
void minkowskiSubtraction(struct portImage *pi, int maskwidth, int maskheight);
void minkowskiOpening(struct portImage *pi, int maskwidth, int maskheight);
void minkowskiClosing(struct portImage *pi, int maskwidth, int maskheight);
double vectorMag(struct vector *v);
struct vector vectorInit(struct portImage *pi, int point);
void insertionSort(struct pair *a, int length, double value, int p, int q);
void FFT2D(struct portImage *pi);
void IFFT2D(struct portImage *pi);
void writeImage(struct portImage *pi, FILE* out);
void writeComplex(struct portImage *pi, FILE* out);
void graph_fftlogplot(struct portImage *pi, char *filename);
void lowpass(struct portImage *pi, double percent);
void highpass(struct portImage *pi, double percent);
double fftSimilarity(struct portImage* a, struct portImage* b, double percent);
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
struct portImage* readImage(FILE* fid);
void histogram(struct portImage *pi, int *hist);
void cumulative_distribution(int scales, int *hist, int *cdf);
void graph_histogram(struct portImage *pi, char *filename);
void graph_cdf(struct portImage *pi, char *filename);
void level_slice(struct portImage *pi, int level);
void equalize(struct portImage *pi);
void contrast_stretching(struct portImage *pi, float low_bound, float high_bound);

#endif
