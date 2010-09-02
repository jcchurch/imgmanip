#ifndef __JAMES_IMAGE_TYPES
#define __JAMES_IMAGE_TYPES

#define MAX_LINE    100
#define TWO_PI      6.28318531

enum Colorspace { 
    csGREYSCALE,
    csRGB,
    csHSV
};

struct Complex { double real, imag; };

struct portImage {
    int height;               // Image Height
    int width;                // Image Width
    int size;                 // Width * Height * Color Scales (so I don't have to keep recalculating it.)
    int scale;                // Max Color Scale value (usually 255 for full color or full grey scale images)
    enum Colorspace color;    // Boolean. 1 = COLOR IMAGE, 0 = GREY SCALE IMAGE
    int type;                 // P-type of the image. More on this here: http://en.wikipedia.org/wiki/Portable_pixmap
    unsigned char  *f;        // Stores the raw RED, GREEN, BLUE scalers for each pixel in image. Very traditional.
    struct Complex *c;        // Stores the complex values after an FFT of the image is compute. Only used when FFT operations are needed.
    struct hsv     *p;        // Perceptual Color Space values for each pixel. HSV stands for Hue, Saturation, and Value.
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

#endif
