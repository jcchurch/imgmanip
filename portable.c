#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include "imgio.h"
#include "portable.h"

int compare (const void * a, const void * b) {
    double d = (*(struct pair*)a).v - (*(struct pair*)b).v;
    if      (d > 0.0)
        return 1;
    else if (d < 0.0)
        return -1;

    return 0;
}

void thresholdImage(struct portImage *pi) {
    int half = pi->scale >> 1;
    int i;
    for (i = 0; i < pi->size; i++)
        if (pi->f[i] > half)
            pi->f[i] = pi->scale;
        else
            pi->f[i] = 0;
}

double vectorMag(struct vector *v) {
    v->mag = sqrt( v->r*v->r + v->g*v->g + v->b*v->b );
    if (v->mag != 0.0) {
        v->ur  = v->r / v->mag;
        v->ug  = v->g / v->mag;
        v->ub  = v->b / v->mag;
    }
    else {
        v->ur  = 0.0;
        v->ug  = 0.0;
        v->ub  = 0.0;
    }

    return v->mag;
}

struct vector vectorInit(struct portImage *pi, int point) {
    struct vector v;
    v.r = pi->f[point+0];
    v.g = pi->f[point+1];
    v.b = pi->f[point+2];

    vectorMag(&v);
    return v;
}

void insertionSort(struct pair *a, int length, double value, int p, int q) {
    struct pair x = {value,p,q};
    int i = length - 1;
    while (i >= 0 && a[i].v > value) {
        a[i+1] = a[i];
        i--;
    }
    a[i+1]=x;
}

/** RGBtoHSV
 * Given a rgb data structure with all values in the range [0..255],
 * returns an hsv data structure with the Hue in the range of [0..300], Saturation [0..1], and Value [0..1]
 */
struct hsv RGBtoHSV(struct rgb p) {
    double dr = (double)p.r / 255.0;
    double dg = (double)p.g / 255.0;
    double db = (double)p.b / 255.0;
    double max = dr;
    double min = dr;
    double delta;
    struct hsv x;
    if (dg > max) max = dg;
    if (db > max) max = db;
    if (dg < min) min = dg;
    if (db < min) min = db;
    x.v = max;
    delta = max - min;

    if (max != 0) {
        x.s = delta / max;       /* s */
    }
    else {
        /* r = g = b = 0        // s = 0, v is undefined */
        x.s = 0;
        x.h = 0;
        return x;
    }

    if (x.s != 0) {
        if( dr == max )
            x.h = ( dg - db ) / delta;     /* between yellow & magenta */
        else if( dg == max )
            x.h = 2.0 + ( db - dr ) / delta; /* between cyan & yellow */
        else
            x.h = 4.0 + ( dr - dg ) / delta; /* between magenta & cyan */
    }
    else {
        x.h = 0;
    }

    x.h *= 60;               /* degrees */
    if( x.h < 0 )
        x.h += 360;

    return x;
}

struct rgb HSVtoRGB(struct hsv x) {
    int i;
    struct rgb y;
    double f, p, q, t, v, h;

    if( x.s == 0 ) {
        /* achromatic (grey) */
        y.r = y.g = y.b = (int) (x.v * 255.0);
        return y;
    }

    h =  x.h / 60.0;      /* sector 0 to 5 */
    i = (int) h;
    f = h - i;            /* factorial part of h */
    p = x.v * ( 1 - x.s );
    q = x.v * ( 1 - x.s * f );
    t = x.v * ( 1 - x.s * ( 1 - f ) );

    v  = x.v * 255.0;
    t *= 255.0;
    q *= 255.0;
    p *= 255.0;

    switch( i ) {
        case 0:
            y.r = v;
            y.g = t;
            y.b = p;
            break;
        case 1:
            y.r = q;
            y.g = v;
            y.b = p;
            break;
        case 2:
            y.r = p;
            y.g = v;
            y.b = t;
            break;
        case 3:
            y.r = p;
            y.g = q;
            y.b = v;
            break;
        case 4:
            y.r = t;
            y.g = p;
            y.b = v;
            break;
        default:        /* case 5: */
            y.r = v;
            y.g = p;
            y.b = q;
            break;
    }

    return y;
}

void RGB(struct portImage *pi) {
    int i;
    int pixels = pi->width * pi->height;
    struct rgb c;

    if (pi->color != csHSV)
        return;

    if (!pi->p)
        return;

    for (i = 0; i < pixels; i++) {
        c = HSVtoRGB(pi->p[i]);
        pi->f[i * 3 + 0] = c.r;
        pi->f[i * 3 + 1] = c.g;
        pi->f[i * 3 + 2] = c.b;
    }
}

void HSV(struct portImage *pi) {
    int i;
    int pixels = pi->width * pi->height;
    struct rgb c;

    if (pi->color != csRGB)
        return;

    if (!pi->p) {
        pi->p = (struct hsv *)malloc(pi->size * sizeof(struct hsv));
    }

    for (i = 0; i < pixels; i++) {
        c.r = pi->f[i * 3 + 0];
        c.g = pi->f[i * 3 + 1];
        c.b = pi->f[i * 3 + 2];
        pi->p[i] = RGBtoHSV(c);
    }
}

void setBrightness(struct portImage *pi, float brightness) {
    int i;
    HSV(pi);
    for (i = 0; i < pi->size; i++)
        pi->p[i].v = (double) brightness;
    RGB(pi);
}

void setSaturation(struct portImage *pi, float sat) {
    int i;
    HSV(pi);
    for (i = 0; i < pi->size; i++)
        pi->p[i].s = (double) sat;
    RGB(pi);
}

void simpleMeanFilter(struct portImage *pi, int window) {
    int            i,j,k,l,n;
    int            w    = pi->width;
    unsigned char *buffer;
    int            window_screen = window*window;
    struct pair   *neighbors;
    int            bound = window >> 1;
    int            index_pixel, index_pixel_window;

    if (pi->color != csGREYSCALE)
        w *= 3;

    neighbors = (struct pair*)malloc( window_screen * sizeof(struct pair) );
    buffer = (unsigned char*)malloc(pi->size * sizeof(unsigned char));

    if (!buffer || !neighbors) {
        fprintf(stderr, "ERROR: Cannot allocate memory for mean filter.\n");
        return;
    }

    if (pi->color == csGREYSCALE) {
        for (j = 0; j < pi->height; j++) {
            for (i = 0; i < pi->width; i++) {
                if (j < bound || j+bound >= pi->height || i < bound || i+bound >= pi->width)
                    buffer[j*w+i] = pi->f[j*w+i];
                else {
                    double sum = 0.0;
                    n = 0;
                    for (k = -bound; k <= bound; k++) {
                        for (l = -bound; l <= bound; l++) {
                            sum += (double)pi->f[(j+k)*w+(i+l)];
                        }
                    }
                    buffer[j*w+i] = sum / (double)window_screen;
                }
            }
        }
    }
    else {
        for (j = 0; j < pi->height; j++) {
            for (i = 0; i < pi->width; i++) {
                index_pixel = j*w+(i*3);
                if (j < bound || j+bound >= pi->height || i < bound || i+bound >= pi->width) {
                    buffer[index_pixel+0] = pi->f[index_pixel+0];
                    buffer[index_pixel+1] = pi->f[index_pixel+1];
                    buffer[index_pixel+2] = pi->f[index_pixel+2];
                }
                else {
                    double sum_red = 0.0;
                    double sum_green = 0.0;
                    double sum_blue = 0.0;
                    n = 0;
                    for (k = -bound; k <= bound; k++) {
                        for (l = -bound; l <= bound; l++) {
                            index_pixel_window = (j+k)*w+(i+l)*3;
                            sum_red   += (double) pi->f[index_pixel_window+0];
                            sum_green += (double) pi->f[index_pixel_window+1];
                            sum_blue  += (double) pi->f[index_pixel_window+2];
                        }
                    }
                    buffer[index_pixel+0] = sum_red / (double)window_screen;
                    buffer[index_pixel+1] = sum_green / (double)window_screen;
                    buffer[index_pixel+2] = sum_blue / (double)window_screen;
                }
            }
        }
    }

    for (i = 0; i < pi->size; i++)
        pi->f[i] = buffer[i];

    free(neighbors);
    free(buffer);
}

void vectorOrderStatistic(struct portImage *pi, int window, int value) {
    int            i,j,k,l,n;
    int            w    = pi->width;
    unsigned char *buffer;
    struct vector        *v;
    struct pair          *neighbors;
    int            window_screen = window*window;
    int            bound = window >> 1;

    if (value >= window_screen)
        value = window_screen - 1;

    if (pi->color != csGREYSCALE)
        w *= 3;

    neighbors = (struct pair*)         malloc(window_screen * sizeof(struct pair));
    buffer    = (unsigned char*)malloc(pi->size * sizeof(unsigned char));
    v         = (struct vector*)       malloc(pi->height * pi->width * sizeof(struct vector));

    if (!buffer || !v || !neighbors) {
        fprintf(stderr, "ERROR: Cannot allocate memory for median filter.\n");
        return;
    }

    if (pi->color == csGREYSCALE) {
        for (j = 0; j < pi->height; j++) {
            for (i = 0; i < pi->width; i++) {
                if (j < bound || j+bound >= pi->height || i < bound || i+bound >= pi->width)
                    buffer[j*w+i] = pi->f[j*w+i];
                else {
                    n = 0;
                    for (k = j-bound; k <= j+bound; k++) {
                        for (l = i-bound; l <= i+bound; l++) {
                            insertionSort(neighbors, n++, (double)pi->f[k*w+l], l, k);
                        }
                    }
                    buffer[j*w+i] = (unsigned char)neighbors[value].v;
                }
            }
        }
    }
    else {
        for (j = 0; j < pi->height; j++)
            for (i = 0; i < pi->width; i++)
                v[j*pi->width+i] = vectorInit(pi, j*w+i*3);

        for (j = 0; j < pi->height; j++) {
            for (i = 0; i < pi->width; i++) {
                if (j < bound || j+bound >= pi->height || i < bound || i+bound >= pi->width) {
                    buffer[j*w+(i*3)+0] = pi->f[j*w+(i*3)+0];
                    buffer[j*w+(i*3)+1] = pi->f[j*w+(i*3)+1];
                    buffer[j*w+(i*3)+2] = pi->f[j*w+(i*3)+2];
                }
                else {
                    n = 0;
                    for (k = j-bound; k <= j+bound; k++) {
                        for (l = i-bound; l <= i+bound; l++) {
                            insertionSort(neighbors, n++, v[k*pi->width+l].mag, l, k);
                        }
                    }
                    buffer[j*w+(i*3)+0] = pi->f[neighbors[value].j*w+(neighbors[value].i*3)+0];
                    buffer[j*w+(i*3)+1] = pi->f[neighbors[value].j*w+(neighbors[value].i*3)+1];
                    buffer[j*w+(i*3)+2] = pi->f[neighbors[value].j*w+(neighbors[value].i*3)+2];
                }
            }
        }
    }

    for (i = 0; i < pi->size; i++)
        pi->f[i] = buffer[i];

    free(neighbors);
    free(buffer);
    free(v);
}

void simpleMedianFilter(struct portImage *pi, int window) {
    vectorOrderStatistic(pi, window, (window*window)/2);
}

int imageDifferentPixels(struct portImage *pa, struct portImage *pb) {
    int i;
    int count = 0;
    if (pa->color != pb->color && pa->height != pb->height && pa->width != pb->width)
        return -1;

    if (pa->color != csGREYSCALE) {
        for (i = 0; i < pa->height * pa->width; i++)
            if (pa->f[i*3+0] != pb->f[i*3+0] || pa->f[i*3+1] != pb->f[i*3+1] || pa->f[i*3+2] != pb->f[i*3+2])
                count++;
    }
    else {
        for (i = 0; i < pa->size; i++)
            if (pa->f[i] != pb->f[i])
                count++;
    }

    return count;
}

double imageAverageError(struct portImage *pa, struct portImage *pb) {
    int i;
    double score = 0.0;
    if (pa->color != pb->color && pa->height != pb->height && pa->width != pb->width)
        return -1.0;

    if (pa->color != csGREYSCALE) {
        for (i = 0; i < pa->size; i+=3) {
            score += sqrt((pa->f[i+0]-pb->f[i+0])*(pa->f[i+0]-pb->f[i+0])+
                          (pa->f[i+1]-pb->f[i+1])*(pa->f[i+1]-pb->f[i+1])+
                          (pa->f[i+2]-pb->f[i+2])*(pa->f[i+2]-pb->f[i+2]));

        }
    }
    else {
        for (i = 0; i < pa->size; i++)
            score += (double) abs(pa->f[i]-pb->f[i]);
    }

    return score / (double)(pa->height * pa->width);
}

double imageCompare(struct portImage *pa, struct portImage *pb) {
    int i;
    double score = 0.0;
    if (pa->color != pb->color && pa->height != pb->height && pa->width != pb->width)
        return -1.0;

    if (pa->color != csGREYSCALE) {
        for (i = 0; i < pa->size; i+=3) {
            score += ((pa->f[i+0]-pb->f[i+0])*(pa->f[i+0]-pb->f[i+0])+
                      (pa->f[i+1]-pb->f[i+1])*(pa->f[i+1]-pb->f[i+1])+
                      (pa->f[i+2]-pb->f[i+2])*(pa->f[i+2]-pb->f[i+2]));
            
        }
    }
    else {
        for (i = 0; i < pa->size; i++)
            score += (double) abs(pa->f[i]-pb->f[i]);
    }

    return sqrt( score / (double)(pa->height * pa->width) );
}

void componentOrderStatistic(struct portImage *pi, int window, int value) {
    int            i,j,k,l,n,m, point;
    int            w = pi->width;
    unsigned char *buffer;
    int            window_screen = window*window;
    int            bound = window >> 1;
    struct pair          *neighbors;

    if (pi->color != csGREYSCALE)
        w *= 3;

    neighbors = (struct pair*)malloc(window * window * sizeof(struct pair));
    buffer = (unsigned char*)malloc(pi->size * sizeof(unsigned char));

    if (!buffer || !neighbors) {
        fprintf(stderr, "ERROR: Cannot allocate memory for spacial filter.\n");
        return;
    }

    if (pi->color != csGREYSCALE) {
        for (j = 0; j < pi->height; j++) {
            for (i = 0; i < pi->width; i++) {
                point = j*w+(i*3);
                if (j < bound || j+bound >= pi->height || i < bound || i+bound >= pi->width) {
                    buffer[point+0] = pi->f[point+0];
                    buffer[point+1] = pi->f[point+1];
                    buffer[point+2] = pi->f[point+2];
                }
                else {
                    for (m = 0; m < 3; m++) {

                        n = 0;
                        for (k = j-bound; k <= j+bound; k++)
                            for (l = i-bound; l <= i+bound; l++)
                                insertionSort(neighbors, n++, (double) pi->f[k*w+(l*3)+m], l, k);

                        buffer[point+m] = (unsigned char)neighbors[value].v;
                    }
                }
            }
        }
    }

    for (i = 0; i < pi->size; i++)
        pi->f[i] = buffer[i];

    free(buffer);
    free(neighbors);
}

void componentMedianFilter(struct portImage *pi, int window) {
    componentOrderStatistic(pi, window, (window*window)/2);
}

int vectorMedianSinglePixel(struct portImage *pi, int i, int j, int window, struct pair *neighbors) {
    int n, p, q, k, l;
    int a, b;
    double diff_sum = 0.0;
    int w = pi->width * 3;
    int bound = window >> 1;
    int lowbound, rightbound;

    if (i-bound < 0 || i+bound >= pi->width || j-bound < 0 || j+bound >= pi->height)
        return 0;

    n = 0;
    lowbound = j + bound;
    rightbound = i + bound;
    for (q = j-bound; q <= lowbound; q++) {
        for (p = i-bound; p <= rightbound; p++) {
            diff_sum = 0.0;

            a = q*w+(p*3);
            for (k = j-bound; k <= j+bound; k++) {
                for (l = i-bound; l <= i+bound; l++) {

                    if (k == q && l == p)
                        continue;

                    b = k*w+(l*3);
                    diff_sum += (pi->f[a+0]-pi->f[b+0])*(pi->f[a+0]-pi->f[b+0]) +
                                (pi->f[a+1]-pi->f[b+1])*(pi->f[a+1]-pi->f[b+2]) +
                                (pi->f[a+2]-pi->f[b+2])*(pi->f[a+2]-pi->f[b+2]);
                }
            }

            neighbors[n].v = diff_sum;
            neighbors[n].i = p;
            neighbors[n].j = q;
            n++;
        }
    }

    qsort(neighbors, n, sizeof(struct pair), compare);
    return n;
}

void vectorMedianOrderStatistic(struct portImage *pi, int window, int m, int replaceWithM) {
    int            i,j,k,l,p,q,n;
    int            a, b;
    int            rc, gc, bc, ro, go, bo;
    int            w    = pi->width;
    unsigned char *buffer;
    struct pair   *neighbors;
    int            window_screen = window*window;
    int            bound = window >> 1;

    if (pi->color != csGREYSCALE)
        w *= 3;

    neighbors = (struct pair*)malloc(window_screen * sizeof(struct pair));
    buffer = (unsigned char*)malloc(pi->size * sizeof(unsigned char));

    if (!buffer || !neighbors) {
        fprintf(stderr, "ERROR: Cannot allocate memory for spacial filter.\n");
        return;
    }

    if (pi->color != csGREYSCALE) {
        for (j = 0; j < pi->height; j++) {
            for (i = 0; i < pi->width; i++) {

                n = vectorMedianSinglePixel(pi, i, j, window, neighbors);

                if (n == 0) {
                    buffer[j*w+(i*3)+0] = pi->f[j*w+(i*3)+0];
                    buffer[j*w+(i*3)+1] = pi->f[j*w+(i*3)+1];
                    buffer[j*w+(i*3)+2] = pi->f[j*w+(i*3)+2];
                }
                else {

                    /* Search for i and j in neighbors listing */
                    for (p = 0; p < n; p++)
                        if (neighbors[p].j == j && neighbors[p].i == i)
                            break;

                    if (replaceWithM) {
                        buffer[j*w+(i*3)+0] = pi->f[neighbors[m].j*w+(neighbors[m].i*3)+0];
                        buffer[j*w+(i*3)+1] = pi->f[neighbors[m].j*w+(neighbors[m].i*3)+1];
                        buffer[j*w+(i*3)+2] = pi->f[neighbors[m].j*w+(neighbors[m].i*3)+2];
                    }
                    else {
                        if (p <= m) {
                            buffer[j*w+(i*3)+0] = pi->f[j*w+(i*3)+0];
                            buffer[j*w+(i*3)+1] = pi->f[j*w+(i*3)+1];
                            buffer[j*w+(i*3)+2] = pi->f[j*w+(i*3)+2];
                        }
                        else {
                            buffer[j*w+(i*3)+0] = pi->f[neighbors[0].j*w+(neighbors[0].i*3)+0];
                            buffer[j*w+(i*3)+1] = pi->f[neighbors[0].j*w+(neighbors[0].i*3)+1];
                            buffer[j*w+(i*3)+2] = pi->f[neighbors[0].j*w+(neighbors[0].i*3)+2];
                        }
                    }
                }
            }
        }
    }

    for (i = 0; i < pi->size; i++)
        pi->f[i] = buffer[i];

    free(neighbors);
    free(buffer);
}

void vectorMedianFilter(struct portImage *pi, int window) {
    vectorMedianOrderStatistic(pi, window, 0, 1);
}

void HSV_ValueFilter(struct portImage *pi, int window) {
    int            i,j,k,l,n, point, better_point;
    int            value = (window * window) / 2;
    int            w = pi->width;
    unsigned char *buffer;
    int            window_screen = window*window;
    int            bound = window >> 1;
    struct pair   *neighbors;

    if (pi->color != csGREYSCALE)
        w *= 3;

    HSV(pi);

    neighbors = (struct pair*)malloc(window * window * sizeof(struct pair));
    buffer = (unsigned char*)malloc(pi->size * sizeof(unsigned char));

    if (!buffer || !neighbors) {
        fprintf(stderr, "ERROR: Cannot allocate memory for spacial filter.\n");
        return;
    }

    if (pi->color != csGREYSCALE) {
        for (j = 0; j < pi->height; j++) {
            for (i = 0; i < pi->width; i++) {
                point = j*w+(i*3);
                if (j < bound || j+bound >= pi->height || i < bound || i+bound >= pi->width) {
                    buffer[point+0] = pi->f[point+0];
                    buffer[point+1] = pi->f[point+1];
                    buffer[point+2] = pi->f[point+2];
                }
                else {
                    n = 0;
                    for (k = j-bound; k <= j+bound; k++)
                        for (l = i-bound; l <= i+bound; l++)
                            insertionSort(neighbors, n++, (double) pi->p[k*w+(l*3)].v, l, k);

                    better_point = neighbors[value].j*w+(3 * neighbors[value].i);
                    buffer[point+0] = pi->f[better_point+0];
                    buffer[point+1] = pi->f[better_point+1];
                    buffer[point+2] = pi->f[better_point+2];
                }
            }
        }
    }

    for (i = 0; i < pi->size; i++)
        pi->f[i] = buffer[i];

    free(buffer);
    free(neighbors);
}

double randDouble() {
    return rand() / (double) RAND_MAX;
}

void addSaltPepperNoise(struct portImage *pi, double percentCorrupt) {
    int i;
    int x;
    int noise;

    for (i = 0; i < pi->width * pi->height; i++) {
        if (randDouble() < percentCorrupt) {
            x     = (int)(randDouble() * 3.0);                     /* Random number from 0 to 2 */
            noise = (int)(randDouble() * (double)(pi->scale + 1)); /* 0 to 255 */
            pi->f[i*3+x] = noise;
        }
    }
}

/* This is a bad noise model. Do not use it. */
void addGaussianNoise(struct portImage *pi, double percentCorrupt, double sigma) {
    int i;
    int pixelsToCorrupt = (int) ((double)pi->height * pi->width * percentCorrupt);
    double rand_num1, rand_num2;
    int new_value, pos;

    for (i = 0; i < pixelsToCorrupt; i++) {
        rand_num1 = randDouble();
        rand_num2 = randDouble();
        pos = rand() % pi->size;
        new_value = pi->f[ pos ] + (int) (sqrt( -2*log(rand_num1) * sin(TWO_PI*rand_num2) ) * sigma);
        if (new_value < 0)
            pi->f[pos] = 0;
        else if (new_value > pi->scale)
            pi->f[pos] = 255;
        else
            pi->f[pos] = new_value;
    }
}

int vectorSpacialSinglePixel(struct portImage *pi, int i, int j, int window, struct pair *neighbors) {
    int n, p, q, k, l;
    int rc, gc, bc, ro, go, bo;
    int a, b;
    struct vector diff_vector;
    double diff_mag = 0.0;
    double diff_mag_denom = 0.0;
    int w = pi->width * 3;
    int bound = window >> 1;
    int lowbound, rightbound;

    if (i-bound < 0 || i+bound >= pi->width || j-bound < 0 || j+bound >= pi->height)
        return 0;

    n = 0;
    lowbound = j + bound;
    rightbound = i + bound;
    for (q = j-bound; q <= lowbound; q++) {
        for (p = i-bound; p <= rightbound; p++) {

            diff_vector.r = 0.0;
            diff_vector.g = 0.0;
            diff_vector.b = 0.0;

            a = q*w+(p*3);

            rc = pi->f[a];
            gc = pi->f[a+1];
            bc = pi->f[a+2];

            for (k = j-bound; k <= lowbound; k++) {
                for (l = i-bound; l <= rightbound; l++) {

                    if (k == q && l == p)
                        continue;

                    b = k*w+(l*3);
                    ro = rc - pi->f[b];
                    go = gc - pi->f[b+1];
                    bo = bc - pi->f[b+2];

                    if( ro || go || bo ) {
                        diff_mag_denom = sqrt( ro*ro + go*go + bo*bo );

                        diff_vector.r += ro / diff_mag_denom;
                        diff_vector.g += go / diff_mag_denom;
                        diff_vector.b += bo / diff_mag_denom;
                    }
                }
            }

            diff_mag = /*sqrt*/(diff_vector.r*diff_vector.r + diff_vector.g*diff_vector.g + diff_vector.b*diff_vector.b);

            neighbors[n].v = diff_mag;
            neighbors[n].i = p;
            neighbors[n].j = q;
            n++;
        }
    }

    qsort(neighbors, n, sizeof(struct pair), compare);
    return n;
}

void vectorSpacialOrderStatistic(struct portImage *pi, int window, int m, int replaceWithM) {
    int            i,j,n,p, index_pixel;
    int            w    = pi->width;
    unsigned char *buffer;
    struct pair          *neighbors;
    int            window_screen = window*window;

    if (pi->color != csGREYSCALE)
        w *= 3;

    if (m >= window_screen) {
        m = window_screen - 1;
    }

    neighbors = (struct pair*)malloc(window_screen * sizeof(struct pair));
    buffer = (unsigned char*)malloc(pi->size * sizeof(unsigned char));
    if (!buffer || !neighbors) {
        fprintf(stderr, "ERROR: Cannot allocate memory for spacial filter.\n");
        return;
    }

    if (pi->color != csGREYSCALE) {
        for (j = 0; j < pi->height; j++) {
            for (i = 0; i < pi->width; i++) {

                index_pixel = j*w+(i*3);

                n = vectorSpacialSinglePixel(pi, i, j, window, neighbors);

                /* Search for i and j in neighbors listing */
                for (p = 0; p < n; p++)
                    if (neighbors[p].j == j && neighbors[p].i == i)
                        break;

                if (n < 0) {
                    /* Center Pixel (i.e. Do not change pixel, for we are on edge of image.) */
                    buffer[index_pixel+0] = pi->f[index_pixel+0];
                    buffer[index_pixel+1] = pi->f[index_pixel+1];
                    buffer[index_pixel+2] = pi->f[index_pixel+2];
                    continue;
                }

                if (replaceWithM) {
                    /* Threshold Pixel */
                    buffer[index_pixel+0] = pi->f[neighbors[m].j*w+(neighbors[m].i*3)+0];
                    buffer[index_pixel+1] = pi->f[neighbors[m].j*w+(neighbors[m].i*3)+1];
                    buffer[index_pixel+2] = pi->f[neighbors[m].j*w+(neighbors[m].i*3)+2];
                }
                else {
                    if (p <= m) {
                        /* Center Pixel (i.e. Do not change pixel) */
                        buffer[index_pixel+0] = pi->f[index_pixel+0];
                        buffer[index_pixel+1] = pi->f[index_pixel+1];
                        buffer[index_pixel+2] = pi->f[index_pixel+2];
                    }

                    else {
                        /* Best Pixel */
                        buffer[index_pixel+0] = pi->f[neighbors[0].j*w+(neighbors[0].i*3)+0];
                        buffer[index_pixel+1] = pi->f[neighbors[0].j*w+(neighbors[0].i*3)+1];
                        buffer[index_pixel+2] = pi->f[neighbors[0].j*w+(neighbors[0].i*3)+2];
                    }
                }
            }
        }
    }

    for (i = 0; i < pi->size; i++)
        pi->f[i] = buffer[i];

    free(neighbors);
    free(buffer);
}

void spacialFilter(struct portImage *pi) {
    vectorSpacialOrderStatistic(pi, 3, 3, 0);
}

void spacialEdgeDetect(struct portImage *pi, int window) {
}

int findBestThreshold(struct portImage *pi, int window, double noiseCorruption) {
    struct portImage *noise = copyImage(pi);
    int window_screen = window * window;
    struct pair *neighbors;
    double *scores;
    double minValue = (double)(pi->scale * 2);
    int minThreshold = -1;
    int i, j, k, p, n, a, b;
    int w = 3*pi->width;

    neighbors = (struct pair*)  malloc( sizeof(struct pair)   * window_screen );
    scores    = (double*)malloc( sizeof(double) * window_screen );

    if (!neighbors || !scores) {
        fprintf(stderr, "Memory allocation failed when attempting to find best threshold.");
        return -1;
    }

    addSaltPepperNoise(noise, noiseCorruption);

    for (i = 0; i < window_screen; i++)
        scores[i] = 0.0;

    for (j = 0; j < pi->height; j++) {
        for (i = 0; i < pi->width; i++) {
            n = vectorSpacialSinglePixel(noise, i, j, window, neighbors);

            /* If n is 0, then p will also be 0. */
            for (p = 0; p < n; p++)
                if (neighbors[p].j == j && neighbors[p].i == i)
                    break;

            a = (j*w)+(i*3);
            b = neighbors[0].j*w+(neighbors[0].i*3);

            for (k = 0; k < window_screen; k++) {
                if (p <= k) {
                    /* Compute as if using Center Pixel */
                    scores[k] += ((pi->f[a  ]-noise->f[a  ])*(pi->f[a  ]-noise->f[a  ]) +
                                  (pi->f[a+1]-noise->f[a+1])*(pi->f[a+1]-noise->f[a+1]) +
                                  (pi->f[a+2]-noise->f[a+2])*(pi->f[a+2]-noise->f[a+2]));
                }
                else {
                    /* Compute as if using Best Pixel */
                    scores[k] += ((pi->f[a  ]-noise->f[b  ])*(pi->f[a  ]-noise->f[b  ]) +
                                  (pi->f[a+1]-noise->f[b+1])*(pi->f[a+1]-noise->f[b+1]) +
                                  (pi->f[a+2]-noise->f[b+2])*(pi->f[a+2]-noise->f[b+2]));
                }
            }
        }
    }

    minValue = scores[0];
    minThreshold = 0;
    for (i = 1; i < window_screen; i++) {
        if (scores[i] < minValue) {
            minValue = scores[i];
            minThreshold = i;
        }
    }

    free(neighbors);
    free(scores);
    freeImage(noise);
    return minThreshold;
}

double calculatePotentialRMSE_mean(struct portImage *pi, struct portImage *noise, int window) {
    int            i,j,k,l,n;
    int            w    = pi->width;
    int            window_screen = window*window;
    struct pair          *neighbors;
    int            bound = window >> 1;
    int            index_pixel, index_pixel_window;
    double         score = 0.0;

    if (pi->color != csGREYSCALE)
        w *= 3;

    neighbors = (struct pair*)malloc( window_screen * sizeof(struct pair) );

    if (!neighbors) {
        fprintf(stderr, "ERROR: Cannot allocate memory for mean filter.\n");
        return;
    }

    for (j = 0; j < pi->height; j++) {
        for (i = 0; i < pi->width; i++) {
            index_pixel = j*w+(i*3);
            if (j >= bound && j+bound < pi->height && i >= bound && i+bound < pi->width) {
                double sum_red = 0.0;
                double sum_green = 0.0;
                double sum_blue = 0.0;
                n = 0;
                for (k = -bound; k <= bound; k++) {
                    for (l = -bound; l <= bound; l++) {
                        index_pixel_window = (j+k)*w+(i+l)*3;
                        sum_red   += (double) noise->f[index_pixel_window+0];
                        sum_green += (double) noise->f[index_pixel_window+1];
                        sum_blue  += (double) noise->f[index_pixel_window+2];
                    }
                }
                sum_red /= (double)window_screen;
                sum_green /= (double)window_screen;
                sum_blue /= (double)window_screen;
                score += (pi->f[index_pixel  ]-sum_red  )*(pi->f[index_pixel  ]-sum_red  ) +
                         (pi->f[index_pixel+1]-sum_green)*(pi->f[index_pixel+1]-sum_green) +
                         (pi->f[index_pixel+2]-sum_blue )*(pi->f[index_pixel+2]-sum_blue );
            }
        }
    }

    free(neighbors);
    return sqrt( score / (double)(pi->height * pi->width) );
}

double calculatePotentialRMSE_smf(struct portImage *pi, struct portImage *noise, int window) {
    int            i,j,k,l,n, r,g,b, indexa, indexb;
    int            w    = pi->width;
    struct vector        *v;
    struct pair          *neighbors;
    int            window_screen = window*window;
    int            bound = window >> 1;
    int            value = window_screen >> 1;
    double         score;

    if (pi->color != csGREYSCALE)
        w *= 3;

    neighbors = (struct pair*)   malloc(window_screen * sizeof(struct pair));
    v         = (struct vector*) malloc(pi->height * pi->width * sizeof(struct vector));

    if (!v || !neighbors) {
        fprintf(stderr, "ERROR: Cannot allocate memory for median filter.\n");
        return;
    }

    for (j = 0; j < pi->height; j++)
        for (i = 0; i < pi->width; i++)
            v[j*pi->width+i] = vectorInit(noise, j*w+i*3);

    for (j = 0; j < pi->height; j++) {
        for (i = 0; i < pi->width; i++) {
            if (j >= bound && j+bound < pi->height && i >= bound && i+bound < pi->width) {
                n = 0;
                for (k = j-bound; k <= j+bound; k++) {
                    for (l = i-bound; l <= i+bound; l++) {
                        insertionSort(neighbors, n++, v[k*pi->width+l].mag, l, k);
                    }
                }
                indexa = j*w+(i*3);
                indexb = neighbors[value].j*w+(neighbors[value].i*3);

                r = pi->f[indexb+0];
                g = pi->f[indexb+1];
                b = pi->f[indexb+2];
                score += (pi->f[indexa+0]-r)*(pi->f[indexa+0]-r) +
                         (pi->f[indexa+1]-g)*(pi->f[indexa+1]-g) +
                         (pi->f[indexa+2]-b)*(pi->f[indexa+2]-b);
            }
        }
    }

    free(neighbors);
    free(v);
    return sqrt( score / (double)(pi->height * pi->width) );
}

double calculatePotentialRMSE_cmf(struct portImage *pi, struct portImage *noise, int window) {
    int            i,j,k,l,n,m, thisComponent;
    int            w = pi->width;
    struct pair          *neighbors;
    int            window_screen = window*window;
    int            value = window_screen >> 1;
    int            bound = window >> 1;
    double         score = 0;
    int            index;

    if (pi->color != csGREYSCALE)
        w *= 3;

    neighbors = (struct pair*)malloc(window_screen * sizeof(struct pair));

    if (!neighbors) {
        fprintf(stderr, "ERROR: Cannot allocate memory for spacial filter.\n");
        return;
    }

    if (pi->color != csGREYSCALE) {
        for (j = 0; j < pi->height; j++) {
            for (i = 0; i < pi->width; i++) {
                if (j >= bound && j+bound < pi->height && i >= bound && i+bound < pi->width) {
                    for (m = 0; m <= 2; m++) {
                        n = 0;
                        for (k = j-bound; k <= j+bound; k++) {
                            for (l = i-bound; l <= i+bound; l++) {
                                insertionSort(neighbors, n++, (double) noise->f[k*w+(l*3)+m], l, k);
                            }
                        }
                        index = j*w+(i*3)+m;
                        thisComponent = noise->f[neighbors[value].j*w+(neighbors[value].i*3)+m];
                        score += (pi->f[index]-thisComponent)*(pi->f[index]-thisComponent);
                    }
                }
            }
        }
    }

    free(neighbors);
    return sqrt( score / (double)(pi->height * pi->width) );
}

void calculatePotentialRMSE_vmf(struct portImage *pi, struct portImage *noise, int window, double *scores) {
    int window_screen = window * window;
    struct pair *neighbors;
    int i, j, p, n, a, b;
    int w = 3*pi->width;
    int threshold;

    neighbors = (struct pair*)malloc(sizeof(struct pair)*window_screen);
    if (!neighbors) {
        fprintf(stderr, "Memory allocation failed when calculating potential RMSE.\n");
        return;
    }

    for (threshold = 0; threshold < window_screen; threshold++)
        scores[threshold] = 0.0;

    for (j = 0; j < pi->height; j++) {
        for (i = 0; i < pi->width; i++) {
            n = vectorMedianSinglePixel(noise, i, j, window, neighbors);

            /* If n is 0, p will be 0 also */
            for (p = 0; p < n; p++)
                if (neighbors[p].j == j && neighbors[p].i == i)
                    break;

            a = (j*w)+(i*3);
            b = neighbors[0].j*w+(neighbors[0].i*3);

            for (threshold = 0; threshold < window_screen; threshold++) {
                if (p <= threshold) {
                    /* Calculate using Center Pixel */
                    scores[threshold] += ((pi->f[a  ]-noise->f[a  ])*(pi->f[a  ]-noise->f[a  ]) +
                                          (pi->f[a+1]-noise->f[a+1])*(pi->f[a+1]-noise->f[a+1]) +
                                          (pi->f[a+2]-noise->f[a+2])*(pi->f[a+2]-noise->f[a+2]));
                }
                else {
                    /* Calculate using Best Pixel */
                    scores[threshold] += ((pi->f[a  ]-noise->f[b  ])*(pi->f[a  ]-noise->f[b  ]) +
                                          (pi->f[a+1]-noise->f[b+1])*(pi->f[a+1]-noise->f[b+1]) +
                                          (pi->f[a+2]-noise->f[b+2])*(pi->f[a+2]-noise->f[b+2]));
                }
            }
        }
    }

    for (threshold = 0; threshold < window_screen; threshold++)
        scores[threshold] = sqrt( scores[threshold] / (double)(pi->height * pi->width) );

    free(neighbors);
}

void calculatePotentialRMSE(struct portImage *pi, struct portImage *noise, int window, double *scores) {
    int window_screen = window * window;
    struct pair *neighbors;
    int i, j, p, n, a, b;
    int w = 3*pi->width;
    int threshold;

    neighbors = (struct pair*)malloc(sizeof(struct pair)*window_screen);
    if (!neighbors) {
        fprintf(stderr, "Memory allocation failed when calculating potential RMSE.\n");
        return;
    }

    for (threshold = 0; threshold < window_screen; threshold++)
        scores[threshold] = 0.0;

    for (j = 0; j < pi->height; j++) {
        for (i = 0; i < pi->width; i++) {
            n = vectorSpacialSinglePixel(noise, i, j, window, neighbors);

            /* If n is 0, p will be 0 also */
            for (p = 0; p < n; p++)
                if (neighbors[p].j == j && neighbors[p].i == i)
                    break;

            a = (j*w)+(i*3);
            b = neighbors[0].j*w+(neighbors[0].i*3);

            for (threshold = 0; threshold < window_screen; threshold++) {
                if (p <= threshold) {
                    /* Calculate using Center Pixel */
                    scores[threshold] += ((pi->f[a  ]-noise->f[a  ])*(pi->f[a  ]-noise->f[a  ]) +
                                          (pi->f[a+1]-noise->f[a+1])*(pi->f[a+1]-noise->f[a+1]) +
                                          (pi->f[a+2]-noise->f[a+2])*(pi->f[a+2]-noise->f[a+2]));
                }
                else {
                    /* Calculate using Best Pixel */
                    scores[threshold] += ((pi->f[a  ]-noise->f[b  ])*(pi->f[a  ]-noise->f[b  ]) +
                                          (pi->f[a+1]-noise->f[b+1])*(pi->f[a+1]-noise->f[b+1]) +
                                          (pi->f[a+2]-noise->f[b+2])*(pi->f[a+2]-noise->f[b+2]));
                }
            }
        }
    }

    for (threshold = 0; threshold < window_screen; threshold++)
        scores[threshold] = sqrt( scores[threshold] / (double)(pi->height * pi->width) );

    free(neighbors);
}

void calculateAverageDifference(struct portImage *pi, int window, double noiseCorruption, double *scores) {
    struct portImage *noise = copyImage(pi);
    int window_screen = window * window;
    struct pair *neighbors;
    double score = 0.0;
    double center, best;
    int i, j, p, n, a, b;
    int w = 3*pi->width;
    int threshold;

    neighbors = (struct pair*)malloc(sizeof(struct pair)*window_screen);
    if (!neighbors) {
        fprintf(stderr, "Memory allocation failed when calculating potential RMSE.\n");
        return;
    }

    addSaltPepperNoise(noise, noiseCorruption);

    for (threshold = 0; threshold < window_screen; threshold++)
        scores[threshold] = 0.0;

    for (j = 0; j < pi->height; j++) {
        for (i = 0; i < pi->width; i++) {
            n = vectorSpacialSinglePixel(noise, i, j, window, neighbors);

            /* If n is 0, p will be 0 also */
            for (p = 0; p < n; p++)
                if (neighbors[p].j == j && neighbors[p].i == i)
                    break;

            a = (j*w)+(i*3);
            b = neighbors[0].j*w+(neighbors[0].i*3);

            center = sqrt ((pi->f[a  ]-noise->f[a  ])*(pi->f[a  ]-noise->f[a  ]) +
                           (pi->f[a+1]-noise->f[a+1])*(pi->f[a+1]-noise->f[a+1]) +
                           (pi->f[a+2]-noise->f[a+2])*(pi->f[a+2]-noise->f[a+2]));

            best   = sqrt((pi->f[a  ]-noise->f[b  ])*(pi->f[a  ]-noise->f[b  ]) +
                          (pi->f[a+1]-noise->f[b+1])*(pi->f[a+1]-noise->f[b+1]) +
                          (pi->f[a+2]-noise->f[b+2])*(pi->f[a+2]-noise->f[b+2]));

            for (threshold = 0; threshold < window_screen; threshold++) {
                if (p <= threshold) {
                    /* Calculate using Center Pixel */
                    scores[threshold] += center;
                }
                else {
                    /* Calculate using Best Pixel */
                    scores[threshold] += best;
                }
            }
        }
    }

    for (threshold = 0; threshold < window_screen; threshold++)
        scores[threshold] = scores[threshold] / (double)(pi->height * pi->width);

    free(neighbors);
    freeImage(noise);
}

void scale_reduce(struct portImage *pi, int newscale) {

    int i;
    int reduction = (pi->scale+1) / (newscale);

    for (i = 0; i < pi->size; i++)
        pi->f[i] /= reduction;

    pi->scale = newscale-1;
}

void spacial_reduce(struct portImage *pi, int resWidth, int resHeigth) {

    int reduce_h = pi->height / resWidth;
    int reduce_w = pi->width  / resHeigth;
    int w    = pi->width;
    if (pi->color == csRGB) {
        int i, j, x, y, average_red, average_green, average_blue;
        w *= 3;

        for (j = 0; j < pi->height; j += reduce_h) {
            for (i = 0; i < w; i += reduce_w*3) {

                average_red = 0;
                average_green = 0;
                average_blue = 0;
                for (y = 0; y < reduce_h; y++)
                    for (x = 0; x < reduce_w; x++) {
                        average_red   += pi->f[(j+y)*w + (i+x) + 0];
                        average_green += pi->f[(j+y)*w + (i+x) + 1];
                        average_blue  += pi->f[(j+y)*w + (i+x) + 2];
                    }

                average_red   /= reduce_h*reduce_w;
                average_green /= reduce_h*reduce_w;
                average_blue  /= reduce_h*reduce_w;

                for (y = 0; y < reduce_h; y++)
                    for (x = 0; x < reduce_w; x++) {
                        pi->f[(j+y)*w + (i+x) + 0] = average_red;
                        pi->f[(j+y)*w + (i+x) + 1] = average_green;
                        pi->f[(j+y)*w + (i+x) + 2] = average_blue;
                    }
            }
        }

    }
    else if (pi->color == csGREYSCALE) {
        int i, j, x, y, average = 0;
        for (j = 0; j < pi->height; j += reduce_h) {
            for (i = 0; i < w; i += reduce_w) {

                average = 0;
                for (y = 0; y < reduce_h; y++)
                    for (x = 0; x < reduce_w; x++)
                        average += pi->f[(j+y)*w + (i+x)];

                average /= reduce_h*reduce_w;

                for (y = 0; y < reduce_h; y++)
                    for (x = 0; x < reduce_w; x++)
                        pi->f[(j+y)*w + (i+x)] = average;
            }
        }
    }
}

void histogram(struct portImage *pi, int *hist) {
    int i;
    if (pi->color != csGREYSCALE)
        return;

    for (i = 0; i <= pi->scale; i++)
        hist[i] = 0;

    for (i = 0; i < (pi->height * pi->width); i++)
        hist[ (int)pi->f[i] ]++;
}

void cumulative_distribution(int scales, int *hist, int *cdf) {
    int i;
    int total = 0;
    for (i = 0; i < scales; i++) {
        total += hist[i];
        cdf[i] = total;
    }
}

void graph_histogram(struct portImage *pi, char *filename) {
    int i, j;
    int max;
    int width = pi->scale + 1;
    int height;
    int *hist;
    double scale;
    FILE *out;
    struct portImage histogram_out;

    histogram_out.height = width;
    histogram_out.width  = width;
    histogram_out.scale  = 1;
    histogram_out.color  = csGREYSCALE;
    histogram_out.type   = 2;
    histogram_out.c      = 0;

    histogram_out.f = (unsigned char*)malloc(sizeof(unsigned char) * width * width);
    hist = (int*)malloc(width * sizeof(int));

    if (!hist || !histogram_out.f) {
        fprintf(stderr, "ERROR: Memory for memory allocation of histogram failed.\n");
        return;
    }

    histogram(pi, hist);

    max = hist[0];
    for (i = 1; i < width; i++) /* For each bar. */
        if (hist[i] > max)
            max = hist[i];

    scale = (double)width / (double)max;

    for (i = 0; i < width; i++) { /* For each bar. */
        height = (int) ( scale * (double)hist[i] );

        for (j = width-height; j <= pi->scale; j++)
            histogram_out.f[j*width+i] = 1;
    }

    out = fopen(filename, "w");
    if (out == 0) {
        fprintf(stderr, "ERROR: Failed to open %s for writing histogram.\n", filename);
        exit(1);
    }

    writeImage(&histogram_out, out);
    fclose(out);
    free(hist);
    free(histogram_out.f);
}

void graph_cdf(struct portImage *pi, char *filename) {
    int i, j;
    int max;
    int width = pi->scale + 1;
    int height;
    int *hist;
    int *cdf;
    double scale;
    FILE *out;
    struct portImage cdf_out;

    cdf_out.height = width;
    cdf_out.width  = width;
    cdf_out.scale  = 1;
    cdf_out.color  = csGREYSCALE;
    cdf_out.type   = 2;
    cdf_out.c      = 0;

    cdf_out.f = (unsigned char*)malloc(sizeof(unsigned char) * width * width);
    hist = (int*)malloc(width * sizeof(int));
    cdf  = (int*)malloc(width * sizeof(int));

    if (!hist || !cdf || !cdf_out.f) {
        fprintf(stderr, "Memory allocation for CDF failed.\n");
        return;
    }

    histogram(pi, hist);
    cumulative_distribution(width, hist, cdf);

    max = cdf[width-1];
    scale = (double)width / (double)max;

    for (i = 0; i < width; i++) { /* For each bar. */
        height = (int) ( scale * (double)cdf[i] );

        for (j = width-height; j <= pi->scale; j++)
            cdf_out.f[j*width+i] = 1;
    }

    out = fopen(filename, "w");
    if (out == 0) {
        fprintf(stderr, "ERROR: Failed to open %s for writing histogram.\n", filename);
        exit(1);
    }

    writeImage(&cdf_out, out);
    fclose(out);
    free(cdf_out.f);
    free(hist);
    free(cdf);
}

void level_slice(struct portImage *pi, int level) {
    int i;
    for (i = 0; i < (pi->height * pi->width); i++)
        if  (level <= pi->f[i] && pi->f[i] <= level+10)
            pi->f[i] = pi->scale;
        else
            pi->f[i] = 0;
}

void equalize(struct portImage *pi) {
    int i;
    int *hist;
    int *cdf;
    int *nhist;

    hist  = (int*)malloc( (pi->scale + 1) * sizeof(int) );
    cdf   = (int*)malloc( (pi->scale + 1) * sizeof(int) );
    nhist = (int*)malloc( (pi->scale + 1) * sizeof(int) );

    if (!hist || !cdf || !nhist) {
        fprintf(stderr, "Memory allocation for image equalization failed.\n");
        return;
    }

    if (pi->color == csGREYSCALE) {
        histogram(pi, hist);
        cumulative_distribution(pi->scale+1, hist, cdf);

        for (i = 0; i <= pi->scale; i++) {
            double Crr = ( (double)cdf[i] / (double)(pi->height * pi->width) ) * (double)pi->scale;

            double Crr_dec = Crr - (double)(int)Crr;

            nhist[i] = (int)Crr;

            if (Crr_dec >= 0.5)
                nhist[i] += 1;
        }

        for (i = 0; i < (pi->height * pi->width); i++)
            pi->f[i] = nhist[ pi->f[i] ];
    }

    free(hist);
    free(cdf);
    free(nhist);
}

void contrast_stretching(struct portImage *pi, float low_bound, float high_bound) {
    int i;
    int *hist;
    int *cdf;
    int low_scale  = -1;
    int high_scale = -1;
    int low_pixels  = (int) ( (float)(pi->height * pi->width) * low_bound  );
    int high_pixels = (int) ( (float)(pi->height * pi->width) * high_bound );
    float factor;

    hist  = (int*)malloc( (pi->scale + 1) * sizeof(int) );
    cdf   = (int*)malloc( (pi->scale + 1) * sizeof(int) );

    if (!hist || !cdf) {
        fprintf(stderr, "Memory allocation for image equalization failed.\n");
        return;
    }

    histogram(pi, hist);
    cumulative_distribution(pi->scale+1, hist, cdf);

    for (i = 0; i <= pi->scale; i++) {
        if (low_scale  == -1 && cdf[i] > low_pixels)  low_scale  = i-1;
        if (high_scale == -1 && cdf[i] > high_pixels) high_scale = i;
    }

    factor = (float)(pi->scale) / (float)(high_scale - low_scale);

    for (i = 0; i < (pi->height * pi->width); i++) {
        /* Set all pixels in the bottom lowbound to black */
        if      (pi->f[i] < low_scale)
            pi->f[i] = 0;

        /* Set all of the pixels in the upper highbound to white */
        else if (pi->f[i] > high_scale)
            pi->f[i] = pi->scale;

        /* Scale whatever is left. */
        else
            pi->f[i] = (int)(factor * (float)(pi->f[i] - low_scale));
    }
    free(hist);
    free(cdf);
}
