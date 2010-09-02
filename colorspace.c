#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include "imgio.h"

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

    pi->color = csRGB;
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

    pi->color = csHSV;
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

void HSV_ValueFilter(struct portImage *pi, int window) {
    int            i,j,k,l,n, point, better_point;
    int            value = (window * window) / 2;
    int            w = pi->width * 3;
    unsigned char *buffer;
    int            window_screen = window*window;
    int            bound = window >> 1;
    struct pair   *neighbors;

    if (pi->color == csGREYSCALE)
        return;

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
