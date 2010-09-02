#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include "imgio.h"
#include "imgfilter.h"

int compare (const void * a, const void * b) {
    double d = (*(struct pair*)a).v - (*(struct pair*)b).v;
    if      (d > 0.0)
        return 1;
    else if (d < 0.0)
        return -1;

    return 0;
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
