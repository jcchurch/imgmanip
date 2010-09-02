#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include "imgio.h"

void freeImage(struct portImage *pi) {
    free(pi->f);
    if (pi->c) free(pi->c);
    if (pi->p) free(pi->p);
    free(pi);
}

void writeImage(struct portImage *pi, FILE* out) {
    int i, j;
    int w = pi->width;

    if (pi->color)
        w *= 3;

    /* Header: */
    fprintf(out, "P%d\n#\n%d %d\n%d\n", pi->type, pi->width, pi->height, pi->scale);

    /* ASCII Format */
    if (pi->type == 1 || pi->type == 2 || pi->type == 3)
        for (i = 0; i < pi->height; i++) {
            for (j = 0; j < w; j++)
                fprintf(out, "%3d ", pi->f[i*w+j]);

            fprintf(out, "\n");
        }

    /* Binary Format */
    else
        fwrite(pi->f, sizeof(unsigned char), pi->size, out);
}

struct portImage* readImage(FILE* fid) {
    char line[MAX_LINE+1];
    char p;
    int loop;
    int i;

    struct portImage *pi = (struct portImage*)malloc(sizeof(struct portImage));

    fgets(line, MAX_LINE, fid);
    sscanf(line, "%c%d", &p, &(pi->type) );

    loop = 1;
    while (loop == 1) {
        fgets(line, MAX_LINE, fid); /* Should be the comment line, but not always */

        if (line[0] != '#') {
            break;
        }
    }

    sscanf(line, "%d %d", &(pi->width), &(pi->height));
    fgets(line, MAX_LINE, fid);
    sscanf(line, "%d", &(pi->scale));

    pi->size = pi->width * pi->height;

    if (pi->type == 3 || pi->type == 6) {
        pi->color = 1;
        pi->size *= 3;
    }
    else
        pi->color = 0;

    pi->c = 0; /* This is the complex image buffer */
    pi->p = 0; /* This is the HSV image buffer */
    pi->f = (unsigned char*)malloc(pi->size * sizeof(unsigned char));

    /* ASCII Format */
    if (pi->type == 1 || pi->type == 2 || pi->type == 3)
        for (i = 0; i < pi->size; i++)
            fscanf(fid, "%d", &(pi->f[i]));

    /* Binary Format */
    else
        fread(pi->f, sizeof(char), pi->size, fid);

    return pi;
}

struct portImage* copyImage(struct portImage *pi) {
    int i;
    struct portImage *pc;
    pc = (struct portImage*)malloc(sizeof(struct portImage));
    if (!pc) {
        fprintf(stderr, "ERROR: memeory allocation failed for copying.\n");
        return pc;
    }

    pc->height = pi->height;
    pc->width  = pi->width;
    pc->size   = pi->size;
    pc->scale  = pi->scale;
    pc->color  = pi->color;
    pc->type   = pi->type;

    if (pi->f) {
        pc->f = (unsigned char*)malloc( pi->size * sizeof(unsigned char) );
        if (!pc->f) {
            fprintf(stderr, "ERROR: memory allocation failed for copying.\n");
            return pc;
        }

        for (i = 0; i < pi->size; i++)
            pc->f[i] = pi->f[i];
    }

    if (pi->c) {
        pc->c = (struct Complex*)malloc( pi->size * sizeof(struct Complex) );
        if (!pc->f) {
            fprintf(stderr, "ERROR: memory allocation failed for copying.\n");
            return pc;
        }

        for (i = 0; i < pi->size; i++)
            pc->c[i] = pi->c[i];
    }

    return pc;
}
