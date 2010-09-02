#ifndef __JAMES_IMAGE_IO
#define __JAMES_IMAGE_IO
#include "imgtypes.h"

void freeImage(struct portImage *pi);
void writeImage(struct portImage *pi, FILE* out);
struct portImage* readImage(FILE* fid);
struct portImage* copyImage(struct portImage *pi);

#endif
