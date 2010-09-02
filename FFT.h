#ifndef __JAMES_FFT
#define __JAMES_FFT
#include "imgtypes.h"

/*
  Fast Fourier Transform (FFT) and Inverse Fast Fourier Transform (IFFT)
  routines; both routines overwrite the input Complex array with the result;
  n is the number of elements in the input array and must be a power of 2
*/

void FFT(struct Complex *x, int n);
void IFFT(struct Complex *x, int n);

#endif
