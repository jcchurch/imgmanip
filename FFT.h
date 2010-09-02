/*
  Fast Fourier Transform (FFT) and Inverse Fast Fourier Transform (IFFT)
  routines; both routines overwrite the input Complex array with the result;
  n is the number of elements in the input array and must be a power of 2
*/

typedef struct { double real, imag; } Complex;

void FFT(Complex *x, int n);
void IFFT(Complex *x, int n);
