/* from "C Language Algorithms for Digital Signal Processing"
   by Paul M. Embree and Bruce Kimble, Prentice Hall, 1991 */

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include "imgtypes.h"
#include "FFT.h"

static void FFT_process(struct Complex *x, int n, int sign)
{
   struct Complex u, temp, tm, *w, *xi, *xip, *xj, *wptr;
   int i, j, k, l, le, m, windex = 1;
   double arg, w_real, w_imag, wrecur_real, wrecur_imag, wtemp_real;

   for (m = 0; (1 << m) < n; m++);
   le = n / 2;
   w = calloc(le - 1, sizeof(struct Complex));

   arg = 4.0 * atan(1.0) / le;
   wrecur_real = w_real = cos(arg);
   wrecur_imag = w_imag = sin(arg) * sign;
   xj = w;

   for (j = 1; j < le; j++)
   {
      xj->real = wrecur_real;
      xj->imag = wrecur_imag;
      xj++;
      wtemp_real  = wrecur_real * w_real - wrecur_imag * w_imag;
      wrecur_imag = wrecur_real * w_imag + wrecur_imag * w_real;
      wrecur_real = wtemp_real;
   }

   le = n;
   for (l = 0; l < m; l++)
   {
      le /= 2;
      for (i = 0; i < n; i += 2 * le)
      {
         xi = x + i;
         xip = xi + le;
         temp.real = xi->real + xip->real;
         temp.imag = xi->imag + xip->imag;
         xip->real = xi->real - xip->real;
         xip->imag = xi->imag - xip->imag;
         *xi = temp;
      }
      wptr = w + windex - 1;
      for (j = 1; j < le; j++)
      {
         u = *wptr;
         for (i = j; i < n; i += 2 * le)
         {
            xi = x + i;
            xip = xi + le;
            temp.real = xi->real + xip->real;
            temp.imag = xi->imag + xip->imag;
            tm.real = xi->real - xip->real;
            tm.imag = xi->imag - xip->imag;
            xip->real = tm.real * u.real - tm.imag * u.imag;
            xip->imag = tm.real * u.imag + tm.imag * u.real;
            *xi = temp;
         }
         wptr += windex;
      }
      windex *= 2;
   }
   j = 0;
   for (i = 1; i < n - 1; i++)
   {
      for (k = n / 2; k <= j; k /= 2)
         j -= k;
      j += k;
      if (i < j)
      {
         xi = x + i;
         xj = x + j;
         temp = *xj;
         *xj = *xi;
         *xi = temp;
      }
   }
}

void FFT(struct Complex *x, int n)
{
   FFT_process(x, n, -1);
}

void IFFT(struct Complex *x, int n)
{
   int i;
   FFT_process(x, n, 1);
   for (i = 0; i < n; i++)
   {
      x[i].real /= n;
      x[i].imag /= n;
   }
}

