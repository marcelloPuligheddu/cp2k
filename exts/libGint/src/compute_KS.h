/*
Copyright (c) 2023 Science and Technology Facilities Council

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.
*/

#ifndef COMPUTE_KS_H_COMPILE_GUARD
#define COMPUTE_KS_H_COMPILE_GUARD

void compute_KS(
      const int Nqrtt,
      const unsigned int * const __restrict__ KS,
      const int la, const int lb, const int lc, const int ld,
      const double * const __restrict__ P,
      const double * const __restrict__ I,
      double * const __restrict__ K,
      const double * const __restrict__ data,
      double hf_fac  );
 
__global__ void compute_KS_gpu(
      const int Nqrtt,
      const unsigned int * const __restrict__ KS,
      const int la, const int lb, const int lc, const int ld,
      const double * const __restrict__ P,
      const double * const __restrict__ I,
      double * const __restrict__ K,
      const double * const __restrict__ data,
      double hf_fac );


#endif
