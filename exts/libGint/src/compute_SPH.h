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

#ifndef COMPUTE_SPH_H_COMPILE_GUARD
#define COMPUTE_SPH_H_COMPILE_GUARD

#include <vector>
#include "cublas_v2.h"

void compute_SPH_batched(
      const int Nqrtt, const int la, const int lb, const int lc, const int ld,
      std::vector<double>& ABCD0, std::vector<double>& SPHER, std::vector<double>& tmp_scratch );

void compute_SPH_batched_low(
      const int Nqrtt, const int la, const int lb, const int lc, const int ld,
      double* const __restrict__ ABCD0,
      double* const __restrict__ SPHER,
      double* const __restrict__ tmp_scratch );


__global__ void compute_SPH_batched_gpu_low(
      const int Nqrtt, const int la, const int lb, const int lc, const int ld,
      double* const __restrict__ ABCD0,
      double* const __restrict__ SPHER,
      double* const __restrict__ tmp_scratch );

void compute_SPH_batched_gpu_alt(
      const int Nqrtt, const int la, const int lb, const int lc, const int ld,
      double* const __restrict__ ABCD0,
      double* const __restrict__ SPHER,
      double* const __restrict__ tmp_scratch, 
      double* const __restrict__ C2S_dev, cublasHandle_t handle );



#endif // #ifndef COMPUTE_SPH_H_COMPILE_GUARD 
