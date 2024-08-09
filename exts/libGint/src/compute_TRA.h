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

#ifndef COMPUTE_TRA_H_COMPILE_GUARD
#define COMPUTE_TRA_H_COMPILE_GUARD

#include <vector>


void compute_TRA_batched(
   const int Nshells, const int la, const int lb, const int lc, const int ld,
   const std::vector<unsigned int>& TRA, const std::vector<double>& SPHER, std::vector<double>& OUT );

void compute_TRA_batched_low(
   const int Nshells, const int la, const int lb, const int lc, const int ld,
   const unsigned int* const __restrict__ TRA,
   const double* const __restrict__ SPHER,
   double* const __restrict__ OUT );


__global__ void compute_TRA_batched_gpu_low(
   const int Nshells, const int la, const int lb, const int lc, const int ld,
   const unsigned int* const __restrict__ TRA,
   const double* const __restrict__ SPHER,
   double* const __restrict__ OUT );

#endif
