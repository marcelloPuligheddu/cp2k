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

#ifndef COMPUTE_HRR_H_COMPILE_GUARD
#define COMPUTE_HRR_H_COMPILE_GUARD

#include <vector>

void compute_HRR_batched(
      const int Ncells, const std::vector<int>& Plan, const std::vector<unsigned int>& FVH, const std::vector<double>& data,
      std::vector<double>& ABCD, std::vector<double>& ABCD0, int hrr_blocksize, int Nc, int numVC, int numVCH );

__global__ void compute_HRR_batched_gpu_low(
      const int Ncells,
      const int* const __restrict__ plan,
      const unsigned int* const __restrict__ FVH,
      const double* const __restrict__ data,
      double* const __restrict__ ABCD,
      double* const __restrict__ ABCD0,
      bool periodic,
      const double* const __restrict__ cell,
      const double* const __restrict__ neighs,
      int hrr_blocksize, int Nc, int numVC, int numVCH );

void compute_HRR_batched_low(
      const int Ncells,
      const int* const __restrict__ plan,
      const unsigned int* const __restrict__ FVH,
      const double* const __restrict__ data,
      double* const __restrict__ ABCD,
      double* const __restrict__ ABCD0,
      int hrr_blocksize, int Nc, int numVC, int numVCH );



#endif // COMPUTE_HRR_H_COMPILE_GUARD
