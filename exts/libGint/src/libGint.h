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

#ifndef LIBGINT_H_COMPILE_GUARD
#define LIBGINT_H_COMPILE_GUARD


#include <vector>
#include <iostream>
#include <cassert>
#include <unordered_set>
#include "plan.h"
#include "util.h"
#include "timer.h"
#include "UniqueArray.h"
#include "compute_Fm.h"
#include "compute_VRR.h"
#include "compute_HRR.h"
#include "compute_SPH.h"
#include "compute_TRA.h"
#include "compute_KS.h"
#include "fgamma.h"
#include <omp.h>

using std::max;

class libGint {
   public:
   libGint(){ 
      my_thr = omp_get_thread_num() ;
   }
   void add_prm ( const int ipa, const int ipb, const int ipc, const int ipd ) ;
   void add_shell (int i, int j , int k, int l, int n1, int n2);
   void add_cell();
   int add_qrt( int la, int lb, int lc, int ld, int nla, int nlb, int nlc, int nld );
   void add_qrtt(
         double symm_fac, int la, int lb, int lc, int ld, 
         int inla, int inlb, int inlc, int inld,
         int ld_ac, int ld_ad, int ld_bc, int ld_bd, 
         unsigned int offset_ac_L_set, unsigned int offset_ad_L_set, 
         unsigned int offset_bc_L_set, unsigned int offset_bd_L_set, 
         bool Tac, bool Tad, bool Tbc, bool Tbd );
   void add_set();
   void set_Atom( int i, double* R_, double*Z_, int np_ );
   void set_Atom_L( int i, int l_, int nl_, double* K_  );
   void set_max_n_prm( int max_n3 );
   void init();
   void set_Potential_Truncated( double R_cut, double * C0, int ld_C0, int C0_size );
   void compute_max_vector_size();
   size_t memory_needed();

   void dispatch(bool skip_cpu);
   size_t out_size = 0;

   std::vector<double> OUT;
   bool periodic = false;
   void show_state();
   void report_througput(bool skip_gpu);
//   private:

   void reset_indices();
   int my_thr = 0;
   std::vector<int> np;
   std::vector<unsigned int> idx_R;
   std::vector<unsigned int> idx_Z;
   std::vector<std::vector<int> > all_l;
   std::vector<std::vector<int> > all_nl;
   std::vector<std::vector< unsigned int>> all_idx_K;
   cublasHandle_t cublas_handle;
   cudaStream_t cuda_stream;

   int nspin = 0 ;
   double * K_a; // not owned 
   double * P_a; // not owned  
   double * K_a_dev = 0; // owned and managed 
   double * P_a_dev = 0; // owned and managed 

   double * K_b; // not owned
   double * P_b; // not owned
   double * K_b_dev = 0; // owned and managed
   double * P_b_dev = 0; // owned and managed

   double hf_fac; // K += fac * I @@ P
 
   size_t FP_size;
   void set_P( double * P, int P_size );
   void set_P( double * Pa, double * Pb, int P_size );
   void set_K( double * K, int K_size );
   void set_K( double * Ka, double * Kb, int K_size );

   void set_P( std::vector<double> & P );
   void set_K( std::vector<double> & K );
   void set_P( std::vector<double> & P_a, std::vector<double> & P_b );
   void set_K( std::vector<double> & K_a, std::vector<double> & K_b );

   void get_K( double * K );
   void get_K( double * K_a, double * K_b ); 
   void get_K( std::vector<double> & K );
   void get_K( std::vector<double> & K_a, std::vector<double> & K_b );

   void set_cell( bool periodic, double * cell_, double * cell_inv_);
   void set_neighs( double * neighs_, int nneighs );

   size_t max_integral_scratch_size = 0;
   size_t max_plan_size = 0;
   size_t max_OF_size = 0;
   size_t max_PMX_size = 0;
   size_t max_FVH_size = 0;
   size_t max_SPH_size = 0;
   size_t max_KS_size = 0;
   size_t max_TRA_size = 0;

   unsigned int offset_F[NL4] = {0};
   unsigned int offset_V[NL4] = {0};
   unsigned int offset_G[NL4] = {0};
   unsigned int offset_Q[NL4] = {0};
   unsigned int offset_T[NL4] = {0};
   int all_vrr_blocksize[NL4] = {0};
   int all_hrr_blocksize[NL4] = {0};

   std::vector<unsigned int> FVH[NL4];
   std::vector<unsigned int> OF[NL4];
   std::vector<unsigned int> PMX[NL4];
   std::vector<unsigned int> SPH[NL4];
   std::vector<unsigned int> TRA[NL4];
   std::vector<unsigned int> KS[NL4];

   unsigned int dest=0;
   std::vector<unsigned int> prm_tmp_list;
   UniqueArray ua;
//   bool is_gamma = true;
   unsigned int n_set = 0;
   unsigned int prm_in_set = 0;
   unsigned int n_prm = 0;
   int max_n_prm;
   int max_ncells;
   unsigned int p0 = 0;
   unsigned int cell_in_set = 0;
   PlanCollection plans;

   double R_cut;
   double * C0, * C0_dev;
   int ld_C0, C0_size;
   int potential_type = COULOMB; // default

   unsigned int Fm_size[NL4] = {0};
   unsigned int AC_size[NL4] = {0};
   unsigned int ABCD_size[NL4] = {0};
   unsigned int ABCD0_size[NL4] = {0};
   unsigned int SPHER_size[NL4] = {0};
   unsigned int OUT_size[NL4] = {0};

   double cell_h[9] = {0};
   double cell_inv_h[9] = {0};
   std::vector<double> neighs;

   std::unordered_set<unsigned int> encoded_moments ;
   bool first = true;
   bool all_moments[NL4] = {false} ;
   std::vector<size_t> record_of_out_sizes[NL4];
   std::vector<double> record_of_times_cpu[NL4];
   std::vector<double> record_of_times_gpu[NL4];

};





#endif // #ifndef LIBGINT_H_COMPILE_GUARD
