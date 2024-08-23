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
#include "fgamma.h"
#include "libGint.h"
#include "c2s.h"
#include <omp.h>

#include <iomanip>      // std::setprecision

using std::max;





void libGint::show_state(){

   for (unsigned int L : encoded_moments ){
      int la,lb,lc,ld;
      decodeL(L,&la,&lb,&lc,&ld);
      cout << " Moments : " << la << " " << lb << " " << lc << " " << ld << endl;
   }

}

inline unsigned int encode_ipabcd_n123( const int ipa, const int ipb, const int ipc, const int ipd, const int n1, const int n2, const int n3 ){
    unsigned int ret = 0;
    ret +=  n3;
    ret +=  n2 * 16;
    ret +=  n1 * 16 * 16;
    ret += ipd * 16 * 16 * 16;
    ret += ipc * 16 * 16 * 16 * 32;
    ret += ipb * 16 * 16 * 16 * 32 * 32;
    ret += ipa * 16 * 16 * 16 * 32 * 32 * 32;    
    return ret;
}

int max( std::vector<int> x ){
   if ( x.size() == 0 ){ return 0; };
   int ret = x[0];
   for( int idx=1; idx<x.size(); idx++ ){
      ret = max(ret, x[idx]);
   }
   return ret;
}

// TODO max_n3 is set by nneighs, which is also max_ncells
void libGint::set_max_n_prm( int max_n3 ){
   max_n_prm  = max_n3;
   max_n_prm *= max( np );
   max_n_prm *= max( np );
   max_n_prm *= max( np );
   max_n_prm *= max( np );
   prm_tmp_list.resize( max_n_prm );
   n_prm = 0;
   cout << " resized prm to accept up to " << max_n_prm << " prm from " << max_n3 << " cells and " << max(np) << " prm " << endl;
}

void libGint::init(){
   PUSH_RANGE("libGint init",1);

   Timer timer;
   timer.start();

   for( int la=0; la <= 1; la++ ){
   for( int lb=0; lb <= 1; lb++ ){
   for( int lc=0; lc <= 1; lc++ ){
   for( int ld=0; ld <= 1; ld++ ){
      int L = encodeL( la,lb,lc,ld );
      PMX[L].reserve( 10000 );
       OF[L].reserve( 10000 );
   }}}}

   CUDA_GPU_ERR_CHECK( cudaStreamCreate( &cuda_stream ));
   CUBLAS_GPU_ERR_CHECK( cublasCreate(&cublas_handle) );
   CUBLAS_GPU_ERR_CHECK( cublasSetStream( cublas_handle, cuda_stream ));

   potential_type = COULOMB; // default

   timer.stop();
   { cout << "Cuda create stream " << cuda_stream << " @ " << &cuda_stream << " \n" ; cout.flush(); }
   POP_RANGE;
}

void libGint::set_Potential_Truncated( double R_cut_, double * C0_, int ld_C0_, int C0_size_ ){
   C0_size = C0_size_;
   ld_C0 = ld_C0_;
   R_cut = R_cut_;
   C0 = C0_;
   potential_type = TRUNCATED;

   CUDA_GPU_ERR_CHECK( cudaMalloc( (void**)&C0_dev, C0_size * sizeof(double) ) );
   CUDA_GPU_ERR_CHECK( cudaMemcpy( C0_dev, C0, C0_size * sizeof(double), cudaMemcpyHostToDevice ));   

//   cout << " Setting C0 " << C0_size << " | " << ld_C0 << endl;
//   for ( int ic=0; ic < C0_size; ic++ ){
//      cout << " " << C0[ic] << " " ;
//      if ( ic % ld_C0 == ld_C0-1 ){ cout << endl ; }
//   } // cout << endl;
}

void libGint::set_Atom( int i, double* R_, double* Z_, int np_ ){
   if ( i >= idx_R.size() ){ idx_R.resize(i+1); }
   if ( i >= idx_Z.size() ){ idx_Z.resize(i+1); }
   if ( i >= np.size() ){ np.resize(i+1); }

   unsigned int tmp_idx_R = ua.add( R_, 3 );
   unsigned int tmp_idx_Z = ua.add( Z_, np_ );
   idx_R[i] = tmp_idx_R;
   idx_Z[i] = tmp_idx_Z;
   np[i] = np_;
   if ( i < all_l.size() ) {
      all_l[i].clear();
      all_nl[i].clear();
      all_idx_K[i].clear();
   }
}

void libGint::set_Atom_L( int i, int l_, int nl_, double* K_ ){
   if ( i >= all_l.size()     ){ all_l.resize( i+1 ); }
   if ( i >= all_nl.size()    ){ all_nl.resize( i+1 ); }
   if ( i >= all_idx_K.size() ){ all_idx_K.resize( i+1 ); }
   all_l[i].push_back(l_);
   all_nl[i].push_back(nl_);
   all_idx_K[i].push_back( ua.add( K_, np[i]*nl_ ) );
}

void libGint::add_prm( const int ipa, const int ipb, const int ipc, const int ipd, const int n1, const int n2, const int n3 ){
//   cout << "|" << ipa << ipb << ipc << ipd << '|' << n1 << n2 << n3 ;
//   cout.flush();
   unsigned int piabcdxyz = encode_ipabcd_n123(ipa,ipb,ipc,ipd,n1,n2,n3);
   prm_tmp_list[ n_prm ] = piabcdxyz;
   n_prm++;
}

void libGint::add_shell ( int i, int j, int k, int l, int n1, int n2 ){
   
//   Timer timer;
//   Timer timer2;

   if (n_prm == 0){
      return;
   }

   int nnla = all_l[i].size();
   int nnlb = all_l[j].size();
   int nnlc = all_l[k].size();
   int nnld = all_l[l].size();

//   cout << endl;

   for( int idx_la=0; idx_la < nnla; idx_la++ ){
      unsigned int idx_Ka = all_idx_K[i][idx_la];
      int nla = all_nl[i][idx_la];
      int la = all_l[i][idx_la];
      for( int idx_lb=0; idx_lb < nnlb; idx_lb++ ){
         unsigned int idx_Kb = all_idx_K[j][idx_lb];
         int nlb = all_nl[j][idx_lb];
         int lb = all_l[j][idx_lb];
         for( int idx_lc=0; idx_lc < nnlc; idx_lc++ ){
            unsigned int idx_Kc = all_idx_K[k][idx_lc];
            int nlc = all_nl[k][idx_lc];
            int lc = all_l[k][idx_lc];
            for( int idx_ld=0; idx_ld < nnld; idx_ld++ ){
               unsigned int idx_Kd = all_idx_K[l][idx_ld];
               int nld = all_nl[l][idx_ld];
               int ld = all_l[l][idx_ld];

               const unsigned int N_cc = nla*nlb*nlc*nld;
               const unsigned int L = encodeL(la,lb,lc,ld);

               const unsigned int Of = offset_F[L];
               const unsigned int Ov = offset_V[L];
               const unsigned int Og = offset_G[L];
               const unsigned int Oq = offset_Q[L];

               const unsigned int encoded_nlabcd = encode_ipabcd_n123(nla,nlb,nlc,nld,n1,n2,0);
               const unsigned int encoded_npabcd = encode4(np[i],np[j],np[k],np[l]);

//               cout << " TH " << my_thr << " Inserting " << Of << " repeated " << n_prm << " times into OF[" << la<<lb<<lc<<ld << "] at " << OF[L].size() << " | " << PMX[L].size() << endl;

               OF[L].insert(  OF[L].end(), n_prm, Of );

//               for ( int prm_idx=0; prm_idx<n_prm; prm_idx++ ){ cout << prm_tmp_list[prm_idx] << " "; } cout << endl;

               PMX[L].insert( PMX[L].end(), &prm_tmp_list[0], &prm_tmp_list[n_prm] );

               const unsigned int tmp[FVH_SIZE] = {
                  Ov, Og, Oq, n_prm, idx_R[i], idx_R[j], idx_R[k], idx_R[l],
                  idx_Z[i], idx_Z[j], idx_Z[k], idx_Z[l], idx_Ka, idx_Kb, idx_Kc, idx_Kd,
                  encoded_nlabcd, encoded_npabcd
               };

               cout << " Adding to FVH | " << la << lb << lc << ld << " | " << i << " " << j << " " << k << " " << l << " " << n1 << " " << n2 << " | " ;
               for ( int tmp_idx=0; tmp_idx<FVH_SIZE; tmp_idx++ ){
                  cout << tmp[tmp_idx] << " " ;
                  if ( (tmp_idx==2) or (tmp_idx==3) or (tmp_idx==7) or (tmp_idx==11) or (tmp_idx==15) ){ cout << "|" ; }
               } cout << endl;

               FVH[L].insert( FVH[L].end(), tmp, tmp+FVH_SIZE );

               const int labcd = la+lb+lc+ld;
               Fm_size[L] += (1+labcd) * n_prm;
               if ( labcd > 0 ){
                  Fm_size[L] += (4*3+5) * n_prm;
               }

               if ( all_moments[L] == false ){

//                  cout << " planning " << la << " " << lb << " " << lc << " " << ld << endl;

                  std::vector<int> * plan = NULL ;
                  unsigned int vrr_blocksize, hrr_blocksize, numV, numVC, numVCH;
                  plans.get( la, lb, lc, ld, &plan, &vrr_blocksize, &hrr_blocksize, &numV, &numVC, &numVCH );
                  all_vrr_blocksize[L] = vrr_blocksize;
                  all_hrr_blocksize[L] = hrr_blocksize;
                  all_moments[L] = true;
               }

               cout << " AC size " << la<<lb<<lc<<ld << " += " << all_vrr_blocksize[L] << " " << n_prm << endl;
               AC_size[L] += all_vrr_blocksize[L] * n_prm;
               ABCD_size[L] += all_hrr_blocksize[L] * N_cc;

               offset_G[L] += N_cc;
               offset_V[L] += n_prm;
               offset_F[L] ++ ;

               encoded_moments.insert(L);
}}}}

}

void libGint::add_cell() {

//   cout << " TH " << my_thr << " Adding cell with " << n_prm << endl;
   if (n_prm == 0){
      return;
   }

   cell_in_set++;
   n_prm = 0;

}

int libGint::add_qrt( int la, int lb, int lc, int ld, int nla, int nlb, int nlc, int nld ){
   unsigned int nlabcd = nla*nlb*nlc*nld;
   unsigned int L = encodeL(la,lb,lc,ld);

   unsigned int SPH_idxs[SPH_SIZE] = {0} ;
   SPH_idxs[SPH_OFFSET_Q     ] = offset_Q[L];
   SPH_idxs[SPH_OFFSET_NLABCD] = nlabcd;
   SPH[L].insert(SPH[L].end(), SPH_idxs, SPH_idxs+SPH_SIZE);

   unsigned int nla_as_uint = (unsigned int) nla ;
   unsigned int nlb_as_uint = (unsigned int) nlb ;
   unsigned int nlc_as_uint = (unsigned int) nlc ;
   unsigned int nld_as_uint = (unsigned int) nld ;
   unsigned int TRA_idxs[TRA_SIZE] = {0};
   TRA_idxs[TRA_OFFSET_Q   ] = offset_Q[L];
   TRA_idxs[TRA_OFFSET_NLA ] = nla_as_uint;
   TRA_idxs[TRA_OFFSET_NLB ] = nlb_as_uint;
   TRA_idxs[TRA_OFFSET_NLC ] = nlc_as_uint;
   TRA_idxs[TRA_OFFSET_NLD ] = nld_as_uint;
   TRA_idxs[TRA_OFFSET_DEST] = dest;
   TRA[L].insert(TRA[L].end(), TRA_idxs, TRA_idxs+TRA_SIZE );

   ABCD0_size[L] += compute_Nc(la,lb,lc,ld) * nlabcd ;
   SPHER_size[L] += compute_Ns(la,lb,lc,ld) * nlabcd ;
   OUT_size[L] += compute_Ns(la,lb,lc,ld) * nlabcd ; 
   offset_Q[L] += nlabcd ;
   offset_T[L] += 1 ;

   unsigned int prev_dest = dest;
   dest += compute_Ns(la,lb,lc,ld) * nlabcd ;
   return prev_dest;


}

void libGint::add_qrtt(
      double symm_fac, int la, int lb, int lc, int ld, 
      int inla, int inlb, int inlc, int inld,
      int ld_ac, int ld_ad, int ld_bc, int ld_bd, 
      unsigned int offset_ac_L_set, unsigned int offset_ad_L_set, 
      unsigned int offset_bc_L_set, unsigned int offset_bd_L_set, 
      bool Tac, bool Tad, bool Tbc, bool Tbd ){

  

   unsigned int idx_fac ;
   // The symm_fac can be only 2,1,0.5 or 0.25 (or 0 in dev version). We assign them to 0,1,2 and 3 (and 4) respectively
   // We encode it here, decode it inside compute_KS

   if ( symm_fac < 0.1  ){ idx_fac = 4 ; }
   else if ( symm_fac < 0.3  ){ idx_fac = 3 ; }
   else if ( symm_fac < 0.75 ){ idx_fac = 2 ; }
   else if ( symm_fac < 1.5  ){ idx_fac = 1 ; }
   else { idx_fac = 0 ; }

   unsigned int L = encodeL(la,lb,lc,ld);
   unsigned int KS_idxs[KS_SIZE] = {0};
   KS_idxs[KS_OFFSET_IDXFAC ] = idx_fac;
   KS_idxs[KS_OFFSET_INLABCD] = encode4(  inla, inlb, inlc, inld );
   KS_idxs[KS_OFFSET_LDABCD ] = encode4( ld_ac,ld_ad,ld_bc,ld_bd );
   KS_idxs[KS_OFFSET_OFFAC  ] = offset_ac_L_set;
   KS_idxs[KS_OFFSET_OFFAD  ] = offset_ad_L_set;
   KS_idxs[KS_OFFSET_OFFBC  ] = offset_bc_L_set;
   KS_idxs[KS_OFFSET_OFFBD  ] = offset_bd_L_set;
   KS_idxs[KS_OFFSET_TALL   ] = encode4(  (int)Tac, (int)Tad, (int)Tbc, (int)Tbd );

//   cout << endl;
//   cout << " KS adding qrtt at weight " << symm_fac << " i: |" << inla << " " << inlb << " " << inlc << " " << inld << "|" ;
//   cout << "     l: |" << la << " " << lb << " " << lc << " " << ld << "|" << endl;
//   cout << "     offsets are " << offset_ac_L_set << " " << offset_ad_L_set << " " << offset_bc_L_set << " " << offset_bd_L_set << endl;
//   cout << "     Ts are " << Tac << Tad << Tbc << Tbd << endl;
//   cout << endl;

   KS[L].insert( KS[L].end(), KS_idxs, KS_idxs+KS_SIZE );

//   cout << " KS done " << endl ; cout.flush();

}

void libGint::add_set(){
   prm_in_set = 0;
   n_set += 1;
}

double Lcorrection(int L){
   if (L == 0){
       return 1.0;
   } else if ( L == 1 ) {
       return sqrt(3.);
   } else {
        return 2. * sqrt(M_PI);
   }
}

__global__ void apply_correction( int N, double* mem_dev, double corr ){
   for( int i= blockIdx.x*blockDim.x + threadIdx.x ; i < N ; i += blockDim.x*gridDim.x ){
      mem_dev[i] *= corr;
   }
}


void libGint::compute_max_vector_size(){

   max_integral_scratch_size = 0;
   max_plan_size = 0;
   max_OF_size = 0;
   max_PMX_size = 0;
   max_FVH_size = 0;
   max_SPH_size = 0;
   max_KS_size  = 0;
   max_TRA_size = 0;
   out_size = 0;

   for ( unsigned int L : encoded_moments ){

//      int la,lb,lc,ld;
//      decodeL(L,&la,&lb,&lc,&ld);
//      std::vector<int> * plan = NULL ;
//      unsigned int vrr_blocksize, hrr_blocksize, numV, numVC, numVCH;
//      plans.get( la, lb, lc, ld, &plan, &vrr_blocksize, &hrr_blocksize, &numV, &numVC, &numVCH );
//      max_plan_size  = max(max_plan_size,  plan->size());
      max_plan_size = 1e6;
      size_t integral_scratch_size = Fm_size[L] + AC_size[L] + ABCD_size[L] + ABCD0_size[L] + SPHER_size[L] ;
      max_integral_scratch_size = max( max_integral_scratch_size, integral_scratch_size );

      max_OF_size    = max(max_OF_size,    OF[L].size());
      max_PMX_size   = max(max_PMX_size,  PMX[L].size());
      max_FVH_size   = max(max_FVH_size,  FVH[L].size());
      max_SPH_size   = max(max_SPH_size,  SPH[L].size());
      max_KS_size    = max(max_KS_size,    KS[L].size());
      max_TRA_size   = max(max_TRA_size,  TRA[L].size());

      out_size += OUT_size[L];

   }

}

size_t libGint::memory_needed( ){
   compute_max_vector_size();
   size_t tmp = (max_plan_size + max_OF_size + max_PMX_size + max_FVH_size + max_SPH_size + max_KS_size + max_TRA_size);
   size_t add1 = tmp * sizeof(unsigned int);
   size_t add2 = (max_integral_scratch_size + out_size) * sizeof(double);
   size_t add3 = 2 * FP_size * sizeof(double);
   return add1+add2+add3;
}

void libGint::set_cell( bool periodic_, double * cell_h_,  double * cell_i_ ){
   periodic = periodic_; 
   for (int i=0;i<9;i++){ cell_h[CELL_HMAT_OFF+i] = cell_h_[i]; }
   for (int i=0;i<9;i++){ cell_h[CELL_HINV_OFF+i] = cell_i_[i]; }
}

void libGint::set_neighs( double * neighs_, int nneighs_ ){
   max_ncells = nneighs_;
   neighs.resize(3*nneighs_);
   for (int i=0;i<3*nneighs_;i++){ neighs[i] = neighs_[i]; }
}

void libGint::set_P( double * P_, int P_size ){
   nspin = 1;
   P_a = P_ ;
   FP_size = P_size;
   CUDA_GPU_ERR_CHECK( cudaMalloc( (void**)&P_a_dev, sizeof(double)*FP_size ));
   CUDA_GPU_ERR_CHECK( cudaMemcpy( P_a_dev, P_a, sizeof(double)*FP_size, cudaMemcpyHostToDevice ));  
}
void libGint::set_P( std::vector<double> & P_ ){ set_P(P_.data(), P_.size()); }

void libGint::set_P( double * P_a_, double * P_b_, int P_size ){
//   cout << " Setting P Polarized with size " << P_a_.size() << " " << P_b_.size() << endl ;
   nspin = 2;
   P_a = P_a_;
   P_b = P_b_;
   FP_size = P_size;
   CUDA_GPU_ERR_CHECK( cudaMalloc( (void**)&P_a_dev, sizeof(double)*FP_size ));
   CUDA_GPU_ERR_CHECK( cudaMemcpy( P_a_dev, P_a_, sizeof(double)*FP_size, cudaMemcpyHostToDevice ));
   CUDA_GPU_ERR_CHECK( cudaMalloc( (void**)&P_b_dev, sizeof(double)*FP_size ));
   CUDA_GPU_ERR_CHECK( cudaMemcpy( P_b_dev, P_b_, sizeof(double)*FP_size, cudaMemcpyHostToDevice ));
}
void libGint::set_P( std::vector<double> & P_a_, std::vector<double> & P_b_ ){ set_P( P_a_.data(), P_b_.data(), P_a_.size()); }

void libGint::set_K( double * K_ , int K_size ){
   assert( nspin == 1 );
   K_a = K_;
   FP_size = K_size;
   CUDA_GPU_ERR_CHECK( cudaMalloc( (void**)&K_a_dev, sizeof(double)*FP_size ));
   CUDA_GPU_ERR_CHECK( cudaMemcpy( K_a_dev, K_, sizeof(double)*FP_size, cudaMemcpyHostToDevice ));

//   cout << " Setting K with size " << K_size << " @ " << K_a << " and " << K_a_dev << endl ;

}
void libGint::set_K( std::vector<double> & K_ ){ set_K( K_.data(), K_.size()); }

void libGint::set_K( double * K_a_, double * K_b_, int K_size ){
//   cout << " Setting K Polarized with size " << K_a_.size() << " " << K_b_.size() << endl ;
   assert( nspin == 2 );
   K_a = K_a_;
   K_b = K_b_;
   FP_size = K_size;
   CUDA_GPU_ERR_CHECK( cudaMalloc( (void**)&K_a_dev, sizeof(double)*FP_size ));
   CUDA_GPU_ERR_CHECK( cudaMemcpy( K_a_dev, K_a_, sizeof(double)*FP_size, cudaMemcpyHostToDevice ));
   CUDA_GPU_ERR_CHECK( cudaMalloc( (void**)&K_b_dev, sizeof(double)*FP_size ));
   CUDA_GPU_ERR_CHECK( cudaMemcpy( K_b_dev, K_b_, sizeof(double)*FP_size, cudaMemcpyHostToDevice ));
}
void libGint::set_K( std::vector<double> & K_a_ , std::vector<double> & K_b_ ){ set_K( K_a_.data(), K_b_.data(), K_a_.size()); }

void libGint::get_K( double * K_ ){
   assert( nspin == 1 );
   CUDA_GPU_ERR_CHECK( cudaMemcpy( K_, K_a_dev, sizeof(double)*FP_size, cudaMemcpyDeviceToHost ));
}
void libGint::get_K( std::vector<double> & K_a_ ){ get_K( K_a_.data()); }

void libGint::get_K( double * K_a_,  double * K_b_ ){
   assert( nspin == 2 );
   CUDA_GPU_ERR_CHECK( cudaMemcpy( K_a_, K_a_dev, sizeof(double)*FP_size, cudaMemcpyDeviceToHost ));
   CUDA_GPU_ERR_CHECK( cudaMemcpy( K_b_, K_b_dev, sizeof(double)*FP_size, cudaMemcpyDeviceToHost ));
}
void libGint::get_K( std::vector<double> & K_a_,  std::vector<double> & K_b_ ){ get_K( K_a_.data(), K_b_.data()); }

void libGint::reset_indices(){
   max_integral_scratch_size = 0;
   max_plan_size = 0;
   max_OF_size = 0;
   max_PMX_size = 0;
   max_FVH_size = 0;
   max_SPH_size = 0;
   max_TRA_size = 0;
   out_size = 0;

   for ( unsigned int L : encoded_moments ){
       OF[L].clear();
      PMX[L].clear();
      FVH[L].clear();
      SPH[L].clear();
      TRA[L].clear();
       KS[L].clear();
      offset_F[L] = 0;
      offset_V[L] = 0;
      offset_G[L] = 0;
      offset_Q[L] = 0;
      offset_T[L] = 0;
      Fm_size[L] = 0;
      AC_size[L] = 0;
      ABCD_size[L] = 0;
      ABCD0_size[L] = 0;
      SPHER_size[L] = 0;
      OUT_size[L] = 0;
   }

   dest=0;
   n_set = 0;
   prm_in_set = 0;
   n_prm = 0;
   p0 = 0;
   cell_in_set = 0;

   encoded_moments.clear();

}

void libGint::dispatch( bool skip_cpu ){

   { cout << "Dispatch on stream " << cuda_stream << " @ " << &cuda_stream << endl; cout.flush(); }
   CUDA_GPU_ERR_CHECK( cudaStreamSynchronize(cuda_stream) );
   CUDA_GPU_ERR_CHECK( cudaPeekAtLastError() );
 
   compute_max_vector_size();

   #pragma omp critical
   { cout << "Dispatch: Will compute " << out_size << " values " << endl; cout.flush(); }
   OUT.resize(out_size);


   int ftable_ld = 0; // ld of table for fgamma
   int nmax = 21;
   double tdelta = 0.1;  // hardcoded in source code
   double tmin = 0.; // hardcoded in source code
   double tmax = 12.;  // hardcoded in source code
   int n = nmax + 6; // 6 is the max taylor expansion. Also hardcoded
   int itabmin = floor(tmin/tdelta); // === 0
   int itabmax = int((tmax - tmin)/tdelta); // === 120
   int nelem = (itabmax - itabmin + 1 ) * (n+1); // === 121*(n+1) == 121*ftable_ld
   double* ftable = create_md_ftable( nmax, tmin, tmax, tdelta, &ftable_ld);

   size_t tot_mem = 0;
   tot_mem += OUT.size()*sizeof(double);
   tot_mem += ua.internal_buffer.size()*sizeof(double) ;
   tot_mem += sizeof(double)*max_integral_scratch_size;
   tot_mem += sizeof(int)*max_plan_size          ;
   tot_mem += sizeof(unsigned int)*max_OF_size  ; 
   tot_mem += sizeof(unsigned int)*max_PMX_size  ; 
   tot_mem += sizeof(unsigned int)*max_FVH_size  ;
   tot_mem += sizeof(unsigned int)*max_SPH_size  ;
   tot_mem += sizeof(unsigned int)*max_KS_size  ;
   tot_mem += sizeof(unsigned int)*max_TRA_size  ;
   tot_mem += sizeof(unsigned int)*(3*max_ncells+2*9+nelem+245) ;
   tot_mem += 2 * sizeof(double)*(FP_size) ;

   #pragma omp single
   if ( first ){
     cout << "Memory use: (B)" << endl;
     cout << " OUT DAT SCRT PLAN OF PMX FVH SPH KS TRA AUX FK TOT" << endl;
     cout.flush();
   }


   first = false;
   #pragma omp barrier

   #pragma omp critical
   {
   double scale = 1.0;
   cout << int( OUT.size()*sizeof(double) *scale ) << " ";
   cout << int( ua.internal_buffer.size()*sizeof(double) *scale ) << " ";
   cout << int( sizeof(double)*max_integral_scratch_size *scale ) << " ";
   cout << int( sizeof(int)*max_plan_size          *scale ) << "  ";
   cout << int( sizeof(unsigned int)*max_OF_size   *scale ) << "  "; 
   cout << int( sizeof(unsigned int)*max_PMX_size  *scale ) << "  "; 
   cout << int( sizeof(unsigned int)*max_FVH_size  *scale ) << "  ";
   cout << int( sizeof(unsigned int)*max_SPH_size  *scale ) << "  ";
   cout << int( sizeof(unsigned int)*max_KS_size   *scale ) << "  ";
   cout << int( sizeof(unsigned int)*max_TRA_size  *scale ) << "  ";
   cout << int( sizeof(unsigned int)*(9+nelem+245) *scale ) << "  "; 
   cout << int( 2 * sizeof(double)*(FP_size)       *scale ) << "  "; 
   cout << int( tot_mem *scale ) << endl; cout.flush();
   }

   double *data_dev, *cell_h_dev, *neighs_dev, *ftable_dev, *OUT_dev, *C2S_dev;
   double *integral_scratch_dev;
   unsigned int *OF_dev, *PMX_dev, *FVH_dev, *SPH_dev, *KS_dev, *TRA_dev;
   int *plan_dev;

   PUSH_RANGE("dispatch malloc",1);
   CUDA_GPU_ERR_CHECK( cudaMalloc( (void**)&OUT_dev, sizeof(double)*OUT.size() ));
   CUDA_GPU_ERR_CHECK( cudaMalloc( (void**)&data_dev, sizeof(double)*(ua.internal_buffer.size()) ));
   CUDA_GPU_ERR_CHECK( cudaMalloc( (void**)&cell_h_dev, sizeof(double)*(2*9) ));
   CUDA_GPU_ERR_CHECK( cudaMalloc( (void**)&neighs_dev, sizeof(double)*(3*max_ncells) ));
   CUDA_GPU_ERR_CHECK( cudaMalloc( (void**)&ftable_dev, sizeof(double)*(nelem) ));
   CUDA_GPU_ERR_CHECK( cudaMalloc( (void**)&C2S_dev, sizeof(double)*245 ));
   CUDA_GPU_ERR_CHECK( cudaMalloc( (void**)&integral_scratch_dev, sizeof(double)*max_integral_scratch_size ));
   CUDA_GPU_ERR_CHECK( cudaMalloc( (void**)&plan_dev,sizeof(int)*max_plan_size ));
   CUDA_GPU_ERR_CHECK( cudaMalloc( (void**)&OF_dev , sizeof(unsigned int)*max_OF_size )); 
   CUDA_GPU_ERR_CHECK( cudaMalloc( (void**)&PMX_dev, sizeof(unsigned int)*max_PMX_size )); 
   CUDA_GPU_ERR_CHECK( cudaMalloc( (void**)&FVH_dev, sizeof(unsigned int)*max_FVH_size ));
   CUDA_GPU_ERR_CHECK( cudaMalloc( (void**)&SPH_dev, sizeof(unsigned int)*max_SPH_size )); 
   CUDA_GPU_ERR_CHECK( cudaMalloc( (void**)&KS_dev , sizeof(unsigned int)*max_KS_size  )); 
   CUDA_GPU_ERR_CHECK( cudaMalloc( (void**)&TRA_dev, sizeof(unsigned int)*max_TRA_size )); 
   CUDA_GPU_ERR_CHECK( cudaDeviceSynchronize() );
   CUDA_GPU_ERR_CHECK( cudaPeekAtLastError() );
   POP_RANGE;

   PUSH_RANGE("dispatch memcy",1);
   CUDA_GPU_ERR_CHECK( cudaMemcpyAsync(
      data_dev, ua.internal_buffer.data(), sizeof(double)*(ua.internal_buffer.size()), cudaMemcpyHostToDevice, cuda_stream ));
   CUDA_GPU_ERR_CHECK( cudaMemcpyAsync(
      cell_h_dev, cell_h, sizeof(double)*(2*9), cudaMemcpyHostToDevice, cuda_stream ));
   CUDA_GPU_ERR_CHECK( cudaMemcpyAsync(
      neighs_dev, neighs.data(), sizeof(double)*(3*max_ncells), cudaMemcpyHostToDevice, cuda_stream )); 
   CUDA_GPU_ERR_CHECK( cudaMemcpyAsync(
      ftable_dev, ftable, sizeof(double)*(nelem), cudaMemcpyHostToDevice, cuda_stream ));
   CUDA_GPU_ERR_CHECK( cudaMemcpyAsync(
      C2S_dev, c2s, sizeof(double)*245, cudaMemcpyHostToDevice, cuda_stream ));
   POP_RANGE;
   CUDA_GPU_ERR_CHECK( cudaStreamSynchronize(cuda_stream) );
   CUDA_GPU_ERR_CHECK( cudaPeekAtLastError() );

   // Main cycle. 
   // 1) Get the plan
   // 2) Copy the input vectors to device memory
   // 3) Run

   PUSH_RANGE("dispatch all L",3);
   for ( unsigned int L : encoded_moments ){

      int la,lb,lc,ld,labcd;
      decodeL(L,&la,&lb,&lc,&ld);
      labcd = la+lb+lc+ld;
      int Nc = compute_Nc(la,lb,lc,ld);
      int Ns = compute_Ns(la,lb,lc,ld);
      double corr = Lcorrection(la)*Lcorrection(lb)*Lcorrection(lc)*Lcorrection(ld);

      std::vector<int> * plan = NULL ;
      unsigned int vrr_blocksize, hrr_blocksize, numV, numVC, numVCH;   

      cout << " L " << la << "" << lb << "" << lc << "" << ld << " | "; cout.flush();

      unsigned int Nprm   = offset_V[L];
      unsigned int Ncells = offset_F[L];
      unsigned int Nqrtt  = offset_Q[L];

      cout << "Computing " << Nprm << " prms " << Ncells << " cells " << Nqrtt << " qrtts " << endl; cout.flush();

      double* Fm_dev    = &integral_scratch_dev[0];
      double* AC_dev    = Fm_dev    + Fm_size[L];
      double* ABCD_dev  = AC_dev    + AC_size[L];
      double* ABCD0_dev = ABCD_dev  + ABCD_size[L];
      double* SPHER_dev = ABCD0_dev + ABCD0_size[L];

      std::string Lname = std::to_string(la) + "_" + std::to_string(lb) + "_" + std::to_string(lc) + "_" + std::to_string(ld);
      PUSH_RANGE(Lname.c_str(),3);

      plans.get( la, lb, lc, ld, &plan, &vrr_blocksize, &hrr_blocksize, &numV, &numVC, &numVCH );

      PUSH_RANGE("transfer indeces",4);
      CUDA_GPU_ERR_CHECK( cudaMemcpyAsync(
         plan_dev, plan->data(), sizeof(int)*(plan->size()), cudaMemcpyHostToDevice, cuda_stream ));
      CUDA_GPU_ERR_CHECK( cudaMemcpyAsync(
          OF_dev,  OF[L].data(), sizeof(unsigned int)*(Nprm), cudaMemcpyHostToDevice, cuda_stream ));  
      CUDA_GPU_ERR_CHECK( cudaMemcpyAsync(
         PMX_dev, PMX[L].data(), sizeof(unsigned int)*(Nprm), cudaMemcpyHostToDevice, cuda_stream )); 
      CUDA_GPU_ERR_CHECK( cudaMemcpyAsync(
         FVH_dev, FVH[L].data(), sizeof(unsigned int)*(FVH[L].size()), cudaMemcpyHostToDevice, cuda_stream ));
      CUDA_GPU_ERR_CHECK( cudaMemcpyAsync(
          KS_dev,  KS[L].data(), sizeof(unsigned int)*( KS[L].size()), cudaMemcpyHostToDevice, cuda_stream )); 

//      CUDA_GPU_ERR_CHECK( cudaMemcpyAsync(
//         TRA_dev, TRA[L].data(), sizeof(unsigned int)*(TRA[L].size()), cudaMemcpyHostToDevice, cuda_stream ));

      POP_RANGE; // transfer indeces

      CUDA_GPU_ERR_CHECK( cudaStreamSynchronize(cuda_stream) );
      CUDA_GPU_ERR_CHECK( cudaPeekAtLastError() );
      CUDA_GPU_ERR_CHECK( cudaDeviceSynchronize() );

      int Fm_blocksize = 256;
      int Fm_numblocks = (Nprm+Fm_blocksize-1)/Fm_blocksize;
      PUSH_RANGE("compute",5);
      compute_Fm_batched_low_gpu<<<Fm_numblocks,Fm_blocksize,0,cuda_stream>>>(
         FVH_dev, OF_dev, PMX_dev, data_dev, Fm_dev, Nprm, labcd,
         periodic, cell_h_dev, neighs_dev, ftable_dev, ftable_ld,R_cut,C0_dev,ld_C0,potential_type );
      CUDA_GPU_ERR_CHECK( cudaPeekAtLastError() );
      CUDA_GPU_ERR_CHECK( cudaDeviceSynchronize() );

      std::vector<double> FM_on_cpu(Fm_size[L]);
      CUDA_GPU_ERR_CHECK( cudaMemcpy( FM_on_cpu.data(),  Fm_dev, sizeof(double)*(Fm_size[L]), cudaMemcpyDeviceToHost) );
      cout << " FM " << endl;
      for( int ifm=0; ifm < Fm_size[L]; ifm++ ){
         cout << ifm << " " << std::setprecision(16) << FM_on_cpu[ifm] << endl;
      } cout << endl;

      std::vector<unsigned int> FVH_on_cpu(FVH[L].size());
      CUDA_GPU_ERR_CHECK( cudaMemcpy( FVH_on_cpu.data(),  FVH_dev, sizeof(unsigned int )*(FVH[L].size()), cudaMemcpyDeviceToHost) );
      cout << " FVH <- " << endl;
      for( int ifm=0; ifm < FVH[L].size(); ifm++ ){
         cout << " " << FVH_on_cpu[ifm] ;
         if (ifm % FVH_SIZE == FVH_SIZE-1 ){ cout << endl ; }
      } cout << endl;

      CUDA_GPU_ERR_CHECK( cudaPeekAtLastError() );
      CUDA_GPU_ERR_CHECK( cudaDeviceSynchronize() );


      compute_VRR_batched_gpu_low<<<Ncells,1,0,cuda_stream>>>(
         Ncells, plan_dev, PMX_dev, FVH_dev, Fm_dev, data_dev,
         AC_dev, ABCD_dev, vrr_blocksize, hrr_blocksize, labcd, numV, numVC ); 
      CUDA_GPU_ERR_CHECK( cudaPeekAtLastError() );
      CUDA_GPU_ERR_CHECK( cudaDeviceSynchronize() );

      std::vector<double> AC_on_cpu(AC_size[L]);
      CUDA_GPU_ERR_CHECK( cudaMemcpy( AC_on_cpu.data(),  AC_dev, sizeof(double)*(AC_size[L]), cudaMemcpyDeviceToHost) );
      cout << " AC " << AC_size[L] << endl;
      for( int ifm=0; ifm < AC_size[L]; ifm++ ){
         cout << ifm << " " << std::setprecision(16) << AC_on_cpu[ifm] << endl;
      } cout << endl;

      std::vector<unsigned int> FVH2_on_cpu(FVH[L].size());
      CUDA_GPU_ERR_CHECK( cudaMemcpy( FVH2_on_cpu.data(),  FVH_dev, sizeof(unsigned int )*(FVH[L].size()), cudaMemcpyDeviceToHost) );
      cout << " FVH 2 <- " << endl;
      for( int ifm=0; ifm < FVH[L].size(); ifm++ ){
         cout << " " << FVH2_on_cpu[ifm] ;
         if (ifm % FVH_SIZE == FVH_SIZE-1 ){ cout << endl ; }
      } cout << endl;


      for ( int ii = 0 ; ii < FVH[L].size() ; ii++ ){
         cout << FVH[L][ii] << " " ;
         if ( ii % FVH_SIZE == FVH_SIZE-1 ){ cout << endl ; }
      }
    
      compute_HRR_batched_gpu_low<<<Ncells,128,0,cuda_stream>>>(
         Ncells, plan_dev, FVH_dev, data_dev, ABCD_dev, ABCD0_dev,
         periodic, cell_h_dev, neighs_dev,
         hrr_blocksize, Nc, numVC, numVCH );

      std::vector<double> ABCD0_on_cpu(ABCD0_size[L]);
      CUDA_GPU_ERR_CHECK( cudaMemcpy( ABCD0_on_cpu.data(),  ABCD0_dev, sizeof(double)*(ABCD0_size[L]), cudaMemcpyDeviceToHost) );
      cout << " ABCD0 " << endl;
      for( int ifm=0; ifm < ABCD0_size[L]; ifm++ ){
         cout << ifm << " " << std::setprecision(16) << ABCD0_on_cpu[ifm] << endl;
      } cout << endl;


      // Note: we need to DeviceSynchronize before going from kernels to cublas. TODO actually check it is true
      CUDA_GPU_ERR_CHECK( cudaPeekAtLastError() );
      CUDA_GPU_ERR_CHECK( cudaDeviceSynchronize() );

      compute_SPH_batched_gpu_alt ( Nqrtt, la, lb, lc, ld, ABCD0_dev, SPHER_dev, ABCD_dev, C2S_dev, cublas_handle );
      CUDA_GPU_ERR_CHECK( cudaPeekAtLastError() );
      CUDA_GPU_ERR_CHECK( cudaDeviceSynchronize() );

//      int corrBS = 64;
//      int corrNB = (Nqrtt*Ns+corrBS-1)/corrBS;
//      apply_correction<<<corrNB,corrBS>>>( Nqrtt*Ns, SPHER_dev, corr );
//      CUDA_GPU_ERR_CHECK( cudaPeekAtLastError() );
//      CUDA_GPU_ERR_CHECK( cudaDeviceSynchronize() );

      compute_KS_gpu<<<Nqrtt,128,0,cuda_stream>>>( Nqrtt, KS_dev, la,lb,lc,ld, P_a_dev, SPHER_dev, K_a_dev, data_dev, hf_fac );
      if ( nspin == 2 ){
         compute_KS_gpu<<<Nqrtt,128,0,cuda_stream>>>( Nqrtt, KS_dev, la,lb,lc,ld, P_b_dev, SPHER_dev, K_b_dev, data_dev, hf_fac );
      }
      CUDA_GPU_ERR_CHECK( cudaPeekAtLastError() );
      CUDA_GPU_ERR_CHECK( cudaDeviceSynchronize() );

//      compute_TRA_batched_gpu_low<<<Nshell,128>>>( Nshell, la, lb, lc, ld, TRA_dev, SPHER_dev, OUT_dev );

      POP_RANGE; // compute
      POP_RANGE; // Lname
   }
   POP_RANGE; // compute all L

//   std::vector<double> OUT_from_gpu( OUT.size() );
//   CUDA_GPU_ERR_CHECK( cudaMemcpy( OUT_from_gpu.data(), OUT_dev, sizeof(double)*(OUT.size()), cudaMemcpyDeviceToHost ));



   std::vector<double> F_a_from_gpu( FP_size );

//   cout << endl;
//   cout << " Copying K back from gpu to cpu. " << endl;
//   cout << " Gpu address is " << K_a_dev << endl;
//   cout << " Cpu address is " << F_a_from_gpu.data() << endl;
//   cout << " Copying " << FP_size << " elements " << endl;
//   cout << endl;

   CUDA_GPU_ERR_CHECK( cudaMemcpy( F_a_from_gpu.data(), K_a_dev, sizeof(double)*(FP_size), cudaMemcpyDeviceToHost ));
   if ( nspin == 2 ){
      std::vector<double> F_b_from_gpu( FP_size );
      CUDA_GPU_ERR_CHECK( cudaMemcpy( F_b_from_gpu.data(), K_b_dev, sizeof(double)*(FP_size), cudaMemcpyDeviceToHost ));
   }

//   for ( int ipf=0; ipf < FP_size; ipf++ ){ cout <<  ipf << " " << F_a_from_gpu[ipf] << endl ; } cout << endl;
  
   CUDA_GPU_ERR_CHECK( cudaFree(OUT_dev) );
   CUDA_GPU_ERR_CHECK( cudaFree(data_dev) );
   CUDA_GPU_ERR_CHECK( cudaFree(cell_h_dev) );
   CUDA_GPU_ERR_CHECK( cudaFree(neighs_dev) );
   CUDA_GPU_ERR_CHECK( cudaFree(ftable_dev) );
   CUDA_GPU_ERR_CHECK( cudaFree(C2S_dev) );
   CUDA_GPU_ERR_CHECK( cudaFree(integral_scratch_dev) );
   CUDA_GPU_ERR_CHECK( cudaFree(plan_dev) );
   CUDA_GPU_ERR_CHECK( cudaFree( OF_dev) );
   CUDA_GPU_ERR_CHECK( cudaFree(PMX_dev) );
   CUDA_GPU_ERR_CHECK( cudaFree(FVH_dev) );
   CUDA_GPU_ERR_CHECK( cudaFree(SPH_dev) );
   CUDA_GPU_ERR_CHECK( cudaFree( KS_dev) );
   CUDA_GPU_ERR_CHECK( cudaFree(TRA_dev) );

   reset_indices();

}


template < class T > 
T sum( std::vector< T > x ){
   T ret = (T)0;
   for( auto i : x ){
      ret += i;
   }
   return ret;
}

void libGint::report_througput(bool skip_cpu){
   cout << " la lb lc ld L OUT(MB) ";
   if ( not skip_cpu ) { cout << "CPU_Throughput(GB/s) " ; }
   cout << " GPU_Throughput(GB/s) " << endl;
/*
   for ( auto L : all_moments ){

      int la,lb,lc,ld,labcd;
      decodeL(L,&la,&lb,&lc,&ld);
      labcd = la+lb+lc+ld;
 
      size_t sum_output_size = sum( record_of_out_sizes[L] );
      double sum_times_cpu   = sum( record_of_times_cpu[L] );
      double sum_times_gpu   = sum( record_of_times_gpu[L] );

      double avg_thr_cpu = sum_output_size / sum_times_cpu * sizeof(double) / 1.e3;
      double avg_thr_gpu = sum_output_size / sum_times_gpu * sizeof(double) / 1.e3;

      cout << la << " " << lb << " " << lc << " " << ld << " " << labcd << " " ;
      cout << sum_output_size / 1.e6 << " " ;
      if ( not skip_cpu ) { cout << avg_thr_cpu << " " ; }
      cout << avg_thr_gpu ;
      cout << endl;

   }
*/
}

