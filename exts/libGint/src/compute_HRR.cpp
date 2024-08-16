#include "compute_HRR.h"
#include "define.h"
#include "util.h"


//#pragma omp declare target
__device__ void execute_HRR1_gpu(
      const int AL, const int CL, const int DL,
      double* __restrict__ abcd,
      const double* const __restrict__ abpm,
      const double* const __restrict__ abcm,
      const double CD[3], const int hrr_blocksize, const int nlabcd,
      const int tid=0, const int Nt=1 ){

   const int NcoA  = NLco_dev(AL);
   const int NcoC  = NLco_dev(CL);
   const int NcoCp = NLco_dev(CL+1);
   const int NcoD  = NLco_dev(DL);
   const int NcoDm = NLco_dev(DL-1);
   
//#pragma omp parallel for 
   for( int ikl=tid ; ikl < NcoA*NcoC*NcoD ; ikl+=Nt ){
      int l = ikl%NcoD;
      int k = (ikl/NcoD)%NcoC;
      int i = ikl/NcoD/NcoC;
      int kp, lm, d ;
      int fx = lx_dev(k,CL);

      if ( l < NcoDm ){
         kp = k;
         lm = l;
         d = 0;
      } else if ( l < NcoD-1 ){
         kp = k + (CL-fx) + 1;
         lm = l - DL;
         d = 1;
      } else {
         lm = l - DL - 1;
         kp = k + (CL-fx) + 2;
         d = 2;
      }
      for ( int nl_idx = 0; nl_idx < nlabcd; nl_idx++ ){
         int idx_off = nl_idx * hrr_blocksize;
         int idx_00 = idx_off + ikl ; // (i*NcoC +k )*NcoD  + l;
         int idx_pm = idx_off + (i*NcoCp+kp)*NcoDm + lm;
         int idx_0m = idx_off + (i*NcoC +k )*NcoDm + lm;
         abcd[ idx_00 ] = abpm[idx_pm] + CD[d] * abcm[idx_0m];
      }
   }
}
//#pragma omp end declare target

void execute_HRR1_cpu(
      const int AL, const int CL, const int DL,
      double* __restrict__ abcd,
      const double* const __restrict__ abpm,
      const double* const __restrict__ abcm,
      const double CD[3], const int hrr_blocksize, const int nlabcd,
      const int tid=0, const int Nt=1 ){

   const int NcoA  = NLco(AL);
   const int NcoC  = NLco(CL);
   const int NcoCp = NLco(CL+1);
   const int NcoD  = NLco(DL);
   const int NcoDm = NLco(DL-1);
   
//#pragma omp parallel for 
   for( int ikl=tid ; ikl < NcoA*NcoC*NcoD ; ikl+=Nt ){
      int l = ikl%NcoD;
      int k = (ikl/NcoD)%NcoC;
      int i = ikl/NcoD/NcoC;
      int kp, lm, d ;
      int fx = lx(k,CL);

      if ( l < NcoDm ){
         kp = k;
         lm = l;
         d = 0;
      } else if ( l < NcoD-1 ){
         kp = k + (CL-fx) + 1;
         lm = l - DL;
         d = 1;
      } else {
         lm = l - DL - 1;
         kp = k + (CL-fx) + 2;
         d = 2;
      }
      for ( int nl_idx = 0; nl_idx < nlabcd; nl_idx++ ){
         int idx_off = nl_idx * hrr_blocksize;
         int idx_00 = idx_off + ikl ; // (i*NcoC +k )*NcoD  + l;
         int idx_pm = idx_off + (i*NcoCp+kp)*NcoDm + lm;
         int idx_0m = idx_off + (i*NcoC +k )*NcoDm + lm;
         abcd[ idx_00 ] = abpm[idx_pm] + CD[d] * abcm[idx_0m];
      }
   }
}


//#pragma omp declare target
__device__ void execute_HRR2_gpu(
      const int AL, const int BL, const int CL, const int DL,
      double* const __restrict__ abcd,
      const double* const __restrict__ pmcd,
      const double* const __restrict__ amcd,
      const double AB[3], const int hrr_blocksize, const int nlabcd,
      const int tid=0, const int Nt=1 ){

   const int NcoA = NLco_dev(AL);
   const int NcoB = NLco_dev(BL);
   const int NcoBm = NLco_dev(BL-1);
   const int NcoC = NLco_dev(CL);
   const int NcoD = NLco_dev(DL);
   const int NcoABCD = NcoA*NcoB*NcoC*NcoD;
//#pragma omp parallel for
   for( int ijkl=tid ; ijkl < NcoABCD ; ijkl+=Nt ){

      int l = ijkl%NcoD;
      int k = (ijkl/NcoD)%NcoC;
      int j = (ijkl/NcoD/NcoC)%NcoB;
      int i = ijkl/NcoD/NcoC/NcoB;// % NcoA
      int ex = lx_dev(i,AL);
      int ip, jm, d ;

      if ( j < NcoBm ){
         ip = i;
         jm = j;
         d = 0;
      } else if ( j<NcoBm+BL ) {
         ip = i + (AL-ex) + 1;
         jm = j - BL;
         d = 1;
      } else {
         ip = i + (AL-ex) + 2;
         jm = j - BL - 1;
         d = 2;
      }
      for ( int nl_idx = 0; nl_idx < nlabcd; nl_idx++ ){ // TODO Maybe fuse ?
         int idx_off = nl_idx * hrr_blocksize;
         int idx_00 = idx_off + ijkl ; // (i *NcoB +j )*NcoC*NcoD+k*NcoD+l;
         int idx_pm = idx_off + (ip*NcoBm+jm)*NcoC*NcoD+k*NcoD+l;
         int idx_0m = idx_off + (i *NcoBm+jm)*NcoC*NcoD+k*NcoD+l;
         abcd[ idx_00 ] = pmcd[ idx_pm ] + AB[d] * amcd[ idx_0m ];
      }
   }
}
//#pragma omp end declare target

void execute_HRR2_cpu(
      const int AL, const int BL, const int CL, const int DL,
      double* const __restrict__ abcd,
      const double* const __restrict__ pmcd,
      const double* const __restrict__ amcd,
      const double AB[3], const int hrr_blocksize, const int nlabcd,
      const int tid=0, const int Nt=1 ){

   const int NcoA = NLco(AL);
   const int NcoB = NLco(BL);
   const int NcoBm = NLco(BL-1);
   const int NcoC = NLco(CL);
   const int NcoD = NLco(DL);
   const int NcoABCD = NcoA*NcoB*NcoC*NcoD;
//#pragma omp parallel for
   for( int ijkl=tid ; ijkl < NcoABCD ; ijkl+=Nt ){

      int l = ijkl%NcoD;
      int k = (ijkl/NcoD)%NcoC;
      int j = (ijkl/NcoD/NcoC)%NcoB;
      int i = ijkl/NcoD/NcoC/NcoB;// % NcoA
      int ex = lx(i,AL);
      int ip, jm, d ;

      if ( j < NcoBm ){
         ip = i;
         jm = j;
         d = 0;
      } else if ( j<NcoBm+BL ) {
         ip = i + (AL-ex) + 1;
         jm = j - BL;
         d = 1;
      } else {
         ip = i + (AL-ex) + 2;
         jm = j - BL - 1;
         d = 2;
      }
      for ( int nl_idx = 0; nl_idx < nlabcd; nl_idx++ ){ // TODO Maybe fuse ?
         int idx_off = nl_idx * hrr_blocksize;
         int idx_00 = idx_off + ijkl ; // (i *NcoB +j )*NcoC*NcoD+k*NcoD+l;
         int idx_pm = idx_off + (ip*NcoBm+jm)*NcoC*NcoD+k*NcoD+l;
         int idx_0m = idx_off + (i *NcoBm+jm)*NcoC*NcoD+k*NcoD+l;
         abcd[ idx_00 ] = pmcd[ idx_pm ] + AB[d] * amcd[ idx_0m ];
      }
   }
}

__global__ void compute_HRR_batched_gpu_low(
      const int Ncells,
      const int* const __restrict__ plan,
      const unsigned int* const __restrict__ FVH,
      const double* const __restrict__ data,
      double* const __restrict__ ABCD,
      double* const __restrict__ ABCD0,
      int hrr_blocksize, int Nc, int numVC, int numVCH ){

   for( int block=blockIdx.x; block < Ncells ; block+=gridDim.x ){

      unsigned int Og     = FVH[block*FVH_SIZE+FVH_OFFSET_OG];
      unsigned int Oq     = FVH[block*FVH_SIZE+FVH_OFFSET_OQ];
      unsigned int idx_A  = FVH[block*FVH_SIZE+FVH_OFFSET_IDX_A];
      unsigned int idx_B  = FVH[block*FVH_SIZE+FVH_OFFSET_IDX_B];
      unsigned int idx_C  = FVH[block*FVH_SIZE+FVH_OFFSET_IDX_C];
      unsigned int idx_D  = FVH[block*FVH_SIZE+FVH_OFFSET_IDX_D];
      unsigned int encoded_nlabcd = FVH[block*FVH_SIZE+FVH_OFFSET_NLABCD];

      unsigned int nla,nlb,nlc,nld;
      decode4( encoded_nlabcd, &nla,&nlb,&nlc,&nld );
      int nlabcd = nla*nlb*nlc*nld;

      const double * const A =  &data[idx_A];
      const double * const B =  &data[idx_B];
      const double AB[3] = { A[0]-B[0], A[1]-B[1], A[2]-B[2] };
      const double * const C =  &data[idx_C];
      const double * const D =  &data[idx_D];
      const double CD[3] = { C[0]-D[0], C[1]-D[1], C[2]-D[2] };

      // TODO // ABSOLUTELY
      // add CELL[n1] to pbc(AB) and CELL[n2] to pbc(CD)

      double* sh_mem = &ABCD[ Og * hrr_blocksize ];

      for ( unsigned int ilabcd=0 ; ilabcd < nlabcd; ilabcd++ ){
         int s0_st = (Oq + ilabcd) * Nc;
         for ( int i=threadIdx.x ; i < Nc; i+=blockDim.x ){
            ABCD0[ s0_st + i ] = 0. ;
         }
      }

      __syncthreads();

      const int tid = threadIdx.x;
      const int Nt = blockDim.x;

      for ( int op=numVC; op < numVCH; op++ ){
         int t  = plan[ op*OP_SIZE + T__OFFSET ];
         int la = plan[ op*OP_SIZE + LA_OFFSET ];
         int lb = plan[ op*OP_SIZE + LB_OFFSET ];
         int lc = plan[ op*OP_SIZE + LC_OFFSET ];
         int ld = plan[ op*OP_SIZE + LD_OFFSET ];
         int off_m1 = plan[ op*OP_SIZE + M1_OFFSET ];
         int off_m2 = plan[ op*OP_SIZE + M2_OFFSET ];
         int off_m3 = plan[ op*OP_SIZE + M3_OFFSET ];

         if ( t == HRR1 ){
            execute_HRR1_gpu(
               la, lc, ld, &sh_mem[off_m1], &sh_mem[off_m2], &sh_mem[off_m3],
               CD, hrr_blocksize, nlabcd, tid, Nt );
         } else if ( t == HRR2 ){
            execute_HRR2_gpu(
               la, lb, lc, ld, &sh_mem[off_m1], &sh_mem[off_m2], &sh_mem[off_m3],
               AB, hrr_blocksize, nlabcd, tid, Nt );
         } else if ( t == SYBL ){
//            #pragma omp team barrier (?)
            __syncthreads();
         }
      }

      __syncthreads();
      // sums over cell from ABCD[Og:Og+nlabcd*Nc] to ABCD0[Oq:Oq+nlabcd*Nc]
      for ( unsigned int ilabcd=threadIdx.x ; ilabcd < nlabcd; ilabcd+=blockDim.x ){
         int sh_st = ilabcd * hrr_blocksize;
         int s0_st = (Oq + ilabcd) * Nc;
         for ( int i=0 ; i < Nc; i++ ){
            // must be atomic
            atomicAdd( &ABCD0[ s0_st + i ] , sh_mem[ sh_st + i ]);
         }
      } 
   }
}



void compute_HRR_batched_low(
      const int Ncells,
      const int* const __restrict__ plan,
      const unsigned int* const __restrict__ FVH,
      const double* const __restrict__ data,
      double* const __restrict__ ABCD,
      double* const __restrict__ ABCD0,
      int hrr_blocksize, int Nc, int numVC, int numVCH ){

//#pragma omp target teams num_teams(Ncells) distribute // map(to:plan,HRR,data ) map(tofrom:ABCD) map(from:ABCD0)
   for( int block=0; block < Ncells ; block++ ){


      unsigned int Og     = FVH[block*FVH_SIZE+FVH_OFFSET_OG];
      unsigned int Oq     = FVH[block*FVH_SIZE+FVH_OFFSET_OQ];
      unsigned int idx_A  = FVH[block*FVH_SIZE+FVH_OFFSET_IDX_A];
      unsigned int idx_B  = FVH[block*FVH_SIZE+FVH_OFFSET_IDX_B];
      unsigned int idx_C  = FVH[block*FVH_SIZE+FVH_OFFSET_IDX_C];
      unsigned int idx_D  = FVH[block*FVH_SIZE+FVH_OFFSET_IDX_D];
      unsigned int encoded_nlabcd = FVH[block*FVH_SIZE+FVH_OFFSET_NLABCD];

      unsigned int nla,nlb,nlc,nld;
      decode4( encoded_nlabcd, &nla,&nlb,&nlc,&nld );
      int nlabcd = nla*nlb*nlc*nld;


      const double * const A =  &data[idx_A];
      const double * const B =  &data[idx_B];
      const double AB[3] = { A[0]-B[0], A[1]-B[1], A[2]-B[2] };
      const double * const C =  &data[idx_C];
      const double * const D =  &data[idx_D];
      const double CD[3] = { C[0]-D[0], C[1]-D[1], C[2]-D[2] };

      double* sh_mem = &ABCD[ Og * hrr_blocksize ];

      for ( unsigned int ilabcd=0 ; ilabcd < nlabcd; ilabcd++ ){
         int s0_st = (Oq + ilabcd) * Nc;
         for ( int i=0 ; i < Nc; i++ ){
            ABCD0[ s0_st + i ] = 0. ;
         }
      }

      for ( int op=numVC; op < numVCH; op++ ){
         int t  = plan[ op*OP_SIZE + T__OFFSET ];
         int la = plan[ op*OP_SIZE + LA_OFFSET ];
         int lb = plan[ op*OP_SIZE + LB_OFFSET ];
         int lc = plan[ op*OP_SIZE + LC_OFFSET ];
         int ld = plan[ op*OP_SIZE + LD_OFFSET ];
         int off_m1 = plan[ op*OP_SIZE + M1_OFFSET ];
         int off_m2 = plan[ op*OP_SIZE + M2_OFFSET ];
         int off_m3 = plan[ op*OP_SIZE + M3_OFFSET ];

         if ( t == HRR1 ){
            execute_HRR1_cpu(
               la, lc, ld, &sh_mem[off_m1], &sh_mem[off_m2], &sh_mem[off_m3],
               CD, hrr_blocksize, nlabcd );
         } else if ( t == HRR2 ){
            execute_HRR2_cpu(
               la, lb, lc, ld, &sh_mem[off_m1], &sh_mem[off_m2], &sh_mem[off_m3],
               AB, hrr_blocksize, nlabcd );
         } else if ( t == SYBL ){
//            #pragma omp barrier
         }
      }

      // sums over cell from ABCD[Og:Og+nlabcd*Nc] to ABCD0[Oq:Oq+nlabcd*Nc]
      for ( unsigned int ilabcd=0 ; ilabcd < nlabcd; ilabcd++ ){
         int sh_st = ilabcd * hrr_blocksize;
         int s0_st = (Oq + ilabcd) * Nc;
         for ( int i=0 ; i < Nc; i++ ){
            // must be atomic
            ABCD0[ s0_st + i ] += sh_mem[ sh_st + i ];
         }
      }

   }
}



void compute_HRR_batched(
      const int Ncells, const std::vector<int>& Plan, const std::vector<unsigned int>& FVH, const std::vector<double>& data,
      std::vector<double>& ABCD, std::vector<double>& ABCD0, int hrr_blocksize, int Nc, int numVC, int numVCH ){
   compute_HRR_batched_low(
      Ncells, Plan.data(), FVH.data(), data.data(),
      ABCD.data(), ABCD0.data(), hrr_blocksize, Nc, numVC, numVCH );
}






