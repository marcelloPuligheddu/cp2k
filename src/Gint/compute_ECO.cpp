#include <iostream>
#include "compute_VRR.h"
#include "util.h"
#include "define.h"

using std::cout;
using std::endl;



__global__ void compute_ECO_batched_gpu_low(
      const int Ncells, const int* __restrict__ plan,
      const unsigned int* const __restrict__ PMX,
      const unsigned int* const __restrict__ FVH,
      const double* const __restrict__ Fm, // unused
      const double* const __restrict__ data,
      double* const __restrict__ AC,
      double* const __restrict__ ABCD,
      int vrr_blocksize, int hrr_blocksize, int L, int numV, int numVC, const int Ng ){


   for( int block=blockIdx.x; block < Ncells*Ng ; block += gridDim.x ){

      unsigned int p      = block / Ng; 
      int n3 = block % Ng;

      unsigned int Ov     = FVH[p*FVH_SIZE+FVH_OFFSET_OV];
      unsigned int Og     = FVH[p*FVH_SIZE+FVH_OFFSET_OG];
      unsigned int n_prm  = FVH[p*FVH_SIZE+FVH_OFFSET_NPRM];
      unsigned int nlabcd = FVH[p*FVH_SIZE+FVH_OFFSET_NLABCD];
      unsigned int npabcd = FVH[p*FVH_SIZE+FVH_OFFSET_NPABCD];
      unsigned int idx_Ka = FVH[p*FVH_SIZE+FVH_OFFSET_IDX_KA];
      unsigned int idx_Kb = FVH[p*FVH_SIZE+FVH_OFFSET_IDX_KB];
      unsigned int idx_Kc = FVH[p*FVH_SIZE+FVH_OFFSET_IDX_KC];
      unsigned int idx_Kd = FVH[p*FVH_SIZE+FVH_OFFSET_IDX_KD];

      const double* Ka = &data[idx_Ka];
      const double* Kb = &data[idx_Kb];
      const double* Kc = &data[idx_Kc];
      const double* Kd = &data[idx_Kd];

      unsigned int nla,nlb,nlc,nld,npa,npb,npc,npd;
      decode_shell( nlabcd, &nla,&nlb,&nlc,&nld, &npa,&npb);
      decode4( npabcd, &npa,&npb,&npc,&npd );
      const unsigned int nl___d = nld;
      const unsigned int nl__cd = nlc*nl___d;
      const unsigned int nl_bcd = nlb*nl__cd;
      nlabcd = nla*nl_bcd;

      double* sh_mem = &ABCD[ Og * hrr_blocksize ];

      // arguable
      const int best_eco_team_size = 1; // (L+1) * (L+2) ; // max_NcoAC ;
      int eco_team_size = blockDim.x;
      while ( eco_team_size > best_eco_team_size ){ eco_team_size /= 2; }

      int num_eco_teams = blockDim.x / eco_team_size;
      int my_eco_team = threadIdx.x / eco_team_size;
      int my_eco_rank = threadIdx.x % eco_team_size;

      // Contraction coefficients needed by this block
      __shared__ double sKa[MAX_N_L * MAX_N_PRM];
      __shared__ double sKb[MAX_N_L * MAX_N_PRM];
      __shared__ double sKc[MAX_N_L * MAX_N_PRM];
      __shared__ double sKd[MAX_N_L * MAX_N_PRM];
      for( unsigned int idx=threadIdx.x; idx < nla * npa ; idx += blockDim.x ){ sKa[idx] = Ka[idx]; }
      for( unsigned int idx=threadIdx.x; idx < nlb * npb ; idx += blockDim.x ){ sKb[idx] = Kb[idx]; }
      for( unsigned int idx=threadIdx.x; idx < nlc * npc ; idx += blockDim.x ){ sKc[idx] = Kc[idx]; }
      for( unsigned int idx=threadIdx.x; idx < nld * npd ; idx += blockDim.x ){ sKd[idx] = Kd[idx]; }

      __syncthreads();

      // Whether a primitive is screened. Since saving this info for all primitives
      // could require too much shared memory, we tile over the them
      const int BS = 1;
      __shared__ bool passed_screening[ BS ] ; // MAX_N_PRM * MAX_N_PRM * MAX_N_PRM * MAX_N_PRM * Ng ];

      for ( unsigned int blk_prm = 0; blk_prm < n_prm ; blk_prm += BS ){

         unsigned int idx_prm_lim = min( BS, n_prm - blk_prm );

         // save screening
         for ( unsigned idx_prm = threadIdx.x; idx_prm < idx_prm_lim ; idx_prm += blockDim.x ){
            unsigned int i_prm = blk_prm+idx_prm;
            double * pr_mem = &AC[ ((Ov+i_prm) * Ng + n3) * vrr_blocksize ];
//            printf("Screening %d.%d %d.%d | %d %d -> %d @%d = %lg \n", blockIdx.x, threadIdx.x, p, n3, blk_prm, idx_prm, i_prm , ((Ov+i_prm) * Ng + n3), pr_mem[0] );
            passed_screening[idx_prm] = (pr_mem[0] > 1.e-20) ;
         }

         __syncthreads();


         for ( int op = /*numV + */ my_eco_team ; op < numVC; op += num_eco_teams ){
            
            // Find the contraction we are doing
            const int t      = plan[ op*OP_SIZE + T__OFFSET ];
            if ( t != CP2S){ continue; }
            const int off_m1 = plan[ op*OP_SIZE + M1_OFFSET ];
            const int off_m2 = plan[ op*OP_SIZE + M2_OFFSET ];
            const int la     = plan[ op*OP_SIZE + LA_OFFSET ];
            const int lc     = plan[ op*OP_SIZE + LC_OFFSET ];

            const int NcoA = NLco_dev(la);
            const int NcoC = NLco_dev(lc);
            const int NcoAC = NcoA*NcoC;

            double* m2 = &sh_mem[off_m2];

            // Each thread is assigned a unique combination of op, il and j
            // This mean each thread con compute its sum independently of the others
            for ( unsigned int i = my_eco_rank; i < nlabcd * NcoAC ; i+= eco_team_size ){

               unsigned int ilabcd    = i / NcoAC;
               unsigned int j         = i % NcoAC;

               unsigned int a = (ilabcd / nl_bcd ) % nla;
               unsigned int b = (ilabcd / nl__cd ) % nlb ;
               unsigned int c = (ilabcd / nl___d ) % nlc ;
               unsigned int d =  ilabcd            % nld ;

               // Each threads computes its sum over non screened primitives 
//               double contr_int = 0.;
               for ( unsigned idx_prm = 0; idx_prm < idx_prm_lim ;  idx_prm ++ ){

                  if ( not passed_screening[idx_prm] ){ continue; }

                  unsigned int i_prm = blk_prm+idx_prm;
                  double * pr_mem = &AC[ ((Ov+i_prm) * Ng + n3) * vrr_blocksize ];
                  double * m1 = &pr_mem[off_m1];

                  unsigned int ipzn = PMX[Ov+i_prm];
                  unsigned int ipa,ipb,ipc,ipd;
                  decode4( ipzn, &ipa,&ipb,&ipc,&ipd );

                  double K = sKa[ a*npa + ipa ] * sKb[ b*npb + ipb ] * sKc[ c*npc + ipc ] * sKd[ d*npd + ipd ];

//                  contr_int += K * m1[j];
                  atomicAdd( &m2[ ilabcd*hrr_blocksize + j ] , K * m1[j] );
               }

               // must be atomic if different warps/blocks share the ABCD [i + j] array
               // printf("CP2S %d %d adding %lg %lg @ %p : %lg \n", blockIdx.x, threadIdx.x, K , m1[j], &m2[ilabcd*hrr_blocksize+j], m2[ilabcd*hrr_blocksize+j] );
               // We need to sum over n3. Unfortunately, different n3 belonging to the same cell are assigned to different blocks
//               atomicAdd( &m2[ ilabcd*hrr_blocksize + j ] , contr_int );
               // m2[ ilabcd*hrr_blocksize + j ] += contr_int; // K * m1[j]; 
            }
         }
         __syncthreads();
      }
   }
}

