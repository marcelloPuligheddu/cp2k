#include "compute_KS.h"
#include "define.h"
#include "util.h"


__constant__ double symm_factors[5] = {2.0, 1.0, 0.5, 0.25, 0.0 };

__global__ void compute_KS_gpu(
      const int Nqrtt,
      const unsigned int * const __restrict__ KS,
      const int la, const int lb, const int lc, const int ld,
      const double * const __restrict__ P,
      const double * const __restrict__ I,
      double * const __restrict__ K,
      const double * const __restrict__ data,
      double hf_fac ){

   const int nsa = 2*la+1;
   const int nsb = 2*lb+1;
   const int nsc = 2*lc+1;
   const int nsd = 2*ld+1;
   const int nscd = nsc*nsd;
   const int nsbcd = nsb*nscd;
   const int nsabcd = nsa*nsbcd;

   for ( int block = blockIdx.x ; block < Nqrtt ; block += gridDim.x ){
      const unsigned int idx_fac = KS[block * KS_SIZE + KS_OFFSET_IDXFAC ];
      const unsigned int inlabcd = KS[block * KS_SIZE + KS_OFFSET_INLABCD];
      const unsigned int ld_abcd = KS[block * KS_SIZE + KS_OFFSET_LDABCD ];
      const unsigned int offset_L_set_atom_ac = KS[block * KS_SIZE + KS_OFFSET_OFFAC  ];
      const unsigned int offset_L_set_atom_ad = KS[block * KS_SIZE + KS_OFFSET_OFFAD  ];
      const unsigned int offset_L_set_atom_bc = KS[block * KS_SIZE + KS_OFFSET_OFFBC  ];
      const unsigned int offset_L_set_atom_bd = KS[block * KS_SIZE + KS_OFFSET_OFFBD  ];
      const unsigned int Tall    = KS[block * KS_SIZE + KS_OFFSET_TALL   ];

      unsigned int inla,inlb,inlc,inld;
      unsigned int ldac,ldad,ldbc,ldbd;
      unsigned int  Tac, Tad, Tbc, Tbd;

      decode4( inlabcd, &inla, &inlb, &inlc, &inld );
      decode4( ld_abcd, &ldac, &ldad, &ldbc, &ldbd );
      decode4( Tall   , &Tac,  &Tad,  &Tbc,  &Tbd  );

      const int Oq = block * nsabcd ;
      // MAYBE move the hf_fac multiplication outside this kernel, since K initial value is zero, only multiply once at the end
      // Cost would go from N(integral computeted) / n_gpu to SIZEOF( K ) for every gpu, so it may not be always useful
      // OR move the hf_fac multiplication to the GCC factors
      const double fac = symm_factors[idx_fac] * hf_fac;
      
      
      for ( int t = threadIdx.x ; t < nsabcd; t += blockDim.x ){


//         if ( block + t == 0 ){ printf( " ! P ! %p \n", P ); }
//         if ( block + t == 0 ){ printf( " ! P ! %p %lg         \n", P , P[0] ); }
//         if ( block + t == 0 ){ printf( " ! P ! %p %lg %lg     \n", P , P[0], P[1] ); }
//         if ( block + t == 0 ){ printf( " ! P ! %p %lg %lg %lg \n", P , P[0], P[1], P[2] ); }

         const int isa =   t / nsbcd;
         const int isb = ( t / nscd ) % nsb;
         const int isc = ( t / nsd  ) % nsc;
         const int isd =   t          % nsd;


//         const int idx_a = inla*nsa + isa;
//         const int idx_b = inlb*nsb + isb;
//         const int idx_c = inlc*nsc + isc;
//         const int idx_d = inld*nsd + isd;

         const int idx_a = isa;
         const int idx_b = isb;
         const int idx_c = isc;
         const int idx_d = isd;


//         const int idx_ac = Tac ?  idx_c * ldac + idx_a :  idx_a * ldac + idx_c;
//         const int idx_ad = Tad ?  idx_d * ldad + idx_a :  idx_a * ldad + idx_d;
//         const int idx_bc = Tbc ?  idx_c * ldbc + idx_b :  idx_b * ldbc + idx_c;
//         const int idx_bd = Tbd ?  idx_d * ldbd + idx_b :  idx_b * ldbd + idx_d;

         const int idx_ac_T = not Tac ?  idx_c * ldac + idx_a :  idx_a * ldac + idx_c;
         const int idx_ad_T = not Tad ?  idx_d * ldad + idx_a :  idx_a * ldad + idx_d;
         const int idx_bc_T = not Tbc ?  idx_c * ldbc + idx_b :  idx_b * ldbc + idx_c;
         const int idx_bd_T = not Tbd ?  idx_d * ldbd + idx_b :  idx_b * ldbd + idx_d;

//         const int Iac = offset_L_set_atom_ac + idx_ac;
//         const int Iad = offset_L_set_atom_ad + idx_ad;
//         const int Ibc = offset_L_set_atom_bc + idx_bc;
//         const int Ibd = offset_L_set_atom_bd + idx_bd;

         const int Iac_T = offset_L_set_atom_ac + idx_ac_T;
         const int Iad_T = offset_L_set_atom_ad + idx_ad_T;
         const int Ibc_T = offset_L_set_atom_bc + idx_bc_T;
         const int Ibd_T = offset_L_set_atom_bd + idx_bd_T;

         const double iabcd = - fac * I[ Oq + t ];
         const double kbd = iabcd * P [Iac_T];
         const double kbc = iabcd * P [Iad_T];
         const double kad = iabcd * P [Ibc_T];
         const double kac = iabcd * P [Ibd_T];

//         printf( " ---------- KS GPU %lg %lg %lg %lg %lg \n", iabcd, kbd,kbc,kad,kac );
//         printf("KS GPU %d.%d.0: Adding %4.10lg ( - %lg * %lg * %lg ) to %lg from P %d [%p] @ K %d [%p] \n ", block, t, kbd, fac, I[ Oq + t ], P[Iac_T], K[Ibd_T], Iac_T, &P[Iac_T], Ibd_T, &K[Ibd_T] );
//         printf("KS GPU %d.0: Adding %4.10lg ( - %lg * %lg * %lg ) to %lg from P %d @ K %d \n ", block, kbd, fac, I[ Oq + t ], P[Iac_T], K[Ibd_T], Iac_T, Ibd_T );
//         printf("KS GPU %d.1: Adding %4.10lg ( - %lg * %lg * %lg ) to %lg from P %d @ K %d \n ", block, kbc, fac, I[ Oq + t ], P[Iad_T], K[Ibc_T], Iad_T, Ibc_T );
//         printf("KS GPU %d.2: Adding %4.10lg ( - %lg * %lg * %lg ) to %lg from P %d @ K %d \n ", block, kad, fac, I[ Oq + t ], P[Ibc_T], K[Iad_T], Ibc_T, Iad_T );
//         printf("KS GPU %d.3: Adding %4.10lg ( - %lg * %lg * %lg ) to %lg from P %d @ K %d \n ", block, kac, fac, I[ Oq + t ], P[Ibd_T], K[Iac_T], Ibd_T, Iac_T );

         // MUST be atomics on device, or however K is distributed
         atomicAdd( &K[ Ibd_T ] , kbd);
         atomicAdd( &K[ Ibc_T ] , kbc);
         atomicAdd( &K[ Iad_T ] , kad);
         atomicAdd( &K[ Iac_T ] , kac);

      }
   }
}

/////////////// OUTDATED ////////////////////////

void compute_KS(
      const int Nqrtt,
      const unsigned int * const __restrict__ KS,
      const int la, const int lb, const int lc, const int ld,
      const double * const __restrict__ P,
      const double * const __restrict__ I,
      double * const __restrict__ K,
      const double * const __restrict__ data,
      double hf_fac ){

   const int nsa = 2*la+1;
   const int nsb = 2*lb+1;
   const int nsc = 2*lc+1;
   const int nsd = 2*ld+1;
   const int nscd = nsc*nsd;
   const int nsbcd = nsb*nscd;
   const int nsabcd = nsa*nsbcd;

   for ( int block = 0 ; block < Nqrtt ; block++ ){
      const unsigned int idx_fac = KS[block * KS_SIZE + KS_OFFSET_IDXFAC ];
      const unsigned int inlabcd = KS[block * KS_SIZE + KS_OFFSET_INLABCD];
      const unsigned int ld_abcd = KS[block * KS_SIZE + KS_OFFSET_LDABCD ];
      const unsigned int offset_L_set_atom_ac = KS[block * KS_SIZE + KS_OFFSET_OFFAC  ];
      const unsigned int offset_L_set_atom_ad = KS[block * KS_SIZE + KS_OFFSET_OFFAD  ];
      const unsigned int offset_L_set_atom_bc = KS[block * KS_SIZE + KS_OFFSET_OFFBC  ];
      const unsigned int offset_L_set_atom_bd = KS[block * KS_SIZE + KS_OFFSET_OFFBD  ];
      const unsigned int Tall    = KS[block * KS_SIZE + KS_OFFSET_TALL   ];

      unsigned int inla,inlb,inlc,inld;
      unsigned int ldac,ldad,ldbc,ldbd;
      unsigned int  Tac, Tad, Tbc, Tbd;

      decode4( inlabcd, &inla, &inlb, &inlc, &inld );
      decode4( ld_abcd, &ldac, &ldad, &ldbc, &ldbd );
      decode4( Tall   , &Tac,  &Tad,  &Tbc,  &Tbd  );

      const unsigned int Oq = block * nsabcd ;
      const double fac = data[idx_fac];
      
      for ( int t = 0 ; t < nsabcd; t++ ){
         const int isa =   t / nsbcd;
         const int isb = ( t / nscd ) % nsb;
         const int isc = ( t / nsd  ) % nsc;
         const int isd =   t          % nsd;

         const int idx_a = inla*nsa + isa;
         const int idx_b = inlb*nsb + isb;
         const int idx_c = inlc*nsc + isc;
         const int idx_d = inld*nsd + isd;

         const int idx_ac = Tac ?  idx_c * ldac + idx_a :  idx_a * ldac + idx_c;
         const int idx_ad = Tad ?  idx_d * ldad + idx_a :  idx_a * ldad + idx_d;
         const int idx_bc = Tbc ?  idx_c * ldbc + idx_b :  idx_b * ldbc + idx_c;
         const int idx_bd = Tbd ?  idx_d * ldbd + idx_b :  idx_b * ldbd + idx_d;

         const int Iac = offset_L_set_atom_ac + idx_ac;
         const int Iad = offset_L_set_atom_ad + idx_ad;
         const int Ibc = offset_L_set_atom_bc + idx_bc;
         const int Ibd = offset_L_set_atom_bd + idx_bd;
         const double iabcd = - fac * I[ Oq + t ];
         const double pac = P [Iac];
         const double pad = P [Iad];
         const double pbc = P [Ibc];
         const double pbd = P [Ibd];
         // must be atomics
         // MAYBE move the hf_fac multiplication outside this kernel, since K initial value is zero, only multiply once at the end
         // Cost would go from N(integral computeted) / n_gpu to SIZEOF( K ) for every gpu, so it may not be always useful
         K[ Ibd ] += hf_fac * pac * iabcd;
         K[ Ibc ] += hf_fac * pad * iabcd;
         K[ Iad ] += hf_fac * pbc * iabcd;
         K[ Iac ] += hf_fac * pbd * iabcd;
      }
   }
}



