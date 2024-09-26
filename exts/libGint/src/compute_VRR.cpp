#include <iostream>
#include "compute_VRR.h"
#include "util.h"
#include "define.h"

using std::cout;
using std::endl;


void execute_VRR1( // int AL, int CL, int m, 
      const int m,
      double* __restrict__ a0c0m0,
      const double* __restrict__ amc0m0,
      const double* __restrict__ amc0mp,
      const double PA[3], const double WP[3], const int my_vrr_rank=0, const int vrr_team_size=1 ){

   for( int imm=my_vrr_rank; imm < 3*m; imm+=vrr_team_size ){

      int mm = imm / 3 ;
      int i  = imm % 3 ;

      int d  = i ;
      int im = 0 ; 

      int idx_000 = i  + 3*mm;
      int idx_m00 = im + 1*mm;

      a0c0m0[idx_000] = PA[d] * amc0m0[idx_m00] + WP[d] * amc0mp[idx_m00];
   }
}

__device__ void execute_VRR1_gpu( // int AL, int CL, int m, 
      const int m,
      double* __restrict__ a0c0m0,
      const double* __restrict__ amc0m0,
      const double* __restrict__ amc0mp,
      const double PA[3], const double WP[3], const int my_vrr_rank=0, const int vrr_team_size=1 ){

   for( int imm=my_vrr_rank; imm < 3*m; imm+=vrr_team_size ){

      int mm = imm / 3 ;
      int i  = imm % 3 ;

      int d  = i ;
      int im = 0 ; 

      int idx_000 = i  + 3*mm;
      int idx_m00 = im + 1*mm;

      a0c0m0[idx_000] = PA[d] * amc0m0[idx_m00] + WP[d] * amc0mp[idx_m00];
   }
}




// Specialization of VRR for a>1, c=0. Applies
// <a+1i,sss>(m) = PA_i <asss>(m) + WP_i <asss>(m+1) 
//               + 1/2zab N_i(a) ( <a-1,sss>(m) - rho/zab <a-1,sss>(m+1) )
void execute_VRR2(
      const int AL, const int m,
      double* __restrict__ a0c0m0,
      const double* __restrict__ amc0m0,
      const double* __restrict__ amc0mp,
      const double* __restrict__ awc0m0,
      const double* __restrict__ awc0mp,
      const double PA[3], const double WP[3], const double inv_2zab, const double min_rho_zab2,
      const int my_vrr_rank=0, const int vrr_team_size=1 ){
   const int NcoA   = (AL+1)*(AL+2)/2;
   const int NcoAxx = (AL+0)*(AL-1)/2;
   const int NcoAx  = (AL+1)*(AL+0)/2;
   const int NcoAyy = (AL+3)*(AL+0)/2 - 1;
   const int NcoAy  = (AL+1)*(AL+2)/2 - 2 ;
   const int NcoAm  = (AL+1)*(AL+0)/2;
   const int NcoAw  = (AL-1)*(AL+0)/2;


   for ( int imm = my_vrr_rank; imm < m*NcoA; imm+=vrr_team_size ){

      int i  = imm % NcoA ;
      int mm = imm / NcoA ;

      int ex = lx(i,AL);
      int ey = ly(i,AL);
      int ez = AL - ex - ey ;
      int d, im, iw, e2 ;
      if ( i < NcoAxx ){
         d  = 0;
         im = i;
         iw = i;
         e2 = ex - 1;
      } else if ( i < NcoAx ) {
         d  = 0;
         im = i;
         iw = 0;
         e2 = 0;
      } else if ( i < NcoAyy ){
         d  = 1;
         im = i - AL;
         iw = i - 2 * AL + 1;
         e2 = ey - 1;
      } else if ( i == NcoAy ){
         d  = 1;
         im = i - AL;
         iw = 0;
         e2 = 0;
      } else {
         d  = 2;
         im = i - AL - 1;
         iw = i - 2 * AL - 1;
         e2 = ez - 1;
      }
      int idx_000 = imm ; // i  + NcoA *mm ;
      int idx_m00 = im + NcoAm*mm ;
      int idx_w00 = iw + NcoAw*mm ;

      a0c0m0[ idx_000 ] = PA[d] * amc0m0[idx_m00] + 
                          WP[d] * amc0mp[idx_m00] + 
                          e2*( inv_2zab * awc0m0[idx_w00] + min_rho_zab2 * awc0mp[idx_w00] );
   }
}

__device__ void execute_VRR2_gpu(
      const int AL, const int m,
      double* __restrict__ a0c0m0,
      const double* __restrict__ amc0m0,
      const double* __restrict__ amc0mp,
      const double* __restrict__ awc0m0,
      const double* __restrict__ awc0mp,
      const double PA[3], const double WP[3], const double inv_2zab, const double min_rho_zab2,
      const int my_vrr_rank=0, const int vrr_team_size=1 ){
   const int NcoA   = (AL+1)*(AL+2)/2;
   const int NcoAxx = (AL+0)*(AL-1)/2;
   const int NcoAx  = (AL+1)*(AL+0)/2;
   const int NcoAyy = (AL+3)*(AL+0)/2 - 1;
   const int NcoAy  = (AL+1)*(AL+2)/2 - 2 ;
   const int NcoAm  = (AL+1)*(AL+0)/2;
   const int NcoAw  = (AL-1)*(AL+0)/2;


   for ( int imm = my_vrr_rank; imm < m*NcoA; imm+=vrr_team_size ){

      int i  = imm % NcoA ;
      int mm = imm / NcoA ;

      int ex = lx_dev(i,AL);
      int ey = ly_dev(i,AL);
      int ez = AL - ex - ey ;
      int d, im, iw, e2 ;
      if ( i < NcoAxx ){
         d  = 0;
         im = i;
         iw = i;
         e2 = ex - 1;
      } else if ( i < NcoAx ) {
         d  = 0;
         im = i;
         iw = 0;
         e2 = 0;
      } else if ( i < NcoAyy ){
         d  = 1;
         im = i - AL;
         iw = i - 2 * AL + 1;
         e2 = ey - 1;
      } else if ( i == NcoAy ){
         d  = 1;
         im = i - AL;
         iw = 0;
         e2 = 0;
      } else {
         d  = 2;
         im = i - AL - 1;
         iw = i - 2 * AL - 1;
         e2 = ez - 1;
      }
      int idx_000 = imm ; // i  + NcoA *mm ;
      int idx_m00 = im + NcoAm*mm ;
      int idx_w00 = iw + NcoAw*mm ;

      a0c0m0[ idx_000 ] = PA[d] * amc0m0[idx_m00] + 
                          WP[d] * amc0mp[idx_m00] + 
                          e2*( inv_2zab * awc0m0[idx_w00] + min_rho_zab2 * awc0mp[idx_w00] );
   }
}

// more general case, for al>0 and cl>1
void execute_VRR5(
      const int AL, const int CL, const int m,
      double* __restrict__ a0c0m0,
      const double* __restrict__ a0cmm0,
      const double* __restrict__ a0cmmp,
      const double* __restrict__ a0cwm0,
      const double* __restrict__ a0cwmp,
      const double* __restrict__ amcmmp,
      const double QC[3], const double WQ[3],
      const double inv_2zcd, const double min_rho_zcd2, const double inv_2z,
      const int my_vrr_rank=0, const int vrr_team_size=1 ){
   //
   const int NcoA = NLco(AL);
   const int NcoAm= NLco(AL-1);
   const int NcoC = NLco(CL);
   const int NcoCm = NLco(CL-1);
   const int NcoCw = NLco(CL-2);
   const int NcoCy = (CL+3)*CL/2-1;
   // k is the faster variable, followed by i, then m
   for ( int mik=my_vrr_rank; mik < NcoA*NcoC*m; mik+=vrr_team_size ){
      int ex, ey, ez, fx, fy, fz;
      int i, k, mm, d, km, kw, f2, im, e2 ;
      int idx_000, idx_0m0, idx_0w0, idx_mmp;
      double i_0m0, i_0mp, i_0w0, i_0wp, i_mmp;


      k  = mik % NcoC;
      i  = (mik/NcoC)%NcoA;
      mm = mik/NcoC/NcoA;

      ex = lx(i,AL);
      ey = ly(i,AL);
      ez = AL - ex - ey ; // lz(i,AL);
      fx = lx(k,CL);
      fy = ly(k,CL);
      fz = CL- fx - fy ; // lz(k,CL);
      // if the x moment of k is more than zero, than apply the vrr along x
      if (k < NcoCw){ 
         d  = 0;
         km = k;
         kw = k;
         f2 = fx-1;
         im = i;
         e2 = ex;
      // x[....]
      } else if (k < NcoCm){
         d  = 0;
         km = k;
         kw = 0; // not used since f2 = fx-1 = 0 
         f2 = fx-1; // == 0
         im = i;
         e2 = ex;
      // yy[...]
      } else if (k < NcoCy){
         d  = 1;
         km = k - CL; // k - (fy+fz) = k - (CL-fx) = k - CL
         kw = k - 2 * CL + 1;
         f2 = fy-1;
         im = i - (ey+ez); // we need the general version since ex may not be zero
         e2 = ey;
      // y[...]
      } else if (k == NcoC - 2){
          d  = 1;
          km = k - CL;
          kw = k - 2 * CL + 1;
          f2 = fy-1;
          im = i - (ey+ez);
          e2 = ey;
      // zzzzzzz
      } else {
          d  = 2;
          km = k - CL - 1;
          kw = k - 2 * CL - 1;
          f2 = fz-1;
          im = i - (ey+ez) - 1;
          e2 = ez;
      }
      idx_000 = mik ; // k +i *NcoC + NcoA *NcoC *(mm)
      idx_0m0 = km+i *NcoCm + NcoA *NcoCm*mm;
      idx_0w0 = kw+i *NcoCw + NcoA *NcoCw*mm;
      idx_mmp = km+im*NcoCm + NcoAm*NcoCm*mm;
      i_0m0 = a0cmm0[ idx_0m0 ];
      i_0mp = a0cmmp[ idx_0m0 ];
      i_0w0 = a0cwm0[ idx_0w0 ];
      i_0wp = a0cwmp[ idx_0w0 ];
      i_mmp = amcmmp[ idx_mmp ];

      a0c0m0[ idx_000 ] = QC[d] * i_0m0 + WQ[d] * i_0mp + f2*( inv_2zcd * i_0w0 + min_rho_zcd2 * i_0wp ) + e2*inv_2z*i_mmp;
   }
}

// more general case, for al>0 and cl>1
__device__ void execute_VRR5_gpu(
      const int AL, const int CL, const int m,
      double* __restrict__ a0c0m0,
      const double* __restrict__ a0cmm0,
      const double* __restrict__ a0cmmp,
      const double* __restrict__ a0cwm0,
      const double* __restrict__ a0cwmp,
      const double* __restrict__ amcmmp,
      const double QC[3], const double WQ[3],
      const double inv_2zcd, const double min_rho_zcd2, const double inv_2z,
      const int my_vrr_rank=0, const int vrr_team_size=1 ){
   //
   const int NcoA = NLco_dev(AL);
   const int NcoAm= NLco_dev(AL-1);
   const int NcoC = NLco_dev(CL);
   const int NcoCm = NLco_dev(CL-1);
   const int NcoCw = NLco_dev(CL-2);
   const int NcoCy = (CL+3)*CL/2-1;
   // k is the faster variable, followed by i, then m
   for ( int mik=my_vrr_rank; mik < NcoA*NcoC*m; mik+=vrr_team_size ){
      int ex, ey, ez, fx, fy, fz;
      int i, k, mm, d, km, kw, f2, im, e2 ;
      int idx_000, idx_0m0, idx_0w0, idx_mmp;
      double i_0m0, i_0mp, i_0w0, i_0wp, i_mmp;


      k  = mik % NcoC;
      i  = (mik/NcoC)%NcoA;
      mm = mik/NcoC/NcoA;

      ex = lx_dev(i,AL);
      ey = ly_dev(i,AL);
      ez = AL - ex - ey ; // lz(i,AL);
      fx = lx_dev(k,CL);
      fy = ly_dev(k,CL);
      fz = CL- fx - fy ; // lz(k,CL);
      // if the x moment of k is more than zero, than apply the vrr along x
      if (k < NcoCw){ 
         d  = 0;
         km = k;
         kw = k;
         f2 = fx-1;
         im = i;
         e2 = ex;
      // x[....]
      } else if (k < NcoCm){
         d  = 0;
         km = k;
         kw = 0; // not used since f2 = fx-1 = 0 
         f2 = fx-1; // == 0
         im = i;
         e2 = ex;
      // yy[...]
      } else if (k < NcoCy){
         d  = 1;
         km = k - CL; // k - (fy+fz) = k - (CL-fx) = k - CL
         kw = k - 2 * CL + 1;
         f2 = fy-1;
         im = i - (ey+ez); // we need the general version since ex may not be zero
         e2 = ey;
      // y[...]
      } else if (k == NcoC - 2){
          d  = 1;
          km = k - CL;
          kw = k - 2 * CL + 1;
          f2 = fy-1;
          im = i - (ey+ez);
          e2 = ey;
      // zzzzzzz
      } else {
          d  = 2;
          km = k - CL - 1;
          kw = k - 2 * CL - 1;
          f2 = fz-1;
          im = i - (ey+ez) - 1;
          e2 = ez;
      }
      idx_000 = mik ; // k +i *NcoC + NcoA *NcoC *(mm)
      idx_0m0 = km+i *NcoCm + NcoA *NcoCm*mm;
      idx_0w0 = kw+i *NcoCw + NcoA *NcoCw*mm;
      idx_mmp = km+im*NcoCm + NcoAm*NcoCm*mm;
      i_0m0 = a0cmm0[ idx_0m0 ];
      i_0mp = a0cmmp[ idx_0m0 ];
      i_0w0 = a0cwm0[ idx_0w0 ];
      i_0wp = a0cwmp[ idx_0w0 ];
      i_mmp = amcmmp[ idx_mmp ];

      a0c0m0[ idx_000 ] = QC[d] * i_0m0 + WQ[d] * i_0mp + f2*( inv_2zcd * i_0w0 + min_rho_zcd2 * i_0wp ) + e2*inv_2z*i_mmp;
   }
}


// case for a>0, c=1
void execute_VRR6( 
      const int AL, const int m,
      double* __restrict__ a0c0m0,
      const double* __restrict__ a0cmm0,
      const double* __restrict__ a0cmmp,
      const double* __restrict__ amcmmp,
      const double QC[3], const double WQ[3], const double inv_2z,
      const int my_vrr_rank=0, const int vrr_team_size=1 ){
   const int NcoA = NLco(AL);
   const int NcoAm= NLco(AL-1);
   const int NcoC = 3 ; // NLco(CL);
   // k is the faster variable, followed by i, then m
   for ( int mik=my_vrr_rank; mik < NcoA*NcoC*m; mik+=vrr_team_size ){
      int ex, ey, ez ;
      int i, k, mm, d, im, e2 ;
      int idx_000, idx_0m0, idx_mmp;
      double i_0m0, i_0mp, i_mmp;

      k  = mik % NcoC;
      i  = (mik/NcoC)%NcoA;
      mm = mik/NcoC/NcoA;

      ex = lx(i,AL);
      ey = ly(i,AL);
      ez = AL - ex - ey ; // lz(i,AL);

      if (k == 0 ){
         d  = 0;
         im = i;
         e2 = ex;
      } else if (k == 1){
         d  = 1;
         im = i - (ey+ez);
         e2 = ey;
      } else {
          d  = 2;
          im = i - (ey+ez) - 1;
          e2 = ez;
      }
      idx_000 = mik ; // k +i *NcoC + NcoA *NcoC *(mm)
      idx_0m0 = i  + NcoA *mm;
      idx_mmp = im + NcoAm*mm;
      i_0m0 = a0cmm0[ idx_0m0 ];
      i_0mp = a0cmmp[ idx_0m0 ];
      i_mmp = amcmmp[ idx_mmp ];

      a0c0m0[ idx_000 ] = QC[d] * i_0m0 + WQ[d] * i_0mp + e2*inv_2z*i_mmp;
   }
}


// case for a>0, c=1
__device__ void execute_VRR6_gpu( 
      const int AL, const int m,
      double* __restrict__ a0c0m0,
      const double* __restrict__ a0cmm0,
      const double* __restrict__ a0cmmp,
      const double* __restrict__ amcmmp,
      const double QC[3], const double WQ[3], const double inv_2z,
      const int my_vrr_rank=0, const int vrr_team_size=1 ){
   const int NcoA = NLco_dev(AL);
   const int NcoAm= NLco_dev(AL-1);
   const int NcoC = 3 ; // NLco_dev(CL);
   // k is the faster variable, followed by i, then m
   for ( int mik=my_vrr_rank; mik < NcoA*NcoC*m; mik+=vrr_team_size ){
      int ex, ey, ez ;
      int i, k, mm, d, im, e2 ;
      int idx_000, idx_0m0, idx_mmp;
      double i_0m0, i_0mp, i_mmp;

      k  = mik % NcoC;
      i  = (mik/NcoC)%NcoA;
      mm = mik/NcoC/NcoA;

      ex = lx_dev(i,AL);
      ey = ly_dev(i,AL);
      ez = AL - ex - ey ; // lz(i,AL);

      if (k == 0 ){
         d  = 0;
         im = i;
         e2 = ex;
      } else if (k == 1){
         d  = 1;
         im = i - (ey+ez);
         e2 = ey;
      } else {
          d  = 2;
          im = i - (ey+ez) - 1;
          e2 = ez;
      }
      idx_000 = mik ; // k +i *NcoC + NcoA *NcoC *(mm)
      idx_0m0 = i  + NcoA *mm;
      idx_mmp = im + NcoAm*mm;
      i_0m0 = a0cmm0[ idx_0m0 ];
      i_0mp = a0cmmp[ idx_0m0 ];
      i_mmp = amcmmp[ idx_mmp ];

      a0c0m0[ idx_000 ] = QC[d] * i_0m0 + WQ[d] * i_0mp + e2*inv_2z*i_mmp;
   }
}


__device__ void execute_CP2S_gpu( 
      const int AL, const int CL, 
      const double* __restrict__ pr_mem,
      double* sh_mem, // sh is -very much- not restrict-ed.
      const int my_vrr_rank, const int vrr_team_size, const int hrr_blocksize,
      const unsigned int ipa, const unsigned int ipb, const unsigned int ipc, const unsigned int ipd,
      const unsigned int nla, const int unsigned nlb, const int unsigned nlc, const int unsigned nld, 
      const int unsigned nga, const int unsigned ngb, const int unsigned ngc, const int unsigned ngd,
      const double* const __restrict__ Ka, const double* const __restrict__ Kb,
      const double* const __restrict__ Kc, const double* const __restrict__ Kd ){

   const int NcoA = NLco_dev(AL);
   const int NcoC = NLco_dev(CL);
   const int NcoAC = NcoA*NcoC;
  
   const unsigned int nl___d = nld;
   const unsigned int nl__cd = nlc*nl___d;
   const unsigned int nl_bcd = nlb*nl__cd;
   const unsigned int nlabcd = nla*nl_bcd;

   for( unsigned int ilabcd = 0; ilabcd < nlabcd; ilabcd++ ){
      unsigned int a = (ilabcd / nl_bcd) ;
      unsigned int b = (ilabcd / nl__cd) % nlb ;
      unsigned int c = (ilabcd / nl___d) % nlc ;
      unsigned int d = ilabcd % nld ;

      double K = Ka[ a*nga + ipa ] * Kb[ b*ngb + ipb ] * Kc[ c*ngc + ipc ] * Kd[ d*ngd + ipd ];

      for( int i=my_vrr_rank; i < NcoAC; i+=vrr_team_size){
        // must be atomic if different warps/blocks share the ABCD array
//        printf("CP2S %d %d adding %lg %lg @ %p : %lg \n", blockIdx.x, threadIdx.x, K , pr_mem[i], &sh_mem[ilabcd*hrr_blocksize+i], sh_mem[ilabcd*hrr_blocksize+i] );
        atomicAdd( &sh_mem[ ilabcd*hrr_blocksize + i ] , K * pr_mem[i]);
//        sh_mem[ ilabcd*hrr_blocksize + i ] += K * pr_mem[i];
      }
   }
}


__global__ void compute_VRR_batched_gpu_low(
      const int Ncells, const int* __restrict__ plan,
      const unsigned int* const __restrict__ PMX,
      const unsigned int* const __restrict__ FVH,
      const double* const __restrict__ Fm,
      const double* const __restrict__ data,
      double* const __restrict__ AC,
      double* const __restrict__ ABCD,
      int vrr_blocksize, int hrr_blocksize, int L, int numV, int numVC, const int Ng ){
   
   int F_size = L+1;
   if (L > 0){ F_size += 4*3+5; }

   for( int block=blockIdx.x; block < Ncells*Ng ; block += gridDim.x ){

      unsigned int p      = block / Ng; 
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

       // Note: VRR does not need to know the n1 n2 cell, since it is already in the PA vectors     
      unsigned int nla,nlb,nlc,nld,npa,npb,npc,npd;
      unsigned int n1,n2;
      decode_shell( nlabcd, &nla,&nlb,&nlc,&nld, &n1,&n2);
      decode4( npabcd, &npa,&npb,&npc,&npd );

      double* sh_mem = &ABCD[ Og * hrr_blocksize ];

      int best_vrr_team_size = max( 1, (L*L+1) / n_prm );
      int vrr_team_size = blockDim.x;
      while ( vrr_team_size > best_vrr_team_size ){ vrr_team_size /= 2; }

//      if ( blockIdx.x + threadIdx.x == 0  ){ printf( " VTS calc: L = %d N = %d  -> %d \n" , L, n_prm, vrr_team_size); }

      int num_vrr_teams = blockDim.x / vrr_team_size;
      int my_vrr_team = threadIdx.x / vrr_team_size;
      int my_vrr_rank = threadIdx.x % vrr_team_size;
      int n3 = block % Ng; 

      for ( unsigned i = my_vrr_team; i < n_prm ;  i += num_vrr_teams ){

         unsigned int Of   = ((Ov+i) * Ng + n3 ) * F_size;
         unsigned int ipzn = PMX[Ov+i];
         unsigned int ipa,ipb,ipc,ipd;

         // We need to know the index of the pgfs to find the K coefficents
         decode4( ipzn, &ipa,&ipb,&ipc,&ipd );

         double* pr_mem = &AC[ ((Ov+i) * Ng + n3) * vrr_blocksize ];
         for( int il=0; il < L+1; il++ ){
            pr_mem[il] = Fm[Of+il];
         }

         if (pr_mem[0] < 1.e-20 ){ continue; }

         const double* PA = nullptr;
         const double* WP = nullptr;
         const double* QC = nullptr;
         const double* WQ = nullptr;
         double inv_2zab = 0.;
         double min_rho_zab2 = 0.;
         double inv_2zcd = 0.;
         double min_rho_zcd2 = 0.;
         double inv_2z = 0.;

         if (L > 0){
            PA = &Fm[Of+L+1];
            WP = &Fm[Of+L+4];
            QC = &Fm[Of+L+7];
            WQ = &Fm[Of+L+10];
            inv_2zab = Fm[Of+L+13];
            min_rho_zab2  = Fm[Of+L+14]; // - rho/zab**2
            inv_2zcd = Fm[Of+L+15];
            min_rho_zcd2  = Fm[Of+L+16]; // - rho/zcd**2
            inv_2z = Fm[Of+L+17];
         }

         for ( int op=0; op < numVC; op++ ){
            const int t  = plan[ op*OP_SIZE + T__OFFSET ];
            const int la = plan[ op*OP_SIZE + LA_OFFSET ];
            const int lc = plan[ op*OP_SIZE + LC_OFFSET ];
            const int max_m  = plan[ op*OP_SIZE + M__OFFSET ];
            const int min_m  = plan[ op*OP_SIZE + H__OFFSET ];
            const int m = max_m - min_m + 1 ; 
            const int off_m1 = plan[ op*OP_SIZE + M1_OFFSET ];
            const int off_m2 = plan[ op*OP_SIZE + M2_OFFSET ];
            const int off_m3 = plan[ op*OP_SIZE + M3_OFFSET ];
            const int off_m4 = plan[ op*OP_SIZE + M4_OFFSET ];
            const int off_m5 = plan[ op*OP_SIZE + M5_OFFSET ];
            const int off_m6 = plan[ op*OP_SIZE + M6_OFFSET ];
            double* m1 = &pr_mem[off_m1];
            double* m2 = &pr_mem[off_m2];
            double* m3 = &pr_mem[off_m3];
            double *m4, *m5, *m6 ;

            // TODO is it worth it to unify the VRRs types as a single, maybe individually less efficient, call ?
            if ( t == VRR1 ){ 
               execute_VRR1_gpu( m, m1, m2, m3, PA, WP, my_vrr_rank, vrr_team_size );
            } else if ( t == VRR2 ){
               m4 = &pr_mem[off_m4]; m5 = &pr_mem[off_m5]; 
               execute_VRR2_gpu( la, m, m1, m2, m3, m4, m5, PA, WP, inv_2zab, min_rho_zab2, my_vrr_rank, vrr_team_size);
            } else if ( t == VRR3 ){ 
               execute_VRR1_gpu( m, m1, m2, m3, QC, WQ, my_vrr_rank, vrr_team_size );
            } else if ( t == VRR4 ){ 
               m4 = &pr_mem[off_m4]; m5 = &pr_mem[off_m5];
               execute_VRR2_gpu( lc, m, m1, m2, m3, m4, m5, QC, WQ, inv_2zcd, min_rho_zcd2, my_vrr_rank, vrr_team_size );
            } else if ( t == VRR5 ){ 
               m4 = &pr_mem[off_m4]; m5 = &pr_mem[off_m5]; m6 = &pr_mem[off_m6];
               execute_VRR5_gpu(
                  la, lc, m, m1, m2, m3, m4, m5, m6, 
                  QC, WQ, inv_2zcd, min_rho_zcd2, inv_2z, my_vrr_rank, vrr_team_size );
            } else if ( t == VRR6 ){ 
               m4 = &pr_mem[off_m4];
               execute_VRR6_gpu( la, m, m1, m2, m3, m4, QC, WQ, inv_2z, my_vrr_rank, vrr_team_size);
            } else if ( t == CP2S){
               m2 = &sh_mem[off_m2];
               execute_CP2S_gpu( 
                  la, lc, m1, m2, my_vrr_rank, vrr_team_size, hrr_blocksize,
                  ipa, ipb, ipc, ipd, nla, nlb, nlc, nld, npa, npb, npc, npd, Ka, Kb, Kc, Kd );
            } else if ( t == SYTM ){
               __syncthreads();
            }
         } // end of loop over op  
      }
   }
}




