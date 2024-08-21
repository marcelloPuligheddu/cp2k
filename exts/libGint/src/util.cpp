#include <cassert>
#include "util.h"
#include "define.h"



__device__ __host__ void decode4(
   unsigned int abcd, unsigned int* a, unsigned int* b,
   unsigned int* c, unsigned int * d ){
(*d) = abcd % 256;
(*c) = abcd / 256 % 256 ;
(*b) = abcd / (256*256) % 256 ;
(*a) = abcd / (256*256*256) ;
}


unsigned int encodeL( int la, int lb, int lc, int ld ){
   return la * NL3 + lb * NL2 + lc * NL + ld;
}


__device__ __host__ void decodeL( unsigned int L, int* la, int* lb, int* lc, int* ld ){
   (*ld) = L % NL;
   (*lc) = L / NL % NL ;
   (*lb) = L / NL2 % NL ;
   (*la) = L / NL3 ;
}

/*
unsigned int encode_ipabcd_n123( const int ipa, const int ipb, const int ipc, const int ipd, const int n1, const int n2, const int n3 ){
    unsigned int ret = 0;
    ret += n3+8;
    ret += (n2+8) * 16;
    ret += (n1+8) * 16 * 16;
    ret += (ipd) * 16 * 16 * 16;
    ret += (ipc) * 16 * 16 * 16 * 32;
    ret += (ipb) * 16 * 16 * 16 * 32 * 32;
    ret += (ipa) * 16 * 16 * 16 * 32 * 32 * 32;    
    return ret;
}
*/


__device__ __host__ void compute_pbc_shift( double A[3], double B[3], double * cell, double * shift ){
   double * h_mat = &cell[CELL_HMAT_OFF];
   double * h_inv = &cell[CELL_HINV_OFF];
   double AB[3];
   AB[0] = A[0]-B[0];
   AB[1] = A[1]-B[1];
   AB[2] = A[2]-B[2];
   // note it is a 3x3 by 3 matrix vector product
   double s0 = h_inv[0*3+0] * AB[0] + h_inv[1*3+0] * AB[1] + h_inv[2*3+0] * AB[2] ;
   double s1 = h_inv[0*3+1] * AB[0] + h_inv[1*3+1] * AB[1] + h_inv[2*3+1] * AB[2] ;
   double s2 = h_inv[0*3+2] * AB[0] + h_inv[1*3+2] * AB[1] + h_inv[2*3+2] * AB[2] ;
   s0 = rint( s0 );
   s1 = rint( s1 );
   s2 = rint( s2 );
   // note it is a 3x3 by 3 matrix vector product
   shift[0] = h_mat[0*3+0] * s0 + h_mat[1*3+0] * s1 + h_mat[2*3+0] * s2;
   shift[1] = h_mat[0*3+1] * s0 + h_mat[1*3+1] * s1 + h_mat[2*3+1] * s2;
   shift[2] = h_mat[0*3+2] * s0 + h_mat[1*3+2] * s1 + h_mat[2*3+2] * s2;
}


__device__ __host__ void decode_ipabcd_none(
      unsigned int ipzn,
      unsigned int* __restrict__ ipa, unsigned int* __restrict__ ipb,
      unsigned int* __restrict__ ipc, unsigned int* __restrict__ ipd ){
   (*ipa) = ipzn / (16*16*16*32*32*32);
   (*ipb) = ipzn / (16*16*16*32*32)%32;
   (*ipc) = ipzn / (16*16*16*32)%32;
   (*ipd) = ipzn / (16*16*16)%32;
}


__host__ __device__ void decode_ipabcd_n123(
      unsigned int ipzn,
      unsigned int* __restrict__ ipa, unsigned int* __restrict__ ipb,
      unsigned int* __restrict__ ipc, unsigned int* __restrict__ ipd,
      int* __restrict__ n1, int* __restrict__ n2, int* __restrict__ n3 ){
   (*ipa) = ipzn / (16*16*16*32*32*32);
   (*ipb) = ipzn / (16*16*16*32*32)%32;
   (*ipc) = ipzn / (16*16*16*32)%32;
   (*ipd) = ipzn / (16*16*16)%32;
   (*n1)  = ipzn / (16*16)%16-8;
   (*n2)  = ipzn / (16)%16-8;
   (*n3)  = ipzn %16-8;
}




// #### device L ####

__constant__ int _NLco_lut_dev[35] = { 0, 0, 1, 3, 6, 10, 15, 21, 28, 36, 45, 55, 66, 78, 91, 105, 120, 136, 153, 171, 190, 210, 231, 253, 276, 300, 325, 351, 378, 406, 435, 465, 496, 528, 561 };
__device__ int NLco_dev( int L ){ return _NLco_lut_dev[L+2]; }

// essentialy is using the pattern:
// s = 0 0 0
// p = 1 0 0, 0 1 0, 0 0 1
// d = 2 0 0, 1 1 0, 1 0 1, 0 2 0, 0 1 1, 0 0 2
// and noting that L-lx does not really depend on L
__constant__ short int lx_lut_dev[45] = { 0, 1,1, 2,2,2, 3,3,3,3, 4,4,4,4,4, 5,5,5,5,5,5, 6,6,6,6,6,6,6, 7,7,7,7,7,7,7,7, 8,8,8,8,8,8,8,8,8 };

// compute (cartesian) moment on x axis for a given total moment.
__device__ int lx_dev( const int i, const int L ){
   return L - lx_lut_dev[i];
}

// 
__device__ int lz_dev( const int i, const int L ){
   int i0 = NLco_dev(lx_lut_dev[i]-1);
   int lz_ = i - i0;
   return lz_;
}

// computes ly as L-lx-lz
__device__ int ly_dev( const int i, const int L ){
   int lx_ = lx_dev(i,L);
   int i0 = NLco_dev(lx_lut_dev[i]-1);
   int lz_ = i - i0;
   int ly_ = L-lx_-lz_;
   return ly_;
}


// #### host L ####

int NLco( int L ){
   const int _NLco_lut[35] = { 0, 0, 1, 3, 6, 10, 15, 21, 28, 36, 45, 55, 66, 78, 91, 105, 120, 136, 153, 171, 190, 210, 231, 253, 276, 300, 325, 351, 378, 406, 435, 465, 496, 528, 561 };
   return _NLco_lut[L+2];
}

int L_lx(const int i ){
   const short int lx_lut[45] = { 0, 1,1, 2,2,2, 3,3,3,3, 4,4,4,4,4, 5,5,5,5,5,5, 6,6,6,6,6,6,6, 7,7,7,7,7,7,7,7, 8,8,8,8,8,8,8,8,8 };
   return lx_lut[i];
}

int lx( const int i, const int L ){
   return L - L_lx(i);
}

int lz( const int i, const int L ){
   int i0 = NLco(L_lx(i)-1);
   int lz_ = i - i0;
   return lz_;
}

// computes ly as L-lx-lz
int ly( const int i, const int L ){
   int i0 = NLco(L_lx(i)-1);
   int lz_ = i - i0;
   int ly_ = L_lx(i)-lz_;
   return ly_;
}

// #### so far, both __device__ and __host__ are ok ####

__device__ __host__ int compute_Nc( int la, int lb, int lc, int ld ){
   return (la+1)*(la+2) * (lb+1)*(lb+2) * (lc+1)*(lc+2) * (ld+1)*(ld+2) / 16 ;
}


__device__ __host__ int compute_Ns( int la, int lb, int lc, int ld ){
   return (2*la+1) * (2*lb+1) * (2*lc+1) * (2*ld+1) ;
}


//__host__ __device__ void compute_weighted_distance(
//      double X12[3], const double X1[3], const double X2[3],
//      const double c1, const double c2, const double c12 ){
//}

