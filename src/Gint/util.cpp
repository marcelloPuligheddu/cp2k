#include <cassert>
#include <vector>
#include "util.h"
#include "define.h"


unsigned int encode4( int a, int b, int c, int d ){
   return a*256*256*256 + b*256*256 + c*256 + d;
}

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

unsigned int encode_prm( const int ipa, const int ipb, const int ipc, const int ipd, const int n3 ){
    assert(ipa >= 0);
    assert(ipa < MAX_N_PRM);
    assert(ipb >= 0);
    assert(ipb < MAX_N_PRM);
    assert(ipc >= 0);
    assert(ipc < MAX_N_PRM);
    assert(ipd >= 0);
    assert(ipd < MAX_N_PRM);
    assert(n3 >= 0);
    assert(n3 < MAX_N_CELL);

    unsigned int ret = 0;
    ret +=  n3;
    ret += ipd * MAX_N_CELL;
    ret += ipc * MAX_N_CELL * MAX_N_PRM;
    ret += ipb * MAX_N_CELL * MAX_N_PRM * MAX_N_PRM;
    ret += ipa * MAX_N_CELL * MAX_N_PRM * MAX_N_PRM * MAX_N_PRM;
    return ret;
}

__host__ __device__ void decode_prm( 
      const unsigned int prm,
      unsigned int* __restrict__ ipa, unsigned int* __restrict__ ipb,
      unsigned int* __restrict__ ipc, unsigned int* __restrict__ ipd,
      unsigned int* __restrict__ n3 ){
   (*ipa) = (prm / (MAX_N_CELL * MAX_N_PRM * MAX_N_PRM * MAX_N_PRM));
   (*ipb) = (prm / (MAX_N_CELL * MAX_N_PRM * MAX_N_PRM)) % MAX_N_PRM;
   (*ipc) = (prm / (MAX_N_CELL * MAX_N_PRM)) % MAX_N_PRM;
   (*ipd) = (prm / (MAX_N_CELL)) % MAX_N_PRM;
   (*n3)  = (prm) % MAX_N_CELL;
}

unsigned int encode_shell( const int nla, const int nlb, const int nlc, const int nld, const int n1, const int n2 ){
    assert(nla >= 0);
    assert(nla < MAX_N_L);
    assert(nlb >= 0);
    assert(nlb < MAX_N_L);
    assert(nlc >= 0);
    assert(nlc < MAX_N_L);
    assert(nld >= 0);
    assert(nld < MAX_N_L);
    assert(n1 >= 0);
    assert(n1 < MAX_N_CELL);
    assert(n2 >= 0);
    assert(n2 < MAX_N_CELL);

    unsigned int ret = 0;
    ret +=  n2;
    ret +=  n1 * MAX_N_CELL;
    ret += nld * MAX_N_CELL * MAX_N_CELL;
    ret += nlc * MAX_N_CELL * MAX_N_CELL * MAX_N_L;
    ret += nlb * MAX_N_CELL * MAX_N_CELL * MAX_N_L * MAX_N_L;
    ret += nla * MAX_N_CELL * MAX_N_CELL * MAX_N_L * MAX_N_L * MAX_N_L;
    return ret;
}

__host__ __device__ void decode_shell(
      const unsigned int shell,
      unsigned int* __restrict__ nla, unsigned int* __restrict__ nlb,
      unsigned int* __restrict__ nlc, unsigned int* __restrict__ nld,
      unsigned int* __restrict__ n1,  unsigned int* __restrict__ n2 ){
   (*nla) = (shell / (MAX_N_CELL * MAX_N_CELL * MAX_N_L * MAX_N_L * MAX_N_L));
   (*nlb) = (shell / (MAX_N_CELL * MAX_N_CELL * MAX_N_L * MAX_N_L)) % MAX_N_L;
   (*nlc) = (shell / (MAX_N_CELL * MAX_N_CELL * MAX_N_L)) % MAX_N_L;
   (*nld) = (shell / (MAX_N_CELL * MAX_N_CELL)) % MAX_N_L;
   (*n1)  = (shell / MAX_N_CELL) % MAX_N_CELL;
   (*n2)  = (shell) % MAX_N_CELL;
}

int max( std::vector<int> x ){
   if ( x.size() == 0 ){ return 0; };
   int ret = x[0];
   for( int idx=1; idx<x.size(); idx++ ){ ret = max(ret, x[idx]); }
   return ret;
}

__device__ __host__ void compute_pbc_shift( const double A[3], const double B[3], const double * cell, double * shift ){
   const double * h_mat = &cell[CELL_HMAT_OFF];
   const double * h_inv = &cell[CELL_HINV_OFF];
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

