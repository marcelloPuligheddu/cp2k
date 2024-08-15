#include <vector>
#include <iostream>
#include "define.h"
#include "util.h"

using std::cout;
using std::endl;

__global__ void compute_TRA_batched_gpu_low(
   const int Nshells, const int la, const int lb, const int lc, const int ld,
   const unsigned int* const __restrict__ TRA,
   const double* const __restrict__ SPHER,
   double* const __restrict__ OUT ){

   for ( int block=blockIdx.x ; block < Nshells; block+=gridDim.x ){
      int Oq  = TRA[block*TRA_SIZE+TRA_OFFSET_Q];
      int nla = TRA[block*TRA_SIZE+TRA_OFFSET_NLA];
      int nlb = TRA[block*TRA_SIZE+TRA_OFFSET_NLB];
      int nlc = TRA[block*TRA_SIZE+TRA_OFFSET_NLC];
      int nld = TRA[block*TRA_SIZE+TRA_OFFSET_NLD];
      int Od  = TRA[block*TRA_SIZE+TRA_OFFSET_DEST];

      int nsoa = compute_Ns(la);
      int nsob = compute_Ns(lb);
      int nsoc = compute_Ns(lc);
      int nsod = compute_Ns(ld);
      int Ns = compute_Ns(la,lb,lc,ld);
      int n = nla*nlb*nlc*nld * Ns;
      int spher_off = Oq * Ns;

      if (nla*nlb*nlc*nld == 1){
         for( int is=threadIdx.x; is < Ns; is+=blockDim.x ){
            OUT[ Od + is ] = SPHER[ spher_off + is ];
         }
      } else {
         for ( int thread=threadIdx.x; thread < n; thread+=blockDim.x ){
            // SPHER is ordered as ila,ilb,ilc,ild,isa,isb,isc,isd
            // OUT is ordered as ila,isa,ilb,isb,ilc,isc,ild,isd
            int d = 0;
            int i = thread;
            int nel = 1;
            nel *=  1   ; int isd = ( (i-d) / nel ) % nsod ; d += isd * nel ;
            nel *= nsod ; int isc = ( (i-d) / nel ) % nsoc ; d += isc * nel ;
            nel *= nsoc ; int isb = ( (i-d) / nel ) % nsob ; d += isb * nel ;
            nel *= nsob ; int isa = ( (i-d) / nel ) % nsoa ; d += isa * nel ;
            nel *= nsoa ; int ild = ( (i-d) / nel ) % nld  ; d += ild * nel ;
            nel *= nld  ; int ilc = ( (i-d) / nel ) % nlc  ; d += ilc * nel ;
            nel *= nlc  ; int ilb = ( (i-d) / nel ) % nlb  ; d += ilb * nel ;
            nel *= nlb  ; int ila = ( (i-d) / nel ) % nla  ; d += ila * nel ;
            int j = ((((((ila*nsoa+isa)*nlb+ilb)*nsob+isb)*nlc+ilc)*nsoc+isc)*nld+ild)*nsod+isd;
            double val = SPHER[ spher_off + i ];
            OUT[ Od + j ] = val;
         }
      }
   }
}


void compute_TRA_batched_low(
   const int Nshells, const int la, const int lb, const int lc, const int ld,
   const unsigned int* const __restrict__ TRA,
   const double* const __restrict__ SPHER,
   double* const __restrict__ OUT ){

   for ( int block=0 ; block < Nshells; block++ ){
      int Oq  = TRA[block*TRA_SIZE+TRA_OFFSET_Q];
      int nla = TRA[block*TRA_SIZE+TRA_OFFSET_NLA];
      int nlb = TRA[block*TRA_SIZE+TRA_OFFSET_NLB];
      int nlc = TRA[block*TRA_SIZE+TRA_OFFSET_NLC];
      int nld = TRA[block*TRA_SIZE+TRA_OFFSET_NLD];
      int Od  = TRA[block*TRA_SIZE+TRA_OFFSET_DEST];

      int nsoa = compute_Ns(la);
      int nsob = compute_Ns(lb);
      int nsoc = compute_Ns(lc);
      int nsod = compute_Ns(ld);
      int Ns = compute_Ns(la,lb,lc,ld);
      int n = nla*nlb*nlc*nld * Ns;
      int spher_off = Oq * Ns;

//      cout << "TRA " << block << " " << Oq << " " << nla << " " << nlb << " " << nlc << " " << nld << " " << Od << endl;
      if (nla*nlb*nlc*nld == 1){
         for( int is=0; is < Ns; is++ ){
//            cout << "OUT["<<Od+is<<"] = " << SPHER[ spher_off + is ] << " @ " << spher_off + is << endl;
            OUT[ Od + is ] = SPHER[ spher_off + is ];
         }
      } else {
         for ( int thread=0; thread < n; thread++ ){
            // SPHER is ordered as ila,ilb,ilc,ild,isa,isb,isc,isd
            // OUT is ordered as ila,isa,ilb,isb,ilc,isc,ild,isd
            int d = 0;
            int i = thread;
            int nel = 1;
            nel *=  1   ; int isd = ( (i-d) / nel ) % nsod ; d += isd * nel ;
            nel *= nsod ; int isc = ( (i-d) / nel ) % nsoc ; d += isc * nel ;
            nel *= nsoc ; int isb = ( (i-d) / nel ) % nsob ; d += isb * nel ;
            nel *= nsob ; int isa = ( (i-d) / nel ) % nsoa ; d += isa * nel ;
            nel *= nsoa ; int ild = ( (i-d) / nel ) % nld  ; d += ild * nel ;
            nel *= nld  ; int ilc = ( (i-d) / nel ) % nlc  ; d += ilc * nel ;
            nel *= nlc  ; int ilb = ( (i-d) / nel ) % nlb  ; d += ilb * nel ;
            nel *= nlb  ; int ila = ( (i-d) / nel ) % nla  ; d += ila * nel ;
            int j = ((((((ila*nsoa+isa)*nlb+ilb)*nsob+isb)*nlc+ilc)*nsoc+isc)*nld+ild)*nsod+isd;
            double val = SPHER[ spher_off + i ];
            OUT[ Od + j ] = val;
         }
      }
   }
}

void compute_TRA_batched(
      const int Nshells, const int la, const int lb, const int lc, const int ld,
      const std::vector<unsigned int>& TRA, const std::vector<double>& SPHER, std::vector<double>& OUT ){
   compute_TRA_batched_low(Nshells, la, lb, lc, ld, TRA.data(), SPHER.data(), OUT.data() );
}






