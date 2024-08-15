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


#include "fgamma.h"

#define Teps 1.e-13
#define NDERIV_MAX 21
#define DEGREE 13

/*
double* read_c0( int Nder, FILE* stream, int* ldc0_ptr ){
   const int patches = 207;
   assert( Nder <= NDERIV_MAX );
   int size1 = ((Nder       + 1)*(DEGREE + 1)*(DEGREE + 2))/2;
   int size2 = ((NDERIV_MAX + 1)*(DEGREE + 1)*(DEGREE + 2))/2;
   int ldc0 = size1;
   double* C0_ptr = (double*)malloc( patches*size1*sizeof(double) );
   for( int i=0; i < patches; i++ ){
      for ( int j=0; j < size2 ; j++ ){
         double temp;
         int nscanned = fscanf( stream , "%lE", &temp);
         assert ( nscanned == 1 );
         if ( j < ldc0 ){
            C0_ptr[i*ldc0+j] = temp;
         }
      }
   }
   (*ldc0_ptr) = ldc0;
   return C0_ptr;
}
*/

void fgamma_ref( int nmax , double T, double* f){
   const int kmax = 50;
   double r[kmax+10];

   if (T < Teps) {
      for( int n = 0 ; n <= nmax ; n++ ){
         f[n] = 1.0/(2.*n + 1.);
      }
   } else if (T <= 50.0) {
      // Initialize ratios of Bessel functions
      r[kmax + 9] = 0.0;
      for ( int j = kmax+9; j >= 1; j--){
         r[j-1] = -T/((4.*j + 2.) - T*r[j+1-1]);
      }
      double factor = 2.0*sinh(0.5*T)*exp(-0.5*T)/T;
      for( int n = 0; n <= nmax ; n++ ){
      // Initialize iteration
         double sumtot = factor/(2.*n + 1.);
         double term = 1.0;

         // Begin the summation and recursion
         for ( int k = 1; k <= kmax ; k++ ){
            term = term*(2.*n - 2.*k + 1.)/(2.*n + 2.*k + 1);

            // Product of Bessel function quotients
            double p = 1.0;
            for ( int j = 1 ; j <= k ; j++ ){
               p = p*r[j-1];
            }
            double sumterm = factor*term*p*(2.*k + 1.)/(2.*n + 1.);
            //
            if (abs(sumterm) < DBL_EPSILON) {
               // Iteration converged
               break;
            } else if (k == kmax) {
               // No convergence with kmax iterations
               assert(k<kmax);
            } else {
               // Add the current term to the sum and continue the iteration
               sumtot = sumtot + sumterm;
            }
         }
         f[n] = sumtot;
      }
   } else {
      //Use asymptotic formula for T > 50
      f[0] = 0.5*sqrt(M_PI/T);
      // Use the upward recursion relation to
      // generate the remaining F_n(t) values
      double expt = exp(-T);
      for( int n = 1 ; n <= nmax ; n++ ){
         f[n] = 0.5*((2.*n-1.)*f[n-1] - expt)/T;
      }
   }
}

double* create_md_ftable( int nmax, double tmin, double tmax, double tdelta, int* ld ){
   int n = nmax + 6;
   int itabmin = int(floor(tmin/tdelta));
   int itabmax = int((tmax - tmin)/tdelta);
   int nelem = (itabmax - itabmin + 1 ) * (n+1);
   double * ftable = (double*) malloc(nelem*sizeof(double) );
   // Fill table
   for ( int itab = itabmin; itab <= itabmax; itab++ ){
      double T = double(itab)*tdelta;
      // fill ftable from itab*(n+1) to (itab+1)*(n+1) with f values
      fgamma_ref(n, T, &ftable[itab*(n+1)]);
   }
   *ld = n+1;
   return ftable;
}

void fgamma0_ref( int nmax, double T, double* f, const double* ftable, int ftable_ld ){
   if ( T < Teps ){
      // eps < T -> T=0
      for( int n = 0 ; n <= nmax ; n++ ){
         f[n] = 1.0/(2.*n+1.);
      }
   } else if ( T < 12.  ){
      // eps < T < 12 -> Taylor expansion + downward
      double tdelta = 0.1;
      int itab = int(round(T/tdelta));
      double ttab = double(itab)*tdelta;
      double tmp = 1.0;
      double inv_fac = 1.0; // 1/i!
      f[nmax] = ftable[itab*ftable_ld+nmax];
      // 6th order Taylor expansion
      tmp *= (ttab-T); inv_fac /= 1.0 ; f[nmax] += ftable[itab*ftable_ld+nmax+1]*tmp*inv_fac;
      tmp *= (ttab-T); inv_fac /= 2.0 ; f[nmax] += ftable[itab*ftable_ld+nmax+2]*tmp*inv_fac;
      tmp *= (ttab-T); inv_fac /= 3.0 ; f[nmax] += ftable[itab*ftable_ld+nmax+3]*tmp*inv_fac;
      tmp *= (ttab-T); inv_fac /= 4.0 ; f[nmax] += ftable[itab*ftable_ld+nmax+4]*tmp*inv_fac;
      tmp *= (ttab-T); inv_fac /= 5.0 ; f[nmax] += ftable[itab*ftable_ld+nmax+5]*tmp*inv_fac;
      tmp *= (ttab-T); inv_fac /= 6.0 ; f[nmax] += ftable[itab*ftable_ld+nmax+6]*tmp*inv_fac;
//      for ( int k=1; k <= 6; k++ ){
//         tmp *= (ttab-T);
//         f[nmax] += ftable[itab*ftable_ld+nmax+k]*tmp*ifac[k];
//      }
      // Use the downward recursion relation to 
      // generate the remaining F_n(t) values   
      double expt = exp(-T);
      for( int n = nmax-1 ; n >= 0 ; n--){
         f[n] = (2.0*T*f[n + 1] + expt)/(2.*n + 1.);
      }
   } else {
   // T > 12 -> segmented Taylor + upward
      double tmp = 1.0/T; // reciprocal value
      if (T <= 15.0) {
         // 12 < T <= 15 -> Four term polynom expansion
         double tmp2 = tmp*tmp;
         double tmp3 = tmp*tmp2;
         double g = 0.4999489092 - 0.2473631686*tmp + 0.321180909*tmp2 - 0.3811559346*tmp3;
         f[0] = 0.5*sqrt(M_PI*tmp) - g*exp(-T)*tmp;
      } else if (T <= 18.0) {
         // 15 < T <= 18 -> Three term polynom expansion
         double g = 0.4998436875 - 0.24249438*tmp + 0.24642845*tmp*tmp;
         f[0] = 0.5*sqrt(M_PI*tmp) - g*exp(-T)*tmp;
      } else if (T <= 24.0) {
         // 18 < T <= 24 -> Two term polynom expansion
         double g = 0.499093162 - 0.2152832*tmp;
         f[0] = 0.5*sqrt(M_PI*tmp) - g*exp(-T)*tmp;
      } else if (T <= 30.0) {
         // 24 < T <= 30 -> One term polynom expansion
         double g = 0.49;
         f[0] = 0.5*sqrt(M_PI*tmp) - g*exp(-T)*tmp;
      } else {
         // T > 30 -> Asymptotic formula
         f[0] = 0.5*sqrt(M_PI*tmp);
      }
      // Use the upward recursion relation to
      // generate the remaining F_n(t) values
      double expt = 0.0;
      if ( T > 2.*nmax + 36.){
         expt = 0.0;
      } else {
         expt = exp(-T);
      }
      for ( int n=1 ; n <= nmax ; n++ ){
         f[n] = (0.5*tmp)*( (2.*n - 1.)*f[n - 1] - expt );
      }
   }
}


__host__ __device__ void fgamma0( int nmax, double T, double* f, const double* ftable, int ftable_ld ){
   if ( T < Teps ){
      // eps < T -> T=0
      for( int n = 0 ; n <= nmax ; n++ ){
         f[n] = 1.0/(2.*n+1.);
      }
   } else if ( T < 12.  ){
      // eps < T < 12 -> Taylor expansion + downward
      double tdelta = 0.1;
      int itab = int(round(T/tdelta));
      double ttab = double(itab)*tdelta;
      double tmp = 1.0;
      double inv_fac = 1.0; // 1/i!
      f[nmax] = ftable[itab*ftable_ld+nmax];
      // 6th order Taylor expansion
      tmp *= (ttab-T); inv_fac = 1.0    ; f[nmax] += ftable[itab*ftable_ld+nmax+1]*tmp*inv_fac;
      tmp *= (ttab-T); inv_fac = 0.5    ; f[nmax] += ftable[itab*ftable_ld+nmax+2]*tmp*inv_fac;
      tmp *= (ttab-T); inv_fac = 1./6.  ; f[nmax] += ftable[itab*ftable_ld+nmax+3]*tmp*inv_fac;
      tmp *= (ttab-T); inv_fac = 1./24. ; f[nmax] += ftable[itab*ftable_ld+nmax+4]*tmp*inv_fac;
      tmp *= (ttab-T); inv_fac = 1./120.; f[nmax] += ftable[itab*ftable_ld+nmax+5]*tmp*inv_fac;
      tmp *= (ttab-T); inv_fac = 1./720.; f[nmax] += ftable[itab*ftable_ld+nmax+6]*tmp*inv_fac;
      // Use the downward recursion relation to 
      // generate the remaining F_n(t) values   
      double expt = exp(-T);
      for( int n = nmax-1 ; n >= 0 ; n--){
         f[n] = (2.0*T*f[n + 1] + expt)/(2.*n + 1.);
      }
   } else {
   // T > 12 -> 4 Taylor + upward
      double tmp = 1.0/T; // reciprocal value
      double tmp2 = tmp*tmp;
      double tmp3 = tmp*tmp2;
      double g = 0.4999489092 - 0.2473631686*tmp + 0.321180909*tmp2 - 0.3811559346*tmp3;
      f[0] = 0.5*sqrt(M_PI*tmp) - g*exp(-T)*tmp;
      // Use the upward recursion relation to
      // generate the remaining F_n(t) values
      double expt = exp(-T);

      for ( int n=1 ; n <= nmax ; n++ ){
         f[n] = (0.5*tmp)*( (2.*n - 1.)*f[n - 1] - expt );
      }
   }
}

#undef Teps
#undef NDERIV_MAX
#undef DEGREE
