#include "libGint.h"

extern "C" {

void * libgint_create_handle () {
   libGint * handle = new libGint() ;
   return (void*) handle ;
}

void libgint_init( void * handle ){ 
   libGint * g_handle = ( libGint * ) handle ;
   g_handle -> init();
}

void libgint_set_Potential_Truncated( void * handle, double R_cut_, double * C0_, int ld_C0_, int C0_size_ ){
   libGint * g_handle = ( libGint * ) handle ;
   g_handle -> set_Potential_Truncated( R_cut_, C0_, ld_C0_, C0_size_ );
}

void libgint_set_P( void * handle, double * P, int P_size ){
   libGint * g_handle = (libGint *) handle ;
   cout << " Setting P with handle " << handle << "|" << g_handle << " @ " << P << " x " << P_size << endl;
   g_handle -> set_P( P, P_size );
}

void libgint_set_P_polarized( void * handle, double * Pa, double * Pb, int P_size ){
   libGint * g_handle = (libGint *) handle ;
   g_handle -> set_P( Pa, Pb, P_size ); 
}

void libgint_set_K( void * handle, double * K, int K_size, double fac ){
   libGint * g_handle = (libGint *) handle ;
   cout << " Setting K with handle " << handle << "|" << g_handle << " @ " << K << " x " << K_size << " fac: " << fac << endl;
   g_handle -> set_K( K, K_size ); 
   g_handle -> hf_fac = fac ;
}

void libgint_set_K_polarized( void * handle, double * Ka, double * Kb, int K_size ){
   libGint * g_handle = (libGint *) handle ;
   g_handle -> set_K( Ka, Kb, K_size ); 
}

void libgint_get_K( void * handle, double * K ){
   libGint * g_handle = (libGint *) handle ;
   g_handle -> get_K( K ); 
}

void libgint_get_K_polarized( void * handle, double * Ka, double * Kb ){
   libGint * g_handle = (libGint *) handle ;
   g_handle -> get_K( Ka, Kb ); 
}

void libgint_set_cell( void * handle, bool periodic, double * cell_h, double * cell_i ){
   libGint * g_handle = (libGint *) handle ;
   g_handle -> set_cell( periodic, cell_h, cell_i );
}

void libgint_set_neighs( void * handle, double * neighs_, int nneighs ){
   libGint * g_handle = (libGint *) handle ;
   g_handle -> set_neighs( neighs_, nneighs );
}

void libgint_set_Atom( void * handle, int i, double * R, double * Z, int np ){
   libGint * g_handle = (libGint *) handle ;

   cout << " Setting Atom(set) " << i << " at " << R[0] << " " << R[1] << " " << R[2] << " w Z: | " ;
   for ( int i = 0 ; i < np ; i++ ){ cout << Z[i] << " " ; }
   cout << " | " << np << endl;

   g_handle -> set_Atom( i, R, Z, np );
}



double compute_norm_psi( double * K, double * z, int np, int l ){
   double ans = 0.0;
   for ( int i=0; i < np; i++ ){
      for ( int j=0; j < np; j++ ){
         ans += K[i] * K[j] / pow(z[i]+z[j],1.5+l);
      }
   }
   ans *= pow(M_PI,1.5) / pow(2.0,l);
//   cout << " norm at l " << l << " of " << K[0] << " ^ " << z[0] << " [." << np << ".] = " << K[0] / sqrt(ans) << endl ;
   return sqrt(ans);

}

#include <vector>
void libgint_set_Atom_L( void * handle, int i, int l, int nl, double * K ){
   libGint * g_handle = (libGint *) handle ;

   int np =  g_handle -> np[i];
   double * z = &( g_handle -> ua.internal_buffer[g_handle -> idx_Z[i]] );

   std::vector<double> K_normalized( nl * np );
   for ( int inl=0; inl < nl ; inl++ ){
      double norm = compute_norm_psi( &K[inl * np], z, np, l );
      for ( int ip = 0; ip < np ; ip ++ ){
         K_normalized[inl * np + ip] = K[inl * np + ip] / norm ;
      }
   }

   cout << " Setting Atom L (set) " << " as " << i << " " << l << " " << nl << "x" << np << " K: " << endl;
   for ( int inl=0; inl < nl ; inl++ ){
      for ( int ip = 0; ip < np ; ip ++ ){
         cout << K[ inl * np + ip ] << " ( " <<  K_normalized[inl * np + ip] << " ) " ;
      } cout << endl ;
   } cout << endl;



   g_handle -> set_Atom_L( i, l, nl, K_normalized.data() );
}

void libgint_set_AtomInfo( void * handle, int i, double * R, double * Z, int np, int lmin, int Lmax, int * nl, double * K ){
   libGint * g_handle = (libGint *) handle ;
   g_handle -> set_Atom( i, R, Z, np );
   int k_offset = 0 ;
   for ( int l = lmin; l < Lmax; l++ ){
      g_handle -> set_Atom_L( i, l, nl[i-lmin], &K[k_offset] );
      k_offset += nl[i-lmin] * np;
   }
}

void libgint_set_max_n_cell( void * handle , int max_n_cell_ ){
   libGint * g_handle = (libGint *) handle ;
   g_handle -> set_max_n_prm( max_n_cell_ );
}

void libgint_add_prm( void * handle, int ipa, int ipb, int ipc, int ipd, int n3 ){
   libGint * g_handle = (libGint *) handle ;
   g_handle -> add_prm( ipa, ipb, ipc, ipd, n3 );
}

void libgint_add_shell( void * handle, int i, int j ,int k, int l, int n1, int n2 ){
   libGint * g_handle = (libGint *) handle ;
   g_handle -> add_shell( i,j,k,l,n1,n2 );
}

void libgint_add_cell( void * handle ){
   libGint * g_handle = (libGint *) handle ;
   g_handle -> add_cell();
}

void libgint_add_qrt( void * handle, int la, int lb, int lc, int ld, int nla, int nlb, int nlc, int nld ){
   libGint * g_handle = (libGint *) handle ;
   g_handle -> add_qrt( la, lb, lc, ld, nla, nlb, nlc, nld );
}

void libgint_add_qrtt(
      void * handle, double symm_fac, int la, int lb, int lc, int ld,
      int inla, int inlb, int inlc, int inld, 
      int ld_ac, int ld_ad, int ld_bc, int ld_bd, 
      unsigned int offset_ac_L_set, unsigned int offset_ad_L_set,
      unsigned int offset_bc_L_set, unsigned int offset_bd_L_set,
      bool Tac, bool Tad, bool Tbc, bool Tbd ){
   libGint * g_handle = (libGint *) handle ;
   g_handle -> add_qrtt( symm_fac, la,lb,lc,ld, inla,inlb,inlc,inld, ld_ac,ld_ad,ld_bc,ld_bd,
                         offset_ac_L_set,offset_ad_L_set,offset_bc_L_set,offset_bd_L_set,
                         Tac, Tad, Tbc, Tbd );
}

void libgint_add_set( void * handle ){
   libGint * g_handle = (libGint *) handle ;
   g_handle -> add_set();
}

void libgint_memory_needed( void * handle, int * mem ){
   libGint * g_handle = (libGint *) handle ;
   (*mem) = ( g_handle -> memory_needed() );
}

void libgint_dispatch( void * handle ){
   libGint * g_handle = (libGint *) handle ;
   g_handle -> dispatch( true );
}

} // end of extern "C"

