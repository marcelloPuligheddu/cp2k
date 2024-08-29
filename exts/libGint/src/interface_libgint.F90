module libgint
   use, intrinsic :: iso_c_binding, only : c_loc, c_ptr, c_null_ptr, c_int, c_bool
   implicit none
   private
   interface
   function libgint_internal_create_handle() result ( libgint_handle ) bind(C,name="libgint_create_handle")
      import :: c_ptr
      implicit none
      type(c_ptr) :: libgint_handle 
   end function libgint_internal_create_handle
   end interface

   interface
   subroutine libgint_internal_init( libgint_handle ) bind(C,name="libgint_init")
      import :: c_ptr
      implicit none
      type(c_ptr), value :: libgint_handle
   end subroutine libgint_internal_init
   end interface

   interface
   subroutine libgint_internal_set_P( libgint_handle, P, P_size ) bind(C,name="libgint_set_P")
      import :: c_ptr, c_int
      implicit none
      type(c_ptr), value :: libgint_handle, P
      integer(kind=c_int), value :: P_size
   end subroutine libgint_internal_set_P
   end interface

   interface
   subroutine libgint_internal_set_P_polarized( libgint_handle, Pa, Pb, P_size ) bind(C,name="libgint_set_P_polarized")
      import :: c_ptr, c_int
      implicit none
      type(c_ptr), value :: libgint_handle, Pa, Pb
      integer(kind=c_int), value :: P_size
   end subroutine libgint_internal_set_P_polarized
   end interface

   interface
   subroutine libgint_internal_set_K( libgint_handle, K, K_size, fac ) bind(C,name="libgint_set_K")
      import :: c_ptr, c_int
      implicit none
      type(c_ptr), value :: libgint_handle, K
      integer(kind=c_int), value :: K_size
      real(kind=8), value :: fac
   end subroutine libgint_internal_set_K
   end interface

   interface
   subroutine libgint_internal_set_K_polarized( libgint_handle, Ka, Kb, K_size ) bind(C,name="libgint_set_K_polarized")
      import :: c_ptr, c_int
      implicit none
      type(c_ptr), value :: libgint_handle, Ka, Kb
      integer(kind=c_int), value :: K_size
   end subroutine libgint_internal_set_K_polarized
   end interface

   interface
   subroutine libgint_internal_get_K( libgint_handle, K ) bind(C,name="libgint_get_K")
      import :: c_ptr, c_int
      implicit none
      type(c_ptr), value :: libgint_handle, K
   end subroutine libgint_internal_get_K
   end interface

   interface
   subroutine libgint_internal_get_K_polarized( libgint_handle, Ka, Kb ) bind(C,name="libgint_get_K_polarized")
      import :: c_ptr, c_int
      implicit none
      type(c_ptr), value :: libgint_handle, Ka, Kb
   end subroutine libgint_internal_get_K_polarized
   end interface

   interface 
   subroutine libgint_internal_set_Atom( libgint_handle, i, R, Z, np ) bind(C,name="libgint_set_Atom")
      import :: c_ptr, c_int
      implicit none
      integer(kind=c_int), value :: i, np
      type (c_ptr), value :: libgint_handle, R ,Z
   end subroutine libgint_internal_set_Atom
   end interface

   interface 
   subroutine libgint_internal_set_Atom_L( libgint_handle, i, l, nl, K ) bind(C,name="libgint_set_Atom_L")
      import :: c_ptr, c_int
      implicit none
      integer(kind=c_int), value :: i, l, nl
      type (c_ptr), value :: libgint_handle, K
   end subroutine libgint_internal_set_Atom_L
   end interface

   interface
   subroutine libgint_internal_set_AtomInfo( libgint_handle, i, R, Z, np, lmin, Lmax, nl, K ) bind(C,name="libgint_set_AtomInfo")
      import c_ptr, c_int
      implicit none
      integer(kind=c_int), value :: i, np, lmin, Lmax
      type(c_ptr), value :: libgint_handle, R,Z,nl,K
   end subroutine libgint_internal_set_AtomInfo
   end interface

   interface
   subroutine libgint_internal_set_cell( libgint_handle, periodic, cell_h, cell_i ) bind(C,name="libgint_set_cell")
      import :: c_ptr, c_int, c_bool
      implicit none
      logical(kind=c_bool) :: periodic
      type(c_ptr), value :: libgint_handle, cell_h,cell_i
   end subroutine libgint_internal_set_cell
   end interface

! void libgint_set_neighs( void * handle, double * neighs_, int nneighs ){
   interface
   subroutine libgint_internal_set_neighs( libgint_handle, neighs, nneighs ) bind(C,name="libgint_set_neighs")
      import :: c_ptr, c_int
      implicit none
      type(c_ptr), value :: libgint_handle, neighs
      integer(kind=c_int), value :: nneighs
   end subroutine libgint_internal_set_neighs
   end interface

   interface
   subroutine libgint_internal_set_max_n_cell( libgint_handle, max_n ) bind(C,name="libgint_set_max_n_cell")
      import :: c_ptr, c_int
      implicit none
      integer(kind=c_int), value :: max_n
      type (c_ptr), value :: libgint_handle
   end subroutine libgint_internal_set_max_n_cell
   end interface

   interface
   subroutine libgint_internal_add_prm( libgint_handle, ipa,ipb,ipc,ipd,n3 ) bind (C,name="libgint_add_prm")
      import :: c_ptr, c_int
      implicit none
      integer(kind=c_int), value :: ipa,ipb,ipc,ipd,n3
      type (c_ptr), value :: libgint_handle
   end subroutine libgint_internal_add_prm
   end interface

   interface   
   subroutine libgint_internal_add_shell( libgint_handle, i,j,k,l,n1,n2 ) bind (C,name="libgint_add_shell")
      import :: c_ptr, c_int
      implicit none
      integer(kind=c_int), value :: i,j,k,l,n1,n2
      type (c_ptr), value :: libgint_handle
   end subroutine libgint_internal_add_shell
   end interface

   interface
   subroutine libgint_internal_add_cell( libgint_handle ) bind (C,name="libgint_add_cell")
      import :: c_ptr
      implicit none
      type (c_ptr), value :: libgint_handle
   end subroutine libgint_internal_add_cell
   end interface

   interface
   subroutine libgint_internal_add_qrt( libgint_handle, la,lb,lc,ld,nla,nlb,nlc,nld ) bind (C,name="libgint_add_qrt")
      import :: c_ptr, c_int
      implicit none
      integer(kind=c_int), value :: la,lb,lc,ld,nla,nlb,nlc,nld
      type (c_ptr), value :: libgint_handle
   end subroutine libgint_internal_add_qrt
   end interface

   interface
   subroutine libgint_internal_add_qrtt( libgint_handle, symm_fac, la,lb,lc,ld, inla,inlb,inlc,inld, &
         ld_ac,ld_ad,ld_bc,ld_bd, offset_ac_L_set,offset_ad_L_set, &
         offset_bc_L_set,offset_bd_L_set, Tac,Tad,Tbc,Tbd ) bind (C,name="libgint_add_qrtt")
      import :: c_ptr, c_int, c_bool
      implicit none
      real( kind=8), value :: symm_fac
      integer(kind=c_int), value :: la,lb,lc,ld, inla,inlb,inlc,inld
      integer(kind=c_int), value :: ld_ac,ld_ad,ld_bc,ld_bd
      integer(kind=c_int), value :: offset_ac_L_set, offset_ad_L_set
      integer(kind=c_int), value :: offset_bc_L_set, offset_bd_L_set
      logical(kind=c_bool), value :: Tac,Tad,Tbc,Tbd
      type (c_ptr), value :: libgint_handle
   end subroutine libgint_internal_add_qrtt
   end interface

   interface
   subroutine libgint_internal_add_set( libgint_handle ) bind (C,name="libgint_add_set")
      import :: c_ptr
      implicit none
      type (c_ptr), value :: libgint_handle
   end subroutine libgint_internal_add_set
   end interface

   interface
   subroutine libgint_internal_memory_needed( libgint_handle, mem ) bind (C,name="libgint_memory_needed")
      import :: c_ptr, c_int
      implicit none
      integer(kind=c_int) :: mem
      type (c_ptr), value :: libgint_handle
   end subroutine libgint_internal_memory_needed
   end interface

   interface
   subroutine libgint_internal_dispatch ( libgint_handle ) bind (C,name="libgint_dispatch")
      import :: c_ptr
      implicit none
      type (c_ptr), value :: libgint_handle
   end subroutine libgint_internal_dispatch
   end interface

   interface 
   subroutine libgint_internal_set_Potential_Truncated( libgint_handle, R_cut, C0, ld_C0, C0_size ) &
                                                      bind (C,name="libgint_set_Potential_Truncated")
      import :: c_ptr, c_int
      implicit none
      type (c_ptr), value :: libgint_handle, C0
      real(kind=8), value :: R_cut
      integer(kind=c_int), value :: ld_C0, C0_size
   end subroutine libgint_internal_set_Potential_Truncated
   end interface



   public :: libgint_init, libgint_set_Potential_Truncated
   public :: libgint_set_P, libgint_set_P_polarized, libgint_set_K, libgint_set_K_polarized
   public :: libgint_get_K, libgint_get_K_polarized, libgint_set_Atom, libgint_set_Atom_L, libgint_set_cell
   public :: libgint_set_neighs
   public :: libgint_set_max_n_cell, libgint_add_prm, libgint_add_shell, libgint_add_cell, libgint_add_qrt
   public :: libgint_add_qrtt, libgint_add_set, libgint_memory_needed, libgint_dispatch

contains
   subroutine libgint_init ( handle )
      type (c_ptr) :: handle
      handle = libgint_internal_create_handle()
      call libgint_internal_init( handle )
   end subroutine libgint_init

   subroutine libgint_set_Potential_Truncated( handle, R_cut, C0 )
      type (c_ptr) :: handle
      real(kind=8) :: R_cut
      integer(kind=c_int) :: ld_C0, C0_size
      real(kind=8), dimension(:,:), target :: C0
      C0_size = size(C0,kind=c_int)
      ld_C0 = size(C0, dim=1, kind=c_int)
      call libgint_internal_set_Potential_Truncated( handle, %val(R_cut), c_loc(C0), %val(ld_C0), %val(C0_size) )
   end subroutine libgint_set_Potential_Truncated

   subroutine libgint_set_P ( handle, P )
      type (c_ptr) :: handle
      integer(kind=c_int) :: P_size
      real(kind=8), dimension(:), target :: P
      P_size = size(P, kind=c_int )
      call libgint_internal_set_P( handle, c_loc(P), %val(P_size) )
   end subroutine libgint_set_P
   
   subroutine libgint_set_P_polarized ( handle, Pa, Pb )
      type (c_ptr) :: handle
      integer(kind=c_int) :: P_size
      real(kind=8), dimension(:), target :: Pa, Pb
      P_size = size(Pa, kind=c_int )
      call libgint_internal_set_P_polarized( handle, c_loc(Pa), c_loc(Pb), %val(P_size) )
   end subroutine libgint_set_P_polarized

   subroutine libgint_set_K ( handle, K, fac )
      type (c_ptr) :: handle
      integer(kind=c_int) :: K_size
      real(kind=8), dimension(:), target :: K
      real(kind=8), value :: fac
      K_size = size(K, kind=c_int )
      call libgint_internal_set_K( handle, c_loc(K), %val(K_size), fac )
   end subroutine libgint_set_K
   
   subroutine libgint_set_K_polarized ( handle, Ka, Kb )
      type (c_ptr) :: handle
      integer(kind=c_int) :: K_size
      real(kind=8), dimension(:), target :: Ka, Kb
      K_size = size(Ka, kind=c_int )
      call libgint_internal_set_K_polarized( handle, c_loc(Ka), c_loc(Kb), %val(K_size) )
   end subroutine libgint_set_K_polarized

   subroutine libgint_get_K ( handle, K )
      type (c_ptr) :: handle
      real(kind=8), dimension(:), target :: K
      call libgint_internal_get_K( handle, c_loc(K) )
   end subroutine libgint_get_K

   subroutine libgint_get_K_polarized ( handle, Ka, Kb )
      type (c_ptr) :: handle
      real(kind=8), dimension(:), target :: Ka, Kb
      call libgint_internal_get_K_polarized( handle, c_loc(Ka), c_loc(Kb) )
   end subroutine libgint_get_K_polarized

   subroutine libgint_set_cell( handle, periodic, cell_h, cell_i )
      type (c_ptr) :: handle
      logical(kind=c_bool), value :: periodic
      real(kind=8), dimension(:,:), target :: cell_h, cell_i
      call libgint_internal_set_cell( handle, periodic, c_loc(cell_h), c_loc(cell_i) )
   end subroutine libgint_set_cell

   subroutine libgint_set_neighs( handle, neighs, nneighs )
      type (c_ptr) :: handle
      integer(kind=c_int) :: nneighs
      real(kind=8), dimension(:,:), target :: neighs
      call libgint_internal_set_neighs( handle, c_loc(neighs), nneighs )
   end subroutine libgint_set_neighs

   subroutine libgint_set_Atom( handle, i, R, Z, np )
      type (c_ptr) :: handle
      integer(kind=c_int), value :: i, np
      real(kind=8), dimension(:), target :: R, Z
!      write(*,*) " calling set Atom(set) ", i, R, Z, np
      call libgint_internal_set_Atom( handle, i, c_loc(R), c_loc(Z), np )
   end subroutine libgint_set_Atom

   subroutine libgint_set_Atom_L( handle, i, l, nl, K )
      type (c_ptr) :: handle
      integer(kind=c_int), value :: i, l , nl
      real(kind=8), dimension(:,:), target :: K
!      write(*,*) " calling set Atom L ", i, l, nl , K
      call libgint_internal_set_Atom_L( handle, i, l, nl, c_loc(K) )
   end subroutine libgint_set_Atom_L

   subroutine ligbint_set_AtomInfo( handle, i, R, Z, np, lmin, Lmax, nl, K )
      type (c_ptr) :: handle
      integer(kind=c_int), value :: i, np, lmin, Lmax
      integer(kind=4), dimension(:), target :: nl
      real(kind=8), dimension(:), target :: R,Z,K
      call libgint_internal_set_AtomInfo( handle, i, c_loc(R), c_loc(Z), np, lmin, Lmax, c_loc(nl) , c_loc(K) )
   end subroutine ligbint_set_AtomInfo

   subroutine libgint_set_max_n_cell( handle, max_periodic_cells )
      type (c_ptr) :: handle
      integer(kind=c_int), value :: max_periodic_cells
      call libgint_internal_set_max_n_cell( handle, max_periodic_cells )
   end subroutine libgint_set_max_n_cell

   subroutine libgint_add_prm( handle, ipa, ipb, ipc, ipd, n3 )
      type (c_ptr) :: handle
      integer(kind=c_int), value :: ipa,ipb,ipc,ipd,n3
      call libgint_internal_add_prm( handle, ipa,ipb,ipc,ipd,n3 )
   end subroutine libgint_add_prm

   subroutine libgint_add_shell( handle, i, j, k, l, n1, n2 )
      type (c_ptr) :: handle
      integer(kind=c_int), value :: i,j,k,l,n1,n2
      call libgint_internal_add_shell( handle, i,j,k,l,n1,n2 )
   end subroutine libgint_add_shell

   subroutine libgint_add_cell( handle )
      type (c_ptr) :: handle
      call libgint_internal_add_cell( handle )
   end subroutine libgint_add_cell

   subroutine libgint_add_qrt( handle, la,lb,lc,ld, nla,nlb,nlc,nld )
      type (c_ptr) :: handle
      integer(kind=c_int), value :: la,lb,lc,ld, nla,nlb,nlc,nld
      call libgint_internal_add_qrt( handle, la,lb,lc,ld,nla,nlb,nlc,nld )
   end subroutine libgint_add_qrt

   subroutine libgint_add_qrtt( handle, symm_fac, la,lb,lc,ld, inla,inlb,inlc,inld, &
         ld_ac,ld_ad,ld_bc,ld_bd, offset_ac_L_set,offset_ad_L_set, &
         offset_bc_L_set,offset_bd_L_set, Tac,Tad,Tbc,Tbd )
      type (c_ptr) :: handle
      real( kind=8), value :: symm_fac
      integer(kind=c_int), value :: la,lb,lc,ld, inla,inlb,inlc,inld
      integer(kind=c_int), value :: ld_ac,ld_ad,ld_bc,ld_bd
      integer(kind=c_int), value :: offset_ac_L_set, offset_ad_L_set 
      integer(kind=c_int), value :: offset_bc_L_set, offset_bd_L_set
      logical(kind=c_bool), value :: Tac,Tad,Tbc,Tbd
      call libgint_internal_add_qrtt( handle, symm_fac, la,lb,lc,ld,inla,inlb,inlc,inld, &
         ld_ac,ld_ad,ld_bc,ld_bd, offset_ac_L_set,offset_ad_L_set, &
         offset_bc_L_set,offset_bd_L_set, Tac,Tad,Tbc,Tbd )
   end subroutine libgint_add_qrtt

   subroutine libgint_add_set( handle )
      type (c_ptr) :: handle
      call libgint_internal_add_set( handle )
   end subroutine libgint_add_set

   function libgint_memory_needed( handle ) result( mem )
      type (c_ptr) :: handle
      integer(kind=c_int) :: mem
      call libgint_internal_memory_needed( handle, mem )
      return
   end function libgint_memory_needed

   subroutine libgint_dispatch( handle )
      type (c_ptr) :: handle
      call libgint_internal_dispatch( handle )
   end subroutine libgint_dispatch

end module libgint























