!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2024 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief Interface to the Libint-Library or a c++ wrapper.
!> \par History
!>      10.2024 Created
!> \author Marcello Puligheddu
! **************************************************************************************************
MODULE libGint_wrapper

   USE kinds, ONLY: dp
#if(__LIBINT)
   use libGint, only : libgint_init, libgint_set_Potential_Truncated, libgint_set_hf_fac, libgint_set_max_mem, &
                       libgint_set_P, libgint_set_P_polarized, libgint_set_K, libgint_set_K_polarized, &
                       libgint_get_K, libgint_get_K_polarized, libgint_set_Atom, libgint_set_Atom_L, &
                       libgint_set_cell, libgint_set_neighs, &
                       libgint_add_prm, libgint_add_shell, libgint_add_cell, libgint_add_qrt, &
                       libgint_add_qrtt, libgint_add_set
#endif
    use hfx_types, only : hfx_type
    use cell_types, only: cell_type
    USE hfx_pair_list_methods, ONLY: build_pair_list_pbc_pgf, bra_t, allocate_bra
    USE particle_types, ONLY: particle_type

#include "./base/base_uses.f90"

   IMPLICIT NONE
   PRIVATE
   INTEGER, ALLOCATABLE, DIMENSION(:) :: first_set_of_atom
   type(bra_t) :: bra, ket
   logical, save :: first_call = .true.
   type(c_ptr), save :: libGint_handle
   !$OMP threadprivate( libGint_handle )
   
   PUBLIC :: cp_libGint_init, libGint_update_env, libGint_set_density, libGint_coulomb4, libGint_update_fock_matrix, libGint_get_fock_matrix
   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'libGint_wrapper'

CONTAINS


subroutine cp_libGint_init( actual_x_data )
  TYPE(hfx_type), POINTER :: actual_x_data
#if(__LINGINT)
   if( first_call) then
      ! Init the offload library, creates an handle unique to this omp thread
      call libgint_init(libGint_handle)
      ! Comunicate the chosen potential and its parameters to libGint
      potential_parameter = actual_x_data%potential_parameter
      libGint_select_potential: select case(potential_parameter%potential_type)
         ! coulomb is the default assumed by libGint, so no more info is needed
         case (do_potential_coulomb)
            exit libGint_select_potential
         case (do_potential_truncated) 
            call libgint_set_Potential_Truncated( libGint_handle, potential_parameter%cutoff_radius, C0(:,:) )
         case default
            CPASSERT(.false.) ! not (yet) implemented
      end select libGint_select_potential
      first_call = .false.
   end if
#else
   MARK_USED(actual_x_data)
   CPABORT("This CP2K executable has not been linked against the required library libGint.")
#endif
end subroutine cp_libGint_init()

subroutine libGint_update_env( fac, memory_parameter, do_periodic, cell, actual_x_data, nneighbors, max_pgf, &
                               natoms, kind_of, particle_set, basis_parameter )

   real(dp) :: fac
   TYPE(hfx_memory_type), POINTER :: memory_parameter
   logical :: do_periodic
   TYPE(cell_type), POINTER :: cell
   TYPE(hfx_type), POINTER :: actual_x_data
   integer :: nneighbors, max_pgf, natoms
   INTEGER, ALLOCATABLE, DIMENSION(:) :: kind_of
   TYPE(particle_type), DIMENSION(:), POINTER         :: particle_set

   logical(1) :: do_pbc
   real(dp), dimension(:,:), allocatable :: cell_r
   real(dp), dimension(:), allocatable :: flat_gcc
   INTEGER, DIMENSION(:), POINTER :: la_mim, la_max, npgfa, nsgfa
   INTEGER, DIMENSION(:, :), POINTER :: nsgfl_a
   integer :: iset,jest, iatom, ikind, nseta, inla
   real(dp) :: ra(3)
   REAL(dp), DIMENSION(:, :), POINTER :: zeta
   REAL(dp), DIMENSION(:, :, : ), POINTER :: gcc

#if(__LINGINT)
   ! Set the multiplicative factor fac (the fraction of exact exchange times the spin factor) for libGint
   call libGint_set_hf_fac( libGint_handle, fac )
   !
   call libGint_set_max_mem( libGint_handle, memory_parameter%max_memory )
   ! Info about periodic cells and neighbouring cells
   do_pbc = do_periodic
   allocate(cell_r(3,nneighbors))
   do i=1, nneighbors
      cell_r(:,i) = actual_x_data%neighbor_cells(i)%cell_r(:)
   end do
   call libgint_set_cell( libGint_handle, do_pbc, cell%hmat, cell%h_inv)
   call libgint_set_neighs( libGint_handle, cell_r, nneighbors) 
   ! CPU side temporary arrays for info about the <AB(g) and CD(n)> pairs
   call allocate_bra( bra, max_pgf, nneighbors )
   call allocate_bra( ket, max_pgf, nneighbors )

   ! Comunicate atomset and atom info to the engine
   jset = 1
   allocate(first_set_of_atom(natom) )
   do iatom = 1, natom
      ikind = kind_of(iatom)
      ra = particle_set(iatom)%r(:)
      nseta = basis_parameter(ikind)%nset
      npgfa => basis_parameter(ikind)%npgf
      la_min => basis_parameter(ikind)%lmin
      la_max => basis_parameter(ikind)%lmax
      zeta => basis_parameter(ikind)%zet
      nsgfa => basis_parameter(ikind)%nsgf
      nsgfl_a => basis_parameter(ikind)%nsgfl
      gcc => basis_parameter(ikind)%gcc
      first_set_of_atom(iatom) = jset
      do iset = 1, nseta
         call libgint_set_Atom( libGint_handle, jset-1 , ra, zeta(:,iset), npgfa(iset) )
         inla = 1
         do l = la_min( iset ) , la_max(iset)
            nla = nsgfl_a(l,iset)
  !                  write (*,*) " CP2K setting atom set L ", jset ,l, nla,  gcc(1:npgfa(iset),inla:inla+nla-1,iset)
  !                  write (*,*) "    note sphi is : " ,  basis_parameter(ikind)%sphi_ext(:, :, :, :)
            IF( ALLOCATED(flat_gcc) )  DEALLOCATE( flat_gcc ) 
            allocate(flat_gcc(npgfa(iset)*nla))
            flat_gcc = pack(gcc(1:npgfa(iset),inla:inla+nla-1,iset),.true.)
            call libgint_set_Atom_L( libGint_handle, jset-1, l, nla, flat_gcc )
            inla = inla + nla
         end do
         jset = jset + 1
      end do 
   end do
#else
   MARK_USED(ALL)
   CPABORT("This CP2K executable has not been linked against the required library libGint.")
#endif
end subroutine libGint_update_env

! Comunicates the current density to the libGint engine.
! Implicitly sets the Fock matrix on the device to zero
interface libGint_set_density
   module subroutine libGint_set_density_A
   module subroutine libGint_set_density_AB
end interface

subroutine libGint_set_density_A( full_density_alpha )
#if(__LIBGINT)
   call libgint_set_P( libGint_handle, full_density_alpha(:, 1) )
#else
   MARK_USED(ALL)
   CPABORT("This CP2K executable has not been linked against the required library libGint.")
#endif

subroutine libGint_set_density_AB( full_density_alpha, full_density_beta )
#if(__LIBGINT)
   call libgint_set_P_polarized( libGint_handle, full_density_alpha(:, 1), full_density_beta(:, 1) )
#else
   MARK_USED(ALL)
   CPABORT("This CP2K executable has not been linked against the required library libGint.")
#endif

interface libGint_get_fock_matrix
   module subroutine libGint_get_fock_matrix_A
   module subroutine libGint_get_fock_matrix_AB
end interface

   subroutine libGint_get_fock_matrix_A( full_ks_alpha_from_gpu(:,1) )
#if(__LIBGINT)
  call libGint_get_K( full_ks_alpha_from_gpu(:,1) )
#else
   MARK_USED(ALL)
   CPABORT("This CP2K executable has not been linked against the required library libGint.")
#endif
   end subroutine update_fock_matrix_gpu

   subroutine libGint_get_fock_matrix_AB( full_ks_alpha_from_gpu(:,1), full_ks_beta_from_gpu(:,1) )
#if(__LIBGINT)
  call libgint_get_K_polarized( libGint_handle, full_ks_alpha_from_gpu(:, 1), full_ks_beta_from_gpu(:, 1) )
#else
   MARK_USED(ALL)
   CPABORT("This CP2K executable has not been linked against the required library libGint.")
#endif
   end subroutine update_fock_matrix_gpu


!> \brief Assign two-electron integrals of a quartet/shell to libGint
   subroutine libGint_coulomb4( iatom, jatom, katom, latom, iset, jset, kset, lset, &
                     ra, rb, rc, rd, npgfa, npgfb, npgfc, npgfd, &
                     potential_parameter, &
                     screen1, screen2, log10_pmax, log10_eps_schwarz, &
                     pgf1, pgf2, &
                     bra, ket, &
                     neighbor_cells, cell, do_periodic, screened )
#if(__LIBGINT)
      type(c_ptr) :: libGint_handle
      integer, intent(in) :: iatom, jatom, katom, latom, iset, jset, kset, lset
      REAL(dp), INTENT(IN)                               :: ra(3), rb(3), rc(3), rd(3)
      INTEGER, INTENT(IN)                                :: npgfa, npgfb, npgfc, npgfd 
      TYPE(hfx_potential_type)                           :: potential_parameter
      TYPE(hfx_cell_type), DIMENSION(:), POINTER         :: neighbor_cells
      REAL(dp), INTENT(IN)                               :: screen1(2), screen2(2)
      TYPE(cell_type), POINTER                           :: cell
      REAL(dp), INTENT(IN)                               :: log10_pmax, log10_eps_schwarz
      TYPE(hfx_screen_coeff_type), DIMENSION(:, :), &
         POINTER                                         :: pgf1, pgf2
      TYPE(bra_t)                                        :: bra, ket
      LOGICAL, INTENT(IN)                                :: do_periodic
      logical, intent(out) :: screened
      logical :: cell_was_screened
      

      integer  :: idx_n1,idx_n2,n1,n2,idx_ij,idx_kl,o_ij,n_ij,o_kl,n_kl !, n3
      integer  :: ipgf,jpgf,kpgf,lpgf,iatom_set,jatom_set,katom_set,latom_set
      real(dp) :: pgf_max_1, pgf_max_2 ! , R1, R2, rpq2

      integer :: nelements_ij, nelements_kl

      cell = cell
      potential_parameter = potential_parameter

      screened = .true.
      iatom_set = first_set_of_atom(iatom)+iset-2
      jatom_set = first_set_of_atom(jatom)+jset-2
      katom_set = first_set_of_atom(katom)+kset-2
      latom_set = first_set_of_atom(latom)+lset-2

      CALL build_pair_list_pbc_pgf( npgfa, npgfb, bra, screen1, screen2, &
                               pgf1, log10_pmax, log10_eps_schwarz, ra, rb, &
                               nelements_ij, neighbor_cells, do_periodic )

      CALL build_pair_list_pbc_pgf( npgfc, npgfd, ket, screen2, screen1, &
                               pgf2, log10_pmax, log10_eps_schwarz, rc, rd, &
                               nelements_kl, neighbor_cells, do_periodic)

      ! Note: we use 3 numbers n1 n2 and n3 as indices for the lattice traslantion vectors
      ! n1 for the AB pair, n2 for the CD pair and n3 for the PQ pair
      ! so that e.g. B = B0 + ua(n1) means B.x = B0.x + ua(n1).x (and same for y and z)
      ! the ua, saved in neighbor_cells(:)%cell_r(:), are already computed 
      ! lattice translation vectors T = i a1 + j a2 + k a3
      ! where a1,a2 and a3 are lattice vectors and i,j and k integers.
      ! So, B.x = B0.x + n1.i * a1.x + n1.j * a2.x + n1.k * a3.x (and same for y and z)
      !  

      ! BIGGER TODO move this to GPU

      do idx_n1 = 1, bra%cell_cnt
      ! TODO use define instead of number for array access
         n1   = bra%cell_idx(1, idx_n1)
         n_ij = bra%cell_idx(2, idx_n1)
         o_ij = bra%cell_idx(3, idx_n1)

         do idx_n2 = 1, ket%cell_cnt

            n2   = ket%cell_idx(1, idx_n2)
            n_kl = ket%cell_idx(2, idx_n2)
            o_kl = ket%cell_idx(3, idx_n2)

            cell_was_screened = .true.
            do idx_ij = o_ij+1, o_ij+n_ij

               ipgf = bra%pgf_idx(1, idx_ij)
               jpgf = bra%pgf_idx(2, idx_ij)
!               n1b  = bra%pgf_idx(3, idx_ij)
               pgf_max_1 = bra%pgf_scr(1, idx_ij)

               do idx_kl = o_kl+1, o_kl+n_kl
                  kpgf = ket%pgf_idx(1,idx_kl)
                  lpgf = ket%pgf_idx(2,idx_kl)
!                  n2b  = ket%pgf_idx(3,idx_kl)
                  pgf_max_2 = ket%pgf_scr(1,idx_kl)

                  IF (pgf_max_1 + pgf_max_2 + log10_pmax < log10_eps_schwarz) CYCLE

                  call libGint_add_prm( libGint_handle, ipgf-1,jpgf-1,kpgf-1,lpgf-1 )
                  cell_was_screened = .false.

               end do ! ket pgf
            end do ! bra pgf

            if ( .not. cell_was_screened ) then
!               write(*,*) 'Adding ', iatom_set, jatom_set, katom_set, latom_set, n1-1, n2-1, 'shell'
               call libgint_add_shell( libGint_handle, iatom_set, jatom_set, katom_set, latom_set, n1-1, n2-1 )
               screened = .false.
            end if

         end do ! ket n2
      end do ! bra n1
#else
   MARK_USED(ALL)
   CPABORT("This CP2K executable has not been linked against the required library libGint.")
#endif
   end subroutine coulomb_4_gpu

   subroutine update_fock_matrix_gpu( &
         libGint_handle, symm_fac, &
         iatom, jatom, katom, latom, &
         iset, jset, kset, lset, &
         atomic_offset_ac, atomic_offset_ad, atomic_offset_bc, atomic_offset_bd, &
         offset_ac_set, offset_ad_set, offset_bc_set, offset_bd_set, &
         nsgfa, nsgfb, nsgfc, nsgfd, &
         la_min, la_max, lb_min, lb_max, &
         lc_min, lc_max, ld_min, ld_max, &
         nsgfl_a, nsgfl_b, nsgfl_c, nsgfl_d)
   type(c_ptr) :: libGint_handle
   real(dp) :: symm_fac
   integer :: iatom, jatom, katom, latom
   integer :: iset, jset, kset, lset
   integer :: atomic_offset_ac, atomic_offset_ad, atomic_offset_bc, atomic_offset_bd
   INTEGER, DIMENSION(:, :), POINTER :: offset_ac_set, offset_ad_set
   INTEGER, DIMENSION(:, :), POINTER :: offset_bc_set, offset_bd_set
   integer :: nsgfa, nsgfb, nsgfc, nsgfd
   integer :: la_min, la_max, lb_min, lb_max
   integer :: lc_min, lc_max, ld_min, ld_max
   INTEGER, DIMENSION(:, :), POINTER :: nsgfl_a, nsgfl_b, nsgfl_c, nsgfl_d

#if(__LIBGINT)
   !! (Hyp)
   ! Let a be a set composed of 2 s and 1 p function.
   ! Let c be a set composed of 1 s and 2 p function.
   ! (1) The density matrix for the ac pair is a 5 x 7 matrix organized as
   !
   !    / -------------------------------------------------------------------------------------------------------------------\
   !    |   a_s1_0@c_s1_0  || a_s1_0@c_p1_0 | a_s1_0@c_p1_1 | a_s1_0@c_p1_2 || a_s1_0@c_p2_0 | a_s1_0@c_p2_1 | a_s1_0@c_p2_2 |
   !    |   a_s2_0@c_s1_0  || a_s2_0@c_p1_0 | a_s2_0@c_p1_1 | a_s2_0@c_p1_2 || a_s2_0@c_p2_0 | a_s2_0@c_p2_1 | a_s2_0@c_p2_2 |
   !    |   a_p1_0@c_s1_0  || a_p1_0@c_p1_0 | a_p1_0@c_p1_1 | a_p1_0@c_p1_2 || a_p1_0@c_p2_0 | a_p1_0@c_p2_1 | a_p1_0@c_p2_2 |
   !    |   a_p1_1@c_s1_0  || a_p1_1@c_p1_0 | a_p1_1@c_p1_1 | a_p1_1@c_p1_2 || a_p1_1@c_p2_0 | a_p1_1@c_p2_1 | a_p1_1@c_p2_2 |
   !    |   a_p1_2@c_s1_0  || a_p1_2@c_p1_0 | a_p1_2@c_p1_1 | a_p1_2@c_p1_2 || a_p1_2@c_p2_0 | a_p1_2@c_p2_1 | a_p1_2@c_p2_2 |
   !    \ -------------------------------------------------------------------------------------------------------------------/
   !
   ! where A_LX_Y means the (Y+1) component of the Xth linear combination of the L angular moment for atom A
   !
   ! (2) This matrix is dense, rectangular and contigous in memory, in fortran column major order.
   ! (3) The big matrix with all pairs is block sparse triangular, only the lower part is valid.


      logical(1) :: Tac,Tad,Tbc,Tbd
      integer :: offset_ac_L_set, offset_ad_L_set, offset_bc_L_set, offset_bd_L_set
      integer :: s_offset_a, s_offset_b, s_offset_c, s_offset_d
      integer :: s_offset_a_l, s_offset_b_l, s_offset_c_l, s_offset_d_l
      integer :: s_offset_ac, s_offset_ad, s_offset_bc, s_offset_bd
      integer :: ld_ac_set, ld_ad_set, ld_bc_set, ld_bd_set
      integer :: la,lb,lc,ld,nla,nlb,nlc,nld,inla,inlb,inlc,inld

      ! TODO rewrite as update_fock_matrix_gpu(libGint_handle, iatomset,jatomset,katomset,latomset )
      ! AFTER TODO communicate (the pointer to) atomic_offset to libGint AND
      ! AFTER TODO communicate (the pointer to)    set_offset to libGint AND
      ! AFTER TODO check if this idea makes sense in general for other codes
      !
      ! Note: this would not change the need to compute sub offsets and
      ! the 8 loops, it would just transfer them to libGint
      ! Except, if libGint can be sure every set has 1  l, it can collapse the l loops
      ! and/or, if libGint can be sure every l   has 1 nl, it can collapse the n loops

      if ( jatom >= latom ) then 
        offset_bd_L_set = offset_bd_set(jset, lset) + atomic_offset_bd - 2
        ld_bd_set = nsgfb
        Tbd = .false.
      else
        offset_bd_L_set = offset_bd_set(lset, jset) + atomic_offset_bd - 2
        ld_bd_set = nsgfd
        Tbd = .true.
      end if
      if ( jatom >= katom ) then
        offset_bc_L_set = offset_bc_set(jset, kset) + atomic_offset_bc - 2
        ld_bc_set = nsgfb
        Tbc = .false.
      else
        offset_bc_L_set = offset_bc_set(kset, jset) + atomic_offset_bc - 2
        ld_bc_set = nsgfc
        Tbc = .true.
      end if

      if ( iatom >= latom ) then
        offset_ad_L_set = offset_ad_set(iset, lset) + atomic_offset_ad - 2
        ld_ad_set = nsgfa
        Tad = .false.
      else
        offset_ad_L_set = offset_ad_set(lset, iset) + atomic_offset_ad - 2
        ld_ad_set = nsgfd
        Tad = .true.
      end if

      if ( iatom >= katom ) then
        offset_ac_L_set = offset_ac_set(iset, kset) + atomic_offset_ac - 2
        ld_ac_set = nsgfa
        Tac = .false.
      else
        offset_ac_L_set = offset_ac_set(kset, iset) + atomic_offset_ac - 2
        ld_ac_set = nsgfc
        Tac = .true.
      end if

!!$omp critical
!      write (*,*) " Atoms ", iatom,jatom,katom,latom,iset,jset,kset,lset
!      write (*,*) " General Offsets: " , atomic_offset_ac, atomic_offset_ad, atomic_offset_bc, atomic_offset_bd
!      write (*,*) "                  " , offset_ac_set, offset_ad_set, offset_bc_set, offset_bd_set
!      write (*,*) "     ->           " , offset_ac_L_set, offset_ad_L_set, offset_bc_L_set, offset_bd_L_set
!!$omp end critical
      
      s_offset_a_l = 0
      DO la = la_min, la_max
         nla = nsgfl_a(la,iset)
         s_offset_b_l = 0
         DO lb = lb_min, lb_max
            nlb = nsgfl_b(lb,jset)
            s_offset_c_l = 0
            DO lc = lc_min, lc_max
               nlc = nsgfl_c(lc,kset)
               s_offset_d_l = 0
               ld_loop: DO ld = ld_min, ld_max
                  nld = nsgfl_d(ld,lset)
                  call libgint_add_qrt(libGint_handle, la,lb,lc,ld,  nla,nlb,nlc,nld )
                  do inla=1, nla
                     s_offset_a = s_offset_a_l + (inla-1) * ( 2 * la + 1 )
                     do inlb=1, nlb
                        s_offset_b = s_offset_b_l + (inlb-1) * ( 2 * lb + 1 )
                        do inlc=1, nlc
                           s_offset_c = s_offset_c_l + (inlc-1) * ( 2 * lc + 1 )
                           do inld=1, nld
                              s_offset_d = s_offset_d_l + (inld-1) * ( 2 * ld + 1 )
                              if (.not. Tac) then
                                 s_offset_ac = offset_ac_L_set + s_offset_c * ld_ac_set + s_offset_a
                              else
                                 s_offset_ac = offset_ac_L_set + s_offset_a * ld_ac_set + s_offset_c
                              end if

                              if (.not. Tad) then
                                 s_offset_ad = offset_ad_L_set + s_offset_d * ld_ad_set + s_offset_a
                              else
                                 s_offset_ad = offset_ad_L_set + s_offset_a * ld_ad_set + s_offset_d
                              end if

                              if (.not. Tbc) then
                                 s_offset_bc = offset_bc_L_set + s_offset_c * ld_bc_set + s_offset_b
                              else
                                 s_offset_bc = offset_bc_L_set + s_offset_b * ld_bc_set + s_offset_c
                              end if

                              if (.not. Tbd) then
                                 s_offset_bd = offset_bd_L_set + s_offset_d * ld_bd_set + s_offset_b
                              else
                                 s_offset_bd = offset_bd_L_set + s_offset_b * ld_bd_set + s_offset_d
                              end if
!!$omp critical
!                              write (*,*) " Atom ", iatom,jatom,katom,latom, iset,jset,kset,lset
!                              write (*,*) " Adding qrtt ", la,lb,lc,ld, inla-1,inlb-1,inlc-1,inld-1
!                              write (*,*) " At offset ",  s_offset_ac,s_offset_ad,s_offset_bc,s_offset_bd
!!$omp end critical
                              call libgint_add_qrtt(libGint_handle, symm_fac, &
                                    la,lb,lc,ld, inla-1,inlb-1,inlc-1,inld-1, &
                                    ld_ac_set, ld_ad_set,ld_bc_set,ld_bd_set, &
                                    s_offset_ac,s_offset_ad,s_offset_bc,s_offset_bd,&
                                    Tac,Tad,Tbc,Tbd )

                           end do
                        end do
                     end do
                  end do
                  s_offset_d_l = s_offset_d_l + nld * ( 2 * ld + 1 )
               end do ld_loop
               s_offset_c_l = s_offset_c_l + nlc * ( 2 * lc + 1 )
            end do
            s_offset_b_l = s_offset_b_l + nlb * ( 2 * lb + 1 )
         end do
         s_offset_a_l = s_offset_a_l + nla * ( 2 * la + 1 )
      end do

      call libgint_add_set( libGint_handle )

#else
   MARK_USED(ALL)
   CPABORT("This CP2K executable has not been linked against the required library libGint.")
#endif
   end subroutine update_fock_matrix_gpu




END MODULE libGint_wrapper
