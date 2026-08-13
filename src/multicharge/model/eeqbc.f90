! This file is part of multicharge.
! SPDX-Identifier: Apache-2.0
!
! Licensed under the Apache License, Version 2.0 (the "License");
! you may not use this file except in compliance with the License.
! You may obtain a copy of the License at
!
!     http://www.apache.org/licenses/LICENSE-2.0
!
! Unless required by applicable law or agreed to in writing, software
! distributed under the License is distributed on an "AS IS" BASIS,
! WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
! See the License for the specific language governing permissions and
! limitations under the License.

!> @file multicharge/model/eeqbc.f90
!> Provides implementation of the bond capacitor electronegativity equilibration model (EEQ_BC)

!> Bond capacitor electronegativity equilibration charge model published in
!>
!> Thomas Froitzheim, Marcel Müller, Andreas Hansen, and Stefan Grimme,
!> *J. Chem. Phys.*, **2025**, 162, 214109.
!> DOI: [10.1063/5.0268978](https://dx.doi.org/10.1063/5.0268978)
module multicharge_model_eeqbc
   use mctc_env, only: timer_type, format_time, error_type, wp
   use mctc_io, only: structure_type
   use mctc_io_constants, only: pi
   use mctc_ncoord, only: new_ncoord, cn_count, ncoord_type
   use mctc_csrlist, only: csr_list
   use multicharge_wignerseitz, only: new_wignerseitz_cell, wignerseitz_cell_type
   use multicharge_blascomp, only: gemv_cmp, gemm_cmp, gemm_cmp_212
   use multicharge_model_type, only: mchrg_model_type, get_dir_trans
   use multicharge_blas, only: gemv, gemm, symv
   use multicharge_model_cache, only: mchrg_cache
   implicit none
   private

   public :: eeqbc_model, new_eeqbc_model

   !> EEQBC model type, extends base mchrg_model_type.
   type, extends(mchrg_model_type) :: eeqbc_model
      !> Bond capacitance parameters
      real(wp), allocatable :: cap(:)
      !> Average coordination number
      real(wp), allocatable :: avg_cn(:)
      !> Exponent of error function in bond capacitance
      real(wp) :: kbc
      !> Exponent of the distance/CN normalization
      real(wp) :: norm_exp
      !> Van der Waals radii matrix (nat × nat)
      real(wp), allocatable :: rvdw(:, :)
   contains
      !> Update and allocate cache
      procedure :: update
      !> Calculate capacitance matrix
      procedure :: get_capacitance_matrix
      !> Calculate Coulomb matrix
      procedure :: get_coulomb_matrix
      !> Calculate derivatives of Coulomb matrix multiplied by charge
      procedure :: get_coulomb_derivs
      !> Calculate right-hand side (electronegativity vector)
      procedure :: get_xvec
      !> Calculate derivatives of EN vector
      procedure :: get_xvec_derivs
      !> Calculate constraint matrix (molecular)
      procedure :: get_cmat_0d
      !> Calculate constraint matrix (molecular) using neighbour list
      procedure :: get_cmat_0d_list
      !> Calculate full constraint matrix (periodic)
      procedure :: get_cmat_3d
      !> Calculate constraint matrix derivatives (molecular)
      procedure :: get_dcmat_0d
      !> Calculate constraint matrix derivatives (molecular) using neighbour list
      procedure :: get_dcmat_0d_list
      !> Calculate constraint matrix derivatives (periodic)
      procedure :: get_dcmat_3d
      !> Calculate gradient
      procedure :: get_grad
   end type eeqbc_model

   real(wp), parameter :: sqrtpi = sqrt(pi)
   real(wp), parameter :: sqrt2pi = sqrt(2.0_wp / pi)
   real(wp), parameter :: eps = sqrt(epsilon(0.0_wp))

   !> Default exponent of distance/CN normalization
   real(wp), parameter :: default_norm_exp = 1.0_wp

   !> Default exponent of error function in bond capacitance
   real(wp), parameter :: default_kbc = 0.65_wp

   !> Default cutoff radius
   real(wp), parameter :: cutoff = 25.0_wp

contains

!> Constructor for the EEQBC model.
   subroutine new_eeqbc_model(self, mol, error, chi, rad, &
   & eta, kcnchi, kqchi, kqeta, kcnrad, cap, avg_cn, rvdw, &
   & kbc, cutoff, cn_exp, rcov, en, cn_max, norm_exp)
      !> Bond capacitor electronegativity equilibration model
      type(eeqbc_model), intent(out) :: self
      !> Molecular structure data
      type(structure_type), intent(in) :: mol
      !> Error handling
      type(error_type), allocatable, intent(out) :: error
      !> Electronegativity
      real(wp), intent(in) :: chi(:)
      !> Exponent gaussian charge
      real(wp), intent(in) :: rad(:)
      !> Chemical hardness
      real(wp), intent(in) :: eta(:)
      !> CN scaling factor for electronegativity
      real(wp), intent(in) :: kcnchi(:)
      !> Local charge scaling factor for electronegativity
      real(wp), intent(in) :: kqchi(:)
      !> Local charge scaling factor for chemical hardness
      real(wp), intent(in) :: kqeta(:)
      !> CN scaling factor for charge width
      real(wp), intent(in) :: kcnrad
      !> Bond capacitance
      real(wp), intent(in) :: cap(:)
      !> Average coordination number
      real(wp), intent(in) :: avg_cn(:)
      !> Van-der-Waals radii
      real(wp), intent(in) :: rvdw(:, :)
      !> Exponent of error function in bond capacitance
      real(wp), intent(in), optional :: kbc
      !> Exponent of the distance normalization
      real(wp), intent(in), optional :: norm_exp
      !> Cutoff radius for coordination number
      real(wp), intent(in), optional :: cutoff
      !> Steepness of the CN counting function
      real(wp), intent(in), optional :: cn_exp
      !> Covalent radii for CN
      real(wp), intent(in), optional :: rcov(:)
      !> Maximum CN cutoff for CN
      real(wp), intent(in), optional :: cn_max
      !> Pauling electronegativities normalized to fluorine
      real(wp), intent(in), optional :: en(:)

      self%chi = chi
      self%rad = rad
      self%eta = eta
      self%kcnchi = kcnchi
      self%kqchi = kqchi
      self%kqeta = kqeta
      self%kcnrad = kcnrad
      self%cap = cap
      self%avg_cn = avg_cn
      self%rvdw = rvdw

      if (present(kbc)) then
         self%kbc = kbc
      else
         self%kbc = default_kbc
      end if

      if (present(norm_exp)) then
         self%norm_exp = norm_exp
      else
         self%norm_exp = default_norm_exp
      end if

      ! Coordination number
      call new_ncoord(self%ncoord, mol, cn_count%erf, error, &
      & cutoff=cutoff, kcn=cn_exp, rcov=rcov, cut=cn_max, &
      & norm_exp=self%norm_exp)
      ! Electronegativity weighted coordination number for local charge
      call new_ncoord(self%ncoord_en, mol, cn_count%erf_en, error, &
      & cutoff=cutoff, kcn=cn_exp, rcov=rcov, en=en, cut=cn_max, &
      & norm_exp=self%norm_exp)

   end subroutine new_eeqbc_model

!> Update coordination numbers and local charges, and set up Wigner–Seitz cell if periodic.
   subroutine update(self, mol, cache, trans, grad, list)
      !> EEQBC model type
      class(eeqbc_model), intent(in) :: self
      !> Structure type
      type(structure_type), intent(in) :: mol
      !> Multicharge cache
      type(mchrg_cache), intent(inout) :: cache
      !> Lattice vectors
      real(wp), intent(in) :: trans(:, :)
      !> Flag to compute derivatives
      logical, intent(in) :: grad
      !> Multicharge neighbourlist type
      type(csr_list), intent(in), optional :: list

      cache%trans = trans

      ! Refer CN and local charge arrays in cache
      if (.not. allocated(cache%cn)) then
         allocate(cache%cn(mol%nat))
      end if
      if (.not. allocated(cache%qloc)) then
         allocate(cache%qloc(mol%nat))
      end if

      if (grad) then
         cache%grad = .true.
      else
         cache%grad = .false.
      end if

      if (present(list)) then
         call self%ncoord%get_coordination_number(mol, trans, cache%cn, list=list)
         call self%local_charge(mol, trans, cache%qloc, list=list)
      else
         call self%ncoord%get_coordination_number(mol, trans, cache%cn)
         call self%local_charge(mol, trans, qloc=cache%qloc)
      end if

      if (any(mol%periodic) .and. .not. present(list)) then
         ! Create WSC
         call new_wignerseitz_cell(cache%wsc, mol)
      end if

   end subroutine update

!> Compute the capacitance matrix (and its derivatives if needed).
   subroutine get_capacitance_matrix(self, mol, ndim, cache, list)
      !> EEQBC model type
      class(eeqbc_model), intent(in) :: self
      !> Structure type
      type(structure_type), intent(in) :: mol
      !> System size
      integer, intent(in) :: ndim
      !> Multicharge cache
      type(mchrg_cache), intent(inout) :: cache
      !> Multicharge neighbourlist type
      type(csr_list), intent(in), optional :: list

      if (cache%grad) then
         if (present(list)) then
            if (.not. allocated(cache%dcdrdiag)) allocate(cache%dcdrdiag(3, mol%nat))
            if (.not. allocated(cache%dcdL)) allocate(cache%dcdL(3, 3, mol%nat))
         else
            if (.not. allocated(cache%dcdr)) allocate(cache%dcdr(3, mol%nat, ndim))
            if (.not. allocated(cache%dcdL)) allocate(cache%dcdL(3, 3, ndim))
         end if
      end if

      if (present(list)) then
         ! Allocate cmat
         if (.not. allocated(cache%clist)) then
            allocate(cache%clist(size(list%nlat)))
         end if
         if (.not. allocated(cache%cdiag)) then
            allocate(cache%cdiag(mol%nat))
         end if
         ! Neighbour list routines
         if (any(mol%periodic)) then
            call get_cmat_3d_list(self, mol, list, cache%clist, cache%cdiag)
            ! cmat gradients
            if (cache%grad) then
               call get_dcmat_3d_list(self, mol, list, cache)
            end if
         else
            call get_cmat_0d_list(self, mol, list, cache%clist, cache%cdiag)
            ! cmat gradients
            if (cache%grad) then
               call get_dcmat_0d_list(self, mol, list, cache)
            end if
         end if
      else
         ! Allocate cmat
         if (.not. allocated(cache%cmat)) then
            allocate(cache%cmat(ndim, ndim))
         end if
         ! Direct routines
         if (any(mol%periodic)) then
            ! Get full cmat sum over all WSC images (for get_xvec and xvec_derivs)
            call get_cmat_3d(self, mol, cache%wsc, cache%cmat)
            if (cache%grad) then
               call get_dcmat_3d(self, mol, cache%wsc, cache%dcdr, cache%dcdL)
            end if
         else
            call get_cmat_0d(self, mol, cache%cmat)
            ! cmat gradients
            if (cache%grad) then
               call get_dcmat_0d(self, mol, cache%dcdr, cache%dcdL)
            end if
         end if
      end if

   end subroutine get_capacitance_matrix

   subroutine get_xvec(self, mol, ndim, cache, list)
      !> EEQBC model type
      class(eeqbc_model), intent(in) :: self
      !> Structure type
      type(structure_type), intent(in) :: mol
      !> System size (number of atoms or atoms+1 if Lagrange multiplier used)
      integer, intent(in) :: ndim
      !> Multicharge cache
      type(mchrg_cache), intent(inout) :: cache
      !> Multicharge neighbourlist type
      type(csr_list), intent(in), optional :: list

      integer :: iat, izp, img, idx
      real(wp) :: ctmp, vec(3), rvdw, capi, wsw
      real(wp), allocatable :: dtrans(:, :)

      if (.not. allocated(cache%xtmp)) then
         allocate(cache%xtmp(ndim))
      end if

      if (.not. allocated(cache%xvec)) then
         allocate(cache%xvec(ndim))
      end if

      cache%xvec(:) = 0.0_wp

      ! Step 1: Perfectly parallelized
      !$omp parallel do default(none) schedule(runtime) &
      !$omp shared(mol, self, cache) &
      !$omp private(iat, izp)
      do iat = 1, mol%nat
         izp = mol%id(iat)
         cache%xtmp(iat) = -self%chi(izp) + self%kcnchi(izp) * cache%cn(iat) &
         & + self%kqchi(izp) * cache%qloc(iat)
      end do
      !$omp end parallel do

      ! Only write the extra element if xtmp has room for it
      if (size(cache%xtmp) == mol%nat + 1) then
         cache%xtmp(mol%nat + 1) = mol%charge
      end if

      ! Step 2: Matrix-vector multiplication
      if (present(list)) then
         call gemv_cmp(list, cache%clist, cache%cdiag, cache%xtmp, cache%xvec, alpha=1.0_wp, beta=0.0_wp)
      else
         call gemv(cache%cmat, cache%xtmp, cache%xvec)
      end if

      ! Step 3: Periodic corrections
      if (any(mol%periodic)) then
         call get_dir_trans(mol, dtrans, cutoff)

         if (present(list)) then
            !$omp parallel do default(none) schedule(runtime) &
            !$omp shared(mol, self, list, cache, dtrans) &
            !$omp private(iat, izp, img, wsw, capi, vec, rvdw, ctmp)
            do iat = 1, mol%nat
               izp = mol%id(iat)
               capi = self%cap(izp)
               rvdw = self%rvdw(izp, izp)

               ! Execute only if periodic self-images exist (nimg > 1, since index 1 is identity R=0)
               if (list%wsc%nimg_list(list%inl(iat)) > 1) then
                  wsw = 1.0_wp / real(list%wsc%nimg_list(list%inl(iat)) - 1, wp)

                  ! Skip identity image at list%wsc%itr_list(kat) and loop through R /= 0 images
                  do img = list%wsc%itr_list(list%inl(iat)) + 1, list%wsc%itr_list(list%inl(iat) + 1) - 1
                     vec = list%wsc%trans(:, list%wsc%tridx_list(img))
                     call get_cpair_dir(self%kbc, vec, dtrans, rvdw, capi, capi, ctmp)

                     ! Direct write is safe: 'iat' is private to this specific thread
                     cache%xvec(iat) = cache%xvec(iat) - wsw * ctmp * cache%xtmp(iat)
                  end do
               end if
            end do
            !$omp end parallel do


         else
            !$omp parallel do default(none) schedule(runtime) &
            !$omp shared(mol, self, cache, dtrans) &
            !$omp private(iat, izp, img, wsw, capi, vec, rvdw, ctmp)
            do iat = 1, mol%nat
               izp = mol%id(iat)
               capi = self%cap(izp)
               rvdw = self%rvdw(izp, izp)

               wsw = 1.0_wp / real(cache%wsc%nimg(iat, iat), wp)
               do img = 1, cache%wsc%nimg(iat, iat)
                  vec = cache%wsc%trans(:, cache%wsc%tridx(img, iat, iat))
                  call get_cpair_dir(self%kbc, vec, dtrans, rvdw, capi, capi, ctmp)

                  ! Direct write is safe: 'iat' is private to this specific thread
                  cache%xvec(iat) = cache%xvec(iat) - wsw * ctmp * cache%xtmp(iat)
               end do
            end do
            !$omp end parallel do
         end if
      end if

   end subroutine get_xvec

!> Compute derivatives of the electronegativity vector with respect to atomic positions and lattice parameters.
   subroutine get_xvec_derivs(self, mol, ndim, cache, list)
      !> EEQBC model type
      class(eeqbc_model), intent(in) :: self
      !> Structure type
      type(structure_type), intent(in) :: mol
      !> System size
      integer, intent(in) :: ndim
      !> Multicharge cache
      type(mchrg_cache), intent(inout) :: cache
      !> Multicharge neighbourlist type
      type(csr_list), intent(in), optional :: list

      if (.not. allocated(cache%dcndr)) then
         allocate(cache%dcndr(3, mol%nat, mol%nat))
      end if
      if (.not. allocated(cache%dcndL)) then
         allocate(cache%dcndL(3, 3, mol%nat))
         call self%ncoord%get_coordination_number(mol, cache%trans, cache%cn, dcndr=cache%dcndr, dcndL=cache%dcndL)
      end if
      if (.not. allocated(cache%dqlocdr)) then
         allocate(cache%dqlocdr(3, mol%nat, mol%nat))
      end if
      if (.not. allocated(cache%dqlocdL)) then
         allocate(cache%dqlocdL(3, 3, mol%nat))
         call self%local_charge(mol, cache%trans, qloc=cache%qloc, dqlocdr=cache%dqlocdr, dqlocdL=cache%dqlocdL)
      end if

      if (any(mol%periodic)) then
         if (present(list)) then
            if (.not. allocated(cache%dxdrij)) allocate(cache%dxdrij(3, size(list%nlat)))
            if (.not. allocated(cache%dxdrji)) allocate(cache%dxdrji(3, size(list%nlat)))
            if (.not. allocated(cache%dxdrdiag)) allocate(cache%dxdrdiag(3, mol%nat))
            if (.not. allocated(cache%dxdL)) allocate(cache%dxdL(3, 3, mol%nat))
            call get_xvec_derivs_3d_list(self, mol, cache, list)
         else
            if (.not. allocated(cache%dxdr)) allocate(cache%dxdr(3, mol%nat, ndim))
            if (.not. allocated(cache%dxdL)) allocate(cache%dxdL(3, 3, ndim))
            call get_xvec_derivs_3d(self, mol, ndim, cache)
         end if
      else
         if (present(list)) then
            if (.not. allocated(cache%dxdrij)) allocate(cache%dxdrij(3, size(list%nlat)))
            if (.not. allocated(cache%dxdrji)) allocate(cache%dxdrji(3, size(list%nlat)))
            if (.not. allocated(cache%dxdrdiag)) allocate(cache%dxdrdiag(3, mol%nat))
            if (.not. allocated(cache%dxdL)) allocate(cache%dxdL(3, 3, mol%nat))
            call get_xvec_derivs_0d_list(self, mol, cache, list)
         else
            if (.not. allocated(cache%dxdr)) allocate(cache%dxdr(3, mol%nat, ndim))
            if (.not. allocated(cache%dxdL)) allocate(cache%dxdL(3, 3, ndim))
            call get_xvec_derivs_0d(self, mol, ndim, cache)
         end if
      end if

   end subroutine get_xvec_derivs

!> Compute derivatives of the electronegativity vector for non-periodic system.
   subroutine get_xvec_derivs_0d(self, mol, ndim, cache)
      !> EEQBC model type
      class(eeqbc_model), intent(in) :: self
      !> Structure type
      type(structure_type), intent(in) :: mol
      !> System size
      integer, intent(in) :: ndim
      !> Multicharge cache
      type(mchrg_cache), intent(inout) :: cache

      integer :: iat, izp, jat, jzp
      real(wp) :: capi, capj, vec(3), ctmp, rvdw, dG(3), dS(3, 3)
      real(wp), allocatable :: dtmpdr(:, :, :), dtmpdL(:, :, :)

      ! Thread-private arrays for reduction
      real(wp), allocatable :: dxdr_local(:, :, :), dxdL_local(:, :, :)
      real(wp), allocatable :: dtmpdr_local(:, :, :), dtmpdL_local(:, :, :)

      allocate(dtmpdr(3, mol%nat, ndim), dtmpdL(3, 3, ndim))

      cache%dxdr(:, :, :) = 0.0_wp
      cache%dxdL(:, :, :) = 0.0_wp
      dtmpdr(:, :, :) = 0.0_wp
      dtmpdL(:, :, :) = 0.0_wp

      !$omp parallel default(none) &
      !$omp shared(mol, self, cache, dtmpdr, dtmpdL) &
      !$omp private(iat, izp, dtmpdr_local, dtmpdL_local)
      allocate(dtmpdr_local, source=dtmpdr)
      allocate(dtmpdL_local, source=dtmpdL)
      !$omp do schedule(runtime)
      do iat = 1, mol%nat
         izp = mol%id(iat)
         ! CN and effective charge derivative
         dtmpdr_local(:, :, iat) = self%kcnchi(izp) * cache%dcndr(:, :, iat) + dtmpdr_local(:, :, iat)
         dtmpdL_local(:, :, iat) = self%kcnchi(izp) * cache%dcndL(:, :, iat) + dtmpdL_local(:, :, iat)
         dtmpdr_local(:, :, iat) = self%kqchi(izp) * cache%dqlocdr(:, :, iat) + dtmpdr_local(:, :, iat)
         dtmpdL_local(:, :, iat) = self%kqchi(izp) * cache%dqlocdL(:, :, iat) + dtmpdL_local(:, :, iat)
      end do
      !$omp end do
      !$omp critical (get_xvec_derivs_0d_)
      dtmpdr(:, :, :) = dtmpdr + dtmpdr_local
      dtmpdL(:, :, :) = dtmpdL + dtmpdL_local
      !$omp end critical (get_xvec_derivs_0d_)
      deallocate(dtmpdL_local, dtmpdr_local)
      !$omp end parallel

      call gemm(dtmpdr, cache%cmat, cache%dxdr)
      call gemm(dtmpdL, cache%cmat, cache%dxdL)

      !$omp parallel default(none) &
      !$omp shared(mol, self, cache) &
      !$omp private(iat, izp, jat, jzp, vec, dxdr_local, dxdL_local)
      allocate(dxdr_local, mold=cache%dxdr)
      allocate(dxdL_local, mold=cache%dxdL)
      dxdr_local(:, :, :) = 0.0_wp
      dxdL_local(:, :, :) = 0.0_wp
      !$omp do schedule(runtime)
      do iat = 1, mol%nat
         do jat = 1, iat - 1
            ! Diagonal elements
            dxdr_local(:, iat, iat) = dxdr_local(:, iat, iat) + cache%xtmp(jat) * cache%dcdr(:, iat, jat)
            dxdr_local(:, jat, jat) = dxdr_local(:, jat, jat) + cache%xtmp(iat) * cache%dcdr(:, jat, iat)

            ! Derivative of capacitance matrix
            dxdr_local(:, iat, jat) = (cache%xtmp(iat) - cache%xtmp(jat)) * cache%dcdr(:, iat, jat) + dxdr_local(:, iat, jat)
            dxdr_local(:, jat, iat) = (cache%xtmp(jat) - cache%xtmp(iat)) * cache%dcdr(:, jat, iat) + dxdr_local(:, jat, iat)

            vec = mol%xyz(:, iat) - mol%xyz(:, jat)
            dxdL_local(:, :, iat) = dxdL_local(:, :, iat) + cache%xtmp(jat) * &
            & spread(cache%dcdr(:, iat, jat), 1, 3) * spread(vec, 2, 3)
            dxdL_local(:, :, jat) = dxdL_local(:, :, jat) + cache%xtmp(iat) * &
            & spread(cache%dcdr(:, jat, iat), 1, 3) * spread(-vec, 2, 3)
         end do
         dxdr_local(:, iat, iat) = dxdr_local(:, iat, iat) + cache%xtmp(iat) * cache%dcdr(:, iat, iat)
         dxdL_local(:, :, iat) = dxdL_local(:, :, iat) + cache%xtmp(iat) * cache%dcdL(:, :, iat)
      end do
      !$omp end do
      !$omp critical (get_xvec_derivs_0d_)
      cache%dxdr(:, :, :) = cache%dxdr + dxdr_local
      cache%dxdL(:, :, :) = cache%dxdL + dxdL_local
      !$omp end critical (get_xvec_derivs_0d_)
      deallocate(dxdL_local, dxdr_local)
      !$omp end parallel

   end subroutine get_xvec_derivs_0d

   !> Compute derivatives of the electronegativity vector for non-periodic system using neighbour list.
   subroutine get_xvec_derivs_0d_list(self, mol, cache, list)
      class(eeqbc_model), intent(in) :: self
      type(structure_type), intent(in) :: mol
      type(mchrg_cache), intent(inout) :: cache
      type(csr_list), intent(in) :: list

      integer :: iat, izp, jat, kat, jzp
      real(wp) :: vec(3)
      real(wp), allocatable :: dtmpdrij(:, :), dtmpdrji(:, :)
      real(wp), allocatable :: dtmpdrdiag(:, :), dtmpdL(:, :, :)

      ! Thread-private arrays for reduction
      real(wp), allocatable :: dxdrij_local(:, :), dxdrji_local(:, :), dxdrdiag_local(:, :)
      real(wp), allocatable :: dxdL_local(:, :, :)
      real(wp), allocatable :: dtmpdrij_local(:, :), dtmpdrji_local(:, :), dtmpdrdiag_local(:, :)
      real(wp), allocatable :: dtmpdL_local(:, :, :)

      allocate(dtmpdrij(3, size(list%nlat)), dtmpdrji(3, size(list%nlat)),&
      & dtmpdrdiag(3, mol%nat), dtmpdL(3, 3, mol%nat))

      cache%dxdrij(:, :) = 0.0_wp
      cache%dxdrji(:, :) = 0.0_wp
      cache%dxdrdiag(:, :) = 0.0_wp
      cache%dxdL(:, :, :) = 0.0_wp
      dtmpdrij(:, :) = 0.0_wp
      dtmpdrji(:, :) = 0.0_wp
      dtmpdrdiag(:, :) = 0.0_wp
      dtmpdL(:, :, :) = 0.0_wp

      !$omp parallel default(none) &
      !$omp shared(mol, self, list, cache, dtmpdrij, dtmpdrji, dtmpdrdiag, dtmpdL) &
      !$omp private(iat, jat, jzp, izp, kat)

      !$omp do schedule(runtime)
      do iat = 1, mol%nat
         izp = mol%id(iat)
         ! Standard diagonal contributions
         dtmpdrdiag(:, iat) = self%kcnchi(izp) * cache%dcndrdiag(:, iat) + dtmpdrdiag(:, iat)
         dtmpdL(:, :, iat) = self%kcnchi(izp) * cache%dcndL(:, :, iat) + dtmpdL(:, :, iat)
         dtmpdrdiag(:, iat) = self%kqchi(izp) * cache%dqlocdrdiag(:, iat) + dtmpdrdiag(:, iat)
         dtmpdL(:, :, iat) = self%kqchi(izp) * cache%dqlocdL(:, :, iat) + dtmpdL(:, :, iat)

         do kat = list%inl(iat) + 1, list%inl(iat+1) - 1
            jat = list%nlat(kat)
            jzp = mol%id(jat)
            ! Since jat > iat, kat represents the pair (iat, jat)
            dtmpdrij(:, kat) = self%kcnchi(jzp) * cache%dcndrij(:, kat) + dtmpdrij(:, kat)
            dtmpdrji(:, kat) = self%kcnchi(izp) * cache%dcndrji(:, kat) + dtmpdrji(:, kat)
            dtmpdrij(:, kat) = self%kqchi(jzp) * cache%dqlocdrij(:, kat) + dtmpdrij(:, kat)
            dtmpdrji(:, kat) = self%kqchi(izp) * cache%dqlocdrji(:, kat) + dtmpdrji(:, kat)
         end do
      end do
      !$omp end do
      !$omp end parallel


      call gemm_cmp_212(list, cache%clist, cache%cdiag, dtmpdrij, dtmpdrji, dtmpdrdiag, &
      & cache%dxdrij, cache%dxdrji, cache%dxdrdiag, 1.0_wp, 0.0_wp)
      call gemm_cmp(list, cache%clist, cache%cdiag, dtmpdL, cache%dxdL, 1.0_wp, 0.0_wp)

      !$omp parallel default(none) &
      !$omp shared(mol, self, cache, list) &
      !$omp private(iat, jat, kat, vec)
      !$omp do schedule(runtime)
      do iat = 1, mol%nat
         do kat = list%inl(iat) + 1, list%inl(iat+1) - 1
            jat = list%nlat(kat)

            ! Diagonal updates (Correct for symmetric list jat > iat)
            cache%dxdrdiag(:, iat) = cache%dxdrdiag(:, iat) + cache%xtmp(jat) * cache%dcdrij(:, kat)
            cache%dxdrdiag(:, jat) = cache%dxdrdiag(:, jat) + cache%xtmp(iat) * cache%dcdrji(:, kat)

            ! Off-diagonal updates
            cache%dxdrij(:, kat) = (cache%xtmp(iat) - cache%xtmp(jat)) * cache%dcdrij(:, kat) + cache%dxdrij(:, kat)
            cache%dxdrji(:, kat) = (cache%xtmp(jat) - cache%xtmp(iat)) * cache%dcdrji(:, kat) + cache%dxdrji(:, kat)

            ! Cell parameter updates
            vec = mol%xyz(:, iat) - mol%xyz(:, jat)
            cache%dxdL(:, :, iat) = cache%dxdL(:, :, iat) + cache%xtmp(jat) * &
            & spread(cache%dcdrij(:, kat), 1, 3) * spread(vec, 2, 3)
            cache%dxdL(:, :, jat) = cache%dxdL(:, :, jat) + cache%xtmp(iat) * &
            & spread(cache%dcdrji(:, kat), 1, 3) * spread(-vec, 2, 3)
         end do
         ! Add pure self-diagonal contribution
         cache%dxdrdiag(:, iat) = cache%dxdrdiag(:, iat) + cache%xtmp(iat) * cache%dcdrdiag(:, iat)
         cache%dxdL(:, :, iat) = cache%dxdL(:, :, iat) + cache%xtmp(iat) * cache%dcdL(:, :, iat)
      end do
      !$omp end do
      !$omp end parallel


   end subroutine get_xvec_derivs_0d_list

!> Compute derivatives of the electronegativity vector for peridoic system.
   subroutine get_xvec_derivs_3d(self, mol, ndim, cache)
      !> EEQBC model type
      class(eeqbc_model), intent(in) :: self
      !> Structure type
      type(structure_type), intent(in) :: mol
      !> System size
      integer, intent(in) :: ndim
      !> Multicharge cache
      type(mchrg_cache), intent(inout) :: cache

      integer :: iat, izp, jat, jzp, img
      real(wp) :: capi, capj, wsw, rvdw, vec(3), ctmp, dG(3), dS(3, 3)
      real(wp), allocatable :: dtmpdr(:, :, :), dtmpdL(:, :, :)
      real(wp), allocatable :: dtrans(:, :)

      ! Thread-private arrays for reduction
      real(wp), allocatable :: dxdr_local(:, :, :), dxdL_local(:, :, :)
      real(wp), allocatable :: dtmpdr_local(:, :, :), dtmpdL_local(:, :, :)

      allocate(dtmpdr(3, mol%nat, ndim), dtmpdL(3, 3, ndim))

      cache%dxdr(:, :, :) = 0.0_wp
      cache%dxdL(:, :, :) = 0.0_wp
      dtmpdr(:, :, :) = 0.0_wp
      dtmpdL(:, :, :) = 0.0_wp

      !$omp parallel default(none) &
      !$omp shared(mol, self, cache, dtmpdr, dtmpdL) &
      !$omp private(iat, izp, dtmpdr_local, dtmpdL_local)
      allocate(dtmpdr_local, source=dtmpdr)
      allocate(dtmpdL_local, source=dtmpdL)
      !$omp do schedule(runtime)
      do iat = 1, mol%nat
         izp = mol%id(iat)
         ! CN and effective charge derivative
         dtmpdr_local(:, :, iat) = self%kcnchi(izp) * cache%dcndr(:, :, iat) + dtmpdr_local(:, :, iat)
         dtmpdL_local(:, :, iat) = self%kcnchi(izp) * cache%dcndL(:, :, iat) + dtmpdL_local(:, :, iat)
         dtmpdr_local(:, :, iat) = self%kqchi(izp) * cache%dqlocdr(:, :, iat) + dtmpdr_local(:, :, iat)
         dtmpdL_local(:, :, iat) = self%kqchi(izp) * cache%dqlocdL(:, :, iat) + dtmpdL_local(:, :, iat)
      end do
      !$omp end do
      !$omp critical (get_xvec_derivs_3d_)
      dtmpdr(:, :, :) = dtmpdr + dtmpdr_local
      dtmpdL(:, :, :) = dtmpdL + dtmpdL_local
      !$omp end critical (get_xvec_derivs_3d_)
      deallocate(dtmpdL_local, dtmpdr_local)
      !$omp end parallel

      call gemm(dtmpdr, cache%cmat, cache%dxdr)
      call gemm(dtmpdL, cache%cmat, cache%dxdL)

      call get_dir_trans(mol, dtrans, cutoff)
      !$omp parallel default(none) &
      !$omp shared(mol, self, cache, dtrans) &
      !$omp private(iat, izp, jat, jzp, img, wsw) &
      !$omp private(capi, capj, vec, rvdw, ctmp, dG, dS) &
      !$omp private(dxdr_local, dxdL_local)
      allocate(dxdr_local, mold=cache%dxdr)
      allocate(dxdL_local, mold=cache%dxdL)
      dxdr_local(:, :, :) = 0.0_wp
      dxdL_local(:, :, :) = 0.0_wp
      !$omp do schedule(runtime)
      do iat = 1, mol%nat
         izp = mol%id(iat)
         capi = self%cap(izp)
         do jat = 1, mol%nat
            jzp = mol%id(jat)
            rvdw = self%rvdw(izp, jzp)
            capj = self%cap(jzp)

            ! Diagonal elements
            dxdr_local(:, iat, iat) = dxdr_local(:, iat, iat) + cache%xtmp(jat) * cache%dcdr(:, iat, jat)

            ! Derivative of capacitance matrix
            dxdr_local(:, iat, jat) = dxdr_local(:, iat, jat)  &
            & + (cache%xtmp(iat) - cache%xtmp(jat)) * cache%dcdr(:, iat, jat)

            wsw = 1.0_wp / real(cache%wsc%nimg(iat, jat), wp)
            do img = 1, cache%wsc%nimg(iat, jat)
               vec = mol%xyz(:, jat) - mol%xyz(:, iat) + cache%wsc%trans(:, cache%wsc%tridx(img, jat, iat))
               call get_dcpair_dir(self%kbc, vec, dtrans, rvdw, capi, capj, dG, dS)
               dxdL_local(:, :, iat) = dxdL_local(:, :, iat) - wsw * dS * cache%xtmp(jat)
            end do
         end do
         dxdL_local(:, :, iat) = dxdL_local(:, :, iat) + cache%xtmp(iat) * cache%dcdL(:, :, iat)

         ! Capacitance terms for i = j, T != 0
         rvdw = self%rvdw(izp, izp)
         wsw = 1.0_wp / real(cache%wsc%nimg(iat, iat), wp)
         do img = 1, cache%wsc%nimg(iat, iat)
            vec = cache%wsc%trans(:, cache%wsc%tridx(img, iat, iat))

            call get_cpair_dir(self%kbc, vec, dtrans, rvdw, capi, capi, ctmp)
            ctmp = ctmp * wsw
            ! EN derivative
            dxdr_local(:, :, iat) = dxdr_local(:, :, iat) - ctmp * self%kcnchi(izp) * cache%dcndr(:, :, iat)
            dxdL_local(:, :, iat) = dxdL_local(:, :, iat) - ctmp * self%kcnchi(izp) * cache%dcndL(:, :, iat)
            dxdr_local(:, :, iat) = dxdr_local(:, :, iat) - ctmp * self%kqchi(izp) * cache%dqlocdr(:, :, iat)
            dxdL_local(:, :, iat) = dxdL_local(:, :, iat) - ctmp * self%kqchi(izp) * cache%dqlocdL(:, :, iat)
         end do
      end do
      !$omp end do
      !$omp critical (get_xvec_derivs_3d_)
      cache%dxdr(:, :, :) = cache%dxdr + dxdr_local
      cache%dxdL(:, :, :) = cache%dxdL + dxdL_local
      !$omp end critical (get_xvec_derivs_3d_)
      deallocate(dxdL_local, dxdr_local)
      !$omp end parallel

   end subroutine get_xvec_derivs_3d

!> Compute derivatives of the electronegativity vector for peridoic system.
   subroutine get_xvec_derivs_3d_list(self, mol, cache, list)
      !> EEQBC model type
      class(eeqbc_model), intent(in) :: self
      !> Molecular structure data
      type(structure_type), intent(in) :: mol
      !> Multicharge cache (contains xtmp, dcdr, dcdL, etc.)
      type(mchrg_cache), intent(inout) :: cache
      !> Multicharge neighbourlist type
      type(csr_list), intent(in) :: list

      integer :: iat, izp, jat, kat, jzp, img
      real(wp) :: capi, capj, vec(3), ctmp, dG(3), dS(3, 3), rvdw, wsw
      real(wp), allocatable :: dtrans(:, :)
      real(wp), allocatable :: dtmpdrij(:, :), dtmpdrji(:, :)
      real(wp), allocatable :: dtmpdrdiag(:, :), dtmpdL(:, :, :)

      ! Thread-private arrays for reduction
      real(wp), allocatable :: dxdrij_local(:, :), dxdrji_local(:, :), dxdrdiag_local(:, :)
      real(wp), allocatable :: dxdL_local(:, :, :)
      real(wp), allocatable :: dtmpdrij_local(:, :), dtmpdrji_local(:, :), dtmpdrdiag_local(:, :)
      real(wp), allocatable :: dtmpdL_local(:, :, :)

      allocate(dtmpdrij(3, size(list%nlat)), dtmpdrji(3, size(list%nlat)),&
      & dtmpdrdiag(3, mol%nat), dtmpdL(3, 3, mol%nat))

      cache%dxdrij(:, :) = 0.0_wp
      cache%dxdrji(:, :) = 0.0_wp
      cache%dxdrdiag(:, :) = 0.0_wp
      cache%dxdL(:, :, :) = 0.0_wp
      dtmpdrij(:, :) = 0.0_wp
      dtmpdrji(:, :) = 0.0_wp
      dtmpdrdiag(:, :) = 0.0_wp
      dtmpdL(:, :, :) = 0.0_wp

      !$omp parallel default(none) &
      !$omp shared(mol, self, list, cache, dtmpdrij, dtmpdrji, dtmpdrdiag, dtmpdL) &
      !$omp private(iat, jat, jzp, izp, dtmpdrij_local, dtmpdrji_local, dtmpdrdiag_local, dtmpdL_local)
      allocate(dtmpdrdiag_local, source=dtmpdrdiag)
      allocate(dtmpdrij_local, source=dtmpdrij)
      allocate(dtmpdrji_local, source=dtmpdrji)
      allocate(dtmpdL_local, source=dtmpdL)
      !$omp do schedule(runtime)
      do iat = 1, mol%nat
         izp = mol%id(iat)
         ! CN and effective charge derivative
         dtmpdrdiag_local(:, iat) = self%kcnchi(izp) * cache%dcndrdiag(:, iat) + dtmpdrdiag_local(:, iat)
         dtmpdL_local(:, :, iat) = self%kcnchi(izp) * cache%dcndL(:, :, iat) + dtmpdL_local(:, :, iat)
         dtmpdrdiag_local(:, iat) = self%kqchi(izp) * cache%dqlocdrdiag(:, iat) + dtmpdrdiag_local(:, iat)
         dtmpdL_local(:, :, iat) = self%kqchi(izp) * cache%dqlocdL(:, :, iat) + dtmpdL_local(:, :, iat)
      end do
      !$omp end do
      !$omp do schedule(runtime)
      do kat = 1, size(list%nlat)
         jat = list%nlat(kat)
         jzp = mol%id(jat)
         dtmpdrij_local(:, kat) = self%kcnchi(jzp) * cache%dcndrij(:, kat) + dtmpdrij_local(:, kat)
         dtmpdrji_local(:, kat) = self%kcnchi(jzp) * cache%dcndrji(:, kat) + dtmpdrji_local(:, kat)
         dtmpdrij_local(:, kat) = self%kqchi(jzp) * cache%dqlocdrij(:, kat) + dtmpdrij_local(:, kat)
         dtmpdrji_local(:, kat) = self%kqchi(jzp) * cache%dqlocdrji(:, kat) + dtmpdrji_local(:, kat)
      end do
      !$omp end do
      !$omp critical (get_xvec_derivs_3d_list_)
      dtmpdrij(:, :) = dtmpdrij + dtmpdrij_local
      dtmpdrji(:, :) = dtmpdrji + dtmpdrji_local
      dtmpdrdiag(:, :) = dtmpdrdiag + dtmpdrdiag_local
      dtmpdL(:, :, :) = dtmpdL + dtmpdL_local
      !$omp end critical (get_xvec_derivs_3d_list_)
      deallocate(dtmpdL_local, dtmpdrij_local, dtmpdrji_local, dtmpdrdiag_local)
      !$omp end parallel

      call gemm_cmp(list, cache%clist, cache%cdiag, dtmpdrij, dtmpdrji, dtmpdrdiag, &
      & cache%dxdrij, cache%dxdrji, cache%dxdrdiag, 1.0_wp, 0.0_wp)
      call gemm_cmp(list, cache%clist, cache%cdiag, dtmpdL, cache%dxdL, 1.0_wp, 0.0_wp)

      call get_dir_trans(mol, dtrans, cutoff)

      !$omp parallel default(none) &
      !$omp shared(mol, self, cache, list, dtrans) &
      !$omp private(iat, izp, jat, kat, jzp, img, vec, dxdrdiag_local) &
      !$omp private(dxdrij_local, dxdrji_local, dxdL_local, capi, capj, rvdw, wsw, dG, dS, ctmp)
      allocate(dxdrdiag_local, mold=cache%dxdrdiag)
      allocate(dxdrij_local, mold=cache%dxdrij)
      allocate(dxdrji_local, mold=cache%dxdrji)
      allocate(dxdL_local, mold=cache%dxdL)
      dxdrdiag_local(:, :) = 0.0_wp
      dxdrij_local(:, :) = 0.0_wp
      dxdrji_local(:, :) = 0.0_wp
      dxdL_local(:, :, :) = 0.0_wp

      !$omp do schedule(runtime)
      do iat = 1, mol%nat
         izp = mol%id(iat)
         capi = self%cap(izp)

         ! Loop over all neighbours, including self‑images
         do kat = list%inl(iat) + 1, list%inl(iat+1) - 1
            jat = list%nlat(kat)
            jzp = mol%id(jat)
            capj = self%cap(jzp)
            rvdw = self%rvdw(izp, jzp)
            wsw = 1.0_wp / real(list%wsc%nimg_list(kat)-1, wp)

            if (jat /= iat) then
               ! Off‑diagonal pair (iat ≠ jat)
               ! Derivative of capacitance matrix (dcdr)
               dxdrdiag_local(:, iat) = dxdrdiag_local(:, iat) + cache%xtmp(jat) * cache%dcdrij(:, kat)
               dxdrdiag_local(:, jat) = dxdrdiag_local(:, jat) + cache%xtmp(iat) * cache%dcdrji(:, kat)
               dxdrij_local(:, kat) = dxdrij_local(:, kat) + (cache%xtmp(iat) - cache%xtmp(jat)) * cache%dcdrij(:, kat)
               dxdrji_local(:, kat) = dxdrji_local(:, kat) + (cache%xtmp(jat) - cache%xtmp(iat)) * cache%dcdrji(:, kat)

               ! Periodic images: lattice derivative contributions from capacitance matrix
               do img = list%wsc%itr_list(kat), list%wsc%itr_list(kat+1) - 1
                  vec = mol%xyz(:, jat) - mol%xyz(:, iat) + list%wsc%trans(:, list%wsc%tridx_list(img))
                  call get_dcpair_dir(self%kbc, vec, dtrans, rvdw, capi, capj, dG, dS)
                  dxdL_local(:, :, iat) = dxdL_local(:, :, iat) - wsw * dS * cache%xtmp(jat)
               end do
            else
               ! Execute only if periodic self-images exist (nimg > 1, since index 1 is identity R=0)
               if (list%wsc%itr_list(list%inl(iat+1)) - list%wsc%itr_list(list%inl(iat)) > 1) then
                  wsw = 1.0_wp / real(list%wsc%nimg_list(list%inl(iat)) - 1, wp)

                  ! Skip identity image at list%wsc%itr_list(list%inl(iat)) and loop through R /= 0 images
                  do img = list%wsc%itr_list(list%inl(iat)), list%wsc%itr_list(list%inl(iat+1)) - 1
                     vec = list%wsc%trans(:, list%wsc%tridx_list(img))
                     call get_cpair_dir(self%kbc, vec, dtrans, rvdw, capi, capi, ctmp)
                     ctmp = ctmp * wsw
                     ! EN derivative contributions
                     dxdrdiag_local(:, iat) = dxdrdiag_local(:, iat) - ctmp * self%kcnchi(izp) * cache%dcndrdiag(:, iat)
                     dxdL_local(:, :, iat) = dxdL_local(:, :, iat) - ctmp * self%kcnchi(izp) * cache%dcndL(:, :, iat)
                     ! Local charge derivative contributions
                     dxdrdiag_local(:, iat) = dxdrdiag_local(:, iat) - ctmp * self%kqchi(izp) * cache%dqlocdrdiag(:, iat)
                     dxdL_local(:, :, iat) = dxdL_local(:, :, iat) - ctmp * self%kqchi(izp) * cache%dqlocdL(:, :, iat)
                  end do
               end if
            end if
         end do

         ! Add the direct lattice derivative contribution from the diagonal of dcdL
         dxdL_local(:, :, iat) = dxdL_local(:, :, iat) + cache%xtmp(iat) * cache%dcdL(:, :, iat)
      end do
      !$omp end do

      !$omp critical (get_xvec_derivs_3d_list_)
      cache%dxdrdiag(:, :) = cache%dxdrdiag + dxdrdiag_local
      cache%dxdrij(:, :) = cache%dxdrij + dxdrij_local
      cache%dxdrji(:, :) = cache%dxdrji + dxdrji_local
      cache%dxdL(:, :, :) = cache%dxdL + dxdL_local
      !$omp end critical (get_xvec_derivs_3d_list_)
      deallocate(dxdL_local, dxdrij_local, dxdrji_local, dxdrdiag_local)
      !$omp end parallel

   end subroutine get_xvec_derivs_3d_list

!> Assemble the Coulomb matrix (periodic or non‑periodic) including bond capacitance contributions.
   subroutine get_coulomb_matrix(self, mol, ndim, cache, list)
      !> EEQBC model type
      class(eeqbc_model), intent(in) :: self
      !> Structure type
      type(structure_type), intent(in) :: mol
      !> System size
      integer, intent(in) :: ndim
      !> Multicharge cache
      type(mchrg_cache), intent(inout) :: cache
      !> Multicharge neighbourlist type
      type(csr_list), intent(in), optional :: list


      if (present(list)) then
         ! Allocate amat
         if (.not. allocated(cache%alist)) then
            allocate(cache%alist(size(list%nlat)))
         end if
         if (.not. allocated(cache%adiag)) then
            allocate(cache%adiag(mol%nat))
         end if
         if (any(mol%periodic)) then
            call get_amat_3d_list(self, mol, list, cache)
         else
            call get_amat_0d_list(self, mol, list, cache)
         end if
      else
         ! Allocate amat
         if (.not. allocated(cache%amat)) then
            allocate(cache%amat(ndim, ndim))
         end if
         if (any(mol%periodic)) then
            call get_amat_3d(self, mol, cache)
         else
            call get_amat_0d(self, mol, cache)
         end if
      end if
   end subroutine get_coulomb_matrix

!> Build the Coulomb matrix for a non‑periodic system (0D).
   subroutine get_amat_0d(self, mol, cache)
      !> EEQBC model type
      class(eeqbc_model), intent(in) :: self
      !> Molecular structure data
      type(structure_type), intent(in) :: mol
      !> Multicharge cache
      type(mchrg_cache), intent(inout) :: cache

      integer :: iat, jat, izp, jzp
      real(wp) :: vec(3), r2, gam2, tmp, norm_cn, radi, radj

      ! Thread-private array for reduction
      real(wp), allocatable :: amat_local(:, :)

      cache%amat(:, :) = 0.0_wp

      !$omp parallel default(none) &
      !$omp shared(cache, mol, self) &
      !$omp private(iat, izp, jat, jzp, gam2, vec, r2, tmp) &
      !$omp private(norm_cn, radi, radj, amat_local)
      allocate(amat_local, source=cache%amat)
      !$omp do schedule(runtime)
      do iat = 1, mol%nat
         izp = mol%id(iat)
         ! Effective charge width of i
         norm_cn = 1.0_wp / self%avg_cn(izp)**self%norm_exp
         radi = self%rad(izp) * (1.0_wp - self%kcnrad * cache%cn(iat) * norm_cn)
         do jat = 1, iat - 1
            jzp = mol%id(jat)
            vec = mol%xyz(:, jat) - mol%xyz(:, iat)
            r2 = vec(1)**2 + vec(2)**2 + vec(3)**2
            ! Effective charge width of j
            norm_cn = cache%cn(jat) / self%avg_cn(jzp)**self%norm_exp
            radj = self%rad(jzp) * (1.0_wp - self%kcnrad * norm_cn)
            ! Coulomb interaction of Gaussian charges
            gam2 = 1.0_wp / (radi**2 + radj**2)
            tmp = erf(sqrt(r2 * gam2)) / sqrt(r2) * cache%cmat(jat, iat)
            amat_local(jat, iat) = tmp
            amat_local(iat, jat) = tmp
         end do
         ! Effective hardness
         tmp = self%eta(izp) + self%kqeta(izp) * cache%qloc(iat) + sqrt2pi / radi
         amat_local(iat, iat) = amat_local(iat, iat) + tmp * cache%cmat(iat, iat) + 1.0_wp
      end do
      !$omp end do
      !$omp critical (get_amat_0d_)
      cache%amat(:, :) = cache%amat + amat_local
      !$omp end critical (get_amat_0d_)
      deallocate(amat_local)
      !$omp end parallel

      if (size(cache%amat, 1) == mol%nat + 1) then
         cache%amat(mol%nat + 1, 1:mol%nat + 1) = 1.0_wp
         cache%amat(1:mol%nat + 1, mol%nat + 1) = 1.0_wp
         cache%amat(mol%nat + 1, mol%nat + 1) = 0.0_wp
      end if

   end subroutine get_amat_0d

!> Build the Coulomb matrix for a non‑periodic system (0D).
   subroutine get_amat_0d_list(self, mol, list, cache)
      !> EEQBC model type
      class(eeqbc_model), intent(in) :: self
      !> Molecular structure data
      type(structure_type), intent(in) :: mol
      !> Multicharge neighbourlist type
      type(csr_list), intent(in) :: list
      !> Multicharge cache
      type(mchrg_cache), intent(inout) :: cache

      integer :: iat, jat, kat, izp, jzp
      real(wp) :: vec(3), r2, gam2, tmp, norm_cn, radi, radj

      ! Zero out global shared target arrays upfront
      cache%alist(:) = 0.0_wp
      cache%adiag(:) = 0.0_wp

      !$omp parallel default(none) &
      !$omp shared(cache, mol, self, list) &
      !$omp private(iat, izp, jat, kat, jzp, gam2, vec, r2, tmp, norm_cn, radi, radj)

      !$omp do schedule(runtime)
      do iat = 1, mol%nat
         izp = mol%id(iat)
         ! Effective charge width of i
         norm_cn = 1.0_wp / self%avg_cn(izp)**self%norm_exp
         radi = self%rad(izp) * (1.0_wp - self%kcnrad * cache%cn(iat) * norm_cn)

         do kat = list%inl(iat) + 1, list%inl(iat + 1) - 1
            jat = list%nlat(kat)
            jzp = mol%id(jat)
            vec = mol%xyz(:, jat) - mol%xyz(:, iat)
            r2 = vec(1)**2 + vec(2)**2 + vec(3)**2
            ! Effective charge width of j
            norm_cn = cache%cn(jat) / self%avg_cn(jzp)**self%norm_exp
            radj = self%rad(jzp) * (1.0_wp - self%kcnrad * norm_cn)
            ! Coulomb interaction of Gaussian charges
            gam2 = 1.0_wp / (radi**2 + radj**2)
            tmp = erf(sqrt(r2 * gam2)) / sqrt(r2) * cache%clist(kat)
            cache%alist(kat) = tmp
         end do

         ! Effective hardness
         tmp = self%eta(izp) + self%kqeta(izp) * cache%qloc(iat) + sqrt2pi / radi
         cache%adiag(iat) = tmp * cache%cdiag(iat) + 1.0_wp
      end do
      !$omp end do
      !$omp end parallel

   end subroutine get_amat_0d_list

!> Build the Coulomb matrix for a periodic system (3D) using Ewald summation and bond capacitance.
   subroutine get_amat_3d(self, mol, cache)
      !> EEQBC model type
      class(eeqbc_model), intent(in) :: self
      !> Molecular structure data
      type(structure_type), intent(in) :: mol
      !> Multicharge cache
      type(mchrg_cache), intent(inout) :: cache

      integer :: iat, jat, izp, jzp, img
      real(wp) :: vec(3), r1, gam, dtmp, ctmp, capi, capj, radi, radj, norm_cn, rvdw, wsw
      real(wp), allocatable :: dtrans(:, :)

      ! Thread-private array for reduction
      real(wp), allocatable :: amat_local(:, :)

      call get_dir_trans(mol, dtrans, cutoff)

      cache%amat(:, :) = 0.0_wp


      !$omp parallel default(none) &
      !$omp shared(cache, mol, self, dtrans)  &
      !$omp private(iat, izp, jat, jzp, gam, vec, dtmp, ctmp, norm_cn) &
      !$omp private(radi, radj, capi, capj, rvdw, r1, wsw, amat_local)
      allocate(amat_local, source=cache%amat)
      !$omp do schedule(runtime)
      do iat = 1, mol%nat
         izp = mol%id(iat)
         ! Effective charge width of i
         norm_cn = 1.0_wp / self%avg_cn(izp)**self%norm_exp
         radi = self%rad(izp) * (1.0_wp - self%kcnrad * cache%cn(iat) * norm_cn)
         capi = self%cap(izp)
         do jat = 1, iat - 1
            jzp = mol%id(jat)
            ! vdw distance in Angstrom (approximate factor 2)
            rvdw = self%rvdw(izp, jzp)
            ! Effective charge width of j
            norm_cn = cache%cn(jat) / self%avg_cn(jzp)**self%norm_exp
            radj = self%rad(jzp) * (1.0_wp - self%kcnrad * norm_cn)
            capj = self%cap(jzp)
            ! Coulomb interaction of Gaussian charges
            gam = 1.0_wp / sqrt(radi**2 + radj**2)
            wsw = 1.0_wp / real(cache%wsc%nimg(jat, iat), wp)
            do img = 1, cache%wsc%nimg(jat, iat)
               vec = mol%xyz(:, jat) - mol%xyz(:, iat) + cache%wsc%trans(:, cache%wsc%tridx(img, jat, iat))
               call get_amat_dir_3d(vec, gam, dtrans, self%kbc, rvdw, capi, capj, dtmp)
               amat_local(jat, iat) = amat_local(jat, iat) + dtmp * wsw
               amat_local(iat, jat) = amat_local(iat, jat) + dtmp * wsw
            end do
         end do

         ! diagonal Coulomb interaction terms
         gam = 1.0_wp / sqrt(2.0_wp * radi**2)
         rvdw = self%rvdw(izp, izp)
         wsw = 1.0_wp / real(cache%wsc%nimg(iat, iat), wp)
         do img = 1, cache%wsc%nimg(iat, iat)
            vec = cache%wsc%trans(:, cache%wsc%tridx(img, iat, iat))
            call get_amat_dir_3d(vec, gam, dtrans, self%kbc, rvdw, capi, capi, dtmp)
            amat_local(iat, iat) = amat_local(iat, iat) + dtmp * wsw
         end do

         ! Effective hardness
         dtmp = self%eta(izp) + self%kqeta(izp) * cache%qloc(iat) + sqrt2pi / radi
         amat_local(iat, iat) = amat_local(iat, iat) + cache%cmat(iat, iat) * dtmp + 1.0_wp
      end do
      !$omp end do
      !$omp critical (get_amat_3d_)
      cache%amat(:, :) = cache%amat + amat_local
      !$omp end critical (get_amat_3d_)
      deallocate(amat_local)
      !$omp end parallel

      if (size(cache%amat, 1) == mol%nat + 1) then
         cache%amat(mol%nat + 1, 1:mol%nat + 1) = 1.0_wp
         cache%amat(1:mol%nat + 1, mol%nat + 1) = 1.0_wp
         cache%amat(mol%nat + 1, mol%nat + 1) = 0.0_wp
      end if

   end subroutine get_amat_3d

!> Build the Coulomb matrix for a periodic system using CSR adjacency list.
   subroutine get_amat_3d_list(self, mol, list, cache)
      class(eeqbc_model), intent(in) :: self
      type(structure_type), intent(in) :: mol
      type(csr_list), intent(in) :: list
      !> Multicharge cache
      type(mchrg_cache), intent(inout) :: cache

      integer :: iat, jat, izp, jzp, img, kat
      real(wp) :: vec(3), gam, dtmp, capi, capj, radi, radj, norm_cn, rvdw, wsw
      real(wp) :: atmp, adiag_tmp
      real(wp), allocatable :: dtrans(:, :)

      call get_dir_trans(mol, dtrans, cutoff)

      ! Zero out global shared target arrays upfront
      cache%alist(:) = 0.0_wp
      cache%adiag(:) = 0.0_wp

      !$omp parallel default(none) &
      !$omp shared(cache, mol, self, list, dtrans) &
      !$omp private(iat, izp, kat, jat, jzp, gam, vec, dtmp, norm_cn, radi, radj) &
      !$omp private(capi, capj, rvdw, wsw, img, atmp, adiag_tmp)

      !$omp do schedule(runtime)
      do iat = 1, mol%nat
         izp = mol%id(iat)
         norm_cn = 1.0_wp / self%avg_cn(izp)**self%norm_exp
         radi = self%rad(izp) * (1.0_wp - self%kcnrad * cache%cn(iat) * norm_cn)
         capi = self%cap(izp)

         ! Initialize scalar accumulator for the diagonal of this atom
         adiag_tmp = 0.0_wp

         do kat = list%inl(iat) + 1, list%inl(iat+1) - 1
            if (list%wsc%nimg_list(kat) == 0) cycle

            jat = list%nlat(kat)
            jzp = mol%id(jat)
            capj = self%cap(jzp)
            rvdw = self%rvdw(izp, jzp)
            wsw = 1.0_wp / real(list%wsc%nimg_list(kat), wp)

            norm_cn = cache%cn(jat) / self%avg_cn(jzp)**self%norm_exp
            radj = self%rad(jzp) * (1.0_wp - self%kcnrad * norm_cn)
            gam = 1.0_wp / sqrt(radi**2 + radj**2)

            ! Accumulate image contributions in a local scalar
            atmp = 0.0_wp
            do img = list%wsc%itr_list(kat), list%wsc%itr_list(kat+1) - 1
               vec = mol%xyz(:, jat) - mol%xyz(:, iat) + list%wsc%trans(:, list%wsc%tridx_list(img))
               call get_amat_dir_3d(vec, gam, dtrans, self%kbc, rvdw, capi, capj, dtmp)
               atmp = atmp + dtmp * wsw
            end do

            cache%alist(kat) = atmp
         end do

         ! Diagonal Coulomb interaction terms (Self-Image)
         gam = 1.0_wp / sqrt(2.0_wp * radi**2)
         rvdw = self%rvdw(izp, izp)

         ! Execute only if periodic self-images exist (nimg > 1, since index 1 is identity R=0)
         if (list%wsc%nimg_list(list%inl(iat)) > 1) then
            wsw = 1.0_wp / real(list%wsc%nimg_list(list%inl(iat)) - 1, wp)

            ! Skip identity image at list%wsc%itr_list(list%inl(iat)) and loop through R /= 0 images
            do img = list%wsc%itr_list(list%inl(iat)) + 1, list%wsc%itr_list(list%inl(iat) + 1) - 1
               vec = list%wsc%trans(:, list%wsc%tridx_list(img))
               call get_amat_dir_3d(vec, gam, dtrans, self%kbc, rvdw, capi, capi, dtmp)
               adiag_tmp = adiag_tmp + dtmp * wsw
            end do
         end if

         ! Effective hardness
         dtmp = self%eta(izp) + self%kqeta(izp) * cache%qloc(iat) + sqrt2pi / radi

         ! Single safe direct write (iat is thread-exclusive)
         cache%adiag(iat) = adiag_tmp + cache%cdiag(iat) * dtmp + 1.0_wp
      end do
      !$omp end do
      !$omp end parallel

   end subroutine get_amat_3d_list

!> Real-space contribution to the Coulomb matrix for the EEQBC model.
   subroutine get_amat_dir_3d(rij, gam, trans, kbc, rvdw, capi, capj, amat)
      !> Distance vector between two atoms (including lattice translation)
      real(wp), intent(in) :: rij(3)
      !> Gaussian width parameter
      real(wp), intent(in) :: gam
      !> Direct lattice translation vectors (3 × N)
      real(wp), intent(in) :: trans(:, :)
      !> Bond capacitance exponent
      real(wp), intent(in) :: kbc
      !> Van der Waals distance
      real(wp), intent(in) :: rvdw
      !> Bond capacitance of atom i
      real(wp), intent(in) :: capi
      !> Bond capacitance of atom j
      real(wp), intent(in) :: capj
      !> Output contribution to the Coulomb matrix
      real(wp), intent(out) :: amat

      integer :: itr
      real(wp) :: vec(3), r1, tmp, ctmp

      amat = 0.0_wp

      do itr = 1, size(trans, 2)
         vec(:) = rij + trans(:, itr)
         r1 = norm2(vec)
         if (r1 < eps) cycle
         call get_cpair(kbc, ctmp, r1, rvdw, capi, capj)
         tmp = -ctmp * erf(gam * r1) / r1
         amat = amat + tmp
      end do

   end subroutine get_amat_dir_3d

!> Compute derivatives of the Coulomb matrix (multiplied by the charge vector) for the EEQBC model.
   subroutine get_coulomb_derivs(self, mol, ndim, cache, list)
      !> EEQBC model type
      class(eeqbc_model), intent(in) :: self
      !> Structure type
      type(structure_type), intent(in) :: mol
      !> System size
      integer, intent(in) :: ndim
      !> Multicharge cache
      type(mchrg_cache), intent(inout) :: cache
      !> Multicharge neighbourlist type
      type(csr_list), intent(in), optional :: list

      real(wp), allocatable :: atrace(:,:)

      integer :: iat

      if (.not. allocated(cache%dcndr)) then
         allocate(cache%dcndr(3, mol%nat, mol%nat))
      end if
      if (.not. allocated(cache%dcndL)) then
         allocate(cache%dcndL(3, 3, mol%nat))
         call self%ncoord%get_coordination_number(mol, cache%trans, cache%cn, dcndr=cache%dcndr, dcndL=cache%dcndL)
      end if
      if (.not. allocated(cache%dqlocdr)) then
         allocate(cache%dqlocdr(3, mol%nat, mol%nat))
      end if
      if (.not. allocated(cache%dqlocdL)) then
         allocate(cache%dqlocdL(3, 3, mol%nat))
         call self%local_charge(mol, cache%trans, qloc=cache%qloc, dqlocdr=cache%dqlocdr, dqlocdL=cache%dqlocdL)
      end if


      if (present(list)) then
         if (.not. allocated(cache%dadrij)) then
            allocate(cache%dadrij(3, size(list%nlat)))
         end if
         if (.not. allocated(cache%dadrji)) then
            allocate(cache%dadrji(3, size(list%nlat)))
         end if
         if (.not. allocated(cache%dadrdiag)) then
            allocate(cache%dadrdiag(3, mol%nat))
         end if
         if (.not. allocated(cache%dadL)) then
            allocate(cache%dadL(3, 3, mol%nat))
         end if

         if (any(mol%periodic)) then
            call get_damat_3d_list(self, mol, list, cache)
         else
            call get_damat_0d_list(self, mol, list, cache)
         end if
      else
         allocate(atrace(3, mol%nat), source = 0.0_wp)
         if (.not. allocated(cache%dadr)) allocate(cache%dadr(3, mol%nat, ndim))
         if (.not. allocated(cache%dadL)) allocate(cache%dadL(3, 3, ndim))

         if (any(mol%periodic)) then
            call get_damat_3d(self, mol, cache, atrace)
         else
            call get_damat_0d(self, mol, cache, atrace)
         end if
         do iat = 1, mol%nat
            cache%dadr(:, iat, iat) = atrace(:, iat) + cache%dadr(:, iat, iat)
         end do
      end if

   end subroutine get_coulomb_derivs

!> Build derivatives of the Coulomb matrix for a non-periodic system.
   subroutine get_damat_0d(self, mol, cache, atrace)
      !> EEQBC model type
      class(eeqbc_model), intent(in) :: self
      !> Molecular structure data
      type(structure_type), intent(in) :: mol
      !> Multicharge cache
      type(mchrg_cache), intent(inout) :: cache
      !> Trace-like array for diagonal contributions
      real(wp), intent(out) :: atrace(:, :)

      integer :: iat, jat, izp, jzp
      real(wp) :: vec(3), r2, gam, arg, dtmp, norm_cn
      real(wp) :: radi, radj, dradi, dradj, dG(3), dS(3, 3), dgamdL(3, 3)
      real(wp), allocatable :: dgamdr(:, :)

      ! Thread-private arrays for reduction
      real(wp), allocatable :: atrace_local(:, :)
      real(wp), allocatable :: dadr_local(:, :, :), dadL_local(:, :, :)

      allocate(dgamdr(3, mol%nat))

      atrace(:, :) = 0.0_wp
      cache%dadr(:, :, :) = 0.0_wp
      cache%dadL(:, :, :) = 0.0_wp

      !$omp parallel default(none) &
      !$omp shared(atrace, cache, mol, self) &
      !$omp private(iat, izp, jat, jzp, gam, vec, r2, dtmp, norm_cn, arg) &
      !$omp private(radi, radj, dradi, dradj, dgamdr, dgamdL, dG, dS) &
      !$omp private(atrace_local, dadr_local, dadL_local)
      allocate(atrace_local, source=atrace)
      allocate(dadr_local, source=cache%dadr)
      allocate(dadL_local, source=cache%dadL)
      !$omp do schedule(runtime)
      do iat = 1, mol%nat
         izp = mol%id(iat)
         ! Effective charge width of i
         norm_cn = 1.0_wp / self%avg_cn(izp)**self%norm_exp
         radi = self%rad(izp) * (1.0_wp - self%kcnrad * cache%cn(iat) * norm_cn)
         dradi = -self%rad(izp) * self%kcnrad * norm_cn
         do jat = 1, iat - 1
            jzp = mol%id(jat)
            vec = mol%xyz(:, jat) - mol%xyz(:, iat)
            r2 = vec(1)**2 + vec(2)**2 + vec(3)**2
            ! Effective charge width of j
            norm_cn = 1.0_wp / self%avg_cn(jzp)**self%norm_exp
            radj = self%rad(jzp) * (1.0_wp - self%kcnrad * cache%cn(jat) * norm_cn)
            dradj = -self%rad(jzp) * self%kcnrad * norm_cn

            ! Coulomb interaction of Gaussian charges
            gam = 1.0_wp / sqrt(radi**2 + radj**2)
            dgamdr(:, :) = -(radi * dradi * cache%dcndr(:, :, iat) + radj * dradj * cache%dcndr(:, :, jat)) &
            & * gam**3.0_wp
            dgamdL(:, :) = -(radi * dradi * cache%dcndL(:, :, iat) + radj * dradj * cache%dcndL(:, :, jat)) &
            & * gam**3.0_wp

            ! Explicit derivative
            arg = gam * gam * r2
            dtmp = 2.0_wp * gam * exp(-arg) / (sqrtpi * r2) &
            & - erf(sqrt(arg)) / (r2 * sqrt(r2))
            dG(:) = dtmp * vec
            dS(:, :) = spread(dG, 1, 3) * spread(vec, 2, 3)
            atrace_local(:, iat) = -dG * cache%vrhs(jat) * cache%cmat(jat, iat) + atrace_local(:, iat)
            atrace_local(:, jat) = +dG * cache%vrhs(iat) * cache%cmat(iat, jat) + atrace_local(:, jat)
            dadr_local(:, iat, jat) = -dG * cache%vrhs(iat) * cache%cmat(iat, jat) + dadr_local(:, iat, jat)
            dadr_local(:, jat, iat) = +dG * cache%vrhs(jat) * cache%cmat(jat, iat) + dadr_local(:, jat, iat)
            dadL_local(:, :, iat) = +dS * cache%vrhs(jat) * cache%cmat(jat, iat) + dadL_local(:, :, iat)
            dadL_local(:, :, jat) = +dS * cache%vrhs(iat) * cache%cmat(iat, jat) + dadL_local(:, :, jat)

            ! Effective charge width derivative
            dtmp = 2.0_wp * exp(-arg) / (sqrtpi)
            atrace_local(:, iat) = -dtmp * cache%vrhs(jat) * dgamdr(:, jat) * cache%cmat(jat, iat) + atrace_local(:, iat)
            atrace_local(:, jat) = -dtmp * cache%vrhs(iat) * dgamdr(:, iat) * cache%cmat(iat, jat) + atrace_local(:, jat)
            dadr_local(:, iat, jat) = +dtmp * cache%vrhs(iat) * dgamdr(:, iat) * cache%cmat(iat, jat) + dadr_local(:, iat, jat)
            dadr_local(:, jat, iat) = +dtmp * cache%vrhs(jat) * dgamdr(:, jat) * cache%cmat(jat, iat) + dadr_local(:, jat, iat)
            dadL_local(:, :, iat) = +dtmp * cache%vrhs(jat) * dgamdL(:, :) * cache%cmat(jat, iat) + dadL_local(:, :, iat)
            dadL_local(:, :, jat) = +dtmp * cache%vrhs(iat) * dgamdL(:, :) * cache%cmat(iat, jat) + dadL_local(:, :, jat)

            ! Capacitance derivative off-diagonal
            dtmp = erf(sqrt(r2) * gam) / (sqrt(r2))
            atrace_local(:, iat) = -dtmp * cache%vrhs(jat) * cache%dcdr(:, jat, iat) + atrace_local(:, iat)
            atrace_local(:, jat) = -dtmp * cache%vrhs(iat) * cache%dcdr(:, iat, jat) + atrace_local(:, jat)
            dadr_local(:, iat, jat) = +dtmp * cache%vrhs(iat) * cache%dcdr(:, iat, jat) + dadr_local(:, iat, jat)
            dadr_local(:, jat, iat) = +dtmp * cache%vrhs(jat) * cache%dcdr(:, jat, iat) + dadr_local(:, jat, iat)
            dadL_local(:, :, iat) = -dtmp * cache%vrhs(jat) * spread(cache%dcdr(:, iat, jat), 2, 3) * spread(vec, 1, 3) &
            & + dadL_local(:, :, iat)
            dadL_local(:, :, jat) = -dtmp * cache%vrhs(iat) * spread(cache%dcdr(:, iat, jat), 2, 3) * spread(vec, 1, 3) &
            & + dadL_local(:, :, jat)

            ! Capacitance derivative diagonal
            dtmp = (self%eta(izp) + self%kqeta(izp) * cache%qloc(iat) + sqrt2pi / radi) * cache%vrhs(iat)
            dadr_local(:, jat, iat) = -dtmp * cache%dcdr(:, jat, iat) + dadr_local(:, jat, iat)

            dtmp = (self%eta(jzp) + self%kqeta(jzp) * cache%qloc(jat) + sqrt2pi / radj) * cache%vrhs(jat)
            dadr_local(:, iat, jat) = -dtmp * cache%dcdr(:, iat, jat) + dadr_local(:, iat, jat)
         end do

         ! Hardness derivative
         dtmp = self%kqeta(izp) * cache%vrhs(iat) * cache%cmat(iat, iat)
         dadr_local(:, :, iat) = +dtmp * cache%dqlocdr(:, :, iat) + dadr_local(:, :, iat)
         dadL_local(:, :, iat) = +dtmp * cache%dqlocdL(:, :, iat) + dadL_local(:, :, iat)

         ! Effective charge width derivative
         dtmp = -sqrt2pi * dradi / (radi**2) * cache%vrhs(iat) * cache%cmat(iat, iat)
         dadr_local(:, :, iat) = +dtmp * cache%dcndr(:, :, iat) + dadr_local(:, :, iat)
         dadL_local(:, :, iat) = +dtmp * cache%dcndL(:, :, iat) + dadL_local(:, :, iat)

         ! Capacitance derivative
         dtmp = (self%eta(izp) + self%kqeta(izp) * cache%qloc(iat) + sqrt2pi / radi) * cache%vrhs(iat)
         dadr_local(:, iat, iat) = +dtmp * cache%dcdr(:, iat, iat) + dadr_local(:, iat, iat)
         dadL_local(:, :, iat) = +dtmp * cache%dcdL(:, :, iat) + dadL_local(:, :, iat)

      end do
      !$omp end do
      !$omp critical (get_damat_0d_)
      atrace(:, :) = atrace + atrace_local
      cache%dadr(:, :, :) = cache%dadr + dadr_local
      cache%dadL(:, :, :) = cache%dadL + dadL_local
      !$omp end critical (get_damat_0d_)
      deallocate(dadL_local, dadr_local, atrace_local)
      !$omp end parallel

   end subroutine get_damat_0d

!> Build derivatives of the Coulomb matrix for a non‑periodic system using neighbour list.
   subroutine get_damat_0d_list(self, mol, list, cache)
      !> EEQBC model type
      class(eeqbc_model), intent(in) :: self
      !> Molecular structure data
      type(structure_type), intent(in) :: mol
      !> Multicharge neighbourlist type
      type(csr_list), intent(in) :: list
      !> Multicharge cache
      type(mchrg_cache), intent(inout) :: cache

      integer :: iat, jat, kat, izp, jzp, start_kat, finish_kat
      real(wp) :: vec(3), r2, gam, arg, dtmp, norm_cn
      real(wp) :: radi, radj, dradi, dradj, dG(3), dS(3, 3), dgamdL(3, 3)
      real(wp) :: pre_i, pre_j, dgam_pre, dgami(3), dgamj(3)

      ! Thread-private arrays for reduction
      real(wp), allocatable :: dadrij_local(:, :), dadrji_local(:, :)
      real(wp), allocatable :: dadrdiag_local(:, :), dadL_local(:, :, :)

      cache%dadrdiag = 0.0_wp
      cache%dadrij = 0.0_wp
      cache%dadrji = 0.0_wp
      cache%dadL = 0.0_wp

      !$omp parallel default(none) &
      !$omp shared(cache, mol, list, self) &
      !$omp private(iat, kat, izp, jat, jzp, gam, dgami, dgamj, dgamdL, vec, r2, dtmp, norm_cn, arg) &
      !$omp private(start_kat, finish_kat, radi, radj, dradi, dradj, dG, dS, pre_i, pre_j, dgam_pre) &
      !$omp private(dadrij_local, dadrji_local, dadrdiag_local, dadL_local)

      allocate(dadrij_local, source=cache%dadrij)
      allocate(dadrji_local, source=cache%dadrji)
      allocate(dadrdiag_local, source=cache%dadrdiag)
      allocate(dadL_local, source=cache%dadL)

      !$omp do schedule(runtime)
      do iat = 1, mol%nat
         izp = mol%id(iat)
         norm_cn = 1.0_wp / self%avg_cn(izp)**self%norm_exp
         radi = self%rad(izp) * (1.0_wp - self%kcnrad * cache%cn(iat) * norm_cn)
         dradi = -self%rad(izp) * self%kcnrad * norm_cn

         start_kat = list%inl(iat) + 1
         finish_kat = list%inl(iat+1) - 1

         do kat = start_kat, finish_kat
            jat = list%nlat(kat)
            jzp = mol%id(jat)
            vec = mol%xyz(:, jat) - mol%xyz(:, iat)
            r2 = dot_product(vec, vec)

            norm_cn = 1.0_wp / self%avg_cn(jzp)**self%norm_exp
            radj = self%rad(jzp) * (1.0_wp - self%kcnrad * cache%cn(jat) * norm_cn)
            dradj = -self%rad(jzp) * self%kcnrad * norm_cn

            gam = 1.0_wp / sqrt(radi**2 + radj**2)
            dgami(:) = -(radi * dradi * cache%dcndrdiag(:, iat) + radj * dradj * cache%dcndrij(:, kat)) &
            & * gam**3.0_wp
            dgamj(:) = -(radi * dradi * cache%dcndrji(:, kat) + radj * dradj * cache%dcndrdiag(:, jat)) &
            & * gam**3.0_wp
            dgamdL(:, :) = -(radi * dradi * cache%dcndL(:, :, iat) + radj * dradj * cache%dcndL(:, :, jat)) &
            & * gam**3.0_wp
            arg = gam * gam * r2

            ! 1. Explicit Geometry Derivative (Coulomb kernel)
            dtmp = 2.0_wp * gam * exp(-arg) / (sqrtpi * r2) - erf(sqrt(arg)) / (r2 * sqrt(r2))
            dG = dtmp * vec
            dS = spread(dG, 1, 3) * spread(vec, 2, 3)

            ! Equivalent to atrace logic (diagonal updates)
            dadrdiag_local(:, iat) = dadrdiag_local(:, iat) - dG * cache%vrhs(jat) * cache%clist(kat)
            dadrdiag_local(:, jat) = dadrdiag_local(:, jat) + dG * cache%vrhs(iat) * cache%clist(kat)

            ! Off-diagonal updates
            dadrij_local(:, kat) = dadrij_local(:, kat) - dG * cache%vrhs(iat) * cache%clist(kat)
            dadrji_local(:, kat) = dadrji_local(:, kat) + dG * cache%vrhs(jat) * cache%clist(kat)

            dadL_local(:, :, iat)  = dadL_local(:, :, iat)  + dS * cache%vrhs(jat) * cache%clist(kat)
            dadL_local(:, :, jat)  = dadL_local(:, :, jat)  + dS * cache%vrhs(iat) * cache%clist(kat)

            ! 2. Effective charge width derivative

            ! Update dadrdiag (atrace equivalent)
            dtmp = 2.0_wp * exp(-arg) / (sqrtpi)
            dadrdiag_local(:, iat) = dadrdiag_local(:, iat) - dtmp * dgamj(:) * cache%vrhs(jat) * cache%clist(kat)
            dadrdiag_local(:, jat) = dadrdiag_local(:, jat) - dtmp * dgami(:) * cache%vrhs(iat) * cache%clist(kat)

            ! Update dadrij/ji (off-diagonal)
            dadrij_local(:, kat) = dadrij_local(:, kat) + dtmp * dgami(:) * cache%vrhs(iat) * cache%clist(kat)
            dadrji_local(:, kat) = dadrji_local(:, kat) + dtmp * dgamj(:) * cache%vrhs(jat) * cache%clist(kat)

            ! Lattice derivative
            dadL_local(:, :, iat) = +dtmp * cache%vrhs(jat) * dgamdL(:, :) * cache%clist(kat) + dadL_local(:, :, iat)
            dadL_local(:, :, jat) = +dtmp * cache%vrhs(iat) * dgamdL(:, :) * cache%clist(kat) + dadL_local(:, :, jat)


            ! 3. Capacitance derivative off-diagonal
            dtmp = erf(sqrt(r2) * gam) / sqrt(r2)
            dadrdiag_local(:, iat) = dadrdiag_local(:, iat) - dtmp * cache%vrhs(jat) * cache%dcdrji(:, kat)
            dadrdiag_local(:, jat) = dadrdiag_local(:, jat) - dtmp * cache%vrhs(iat) * cache%dcdrij(:, kat)

            dadrij_local(:, kat) = dadrij_local(:, kat) + dtmp * cache%vrhs(iat) * cache%dcdrij(:, kat)
            dadrji_local(:, kat) = dadrji_local(:, kat) + dtmp * cache%vrhs(jat) * cache%dcdrji(:, kat)

            dadL_local(:, :, iat) = dadL_local(:, :, iat) - dtmp * cache%vrhs(jat) * &
               spread(cache%dcdrij(:, kat), 2, 3) * spread(vec, 1, 3)
            dadL_local(:, :, jat) = dadL_local(:, :, jat) - dtmp * cache%vrhs(iat) * &
               spread(cache%dcdrij(:, kat), 2, 3) * spread(vec, 1, 3)

            ! 4. Capacitance derivative diagonal contribution
            dtmp = (self%eta(izp) + self%kqeta(izp) * cache%qloc(iat) + sqrt2pi / radi) * cache%vrhs(iat)
            dadrji_local(:, kat) = -dtmp * cache%dcdrji(:, kat) + dadrji_local(:, kat)

            dtmp = (self%eta(jzp) + self%kqeta(jzp) * cache%qloc(jat) + sqrt2pi / radj) * cache%vrhs(jat)
            dadrij_local(:, kat) = -dtmp * cache%dcdrij(:, kat) + dadrij_local(:, kat)

            ! 5. Hardness and coordination-dependent diagonal corrections
            pre_i = self%kqeta(izp) * cache%vrhs(iat) * cache%cdiag(iat)
            pre_j = self%kqeta(jzp) * cache%vrhs(jat) * cache%cdiag(jat)

            dadrij_local(:, kat) = dadrij_local(:, kat) + pre_j * cache%dqlocdrij(:, kat)
            dadrji_local(:, kat) = dadrji_local(:, kat) + pre_i * cache%dqlocdrji(:, kat)

            ! Effective charge width diagonal derivative
            pre_i = -sqrt2pi * dradi / (radi**2) * cache%vrhs(iat) * cache%cdiag(iat)
            pre_j = -sqrt2pi * dradj / (radj**2) * cache%vrhs(jat) * cache%cdiag(jat)

            dadrij_local(:, kat) = dadrij_local(:, kat) + pre_j * cache%dcndrij(:, kat)
            dadrji_local(:, kat) = dadrji_local(:, kat) + pre_i * cache%dcndrji(:, kat)

         end do

         ! 5. Hardness and coordination-dependent diagonal corrections (Diagonal-only)
         dtmp = self%kqeta(izp) * cache%vrhs(iat) * cache%cdiag(iat)
         dadrdiag_local(:, iat) = dadrdiag_local(:, iat) + dtmp * cache%dqlocdrdiag(:, iat)
         dadL_local(:, :, iat) = dadL_local(:, :, iat) + dtmp * cache%dqlocdL(:, :, iat)

         ! Effective charge width diagonal derivative (Diagonal-only)
         dtmp = -sqrt2pi * dradi / (radi**2) * cache%vrhs(iat) * cache%cdiag(iat)
         dadrdiag_local(:, iat) = dadrdiag_local(:, iat) + dtmp * cache%dcndrdiag(:, iat)
         dadL_local(:, :, iat) = dadL_local(:, :, iat) + dtmp * cache%dcndL(:, :, iat)

         ! 6. Intrinsic capacitance derivative
         dtmp = (self%eta(izp) + self%kqeta(izp) * cache%qloc(iat) + sqrt2pi / radi) * cache%vrhs(iat)
         dadrdiag_local(:, iat) = dadrdiag_local(:, iat) + dtmp * cache%dcdrdiag(:, iat)
         dadL_local(:, :, iat) = dadL_local(:, :, iat) + dtmp * cache%dcdL(:, :, iat)
      end do
      !$omp end do

      !$omp critical (get_damat_0d_list_)
      cache%dadrdiag = cache%dadrdiag + dadrdiag_local
      cache%dadrij = cache%dadrij + dadrij_local
      cache%dadrji = cache%dadrji + dadrji_local
      cache%dadL = cache%dadL + dadL_local
      !$omp end critical (get_damat_0d_list_)

      deallocate(dadL_local, dadrij_local, dadrji_local, dadrdiag_local)
      !$omp end parallel

   end subroutine get_damat_0d_list

!> Build derivatives of the Coulomb matrix for a periodic system.
   subroutine get_damat_3d(self, mol, cache, atrace)
      !> EEQBC model type
      class(eeqbc_model), intent(in) :: self
      !> Molecular structure data
      type(structure_type), intent(in) :: mol
      !> Multicharge cache
      type(mchrg_cache), intent(inout) :: cache
      !> Trace-like array for diagonal contributions
      real(wp), intent(out) :: atrace(:, :)

      integer :: iat, jat, izp, jzp, img
      real(wp) :: vec(3), r2, gam, arg, dtmp, norm_cn, rvdw, wsw, dgam
      real(wp) :: radi, radj, dradi, dradj, dG(3), dS(3, 3)
      real(wp) :: dgamdL(3, 3), capi, capj
      real(wp), allocatable :: dgamdr(:, :), dtrans(:, :)

      ! Thread-private arrays for reduction
      real(wp), allocatable :: atrace_local(:, :)
      real(wp), allocatable :: dadr_local(:, :, :), dadL_local(:, :, :)

      call get_dir_trans(mol, dtrans, cutoff)

      allocate(dgamdr(3, mol%nat))

      atrace(:, :) = 0.0_wp
      cache%dadr(:, :, :) = 0.0_wp
      cache%dadL(:, :, :) = 0.0_wp

      !$omp parallel default(none) &
      !$omp shared(self, mol, cache, atrace, dtrans) &
      !$omp private(iat, izp, jat, jzp, img, gam, vec, r2, dtmp, norm_cn, arg, rvdw) &
      !$omp private(radi, radj, dradi, dradj, capi, capj, dgamdr, dgamdL, dG, dS, wsw) &
      !$omp private(dgam, dadr_local, dadL_local, atrace_local)
      allocate(atrace_local, source=atrace)
      allocate(dadr_local, source=cache%dadr)
      allocate(dadL_local, source=cache%dadL)
      !$omp do schedule(runtime)
      do iat = 1, mol%nat
         izp = mol%id(iat)
         ! Effective charge width of i
         norm_cn = 1.0_wp / self%avg_cn(izp)**self%norm_exp
         radi = self%rad(izp) * (1.0_wp - self%kcnrad * cache%cn(iat) * norm_cn)
         dradi = -self%rad(izp) * self%kcnrad * norm_cn
         capi = self%cap(izp)
         do jat = 1, iat - 1
            jzp = mol%id(jat)
            capj = self%cap(jzp)
            rvdw = self%rvdw(izp, jzp)

            ! Effective charge width of j
            norm_cn = 1.0_wp / self%avg_cn(jzp)**self%norm_exp
            radj = self%rad(jzp) * (1.0_wp - self%kcnrad * cache%cn(jat) * norm_cn)
            dradj = -self%rad(jzp) * self%kcnrad * norm_cn

            ! Coulomb interaction of Gaussian charges
            gam = 1.0_wp / sqrt(radi**2 + radj**2)
            dgamdr(:, :) = -(radi * dradi * cache%dcndr(:, :, iat) + radj * dradj * cache%dcndr(:, :, jat)) &
            & * gam**3.0_wp
            dgamdL(:, :) = -(radi * dradi * cache%dcndL(:, :, iat) + radj * dradj * cache%dcndL(:, :, jat)) &
            & * gam**3.0_wp

            wsw = 1.0_wp / real(cache%wsc%nimg(jat, iat), wp)
            do img = 1, cache%wsc%nimg(jat, iat)
               vec = mol%xyz(:, jat) - mol%xyz(:, iat) + cache%wsc%trans(:, cache%wsc%tridx(img, jat, iat))

               call get_damat_dir(vec, dtrans, capi, capj, rvdw, self%kbc, gam, dG, dS, dgam)
               dG = dG * wsw
               dS = dS * wsw
               dgam = dgam * wsw

               ! Explicit derivative
               atrace_local(:, iat) = -dG * cache%vrhs(jat) + atrace_local(:, iat)
               atrace_local(:, jat) = +dG * cache%vrhs(iat) + atrace_local(:, jat)
               dadr_local(:, iat, jat) = -dG * cache%vrhs(iat) + dadr_local(:, iat, jat)
               dadr_local(:, jat, iat) = +dG * cache%vrhs(jat) + dadr_local(:, jat, iat)
               dadL_local(:, :, jat) = +dS * cache%vrhs(iat) + dadL_local(:, :, jat)
               dadL_local(:, :, iat) = +dS * cache%vrhs(jat) + dadL_local(:, :, iat)

               ! Effective charge width derivative
               atrace_local(:, iat) = +dgam * cache%vrhs(jat) * dgamdr(:, jat) + atrace_local(:, iat)
               atrace_local(:, jat) = +dgam * cache%vrhs(iat) * dgamdr(:, iat) + atrace_local(:, jat)
               dadr_local(:, iat, jat) = -dgam * cache%vrhs(iat) * dgamdr(:, iat) + dadr_local(:, iat, jat)
               dadr_local(:, jat, iat) = -dgam * cache%vrhs(jat) * dgamdr(:, jat) + dadr_local(:, jat, iat)
               dadL_local(:, :, iat) = -dgam * cache%vrhs(jat) * dgamdL(:, :) + dadL_local(:, :, iat)
               dadL_local(:, :, jat) = -dgam * cache%vrhs(iat) * dgamdL(:, :) + dadL_local(:, :, jat)

               call get_damat_dc_dir(vec, dtrans, capi, capj, rvdw, self%kbc, gam, dG, dS)
               dG = dG * wsw
               dS = dS * wsw

               ! Capacitance derivative off-diagonal
               atrace_local(:, iat) = +cache%vrhs(jat) * dG(:) + atrace_local(:, iat)
               atrace_local(:, jat) = -cache%vrhs(iat) * dG(:) + atrace_local(:, jat)
               dadr_local(:, jat, iat) = -cache%vrhs(jat) * dG(:) + dadr_local(:, jat, iat)
               dadr_local(:, iat, jat) = +cache%vrhs(iat) * dG(:) + dadr_local(:, iat, jat)
               dadL_local(:, :, jat) = -cache%vrhs(iat) * dS(:, :) + dadL_local(:, :, jat)
               dadL_local(:, :, iat) = -cache%vrhs(jat) * dS(:, :) + dadL_local(:, :, iat)

               call get_dcpair_dir(self%kbc, vec, dtrans, rvdw, capi, capj, dG, dS)
               dG = dG * wsw

               ! Capacitance derivative diagonal
               dtmp = (self%eta(izp) + self%kqeta(izp) * cache%qloc(iat) + sqrt2pi / radi) * cache%vrhs(iat)
               dadr_local(:, jat, iat) = +dtmp * dG(:) + dadr_local(:, jat, iat)
               dtmp = (self%eta(jzp) + self%kqeta(jzp) * cache%qloc(jat) + sqrt2pi / radj) * cache%vrhs(jat)
               dadr_local(:, iat, jat) = -dtmp * dG(:) + dadr_local(:, iat, jat)
            end do
         end do

         ! diagonal explicit, charge width, and capacitance derivative terms
         gam = 1.0_wp / sqrt(2.0_wp * radi**2)
         dtmp = -sqrt2pi * dradi / (radi**2) * cache%vrhs(iat)
         rvdw = self%rvdw(izp, izp)
         wsw = 1.0_wp / real(cache%wsc%nimg(iat, iat), wp)
         do img = 1, cache%wsc%nimg(iat, iat)
            vec = cache%wsc%trans(:, cache%wsc%tridx(img, iat, iat))
            call get_damat_dir(vec, dtrans, capi, capi, rvdw, self%kbc, gam, dG, dS, dgam)
            dgam = dgam * wsw

            ! Explicit derivative
            dadL_local(:, :, iat) = +dS * wsw * cache%vrhs(iat) + dadL_local(:, :, iat)

            ! Effective charge width derivative
            atrace_local(:, iat) = +dtmp * cache%dcndr(:, iat, iat) * dgam + atrace_local(:, iat)
            dadr_local(:, iat, iat) = -dtmp * cache%dcndr(:, iat, iat) * dgam + dadr_local(:, iat, iat)
            dadL_local(:, :, iat) = -dtmp * cache%dcndL(:, :, iat) * dgam + dadL_local(:, :, iat)

            ! Capacitance derivative
            call get_damat_dc_dir(vec, dtrans, capi, capi, rvdw, self%kbc, gam, dG, dS)
            dadL_local(:, :, iat) = -cache%vrhs(iat) * dS * wsw + dadL_local(:, :, iat)
         end do

         ! Hardness derivative
         dtmp = self%kqeta(izp) * cache%vrhs(iat) * cache%cmat(iat, iat)
         dadr_local(:, :, iat) = +dtmp * cache%dqlocdr(:, :, iat) + dadr_local(:, :, iat)
         dadL_local(:, :, iat) = +dtmp * cache%dqlocdL(:, :, iat) + dadL_local(:, :, iat)

         ! Effective charge width derivative
         dtmp = -sqrt2pi * dradi / (radi**2) * cache%vrhs(iat) * cache%cmat(iat, iat)
         dadr_local(:, :, iat) = +dtmp * cache%dcndr(:, :, iat) + dadr_local(:, :, iat)
         dadL_local(:, :, iat) = +dtmp * cache%dcndL(:, :, iat) + dadL_local(:, :, iat)

         dtmp = (self%eta(izp) + self%kqeta(izp) * cache%qloc(iat) + sqrt2pi / radi) * cache%vrhs(iat)
         dadr_local(:, iat, iat) = +dtmp * cache%dcdr(:, iat, iat) + dadr_local(:, iat, iat)
         dadL_local(:, :, iat) = +dtmp * cache%dcdL(:, :, iat) + dadL_local(:, :, iat)

      end do
      !$omp end do
      !$omp critical (get_damat_3d_)
      atrace(:, :) = atrace + atrace_local
      cache%dadr(:, :, :) = cache%dadr + dadr_local
      cache%dadL(:, :, :) = cache%dadL + dadL_local
      !$omp end critical (get_damat_3d_)
      deallocate(dadL_local, dadr_local, atrace_local)
      !$omp end parallel

   end subroutine get_damat_3d

!> Build derivatives of the Coulomb matrix for a periodic system.
   subroutine get_damat_3d_list(self, mol, list, cache)
      !> EEQBC model type
      class(eeqbc_model), intent(in) :: self
      !> Molecular structure data
      type(structure_type), intent(in) :: mol
      !> Multicharge neighbourlist type
      type(csr_list), intent(in) :: list
      !> Multicharge cache
      type(mchrg_cache), intent(inout) :: cache

      integer :: iat, jat, izp, jzp, img, kat
      integer :: start_kat, finish_kat
      real(wp) :: vec(3), gam, norm_cn, rvdw, wsw, dgam, dtmp
      real(wp) :: radi, radj, dradi, dradj, dG(3), dS(3, 3)
      real(wp) :: dgamdL(3, 3), capi, capj
      real(wp) :: pre_i, pre_j
      real(wp), allocatable :: dtrans(:, :)

      ! Thread-private local arrays
      real(wp), allocatable :: dadrij_local(:, :), dadrji_local(:, :)
      real(wp), allocatable :: dadrdiag_local(:, :), dadL_local(:, :, :)

      call get_dir_trans(mol, dtrans, cutoff)

      cache%dadrdiag(:, :) = 0.0_wp
      cache%dadrij(:, :) = 0.0_wp
      cache%dadrji(:, :) = 0.0_wp
      cache%dadL(:, :, :) = 0.0_wp

      !$omp parallel default(none) &
      !$omp shared(self, mol, list, cache, dtrans) &
      !$omp private(iat, izp, jat, kat, start_kat, finish_kat, jzp, img, gam, vec, dtmp, norm_cn, rvdw) &
      !$omp private(radi, radj, dradi, dradj, capi, capj, dgamdL, dG, dS, wsw) &
      !$omp private(dgam, dadrij_local, dadrji_local, dadL_local, dadrdiag_local, pre_i, pre_j)

      allocate(dadrdiag_local, source=cache%dadrdiag)
      allocate(dadrij_local, source=cache%dadrij)
      allocate(dadrji_local, source=cache%dadrji)
      allocate(dadL_local, source=cache%dadL)

      !$omp do schedule(runtime)
      do iat = 1, mol%nat
         izp = mol%id(iat)
         norm_cn = 1.0_wp / self%avg_cn(izp)**self%norm_exp
         radi = self%rad(izp) * (1.0_wp - self%kcnrad * cache%cn(iat) * norm_cn)
         dradi = -self%rad(izp) * self%kcnrad * norm_cn
         capi = self%cap(izp)

         start_kat = list%inl(iat) + 1
         finish_kat = list%inl(iat+1) - 1

         do kat = start_kat, finish_kat
            jat = list%nlat(kat)
            jzp = mol%id(jat)
            capj = self%cap(jzp)
            rvdw = self%rvdw(izp, jzp)
            wsw = 1.0_wp / real(list%wsc%itr_list(kat+1) - list%wsc%itr_list(kat), wp)

            norm_cn = 1.0_wp / self%avg_cn(jzp)**self%norm_exp
            radj = self%rad(jzp) * (1.0_wp - self%kcnrad * cache%cn(jat) * norm_cn)
            dradj = -self%rad(jzp) * self%kcnrad * norm_cn

            gam = 1.0_wp / sqrt(radi**2 + radj**2)
            pre_i = -radi * dradi * gam**3.0_wp
            pre_j = -radj * dradj * gam**3.0_wp

            ! Lattice derivative of interaction width
            dgamdL(:, :) = (pre_i * cache%dcndL(:, :, iat) + pre_j * cache%dcndL(:, :, jat))

            do img = list%wsc%itr_list(kat), list%wsc%itr_list(kat+1) - 1
               vec = mol%xyz(:, jat) - mol%xyz(:, iat) + list%wsc%trans(:, list%wsc%tridx_list(img))

               call get_damat_dir(vec, dtrans, capi, capj, rvdw, self%kbc, gam, dG, dS, dgam)
               dG = dG * wsw
               dS = dS * wsw
               dgam = dgam * wsw

               ! 1. Geometry-only derivatives
               dadrdiag_local(:, iat) = dadrdiag_local(:, iat) - dG * cache%vrhs(jat)
               dadrdiag_local(:, jat) = dadrdiag_local(:, jat) + dG * cache%vrhs(iat)
               dadrij_local(:, kat) = dadrij_local(:, kat) - dG * cache%vrhs(iat)
               dadrji_local(:, kat) = dadrji_local(:, kat) + dG * cache%vrhs(jat)
               dadL_local(:, :, jat) = dadL_local(:, :, jat) + dS * cache%vrhs(iat)
               dadL_local(:, :, iat) = dadL_local(:, :, iat) + dS * cache%vrhs(jat)

               ! 2. CN-dependent interaction width (dgam/dr)
               ! Deriv w.r.t iat: pre_i * dcndrdiag(iat) + pre_j * dcndrlist(kat_rev)
               ! Note: We only add terms where atom iat or jat is the center of the derivative
               dadrdiag_local(:, iat) = dadrdiag_local(:, iat) + &
               & dgam * cache%vrhs(jat) * (pre_i * cache%dcndrdiag(:, iat))
               dadrdiag_local(:, jat) = dadrdiag_local(:, jat) + &
               & dgam * cache%vrhs(iat) * (pre_j * cache%dcndrdiag(:, jat))
               dadrij_local(:, kat) = dadrij_local(:, kat) - &
               & dgam * cache%vrhs(iat) * (pre_i * cache%dcndrij(:, kat))
               dadrji_local(:, kat) = dadrji_local(:, kat) + &
               & dgam * cache%vrhs(jat) * (pre_j * cache%dcndrji(:, kat))

               dadL_local(:, :, iat) = dadL_local(:, :, iat) - &
               & dgam * cache%vrhs(jat) * dgamdL
               dadL_local(:, :, jat) = dadL_local(:, :, jat) - &
               & dgam * cache%vrhs(iat) * dgamdL

               ! 3. Bond capacitance derivative (Off-diagonal)
               call get_damat_dc_dir(vec, dtrans, capi, capj, rvdw, self%kbc, gam, dG, dS)
               dG = dG * wsw
               dS = dS * wsw
               dadrdiag_local(:, iat) = dadrdiag_local(:, iat) + cache%vrhs(jat) * dG
               dadrdiag_local(:, jat) = dadrdiag_local(:, jat) - cache%vrhs(iat) * dG
               dadrij_local(:, kat) = dadrij_local(:, kat) + cache%vrhs(iat) * dG
               dadrji_local(:, kat) = dadrji_local(:, kat) - cache%vrhs(jat) * dG
               dadL_local(:, :, jat) = dadL_local(:, :, jat) - cache%vrhs(iat) * dS
               dadL_local(:, :, iat) = dadL_local(:, :, iat) - cache%vrhs(jat) * dS

               ! 4. Bond capacitance derivative (Diagonal contribution)
               call get_dcpair_dir(self%kbc, vec, dtrans, rvdw, capi, capj, dG, dS)
               dG = dG * wsw
               dtmp = (self%eta(jzp) + self%kqeta(jzp) * &
               & cache%qloc(jat) + sqrt2pi / radj) * cache%vrhs(jat)
               dadrij_local(:, kat) = dadrij_local(:, kat) - dtmp * dG * wsw
               dtmp = (self%eta(izp) + self%kqeta(izp) * &
               & cache%qloc(iat) + sqrt2pi / radj) * cache%vrhs(iat)
               dadrji_local(:, kat) = dadrji_local(:, kat) - dtmp * dG * wsw

            end do

            ! Self-image interaction (jat == iat)
            gam = 1.0_wp / sqrt(2.0_wp * radi**2)
            pre_i = -2.0_wp * radi * dradi * gam**3.0_wp
            dtmp = -sqrt2pi * dradi / (radi**2) * cache%vrhs(iat)

            wsw = 1.0_wp / real(list%wsc%itr_list(kat+1) - list%wsc%itr_list(kat), wp)
            do img = list%wsc%itr_list(kat), list%wsc%itr_list(kat+1) - 1
               vec = list%wsc%trans(:, list%wsc%tridx_list(img))
               call get_damat_dir(vec, dtrans, capi, capi, rvdw, self%kbc, gam, dG, dS, dgam)
               dgam = dgam * wsw

               dadL_local(:, :, iat) = dadL_local(:, :, iat) + (dS * cache%vrhs(iat)) * wsw
               dadL_local(:, :, iat) = dadL_local(:, :, iat) - (dtmp * cache%dcndL(:, :, iat) * dgam)

               ! Self-image width change w.r.t position i
               dadrdiag_local(:, iat) = dadrdiag_local(:, iat) + &
               & dgam * cache%vrhs(iat) * (pre_i * cache%dcndrdiag(:, iat))

               call get_damat_dc_dir(vec, dtrans, capi, capi, rvdw, self%kbc, gam, dG, dS)
               dadL_local(:, :, iat) = dadL_local(:, :, iat) - (cache%vrhs(iat) * dS) * wsw
            end do

         end do

         ! 5. Hardness and coordination-dependent diagonal corrections
         dtmp = self%kqeta(izp) * cache%vrhs(iat) * cache%cdiag(iat)
         dadrdiag_local(:, iat) = dadrdiag_local(:, iat) + dtmp * cache%dqlocdrdiag(:, iat)
         dadrij_local(:, start_kat:finish_kat) = dadrij_local(:, start_kat:finish_kat) + &
         &  dtmp * cache%dqlocdrij(:, start_kat:finish_kat)
         dadrji_local(:, start_kat:finish_kat) = dadrji_local(:, start_kat:finish_kat) + &
         &  dtmp * cache%dqlocdrji(:, start_kat:finish_kat)
         dadL_local(:, :, iat) = dadL_local(:, :, iat) + dtmp * cache%dqlocdL(:, :, iat)

         dtmp = -sqrt2pi * dradi / (radi**2) * cache%vrhs(iat) * cache%cdiag(iat)
         dadrdiag_local(:, iat) = dadrdiag_local(:, iat) + dtmp * cache%dcndrdiag(:, iat)
         dadrij_local(:, start_kat:finish_kat) = dadrij_local(:, start_kat:finish_kat) + &
         &  dtmp * cache%dcndrij(:, start_kat:finish_kat)
         dadrji_local(:, start_kat:finish_kat) = dadrji_local(:, start_kat:finish_kat) + &
         &  dtmp * cache%dcndrji(:, start_kat:finish_kat)
         dadL_local(:, :, iat) = dadL_local(:, :, iat) + dtmp * cache%dcndL(:, :, iat)

         ! 6. Intrinsic capacitance derivative
         dtmp = (self%eta(izp) + self%kqeta(izp) * cache%qloc(iat) + sqrt2pi / radi) * cache%vrhs(iat)
         dadrdiag_local(:, iat) = dadrdiag_local(:, iat) + dtmp * cache%dcdrdiag(:, iat)
         dadL_local(:, :, iat) = dadL_local(:, :, iat) + dtmp * cache%dcdL(:, :, iat)

      end do
      !$omp end do

      !$omp critical (get_damat_3d_list_)
      cache%dadrdiag = cache%dadrdiag + dadrdiag_local
      cache%dadrij = cache%dadrij + dadrij_local
      cache%dadrji = cache%dadrji + dadrji_local
      cache%dadL = cache%dadL + dadL_local
      !$omp end critical (get_damat_3d_list_)
      deallocate(dadL_local, dadrij_local, dadrji_local, dadrdiag_local)
      !$omp end parallel

   end subroutine get_damat_3d_list

!> Real-space contribution to the derivative of the Coulomb matrix for the EEQBC model.
   subroutine get_damat_dir(rij, trans, capi, capj, rvdw, kbc, gam, dG, dS, dgam)
      !> Distance vector between two atoms (including lattice translation)
      real(wp), intent(in) :: rij(3)
      !> Direct lattice translation vectors (3 × N)
      real(wp), intent(in) :: trans(:, :)
      !> Bond capacitance of atom i
      real(wp), intent(in) :: capi
      !> Bond capacitance of atom j
      real(wp), intent(in) :: capj
      !> Van der Waals distance
      real(wp), intent(in) :: rvdw
      !> Bond capacitance exponent
      real(wp), intent(in) :: kbc
      !> Gaussian width parameter
      real(wp), intent(in) :: gam
      !> Derivative of Coulomb matrix element w.r.t. atomic position (3)
      real(wp), intent(out) :: dG(3)
      !> Derivative of Coulomb matrix element w.r.t. lattice parameters (3×3)
      real(wp), intent(out) :: dS(3, 3)
      !> Contribution to the derivative w.r.t. gam
      real(wp), intent(out) :: dgam

      integer :: itr
      real(wp) :: vec(3), r1, r2, gtmp, gam2, cmat

      dG(:) = 0.0_wp
      dS(:, :) = 0.0_wp
      dgam = 0.0_wp

      gam2 = gam * gam

      do itr = 1, size(trans, 2)
         vec(:) = rij(:) + trans(:, itr)
         r1 = norm2(vec)
         if (r1 < eps) cycle
         r2 = r1 * r1
         call get_cpair(kbc, cmat, r1, rvdw, capi, capj)
         gtmp = 2.0_wp * gam * exp(-r2 * gam2) / (sqrtpi * r2) - erf(r1 * gam) / (r2 * r1)
         dG(:) = dG - cmat * gtmp * vec
         dS(:, :) = dS - cmat * gtmp * spread(vec, 1, 3) * spread(vec, 2, 3)
         dgam = dgam + cmat * 2.0_wp * exp(-gam2 * r2) / sqrtpi
      end do

   end subroutine get_damat_dir

!> Contribution to the derivative of the Coulomb matrix from the derivative of the bond capacitance (direct part).
   subroutine get_damat_dc_dir(rij, trans, capi, capj, rvdw, kbc, gam, dG, dS)
      !> Distance vector between two atoms (including lattice translation)
      real(wp), intent(in) :: rij(3)
      !> Direct lattice translation vectors (3 × N)
      real(wp), intent(in) :: trans(:, :)
      !> Bond capacitance of atom i
      real(wp), intent(in) :: capi
      !> Bond capacitance of atom j
      real(wp), intent(in) :: capj
      !> Van der Waals distance
      real(wp), intent(in) :: rvdw
      !> Bond capacitance exponent
      real(wp), intent(in) :: kbc
      !> Gaussian width parameter
      real(wp), intent(in) :: gam
      !> Derivative of Coulomb matrix element w.r.t. atomic position (3)
      real(wp), intent(out) :: dG(3)
      !> Derivative of Coulomb matrix element w.r.t. lattice parameters (3×3)
      real(wp), intent(out) :: dS(3, 3)

      integer :: itr
      real(wp) :: vec(3), r1, gtmp(3), stmp(3, 3), tmp

      dG(:) = 0.0_wp
      dS(:, :) = 0.0_wp

      do itr = 1, size(trans, 2)
         vec(:) = rij(:) + trans(:, itr)
         r1 = norm2(vec)
         if (r1 < eps) cycle
         call get_dcpair(kbc, vec, rvdw, capi, capj, gtmp, stmp)
         tmp = erf(gam * r1) / r1
         dG(:) = dG(:) + tmp * gtmp
         dS(:, :) = dS(:, :) + tmp * stmp
      end do

   end subroutine get_damat_dc_dir

!> Build the bond capacitance matrix for a non‑periodic system.
   subroutine get_cmat_0d(self, mol, cmat)
      !> EEQBC model type
      class(eeqbc_model), intent(in) :: self
      !> Molecular structure data
      type(structure_type), intent(in) :: mol
      !> Output capacitance matrix (size ndim × ndim)
      real(wp), intent(out) :: cmat(:, :)

      integer :: iat, jat, izp, jzp
      real(wp) :: vec(3), rvdw, tmp, capi, capj, r1

      ! Thread-private array for reduction
      real(wp), allocatable :: cmat_local(:, :)

      cmat(:, :) = 0.0_wp

      !$omp parallel default(none) &
      !$omp shared(cmat, mol, self) &
      !$omp private(iat, izp, jat, jzp) &
      !$omp private(vec, r1, rvdw, tmp, capi, capj, cmat_local)
      allocate(cmat_local, source=cmat)
      !$omp do schedule(runtime)
      do iat = 1, mol%nat
         izp = mol%id(iat)
         capi = self%cap(izp)
         do jat = 1, iat - 1
            jzp = mol%id(jat)
            vec = mol%xyz(:, jat) - mol%xyz(:, iat)
            r1 = norm2(vec)
            rvdw = self%rvdw(izp, jzp)
            capj = self%cap(jzp)

            call get_cpair(self%kbc, tmp, r1, rvdw, capi, capj)

            ! Off-diagonal elements
            cmat_local(jat, iat) = -tmp
            cmat_local(iat, jat) = -tmp
            ! Diagonal elements
            cmat_local(iat, iat) = cmat_local(iat, iat) + tmp
            cmat_local(jat, jat) = cmat_local(jat, jat) + tmp
         end do
      end do
      !$omp end do
      !$omp critical (get_cmat_0d_)
      cmat(:, :) = cmat + cmat_local
      !$omp end critical (get_cmat_0d_)
      deallocate(cmat_local)
      !$omp end parallel

      if (size(cmat, 1) == mol%nat + 1) then
         cmat(mol%nat + 1, mol%nat + 1) = 1.0_wp
      end if

   end subroutine get_cmat_0d


   !> Build the bond capacitance matrix for a non‑periodic system.
   subroutine get_cmat_0d_list(self, mol, list, clist, cdiag)
      !> EEQBC model type
      class(eeqbc_model), intent(in) :: self
      !> Molecular structure data
      type(structure_type), intent(in) :: mol
      !> Multicharge neighbourlist type
      type(csr_list), intent(in) :: list
      !> Output capacitance matrix in compressed format, size of list%nlat
      real(wp), intent(out) :: clist(:)
      !> Output diagonal elements capacitance matrix in compressed format, size of mol%nat
      real(wp), intent(out) :: cdiag(:)

      integer :: iat, jat, kat, izp, jzp
      real(wp) :: vec(3), rvdw, tmp, capi, capj, r1

      ! Zero out global shared target arrays upfront
      clist(:) = 0.0_wp
      cdiag(:) = 0.0_wp

      !$omp parallel default(none) &
      !$omp shared(clist, mol, list, self) &
      !$omp private(iat, kat, izp, jat, jzp, vec, r1, rvdw, tmp, capi, capj) &
      !$omp reduction(+:cdiag)

      !$omp do schedule(runtime)
      do iat = 1, mol%nat
         izp = mol%id(iat)
         capi = self%cap(izp)
         do kat = list%inl(iat) + 1, list%inl(iat + 1) - 1
            jat = list%nlat(kat)
            jzp = mol%id(jat)
            vec = mol%xyz(:, jat) - mol%xyz(:, iat)
            r1 = norm2(vec)
            rvdw = self%rvdw(izp, jzp)
            capj = self%cap(jzp)

            call get_cpair(self%kbc, tmp, r1, rvdw, capi, capj)

            ! Safe direct write (kat is unique per thread)
            clist(kat) = -tmp

            ! Safe reduction writes (tracked in thread-local array copies)
            cdiag(iat) = cdiag(iat) + tmp
            cdiag(jat) = cdiag(jat) + tmp
         end do
      end do
      !$omp end do
      !$omp end parallel

   end subroutine get_cmat_0d_list

!> Build the bond capacitance matrix for a periodic system.
   subroutine get_cmat_3d(self, mol, wsc, cmat)
      !> EEQBC model type
      class(eeqbc_model), intent(in) :: self
      !> Molecular structure data
      type(structure_type), intent(in) :: mol
      !> Wigner–Seitz cell
      type(wignerseitz_cell_type), intent(in) :: wsc
      !> Output capacitance matrix (size ndim × ndim)
      real(wp), intent(out) :: cmat(:, :)

      integer :: iat, jat, izp, jzp, img
      real(wp) :: vec(3), rvdw, tmp, capi, capj, wsw
      real(wp), allocatable :: dtrans(:, :)

      ! Thread-private array for reduction
      real(wp), allocatable :: cmat_local(:, :)

      call get_dir_trans(mol, dtrans, cutoff)

      cmat(:, :) = 0.0_wp

      !$omp parallel default(none) &
      !$omp shared(cmat, mol, self, wsc, dtrans) &
      !$omp private(iat, izp, jat, jzp, img) &
      !$omp private(vec, rvdw, tmp, capi, capj, wsw, cmat_local)
      allocate(cmat_local, source=cmat)
      !$omp do schedule(runtime)
      do iat = 1, mol%nat
         izp = mol%id(iat)
         capi = self%cap(izp)
         do jat = 1, iat - 1
            jzp = mol%id(jat)
            rvdw = self%rvdw(izp, jzp)
            capj = self%cap(jzp)
            wsw = 1.0_wp / real(wsc%nimg(jat, iat), wp)
            do img = 1, wsc%nimg(jat, iat)
               vec = mol%xyz(:, iat) - mol%xyz(:, jat) - wsc%trans(:, wsc%tridx(img, jat, iat))

               call get_cpair_dir(self%kbc, vec, dtrans, rvdw, capi, capj, tmp)

               ! Off-diagonal elements
               cmat_local(jat, iat) = cmat_local(jat, iat) - tmp * wsw
               cmat_local(iat, jat) = cmat_local(iat, jat) - tmp * wsw
               ! Diagonal elements
               cmat_local(iat, iat) = cmat_local(iat, iat) + tmp * wsw
               cmat_local(jat, jat) = cmat_local(jat, jat) + tmp * wsw
            end do
         end do

         ! diagonal capacitance (interaction with images)
         rvdw = self%rvdw(izp, izp)
         wsw = 1.0_wp / real(wsc%nimg(iat, iat), wp)
         do img = 1, wsc%nimg(iat, iat)
            vec = wsc%trans(:, wsc%tridx(img, iat, iat))
            call get_cpair_dir(self%kbc, vec, dtrans, rvdw, capi, capi, tmp)
            cmat_local(iat, iat) = cmat_local(iat, iat) + tmp * wsw
         end do
      end do
      !$omp end do
      !$omp critical (get_cmat_3d_)
      cmat(:, :) = cmat + cmat_local
      !$omp end critical (get_cmat_3d_)
      deallocate(cmat_local)
      !$omp end parallel
      !
      if (size(cmat, 1) == mol%nat + 1) then
         cmat(mol%nat + 1, mol%nat + 1) = 1.0_wp
      end if

   end subroutine get_cmat_3d

!> Build the bond capacitance matrix for a periodic system using CSR adjacency list.
   subroutine get_cmat_3d_list(self, mol, list, clist, cdiag)
      !> EEQBC model type
      class(eeqbc_model), intent(in) :: self
      !> Molecular structure data
      type(structure_type), intent(in) :: mol
      !> Multicharge neighbourlist type (CSR format)
      type(csr_list), intent(in) :: list
      !> Output capacitance matrix in compressed format, size of list%nlat
      real(wp), intent(out) :: clist(:)
      !> Output diagonal elements capacitance matrix in compressed format, size of mol%nat
      real(wp), intent(out) :: cdiag(:)

      integer :: iat, jat, izp, jzp, img, kat
      real(wp) :: vec(3), rvdw, tmp, capi, capj, wsw, ctmp
      real(wp), allocatable :: dtrans(:, :)

      call get_dir_trans(mol, dtrans, cutoff)

      ! Zero out global shared target arrays upfront
      clist(:) = 0.0_wp
      cdiag(:) = 0.0_wp

      !$omp parallel default(none) &
      !$omp shared(clist, mol, list, self, dtrans) &
      !$omp private(iat, izp, jat, kat, jzp, img, vec, rvdw, tmp, capi, capj, wsw, ctmp) &
      !$omp reduction(+:cdiag)
      !$omp do schedule(runtime)
      do iat = 1, mol%nat
         izp = mol%id(iat)
         capi = self%cap(izp)

         ! 1. Off-diagonal neighbour pairs (jat /= iat)
         do kat = list%inl(iat) + 1, list%inl(iat+1) - 1
            if (list%wsc%nimg_list(kat) == 0) cycle

            jat = list%nlat(kat)
            jzp = mol%id(jat)
            capj = self%cap(jzp)
            rvdw = self%rvdw(izp, jzp)
            wsw = 1.0_wp / real(list%wsc%nimg_list(kat), wp)

            ! Initialize scalar accumulator for off-diagonal element
            ctmp = 0.0_wp

            ! Loop over translation images for pair kat
            do img = list%wsc%itr_list(kat), list%wsc%itr_list(kat+1) - 1
               vec = mol%xyz(:, iat) - mol%xyz(:, jat) - list%wsc%trans(:, list%wsc%tridx_list(img))

               call get_cpair_dir(self%kbc, vec, dtrans, rvdw, capi, capj, tmp)

               ! Accumulate off-diagonal elements in local scalar
               ctmp = ctmp - tmp * wsw

               ! Diagonal elements (Safe reduction writes)
               cdiag(iat) = cdiag(iat) + tmp * wsw
               cdiag(jat) = cdiag(jat) + tmp * wsw
            end do

            ! Safe direct write (kat is unique per thread)
            clist(kat) = ctmp

         end do

         ! 2. Self-interaction with periodic images R /= 0 (j = iat)
         if (list%wsc%nimg_list(list%inl(iat)) > 1) then
            rvdw = self%rvdw(izp, izp)
            wsw  = 1.0_wp / real(list%wsc%nimg_list(list%inl(iat)) - 1, wp)

            ! Skip identity image at list%wsc%itr_list(list%inl(iat)) and loop through R /= 0 images
            do img = list%wsc%itr_list(list%inl(iat)) + 1, list%wsc%itr_list(list%inl(iat) + 1) - 1
               vec = list%wsc%trans(:, list%wsc%tridx_list(img))

               call get_cpair_dir(self%kbc, vec, dtrans, rvdw, capi, capi, tmp)
               cdiag(iat) = cdiag(iat) + tmp * wsw
            end do
         end if

      end do
      !$omp end do
      !$omp end parallel

   end subroutine get_cmat_3d_list


!> Compute the bond capacitance between two atoms based on distance.
   subroutine get_cpair(kbc, cpair, r1, rvdw, capi, capj)
      !> Bond capacitance exponent
      real(wp), intent(in) :: kbc
      !> Distance between atoms
      real(wp), intent(in) :: r1
      !> Bond capacitance of atom i
      real(wp), intent(in) :: capi
      !> Bond capacitance of atom j
      real(wp), intent(in) :: capj
      !> Van der Waals distance
      real(wp), intent(in) :: rvdw
      !> Output bond capacitance value
      real(wp), intent(out) :: cpair

      real(wp) :: arg

      ! Capacitance of bond between atom i and j
      arg = -kbc * (r1 - rvdw) / rvdw
      cpair = sqrt(capi * capj) * 0.5_wp * (1.0_wp + erf(arg))
   end subroutine get_cpair

!> Sum the bond capacitance contributions over all lattice translations.
   subroutine get_cpair_dir(kbc, rij, trans, rvdw, capi, capj, cpair)
      !> Bond capacitance exponent
      real(wp), intent(in) :: kbc
      !> Distance vector between atoms
      real(wp), intent(in) :: rij(3)
      !> Direct lattice translation vectors (3 × N)
      real(wp), intent(in) :: trans(:, :)
      !> Van der Waals distance
      real(wp), intent(in) :: rvdw
      !> Bond capacitance of atom i
      real(wp), intent(in) :: capi
      !> Bond capacitance of atom j
      real(wp), intent(in) :: capj
      !> Output total bond capacitance
      real(wp), intent(out) :: cpair

      integer :: itr
      real(wp) :: vec(3), r1, tmp

      cpair = 0.0_wp
      do itr = 1, size(trans, 2)
         vec(:) = rij + trans(:, itr)
         r1 = norm2(vec)
         if (r1 < eps) cycle
         call get_cpair(kbc, tmp, r1, rvdw, capi, capj)
         cpair = cpair + tmp
      end do
   end subroutine get_cpair_dir

!> Compute the derivative of the bond capacitance with respect to atomic positions and lattice parameters.
   subroutine get_dcpair(kbc, vec, rvdw, capi, capj, dgpair, dspair)
      !> Bond capacitance exponent
      real(wp), intent(in) :: kbc
      !> Vector between atoms
      real(wp), intent(in) :: vec(3)
      !> Van der Waals distance
      real(wp), intent(in) :: rvdw
      !> Bond capacitance of atom i
      real(wp), intent(in) :: capi
      !> Bond capacitance of atom j
      real(wp), intent(in) :: capj
      !> Derivative w.r.t. atomic position (3)
      real(wp), intent(out) :: dgpair(3)
      !> Derivative w.r.t. lattice parameters (3��3)
      real(wp), intent(out) :: dspair(3, 3)

      real(wp) :: r1, arg, dtmp

      dgpair(:) = 0.0_wp
      dspair(:, :) = 0.0_wp

      r1 = norm2(vec)
      ! Capacitance of bond between atom i and j
      arg = -(kbc * (r1 - rvdw) / rvdw)**2
      dtmp = -sqrt(capi * capj) * kbc * exp(arg) / (sqrtpi * rvdw)
      dgpair = dtmp * vec / r1
      dspair = spread(dgpair, 1, 3) * spread(vec, 2, 3)
   end subroutine get_dcpair

!> Build the derivative of the bond capacitance matrix for a non‑periodic system.
   subroutine get_dcmat_0d(self, mol, dcdr, dcdL)
      !> EEQBC model type
      class(eeqbc_model), intent(in) :: self
      !> Molecular structure data
      type(structure_type), intent(in) :: mol
      !> Derivative of capacitance matrix w.r.t. atomic positions (3 × nat × ndim)
      real(wp), intent(out) :: dcdr(:, :, :)
      !> Derivative of capacitance matrix w.r.t. lattice parameters (3 × 3 �� ndim)
      real(wp), intent(out) :: dcdL(:, :, :)

      integer :: iat, jat, izp, jzp
      real(wp) :: vec(3), r2, rvdw, dtmp, arg, dG(3), dS(3, 3), capi, capj

      ! Thread-private arrays for reduction
      real(wp), allocatable :: dcdr_local(:, :, :), dcdL_local(:, :, :)

      dcdr(:, :, :) = 0.0_wp
      dcdL(:, :, :) = 0.0_wp

      allocate(dcdr_local, source=dcdr)
      allocate(dcdL_local, source=dcdL)

      !$omp parallel default(none) &
      !$omp shared(dcdr, dcdL, mol, self) &
      !$omp private(iat, izp, jat, jzp, r2, vec, rvdw) &
      !$omp private(dG, dS, dtmp, arg, capi, capj) &
      !$omp reduction(+:dcdr_local, dcdL_local)

      !$omp do schedule(runtime)
      do iat = 1, mol%nat
         izp = mol%id(iat)
         capi = self%cap(izp)
         do jat = 1, iat - 1
            jzp = mol%id(jat)
            capj = self%cap(jzp)
            rvdw = self%rvdw(izp, jzp)
            vec = mol%xyz(:, jat) - mol%xyz(:, iat)

            call get_dcpair(self%kbc, vec, rvdw, capi, capj, dG, dS)

            ! Off-diagonal elements
            dcdr_local(:, iat, jat) = +dG
            dcdr_local(:, jat, iat) = -dG
            ! Diagonal elements
            dcdr_local(:, iat, iat) = -dG + dcdr_local(:, iat, iat)
            dcdr_local(:, jat, jat) = +dG + dcdr_local(:, jat, jat)
            dcdL_local(:, :, iat) = +dS + dcdL_local(:, :, iat)
            dcdL_local(:, :, jat) = +dS + dcdL_local(:, :, jat)
         end do
      end do
      !$omp end do
      !$omp end parallel

      dcdr(:, :, :) = dcdr + dcdr_local
      dcdL(:, :, :) = dcdL + dcdL_local

      deallocate(dcdL_local, dcdr_local)

   end subroutine get_dcmat_0d

!> Build the derivative of the bond capacitance matrix for a non‑periodic system.
   subroutine get_dcmat_0d_list(self, mol, list, cache)
      !> EEQBC model type
      class(eeqbc_model), intent(in) :: self
      !> Molecular structure data
      type(structure_type), intent(in) :: mol
      !> Multicharge neighbourlist type
      type(csr_list), intent(in) :: list
      !> Multicharge cache
      type(mchrg_cache), intent(inout) :: cache

      integer :: iat, jat, kat, izp, jzp
      real(wp) :: vec(3), rvdw, dG(3), dS(3, 3), capi, capj
      real(wp), allocatable :: dcdrdiag_local(:, :), dcdL_local(:, :, :)

      cache%dcdrdiag(:, :) = 0.0_wp
      cache%dcdL(:, :, :) = 0.0_wp

      allocate(dcdrdiag_local, source=cache%dcdrdiag)
      allocate(dcdL_local, source=cache%dcdL)


      !$omp parallel default(none) &
      !$omp shared(cache, mol, list, self) &
      !$omp private(iat, izp, jat, kat, jzp, vec, rvdw) &
      !$omp private(dG, dS, capi, capj) &
      !$omp reduction(+:dcdrdiag_local, dcdL_local)

      !$omp do schedule(runtime)
      do iat = 1, mol%nat
         izp = mol%id(iat)
         capi = self%cap(izp)
         do kat = list%inl(iat) + 1, list%inl(iat + 1) - 1
            jat = list%nlat(kat)
            jzp = mol%id(jat)
            capj = self%cap(jzp)
            rvdw = self%rvdw(izp, jzp)
            vec = mol%xyz(:, jat) - mol%xyz(:, iat)

            call get_dcpair(self%kbc, vec, rvdw, capi, capj, dG, dS)

            ! Diagonal elements (Matches reference dcdr(i,i) and dcdr(j,j))
            dcdrdiag_local(:, iat) = dcdrdiag_local(:, iat) - dG
            dcdrdiag_local(:, jat) = dcdrdiag_local(:, jat) + dG

            ! Lattice derivatives
            dcdL_local(:, :, iat) = dcdL_local(:, :, iat) + dS
            dcdL_local(:, :, jat) = dcdL_local(:, :, jat) + dS
         end do
      end do
      !$omp end do
      !$omp end parallel

      cache%dcdrdiag(:, :) = cache%dcdrdiag + dcdrdiag_local
      cache%dcdL(:, :, :) = cache%dcdL + dcdL_local

      deallocate(dcdL_local, dcdrdiag_local)
   end subroutine get_dcmat_0d_list

!> Build the derivative of the bond capacitance matrix for a periodic system.
   subroutine get_dcmat_3d(self, mol, wsc, dcdr, dcdL)
      !> EEQBC model type
      class(eeqbc_model), intent(in) :: self
      !> Molecular structure data
      type(structure_type), intent(in) :: mol
      !> Wigner–Seitz cell
      type(wignerseitz_cell_type), intent(in) :: wsc
      !> Derivative of capacitance matrix w.r.t. atomic positions (3 × nat × ndim)
      real(wp), intent(out) :: dcdr(:, :, :)
      !> Derivative of capacitance matrix w.r.t. lattice parameters (3 × 3 × ndim)
      real(wp), intent(out) :: dcdL(:, :, :)

      integer :: iat, jat, izp, jzp, img
      real(wp) :: vec(3), r2, rvdw, dtmp, arg, dG(3), dS(3, 3), capi, capj, wsw
      real(wp), allocatable :: dtrans(:, :)

      ! Thread-private arrays for reduction
      real(wp), allocatable :: dcdr_local(:, :, :), dcdL_local(:, :, :)

      call get_dir_trans(mol, dtrans, cutoff)

      dcdr(:, :, :) = 0.0_wp
      dcdL(:, :, :) = 0.0_wp

      allocate(dcdr_local, source=dcdr)
      allocate(dcdL_local, source=dcdL)

      !$omp parallel default(none) &
      !$omp shared(dcdr, dcdL, mol, self, dtrans, wsc) &
      !$omp private(iat, izp, jat, jzp, r2, vec, rvdw) &
      !$omp private(dG, dS, dtmp, arg, capi, capj, wsw) &
      !$omp reduction(+:dcdr_local, dcdL_local)

      !$omp do schedule(runtime)
      do iat = 1, mol%nat
         izp = mol%id(iat)
         capi = self%cap(izp)
         do jat = 1, iat - 1
            jzp = mol%id(jat)
            capj = self%cap(jzp)
            rvdw = self%rvdw(izp, jzp)
            wsw = 1.0_wp / real(wsc%nimg(jat, iat), wp)
            do img = 1, wsc%nimg(jat, iat)
               vec = mol%xyz(:, jat) - mol%xyz(:, iat) + wsc%trans(:, wsc%tridx(img, jat, iat))

               call get_dcpair_dir(self%kbc, vec, dtrans, rvdw, capi, capj, dG, dS)

               ! Off-diagonal elements
               dcdr_local(:, iat, jat) = +dG * wsw + dcdr_local(:, iat, jat)
               dcdr_local(:, jat, iat) = -dG * wsw + dcdr_local(:, jat, iat)
               ! Diagonal elements
               dcdr_local(:, iat, iat) = -dG * wsw + dcdr_local(:, iat, iat)
               dcdr_local(:, jat, jat) = +dG * wsw + dcdr_local(:, jat, jat)
               dcdL_local(:, :, jat) = +dS * wsw + dcdL_local(:, :, jat)
               dcdL_local(:, :, iat) = +dS * wsw + dcdL_local(:, :, iat)
            end do
         end do

         rvdw = self%rvdw(izp, izp)
         wsw = 1.0_wp / real(wsc%nimg(iat, iat), wp)
         do img = 1, wsc%nimg(iat, iat)
            vec = wsc%trans(:, wsc%tridx(img, iat, iat))

            call get_dcpair_dir(self%kbc, vec, dtrans, rvdw, capi, capi, dG, dS)

            ! Positive diagonal elements
            dcdL_local(:, :, iat) = +dS * wsw + dcdL_local(:, :, iat)
         end do
      end do
      !$omp end do
      !$omp end parallel

      dcdr(:, :, :) = dcdr + dcdr_local
      dcdL(:, :, :) = dcdL + dcdL_local

      deallocate(dcdL_local, dcdr_local)

   end subroutine get_dcmat_3d

!> Build the derivative of the bond capacitance matrix for a periodic system.
   subroutine get_dcmat_3d_list(self, mol, list, cache)
      class(eeqbc_model), intent(in) :: self
      type(structure_type), intent(in) :: mol
      type(csr_list), intent(in) :: list
      type(mchrg_cache), intent(inout) :: cache

      integer :: iat, jat, izp, jzp, img, kat
      real(wp) :: vec(3), rvdw, dG(3), dS(3, 3), capi, capj, wsw
      real(wp), allocatable :: dtrans(:, :)
      real(wp), allocatable :: dcdrdiag_local(:, :), dcdL_local(:, :, :)

      call get_dir_trans(mol, dtrans, cutoff)

      cache%dcdrdiag(:, :) = 0.0_wp
      cache%dcdL(:, :, :) = 0.0_wp

      allocate(dcdrdiag_local, source=cache%dcdrdiag)
      allocate(dcdL_local, source=cache%dcdL)

      !$omp parallel default(none) &
      !$omp shared(cache, mol, list, self, dtrans) &
      !$omp private(iat, izp, jat, kat, jzp, vec, rvdw, dG, dS, capi, capj, wsw, img) &
      !$omp reduction(+:dcdrdiag_local, dcdL_local)

      !$omp do schedule(runtime)
      do iat = 1, mol%nat
         izp = mol%id(iat)
         capi = self%cap(izp)

         ! ------------------------------------------------------------------
         ! 1. Off-diagonal neighbor pairs (jat /= iat)
         ! ------------------------------------------------------------------
         do kat = list%inl(iat) + 1, list%inl(iat+1) - 1
            if (list%wsc%nimg_list(kat) == 0) cycle

            jat = list%nlat(kat)
            jzp = mol%id(jat)
            capj = self%cap(jzp)
            rvdw = self%rvdw(izp, jzp)
            wsw = 1.0_wp / real(list%wsc%nimg_list(kat), wp)

            do img = list%wsc%itr_list(kat), list%wsc%itr_list(kat+1) - 1
               vec = mol%xyz(:, jat) - mol%xyz(:, iat) + list%wsc%trans(:, list%wsc%tridx_list(img))
               call get_dcpair_dir(self%kbc, vec, dtrans, rvdw, capi, capj, dG, dS)

               ! Diagonal elements
               dcdrdiag_local(:, iat) = -dG * wsw + dcdrdiag_local(:, iat)
               dcdrdiag_local(:, jat) = +dG * wsw + dcdrdiag_local(:, jat)
               dcdL_local(:, :, jat) = +dS * wsw + dcdL_local(:, :, jat)
               dcdL_local(:, :, iat) = +dS * wsw + dcdL_local(:, :, iat)
            end do
         end do

         ! ------------------------------------------------------------------
         ! 2. Self-interaction with periodic images R /= 0 (j = iat)
         ! ------------------------------------------------------------------
         rvdw = self%rvdw(izp, izp)
         if (list%wsc%nimg_list(list%inl(iat)) > 1) then
            wsw = 1.0_wp / real(list%wsc%nimg_list(list%inl(iat)) - 1, wp)

            ! Skip identity image at list%wsc%itr_list(list%inl(iat)) and loop through R /= 0 images
            do img = list%wsc%itr_list(list%inl(iat)) + 1, list%wsc%itr_list(list%inl(iat) + 1) - 1
               vec = list%wsc%trans(:, list%wsc%tridx_list(img))
               call get_dcpair_dir(self%kbc, vec, dtrans, rvdw, capi, capi, dG, dS)
               dcdL_local(:, :, iat) = dcdL_local(:, :, iat) + dS * wsw
            end do
         end if

      end do
      !$omp end do
      !$omp end parallel

      cache%dcdrdiag(:, :) = cache%dcdrdiag + dcdrdiag_local
      cache%dcdL(:, :, :) = cache%dcdL + dcdL_local

      deallocate(dcdL_local, dcdrdiag_local)

   end subroutine get_dcmat_3d_list

!> Sum the derivative of bond capacitance over lattice translations.
   subroutine get_dcpair_dir(kbc, rij, trans, rvdw, capi, capj, dgpair, dspair)
      !> Bond capacitance exponent
      real(wp), intent(in) :: kbc
      !> Distance vector between atoms
      real(wp), intent(in) :: rij(3)
      !> Direct lattice translation vectors (3 × N)
      real(wp), intent(in) :: trans(:, :)
      !> Van der Waals distance
      real(wp), intent(in) :: rvdw
      !> Bond capacitance of atom i
      real(wp), intent(in) :: capi
      !> Bond capacitance of atom j
      real(wp), intent(in) :: capj
      !> Summed derivative w.r.t. atomic position (3)
      real(wp), intent(out) :: dgpair(3)
      !> Summed derivative w.r.t. lattice parameters (3×3)
      real(wp), intent(out) :: dspair(3, 3)

      integer :: itr
      real(wp) :: r1, arg, dtmp, dgtmp(3), dstmp(3, 3), vec(3)

      dgpair(:) = 0.0_wp
      dspair(:, :) = 0.0_wp
      do itr = 1, size(trans, 2)
         vec(:) = rij + trans(:, itr)
         r1 = norm2(vec)
         if (r1 < eps) cycle
         call get_dcpair(kbc, vec, rvdw, capi, capj, dgtmp, dstmp)
         dgpair(:) = dgpair + dgtmp
         dspair(:, :) = dspair + dstmp
      end do
   end subroutine get_dcpair_dir

   subroutine get_dcnpair(self, mol, iat, jat, rij, dG_ij, dG_ji)
      !> Coordination number container
      class(ncoord_type), intent(in) :: self
      !> Molecular structure data
      type(structure_type), intent(in) :: mol
      !> Indices of the interacting pair
      integer, intent(in) :: iat, jat
      !> Distance vector and scalar (rij = r_i - r_j - trans)
      real(wp), intent(in) :: rij(3)
      !> Derivatives: dG_ij = d(CN_j)/d(r_i) and dG_ji = d(CN_i)/d(r_j)
      real(wp), intent(out) :: dG_ij(3), dG_ji(3)

      ! Atomic numbers / species indices
      integer :: izp, jzp

      real(wp) :: den, countf, countd(3), r1, r2

      izp = mol%id(iat)
      jzp = mol%id(jat)

      r2 = sum(rij**2)
      r1 = sqrt(r2)

      den = self%get_en_factor(izp, jzp)
      countd = den * self%ncoord_dcount(izp, jzp, r1) * rij / r1


      ! Pay attention to the case when atoms are the same
      ! (e.g., self-interaction through periodic boundaries)
      if (iat == jat) then
         ! Avoid double counting for the same atom
         dG_ij(:) = 0.0_wp
         dG_ji(:) = 0.0_wp
      else

         ! dG_ij corresponds to the off-diagonal 'dcndrij'
         dG_ij(:) = countd * self%directed_factor

         ! dG_ji corresponds to the off-diagonal 'dcndrji'
         dG_ji(:) = -countd
      end if

   end subroutine get_dcnpair

   !> Sum the derivative of bond capacitance over lattice translations.
   subroutine get_dcnpair_dir(self, mol, trans, iat, jat, rij, dG_ij, dG_ji)
      !> Coordination number container
      class(ncoord_type), intent(in) :: self
      !> Molecular structure data
      type(structure_type), intent(in) :: mol
      !> Lattice points
      real(wp), intent(in) :: trans(:, :)
      !> Indices of the interacting pair
      integer, intent(in) :: iat, jat
      !> Distance vector and scalar (rij = r_i - r_j - trans)
      real(wp), intent(in) :: rij(3)
      !> Derivatives: dG_ij = d(CN_j)/d(r_i) and dG_ji = d(CN_i)/d(r_j)
      real(wp), intent(out) :: dG_ij(3), dG_ji(3)

      integer :: itr
      real(wp) :: r1, arg, dtmp, dgtmpij(3), dgtmpji(3), vec(3)

      dG_ij(:) = 0.0_wp
      dG_ji(:) = 0.0_wp
      do itr = 1, size(trans, 2)
         vec(:) = rij + trans(:, itr)
         r1 = norm2(vec)
         if (r1 < eps) cycle
         call get_dcnpair(self, mol, iat, jat, rij, dgtmpij, dgtmpji)
         dG_ij(:) = dG_ij + dgtmpij
         dG_ji(:) = dG_ji + dgtmpji

      end do
   end subroutine get_dcnpair_dir

   subroutine get_dcndiag_list(self, mol, trans, cn, dcndrdiag, dcndL, list)
      !> Coordination number container
      class(ncoord_type), intent(in) :: self
      !> Molecular structure data
      type(structure_type), intent(in) :: mol
      !> Lattice points
      real(wp), intent(in) :: trans(:, :)
      !> Error function coordination number.
      real(wp), intent(out) :: cn(:)
      !> Diagonal derivative of the CN with respect to the Cartesian coordinates.
      real(wp), intent(out) :: dcndrdiag(:, :)
      !> Derivative of the CN with respect to strain deformations.
      real(wp), intent(out) :: dcndL(:, :, :)
      !> Adjacency list for neighbourlist-based CN evaluation
      type(csr_list), intent(in) :: list

      integer :: iat, jat, kat, izp, jzp, itr
      real(wp) :: r2, r1, rij(3), countf, countd(3), sigma(3, 3), cutoff2, den

      ! Thread-private arrays for reduction
      real(wp), allocatable :: cn_local(:)
      real(wp), allocatable :: dcndrdiag_local(:, :),  dcndL_local(:, :, :)

      cn(:) = 0.0_wp
      dcndrdiag(:, :)  = 0.0_wp
      dcndL(:, :, :) = 0.0_wp
      cutoff2 = self%cutoff**2

      !$omp parallel default(none) &
      !$omp shared(self, mol, list, trans, cutoff2, cn, dcndrdiag, dcndL) &
      !$omp private(jat, kat, itr, izp, jzp, r2, rij, r1, den, countf, countd) &
      !$omp private(sigma, cn_local, dcndrdiag_local, dcndL_local)
      allocate(cn_local, source=cn)
      allocate(dcndrdiag_local, source=dcndrdiag)
      allocate(dcndL_local, source=dcndL)

      !$omp do schedule(runtime)
      do iat = 1, mol%nat
         izp = mol%id(iat)
         do kat = list%inl(iat), list%inl(iat+1) - 1
            jat = list%nlat(kat)
            jzp = mol%id(jat)
            den = self%get_en_factor(izp, jzp)

            do itr = 1, size(trans, dim=2)
               rij = mol%xyz(:, iat) - (mol%xyz(:, jat) + trans(:, itr))
               r2 = sum(rij**2)
               r1 = sqrt(r2)

               countf = den * self%ncoord_count(izp, jzp, r1)
               countd = den * self%ncoord_dcount(izp, jzp, r1) * rij/r1
               sigma = spread(countd, 1, 3) * spread(rij, 2, 3)

               ! Accumulate terms for the current atom (i)
               cn_local(iat) = cn_local(iat) + countf
               dcndrdiag_local(:,iat) = dcndrdiag_local(:,iat) + countd
               dcndL_local(:, :, iat) = dcndL_local(:, :, iat) + sigma

               ! Accumulate terms for the neighbour atom (j), avoiding double counting for self-images
               if (iat /= jat) then
                  cn_local(jat) = cn_local(jat) + countf * self%directed_factor
                  dcndrdiag_local(:,jat) = dcndrdiag_local(:,jat) - countd * self%directed_factor
                  dcndL_local(:, :, jat) = dcndL_local(:, :, jat) + sigma * self%directed_factor
               end if

            end do
         end do
      end do
      !$omp end do

      !$omp critical (ncoord_d_diag_list_)
      cn(:)            = cn(:)            + cn_local(:)
      dcndrdiag(:, :)  = dcndrdiag(:, :)  + dcndrdiag_local(:, :)
      dcndL(:, :, :)   = dcndL(:, :, :)   + dcndL_local(:, :, :)
      !$omp end critical (ncoord_d_diag_list_)

      deallocate(cn_local, dcndrdiag_local, dcndL_local)
      !$omp end parallel

   end subroutine get_dcndiag_list

   subroutine get_grad(self, mol, cache, p, gradient, sigma, alpha, beta, list)
      class(eeqbc_model), intent(in) :: self
      type(structure_type), intent(in) :: mol
      type(mchrg_cache), intent(in) :: cache
      real(wp), intent(in) :: p(:)
      real(wp), intent(inout) :: gradient(:, :)
      real(wp), intent(inout) :: sigma(:, :)
      real(wp), optional, intent(in) :: alpha
      real(wp), optional, intent(in) :: beta
      type(csr_list), intent(in), optional :: list

      if (.not. present(list)) then
         if (any(mol%periodic)) then
            call get_grad_3d(self, mol, cache, p, gradient, sigma, alphain=alpha, betain=beta)
         else
            call get_grad_0d(self, mol, cache, p, gradient, sigma, alphain=alpha, betain=beta)
         end if
      else
         if (any(mol%periodic)) then
            call get_grad_3d_list(self, mol, list, cache, p, gradient, sigma, alphain=alpha, betain=beta)
         else
            call get_grad_0d_list(self, mol, list, cache, p, gradient, sigma, alphain=alpha, betain=beta)
         end if
      end if

   end subroutine get_grad


   subroutine get_grad_0d_list(self, mol, list, cache, p, gradient, sigma, alphain, betain)
      class(eeqbc_model), intent(in) :: self
      type(structure_type), intent(in) :: mol
      type(csr_list), intent(in) :: list
      type(mchrg_cache), intent(in) :: cache
      real(wp), intent(in) :: p(:)
      real(wp), intent(inout) :: gradient(:, :)
      real(wp), intent(inout) :: sigma(:, :)
      real(wp), optional, intent(in) :: alphain
      real(wp), optional, intent(in) :: betain

      real(wp) :: alpha, beta
      integer :: iat, jat, kat, izp, jzp, start_kat, finish_kat
      real(wp) :: vec(3), r2, gam, arg, dtmp, norm_cn, factor
      real(wp) :: radi, radj, dradi, dradj, dG(3), dS(3, 3)
      real(wp) :: W_ii, W_jj, W_ij
      real(wp), allocatable :: gradient_local(:, :), sigma_local(:, :)
      real(wp), allocatable :: v(:), qlocacc(:), cnacc(:)
      real(wp), allocatable :: dtrans(:, :)

      alpha = 1.0_wp
      if (present(alphain)) alpha = alphain
      beta = 1.0_wp
      if (present(betain)) beta = betain

      allocate(dtrans, source=list%trans)
      allocate(qlocacc(mol%nat), cnacc(mol%nat))
      qlocacc = 0.0_wp
      cnacc = 0.0_wp

      allocate(v(mol%nat))
      call gemv_cmp(list, cache%clist, cache%cdiag, p, v, alpha=1.0_wp, beta=0.0_wp)

      allocate(gradient_local(3, mol%nat), source=0.0_wp)
      allocate(sigma_local(3, 3), source=0.0_wp)

      !$omp parallel do default(none) schedule(runtime) &
      !$omp shared(cache, mol, self, p, v, list, alpha, beta) &
      !$omp private(iat, jat, izp, jzp, gam, vec, r2, dtmp, norm_cn, arg) &
      !$omp private(radi, radj, dradi, dradj, dG, dS, W_ii, W_jj, W_ij) &
      !$omp reduction(+:cnacc, qlocacc, gradient_local, sigma_local)

      do iat = 1, mol%nat
         izp = mol%id(iat)
         norm_cn = 1.0_wp / self%avg_cn(izp)**self%norm_exp
         radi = self%rad(izp) * (1.0_wp - self%kcnrad * cache%cn(iat) * norm_cn)
         dradi = -self%rad(izp) * self%kcnrad * norm_cn

         W_ii = p(iat) * cache%vrhs(iat)

         do kat = list%inl(iat) + 1, list%inl(iat + 1) - 1
            jat = list%nlat(kat)
            jzp = mol%id(jat)
            vec = mol%xyz(:, jat) - mol%xyz(:, iat)
            r2 = dot_product(vec, vec)

            W_jj = p(jat) * cache%vrhs(jat)
            W_ij = p(iat) * cache%vrhs(jat) + p(jat) * cache%vrhs(iat)

            norm_cn = 1.0_wp / self%avg_cn(jzp)**self%norm_exp
            radj = self%rad(jzp) * (1.0_wp - self%kcnrad * cache%cn(jat) * norm_cn)
            dradj = -self%rad(jzp) * self%kcnrad * norm_cn

            gam = 1.0_wp / sqrt(radi**2 + radj**2)
            arg = gam * gam * r2
            dtmp = 2.0_wp * exp(-arg) / (sqrtpi)

            ! Accumulate scalar weights for CN derivatives
            cnacc(iat) = cnacc(iat) - (dtmp * radi * dradi * gam**3.0_wp) * cache%clist(kat) * W_ij * alpha
            cnacc(jat) = cnacc(jat) - (dtmp * radj * dradj * gam**3.0_wp) * cache%clist(kat) * W_ij * alpha

            ! 1. Explicit Geometry Derivative
            dtmp = 2.0_wp * gam * exp(-arg) / (sqrtpi * r2) - erf(sqrt(arg)) / (r2 * sqrt(r2))
            dG = dtmp * vec
            dS = spread(dG, 1, 3) * spread(vec, 2, 3)

            gradient_local(:, iat) = gradient_local(:, iat) - dG * cache%clist(kat) * W_ij * alpha
            gradient_local(:, jat) = gradient_local(:, jat) + dG * cache%clist(kat) * W_ij * alpha
            sigma_local(:, :)   = sigma_local(:, :)   + dS * cache%clist(kat) * W_ij * alpha

            ! 3 & 4. Capacitance derivatives
            call get_dcpair(self%kbc, vec, self%rvdw(izp, jzp), self%cap(izp), self%cap(jzp), dG, dS)
            dtmp = erf(sqrt(r2) * gam) / sqrt(r2)
            gradient_local(:, iat) = gradient_local(:, iat) + dtmp * dG * W_ij * alpha &
            & + p(iat) * cache%xtmp(jat) * dG * beta &
            & - p(jat) * (cache%xtmp(jat) - cache%xtmp(iat)) * dG * beta

            gradient_local(:, jat) = gradient_local(:, jat) - dtmp * dG * W_ij * alpha &
            & - p(jat) * cache%xtmp(iat) * dG * beta &
            & + p(iat) * (cache%xtmp(iat) - cache%xtmp(jat)) * dG * beta

            sigma_local(:, :)   = sigma_local(:, :)   - dtmp * W_ij * dS * alpha &
            & - p(iat) * cache%xtmp(jat) * dS * beta &
            & - p(jat) * cache%xtmp(iat) * dS * beta

            dtmp = (self%eta(jzp) + self%kqeta(jzp) * cache%qloc(jat) + sqrt2pi / radj)
            gradient_local(:, iat) = gradient_local(:, iat) - dtmp * dG * W_jj * alpha
            dtmp = (self%eta(izp) + self%kqeta(izp) * cache%qloc(iat) + sqrt2pi / radi)
            gradient_local(:, jat) = gradient_local(:, jat) + dtmp * dG * W_ii * alpha
         end do

         ! 5. Diagonal weight accumulation
         qlocacc(iat) = qlocacc(iat) + self%kqeta(izp) * W_ii * cache%cdiag(iat) * alpha &
         & + v(iat) * self%kqchi(mol%id(iat)) * beta
         cnacc(iat) = cnacc(iat) - sqrt2pi * dradi / (radi**2) * W_ii * cache%cdiag(iat) * alpha &
         & + v(iat) * self%kcnchi(mol%id(iat)) * beta

         ! 6. Intrinsic capacitance
         dtmp = (self%eta(izp) + self%kqeta(izp) * cache%qloc(iat) + sqrt2pi / radi) * W_ii
         gradient_local(:, iat) = gradient_local(:, iat) + dtmp * cache%dcdrdiag(:, iat) * alpha &
         & + p(iat) * cache%xtmp(iat) * cache%dcdrdiag(:, iat) * beta

         sigma_local(:, :)   = sigma_local(:, :)   + dtmp * cache%dcdL(:, :, iat) * alpha &
         & + p(iat) * cache%xtmp(iat) * cache%dcdL(:, :, iat) * beta
      end do
      !$omp end parallel do

      gradient = gradient + gradient_local
      sigma = sigma + sigma_local

      gradient_local = 0.0_wp
      sigma_local = 0.0_wp

      call self%ncoord%add_coordination_number_derivs_list(mol, dtrans, cnacc, gradient_local, sigma_local, list)
      call self%ncoord_en%add_coordination_number_derivs_list(mol, dtrans, qlocacc, gradient_local, sigma_local, list)

      gradient = gradient + gradient_local
      sigma = sigma + sigma_local

      deallocate(gradient_local, sigma_local, qlocacc, cnacc, dtrans, v)
   end subroutine get_grad_0d_list

   subroutine get_grad_0d(self, mol, cache, p, gradient, sigma, alphain, betain)
      class(eeqbc_model), intent(in) :: self
      type(structure_type), intent(in) :: mol
      type(mchrg_cache), intent(in) :: cache
      real(wp), intent(in) :: p(:)
      real(wp), intent(inout) :: gradient(:, :)
      real(wp), intent(inout) :: sigma(:, :)
      real(wp), optional, intent(in) :: alphain
      real(wp), optional, intent(in) :: betain

      real(wp) :: alpha, beta
      integer :: iat, jat, kat, izp, jzp, start_kat, finish_kat
      real(wp) :: vec(3), r2, gam, arg, dtmp, norm_cn, factor
      real(wp) :: radi, radj, dradi, dradj, dG(3), dS(3, 3)
      real(wp) :: W_ii, W_jj, W_ij
      real(wp), allocatable :: gradient_local(:, :), sigma_local(:, :)
      real(wp), allocatable :: v(:), qlocacc(:), cnacc(:)
      real(wp), allocatable :: dtrans(:, :)

      alpha = 1.0_wp
      if (present(alphain)) alpha = alphain
      beta = 1.0_wp
      if (present(betain)) beta = betain

      allocate(dtrans, source=cache%trans)
      allocate(qlocacc(mol%nat), cnacc(mol%nat))
      qlocacc = 0.0_wp
      cnacc = 0.0_wp

      allocate(v(mol%nat))
      call symv(cache%cmat(:mol%nat, :mol%nat), p, v, alpha=1.0_wp, beta=0.0_wp, uplo='l')

      allocate(gradient_local(3, mol%nat), source=0.0_wp)
      allocate(sigma_local(3, 3), source=0.0_wp)

      !$omp parallel do default(none) schedule(runtime) &
      !$omp shared(cache, mol, self, p, v, alpha, beta) &
      !$omp private(iat, jat, izp, jzp, gam, vec, r2, dtmp, norm_cn, arg) &
      !$omp private(radi, radj, dradi, dradj, dG, dS, W_ii, W_jj, W_ij) &
      !$omp reduction(+:cnacc, qlocacc, gradient_local, sigma_local)

      do iat = 1, mol%nat
         izp = mol%id(iat)
         norm_cn = 1.0_wp / self%avg_cn(izp)**self%norm_exp
         radi = self%rad(izp) * (1.0_wp - self%kcnrad * cache%cn(iat) * norm_cn)
         dradi = -self%rad(izp) * self%kcnrad * norm_cn

         W_ii = p(iat) * cache%vrhs(iat)

         do jat = 1, iat - 1
            jzp = mol%id(jat)
            vec = mol%xyz(:, jat) - mol%xyz(:, iat)
            r2 = dot_product(vec, vec)

            W_jj = p(jat) * cache%vrhs(jat)
            W_ij = p(iat) * cache%vrhs(jat) + p(jat) * cache%vrhs(iat)

            norm_cn = 1.0_wp / self%avg_cn(jzp)**self%norm_exp
            radj = self%rad(jzp) * (1.0_wp - self%kcnrad * cache%cn(jat) * norm_cn)
            dradj = -self%rad(jzp) * self%kcnrad * norm_cn

            gam = 1.0_wp / sqrt(radi**2 + radj**2)
            arg = gam * gam * r2
            dtmp = 2.0_wp * exp(-arg) / (sqrtpi)

            ! Accumulate scalar weights for CN derivatives
            cnacc(iat) = cnacc(iat) - (dtmp * radi * dradi * gam**3.0_wp) * cache%cmat(iat, jat) * W_ij * alpha
            cnacc(jat) = cnacc(jat) - (dtmp * radj * dradj * gam**3.0_wp) * cache%cmat(iat, jat) * W_ij * alpha

            ! 1. Explicit Geometry Derivative
            dtmp = 2.0_wp * gam * exp(-arg) / (sqrtpi * r2) - erf(sqrt(arg)) / (r2 * sqrt(r2))
            dG = dtmp * vec
            dS = spread(dG, 1, 3) * spread(vec, 2, 3)

            gradient_local(:, iat) = gradient_local(:, iat) - dG * cache%cmat(iat, jat) * W_ij * alpha
            gradient_local(:, jat) = gradient_local(:, jat) + dG * cache%cmat(iat, jat) * W_ij * alpha
            sigma_local(:, :)   = sigma_local(:, :)   + dS * cache%cmat(iat, jat) * W_ij * alpha

            ! 3 & 4. Capacitance derivatives
            call get_dcpair(self%kbc, vec, self%rvdw(izp, jzp), self%cap(izp), self%cap(jzp), dG, dS)
            dtmp = erf(sqrt(r2) * gam) / sqrt(r2)
            gradient_local(:, iat) = gradient_local(:, iat) + dtmp * dG * W_ij * alpha &
            & + p(iat) * cache%xtmp(jat) * dG * beta &
            & - p(jat) * (cache%xtmp(jat) - cache%xtmp(iat)) * dG * beta

            gradient_local(:, jat) = gradient_local(:, jat) - dtmp * dG * W_ij * alpha &
            & - p(jat) * cache%xtmp(iat) * dG * beta &
            & + p(iat) * (cache%xtmp(iat) - cache%xtmp(jat)) * dG * beta

            sigma_local(:, :)   = sigma_local(:, :)   - dtmp * W_ij * dS * alpha &
            & - p(iat) * cache%xtmp(jat) * dS * beta &
            & - p(jat) * cache%xtmp(iat) * dS * beta

            dtmp = (self%eta(jzp) + self%kqeta(jzp) * cache%qloc(jat) + sqrt2pi / radj)
            gradient_local(:, iat) = gradient_local(:, iat) - dtmp * dG * W_jj * alpha
            dtmp = (self%eta(izp) + self%kqeta(izp) * cache%qloc(iat) + sqrt2pi / radi)
            gradient_local(:, jat) = gradient_local(:, jat) + dtmp * dG * W_ii * alpha
         end do

         ! 5. Diagonal weight accumulation
         qlocacc(iat) = qlocacc(iat) + self%kqeta(izp) * W_ii * cache%cmat(iat, iat) * alpha &
         & + v(iat) * self%kqchi(mol%id(iat)) * beta
         cnacc(iat) = cnacc(iat) - sqrt2pi * dradi / (radi**2) * W_ii * cache%cmat(iat, iat) * alpha &
         & + v(iat) * self%kcnchi(mol%id(iat)) * beta

         ! 6. Intrinsic capacitance
         dtmp = (self%eta(izp) + self%kqeta(izp) * cache%qloc(iat) + sqrt2pi / radi) * W_ii
         gradient_local(:, iat) = gradient_local(:, iat) + dtmp * cache%dcdr(:, iat, iat) * alpha &
         & + p(iat) * cache%xtmp(iat) * cache%dcdr(:, iat, iat) * beta

         sigma_local(:, :)   = sigma_local(:, :)   + dtmp * cache%dcdL(:, :, iat) * alpha &
         & + p(iat) * cache%xtmp(iat) * cache%dcdL(:, :, iat) * beta
      end do
      !$omp end parallel do

      gradient = gradient + gradient_local
      sigma = sigma + sigma_local

      gradient_local = 0.0_wp
      sigma_local = 0.0_wp

      call self%ncoord%add_coordination_number_derivs(mol, dtrans, cnacc, gradient_local, sigma_local)
      call self%ncoord_en%add_coordination_number_derivs(mol, dtrans, qlocacc, gradient_local, sigma_local)

      gradient = gradient + gradient_local
      sigma = sigma + sigma_local

      deallocate(gradient_local, sigma_local, qlocacc, cnacc, dtrans, v)
   end subroutine get_grad_0d


   subroutine get_grad_3d_list(self, mol, list, cache, p, gradient, sigma, alphain, betain)
      class(eeqbc_model), intent(in) :: self
      type(structure_type), intent(in) :: mol
      type(csr_list), intent(in) :: list
      type(mchrg_cache), intent(in) :: cache
      real(wp), intent(in) :: p(:)
      real(wp), intent(inout) :: gradient(:, :)
      real(wp), intent(inout) :: sigma(:, :)
      real(wp), optional, intent(in) :: alphain
      real(wp), optional, intent(in) :: betain

      real(wp) :: alpha, beta
      integer :: iat, jat, kat, izp, jzp, img
      real(wp) :: vec(3), r2, gam, dtmp, capi, capj, dgam, rvdw, wsw, factor
      real(wp) :: radi, radj, dradi, dradj, dG(3), dS(3, 3), rij(3), ctmp
      real(wp) :: W_ii, W_jj, W_ij, norm_cn
      real(wp), allocatable :: gradient_local(:, :), sigma_local(:, :)
      real(wp), allocatable :: v(:), qlocacc(:), cnacc(:), dtrans(:, :)

      alpha = 1.0_wp
      if (present(alphain)) alpha = alphain
      beta = 1.0_wp
      if (present(betain)) beta = betain

      call get_dir_trans(mol, dtrans, cutoff)
      allocate(qlocacc(mol%nat), cnacc(mol%nat))
      qlocacc = 0.0_wp
      cnacc = 0.0_wp

      ! Allocate reduction targets
      allocate(gradient_local(3, mol%nat), source=0.0_wp)
      allocate(sigma_local(3, 3), source=0.0_wp)

      allocate(v(mol%nat))
      call gemv_cmp(list, cache%clist, cache%cdiag, p, v, alpha=1.0_wp, beta=0.0_wp)

      !$omp parallel do default(none) schedule(runtime) &
      !$omp shared(cache, mol, list, self, p, v, dtrans, alpha, beta) &
      !$omp private(iat, kat, izp, jat, jzp, gam, vec, r2, dtmp, radi, radj, dradi, dradj, dG, dS, norm_cn) &
      !$omp private(W_ii, W_jj, W_ij, wsw, capi, capj, rvdw, img, dgam, ctmp) &
      !$omp reduction(+:cnacc, qlocacc, gradient_local, sigma_local)
      do iat = 1, mol%nat
         izp = mol%id(iat)
         norm_cn = 1.0_wp / self%avg_cn(izp)**self%norm_exp
         radi = self%rad(izp) * (1.0_wp - self%kcnrad * cache%cn(iat) * norm_cn)
         dradi = -self%rad(izp) * self%kcnrad * norm_cn
         capi = self%cap(izp)
         W_ii = p(iat) * cache%vrhs(iat)

         do kat = list%inl(iat) + 1, list%inl(iat + 1) - 1
            jat = list%nlat(kat)
            jzp = mol%id(jat)
            capj = self%cap(jzp)
            rvdw = self%rvdw(izp, jzp)
            if (list%wsc%nimg_list(kat) == 0) cycle
            wsw = 1.0_wp / real(list%wsc%nimg_list(kat), wp)
            W_jj = p(jat) * cache%vrhs(jat)
            W_ij = p(iat) * cache%vrhs(jat) + p(jat) * cache%vrhs(iat)

            norm_cn = 1.0_wp / self%avg_cn(jzp)**self%norm_exp
            radj = self%rad(jzp) * (1.0_wp - self%kcnrad * cache%cn(jat) * norm_cn)
            dradj = -self%rad(jzp) * self%kcnrad * norm_cn
            gam = 1.0_wp / sqrt(radi**2 + radj**2)

            do img = list%wsc%itr_list(kat), list%wsc%itr_list(kat+1) - 1
               vec = mol%xyz(:, jat) - mol%xyz(:, iat) + list%wsc%trans(:, list%wsc%tridx_list(img))
               call get_damat_dir(vec, dtrans, capi, capj, rvdw, self%kbc, gam, dG, dS, dgam)

               cnacc(iat) = cnacc(iat) - (dgam * wsw * radi * dradi * gam**3.0_wp) * W_ij * alpha
               cnacc(jat) = cnacc(jat) - (dgam * wsw * radj * dradj * gam**3.0_wp) * W_ij * alpha

               gradient_local(:, iat) = gradient_local(:, iat) - dG * wsw * W_ij * alpha
               gradient_local(:, jat) = gradient_local(:, jat) + dG * wsw * W_ij * alpha
               sigma_local(:, :) = sigma_local(:, :) + dS * wsw * W_ij * alpha

               call get_damat_dc_dir(vec, dtrans, capi, capj, rvdw, self%kbc, gam, dG, dS)
               gradient_local(:, iat) = gradient_local(:, iat) + dG * wsw * W_ij * alpha
               gradient_local(:, jat) = gradient_local(:, jat) - dG * wsw * W_ij * alpha
               sigma_local(:, :)   = sigma_local(:, :) - W_ij * dS * wsw * alpha

               call get_dcpair_dir(self%kbc, vec, dtrans, rvdw, capi, capj, dG, dS)
               dtmp = (self%eta(jzp) + self%kqeta(jzp) * cache%qloc(jat) + sqrt2pi / radj)
               gradient_local(:, iat) = gradient_local(:, iat) - dtmp * dG * wsw * W_jj * alpha &
               & + p(iat) * cache%xtmp(jat) * dG * wsw * beta &
               & - p(jat) * (cache%xtmp(jat) - cache%xtmp(iat)) * dG * wsw * beta

               dtmp = (self%eta(izp) + self%kqeta(izp) * cache%qloc(iat) + sqrt2pi / radi)
               gradient_local(:, jat) = gradient_local(:, jat) + dtmp * dG * wsw * W_ii * alpha &
               & - p(jat) * cache%xtmp(iat) * dG * wsw * beta &
               & + p(iat) * (cache%xtmp(iat) - cache%xtmp(jat)) * dG * wsw * beta

               ! Lattice Sigma Updates
               ! This ensures the (xi - xj) term is correctly formed against the diagonal dcdL
               sigma_local(:, :) = sigma_local(:, :) - p(iat) * cache%xtmp(jat) * dS * wsw * beta
               sigma_local(:, :) = sigma_local(:, :) - p(jat) * cache%xtmp(iat) * dS * wsw * beta

            end do
         end do

         ! Diagonal corrections
         gam = 1.0_wp / sqrt(2.0_wp * radi**2)
         rvdw = self%rvdw(izp, izp)

         if (list%wsc%nimg_list(list%inl(iat)) > 1) then
            wsw = 1.0_wp / real(list%wsc%nimg_list(list%inl(iat)) - 1, wp)

            do img = list%wsc%itr_list(list%inl(iat)) + 1, list%wsc%itr_list(list%inl(iat) + 1) - 1
               vec = list%wsc%trans(:, list%wsc%tridx_list(img))
               call get_damat_dir(vec, dtrans, capi, capi, rvdw, self%kbc, gam, dG, dS, dgam)
               sigma_local(:, :) = sigma_local(:, :) + dS * wsw * W_ii * alpha
               cnacc(iat) = cnacc(iat) - (dgam * wsw * radi * dradi * gam**3.0_wp) * W_ii * alpha

               call get_damat_dc_dir(vec, dtrans, capi, capi, rvdw, self%kbc, gam, dG, dS)
               sigma_local(:, :) = sigma_local(:, :) - W_ii * dS * wsw * alpha

               call get_dcpair_dir(self%kbc, vec, dtrans, rvdw, capi, capi, dG, dS)
               sigma_local(:, :) = sigma_local(:, :) - p(iat) * cache%xtmp(iat) * dS * wsw * beta

               call get_cpair_dir(self%kbc, vec, dtrans, rvdw, capi, capi, ctmp)
               ctmp = ctmp * wsw

               v(iat) = v(iat) - ctmp * p(iat)
            end do
         end if

         qlocacc(iat) = qlocacc(iat) + self%kqeta(izp) * W_ii * cache%cdiag(iat) * alpha &
         & + v(iat) * self%kqchi(izp) * beta
         cnacc(iat) = cnacc(iat) - sqrt2pi * dradi / (radi**2) * W_ii * cache%cdiag(iat) * alpha &
         & + v(iat) * self%kcnchi(izp) * beta

         dtmp = (self%eta(izp) + self%kqeta(izp) * cache%qloc(iat) + sqrt2pi / radi) * W_ii
         gradient_local(:, iat) = gradient_local(:, iat) + dtmp * cache%dcdrdiag(:, iat) * alpha &
         & + p(iat) * cache%xtmp(iat) * cache%dcdrdiag(:, iat) * beta
         sigma_local(:, :) = sigma_local(:, :) + dtmp * cache%dcdL(:, :, iat) * alpha &
         & + p(iat) * cache%xtmp(iat) * cache%dcdL(:, :, iat) * beta
      end do
      !$omp end parallel do

      gradient = gradient + gradient_local
      sigma = sigma + sigma_local

      gradient_local = 0.0_wp
      sigma_local = 0.0_wp

      call self%ncoord%add_coordination_number_derivs_list(mol, dtrans, cnacc, gradient_local, sigma_local, list)
      call self%ncoord_en%add_coordination_number_derivs_list(mol, dtrans, qlocacc, gradient_local, sigma_local, list)

      gradient = gradient + gradient_local
      sigma = sigma + sigma_local

      deallocate(gradient_local, sigma_local, qlocacc, cnacc, dtrans, v)
   end subroutine get_grad_3d_list

   subroutine get_grad_3d(self, mol, cache, p, gradient, sigma, alphain, betain)
      class(eeqbc_model), intent(in) :: self
      type(structure_type), intent(in) :: mol
      type(mchrg_cache), intent(in) :: cache
      real(wp), intent(in) :: p(:)
      real(wp), intent(inout) :: gradient(:, :)
      real(wp), intent(inout) :: sigma(:, :)
      real(wp), optional, intent(in) :: alphain
      real(wp), optional, intent(in) :: betain

      real(wp) :: alpha, beta
      integer :: iat, jat, izp, jzp, img
      real(wp) :: vec(3), gam, dtmp, capi, capj, dgam, rvdw, wsw, ctmp
      real(wp) :: radi, radj, dradi, dradj, dG(3), dS(3, 3), norm_cn
      real(wp) :: W_ii, W_jj, W_ij
      real(wp), allocatable :: gradient_local(:, :), sigma_local(:, :)
      real(wp), allocatable :: qlocacc(:), cnacc(:), dtrans(:, :), v(:)

      alpha = 1.0_wp
      if (present(alphain)) alpha = alphain
      beta = 1.0_wp
      if (present(betain)) beta = betain

      call get_dir_trans(mol, dtrans, cutoff)

      allocate(qlocacc(mol%nat), cnacc(mol%nat))
      qlocacc = 0.0_wp
      cnacc = 0.0_wp

      allocate(v(mol%nat))
      call symv(cache%cmat(:mol%nat, :mol%nat), p, v, alpha=1.0_wp, beta=0.0_wp, uplo='l')

      allocate(gradient_local(3, mol%nat), source=0.0_wp)
      allocate(sigma_local(3, 3), source=0.0_wp)

      !$omp parallel do default(none) schedule(runtime) &
      !$omp shared(cache, mol, self, p, dtrans, v, alpha, beta) &
      !$omp private(iat, izp, jat, jzp, gam, vec, dtmp, radi, radj, dradi, dradj, dG, dS, norm_cn) &
      !$omp private(W_ii, W_jj, W_ij, wsw, capi, capj, rvdw, img, dgam, ctmp) &
      !$omp reduction(+:cnacc, qlocacc, gradient_local, sigma_local)
      do iat = 1, mol%nat
         izp = mol%id(iat)
         norm_cn = 1.0_wp / self%avg_cn(izp)**self%norm_exp
         radi = self%rad(izp) * (1.0_wp - self%kcnrad * cache%cn(iat) * norm_cn)
         dradi = -self%rad(izp) * self%kcnrad * norm_cn
         capi = self%cap(izp)
         W_ii = p(iat) * cache%vrhs(iat)

         do jat = 1, iat - 1
            jzp = mol%id(jat)
            capj = self%cap(jzp)
            rvdw = self%rvdw(izp, jzp)
            wsw = 1.0_wp / real(cache%wsc%nimg(iat, jat), wp)
            W_jj = p(jat) * cache%vrhs(jat)
            W_ij = p(iat) * cache%vrhs(jat) + p(jat) * cache%vrhs(iat)

            norm_cn = 1.0_wp / self%avg_cn(jzp)**self%norm_exp
            radj = self%rad(jzp) * (1.0_wp - self%kcnrad * cache%cn(jat) * norm_cn)
            dradj = -self%rad(jzp) * self%kcnrad * norm_cn
            gam = 1.0_wp / sqrt(radi**2 + radj**2)

            do img = 1, cache%wsc%nimg(iat, jat)
               vec = mol%xyz(:, jat) - mol%xyz(:, iat) + cache%wsc%trans(:, cache%wsc%tridx(img, jat, iat))

               ! 1. A-matrix explicit kernel derivative
               call get_damat_dir(vec, dtrans, capi, capj, rvdw, self%kbc, gam, dG, dS, dgam)

               cnacc(iat) = cnacc(iat) - (dgam * wsw * radi * dradi * gam**3.0_wp) * W_ij * alpha
               cnacc(jat) = cnacc(jat) - (dgam * wsw * radj * dradj * gam**3.0_wp) * W_ij * alpha

               gradient_local(:, iat) = gradient_local(:, iat) - dG * wsw * W_ij * alpha
               gradient_local(:, jat) = gradient_local(:, jat) + dG * wsw * W_ij * alpha
               sigma_local(:, :)   = sigma_local(:, :)   + dS * wsw * W_ij * alpha

               ! 2. A-matrix capacitance terms
               call get_damat_dc_dir(vec, dtrans, capi, capj, rvdw, self%kbc, gam, dG, dS)
               gradient_local(:, iat) = gradient_local(:, iat) + dG * wsw * W_ij * alpha
               gradient_local(:, jat) = gradient_local(:, jat) - dG * wsw * W_ij * alpha
               sigma_local(:, :)   = sigma_local(:, :)   - W_ij * dS * wsw * alpha

               ! 3. Pair capacitance (b-vector) and hardness terms
               call get_dcpair_dir(self%kbc, vec, dtrans, rvdw, capi, capj, dG, dS)

               dtmp = (self%eta(jzp) + self%kqeta(jzp) * cache%qloc(jat) + sqrt2pi / radj)
               gradient_local(:, iat) = gradient_local(:, iat) - dtmp * dG * wsw * W_jj * alpha &
               & + p(iat) * cache%xtmp(jat) * dG * wsw * beta &
               & - p(jat) * (cache%xtmp(jat) - cache%xtmp(iat)) * dG * wsw * beta

               dtmp = (self%eta(izp) + self%kqeta(izp) * cache%qloc(iat) + sqrt2pi / radi)
               gradient_local(:, jat) = gradient_local(:, jat) + dtmp * dG * wsw * W_ii * alpha &
               & - p(jat) * cache%xtmp(iat) * dG * wsw * beta &
               & + p(iat) * (cache%xtmp(iat) - cache%xtmp(jat)) * dG * wsw * beta

               sigma_local(:, :)   = sigma_local(:, :)   - p(iat) * cache%xtmp(jat) * dS * wsw * beta &
               & - p(jat) * cache%xtmp(iat) * dS * wsw * beta
            end do
         end do

         ! Diagonal corrections & Self-images
         gam = 1.0_wp / sqrt(2.0_wp * radi**2)
         rvdw = self%rvdw(izp, izp)
         wsw = 1.0_wp / real(cache%wsc%nimg(iat, iat), wp)

         do img = 1, cache%wsc%nimg(iat, iat)
            vec = cache%wsc%trans(:, cache%wsc%tridx(img, iat, iat))

            call get_damat_dir(vec, dtrans, capi, capi, rvdw, self%kbc, gam, dG, dS, dgam)
            sigma_local(:, :) = sigma_local(:, :) + dS * wsw * W_ii * alpha
            cnacc(iat) = cnacc(iat) - (dgam * wsw * radi * dradi * gam**3.0_wp) * W_ii * alpha

            call get_damat_dc_dir(vec, dtrans, capi, capi, rvdw, self%kbc, gam, dG, dS)
            sigma_local(:, :) = sigma_local(:, :) - W_ii * dS * wsw * alpha

            call get_dcpair_dir(self%kbc, vec, dtrans, rvdw, capi, capi, dG, dS)
            sigma_local(:, :) = sigma_local(:, :) - p(iat) * cache%xtmp(iat) * dS * wsw * beta

            ! Apply self-image correction to potential v(iat) before weight calculation
            call get_cpair_dir(self%kbc, vec, dtrans, rvdw, capi, capi, ctmp)
            v(iat) = v(iat) - ctmp * wsw * p(iat)
         end do

         ! Diagonal contributions to weights
         qlocacc(iat) = qlocacc(iat) + self%kqeta(izp) * W_ii * cache%cmat(iat, iat) * alpha &
         & + v(iat) * self%kqchi(izp) * beta
         cnacc(iat) = cnacc(iat) - sqrt2pi * dradi / (radi**2) * W_ii * cache%cmat(iat, iat) * alpha &
         & + v(iat) * self%kcnchi(izp) * beta

         ! Intrinsic capacitance
         dtmp = (self%eta(izp) + self%kqeta(izp) * cache%qloc(iat) + sqrt2pi / radi) * W_ii
         gradient_local(:, iat) = gradient_local(:, iat) + dtmp * cache%dcdr(:, iat, iat) * alpha &
         & + p(iat) * cache%xtmp(iat) * cache%dcdr(:, iat, iat) * beta

         sigma_local(:, :)   = sigma_local(:, :)   + dtmp * cache%dcdL(:, :, iat) * alpha &
         & + p(iat) * cache%xtmp(iat) * cache%dcdL(:, :, iat) * beta
      end do
      !$omp end parallel do

      gradient = gradient + gradient_local
      sigma = sigma + sigma_local

      gradient_local = 0.0_wp
      sigma_local = 0.0_wp

      call self%ncoord%add_coordination_number_derivs(mol, dtrans, cnacc, gradient_local, sigma_local)
      call self%ncoord_en%add_coordination_number_derivs(mol, dtrans, qlocacc, gradient_local, sigma_local)

      gradient = gradient + gradient_local
      sigma = sigma + sigma_local

      deallocate(gradient_local, sigma_local, qlocacc, cnacc, dtrans, v)
   end subroutine get_grad_3d

end module multicharge_model_eeqbc
