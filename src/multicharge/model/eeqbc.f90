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
!> Provides implementation of the bond capacitor electronegativity
!> equilibration model (EEQ_BC)

!> Bond capacitor electronegativity equilibration charge model published in
!>
!> Thomas Froitzheim, Marcel Müller, Andreas Hansen, and Stefan Grimme,
!> *J. Chem. Phys.*, **2025**, 162, 214109.
!> DOI: [10.1063/5.0268978](https://dx.doi.org/10.1063/5.0268978)
!> Updated from of the parametrization and minor model changes published in
!>
!> Thomas Froitzheim, Marcel Müller, Andreas Hansen, and Stefan Grimme,
!> *ChemRxiv*, **2025**.
!> DOI: [10.26434/chemrxiv-2025-bjxvt](https://doi.org/10.26434/chemrxiv-2025-bjxvt)

module multicharge_model_eeqbc
   use mctc_env, only : error_type, wp, i8
   use mctc_io, only : structure_type
   use mctc_io_constants, only : pi
   use mctc_ncoord, only : cn_count, new_ncoord, ncoord_type
   use mctc_csrlist, only : csr_list, spgemv_csr, spsymv_csr
   use mctc_wignerseitz, only : wignerseitz_cell
   use multicharge_wignerseitz, only : new_wignerseitz_cell
   use multicharge_model_type, only : get_dir_trans, mchrg_model_type, mchrg_hessian
   use multicharge_blas, only : gemm, gemv, symv
   use multicharge_model_cache, only : mchrg_cache
   implicit none
   private

   public :: eeqbc_model, new_eeqbc_model

   !> EEQBC model type extending the shared charge-model state
   type, extends(mchrg_model_type) :: eeqbc_model
      !> Bond capacitance parameters
      real(wp), allocatable :: cap(:)

      !> Local charge scaling factor prefactor for chemical hardness
      real(wp) :: kqeta_pre

      !> Average coordination number
      real(wp), allocatable :: avg_cn(:)

      !> Exponent of error function in bond capacitance
      real(wp) :: kbc

      !> Van der Waals radii matrix (nat × nat)
      real(wp), allocatable :: rvdw(:, :)
   contains
      !> Update and allocate cache
      procedure :: update
      !> Calculate capacitance matrix
      procedure :: get_capacitance_matrix
      !> Calculate Coulomb matrix
      procedure :: get_coulomb_matrix
      !> Calculate alpha * dA/dR*q + beta * dX/dR
      procedure :: get_partial_derivs
      !> Calculate right-hand side (electronegativity vector)
      procedure :: get_xvec
      !> Calculate constraint matrix (molecular)
      procedure :: get_cmat_0d
      !> Calculate constraint matrix (molecular) using neighborlist
      procedure :: get_cmat_0d_list
      !> Calculate full constraint matrix (periodic)
      procedure :: get_cmat_3d
      !> Calculate constraint matrix derivatives (molecular)
      procedure :: get_dcmat_0d
      !> Calculate constraint matrix derivatives (molecular) using neighborlist
      procedure :: get_dcmat_0d_list
      !> Calculate constraint matrix derivatives (periodic)
      procedure :: get_dcmat_3d
      !> Calculate diagonal blocks of the constraint matrix second derivatives
      !> (molecular)
      procedure :: get_d2cmat_0d
      !> Calculate diagonal blocks of the constraint matrix second derivatives
      !> (molecular) using neighborlist
      procedure :: get_d2cmat_0d_list
      !> Calculate the explicit Hessian alpha * p^T d2A/dR2 q + beta * p^T d2X/dR2
      !> at fixed p and q (molecular)
      procedure :: get_explicit_hessian_dense
      !> Calculate the explicit Hessian at fixed p and q on a complete
      !> neighborlist in compressed storage (molecular)
      procedure :: get_explicit_hessian_list
      !> Calculate gradient
      procedure :: get_grad
   end type eeqbc_model

   real(wp), parameter :: sqrtpi = sqrt(pi)
   real(wp), parameter :: sqrt2pi = sqrt(2.0_wp / pi)

   !> Threshold of the CN and local charge derivatives on the neighborlist below
   !> which their products in the partial derivatives and the Hessian are skipped
   real(wp), parameter :: cnthr = epsilon(1.0_wp)**2
   real(wp), parameter :: eps = sqrt(epsilon(0.0_wp))

   !> Default exponent of error function in bond capacitance
   real(wp), parameter :: default_kbc = 0.65_wp

   !> Default cutoff radius
   real(wp), parameter :: cutoff = 25.0_wp

   !> Default scaling factor for the external electric field
   real(wp), parameter :: default_efield_scale = 10.0_wp


contains


!> Construct an EEQBC model from element-wise parameters
subroutine new_eeqbc_model(self, mol, error, chi, rad, &
   & eta, kcnchi, kqchi, kqeta, kqeta_pre, kcnrad, cap, avg_cn, rvdw, &
   & kbc, cutoff, cn_exp, rcov, en, cn_max, efield_scale)
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
   !> Local charge scaling factor with tanh-clipping for chemical hardness
   real(wp), intent(in) :: kqeta(:)
   !> Local charge scaling factor prefactor for chemical hardness
   real(wp), intent(in) :: kqeta_pre
   !> CN scaling factor for charge width
   real(wp), intent(in) :: kcnrad(:)
   !> Bond capacitance
   real(wp), intent(in) :: cap(:)
   !> Average coordination number
   real(wp), intent(in) :: avg_cn(:)
   !> Van-der-Waals radii
   real(wp), intent(in) :: rvdw(:, :)
   !> Exponent of error function in bond capacitance
   real(wp), intent(in), optional :: kbc
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
   !> Scaling factor for external electric field
   real(wp), intent(in), optional :: efield_scale

   self%chi = chi
   self%rad = rad
   self%eta = eta
   self%kcnchi = kcnchi
   self%kqchi = kqchi
   self%kqeta = kqeta
   self%kqeta_pre = kqeta_pre
   self%kcnrad = kcnrad
   self%cap = cap
   self%avg_cn = avg_cn
   self%rvdw = rvdw

   if (present(kbc)) then
      self%kbc = kbc
   else
      self%kbc = default_kbc
   end if

   if (present(efield_scale)) then
      self%efield_scale = efield_scale
   else
      self%efield_scale = default_efield_scale
   end if

   ! Coordination number
   call new_ncoord(self%ncoord, mol, cn_count%erf, error, &
      & cutoff=cutoff, kcn=cn_exp, rcov=rcov, cut=cn_max)
   ! Electronegativity weighted coordination number for local charge
   call new_ncoord(self%ncoord_en, mol, cn_count%erf_en, error, &
      & cutoff=cutoff, kcn=cn_exp, rcov=rcov, en=en, cut=cn_max)

end subroutine new_eeqbc_model

!> Update coordination numbers and local charges, their derivatives if the
!> charge derivatives or the Hessian are requested, and set up the Wigner-Seitz
!> cell if periodic
subroutine update(self, mol, cache, trans, grad, list, dqdr, hess)
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
   !> Multicharge neighborlist type
   type(csr_list), intent(in), optional :: list
   !> Flag for the charge derivatives w.r.t. positions and strain
   logical, intent(in), optional :: dqdr
   !> Flag for the Hessian
   logical, intent(in), optional :: hess

   logical :: deriv

   ! The charge derivatives and the Hessian need the position and strain
   ! derivatives of the CN and the local charges
   deriv = .false.
   if (present(dqdr)) deriv = dqdr
   if (present(hess)) deriv = deriv .or. hess

   cache%trans = trans

   ! Refer CN and local charge arrays in cache
   if (.not. allocated(cache%cn)) then
      allocate(cache%cn(mol%nat))
   end if
   if (.not. allocated(cache%qloc)) then
      allocate(cache%qloc(mol%nat))
   end if

   ! The charge derivatives and the Hessian also need the capacitance derivatives
   cache%grad = grad .or. deriv

   if (deriv) then
      if (.not. allocated(cache%dcndL)) allocate(cache%dcndL(3, 3, mol%nat))
      if (.not. allocated(cache%dqlocdL)) allocate(cache%dqlocdL(3, 3, mol%nat))
   end if

   if (present(list)) then
      if (deriv) then
         ! The pattern of the position derivatives follows the neighborlist
         if (allocated(cache%dcndrlist)) deallocate(cache%dcndrlist)
         if (allocated(cache%dqlocdrlist)) deallocate(cache%dqlocdrlist)
         allocate(cache%dcndrlist(3, size(list%nlat, kind=i8)))
         allocate(cache%dqlocdrlist(3, size(list%nlat, kind=i8)))
         call self%ncoord%get_coordination_number(mol, trans, cache%cn, &
            & dcndL=cache%dcndL, list=list, dcndrlist=cache%dcndrlist)
         call self%local_charge(mol, trans, cache%qloc, dqlocdL=cache%dqlocdL, &
            & list=list, dqlocdrlist=cache%dqlocdrlist)
      else
         call self%ncoord%get_coordination_number(mol, trans, cache%cn, list=list)
         call self%local_charge(mol, trans, cache%qloc, list=list)
      end if
   else
      if (deriv) then
         if (.not. allocated(cache%dcndr)) allocate(cache%dcndr(3, mol%nat, mol%nat))
         if (.not. allocated(cache%dqlocdr)) allocate(cache%dqlocdr(3, mol%nat, mol%nat))
         call self%ncoord%get_coordination_number(mol, trans, cache%cn, &
            & dcndr=cache%dcndr, dcndL=cache%dcndL)
         call self%local_charge(mol, trans, cache%qloc, dqlocdr=cache%dqlocdr, &
            & dqlocdL=cache%dqlocdL)
      else
         call self%ncoord%get_coordination_number(mol, trans, cache%cn)
         call self%local_charge(mol, trans, qloc=cache%qloc)
      end if
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
   !> Multicharge neighborlist type
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
         allocate(cache%clist(size(list%nlat, kind=i8)))
      end if
      ! neighborlist routines
      if (any(mol%periodic)) then
         call get_cmat_3d_list(self, mol, list, cache)
         ! Capacitance-matrix gradients
         if (cache%grad) then
            call get_dcmat_3d_list(self, mol, list, cache)
         end if
      else
         call get_cmat_0d_list(self, mol, list, cache)
         ! Capacitance-matrix gradients
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
         ! Capacitance-matrix gradients
         if (cache%grad) then
            call get_dcmat_0d(self, mol, cache%dcdr, cache%dcdL)
         end if
      end if
   end if

end subroutine get_capacitance_matrix

!> Build the electronegativity vector, including local-charge corrections
subroutine get_xvec(self, mol, ndim, cache, list, efield)
   !> EEQBC model type
   class(eeqbc_model), intent(in) :: self
   !> Structure type
   type(structure_type), intent(in) :: mol
   !> System size (number of atoms or atoms+1 if Lagrange multiplier used)
   integer, intent(in) :: ndim
   !> Multicharge cache
   type(mchrg_cache), intent(inout) :: cache
   !> Multicharge neighborlist type
   type(csr_list), intent(in), optional :: list
   !> External electric field
   real(wp), intent(in), optional :: efield(:)

   integer :: iat, izp
   integer(i8) :: img
   real(wp) :: ctmp, vec(3), rvdw, capi, wsw
   real(wp), allocatable :: dtrans(:, :)

   if (.not. allocated(cache%xtmp)) then
      allocate(cache%xtmp(ndim))
   end if

   if (.not. allocated(cache%xvec)) then
      allocate(cache%xvec(ndim))
   end if

   cache%xvec(:) = 0.0_wp

   !$omp parallel do default(none) schedule(runtime) &
   !$omp shared(mol, self, cache) &
   !$omp private(iat, izp)
   do iat = 1, mol%nat
      izp = mol%id(iat)
      cache%xtmp(iat) = -self%chi(izp) + self%kcnchi(izp) * cache%cn(iat) &
         & + self%kqchi(izp) * cache%qloc(iat)
   end do
   !$omp end parallel do

   ! Add external electric field to the RHS if present
   if (present(efield)) then
      !$omp parallel do default(none) schedule(runtime) &
      !$omp shared(mol, self, cache, efield) private(iat)
      do iat = 1, mol%nat
         cache%xtmp(iat) = cache%xtmp(iat) + self%efield_scale &
            & * dot_product(mol%xyz(:, iat), efield)
      end do
      !$omp end parallel do
   end if

   if (size(cache%xtmp) == mol%nat + 1) then
      cache%xtmp(mol%nat + 1) = mol%charge
   end if

   if (present(list)) then
      if (list%complete) then
         call spgemv_csr(mol%nat, cache%clist, list%inl, list%nlat, &
            & cache%xtmp, cache%xvec)
      else
         call spsymv_csr(mol%nat, cache%clist, list%inl, list%nlat, &
            & cache%xtmp, cache%xvec)
      end if
   else
      call gemv(cache%cmat, cache%xtmp, cache%xvec)
   end if

   ! Periodic contributions
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

            if (cache%wsc%nimg_list(list%inl(iat)) > 0) then
               wsw = 1.0_wp / real(cache%wsc%nimg_list(list%inl(iat)), wp)
               do img = cache%wsc%itr_list(list%inl(iat)), &
               & cache%wsc%itr_list(list%inl(iat) + 1) - 1
                  vec = cache%wsc%trans(:, cache%wsc%tridx_list(img))
                  call get_cpair_dir(self%kbc, vec, dtrans, rvdw, capi, capi, ctmp)
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

!> Compute the linear combination of the Coulomb matrix derivatives (multiplied
!> by the charge vector) and the electronegativity vector derivatives,
!> alpha * dA/dR*q + beta * dX/dR, w.r.t. positions and strain
subroutine get_partial_derivs(self, mol, ndim, cache, alpha, beta, list)
   !> EEQBC model type
   class(eeqbc_model), intent(in) :: self
   !> Structure type
   type(structure_type), intent(in) :: mol
   !> System size
   integer, intent(in) :: ndim
   !> Multicharge cache
   type(mchrg_cache), intent(inout) :: cache
   !> dA/dR*q multiplier
   real(wp), intent(in) :: alpha
   !> dX/dR multiplier
   real(wp), intent(in) :: beta
   !> Multicharge neighborlist type (complete)
   type(csr_list), intent(in), optional :: list

   if (.not. allocated(cache%dabdL)) allocate(cache%dabdL(3, 3, ndim))
   cache%dabdL = 0.0_wp

   if (present(list)) then
      if (.not. allocated(cache%dabdrlist)) then
         allocate(cache%dabdrlist(3, size(list%nlat, kind=i8)))
      end if
      cache%dabdrlist = 0.0_wp
      if (any(mol%periodic)) then
         call get_partial_derivs_3d_list(self, mol, list, cache, alpha, beta)
      else
         call get_partial_derivs_0d_list(self, mol, list, cache, alpha, beta)
      end if
   else
      if (.not. allocated(cache%dabdr)) allocate(cache%dabdr(3, mol%nat, ndim))
      cache%dabdr = 0.0_wp
      if (any(mol%periodic)) then
         call get_partial_derivs_3d(self, mol, ndim, cache, alpha, beta)
      else
         call get_partial_derivs_0d(self, mol, ndim, cache, alpha, beta)
      end if
   end if

end subroutine get_partial_derivs

!> Compute the alpha * dA/dR*q + beta * dX/dR, for a non-periodic system.
subroutine get_partial_derivs_0d(self, mol, ndim, cache, alpha, beta)
   !> EEQBC model type
   class(eeqbc_model), intent(in) :: self
   !> Structure type
   type(structure_type), intent(in) :: mol
   !> System size
   integer, intent(in) :: ndim
   !> Multicharge cache
   type(mchrg_cache), intent(inout) :: cache
   !> dA/dR*q multiplier
   real(wp), intent(in) :: alpha
   !> dX/dR multiplier
   real(wp), intent(in) :: beta

   integer :: iat, izp, jat, jzp
   real(wp) :: vec(3), r2, gam, arg, dtmp, norm_cn, hardi, wtmp, wcni
   real(wp) :: radi, radj, dradi, dradj, dG(3), dS(3, 3)
   real(wp), allocatable :: dtmpdr(:, :, :), dtmpdL(:, :, :)

   allocate(dtmpdr(3, mol%nat, ndim), dtmpdL(3, 3, ndim), source=0.0_wp)

   ! CN and effective charge derivative of the electronegativity
   !$omp parallel do default(none) schedule(runtime) &
   !$omp shared(self, mol, cache, dtmpdr, dtmpdL) private(iat, izp)
   do iat = 1, mol%nat
      izp = mol%id(iat)
      dtmpdr(:, :, iat) = self%kcnchi(izp) * cache%dcndr(:, :, iat) &
         & + self%kqchi(izp) * cache%dqlocdr(:, :, iat)
      dtmpdL(:, :, iat) = self%kcnchi(izp) * cache%dcndL(:, :, iat) &
         & + self%kqchi(izp) * cache%dqlocdL(:, :, iat)
   end do
   !$omp end parallel do

   ! EN derivative through the capacitance matrix, initializes dabdr and dabdL
   call gemm(dtmpdr, cache%cmat, cache%dabdr, alpha=beta)
   call gemm(dtmpdL, cache%cmat, cache%dabdL, alpha=beta)

   ! Each iteration only updates the derivatives of column iat (last index)
   !$omp parallel do default(none) schedule(runtime) &
   !$omp shared(cache, mol, self, alpha, beta) &
   !$omp private(iat, izp, jat, jzp, gam, vec, r2, dtmp, norm_cn, arg, hardi) &
   !$omp private(wtmp, wcni, radi, radj, dradi, dradj, dG, dS)
   do iat = 1, mol%nat
      izp = mol%id(iat)
      ! Effective charge width of i
      norm_cn = 1.0_wp / self%avg_cn(izp)
      radi = self%rad(izp) * exp(-self%kcnrad(izp) * cache%cn(iat) * norm_cn)
      dradi = -self%kcnrad(izp) * norm_cn * radi
      ! Effective hardness of i
      hardi = self%eta(izp) + self%kqeta_pre &
      & * tanh(self%kqeta(izp) * cache%qloc(iat)) + sqrt2pi / radi
      wcni = 0.0_wp
      do jat = 1, mol%nat
         if (jat == iat) cycle
         jzp = mol%id(jat)
         vec = mol%xyz(:, jat) - mol%xyz(:, iat)
         r2 = vec(1)**2 + vec(2)**2 + vec(3)**2
         ! Effective charge width of j
         norm_cn = 1.0_wp / self%avg_cn(jzp)
         radj = self%rad(jzp) * exp(-self%kcnrad(jzp) * cache%cn(jat) * norm_cn)
         dradj = -self%kcnrad(jzp) * norm_cn * radj

         ! EN derivative: capacitance matrix derivative
         cache%dabdr(:, iat, iat) = +beta * cache%xtmp(jat) * cache%dcdr(:, iat, &
         & jat) + cache%dabdr(:, iat, iat)
         cache%dabdr(:, jat, iat) = +beta * (cache%xtmp(jat) - cache%xtmp(iat)) &
         & * cache%dcdr(:, jat, iat) + cache%dabdr(:, jat, iat)
         cache%dabdL(:, :, iat) = +beta * cache%xtmp(jat) * spread(cache%dcdr(:, &
         & iat, jat), 1, 3) * spread(-vec, 2, 3) + cache%dabdL(:, :, iat)

         ! Coulomb interaction of Gaussian charges
         gam = 1.0_wp / sqrt(radi**2 + radj**2)
         arg = gam * gam * r2

         ! Explicit derivative
         dtmp = 2.0_wp * gam * exp(-arg) / (sqrtpi * r2) &
            & - erf(sqrt(arg)) / (r2 * sqrt(r2))
         dG(:) = alpha * dtmp * vec * cache%vrhs(jat) * cache%cmat(jat, iat)
         dS(:, :) = spread(dG, 1, 3) * spread(vec, 2, 3)
         cache%dabdr(:, iat, iat) = cache%dabdr(:, iat, iat) - dG
         cache%dabdr(:, jat, iat) = cache%dabdr(:, jat, iat) + dG
         cache%dabdL(:, :, iat) = cache%dabdL(:, :, iat) + dS

         ! Effective charge width derivative, CN of i is collected in wcni
         wtmp = alpha * 2.0_wp * exp(-arg) / sqrtpi * cache%vrhs(jat) &
         & * cache%cmat(jat, iat)
         wcni = wcni - wtmp * radi * dradi * gam**3
         dtmp = -wtmp * radj * dradj * gam**3
         cache%dabdr(:, :, iat) = +dtmp * cache%dcndr(:, :, jat) &
            & + cache%dabdr(:, :, iat)
         cache%dabdL(:, :, iat) = +dtmp * cache%dcndL(:, :, jat) &
            & + cache%dabdL(:, :, iat)

         ! Capacitance derivative off-diagonal
         dtmp = alpha * erf(sqrt(r2) * gam) / (sqrt(r2)) * cache%vrhs(jat)
         cache%dabdr(:, iat, iat) = -dtmp * cache%dcdr(:, jat, iat) &
         & + cache%dabdr(:, iat, iat)
         cache%dabdr(:, jat, iat) = +dtmp * cache%dcdr(:, jat, iat) &
         & + cache%dabdr(:, jat, iat)
         cache%dabdL(:, :, iat) = -dtmp * spread(cache%dcdr(:, iat, jat), 2, 3) &
         & * spread(vec, 1, 3) + cache%dabdL(:, :, iat)

         ! Capacitance derivative diagonal
         dtmp = alpha * hardi * cache%vrhs(iat)
         cache%dabdr(:, jat, iat) = -dtmp * cache%dcdr(:, jat, iat) &
         & + cache%dabdr(:, jat, iat)
      end do

      ! EN derivative: capacitance matrix derivative diagonal
      cache%dabdr(:, iat, iat) = +beta * cache%xtmp(iat) * cache%dcdr(:, iat, &
      & iat) + cache%dabdr(:, iat, iat)
      cache%dabdL(:, :, iat) = +beta * cache%xtmp(iat) * cache%dcdL(:, :, &
      & iat) + cache%dabdL(:, :, iat)

      ! Hardness derivative
      dtmp = alpha * self%kqeta_pre * self%kqeta(izp) &
      & / cosh(self%kqeta(izp) * cache%qloc(iat))**2 * cache%vrhs(iat) &
      & * cache%cmat(iat, iat)
      cache%dabdr(:, :, iat) = +dtmp * cache%dqlocdr(:, :, iat) &
         & + cache%dabdr(:, :, iat)
      cache%dabdL(:, :, iat) = +dtmp * cache%dqlocdL(:, :, iat) &
         & + cache%dabdL(:, :, iat)

      ! Effective charge width derivative
      dtmp = -alpha * sqrt2pi * dradi / (radi**2) * cache%vrhs(iat) * cache%cmat(iat, iat)
      cache%dabdr(:, :, iat) = +(dtmp + wcni) * cache%dcndr(:, :, iat) &
         & + cache%dabdr(:, :, iat)
      cache%dabdL(:, :, iat) = +(dtmp + wcni) * cache%dcndL(:, :, iat) &
         & + cache%dabdL(:, :, iat)

      ! Capacitance derivative
      dtmp = alpha * hardi * cache%vrhs(iat)
      cache%dabdr(:, iat, iat) = +dtmp * cache%dcdr(:, iat, iat) &
      & + cache%dabdr(:, iat, iat)
      cache%dabdL(:, :, iat) = +dtmp * cache%dcdL(:, :, iat) + cache%dabdL(:, :, iat)
   end do
   !$omp end parallel do

end subroutine get_partial_derivs_0d

!> Compute the alpha * dA/dR*q + beta * dX/dR, for a periodic system.
subroutine get_partial_derivs_3d(self, mol, ndim, cache, alpha, beta)
   !> EEQBC model type
   class(eeqbc_model), intent(in) :: self
   !> Structure type
   type(structure_type), intent(in) :: mol
   !> System size
   integer, intent(in) :: ndim
   !> Multicharge cache
   type(mchrg_cache), intent(inout) :: cache
   !> dA/dR*q multiplier
   real(wp), intent(in) :: alpha
   !> dX/dR multiplier
   real(wp), intent(in) :: beta

   integer :: iat, jat, izp, jzp, img
   real(wp) :: vec(3), gam, dtmp, norm_cn, rvdw, wsw, dgam, wgam, ctmp, hardi, hardj
   real(wp) :: radi, radj, dradi, dradj, dG(3), dS(3, 3), dgamdL(3, 3)
   real(wp) :: capi, capj
   real(wp), allocatable :: dtrans(:, :), dgamdr(:, :)
   real(wp), allocatable :: dtmpdr(:, :, :), dtmpdL(:, :, :)

   ! Thread-private arrays for reduction
   real(wp), allocatable :: dabdr_local(:, :, :), dabdL_local(:, :, :)

   call get_dir_trans(mol, dtrans, cutoff)

   allocate(dtmpdr(3, mol%nat, ndim), dtmpdL(3, 3, ndim), source=0.0_wp)

   ! CN and effective charge derivative of the electronegativity
   !$omp parallel do default(none) schedule(runtime) &
   !$omp shared(self, mol, cache, dtmpdr, dtmpdL) private(iat, izp)
   do iat = 1, mol%nat
      izp = mol%id(iat)
      dtmpdr(:, :, iat) = self%kcnchi(izp) * cache%dcndr(:, :, iat) &
         & + self%kqchi(izp) * cache%dqlocdr(:, :, iat)
      dtmpdL(:, :, iat) = self%kcnchi(izp) * cache%dcndL(:, :, iat) &
         & + self%kqchi(izp) * cache%dqlocdL(:, :, iat)
   end do
   !$omp end parallel do

   ! EN derivative through the capacitance matrix, initializes dabdr and dabdL
   call gemm(dtmpdr, cache%cmat, cache%dabdr, alpha=beta)
   call gemm(dtmpdL, cache%cmat, cache%dabdL, alpha=beta)

   !$omp parallel default(none) &
   !$omp shared(self, mol, cache, dtrans, alpha, beta, dtmpdr, dtmpdL) &
   !$omp private(iat, izp, jat, jzp, img, gam, vec, dtmp, norm_cn, rvdw, wsw) &
   !$omp private(radi, radj, dradi, dradj, capi, capj, dgam, wgam, dgamdr, dgamdL) &
   !$omp private(ctmp, hardi, hardj, dG, dS) &
   !$omp private(dabdr_local, dabdL_local)
   allocate(dgamdr(3, mol%nat))
   allocate(dabdr_local(3, mol%nat, size(cache%dabdr, 3)), source=0.0_wp)
   allocate(dabdL_local(3, 3, size(cache%dabdL, 3)), source=0.0_wp)
   !$omp do schedule(runtime)
   do iat = 1, mol%nat
      izp = mol%id(iat)
      ! Effective charge width of i
      norm_cn = 1.0_wp / self%avg_cn(izp)
      radi = self%rad(izp) * exp(-self%kcnrad(izp) * cache%cn(iat) * norm_cn)
      dradi = -self%kcnrad(izp) * norm_cn * radi
      capi = self%cap(izp)
      ! Effective hardness of i
      hardi = self%eta(izp) + self%kqeta_pre &
      & * tanh(self%kqeta(izp) * cache%qloc(iat)) + sqrt2pi / radi

      do jat = 1, mol%nat
         jzp = mol%id(jat)
         capj = self%cap(jzp)
         rvdw = self%rvdw(izp, jzp)

         ! EN derivative: capacitance matrix derivative
         dabdr_local(:, iat, iat) = dabdr_local(:, iat, iat) &
            & + beta * cache%xtmp(jat) * cache%dcdr(:, iat, jat)
         dabdr_local(:, iat, jat) = dabdr_local(:, iat, jat) &
            & + beta * (cache%xtmp(iat) - cache%xtmp(jat)) * cache%dcdr(:, iat, jat)

         wsw = 1.0_wp / real(cache%wsc%nimg(iat, jat), wp)
         do img = 1, cache%wsc%nimg(iat, jat)
            vec = mol%xyz(:, jat) - mol%xyz(:, iat) + cache%wsc%trans(:, &
            & cache%wsc%tridx(img, jat, iat))
            call get_dcpair_dir(self%kbc, vec, dtrans, rvdw, capi, capj, dG, dS)
            dabdL_local(:, :, iat) = dabdL_local(:, :, iat) &
               & - beta * wsw * dS * cache%xtmp(jat)
         end do

         ! Coulomb matrix derivatives for each unordered pair
         if (jat >= iat) cycle

         ! Effective charge width of j
         norm_cn = 1.0_wp / self%avg_cn(jzp)
         radj = self%rad(jzp) * exp(-self%kcnrad(jzp) * cache%cn(jat) * norm_cn)
         dradj = -self%kcnrad(jzp) * norm_cn * radj
         ! Effective hardness of j
         hardj = self%eta(jzp) + self%kqeta_pre &
         & * tanh(self%kqeta(jzp) * cache%qloc(jat)) + sqrt2pi / radj

         ! Coulomb interaction of Gaussian charges
         gam = 1.0_wp / sqrt(radi**2 + radj**2)
         dgamdr(:, :) = -(radi * dradi * cache%dcndr(:, :, iat) &
         & + radj * dradj * cache%dcndr(:, :, jat)) * gam**3
         dgamdL(:, :) = -(radi * dradi * cache%dcndL(:, :, iat) &
         & + radj * dradj * cache%dcndL(:, :, jat)) * gam**3

         wgam = 0.0_wp
         wsw = 1.0_wp / real(cache%wsc%nimg(jat, iat), wp)
         do img = 1, cache%wsc%nimg(jat, iat)
            vec = mol%xyz(:, jat) - mol%xyz(:, iat) + cache%wsc%trans(:, &
            & cache%wsc%tridx(img, jat, iat))

            call get_damat_dir(vec, dtrans, capi, capj, rvdw, self%kbc, gam, dG, dS, dgam)
            dG = alpha * dG * wsw
            dS = alpha * dS * wsw
            wgam = wgam + alpha * dgam * wsw

            ! Explicit derivative
            dabdr_local(:, iat, iat) = dabdr_local(:, iat, iat) - dG * cache%vrhs(jat)
            dabdr_local(:, jat, jat) = dabdr_local(:, jat, jat) + dG * cache%vrhs(iat)
            dabdr_local(:, iat, jat) = dabdr_local(:, iat, jat) - dG * cache%vrhs(iat)
            dabdr_local(:, jat, iat) = dabdr_local(:, jat, iat) + dG * cache%vrhs(jat)
            dabdL_local(:, :, jat) = dabdL_local(:, :, jat) + dS * cache%vrhs(iat)
            dabdL_local(:, :, iat) = dabdL_local(:, :, iat) + dS * cache%vrhs(jat)

            call get_damat_dc_dir(vec, dtrans, capi, capj, rvdw, self%kbc, gam, dG, dS)
            dG = alpha * dG * wsw
            dS = alpha * dS * wsw

            ! Capacitance derivative off-diagonal
            dabdr_local(:, iat, iat) = dabdr_local(:, iat, iat) + cache%vrhs(jat) * dG
            dabdr_local(:, jat, jat) = dabdr_local(:, jat, jat) - cache%vrhs(iat) * dG
            dabdr_local(:, jat, iat) = dabdr_local(:, jat, iat) - cache%vrhs(jat) * dG
            dabdr_local(:, iat, jat) = dabdr_local(:, iat, jat) + cache%vrhs(iat) * dG
            dabdL_local(:, :, jat) = dabdL_local(:, :, jat) - cache%vrhs(iat) * dS
            dabdL_local(:, :, iat) = dabdL_local(:, :, iat) - cache%vrhs(jat) * dS

            call get_dcpair_dir(self%kbc, vec, dtrans, rvdw, capi, capj, dG, dS)
            dG = alpha * dG * wsw

            ! Capacitance derivative diagonal
            dabdr_local(:, jat, iat) = dabdr_local(:, jat, iat) &
               & + hardi * cache%vrhs(iat) * dG
            dabdr_local(:, iat, jat) = dabdr_local(:, iat, jat) &
               & - hardj * cache%vrhs(jat) * dG
         end do

         ! Effective charge width derivative
         dabdr_local(:, :, iat) = dabdr_local(:, :, iat) - wgam * cache%vrhs(jat) * dgamdr
         dabdr_local(:, :, jat) = dabdr_local(:, :, jat) - wgam * cache%vrhs(iat) * dgamdr
         dabdL_local(:, :, iat) = dabdL_local(:, :, iat) - wgam * cache%vrhs(jat) * dgamdL
         dabdL_local(:, :, jat) = dabdL_local(:, :, jat) - wgam * cache%vrhs(iat) * dgamdL
      end do

      ! EN derivative: capacitance matrix derivative diagonal
      dabdL_local(:, :, iat) = dabdL_local(:, :, iat) &
         & + beta * cache%xtmp(iat) * cache%dcdL(:, :, iat)

      ! Self-images, i = j and T != 0
      gam = 1.0_wp / sqrt(2.0_wp * radi**2)
      dgamdr(:, :) = -2.0_wp * radi * dradi * cache%dcndr(:, :, iat) * gam**3
      dgamdL(:, :) = -2.0_wp * radi * dradi * cache%dcndL(:, :, iat) * gam**3
      rvdw = self%rvdw(izp, izp)
      wsw = 1.0_wp / real(cache%wsc%nimg(iat, iat), wp)
      do img = 1, cache%wsc%nimg(iat, iat)
         vec = cache%wsc%trans(:, cache%wsc%tridx(img, iat, iat))

         ! EN derivative through the self-image capacitance
         call get_cpair_dir(self%kbc, vec, dtrans, rvdw, capi, capi, ctmp)
         ctmp = beta * ctmp * wsw
         dabdr_local(:, :, iat) = dabdr_local(:, :, iat) - ctmp * dtmpdr(:, :, iat)
         dabdL_local(:, :, iat) = dabdL_local(:, :, iat) - ctmp * dtmpdL(:, :, iat)

         call get_damat_dir(vec, dtrans, capi, capi, rvdw, self%kbc, gam, dG, dS, dgam)
         dS = alpha * dS * wsw
         dgam = alpha * dgam * wsw

         ! Explicit derivative
         dabdL_local(:, :, iat) = dabdL_local(:, :, iat) + dS * cache%vrhs(iat)

         ! Effective charge width derivative
         dabdr_local(:, :, iat) = dabdr_local(:, :, iat) - cache%vrhs(iat) * dgam * dgamdr
         dabdL_local(:, :, iat) = dabdL_local(:, :, iat) - cache%vrhs(iat) * dgam * dgamdL

         ! Capacitance derivative
         call get_damat_dc_dir(vec, dtrans, capi, capi, rvdw, self%kbc, gam, dG, dS)
         dabdL_local(:, :, iat) = dabdL_local(:, :, iat) &
            & - alpha * cache%vrhs(iat) * dS * wsw
      end do

      ! Hardness derivative
      dtmp = alpha * self%kqeta_pre * self%kqeta(izp) &
      & / cosh(self%kqeta(izp) * cache%qloc(iat))**2 * cache%vrhs(iat) &
      & * cache%cmat(iat, iat)
      dabdr_local(:, :, iat) = dabdr_local(:, :, iat) + dtmp * cache%dqlocdr(:, :, iat)
      dabdL_local(:, :, iat) = dabdL_local(:, :, iat) + dtmp * cache%dqlocdL(:, :, iat)

      ! Effective charge width derivative
      dtmp = -alpha * sqrt2pi * dradi / (radi**2) * cache%vrhs(iat) * cache%cmat(iat, iat)
      dabdr_local(:, :, iat) = dabdr_local(:, :, iat) + dtmp * cache%dcndr(:, :, iat)
      dabdL_local(:, :, iat) = dabdL_local(:, :, iat) + dtmp * cache%dcndL(:, :, iat)

      ! Capacitance derivative
      dtmp = alpha * hardi * cache%vrhs(iat)
      dabdr_local(:, iat, iat) = dabdr_local(:, iat, iat) + dtmp * cache%dcdr(:, iat, iat)
      dabdL_local(:, :, iat) = dabdL_local(:, :, iat) + dtmp * cache%dcdL(:, :, iat)
   end do
   !$omp end do
   !$omp critical (get_partial_derivs_3d_)
   cache%dabdr(:, :, :) = cache%dabdr + dabdr_local
   cache%dabdL(:, :, :) = cache%dabdL + dabdL_local
   !$omp end critical (get_partial_derivs_3d_)
   deallocate(dgamdr, dabdr_local, dabdL_local)
   !$omp end parallel

end subroutine get_partial_derivs_3d

!> Compute alpha * dA/dR*q + beta * dX/dR for a non-periodic system using a
!> complete neighborlist, the products with CN or local charge derivatives
!> below cnthr are skipped
subroutine get_partial_derivs_0d_list(self, mol, list, cache, alpha, beta)
   !> EEQBC model type
   class(eeqbc_model), intent(in) :: self
   !> Molecular structure data
   type(structure_type), intent(in) :: mol
   !> Multicharge neighborlist type (complete)
   type(csr_list), intent(in) :: list
   !> Multicharge cache
   type(mchrg_cache), intent(inout) :: cache
   !> dA/dR*q multiplier
   real(wp), intent(in) :: alpha
   !> dX/dR multiplier
   real(wp), intent(in) :: beta

   integer :: iat, izp, jat, jzp
   integer(i8) :: kat, lat, ist, ien, cat, xpos
   real(wp) :: vec(3), r1, r2, gam, arg, dtmp, norm_cn, hardi, wtmp, wcni, cij
   real(wp) :: radi, radj, dradi, dradj, dG(3), dS(3, 3), dGc(3), dSc(3, 3)
   real(wp), allocatable :: dtmpdr(:, :), dtmpdL(:, :, :)
   integer(i8), allocatable :: cninl(:), cnidx(:)

   ! Thread-private position of each atom in the row of atom iat
   integer(i8), allocatable :: pos(:)

   allocate(dtmpdr(3, size(list%nlat, kind=i8)), dtmpdL(3, 3, mol%nat))

   ! Entries with CN or local charge derivatives, the products with the
   ! derivatives of the second neighbors only run over these entries
   call get_cn_entries(mol, list, cache, cninl, cnidx)

   ! CN and effective charge derivative of the electronegativity in the row of
   ! atom iat, without a maximum CN cutoff dCN_i/dR_j = -dCN_j/dR_i and
   ! dqloc_i/dR_j = dqloc_j/dR_i give them from the entries of the row
   !$omp parallel do default(none) schedule(runtime) &
   !$omp shared(self, mol, list, cache, dtmpdr, dtmpdL) &
   !$omp private(iat, izp, ist, ien)
   do iat = 1, mol%nat
      izp = mol%id(iat)
      ist = list%inl(iat)
      ien = list%inl(iat + 1) - 1
      dtmpdr(:, ist) = self%kcnchi(izp) * cache%dcndrlist(:, ist) &
         & + self%kqchi(izp) * cache%dqlocdrlist(:, ist)
      dtmpdr(:, ist+1:ien) = -self%kcnchi(izp) * cache%dcndrlist(:, ist+1:ien) &
         & + self%kqchi(izp) * cache%dqlocdrlist(:, ist+1:ien)
      dtmpdL(:, :, iat) = self%kcnchi(izp) * cache%dcndL(:, :, iat) &
         & + self%kqchi(izp) * cache%dqlocdL(:, :, iat)
   end do
   !$omp end parallel do

   cache%dabdL(:, :, :) = 0.0_wp

   ! Each iteration only updates the derivatives of component iat (row iat),
   ! derivatives w.r.t. atoms outside of the row are dropped
   !$omp parallel default(none) &
   !$omp shared(cache, mol, list, self, alpha, beta, dtmpdr, dtmpdL, cninl, cnidx) &
   !$omp private(iat, izp, jat, jzp, kat, lat, ist, cat, xpos, gam, vec, r1, r2) &
   !$omp private(dtmp, norm_cn, arg, hardi, wtmp, wcni, cij, radi, radj, dradi) &
   !$omp private(dradj, dG, dS, dGc, dSc, pos)
   allocate(pos(mol%nat), source=0_i8)
   !$omp do schedule(runtime)
   do iat = 1, mol%nat
      izp = mol%id(iat)
      ist = list%inl(iat)
      do kat = list%inl(iat), list%inl(iat + 1) - 1
         pos(list%nlat(kat)) = kat
      end do

      ! EN derivative through the capacitance matrix
      do kat = list%inl(iat), list%inl(iat + 1) - 1
         jat = list%nlat(kat)
         do cat = cninl(jat), cninl(jat + 1) - 1
            lat = cnidx(cat)
            xpos = pos(list%nlat(lat))
            if (xpos == 0) cycle
            cache%dabdrlist(:, xpos) = cache%dabdrlist(:, xpos) &
               & + beta * cache%clist(kat) * dtmpdr(:, lat)
         end do
         cache%dabdL(:, :, iat) = cache%dabdL(:, :, iat) &
            & + beta * cache%clist(kat) * dtmpdL(:, :, jat)
      end do

      ! Effective charge width of i
      norm_cn = 1.0_wp / self%avg_cn(izp)
      radi = self%rad(izp) * exp(-self%kcnrad(izp) * cache%cn(iat) * norm_cn)
      dradi = -self%kcnrad(izp) * norm_cn * radi
      ! Effective hardness of i
      hardi = self%eta(izp) + self%kqeta_pre &
      & * tanh(self%kqeta(izp) * cache%qloc(iat)) + sqrt2pi / radi
      wcni = 0.0_wp
      do kat = list%inl(iat) + 1, list%inl(iat + 1) - 1
         jat = list%nlat(kat)
         jzp = mol%id(jat)
         vec = mol%xyz(:, jat) - mol%xyz(:, iat)
         r2 = vec(1)**2 + vec(2)**2 + vec(3)**2
         r1 = sqrt(r2)
         cij = cache%clist(kat)
         ! Effective charge width of j
         norm_cn = 1.0_wp / self%avg_cn(jzp)
         radj = self%rad(jzp) * exp(-self%kcnrad(jzp) * cache%cn(jat) * norm_cn)
         dradj = -self%kcnrad(jzp) * norm_cn * radj

         ! Capacitance derivative, dC(i, j)/dR(i) = dGc = -dC(i, j)/dR(j)
         call get_dcpair(self%kbc, vec, self%rvdw(izp, jzp), self%cap(izp), &
            & self%cap(jzp), dGc, dSc)

         ! EN derivative: capacitance matrix derivative
         cache%dabdrlist(:, ist) = cache%dabdrlist(:, ist) + beta * cache%xtmp(jat) * dGc
         cache%dabdrlist(:, kat) = cache%dabdrlist(:, kat) &
            & - beta * (cache%xtmp(jat) - cache%xtmp(iat)) * dGc
         cache%dabdL(:, :, iat) = cache%dabdL(:, :, iat) - beta * cache%xtmp(jat) * dSc

         ! Coulomb interaction of Gaussian charges
         gam = 1.0_wp / sqrt(radi**2 + radj**2)
         arg = gam * gam * r2

         ! Explicit derivative
         dtmp = 2.0_wp * gam * exp(-arg) / (sqrtpi * r2) - erf(sqrt(arg)) / (r2 * r1)
         dG(:) = alpha * dtmp * vec * cache%vrhs(jat) * cij
         dS(:, :) = spread(dG, 1, 3) * spread(vec, 2, 3)
         cache%dabdrlist(:, ist) = cache%dabdrlist(:, ist) - dG
         cache%dabdrlist(:, kat) = cache%dabdrlist(:, kat) + dG
         cache%dabdL(:, :, iat) = cache%dabdL(:, :, iat) + dS

         ! Effective charge width derivative, CN of i is collected in wcni
         wtmp = alpha * 2.0_wp * exp(-arg) / sqrtpi * cache%vrhs(jat) * cij
         wcni = wcni - wtmp * radi * dradi * gam**3
         dtmp = -wtmp * radj * dradj * gam**3
         lat = list%inl(jat)
         cache%dabdrlist(:, kat) = cache%dabdrlist(:, kat) &
            & + dtmp * cache%dcndrlist(:, lat)
         do cat = cninl(jat) + 1, cninl(jat + 1) - 1
            lat = cnidx(cat)
            xpos = pos(list%nlat(lat))
            if (xpos == 0) cycle
            cache%dabdrlist(:, xpos) = cache%dabdrlist(:, xpos) - &
               & dtmp * cache%dcndrlist(:, lat)
         end do
         cache%dabdL(:, :, iat) = cache%dabdL(:, :, iat) &
            & + dtmp * cache%dcndL(:, :, jat)

         ! Capacitance derivative off-diagonal
         dtmp = alpha * erf(r1 * gam) / r1 * cache%vrhs(jat)
         cache%dabdrlist(:, ist) = cache%dabdrlist(:, ist) + dtmp * dGc
         cache%dabdrlist(:, kat) = cache%dabdrlist(:, kat) - dtmp * dGc
         cache%dabdL(:, :, iat) = cache%dabdL(:, :, iat) - dtmp * dSc

         ! Capacitance derivative diagonal
         dtmp = alpha * hardi * cache%vrhs(iat)
         cache%dabdrlist(:, kat) = cache%dabdrlist(:, kat) + dtmp * dGc
      end do

      ! EN derivative: capacitance matrix derivative diagonal
      cache%dabdrlist(:, ist) = cache%dabdrlist(:, ist) &
         & + beta * cache%xtmp(iat) * cache%dcdrdiag(:, iat)
      cache%dabdL(:, :, iat) = cache%dabdL(:, :, iat) &
         & + beta * cache%xtmp(iat) * cache%dcdL(:, :, iat)

      ! Hardness derivative
      dtmp = alpha * self%kqeta_pre * self%kqeta(izp) &
      & / cosh(self%kqeta(izp) * cache%qloc(iat))**2 * cache%vrhs(iat) &
      & * cache%clist(list%inl(iat))
      do kat = list%inl(iat), list%inl(iat + 1) - 1
         cache%dabdrlist(:, kat) = cache%dabdrlist(:, kat) + &
            &dtmp * cache%dqlocdrlist(:, kat)
      end do
      cache%dabdL(:, :, iat) = cache%dabdL(:, :, iat) &
         & + dtmp * cache%dqlocdL(:, :, iat)

      ! Effective charge width derivative
      dtmp = -alpha * sqrt2pi * dradi / (radi**2) * cache%vrhs(iat) &
      & * cache%clist(list%inl(iat)) + wcni
      cache%dabdrlist(:, ist) = cache%dabdrlist(:, ist) &
         & + dtmp * cache%dcndrlist(:, ist)
      do kat = list%inl(iat) + 1, list%inl(iat + 1) - 1
         cache%dabdrlist(:, kat) = cache%dabdrlist(:, kat) &
            & - dtmp * cache%dcndrlist(:, kat)
      end do
      cache%dabdL(:, :, iat) = cache%dabdL(:, :, iat) + dtmp * cache%dcndL(:, :, iat)

      ! Capacitance derivative
      dtmp = alpha * hardi * cache%vrhs(iat)
      cache%dabdrlist(:, ist) = cache%dabdrlist(:, ist) &
         & + dtmp * cache%dcdrdiag(:, iat)
      cache%dabdL(:, :, iat) = cache%dabdL(:, :, iat) + dtmp * cache%dcdL(:, :, iat)

      do kat = list%inl(iat), list%inl(iat + 1) - 1
         pos(list%nlat(kat)) = 0
      end do
   end do
   !$omp end do
   deallocate(pos)
   !$omp end parallel

end subroutine get_partial_derivs_0d_list

!> Compute alpha * dA/dR*q + beta * dX/dR for a periodic system using a
!> complete neighborlist with Wigner-Seitz images
subroutine get_partial_derivs_3d_list(self, mol, list, cache, alpha, beta)
   !> EEQBC model type
   class(eeqbc_model), intent(in) :: self
   !> Molecular structure data
   type(structure_type), intent(in) :: mol
   !> Multicharge neighborlist type (complete)
   type(csr_list), intent(in) :: list
   !> Multicharge cache
   type(mchrg_cache), intent(inout) :: cache
   !> dA/dR*q multiplier
   real(wp), intent(in) :: alpha
   !> dX/dR multiplier
   real(wp), intent(in) :: beta

   integer :: iat, izp, jat, jzp, xat
   integer(i8) :: kat, lat, img, ist, ien
   real(wp) :: vec(3), gam, dtmp, norm_cn, rvdw, wsw, dgam, wgam, ctmp, cself, hardi
   real(wp) :: radi, radj, dradi, dradj, dG(3), dS(3, 3), dGc(3), dSc(3, 3), wcni
   real(wp) :: capi, capj
   real(wp), allocatable :: dtrans(:, :)
   real(wp), allocatable :: dtmpdr(:, :), dtmpdL(:, :, :)

   ! Thread-private derivatives of component iat w.r.t. all atoms
   real(wp), allocatable :: dabdr_local(:, :)

   call get_dir_trans(mol, dtrans, cutoff)

   allocate(dtmpdr(3, size(list%nlat, kind=i8)), dtmpdL(3, 3, mol%nat))

   ! CN and effective charge derivative of the electronegativity in the row of
   ! atom iat, without a maximum CN cutoff dCN_i/dR_j = -dCN_j/dR_i and
   ! dqloc_i/dR_j = dqloc_j/dR_i give them from the entries of the row
   !$omp parallel do default(none) schedule(runtime) &
   !$omp shared(self, mol, list, cache, dtmpdr, dtmpdL) &
   !$omp private(iat, izp, ist, ien)
   do iat = 1, mol%nat
      izp = mol%id(iat)
      ist = list%inl(iat)
      ien = list%inl(iat + 1) - 1
      dtmpdr(:, ist) = self%kcnchi(izp) * cache%dcndrlist(:, ist) &
         & + self%kqchi(izp) * cache%dqlocdrlist(:, ist)
      dtmpdr(:, ist+1:ien) = -self%kcnchi(izp) * cache%dcndrlist(:, ist+1:ien) &
         & + self%kqchi(izp) * cache%dqlocdrlist(:, ist+1:ien)
      dtmpdL(:, :, iat) = self%kcnchi(izp) * cache%dcndL(:, :, iat) &
         & + self%kqchi(izp) * cache%dqlocdL(:, :, iat)
   end do
   !$omp end parallel do

   cache%dabdL(:, :, :) = 0.0_wp

   ! Each iteration only updates the derivatives of component iat (row iat)
   !$omp parallel default(none) &
   !$omp shared(self, mol, list, cache, dtrans, alpha, beta, dtmpdr, dtmpdL) &
   !$omp private(iat, izp, jat, jzp, xat, kat, lat, img, gam, vec, dtmp, norm_cn) &
   !$omp private(rvdw, wsw, dgam, wgam, ctmp, cself, hardi, wcni, radi, radj) &
   !$omp private(dradi, dradj, capi, capj, dG, dS, dGc, dSc, dabdr_local)
   allocate(dabdr_local(3, mol%nat))
   !$omp do schedule(runtime)
   do iat = 1, mol%nat
      izp = mol%id(iat)
      dabdr_local(:, :) = 0.0_wp

      ! EN derivative through the capacitance matrix
      do kat = list%inl(iat), list%inl(iat + 1) - 1
         jat = list%nlat(kat)
         do lat = list%inl(jat), list%inl(jat + 1) - 1
            xat = list%nlat(lat)
            dabdr_local(:, xat) = dabdr_local(:, xat) &
               & + beta * cache%clist(kat) * dtmpdr(:, lat)
         end do
         cache%dabdL(:, :, iat) = cache%dabdL(:, :, iat) &
            & + beta * cache%clist(kat) * dtmpdL(:, :, jat)
      end do

      ! Effective charge width of i
      norm_cn = 1.0_wp / self%avg_cn(izp)
      radi = self%rad(izp) * exp(-self%kcnrad(izp) * cache%cn(iat) * norm_cn)
      dradi = -self%kcnrad(izp) * norm_cn * radi
      capi = self%cap(izp)
      ! Effective hardness of i
      hardi = self%eta(izp) + self%kqeta_pre &
      & * tanh(self%kqeta(izp) * cache%qloc(iat)) + sqrt2pi / radi
      wcni = 0.0_wp
      do kat = list%inl(iat) + 1, list%inl(iat + 1) - 1
         if (cache%wsc%nimg_list(kat) == 0) cycle
         jat = list%nlat(kat)
         jzp = mol%id(jat)
         capj = self%cap(jzp)
         rvdw = self%rvdw(izp, jzp)
         ! Effective charge width of j
         norm_cn = 1.0_wp / self%avg_cn(jzp)
         radj = self%rad(jzp) * exp(-self%kcnrad(jzp) * cache%cn(jat) * norm_cn)
         dradj = -self%kcnrad(jzp) * norm_cn * radj

         ! Coulomb interaction of Gaussian charges
         gam = 1.0_wp / sqrt(radi**2 + radj**2)

         wgam = 0.0_wp
         wsw = 1.0_wp / real(cache%wsc%nimg_list(kat), wp)
         do img = cache%wsc%itr_list(kat), cache%wsc%itr_list(kat + 1) - 1
            vec = mol%xyz(:, jat) - mol%xyz(:, iat) + cache%wsc%trans(:, &
            & cache%wsc%tridx_list(img))

            ! Capacitance derivative, dC(i, j)/dR(i) = dGc = -dC(i, j)/dR(j)
            call get_dcpair_dir(self%kbc, vec, dtrans, rvdw, capi, capj, dGc, dSc)
            dGc = dGc * wsw
            dSc = dSc * wsw

            ! EN derivative: capacitance matrix derivative
            dabdr_local(:, iat) = dabdr_local(:, iat) + beta * cache%xtmp(jat) * dGc
            dabdr_local(:, jat) = dabdr_local(:, jat) &
               & - beta * (cache%xtmp(jat) - cache%xtmp(iat)) * dGc
            cache%dabdL(:, :, iat) = cache%dabdL(:, :, iat) &
               & - beta * cache%xtmp(jat) * dSc

            ! Explicit derivative
            call get_damat_dir(vec, dtrans, capi, capj, rvdw, self%kbc, gam, dG, dS, dgam)
            dG = alpha * dG * wsw
            dS = alpha * dS * wsw
            wgam = wgam + alpha * dgam * wsw
            dabdr_local(:, iat) = dabdr_local(:, iat) - dG * cache%vrhs(jat)
            dabdr_local(:, jat) = dabdr_local(:, jat) + dG * cache%vrhs(jat)
            cache%dabdL(:, :, iat) = cache%dabdL(:, :, iat) + dS * cache%vrhs(jat)

            ! Capacitance derivative off-diagonal
            call get_damat_dc_dir(vec, dtrans, capi, capj, rvdw, self%kbc, gam, dG, dS)
            dG = alpha * dG * wsw
            dS = alpha * dS * wsw
            dabdr_local(:, iat) = dabdr_local(:, iat) + dG * cache%vrhs(jat)
            dabdr_local(:, jat) = dabdr_local(:, jat) - dG * cache%vrhs(jat)
            cache%dabdL(:, :, iat) = cache%dabdL(:, :, iat) - dS * cache%vrhs(jat)

            ! Capacitance derivative diagonal
            dabdr_local(:, jat) = dabdr_local(:, jat) &
               & + alpha * hardi * cache%vrhs(iat) * dGc
         end do

         ! Effective charge width derivative, CN of i is collected in wcni
         wcni = wcni + wgam * cache%vrhs(jat) * radi * dradi * gam**3
         dtmp = wgam * cache%vrhs(jat) * radj * dradj * gam**3
         lat = list%inl(jat)
         dabdr_local(:, jat) = dabdr_local(:, jat) + dtmp * cache%dcndrlist(:, lat)
         do lat = list%inl(jat) + 1, list%inl(jat + 1) - 1
            xat = list%nlat(lat)
            dabdr_local(:, xat) = dabdr_local(:, xat) - dtmp * cache%dcndrlist(:, lat)
         end do
         cache%dabdL(:, :, iat) = cache%dabdL(:, :, iat) &
            & + dtmp * cache%dcndL(:, :, jat)
      end do

      ! Self-images, i = j and T != 0
      gam = 1.0_wp / sqrt(2.0_wp * radi**2)
      rvdw = self%rvdw(izp, izp)
      cself = 0.0_wp
      if (cache%wsc%nimg_list(list%inl(iat)) > 0) then
         wgam = 0.0_wp
         wsw = 1.0_wp / real(cache%wsc%nimg_list(list%inl(iat)), wp)
         do img = cache%wsc%itr_list(list%inl(iat)), &
         & cache%wsc%itr_list(list%inl(iat) + 1) - 1
            vec = cache%wsc%trans(:, cache%wsc%tridx_list(img))

            ! Self-image capacitance in the EN
            call get_cpair_dir(self%kbc, vec, dtrans, rvdw, capi, capi, ctmp)
            cself = cself + ctmp * wsw

            ! EN derivative: self-image capacitance derivative
            call get_dcpair_dir(self%kbc, vec, dtrans, rvdw, capi, capi, dGc, dSc)
            cache%dabdL(:, :, iat) = cache%dabdL(:, :, iat) &
               & - beta * cache%xtmp(iat) * dSc * wsw

            ! Explicit derivative
            call get_damat_dir(vec, dtrans, capi, capi, rvdw, self%kbc, gam, dG, dS, dgam)
            cache%dabdL(:, :, iat) = cache%dabdL(:, :, iat) &
               & + alpha * cache%vrhs(iat) * dS * wsw
            wgam = wgam + alpha * dgam * wsw

            ! Capacitance derivative
            call get_damat_dc_dir(vec, dtrans, capi, capi, rvdw, self%kbc, gam, dG, dS)
            cache%dabdL(:, :, iat) = cache%dabdL(:, :, iat) &
               & - alpha * cache%vrhs(iat) * dS * wsw
         end do

         ! Effective charge width derivative
         wcni = wcni + wgam * cache%vrhs(iat) * 2.0_wp * radi * dradi * gam**3
      end if

      ! EN derivative through the self-image capacitance
      do kat = list%inl(iat), list%inl(iat + 1) - 1
         xat = list%nlat(kat)
         dabdr_local(:, xat) = dabdr_local(:, xat) - beta * cself * dtmpdr(:, kat)
      end do
      cache%dabdL(:, :, iat) = cache%dabdL(:, :, iat) - beta * cself * dtmpdL(:, :, iat)

      ! EN derivative: capacitance matrix derivative diagonal
      dabdr_local(:, iat) = dabdr_local(:, iat) &
         & + beta * cache%xtmp(iat) * cache%dcdrdiag(:, iat)
      cache%dabdL(:, :, iat) = cache%dabdL(:, :, iat) &
         & + beta * cache%xtmp(iat) * cache%dcdL(:, :, iat)

      ! Hardness derivative
      dtmp = alpha * self%kqeta_pre * self%kqeta(izp) &
      & / cosh(self%kqeta(izp) * cache%qloc(iat))**2 * cache%vrhs(iat) &
      & * cache%clist(list%inl(iat))
      do kat = list%inl(iat), list%inl(iat + 1) - 1
         xat = list%nlat(kat)
         dabdr_local(:, xat) = dabdr_local(:, xat) + dtmp * cache%dqlocdrlist(:, kat)
      end do
      cache%dabdL(:, :, iat) = cache%dabdL(:, :, iat) &
         & + dtmp * cache%dqlocdL(:, :, iat)

      ! Effective charge width derivative
      dtmp = -alpha * sqrt2pi * dradi / (radi**2) * cache%vrhs(iat) &
      & * cache%clist(list%inl(iat)) + wcni
      kat = list%inl(iat)
      dabdr_local(:, iat) = dabdr_local(:, iat) + dtmp * cache%dcndrlist(:, kat)
      do kat = list%inl(iat) + 1, list%inl(iat + 1) - 1
         xat = list%nlat(kat)
         dabdr_local(:, xat) = dabdr_local(:, xat) - dtmp * cache%dcndrlist(:, kat)
      end do
      cache%dabdL(:, :, iat) = cache%dabdL(:, :, iat) + dtmp * cache%dcndL(:, :, iat)

      ! Capacitance derivative
      dtmp = alpha * hardi * cache%vrhs(iat)
      dabdr_local(:, iat) = dabdr_local(:, iat) + dtmp * cache%dcdrdiag(:, iat)
      cache%dabdL(:, :, iat) = cache%dabdL(:, :, iat) + dtmp * cache%dcdL(:, :, iat)

      ! Derivatives of component iat w.r.t. the atoms in its row
      do kat = list%inl(iat), list%inl(iat + 1) - 1
         cache%dabdrlist(:, kat) = dabdr_local(:, list%nlat(kat))
      end do
   end do
   !$omp end do
   deallocate(dabdr_local)
   !$omp end parallel

end subroutine get_partial_derivs_3d_list

!> Assemble the Coulomb matrix, including bond-capacitance contributions
subroutine get_coulomb_matrix(self, mol, ndim, cache, list)
   !> EEQBC model type
   class(eeqbc_model), intent(in) :: self
   !> Structure type
   type(structure_type), intent(in) :: mol
   !> System size
   integer, intent(in) :: ndim
   !> Multicharge cache
   type(mchrg_cache), intent(inout) :: cache
   !> Multicharge neighborlist type
   type(csr_list), intent(in), optional :: list

   if (present(list)) then
      ! Allocate amat
      if (.not. allocated(cache%alist)) then
         allocate(cache%alist(size(list%nlat, kind=i8)))
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
      norm_cn = cache%cn(iat) / self%avg_cn(izp)
      radi = self%rad(izp) * exp(-self%kcnrad(izp) * norm_cn)
      do jat = 1, iat - 1
         jzp = mol%id(jat)
         vec = mol%xyz(:, jat) - mol%xyz(:, iat)
         r2 = vec(1)**2 + vec(2)**2 + vec(3)**2
         ! Effective charge width of j
         norm_cn = cache%cn(jat) / self%avg_cn(jzp)
         radj = self%rad(jzp) * exp(-self%kcnrad(jzp) * norm_cn)
         ! Coulomb interaction of Gaussian charges
         gam2 = 1.0_wp / (radi**2 + radj**2)
         tmp = erf(sqrt(r2 * gam2)) / sqrt(r2) * cache%cmat(jat, iat)
         amat_local(jat, iat) = tmp
         amat_local(iat, jat) = tmp
      end do
      ! Effective hardness
      tmp = self%eta(izp) + self%kqeta_pre * tanh(self%kqeta(izp) * cache%qloc(iat)) &
      & + sqrt2pi / radi
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
   !> Multicharge neighborlist type
   type(csr_list), intent(in) :: list
   !> Multicharge cache
   type(mchrg_cache), intent(inout) :: cache

   integer(i8) :: kat
   integer :: iat, jat, izp, jzp
   real(wp) :: vec(3), r2, gam2, tmp, norm_cn, radi, radj

   ! Zero out global shared target arrays upfront
   cache%alist(:) = 0.0_wp

   !$omp parallel default(none) &
   !$omp shared(cache, mol, self, list) &
   !$omp private(iat, izp, jat, kat, jzp, gam2, vec, r2, tmp, norm_cn, radi, radj)

   !$omp do schedule(runtime)
   do iat = 1, mol%nat
      izp = mol%id(iat)
      ! Effective charge width of i
      norm_cn = cache%cn(iat) / self%avg_cn(izp)
      radi = self%rad(izp) * exp(-self%kcnrad(izp) * norm_cn)

      do kat = list%inl(iat) + 1, list%inl(iat + 1) - 1
         jat = list%nlat(kat)
         jzp = mol%id(jat)
         vec = mol%xyz(:, jat) - mol%xyz(:, iat)
         r2 = vec(1)**2 + vec(2)**2 + vec(3)**2
         ! Effective charge width of j
         norm_cn = cache%cn(jat) / self%avg_cn(jzp)
         radj = self%rad(jzp) * exp(-self%kcnrad(jzp) * norm_cn)
         ! Coulomb interaction of Gaussian charges
         gam2 = 1.0_wp / (radi**2 + radj**2)
         tmp = erf(sqrt(r2 * gam2)) / sqrt(r2) * cache%clist(kat)
         cache%alist(kat) = tmp
      end do

      ! Effective hardness
      tmp = self%eta(izp) + self%kqeta_pre * tanh(self%kqeta(izp) * cache%qloc(iat)) &
      & + sqrt2pi / radi
      cache%alist(list%inl(iat)) = tmp * cache%clist(list%inl(iat)) + 1.0_wp
   end do
   !$omp end do
   !$omp end parallel

end subroutine get_amat_0d_list

!> Build the Coulomb matrix for a periodic system (3D) using Ewald summation
!> and bond capacitance.
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

   real(wp), allocatable :: amat_local(:, :)

   call get_dir_trans(mol, dtrans, cutoff)

   cache%amat(:, :) = 0.0_wp

   !$omp parallel default(none) &
   !$omp shared(cache, mol, self, dtrans) &
   !$omp private(iat, izp, jat, jzp, gam, vec, dtmp, ctmp, norm_cn) &
   !$omp private(radi, radj, capi, capj, rvdw, r1, wsw, amat_local)
   allocate(amat_local, source=cache%amat)
   !$omp do schedule(runtime)
   do iat = 1, mol%nat
      izp = mol%id(iat)
      ! Effective charge width of i
      norm_cn = cache%cn(iat) / self%avg_cn(izp)
      radi = self%rad(izp) * exp(-self%kcnrad(izp) * norm_cn)
      capi = self%cap(izp)
      do jat = 1, iat - 1
         jzp = mol%id(jat)
         ! Van der Waals distance in Angstrom (approximate factor 2)
         rvdw = self%rvdw(izp, jzp)
         ! Effective charge width of j
         norm_cn = cache%cn(jat) / self%avg_cn(jzp)
         radj = self%rad(jzp) * exp(-self%kcnrad(jzp) * norm_cn)
         capj = self%cap(jzp)
         ! Coulomb interaction of Gaussian charges
         gam = 1.0_wp / sqrt(radi**2 + radj**2)
         wsw = 1.0_wp / real(cache%wsc%nimg(jat, iat), wp)
         do img = 1, cache%wsc%nimg(jat, iat)
            vec = mol%xyz(:, jat) - mol%xyz(:, iat) + cache%wsc%trans(:, &
            & cache%wsc%tridx(img, jat, iat))
            call get_amat_dir_3d(vec, gam, dtrans, self%kbc, rvdw, capi, capj, dtmp)
            amat_local(jat, iat) = amat_local(jat, iat) + dtmp * wsw
            amat_local(iat, jat) = amat_local(iat, jat) + dtmp * wsw
         end do
      end do

      ! Diagonal Coulomb interaction terms
      gam = 1.0_wp / sqrt(2.0_wp * radi**2)
      rvdw = self%rvdw(izp, izp)
      wsw = 1.0_wp / real(cache%wsc%nimg(iat, iat), wp)
      do img = 1, cache%wsc%nimg(iat, iat)
         vec = cache%wsc%trans(:, cache%wsc%tridx(img, iat, iat))
         call get_amat_dir_3d(vec, gam, dtrans, self%kbc, rvdw, capi, capi, dtmp)
         amat_local(iat, iat) = amat_local(iat, iat) + dtmp * wsw
      end do

      ! Effective hardness
      dtmp = self%eta(izp) + self%kqeta_pre &
      & * tanh(self%kqeta(izp) * cache%qloc(iat)) + sqrt2pi / radi
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

!> Build the Coulomb matrix for a periodic system using a CSR adjacency list
subroutine get_amat_3d_list(self, mol, list, cache)
   !> EEQBC model type
   class(eeqbc_model), intent(in) :: self
   !> Molecular structure data
   type(structure_type), intent(in) :: mol
   !> Multicharge neighborlist type
   type(csr_list), intent(in) :: list
   !> Multicharge cache
   type(mchrg_cache), intent(inout) :: cache

   integer :: iat, jat, izp, jzp
   integer(i8) :: kat, img
   real(wp) :: vec(3), gam, dtmp, capi, capj, radi, radj, norm_cn, rvdw, wsw
   real(wp) :: atmp, adiag_tmp
   real(wp), allocatable :: dtrans(:, :)

   call get_dir_trans(mol, dtrans, cutoff)

   ! Zero out global shared target arrays upfront
   cache%alist(:) = 0.0_wp

   !$omp parallel default(none) &
   !$omp shared(cache, mol, self, list, dtrans) &
   !$omp private(iat, izp, kat, jat, jzp, gam, vec, dtmp, norm_cn, radi, radj) &
   !$omp private(capi, capj, rvdw, wsw, img, atmp, adiag_tmp)

   !$omp do schedule(runtime)
   do iat = 1, mol%nat
      izp = mol%id(iat)
      norm_cn = cache%cn(iat) / self%avg_cn(izp)
      radi = self%rad(izp) * exp(-self%kcnrad(izp) * norm_cn)
      capi = self%cap(izp)

      ! Initialize scalar accumulator for the diagonal of this atom
      adiag_tmp = 0.0_wp

      do kat = list%inl(iat) + 1, list%inl(iat+1) - 1
         if (cache%wsc%nimg_list(kat) == 0) cycle

         jat = list%nlat(kat)
         jzp = mol%id(jat)
         capj = self%cap(jzp)
         rvdw = self%rvdw(izp, jzp)
         wsw = 1.0_wp / real(cache%wsc%nimg_list(kat), wp)

         norm_cn = cache%cn(jat) / self%avg_cn(jzp)
         radj = self%rad(jzp) * exp(-self%kcnrad(jzp) * norm_cn)
         gam = 1.0_wp / sqrt(radi**2 + radj**2)

         ! Accumulate image contributions in a local scalar
         atmp = 0.0_wp
         do img = cache%wsc%itr_list(kat), cache%wsc%itr_list(kat+1) - 1
            vec = mol%xyz(:, jat) - mol%xyz(:, iat) + cache%wsc%trans(:, &
            & cache%wsc%tridx_list(img))
            call get_amat_dir_3d(vec, gam, dtrans, self%kbc, rvdw, capi, capj, dtmp)
            atmp = atmp + dtmp * wsw
         end do

         cache%alist(kat) = atmp
      end do

      ! Diagonal Coulomb interaction terms
      gam = 1.0_wp / sqrt(2.0_wp * radi**2)
      rvdw = self%rvdw(izp, izp)
      if (cache%wsc%nimg_list(list%inl(iat)) > 0) then
         wsw = 1.0_wp / real(cache%wsc%nimg_list(list%inl(iat)), wp)
         do img = cache%wsc%itr_list(list%inl(iat)), &
         & cache%wsc%itr_list(list%inl(iat) + 1) - 1
            vec = cache%wsc%trans(:, cache%wsc%tridx_list(img))
            call get_amat_dir_3d(vec, gam, dtrans, self%kbc, rvdw, capi, capi, dtmp)
            adiag_tmp = adiag_tmp + dtmp * wsw
         end do
      end if

      ! Effective hardness
      dtmp = self%eta(izp) + self%kqeta_pre &
      & * tanh(self%kqeta(izp) * cache%qloc(iat)) + sqrt2pi / radi

      ! Single safe direct write (iat is thread-exclusive)
      cache%alist(list%inl(iat)) = adiag_tmp + cache%clist(list%inl(iat)) * dtmp + 1.0_wp
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

!> Contribution to the derivative of the Coulomb matrix from the derivative of
!> the bond capacitance (direct part).
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

   cmat(:, :) = 0.0_wp

   !$omp parallel default(none) &
   !$omp shared(cmat, mol, self) &
   !$omp private(iat, izp, jat, jzp, vec, r1, rvdw, tmp, capi, capj)
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

         cmat(jat, iat) = -tmp
         cmat(iat, jat) = -tmp
      end do
   end do
   !$omp end do

   !$omp do schedule(static)
   do iat = 1, mol%nat
      cmat(iat, iat) = - sum(cmat(iat, :))
   end do
   !$omp end do
   !$omp end parallel

   if (size(cmat, 1) == mol%nat + 1) then
      cmat(mol%nat + 1, mol%nat + 1) = 1.0_wp
   end if

end subroutine get_cmat_0d

!> Build the bond capacitance matrix for a non‑periodic system.
subroutine get_cmat_0d_list(self, mol, list, cache)
   !> EEQBC model type
   class(eeqbc_model), intent(in) :: self
   !> Molecular structure data
   type(structure_type), intent(in) :: mol
   !> CSR list of neighbors
   type(csr_list), intent(in) :: list
   !> EEQBC cache
   type(mchrg_cache), intent(inout) :: cache

   integer :: iat, jat, izp, jzp
   integer(i8) :: kat
   real(wp) :: vec(3), rvdw, tmp, capi, capj, r1
   real(wp) :: diag(mol%nat)

   cache%clist(:) = 0.0_wp
   diag(:) = 0.0_wp

   !$omp parallel default(none) &
   !$omp shared(cache, mol, list, self, diag) &
   !$omp private(iat, kat, izp, jat, jzp, vec, r1, rvdw, tmp, capi, capj)

   !$omp do schedule(runtime) reduction(+:diag)
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

         cache%clist(kat) = -tmp
         diag(iat) = diag(iat) + tmp
         ! A complete list also holds the pair in the row of jat
         if (.not. list%complete) diag(jat) = diag(jat) + tmp
      end do
   end do
   !$omp end do

   !$omp do schedule(static)
   do iat = 1, mol%nat
      cache%clist(list%inl(iat)) = diag(iat)
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
   !> Wigner-Seitz cell
   type(wignerseitz_cell), intent(in) :: wsc
   !> Output capacitance matrix
   real(wp), intent(out) :: cmat(:, :)

   integer :: iat, jat, izp, jzp, img
   real(wp) :: vec(3), rvdw, tmp, capi, capj, wsw
   real(wp), allocatable :: dtrans(:, :)

   call get_dir_trans(mol, dtrans, cutoff)

   cmat(:, :) = 0.0_wp

   !$omp parallel default(none) &
   !$omp shared(cmat, mol, self, wsc, dtrans) &
   !$omp private(iat, izp, jat, jzp, img) &
   !$omp private(vec, rvdw, tmp, capi, capj, wsw)

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
            vec = mol%xyz(:, iat) - mol%xyz(:, jat) - wsc%trans(:, wsc%tridx(img, jat, &
            & iat))
            call get_cpair_dir(self%kbc, vec, dtrans, rvdw, capi, capj, tmp)

            cmat(jat, iat) = cmat(jat, iat) - tmp * wsw
            cmat(iat, jat) = cmat(iat, jat) - tmp * wsw
         end do
      end do

      ! Self-image diagonal term - also race-free, touches only own iat
      rvdw = self%rvdw(izp, izp)
      wsw = 1.0_wp / real(wsc%nimg(iat, iat), wp)
      do img = 1, wsc%nimg(iat, iat)
         vec = wsc%trans(:, wsc%tridx(img, iat, iat))
         call get_cpair_dir(self%kbc, vec, dtrans, rvdw, capi, capi, tmp)
         cmat(iat, iat) = cmat(iat, iat) + tmp * wsw
      end do
   end do
   !$omp end do

   !$omp do schedule(static)
   do iat = 1, mol%nat
      cmat(iat, iat) = cmat(iat, iat) &
         - sum(cmat(iat, 1:iat-1)) &
         - sum(cmat(iat, iat+1:mol%nat))
   end do
   !$omp end do
   !$omp end parallel

   if (size(cmat, 1) == mol%nat + 1) then
      cmat(mol%nat + 1, mol%nat + 1) = 1.0_wp
   end if

end subroutine get_cmat_3d

!> Build the bond capacitance matrix for a periodic system using CSR adjacency list.
subroutine get_cmat_3d_list(self, mol, list, cache)
   !> EEQBC model type
   class(eeqbc_model), intent(in) :: self
   !> Molecular structure data
   type(structure_type), intent(in) :: mol
   !> Multicharge neighborlist type (CSR format)
   type(csr_list), intent(in) :: list
   !> EEQBC cache
   type(mchrg_cache), intent(inout) :: cache

   integer :: iat, jat, izp, jzp
   integer(i8) :: kat, img
   real(wp) :: vec(3), rvdw, tmp, capi, capj, wsw, ctmp
   real(wp), allocatable :: dtrans(:, :)
   real(wp) :: diag(mol%nat)

   call get_dir_trans(mol, dtrans, cutoff)

   cache%clist(:) = 0.0_wp
   diag(:) = 0.0_wp

   !$omp parallel default(none) &
   !$omp shared(cache, mol, list, self, dtrans, diag) &
   !$omp private(iat, izp, jat, kat, jzp, img) &
   !$omp private(vec, rvdw, tmp, capi, capj, wsw, ctmp)

   !$omp do schedule(runtime) reduction(+:diag)
   do iat = 1, mol%nat
      izp = mol%id(iat)
      capi = self%cap(izp)

      ! 1. Off-diagonal neighbor pairs (jat /= iat)
      do kat = list%inl(iat) + 1, list%inl(iat+1) - 1
         if (cache%wsc%nimg_list(kat) == 0) cycle

         jat = list%nlat(kat)
         jzp = mol%id(jat)
         capj = self%cap(jzp)
         rvdw = self%rvdw(izp, jzp)
         wsw = 1.0_wp / real(cache%wsc%nimg_list(kat), wp)

         ctmp = 0.0_wp

         do img = cache%wsc%itr_list(kat), cache%wsc%itr_list(kat+1) - 1
            vec = mol%xyz(:, iat) - mol%xyz(:, jat) - cache%wsc%trans(:, &
            & cache%wsc%tridx_list(img))

            call get_cpair_dir(self%kbc, vec, dtrans, rvdw, capi, capj, tmp)

            ctmp = ctmp - tmp * wsw
            diag(iat) = diag(iat) + tmp * wsw
            ! A complete list also holds the pair in the row of jat
            if (.not. list%complete) diag(jat) = diag(jat) + tmp * wsw
         end do

         cache%clist(kat) = ctmp

      end do

      ! 2. Self-interaction with periodic images R /= 0 (j = iat)
      rvdw = self%rvdw(izp, izp)
      wsw = 1.0_wp / real(cache%wsc%nimg_list(list%inl(iat)), wp)

      do img = cache%wsc%itr_list(list%inl(iat)), &
      & cache%wsc%itr_list(list%inl(iat) + 1) - 1
         vec = cache%wsc%trans(:, cache%wsc%tridx_list(img))
         call get_cpair_dir(self%kbc, vec, dtrans, rvdw, capi, capi, tmp)
         ! Direct write is safe: touches only own atom, not even a cross-thread hazard
         diag(iat) = diag(iat) + tmp * wsw
      end do

   end do
   !$omp end do

   !$omp do schedule(static)
   do iat = 1, mol%nat
      cache%clist(list%inl(iat)) = diag(iat)
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

!> Compute the derivative of the bond capacitance with respect to atomic
!> positions and lattice parameters.
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
   !> Derivative w.r.t. lattice parameters (3×3)
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
   !> Derivative of capacitance matrix w.r.t. lattice parameters (3 × 3 × ndim)
   real(wp), intent(out) :: dcdL(:, :, :)

   integer :: iat, jat, izp, jzp
   real(wp) :: vec(3), rvdw, dG(3), dS(3, 3), capi, capj

   real(wp), allocatable :: dcdL_acc(:, :, :)

   dcdr(:, :, :) = 0.0_wp
   dcdL(:, :, :) = 0.0_wp

   allocate(dcdL_acc(3, 3, mol%nat), source=0.0_wp)

   !$omp parallel default(none) &
   !$omp shared(dcdr, mol, self) &
   !$omp private(iat, izp, jat, jzp, vec, rvdw, dG, dS, capi, capj) &
   !$omp reduction(+:dcdL_acc)

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

         dcdr(:, iat, jat) = +dG
         dcdr(:, jat, iat) = -dG

         dcdL_acc(:, :, iat) = dcdL_acc(:, :, iat) + dS
         dcdL_acc(:, :, jat) = dcdL_acc(:, :, jat) + dS
      end do
   end do
   !$omp end do

   !$omp do schedule(static)
   do iat = 1, mol%nat
      dcdr(:, iat, iat) = sum(dcdr(:, 1:iat-1, iat), dim=2) &
         + sum(dcdr(:, iat+1:mol%nat, iat), dim=2)
   end do
   !$omp end do

   !$omp end parallel

   dcdL(:, :, 1:mol%nat) = dcdL_acc
   deallocate(dcdL_acc)

end subroutine get_dcmat_0d

!> Build the derivative of the bond capacitance matrix for a non‑periodic system.
subroutine get_dcmat_0d_list(self, mol, list, cache)
   !> EEQBC model type
   class(eeqbc_model), intent(in) :: self
   !> Molecular structure data
   type(structure_type), intent(in) :: mol
   !> Multicharge neighborlist type
   type(csr_list), intent(in) :: list
   !> Multicharge cache
   type(mchrg_cache), intent(inout) :: cache

   integer :: iat, jat, izp, jzp, ic, i, j
   integer(i8) :: kat

   real(wp) :: vec(3), rvdw, dG(3), dS(3, 3), capi, capj

   real(wp), allocatable :: dcdrdiag(:, :), dcdL(:, :, :)

   ! Zero out global shared target arrays upfront
   cache%dcdrdiag(:, :) = 0.0_wp
   cache%dcdL(:, :, :) = 0.0_wp

   allocate(dcdrdiag, source=cache%dcdrdiag)
   allocate(dcdL, source=cache%dcdL)

   !$omp parallel default(none) &
   !$omp shared(cache, mol, list, self) &
   !$omp private(iat, izp, jat, kat, jzp, vec, rvdw) &
   !$omp private(dG, dS, capi, capj, ic, i, j) &
   !$omp reduction(+:dcdrdiag, dcdL)

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

         dcdrdiag(:, iat) = dcdrdiag(:, iat) - dG(:)
         dcdL(:, :, iat) = dcdL(:, :, iat) + dS(:, :)

         ! A complete list also holds the pair in the row of jat
         if (.not. list%complete) then
            dcdrdiag(:, jat) = dcdrdiag(:, jat) + dG(:)
            dcdL(:, :, jat) = dcdL(:, :, jat) + dS(:, :)
         end if
      end do
   end do
   !$omp end do
   !$omp end parallel

   cache%dcdrdiag(:, :) = cache%dcdrdiag(:, :) + dcdrdiag(:, :)
   cache%dcdL(:, :, :) = cache%dcdL(:, :, :) + dcdL(:, :, :)

end subroutine get_dcmat_0d_list

!> Build the derivative of the bond capacitance matrix for a periodic system.
subroutine get_dcmat_3d(self, mol, wsc, dcdr, dcdL)
   !> EEQBC model type
   class(eeqbc_model), intent(in) :: self
   !> Molecular structure data
   type(structure_type), intent(in) :: mol
   !> Wigner-Seitz cell
   type(wignerseitz_cell), intent(in) :: wsc
   !> Derivative of capacitance matrix w.r.t. atomic positions (3 × nat × ndim)
   real(wp), intent(out) :: dcdr(:, :, :)
   !> Derivative of capacitance matrix w.r.t. lattice parameters (3 × 3 × ndim)
   real(wp), intent(out) :: dcdL(:, :, :)

   integer :: iat, jat, izp, jzp, img
   real(wp) :: vec(3), rvdw, dG(3), dS(3, 3), capi, capj, wsw
   real(wp), allocatable :: dtrans(:, :)
   real(wp), allocatable :: dcdL_acc(:, :, :)

   dcdr(:, :, :) = 0.0_wp
   dcdL(:, :, :) = 0.0_wp

   allocate(dcdL_acc(3, 3, mol%nat), source=0.0_wp)

   call get_dir_trans(mol, dtrans, cutoff)

   !$omp parallel default(none) &
   !$omp shared(dcdr, mol, self, dtrans, wsc) &
   !$omp private(iat, izp, jat, jzp, vec, rvdw, dG, dS, capi, capj, wsw, img) &
   !$omp reduction(+:dcdL_acc)

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
            vec = mol%xyz(:, jat) - mol%xyz(:, iat) + wsc%trans(:, wsc%tridx(img, jat, &
            & iat))
            call get_dcpair_dir(self%kbc, vec, dtrans, rvdw, capi, capj, dG, dS)

            dcdr(:, iat, jat) = dcdr(:, iat, jat) + dG * wsw
            dcdr(:, jat, iat) = dcdr(:, jat, iat) - dG * wsw

            dcdL_acc(:, :, iat) = dcdL_acc(:, :, iat) + dS * wsw
            dcdL_acc(:, :, jat) = dcdL_acc(:, :, jat) + dS * wsw
         end do
      end do

      rvdw = self%rvdw(izp, izp)
      wsw = 1.0_wp / real(wsc%nimg(iat, iat), wp)
      do img = 1, wsc%nimg(iat, iat)
         vec = wsc%trans(:, wsc%tridx(img, iat, iat))
         call get_dcpair_dir(self%kbc, vec, dtrans, rvdw, capi, capi, dG, dS)
         dcdL_acc(:, :, iat) = dcdL_acc(:, :, iat) + dS * wsw
      end do
   end do
   !$omp end do

   !$omp do schedule(static)
   do iat = 1, mol%nat
      dcdr(:, iat, iat) = sum(dcdr(:, 1:iat-1, iat), dim=2) &
         + sum(dcdr(:, iat+1:mol%nat, iat), dim=2)
   end do
   !$omp end do

   !$omp end parallel

   dcdL(:, :, 1:mol%nat) = dcdL_acc
   deallocate(dcdL_acc)

end subroutine get_dcmat_3d

!> Build the derivative of the bond capacitance matrix for a periodic system using a
!> CSR adjacency list.
subroutine get_dcmat_3d_list(self, mol, list, cache)
   !> EEQBC model type
   class(eeqbc_model), intent(in) :: self
   !> Molecular structure data
   type(structure_type), intent(in) :: mol
   !> Multicharge neighborlist type
   type(csr_list), intent(in) :: list
   !> Multicharge cache
   type(mchrg_cache), intent(inout) :: cache

   integer :: iat, jat, izp, jzp
   integer(i8) :: kat, img
   real(wp) :: vec(3), rvdw, dG(3), dS(3, 3), capi, capj, wsw
   real(wp), allocatable :: dtrans(:, :)
   real(wp), allocatable :: dcdrdiag_acc(:, :), dcdL_acc(:, :, :)

   call get_dir_trans(mol, dtrans, cutoff)

   cache%dcdrdiag(:, :) = 0.0_wp
   cache%dcdL(:, :, :) = 0.0_wp

   allocate(dcdrdiag_acc(3, mol%nat), source=0.0_wp)
   allocate(dcdL_acc(3, 3, mol%nat), source=0.0_wp)

   !$omp parallel default(none) &
   !$omp shared(mol, cache, list, self, dtrans) &
   !$omp private(iat, izp, jat, kat, jzp, vec, rvdw, dG, dS, capi, capj, wsw, img) &
   !$omp reduction(+:dcdrdiag_acc, dcdL_acc)

   !$omp do schedule(runtime)
   do iat = 1, mol%nat
      izp = mol%id(iat)
      capi = self%cap(izp)

      do kat = list%inl(iat) + 1, list%inl(iat+1) - 1
         if (cache%wsc%nimg_list(kat) == 0) cycle

         jat = list%nlat(kat)
         jzp = mol%id(jat)
         capj = self%cap(jzp)
         rvdw = self%rvdw(izp, jzp)
         wsw = 1.0_wp / real(cache%wsc%nimg_list(kat), wp)

         do img = cache%wsc%itr_list(kat), cache%wsc%itr_list(kat+1) - 1
            vec = mol%xyz(:, jat) - mol%xyz(:, iat) + &
               & cache%wsc%trans(:, cache%wsc%tridx_list(img))
            call get_dcpair_dir(self%kbc, vec, dtrans, rvdw, capi, capj, dG, dS)

            dcdrdiag_acc(:, iat) = dcdrdiag_acc(:, iat) - dG(:) * wsw
            dcdL_acc(:, :, iat) = dcdL_acc(:, :, iat) + dS(:, :) * wsw

            ! A complete list also holds the pair in the row of jat
            if (.not. list%complete) then
               dcdrdiag_acc(:, jat) = dcdrdiag_acc(:, jat) + dG(:) * wsw
               dcdL_acc(:, :, jat) = dcdL_acc(:, :, jat) + dS(:, :) * wsw
            end if
         end do
      end do

      rvdw = self%rvdw(izp, izp)
      if (cache%wsc%nimg_list(list%inl(iat)) > 0) then
         wsw = 1.0_wp / real(cache%wsc%nimg_list(list%inl(iat)), wp)
         do img = cache%wsc%itr_list(list%inl(iat)), &
         & cache%wsc%itr_list(list%inl(iat) + 1) - 1
            vec = cache%wsc%trans(:, cache%wsc%tridx_list(img))
            call get_dcpair_dir(self%kbc, vec, dtrans, rvdw, capi, capi, dG, dS)
            dcdL_acc(:, :, iat) = dcdL_acc(:, :, iat) + dS(:, :) * wsw
         end do
      end if

   end do
   !$omp end do
   !$omp end parallel

   cache%dcdrdiag(:, :) = dcdrdiag_acc
   cache%dcdL(:, :, :) = dcdL_acc
   deallocate(dcdrdiag_acc, dcdL_acc)

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
   real(wp) :: r1, dgtmp(3), dstmp(3, 3), vec(3)

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

!> Accumulate the gradient and stress contributions of the EEQBC model
subroutine get_grad(self, mol, cache, p, gradient, sigma, alpha, beta, list)
   !> EEQBC model type
   class(eeqbc_model), intent(in) :: self
   !> Molecular structure data
   type(structure_type), intent(in) :: mol
   !> Multicharge cache
   type(mchrg_cache), intent(in) :: cache
   !> Solution vector
   real(wp), intent(in) :: p(:)
   !> Cartesian gradient
   real(wp), intent(inout) :: gradient(:, :)
   !> Lattice stress contribution
   real(wp), intent(inout) :: sigma(:, :)
   !> Optional gradient prefactor
   real(wp), optional, intent(in) :: alpha
   !> Optional electronegativity prefactor
   real(wp), optional, intent(in) :: beta
   !> Optional neighborlist
   type(csr_list), intent(in), optional :: list

   if (.not. present(list)) then
      if (any(mol%periodic)) then
         call get_grad_3d(self, mol, cache, p, gradient, sigma, alphain=alpha, &
         & betain=beta)
      else
         call get_grad_0d(self, mol, cache, p, gradient, sigma, alphain=alpha, &
         & betain=beta)
      end if
   else
      if (any(mol%periodic)) then
         call get_grad_3d_list(self, mol, list, cache, p, gradient, sigma, &
         & alphain=alpha, betain=beta)
      else
         call get_grad_0d_list(self, mol, list, cache, p, gradient, sigma, &
         & alphain=alpha, betain=beta)
      end if
   end if

end subroutine get_grad

!> Accumulate gradient and stress contributions for a non-periodic CSR system
subroutine get_grad_0d_list(self, mol, list, cache, p, gradient, sigma, alphain, betain)
   !> EEQBC model type
   class(eeqbc_model), intent(in) :: self
   !> Molecular structure data
   type(structure_type), intent(in) :: mol
   !> Multicharge neighborlist type
   type(csr_list), intent(in) :: list
   !> Multicharge cache
   type(mchrg_cache), intent(in) :: cache
   !> Solution vector
   real(wp), intent(in) :: p(:)
   !> Cartesian gradient
   real(wp), intent(inout) :: gradient(:, :)
   !> Lattice stress contribution
   real(wp), intent(inout) :: sigma(:, :)
   !> Optional gradient prefactor
   real(wp), optional, intent(in) :: alphain
   !> Optional electronegativity prefactor
   real(wp), optional, intent(in) :: betain

   real(wp) :: alpha, beta
   integer :: iat, jat, izp, jzp
   integer(i8) :: kat
   real(wp) :: vec(3), r2, gam, arg, dtmp, norm_cn
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
   if (list%complete) then
      call spgemv_csr(mol%nat, cache%clist, list%inl, list%nlat, p, v)
   else
      call spsymv_csr(mol%nat, cache%clist, list%inl, list%nlat, p, v)
   end if

   allocate(gradient_local(3, mol%nat), source=0.0_wp)
   allocate(sigma_local(3, 3), source=0.0_wp)

   !$omp parallel do default(none) schedule(runtime) &
   !$omp shared(cache, mol, self, p, v, list, alpha, beta) &
   !$omp private(iat, jat, izp, jzp, gam, vec, r2, dtmp, norm_cn, arg) &
   !$omp private(radi, radj, dradi, dradj, dG, dS, W_ii, W_jj, W_ij) &
   !$omp reduction(+:cnacc, qlocacc, gradient_local, sigma_local)

   do iat = 1, mol%nat
      izp = mol%id(iat)
      norm_cn = 1.0_wp / self%avg_cn(izp)
      radi = self%rad(izp) * exp(-self%kcnrad(izp) * cache%cn(iat) * norm_cn)
      dradi = -self%kcnrad(izp) * norm_cn * radi

      W_ii = p(iat) * cache%vrhs(iat)

      do kat = list%inl(iat) + 1, list%inl(iat + 1) - 1
         jat = list%nlat(kat)
         ! Each pair once, a complete list also holds it in the row of jat
         if (list%complete .and. jat < iat) cycle
         jzp = mol%id(jat)
         vec = mol%xyz(:, jat) - mol%xyz(:, iat)
         r2 = dot_product(vec, vec)

         W_jj = p(jat) * cache%vrhs(jat)
         W_ij = p(iat) * cache%vrhs(jat) + p(jat) * cache%vrhs(iat)

         norm_cn = 1.0_wp / self%avg_cn(jzp)
         radj = self%rad(jzp) * exp(-self%kcnrad(jzp) * cache%cn(jat) * norm_cn)
         dradj = -self%kcnrad(jzp) * norm_cn * radj

         gam = 1.0_wp / sqrt(radi**2 + radj**2)
         arg = gam * gam * r2
         dtmp = 2.0_wp * exp(-arg) / (sqrtpi)

         ! Accumulate scalar weights for CN derivatives
         cnacc(iat) = cnacc(iat) - (dtmp * radi * dradi * gam**3.0_wp) &
         & * cache%clist(kat) * W_ij * alpha
         cnacc(jat) = cnacc(jat) - (dtmp * radj * dradj * gam**3.0_wp) &
         & * cache%clist(kat) * W_ij * alpha

         ! 1. Explicit Geometry Derivative
         dtmp = 2.0_wp * gam * exp(-arg) / (sqrtpi * r2) - erf(sqrt(arg)) &
         & / (r2 * sqrt(r2))
         dG = dtmp * vec
         dS = spread(dG, 1, 3) * spread(vec, 2, 3)

         gradient_local(:, iat) = gradient_local(:, &
         & iat) - dG * cache%clist(kat) * W_ij * alpha
         gradient_local(:, jat) = gradient_local(:, &
         & jat) + dG * cache%clist(kat) * W_ij * alpha
         sigma_local(:, :)   = sigma_local(:, :)   + dS * cache%clist(kat) * W_ij * alpha

         ! 3 & 4. Capacitance derivatives
         call get_dcpair(self%kbc, vec, self%rvdw(izp, jzp), self%cap(izp), &
         & self%cap(jzp), dG, dS)
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

         dtmp = (self%eta(jzp) + self%kqeta_pre &
         & * tanh(self%kqeta(jzp) * cache%qloc(jat)) + sqrt2pi / radj)
         gradient_local(:, iat) = gradient_local(:, iat) - dtmp * dG * W_jj * alpha
         dtmp = (self%eta(izp) + self%kqeta_pre &
         & * tanh(self%kqeta(izp) * cache%qloc(iat)) + sqrt2pi / radi)
         gradient_local(:, jat) = gradient_local(:, jat) + dtmp * dG * W_ii * alpha
      end do

      ! 5. Diagonal weight accumulation
      qlocacc(iat) = qlocacc(iat) + self%kqeta_pre * self%kqeta(izp) &
         & / cosh(self%kqeta(izp) * cache%qloc(iat))**2 &
         & * W_ii * cache%clist(list%inl(iat)) * alpha &
         & + v(iat) * self%kqchi(mol%id(iat)) * beta
      cnacc(iat) = cnacc(iat) - sqrt2pi * dradi / (radi**2) * W_ii &
      & * cache%clist(list%inl(iat)) * alpha &
         & + v(iat) * self%kcnchi(mol%id(iat)) * beta

      ! 6. Intrinsic capacitance
      dtmp = (self%eta(izp) + self%kqeta_pre &
      & * tanh(self%kqeta(izp) * cache%qloc(iat)) + sqrt2pi / radi) * W_ii
      gradient_local(:, iat) = gradient_local(:, iat) + dtmp * cache%dcdrdiag(:, &
      & iat) * alpha &
         & + p(iat) * cache%xtmp(iat) * cache%dcdrdiag(:, iat) * beta

      sigma_local(:, :)   = sigma_local(:, :)   + dtmp * cache%dcdL(:, :, iat) * alpha &
         & + p(iat) * cache%xtmp(iat) * cache%dcdL(:, :, iat) * beta
   end do
   !$omp end parallel do

   gradient = gradient + gradient_local
   sigma = sigma + sigma_local

   gradient_local = 0.0_wp
   sigma_local = 0.0_wp

   call self%ncoord%add_coordination_number_derivs_list(mol, dtrans, cnacc, &
   & gradient_local, sigma_local, list)
   call self%ncoord_en%add_coordination_number_derivs_list(mol, dtrans, qlocacc, &
   & gradient_local, sigma_local, list)

   gradient = gradient + gradient_local
   sigma = sigma + sigma_local

   deallocate(gradient_local, sigma_local, qlocacc, cnacc, dtrans, v)
end subroutine get_grad_0d_list

!> Accumulate gradient and stress contributions for a non-periodic system
subroutine get_grad_0d(self, mol, cache, p, gradient, sigma, alphain, betain)
   !> EEQBC model type
   class(eeqbc_model), intent(in) :: self
   !> Molecular structure data
   type(structure_type), intent(in) :: mol
   !> Multicharge cache
   type(mchrg_cache), intent(in) :: cache
   !> Solution vector
   real(wp), intent(in) :: p(:)
   !> Cartesian gradient
   real(wp), intent(inout) :: gradient(:, :)
   !> Lattice stress contribution
   real(wp), intent(inout) :: sigma(:, :)
   !> Optional gradient prefactor
   real(wp), optional, intent(in) :: alphain
   !> Optional electronegativity prefactor
   real(wp), optional, intent(in) :: betain

   real(wp) :: alpha, beta
   integer :: iat, jat, izp, jzp
   real(wp) :: vec(3), r2, gam, arg, dtmp, norm_cn
   real(wp) :: radi, radj, dradi, dradj, dG(3), dS(3, 3)
   real(wp) :: W_ii, W_jj, W_ij
   real(wp), allocatable :: gradient_local(:, :), sigma_local(:, :)
   real(wp), allocatable :: v(:), ptmp(:), qlocacc(:), cnacc(:)
   real(wp), allocatable :: dtrans(:, :)

   alpha = 1.0_wp
   if (present(alphain)) alpha = alphain
   beta = 1.0_wp
   if (present(betain)) beta = betain

   allocate(dtrans, source=cache%trans)
   allocate(qlocacc(mol%nat), cnacc(mol%nat))
   qlocacc = 0.0_wp
   cnacc = 0.0_wp

   ! Zero-pad p to the full cmat instead of copying its atomic block, the
   ! constraint row and column of cmat vanish apart from the diagonal
   allocate(ptmp(size(cache%cmat, 1)), source=0.0_wp)
   ptmp(:mol%nat) = p
   allocate(v(size(cache%cmat, 1)))
   call symv(cache%cmat, ptmp, v, alpha=1.0_wp, beta=0.0_wp, uplo='l')

   allocate(gradient_local(3, mol%nat), source=0.0_wp)
   allocate(sigma_local(3, 3), source=0.0_wp)

   !$omp parallel do default(none) schedule(runtime) &
   !$omp shared(cache, mol, self, p, v, alpha, beta) &
   !$omp private(iat, jat, izp, jzp, gam, vec, r2, dtmp, norm_cn, arg) &
   !$omp private(radi, radj, dradi, dradj, dG, dS, W_ii, W_jj, W_ij) &
   !$omp reduction(+:cnacc, qlocacc, gradient_local, sigma_local)

   do iat = 1, mol%nat
      izp = mol%id(iat)
      norm_cn = 1.0_wp / self%avg_cn(izp)
      radi = self%rad(izp) * exp(-self%kcnrad(izp) * cache%cn(iat) * norm_cn)
      dradi = -self%kcnrad(izp) * norm_cn * radi

      W_ii = p(iat) * cache%vrhs(iat)

      do jat = 1, iat - 1
         jzp = mol%id(jat)
         vec = mol%xyz(:, jat) - mol%xyz(:, iat)
         r2 = dot_product(vec, vec)

         W_jj = p(jat) * cache%vrhs(jat)
         W_ij = p(iat) * cache%vrhs(jat) + p(jat) * cache%vrhs(iat)

         norm_cn = 1.0_wp / self%avg_cn(jzp)
         radj = self%rad(jzp) * exp(-self%kcnrad(jzp) * cache%cn(jat) * norm_cn)
         dradj = -self%kcnrad(jzp) * norm_cn * radj

         gam = 1.0_wp / sqrt(radi**2 + radj**2)
         arg = gam * gam * r2
         dtmp = 2.0_wp * exp(-arg) / (sqrtpi)

         ! Accumulate scalar weights for CN derivatives
         cnacc(iat) = cnacc(iat) - (dtmp * radi * dradi * gam**3.0_wp) &
         & * cache%cmat(iat, jat) * W_ij * alpha
         cnacc(jat) = cnacc(jat) - (dtmp * radj * dradj * gam**3.0_wp) &
         & * cache%cmat(iat, jat) * W_ij * alpha

         ! 1. Explicit Geometry Derivative
         dtmp = 2.0_wp * gam * exp(-arg) / (sqrtpi * r2) - erf(sqrt(arg)) &
         & / (r2 * sqrt(r2))
         dG = dtmp * vec
         dS = spread(dG, 1, 3) * spread(vec, 2, 3)

         gradient_local(:, iat) = gradient_local(:, iat) - dG * cache%cmat(iat, &
         & jat) * W_ij * alpha
         gradient_local(:, jat) = gradient_local(:, jat) + dG * cache%cmat(iat, &
         & jat) * W_ij * alpha
         sigma_local(:, :)   = sigma_local(:, :)   + dS * cache%cmat(iat, &
         & jat) * W_ij * alpha

         ! 3 & 4. Capacitance derivatives
         call get_dcpair(self%kbc, vec, self%rvdw(izp, jzp), self%cap(izp), &
         & self%cap(jzp), dG, dS)
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

         dtmp = (self%eta(jzp) + self%kqeta_pre &
         & * tanh(self%kqeta(jzp) * cache%qloc(jat)) + sqrt2pi / radj)
         gradient_local(:, iat) = gradient_local(:, iat) - dtmp * dG * W_jj * alpha
         dtmp = (self%eta(izp) + self%kqeta_pre &
         & * tanh(self%kqeta(izp) * cache%qloc(iat)) + sqrt2pi / radi)
         gradient_local(:, jat) = gradient_local(:, jat) + dtmp * dG * W_ii * alpha
      end do

      ! 5. Diagonal weight accumulation
      qlocacc(iat) = qlocacc(iat) + self%kqeta_pre * self%kqeta(izp) &
         & / cosh(self%kqeta(izp) * cache%qloc(iat))**2 &
         & * W_ii * cache%cmat(iat, iat) * alpha &
         & + v(iat) * self%kqchi(mol%id(iat)) * beta
      cnacc(iat) = cnacc(iat) - sqrt2pi * dradi / (radi**2) * W_ii * cache%cmat(iat, &
      & iat) * alpha &
         & + v(iat) * self%kcnchi(mol%id(iat)) * beta

      ! 6. Intrinsic capacitance
      dtmp = (self%eta(izp) + self%kqeta_pre &
      & * tanh(self%kqeta(izp) * cache%qloc(iat)) + sqrt2pi / radi) * W_ii
      gradient_local(:, iat) = gradient_local(:, iat) + dtmp * cache%dcdr(:, iat, &
      & iat) * alpha &
         & + p(iat) * cache%xtmp(iat) * cache%dcdr(:, iat, iat) * beta

      sigma_local(:, :)   = sigma_local(:, :)   + dtmp * cache%dcdL(:, :, iat) * alpha &
         & + p(iat) * cache%xtmp(iat) * cache%dcdL(:, :, iat) * beta
   end do
   !$omp end parallel do

   gradient = gradient + gradient_local
   sigma = sigma + sigma_local

   gradient_local = 0.0_wp
   sigma_local = 0.0_wp

   call self%ncoord%add_coordination_number_derivs(mol, dtrans, cnacc, gradient_local, &
   & sigma_local)
   call self%ncoord_en%add_coordination_number_derivs(mol, dtrans, qlocacc, &
   & gradient_local, sigma_local)

   gradient = gradient + gradient_local
   sigma = sigma + sigma_local

   deallocate(gradient_local, sigma_local, qlocacc, cnacc, dtrans, v, ptmp)
end subroutine get_grad_0d

!> Accumulate gradient and stress contributions for a periodic CSR system
subroutine get_grad_3d_list(self, mol, list, cache, p, gradient, sigma, alphain, betain)
   !> EEQBC model type
   class(eeqbc_model), intent(in) :: self
   !> Molecular structure data
   type(structure_type), intent(in) :: mol
   !> Multicharge neighborlist type
   type(csr_list), intent(in) :: list
   !> Multicharge cache
   type(mchrg_cache), intent(in) :: cache
   !> Solution vector
   real(wp), intent(in) :: p(:)
   !> Cartesian gradient
   real(wp), intent(inout) :: gradient(:, :)
   !> Lattice stress contribution
   real(wp), intent(inout) :: sigma(:, :)
   !> Optional gradient prefactor
   real(wp), optional, intent(in) :: alphain
   !> Optional electronegativity prefactor
   real(wp), optional, intent(in) :: betain

   real(wp) :: alpha, beta
   integer :: iat, jat, izp, jzp
   integer(i8) :: kat, img
   real(wp) :: vec(3), r2, gam, dtmp, capi, capj, dgam, rvdw, wsw
   real(wp) :: radi, radj, dradi, dradj, dG(3), dS(3, 3), ctmp
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

   allocate(gradient_local(3, mol%nat), source=0.0_wp)
   allocate(sigma_local(3, 3), source=0.0_wp)

   allocate(v(mol%nat))
   if (list%complete) then
      call spgemv_csr(mol%nat, cache%clist, list%inl, list%nlat, p, v)
   else
      call spsymv_csr(mol%nat, cache%clist, list%inl, list%nlat, p, v)
   end if

   !$omp parallel do default(none) schedule(runtime) &
   !$omp shared(cache, mol, list, self, p, v, dtrans, alpha, beta) &
   !$omp private(iat, kat, izp, jat, jzp, gam, vec, r2) &
   !$omp private(dtmp, radi, radj, dradi, dradj, dG, dS, norm_cn) &
   !$omp private(W_ii, W_jj, W_ij, wsw, capi, capj, rvdw, img, dgam, ctmp) &
   !$omp reduction(+:cnacc, qlocacc, gradient_local, sigma_local)
   do iat = 1, mol%nat
      izp = mol%id(iat)
      norm_cn = 1.0_wp / self%avg_cn(izp)
      radi = self%rad(izp) * exp(-self%kcnrad(izp) * cache%cn(iat) * norm_cn)
      dradi = -self%kcnrad(izp) * norm_cn * radi
      capi = self%cap(izp)
      W_ii = p(iat) * cache%vrhs(iat)

      do kat = list%inl(iat) + 1, list%inl(iat + 1) - 1
         jat = list%nlat(kat)
         ! Each pair once, a complete list also holds it in the row of jat
         if (list%complete .and. jat < iat) cycle
         jzp = mol%id(jat)
         capj = self%cap(jzp)
         rvdw = self%rvdw(izp, jzp)
         if (cache%wsc%nimg_list(kat) == 0) cycle
         wsw = 1.0_wp / real(cache%wsc%nimg_list(kat), wp)
         W_jj = p(jat) * cache%vrhs(jat)
         W_ij = p(iat) * cache%vrhs(jat) + p(jat) * cache%vrhs(iat)

         norm_cn = 1.0_wp / self%avg_cn(jzp)
         radj = self%rad(jzp) * exp(-self%kcnrad(jzp) * cache%cn(jat) * norm_cn)
         dradj = -self%kcnrad(jzp) * norm_cn * radj
         gam = 1.0_wp / sqrt(radi**2 + radj**2)

         do img = cache%wsc%itr_list(kat), cache%wsc%itr_list(kat+1) - 1
            vec = mol%xyz(:, jat) - mol%xyz(:, iat) + cache%wsc%trans(:, &
            & cache%wsc%tridx_list(img))
            call get_damat_dir(vec, dtrans, capi, capj, rvdw, self%kbc, gam, dG, dS, dgam)

            cnacc(iat) = cnacc(iat) + (dgam * wsw * radi * dradi * gam**3.0_wp) &
            & * W_ij * alpha
            cnacc(jat) = cnacc(jat) + (dgam * wsw * radj * dradj * gam**3.0_wp) &
            & * W_ij * alpha

            gradient_local(:, iat) = gradient_local(:, iat) - dG * wsw * W_ij * alpha
            gradient_local(:, jat) = gradient_local(:, jat) + dG * wsw * W_ij * alpha
            sigma_local(:, :) = sigma_local(:, :) + dS * wsw * W_ij * alpha

            call get_damat_dc_dir(vec, dtrans, capi, capj, rvdw, self%kbc, gam, dG, dS)
            gradient_local(:, iat) = gradient_local(:, iat) + dG * wsw * W_ij * alpha
            gradient_local(:, jat) = gradient_local(:, jat) - dG * wsw * W_ij * alpha
            sigma_local(:, :)   = sigma_local(:, :) - W_ij * dS * wsw * alpha

            call get_dcpair_dir(self%kbc, vec, dtrans, rvdw, capi, capj, dG, dS)
            dtmp = (self%eta(jzp) + self%kqeta_pre &
            & * tanh(self%kqeta(jzp) * cache%qloc(jat)) + sqrt2pi / radj)
            gradient_local(:, iat) = gradient_local(:, &
            & iat) - dtmp * dG * wsw * W_jj * alpha &
               & + p(iat) * cache%xtmp(jat) * dG * wsw * beta &
               & - p(jat) * (cache%xtmp(jat) - cache%xtmp(iat)) * dG * wsw * beta

            dtmp = (self%eta(izp) + self%kqeta_pre &
            & * tanh(self%kqeta(izp) * cache%qloc(iat)) + sqrt2pi / radi)
            gradient_local(:, jat) = gradient_local(:, &
            & jat) + dtmp * dG * wsw * W_ii * alpha &
               & - p(jat) * cache%xtmp(iat) * dG * wsw * beta &
               & + p(iat) * (cache%xtmp(iat) - cache%xtmp(jat)) * dG * wsw * beta

            ! Lattice Sigma Updates
            ! This ensures the (xi - xj) term is correctly formed
            ! against the diagonal dcdL
            sigma_local(:, :) = sigma_local(:, &
            & :) - p(iat) * cache%xtmp(jat) * dS * wsw * beta
            sigma_local(:, :) = sigma_local(:, &
            & :) - p(jat) * cache%xtmp(iat) * dS * wsw * beta

         end do
      end do

      ! Diagonal corrections
      gam = 1.0_wp / sqrt(2.0_wp * radi**2)
      rvdw = self%rvdw(izp, izp)

      if (cache%wsc%nimg_list(list%inl(iat)) > 0) then
         wsw = 1.0_wp / real(cache%wsc%nimg_list(list%inl(iat)), wp)

         do img = cache%wsc%itr_list(list%inl(iat)), &
         & cache%wsc%itr_list(list%inl(iat) + 1) - 1
            vec = cache%wsc%trans(:, cache%wsc%tridx_list(img))
            call get_damat_dir(vec, dtrans, capi, capi, rvdw, self%kbc, gam, dG, dS, dgam)
            sigma_local(:, :) = sigma_local(:, :) + dS * wsw * W_ii * alpha
            cnacc(iat) = cnacc(iat) + (dgam * wsw * 2.0_wp * radi * dradi &
            & * gam**3.0_wp) * W_ii * alpha

            call get_damat_dc_dir(vec, dtrans, capi, capi, rvdw, self%kbc, gam, dG, dS)
            sigma_local(:, :) = sigma_local(:, :) - W_ii * dS * wsw * alpha

            call get_dcpair_dir(self%kbc, vec, dtrans, rvdw, capi, capi, dG, dS)
            sigma_local(:, :) = sigma_local(:, &
            & :) - p(iat) * cache%xtmp(iat) * dS * wsw * beta

            call get_cpair_dir(self%kbc, vec, dtrans, rvdw, capi, capi, ctmp)
            ctmp = ctmp * wsw

            v(iat) = v(iat) - ctmp * p(iat)
         end do
      end if

      qlocacc(iat) = qlocacc(iat) + self%kqeta_pre * self%kqeta(izp) / &
         & cosh(self%kqeta(izp) * cache%qloc(iat))**2 * W_ii &
         & * cache%clist(list%inl(iat)) * alpha &
         & + v(iat) * self%kqchi(izp) * beta
      cnacc(iat) = cnacc(iat) - sqrt2pi * dradi / (radi**2) * W_ii &
      & * cache%clist(list%inl(iat)) * alpha &
         & + v(iat) * self%kcnchi(izp) * beta

      dtmp = (self%eta(izp) + self%kqeta_pre &
      & * tanh(self%kqeta(izp) * cache%qloc(iat)) + sqrt2pi / radi) * W_ii
      gradient_local(:, iat) = gradient_local(:, iat) + dtmp * cache%dcdrdiag(:, &
      & iat) * alpha &
         & + p(iat) * cache%xtmp(iat) * cache%dcdrdiag(:, iat) * beta
      sigma_local(:, :) = sigma_local(:, :) + dtmp * cache%dcdL(:, :, iat) * alpha &
         & + p(iat) * cache%xtmp(iat) * cache%dcdL(:, :, iat) * beta
   end do
   !$omp end parallel do

   gradient = gradient + gradient_local
   sigma = sigma + sigma_local

   gradient_local = 0.0_wp
   sigma_local = 0.0_wp

   call self%ncoord%add_coordination_number_derivs_list(mol, dtrans, cnacc, &
   & gradient_local, sigma_local, list)
   call self%ncoord_en%add_coordination_number_derivs_list(mol, dtrans, qlocacc, &
   & gradient_local, sigma_local, list)

   gradient = gradient + gradient_local
   sigma = sigma + sigma_local

   deallocate(gradient_local, sigma_local, qlocacc, cnacc, dtrans, v)
end subroutine get_grad_3d_list

!> Accumulate gradient and stress contributions for a periodic system
subroutine get_grad_3d(self, mol, cache, p, gradient, sigma, alphain, betain)
   !> EEQBC model type
   class(eeqbc_model), intent(in) :: self
   !> Molecular structure data
   type(structure_type), intent(in) :: mol
   !> Multicharge cache
   type(mchrg_cache), intent(in) :: cache
   !> Solution vector
   real(wp), intent(in) :: p(:)
   !> Cartesian gradient
   real(wp), intent(inout) :: gradient(:, :)
   !> Lattice stress contribution
   real(wp), intent(inout) :: sigma(:, :)
   !> Optional gradient prefactor
   real(wp), optional, intent(in) :: alphain
   !> Optional electronegativity prefactor
   real(wp), optional, intent(in) :: betain

   real(wp) :: alpha, beta
   integer :: iat, jat, izp, jzp, img
   real(wp) :: vec(3), gam, dtmp, capi, capj, dgam, rvdw, wsw, ctmp
   real(wp) :: radi, radj, dradi, dradj, dG(3), dS(3, 3), norm_cn
   real(wp) :: W_ii, W_jj, W_ij
   real(wp), allocatable :: gradient_local(:, :), sigma_local(:, :)
   real(wp), allocatable :: qlocacc(:), cnacc(:), dtrans(:, :), v(:), ptmp(:)

   alpha = 1.0_wp
   if (present(alphain)) alpha = alphain
   beta = 1.0_wp
   if (present(betain)) beta = betain

   call get_dir_trans(mol, dtrans, cutoff)

   allocate(qlocacc(mol%nat), cnacc(mol%nat))
   qlocacc = 0.0_wp
   cnacc = 0.0_wp

   ! Zero-pad p to the full cmat instead of copying its atomic block, the
   ! constraint row and column of cmat vanish apart from the diagonal
   allocate(ptmp(size(cache%cmat, 1)), source=0.0_wp)
   ptmp(:mol%nat) = p
   allocate(v(size(cache%cmat, 1)))
   call symv(cache%cmat, ptmp, v, alpha=1.0_wp, beta=0.0_wp, uplo='l')

   allocate(gradient_local(3, mol%nat), source=0.0_wp)
   allocate(sigma_local(3, 3), source=0.0_wp)

   !$omp parallel do default(none) schedule(runtime) &
   !$omp shared(cache, mol, self, p, dtrans, v, alpha, beta) &
   !$omp private(iat, izp, jat, jzp, gam, vec, dtmp) &
   !$omp private(radi, radj, dradi, dradj, dG, dS, norm_cn) &
   !$omp private(W_ii, W_jj, W_ij, wsw, capi, capj, rvdw, img, dgam, ctmp) &
   !$omp reduction(+:cnacc, qlocacc, gradient_local, sigma_local)
   do iat = 1, mol%nat
      izp = mol%id(iat)
      norm_cn = 1.0_wp / self%avg_cn(izp)
      radi = self%rad(izp) * exp(-self%kcnrad(izp) * cache%cn(iat) * norm_cn)
      dradi = -self%kcnrad(izp) * norm_cn * radi
      capi = self%cap(izp)
      W_ii = p(iat) * cache%vrhs(iat)

      do jat = 1, iat - 1
         jzp = mol%id(jat)
         capj = self%cap(jzp)
         rvdw = self%rvdw(izp, jzp)
         wsw = 1.0_wp / real(cache%wsc%nimg(iat, jat), wp)
         W_jj = p(jat) * cache%vrhs(jat)
         W_ij = p(iat) * cache%vrhs(jat) + p(jat) * cache%vrhs(iat)

         norm_cn = 1.0_wp / self%avg_cn(jzp)
         radj = self%rad(jzp) * exp(-self%kcnrad(jzp) * cache%cn(jat) * norm_cn)
         dradj = -self%kcnrad(jzp) * norm_cn * radj
         gam = 1.0_wp / sqrt(radi**2 + radj**2)

         do img = 1, cache%wsc%nimg(iat, jat)
            vec = mol%xyz(:, jat) - mol%xyz(:, iat) + cache%wsc%trans(:, &
            & cache%wsc%tridx(img, jat, iat))

            ! 1. A-matrix explicit kernel derivative
            call get_damat_dir(vec, dtrans, capi, capj, rvdw, self%kbc, gam, dG, dS, dgam)

            cnacc(iat) = cnacc(iat) + (dgam * wsw * radi * dradi * gam**3.0_wp) &
            & * W_ij * alpha
            cnacc(jat) = cnacc(jat) + (dgam * wsw * radj * dradj * gam**3.0_wp) &
            & * W_ij * alpha

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

            dtmp = (self%eta(jzp) + self%kqeta_pre &
            & * tanh(self%kqeta(jzp) * cache%qloc(jat)) + sqrt2pi / radj)
            gradient_local(:, iat) = gradient_local(:, &
            & iat) - dtmp * dG * wsw * W_jj * alpha &
               & + p(iat) * cache%xtmp(jat) * dG * wsw * beta &
               & - p(jat) * (cache%xtmp(jat) - cache%xtmp(iat)) * dG * wsw * beta

            dtmp = (self%eta(izp) + self%kqeta_pre &
            & * tanh(self%kqeta(izp) * cache%qloc(iat)) + sqrt2pi / radi)
            gradient_local(:, jat) = gradient_local(:, &
            & jat) + dtmp * dG * wsw * W_ii * alpha &
               & - p(jat) * cache%xtmp(iat) * dG * wsw * beta &
               & + p(iat) * (cache%xtmp(iat) - cache%xtmp(jat)) * dG * wsw * beta

            sigma_local(:, :)   = sigma_local(:, &
            & :)   - p(iat) * cache%xtmp(jat) * dS * wsw * beta &
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
         cnacc(iat) = cnacc(iat) + (dgam * wsw * 2.0_wp * radi * dradi &
         & * gam**3.0_wp) * W_ii * alpha

         call get_damat_dc_dir(vec, dtrans, capi, capi, rvdw, self%kbc, gam, dG, dS)
         sigma_local(:, :) = sigma_local(:, :) - W_ii * dS * wsw * alpha

         call get_dcpair_dir(self%kbc, vec, dtrans, rvdw, capi, capi, dG, dS)
         sigma_local(:, :) = sigma_local(:, &
         & :) - p(iat) * cache%xtmp(iat) * dS * wsw * beta

         ! Apply self-image correction to potential v(iat) before weight calculation
         call get_cpair_dir(self%kbc, vec, dtrans, rvdw, capi, capi, ctmp)
         v(iat) = v(iat) - ctmp * wsw * p(iat)
      end do

      ! Diagonal contributions to weights
      qlocacc(iat) = qlocacc(iat) + self%kqeta_pre * self%kqeta(izp) / &
         & cosh(self%kqeta(izp) * cache%qloc(iat))**2 * W_ii * cache%cmat(iat, &
         & iat) * alpha &
         & + v(iat) * self%kqchi(izp) * beta
      cnacc(iat) = cnacc(iat) - sqrt2pi * dradi / (radi**2) * W_ii * cache%cmat(iat, &
      & iat) * alpha &
         & + v(iat) * self%kcnchi(izp) * beta

      ! Intrinsic capacitance
      dtmp = (self%eta(izp) + self%kqeta_pre &
      & * tanh(self%kqeta(izp) * cache%qloc(iat)) + sqrt2pi / radi) * W_ii
      gradient_local(:, iat) = gradient_local(:, iat) + dtmp * cache%dcdr(:, iat, &
      & iat) * alpha &
         & + p(iat) * cache%xtmp(iat) * cache%dcdr(:, iat, iat) * beta

      sigma_local(:, :)   = sigma_local(:, :)   + dtmp * cache%dcdL(:, :, iat) * alpha &
         & + p(iat) * cache%xtmp(iat) * cache%dcdL(:, :, iat) * beta
   end do
   !$omp end parallel do

   gradient = gradient + gradient_local
   sigma = sigma + sigma_local

   gradient_local = 0.0_wp
   sigma_local = 0.0_wp

   call self%ncoord%add_coordination_number_derivs(mol, dtrans, cnacc, gradient_local, &
   & sigma_local)
   call self%ncoord_en%add_coordination_number_derivs(mol, dtrans, qlocacc, &
   & gradient_local, sigma_local)

   gradient = gradient + gradient_local
   sigma = sigma + sigma_local

   deallocate(gradient_local, sigma_local, qlocacc, cnacc, dtrans, v, ptmp)
end subroutine get_grad_3d

!> Compute the second derivative of the bond capacitance with respect to the
!> distance vector between the atoms.
subroutine get_d2cpair(kbc, vec, rvdw, capi, capj, d2cpair)
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
   !> Second derivative w.r.t. the distance vector (3 × 3)
   real(wp), intent(out) :: d2cpair(3, 3)

   real(wp) :: r1, arg, dtmp, dtmp1, r2
   integer :: ic, jc

   d2cpair(:, :) = 0.0_wp

   r1 = norm2(vec)
   r2 = r1 * r1
   ! Capacitance of bond between atom i and j
   arg = -(kbc * (r1 - rvdw) / rvdw)**2
   dtmp = -sqrt(capi * capj) * kbc * exp(arg) / (sqrtpi * rvdw)
   dtmp1 = - 2.0_wp * kbc * kbc * (r1 - rvdw) / (rvdw * rvdw)

   do ic = 1, 3
      do jc = 1, 3
         d2cpair(jc, ic) = dtmp * (dtmp1 - 1.0_wp / r1) * vec(ic) * vec(jc) / r2
      end do
      d2cpair(ic, ic) = d2cpair(ic, ic) + dtmp / r1
   end do

end subroutine get_d2cpair

!> Build the diagonal blocks d2C_ii/dR_i dR_i of the second derivative of the
!> bond capacitance matrix for a non‑periodic system.
subroutine get_d2cmat_0d(self, mol, cache)
   !> EEQBC model type
   class(eeqbc_model), intent(in) :: self
   !> Molecular structure data
   type(structure_type), intent(in) :: mol
   !> EEQBC cache
   type(mchrg_cache), intent(inout) :: cache

   integer :: iat, jat, izp, jzp, npair, ipair
   real(wp) :: vec(3), rvdw, d2cpair(3, 3), capi, capj
   real(wp), allocatable :: diagonal_local(:, :, :)

   cache%d2cdr2diag(:, :, :) = 0.0_wp
   npair = mol%nat*(mol%nat - 1)/2

   !$omp parallel default(none) &
   !$omp shared(self, mol, cache, npair) &
   !$omp private(ipair, iat, jat, izp, jzp, vec) &
   !$omp private(rvdw, d2cpair, capi, capj, diagonal_local)
   allocate(diagonal_local(3, 3, mol%nat), source=0.0_wp)
   !$omp do schedule(runtime)
   do ipair = 1, npair
      iat = int(0.5_wp*(1.0_wp + sqrt(8.0_wp*real(ipair, wp) + 1.0_wp)))
      if (iat*(iat - 1)/2 < ipair) iat = iat + 1
      jat = ipair - (iat - 1)*(iat - 2)/2
      izp = mol%id(iat)
      jzp = mol%id(jat)
      capi = self%cap(izp)
      capj = self%cap(jzp)
      rvdw = self%rvdw(izp, jzp)
      vec = mol%xyz(:, jat) - mol%xyz(:, iat)

      call get_d2cpair(self%kbc, vec, rvdw, capi, capj, d2cpair)
      diagonal_local(:, :, iat) = diagonal_local(:, :, iat) + d2cpair
      diagonal_local(:, :, jat) = diagonal_local(:, :, jat) + d2cpair

   end do
   !$omp end do

   !$omp critical
   cache%d2cdr2diag(:, :, :) = cache%d2cdr2diag(:, :, :) + &
      & diagonal_local(:, :, :)
   !$omp end critical
   deallocate(diagonal_local)
   !$omp end parallel

end subroutine get_d2cmat_0d

!> Build the diagonal blocks d2C_ii/dR_i dR_i of the second derivative of the
!> bond capacitance matrix for a non‑periodic system using CSR adjacency list.
subroutine get_d2cmat_0d_list(self, mol, list, cache)
   !> EEQBC model type
   class(eeqbc_model), intent(in) :: self
   !> Molecular structure data
   type(structure_type), intent(in) :: mol
   !> CSR neighbor list
   type(csr_list), intent(in) :: list
   !> EEQBC cache
   type(mchrg_cache), intent(inout) :: cache

   integer :: iat, jat, izp, jzp
   integer(i8) :: kat
   real(wp) :: vec(3), rvdw, d2cpair(3, 3), capi, capj
   real(wp), allocatable :: diagonal_local(:, :, :)
   logical :: halflist

   cache%d2cdr2diag(:, :, :) = 0.0_wp
   halflist = .not. list%complete

   !$omp parallel default(none) &
   !$omp shared(self, mol, cache, halflist, list) &
   !$omp private(iat, jat, kat, izp, jzp, vec) &
   !$omp private(rvdw, d2cpair, capi, capj, diagonal_local)
   allocate(diagonal_local(3, 3, mol%nat), source=0.0_wp)
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

         call get_d2cpair(self%kbc, vec, rvdw, capi, capj, d2cpair)
         diagonal_local(:, :, iat) = diagonal_local(:, :, iat) + d2cpair
         if (halflist) then
            diagonal_local(:, :, jat) = diagonal_local(:, :, jat) + d2cpair
         end if
      end do
   end do
   !$omp end do

   !$omp critical
   cache%d2cdr2diag(:, :, :) = cache%d2cdr2diag(:, :, :) + &
      & diagonal_local(:, :, :)
   !$omp end critical
   deallocate(diagonal_local)
   !$omp end parallel

end subroutine get_d2cmat_0d_list

!> Compute the explicit part of the molecular nuclear Hessian,
!> alpha * p^T d2A/dR2 q + beta * p^T d2X/dR2 at fixed p and charges q. The
!> cache has to hold the position derivatives of the CN and the local charges,
!> computed by the update for the Hessian. The diagonal blocks of the second
!> derivative of the capacitance matrix are computed into the cache.
subroutine get_explicit_hessian_dense(self, mol, cache, p, hessian, error, alpha, &
   & beta)
   !> EEQBC model type
   class(eeqbc_model), intent(in) :: self
   !> Molecular structure data
   type(structure_type), intent(in) :: mol
   !> Multicharge cache
   type(mchrg_cache), intent(inout) :: cache
   !> Charge-like contraction vector
   real(wp), intent(in) :: p(:)
   !> Molecular nuclear Hessian
   real(wp), intent(inout) :: hessian(:, :)
   !> Error handling
   type(error_type), allocatable, intent(out) :: error
   !> Optional Coulomb matrix prefactor
   real(wp), intent(in), optional :: alpha
   !> Optional electronegativity prefactor
   real(wp), intent(in), optional :: beta

   if (.not. allocated(cache%d2cdr2diag)) allocate(cache%d2cdr2diag(3, 3, mol%nat))
   call self%get_d2cmat_0d(mol, cache)
   call get_explicit_hessian_0d(self, mol, cache, p, hessian, alphain=alpha, &
      & betain=beta)

end subroutine get_explicit_hessian_dense

!> Compute the explicit part of the molecular nuclear Hessian at fixed p and
!> charges q on a complete neighborlist. The cache has to hold the position
!> derivatives of the CN and the local charges on the neighborlist, computed by
!> the update for the Hessian. The blocks are added to the compressed Hessian
!> on the pattern of its own list, blocks outside of the pattern are dropped.
subroutine get_explicit_hessian_list(self, mol, cache, p, hessian, error, list, &
   & alpha, beta)
   !> EEQBC model type
   class(eeqbc_model), intent(in) :: self
   !> Molecular structure data
   type(structure_type), intent(in) :: mol
   !> Multicharge cache
   type(mchrg_cache), intent(inout) :: cache
   !> Charge-like contraction vector
   real(wp), intent(in) :: p(:)
   !> Molecular nuclear Hessian in compressed storage with allocated blocks
   type(mchrg_hessian), intent(inout) :: hessian
   !> Error handling
   type(error_type), allocatable, intent(out) :: error
   !> Multicharge neighborlist type (complete)
   type(csr_list), intent(in) :: list
   !> Optional Coulomb matrix prefactor
   real(wp), intent(in), optional :: alpha
   !> Optional electronegativity prefactor
   real(wp), intent(in), optional :: beta

   if (.not. allocated(cache%d2cdr2diag)) allocate(cache%d2cdr2diag(3, 3, mol%nat))
   call self%get_d2cmat_0d_list(mol, list, cache)
   call get_explicit_hessian_0d_list(self, mol, list, cache, p, hessian, &
      & alphain=alpha, betain=beta)

end subroutine get_explicit_hessian_list

!> Compute the explicit part of the molecular nuclear Hessian,
!> alpha * p^T d2A/dR2 q + beta * p^T d2X/dR2 at fixed p and charges q
subroutine get_explicit_hessian_0d(self, mol, cache, p, hessian, alphain, betain)
   !> EEQBC model type
   class(eeqbc_model), intent(in) :: self

   !> Molecular structure data
   type(structure_type), intent(in) :: mol

   !> Cache for intermediate quantities
   type(mchrg_cache), intent(in) :: cache

   !> Charge-like contraction vector
   real(wp), intent(in) :: p(:)

   !> Molecular nuclear Hessian
   real(wp), intent(inout) :: hessian(:,:)

   !> Alpha parameter
   real(wp), intent(in), optional :: alphain

   !> Beta parameter
   real(wp), intent(in), optional :: betain

   real(wp) :: alpha, beta
   integer :: iat, jat, kat, izp, jzp, ic, jc, ii, jj
   real(wp) :: vec(3), evec(3), r1, r2, gam, arg, dtmp, norm_cn_i, norm_cn_j
   real(wp) :: radi, radj, dradi, dradj, d2radi, d2radj, capi, capj, rvdw
   real(wp) :: tqeta, hardi, hardj, dhdcni, dhdcnj, d2hdcni2, detai, detaj, d2etai
   real(wp) :: W_ii, W_jj, W_ij, wcapi, wcapj, pchi, cij
   real(wp) :: gtmp, dgdgam, dgdr, d2gdr2, d2gdrdgam, d2gdgam2
   real(wp) :: dgamdcni, dgamdcnj, d2gamdcni2, d2gamdcnj2, d2gamdcnij
   real(wp) :: kii, kij, kjj, dvec(3), dcpair(3), dS(3, 3)
   real(wp) :: d2cpair(3, 3), d2J(3, 3), d2pair(3, 3)
   real(wp), allocatable :: hessian_local(:, :)
   real(wp), allocatable :: v(:), ptmp(:), qlocacc(:), cnacc(:)
   real(wp), allocatable :: dcnacc(:, :, :), dqlocacc(:, :, :)
   real(wp), allocatable :: dtrans(:, :)

   alpha = 1.0_wp
   if (present(alphain)) alpha = alphain
   beta = 1.0_wp
   if (present(betain)) beta = betain

   allocate(dtrans, source=cache%trans)
   allocate(qlocacc(mol%nat), cnacc(mol%nat))
   qlocacc = 0.0_wp
   cnacc = 0.0_wp

   ! Zero-pad p to the full cmat instead of copying its atomic block, the
   ! constraint row and column of cmat vanish apart from the diagonal
   allocate(ptmp(size(cache%cmat, 1)), source=0.0_wp)
   ptmp(:mol%nat) = p
   allocate(v(size(cache%cmat, 1)))
   call symv(cache%cmat, ptmp, v, alpha=1.0_wp, beta=0.0_wp, uplo='l')

   allocate(hessian_local(3 * mol%nat, 3 * mol%nat), source=0.0_wp)

   ! Vector weights of the CN and local charge derivatives in the products of
   ! first derivatives, dcnacc(:, :, iat) multiplies dCN_iat/dR
   allocate(dcnacc(3, mol%nat, mol%nat), dqlocacc(3, mol%nat, mol%nat), &
      & source=0.0_wp)

   !$omp parallel do default(none) schedule(runtime) &
   !$omp shared(cache, mol, self, p, v, alpha, beta) &
   !$omp private(iat, jat, izp, jzp, ic, ii, jj, vec, evec, r1, r2, gam, arg) &
   !$omp private(dtmp, norm_cn_i, norm_cn_j, radi, radj, dradi, dradj, d2radi) &
   !$omp private(d2radj, capi, capj, rvdw, tqeta, hardi, hardj, dhdcni, dhdcnj) &
   !$omp private(d2hdcni2, detai, detaj, d2etai, W_ii, W_jj, W_ij, wcapi, wcapj) &
   !$omp private(pchi, cij, gtmp, dgdgam, dgdr, d2gdr2, d2gdrdgam, d2gdgam2) &
   !$omp private(dgamdcni, dgamdcnj, d2gamdcni2, d2gamdcnj2, d2gamdcnij, kii) &
   !$omp private(kij, kjj, dvec, dcpair, dS, d2cpair, d2J, d2pair) &
   !$omp reduction(+:cnacc, qlocacc, dcnacc, dqlocacc, hessian_local)

   do iat = 1, mol%nat
      izp = mol%id(iat)
      ii = 3 * (iat - 1)
      norm_cn_i = 1.0_wp / self%avg_cn(izp)
      radi = self%rad(izp) * exp(-self%kcnrad(izp) * cache%cn(iat) * norm_cn_i)
      dradi = -self%kcnrad(izp) * norm_cn_i * radi
      d2radi = -self%kcnrad(izp) * norm_cn_i * dradi
      capi = self%cap(izp)

      ! Effective hardness of i and its derivatives w.r.t. CN and local charge
      tqeta = tanh(self%kqeta(izp) * cache%qloc(iat))
      hardi = self%eta(izp) + self%kqeta_pre * tqeta + sqrt2pi / radi
      dhdcni = -sqrt2pi * dradi / radi**2
      d2hdcni2 = sqrt2pi * (2.0_wp * dradi**2 / radi**3 - d2radi / radi**2)
      detai = self%kqeta_pre * self%kqeta(izp) * (1.0_wp - tqeta**2)
      d2etai = -2.0_wp * self%kqeta(izp) * tqeta * detai

      W_ii = p(iat) * cache%vrhs(iat)

      ! Weight of the diagonal capacitance C_ii
      wcapi = alpha * W_ii * hardi + beta * p(iat) * cache%xtmp(iat)

      do jat = 1, iat - 1
         jzp = mol%id(jat)
         jj = 3 * (jat - 1)
         vec = mol%xyz(:, jat) - mol%xyz(:, iat)
         r2 = dot_product(vec, vec)
         r1 = sqrt(r2)
         evec = vec / r1
         norm_cn_j = 1.0_wp / self%avg_cn(jzp)
         radj = self%rad(jzp) * exp(-self%kcnrad(jzp) * cache%cn(jat) * norm_cn_j)
         dradj = -self%kcnrad(jzp) * norm_cn_j * radj
         d2radj = -self%kcnrad(jzp) * norm_cn_j * dradj
         capj = self%cap(jzp)
         rvdw = self%rvdw(izp, jzp)

         ! Effective hardness of j and its derivatives w.r.t. CN and local charge
         tqeta = tanh(self%kqeta(jzp) * cache%qloc(jat))
         hardj = self%eta(jzp) + self%kqeta_pre * tqeta + sqrt2pi / radj
         dhdcnj = -sqrt2pi * dradj / radj**2
         detaj = self%kqeta_pre * self%kqeta(jzp) * (1.0_wp - tqeta**2)

         W_jj = p(jat) * cache%vrhs(jat)
         W_ij = p(iat) * cache%vrhs(jat) + p(jat) * cache%vrhs(iat)

         ! Weights of the diagonal capacitance C_jj and of the off-diagonal
         ! capacitance C_ij in the electronegativity term
         wcapj = alpha * W_jj * hardj + beta * p(jat) * cache%xtmp(jat)
         pchi = beta * (p(iat) * cache%xtmp(jat) + p(jat) * cache%xtmp(iat))
         cij = cache%cmat(iat, jat)

         ! Coulomb interaction of Gaussian charges g = erf(gam*r)/r and its
         ! derivatives w.r.t. the distance and gam
         gam = 1.0_wp / sqrt(radi**2 + radj**2)
         arg = gam * gam * r2
         gtmp = erf(gam * r1) / r1
         dgdgam = 2.0_wp * exp(-arg) / sqrtpi
         dgdr = (gam * dgdgam - gtmp) / r1
         d2gdr2 = -2.0_wp * (dgdr / r1 + gam**3 * dgdgam)
         d2gdrdgam = -2.0_wp * gam**2 * r1 * dgdgam
         d2gdgam2 = -2.0_wp * gam * r2 * dgdgam

         ! Derivatives of gam w.r.t. the CNs of i and j
         dgamdcni = -radi * dradi * gam**3
         dgamdcnj = -radj * dradj * gam**3
         d2gamdcni2 = 3.0_wp * (radi * dradi)**2 * gam**5 &
            & - (dradi**2 + radi * d2radi) * gam**3
         d2gamdcnj2 = 3.0_wp * (radj * dradj)**2 * gam**5 &
            & - (dradj**2 + radj * d2radj) * gam**3
         d2gamdcnij = 3.0_wp * radi * dradi * radj * dradj * gam**5

         ! Bond capacitance derivatives w.r.t. R_j, dC_ij/dR_j = -dcpair
         call get_dcpair(self%kbc, vec, rvdw, capi, capj, dcpair, dS)
         call get_d2cpair(self%kbc, vec, rvdw, capi, capj, d2cpair)

         ! Second derivative of g w.r.t. R_j at fixed CNs
         do ic = 1, 3
            d2J(:, ic) = (d2gdr2 - dgdr / r1) * evec * evec(ic)
            d2J(ic, ic) = d2J(ic, ic) + dgdr / r1
         end do

         ! Second derivative of the pair w.r.t. R_j at fixed CNs and local charges
         ! without the diagonal capacitances,
         ! d2(C_ij g) = d2C_ij g + dC_ij dg^T + dg dC_ij^T + C_ij d2g
         do ic = 1, 3
            d2pair(:, ic) = alpha * W_ij * (-gtmp * d2cpair(:, ic) &
               & - dgdr * (dcpair * evec(ic) + evec * dcpair(ic)) &
               & + cij * d2J(:, ic)) - pchi * d2cpair(:, ic)
         end do

         ! The diagonal capacitances contribute to all blocks of the pair apart
         ! from their own diagonal block, which is added from d2cdr2diag
         hessian_local(ii+1:ii+3, ii+1:ii+3) = hessian_local(ii+1:ii+3, ii+1:ii+3) &
            & + d2pair + wcapj * d2cpair
         hessian_local(jj+1:jj+3, jj+1:jj+3) = hessian_local(jj+1:jj+3, jj+1:jj+3) &
            & + d2pair + wcapi * d2cpair
         hessian_local(ii+1:ii+3, jj+1:jj+3) = hessian_local(ii+1:ii+3, jj+1:jj+3) &
            & - d2pair - (wcapi + wcapj) * d2cpair
         hessian_local(jj+1:jj+3, ii+1:ii+3) = hessian_local(jj+1:jj+3, ii+1:ii+3) &
            & - d2pair - (wcapi + wcapj) * d2cpair

         ! Accumulate scalar weights for the CN second derivatives
         cnacc(iat) = cnacc(iat) + alpha * W_ij * cij * dgdgam * dgamdcni
         cnacc(jat) = cnacc(jat) + alpha * W_ij * cij * dgdgam * dgamdcnj

         ! Mixed derivatives of the Coulomb interaction w.r.t. the positions and
         ! the CNs, dvec is the derivative w.r.t. R_j and its negative w.r.t. R_i
         dvec = alpha * W_ij * (cij * d2gdrdgam * evec - dgdgam * dcpair)
         dcnacc(:, jat, iat) = dcnacc(:, jat, iat) + dgamdcni * dvec
         dcnacc(:, iat, iat) = dcnacc(:, iat, iat) - dgamdcni * dvec
         dcnacc(:, jat, jat) = dcnacc(:, jat, jat) + dgamdcnj * dvec
         dcnacc(:, iat, jat) = dcnacc(:, iat, jat) - dgamdcnj * dvec

         ! Mixed derivatives of the capacitance weights w.r.t. the positions and
         ! the CNs or local charges, dc_ij/dR_j = dcpair
         dtmp = alpha * W_ii * dhdcni + beta * (p(iat) - p(jat)) * self%kcnchi(izp)
         dcnacc(:, jat, iat) = dcnacc(:, jat, iat) + dtmp * dcpair
         dcnacc(:, iat, iat) = dcnacc(:, iat, iat) - dtmp * dcpair
         dtmp = alpha * W_jj * dhdcnj + beta * (p(jat) - p(iat)) * self%kcnchi(jzp)
         dcnacc(:, jat, jat) = dcnacc(:, jat, jat) + dtmp * dcpair
         dcnacc(:, iat, jat) = dcnacc(:, iat, jat) - dtmp * dcpair
         dtmp = alpha * W_ii * detai + beta * (p(iat) - p(jat)) * self%kqchi(izp)
         dqlocacc(:, jat, iat) = dqlocacc(:, jat, iat) + dtmp * dcpair
         dqlocacc(:, iat, iat) = dqlocacc(:, iat, iat) - dtmp * dcpair
         dtmp = alpha * W_jj * detaj + beta * (p(jat) - p(iat)) * self%kqchi(jzp)
         dqlocacc(:, jat, jat) = dqlocacc(:, jat, jat) + dtmp * dcpair
         dqlocacc(:, iat, jat) = dqlocacc(:, iat, jat) - dtmp * dcpair

         ! Products of the CN derivatives of i and j from the Coulomb interaction,
         ! each accumulator takes half of the symmetric product
         kii = 0.5_wp * alpha * W_ij * cij &
            & * (d2gdgam2 * dgamdcni**2 + dgdgam * d2gamdcni2)
         kij = 0.5_wp * alpha * W_ij * cij &
            & * (d2gdgam2 * dgamdcni * dgamdcnj + dgdgam * d2gamdcnij)
         kjj = 0.5_wp * alpha * W_ij * cij &
            & * (d2gdgam2 * dgamdcnj**2 + dgdgam * d2gamdcnj2)
         dcnacc(:, :, iat) = dcnacc(:, :, iat) + kii * cache%dcndr(:, :, iat) &
            & + kij * cache%dcndr(:, :, jat)
         dcnacc(:, :, jat) = dcnacc(:, :, jat) + kij * cache%dcndr(:, :, iat) &
            & + kjj * cache%dcndr(:, :, jat)
      end do

      ! Own diagonal block of the diagonal capacitance C_ii
      hessian_local(ii+1:ii+3, ii+1:ii+3) = hessian_local(ii+1:ii+3, ii+1:ii+3) &
         & + wcapi * cache%d2cdr2diag(:, :, iat)

      ! Scalar weights for the CN and local charge second derivatives from the
      ! hardness and the electronegativity
      cnacc(iat) = cnacc(iat) + alpha * W_ii * cache%cmat(iat, iat) * dhdcni &
         & + beta * self%kcnchi(izp) * v(iat)
      qlocacc(iat) = qlocacc(iat) + alpha * W_ii * cache%cmat(iat, iat) * detai &
         & + beta * self%kqchi(izp) * v(iat)

      ! Products of the CN and local charge derivatives from the hardness
      dcnacc(:, :, iat) = dcnacc(:, :, iat) + 0.5_wp * alpha * W_ii &
         & * cache%cmat(iat, iat) * d2hdcni2 * cache%dcndr(:, :, iat)
      dqlocacc(:, :, iat) = dqlocacc(:, :, iat) + 0.5_wp * alpha * W_ii &
         & * cache%cmat(iat, iat) * d2etai * cache%dqlocdr(:, :, iat)
   end do
   !$omp end parallel do

   ! Products of the CN and local charge derivatives with their vector weights,
   ! each iteration only updates the columns of atom jat
   !$omp parallel do default(none) schedule(runtime) &
   !$omp shared(mol, cache, dcnacc, dqlocacc, hessian_local) &
   !$omp private(iat, jat, kat, jc, ii, jj)
   do jat = 1, mol%nat
      jj = 3 * (jat - 1)
      do jc = 1, 3
         do kat = 1, mol%nat
            do iat = 1, mol%nat
               ii = 3 * (iat - 1)
               hessian_local(ii+1:ii+3, jj+jc) = hessian_local(ii+1:ii+3, jj+jc) &
                  & + cache%dcndr(:, iat, kat) * dcnacc(jc, jat, kat) &
                  & + dcnacc(:, iat, kat) * cache%dcndr(jc, jat, kat) &
                  & + cache%dqlocdr(:, iat, kat) * dqlocacc(jc, jat, kat) &
                  & + dqlocacc(:, iat, kat) * cache%dqlocdr(jc, jat, kat)
            end do
         end do
      end do
   end do
   !$omp end parallel do

   hessian = hessian + hessian_local
   hessian_local = 0.0_wp

   call self%ncoord%add_coordination_number_hessian(mol, dtrans, cnacc, hessian_local)
   call self%ncoord_en%add_coordination_number_hessian(mol, dtrans, qlocacc, &
      & hessian_local)

   hessian = hessian + hessian_local

   deallocate(hessian_local, qlocacc, cnacc, dcnacc, dqlocacc, dtrans, v, ptmp)

end subroutine get_explicit_hessian_0d

!> Compute the explicit part of the molecular nuclear Hessian,
!> alpha * p^T d2A/dR2 q + beta * p^T d2X/dR2 at fixed p and charges q, using a
!> complete neighborlist. Each pair is visited from the rows of both atoms, a
!> row only updates its own row of the compressed Hessian and its own entries
!> of the weights. Blocks outside of the pattern of the list of the Hessian are
!> dropped, and the products with CN or local charge derivatives below cnthr
!> are skipped. Without a maximum CN cutoff the pair derivatives give
!> dCN_i/dR_j = -dCN_j/dR_i and dqloc_i/dR_j = dqloc_j/dR_i, so the derivatives
!> of the CN and local charge of atom k are read from the row of atom k.
subroutine get_explicit_hessian_0d_list(self, mol, list, cache, p, hessian, &
   & alphain, betain)
   !> EEQBC model type
   class(eeqbc_model), intent(in) :: self

   !> Molecular structure data
   type(structure_type), intent(in) :: mol

   !> Multicharge neighborlist type (complete)
   type(csr_list), intent(in) :: list

   !> Cache for intermediate quantities
   type(mchrg_cache), intent(in) :: cache

   !> Charge-like contraction vector
   real(wp), intent(in) :: p(:)

   !> Molecular nuclear Hessian in compressed storage with allocated blocks
   type(mchrg_hessian), intent(inout) :: hessian

   !> Alpha parameter
   real(wp), intent(in), optional :: alphain

   !> Beta parameter
   real(wp), intent(in), optional :: betain

   real(wp) :: alpha, beta
   integer :: iat, jat, izp, jzp, ic, it, ntouch
   integer(i8) :: kat, lat, ist, hat, cat
   real(wp) :: vec(3), evec(3), r1, r2, gam, arg, dtmp, norm_cn_i, norm_cn_j
   real(wp) :: radi, radj, dradi, dradj, d2radi, capi, capj, rvdw
   real(wp) :: tqeta, hardi, hardj, dhdcni, dhdcnj, d2hdcni2, detai, detaj, d2etai
   real(wp) :: W_ii, W_jj, W_ij, wcapi, wcapj, pchi, cij
   real(wp) :: gtmp, dgdgam, dgdr, d2gdr2, d2gdrdgam, d2gdgam2
   real(wp) :: dgamdcni, dgamdcnj, d2gamdcni2, d2gamdcnij
   real(wp) :: dvec(3), dcpair(3), dS(3, 3)
   real(wp) :: d2cpair(3, 3), d2J(3, 3), d2pair(3, 3)
   real(wp) :: dcnj(3), dqlocj(3), dcnx(3), gcni(3), gqloci(3)
   integer, allocatable :: touch(:)
   integer(i8), allocatable :: hpos(:), cninl(:), cnidx(:)
   logical, allocatable :: mark(:)
   real(wp), allocatable :: v(:), qlocacc(:), cnacc(:), kcnacc(:), kqlocacc(:)
   real(wp), allocatable :: dcnacc(:, :, :), dqlocacc(:, :, :)
   real(wp), allocatable :: gcn(:, :), gqloc(:, :)

   alpha = 1.0_wp
   if (present(alphain)) alpha = alphain
   beta = 1.0_wp
   if (present(betain)) beta = betain

   allocate(qlocacc(mol%nat), cnacc(mol%nat))
   qlocacc = 0.0_wp
   cnacc = 0.0_wp

   allocate(v(mol%nat))
   call spgemv_csr(mol%nat, cache%clist, list%inl, list%nlat, p, v)

   ! Vector weights of the CN and local charge derivatives on the list, for the
   ! entry of atom j in the row of atom i dcnacc(:, 1, kat) multiplies
   ! dCN_i/dR_j and dcnacc(:, 2, kat) multiplies dCN_j/dR_i. The products of
   ! the CN derivatives are kept as the second derivatives w.r.t. the CNs on
   ! the list in kcnacc and the local charges in kqlocacc.
   allocate(dcnacc(3, 2, size(list%nlat, kind=i8)), &
      & dqlocacc(3, 2, size(list%nlat, kind=i8)), kcnacc(size(list%nlat, kind=i8)), &
      & kqlocacc(mol%nat), source=0.0_wp)

   ! Each iteration only updates row iat of the list and of the Hessian, hpos
   ! holds the position of each atom in the row of the Hessian
   !$omp parallel default(none) &
   !$omp shared(cache, mol, self, list, p, v, alpha, beta, hessian) &
   !$omp shared(cnacc, qlocacc, dcnacc, dqlocacc, kcnacc, kqlocacc) &
   !$omp private(iat, jat, izp, jzp, ic, kat, ist, hat, vec, evec, r1, r2) &
   !$omp private(gam, arg, dtmp, norm_cn_i, norm_cn_j, radi, radj, dradi, dradj) &
   !$omp private(d2radi, capi, capj, rvdw, tqeta, hardi, hardj, dhdcni, dhdcnj) &
   !$omp private(d2hdcni2, detai, detaj, d2etai, W_ii, W_jj, W_ij, wcapi, wcapj) &
   !$omp private(pchi, cij, gtmp, dgdgam, dgdr, d2gdr2, d2gdrdgam, d2gdgam2) &
   !$omp private(dgamdcni, dgamdcnj, d2gamdcni2, d2gamdcnij, dvec, dcpair, dS) &
   !$omp private(d2cpair, d2J, d2pair, hpos)
   allocate(hpos(mol%nat), source=0_i8)
   !$omp do schedule(runtime)
   do iat = 1, mol%nat
      izp = mol%id(iat)
      ist = list%inl(iat)
      do hat = hessian%list%inl(iat), hessian%list%inl(iat + 1) - 1
         hpos(hessian%list%nlat(hat)) = hat
      end do
      norm_cn_i = 1.0_wp / self%avg_cn(izp)
      radi = self%rad(izp) * exp(-self%kcnrad(izp) * cache%cn(iat) * norm_cn_i)
      dradi = -self%kcnrad(izp) * norm_cn_i * radi
      d2radi = -self%kcnrad(izp) * norm_cn_i * dradi
      capi = self%cap(izp)

      ! Effective hardness of i and its derivatives w.r.t. CN and local charge
      tqeta = tanh(self%kqeta(izp) * cache%qloc(iat))
      hardi = self%eta(izp) + self%kqeta_pre * tqeta + sqrt2pi / radi
      dhdcni = -sqrt2pi * dradi / radi**2
      d2hdcni2 = sqrt2pi * (2.0_wp * dradi**2 / radi**3 - d2radi / radi**2)
      detai = self%kqeta_pre * self%kqeta(izp) * (1.0_wp - tqeta**2)
      d2etai = -2.0_wp * self%kqeta(izp) * tqeta * detai

      W_ii = p(iat) * cache%vrhs(iat)

      ! Weight of the diagonal capacitance C_ii
      wcapi = alpha * W_ii * hardi + beta * p(iat) * cache%xtmp(iat)

      do kat = ist + 1, list%inl(iat + 1) - 1
         jat = list%nlat(kat)
         jzp = mol%id(jat)
         vec = mol%xyz(:, jat) - mol%xyz(:, iat)
         r2 = dot_product(vec, vec)
         r1 = sqrt(r2)
         evec = vec / r1
         norm_cn_j = 1.0_wp / self%avg_cn(jzp)
         radj = self%rad(jzp) * exp(-self%kcnrad(jzp) * cache%cn(jat) * norm_cn_j)
         dradj = -self%kcnrad(jzp) * norm_cn_j * radj
         capj = self%cap(jzp)
         rvdw = self%rvdw(izp, jzp)

         ! Effective hardness of j and its derivatives w.r.t. CN and local charge
         tqeta = tanh(self%kqeta(jzp) * cache%qloc(jat))
         hardj = self%eta(jzp) + self%kqeta_pre * tqeta + sqrt2pi / radj
         dhdcnj = -sqrt2pi * dradj / radj**2
         detaj = self%kqeta_pre * self%kqeta(jzp) * (1.0_wp - tqeta**2)

         W_jj = p(jat) * cache%vrhs(jat)
         W_ij = p(iat) * cache%vrhs(jat) + p(jat) * cache%vrhs(iat)

         ! Weights of the diagonal capacitance C_jj and of the off-diagonal
         ! capacitance C_ij in the electronegativity term
         wcapj = alpha * W_jj * hardj + beta * p(jat) * cache%xtmp(jat)
         pchi = beta * (p(iat) * cache%xtmp(jat) + p(jat) * cache%xtmp(iat))
         cij = cache%clist(kat)

         ! Coulomb interaction of Gaussian charges g = erf(gam*r)/r and its
         ! derivatives w.r.t. the distance and gam
         gam = 1.0_wp / sqrt(radi**2 + radj**2)
         arg = gam * gam * r2
         gtmp = erf(gam * r1) / r1
         dgdgam = 2.0_wp * exp(-arg) / sqrtpi
         dgdr = (gam * dgdgam - gtmp) / r1
         d2gdr2 = -2.0_wp * (dgdr / r1 + gam**3 * dgdgam)
         d2gdrdgam = -2.0_wp * gam**2 * r1 * dgdgam
         d2gdgam2 = -2.0_wp * gam * r2 * dgdgam

         ! Derivatives of gam w.r.t. the CNs of i and j
         dgamdcni = -radi * dradi * gam**3
         dgamdcnj = -radj * dradj * gam**3
         d2gamdcni2 = 3.0_wp * (radi * dradi)**2 * gam**5 &
            & - (dradi**2 + radi * d2radi) * gam**3
         d2gamdcnij = 3.0_wp * radi * dradi * radj * dradj * gam**5

         ! Bond capacitance derivatives w.r.t. R_j, dC_ij/dR_j = -dcpair
         call get_dcpair(self%kbc, vec, rvdw, capi, capj, dcpair, dS)
         call get_d2cpair(self%kbc, vec, rvdw, capi, capj, d2cpair)

         ! Second derivative of g w.r.t. R_j at fixed CNs
         do ic = 1, 3
            d2J(:, ic) = (d2gdr2 - dgdr / r1) * evec * evec(ic)
            d2J(ic, ic) = d2J(ic, ic) + dgdr / r1
         end do

         ! Second derivative of the pair w.r.t. R_j at fixed CNs and local charges
         ! without the diagonal capacitances,
         ! d2(C_ij g) = d2C_ij g + dC_ij dg^T + dg dC_ij^T + C_ij d2g
         do ic = 1, 3
            d2pair(:, ic) = alpha * W_ij * (-gtmp * d2cpair(:, ic) &
               & - dgdr * (dcpair * evec(ic) + evec * dcpair(ic)) &
               & + cij * d2J(:, ic)) - pchi * d2cpair(:, ic)
         end do

         ! The diagonal capacitances contribute to all blocks of the pair apart
         ! from their own diagonal block, which is added from d2cdr2diag, the
         ! row of atom j is updated from the row of atom j
         hat = hessian%list%inl(iat)
         hessian%hesslist(:, :, hat) = hessian%hesslist(:, :, hat) &
            & + d2pair + wcapj * d2cpair
         hat = hpos(jat)
         if (hat > 0) hessian%hesslist(:, :, hat) = hessian%hesslist(:, :, hat) &
            & - d2pair - (wcapi + wcapj) * d2cpair

         ! Accumulate scalar weights for the CN second derivatives
         cnacc(iat) = cnacc(iat) + alpha * W_ij * cij * dgdgam * dgamdcni

         ! Mixed derivatives of the Coulomb interaction w.r.t. the positions and
         ! the CNs, dvec is the derivative w.r.t. R_j and its negative w.r.t. R_i
         dvec = alpha * W_ij * (cij * d2gdrdgam * evec - dgdgam * dcpair)
         dcnacc(:, 1, kat) = dcnacc(:, 1, kat) + dgamdcni * dvec
         dcnacc(:, 1, ist) = dcnacc(:, 1, ist) - dgamdcni * dvec
         dcnacc(:, 2, kat) = dcnacc(:, 2, kat) - dgamdcnj * dvec

         ! Mixed derivatives of the capacitance weights w.r.t. the positions and
         ! the CNs or local charges, dc_ij/dR_j = dcpair
         dtmp = alpha * W_ii * dhdcni + beta * (p(iat) - p(jat)) * self%kcnchi(izp)
         dcnacc(:, 1, kat) = dcnacc(:, 1, kat) + dtmp * dcpair
         dcnacc(:, 1, ist) = dcnacc(:, 1, ist) - dtmp * dcpair
         dtmp = alpha * W_jj * dhdcnj + beta * (p(jat) - p(iat)) * self%kcnchi(jzp)
         dcnacc(:, 2, kat) = dcnacc(:, 2, kat) - dtmp * dcpair
         dtmp = alpha * W_ii * detai + beta * (p(iat) - p(jat)) * self%kqchi(izp)
         dqlocacc(:, 1, kat) = dqlocacc(:, 1, kat) + dtmp * dcpair
         dqlocacc(:, 1, ist) = dqlocacc(:, 1, ist) - dtmp * dcpair
         dtmp = alpha * W_jj * detaj + beta * (p(jat) - p(iat)) * self%kqchi(jzp)
         dqlocacc(:, 2, kat) = dqlocacc(:, 2, kat) - dtmp * dcpair

         ! Second derivatives w.r.t. the CNs of i and j from the Coulomb
         ! interaction
         kcnacc(ist) = kcnacc(ist) + alpha * W_ij * cij &
            & * (d2gdgam2 * dgamdcni**2 + dgdgam * d2gamdcni2)
         kcnacc(kat) = alpha * W_ij * cij &
            & * (d2gdgam2 * dgamdcni * dgamdcnj + dgdgam * d2gamdcnij)
      end do

      ! Own diagonal block of the diagonal capacitance C_ii
      hat = hessian%list%inl(iat)
      hessian%hesslist(:, :, hat) = hessian%hesslist(:, :, hat) &
         & + wcapi * cache%d2cdr2diag(:, :, iat)
      do hat = hessian%list%inl(iat), hessian%list%inl(iat + 1) - 1
         hpos(hessian%list%nlat(hat)) = 0
      end do

      ! Scalar weights for the CN and local charge second derivatives from the
      ! hardness and the electronegativity
      cnacc(iat) = cnacc(iat) + alpha * W_ii * cache%clist(ist) * dhdcni &
         & + beta * self%kcnchi(izp) * v(iat)
      qlocacc(iat) = qlocacc(iat) + alpha * W_ii * cache%clist(ist) * detai &
         & + beta * self%kqchi(izp) * v(iat)

      ! Second derivatives w.r.t. the CN and the local charge from the hardness
      kcnacc(ist) = kcnacc(ist) + alpha * W_ii * cache%clist(ist) * d2hdcni2
      kqlocacc(iat) = alpha * W_ii * cache%clist(ist) * d2etai

      ! The diagonal weights multiply dCN_i/dR_i in both roles
      dcnacc(:, 2, ist) = dcnacc(:, 1, ist)
      dqlocacc(:, 2, ist) = dqlocacc(:, 1, ist)
   end do
   !$omp end do
   deallocate(hpos)
   !$omp end parallel

   ! Entries with CN or local charge derivatives, the products with the
   ! derivatives of the second neighbors only run over these entries
   call get_cn_entries(mol, list, cache, cninl, cnidx)

   ! Products of the CN and local charge derivatives with their vector weights,
   ! as in the dense variant the block of the atoms xat and jat collects
   ! dCN_iat/dR_xat gcn_iat(jat)^T + gcn_iat(xat) dCN_iat/dR_jat^T over the
   ! atoms iat, and the same for the local charges. Each iteration only updates
   ! the row of atom jat of the Hessian, hpos holds the position of each atom in
   ! this row.
   !$omp parallel default(none) &
   !$omp shared(mol, list, cache, dcnacc, dqlocacc, kcnacc, kqlocacc, hessian) &
   !$omp shared(cninl, cnidx) &
   !$omp private(iat, jat, ic, it, ntouch, kat, lat, hat, cat) &
   !$omp private(dcnj, dqlocj, dcnx, gcni, gqloci, gcn, gqloc, mark, touch, hpos)
   allocate(gcn(3, mol%nat), gqloc(3, mol%nat), source=0.0_wp)
   allocate(mark(mol%nat), source=.false.)
   allocate(touch(mol%nat), hpos(mol%nat))
   hpos(:) = 0_i8
   !$omp do schedule(runtime)
   do jat = 1, mol%nat
      do hat = hessian%list%inl(jat), hessian%list%inl(jat + 1) - 1
         hpos(hessian%list%nlat(hat)) = hat
      end do

      ! Weights gcn(:, iat) of dCN_iat/dR_jat and gqloc(:, iat) of
      ! dqloc_iat/dR_jat: the second derivatives w.r.t. the CNs give
      ! K(iat, k) dCN_k/dR_jat, which reaches the neighbors iat of the atoms k
      ! with CN derivatives w.r.t. R_jat, and the mixed weights of the pairs of
      ! jat are stored in the entries of the row of jat. The second derivatives
      ! w.r.t. the local charges are diagonal. The diagonal entry of K touches
      ! every neighbor of jat.
      ntouch = 0
      do cat = cninl(jat), cninl(jat + 1) - 1
         kat = cnidx(cat)
         dcnj = cache%dcndrlist(:, kat)
         do lat = list%inl(list%nlat(kat)), list%inl(list%nlat(kat) + 1) - 1
            iat = list%nlat(lat)
            if (.not. mark(iat)) then
               mark(iat) = .true.
               ntouch = ntouch + 1
               touch(ntouch) = iat
            end if
            gcn(:, iat) = gcn(:, iat) + kcnacc(lat) * dcnj
         end do
      end do
      do kat = list%inl(jat), list%inl(jat + 1) - 1
         iat = list%nlat(kat)
         gcn(:, iat) = gcn(:, iat) + dcnacc(:, 2, kat)
         gqloc(:, iat) = dqlocacc(:, 2, kat) &
            & + kqlocacc(iat) * cache%dqlocdrlist(:, kat)
      end do

      ! First product, dCN_iat/dR_xat gcn(:, iat)^T and dqloc_iat/dR_xat
      ! gqloc(:, iat)^T, over the entries of xat in the row of iat with CN or
      ! local charge derivatives. Apart from the diagonal, dCN_iat/dR_xat is
      ! the negative of the entry of xat in the row of iat.
      do it = 1, ntouch
         iat = touch(it)
         gcni = gcn(:, iat)
         gqloci = gqloc(:, iat)
         do cat = cninl(iat), cninl(iat + 1) - 1
            lat = cnidx(cat)
            hat = hpos(list%nlat(lat))
            if (hat == 0) cycle
            dcnx = cache%dcndrlist(:, lat)
            if (cat > cninl(iat)) dcnx = -dcnx
            do ic = 1, 3
               hessian%hesslist(:, ic, hat) = hessian%hesslist(:, ic, hat) &
                  & + dcnx * gcni(ic) + cache%dqlocdrlist(:, lat) * gqloci(ic)
            end do
         end do
         gcn(:, iat) = 0.0_wp
         gqloc(:, iat) = 0.0_wp
         mark(iat) = .false.
      end do

      ! Second product, gcn_iat(xat) dCN_iat/dR_jat^T and the same for the local
      ! charges, over the atoms iat with CN or local charge derivatives w.r.t.
      ! R_jat. The weights of dCN_iat/dR_xat of the pairs of iat are stored in
      ! the entries of the row of iat.
      do cat = cninl(jat), cninl(jat + 1) - 1
         kat = cnidx(cat)
         iat = list%nlat(kat)
         dcnj = cache%dcndrlist(:, kat)
         dqlocj = cache%dqlocdrlist(:, kat)
         do lat = list%inl(iat), list%inl(iat + 1) - 1
            hat = hpos(list%nlat(lat))
            if (hat == 0) cycle
            do ic = 1, 3
               hessian%hesslist(:, ic, hat) = hessian%hesslist(:, ic, hat) &
                  & + dcnacc(:, 1, lat) * dcnj(ic) + dqlocacc(:, 1, lat) * dqlocj(ic)
            end do
         end do
      end do

      do hat = hessian%list%inl(jat), hessian%list%inl(jat + 1) - 1
         hpos(hessian%list%nlat(hat)) = 0
      end do
   end do
   !$omp end do
   deallocate(gcn, gqloc, mark, touch, hpos)
   !$omp end parallel

   call add_ncoord_hessian_list(self%ncoord, mol, list, cnacc, hessian)
   call add_ncoord_hessian_list(self%ncoord_en, mol, list, qlocacc, hessian)

   deallocate(qlocacc, cnacc, kcnacc, kqlocacc, dcnacc, dqlocacc, v, cninl, cnidx)

end subroutine get_explicit_hessian_0d_list

!> Add dE/dCN contracted with the Cartesian Hessian of the coordination numbers
!> to the compressed Hessian, as add_coordination_number_hessian_list of
!> mctc-lib for a molecular system with a complete neighborlist and without a
!> maximum CN cutoff. Each iteration only updates the row of atom iat of the
!> Hessian, blocks outside of the pattern of its list are dropped.
subroutine add_ncoord_hessian_list(ncoord, mol, list, dEdcn, hessian)
   !> Coordination number container
   class(ncoord_type), intent(in) :: ncoord

   !> Molecular structure data
   type(structure_type), intent(in) :: mol

   !> Multicharge neighborlist type (complete)
   type(csr_list), intent(in) :: list

   !> Derivative of the expression w.r.t. the coordination numbers
   real(wp), intent(in) :: dEdcn(:)

   !> Molecular nuclear Hessian in compressed storage with allocated blocks
   type(mchrg_hessian), intent(inout) :: hessian

   integer :: iat, jat, izp, jzp, ic
   integer(i8) :: kat, hat
   real(wp) :: rij(3), r1, r2, cutoff2, den, dEdcnij, countd, countd2
   real(wp) :: box(3, 3)
   integer(i8), allocatable :: hpos(:)

   cutoff2 = ncoord%cutoff**2

   !$omp parallel default(none) &
   !$omp shared(ncoord, mol, list, dEdcn, hessian, cutoff2) &
   !$omp private(iat, jat, izp, jzp, ic, kat, hat, rij, r1, r2, den, dEdcnij) &
   !$omp private(countd, countd2, box, hpos)
   allocate(hpos(mol%nat), source=0_i8)
   !$omp do schedule(runtime)
   do iat = 1, mol%nat
      izp = mol%id(iat)
      do hat = hessian%list%inl(iat), hessian%list%inl(iat + 1) - 1
         hpos(hessian%list%nlat(hat)) = hat
      end do
      do kat = list%inl(iat) + 1, list%inl(iat + 1) - 1
         jat = list%nlat(kat)
         jzp = mol%id(jat)
         rij = mol%xyz(:, iat) - mol%xyz(:, jat)
         r2 = sum(rij**2)
         if (r2 > cutoff2 .or. r2 < 1.0e-12_wp) cycle
         r1 = sqrt(r2)
         den = ncoord%get_en_factor(izp, jzp)
         dEdcnij = dEdcn(iat) + dEdcn(jat) * ncoord%directed_factor
         countd = den * ncoord%ncoord_dcount(izp, jzp, r1)
         countd2 = den * ncoord%ncoord_d2count(izp, jzp, r1)

         ! Second derivative of the pair count w.r.t. R_i, R_j gives its negative
         do ic = 1, 3
            box(:, ic) = dEdcnij * (countd2 / r2 - countd / (r2 * r1)) &
               & * rij * rij(ic)
            box(ic, ic) = box(ic, ic) + dEdcnij * countd / r1
         end do
         hat = hessian%list%inl(iat)
         hessian%hesslist(:, :, hat) = hessian%hesslist(:, :, hat) + box
         hat = hpos(jat)
         if (hat > 0) hessian%hesslist(:, :, hat) = hessian%hesslist(:, :, hat) - box
      end do
      do hat = hessian%list%inl(iat), hessian%list%inl(iat + 1) - 1
         hpos(hessian%list%nlat(hat)) = 0
      end do
   end do
   !$omp end do
   deallocate(hpos)
   !$omp end parallel

end subroutine add_ncoord_hessian_list

!> Entries of a complete neighborlist with CN or local charge derivatives above
!> cnthr, given in the rows cninl as positions in the list, the diagonal entry
!> first. The derivatives of the counting functions decay as Gaussians and only
!> reach the close neighbors in the list.
subroutine get_cn_entries(mol, list, cache, cninl, cnidx)
   !> Molecular structure data
   type(structure_type), intent(in) :: mol

   !> Multicharge neighborlist type (complete)
   type(csr_list), intent(in) :: list

   !> Cache with the CN and local charge derivatives on the list
   type(mchrg_cache), intent(in) :: cache

   !> Offsets of the rows of the entries
   integer(i8), allocatable, intent(out) :: cninl(:)

   !> Positions of the entries in the list
   integer(i8), allocatable, intent(out) :: cnidx(:)

   integer :: iat
   integer(i8) :: kat, cnt

   allocate(cninl(mol%nat + 1))
   cninl(1) = 1
   do iat = 1, mol%nat
      cnt = 1
      do kat = list%inl(iat) + 1, list%inl(iat + 1) - 1
         if (any(abs(cache%dcndrlist(:, kat)) > cnthr) &
            & .or. any(abs(cache%dqlocdrlist(:, kat)) > cnthr)) cnt = cnt + 1
      end do
      cninl(iat + 1) = cninl(iat) + cnt
   end do

   allocate(cnidx(cninl(mol%nat + 1) - 1))
   !$omp parallel do default(none) schedule(runtime) &
   !$omp shared(mol, list, cache, cninl, cnidx) private(iat, kat, cnt)
   do iat = 1, mol%nat
      cnt = cninl(iat)
      cnidx(cnt) = list%inl(iat)
      do kat = list%inl(iat) + 1, list%inl(iat + 1) - 1
         if (any(abs(cache%dcndrlist(:, kat)) > cnthr) &
            & .or. any(abs(cache%dqlocdrlist(:, kat)) > cnthr)) then
            cnt = cnt + 1
            cnidx(cnt) = kat
         end if
      end do
   end do
   !$omp end parallel do

end subroutine get_cn_entries

end module multicharge_model_eeqbc
