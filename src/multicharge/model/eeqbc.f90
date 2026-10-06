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
   use multicharge_model_type, only : get_dir_trans, mchrg_model_type
   use multicharge_blas, only : gemm, gemv, symv
   use multicharge_model_cache, only : mchrg_cache
   use multicharge_ncoord, only : add_dcndr, add_dqlocdr, get_dcndr_pair, &
      & get_dqlocdr_pair
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
      !> Calculate gradient
      procedure :: get_grad
   end type eeqbc_model

   real(wp), parameter :: sqrtpi = sqrt(pi)
   real(wp), parameter :: sqrt2pi = sqrt(2.0_wp / pi)
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

!> Update coordination numbers and local charges, and set up the Wigner-Seitz
!> cell if periodic
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
   !> Multicharge neighborlist type
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

   integer :: iat, izp, img
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

!> Derivative of the electronegativity of each atom through its CN and local
!> charge, dtmp(:, :, a) = kcnchi(a) * dCN(a)/dr + kqchi(a) * dqloc(a)/dr
subroutine get_dtmp_derivs(self, mol, cache, dtmpdr, dtmpdL)
   !> EEQBC model type
   class(eeqbc_model), intent(in) :: self
   !> Structure type
   type(structure_type), intent(in) :: mol
   !> Multicharge cache
   type(mchrg_cache), intent(in) :: cache
   !> Derivative w.r.t. the Cartesian coordinates
   real(wp), intent(out) :: dtmpdr(:, :, :)
   !> Derivative w.r.t. strain deformations
   real(wp), intent(out) :: dtmpdL(:, :, :)

   dtmpdr(:, :, :) = 0.0_wp
   dtmpdL(:, :, :) = 0.0_wp
   call add_dcndr(self%ncoord, mol, cache%trans, cache%cn, &
      & self%kcnchi(mol%id), dtmpdr, dtmpdL)
   call add_dqlocdr(self%ncoord_en, mol, cache%trans, cache%qloc, &
      & self%kqchi(mol%id), dtmpdr, dtmpdL)

end subroutine get_dtmp_derivs

!> Contract weights with the on-the-fly CN and local charge derivatives,
!> dxdr(:, :, b) += sum_a wcn(a, b) * dCN(a)/dr + wq(b) * dqloc(b)/dr
subroutine add_cn_weighted(self, mol, cache, wcn, wq, dxdr, dxdL)
   !> EEQBC model type
   class(eeqbc_model), intent(in) :: self
   !> Structure type
   type(structure_type), intent(in) :: mol
   !> Multicharge cache
   type(mchrg_cache), intent(in) :: cache
   !> Weights of the CN derivatives
   real(wp), intent(in) :: wcn(:, :)
   !> Weights of the local charge derivatives
   real(wp), intent(in) :: wq(:)
   !> Derivative w.r.t. the Cartesian coordinates
   real(wp), intent(inout), contiguous :: dxdr(:, :, :)
   !> Derivative w.r.t. strain deformations
   real(wp), intent(inout), contiguous :: dxdL(:, :, :)

   real(wp), allocatable :: dcndr(:, :, :), dcndL(:, :, :), unity(:)

   allocate(dcndr(3, mol%nat, mol%nat), dcndL(3, 3, mol%nat), source=0.0_wp)
   allocate(unity(mol%nat), source=1.0_wp)
   call add_dcndr(self%ncoord, mol, cache%trans, cache%cn, unity, dcndr, dcndL)
   call gemm(dcndr, wcn, dxdr, beta=1.0_wp)
   call gemm(dcndL, wcn, dxdL, beta=1.0_wp)
   call add_dqlocdr(self%ncoord_en, mol, cache%trans, cache%qloc, wq, dxdr, dxdL)

end subroutine add_cn_weighted

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
   !> Multicharge neighborlist type (upper triangle)
   type(csr_list), intent(in), optional :: list

   if (.not. allocated(cache%dabdL)) allocate(cache%dabdL(3, 3, ndim))

   if (present(list)) then
      if (.not. allocated(cache%dabdrij)) allocate(cache%dabdrij(3, size(list%nlat)))
      if (.not. allocated(cache%dabdrji)) allocate(cache%dabdrji(3, size(list%nlat)))
      if (.not. allocated(cache%dabdrdiag)) allocate(cache%dabdrdiag(3, mol%nat))
      if (any(mol%periodic)) then
         call get_partial_derivs_3d_list(self, mol, list, cache, alpha, beta)
      else
         call get_partial_derivs_0d_list(self, mol, list, cache, alpha, beta)
      end if
   else
      if (.not. allocated(cache%dabdr)) allocate(cache%dabdr(3, mol%nat, ndim))
      if (any(mol%periodic)) then
         call get_partial_derivs_3d(self, mol, ndim, cache, alpha, beta)
      else
         call get_partial_derivs_0d(self, mol, ndim, cache, alpha, beta)
      end if
   end if

end subroutine get_partial_derivs

!> Compute the linear combination of the Coulomb matrix derivatives (multiplied
!> by the charge vector) and the electronegativity vector derivatives,
!> alpha * dA/dR*q + beta * dX/dR, for a non-periodic system.
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
   real(wp), allocatable :: wcn(:, :), wq(:)

   allocate(dtmpdr(3, mol%nat, ndim), dtmpdL(3, 3, ndim))
   ! Weights of the CN derivatives, column b collects d(Aq)_b/dCN(a),
   ! and of the local charge derivatives (diagonal)
   allocate(wcn(mol%nat, ndim), wq(mol%nat), source=0.0_wp)

   ! CN and effective charge derivative
   call get_dtmp_derivs(self, mol, cache, dtmpdr, dtmpdL)

   ! EN derivative through the capacitance matrix, initializes dabdr and dabdL
   call gemm(dtmpdr, cache%cmat, cache%dabdr, alpha=beta)
   call gemm(dtmpdL, cache%cmat, cache%dabdL, alpha=beta)

   ! Each iteration only updates the derivatives of row iat (last index)
   !$omp parallel do default(none) schedule(runtime) &
   !$omp shared(cache, mol, self, alpha, beta, wcn, wq) &
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
         wcn(jat, iat) = wcn(jat, iat) - wtmp * radj * dradj * gam**3

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
      wq(iat) = alpha * self%kqeta_pre * self%kqeta(izp) &
      & / cosh(self%kqeta(izp) * cache%qloc(iat))**2 * cache%vrhs(iat) &
      & * cache%cmat(iat, iat)

      ! Effective charge width derivative
      dtmp = -alpha * sqrt2pi * dradi / (radi**2) * cache%vrhs(iat) * cache%cmat(iat, iat)
      wcn(iat, iat) = wcn(iat, iat) + dtmp + wcni

      ! Capacitance derivative
      dtmp = alpha * hardi * cache%vrhs(iat)
      cache%dabdr(:, iat, iat) = +dtmp * cache%dcdr(:, iat, iat) &
      & + cache%dabdr(:, iat, iat)
      cache%dabdL(:, :, iat) = +dtmp * cache%dcdL(:, :, iat) + cache%dabdL(:, :, iat)
   end do
   !$omp end parallel do

   call add_cn_weighted(self, mol, cache, wcn, wq, cache%dabdr, cache%dabdL)

end subroutine get_partial_derivs_0d

!> Compute the linear combination of the Coulomb matrix derivatives (multiplied
!> by the charge vector) and the electronegativity vector derivatives,
!> alpha * dA/dR*q + beta * dX/dR, for a periodic system.
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
   real(wp) :: vec(3), gam, dtmp, norm_cn, rvdw, wsw, dgam, ctmp, hardi, hardj
   real(wp) :: radi, radj, dradi, dradj, dG(3), dS(3, 3)
   real(wp) :: dgami, dgamj, capi, capj, wqi
   real(wp), allocatable :: dtrans(:, :), dtmpdr(:, :, :), dtmpdL(:, :, :)
   real(wp), allocatable :: wcn(:, :), wq(:)

   ! Thread-private arrays for reduction
   real(wp), allocatable :: dabdr_local(:, :, :), dabdL_local(:, :, :)
   real(wp), allocatable :: wcn_local(:, :)

   call get_dir_trans(mol, dtrans, cutoff)

   allocate(dtmpdr(3, mol%nat, ndim), dtmpdL(3, 3, ndim))
   ! Weights of the CN derivatives, column b collects d(Aq)_b/dCN(a),
   ! and of the local charge derivatives (diagonal)
   allocate(wcn(mol%nat, ndim), wq(mol%nat), source=0.0_wp)

   ! CN and effective charge derivative
   call get_dtmp_derivs(self, mol, cache, dtmpdr, dtmpdL)

   ! EN derivative through the capacitance matrix, initializes dabdr and dabdL
   call gemm(dtmpdr, cache%cmat, cache%dabdr, alpha=beta)
   call gemm(dtmpdL, cache%cmat, cache%dabdL, alpha=beta)

   !$omp parallel default(none) &
   !$omp shared(self, mol, cache, dtrans, alpha, beta, wcn, wq) &
   !$omp private(iat, izp, jat, jzp, img, gam, vec, dtmp, norm_cn, rvdw, wsw) &
   !$omp private(radi, radj, dradi, dradj, capi, capj, dgami, dgamj, dgam) &
   !$omp private(ctmp, hardi, hardj, wqi, dG, dS) &
   !$omp private(dabdr_local, dabdL_local, wcn_local)
   allocate(dabdr_local(3, mol%nat, size(cache%dabdr, 3)), source=0.0_wp)
   allocate(dabdL_local(3, 3, size(cache%dabdL, 3)), source=0.0_wp)
   allocate(wcn_local(mol%nat, size(wcn, 2)), source=0.0_wp)
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
         ! Derivative of gamma w.r.t. CN(i) and CN(j)
         dgami = -radi * dradi * gam**3
         dgamj = -radj * dradj * gam**3

         wsw = 1.0_wp / real(cache%wsc%nimg(jat, iat), wp)
         do img = 1, cache%wsc%nimg(jat, iat)
            vec = mol%xyz(:, jat) - mol%xyz(:, iat) + cache%wsc%trans(:, &
            & cache%wsc%tridx(img, jat, iat))

            call get_damat_dir(vec, dtrans, capi, capj, rvdw, self%kbc, gam, dG, dS, dgam)
            dG = alpha * dG * wsw
            dS = alpha * dS * wsw
            dgam = alpha * dgam * wsw

            ! Explicit derivative
            dabdr_local(:, iat, iat) = dabdr_local(:, iat, iat) - dG * cache%vrhs(jat)
            dabdr_local(:, jat, jat) = dabdr_local(:, jat, jat) + dG * cache%vrhs(iat)
            dabdr_local(:, iat, jat) = dabdr_local(:, iat, jat) - dG * cache%vrhs(iat)
            dabdr_local(:, jat, iat) = dabdr_local(:, jat, iat) + dG * cache%vrhs(jat)
            dabdL_local(:, :, jat) = dabdL_local(:, :, jat) + dS * cache%vrhs(iat)
            dabdL_local(:, :, iat) = dabdL_local(:, :, iat) + dS * cache%vrhs(jat)

            ! Effective charge width derivative
            wcn_local(iat, iat) = wcn_local(iat, iat) - dgam * cache%vrhs(jat) * dgami
            wcn_local(jat, iat) = wcn_local(jat, iat) - dgam * cache%vrhs(jat) * dgamj
            wcn_local(iat, jat) = wcn_local(iat, jat) - dgam * cache%vrhs(iat) * dgami
            wcn_local(jat, jat) = wcn_local(jat, jat) - dgam * cache%vrhs(iat) * dgamj

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
      end do

      ! EN derivative: capacitance matrix derivative diagonal
      dabdL_local(:, :, iat) = dabdL_local(:, :, iat) &
         & + beta * cache%xtmp(iat) * cache%dcdL(:, :, iat)

      ! Self-images, i = j and T != 0
      gam = 1.0_wp / sqrt(2.0_wp * radi**2)
      dgami = -(2.0_wp * radi * dradi) * gam**3
      rvdw = self%rvdw(izp, izp)
      wqi = 0.0_wp
      wsw = 1.0_wp / real(cache%wsc%nimg(iat, iat), wp)
      do img = 1, cache%wsc%nimg(iat, iat)
         vec = cache%wsc%trans(:, cache%wsc%tridx(img, iat, iat))

         ! EN derivative through the self-image capacitance
         call get_cpair_dir(self%kbc, vec, dtrans, rvdw, capi, capi, ctmp)
         ctmp = beta * ctmp * wsw
         wcn_local(iat, iat) = wcn_local(iat, iat) - ctmp * self%kcnchi(izp)
         wqi = wqi - ctmp * self%kqchi(izp)

         call get_damat_dir(vec, dtrans, capi, capi, rvdw, self%kbc, gam, dG, dS, dgam)
         dS = alpha * dS * wsw
         dgam = alpha * dgam * wsw

         ! Explicit derivative
         dabdL_local(:, :, iat) = dabdL_local(:, :, iat) + dS * cache%vrhs(iat)

         ! Effective charge width derivative
         wcn_local(iat, iat) = wcn_local(iat, iat) - cache%vrhs(iat) * dgam * dgami

         ! Capacitance derivative
         call get_damat_dc_dir(vec, dtrans, capi, capi, rvdw, self%kbc, gam, dG, dS)
         dabdL_local(:, :, iat) = dabdL_local(:, :, iat) &
            & - alpha * cache%vrhs(iat) * dS * wsw
      end do

      ! Hardness derivative, wq(iat) is exclusive to this iteration
      wq(iat) = wqi + alpha * self%kqeta_pre * self%kqeta(izp) &
      & / cosh(self%kqeta(izp) * cache%qloc(iat))**2 * cache%vrhs(iat) &
      & * cache%cmat(iat, iat)

      ! Effective charge width derivative
      dtmp = -alpha * sqrt2pi * dradi / (radi**2) * cache%vrhs(iat) * cache%cmat(iat, iat)
      wcn_local(iat, iat) = wcn_local(iat, iat) + dtmp

      ! Capacitance derivative
      dtmp = alpha * hardi * cache%vrhs(iat)
      dabdr_local(:, iat, iat) = dabdr_local(:, iat, iat) + dtmp * cache%dcdr(:, iat, iat)
      dabdL_local(:, :, iat) = dabdL_local(:, :, iat) + dtmp * cache%dcdL(:, :, iat)
   end do
   !$omp end do
   !$omp critical (get_partial_derivs_3d_)
   cache%dabdr(:, :, :) = cache%dabdr + dabdr_local
   cache%dabdL(:, :, :) = cache%dabdL + dabdL_local
   wcn(:, :) = wcn + wcn_local
   !$omp end critical (get_partial_derivs_3d_)
   deallocate(dabdr_local, dabdL_local, wcn_local)
   !$omp end parallel

   call add_cn_weighted(self, mol, cache, wcn, wq, cache%dabdr, cache%dabdL)

end subroutine get_partial_derivs_3d

!> Compute alpha * dA/dR*q + beta * dX/dR for a non-periodic system using an
!> upper-triangle neighborlist.
!>
!> For the list entry kat of row i with j = list%nlat(kat) the derivatives are
!> stored as dabdrij(:, kat) = dB(j)/dR(i), dabdrji(:, kat) = dB(i)/dR(j) and
!> dabdrdiag(:, i) = dB(i)/dR(i). Contributions outside the sparsity pattern of
!> the list are neglected.
subroutine get_partial_derivs_0d_list(self, mol, list, cache, alpha, beta)
   !> EEQBC model type
   class(eeqbc_model), intent(in) :: self
   !> Molecular structure data
   type(structure_type), intent(in) :: mol
   !> Multicharge neighborlist type (upper triangle)
   type(csr_list), intent(in) :: list
   !> Multicharge cache
   type(mchrg_cache), intent(inout) :: cache
   !> dA/dR*q multiplier
   real(wp), intent(in) :: alpha
   !> dX/dR multiplier
   real(wp), intent(in) :: beta

   integer :: iat, jat, izp, jzp
   integer(i8) :: kat
   real(wp) :: vec(3), r1, r2, gam, arg, dtmp, norm_cn, hardi, hardj, cij
   real(wp) :: radi, radj, dradi, dradj, dG(3), dS(3, 3), ti, tj
   real(wp), allocatable :: wdiag(:), wij(:), wji(:), qdiag(:), qij(:), qji(:)

   ! Thread-private arrays for reduction
   real(wp), allocatable :: bdiag_local(:, :), bL_local(:, :, :), wdiag_local(:)

   ! Weights of the CN and local charge derivatives in list storage
   allocate(wij(size(list%nlat)), wji(size(list%nlat)), source=0.0_wp)
   allocate(qij(size(list%nlat)), qji(size(list%nlat)), source=0.0_wp)
   allocate(wdiag(mol%nat), qdiag(mol%nat), source=0.0_wp)

   cache%dabdrij(:, :) = 0.0_wp
   cache%dabdrji(:, :) = 0.0_wp
   cache%dabdrdiag(:, :) = 0.0_wp
   cache%dabdL(:, :, :) = 0.0_wp

   !$omp parallel default(none) &
   !$omp shared(self, mol, list, cache, alpha, beta) &
   !$omp shared(wdiag, wij, wji, qdiag, qij, qji) &
   !$omp private(iat, jat, izp, jzp, kat, vec, r1, r2, gam, arg, dtmp, norm_cn) &
   !$omp private(hardi, hardj, cij, radi, radj, dradi, dradj, dG, dS, ti, tj) &
   !$omp private(bdiag_local, bL_local, wdiag_local)
   allocate(bdiag_local(3, mol%nat), bL_local(3, 3, mol%nat), source=0.0_wp)
   allocate(wdiag_local(mol%nat), source=0.0_wp)
   !$omp do schedule(runtime)
   do iat = 1, mol%nat
      izp = mol%id(iat)
      ! Effective charge width of i
      norm_cn = 1.0_wp / self%avg_cn(izp)
      radi = self%rad(izp) * exp(-self%kcnrad(izp) * cache%cn(iat) * norm_cn)
      dradi = -self%kcnrad(izp) * norm_cn * radi
      ! Effective hardness of i
      hardi = self%eta(izp) + self%kqeta_pre &
      & * tanh(self%kqeta(izp) * cache%qloc(iat)) + sqrt2pi / radi

      ! Entries kat of row iat are exclusive to this iteration
      do kat = list%inl(iat) + 1, list%inl(iat + 1) - 1
         jat = list%nlat(kat)
         jzp = mol%id(jat)
         vec = mol%xyz(:, jat) - mol%xyz(:, iat)
         r2 = dot_product(vec, vec)
         r1 = sqrt(r2)
         cij = cache%clist(kat)
         ! Effective charge width of j
         norm_cn = 1.0_wp / self%avg_cn(jzp)
         radj = self%rad(jzp) * exp(-self%kcnrad(jzp) * cache%cn(jat) * norm_cn)
         dradj = -self%kcnrad(jzp) * norm_cn * radj
         ! Effective hardness of j
         hardj = self%eta(jzp) + self%kqeta_pre &
         & * tanh(self%kqeta(jzp) * cache%qloc(jat)) + sqrt2pi / radj

         ! Coulomb interaction of Gaussian charges
         gam = 1.0_wp / sqrt(radi**2 + radj**2)
         arg = gam * gam * r2

         ! EN derivative through the capacitance matrix
         wij(kat) = beta * self%kcnchi(izp) * cij
         wji(kat) = beta * self%kcnchi(jzp) * cij
         qij(kat) = beta * self%kqchi(izp) * cij
         qji(kat) = beta * self%kqchi(jzp) * cij

         ! Effective charge width derivative
         dtmp = alpha * 2.0_wp * exp(-arg) / sqrtpi * cij
         ti = -dtmp * radi * dradi * gam**3
         tj = -dtmp * radj * dradj * gam**3
         wdiag_local(iat) = wdiag_local(iat) + ti * cache%vrhs(jat)
         wij(kat) = wij(kat) + ti * cache%vrhs(iat)
         wji(kat) = wji(kat) + tj * cache%vrhs(jat)
         wdiag_local(jat) = wdiag_local(jat) + tj * cache%vrhs(iat)

         ! Explicit derivative
         dtmp = 2.0_wp * gam * exp(-arg) / (sqrtpi * r2) - erf(sqrt(arg)) / (r2 * r1)
         dG(:) = alpha * dtmp * cij * vec
         dS(:, :) = spread(dG, 1, 3) * spread(vec, 2, 3)
         bdiag_local(:, iat) = bdiag_local(:, iat) - cache%vrhs(jat) * dG
         cache%dabdrij(:, kat) = cache%dabdrij(:, kat) - cache%vrhs(iat) * dG
         cache%dabdrji(:, kat) = cache%dabdrji(:, kat) + cache%vrhs(jat) * dG
         bdiag_local(:, jat) = bdiag_local(:, jat) + cache%vrhs(iat) * dG
         bL_local(:, :, iat) = bL_local(:, :, iat) + cache%vrhs(jat) * dS
         bL_local(:, :, jat) = bL_local(:, :, jat) + cache%vrhs(iat) * dS

         ! Capacitance derivatives
         call get_dcpair(self%kbc, vec, self%rvdw(izp, jzp), self%cap(izp), &
         & self%cap(jzp), dG, dS)
         dtmp = alpha * erf(r1 * gam) / r1
         bdiag_local(:, iat) = bdiag_local(:, iat) &
            & + (dtmp * cache%vrhs(jat) + beta * cache%xtmp(jat)) * dG
         cache%dabdrij(:, kat) = cache%dabdrij(:, kat) &
            & + (dtmp * cache%vrhs(iat) - alpha * hardj * cache%vrhs(jat) &
            & - beta * (cache%xtmp(jat) - cache%xtmp(iat))) * dG
         cache%dabdrji(:, kat) = cache%dabdrji(:, kat) &
            & + (-dtmp * cache%vrhs(jat) + alpha * hardi * cache%vrhs(iat) &
            & + beta * (cache%xtmp(iat) - cache%xtmp(jat))) * dG
         bdiag_local(:, jat) = bdiag_local(:, jat) &
            & - (dtmp * cache%vrhs(iat) + beta * cache%xtmp(iat)) * dG
         bL_local(:, :, iat) = bL_local(:, :, iat) &
            & - (dtmp * cache%vrhs(jat) + beta * cache%xtmp(jat)) * dS
         bL_local(:, :, jat) = bL_local(:, :, jat) &
            & - (dtmp * cache%vrhs(iat) + beta * cache%xtmp(iat)) * dS
      end do

      ! EN derivative through the capacitance matrix and hardness derivative,
      ! qdiag(iat) is exclusive to this iteration
      cij = cache%clist(list%inl(iat))
      wdiag_local(iat) = wdiag_local(iat) + beta * self%kcnchi(izp) * cij
      qdiag(iat) = beta * self%kqchi(izp) * cij + alpha * self%kqeta_pre &
      & * self%kqeta(izp) / cosh(self%kqeta(izp) * cache%qloc(iat))**2 &
      & * cache%vrhs(iat) * cij

      ! Effective charge width derivative
      wdiag_local(iat) = wdiag_local(iat) &
         & - alpha * sqrt2pi * dradi / (radi**2) * cache%vrhs(iat) * cij

      ! Intrinsic capacitance derivative
      dtmp = alpha * hardi * cache%vrhs(iat) + beta * cache%xtmp(iat)
      bdiag_local(:, iat) = bdiag_local(:, iat) + dtmp * cache%dcdrdiag(:, iat)
      bL_local(:, :, iat) = bL_local(:, :, iat) + dtmp * cache%dcdL(:, :, iat)
   end do
   !$omp end do
   !$omp critical (get_partial_derivs_0d_list_)
   cache%dabdrdiag(:, :) = cache%dabdrdiag + bdiag_local
   cache%dabdL(:, :, :mol%nat) = cache%dabdL(:, :, :mol%nat) + bL_local
   wdiag(:) = wdiag + wdiag_local
   !$omp end critical (get_partial_derivs_0d_list_)
   deallocate(bdiag_local, bL_local, wdiag_local)
   !$omp end parallel

   call add_cn_weighted_list(self, mol, list, cache, wdiag, wij, wji, &
      & qdiag, qij, qji, cache%dabdrdiag, cache%dabdrij, cache%dabdrji, cache%dabdL)

end subroutine get_partial_derivs_0d_list

!> Compute alpha * dA/dR*q + beta * dX/dR for a periodic system using an
!> upper-triangle neighborlist with Wigner-Seitz images.
!>
!> The storage follows get_partial_derivs_0d_list, contributions outside the
!> sparsity pattern of the list are neglected.
subroutine get_partial_derivs_3d_list(self, mol, list, cache, alpha, beta)
   !> EEQBC model type
   class(eeqbc_model), intent(in) :: self
   !> Molecular structure data
   type(structure_type), intent(in) :: mol
   !> Multicharge neighborlist type (upper triangle)
   type(csr_list), intent(in) :: list
   !> Multicharge cache
   type(mchrg_cache), intent(inout) :: cache
   !> dA/dR*q multiplier
   real(wp), intent(in) :: alpha
   !> dX/dR multiplier
   real(wp), intent(in) :: beta

   integer :: iat, jat, izp, jzp
   integer(i8) :: kat, img
   real(wp) :: vec(3), gam, dtmp, norm_cn, rvdw, wsw, dgam, hardi, hardj, cij
   real(wp) :: radi, radj, dradi, dradj, dG(3), dS(3, 3), ti, tj, capi, capj
   real(wp) :: ctmp, cself
   real(wp), allocatable :: dtrans(:, :)
   real(wp), allocatable :: wdiag(:), wij(:), wji(:), qdiag(:), qij(:), qji(:)

   ! Thread-private arrays for reduction
   real(wp), allocatable :: bdiag_local(:, :), bL_local(:, :, :), wdiag_local(:)

   call get_dir_trans(mol, dtrans, cutoff)

   ! Weights of the CN and local charge derivatives in list storage
   allocate(wij(size(list%nlat)), wji(size(list%nlat)), source=0.0_wp)
   allocate(qij(size(list%nlat)), qji(size(list%nlat)), source=0.0_wp)
   allocate(wdiag(mol%nat), qdiag(mol%nat), source=0.0_wp)

   cache%dabdrij(:, :) = 0.0_wp
   cache%dabdrji(:, :) = 0.0_wp
   cache%dabdrdiag(:, :) = 0.0_wp
   cache%dabdL(:, :, :) = 0.0_wp

   !$omp parallel default(none) &
   !$omp shared(self, mol, list, cache, dtrans, alpha, beta) &
   !$omp shared(wdiag, wij, wji, qdiag, qij, qji) &
   !$omp private(iat, jat, izp, jzp, kat, img, vec, gam, dtmp, norm_cn, rvdw, wsw) &
   !$omp private(dgam, hardi, hardj, cij, radi, radj, dradi, dradj, dG, dS, ti, tj) &
   !$omp private(capi, capj, ctmp, cself) &
   !$omp private(bdiag_local, bL_local, wdiag_local)
   allocate(bdiag_local(3, mol%nat), bL_local(3, 3, mol%nat), source=0.0_wp)
   allocate(wdiag_local(mol%nat), source=0.0_wp)
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

      ! Entries kat of row iat are exclusive to this iteration
      do kat = list%inl(iat) + 1, list%inl(iat + 1) - 1
         jat = list%nlat(kat)
         jzp = mol%id(jat)
         capj = self%cap(jzp)
         rvdw = self%rvdw(izp, jzp)
         cij = cache%clist(kat)

         ! EN derivative through the capacitance matrix
         wij(kat) = beta * self%kcnchi(izp) * cij
         wji(kat) = beta * self%kcnchi(jzp) * cij
         qij(kat) = beta * self%kqchi(izp) * cij
         qji(kat) = beta * self%kqchi(jzp) * cij

         if (cache%wsc%nimg_list(kat) == 0) cycle
         wsw = 1.0_wp / real(cache%wsc%nimg_list(kat), wp)

         ! Effective charge width of j
         norm_cn = 1.0_wp / self%avg_cn(jzp)
         radj = self%rad(jzp) * exp(-self%kcnrad(jzp) * cache%cn(jat) * norm_cn)
         dradj = -self%kcnrad(jzp) * norm_cn * radj
         ! Effective hardness of j
         hardj = self%eta(jzp) + self%kqeta_pre &
         & * tanh(self%kqeta(jzp) * cache%qloc(jat)) + sqrt2pi / radj
         gam = 1.0_wp / sqrt(radi**2 + radj**2)

         do img = cache%wsc%itr_list(kat), cache%wsc%itr_list(kat + 1) - 1
            vec = mol%xyz(:, jat) - mol%xyz(:, iat) + cache%wsc%trans(:, &
            & cache%wsc%tridx_list(img))
            call get_damat_dir(vec, dtrans, capi, capj, rvdw, self%kbc, gam, dG, dS, dgam)

            ! Effective charge width derivative
            ti = alpha * dgam * wsw * radi * dradi * gam**3
            tj = alpha * dgam * wsw * radj * dradj * gam**3
            wdiag_local(iat) = wdiag_local(iat) + ti * cache%vrhs(jat)
            wij(kat) = wij(kat) + ti * cache%vrhs(iat)
            wji(kat) = wji(kat) + tj * cache%vrhs(jat)
            wdiag_local(jat) = wdiag_local(jat) + tj * cache%vrhs(iat)

            ! Explicit derivative
            dG = alpha * dG * wsw
            dS = alpha * dS * wsw
            bdiag_local(:, iat) = bdiag_local(:, iat) - cache%vrhs(jat) * dG
            cache%dabdrij(:, kat) = cache%dabdrij(:, kat) - cache%vrhs(iat) * dG
            cache%dabdrji(:, kat) = cache%dabdrji(:, kat) + cache%vrhs(jat) * dG
            bdiag_local(:, jat) = bdiag_local(:, jat) + cache%vrhs(iat) * dG
            bL_local(:, :, iat) = bL_local(:, :, iat) + cache%vrhs(jat) * dS
            bL_local(:, :, jat) = bL_local(:, :, jat) + cache%vrhs(iat) * dS

            ! Capacitance derivative of the Coulomb interaction
            call get_damat_dc_dir(vec, dtrans, capi, capj, rvdw, self%kbc, gam, dG, dS)
            dG = alpha * dG * wsw
            dS = alpha * dS * wsw
            bdiag_local(:, iat) = bdiag_local(:, iat) + cache%vrhs(jat) * dG
            cache%dabdrij(:, kat) = cache%dabdrij(:, kat) + cache%vrhs(iat) * dG
            cache%dabdrji(:, kat) = cache%dabdrji(:, kat) - cache%vrhs(jat) * dG
            bdiag_local(:, jat) = bdiag_local(:, jat) - cache%vrhs(iat) * dG
            bL_local(:, :, iat) = bL_local(:, :, iat) - cache%vrhs(jat) * dS
            bL_local(:, :, jat) = bL_local(:, :, jat) - cache%vrhs(iat) * dS

            ! Capacitance derivatives of the hardness and the EN
            call get_dcpair_dir(self%kbc, vec, dtrans, rvdw, capi, capj, dG, dS)
            dG = dG * wsw
            dS = dS * wsw
            bdiag_local(:, iat) = bdiag_local(:, iat) + beta * cache%xtmp(jat) * dG
            cache%dabdrij(:, kat) = cache%dabdrij(:, kat) &
               & - (alpha * hardj * cache%vrhs(jat) &
               & + beta * (cache%xtmp(jat) - cache%xtmp(iat))) * dG
            cache%dabdrji(:, kat) = cache%dabdrji(:, kat) &
               & + (alpha * hardi * cache%vrhs(iat) &
               & + beta * (cache%xtmp(iat) - cache%xtmp(jat))) * dG
            bdiag_local(:, jat) = bdiag_local(:, jat) - beta * cache%xtmp(iat) * dG
            bL_local(:, :, iat) = bL_local(:, :, iat) - beta * cache%xtmp(jat) * dS
            bL_local(:, :, jat) = bL_local(:, :, jat) - beta * cache%xtmp(iat) * dS
         end do
      end do

      ! Self-images, i = j and T != 0
      gam = 1.0_wp / sqrt(2.0_wp * radi**2)
      rvdw = self%rvdw(izp, izp)
      cself = 0.0_wp
      if (cache%wsc%nimg_list(list%inl(iat)) > 0) then
         wsw = 1.0_wp / real(cache%wsc%nimg_list(list%inl(iat)), wp)
         do img = cache%wsc%itr_list(list%inl(iat)), &
         & cache%wsc%itr_list(list%inl(iat) + 1) - 1
            vec = cache%wsc%trans(:, cache%wsc%tridx_list(img))

            ! Explicit and effective charge width derivative
            call get_damat_dir(vec, dtrans, capi, capi, rvdw, self%kbc, gam, dG, dS, dgam)
            bL_local(:, :, iat) = bL_local(:, :, iat) &
               & + alpha * cache%vrhs(iat) * dS * wsw
            wdiag_local(iat) = wdiag_local(iat) + alpha * cache%vrhs(iat) &
               & * dgam * wsw * 2.0_wp * radi * dradi * gam**3

            ! Capacitance derivative of the Coulomb interaction
            call get_damat_dc_dir(vec, dtrans, capi, capi, rvdw, self%kbc, gam, dG, dS)
            bL_local(:, :, iat) = bL_local(:, :, iat) &
               & - alpha * cache%vrhs(iat) * dS * wsw

            ! Capacitance derivative of the EN
            call get_dcpair_dir(self%kbc, vec, dtrans, rvdw, capi, capi, dG, dS)
            bL_local(:, :, iat) = bL_local(:, :, iat) &
               & - beta * cache%xtmp(iat) * dS * wsw

            ! Self-image capacitance in the EN
            call get_cpair_dir(self%kbc, vec, dtrans, rvdw, capi, capi, ctmp)
            cself = cself + ctmp * wsw
         end do
      end if

      ! EN derivative through the capacitance matrix and hardness derivative,
      ! qdiag(iat) is exclusive to this iteration
      cij = cache%clist(list%inl(iat))
      wdiag_local(iat) = wdiag_local(iat) + beta * self%kcnchi(izp) * (cij - cself)
      qdiag(iat) = beta * self%kqchi(izp) * (cij - cself) + alpha * self%kqeta_pre &
      & * self%kqeta(izp) / cosh(self%kqeta(izp) * cache%qloc(iat))**2 &
      & * cache%vrhs(iat) * cij

      ! Effective charge width derivative
      wdiag_local(iat) = wdiag_local(iat) &
         & - alpha * sqrt2pi * dradi / (radi**2) * cache%vrhs(iat) * cij

      ! Intrinsic capacitance derivative
      dtmp = alpha * hardi * cache%vrhs(iat) + beta * cache%xtmp(iat)
      bdiag_local(:, iat) = bdiag_local(:, iat) + dtmp * cache%dcdrdiag(:, iat)
      bL_local(:, :, iat) = bL_local(:, :, iat) + dtmp * cache%dcdL(:, :, iat)
   end do
   !$omp end do
   !$omp critical (get_partial_derivs_3d_list_)
   cache%dabdrdiag(:, :) = cache%dabdrdiag + bdiag_local
   cache%dabdL(:, :, :mol%nat) = cache%dabdL(:, :, :mol%nat) + bL_local
   wdiag(:) = wdiag + wdiag_local
   !$omp end critical (get_partial_derivs_3d_list_)
   deallocate(bdiag_local, bL_local, wdiag_local)
   !$omp end parallel

   call add_cn_weighted_list(self, mol, list, cache, wdiag, wij, wji, &
      & qdiag, qij, qji, cache%dabdrdiag, cache%dabdrij, cache%dabdrji, cache%dabdL)

end subroutine get_partial_derivs_3d_list

!> Contract list weights with the on-the-fly CN and local charge derivatives,
!> B(:, x, b) += sum_a W(a, b) * dCN(a)/dR(x) + Q(a, b) * dqloc(a)/dR(x),
!> restricted to the sparsity pattern of the upper-triangle neighborlist.
!>
!> Weights and derivatives are stored as diag(b) for (b, b), ij(kat) for (i, j)
!> and ji(kat) for (j, i), where kat is an entry of row i with j = list%nlat(kat).
!> The list has to contain every pair within the CN cutoff.
subroutine add_cn_weighted_list(self, mol, list, cache, wdiag, wij, wji, &
   & qdiag, qij, qji, bdiag, bij, bji, bL)
   !> EEQBC model type
   class(eeqbc_model), intent(in) :: self
   !> Molecular structure data
   type(structure_type), intent(in) :: mol
   !> Multicharge neighborlist type (upper triangle)
   type(csr_list), intent(in) :: list
   !> Multicharge cache
   type(mchrg_cache), intent(in) :: cache
   !> Diagonal weights of the CN derivatives, W(b, b)
   real(wp), intent(in) :: wdiag(:)
   !> Weights of the CN derivatives, W(i, j)
   real(wp), intent(in) :: wij(:)
   !> Weights of the CN derivatives, W(j, i)
   real(wp), intent(in) :: wji(:)
   !> Diagonal weights of the local charge derivatives, Q(b, b)
   real(wp), intent(in) :: qdiag(:)
   !> Weights of the local charge derivatives, Q(i, j)
   real(wp), intent(in) :: qij(:)
   !> Weights of the local charge derivatives, Q(j, i)
   real(wp), intent(in) :: qji(:)
   !> Derivative of each atom w.r.t. its own position, dB(b)/dR(b)
   real(wp), intent(inout) :: bdiag(:, :)
   !> Derivative of the neighbor w.r.t. the central atom, dB(j)/dR(i)
   real(wp), intent(inout) :: bij(:, :)
   !> Derivative of the central atom w.r.t. the neighbor, dB(i)/dR(j)
   real(wp), intent(inout) :: bji(:, :)
   !> Derivative w.r.t. strain deformations
   real(wp), intent(inout) :: bL(:, :, :)

   integer :: iat, jat, bat, itr
   integer(i8) :: kat, ent, mq, ip
   real(wp) :: vec(3), dpair(3, 2), qpair(3, 2), dcn(3, 2), dql(3, 2)
   real(wp) :: scn(3, 3, 2), sql(3, 3, 2), wpb, qpb, wqb, qqb, g(3)
   integer, allocatable :: row(:)
   integer(i8), allocatable :: adj_ptr(:), adj_ent(:)

   ! Thread-private arrays for reduction and neighborhood marks
   real(wp), allocatable :: bdiag_local(:, :), bij_local(:, :), bji_local(:, :)
   real(wp), allocatable :: bL_local(:, :, :)
   integer(i8), allocatable :: markp(:), markq(:)

   call get_list_adjacency(list, mol%nat, row, adj_ptr, adj_ent)

   !$omp parallel default(none) &
   !$omp shared(self, mol, list, cache, wdiag, wij, wji, qdiag, qij, qji) &
   !$omp shared(row, adj_ptr, adj_ent) &
   !$omp private(jat, bat, itr, kat, ent, mq, ip, vec, dpair, qpair, dcn, dql) &
   !$omp private(scn, sql, wpb, qpb, wqb, qqb, g, markp, markq) &
   !$omp shared(bdiag, bij, bji, bL) &
   !$omp private(bdiag_local, bij_local, bji_local, bL_local)
   allocate(bdiag_local(3, size(bdiag, 2)), source=0.0_wp)
   allocate(bij_local(3, size(bij, 2)), bji_local(3, size(bji, 2)), source=0.0_wp)
   allocate(bL_local(3, 3, size(bL, 3)), source=0.0_wp)
   allocate(markp(mol%nat), markq(mol%nat), source=0_i8)
   !$omp do schedule(runtime)
   do iat = 1, mol%nat
      ! Mark the neighborhood of iat with its list entries
      do ip = adj_ptr(iat), adj_ptr(iat + 1) - 1
         markp(list_partner(list, row, adj_ent(ip))) = adj_ent(ip)
      end do

      do kat = list%inl(iat), list%inl(iat + 1) - 1
         jat = list%nlat(kat)

         ! Pair derivatives summed over all lattice points
         dcn(:, :) = 0.0_wp
         dql(:, :) = 0.0_wp
         scn(:, :, :) = 0.0_wp
         sql(:, :, :) = 0.0_wp
         do itr = 1, size(cache%trans, 2)
            vec(:) = mol%xyz(:, iat) - (mol%xyz(:, jat) + cache%trans(:, itr))
            dpair(:, :) = get_dcndr_pair(self%ncoord, mol, iat, jat, vec, cache%cn)
            qpair(:, :) = get_dqlocdr_pair(self%ncoord_en, mol, iat, jat, vec, &
               & cache%qloc)
            dcn(:, :) = dcn + dpair
            dql(:, :) = dql + qpair
            scn(:, :, 1) = scn(:, :, 1) + spread(dpair(:, 1), 2, 3) * spread(vec, 1, 3)
            scn(:, :, 2) = scn(:, :, 2) + spread(dpair(:, 2), 2, 3) * spread(vec, 1, 3)
            sql(:, :, 1) = sql(:, :, 1) + spread(qpair(:, 1), 2, 3) * spread(vec, 1, 3)
            sql(:, :, 2) = sql(:, :, 2) + spread(qpair(:, 2), 2, 3) * spread(vec, 1, 3)
         end do

         ! Self-images only contribute to the strain derivative
         if (jat == iat) then
            bL_local(:, :, iat) = bL_local(:, :, iat) &
               & + wdiag(iat) * scn(:, :, 1) + qdiag(iat) * sql(:, :, 1)
            do ip = adj_ptr(iat), adj_ptr(iat + 1) - 1
               ent = adj_ent(ip)
               bat = list_partner(list, row, ent)
               bL_local(:, :, bat) = bL_local(:, :, bat) &
                  & + list_weight(ent, wij, wji) * scn(:, :, 1) &
                  & + list_weight(ent, qij, qji) * sql(:, :, 1)
            end do
            cycle
         end if

         ! Mark the neighborhood of jat with its list entries
         do ip = adj_ptr(jat), adj_ptr(jat + 1) - 1
            markq(list_partner(list, row, adj_ent(ip))) = adj_ent(ip)
         end do

         ! Columns b in the neighborhood of iat, starting with b = iat
         do ip = adj_ptr(iat) - 1, adj_ptr(iat + 1) - 1
            if (ip < adj_ptr(iat)) then
               bat = iat
               wpb = wdiag(iat)
               qpb = qdiag(iat)
            else
               ent = adj_ent(ip)
               bat = list_partner(list, row, ent)
               wpb = list_weight(ent, wij, wji)
               qpb = list_weight(ent, qij, qji)
            end if

            mq = 0_i8
            if (bat == jat) then
               wqb = wdiag(jat)
               qqb = qdiag(jat)
            else
               mq = markq(bat)
               wqb = list_weight(mq, wij, wji)
               qqb = list_weight(mq, qij, qji)
            end if

            g(:) = wpb * dcn(:, 1) + wqb * dcn(:, 2) + qpb * dql(:, 1) + qqb * dql(:, 2)

            ! d/dR(iat), always part of the pattern
            if (ip < adj_ptr(iat)) then
               bdiag_local(:, iat) = bdiag_local(:, iat) + g
            else if (ent > 0) then
               bij_local(:, ent) = bij_local(:, ent) + g
            else
               bji_local(:, -ent) = bji_local(:, -ent) + g
            end if

            ! d/dR(jat), only within the pattern
            if (bat == jat) then
               bdiag_local(:, jat) = bdiag_local(:, jat) - g
            else if (mq > 0) then
               bij_local(:, mq) = bij_local(:, mq) - g
            else if (mq < 0) then
               bji_local(:, -mq) = bji_local(:, -mq) - g
            end if

            bL_local(:, :, bat) = bL_local(:, :, bat) &
               & + wpb * scn(:, :, 1) + wqb * scn(:, :, 2) &
               & + qpb * sql(:, :, 1) + qqb * sql(:, :, 2)
         end do

         ! Columns b in the neighborhood of jat but not of iat
         do ip = adj_ptr(jat), adj_ptr(jat + 1) - 1
            ent = adj_ent(ip)
            bat = list_partner(list, row, ent)
            if (bat == iat .or. markp(bat) /= 0) cycle
            wqb = list_weight(ent, wij, wji)
            qqb = list_weight(ent, qij, qji)
            g(:) = wqb * dcn(:, 2) + qqb * dql(:, 2)
            if (ent > 0) then
               bij_local(:, ent) = bij_local(:, ent) - g
            else
               bji_local(:, -ent) = bji_local(:, -ent) - g
            end if
            bL_local(:, :, bat) = bL_local(:, :, bat) &
               & + wqb * scn(:, :, 2) + qqb * sql(:, :, 2)
         end do

         do ip = adj_ptr(jat), adj_ptr(jat + 1) - 1
            markq(list_partner(list, row, adj_ent(ip))) = 0_i8
         end do
      end do

      do ip = adj_ptr(iat), adj_ptr(iat + 1) - 1
         markp(list_partner(list, row, adj_ent(ip))) = 0_i8
      end do
   end do
   !$omp end do
   !$omp critical (add_cn_weighted_list_)
   bdiag(:, :) = bdiag + bdiag_local
   bij(:, :) = bij + bij_local
   bji(:, :) = bji + bji_local
   bL(:, :, :size(bL_local, 3)) = bL(:, :, :size(bL_local, 3)) + bL_local
   !$omp end critical (add_cn_weighted_list_)
   deallocate(bdiag_local, bij_local, bji_local, bL_local, markp, markq)
   !$omp end parallel

end subroutine add_cn_weighted_list

!> Partner atom of a signed list entry, +kat for the row atom of the entry
!> (partner is list%nlat(kat)) and -kat for the neighbor (partner is the row atom)
pure function list_partner(list, row, ent) result(bat)
   !> Multicharge neighborlist type
   type(csr_list), intent(in) :: list
   !> Row atom of each list entry
   integer, intent(in) :: row(:)
   !> Signed list entry
   integer(i8), intent(in) :: ent
   !> Partner atom
   integer :: bat

   if (ent > 0) then
      bat = list%nlat(ent)
   else
      bat = row(-ent)
   end if

end function list_partner

!> Weight of a signed list entry, ij(kat) for +kat, ji(kat) for -kat, zero for 0
pure function list_weight(ent, wij, wji) result(weight)
   !> Signed list entry
   integer(i8), intent(in) :: ent
   !> Weights of the entries (i, j)
   real(wp), intent(in) :: wij(:)
   !> Weights of the entries (j, i)
   real(wp), intent(in) :: wji(:)
   !> Weight of the entry
   real(wp) :: weight

   if (ent > 0) then
      weight = wij(ent)
   else if (ent < 0) then
      weight = wji(-ent)
   else
      weight = 0.0_wp
   end if

end function list_weight

!> Per-atom adjacency of an upper-triangle neighborlist, each off-diagonal entry
!> kat appears as +kat for its row atom and as -kat for its neighbor
subroutine get_list_adjacency(list, nat, row, adj_ptr, adj_ent)
   !> Multicharge neighborlist type (upper triangle)
   type(csr_list), intent(in) :: list
   !> Number of atoms
   integer, intent(in) :: nat
   !> Row atom of each list entry
   integer, allocatable, intent(out) :: row(:)
   !> Offset of each atom in the adjacency
   integer(i8), allocatable, intent(out) :: adj_ptr(:)
   !> Signed list entries of each atom
   integer(i8), allocatable, intent(out) :: adj_ent(:)

   integer :: iat, jat
   integer(i8) :: kat
   integer(i8), allocatable :: fill(:)

   allocate(row(size(list%nlat)), source=0)
   allocate(adj_ptr(nat + 1), source=0_i8)
   do iat = 1, nat
      do kat = list%inl(iat), list%inl(iat + 1) - 1
         row(kat) = iat
         jat = list%nlat(kat)
         if (jat == iat) cycle
         adj_ptr(iat + 1) = adj_ptr(iat + 1) + 1_i8
         adj_ptr(jat + 1) = adj_ptr(jat + 1) + 1_i8
      end do
   end do
   adj_ptr(1) = 1_i8
   do iat = 1, nat
      adj_ptr(iat + 1) = adj_ptr(iat + 1) + adj_ptr(iat)
   end do

   allocate(adj_ent(adj_ptr(nat + 1) - 1_i8))
   allocate(fill, source=adj_ptr(1:nat))
   do iat = 1, nat
      do kat = list%inl(iat), list%inl(iat + 1) - 1
         jat = list%nlat(kat)
         if (jat == iat) cycle
         adj_ent(fill(iat)) = kat
         fill(iat) = fill(iat) + 1_i8
         adj_ent(fill(jat)) = -kat
         fill(jat) = fill(jat) + 1_i8
      end do
   end do

end subroutine get_list_adjacency

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
         diag(jat) = diag(jat) + tmp
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
            diag(jat) = diag(jat) + tmp * wsw
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

         ! Updates for coordinate gradients (dcdrdiag)
         dcdrdiag(:, iat) = dcdrdiag(:, iat) - dG(:)
         dcdrdiag(:, jat) = dcdrdiag(:, jat) + dG(:)

         ! Updates for strain/lattice derivatives (dcdL)
         dcdL(:, :, iat) = dcdL(:, :, iat) + dS(:, :)
         dcdL(:, :, jat) = dcdL(:, :, jat) + dS(:, :)

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
            dcdrdiag_acc(:, jat) = dcdrdiag_acc(:, jat) + dG(:) * wsw

            dcdL_acc(:, :, iat) = dcdL_acc(:, :, iat) + dS(:, :) * wsw
            dcdL_acc(:, :, jat) = dcdL_acc(:, :, jat) + dS(:, :) * wsw
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

!> Compute the derivative of the coordination-number pair contribution for a single image.
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

   real(wp) :: den, countd(3), r1, r2

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
   real(wp) :: r1, dgtmpij(3), dgtmpji(3), vec(3)

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

!> Compute the coordination number and its diagonal derivatives using a neighborlist.
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
   !> Adjacency list for neighborlist-based CN evaluation
   type(csr_list), intent(in) :: list

   integer :: iat, jat, izp, jzp, itr
   integer(i8) :: kat
   real(wp) :: r2, r1, rij(3), countf, countd(3), sigma(3, 3), cutoff2, den

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

            ! Accumulate terms for the neighbor atom (j), avoiding
            ! double counting for self-images
            if (iat /= jat) then
               cn_local(jat) = cn_local(jat) + countf * self%directed_factor
               dcndrdiag_local(:,jat) = dcndrdiag_local(:, &
               & jat) - countd * self%directed_factor
               dcndL_local(:, :, jat) = dcndL_local(:, :, &
               & jat) + sigma * self%directed_factor
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

   deallocate(gradient_local, sigma_local, qlocacc, cnacc, dtrans, v)
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

   deallocate(gradient_local, sigma_local, qlocacc, cnacc, dtrans, v)
end subroutine get_grad_3d

end module multicharge_model_eeqbc
