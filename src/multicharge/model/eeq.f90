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

!> @file multicharge/model/eeq.f90
!> Provides implementation of the electronegativity equilibration model (EEQ)

!> Electronegativity equilibration charge model published in
!>
!> E. Caldeweyher, S. Ehlert, A. Hansen, H. Neugebauer, S. Spicher, C. Bannwarth
!> and S. Grimme, *J. Chem. Phys.*, **2019**, 150, 154122.
!> DOI: [10.1063/1.5090222](https://dx.doi.org/10.1063/1.5090222)

module multicharge_model_eeq
   use mctc_env, only : error_type, fatal_error, wp
   use mctc_io, only : structure_type
   use mctc_io_constants, only : pi
   use mctc_io_math, only : matdet_3x3, outerprod
   use mctc_ncoord, only : new_ncoord, cn_count
   use mctc_wignerseitz, only : wignerseitz_cell
   use multicharge_wignerseitz, only : new_wignerseitz_cell
   use mctc_csrlist, only : csr_list
   use multicharge_blas, only : gemm
   use multicharge_linal, only : add_strain_block, add_strain_virial
   use multicharge_ewald, only : get_alpha
   use multicharge_model_type, only : mchrg_model_type, mchrg_hessian, get_dir_trans, &
      & get_rec_trans
   use multicharge_model_cache, only : mchrg_cache, clear_cache_derivs
   implicit none
   private

   public :: eeq_model, new_eeq_model

   !> EEQ model operations on shared charge-model parameters
   type, extends(mchrg_model_type) :: eeq_model
   contains
      !> Update and allocate cache
      procedure :: update
      !> Calculate capacitance matrix
      procedure :: get_capacitance_matrix
      !> Calculate capacitance matrix derivatives
      procedure :: get_capacitance_derivs
      !> Calculate Coulomb matrix
      procedure :: get_coulomb_matrix
      !> Add alpha * dA/dR*p
      procedure :: get_coulomb_derivs
      !> Calculate right-hand side (electronegativity vector)
      procedure :: get_xvec
      !> Add beta * dX/dR
      procedure :: get_xvec_derivs
      !> Calculate gradient
      procedure :: get_grad
      !> Explicit Hessian at a fixed lattice
      procedure :: get_explicit_hessian_dense
      !> Explicit Hessian on a neighborlist, not available for the EEQ model
      procedure :: get_explicit_hessian_list
   end type eeq_model

   real(wp), parameter :: sqrtpi = sqrt(pi)
   real(wp), parameter :: sqrt2pi = sqrt(2.0_wp / pi)
   real(wp), parameter :: eps = sqrt(epsilon(0.0_wp))

   real(wp), parameter :: cutoff = 15.0_wp

   !> Default scaling factor for the external electric field
   real(wp), parameter :: default_efield_scale = 1.0_wp


contains


!> Construct an EEQ model from element-wise parameters
subroutine new_eeq_model(self, mol, error, chi, rad, eta, kcnchi, &
   & cutoff, cn_exp, rcov, cn_max, efield_scale)
   !> Electronegativity equilibration model
   type(eeq_model), intent(out) :: self
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
   !> Cutoff radius for coordination number
   real(wp), intent(in), optional :: cutoff
   !> Steepness of the CN counting function
   real(wp), intent(in), optional :: cn_exp
   !> Covalent radii for CN
   real(wp), intent(in), optional :: rcov(:)
   !> Maximum CN cutoff for CN
   real(wp), intent(in), optional :: cn_max
   !> Scaling factor for the external electric field
   real(wp), intent(in), optional :: efield_scale


   self%chi = chi
   self%rad = rad
   self%eta = eta
   self%kcnchi = kcnchi

   if (present(efield_scale)) then
      self%efield_scale = efield_scale
   else
      self%efield_scale = default_efield_scale
   end if

   call new_ncoord(self%ncoord, mol, cn_count%erf, error, &
      & cutoff=cutoff, kcn=cn_exp, rcov=rcov, cut=cn_max)

end subroutine new_eeq_model


!> Update coordination numbers, qcp adds their missing derivatives
subroutine update(self, mol, cache, trans, list, qcp)
   !> EEQ model type
   class(eeq_model), intent(in) :: self
   !> Structure type
   type(structure_type), intent(in) :: mol
   !> Multicharge cache
   type(mchrg_cache), intent(inout) :: cache
   !> Lattice vectors
   real(wp), intent(in) :: trans(:, :)
   !> Multicharge neighborlist type
   type(csr_list), intent(in), optional :: list
   !> Flag for the CN derivatives of the charge response
   logical, intent(in), optional :: qcp

   logical :: deriv, fresh

   deriv = .false.
   if (present(qcp)) deriv = qcp

   ! A plain update starts over for a new structure
   if (.not. deriv) call clear_cache_derivs(cache)
   fresh = .not. (deriv .and. allocated(cache%cn))

   if (.not. allocated(cache%cn)) allocate(cache%cn(mol%nat))
   cache%trans = trans

   if (deriv .and. .not. allocated(cache%dcndr)) then
      allocate(cache%dcndr(3, mol%nat, mol%nat), cache%dcndL(3, 3, mol%nat))
      call self%ncoord%get_coordination_number(mol, trans, cache%cn, &
         & dcndr=cache%dcndr, dcndL=cache%dcndL)
   else if (fresh) then
      call self%ncoord%get_coordination_number(mol, trans, cache%cn)
   end if

   if (fresh .and. any(mol%periodic)) then
      call new_wignerseitz_cell(cache%wsc, mol)
      call get_alpha(mol%lattice, cache%alpha)
   end if

end subroutine update


!> No-op override: the EEQ model does not use a capacitance matrix.
subroutine get_capacitance_matrix(self, mol, ndim, cache, list)
   !> EEQ model type
   class(eeq_model), intent(in) :: self
   !> Structure type
   type(structure_type), intent(in) :: mol
   !> System size
   integer, intent(in) :: ndim
   !> Multicharge cache
   type(mchrg_cache), intent(inout) :: cache
   !> Multicharge neighborlist type
   type(csr_list), intent(in), optional :: list
end subroutine get_capacitance_matrix


!> No-op override: the EEQ model does not use a capacitance matrix.
subroutine get_capacitance_derivs(self, mol, cache, list)
   !> EEQ model type
   class(eeq_model), intent(in) :: self
   !> Structure type
   type(structure_type), intent(in) :: mol
   !> Multicharge cache
   type(mchrg_cache), intent(inout) :: cache
   !> Multicharge neighborlist type
   type(csr_list), intent(in), optional :: list
end subroutine get_capacitance_derivs


!> Build the electronegativity vector with CN correction.
subroutine get_xvec(self, mol, ndim, cache, list, efield)
   !> EEQ model type
   class(eeq_model), intent(in) :: self
   !> Structure type
   type(structure_type), intent(in) :: mol
   !> System size (number of atoms or atoms+1 if a Lagrange multiplier is used)
   integer, intent(in) :: ndim
   !> Multicharge cache (provides CN and will store the vector)
   type(mchrg_cache), intent(inout) :: cache
   !> Multicharge neighborlist type
   type(csr_list), intent(in), optional :: list
   !> External electric field
   real(wp), intent(in), optional :: efield(:)

   real(wp), parameter :: reg = 1.0e-14_wp

   integer :: iat, izp
   real(wp) :: tmp

   if (.not. allocated(cache%xvec)) then
      allocate(cache%xvec(ndim))
   else if (size(cache%xvec) /= ndim) then
      deallocate(cache%xvec)
      allocate(cache%xvec(ndim))
   end if

   !$omp parallel do default(none) schedule(runtime) &
   !$omp shared(mol, self, cache) private(iat, izp, tmp)
   do iat = 1, mol%nat
      izp = mol%id(iat)
      tmp = self%kcnchi(izp) / sqrt(cache%cn(iat) + reg)
      cache%xvec(iat) = -self%chi(izp) + tmp * cache%cn(iat)
   end do

   ! Add external electric field to the RHS if present
   if (present(efield)) then
      !$omp parallel do default(none) schedule(runtime) &
      !$omp shared(mol, self, cache, efield) private(iat)
      do iat = 1, mol%nat
         cache%xvec(iat) = cache%xvec(iat) + self%efield_scale &
            & * dot_product(mol%xyz(:, iat), efield)
      end do
   end if

   if (ndim == mol%nat + 1) then
      cache%xvec(mol%nat + 1) = mol%charge
   end if

end subroutine get_xvec


!> Assemble the Coulomb matrix (periodic or non‑periodic).
subroutine get_coulomb_matrix(self, mol, ndim, cache, list)
   !> EEQ model type
   class(eeq_model), intent(in) :: self
   !> Structure type
   type(structure_type), intent(in) :: mol
   !> System size
   integer, intent(in) :: ndim
   !> Multicharge cache (will hold the Coulomb matrix)
   type(mchrg_cache), intent(inout) :: cache
   !> Multicharge neighborlist type
   type(csr_list), intent(in), optional :: list

   if (.not. allocated(cache%amat)) then
      allocate(cache%amat(ndim, ndim))
   else if (size(cache%amat, 1) /= ndim) then
      deallocate(cache%amat)
      allocate(cache%amat(ndim, ndim))
   end if

   if (any(mol%periodic)) then
      call get_amat_3d(self, mol, cache%wsc, cache%alpha, cache%amat)
   else
      call get_amat_0d(self, mol, cache%amat)
   end if
end subroutine get_coulomb_matrix


!> Build the Coulomb matrix for a non‑periodic system (0D).
subroutine get_amat_0d(self, mol, amat)
   !> EEQ model type
   class(eeq_model), intent(in) :: self
   !> Molecular structure data
   type(structure_type), intent(in) :: mol
   !> Output Coulomb matrix (size ndim × ndim)
   real(wp), intent(out) :: amat(:, :)

   integer :: iat, jat, izp, jzp
   real(wp) :: vec(3), r2, gam, tmp

   real(wp), allocatable :: amat_local(:, :)

   amat(:, :) = 0.0_wp

   !$omp parallel default(none) &
   !$omp shared(amat, mol, self) &
   !$omp private(iat, izp, jat, jzp, gam, vec, r2, tmp, amat_local)
   allocate(amat_local, source=amat)
   !$omp do schedule(runtime)
   do iat = 1, mol%nat
      izp = mol%id(iat)
      do jat = 1, iat - 1
         jzp = mol%id(jat)
         vec = mol%xyz(:, jat) - mol%xyz(:, iat)
         r2 = vec(1)**2 + vec(2)**2 + vec(3)**2
         gam = 1.0_wp / (self%rad(izp)**2 + self%rad(jzp)**2)
         tmp = erf(sqrt(r2 * gam)) / sqrt(r2)
         amat_local(jat, iat) = amat_local(jat, iat) + tmp
         amat_local(iat, jat) = amat_local(iat, jat) + tmp
      end do
      tmp = self%eta(izp) + sqrt2pi / self%rad(izp)
      amat_local(iat, iat) = amat_local(iat, iat) + tmp
   end do
   !$omp end do
   !$omp critical (get_amat_0d_)
   amat(:, :) = amat(:, :) + amat_local(:, :)
   !$omp end critical (get_amat_0d_)
   deallocate(amat_local)
   !$omp end parallel

   if (size(amat, 1) == mol%nat + 1) then
      amat(mol%nat + 1, 1:mol%nat + 1) = 1.0_wp
      amat(1:mol%nat + 1, mol%nat + 1) = 1.0_wp
      amat(mol%nat + 1, mol%nat + 1) = 0.0_wp
   end if

end subroutine get_amat_0d


!> Build the Coulomb matrix for a periodic system (3D) using Ewald summation.
subroutine get_amat_3d(self, mol, wsc, alpha, amat)
   !> EEQ model type
   class(eeq_model), intent(in) :: self
   !> Molecular structure data
   type(structure_type), intent(in) :: mol
   !> Wigner–Seitz cell for the given structure
   type(wignerseitz_cell), intent(in) :: wsc
   !> Ewald splitting parameter
   real(wp), intent(in) :: alpha
   !> Output Coulomb matrix (size ndim × ndim)
   real(wp), intent(out) :: amat(:, :)

   integer :: iat, jat, izp, jzp, img
   real(wp) :: vec(3), gam, wsw, dtmp, rtmp, vol
   real(wp), allocatable :: dtrans(:, :), rtrans(:, :)

   real(wp), allocatable :: amat_local(:, :)

   amat(:, :) = 0.0_wp

   vol = abs(matdet_3x3(mol%lattice))
   call get_dir_trans(mol, dtrans)
   call get_rec_trans(mol, rtrans)

   !$omp parallel default(none) &
   !$omp shared(amat, mol, self, wsc, dtrans, rtrans, alpha, vol) &
   !$omp private(iat, izp, jat, jzp, gam, wsw, vec, dtmp, rtmp, amat_local)
   allocate(amat_local, source=amat)
   !$omp do schedule(runtime)
   do iat = 1, mol%nat
      izp = mol%id(iat)
      do jat = 1, iat - 1
         jzp = mol%id(jat)
         gam = 1.0_wp / sqrt(self%rad(izp)**2 + self%rad(jzp)**2)
         wsw = 1.0_wp / real(wsc%nimg(jat, iat), wp)
         do img = 1, wsc%nimg(jat, iat)
            vec = mol%xyz(:, jat) - mol%xyz(:, iat) + wsc%trans(:, wsc%tridx(img, jat, &
               & iat))
            call get_amat_dir_3d(vec, gam, alpha, dtrans, dtmp)
            call get_amat_rec_3d(vec, vol, alpha, rtrans, rtmp)
            amat_local(jat, iat) = amat_local(jat, iat) + (dtmp + rtmp) * wsw
            amat_local(iat, jat) = amat_local(iat, jat) + (dtmp + rtmp) * wsw
         end do
      end do

      gam = 1.0_wp / sqrt(2.0_wp * self%rad(izp)**2)
      wsw = 1.0_wp / real(wsc%nimg(iat, iat), wp)
      do img = 1, wsc%nimg(iat, iat)
         vec = wsc%trans(:, wsc%tridx(img, iat, iat))
         call get_amat_dir_3d(vec, gam, alpha, dtrans, dtmp)
         call get_amat_rec_3d(vec, vol, alpha, rtrans, rtmp)
         amat_local(iat, iat) = amat_local(iat, iat) + (dtmp + rtmp) * wsw
      end do

      dtmp = self%eta(izp) + sqrt2pi / self%rad(izp) - 2.0_wp * alpha / sqrtpi
      amat_local(iat, iat) = amat_local(iat, iat) + dtmp
   end do
   !$omp end do
   !$omp critical (get_amat_3d_)
   amat(:, :) = amat(:, :) + amat_local(:, :)
   !$omp end critical (get_amat_3d_)
   deallocate(amat_local)
   !$omp end parallel

   if (size(amat, 1) == mol%nat + 1) then
      amat(mol%nat + 1, 1:mol%nat + 1) = 1.0_wp
      amat(1:mol%nat + 1, mol%nat + 1) = 1.0_wp
      amat(mol%nat + 1, mol%nat + 1) = 0.0_wp
   end if

end subroutine get_amat_3d


!> Real‑space contribution to the Coulomb matrix (direct sum).
subroutine get_amat_dir_3d(rij, gam, alp, trans, amat)
   !> Distance vector between two atoms (including lattice translation)
   real(wp), intent(in) :: rij(3)
   !> Gaussian width parameter
   real(wp), intent(in) :: gam
   !> Ewald splitting parameter
   real(wp), intent(in) :: alp
   !> Direct lattice translation vectors (3 × N)
   real(wp), intent(in) :: trans(:, :)
   !> Output contribution to the Coulomb matrix
   real(wp), intent(out) :: amat

   integer :: itr
   real(wp) :: vec(3), r1, tmp

   amat = 0.0_wp

   do itr = 1, size(trans, 2)
      vec(:) = rij + trans(:, itr)
      r1 = norm2(vec)
      if (r1 < eps) cycle
      tmp = erf(gam * r1) / r1 - erf(alp * r1) / r1
      amat = amat + tmp
   end do

end subroutine get_amat_dir_3d


!> Reciprocal‑space contribution to the Coulomb matrix (Ewald sum).
subroutine get_amat_rec_3d(rij, vol, alp, trans, amat)
   !> Distance vector between two atoms (including lattice translation)
   real(wp), intent(in) :: rij(3)
   !> Unit cell volume
   real(wp), intent(in) :: vol
   !> Ewald splitting parameter
   real(wp), intent(in) :: alp
   !> Reciprocal lattice translation vectors (3 × N)
   real(wp), intent(in) :: trans(:, :)
   !> Output contribution to the Coulomb matrix
   real(wp), intent(out) :: amat

   integer :: itr
   real(wp) :: fac, vec(3), g2, tmp

   amat = 0.0_wp
   fac = 4.0_wp * pi / vol

   do itr = 1, size(trans, 2)
      vec(:) = trans(:, itr)
      g2 = dot_product(vec, vec)
      if (g2 < eps) cycle
      tmp = cos(dot_product(rij, vec)) * fac * exp(-0.25_wp * g2 / (alp * alp)) / g2
      amat = amat + tmp
   end do

end subroutine get_amat_rec_3d


!> Add the electronegativity vector derivatives, beta * dX/dR
subroutine get_xvec_derivs(self, mol, ndim, cache, beta, list)
   !> EEQ model type
   class(eeq_model), intent(in) :: self
   !> Structure type
   type(structure_type), intent(in) :: mol
   !> System size
   integer, intent(in) :: ndim
   !> Multicharge cache (provides CNs, stores the derivatives)
   type(mchrg_cache), intent(inout) :: cache
   !> dX/dR multiplier
   real(wp), intent(in) :: beta
   !> Multicharge neighborlist type
   type(csr_list), intent(in), optional :: list

   real(wp), parameter :: reg = 1.0e-14_wp

   integer :: iat
   real(wp) :: dxdcn

   ! EN derivative through the CN, each iteration only updates column iat
   !$omp parallel do default(none) schedule(runtime) &
   !$omp shared(self, mol, cache, beta) private(iat, dxdcn)
   do iat = 1, mol%nat
      dxdcn = beta * 0.5_wp * self%kcnchi(mol%id(iat)) / sqrt(cache%cn(iat) + reg)
      cache%dabdr(:, :, iat) = cache%dabdr(:, :, iat) &
         & + dxdcn * cache%dcndr(:, :, iat)
      cache%dabdL(:, :, iat) = cache%dabdL(:, :, iat) &
         & + dxdcn * cache%dcndL(:, :, iat)
   end do
   !$omp end parallel do

end subroutine get_xvec_derivs


!> Add the Coulomb matrix derivatives contracted with a vector
subroutine get_coulomb_derivs(self, mol, ndim, cache, alpha, p, list, dadr, &
   & dadrlist, dadL)
   !> EEQ model type
   class(eeq_model), intent(in) :: self
   !> Structure type
   type(structure_type), intent(in) :: mol
   !> System size
   integer, intent(in) :: ndim
   !> Multicharge cache
   type(mchrg_cache), intent(inout) :: cache
   !> dA/dR*p multiplier
   real(wp), intent(in) :: alpha
   !> Vector contracted with the Coulomb matrix derivatives
   real(wp), intent(in) :: p(:)
   !> Geometrical derivatives of the Coulomb matrix multiplied p
   real(wp), intent(inout), optional :: dadr(:, :, :)
   !> Compressed geometrical derivatives of the Coulomb matrix multiplied by p
   real(wp), intent(inout), optional :: dadrlist(:, :)
   !> Lattice derivatives of the Coulomb matrix multiplied p
   real(wp), intent(inout) :: dadL(:, :, :)
   !> Multicharge neighborlist type (sorted complete)
   type(csr_list), intent(in), optional :: list

   if (any(mol%periodic)) then
      call get_damat_3d(self, mol, cache, alpha, p, dadr, dadL)
   else
      call get_damat_0d(self, mol, cache, alpha, p, dadr, dadL)
   end if

end subroutine get_coulomb_derivs


!> Add alpha * dA/dR*p for a non-periodic system
subroutine get_damat_0d(self, mol, cache, alpha, p, dadr, dadL)
   !> EEQ model type
   class(eeq_model), intent(in) :: self
   !> Molecular structure data
   type(structure_type), intent(in) :: mol
   !> Multicharge cache
   type(mchrg_cache), intent(in) :: cache
   !> dA/dR*p multiplier
   real(wp), intent(in) :: alpha
   !> Vector contracted with the Coulomb matrix derivatives
   real(wp), intent(in) :: p(:)
   !> Derivatives w.r.t. positions, dadr(:, x, b) of component b w.r.t. atom x
   real(wp), intent(inout) :: dadr(:, :, :)
   !> Derivatives w.r.t. strain
   real(wp), intent(inout) :: dadL(:, :, :)

   integer :: iat, jat, izp, jzp
   real(wp) :: vec(3), r2, gam, arg, dtmp, dG(3), dS(3, 3)

   ! Thread-private arrays for reduction
   real(wp), allocatable :: dadr_local(:, :, :), dadL_local(:, :, :)

   !$omp parallel default(none) &
   !$omp shared(dadr, dadL, cache, mol, self, alpha, p) &
   !$omp private(iat, izp, jat, jzp, gam, r2, vec, dG, dS, dtmp, arg) &
   !$omp private(dadr_local, dadL_local)
   allocate(dadr_local(3, mol%nat, size(dadr, 3)), source=0.0_wp)
   allocate(dadL_local(3, 3, size(dadL, 3)), source=0.0_wp)
   !$omp do schedule(runtime)
   do iat = 1, mol%nat
      izp = mol%id(iat)
      do jat = 1, iat - 1
         jzp = mol%id(jat)
         vec = mol%xyz(:, jat) - mol%xyz(:, iat)
         r2 = vec(1)**2 + vec(2)**2 + vec(3)**2
         gam = 1.0_wp / sqrt(self%rad(izp)**2 + self%rad(jzp)**2)
         arg = gam * gam * r2
         dtmp = 2.0_wp * gam * exp(-arg) / (sqrtpi * r2) - erf(sqrt(arg)) &
            & / (r2 * sqrt(r2))
         dG = alpha * dtmp * vec
         dS = outerprod(vec, dG)
         dadr_local(:, iat, iat) = dadr_local(:, iat, iat) - dG * p(jat)
         dadr_local(:, jat, jat) = dadr_local(:, jat, jat) + dG * p(iat)
         dadr_local(:, iat, jat) = dadr_local(:, iat, jat) - dG * p(iat)
         dadr_local(:, jat, iat) = dadr_local(:, jat, iat) + dG * p(jat)
         dadL_local(:, :, jat) = dadL_local(:, :, jat) + dS * p(iat)
         dadL_local(:, :, iat) = dadL_local(:, :, iat) + dS * p(jat)
      end do
   end do
   !$omp end do
   !$omp critical (get_damat_0d_)
   dadr(:, :, :) = dadr + dadr_local
   dadL(:, :, :) = dadL + dadL_local
   !$omp end critical (get_damat_0d_)
   deallocate(dadr_local, dadL_local)
   !$omp end parallel

end subroutine get_damat_0d


!> Add alpha * dA/dR*p for a periodic system
subroutine get_damat_3d(self, mol, cache, alpha, p, dadr, dadL)
   !> EEQ model type
   class(eeq_model), intent(in) :: self
   !> Molecular structure data
   type(structure_type), intent(in) :: mol
   !> Multicharge cache
   type(mchrg_cache), intent(in) :: cache
   !> dA/dR*p multiplier
   real(wp), intent(in) :: alpha
   !> Vector contracted with the Coulomb matrix derivatives
   real(wp), intent(in) :: p(:)
   !> Derivatives w.r.t. positions, dadr(:, x, b) of component b w.r.t. atom x
   real(wp), intent(inout) :: dadr(:, :, :)
   !> Derivatives w.r.t. strain
   real(wp), intent(inout) :: dadL(:, :, :)

   integer :: iat, jat, izp, jzp, img
   real(wp) :: vol, gam, wsw, vec(3), dG(3), dS(3, 3)
   real(wp) :: dGd(3), dSd(3, 3), dGr(3), dSr(3, 3)
   real(wp), allocatable :: dtrans(:, :), rtrans(:, :)

   ! Thread-private arrays for reduction
   real(wp), allocatable :: dadr_local(:, :, :), dadL_local(:, :, :)

   vol = abs(matdet_3x3(mol%lattice))
   call get_dir_trans(mol, dtrans)
   call get_rec_trans(mol, rtrans)

   !$omp parallel default(none) &
   !$omp shared(dadr, dadL, mol, self, cache, alpha, p, vol, dtrans, rtrans) &
   !$omp private(iat, izp, jat, jzp, img, gam, wsw, vec, dG, dS) &
   !$omp private(dGr, dSr, dGd, dSd, dadr_local, dadL_local)
   allocate(dadr_local(3, mol%nat, size(dadr, 3)), source=0.0_wp)
   allocate(dadL_local(3, 3, size(dadL, 3)), source=0.0_wp)
   !$omp do schedule(runtime)
   do iat = 1, mol%nat
      izp = mol%id(iat)
      do jat = 1, iat - 1
         jzp = mol%id(jat)
         dG(:) = 0.0_wp
         dS(:, :) = 0.0_wp
         gam = 1.0_wp / sqrt(self%rad(izp)**2 + self%rad(jzp)**2)
         wsw = 1.0_wp / real(cache%wsc%nimg(jat, iat), wp)
         do img = 1, cache%wsc%nimg(jat, iat)
            vec = mol%xyz(:, jat) - mol%xyz(:, iat) + cache%wsc%trans(:, &
               & cache%wsc%tridx(img, jat, iat))
            call get_damat_dir_3d(vec, gam, cache%alpha, dtrans, dGd, dSd)
            call get_damat_rec_3d(vec, vol, cache%alpha, rtrans, dGr, dSr)
            dG = dG + (dGd + dGr) * wsw
            dS = dS + (dSd + dSr) * wsw
         end do
         dG = alpha * dG
         dS = alpha * dS
         dadr_local(:, iat, iat) = dadr_local(:, iat, iat) - dG * p(jat)
         dadr_local(:, jat, jat) = dadr_local(:, jat, jat) + dG * p(iat)
         dadr_local(:, iat, jat) = dadr_local(:, iat, jat) - dG * p(iat)
         dadr_local(:, jat, iat) = dadr_local(:, jat, iat) + dG * p(jat)
         dadL_local(:, :, jat) = dadL_local(:, :, jat) + dS * p(iat)
         dadL_local(:, :, iat) = dadL_local(:, :, iat) + dS * p(jat)
      end do

      dS(:, :) = 0.0_wp
      gam = 1.0_wp / sqrt(2.0_wp * self%rad(izp)**2)
      wsw = 1.0_wp / real(cache%wsc%nimg(iat, iat), wp)
      do img = 1, cache%wsc%nimg(iat, iat)
         vec = cache%wsc%trans(:, cache%wsc%tridx(img, iat, iat))
         call get_damat_dir_3d(vec, gam, cache%alpha, dtrans, dGd, dSd)
         call get_damat_rec_3d(vec, vol, cache%alpha, rtrans, dGr, dSr)
         dS = dS + (dSd + dSr) * wsw
      end do
      dadL_local(:, :, iat) = dadL_local(:, :, iat) + alpha * dS * p(iat)
   end do
   !$omp end do
   !$omp critical (get_damat_3d_)
   dadr(:, :, :) = dadr + dadr_local
   dadL(:, :, :) = dadL + dadL_local
   !$omp end critical (get_damat_3d_)
   deallocate(dadr_local, dadL_local)
   !$omp end parallel

end subroutine get_damat_3d


!> Real‑space contribution to the Coulomb matrix derivatives.
subroutine get_damat_dir_3d(rij, gam, alp, trans, dg, ds)
   !> Distance vector between two atoms (including lattice translation)
   real(wp), intent(in) :: rij(3)
   !> Gaussian width parameter
   real(wp), intent(in) :: gam
   !> Ewald splitting parameter
   real(wp), intent(in) :: alp
   !> Direct lattice translation vectors (3 × N)
   real(wp), intent(in) :: trans(:, :)
   !> Derivative of the Coulomb matrix element w.r.t. atomic position (3)
   real(wp), intent(out) :: dg(3)
   !> Derivative of the Coulomb matrix element w.r.t. lattice parameters (3×3)
   real(wp), intent(out) :: ds(3, 3)

   integer :: itr
   real(wp) :: vec(3), r1, r2, gtmp, atmp, gam2, alp2

   dg(:) = 0.0_wp
   ds(:, :) = 0.0_wp

   gam2 = gam * gam
   alp2 = alp * alp

   do itr = 1, size(trans, 2)
      vec(:) = rij + trans(:, itr)
      r1 = norm2(vec)
      if (r1 < eps) cycle
      r2 = r1 * r1
      gtmp = +2.0_wp * gam * exp(-r2 * gam2) / (sqrtpi * r2) - erf(r1 * gam) / (r2 * r1)
      atmp = -2.0_wp * alp * exp(-r2 * alp2) / (sqrtpi * r2) + erf(r1 * alp) / (r2 * r1)
      dg(:) = dg + (gtmp + atmp) * vec
      ds(:, :) = ds + (gtmp + atmp) * outerprod(vec, vec)
   end do

end subroutine get_damat_dir_3d


!> Reciprocal‑space contribution to the Coulomb matrix derivatives.
subroutine get_damat_rec_3d(rij, vol, alp, trans, dg, ds)
   !> Distance vector between two atoms (including lattice translation)
   real(wp), intent(in) :: rij(3)
   !> Unit cell volume
   real(wp), intent(in) :: vol
   !> Ewald splitting parameter
   real(wp), intent(in) :: alp
   !> Reciprocal lattice translation vectors (3 × N)
   real(wp), intent(in) :: trans(:, :)
   !> Derivative of the Coulomb matrix element w.r.t. atomic position (3)
   real(wp), intent(out) :: dg(3)
   !> Derivative of the Coulomb matrix element w.r.t. lattice parameters (3×3)
   real(wp), intent(out) :: ds(3, 3)

   integer :: itr
   real(wp) :: fac, vec(3), g2, gv, etmp, dtmp, alp2
   real(wp), parameter :: unity(3, 3) = reshape(&
      & [1.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, 1.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, 1.0_wp], [3, 3])

   dg(:) = 0.0_wp
   ds(:, :) = 0.0_wp
   fac = 4.0_wp * pi / vol
   alp2 = alp * alp

   do itr = 1, size(trans, 2)
      vec(:) = trans(:, itr)
      g2 = dot_product(vec, vec)
      if (g2 < eps) cycle
      gv = dot_product(rij, vec)
      etmp = fac * exp(-0.25_wp * g2 / alp2) / g2
      dtmp = -sin(gv) * etmp
      dg(:) = dg + dtmp * vec
      ds(:, :) = ds + etmp * cos(gv) &
         & * ((2.0_wp / g2 + 0.5_wp / alp2) * outerprod(vec, vec) - unity)
   end do

end subroutine get_damat_rec_3d


!> Accumulate the gradient and stress contributions of the EEQ model
subroutine get_grad(self, mol, cache, p, gradient, sigma, alpha, beta, list)
   !> EEQ model type
   class(eeq_model), intent(in) :: self
   !> Structure type
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
   real(wp), intent(in), optional :: alpha
   !> Optional electronegativity prefactor
   real(wp), intent(in), optional :: beta
   !> neighborlist (each unordered pair appears once)
   type(csr_list), optional, intent(in) :: list

   if (any(mol%periodic)) then
      call get_grad_3d(self, mol, cache, p, gradient, sigma, alphain=alpha, betain=beta)
   else
      call get_grad_0d(self, mol, cache, p, gradient, sigma, alphain=alpha, betain=beta)
   end if

end subroutine get_grad


!> Explicit part of the nuclear Hessian at fixed p and charges q
subroutine get_explicit_hessian_dense(self, mol, cache, p, hessian, elst, error, &
   & alpha, beta)
   !> EEQ model type
   class(eeq_model), intent(in) :: self
   !> Structure type
   type(structure_type), intent(in) :: mol
   !> Multicharge cache
   type(mchrg_cache), intent(inout) :: cache
   !> Charge-like contraction vector
   real(wp), intent(in) :: p(:)
   !> Nuclear Hessian
   real(wp), intent(inout) :: hessian(:, :)
   !> Strain Hessian x -> (1 + eps) x, element (a, b, c, d) is eps_ab, eps_cd
   real(wp), intent(inout) :: elst(:, :, :, :)
   !> Error handling
   type(error_type), allocatable, intent(out) :: error
   !> Optional Coulomb matrix prefactor
   real(wp), intent(in), optional :: alpha
   !> Optional electronegativity prefactor
   real(wp), intent(in), optional :: beta

   real(wp), parameter :: reg = 1.0e-14_wp

   integer :: iat, izp, ii, nat
   real(wp) :: alp, bet, sqcn, sigma(3, 3), el4(3, 3, 3, 3), dl9(9)
   real(wp), allocatable :: dEdcn(:), d2Edcn2(:), work(:, :), diagonal(:, :, :)
   real(wp), allocatable :: gtmp(:, :)

   alp = 1.0_wp
   if (present(alpha)) alp = alpha
   bet = 1.0_wp
   if (present(beta)) bet = beta
   nat = mol%nat

   if (.not. allocated(cache%vrhs) .or. .not. allocated(cache%dcndr)) then
      call fatal_error(error, "Charges and CN derivatives are not in the cache, "// &
         & "solve and update the cache with qcp first")
      return
   end if

   allocate(dEdcn(nat), d2Edcn2(nat), work(3 * nat, nat))
   !$omp parallel do default(none) schedule(runtime) &
   !$omp shared(self, mol, cache, p, bet, dEdcn, d2Edcn2, work, nat) &
   !$omp private(iat, izp, sqcn)
   do iat = 1, nat
      izp = mol%id(iat)
      sqcn = sqrt(cache%cn(iat) + reg)
      dEdcn(iat) = bet * p(iat) * 0.5_wp * self%kcnchi(izp) / sqcn
      d2Edcn2(iat) = -bet * p(iat) * 0.25_wp * self%kcnchi(izp) / sqcn**3
      work(:, iat) = reshape(d2Edcn2(iat) * cache%dcndr(:, :, iat), [3 * nat])
   end do
   !$omp end parallel do
   el4(:, :, :, :) = 0.0_wp
   call self%ncoord%add_coordination_number_hessian(mol, cache%trans, dEdcn, &
      & hessian, el4)
   call gemm(work, reshape(cache%dcndr, [3 * nat, nat]), hessian, transb='t', &
      & beta=1.0_wp)
   deallocate(work)

   ! CN-CN strain term
   do iat = 1, nat
      dl9 = reshape(cache%dcndL(:, :, iat), [9])
      el4(:, :, :, :) = el4 + d2Edcn2(iat) * reshape(spread(dl9, 2, 9) &
         & * spread(dl9, 1, 9), [3, 3, 3, 3])
   end do

   allocate(diagonal(3, 3, nat), source=0.0_wp)
   if (any(mol%periodic)) then
      call get_d2amat_3d(self, mol, cache, alp, p, hessian, diagonal, el4)
   else
      call get_d2amat_0d(self, mol, cache, alp, p, hessian, diagonal, el4)
   end if

   !$omp parallel do default(none) schedule(static) &
   !$omp shared(nat, hessian, diagonal) private(iat, ii)
   do iat = 1, nat
      ii = 3 * (iat - 1)
      hessian(ii + 1:ii + 3, ii + 1:ii + 3) = hessian(ii + 1:ii + 3, ii + 1:ii + 3) &
         & + diagonal(:, :, iat)
   end do
   !$omp end parallel do

   ! Virial term of the strain Hessian, the CN part is set by the CN Hessian
   allocate(gtmp(3, nat), source=0.0_wp)
   sigma(:, :) = 0.0_wp
   call self%get_grad(mol, cache, p, gtmp, sigma, alpha=alp, beta=bet)
   do iat = 1, nat
      sigma(:, :) = sigma - dEdcn(iat) * cache%dcndL(:, :, iat)
   end do
   call add_strain_virial(sigma, el4)
   elst(:, :, :, :) = elst + el4

end subroutine get_explicit_hessian_dense


!> Add alpha * p^T d2A/dR2 q of the pair terms (molecular)
subroutine get_d2amat_0d(self, mol, cache, alpha, p, hessian, diagonal, elst)
   !> EEQ model type
   class(eeq_model), intent(in) :: self
   !> Molecular structure data
   type(structure_type), intent(in) :: mol
   !> Multicharge cache with the charges
   type(mchrg_cache), intent(in) :: cache
   !> d2A/dR2 multiplier
   real(wp), intent(in) :: alpha
   !> Charge-like contraction vector
   real(wp), intent(in) :: p(:)
   !> Nuclear Hessian, receives the off-diagonal blocks
   real(wp), intent(inout) :: hessian(:, :)
   !> Diagonal blocks of the Hessian
   real(wp), intent(inout) :: diagonal(:, :, :)
   !> Strain Hessian
   real(wp), intent(inout) :: elst(3, 3, 3, 3)

   integer :: iat, jat, izp, jzp
   real(wp) :: vec(3), r1, gam, dr, drr, weight, box(3, 3)
   real(wp) :: elst_local(3, 3, 3, 3)
   real(wp), allocatable :: diagonal_local(:, :, :)

   !$omp parallel default(none) &
   !$omp shared(self, mol, cache, alpha, p, hessian) &
   !$omp private(iat, jat, izp, jzp, vec, r1, gam, dr, drr, weight, box) &
   !$omp shared(diagonal, elst) &
   !$omp private(diagonal_local, elst_local)
   allocate(diagonal_local(3, 3, size(diagonal, 3)), source=0.0_wp)
   elst_local(:, :, :, :) = 0.0_wp
   !$omp do schedule(runtime)
   do iat = 1, mol%nat
      izp = mol%id(iat)
      do jat = 1, iat - 1
         jzp = mol%id(jat)
         vec = mol%xyz(:, iat) - mol%xyz(:, jat)
         r1 = norm2(vec)
         gam = 1.0_wp / sqrt(self%rad(izp)**2 + self%rad(jzp)**2)
         call get_coulomb_kernel(r1, gam, dr, drr)
         weight = alpha * (p(iat) * cache%vrhs(jat) + p(jat) * cache%vrhs(iat))
         box = weight * radial_hessian(vec, dr, drr)
         call add_pair_block(iat, jat, box, hessian, diagonal_local)
         call add_strain_block(box, vec, elst_local)
      end do
   end do
   !$omp end do

   !$omp critical (get_d2amat_0d_)
   diagonal(:, :, :) = diagonal + diagonal_local
   elst(:, :, :, :) = elst + elst_local
   !$omp end critical (get_d2amat_0d_)
   deallocate(diagonal_local)
   !$omp end parallel

end subroutine get_d2amat_0d


!> Add alpha * p^T d2A/dR2 q of the pair terms (periodic)
subroutine get_d2amat_3d(self, mol, cache, alpha, p, hessian, diagonal, elst)
   !> EEQ model type
   class(eeq_model), intent(in) :: self
   !> Molecular structure data
   type(structure_type), intent(in) :: mol
   !> Multicharge cache with the charges, the Wigner-Seitz cell and Ewald alpha
   type(mchrg_cache), intent(in) :: cache
   !> d2A/dR2 multiplier
   real(wp), intent(in) :: alpha
   !> Charge-like contraction vector
   real(wp), intent(in) :: p(:)
   !> Nuclear Hessian, receives the off-diagonal blocks
   real(wp), intent(inout) :: hessian(:, :)
   !> Diagonal blocks of the Hessian
   real(wp), intent(inout) :: diagonal(:, :, :)
   !> Strain Hessian
   real(wp), intent(inout) :: elst(3, 3, 3, 3)

   integer :: iat, jat, izp, jzp, img, itr
   real(wp) :: rij(3), vec(3), box(3, 3), rbox(3, 3), gam, r1, g2, wsw, weight, fac
   real(wp) :: dr, drr, ar, arr
   real(wp) :: elst_local(3, 3, 3, 3)
   real(wp), allocatable :: dtrans(:, :), rtrans(:, :), rweight(:)
   real(wp), allocatable :: diagonal_local(:, :, :), rcoef(:), rcoef_local(:)

   call get_dir_trans(mol, dtrans)
   call get_rec_trans(mol, rtrans)

   ! Reciprocal-space weights independent of the atom pair
   fac = 4.0_wp * pi / abs(matdet_3x3(mol%lattice))
   allocate(rweight(size(rtrans, 2)), rcoef(size(rtrans, 2)), source=0.0_wp)
   do itr = 1, size(rtrans, 2)
      g2 = dot_product(rtrans(:, itr), rtrans(:, itr))
      if (g2 < eps) cycle
      rweight(itr) = fac * exp(-0.25_wp * g2 / cache%alpha**2) / g2
   end do

   !$omp parallel default(none) &
   !$omp shared(self, mol, cache, alpha, p, hessian, dtrans, rtrans, rweight) &
   !$omp private(iat, jat, izp, jzp, img, itr, rij, vec, box, rbox, gam, r1, wsw) &
   !$omp private(weight, dr, drr, ar, arr) &
   !$omp shared(diagonal, elst, rcoef) &
   !$omp private(diagonal_local, elst_local, rcoef_local)
   allocate(diagonal_local(3, 3, size(diagonal, 3)), source=0.0_wp)
   allocate(rcoef_local(size(rtrans, 2)), source=0.0_wp)
   elst_local(:, :, :, :) = 0.0_wp
   !$omp do schedule(runtime)
   do iat = 1, mol%nat
      izp = mol%id(iat)
      do jat = 1, iat - 1
         jzp = mol%id(jat)
         gam = 1.0_wp / sqrt(self%rad(izp)**2 + self%rad(jzp)**2)
         wsw = 1.0_wp / real(cache%wsc%nimg(jat, iat), wp)
         weight = alpha * (p(iat) * cache%vrhs(jat) + p(jat) * cache%vrhs(iat))
         box(:, :) = 0.0_wp
         do img = 1, cache%wsc%nimg(jat, iat)
            rij = mol%xyz(:, jat) - mol%xyz(:, iat) &
               & + cache%wsc%trans(:, cache%wsc%tridx(img, jat, iat))
            do itr = 1, size(dtrans, 2)
               vec = rij + dtrans(:, itr)
               r1 = norm2(vec)
               if (r1 < eps) cycle
               call get_coulomb_kernel(r1, gam, dr, drr)
               call get_coulomb_kernel(r1, cache%alpha, ar, arr)
               rbox = wsw * radial_hessian(vec, dr - ar, drr - arr)
               box = box + rbox
               call add_strain_block(weight * rbox, vec, elst_local)
            end do
            do itr = 1, size(rtrans, 2)
               vec = rtrans(:, itr)
               box = box - wsw * rweight(itr) * cos(dot_product(rij, vec)) &
                  & * outerprod(vec, vec)
               rcoef_local(itr) = rcoef_local(itr) &
                  & + weight * wsw * cos(dot_product(rij, vec))
            end do
         end do
         call add_pair_block(iat, jat, weight * box, hessian, diagonal_local)
      end do

      ! Interaction with the own images depends on the strain only
      gam = 1.0_wp / sqrt(2.0_wp * self%rad(izp)**2)
      wsw = 1.0_wp / real(cache%wsc%nimg(iat, iat), wp)
      weight = alpha * p(iat) * cache%vrhs(iat)
      do img = 1, cache%wsc%nimg(iat, iat)
         rij = cache%wsc%trans(:, cache%wsc%tridx(img, iat, iat))
         do itr = 1, size(dtrans, 2)
            vec = rij + dtrans(:, itr)
            r1 = norm2(vec)
            if (r1 < eps) cycle
            call get_coulomb_kernel(r1, gam, dr, drr)
            call get_coulomb_kernel(r1, cache%alpha, ar, arr)
            call add_strain_block(weight * wsw * radial_hessian(vec, dr - ar, &
               & drr - arr), vec, elst_local)
         end do
         do itr = 1, size(rtrans, 2)
            rcoef_local(itr) = rcoef_local(itr) &
               & + weight * wsw * cos(dot_product(rij, rtrans(:, itr)))
         end do
      end do
   end do
   !$omp end do

   !$omp critical (get_d2amat_3d_)
   diagonal(:, :, :) = diagonal + diagonal_local
   elst(:, :, :, :) = elst + elst_local
   rcoef(:) = rcoef + rcoef_local
   !$omp end critical (get_d2amat_3d_)
   deallocate(diagonal_local, rcoef_local)
   !$omp end parallel

   ! Reciprocal terms depend on the strain through the volume and |G|^2 only
   do itr = 1, size(rtrans, 2)
      if (rweight(itr) == 0.0_wp) cycle
      call add_rec_strain(rtrans(:, itr), cache%alpha, rcoef(itr) * rweight(itr), elst)
   end do

end subroutine get_d2amat_3d


!> Strain Hessian of a reciprocal term w(G) = 4 pi / V exp(-G^2 / 4 alp^2) / G^2
pure subroutine add_rec_strain(gvec, alp, weight, elst)
   !> Reciprocal lattice vector
   real(wp), intent(in) :: gvec(3)
   !> Ewald splitting parameter
   real(wp), intent(in) :: alp
   !> Prefactor times w(G)
   real(wp), intent(in) :: weight
   !> Strain Hessian
   real(wp), intent(inout) :: elst(3, 3, 3, 3)

   integer :: ia, ib, ic, id
   real(wp) :: g2, u, c4, delta(3, 3)

   delta(:, :) = 0.0_wp
   do ia = 1, 3
      delta(ia, ia) = 1.0_wp
   end do
   g2 = dot_product(gvec, gvec)
   u = 0.25_wp / alp**2 + 1.0_wp / g2
   c4 = 4.0_wp * (u**2 + 1.0_wp / g2**2)
   do id = 1, 3
      do ic = 1, 3
         do ib = 1, 3
            do ia = 1, 3
               elst(ia, ib, ic, id) = elst(ia, ib, ic, id) + weight * ( &
                  & c4 * gvec(ia) * gvec(ib) * gvec(ic) * gvec(id) &
                  & - 2.0_wp * u * (delta(ib, id) * gvec(ia) * gvec(ic) &
                  & + delta(ia, id) * gvec(ib) * gvec(ic) &
                  & + delta(ib, ic) * gvec(ia) * gvec(id) &
                  & + delta(ic, id) * gvec(ia) * gvec(ib) &
                  & + delta(ia, ib) * gvec(ic) * gvec(id)) &
                  & + delta(ia, ib) * delta(ic, id) + delta(ia, id) * delta(ib, ic))
            end do
         end do
      end do
   end do

end subroutine add_rec_strain


!> First and second derivative of the Gaussian Coulomb kernel
pure subroutine get_coulomb_kernel(r1, gam, dr, drr)
   !> Distance
   real(wp), intent(in) :: r1
   !> Gaussian width parameter
   real(wp), intent(in) :: gam
   !> First derivative w.r.t. the distance
   real(wp), intent(out) :: dr
   !> Second derivative w.r.t. the distance
   real(wp), intent(out) :: drr

   real(wp) :: erfval, ex

   erfval = erf(gam * r1)
   ex = 2.0_wp / sqrtpi * exp(-(gam * r1)**2)
   dr = gam * ex / r1 - erfval / r1**2
   drr = -2.0_wp * gam**3 * ex - 2.0_wp * gam * ex / r1**2 + 2.0_wp * erfval / r1**3

end subroutine get_coulomb_kernel


!> Cartesian Hessian of a radial function
pure function radial_hessian(vec, dr, drr) result(box)
   !> Distance vector
   real(wp), intent(in) :: vec(3)
   !> First derivative w.r.t. the distance
   real(wp), intent(in) :: dr
   !> Second derivative w.r.t. the distance
   real(wp), intent(in) :: drr
   !> Cartesian Hessian
   real(wp) :: box(3, 3)

   integer :: ic
   real(wp) :: r1, evec(3)

   r1 = norm2(vec)
   evec = vec / r1
   box = (drr - dr / r1) * outerprod(evec, evec)
   do ic = 1, 3
      box(ic, ic) = box(ic, ic) + dr / r1
   end do

end function radial_hessian


!> Add the Hessian blocks of a pair term
pure subroutine add_pair_block(iat, jat, box, hessian, diagonal)
   !> First atom of the pair
   integer, intent(in) :: iat
   !> Second atom of the pair
   integer, intent(in) :: jat
   !> Cartesian Hessian of the pair term w.r.t. the distance vector
   real(wp), intent(in) :: box(3, 3)
   !> Nuclear Hessian, receives the off-diagonal blocks
   real(wp), intent(inout) :: hessian(:, :)
   !> Diagonal blocks of the Hessian
   real(wp), intent(inout) :: diagonal(:, :, :)

   integer :: ii, jj

   ii = 3 * (iat - 1)
   jj = 3 * (jat - 1)
   hessian(ii + 1:ii + 3, jj + 1:jj + 3) = hessian(ii + 1:ii + 3, jj + 1:jj + 3) - box
   hessian(jj + 1:jj + 3, ii + 1:ii + 3) = hessian(jj + 1:jj + 3, ii + 1:ii + 3) - box
   diagonal(:, :, iat) = diagonal(:, :, iat) + box
   diagonal(:, :, jat) = diagonal(:, :, jat) + box

end subroutine add_pair_block


!> The explicit part of the Hessian is not implemented for the EEQ model
subroutine get_explicit_hessian_list(self, mol, cache, p, hessian, elst, error, list, &
   & alpha, beta)
   !> EEQ model type
   class(eeq_model), intent(in) :: self
   !> Structure type
   type(structure_type), intent(in) :: mol
   !> Multicharge cache
   type(mchrg_cache), intent(inout) :: cache
   !> Charge-like contraction vector
   real(wp), intent(in) :: p(:)
   !> Molecular nuclear Hessian in compressed storage
   type(mchrg_hessian), intent(inout) :: hessian
   !> Strain Hessian
   real(wp), intent(inout) :: elst(:, :, :, :)
   !> Error handling
   type(error_type), allocatable, intent(out) :: error
   !> Multicharge neighborlist type (sorted complete)
   type(csr_list), intent(in) :: list
   !> Optional Coulomb matrix prefactor
   real(wp), intent(in), optional :: alpha
   !> Optional electronegativity prefactor
   real(wp), intent(in), optional :: beta

   call fatal_error(error, "Analytical Hessian is not available for the EEQ model")

end subroutine get_explicit_hessian_list


!> Accumulate gradient and stress contributions for a non-periodic system
subroutine get_grad_0d(self, mol, cache, p, gradient, sigma, alphain, betain)
   !> EEQ model type
   class(eeq_model), intent(in) :: self
   !> Structure type
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
   real(wp), intent(in), optional :: alphain
   !> Optional electronegativity prefactor
   real(wp), intent(in), optional :: betain

   integer :: iat, jat, izp, jzp
   real(wp) :: alpha, beta
   real(wp) :: vec(3), r2, gam, arg, dtmp, tmp
   real(wp) :: dG(3), dS(3, 3)
   real(wp) :: W_ij
   real(wp), allocatable :: cnacc(:), gradient_local(:, :), sigma_local(:, :)

   real(wp), parameter :: reg = 1.0e-14_wp

   alpha = 1.0_wp
   if (present(alphain)) alpha = alphain
   beta = 1.0_wp
   if (present(betain)) beta = betain
   allocate(cnacc(mol%nat), source=0.0_wp)

   !$omp parallel default(none) &
   !$omp shared(cache, mol, self, p, cnacc, gradient, sigma, alpha) &
   !$omp private(iat, izp, jat, jzp, gam, vec, r2, dtmp, arg, tmp) &
   !$omp private(W_ij, gradient_local, sigma_local, dG, dS)
   allocate(gradient_local(3, mol%nat))
   allocate(sigma_local(3, 3))
   gradient_local = 0.0_wp
   sigma_local = 0.0_wp

   !$omp do schedule(runtime)
   do iat = 1, mol%nat
      izp = mol%id(iat)
      do jat = 1, iat - 1
         jzp = mol%id(jat)
         vec = mol%xyz(:, jat) - mol%xyz(:, iat)
         r2 = dot_product(vec, vec)
         gam = 1.0_wp / sqrt(self%rad(izp)**2 + self%rad(jzp)**2)
         arg = gam * gam * r2
         dtmp = 2.0_wp * gam * exp(-arg) / (sqrtpi * r2) - erf(sqrt(arg)) &
            & / (r2 * sqrt(r2))
         dG = dtmp * vec
         dS = outerprod(vec, dG)

         ! Weights for the pair
         W_ij = p(iat) * cache%vrhs(jat) + p(jat) * cache%vrhs(iat)

         gradient_local(:, iat) = gradient_local(:, iat) - dG * W_ij * alpha
         gradient_local(:, jat) = gradient_local(:, jat) + dG * W_ij * alpha
         sigma_local(:, :)   = sigma_local(:, :)   + dS * W_ij * alpha

      end do

      tmp = self%kcnchi(izp) / sqrt(cache%cn(iat) + reg)
      cnacc(iat) = 0.5_wp * tmp * p(iat)

   end do
   !$omp end do

   !$omp critical
   gradient = gradient + gradient_local
   sigma = sigma + sigma_local
   !$omp end critical

   deallocate(gradient_local, sigma_local)
   !$omp end parallel

   allocate(gradient_local(3, mol%nat), source = 0.0_wp)
   allocate(sigma_local(3, 3), source = 0.0_wp)

   call self%ncoord%add_coordination_number_derivs(mol, cache%trans, cnacc, &
      & gradient_local, sigma_local)

   gradient = gradient + beta * gradient_local
   sigma = sigma + beta * sigma_local

end subroutine get_grad_0d


!> Accumulate gradient and stress contributions for a periodic system
subroutine get_grad_3d(self, mol, cache, p, gradient, sigma, alphain, betain)
   !> EEQ model type
   class(eeq_model), intent(in) :: self
   !> Structure type
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
   real(wp), intent(in), optional :: alphain
   !> Optional electronegativity prefactor
   real(wp), intent(in), optional :: betain

   integer :: iat, jat, izp, jzp, img
   real(wp) :: alpha, beta
   real(wp) :: vec(3), r2, gam, arg, dtmp, tmp
   real(wp) :: dG(3), dS(3, 3), dGd(3), dSd(3, 3), dGr(3), dSr(3, 3)
   real(wp) :: W_ij, W_ii, wsw, vol
   real(wp), allocatable :: cnacc(:), gradient_local(:, :), sigma_local(:, :)
   real(wp), allocatable :: dtrans(:, :), rtrans(:, :)

   real(wp), parameter :: reg = 1.0e-14_wp

   vol = abs(matdet_3x3(mol%lattice))
   call get_dir_trans(mol, dtrans)
   call get_rec_trans(mol, rtrans)

   alpha = 1.0_wp
   if (present(alphain)) alpha = alphain
   beta = 1.0_wp
   if (present(betain)) beta = betain
   allocate(cnacc(mol%nat), source=0.0_wp)

   !$omp parallel default(none) &
   !$omp shared(cache, mol, self, p, gradient, sigma, vol, dtrans, rtrans, alpha, cnacc) &
   !$omp private(iat, izp, jat, jzp, gam, vec, r2, dtmp, arg, wsw, tmp) &
   !$omp private(W_ii, W_ij, gradient_local, sigma_local) &
   !$omp private(dG, dS, dGd, dSd, dGr, dSr)
   allocate(gradient_local(3, mol%nat))
   allocate(sigma_local(3, 3))
   gradient_local = 0.0_wp
   sigma_local = 0.0_wp

   !$omp do schedule(runtime)
   do iat = 1, mol%nat
      izp = mol%id(iat)
      do jat = 1, iat - 1
         jzp = mol%id(jat)
         dG(:) = 0.0_wp
         dS(:, :) = 0.0_wp
         gam = 1.0_wp / sqrt(self%rad(izp)**2 + self%rad(jzp)**2)
         wsw = 1.0_wp / real(cache%wsc%nimg(jat, iat), wp)
         do img = 1, cache%wsc%nimg(jat, iat)
            vec = mol%xyz(:, jat) - mol%xyz(:, iat) + cache%wsc%trans(:, &
               & cache%wsc%tridx(img, jat, iat))
            call get_damat_dir_3d(vec, gam, cache%alpha, dtrans, dGd, dSd)
            call get_damat_rec_3d(vec, vol, cache%alpha, rtrans, dGr, dSr)
            dG = dG + (dGd + dGr) * wsw
            dS = dS + (dSd + dSr) * wsw
         end do
         ! Weights for the pair
         W_ij = p(iat) * cache%vrhs(jat) + p(jat) * cache%vrhs(iat)

         gradient_local(:, iat) = gradient_local(:, iat) - dG * W_ij * alpha
         gradient_local(:, jat) = gradient_local(:, jat) + dG * W_ij * alpha
         sigma_local(:, :) = sigma_local(:, :) + dS * W_ij * alpha
      end do
      dS(:, :) = 0.0_wp
      W_ii = p(iat) * cache%vrhs(iat)
      gam = 1.0_wp / sqrt(2.0_wp * self%rad(izp)**2)
      wsw = 1.0_wp / real(cache%wsc%nimg(iat, iat), wp)
      do img = 1, cache%wsc%nimg(iat, iat)
         vec = cache%wsc%trans(:, cache%wsc%tridx(img, iat, iat))
         call get_damat_dir_3d(vec, gam, cache%alpha, dtrans, dGd, dSd)
         call get_damat_rec_3d(vec, vol, cache%alpha, rtrans, dGr, dSr)
         dS = dS + (dSd + dSr) * wsw
      end do
      sigma_local(:, :) = sigma_local(:, :) + dS * W_ii * alpha

      tmp = self%kcnchi(izp) / sqrt(cache%cn(iat) + reg)
      cnacc(iat) = 0.5_wp * tmp * p(iat)

   end do
   !$omp end do

   !$omp critical
   gradient = gradient + gradient_local
   sigma = sigma + sigma_local
   !$omp end critical

   deallocate(gradient_local, sigma_local)
   !$omp end parallel

   allocate(gradient_local(3, mol%nat), source = 0.0_wp)
   allocate(sigma_local(3, 3), source = 0.0_wp)

   call self%ncoord%add_coordination_number_derivs(mol, cache%trans, cnacc, &
      & gradient_local, sigma_local)

   gradient = gradient + beta * gradient_local
   sigma = sigma + beta * sigma_local

end subroutine get_grad_3d

end module multicharge_model_eeq
