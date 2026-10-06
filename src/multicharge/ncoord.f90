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

!> @file multicharge/ncoord.f90
!> Provides on-the-fly pair derivatives of the coordination number and local charges

!> Wrapper around the coordination number evaluators returning the derivative
!> of a single pair contribution, to be called inside a loop over pairs.
!>
!> For a pair (iat, jat) and its distance vector vec = r(iat) - r(jat) - T the
!> result dcndr(3, 2) contains
!>
!> - dcndr(:, 1) = dCN(iat)/dvec
!> - dcndr(:, 2) = dCN(jat)/dvec
!>
!> The derivative w.r.t. r(iat) equals the one w.r.t. vec, the derivative w.r.t.
!> r(jat) is its negative, and the strain derivative is dcndr(:, k) x vec.
!> Both columns are needed since the partners are not symmetric: the directed
!> factor of the electronegativity weighted CN (local charges) gives them
!> opposite signs, and the smooth maximum CN cutoff scales each partner with the
!> slope at its own CN. For a self-image (iat == jat) the pair contributes to
!> CN(iat) only, the second column is zero.
!>
!> The pair derivatives are accumulated over all pairs and lattice points by
!> add_dcndr / add_dqlocdr (weighted dense derivatives) and get_dcndr_diag /
!> get_dqlocdr_diag (self derivatives and strain derivatives of each atom).
module multicharge_ncoord
   use mctc_env, only : wp
   use mctc_io, only : structure_type
   use mctc_ncoord, only : ncoord_type
   implicit none
   private

   public :: get_dcndr_pair, get_dqlocdr_pair
   public :: add_dcndr, add_dqlocdr, get_dcndr_diag, get_dqlocdr_diag


contains


!> Derivative of a single pair contribution to the coordination number
pure function get_dcndr_pair(ncoord, mol, iat, jat, vec, cn) result(dcndr)

   !> Coordination number container
   class(ncoord_type), intent(in) :: ncoord

   !> Molecular structure data
   type(structure_type), intent(in) :: mol

   !> Index of the first atom of the pair
   integer, intent(in) :: iat

   !> Index of the second atom of the pair
   integer, intent(in) :: jat

   !> Distance vector, vec = r(iat) - r(jat) - T
   real(wp), intent(in) :: vec(3)

   !> Coordination numbers including the maximum CN cutoff
   real(wp), intent(in) :: cn(:)

   !> Derivatives of CN(iat) and CN(jat) w.r.t. the distance vector
   real(wp) :: dcndr(3, 2)

   dcndr = get_pair_derivs(ncoord, mol, iat, jat, vec, cn(iat), cn(jat))

end function get_dcndr_pair


!> Derivative of a single pair contribution to the local charges
pure function get_dqlocdr_pair(ncoord_en, mol, iat, jat, vec, qloc) result(dqlocdr)

   !> Electronegativity weighted coordination number container
   class(ncoord_type), intent(in) :: ncoord_en

   !> Molecular structure data
   type(structure_type), intent(in) :: mol

   !> Index of the first atom of the pair
   integer, intent(in) :: iat

   !> Index of the second atom of the pair
   integer, intent(in) :: jat

   !> Distance vector, vec = r(iat) - r(jat) - T
   real(wp), intent(in) :: vec(3)

   !> Local charges including the equally distributed total charge
   real(wp), intent(in) :: qloc(:)

   !> Derivatives of qloc(iat) and qloc(jat) w.r.t. the distance vector
   real(wp) :: dqlocdr(3, 2)

   real(wp) :: qshift

   ! The cutoff acts on the weighted CN before the total charge is distributed
   qshift = mol%charge / real(mol%nat, wp)
   dqlocdr = get_pair_derivs(ncoord_en, mol, iat, jat, vec, &
      & qloc(iat) - qshift, qloc(jat) - qshift)

end function get_dqlocdr_pair


!> Add weighted CN derivatives, dxdr(:, :, a) += weight(a) * dCN(a)/dr
subroutine add_dcndr(ncoord, mol, trans, cn, weight, dxdr, dxdL)

   !> Coordination number container
   class(ncoord_type), intent(in) :: ncoord

   !> Molecular structure data
   type(structure_type), intent(in) :: mol

   !> Lattice points
   real(wp), intent(in) :: trans(:, :)

   !> Coordination numbers including the maximum CN cutoff
   real(wp), intent(in) :: cn(:)

   !> Weight of the CN derivative of each atom
   real(wp), intent(in) :: weight(:)

   !> Weighted derivative w.r.t. the Cartesian coordinates
   real(wp), intent(inout) :: dxdr(:, :, :)

   !> Weighted derivative w.r.t. strain deformations
   real(wp), intent(inout) :: dxdL(:, :, :)

   call add_pair_derivs(ncoord, mol, trans, cn, 0.0_wp, weight, dxdr, dxdL)

end subroutine add_dcndr


!> Add weighted local charge derivatives, dxdr(:, :, a) += weight(a) * dqloc(a)/dr
subroutine add_dqlocdr(ncoord_en, mol, trans, qloc, weight, dxdr, dxdL)

   !> Electronegativity weighted coordination number container
   class(ncoord_type), intent(in) :: ncoord_en

   !> Molecular structure data
   type(structure_type), intent(in) :: mol

   !> Lattice points
   real(wp), intent(in) :: trans(:, :)

   !> Local charges including the equally distributed total charge
   real(wp), intent(in) :: qloc(:)

   !> Weight of the local charge derivative of each atom
   real(wp), intent(in) :: weight(:)

   !> Weighted derivative w.r.t. the Cartesian coordinates
   real(wp), intent(inout) :: dxdr(:, :, :)

   !> Weighted derivative w.r.t. strain deformations
   real(wp), intent(inout) :: dxdL(:, :, :)

   call add_pair_derivs(ncoord_en, mol, trans, qloc, &
      & mol%charge / real(mol%nat, wp), weight, dxdr, dxdL)

end subroutine add_dqlocdr


!> CN derivative of each atom w.r.t. its own position and w.r.t. strain
subroutine get_dcndr_diag(ncoord, mol, trans, cn, dcndrdiag, dcndL)

   !> Coordination number container
   class(ncoord_type), intent(in) :: ncoord

   !> Molecular structure data
   type(structure_type), intent(in) :: mol

   !> Lattice points
   real(wp), intent(in) :: trans(:, :)

   !> Coordination numbers including the maximum CN cutoff
   real(wp), intent(in) :: cn(:)

   !> Derivative of CN(a) w.r.t. r(a)
   real(wp), intent(out) :: dcndrdiag(:, :)

   !> Derivative of CN(a) w.r.t. strain deformations
   real(wp), intent(out) :: dcndL(:, :, :)

   call get_diag_derivs(ncoord, mol, trans, cn, 0.0_wp, dcndrdiag, dcndL)

end subroutine get_dcndr_diag


!> Local charge derivative of each atom w.r.t. its own position and w.r.t. strain
subroutine get_dqlocdr_diag(ncoord_en, mol, trans, qloc, dqlocdrdiag, dqlocdL)

   !> Electronegativity weighted coordination number container
   class(ncoord_type), intent(in) :: ncoord_en

   !> Molecular structure data
   type(structure_type), intent(in) :: mol

   !> Lattice points
   real(wp), intent(in) :: trans(:, :)

   !> Local charges including the equally distributed total charge
   real(wp), intent(in) :: qloc(:)

   !> Derivative of qloc(a) w.r.t. r(a)
   real(wp), intent(out) :: dqlocdrdiag(:, :)

   !> Derivative of qloc(a) w.r.t. strain deformations
   real(wp), intent(out) :: dqlocdL(:, :, :)

   call get_diag_derivs(ncoord_en, mol, trans, qloc, &
      & mol%charge / real(mol%nat, wp), dqlocdrdiag, dqlocdL)

end subroutine get_dqlocdr_diag


!> Accumulate weighted pair derivatives over all pairs and lattice points
subroutine add_pair_derivs(ncoord, mol, trans, val, shift, weight, dxdr, dxdL)

   !> Coordination number container
   class(ncoord_type), intent(in) :: ncoord

   !> Molecular structure data
   type(structure_type), intent(in) :: mol

   !> Lattice points
   real(wp), intent(in) :: trans(:, :)

   !> Coordination numbers including the cutoff and the shift
   real(wp), intent(in) :: val(:)

   !> Constant shift added to the coordination numbers after the cutoff
   real(wp), intent(in) :: shift

   !> Weight of the derivative of each atom
   real(wp), intent(in) :: weight(:)

   !> Weighted derivative w.r.t. the Cartesian coordinates
   real(wp), intent(inout) :: dxdr(:, :, :)

   !> Weighted derivative w.r.t. strain deformations
   real(wp), intent(inout) :: dxdL(:, :, :)

   integer :: iat, jat, itr
   real(wp) :: vec(3), dpair(3, 2), gi(3), gj(3)

   ! Thread-private arrays for reduction
   ! Set to zero explicitly as the shared variants are potentially non-zero (inout)
   real(wp), allocatable :: dxdr_local(:, :, :), dxdL_local(:, :, :)

   !$omp parallel default(none) &
   !$omp shared(ncoord, mol, trans, val, shift, weight) &
   !$omp private(jat, itr, vec, dpair, gi, gj) &
   !$omp shared(dxdr, dxdL) &
   !$omp private(dxdr_local, dxdL_local)
   allocate(dxdr_local(size(dxdr, 1), size(dxdr, 2), size(dxdr, 3)), source=0.0_wp)
   allocate(dxdL_local(size(dxdL, 1), size(dxdL, 2), size(dxdL, 3)), source=0.0_wp)
   !$omp do schedule(runtime)
   do iat = 1, mol%nat
      do jat = 1, iat
         if (abs(weight(iat)) + abs(weight(jat)) <= 0.0_wp) cycle
         do itr = 1, size(trans, 2)
            vec(:) = mol%xyz(:, iat) - (mol%xyz(:, jat) + trans(:, itr))
            dpair(:, :) = get_pair_derivs(ncoord, mol, iat, jat, vec, &
               & val(iat) - shift, val(jat) - shift)
            gi(:) = weight(iat) * dpair(:, 1)
            gj(:) = weight(jat) * dpair(:, 2)

            dxdr_local(:, iat, iat) = dxdr_local(:, iat, iat) + gi
            dxdr_local(:, jat, iat) = dxdr_local(:, jat, iat) - gi
            dxdr_local(:, iat, jat) = dxdr_local(:, iat, jat) + gj
            dxdr_local(:, jat, jat) = dxdr_local(:, jat, jat) - gj

            dxdL_local(:, :, iat) = dxdL_local(:, :, iat) &
               & + spread(gi, 2, 3) * spread(vec, 1, 3)
            dxdL_local(:, :, jat) = dxdL_local(:, :, jat) &
               & + spread(gj, 2, 3) * spread(vec, 1, 3)
         end do
      end do
   end do
   !$omp end do
   !$omp critical (add_pair_derivs_)
   dxdr(:, :, :) = dxdr(:, :, :) + dxdr_local(:, :, :)
   dxdL(:, :, :) = dxdL(:, :, :) + dxdL_local(:, :, :)
   !$omp end critical (add_pair_derivs_)
   deallocate(dxdr_local, dxdL_local)
   !$omp end parallel

end subroutine add_pair_derivs


!> Accumulate the self derivatives and strain derivatives over all pairs
subroutine get_diag_derivs(ncoord, mol, trans, val, shift, ddiag, dL)

   !> Coordination number container
   class(ncoord_type), intent(in) :: ncoord

   !> Molecular structure data
   type(structure_type), intent(in) :: mol

   !> Lattice points
   real(wp), intent(in) :: trans(:, :)

   !> Coordination numbers including the cutoff and the shift
   real(wp), intent(in) :: val(:)

   !> Constant shift added to the coordination numbers after the cutoff
   real(wp), intent(in) :: shift

   !> Derivative of each atom w.r.t. its own position
   real(wp), intent(out) :: ddiag(:, :)

   !> Derivative of each atom w.r.t. strain deformations
   real(wp), intent(out) :: dL(:, :, :)

   integer :: iat, jat, itr
   real(wp) :: vec(3), dpair(3, 2)

   ! Thread-private arrays for reduction
   real(wp), allocatable :: ddiag_local(:, :), dL_local(:, :, :)

   ddiag(:, :) = 0.0_wp
   dL(:, :, :) = 0.0_wp

   !$omp parallel default(none) &
   !$omp shared(ncoord, mol, trans, val, shift) &
   !$omp private(jat, itr, vec, dpair) &
   !$omp shared(ddiag, dL) &
   !$omp private(ddiag_local, dL_local)
   allocate(ddiag_local(size(ddiag, 1), size(ddiag, 2)), source=0.0_wp)
   allocate(dL_local(size(dL, 1), size(dL, 2), size(dL, 3)), source=0.0_wp)
   !$omp do schedule(runtime)
   do iat = 1, mol%nat
      do jat = 1, iat
         do itr = 1, size(trans, 2)
            vec(:) = mol%xyz(:, iat) - (mol%xyz(:, jat) + trans(:, itr))
            dpair(:, :) = get_pair_derivs(ncoord, mol, iat, jat, vec, &
               & val(iat) - shift, val(jat) - shift)

            ! Self-images cancel in the position derivative
            if (iat /= jat) then
               ddiag_local(:, iat) = ddiag_local(:, iat) + dpair(:, 1)
               ddiag_local(:, jat) = ddiag_local(:, jat) - dpair(:, 2)
            end if

            dL_local(:, :, iat) = dL_local(:, :, iat) &
               & + spread(dpair(:, 1), 2, 3) * spread(vec, 1, 3)
            dL_local(:, :, jat) = dL_local(:, :, jat) &
               & + spread(dpair(:, 2), 2, 3) * spread(vec, 1, 3)
         end do
      end do
   end do
   !$omp end do
   !$omp critical (get_diag_derivs_)
   ddiag(:, :) = ddiag(:, :) + ddiag_local(:, :)
   dL(:, :, :) = dL(:, :, :) + dL_local(:, :, :)
   !$omp end critical (get_diag_derivs_)
   deallocate(ddiag_local, dL_local)
   !$omp end parallel

end subroutine get_diag_derivs


!> Derivative of a single pair contribution with the directed factor and the
!> chain rule of the maximum CN cutoff
pure function get_pair_derivs(ncoord, mol, iat, jat, vec, cni, cnj) result(dcndr)

   !> Coordination number container
   class(ncoord_type), intent(in) :: ncoord

   !> Molecular structure data
   type(structure_type), intent(in) :: mol

   !> Index of the first atom of the pair
   integer, intent(in) :: iat

   !> Index of the second atom of the pair
   integer, intent(in) :: jat

   !> Distance vector, vec = r(iat) - r(jat) - T
   real(wp), intent(in) :: vec(3)

   !> Coordination number of the first atom including the cutoff
   real(wp), intent(in) :: cni

   !> Coordination number of the second atom including the cutoff
   real(wp), intent(in) :: cnj

   !> Derivatives of CN(iat) and CN(jat) w.r.t. the distance vector
   real(wp) :: dcndr(3, 2)

   integer :: izp, jzp
   real(wp) :: r1, r2, countd(3), dcuti, dcutj

   dcndr(:, :) = 0.0_wp

   r2 = sum(vec**2)
   if (r2 > ncoord%cutoff**2 .or. r2 < 1.0e-12_wp) return
   r1 = sqrt(r2)

   izp = mol%id(iat)
   jzp = mol%id(jat)
   countd = ncoord%get_en_factor(izp, jzp) * ncoord%ncoord_dcount(izp, jzp, r1) &
      & * vec / r1

   dcuti = 1.0_wp
   dcutj = 1.0_wp
   if (ncoord%cut > 0.0_wp) then
      dcuti = dlog_cn_cut(cni, ncoord%cut)
      dcutj = dlog_cn_cut(cnj, ncoord%cut)
   end if

   dcndr(:, 1) = countd * dcuti
   if (iat /= jat) then
      dcndr(:, 2) = countd * ncoord%directed_factor * dcutj
   end if

end function get_pair_derivs


!> Slope of the smooth maximum CN cutoff, evaluated from the cut CN
!>
!> With cnp = log(1 + exp(cnmax)) - log(1 + exp(cnmax - cn)) the slope
!> exp(cnmax) / (exp(cnmax) + exp(cn)) equals 1 - exp(cnp) / (1 + exp(cnmax)).
elemental function dlog_cn_cut(cnp, cnmax) result(dcnpdcn)

   !> Coordination number including the cutoff
   real(wp), intent(in) :: cnp

   !> Maximum coordination number
   real(wp), intent(in) :: cnmax

   real(wp) :: dcnpdcn

   dcnpdcn = 1.0_wp - exp(cnp - log(1.0_wp + exp(cnmax)))

end function dlog_cn_cut


end module multicharge_ncoord
