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
!> Provides on-the-fly pair derivatives of the coordination number

!> Derivative of a single pair contribution to the coordination number, to be
!> called inside a loop over pairs.
!>
!> For a pair (iat, jat) and its distance vector vec = r(iat) - r(jat) - T the
!> results contain
!>
!> - dcndr(:, 1) = dCN(iat)/dvec, dcndL(:, :, 1) = dCN(iat)/dL
!> - dcndr(:, 2) = dCN(jat)/dvec, dcndL(:, :, 2) = dCN(jat)/dL
!>
!> The derivative w.r.t. r(iat) equals the one w.r.t. vec, the derivative w.r.t.
!> r(jat) is its negative. Both columns are needed since the partners are not
!> symmetric: the directed factor of the electronegativity weighted CN (local
!> charges) gives them opposite signs, and the smooth maximum CN cutoff scales
!> each partner with the slope at its own CN. For a self-image (iat == jat) the
!> pair contributes to CN(iat) only, the second column is zero.
!>
!> For the local charges the CN of the electronegativity weighted container
!> has to be passed, i.e. the local charges without the distributed total charge.
module multicharge_ncoord
   use mctc_env, only : wp
   use mctc_io, only : structure_type
   use mctc_ncoord, only : ncoord_type
   implicit none
   private

   public :: get_pair_derivs


contains


!> Derivative of a single pair contribution with the directed factor and the
!> chain rule of the maximum CN cutoff
pure subroutine get_pair_derivs(ncoord, mol, iat, jat, vec, cni, cnj, dcndr, dcndL)

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
   real(wp), intent(out) :: dcndr(3, 2)

   !> Derivatives of CN(iat) and CN(jat) w.r.t. strain deformations
   real(wp), intent(out) :: dcndL(3, 3, 2)

   integer :: izp, jzp
   real(wp) :: r1, r2, countd(3), sigma(3, 3), dcuti, dcutj

   dcndr(:, :) = 0.0_wp
   dcndL(:, :, :) = 0.0_wp

   r2 = sum(vec**2)
   if (r2 > ncoord%cutoff**2 .or. r2 < 1.0e-12_wp) return
   r1 = sqrt(r2)

   izp = mol%id(iat)
   jzp = mol%id(jat)
   countd = ncoord%get_en_factor(izp, jzp) * ncoord%ncoord_dcount(izp, jzp, r1) &
      & * vec / r1
   sigma = spread(countd, 1, 3) * spread(vec, 2, 3)

   dcuti = 1.0_wp
   dcutj = 1.0_wp
   if (ncoord%cut > 0.0_wp) then
      dcuti = dlog_cn_cut(cni, ncoord%cut)
      dcutj = dlog_cn_cut(cnj, ncoord%cut)
   end if

   dcndr(:, 1) = countd * dcuti
   dcndL(:, :, 1) = sigma * dcuti
   if (iat /= jat) then
      dcndr(:, 2) = countd * ncoord%directed_factor * dcutj
      dcndL(:, :, 2) = sigma * ncoord%directed_factor * dcutj
   end if

end subroutine get_pair_derivs


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
