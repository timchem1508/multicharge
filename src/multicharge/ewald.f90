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

module multicharge_ewald
   use mctc_env, only : wp
   use mctc_io_constants, only : pi
   use mctc_io_math, only : matdet_3x3
   implicit none
   private

   public :: get_alpha, get_dir_cutoff, get_rec_cutoff

   real(wp), parameter :: sqrtpi = sqrt(pi)
   real(wp), parameter :: eps = sqrt(epsilon(0.0_wp))

   !> Scaling of the balanced Ewald parameter, shifting work to the reciprocal
   !> sum, whose terms are cheaper to evaluate than the real-space terms
   real(wp), parameter :: alpha_scale = 1.25_wp

   abstract interface
      !> Returns the max. value of a term in the real or reciprocal space part
      !> of the Ewald summation for a given vector length.
      pure function get_term_gen(dist, alpha, vol) result(val)
         import :: wp

         !> Length of the real or reciprocal space vector
         real(wp), intent(in) :: dist

         !> Parameter of the Ewald summation
         real(wp), intent(in) :: alpha

         !> Volume of the real space unit cell
         real(wp), intent(in) :: vol

         !> Value of the term
         real(wp) :: val

      end function get_term_gen
   end interface

contains


!> Ewald parameter balancing the cost of the real-space and reciprocal sums.
!> With both sums converged to a fixed tolerance per term, the number of
!> real-space terms scales as 1/(alpha^3 V) and the number of reciprocal terms
!> as alpha^3 V, which balance for alpha = sqrt(pi) / V^(1/3).
subroutine get_alpha(lattice, alpha)
   !> Lattice vectors
   real(wp), intent(in) :: lattice(:, :)
   !> Parameter of the Ewald summation
   real(wp), intent(out) :: alpha

   alpha = alpha_scale * sqrtpi / abs(matdet_3x3(lattice))**(1.0_wp / 3.0_wp)

end subroutine get_alpha


!> Real-space cutoff beyond which every term of the direct Ewald sum stays
!> below the given tolerance
pure function get_dir_cutoff(alpha, conv) result(cutoff)
   !> Parameter of the Ewald summation
   real(wp), intent(in) :: alpha
   !> Tolerance for a single real-space term
   real(wp), intent(in) :: conv
   !> Real-space cutoff
   real(wp) :: cutoff

   cutoff = search_cutoff(get_dir_term, alpha, 0.0_wp, conv)

end function get_dir_cutoff


!> Reciprocal-space cutoff beyond which every term of the reciprocal Ewald sum
!> stays below the given tolerance
pure function get_rec_cutoff(alpha, volume, conv) result(cutoff)
   !> Parameter of the Ewald summation
   real(wp), intent(in) :: alpha
   !> Volume of the real space unit cell
   real(wp), intent(in) :: volume
   !> Tolerance for a single reciprocal-space term
   real(wp), intent(in) :: conv
   !> Reciprocal-space cutoff
   real(wp) :: cutoff

   cutoff = search_cutoff(get_rec_term_3d, alpha, volume, conv)

end function get_rec_cutoff


!> Find the vector length beyond which a monotonically decaying Ewald term
!> stays below the given tolerance. The upper bound of the bisection is
!> returned, such that the term is guaranteed to be converged at the cutoff.
pure function search_cutoff(get_term, alpha, volume, conv) result(cutoff)
   !> Procedure to evaluate the Ewald term
   procedure(get_term_gen) :: get_term
   !> Parameter of the Ewald summation
   real(wp), intent(in) :: alpha
   !> Volume of the real space unit cell
   real(wp), intent(in) :: volume
   !> Tolerance for a single term
   real(wp), intent(in) :: conv
   !> Cutoff for the vector length
   real(wp) :: cutoff

   integer, parameter :: niter = 30
   real(wp) :: xl, xr, x
   integer :: iter

   xr = eps
   do while (get_term(xr, alpha, volume) > conv .and. xr <= huge(1.0_wp))
      xr = 2.0_wp * xr
   end do
   xl = 0.5_wp * xr

   do iter = 1, niter
      if (xr - xl <= eps * xr) exit
      x = 0.5_wp * (xl + xr)
      if (get_term(x, alpha, volume) > conv) then
         xl = x
      else
         xr = x
      end if
   end do
   cutoff = xr

end function search_cutoff


!> Returns the max. value of a term in the real space part of the Ewald summation
!> for a given vector length.
pure function get_dir_term(rr, alpha, vol) result(dval)

   !> Length of the real space vector
   real(wp), intent(in) :: rr

   !> Parameter of the Ewald summation
   real(wp), intent(in) :: alpha

   !> Volume of the real space unit cell, unused but required by the interface
   real(wp), intent(in) :: vol

   !> Real space term
   real(wp) :: dval

   dval = erfc(alpha*rr)/rr

end function get_dir_term


!> Returns the max. value of a term in the reciprocal space part of the Ewald
!> summation for a given vector length.
pure function get_rec_term_3d(gg, alpha, vol) result(rval)

   !> Length of the reciprocal space vector
   real(wp), intent(in) :: gg

   !> Parameter of the Ewald summation
   real(wp), intent(in) :: alpha

   !> Volume of the real space unit cell
   real(wp), intent(in) :: vol

   !> Reciprocal term
   real(wp) :: rval

   rval = 4.0_wp*pi*(exp(-0.25_wp*gg*gg/(alpha**2))/(vol*gg*gg))

end function get_rec_term_3d

end module multicharge_ewald
