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

!> @file multicharge/coulomb/shortrange.f90
!> Provides the sparse real-space part of the particle mesh Ewald operator

!> Real-space pair interactions of Gaussian charges in the CSR list layout.
!>
!> The pair kernel is erf(gam r)/r minus the smooth Ewald kernel
!> erf(alpha r)/r. Pair vectors follow the convention of the CSR list,
!> r_i - r_j - T, and the first entry of each row holds the diagonal.
module multicharge_coulomb_shortrange
   use mctc_csrlist, only : csr_list
   use mctc_env, only : wp, i8
   use mctc_io_constants, only : pi
   implicit none
   private

   public :: pair_kernel, get_gamma
   public :: get_pair_values, add_pair_gradient

   real(wp), parameter :: tsqrtpi = 2.0_wp / sqrt(pi)


contains


!> Pair-wise Gaussian width parameters gam_ij = 1/sqrt(rad_i^2 + rad_j^2)
pure subroutine get_gamma(rad, gam)
   !> Species-wise charge widths
   real(wp), intent(in) :: rad(:)
   !> Pair-wise width parameters
   real(wp), allocatable, intent(out) :: gam(:, :)

   integer :: isp, jsp

   allocate(gam(size(rad), size(rad)))
   do isp = 1, size(rad)
      do jsp = 1, size(rad)
         gam(jsp, isp) = 1.0_wp / sqrt(rad(isp)**2 + rad(jsp)**2)
      end do
   end do

end subroutine get_gamma


!> Value of the real-space kernel and its radial derivative divided by r
elemental subroutine pair_kernel(r, gam, alpha, f, df)
   !> Pair distance
   real(wp), intent(in) :: r
   !> Gaussian width parameter of the pair
   real(wp), intent(in) :: gam
   !> Ewald splitting parameter
   real(wp), intent(in) :: alpha
   !> Kernel value
   real(wp), intent(out) :: f
   !> Derivative of the kernel divided by the distance, f'(r)/r
   real(wp), intent(out) :: df

   real(wp) :: r2, egam, ealp

   r2 = r * r
   egam = tsqrtpi * gam * exp(-gam * gam * r2)
   ealp = tsqrtpi * alpha * exp(-alpha * alpha * r2)
   f = (erfc(alpha * r) - erfc(gam * r)) / r
   df = (egam - ealp - f) / r2

end subroutine pair_kernel


!> Evaluate the real-space matrix in the sparsity pattern of the CSR list
subroutine get_pair_values(list, xyz, id, gam, alpha, diag, val)
   !> Upper triangular CSR list, diagonal first in each row
   type(csr_list), intent(in) :: list
   !> Cartesian coordinates
   real(wp), intent(in) :: xyz(:, :)
   !> Species index of each atom
   integer, intent(in) :: id(:)
   !> Pair-wise Gaussian width parameters
   real(wp), intent(in) :: gam(:, :)
   !> Ewald splitting parameter
   real(wp), intent(in) :: alpha
   !> Diagonal elements
   real(wp), intent(in) :: diag(:)
   !> Matrix values in the CSR list layout
   real(wp), intent(out) :: val(:)

   integer :: iat, jat
   integer(i8) :: kat
   logical :: trlist
   real(wp) :: vec(3), r, f, df

   trlist = allocated(list%nltr)

   ! Every list entry belongs to exactly one row, no reduction is needed
   !$omp parallel do default(none) schedule(runtime) &
   !$omp shared(list, xyz, id, gam, alpha, diag, val, trlist) &
   !$omp private(iat, jat, kat, vec, r, f, df)
   do iat = 1, size(diag)
      val(list%inl(iat)) = diag(iat)
      do kat = list%inl(iat) + 1, list%inl(iat + 1) - 1
         jat = list%nlat(kat)
         vec = xyz(:, iat) - xyz(:, jat)
         if (trlist) vec = vec - list%trans(:, list%nltr(kat))
         r = norm2(vec)
         call pair_kernel(r, gam(id(jat), id(iat)), alpha, f, df)
         val(kat) = f
      end do
   end do

end subroutine get_pair_values


!> Add scale * sum_ij p_i q_j dS_ij/dR of the real-space matrix S
subroutine add_pair_gradient(list, xyz, id, gam, alpha, pvec, qvec, &
   & scale, gradient, sigma)
   !> Upper triangular CSR list, diagonal first in each row
   type(csr_list), intent(in) :: list
   !> Cartesian coordinates
   real(wp), intent(in) :: xyz(:, :)
   !> Species index of each atom
   integer, intent(in) :: id(:)
   !> Pair-wise Gaussian width parameters
   real(wp), intent(in) :: gam(:, :)
   !> Ewald splitting parameter
   real(wp), intent(in) :: alpha
   !> First charge-like vector
   real(wp), intent(in) :: pvec(:)
   !> Second charge-like vector
   real(wp), intent(in) :: qvec(:)
   !> Scaling factor of the contraction
   real(wp), intent(in) :: scale
   !> Cartesian gradient
   real(wp), intent(inout) :: gradient(:, :)
   !> Strain derivative
   real(wp), intent(inout) :: sigma(:, :)

   integer :: iat, jat, ic
   integer(i8) :: kat
   logical :: trlist
   real(wp) :: vec(3), r, f, df, weight, dG(3)

   ! Thread-private arrays for reduction
   ! Set to 0 explicitly as the shared variants are potentially non-zero (inout)
   real(wp), allocatable :: gradient_local(:, :), sigma_local(:, :)

   trlist = allocated(list%nltr)

   !$omp parallel default(none) &
   !$omp shared(list, xyz, id, gam, alpha, pvec, qvec, scale, trlist) &
   !$omp private(iat, jat, kat, ic, vec, r, f, df, weight, dG) &
   !$omp shared(gradient, sigma) &
   !$omp private(gradient_local, sigma_local)
   allocate(gradient_local(size(gradient, 1), size(gradient, 2)), source=0.0_wp)
   allocate(sigma_local(3, 3), source=0.0_wp)
   !$omp do schedule(runtime)
   do iat = 1, size(pvec)
      do kat = list%inl(iat) + 1, list%inl(iat + 1) - 1
         jat = list%nlat(kat)
         vec = xyz(:, iat) - xyz(:, jat)
         if (trlist) vec = vec - list%trans(:, list%nltr(kat))
         r = norm2(vec)
         call pair_kernel(r, gam(id(jat), id(iat)), alpha, f, df)
         if (iat == jat) then
            weight = pvec(iat) * qvec(iat)
         else
            weight = pvec(iat) * qvec(jat) + pvec(jat) * qvec(iat)
         end if
         dG = scale * weight * df * vec
         if (iat /= jat) then
            gradient_local(:, iat) = gradient_local(:, iat) + dG
            gradient_local(:, jat) = gradient_local(:, jat) - dG
         end if
         do ic = 1, 3
            sigma_local(:, ic) = sigma_local(:, ic) + vec * dG(ic)
         end do
      end do
   end do
   !$omp end do
   !$omp critical (add_pair_gradient_)
   gradient(:, :) = gradient(:, :) + gradient_local(:, :)
   sigma(:, :) = sigma(:, :) + sigma_local(:, :)
   !$omp end critical (add_pair_gradient_)
   deallocate(gradient_local, sigma_local)
   !$omp end parallel

end subroutine add_pair_gradient

end module multicharge_coulomb_shortrange
