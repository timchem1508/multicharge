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

!> @file multicharge/linal.f90
!> Provides dense and sparse linear algebra kernels of the charge models

!> Matrix-vector and matrix-matrix products on dense arrays and sorted lists
module multicharge_linal
   use mctc_env, only : wp, i8
   use mctc_csrlist, only : csr_list, find_lower_start
   use multicharge_blas, only : gemm
   implicit none
   private

   public :: max_group, add_transposed_product, gather_group_list, &
      & multiply_group_list, add_chain_products, get_transposed_positions


   !> Number of atoms in a group of the charge derivatives on a neighborlist
   integer, parameter :: max_group = 8


contains


!> Add the transposed product of partial derivatives on a list with a vector
subroutine add_transposed_product(list, dlist, vec, wvec)

   !> Sorted complete neighborlist
   type(csr_list), intent(in) :: list

   !> Partial derivatives on the list
   real(wp), intent(in) :: dlist(:, :)

   !> Vector contracted with the components
   real(wp), intent(in) :: vec(:)

   !> Product for each atom, dimension (3, nat)
   real(wp), intent(inout) :: wvec(:, :)

   integer :: xat
   integer(i8) :: kat
   real(wp) :: wx(3)

   !$omp parallel do default(none) schedule(runtime) &
   !$omp shared(list, dlist, vec, wvec) private(xat, kat, wx)
   do xat = 1, size(wvec, 2)
      wx(:) = 0.0_wp
      do kat = list%inl(xat), list%inl(xat + 1) - 1
         wx(:) = wx + dlist(:, kat) * vec(list%nlat(kat))
      end do
      wvec(:, xat) = wvec(:, xat) + wx
   end do

end subroutine add_transposed_product


!> Union of the neighbors of a group of atoms on a sorted complete list
subroutine gather_group_list(list, dlist, group, nb, bloc, blist, qpos, pan)

   !> Sorted complete neighborlist
   type(csr_list), intent(in) :: list

   !> Partial derivatives on the list, the row of x holds the column x
   real(wp), intent(in), contiguous :: dlist(:, :)

   !> Atoms of the group
   integer, intent(in), contiguous :: group(:)

   !> Number of atoms in the union of the neighbors of the group
   integer, intent(out) :: nb

   !> Index of each atom in the union, zero outside
   integer, intent(inout), contiguous :: bloc(:)

   !> Atoms of the union, sorted in ascending order
   integer, intent(inout), contiguous :: blist(:)

   !> Position of the entry of each union atom in the row of x
   integer(i8), intent(inout), contiguous :: qpos(:, :)

   !> Partial derivatives dB_c/dR_x of the union
   real(wp), intent(inout), contiguous :: pan(:, :)

   integer :: i, ic, xat, bat, n
   integer(i8) :: kat

   nb = 0
   do i = 1, size(group)
      xat = group(i)
      ic = 3 * (i - 1)
      do kat = list%inl(xat), list%inl(xat + 1) - 1
         bat = list%nlat(kat)
         n = bloc(bat)
         if (n == 0) then
            nb = nb + 1
            n = nb
            bloc(bat) = n
            blist(n) = bat
            qpos(:, n) = 0_i8
            pan(:, n) = 0.0_wp
         end if
         pan(ic + 1:ic + 3, n) = dlist(:, kat)
         qpos(i, n) = kat
      end do
   end do

end subroutine gather_group_list


!> Product of a symmetric matrix on a list with the derivatives of a union
pure subroutine multiply_group_list(list, ufin, smat, nb, bloc, blist, pan, qcol)

   !> Sorted complete list of the symmetric matrix
   type(csr_list), intent(in) :: list

   !> End of the upper triangle in each row of the list
   integer(i8), intent(in), contiguous :: ufin(:)

   !> Symmetric matrix on the upper triangle of the list
   real(wp), intent(in), contiguous :: smat(:)

   !> Number of atoms in the union
   integer, intent(in) :: nb

   !> Index of each atom in the union, zero outside
   integer, intent(in), contiguous :: bloc(:)

   !> Atoms of the union
   integer, intent(in), contiguous :: blist(:)

   !> Derivatives of the union with fixed leading dimension 3 * max_group
   real(wp), intent(in) :: pan(3 * max_group, *)

   !> Product with the symmetric matrix on the union
   real(wp), intent(inout) :: qcol(3 * max_group, *)

   integer, parameter :: nvec = 3 * max_group
   integer :: n, m, bat
   integer(i8) :: hat
   real(wp) :: s, pn(nvec), acc(nvec)

   qcol(:, :nb) = 0.0_wp
   do n = 1, nb
      bat = blist(n)
      pn(:) = pan(:, n)
      acc(:) = qcol(:, n) + smat(list%inl(bat)) * pn(:)
      do hat = list%inl(bat) + 1, ufin(bat)
         m = bloc(list%nlat(hat))
         if (m == 0) cycle
         s = smat(hat)
         acc(:) = acc(:) + s * pan(:, m)
         qcol(:, m) = qcol(:, m) + s * pn(:)
      end do
      qcol(:, n) = acc(:)
   end do

end subroutine multiply_group_list


!> Add D M^T + M D^T of derivatives and their weights to the Hessian
subroutine add_chain_products(dzdr, weight, hessian)

   !> Position derivatives of the atomic quantities, dzdr(:, x, k) = dz_k/dR_x
   real(wp), intent(in), contiguous, target :: dzdr(:, :, :)

   !> Vector weights of the position derivatives, in the layout of dzdr
   real(wp), intent(in), contiguous, target :: weight(:, :, :)

   !> Molecular nuclear Hessian
   real(wp), intent(inout) :: hessian(:, :)

   real(wp), pointer :: dptr(:, :), wptr(:, :)

   dptr(1:size(dzdr, 1) * size(dzdr, 2), 1:size(dzdr, 3)) => dzdr
   wptr(1:size(weight, 1) * size(weight, 2), 1:size(weight, 3)) => weight
   call gemm(dptr, wptr, hessian, transb='t', beta=1.0_wp)
   call gemm(wptr, dptr, hessian, transb='t', beta=1.0_wp)

end subroutine add_chain_products


!> Transposed positions of the entries of a sorted complete neighborlist
subroutine get_transposed_positions(list, tpos)

   !> Sorted complete neighborlist
   type(csr_list), intent(in) :: list

   !> Position of the atom in the row of the neighbor for each entry
   integer(i8), allocatable, intent(out) :: tpos(:)

   integer :: iat, jat, nat
   integer(i8) :: kat, lo, hi, mid, nfin

   nat = size(list%inl) - 1
   allocate(tpos(list%inl(nat + 1) - 1))

   !$omp parallel do default(none) schedule(runtime) &
   !$omp shared(list, tpos, nat) private(iat, jat, kat, lo, hi, mid, nfin)
   do iat = 1, nat
      tpos(list%inl(iat)) = list%inl(iat)
      nfin = find_lower_start(list, iat) - 1
      do kat = list%inl(iat) + 1, nfin
         jat = list%nlat(kat)
         ! First entry in the row of jat with a neighbor not above iat
         lo = list%inl(jat) + 1
         hi = list%inl(jat + 1) - 1
         do while (lo < hi)
            mid = (lo + hi) / 2
            if (list%nlat(mid) > iat) then
               lo = mid + 1
            else
               hi = mid
            end if
         end do
         tpos(kat) = lo
         tpos(lo) = kat
      end do
   end do
   !$omp end parallel do

end subroutine get_transposed_positions


end module multicharge_linal
