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

!> @file multicharge/lowerlist.f90
!> Provides the lower triangle of a half neighborlist

!> Lower triangle of a half neighborlist, which stores each atom pair once in
!> the row of the atom with the lower index
module multicharge_lowerlist
   use mctc_env, only : i8
   use mctc_csrlist, only : csr_list
   !$ use omp_lib, only : omp_get_num_threads, omp_get_thread_num
   implicit none
   private

   public :: lower_list, new_lower_list, count_neighbors, get_neighbor


   !> Lower triangle of a half neighborlist. Row x holds the atoms b < x whose
   !> row of the list holds x, in ascending order, with the position of the
   !> entry of x in the row of b. Together with row x of the list it gives all
   !> neighbors of x.
   type :: lower_list

      !> Offset of each row, the last entry is the number of entries plus one
      integer(i8), allocatable :: inl(:)

      !> Atom b of each entry
      integer, allocatable :: nlat(:)

      !> Position of the entry of atom x in the row of atom b of the list
      integer(i8), allocatable :: pos(:)

   end type lower_list

contains


!> Create the lower triangle of a half neighborlist. Each thread counts the
!> entries of its rows of the list in each column, the entries of a thread
!> follow the ones of the threads with lower rows.
subroutine new_lower_list(self, nat, list)

   !> Lower triangle of the list
   type(lower_list), intent(out) :: self

   !> Number of atoms
   integer, intent(in) :: nat

   !> Half neighborlist with the diagonal and the atoms of higher index
   type(csr_list), intent(in) :: list

   integer :: bat, xat, ith, jth, nth, num
   integer(i8) :: kat, lat
   integer, allocatable :: cnt(:, :)

   allocate(self%inl(nat + 1))

   !$omp parallel default(none) &
   !$omp shared(self, nat, list, cnt, nth) &
   !$omp private(bat, xat, ith, jth, num, kat, lat)
   !$omp single
   nth = 1
   !$ nth = omp_get_num_threads()
   allocate(cnt(nat, 0:nth), source=0)
   !$omp end single
   ith = 1
   !$ ith = omp_get_thread_num() + 1
   !$omp do schedule(static)
   do bat = 1, nat
      do kat = list%inl(bat), list%inl(bat + 1) - 1
         xat = list%nlat(kat)
         if (xat == bat) cycle
         cnt(xat, ith) = cnt(xat, ith) + 1
      end do
   end do
   !$omp end do

   ! Offset of each thread in each row, cnt(:, 0) receives the row lengths
   !$omp do schedule(static)
   do xat = 1, nat
      do jth = 1, nth
         num = cnt(xat, jth)
         cnt(xat, jth) = cnt(xat, 0)
         cnt(xat, 0) = cnt(xat, 0) + num
      end do
   end do
   !$omp end do
   !$omp single
   self%inl(1) = 1
   do xat = 1, nat
      self%inl(xat + 1) = self%inl(xat) + cnt(xat, 0)
   end do
   allocate(self%nlat(self%inl(nat + 1) - 1), self%pos(self%inl(nat + 1) - 1))
   !$omp end single

   !$omp do schedule(static)
   do bat = 1, nat
      do kat = list%inl(bat), list%inl(bat + 1) - 1
         xat = list%nlat(kat)
         if (xat == bat) cycle
         lat = self%inl(xat) + cnt(xat, ith)
         cnt(xat, ith) = cnt(xat, ith) + 1
         self%nlat(lat) = bat
         self%pos(lat) = kat
      end do
   end do
   !$omp end do
   !$omp end parallel

end subroutine new_lower_list


!> Number of neighbors of an atom in a half neighborlist with its lower
!> triangle, including the atom itself
pure function count_neighbors(list, low, iat) result(nrow)

   !> Half neighborlist
   type(csr_list), intent(in) :: list

   !> Lower triangle of the list
   type(lower_list), intent(in) :: low

   !> Atom of the row
   integer, intent(in) :: iat

   !> Number of neighbors
   integer :: nrow

   nrow = int(list%inl(iat + 1) - list%inl(iat) + low%inl(iat + 1) - low%inl(iat))

end function count_neighbors


!> Neighbor of an atom in a half neighborlist with its lower triangle, the
!> entries of the row of the list come first, starting with the atom itself,
!> followed by the entries of the lower triangle
pure subroutine get_neighbor(list, low, iat, ient, jat, kat)

   !> Half neighborlist
   type(csr_list), intent(in) :: list

   !> Lower triangle of the list
   type(lower_list), intent(in) :: low

   !> Atom of the row
   integer, intent(in) :: iat

   !> Index of the neighbor in the row
   integer, intent(in) :: ient

   !> Neighbor atom
   integer, intent(out) :: jat

   !> Position of the pair in the list
   integer(i8), intent(out) :: kat

   integer(i8) :: lat

   kat = list%inl(iat) + ient - 1
   if (kat < list%inl(iat + 1)) then
      jat = list%nlat(kat)
   else
      lat = low%inl(iat) + (kat - list%inl(iat + 1))
      jat = low%nlat(lat)
      kat = low%pos(lat)
   end if

end subroutine get_neighbor

end module multicharge_lowerlist
