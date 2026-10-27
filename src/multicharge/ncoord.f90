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
!> Provides the coordination number terms of the charge derivatives and the
!> Hessian on a sorted complete neighborlist

!> Coordination number terms on a sorted complete neighborlist, whose CN and
!> local charge derivatives are kept in the half storage of the upper
!> triangles of the rows. Without a maximum CN cutoff the entry of a pair holds
!> the derivatives of the CN of the atom of higher index, the ones of the atom
!> of lower index follow from dCN_i/dR_j = -dCN_j/dR_i and
!> dqloc_i/dR_j = dqloc_j/dR_i.
module multicharge_ncoord
   use mctc_csrlist, only : csr_list
   use multicharge_atomsort, only : get_pair_positions
   use multicharge_model_cache, only : mchrg_cache
   use multicharge_model_type, only : mchrg_hessian
   use mctc_env, only : wp, i8
   use mctc_io, only : structure_type
   use mctc_ncoord, only : ncoord_type
   implicit none
   private

   public :: add_ncoord_hessian_list, get_cn_entries_list

   !> Threshold of the CN and local charge derivatives on the neighborlist below
   !> which their products in the partial derivatives and the Hessian are skipped
   real(wp), parameter :: cnthr = epsilon(1.0_wp)**2

contains


!> Add dE/dCN contracted with the Cartesian Hessian of the coordination numbers
!> to the compressed Hessian, as add_coordination_number_hessian of mctc-lib
!> for a molecular system with a sorted complete neighborlist and without a
!> maximum CN cutoff. Each pair is evaluated once in the row of the atom of
!> lower index, which writes the block of the pair in the upper triangle of its
!> row of the Hessian, only the diagonal blocks are reduced over the threads.
subroutine add_ncoord_hessian_list(ncoord, mol, list, dEdcn, hessian)
   !> Coordination number container
   class(ncoord_type), intent(in) :: ncoord

   !> Molecular structure data
   type(structure_type), intent(in) :: mol

   !> Sorted complete neighborlist
   type(csr_list), intent(in) :: list

   !> Derivative of the expression w.r.t. the coordination numbers
   real(wp), intent(in) :: dEdcn(:)

   !> Molecular nuclear Hessian with allocated blocks
   type(mchrg_hessian), intent(inout) :: hessian

   integer :: iat, jat, izp, jzp, ic
   integer(i8) :: kat, hat, hoff
   real(wp) :: rij(3), r1, r2, cutoff2, den, dEdcnij, countd, countd2
   real(wp) :: box(3, 3)
   integer(i8), allocatable :: hmap(:)
   real(wp), allocatable :: diagonal(:, :, :)

   cutoff2 = ncoord%cutoff**2
   allocate(diagonal(3, 3, mol%nat), source=0.0_wp)

   ! The neighbors of higher index lead each row of the sorted complete list,
   ! hmap holds the position of the block of each atom in the upper triangle
   ! of the row of iat of the Hessian
   !$omp parallel default(none) &
   !$omp shared(ncoord, mol, list, dEdcn, hessian, cutoff2) &
   !$omp private(iat, jat, izp, jzp, ic, kat, hat, hoff, rij, r1, r2, den) &
   !$omp private(dEdcnij, countd, countd2, box, hmap) &
   !$omp reduction(+:diagonal)
   allocate(hmap(mol%nat), source=0_i8)
   !$omp do schedule(runtime)
   do iat = 1, mol%nat
      izp = mol%id(iat)
      hoff = hessian%list%inl(iat) - hessian%hinl(iat)
      do hat = hessian%hinl(iat), hessian%hinl(iat + 1) - 1
         hmap(hessian%list%nlat(hat + hoff)) = hat
      end do
      do kat = list%inl(iat) + 1, list%inl(iat + 1) - 1
         jat = list%nlat(kat)
         if (jat < iat) exit
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
         diagonal(:, :, iat) = diagonal(:, :, iat) + box
         diagonal(:, :, jat) = diagonal(:, :, jat) + box
         hat = hmap(jat)
         if (hat > 0) hessian%hesslist(:, :, hat) = hessian%hesslist(:, :, hat) - box
      end do
      do hat = hessian%hinl(iat), hessian%hinl(iat + 1) - 1
         hmap(hessian%list%nlat(hat + hoff)) = 0
      end do
   end do
   !$omp end do
   deallocate(hmap)
   !$omp end parallel

   !$omp parallel do default(none) schedule(static) &
   !$omp shared(mol, hessian, diagonal) &
   !$omp private(iat, hat)
   do iat = 1, mol%nat
      hat = hessian%hinl(iat)
      hessian%hesslist(:, :, hat) = hessian%hesslist(:, :, hat) + diagonal(:, :, iat)
   end do
   !$omp end parallel do

   deallocate(diagonal)

end subroutine add_ncoord_hessian_list

!> Neighbors of each atom with CN or local charge derivatives above cnthr on a
!> sorted complete neighborlist, the atom itself first. The rows cninl give
!> the neighbor atoms in cnat and the positions of the pairs in the half
!> storage in cnpos. Apart from the diagonal, the entry of a pair holds the
!> derivatives of the CN of the atom with the higher index.
subroutine get_cn_entries_list(mol, list, cache, cninl, cnpos, cnat)
   !> Molecular structure data
   type(structure_type), intent(in) :: mol

   !> Sorted complete neighborlist
   type(csr_list), intent(in) :: list

   !> Cache with the CN and local charge derivatives in the half storage
   type(mchrg_cache), intent(in) :: cache

   !> Offsets of the rows of the entries
   integer(i8), allocatable, intent(out) :: cninl(:)

   !> Positions of the pairs in the half storage
   integer(i8), allocatable, intent(out) :: cnpos(:)

   !> Neighbor atoms of the entries
   integer, allocatable, intent(out) :: cnat(:)

   integer :: iat, ient, nrow, kmax
   integer(i8) :: cnt
   integer(i8), allocatable :: hpos(:)

   kmax = int(maxval(list%inl(2:mol%nat + 1) - list%inl(:mol%nat)))
   allocate(cninl(mol%nat + 1))
   cninl(1) = 1

   !$omp parallel default(none) &
   !$omp shared(mol, list, cache, cninl, kmax) private(iat, ient, nrow, cnt, hpos)
   allocate(hpos(kmax))
   !$omp do schedule(runtime)
   do iat = 1, mol%nat
      nrow = int(list%inl(iat + 1) - list%inl(iat))
      call get_pair_positions(list, cache%hinl, iat, hpos(:nrow))
      cnt = 1
      do ient = 2, nrow
         if (any(abs(cache%dcndrlist(:, hpos(ient))) > cnthr) &
            & .or. any(abs(cache%dqlocdrlist(:, hpos(ient))) > cnthr)) cnt = cnt + 1
      end do
      cninl(iat + 1) = cnt
   end do
   !$omp end do
   deallocate(hpos)
   !$omp end parallel
   do iat = 1, mol%nat
      cninl(iat + 1) = cninl(iat) + cninl(iat + 1)
   end do

   allocate(cnpos(cninl(mol%nat + 1) - 1), cnat(cninl(mol%nat + 1) - 1))
   !$omp parallel default(none) &
   !$omp shared(mol, list, cache, cninl, cnpos, cnat, kmax) &
   !$omp private(iat, ient, nrow, cnt, hpos)
   allocate(hpos(kmax))
   !$omp do schedule(runtime)
   do iat = 1, mol%nat
      nrow = int(list%inl(iat + 1) - list%inl(iat))
      call get_pair_positions(list, cache%hinl, iat, hpos(:nrow))
      cnt = cninl(iat)
      cnpos(cnt) = hpos(1)
      cnat(cnt) = iat
      do ient = 2, nrow
         if (any(abs(cache%dcndrlist(:, hpos(ient))) > cnthr) &
            & .or. any(abs(cache%dqlocdrlist(:, hpos(ient))) > cnthr)) then
            cnt = cnt + 1
            cnpos(cnt) = hpos(ient)
            cnat(cnt) = list%nlat(list%inl(iat) + ient - 1)
         end if
      end do
   end do
   !$omp end do
   deallocate(hpos)
   !$omp end parallel

end subroutine get_cn_entries_list


end module multicharge_ncoord
