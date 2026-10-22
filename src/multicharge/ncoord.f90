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
   use mctc_csrlist, only : csr_list
   use multicharge_lowerlist, only : count_neighbors, get_neighbor
   use multicharge_model_cache, only : mchrg_cache
   use multicharge_model_type, only : mchrg_hessian
   use mctc_env, only : wp, i8
   use mctc_io, only : structure_type
   use mctc_ncoord, only : ncoord_type
   implicit none
   private

   public :: add_ncoord_hessian_list, get_cn_entries_list

   real(wp), parameter :: cnthr = epsilon(1.0_wp)**2

contains


!> Add dE/dCN contracted with the Cartesian Hessian of the coordination numbers
!> to the compressed Hessian on a half list, as add_coordination_number_hessian
!> of mctc-lib for a molecular system with a half neighborlist with its lower
!> triangle in the cache and without a maximum CN cutoff. Each iteration only
!> updates the row of atom iat of the Hessian, the diagonal block from all
!> pairs of iat and the blocks of the atoms of higher index.
subroutine add_ncoord_hessian_list(ncoord, mol, list, cache, dEdcn, hessian)
   !> Coordination number container
   class(ncoord_type), intent(in) :: ncoord

   !> Molecular structure data
   type(structure_type), intent(in) :: mol

   !> Multicharge neighborlist type (half)
   type(csr_list), intent(in) :: list

   !> Cache with the lower triangle of the list
   type(mchrg_cache), intent(in) :: cache

   !> Derivative of the expression w.r.t. the coordination numbers
   real(wp), intent(in) :: dEdcn(:)

   !> Molecular nuclear Hessian on a half list with allocated blocks
   type(mchrg_hessian), intent(inout) :: hessian

   integer :: iat, jat, izp, jzp, ic, ient
   integer(i8) :: kat, hat
   real(wp) :: rij(3), r1, r2, cutoff2, den, dEdcnij, countd, countd2
   real(wp) :: box(3, 3)
   integer(i8), allocatable :: hpos(:)

   cutoff2 = ncoord%cutoff**2

   !$omp parallel default(none) &
   !$omp shared(ncoord, mol, list, cache, dEdcn, hessian, cutoff2) &
   !$omp private(iat, jat, izp, jzp, ic, ient, kat, hat, rij, r1, r2, den) &
   !$omp private(dEdcnij, countd, countd2, box, hpos)
   allocate(hpos(mol%nat), source=0_i8)
   !$omp do schedule(runtime)
   do iat = 1, mol%nat
      izp = mol%id(iat)
      do hat = hessian%list%inl(iat), hessian%list%inl(iat + 1) - 1
         hpos(hessian%list%nlat(hat)) = hat
      end do
      do ient = 2, count_neighbors(list, cache%lowlist, iat)
         call get_neighbor(list, cache%lowlist, iat, ient, jat, kat)
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

!> Neighbors of each atom with CN or local charge derivatives above cnthr in a
!> half neighborlist with its lower triangle in the cache, the atom itself
!> first. The rows cninl give the neighbor atoms in cnat and the positions of
!> the pairs in the list in cnpos. Apart from the diagonal, the entry of a pair
!> holds the derivatives of the CN of the atom with the higher index.
subroutine get_cn_entries_list(mol, list, cache, cninl, cnpos, cnat)
   !> Molecular structure data
   type(structure_type), intent(in) :: mol

   !> Multicharge neighborlist type (half)
   type(csr_list), intent(in) :: list

   !> Cache with the CN and local charge derivatives and the lower triangle
   type(mchrg_cache), intent(in) :: cache

   !> Offsets of the rows of the entries
   integer(i8), allocatable, intent(out) :: cninl(:)

   !> Positions of the pairs in the list
   integer(i8), allocatable, intent(out) :: cnpos(:)

   !> Neighbor atoms of the entries
   integer, allocatable, intent(out) :: cnat(:)

   integer :: iat, jat, ient
   integer(i8) :: kat, cnt

   allocate(cninl(mol%nat + 1))
   cninl(1) = 1
   !$omp parallel do default(none) schedule(runtime) &
   !$omp shared(mol, list, cache, cninl) private(iat, jat, ient, kat, cnt)
   do iat = 1, mol%nat
      cnt = 1
      do ient = 2, count_neighbors(list, cache%lowlist, iat)
         call get_neighbor(list, cache%lowlist, iat, ient, jat, kat)
         if (any(abs(cache%dcndrlist(:, kat)) > cnthr) &
            & .or. any(abs(cache%dqlocdrlist(:, kat)) > cnthr)) cnt = cnt + 1
      end do
      cninl(iat + 1) = cnt
   end do
   !$omp end parallel do
   do iat = 1, mol%nat
      cninl(iat + 1) = cninl(iat) + cninl(iat + 1)
   end do

   allocate(cnpos(cninl(mol%nat + 1) - 1), cnat(cninl(mol%nat + 1) - 1))
   !$omp parallel do default(none) schedule(runtime) &
   !$omp shared(mol, list, cache, cninl, cnpos, cnat) private(iat, jat, ient, kat, cnt)
   do iat = 1, mol%nat
      cnt = cninl(iat)
      cnpos(cnt) = list%inl(iat)
      cnat(cnt) = iat
      do ient = 2, count_neighbors(list, cache%lowlist, iat)
         call get_neighbor(list, cache%lowlist, iat, ient, jat, kat)
         if (any(abs(cache%dcndrlist(:, kat)) > cnthr) &
            & .or. any(abs(cache%dqlocdrlist(:, kat)) > cnthr)) then
            cnt = cnt + 1
            cnpos(cnt) = kat
            cnat(cnt) = jat
         end if
      end do
   end do
   !$omp end parallel do

end subroutine get_cn_entries_list


end module multicharge_ncoord
