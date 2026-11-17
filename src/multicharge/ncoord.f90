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
!> Provides the coordination numbers and local charges of the charge models

!> Coordination numbers, local charges and their derivatives on neighborlists
module multicharge_ncoord
   use mctc_env, only : wp, i8
   use mctc_io, only : structure_type
   use mctc_ncoord, only : ncoord_type
   use mctc_csrlist, only : csr_list, find_lower_start
   use multicharge_model_cache, only : mchrg_cache
   implicit none
   private

   public :: local_charge, update_ncoord_list, get_cn_entries_list


   !> Threshold of the CN and local charge derivatives on the neighborlist below
   !> which their products in the partial derivatives and the Hessian are skipped
   real(wp), parameter :: cnthr = epsilon(1.0_wp)**2


contains


!> Local charges from the electronegativity weighted CN
subroutine local_charge(ncoord_en, mol, trans, qloc, dqlocdr, dqlocdL, &
   & list, dqlocdrlist)

   !> Electronegativity weighted CN, the total charge is distributed equally
   !> if not allocated
   class(ncoord_type), allocatable, intent(in) :: ncoord_en

   !> Molecular structure data
   type(structure_type), intent(in) :: mol

   !> Lattice translation vectors
   real(wp), intent(in) :: trans(:, :)

   !> Local atomic partial charges
   real(wp), intent(out) :: qloc(:)

   !> Optional derivative of local atomic partial charges w.r.t. atomic positions
   real(wp), intent(out), optional :: dqlocdr(3, mol%nat, mol%nat)

   !> Optional derivative of local atomic partial charges w.r.t. lattice vectors
   real(wp), intent(out), optional :: dqlocdL(3, 3, mol%nat)

   !> Multicharge neighborlist type
   type(csr_list), intent(in), optional :: list

   !> Optional derivative of local charges w.r.t. positions on the list
   real(wp), intent(out), optional :: dqlocdrlist(:, :)

   qloc = 0.0_wp
   if (present(dqlocdr) .and. present(dqlocdL)) then
      dqlocdr = 0.0_wp
      dqlocdL = 0.0_wp
   end if
   if (present(list) .and. present(dqlocdrlist) .and. present(dqlocdL)) then
      dqlocdrlist = 0.0_wp
      dqlocdL = 0.0_wp
   end if
   ! Get the electronegativity weighted CN for local charge
   if (allocated(ncoord_en)) then
      call ncoord_en%get_coordination_number(mol, trans, qloc, &
         & dcndr=dqlocdr, dcndrlist=dqlocdrlist, dcndL=dqlocdL, list=list)
   end if

   ! Distribute the total charge equally
   qloc = qloc + mol%charge / real(mol%nat, wp)

end subroutine local_charge


!> Coordination numbers and local charges on a sorted neighborlist
subroutine update_ncoord_list(ncoord, ncoord_en, mol, cache, trans, deriv, list)

   !> Coordination number
   class(ncoord_type), intent(in) :: ncoord

   !> Electronegativity weighted CN for local charge
   class(ncoord_type), allocatable, intent(in) :: ncoord_en

   !> Structure type
   type(structure_type), intent(in) :: mol

   !> Multicharge cache
   type(mchrg_cache), intent(inout) :: cache

   !> Lattice vectors
   real(wp), intent(in) :: trans(:, :)

   !> Flag for the position and strain derivatives
   logical, intent(in) :: deriv

   !> Sorted neighborlist, half or complete, holds the pair arrays of the cache
   type(csr_list), intent(in) :: list

   integer(i8) :: nnz

   ! Derivatives of a previous update may belong to another list, the transposed
   ! positions and CN entries on the complete list are formed on first use
   if (allocated(cache%dcndrlist)) deallocate(cache%dcndrlist)
   if (allocated(cache%dqlocdrlist)) deallocate(cache%dqlocdrlist)
   if (allocated(cache%cnmap)) deallocate(cache%cnmap)
   if (allocated(cache%tpos)) deallocate(cache%tpos)
   if (allocated(cache%cninl)) deallocate(cache%cninl)
   if (allocated(cache%cnpos)) deallocate(cache%cnpos)
   cache%hinl = list%inl

   if (deriv) then
      nnz = size(list%nlat, kind=i8)
      allocate(cache%dcndrlist(3, nnz), cache%dqlocdrlist(3, nnz))
      call ncoord%get_coordination_number(mol, trans, cache%cn, &
         & dcndL=cache%dcndL, list=list, dcndrlist=cache%dcndrlist)
      call local_charge(ncoord_en, mol, trans, cache%qloc, dqlocdL=cache%dqlocdL, &
         & list=list, dqlocdrlist=cache%dqlocdrlist)
      call compress_cn_derivs(cache%dcndrlist, cache%dqlocdrlist, cache%cnmap)
   else
      call ncoord%get_coordination_number(mol, trans, cache%cn, list=list)
      call local_charge(ncoord_en, mol, trans, cache%qloc, list=list)
   end if

end subroutine update_ncoord_list


!> Compress the CN and local charge derivatives on a list to the entries above
!> cnthr, the derivatives decay much faster than the pair distances of the list
subroutine compress_cn_derivs(dcndrlist, dqlocdrlist, cnmap)

   !> CN derivatives on the list, replaced by the compressed columns
   real(wp), allocatable, intent(inout) :: dcndrlist(:, :)

   !> Local charge derivatives on the list, replaced by the compressed columns
   real(wp), allocatable, intent(inout) :: dqlocdrlist(:, :)

   !> Column of each entry of the list, 0 for the zero column
   integer, allocatable, intent(out) :: cnmap(:)

   integer(i8) :: kat, nnz, ncn
   real(wp), allocatable :: dcn(:, :), dqloc(:, :)

   nnz = size(dcndrlist, 2, kind=i8)
   allocate(cnmap(nnz))
   !$omp parallel do default(none) schedule(static) &
   !$omp shared(nnz, dcndrlist, dqlocdrlist, cnmap) private(kat)
   do kat = 1, nnz
      cnmap(kat) = merge(1, 0, any(abs(dcndrlist(:, kat)) > cnthr) &
         & .or. any(abs(dqlocdrlist(:, kat)) > cnthr))
   end do
   !$omp end parallel do

   ncn = 0
   do kat = 1, nnz
      if (cnmap(kat) > 0) then
         ncn = ncn + 1
         if (ncn > huge(cnmap)) &
            & error stop "Too many CN derivative entries for the int32 cnmap"
         cnmap(kat) = int(ncn)
      end if
   end do

   allocate(dcn(3, 0:ncn), dqloc(3, 0:ncn))
   dcn(:, 0) = 0.0_wp
   dqloc(:, 0) = 0.0_wp
   !$omp parallel do default(none) schedule(static) &
   !$omp shared(nnz, dcndrlist, dqlocdrlist, cnmap, dcn, dqloc) private(kat)
   do kat = 1, nnz
      if (cnmap(kat) > 0) then
         dcn(:, cnmap(kat)) = dcndrlist(:, kat)
         dqloc(:, cnmap(kat)) = dqlocdrlist(:, kat)
      end if
   end do
   !$omp end parallel do
   call move_alloc(dcn, dcndrlist)
   call move_alloc(dqloc, dqlocdrlist)

end subroutine compress_cn_derivs


!> Entries of a sorted complete neighborlist with CN or local charge
subroutine get_cn_entries_list(mol, list, cache, tpos, cninl, cnpos)

   !> Molecular structure data
   type(structure_type), intent(in) :: mol

   !> Sorted complete neighborlist
   type(csr_list), intent(in) :: list

   !> Cache with the CN and local charge derivatives on the list of the update
   type(mchrg_cache), intent(in) :: cache

   !> Position of the atom in the row of the neighbor for each entry
   integer(i8), intent(in) :: tpos(:)

   !> Offsets of the rows of the entries
   integer(i8), allocatable, intent(out) :: cninl(:)

   !> Positions of the entries in the sorted complete list
   integer(i8), allocatable, intent(out) :: cnpos(:)

   integer :: iat
   integer(i8) :: kat, tat, hat, hoff, nfin, cnt
   logical :: half
   logical, allocatable :: keep(:)

   ! On a half list both entries of a pair share the position of the upper one
   half = cache%hinl(mol%nat + 1) /= list%inl(mol%nat + 1)

   allocate(cninl(mol%nat + 1))
   allocate(keep(size(list%nlat, kind=i8)))
   cninl(1) = 1

   !$omp parallel default(none) &
   !$omp shared(mol, list, cache, tpos, keep, cninl, half) &
   !$omp private(iat, kat, tat, hat, hoff, nfin, cnt)
   !$omp do schedule(runtime)
   do iat = 1, mol%nat
      nfin = find_lower_start(list, iat) - 1
      hoff = list%inl(iat) - cache%hinl(iat)
      do kat = list%inl(iat) + 1, nfin
         tat = tpos(kat)
         hat = kat - hoff
         keep(kat) = cache%cnmap(hat) > 0 .or. cache%cnmap(merge(hat, tat, half)) > 0
         keep(tat) = keep(kat)
      end do
   end do
   !$omp end do

   !$omp do schedule(runtime)
   do iat = 1, mol%nat
      cnt = 1
      do kat = list%inl(iat) + 1, list%inl(iat + 1) - 1
         if (keep(kat)) cnt = cnt + 1
      end do
      cninl(iat + 1) = cnt
   end do
   !$omp end do
   !$omp end parallel

   do iat = 1, mol%nat
      cninl(iat + 1) = cninl(iat) + cninl(iat + 1)
   end do

   allocate(cnpos(cninl(mol%nat + 1) - 1))
   !$omp parallel do default(none) schedule(runtime) &
   !$omp shared(mol, list, keep, cninl, cnpos) private(iat, kat, cnt)
   do iat = 1, mol%nat
      cnt = cninl(iat)
      cnpos(cnt) = list%inl(iat)
      do kat = list%inl(iat) + 1, list%inl(iat + 1) - 1
         if (keep(kat)) then
            cnt = cnt + 1
            cnpos(cnt) = kat
         end if
      end do
   end do
   !$omp end parallel do
   deallocate(keep)

end subroutine get_cn_entries_list


end module multicharge_ncoord
