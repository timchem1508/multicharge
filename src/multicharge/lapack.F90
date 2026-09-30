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

#ifndef IK
#define IK i4
#endif

!> Interface to LAPACK library for symmetric indefinite factorization, solve and
!> inversion
module multicharge_lapack
   use mctc_env, only : sp, dp, ik => IK
   implicit none
   private

   public :: sytrf, sytrs, sytri


   !> Bunch-Kaufman factorization of a symmetric matrix
   interface sytrf
      module procedure :: mchrg_ssytrf
      module procedure :: mchrg_dsytrf
   end interface sytrf

   !> Solve a linear system using the sytrf factorization
   interface sytrs
      module procedure :: mchrg_ssytrs
      module procedure :: mchrg_ssytrs1
      module procedure :: mchrg_ssytrs3
      module procedure :: mchrg_dsytrs
      module procedure :: mchrg_dsytrs1
      module procedure :: mchrg_dsytrs3
   end interface sytrs

   !> Invert a symmetric matrix using the sytrf factorization
   interface sytri
      module procedure :: mchrg_ssytri
      module procedure :: mchrg_dsytri
   end interface sytri


   !> Symmetric indefinite factorization (LAPACK)
   interface lapack_sytrf
      pure subroutine ssytrf(uplo, n, a, lda, ipiv, work, lwork, info)
         import :: sp, ik
         integer(ik), intent(in) :: lda
         real(sp), intent(inout) :: a(lda, *)
         character(len=1), intent(in) :: uplo
         integer(ik), intent(out) :: ipiv(*)
         integer(ik), intent(out) :: info
         integer(ik), intent(in) :: n
         real(sp), intent(inout) :: work(*)
         integer(ik), intent(in) :: lwork
      end subroutine ssytrf
      pure subroutine dsytrf(uplo, n, a, lda, ipiv, work, lwork, info)
         import :: dp, ik
         integer(ik), intent(in) :: lda
         real(dp), intent(inout) :: a(lda, *)
         character(len=1), intent(in) :: uplo
         integer(ik), intent(out) :: ipiv(*)
         integer(ik), intent(out) :: info
         integer(ik), intent(in) :: n
         real(dp), intent(inout) :: work(*)
         integer(ik), intent(in) :: lwork
      end subroutine dsytrf
   end interface lapack_sytrf

   !> Solve with a symmetric indefinite factorization (LAPACK)
   interface lapack_sytrs
      pure subroutine ssytrs(uplo, n, nrhs, a, lda, ipiv, b, ldb, info)
         import :: sp, ik
         integer(ik), intent(in) :: lda
         integer(ik), intent(in) :: ldb
         real(sp), intent(in) :: a(lda, *)
         real(sp), intent(inout) :: b(ldb, *)
         integer(ik), intent(in) :: ipiv(*)
         character(len=1), intent(in) :: uplo
         integer(ik), intent(out) :: info
         integer(ik), intent(in) :: n
         integer(ik), intent(in) :: nrhs
      end subroutine ssytrs
      pure subroutine dsytrs(uplo, n, nrhs, a, lda, ipiv, b, ldb, info)
         import :: dp, ik
         integer(ik), intent(in) :: lda
         integer(ik), intent(in) :: ldb
         real(dp), intent(in) :: a(lda, *)
         real(dp), intent(inout) :: b(ldb, *)
         integer(ik), intent(in) :: ipiv(*)
         character(len=1), intent(in) :: uplo
         integer(ik), intent(out) :: info
         integer(ik), intent(in) :: n
         integer(ik), intent(in) :: nrhs
      end subroutine dsytrs
   end interface lapack_sytrs

   !> Symmetric indefinite inversion (LAPACK)
   interface lapack_sytri
      pure subroutine ssytri(uplo, n, a, lda, ipiv, work, info)
         import :: sp, ik
         integer(ik), intent(in) :: lda
         real(sp), intent(inout) :: a(lda, *)
         integer(ik), intent(in) :: ipiv(*)
         character(len=1), intent(in) :: uplo
         integer(ik), intent(out) :: info
         integer(ik), intent(in) :: n
         real(sp), intent(in) :: work(*)
      end subroutine ssytri
      pure subroutine dsytri(uplo, n, a, lda, ipiv, work, info)
         import :: dp, ik
         integer(ik), intent(in) :: lda
         real(dp), intent(inout) :: a(lda, *)
         integer(ik), intent(in) :: ipiv(*)
         character(len=1), intent(in) :: uplo
         integer(ik), intent(out) :: info
         integer(ik), intent(in) :: n
         real(dp), intent(in) :: work(*)
      end subroutine dsytri
   end interface lapack_sytri


contains


!> Symmetric indefinite factorization (single)
subroutine mchrg_ssytrf(amat, ipiv, uplo, info)

   !> Matrix A
   real(sp), intent(inout) :: amat(:, :)

   !> Pivot indices
   integer(ik), intent(out) :: ipiv(:)

   !> Optional triangle of A to reference ('u' or 'l')
   character(len=1), intent(in), optional :: uplo

   !> Optional status flag, absent aborts on failure
   integer(ik), intent(out), optional :: info

   character(len=1) :: ula
   integer(ik) :: stat, n, lda, lwork, stat_alloc, stat_dealloc
   real(sp), allocatable :: work(:)
   real(sp) :: test(1)

   if (present(uplo)) then
      ula = uplo
   else
      ula = 'u'
   end if
   lda = max(1, size(amat, 1))
   n = size(amat, 2)
   stat_alloc = 0_ik
   lwork = -1_ik
   call lapack_sytrf(ula, n, amat, lda, ipiv, test, lwork, stat)
   if (stat == 0) then
      lwork = nint(test(1))
      if (stat_alloc == 0) then
         allocate(work(lwork), stat=stat_alloc)
      end if
      if (stat_alloc == 0) then
         call lapack_sytrf(ula, n, amat, lda, ipiv, work, lwork, stat)
      else
         stat = -1000_ik
      end if
      deallocate(work, stat=stat_dealloc)
   end if
   if (present(info)) then
      info = stat
   else
      if (stat /= 0) error stop "[multicharge_lapack] ssytrf failed"
   end if

end subroutine mchrg_ssytrf


!> Symmetric indefinite factorization (double)
subroutine mchrg_dsytrf(amat, ipiv, uplo, info)

   !> Matrix A
   real(dp), intent(inout) :: amat(:, :)

   !> Pivot indices
   integer(ik), intent(out) :: ipiv(:)

   !> Optional triangle of A to reference ('u' or 'l')
   character(len=1), intent(in), optional :: uplo

   !> Optional status flag, absent aborts on failure
   integer(ik), intent(out), optional :: info

   character(len=1) :: ula
   integer(ik) :: stat, n, lda, lwork, stat_alloc, stat_dealloc
   real(dp), allocatable :: work(:)
   real(dp) :: test(1)

   if (present(uplo)) then
      ula = uplo
   else
      ula = 'u'
   end if
   lda = max(1, size(amat, 1))
   n = size(amat, 2)
   stat_alloc = 0_ik
   lwork = -1_ik
   call lapack_sytrf(ula, n, amat, lda, ipiv, test, lwork, stat)
   if (stat == 0) then
      lwork = nint(test(1))
      if (stat_alloc == 0) then
         allocate(work(lwork), stat=stat_alloc)
      end if
      if (stat_alloc == 0) then
         call lapack_sytrf(ula, n, amat, lda, ipiv, work, lwork, stat)
      else
         stat = -1000_ik
      end if
      deallocate(work, stat=stat_dealloc)
   end if
   if (present(info)) then
      info = stat
   else
      if (stat /= 0) error stop "[multicharge_lapack] dsytrf failed"
   end if

end subroutine mchrg_dsytrf


!> Solve with a symmetric indefinite factorization (single)
subroutine mchrg_ssytrs(amat, bmat, ipiv, uplo, info)

   !> Matrix A
   real(sp), intent(in) :: amat(:, :)

   !> Matrix B
   real(sp), intent(inout) :: bmat(:, :)

   !> Pivot indices
   integer(ik), intent(in) :: ipiv(:)

   !> Optional triangle of A to reference ('u' or 'l')
   character(len=1), intent(in), optional :: uplo

   !> Optional status flag, absent aborts on failure
   integer(ik), intent(out), optional :: info

   character(len=1) :: ula
   integer(ik) :: stat, n, nrhs, lda, ldb

   if (present(uplo)) then
      ula = uplo
   else
      ula = 'u'
   end if
   lda = max(1, size(amat, 1))
   ldb = max(1, size(bmat, 1))
   n = size(amat, 2)
   nrhs = size(bmat, 2)
   call lapack_sytrs(ula, n, nrhs, amat, lda, ipiv, bmat, ldb, stat)
   if (present(info)) then
      info = stat
   else
      if (stat /= 0) error stop "[multicharge_lapack] ssytrs failed"
   end if

end subroutine mchrg_ssytrs


!> Solve with a symmetric indefinite factorization (double)
subroutine mchrg_dsytrs(amat, bmat, ipiv, uplo, info)

   !> Matrix A
   real(dp), intent(in) :: amat(:, :)

   !> Matrix B
   real(dp), intent(inout) :: bmat(:, :)

   !> Pivot indices
   integer(ik), intent(in) :: ipiv(:)

   !> Optional triangle of A to reference ('u' or 'l')
   character(len=1), intent(in), optional :: uplo

   !> Optional status flag, absent aborts on failure
   integer(ik), intent(out), optional :: info

   character(len=1) :: ula
   integer(ik) :: stat, n, nrhs, lda, ldb

   if (present(uplo)) then
      ula = uplo
   else
      ula = 'u'
   end if
   lda = max(1, size(amat, 1))
   ldb = max(1, size(bmat, 1))
   n = size(amat, 2)
   nrhs = size(bmat, 2)
   call lapack_sytrs(ula, n, nrhs, amat, lda, ipiv, bmat, ldb, stat)
   if (present(info)) then
      info = stat
   else
      if (stat /= 0) error stop "[multicharge_lapack] dsytrs failed"
   end if

end subroutine mchrg_dsytrs


!> Solve with a symmetric indefinite factorization (single, array ranks 1)
subroutine mchrg_ssytrs1(amat, bvec, ipiv, uplo, info)

   !> Matrix A
   real(sp), intent(in) :: amat(:, :)

   !> Right-hand side vector b
   real(sp), intent(inout), target :: bvec(:)

   !> Pivot indices
   integer(ik), intent(in) :: ipiv(:)

   !> Optional triangle of A to reference ('u' or 'l')
   character(len=1), intent(in), optional :: uplo

   !> Optional status flag, absent aborts on failure
   integer(ik), intent(out), optional :: info

   real(sp), pointer :: bptr(:, :)

   bptr(1:size(bvec), 1:1) => bvec
   call sytrs(amat, bptr, ipiv, uplo, info)

end subroutine mchrg_ssytrs1


!> Solve with a symmetric indefinite factorization (single, array ranks 3)
subroutine mchrg_ssytrs3(amat, bmat, ipiv, uplo, info)

   !> Matrix A
   real(sp), intent(in) :: amat(:, :)

   !> Matrix B
   real(sp), intent(inout), contiguous, target :: bmat(:, :, :)

   !> Pivot indices
   integer(ik), intent(in) :: ipiv(:)

   !> Optional triangle of A to reference ('u' or 'l')
   character(len=1), intent(in), optional :: uplo

   !> Optional status flag, absent aborts on failure
   integer(ik), intent(out), optional :: info

   real(sp), pointer :: bptr(:, :)

   bptr(1:size(bmat, 1), 1:size(bmat, 2)*size(bmat, 3)) => bmat
   call sytrs(amat, bptr, ipiv, uplo, info)

end subroutine mchrg_ssytrs3


!> Solve with a symmetric indefinite factorization (double, array ranks 1)
subroutine mchrg_dsytrs1(amat, bvec, ipiv, uplo, info)

   !> Matrix A
   real(dp), intent(in) :: amat(:, :)

   !> Right-hand side vector b
   real(dp), intent(inout), target :: bvec(:)

   !> Pivot indices
   integer(ik), intent(in) :: ipiv(:)

   !> Optional triangle of A to reference ('u' or 'l')
   character(len=1), intent(in), optional :: uplo

   !> Optional status flag, absent aborts on failure
   integer(ik), intent(out), optional :: info

   real(dp), pointer :: bptr(:, :)

   bptr(1:size(bvec), 1:1) => bvec
   call sytrs(amat, bptr, ipiv, uplo, info)

end subroutine mchrg_dsytrs1


!> Solve with a symmetric indefinite factorization (double, array ranks 3)
subroutine mchrg_dsytrs3(amat, bmat, ipiv, uplo, info)

   !> Matrix A
   real(dp), intent(in) :: amat(:, :)

   !> Matrix B
   real(dp), intent(inout), contiguous, target :: bmat(:, :, :)

   !> Pivot indices
   integer(ik), intent(in) :: ipiv(:)

   !> Optional triangle of A to reference ('u' or 'l')
   character(len=1), intent(in), optional :: uplo

   !> Optional status flag, absent aborts on failure
   integer(ik), intent(out), optional :: info

   real(dp), pointer :: bptr(:, :)

   bptr(1:size(bmat, 1), 1:size(bmat, 2)*size(bmat, 3)) => bmat
   call sytrs(amat, bptr, ipiv, uplo, info)

end subroutine mchrg_dsytrs3


!> Symmetric indefinite inversion from its factorization (single)
subroutine mchrg_ssytri(amat, ipiv, uplo, info)

   !> Matrix A
   real(sp), intent(inout) :: amat(:, :)

   !> Pivot indices
   integer(ik), intent(in) :: ipiv(:)

   !> Optional triangle of A to reference ('u' or 'l')
   character(len=1), intent(in), optional :: uplo

   !> Optional status flag, absent aborts on failure
   integer(ik), intent(out), optional :: info

   character(len=1) :: ula
   integer(ik) :: stat, n, lda, stat_alloc, stat_dealloc
   real(sp), allocatable :: work(:)

   if (present(uplo)) then
      ula = uplo
   else
      ula = 'u'
   end if
   lda = max(1, size(amat, 1))
   n = size(amat, 2)
   stat_alloc = 0_ik
   allocate(work(n), stat=stat_alloc)
   if (stat_alloc == 0) then
      call lapack_sytri(ula, n, amat, lda, ipiv, work, stat)
   else
      stat = -1000_ik
   end if
   deallocate(work, stat=stat_dealloc)
   if (present(info)) then
      info = stat
   else
      if (stat /= 0) error stop "[multicharge_lapack] ssytri failed"
   end if

end subroutine mchrg_ssytri


!> Symmetric indefinite inversion from its factorization (double)
subroutine mchrg_dsytri(amat, ipiv, uplo, info)

   !> Matrix A
   real(dp), intent(inout) :: amat(:, :)

   !> Pivot indices
   integer(ik), intent(in) :: ipiv(:)

   !> Optional triangle of A to reference ('u' or 'l')
   character(len=1), intent(in), optional :: uplo

   !> Optional status flag, absent aborts on failure
   integer(ik), intent(out), optional :: info

   character(len=1) :: ula
   integer(ik) :: stat, n, lda, stat_alloc, stat_dealloc
   real(dp), allocatable :: work(:)

   if (present(uplo)) then
      ula = uplo
   else
      ula = 'u'
   end if
   lda = max(1, size(amat, 1))
   n = size(amat, 2)
   stat_alloc = 0_ik
   allocate(work(n), stat=stat_alloc)
   if (stat_alloc == 0) then
      call lapack_sytri(ula, n, amat, lda, ipiv, work, stat)
   else
      stat = -1000_ik
   end if
   deallocate(work, stat=stat_dealloc)
   if (present(info)) then
      info = stat
   else
      if (stat /= 0) error stop "[multicharge_lapack] dsytri failed"
   end if

end subroutine mchrg_dsytri


end module multicharge_lapack
