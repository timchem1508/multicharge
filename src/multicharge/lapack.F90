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

module multicharge_lapack
   use mctc_env, only : sp, dp, ik => IK
   implicit none
   private

   public :: sytrf, sytrs, sytri, potrf, potrs, syevd

   interface sytrf
      module procedure :: mchrg_ssytrf
      module procedure :: mchrg_dsytrf
   end interface sytrf

   interface sytrs
      module procedure :: mchrg_ssytrs
      module procedure :: mchrg_ssytrs1
      module procedure :: mchrg_ssytrs3
      module procedure :: mchrg_dsytrs
      module procedure :: mchrg_dsytrs1
      module procedure :: mchrg_dsytrs3
   end interface sytrs

   interface sytri
      module procedure :: mchrg_ssytri
      module procedure :: mchrg_dsytri
   end interface sytri

   interface potrf
      module procedure :: mchrg_spotrf
      module procedure :: mchrg_dpotrf
   end interface potrf

   interface potrs
      module procedure :: mchrg_spotrs
      module procedure :: mchrg_dpotrs
   end interface potrs

   interface syevd
      module procedure :: mchrg_ssyevd
      module procedure :: mchrg_dsyevd
   end interface syevd


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

   interface lapack_potrf
      pure subroutine spotrf(uplo, n, a, lda, info)
         import :: sp, ik
         integer(ik), intent(in) :: lda
         real(sp), intent(inout) :: a(lda, *)
         character(len=1), intent(in) :: uplo
         integer(ik), intent(out) :: info
         integer(ik), intent(in) :: n
      end subroutine spotrf
      pure subroutine dpotrf(uplo, n, a, lda, info)
         import :: dp, ik
         integer(ik), intent(in) :: lda
         real(dp), intent(inout) :: a(lda, *)
         character(len=1), intent(in) :: uplo
         integer(ik), intent(out) :: info
         integer(ik), intent(in) :: n
      end subroutine dpotrf
   end interface lapack_potrf

   interface lapack_potrs
      pure subroutine spotrs(uplo, n, nrhs, a, lda, b, ldb, info)
         import :: sp, ik
         integer(ik), intent(in) :: lda
         integer(ik), intent(in) :: ldb
         real(sp), intent(in) :: a(lda, *)
         real(sp), intent(inout) :: b(ldb, *)
         character(len=1), intent(in) :: uplo
         integer(ik), intent(out) :: info
         integer(ik), intent(in) :: n
         integer(ik), intent(in) :: nrhs
      end subroutine spotrs
      pure subroutine dpotrs(uplo, n, nrhs, a, lda, b, ldb, info)
         import :: dp, ik
         integer(ik), intent(in) :: lda
         integer(ik), intent(in) :: ldb
         real(dp), intent(in) :: a(lda, *)
         real(dp), intent(inout) :: b(ldb, *)
         character(len=1), intent(in) :: uplo
         integer(ik), intent(out) :: info
         integer(ik), intent(in) :: n
         integer(ik), intent(in) :: nrhs
      end subroutine dpotrs
   end interface lapack_potrs

   interface lapack_syevd
      pure subroutine ssyevd(jobz, uplo, n, a, lda, w, work, lwork, iwork, &
            & liwork, info)
         import :: sp, ik
         character(len=1), intent(in) :: jobz
         character(len=1), intent(in) :: uplo
         integer(ik), intent(in) :: n
         integer(ik), intent(in) :: lda
         real(sp), intent(inout) :: a(lda, *)
         real(sp), intent(out) :: w(*)
         real(sp), intent(inout) :: work(*)
         integer(ik), intent(in) :: lwork
         integer(ik), intent(inout) :: iwork(*)
         integer(ik), intent(in) :: liwork
         integer(ik), intent(out) :: info
      end subroutine ssyevd
      pure subroutine dsyevd(jobz, uplo, n, a, lda, w, work, lwork, iwork, &
            & liwork, info)
         import :: dp, ik
         character(len=1), intent(in) :: jobz
         character(len=1), intent(in) :: uplo
         integer(ik), intent(in) :: n
         integer(ik), intent(in) :: lda
         real(dp), intent(inout) :: a(lda, *)
         real(dp), intent(out) :: w(*)
         real(dp), intent(inout) :: work(*)
         integer(ik), intent(in) :: lwork
         integer(ik), intent(inout) :: iwork(*)
         integer(ik), intent(in) :: liwork
         integer(ik), intent(out) :: info
      end subroutine dsyevd
   end interface lapack_syevd


contains


subroutine mchrg_ssytrf(amat, ipiv, uplo, info)
   real(sp), intent(inout) :: amat(:, :)
   integer(ik), intent(out) :: ipiv(:)
   character(len=1), intent(in), optional :: uplo
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
      if (stat_alloc==0) then
         allocate(work(lwork), stat=stat_alloc)
      end if
      if (stat_alloc==0) then
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


subroutine mchrg_dsytrf(amat, ipiv, uplo, info)
   real(dp), intent(inout) :: amat(:, :)
   integer(ik), intent(out) :: ipiv(:)
   character(len=1), intent(in), optional :: uplo
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
      if (stat_alloc==0) then
         allocate(work(lwork), stat=stat_alloc)
      end if
      if (stat_alloc==0) then
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


subroutine mchrg_ssytrs(amat, bmat, ipiv, uplo, info)
   real(sp), intent(in) :: amat(:, :)
   real(sp), intent(inout) :: bmat(:, :)
   integer(ik), intent(in) :: ipiv(:)
   character(len=1), intent(in), optional :: uplo
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


subroutine mchrg_dsytrs(amat, bmat, ipiv, uplo, info)
   real(dp), intent(in) :: amat(:, :)
   real(dp), intent(inout) :: bmat(:, :)
   integer(ik), intent(in) :: ipiv(:)
   character(len=1), intent(in), optional :: uplo
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


subroutine mchrg_ssytrs1(amat, bvec, ipiv, uplo, info)
   real(sp), intent(in) :: amat(:, :)
   real(sp), intent(inout), target :: bvec(:)
   integer(ik), intent(in) :: ipiv(:)
   character(len=1), intent(in), optional :: uplo
   integer(ik), intent(out), optional :: info
   real(sp), pointer :: bptr(:, :)
   bptr(1:size(bvec), 1:1) => bvec
   call sytrs(amat, bptr, ipiv, uplo, info)
end subroutine mchrg_ssytrs1


subroutine mchrg_ssytrs3(amat, bmat, ipiv, uplo, info)
   real(sp), intent(in) :: amat(:, :)
   real(sp), intent(inout), contiguous, target :: bmat(:, :, :)
   integer(ik), intent(in) :: ipiv(:)
   character(len=1), intent(in), optional :: uplo
   integer(ik), intent(out), optional :: info
   real(sp), pointer :: bptr(:, :)
   bptr(1:size(bmat, 1), 1:size(bmat, 2)*size(bmat, 3)) => bmat
   call sytrs(amat, bptr, ipiv, uplo, info)
end subroutine mchrg_ssytrs3


subroutine mchrg_dsytrs1(amat, bvec, ipiv, uplo, info)
   real(dp), intent(in) :: amat(:, :)
   real(dp), intent(inout), target :: bvec(:)
   integer(ik), intent(in) :: ipiv(:)
   character(len=1), intent(in), optional :: uplo
   integer(ik), intent(out), optional :: info
   real(dp), pointer :: bptr(:, :)
   bptr(1:size(bvec), 1:1) => bvec
   call sytrs(amat, bptr, ipiv, uplo, info)
end subroutine mchrg_dsytrs1


subroutine mchrg_dsytrs3(amat, bmat, ipiv, uplo, info)
   real(dp), intent(in) :: amat(:, :)
   real(dp), intent(inout), contiguous, target :: bmat(:, :, :)
   integer(ik), intent(in) :: ipiv(:)
   character(len=1), intent(in), optional :: uplo
   integer(ik), intent(out), optional :: info
   real(dp), pointer :: bptr(:, :)
   bptr(1:size(bmat, 1), 1:size(bmat, 2)*size(bmat, 3)) => bmat
   call sytrs(amat, bptr, ipiv, uplo, info)
end subroutine mchrg_dsytrs3


subroutine mchrg_ssytri(amat, ipiv, uplo, info)
   real(sp), intent(inout) :: amat(:, :)
   integer(ik), intent(in) :: ipiv(:)
   character(len=1), intent(in), optional :: uplo
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
   if (stat_alloc==0) then
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


subroutine mchrg_dsytri(amat, ipiv, uplo, info)
   real(dp), intent(inout) :: amat(:, :)
   integer(ik), intent(in) :: ipiv(:)
   character(len=1), intent(in), optional :: uplo
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
   if (stat_alloc==0) then
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


subroutine mchrg_spotrf(amat, uplo, info)
   real(sp), intent(inout) :: amat(:, :)
   character(len=1), intent(in), optional :: uplo
   integer(ik), intent(out), optional :: info
   character(len=1) :: ula
   integer(ik) :: stat, n, lda
   if (present(uplo)) then
      ula = uplo
   else
      ula = 'u'
   end if
   lda = max(1, size(amat, 1))
   n = size(amat, 2)
   call lapack_potrf(ula, n, amat, lda, stat)
   if (present(info)) then
      info = stat
   else
      if (stat /= 0) error stop "[multicharge_lapack] spotrf failed"
   end if
end subroutine mchrg_spotrf


subroutine mchrg_spotrs(amat, bmat, uplo, info)
   real(sp), intent(in) :: amat(:, :)
   real(sp), intent(inout) :: bmat(:, :)
   character(len=1), intent(in), optional :: uplo
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
   call lapack_potrs(ula, n, nrhs, amat, lda, bmat, ldb, stat)
   if (present(info)) then
      info = stat
   else
      if (stat /= 0) error stop "[multicharge_lapack] spotrs failed"
   end if
end subroutine mchrg_spotrs


subroutine mchrg_ssyevd(amat, eval, jobz, uplo, info)
   real(sp), intent(inout) :: amat(:, :)
   real(sp), intent(out) :: eval(:)
   character(len=1), intent(in), optional :: jobz
   character(len=1), intent(in), optional :: uplo
   integer(ik), intent(out), optional :: info
   character(len=1) :: job, ula
   integer(ik) :: stat, n, lda, lwork, liwork, stat_alloc, stat_dealloc
   real(sp), allocatable :: work(:)
   integer(ik), allocatable :: iwork(:)
   real(sp) :: test(1)
   integer(ik) :: itest(1)
   if (present(jobz)) then
      job = jobz
   else
      job = 'v'
   end if
   if (present(uplo)) then
      ula = uplo
   else
      ula = 'u'
   end if
   lda = max(1, size(amat, 1))
   n = size(amat, 2)
   stat_alloc = 0_ik
   lwork = -1_ik
   liwork = -1_ik
   call lapack_syevd(job, ula, n, amat, lda, eval, test, lwork, itest, liwork, &
      & stat)
   if (stat == 0) then
      lwork = nint(test(1))
      liwork = itest(1)
      allocate(work(lwork), iwork(liwork), stat=stat_alloc)
      if (stat_alloc==0) then
         call lapack_syevd(job, ula, n, amat, lda, eval, work, lwork, iwork, &
            & liwork, stat)
      else
         stat = -1000_ik
      end if
      deallocate(work, iwork, stat=stat_dealloc)
   end if
   if (present(info)) then
      info = stat
   else
      if (stat /= 0) error stop "[multicharge_lapack] ssyevd failed"
   end if
end subroutine mchrg_ssyevd


subroutine mchrg_dpotrf(amat, uplo, info)
   real(dp), intent(inout) :: amat(:, :)
   character(len=1), intent(in), optional :: uplo
   integer(ik), intent(out), optional :: info
   character(len=1) :: ula
   integer(ik) :: stat, n, lda
   if (present(uplo)) then
      ula = uplo
   else
      ula = 'u'
   end if
   lda = max(1, size(amat, 1))
   n = size(amat, 2)
   call lapack_potrf(ula, n, amat, lda, stat)
   if (present(info)) then
      info = stat
   else
      if (stat /= 0) error stop "[multicharge_lapack] dpotrf failed"
   end if
end subroutine mchrg_dpotrf


subroutine mchrg_dpotrs(amat, bmat, uplo, info)
   real(dp), intent(in) :: amat(:, :)
   real(dp), intent(inout) :: bmat(:, :)
   character(len=1), intent(in), optional :: uplo
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
   call lapack_potrs(ula, n, nrhs, amat, lda, bmat, ldb, stat)
   if (present(info)) then
      info = stat
   else
      if (stat /= 0) error stop "[multicharge_lapack] dpotrs failed"
   end if
end subroutine mchrg_dpotrs


subroutine mchrg_dsyevd(amat, eval, jobz, uplo, info)
   real(dp), intent(inout) :: amat(:, :)
   real(dp), intent(out) :: eval(:)
   character(len=1), intent(in), optional :: jobz
   character(len=1), intent(in), optional :: uplo
   integer(ik), intent(out), optional :: info
   character(len=1) :: job, ula
   integer(ik) :: stat, n, lda, lwork, liwork, stat_alloc, stat_dealloc
   real(dp), allocatable :: work(:)
   integer(ik), allocatable :: iwork(:)
   real(dp) :: test(1)
   integer(ik) :: itest(1)
   if (present(jobz)) then
      job = jobz
   else
      job = 'v'
   end if
   if (present(uplo)) then
      ula = uplo
   else
      ula = 'u'
   end if
   lda = max(1, size(amat, 1))
   n = size(amat, 2)
   stat_alloc = 0_ik
   lwork = -1_ik
   liwork = -1_ik
   call lapack_syevd(job, ula, n, amat, lda, eval, test, lwork, itest, liwork, &
      & stat)
   if (stat == 0) then
      lwork = nint(test(1))
      liwork = itest(1)
      allocate(work(lwork), iwork(liwork), stat=stat_alloc)
      if (stat_alloc==0) then
         call lapack_syevd(job, ula, n, amat, lda, eval, work, lwork, iwork, &
            & liwork, stat)
      else
         stat = -1000_ik
      end if
      deallocate(work, iwork, stat=stat_dealloc)
   end if
   if (present(info)) then
      info = stat
   else
      if (stat /= 0) error stop "[multicharge_lapack] dsyevd failed"
   end if
end subroutine mchrg_dsyevd


end module multicharge_lapack
