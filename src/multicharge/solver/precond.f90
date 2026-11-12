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

!> @file multicharge/solver/precond.f90
!> Provides preconditioners for the iterative linear system solvers

!> Preconditioners of the linear system for the conjugate gradient solver.
!>
!> The factorized sparse approximate inverse approximates the inverse of the
!> linear system matrix by M^-1 = G^T G with a sparse triangular factor G.
!> Each row of G is obtained independently from a small dense system in the
!> neighbourhood of an atom, which keeps the setup local and parallel, and the
!> application only requires two sparse matrix-vector products.
!>
!> L. Yu. Kolotilina and A. Yu. Yeremin, *SIAM J. Matrix Anal. Appl.*,
!> **1993**, 14, 45.
!> DOI: [10.1137/0614004](https://doi.org/10.1137/0614004)
module multicharge_solver_precond
   use mctc_csrlist, only : csr_list
   use mctc_env, only : error_type, fatal_error, i8, wp
   use multicharge_solver_trialprod, only : trialprod
   implicit none
   private

   public :: precond_type, identity_precond, jacobi_precond, fsai_precond
   public :: new_jacobi_precond, new_fsai_precond, precond_kind


   !> Available preconditioners
   type :: enum_precond

      !> Preconditioner preferred by the charge model
      integer :: auto = -1

      !> No preconditioning
      integer :: none = 0

      !> Inverse diagonal of the matrix
      integer :: jacobi = 1

      !> Factorized sparse approximate inverse
      integer :: fsai = 2

   end type enum_precond

   !> Actual enumerator of the preconditioners
   type(enum_precond), parameter :: precond_kind = enum_precond()


   !> Abstract preconditioner, applies an approximate inverse of the matrix
   type, abstract :: precond_type
   contains

      !> Apply the preconditioner to a residual vector
      procedure(apply_precond), deferred :: apply

      !> Name of the preconditioner for printout
      procedure(get_name), deferred :: get_name

   end type precond_type


   abstract interface
      !> Apply the preconditioner to a residual vector
      subroutine apply_precond(self, rvec, zvec)
         import :: precond_type, wp

         !> Preconditioner instance
         class(precond_type), intent(in) :: self

         !> Residual vector
         real(wp), intent(in) :: rvec(:)

         !> Preconditioned residual vector
         real(wp), intent(out) :: zvec(:)

      end subroutine apply_precond

      !> Name of the preconditioner for printout
      pure function get_name(self) result(name)
         import :: precond_type

         !> Preconditioner instance
         class(precond_type), intent(in) :: self

         !> Name of the preconditioner
         character(len=:), allocatable :: name

      end function get_name
   end interface


   !> Identity, leaves the residual unchanged
   type, extends(precond_type) :: identity_precond
   contains

      !> Copy the residual vector
      procedure :: apply => identity_apply

      !> Name of the preconditioner
      procedure :: get_name => identity_get_name

   end type identity_precond


   !> Jacobi preconditioner using the inverse diagonal of the matrix
   type, extends(precond_type) :: jacobi_precond

      !> Inverse diagonal elements
      real(wp), allocatable :: dinv(:)
   contains

      !> Scale the residual by the inverse diagonal
      procedure :: apply => jacobi_apply

      !> Name of the preconditioner
      procedure :: get_name => jacobi_get_name

   end type jacobi_precond


   !> Factorized sparse approximate inverse M^-1 = G^T G
   type, extends(precond_type) :: fsai_precond

      !> Row offsets of the sparse factor G
      integer(i8), allocatable :: gptr(:)

      !> Column indices of the sparse factor G
      integer, allocatable :: gidx(:)

      !> Values of the sparse factor G
      real(wp), allocatable :: gval(:)

      !> Row offsets of the transposed factor G^T
      integer(i8), allocatable :: tptr(:)

      !> Column indices of the transposed factor G^T
      integer, allocatable :: tidx(:)

      !> Values of the transposed factor G^T
      real(wp), allocatable :: tval(:)

      !> Neighbour cutoff of the sparsity pattern in Bohr
      real(wp) :: cutoff = 0.0_wp
   contains

      !> Apply G^T G to the residual vector
      procedure :: apply => fsai_apply

      !> Name of the preconditioner
      procedure :: get_name => fsai_get_name

   end type fsai_precond

   !> Positive number used to prevent division by zero
   real(wp), parameter :: eps = tiny(1.0_wp)

   !> Minimal number of rows to set up the factor in parallel
   integer, parameter :: min_rows_setup = 256

   !> Minimal number of rows to apply the factor in parallel, below the
   !> overhead of a parallel region exceeds the sparse products
   integer, parameter :: min_rows_apply = 4096


contains


!> Create a Jacobi preconditioner from the diagonal of the matrix
subroutine new_jacobi_precond(self, tprod, ndim)

   !> Jacobi preconditioner instance
   type(jacobi_precond), intent(out) :: self

   !> Trial product operator representing the linear system matrix
   class(trialprod), intent(in) :: tprod

   !> Dimension of the linear system
   integer, intent(in) :: ndim

   allocate(self%dinv(ndim))
   call tprod%get_diagonal(self%dinv)
   self%dinv(:) = 1.0_wp / (self%dinv + eps)

end subroutine new_jacobi_precond


!> Create a factorized sparse approximate inverse.
!>
!> Row i of the factor G couples atom i with its neighbours j > i in the
!> pattern, its values solve A(J, J) g = e_i normalized to g_i = 1/sqrt(g_i).
!> Rows with a block that is not positive definite fall back to the inverse
!> square root of the diagonal.
subroutine new_fsai_precond(self, tprod, pattern, error)

   !> Factorized sparse approximate inverse instance
   type(fsai_precond), intent(out) :: self

   !> Trial product operator representing the linear system matrix
   class(trialprod), intent(in) :: tprod

   !> Neighbour list defining the sparsity pattern of the factor
   type(csr_list), intent(in) :: pattern

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   integer :: ndim, iat, nrow, info, kk
   integer(i8) :: ij
   integer, allocatable :: cnt(:), mark(:), idx(:)
   real(wp), allocatable :: block(:, :), gvec(:)

   if (.not. allocated(pattern%inl) .or. .not. allocated(pattern%nlat)) then
      call fatal_error(error, "Sparsity pattern of the FSAI preconditioner "//&
         & "is not initialized.")
      return
   end if
   ndim = size(pattern%inl) - 1
   if (allocated(pattern%cutoff)) self%cutoff = pattern%cutoff

   ! Count the unique neighbours of each row, the diagonal comes first
   allocate(cnt(ndim))
   !$omp parallel default(none) shared(pattern, ndim, cnt) &
   !$omp private(iat, mark, idx, nrow) if(ndim >= min_rows_setup)
   allocate(mark(ndim), source=0)
   !$omp do schedule(dynamic, 64)
   do iat = 1, ndim
      call get_row_pattern(pattern, iat, mark, idx, nrow)
      cnt(iat) = nrow
   end do
   !$omp end do
   !$omp end parallel

   allocate(self%gptr(ndim + 1))
   self%gptr(1) = 1_i8
   do iat = 1, ndim
      self%gptr(iat + 1) = self%gptr(iat) + int(cnt(iat), i8)
   end do
   allocate(self%gidx(self%gptr(ndim + 1) - 1), self%gval(self%gptr(ndim + 1) - 1))

   ! Solve the local systems of each row independently
   !$omp parallel default(none) shared(self, tprod, pattern, ndim) &
   !$omp private(iat, mark, idx, nrow, block, gvec, info, kk, ij) &
   !$omp if(ndim >= min_rows_setup)
   allocate(mark(ndim), source=0)
   !$omp do schedule(dynamic, 16)
   do iat = 1, ndim
      call get_row_pattern(pattern, iat, mark, idx, nrow)
      if (allocated(block)) then
         if (size(block, 1) < nrow) deallocate(block, gvec)
      end if
      if (.not. allocated(block)) allocate(block(nrow, nrow), gvec(nrow))
      call tprod%get_block(idx(:nrow), block(:nrow, :nrow))
      gvec(:nrow) = 0.0_wp
      gvec(1) = 1.0_wp
      call solve_spd(block(:nrow, :nrow), gvec(:nrow), info)
      if (info /= 0 .or. .not. gvec(1) > 0.0_wp) then
         call tprod%get_block(idx(:1), block(:1, :1))
         gvec(:nrow) = 0.0_wp
         gvec(1) = 1.0_wp / max(block(1, 1), eps)
      end if
      gvec(:nrow) = gvec(:nrow) / sqrt(gvec(1))
      do kk = 1, nrow
         ij = self%gptr(iat) + int(kk - 1, i8)
         self%gidx(ij) = idx(kk)
         self%gval(ij) = gvec(kk)
      end do
   end do
   !$omp end do
   !$omp end parallel

   call transpose_factor(self, ndim)

end subroutine new_fsai_precond


!> Unique neighbours of a row in the pattern, starting with the diagonal and
!> restricted to the upper triangle
subroutine get_row_pattern(pattern, iat, mark, idx, nrow)

   !> Neighbour list defining the sparsity pattern
   type(csr_list), intent(in) :: pattern

   !> Row index
   integer, intent(in) :: iat

   !> Work array marking visited columns with the current row, size ndim
   integer, intent(inout) :: mark(:)

   !> Column indices of the row
   integer, allocatable, intent(inout) :: idx(:)

   !> Number of unique columns in the row
   integer, intent(out) :: nrow

   integer :: jat
   integer(i8) :: ij, nmax

   nmax = pattern%inl(iat + 1) - pattern%inl(iat) + 1_i8
   if (allocated(idx)) then
      if (size(idx, kind=i8) < nmax) deallocate(idx)
   end if
   if (.not. allocated(idx)) allocate(idx(nmax))

   nrow = 1
   idx(1) = iat
   mark(iat) = iat
   do ij = pattern%inl(iat), pattern%inl(iat + 1) - 1
      jat = pattern%nlat(ij)
      if (jat <= iat .or. mark(jat) == iat) cycle
      mark(jat) = iat
      nrow = nrow + 1
      idx(nrow) = jat
   end do

end subroutine get_row_pattern


!> Solve a small symmetric positive definite system by Cholesky factorization
pure subroutine solve_spd(amat, bvec, info)

   !> Matrix, overwritten by its Cholesky factor
   real(wp), intent(inout) :: amat(:, :)

   !> Right-hand side, overwritten by the solution
   real(wp), intent(inout) :: bvec(:)

   !> Zero on success, otherwise the first non-positive pivot
   integer, intent(out) :: info

   integer :: ii, jj, nn
   real(wp) :: piv

   nn = size(bvec)
   info = 0
   do jj = 1, nn
      piv = amat(jj, jj) - dot_product(amat(jj, :jj - 1), amat(jj, :jj - 1))
      if (.not. piv > 0.0_wp) then
         info = jj
         return
      end if
      amat(jj, jj) = sqrt(piv)
      do ii = jj + 1, nn
         amat(ii, jj) = (amat(ii, jj) - dot_product(amat(ii, :jj - 1), &
            & amat(jj, :jj - 1))) / amat(jj, jj)
      end do
   end do

   ! Forward substitution L y = b
   do ii = 1, nn
      bvec(ii) = (bvec(ii) - dot_product(amat(ii, :ii - 1), bvec(:ii - 1))) &
         & / amat(ii, ii)
   end do
   ! Backward substitution L^T x = y
   do ii = nn, 1, -1
      bvec(ii) = (bvec(ii) - dot_product(amat(ii + 1:, ii), bvec(ii + 1:))) &
         & / amat(ii, ii)
   end do

end subroutine solve_spd


!> Store the transposed factor for a gather-only application
subroutine transpose_factor(self, ndim)

   !> Factorized sparse approximate inverse instance
   type(fsai_precond), intent(inout) :: self

   !> Dimension of the linear system
   integer, intent(in) :: ndim

   integer :: iat, jat
   integer(i8) :: ij, pos
   integer(i8), allocatable :: fill(:)

   allocate(self%tptr(ndim + 1), source=0_i8)
   do ij = 1, size(self%gidx, kind=i8)
      jat = self%gidx(ij)
      self%tptr(jat + 1) = self%tptr(jat + 1) + 1_i8
   end do
   self%tptr(1) = 1_i8
   do iat = 1, ndim
      self%tptr(iat + 1) = self%tptr(iat + 1) + self%tptr(iat)
   end do

   allocate(self%tidx(size(self%gidx)), self%tval(size(self%gval)))
   allocate(fill(ndim))
   fill(:) = self%tptr(:ndim)
   do iat = 1, ndim
      do ij = self%gptr(iat), self%gptr(iat + 1) - 1
         jat = self%gidx(ij)
         pos = fill(jat)
         self%tidx(pos) = iat
         self%tval(pos) = self%gval(ij)
         fill(jat) = pos + 1_i8
      end do
   end do

end subroutine transpose_factor


!> Copy the residual vector
subroutine identity_apply(self, rvec, zvec)

   !> Identity preconditioner instance
   class(identity_precond), intent(in) :: self

   !> Residual vector
   real(wp), intent(in) :: rvec(:)

   !> Preconditioned residual vector
   real(wp), intent(out) :: zvec(:)

   zvec(:) = rvec

end subroutine identity_apply


!> Scale the residual by the inverse diagonal
subroutine jacobi_apply(self, rvec, zvec)

   !> Jacobi preconditioner instance
   class(jacobi_precond), intent(in) :: self

   !> Residual vector
   real(wp), intent(in) :: rvec(:)

   !> Preconditioned residual vector
   real(wp), intent(out) :: zvec(:)

   zvec(:) = self%dinv * rvec

end subroutine jacobi_apply


!> Apply G^T G to the residual vector
subroutine fsai_apply(self, rvec, zvec)

   !> Factorized sparse approximate inverse instance
   class(fsai_precond), intent(in) :: self

   !> Residual vector
   real(wp), intent(in) :: rvec(:)

   !> Preconditioned residual vector
   real(wp), intent(out) :: zvec(:)

   integer :: iat, ndim
   integer(i8) :: ij
   real(wp) :: tmp
   real(wp), allocatable :: gres(:)

   ndim = size(self%gptr) - 1
   allocate(gres(ndim))

   !$omp parallel default(none) shared(self, rvec, zvec, gres, ndim) &
   !$omp private(iat, ij, tmp) if(ndim >= min_rows_apply)
   !$omp do schedule(static)
   do iat = 1, ndim
      tmp = 0.0_wp
      do ij = self%gptr(iat), self%gptr(iat + 1) - 1
         tmp = tmp + self%gval(ij) * rvec(self%gidx(ij))
      end do
      gres(iat) = tmp
   end do
   !$omp end do
   !$omp do schedule(static)
   do iat = 1, ndim
      tmp = 0.0_wp
      do ij = self%tptr(iat), self%tptr(iat + 1) - 1
         tmp = tmp + self%tval(ij) * gres(self%tidx(ij))
      end do
      zvec(iat) = tmp
   end do
   !$omp end do
   !$omp end parallel

end subroutine fsai_apply


!> Name of the identity preconditioner
pure function identity_get_name(self) result(name)

   !> Identity preconditioner instance
   class(identity_precond), intent(in) :: self

   !> Name of the preconditioner
   character(len=:), allocatable :: name

   name = "None"

end function identity_get_name


!> Name of the Jacobi preconditioner
pure function jacobi_get_name(self) result(name)

   !> Jacobi preconditioner instance
   class(jacobi_precond), intent(in) :: self

   !> Name of the preconditioner
   character(len=:), allocatable :: name

   name = "Jacobi (Diagonal)"

end function jacobi_get_name


!> Name of the factorized sparse approximate inverse with its cutoff
pure function fsai_get_name(self) result(name)

   !> Factorized sparse approximate inverse instance
   class(fsai_precond), intent(in) :: self

   !> Name of the preconditioner
   character(len=:), allocatable :: name

   character(len=32) :: buffer

   write(buffer, "(f0.2)") self%cutoff
   name = "FSAI (cutoff " // trim(buffer) // " Bohr)"

end function fsai_get_name


end module multicharge_solver_precond
