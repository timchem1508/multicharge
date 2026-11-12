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

!> @file multicharge/solver/trialprod.f90
!> Provides the product of the linear system matrix with a trial vector

!> Trial products of the linear system matrix in its different representations.
!>
!> Iterative solvers only access the coefficient matrix through the product
!> with a trial vector and its diagonal, such that the dense matrix, the
!> neighbour list and the matrix-free operators share a single interface.
!> The dense and the sparse trial products reference the matrix of the caller,
!> which must outlive the trial product.
module multicharge_solver_trialprod
   use mctc_csrlist, only : csr_list, spgemv_csr, spsymv_csr
   use mctc_env, only : i8, wp
   use multicharge_blas, only : symv
   implicit none
   private

   public :: trialprod, dense_trialprod, csr_trialprod
   public :: new_dense_trialprod, new_csr_trialprod


   !> Abstract product of the linear system matrix with a trial vector
   type, abstract :: trialprod
   contains

      !> Multiply the linear system matrix with a trial vector
      procedure(multiply), deferred :: multiply

      !> Diagonal of the linear system matrix
      procedure(get_diagonal), deferred :: get_diagonal

      !> Dense block of the linear system matrix for a set of indices
      procedure(get_block), deferred :: get_block

   end type trialprod


   abstract interface
      !> Multiply the linear system matrix with a trial vector
      subroutine multiply(self, xvec, yvec)
         import :: trialprod, wp

         !> Trial product operator
         class(trialprod), intent(in) :: self

         !> Trial vector
         real(wp), intent(in) :: xvec(:)

         !> Product of the matrix with the trial vector
         real(wp), intent(out) :: yvec(:)

      end subroutine multiply

      !> Diagonal of the linear system matrix
      subroutine get_diagonal(self, diag)
         import :: trialprod, wp

         !> Trial product operator
         class(trialprod), intent(in) :: self

         !> Diagonal elements of the matrix
         real(wp), intent(out) :: diag(:)

      end subroutine get_diagonal

      !> Dense block of the linear system matrix for a set of indices, used
      !> to set up local preconditioners. Matrix-free operators may return
      !> an approximation, the block must remain symmetric positive definite.
      subroutine get_block(self, idx, block)
         import :: trialprod, wp

         !> Trial product operator
         class(trialprod), intent(in) :: self

         !> Row and column indices of the block
         integer, intent(in) :: idx(:)

         !> Matrix elements of the block
         real(wp), intent(out) :: block(:, :)

      end subroutine get_block
   end interface


   !> Trial product with a dense symmetric matrix
   type, extends(trialprod) :: dense_trialprod

      !> Dense symmetric matrix, referenced but not owned
      real(wp), pointer, contiguous :: amat(:, :) => null()
   contains

      !> Dense symmetric matrix-vector product
      procedure :: multiply => dense_multiply

      !> Diagonal of the dense matrix
      procedure :: get_diagonal => dense_get_diagonal

      !> Dense block of the matrix
      procedure :: get_block => dense_get_block

   end type dense_trialprod


   !> Trial product with a sparse matrix in the CSR neighbour list layout
   type, extends(trialprod) :: csr_trialprod

      !> Matrix values in the sparsity pattern of the list, referenced
      real(wp), pointer, contiguous :: val(:) => null()

      !> CSR neighbour list with the diagonal first in each row, referenced
      type(csr_list), pointer :: list => null()
   contains

      !> Sparse matrix-vector product
      procedure :: multiply => csr_multiply

      !> Diagonal of the sparse matrix
      procedure :: get_diagonal => csr_get_diagonal

      !> Dense block of the sparse matrix
      procedure :: get_block => csr_get_block

   end type csr_trialprod


contains


!> Create a trial product referencing a dense symmetric matrix
subroutine new_dense_trialprod(self, amat)

   !> Dense trial product operator
   type(dense_trialprod), intent(out) :: self

   !> Dense symmetric matrix, must outlive the trial product
   real(wp), intent(in), target, contiguous :: amat(:, :)

   self%amat => amat

end subroutine new_dense_trialprod


!> Create a trial product referencing a sparse matrix in the CSR layout
subroutine new_csr_trialprod(self, val, list)

   !> Sparse trial product operator
   type(csr_trialprod), intent(out) :: self

   !> Matrix values in the sparsity pattern of the list, must outlive the
   !> trial product
   real(wp), intent(in), target, contiguous :: val(:)

   !> CSR neighbour list, must outlive the trial product
   type(csr_list), intent(in), target :: list

   self%val => val
   self%list => list

end subroutine new_csr_trialprod


!> Dense symmetric matrix-vector product
subroutine dense_multiply(self, xvec, yvec)

   !> Dense trial product operator
   class(dense_trialprod), intent(in) :: self

   !> Trial vector
   real(wp), intent(in) :: xvec(:)

   !> Product of the matrix with the trial vector
   real(wp), intent(out) :: yvec(:)

   call symv(self%amat, xvec, yvec, alpha=1.0_wp, beta=0.0_wp)

end subroutine dense_multiply


!> Diagonal of the dense matrix
subroutine dense_get_diagonal(self, diag)

   !> Dense trial product operator
   class(dense_trialprod), intent(in) :: self

   !> Diagonal elements of the matrix
   real(wp), intent(out) :: diag(:)

   integer :: iat

   do iat = 1, size(diag)
      diag(iat) = self%amat(iat, iat)
   end do

end subroutine dense_get_diagonal


!> Dense block of the dense matrix
subroutine dense_get_block(self, idx, block)

   !> Dense trial product operator
   class(dense_trialprod), intent(in) :: self

   !> Row and column indices of the block
   integer, intent(in) :: idx(:)

   !> Matrix elements of the block
   real(wp), intent(out) :: block(:, :)

   integer :: ii, jj

   do jj = 1, size(idx)
      do ii = 1, size(idx)
         block(ii, jj) = self%amat(idx(ii), idx(jj))
      end do
   end do

end subroutine dense_get_block


!> Sparse matrix-vector product, the upper triangular list is symmetrized
subroutine csr_multiply(self, xvec, yvec)

   !> Sparse trial product operator
   class(csr_trialprod), intent(in) :: self

   !> Trial vector
   real(wp), intent(in) :: xvec(:)

   !> Product of the matrix with the trial vector
   real(wp), intent(out) :: yvec(:)

   integer :: ndim

   ndim = size(self%list%inl) - 1
   if (self%list%complete) then
      call spgemv_csr(ndim, self%val, self%list%inl, self%list%nlat, xvec, yvec)
   else
      call spsymv_csr(ndim, self%val, self%list%inl, self%list%nlat, xvec, yvec)
   end if

end subroutine csr_multiply


!> Diagonal of the sparse matrix, stored first in each row of the list
subroutine csr_get_diagonal(self, diag)

   !> Sparse trial product operator
   class(csr_trialprod), intent(in) :: self

   !> Diagonal elements of the matrix
   real(wp), intent(out) :: diag(:)

   integer :: iat

   do iat = 1, size(diag)
      diag(iat) = self%val(self%list%inl(iat))
   end do

end subroutine csr_get_diagonal


!> Dense block of the sparse matrix, elements outside of the list vanish
subroutine csr_get_block(self, idx, block)

   !> Sparse trial product operator
   class(csr_trialprod), intent(in) :: self

   !> Row and column indices of the block
   integer, intent(in) :: idx(:)

   !> Matrix elements of the block
   real(wp), intent(out) :: block(:, :)

   integer :: ii, jj, iat, jat
   integer(i8) :: ij

   do jj = 1, size(idx)
      do ii = 1, jj
         ! The upper triangular list stores each pair in the row of the
         ! smaller index, the complete list stores it in both rows
         iat = min(idx(ii), idx(jj))
         jat = max(idx(ii), idx(jj))
         block(ii, jj) = 0.0_wp
         do ij = self%list%inl(iat), self%list%inl(iat + 1) - 1
            if (self%list%nlat(ij) == jat) block(ii, jj) = block(ii, jj) + self%val(ij)
         end do
         block(jj, ii) = block(ii, jj)
      end do
   end do

end subroutine csr_get_block


end module multicharge_solver_trialprod
