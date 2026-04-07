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

!> @file blascomp.f90
!> Matrix-vector and Matrix-matrix routines for CSR compressed matrices.
module multicharge_blascomp
   use mctc_env, only : wp
   use mctc_ncoord, only: adjacency_list
   implicit none
   private

   public :: gemv_cmp, gemm_cmp

   interface gemv_cmp
      module procedure gemv_cmp_111
      module procedure gemv_cmp_212
   end interface gemv_cmp

   interface gemm_cmp
      module procedure gemm_cmp_122
      module procedure gemm_cmp_222
      module procedure gemm_cmp_133
      module procedure gemm_cmp_211
      module procedure gemm_cmp_211_dir
   end interface gemm_cmp

contains

!=========================================================
! GEMV 111
!=========================================================
   pure subroutine gemv_cmp_111(list, mlist, mdiag, x, y, alpha, beta, symmetric)
      type(adjacency_list), intent(in) :: list
      real(wp), intent(in)  :: mlist(:)
      real(wp), intent(in)  :: mdiag(:)
      real(wp), intent(in)  :: x(:)
      real(wp), intent(inout) :: y(:)
      real(wp), intent(in)  :: alpha, beta
      logical, intent(in), optional :: symmetric

      integer :: i, k, j
      logical :: is_sym

      is_sym = .true.
      if (present(symmetric)) is_sym = symmetric

      if (size(mlist) /= size(list%nlat)) return

      if (beta == 0.0_wp) then
         y(:) = 0.0_wp
      else if (beta /= 1.0_wp) then
         y(:) = beta * y(:)
      end if

      do i = 1, size(list%nnl)

         do k = list%inl(i) + 1, list%inl(i) + list%nnl(i)
            j = list%nlat(k)

            ! A(i,j)
            y(i) = y(i) + alpha * mlist(k) * x(j)

            ! Mirror entry
            if (is_sym) then
               y(j) = y(j) + alpha * mlist(k) * x(i)
            else
               y(j) = y(j) - alpha * mlist(k) * x(i)
            end if
         end do

         ! Diagonal only for symmetric matrices
         if (is_sym) then
            y(i) = y(i) + alpha * mdiag(i) * x(i)
         end if

      end do
   end subroutine gemv_cmp_111


!=========================================================
! GEMV 212
!=========================================================
   pure subroutine gemv_cmp_212(list, mlist, mdiag, x, y, alpha, beta, symmetric)
      type(adjacency_list), intent(in) :: list
      real(wp), intent(in)  :: mlist(:,:)
      real(wp), intent(in)  :: mdiag(:,:)
      real(wp), intent(in)  :: x(:)
      real(wp), intent(inout) :: y(:,:)
      real(wp), intent(in)  :: alpha, beta
      logical, intent(in), optional :: symmetric

      integer :: i, k, j
      logical :: is_sym

      is_sym = .true.
      if (present(symmetric)) is_sym = symmetric

      if (size(mlist, 2) /= size(list%nlat)) return

      if (beta == 0.0_wp) then
         y(:,:) = 0.0_wp
      else if (beta /= 1.0_wp) then
         y(:,:) = beta * y(:,:)
      end if

      do i = 1, size(list%nnl)

         do k = list%inl(i) + 1, list%inl(i) + list%nnl(i)
            j = list%nlat(k)

            y(:, i) = y(:, i) + alpha * mlist(:, k) * x(j)

            if (is_sym) then
               y(:, j) = y(:, j) + alpha * mlist(:, k) * x(i)
            else
               y(:, j) = y(:, j) - alpha * mlist(:, k) * x(i)
            end if
         end do
         if (is_sym) then
            y(:, i) = y(:, i) + alpha * mdiag(:, i) * x(i)
         end if
      end do
   end subroutine gemv_cmp_212


!=========================================================
! GEMM 122
!=========================================================
   pure subroutine gemm_cmp_122(list, mlist, mdiag, x, y, alpha, beta, symmetric)
      type(adjacency_list), intent(in) :: list
      real(wp), intent(in)  :: mlist(:)
      real(wp), intent(in)  :: mdiag(:)
      real(wp), intent(in)  :: x(:,:)
      real(wp), intent(inout) :: y(:,:)
      real(wp), intent(in)  :: alpha, beta
      logical, intent(in), optional :: symmetric

      integer :: i, k, j
      logical :: is_sym

      is_sym = .true.
      if (present(symmetric)) is_sym = symmetric

      if (size(mlist) /= size(list%nlat)) return

      if (beta == 0.0_wp) then
         y(:,:) = 0.0_wp
      else if (beta /= 1.0_wp) then
         y(:,:) = beta * y(:,:)
      end if

      do i = 1, size(list%nnl)

         if (is_sym) then
            y(:, i) = y(:, i) + alpha * mdiag(i) * x(:, i)
         end if

         do k = list%inl(i) + 1, list%inl(i) + list%nnl(i)
            j = list%nlat(k)

            y(:, i) = y(:, i) + alpha * mlist(k) * x(:, j)

            if (is_sym) then
               y(:, j) = y(:, j) + alpha * mlist(k) * x(:, i)
            else
               y(:, j) = y(:, j) - alpha * mlist(k) * x(:, i)
            end if
         end do
      end do
   end subroutine gemm_cmp_122


!=========================================================
! GEMM 222
!=========================================================
   pure subroutine gemm_cmp_222(list, mlist, mdiag, x, y, alpha, beta, symmetric)
      type(adjacency_list), intent(in) :: list
      real(wp), intent(in)  :: mlist(:,:)
      real(wp), intent(in)  :: mdiag(:,:)
      real(wp), intent(in)  :: x(:,:)
      real(wp), intent(inout) :: y(:,:)
      real(wp), intent(in)  :: alpha, beta
      logical, intent(in), optional :: symmetric

      integer :: i, k, j
      logical :: is_sym

      is_sym = .true.
      if (present(symmetric)) is_sym = symmetric

      if (size(mlist, 2) /= size(list%nlat)) return

      if (beta == 0.0_wp) then
         y(:,:) = 0.0_wp
      else if (beta /= 1.0_wp) then
         y(:,:) = beta * y(:,:)
      end if

      do i = 1, size(list%nnl)

         if (is_sym) then
            y(:, i) = y(:, i) + alpha * mdiag(:, i) * x(:, i)
         end if

         do k = list%inl(i) + 1, list%inl(i) + list%nnl(i)
            j = list%nlat(k)

            y(:, i) = y(:, i) + alpha * mlist(:, k) * x(:, j)

            if (is_sym) then
               y(:, j) = y(:, j) + alpha * mlist(:, k) * x(:, i)
            else
               y(:, j) = y(:, j) - alpha * mlist(:, k) * x(:, i)
            end if
         end do
      end do
   end subroutine gemm_cmp_222

!=========================================================
! RIGHT GEMM 133
! Computes Y(:,:,i) = sum_j X(:,:,j) * C(j,i)
! where C is symmetric and stored via neighbour list
!=========================================================
   pure subroutine gemm_cmp_133(list, clist, cdiag, x, y, alpha, beta, symmetric)
      type(adjacency_list), intent(in) :: list
      real(wp), intent(in)  :: clist(:)
      real(wp), intent(in)  :: cdiag(:)
      real(wp), intent(in)  :: x(:,:,:)
      real(wp), intent(inout) :: y(:,:,:)
      real(wp), intent(in)  :: alpha, beta
      logical, intent(in), optional :: symmetric

      integer :: i, j, k
      logical :: is_sym
      is_sym = .true.
      if (present(symmetric)) is_sym = symmetric
      if (size(clist) /= size(list%nlat)) return

      ! Scale output
      if (beta == 0.0_wp) then
         y(:,:,:) = 0.0_wp
      else if (beta /= 1.0_wp) then
         y(:,:,:) = beta * y(:,:,:)
      end if

      do i = 1, size(list%nnl)

         ! Diagonal contribution
         y(:,:,i) = y(:,:,i) + alpha * cdiag(i) * x(:,:,i)

         ! Off-diagonal neighbours
         do k = list%inl(i) + 1, list%inl(i) + list%nnl(i)
            j = list%nlat(k)

            ! Y_i += X_j * C(j,i)
            y(:,:,i) = y(:,:,i) + alpha * clist(k) * x(:,:,j)

            if (is_sym) then
               ! symmetric partner
               y(:,:,j) = y(:,:,j) + alpha * clist(k) * x(:,:,i)
            else
               ! antisymmetric partner
               y(:,:,j) = y(:,:,j) - alpha * clist(k) * x(:,:,i)
            end if

         end do
      end do

   end subroutine gemm_cmp_133

!=========================================================
! GEMM COMPRESSED 211
!
! DX = DTMP * C
!
! DTMP may be symmetric OR antisymmetric
! C is symmetric
!=========================================================
   pure subroutine gemm_cmp_211(list, clist, cdiag, &
      dtmp_list, dtmp_diag, &
      dx_list, dx_diag, alpha, beta, &
      symmetric)

      type(adjacency_list), intent(in) :: list
      real(wp), intent(in)  :: clist(:)
      real(wp), intent(in)  :: cdiag(:)
      real(wp), intent(in)  :: dtmp_list(:,:)   ! (ncomp, nnz)
      real(wp), intent(in)  :: dtmp_diag(:,:)   ! (ncomp, nat)
      real(wp), intent(inout) :: dx_list(:,:)   ! (ncomp, nnz)
      real(wp), intent(inout) :: dx_diag(:,:)   ! (ncomp, nat)
      real(wp), intent(in) :: alpha, beta
      logical, intent(in), optional :: symmetric

      logical :: is_sym
      integer :: i, j, k

      is_sym = .true.
      if (present(symmetric)) is_sym = symmetric

      ! Scale outputs
      if (beta == 0.0_wp) then
         dx_list(:,:) = 0.0_wp
         dx_diag(:,:) = 0.0_wp
      else if (beta /= 1.0_wp) then
         dx_list(:,:) = beta * dx_list(:,:)
         dx_diag(:,:) = beta * dx_diag(:,:)
      end if

      do i = 1, size(list%nnl)

         !==================================================
         ! DIAGONAL
         !==================================================
         if (is_sym) then
            ! symmetric result
            dx_diag(:,i) = dx_diag(:,i) + alpha * dtmp_diag(:,i) * cdiag(i)

            do k = list%inl(i) + 1, list%inl(i) + list%nnl(i)
               j = list%nlat(k)
               dx_diag(:,i) = dx_diag(:,i) + alpha * dtmp_list(:,k) * clist(k)
               dx_diag(:,j) = dx_diag(:,j) + alpha * dtmp_list(:,k) * clist(k)
            end do
         else
            ! antisymmetric result → diagonal must remain ZERO
            dx_diag(:,i) = 0.0_wp
         end if

         !==================================================
         ! OFF-DIAGONAL
         !==================================================
         do k = list%inl(i) + 1, list%inl(i) + list%nnl(i)
            j = list%nlat(k)

            if (is_sym) then
               ! symmetric × symmetric → symmetric
               dx_list(:,k) = dx_list(:,k) + alpha * ( &
                  dtmp_diag(:,i) * clist(k) &
                  + cdiag(i)      * dtmp_list(:,k) &
                  + dtmp_list(:,k)* cdiag(j) &
                  + clist(k)      * dtmp_diag(:,j) )
            else
               ! antisymmetric × symmetric → antisymmetric
               dx_list(:,k) = dx_list(:,k) + alpha * &
                  dtmp_list(:,k) * (cdiag(j) - cdiag(i))
            end if

         end do

      end do

   end subroutine gemm_cmp_211

!=========================================================
! GEMM COMPRESSED 211 DIRECTED
!
! DX = DTMP * C
! Matches the 12-argument signature with explicit i->j
! (drij) and j->i (drji) directed edge dependencies.
!=========================================================
   pure subroutine gemm_cmp_211_dir(list, clist, cdiag, &
      dtmp_drij, dtmp_drji, dtmp_diag, &
      dx_drij, dx_drji, dx_diag, alpha, beta, &
      symmetric)

      type(adjacency_list), intent(in) :: list
      real(wp), intent(in)  :: clist(:)
      real(wp), intent(in)  :: cdiag(:)
      real(wp), intent(in)  :: dtmp_drij(:,:)   ! (ncomp, nnz)
      real(wp), intent(in)  :: dtmp_drji(:,:)   ! (ncomp, nnz)
      real(wp), intent(in)  :: dtmp_diag(:,:)   ! (ncomp, nat)
      real(wp), intent(inout) :: dx_drij(:,:)   ! (ncomp, nnz)
      real(wp), intent(inout) :: dx_drji(:,:)   ! (ncomp, nnz)
      real(wp), intent(inout) :: dx_diag(:,:)   ! (ncomp, nat)
      real(wp), intent(in) :: alpha, beta
      logical, intent(in), optional :: symmetric

      logical :: is_sym
      integer :: i, j, k

      is_sym = .true.
      if (present(symmetric)) is_sym = symmetric

      ! Scale outputs
      if (beta == 0.0_wp) then
         dx_drij(:,:) = 0.0_wp
         dx_drji(:,:) = 0.0_wp
         dx_diag(:,:) = 0.0_wp
      else if (beta /= 1.0_wp) then
         dx_drij(:,:) = beta * dx_drij(:,:)
         dx_drji(:,:) = beta * dx_drji(:,:)
         dx_diag(:,:) = beta * dx_diag(:,:)
      end if

      do i = 1, size(list%nnl)

         !==================================================
         ! DIAGONAL
         !==================================================
         if (is_sym) then
            ! symmetric result
            dx_diag(:,i) = dx_diag(:,i) + alpha * dtmp_diag(:,i) * cdiag(i)

            do k = list%inl(i) + 1, list%inl(i) + list%nnl(i)
               j = list%nlat(k)
               dx_diag(:,i) = dx_diag(:,i) + alpha * dtmp_drij(:,k) * clist(k)
               dx_diag(:,j) = dx_diag(:,j) + alpha * dtmp_drji(:,k) * clist(k)
            end do
         else
            ! antisymmetric result → diagonal must remain ZERO
            dx_diag(:,i) = 0.0_wp
         end if

         !==================================================
         ! OFF-DIAGONAL
         !==================================================
         do k = list%inl(i) + 1, list%inl(i) + list%nnl(i)
            j = list%nlat(k)

            if (is_sym) then
               ! symmetric × symmetric → symmetric (forward edge i->j)
               dx_drij(:,k) = dx_drij(:,k) + alpha * ( &
                  dtmp_diag(:,i) * clist(k) &
                  + cdiag(i)      * dtmp_drij(:,k) &
                  + dtmp_drij(:,k)* cdiag(j) &
                  + clist(k)      * dtmp_diag(:,j) )

               ! symmetric × symmetric → symmetric (backward edge j->i)
               dx_drji(:,k) = dx_drji(:,k) + alpha * ( &
                  dtmp_diag(:,j) * clist(k) &
                  + cdiag(j)      * dtmp_drji(:,k) &
                  + dtmp_drji(:,k)* cdiag(i) &
                  + clist(k)      * dtmp_diag(:,i) )
            else
               ! antisymmetric × symmetric → antisymmetric
               ! swapping nodes i and j shifts the target cdiag value for drji
               dx_drij(:,k) = dx_drij(:,k) + alpha * &
                  dtmp_drij(:,k) * (cdiag(j) - cdiag(i))

               dx_drji(:,k) = dx_drji(:,k) + alpha * &
                  dtmp_drji(:,k) * (cdiag(i) - cdiag(j))
            end if

         end do

      end do

   end subroutine gemm_cmp_211_dir

end module multicharge_blascomp
