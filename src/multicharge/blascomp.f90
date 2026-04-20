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

   public :: gemv_cmp, gemm_cmp, gemv_cmp_212, gemm_cmp_212

   interface gemv_cmp
      module procedure gemv_cmp_111
      module procedure gemv_cmp_212
   end interface gemv_cmp

   interface gemv_cmp_212
      module procedure gemv_cmp_212
   end interface gemv_cmp_212

   interface gemm_cmp
      module procedure gemm_cmp_122
      module procedure gemm_cmp_222
      module procedure gemm_cmp_133
      module procedure gemm_cmp_211
      module procedure gemm_cmp_211_dir
   end interface gemm_cmp

   interface gemm_cmp_211_dir
      module procedure gemm_cmp_211_dir
   end interface gemm_cmp_211_dir

   interface gemm_cmp_212
      module procedure gemm_cmp_212
   end interface gemm_cmp_212


contains

!=========================================================
! GEMV 111
!=========================================================
   subroutine gemv_cmp_111(list, mlist, mdiag, x, y, alpha, beta, symmetric)
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

      !$omp parallel do default(shared) private(i, k, j) reduction(+:y) schedule(static)
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
      !$omp end parallel do
   end subroutine gemv_cmp_111

!=========================================================
! GEMV 212
!=========================================================
   subroutine gemv_cmp_212(list, mdrij, mdrji, mdrdiag, x, y, alpha, beta)
      type(adjacency_list), intent(in) :: list
      real(wp), intent(in)  :: mdrij(:,:)
      real(wp), intent(in)  :: mdrji(:,:)
      real(wp), intent(in)  :: mdrdiag(:,:)
      real(wp), intent(in)  :: x(:)
      real(wp), intent(inout) :: y(:,:)
      real(wp), intent(in), optional :: alpha
      real(wp), intent(in), optional :: beta

      real(wp) :: a, b
      integer  :: i, j, k

      a = 1.0_wp
      if (present(alpha)) a = alpha
      b = 0.0_wp
      if (present(beta)) b = beta

      if (b == 0.0_wp) then
         y(:, :) = 0.0_wp
      else if (b /= 1.0_wp) then
         y(:, :) = b * y(:, :)
      end if

      !$omp parallel do default(shared) private(i, k, j) reduction(+:y) schedule(static)
      do i = 1, size(list%nnl)

         do k = list%inl(i) + 1, list%inl(i) + list%nnl(i)
            j = list%nlat(k)

            ! Forward edge A(:, i, j)
            y(:, i) = y(:, i) + a * mdrij(:, k) * x(j)

            ! Backward edge A(:, j, i)
            y(:, j) = y(:, j) + a * mdrji(:, k) * x(i)

         end do

         ! Diagonal contribution
         y(:, i) = y(:, i) + a * mdrdiag(:, i) * x(i)

      end do
      !$omp end parallel do

   end subroutine gemv_cmp_212

!=========================================================
! GEMV 212 DIRECTED
!
! Y = Beta * Y + Alpha * M * X
! Matches the 9-argument signature with explicit i->j
! (drij) and j->i (drji) directed edge dependencies.
!=========================================================
   subroutine gemv_cmp_212_dir(list, mlist_drij, mlist_drji, mdiag, x, y, alpha, beta, symmetric)
      type(adjacency_list), intent(in) :: list
      real(wp), intent(in)  :: mlist_drij(:,:)
      real(wp), intent(in)  :: mlist_drji(:,:)
      real(wp), intent(in)  :: mdiag(:,:)
      real(wp), intent(in)  :: x(:)
      real(wp), intent(inout) :: y(:,:)
      real(wp), intent(in)  :: alpha, beta
      logical, intent(in), optional :: symmetric

      integer :: i, k, j
      logical :: is_sym

      is_sym = .true.
      if (present(symmetric)) is_sym = symmetric

      if (size(mlist_drij, 2) /= size(list%nlat)) return
      if (size(mlist_drji, 2) /= size(list%nlat)) return

      if (beta == 0.0_wp) then
         y(:,:) = 0.0_wp
      else if (beta /= 1.0_wp) then
         y(:,:) = beta * y(:,:)
      end if

      !$omp parallel do default(shared) private(i, k, j) reduction(+:y) schedule(static)
      do i = 1, size(list%nnl)

         do k = list%inl(i) + 1, list%inl(i) + list%nnl(i)
            j = list%nlat(k)

            ! Contribution to node i
            y(:, i) = y(:, i) + alpha * mlist_drij(:, k) * x(j)

            ! Contribution to node j
            if (is_sym) then
               y(:, j) = y(:, j) + alpha * mlist_drji(:, k) * x(i)
            else
               y(:, j) = y(:, j) - alpha * mlist_drji(:, k) * x(i)
            end if
         end do

         if (is_sym) then
            y(:, i) = y(:, i) + alpha * mdiag(:, i) * x(i)
         end if
      end do
      !$omp end parallel do
   end subroutine gemv_cmp_212_dir

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
   pure subroutine gemm_cmp_211_dir(list, clist, cdiag, drij, drji, ddiag, &
      xrij, xrji, xdiag, alpha, beta)
      use, intrinsic :: iso_fortran_env, only: wp => real64
      type(adjacency_list), intent(in) :: list
      real(wp), intent(in)    :: clist(:)          ! B off‑diagonal (i<j)
      real(wp), intent(in)    :: cdiag(:)          ! B diagonal
      real(wp), intent(in)    :: drij(:, :)        ! A(:,i,j) for i<j
      real(wp), intent(in)    :: drji(:, :)        ! A(:,j,i) for i<j
      real(wp), intent(in)    :: ddiag(:, :)       ! A(:,i,i)
      real(wp), intent(inout) :: xrij(:, :)        ! C(:,i,j) for i<j
      real(wp), intent(inout) :: xrji(:, :)        ! C(:,j,i) for i<j
      real(wp), intent(inout) :: xdiag(:, :)       ! C(:,i,i)
      real(wp), intent(in)    :: alpha, beta

      integer :: i, j, k, idx_ij, idx_ik, idx_jk
      integer :: ncomp
      real(wp) :: A_ij(3), A_ji(3), B_ij, B_ik, B_jk, B_ii, B_jj

      ncomp = size(drij, 1)   ! number of components (here 3)

      ! ----- beta scaling -----
      if (beta == 0.0_wp) then
         xrij  = 0.0_wp
         xrji  = 0.0_wp
         xdiag = 0.0_wp
      else if (beta /= 1.0_wp) then
         xrij  = beta * xrij
         xrji  = beta * xrji
         xdiag = beta * xdiag
      end if

      ! ----- Off‑diagonal contributions from A(i,j) and A(j,i) -----
      do i = 1, size(list%nnl)
         B_ii = cdiag(i)
         do idx_ij = list%inl(i) + 1, list%inl(i) + list%nnl(i)
            j = list%nlat(idx_ij)

            A_ij = drij(:, idx_ij)
            A_ji = drji(:, idx_ij)
            B_ij = clist(idx_ij)
            B_jj = cdiag(j)

            ! 1) C(i,j) += alpha * A(i,j) * B(j,j)
            xrij(:, idx_ij) = xrij(:, idx_ij) + alpha * A_ij * B_jj
            !    C(j,i) += alpha * A(j,i) * B(i,i)
            xrji(:, idx_ij) = xrji(:, idx_ij) + alpha * A_ji * B_ii

            ! 2) Diagonal contributions from this pair
            !    C(i,i) += alpha * A(i,j) * B(j,i)  (B(j,i)=B_ij)
            xdiag(:, i) = xdiag(:, i) + alpha * A_ij * B_ij
            !    C(j,j) += alpha * A(j,i) * B(i,j)  (B(i,j)=B_ij)
            xdiag(:, j) = xdiag(:, j) + alpha * A_ji * B_ij

            ! 3) Coupling to other neighbours:
            !    a) For every neighbour k of j (k /= i):
            !       C(i,k) += alpha * A(i,j) * B(j,k)
            do idx_jk = list%inl(j) + 1, list%inl(j) + list%nnl(j)
               k = list%nlat(idx_jk)
               if (k == i) cycle
               B_jk = clist(idx_jk)
               if (i < k) then
                  idx_ik = find_index(list, i, k)   ! helper function (see below)
                  xrij(:, idx_ik) = xrij(:, idx_ik) + alpha * A_ij * B_jk
               else if (i > k) then
                  idx_ik = find_index(list, k, i)   ! stored as (k,i)
                  xrji(:, idx_ik) = xrji(:, idx_ik) + alpha * A_ij * B_jk
               end if
            end do

            !    b) For every neighbour k of i (k /= j):
            !       C(j,k) += alpha * A(j,i) * B(i,k)
            do idx_ik = list%inl(i) + 1, list%inl(i) + list%nnl(i)
               k = list%nlat(idx_ik)
               if (k == j) cycle
               B_ik = clist(idx_ik)
               if (j < k) then
                  idx_jk = find_index(list, j, k)
                  xrij(:, idx_jk) = xrij(:, idx_jk) + alpha * A_ji * B_ik
               else if (j > k) then
                  idx_jk = find_index(list, k, j)
                  xrji(:, idx_jk) = xrji(:, idx_jk) + alpha * A_ji * B_ik
               end if
            end do
         end do
      end do

      ! ----- Diagonal contributions from A(i,i) -----
      do i = 1, size(list%nnl)
         A_ij = ddiag(:, i)          ! actually A_ii
         B_ii = cdiag(i)

         ! C(i,i) += alpha * A(i,i) * B(i,i)
         xdiag(:, i) = xdiag(:, i) + alpha * A_ij * B_ii

         ! C(i,j) += alpha * A(i,i) * B(i,j)   for j > i
         do idx_ij = list%inl(i) + 1, list%inl(i) + list%nnl(i)
            B_ij = clist(idx_ij)
            xrij(:, idx_ij) = xrij(:, idx_ij) + alpha * A_ij * B_ij
            ! Note: C(j,i) does *not* get a contribution from A(i,i)
         end do
      end do

   contains
      ! Helper: find compressed index for pair (i,j) with i < j
      pure function find_index(list, i, j) result(idx)
         type(adjacency_list), intent(in) :: list
         integer, intent(in) :: i, j
         integer :: idx, k
         do idx = list%inl(i) + 1, list%inl(i) + list%nnl(i)
            if (list%nlat(idx) == j) return
         end do
         idx = 0   ! should never happen if (i,j) is a neighbour pair
      end function
   end subroutine gemm_cmp_211_dir

   pure subroutine gemm_cmp_212(list, clist, cdiag, dtmpdrij, dtmpdrji, dtmpdrdiag, &
   & dxdrij, dxdrji, dxdrdiag, alpha, beta)
      !> Assumes wp (working precision) is accessible via module or host association.
      !> import :: wp
      type(adjacency_list), intent(in) :: list
      real(wp), intent(in) :: clist(:), cdiag(:)
      real(wp), intent(in) :: dtmpdrij(:,:), dtmpdrji(:,:), dtmpdrdiag(:,:)
      real(wp), intent(inout) :: dxdrij(:,:), dxdrji(:,:), dxdrdiag(:,:)
      real(wp), intent(in) :: alpha, beta

      integer :: nat, n_edges, i, j, k, e
      integer :: idx, idx_k, idx_y, y, e_k, e_y, m
      integer, allocatable :: deg(:), head(:)
      integer, allocatable :: adj_node(:), adj_edge(:)
      logical, allocatable :: adj_fw(:)
      real(wp) :: A_xk(3), B_ky
      real(wp), allocatable :: C_row(:,:)
      logical, allocatable :: modified(:)
      integer, allocatable :: mod_list(:)
      integer :: num_mod

      nat = size(list%nnl)
      n_edges = size(list%nlat)

      ! 1. Apply beta scaling upfront
      if (beta == 0.0_wp) then
         dxdrdiag = 0.0_wp
         dxdrij   = 0.0_wp
         dxdrji   = 0.0_wp
      else if (beta /= 1.0_wp) then
         dxdrdiag = beta * dxdrdiag
         dxdrij   = beta * dxdrij
         dxdrji   = beta * dxdrji
      end if

      ! 2. Build full symmetric neighbor list for O(1) reverse lookups
      allocate(deg(nat))
      deg = list%nnl
      do e = 1, n_edges
         j = list%nlat(e)
         deg(j) = deg(j) + 1
      end do

      allocate(head(nat + 1))
      head(1) = 1
      do i = 1, nat
         head(i+1) = head(i) + deg(i)
      end do

      allocate(adj_node(2 * n_edges))
      allocate(adj_edge(2 * n_edges))
      allocate(adj_fw(2 * n_edges))

      deg = head(1:nat) ! Reuse deg array as insertion pointers
      do i = 1, nat
         do e = list%inl(i) + 1, list%inl(i) + list%nnl(i)
            j = list%nlat(e)

            ! Forward edge (i -> j)
            idx = deg(i)
            adj_node(idx) = j
            adj_edge(idx) = e
            adj_fw(idx) = .true.
            deg(i) = deg(i) + 1

            ! Backward edge (j -> i)
            idx = deg(j)
            adj_node(idx) = i
            adj_edge(idx) = e
            adj_fw(idx) = .false.
            deg(j) = deg(j) + 1
         end do
      end do
      deallocate(deg)

      ! 3. Accumulate matrix product C = A * B row-by-row
      allocate(C_row(3, nat))
      C_row = 0.0_wp
      allocate(modified(nat))
      modified = .false.
      allocate(mod_list(nat))
      num_mod = 0

      do i = 1, nat

         ! -- Process k = i --
         A_xk(:) = dtmpdrdiag(:, i)

         ! y = i
         B_ky = cdiag(i)
         C_row(:, i) = C_row(:, i) + A_xk(:) * B_ky
         if (.not. modified(i)) then
            num_mod = num_mod + 1
            mod_list(num_mod) = i
            modified(i) = .true.
         end if

         ! y in N(i)
         do idx_y = head(i), head(i+1) - 1
            y = adj_node(idx_y)
            e_y = adj_edge(idx_y)
            B_ky = clist(e_y)
            C_row(:, y) = C_row(:, y) + A_xk(:) * B_ky
            if (.not. modified(y)) then
               num_mod = num_mod + 1
               mod_list(num_mod) = y
               modified(y) = .true.
            end if
         end do

         ! -- Process k in N(i) --
         do idx_k = head(i), head(i+1) - 1
            k = adj_node(idx_k)
            e_k = adj_edge(idx_k)

            ! Determine if we use the forward or backward asymmetric element
            if (adj_fw(idx_k)) then
               A_xk(:) = dtmpdrij(:, e_k)
            else
               A_xk(:) = dtmpdrji(:, e_k)
            end if

            ! y = k
            B_ky = cdiag(k)
            C_row(:, k) = C_row(:, k) + A_xk(:) * B_ky
            if (.not. modified(k)) then
               num_mod = num_mod + 1
               mod_list(num_mod) = k
               modified(k) = .true.
            end if

            ! y in N(k)
            do idx_y = head(k), head(k+1) - 1
               y = adj_node(idx_y)
               e_y = adj_edge(idx_y)
               B_ky = clist(e_y)
               C_row(:, y) = C_row(:, y) + A_xk(:) * B_ky
               if (.not. modified(y)) then
                  num_mod = num_mod + 1
                  mod_list(num_mod) = y
                  modified(y) = .true.
               end if
            end do
         end do

         ! 4. Scatter computed row back to compressed layout formats
         dxdrdiag(:, i) = dxdrdiag(:, i) + alpha * C_row(:, i)
         do idx_y = head(i), head(i+1) - 1
            y = adj_node(idx_y)
            e_y = adj_edge(idx_y)
            if (adj_fw(idx_y)) then
               dxdrij(:, e_y) = dxdrij(:, e_y) + alpha * C_row(:, y)
            else
               dxdrji(:, e_y) = dxdrji(:, e_y) + alpha * C_row(:, y)
            end if
         end do

         ! 5. Clear accumulator array sparse-ly
         do m = 1, num_mod
            y = mod_list(m)
            C_row(:, y) = 0.0_wp
            modified(y) = .false.
         end do
         num_mod = 0

      end do

   end subroutine gemm_cmp_212

end module multicharge_blascomp
