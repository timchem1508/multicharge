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

module multicharge_adjlist
    use mctc_env, only : wp
    use mctc_io, only : structure_type
    use mctc_io_resize, only : resize
    use mctc_ncoord_adjlist_type, only: adjacency_list
    implicit none
    private
    
    public :: symv_sparse, gemv_sparse, gemv_cmp

contains

    subroutine symv_sparse(adjlist, amat, x, y, alpha, beta)
        type(adjacency_list), intent(in) :: adjlist
        real(wp), intent(in) :: amat(:, :)
        real(wp), intent(in) :: x(:)
        real(wp), intent(inout) :: y(:)
        real(wp), intent(in) :: alpha, beta

        integer :: i, j, k, nat

        nat = size(adjlist%inl)

        ! Scale y by beta
        if (beta == 0.0_wp) then
            y(:) = 0.0_wp
        else if (beta /= 1.0_wp) then
            y(:) = beta * y(:)
        end if

        ! Loop over all atoms (rows)
        do i = 1, nat
            ! Add diagonal contribution
            y(i) = y(i) + alpha * amat(i, i) * x(i)

            ! Loop over off‑diagonal neighbors (stored in the reduced map)
            do k = adjlist%inl(i) + 1, adjlist%inl(i) + adjlist%nnl(i)
                j = adjlist%nlat(k)

                ! Update y(i) using the off‑diagonal element
                y(i) = y(i) + alpha * amat(i, j)  * x(j)

                ! Because the matrix is symmetric, also update y(j) using the same element
                if (i /= j) then
                    y(j) = y(j) + alpha * amat(i, j) * x(i)
                end if
            end do
        end do
    end subroutine symv_sparse

    subroutine gemv_sparse(adjlist, amat, x, y, alpha, beta)
        type(adjacency_list), intent(in) :: adjlist
        real(wp), intent(in) :: amat(:, :)
        real(wp), intent(in) :: x(:)
        real(wp), intent(inout) :: y(:)
        real(wp), intent(in) :: alpha, beta

        integer :: i, j, k, nat

        nat = size(adjlist%inl)

        !--------------------------------------------------
        ! Scale y by beta   (BLAS behaviour)
        !--------------------------------------------------
        if (beta == 0.0_wp) then
            y(:) = 0.0_wp
        else if (beta /= 1.0_wp) then
            y(:) = beta * y(:)
        end if

        !--------------------------------------------------
        ! y = alpha * A * x + y
        !--------------------------------------------------
        do i = 1, nat

            ! Diagonal element A(i,i)
            y(i) = y(i) + alpha * amat(i,i) * x(i)

            ! Off-diagonal nonzeros of row i
            do k = adjlist%inl(i) + 1, adjlist%inl(i) + adjlist%nnl(i)
                j = adjlist%nlat(k)

                y(i) = y(i) + alpha * amat(i,j) * x(j)

            end do
        end do

    end subroutine gemv_sparse

    subroutine gemv_cmp(list, mlist, mdiag, x, y, alpha, beta)
        type(adjacency_list), intent(in) :: list
        real(wp), intent(in)  :: mlist(:)   ! same size as list%nlat
        real(wp), intent(in)  :: mdiag(:)
        real(wp), intent(in)  :: x(:)
        real(wp), intent(out) :: y(:)
        real(wp), intent(in) :: alpha, beta

        integer :: i, k, j

        if (size(mlist) /= size(list%nlat)) return

        ! Scale y by beta
        if (beta == 0.0_wp) then
            y(:) = 0.0_wp
        else if (beta /= 1.0_wp) then
            y(:) = beta * y(:)
        end if

        do i = 1, size(list%nnl)

            ! Diagonal element A(i,i)
            y(i) = y(i) + alpha * mdiag(i) * x(i)

            ! Off-diagonal nonzeros of row i
            do k = list%inl(i) + 1, list%inl(i) + list%nnl(i)
                j = list%nlat(k)

                y(i) = y(i) + alpha * mlist(k) * x(j)
                y(j) = y(j) + alpha * mlist(k) * x(i) 
            end do
        end do
    end subroutine gemv_cmp

    subroutine gemm_sparse(adjlist, amat, B, C, alpha, beta)
        use omp_lib
        type(adjacency_list), intent(in) :: adjlist
        real(wp), intent(in) :: amat(:, :)
        real(wp), intent(in) :: B(:, :)
        real(wp), intent(inout) :: C(:, :)
        real(wp), intent(in) :: alpha, beta

        integer :: i, j, k, col, nat, ncol
        real(wp) :: aij

        nat  = size(adjlist%inl)
        ncol = size(B,2)

        !--------------------------------------------------
        ! Scale C by beta (parallel friendly)
        !--------------------------------------------------
        if (beta == 0.0_wp) then
            !$omp parallel do collapse(2) schedule(static)
            do i = 1, nat
                do col = 1, ncol
                    C(i,col) = 0.0_wp
                end do
            end do
        else if (beta /= 1.0_wp) then
            !$omp parallel do collapse(2) schedule(static)
            do i = 1, nat
                do col = 1, ncol
                    C(i,col) = beta * C(i,col)
                end do
            end do
        end if

        !--------------------------------------------------
        ! C = alpha * A * B + C
        ! Parallel over rows of A/C
        !--------------------------------------------------
        !$omp parallel do default(none) &
        !$omp private(i,j,k,col,aij) &
        !$omp shared(adjlist,amat,B,C,alpha,nat,ncol) &
        !$omp schedule(static)
        do i = 1, nat

            ! ----- diagonal -----
            aij = alpha * amat(i,i)
            if (aij /= 0.0_wp) then
                do col = 1, ncol
                    C(i,col) = C(i,col) + aij * B(i,col)
                end do
            end if

            ! ----- neighbour list (CSR row) -----
            do k = adjlist%inl(i) + 1, adjlist%inl(i) + adjlist%nnl(i)
                j   = adjlist%nlat(k)
                aij = alpha * amat(i,j)

                do col = 1, ncol
                    C(i,col) = C(i,col) + aij * B(j,col)
                end do
            end do

        end do
        !$omp end parallel do

    end subroutine gemm_sparse

end module multicharge_adjlist