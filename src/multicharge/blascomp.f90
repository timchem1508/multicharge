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
    use multicharge_adjlist, only : adjacency_list
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

        ! Diagonal only for symmetric matrices
        if (is_sym) then
            y(i) = y(i) + alpha * mdiag(i) * x(i)
        end if

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

        if (is_sym) then
            y(:, i) = y(:, i) + alpha * mdiag(:, i) * x(i)
        end if

        do k = list%inl(i) + 1, list%inl(i) + list%nnl(i)
            j = list%nlat(k)

            y(:, i) = y(:, i) + alpha * mlist(:, k) * x(j)

            if (is_sym) then
                y(:, j) = y(:, j) + alpha * mlist(:, k) * x(i)
            else
                y(:, j) = y(:, j) - alpha * mlist(:, k) * x(i)
            end if
        end do
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

end module multicharge_blascomp