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

!> @file multicharge/adjlist.f90
!> Sparse neighbour map / adjacency list implementation.

!> Implementation of a sparse neighbour map in compressed sparse row format.
!>
!> A symmetric neighbour map given in dense format like
!>
!>   |   | 1 | 2 | 3 | 4 | 5 | 6 |
!>   |---|---|---|---|---|---|---|
!>   | 1 |   | x |   | x | x |   |
!>   | 2 | x |   | x |   | x | x |
!>   | 3 |   | x |   | x |   | x |
!>   | 4 | x |   | x |   | x | x |
!>   | 5 | x | x |   | x |   |   |
!>   | 6 |   | x | x | x |   |   |
!>
!> Is stored in two compressed array identifying the neighbouring atom `nlat`
!> and its cell index `nltr`. Two index arrays `inl` for the offset
!> and `nnl` for the number of entries map the atomic index to the row index.
!>
!> ```
!> inl   =  0,       3,          7,      10,         14,      17, 20
!> nnl   =  |  2 ->  |  3 ->     |  2 ->  |  3 ->     |  2 ->  |  |
!> nlat  =     2, 4, 5, 1, 3, 5, 6, 2, 4, 6, 1, 3, 5, 6, 1, 2, 4, 2, 3, 4
!> nltr  =     1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1
!> ```
!>
!> An alternative representation would be to store just the offsets in `inl` and
!> additional beyond the last element the total number of neighbors. However,
!> the indexing is from inl(i) to inl(i+1)-1 could be confusing, therefore
!> two arrays are used for clarity.
module multicharge_adjlist
    use mctc_env, only : wp
    use mctc_io, only : structure_type
    use mctc_io_resize, only : resize
    implicit none
    private

    public :: adjacency_list, mchrg_adjlist_input, new_adjacency_list
    public :: symv_sparse, gemv_sparse, gemv_cmp, gemm_sparse

    !> @class adjacency_list
    !> Neighbourlist in CSR format
    type :: adjacency_list
        !> Realspace cutoff for neighbourlist generation
        real(wp), allocatable :: cutoff
        !> Whether a complete or a symmetrical reduced map should be generated
        logical, allocatable :: complete
        !> Offset index in the neighbour map
        integer, allocatable :: inl(:)
        !> Number of neighbours for each atom
        integer, allocatable :: nnl(:)
        !> Index of the neighbouring atom
        integer, allocatable :: nlat(:)
        !> Cell index of the neighbouring atom
        integer, allocatable :: nltr(:)
    end type adjacency_list

    !> Solver input abstract type
    type, public :: mchrg_adjlist_input
        !> Realspace cutoff for neighbourlist generation
        real(wp), allocatable :: cutoff
        !> Whether a complete or a symmetrical reduced map should be generated
        logical, allocatable :: complete
    end type mchrg_adjlist_input

    ! Default input 
    real(wp), parameter :: cutoff_def = 28.0_wp
    logical, parameter :: complete_def = .false.
    integer, parameter :: init_size = 10
    real(wp), parameter :: buffer = 0.1_wp

    real(wp), parameter :: eps = tiny(1.0_wp)

contains

    !> Create new neighbourlist for a given geometry and cutoff
    subroutine new_adjacency_list(self, input, mol, trans)
        !> Instance of the neighbourlist
        type(adjacency_list), intent(out) :: self
        !> Neighbour list input
        type(mchrg_adjlist_input), intent(in) :: input
        !> Molecular structure data
        type(structure_type), intent(in) :: mol
        !> Translation vectors for all images
        real(wp), intent(in) :: trans(:, :)

        if (allocated(input%cutoff)) then
            self%cutoff = input%cutoff
        else 
            self%cutoff = cutoff_def
        end if
        if (allocated(input%complete)) then
            self%complete = input%complete
        else 
            self%complete = complete_def
        end if

        allocate(self%inl(mol%nat), source=0)
        allocate(self%nnl(mol%nat), source=0)
        call generate(self, mol, trans)
    end subroutine new_adjacency_list

    !> Generator for neighbourlist using a Linked Cell List approach (O(N) scaling)
    subroutine generate(self, mol, trans)
        !> Instance of the neighbourlist
        type(adjacency_list), intent(inout) :: self
        !> Molecular structure data
        type(structure_type), intent(in) :: mol
        !> Translation vectors for all images
        real(wp), intent(in) :: trans(:, :)

        integer :: iat, jat, itr, img, ic, jc
        integer :: ix, iy, iz, jx, jy, jz, di, dj, dk
        integer, allocatable :: head(:), nxt(:)
        integer :: n_xyz(3)
        real(wp) :: r2, vec(3), cutoff2, cell_w(3), min_xyz(3), max_xyz(3)

        img = 0
        cutoff2 = self%cutoff**2

        ! 1. Define the grid boundaries and dimensions
        ! We add a small buffer to the bounding box to ensure all atoms are contained
        min_xyz = minval(mol%xyz, dim=2) - buffer
        max_xyz = maxval(mol%xyz, dim=2) + buffer
        
        ! Number of cells: must be at least 1, and cell width >= cutoff
        n_xyz = max(1, floor((max_xyz - min_xyz) / (self%cutoff + eps)))
        cell_w = (max_xyz - min_xyz) / (real(n_xyz, wp) + eps) + eps

        ! 2. Build the Linked List
        allocate(head(product(n_xyz)), source=0)
        allocate(nxt(mol%nat), source=0)
        
        do iat = 1, mol%nat
            ix = min(n_xyz(1), max(1, int((mol%xyz(1, iat) - min_xyz(1)) / cell_w(1)) + 1))
            iy = min(n_xyz(2), max(1, int((mol%xyz(2, iat) - min_xyz(2)) / cell_w(2)) + 1))
            iz = min(n_xyz(3), max(1, int((mol%xyz(3, iat) - min_xyz(3)) / cell_w(3)) + 1))
            
            ic = ix + n_xyz(1)*(iy-1) + n_xyz(1)*n_xyz(2)*(iz-1)
            nxt(iat) = head(ic)
            head(ic) = iat
        end do

        ! Pre-allocate neighbor arrays
        call resize(self%nlat, init_size*mol%nat)
        call resize(self%nltr, init_size*mol%nat)

        ! 3. Triple loop search over nearby cells (O(N) time)
        do iat = 1, mol%nat
            self%inl(iat) = img
            
            ix = min(n_xyz(1), max(1, int((mol%xyz(1, iat) - min_xyz(1)) / cell_w(1)) + 1))
            iy = min(n_xyz(2), max(1, int((mol%xyz(2, iat) - min_xyz(2)) / cell_w(2)) + 1))
            iz = min(n_xyz(3), max(1, int((mol%xyz(3, iat) - min_xyz(3)) / cell_w(3)) + 1))

            ! Check 27 neighboring cells (3x3x3 block)
            do dk = -1, 1; do dj = -1, 1; do di = -1, 1
                jx = ix + di; jy = iy + dj; jz = iz + dk
                
                ! Skip cells outside the defined grid
                if (jx < 1 .or. jx > n_xyz(1) .or. &
                    jy < 1 .or. jy > n_xyz(2) .or. &
                    jz < 1 .or. jz > n_xyz(3)) cycle
                
                jc = jx + n_xyz(1)*(jy-1) + n_xyz(1)*n_xyz(2)*(jz-1)
                jat = head(jc)
                
                do while (jat > 0)

                ! Symmetrical optimization: skip if jat > iat and complete is false
                if (.not. self%complete .and. jat > iat) then
                    jat = nxt(jat)
                    cycle
                end if

                ! Check all translation images for this atom pair
                do itr = 1, size(trans, 2)
                    vec(:) = mol%xyz(:, iat) - mol%xyz(:, jat) - trans(:, itr)
                    r2 = sum(vec**2)
                    
                    ! Standard distance check and self-interaction exclusion
                    if (r2 < epsilon(cutoff2) .or. r2 > cutoff2) cycle
                    
                    img = img + 1
                    if (size(self%nlat) < img) call resize(self%nlat)
                    if (size(self%nltr) < img) call resize(self%nltr)
                    self%nlat(img) = jat
                    self%nltr(img) = itr
                end do
                jat = nxt(jat)
                end do
            end do; end do; end do
            self%nnl(iat) = img - self%inl(iat)
        end do

        ! Cleanup and final sizing
        if (allocated(head)) deallocate(head)
        if (allocated(nxt)) deallocate(nxt)
        call resize(self%nlat, img)
        call resize(self%nltr, img)

    end subroutine generate

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