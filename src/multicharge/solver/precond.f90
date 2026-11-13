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
   use mctc_io, only : structure_type
   use mctc_io_constants, only : pi
   use mctc_io_math, only : matdet_3x3, matinv_3x3
   use multicharge_blas, only : symv
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

      !> Number of coarse cells in the two-level correction
      integer :: ncoarse = 0

      !> Coarse cell index of each atom
      integer, allocatable :: cidx(:)

      !> Coarse correction matrix W^-1 - B^-1
      real(wp), allocatable :: cmat(:, :)
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

   !> Maximum number of coarse cells per axis to bound coarse setup cost
   integer, parameter :: max_coarse_grid = 256

   !> Target maximum total number of coarse cells in 3D
   integer, parameter :: max_coarse_cells = 300

   real(wp), parameter :: sqrtpi = sqrt(pi)


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
subroutine new_fsai_precond(self, tprod, pattern, error, mol)

   !> Factorized sparse approximate inverse instance
   type(fsai_precond), intent(out) :: self

   !> Trial product operator representing the linear system matrix
   class(trialprod), intent(in) :: tprod

   !> Neighbour list defining the sparsity pattern of the factor
   type(csr_list), intent(in) :: pattern

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   !> Optional molecular structure for the two-level coarse correction
   type(structure_type), intent(in), optional :: mol

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
   deallocate(mark)
   if (allocated(idx)) deallocate(idx)
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
   deallocate(mark)
   if (allocated(idx)) deallocate(idx)
   if (allocated(block)) deallocate(block, gvec)
   !$omp end parallel

   call transpose_factor(self, ndim)
   if (present(mol) .and. self%cutoff > 0.0_wp) then
      call setup_coarse_correction(self, mol, ndim)
   end if

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


!> Construct the two-level coarse-space Woodbury correction W^-1 - B^-1
subroutine setup_coarse_correction(self, mol, ndim)

   !> Factorized sparse approximate inverse instance
   type(fsai_precond), intent(inout) :: self

   !> Molecular structure data
   type(structure_type), intent(in) :: mol

   !> Dimension of the linear system
   integer, intent(in) :: ndim

   integer :: iat, idir, ic, jc, nc, info
   integer :: ngrid(3), ix(3), m1, m2, m3, mmax(3)
   integer(i8) :: ij, kl
   integer, allocatable :: raw_cid(:), map(:)
   real(wp) :: hc, xmin(3), xmax(3), span(3), u(3), inv_lat(3, 3), rec_lat(3, 3)
   real(wp) :: vol, alpha_c, r, gvec(3), g2, wg, val, dxyz(3), scale_grid
   real(wp), allocatable :: gmat(:, :), wmat(:, :), bmat(:, :), wvec(:)
   real(wp), allocatable :: rcell(:, :), ref_xyz(:, :), cphase(:), sphase(:)
   integer, allocatable :: cnt_cell(:)

   if (ndim < 16 .or. mol%nat /= ndim) return
   if (any(mol%periodic) .and. .not. all(mol%periodic)) return

   hc = 0.8_wp * self%cutoff
   if (all(mol%periodic)) then
      inv_lat = matinv_3x3(mol%lattice)
      do idir = 1, 3
         span(idir) = norm2(mol%lattice(:, idir))
         ngrid(idir) = max(1, nint(span(idir) / hc))
      end do
   else
      xmin = minval(mol%xyz, dim=2)
      xmax = maxval(mol%xyz, dim=2)
      span = xmax - xmin
      do idir = 1, 3
         ngrid(idir) = max(1, ceiling(span(idir) / hc))
      end do
   end if

   if (product(ngrid) > max_coarse_cells) then
      scale_grid = (real(max_coarse_cells, wp) / real(product(ngrid), wp))**(1.0_wp / 3.0_wp)
      do idir = 1, 3
         ngrid(idir) = max(1, min(max_coarse_grid, nint(real(ngrid(idir), wp) * scale_grid)))
      end do
   else
      do idir = 1, 3
         ngrid(idir) = min(max_coarse_grid, ngrid(idir))
      end do
   end if

   if (product(ngrid) <= 1) return

   allocate(raw_cid(ndim))
   do iat = 1, ndim
      if (all(mol%periodic)) then
         u = matmul(inv_lat, mol%xyz(:, iat))
         u = u - floor(u)
         do idir = 1, 3
            ix(idir) = modulo(int(u(idir) * real(ngrid(idir), wp)), ngrid(idir))
         end do
      else
         do idir = 1, 3
            if (span(idir) > eps) then
               u(idir) = (mol%xyz(idir, iat) - xmin(idir)) / span(idir)
               ix(idir) = min(ngrid(idir) - 1, max(0, int(u(idir) * real(ngrid(idir), wp))))
            else
               ix(idir) = 0
            end if
         end do
      end if
      raw_cid(iat) = 1 + ix(1) + ngrid(1) * (ix(2) + ngrid(2) * ix(3))
   end do

   allocate(map(product(ngrid)), source=0)
   nc = 0
   do iat = 1, ndim
      ic = raw_cid(iat)
      if (map(ic) == 0) then
         nc = nc + 1
         map(ic) = nc
      end if
   end do

   if (nc <= 1) return

   self%ncoarse = nc
   allocate(self%cidx(ndim))
   do iat = 1, ndim
      self%cidx(iat) = map(raw_cid(iat))
   end do

   ! Compute coarse cell centers and W = (G P_0)^T (G P_0)
   allocate(rcell(3, nc), ref_xyz(3, nc), source=0.0_wp)
   allocate(cnt_cell(nc), source=0)
   do iat = 1, ndim
      ic = self%cidx(iat)
      if (cnt_cell(ic) == 0) then
         ref_xyz(:, ic) = mol%xyz(:, iat)
         rcell(:, ic) = mol%xyz(:, iat)
      else if (all(mol%periodic)) then
         dxyz = mol%xyz(:, iat) - ref_xyz(:, ic)
         u = matmul(inv_lat, dxyz)
         dxyz = matmul(mol%lattice, u - real(nint(u), wp))
         rcell(:, ic) = rcell(:, ic) + ref_xyz(:, ic) + dxyz
      else
         rcell(:, ic) = rcell(:, ic) + mol%xyz(:, iat)
      end if
      cnt_cell(ic) = cnt_cell(ic) + 1
   end do
   do ic = 1, nc
      rcell(:, ic) = rcell(:, ic) / real(cnt_cell(ic), wp)
   end do

   allocate(wmat(nc, nc), source=0.0_wp)
   do iat = 1, ndim
      do ij = self%gptr(iat), self%gptr(iat + 1) - 1
         ic = self%cidx(self%gidx(ij))
         val = self%gval(ij)
         wmat(ic, ic) = wmat(ic, ic) + val * val
         do kl = ij + 1, self%gptr(iat + 1) - 1
            jc = self%cidx(self%gidx(kl))
            if (ic == jc) then
               wmat(ic, ic) = wmat(ic, ic) + 2.0_wp * val * self%gval(kl)
            else
               wmat(ic, jc) = wmat(ic, jc) + val * self%gval(kl)
               wmat(jc, ic) = wmat(jc, ic) + val * self%gval(kl)
            end if
         end do
      end do
   end do

   allocate(wvec(nc))
   do ic = 1, nc
      wvec(ic) = sum(wmat(:, ic))
   end do

   ! Build analytical coarse matrix B = W + D_w K_long D_w
   allocate(bmat(nc, nc))
   bmat(:, :) = wmat
   alpha_c = 1.2_wp / self%cutoff

   if (all(mol%periodic)) then
      vol = abs(matdet_3x3(mol%lattice))
      rec_lat = 2.0_wp * pi * transpose(inv_lat)
      do idir = 1, 3
         mmax(idir) = max(2, min(6, ceiling(2.5_wp * alpha_c * span(idir) / pi)))
      end do
      allocate(cphase(nc), sphase(nc))
      do m3 = -mmax(3), mmax(3)
         do m2 = -mmax(2), mmax(2)
            do m1 = -mmax(1), mmax(1)
               if (m1 == 0 .and. m2 == 0 .and. m3 == 0) cycle
               gvec = matmul(rec_lat, real([m1, m2, m3], wp))
               g2 = dot_product(gvec, gvec)
               if (0.25_wp * g2 / (alpha_c * alpha_c) > 10.0_wp) cycle
               wg = (4.0_wp * pi / vol) * exp(-0.25_wp * g2 / (alpha_c * alpha_c)) / g2
               do ic = 1, nc
                  val = dot_product(gvec, rcell(:, ic))
                  cphase(ic) = wvec(ic) * cos(val)
                  sphase(ic) = wvec(ic) * sin(val)
               end do
               do jc = 1, nc
                  bmat(:, jc) = bmat(:, jc) + wg * (cphase(jc) * cphase(:) + sphase(jc) * sphase(:))
               end do
            end do
         end do
      end do
   else
      do jc = 1, nc
         bmat(jc, jc) = bmat(jc, jc) + wvec(jc) * wvec(jc) * (2.0_wp * alpha_c / sqrtpi)
         do ic = 1, jc - 1
            r = norm2(rcell(:, ic) - rcell(:, jc))
            if (r > 1.0e-10_wp) then
               val = wvec(ic) * wvec(jc) * erf(alpha_c * r) / r
            else
               val = wvec(ic) * wvec(jc) * (2.0_wp * alpha_c / sqrtpi)
            end if
            bmat(ic, jc) = bmat(ic, jc) + val
            bmat(jc, ic) = bmat(ic, jc)
         end do
      end do
   end if

   ! Invert W and B by Cholesky factorization: cmat = W^-1 - B^-1
   allocate(self%cmat(nc, nc), gmat(nc, nc))
   call invert_spd(wmat, self%cmat, info)
   if (info /= 0) then
      self%ncoarse = 0
      deallocate(self%cidx, self%cmat)
      return
   end if
   call invert_spd(bmat, gmat, info)
   if (info /= 0) then
      self%ncoarse = 0
      deallocate(self%cidx, self%cmat)
      return
   end if
   self%cmat(:, :) = self%cmat - gmat

end subroutine setup_coarse_correction


!> Invert a symmetric positive definite matrix by Cholesky factorization
pure subroutine invert_spd(amat, ainv, info)

   !> Matrix, overwritten by its Cholesky factor
   real(wp), intent(inout) :: amat(:, :)

   !> Inverse matrix
   real(wp), intent(out) :: ainv(:, :)

   !> Zero on success, otherwise the first non-positive pivot
   integer, intent(out) :: info

   integer :: ii, jj, nn
   real(wp) :: piv

   nn = size(amat, 1)
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

   ainv(:, :) = 0.0_wp
   do jj = 1, nn
      ainv(jj, jj) = 1.0_wp
      do ii = jj, nn
         ainv(ii, jj) = (ainv(ii, jj) - dot_product(amat(ii, jj:ii - 1), &
            & ainv(jj:ii - 1, jj))) / amat(ii, ii)
      end do
      do ii = nn, 1, -1
         ainv(ii, jj) = (ainv(ii, jj) - dot_product(amat(ii + 1:, ii), &
            & ainv(ii + 1:, jj))) / amat(ii, ii)
      end do
      do ii = 1, jj - 1
         ainv(ii, jj) = ainv(jj, ii)
      end do
   end do

end subroutine invert_spd


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


!> Apply G^T G and the optional two-level coarse correction to the residual vector
subroutine fsai_apply(self, rvec, zvec)

   !> Factorized sparse approximate inverse instance
   class(fsai_precond), intent(in) :: self

   !> Residual vector
   real(wp), intent(in) :: rvec(:)

   !> Preconditioned residual vector
   real(wp), intent(out) :: zvec(:)

   integer :: ndim
   real(wp), allocatable :: gres(:), cvec(:), dvec(:), pvec(:), corr(:)

   ndim = size(self%gptr) - 1
   allocate(gres(ndim))

   call fsai_raw_apply(self, ndim, rvec, gres, zvec)

   if (self%ncoarse > 1) then
      allocate(cvec(self%ncoarse), dvec(self%ncoarse), pvec(ndim), corr(ndim))
      call restrict_coarse(self%cidx, ndim, self%ncoarse, zvec, cvec)
      call symv(self%cmat, cvec, dvec)
      call prolong_coarse(self%cidx, ndim, dvec, pvec)
      call fsai_raw_apply(self, ndim, pvec, gres, corr)
      zvec(:) = zvec - corr
   end if

end subroutine fsai_apply


!> Apply G^T G to a vector
subroutine fsai_raw_apply(self, ndim, rvec, gres, zvec)

   !> Factorized sparse approximate inverse instance
   type(fsai_precond), intent(in) :: self

   !> Dimension of the linear system
   integer, intent(in) :: ndim

   !> Input vector
   real(wp), intent(in) :: rvec(:)

   !> Intermediate vector G * rvec
   real(wp), intent(inout) :: gres(:)

   !> Output vector G^T * gres
   real(wp), intent(out) :: zvec(:)

   integer :: iat
   integer(i8) :: ij
   real(wp) :: tmp

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

end subroutine fsai_raw_apply


!> Restrict a fine-grid vector to coarse cells by summation
pure subroutine restrict_coarse(cidx, ndim, ncoarse, fvec, cvec)

   !> Coarse cell index of each atom
   integer, intent(in) :: cidx(:)

   !> Number of atoms
   integer, intent(in) :: ndim

   !> Number of coarse cells
   integer, intent(in) :: ncoarse

   !> Fine-grid vector
   real(wp), intent(in) :: fvec(:)

   !> Coarse-grid vector
   real(wp), intent(out) :: cvec(:)

   integer :: iat, ic

   cvec(:ncoarse) = 0.0_wp
   do iat = 1, ndim
      ic = cidx(iat)
      cvec(ic) = cvec(ic) + fvec(iat)
   end do

end subroutine restrict_coarse


!> Prolongate a coarse-grid vector to atoms
pure subroutine prolong_coarse(cidx, ndim, cvec, fvec)

   !> Coarse cell index of each atom
   integer, intent(in) :: cidx(:)

   !> Number of atoms
   integer, intent(in) :: ndim

   !> Coarse-grid vector
   real(wp), intent(in) :: cvec(:)

   !> Fine-grid vector
   real(wp), intent(out) :: fvec(:)

   integer :: iat

   do iat = 1, ndim
      fvec(iat) = cvec(cidx(iat))
   end do

end subroutine prolong_coarse


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
