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

!> @file multicharge/coulomb/pme.f90
!> Provides the particle mesh Ewald representation of the Coulomb matrix

!> Smooth particle mesh Ewald summation for Gaussian charges.
!>
!> The Coulomb matrix is split into a sparse real-space part, stored in the
!> CSR list layout, and the smooth kernel erf(alpha r)/r, which is applied
!> on a grid with cardinal B-spline charge assignment.
!>
!> U. Essmann, L. Perera, M. L. Berkowitz, T. Darden, H. Lee and
!> L. G. Pedersen, *J. Chem. Phys.*, **1995**, 103, 8577.
!> DOI: [10.1063/1.470117](https://doi.org/10.1063/1.470117)
!>
!> Non-periodic systems use the smooth Fourier extension of the long-range
!> kernel on a doubled domain.
!>
!> D. M. Elking, *Comput. Phys. Commun.*, **2025**, 315, 109739.
!> DOI: [10.1016/j.cpc.2025.109739](https://doi.org/10.1016/j.cpc.2025.109739)
module multicharge_coulomb_pme
   use, intrinsic :: iso_fortran_env, only : output_unit
   use mctc_csrlist, only : csr_list, new_csr_list, spsymv_csr
   use mctc_env, only : error_type, fatal_error, wp, i8
   use mctc_io, only : structure_type
   use mctc_io_constants, only : pi
   use mctc_io_math, only : matdet_3x3, matinv_3x3
   use multicharge_coulomb_shortrange, only : get_gamma, get_pair_values, &
      & add_pair_gradient
   use multicharge_coulomb_type, only : coulomb_operator, coulomb_input
   use multicharge_ewald, only : get_dir_cutoff
   use multicharge_fft, only : fft_plan, new_fft_plan, fft_3d, good_fft_size
   use multicharge_fft_real, only : fft_rplan, new_fft_rplan, fft_rwork, &
      & fft_r2c_3d, fft_c2r_3d
   implicit none
   private

   public :: pme_operator, new_pme_operator


   !> Particle mesh Ewald representation of the Coulomb matrix
   type, extends(coulomb_operator) :: pme_operator
      !> Periodic boundary conditions
      logical :: periodic = .false.

      !> Order of the cardinal B-splines
      integer :: order = 8

      !> Ewald splitting parameter
      real(wp) :: alpha = 0.0_wp

      !> Real-space cutoff
      real(wp) :: cutoff = 0.0_wp

      !> Number of grid points along each grid axis
      integer :: ngrid(3) = 1

      !> Leading block of the grid holding the charges
      integer :: active(3) = 1

      !> Transformation from Cartesian to grid coordinates
      real(wp) :: frac(3, 3) = 0.0_wp

      !> Cartesian origin of the grid
      real(wp) :: origin(3) = 0.0_wp

      !> Reciprocal lattice vectors without the factor 2 pi
      real(wp) :: rec_lat(3, 3) = 0.0_wp

      !> Volume of the unit cell
      real(wp) :: volume = 0.0_wp

      !> Constant added to all elements to match the reference Ewald parameter
      real(wp) :: shift = 0.0_wp

      !> Plan of the complex grid transformation
      type(fft_plan) :: plan

      !> Plan of the real grid transformation applying the operator
      type(fft_rplan) :: rplan

      !> Influence function on the grid, including the B-spline moduli
      real(wp), allocatable :: influence(:, :, :)

      !> Cartesian coordinates
      real(wp), allocatable :: xyz(:, :)

      !> Species index of each atom
      integer, allocatable :: id(:)

      !> Pair-wise Gaussian width parameters
      real(wp), allocatable :: gam(:, :)

      !> Diagonal elements of the Coulomb matrix without lattice sums
      real(wp), allocatable :: diag(:)

      !> First grid point covered by the B-splines of each atom
      integer, allocatable :: gidx(:, :)

      !> B-spline weights of each atom along the grid axes
      real(wp), allocatable :: theta(:, :, :)

      !> Number of blocks of the colored charge assignment
      integer :: nblock(2) = 1

      !> Offset of each grid block into the block-ordered atom list
      integer, allocatable :: block_start(:)

      !> Atom indices ordered by grid block
      integer, allocatable :: block_atom(:)

      !> Color of each grid block, blocks of one color do not overlap
      integer, allocatable :: block_color(:)

      !> Neighbour list of the real-space part
      type(csr_list) :: list

      !> Real-space part in the CSR list layout
      real(wp), allocatable :: val(:)
   contains
      !> Multiply the Coulomb matrix with a trial vector
      procedure :: multiply
      !> Diagonal of the Coulomb matrix for preconditioning
      procedure :: get_diagonal
      !> Approximate dense block of the Coulomb matrix for preconditioning
      procedure :: get_block
      !> Contract the operator derivatives with two charge-like vectors
      procedure :: add_gradient
      !> Write a summary of the operator setup
      procedure :: write_info
   end type pme_operator

   real(wp), parameter :: sqrtpi = sqrt(pi)
   real(wp), parameter :: sqrt2pi = sqrt(2.0_wp / pi)

   !> Steepness of the ramp function of the Fourier extension
   real(wp), parameter :: ramp_steepness = 4.0_wp

   !> Default real-space cutoff in Bohr
   real(wp), parameter :: default_cutoff = 12.0_wp


contains


!> Construct the particle mesh Ewald representation of the Coulomb matrix
subroutine new_pme_operator(self, mol, error, input, rad, eta, alpha_ref)
   !> Instance of the operator
   type(pme_operator), intent(out) :: self
   !> Molecular structure data
   type(structure_type), intent(in) :: mol
   !> Error handling
   type(error_type), allocatable, intent(out) :: error
   !> Input configuration
   type(coulomb_input), intent(in) :: input
   !> Species-wise charge widths
   real(wp), intent(in) :: rad(:)
   !> Species-wise chemical hardness
   real(wp), intent(in) :: eta(:)
   !> Ewald parameter of the reference dense matrix in periodic systems
   real(wp), intent(in) :: alpha_ref

   integer :: iat
   real(wp) :: tol, rgauss, spacing

   if (any(mol%periodic) .and. .not. all(mol%periodic)) then
      call fatal_error(error, "Particle mesh Ewald requires molecular or 3D "//&
         & "periodic boundary conditions")
      return
   end if

   self%nat = mol%nat
   self%periodic = any(mol%periodic)
   self%xyz = mol%xyz
   self%id = mol%id
   call get_gamma(rad, self%gam)

   tol = input%accuracy
   if (allocated(input%order)) then
      self%order = input%order
   else
      self%order = default_order(tol)
   end if
   if (self%order < 4 .or. mod(self%order, 2) /= 0) then
      call fatal_error(error, "B-spline order of the PME must be even and at least 4")
      return
   end if

   ! The real-space kernel decays with the smaller of the Ewald parameter and
   ! the widest Gaussian, the cutoff cannot be smaller than the latter range
   rgauss = get_dir_cutoff(minval(self%gam), tol)
   if (allocated(input%cutoff)) then
      self%cutoff = max(input%cutoff, rgauss)
   else
      self%cutoff = max(default_cutoff, rgauss)
   end if
   self%alpha = get_alpha_cutoff(self%cutoff, tol)

   if (allocated(input%spacing)) then
      spacing = input%spacing
   else
      spacing = default_spacing(self%alpha, tol, self%order)
   end if

   if (self%periodic) then
      call setup_periodic_grid(self, mol, spacing, alpha_ref)
   else
      call setup_free_grid(self, mol, spacing)
   end if

   ! The work grids of the real transforms are allocated after the setup
   ! released its temporary grids, which lowers the peak memory
   call new_fft_rplan(self%rplan, self%ngrid, self%active)
   call setup_splines(self)

   ! Real-space part with the self-interaction of the smooth kernel removed
   allocate(self%diag(mol%nat))
   do iat = 1, mol%nat
      self%diag(iat) = eta(mol%id(iat)) + sqrt2pi / rad(mol%id(iat))
   end do
   call new_csr_list(self%list, mol, error, cutoff=self%cutoff)
   if (allocated(error)) return
   allocate(self%val(size(self%list%nlat)))
   call get_pair_values(self%list, self%xyz, self%id, self%gam, self%alpha, &
      & self%diag - 2.0_wp * self%alpha / sqrtpi, self%val)

end subroutine new_pme_operator


!> Default B-spline order for a target accuracy
pure function default_order(tol) result(order)
   !> Target accuracy
   real(wp), intent(in) :: tol
   !> B-spline order
   integer :: order

   if (tol >= 1.0e-5_wp) then
      order = 6
   else if (tol >= 1.0e-9_wp) then
      order = 8
   else
      order = 10
   end if

end function default_order


!> Default grid spacing for the Ewald parameter, target accuracy and order.
!> The Gaussian of the reciprocal kernel must have decayed at the Nyquist
!> frequency, the remaining factor accounts for the B-spline aliasing.
pure function default_spacing(alpha, tol, order) result(spacing)
   !> Ewald splitting parameter
   real(wp), intent(in) :: alpha
   !> Target accuracy
   real(wp), intent(in) :: tol
   !> B-spline order
   integer, intent(in) :: order
   !> Grid spacing
   real(wp) :: spacing

   real(wp) :: aliasing

   aliasing = min(0.95_wp, 0.55_wp + 0.05_wp * real(order, wp))
   spacing = aliasing * pi / (2.0_wp * alpha * sqrt(log(1.0_wp / tol)))

end function default_spacing


!> Ewald parameter for which the real-space kernel decays below the
!> tolerance at the cutoff
pure function get_alpha_cutoff(cutoff, tol) result(alpha)
   !> Real-space cutoff
   real(wp), intent(in) :: cutoff
   !> Target accuracy
   real(wp), intent(in) :: tol
   !> Ewald splitting parameter
   real(wp) :: alpha

   integer :: iter
   real(wp) :: lo, hi

   lo = 0.0_wp
   hi = 1.0_wp
   do while (erfc(hi * cutoff) / cutoff > tol)
      hi = 2.0_wp * hi
   end do
   do iter = 1, 100
      alpha = 0.5_wp * (lo + hi)
      if (erfc(alpha * cutoff) / cutoff > tol) then
         lo = alpha
      else
         hi = alpha
      end if
      if (hi - lo < 1.0e-12_wp * hi) exit
   end do
   alpha = hi

end function get_alpha_cutoff


!> Smallest FFT friendly even grid size not below n
pure function good_grid_size(n) result(m)
   !> Requested minimal size
   integer, intent(in) :: n
   !> Even FFT friendly grid size
   integer :: m

   m = good_fft_size(n)
   do while (mod(m, 2) /= 0)
      m = good_fft_size(m + 1)
   end do

end function good_grid_size


!> Grid and influence function for a periodic system
subroutine setup_periodic_grid(self, mol, spacing, alpha_ref)
   !> Instance of the operator
   type(pme_operator), intent(inout) :: self
   !> Molecular structure data
   type(structure_type), intent(in) :: mol
   !> Maximal grid spacing
   real(wp), intent(in) :: spacing
   !> Ewald parameter of the reference dense matrix
   real(wp), intent(in) :: alpha_ref

   integer :: idir, i1, i2, i3, m(3)
   real(wp) :: inv_lat(3, 3), gvec(3), g2, fac
   real(wp), allocatable :: bsp(:, :)

   do idir = 1, 3
      self%ngrid(idir) = good_grid_size(max(2 * self%order, &
         & ceiling(norm2(mol%lattice(:, idir)) / spacing)))
   end do
   self%active = self%ngrid
   self%volume = abs(matdet_3x3(mol%lattice))
   inv_lat = matinv_3x3(mol%lattice)
   self%rec_lat = transpose(inv_lat)
   do idir = 1, 3
      self%frac(idir, :) = real(self%ngrid(idir), wp) * inv_lat(idir, :)
   end do
   self%origin = 0.0_wp
   self%shift = pi / self%volume * (1.0_wp / alpha_ref**2 - 1.0_wp / self%alpha**2)

   call new_fft_plan(self%plan, self%ngrid)
   call get_bspline_moduli(self%order, self%ngrid, bsp)

   fac = 4.0_wp * pi / self%volume
   allocate(self%influence(self%ngrid(1), self%ngrid(2), self%ngrid(3)))
   !$omp parallel do default(none) collapse(2) schedule(static) &
   !$omp shared(self, bsp, fac) private(i1, i2, i3, m, gvec, g2)
   do i3 = 1, self%ngrid(3)
      do i2 = 1, self%ngrid(2)
         do i1 = 1, self%ngrid(1)
            m = alias_index([i1, i2, i3], self%ngrid)
            gvec = 2.0_wp * pi * matmul(self%rec_lat, real(m, wp))
            g2 = dot_product(gvec, gvec)
            if (all(m == 0)) then
               self%influence(i1, i2, i3) = 0.0_wp
            else
               self%influence(i1, i2, i3) = fac * exp(-0.25_wp * g2 / self%alpha**2) &
                  & / g2 * bsp(i1, 1) * bsp(i2, 2) * bsp(i3, 3)
            end if
         end do
      end do
   end do

end subroutine setup_periodic_grid


!> Grid and influence function from the Fourier extension for a molecule
subroutine setup_free_grid(self, mol, spacing)
   !> Instance of the operator
   type(pme_operator), intent(inout) :: self
   !> Molecular structure data
   type(structure_type), intent(in) :: mol
   !> Maximal grid spacing
   real(wp), intent(in) :: spacing

   integer :: idir, i1, i2, i3, j1, j2, j3
   real(wp) :: xmin(3), xmax(3), exact(3), smooth, length(3), step(3), vec(3)
   real(wp) :: r, f, wgt
   real(wp), allocatable :: pos(:, :, :), wgt_ramp(:, :, :), bsp(:, :)
   complex(wp), allocatable :: grid(:, :, :)

   xmin = minval(mol%xyz, dim=2)
   xmax = maxval(mol%xyz, dim=2)

   ! The kernel must be exact for all differences between grid points
   ! carrying charge, including the B-spline support and the sharpening by
   ! the B-spline moduli, the smoothing zone transitions to the periodic image
   smooth = self%cutoff
   do idir = 1, 3
      exact(idir) = xmax(idir) - xmin(idir) + real(2 * self%order + 2, wp) * spacing
      length(idir) = 2.0_wp * (exact(idir) + smooth)
      self%ngrid(idir) = good_grid_size(max(2 * self%order, &
         & ceiling(length(idir) / spacing)))
      step(idir) = length(idir) / real(self%ngrid(idir), wp)
      self%origin(idir) = xmin(idir) - real(self%order, wp) * step(idir)
      self%active(idir) = min(self%ngrid(idir), &
         & int((xmax(idir) - self%origin(idir)) / step(idir)) + 2)
   end do
   self%frac = 0.0_wp
   do idir = 1, 3
      self%frac(idir, idir) = 1.0_wp / step(idir)
   end do
   self%volume = product(length)
   self%shift = 0.0_wp

   call new_fft_plan(self%plan, self%ngrid)
   call get_bspline_moduli(self%order, self%ngrid, bsp)

   ! Positions and ramp weights of the kernel and its periodic image
   allocate(pos(maxval(self%ngrid), 2, 3), source=0.0_wp)
   allocate(wgt_ramp(maxval(self%ngrid), 2, 3), source=0.0_wp)
   do idir = 1, 3
      call get_extension(self%ngrid(idir), length(idir), exact(idir), &
         & pos(:, :, idir), wgt_ramp(:, :, idir))
   end do

   ! Sample the smooth extension of erf(alpha r)/r and transform it
   allocate(grid(self%ngrid(1), self%ngrid(2), self%ngrid(3)))
   !$omp parallel do default(none) collapse(2) schedule(static) &
   !$omp shared(self, grid, pos, wgt_ramp) &
   !$omp private(i1, i2, i3, j1, j2, j3, vec, r, f, wgt)
   do i3 = 1, self%ngrid(3)
      do i2 = 1, self%ngrid(2)
         do i1 = 1, self%ngrid(1)
            f = 0.0_wp
            do j3 = 1, 2
               if (wgt_ramp(i3, j3, 3) == 0.0_wp) cycle
               do j2 = 1, 2
                  if (wgt_ramp(i2, j2, 2) == 0.0_wp) cycle
                  do j1 = 1, 2
                     wgt = wgt_ramp(i1, j1, 1) * wgt_ramp(i2, j2, 2) * wgt_ramp(i3, j3, 3)
                     if (wgt == 0.0_wp) cycle
                     vec = [pos(i1, j1, 1), pos(i2, j2, 2), pos(i3, j3, 3)]
                     r = norm2(vec)
                     if (r > 1.0e-10_wp) then
                        f = f + wgt * erf(self%alpha * r) / r
                     else
                        f = f + wgt * 2.0_wp * self%alpha / sqrtpi
                     end if
                  end do
               end do
            end do
            grid(i1, i2, i3) = cmplx(f, 0.0_wp, wp)
         end do
      end do
   end do
   call fft_3d(self%plan, grid, .false.)

   allocate(self%influence(self%ngrid(1), self%ngrid(2), self%ngrid(3)))
   !$omp parallel do default(none) collapse(2) schedule(static) &
   !$omp shared(self, grid, bsp) private(i1, i2, i3)
   do i3 = 1, self%ngrid(3)
      do i2 = 1, self%ngrid(2)
         do i1 = 1, self%ngrid(1)
            self%influence(i1, i2, i3) = real(grid(i1, i2, i3), wp) &
               & / real(product(self%ngrid), wp) * bsp(i1, 1) * bsp(i2, 2) * bsp(i3, 3)
         end do
      end do
   end do

end subroutine setup_free_grid


!> Positions and weights of the overlapping smooth periodic extension along
!> one grid axis. Inside the exact region the kernel is sampled directly, in
!> the smoothing zone it is blended with its periodic image by a ramp.
pure subroutine get_extension(ngrid, length, exact, pos, wgt)
   !> Number of grid points
   integer, intent(in) :: ngrid
   !> Period of the grid
   real(wp), intent(in) :: length
   !> Half width of the exact region
   real(wp), intent(in) :: exact
   !> Position of the kernel and of its periodic image
   real(wp), intent(out) :: pos(:, :)
   !> Weight of the kernel and of its periodic image
   real(wp), intent(out) :: wgt(:, :)

   integer :: ip
   real(wp) :: theta, x, width

   theta = 0.5_wp * length
   width = theta - exact
   pos(:, :) = 0.0_wp
   wgt(:, :) = 0.0_wp
   do ip = 1, ngrid
      if (ip - 1 < ngrid / 2) then
         x = real(ip - 1, wp) * length / real(ngrid, wp)
      else
         x = real(ip - 1 - ngrid, wp) * length / real(ngrid, wp)
      end if
      pos(ip, 1) = x
      if (abs(x) <= exact) then
         wgt(ip, 1) = 1.0_wp
      else
         pos(ip, 2) = x - sign(length, x)
         wgt(ip, 1) = ramp((exact + width - abs(x)) / width)
         wgt(ip, 2) = ramp((exact + width - abs(pos(ip, 2))) / width)
      end if
   end do

end subroutine get_extension


!> Smooth ramp rising from zero at -1 to one at +1 with all derivatives
!> vanishing at the end points
elemental function ramp(x) result(val)
   !> Argument of the ramp
   real(wp), intent(in) :: x
   !> Value of the ramp
   real(wp) :: val

   if (x <= -1.0_wp) then
      val = 0.0_wp
   else if (x >= 1.0_wp) then
      val = 1.0_wp
   else
      val = 0.5_wp * (1.0_wp + erf(ramp_steepness * x / sqrt(1.0_wp - x * x)))
   end if

end function ramp


!> Frequency index of a grid point in the range [-n/2, n/2)
pure function alias_index(idx, ngrid) result(m)
   !> One-based grid index
   integer, intent(in) :: idx(3)
   !> Number of grid points
   integer, intent(in) :: ngrid(3)
   !> Frequency index
   integer :: m(3)

   m = idx - 1
   where (m >= ngrid / 2) m = m - ngrid

end function alias_index


!> Inverse squared moduli of the B-spline structure factor approximation
subroutine get_bspline_moduli(order, ngrid, bsp)
   !> B-spline order
   integer, intent(in) :: order
   !> Number of grid points
   integer, intent(in) :: ngrid(3)
   !> Inverse squared moduli along each grid axis
   real(wp), allocatable, intent(out) :: bsp(:, :)

   integer :: idir, im, k
   real(wp) :: theta(order), arg, sre, sim

   call fill_bspline(0.0_wp, order, theta)
   allocate(bsp(maxval(ngrid), 3), source=0.0_wp)
   do idir = 1, 3
      do im = 1, ngrid(idir)
         sre = 0.0_wp
         sim = 0.0_wp
         do k = 0, order - 2
            arg = 2.0_wp * pi * real((im - 1) * k, wp) / real(ngrid(idir), wp)
            sre = sre + theta(order - k - 1) * cos(arg)
            sim = sim + theta(order - k - 1) * sin(arg)
         end do
         bsp(im, idir) = 1.0_wp / (sre**2 + sim**2)
      end do
   end do

end subroutine get_bspline_moduli


!> Cardinal B-spline weights and derivatives for the fractional offset w,
!> theta(i) = M_n(w + n - i)
pure subroutine fill_bspline(w, order, theta, dtheta)
   !> Fractional part of the grid coordinate
   real(wp), intent(in) :: w
   !> B-spline order
   integer, intent(in) :: order
   !> B-spline weights
   real(wp), intent(out) :: theta(:)
   !> Derivatives of the B-spline weights
   real(wp), intent(out), optional :: dtheta(:)

   integer :: k, j
   real(wp) :: div

   theta(order) = 0.0_wp
   theta(2) = w
   theta(1) = 1.0_wp - w
   do k = 3, order - 1
      div = 1.0_wp / real(k - 1, wp)
      theta(k) = div * w * theta(k - 1)
      do j = 1, k - 2
         theta(k - j) = div * ((w + j) * theta(k - j - 1) + (k - j - w) * theta(k - j))
      end do
      theta(1) = div * (1.0_wp - w) * theta(1)
   end do
   if (present(dtheta)) then
      dtheta(1) = -theta(1)
      do j = 2, order
         dtheta(j) = theta(j - 1) - theta(j)
      end do
   end if
   div = 1.0_wp / real(order - 1, wp)
   theta(order) = div * w * theta(order - 1)
   do j = 1, order - 2
      theta(order - j) = div * ((w + j) * theta(order - j - 1) &
         & + (order - j - w) * theta(order - j))
   end do
   theta(1) = div * (1.0_wp - w) * theta(1)

end subroutine fill_bspline


!> Grid coordinates of an atom
pure function grid_coord(self, xyz) result(u)
   !> Instance of the operator
   type(pme_operator), intent(in) :: self
   !> Cartesian position
   real(wp), intent(in) :: xyz(3)
   !> Grid coordinate
   real(wp) :: u(3)

   u = matmul(self%frac, xyz - self%origin)
   if (self%periodic) then
      u = u - real(self%ngrid, wp) * floor(u / real(self%ngrid, wp))
      where (u >= real(self%ngrid, wp)) u = 0.0_wp
   end if

end function grid_coord


!> B-spline weights of all atoms and the colored blocks for charge assignment
subroutine setup_splines(self)
   !> Instance of the operator
   type(pme_operator), intent(inout) :: self

   integer :: iat, idir, ib, nb, ic(2), bsize(2)
   integer, allocatable :: atom_block(:), cnt(:)
   real(wp) :: u(3)

   allocate(self%gidx(3, self%nat), self%theta(self%order, 3, self%nat))
   !$omp parallel do default(none) schedule(static) &
   !$omp shared(self) private(iat, idir, u)
   do iat = 1, self%nat
      u = grid_coord(self, self%xyz(:, iat))
      do idir = 1, 3
         self%gidx(idir, iat) = int(u(idir)) - self%order + 1
         call fill_bspline(u(idir) - int(u(idir)), self%order, &
            & self%theta(:, idir, iat))
      end do
   end do

   ! Blocks along the second and third grid axis, each at least as thick as
   ! the B-spline support, an even number of blocks avoids wrap-around overlap
   do idir = 1, 2
      nb = 2 * (self%ngrid(idir + 1) / (2 * self%order))
      if (nb < 2) nb = 1
      self%nblock(idir) = nb
      bsize(idir) = self%ngrid(idir + 1) / nb
   end do
   allocate(atom_block(self%nat), cnt(product(self%nblock)), source=0)
   do iat = 1, self%nat
      do idir = 1, 2
         ic(idir) = min(modulo(self%gidx(idir + 1, iat), self%ngrid(idir + 1)) &
            & / bsize(idir), self%nblock(idir) - 1)
      end do
      atom_block(iat) = ic(1) + self%nblock(1) * ic(2) + 1
      cnt(atom_block(iat)) = cnt(atom_block(iat)) + 1
   end do
   allocate(self%block_start(product(self%nblock) + 1), &
      & self%block_color(product(self%nblock)), self%block_atom(self%nat))
   self%block_start(1) = 1
   do ib = 1, product(self%nblock)
      self%block_start(ib + 1) = self%block_start(ib) + cnt(ib)
      ic = [mod(ib - 1, self%nblock(1)), (ib - 1) / self%nblock(1)]
      self%block_color(ib) = mod(ic(1), 2) + 2 * mod(ic(2), 2)
   end do
   cnt(:) = 0
   do iat = 1, self%nat
      ib = atom_block(iat)
      self%block_atom(self%block_start(ib) + cnt(ib)) = iat
      cnt(ib) = cnt(ib) + 1
   end do

end subroutine setup_splines


!> Assign a charge-like vector to the grid
subroutine spread_charges(self, qvec, grid)
   !> Instance of the operator
   type(pme_operator), intent(in) :: self
   !> Charge-like vector
   real(wp), intent(in) :: qvec(:)
   !> Charge grid
   real(wp), intent(out) :: grid(:, :, :)

   integer :: icolor, ib, ia, iat, i1, i2, i3, n, i3w, i2w
   integer :: idx(self%order, 3)
   real(wp) :: w3, w23, chrg

   n = self%order
   !$omp parallel default(none) &
   !$omp shared(self, qvec, grid, n) &
   !$omp private(icolor, ib, ia, iat, i1, i2, i3, i2w, i3w, idx, w3, w23, chrg)
   !$omp do collapse(2) schedule(static)
   do i3 = 1, size(grid, 3)
      do i2 = 1, size(grid, 2)
         grid(:, i2, i3) = 0.0_wp
      end do
   end do
   !$omp end do
   do icolor = 0, 3
      !$omp do schedule(dynamic)
      do ib = 1, size(self%block_color)
         if (self%block_color(ib) /= icolor) cycle
         do ia = self%block_start(ib), self%block_start(ib + 1) - 1
            iat = self%block_atom(ia)
            chrg = qvec(iat)
            call get_indices(self, iat, idx)
            do i3 = 1, n
               i3w = idx(i3, 3)
               w3 = self%theta(i3, 3, iat)
               do i2 = 1, n
                  i2w = idx(i2, 2)
                  w23 = w3 * self%theta(i2, 2, iat)
                  do i1 = 1, n
                     grid(idx(i1, 1), i2w, i3w) = grid(idx(i1, 1), i2w, i3w) &
                        & + chrg * (w23 * self%theta(i1, 1, iat))
                  end do
               end do
            end do
         end do
      end do
      !$omp end do
   end do
   !$omp end parallel

end subroutine spread_charges


!> Wrapped grid indices covered by the B-splines of an atom
pure subroutine get_indices(self, iat, idx)
   !> Instance of the operator
   type(pme_operator), intent(in) :: self
   !> Atom index
   integer, intent(in) :: iat
   !> One-based grid indices along each axis
   integer, intent(out) :: idx(:, :)

   integer :: idir, k

   do idir = 1, 3
      do k = 1, self%order
         idx(k, idir) = modulo(self%gidx(idir, iat) + k - 1, self%ngrid(idir)) + 1
      end do
   end do

end subroutine get_indices


!> Convolve the charge grid with the long-range kernel
subroutine convolve(self, work)
   !> Instance of the operator
   type(pme_operator), intent(in) :: self
   !> Charge grid on input, potential grid on output
   type(fft_rwork), intent(inout) :: work

   call fft_r2c_3d(self%rplan, work)
   call apply_influence(self, work%cgrid)
   call fft_c2r_3d(self%rplan, work)

end subroutine convolve


!> Multiply the spectrum with the influence function. The influence function
!> is even, the half spectrum of a real transform uses its leading part along
!> the first dimension.
subroutine apply_influence(self, cgrid)
   !> Instance of the operator
   type(pme_operator), intent(in) :: self
   !> Spectrum of the charge grid
   complex(wp), intent(inout) :: cgrid(:, :, :)

   integer :: i2, i3, nspec

   nspec = size(cgrid, 1)
   !$omp parallel do default(none) collapse(2) schedule(static) &
   !$omp shared(self, cgrid, nspec) private(i2, i3)
   do i3 = 1, self%ngrid(3)
      do i2 = 1, self%ngrid(2)
         cgrid(:, i2, i3) = cgrid(:, i2, i3) * self%influence(:nspec, i2, i3)
      end do
   end do

end subroutine apply_influence


!> Set the real or the imaginary part of a complex grid
subroutine set_grid_part(grid, part, imaginary)
   !> Complex grid
   complex(wp), intent(inout) :: grid(:, :, :)
   !> Values of the part
   real(wp), intent(in) :: part(:, :, :)
   !> Set the imaginary instead of the real part
   logical, intent(in) :: imaginary

   integer :: i2, i3

   !$omp parallel do default(none) collapse(2) schedule(static) &
   !$omp shared(grid, part, imaginary) private(i2, i3)
   do i3 = 1, size(grid, 3)
      do i2 = 1, size(grid, 2)
         if (imaginary) then
            grid(:, i2, i3) = cmplx(real(grid(:, i2, i3), wp), part(:, i2, i3), wp)
         else
            grid(:, i2, i3) = cmplx(part(:, i2, i3), 0.0_wp, wp)
         end if
      end do
   end do

end subroutine set_grid_part


!> Add the potential of a grid at the atoms
subroutine gather_potential(self, grid, yvec)
   !> Instance of the operator
   type(pme_operator), intent(in) :: self
   !> Potential grid
   real(wp), intent(in) :: grid(:, :, :)
   !> Potential at the atoms, the grid potential is added
   real(wp), intent(inout) :: yvec(:)

   integer :: iat, i1, i2, i3, n
   integer :: idx(self%order, 3)
   real(wp) :: phi, w3, w23

   n = self%order
   !$omp parallel do default(none) schedule(static) &
   !$omp shared(self, grid, yvec, n) private(iat, i1, i2, i3, idx, phi, w3, w23)
   do iat = 1, self%nat
      call get_indices(self, iat, idx)
      phi = 0.0_wp
      do i3 = 1, n
         w3 = self%theta(i3, 3, iat)
         do i2 = 1, n
            w23 = w3 * self%theta(i2, 2, iat)
            do i1 = 1, n
               phi = phi + w23 * self%theta(i1, 1, iat) &
                  & * grid(idx(i1, 1), idx(i2, 2), idx(i3, 3))
            end do
         end do
      end do
      yvec(iat) = yvec(iat) + phi
   end do

end subroutine gather_potential


!> Multiply the Coulomb matrix with a trial vector
subroutine multiply(self, xvec, yvec)
   !> Instance of the operator
   class(pme_operator), intent(in) :: self
   !> Input vector
   real(wp), intent(in) :: xvec(:)
   !> Operator-vector product
   real(wp), intent(out) :: yvec(:)

   call spsymv_csr(self%nat, self%val, self%list%inl, self%list%nlat, xvec, yvec)

   ! The work arrays of the plan are reused, products of the same operator
   ! must not run concurrently
   call spread_charges(self, xvec, self%rplan%work%rgrid)
   call convolve(self, self%rplan%work)
   call gather_potential(self, self%rplan%work%rgrid, yvec)

   if (self%periodic) yvec(:) = yvec + self%shift * sum(xvec)

end subroutine multiply


!> Diagonal of the Coulomb matrix without the small lattice sums
subroutine get_diagonal(self, diag)
   !> Instance of the operator
   class(pme_operator), intent(in) :: self
   !> Diagonal of the operator
   real(wp), intent(out) :: diag(:)

   diag(:) = self%diag

end subroutine get_diagonal


!> Dense block of the Coulomb matrix without lattice sums.
!>
!> The pairs interact with erf(gam_ij r_ij)/r_ij, for periodic systems at the
!> minimum image distance. The block is exact for molecules and an
!> approximation for periodic systems, which remains positive definite as the
!> Gaussian Coulomb kernel is positive definite.
subroutine get_block(self, idx, block)
   !> Instance of the operator
   class(pme_operator), intent(in) :: self
   !> Row and column indices of the block
   integer, intent(in) :: idx(:)
   !> Matrix elements of the block
   real(wp), intent(out) :: block(:, :)

   integer :: ii, jj, iat, jat
   real(wp) :: lattice(3, 3), inv_lat(3, 3), vec(3), r

   if (self%periodic) then
      inv_lat = transpose(self%rec_lat)
      lattice = matinv_3x3(inv_lat)
   end if
   do jj = 1, size(idx)
      jat = idx(jj)
      block(jj, jj) = self%diag(jat)
      do ii = 1, jj - 1
         iat = idx(ii)
         vec = self%xyz(:, iat) - self%xyz(:, jat)
         if (self%periodic) then
            vec = matmul(inv_lat, vec)
            vec = matmul(lattice, vec - real(nint(vec), wp))
         end if
         r = norm2(vec)
         if (r > 1.0e-10_wp) then
            block(ii, jj) = erf(self%gam(self%id(iat), self%id(jat)) * r) / r
         else
            block(ii, jj) = 2.0_wp * self%gam(self%id(iat), self%id(jat)) / sqrtpi
         end if
         block(jj, ii) = block(ii, jj)
      end do
   end do

end subroutine get_block


!> Add scale * sum_ij p_i q_j dA_ij/dR to the gradient and the strain
!> derivative to sigma
subroutine add_gradient(self, pvec, qvec, gradient, sigma, scale, error)
   !> Instance of the operator
   class(pme_operator), intent(in) :: self
   !> First charge-like vector
   real(wp), intent(in) :: pvec(:)
   !> Second charge-like vector
   real(wp), intent(in) :: qvec(:)
   !> Cartesian gradient
   real(wp), intent(inout) :: gradient(:, :)
   !> Strain derivative
   real(wp), intent(inout) :: sigma(:, :)
   !> Scaling factor of the contraction
   real(wp), intent(in) :: scale
   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   integer :: iat, i1, i2, i3, n, idir, ic, m(3), mm(3)
   integer :: idx(self%order, 3)
   real(wp) :: u(3), dphi_p(3), dphi_q(3), val, wgt(3), dwgt(3)
   real(wp) :: gvec(3), g2, factor, pq, sigma_rec(3, 3)
   real(wp) :: theta(self%order, 3), dtheta(self%order, 3)
   real(wp), allocatable :: grad_loc(:, :), sigma_loc(:, :)
   complex(wp) :: zp, zm, qp, qq
   complex(wp), allocatable :: grid(:, :, :)

   allocate(grad_loc(3, self%nat), source=0.0_wp)
   allocate(sigma_loc(3, 3), source=0.0_wp)

   ! Real-space part
   call add_pair_gradient(self%list, self%xyz, self%id, self%gam, self%alpha, &
      & pvec, qvec, scale, grad_loc, sigma_loc)

   ! Both vectors are transformed together, Z = Q_q + i Q_p, the real work
   ! grid of the operator holds each charge grid in turn
   allocate(grid(self%ngrid(1), self%ngrid(2), self%ngrid(3)))
   call spread_charges(self, qvec, self%rplan%work%rgrid)
   call set_grid_part(grid, self%rplan%work%rgrid, .false.)
   call spread_charges(self, pvec, self%rplan%work%rgrid)
   call set_grid_part(grid, self%rplan%work%rgrid, .true.)
   call fft_3d(self%plan, grid, .false., self%active)

   ! Reciprocal strain derivative from the structure factors of both vectors
   if (self%periodic) then
      sigma_rec(:, :) = 0.0_wp
      !$omp parallel do default(none) collapse(2) schedule(static) &
      !$omp shared(self, grid) private(i1, i2, i3, m, mm, zp, zm, qp, qq, pq) &
      !$omp private(gvec, g2, factor, ic) reduction(+:sigma_rec)
      do i3 = 1, self%ngrid(3)
         do i2 = 1, self%ngrid(2)
            do i1 = 1, self%ngrid(1)
               m = alias_index([i1, i2, i3], self%ngrid)
               if (all(m == 0)) cycle
               mm = modulo(-m, self%ngrid) + 1
               zp = grid(i1, i2, i3)
               zm = conjg(grid(mm(1), mm(2), mm(3)))
               qq = 0.5_wp * (zp + zm)
               qp = cmplx(0.0_wp, -0.5_wp, wp) * (zp - zm)
               pq = real(qp * conjg(qq), wp)
               gvec = 2.0_wp * pi * matmul(self%rec_lat, real(m, wp))
               g2 = dot_product(gvec, gvec)
               factor = self%influence(i1, i2, i3) * pq
               do ic = 1, 3
                  sigma_rec(:, ic) = sigma_rec(:, ic) + factor &
                     & * (2.0_wp / g2 + 0.5_wp / self%alpha**2) * gvec * gvec(ic)
                  sigma_rec(ic, ic) = sigma_rec(ic, ic) - factor
               end do
            end do
         end do
      end do
      sigma_loc(:, :) = sigma_loc + scale * sigma_rec
      ! Strain derivative of the constant shift
      do ic = 1, 3
         sigma_loc(ic, ic) = sigma_loc(ic, ic) - scale * self%shift * sum(pvec) * sum(qvec)
      end do
   end if

   !$omp parallel do default(none) collapse(2) schedule(static) &
   !$omp shared(self, grid) private(i2, i3)
   do i3 = 1, self%ngrid(3)
      do i2 = 1, self%ngrid(2)
         grid(:, i2, i3) = grid(:, i2, i3) * self%influence(:, i2, i3)
      end do
   end do
   call fft_3d(self%plan, grid, .true., self%active)

   ! Gradients of both potentials at the atoms
   n = self%order
   !$omp parallel do default(none) schedule(static) &
   !$omp shared(self, grid, grad_loc, pvec, qvec, scale, n) &
   !$omp private(iat, i1, i2, i3, idir, idx, u, theta, dtheta, dphi_p, dphi_q) &
   !$omp private(val, wgt, dwgt)
   do iat = 1, self%nat
      u = grid_coord(self, self%xyz(:, iat))
      do idir = 1, 3
         call fill_bspline(u(idir) - int(u(idir)), n, theta(:, idir), dtheta(:, idir))
      end do
      call get_indices(self, iat, idx)
      dphi_p(:) = 0.0_wp
      dphi_q(:) = 0.0_wp
      do i3 = 1, n
         do i2 = 1, n
            do i1 = 1, n
               wgt = [theta(i1, 1), theta(i2, 2), theta(i3, 3)]
               dwgt = [dtheta(i1, 1) * wgt(2) * wgt(3), &
                  & wgt(1) * dtheta(i2, 2) * wgt(3), wgt(1) * wgt(2) * dtheta(i3, 3)]
               dphi_q = dphi_q + dwgt * real(grid(idx(i1, 1), idx(i2, 2), idx(i3, 3)), wp)
               dphi_p = dphi_p + dwgt * aimag(grid(idx(i1, 1), idx(i2, 2), idx(i3, 3)))
            end do
         end do
      end do
      ! Chain rule from grid to Cartesian coordinates
      dphi_p = matmul(dphi_p, self%frac)
      dphi_q = matmul(dphi_q, self%frac)
      grad_loc(:, iat) = grad_loc(:, iat) &
         & + scale * (pvec(iat) * dphi_q + qvec(iat) * dphi_p)
   end do

   ! The strain derivative of a molecule follows from the virial of the forces
   if (.not. self%periodic) then
      sigma_loc(:, :) = 0.0_wp
      do iat = 1, self%nat
         do ic = 1, 3
            sigma_loc(:, ic) = sigma_loc(:, ic) + self%xyz(:, iat) * grad_loc(ic, iat)
         end do
      end do
   end if

   gradient(:, :) = gradient + grad_loc
   sigma(:, :) = sigma + sigma_loc

end subroutine add_gradient


!> Write a summary of the operator setup
subroutine write_info(self, unit)
   !> Instance of the operator
   class(pme_operator), intent(in) :: self
   !> Output unit
   integer, intent(in) :: unit

   write(unit, "(a, t32, a)") "Coulomb operator", &
      & merge("particle mesh Ewald  ", "PME, Fourier extended", self%periodic)
   write(unit, "(a, t32, f12.6)") "Ewald parameter / bohr^-1", self%alpha
   write(unit, "(a, t32, f12.6)") "Real-space cutoff / bohr", self%cutoff
   write(unit, "(a, t32, i12)") "B-spline order", self%order
   write(unit, '(a, t32, 3(i0, :, " x "))') "Grid points", self%ngrid
   write(unit, '(a, t32, 3(i0, :, " x "))') "Grid points with charge", self%active
   write(unit, "(a, t32, a)") "Fourier transform", self%rplan%backend
   write(unit, "(a, t32, i12)") "Real-space pairs", size(self%list%nlat) - self%nat

end subroutine write_info

end module multicharge_coulomb_pme
