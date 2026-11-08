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

!> @file multicharge/fft.f90
!> Provides a self-contained mixed-radix fast Fourier transform

!> Three-dimensional complex fast Fourier transform.
!>
!> The one-dimensional transforms use a self-sorting Stockham formulation with
!> radix 2, 3, 4 and 5 stages, which is applied to a block of lines at once,
!> such that the innermost loops run over contiguous memory. Transforms are
!> unnormalized, the forward transform uses exp(-2 pi i jk/n).
!>
!> P. N. Swarztrauber, *Parallel Comput.*, **1984**, 1, 45.
!> DOI: [10.1016/S0167-8191(84)90413-7](https://doi.org/10.1016/S0167-8191(84)90413-7)
module multicharge_fft
   use mctc_env, only : wp
   use mctc_io_constants, only : pi
   implicit none
   private

   public :: fft_plan, new_fft_plan, fft_3d, good_fft_size


   !> Plan for the one-dimensional transforms of a single length
   type :: fft_plan_1d
      !> Transform length
      integer :: n = 1

      !> Radix of each stage
      integer, allocatable :: radix(:)

      !> Offset of each stage into the twiddle tables
      integer, allocatable :: offset(:)

      !> Twiddle factors of the forward transform
      complex(wp), allocatable :: twf(:)

      !> Twiddle factors of the backward transform
      complex(wp), allocatable :: twb(:)
   end type fft_plan_1d

   !> Plan for three-dimensional complex transforms
   type :: fft_plan
      !> Grid dimensions
      integer :: n(3) = 1

      !> One-dimensional plans along each grid dimension
      type(fft_plan_1d) :: line(3)
   end type fft_plan

   !> Number of lines transformed together
   integer, parameter :: nblock = 16


contains


!> Smallest integer not below n with only the prime factors 2, 3 and 5
pure function good_fft_size(n) result(m)
   !> Requested minimal size
   integer, intent(in) :: n
   !> FFT friendly size
   integer :: m

   integer :: k

   m = max(n, 1)
   do
      k = m
      do while (mod(k, 2) == 0)
         k = k / 2
      end do
      do while (mod(k, 3) == 0)
         k = k / 3
      end do
      do while (mod(k, 5) == 0)
         k = k / 5
      end do
      if (k == 1) exit
      m = m + 1
   end do

end function good_fft_size


!> Create a plan for three-dimensional transforms of the given grid
subroutine new_fft_plan(self, n)
   !> Instance of the plan
   type(fft_plan), intent(out) :: self
   !> Grid dimensions
   integer, intent(in) :: n(3)

   integer :: idim

   self%n = n
   do idim = 1, 3
      call new_fft_plan_1d(self%line(idim), n(idim))
   end do

end subroutine new_fft_plan


!> Factorize the transform length and tabulate the twiddle factors
subroutine new_fft_plan_1d(self, n)
   !> Instance of the plan
   type(fft_plan_1d), intent(out) :: self
   !> Transform length
   integer, intent(in) :: n

   integer :: rest, nstage, istage, p, ns, m, j, t, pos, ntw
   integer :: stages(64)
   real(wp) :: arg

   self%n = n
   rest = n
   nstage = 0
   do while (mod(rest, 4) == 0)
      nstage = nstage + 1
      stages(nstage) = 4
      rest = rest / 4
   end do
   do p = 2, n
      if (rest == 1) exit
      do while (mod(rest, p) == 0)
         nstage = nstage + 1
         stages(nstage) = p
         rest = rest / p
      end do
   end do
   self%radix = stages(:nstage)

   ! Twiddle table w^(j*t) of each stage, with w = exp(-2 pi i / ns)
   allocate(self%offset(nstage))
   ntw = 0
   ns = n
   do istage = 1, nstage
      p = self%radix(istage)
      m = ns / p
      self%offset(istage) = ntw
      ntw = ntw + m * (p - 1)
      ns = m
   end do
   allocate(self%twf(max(ntw, 1)), self%twb(max(ntw, 1)))
   ns = n
   do istage = 1, nstage
      p = self%radix(istage)
      m = ns / p
      do t = 1, p - 1
         do j = 0, m - 1
            pos = self%offset(istage) + j + 1 + m * (t - 1)
            arg = 2.0_wp * pi * real(j * t, wp) / real(ns, wp)
            self%twf(pos) = cmplx(cos(arg), -sin(arg), wp)
            self%twb(pos) = cmplx(cos(arg), sin(arg), wp)
         end do
      end do
      ns = m
   end do

end subroutine new_fft_plan_1d


!> In-place three-dimensional transform of a complex grid.
!>
!> The optional active extents restrict the work for zero-padded grids. In a
!> forward transform the input must vanish outside of the leading active block,
!> in a backward transform only the leading active block of the output is
!> computed, the rest of the grid is left in an intermediate state.
subroutine fft_3d(self, data, backward, active)
   !> Instance of the plan
   type(fft_plan), intent(in) :: self
   !> Grid data, transformed in place
   complex(wp), intent(inout) :: data(:, :, :)
   !> Select the backward transform
   logical, intent(in) :: backward
   !> Leading active extent of the grid in each dimension
   integer, intent(in), optional :: active(3)

   integer :: act(3)

   act = self%n
   if (present(active)) act = min(active, self%n)

   if (backward) then
      call transform_dim3(self, data, backward, self%n(2))
      call transform_dim2(self, data, backward, act(3))
      call transform_dim1(self, data, backward, act(2), act(3))
   else
      call transform_dim1(self, data, backward, act(2), act(3))
      call transform_dim2(self, data, backward, act(3))
      call transform_dim3(self, data, backward, self%n(2))
   end if

end subroutine fft_3d


!> Transform along the first, contiguous dimension for the leading planes
subroutine transform_dim1(self, data, backward, n2, n3)
   !> Instance of the plan
   type(fft_plan), intent(in) :: self
   !> Grid data, transformed in place
   complex(wp), intent(inout) :: data(:, :, :)
   !> Select the backward transform
   logical, intent(in) :: backward
   !> Number of lines along the second dimension to transform
   integer, intent(in) :: n2
   !> Number of planes along the third dimension to transform
   integer, intent(in) :: n3

   integer :: n1, nb2, ib, jb, j0, nb, k, j
   complex(wp), allocatable :: xbuf(:, :), ybuf(:, :)

   n1 = self%n(1)
   nb2 = (n2 + nblock - 1) / nblock

   !$omp parallel default(none) &
   !$omp shared(self, data, backward, n1, n2, n3, nb2) &
   !$omp private(ib, jb, j0, nb, k, j, xbuf, ybuf)
   allocate(xbuf(nblock, n1), ybuf(nblock, n1))
   !$omp do schedule(static)
   do ib = 1, nb2 * n3
      k = (ib - 1) / nb2 + 1
      jb = mod(ib - 1, nb2)
      j0 = jb * nblock
      nb = min(nblock, n2 - j0)
      do j = 1, nb
         xbuf(j, :) = data(:, j0 + j, k)
      end do
      call fft_lines(self%line(1), nb, xbuf, ybuf, backward)
      do j = 1, nb
         data(:, j0 + j, k) = xbuf(j, :)
      end do
   end do
   !$omp end do
   deallocate(xbuf, ybuf)
   !$omp end parallel

end subroutine transform_dim1


!> Transform along the second dimension for the leading planes
subroutine transform_dim2(self, data, backward, n3)
   !> Instance of the plan
   type(fft_plan), intent(in) :: self
   !> Grid data, transformed in place
   complex(wp), intent(inout) :: data(:, :, :)
   !> Select the backward transform
   logical, intent(in) :: backward
   !> Number of planes along the third dimension to transform
   integer, intent(in) :: n3

   integer :: n1, n2, nb1, ib, i0, nb, k
   complex(wp), allocatable :: xbuf(:, :), ybuf(:, :)

   n1 = self%n(1)
   n2 = self%n(2)
   nb1 = (n1 + nblock - 1) / nblock

   !$omp parallel default(none) &
   !$omp shared(self, data, backward, n1, n2, n3, nb1) &
   !$omp private(ib, i0, nb, k, xbuf, ybuf)
   allocate(xbuf(nblock, n2), ybuf(nblock, n2))
   !$omp do schedule(static)
   do ib = 1, nb1 * n3
      k = (ib - 1) / nb1 + 1
      i0 = mod(ib - 1, nb1) * nblock
      nb = min(nblock, n1 - i0)
      xbuf(:nb, :) = data(i0 + 1:i0 + nb, :, k)
      call fft_lines(self%line(2), nb, xbuf, ybuf, backward)
      data(i0 + 1:i0 + nb, :, k) = xbuf(:nb, :)
   end do
   !$omp end do
   deallocate(xbuf, ybuf)
   !$omp end parallel

end subroutine transform_dim2


!> Transform along the third dimension for the leading planes of the second one
subroutine transform_dim3(self, data, backward, n2)
   !> Instance of the plan
   type(fft_plan), intent(in) :: self
   !> Grid data, transformed in place
   complex(wp), intent(inout) :: data(:, :, :)
   !> Select the backward transform
   logical, intent(in) :: backward
   !> Number of planes along the second dimension to transform
   integer, intent(in) :: n2

   integer :: n1, n3, nb1, ib, i0, nb, j
   complex(wp), allocatable :: xbuf(:, :), ybuf(:, :)

   n1 = self%n(1)
   n3 = self%n(3)
   nb1 = (n1 + nblock - 1) / nblock

   !$omp parallel default(none) &
   !$omp shared(self, data, backward, n1, n2, n3, nb1) &
   !$omp private(ib, i0, nb, j, xbuf, ybuf)
   allocate(xbuf(nblock, n3), ybuf(nblock, n3))
   !$omp do schedule(static)
   do ib = 1, nb1 * n2
      j = (ib - 1) / nb1 + 1
      i0 = mod(ib - 1, nb1) * nblock
      nb = min(nblock, n1 - i0)
      xbuf(:nb, :) = data(i0 + 1:i0 + nb, j, :)
      call fft_lines(self%line(3), nb, xbuf, ybuf, backward)
      data(i0 + 1:i0 + nb, j, :) = xbuf(:nb, :)
   end do
   !$omp end do
   deallocate(xbuf, ybuf)
   !$omp end parallel

end subroutine transform_dim3


!> Transform a block of lines stored with the line index running fastest
subroutine fft_lines(plan, nb, x, y, backward)
   !> One-dimensional plan
   type(fft_plan_1d), intent(in) :: plan
   !> Number of active lines in the block
   integer, intent(in) :: nb
   !> Lines to transform, the result is returned in place
   complex(wp), intent(inout) :: x(nblock, plan%n)
   !> Work space of the same shape as the lines
   complex(wp), intent(inout) :: y(nblock, plan%n)
   !> Select the backward transform
   logical, intent(in) :: backward

   integer :: istage, p, ns, m, s, ld, off
   logical :: in_x

   ld = nblock
   ns = plan%n
   s = 1
   in_x = .true.
   do istage = 1, size(plan%radix)
      p = plan%radix(istage)
      m = ns / p
      off = plan%offset(istage)
      if (backward) then
         if (in_x) then
            call fft_pass(p, ld, nb, s, m, x, y, plan%twb(off + 1:), .true.)
         else
            call fft_pass(p, ld, nb, s, m, y, x, plan%twb(off + 1:), .true.)
         end if
      else
         if (in_x) then
            call fft_pass(p, ld, nb, s, m, x, y, plan%twf(off + 1:), .false.)
         else
            call fft_pass(p, ld, nb, s, m, y, x, plan%twf(off + 1:), .false.)
         end if
      end if
      in_x = .not. in_x
      ns = m
      s = s * p
   end do
   if (.not. in_x) x(:nb, :) = y(:nb, :)

end subroutine fft_lines


!> Dispatch a single Stockham stage to the radix-specific kernel
subroutine fft_pass(p, ld, nb, s, m, x, y, tw, backward)
   !> Radix of the stage
   integer, intent(in) :: p
   !> Leading dimension of the line block
   integer, intent(in) :: ld
   !> Number of active lines
   integer, intent(in) :: nb
   !> Number of interleaved subsequences
   integer, intent(in) :: s
   !> Length of the subsequences after this stage
   integer, intent(in) :: m
   !> Input of the stage
   complex(wp), intent(in) :: x(ld, s, m, p)
   !> Output of the stage
   complex(wp), intent(out) :: y(ld, s, p, m)
   !> Twiddle factors of the stage
   complex(wp), intent(in) :: tw(m, p - 1)
   !> Select the backward transform
   logical, intent(in) :: backward

   select case(p)
   case(2)
      call pass2(ld, nb, s, m, x, y, tw)
   case(3)
      call pass3(ld, nb, s, m, x, y, tw, backward)
   case(4)
      call pass4(ld, nb, s, m, x, y, tw, backward)
   case(5)
      call pass5(ld, nb, s, m, x, y, tw, backward)
   case default
      call passg(p, ld, nb, s, m, x, y, tw, backward)
   end select

end subroutine fft_pass


!> Radix-2 Stockham stage
pure subroutine pass2(ld, nb, s, m, x, y, tw)
   !> Leading dimension of the line block
   integer, intent(in) :: ld
   !> Number of active lines
   integer, intent(in) :: nb
   !> Number of interleaved subsequences
   integer, intent(in) :: s
   !> Length of the subsequences after this stage
   integer, intent(in) :: m
   !> Input of the stage
   complex(wp), intent(in) :: x(ld, s, m, 2)
   !> Output of the stage
   complex(wp), intent(out) :: y(ld, s, 2, m)
   !> Twiddle factors of the stage
   complex(wp), intent(in) :: tw(m, 1)

   integer :: j, q, b
   complex(wp) :: w1, a0, a1

   do j = 1, m
      w1 = tw(j, 1)
      do q = 1, s
         do b = 1, nb
            a0 = x(b, q, j, 1)
            a1 = x(b, q, j, 2)
            y(b, q, 1, j) = a0 + a1
            y(b, q, 2, j) = (a0 - a1) * w1
         end do
      end do
   end do

end subroutine pass2


!> Radix-3 Stockham stage
pure subroutine pass3(ld, nb, s, m, x, y, tw, backward)
   !> Leading dimension of the line block
   integer, intent(in) :: ld
   !> Number of active lines
   integer, intent(in) :: nb
   !> Number of interleaved subsequences
   integer, intent(in) :: s
   !> Length of the subsequences after this stage
   integer, intent(in) :: m
   !> Input of the stage
   complex(wp), intent(in) :: x(ld, s, m, 3)
   !> Output of the stage
   complex(wp), intent(out) :: y(ld, s, 3, m)
   !> Twiddle factors of the stage
   complex(wp), intent(in) :: tw(m, 2)
   !> Select the backward transform
   logical, intent(in) :: backward

   real(wp), parameter :: s3 = 0.5_wp * sqrt(3.0_wp)
   integer :: j, q, b
   real(wp) :: sg
   complex(wp) :: w1, w2, a0, a1, a2, t1, t2, t3

   sg = merge(s3, -s3, backward)
   do j = 1, m
      w1 = tw(j, 1)
      w2 = tw(j, 2)
      do q = 1, s
         do b = 1, nb
            a0 = x(b, q, j, 1)
            a1 = x(b, q, j, 2)
            a2 = x(b, q, j, 3)
            t1 = a1 + a2
            t2 = a0 - 0.5_wp * t1
            t3 = a1 - a2
            t3 = cmplx(-sg * aimag(t3), sg * real(t3), wp)
            y(b, q, 1, j) = a0 + t1
            y(b, q, 2, j) = (t2 + t3) * w1
            y(b, q, 3, j) = (t2 - t3) * w2
         end do
      end do
   end do

end subroutine pass3


!> Radix-4 Stockham stage
pure subroutine pass4(ld, nb, s, m, x, y, tw, backward)
   !> Leading dimension of the line block
   integer, intent(in) :: ld
   !> Number of active lines
   integer, intent(in) :: nb
   !> Number of interleaved subsequences
   integer, intent(in) :: s
   !> Length of the subsequences after this stage
   integer, intent(in) :: m
   !> Input of the stage
   complex(wp), intent(in) :: x(ld, s, m, 4)
   !> Output of the stage
   complex(wp), intent(out) :: y(ld, s, 4, m)
   !> Twiddle factors of the stage
   complex(wp), intent(in) :: tw(m, 3)
   !> Select the backward transform
   logical, intent(in) :: backward

   integer :: j, q, b
   real(wp) :: sg
   complex(wp) :: w1, w2, w3, a0, a1, a2, a3, t0, t1, t2, t3

   sg = merge(1.0_wp, -1.0_wp, backward)
   do j = 1, m
      w1 = tw(j, 1)
      w2 = tw(j, 2)
      w3 = tw(j, 3)
      do q = 1, s
         do b = 1, nb
            a0 = x(b, q, j, 1)
            a1 = x(b, q, j, 2)
            a2 = x(b, q, j, 3)
            a3 = x(b, q, j, 4)
            t0 = a0 + a2
            t1 = a0 - a2
            t2 = a1 + a3
            t3 = a1 - a3
            t3 = cmplx(-sg * aimag(t3), sg * real(t3), wp)
            y(b, q, 1, j) = t0 + t2
            y(b, q, 2, j) = (t1 + t3) * w1
            y(b, q, 3, j) = (t0 - t2) * w2
            y(b, q, 4, j) = (t1 - t3) * w3
         end do
      end do
   end do

end subroutine pass4


!> Radix-5 Stockham stage
pure subroutine pass5(ld, nb, s, m, x, y, tw, backward)
   !> Leading dimension of the line block
   integer, intent(in) :: ld
   !> Number of active lines
   integer, intent(in) :: nb
   !> Number of interleaved subsequences
   integer, intent(in) :: s
   !> Length of the subsequences after this stage
   integer, intent(in) :: m
   !> Input of the stage
   complex(wp), intent(in) :: x(ld, s, m, 5)
   !> Output of the stage
   complex(wp), intent(out) :: y(ld, s, 5, m)
   !> Twiddle factors of the stage
   complex(wp), intent(in) :: tw(m, 4)
   !> Select the backward transform
   logical, intent(in) :: backward

   real(wp), parameter :: c1 = cos(0.4_wp * pi), c2 = cos(0.8_wp * pi)
   real(wp), parameter :: s1 = sin(0.4_wp * pi), s2 = sin(0.8_wp * pi)
   integer :: j, q, b
   real(wp) :: sg
   complex(wp) :: w1, w2, w3, w4, a0, a1, a2, a3, a4
   complex(wp) :: b1, b2, d1, d2, e1, e2, f1, f2

   sg = merge(1.0_wp, -1.0_wp, backward)
   do j = 1, m
      w1 = tw(j, 1)
      w2 = tw(j, 2)
      w3 = tw(j, 3)
      w4 = tw(j, 4)
      do q = 1, s
         do b = 1, nb
            a0 = x(b, q, j, 1)
            a1 = x(b, q, j, 2)
            a2 = x(b, q, j, 3)
            a3 = x(b, q, j, 4)
            a4 = x(b, q, j, 5)
            b1 = a1 + a4
            b2 = a2 + a3
            d1 = a1 - a4
            d2 = a2 - a3
            e1 = a0 + c1 * b1 + c2 * b2
            e2 = a0 + c2 * b1 + c1 * b2
            f1 = s1 * d1 + s2 * d2
            f2 = s2 * d1 - s1 * d2
            f1 = cmplx(-sg * aimag(f1), sg * real(f1), wp)
            f2 = cmplx(-sg * aimag(f2), sg * real(f2), wp)
            y(b, q, 1, j) = a0 + b1 + b2
            y(b, q, 2, j) = (e1 + f1) * w1
            y(b, q, 3, j) = (e2 + f2) * w2
            y(b, q, 4, j) = (e2 - f2) * w3
            y(b, q, 5, j) = (e1 - f1) * w4
         end do
      end do
   end do

end subroutine pass5


!> Generic Stockham stage for an arbitrary radix
pure subroutine passg(p, ld, nb, s, m, x, y, tw, backward)
   !> Radix of the stage
   integer, intent(in) :: p
   !> Leading dimension of the line block
   integer, intent(in) :: ld
   !> Number of active lines
   integer, intent(in) :: nb
   !> Number of interleaved subsequences
   integer, intent(in) :: s
   !> Length of the subsequences after this stage
   integer, intent(in) :: m
   !> Input of the stage
   complex(wp), intent(in) :: x(ld, s, m, p)
   !> Output of the stage
   complex(wp), intent(out) :: y(ld, s, p, m)
   !> Twiddle factors of the stage
   complex(wp), intent(in) :: tw(m, p - 1)
   !> Select the backward transform
   logical, intent(in) :: backward

   integer :: j, q, b, r, t
   real(wp) :: sg, arg
   complex(wp) :: acc, wt
   complex(wp) :: omega(0:p - 1)

   sg = merge(1.0_wp, -1.0_wp, backward)
   do r = 0, p - 1
      arg = 2.0_wp * pi * real(r, wp) / real(p, wp)
      omega(r) = cmplx(cos(arg), sg * sin(arg), wp)
   end do
   do j = 1, m
      do q = 1, s
         do t = 0, p - 1
            wt = (1.0_wp, 0.0_wp)
            if (t > 0) wt = tw(j, t)
            do b = 1, nb
               acc = x(b, q, j, 1)
               do r = 1, p - 1
                  acc = acc + x(b, q, j, r + 1) * omega(mod(r * t, p))
               end do
               y(b, q, t + 1, j) = acc * wt
            end do
         end do
      end do
   end do

end subroutine passg

end module multicharge_fft
