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

!> @file multicharge/fft_real.F90
!> Provides three-dimensional transforms of real grids

!> Real-to-complex and complex-to-real three-dimensional transforms.
!>
!> With FFTW the transforms of real grids only form the half spectrum
!> n(1)/2+1 x n(2) x n(3), which halves the work compared to complex
!> transforms. Lines and planes outside of the leading active block are
!> skipped, the input must vanish outside of the active block and only the
!> active block of the backward transform is formed. Without FFTW the complex
!> in-house transform is used with the full spectrum. Transforms are
!> unnormalized, the forward transform uses exp(-2 pi i jk/n).
!>
!> The FFTW plans are single-threaded and applied to planes and column blocks
!> in OpenMP loops, which keeps the threading under control of the caller.
!> FFTW plans depend on the alignment of the arrays, separate plans are kept
!> for each alignment class encountered in the grid, with an unaligned plan
!> as fallback.
!>
!> M. Frigo and S. G. Johnson, *Proc. IEEE*, **2005**, 93, 216.
!> DOI: [10.1109/JPROC.2004.840301](https://doi.org/10.1109/JPROC.2004.840301)
module multicharge_fft_real
#ifdef WITH_FFTW
   ! The FFTW interface declarations require the complete C binding module
   use, intrinsic :: iso_c_binding
   use mctc_env, only : wp
#else
   use mctc_env, only : wp
   use multicharge_fft, only : fft_plan, new_fft_plan, fft_3d
#endif
   implicit none
   private

#ifdef WITH_FFTW
#include "fftw3.f03"
#endif

   public :: fft_rplan, new_fft_rplan, fft_rwork, new_fft_rwork
   public :: fft_r2c_3d, fft_c2r_3d


#ifdef WITH_FFTW
   !> Number of alignment classes of 16 byte within the widest SIMD alignment
   integer, parameter :: nalign = 4

   !> Single-threaded FFTW plans for each alignment class, the last entry is
   !> the unaligned fallback
   type :: aligned_plans

      !> Plans for each alignment class
      type(c_ptr) :: plan(0:nalign) = c_null_ptr

   end type aligned_plans
#endif


   !> Plan for transforms of real grids
   type :: fft_rplan

      !> Grid dimensions
      integer :: n(3) = 1

      !> Leading active extent of the grid in each dimension
      integer :: active(3) = 1

      !> Extent of the spectrum along the first dimension
      integer :: nspec = 1

      !> Number of threads used for the transforms
      integer :: nthreads = 1

      !> Name of the transform backend
      character(len=:), allocatable :: backend

      !> Work arrays reused by all transforms of the plan, the transforms of a
      !> plan must therefore not run concurrently
      type(fft_rwork), pointer :: work => null()

#ifdef WITH_FFTW
      !> Real-to-complex transforms of the active lines of a plane
      type(aligned_plans) :: r2c_plane

      !> Complex-to-real transforms of the active lines of a plane
      type(aligned_plans) :: c2r_plane

      !> Forward transforms along the second dimension of a plane
      type(aligned_plans) :: fwd_plane

      !> Backward transforms along the second dimension of a plane
      type(aligned_plans) :: bwd_plane

      !> Forward transforms along the third dimension of a column block
      type(aligned_plans) :: fwd_column

      !> Backward transforms along the third dimension of a column block
      type(aligned_plans) :: bwd_column
#else
      !> Plan of the complex in-house transform
      type(fft_plan) :: plan
#endif
   contains

      !> Release the plans and the work arrays
      final :: delete_fft_rplan

   end type fft_rplan


   !> Work arrays of the transforms, aligned as required by the plans
   type :: fft_rwork

      !> Real grid of dimension n(1) x n(2) x n(3)
      real(wp), pointer, contiguous :: rgrid(:, :, :) => null()

      !> Spectrum of dimension nspec x n(2) x n(3)
      complex(wp), pointer, contiguous :: cgrid(:, :, :) => null()

#ifdef WITH_FFTW
      !> Memory of the real grid allocated by FFTW
      type(c_ptr) :: rmem = c_null_ptr

      !> Memory of the spectrum allocated by FFTW
      type(c_ptr) :: cmem = c_null_ptr
#endif
   contains

      !> Release the work arrays
      final :: delete_fft_rwork

   end type fft_rwork

   !> Minimal number of grid points per thread, smaller grids use fewer
   !> threads as starting the threads outweighs the transform
   integer, parameter :: points_per_thread = 32768


contains


!> Create a plan for transforms of real grids.
!>
!> Planning is not thread-safe and happens in a critical section, the
!> execution of the plans is thread-safe.
subroutine new_fft_rplan(self, n, active)

   !> Instance of the plan
   type(fft_rplan), intent(out) :: self

   !> Grid dimensions, the first one must be even
   integer, intent(in) :: n(3)

   !> Leading active extent of the grid in each dimension
   integer, intent(in), optional :: active(3)

   integer :: nthreads

   self%n = n
   self%active = n
   if (present(active)) self%active = max(1, min(active, n))

   nthreads = 1
   !$ nthreads = omp_max_threads()
   self%nthreads = max(1, min(nthreads, product(n) / points_per_thread))

#ifdef WITH_FFTW
   self%backend = "FFTW (real-to-complex)"
   self%nspec = n(1) / 2 + 1
   call plan_fftw(self)
#else
   self%backend = "in-house (complex)"
   self%nspec = n(1)
   call new_fft_plan(self%plan, n)
#endif

   allocate(self%work)
   call new_fft_rwork(self%work, self)

contains

   !> Number of OpenMP threads available to the transforms
   function omp_max_threads() result(nthreads)
      !$ use omp_lib, only : omp_get_max_threads
      !> Number of threads
      integer :: nthreads

      nthreads = 1
      !$ nthreads = omp_get_max_threads()
   end function omp_max_threads

end subroutine new_fft_rplan


!> Allocate the work arrays for the transforms of a plan
subroutine new_fft_rwork(self, plan)

   !> Work arrays of the transforms
   type(fft_rwork), intent(out) :: self

   !> Plan the work arrays are used with
   type(fft_rplan), intent(in) :: plan

#ifdef WITH_FFTW
   self%rmem = fftw_alloc_real(int(product(plan%n), c_size_t))
   self%cmem = fftw_alloc_complex(int(plan%nspec, c_size_t) &
      & * int(plan%n(2), c_size_t) * int(plan%n(3), c_size_t))
   call c_f_pointer(self%rmem, self%rgrid, plan%n)
   call c_f_pointer(self%cmem, self%cgrid, [plan%nspec, plan%n(2), plan%n(3)])
#else
   allocate(self%rgrid(plan%n(1), plan%n(2), plan%n(3)))
   allocate(self%cgrid(plan%nspec, plan%n(2), plan%n(3)))
#endif

end subroutine new_fft_rwork


!> Release the work arrays
subroutine delete_fft_rwork(self)

   !> Work arrays of the transforms
   type(fft_rwork), intent(inout) :: self

#ifdef WITH_FFTW
   if (c_associated(self%rmem)) call fftw_free(self%rmem)
   if (c_associated(self%cmem)) call fftw_free(self%cmem)
   self%rmem = c_null_ptr
   self%cmem = c_null_ptr
#else
   if (associated(self%rgrid)) deallocate(self%rgrid)
   if (associated(self%cgrid)) deallocate(self%cgrid)
#endif
   nullify(self%rgrid, self%cgrid)

end subroutine delete_fft_rwork


!> Forward transform of a real grid vanishing outside of the active block
subroutine fft_r2c_3d(self, work)

   !> Instance of the plan
   type(fft_rplan), intent(in) :: self

   !> Work arrays, the real grid is transformed into the spectrum
   type(fft_rwork), intent(inout) :: work

#ifdef WITH_FFTW
   call forward_planes(self, work%rgrid, work%cgrid)
   call transform_columns(self, self%fwd_column, work%cgrid)
#else
   call pack_real(self, work%rgrid, work%cgrid)
   call fft_3d(self%plan, work%cgrid, .false., self%active)
#endif

end subroutine fft_r2c_3d


!> Backward transform of the spectrum, only the active block of the real grid
!> is formed and the spectrum is destroyed
subroutine fft_c2r_3d(self, work)

   !> Instance of the plan
   type(fft_rplan), intent(in) :: self

   !> Work arrays, the spectrum is transformed into the real grid
   type(fft_rwork), intent(inout) :: work

#ifdef WITH_FFTW
   call transform_columns(self, self%bwd_column, work%cgrid)
   call backward_planes(self, work%cgrid, work%rgrid)
#else
   call fft_3d(self%plan, work%cgrid, .true., self%active)
   call unpack_real(self, work%cgrid, work%rgrid)
#endif

end subroutine fft_c2r_3d


#ifdef WITH_FFTW
!> Create the FFTW plans for all alignment classes of the planes and columns
subroutine plan_fftw(self)

   !> Instance of the plan
   type(fft_rplan), intent(inout) :: self

   type(fft_rwork) :: work
   integer :: ip, ia
   integer(c_int) :: n1, n2, n3, h1, a2, nn(1), flags

   n1 = int(self%n(1), c_int)
   n2 = int(self%n(2), c_int)
   n3 = int(self%n(3), c_int)
   h1 = int(self%nspec, c_int)
   a2 = int(self%active(2), c_int)

   ! Planning with FFTW_ESTIMATE leaves the arrays untouched and is cheap
   call new_fft_rwork(work, self)

   !$omp critical (multicharge_fftw_planner)
   do ia = 0, nalign
      flags = FFTW_ESTIMATE
      if (ia == nalign) flags = ior(flags, FFTW_UNALIGNED)

      ! Planes along the third dimension
      do ip = 1, int(min(n3, int(nalign, c_int)))
         if (ia < nalign) then
            if (alignment_class(c_loc(work%cgrid(1, 1, ip))) /= ia) cycle
         end if
         nn = n2
         if (.not. c_associated(self%fwd_plane%plan(ia))) then
            self%fwd_plane%plan(ia) = fftw_plan_many_dft(1_c_int, nn, h1, &
               & work%cgrid(:, :, ip), nn, h1, 1_c_int, work%cgrid(:, :, ip), nn, h1, &
               & 1_c_int, FFTW_FORWARD, flags)
            self%bwd_plane%plan(ia) = fftw_plan_many_dft(1_c_int, nn, h1, &
               & work%cgrid(:, :, ip), nn, h1, 1_c_int, work%cgrid(:, :, ip), nn, h1, &
               & 1_c_int, FFTW_BACKWARD, flags)
         end if
         exit
      end do

      ! Real planes, classified by the alignment of the real plane
      do ip = 1, int(min(n3, int(nalign, c_int)))
         if (ia < nalign) then
            if (alignment_class(c_loc(work%rgrid(1, 1, ip))) /= ia) cycle
         end if
         nn = n1
         if (.not. c_associated(self%r2c_plane%plan(ia))) then
            self%r2c_plane%plan(ia) = fftw_plan_many_dft_r2c(1_c_int, nn, a2, &
               & work%rgrid(:, :, ip), nn, 1_c_int, n1, &
               & work%cgrid(:, :, ip), nn, 1_c_int, h1, flags)
            self%c2r_plane%plan(ia) = fftw_plan_many_dft_c2r(1_c_int, nn, a2, &
               & work%cgrid(:, :, ip), nn, 1_c_int, h1, &
               & work%rgrid(:, :, ip), nn, 1_c_int, n1, flags)
         end if
         exit
      end do

      ! Column blocks along the second dimension
      do ip = 1, int(min(n2, int(nalign, c_int)))
         if (ia < nalign) then
            if (alignment_class(c_loc(work%cgrid(1, ip, 1))) /= ia) cycle
         end if
         nn = n3
         if (.not. c_associated(self%fwd_column%plan(ia))) then
            self%fwd_column%plan(ia) = plan_column(work%cgrid, ip, FFTW_FORWARD, flags)
            self%bwd_column%plan(ia) = plan_column(work%cgrid, ip, FFTW_BACKWARD, flags)
         end if
         exit
      end do
   end do
   !$omp end critical (multicharge_fftw_planner)

contains

   !> Plan the transform along the third dimension of the column block ip
   function plan_column(cgrid, ip, sign, flags) result(plan)
      !> Spectrum used for planning
      complex(wp), intent(inout), target :: cgrid(:, :, :)
      !> Index of the column block along the second dimension
      integer, intent(in) :: ip
      !> Direction of the transform
      integer(c_int), intent(in) :: sign
      !> Planner flags
      integer(c_int), intent(in) :: flags
      !> FFTW plan
      type(c_ptr) :: plan

      complex(c_double_complex), pointer :: column(:)

      call c_f_pointer(c_loc(cgrid(1, ip, 1)), column, [h1 * n2 * n3])
      plan = fftw_plan_many_dft(1_c_int, [n3], h1, column, [n3], h1 * n2, 1_c_int, &
         & column, [n3], h1 * n2, 1_c_int, sign, flags)
   end function plan_column

end subroutine plan_fftw


!> Alignment class of an address within the widest SIMD alignment
pure function alignment_class(ptr) result(ia)

   !> Address
   type(c_ptr), intent(in) :: ptr

   !> Alignment class in units of 16 byte
   integer :: ia

   ia = int(modulo(transfer(ptr, 0_c_intptr_t), int(16 * nalign, c_intptr_t)) / 16)

end function alignment_class


!> Select the plan matching the alignment of an address
function select_plan(plans, ptr) result(plan)

   !> Plans for each alignment class
   type(aligned_plans), intent(in) :: plans

   !> Address the plan is executed on
   type(c_ptr), intent(in) :: ptr

   !> Matching plan, the unaligned plan if no plan matches
   type(c_ptr) :: plan

   plan = plans%plan(alignment_class(ptr))
   if (.not. c_associated(plan)) plan = plans%plan(nalign)

end function select_plan


!> Transform the active planes along the first two dimensions
subroutine forward_planes(self, rgrid, cgrid)

   !> Instance of the plan
   type(fft_rplan), intent(in) :: self

   !> Real grid
   real(wp), intent(inout), contiguous, target :: rgrid(:, :, :)

   !> Spectrum
   complex(wp), intent(inout), contiguous, target :: cgrid(:, :, :)

   integer :: i3

   !$omp parallel do default(none) schedule(static) shared(self, rgrid, cgrid) &
   !$omp private(i3) num_threads(self%nthreads)
   do i3 = 1, self%n(3)
      if (i3 <= self%active(3)) then
         cgrid(:, self%active(2) + 1:, i3) = (0.0_wp, 0.0_wp)
         call fftw_execute_dft_r2c(select_plan(self%r2c_plane, c_loc(rgrid(1, 1, i3))), &
            & rgrid(:, :, i3), cgrid(:, :, i3))
         call fftw_execute_dft(select_plan(self%fwd_plane, c_loc(cgrid(1, 1, i3))), &
            & cgrid(:, :, i3), cgrid(:, :, i3))
      else
         cgrid(:, :, i3) = (0.0_wp, 0.0_wp)
      end if
   end do

end subroutine forward_planes


!> Back transform the active planes along the second and first dimension
subroutine backward_planes(self, cgrid, rgrid)

   !> Instance of the plan
   type(fft_rplan), intent(in) :: self

   !> Spectrum, destroyed on output
   complex(wp), intent(inout), contiguous, target :: cgrid(:, :, :)

   !> Real grid
   real(wp), intent(inout), contiguous, target :: rgrid(:, :, :)

   integer :: i3

   !$omp parallel do default(none) schedule(static) shared(self, rgrid, cgrid) &
   !$omp private(i3) num_threads(self%nthreads)
   do i3 = 1, self%active(3)
      call fftw_execute_dft(select_plan(self%bwd_plane, c_loc(cgrid(1, 1, i3))), &
         & cgrid(:, :, i3), cgrid(:, :, i3))
      call fftw_execute_dft_c2r(select_plan(self%c2r_plane, c_loc(rgrid(1, 1, i3))), &
         & cgrid(:, :, i3), rgrid(:, :, i3))
   end do

end subroutine backward_planes


!> Transform the column blocks along the third dimension
subroutine transform_columns(self, plans, cgrid)

   !> Instance of the plan
   type(fft_rplan), intent(in) :: self

   !> Plans of the column blocks in the direction of the transform
   type(aligned_plans), intent(in) :: plans

   !> Spectrum, transformed in place
   complex(wp), intent(inout), contiguous, target :: cgrid(:, :, :)

   integer :: i2
   complex(c_double_complex), pointer :: column(:)

   !$omp parallel do default(none) schedule(static) shared(self, plans, cgrid) &
   !$omp private(i2, column) num_threads(self%nthreads)
   do i2 = 1, self%n(2)
      call c_f_pointer(c_loc(cgrid(1, i2, 1)), column, &
         & [self%nspec * self%n(2) * self%n(3)])
      call fftw_execute_dft(select_plan(plans, c_loc(cgrid(1, i2, 1))), column, column)
   end do

end subroutine transform_columns
#endif


#ifndef WITH_FFTW
!> Copy the active lines of a real grid into the complex grid of the in-house
!> transform, the remaining lines vanish
subroutine pack_real(self, rgrid, cgrid)

   !> Instance of the plan
   type(fft_rplan), intent(in) :: self

   !> Real grid
   real(wp), intent(in) :: rgrid(:, :, :)

   !> Complex grid
   complex(wp), intent(out) :: cgrid(:, :, :)

   integer :: i2, i3

   !$omp parallel do default(none) collapse(2) schedule(static) &
   !$omp shared(self, rgrid, cgrid) private(i2, i3) num_threads(self%nthreads)
   do i3 = 1, self%n(3)
      do i2 = 1, self%n(2)
         if (i2 <= self%active(2) .and. i3 <= self%active(3)) then
            cgrid(:, i2, i3) = cmplx(rgrid(:, i2, i3), 0.0_wp, wp)
         else
            cgrid(:, i2, i3) = (0.0_wp, 0.0_wp)
         end if
      end do
   end do

end subroutine pack_real


!> Copy the real part of the active block of the complex grid
subroutine unpack_real(self, cgrid, rgrid)

   !> Instance of the plan
   type(fft_rplan), intent(in) :: self

   !> Complex grid
   complex(wp), intent(in) :: cgrid(:, :, :)

   !> Real grid, only the active block is formed
   real(wp), intent(inout) :: rgrid(:, :, :)

   integer :: i2, i3

   !$omp parallel do default(none) collapse(2) schedule(static) &
   !$omp shared(self, rgrid, cgrid) private(i2, i3) num_threads(self%nthreads)
   do i3 = 1, self%active(3)
      do i2 = 1, self%active(2)
         rgrid(:, i2, i3) = real(cgrid(:, i2, i3), wp)
      end do
   end do

end subroutine unpack_real
#endif


!> Release the plans and the work arrays
subroutine delete_fft_rplan(self)

   !> Instance of the plan
   type(fft_rplan), intent(inout) :: self

   if (associated(self%work)) deallocate(self%work)
   self%work => null()

#ifdef WITH_FFTW
   !$omp critical (multicharge_fftw_planner)
   call destroy_plans(self%r2c_plane)
   call destroy_plans(self%c2r_plane)
   call destroy_plans(self%fwd_plane)
   call destroy_plans(self%bwd_plane)
   call destroy_plans(self%fwd_column)
   call destroy_plans(self%bwd_column)
   !$omp end critical (multicharge_fftw_planner)

contains

   !> Destroy the plans of all alignment classes
   subroutine destroy_plans(plans)
      !> Plans for each alignment class
      type(aligned_plans), intent(inout) :: plans

      integer :: ia

      do ia = 0, nalign
         if (c_associated(plans%plan(ia))) call fftw_destroy_plan(plans%plan(ia))
         plans%plan(ia) = c_null_ptr
      end do
   end subroutine destroy_plans
#endif

end subroutine delete_fft_rplan


end module multicharge_fft_real
