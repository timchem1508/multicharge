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

!> @file multicharge/solver/direct.f90
!> Provides implementation of the direct solver using LAPACK for symmetric indefinite systems.

module multicharge_solver_direct
   use iso_fortran_env, only : output_unit
   use mctc_env, only : error_type, fatal_error, format_time, timer_type, wp
   use mctc_csrlist, only : csr_list
   use multicharge_blas, only : symv
   use multicharge_lapack, only : sytrf, sytri, sytrs
   use multicharge_solver_type, only : mchrg_solver_input, mchrg_solver_type
   implicit none
   private

   public :: direct_solver, direct_input, new_direct_solver

   !> Input configuration for the direct solver
   type, extends(mchrg_solver_input) :: direct_input
      !> Whether to use the direct solver
      logical :: direct = .true.

      !> Optional output verbosity
      integer, allocatable :: verbosity
   end type direct_input

   !> Direct LAPACK solver for symmetric indefinite systems
   type, extends(mchrg_solver_type) :: direct_solver
      !> Output verbosity
      integer, allocatable :: verbosity
   contains
      !> Solve the linear system directly
      procedure :: solve
   end type direct_solver

   !> Default verbosity level
   integer, parameter :: verbosity_def = 0

contains

   !> Construct a direct solver from its input configuration
   subroutine new_direct_solver(self, input)
      !> Direct solver instance
      class(direct_solver), intent(out) :: self

      !> Direct solver configuration
      type(direct_input), intent(in) :: input

      self%need_pos_def = .false.
      if (allocated(input%verbosity)) then
         self%verbosity = input%verbosity
      else
         self%verbosity = verbosity_def
      end if

   end subroutine new_direct_solver

   !> Solve a dense symmetric linear system using LAPACK
   subroutine solve(self, amat, alist, xvec, vrhs, ainv, cpq, list, new_unit, error)
      !> Direct solver instance
      class(direct_solver), intent(in) :: self

      !> Dense coefficient matrix of the linear system
      real(wp), intent(in), optional :: amat(:, :)

      !> Coefficient matrix values in compressed-row storage
      real(wp), intent(in), optional :: alist(:)

      !> Right-hand side vector
      real(wp), intent(in) :: xvec(:)

      !> On input: initial guess; on output: solution
      real(wp), intent(inout), contiguous :: vrhs(:)

      !> Inverse coefficient matrix
      real(wp), intent(out), optional :: ainv(:, :)

      !> Whether to solve coupled-perturbed equations
      logical, intent(in), optional :: cpq

      !> Optional neighbour-list representation of the matrix
      type(csr_list), intent(in), optional :: list

      !> Output unit (optional)
      integer, intent(in), optional :: new_unit

      !> Error handling
      type(error_type), allocatable, intent(out) :: error

      real(wp), allocatable :: invmat(:, :)
      integer, allocatable :: ipiv(:)

      integer :: local_info
      integer :: ndim, ic, jc

      integer :: unit

      logical :: want_cpq

      type(timer_type) :: timer

      if (self%verbosity > 1) call timer%push("total")

      if (present(new_unit)) then
         unit = new_unit
      else
         unit = output_unit
      end if

      if (self%verbosity > 0) then
         call write_direct_solver(unit)
      end if

      ! Dimensions match check
      ndim = size(xvec)

      if (present(amat)) then
         if (size(amat,1) /= ndim .or. size(amat,2) /= ndim) then
            call fatal_error(error, "solve_direct: dimension mismatch.")
            return
         end if

         allocate(invmat(ndim, ndim))
         invmat = amat
         vrhs = xvec

         want_cpq = .false.
         if (present(cpq)) want_cpq = cpq

         ! Factorize the Coulomb matrix
         allocate(ipiv(ndim))
         call sytrf(invmat, ipiv, info=local_info, uplo='l')
         if (local_info /= 0) then
            call fatal_error(error, "Bunch-Kaufman factorization failed.")
            return
         end if

         if (want_cpq) then
            ! Inverted matrix is needed for coupled-perturbed equations
            call sytri(invmat, ipiv, info=local_info, uplo='l')
            if (local_info /= 0) then
               call fatal_error(error, "Inversion of factorized matrix failed.")
               return
            end if
            ! Solve the linear system
            call symv(invmat, xvec, vrhs, uplo='l')
            do ic = 1, ndim
               do jc = ic + 1, ndim
                  invmat(ic, jc) = invmat(jc, ic)
               end do
            end do
         else
            ! Solve the linear system
            call sytrs(invmat, vrhs, ipiv, info=local_info, uplo='l')
            if (local_info /= 0) then
               call fatal_error(error, "Solution of linear system failed.")
               return
            end if

         end if

         if (present(ainv)) ainv = invmat

      end if

      ! pop solve timer
      call timer%pop
      call print_direct_final(unit, timer, self%verbosity)

   end subroutine solve

!> Print the direct solver banner
   subroutine write_direct_solver(unit)
      !> Output unit
      integer, intent(in) :: unit

      write(unit, '(a)') "Using Direct Solver"
      write(unit, '(a)')

   end subroutine write_direct_solver

   !> Print final summary
   subroutine print_direct_final(unit, timer, verbosity)
      !> Output unit
      integer, intent(in) :: unit, verbosity

      !> Timer holding the accumulated execution time
      type(timer_type), intent(in) :: timer

      if (verbosity > 1) then
         write(unit, '(a, 1x, a)') "Direct solver time : ", format_time(timer%get("total"))
         write(unit, '(a)') ''
      end if
   end subroutine print_direct_final

end module multicharge_solver_direct
