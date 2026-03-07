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
    use iso_fortran_env, only: output_unit
    use mctc_env, only: error_type, fatal_error, wp
    use multicharge_blas, only: symv
    use multicharge_lapack, only: sytrf, sytrs, sytri
    use multicharge_solver_type, only: mchrg_solver_type, mchrg_solver_input
    implicit none
    private

    public :: direct_solver, direct_input, new_direct_solver

    !> Input for Direct solver
    type, extends(mchrg_solver_input) :: direct_input
        ! Use Direct solver
        logical :: direct = .true.
        ! Verbosity
        integer, allocatable :: verbosity
    end type direct_input

    !> Direct solver using LAPACK
    type, extends(mchrg_solver_type) :: direct_solver
        integer, allocatable :: verbosity
    contains
        procedure :: solve
    end type direct_solver

    !> Default verbosity level
    integer, parameter :: verbosity_def = 0

contains

    subroutine new_direct_solver(self, input)
        class(direct_solver), intent(out) :: self
        type(direct_input), intent(in) :: input 
        
        self%need_pos_def = .false.
        if (allocated(input%verbosity)) then
            self%verbosity = input%verbosity
        else 
            self%verbosity = verbosity_def
        end if

    end subroutine new_direct_solver

    !> Solve method for direct solver
    subroutine solve(self, amat, xvec, vrhs, ainv, cpq, new_unit, error)
        class(direct_solver), intent(in) :: self
        !> A matrix of Ax=b system
        real(wp), intent(in)  :: amat(:, :)
        !> Initial search direction (b)
        real(wp), intent(in)  :: xvec(:)
        !> Initial guess and solution
        real(wp), intent(inout) :: vrhs(:)
        !> Inverse A-matrix and coupled perturbed logical
        real(wp), intent(out), optional :: ainv(:, :)
        !> Coupled-perturbed equations flag (optional)
        logical, intent(in), optional :: cpq
        !> Output unit (optional)
        integer, intent(in), optional :: new_unit
        !> Error handling
        type(error_type), allocatable, intent(out) :: error
    
        ! Local inverse matrix required for calculations
        real(wp), allocatable :: invmat(:,:)
        integer  :: local_info
        integer :: ndim, ic, jc
        integer, allocatable :: ipiv(:)
        logical :: want_cpq
        integer :: unit

        if (present(new_unit)) then
            unit = new_unit
        else
            unit = output_unit
        end if

        if (self%verbosity > 0) then
            call write_direct_solver(unit, self)
        end if 
    
        ! Dimensions match check
        ndim = size(xvec)
        if (size(amat,1) /= ndim .or. size(amat,2) /= ndim) then
            call fatal_error(error, "solve_direct: dimension mismatch.")
            return
        end if

        allocate(invmat(ndim, ndim))
        invmat = amat
        vrhs = xvec
    
        ! Logical: solve coupled-perturbed equations flag
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

        if (present(ainv)) ainv=invmat
    
    end subroutine solve

subroutine write_direct_solver(unit, solver)
   integer, intent(in) :: unit
   class(direct_solver), intent(in) :: solver

   write(unit, '(a)') "Using Direct Solver"
   write(unit, '(a)')

end subroutine write_direct_solver

end module multicharge_solver_direct