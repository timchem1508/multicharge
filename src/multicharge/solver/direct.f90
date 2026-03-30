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
    use mctc_env, only: error_type, fatal_error, wp, timer_type, timer_type, format_time
    use mctc_ncoord, only: adjacency_list
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

    !> New direct solver construction
    subroutine new_direct_solver(self, input)
        !> Direct solver type
        class(direct_solver), intent(out) :: self
        !> Direct input type
        type(direct_input), intent(in) :: input 
        
        self%need_pos_def = .false.
        if (allocated(input%verbosity)) then
            self%verbosity = input%verbosity
        else 
            self%verbosity = verbosity_def
        end if

    end subroutine new_direct_solver

    !> Solve method for direct solver
    subroutine solve(self, amat, alist, adiag, xvec, vrhs, ainv, cpq, list, new_unit, error)
        class(direct_solver), intent(in) :: self
        !> A matrix of Ax=b system
        real(wp), intent(in), optional  :: amat(:, :)
        !> Off-diagonall elements of matrix for in compressed
        real(wp), intent(in), optional  :: alist(:)
        !> Diagonall elements of tmatrix for in compressed
        real(wp), intent(in), optional  :: adiag(:)
        !> Right-hand side vector 
        real(wp), intent(in)  :: xvec(:)
        !> On input: initial guess; on output: solution
        real(wp), intent(inout), contiguous :: vrhs(:)
        !> Inverse matrix – not computed by CG, but required by interface
        real(wp), intent(out), optional :: ainv(:, :)
        !> Flag for coupled-perturbed equations
        logical, intent(in), optional :: cpq
        !> Neighbour list optional type
        type(adjacency_list), intent(in), optional :: list
        !> Output unit (optional)
        integer, intent(in), optional :: new_unit
        !> Error handling
        type(error_type), allocatable, intent(out) :: error
    
        real(wp), allocatable :: invmat(:,:)
        integer, allocatable :: ipiv(:)
        integer  :: local_info
        integer :: ndim, ic, jc
        logical :: want_cpq
        integer :: unit
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

    subroutine write_direct_solver(unit)
    integer, intent(in) :: unit

    write(unit, '(a)') "Using Direct Solver"
    write(unit, '(a)')

    end subroutine write_direct_solver

    !> Print final summary
    subroutine print_direct_final(unit, timer, verbosity)
        integer, intent(in) :: unit, verbosity
        type(timer_type), intent(in) :: timer

        if (verbosity > 1) then
            write(unit, '(a, 1x, a)') "Direct solver time : ", format_time(timer%get("total"))
            write(unit, '(a)') ''
        end if
    end subroutine print_direct_final

end module multicharge_solver_direct