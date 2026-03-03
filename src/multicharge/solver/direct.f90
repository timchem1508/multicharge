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

module direct_solver
    use iso_fortran_env, only: output_unit
    use mctc_env, only: error_type, fatal_error, wp
    use multicharge_blas, only: symv
    use multicharge_lapack, only: sytrf, sytrs, sytri
    use solver_type, only: mchrg_solver_type, mchrg_solver_input
    use solver_cache, only: cache_container, mchrg_solver_cache

    implicit none
    private

    public :: mchrg_solver_direct, direct_input, new_direct_solver
    type, extends(mchrg_solver_cache), public :: direct_cache
    end type direct_cache

    !> Input for Direct solver
    type, extends(mchrg_solver_input) :: direct_input
        ! Use Direct solver
        logical :: direct = .true.
    end type direct_input

    !> Direct solver using LAPACK
    type, extends(mchrg_solver_type) :: mchrg_solver_direct
    contains
       procedure :: solve
       procedure :: update 
    end type mchrg_solver_direct

contains

    subroutine new_direct_solver(self, input)
        class(mchrg_solver_type), intent(out) :: self
        type(direct_input), intent(in) :: input 

        self%need_pos_def = .false.

    end subroutine new_direct_solver

    !> Update method for direct solver
    subroutine update(self, cache, vrhs, ainv, cpq)
        class(mchrg_solver_direct), intent(in) :: self
        type(cache_container), intent(inout) :: cache
        real(wp), intent(inout) :: vrhs(:)
        real(wp), intent(out) :: ainv(:, :)
        logical, intent(in), optional :: cpq
        
    end subroutine update

    !> Solve method for direct solver
    subroutine solve(self, amat, xvec, vrhs, ainv, cpq, new_unit, error)
        class(mchrg_solver_direct), intent(in) :: self
        ! A matrix of Ax=b system
        real(wp), intent(in)  :: amat(:, :)
        ! Initial search direction (b)
        real(wp), intent(in)  :: xvec(:)
        ! Initial guess and solution
        real(wp), intent(inout) :: vrhs(:)
        ! Inverse A-matrix and coupled perturbed logical
        real(wp), intent(out) :: ainv(:, :)
        ! Coupled-perturbed equations flag (optional)
        logical, intent(in), optional :: cpq
        ! Output unit (optional)
        integer, intent(in), optional :: new_unit
        !> Error handling
        type(error_type), allocatable, intent(out) :: error
    
        integer  :: local_info
        integer :: ndim, ic, jc
        integer, allocatable :: ipiv(:)
        logical :: want_cpq
        integer :: unit
        type(cache_container), allocatable :: cache

        if (present(new_unit)) then
            unit = new_unit
        else
            unit = output_unit
        end if
    
        ! Dimensions match check
        ndim = size(xvec)
        if (size(amat,1) /= ndim .or. size(amat,2) /= ndim) then
            call fatal_error(error, "solve_direct: dimension mismatch.")
            return
        end if

        ainv = amat
        vrhs = xvec
        
        ! Update cache and prepare vrhs and ainv
        allocate(cache)
        call self%update(cache, vrhs, ainv, cpq)
    
        ! Logical: solve coupled-perturbed equations flag
        want_cpq = .false.
        if (present(cpq)) want_cpq = cpq
    
        allocate(ipiv(ndim))
        call sytrf(ainv, ipiv, info=local_info, uplo='l')
        if (local_info /= 0) then
            call fatal_error(error, "solve_direct: Bunch-Kaufman factorization failed.")
            return
        end if
    
        if (want_cpq) then
            call sytri(ainv, ipiv, info=local_info, uplo='l')
            if (local_info /= 0) then
                call fatal_error(error, "solve_direct: Inversion of factorized matrix failed.")
                return
            end if
            call symv(ainv, xvec, vrhs, uplo='l')
            do ic = 1, ndim
                do jc = ic + 1, ndim
                    ainv(ic, jc) = ainv(jc, ic)
                end do
            end do
        else
            call sytrs(ainv, vrhs, ipiv, info=local_info, uplo='l')
            if (local_info /= 0) then
                call fatal_error(error, "solve_direct: Solving factorized system failed.")
                return
            end if
        end if
    
    end subroutine solve

end module direct_solver