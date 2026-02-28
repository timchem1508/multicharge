!> @dir multicharge/solver
!> Contains the implementation of the linear equations solvers

!> @file multicharge/solver.f90
!> Provides a reexport of the solvers implementations

!> Proxy module to reexport the solver implementations
module multicharge_solver
   use mctc_env, only : error_type, fatal_error
   use solver_type, only : mchrg_solver_type, mchrg_solver_input
   use direct_solver, only : mchrg_solver_direct, new_direct_solver, & 
                                        & direct_input, direct_cache
   use cg_solver, only : mchrg_solver_cg, new_cg_solver, cg_input, &
                                    & cg_cache
   use solver_cache, only: mchrg_solver_cache
   implicit none
   private

   public :: mchrg_solver_type, mchrg_solver_input
   public :: mchrg_solver_direct, new_direct_solver, direct_input
   public :: mchrg_solver_cg, new_cg_solver, cg_input
   public :: new_mchrg_solver

contains 

subroutine new_mchrg_solver(solver, input, error)
    class(mchrg_solver_type), intent(out), allocatable :: solver
    class(mchrg_solver_input), intent(in) :: input
    type(error_type), allocatable, intent(out) :: error

    select type (input)
    type is (cg_input)
        block
            class(mchrg_solver_cg), allocatable :: tmp
            allocate(tmp)
            call new_cg_solver(tmp, input)
            call move_alloc(tmp, solver)
        end block
    type is (direct_input)
        block
            class(mchrg_solver_direct), allocatable :: tmp
            allocate(tmp)
            call new_direct_solver(tmp, input)
            call move_alloc(tmp, solver)
        end block
    class default 
        call fatal_error(error, "multicharge/solver.f90: Unknown solver input type")
        return
    end select
    
end subroutine new_mchrg_solver

end module multicharge_solver