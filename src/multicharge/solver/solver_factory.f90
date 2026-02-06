module solver_factory
    use mctc_env, only: wp
    use solver, only: mchrg_solver_type
    use direct_solver, only: direct_solver_type
    use iterative_solver, only: cg_solver_type, cg_input
    use mchrg_solver_input, only: solver_input
    implicit none
    private

    public :: new_mchrg_solver

contains

    !> Create a new solver instance
    !> solver_type String identifying the solver ("DIRECT" or "CG"/"ITERATIVE")
    !> cgmiter Optional maximum iterations (for iterative solver)
    !> cgtol Optional tolerance (for iterative solver)
    !> cgmode Mode of the iterative solver (for the following benchmarks)
    subroutine new_mchrg_solver(solver_type, cgmiter, cgtol, cgmode, solver)
        character(len=*), intent(in) :: solver_type
        integer, intent(in), optional :: cgmiter
        real(wp), intent(in), optional :: cgtol
        character(len=32), optional :: cgmode
        class(mchrg_solver_type), intent(out), allocatable :: solver

        class(solver_input), allocatable :: sinput
        logical, allocatable :: cg

        character(len=10) :: type_upper
        
        type_upper = trim(solver_type)
        ! Simple uppercase conversion (assuming ASCII)
        ! In a robust app, use a dedicated utility
        
        if (type_upper == "DIRECT" .or. type_upper == "direct") then
           allocate(direct_solver_type :: solver)
        else
            ! Default to CG
            allocate(cg_solver_type :: solver)
            cg = .true.
            allocate(sinput)
            sinput%cg = cg_input(cgmiter, cgtol, cgmode)
            ! Apply settings if allocated as CG
            select type(solver)
            type is (cg_solver_type)
                solver%cgmiter = sinput%cg%cgmiter
                solver%cgtol = sinput%cg%cgtol
                solver%cgmode = sinput%cg%cgmode
            end select

        end if
    end subroutine new_mchrg_solver

end module solver_factory