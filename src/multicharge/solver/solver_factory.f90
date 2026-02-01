module solver_factory
    use mctc_env, only: wp
    use solver, only: mchrg_solver_type
    use direct_solver, only: direct_solver_type
    use iterative_solver, only: cg_solver_type
    implicit none
    private

    public :: new_mchrg_solver

contains

    !> Create a new solver instance
    !> @param solver_type String identifying the solver ("DIRECT" or "CG"/"ITERATIVE")
    !> @param cgmiter Optional maximum iterations (for iterative solver)
    !> @param cgtol Optional tolerance (for iterative solver)
    !> @param cgmode Mode of the iterative solver (for the following benchmarks)
    function new_mchrg_solver(solver_type, cgmiter, cgtol, cgmode) result(solver)
        character(len=*), intent(in) :: solver_type
        integer, intent(in) :: cgmiter
        real(wp), intent(in) :: cgtol
        character(len=32), intent(in) :: cgmode
        class(mchrg_solver_type), allocatable :: solver
        
        character(len=10) :: type_upper

        write(*,*) "Solver type:", solver_type
        
        type_upper = trim(solver_type)
        ! Simple uppercase conversion (assuming ASCII)
        ! In a robust app, use a dedicated utility
        
        if (type_upper == "DIRECT" .or. type_upper == "direct") then
           allocate(direct_solver_type :: solver)
        else
           ! Default to CG
           allocate(cg_solver_type :: solver)
           
           ! Apply settings if allocated as CG
           select type(slv => solver)
           type is (cg_solver_type)
               slv%cgmiter = cgmiter
               !write(*,*) "CG max iterations:", cgmiter
               slv%cgtol = cgtol
               !write(*,*) "CG tolerance", cgtol
               slv%cgmode = cgmode
               !write(*,*) "CG mode", cgmode
           end select
        end if
     end function new_mchrg_solver

end module solver_factory