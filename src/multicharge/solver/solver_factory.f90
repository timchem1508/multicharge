module solver_factory
    use solver, only: mchrg_solver_type
    use direct_solver, only: direct_solver_type
    use iterative_solver, only: cg_solver_type
    implicit none
    private

    public :: new_mchrg_solver

contains

    function new_mchrg_solver(use_cg) result(solver)
        logical, intent(in), optional :: use_cg
        class(mchrg_solver_type), allocatable :: solver
        logical :: cg
        character(len=32) :: env
    
        cg = .true.
        if (present(use_cg)) then
           cg = use_cg
        else
           call get_environment_variable("MCHARGE_SOLVER", env)
           if (trim(env) == "DIRECT") cg = .false.
        end if
    
        if (cg) then
           allocate(cg_solver_type :: solver)
        else
           allocate(direct_solver_type :: solver)
        end if
     end function new_mchrg_solver

end module solver_factory