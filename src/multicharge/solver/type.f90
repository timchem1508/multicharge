module solver
    use mctc_env, only: error_type, wp
    use solver_type_cache, only: cache_container
    implicit none
    private

    public :: mchrg_solver_type

    !> Abstract base type for Multi-charge solvers
    type, abstract :: mchrg_solver_type
    contains
       procedure(solve), deferred :: solve
       procedure(update), deferred :: update
    end type mchrg_solver_type

    abstract interface
        subroutine solve(self, amat, xvec, vrhs, ainv, cpq, error)
            import :: mchrg_solver_type, error_type, wp
            class(mchrg_solver_type), intent(in) :: self
            real(wp), intent(in)  :: amat(:, :)
            real(wp), intent(in)  :: xvec(:)
            real(wp), intent(inout) :: vrhs(:)
            real(wp), intent(out) :: ainv(:, :)
            logical, intent(in), optional :: cpq
            type(error_type), allocatable, intent(out) :: error
        end subroutine solve

        subroutine update(self, cache, vrhs, ainv, cpq)
            import :: mchrg_solver_type, cache_container, wp
            class(mchrg_solver_type), intent(in) :: self
            type(cache_container), intent(inout) :: cache
            real(wp), intent(inout) :: vrhs(:)
            real(wp), intent(out) :: ainv(:, :)
            logical, intent(in), optional :: cpq
        end subroutine update
    end interface

end module solver