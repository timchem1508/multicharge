# 1 "/home/zakharov/Documents/prog/mcharge/multicharge/src/multicharge/solver/type.f90"
# 1 "<built-in>"
# 1 "<command-line>"
# 1 "/home/zakharov/Documents/prog/mcharge/multicharge/cmake_build//"
# 1 "/home/zakharov/Documents/prog/mcharge/multicharge/src/multicharge/solver/type.f90"
module multicharge_solver_type
    use mctc_env, only: error_type, wp
    use multicharge_solver_cache, only: cache_container
    implicit none
    private

    public :: mchrg_solver_type, mchrg_solver_input

    !> Abstract base type for Multi-charge solvers
    type, abstract :: mchrg_solver_type
        !> Type of matrix availiable for solver
        !> CG solver can use only positive definite matrices
        !> Direct solver can use either type of matrix
        logical, allocatable :: need_pos_def
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

    !> Solver input abstract type
    type, abstract, public :: mchrg_solver_input
    end type mchrg_solver_input

end module multicharge_solver_type
