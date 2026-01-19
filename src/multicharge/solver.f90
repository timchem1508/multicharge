#ifndef IK
#define IK i4
#endif


module solver
    use mctc_env, only: wp, ik => IK, fatal_error
    use multicharge_blas, only: symv, gemv
    use multicharge_lapack, only: sytrf, sytrs, sytri
    use multicharge_model_cache, only: model_cache, cache_container
    use multicharge_model_type, only: local_charge 
    implicit none
    private

    public :: mchrg_solver_type, new_mchrg_solver
 
    type, abstract :: mchrg_solver_type
    contains
       procedure(solve_if), deferred :: solve
       procedure(update), deferred :: update
    end type mchrg_solver_type
 

    abstract interface
        subroutine update(self, cache, amat, xvec, vrhs, ainv, cpq)
            import :: mchrg_model_type, structure_type, cache_container, wp
            class(mchrg_solver_type), intent(in) :: self
            type(cache_container), intent(inout) :: cache
            real(wp), intent(in)  :: amat(:, :)
            real(wp), intent(in)  :: xvec(:)
            real(wp), intent(out) :: vrhs(:)
            real(wp), intent(out), optional :: ainv(:, :)    ! <-- made optional for consistency
            logical, intent(in), optional :: cpq
            integer(ik), intent(out), optional :: info
        end subroutine update   
    end interface

contains

    ! direct solver concrete type (was referenced but missing)
    type, extends(mchrg_solver_type) :: direct_solver_type
    contains
       procedure :: solve => solve_direct
       procedure :: update   ! note: concrete update must be provided elsewhere
    end type direct_solver_type

subroutine solve_direct(self, amat, xvec, vrhs, ainv, cpq, error, info)
       class(mchrg_solver_type), intent(in) :: self
       type(error_type), allocatable, intent(out) :: error
       real(wp), intent(in)  :: amat(:, :)
       real(wp), intent(in)  :: xvec(:)
       real(wp), intent(out) :: vrhs(:)
       real(wp), intent(out) :: ainv(:, :)
       logical, intent(in), optional :: cpq
       integer(ik), intent(out), optional :: info

       integer(ik) :: local_info
       integer :: ndim, ic, jc
       integer(ik), allocatable :: ipiv(:)
       logical :: want_cpq
       type(cache_container), allocatable :: cache

       ! Dimensions match check (do this before calling update)
       ndim = size(xvec)
       if (size(amat,1) /= ndim .or. size(amat,2) /= ndim .or. size(vrhs) /= ndim) then
          call fatal_error(local_info, "solve_direct: dimension mismatch.")
          if (present(info)) info = -1_ik
          return
       end if

       vrhs = xvec
       ainv = amat

       ! Update cache
       allocate(cache)
       call self%update(cache, amat, xvec, vrhs, ainv, cpq, info)    ! <-- call bound procedure without explicit self

       ! Logical: compute inverse flag
       want_cpq = .false.
       if (present(cpq)) want_cpq = cpq

       call sytrf(ainv, ipiv, info=local_info, uplo='l')
       if (local_info /= 0) then
          call fatal_error(error, "solve_direct: Bunch-Kaufman factorization failed.")
          if (present(info)) info = local_info
          return
       end if

       if (want_cpq) then
          call sytri(ainv, ipiv, info=local_info, uplo='l')
          if (local_info /= 0) then
             call fatal_error(error, "solve_direct: Inversion of factorized matrix failed.")
             if (present(info)) info = local_info
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
             if (present(info)) info = local_info
             return
          end if
       end if

       if (present(info)) info = local_info
end subroutine solve_direct

    ! CG solver with Jacobi preconditioner
    type, extends(mchrg_solver_type) :: cg_solver_type
    contains
       procedure :: solve => solve_cg
    end type cg_solver_type

subroutine solve_cg(self, amat, xvec, vrhs, ainv, cpq, error, info)
       class(mchrg_solver_type), intent(in) :: self
       type(error_type), allocatable, intent(out) :: error
       real(wp), intent(in)  :: amat(:, :)
       real(wp), intent(in)  :: xvec(:)
       real(wp), intent(out) :: vrhs(:)
       real(wp), intent(out), optional :: ainv(:, :)
       logical, intent(in), optional :: cpq
       integer(ik), intent(out), optional :: info
       
       integer :: ndim, it, maxit
       real(wp) :: tol, bnorm, rnorm, alpha, beta, denom
       real(wp), allocatable :: r(:), p(:), z(:), Ap(:), Mdiag(:)
       integer(ik) :: local_info
       real(wp) :: rz_old, rz_new
       type(cache_container), allocatable :: cache
       
       ! Dimensions match check (do this before calling update)
       ndim = size(xvec)
       if (size(amat,1) /= ndim .or. size(amat,2) /= ndim .or. size(vrhs) /= ndim) then
          call fatal_error(local_info, "solve_direct: dimension mismatch.")
          if (present(info)) info = -1_ik
          return
       end if

       ! Prepare/cache and allow update to modify vrhs/ainv
       allocate(cache)
       call self%update(cache, amat, xvec, vrhs, ainv, cpq, info)   ! <-- no explicit self
 
       ! Global thresholds
       tol = 1.0e-8_wp
       maxit = max(10, ndim*10)

       allocate(r(ndim), p(ndim), z(ndim), Ap(ndim), Mdiag(ndim))


         ! Jacobi preconditioner (inverse of diagonal)
       do it = 1, ndim
           Mdiag(it)=amat(it,it)
           if (abs(Mdiag(it)) < tol**3) Mdiag(it) = tol**3
           Mdiag(it) = 1.0_wp / Mdiag(it)
       end do

       vrhs = 0.0_wp  ! initial guess zero, later the local charge vector will be added

         ! Initial residual r = b - A*x (x=0)
       call gemv(amat, vrhs, Ap, alpha=1.0_wp, beta=0.0_wp, trans='n')
       ! Residual compute
       r = xvec - Ap
       ! Apply preconditioner z = M * r
       z = r * Mdiag
       ! Initial search direction
       p = z
       ! Initial direction udate factor
       bnorm = sqrt(sum(xvec*xvec))
         if (bnorm < tol**3) bnorm = 1.0_wp
         rnorm = sqrt(sum(r*r))
         if (rnorm / bnorm <= tol) then
            if (present(info)) info = 0_ik
            return
         end if

       ! Dynamical residual
       rz_old = sum(r*z)
       local_info = -1_ik

       ! Conjugate Gradient iterations
       do it = 1, maxit
           call gemv(amat, p, Ap, alpha=1.0_wp, beta=0.0_wp, trans='n')
           ! Compute step size alpha
           denom = sum(p * Ap)
              if (abs(denom) < tol**4) then
                 local_info = 0_ik
                 exit
              end if
           alpha = rz_old / denom
           ! Update solution and residual
           vrhs = vrhs + alpha * p
           r = r - alpha * Ap
           ! Check convergence
           rnorm = sqrt(sum(r*r))
              if (rnorm / bnorm <= tol) then
                 local_info = 0_ik
                 exit
              end if
           ! Apply preconditioner z = M * r
           z = r * Mdiag
           rz_new = sum(r * z)
           beta = rz_new / rz_old
           p = z + beta * p
           rz_old = rz_new
           if (it == maxit) then
              local_info = 1_ik  ! did not converge
              call fatal_error(error, "solve_cg: CG did not converge within max iterations.")
           end if
       end do

       if (present(info)) info = local_info
end subroutine solve_cg

function new_mchrg_solver(use_cg) result(solver)
    logical, intent(in), optional :: use_cg
    class(mchrg_solver_type), allocatable :: solver
    logical :: cg
    character(len=32) :: env

    cg = .false.
    if (present(use_cg)) then
       cg = use_cg
    else
       call get_environment_variable("MCHARGE_SOLVER", env)
       if (trim(env) == "CG") cg = .true.
    end if

    if (cg) then
       allocate(cg_solver_type :: solver)
    else
       allocate(direct_solver_type :: solver)
    end if
 end function new_mchrg_solver

end module solver