! Abstract solver type with direct and CG (Jacobi preconditioned) implementations
module solver
    use mctc_env, only: wp, ik => IK, fatal_error
    use multicharge_blas, only: symv, gemv
    use multicharge_lapack, only: sytrf, sytrs, sytri
    implicit none
    private
    public :: mchrg_solver_type, new_mchrg_solver
 
    type, abstract :: mchrg_solver_type
    contains
       procedure(solve_if), deferred :: solve
    end type mchrg_solver_type
 
    abstract interface
       subroutine solve_if(self, amat, xvec, vrhs, compute_inverse, ainv, info)
          import :: mchrg_solver_type, wp, ik
          class(mchrg_solver_type), intent(in) :: self
          real(wp), intent(in)  :: amat(:, :)
          real(wp), intent(in)  :: xvec(:)
          real(wp), intent(out) :: vrhs(:)
          logical, intent(in), optional :: compute_inverse
          real(wp), intent(out), optional :: ainv(:, :)
          integer(ik), intent(out), optional :: info
       end subroutine solve_if
    end interface
 
 contains
 
    ! Concrete direct solver
    type, extends(mchrg_solver_type) :: direct_solver_type
    contains
       procedure :: solve => solve_direct
    end type direct_solver_type
 
    subroutine solve_direct(self, amat, xvec, vrhs, compute_inverse, ainv, info)
       class(direct_solver_type), intent(in) :: self
       real(wp), intent(in)  :: amat(:, :)
       real(wp), intent(in)  :: xvec(:)
       real(wp), intent(out) :: vrhs(:)
       logical, intent(in), optional :: compute_inverse
       real(wp), intent(out), optional :: ainv(:, :)
       integer(ik), intent(out), optional :: info
 
       integer(ik) :: local_info
       integer :: n, i, j
       integer(ik), allocatable :: ipiv(:)
       real(wp), allocatable :: a(:,:)
       logical :: want_inv
 
       n = size(xvec)
       if (size(amat,1) /= n .or. size(amat,2) /= n .or. size(vrhs) /= n) then
          call fatal_error(local_info, "solve_direct: dimension mismatch.")
          if (present(info)) info = -1_ik
          return
       end if
 
       want_inv = .false.
       if (present(compute_inverse)) want_inv = compute_inverse
 
       allocate(a(n,n))
       a = amat
       allocate(ipiv(n))
 
       call sytrf(a, ipiv, info=local_info, uplo='l')
       if (local_info /= 0) then
          call fatal_error(local_info, "solve_direct: sytrf failed.")
          if (present(info)) info = local_info
          return
       end if
 
       vrhs = xvec
 
       if (want_inv) then
          call sytri(a, ipiv, info=local_info, uplo='l')
          if (local_info /= 0) then
             call fatal_error(local_info, "solve_direct: sytri failed.")
             if (present(info)) info = local_info
             return
          end if
          call symv(a, xvec, vrhs, uplo='l')
          ! mirror lower to upper if ainv requested
          if (present(ainv)) then
             ainv = a
             do i = 1, n
                do j = i+1, n
                   ainv(i,j) = ainv(j,i)
                end do
             end do
          end if
       else
          call sytrs(a, vrhs, ipiv, info=local_info, uplo='l')
          if (local_info /= 0) then
             call fatal_error(local_info, "solve_direct: sytrs failed.")
             if (present(info)) info = local_info
             return
          end if
       end if
 
       if (present(info)) info = local_info
    end subroutine solve_direct
 
    ! Concrete CG solver with Jacobi preconditioner
    type, extends(mchrg_solver_type) :: cg_solver_type
    contains
       procedure :: solve => solve_cg
    end type cg_solver_type
 
    subroutine solve_cg(self, amat, xvec, vrhs, compute_inverse, ainv, info)
       class(cg_solver_type), intent(in) :: self
       real(wp), intent(in)  :: amat(:, :)
       real(wp), intent(in)  :: xvec(:)
       real(wp), intent(out) :: vrhs(:)
       logical, intent(in), optional :: compute_inverse  ! ignored for CG
       real(wp), intent(out), optional :: ainv(:, :)     ! not produced
       integer(ik), intent(out), optional :: info
 
       integer :: n, it, maxit
       real(wp) :: tol, bnorm, rnorm, alpha, beta, denom
       real(wp), allocatable :: r(:), p(:), z(:), Ap(:), Mdiag(:)
       integer(ik) :: local_info
 
       n = size(xvec)
       if (size(amat,1) /= n .or. size(amat,2) /= n .or. size(vrhs) /= n) then
          call fatal_error(local_info, "solve_cg: dimension mismatch.")
          if (present(info)) info = -1_ik
          return
       end if
 
       tol = 1.0e-8_wp
       maxit = max(10, n*10)
 
       allocate(r(n), p(n), z(n), Ap(n), Mdiag(n))
 
       ! Jacobi preconditioner: inverse diagonal
       do it = 1, n
          Mdiag(it) = amat(it,it)
          if (abs(Mdiag(it)) < 1.0e-24_wp) Mdiag(it) = 1.0e-24_wp
          Mdiag(it) = 1.0_wp / Mdiag(it)
       end do
 
       ! initial guess zero
       vrhs = 0.0_wp
 
       ! r = b - A x
       call gemv(amat, vrhs, Ap, alpha=1.0_wp, beta=0.0_wp, trans='n')
       r = xvec - Ap
 
       ! z = M^{-1} r
       z = r * Mdiag
       p = z
       bnorm = sqrt(sum(xvec*xvec))
       if (bnorm == 0.0_wp) bnorm = 1.0_wp
       rnorm = sqrt(sum(r*r))
       if (rnorm / bnorm <= tol) then
          if (present(info)) info = 0_ik
          return
       end if
 
       real(wp) :: rz_old, rz_new
       rz_old = sum(r*z)
       local_info = -1_ik
 
       do it = 1, maxit
          call gemv(amat, p, Ap, alpha=1.0_wp, beta=0.0_wp, trans='n')
          denom = sum(p * Ap)
          if (abs(denom) < 1.0e-30_wp) then
             local_info = -2_ik
             exit
          end if
          alpha = rz_old / denom
 
          vrhs = vrhs + alpha * p
          r = r - alpha * Ap
 
          rnorm = sqrt(sum(r*r))
          if (rnorm / bnorm <= tol) then
             local_info = 0_ik
             exit
          end if
 
          z = r * Mdiag
          rz_new = sum(r * z)
          beta = rz_new / rz_old
          p = z + beta * p
          rz_old = rz_new
       end do
 
       if (present(info)) info = local_info
    end subroutine solve_cg
 
    ! Factory: create solver (use_cg optional). If not present, env var MCHARGE_SOLVER="CG" selects cg.
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