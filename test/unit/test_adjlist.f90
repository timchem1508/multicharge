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

module test_adjlist
   use iso_fortran_env, only: output_unit
   use mctc_env, only: wp
   use mctc_env_testing, only: new_unittest, unittest_type, error_type, test_failed
   use mctc_cutoff, only: get_lattice_points
   use mctc_io_structure, only: structure_type, new
   use mctc_ncoord, only: adjacency_list, new_adjacency_list
   use mstore, only: get_structure
   use multicharge_blas, only: gemv
   use multicharge_wignerseitz, only: new_wignerseitz_cell
   use multicharge_model_type, only: mchrg_model_type
   use multicharge_model_eeqbc, only: eeqbc_model
   use multicharge_param, only: new_eeq2019_model, new_eeqbc2025_model
   use multicharge_model_cache, only: mchrg_cache
   use multicharge_charge, only: get_charges, get_eeq_charges, get_eeqbc_charges
   use multicharge_solver_type, only: mchrg_solver_type, mchrg_solver_input
   use multicharge_solver_direct, only : direct_solver, new_direct_solver, direct_input
   use multicharge_solver_cg, only : cg_solver, new_cg_solver, cg_input
   use multicharge_adjlist, only: symv_sparse
   implicit none
   private

   public :: collect_adjlist

   real(wp), parameter :: thr = 100 * epsilon(1.0_wp)
   real(wp), parameter :: thr2 = sqrt(epsilon(1.0_wp))
   real(wp), parameter :: thr3 = 100 * thr2
   real(wp), parameter :: cutoff = 29.0_wp

contains

!> Collect all exported unit tests
   subroutine collect_adjlist(testsuite)

      !> Collection of tests
      type(unittest_type), allocatable, intent(out) :: testsuite(:)

      testsuite = [ &
      & new_unittest("eeqbc-charges-mb01", test_eeqbc_q_mb01), &
      & new_unittest("eeqbc-charges-mb02", test_eeqbc_q_mb02), &
      & new_unittest("eeqbc-charges-actinides", test_eeqbc_q_actinides), &
      & new_unittest("eeqbc-energy-mb03", test_eeqbc_e_mb03), &
      & new_unittest("eeqbc-energy-mb04", test_eeqbc_e_mb04), &
      & new_unittest("eeqbc-gradient-mb05", test_eeqbc_g_mb05), &
      & new_unittest("eeqbc-gradient-mb06", test_eeqbc_g_mb06), &
      & new_unittest("eeqbc-energy-co2", test_eeqbc_e_co2) &
      !& new_unittest("eeqbc-gradient-co2", test_eeqbc_g_co2) &
      & ]

   end subroutine collect_adjlist

   subroutine solver_maker(solver, input, error)
      !> Solver type
      class(mchrg_solver_type), intent(out), allocatable :: solver
      !> Solver input
      class(mchrg_solver_input), intent(in) :: input
      !> Error handling
      type(error_type), allocatable, intent(out) :: error

      select type (input)
       type is (cg_input)
         block
            class(cg_solver), allocatable :: tmp
            allocate(tmp)
            call new_cg_solver(tmp, input)
            call move_alloc(tmp, solver)
         end block
       type is (direct_input)
         block
            class(direct_solver), allocatable :: tmp
            allocate(tmp)
            call new_direct_solver(tmp, input)
            call move_alloc(tmp, solver)
         end block
       class default
         allocate(error)
         return
      end select

   end subroutine solver_maker

!------------------------------------------------------------------------
! General helper routines – now they accept a pre‑built adjacency list.
!------------------------------------------------------------------------
   subroutine gen_test_molecular(error, mol, model, qref, eref)

      !> Error handling
      type(error_type), allocatable, intent(out) :: error

      !> Molecular structure data
      type(structure_type), intent(in) :: mol

      !> Electronegativity equilibration model
      class(mchrg_model_type), intent(in) :: model

      !> Reference charges
      real(wp), intent(in), optional :: qref(:)

      !> Reference energies
      real(wp), intent(in), optional :: eref(:)

      type(mchrg_cache), allocatable :: cache

      ! Solver variables
      class(mchrg_solver_type), allocatable :: solver
      class(mchrg_solver_input), allocatable :: solver_input
      real(wp) :: tol = 1.0e-15_wp
      integer :: maxiter = 1000
      integer :: verbosity = 0

      type(adjacency_list), allocatable :: list

      real(wp) :: trans(3, 1) = 0.0_wp
      real(wp), allocatable :: energy(:)
      real(wp), allocatable :: qvec(:)

      allocate(cg_input :: solver_input)
      select type (solver_input)
       type is (cg_input)
         solver_input%cgtol = tol
         solver_input%cgmiter = maxiter
         solver_input%verbosity = verbosity
         solver_input%use_nlist = .true.
      end select
      call solver_maker(solver, solver_input, error)

      allocate(cache)

      if (present(eref)) then
         allocate (energy(mol%nat))
         energy(:) = 0.0_wp
      end if
      if (present(qref)) then
         allocate (qvec(mol%nat))
      end if

      ! Build adjacency list
      allocate(list)
      call new_adjacency_list(list, mol, cutoff, .false.)

      call model%update(mol, cache, trans, grad=.false., list=list)
      call model%solve(mol, solver, cache, error, energy=energy, qvec=qvec, list=list, unit=output_unit)
      if (allocated(error)) return

      if (present(qref)) then
         if (any(abs(qvec - qref) > thr)) then
            call test_failed(error, "Partial charges do not match")
            print'(a)', "Charges:"
            print'(3es21.14)', qvec
            print'(a)', "diff:"
            print'(3es21.14)', qvec - qref
         end if
      end if
      if (allocated(error)) return

      if (present(eref)) then
         if (any(abs(energy - eref) > thr)) then
            call test_failed(error, "Energies do not match")
            print'(a)', "Energy:"
            print'(3es21.14)', energy
            print'(a)', "diff:"
            print'(3es21.14)', energy - eref
         end if
      end if

   end subroutine gen_test_molecular

   subroutine gen_test_periodic(error, mol, model)

      !> Error handling
      type(error_type), allocatable, intent(out) :: error

      !> Molecular structure data
      type(structure_type), intent(in) :: mol

      !> Electronegativity equilibration model
      class(mchrg_model_type), intent(in) :: model

      !> Reference charges
      real(wp), allocatable :: qref(:)

      !> Reference energies
      real(wp), allocatable :: eref(:)

      type(mchrg_cache), allocatable :: cache

      type(adjacency_list), allocatable :: list

      ! Solver variables
      class(mchrg_solver_type), allocatable :: solver
      class(mchrg_solver_input), allocatable :: solver_input
      real(wp) :: tol = 1.0e-15_wp
      integer :: maxiter = 1000
      integer :: verbosity = 0

      integer :: ndim
      real(wp), allocatable :: trans(:, :)
      real(wp), allocatable :: energy(:)
      real(wp), allocatable :: qvec(:)
      real(wp), allocatable :: amat_dir(:,:)
      real(wp), allocatable :: amat_list(:,:)

      integer :: iat, jat, kat

      allocate(cg_input :: solver_input)
      select type (solver_input)
       type is (cg_input)
         solver_input%cgtol = tol
         solver_input%cgmiter = maxiter
         solver_input%verbosity = verbosity
         ndim = mol%nat
         solver_input%use_nlist = .true.
      end select
      call solver_maker(solver, solver_input, error)

      allocate(cache)

      call get_lattice_points(mol%periodic, mol%lattice, cutoff, trans)

      allocate(eref(mol%nat))
      eref(:) = 0.0_wp
      allocate(energy(mol%nat))
      energy(:) = 0.0_wp
      allocate (qref(mol%nat))
      allocate (qvec(mol%nat))

      allocate(amat_dir(mol%nat, mol%nat), amat_list(mol%nat, mol%nat))

      call model%update(mol, cache, trans, grad=.false.)
      !call model%get_capacitance_matrix(mol, mol%nat, cache)
      !call model%get_xvec( mol, mol%nat, cache)
      !call model%get_coulomb_matrix(mol, mol%nat, cache)
      !amat_dir(:,:) = cache%amat
      !write(*,*) "DIRECT AMAT"
      !print'(12es21.14)', amat_dir
      !write(*,*) "DIRECT XVEC"
      !write(*,'(es21.14)') cache%xvec
      !write(*,'(50("-"))')
      call model%solve(mol, solver, cache, error, energy=eref, qvec=qref, unit=output_unit)


      if (allocated(error)) return

      ! Build adjacency list
      deallocate(cache)
      allocate(cache)
      allocate(list)
      call new_adjacency_list(list, mol, cutoff , .false.)
      call model%update(mol, cache, trans, grad=.false., list=list)
      call model%solve(mol, solver, cache, error, energy=energy, qvec=qvec, list=list, unit=output_unit)

      if (any(abs(qvec - qref) > thr)) then
         call test_failed(error, "Partial charges do not match")
         print'(a)', "Charges:"
         print'(3es21.14)', qvec
         print'(a)', "diff:"
         print'(3es21.14)', qvec - qref
      end if

      if (any(abs(energy - eref) > thr)) then
         call test_failed(error, "Energies do not match")
         print'(a)', "Energy:"
         print'(3es21.14)', energy
         print'(a)', "diff:"
         print'(3es21.14)', energy - eref
      end if
      if (allocated(error)) return


   end subroutine gen_test_periodic

   subroutine test_numgrad(error, mol, model)

      !> Error handling
      type(error_type), allocatable, intent(out) :: error

      !> Molecular structure data
      type(structure_type), intent(inout) :: mol

      !> Electronegativity equilibration model
      class(mchrg_model_type), intent(in) :: model

      type(mchrg_cache), allocatable :: cache1, cache2

      ! Solver variables
      class(mchrg_solver_type), allocatable :: solver
      class(mchrg_solver_input), allocatable :: solver_input
      real(wp) :: tol = 1.0e-15_wp
      integer :: maxiter = 1000
      integer :: verbosity = 0

      type(adjacency_list), allocatable :: list

      integer :: iat, jat, kat, ic, ndim
      real(wp), parameter :: trans(3, 1) = 0.0_wp   ! dummy for non‑periodic systems
      real(wp), allocatable :: energy(:), gradient(:, :), sigma(:, :)
      real(wp), allocatable :: numgrad(:, :), numsigma(:, :)
      real(wp) :: er, el
      real(wp), allocatable :: dcmat_list(:, :, :), damat_list(:, :, :)
      real(wp), allocatable :: dxvec_list(:, :, :), dcn_list(:, :, :), dqloc_list(:, :, :)
      logical :: grad = .true.

      allocate(cg_input :: solver_input)
      select type (solver_input)
       type is (cg_input)
         solver_input%cgtol = tol
         solver_input%cgmiter = maxiter
         solver_input%verbosity = verbosity
         solver_input%use_nlist = .true.
         ndim = mol%nat
      end select
      call solver_maker(solver, solver_input, error)

      allocate(cache1)

      allocate (energy(mol%nat), gradient(3, mol%nat), sigma(3, 3))
      energy(:) = 0.0_wp
      gradient(:, :) = 0.0_wp
      sigma(:, :) = 0.0_wp

      allocate(damat_list(3, mol%nat, mol%nat))
      allocate(dcmat_list(3, mol%nat, mol%nat))
      allocate(dxvec_list(3, mol%nat, mol%nat))
      allocate(dcn_list(3, mol%nat, mol%nat), dqloc_list(3, mol%nat, mol%nat))
      damat_list(:,:,:) = 0.0_wp
      dcmat_list(:,:,:) = 0.0_wp
      dxvec_list(:,:,:) = 0.0_wp
      dcn_list(:,:,:) = 0.0_wp
      dqloc_list(:,:,:) = 0.0_wp
      ! Build adjacency list
      allocate(numsigma(3, 3), source=0.0_wp)
      allocate(numgrad(3, mol%nat), source=0.0_wp)

      call model%update(mol, cache1, trans, grad=.true.)
      call model%solve(mol, solver, cache1, error, &
      & gradient=numgrad, sigma=numsigma, unit=output_unit)
      print'(16es21.14)', cache1%dcndr(3, :, :)
      write(*, *)
      if (allocated(error)) return

      allocate(cache2)
      allocate(list)
      gradient = 0.0_wp
      sigma(:, :) = 0.0_wp
      call new_adjacency_list(list, mol, cutoff, .false.)
      call model%update(mol, cache2, trans, grad=.true., list=list)

      !call model%get_capacitance_matrix(mol, mol%nat, cache2, list=list)
      !call model%get_xvec( mol, mol%nat, cache2, list=list)
      !call model%get_coulomb_matrix(mol, mol%nat, cache2, list=list)

      call model%solve(mol, solver, cache2, error, &
      & gradient=gradient, sigma=sigma, list=list, unit=output_unit)

      !do iat = 1, mol%nat
!
      !   do kat = list%inl(iat) + 1, list%inl(iat) + list%nnl(iat)
      !      jat = list%nlat(kat)
      !      dcmat_list(:, iat, jat) = cache2%dcdrij(:, kat)
      !      dcmat_list(:, jat, iat) = cache2%dcdrji(:, kat)
      !      dqloc_list(:, iat, jat) = cache2%dqlocdrij(:, kat)
      !      dqloc_list(:, jat, iat) = cache2%dqlocdrji(:, kat)
      !      dcn_list(:, iat, jat) = cache2%dcndrij(:, kat)
      !      dcn_list(:, jat, iat) = cache2%dcndrji(:, kat)
      !      dxvec_list(:, iat, jat) = cache2%dxdrij(:, kat)
      !      dxvec_list(:, jat, iat) = cache2%dxdrji(:, kat)
      !      damat_list(:, iat, jat) = cache2%dadrij(:, kat)
      !      damat_list(:, jat, iat) = cache2%dadrji(:, kat)
      !   end do
!
      !   dcmat_list(:, iat, iat) = cache2%dcdrdiag(:, iat)
      !   damat_list(:, iat, iat) = cache2%dadrdiag(:, iat)
      !   dqloc_list(:, iat, iat) = cache2%dqlocdrdiag(:, iat)
      !   dxvec_list(:, iat, iat) = cache2%dxdrdiag(:, iat)
      !   dcn_list(:, iat, iat) = cache2%dcndrdiag(:, iat)
      !end do

      if (allocated(error)) return

      !if (any(abs(cache1%dqlocdr(:, :, :) - dqloc_list(:, :, :)) > thr3)) then
      !   call test_failed(error, "Derivative of local charge does not match")
      !   print'(a)', "Local charge derivative:"
      !   print'(3es21.14)', cache1%dqlocdr
      !   print'(a)', "Nlist local charge derivative:"
      !   print'(3es21.14)', dqloc_list
      !   print'(a)', "diff:"
      !   print'(3es21.14)', dqloc_list - cache1%dqlocdr
      !end if
!
      !if (any(abs(cache1%dcndr(:, :, :) - dcn_list(:, :, :)) > thr3)) then
      !   call test_failed(error, "Derivative of CN does not match")
      !   print'(a)', "CN derivative:"
      !   print'(3es21.14)', cache1%dcndr
      !   print'(a)', "Nlist CN derivative:"
      !   print'(3es21.14)', dcn_list
      !   print'(a)', "diff:"
      !   print'(3es21.14)', dcn_list - cache1%dcndr
      !end if
!
      !if (any(abs(cache1%dcdr(:, :, :) - dcmat_list(:, :, :)) > thr3)) then
      !   call test_failed(error, "Derivative of Capacitance matrix does not match")
      !   print'(a)', "Capacitance derivative:"
      !   print'(16es21.14)', cache1%dcdr
      !   print'(a)', "Nlist Capacitance derivative:"
      !   print'(16es21.14)', dcmat_list
      !   print'(a)', "diff:"
      !   print'(16es21.14)', dcmat_list - cache1%dcdr
      !end if
!
      !if (any(abs(cache1%dxdr(:, :, :) - dxvec_list(:, :, :)) > thr3)) then
      !   call test_failed(error, "Derivative of electronegativity does not match")
      !   print'(a)', "Electronegativity derivative:"
      !   print'(16es21.14)', cache1%dxdr(3, :, :)
      !   print'(a)', "Nlist electronegativity derivative:"
      !   print'(16es21.14)', dxvec_list(3, :, :)
      !   print'(a)', "diff:"
      !   print'(16es21.14)', dxvec_list(3, :, :) - cache1%dxdr(3, :, :)
      !end if
!
      !if (any(abs(cache1%dadr(:, :, :) - damat_list(:, :, :)) > thr3)) then
      !   call test_failed(error, "Derivative of Coulomb matrix does not match")
      !   print'(a)', "Coulomb matrix derivative:"
      !   print'(16es21.14)', cache1%dadr(3, :, :)
      !   print'(a)', "Nlist Coulomb matrix derivative:"
      !   print'(16es21.14)', damat_list( 3, :, :)
      !   print'(a)', "diff:"
      !   print'(16es21.14)', damat_list( 3, :, :) - cache1%dadr(3, :, :)
      !end if

      if (any(abs(gradient(:, :) - numgrad(:, :)) > thr3)) then
         call test_failed(error, "Derivative of energy does not match")
         print'(a)', "Energy gradient:"
         print'(3es21.14)', gradient
         print'(a)', "numgrad:"
         print'(3es21.14)', numgrad
         print'(a)', "diff:"
         print'(3es21.14)', gradient - numgrad
      end if

      if (any(abs(sigma(:, :) - numsigma(:, :)) > thr3)) then
         call test_failed(error, "Derivative of energy does not match")
         print'(a)', "Energy sigma:"
         print'(3es21.14)', sigma
         print'(a)', "numsigma:"
         print'(3es21.14)', numsigma
         print'(a)', "diff:"
         print'(3es21.14)', sigma - numsigma
      end if

   end subroutine test_numgrad

   subroutine test_numgrad_periodic(error, mol, model)

      !> Error handling
      type(error_type), allocatable, intent(out) :: error

      !> Molecular structure data
      type(structure_type), intent(inout) :: mol

      !> Electronegativity equilibration model
      class(mchrg_model_type), intent(in) :: model

      type(mchrg_cache), allocatable :: cache

      ! Solver variables
      class(mchrg_solver_type), allocatable :: solver
      class(mchrg_solver_input), allocatable :: solver_input
      real(wp) :: tol = 1.0e-15_wp
      integer :: maxiter = 1000
      integer :: verbosity = 0

      type(adjacency_list), allocatable :: list

      integer :: iat, ic, ndim
      real(wp), allocatable :: trans(:, :)
      real(wp), parameter :: step = 1.0e-6_wp
      real(wp), allocatable :: energy(:), gradient(:, :), sigma(:, :)
      real(wp), allocatable :: numgrad(:, :)
      real(wp) :: er, el
      logical :: grad = .true.

      allocate(cg_input :: solver_input)
      select type (solver_input)
       type is (cg_input)
         solver_input%cgtol = tol
         solver_input%cgmiter = maxiter
         solver_input%verbosity = verbosity
         solver_input%use_nlist = .true.
         ndim = mol%nat
      end select
      call solver_maker(solver, solver_input, error)

      allocate(cache)

      allocate (energy(mol%nat), gradient(3, mol%nat), sigma(3, 3), numgrad(3, mol%nat))
      energy(:) = 0.0_wp
      gradient(:, :) = 0.0_wp
      sigma(:, :) = 0.0_wp

      call get_lattice_points(mol%periodic, mol%lattice, 25.0_wp, trans)

      lp: do iat = 1, mol%nat
         do ic = 1, 3
            energy(:) = 0.0_wp
            mol%xyz(ic, iat) = mol%xyz(ic, iat) + step
            call model%update(mol, cache, trans, grad=.false.)
            call model%solve(mol, solver, cache, error, energy=energy, unit=output_unit)
            if (allocated(error)) exit lp
            er = sum(energy)

            energy(:) = 0.0_wp
            mol%xyz(ic, iat) = mol%xyz(ic, iat) - 2*step
            call model%update(mol, cache, trans, grad=.false.)
            call model%solve(mol, solver, cache, error, energy=energy, unit=output_unit)
            if (allocated(error)) exit lp
            el = sum(energy)

            mol%xyz(ic, iat) = mol%xyz(ic, iat) + step
            numgrad(ic, iat) = 0.5_wp*(er - el)/step
         end do
      end do lp
      if (allocated(error)) return

      ! Build adjacency list
      allocate(list)
      call new_adjacency_list(list, mol, cutoff, .false.)

      call model%update(mol, cache, trans, grad, list=list)
      call model%solve(mol, solver, cache, error, &
      & gradient=gradient, sigma=sigma, list=list, unit=output_unit)
      if (allocated(error)) return

      if (any(abs(gradient(:, :) - numgrad(:, :)) > thr3)) then
         call test_failed(error, "Derivative of energy does not match")
         print'(a)', "Energy gradient:"
         print'(3es21.14)', gradient
         print'(a)', "numgrad:"
         print'(3es21.14)', numgrad
         print'(a)', "diff:"
         print'(3es21.14)', gradient - numgrad
      end if

   end subroutine test_numgrad_periodic

   subroutine test_dadr(error, mol, model)


      !> Molecular structure data
      type(structure_type), intent(inout) :: mol

      !> Electronegativity equilibration model
      class(mchrg_model_type), intent(in) :: model

      !> Error handling
      type(error_type), allocatable, intent(out) :: error

      ! Solver variables
      class(mchrg_solver_type), allocatable :: solver
      class(mchrg_solver_input), allocatable :: solver_input
      real(wp) :: tol = 1.0e-15_wp
      integer :: maxiter = 1000
      integer :: verbosity = 0

      integer :: iat, ic, jat, kat, ndim
      real(wp) :: thr2_local
      real(wp), parameter :: step = 1.0e-6_wp
      real(wp), allocatable :: trans(:, :)
      real(wp), allocatable :: qvec(:), numgrad(:, :, :),  numtrace(:, :)
      real(wp), allocatable :: amatr1(:, :), amatr2(:, :), amatl1(:, :), amatl2(:, :)
      type(mchrg_cache), allocatable :: cache
      logical :: grad = .true.

      allocate(cg_input :: solver_input)
      select type (solver_input)
       type is (cg_input)
         solver_input%cgtol = tol
         solver_input%cgmiter = maxiter
         solver_input%verbosity = verbosity
         ndim = mol%nat
      end select
      call solver_maker(solver, solver_input, error)

      allocate (cache)

      allocate (amatr1(ndim, ndim), amatl1(ndim, ndim), amatr2(ndim, ndim), amatl2(ndim, ndim), &
      & numtrace(3, mol%nat), numgrad(3, mol%nat, ndim), qvec(mol%nat))

      ! Set tolerance higher if testing eeqbc model
      select type (model)
       type is (eeqbc_model)
         thr2_local = 3.0_wp*thr2
       class default
         thr2_local = thr2
      end select

      call get_lattice_points(mol%periodic, mol%lattice, cutoff, trans)

      ! Obtain the vector of charges
      call model%update(mol, cache, trans, grad=.false.)
      call model%solve(mol, solver, cache, error, qvec=qvec, unit=output_unit)
      if (allocated(error)) return

      ! Analytical gradient
      call model%update(mol, cache, trans, grad)
      call model%get_capacitance_matrix(mol, ndim, cache)
      call model%get_coulomb_derivs(mol, ndim, cache)

      if (any(abs(cache%dadr(:, :, :) - numgrad(:, :, :)) > thr2_local)) then
         call test_failed(error, "Derivative of the A matrix does not match")
         print'(a)', "dadr:"
         print'(3es21.12)', cache%dadr
         print'(a)', "numgrad:"
         print'(3es21.12)', numgrad
         print'(a)', "diff:"
         print'(3es21.12)', cache%dadr - numgrad
      end if

   end subroutine test_dadr

!------------------------------------------------------------------------
! Test routines – now each builds its own adjacency list.
!------------------------------------------------------------------------
   subroutine test_eeqbc_q_mb01(error)

      !> Error handling
      type(error_type), allocatable, intent(out) :: error

      !> Molecular structure data
      type(structure_type) :: mol
      class(mchrg_model_type), allocatable :: model

      real(wp), parameter :: ref(16) = [&
      & 4.75783090912440E-1_wp, -4.26540500638442E-2_wp, -3.77871226005535E-1_wp, &
      &-9.67376090029522E-2_wp, -1.73364116997142E-1_wp, 1.08660101025683E-1_wp, &
      &-1.13628448410420E-1_wp, -3.17939699645693E-1_wp, -2.45655524697400E-1_wp, &
      & 1.76106572419156E-1_wp, 1.14510850652006E-1_wp, -1.22241025474265E-1_wp, &
      &-1.44595425453640E-2_wp, 2.57782082780412E-1_wp, -1.11777579535162E-1_wp, &
      & 4.83486124588080E-1_wp]

      real(wp), allocatable :: qvec(:)
      real(wp), parameter :: trans(3, 1) = 0.0_wp

      call get_structure(mol, "MB16-43", "01")

      call new_eeqbc2025_model(mol, model, error)
      if (allocated(error)) return
      call gen_test_molecular(error, mol, model, qref=ref)

   end subroutine test_eeqbc_q_mb01

   subroutine test_eeqbc_q_mb02(error)

      !> Error handling
      type(error_type), allocatable, intent(out) :: error

      !> Molecular structure data
      type(structure_type) :: mol
      class(mchrg_model_type), allocatable :: model

      real(wp), parameter :: ref(16) = [&
      &-7.89571755894845E-2_wp, -1.84724587297173E-1_wp, -1.63060175795952E-2_wp, &
      &-2.36115890461711E-1_wp, 5.05729582512203E-1_wp, 1.37556939519704E-1_wp, &
      &-2.29048340967271E-2_wp, -4.31722346626804E-2_wp, 2.26466952977883E-1_wp, &
      & 1.25047857913714E-1_wp, 6.72899182661252E-3_wp, 3.08986208662492E-1_wp, &
      &-3.34344661086462E-1_wp, -3.16758668376149E-2_wp, -5.24170403450005E-2_wp, &
      &-3.09898225456160E-1_wp]

      real(wp), parameter :: trans(3, 1) = 0.0_wp

      call get_structure(mol, "MB16-43", "02")

      call new_eeqbc2025_model(mol, model, error)
      if (allocated(error)) return
      call gen_test_molecular(error, mol, model, qref=ref)

   end subroutine test_eeqbc_q_mb02

   subroutine test_eeqbc_q_actinides(error)

      !> Error handling
      type(error_type), allocatable, intent(out) :: error

      type(structure_type) :: mol
      class(mchrg_model_type), allocatable :: model

      real(wp), parameter :: ref(17) = [&
      & 9.27195802124755E-2_wp, -2.78358027117801E-1_wp, 1.71815557281178E-1_wp, &
      & 7.85579953672371E-2_wp, -1.08186262417305E-2_wp, -4.81860290986309E-2_wp, &
      & 1.57794666483371E-1_wp, -1.61830258916072E-1_wp, -2.76569765724910E-1_wp, &
      & 2.99654899926371E-1_wp, -5.24433579322476E-1_wp, -1.99523360511699E-1_wp, &
      &-3.42285450387671E-2_wp, -3.15076271542101E-2_wp, 1.49700940990172E-1_wp, &
      & 1.45447393911445E-1_wp, 4.69764784954047E-1_wp]

      real(wp), parameter :: trans(3, 1) = 0.0_wp

      ! Molecular structure data
      mol%nat = 17
      mol%nid = 17
      mol%id = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, &
      & 12, 13, 14, 15, 16, 17]
      mol%num = [87, 88, 89, 90, 91, 92, 93, 94, 95, 96, 97, &
      & 98, 99, 100, 101, 102, 103]
      mol%xyz = reshape([ &
      & 0.98692316414074_wp, 6.12727238368797_wp, -6.67861597188102_wp, &
      & 3.63898862390869_wp, 5.12109301182962_wp, 3.01908613326278_wp, &
      & 5.14503571563551_wp, -3.97172984617710_wp, 3.82011791828867_wp, &
      & 6.71986847575494_wp, 1.71382138402812_wp, 3.92749159076307_wp, &
      & 4.13783589704826_wp, -2.10695793491818_wp, 0.19753203068899_wp, &
      & 8.97685097698326_wp, -3.08813636191844_wp, -4.45568615593938_wp, &
      & 12.5486412940776_wp, -1.77128765259458_wp, 0.59261498922861_wp, &
      & 7.82051475868325_wp, -3.97159756604558_wp, -0.53637703616916_wp, &
      &-0.43444574624893_wp, -1.69696511583960_wp, -1.65898182093050_wp, &
      &-4.71270645149099_wp, -0.11534827468942_wp, 2.84863373521297_wp, &
      &-2.52061680335614_wp, 1.82937752749537_wp, -2.10366982879172_wp, &
      & 0.13551154616576_wp, 7.99805359235043_wp, -1.55508522619903_wp, &
      & 3.91594542499717_wp, -1.72975169129597_wp, -5.07944366756113_wp, &
      &-1.03393930231679_wp, 4.69307230054046_wp, 0.02656940927472_wp, &
      & 6.20675384557240_wp, 4.24490721493632_wp, -0.71004195169885_wp, &
      & 7.04586341131562_wp, 5.20053667939076_wp, -7.51972863675876_wp, &
      & 2.01082807362334_wp, 1.34838807211157_wp, -4.70482633508447_wp],&
      & [3, 17])
      mol%periodic = [.false.]

      call new_eeqbc2025_model(mol, model, error)
      if (allocated(error)) return
      call gen_test_molecular(error, mol, model, qref=ref)

   end subroutine test_eeqbc_q_actinides

   subroutine test_eeqbc_e_mb03(error)

      !> Error handling
      type(error_type), allocatable, intent(out) :: error

      type(structure_type) :: mol
      class(mchrg_model_type), allocatable :: model

      real(wp), parameter :: ref(16) = [&
      &-6.96992195046228E-2_wp, -1.62155815983893E+0_wp, -1.38060751929644E-3_wp, &
      &-9.06342279911342E-1_wp, -1.83281566961757E+0_wp, -1.20333262207652E-1_wp, &
      &-6.51187555181622E-1_wp, -3.27410111288548E-3_wp, -8.00565881078213E-3_wp, &
      &-2.60385867643294E-2_wp, -9.33285940415006E-1_wp, -1.48859947660327E-1_wp, &
      &-7.19456827995756E-1_wp, -9.58311834831915E-2_wp, -1.54672086637309E+0_wp, &
      &-1.03483694342593E-5_wp]

      real(wp), parameter :: trans(3, 1) = 0.0_wp

      call get_structure(mol, "MB16-43", "03")

      call new_eeqbc2025_model(mol, model, error)
      if (allocated(error)) return
      call gen_test_molecular(error, mol, model, eref=ref)

   end subroutine test_eeqbc_e_mb03

   subroutine test_eeqbc_e_mb04(error)

      !> Error handling
      type(error_type), allocatable, intent(out) :: error

      type(structure_type) :: mol
      class(mchrg_model_type), allocatable :: model

      real(wp), parameter :: ref(16) = [&
      &-3.91054587109712E-2_wp, -8.21933095021462E-4_wp, -1.28550631772418E-2_wp, &
      &-8.95571658260288E-2_wp, -4.94655224590082E-1_wp, -3.34598696522549E-2_wp, &
      &-3.75768676247744E-2_wp, -1.36087478076862E-2_wp, -2.07985587717960E-3_wp, &
      &-1.17711924662077E-2_wp, -2.68707428024071E-1_wp, -1.00650791933494E+0_wp, &
      &-5.64487253409848E-2_wp, -4.89693252471477E-1_wp, -3.74734977139679E-2_wp, &
      &-9.22642011641358E-3_wp]

      real(wp), parameter :: trans(3, 1) = 0.0_wp

      call get_structure(mol, "MB16-43", "04")

      call new_eeqbc2025_model(mol, model, error)
      if (allocated(error)) return
      call gen_test_molecular(error, mol, model, eref=ref)

   end subroutine test_eeqbc_e_mb04

   subroutine test_eeqbc_g_mb05(error)

      !> Error handling
      type(error_type), allocatable, intent(out) :: error

      type(structure_type) :: mol
      class(mchrg_model_type), allocatable :: model

      real(wp), parameter :: trans(3, 1) = 0.0_wp

      call get_structure(mol, "MB16-43", "05")

      call new_eeqbc2025_model(mol, model, error)
      if (allocated(error)) return
      call test_numgrad(error, mol, model)

   end subroutine test_eeqbc_g_mb05

   subroutine test_eeqbc_g_mb06(error)

      !> Error handling
      type(error_type), allocatable, intent(out) :: error

      type(structure_type) :: mol
      class(mchrg_model_type), allocatable :: model

      real(wp), parameter :: trans(3, 1) = 0.0_wp

      call get_structure(mol, "MB16-43", "06")

      call new_eeqbc2025_model(mol, model, error)
      if (allocated(error)) return
      call test_numgrad(error, mol, model)

   end subroutine test_eeqbc_g_mb06

   subroutine test_eeqbc_g_co2(error)

      !> Error handling
      type(error_type), allocatable, intent(out) :: error

      type(structure_type) :: mol
      class(mchrg_model_type), allocatable :: model

      call get_structure(mol, "X23", "CO2")
      call new_eeqbc2025_model(mol, model, error)
      if (allocated(error)) return
      call test_numgrad_periodic(error, mol, model)

   end subroutine test_eeqbc_g_co2

   subroutine test_eeqbc_e_co2(error)

      !> Error handling
      type(error_type), allocatable, intent(out) :: error

      type(structure_type) :: mol
      class(mchrg_model_type), allocatable :: model

      call get_structure(mol, "X23", "CO2")
      call new_eeqbc2025_model(mol, model, error)
      if (allocated(error)) return
      call gen_test_periodic(error, mol, model)

   end subroutine test_eeqbc_e_co2

end module test_adjlist
