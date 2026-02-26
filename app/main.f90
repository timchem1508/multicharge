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

program main
   use, intrinsic :: iso_fortran_env, only: output_unit, error_unit, input_unit
   use mctc_env, only: error_type, fatal_error, get_argument, wp
   use mctc_env_timer, only : timer_type, format_time
   use mctc_io, only: structure_type, read_structure, filetype, get_filetype
   use mctc_cutoff, only: get_lattice_points
   use multicharge, only: mchrg_model_type, mchrg_model, new_eeq2019_model, &
      & new_eeqbc2025_model, get_multicharge_version, &
      & write_ascii_model, write_ascii_properties, write_ascii_results
   use multicharge_output, only: json_results, write_cg_solver, write_direct_solver
   use multicharge_solver, only: new_mchrg_solver, mchrg_solver_type, mchrg_solver_direct, &
      & mchrg_solver_cg, mchrg_solver_input, cg_input, direct_input


   implicit none
   character(len=*), parameter :: prog_name = "multicharge"
   character(len=*), parameter :: json_output = "multicharge.json"

   character(len=:), allocatable :: input, chargeinput
   integer, allocatable :: input_format
   integer :: stat, unit, model_id
   type(error_type), allocatable :: error
   type(structure_type) :: mol
   class(mchrg_model_type), allocatable :: model
   class(mchrg_solver_type), allocatable :: solver
   class(mchrg_solver_input), allocatable :: solver_input
   logical :: grad, json, exist
   real(wp), parameter :: cn_max = 8.0_wp, cutoff = 25.0_wp
   real(wp), allocatable :: cn(:), rcov(:), trans(:, :)
   real(wp), allocatable :: qloc(:)
   real(wp), allocatable :: dcndr(:, :, :), dcndL(:, :, :), dqlocdr(:, :, :), dqlocdL(:, :, :)
   real(wp), allocatable :: energy(:), gradient(:, :), sigma(:, :)
   real(wp), allocatable :: qvec(:)
   real(wp), allocatable :: dqdr(:, :, :), dqdL(:, :, :)
   real(wp), allocatable :: charge

   ! Solver configuration
   integer, allocatable :: maxiter
   real(wp), allocatable :: tol

   ! Verbose level
   integer, allocatable :: verbose

   ! Timer 
   type(timer_type) :: timer

   call timer%push("total")

   ! 1. Parse Arguments
   call get_arguments(input, model_id, input_format, grad, charge, json, &
                      solver_input, verbose, error)
   if (allocated(error)) then
      write(error_unit, '(a)') error%message
      error stop
   end if

   ! 2. Initialize Solver using the solver_input

   call new_mchrg_solver(solver, solver_input,  error)

   ! 3. Load Structure
   if (input == "-") then
      if (.not. allocated(input_format)) input_format = filetype%xyz
      call read_structure(mol, input_unit, input_format, error)
   else
      call read_structure(mol, input, error, input_format)
   end if
   if (allocated(error)) then
      write(error_unit, '(a)') error%message
      error stop
   end if

   if (allocated(charge)) then
      mol%charge = charge
   else
      chargeinput = ".CHRG"
      inquire(file=chargeinput, exist=exist)
      if (exist) then
         open(file=chargeinput, newunit=unit)
         allocate(charge)
         read(unit, *, iostat=stat) charge
         if (stat == 0) then
            mol%charge = charge
            write(output_unit, '(a,/)') &
               "[Info] Molecular charge read from '"//chargeinput//"'"
            else
               write(output_unit, '(a,/)') &
                  "[Warn] Could not read molecular charge read from '"//chargeinput//"'"
         end if
         close(unit)
      end if
   end if

   ! 4. Initialize Model
   if (model_id == mchrg_model%eeq2019) then
      call new_eeq2019_model(mol, model, error)
   else if (model_id == mchrg_model%eeqbc2025) then
      call new_eeqbc2025_model(mol, model, error)
   else
      call fatal_error(error, "Invalid model was choosen.")
   end if
   if (allocated(error)) then
      write(error_unit, '(a)') error%message
      error stop
   end if

   call write_ascii_model(output_unit, mol, model)

   allocate(energy(mol%nat), qvec(mol%nat))
   energy(:) = 0.0_wp

   allocate(cn(mol%nat), qloc(mol%nat))
   if (grad) then
      allocate(gradient(3, mol%nat), sigma(3, 3))
      gradient(:, :) = 0.0_wp
      sigma(:, :) = 0.0_wp
      allocate(dqdr(3, mol%nat, mol%nat), dqdL(3, 3, mol%nat))
      dqdr = 0.0_wp
      dqdL = 0.0_wp
      allocate(dcndr(3, mol%nat, mol%nat), dcndL(3, 3, mol%nat))
      allocate(dqlocdr(3, mol%nat, mol%nat), dqlocdL(3, 3, mol%nat))
   end if

   call get_lattice_points(mol%periodic, mol%lattice, model%ncoord%cutoff, trans)
   call model%ncoord%get_coordination_number(mol, trans, cn, dcndr, dcndL)
   call model%local_charge(mol, trans, qloc, dqlocdr, dqlocdL)
   
   ! 5. Run Solve (Solver instance passed implicitly via argument or model)

   if (verbose > 0) then
      select type (solver)
      type is (mchrg_solver_cg)
         call write_cg_solver(output_unit, solver)
      type is (mchrg_solver_direct)
         call write_direct_solver(output_unit, solver)
      end select
   end if 

   call model%solve(mol, solver, error, cn, qloc, dcndr, dcndL, dqlocdr, dqlocdL, &
      & energy, gradient, sigma, qvec, dqdr, dqdL, verbose=verbose)

   if (allocated(error)) then
      write(error_unit, '(a)') error%message
      error stop
   end if

   call write_ascii_properties(output_unit, mol, model, cn, qvec)
   call write_ascii_results(output_unit, mol, energy, gradient, sigma)

   if (json) then
      open(file=json_output, newunit=unit)
      call json_results(unit, "  ", energy=sum(energy), gradient=gradient, charges=qvec, cn=cn)
      close(unit)
      write(output_unit, '(a)') &
         "[Info] JSON dump of results written to '"//json_output//"'"
   end if

   call timer%pop
   if (verbose > 1) then
      write(*, '(a, 1x, a)') "Total execution time : ", format_time(timer%get("total"))
      write(*,*)
   end if

contains

subroutine help(unit)
   integer, intent(in) :: unit

   write(unit, '(a, *(1x, a))') &
      "Usage: "//prog_name//" [options] <input>"

   write(unit, '(a)') &
      "", &
      "Electronegativity equilibration model for atomic charges and", &
      "higher multipole moments", &
      ""

   write(unit, '(2x, a, t35, a)') &
      "-m, -model, --model <model>", "Choose the charge model", &
      "-i, -input, --input <format>", "Hint for the format of the input file", &
      "-c, -charge, --charge <value>", "Set the molecular charge", &
      "-g, -grad, --grad", "Evaluate molecular gradient and virial. Only for the direct solver.", &
      "-j, -json, --json", "Provide output in JSON format to the file 'multicharge.json'", &
      "-version, --version", "Print program version and exit", &
      "-v, -verbosity, --verbosity <int>", "Size of an output and exit", &
      "-h, -help, --help", "Show this help message", &
      "-s, -solver, --solver <type>", "Solver: 'CG' or 'DIRECT'", &
      "-it, -maxiter, --maxiter <int>", "Max iterations", &
      "-tol, -tolerance, --tolerance <real>", "Tolerance"
   write(unit, '(a)')

end subroutine help

subroutine version(unit)
   integer, intent(in) :: unit
   character(len=:), allocatable :: version_string

   call get_multicharge_version(string=version_string)
   write(unit, '(a, *(1x, a))') &
      & prog_name, "version", version_string

end subroutine version

subroutine get_arguments(input, model_id, input_format, grad, charge, &
   & json, solver_input, verbose, error)

   !> Input file name
   character(len=:), allocatable :: input
   !> ID of choosen model type
   integer, intent(out) :: model_id
   !> Input file format
   integer, allocatable, intent(out) :: input_format
   !> Evaluate gradient
   logical, intent(out) :: grad
   !> Provide JSON output
   logical, intent(out) :: json
   !> Charge
   real(wp), allocatable, intent(out) :: charge
   !> Solver args
   class(mchrg_solver_input), allocatable, intent(out) :: solver_input
   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   integer :: iarg, narg, iostat
   character(len=:), allocatable :: arg

   !> Solver type: CG or DIRECT
   character(len=:), allocatable :: solver_name
   !> Maximal number of the cg solver iterations
   integer, allocatable :: maxiter
   !> CG solver tolerance
   real(wp), allocatable :: tol

   ! Verbose number
   integer, allocatable :: verbose

   model_id = mchrg_model%eeq2019
   grad = .false.
   json = .false.
   iarg = 0
   verbose = 0
   narg = command_argument_count()

   do while(iarg < narg)
      iarg = iarg + 1
      call get_argument(iarg, arg)
      select case(arg)
      case("-h", "-help", "--help")
         call help(output_unit)
         stop
      case("-version", "--version")
         call version(output_unit)
         stop
      case("-v", "-verbosity", "--verbosity")
         iarg = iarg + 1
         call get_argument(iarg, arg)
         read(arg, *, iostat=iostat) verbose
         if (iostat >  1) then
            call fatal_error(error, "Invalid verbose level")
            exit
         end if
      case default
         if (.not. allocated(input)) then
            call move_alloc(arg, input)
            cycle
         end if
         call fatal_error(error, "Too many positional arguments present")
         exit
      case("-m", "-model", "--model")
         iarg = iarg + 1
         call get_argument(iarg, arg)
         if (.not. allocated(arg)) then
            call fatal_error(error, "Missing argument for model")
            exit
         end if
         if (arg == "eeq2019" .or. arg == "eeq") then
            model_id = mchrg_model%eeq2019
         else if (arg == "eeqbc2025" .or. arg == "eeqbc") then
            model_id = mchrg_model%eeqbc2025
         else
            call fatal_error(error, "Invalid model")
            exit
         end if
      case("-i", "-input", "--input")
         iarg = iarg + 1
         call get_argument(iarg, arg)
         if (.not. allocated(arg)) then
            call fatal_error(error, "Missing argument for input format")
            exit
         end if
         input_format = get_filetype("."//arg)
      case("-c", "-charge", "--charge")
         iarg = iarg + 1
         call get_argument(iarg, arg)
         if (.not. allocated(arg)) then
            call fatal_error(error, "Missing argument for charge")
            exit
         end if
         allocate(charge)
         read(arg, *, iostat=iostat) charge
         if (iostat /= 0) then
            call fatal_error(error, "Invalid charge value")
            exit
         end if
      case("-g", "-grad", "--grad")
         grad = .true.
      case("-j", "-json", "--json")
         json = .true.
      case("-s", "-solver", "--solver")
         if (allocated(solver_name)) then
            call fatal_error(error, "Cannot use multiple solvers")
            exit
         end if
         iarg = iarg + 1; call get_argument(iarg, solver_name)
         if (solver_name == "DIRECT" .or. solver_name == "direct" .or. solver_name == "LAPACK" ) then
            allocate(direct_input :: solver_input)
         end if
         if (solver_name == "CG" .or. solver_name == "cg" .or. solver_name == "iterative") then
              allocate(cg_input :: solver_input)
         end if
      case("-it", "-maxiter", "--maxiter")
         allocate(maxiter)
         iarg = iarg + 1; call get_argument(iarg, arg)
         read(arg, *, iostat=iostat) maxiter
         if (iostat /= 0) call fatal_error(error, "Invalid maximal number of iterations")
      case("-tol", "-tolerance", "--tolerance")
         allocate(tol)
         iarg = iarg + 1; call get_argument(iarg, arg)
         read(arg, *, iostat=iostat) tol
         if (iostat /= 0) call fatal_error(error, "Invalid tolerance")
      end select
   end do

   if (.not. allocated(solver_name)) then
      allocate(direct_input :: solver_input)
   end if

   if (grad .eqv. .true.) then
      if (allocated(solver_input)) deallocate(solver_input)
      allocate(direct_input :: solver_input)
   end if

   select type(solver_input)
   type is (cg_input)
   if (allocated(maxiter)) then
      solver_input%cgmiter = maxiter
   end if
   if (allocated(tol)) then
      solver_input%cgtol = tol
   end if
   if (allocated(verbose)) then
      solver_input%verbose = verbose
   end if
   end select

   
   if (.not. allocated(input)) then
      if (.not. allocated(error)) then
         call help(output_unit)
         error stop
      end if
   end if

end subroutine get_arguments

end program main