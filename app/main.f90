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
   use mctc_io, only: structure_type, read_structure, filetype, get_filetype
   use mctc_cutoff, only: get_lattice_points
   use multicharge, only: mchrg_model_type, mchrg_model, new_eeq2019_model, &
      & new_eeqbc2025_model, get_multicharge_version, &
      & write_ascii_model, write_ascii_properties, write_ascii_results
   use multicharge_output, only: json_results
   use solver, only: mchrg_solver_type
   use solver_factory, only: new_mchrg_solver

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
   logical :: grad, json, exist
   real(wp), parameter :: cn_max = 8.0_wp, cutoff = 25.0_wp
   real(wp), allocatable :: cn(:), trans(:, :)
   real(wp), allocatable :: qloc(:)
   real(wp), allocatable :: dcndr(:, :, :), dcndL(:, :, :), dqlocdr(:, :, :), dqlocdL(:, :, :)
   real(wp), allocatable :: energy(:), gradient(:, :), sigma(:, :)
   real(wp), allocatable :: qvec(:)
   real(wp), allocatable :: dqdr(:, :, :), dqdL(:, :, :)
   real(wp), allocatable :: charge

   ! Solver configuration
   character(len=:), allocatable :: solver_choice
   integer, allocatable :: cgmiter
   real(wp), allocatable :: cgtol
   character(len=32), allocatable :: cgmode

   ! 1. Parse Arguments
   call get_arguments(input, model_id, input_format, grad, charge, json, &
                      solver_choice, cgmiter, cgtol, cgmode, error)
   if (allocated(error)) then
      write(error_unit, '(a)') error%message
      error stop
   end if

   ! 2. Initialize Solver using factory with parsed arguments
   if (.not. allocated(solver_choice)) solver_choice = "CG" ! Default
   if (.not. allocated(cgmiter)) cgmiter = 1000 ! Default
   if (.not. allocated(cgtol)) cgtol = 1.0e-11_wp ! Default
   if (.not. allocated(cgmode)) cgmode = "default" ! Default
   solver = new_mchrg_solver(solver_choice, cgmiter, cgtol, cgmode)

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
            write(output_unit, '(a,/)') "[Info] Molecular charge read from '"//chargeinput//"'"
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
   if (allocated(error)) error stop

   call write_ascii_model(output_unit, mol, model)

   allocate(energy(mol%nat), qvec(mol%nat))
   energy(:) = 0.0_wp

   allocate(cn(mol%nat), qloc(mol%nat))
   if (grad) then
      allocate(gradient(3, mol%nat), sigma(3, 3))
      gradient(:, :) = 0.0_wp; sigma(:, :) = 0.0_wp
      allocate(dqdr(3, mol%nat, mol%nat), dqdL(3, 3, mol%nat))
      dqdr = 0.0_wp; dqdL = 0.0_wp
      allocate(dcndr(3, mol%nat, mol%nat), dcndL(3, 3, mol%nat))
      allocate(dqlocdr(3, mol%nat, mol%nat), dqlocdL(3, 3, mol%nat))
   end if

   call get_lattice_points(mol%periodic, mol%lattice, model%ncoord%cutoff, trans)
   call model%ncoord%get_coordination_number(mol, trans, cn, dcndr, dcndL)
   call model%local_charge(mol, trans, qloc, dqlocdr, dqlocdL)
   
   ! 5. Run Solve (Solver instance passed implicitly via argument or model)
   call model%solve(mol, solver, error, cn, qloc, dcndr, dcndL, dqlocdr, dqlocdL, &
      & energy, gradient, sigma, qvec, dqdr, dqdL)

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
      write(output_unit, '(a)') "[Info] JSON dump written to '"//json_output//"'"
   end if

contains

subroutine help(unit)
   integer, intent(in) :: unit
   write(unit, '(a, *(1x, a))') "Usage: "//prog_name//" [options] <input>"
   write(unit, '(a)') "", "Electronegativity equilibration model", ""
   write(unit, '(2x, a, t35, a)') &
      "-m, --model <model>", "Choose charge model", &
      "-s, --solver <type>", "Solver: 'CG' (Iterative) or 'DIRECT'", &
      "--max-iter <int>", "Max iterations (for CG)", &
      "--tol <real>", "Tolerance (for CG)", &
      "-i, --input <format>", "Input format hint", &
      "-c, --charge <value>", "Molecular charge", &
      "-g, --grad", "Evaluate gradient", &
      "-j, --json", "Output JSON", &
      "-h, --help", "Show help"
   write(unit, '(a)')
end subroutine help

subroutine get_arguments(input, model_id, input_format, grad, charge, &
   & json, solver_choice, cgmiter, cgtol, cgmode, error)

   character(len=:), allocatable, intent(out) :: input
   integer, intent(out) :: model_id
   integer, allocatable, intent(out) :: input_format
   logical, intent(out) :: grad, json
   real(wp), allocatable, intent(out) :: charge
   ! Solver args
   character(len=:), allocatable, intent(out) :: solver_choice
   integer, allocatable, intent(out) :: cgmiter
   real(wp), allocatable, intent(out) :: cgtol
   character(len=32), allocatable, intent(out) :: cgmode
   type(error_type), allocatable, intent(out) :: error

   integer :: iarg, narg, iostat
   character(len=:), allocatable :: arg

   model_id = mchrg_model%eeq2019
   grad = .false.; json = .false.
   iarg = 0; narg = command_argument_count()

   do while(iarg < narg)
      iarg = iarg + 1
      call get_argument(iarg, arg)
      select case(arg)
      case("-h", "-help", "--help")
         call help(output_unit)
         stop
      case("-v", "--version")
         ! (version logic omitted for brevity, similar to before)
         stop
      case default
         if (.not. allocated(input)) then
            call move_alloc(arg, input)
            cycle
         end if
         call fatal_error(error, "Too many positional arguments")
         exit
      case("-m", "--model")
         iarg = iarg + 1; call get_argument(iarg, arg)
         if (.not. allocated(arg)) then; call fatal_error(error, "Missing model"); exit; end if
         if (arg == "eeq2019") then; model_id = mchrg_model%eeq2019
         else if (arg == "eeqbc2025") then; model_id = mchrg_model%eeqbc2025
         else; call fatal_error(error, "Invalid model"); exit; end if
      ! --- Solver Options ---
      case("-s", "--solver")
         iarg = iarg + 1; call get_argument(iarg, solver_choice)
         if (.not. allocated(solver_choice)) then
            call fatal_error(error, "Missing solver type")
            exit
         end if
      case("--max-iter")
         iarg = iarg + 1; call get_argument(iarg, arg)
         allocate(cgmiter)
         read(arg, *, iostat=iostat) cgmiter
         if (iostat /= 0) call fatal_error(error, "Invalid max-iter")
      case("--cgtol")
         iarg = iarg + 1; call get_argument(iarg, arg)
         allocate(cgtol)
         read(arg, *, iostat=iostat) cgtol
         if (iostat /= 0) call fatal_error(error, "Invalid tolerance")
      case("--cgmode")
         iarg = iarg + 1; call get_argument(iarg, arg)
         allocate(cgmode)
         read(arg, *, iostat=iostat) cgmode
         if (iostat /= 0) call fatal_error(error, "Invalid iterative solver mode")
      ! ----------------------
      case("-i", "--input")
         iarg = iarg + 1; call get_argument(iarg, arg)
         input_format = get_filetype("."//arg)
      case("-c", "--charge")
         iarg = iarg + 1; call get_argument(iarg, arg)
         allocate(charge)
         read(arg, *, iostat=iostat) charge
      case("-g", "--grad")
         grad = .true.
      case("-j", "--json")
         json = .true.
      end select
   end do
   
   if (.not. allocated(input) .and. .not. allocated(error)) then
       call help(output_unit); error stop
   end if

end subroutine get_arguments

end program main