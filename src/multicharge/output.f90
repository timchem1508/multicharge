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

!> Formatted output routines for the charge model, properties, and results
module multicharge_output
   use mctc_env, only : wp, timer_type, format_time
   use mctc_io, only : structure_type
   use mctc_io_convert, only : autoaa
   use mctc_io_constants, only : pi
   use multicharge_model, only : mchrg_model_type
   use multicharge_model_cache, only : mchrg_cache
   use multicharge_model_type, only : hess_index
   use multicharge_version, only : get_multicharge_version
   implicit none
   private

   public :: write_ascii_model, write_ascii_properties, write_ascii_results, json_results


contains


!> Write the charge model parameters in tabular form
subroutine write_ascii_model(unit, mol, model, verbosity, timer)

   !> Formatted unit
   integer, intent(in) :: unit

   !> Molecular structure data
   class(structure_type), intent(in) :: mol

   !> Electronegativity equilibration model
   class(mchrg_model_type), intent(in) :: model

   !> Verbosity level for output
   integer, intent(in) :: verbosity

   !> Timer for performance measurement
   real(wp), intent(in) :: timer

   integer :: isp
   real(wp), parameter :: sqrt2pi = sqrt(2.0_wp/pi)

   write(unit, '(a, ":")') "Charge model parameter"
   write(unit, '(54("-"))')
   write(unit, '(a4,5x,*(1x,a10))') "Z", "chi/Eh", "kcn_chi/Eh", "eta/Eh", "rad/AA"
   write(unit, '(54("-"))')
   do isp = 1, mol%nid
      write(unit, '(i4, 1x, a4, *(1x,f10.4))') &
      & mol%num(isp), mol%sym(isp), model%chi(isp), model%kcnchi(isp), &
      & model%eta(isp) + sqrt2pi/model%rad(isp), model%rad(isp) * autoaa
   end do
   write(unit, '(54("-"),/)')

   if (verbosity > 1) then
      write(unit, '(a, 1x, a)') "Model setup time :", format_time(timer)
      write(unit, '(a)') ""
   end if

end subroutine write_ascii_model


!> Write coordination numbers and partial charges in tabular form
subroutine write_ascii_properties(unit, mol, model, cn, qvec)

   !> Unit for output
   integer, intent(in) :: unit

   !> Molecular structure data
   class(structure_type), intent(in) :: mol

   !> Electronegativity equilibration model
   class(mchrg_model_type), intent(in) :: model

   !> Coordination numbers
   real(wp), intent(in) :: cn(:)

   !> Atomic partial charges
   real(wp), intent(in) :: qvec(:)

   integer :: iat, isp

   write(unit, '(54("-"))')
   write(unit, '(24x,a)') "Results"
   write(unit, '(54("-"))')
   write(unit, '(a)') ''
   write(unit, '(a,":")') "Electrostatic properties (in atomic units)"
   write(unit, '(54("-"))')
   write(unit, '(a10,1x,a4,5x,*(1x,a10))') "#", "Z", "CN", "q", "chi"
   write(unit, '(54("-"))')
   do iat = 1, mol%nat
      isp = mol%id(iat)
      write(unit, '(i10,1x,i4,1x,a4,*(1x,f10.4))') &
      & iat, mol%num(isp), mol%sym(isp), cn(iat), qvec(iat), &
      & model%chi(isp) - model%kcnchi(isp) * sqrt(cn(iat))
   end do
   write(unit, '(54("-"))')
   write(unit, '(a7,26x,f10.4)') &
   & "Σ", sum(qvec)
   write(unit, '(54("-"),/)')

end subroutine write_ascii_properties


!> Write the electrostatic energy, gradient, charge derivatives, and Hessian
subroutine write_ascii_results(unit, mol, energy, gradient, sigma, dqdr, dqdL, hess, &
   & press, cache)

   !> Unit for output
   integer, intent(in) :: unit

   !> Molecular structure data
   class(structure_type), intent(in) :: mol

   !> Electrostatic energy contributions
   real(wp), intent(in) :: energy(:)

   !> Derivative of the energy w.r.t. the Cartesian coordinates
   real(wp), intent(in), optional :: gradient(:, :)

   !> Derivative of the energy w.r.t. strain deformations
   real(wp), intent(in), optional :: sigma(:, :)

   !> Derivative of the partial charges w.r.t. the Cartesian coordinates
   real(wp), intent(in), optional :: dqdr(:, :, :)

   !> Derivative of the partial charges w.r.t. strain deformations
   real(wp), intent(in), optional :: dqdL(:, :, :)

   !> Second derivative of the energy w.r.t. the Cartesian coordinates in
   !> packed lower-triangle storage
   real(wp), intent(in), optional :: hess(:)

   !> Derivative of the virial w.r.t. positions (3, 3, 3, nat) or strain (3, 3, 3, 3)
   real(wp), intent(in), optional :: press(:, :, :, :)

   !> Cache with the derivatives of the partial charges w.r.t. the Cartesian
   !> coordinates in compressed storage, used if dqdr is absent
   type(mchrg_cache), intent(in), optional :: cache

   integer :: iat, jat, isp, jsp, ic, jc
   real(wp) :: hrow(3), hnorm, dqsum
   real(wp), allocatable :: dqrow(:, :)
   logical :: grad, qgrad, sparse
   character(len=1), parameter :: comp(3) = ["x", "y", "z"]

   grad = present(gradient) .and. present(sigma)
   sparse = .false.
   if (.not. present(dqdr) .and. present(cache)) sparse = allocated(cache%dqdrlist)
   qgrad = present(dqdL) .and. (present(dqdr) .or. sparse)

   write(unit, '(a,":", t25, es20.13, 1x, a)') &
   & "Electrostatic energy", sum(energy), "Eh"
   write(unit, '(a)')

   if (grad) then
      write(unit, '(a,":", t25, es20.13, 1x, a)') &
      & "Energy gradient norm", norm2(gradient), "Eh/a0"
      write(unit, '(54("-"))')
      write(unit, '(a10,1x,a4,5x,*(1x,a10))') "#", "Z", "dE/dx", "dE/dy", "dE/dz"
      write(unit, '(54("-"))')
      do iat = 1, mol%nat
         isp = mol%id(iat)
         write(unit, '(i10,1x,i4,1x,a4,*(es11.3))') &
         & iat, mol%num(isp), mol%sym(isp), gradient(:, iat)
      end do
      write(unit, '(54("-"))')
      write(unit, '(a)')

      write(unit, '(a,":")') &
      & "Energy virial"
      write(unit, '(50("-"))')
      write(unit, '(a15,1x,*(1x,a10))') "component", "x", "y", "z"
      write(unit, '(50("-"))')
      do iat = 1, 3
         write(unit, '(2x,4x,1x,a4,1x,4x,*(es11.3))') &
         & comp(iat), sigma(:, iat)
      end do
      write(unit, '(50("-"))')
      write(unit, '(a)')
   end if

   if (qgrad) then
      ! The compressed derivatives are expanded one charge at a time
      if (sparse) then
         dqsum = sum(cache%dqdrlist) - sum(cache%uvec) * sum(cache%dqdrscal)
      else
         dqsum = sum(dqdr)
      end if
      write(unit, '(a,":", t25, es20.13, 1x, a)') &
      & "Sum over dqdr elements", dqsum, "a.u./a0"
      write(unit, '(72("-"))')
      write(unit, '(a10,1x,a4,3x,a6,1x,a4,3x,*(1x,a12))') "#", "Z", "#", "A", &
      & "dQ(Z)/dx(A)", "dQ(Z)/dy(A)", "dQ(Z)/dz(A)"
      write(unit, '(72("-"))')
      allocate(dqrow(3, mol%nat))
      do iat = 1, mol%nat
         isp = mol%id(iat)
         if (sparse) then
            call cache%get_dqdr_row(iat, dqrow)
         else
            dqrow(:, :) = dqdr(:, :, iat)
         end if
         do jat = 1, mol%nat
            jsp = mol%id(jat)
            write(unit, '(i10,1x,i3,1x,a2,1x,i6,1x,i3,1x,a2,*(2x ,es11.3))') &
            & iat, mol%num(isp), mol%sym(isp), jat, mol%num(jsp), mol%sym(jsp), &
            & dqrow(:, jat)
         end do
      end do
      write(unit, '(72("-"))')
      write(unit, '(a)')

      write(unit, '(a,":")') &
      & "Charge virial"
      write(unit, '(62("-"))')
      write(unit, '(a10,1x,a4,3x,a9,1x,*(1x, a10))')  "#", "Z", "component", "x", "y", "z"
      write(unit, '(62("-"))')
      do iat = 1, mol%nat
         isp = mol%id(iat)
         do jat = 1, 3
            write(unit, '(i10,1x,i3,1x,a2, 2x, a4, 5x,*(es11.3))') &
            & iat, mol%num(isp), mol%sym(isp), comp(jat), dqdL(:, jat, iat)
         end do
      end do
      write(unit, '(62("-"))')
      write(unit, '(a)')
   end if

   if (present(hess)) then
      hnorm = 2.0_wp * sum(hess**2)
      do ic = 1, 3 * mol%nat
         hnorm = hnorm - hess(hess_index(ic, ic))**2
      end do
      write(unit, '(a,":", t25, es20.13, 1x, a)') &
      & "Hessian matrix norm", sqrt(hnorm), "Eh/a0^2"
      write(unit, '(78("-"))')
      write(unit, '(a10,1x,a4,3x,a6,1x,a4,3x,a9,1x,*(1x,a12))') &
      & "#", "Z", "#", "A", "component", "d2E/dxdR", "d2E/dydR", "d2E/dzdR"
      write(unit, '(78("-"))')
      do iat = 1, mol%nat
         isp = mol%id(iat)
         do jat = 1, mol%nat
            jsp = mol%id(jat)
            do ic = 1, 3
               do jc = 1, 3
                  hrow(jc) = hess(hess_index(3 * (jat - 1) + jc, 3 * (iat - 1) + ic))
               end do
               write(unit, '(i10,1x,i3,1x,a2,1x,i6,1x,i3,1x,a2,2x,a4,5x,*(2x,es11.3))') &
               & iat, mol%num(isp), mol%sym(isp), jat, mol%num(jsp), mol%sym(jsp), &
               & comp(ic), hrow
            end do
         end do
      end do
      write(unit, '(78("-"))')
      write(unit, '(a)')
   end if

   if (present(press)) then
      if (size(press, 4) == mol%nat) then
         write(unit, '(a,":")') "Stress gradient (d sigma / d R)"
         write(unit, '(72("-"))')
         write(unit, '(a10,1x,a4,3x,a9,1x,a4,5x,*(1x,a10))') &
         & "#", "Z", "component", "dR", "x", "y", "z"
         write(unit, '(72("-"))')
         do jat = 1, mol%nat
            jsp = mol%id(jat)
            do jc = 1, 3
               do ic = 1, 3
                  write(unit, '(i10,1x,i3,1x,a2,2x,a4,5x,a4,5x,*(es11.3))') &
                  & jat, mol%num(jsp), mol%sym(jsp), comp(ic), comp(jc), press(:, ic, jc, jat)
               end do
            end do
         end do
         write(unit, '(72("-"))')
         write(unit, '(a)')
      else if (size(press, 4) == 3) then
         write(unit, '(a,":")') "Elastic tensor (d sigma / d eps)"
         write(unit, '(50("-"))')
         write(unit, '(a15,1x,a10,1x,*(1x,a10))') "component", "strain", "x", "y", "z"
         write(unit, '(50("-"))')
         do ic = 1, 3
            do jc = 1, 3
               write(unit, '(2x,4x,1x,a4,5x,a4,5x,*(es11.3))') &
               & comp(ic), comp(jc), press(:, jc, ic, 1)
            end do
         end do
         write(unit, '(50("-"))')
         write(unit, '(a)')
      end if
   end if

end subroutine write_ascii_results


!> Write computed properties as a JSON object
subroutine json_results(unit, indentation, energy, gradient, dqdr, charges, cn, cache)

   !> Unit for output
   integer, intent(in) :: unit

   !> Indentation string for pretty printing
   character(len=*), intent(in), optional :: indentation

   !> Electrostatic energy
   real(wp), intent(in), optional :: energy

   !> Derivative of the energy w.r.t. the Cartesian coordinates
   real(wp), intent(in), optional :: gradient(:, :)

   !> Derivative of the partial charges w.r.t. the Cartesian coordinates
   real(wp), intent(in), optional :: dqdr(:, :, :)

   !> Atomic partial charges
   real(wp), intent(in), optional :: charges(:)

   !> Coordination numbers
   real(wp), intent(in), optional :: cn(:)

   !> Cache with the derivatives of the partial charges w.r.t. the Cartesian
   !> coordinates in compressed storage, used if dqdr is absent
   type(mchrg_cache), intent(in), optional :: cache

   character(len=:), allocatable :: indent, version_string
   character(len=*), parameter :: jsonkey = "('""',a,'"":',1x)"
   real(wp), allocatable :: array(:)
   logical :: sparse

   call get_multicharge_version(string=version_string)

   if (present(indentation)) then
      indent = indentation
   else
      indent = ""
   end if

   write(unit, '("{")', advance='no')
   if (allocated(indent)) write(unit, '(/,a)', advance='no') repeat(indent, 1)
   write(unit, jsonkey, advance='no') 'version'
   write(unit, '(1x,a)', advance='no') '"'//version_string//'"'
   if (present(energy)) then
      write(unit, '(",")', advance='no')
      if (allocated(indent)) write(unit, '(/,a)', advance='no') repeat(indent, 1)
      write(unit, jsonkey, advance='no') 'energy'
      write(unit, '(1x,es25.16)', advance='no') energy
   end if
   if (present(gradient)) then
      write(unit, '(",")', advance='no')
      if (allocated(indent)) write(unit, '(/,a)', advance='no') repeat(indent, 1)
      write(unit, jsonkey, advance='no') 'gradient'
      array = reshape(gradient, [size(gradient)])
      call write_json_array(unit, array, indent)
   end if
   sparse = .false.
   if (.not. present(dqdr) .and. present(cache)) sparse = allocated(cache%dqdrlist)
   if (present(dqdr) .or. sparse) then
      write(unit, '(",")', advance='no')
      if (allocated(indent)) write(unit, '(/,a)', advance='no') repeat(indent, 1)
      write(unit, jsonkey, advance='no') 'dq/dr'
      if (sparse) then
         call write_json_dqdr(unit, cache, indent)
      else
         array = reshape(dqdr, [size(dqdr)])
         call write_json_array(unit, array, indent)
      end if
   end if
   if (present(charges)) then
      write(unit, '(",")', advance='no')
      if (allocated(indent)) write(unit, '(/,a)', advance='no') repeat(indent, 1)
      write(unit, jsonkey, advance='no') 'charges'
      array = reshape(charges, [size(charges)])
      call write_json_array(unit, array, indent)
   end if
   if (present(cn)) then
      write(unit, '(",")', advance='no')
      if (allocated(indent)) write(unit, '(/,a)', advance='no') repeat(indent, 1)
      write(unit, jsonkey, advance='no') 'coordination numbers'
      array = reshape(cn, [size(cn)])
      call write_json_array(unit, array, indent)
   end if
   if (allocated(indent)) write(unit, '(/)', advance='no')
   write(unit, '("}")')

end subroutine json_results


!> Write a real array as a JSON array, optionally pretty-printed
subroutine write_json_array(unit, array, indent)

   !> Unit for output
   integer, intent(in) :: unit

   !> Array of values to write
   real(wp), intent(in) :: array(:)

   !> Indentation string for pretty printing
   character(len=:), allocatable, intent(in) :: indent

   integer :: i

   write(unit, '("[")', advance='no')
   do i = 1, size(array)
      call write_json_element(unit, array(i), i == size(array), indent)
   end do
   if (allocated(indent)) write(unit, '(/,a)', advance='no') repeat(indent, 1)
   write(unit, '("]")', advance='no')
end subroutine write_json_array


!> Write the derivatives of the partial charges w.r.t. the Cartesian coordinates
!> as a JSON array in the order of a dense dqdr, expanded from the compressed
!> storage one charge at a time
subroutine write_json_dqdr(unit, cache, indent)

   !> Unit for output
   integer, intent(in) :: unit

   !> Cache with the compressed derivatives of the partial charges
   type(mchrg_cache), intent(in) :: cache

   !> Indentation string for pretty printing
   character(len=:), allocatable, intent(in) :: indent

   integer :: nat, iat, jat, ic
   real(wp), allocatable :: dqrow(:, :)

   nat = size(cache%dqdrscal, 2)
   allocate(dqrow(3, nat))
   write(unit, '("[")', advance='no')
   do jat = 1, nat
      call cache%get_dqdr_row(jat, dqrow)
      do iat = 1, nat
         do ic = 1, 3
            call write_json_element(unit, dqrow(ic, iat), &
               & jat == nat .and. iat == nat .and. ic == 3, indent)
         end do
      end do
   end do
   if (allocated(indent)) write(unit, '(/,a)', advance='no') repeat(indent, 1)
   write(unit, '("]")', advance='no')
end subroutine write_json_dqdr


!> Write one element of a JSON array, followed by a separator unless it is last
subroutine write_json_element(unit, value, last, indent)

   !> Unit for output
   integer, intent(in) :: unit

   !> Value of the element
   real(wp), intent(in) :: value

   !> Whether this is the last element of the array
   logical, intent(in) :: last

   !> Indentation string for pretty printing
   character(len=:), allocatable, intent(in) :: indent

   if (allocated(indent)) write(unit, '(/,a)', advance='no') repeat(indent, 2)
   write(unit, '(es23.16)', advance='no') value
   if (.not. last) write(unit, '(",")', advance='no')
end subroutine write_json_element


end module multicharge_output
