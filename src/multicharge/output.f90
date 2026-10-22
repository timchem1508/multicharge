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
   use mctc_env, only : wp, i8, timer_type, format_time
   use mctc_io, only : structure_type
   use mctc_csrlist, only : csr_list
   use mctc_io_convert, only : autoaa
   use mctc_io_constants, only : pi
   use multicharge_model, only : mchrg_model_type
   use multicharge_model_cache, only : mchrg_cache
   use multicharge_model_type, only : mchrg_hessian, hess_index, get_dqdr_row
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
   & press, cache, hessian, list)

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

   !> Second derivative of the energy w.r.t. the Cartesian coordinates in
   !> compressed storage with the factored charge constraint term, used if hess
   !> is absent
   type(mchrg_hessian), intent(in), optional :: hessian

   !> Neighborlist of the compressed derivatives of the partial charges in the
   !> cache
   type(csr_list), intent(in), optional :: list

   integer :: iat, jat, isp, jsp, ic, jc
   real(wp) :: hrow(3), hnorm, dqsum
   real(wp), allocatable :: dqrow(:, :)
   logical :: grad, qgrad, sparse
   character(len=1), parameter :: comp(3) = ["x", "y", "z"]

   grad = present(gradient) .and. present(sigma)
   sparse = .false.
   if (.not. present(dqdr) .and. present(cache) .and. present(list)) &
      & sparse = allocated(cache%dqdrlist_ut)
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
         dqsum = sum(cache%dqdrlist_ut) + sum(cache%dqdrlist_lt) &
            & - sum(cache%uvec) * sum(cache%dlambdadr)
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
            call get_dqdr_row(cache, list, iat, dqrow)
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
   else if (present(hessian)) then
      call write_ascii_hessian_list(unit, mol, hessian)
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


!> Write the Hessian in compressed storage, the blocks on the pattern of its list
!> including the charge constraint term and the factors of this term, which
!> give the blocks outside of the pattern. The rows follow the order of the
!> half list, which only holds the blocks of the atoms y >= x of each row x.
!> The lines of each row are formatted in parallel, as the format of the dense
!> Hessian, and the rows are written in order.
subroutine write_ascii_hessian_list(unit, mol, hessian)

   !> Unit for output
   integer, intent(in) :: unit

   !> Molecular structure data
   class(structure_type), intent(in) :: mol

   !> Second derivative of the energy w.r.t. the Cartesian coordinates in
   !> compressed storage with the factored charge constraint term
   type(mchrg_hessian), intent(in) :: hessian

   integer :: iat, jat, isp, jsp, ic, pos
   integer(i8) :: hat
   real(wp) :: hblk(3, 3), dnorm2, hnorm, usum
   real(wp) :: hrow(3)
   character(len=:), allocatable :: buffer
   character(len=1), parameter :: comp(3) = ["x", "y", "z"]
   ! Length of a line of the block table
   integer, parameter :: hlen = 81

   ! Frobenius norm from the blocks on the pattern, the charge constraint term
   ! of all other blocks is sum(u)^2 |dlambda|^4 without its blocks on the
   ! pattern. The half list holds the block of y and x as the transpose of the
   ! block of x and y.
   usum = hessian%weight
   dnorm2 = sum(hessian%dlambda**2)
   hnorm = usum**2 * dnorm2**2
   !$omp parallel do default(none) schedule(runtime) reduction(+:hnorm) &
   !$omp shared(mol, hessian, usum) private(iat, jat, ic, hat, hblk)
   do iat = 1, mol%nat
      do hat = hessian%list%inl(iat), hessian%list%inl(iat + 1) - 1
         jat = hessian%list%nlat(hat)
         do ic = 1, 3
            hblk(:, ic) = hessian%hesslist(:, ic, hat) &
               & + usum * hessian%dlambda(ic, iat) * hessian%dlambda(:, jat)
         end do
         hnorm = hnorm + sum(hblk**2) - usum**2 &
            & * sum(hessian%dlambda(:, iat)**2) * sum(hessian%dlambda(:, jat)**2)
         if (jat == iat) cycle
         do ic = 1, 3
            hblk(:, ic) = hessian%hesslist(ic, :, hat) &
               & + usum * hessian%dlambda(ic, jat) * hessian%dlambda(:, iat)
         end do
         hnorm = hnorm + sum(hblk**2) - usum**2 &
            & * sum(hessian%dlambda(:, jat)**2) * sum(hessian%dlambda(:, iat)**2)
      end do
   end do

   write(unit, '(a,":", t25, es20.13, 1x, a)') &
   & "Hessian matrix norm", sqrt(max(hnorm, 0.0_wp)), "Eh/a0^2"
   write(unit, '(a, 1x, i0, a, f0.1, a)') "Hessian blocks on the pattern:", &
   & size(hessian%list%nlat, kind=i8), " (", &
   & real(size(hessian%list%nlat, kind=i8), wp) / mol%nat, " per atom)"
   write(unit, '(a)') &
   & "Upper triangle, the block of A and # is the transpose of the block of # and A"
   write(unit, '(78("-"))')
   write(unit, '(a10,1x,a4,3x,a6,1x,a4,3x,a9,1x,*(1x,a12))') &
   & "#", "Z", "#", "A", "component", "d2E/dxdR", "d2E/dydR", "d2E/dzdR"
   write(unit, '(78("-"))')
   !$omp parallel do default(none) ordered schedule(dynamic) &
   !$omp shared(unit, mol, hessian) &
   !$omp private(iat, jat, isp, jsp, ic, hat, pos, hrow, buffer)
   do iat = 1, mol%nat
      isp = mol%id(iat)
      allocate(character(len=3 * (hlen + 1) &
         & * (hessian%list%inl(iat + 1) - hessian%list%inl(iat))) :: buffer)
      pos = 0
      do hat = hessian%list%inl(iat), hessian%list%inl(iat + 1) - 1
         jat = hessian%list%nlat(hat)
         jsp = mol%id(jat)
         do ic = 1, 3
            hrow(:) = hessian%hesslist(:, ic, hat) + hessian%weight &
               & * hessian%dlambda(ic, iat) * hessian%dlambda(:, jat)
            ! (i10,1x,i3,1x,a2,1x,i6,1x,i3,1x,a2,2x,a4,5x,3(2x,es11.3))
            buffer(pos+1:pos+hlen+1) = format_int(iat, 10) // " " &
               & // format_int(mol%num(isp), 3) // " " // mol%sym(isp)(1:2) // " " &
               & // format_int(jat, 6) // " " // format_int(mol%num(jsp), 3) // " " &
               & // mol%sym(jsp)(1:2) // "     " // comp(ic) // "     " &
               & // "  " // format_es(hrow(1)) // "  " // format_es(hrow(2)) &
               & // "  " // format_es(hrow(3)) // new_line('a')
            pos = pos + hlen + 1
         end do
      end do
      ! The record ends with the last line of the row
      !$omp ordered
      write(unit, '(a)') buffer(:pos-1)
      !$omp end ordered
      deallocate(buffer)
   end do
   write(unit, '(78("-"))')
   write(unit, '(a)')

   write(unit, '(a,":")') &
   & "Charge constraint term, outside of the pattern d2E/dR(#)dR(A) = w(#) dlambda(A)^T"
   write(unit, '(62("-"))')
   write(unit, '(a10,1x,a4,3x,a9,1x,*(1x,a12))') "#", "Z", "component", "w", &
   & "dlambda"
   write(unit, '(62("-"))')
   do iat = 1, mol%nat
      isp = mol%id(iat)
      do ic = 1, 3
         write(unit, '(i10,1x,i3,1x,a2,2x,a4,5x,*(2x,es11.3))') &
         & iat, mol%num(isp), mol%sym(isp), comp(ic), &
         & hessian%weight * hessian%dlambda(ic, iat), hessian%dlambda(ic, iat)
      end do
   end do
   write(unit, '(62("-"))')
   write(unit, '(a)')

end subroutine write_ascii_hessian_list


!> Format a real number as the ES11.3 edit descriptor, without the overhead of
!> formatted I/O. Values at an exact decimal tie of the last digit may round to
!> the other side.
function format_es(val) result(str)

   !> Value to format
   real(wp), intent(in) :: val

   !> Formatted value
   character(len=11) :: str

   real(wp) :: aval, mant
   integer :: expo, digits, eabs

   aval = abs(val)
   if (.not. (aval <= huge(1.0_wp))) then
      write(str, '(es11.3)') val
      return
   end if
   if (aval > 0.0_wp) then
      ! Mantissa in [1, 10) rounded to four digits
      expo = floor(log10(aval))
      mant = scale_decimal(aval, expo)
      if (mant < 1.0_wp) then
         expo = expo - 1
         mant = scale_decimal(aval, expo)
      else if (mant >= 10.0_wp) then
         expo = expo + 1
         mant = scale_decimal(aval, expo)
      end if
      digits = nint(mant * 1000.0_wp)
      if (digits >= 10000) then
         digits = digits / 10
         expo = expo + 1
      end if
   else
      expo = 0
      digits = 0
   end if

   str = "  0.000E+00"
   if (sign(1.0_wp, val) < 0.0_wp) str(2:2) = "-"
   str(3:3) = digit(digits / 1000)
   str(5:5) = digit(mod(digits / 100, 10))
   str(6:6) = digit(mod(digits / 10, 10))
   str(7:7) = digit(mod(digits, 10))
   eabs = abs(expo)
   if (eabs < 100) then
      if (expo < 0) str(9:9) = "-"
      str(10:10) = digit(eabs / 10)
      str(11:11) = digit(mod(eabs, 10))
   else
      ! Three-digit exponents drop the exponent letter
      str(8:8) = merge("-", "+", expo < 0)
      str(9:9) = digit(eabs / 100)
      str(10:10) = digit(mod(eabs / 10, 10))
      str(11:11) = digit(mod(eabs, 10))
   end if

end function format_es


!> Divide a positive number by a power of ten without overflow of the power
pure function scale_decimal(aval, expo) result(mant)

   !> Positive number
   real(wp), intent(in) :: aval

   !> Decimal exponent
   integer, intent(in) :: expo

   !> Number divided by 10**expo
   real(wp) :: mant

   if (expo < -300) then
      mant = (aval * 1.0e30_wp) * 10.0_wp**(-expo - 30)
   else
      mant = aval / 10.0_wp**expo
   end if

end function scale_decimal


!> Format an integer as the Iw edit descriptor, asterisks if it does not fit
pure function format_int(val, width) result(str)

   !> Value to format
   integer, intent(in) :: val

   !> Field width
   integer, intent(in) :: width

   !> Formatted value
   character(len=width) :: str

   integer :: rest, pos

   str = repeat(" ", width)
   rest = abs(val)
   pos = width
   do
      if (pos < 1) then
         str = repeat("*", width)
         return
      end if
      str(pos:pos) = digit(mod(rest, 10))
      rest = rest / 10
      pos = pos - 1
      if (rest == 0) exit
   end do
   if (val < 0) then
      if (pos < 1) then
         str = repeat("*", width)
         return
      end if
      str(pos:pos) = "-"
   end if

end function format_int


!> Character of a decimal digit
pure function digit(ival) result(chr)

   !> Digit from 0 to 9
   integer, intent(in) :: ival

   !> Character of the digit
   character(len=1) :: chr

   chr = achar(iachar("0") + ival)

end function digit


!> Write computed properties as a JSON object
subroutine json_results(unit, indentation, energy, gradient, dqdr, charges, cn, cache, &
   & list)

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

   !> Neighborlist of the compressed derivatives of the partial charges in the
   !> cache
   type(csr_list), intent(in), optional :: list

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
   if (.not. present(dqdr) .and. present(cache) .and. present(list)) &
      & sparse = allocated(cache%dqdrlist_ut)
   if (present(dqdr) .or. sparse) then
      write(unit, '(",")', advance='no')
      if (allocated(indent)) write(unit, '(/,a)', advance='no') repeat(indent, 1)
      write(unit, jsonkey, advance='no') 'dq/dr'
      if (sparse) then
         call write_json_dqdr(unit, cache, list, indent)
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
subroutine write_json_dqdr(unit, cache, list, indent)

   !> Unit for output
   integer, intent(in) :: unit

   !> Cache with the compressed derivatives of the partial charges
   type(mchrg_cache), intent(in) :: cache

   !> Neighborlist of the compressed derivatives
   type(csr_list), intent(in) :: list

   !> Indentation string for pretty printing
   character(len=:), allocatable, intent(in) :: indent

   integer :: nat, iat, jat, ic
   real(wp), allocatable :: dqrow(:, :)

   nat = size(cache%dlambdadr, 2)
   allocate(dqrow(3, nat))
   write(unit, '("[")', advance='no')
   do jat = 1, nat
      call get_dqdr_row(cache, list, jat, dqrow)
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
