# 1 "/home/zakharov/Documents/prog/mcharge/multicharge/subprojects/mctc-lib/src/mctc/io/read/qcschema.F90"
# 1 "<built-in>"
# 1 "<command-line>"
# 1 "/home/zakharov/Documents/prog/mcharge/multicharge/cmake_build//"
# 1 "/home/zakharov/Documents/prog/mcharge/multicharge/subprojects/mctc-lib/src/mctc/io/read/qcschema.F90"
! This file is part of mctc-lib.
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


# 1 "/home/zakharov/Documents/prog/mcharge/multicharge/subprojects/mctc-lib/include/mctc/defs.h" 1







# 16 "/home/zakharov/Documents/prog/mcharge/multicharge/subprojects/mctc-lib/src/mctc/io/read/qcschema.F90" 2

module mctc_io_read_qcschema
   use mctc_env_accuracy, only : wp
   use mctc_env_error, only : error_type, fatal_error
   use mctc_io_structure, only : structure_type, new
   use mctc_io_symbols, only : to_number, symbol_length
   use mctc_io_utils, only : to_string





   implicit none
   private

   public :: read_qcschema








contains


subroutine read_qcschema(self, unit, error)

   !> Instance of the molecular structure data
   type(structure_type), intent(out) :: self

   !> File handle
   integer, intent(in) :: unit

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

# 75 "/home/zakharov/Documents/prog/mcharge/multicharge/subprojects/mctc-lib/src/mctc/io/read/qcschema.F90"
   call fatal_error(error, "JSON support not enabled")

end subroutine read_qcschema

# 304 "/home/zakharov/Documents/prog/mcharge/multicharge/subprojects/mctc-lib/src/mctc/io/read/qcschema.F90"


end module mctc_io_read_qcschema
