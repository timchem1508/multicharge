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

!> @file multicharge/solver/solver_cache.f90
!> Contains the cache baseclass for the linear equations solvers and a container for mutable cache data

!> Cache for charge models
module solver_type_cache
   use mctc_env, only: wp
   implicit none
   private

   type, public :: cache_container
      !> Mutable data attribute
      class(*), allocatable :: raw
   end type cache_container

   !> Cache for the solvers
   type, abstract, public :: solver_cache

    logical :: cpq
   end type solver_cache

end module solver_type_cache