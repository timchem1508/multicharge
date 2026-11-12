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

!> @file multicharge/model/hessian.f90
!> Provides the nuclear Hessians of the charge models

!> Analytical and semi-numerical nuclear Hessians of the charge models
submodule (multicharge_model_type) multicharge_model_hessian
   use mctc_csrlist, only : new_csr_list, distance_screening, sort_neighborlist
   implicit none


contains


!> Compute the molecular nuclear Hessian of the electrostatic energy
module procedure get_hessian_dense

   real(wp), allocatable :: dqdr(:, :, :), dqdL(:, :, :)

   if (any(mol%periodic)) then
      call fatal_error(error, &
         & "Analytical Hessian is not available for periodic systems")
      return
   end if

   ! Get variables from the model cache
   if (.not. allocated(cache%vrhs) .or. .not. allocated(cache%dcndr)) then
      call fatal_error(error, "Charges and CN derivatives are not in the cache, "// &
         & "update the cache with hess and solve first")
      return
   end if

   ! Charge derivatives
   call new_partial_derivs(mol, mol%nat, cache)
   call self%get_xvec_derivs(mol, mol%nat, cache, beta=1.0_wp)
   call self%get_coulomb_derivs(mol, mol%nat, cache, alpha=-1.0_wp, &
      & p=cache%vrhs(:mol%nat))
   allocate(dqdL(3, 3, mol%nat), dqdr(3, mol%nat, mol%nat))
   call get_q_derivs(mol, solver, cache, error, mol%nat, dqdr, dqdL, &
      & unit=unit, verbosity=verbosity)
   if (allocated(error)) return

   ! Second derivatives of 0.5 * q^T A q - q^T X at fixed charges
   hessian(:, :) = 0.0_wp
   call self%get_explicit_hessian(mol, cache, cache%vrhs(:mol%nat), hessian, &
      & error, alpha=0.5_wp, beta=-1.0_wp)
   if (allocated(error)) return
   call self%get_implicit_hessian(mol, cache, hessian, dqdr=dqdr)

end procedure get_hessian_dense


!> Compute the molecular nuclear Hessian of the electrostatic energy on a list
module procedure get_hessian_list

   real(wp) :: tgroup(2)
   type(timer_type) :: timer
   integer :: print_unit, verbosity_hess

   ! Without a list of its own the blocks stay on the complete list
   if (.not. allocated(hessian%list)) hessian%list = list
   call check_hessian_setup(self, mol, solver, cache, hessian, list, error)
   if (allocated(error)) return

   print_unit = output_unit
   if (present(unit)) print_unit = unit
   verbosity_hess = 0
   if (present(verbosity)) verbosity_hess = verbosity

   call timer%push("total")

   ! Partial derivatives dX/dR - dA/dR*q and the inverse on the list
   call print_hessian_header(print_unit, verbosity_hess)
   call timer%push("dabdr")
   call new_partial_derivs(mol, mol%nat, cache, list)
   call self%get_xvec_derivs(mol, mol%nat, cache, beta=1.0_wp, list=list)
   call self%get_coulomb_derivs(mol, mol%nat, cache, alpha=-1.0_wp, &
      & p=cache%vrhs(:mol%nat), list=list)
   call timer%pop
   if (allocated(cache%tpos)) deallocate(cache%tpos)
   call timer%push("inverse")
   select type (solver)
   class is (cg_solver)
      call get_inverse_list(mol, solver, cache, list, error, unit, hessian%list)
   end select
   call timer%pop
   if (allocated(error)) return
   deallocate(cache%alist)

   call new_hessian_blocks(hessian)

   ! Response of the charges, -(dX/dR - dA/dR*q)^T dq/dR
   call timer%push("implicit")
   call self%get_implicit_hessian(mol, cache, hessian, list, times=tgroup)
   call timer%pop
   deallocate(cache%ainvlist, cache%dabdrlist)

   ! Second derivatives of 0.5 * q^T A q - q^T X at fixed charges
   call timer%push("explicit")
   call self%get_explicit_hessian(mol, cache, cache%vrhs(:mol%nat), hessian, &
      & error, list, alpha=0.5_wp, beta=-1.0_wp)
   call timer%pop
   call timer%pop
   if (allocated(error)) return

   if (verbosity_hess > 0) call print_hessian_stages(print_unit, timer, tgroup)

end procedure get_hessian_list


!> Allocate the zero blocks of a compressed Hessian on the upper triangle of its list
module procedure new_hessian_blocks

   if (present(list)) hessian%list = list
   if (allocated(hessian%hesslist)) deallocate(hessian%hesslist)
   allocate(hessian%hesslist(3, 3, size(hessian%list%nlat, kind=i8)), source=0.0_wp)

end procedure new_hessian_blocks


!> Check the setup of the Hessian on a neighborlist
subroutine check_hessian_setup(self, mol, solver, cache, hessian, list, error)

   !> Electronegativity-equilibration model
   class(mchrg_model_type), intent(in) :: self

   !> Molecular structure data
   type(structure_type), intent(in) :: mol

   !> Solver instance
   class(mchrg_solver_type), intent(in) :: solver

   !> Cache updated with the list of the Hessian
   type(mchrg_cache), intent(in) :: cache

   !> Molecular nuclear Hessian with its list
   type(mchrg_hessian), intent(in) :: hessian

   !> Sorted complete neighborlist
   type(csr_list), intent(in) :: list

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   integer :: iat
   integer(i8) :: nup
   logical :: same

   if (any(mol%periodic)) then
      call fatal_error(error, &
         & "Analytical Hessian is not available for periodic systems")
      return
   end if
   if (.not. list%complete) then
      call fatal_error(error, "The Hessian on a neighborlist requires a sorted "// &
         & "complete list, see sort_neighborlist")
      return
   end if
   select type (solver)
   class is (cg_solver)
   class default
      call fatal_error(error, "The Hessian on a neighborlist requires the CG solver")
      return
   end select
   if (.not. allocated(cache%vrhs) .or. .not. allocated(cache%dcndrlist)) then
      call fatal_error(error, "Charges and CN derivatives on the list are not "// &
         & "in the cache, update the cache with hess and solve first")
      return
   end if

   ! The mirrored CN derivatives of the pair arrays hold only without CN cap
   if (self%ncoord%cut > 0.0_wp) then
      call fatal_error(error, "The Hessian on a neighborlist requires the "// &
         & "coordination numbers without cap")
      return
   end if
   if (allocated(self%ncoord_en)) then
      if (self%ncoord_en%cut > 0.0_wp) then
         call fatal_error(error, "The Hessian on a neighborlist requires the "// &
            & "local charges without cap")
         return
      end if
   end if

   ! The cache holds the pair arrays on the rows of the list of the Hessian
   if (size(list%inl) /= mol%nat + 1 .or. .not. allocated(cache%hinl)) then
      call fatal_error(error, "The cache was updated with a different neighborlist")
      return
   end if
   if (size(hessian%list%inl) /= mol%nat + 1 .or. size(cache%hinl) /= mol%nat + 1) then
      call fatal_error(error, "The cache was updated with a different neighborlist")
      return
   end if
   if (any(cache%hinl /= hessian%list%inl) &
      & .or. size(cache%cnmap, kind=i8) /= size(hessian%list%nlat, kind=i8)) &
      & then
      call fatal_error(error, "The cache was updated with a different neighborlist")
      return
   end if

   ! Each row of the list of the Hessian starts with the upper triangle of the
   ! same row of the complete list, a half list holds nothing else
   same = .true.
   !$omp parallel do default(none) schedule(runtime) reduction(.and.:same) &
   !$omp shared(mol, list, hessian) private(iat, nup)
   do iat = 1, mol%nat
      nup = find_lower_start(list, iat) - list%inl(iat)
      if (find_lower_start(hessian%list, iat) - hessian%list%inl(iat) /= nup) then
         same = .false.
      else
         same = same .and. all(hessian%list%nlat(hessian%list%inl(iat): &
            & hessian%list%inl(iat) + nup - 1) &
            & == list%nlat(list%inl(iat):list%inl(iat) + nup - 1))
      end if
   end do
   if (.not. same) then
      call fatal_error(error, "The list of the Hessian is not the upper triangle "// &
         & "of the complete neighborlist")
      return
   end if

end subroutine check_hessian_setup


!> Compute the implicit Hessian of the electrostatic energy
module procedure get_implicit_hessian_dense

   call gemm(reshape(cache%dabdr(:, :, :mol%nat), [3 * mol%nat, mol%nat]), &
      & reshape(dqdr, [3 * mol%nat, mol%nat]), hessian, transb='t', &
      & alpha=-1.0_wp, beta=1.0_wp)

end procedure get_implicit_hessian_dense


!> Compute the implicit Hessian of the electrostatic energy on a list
module procedure get_implicit_hessian_list

   real(wp) :: tgroup(2)
   integer, allocatable :: order(:)

   call get_spatial_order(mol%xyz, order)
   call add_implicit_blocks_list(mol, cache, hessian, list, order, tgroup)
   if (present(times)) times(:) = tgroup

   call get_constraint_factors(mol, cache, list, hessian%dlambda)
   hessian%weight = sum(cache%uvec)
   if (allocated(hessian%wvec)) deallocate(hessian%wvec)

end procedure get_implicit_hessian_list


!> Add the implicit blocks for groups of atoms on a sorted complete list
subroutine add_implicit_blocks_list(mol, cache, hessian, list, order, times, &
   & d2fdq2, d2fdq2_list)

   !> Molecular structure data
   type(structure_type), intent(in) :: mol

   !> Cache with the inverse, partial and mixed derivatives
   type(mchrg_cache), intent(in) :: cache

   !> Molecular nuclear Hessian with allocated blocks on its list, the list of
   !> the update with the inverse on its upper triangle
   type(mchrg_hessian), intent(inout) :: hessian

   !> Sorted complete neighborlist of the partial derivatives
   type(csr_list), intent(in) :: list

   !> Atoms in the spatial order of the groups
   integer, intent(in), contiguous :: order(:)

   !> Times of the group derivatives and their contraction
   real(wp), intent(out) :: times(2)

   !> Second derivatives of an external function w.r.t. the charges
   real(wp), intent(in), contiguous, optional :: d2fdq2(:)

   !> Sorted complete list of the second derivatives w.r.t. the charges
   type(csr_list), intent(in), optional :: d2fdq2_list

   integer :: igrp, ngrp, nbmax, numax, x0, ng, i, a, na, xat, yat, nb, nu, n, u
   integer :: act(max_group)
   integer(i8) :: hat, tick(3), rate
   real(wp) :: scale, tpan, tcon
   real(wp) :: blk(3, 3, max_group)
   logical :: ext
   integer, allocatable :: bloc(:), hloc(:), blist(:), ulist(:)
   integer(i8), allocatable :: hfin(:), ufin_d2f(:), qpos(:, :), jmap(:, :)
   real(wp), allocatable :: pan(:, :), qcol(:, :), pcol(:, :)

   times(:) = 0.0_wp
   ngrp = (mol%nat + max_group - 1) / max_group
   if (ngrp < 1) return
   nbmax = min(mol%nat, max_group * get_max_neighbors(mol%nat, list))
   numax = nbmax
   ext = present(d2fdq2) .and. present(d2fdq2_list)
   scale = -1.0_wp
   if (ext) scale = 1.0_wp

   ! Rows of the blocks and of the inverse
   allocate(hfin(mol%nat))
   do i = 1, mol%nat
      hfin(i) = find_lower_start(hessian%list, i) - 1_i8
   end do
   if (ext) then
      allocate(ufin_d2f(mol%nat))
      do i = 1, mol%nat
         ufin_d2f(i) = find_lower_start(d2fdq2_list, i) - 1_i8
      end do
   end if

   ! bloc and hloc hold the position of each atom in the unions of the group
   tpan = 0.0_wp
   tcon = 0.0_wp
   !$omp parallel default(none) &
   !$omp shared(mol, cache, hessian, list, hfin, ufin_d2f, order, ngrp, nbmax) &
   !$omp shared(numax, ext, scale, d2fdq2, d2fdq2_list) &
   !$omp private(igrp, x0, ng, i, a, na, act, xat, yat, nb, nu, n, u, hat, tick, rate) &
   !$omp private(bloc, hloc, blist, ulist, qpos, jmap, pan, qcol, pcol, blk) &
   !$omp reduction(+:tpan, tcon)
   allocate(bloc(mol%nat), hloc(mol%nat), source=0)
   allocate(blist(nbmax), ulist(numax), qpos(max_group, nbmax), &
      & jmap(max_group, numax), pan(3 * max_group, nbmax), &
      & qcol(3 * max_group, nbmax))
   if (ext) allocate(pcol(3 * max_group, nbmax))
   !$omp do schedule(runtime)
   do igrp = 1, ngrp
      call system_clock(tick(1), rate)
      x0 = max_group * (igrp - 1)
      ng = min(max_group, mol%nat - x0)

      call get_group_dqdr(cache, list, hessian%list, hfin, order(x0 + 1:x0 + ng), &
         & nb, bloc, blist, qpos, pan, qcol)
      if (ext) call get_group_dpdr(cache, hessian%list, hfin, d2fdq2, d2fdq2_list, &
         & ufin_d2f, ng, nb, bloc, blist, qpos, qcol, pan, pcol)
      call system_clock(tick(2))

      nu = 0
      do i = 1, ng
         xat = order(x0 + i)
         do hat = hessian%list%inl(xat), hfin(xat)
            yat = hessian%list%nlat(hat)
            if (hloc(yat) == 0) then
               nu = nu + 1
               hloc(yat) = nu
               ulist(nu) = yat
               jmap(:ng, nu) = 0
            end if
            jmap(i, hloc(yat)) = hat
         end do
      end do

      do u = 1, nu
         yat = ulist(u)
         hloc(yat) = 0
         na = 0
         do i = 1, ng
            hat = jmap(i, u)
            if (hat == 0) cycle
            na = na + 1
            act(na) = i
         end do
         if (ext) then
            call add_column_products_ext(list, cache%dabdrlist, cache%dafdrlist, &
               & yat, bloc, pcol, qcol, na, act, blk)
         else
            call add_column_products(list, cache%dabdrlist, yat, bloc, qcol, &
               & na, act, blk)
         end if
         do a = 1, na
            i = act(a)
            hat = jmap(i, u)
            if (yat == order(x0 + i)) then
               hessian%hesslist(:, :, hat) = hessian%hesslist(:, :, hat) &
                  & + 0.5_wp * scale * (blk(:, :, a) + transpose(blk(:, :, a)))
            else
               hessian%hesslist(:, :, hat) = hessian%hesslist(:, :, hat) &
                  & + scale * blk(:, :, a)
            end if
         end do
      end do
      do n = 1, nb
         bloc(blist(n)) = 0
      end do
      call system_clock(tick(3))
      tpan = tpan + real(tick(2) - tick(1), wp) / real(rate, wp)
      tcon = tcon + real(tick(3) - tick(2), wp) / real(rate, wp)
   end do
   !$omp end do
   deallocate(bloc, hloc, blist, ulist, qpos, jmap, pan, qcol)
   if (ext) deallocate(pcol)
   !$omp end parallel
   times(:) = [tpan, tcon]

end subroutine add_implicit_blocks_list


!> Derivatives of the adjoint vector for a group of atoms
subroutine get_group_dpdr(cache, half, hfin, d2fdq2, d2fdq2_list, ufin_d2f, &
   & ng, nb, bloc, blist, qpos, qcol, pan, pcol)

   !> Cache with the inverse and the mixed derivatives on the list
   type(mchrg_cache), intent(in) :: cache

   !> Sorted list of the update with the inverse on its upper triangle
   type(csr_list), intent(in) :: half

   !> End of the upper triangle in each row of the list of the inverse
   integer(i8), intent(in), contiguous :: hfin(:)

   !> Second derivatives w.r.t. the charges on the upper triangle of their list
   real(wp), intent(in), contiguous :: d2fdq2(:)

   !> Sorted complete list of the second derivatives w.r.t. the charges
   type(csr_list), intent(in) :: d2fdq2_list

   !> End of the upper triangle in each row of d2fdq2_list
   integer(i8), intent(in), contiguous :: ufin_d2f(:)

   !> Number of atoms in the group
   integer, intent(in) :: ng

   !> Number of atoms in the union of the neighbors of the group
   integer, intent(in) :: nb

   !> Index of each atom in the union, zero outside
   integer, intent(in), contiguous :: bloc(:)

   !> Atoms of the union
   integer, intent(in), contiguous :: blist(:)

   !> Position of the entry of each union atom in the row of x
   integer(i8), intent(in), contiguous :: qpos(:, :)

   !> Charge derivatives dq_b/dR_x of the union
   real(wp), intent(in), contiguous :: qcol(:, :)

   !> Workspace for G + d2f/dq2 dq/dR of the union
   real(wp), intent(inout), contiguous :: pan(:, :)

   !> Derivatives dp_b/dR_x of the adjoint vector of the union
   real(wp), intent(inout), contiguous :: pcol(:, :)

   integer :: n, i
   integer(i8) :: kat

   ! G + d2f/dq2 dq/dR on the neighbors of x
   call multiply_group_list(d2fdq2_list, ufin_d2f, d2fdq2, nb, bloc, blist, &
      & qcol, pan)
   do n = 1, nb
      do i = 1, ng
         kat = qpos(i, n)
         if (kat > 0) then
            pan(3 * i - 2:3 * i, n) = pan(3 * i - 2:3 * i, n) &
               & + cache%dafdrlist(:, kat)
         else
            pan(3 * i - 2:3 * i, n) = 0.0_wp
         end if
      end do
   end do
   ! Product with the inverse on the neighbors of x
   call multiply_group_list(half, hfin, cache%ainvlist, nb, bloc, blist, pan, &
      & pcol)
   do n = 1, nb
      do i = 1, ng
         if (qpos(i, n) == 0) pcol(3 * i - 2:3 * i, n) = 0.0_wp
      end do
   end do

end subroutine get_group_dpdr


!> Compute the products of a column of partial derivatives with group derivatives
pure subroutine add_column_products(list, dlist, yat, bloc, col, na, act, blk)

   !> Sorted complete neighborlist
   type(csr_list), intent(in) :: list

   !> Partial derivatives on the list, the row of y holds the column y
   real(wp), intent(in), contiguous :: dlist(:, :)

   !> Atom of the column
   integer, intent(in) :: yat

   !> Index of each atom in the union, zero outside
   integer, intent(in), contiguous :: bloc(:)

   !> Derivatives of the group for the union with leading dimension 3 * max_group
   real(wp), intent(in) :: col(3 * max_group, *)

   !> Number of group atoms with an upper-triangle block with y
   integer, intent(in) :: na

   !> Indices of the active group atoms
   integer, intent(in) :: act(max_group)

   !> 3x3 blocks for the active group atoms
   real(wp), intent(out) :: blk(3, 3, max_group)

   integer :: a, ic, n
   integer(i8) :: kat
   real(wp) :: dy(3)

   blk(:, :, :na) = 0.0_wp
   if (na == max_group) then
      do kat = list%inl(yat), list%inl(yat + 1) - 1
         n = bloc(list%nlat(kat))
         if (n == 0) cycle
         dy(:) = dlist(:, kat)
         do a = 1, max_group
            ic = 3 * (a - 1)
            blk(:, 1, a) = blk(:, 1, a) + dy * col(ic + 1, n)
            blk(:, 2, a) = blk(:, 2, a) + dy * col(ic + 2, n)
            blk(:, 3, a) = blk(:, 3, a) + dy * col(ic + 3, n)
         end do
      end do
   else
      do kat = list%inl(yat), list%inl(yat + 1) - 1
         n = bloc(list%nlat(kat))
         if (n == 0) cycle
         dy(:) = dlist(:, kat)
         do a = 1, na
            ic = 3 * (act(a) - 1)
            blk(:, 1, a) = blk(:, 1, a) + dy * col(ic + 1, n)
            blk(:, 2, a) = blk(:, 2, a) + dy * col(ic + 2, n)
            blk(:, 3, a) = blk(:, 3, a) + dy * col(ic + 3, n)
         end do
      end do
   end if

end subroutine add_column_products


!> Compute the fused products of two columns of partial derivatives with group derivatives
pure subroutine add_column_products_ext(list, dlist1, dlist2, yat, bloc, col1, &
   & col2, na, act, blk)

   !> Sorted complete neighborlist
   type(csr_list), intent(in) :: list

   !> First partial derivatives on the list
   real(wp), intent(in), contiguous :: dlist1(:, :)

   !> Second partial derivatives on the list
   real(wp), intent(in), contiguous :: dlist2(:, :)

   !> Atom of the column
   integer, intent(in) :: yat

   !> Index of each atom in the union, zero outside
   integer, intent(in), contiguous :: bloc(:)

   !> First derivatives of the group for the union
   real(wp), intent(in) :: col1(3 * max_group, *)

   !> Second derivatives of the group for the union
   real(wp), intent(in) :: col2(3 * max_group, *)

   !> Number of group atoms with an upper-triangle block with y
   integer, intent(in) :: na

   !> Indices of the active group atoms
   integer, intent(in) :: act(max_group)

   !> 3x3 blocks for the active group atoms
   real(wp), intent(out) :: blk(3, 3, max_group)

   integer :: a, ic, n
   integer(i8) :: kat
   real(wp) :: dy1(3), dy2(3)

   blk(:, :, :na) = 0.0_wp
   do kat = list%inl(yat), list%inl(yat + 1) - 1
      n = bloc(list%nlat(kat))
      if (n == 0) cycle
      dy1(:) = dlist1(:, kat)
      dy2(:) = dlist2(:, kat)
      do a = 1, na
         ic = 3 * (act(a) - 1)
         blk(:, 1, a) = blk(:, 1, a) + dy1 * col1(ic + 1, n) + dy2 * col2(ic + 1, n)
         blk(:, 2, a) = blk(:, 2, a) + dy1 * col1(ic + 2, n) + dy2 * col2(ic + 2, n)
         blk(:, 3, a) = blk(:, 3, a) + dy1 * col1(ic + 3, n) + dy2 * col2(ic + 3, n)
      end do
   end do

end subroutine add_column_products_ext


!> Compute the nuclear Hessian of an external function of the charges
module procedure get_external_hessian_dense

   real(wp), allocatable :: dqdr(:, :, :), dqdL(:, :, :), padj(:)
   integer :: ndim

   ndim = 3 * mol%nat
   if (size(dfdq) /= mol%nat .or. any(shape(d2fdq2) /= [mol%nat, mol%nat]) &
      & .or. any(shape(dfdrdq) /= [3, mol%nat, mol%nat]) &
      & .or. any(shape(hessian) /= [ndim, ndim])) then
      call fatal_error(error, "External partial derivatives are wrong size")
      return
   end if

   ! Get variables from the model cache
   if (.not. allocated(cache%vrhs) .or. .not. allocated(cache%dcndr)) then
      call fatal_error(error, "Charges and CN derivatives are not in the cache, "// &
         & "update the cache with hess and solve first")
      return
   end if

   ! Charge derivatives
   call new_partial_derivs(mol, mol%nat, cache)
   call self%get_xvec_derivs(mol, mol%nat, cache, beta=1.0_wp)
   call self%get_coulomb_derivs(mol, mol%nat, cache, alpha=-1.0_wp, &
      & p=cache%vrhs(:mol%nat))
   allocate(dqdL(3, 3, mol%nat), dqdr(3, mol%nat, mol%nat))
   call get_q_derivs(mol, solver, cache, error, mol%nat, dqdr, dqdL, &
      & unit=unit, verbosity=verbosity)
   if (allocated(error)) return

   call get_adjoint_vector(mol, solver, cache, error, dfdq, padj, unit=unit, &
      & verbosity=verbosity)
   if (allocated(error)) return

   ! Mixed derivatives of the Lagrangian, G^T = d2f/dRdq - dA/dR*p
   call new_partial_derivs(mol, size(padj), cache, dafdr=.true.)
   call self%get_coulomb_derivs(mol, size(padj), cache, alpha=-1.0_wp, &
      & p=padj(:mol%nat), dafdr=.true.)
   cache%dafdr = cache%dafdr + dfdrdq

   call gemm(reshape(cache%dafdr, [ndim, mol%nat]), &
      & reshape(dqdr, [ndim, mol%nat]), hessian, transb='t')
   call gemm(dqdr, d2fdq2, cache%dafdr, beta=1.0_wp)
   call gemm(reshape(dqdr, [ndim, mol%nat]), &
      & reshape(cache%dafdr, [ndim, mol%nat]), hessian, transb='t', &
      & beta=1.0_wp)
   deallocate(cache%dafdr)

   ! Second derivatives of p^T (X - A*q) at fixed charges and adjoint vector
   call self%get_explicit_hessian(mol, cache, padj(:mol%nat), hessian, &
      & error, alpha=-1.0_wp, beta=1.0_wp)
   if (allocated(error)) return

end procedure get_external_hessian_dense


!> Compute the nuclear Hessian of an external function on a list
module procedure get_external_hessian_list

   integer :: print_unit, verbosity_hess
   real(wp) :: tgroup(2)
   integer, allocatable :: order(:)
   real(wp), allocatable :: padj(:), fuvec(:), tvec(:)
   type(timer_type) :: timer

   ! Without a list of its own the blocks stay on the complete list
   if (.not. allocated(hessian%list)) hessian%list = list
   call check_hessian_setup(self, mol, solver, cache, hessian, list, error)
   if (allocated(error)) return
   call check_external_derivs(mol, dfdq, d2fdq2, d2fdq2_list, dfdrdq, &
      & dfdrdq_list, error)
   if (allocated(error)) return

   print_unit = output_unit
   if (present(unit)) print_unit = unit
   verbosity_hess = 0
   if (present(verbosity)) verbosity_hess = verbosity

   call timer%push("total")

   ! Adjoint vector p = S df/dq and the response of the charge constraint
   call timer%push("adjoint")
   ! The solves with J only read its upper triangle on the list of the update
   call get_adjoint_vector(mol, solver, cache, error, dfdq, padj, hessian%list, &
      & unit, verbosity)
   if (.not. allocated(error)) call get_constraint_response(mol, solver, cache, &
      & hessian%list, d2fdq2, d2fdq2_list, fuvec, tvec, print_unit, error)
   call timer%pop
   if (allocated(error)) return

   call print_hessian_header(print_unit, verbosity_hess)
   call timer%push("dabdr")
   call new_partial_derivs(mol, mol%nat, cache, list)
   call self%get_xvec_derivs(mol, mol%nat, cache, beta=1.0_wp, list=list)
   call self%get_coulomb_derivs(mol, mol%nat, cache, alpha=-1.0_wp, &
      & p=cache%vrhs(:mol%nat), list=list)
   call new_partial_derivs(mol, mol%nat, cache, list, dafdr=.true.)
   call self%get_coulomb_derivs(mol, mol%nat, cache, alpha=-1.0_wp, p=padj, &
      & list=list, dafdr=.true.)
   call add_list_derivs(mol%nat, list, cache%dafdrlist, dfdrdq_list, dfdrdq)
   call timer%pop
   if (allocated(cache%tpos)) deallocate(cache%tpos)
   call timer%push("inverse")
   select type (solver)
   class is (cg_solver)
      call get_inverse_list(mol, solver, cache, list, error, unit, hessian%list)
   end select
   call timer%pop
   if (allocated(error)) return
   deallocate(cache%alist)

   call new_hessian_blocks(hessian)
   call get_external_constraint_factors(mol, cache, list, fuvec, tvec, hessian)

   call timer%push("implicit")
   call get_spatial_order(mol%xyz, order)
   call add_implicit_blocks_list(mol, cache, hessian, list, order, tgroup, &
      & d2fdq2, d2fdq2_list)
   call timer%pop
   deallocate(cache%ainvlist, cache%dabdrlist, cache%dafdrlist)

   ! Second derivatives of p^T (X - A*q) at fixed charges and adjoint vector
   call timer%push("explicit")
   call self%get_explicit_hessian(mol, cache, padj, hessian, error, list, &
      & alpha=-1.0_wp, beta=1.0_wp)
   call timer%pop
   call timer%pop
   if (allocated(error)) return

   if (verbosity_hess > 0) call print_hessian_stages(print_unit, timer, tgroup, &
      & adjoint=.true.)

end procedure get_external_hessian_list


!> Check the derivatives of an external function on their lists
subroutine check_external_derivs(mol, dfdq, d2fdq2, d2fdq2_list, dfdrdq, &
   & dfdrdq_list, error)

   !> Molecular structure data
   type(structure_type), intent(in) :: mol

   !> Derivative of the external function w.r.t. atomic charges
   real(wp), intent(in) :: dfdq(:)

   !> Second derivative w.r.t. atomic charges on the upper triangle of its list
   real(wp), intent(in) :: d2fdq2(:)

   !> Sorted complete list of the second derivative w.r.t. atomic charges
   type(csr_list), intent(in) :: d2fdq2_list

   !> Mixed derivative w.r.t. positions and charges on its list
   real(wp), intent(in) :: dfdrdq(:, :)

   !> Sorted complete list of the mixed derivative w.r.t. positions and charges
   type(csr_list), intent(in) :: dfdrdq_list

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   if (.not. (d2fdq2_list%complete .and. dfdrdq_list%complete)) then
      call fatal_error(error, "The external Hessian on a neighborlist requires "// &
         & "sorted complete lists")
      return
   end if

   if (size(dfdq) /= mol%nat &
      & .or. size(d2fdq2, kind=i8) /= size(d2fdq2_list%nlat, kind=i8) &
      & .or. size(dfdrdq, 1) /= 3 &
      & .or. size(dfdrdq, 2, kind=i8) /= size(dfdrdq_list%nlat, kind=i8)) then
      call fatal_error(error, "External partial derivatives are wrong size")
      return
   end if

end subroutine check_external_derivs


!> Add partial derivatives on another sorted complete list
subroutine add_list_derivs(nat, list, dlist, plist, pval)

   !> Number of atoms
   integer, intent(in) :: nat

   !> Sorted complete neighborlist
   type(csr_list), intent(in) :: list

   !> Partial derivatives on the neighborlist
   real(wp), intent(inout) :: dlist(:, :)

   !> Sorted complete list of the added partial derivatives
   type(csr_list), intent(in) :: plist

   !> Added partial derivatives on their list
   real(wp), intent(in) :: pval(:, :)

   integer :: xat
   integer(i8) :: kat, pat
   integer(i8), allocatable :: loc(:)

   ! loc holds the position of each atom in the row of the neighborlist
   !$omp parallel default(none) shared(nat, list, dlist, plist, pval) &
   !$omp private(xat, kat, pat, loc)
   allocate(loc(nat), source=0_i8)
   !$omp do schedule(runtime)
   do xat = 1, nat
      do kat = list%inl(xat), list%inl(xat + 1) - 1
         loc(list%nlat(kat)) = kat
      end do
      do pat = plist%inl(xat), plist%inl(xat + 1) - 1
         kat = loc(plist%nlat(pat))
         if (kat == 0) cycle
         dlist(:, kat) = dlist(:, kat) + pval(:, pat)
      end do
      do kat = list%inl(xat), list%inl(xat + 1) - 1
         loc(list%nlat(kat)) = 0
      end do
   end do
   !$omp end do
   deallocate(loc)
   !$omp end parallel

end subroutine add_list_derivs


!> Response of the charge constraint to an external function
subroutine get_constraint_response(mol, solver, cache, list, d2fdq2, &
   & d2fdq2_list, fuvec, tvec, unit, error)

   !> Molecular structure data
   type(structure_type), intent(in) :: mol

   !> Solver instance
   class(mchrg_solver_type), intent(in) :: solver

   !> Cache with the coefficient matrix on the list and the constraint response
   type(mchrg_cache), intent(in) :: cache

   !> Sorted list of the update with the coefficient matrix on its upper triangle
   type(csr_list), intent(in) :: list

   !> Second derivative w.r.t. atomic charges on the upper triangle of its list
   real(wp), intent(in), contiguous :: d2fdq2(:)

   !> Sorted complete list of the second derivative w.r.t. atomic charges
   type(csr_list), intent(in) :: d2fdq2_list

   !> Product F u
   real(wp), allocatable, intent(out) :: fuvec(:)

   !> Response t = A^-1 F u
   real(wp), allocatable, intent(out) :: tvec(:)

   !> Output unit
   integer, intent(in) :: unit

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   integer :: iat

   allocate(fuvec(mol%nat), tvec(mol%nat))
   call spsymv_csr(mol%nat, d2fdq2, d2fdq2_list%inl, d2fdq2_list%nlat, cache%uvec, &
      & fuvec, uplo='U')

   ! Diagonal guess of the CG solve
   do iat = 1, mol%nat
      tvec(iat) = fuvec(iat) / (cache%alist(list%inl(iat)) + eps)
   end do
   call solver%solve(alist=cache%alist, xvec=fuvec, vrhs=tvec, list=list, &
      & new_unit=unit, error=error)

end subroutine get_constraint_response


!> Charge constraint factors of the Hessian of an external function
subroutine get_external_constraint_factors(mol, cache, list, fuvec, tvec, &
   & hessian)

   !> Molecular structure data
   type(structure_type), intent(in) :: mol

   !> Cache with the partial and mixed derivatives
   type(mchrg_cache), intent(in) :: cache

   !> Sorted complete neighborlist
   type(csr_list), intent(in) :: list

   !> Product F u of the second derivatives w.r.t. the charges
   real(wp), intent(in) :: fuvec(:)

   !> Response t = A^-1 F u
   real(wp), intent(in) :: tvec(:)

   !> Molecular nuclear Hessian, receives the factors
   type(mchrg_hessian), intent(inout) :: hessian

   call get_constraint_factors(mol, cache, list, hessian%dlambda)
   if (allocated(hessian%wvec)) deallocate(hessian%wvec)
   allocate(hessian%wvec(3, mol%nat), source=0.0_wp)
   call add_transposed_product(list, cache%dafdrlist, cache%uvec, hessian%wvec)
   call add_transposed_product(list, cache%dabdrlist, tvec, hessian%wvec)
   hessian%weight = dot_product(cache%uvec, fuvec)

end subroutine get_external_constraint_factors


!> Semi-numerical Hessian and pressure tensor
module procedure get_numhess_packed

   real(wp), parameter :: step = 1.0e-6_wp
   type(structure_type) :: mol_work
   real(wp), allocatable :: trans(:, :)
   real(wp), allocatable :: g_plus(:, :), g_minus(:, :)
   real(wp) :: s_plus(3, 3), s_minus(3, 3)
   real(wp), allocatable :: xyz_orig(:, :)
   real(wp) :: lattice_orig(3, 3), eps_mat(3, 3)
   integer :: ic, jc, kc, lc, jat
   logical :: periodic

   hess(:) = 0.0_wp
   press(:, :, :, :) = 0.0_wp

   ! Mutable local copy of the structure
   periodic = any(mol%periodic)
   mol_work = mol
   allocate(xyz_orig(3, mol%nat))
   xyz_orig(:, :) = mol%xyz
   lattice_orig(:, :) = 0.0_wp
   if (periodic) lattice_orig(:, :) = mol%lattice

   ! Evaluate the unperturbed system
   call get_lattice_points(mol%periodic, mol%lattice, self%ncoord%cutoff, trans)
   call self%update(mol, cache, trans, .true., list)
   call self%solve(mol, solver, cache, error, qvec=qvec, energy=energy, &
      & gradient=grad, sigma=sigma, list=list, unit=unit, verbosity=verbosity)
   if (allocated(error)) return

   allocate(g_plus(3, mol%nat), g_minus(3, mol%nat))

   ! Cartesian displacements
   do jat = 1, mol%nat
      do jc = 1, 3
         mol_work%xyz(jc, jat) = xyz_orig(jc, jat) + step
         call get_displaced_gradient(self, mol_work, solver, error, g_plus, s_plus, &
            & list, unit, verbosity)
         if (allocated(error)) return

         mol_work%xyz(jc, jat) = xyz_orig(jc, jat) - step
         call get_displaced_gradient(self, mol_work, solver, error, g_minus, s_minus, &
            & list, unit, verbosity)
         if (allocated(error)) return

         mol_work%xyz(jc, jat) = xyz_orig(jc, jat)

         call add_hess_column(3 * (jat - 1) + jc, &
            & 0.5_wp * (g_plus - g_minus) / step, hess)
         if (size(press, 4) == mol%nat) then
            press(:, :, jc, jat) = 0.5_wp * (s_plus - s_minus) / step
         end if
      end do
   end do

   ! Strain deformations
   if (size(press, 4) == 3) then
      eps_mat(:, :) = 0.0_wp
      do ic = 1, 3
         eps_mat(ic, ic) = 1.0_wp
      end do

      do kc = 1, 3
         do lc = 1, 3
            eps_mat(lc, kc) = eps_mat(lc, kc) + step
            mol_work%xyz(:, :) = matmul(eps_mat, xyz_orig)
            if (periodic) mol_work%lattice(:, :) = matmul(eps_mat, lattice_orig)
            call get_displaced_gradient(self, mol_work, solver, error, g_plus, s_plus, &
               & list, unit, verbosity)
            if (allocated(error)) return

            eps_mat(lc, kc) = eps_mat(lc, kc) - 2.0_wp * step
            mol_work%xyz(:, :) = matmul(eps_mat, xyz_orig)
            if (periodic) mol_work%lattice(:, :) = matmul(eps_mat, lattice_orig)
            call get_displaced_gradient(self, mol_work, solver, error, g_minus, s_minus, &
               & list, unit, verbosity)
            if (allocated(error)) return

            eps_mat(lc, kc) = eps_mat(lc, kc) + step
            mol_work%xyz(:, :) = xyz_orig
            if (periodic) mol_work%lattice(:, :) = lattice_orig

            press(:, :, lc, kc) = 0.5_wp * (s_plus - s_minus) / step
         end do
      end do
   end if

end procedure get_numhess_packed


!> Semi-numerical Hessian in dense symmetric storage
module procedure get_numhess_dense

   real(wp), allocatable :: hess_packed(:)
   integer :: idx, jdx, ij, ndim

   hess(:, :) = 0.0_wp

   ndim = 3 * mol%nat
   allocate(hess_packed(ndim * (ndim + 1) / 2))
   call self%get_numhess_packed(mol, solver, cache, error, qvec, energy, grad, sigma, &
      & hess_packed, press, list, unit, verbosity)
   if (allocated(error)) return

   ij = 0
   do idx = 1, ndim
      do jdx = 1, idx
         ij = ij + 1
         hess(jdx, idx) = hess_packed(ij)
         hess(idx, jdx) = hess_packed(ij)
      end do
   end do

end procedure get_numhess_dense


!> Accumulate one column of the Hessian in packed storage
subroutine add_hess_column(jdx, dgdr, hess)

   !> Index of the displaced Cartesian coordinate
   integer, intent(in) :: jdx

   !> Derivative of the energy gradient w.r.t. the displaced coordinate
   real(wp), intent(in) :: dgdr(:, :)

   !> Packed lower-triangle Hessian
   real(wp), intent(inout) :: hess(:)

   integer :: iat, ic, idx, ij

   do iat = 1, size(dgdr, 2)
      do ic = 1, 3
         idx = 3 * (iat - 1) + ic
         ij = hess_index(idx, jdx)
         if (idx == jdx) then
            hess(ij) = dgdr(ic, iat)
         else
            hess(ij) = hess(ij) + 0.5_wp * dgdr(ic, iat)
         end if
      end do
   end do

end subroutine add_hess_column


!> Position of a symmetric-matrix element in packed lower-triangle storage
module procedure hess_index

   ij = max(idx, jdx) * (max(idx, jdx) - 1) / 2 + min(idx, jdx)

end procedure hess_index


!> Energy gradient and virial for a displaced structure
subroutine get_displaced_gradient(self, mol, solver, error, gradient, sigma, &
   & list, unit, verbosity)

   !> Electronegativity-equilibration model
   class(mchrg_model_type), intent(in) :: self

   !> Displaced molecular structure data
   type(structure_type), intent(in) :: mol

   !> The solver instance
   class(mchrg_solver_type), intent(in) :: solver

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   !> Energy gradient
   real(wp), intent(out), contiguous :: gradient(:, :)

   !> Virial
   real(wp), intent(out), contiguous :: sigma(:, :)

   !> Reference neighborlist
   type(csr_list), intent(in), optional :: list

   !> Output unit
   integer, intent(in), optional :: unit

   !> Verbosity level
   integer, intent(in), optional :: verbosity

   type(mchrg_cache) :: cache
   type(csr_list), allocatable :: list_work
   type(distance_screening) :: screen
   real(wp), allocatable :: trans(:, :)

   gradient(:, :) = 0.0_wp
   sigma(:, :) = 0.0_wp

   call get_lattice_points(mol%periodic, mol%lattice, self%ncoord%cutoff, trans)

   if (present(list)) then
      allocate(list_work)
      if (any(mol%periodic)) then
         call new_csr_list(list_work, screen, mol, error, cache%wsc, list%cutoff, &
            & complete=.true.)
         if (allocated(error)) return
         call sort_neighborlist(list_work, cache%wsc)
      else
         call new_csr_list(list_work, screen, mol, error, cutoff=list%cutoff, &
            & complete=.true.)
         if (allocated(error)) return
         call sort_neighborlist(list_work)
      end if
   end if

   call self%update(mol, cache, trans, .true., list_work)
   call self%solve(mol, solver, cache, error, gradient=gradient, sigma=sigma, &
      & list=list_work, unit=unit, verbosity=verbosity)

end subroutine get_displaced_gradient


!> Print header for the analytical Hessian
subroutine print_hessian_header(unit, verbosity)
   !> Output unit
   integer, intent(in) :: unit

   !> Verbosity level
   integer, intent(in) :: verbosity

   if (verbosity > 0) then
      write(unit, '(54("-"))')
      write(unit, '(18x, a)') "Analytical Hessian"
      write(unit, '(54("-"))')
      write(unit, '(a)') ''
   end if
end subroutine print_hessian_header


!> Print the elapsed time of the stages of the Hessian on a neighborlist
subroutine print_hessian_stages(unit, timer, tgroup, adjoint)
   !> Output unit
   integer, intent(in) :: unit

   !> Timer of the stages
   type(timer_type), intent(in) :: timer

   !> Times of the group derivatives and their contraction
   real(wp), intent(in) :: tgroup(2)

   !> Print the stage of the adjoint solves of an external function
   logical, intent(in), optional :: adjoint

   if (present(adjoint)) then
      if (adjoint) call print_stage_time(unit, "Adjoint vector", &
         & timer%get("adjoint"))
   end if
   call print_stage_time(unit, "Partial derivatives dB/dR", timer%get("dabdr"))
   call print_stage_time(unit, "Inverse on the neighborlist", &
      & timer%get("inverse"))
   call print_stage_time(unit, "Implicit Hessian", timer%get("implicit"))
   write(unit, '(a, t34, ":", 2(1x, f0.3))') &
      & "  S*dB/dR/contraction (thread s)", tgroup
   call print_stage_time(unit, "Explicit Hessian", timer%get("explicit"))
   call print_stage_time(unit, "Total Hessian", timer%get("total"))
   write(unit, '(a)') ''
end subroutine print_hessian_stages


end submodule multicharge_model_hessian
