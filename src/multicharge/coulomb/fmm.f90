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

!> @file multicharge/coulomb/fmm.f90
!> Provides the fast multipole representation of the Coulomb matrix

!> Fast multipole method for Gaussian charges of molecules.
!>
!> The Coulomb matrix is split as A = D + C + K, where D is the diagonal, K the
!> kernel 1/r between all pairs of distinct atoms, and C the short-ranged
!> correction -erfc(gam_ij r)/r of the Gaussian charges, which is stored in
!> the CSR list layout. Only K is evaluated by the fast multipole method, such
!> that all expansions are those of point charges and independent of the
!> species and of the charge model.
!>
!> The expansions use complex solid harmonics with the normalization
!> R_l^m = r^l P_l^m exp(i m phi) / (l+m)!, I_l^m = (l-m)! P_l^m exp(i m phi) / r^(l+1),
!> for which all translations are plain convolutions. The geometry of a
!> product is a sparse uniform octree in Morton order, whose boxes at all levels
!> have the same translation operators in scaled coordinates. The operators of
!> the multipole-to-multipole, multipole-to-local and local-to-local translations
!> are therefore precomputed as dense real matrices and applied to all boxes of
!> one relative position at once by a matrix-matrix product.
!>
!> The expansion order follows from the target accuracy and the number of levels
!> from a cost model of the actual tree.
!>
!> L. Greengard and V. Rokhlin, *J. Comput. Phys.*, **1987**, 73, 325.
!> DOI: [10.1016/0021-9991(87)90140-9](https://doi.org/10.1016/0021-9991(87)90140-9)
!>
!> C. A. White and M. Head-Gordon, *J. Chem. Phys.*, **1994**, 101, 6593.
!> DOI: [10.1063/1.468354](https://doi.org/10.1063/1.468354)
module multicharge_coulomb_fmm
   use mctc_csrlist, only : csr_list, new_csr_list, spsymv_csr
   use mctc_env, only : error_type, fatal_error, sp, wp, i8
   use mctc_io, only : structure_type
   use mctc_io_constants, only : pi
   use multicharge_blas, only : gemm
   use multicharge_coulomb_shortrange, only : get_gamma, get_pair_values, &
      & add_pair_gradient
   use multicharge_coulomb_type, only : coulomb_operator, coulomb_input
   use multicharge_ewald, only : get_dir_cutoff
   !$ use omp_lib, only : omp_get_max_threads
   implicit none
   private

   public :: fmm_operator, new_fmm_operator


   !> Tree level of the fast multipole method
   type :: fmm_level
      !> Number of occupied boxes
      integer :: nbox = 0

      !> Index of the first box of this level in the global expansion arrays
      integer :: base = 0

      !> Morton key of each box
      integer(i8), allocatable :: key(:)

      !> Boxes grouped by the octant in the parent, ranges of the lists below
      integer :: oct_ptr(0:8) = 0

      !> Box of this level, ordered by octant
      integer, allocatable :: oct_box(:)

      !> Parent box of the box in the list above
      integer, allocatable :: oct_par(:)

   end type fmm_level


   !> Fast multipole representation of the Coulomb matrix
   type, extends(coulomb_operator) :: fmm_operator
      !> Order of the multipole expansions
      integer :: order = 10

      !> Number of real coefficients of an expansion
      integer :: ncoef = 0

      !> Finest tree level
      integer :: nlevel = 0

      !> Edge length of the root box
      real(wp) :: edge = 0.0_wp

      !> Cartesian origin of the root box
      real(wp) :: origin(3) = 0.0_wp

      !> Cutoff of the correction for the Gaussian charge distribution
      real(wp) :: cutoff = 0.0_wp

      !> Total number of boxes with expansions
      integer :: ntotal = 0

      !> Original index of each atom in the Morton order
      integer, allocatable :: perm(:)

      !> Cartesian coordinates in the input order
      real(wp), allocatable :: xyz(:, :)

      !> Cartesian coordinates in the Morton order
      real(wp), allocatable :: sx(:), sy(:), sz(:)

      !> Species index of each atom
      integer, allocatable :: id(:)

      !> Pair-wise Gaussian width parameters
      real(wp), allocatable :: gam(:, :)

      !> Diagonal elements of the Coulomb matrix
      real(wp), allocatable :: diag(:)

      !> Levels of the tree from the second to the finest
      type(fmm_level), allocatable :: lev(:)

      !> First sorted atom of each leaf box
      integer, allocatable :: astart(:)

      !> Leaf box and range of the sorted atoms of each work item of the near field
      integer, allocatable :: task_leaf(:), task_a0(:), task_a1(:)

      !> Range of the near field boxes of each leaf box
      integer, allocatable :: near_ptr(:)

      !> Near field boxes of the leaf boxes, including the box itself
      integer, allocatable :: near_box(:)

      !> Compact index of the multipole-to-local matrix of each relative position
      integer, allocatable :: offmap(:)

      !> Range of the interaction pairs of each relative box position, the
      !> pairs of all levels are grouped together as the matrices are shared
      integer, allocatable :: m2l_ptr(:)

      !> Target expansions of the interaction pairs, sorted within each position
      integer, allocatable :: m2l_dst(:)

      !> Source expansions of the interaction pairs
      integer, allocatable :: m2l_src(:)

      !> Multipole-to-multipole matrices of the eight octants
      real(wp), allocatable :: tm2m(:, :, :)

      !> Local-to-local matrices of the eight octants
      real(wp), allocatable :: tl2l(:, :, :)

      !> Multipole-to-local matrices of the relative box positions
      real(wp), allocatable :: tm2l(:, :, :)

      !> Neighbour list of the correction for the Gaussian charge distribution
      type(csr_list) :: list

      !> Diagonal and correction in the CSR list layout
      real(wp), allocatable :: val(:)

      !> Number of near field atom pairs and multipole-to-local translations
      integer(i8) :: npair_near = 0_i8, npair_m2l = 0_i8
   contains
      !> Multiply the Coulomb matrix with a trial vector
      procedure :: multiply
      !> Diagonal of the Coulomb matrix for preconditioning
      procedure :: get_diagonal
      !> Dense block of the Coulomb matrix for preconditioning
      procedure :: get_block
      !> Contract the operator derivatives with two charge-like vectors
      procedure :: add_gradient
      !> Write a summary of the operator setup
      procedure :: write_info
   end type fmm_operator

   real(wp), parameter :: sqrtpi = sqrt(pi)
   real(wp), parameter :: sqrt2pi = sqrt(2.0_wp / pi)

   !> Finest level of the keys of the atoms
   integer, parameter :: key_level = 20

   !> Splitting parameter replacing the Ewald parameter in the real-space
   !> kernels, large enough to remove the smooth kernel entirely
   real(wp), parameter :: alpha_big = 100.0_wp

   !> Largest number of atoms of a work item in the near field
   integer, parameter :: task_size = 32

   !> Largest and smallest number of expansions handled by one matrix product
   integer, parameter :: chunk = 96, min_chunk = 12

   !> Number of distinct relative positions of boxes with the same parent
   !> neighbourhood, offsets run from -3 to 3 along each axis
   integer, parameter :: noff = 343

   !> Estimated cost in nanoseconds of a pair interaction in the near field
   real(wp), parameter :: cost_pair = 0.6_wp

   !> Estimated cost in nanoseconds per coefficient pair of a translation
   real(wp), parameter :: cost_coef = 0.10_wp

   !> Estimated cost in nanoseconds per atom and coefficient for the
   !> assignment of the charges and the evaluation of the fields
   real(wp), parameter :: cost_atom = 4.5_wp


contains


!> Construct the fast multipole representation of the Coulomb matrix
subroutine new_fmm_operator(self, mol, error, input, rad, eta)
   !> Instance of the operator
   type(fmm_operator), intent(out) :: self
   !> Molecular structure data
   type(structure_type), intent(in) :: mol
   !> Error handling
   type(error_type), allocatable, intent(out) :: error
   !> Input configuration
   type(coulomb_input), intent(in) :: input
   !> Species-wise charge widths
   real(wp), intent(in) :: rad(:)
   !> Species-wise chemical hardness
   real(wp), intent(in) :: eta(:)

   integer :: iat
   real(wp) :: tol, rgauss

   if (any(mol%periodic)) then
      call fatal_error(error, "The fast multipole method is only available "//&
         & "for molecules")
      return
   end if

   self%nat = mol%nat
   self%xyz = mol%xyz
   self%id = mol%id
   call get_gamma(rad, self%gam)

   tol = input%accuracy
   if (allocated(input%order)) then
      self%order = input%order
   else
      self%order = default_order(tol)
   end if
   if (self%order < 2 .or. self%order > 40) then
      call fatal_error(error, "Order of the multipole expansion must be "//&
         & "between 2 and 40")
      return
   end if
   self%ncoef = (self%order + 1)**2

   ! The correction decays with the widest Gaussian charge distribution
   rgauss = get_dir_cutoff(minval(self%gam), tol)
   if (allocated(input%cutoff)) then
      self%cutoff = max(input%cutoff, rgauss)
   else
      self%cutoff = rgauss
   end if

   allocate(self%diag(mol%nat))
   do iat = 1, mol%nat
      self%diag(iat) = eta(mol%id(iat)) + sqrt2pi / rad(mol%id(iat))
   end do

   call new_csr_list(self%list, mol, error, cutoff=self%cutoff)
   if (allocated(error)) return
   allocate(self%val(size(self%list%nlat)))
   call get_pair_values(self%list, self%xyz, self%id, self%gam, alpha_big, &
      & self%diag, self%val)

   call setup_tree(self, input)
   call setup_translations(self)

end subroutine new_fmm_operator


!> Default order of the multipole expansion for a target accuracy.
!>
!> The maximal error of a product falls by a factor of 2.3 per order down to
!> 1e-7, which is reached at order 18. Beyond that the pairs of the closest well
!> separated boxes of deep trees dominate and the error falls by only 1.55 per
!> order, which sets the slope for tight accuracies.
pure function default_order(tol) result(order)
   !> Target accuracy
   real(wp), intent(in) :: tol
   !> Order of the multipole expansion
   integer :: order

   real(wp) :: digits

   digits = log10(1.0_wp / tol)
   if (digits <= 7.0_wp) then
      order = ceiling(12.0_wp + (digits - 5.0_wp) / 0.36_wp)
   else
      order = 18 + ceiling((digits - 7.0_wp) / 0.19_wp)
   end if
   order = max(4, min(40, order))

end function default_order


!> Spread the lowest 21 bits of an integer to every third bit
pure function spread_bits(x) result(y)
   !> Integer with at most 21 significant bits
   integer(i8), intent(in) :: x
   !> Spread bits
   integer(i8) :: y

   y = iand(x, int(z'1fffff', i8))
   y = iand(ior(y, ishft(y, 32)), int(z'1f00000000ffff', i8))
   y = iand(ior(y, ishft(y, 16)), int(z'1f0000ff0000ff', i8))
   y = iand(ior(y, ishft(y, 8)), int(z'100f00f00f00f00f', i8))
   y = iand(ior(y, ishft(y, 4)), int(z'10c30c30c30c30c3', i8))
   y = iand(ior(y, ishft(y, 2)), int(z'1249249249249249', i8))

end function spread_bits


!> Collect every third bit of an integer, inverse of the spreading
pure function compact_bits(key) result(x)
   !> Spread bits
   integer(i8), intent(in) :: key
   !> Compacted integer
   integer(i8) :: x

   x = iand(key, int(z'1249249249249249', i8))
   x = iand(ior(x, ishft(x, -2)), int(z'10c30c30c30c30c3', i8))
   x = iand(ior(x, ishft(x, -4)), int(z'100f00f00f00f00f', i8))
   x = iand(ior(x, ishft(x, -8)), int(z'1f0000ff0000ff', i8))
   x = iand(ior(x, ishft(x, -16)), int(z'1f00000000ffff', i8))
   x = iand(ior(x, ishft(x, -32)), int(z'1fffff', i8))

end function compact_bits


!> Morton key of a box from its integer coordinates
pure function encode_key(ix) result(key)
   !> Integer coordinates of the box
   integer, intent(in) :: ix(3)
   !> Morton key
   integer(i8) :: key

   key = ior(spread_bits(int(ix(1), i8)), &
      & ior(ishft(spread_bits(int(ix(2), i8)), 1), &
      & ishft(spread_bits(int(ix(3), i8)), 2)))

end function encode_key


!> Integer coordinates of a box from its Morton key
pure function decode_key(key) result(ix)
   !> Morton key
   integer(i8), intent(in) :: key
   !> Integer coordinates of the box
   integer :: ix(3)

   ix(1) = int(compact_bits(key))
   ix(2) = int(compact_bits(ishft(key, -1)))
   ix(3) = int(compact_bits(ishft(key, -2)))

end function decode_key


!> Index of a key in a sorted list of keys, zero if it is not present
pure function find_key(keys, key) result(idx)
   !> Sorted keys
   integer(i8), intent(in) :: keys(:)
   !> Key to look up
   integer(i8), intent(in) :: key
   !> Index of the key
   integer :: idx

   integer :: lo, hi, mid

   idx = 0
   lo = 1
   hi = size(keys)
   do while (lo <= hi)
      mid = (lo + hi) / 2
      if (keys(mid) < key) then
         lo = mid + 1
      else if (keys(mid) > key) then
         hi = mid - 1
      else
         idx = mid
         return
      end if
   end do

end function find_key


!> Sort keys in ascending order and return the permutation, by a stable
!> bottom-up merge sort
subroutine sort_keys(keys, perm)
   !> Keys to sort, sorted on exit
   integer(i8), intent(inout) :: keys(:)
   !> Permutation from the sorted to the input order
   integer, intent(out) :: perm(:)

   integer :: n, i, width, lo, mid, hi, a, b, k
   integer(i8), allocatable :: tkey(:)
   integer, allocatable :: tperm(:)

   n = size(keys)
   allocate(tkey(n), tperm(n))
   do i = 1, n
      perm(i) = i
   end do

   width = 1
   do while (width < n)
      do lo = 1, n, 2 * width
         mid = min(lo + width - 1, n)
         hi = min(lo + 2 * width - 1, n)
         a = lo
         b = mid + 1
         k = lo
         do while (a <= mid .and. b <= hi)
            if (keys(b) < keys(a)) then
               tkey(k) = keys(b)
               tperm(k) = perm(b)
               b = b + 1
            else
               tkey(k) = keys(a)
               tperm(k) = perm(a)
               a = a + 1
            end if
            k = k + 1
         end do
         do while (a <= mid)
            tkey(k) = keys(a)
            tperm(k) = perm(a)
            a = a + 1
            k = k + 1
         end do
         do while (b <= hi)
            tkey(k) = keys(b)
            tperm(k) = perm(b)
            b = b + 1
            k = k + 1
         end do
      end do
      keys(:) = tkey
      perm(:) = tperm
      width = 2 * width
   end do

end subroutine sort_keys


!> Occupied boxes of one level from the sorted keys of the atoms
subroutine get_boxes(skey, level, ukey, astart)
   !> Sorted keys of the atoms on the finest key level
   integer(i8), intent(in) :: skey(:)
   !> Level of the boxes
   integer, intent(in) :: level
   !> Keys of the occupied boxes
   integer(i8), allocatable, intent(out) :: ukey(:)
   !> First atom of each box
   integer, allocatable, intent(out) :: astart(:)

   integer :: i, nb, shift
   integer(i8) :: last

   shift = 3 * (key_level - level)
   nb = 0
   last = -1_i8
   do i = 1, size(skey)
      if (ishft(skey(i), -shift) /= last) then
         nb = nb + 1
         last = ishft(skey(i), -shift)
      end if
   end do
   allocate(ukey(nb), astart(nb + 1))
   nb = 0
   last = -1_i8
   do i = 1, size(skey)
      if (ishft(skey(i), -shift) /= last) then
         nb = nb + 1
         last = ishft(skey(i), -shift)
         ukey(nb) = last
         astart(nb) = i
      end if
   end do
   astart(nb + 1) = size(skey) + 1

end subroutine get_boxes


!> Occupied boxes of a level that are neighbours of a box, including itself
pure subroutine get_near_boxes(ukey, level, ix, nnear, near)
   !> Keys of the occupied boxes of the level
   integer(i8), intent(in) :: ukey(:)
   !> Level of the boxes
   integer, intent(in) :: level
   !> Integer coordinates of the box
   integer, intent(in) :: ix(3)
   !> Number of occupied neighbours
   integer, intent(out) :: nnear
   !> Indices of the occupied neighbours
   integer, intent(out) :: near(27)

   integer :: dx, dy, dz, jx(3), idx, nmax

   nmax = 2**level
   nnear = 0
   do dz = -1, 1
      jx(3) = ix(3) + dz
      if (jx(3) < 0 .or. jx(3) >= nmax) cycle
      do dy = -1, 1
         jx(2) = ix(2) + dy
         if (jx(2) < 0 .or. jx(2) >= nmax) cycle
         do dx = -1, 1
            jx(1) = ix(1) + dx
            if (jx(1) < 0 .or. jx(1) >= nmax) cycle
            idx = find_key(ukey, encode_key(jx))
            if (idx > 0) then
               nnear = nnear + 1
               near(nnear) = idx
            end if
         end do
      end do
   end do

end subroutine get_near_boxes


!> Occupied boxes that interact with a box by a multipole-to-local translation,
!> the children of the neighbours of the parent that are no neighbours
pure subroutine get_far_boxes(ukey, level, ix, nfar, far, offid)
   !> Keys of the occupied boxes of the level
   integer(i8), intent(in) :: ukey(:)
   !> Level of the boxes
   integer, intent(in) :: level
   !> Integer coordinates of the box
   integer, intent(in) :: ix(3)
   !> Number of interacting boxes
   integer, intent(out) :: nfar
   !> Indices of the interacting boxes
   integer, intent(out) :: far(189)
   !> Relative position of the target to each source box
   integer, intent(out) :: offid(189)

   integer :: cx, cy, cz, jx(3), idx, nmax, lo(3)

   nmax = 2**level
   nfar = 0
   lo = 2 * (ix / 2) - 2
   do cz = lo(3), lo(3) + 5
      if (cz < 0 .or. cz >= nmax) cycle
      do cy = lo(2), lo(2) + 5
         if (cy < 0 .or. cy >= nmax) cycle
         do cx = lo(1), lo(1) + 5
            if (cx < 0 .or. cx >= nmax) cycle
            if (abs(cx - ix(1)) <= 1 .and. abs(cy - ix(2)) <= 1 &
               & .and. abs(cz - ix(3)) <= 1) cycle
            jx = [cx, cy, cz]
            idx = find_key(ukey, encode_key(jx))
            if (idx > 0) then
               nfar = nfar + 1
               far(nfar) = idx
               offid(nfar) = offset_id(ix - jx)
            end if
         end do
      end do
   end do

end subroutine get_far_boxes


!> Index of a relative box position with components from -3 to 3
pure function offset_id(d) result(id)
   !> Relative position
   integer, intent(in) :: d(3)
   !> Index of the position
   integer :: id

   id = (d(1) + 3) + 7 * (d(2) + 3) + 49 * (d(3) + 3) + 1

end function offset_id


!> Relative box position from its index
pure function offset_vec(id) result(d)
   !> Index of the position
   integer, intent(in) :: id
   !> Relative position
   integer :: d(3)

   d(1) = mod(id - 1, 7) - 3
   d(2) = mod((id - 1) / 7, 7) - 3
   d(3) = (id - 1) / 49 - 3

end function offset_vec


!> Estimated time in nanoseconds of a product for a finest level
function get_level_cost(skey, level, ncoef) result(cost)
   !> Sorted keys of the atoms on the finest key level
   integer(i8), intent(in) :: skey(:)
   !> Candidate for the finest level
   integer, intent(in) :: level
   !> Number of coefficients of an expansion
   integer, intent(in) :: ncoef
   !> Estimated cost
   real(wp) :: cost

   integer :: ib, jb, nb, nnear, nfar, near(27), far(189), offid(189), ix(3)
   integer(i8), allocatable :: ukey(:)
   integer, allocatable :: astart(:)
   real(wp) :: pairs, trans

   call get_boxes(skey, level, ukey, astart)
   nb = size(ukey)
   pairs = 0.0_wp
   trans = 0.0_wp
   !$omp parallel do default(none) schedule(static) &
   !$omp shared(ukey, astart, level, nb) private(ib, jb, nnear, nfar, near) &
   !$omp private(far, offid, ix) reduction(+:pairs, trans)
   do ib = 1, nb
      ix = decode_key(ukey(ib))
      call get_near_boxes(ukey, level, ix, nnear, near)
      do jb = 1, nnear
         pairs = pairs + real(astart(ib + 1) - astart(ib), wp) &
            & * real(astart(near(jb) + 1) - astart(near(jb)), wp)
      end do
      if (level >= 2) then
         call get_far_boxes(ukey, level, ix, nfar, far, offid)
         trans = trans + real(nfar, wp)
      end if
   end do

   ! Coarser levels add roughly one seventh of the translations
   cost = cost_pair * pairs + cost_coef * real(ncoef, wp)**2 * trans * 8.0_wp / 7.0_wp &
      & + cost_atom * real(size(skey), wp) * real(ncoef, wp)

end function get_level_cost


!> Build the octree, the interaction lists and the sorted coordinates
subroutine setup_tree(self, input)
   !> Instance of the operator
   type(fmm_operator), intent(inout) :: self
   !> Input configuration
   type(coulomb_input), intent(in) :: input

   integer :: iat, level, lc, best, nlev, ibox, ilev, nb, ix(3), ncell, nused
   integer(i8), allocatable :: skey(:), ukey(:)
   integer, allocatable :: astart(:)
   real(wp) :: lo(3), hi(3), cost, bestcost, scale

   lo = minval(self%xyz, 2)
   hi = maxval(self%xyz, 2)
   self%edge = maxval(hi - lo) * 1.0001_wp + 1.0e-3_wp
   self%origin = 0.5_wp * (lo + hi) - 0.5_wp * self%edge

   allocate(skey(self%nat), self%perm(self%nat))
   ncell = 2**key_level
   scale = real(ncell, wp) / self%edge
   do iat = 1, self%nat
      ix = max(0, min(ncell - 1, int((self%xyz(:, iat) - self%origin) * scale)))
      skey(iat) = encode_key(ix)
   end do
   call sort_keys(skey, self%perm)
   allocate(self%sx(self%nat), self%sy(self%nat), self%sz(self%nat))
   do iat = 1, self%nat
      self%sx(iat) = self%xyz(1, self%perm(iat))
      self%sy(iat) = self%xyz(2, self%perm(iat))
      self%sz(iat) = self%xyz(3, self%perm(iat))
   end do

   ! Choose the finest level
   if (allocated(input%leaf_size)) then
      best = 0
      do level = 0, key_level - 1
         call get_boxes(skey, level, ukey, astart)
         best = level
         if (size(ukey) * input%leaf_size >= self%nat) exit
      end do
   else
      lc = 0
      do level = 0, key_level - 1
         call get_boxes(skey, level, ukey, astart)
         lc = level
         if (size(ukey) * 256 >= self%nat) exit
      end do
      best = lc
      bestcost = huge(1.0_wp)
      do level = max(0, lc - 2), min(key_level - 1, lc + 1)
         cost = get_level_cost(skey, level, self%ncoef)
         if (cost < bestcost) then
            bestcost = cost
            best = level
         end if
      end do
   end if
   self%nlevel = best

   ! Occupied boxes of all levels with expansions, from the finest upwards
   nlev = self%nlevel
   allocate(self%lev(2:max(nlev, 1)))
   if (nlev >= 2) then
      call get_boxes(skey, nlev, ukey, astart)
      call move_alloc(ukey, self%lev(nlev)%key)
      self%lev(nlev)%nbox = size(self%lev(nlev)%key)
      do ilev = nlev - 1, 2, -1
         call coarsen_level(self%lev(ilev + 1)%key, self%lev(ilev)%key)
         self%lev(ilev)%nbox = size(self%lev(ilev)%key)
      end do
      self%ntotal = 0
      do ilev = 2, nlev
         self%lev(ilev)%base = self%ntotal
         self%ntotal = self%ntotal + self%lev(ilev)%nbox
      end do
   end if

   ! Leaf boxes and their near field
   call get_boxes(skey, nlev, ukey, astart)
   nb = size(ukey)
   call move_alloc(astart, self%astart)
   allocate(self%near_ptr(nb + 1))
   !$omp parallel do default(none) schedule(static) shared(ukey, nlev, nb, self) &
   !$omp private(ibox, ix)
   do ibox = 1, nb
      block
         integer :: nnear, near(27)
         ix = decode_key(ukey(ibox))
         call get_near_boxes(ukey, nlev, ix, nnear, near)
         self%near_ptr(ibox) = nnear
      end block
   end do
   block
      integer :: ptr, cnt
      ptr = 1
      do ibox = 1, nb
         cnt = self%near_ptr(ibox)
         self%near_ptr(ibox) = ptr
         ptr = ptr + cnt
      end do
      self%near_ptr(nb + 1) = ptr
      allocate(self%near_box(ptr - 1))
   end block
   !$omp parallel do default(none) schedule(static) shared(ukey, nlev, nb, self) &
   !$omp private(ibox, ix)
   do ibox = 1, nb
      block
         integer :: nnear, near(27)
         ix = decode_key(ukey(ibox))
         call get_near_boxes(ukey, nlev, ix, nnear, near)
         self%near_box(self%near_ptr(ibox):self%near_ptr(ibox) + nnear - 1) &
            & = near(:nnear)
      end block
   end do
   ! Work items of at most task_size atoms of one leaf box balance the near field
   block
      integer :: ntask, ia0, ia1
      ntask = 0
      do ibox = 1, nb
         ntask = ntask + (self%astart(ibox + 1) - self%astart(ibox) + task_size - 1) &
            & / task_size
      end do
      allocate(self%task_leaf(ntask), self%task_a0(ntask), self%task_a1(ntask))
      ntask = 0
      do ibox = 1, nb
         do ia0 = self%astart(ibox), self%astart(ibox + 1) - 1, task_size
            ia1 = min(ia0 + task_size - 1, self%astart(ibox + 1) - 1)
            ntask = ntask + 1
            self%task_leaf(ntask) = ibox
            self%task_a0(ntask) = ia0
            self%task_a1(ntask) = ia1
         end do
      end do
   end block
   self%npair_near = 0_i8
   do ibox = 1, nb
      do iat = self%near_ptr(ibox), self%near_ptr(ibox + 1) - 1
         self%npair_near = self%npair_near &
            & + int(self%astart(ibox + 1) - self%astart(ibox), i8) &
            & * int(self%astart(self%near_box(iat) + 1) &
            & - self%astart(self%near_box(iat)), i8)
      end do
   end do
   if (nlev < 2) then
      allocate(self%offmap(noff), source=0)
      return
   end if

   ! Interaction lists of the translations, merged over all levels
   block
      integer, allocatable :: lptr(:, :), pdst(:), psrc(:)
      integer, allocatable :: cnt(:), pos(:)
      integer :: id, k
      type :: pair_list
         integer, allocatable :: dst(:), src(:)
      end type pair_list
      type(pair_list), allocatable :: pl(:)

      allocate(lptr(noff + 1, 2:nlev), pl(2:nlev), cnt(noff), pos(noff))
      cnt(:) = 0
      do ilev = 2, nlev
         call build_far_lists(self%lev(ilev), ilev, pdst, pl(ilev)%dst, pl(ilev)%src)
         lptr(:, ilev) = pdst
         cnt(:) = cnt(:) + (lptr(2:, ilev) - lptr(:noff, ilev))
         if (ilev > 2) call build_octant_lists(self%lev(ilev), self%lev(ilev - 1))
      end do
      allocate(self%m2l_ptr(noff + 1), self%offmap(noff), source=0)
      self%m2l_ptr(1) = 1
      nused = 0
      do id = 1, noff
         self%m2l_ptr(id + 1) = self%m2l_ptr(id) + cnt(id)
         if (cnt(id) > 0) then
            nused = nused + 1
            self%offmap(id) = nused
         end if
      end do
      allocate(self%m2l_dst(self%m2l_ptr(noff + 1) - 1), &
         & self%m2l_src(self%m2l_ptr(noff + 1) - 1))
      pos(:) = self%m2l_ptr(:noff)
      do ilev = 2, nlev
         do id = 1, noff
            k = lptr(id + 1, ilev) - lptr(id, ilev)
            if (k == 0) cycle
            self%m2l_dst(pos(id):pos(id) + k - 1) = self%lev(ilev)%base &
               & + pl(ilev)%dst(lptr(id, ilev):lptr(id + 1, ilev) - 1)
            self%m2l_src(pos(id):pos(id) + k - 1) = self%lev(ilev)%base &
               & + pl(ilev)%src(lptr(id, ilev):lptr(id + 1, ilev) - 1)
            pos(id) = pos(id) + k
         end do
      end do
      self%npair_m2l = size(self%m2l_dst)
   end block

end subroutine setup_tree


!> Keys of the parent boxes from the sorted keys of the boxes
subroutine coarsen_level(ckey, pkey)
   !> Keys of the child boxes
   integer(i8), intent(in) :: ckey(:)
   !> Keys of the parent boxes
   integer(i8), allocatable, intent(out) :: pkey(:)

   integer :: i, n
   integer(i8) :: last

   n = 0
   last = -1_i8
   do i = 1, size(ckey)
      if (ishft(ckey(i), -3) /= last) then
         n = n + 1
         last = ishft(ckey(i), -3)
      end if
   end do
   allocate(pkey(n))
   n = 0
   last = -1_i8
   do i = 1, size(ckey)
      if (ishft(ckey(i), -3) /= last) then
         n = n + 1
         last = ishft(ckey(i), -3)
         pkey(n) = last
      end if
   end do

end subroutine coarsen_level


!> Interaction lists of one level, grouped by the relative position
subroutine build_far_lists(lv, level, ptr, dst, src)
   !> Level of the tree
   type(fmm_level), intent(in) :: lv
   !> Number of the level
   integer, intent(in) :: level
   !> Range of the pairs of each relative position
   integer, allocatable, intent(out) :: ptr(:)
   !> Target boxes of the pairs
   integer, allocatable, intent(out) :: dst(:)
   !> Source boxes of the pairs
   integer, allocatable, intent(out) :: src(:)

   integer :: ib, k, nfar, far(189), offid(189), ix(3), npair, id
   integer, allocatable :: bptr(:), bsrc(:), boff(:), cnt(:), pos(:)

   allocate(bptr(lv%nbox + 1))
   !$omp parallel do default(none) schedule(static) shared(lv, level, bptr) &
   !$omp private(ib, ix, nfar, far, offid)
   do ib = 1, lv%nbox
      ix = decode_key(lv%key(ib))
      call get_far_boxes(lv%key, level, ix, nfar, far, offid)
      bptr(ib) = nfar
   end do
   npair = 1
   do ib = 1, lv%nbox
      k = bptr(ib)
      bptr(ib) = npair
      npair = npair + k
   end do
   bptr(lv%nbox + 1) = npair
   npair = npair - 1
   allocate(bsrc(npair), boff(npair))
   !$omp parallel do default(none) schedule(static) shared(lv, level, bptr, bsrc, boff) &
   !$omp private(ib, ix, nfar, far, offid)
   do ib = 1, lv%nbox
      ix = decode_key(lv%key(ib))
      call get_far_boxes(lv%key, level, ix, nfar, far, offid)
      bsrc(bptr(ib):bptr(ib) + nfar - 1) = far(:nfar)
      boff(bptr(ib):bptr(ib) + nfar - 1) = offid(:nfar)
   end do

   ! Group the pairs by their relative position
   allocate(cnt(noff), pos(noff), source=0)
   do k = 1, npair
      cnt(boff(k)) = cnt(boff(k)) + 1
   end do
   allocate(ptr(noff + 1))
   ptr(1) = 1
   do id = 1, noff
      ptr(id + 1) = ptr(id) + cnt(id)
   end do
   allocate(dst(npair), src(npair))
   pos(:) = ptr(:noff)
   do ib = 1, lv%nbox
      do k = bptr(ib), bptr(ib + 1) - 1
         dst(pos(boff(k))) = ib
         src(pos(boff(k))) = bsrc(k)
         pos(boff(k)) = pos(boff(k)) + 1
      end do
   end do

end subroutine build_far_lists


!> Group the boxes of a level by their octant in the parent
subroutine build_octant_lists(lv, parent)
   !> Level of the tree
   type(fmm_level), intent(inout) :: lv
   !> Parent level of the tree
   type(fmm_level), intent(in) :: parent

   integer :: ib, o, pos(0:7)

   lv%oct_ptr(:) = 0
   do ib = 1, lv%nbox
      o = int(iand(lv%key(ib), 7_i8))
      lv%oct_ptr(o + 1) = lv%oct_ptr(o + 1) + 1
   end do
   lv%oct_ptr(0) = 1
   do o = 0, 7
      lv%oct_ptr(o + 1) = lv%oct_ptr(o) + lv%oct_ptr(o + 1)
   end do
   allocate(lv%oct_box(lv%nbox), lv%oct_par(lv%nbox))
   pos(:) = lv%oct_ptr(0:7)
   do ib = 1, lv%nbox
      o = int(iand(lv%key(ib), 7_i8))
      lv%oct_box(pos(o)) = ib
      lv%oct_par(pos(o)) = find_key(parent%key, ishft(lv%key(ib), -3))
      pos(o) = pos(o) + 1
   end do

end subroutine build_octant_lists


!> Value of a table of regular or irregular harmonics for any sign of m
pure function hfull(tab, l, m) result(val)
   !> Table of harmonics for non-negative m
   complex(wp), intent(in) :: tab(0:, 0:)
   !> Degree
   integer, intent(in) :: l
   !> Order
   integer, intent(in) :: m
   !> Harmonic
   complex(wp) :: val

   if (abs(m) > l) then
      val = (0.0_wp, 0.0_wp)
   else if (m >= 0) then
      val = tab(l, m)
   else
      val = conjg(tab(l, -m))
      if (mod(-m, 2) /= 0) val = -val
   end if

end function hfull


!> Regular solid harmonics R_l^m(r) = r^l P_l^m exp(i m phi) / (l+m)!
pure subroutine regular_harmonics(vec, p, rt)
   !> Cartesian position
   real(wp), intent(in) :: vec(3)
   !> Maximal degree
   integer, intent(in) :: p
   !> Harmonics for non-negative orders
   complex(wp), intent(out) :: rt(0:, 0:)

   integer :: l, m
   real(wp) :: r2
   complex(wp) :: xy

   r2 = vec(1)**2 + vec(2)**2 + vec(3)**2
   xy = cmplx(vec(1), vec(2), wp)
   rt(0, 0) = (1.0_wp, 0.0_wp)
   do m = 0, p
      if (m > 0) rt(m, m) = -xy / real(2 * m, wp) * rt(m - 1, m - 1)
      if (m < p) rt(m + 1, m) = vec(3) * rt(m, m)
      do l = m + 2, p
         rt(l, m) = (real(2 * l - 1, wp) * vec(3) * rt(l - 1, m) &
            & - r2 * rt(l - 2, m)) / real((l + m) * (l - m), wp)
      end do
   end do

end subroutine regular_harmonics


!> Irregular solid harmonics I_l^m(r) = (l-m)! P_l^m exp(i m phi) / r^(l+1)
pure subroutine irregular_harmonics(vec, p, it)
   !> Cartesian position
   real(wp), intent(in) :: vec(3)
   !> Maximal degree
   integer, intent(in) :: p
   !> Harmonics for non-negative orders
   complex(wp), intent(out) :: it(0:, 0:)

   integer :: l, m
   real(wp) :: r2
   complex(wp) :: xy

   r2 = vec(1)**2 + vec(2)**2 + vec(3)**2
   xy = cmplx(vec(1), vec(2), wp)
   it(0, 0) = cmplx(1.0_wp / sqrt(r2), 0.0_wp, wp)
   do m = 0, p
      if (m > 0) it(m, m) = -real(2 * m - 1, wp) * xy / r2 * it(m - 1, m - 1)
      if (m < p) it(m + 1, m) = real(2 * m + 1, wp) * vec(3) / r2 * it(m, m)
      do l = m + 2, p
         it(l, m) = (real(2 * l - 1, wp) * vec(3) * it(l - 1, m) &
            & - real((l - 1 + m) * (l - 1 - m), wp) * it(l - 2, m)) / r2
      end do
   end do

end subroutine irregular_harmonics


!> Index of the real coefficient of degree l and order m, the imaginary part
!> of an order above zero is stored behind the real part
pure function coef_index(l, m) result(idx)
   !> Degree
   integer, intent(in) :: l
   !> Non-negative order
   integer, intent(in) :: m
   !> Index of the real part
   integer :: idx

   if (m == 0) then
      idx = l * l + 1
   else
      idx = l * l + 2 * m
   end if

end function coef_index


!> Coefficient of a translation between the complex expansions
pure function translation_coef(mode, lo, mo, li, mi, tab) result(val)
   !> Kind of translation, multipole-to-local (1), multipole-to-multipole (2)
   !> or local-to-local (3)
   integer, intent(in) :: mode
   !> Degree and order of the output coefficient
   integer, intent(in) :: lo, mo
   !> Degree and order of the input coefficient
   integer, intent(in) :: li, mi
   !> Irregular harmonics of the shift for the first, regular harmonics otherwise
   complex(wp), intent(in) :: tab(0:, 0:)
   !> Coefficient of the translation
   complex(wp) :: val

   integer :: l, m

   val = (0.0_wp, 0.0_wp)
   select case(mode)
   case(1)
      val = hfull(tab, li + lo, mi - mo)
      if (mod(lo + mo, 2) /= 0) val = -val
   case(2)
      if (li > lo) return
      l = lo - li
      m = mo - mi
      if (abs(m) > l) return
      val = conjg(hfull(tab, l, m)) * 0.5_wp**li
   case(3)
      if (li < lo) return
      l = li - lo
      m = mi - mo
      if (abs(m) > l) return
      val = hfull(tab, l, m) * 0.5_wp**(lo + 1)
   end select

end function translation_coef


!> Real matrix of a translation acting on the real coefficients of the
!> expansions. Coefficients of negative order follow from the symmetry of a
!> real field, which maps the matrix of the complex translation to the
!> real and imaginary parts of the non-negative orders.
subroutine build_matrix(p, mode, tab, tmat)
   !> Order of the expansions
   integer, intent(in) :: p
   !> Kind of translation
   integer, intent(in) :: mode
   !> Table of harmonics of the shift
   complex(wp), intent(in) :: tab(0:, 0:)
   !> Matrix of the translation
   real(wp), intent(out) :: tmat(:, :)

   integer :: lo, mo, li, mi, ro, ci
   real(wp) :: s
   complex(wp) :: c1, c2

   tmat(:, :) = 0.0_wp
   do lo = 0, p
      do mo = 0, lo
         ro = coef_index(lo, mo)
         do li = 0, p
            do mi = 0, li
               ci = coef_index(li, mi)
               c1 = translation_coef(mode, lo, mo, li, mi, tab)
               if (mi == 0) then
                  tmat(ro, ci) = tmat(ro, ci) + real(c1, wp)
                  if (mo > 0) tmat(ro + 1, ci) = tmat(ro + 1, ci) + aimag(c1)
               else
                  c2 = translation_coef(mode, lo, mo, li, -mi, tab)
                  s = merge(-1.0_wp, 1.0_wp, mod(mi, 2) /= 0)
                  tmat(ro, ci) = tmat(ro, ci) + real(c1, wp) + s * real(c2, wp)
                  tmat(ro, ci + 1) = tmat(ro, ci + 1) - aimag(c1) + s * aimag(c2)
                  if (mo > 0) then
                     tmat(ro + 1, ci) = tmat(ro + 1, ci) + aimag(c1) + s * aimag(c2)
                     tmat(ro + 1, ci + 1) = tmat(ro + 1, ci + 1) &
                        & + real(c1, wp) - s * real(c2, wp)
                  end if
               end if
            end do
         end do
      end do
   end do

end subroutine build_matrix


!> Precompute the translation matrices that occur in the tree
subroutine setup_translations(self)
   !> Instance of the operator
   type(fmm_operator), intent(inout) :: self

   integer :: p, nc, o, id, k, nused
   integer, allocatable :: ids(:)
   real(wp) :: delta(3)
   complex(wp), allocatable :: tab(:, :), tabi(:, :)

   if (self%nlevel < 2) return
   p = self%order
   nc = self%ncoef

   allocate(self%tm2m(nc, nc, 0:7), self%tl2l(nc, nc, 0:7))
   allocate(tab(0:p, 0:p))
   do o = 0, 7
      delta = 0.5_wp * (real([iand(o, 1), iand(ishft(o, -1), 1), &
         & iand(ishft(o, -2), 1)], wp) - 0.5_wp)
      call regular_harmonics(delta, p, tab)
      call build_matrix(p, 2, tab, self%tm2m(:, :, o))
      call build_matrix(p, 3, tab, self%tl2l(:, :, o))
   end do

   deallocate(tab)

   nused = maxval(self%offmap)
   allocate(self%tm2l(nc, nc, nused), ids(nused))
   do id = 1, noff
      if (self%offmap(id) > 0) ids(self%offmap(id)) = id
   end do
   !$omp parallel default(none) shared(self, ids, nused, p) private(k, tabi, delta)
   allocate(tabi(0:2 * p, 0:2 * p))
   !$omp do schedule(dynamic)
   do k = 1, nused
      delta = real(offset_vec(ids(k)), wp)
      call irregular_harmonics(delta, 2 * p, tabi)
      call build_matrix(p, 1, tabi, self%tm2l(:, :, k))
   end do
   !$omp end do
   deallocate(tabi)
   !$omp end parallel

end subroutine setup_translations


!> Apply a translation matrix to selected expansions and add the results:
!> y(:, dst(k)) = y(:, dst(k)) + T x(:, src(k)). Destinations must be unique.
!> Without y the results are added to the expansions x themselves, the sources
!> and the destinations must then not overlap.
subroutine apply_translation(tmat, src, dst, x, y)
   !> Translation matrix
   real(wp), intent(in) :: tmat(:, :)
   !> Source expansions
   integer, intent(in) :: src(:)
   !> Destination expansions
   integer, intent(in) :: dst(:)
   !> Expansions the translation acts on
   real(wp), intent(inout), target :: x(:, :)
   !> Expansions the results are added to
   real(wp), intent(inout), target, optional :: y(:, :)

   integer :: n, nc, ic, jc, m, k, cs, nthr
   real(wp), allocatable :: xb(:, :), yb(:, :)
   real(wp), pointer :: out(:, :)

   n = size(src)
   nc = size(tmat, 1)
   if (n == 0) return

   ! Spread few pairs over all threads, many pairs favour large products
   nthr = 1
   !$ nthr = omp_get_max_threads()
   cs = max(min_chunk, min(chunk, (n + 2 * nthr - 1) / (2 * nthr)))

   if (present(y)) then
      out => y
   else
      out => x
   end if

   !$omp parallel default(none) shared(tmat, src, dst, x, out, n, nc, cs) &
   !$omp private(ic, jc, m, k, xb, yb)
   allocate(xb(nc, cs), yb(nc, cs))
   !$omp do schedule(dynamic)
   do ic = 1, n, cs
      m = min(cs, n - ic + 1)
      do k = 1, m
         xb(:, k) = x(:, src(ic + k - 1))
      end do
      call gemm(tmat, xb(:, :m), yb(:, :m))
      do k = 1, m
         jc = dst(ic + k - 1)
         out(:, jc) = out(:, jc) + yb(:, k)
      end do
   end do
   !$omp end do
   deallocate(xb, yb)
   !$omp end parallel

end subroutine apply_translation


!> Multipole-to-local translations of all interaction pairs. The work items
!> are chunks of the pairs of one relative position, which are spread over
!> all threads together, such that few pairs per position do not limit the
!> parallelism. The results of a chunk are added to the shared local expansions
!> in a short critical section.
subroutine apply_m2l(self, mult, loc)
   !> Instance of the operator
   type(fmm_operator), intent(in) :: self
   !> Multipole expansions
   real(wp), intent(in) :: mult(:, :)
   !> Local expansions
   real(wp), intent(inout) :: loc(:, :)

   integer :: id, it, ntask, nc, cs, nthr, i0, i1, m, k, jc, npair
   integer, allocatable :: t_id(:), t_i0(:), t_i1(:)
   real(wp), allocatable :: xb(:, :), yb(:, :)

   nc = self%ncoef
   npair = size(self%m2l_dst)
   if (npair == 0) return

   nthr = 1
   !$ nthr = omp_get_max_threads()
   cs = max(min_chunk, min(chunk, npair / (8 * nthr)))

   ntask = 0
   do id = 1, noff
      ntask = ntask + (self%m2l_ptr(id + 1) - self%m2l_ptr(id) + cs - 1) / cs
   end do
   allocate(t_id(ntask), t_i0(ntask), t_i1(ntask))
   ntask = 0
   do id = 1, noff
      do i0 = self%m2l_ptr(id), self%m2l_ptr(id + 1) - 1, cs
         ntask = ntask + 1
         t_id(ntask) = id
         t_i0(ntask) = i0
         t_i1(ntask) = min(i0 + cs - 1, self%m2l_ptr(id + 1) - 1)
      end do
   end do

   !$omp parallel default(none) &
   !$omp shared(self, mult, loc, ntask, nc, cs, t_id, t_i0, t_i1) &
   !$omp private(it, i0, i1, m, k, jc, id, xb, yb)
   allocate(xb(nc, cs), yb(nc, cs))
   !$omp do schedule(dynamic)
   do it = 1, ntask
      id = t_id(it)
      i0 = t_i0(it)
      i1 = t_i1(it)
      m = i1 - i0 + 1
      do k = 1, m
         xb(:, k) = mult(:, self%m2l_src(i0 + k - 1))
      end do
      call gemm(self%tm2l(:, :, self%offmap(id)), xb(:, :m), yb(:, :m))
      !$omp critical (m2l_scatter)
      do k = 1, m
         jc = self%m2l_dst(i0 + k - 1)
         loc(:, jc) = loc(:, jc) + yb(:, k)
      end do
      !$omp end critical (m2l_scatter)
   end do
   !$omp end do
   deallocate(xb, yb)
   !$omp end parallel

end subroutine apply_m2l


!> Center of a box and its edge length
pure subroutine box_geometry(self, level, key, center, h)
   !> Instance of the operator
   type(fmm_operator), intent(in) :: self
   !> Level of the box
   integer, intent(in) :: level
   !> Morton key of the box
   integer(i8), intent(in) :: key
   !> Center of the box
   real(wp), intent(out) :: center(3)
   !> Edge length of the box
   real(wp), intent(out) :: h

   integer :: ix(3)

   h = self%edge / real(2**level, wp)
   ix = decode_key(key)
   center = self%origin + (real(ix, wp) + 0.5_wp) * h

end subroutine box_geometry


!> Potential and field of the point charges 1/r in the Morton order. The
!> potential is added to phi and the field to grad, either may be omitted.
subroutine far_field(self, qs, phi, grad)
   !> Instance of the operator
   class(fmm_operator), intent(in), target :: self
   !> Charges in the Morton order
   real(wp), intent(in) :: qs(:)
   !> Potential at the atoms in the Morton order
   real(wp), intent(inout), optional :: phi(:)
   !> Gradient of the potential at the atoms in the Morton order
   real(wp), intent(inout), optional :: grad(:, :)

   integer :: nlev, nc, p, ilev, ib, ia, o, id, k, nleaf
   real(wp) :: center(3), h, u(3)
   real(wp), allocatable :: mult(:, :), loc(:, :)
   complex(wp), allocatable :: rt(:, :)
   type(fmm_level), pointer :: lv, plv

   nlev = self%nlevel
   p = self%order
   nc = self%ncoef

   if (nlev < 2) then
      call near_field(self, qs, phi, grad)
      return
   end if
   nleaf = self%lev(nlev)%nbox

   allocate(mult(nc, self%ntotal), loc(nc, self%ntotal))
   !$omp parallel do default(none) schedule(static) shared(mult, loc, self, nc)
   do ib = 1, self%ntotal
      mult(:, ib) = 0.0_wp
      loc(:, ib) = 0.0_wp
   end do

   ! Multipole expansions of the leaf boxes
   lv => self%lev(nlev)
   !$omp parallel default(none) &
   !$omp shared(self, lv, qs, mult, nlev, nleaf, p) &
   !$omp private(ib, ia, center, h, u, rt)
   allocate(rt(0:p, 0:p))
   !$omp do schedule(dynamic, 4)
   do ib = 1, nleaf
      call box_geometry(self, nlev, lv%key(ib), center, h)
      do ia = self%astart(ib), self%astart(ib + 1) - 1
         u = ([self%sx(ia), self%sy(ia), self%sz(ia)] - center) / h
         call regular_harmonics(u, p, rt)
         call add_moment(qs(ia), rt, p, mult(:, lv%base + ib))
      end do
   end do
   !$omp end do
   deallocate(rt)
   !$omp end parallel

   ! Upward pass
   do ilev = nlev, 3, -1
      lv => self%lev(ilev)
      plv => self%lev(ilev - 1)
      do o = 0, 7
         call apply_translation(self%tm2m(:, :, o), &
            & lv%base + lv%oct_box(lv%oct_ptr(o):lv%oct_ptr(o + 1) - 1), &
            & plv%base + lv%oct_par(lv%oct_ptr(o):lv%oct_ptr(o + 1) - 1), mult)
      end do
   end do

   ! Interactions between well separated boxes of all levels
   call apply_m2l(self, mult, loc)

   ! Downward pass
   do ilev = 3, nlev
      lv => self%lev(ilev)
      plv => self%lev(ilev - 1)
      do o = 0, 7
         call apply_translation(self%tl2l(:, :, o), &
            & plv%base + lv%oct_par(lv%oct_ptr(o):lv%oct_ptr(o + 1) - 1), &
            & lv%base + lv%oct_box(lv%oct_ptr(o):lv%oct_ptr(o + 1) - 1), loc)
      end do
   end do

   ! Evaluate the local expansions of the leaf boxes
   lv => self%lev(nlev)
   !$omp parallel default(none) &
   !$omp shared(self, lv, loc, nlev, nleaf, p, phi, grad) &
   !$omp private(ib, ia, center, h, u, rt)
   allocate(rt(0:p, 0:p))
   !$omp do schedule(dynamic, 4)
   do ib = 1, nleaf
      call box_geometry(self, nlev, lv%key(ib), center, h)
      do ia = self%astart(ib), self%astart(ib + 1) - 1
         u = ([self%sx(ia), self%sy(ia), self%sz(ia)] - center) / h
         call regular_harmonics(u, p, rt)
         if (present(phi)) then
            phi(ia) = phi(ia) + eval_local(loc(:, lv%base + ib), rt, p) / h
         end if
         if (present(grad)) then
            grad(:, ia) = grad(:, ia) &
               & + eval_local_gradient(loc(:, lv%base + ib), rt, p) / (h * h)
         end if
      end do
   end do
   !$omp end do
   deallocate(rt)
   !$omp end parallel

   call near_field(self, qs, phi, grad)

end subroutine far_field


!> Add the multipole expansion of a point charge to the real coefficients
pure subroutine add_moment(q, rt, p, vec)
   !> Charge
   real(wp), intent(in) :: q
   !> Regular harmonics of the position
   complex(wp), intent(in) :: rt(0:, 0:)
   !> Order of the expansion
   integer, intent(in) :: p
   !> Real coefficients of the expansion
   real(wp), intent(inout) :: vec(:)

   integer :: l, m, i

   do l = 0, p
      vec(l * l + 1) = vec(l * l + 1) + q * real(rt(l, 0), wp)
      do m = 1, l
         i = l * l + 2 * m
         vec(i) = vec(i) + q * real(rt(l, m), wp)
         vec(i + 1) = vec(i + 1) - q * aimag(rt(l, m))
      end do
   end do

end subroutine add_moment


!> Potential of a local expansion, in the scaled units of the box
pure function eval_local(vec, rt, p) result(phi)
   !> Real coefficients of the local expansion
   real(wp), intent(in) :: vec(:)
   !> Regular harmonics of the position
   complex(wp), intent(in) :: rt(0:, 0:)
   !> Order of the expansion
   integer, intent(in) :: p
   !> Potential
   real(wp) :: phi

   integer :: l, m, i

   phi = 0.0_wp
   do l = 0, p
      phi = phi + vec(l * l + 1) * real(rt(l, 0), wp)
      do m = 1, l
         i = l * l + 2 * m
         phi = phi + 2.0_wp * (vec(i) * real(rt(l, m), wp) &
            & - vec(i + 1) * aimag(rt(l, m)))
      end do
   end do

end function eval_local


!> Gradient of a local expansion with respect to the scaled position
pure function eval_local_gradient(vec, rt, p) result(grad)
   !> Real coefficients of the local expansion
   real(wp), intent(in) :: vec(:)
   !> Regular harmonics of the position
   complex(wp), intent(in) :: rt(0:, 0:)
   !> Order of the expansion
   integer, intent(in) :: p
   !> Gradient
   real(wp) :: grad(3)

   integer :: j, k, l
   real(wp) :: fz, w
   complex(wp) :: fp, lk

   fz = 0.0_wp
   fp = (0.0_wp, 0.0_wp)
   do j = 1, p
      l = j - 1
      do k = 0, j
         if (k == 0) then
            lk = cmplx(vec(j * j + 1), 0.0_wp, wp)
            w = 1.0_wp
         else
            lk = cmplx(vec(j * j + 2 * k), vec(j * j + 2 * k + 1), wp)
            w = 2.0_wp
         end if
         if (k <= l) fz = fz + w * real(lk * rt(l, k), wp)
         if (k + 1 <= l) fp = fp + lk * rt(l, k + 1)
         if (k >= 1) fp = fp - conjg(lk * rt(l, k - 1))
      end do
   end do
   grad = [real(fp, wp), aimag(fp), fz]

end function eval_local_gradient


!> Inverse square root from a single precision estimate and one Newton step,
!> which vectorizes better than the double precision square root and division
!> and recovers the full double precision accuracy
elemental function inv_sqrt(r2) result(ri)
   !> Squared distance
   real(wp), intent(in) :: r2
   !> Inverse of the distance
   real(wp) :: ri

   real(wp) :: y0

   y0 = real(1.0_sp / sqrt(real(r2, sp)), wp)
   ri = y0 * (1.5_wp - 0.5_wp * r2 * y0 * y0)

end function inv_sqrt


!> Direct evaluation of the point charge interactions between neighbouring
!> leaf boxes
subroutine near_field(self, qs, phi, grad)
   !> Instance of the operator
   class(fmm_operator), intent(in) :: self
   !> Charges in the Morton order
   real(wp), intent(in) :: qs(:)
   !> Potential at the atoms in the Morton order
   real(wp), intent(inout), optional :: phi(:)
   !> Gradient of the potential at the atoms in the Morton order
   real(wp), intent(inout), optional :: grad(:, :)

   integer :: it, ib, jb, kb, ia, ja, a0, a1, b0, b1, ntask
   real(wp) :: xi, yi, zi, dx, dy, dz, r2, ri, acc, gx, gy, gz, w

   ntask = size(self%task_leaf)

   !$omp parallel do default(none) schedule(dynamic) &
   !$omp shared(self, qs, phi, grad, ntask) &
   !$omp private(it, ib, jb, kb, ia, ja, a0, a1, b0, b1, xi, yi, zi, dx, dy, dz) &
   !$omp private(r2, ri, acc, gx, gy, gz, w)
   do it = 1, ntask
      ib = self%task_leaf(it)
      a0 = self%task_a0(it)
      a1 = self%task_a1(it)
      do ia = a0, a1
         xi = self%sx(ia)
         yi = self%sy(ia)
         zi = self%sz(ia)
         acc = 0.0_wp
         gx = 0.0_wp
         gy = 0.0_wp
         gz = 0.0_wp
         do kb = self%near_ptr(ib), self%near_ptr(ib + 1) - 1
            jb = self%near_box(kb)
            b0 = self%astart(jb)
            b1 = self%astart(jb + 1) - 1
            if (jb == ib) then
               ! Skip the self interaction
               if (present(phi)) then
                  !$omp simd private(dx, dy, dz, r2, ri) reduction(+:acc)
                  do ja = b0, ia - 1
                     dx = xi - self%sx(ja)
                     dy = yi - self%sy(ja)
                     dz = zi - self%sz(ja)
                     r2 = dx * dx + dy * dy + dz * dz
                     ri = inv_sqrt(r2)
                     acc = acc + qs(ja) * ri
                  end do
                  !$omp simd private(dx, dy, dz, r2, ri) reduction(+:acc)
                  do ja = ia + 1, b1
                     dx = xi - self%sx(ja)
                     dy = yi - self%sy(ja)
                     dz = zi - self%sz(ja)
                     r2 = dx * dx + dy * dy + dz * dz
                     ri = inv_sqrt(r2)
                     acc = acc + qs(ja) * ri
                  end do
               end if
               if (present(grad)) then
                  !$omp simd private(dx, dy, dz, r2, ri, w) reduction(+:gx, gy, gz)
                  do ja = b0, ia - 1
                     dx = xi - self%sx(ja)
                     dy = yi - self%sy(ja)
                     dz = zi - self%sz(ja)
                     r2 = dx * dx + dy * dy + dz * dz
                     ri = inv_sqrt(r2)
                     w = qs(ja) * ri * ri * ri
                     gx = gx - w * dx
                     gy = gy - w * dy
                     gz = gz - w * dz
                  end do
                  !$omp simd private(dx, dy, dz, r2, ri, w) reduction(+:gx, gy, gz)
                  do ja = ia + 1, b1
                     dx = xi - self%sx(ja)
                     dy = yi - self%sy(ja)
                     dz = zi - self%sz(ja)
                     r2 = dx * dx + dy * dy + dz * dz
                     ri = inv_sqrt(r2)
                     w = qs(ja) * ri * ri * ri
                     gx = gx - w * dx
                     gy = gy - w * dy
                     gz = gz - w * dz
                  end do
               end if
            else
               if (present(phi)) then
                  !$omp simd private(dx, dy, dz, r2, ri) reduction(+:acc)
                  do ja = b0, b1
                     dx = xi - self%sx(ja)
                     dy = yi - self%sy(ja)
                     dz = zi - self%sz(ja)
                     r2 = dx * dx + dy * dy + dz * dz
                     ri = inv_sqrt(r2)
                     acc = acc + qs(ja) * ri
                  end do
               end if
               if (present(grad)) then
                  !$omp simd private(dx, dy, dz, r2, ri, w) reduction(+:gx, gy, gz)
                  do ja = b0, b1
                     dx = xi - self%sx(ja)
                     dy = yi - self%sy(ja)
                     dz = zi - self%sz(ja)
                     r2 = dx * dx + dy * dy + dz * dz
                     ri = inv_sqrt(r2)
                     w = qs(ja) * ri * ri * ri
                     gx = gx - w * dx
                     gy = gy - w * dy
                     gz = gz - w * dz
                  end do
               end if
            end if
         end do
         if (present(phi)) phi(ia) = phi(ia) + acc
         if (present(grad)) then
            grad(1, ia) = grad(1, ia) + gx
            grad(2, ia) = grad(2, ia) + gy
            grad(3, ia) = grad(3, ia) + gz
         end if
      end do
   end do

end subroutine near_field


!> Multiply the Coulomb matrix with a trial vector
subroutine multiply(self, xvec, yvec)
   !> Instance of the operator
   class(fmm_operator), intent(in) :: self
   !> Input vector
   real(wp), intent(in) :: xvec(:)
   !> Operator-vector product
   real(wp), intent(out) :: yvec(:)

   integer :: iat
   real(wp), allocatable :: qs(:), phi(:)

   ! Diagonal and correction for the Gaussian charge distribution
   call spsymv_csr(self%nat, self%val, self%list%inl, self%list%nlat, xvec, yvec)

   allocate(qs(self%nat), phi(self%nat))
   do iat = 1, self%nat
      qs(iat) = xvec(self%perm(iat))
      phi(iat) = 0.0_wp
   end do
   call far_field(self, qs, phi=phi)
   do iat = 1, self%nat
      yvec(self%perm(iat)) = yvec(self%perm(iat)) + phi(iat)
   end do

end subroutine multiply


!> Diagonal of the Coulomb matrix
subroutine get_diagonal(self, diag)
   !> Instance of the operator
   class(fmm_operator), intent(in) :: self
   !> Diagonal of the operator
   real(wp), intent(out) :: diag(:)

   diag(:) = self%diag

end subroutine get_diagonal


!> Dense block of the Coulomb matrix, exact for molecules
subroutine get_block(self, idx, block)
   !> Instance of the operator
   class(fmm_operator), intent(in) :: self
   !> Row and column indices of the block
   integer, intent(in) :: idx(:)
   !> Matrix elements of the block
   real(wp), intent(out) :: block(:, :)

   integer :: ii, jj, iat, jat
   real(wp) :: vec(3), r

   do jj = 1, size(idx)
      jat = idx(jj)
      block(jj, jj) = self%diag(jat)
      do ii = 1, jj - 1
         iat = idx(ii)
         vec = self%xyz(:, iat) - self%xyz(:, jat)
         r = norm2(vec)
         if (r > 1.0e-10_wp) then
            block(ii, jj) = erf(self%gam(self%id(iat), self%id(jat)) * r) / r
         else
            block(ii, jj) = 2.0_wp * self%gam(self%id(iat), self%id(jat)) / sqrtpi
         end if
         block(jj, ii) = block(ii, jj)
      end do
   end do

end subroutine get_block


!> Add scale * sum_ij p_i q_j dA_ij/dR to the gradient and the strain
!> derivative to sigma
subroutine add_gradient(self, pvec, qvec, gradient, sigma, scale, error)
   !> Instance of the operator
   class(fmm_operator), intent(in) :: self
   !> First charge-like vector
   real(wp), intent(in) :: pvec(:)
   !> Second charge-like vector
   real(wp), intent(in) :: qvec(:)
   !> Cartesian gradient
   real(wp), intent(inout) :: gradient(:, :)
   !> Strain derivative
   real(wp), intent(inout) :: sigma(:, :)
   !> Scaling factor of the contraction
   real(wp), intent(in) :: scale
   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   integer :: iat, ic
   real(wp), allocatable :: grad_loc(:, :), sigma_loc(:, :), qs(:), ps(:)
   real(wp), allocatable :: gq(:, :), gp(:, :)

   allocate(grad_loc(3, self%nat), source=0.0_wp)
   allocate(sigma_loc(3, 3), source=0.0_wp)

   ! Correction for the Gaussian charge distribution
   call add_pair_gradient(self%list, self%xyz, self%id, self%gam, alpha_big, &
      & pvec, qvec, scale, grad_loc, sigma_loc)

   ! The field of both charge vectors at the atoms
   allocate(qs(self%nat), ps(self%nat), source=0.0_wp)
   allocate(gq(3, self%nat), gp(3, self%nat), source=0.0_wp)
   do iat = 1, self%nat
      qs(iat) = qvec(self%perm(iat))
      ps(iat) = pvec(self%perm(iat))
   end do
   call far_field(self, qs, grad=gq)
   call far_field(self, ps, grad=gp)
   do iat = 1, self%nat
      grad_loc(:, self%perm(iat)) = grad_loc(:, self%perm(iat)) + scale &
         & * (pvec(self%perm(iat)) * gq(:, iat) + qvec(self%perm(iat)) * gp(:, iat))
   end do

   ! The strain derivative of a molecule follows from the virial of the forces
   sigma_loc(:, :) = 0.0_wp
   do iat = 1, self%nat
      do ic = 1, 3
         sigma_loc(:, ic) = sigma_loc(:, ic) + self%xyz(:, iat) * grad_loc(ic, iat)
      end do
   end do

   gradient(:, :) = gradient + grad_loc
   sigma(:, :) = sigma + sigma_loc

end subroutine add_gradient


!> Write a summary of the operator setup
subroutine write_info(self, unit)
   !> Instance of the operator
   class(fmm_operator), intent(in) :: self
   !> Output unit
   integer, intent(in) :: unit

   integer :: nleaf
   real(wp) :: mem

   write(unit, "(a, t32, a)") "Coulomb operator", "fast multipole method"
   write(unit, "(a, t32, i12)") "Expansion order", self%order
   write(unit, "(a, t32, f12.6)") "Correction cutoff / bohr", self%cutoff
   write(unit, "(a, t32, i12)") "Finest tree level", self%nlevel
   nleaf = size(self%near_ptr) - 1
   write(unit, "(a, t32, i12)") "Leaf boxes", nleaf
   write(unit, "(a, t32, f12.2)") "Atoms per leaf box", real(self%nat, wp) / nleaf
   write(unit, "(a, t32, f12.6)") "Leaf box edge / bohr", &
      & self%edge / real(2**self%nlevel, wp)
   write(unit, "(a, t32, i12)") "Near field pairs", self%npair_near
   write(unit, "(a, t32, i12)") "Multipole translations", self%npair_m2l
   write(unit, "(a, t32, i12)") "Correction pairs", size(self%list%nlat) - self%nat
   mem = 0.0_wp
   if (allocated(self%tm2l)) mem = mem + 8.0_wp * size(self%tm2l)
   if (allocated(self%tm2m)) mem = mem + 8.0_wp * (size(self%tm2m) + size(self%tl2l))
   write(unit, "(a, t32, f12.2)") "Translation matrices / MB", mem * 1.0e-6_wp

end subroutine write_info

end module multicharge_coulomb_fmm
