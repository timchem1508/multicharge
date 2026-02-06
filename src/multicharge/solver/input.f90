module mchrg_solver_input
   use iterative_solver, only : cg_input

   implicit none
   private


   !> Collection of possible solvation models
   type, public :: solver_input
      !> Input for Conjugate gradient iterative solver
      type(cg_input), allocatable :: cg
   end type solver_input

end module mchrg_solver_input