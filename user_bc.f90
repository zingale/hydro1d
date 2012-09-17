module user_bc_module

  use grid_module
  use variables_module
  use datatypes_module

  implicit none

  private

  public :: user_bc_xm, user_bc_xp

contains

  subroutine user_bc_xm(U)

    type(gridvar_t), intent(inout) :: U

  end subroutine user_bc_xp

  subroutine user_bc_xp(U)

    type(gridvar_t), intent(inout) :: U

  end subroutine user_bc_xp

end module user_bc_module
