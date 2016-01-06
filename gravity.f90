subroutine gravity(U, g)

  use datatypes_module
  use grid_module
  use variables_module

  implicit none

  type(gridedgevar_t), intent(in) :: U
  type(gridedgevar_t), intent(out) :: g

  ! compute the gravitational acceleration on edges

end subroutine gravity
