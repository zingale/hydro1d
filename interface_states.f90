module interface_states_module

  use datatypes_module
  use grid_module
  use variables_module

  implicit none

contains
  
  subroutine make_interface_states(U, U_l, U_r, dt)

    type(gridvar_t),     intent(in   ) :: U
    type(gridedgevar_t), intent(  out) :: U_l, U_r
    real (kind=dp_t),    intent(in   ) :: dt


  end subroutine make_interface_states

end module interface_states_module
