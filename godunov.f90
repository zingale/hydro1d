module interface_states_godunov_module

  use datatypes_module
  use grid_module
  use variables_module

  implicit none

contains
  
  subroutine make_interface_states_godunov(U, U_l, U_r, dt)

    type(gridvar_t),     intent(in   ) :: U
    type(gridedgevar_t), intent(inout) :: U_l, U_r
    real (kind=dp_t),    intent(in   ) :: dt

    integer :: i

    ! piecewise constant slopes
    
    ! loop over interfaces and fill the states 
    do i = U%grid%lo, U%grid%hi+1
       U_l%data(i,:) = U%data(i-1,:)
       U_r%data(i,:) = U%data(i,:)
    enddo

  end subroutine make_interface_states_godunov

end module interface_states_godunov_module
