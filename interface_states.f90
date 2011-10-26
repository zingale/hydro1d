module interface_states_module

  use datatypes_module
  use grid_module
  use variables_module
  use params_module

  implicit none

contains
  
  subroutine make_interface_states(U, U_l, U_r, dt)

    type(gridvar_t),     intent(in   ) :: U
    type(gridedgevar_t), intent(  out) :: U_l, U_r
    real (kind=dp_t),    intent(in   ) :: dt

    integer :: i

    select case (godunov_type)

    case (0)
       ! piecewise constant slopes
       
       ! loop over interfaces and fill the states 
       do i = U%grid%lo, U%grid%hi+1
          U_l%data(i,:) = U%data(i-1,:)
          U_r%data(i,:) = U%data(i,:)
       enddo

    case (1)
       ! piecewise linear slopes

    case (2)
       ! piecewise parabolic slopes
       
    case default
       print *, "ERROR: invalid godunov_type"
       stop

    end select


  end subroutine make_interface_states

end module interface_states_module
