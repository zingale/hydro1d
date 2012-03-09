module update_module

  use datatypes_module
  use grid_module
  use variables_module

  implicit none

contains
  
  subroutine update(U, fluxes, dt)

    type(gridvar_t),     intent(inout) :: U
    type(gridedgevar_t), intent(in   ) :: fluxes
    real (kind=dp_t),    intent(in   ) :: dt

    integer :: i

    do i = U%grid%lo, U%grid%hi
       U%data(i,:) = U%data(i,:) + &
            (dt/U%grid%dx)*(fluxes%data(i,:) - fluxes%data(i+1,:))
    enddo

  end subroutine update

end module update_module
