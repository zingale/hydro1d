module update_module

  use datatypes_module
  use grid_module
  use variables_module
  use params_module

  implicit none

  private

  public :: update

contains
  
  subroutine update(U, fluxes, dt)

    type(gridvar_t),     intent(inout) :: U
    type(gridedgevar_t), intent(in   ) :: fluxes
    real (kind=dp_t),    intent(in   ) :: dt

    type(gridvar_t) :: Uold

    integer :: i

    ! store the old state
    call build(Uold, U%grid, U%nvar)
    do i = U%grid%lo, U%grid%hi
       Uold%data(i,:) = U%data(i,:)
    enddo
    
    ! update -- this should be (A_l F_l - A_r F_r)/dV + any gradient terms
    do i = U%grid%lo, U%grid%hi
       U%data(i,:) = U%data(i,:) + &
            (dt/U%grid%dx)*(fluxes%data(i,:) - fluxes%data(i+1,:))
    enddo


    ! time-centered source terms 
    do i = U%grid%lo, U%grid%hi
       U%data(i,iumomx) = U%data(i,iumomx) + &
            0.5_dp_t*dt*(U%data(i,iudens) + Uold%data(i,iudens))*grav

       U%data(i,iuener) = U%data(i,iuener) + &
            0.5_dp_t*dt*(U%data(i,iumomx) + Uold%data(i,iumomx))*grav
    enddo


    ! clean-up
    call destroy(Uold)


  end subroutine update

end module update_module
