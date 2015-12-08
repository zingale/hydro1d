module init_module

  use datatypes_module
  use probparams_module
  use params_module
  use grid_module
  use variables_module

  implicit none

  private

  public :: init_data

contains
  
  subroutine init_data(U)

    type(gridvar_t), intent(inout) :: U

    integer :: i
    real (kind=dp_t) :: rhoe

    real (kind=dp_t), parameter :: E_sedov = ONE
    real (kind=dp_t) :: p, p_exp, p_amb, V_exp

    ! we'll put all of the energy in the first zone
    V_exp = 4.0*pi/3.0 * U%grid%dx**3
    p_amb = 1.e-5_dp_t

    do i = U%grid%lo, U%grid%hi
      
       ! store the conserved state
       U%data(i,iudens) = ONE
       U%data(i,iumomx) = ZERO
       
       if (i == U%grid%lo) then
          rhoe = E_sedov/V_exp
       else
          rhoe = p_amb/(gamma - ONE)
       endif

       U%data(i,iuener) = rhoe
          
    enddo

  end subroutine init_data

end module init_module

