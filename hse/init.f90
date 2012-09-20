module init_module

  use datatypes_module
  use params_module
  use probparams_module
  use grid_module
  use eos_module
  use variables_module

  implicit none

  private

  public :: init_data

contains
  
  subroutine init_data(U)

    type(gridvar_t), intent(inout) :: U

    integer :: i
    real (kind=dp_t) :: pres, pres_below, e, H


    ! compute the pressure scale height (isothermal gas)
    H = pres_base / dens_base / abs(grav)

    
    ! first set the density profile
    do i = U%grid%lo, U%grid%hi

       if (do_isentropic) then
          ! we can integrate HSE with p = K rho^gamma analytically
          U%data(i,iudens) = dens_base*(grav*dens_base* &
               (gamma - 1.0_dp_t)*U%grid%x(i)/ &
               (gamma*pres_base) + 1.0_dp_t)**(1.0_dp_t/(gamma - 1.0_dp_t))
          
       else
          ! the density of an isothermal gamma-law atm is exponential
          U%data(i,iudens) = dens_base * exp(-U%grid%x(i)/H)

       endif

       if (U%data(i,iudens) < small_dens) then
          U%data(i:,iudens) = small_dens
          exit
       endif
    enddo

    ! now determine the pressure via HSE and set the state by using the EOS
    pres_below = pres_base
    do i = U%grid%lo, U%grid%hi
       if (U%data(i,iudens) > small_dens) then
          pres = pres_below + &
               0.5_dp_t*U%grid%dx*(U%data(i,iudens) + U%data(i-1,iudens)) * grav
       else
          pres = pres_below
       endif

       ! get the internal energy from the EOS
       call eos(eos_input_p, pres, e, U%data(i,iudens))

       ! store the remainder of the conserved state
       U%data(i,iumomx) = 0.0_dp_t
       U%data(i,iuener) = U%data(i,iudens)*e    ! no kinetic energy

       pres_below = pres

       print *, i, U%data(i,iudens), pres

    enddo

  end subroutine init_data

end module init_module

