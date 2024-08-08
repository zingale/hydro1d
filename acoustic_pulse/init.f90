module init_module

  use datatypes_module
  use probparams_module
  use grid_module
  use eos_module
  use variables_module
  use params_module, only: gamma

  implicit none

  private

  public :: init_data

contains

  subroutine init_data(U)

    type(gridvar_t), intent(inout) :: U

    integer :: i
    real (kind=dp_t) :: rho, ux, p, e

    real (kind=dp_t) :: xcenter, r

    ! here we setup a constant pressure, constant (non-zero) velocity,
    ! and a varying density field that will simply be advected.  The
    ! particular choices here come from the FLASH Code paper (Fryxell
    ! et al. 2000).

    xcenter = 0.5_dp_t*(U%grid%xmin + U%grid%xmax)


    do i = U%grid%lo, U%grid%hi

       ! compute the distance from the center
       r = abs(U%grid%x(i) - xcenter)

       if (r <= 0.5_dp_t) then
          rho = rho0 + drho * exp(-16.0_dp_t * r) * cos(pi * r)**6
       else
          rho = rho0
       endif

       ux  = 0
       p   = (rho / rho0)**gamma

       ! get the internal energy from the EOS
       call eos(eos_input_p, p, e, rho)

       ! store the conserved state
       U%data(i,iudens) = rho
       U%data(i,iumomx) = rho*ux
       U%data(i,iuener) = 0.5_dp_t*rho*ux*ux + rho*e

    enddo

  end subroutine init_data

end module init_module
