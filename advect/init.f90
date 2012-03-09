module init_module

  use datatypes_module
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
    real (kind=dp_t) :: rho, ux, p, e

    real (kind=dp_t) :: xcenter, s

    ! here we setup a constant pressure, constant (non-zero) velocity,
    ! and a varying density field that will simply be advected.  The
    ! particular choices here come from the FLASH Code paper (Fryxell
    ! et al. 2000).

    xcenter = 0.5_dp_t*(U%grid%xmin + U%grid%xmax)


    do i = U%grid%lo, U%grid%hi

       ! compute the distance from the center
       s = U%grid%x(i) - xcenter

       if (perttype == "gaussian") then
          rho = rho1*phi_g(s/width) + rho0*(1.0_dp_t - phi_g(s/width))
       else if (perttype == "tophat") then
          rho = rho1*phi_t(s/width) + rho0*(1.0_dp_t - phi_t(s/width))
       else
          print *, "ERROR: invalid perttype"
          stop
       endif

       ux  = u0
       p   = p0

       ! get the internal energy from the EOS
       call eos(eos_input_p, p, e, rho)
       
       ! store the conserved state
       U%data(i,iudens) = rho
       U%data(i,iumomx) = rho*ux
       U%data(i,iuener) = 0.5_dp_t*rho*ux*ux + rho*e
          
    enddo

  end subroutine init_data

  function phi_g(q)

    real (kind=dp_t) :: phi_g, q

    phi_g = exp(-q**2)

    return
  end function phi_g

  function phi_t(q)

    real (kind=dp_t) :: phi_t, q

    phi_t = 0.0_dp_t
    if (abs(q) < 1.0_dp_t) then
       phi_t = 1.0_dp_t
    endif

    return
  end function phi_t

end module init_module

