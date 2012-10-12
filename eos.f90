! a simple gamma-law equation of state: 
!    p = rho e (gamma - 1)
!
! to define the temperature and entropy, we will use the ideal gas law:
!    p = rho k T / (mu m_amu)
!    
! where mu is the mean molecular weight and m_amu is the atomic mass
! unit

module eos_module

  use datatypes_module
  use params_module, only: gamma
  implicit none

  integer, parameter :: eos_input_p  = 1
  integer, parameter :: eos_input_e  = 2
  integer, parameter :: eos_input_pe = 3
  integer, parameter :: eos_input_ps = 4
  private

  public :: eos_input_p, eos_input_e, eos_input_pe, eos_input_ps
  public :: eos

contains

  subroutine eos(input, p, e, rho, T, s)
    
    integer,          intent(in   ) :: input
    real (kind=dp_t), intent(inout) :: p
    real (kind=dp_t), intent(inout) :: e
    real (kind=dp_t), intent(inout) :: rho
    real (kind=dp_t), intent(inout), optional :: s
    real (kind=dp_t), intent(inout), optional :: T

    ! physical constants
    real (kind=dp_t), parameter :: mu = 1.0_dp_t
    real (kind=dp_t), parameter :: m_amu = 1.66053886e-24_dp_t   ! g
    real (kind=dp_t), parameter :: k = 1.3806503e-16_dp_t        ! erg / K

    real (kind=dp_t) :: c_v 

    ! c_v = de/dT |_rho, so since e = p / (rho (gamma - 1)), and
    ! using p = rho k T / (mu m_amu), we have e = k/(mu m (gamma - 1))
    c_v = k/((gamma - 1.0_dp_t)*mu*m_amu)

    if (input == eos_input_p) then

       if (p < 0.0_dp_t) then
          print *, 'ERROR: input pressure < 0 in EOS'
          stop
       endif

       e = p / (rho * (gamma - 1.0_dp_t))

       if (present(T)) T = mu*m_amu*e*(gamma - 1.0_dp_t)/k
       if (present(s)) s = c_v*log(p/rho**gamma) 

    else if (input == eos_input_e) then

       if (e < 0.0_dp_t) then
          print *, 'ERROR: input internal energy < 0 in EOS'
          stop
       endif

       p = rho*e*(gamma - 1.0_dp_t)

       if (present(T)) T = mu*m_amu*e*(gamma - 1.0_dp_t)/k
       if (present(s)) s = c_v*log(p/rho**gamma) 

    else if (input == eos_input_pe) then

       ! pressure and energy are inputs -- return density
       if (e < 0.0_dp_t .or. p < 0.0_dp_t) then
          print *, 'ERROR: input internal energy or pressure < 0 in EOS'
          stop
       endif

       rho = p/e/(gamma - 1.0_dp_t)

       if (present(T)) T = mu*m_amu*e*(gamma - 1.0_dp_t)/k
       if (present(s)) s = c_v*log(p/rho**gamma) 

    else if (input == eos_input_ps) then

       ! entropy and pressure and inputs, return density and e 
       if (.not. present(s)) then
          print *, 'ERROR: entropy must be input'
          stop
       endif
       
       ! note: our definition of entropy neglects the constant since
       ! we really only care about changes in entropy.
       rho = (p*exp(-s/c_v))**(1.0_dp_t/gamma)

       e = p / (rho * (gamma - 1.0_dp_t))

       if (present(T)) T = mu*m_amu*e*(gamma - 1.0_dp_t)/k

    else
       print *, 'ERROR: invalid input in EOS'
       stop
    endif


  end subroutine eos

end module eos_module
