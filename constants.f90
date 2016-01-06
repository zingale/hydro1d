module constants_module

  use datatypes_module

  implicit none

  ! physical constants
  real (kind=dp_t), parameter :: m_amu = 1.66053886e-24_dp_t   ! g
  real (kind=dp_t), parameter :: k = 1.3806503e-16_dp_t        ! erg / K
  real (kind=dp_t), parameter :: G_newton = 6.67259e-8_dp_t    ! cm^3 / g / s^2

end module constants_module
