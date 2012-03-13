! a simple gamma-law equation of state
! p = rho e (gamma - 1)
module eos_module

  use datatypes_module
  use params_module, only: gamma
  implicit none

  integer, parameter :: eos_input_p = 1
  integer, parameter :: eos_input_e = 2
  integer, parameter :: eos_input_pe = 3

  private

  public :: eos_input_p, eos_input_e, eos_input_pe
  public :: eos

contains

  subroutine eos(input, p, e, rho)
    
    integer,          intent(in   ) :: input
    real (kind=dp_t), intent(inout) :: p
    real (kind=dp_t), intent(inout) :: e
    real (kind=dp_t), intent(inout) :: rho

    if (input == eos_input_p) then

       if (p < 0.0_dp_t) then
          print *, 'ERROR: input pressure < 0 in EOS'
          stop
       endif

       e = p / (rho * (gamma - 1.0_dp_t))

    else if (input == eos_input_e) then

       if (e < 0.0_dp_t) then
          print *, 'ERROR: input internal energy < 0 in EOS'
          stop
       endif

       p = rho*e*(gamma - 1.0_dp_t)

    else if (input == eos_input_pe) then

       ! pressure and energy are inputs -- return density
       if (e < 0.0_dp_t .or. p < 0.0_dp_t) then
          print *, 'ERROR: input internal energy or pressure < 0 in EOS'
          stop
       endif

       rho = p/e/(gamma - 1.0_dp_t)

    else
       print *, 'ERROR: invalid input in EOS'
       stop
    endif

  end subroutine eos

end module eos_module
