! a simple gamma-law equation of state
! p = rho e (gamma - 1)
module eos_module

  use datatypes_module
  use params_module, only: gamma
  implicit none

  integer, parameter :: eos_input_p = 1
  integer, parameter :: eos_input_e = 2

contains

  subroutine eos(input, p, e, rho)
    
    integer,          intent(in   ) :: input
    real (kind=dp_t), intent(inout) :: p
    real (kind=dp_t), intent(inout) :: e
    real (kind=dp_t), intent(in   ) :: rho

    if (input == eos_input_p) then

       e = p / (rho * (gamma - 1.0_dp_t))

    else if (input == eos_input_e) then

       p = rho*e*(gamma - 1.0_dp_t)

    endif

  end subroutine eos

end module eos_module
