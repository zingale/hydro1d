! runtime diagonstics

module runtime_diag_module

  use grid_module
  use datatypes_module
  use params_module
  use variables_module
  use eos_module

  implicit none

  private

  public :: run_diag

contains

  subroutine run_diag(U, t)

    type(gridvar_t),  intent(in) :: U
    real (kind=dp_t), intent(in) :: t



  end subroutine run_diag
end module runtime_diag_module

