module riemann_module

  use datatypes_module
  use grid_module
  use variables_module

  implicit none

contains
  
  subroutine solve_riemann(U_l, U_r, fluxes)

    type(gridedgevar_t), intent(in   ) :: U_l, U_r
    type(gridedgevar_t), intent(in   ) :: fluxes

  end subroutine solve_riemann

end module riemann_module
