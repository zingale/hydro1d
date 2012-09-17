module user_bc_module

  use grid_module
  use variables_module
  use probparams_module
  use datatypes_module

  implicit none

  real (kind=dp_t) :: e_fluff, rho_fluff

contains

  subroutine user_bc_xp(U)

    type(gridvar_t), intent(inout) :: U

    integer :: i
    real (kind=dp_t) :: e, v

    do i = U%grid%hi+1, U%grid%hi+U%grid%ng

       ! give all quantities a zero-gradient
       U%data(i,:) = U%data(i-1,:)

       ! velocity -- outgoing only
       !v = max(0.0_dp_t, U%data(i,iumomx)/U%data(i,iudens))
       v = U%data(i,iumomx)/U%data(i,iudens)

       ! set the density to the small value
       U%data(i,iudens) = rho_fluff

       ! recompute energy
       U%data(i,iuener) = rho_fluff*e_fluff + 0.5_dp_t*rho_fluff*v**2
       U%data(i,iumomx) = rho_fluff*v

    enddo

  end subroutine user_bc_xp
end module user_bc_module
