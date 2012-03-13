module user_bc_module

  use grid_module
  use variables_module
  use probparams_module
  use datatypes_module

  implicit none

contains

  subroutine user_bc_xp(U)

    type(gridvar_t), intent(inout) :: U

    integer :: i
    real (kind=dp_t) :: e

    do i = U%grid%hi+1, U%grid%hi+U%grid%ng

       ! give all quantities a zero-gradient
       U%data(i,:) = U%data(i-1,:)

       ! set the density to the small value
       U%data(i,iudens) = small_dens

       ! store internal energy
       e = (U%data(i,iuener) - &
            0.5_dp_t*U%data(i,iumomx)**2/U%data(i,iudens))/U%data(i,iudens)
       
       U%data(i,iumomx) = max(0.0_dp_t, U%data(i,iumomx))

       ! recompute energy
       U%data(i,iuener) = U%data(i,iudens)*e + &
            0.5_dp_t*U%data(i,iumomx)**2/U%data(i,iudens)
       
    enddo

  end subroutine user_bc_xp
end module user_bc_module
