! dt module
!
! compute the timestep as the minimum dx / (|U| + cs)

module dt_module

  use datatypes_module
  use grid_module
  use params_module
  use eos_module
  use variables_module

  implicit none

contains

  subroutine compute_dt(U, n, dt)

    type(gridvar_t),  intent(in   ) :: U
    integer,          intent(inout) :: n
    real (kind=dp_t), intent(inout) :: dt

    integer :: i
    real (kind=dp_t) :: cs, p, e

    
    dt = huge(0.0_dp_t)

    do i = U%grid%lo, U%grid%hi

       ! compute cs (soundspeed)
       e = (U%data(i,iuener) - 0.5_dp_t*U%data(i,iumomx)**2/U%data(i,iudens)) / &
            U%data(i,iudens)
       
       call eos(eos_input_e, p, e, U%data(i,iudens))

       cs = sqrt(gamma*p/U%data(i,iudens))

       dt = min(dt, U%grid%dx/(abs(U%data(i,iumomx)/U%data(i,iudens)) + cs))

    enddo    

    dt = cfl*dt

    if (n == 0) dt = init_shrink*dt

  end subroutine compute_dt

end module dt_module
