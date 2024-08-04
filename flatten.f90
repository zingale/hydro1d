module flatten_module

  use datatypes_module
  use grid_module
  use variables_module

  implicit none

contains

  subroutine flatten(Q, xi)

    ! flattening kicks in behind strong shocks and reduces the
    ! reconstruction to using piecewise constant slopes, making things
    ! first-order.  See Saltzman (1994) page 159 for this
    ! implementation.

    type(gridvar_t), intent(in) :: Q
    type(gridvar_t), intent(inout) :: xi

    real (kind=dp_t) :: dp, dp2, z
    real (kind=dp_t), parameter :: z0 = 0.75_dp_t
    real (kind=dp_t), parameter :: z1 = 0.85_dp_t
    real (kind=dp_t), parameter :: delta = 0.33_dp_t
    real (kind=dp_t), parameter :: smallp = 1.e-10_dp_t

    type(gridvar_t) :: xi_t

    integer :: i

    call build(xi_t, Q%grid, 1)

    do i = Q%grid%lo-2, Q%grid%hi+2

       dp  = abs(Q%data(i+1,iqpres) - Q%data(i-1,iqpres))
       dp2 = abs(Q%data(i+2,iqpres) - Q%data(i-2,iqpres))
       z = dp/max(dp2,smallp)

       if ( (Q%data(i-1,iqxvel) - Q%data(i+1,iqxvel) > 0.0_dp_t) .and. &
            (dp/min(Q%data(i+1,iqpres),Q%data(i-1,iqpres)) > delta) ) then
          xi_t%data(i,1) = min(1.0_dp_t, &
                             max(0.0_dp_t, 1.0_dp_t - (z - z0)/(z1 - z0)))
       else
          xi_t%data(i,1) = 1.0_dp_t
       endif

    enddo

    do i = Q%grid%lo-1, Q%grid%hi+1

       if (Q%data(i+1,iqpres) - Q%data(i-1,iqpres) > 0.0_dp_t) then
          xi%data(i,1) = min(xi_t%data(i,1), xi_t%data(i-1,1))

       else if (Q%data(i+1,iqpres) - Q%data(i-1,iqpres) < 0.0_dp_t) then
          xi%data(i,1) = min(xi_t%data(i,1), xi_t%data(i+1,1))

       else
          xi%data(i,1) = xi_t%data(i,1)
       endif

    enddo

    call destroy(xi_t)
  end subroutine flatten

end module flatten_module
