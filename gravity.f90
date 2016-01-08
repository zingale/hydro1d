subroutine gravity(U, g)

  ! compute the cell-centered gravitational acceleration

  use datatypes_module
  use grid_module
  use variables_module
  use params_module, only : grav, gravity_monopole
  use constants_module

  implicit none

  type(gridvar_t), intent(in) :: U
  type(gridvar_t), intent(inout) :: g

  type(grid_t) :: grid

  integer :: i, ic, ip, ir
  real (kind=dp_t) :: M_enclosed, g_edge_lo, g_edge_hi
  
  grid = U%grid

  if (gravity_monopole == 0) then
     g%data(:,1) = grav
  else

     ! compute the acceleration on edges and average to the center.
     ! This is the same as doing the potential at centers and doing
     ! a centered difference
     M_enclosed = ZERO

     g_edge_lo = ZERO

     do i = grid%lo, grid%hi
        ! mass enclosed in this zone
        M_enclosed = M_enclosed + FOURTHIRDS*pi* &
             (grid%xr(i)**2 + grid%xl(i)**2 + grid%xl(i)*grid%xr(i)) * &
             grid%dx * &
             U%data(i,iudens)
        
        ! gravitational acceleration on this right edge
        g_edge_hi = -G_newton*M_enclosed/grid%xr(i)**2

        ! cell-centered gravity
        g%data(i,1) = HALF*(g_edge_lo + g_edge_hi)

        g_edge_lo = g_edge_hi
     enddo
  endif

  ! fill boundary conditions -- this is necessary for the interface
  ! state construction

  ! Note: at the moment, we cannot use the bc routines here, since they 
  ! are coded specifically for hydro, so we'll hack something in.  This
  ! needs to be improved.

  ! lower boundary (-x)
  select case (U%grid%xlboundary)
       
  case ("reflect")
       
     ! ir is the index that corresponds to the reflected ghostcell
     ir = U%grid%lo
     do i = U%grid%lo-1, U%grid%lo-U%grid%ng, -1
          
        ! odd-reflect the acceleration
        g%data(i,1) = -g%data(ir,1)
        ir = ir + 1
     enddo
        
    case ("outflow", "hse", "user")
       do i = U%grid%lo-1, U%grid%lo-U%grid%ng, -1
          ! give all quantities a zero-gradient
          g%data(i,1) = g%data(i+1,1)
       enddo

    case ("periodic")
       ip = U%grid%hi
       do i = U%grid%lo-1, U%grid%lo-U%grid%ng, -1
          g%data(i,1) = g%data(ip,1)
          ip = ip - 1
       enddo

    case default
       print *, "ERROR: xlboundary not valid"
       stop
       
    end select
     

    ! upper boundary (+x)
    select case (U%grid%xrboundary)       
    case ("reflect")
       
       ! ir is the index that corresponds to the reflected ghostcell
       ir = U%grid%hi
       do i = U%grid%hi+1, U%grid%hi+U%grid%ng          
          ! odd reflect the acceleration
          g%data(i,1) = -g%data(ir,1)
          ir = ir - 1
       enddo

    case ("outflow", "hse", "user", "diode")
       do i = U%grid%hi+1, U%grid%hi+U%grid%ng
          ! give all quantities a zero-gradient
          g%data(i,1) = g%data(i-1,1)
       enddo
       
    case ("periodic")
       
       ip = U%grid%lo
       do i = U%grid%hi+1, U%grid%hi+U%grid%ng
          g%data(i,1) = g%data(ip,1)
          ip = ip + 1
       enddo
       
    case default
       print *, "ERROR: xrboundary not valid"
       stop
    end select

end subroutine gravity
