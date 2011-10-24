module bcs_module

  use datatypes_module
  use grid_module
  use variables_module
  use params_module
  implicit none

contains

  subroutine fillBCs(U)

    type(gridvar_t), intent(inout) :: U

    integer :: i, ir, ip


    ! sanity check
    if ( (xlboundary == "periodic" .and. xrboundary /= "periodic") .or. &
         (xrboundary == "periodic" .and. xlboundary /= "periodic") ) then
       print *, "ERROR: both boundaries must be periodic"
       stop
    endif


    ! lower boundary (-x)
    select case (xlboundary)

    case ("reflect")

       ! ir is the index that corresponds to the reflected ghostcell
       ir = U%grid%lo
       do i = U%grid%lo-1, -U%grid%ng
          
          ! give all quantities a zero-gradient
          U%data(i,:) = U%data(i+1,:)

          ! reflect the velocities.  
          U%data(i,iumomx) = -U%data(ir,iumomx)

          ir = ir + 1

       enddo

    case ("outflow")

       do i = U%grid%lo-1, -U%grid%ng
          
          ! give all quantities a zero-gradient
          U%data(i,:) = U%data(i+1,:)

       enddo

    case ("periodic")

       ip = U%grid%hi
       do i = U%grid%lo-1, -U%grid%ng
          
          U%data(i,:) = U%data(ip,:)

          ip = ip - 1

       enddo

    case default

       print *, "ERROR: xlboundary not valid"
       stop

    end select


    ! upper boundary (+x)
    select case (xrboundary)

    case ("reflect")

       ! ir is the index that corresponds to the reflected ghostcell
       ir = U%grid%hi
       do i = U%grid%hi+1, U%grid%hi+U%grid%ng
          
          ! give all quantities a zero-gradient
          U%data(i,:) = U%data(i-1,:)

          ! reflect the velocities.  
          U%data(i,iumomx) = -U%data(ir,iumomx)

          ir = ir - 1

       enddo

    case ("outflow")

       do i = U%grid%hi+1, U%grid%hi+U%grid%ng
          
          ! give all quantities a zero-gradient
          U%data(i,:) = U%data(i-1,:)

       enddo

    case ("periodic")

       ip = U%grid%lo
       do i = U%grid%hi+1, U%grid%hi+U%grid%ng
          
          U%data(i,:) = U%data(ip,:)

          ip = ip + 1

       enddo


    case default

       print *, "ERROR: xlboundary not valid"
       stop

    end select

  end subroutine fillBCs

end module bcs_module
