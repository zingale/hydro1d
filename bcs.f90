module bcs_module

  use datatypes_module
  use grid_module
  use variables_module
  use params_module
  use eos_module
  use user_bc_module

  implicit none

  private

  public :: fillBCs
   
contains

  subroutine fillBCs(U)

    type(gridvar_t), intent(inout) :: U
    
    integer :: i, ir, ip, ic
    real (kind=dp_t) :: p, e, p_above, e_above, p_below, e_below
    real (kind=dp_t) :: econst, dens, vel
    
    ! sanity check
    if ( (U%grid%xlboundary == "periodic" .and. &
         U%grid%xrboundary /= "periodic") .or. &
         (U%grid%xrboundary == "periodic" .and. &
         U%grid%xlboundary /= "periodic") ) then
       print *, "ERROR: both boundaries must be periodic"
       stop
    endif
    

    !-------------------------------------------------------------------------
    ! lower boundary (-x)
    !-------------------------------------------------------------------------
    select case (U%grid%xlboundary)
       
    case ("reflect")
       
       ! ir is the index that corresponds to the reflected ghostcell
       ir = U%grid%lo
       do i = U%grid%lo-1, U%grid%lo-U%grid%ng, -1
          
          ! even-reflect density and energy
          U%data(i,iudens) = U%data(ir,iudens)
          U%data(i,iuener) = U%data(ir,iuener)
          
          ! odd-reflect the velocities.  
          U%data(i,iumomx) = -U%data(ir,iumomx)
          
          ir = ir + 1
          
       enddo
        

    case ("hse")
        
       ir = U%grid%lo
       ic = U%grid%lo
             
       do i = U%grid%lo-1, U%grid%lo-U%grid%ng, -1
          
          if (hse_bc_const == "density") then
              
             ! constant density in the ghost cells
              
             ! zero gradient to rho
             U%data(i,iudens) = U%data(ic,iudens)

             ! velocity depends on hse_vel_type
             select case (hse_vel_type)
             case ("outflow")
                vel = U%data(ic,iumomx)/U%data(ic,iudens)

             case ("reflect")
                vel = -U%data(ir,iumomx)/U%data(ir,iudens)
                ir = ir + 1

             case default
                print *, "ERROR: invalid hse_vel_type"
                stop
                
             end select

             U%data(i,iumomx) = U%data(i,iudens)*vel

             ! p via HSE
             e_above = (U%data(i+1,iuener) - &
                  0.5_dp_t*U%data(i+1,iumomx)**2/U%data(i+1,iudens))/ &
                  U%data(i+1,iudens)
             call eos(eos_input_e, p_above, e_above, U%data(i+1,iudens))
             
             p = p_above - 0.5_dp_t*(U%data(i,iudens) + &
                                     U%data(i+1,iudens))*grav*U%grid%dx
             
             call eos(eos_input_p, p, e, U%data(i,iudens))
             
             ! compute E
             U%data(i,iuener) = U%data(i,iudens)*e + &
                  0.5_dp_t*U%data(i,iumomx)**2/U%data(i,iudens)
             
          else if (hse_bc_const == "temperature") then
             
             ! constant temperature (or internal energy) in the ghost
             ! cells
             
             econst = (U%data(ic,iuener) - &
                  0.5_dp_t*U%data(ic,iumomx)**2/U%data(ic,iudens))/ &
                  U%data(ic,iudens)
             
             ! velocity depends on hse_vel_type
             select case (hse_vel_type)
             case ("outflow")
                vel = U%data(ic,iumomx)/U%data(ic,iudens)

             case ("reflect")
                vel = -U%data(ir,iumomx)/U%data(ir,iudens)
                ir = ir + 1

             case default
                print *, "ERROR: invalid hse_vel_type"
                stop
                
             end select

             
             ! get the pressure above (we already know that it has the
             ! internal energy econst)
             call eos(eos_input_e, p_above, econst, U%data(i+1, iudens))

             ! p via HSE.  Here we explicitly make use of the
             ! gamma-law nature of the EOS.  Otherwise, we would need
             ! to iterate to simultaneously satisfy HSE + the EOS.
             p = (p_above - 0.5_dp_t*U%data(i+1,iudens)*grav*U%grid%dx) / &
                  (1.0_dp_t + 0.5_dp_t*grav*U%grid%dx/ &
                   (econst*(gamma - 1.0_dp_t)))

             ! now find the density that corresponds to this p, e
             call eos(eos_input_pe, p, econst, dens)
             
             U%data(i,iudens) = dens
             U%data(i,iumomx) = dens*vel

             ! compute E
             U%data(i,iuener) = dens*econst + &
                  0.5_dp_t*U%data(i,iumomx)**2/dens
             
          else

             print *, "invalid hse_bc_const type"
             stop
             
          endif
       enddo
    

    case ("outflow")

       do i = U%grid%lo-1, U%grid%lo-U%grid%ng, -1

          ! give all quantities a zero-gradient
          U%data(i,:) = U%data(i+1,:)
           
       enddo
       

    case ("user")

       call user_bc_xm(U)


    case ("periodic")
       
       ip = U%grid%hi
       do i = U%grid%lo-1, U%grid%lo-U%grid%ng, -1
          
          U%data(i,:) = U%data(ip,:)
          
          ip = ip - 1
          
       enddo

    case default
       
       print *, "ERROR: xlboundary not valid"
       stop
       
    end select
     

    !-------------------------------------------------------------------------
    ! upper boundary (+x)
    !-------------------------------------------------------------------------
    select case (U%grid%xrboundary)
       
    case ("reflect")
       
       ! ir is the index that corresponds to the reflected ghostcell
       ir = U%grid%hi
       do i = U%grid%hi+1, U%grid%hi+U%grid%ng
          
          ! even reflect density and energy
          U%data(i,iudens) = U%data(ir,iudens)
          U%data(i,iuener) = U%data(ir,iuener)
          
          ! reflect the velocities.  
          U%data(i,iumomx) = -U%data(ir,iumomx)
          
          ir = ir - 1
           
       enddo


    case ("hse")

       ir = U%grid%hi
       ic = U%grid%hi
        
       do i = U%grid%hi+1, U%grid%hi+U%grid%ng
          
          if (hse_bc_const == "density") then
              
             ! constant density in the ghost cells
              
             ! zero gradient to rho
             U%data(i,iudens) = U%data(ic,iudens)

             ! velocity depends on hse_vel_type
             select case (hse_vel_type)
             case ("outflow")
                vel = U%data(ic,iumomx)/U%data(ic,iudens)

             case ("reflect")
                vel = -U%data(ir,iumomx)/U%data(ir,iudens)
                ir = ir - 1

             case default
                print *, "ERROR: invalid hse_vel_type"
                stop
                
             end select

             U%data(i,iumomx) = U%data(i,iudens)*vel

             ! p via HSE
             e_below = (U%data(i-1,iuener) - &
                  0.5_dp_t*U%data(i-1,iumomx)**2/U%data(i-1,iudens))/ &
                  U%data(i-1,iudens)
             call eos(eos_input_e, p_below, e_below, U%data(i-1,iudens))
             
             p = p_below + 0.5_dp_t*(U%data(i,iudens) + &
                                     U%data(i-1,iudens))*grav*U%grid%dx
             
             call eos(eos_input_p, p, e, U%data(i,iudens))
             
             ! compute E
             U%data(i,iuener) = U%data(i,iudens)*e + &
                  0.5_dp_t*U%data(i,iumomx)**2/U%data(i,iudens)
             
          else if (hse_bc_const == "temperature") then
             
             ! constant temperature (or internal energy) in the ghost
             ! cells

             econst = (U%data(ic,iuener) - &
                  0.5_dp_t*U%data(ic,iumomx)**2/U%data(ic,iudens))/ &
                  U%data(ic,iudens)
             
             ! velocity depends on hse_vel_type
             select case(hse_vel_type)
             case ("outflow")
                vel = U%data(ic,iumomx)/U%data(ic,iudens)

             case ("reflect")
                vel = -U%data(ir,iumomx)/U%data(ir,iudens)
                ir = ir - 1

             case default
                print *, "ERROR: invalid hse_vel_type"
                stop
                
             end select

             
             ! get the pressure below (we already know that it has the
             ! internal energy econst)
             call eos(eos_input_e, p_below, econst, U%data(i-1, iudens))

             ! p via HSE.  Here we explicitly make use of the
             ! gamma-law nature of the EOS.  Otherwise, we would need
             ! to iterate to simultaneously satisfy HSE + the EOS.
             p = (p_below + 0.5_dp_t*U%data(i-1,iudens)*grav*U%grid%dx) / &
                  (1.0_dp_t - 0.5_dp_t*grav*U%grid%dx/ &
                   (econst*(gamma - 1.0_dp_t)))

             ! now find the density that corresponds to this p, e
             call eos(eos_input_pe, p, econst, dens)
             
             U%data(i,iudens) = dens
             U%data(i,iumomx) = dens*vel

             ! compute E
             U%data(i,iuener) = dens*econst + &
                  0.5_dp_t*U%data(i,iumomx)**2/dens
             
          else

             print *, "invalid hse_bc_const type"
             stop
             
          endif
       enddo


    case ("outflow")

       do i = U%grid%hi+1, U%grid%hi+U%grid%ng
          
          ! give all quantities a zero-gradient
          U%data(i,:) = U%data(i-1,:)
          
       enddo
       

    case ("diode")
       
       do i = U%grid%hi+1, U%grid%hi+U%grid%ng
          
          ! give all quantities a zero-gradient
          U%data(i,:) = U%data(i-1,:)
          
          ! store internal energy
          e = (U%data(i,iuener) - &
               0.5_dp_t*U%data(i,iumomx)**2/U%data(i,iudens))/U%data(i,iudens)
          
          U%data(i,iumomx) = max(0.0_dp_t, U%data(i,iumomx))
          
          ! recompute energy
          U%data(i,iuener) = U%data(i,iudens)*e + &
               0.5_dp_t*U%data(i,iumomx)**2/U%data(i,iudens)
          
       enddo


    case ("user")

       call user_bc_xp(U)
       
       
    case ("periodic")
       
       ip = U%grid%lo
       do i = U%grid%hi+1, U%grid%hi+U%grid%ng
          
          U%data(i,:) = U%data(ip,:)
          
          ip = ip + 1
          
       enddo
        
       
    case default
       
       print *, "ERROR: xrboundary not valid"
       stop
       
    end select
     
  end subroutine fillBCs

 end module bcs_module
