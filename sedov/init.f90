module init_module

  use datatypes_module
  use probparams_module
  use params_module
  use grid_module
  use variables_module

  implicit none

  private

  public :: init_data

contains
  
  subroutine init_data(U)

    type(gridvar_t), intent(inout) :: U

    integer :: i, ii
    real (kind=dp_t) :: rhoe

    real (kind=dp_t), parameter :: E_sedov = ONE
    real (kind=dp_t), parameter :: p_ambient = 1.d-5

    real (kind=dp_t) :: p_exp, V_exp, p_zone
    real (kind=dp_t) :: xl, xr, xx, vol_pert, vol_ambient
    integer, parameter :: nsub = 4

    ! we'll put all of the energy in the first zone
    if (single_point) then
       V_exp = 4.0*pi/3.0 * U%grid%dx**3

       do i = U%grid%lo, U%grid%hi
      
          ! store the conserved state
          U%data(i,iudens) = ONE
          U%data(i,iumomx) = ZERO
       
          if (i == U%grid%lo) then
             rhoe = E_sedov/V_exp
          else
             rhoe = p_ambient/(gamma - ONE)
          endif

          U%data(i,iuener) = rhoe          
       enddo

    else

       V_exp = 4.0*pi/3.0 * r_init**3
       p_exp = (gamma - 1.d0)*E_sedov/V_exp

       do i = U%grid%lo, U%grid%hi

          ! sub zone
          vol_pert = 0.0d0
          vol_ambient = 0.0d0

          do ii = 0, nsub-1
             xl = U%grid%xl(i) + U%grid%dx*dble(ii)/nsub
             xr = U%grid%xl(i) + U%grid%dx*dble(ii + 1.d0)/nsub
             xx = 0.5d0*(xl + xr)

                                                                                                                                                                          
             ! the volume of a subzone is (4/3) pi (xr^3 - xl^3).                                                                                                           
             ! we can factor this as: (4/3) pi dr (xr^2 + xl*xr + xl^2)                                                                                                     
             ! The (4/3) pi dr factor is common, so we can neglect it.                                                                                                      
             if(xx <= r_init) then                                                                                                                                          
                vol_pert = vol_pert + (xr*xr + xl*xr + xl*xl)                                                                                                               
             else                                                                                                                                                           
                vol_ambient = vol_ambient + (xr*xr + xl*xr + xl*xl)                                                                                                         
             endif
          enddo
 
          p_zone = (vol_pert*p_exp + vol_ambient*p_ambient)/(vol_pert+vol_ambient)                                                                                          

          ! store the conserved state
          U%data(i,iudens) = ONE
          U%data(i,iumomx) = ZERO
     
          rhoe = p_zone/(gamma - 1.d0)

          U%data(i,iuener) = rhoe          
       enddo

    endif

  end subroutine init_data

end module init_module

