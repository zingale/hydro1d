module init_module

  use datatypes_module
  use params_module
  use probparams_module
  use grid_module
  use eos_module
  use variables_module

  implicit none

  private

  public :: init_data

contains
  
  subroutine init_data(U)

    use user_bc_module, only: e_fluff, rho_fluff

    type(gridvar_t), intent(inout) :: U

    integer :: i

    real (kind=dp_t) :: e_const, s_const
    real (kind=dp_t) :: dens_zone, pres_zone, e_zone
    real (kind=dp_t) :: drho, dens_zone_old

    real (kind=dp_t) :: e0, slope, x0

    integer :: iter

    integer :: icutoff

    logical :: first_transition, first_isentropic, cutoff

    real (kind=dp_t), parameter :: tol = 1.e-8_dp_t

    type(gridvar_t) :: p

    

    call build(p, U%grid, 1)

    ! base
    U%data(U%grid%lo,iudens) = dens_base
    p%data(U%grid%lo,1) = pres_base

    ! get the internal energy (temperature) for the isothermal region
    call eos(eos_input_p, pres_base, e_const, dens_base)

    ! store the first zone's state
    U%data(U%grid%lo,iumomx) = 0.0_dp_t    
    U%data(U%grid%lo,iuener) = dens_base*e_const

    first_isentropic = .true.
    first_transition = .true.

    cutoff = .false.

    do i = U%grid%lo+1, U%grid%hi

       ! initial guess for the density
       dens_zone = U%data(i-1,iudens)

       drho = 1.d33
       dens_zone_old = 1.d33

       iter = 0

       ! loop until density converges
       do while (abs(drho) > tol*dens_zone)

          ! what pressure does HSE want?
          pres_zone = p%data(i-1,1) + &
               0.5_dp_t*U%grid%dx*(dens_zone + U%data(i-1,iudens))*grav

          ! use the EOS to get the density for this pressure
          if (U%grid%x(i) < x_jump) then

             ! isothermal region -- constrain to e_const
             call eos(eos_input_pe, pres_zone, e_const, dens_zone)

          elseif (U%grid%x(i) < x_jump + delta) then

             ! transition region -- compute the e we want to constraint to

             if (first_transition) then
                e0 = U%data(i-1,iuener)/U%data(i-1,iudens)

                slope = (e_factor*e0 - e0)/delta
                x0 = U%grid%x(i-1)

                first_transition = .false.
             endif

             e_zone = slope*(U%grid%x(i) - x0) + e0

             call eos(eos_input_pe, pres_zone, e_zone, dens_zone)  
          else

             ! constant entropy region
             
             if (first_isentropic) then

                ! set the entropy to constrain to
                call eos(eos_input_p, p%data(i-1,1), e_zone, &
                         U%data(i-1,iudens), s=s_const)

                first_isentropic = .false.

             endif
             
             call eos(eos_input_ps, pres_zone, e_zone, dens_zone, &
                      s=s_const)

          endif
          
          drho = dens_zone - dens_zone_old
          dens_zone_old = dens_zone

          iter = iter + 1
          if (iter > 100) then
             print *, "ERROR: too many iterations"
             stop
          endif
          
       enddo

       ! store the state
       U%data(i,iudens) = dens_zone
       p%data(i,1) = pres_zone

       ! get the internal energy from the EOS
       call eos(eos_input_p, p%data(i,1), e_zone, U%data(i,iudens))

       ! store the remainder of the conserved state
       U%data(i,iumomx) = 0.0_dp_t
       U%data(i,iuener) = U%data(i,iudens)*e_zone    ! no kinetic energy


       ! check if density dropped below threshold
       if (dens_zone <= small_dens) then
          icutoff = i
          cutoff = .true.
          exit
       endif

    enddo

    if (cutoff) then
       U%data(icutoff+1:,iudens) = U%data(icutoff,iudens)
       U%data(icutoff+1:,iumomx) = U%data(icutoff,iumomx)
       U%data(icutoff+1:,iuener) = U%data(icutoff,iuener)
    endif

    do i = U%grid%lo, U%grid%hi
       print *, i, U%data(i,iudens), U%data(i,iuener)/U%data(i,iudens)
    enddo

    ! top state -- for the user boundary conditions
    rho_fluff = U%data(U%grid%hi,iudens)
    e_fluff = U%data(U%grid%hi,iuener)/rho_fluff


    ! HSE check
    ! do i = U%grid%lo+1, U%grid%hi
    !    if (U%data(i,iudens) <= small_dens) exit
       
    !    hse_err = abs( (p%data(i,1) - p%data(i-1,1))/U%grid%dx - &
    !                   0.5_dp_t*(U%data(i,iudens) + U%data(i-1,iudens))*grav) / &
    !                   abs(0.5_dp_t*(U%data(i,iudens) + U%data(i-1,iudens))*grav)

    !    print *, i, hse_err

    ! enddo


    call destroy(p)

  end subroutine init_data

end module init_module

