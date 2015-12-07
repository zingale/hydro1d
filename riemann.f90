module riemann_module

  use datatypes_module
  use grid_module
  use variables_module

  implicit none

  private

  public :: solve_riemann

contains
  
  subroutine solve_riemann(Uin_l, Uin_r, fluxes, godunov_state)

    use params_module, only: gamma
    use eos_module

    type(gridedgevar_t), intent(in   ) :: Uin_l, Uin_r
    type(gridedgevar_t), intent(inout) :: fluxes
    type(gridedgevar_t), intent(inout) :: godunov_state

    !  Solve riemann shock tube problem for a general equation of
    !  state using the method of Colella, Glaz, and Ferguson.  See
    !  Almgren et al. 2010 (the CASTRO paper) for details.
    !
    !  The Riemann problem for the Euler's equation produces 4 states,
    !  separated by the three characteristics (u - cs, u, u + cs):
    !
    !
    !       u - cs    t    u      u + cs
    !         \       ^   .       /
    !          \  *L  |   . *R   /
    !           \     |  .     /
    !            \    |  .    /
    !        L    \   | .   /    R
    !              \  | .  /
    !               \ |. /
    !                \|./
    !       ----------+----------------> x
    !
    !       only density jumps across the u characteristic
    !

    real (kind=dp_t) :: rho_l, u_l, rhoe_l, p_l
    real (kind=dp_t) :: rho_r, u_r, rhoe_r, p_r
    real (kind=dp_t) :: e_l, e_r

    real (kind=dp_t) :: W_l, W_r, c_l, c_r

    real (kind=dp_t) :: pstar, ustar
    real (kind=dp_t) :: rhostar_l, rhostar_r, rhoestar_l, rhoestar_r
    real (kind=dp_t) :: cstar_l, cstar_r

    real (kind=dp_t) :: lambda_l, lambda_r
    real (kind=dp_t) :: lambdastar_l, lambdastar_r
    real (kind=dp_t) :: sigma, alpha

    real (kind=dp_t) :: p_state, u_state, rhoe_state, rho_state

    real (kind=dp_t), parameter :: smallc   = 1.e-10_dp_t
    real (kind=dp_t), parameter :: smallrho = 1.e-10_dp_t
    real (kind=dp_t), parameter :: smallp   = 1.e-10_dp_t

    integer :: i

    do i = Uin_l%grid%lo, Uin_l%grid%hi+1

       ! left state
       rho_l  = Uin_l%data(i,iudens)
       u_l    = Uin_l%data(i,iumomx)/rho_l
       rhoe_l = Uin_l%data(i,iuener) - 0.5_dp_t*rho_l*u_l**2
       e_l = rhoe_l/rho_l
       call eos(eos_input_e, p_l, e_l, rho_l)   ! get p_l

       ! right state
       rho_r  = Uin_r%data(i,iudens)
       u_r    = Uin_r%data(i,iumomx)/rho_r
       rhoe_r = Uin_r%data(i,iuener) - 0.5_dp_t*rho_r*u_r**2
       e_r = rhoe_r/rho_r
       call eos(eos_input_e, p_r, e_r, rho_r)   ! get p_r


       ! define the Lagrangian sound speed
       W_l = max(smallrho*smallc, sqrt(gamma*p_l*rho_l))
       W_r = max(smallrho*smallc, sqrt(gamma*p_r*rho_r))
       
       ! and the regular sound speeds
       c_l = max(smallc, sqrt(gamma*p_l/rho_l))
       c_r = max(smallc, sqrt(gamma*p_r/rho_r))
       
       ! define the star states
       pstar = (W_l*p_r + W_r*p_l + W_l*W_r*(u_l - u_r))/(W_l + W_r)
       pstar = max(pstar, smallp)
       ustar = (W_l*u_l + W_r*u_r + (p_l - p_r))/(W_l + W_r)
       
       ! now compute the remaining state to the left and right of the
       ! contact (in the star region)
       rhostar_l = rho_l + (pstar - p_l)/c_l**2
       rhostar_r = rho_r + (pstar - p_r)/c_r**2
       
       rhoestar_l = rhoe_l + (pstar - p_l)*(rhoe_l/rho_l + p_l/rho_l)/c_l**2
       rhoestar_r = rhoe_r + (pstar - p_r)*(rhoe_r/rho_r + p_r/rho_r)/c_r**2

       cstar_l = max(smallc,sqrt(gamma*pstar/rhostar_l))
       cstar_r = max(smallc,sqrt(gamma*pstar/rhostar_r))
       
       ! figure out which state we are in, based on the location of
       ! the waves
       if (ustar > 0_dp_t) then

          ! contact is moving right, we need to understand the L and
          ! *L states

          ! define eigenvalues
          lambda_l = u_l - c_l
          lambdastar_l = ustar - cstar_l
          
          if (pstar > p_l) then
             ! the wave is a shock -- find the shock speed
             sigma = (lambda_l + lambdastar_l)/2.0_dp_t
             
             if (sigma > 0.0_dp_t) then
                ! shock is moving to the right -- solution is L state
                rho_state = rho_l
                u_state = u_l
                p_state = p_l
                rhoe_state = rhoe_l
                
             else
                ! solution is *L state
                rho_state = rhostar_l
                u_state = ustar
                p_state = pstar
                rhoe_state = rhoestar_l
                
             endif

          else
             ! the wave is a rarefaction
             if (lambda_l < 0.0_dp_t .and. lambdastar_l < 0.0_dp_t) then
                ! rarefaction fan is moving to the left -- solution is
                ! *L state
                rho_state = rhostar_l
                u_state = ustar
                p_state = pstar
                rhoe_state = rhoestar_l
                
             else if (lambda_l > 0.0_dp_t .and. lambdastar_l > 0.0_dp_t) then
                ! rarefaction fan is moving to the right -- solution
                ! is L state
                rho_state = rho_l
                u_state = u_l
                p_state = p_l
                rhoe_state = rhoe_l
                
             else
                ! rarefaction spans x/t = 0 -- interpolate
                alpha = lambda_l/(lambda_l - lambdastar_l)
                
                rho_state  = alpha*rhostar_l  + (1.0_dp_t - alpha)*rho_l
                u_state    = alpha*ustar      + (1.0_dp_t - alpha)*u_l
                p_state    = alpha*pstar      + (1.0_dp_t - alpha)*p_l
                rhoe_state = alpha*rhoestar_l + (1.0_dp_t - alpha)*rhoe_l

             endif

          endif
                          

       else if (ustar < 0_dp_t) then
          
          ! contact moving left, we need to understand the R and *R
          ! states

          ! define eigenvalues
          lambda_r = u_r + c_r
          lambdastar_r = ustar + cstar_r
          
          if (pstar > p_r) then
             ! the wave if a shock -- find the shock speed
             sigma = (lambda_r + lambdastar_r)/2.0_dp_t
             
             if (sigma > 0.0_dp_t) then
                ! shock is moving to the right -- solution is *R state
                rho_state = rhostar_r
                u_state = ustar
                p_state = pstar
                rhoe_state = rhoestar_r
                
             else
                ! solution is R state
                rho_state = rho_r
                u_state = u_r
                p_state = p_r
                rhoe_state = rhoe_r
                
             endif

          else
             ! the wave is a rarefaction
             if (lambda_r < 0.0_dp_t .and. lambdastar_r < 0.0_dp_t) then
                ! rarefaction fan is moving to the left -- solution is
                ! R state
                rho_state = rho_r
                u_state = u_r
                p_state = p_r
                rhoe_state = rhoe_r
                
             else if (lambda_r > 0.0_dp_t .and. lambdastar_r > 0.0_dp_t) then
                ! rarefaction fan is moving to the right -- solution
                ! is *R state
                rho_state = rhostar_r
                u_state = ustar
                p_state = pstar
                rhoe_state = rhoestar_r
                
             else
                ! rarefaction spans x/t = 0 -- interpolate
                alpha = lambda_r/(lambda_r - lambdastar_r)
                
                rho_state  = alpha*rhostar_r  + (1.0_dp_t - alpha)*rho_r
                u_state    = alpha*ustar      + (1.0_dp_t - alpha)*u_r
                p_state    = alpha*pstar      + (1.0_dp_t - alpha)*p_r
                rhoe_state = alpha*rhoestar_r + (1.0_dp_t - alpha)*rhoe_r

             endif

          endif



       else ! ustar == 0

          rho_state = 0.5_dp_t*(rhostar_l + rhostar_r)
          u_state = ustar
          p_state = pstar
          rhoe_state = 0.5_dp_t*(rhoestar_l + rhoestar_r)

       endif


       ! reflect BC hack
       !if (i == Uin_l%grid%lo .and. Uin_l%grid%xlboundary == "reflect") then
       !   u_state = 0.0_dp_t
       !endif

       ! compute the fluxes
       fluxes%data(i,iudens) = rho_state*u_state

       ! for spherical, don't include the pressure term, since it is
       ! a gradient, not a divergence
       if (Uin_l%grid%geometry == 1) then
          fluxes%data(i,iumomx) = rho_state*u_state*u_state
       else
          fluxes%data(i,iumomx) = rho_state*u_state*u_state + p_state
       endif

       fluxes%data(i,iuener) = rhoe_state*u_state + &
            0.5_dp_t*rho_state*u_state**3 + p_state*u_state

       ! we need the pressure on the interface for Riemann problems
       godunov_state%data(i,iqdens) = rho_state
       godunov_state%data(i,iqxvel) = u_state
       godunov_state%data(i,iqpres) = p_state


    enddo


    
  end subroutine solve_riemann

end module riemann_module
