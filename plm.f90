module interface_states_plm_module

  use datatypes_module
  use grid_module
  use variables_module
  use params_module
  use eos_module

  implicit none

contains
  
  subroutine make_interface_states_plm(U, U_l, U_r, dt)

    type(gridvar_t),     intent(in   ) :: U
    type(gridedgevar_t), intent(inout) :: U_l, U_r
    real (kind=dp_t),    intent(in   ) :: dt

    type(gridvar_t) :: Q
    type(gridedgevar_t) :: Q_l, Q_r
    type(gridvar_t) :: ldelta, temp

    real (kind=dp_t) :: rvec(nprim,nprim), lvec(nprim,nprim), eval(nprim)
    real (kind=dp_t) :: dQ(nprim)

    real (kind=dp_t) :: r, ux, p, cs
    real (kind=dp_t) :: ldr, ldu, ldp
    real (kind=dp_t) :: r_xm, r_xp, u_xm, u_xp, p_xm, p_xp
    real (kind=dp_t) :: beta_xm(nprim), beta_xp(nprim)
    real (kind=dp_t) :: sum_xm, sum_xp
    real (kind=dp_t) :: sum

    real (kind=dp_t) :: dp, dp2, z
    type(gridvar_t) :: xi_t, xi
    real (kind=dp_t), parameter :: z0 = 0.75_dp_t
    real (kind=dp_t), parameter :: z1 = 0.85_dp_t
    real (kind=dp_t), parameter :: delta = 0.33_dp_t
    real (kind=dp_t), parameter :: smallp = 1.e-10_dp_t
    

    real (kind=dp_t) :: dtdx

    real (kind=dp_t) :: e
    real (kind=dp_t) :: test

    integer :: i, m, n

    ! piecewise linear slopes 
    !
    ! This is a 1-d version of the piecewise linear Godunov method
    ! detailed in Colella (1990).  See also Colella & Glaz and
    ! Saltzman (1994).
    !
    ! We wish to solve
    !
    !   U_t + [F(U)]_x = H
    !
    ! we want U_{i+1/2}^{n+1/2} -- the interface values that are input
    ! to the Riemann problem through the faces for each zone.
    !
    ! Taylor expanding yields
    !
    !  n+1/2                     dU           dU  
    ! U          = U   + 0.5 dx  --  + 0.5 dt --  
    !  i+1/2,j,L    i,j          dx           dt   
    ! 
    !   
    !                            dU             dF 
    !            = U   + 0.5 dx  --  - 0.5 dt ( -- - H ) 
    !               i,j          dx             dx  
    !      
    !       
    !                             dU      dF
    !            = U   + 0.5 ( dx -- - dt -- ) + 0.5 dt H 
    !               i,j           dx      dx   
    !  
    !       
    !                                 dt       dU 
    !            = U   + 0.5 dx ( 1 - -- A^x ) -- + 0.5 dt H  
    !               i,j               dx       dx      
    !   
    !   
    !                               dt       _ 
    !            = U   + 0.5  ( 1 - -- A^x ) DU + 0.5 dt H   
    !               i,j             dx     
    !     
    !                    +----------+----------+ +---+---+    
    !                               |                |   
    !               
    !                   this is the monotonized   source term    
    !                   central difference term     
  

    
    ! sanity check
    if (U%grid%ng < 4) then
       print *, "ERROR: ng < 4 in plm states"
       stop
    endif


    !-------------------------------------------------------------------------
    ! convert to primitve variables
    !-------------------------------------------------------------------------
    call build(Q, U%grid, nprim)

    do i = U%grid%lo-U%grid%ng, U%grid%hi+U%grid%ng
       ! density
       Q%data(i,iqdens) = U%data(i,iudens)

       ! velocity
       Q%data(i,iqxvel) = U%data(i,iumomx)/U%data(i,iudens)

       ! pressure
       e = (U%data(i,iuener) - &
            0.5_dp_t*U%data(i,iumomx)**2/U%data(i,iudens))/U%data(i,iudens)
       call eos(eos_input_e, Q%data(i,iqpres), e, Q%data(i,iqdens))
    enddo


    !-------------------------------------------------------------------------
    ! compute the flattening coefficients
    !-------------------------------------------------------------------------

    ! flattening kicks in behind strong shocks.  See Saltzman (1994) page 159
    ! for this implementation.
    call build(xi_t, U%grid, 1)

    do i = U%grid%lo-2, U%grid%hi+2
       
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

    call build(xi, U%grid, 1)

    do i = U%grid%lo-1, U%grid%hi+1

       if (Q%data(i+1,iqpres) - Q%data(i-1,iqpres) > 0.0_dp_t) then
          xi%data(i,1) = min(xi_t%data(i,1), xi_t%data(i-1,1))

       else if (Q%data(i+1,iqpres) - Q%data(i-1,iqpres) < 0.0_dp_t) then
          xi%data(i,1) = min(xi_t%data(i,1), xi_t%data(i+1,1))

       else
          xi%data(i,1) = xi_t%data(i,1)
       endif

    enddo
    
    
    call destroy(xi_t)

    

    !-------------------------------------------------------------------------
    ! compute the monotonized central differences
    !-------------------------------------------------------------------------

    ! 4th order MC limiting.  See Colella (1985) Eq. 2.5 and 2.6,
    ! Colella (1990) page 191 (with the delta a terms all equal) or
    ! Saltzman 1994, page 156
    
    call build(temp, U%grid, 1)
    call build(ldelta, U%grid, nprim)

    do n = 1, nprim

       ! first do the normal MC limiting 
       do i = Q%grid%lo-3, Q%grid%hi+3

          test = (Q%data(i+1,n) - Q%data(i,n))*(Q%data(i,n) - Q%data(i-1,n))

          if (test > 0.0_dp_t) then
             temp%data(i,1) = min(0.5_dp_t*abs(Q%data(i+1,n) - Q%data(i-1,n)), &
                                  min(2.0_dp_t*abs(Q%data(i+1,n) - Q%data(i,n)), &
                                      2.0_dp_t*abs(Q%data(i,n) - Q%data(i-1,n)))) * &
                              sign(1.0_dp_t,Q%data(i+1,n) - Q%data(i-1,n))

          else
             temp%data(i,1) = 0.0_dp_t
          endif

       enddo
    

       ! now do the fourth order part
       do i = Q%grid%lo-2, Q%grid%hi+2
          
          test = (Q%data(i+1,n) - Q%data(i,n))*(Q%data(i,n) - Q%data(i-1,n))

          if (test > 0.0_dp_t) then
             ldelta%data(i,n) = &
                  min((2.0_dp_t/3.0_dp_t)*abs(Q%data(i+1,n) - Q%data(i-1,n) - &
                      0.25_dp_t*(temp%data(i+1,1) + temp%data(i-1,1))), &
                  min(2.0*abs(Q%data(i+1,n) - Q%data(i,n)), &
                      2.0*abs(Q%data(i,n) - Q%data(i-1,n)))) * &
                  sign(1.0_dp_t, Q%data(i+1,n) - Q%data(i-1,n))
          else
             ldelta%data(i,n) = 0.0_dp_t
          endif

       enddo

    enddo


    ! apply flattening to the slopes
    do n = 1, nprim
       do i = Q%grid%lo-2, Q%grid%hi+2
          ldelta%data(i,n) = xi%data(i,1)*ldelta%data(i,n)
       enddo
    enddo

    call destroy(temp)
    call destroy(xi)


    !-------------------------------------------------------------------------
    ! compute left and right primitive variable states
    !-------------------------------------------------------------------------
    call build(Q_l, U%grid, nprim)
    call build(Q_r, U%grid, nprim)

    dtdx = dt/U%grid%dx

    ! We need the left and right eigenvectors and the eigenvalues for
    ! the system
    !
    ! The eigenvalues are u - c, u, u + c
    !
    ! The Jacobian matrix for the primitive variable formulation of the
    ! Euler equations is 
    !
    !       / u   r   0   \  
    !   A = | 0   u   1/r |  
    !       \ 0  rc^2 u  /
    !
    ! With the rows corresponding to rho, u, and p
    !
    ! The right eigenvectors are
    !
    !       /  1  \        / 1 \        /  1  \
    ! r1 =  |-c/r |   r2 = | 0 |   r3 = | c/r |
    !       \ c^2 /        \ 0 /        \ c^2 /
    !
    !
    ! The left eigenvectors are    
    !
    !   l1 =     ( 0,  -r/(2c),  1/(2c^2) )
    !   l2 =     ( 1,     0,     -1/c^2,  )
    !   l3 =     ( 0,   r/(2c),  1/(2c^2) )
    !
    ! The fluxes are going to be defined on the left edge of the
    ! computational zone.
    !
    !          |             |             |             | 
    !          |             |             |             |  
    !         -+------+------+------+------+------+------+--   
    !          |     i-1     |      i      |     i+1     |   
    !                       * *           * 
    !                   V_l,i  V_r,i   V_l,i+1 
    !           
    ! V_l,i+1 are computed using the information in zone i,j.
    !

    ! the basic idea here is that we do a characteristic
    ! decomposition.  The jump in primitive variables (Q)
    ! can be transformed to a jump in characteristic variables
    ! using the left and right eigenvectors.  Then each wave
    ! tells us how much of each characteristic quantity reaches
    ! the interface over dt/2.  We only add the quantity if
    ! it moves toward the interface.
    !
    ! Given a jump dQ, we can express this as:
    !
    !       3
    ! dQ = sum  (l_i dQ) r_i
    !      i=1
    !
    ! Then 
    !         3
    ! A dQ = sum  lambda_i (l_i dQ) r_i
    !        i=1
    !
    ! where lambda_i is the eigenvalue.
    !

    do i = U%grid%lo-1, U%grid%hi+1

       r  = Q%data(i,iqdens)
       ux = Q%data(i,iqxvel)
       p  = Q%data(i,iqpres)

       ldr = ldelta%data(i,iqdens)
       ldu = ldelta%data(i,iqxvel)
       ldp = ldelta%data(i,iqpres)

       dQ(1) = ldr
       dQ(2) = ldu
       dQ(3) = ldp


       ! compute the sound speed
       cs = sqrt(gamma*p/r)

       ! compute the eigenvalues
       eval(1) = ux - cs
       eval(2) = ux
       eval(3) = ux + cs


       ! compute the left eigenvectors
       
       ! u - c eigenvector
       lvec(1,1) = 0.0_dp_t
       lvec(1,2) = -0.5_dp_t*r/cs
       lvec(1,3) = 0.5_dp_t/(cs*cs)

       ! u eigenvector
       lvec(2,1) = 1.0_dp_t
       lvec(2,2) = 0.0_dp_t
       lvec(2,3) = -1.0_dp_t/(cs*cs)

       ! u + c eigenvector
       lvec(3,1) = 0.0_dp_t
       lvec(3,2) = 0.5_dp_t*r/cs
       lvec(3,3) = 0.5_dp_t/(cs*cs)


       ! compute the right eigenvectors

       ! u - c eigenvector
       rvec(1,1) = 1.0_dp_t
       rvec(1,2) = -cs/r
       rvec(1,3) = cs*cs

       ! u eigenvector
       rvec(2,1) = 1.0_dp_t
       rvec(2,2) = 0.0_dp_t
       rvec(2,3) = 0.0_dp_t

       ! u + c eigenvector
       rvec(3,1) = 1.0_dp_t
       rvec(3,2) = cs/r
       rvec(3,3) = cs*cs


       ! Define the reference states (here xp is the right interface
       ! for the current zone and xm is the left interface for the
       ! current zone)
       !                           ~
       ! These expressions are the V_{L,R} in Colella (1990) at the 
       ! bottom of page 191.  They are also in Saltzman (1994) as
       ! V_ref on page 161.
       r_xp = r + 0.5_dp_t*(1.0_dp_t - dtdx*max(eval(3), 0.0_dp_t))*ldr
       r_xm = r - 0.5_dp_t*(1.0_dp_t + dtdx*min(eval(1), 0.0_dp_t))*ldr

       u_xp = ux + 0.5_dp_t*(1.0_dp_t - dtdx*max(eval(3), 0.0_dp_t))*ldu
       u_xm = ux - 0.5_dp_t*(1.0_dp_t + dtdx*min(eval(1), 0.0_dp_t))*ldu

       p_xp = p + 0.5_dp_t*(1.0_dp_t - dtdx*max(eval(3), 0.0_dp_t))*ldp
       p_xm = p - 0.5_dp_t*(1.0_dp_t + dtdx*min(eval(1), 0.0_dp_t))*ldp

       !                                                   ^
       ! Now compute the interface states.   These are the V expressions
       ! in Colella (1990) page 191, and the interface state expressions
       ! -V and +V in Saltzman (1994) on pages 161.
       
       ! first compute beta_xm and beta_xp
       do m = 1, nprim
          sum = 0.0_dp_t

          do n = 1, nprim
             sum = sum + lvec(m,n)*dQ(n)
          enddo

          ! here the sign() function makes sure we only add the right-moving
          ! waves
          beta_xp(m) = 0.25_dp_t*dtdx*(eval(3) - eval(m))* &
               (sign(1.0_dp_t, eval(m)) + 1.0_dp_t)*sum

          ! here the sign() function makes sure we only add the left-moving
          ! waves
          beta_xm(m) = 0.25_dp_t*dtdx*(eval(1) - eval(m))* &
               (1.0_dp_t - sign(1.0_dp_t, eval(m)))*sum

       enddo

       
       ! density
       sum_xm = 0.0_dp_t
       sum_xp = 0.0_dp_t
       do n = 1, nprim
          sum_xm = sum_xm + beta_xm(n)*rvec(n,1)
          sum_xp = sum_xp + beta_xp(n)*rvec(n,1)
       enddo

       Q_l%data(i+1,iqdens) = r_xp + sum_xp
       Q_r%data(i,iqdens) = r_xm + sum_xm


       ! velocity
       sum_xm = 0.0_dp_t
       sum_xp = 0.0_dp_t
       do n = 1, nprim
          sum_xm = sum_xm + beta_xm(n)*rvec(n,2)
          sum_xp = sum_xp + beta_xp(n)*rvec(n,2)
       enddo

       Q_l%data(i+1,iqxvel) = u_xp + sum_xp
       Q_r%data(i,iqxvel) = u_xm + sum_xm


       ! pressure
       sum_xm = 0.0_dp_t
       sum_xp = 0.0_dp_t
       do n = 1, nprim
          sum_xm = sum_xm + beta_xm(n)*rvec(n,3)
          sum_xp = sum_xp + beta_xp(n)*rvec(n,3)
       enddo

       Q_l%data(i+1,iqpres) = p_xp + sum_xp
       Q_r%data(i,iqpres) = p_xm + sum_xm

       
    enddo

    ! clean-up
    call destroy(ldelta)


    !-------------------------------------------------------------------------
    ! transform the states into conserved variables
    !-------------------------------------------------------------------------
    do i = U%grid%lo, U%grid%hi+1

       ! density
       U_l%data(i,iudens) = Q_l%data(i,iqdens) 
       U_r%data(i,iudens) = Q_r%data(i,iqdens) 

       ! momentum
       U_l%data(i,iumomx) = Q_l%data(i,iqdens)*Q_l%data(i,iqxvel)
       U_r%data(i,iumomx) = Q_r%data(i,iqdens)*Q_r%data(i,iqxvel)

       ! total energy
       call eos(eos_input_p, Q_l%data(i,iqpres), e, Q_l%data(i,iqdens))
       U_l%data(i,iuener) = Q_l%data(i,iqdens)*e + &
            0.5_dp_t*Q_l%data(i,iqdens)*Q_l%data(i,iqxvel)**2

       call eos(eos_input_p, Q_r%data(i,iqpres), e, Q_r%data(i,iqdens))
       U_r%data(i,iuener) = Q_r%data(i,iqdens)*e + &
            0.5_dp_t*Q_r%data(i,iqdens)*Q_r%data(i,iqxvel)**2

    enddo
    

    ! clean-up
    call destroy(Q)
    call destroy(Q_l)
    call destroy(Q_r)


    !-------------------------------------------------------------------------
    ! apply the source terms
    !-------------------------------------------------------------------------


  end subroutine make_interface_states_plm

end module interface_states_plm_module
