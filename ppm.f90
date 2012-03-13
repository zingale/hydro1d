module interface_states_ppm_module

  use datatypes_module
  use grid_module
  use variables_module
  use params_module
  use eos_module

  implicit none

  private

  public :: make_interface_states_ppm

contains
  
  subroutine make_interface_states_ppm(U, U_l, U_r, dt)

    type(gridvar_t),     intent(in   ) :: U
    type(gridedgevar_t), intent(inout) :: U_l, U_r
    real (kind=dp_t),    intent(in   ) :: dt

    type(gridvar_t) :: Q
    type(gridedgevar_t) :: Q_l, Q_r

    real (kind=dp_t) :: rvec(nprim,nprim), lvec(nprim,nprim), eval(nprim)

    real (kind=dp_t) :: r, ux, p, cs
    real (kind=dp_t) :: beta_xm(nprim), beta_xp(nprim)

    real (kind=dp_t) :: dq0, dqp
    real (kind=dp_t) :: Iminus(nprim,nprim), Iplus(nprim,nprim)

    type(gridedgevar_t) :: Qminus, Qplus
    type(gridvar_t) :: Q6

    real (kind=dp_t) :: sigma

    real (kind=dp_t) :: dp, dp2, z
    type(gridvar_t) :: xi_t, xi
    real (kind=dp_t), parameter :: z0 = 0.75_dp_t
    real (kind=dp_t), parameter :: z1 = 0.85_dp_t
    real (kind=dp_t), parameter :: delta = 0.33_dp_t
    real (kind=dp_t), parameter :: smallp = 1.e-10_dp_t
    

    real (kind=dp_t) :: dtdx

    real (kind=dp_t) :: e

    integer :: i, m, n

    ! piecewise parabolic slopes 
    !
    ! This is a 1-d version of the piecewise parabolic method detailed
    ! Colella & Woodward (1984).  We follow the description of Almgren
    ! et al. 2010 (the CASTRO paper).
    
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
    ! interpolate the cell-centered data to the edges
    !-------------------------------------------------------------------------
    
    ! For each cell, we will find the Qminus and Qplus states -- these
    ! are the - and + edges of the cell.
    call build(Qminus, U%grid, nprim)
    call build(Qplus,  U%grid, nprim)
    
    do n = 1, nprim
       do i = U%grid%lo-2, U%grid%hi+2

          ! dq (C&W Eq. 1.7)
          dq0 = 0.5_dp_t*(Q%data(i+1,n) - Q%data(i-1,n))
          dqp = 0.5_dp_t*(Q%data(i+2,n) - Q%data(i,n))

          ! limiting (C&W Eq. 1.8)
          if ( (Q%data(i+1,n) - Q%data(i,n))* &
               (Q%data(i,n) - Q%data(i-1,n)) > 0.0_dp_t) then
             dq0 = sign(1.0_dp_t,dq0)* &
                  min( abs(dq0), &
                       2.0_dp_t*abs(Q%data(i,n) - Q%data(i-1,n)), &
                       2.0_dp_t*abs(Q%data(i+1,n) - Q%data(i,n)) )
          else
             dq0 = 0.0_dp_t
          endif

          if ( (Q%data(i+2,n) - Q%data(i+1,n))* &
               (Q%data(i+1,n) - Q%data(i,n)) > 0.0_dp_t) then
             dqp = sign(1.0_dp_t,dqp)* &
                  min( abs(dqp), &
                       2.0_dp_t*abs(Q%data(i+1,n) - Q%data(i,n)), &
                       2.0_dp_t*abs(Q%data(i+2,n) - Q%data(i+1,n)) )
          else
             dqp = 0.0_dp_t
          endif

          ! cubic (C&W Eq. 1.6)
          Qplus%data(i,n) = 0.5_dp_t*(Q%data(i,n) + Q%data(i+1,n)) - &
               (1.0_dp_t/6.0_dp_t)*(dqp - dq0)

          Qminus%data(i+1,n) = Qplus%data(i,n)

          ! make sure that we didn't over or undersoot -- this may not be needed, but
          ! is discussed in Colella & Sekora (2008)
          Qplus%data(i,n) = max(Qplus%data(i,n), min(Q%data(i,n),Q%data(i+1,n)))
          Qplus%data(i,n) = min(Qplus%data(i,n), max(Q%data(i,n),Q%data(i+1,n)))

          Qminus%data(i+1,n) = max(Qminus%data(i+1,n), min(Q%data(i,n),Q%data(i+1,n)))
          Qminus%data(i+1,n) = min(Qminus%data(i+1,n), max(Q%data(i,n),Q%data(i+1,n)))

       enddo
    enddo


    !-------------------------------------------------------------------------
    ! construct the parameters for the parabolic reconstruction polynomials
    !-------------------------------------------------------------------------

    ! Limit (C&W Eq. 1.10).  Here the loop is over cells, and
    ! considers the values on either side of the center of the cell
    ! (Qminus and Qplus).
    do n = 1, nprim
       do i = U%grid%lo-1, U%grid%hi+1

          if ( (Qplus%data(i,n) - Q%data(i,n)) * &
               (Q%data(i,n) - Qminus%data(i,n)) <= 0.0_dp_t) then
             Qminus%data(i,n) = Q%data(i,n)
             Qplus%data(i,n) = Q%data(i,n)

          else if ( (Qplus%data(i,n) - Qminus%data(i,n)) * &
                    (Q%data(i,n) - &
                      0.5_dp_t*(Qminus%data(i,n) + Qplus%data(i,n))) > &
                   (Qplus%data(i,n) - Qminus%data(i,n))**2/6.0_dp_t ) then

          ! alternate test from Colella & Sekora (2008)
          !else if (abs(Qminus%data(i,n) - Q%data(i,n)) >= &
          !     2.0*abs(Qplus%data(i,n) - Q%data(i,n))) then
             Qminus%data(i,n) = 3.0_dp_t*Q%data(i,n) - 2.0_dp_t*Qplus%data(i,n)

          else if (-(Qplus%data(i,n) - Qminus%data(i,n))**2/6.0_dp_t > &
                    (Qplus%data(i,n) - Qminus%data(i,n)) * &
                    (Q%data(i,n) - &
                           0.5_dp_t*(Qminus%data(i,n) + Qplus%data(i,n))) ) then

          !else if (abs(Qplus%data(i,n) - Q%data(i,n)) >= &
          !     2.0*abs(Qminus%data(i,n) - Q%data(i,n))) then
             Qplus%data(i,n) = 3.0_dp_t*Q%data(i,n) - 2.0_dp_t*Qminus%data(i,n)

          endif

       enddo
    enddo

    ! define Q6
    call build(Q6,  U%grid, nprim)

    do n = 1, nprim
       do i = U%grid%lo-1, U%grid%hi+1

          Q6%data(i,n) = 6.0_dp_t*Q%data(i,n) - &
               3.0_dp_t*(Qminus%data(i,n) + Qplus%data(i,n))

       enddo
    enddo



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

    do i = U%grid%lo-1, U%grid%hi+1

       r  = Q%data(i,iqdens)
       ux = Q%data(i,iqxvel)
       p  = Q%data(i,iqpres)


       ! compute the sound speed
       cs = sqrt(gamma*p/r)

       ! compute the eigenvalues
       eval(1) = ux - cs
       eval(2) = ux
       eval(3) = ux + cs


       ! compute the left eigenvectors
       lvec(1,:) = [ 0.0_dp_t, -0.5_dp_t*r/cs, 0.5_dp_t/(cs*cs)  ]   ! u - c
       lvec(2,:) = [ 1.0_dp_t, 0.0_dp_t,       -1.0_dp_t/(cs*cs) ]   ! u
       lvec(3,:) = [ 0.0_dp_t, 0.5_dp_t*r/cs,  0.5_dp_t/(cs*cs)  ]   ! u + c

       ! compute the right eigenvectors
       rvec(1,:) = [ 1.0_dp_t, -cs/r,    cs*cs    ]   ! u - c 
       rvec(2,:) = [ 1.0_dp_t, 0.0_dp_t, 0.0_dp_t ]   ! u  
       rvec(3,:) = [ 1.0_dp_t, cs/r,     cs*cs    ]   ! u + c   


       ! integrate the parabola in the cell from the left interface
       ! (Iminus) over the portion of the cell that each eigenvalue
       ! can reach.  Do the same from the right interface in the
       ! cell, defining Iplus.  See Almgren et al. 2010 (Eq. 30) or
       ! Colella & Sekora (2008)
       do m = 1, nprim
          sigma = abs(eval(m))*dtdx
          do n = 1, nprim

             Iplus(m,n) = Qplus%data(i,n) - 0.5_dp_t*sigma* &
                  (Qplus%data(i,n) - Qminus%data(i,n) - &
                  (1.0_dp_t - (2.0_dp_t/3.0_dp_t)*sigma)*Q6%data(i,n))

             Iminus(m,n) = Qminus%data(i,n) + 0.5_dp_t*sigma* &
                  (Qplus%data(i,n) - Qminus%data(i,n) + &
                  (1.0_dp_t - (2.0_dp_t/3.0_dp_t)*sigma)*Q6%data(i,n))
             
          enddo
       enddo


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
       ! where lambda_i is the eigenvalue.  (See for example LeVeque's
       ! book).
       !
       
       ! compute the dot product of each left eigenvector with (q - I)
       do m = 1, nprim
          beta_xm(m) = 0.0_dp_t
          beta_xp(m) = 0.0_dp_t

          do n = 1, nprim
             beta_xm(m) = beta_xm(m) + lvec(m,n)*(Q%data(i,n) - Iminus(m,n))
             beta_xp(m) = beta_xp(m) + lvec(m,n)*(Q%data(i,n) - Iplus(m,n))
          enddo
       enddo


       ! density
       Q_l%data(i+1,iqdens) = 0.0_dp_t
       Q_r%data(i,iqdens) = 0.0_dp_t

       do n = 1, nprim
          if (eval(n) >= 0.0_dp_t) then
             Q_l%data(i+1,iqdens) = Q_l%data(i+1,iqdens) + &
                  beta_xp(n)*rvec(n,1)
          endif

          if (eval(n) <= 0.0_dp_t) then
             Q_r%data(i,iqdens) = Q_r%data(i,iqdens) + &
                  beta_xm(n)*rvec(n,1)
          endif
       enddo

       Q_l%data(i+1,iqdens) = Q%data(i,iqdens) - &
            xi%data(i,1)*Q_l%data(i+1,iqdens)

       Q_r%data(i,iqdens) = Q%data(i,iqdens) - &
            xi%data(i,1)*Q_r%data(i,iqdens)


       ! velocity
       Q_l%data(i+1,iqxvel) = 0.0_dp_t
       Q_r%data(i,iqxvel) = 0.0_dp_t

       do n = 1, nprim
          if (eval(n) >= 0.0_dp_t) then
             Q_l%data(i+1,iqxvel) = Q_l%data(i+1,iqxvel) + &
                  beta_xp(n)*rvec(n,2)
          endif

          if (eval(n) <= 0.0_dp_t) then
             Q_r%data(i,iqxvel) = Q_r%data(i,iqxvel) + &
                  beta_xm(n)*rvec(n,2)
          endif
       enddo

       Q_l%data(i+1,iqxvel) = Q%data(i,iqxvel) - &
            xi%data(i,1)*Q_l%data(i+1,iqxvel)

       Q_r%data(i,iqxvel) = Q%data(i,iqxvel) - &
            xi%data(i,1)*Q_r%data(i,iqxvel)


       ! pressure
       Q_l%data(i+1,iqpres) = 0.0_dp_t
       Q_r%data(i,iqpres) = 0.0_dp_t

       do n = 1, nprim
          if (eval(n) >= 0.0_dp_t) then
             Q_l%data(i+1,iqpres) = Q_l%data(i+1,iqpres) + &
                  beta_xp(n)*rvec(n,3)
          endif

          if (eval(n) <= 0.0_dp_t) then
             Q_r%data(i,iqpres) = Q_r%data(i,iqpres) + &
                  beta_xm(n)*rvec(n,3)
          endif
       enddo

       Q_l%data(i+1,iqpres) = Q%data(i,iqpres) - &
            xi%data(i,1)*Q_l%data(i+1,iqpres)

       Q_r%data(i,iqpres) = Q%data(i,iqpres) - &
            xi%data(i,1)*Q_r%data(i,iqpres)

       
    enddo

    ! clean-up
    call destroy(Qminus)
    call destroy(Qplus)
    call destroy(Q6)

    call destroy(xi)


    !-------------------------------------------------------------------------
    ! apply the source terms
    !-------------------------------------------------------------------------
    do i = U%grid%lo, U%grid%hi+1
       Q_l%data(i,iqxvel) = Q_l%data(i,iqxvel) + 0.5_dp_t*dt*grav
       Q_r%data(i,iqxvel) = Q_r%data(i,iqxvel) + 0.5_dp_t*dt*grav
    enddo

    ! special fixes at the boundary -- gravity must be reflected (correct the above too)
    if (U%grid%xlboundary == "reflect") then
       Q_l%data(U%grid%lo,iqxvel) = &
            Q_l%data(U%grid%lo,iqxvel) - dt*grav
    endif

    if (U%grid%xrboundary == "reflect") then
       Q_r%data(U%grid%hi+1,iqxvel) = &
            Q_r%data(U%grid%hi+1,iqxvel) + dt*grav
    endif


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

  end subroutine make_interface_states_ppm

end module interface_states_ppm_module
