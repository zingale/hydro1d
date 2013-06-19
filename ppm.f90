module interface_states_ppm_module

  use datatypes_module
  use grid_module
  use variables_module
  use params_module
  use eos_module
  use eigen_module
  use eigen_T_module
  use flatten_module

  implicit none

  private

  public :: make_interface_states_ppm

contains
  
  subroutine make_interface_states_ppm(U, U_l, U_r, dt)

    type(gridvar_t),     intent(in   ) :: U
    type(gridedgevar_t), intent(inout) :: U_l, U_r
    real (kind=dp_t),    intent(in   ) :: dt

    type(gridvar_t) :: Q, Qhat
    type(gridedgevar_t) :: Q_l, Q_r

    real (kind=dp_t) :: rvec(nwaves,nprim), lvec(nwaves,nprim), eval(nwaves)

    real (kind=dp_t) :: r, ux, p, cs, T, p_r, p_T
    real (kind=dp_t) :: beta_xm(nwaves), beta_xp(nwaves)

    real (kind=dp_t) :: dq0, dqp
    real (kind=dp_t) :: Iminus(nwaves,nprim), Iplus(nwaves,nprim)
    real (kind=dp_t) :: Qref_xm(nprim), Qref_xp(nprim)

    type(gridedgevar_t) :: Qminus, Qplus
    type(gridvar_t) :: Q6

    real (kind=dp_t) :: sigma

    type(gridvar_t) :: xi

    real (kind=dp_t) :: dtdx

    real (kind=dp_t) :: e, T_l, T_r

    integer :: i, m, n

    ! piecewise parabolic slopes 
    !
    ! This is a 1-d version of the piecewise parabolic method detailed
    ! Colella & Woodward (1984).  We follow the description of Almgren
    ! et al. 2010 (the CASTRO paper) and Miller and Colella (2002).  
    !
    ! Also note that we do not implement the contact steepening here.
    
    ! sanity check
    if (U%grid%ng < 4) then
       print *, "ERROR: ng < 4 in plm states"
       stop
    endif


    !-------------------------------------------------------------------------
    ! convert to primitve variables
    !-------------------------------------------------------------------------
    call build(Q, U%grid, nprim)
    call build(Qhat, U%grid, nprim)

    do i = U%grid%lo-U%grid%ng, U%grid%hi+U%grid%ng
       ! density
       Q%data(i,iqdens) = U%data(i,iudens)

       ! velocity
       Q%data(i,iqxvel) = U%data(i,iumomx)/U%data(i,iudens)

       Qhat%data(i,iqdens) = Q%data(i,iqdens)
       Qhat%data(i,iqxvel) = Q%data(i,iqxvel)

       ! pressure
       e = (U%data(i,iuener) - &
            HALF*U%data(i,iumomx)**2/U%data(i,iudens))/U%data(i,iudens)
       call eos(eos_input_e, Q%data(i,iqpres), e, Q%data(i,iqdens), T = Qhat%data(i,iqtemp) )
    enddo


    !-------------------------------------------------------------------------
    ! compute the flattening coefficients
    !-------------------------------------------------------------------------
    call build(xi, U%grid, 1)
    call flatten(Q, xi)


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
          dq0 = HALF*(Qhat%data(i+1,n) - Qhat%data(i-1,n))
          dqp = HALF*(Qhat%data(i+2,n) - Qhat%data(i,n))

          ! limiting (C&W Eq. 1.8)
          if ( (Qhat%data(i+1,n) - Qhat%data(i,n))* &
               (Qhat%data(i,n) - Qhat%data(i-1,n)) > ZERO) then
             dq0 = sign(ONE,dq0)* &
                  min( abs(dq0), &
                       2.0_dp_t*abs(Qhat%data(i,n) - Qhat%data(i-1,n)), &
                       2.0_dp_t*abs(Qhat%data(i+1,n) - Qhat%data(i,n)) )
          else
             dq0 = ZERO
          endif

          if ( (Qhat%data(i+2,n) - Qhat%data(i+1,n))* &
               (Qhat%data(i+1,n) - Qhat%data(i,n)) > ZERO) then
             dqp = sign(ONE,dqp)* &
                  min( abs(dqp), &
                       2.0_dp_t*abs(Qhat%data(i+1,n) - Qhat%data(i,n)), &
                       2.0_dp_t*abs(Qhat%data(i+2,n) - Qhat%data(i+1,n)) )
          else
             dqp = ZERO
          endif

          ! cubic (C&W Eq. 1.6)
          Qplus%data(i,n) = HALF*(Qhat%data(i,n) + Qhat%data(i+1,n)) - &
               (ONE/6.0_dp_t)*(dqp - dq0)

          Qminus%data(i+1,n) = Qplus%data(i,n)

          ! make sure that we didn't over or undersoot -- this may not
          ! be needed, but is discussed in Colella & Sekora (2008)
          Qplus%data(i,n) = max(Qplus%data(i,n), min(Qhat%data(i,n),Qhat%data(i+1,n)))
          Qplus%data(i,n) = min(Qplus%data(i,n), max(Qhat%data(i,n),Qhat%data(i+1,n)))

          Qminus%data(i+1,n) = max(Qminus%data(i+1,n), min(Qhat%data(i,n),Qhat%data(i+1,n)))
          Qminus%data(i+1,n) = min(Qminus%data(i+1,n), max(Qhat%data(i,n),Qhat%data(i+1,n)))

       enddo
    enddo


    !-------------------------------------------------------------------------
    ! construct the parameters for the parabolic reconstruction polynomials
    !-------------------------------------------------------------------------

    ! Our parabolic profile has the form:
    !
    !  q(xi) = qminus + xi*(qplus - qminus + q6 * (1-xi) )
    !
    ! with xi = (x - xl)/dx, where xl is the interface of the left
    ! edge of the cell.  qminus and qplus are the values of the 
    ! parabola on the left and right edges of the current cell.

    ! Limit the left and right values of the parabolic interpolant
    ! (C&W Eq. 1.10).  Here the loop is over cells, and considers the
    ! values on either side of the center of the cell (Qminus and
    ! Qplus).
    do n = 1, nprim
       do i = U%grid%lo-1, U%grid%hi+1

          if ( (Qplus%data(i,n) - Qhat%data(i,n)) * &
               (Qhat%data(i,n) - Qminus%data(i,n)) <= ZERO) then
             Qminus%data(i,n) = Qhat%data(i,n)
             Qplus%data(i,n) = Qhat%data(i,n)

          else if ( (Qplus%data(i,n) - Qminus%data(i,n)) * &
                    (Qhat%data(i,n) - &
                      HALF*(Qminus%data(i,n) + Qplus%data(i,n))) > &
                   (Qplus%data(i,n) - Qminus%data(i,n))**2/6.0_dp_t ) then

          ! alternate test from Colella & Sekora (2008)
          !else if (abs(Qminus%data(i,n) - Qhat%data(i,n)) >= &
          !     2.0*abs(Qplus%data(i,n) - Qhat%data(i,n))) then
             Qminus%data(i,n) = 3.0_dp_t*Qhat%data(i,n) - 2.0_dp_t*Qplus%data(i,n)

          else if (-(Qplus%data(i,n) - Qminus%data(i,n))**2/6.0_dp_t > &
                    (Qplus%data(i,n) - Qminus%data(i,n)) * &
                    (Qhat%data(i,n) - &
                           HALF*(Qminus%data(i,n) + Qplus%data(i,n))) ) then

          !else if (abs(Qplus%data(i,n) - Qhat%data(i,n)) >= &
          !     2.0*abs(Qminus%data(i,n) - Qhat%data(i,n))) then
             Qplus%data(i,n) = 3.0_dp_t*Qhat%data(i,n) - 2.0_dp_t*Qminus%data(i,n)

          endif

       enddo
    enddo

    ! define Q6
    call build(Q6,  U%grid, nprim)

    do n = 1, nprim
       do i = U%grid%lo-1, U%grid%hi+1

          Q6%data(i,n) = 6.0_dp_t*Qhat%data(i,n) - &
               3.0_dp_t*(Qminus%data(i,n) + Qplus%data(i,n))

       enddo
    enddo


    !-------------------------------------------------------------------------
    ! compute left and right primitive variable states
    !-------------------------------------------------------------------------

    call build(Q_l, U%grid, nprim)
    call build(Q_r, U%grid, nprim)

    dtdx = dt/U%grid%dx

    ! The fluxes are going to be defined on the left edge of the
    ! computational zone.
    !
    !          |             |             |             | 
    !          |             |             |             |  
    !         -+------+------+------+------+------+------+--   
    !          |     i-1     |      i      |     i+1     |   
    !                       * *           * 
    !                   q_l,i  q_r,i   q_l,i+1 
    !           
    ! q_l,i+1 are computed using the information in zone i.
    !

    do i = U%grid%lo-1, U%grid%hi+1

       r  = Qhat%data(i,iqdens)
       ux = Qhat%data(i,iqxvel)
       T  = Qhat%data(i,iqtemp)

       ! for the EOS call
       p  = Q%datA(i,iqpres)

       call eos(eos_input_p, p, e, r, p_r=p_r, p_T=p_T, c=cs)

       ! get the eigenvalues and eigenvectors
       call eigen_T(r, ux, T, p_r, p_T, cs, lvec, rvec, eval)


       ! integrate the parabola in the cell from the left interface
       ! (Iminus) over the portion of the cell that each eigenvalue
       ! can reach.  Do the same from the right interface in the cell,
       ! defining Iplus.  See Almgren et al. 2010 (Eq. 30) or Colella
       ! & Sekora (2008), or Miller & Colella (2002), Eq. 90.
       do m = 1, nwaves
          sigma = abs(eval(m))*dtdx
          do n = 1, nprim

             ! only integrate if the wave is moving toward the interface
             ! (see Miller & Colella, Eg. 90).  This may not be necessary.
             if (eval(m) >= ZERO) then
                Iplus(m,n) = Qplus%data(i,n) - HALF*sigma* &
                     (Qplus%data(i,n) - Qminus%data(i,n) - &
                     (ONE - (2.0_dp_t/3.0_dp_t)*sigma)*Q6%data(i,n))
             else
                Iplus(m,n) = Qhat%data(i,n)
             endif

             if (eval(m) <= ZERO) then
                Iminus(m,n) = Qminus%data(i,n) + HALF*sigma* &
                     (Qplus%data(i,n) - Qminus%data(i,n) + &
                     (ONE - (2.0_dp_t/3.0_dp_t)*sigma)*Q6%data(i,n))
             else
                Iminus(m,n) = Qhat%data(i,n)
             endif

          enddo
       enddo


       ! the basic idea here is that we do a characteristic
       ! decomposition.  The jump in primitive variables (Q) can be
       ! transformed to a jump in characteristic variables using the
       ! left and right eigenvectors.  Then each wave tells us how
       ! much of each characteristic quantity reaches the interface
       ! over dt/2.  We only add the quantity if it moves toward the
       ! interface.
       !
       ! See Miller & Colella for a good discussion of the
       ! characteristic form.  The basic form is:
       !
       !
       !  n+1/2      ~               ~
       ! q         = q  -  sum l . ( q   - I  ) r
       !  i+1/2,L     L     i   i     L     +    i
       ! 
       !       ~
       ! Where q is the reference state.

       
       ! define the reference states -- Miller & Colella (2002) argue
       ! picking the fastest wave speed.  We follow the convention from
       ! Colella & Glaz, and Colella (1990) -- this is intended to 
       ! minimize the size of the jump that the projection operates on.
       if (.false.) then
          ! CASTRO method
          Qref_xm(:) = Qhat%data(i,:)
          Qref_xp(:) = Qhat%data(i,:)
       else
          ! Miller and Colella method
          if (eval(3) >= ZERO) then
             Qref_xp(:) = Iplus(3,:)
          else
             Qref_xp(:) = Qhat%data(i,:)
          endif

          if (eval(1) <= ZERO) then
             Qref_xm(:) = Iminus(1,:)
          else
             Qref_xm(:) = Qhat%data(i,:)
          endif             
       endif


       ! compute the dot product of each left eigenvector with (qref - I)
       do m = 1, nwaves    ! loop over waves
          beta_xm(m) = ZERO
          beta_xp(m) = ZERO

          beta_xm(m) = beta_xm(m) + dot_product(lvec(m,:),Qref_xm(:) - Iminus(m,:))
          beta_xp(m) = beta_xp(m) + dot_product(lvec(m,:),Qref_xp(:) - Iplus(m,:))
       enddo

       ! finally, sum up all the jumps

       ! density
       Q_l%data(i+1,iqdens) = ZERO
       Q_r%data(i,iqdens) = ZERO

       do n = 1, nwaves
          if (eval(n) >= ZERO) then
             Q_l%data(i+1,iqdens) = Q_l%data(i+1,iqdens) + &
                  beta_xp(n)*rvec(n,iqdens)
          endif

          if (eval(n) <= ZERO) then
             Q_r%data(i,iqdens) = Q_r%data(i,iqdens) + &
                  beta_xm(n)*rvec(n,iqdens)
          endif
       enddo

       Q_l%data(i+1,iqdens) = Qref_xp(iqdens) - Q_l%data(i+1,iqdens)
       Q_r%data(i,iqdens)   = Qref_xm(iqdens) - Q_r%data(i,iqdens)


       ! velocity
       Q_l%data(i+1,iqxvel) = ZERO
       Q_r%data(i,iqxvel) = ZERO

       do n = 1, nwaves
          if (eval(n) >= ZERO) then
             Q_l%data(i+1,iqxvel) = Q_l%data(i+1,iqxvel) + &
                  beta_xp(n)*rvec(n,iqxvel)
          endif

          if (eval(n) <= ZERO) then
             Q_r%data(i,iqxvel) = Q_r%data(i,iqxvel) + &
                  beta_xm(n)*rvec(n,iqxvel)
          endif
       enddo

       Q_l%data(i+1,iqxvel) = Qref_xp(iqxvel) - Q_l%data(i+1,iqxvel)
       Q_r%data(i,iqxvel)   = Qref_xm(iqxvel) - Q_r%data(i,iqxvel)


       ! pressure
       Q_l%data(i+1,iqtemp) = ZERO
       Q_r%data(i,iqtemp) = ZERO

       do n = 1, nwaves
          if (eval(n) >= ZERO) then
             Q_l%data(i+1,iqtemp) = Q_l%data(i+1,iqtemp) + &
                  beta_xp(n)*rvec(n,iqtemp)
          endif

          if (eval(n) <= ZERO) then
             Q_r%data(i,iqtemp) = Q_r%data(i,iqtemp) + &
                  beta_xm(n)*rvec(n,iqtemp)
          endif
       enddo

       Q_l%data(i+1,iqtemp) = Qref_xp(iqtemp) - Q_l%data(i+1,iqtemp)
       Q_r%data(i,iqtemp)   = Qref_xm(iqtemp) - Q_r%data(i,iqtemp)
       
       ! flatten
       Q_l%data(i+1,:) = (ONE - xi%data(i,1))*Qhat%data(i,:) + xi%data(i,1)*Q_l%data(i+1,:)
       Q_r%data(i,:)   = (ONE - xi%data(i,1))*Qhat%data(i,:) + xi%data(i,1)*Q_r%data(i,:)

       ! convert the temperature edge state into a pressure edge state via the EOS
       T_l = Q_l%data(i+1,iqtemp)
       T_r = Q_r%data(i,iqtemp)

       call eos(eos_input_T, Q_l%data(i+1,iqpres), e, Q_l%data(i+1,iqdens), T = T_l)
       call eos(eos_input_T, Q_r%data(i  ,iqpres), e, Q_r%data(i  ,iqdens), T = T_r)

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
       Q_l%data(i,iqxvel) = Q_l%data(i,iqxvel) + HALF*dt*grav
       Q_r%data(i,iqxvel) = Q_r%data(i,iqxvel) + HALF*dt*grav
    enddo

    ! special fixes at the boundary -- gravity must be reflected
    ! (correct the above too)
    if (U%grid%xlboundary == "reflect") then
       Q_l%data(U%grid%lo,iqxvel) = &
            Q_l%data(U%grid%lo,iqxvel) - dt*grav
    endif

    if (U%grid%xrboundary == "reflect") then
       Q_r%data(U%grid%hi+1,iqxvel) = &
            Q_r%data(U%grid%hi+1,iqxvel) - dt*grav
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
            HALF*Q_l%data(i,iqdens)*Q_l%data(i,iqxvel)**2

       call eos(eos_input_p, Q_r%data(i,iqpres), e, Q_r%data(i,iqdens))
       U_r%data(i,iuener) = Q_r%data(i,iqdens)*e + &
            HALF*Q_r%data(i,iqdens)*Q_r%data(i,iqxvel)**2

    enddo
    

    ! clean-up
    call destroy(Q)
    call destroy(Q_l)
    call destroy(Q_r)

  end subroutine make_interface_states_ppm

end module interface_states_ppm_module
