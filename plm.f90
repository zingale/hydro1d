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

    real (kind=dp_t) :: e
    real (kind=dp_t) :: test

    integer :: i, n

    ! piecewise linear slopes -- this is a 1-d version of the Colella
    ! 2nd order unsplit Godunov scheme (Colella 1990).  See also
    ! Colella & Glaz.
    !
    ! we wish to solve
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




    !-------------------------------------------------------------------------
    ! compute the monotonized central differences
    !-------------------------------------------------------------------------

    ! 4th order MC limiting.  See Colella (1985) Eq. 2.5 and 2.6 or 
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


    !-------------------------------------------------------------------------
    ! compute left and right primitive variable states
    !-------------------------------------------------------------------------
    call build(Q_l, U%grid, nprim)
    call build(Q_r, U%grid, nprim)

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
    call destroy(Q_l)
    call destroy(Q_r)


    !-------------------------------------------------------------------------
    ! apply the source terms
    !-------------------------------------------------------------------------




  end subroutine make_interface_states_plm

end module interface_states_plm_module
