module eigen_module

  use datatypes_module
  use variables_module

  implicit none

contains

  subroutine eigen(r, ux, p, cs, lvec, rvec, eval)

    real (kind=dp_t), intent(in) :: r, ux, p, cs
    real (kind=dp_t), intent(out) :: rvec(nwaves,nprim), lvec(nwaves,nprim), eval(nwaves)


    ! Compute the left and right eigenvectors and the eigenvalues for
    ! the system
    !
    ! The eigenvalues are u - c, u, u + c
    !
    ! The Jacobian matrix for the primitive variable formulation of the
    ! Euler equations is 
    !
    !       / u   r   0   \  
    !   A = | 0   u   1/r |  
    !       \ 0  rc^2 u   /
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

    ! compute the eigenvalues
    eval(1) = ux - cs
    eval(2) = ux
    eval(3) = ux + cs


    ! compute the left eigenvectors
    lvec(1,:) = [ ZERO, -HALF*r/cs, HALF/(cs*cs) ]   ! u - c
    lvec(2,:) = [ ONE,  ZERO,       -ONE/(cs*cs) ]   ! u
    lvec(3,:) = [ ZERO, HALF*r/cs,  HALF/(cs*cs) ]   ! u + c

    ! compute the right eigenvectors
    rvec(1,:) = [ ONE, -cs/r, cs*cs ]   ! u - c 
    rvec(2,:) = [ ONE, ZERO,  ZERO  ]   ! u  
    rvec(3,:) = [ ONE, cs/r,  cs*cs ]   ! u + c   

  end subroutine eigen
end module eigen_module
