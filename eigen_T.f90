module eigen_T_module

  use datatypes_module
  use variables_module

  implicit none

contains

  subroutine eigen_T(r, ux, T, p_r, p_T, cs, lvec, rvec, eval)

    real (kind=dp_t), intent(in) :: r, ux, T, cs, p_r, p_T
    real (kind=dp_t), intent(out) :: rvec(nwaves,nprim), lvec(nwaves,nprim), eval(nwaves)

    ! compute the eigenvectors for the rho, u, T system

    ! compute the eigenvalues
    eval(1) = ux - cs
    eval(2) = ux
    eval(3) = ux + cs


    ! compute the left eigenvectors
    lvec(1,:) = [ HALF*p_r/(cs*cs), -HALF*r/cs, HALF*p_T/(cs*cs) ]   ! u - c
    lvec(2,:) = [ r - r*p_r/(cs*cs),  ZERO,       -p_T/(cs*cs) ]   ! u
    lvec(3,:) = [ HALF*p_r/(cs*cs), HALF*r/cs,  HALF*p_T/(cs*cs) ]   ! u + c

    ! compute the right eigenvectors
    rvec(1,:) = [ ONE, -cs/r, -p_r/p_T + cs*cs/p_T ]   ! u - c 
    rvec(2,:) = [ ONE, ZERO,  -p_r/p_T  ]   ! u  
    rvec(3,:) = [ ONE, cs/r,  -p_r/p_T + cs**2/p_T ]   ! u + c   

  end subroutine eigen_T
end module eigen_T_module
