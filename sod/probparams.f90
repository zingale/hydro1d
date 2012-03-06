!=============================================================================
! probparams module
!
! The propparams module controls the runtime parameters for the
! problem-specific initial conditions only.  It reads them in from
! the same inputs file as the params module, but looks for a different
! namelist.
!=============================================================================

module probparams_module

  use datatypes_module
  implicit none

  !---------------------------------------------------------------------------
  ! declare runtime parameters with SAVE attribute and default values
  !---------------------------------------------------------------------------

  real (kind=dp_t), save :: rho_l = 1.0_dp_t
  real (kind=dp_t), save :: u_l   = 0.0_dp_t
  real (kind=dp_t), save :: p_l   = 1.0_dp_t

  real (kind=dp_t), save :: rho_r = 0.125_dp_t
  real (kind=dp_t), save :: u_r   = 0.0_dp_t
  real (kind=dp_t), save :: p_r   = 0.1_dp_t


  ! namelist
  namelist /problem/ rho_l, u_l, p_l, rho_r, u_r, p_r

contains

  !===========================================================================
  ! init_params
  !===========================================================================
  subroutine init_probparams

    use params_module, only: infile

    integer :: lun, status

    if (.not. infile == "") then

       print *, 'reading problem-specific parameters from ', trim(infile)
       
       open(newunit=lun, file=trim(infile), status="old", action="read")
       read(unit=lun, nml=problem, iostat=status)

       if (status < 0) then
          ! namelist doesn't exist.  Use the defaults and carry on.
          continue

       else if (status > 0) then
          ! some problem with the namelist.
          print *, "ERROR: problem namelist invalid"
          stop

       endif

       close(unit=lun)

    endif    

  end subroutine init_probparams

end module probparams_module
