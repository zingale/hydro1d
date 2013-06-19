!=============================================================================
! params module
!
! The params module is responsible for reading in runtime parameters
! and providing them to any routines that need them.
!=============================================================================

module params_module

  use datatypes_module
  implicit none

  !---------------------------------------------------------------------------
  ! declare runtime parameters with SAVE attribute and default values
  !---------------------------------------------------------------------------

  ! grid parameters
  integer, save :: nx = 16
  real (kind=dp_t), save :: xmin = 0.0_dp_t
  real (kind=dp_t), save :: xmax = 1.0_dp_t
  character (len=32), save :: xlboundary = "outflow"
  character (len=32), save :: xrboundary = "outflow"
  character (len=32), save :: hse_bc_const = "temperature"
  character (len=32), save :: hse_vel_type = "reflect"

  ! timestep parameters
  real (kind=dp_t), save :: init_shrink = 1.0_dp_t
  real (kind=dp_t), save :: dt_change = 1.2_dp_t
  real (kind=dp_t), save :: tmax = 1.0_dp_t
  real (kind=dp_t), save :: cfl = 0.8_dp_t

  ! eos parameters
  real (kind=dp_t), save :: gamma = 1.4_dp_t

  ! hydro parameters
  logical, save :: ppm_temp = .true.
  real (kind=dp_t), save :: grav = 0.0_dp_t
  
  ! problem-specific parameters
  character (len=32), save :: problem_name = "sod"

  ! output parameters
  logical, save :: write_ghost = .false.
  real (kind=dp_t), save :: plot_dt = -1.0_dp_t


  ! namelist
  namelist /params/ nx, xmin, xmax, &
                    xlboundary, xrboundary, &
                    hse_bc_const, hse_vel_type, &
                    init_shrink, dt_change, tmax, &
                    cfl, &
                    gamma, &
                    ppm_temp, &
                    grav, &
                    problem_name, &
                    write_ghost, &
                    plot_dt

  character (len=32), save :: infile = ""

contains

  !===========================================================================
  ! init_params
  !===========================================================================
  subroutine init_params()


    integer :: lun

    ! read in the inputs file, if any -- we assume that the only command
    ! line argument is the inputs file name
    if (command_argument_count() == 0) then
       print *, 'no inputs file specified -- using parameter defaults'

    else if (command_argument_count() > 1) then
       print *, 'ERROR: only one command line argument allowed'

    else
       call get_command_argument(1, infile)
       print *, 'reading parameters from ', trim(infile)
       
       open(newunit=lun, file=trim(infile), status="old", action="read")
       read(unit=lun, nml=params)
       close(unit=lun)

    endif
    

  end subroutine init_params

end module params_module
