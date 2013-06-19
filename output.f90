module output_module

  use grid_module
  use datatypes_module
  use params_module
  use variables_module
  use eos_module

  implicit none

  private

  public :: output

contains

  subroutine output(U, t, n)

    type(gridvar_t),  intent(in) :: U
    real (kind=dp_t), intent(in) :: t
    integer,          intent(in) :: n

    character (len=6) :: filenum
    character (len=40) :: outfile
    integer :: i
    integer :: lun
    real (kind=dp_t) :: p, e, rho

    ! construct the filename and open for output
    write(unit=filenum,fmt="(i6.6)") n
    if (ppm_temp) then
       outfile = trim(problem_name) // "_ppmT_" // filenum
    else
       outfile = trim(problem_name) // "_" // filenum
    endif

    open(newunit=lun, file=trim(outfile), status="unknown", action="write")

    
    ! write out the header
1   format("# problem: ", a)
2   format("# time = ", f8.5)
3   format("# n = ", i4)
4   format("# nx = ", i4)
5   format("# ", 9(a16,1x))

    write(unit=lun, fmt=1) trim(problem_name)
    write(unit=lun, fmt=2) t
    write(unit=lun, fmt=3) n
    write(unit=lun, fmt=4) U%grid%nx
    write(unit=lun, fmt=5) "x   ", "rho   ", "rho u   ", "rho E   ", "u   ", "p   ", "e   ", "cs  ", "M   "

    ! write out the data
6   format(2x, 9(g16.10, 1x))

    if (write_ghost) then

       do i = U%grid%lo-U%grid%ng, U%grid%hi+U%grid%ng

          e = (U%data(i,iuener) - &
               0.5_dp_t*U%data(i,iumomx)**2/U%data(i,iudens)) / &
               U%data(i,iudens)
          rho = U%data(i,iudens)
          call eos(eos_input_e, p, e, rho)
          
          write(unit=lun, fmt=6) U%grid%x(i), &
               U%data(i,iudens), U%data(i,iumomx), U%data(i,iuener), &
               U%data(i,iumomx)/U%data(i,iudens), p, e, sqrt(gamma*p/rho), &
               abs(U%data(i,iumomx)/U%data(i,iudens))/sqrt(gamma*p/rho)
          
       enddo

    else

       do i = U%grid%lo, U%grid%hi

          e = (U%data(i,iuener) - &
               0.5_dp_t*U%data(i,iumomx)**2/U%data(i,iudens)) / &
               U%data(i,iudens)
          rho = U%data(i,iudens)
          call eos(eos_input_e, p, e, rho)
          
          write(unit=lun, fmt=6) U%grid%x(i), &
               U%data(i,iudens), U%data(i,iumomx), U%data(i,iuener), &
               U%data(i,iumomx)/U%data(i,iudens), p, e, sqrt(gamma*p/rho)
          
       enddo

    endif

    close(unit=lun)

  end subroutine output
    

end module output_module
