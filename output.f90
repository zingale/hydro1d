module output_module

  use grid_module
  use datatypes_module
  use params_module
  use variables_module
  use eos_module

  implicit none

contains

  subroutine output(U, t, n)

    type(gridvar_t),  intent(in) :: U
    real (kind=dp_t), intent(in) :: t
    integer,          intent(in) :: n

    character (len=4) :: filenum
    character (len=40) :: outfile
    integer :: i
    integer :: lun
    real (kind=dp_t) :: p, e

    ! construct the filename and open for output
    write(unit=filenum,fmt="(i4.4)") n
    outfile = trim(problem_name) // "_" // filenum

    open(newunit=lun, file=trim(outfile), status="unknown", action="write")

    
    ! write out the header
1   format("# problem: ", a)
2   format("# time = ", f8.5)
3   format("# n = ", i4)
4   format("# nx = ", i4)
5   format("#", 1x, a15, 1x, a15, 1x, a15, 1x, a15, 1x, a15, 1x, a15)

    write(unit=lun, fmt=1) trim(problem_name)
    write(unit=lun, fmt=2) t
    write(unit=lun, fmt=3) n
    write(unit=lun, fmt=4) U%grid%nx
    write(unit=lun, fmt=5) "x   ", "rho   ", "rho u   ", "rho E   ", "u   ", "p   ", "e   "

    ! write out the data
6   format(2x, 7(g16.10, 1x))

    do i = U%grid%lo, U%grid%hi

       e = (U%data(i,iuener) - &
            0.5_dp_t*U%data(i,iumomx)**2/U%data(i,iudens)) / &
            U%data(i,iudens)
       call eos(eos_input_e, p, e, U%data(i,iudens))

       write(unit=lun, fmt=6) U%grid%x(i), &
            U%data(i,iudens), U%data(i,iumomx), U%data(i,iuener), &
            U%data(i,iumomx)/U%data(i,iudens), p, e

    enddo

    close(unit=lun)

  end subroutine output
    

end module output_module
