! runtime diagonstics

module runtime_diag_module

  use grid_module
  use datatypes_module
  use params_module
  use variables_module
  use eos_module

  implicit none

  private

  public :: run_diag

contains

  subroutine run_diag(U, t, n)

    type(gridvar_t),  intent(in) :: U
    real (kind=dp_t), intent(in) :: t
    integer, intent(in) :: n

    real (kind=dp_t) :: mach_max

    real (kind=dp_t) :: e, rho, p

    integer :: i, lun

    character (len=256) :: outfile

    mach_max = -1.d33

    do i = U%grid%lo, U%grid%hi

       e = (U%data(i,iuener) - &
            0.5_dp_t*U%data(i,iumomx)**2/U%data(i,iudens)) / &
            U%data(i,iudens)
       rho = U%data(i,iudens)
       call eos(eos_input_e, p, e, rho)

       mach_max = max(mach_max, &
                      abs(U%data(i,iumomx)/U%data(i,iudens))/sqrt(gamma*p/rho))
    enddo
    
    outfile = trim(problem_name) // "_diag.out"

    if (n == 1) then
       open(newunit=lun, file=trim(outfile), status="replace")       
    else
       open(newunit=lun, file=trim(outfile), status="old", position="append")
    endif

100 format(1x, g20.10, 1x, g20.10)
    write(unit=lun,fmt=100) t, mach_max

    close (unit=lun)

  end subroutine run_diag
end module runtime_diag_module

