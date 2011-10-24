module init_module

  use datatypes_module
  use params_module
  use grid_module
  use eos_module
  use variables_module

  implicit none

contains
  
  subroutine init_data(U)

    type(gridvar_t), intent(inout) :: U

    integer :: i
    real (kind=dp_t) :: rho, ux, p, e

    real (kind=dp_t) :: xcenter

    select case (problem_name)

    case ("sod")

       xcenter = 0.5_dp_t*(U%grid%xmin + U%grid%xmax)

       do i = U%grid%lo, U%grid%hi

          if (U%grid%x(i) < xcenter) then

             ! left state
             rho = 1.0_dp_t
             ux  = 0.0_dp_t
             p   = 1.0_dp_t

          else

             ! right state
             rho = 0.125_dp_t
             ux  = 0.0_dp_t
             p   = 0.1_dp_t

          endif

          ! get the internal energy from the EOS
          call eos(eos_input_p, p, e, rho)

          ! store the conserved state
          U%data(i,iudens) = rho
          U%data(i,iumomx) = rho*ux
          U%data(i,iuener) = 0.5_dp_t*rho*ux*ux + rho*e
          
       enddo

    case default
       print *, "ERROR: problem name not found: ", problem_name
       stop

    end select

  end subroutine init_data

end module init_module

