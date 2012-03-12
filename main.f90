program hydro1d

  use datatypes_module
  use grid_module
  use params_module
  use probparams_module, only: init_probparams
  use init_module
  use variables_module
  use output_module
  use bcs_module
  use dt_module
  use interface_states_godunov_module
  use interface_states_plm_module
  use interface_states_ppm_module
  use riemann_module
  use update_module

  implicit none

  type(grid_t) :: grid
  type(gridvar_t) :: U
  type(gridedgevar_t) :: U_l, U_r, fluxes

  integer :: i, n, ng
  real (kind=dp_t) :: t, dt, dt_new

  ! hack
  real (kind=dp_t) :: p, p_below, e

  ! parse the inputs file
  call init_params()
  call init_probparams()

  ! set the number of ghostcells
  if (godunov_type == 0) then
     ng = 1
  else
     ng = 4
  endif


  ! build the grid and storage for grid variables, interface states,
  ! and fluxes
  call build(grid, nx, ng, xmin, xmax, xlboundary, xrboundary)

  call build(U, grid, ncons)
  call build(U_l, grid, ncons)
  call build(U_r, grid, ncons)

  call build(fluxes, grid, ncons)


  ! set the initial conditions
  n = 0
  t = 0.0_dp_t
  call init_data(U)

  ! set the boundary conditions
  call fillBCs(U)



  ! hack -- check HSE
  do i = U%grid%lo+1, U%grid%hi
     e = (U%data(i,iuener) - 0.5_dp_t*U%data(i,iumomx)**2/U%data(i,iudens))/U%data(i,iudens)
     call eos(eos_input_e, p, e, U%data(i,iudens))

     e = (U%data(i-1,iuener) - 0.5_dp_t*U%data(i-1,iumomx)**2/U%data(i-1,iudens))/U%data(i-1,iudens)
     call eos(eos_input_e, p_below, e, U%data(i-1,iudens))

     print *, i, ((p - p_below)/U%grid%dx - 0.5_dp_t*(U%data(i-1,iudens) + U%data(i,iudens))*grav)/ &
          (0.5_dp_t*(U%data(i-1,iudens) +U%data(i,iudens))*grav)
  enddo
!  stop

  ! output the initial model
  call output(U, t, n)


  ! main evolution loop
  do while (t < tmax)

     ! set the boundary conditions
     call fillBCs(U)

     ! compute the timestep
     call compute_dt(U, n, dt_new)
     if (n > 0) then
        dt = min(dt_change*dt, dt_new)
     else
        dt = dt_new
     endif

     if (t + dt > tmax) dt = tmax - t


     ! construct the interface states
     if (godunov_type == 0) then
        call make_interface_states_godunov(U, U_l, U_r, dt)

     else if (godunov_type == 1) then
        call make_interface_states_plm(U, U_l, U_r, dt)

     else if (godunov_type == 2) then
        call make_interface_states_ppm(U, U_l, U_r, dt)
     endif


     ! compute the fluxes
     call solve_riemann(U_l, U_r, fluxes)


     ! do the conservative update
     call update(U, fluxes, dt)

     t = t + dt
     n = n + 1


     ! output (if necessary)
     print *, n, t, dt

  enddo


  ! final output
  call output(U, t, n)


  ! clean-up
  call destroy(U)
  call destroy(grid)

end program hydro1d
