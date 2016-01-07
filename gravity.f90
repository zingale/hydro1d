subroutine gravity(U, g)

  use datatypes_module
  use grid_module
  use variables_module
  use params_module, only : grav, gravity_monopole
  use constants_module

  implicit none

  type(gridedgevar_t), intent(in) :: U
  type(gridedgevar_t), intent(out) :: g

  type(grid_t) :: grid

  integer :: i
  real (kind=dp_t) :: M_enclosed
  
  grid = U%grid

  ! compute the gravitational acceleration on edges
  if (gravity_monopole == 0) then
     g%data(:,1) = grav
  endif

  M_enclosed = ZERO

  g%data(grid%lo,1) = ZERO  

  do i = grid%lo+1, grid%hi+1
     M_enclosed = M_enclosed + FOURTHIRDS*pi* &
          (grid%xr(i-1)**2 + grid%xl(i-1)**2 + grid%xl(i-1)*grid%xr(i-1)) * &
          grid%dx * U%data(i-1,iudens)

     g%data(i,1) = -G_newton*M_enclosed/grid%xl(i)**2
  enddo


  ! fill boundary conditions -- note that g is edge-centered here

end subroutine gravity
