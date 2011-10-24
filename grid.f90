! grid module
!
! The grid module defines datatypes for the grid coordinates (grid_t)
! and variables that live on the grid (gridvar_t).  It is assumed that
! every grid variable will have the same number of ghostcells (ng).
!
! The limits of the valid cell-centered data on the grid are given by
! lo, hi.  For edge-centered data, the limits are lo, hi+1.
!
! For the cell-centered grid variables, the overall grid looks like this:
!
!    
!    |     |      |     X     |     |      |     |     X     |      |     |
!    +--*--+- // -+--*--X--*--+--*--+- // -+--*--+--*--X--*--+- // -+--*--+
!      -ng          -1     0     1           ...  nx-1    nx        nx+ng-1
!                        (lo)                     (hi)    
!    |                  |                              |                  |
!    |<- ng ghostcells->|<---- nx interior zones ----->|<- ng ghostcells->|
!    |                  |                              |                  |
!
! Here we adopt the convention that the first valid zone is 0.1
!
! For the edge-centered variables, the overall grid looks like:
!
!    
!    |     |      |     X     |     |      |     |     X     |      |     |
!    +-----+- // -+-----X-----+-----+- // -+-----+-----X-----+- // -+-----+
!   -ng          -1     0     1           ...  nx-1    nx        nx+ng-1 nx+ng
!                      (lo)                          (hi+1)    
!    |                  |                              |                  |
!    |<- ng ghostcells->|<---- nx interior zones ----->|<- ng ghostcells->|
!    |                  |                              |                  |
!
!
!
! The different datatypes that live on the grid are build from the
! grid type itself.  It is assumed that all quantities on the grid
! will have the same number of ghostcells, for simplicity.

module grid_module

  use datatypes_module

  implicit none

  ! the datatype for the grid coordinate information
  type grid_t 
     integer :: lo = -1
     integer :: hi = -1
     integer :: nx = -1
     integer :: ng = -1
     real (kind=dp_t) :: xmin, xmax, dx
     real (kind=dp_t), allocatable :: x(:), xl(:), xr(:)
  end type grid_t

  ! the datatype for a variable that lives in a zone on the grid 
  type gridvar_t
     integer :: nvar
     real (kind=dp_t), allocatable :: data(:,:)
     type(grid_t) :: grid
  end type gridvar_t

  ! the datatype for a edge-centered variable
  type gridedgevar_t
     integer :: nedgevar
     real (kind=dp_t), allocatable :: data(:,:)
     type(grid_t) :: grid
  end type gridedgevar_t


  interface build
     module procedure build_grid
     module procedure build_gridvar
     module procedure build_gridedgevar
  end interface build

  interface destroy
     module procedure destroy_grid
     module procedure destroy_gridvar
     module procedure destroy_gridedgevar
  end interface destroy
  

contains  

  subroutine build_grid(grid,nx,ng,xmin,xmax)
    
    type(grid_t),     intent(inout) :: grid
    integer,          intent(in   ) :: nx
    integer,          intent(in   ) :: ng
    real (kind=dp_t), intent(in   ) :: xmin, xmax

    integer :: i

    grid%lo = 0
    grid%hi = nx-1

    grid%nx = nx
    grid%ng = ng

    grid%xmin = xmin
    grid%xmax = xmax
    grid%dx = (xmax - xmin)/nx

    allocate(grid%x(-ng:nx+ng-1))
    allocate(grid%xl(-ng:nx+ng-1))
    allocate(grid%xr(-ng:nx+ng-1))

    do i = grid%lo-ng, grid%hi+ng
       grid%xl(i) = dble(i  )*grid%dx + xmin
       grid%xr(i) = dble(i+1)*grid%dx + xmin

       grid%x(i) = 0.5_dp_t*(grid%xl(i) + grid%xr(i))
    enddo

  end subroutine build_grid


  subroutine destroy_grid(grid)

    type(grid_t), intent(inout) :: grid

    deallocate(grid%xl)
    deallocate(grid%xr)
    deallocate(grid%x)
    
  end subroutine destroy_grid


  subroutine build_gridvar(gridvar, grid, nvar)

    type(gridvar_t), intent(inout) :: gridvar
    type(grid_t),    intent(in   ) :: grid
    integer,         intent(in   ) :: nvar

    if (grid%nx == -1) then
       print *, "ERROR: grid not initialized"
    endif

    ! first build gridvar's grid type and then initialize it
    call build(gridvar%grid, grid%nx, grid%ng, grid%xmin, grid%xmax)
    gridvar%grid = grid

    ! now initialize the storage for the grid data
    allocate(gridvar%data(-grid%ng:grid%nx+grid%ng-1,nvar))
    gridvar%data(:,:) = 0.0_dp_t

  end subroutine build_gridvar


  subroutine destroy_gridvar(gridvar)

    type(gridvar_t), intent(inout) :: gridvar

    call destroy(gridvar%grid)
    deallocate(gridvar%data)

  end subroutine destroy_gridvar


  subroutine build_gridedgevar(gridedgevar, grid, nvar)

    type(gridedgevar_t), intent(inout) :: gridedgevar
    type(grid_t),        intent(in   ) :: grid
    integer,             intent(in   ) :: nvar

    if (grid%nx == -1) then
       print *, "ERROR: grid not initialized"
    endif

    ! first build gridvar's grid type and then initialize it
    call build(gridedgevar%grid, grid%nx, grid%ng, grid%xmin, grid%xmax)
    gridedgevar%grid = grid

    ! now initialize the storage for the grid data
    allocate(gridedgevar%data(-grid%ng:grid%nx+grid%ng,nvar))
    gridedgevar%data(:,:) = 0.0_dp_t

  end subroutine build_gridedgevar


  subroutine destroy_gridedgevar(gridedgevar)

    type(gridedgevar_t), intent(inout) :: gridedgevar

    call destroy(gridedgevar%grid)
    deallocate(gridedgevar%data)
    
  end subroutine destroy_gridedgevar
   
end module grid_module
