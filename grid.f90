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
! Here we adopt the convention that the first valid zone is 0
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

  real (kind=dp_t) :: pi = 3.141592653589793238462643383279502884_dp_t

  ! the datatype for the grid coordinate information
  type grid_t 
     integer :: lo = -1
     integer :: hi = -1
     integer :: nx = -1
     integer :: ng = -1

     real (kind=dp_t) :: xmin, xmax, dx

     ! geometry: 0 = Cartesian; 1 = spherical
     integer :: geometry = 0

     real (kind=dp_t), pointer :: x(:) => Null()
     real (kind=dp_t), pointer :: xl(:) => Null()
     real (kind=dp_t), pointer :: xr(:) => Null()  

     real (kind=dp_t), pointer :: Al(:) => Null()
     real (kind=dp_t), pointer :: Ar(:) => Null()

     real (kind=dp_t), pointer :: V(:) => Null()

     character (len=32) :: xlboundary
     character (len=32) :: xrboundary

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
  
  private

  public :: grid_t, gridvar_t, gridedgevar_t
  public :: build, destroy


contains  

  subroutine build_grid(grid,nx,ng,xmin,xmax,xlbc,xrbc,is_spherical)
    
    type(grid_t),       intent(inout) :: grid
    integer,            intent(in   ) :: nx
    integer,            intent(in   ) :: ng
    real (kind=dp_t),   intent(in   ) :: xmin, xmax
    character (len=32), intent(in   ) :: xlbc, xrbc
    integer :: i
    logical, optional,  intent(in   ) :: is_spherical

    if (present(is_spherical)) then
       if (is_spherical) then
          grid%geometry = 1
       else
          grid%geometry = 0
       endif
    endif

    grid%lo = 0
    grid%hi = nx-1

    grid%nx = nx
    grid%ng = ng

    grid%xmin = xmin
    grid%xmax = xmax
    grid%dx = (xmax - xmin)/nx

    grid%xlboundary = xlbc
    grid%xrboundary = xrbc

    allocate(grid%x(-ng:nx+ng-1))
    allocate(grid%xl(-ng:nx+ng-1))
    allocate(grid%xr(-ng:nx+ng-1))

    allocate(grid%Al(-ng:nx+ng-1))
    allocate(grid%Ar(-ng:nx+ng-1))
    allocate(grid%V(-ng:nx+ng-1))

    do i = grid%lo-ng, grid%hi+ng
       grid%xl(i) = dble(i  )*grid%dx + xmin
       grid%xr(i) = dble(i+1)*grid%dx + xmin

       grid%x(i) = 0.5_dp_t*(grid%xl(i) + grid%xr(i))
    enddo

    if (grid%geometry == 1) then
       do i = grid%lo-ng, grid%hi+ng
          grid%Al(i) = grid%xl(i)**2
          grid%Ar(i) = grid%xr(i)**2
          
          grid%dV(i) = grid%x(i)**2 * grid%dx
       enddo
    else
       grid%Al(:) = 1.0_dp_t
       grid%Ar(:) = 1.0_dp_t
       grid%dV(:) = grid%dx
    endif

  end subroutine build_grid


  subroutine destroy_grid(grid)

    type(grid_t), intent(inout) :: grid

    deallocate(grid%xl)
    deallocate(grid%xr)
    deallocate(grid%x)

    deallocate(grid%Al)
    deallocate(grid%Ar)
    deallocate(grid%V)

  end subroutine destroy_grid


  subroutine build_gridvar(gridvar, grid, nvar)

    type(gridvar_t), intent(inout) :: gridvar
    type(grid_t),    intent(in   ) :: grid
    integer,         intent(in   ) :: nvar

    if (grid%nx == -1) then
       print *, "ERROR: grid not initialized"
    endif

    ! gridvar's grid type is simply a copy of the input grid
    gridvar%grid = grid

    ! now initialize the storage for the grid data
    allocate(gridvar%data(-grid%ng:grid%nx+grid%ng-1,nvar))
    gridvar%data(:,:) = 0.0_dp_t

    gridvar%nvar = nvar

  end subroutine build_gridvar


  subroutine destroy_gridvar(gridvar)

    type(gridvar_t), intent(inout) :: gridvar

    deallocate(gridvar%data)

  end subroutine destroy_gridvar


  subroutine build_gridedgevar(gridedgevar, grid, nvar)

    type(gridedgevar_t), intent(inout) :: gridedgevar
    type(grid_t),        intent(in   ) :: grid
    integer,             intent(in   ) :: nvar

    if (grid%nx == -1) then
       print *, "ERROR: grid not initialized"
    endif

    ! gridedgevar's grid type is simply a copy of the input grid
    gridedgevar%grid = grid

    ! now initialize the storage for the grid data
    allocate(gridedgevar%data(-grid%ng:grid%nx+grid%ng,nvar))
    gridedgevar%data(:,:) = 0.0_dp_t

    gridedgevar%nedgevar = nvar

  end subroutine build_gridedgevar


  subroutine destroy_gridedgevar(gridedgevar)

    type(gridedgevar_t), intent(inout) :: gridedgevar

    deallocate(gridedgevar%data)
    
  end subroutine destroy_gridedgevar
   
end module grid_module
