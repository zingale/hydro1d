module interface_states_plm_module

  use datatypes_module
  use grid_module
  use variables_module
  use params_module

  implicit none

contains
  
  subroutine make_interface_states_plm(U, U_l, U_r, dt)

    type(gridvar_t),     intent(in   ) :: U
    type(gridedgevar_t), intent(inout) :: U_l, U_r
    real (kind=dp_t),    intent(in   ) :: dt

    integer :: i

    ! piecewise linear slopes -- this is a 1-d version of the Colella
    ! 2nd order unsplit Godunov scheme (Colella 1990).  See also
    ! Colella & Glaz.
    !
    ! we wish to solve
    !
    !   U_t + [F(U)]_x = H                                                                                    
    !
    ! we want U_{i+1/2}^{n+1/2} -- the interface values that are input
    ! to the Riemann problem through the faces for each zone.
    !
    ! Taylor expanding yields                                                                              
    !
    !  n+1/2                     dU           dU                                                         
    ! U          = U   + 0.5 dx  --  + 0.5 dt --                                                         
    !  i+1/2,j,L    i,j          dx           dt                                                         
    !                                                                                                 
    !                                                                                                 
    !                            dU             dF                                                     
    !            = U   + 0.5 dx  --  - 0.5 dt ( -- - H )                                               
    !               i,j          dx             dx                                                      
    !                                                                                                 
    !                                                                                                 
    !                             dU      dF^x                                                           
    !            = U   + 0.5 ( dx -- - dt ---- ) + 0.5 dt H                                              
    !               i,j           dx       dx                                                            
    !                                                                                                 
    !                                                                                                 
    !                                 dt       dU                                                        
    !            = U   + 0.5 dx ( 1 - -- A^x ) -- + 0.5 dt H                                             
    !               i,j               dx       dx                                                        
    !                                                                                                 
    !                                                                                                 
    !                               dt       _                                                           
    !            = U   + 0.5  ( 1 - -- A^x ) DU + 0.5 dt H                                               
    !               i,j             dx                                                                   
    !                                                                                                 
    !                                                                                                 
    !                    +----------+----------+ +---+---+                                               
    !                               |                |                                                   
    !                                                                                                 
    !                   this is the monotonized   source term                                            
    !                   central difference term                                                          
  

    !-------------------------------------------------------------------------
    ! convert to primitve variables
    !-------------------------------------------------------------------------




    !-------------------------------------------------------------------------
    ! compute the flattening coefficients
    !-------------------------------------------------------------------------




    !-------------------------------------------------------------------------
    ! compute the monotonized central differences
    !-------------------------------------------------------------------------




    !-------------------------------------------------------------------------
    ! compute left and right primitive variable states
    !-------------------------------------------------------------------------




    !-------------------------------------------------------------------------
    ! transform the states into conserved variables
    !-------------------------------------------------------------------------




    !-------------------------------------------------------------------------
    ! apply the source terms
    !-------------------------------------------------------------------------



    ! loop over interfaces and fill the states 
    do i = U%grid%lo, U%grid%hi+1
       U_l%data(i,:) = U%data(i-1,:)
       U_r%data(i,:) = U%data(i,:)
    enddo

  end subroutine make_interface_states_plm

end module interface_states_plm_module
