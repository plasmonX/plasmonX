   !> Subroutine for constructing the accurate Green function 
   !!    In/Out  : bem      -- bem_type
   module subroutine construct_accurate_green_function(bem)

      implicit none
      
      !input/output variables
      class(bem_type), intent(inout) :: bem

   end subroutine construct_accurate_green_function


   !> Subroutine for flat integration diagonal part
   !! This is a subroutine valid for flat tesserae only in polar coordinates.
   !!    In/Out  : bem      -- bem_type
   module subroutine flat_integration_diagonal(bem)
  
      implicit none
       
      !input/output variables
      class(bem_type), intent(inout) :: bem
  
   end subroutine flat_integration_diagonal


   !> Subroutine for quadrature rules setting for an unit triangular mesh
   !! Strang and Fix, "An Analysis of the Finite Element Method"
   !! Prentice Hall, 1973.
   !! Other quadrature rule can be implemented, with different precisions  
   !!    In/Out  : bem      -- bem_type
   module subroutine quadrature_rule(bem)
  
      implicit none
       
      !input/output variables
      class(bem_type), intent(inout) :: bem
  
   end subroutine quadrature_rule


   !> Subroutine for the integration over tesserae in polar coordinates
   !!    In/Out  : bem      -- bem_type
   !!    In/Out  : int_coord_x -- integration coordinates x
   !!    In/Out  : int_coord_y -- integration coordinates y
   module subroutine tesserae_integration_polar_coordinates(bem,int_coord_x, &
                                                                int_coord_y)
  
      implicit none
       
      !input/output variables
      class(bem_type), intent(inout) :: bem
      real(dp), dimension(bem%npol1*bem%npol2*3), intent(inout) :: int_coord_x
      real(dp), dimension(bem%npol1*bem%npol2*3), intent(inout) :: int_coord_y
  
   end subroutine tesserae_integration_polar_coordinates


   !> Subroutine for flat integration off diagonal part
   !!    In/Out  : bem      -- bem_type
   module subroutine flat_integration_off_diagonal(bem)
  
      implicit none
  
      !input/output variables
      class(bem_type), intent(inout) :: bem
  
   end subroutine flat_integration_off_diagonal


   !> Subroutine for selecting the refined elements
   !!    In/Out  : bem      -- bem_type
   module subroutine selection_of_refined_elements(bem)
  
      implicit none
       
      !input/output variables
      class(bem_type), intent(inout) :: bem
  
   end subroutine selection_of_refined_elements


   !> Subroutine for refining the diagonal elements
   !!    In/Out  : bem      -- bem_type
   module subroutine refinement_diagonal_elements(bem)
  
      implicit none
       
      !input/output variables
      class(bem_type), intent(inout) :: bem
  
   end subroutine refinement_diagonal_elements


   !> Subroutine for refining the off diagonal elements
   !!    Input   : bem                  -- bem_type
   !!    Input   : i_row                -- index of the row
   !!    Input   : j_column             -- index of the column
   !!    Output  : off_diagonal_element -- i,j element
   module function refinement_off_diagonal_elements(bem,i_row,j_column) &
                   result(off_diagonal_element)
  
      implicit none
       
      !input/output variables
      class(bem_type), intent(in) :: bem
      integer, intent(in)         :: i_row
      integer, intent(in)         :: j_column
      real(dp)                    :: off_diagonal_element
  
   end function refinement_off_diagonal_elements


   !> Subroutine for the final refinement of the Green function
   !!    In/Out  : bem                  -- bem_type
   !!    Input   : i_row                -- index of the row
   !!    Input   : j_column             -- index of the column
   !!    Output  : off_diagonal_element -- i,j element
   module subroutine final_refinement_green_function(bem, matrix_constant)
  
      implicit none
       
      !input/output variables
      class(bem_type) :: bem
      real(dp), dimension(bem%n_var,bem%n_var), intent(inout) :: matrix_constant
  
   end subroutine final_refinement_green_function
