!> Float manipulation module
!!
!! This module contains the functions for manipulating floats
!!
!! Date         : 2025
!!
module float_manipulation_module

   use parameters_module, only : dp, tol_float
   
     implicit none
     public is_equal

contains

   !> Function for comparing two floats
   !!    Input  : a, b 
   !!    Output : is_equal
   logical function is_equal(a,b)
 
      implicit none
 
      !input variables
      real(dp), intent(in)  :: a, b
 
      is_equal = dabs(a - b) .le. tol_float
 
   end function is_equal

end module float_manipulation_module
