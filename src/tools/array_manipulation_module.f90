!> Array manipulation module
!!
!! This module contains the subroutines and functions for manipulating arrays
!!
!! Date         : 2025
!!
module array_manipulation_module

   use parameters_module, only : zero, dp
   
     implicit none
     public cross_product,       &
            array_clear,         &
            array_clear_complex, &
            array_scale,         &
            array_scale_complex

contains

   !> Function for performing the cross product between two real arrays
   !!    Input  : a, b 
   !!    Output : cross_product
   function cross_product(a, b)
 
      implicit none
 
      !input variables
      real(dp), dimension(3), intent(in)  :: a, b
      real(dp), dimension(3) :: cross_product
 
      cross_product(1) = a(2)*b(3) - a(3)*b(2)
      cross_product(2) = a(3)*b(1) - a(1)*b(3)
      cross_product(3) = a(1)*b(2) - a(2)*b(1)
 
   end function cross_product


   !> Subroutine for clearing an array of length nlen
   !!    Input  : nlen
   !!    In/Out : a  
   subroutine array_clear(nlen,a)
 
      implicit none
       
      !input variables
      integer, intent(in)                      :: nlen
      real(dp), dimension(nlen), intent(inout) :: a
 
      !internal variables
      integer :: i
 
      !$omp parallel do
      do i = 1, nlen
         a(i) = zero
      enddo
      !$omp end parallel do
 
   end subroutine array_clear


   !> Subroutine for clearing an array of length nlen (complex)
   !!    Input  : nlen
   !!    In/Out : a  
   subroutine array_clear_complex(nlen,a)
 
      implicit none
       
      !input variables
      integer, intent(in)                         :: nlen
      complex(dp), dimension(nlen), intent(inout) :: a
 
      !internal variables
      integer :: i
 
      !$omp parallel do
      do i = 1, nlen
         a(i) = dcmplx(zero,zero)
      enddo
      !$omp end parallel do
 
   end subroutine array_clear_complex


   !> Subroutine for scaling an array of length nlen by factor
   !!    Input  : nlen
   !!    Input  : factor
   !!    In/Out : a  
   subroutine array_scale(factor,nlen,a)
  
      implicit none
       
      !input variables
      integer , intent(in)                     :: nlen
      real(dp), intent(in)                     :: factor
      real(dp), dimension(nlen), intent(inout) :: a
      
      !internal variables
      integer :: i 
       
      !$omp parallel do private(i)
      do i = 1, nlen
         a(i) = a(i) * factor
      enddo
      !$omp end parallel do
     
   end subroutine array_scale


   !> Subroutine for scaling an array of length nlen by factor (complex)
   !!    Input  : nlen
   !!    Input  : factor
   !!    In/Out : a  
   subroutine array_scale_complex(factor,nlen,a)
  
      implicit none
       
      !input variables
      integer , intent(in)                        :: nlen
      complex(dp), intent(in)                     :: factor
      complex(dp), dimension(nlen), intent(inout) :: a
      
      !internal variables
      integer :: i 
       
      !$omp parallel do private(i)
      do i = 1, nlen
         a(i) = a(i) * factor
      enddo
      !$omp end parallel do
     
   end subroutine array_scale_complex

end module array_manipulation_module
