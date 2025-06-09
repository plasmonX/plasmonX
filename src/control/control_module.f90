!> Control module
!!
!! This module contains the variables for controlling the calculations
!! and the type control
!!
!! Date         : 2024 
!!
module control_module
      
   implicit none
   public control

   type :: control_type
      logical :: principal_axis = .false. ! rotate axis?
      logical :: restart        = .false. ! restart calculation
      logical :: save_info      = .false. ! if you should save the info file
      ! default: solution for maxima analysis
      character(len=200) :: maxima_analysis = "absorption"  

      contains 

      procedure :: print_info         => print_info_control
   end type control_type
    
   type (control_type), save :: control

contains

   !> Subroutine for writing the info related to control and out_module
   !!    Input  : control      -- control_type
   subroutine print_info_control(control)
  
      use output_module
      
      implicit none
       
      !input/output variables
      class(control_type), intent(in) :: control
  
      write(out_%iunit,'(1x,a)') "What               : "//trim(out_%what)
      write(out_%iunit,'(1x,a,i1)') "Verbose            : ", out_%ivrb
      write(out_%iunit,'(1x,a,l1)') "Principal Axis     : ", &
                                    control%principal_axis
      write(out_%iunit,'(1x,a,l1)') "Restart            : ",control%restart
      write(out_%iunit,'(1x,a,l1)') "Save info          : ",control%save_info
      if (control%save_info) then
         write(out_%iunit,'(1x,a)') "Info File          : "// &
                                    trim(out_%info_file)
         write(out_%iunit,'(1x,a)') "tar.gz File        : "//trim(out_%file_tar)
      else
         write(out_%iunit,'(1x,a)') "Info File          : none"
         write(out_%iunit,'(1x,a)') "tar.gz File        : none"
      endif
      
      write(out_%iunit,out_%sticks) 
      flush(out_%iunit)
         
   end subroutine print_info_control


   !> Subroutine for final deallocation
   subroutine final_deallocation
  
      use field_module
      use target_module
      
      implicit none
       
      call field%dealloc()
      call target_%dealloc()
         
   end subroutine final_deallocation
    
end module control_module
