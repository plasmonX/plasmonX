!> Control analysis module   
!!
!! This module contains the control parameters
!!
!! Date         : 2025
!!
module control_analysis_module

   use control_module
   use parameters_module

   implicit none
   public control_analysis

   type, extends(control_type) :: control_analysis_type
      character (len=100) :: what
      logical             :: plane_requested = .false.
      character (len=2)   :: plane
      integer             :: n_plane
      real(dp)            :: step_plane
      real(dp)            :: start_plane
      real(dp)            :: scale_e0
      logical             :: volume_requested = .false.
      logical             :: separate_q_mu = .false.
      integer             :: n_threads_omp = 1

      contains

      procedure :: print_info         => print_info_control_analysis
   end type control_analysis_type
    
   type (control_analysis_type), save :: control_analysis

contains

   !> Subroutine for writing the info related to control and out_module
   !!    Input  : control      -- control_type
   subroutine print_info_control_analysis(control)
  
      use output_module
      
      implicit none
       
      !input/output variables
      class(control_analysis_type), intent(in) :: control
  
      write(out_%iunit,'(1x,a)') "What               : "//trim(control%what)
      if (trim(control%what).eq.'density'.or.trim(control%what).eq.'field') then
         write(out_%iunit,'(1x,a, l1)') "Plane requested    : ", &
                                     control%plane_requested
         if(control%plane_requested) then
            write(out_%iunit,'(1x,a,i7)') "Number of planes   : ", &
                                        control%n_plane
            write(out_%iunit,'(1x,a,f7.2,a)') "Start plane        : ", &
                                        control%start_plane, " Ang."
            write(out_%iunit,'(1x,a,f7.2,a)') "Step btw planes    : ", &
                                        control%step_plane, " Ang."
         else
            write(out_%iunit,'(1x,a, l1)') "Volume requested   : ", &
                                        control%volume_requested
         endif
      endif
      
      write(out_%iunit,out_%sticks) 
      flush(out_%iunit)
         
   end subroutine print_info_control_analysis

    
end module control_analysis_module
