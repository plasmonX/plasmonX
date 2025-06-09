!> Output analysis module   
!!
!! This module contains the subroutines for printing the output of the analysis
!!
!! Date         : 2024
!!
module output_analysis_module
 
   use parameters_module
   use output_module
   use target_module
   use control_analysis_module
 
   implicit none 
 
   !public variables
   public out_analysis 
 
   type, extends(out_type) :: out_analysis_type
     character(len=1000) :: plane_folder
     character(len=1000) :: root_folder
     character(len=1000) :: root_filename
     character(len=13)   :: tmp_folder = "tmp_analysis/"
 
     contains
 
     procedure :: out_file_fill => out_file_fill_analysis
     procedure :: create_output_folder
     procedure :: create_plane_folder
   end type out_analysis_type
    
   type (out_analysis_type) :: out_analysis
 
contains

   !> Subroutine for define the name of the output file
   !!    In/Out : out_         -- type
   !!    Input  : in_file      -- input_file
   subroutine out_file_fill_analysis(out_,in_file)
  
      implicit none
       
      !input/output variables
      class(out_analysis_type), intent(inout) :: out_
      character(len=*), intent(in) :: in_file
  
      !internal variables
      integer :: nlen
      
      nlen = len_trim(in_file)
      out_%root_filename = trim(in_file(1:nlen-7))
      write(out_%filename,'(a)') trim(out_%root_filename) // '_analysis.log'
      out_%filename = trim(out_%filename)
         
   end subroutine out_file_fill_analysis


   !> Subroutine for creating the output folder
   !!    In/Out : out_         -- type
   subroutine create_output_folder(out_)
  
      implicit none
       
      !input/output variables
      class(out_analysis_type), intent(inout) :: out_
  
      !internal variables
      character(len=200)  :: input_file
      character(len=1000) :: command_folder
       
      write(input_file,'(a)') out_%filename(1:len_trim(out_%filename)-13)
  
      command_folder ='[ ! -d "post_process_'//trim(input_file)// &
                       '" ] && mkdir -p "post_process_'//trim(input_file)//'"'
  
      call execute_command_line(trim(command_folder))
  
      out_%root_folder = trim('post_process_'//trim(input_file)//'/')
       
      write(out_%filename,'(a)') trim(out_%root_folder)//trim(out_%filename)
  
   end subroutine create_output_folder


   !> Subroutine for creating the plane folder
   !!    In/Out : out_         -- type
   subroutine create_plane_folder(out_)
  
      implicit none
       
      !input/output variables
      class(out_analysis_type), intent(inout) :: out_
  
      !internal variables
      character(len=1000) :: command_folder
       
      command_folder ='[ ! -d "post_process_'//trim(out_%root_filename)// &
                       '/planes" ] && mkdir -p "post_process_'// &
                       trim(out_%root_filename)//'/planes"'
  
      call execute_command_line(trim(command_folder))
  
      command_folder ='[ ! -d "post_process_'//trim(out_%root_filename)// &
                       '/planes/'//trim(control_analysis%plane)//&
                       '" ] && mkdir -p "post_process_'// &
                       trim(out_%root_filename)//&
                       '/planes/'//trim(control_analysis%plane)//'"'
  
      call execute_command_line(trim(command_folder))
  
      out_%plane_folder = trim('post_process_'//trim(out_%root_filename)// &
                               '/planes/'//trim(control_analysis%plane)//'/')
       
   end subroutine create_plane_folder

end module output_analysis_module


     !
     ! We SHOULD PRINT THIS
     !
     ! ! Stampa dei parametri letti per verifica
     ! print *, 'what = ', control_analysis%what
     ! print *, 'convolution_what = ', control_analysis%convolution_what
     ! print *, 'fwhm = ', control_analysis%fwhm
     ! print *, 'nofield = ', control_analysis%nofield
     ! print *, 'num_ex_freq = ', field%n_freq
     ! if (allocated(field%freq)) &
     !     print *, 'frequencies = ', field%freq
     ! print *, 'min_freq = ', field%min_freq
     ! print *, 'max_freq = ', field%max_freq
     ! print *, 'step = ', field%step_freq
     ! print *, 'plane = ', control_analysis%plane
     ! print *, 'n_plane = ', control_analysis%n_plane
     ! print *, 'step_plane = ', control_analysis%step_plane
     ! print *, 'start = ', control_analysis%start_plane
     ! print *, 'coord_plane = ', control_analysis%coord_plane
     ! print *, 'field_dir = ', field%polarization
     ! print *, 'point_charge = ', control_analysis%point_charge
     ! print *, 'nx_points = ', control_analysis%nx_points
     ! print *, 'ny_points = ', control_analysis%ny_points
     ! print *, 'nz_points = ', control_analysis%nz_points
     ! print *, 'volume = ', control_analysis%volume
     ! print *, 'min_vol = ', control_analysis%min_vol
     ! print *, 'max_vol = ', control_analysis%max_vol
     ! print *, 'index_atom = ', control_analysis%index_atom

