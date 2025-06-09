!> Main Program plasmonX_analysis
!!
!! Date: 2024
!!
program plasmonX_analysis

use parameters_module
use input_analysis_module
use plot_module
!$ use omp_lib
  
implicit none

!internal variables
integer :: iost

call inp_analysis%check_file()
call inp_analysis%get_output()
call out_analysis%create_output_folder()

open(unit=out_analysis%iunit,file=out_analysis%filename,status="UNKNOWN", &
     iostat=iost,position="append")
   if (iost.ne.0) call out_analysis%error('Error in opening(w) '//&
                                          out_analysis%filename// ' file')
   call inp_analysis%read_file()

   !$ call omp_set_num_threads(control_analysis%n_threads_omp) 

   call inp_analysis%extract_tar_gz()
   call check_tmp_files()
   call read_info_file()
   call inp_analysis%print_input_info()

   !xyz file
   if(control_analysis%what.eq.'xyz') then
      call plot_xyz()
   !pqr file   
   else if (control_analysis%what.eq.'pqr') then
      call plot_pqr()
   !field or density requested   
   else if (control_analysis%what.eq.'field'.or. &
            control_analysis%what.eq.'density') then
      call plot_field_or_density()
   endif
close(out_analysis%iunit)

call execute_command_line("rm -rf tmp_analysis")
  
end program plasmonX_analysis
