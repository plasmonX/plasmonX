!> Main Program plasmonX
!!
!! Date : 2025
program plasmonX

use input_module
use output_module
use target_module
use control_module
use memory_manager_module
use parameters_module
!$ use omp_lib

implicit none

!internal variables
integer :: iost

call inp_%get_arguments()
call inp_%check_input_file()
call inp_%read_general_settings()

!$ call omp_set_num_threads(n_threads_omp) 

open(unit=out_%iunit,file=out_%filename,status="unknown",position="append",&
     iostat=iost)

   if(iost.ne.0) then
      write(*,'(a)') "Error opening output file: "//trim(out_%filename)
      stop
   endif

   call inp_%read_() 
   call inp_%print_input_info()

   if(control%save_info) call initialize_save_info()
   if(control%principal_axis) &
      call target_%rotate_principal_axis(control%save_info)

   if(field%static) then 
      call algorithm%solve_static_field()
   else if(field%dynamic) then
      call algorithm%solve_dynamic_field()
   else 
      call algorithm%solve_ground_state()
   endif

   if(control%save_info) call out_%make_targz()

   call final_deallocation() !inside control
   call mem_man%print_status()

close(out_%iunit)

end program plasmonX
