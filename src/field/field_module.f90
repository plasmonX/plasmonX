!> Field module
!!
!! This module contains the subroutines and variables related to the type
!! field, which is used for all info about the external field 
!!
!! Date         : 2025
!!
module field_module

   !$ use omp_lib
   use output_module
   use parameters_module
   use float_manipulation_module
 
   implicit none
 
   !public type
   public field
 
   !field type
   type field_type
      !parameters of the external frequencies
      integer  :: n_freq    = 0 ! 
      real(dp) :: min_freq  = zero  ! min frequency in au
      real(dp) :: max_freq  = zero  ! max frequency in au
      real(dp) :: step_freq = zero  ! step/graining in au
      
      !list of field frequencies saved in nm
      real(dp), dimension(:), allocatable :: freq
 
      !type of field
      logical  :: static  = .false.
      logical  :: dynamic = .false.
      real(dp) :: e_0     = 0.0d0 ! field intensity in a.u. default: 10^-5
 
      !rhs form: potential or field
      character(len=200) :: rhs_form = "none"
 
      !field polarization
      integer                        :: n_polarization = 3   !default x,y,z
      character(len=4)               :: polarization = 'none' !'x','y','z'
      integer                        :: polarization_index = 0 ! 1, 2, 3
      character(len=2), dimension(3) :: polarization_name 
   
      contains 
  
      procedure :: print_info   => print_info_field
      procedure :: dealloc      => deallocate_field
      procedure :: set_frequencies
      procedure :: set_polarization_variables
      procedure :: index_rhs_polarization
   end type field_type
    
   type (field_type), save :: field

contains

   !> Subroutine for writing the field info in output
   !!    Input  : field        -- type
   subroutine print_info_field(field, analysis)
 
      implicit none
  
      !input/output variables
      class(field_type), intent(in) :: field
      logical, optional, intent(in) :: analysis

      !internal variables
      integer :: i
      real(dp) :: freq_min
      real(dp) :: freq_max
      real(dp) :: freq_step
      real(dp) :: freq_eV
      logical  :: analysis_internal

      call freqautoev(field%min_freq,freq_min)
      call freqautoev(field%max_freq,freq_max)
      call freqautoev(field%step_freq,freq_step)

      if(present(analysis)) then
         analysis_internal = analysis
         if(analysis) then
            freq_min = field%min_freq
            freq_max = field%min_freq
            freq_step = field%step_freq
         endif
      else
         analysis_internal = .false.
      endif
  
      if (field%static) then
         write(out_%iunit,'(1x,a)') "Field type         : static"
      else if (field%dynamic) then
         write(out_%iunit,'(1x,a)') "Field type         : dynamic"
      else
         write(out_%iunit,'(1x,a)') "Field type         : none"
      endif
      write(out_%iunit,'(1x,a,e11.4,a)') "Field Intensity    :",field%e_0," au"
      write(out_%iunit,'(1x,a)') "Field Polarization : "//&
                                 trim(field%polarization)
      write(out_%iunit,'(1x,a,i5)') "Num. Freq.         :", field%n_freq
      if(field%n_freq.gt.0) then
         write(out_%iunit,'(1x,a,e11.4,a)') "Min. Freq.         :", &
                                            freq_min," eV"
         write(out_%iunit,'(1x,a,e11.4,a)') "Max. Freq.         :", &
                                            freq_max," eV"
         write(out_%iunit,'(1x,a,e11.4,a)') "Step Freq.         :", &
                                            freq_step," eV"
         if(field%n_freq.gt.0.and.allocated(field%freq).and.&
            is_equal(field%min_freq,zero).and. &
            is_equal(field%max_freq,zero)) then!external freq
            do i = 1, field%n_freq
               if(analysis_internal) then 
                  freq_eV = field%freq(i)
               else
                  call freqautoev(field%freq(i),freq_eV)
               endif
               write(out_%iunit,'(1x,a,i5,a,e10.3,a)') "Freq ",i,"         :", &
                                                       freq_eV, " eV"
            enddo
         else !defined with min max step
            if (out_%ivrb.ge.2) then
               do i = 1, field%n_freq
                  call freqautoev(field%freq(i), freq_eV)
                  write(out_%iunit,'(1x,a,i5,a,e10.3,a)') "Freq ",i,&
                                                          "         :",&
                                                          freq_eV, " eV"
               enddo
            endif
         endif
      endif
      write(out_%iunit,out_%sticks)
      flush(out_%iunit)
 
   end subroutine print_info_field


   !> Subroutine for deallocating the field variables
   !!    In/Out : field        -- type
   subroutine deallocate_field(field)
 
      use memory_manager_module

      implicit none
  
      !input/output variables
      class(field_type), intent(inout) :: field

      if (allocated(field%freq)) &
         call mem_man%dealloc(field%freq,"field%freq")
 
   end subroutine deallocate_field


   !> Subroutine for setting the frequencies in input
   !!    In/out : field        -- type
   subroutine set_frequencies(field)
 
      implicit none
  
      !input/output variables
      class(field_type), intent(inout) :: field
  
      !internal variables
      integer :: i
      logical  :: explicit = .false.
  
      if(is_equal(field%min_freq,zero).and.is_equal(field%max_freq,zero)) &
         explicit = .true.
      if(.not.explicit) then !frequencies defined within a range
         do i = 1, field%n_freq
            field%freq(i) = field%min_freq + field%step_freq*(i-1)
         enddo 
      endif
 
   end subroutine set_frequencies


   !> Subroutine for setting the polarization variables (number, name)
   !!    In/out : field_       -- type
   subroutine set_polarization_variables(field)
     
      implicit none

      !input variables
      class(field_type), intent(inout) :: field
  
      if(field%polarization.eq.'all') then !default
         field%n_polarization       = 3 
         field%polarization_name(1) = "Ex"
         field%polarization_name(2) = "Ey"
         field%polarization_name(3) = "Ez"
      else if (field%polarization.eq.'x') then
         field%n_polarization       = 1
         field%polarization_name(1) = "Ex"
         field%polarization_index   = 1 
      else if (field%polarization.eq.'y') then
         field%n_polarization       = 1
         field%polarization_name(1) = "Ey"
         field%polarization_index   = 2
      else if (field%polarization.eq.'z') then
         field%n_polarization       = 1
         field%polarization_name(1) = "Ez"
         field%polarization_index   = 3
      endif
  
   end subroutine set_polarization_variables


   !> Function to define the index for the RHS based on polarzation
   !!    In/out : field       -- type
   !!    Input  : int_         -- integer for polarization
   integer function index_rhs_polarization(field, int_) result(index_)

      implicit none
 
      !input variables
      class(field_type), intent(inout) :: field
      integer, intent(in)              :: int_ !Polarization integer 
 
      index_ = int_
      if(field%polarization.eq.'all') then 
         index_ = int_ 
      else if(field%polarization.eq.'x') then
         index_ = 1 
      else if(field%polarization.eq.'y') then
         index_ = 2 
      else if(field%polarization.eq.'z') then
         index_ = 3
      endif
 
   end function index_rhs_polarization

end module field_module
