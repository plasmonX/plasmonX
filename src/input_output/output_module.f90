!> Target module
!!
!! This module contains the subroutines for the output type
!!
!! Date         : 2025
!!
module output_module

   use parameters_module
   use string_manipulation_module
   use array_manipulation_module

   implicit none

   !public variables
   public out_

   type out_type
      integer  :: unit_info = 8 ! unit of the info file
      integer  :: iunit = 12 ! unit of the file
      integer  :: ivrb  = 0  ! Verbose modality 
      logical  :: exists = .false.
      character (len=12)  :: sticks = "(1x,80(1h-))"
      character (len=200) :: filename
      character (len=200) :: info_file =""
      character (len=200) :: file_tar =""
      character (len=200) :: what

      contains

      procedure :: print_matrix
      procedure :: warning
      procedure :: error
      procedure :: clean_up_scratch
      procedure :: make_targz
   end type out_type
    
   type (out_type), Save :: out_

contains

   !> Subroutine for printing a rectangular matrix
   !!    Input  : out_     -- out_type
   !!    Input  : string   -- string defining the matrix
   !!    Input  : matrix   -- matrix to be printed (idim1,idim2)
   !!    Input  : idim1    -- first dimension
   !!    Input  : idim2    -- second dimension
   subroutine print_matrix(out_,string,matrix,idim1,idim2)
  
      implicit none
       
      !input/output variables
      class(out_type), intent(in) :: out_
      integer, intent(in)                             :: idim1
      integer, intent(in)                             :: idim2
      real(dp), dimension(idim1,idim2), intent(in)    :: matrix
      character(len=*)                                :: string
       
      !internal variables
      integer                                 :: i, j, k
      integer                                 :: ierr
      integer                                 :: NVec
      integer                                 :: NElements
      integer                                 :: len_string
      character(len=7)                        :: format_string
      character(len=13)                       :: format_1="(3x,5(8x,i5))"
      real(dp), dimension(:,:,:), allocatable :: dummy
      real(dp), dimension(:,:,:), allocatable :: dummy2
       
      !count how many vectors are there
      NVec      = ceiling(float(idim2)/five)
      NElements = mod(idim2,5)
      if(nelements.eq.0) nelements = 5
      !allocations of dummy matrices
      allocate(dummy(idim1,NVec-1,5),stat=ierr)
      if(ierr.gt.0) call out_%error('not enough space in print_matrix')
      allocate(dummy2(idim1,1,NElements),stat=ierr)
      if(ierr.gt.0) call out_%error('not enough space in print_matrix')
      !definition of dummy matrices
      if(NVec.ne.1) then
         !$omp parallel do collapse(3) 
         do i = 1, NVec-1
            do j = 1, idim1
               do k = 1, 5
                  dummy(j,i,k) = matrix(j,k+(i-1)*5)
               enddo
            enddo
         enddo
        !$omp end parallel do
      endif
      do j = 1, idim1
         do k = 1, nelements
            dummy2(j,1,k) = matrix(j,k+(NVec-1)*5)
         enddo
      enddo    
      !printing of the matrix
      len_string = len(string)
      write(format_string,'(a,i2,a)') "(",40-len_string/2,"X,a)"
      write(out_%iunit,out_%sticks)
      write(out_%iunit,format_string) string
      write(out_%iunit,out_%sticks)
      do i = 1, nvec
         write(out_%iunit,'(a)') ' '
         if(i.ne.nvec) write(out_%iunit,format_1) (k+(i-1)*5,k=1,5)
         if(i.eq.nvec) write(out_%iunit,format_1) (k+(i-1)*5,k=1,NElements)
         do j = 1, idim1
            if(i.ne.NVec) write(out_%iunit,'(i4,4x,5(e11.4,2x))') j, &
                                             (dummy(j,i,k), k = 1, 5 )
            if(i.eq.nvec) write(out_%iunit,'(i4,4x,5(e11.4,2x))') j, &
                                             (dummy2(j,1,k), k = 1,NElements)
         enddo
      enddo
      write(out_%iunit,out_%sticks)
      flush(out_%iunit)
      
   end subroutine print_matrix


   !> Subroutine for printing a warning
   !!    Input  : out_     -- out_type
   !!    Input  : string   -- string defining the warning
   subroutine warning(out_,string)
  
      implicit none
       
      !input/output variables
      class(out_type), intent(in)  :: out_
      character(len=*), intent(in) :: string
  
      !internal variables
      integer          :: unit_
      
      !check if the file is opened with unit_ = out_%iunit
      inquire(file=out_%filename,number=unit_)
      if (unit_.eq.out_%iunit) then
         write(out_%iunit,'(/1x,a/)') "Warning! "//trim(string)
      else
         write(*,'(/1x,a/)') "Warning! "//trim(string)
      endif
         
   end subroutine warning


   !> Subroutine for printing an error
   !!    Input  : out_     -- out_type
   !!    Input  : string   -- string defining the error
   subroutine error(out_,string)
  
      implicit none
       
      !input/output variables
      class(out_type)  :: out_
      character(len=*) :: string
  
      !internal variables
      integer          :: unit_
      
      !check if the file is opened with unit=out_%iunit
      inquire(file=out_%filename,number=unit_)
      if (unit_.eq.out_%iunit) then
         write(out_%iunit,'(/1x,a)') "Error during the execution of plasmonX"
         write(out_%iunit,'(1x,a/)') trim(string)
         flush(out_%iunit)
         stop
      else
         write(*,'(/1x,a)') "Error during the execution of plasmonX"
         write(*,'(1x,a/)') trim(string)
         stop
      endif
         
   end subroutine error


   !> Subroutine for cleaning up the backup files (*.plasmonX.bk)
   !!    Input  : out_     -- out_type
   subroutine clean_up_scratch(out_)
     
      implicit none
  
      class(out_type), intent(in) :: out_
  
      write(out_%iunit,out_%sticks)
      write(out_%iunit,'(/ a /)') ' I am cleaning up the backup (*.plasmonX.bk)'
      write(out_%iunit,out_%sticks)
      call system('rm *.plasmonX.bk')
      
   end subroutine clean_up_scratch


   !> Subroutine for creating the targz file if requested
   !!    Input  : out_     -- out_type
   subroutine make_targz(out_)
     
      implicit none 
  
      !input/output variables
      class(out_type), intent(in) :: out_
       
      !internal variables
      integer :: unit_freq
      integer :: unit_info = 8
      character(len=100)  :: line_file_freq
      character(len=100)  :: file_csv
      character(len=100)  :: file_info
      character(len=100)  :: file_tar
      character(len=318)  :: command
      logical             :: exist_csv
      logical             :: exist_freq =.false.
      logical             :: info_freq_exist =.false.
  
      write(file_info,'(a)') out_%filename(1:len_trim(out_%filename)-4)//'.info'
      write(file_csv,'(a)') out_%filename(1:len_trim(out_%filename)-4)//'.csv'
      Inquire(file=file_csv,exist=exist_csv)     
       
      !check if freq file exists
      call execute_command_line('for i in *.freq ; do test -f "$i" '// &
                                '&& echo "exists one or more files" > &
                                &info_freq.txt && break; done')
      inquire(file="info_freq.txt",exist=info_freq_exist)
      if(info_freq_exist) then
         unit_freq = 13
         Open(unit=unit_freq,file="info_freq.txt",status="OLD",ERR=02)
            read(unit_freq,'(a)') line_file_freq
            if(index(line_file_freq,"exists one or more files").gt.0) &
               exist_freq = .true.
         02 Continue      
         Close(unit_freq)
         call execute_command_line("rm info_freq.txt")
      endif
  
      !open unit_info for writing the correct execution
      Open(unit=unit_info,file=file_info,status="OLD",ACCESS='APPEND',ERR=03)
         write(unit_info,out_%sticks)
         write(unit_info,'(24x,a)') 'Normal termination of plasmonX'
         write(unit_info,out_%sticks)
      03 Continue      
      close(unit_info)
  
      !creation of the tar.gz (with all relevant files)
      write(file_tar,'(a)') out_%filename(1:len_trim(out_%filename)-4) // &
                            '.tar.gz'
      if(exist_csv) then 
         if(exist_freq) then
            command = "tar -czf "//file_tar//" "//file_info//" "//file_csv//&
                      " *.freq"
         else
            command = "tar -czf "//file_tar//" "//file_info//" "//file_csv
         endif
      else
         if(exist_freq) then
            command = "tar -czf "//file_tar//" "//file_info//" *.freq"
         else
            command = "tar -czf "//file_tar//" "//file_info
         endif
      endif
  
      call execute_command_line(trim(command))
      call execute_command_line("rm "//file_info)
      if(exist_freq) call execute_command_line("rm *freq")
  
      If(out_%ivrb.ge.1) then
         write(out_%iunit,'(1x,a)') "Created file tar  : "//file_tar
         write(out_%iunit,out_%sticks) 
      endif
  
   end subroutine make_targz

end module output_module
