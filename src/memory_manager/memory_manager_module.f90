!> Memory manager module
!!
!! This module contains the variables for controlling the memory allocation
!!
!! Date         : 2025
!!
module memory_manager_module

   use output_module
   use parameters_module
   use array_manipulation_module
   use float_manipulation_module

   implicit none
   private

   !parameters
   integer, parameter :: max_arrays = 1000  ! Maximum number of trackable arrays

   !information about the arrays
   type :: array_info
      !name
      character(len=128) :: name_ = ""  
      !size 
      integer            :: n_elements = 0
      !size in byte
      real(dp)           :: size_bytes = zero
      !allocation flag
      logical            :: is_allocated = .false. 
   end type array_info

   !public variables
   public :: mem_man

   !define control_type for memory management
   type :: memory_manager_type
      !memory used when called in bytes
      real(dp)  :: memory_used = zero
      !memory available in GB
      real(dp)  :: memory_available = zero
      !peak memory usage in GB
      real(dp) :: peak_memory_used = zero
      !number of allocated arrays
      integer  :: num_arrays = 0 
      !array info storage
      type(array_info), dimension(max_arrays) :: arrays  
      !List of known OpenMP private arrays
      character(len=15), dimension(7) :: omp_private_arrays = &
                                        [ "matrix_constant", &
                                          "rhs_w          ", &
                                          "rhs_1          ", &
                                          "matrix_w       ", &
                                          "variables_w    ", &
                                          "variables_1    ", &
                                          "polar_w        " ]      

      contains

         procedure :: check_mem 
         procedure :: update_saved_info
         procedure :: check_deallocate
         procedure :: allocate_array_real_1D
         procedure :: allocate_array_real_2D
         procedure :: allocate_array_real_3D
         procedure :: allocate_array_complex_1D
         procedure :: allocate_array_complex_2D
         procedure :: allocate_array_int_1D
         procedure :: allocate_array_int_2D
         procedure :: allocate_array_char_1D
         procedure :: allocate_array_char_2D
         procedure :: print_status
         generic   :: alloc   => allocate_array_real_1D,    &
                                 allocate_array_real_2D,    &
                                 allocate_array_real_3D,    &
                                 allocate_array_complex_1D, &
                                 allocate_array_complex_2D, &
                                 allocate_array_int_1D,     &
                                 allocate_array_int_2D,     &
                                 allocate_array_char_1D,    &
                                 allocate_array_char_2D
         procedure :: deallocate_array_real_1D
         procedure :: deallocate_array_real_2D
         procedure :: deallocate_array_real_3D
         procedure :: deallocate_array_complex_1D
         procedure :: deallocate_array_complex_2D
         procedure :: deallocate_array_int_1D
         procedure :: deallocate_array_int_2D
         procedure :: deallocate_array_char_1D
         procedure :: deallocate_array_char_2D
         generic   :: dealloc   => deallocate_array_real_1D,    &
                                   deallocate_array_real_2D,    &
                                   deallocate_array_real_3D,    &
                                   deallocate_array_complex_1D, &
                                   deallocate_array_complex_2D, &
                                   deallocate_array_int_1D,     &
                                   deallocate_array_int_2D,     &
                                   deallocate_array_char_1D,    &
                                   deallocate_array_char_2D
   end type memory_manager_type

   type(memory_manager_type), save :: mem_man

contains

   !> Subroutine for allocating a generic array real 1D
   !!    In/Out : mem          -- memory_type
   !!    Output : ptr          -- array to allocate
   !!    Input  : name_        -- string identifying the array
   !!    Input  : n1           -- number of elements 1D array
   subroutine allocate_array_real_1D(mem, ptr, n1, name_)

      !input/output variables
      class(memory_manager_type), intent(inout) :: mem
      real(dp), allocatable, intent(out) :: ptr(:)
      integer, intent(in) :: n1
      character(len=*), intent(in) :: name_

      call mem%check_mem(name_,n1,'real')

      !allocate the array
      allocate(ptr(n1))
      call array_clear(n1, ptr)

      !update saved info
      call mem%update_saved_info(name_,n1,'real')

   end subroutine allocate_array_real_1D


   !> Subroutine for allocating a generic array real 2D
   !!    In/Out : mem          -- memory_type
   !!    Output : ptr          -- array to allocate
   !!    Input  : name_        -- string identifying the array
   !!    Input  : n1           -- number of dim1 2D array
   !!    Input  : n2           -- number of dim2 2D array
   subroutine allocate_array_real_2D(mem, ptr, n1, n2, name_)

      !input/output variables
      class(memory_manager_type), intent(inout) :: mem
      real(dp), allocatable, intent(out) :: ptr(:,:)
      integer, intent(in) :: n1
      integer, intent(in) :: n2
      character(len=*), intent(in) :: name_

      call mem%check_mem(name_,n1*n2,'real')

      !allocate the array
      allocate(ptr(n1,n2))
      call array_clear(n1*n2, ptr)

      !update saved info
      call mem%update_saved_info(name_,n1*n2,'real')

   end subroutine allocate_array_real_2D


   !> Subroutine for allocating a generic array real 3D
   !!    In/Out : mem          -- memory_type
   !!    Output : ptr          -- array to allocate
   !!    Input  : name_        -- string identifying the array
   !!    Input  : n1           -- number of dim1 3D array
   !!    Input  : n2           -- number of dim2 3D array
   !!    Input  : n2           -- number of dim3 3D array
   subroutine allocate_array_real_3D(mem, ptr, n1, n2, n3, name_)

      !input/output variables
      class(memory_manager_type), intent(inout) :: mem
      real(dp), allocatable, intent(out) :: ptr(:,:,:)
      integer, intent(in) :: n1
      integer, intent(in) :: n2
      integer, intent(in) :: n3
      character(len=*), intent(in) :: name_

      call mem%check_mem(name_,n1*n2*n3,'real')

      !allocate the array
      allocate(ptr(n1,n2,n3))
      call array_clear(n1*n2*n3, ptr)

      !update saved info
      call mem%update_saved_info(name_,n1*n2*n3,'real')

   end subroutine allocate_array_real_3D


   !> Subroutine for allocating a generic array complex 1D
   !!    In/Out : mem          -- memory_type
   !!    Output : ptr          -- array to allocate
   !!    Input  : name_        -- string identifying the array
   !!    Input  : n1           -- number of elements 1D array
   subroutine allocate_array_complex_1D(mem, ptr, n1, name_)

      !input/output variables
      class(memory_manager_type), intent(inout) :: mem
      complex(dp), allocatable, intent(out) :: ptr(:)
      integer, intent(in) :: n1
      character(len=*), intent(in) :: name_

      call mem%check_mem(name_,n1,'complex')

      !allocate the array
      allocate(ptr(n1))
      call array_clear_complex(n1, ptr)

      !update saved info
      call mem%update_saved_info(name_,n1,'complex')

   end subroutine allocate_array_complex_1D


   !> Subroutine for allocating a generic array complex 2D
   !!    In/Out : mem          -- memory_type
   !!    Output : ptr          -- array to allocate
   !!    Input  : name_        -- string identifying the array
   !!    Input  : n1           -- number of dim1 2D array
   !!    Input  : n2           -- number of dim2 2D array
   subroutine allocate_array_complex_2D(mem, ptr, n1, n2, name_)

      !input/output variables
      class(memory_manager_type), intent(inout) :: mem
      complex(dp), allocatable, intent(out) :: ptr(:,:)
      integer, intent(in) :: n1
      integer, intent(in) :: n2
      character(len=*), intent(in) :: name_

      call mem%check_mem(name_,n1*n2,'complex')

      !allocate the array
      allocate(ptr(n1,n2))
      call array_clear_complex(n1*n2, ptr)

      !update saved info
      call mem%update_saved_info(name_,n1*n2,'complex')

   end subroutine allocate_array_complex_2D


   !> Subroutine for allocating a generic array complex 1D
   !!    In/Out : mem          -- memory_type
   !!    Output : ptr          -- array to allocate
   !!    Input  : name_        -- string identifying the array
   !!    Input  : n1           -- number of elements 1D array
   subroutine allocate_array_int_1D(mem, ptr, n1, name_)

      !input/output variables
      class(memory_manager_type), intent(inout) :: mem
      integer, allocatable, intent(out) :: ptr(:)
      integer, intent(in) :: n1
      character(len=*), intent(in) :: name_

      call mem%check_mem(name_,n1,'integer')

      !allocate the array
      allocate(ptr(n1))
      ptr = 0

      !update saved info
      call mem%update_saved_info(name_,n1,'integer')

   end subroutine allocate_array_int_1D


   !> Subroutine for allocating a generic array complex 2D
   !!    In/Out : mem          -- memory_type
   !!    Output : ptr          -- array to allocate
   !!    Input  : name_        -- string identifying the array
   !!    Input  : n1           -- number of dim1 2D array
   !!    Input  : n2           -- number of dim2 2D array
   subroutine allocate_array_int_2D(mem, ptr, n1, n2, name_)

      !input/output variables
      class(memory_manager_type), intent(inout) :: mem
      integer, allocatable, intent(out) :: ptr(:,:)
      integer, intent(in) :: n1
      integer, intent(in) :: n2
      character(len=*), intent(in) :: name_

      call mem%check_mem(name_,n1*n2,'integer')

      !allocate the array
      allocate(ptr(n1,n2))
      ptr = 0

      !update saved info
      call mem%update_saved_info(name_,n1*n2,'integer')

   end subroutine allocate_array_int_2D


   !> Subroutine for allocating a generic array complex 1D
   !!    In/Out : mem          -- memory_type
   !!    Output : ptr          -- array to allocate
   !!    Input  : name_        -- string identifying the array
   !!    Input  : n1           -- number of elements 1D array
   subroutine allocate_array_char_1D(mem, ptr, n1, name_)

      !input/output variables
      class(memory_manager_type), intent(inout) :: mem
      character(len=*), allocatable, intent(out) :: ptr(:)
      integer, intent(in) :: n1
      character(len=*), intent(in) :: name_

      call mem%check_mem(name_,n1*len(ptr),'character')

      !allocate the array
      allocate(ptr(n1))
      ptr = ""

      !update saved info
      call mem%update_saved_info(name_,n1*len(ptr),'character')

   end subroutine allocate_array_char_1D


   !> Subroutine for allocating a generic array complex 2D
   !!    In/Out : mem          -- memory_type
   !!    Output : ptr          -- array to allocate
   !!    Input  : name_        -- string identifying the array
   !!    Input  : n1           -- number of dim1 2D array
   !!    Input  : n2           -- number of dim2 2D array
   subroutine allocate_array_char_2D(mem, ptr, n1, n2, name_)

      !input/output variables
      class(memory_manager_type), intent(inout) :: mem
      character(len=*), allocatable, intent(out) :: ptr(:,:)
      integer, intent(in) :: n1
      integer, intent(in) :: n2
      character(len=*), intent(in) :: name_

      call mem%check_mem(name_,n1*n2*len(ptr),'character')

      !allocate the array
      allocate(ptr(n1,n2))
      ptr = ""

      !update saved info
      call mem%update_saved_info(name_,n1*n2*len(ptr),'character')

   end subroutine allocate_array_char_2D


   !> Subroutine for checking the memory for array "name_"
   !!    In/Out : mem          -- memory_type
   !!    Input  : name_        -- string identifying the array
   !!    Input  : size_        -- the size of the array [int]
   !!    Input  : type_        -- the type of the array 
   subroutine check_mem(mem, name_, size_, type_)

      !input/output variables
      class(memory_manager_type), intent(inout) :: mem
      character(len=*), intent(in) :: name_
      integer, intent(in) :: size_
      character(len=*), intent(in) :: type_

      !internal variables
      integer :: i
      real(dp) :: size_adj
      real(dp) :: byte_type
      real(dp) :: size_bytes

      if ( type_.eq.'real') then
         byte_type = 8.0d0
      else if (type_.eq.'complex') then
         byte_type = 16.0d0
      else if (type_.eq.'integer') then
         byte_type = 4.0d0
      else if (type_.eq.'character') then
         byte_type = 1.0d0
      else 
         byte_type = zero
         call out_%error("Type not recognised in check_mem")
      endif

      ! check if omp_private_array
      size_adj = real(size_, dp)
      do i = 1, size(mem%omp_private_arrays)
         if (trim(name_).eq.trim(mem%omp_private_arrays(i))) then
            size_adj = real(size_, dp) * n_threads_omp_used
            exit
         end if
      end do

      !size in bytes
      size_bytes = size_adj * byte_type

      !check if the array is already registered
      do i = 1, mem%num_arrays
         if (trim(mem%arrays(i)%name_).eq.trim(name_)) then
            if (mem%arrays(i)%is_allocated) &
               call out_%error("The array '"//name_//"' is already allocated.")
         end if
      end do
      !check memory limit
      if (mem%memory_used + size_bytes.gt. &
          mem%memory_available*bytes_per_gb) &
         call out_%error("Exceeded the available memory limit for array '"// &
                          name_//"'.")

   end subroutine check_mem


   !> Subroutine for updating the saved information of memory allocations
   !!    In/Out : mem          -- memory_type
   !!    Input  : name_        -- string identifying the array
   !!    Input  : size_        -- the size of the array [int]
   !!    Input  : type_        -- the type of the array 
   subroutine update_saved_info(mem, name_, size_, type_)

      !input/output variables
      class(memory_manager_type), intent(inout) :: mem
      character(len=*), intent(in) :: name_
      integer, intent(in) :: size_
      character(len=*), intent(in) :: type_

      !internal variables
      logical :: found
      integer :: i
      real(dp) :: size_adj
      real(dp) :: byte_type
      real(dp) :: size_bytes

      if ( type_.eq.'real') then
         byte_type = 8.0d0
      else if (type_.eq.'complex') then
         byte_type = 16.0d0
      else if (type_.eq.'integer') then
         byte_type = 4.0d0
      else if (type_.eq.'character') then
         byte_type = 1.0d0
      else 
         byte_type = zero
         call out_%error("Type not recognised in update_saved_info")
      endif

      ! check if omp_private_array
      size_adj = real(size_, dp)
      do i = 1, size(mem%omp_private_arrays)
         if (trim(name_).eq.trim(mem%omp_private_arrays(i))) then
            size_adj = real(size_, dp) * n_threads_omp_used
            exit
         end if
      end do

      !size in bytes
      size_bytes = size_adj * byte_type

      found = .false.
      !check if the array is already registered
      do i = 1, mem%num_arrays
         if (trim(mem%arrays(i)%name_).eq.trim(name_)) then
            !update existing entry
            mem%arrays(i)%n_elements = size_
            mem%arrays(i)%size_bytes = size_bytes
            mem%arrays(i)%is_allocated = .true.
            found = .true.
            exit
         end if
      end do

      if (.not. found) then
         ! Add a new entry if not found
         if (mem%num_arrays.lt.max_arrays) then
            mem%num_arrays = mem%num_arrays + 1
            mem%arrays(mem%num_arrays)%name_ = name_
            mem%arrays(mem%num_arrays)%n_elements = size_
            mem%arrays(mem%num_arrays)%size_bytes = size_bytes
            mem%arrays(mem%num_arrays)%is_allocated = .true.
         else
            call out_%error("Exceeded the maximum number of trackable arrays &
                            &for array '"//name_//"'.")
            stop
         end if
      end if      

      !update total memory and peak memory usage
      mem%memory_used = mem%memory_used + size_bytes
      !peak memory
      mem%peak_memory_used = max(mem%peak_memory_used,  &
                                 mem%memory_used)

   end subroutine update_saved_info


   !> Subroutine for checking if it is safe to deallocate
   !!    In/Out : mem          -- memory_type
   !!    Input  : name_        -- string identifying the array
   subroutine check_deallocate(mem, name_)

      !input/output variables
      class(memory_manager_type), intent(inout) :: mem
      character(len=*), intent(in) :: name_

      !internal variables
      integer :: i
      logical :: found

      found = .false.

      !check if the array is registered
      do i = 1, mem%num_arrays
         if (trim(mem%arrays(i)%name_).eq.trim(name_)) then
            if (.not. mem%arrays(i)%is_allocated) &
               call out_%error("The array '"//name_//"' is not allocated.")
            !update the registry
            mem%arrays(i)%is_allocated = .false.
            mem%memory_used = mem%memory_used - mem%arrays(i)%size_bytes
            found = .true.
         end if
      end do
      !the array is not found
      if(.not.found) &
         call out_%error("The array '"//name_//"' was not found.")

   end subroutine check_deallocate


   !> Subroutine for deallocating a generic real array 1D
   !!    In/Out : mem          -- memory_type
   !!    In/Out : ptr          -- array to allocate
   !!    Input  : name_        -- string identifying the array
   subroutine deallocate_array_real_1D(mem, ptr, name_)

      !input/output variables
      class(memory_manager_type), intent(inout) :: mem
      real(dp), allocatable, intent(inout) :: ptr(:)
      character(len=*), intent(in) :: name_

      call mem%check_deallocate(name_)

      deallocate(ptr)

   end subroutine deallocate_array_real_1D


   !> Subroutine for deallocating a generic real array 2D
   !!    In/Out : mem          -- memory_type
   !!    In/Out : ptr          -- array to allocate
   !!    Input  : name_        -- string identifying the array
   subroutine deallocate_array_real_2D(mem, ptr, name_)

      !input/output variables
      class(memory_manager_type), intent(inout) :: mem
      real(dp), allocatable, intent(inout) :: ptr(:,:)
      character(len=*), intent(in) :: name_

      call mem%check_deallocate(name_)

      deallocate(ptr)

   end subroutine deallocate_array_real_2D


   !> Subroutine for deallocating a generic real array 3D
   !!    In/Out : mem          -- memory_type
   !!    In/Out : ptr          -- array to allocate
   !!    Input  : name_        -- string identifying the array
   subroutine deallocate_array_real_3D(mem, ptr, name_)

      !input/output variables
      class(memory_manager_type), intent(inout) :: mem
      real(dp), allocatable, intent(inout) :: ptr(:,:,:)
      character(len=*), intent(in) :: name_

      call mem%check_deallocate(name_)

      deallocate(ptr)

   end subroutine deallocate_array_real_3D


   !> Subroutine for deallocating a generic complex array 1D
   !!    In/Out : mem          -- memory_type
   !!    In/Out : ptr          -- array to allocate
   !!    Input  : name_        -- string identifying the array
   subroutine deallocate_array_complex_1D(mem, ptr, name_)

      !input/output variables
      class(memory_manager_type), intent(inout) :: mem
      complex(dp), allocatable, intent(inout) :: ptr(:)
      character(len=*), intent(in) :: name_

      call mem%check_deallocate(name_)

      deallocate(ptr)

   end subroutine deallocate_array_complex_1D


   !> Subroutine for deallocating a generic complex array 1D
   !!    In/Out : mem          -- memory_type
   !!    In/Out : ptr          -- array to allocate
   !!    Input  : name_        -- string identifying the array
   subroutine deallocate_array_complex_2D(mem, ptr, name_)

      !input/output variables
      class(memory_manager_type), intent(inout) :: mem
      complex(dp), allocatable, intent(inout) :: ptr(:,:)
      character(len=*), intent(in) :: name_

      call mem%check_deallocate(name_)

      deallocate(ptr)

   end subroutine deallocate_array_complex_2D


   !> Subroutine for deallocating a generic int array 1D
   !!    In/Out : mem          -- memory_type
   !!    In/Out : ptr          -- array to allocate
   !!    Input  : name_        -- string identifying the array
   subroutine deallocate_array_int_1D(mem, ptr, name_)

      !input/output variables
      class(memory_manager_type), intent(inout) :: mem
      integer, allocatable, intent(inout) :: ptr(:)
      character(len=*), intent(in) :: name_

      call mem%check_deallocate(name_)

      deallocate(ptr)

   end subroutine deallocate_array_int_1D


   !> Subroutine for deallocating a generic int array 1D
   !!    In/Out : mem          -- memory_type
   !!    In/Out : ptr          -- array to allocate
   !!    Input  : name_        -- string identifying the array
   subroutine deallocate_array_int_2D(mem, ptr, name_)

      !input/output variables
      class(memory_manager_type), intent(inout) :: mem
      integer, allocatable, intent(inout) :: ptr(:,:)
      character(len=*), intent(in) :: name_

      call mem%check_deallocate(name_)

      deallocate(ptr)

   end subroutine deallocate_array_int_2D


   !> Subroutine for deallocating a generic char array 1D
   !!    In/Out : mem          -- memory_type
   !!    In/Out : ptr          -- array to allocate
   !!    Input  : name_        -- string identifying the array
   subroutine deallocate_array_char_1D(mem, ptr, name_)

      !input/output variables
      class(memory_manager_type), intent(inout) :: mem
      character(len=*), allocatable, intent(inout) :: ptr(:)
      character(len=*), intent(in) :: name_

      call mem%check_deallocate(name_)

      deallocate(ptr)

   end subroutine deallocate_array_char_1D


   !> Subroutine for deallocating a generic char array 1D
   !!    In/Out : mem          -- memory_type
   !!    In/Out : ptr          -- array to allocate
   !!    Input  : name_        -- string identifying the array
   subroutine deallocate_array_char_2D(mem, ptr, name_)

      !input/output variables
      class(memory_manager_type), intent(inout) :: mem
      character(len=*), allocatable, intent(inout) :: ptr(:,:)
      character(len=*), intent(in) :: name_

      call mem%check_deallocate(name_)

      deallocate(ptr)

   end subroutine deallocate_array_char_2D


   !> Subroutine for printing the memory status                                   
   !!    Input  : mem          -- memory_type      
   subroutine print_status(mem)             
                                                   
      !input/output variables                      
      class(memory_manager_type), intent(in) :: mem
                                                   
      !internal variables                          
      integer :: i                                 
      real(dp) :: size_in_unit
      real(dp) :: total_memory, peak_memory
      character(len=10) :: unit_total
      character(len=10) :: unit_peak
      character(len=10) :: unit_
      character(len=40) :: array_name      

      ! Determine the unit for total memory used
      if (mem%memory_used.lt.1024.0d0) then
         total_memory = real(mem%memory_used, dp)
         unit_total = " B"
      else if (mem%memory_used.lt.1024.0d0**2) then
         total_memory = real(mem%memory_used, dp) / 1024.0d0
         unit_total = " KB"
      else if (mem%memory_used.lt.1024.0d0**3) then
         total_memory = real(mem%memory_used, dp) / (1024.0d0**2)
         unit_total = " MB"
      else
         total_memory = real(mem%memory_used, dp) / (1024.0d0**3)
         unit_total = " GB"
      end if

      ! Determine the unit for peak memory used
      if (mem%peak_memory_used.lt.1024.0d0) then
         peak_memory = mem%peak_memory_used 
         unit_peak = " B"
      else if (mem%peak_memory_used.lt.1024.0d0**2) then
         peak_memory = mem%peak_memory_used / 1024.0d0
         unit_peak = " KB"
      else if (mem%peak_memory_used.lt.1024.0d0**3) then
         peak_memory = mem%peak_memory_used / 1024.0d0**2
         unit_peak = " MB"
      else
         peak_memory = mem%peak_memory_used / 1024.0d0**3
         unit_peak = " GB"
      end if

      write(out_%iunit, '(a,f10.3,a)') " Peak memory used  : ", peak_memory,  &
                                       unit_peak
      write(out_%iunit, out_%sticks)
      if(.not.is_equal(total_memory,zero)) then
         write(out_%iunit, '(a,f10.3,a)') " Current memory used : ", total_memory, &
                                            unit_total
      endif

      if(any(mem%arrays(:)%is_allocated)) &
         write(out_%iunit,'(a)') " Allocated arrays:"                                 
      do i = 1, mem%num_arrays
         if (mem%arrays(i)%is_allocated) then
            ! determine the unit (kb, mb, gb) based on size
            if (mem%arrays(i)%size_bytes.lt.1024.0d0) then
               size_in_unit = mem%arrays(i)%size_bytes
               unit_ = " B"
            else if (mem%arrays(i)%size_bytes.lt.1024.0d0**2) then
               size_in_unit = mem%arrays(i)%size_bytes / 1024.0d0
               unit_ = " KB"
            else if (mem%arrays(i)%size_bytes.lt.1024.0d0**3) then
               size_in_unit = mem%arrays(i)%size_bytes / (1024.0d0**2)
               unit_ = " MB"
            else
               size_in_unit = mem%arrays(i)%size_bytes / (1024.0d0**3)
               unit_ = " GB"
            end if

            ! Format the array name to fit into 40 characters, left-aligned
            write(array_name, '(a30)') trim(mem%arrays(i)%name_)
            ! Print the formatted information
            write(out_%iunit, '(a,f10.3,a)') ' - '//adjustl(trim(array_name))//&
                                             ": ", size_in_unit, trim(unit_)
         end if
      end do      
                                                                                  
   end subroutine print_status      

end module memory_manager_module
