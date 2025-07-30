!> Input module
!!
!! This module contains the subroutines for the reading the input
!!
!! Date         : 2025
!!
module input_module
      
   !$ use omp_lib
   use output_module
   use parameters_module
   use string_manipulation_module
   use array_manipulation_module
   use target_module
   use bem_module
   use fq_module
   use fqfmu_module
   use wfq_module
   use wfqfmu_module
   use algorithm_module
   use control_module
   use memory_manager_module

   implicit none

   !public variables
   public inp_

   type inp_type
      integer  :: iunit = 10 ! unit of the input file
      character (len=200) :: filename
      logical  :: used_defaults = .false.
      logical  :: yaml_read = .false.
      logical  :: xyz_read  = .false.
      character (len=200) :: xyz_file

      contains

      !get/check the input file
      procedure :: get_arguments
      procedure :: check_input_file
      !read the input file
      procedure :: read_general_settings
      procedure :: read_
      procedure :: get_number_lines_keywords
      procedure :: get_target
      procedure :: get_section_atomtypes
      procedure :: get_section_parameters
      procedure :: get_section_geometry
      !printing
      procedure :: print_input_info
   end type inp_type
    
   type (inp_type), Save :: inp_

   contains

!-------------------------------------------------------------------------------
!  SUBROUTINES FOR GETTING/CHECKING THE INPUT FILE
!-------------------------------------------------------------------------------

   !> Subroutine for getting the input arguments
   !!    In/Out : inp_      -- inp_type 
   subroutine get_arguments(inp_)
     
      implicit none
  
      !input/output variables
      class(inp_type), intent(inout) :: inp_

      !internal variables
      integer :: narg
       
      narg=command_argument_count()
      if (narg .eq. 1) then 
         call get_command_argument(narg,inp_%filename)
      else if (narg .gt. 1) then !more than one variables detected
         write(*,'(/a/)') "More than one input variables to Fortran file"
         stop
      else if (narg .eq. 0) then 
         write(*,'(/a/)') "Type the input file (e.g. filename.tmp)"
         read(*,*) inp_%filename(1:99)
      endIf
      inp_%filename = trim(inp_%filename)
    
   end subroutine get_arguments


   !> Subroutine for checking the existence of the input file
   !!    In/Out : inp_      -- inp_type 
   subroutine check_input_file(inp_)
     
      implicit none
  
      !input/output variables
      class(inp_type)  :: inp_
  
      !internal variables
      logical          :: exists
      
      inquire(file=inp_%filename,exist=exists)
      
      if(.not.exists) &
          write(*, '(a)') "file "//trim(inp_%filename)//" does not exist.&
          & Check the name."
      
   end subroutine check_input_file

!-------------------------------------------------------------------------------
!  SUBROUTINES FOR READING THE INPUT FILE
!-------------------------------------------------------------------------------

   !> Subroutine for reading the general settings from the input file
   !! In the first three lines of the input file, we have:
   !! 1) Output file name
   !! 2) OMP Threads
   !! 3) Available memory for calculation
   !!    In/Out : inp_      -- inp_type 
   subroutine read_general_settings(inp_)
     
      implicit none
  
      !input/output variables
      class(inp_type), intent(inout) :: inp_

      !internal variables
      integer :: iost
      character(len=200) :: line
 
      inp_%iunit = 10! file
      open(unit=inp_%iunit,file=trim(inp_%filename),status="old",iostat=iost,&
           action="read") 
         if (iost.ne.0) then
            write(*,'(a)') 'Error in opening '//trim(inp_%filename)//' file'
            stop
         endif
         rewind(inp_%iunit) 

         !read the first line
         read(inp_%iunit, '(a)') line
         call validate_and_parse_char(line,"output file",out_%filename)
         !read the second line
         read(inp_%iunit, '(a)') line
         call validate_and_parse_int(line,"omp threads",n_threads_omp)
         n_threads_omp_used = n_threads_omp
         !read the third line
         read(inp_%iunit, '(a)') line
         call validate_and_parse_real(line,"memory", mem_man%memory_available)
      close(inp_%iunit) 
      
   end subroutine read_general_settings


   !> Subroutine for reading the input file
   !!    In/Out : inp_      -- inp_type 
   subroutine read_(inp_)
    
      implicit none
  
      !input/output variables
      class(inp_type), intent(inout) :: inp_
  
      !internal variables
      integer :: iost
      integer :: n_lines_keywords
      integer :: i
      character(len=200) :: line, keyword, value_

      open(unit=inp_%iunit,file=trim(inp_%filename),status="old",iostat=iost,&
           action="read") 
         if (iost.ne.0) call out_%error('Error in opening '// &
                                        trim(inp_%filename)//' file')
         !restart reading
         rewind(inp_%iunit) 

         n_lines_keywords = inp_%get_number_lines_keywords()
         call inp_%get_target(n_lines_keywords)

         !read and process each line for the keyword sections
         rewind(inp_%iunit) 
         do i = 1, n_lines_keywords
            read(inp_%iunit,'(a)',iostat=iost) line
            if (iost.ne.0) exit  ! end of file
            
            !cycle if commented line in input file
            if(commented_line(line)) cycle

            !divide line in keyword and value
            call get_keyword_and_value(out_%iunit,line,keyword,value_)

            !determine the section from the first word of the keyword
            select case (nth_word_in_string(keyword, 1)) 
            case ("omp", "memory", "forcefield","bem", "parameters")
               cycle !already read or not to be read here
            case ("what")
               call get_section_what(value_)
            case ("algorithm")
               call get_section_algorithm(line)
            case ("field")
               call get_section_field(line)
            case ("control")
               call get_section_control(line)
            case ("used") 
               call validate_and_parse_bool(line, &
                                 "used defaults for atom_types or parameters", &
                                 inp_%used_defaults)
            case ("number")
               if (trim(keyword).eq.'number of atoms') then
                  call validate_and_parse_int(line,"number of atoms", &
                                                   target_%n_atoms)
               else if (trim(keyword).eq.'number of molecules') then
                  call validate_and_parse_int(line,"number of molecules", &
                                                   target_%n_mol)
               endif
            case ("atom_types")
               if (trim(keyword).eq.'atom_types number') &
                  call validate_and_parse_int(line,"atom_types number", &
                                                   target_%n_atomtypes)
            case ("output")
               call get_section_output(line)
            case default
               call out_%error("Error: Unknown section for keyword: "// &
                                trim(keyword))
            end select
         end do

         if(field%dynamic) call field%set_frequencies()

         if(target_%n_atomtypes.gt.0) then
            call inp_%get_section_atomtypes(n_lines_keywords)
            call inp_%get_section_parameters(n_lines_keywords)
            call inp_%get_section_geometry(n_lines_keywords)
            if(target_%name_(1:3).eq.'wfq') then
               call target_%assign_model_parameters()
               if(target_%name_.eq.'wfqfmu') call read_alpha()
            endif
            !scaling array coordinates
            call array_scale(tobohr,3*target_%n_atoms,target_%coord)
         endif

         if(bem%exists) then 
            call bem%read_gmsh_file()
            call bem%read_permittivity_bem()
            call bem%create_tesserae()
         endif

         !assign nvar and other dimensions
         call target_%assign_model_dimensions()

         !reassign tolerance for RMSE
         algorithm%threshold = algorithm%threshold*sqrt(dble(target_%n_var)) 

         !set correct n_omp_thread_used for memory checks
         if(algorithm%matrix_in_parallel) n_threads_omp_used = 1

      close(inp_%iunit) 
  
   end subroutine read_


   !> Function to get the number of lines containing keywords
   !!    Input  : inp__      -- inp_ type
   integer function get_number_lines_keywords(inp_) result(n_lines)
  
      implicit none
  
      !input/output variables
      class(inp_type), intent(in) :: inp_
  
      !internal variables
      integer :: io
      character (len=200) :: line
        
      n_lines = 0
      do
        read(inp_%iunit,*,iostat=io) line
        if (io.ne.0) exit
        if (trim(nth_word_in_string(line, 1)).eq.'input_geometry:') exit
        n_lines = n_lines + 1
      end do
  
   end function get_number_lines_keywords


   !> Subroutine for getting the target
   !!    Input  : inp__      -- inp_ type
   !!    Input  : nlines     -- number of lines to read
   subroutine get_target(inp_,nlines)
  
      implicit none
  
      !input/output variables
      class(inp_type), intent(in) :: inp_
      integer, intent(in)         :: nlines
  
      !internal variables
      integer :: i 
      integer :: iost
      character(len=200) :: line
      character(len=200) :: keyword
      character(len=200) :: value_
      character(len=200) :: target_name = ""
      character(len=200) :: static_ff 
      character(len=200) :: dynamic_ff
      character(len=200) :: kernel
  
      rewind(inp_%iunit)
   
      do i = 1, nlines
         read(inp_%iunit,'(a)',iostat=iost) line

         !cycle if commented line in input file
         if(commented_line(line)) cycle

         !divide line in keyword and value
         call get_keyword_and_value(out_%iunit,line,keyword,value_)

         !determine the section from the first word of the keyword
         select case (nth_word_in_string(keyword, 1)) 
         case ("forcefield")
            call get_section_forcefield(line, static_ff, dynamic_ff, kernel)
         case ("bem")
            call get_section_bem(line)
         end select
      enddo

      if (trim(bem%gmsh_file).ne.'none') then 
         bem%exists = .true.
      else 
         if(trim(dynamic_ff).ne.'none') then
            target_name = trim(dynamic_ff)
         else
            target_name = trim(static_ff)
         endif
      endif
         
      call assign_target(target_name)

      target_%kernel = trim(kernel)
      target_%forcefield = trim(static_ff)
  
   end subroutine get_target


   !> Subroutine for assigning the target
   !!    Input  : target_name    -- name assigned to the target
   subroutine assign_target(target_name)
 
      implicit none
 
      !input/output variables
      character(len=*) :: target_name 
      
      if(target_name.eq.'fq') then         
         target_ => fq
         target_%name_ = "fq"
      else if(target_name.eq.'fqfmu') then         
         target_ => fqfmu
         target_%name_ = "fqfmu"
      else if(target_name.eq.'wfq') then 
         target_ => wfq
         target_%name_ = "wfq"
      else if(target_name.eq.'wfqfmu') then 
         target_ => wfqfmu
         target_%name_ = "wfqfmu"
      else if (target_name.eq."bem") then
            target_ => bem
            target_%name_ = "bem"
      else if (target_name.eq."") then
         if(bem%exists) then
            target_ => bem
            target_%name_ = "bem"
         endif
      else 
         call out_%error( "No target model defined")
      endif
   
   end subroutine assign_target


   !> Subroutine for getting section what 
   !!    Input  : value_      -- value_ of what keyword
   subroutine get_section_what(value_)
 
      implicit none
 
      !input/output variables
      character(len=200), intent(in) :: value_
 
      !internal variables
      integer :: len_line
      integer :: i
      integer :: index_
      integer :: start
      character(len=200) :: trimmed_line
 
      trimmed_line = trim(value_)
      len_line = len(trimmed_line)
      !just one keyword
      if (index(trimmed_line, ",").eq.0) then
         out_%what = trimmed_line
      else 
         !two keywords 
         start  = 1
         index_ = 0
         do i = 1, len_line
            if (trimmed_line(i:i).eq.",".or.i.eq.len_line) then
               index_ = index_ + 1 
               if(index_.eq.1) out_%what=trim(adjustl(trimmed_line(start:i-1)))
               if(index_.eq.2) control%restart = .true.
               start = i + 1
            endif
         end do
      endif
 
   end subroutine get_section_what


   !> Subroutine for getting section algorithm
   !!    Input  : line       -- line for algorithm section
   subroutine get_section_algorithm(line)
 
      implicit none
 
      !input/output variables
      character(len=200), intent(in) :: line

      !internal variables
      character(len=200) :: keyword
      character(len=200) :: value_
      character(len=200) :: method
      character(len=200) :: parallel_execution
      character(len=200) :: adaptive_tuning

      call get_keyword_and_value(out_%iunit,line,keyword,value_)

      select case (trim(keyword))
      case ("algorithm method")
          call validate_and_parse_char(line, "algorithm method", method)
          if (trim(method).eq."inversion") then 
             algorithm%inversion = .true.
          else if (trim(method).eq."iterative") then 
             algorithm%iterative = .true.
          else if (trim(method).eq."iterative on the fly") then 
             algorithm%iterative = .true.
             algorithm%on_the_fly = .true.
          endif
      case ("algorithm parallel execution")
          call validate_and_parse_char(line, "algorithm parallel execution", &
                                       parallel_execution)
          if (trim(parallel_execution).eq.'frequencies') then 
             algorithm%freq_in_parallel = .true.
          else if (trim(parallel_execution).eq.'matrix') then 
             algorithm%matrix_in_parallel = .true.
          endif
      case ("algorithm number of iterations")
          call validate_and_parse_int(line, "algorithm number of iterations", &
                                      algorithm%n_iter)
      case ("algorithm gmres dimension")
          call validate_and_parse_int(line, "algorithm gmres dimension", &
                                      algorithm%n_dim_gmres)
      case ("algorithm tolerance")
          call validate_and_parse_real(line, "algorithm tolerance", &
                                       algorithm%threshold)
      case ("algorithm adaptive tuning")
          call validate_and_parse_char(line, "algorithm adaptive tuning", &
                                       adaptive_tuning)
          if(trim(adaptive_tuning).eq.'yes') then
             algorithm%adaptive_tuning = .true.
          else if(trim(adaptive_tuning).eq.'no') then
             algorithm%adaptive_tuning = .false.
          endif
      end select

   end subroutine get_section_algorithm


   !> Subroutine for getting section forcefield
   !!    Input  : line       -- line for forcefield
   !!    Output : static_ff  -- static forcefield
   !!    Output : dynamic_ff -- dynamic forcefield
   !!    Output : kernel     -- kernel
   subroutine get_section_forcefield(line, static_ff, dynamic_ff, kernel)
   
      implicit none
   
      !input/output variables
      character(len=200), intent(in) :: line
      character(len=200), intent(out) :: static_ff
      character(len=200), intent(out) :: dynamic_ff
      character(len=200), intent(out) :: kernel
  
      !internal variables
      character(len=200) :: keyword
      character(len=200) :: value_

      call get_keyword_and_value(out_%iunit,line,keyword,value_)

      select case (trim(keyword))
      case ("forcefield static")
          call validate_and_parse_char(line, "forcefield static", &
                                       static_ff)
      case ("forcefield dynamic")
          call validate_and_parse_char(line, "forcefield dynamic", &
                                       dynamic_ff)
      case ("forcefield kernel")
          call validate_and_parse_char(line, "forcefield kernel", &
                                       kernel)
      end select
  
   end subroutine get_section_forcefield


   !> Subroutine for getting section field
   !!    Input  : line       -- line for field section
   subroutine get_section_field(line)
  
      implicit none
  
      !input/output variables
      character(len=200), intent(in) :: line
  
      !internal variables
      character(len=200) :: keyword
      character(len=200) :: value_
      character(len=200) :: type_field

      call get_keyword_and_value(out_%iunit,line,keyword,value_)

      select case (trim(keyword))
      case ("field type")
         call validate_and_parse_char(line,"field type", type_field)
         if(trim(type_field).eq.'static') then
            field%static = .true.
         else if(trim(type_field).eq.'dynamic') then
            field%dynamic = .true.
         endif
      case ("field rhs type")
         call validate_and_parse_char(line,"field rhs type", field%rhs_form)
      case ("field field intensity")
         call validate_and_parse_real(line,"field field intensity", field%e_0)
      case ("field nfreq")
         call validate_and_parse_int(line,"field nfreq", field%n_freq)
         if(field%dynamic) &
            call mem_man%alloc(field%freq, field%n_freq, 'field%freq')
         if(field%n_freq.gt.0.and.field%dynamic &
            .and.field%n_freq.lt.n_threads_omp) &
            n_threads_omp_used = field%n_freq
      case ("field min freq")
         call validate_and_parse_real(line,"field min freq", field%min_freq)
      case ("field max freq")
         call validate_and_parse_real(line,"field max freq", field%max_freq)
      case ("field step freq")
         call validate_and_parse_real(line,"field step freq", field%step_freq)
      case ("field external freq")
         if(field%n_freq.gt.0.and.allocated(field%freq).and.&
            field%min_freq.eq.zero.and.field%max_freq.eq.zero) then
            read(value_,*) field%freq(1:field%n_freq)
         endif
      case ("field polarization")
         call validate_and_parse_char(line,"field polarization", &
                                      field%polarization)
         call field%set_polarization_variables()
      end select
  
   end subroutine get_section_field


   !> Subroutine for getting section control
   !!    Input  : line       -- line for control section
   subroutine get_section_control(line)
  
      implicit none
  
      !input/output variables
      character(len=200), intent(in) :: line
  
      !internal variables
      character(len=200) :: keyword
      character(len=200) :: value_
      logical :: no_save_info = .false.

      call get_keyword_and_value(out_%iunit,line,keyword,value_)

      select case (trim(keyword))
      case ("control no info file")
         call validate_and_parse_bool(line, "control no info file", &
                                      no_save_info)
         control%save_info = .not. no_save_info
         if (control%save_info) then
         write(out_%info_file,'(a)') out_%filename(1:len_trim(out_%filename)-4)&
                                     //'.info'
         write(out_%file_tar,'(a)') out_%filename(1:len_trim(out_%filename)-4) &
                                    //'.tar.gz'
         endif
      case ("control principal axes")
         call validate_and_parse_bool(line, "control principal axes", &
                                      control%principal_axis)
      end select
  
   end subroutine get_section_control


   !> Subroutine for getting section output
   !!    Input  : line       -- line for output section
   subroutine get_section_output(line)
  
      implicit none
  
      !input/output variables
      character(len=200), intent(in) :: line
  
      !internal variables
      character(len=200) :: keyword
      character(len=200) :: value_

      call get_keyword_and_value(out_%iunit,line,keyword,value_)

      select case (trim(keyword))
      case ("output verbose")
         call validate_and_parse_int(line, "output verbose", out_%ivrb)
      case ("output maxima analysis")
         call validate_and_parse_char(line, "output maxima analysis", &
                                      control%maxima_analysis)
      end select
  
   end subroutine get_section_output


   !> Subroutine for getting section bem
   !!    Input  : line       -- line for bem section
   subroutine get_section_bem(line)
  
      implicit none
  
      !input/output variables
      character(len=200), intent(in) :: line
  
      !internal variables
      character(len=200) :: keyword
      character(len=200) :: value_

      call get_keyword_and_value(out_%iunit,line,keyword,value_)

      select case (trim(keyword))
      case ("bem mesh file")
         call validate_and_parse_char(line, "bem mesh file", bem%gmsh_file)
      case ("bem normal scalar factor")
         call validate_and_parse_real(line, "bem normal scalar factor", &
                                      bem%normal_factor)
      case ("bem permittivity file")
         call validate_and_parse_char(line, "bem permittivity file", &
                                      bem%permittivity_file)
      case ("bem permittivity")
         call validate_and_parse_char(line, "bem permittivity", &
                                      bem%permittivity_type)
      case ("bem green function")
         call validate_and_parse_char(line, "bem green function", &
                                      bem%green_function)
      case ("bem sphere radius")
         call validate_and_parse_real(line, "bem sphere radius", bem%sphere_r)
      case ("bem solvent")
         call validate_and_parse_char(line, "bem solvent", bem%solvent)
      case ("bem epsilon solvent")
         call validate_and_parse_real(line, "bem epsilon solvent", &
                                      bem%epsilon_solvent)
      case ("bem charge constraint")
         call validate_and_parse_bool(line, "bem charge constraint", &
                                      bem%charge_constraint)
      case ("bem variant")
         call validate_and_parse_char(line, "bem variant", &
                                      bem%variant)
      end select 
  
   end subroutine get_section_bem


   !> Subroutine for getting section atomtypes
   !!    Input  : inp_       -- inp_ type
   !!    Input  : nlines     -- number of lines to read
   subroutine get_section_atomtypes(inp_,nlines)
  
      implicit none
  
      !input/output variables
      class(inp_type), intent(in) :: inp_
      integer, intent(in)         :: nlines
  
      !internal variables
      character(len=200) :: keyword
      character(len=200) :: value_
      character(len=200) :: line
      integer :: i
      integer :: atom_index

      if(target_%n_atomtypes.gt.1) target_%heterogeneous = .true.
      call mem_man%alloc(target_%atom_type, target_%n_atomtypes, &
                         'target_%atom_type')
      call mem_man%alloc(target_%chi, target_%n_atomtypes, 'target_%chi')
      call mem_man%alloc(target_%eta, target_%n_atomtypes, 'target_%eta')
      call mem_man%alloc(target_%alpha, target_%n_atomtypes, 'target_%alpha')
      call mem_man%alloc(target_%r_q, target_%n_atomtypes, 'target_%r_q')
      call mem_man%alloc(target_%r_mu, target_%n_atomtypes, 'target_%r_mu')

      rewind(inp_%iunit)
      !read data for each atom type
      atom_index = 0
      do i = 1, nlines
         read(inp_%iunit,'(a)') line

         if(commented_line(line)) cycle

         ! Parse keyword and value
         call get_keyword_and_value(out_%iunit,line,keyword,value_)
         if(nth_word_in_string(keyword, 1).eq.'atom_types') then
            select case (trim(keyword))
            case ("atom_types number")
               cycle
            case ("atom_types name")
               atom_index = atom_index + 1
               call validate_and_parse_char(line, "atom_types name", &
                                            target_%atom_type(atom_index))
            case ("atom_types chi")
               call validate_and_parse_real(line, "atom_types chi", &
                                            target_%chi(atom_index))
            case ("atom_types eta")
               call validate_and_parse_real(line, "atom_types eta", &
                                            target_%eta(atom_index))
            case ("atom_types alpha")
               call validate_and_parse_real(line, "atom_types alpha", &
                                            target_%alpha(atom_index))
            case ("atom_types rq")
               call validate_and_parse_real(line, "atom_types rq", &
                                            target_%r_q(atom_index))
            case ("atom_types rmu")
               call validate_and_parse_real(line, "atom_types rmu", &
                                            target_%r_mu(atom_index))
            case default
               call out_%error("Error: Unknown keyword " // trim(keyword))
            end select
         end if
      end do      

   end subroutine get_section_atomtypes


   !> Subroutine for getting section atomtypes
   !!    Input  : inp_       -- inp_ type
   !!    Input  : nlines     -- number of lines to read
   subroutine get_section_parameters(inp_,nlines)
  
      implicit none
  
      !input/output variables
      class(inp_type)  :: inp_
      integer, intent(in)             :: nlines
       
      !internal variables
      character(len=200) :: keyword
      character(len=200) :: value_
      character(len=200) :: line
      character(len=200) :: atom_type_param
      character(len=200) :: interaction_name
      integer :: i
      integer :: atom_index
      integer :: atom_index_1 !for interaction
      integer :: atom_index_2 !for interaction
  
      !these are all zero if static or energy calculation

      call mem_man%alloc(parameters%tau,     target_%n_atomtypes, &
                         "parameters%tau")
      call mem_man%alloc(parameters%sigma_0, target_%n_atomtypes, &
                         "parameters%sigma0")
      call mem_man%alloc(parameters%scaling, target_%n_atomtypes, &
                         "parameters%scaling")
      call mem_man%alloc(parameters%A_ij,    target_%n_atomtypes, &
                         "parameters%A_ij")
      call mem_man%alloc(parameters%fermi_d, target_%n_atomtypes, &
                         target_%n_atomtypes, "parameters%fermi_d")
      call mem_man%alloc(parameters%fermi_s, target_%n_atomtypes, &
                         target_%n_atomtypes, "parameters%fermi_s")
      call mem_man%alloc(parameters%fermi_energy, target_%n_atomtypes, &
                         "parameters%fermi_energy")
      call mem_man%alloc(parameters%permittivity_type, target_%n_atomtypes, &
                         "parameters%permittivity_type")
      call mem_man%alloc(parameters%wfqfmu_file,   target_%n_atomtypes, &
                         "parameters%wfqfmu_file")

      atom_index   = 0
      atom_index_1 = 0
      atom_index_2 = 0
  
      rewind(inp_%iunit)
      !read data for each atom type
      do i = 1, nlines
         read(inp_%iunit,'(a)') line
         if(commented_line(line)) cycle

         ! Parse keyword and value
         call get_keyword_and_value(out_%iunit,line,keyword,value_)
         if(nth_word_in_string(keyword, 1).eq.'parameters') then
            select case (trim(keyword))
            case ("parameters atomtype name")
               call validate_and_parse_char(line, "parameters atomtype name", &
                                            atom_type_param)
               atom_index = get_index_atom_type(atom_type_param)
            case ("parameters atomtype tau")
               call validate_and_parse_real(line, "parameters atomtype tau", &
                                            parameters%tau(atom_index))
            case ("parameters atomtype sigma0")
               call validate_and_parse_real(line, "parameters atomtype sigma0",&
                                            parameters%sigma_0(atom_index))
            case ("parameters atomtype a_ij")
               call validate_and_parse_real(line, "parameters atomtype a_ij", &
                                            parameters%a_ij(atom_index))
            case ("parameters atomtype fermi function d")
               call validate_and_parse_real(line, &
                                       "parameters atomtype fermi function d", &
                                      parameters%fermi_d(atom_index,atom_index))
            case ("parameters atomtype fermi function s")
               call validate_and_parse_real(line, &
                                       "parameters atomtype fermi function s", &
                                      parameters%fermi_s(atom_index,atom_index))
            case ("parameters atomtype fermi energy")
               call validate_and_parse_real(line, &
                                       "parameters atomtype fermi energy", &
                                       parameters%fermi_energy(atom_index))
            case ("parameters atomtype scaling sigma0-tau")
               call validate_and_parse_real(line, &
                                      "parameters atomtype scaling sigma0-tau",&
                                       parameters%scaling(atom_index))
            case ("parameters atomtype wfqfmu file")
               call validate_and_parse_char(line, &
                                            "parameters atomtype wfqfmu file", &
                                            parameters%wfqfmu_file(atom_index))
            case ("parameters atomtype permittivity")
               call validate_and_parse_char(line, &
                                            "parameters atomtype permittivity",&
                                            parameters%permittivity_type(atom_index))
            case ("parameters interaction name")
               call validate_and_parse_char(line,"parameters interaction name",&
                                            interaction_name)
               call get_indeces_atom_type_interaction(interaction_name, &
                                                      atom_index_1,     &
                                                      atom_index_2)
            case ("parameters interaction fermi function d")
               call validate_and_parse_real(line, &
                                    "parameters interaction fermi function d", &
                                  parameters%fermi_d(atom_index_1,atom_index_2))
            case ("parameters interaction fermi function s")
               call validate_and_parse_real(line, &
                                    "parameters interaction fermi function s", &
                                  parameters%fermi_s(atom_index_1,atom_index_2))
            case default
               call out_%error("Error: Unknown keyword " // trim(keyword))
            end select
         end if
      end do      

   end subroutine get_section_parameters


   !> Function to get the atom index given the atom type
   !!    Input  : atomtype      -- atomtype character
   integer function get_index_atom_type(atomtype) result(index_)

      implicit none
  
      !input/output variables
      character(len=*), intent(in) :: atomtype

      !internal variables
      integer :: i
      logical :: found
  
      found = .false.
      index_ = 0

      do i = 1, size(target_%atom_type)
          if (trim(target_%atom_type(i)) == trim(atomtype)) then
              index_ = i
              found = .true.
              exit
          endif
      end do
      
      if(found) then
         return
      else
         call out_%error("Atomtype: "//trim(atomtype)//" not recognised")
      endif 
  
   end function get_index_atom_type


   !> Subroutine to get the atom indeces given the atom type interaction
   !!    Input  : interaction   -- interaction name "Ag->Au", ...
   !!    Output : index_1       -- atomtype 1 index
   !!    Output : index_2       -- atomtype 2 index
   subroutine get_indeces_atom_type_interaction(interaction, index_1, index_2) 

      implicit none
  
      !input/output variables
      character(len=*) :: interaction
      integer, intent(out) :: index_1
      integer, intent(out) :: index_2

      !internal variables
      character(len=200) :: atomtype_1
      character(len=200) :: atomtype_2
      integer :: iend
  
      !format line: atomtype_1->atomtype_2
      iend = move_to_character(interaction, '>') - 2
      atomtype_1 = trim(interaction(1:iend))
      atomtype_2 = trim(interaction(iend+3:)) !->skip
      index_1 = get_index_atom_type(atomtype_1)
      index_2 = get_index_atom_type(atomtype_2)
  
   end subroutine get_indeces_atom_type_interaction


   !> Subroutine for getting section geometry and assign map atoms -> atomtypes
   !!    In/Out : inp_       -- inp_ type
   !!    Input  : nlines     -- number of lines to skip (keywords)
   subroutine get_section_geometry(inp_,nlines)
   
      implicit none
   
      !input/output variables
      class(inp_type), intent(inout) :: inp_
      integer, intent(in) :: nlines
   
      !internal variables
      integer :: i
      character(len=200) :: line
      character(len=200) :: keyword
      character(len=200) :: value_
      character(len=200), dimension(:), allocatable :: atomname
      character(len=24) :: format_geom = "(a6,1x,i10,3(1x,f25.16))"
   
      call mem_man%alloc(atomname, target_%n_atoms, "atomname")
      call mem_man%alloc(target_%i_mol, target_%n_atoms, "target_%i_mol")
      call mem_man%alloc(target_%coord, 3, target_%n_atoms, "target_%coord")
   
      rewind(inp_%iunit)
      !skip nlines
      do i = 1, nlines
         read (inp_%iunit, '(a)') line 
      enddo
      read(inp_%iunit, '(a)') line
      call get_keyword_and_value(out_%iunit,line,keyword,value_)
      if(trim(value_).eq.'read from yaml file') then
         inp_%yaml_read = .true.
      else !xyz file
         inp_%xyz_read = .true.
         inp_%xyz_file = trim(get_xyz_file_from_value(trim(value_)))
      endif

      do i = 1, target_%n_atoms
        read(inp_%iunit,'(a)') line
        !format : atomname 6char, imol i10, coords f25.16
        read(line, format_geom) atomname(i), target_%i_mol(i), &
                                target_%coord(1,i),            &
                                target_%coord(2,i),            &
                                target_%coord(3,i)
      end do

      !create map atomtype - atoms
      call mem_man%alloc(target_%atomic_number, target_%n_atoms, &
                         "target_%atomic_numbers")
      call mem_man%alloc(target_%map_atomtypes, target_%n_atoms, &
                         "target_%map_atomtypes")
      call target_%assign_atomic_numbers(atomname)
      call target_%assign_map_atomtypes(atomname)
      call target_%get_neighbours()
      call create_n_atoms_per_molecule()
      !deallocation
      call mem_man%dealloc(atomname, "atomname")
   
   end subroutine get_section_geometry


   !> Function to get the xyz file from value keyword input geometry
   !!    Input  : value_      -- string
   function get_xyz_file_from_value(value_) result(xyz_file)

      implicit none
  
      !input/output variables
      character(len=*) :: value_
      character(len=200) :: xyz_file

      !internal variables
      integer :: iend

      !format line: read from xyz_file
      iend = index(value_, 'from')
      xyz_file = trim(value_(iend+5:))
  
   end function get_xyz_file_from_value


   !> Subroutine for creating the array: how many atoms per molecules
   subroutine create_n_atoms_per_molecule()

      !internal variables
      integer :: i, k

      call mem_man%alloc(target_%n_atoms_per_molecule, target_%n_mol, &
                         "target_%n_atoms_per_molecule")

      if (target_%n_mol.eq.1) then 
         target_%n_atoms_per_molecule = target_%n_atoms
         return
      endif

      !at least one atom per molecule
      target_%n_atoms_per_molecule = 1
      k = 1
      do i = 1, target_%n_atoms - 1
        if(target_%i_mol(i).eq.target_%i_mol(i+1)) then
           target_%n_atoms_per_molecule(k) = target_%n_atoms_per_molecule(k) + 1
        else
           k = k + 1
        endif
      enddo

   end subroutine create_n_atoms_per_molecule

!-------------------------------------------------------------------------------
!  SUBROUTINES FOR PRINTING
!-------------------------------------------------------------------------------

   !> Subroutine for getting section geometry and assign map atoms -> atomtypes
   !!    In/Out : inp_       -- inp_ type
   subroutine print_input_info(inp_)
  
      implicit none
  
      !input/output variables
      class(inp_type), intent(in) :: inp_
  
      write(out_%iunit,'(1x,a)') "Input File         : "// &
         trim(inp_%filename(1:len_trim(inp_%filename)-4))//".yaml"
      write(out_%iunit,'(1x,a)') "Output File        : "//trim(out_%filename)
      write(out_%iunit,'(1x,a,i3)') "OMP Threads        : ", &
                                       n_threads_omp
      write(out_%iunit,'(1x,a,f7.2)') "Memory (GB)        : ", &
                                       mem_man%memory_available
      write(out_%iunit,out_%sticks) 
      flush(out_%iunit)

      !control & output section
      call control%print_info()
  
      !field informations
      call field%print_info()
       
      !algorithm informations
      call algorithm%print_info()
  
      !Target informations
      call target_%print_info(inp_%used_defaults)
  
      if(bem%exists) call bem%print_parameters()
      if(target_%name_.ne.'bem') call target_%print_coord("Input")
      if(bem%exists) call bem%print_coord_bem("Input")
      flush(out_%iunit)
   
   end subroutine print_input_info


   !> Subroutine to initialize the .info file for post process analysis
   subroutine initialize_save_info()
  
      implicit none
  
      !internal variables 
      integer            :: iost
      integer            :: i
      real(dp)           :: freqeV
  
      If(out_%ivrb.ge.1) then
         write(out_%iunit,'(1x,a)') "Created file info : " // trim(out_%info_file)
         write(out_%iunit,out_%sticks) 
      endif
       
      open(unit=out_%unit_info,file=trim(out_%info_file),status="unknown",iostat=iost)
  
         !write info in the info file
         !first we write all the infos in order
  
         write(out_%unit_info,'(a,i10)')   'NAtoms     : ', target_%n_atoms
         write(out_%unit_info,'(a,i10)')   'NMolec     : ', target_%n_mol
         write(out_%unit_info,'(a,i10)')   'NVar       : ', target_%n_var
         if(field%static) then
            write(out_%unit_info,'(a)')   'Field      :     static'
         else if(field%dynamic) then
            write(out_%unit_info,'(a)')   'Field      :    dynamic'
         else
            write(out_%unit_info,'(a)')   'Field      :     absent'
         endif
         if(target_%kernel.eq.'coulomb') then
            write(out_%unit_info,'(a,i10)')   'Kernel     :    coulomb'
         else if(target_%kernel.eq.'ohno') then
            write(out_%unit_info,'(a,i10)')   'Kernel     :       ohno'
         else if(target_%kernel.eq.'gaussian') then
            write(out_%unit_info,'(a,i10)')   'Kernel     :   gaussian'
         else 
            write(out_%unit_info,'(a,i10)')   'Kernel     :       none'
         endif
         write(out_%unit_info,'(a,1x,a10)') 'Model      :' , trim(target_%name_)
         write(out_%unit_info,'(a,i10)')    'NumExFreq  : ', field%n_freq
         write(out_%unit_info,'(a,E13.6)')  'Input |E|  : ', field%e_0
         write(out_%unit_info,'(a,l10)')    'PrincAxes  : ', &
                                            control%principal_axis
         write(out_%unit_info,'(a,1x,a10)') 'ForceField :' , &
                                            trim(target_%forcefield)
         write(out_%unit_info,'(a)') 'Input Geometry (angstrom)'
      
         if(target_%name_.ne.'bem') then
            do i = 1, target_%n_atoms
               write(out_%unit_info,'(f4.1,2x,i8,3(1x,f25.16))')  &
                                target_%atomic_number(i), &
                                target_%i_mol(i),         &
                                target_%coord(1,i)/tobohr,       &
                                target_%coord(2,i)/tobohr,       &
                                target_%coord(3,i)/tobohr
            enddo 
      
            write(out_%unit_info,'(a,1x,a10)') "Parameters: ", &
                                               trim(target_%forcefield)
            write(out_%unit_info,'(a,i3)') "Num. AtomTypes : ", &
                                            target_%n_atomtypes
            if(target_%forcefield.eq.'fq'.and.target_%kernel.ne.'gaussian') then
               do i = 1, target_%n_atomtypes
                  write(out_%unit_info,'(a2,2(1x,f25.16))')     &
                                    target_%atom_type(i), &
                                    target_%chi(i),            &
                                    target_%eta(i)
               enddo 
            else if(target_%forcefield.eq.'fq'.and. &
                    target_%kernel.eq.'gaussian') then
               do i = 1, target_%n_atomtypes
                  write(out_%unit_info,'(a2,3(1x,f25.16))')     &
                                    target_%atom_type(i), &
                                    target_%chi(i),            &
                                    target_%eta(i),            &
                                    target_%r_q(i)
               enddo 
            else if(target_%forcefield.eq.'fqfmu') then
               do i = 1, target_%n_atomtypes
                  write(out_%unit_info,'(a2,5(1x,f25.16))')    &
                                    target_%atom_type(i),      &
                                    target_%chi(i),            &
                                    target_%eta(i),            &
                                    target_%alpha(i),          &
                                    target_%r_q(i),            &
                                    target_%r_mu(i)
               enddo 
            endif
  
            !NAtoms_per_IMol
            write(out_%unit_info,'(a)') 'Number of atoms per IMol'
            do i = 1, target_%n_mol
               write(out_%unit_info,'(i8,2x,i8)') i, target_%n_atoms_per_molecule(i)
            enddo
  
         else ! BEM 
            do i = 1, target_%n_var
               write(out_%unit_info,'(i8,3(1x,f25.16))')  &
                                i,                           &
                                target_%coord(1,i)/tobohr,   &
                                target_%coord(2,i)/tobohr,   &
                                target_%coord(3,i)/tobohr
            enddo 
         endif
  
         if(target_%name_.eq.'wfq'.or. &
            target_%name_.eq.'wfqfmu') then
  
            write(out_%unit_info,'(a)') "Dynamic parameters: 6"
            do i = 1, target_%n_atomtypes
               write(out_%unit_info,'(i1,3x,f30.16)') i, parameters%tau(i)
               write(out_%unit_info,'(i1,3x,f30.16)') i, parameters%sigma_0(i)      
               write(out_%unit_info,'(i1,3x,f30.16)') i, parameters%A_ij(i)         
               write(out_%unit_info,'(i1,3x,f30.16)') i, parameters%fermi_d(i,i)      
               write(out_%unit_info,'(i1,3x,f30.16)') i, parameters%fermi_s(i,i)      
               write(out_%unit_info,'(i1,3x,f30.16)') i, parameters%fermi_energy(i)
            enddo
  
         endif
                                      
         !Frequenze
         if(field%n_freq.gt.0) then
            write(out_%unit_info,'(a)') 'Frequencies (eV)'
            do i = 1, field%n_freq
               call FreqautoeV(field%freq(i),FreqeV)
               write(out_%unit_info,'(f25.16)') FreqeV
            enddo
         endif
  
      Close(out_%unit_info)
  
   end subroutine initialize_save_info

end module input_module
