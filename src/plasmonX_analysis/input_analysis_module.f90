!> Input analysis module   
!!
!! This module contains the subroutines for reading the input of the analysis
!!
!! Date         : 2025
!!
module input_analysis_module

   use parameters_module
   use input_module
   use output_analysis_module
   use control_analysis_module
   use grid_module

   implicit none 

   !public variables
   public inp_analysis 

   type, extends(inp_type) :: inp_analysis_type
     character(len=977) :: tar_file

     contains

     procedure :: check_file
     procedure :: read_file
     procedure :: parse_line
     procedure :: get_output
     procedure :: extract_tar_gz
     procedure :: print_input_info => print_input_info_inp_analysis
   end type inp_analysis_type
    
   type (inp_analysis_type) :: inp_analysis

contains


   !> Subroutine for check if the file exists and has the right extension
   !!    In/Out : input_analysis -- type
   subroutine check_file(inp_analysis)
     
      implicit none
  
      !input/output variables
      class(inp_analysis_type), intent(inout) :: inp_analysis

      !internal variables
      logical          :: exists
  
      inp_analysis%filename = "parameters.txt"
      
      inquire(file=inp_analysis%filename,exist=exists)
      
      if(.not.exists) call out_analysis%error("File "// &
                                              trim(inp_analysis%filename)// &
                                              " does not exist. Check the name")
      
   end subroutine check_file


   !> Subroutine for getting the name of the output file
   !!    In/Out : input_analysis -- type
   subroutine get_output(inp_analysis)
  
      implicit none 
  
      !input/output variables
      class(inp_analysis_type), intent(inout) :: inp_analysis
  
      !internal variables
      integer :: iost
      character(len=1000) :: line
  
      open(unit=inp_analysis%iunit,file=trim(inp_analysis%filename), &
           status="old",iostat=iost) 
  
         if (iost.ne.0) call out_analysis%error('Error in opening parameters.txt file')
         rewind(inp_analysis%iunit) 
         read(inp_analysis%iunit,'(A)', iostat=iost) line
         call inp_analysis%parse_line(line)
      close(inp_analysis%iunit) 
  
      !extract file name
      call out_analysis%out_file_fill(inp_analysis%tar_file)
  
   end subroutine get_output


   !> Subroutine for reading the parameters file
   !!    In/Out : input_analysis -- type
   subroutine read_file(inp_analysis)
  
      implicit none 
  
      !input/output variables
      class(inp_analysis_type), intent(inout)  :: inp_analysis
  
      !internal variables
      integer :: i
      integer :: iost
      integer :: nlines
      character(len=1000) :: line
  
      open(unit=inp_analysis%iunit,file=trim(inp_analysis%filename), &
           status="old",iostat=iost) 
         if (iost.ne.0) call out_analysis%error('Error in opening &
                                                &parameters.txt file')

         rewind(inp_analysis%iunit) 
          
         call get_number_lines(inp_analysis%iunit,nlines)
  
         rewind(inp_analysis%iunit) 
         do i = 1, nlines
            read(inp_analysis%iunit,'(A)', iostat=iost) line
            if (iost.ne.0) exit
            call inp_analysis%parse_line(line)
         end do
  
      close(inp_analysis%iunit) 
  
   end subroutine read_file


   !> Subroutine for parsing a line of the parameters file
   !!    In/Out : input_analysis -- type
   subroutine parse_line(inp_analysis, line)
     
      implicit none
  
      !input/output variables
      class(inp_analysis_type), intent(inout) :: inp_analysis
       
      !internal variables
      integer :: pos
      integer :: ierr
      character(len=1000), intent(in) :: line
      character(len=20)  :: variable       ! left = 
      character(len=977) :: value_variable !  = right
      character(len=6)  :: extention_input
  
      pos = index(line, '=')
      if (pos .gt. 0) then
          variable       = adjustl(line(1:pos-1))
          value_variable = adjustl(line(pos+1:))
      else
          return
      end if
  
      select case (trim(adjustl(variable)))
      
      case ("input_file")
         inp_analysis%tar_file = trim(adjustl(value_variable))
         extention_input = &
            inp_analysis%tar_file(len_trim(inp_analysis%tar_file)-5:)
         if(extention_input.ne.'tar.gz') &
            call out_analysis%error('Input file: '// &
                                    trim(inp_analysis%tar_file)// &
                                    ' is not a .tar.gz')
             
      case ("what")
         control_analysis%what = trim(adjustl(value_variable))
  
      case ("num_ex_freq")
         read(value_variable, *) field%n_freq
           
      case ("frequencies")
         if(field%n_freq.gt.0) then
            allocate(field%freq(field%n_freq),stat=ierr)
            if(ierr.gt.0) call out_%error('Insufficient space &
                                          &for Freq allocation')
            call array_clear(field%n_freq,field%freq)
            call parse_array(value_variable, field%freq)
         endif
           
      case ("scale_e0")
         read(value_variable, *) control_analysis%scale_e0
           
      case ("plane")
         control_analysis%plane = trim(adjustl(value_variable))
         if(control_analysis%plane.ne.'nu') &
            control_analysis%plane_requested = .true.
          
      case ("n_plane")
         read(value_variable, *) control_analysis%n_plane
           
      case ("step_plane")
         read(value_variable, *) control_analysis%step_plane
           
      case ("start")
         read(value_variable, *) control_analysis%start_plane
           
      case ("field_dir")
         field%polarization = trim(adjustl(value_variable))
         call field%set_polarization_variables()
           
      case ("nx_points")
         read(value_variable, *) grid%nx
           
      case ("ny_points")
         read(value_variable, *) grid%ny
           
      case ("nz_points")
         read(value_variable, *) grid%nz
           
      case ("volume")
         read(value_variable, *) control_analysis%volume_requested
           
      case ("min_grid")
         if(trim(value_variable).ne.'default') then
            grid%user_defined = .true.
            call parse_array(value_variable,grid%min_coord)
         endif
           
      case ("max_grid")
         if(trim(value_variable).ne.'default') then
            grid%user_defined = .true.
            call parse_array(value_variable,grid%max_coord)
         endif
           
      case ("offset_grid")
         read(value_variable, *) grid%offset
           
      case ("format_grid")
         read(value_variable, *) grid%format_

      case ("separate_q_mu")
         read(value_variable, *) control_analysis%separate_q_mu
           
      case ("n_omp_threads")
         read(value_variable, *) control_analysis%n_threads_omp
           
      end select
  
   end subroutine parse_line


   !> Subroutine for extracting files from tar.gz 
   !!    In/Out : input_analysis -- type
   subroutine extract_tar_gz(inp_analysis)
     
      implicit none
  
      class(inp_analysis_type), intent(inout) :: inp_analysis
  
      call execute_command_line( &
          '[ ! -d "tmp_analysis" ] && mkdir -p "tmp_analysis"')
      call execute_command_line( &
          "cp "//trim(inp_analysis%tar_file)//" tmp_analysis")
      call execute_command_line("tar -xf ./tmp_analysis/"//&
                                 trim(inp_analysis%tar_file)//&
                                 " -C tmp_analysis")
      
   end subroutine extract_tar_gz


   !> Subroutine for checking that correct files in tmp_analysis folder
   subroutine check_tmp_files()
     
      implicit none
  
      !internal variables
      integer            :: i
      character(len=300) :: file_to_check
      character(len=12)  :: string_freq
      logical            :: exists
  
      !check file .info
      file_to_check = out_analysis%tmp_folder // &
                      trim(out_analysis%root_filename) // ".info"
      inquire(file=trim(file_to_check),exist=exists)
  
      if(.not.exists) call out_analysis%error("file "//trim(file_to_check)// &
                                              " is not in the tar.gz.")
  
      !check .freq file if n_freq > 0
      if (field%n_freq.gt.0) then
         do i = 1, field%n_freq 
            write(string_freq,'(f12.10)') field%freq(i)
            string_freq = adjustl(string_freq)
            file_to_check = out_analysis%tmp_folder // &
                            trim(out_analysis%root_filename) // &
                            "-"//trim(string_freq)// ".freq"
            inquire(file=trim(file_to_check),exist=exists)
            if(.not.exists) call out_analysis%error("file "//trim(file_to_check)// &
                                                   " is not in the tar.gz.")
         enddo
      endif
      
   end subroutine check_tmp_files


   !> Subroutine for parsing an array
   !!    In/Out : string_in  -- String input
   !!    Output : array      -- array values
   subroutine parse_array(string_in, array)
     
      implicit none
  
      !input/output variables
      character(len=977), intent(inout) :: string_in
      real(dp), dimension(:), intent(out) :: array
  
      !internal variables
      integer :: init, end_
      integer :: i 
      integer :: n_variables
      character(len=977) :: string
  
      init = move_to_character(string_in,'[')+1
      end_ = move_to_character(string_in,']')-1
  
      string = string_in(init:end_)
      string = sweep_blanks(string)
  
      n_variables = 1
      do i = 1, len_trim(string)
          if (string(i:i) == ',') n_variables = n_variables + 1
      end do
  
      if(n_variables.ne.size(array)) & 
         call out_analysis%error("The number of variables is not correct in &
                                 &parse_array")
  
      do i = 1, n_variables
         string = string(1:len(string))
         if(i.eq.n_variables) then
            end_ = len(string)
         else
            end_ = move_to_character(string,',')-1
         endif
         read(string(1:end_),'(f25.16)') array(i)
         string  =  string(end_+2:len(string))
      enddo
         
   end subroutine parse_array


   !> Subroutine for reading info file
   subroutine read_info_file()
     
      implicit none
  
      !internal variables
      integer            :: unit_info = 14
      integer            :: iost
      integer            :: ierr
      integer            :: i
      integer            :: i_dummy
      character(len=300) :: file_info
      character(len=100) :: field_type
      character(len=200) :: kernel_type
      character(len=200) :: target_name
      character(len=200) :: forcefield_type
  
      file_info = out_analysis%tmp_folder // trim(out_analysis%root_filename)//&
                  ".info"
       
      open(unit=unit_info,file=trim(file_info),status="old",iostat=iost)
         if (iost.ne.0) call out_analysis%error('Error in opening '// &
                                                trim(file_info))
         !load info from the info file
         !first we read the info on the target all the infos in order
         do i = 1, 5 
            read(unit_info,*)
         enddo
  
         read(unit_info,'(13x, a)') target_name
         target_name = sweep_blanks(target_name)
  
         call assign_target(trim(target_name))
  
         rewind(unit_info)
         if (target_%name_.ne.'bem') then
             read(unit_info,'(13x, i10)')  target_%n_atoms
             read(unit_info,'(13x, i10)')  target_%n_mol
         else
             read(unit_info,'(13x, i10)')  i_dummy
             read(unit_info,'(13x, i10)')  i_dummy
         endif
         read(unit_info,'(13x, i10)')  target_%n_var
         read(unit_info,'(13x, a)')    field_type
          
         if(trim(adjustl(field_type)) .eq. 'static') then 
            field%static = .true.
         else if(trim(adjustl(field_type)) .eq. 'dynamic') then 
            field%dynamic = .true.
         endif
          
         read(unit_info,'(13x, a)') kernel_type
         kernel_type = sweep_blanks(kernel_type)
         if(target_%name_.ne.'bem') target_%kernel = trim(kernel_type)
          
         read(unit_info,*)  !skip model
         read(unit_info,*)  !skip n_freq 
         read(unit_info,'(13x, E13.6)') field%e_0
         read(unit_info,'(13x, l10)')   control%principal_axis
         read(unit_info,'(13x, a)')     forcefield_type
          
         forcefield_type = trim(forcefield_type)
         forcefield_type = sweep_blanks(forcefield_type)
         
         if(target_%name_.ne.'bem') target_%forcefield = trim(forcefield_type)
         field%e_0 = field%e_0 * control_analysis%scale_e0
  
         !allocation of variables
         if(target_%name_.ne.'bem') then
            allocate(target_%coord(3,target_%n_atoms),stat=ierr)
            if(ierr.ne.0) call out_%error('Not Enough Space for Coord')
         else
            allocate(target_%coord(3,target_%n_var),stat=ierr)
            if(ierr.ne.0) call out_%error('Not Enough Space for Coord')
         endif

         read(unit_info,*) !'Input Geometry (angstrom)'

         ! read info file
         if(target_%name_.ne.'bem') then
            allocate(target_%atomic_number(target_%n_atoms),stat=ierr)
            if(ierr.gt.0) call out_%error('Not Enough Space for atomic numbers')
            allocate(target_%atom_name(target_%n_atoms),stat=ierr)
            if(ierr.gt.0) call out_%error('Not Enough Space for atom name')
            allocate(target_%map_atomtypes(target_%n_atoms),stat=ierr)
            if(ierr.ne.0) call out_%error('Not Enough Space for map_atomtypes')
            allocate(target_%i_mol(target_%n_atoms),stat=ierr)
            if(ierr.ne.0) call out_%error('Not Enough Space for I_Mol')

            do i = 1, target_%n_atoms
               read(unit_info,'(f4.1,2x,i8,3(1x,f25.16))')  &
                                target_%atomic_number(i), &
                                target_%i_mol(i),         &
                                target_%coord(1,i),       &
                                target_%coord(2,i),       &
                                target_%coord(3,i)
            enddo 
            read(unit_info,*) ! "Parameters: ", trim(target_%forcefield)
            read(unit_info,'(17x,i3)') target_%n_atomtypes ! Num. AtomTypes 
  
            !allocation variables for atom types
            allocate(target_%atom_type(target_%n_atomtypes),stat=ierr)
            if(ierr.ne.0) call out_%error('Not Enough Space for atomtype')
            allocate(target_%chi(target_%n_atomtypes),stat=ierr)
            if(ierr.ne.0) call out_%error('Not Enough Space for chi')
            allocate(target_%eta(target_%n_atomtypes),stat=ierr)
            if(ierr.ne.0) call out_%error('Not Enough Space for eta')
            allocate(target_%n_atoms_per_molecule(target_%n_mol),stat=ierr)
            if(ierr.ne.0) call out_%error('Not Enough Space for &
                                          &N_atoms_per_molecule')
            if(target_%forcefield.eq.'fq'.and.target_%kernel.eq.'gaussian') then
               allocate(target_%r_q(target_%n_atomtypes),Stat=IErr)
               if(ierr.ne.0) call out_%error('Not Enough Space for R_q (FQ)')
            else if(target_%forcefield.eq.'fqfmu') then 
               allocate(target_%alpha(target_%n_atomtypes),Stat=IErr)
               if(ierr.ne.0) call out_%error('Not Enough Space for alpha (FMu)')
               allocate(target_%r_q(target_%n_atomtypes),Stat=IErr)
               if(ierr.ne.0) call out_%error('Not Enough Space for R_q (FQFMu)')
               allocate(target_%r_mu(target_%n_atomtypes),Stat=IErr)
               if(ierr.ne.0) call out_%error('Not Enough Space for R_mu(FQFMu)')
            endif
  
            if(target_%forcefield.eq.'fq'.and.target_%kernel.ne.'gaussian') then
               do i = 1, target_%n_atomtypes
                  read(unit_info,'(a2,2(1x,f25.16))')     &
                                    target_%atom_type(i), &
                                    target_%chi(i),       &
                                    target_%eta(i)
               enddo 
            else if(target_%forcefield.eq.'fq'.and. &
                    target_%kernel.eq.'gaussian') then
               do i = 1, target_%n_atomtypes
                  read(unit_info,'(a2,3(1x,f25.16))')     &
                                    target_%atom_type(i), &
                                    target_%chi(i),       &
                                    target_%eta(i),       &
                                    target_%r_q(i)
               enddo 
            else if(target_%forcefield.eq.'fqfmu') then
               do i = 1, target_%n_atomtypes
                  read(unit_info,'(a2,5(1x,f25.16))')     &
                                    target_%atom_type(i), &
                                    target_%chi(i),       &
                                    target_%eta(i),       &
                                    target_%alpha(i),     &
                                    target_%r_q(i),       &
                                    target_%r_mu(i)
               enddo 
            endif
  
            !NAtoms_per_IMol
            read(unit_info,*) !'Number of atoms per IMol'
            do i = 1, target_%n_mol
               read(unit_info,'(10x,i8)') target_%n_atoms_per_molecule(i)
            enddo
            !Assign atom names 
            do i = 1, target_%n_atoms
               target_%atom_name(i) = assign_atom_name(target_%atomic_number(i))
            enddo 
            !Assign map atomtypes
            call target_%assign_map_atomtypes(target_%atom_name)
         else ! READ BEM
            do i = 1, target_%n_var
               read(unit_info,'(i8,3(1x,f25.16))')  &
                               i_dummy,             &
                               target_%coord(1,i),  &
                               target_%coord(2,i),  &
                               target_%coord(3,i)
            enddo 
            allocate(target_%center_of_mass(3), stat=ierr)
            if(ierr.gt.0) call out_%error('not enough space for center of mass')
            target_%center_of_mass = target_%calculate_center_of_mass()
         endif
  
         if(target_%name_.eq.'wfq'.or. &
            target_%name_.eq.'wfqfmu') then
            allocate(parameters%tau           (target_%n_atomtypes), stat=ierr)
            if(ierr.gt.0) call out_%error('not enough space for tau')
            allocate(parameters%sigma_0       (target_%n_atomtypes), stat=ierr)
            if(ierr.gt.0) call out_%error('not enough space for sigma_0')
            allocate(parameters%a_ij          (target_%n_atomtypes), stat=ierr)
            if(ierr.gt.0) call out_%error('not enough space for a_ij')
            allocate(parameters%fermi_d       (target_%n_atomtypes, &
                                               target_%n_atomtypes), stat=ierr)
            if(ierr.gt.0) call out_%error('not enough space for fermi_d')
            allocate(parameters%fermi_s       (target_%n_atomtypes, &
                                               target_%n_atomtypes), stat=ierr)
            if(ierr.gt.0) call out_%error('not enough space for fermi_s')
            allocate(parameters%fermi_energy(target_%n_atomtypes), stat=ierr)
            if(ierr.gt.0) call out_%error('not enough space for Fermi energy')
            read(unit_info,*) !"Dynamic parameters: 6"
            do i = 1, target_%n_atomtypes
               read(unit_info,'(4x,f30.16)') parameters%tau(i)
               read(unit_info,'(4x,f30.16)') parameters%sigma_0(i)      
               read(unit_info,'(4x,f30.16)') parameters%A_ij(i)         
               read(unit_info,'(4x,f30.16)') parameters%fermi_d(i,i)      
               read(unit_info,'(4x,f30.16)') parameters%fermi_s(i,i)      
               read(unit_info,'(4x,f30.16)') parameters%fermi_energy(i)
            enddo
  
         endif
  
      close(unit_info)
      
   end subroutine read_info_file


   !> Subroutine for dynamic variables for a specific frequency
   !!    Input  : inp_analysis   -- type
   !!    Input  : string_freq    -- frequency string
   subroutine read_variables_w_freq_file(string_freq, variables_w)
     
      implicit none
  
      !input/output variables
      character(len=12), intent(in) :: string_freq
      complex(dp), dimension(target_%n_var,3), intent(inout) :: variables_w
  
      !internal variables
      integer            :: unit_freq = 14
      integer            :: iost
      integer            :: j, k
      character(len=300) :: freq_file
      character(len=100) :: line_freq
       
      freq_file = out_analysis%tmp_folder//trim(out_analysis%root_filename)// &
                  '-'//trim(string_freq)//'.freq'
  
      open(unit=unit_freq,file=trim(freq_file),status="old",iostat=iost)
         rewind(unit_freq)
         do j = 1, 3 ! x, y, z
            read(unit_freq,'(a)') line_freq ! read indentation line for x, y, z 
            do k = 1, target_%n_var
               read(unit_freq,'(13x,2(f25.16,3x))') variables_w(k,j)
            enddo
         enddo
      close(unit_freq)
  
      !scale for scale_e0
      variables_w = variables_w * control_analysis%scale_e0
  
   end subroutine read_variables_w_freq_file


   !> Subroutine for input info on output
   !!    Input  : inp_           -- type
   subroutine print_input_info_inp_analysis(inp_)
 
      implicit none
 
      !input/output variables
      class(inp_analysis_type), intent(in) :: inp_
 
      write(out_%iunit,'(1x,a)') "Input File         : "// &
         trim(inp_%tar_file)
      write(out_%iunit,'(1x,a,i3)') "OMP Threads        : ", &
                                    control_analysis%n_threads_omp
      write(out_%iunit,out_%sticks) 
      flush(out_%iunit)
 
      !field informations
      call control_analysis%print_info()
  
      !field informations
      call field%print_info(.true.)
  
      !Target informations
      call target_%print_info()
       
      if(bem%exists) call bem%print_parameters()
      if(bem%exists) call bem%print_coord_bem("Input")
      flush(out_analysis%iunit)
   
   end subroutine print_input_info_inp_analysis

end module input_analysis_module
