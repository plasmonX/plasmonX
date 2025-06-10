!> BEM module
!!
!! This module contains the subroutines for the BEM type
!!
!! Date         : 2025
!!
module bem_module
      
   use output_module
   use parameters_module
   use string_manipulation_module
   use array_manipulation_module
   use field_module
   use target_module
   use rhs_module
   use matrix_module

   implicit none

   !public variables
   public bem

   type, extends(target_type) :: bem_type
      ! number of vertices
      integer                                :: n_verts  
      ! parameters for legendre gauss lobatto integration
      integer                                :: npol1 = 8 
      integer                                :: npol2 = 6 
      ! vertices indeces [nodes]
      integer, dimension(:), allocatable     :: verts_index 
      ! faces [elements] 
      integer, dimension(:,:), allocatable   :: faces       
      ! Sphere radius for approx. Green function
      real(dp)                               :: sphere_r    
      ! factor for the normal vector [ = 1.0d0 or -1.0d0 ]
      ! MNPBEM        : -1.0d0
      ! GMSH(default) : 1.0d0
      real(dp)                               :: normal_factor = 1.0d0  
      ! tesserae area
      real(dp), dimension(:), allocatable    :: area        
      ! vertices coordinates (nodes)
      real(dp), dimension(:,:), allocatable  :: verts_coord 
      ! tangent vector _|_ normal 1direction
      real(dp), dimension(:,:), allocatable  :: tangent_1   
      ! tangent vector _|_ normal 2direction
      real(dp), dimension(:,:), allocatable  :: tangent_2   
      ! normal vector
      real(dp), dimension(:,:), allocatable  :: normal      
      !variable for exact green function
      ! tabulated values x
      real(dp), dimension(:), allocatable    :: x_tab       
      ! tabulated values y
      real(dp), dimension(:), allocatable    :: y_tab       
      ! tabulated values z
      real(dp), dimension(:), allocatable    :: w_tab       
      ! integration weights
      real(dp), dimension(:), allocatable    :: integration_weights  
      ! distance between polar and bem%coord
      real(dp), dimension(:), allocatable    :: distance_polar_coord 
      ! polar coordinates
      real(dp), dimension(:,:), allocatable  :: polar_coord          
      ! integration coordinates
      real(dp), dimension(:,:), allocatable  :: integration_coord    
      ! selected elements to be used for exact BEM 
      real(dp), dimension(:,:), allocatable  :: selected_elements    
      ! exact green function diagonal
      real(dp), dimension(:,:), allocatable  :: green_diagonal_exact 
      ! IEF green derivate
      real(dp), dimension(:,:), allocatable  :: green_der_ief
      ! IEF SA^-1 -- complex for final multiplication
      complex(dp), dimension(:,:), allocatable  :: SA_inv
      !BEM parameters
      character(len=200)                     :: gmsh_file = ""
      character(len=200)                     :: solvent = ""
      ! epsilon solvent -- default vacuum = 1.0d0
      real(dp)                               :: epsilon_solvent = one 
      character(len=200)                     :: green_function = ""
      ! file for the permittivity
      character(len=200)                     :: permittivity_file = ""
      ! permittivity type
      character(len=200)                     :: permittivity_type = ""
      ! variant
      character(len=200)                     :: variant = ""
      logical                                :: exists = .false.
      ! rescale the charges for constraints
      logical                                :: charge_constraint = .false. 
      complex(dp), dimension(:), allocatable :: permittivity

      contains

      !subroutines for defining the BEM system
      procedure :: read_permittivity_bem
      procedure :: set_fitted_permittivity
      procedure :: set_permittivity_from_file
      procedure :: read_gmsh_file
      procedure :: read_nodes_index_coord
      procedure :: read_faces
      procedure :: create_tesserae
      procedure :: print_coord_bem
      procedure :: print_parameters
      procedure :: dealloc => deallocate_bem
      !allocations
      procedure :: allocate_dynamic_field_memory   => &
                   allocate_dynamic_field_memory_bem
      procedure :: deallocate_dynamic_field_memory => &
                   deallocate_dynamic_field_memory_bem
      !construct matrices/vectors
      procedure :: construct_constant_matrix       => &
                   construct_constant_matrix_bem
      procedure :: construct_dynamic_field_rhs     => &
                   construct_dynamic_field_rhs_bem
      procedure :: construct_dynamic_matrix        => &
                   construct_dynamic_matrix_bem
      procedure :: gmres_diagonal_shift            => &
                   gmres_diagonal_shift_bem
      procedure :: new_gmres_diagonal_shift        => &
                   new_gmres_diagonal_shift_bem
      !printing & saving files
      procedure :: print_dynamic_field_variables   => &
                   print_dynamic_field_variables_bem
      !calculate properties
      procedure :: calculate_dynamic_polar         => &
                   calculate_dynamic_polar_bem
      procedure :: calculate_induced_field_at_point => &
                   calculate_induced_field_at_point_bem
      ! exact green function subroutines, defined in the interface
      procedure :: construct_exact_green_function
      procedure :: flat_integration_diagonal
      procedure :: quadrature_rule
      procedure :: tesserae_integration_polar_coordinates
      procedure :: flat_integration_off_diagonal
      procedure :: selection_of_refined_elements
      procedure :: refinement_diagonal_elements
      procedure :: refinement_off_diagonal_elements
      procedure :: final_refinement_green_function
   end type bem_type
    
   type (bem_type), target, save :: bem
 
   interface
 
      include "exact_green_function_interface.f90"
 
   end interface

contains

!-------------------------------------------------------------------------------
!  SUBROUTINES FOR DEFINING THE BEM SYSTEM
!-------------------------------------------------------------------------------

   !> Subroutine for reading the permittivity from file or tabulated data
   !!    In/Out : bem      -- bem_type
   subroutine read_permittivity_bem(bem)
  
      use field_module
       
      implicit none
       
      !input/output variables
      class(bem_type), intent(inout) :: bem
  
      call mem_man%alloc(bem%permittivity, field%n_freq, "bem%permittivity")

      !check if it to be read from tabulated data or file
      if(trim(bem%permittivity_type).ne.'none') then
         call bem%set_fitted_permittivity()
      else if(trim(bem%permittivity_file).ne.'none') then
         call bem%set_permittivity_from_file()
      endif
  
   end subroutine read_permittivity_bem


   !> Subroutine for setting the BEM epsilon from tabulated fitted data
   !!    In/Out : bem      -- bem_type
   subroutine set_fitted_permittivity(bem)
       
      implicit none
       
      !input/output variables
      class(bem_type), intent(inout) :: bem
  
      if(bem%permittivity_type.eq.'silver etchegoin') then
         call silver_etchegoin(out_%iunit,      &
                               'epsilon',       &
                               field%n_freq,   &
                               field%freq,     &
                               bem%permittivity)
      else if(bem%permittivity_type.eq.'gold etchegoin') then
         call gold_etchegoin(out_%iunit,      &
                             'epsilon',       &
                             field%n_freq,   &
                             field%freq,     &
                             bem%permittivity)
      else 
         call out_%error("Fitted permittivity for material: "//&
         trim(bem%permittivity_type)//" not yet implemented")
      endif

   end subroutine set_fitted_permittivity


   !> Subroutine for setting the BEM epsilon from file in input
   !!    In/Out : bem      -- bem_type
   subroutine set_permittivity_from_file(bem)
       
      implicit none
       
      !input/output variables
      class(bem_type), intent(inout) :: bem
  
      call read_permittivity_info_from_file(out_%iunit,              &
                                            bem%permittivity_file,   &
                                            field%n_freq,           &
                                            field%freq,             &
                                            bem%permittivity)
    
   end subroutine set_permittivity_from_file


   !> Subroutine for reading the GMSH file
   !!    In/Out : bem      -- bem_type
   subroutine read_gmsh_file(bem)
  
      implicit none
       
      !input/output variables
      class(bem_type), intent(inout) :: bem
  
      !internal variables
      integer :: i
      integer :: n_lines
      integer :: num_string_initial
      integer :: num_string_end
      integer :: unit_gmsh = 11
      integer :: io
      integer :: n_elements_tot
      logical :: check_file_gmsh
      logical :: found_string
      character(len=200) :: line
  
      !check for file_gmsh in input. 
      inquire(file=bem%gmsh_file, exist=check_file_gmsh)
      If(.not.check_file_gmsh) call out_%error("File "//trim(bem%gmsh_file)// &
                                               " does not exist")
      !read the GMSH file
      open(unit=unit_gmsh,file=bem%gmsh_file,status="OLD",IOSTAT=io,ERR=01)
         call get_number_lines(unit_gmsh,n_lines)
         rewind(unit_gmsh)
         !GMSH nodes
         call go_to_string(unit_gmsh,"$nodes",found_string,n_lines, &
                           num_string_initial)
         if(.not.found_string) call out_%error("Unrecognised .msh file - $Nodes&
                                        & not found - ASCII2 format is required")
         rewind(unit_gmsh)
         call go_to_string(unit_gmsh,"$endnodes",found_string,n_lines, &
                           num_string_end)
         if(.not.found_string) call out_%error("Unrecognised .msh file - &
                               &$EndNodes not found - ASCII2 format is required")
  
         !define where it is the beginning
         !end_string - initial_string -1 -1 (two -1 to skip first line)
         bem%n_verts = num_string_end - num_string_initial - 1 - 1 
  
         !allocate the array containing the nodes and its label number
         call mem_man%alloc(bem%verts_index, bem%n_verts, "bem%verts_index")
         call mem_man%alloc(bem%verts_coord, 3, bem%n_verts, "bem%verts_coord")

         rewind(unit_gmsh)
         call go_to_string(unit_gmsh,"$nodes",found_string,n_lines, &
                           num_string_initial)
  
         !read nodes indeces and xyz coordinates
         call bem%read_nodes_index_coord(unit_gmsh)
  
         !read $elements in .msh
         rewind(unit_gmsh)
         call go_to_string(unit_gmsh,"$elements",found_string,n_lines, &
                           num_string_initial)
         read(unit_gmsh,'(i8)') n_elements_tot
         !assign the number of tesserae
         bem%n_var = 0 
         do i = 1, n_elements_tot
            read(unit_gmsh,'(a)') line
            if(get_n_in_line(line).ge.8) bem%n_var = bem%n_var + 1
         enddo
         !allocate the array containing the elements information
         call mem_man%alloc(bem%faces, 3,bem%n_var, "bem%faces")
         rewind(unit_gmsh)
         call go_to_string(unit_gmsh,"$elements",found_string,n_lines, &
                           num_string_initial)
         !read the faces
         call bem%read_faces(unit_gmsh,n_elements_tot)
      01 continue 
      close(unit_gmsh)
         
   end subroutine read_gmsh_file


   !> Subroutine for reading the nodes indeces and coordinates
   !!    In/Out : bem       -- bem_type
   !!    Input  : unit_gmsh -- Unit of the file
   subroutine read_nodes_index_coord(bem,unit_gmsh)
  
      implicit none
       
      !input/output variables
      class(bem_type), intent(inout) :: bem
      integer, intent(in)            :: unit_gmsh
  
      !internal variables
      integer :: i
      integer :: ierr
      integer :: iend
      character(len=200)  :: line
      character(len=200)  :: line_x,line_y,line_z
      character(len=200)  :: line_index
  
      !Read the first line in dummy
      read(unit_gmsh,*) 
      do i = 1, bem%n_verts
        read(unit_gmsh,'(a)') line
        !eliminate any tab from string
        line = substitute_string(line,char(9)," ")
        !read index
        iend       = move_to_character(line,' ')
        line_index = line(1:iend)
        line_index = sweep_blanks(line_index)
        read(line_index, '(i8)',iostat=ierr) bem%verts_index(i)
        if(ierr.ne.0) call out_%error("node index "//trim(line_index)//&
                                      " is not an integer")
        !read X coord
        line      = line(iend+1:200)
        iend      = move_to_character(line,' ')
        line_X    = line(1:iend)
        line_X    = sweep_blanks(line_X)
        read(line_X, *) bem%verts_coord(1,i)
        !read Y coord
        line      = line(iend+1:200)
        iend      = move_to_character(line,' ')
        line_Y    = line(1:iend)
        line_Y    = sweep_blanks(line_Y)
        read(line_Y, *) bem%verts_coord(2,i)
        !read Z coord
        line      = line(iend+1:200)
        line_Z    = adjustl(line)
        read(line_Z, *) bem%verts_coord(3,i) 
      enddo

      !conversion of bem%verts to_bohr
      call array_scale(ToBohr,3*bem%n_verts,bem%verts_coord)
  
   end subroutine read_nodes_index_coord


   !> Subroutine for reading the faces (elements triangles)
   !!    In/Out : bem       -- bem_type
   !!    Input  : unit_gmsh -- Unit of the file
   !!    Input  : n_lines   -- number of lines to be read
   subroutine read_faces(bem,unit_gmsh,n_lines)
  
      implicit none
       
      !input/output variables
      class(bem_type), intent(inout) :: bem
      integer, intent(in)            :: unit_gmsh
      integer, intent(in)            :: n_lines
  
      !internal variables
      integer :: i,j,k
      integer :: ierr
      integer :: iend
      character(len=200)  :: line
      character(len=200)  :: line_1,line_2,line_3
  
      read(unit_gmsh,*) 
      k = 0 
      do i = 1, n_lines
        read(unit_gmsh,'(a)') line
        if(get_n_in_line(line).ge.8) then
           k = k + 1
           !eliminate any tab from string
           line = substitute_string(line,char(9)," ") 
           !move to start reading the last three numbers
           do j = 1,5
              iend   = move_to_character(line,' ')
              line   = line (iend+1:)
           enddo
           !read first vertex
           iend   = move_to_character(line,' ')
           line_1 = line(1:iend)
           line_1 = sweep_blanks(line_1)
           read(line_1, '(i8)',iostat=ierr) bem%faces(1,k)
           if(ierr.ne.0) call out_%error("element 1 "//trim(line_1)//&
                                         " is not an integer")
           !read second vertex
           line   = line(iend+1:)
           iend   = move_to_character(line,' ')
           line_2 = line(1:iend)
           line_2 = sweep_blanks(line_2)
           read(line_2, '(i8)',iostat=ierr) bem%faces(2,k)
           if(ierr.ne.0) call out_%error("element 2 "//trim(line_2)//&
                                         " is not an integer")
           !read third vertex
           line   = line(iend+1:)
           iend   = move_to_character(line,' ')
           line_3 = line(1:iend)
           line_3 = sweep_blanks(line_3)
           read(line_3, '(i8)',iostat=ierr) bem%faces(3,k)
           if(ierr.ne.0) call out_%error("element 3 "//trim(line_3)//&
                                         " is not an integer")
        endif ! end if 8 elements
      enddo
      
   end subroutine read_faces


   !> Subroutine for printing the BEM parameters
   !!    In/Out : bem      -- bem_type
   subroutine create_tesserae(bem)
  
      implicit none
       
      !input/output variables
      class(bem_type) :: bem
  
      !internal variables
      integer :: i,k
      real(dp), dimension(3)   :: a, b
      real(dp), dimension(3)   :: cross
      real(dp), dimension(3,3) :: tmp_coord
  
      !allocations        
      call mem_man%alloc(bem%coord,3,bem%n_var,"bem%coord")
      call mem_man%alloc(bem%tangent_1,3,bem%n_var,"bem%tangent_1")
      call mem_man%alloc(bem%tangent_2,3,bem%n_var,"bem%tangent_2")
      call mem_man%alloc(bem%normal,3,bem%n_var,"bem%normal")
      call mem_man%alloc(bem%area,bem%n_var,"bem%area")
    
      !fill centroid coordinates array in atomic units
      tmp_coord = zero
      do i = 1, bem%n_var
         do k = 1, bem%n_verts
            if (bem%verts_index(k) .eq. bem%faces(1,i)) then
               tmp_coord(1,1) = bem%verts_coord(1,k)
               tmp_coord(1,2) = bem%verts_coord(2,k)
               tmp_coord(1,3) = bem%verts_coord(3,k)
            else if (bem%verts_index(k) .eq. bem%faces(2,i)) then
               tmp_coord(2,1) = bem%verts_coord(1,k)
               tmp_coord(2,2) = bem%verts_coord(2,k)
               tmp_coord(2,3) = bem%verts_coord(3,k)
            else if (bem%verts_index(k) .eq. bem%faces(3,i)) then
               tmp_coord(3,1) = bem%verts_coord(1,k)
               tmp_coord(3,2) = bem%verts_coord(2,k)
               tmp_coord(3,3) = bem%verts_coord(3,k)
            endif
         enddo
         !bem coordinates 
         bem%coord(1,i) = (tmp_coord(1,1)+tmp_coord(2,1)+tmp_coord(3,1))/three
         bem%coord(2,i) = (tmp_coord(1,2)+tmp_coord(2,2)+tmp_coord(3,2))/three
         bem%coord(3,i) = (tmp_coord(1,3)+tmp_coord(2,3)+tmp_coord(3,3))/three
         !tangential vector 1
         a(1) = tmp_coord(2,1) - tmp_coord(1,1)
         a(2) = tmp_coord(2,2) - tmp_coord(1,2)
         a(3) = tmp_coord(2,3) - tmp_coord(1,3)
         bem%tangent_1(:,i) =  - a(:)/norm2(a)
         !cross vector
         b(1) = tmp_coord(3,1) - tmp_coord(1,1)
         b(2) = tmp_coord(3,2) - tmp_coord(1,2)
         b(3) = tmp_coord(3,3) - tmp_coord(1,3)
         cross = cross_product(a,b)
         !Area of the tesserae
         bem%area(i) = norm2(cross)/two
         !Normal vector
         bem%normal(:,i) = bem%normal_factor * cross/norm2(cross) 
         !tangential vector 2
         bem%tangent_2(:,i) = cross_product(bem%normal(:,i),bem%tangent_1(:,i))
      enddo
      !printing
      if(out_%ivrb.ge.4) then
         call out_%print_matrix('BEM Tangent Vector 1',bem%tangent_1,3,&
                                bem%n_var)
         call out_%print_matrix('BEM Tangent Vector 2',bem%tangent_2,3,&
                                 bem%n_var)
         call out_%print_matrix('BEM Normal Vector',bem%normal,3,bem%n_var)
         call array_scale(1/ToBohr**2,bem%n_var,bem%area)
         call out_%print_matrix('BEM Area (Ang^2)', bem%area,1,bem%n_var)
         call array_scale(ToBohr**2,bem%n_var,bem%area)
      endif
        
   end subroutine create_tesserae


   !> Subroutine for printing the BEM coordinates of the tesserae centroids
   !!    Input  : bem      -- bem_type
   !!    Input  : what     -- what to calculate
   subroutine print_coord_bem(bem,what)
  
      implicit none
       
      !input/output variables
      class(bem_type), intent(in)  :: bem
      character(len=*), intent(in) :: what
  
      !internal variables
      integer :: i 
      character(len=76) :: format_1 = "(7x,'Tess.',21x,'X',19x,'Y',19x,'Z')"
      character(len=50) :: format_2 = "(3x,i9,5x,3(f20.6))"
  
      write(out_%iunit,'(/,27x,a,/)') trim(what)//" Tesserae Centroids (Ang)"
      write(out_%iunit,out_%sticks)
      write(out_%iunit,format_1)
      write(out_%iunit,out_%sticks)
      if(bem%n_var.le.20000.or.out_%ivrb.eq.3) then
         do i = 1, bem%n_var
            write(out_%iunit,format_2) i,                    &
                                       bem%coord(1,i)*ToAng, &
                                       bem%coord(2,i)*ToAng, &
                                       bem%coord(3,i)*ToAng
         enddo
      else
         do i = 1, 2000
            write(out_%iunit,format_2) i,                    &
                                       bem%coord(1,i)*ToAng, &
                                       bem%coord(2,i)*ToAng, &
                                       bem%coord(3,i)*ToAng
         enddo
         write(out_%iunit,'(a)') "           ."
         write(out_%iunit,'(a)') "           ."
         write(out_%iunit,'(a)') "           ."
         do i = bem%n_var-2000, bem%n_var
            write(out_%iunit,format_2) i,                    &
                                       bem%coord(1,i)*ToAng, &
                                       bem%coord(2,i)*ToAng, &
                                       bem%coord(3,i)*ToAng
         enddo
      endif
      write(out_%iunit,out_%sticks)
       
   end subroutine print_coord_bem


   !> Subroutine for printing the BEM parameters
   !!    Input : bem      -- bem_type
   subroutine print_parameters(bem)
  
      implicit none
       
      class(bem_type), intent(in) :: bem
  
      write(out_%iunit,'(1x,a)') "GMSH File          : "//trim(bem%gmsh_file)
      write(out_%iunit,'(1x,a)') "BEM Variant        : "//&
                                   trim(bem%variant)
      write(out_%iunit,'(1x,a)') "Green Function     : "//&
                                   trim(bem%green_function)
      write(out_%iunit,'(1x,a)') "Permittivity       : "//&
                                  trim(bem%permittivity_type)
      write(out_%iunit,'(1x,a)') "Permittivity File  : "//&
                                  trim(bem%permittivity_file)
      write(out_%iunit,'(1x,a)') "Solvent            : "//trim(bem%solvent)
      write(out_%iunit,'(1x,a,f8.3)') "Epsilon Solvent    : ", &
                                        bem%epsilon_solvent
      write(out_%iunit,'(1x,a,f8.3,a)') "Sphere Radius      : ", &
                                         bem%sphere_r, " a.u."
      write(out_%iunit,'(1x,a,f8.3)')   "Normal Factor      : ", &
                                         bem%normal_factor
      write(out_%iunit,'(1x,a,i6)') "Number of Vertices :", bem%n_verts
      write(out_%iunit,'(1x,a,i6)') "Number of Tesserae :", bem%n_var
      write(out_%iunit,'(1x,a,l1)') "Charge Constraint  : ", &
                                     bem%charge_constraint 
      write(out_%iunit,out_%sticks)
       
   end subroutine print_parameters

!-------------------------------------------------------------------------------
!  SUBROUTINES FOR ALLOCATIONS
!-------------------------------------------------------------------------------

   !> Subroutine for deallocating BEM variables
   !!    In/Out : bem      -- bem_type
   subroutine deallocate_bem(target_)
  
      implicit none
       
      !input/output variables
      class(bem_type), intent(inout) :: target_
  
      if(target_%name_.ne.'bem') &
         call out_%error("The target name is not BEM in deallocated BEM")
      call mem_man%dealloc(bem%verts_coord, "bem%verts_coord")
      call mem_man%dealloc(bem%verts_index, "bem%verts_index")
      call mem_man%dealloc(bem%faces, "bem%faces")
      call mem_man%dealloc(bem%area, "bem%area")
      call mem_man%dealloc(bem%coord, "bem%coord")
      call mem_man%dealloc(bem%tangent_1, "bem%tangent_1")
      call mem_man%dealloc(bem%tangent_2, "bem%tangent_2")
      call mem_man%dealloc(bem%normal, "bem%normal")
      call mem_man%dealloc(bem%permittivity, "bem%permittivity")
      if(bem%green_function.eq.'exact') then
         call mem_man%dealloc(bem%x_tab, "bem%x_tab")
         call mem_man%dealloc(bem%y_tab, "bem%y_tab")
         call mem_man%dealloc(bem%w_tab, "bem%w_tab")
         call mem_man%dealloc(bem%integration_weights, &
                              "bem%integration_weights")
         call mem_man%dealloc(bem%polar_coord, "bem%polar_coord")
         call mem_man%dealloc(bem%distance_polar_coord, &
                              "bem%distance_polar_coord")
         call mem_man%dealloc(bem%integration_coord, &
                              "bem%integration_coord")
         call mem_man%dealloc(bem%selected_elements, &
                              "bem%selected_elements")
         call mem_man%dealloc(bem%green_diagonal_exact, &
                              "bem%green_diagonal_exact")
      endif
         
   end subroutine deallocate_bem


   !> Subroutine for allocating dynamic field variables for memory (iter/inv)
   !!    In/Out  : target_             -- bem
   !!    In/Out  : matrix_constant     -- part of matrix_w that is constant
   !!    In/Out  : matrix_w            -- dynamic Matrix
   subroutine allocate_dynamic_field_memory_bem(target_,matrix_constant, &
                                                matrix_w)
  
      implicit none
       
      !input/output variables
      class(bem_type)    :: target_
      real(dp), dimension(:,:), allocatable    :: matrix_constant
      complex(dp), dimension(:,:), allocatable :: matrix_w
  
      target_%n_rhs = 3
      call mem_man%alloc(matrix_w, target_%n_var, target_%n_var, "matrix_w") 
      call mem_man%alloc(matrix_constant, target_%n_var, target_%n_var, &
                        "matrix_constant") 
      if(trim(target_%variant).eq.'iefpcm') then
        call mem_man%alloc(bem%green_der_ief, target_%n_var, target_%n_var, &
                           "green_der_ief")
        call mem_man%alloc(bem%SA_inv, target_%n_var, target_%n_var, "sa_inv")
      endif
       
   end subroutine allocate_dynamic_field_memory_bem


   !> Subroutine for allocating dynamic field variables for memory (iter/inv)
   !!    In/Out  : target_             -- bem
   !!    In/Out  : matrix_constant     -- part of matrix_w that is constant
   !!    In/Out  : matrix_w            -- dynamic Matrix
   subroutine deallocate_dynamic_field_memory_bem(target_,matrix_w, &
                                                  matrix_constant)
  
      implicit none
       
      !input/output variables
      class(bem_type)    :: target_
      complex(dp), dimension(:,:), allocatable :: matrix_w
      real(dp), dimension(:,:), allocatable, optional :: matrix_constant
  
      call mem_man%dealloc(matrix_w, "matrix_w") 
      if(present(matrix_constant)) &
         call mem_man%dealloc(matrix_constant,"matrix_constant") 
      if(trim(target_%variant).eq.'iefpcm') then
        call mem_man%dealloc(bem%green_der_ief, "green_der_ief")
        call mem_man%dealloc(bem%SA_inv, "sa_inv")
      endif
       
   end subroutine deallocate_dynamic_field_memory_bem

!-------------------------------------------------------------------------------
!  SUBROUTINES FOR CONSTRUCTING MATRICES/VECTORS
!-------------------------------------------------------------------------------

   !> Subroutine for constructing the constant matrix
   !! Here, the Green function is constructed 
   !!    In/Out : target_         -- wfqfmu
   !!    In/Out : matrix_constant -- constant part of the dynamic matrix
   subroutine construct_constant_matrix_bem(target_,matrix_constant)
  
      implicit none
      !input/output variables
      class(bem_type), intent(inout) :: target_
      real(dp), dimension(target_%n_var,target_%n_var), intent(inout) :: &
                                                                 matrix_constant
  
      !internal variables
      integer                               :: i, j
      real(dp)                              :: x_ij, y_ij, z_ij
      real(dp)                              :: dist_ij, dist_ij_3
      real(dp), parameter                   :: diag_param = 1.0694d0
      real(dp), dimension(:), allocatable   :: tmp_a_inv
  
      if(bem%green_function.eq.'exact') & 
         call bem%construct_exact_green_function()
      
      !IEFPCM - We need some more matrices
      if(trim(bem%variant).eq."iefpcm") &
         call mem_man%alloc(tmp_a_inv, target_%n_var, "tmp_a_inv")

      do i = 1, target_%n_var
         do j = 1, target_%n_var
            x_ij = target_%coord(1,i) - target_%coord(1,j)   
            y_ij = target_%coord(2,i) - target_%coord(2,j)   
            z_ij = target_%coord(3,i) - target_%coord(3,j)   
            dist_ij   = dsqrt(x_ij**2 + y_ij**2 + z_ij**2)
            dist_ij_3 = dist_ij**3
            !Exact Green function
            if(target_%green_function.eq.'exact') then 
               if(bem%selected_elements(i,j).ge.zero) then  
                  matrix_constant(i,j) = target_%area(j)* &
                                         (x_ij*target_%normal(1,i) + &
                                          y_ij*target_%normal(2,i) + &
                                          z_ij*target_%normal(3,i) )/dist_ij_3
               else
                  !diagonal elements
                  if(i.eq.j) then
                     matrix_constant(i,j) = - (bem%green_diagonal_exact(1,i)* &
                                               bem%normal(1,i)              + & 
                                               bem%green_diagonal_exact(2,i)* &
                                               bem%normal(2,i)              + &
                                               bem%green_diagonal_exact(3,i)* &
                                               bem%normal(3,i) )
                  !off diagonal elements
                  else
                     matrix_constant(i,j) = &
                                    - bem%refinement_off_diagonal_elements(i,j)
                  endif !diagonal elements
               endif ! bem%selected_elements
            !green_function .eq. 'approx'
            else 
               if(i.ne.j) then !off diagonal elements
                  matrix_constant(i,j) = target_%area(j)* &
                                         (x_ij*target_%normal(1,i) + &
                                          y_ij*target_%normal(2,i) + &
                                          z_ij*target_%normal(3,i) )/dist_ij_3
               else 
                  !diag = 1.0694 * sqrt(4 * pi * area) / (two * radius_sphere)
                  matrix_constant(i,i) = diag_param * & 
                            sqrt(four*pi*target_%area(i))/(two*target_%sphere_r)
               endif !diagonal
            endif !green function 'approx' or 'exact'
            !iefpcm
            if (trim(bem%variant).eq.'iefpcm') then
               if(i.eq.j) then
                  tmp_a_inv(i)    = one/target_%area(i)
                  bem%SA_inv(i,j) = dcmplx(diag_param * &
                                    sqrt(four*pi*target_%area(i)),zero)  
                  !save -(2pi*I + D) 
                  bem%green_der_ief(i,j) = matrix_constant(i,j) - two*pi
               else
                  bem%SA_inv(i,j) = dcmplx(target_%area(j)/dist_ij, zero)
                  !save D off diagonal
                  bem%green_der_ief(i,j) = matrix_constant(i,j)
               endif
            endif
         enddo !n_var
      enddo !n_var

      if(trim(bem%variant).eq.'iefpcm') then
         do i = 1, target_%n_var
            do j = 1, target_%n_var
               bem%SA_inv(i,j) = bem%SA_inv(i,j) * tmp_a_inv(j)
            enddo
         enddo
         call mem_man%dealloc(tmp_a_inv, "tmp_a_inv")
         if(out_%ivrb.ge.4) call out_%print_matrix('SA^-1',dble(bem%sa_inv), &
                                                   target_%n_var, target_%n_var)
      endif

      !printing
      if(out_%ivrb.ge.4) & 
         call out_%print_matrix('BEM Green Der', -matrix_constant, &
                                 target_%n_var, target_%n_var)
      !final refinement for Exact Green Function
      if(bem%green_function.eq.'exact') then 
          call bem%final_refinement_green_function(matrix_constant)
          if(out_%ivrb.ge.4) &
             call out_%print_matrix('BEM Green Der Refined', -matrix_constant, &
                                    target_%n_var, target_%n_var)
      endif
  
   end subroutine construct_constant_matrix_bem


   !> Subroutine for constructing the wfqfmu RHS on memory
   !!    In/Out : target_         -- wfqfmu
   !!    In/Out : rhs_w           -- dynamic RHS
   subroutine construct_dynamic_field_rhs_bem(target_,rhs_w)
  
      use field_module
  
      implicit none
  
      !input/output variables
      class(bem_type), intent(inout) :: target_ 
      complex(dp), dimension(target_%n_var,3), intent(inout) :: rhs_w
  
      !internal variables
      integer           :: i
      real(dp), dimension(3,3) :: external_field
      real(dp), dimension(:,:), allocatable :: rhs_static
      real(dp), dimension(:,:), allocatable :: potential
  
      !Polarization Vector
      external_field(1,1) = field%e_0
      external_field(2,1) = Zero
      external_field(3,1) = Zero

      external_field(1,2) = Zero
      external_field(2,2) = field%e_0
      external_field(3,2) = Zero

      external_field(1,3) = Zero
      external_field(2,3) = Zero
      external_field(3,3) = field%e_0

      call mem_man%alloc(rhs_static,target_%n_var,3,"rhs_static")
      if(trim(bem%variant) .eq. 'iefpcm') then
         !-(2pi*I + D) * (-Field * coord_xyz)
         ! rhs_static = matmul(bem%green_der_ief,            &
         !              transpose(matmul(-field, bem%coord)))
         !potential = -external_field * coord
         call mem_man%alloc(potential, 3, target_%n_var,"potential") 
         call dgemm('N', 'N',       &
                     3,             &
                     target_%n_var, &
                     3,             &     
                     -one,          &
                     external_field,&
                     3,             &
                     bem%coord,     &
                     3,             &
                     zero,          &
                     potential,     &
                     3)
         call dgemm('N', 'T',          &
                    target_%n_var,     &
                    3,                 &
                    target_%n_var,     &
                    one,               &
                    bem%green_der_ief, &
                    target_%n_var,     &
                    potential,         &
                    3,                 &
                    zero,              &
                    rhs_static,        &
                    target_%n_var)                     
         call mem_man%dealloc(potential,"potential") 
      else
         !Projection over the normal direction \hat{n}\cdot\vec{E}
         call dgemm('T', 'N',        &
                    target_%n_var,   &
                    3,               &
                    3,               &
                    -one,            &
                    target_%normal,  &
                    3,               &
                    external_field,  &
                    3,               &
                    zero,            &
                    rhs_static,      &
                    target_%n_var)
      endif

      do i=1,target_%n_var
         rhs_w(i,:) = dcmplx(rhs_static(i,:),zero)
      enddo

      if (out_%ivrb.ge.3) call print_dynamic_field_rhs(target_,rhs_w)
  
      call mem_man%dealloc(rhs_static,"rhs_static")
  
   end subroutine construct_dynamic_field_rhs_bem


   !> Subroutine for constructing the BEM dynamic matrix
   !!    In/Out : target_         -- BEM
   !!    In/Out : i_freq          -- index of frequency
   !!    Input  : matrix_constant -- constant part of the dynamic matrix
   !!    In/Out : matrix_w        -- dynamic matrix
   subroutine construct_dynamic_matrix_bem(target_,i_freq,matrix_constant, &
                                           matrix_w)
  
      implicit none
  
      !input/output variables
      class(bem_type), intent(inout) :: target_
      integer, intent(in)            :: i_freq
      real(dp), dimension(target_%n_var,target_%n_var), intent(in) :: &
                                                                 matrix_constant
      complex(dp), dimension(target_%n_var,target_%n_var), intent(inout) :: &
                                                                        matrix_w
  
      !internal variables
      integer           :: i,j
      complex(dp), dimension(target_%n_var,target_%n_var) :: matrix_out
  
      do i = 1, target_%n_var
         do j = 1, target_%n_var
            matrix_w(i,j) = dcmplx(matrix_constant(i,j),zero)
         enddo
         !diagonal = 2 * pi (eps_solv + eps_w) / (eps_solv - eps_w)
         matrix_w(i,i) = matrix_w(i,i) + two*pi* &
                         (bem%epsilon_solvent + bem%permittivity(i_freq))/ &
                         (bem%epsilon_solvent - bem%permittivity(i_freq))
      enddo
      if(trim(bem%variant).eq.'iefpcm') then
         !call mem_man%alloc(matrix_out, target_%n_var, target_%n_var, &
         !                   "matrix_out")
         call zgemm('N', 'N',           &
                     target_%n_var,     &
                     target_%n_var,     &
                     target_%n_var,     &
                     dcmplx(one,zero),  &
                     matrix_w,          &
                     target_%n_var,     &
                     bem%SA_inv,        &
                     target_%n_var,     &
                     dcmplx(zero,zero), &
                     matrix_out,          &
                     target_%n_var)
         matrix_w = -matrix_out
         !call mem_man%dealloc(matrix_out, "matrix_out")
      endif
      !printing
      if(out_%ivrb.ge.3) then
         call out_%print_matrix("BEM Matrix: Real Part",dble(matrix_w), &
                                target_%n_var,target_%n_var)
         call out_%print_matrix("BEM Matrix: IMag Part",dimag(matrix_w), &
                                target_%n_var,target_%n_var)
      endif
  
   end subroutine construct_dynamic_matrix_bem


   !> Subroutine for applying the diagonal shift for GMRES algorithm
   !!    In/Out  : target_      -- model 
   !!    Input   : i_freq       -- index of the frequency
   !!    Input   : colx         -- index of the column x (GMRES)
   !!    Input   : colz         -- index of the column z (GMRES)
   !!    Input   : lwork        -- dimension of the work array
   !!    In/Out  : work         -- array work
   subroutine gmres_diagonal_shift_bem(target_,i_freq,colx,colz,lwork,work)
  
      implicit none
       
      !input/output variables
      class(bem_type)    :: target_
      integer  :: i_freq
      integer  :: colx
      integer  :: colz
      integer  :: lwork
      complex(dp), dimension(lwork), intent(inout) :: work
  
      !internal variables
      integer  :: i
      integer  :: start_x
      integer  :: start_z
      complex(dp) :: z_i
  
      start_z = colz - 1
      start_x = colx - 1
      do i = 1, target_%n_var
         z_i = two*pi*(bem%epsilon_solvent + bem%permittivity(i_freq))/ &
                      (bem%epsilon_solvent - bem%permittivity(i_freq))
         work(start_z+i) = work(start_z + i) + z_i * work(start_x + i)
      enddo
       
   end subroutine gmres_diagonal_shift_bem


   !> Subroutine for applying the diagonal shift for GMRES algorithm
   !!    In/Out  : target_      -- model 
   !!    Input   : i_freq       -- index of the frequency
   !!    Input   : r            -- rhs
   !!    In/Out  : w            -- solution
   subroutine new_gmres_diagonal_shift_bem(target_,i_freq,r,w)
  
      implicit none
       
      !input/output variables
      class(bem_type)    :: target_
      integer  :: i_freq
      complex(dp), dimension(target_%n_var), intent(in) :: r
      complex(dp), dimension(target_%n_var), intent(inout) :: w
  
      !internal variables
      integer  :: i
      complex(dp) :: z_i
  
      do i = 1, target_%n_var
         z_i = two*pi*(bem%epsilon_solvent + bem%permittivity(i_freq))/ &
                      (bem%epsilon_solvent - bem%permittivity(i_freq))
         w(i) = w(i) + z_i * r(i)
      enddo
       
   end subroutine new_gmres_diagonal_shift_bem

!-------------------------------------------------------------------------------
!  SUBROUTINES FOR PRINTING & SAVING FILES
!-------------------------------------------------------------------------------

   !> Subroutine for printing the variables when dynamic field is applied
   !!    Input   : target_      -- model 
   !!    Input   : variables_w  -- w-variables (N_var, 3 [x,y,z]) 
   subroutine print_dynamic_field_variables_bem(target_,variables_w)
  
      implicit none
       
      !input/output variables
      class(bem_type), intent(in) :: target_
      complex(dp), dimension(target_%n_var,3), intent(in) :: variables_w
  
      !internal variables
      integer :: i
      complex(dp) :: sum_X = dcmplx(zero,zero)
      complex(dp) :: sum_Y = dcmplx(zero,zero)
      complex(dp) :: sum_Z = dcmplx(zero,zero)
      character(len=21) :: format_1 = "(38x,'BEM Charges',/)"
      character(len=64) :: format_2 = "(13x,'X Component',14x,'Y Component',&
                                      &14x,'Z Component',/)"
      character(len=60) :: format_3 = "(1x,I4,2x,2(E10.3,2x,E10.3,1X,'|',1x),&
                                      &E10.3,2X,E10.3)"
      character(len=64) :: format_4 = "(/,'ChErr',2x,2(E10.3,2x,E10.3,1X,'|',&
                                      &1x),E10.3,2X,E10.3)"
  
      if(out_%ivrb.ge.2) then
         write(out_%iunit,format_1)
         write(out_%iunit,format_2)
         do i = 1, target_%n_var
            write(out_%iunit,format_3) i,variables_w(i,:)
            sum_x = sum_x + variables_w(i,1)
            sum_y = sum_y + variables_w(i,2)
            sum_z = sum_z + variables_w(i,3)
         enddo
         write(out_%iunit,format_4) sum_x, sum_y, sum_z
         write(out_%iunit,out_%sticks) 
      endif
       
   end subroutine print_dynamic_field_variables_bem

!-------------------------------------------------------------------------------
!  SUBROUTINES FOR CALCULATING PROPERTIES
!-------------------------------------------------------------------------------

   !> Subroutine for calculating the dynamic polar
   !!    In/Out  : target_      -- model 
   !!    Input   : i_freq       -- index of the frequency
   !!    Input   : variables_w  -- w-variables
   !!    Output  : polar_w      -- dynamic polar
   subroutine calculate_dynamic_polar_bem(target_,i_freq,variables_w,polar_w)
  
      implicit none
       
      !input/output variables
      class(bem_type), intent(inout) :: target_
      integer, intent(in) :: i_freq
      complex(dp), dimension(target_%n_var,3), intent(in) :: variables_w
      complex(dp), dimension(3,3), intent(inout)          :: polar_w
  
      !internal variables
      integer :: i
      complex(dp), dimension(3,3) :: internal_polar
  
      internal_polar = dcmplx(zero,zero)
     if(trim(bem%variant).eq.'iefpcm') then
        do i = 1, target_%n_var
          internal_polar(1,1) = internal_polar(1,1) + &
                                variables_w(i,1)*target_%coord(1,i)/field%e_0
          internal_polar(1,2) = internal_polar(1,2) + &
                                variables_w(i,2)*target_%coord(1,i)/field%e_0
          internal_polar(1,3) = internal_polar(1,3) + &
                                variables_w(i,3)*target_%coord(1,i)/field%e_0
          internal_polar(2,1) = internal_polar(2,1) + &
                                variables_w(i,1)*target_%coord(2,i)/field%e_0
          internal_polar(2,2) = internal_polar(2,2) + &
                                variables_w(i,2)*target_%coord(2,i)/field%e_0
          internal_polar(2,3) = internal_polar(2,3) + &
                                variables_w(i,3)*target_%coord(2,i)/field%e_0
          internal_polar(3,1) = internal_polar(3,1) + &
                                variables_w(i,1)*target_%coord(3,i)/field%e_0
          internal_polar(3,2) = internal_polar(3,2) + &
                                variables_w(i,2)*target_%coord(3,i)/field%e_0
          internal_polar(3,3) = internal_polar(3,3) + &
                                variables_w(i,3)*target_%coord(3,i)/field%e_0
        enddo
      else
         do i = 1, target_%n_var
           internal_polar(1,1) = internal_polar(1,1) + &
                               bem%area(i)*variables_w(i,1)*target_%coord(1,i)/&
                               field%e_0
           internal_polar(1,2) = internal_polar(1,2) + &
                               bem%area(i)*variables_w(i,2)*target_%coord(1,i)/&
                               field%e_0
           internal_polar(1,3) = internal_polar(1,3) + &
                               bem%area(i)*variables_w(i,3)*target_%coord(1,i)/&
                               field%e_0
           internal_polar(2,1) = internal_polar(2,1) + &
                               bem%area(i)*variables_w(i,1)*target_%coord(2,i)/&
                               field%e_0
           internal_polar(2,2) = internal_polar(2,2) + &
                               bem%area(i)*variables_w(i,2)*target_%coord(2,i)/&
                               field%e_0
           internal_polar(2,3) = internal_polar(2,3) + &
                               bem%area(i)*variables_w(i,3)*target_%coord(2,i)/&
                               field%e_0
           internal_polar(3,1) = internal_polar(3,1) + &
                               bem%area(i)*variables_w(i,1)*target_%coord(3,i)/&
                               field%e_0
           internal_polar(3,2) = internal_polar(3,2) + &
                               bem%area(i)*variables_w(i,2)*target_%coord(3,i)/&
                               field%e_0
           internal_polar(3,3) = internal_polar(3,3) + &
                               bem%area(i)*variables_w(i,3)*target_%coord(3,i)/&
                               field%e_0
         enddo
      endif
      polar_w = internal_polar
       
      !save data into results vector
      !1) isotropic polar --- real part
      target_%results(i_freq,1) = (dble(polar_w(1,1)) + &
                                   dble(polar_w(2,2)) + &
                                   dble(polar_w(3,3)))/ 3 
      !2) isotropic polar --- imaginary part
      target_%results(i_freq,2) = (dimag(polar_w(1,1)) + &
                                   dimag(polar_w(2,2)) + &
                                   dimag(polar_w(3,3)))/ 3 
      !3) long polar X --- real part                             
      target_%results(i_freq,3) = dble(polar_w(1,1))  
      !4) long polar Y --- real part
      target_%results(i_freq,4) = dble(polar_w(2,2))  
      !5) long polar Z --- real part
      target_%results(i_freq,5) = dble(polar_w(3,3))  
      !6) long polar X --- imag part
      target_%results(i_freq,6) = dimag(polar_w(1,1)) 
      !7) long polar Y --- imag part
      target_%results(i_freq,7) = dimag(polar_w(2,2)) 
      !8) long polar Z --- imag part
      target_%results(i_freq,8) = dimag(polar_w(3,3)) 
  
   end subroutine calculate_dynamic_polar_bem


   !> Function for calculating the induced field at a specific point
   !!    Input   : target_      -- model 
   !!    Input   : i_pol        -- index of the polarization
   !!    Input   : point_coord  -- coordinates of the point
   !!    Input   : variables_w  -- w-variables
   !!    Output  : EField       -- Electric Field
   function calculate_induced_field_at_point_bem(target_,i_pol,point_coord, &
                                                 variables_w) result(EField)
  
      implicit none
       
      !input/output variables
      class(bem_type), intent(in) :: target_
      integer, intent(in) :: i_pol
      real(dp), dimension(3), intent(in) :: point_coord
      complex(dp), dimension(target_%n_var), intent(in) :: variables_w
      complex(dp), dimension(3) :: EField
  
      !internal variables
      integer  :: i, l
      real(dp), dimension(3) :: d_IJ
      real(dp) :: distIJ
      real(dp) :: distIJ_lowest
      real(dp) :: distIJ_BEM
      real(dp) :: distIJ_grid
      logical :: is_out
       
      !we will check if the point is inside or outside bem surface
      l = 1 
      distIJ_lowest = 1.0d6
      do i = 1, target_%n_var
         d_IJ(1)   = (target_%coord(1,i)-point_coord(1))*tobohr
         d_IJ(2)   = (target_%coord(2,i)-point_coord(2))*tobohr
         d_IJ(3)   = (target_%coord(3,i)-point_coord(3))*tobohr
         distIJ = dsqrt(d_IJ(1)**2 + d_IJ(2)**2 + d_IJ(3)**2)
         if (distIJ_lowest.gt.distIJ) then
            distIJ_lowest = distIJ 
            l = i
         endif
      enddo

      !distance tesserae "l" wrt BEM center
      d_IJ(1)   = (target_%coord(1,l)-target_%center_of_mass(1))*tobohr
      d_IJ(2)   = (target_%coord(2,l)-target_%center_of_mass(2))*tobohr
      d_IJ(3)   = (target_%coord(3,l)-target_%center_of_mass(3))*tobohr    
      distIJ_BEM = dsqrt(d_IJ(1)**2 + d_IJ(2)**2 + d_IJ(3)**2) + 0.5d0*tobohr

      !distance point wrt BEM center
      d_IJ(1)   = (point_coord(1)-target_%center_of_mass(1))*tobohr
      d_IJ(2)   = (point_coord(2)-target_%center_of_mass(2))*tobohr
      d_IJ(3)   = (point_coord(3)-target_%center_of_mass(3))*tobohr    
      distIJ_grid = dsqrt(d_IJ(1)**2 + d_IJ(2)**2 + d_IJ(3)**2)

      is_out = .false.
      !check if the point is outside or inside the BEM surface
      if(distIJ_grid.gt.distIJ_bem) is_out = .true.

      EField = dcmplx(zero, zero)
      if (is_out) then
         EField(i_pol) = dcmplx(field%e_0, zero)
         do i = 1, target_%n_var
            d_IJ(1)   = (target_%coord(1,i)-point_coord(1))*tobohr
            d_IJ(2)   = (target_%coord(2,i)-point_coord(2))*tobohr
            d_IJ(3)   = (target_%coord(3,i)-point_coord(3))*tobohr
            call calculate_induced_field_q_bem(d_IJ,variables_w(i),EField)
         enddo
      endif
       
   end function calculate_induced_field_at_point_bem

end module bem_module
