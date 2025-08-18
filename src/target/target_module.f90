!> Target module
!!
!! This module contains the subroutines for the target type
!!
!! Date         : 2025
!!
module target_module
      
   !$ use omp_lib
   use output_module
   use parameters_module
   use float_manipulation_module
   use string_manipulation_module
   use memory_manager_module

   implicit none

   public target_

   type :: target_type
      integer :: n_atoms     = 0  ! number of atoms
      integer :: n_tess      = 0  ! number of BEM tesserae     
      integer :: n_var       = 0  ! number of variables   (n_var, 3)
      integer :: n_rhsre     = 1  ! number of rhs real    (n_var, 3)
      integer :: n_rhs       = 1  ! number of rhs general (n_var, 3)
      integer :: n_atomtypes = 0  ! number of atomtypes
      integer :: n_q         = 0  ! number of charges
      integer :: n_mu        = 0  ! number of dipoles
      integer :: n_mol       = 0  ! number of molecules
      integer, dimension(:), allocatable  :: n_atoms_per_molecule 
      integer, dimension(:), allocatable  :: i_mol ! index of molecule [FQ]

      !center of mass
      real(dp), dimension(:), allocatable :: center_of_mass

      !parameters of the model (wfq - wfqfmu)
      real(dp), dimension(:), allocatable :: chi    ! electronegativies
      real(dp), dimension(:), allocatable :: eta    ! chemical hardnesses
      real(dp), dimension(:), allocatable :: alpha  ! atomic polarizabilities
      real(dp), dimension(:), allocatable :: r_q    ! width gaussian charges
      real(dp), dimension(:), allocatable :: r_mu   ! width gaussian dipoles
      !if alpha is taken from tabulated data [wfqfmu]
      logical :: alpha_inside = .false.

      !heterogeneous 
      logical                              :: heterogeneous = .false. !flag
      !map atoms --> atomtypes
      integer, dimension(:), allocatable   :: map_atomtypes
      !atom indeces of the neighbours
      integer, dimension(:,:), allocatable :: neighbours    
      !number of neighboursh for each atom
      integer, dimension(:), allocatable   :: n_neighbours  

      !generic variables
      real(dp)                              :: energy        ! energy
      real(dp), dimension(:), allocatable   :: atomic_number ! atomic number
      character(4), dimension(:), allocatable :: atom_name   ! atom name

      real(dp), dimension(:,:), allocatable :: polar     ! static polar
      real(dp), dimension(:,:), allocatable :: coord     ! coordinates
      real(dp), dimension(:,:), allocatable :: matrix    ! matrix left hand side
      real(dp), dimension(:,:), allocatable :: rhs       ! right hand side
      real(dp), dimension(:,:), allocatable :: variables ! variables
      real(dp), dimension(:,:), allocatable :: results   ! results
      character(len=200), dimension(:), allocatable :: atom_type !Ag,Au,O-OW,...

      character(len=200) :: forcefield = 'fq' ! forcefield for static part
      character(len=200) :: kernel = 'ohno' ! ohno, gaussian, coulomn
      character(len=200) :: name_      ! fq, fqfmu, wfq, wfqfmu, bem, etc.

      contains

      !assign/get variables
      procedure :: assign_model_dimensions
      procedure :: assign_model_parameters
      procedure :: assign_variables_w
      procedure :: assign_variables_w_on_the_fly
      procedure :: mass_atom
      procedure :: assign_atomic_numbers
      procedure :: assign_map_atomtypes
      procedure :: get_neighbours
      procedure :: nearest_distance
      procedure :: rotate_principal_axis
      procedure :: calculate_center_of_mass
      !allocations
      procedure :: allocate_gs_matrices
      procedure :: allocate_static_field_matrices
      procedure :: allocate_dynamic_field_general
      procedure :: allocate_dynamic_field_memory
      procedure :: allocate_constant_potential_variables
      procedure :: dealloc => deallocate_target
      procedure :: deallocate_atomtypes
      procedure :: deallocate_parameters
      procedure :: deallocate_atomic_quantities
      procedure :: deallocate_gs_matrices
      procedure :: deallocate_static_field_matrices
      procedure :: deallocate_dynamic_field_general
      procedure :: deallocate_dynamic_field_memory
      procedure :: deallocate_constant_potential_variables
      !construct matrices/vectors
      procedure :: construct_static_matrix
      procedure :: construct_constant_matrix
      procedure :: construct_dynamic_matrix
      procedure :: construct_dynamic_matrix_gmres
      procedure :: construct_ground_state_rhs
      procedure :: construct_static_field_rhs
      procedure :: construct_dynamic_field_rhs
      procedure :: construct_dynamic_field_rhs_on_the_fly
      procedure :: product_matrix_vector
      procedure :: gmres_diagonal_shift
      procedure :: new_gmres_diagonal_shift
      !printing & saving files
      procedure :: print_info => print_info_target
      procedure :: print_atomtypes
      procedure :: print_coord
      procedure :: print_gs_variables        
      procedure :: print_static_field_variables        
      procedure :: print_dynamic_field_variables        
      procedure :: print_dynamic_polar
      procedure :: print_dynamic_results
      procedure :: save_intermediate_results
      procedure :: save_csv_file
      !calculate properties
      procedure :: calculate_energy
      procedure :: calculate_static_polar
      procedure :: calculate_dynamic_polar
      procedure :: calculate_cross_section
      procedure :: calculate_induced_field_at_point
      procedure :: calculate_density_at_point
      procedure :: calculate_density_at_point_q_mu_separated
   end type target_type
    
   class (target_type), pointer, save :: target_

contains

!-------------------------------------------------------------------------------
!  SUBROUTINES FOR ASSIGNING & GETTING VARIABLES
!-------------------------------------------------------------------------------

   !> Subroutine for assigning the dimensions of the specific model
   !!    In/Out : target_      -- model 
   subroutine assign_model_dimensions(target_)
  
      implicit none

      !input/output variables
      class(target_type), intent(inout) :: target_
       
      if(target_%name_.ne.'fq'        .and. &
         target_%name_.ne.'fqfmu'     .and. &
         target_%name_.ne.'wfq'       .and. &
         target_%name_.ne.'wfqfmu'    .and. &
         target_%name_.ne.'bem') &
         call out_%error("Target name: "//trim(target_%name_)//" not recognised")
  
   end subroutine assign_model_dimensions


   !> Subroutine for assigning the parameters of the specific model
   !!    In/Out : target_      -- model 
   subroutine assign_model_parameters(target_)
  
      implicit none

      !input/output variables
      class(target_type), intent(inout) :: target_
       
      if(target_%name_.ne.'wfq'       .and. &
         target_%name_.ne.'wfqfmu'    .and. &
         target_%name_.ne.'bem') &
         call out_%error("Target name: "//trim(target_%name_)//" not recognised")
  
   end subroutine assign_model_parameters


   !> Subroutine for assigning the dynamic variables
   !!   Default : variables_w = rhs_w. This changes for heterogeneous
   !!    Input  : target_      -- model 
   !!    In/Out : rhs_w        -- dynamic RHS
   !!    Output : variables_w  -- dynamic variables
   subroutine assign_variables_w(target_, rhs_w, variables_w)
  
      implicit none
       
      !input/output variables
      class(target_type), intent(inout) :: target_
      complex(dp), dimension(target_%n_var,3), intent(inout) :: rhs_w
      complex(dp), dimension(target_%n_var,3), intent(out)   :: variables_w
  
      if(target_%name_.ne.'wfq'       .and. &
         target_%name_.ne.'wfqfmu'    .and. &
         target_%name_.ne.'bem') &
         call out_%error("Target name: "//trim(target_%name_)//" not recognised")

      !variables_w = rhs_w !default
      call zcopy(target_%n_var*3, rhs_w, 1, variables_w, 1)
  
   end subroutine assign_variables_w
   

   !> Subroutine for assigning the dynamic variables on the fly algorithm
   !!   Default : variables_w = rhs_w. This changes for heterogeneous
   !!    Input  : target_      -- model 
   !!    Input  : i_freq       -- index of the frequency
   !!    In/Out : rhs_w        -- dynamic RHS
   !!    Output : variables_w  -- dynamic variables
   subroutine assign_variables_w_on_the_fly(target_, i_freq, rhs_w, variables_w)
  
      implicit none
       
      !input/output variables
      class(target_type), intent(inout) :: target_
      integer, intent(in)               :: i_freq
      complex(dp), dimension(target_%n_var,3), intent(inout) :: rhs_w
      complex(dp), dimension(target_%n_var,3), intent(out)   :: variables_w
  
      if(target_%name_.ne.'wfq'       .and. &
         target_%name_.ne.'wfqfmu'    .and. &
         target_%name_.ne.'bem') &
         call out_%error("Target name: "//trim(target_%name_)//" not recognised")
          
      !variables_w = rhs_w !default
      call zcopy(target_%n_var*3, rhs_w, 1, variables_w, 1)
  
   end subroutine assign_variables_w_on_the_fly


   !> Function for assigning the atomic masses
   !!    Input  : target_      -- model 
   !!    Input  : i            -- index of the atom
   !!    Output : mass         -- mass of the atom
   function mass_atom(target_,i) result(mass)
     
      implicit none
  
      !input/output variables
      class(target_type), intent(in) :: target_
      integer, intent(in)            :: i
      real(dp)                       :: mass
  
      mass = zero !default
      if(is_equal(target_%atomic_number(i),1.0d0)) then !H
         mass = 1.00794d0
      else if(is_equal(target_%atomic_number(i),2.0d0)) then !He
         mass = 4.00260d0
      else if(is_equal(target_%atomic_number(i),3.0d0)) then !Li 
         mass = 6.941d0
      else if(is_equal(target_%atomic_number(i),4.0d0)) then !Be
         mass = 9.01218d0
      else if(is_equal(target_%atomic_number(i),5.0d0)) then !B
         mass = 10.81d0
      else if(is_equal(target_%atomic_number(i),6.0d0)) then !C
         mass = 12.0107d0
      else if(is_equal(target_%atomic_number(i),7.0d0)) then !N
         mass = 14.0067d0 
      else if(is_equal(target_%atomic_number(i),8.0d0)) then !O
         mass = 15.9990d0
      else if(is_equal(target_%atomic_number(i),9.0d0)) then !F
         mass = 18.998403d0
      else if(is_equal(target_%atomic_number(i),10.0d0)) then !Ne
         mass = 20.179d0
      else if(is_equal(target_%atomic_number(i),11.0d0)) then !Na
         mass = 22.989770d0
      else if(is_equal(target_%atomic_number(i),12.0d0)) then !Mg
         mass = 24.305d0
      else if(is_equal(target_%atomic_number(i),13.0d0)) then !Al
         mass = 26.98154d0
      else if(is_equal(target_%atomic_number(i),14.0d0)) then !Si
         mass = 28.0855d0
      else if(is_equal(target_%atomic_number(i),15.0d0)) then !P
         mass = 30.97376d0
      else if(is_equal(target_%atomic_number(i),16.0d0)) then !S
         mass = 32.06d0
      else if(is_equal(target_%atomic_number(i),17.0d0)) then !Cl
         mass = 35.453d0
      else if(is_equal(target_%atomic_number(i),18.0d0)) then !Ar
         mass = 39.948d0
      else if(is_equal(target_%atomic_number(i),29.0d0)) then !Cu
         mass = 63.546d0
      else if(is_equal(target_%atomic_number(i),47.0d0)) then !Ag
         mass = 107.8682d0
      else if(is_equal(target_%atomic_number(i),79.0d0)) then !Au
         mass = 196.9665d0
      else 
         call out_%error('Atomic number not recognized for atomic mass')
      endif
  
   end function mass_atom


   !> Subroutine for assigning the atomic numbers to the atoms
   !!    In/Out : target_      -- model 
   !!    Input  : AtomName     -- the name of the atoms
   subroutine assign_atomic_numbers(target_, atomname) 
  
      implicit none
  
      !input/output variables
      class(target_type), intent(inout) :: target_
      character(len=*), dimension(target_%n_atoms), intent(in) :: atomname
  
      !internal variables
      integer :: i
  
      do i = 1, target_%n_atoms
         target_%atomic_number(i) = assign_atomic_number(lower(atomname(i)))
      enddo
  
   end subroutine assign_atomic_numbers


   !> Function to assign the atomic number based on the atom name
   !!    Input  : AtomName     -- the name of the atoms
   real(dp) function assign_atomic_number(atomname) result(atmnumber)
   
      implicit none
   
      !input/output variables
      character(len=*), intent(in)  :: atomname
   
      atmnumber = 0.0d0
      if(trim(atomname).eq.'h') then
         atmnumber = 1.0d0
      else if(trim(atomname).eq.'hw') then
         atmnumber = 1.0d0
      else if(trim(atomname).eq.'h-hw') then
         atmnumber = 1.0d0
      else if(trim(atomname).eq.'he') then
         atmnumber = 2.0d0
      else if(trim(atomname).eq.'li') then
         atmnumber = 3.0d0
      else if(trim(atomname).eq.'be') then
         atmnumber = 4.0d0
      else if(trim(atomname).eq.'b') then
         atmnumber = 5.0d0
      else if(trim(atomname).eq.'c') then
         atmnumber = 6.0d0
      else if(trim(atomname).eq.'n') then
         atmnumber = 7.0d0
      else if(trim(atomname).eq.'o') then
         atmnumber = 8.0d0
      else if(trim(atomname).eq.'ow') then
         atmnumber = 8.0d0
      else if(trim(atomname).eq.'o-ow') then
         atmnumber = 8.0d0
      else if(trim(atomname).eq.'f') then
         atmnumber = 9.0d0
      else if(trim(atomname).eq.'ne') then
         atmnumber = 10.0d0
      else if(trim(atomname).eq.'na') then
         atmnumber = 11.0d0
      else if(trim(atomname).eq.'mg') then
         atmnumber = 12.0d0
      else if(trim(atomname).eq.'al') then
         atmnumber = 13.0d0
      else if(trim(atomname).eq.'si') then
         atmnumber = 14.0d0
      else if(trim(atomname).eq.'p') then
         atmnumber = 15.0d0
      else if(trim(atomname).eq.'s') then
         atmnumber = 16.0d0
      else if(trim(atomname).eq.'cl') then
         atmnumber = 17.0d0
      else if(trim(atomname).eq.'ar') then
         atmnumber = 18.0d0
      else if(trim(atomname).eq.'ag') then
         atmnumber = 47.0d0
      else if(trim(atomname).eq.'cu') then
         atmnumber = 29.0d0
      else if(trim(atomname).eq.'au') then
         atmnumber = 79.0d0
      else 
         call out_%error('atomname: '//trim(atomname)//&
                         ' not recognized for atomic symbol')
      endif
   
   end function assign_atomic_number


   !> Subroutine for assigning the map of the atomtypes: in this way, we store
   !! less information, by mapping the atomic parameters to the atomtypes.
   !!    In/Out : target_      -- model 
   !!    Input  : AtomName     -- the name of the atoms
   subroutine assign_map_atomtypes(target_, AtomName) 
  
      implicit none
  
      !input/output variables
      class(target_type), intent(inout) :: target_
      character(len=*), dimension(target_%n_atoms), intent(in)  :: AtomName
  
      !internal variables
      integer :: i, j
      logical :: found_atomname
  
      do i = 1, target_%n_atoms
         found_atomname = .false.
         do j = 1, target_%n_atomtypes
            if(trim(atomname(i)).eq.target_%atom_type(j)) then
              target_%map_atomtypes(i) = j
              found_atomname = .true.
              exit
            endif
            if(j.eq.target_%n_atomtypes.and..not.found_atomname) &
               call out_%error("No parameters given for atomtype: "// &
                               trim(AtomName(i)))
         enddo
      enddo
  
   end subroutine assign_map_atomtypes


   !> Subroutine for getting the neighbours of each atom
   !! This is needed for heterostructures to define the polar
   !!    In/Out : target_      -- model 
   subroutine get_neighbours(target_) 
  
      implicit none
  
      !input/output variables
      class(target_type), intent(inout) :: target_
  
      !internal variables
      integer :: i, j
      integer :: n_tot
      integer, dimension(:), allocatable :: factor !tmp vector
      real(dp) :: distij
  
      call mem_man%alloc(target_%neighbours,target_%n_atomtypes, &
                         target_%n_atoms,"target_%neighbours")
      call mem_man%alloc(target_%n_neighbours,target_%n_atoms, &
                         "target_%n_neighbours")
      call mem_man%alloc(factor, target_%n_atomtypes, "factor")
  
      do i = 1, target_%n_atoms
         n_tot  = 1
         factor = 0
         factor(target_%map_atomtypes(i)) = 1
         do j = 1, target_%n_atoms
            if(i.ne.j) then
               distij = sqrt((target_%coord(1,j)-target_%coord(1,i))**2 + &
                             (target_%coord(2,j)-target_%coord(2,i))**2 + &
                             (target_%coord(3,j)-target_%coord(3,i))**2 )*tobohr
               if(distij.le.target_%nearest_distance(i,j)) then !if neighbours
                  n_tot = n_tot + 1
                  if (is_equal(target_%atomic_number(i), &
                               target_%atomic_number(j))) then
                      factor(target_%map_atomtypes(i)) = &
                      factor(target_%map_atomtypes(i)) + 1 !same type
                  else
                      factor(target_%map_atomtypes(j)) = &
                      factor(target_%map_atomtypes(j)) + 1 !different type
                  endif
               endif
            endif
         enddo
         target_%neighbours(:,i) = factor
         target_%n_neighbours(i) = n_tot
      enddo
  
      call mem_man%dealloc(factor,"factor")
  
   end subroutine get_neighbours


   !> Function for defining the nearest neighbour distance
   !! Needed for Fermi function - w models
   !!    Input  : target_      -- model 
   !!    Input  : i,j          -- indices of the atoms
   !!    Output : rij0         -- neighbour distance
   real(dp) function nearest_distance(target_, i,j) result(rij0)
     
      implicit none
               
      !input/output variables
      class(target_type), intent(in) :: target_
      integer, intent(in)            :: i,j 
      
      if( (is_equal(target_%atomic_number(i),8.0d0).and. & !O  -- this is for water
           is_equal(target_%atomic_number(j),1.0d0)).or. & !H
          (is_equal(target_%atomic_number(i),1.0d0).and. & !H
           is_equal(target_%atomic_number(j),8.0d0)).or. & !H
          (is_equal(target_%atomic_number(i),8.0d0).and. & !O
           is_equal(target_%atomic_number(j),8.0d0)).or. & !O
          (is_equal(target_%atomic_number(i),1.0d0).and. & !H
           is_equal(target_%atomic_number(j),1.0d0))) then !H 
         rij0 = rOH0
      else if(is_equal(target_%atomic_number(i),47.0d0).and. & !Ag
              is_equal(target_%atomic_number(j),47.0d0)) then  !Ag
         rij0 = rAg0
      else if(is_equal(target_%atomic_number(i),79.0d0).and. & !Au
              is_equal(target_%atomic_number(j),79.0d0)) then  !Au
         rij0 = rAu0
      else if(is_equal(target_%atomic_number(i),11.0d0).and. & !Na
              is_equal(target_%atomic_number(j),11.0d0)) then  !Na
         rij0 = rNa0
      else if(is_equal(target_%atomic_number(i),6.0d0).and.  & !C 
              is_equal(target_%atomic_number(j),6.0d0)) then   !C
         rij0 = rC0
      else if(is_equal(target_%atomic_number(i),13.0d0).and. & !Al
              is_equal(target_%atomic_number(j),13.0d0)) then  !Al
         rij0 = rAl0  
      else if(is_equal(target_%atomic_number(i),29.0d0).and. & !Cu 
              is_equal(target_%atomic_number(j),29.0d0)) then  !Cu
         rij0 = rCu0
      else if( (is_equal(target_%atomic_number(i),47.0d0).and. & !Ag -- alloys
                is_equal(target_%atomic_number(j),79.0d0)).or. & !Au
               (is_equal(target_%atomic_number(i),79.0d0).and. & !Au
                is_equal(target_%atomic_number(j),47.0d0))) then !Ag
         rij0 = ( rAg0 + rAu0 ) / two
      else !not defined
         rij0 = zero
         call out_%error('I do not recognize atomic couples in &
                         &nearest distance')
      endif
       
   end function nearest_distance


   !> Function for assigning atom name given the atomic number
   !!    Input  : AtmNumber    -- Atomic Number
   !!    Output : AtomName     -- Atom Name
   character(2) function assign_atom_name(AtmNumber) result(AtomName)
 
      implicit none
 
      !input/output variables
      real(dp), intent(in) :: AtmNumber
 
      if(is_equal(AtmNumber,1.0d0)) then
         AtomName = 'H'
      else if(is_equal(AtmNumber,2.0d0)) then
         AtomName = 'He'
      else if(is_equal(AtmNumber,3.0d0)) then
         AtomName = 'Li'
      else if(is_equal(AtmNumber,4.0d0)) then
         AtomName = 'Be'
      else if(is_equal(AtmNumber,5.0d0)) then
         AtomName = 'B'
      else if(is_equal(AtmNumber,6.0d0)) then
         AtomName = 'C'
      else if(is_equal(AtmNumber,7.0d0)) then
         AtomName = 'N'
      else if(is_equal(AtmNumber,8.0d0)) then
         AtomName = 'O'
      else if(is_equal(AtmNumber,9.0d0)) then
         AtomName = 'F'
      else if(is_equal(AtmNumber,10.0d0)) then
         AtomName = 'Ne'
      else if(is_equal(AtmNumber,11.0d0)) then
         AtomName = 'Na'
      else if(is_equal(AtmNumber,12.0d0)) then
         AtomName = 'Mg'
      else if(is_equal(AtmNumber,13.0d0)) then
         AtomName = 'Al'
      else if(is_equal(AtmNumber,14.0d0)) then
         AtomName = 'Si'
      else if(is_equal(AtmNumber,15.0d0)) then
         AtomName = 'P'
      else if(is_equal(AtmNumber,16.0d0)) then
         AtomName = 'S'
      else if(is_equal(AtmNumber,17.0d0)) then
         AtomName = 'Cl'
      else if(is_equal(AtmNumber,18.0d0)) then
         AtomName = 'Ar'
      else if(is_equal(AtmNumber,47.0d0)) then
         AtomName = 'Ag'
      else if(is_equal(AtmNumber,79.0d0)) then
         AtomName = 'Au'
      else 
         call out_%error("AtomName "//AtomName//" not currently implemented.")
      endif
 
   end function assign_atom_name


   !> Function for assigning atom radius given the atom name
   !!    Input  : AtomName     -- Atom Name
   !!    Output : Radius       -- Radius in Angstrom [vdW radii]
   function assign_atomic_radius(AtomName) result(radius)
  
      implicit none
  
      !input/output variables
      character(len=2), intent(in)  :: AtomName
      real(dp) :: radius
  
      radius = 0.0d0 !default
      If(trim(AtomName).eq.'H') then
         radius = 1.200d0
      else If(trim(AtomName).eq.'He') then
         radius = 1.400d0
      else If(trim(AtomName).eq.'Li') then
         radius = 1.820d0
      else If(trim(AtomName).eq.'Be') then
         radius = 1.530d0
      else If(trim(AtomName).eq.'B') then
         radius = 1.920d0
      else If(trim(AtomName).eq.'C') then
         radius = 1.700d0
      else If(trim(AtomName).eq.'N') then
         radius = 1.550d0
      else If(trim(AtomName).eq.'O') then
         radius = 1.520d0
      else If(trim(AtomName).eq.'F') then
         radius = 1.470d0
      else If(trim(AtomName).eq.'Ne') then
         radius = 1.540d0
      else If(trim(AtomName).eq.'Na') then
         radius = 2.270d0
      else If(trim(AtomName).eq.'Mg') then
         radius = 1.730d0
      else If(trim(AtomName).eq.'Al') then
         radius = 1.840d0
      else If(trim(AtomName).eq.'Si') then
         radius = 2.100d0
      else If(trim(AtomName).eq.'P') then
         radius = 1.800d0
      else If(trim(AtomName).eq.'S') then
         radius = 1.800d0
      else If(trim(AtomName).eq.'Cl') then
         radius = 1.750d0
      else If(trim(AtomName).eq.'Ar') then
         radius = 1.880d0
      else If(trim(AtomName).eq.'Ag') then
         radius = 1.720d0
      else If(trim(AtomName).eq.'Au') then
         radius = 1.660d0
      endif
 
   end function assign_atomic_radius


   !> Function for calculating the center of mass
   !!    Input  : target_      -- target_ type
   !!    Output : com          -- center of mass
   function calculate_center_of_mass(target_) result(com)
  
      implicit none
  
      !input/output variables
      class(target_type), intent(in) :: target_
      real(dp), dimension(3) :: com

      !internal variables
      integer  :: i
      real(dp) :: total_mass
  
      call array_clear(3,com)
      if (target_%name_.ne.'bem') then
         !calculation of total mass
         total_mass = zero
         !$omp parallel do reduction(+:total_mass)
         do i = 1, target_%n_atoms
            total_mass = total_mass + target_%mass_atom(i)
         enddo
         !$omp end parallel do
  
         !center of mass = sum_i m(i) * coordinate(i) / total_mass 
         !$omp parallel do reduction(+:com)
         do i = 1,target_%n_atoms
            com(1) = com(1) + target_%mass_atom(i)*target_%coord(1,i)/total_mass
            com(2) = com(2) + target_%mass_atom(i)*target_%coord(2,i)/total_mass
            com(3) = com(3) + target_%mass_atom(i)*target_%coord(3,i)/total_mass
         enddo
         !$omp end parallel do
      else
         !center of mass = sum_i coordinate(i) / target_%n_var
         !$omp parallel do reduction(+:com)
         do i = 1,target_%n_var
            com(1) = com(1) + target_%coord(1,i)/target_%n_var
            com(2) = com(2) + target_%coord(2,i)/target_%n_var
            com(3) = com(3) + target_%coord(3,i)/target_%n_var
         enddo
         !$omp end parallel do
      endif
 
   end function calculate_center_of_mass

   !> Subroutine for rotating the system coordinates to principal axis frame
   !!    In/Out  : target_      -- model 
   !!    Input   : save_info    -- save info file
   subroutine rotate_principal_axis(target_, save_info)
  
      implicit none
  
      !input/output variables
      class(target_type), intent(inout) :: target_
      logical, intent(in) :: save_info
  
      !internal variables
      integer :: i
      integer :: info
      real(dp), dimension(:,:),allocatable :: internal_coordinates
      real(dp), dimension(3)               :: eigenvalues
      real(dp), dimension(3,3)             :: inertial_moment
      integer :: lipiv
      integer, dimension(:), allocatable    :: ipiv
      integer, dimension(1)                 :: ipiv_dummy
      integer :: lwork
      real(dp), dimension(1)                :: work_dummy
      real(dp), dimension(:),   allocatable :: work
      integer  :: iost
  
      if(target_%name_.eq.'bem') return ! no rotations for BEM

      call mem_man%alloc(target_%center_of_mass, 3, "center_of_mass")
      target_%center_of_mass = target_%calculate_center_of_mass()
  
      !put the center of mass in the origin
      do i = 1,target_%n_atoms 
         target_%coord(1,i) = target_%coord(1,i) - target_%center_of_mass(1)
         target_%coord(2,i) = target_%coord(2,i) - target_%center_of_mass(2)
         target_%coord(3,i) = target_%coord(3,i) - target_%center_of_mass(3)
      enddo
      call mem_man%dealloc(target_%center_of_mass, "center_of_mass")
  
      !calculation of the intertia moment
      call array_clear(9,inertial_moment)
      !$omp parallel do reduction(+:inertial_moment)
      do i = 1, target_%n_atoms
         inertial_moment(1,1) = inertial_moment(1,1) + target_%mass_atom(i)* &
                                (target_%coord(2,i)**2+target_%coord(3,i)**2)
         inertial_moment(2,2) = inertial_moment(2,2) + target_%mass_atom(i)* &
                                (target_%coord(1,i)**2+target_%coord(3,i)**2)
         inertial_moment(3,3) = inertial_moment(3,3) + target_%mass_atom(i)* &
                                (target_%coord(2,i)**2+target_%coord(1,i)**2)
         inertial_moment(1,2) = inertial_moment(1,2) - target_%mass_atom(i)* &
                                (target_%coord(1,i)*target_%coord(2,i))
         inertial_moment(1,3) = inertial_moment(1,3) - target_%mass_atom(i)* &
                                (target_%coord(1,i)*target_%coord(3,i))
         inertial_moment(2,3) = inertial_moment(2,3) - target_%mass_atom(i)* &
                                (target_%coord(2,i)*target_%coord(3,i))
      enddo
      !$omp end parallel do
      inertial_moment(2,1)=inertial_moment(1,2)
      inertial_moment(3,1)=inertial_moment(1,3)
      inertial_moment(3,2)=inertial_moment(2,3)
       
      if(out_%ivrb.ge.3) call out_%print_matrix('Inertial momentum', &
                                                inertial_moment,3,3)
       
      !diagonalization of the inertial moment
      call dsyevd('v',             &
                  'u',             &
                  3,               &
                  inertial_moment, &
                  3,               &
                  eigenvalues,     &
                  work_dummy,      &
                  -1,              &
                  ipiv_dummy,      &
                  -1,              &
                  info)
      if(info.ne.0) call out_%error("INFO.ne.0 in DSYEVD (Inertial momentum)")

      lwork = nint(work_dummy(1))
      call mem_man%alloc(work, lwork, "work_rotate_principal_axis")
      lipiv = ipiv_dummy(1)
      call mem_man%alloc(ipiv, lipiv, "ipiv_rotate_principal_axis")
     
      call dsyevd('v',             &
                  'u',             &
                  3,               &
                  inertial_moment, &
                  3,               &
                  eigenvalues,     &
                  work,            &
                  lwork,           &
                  ipiv,            &
                  lipiv,           &
                  info)
      if(info.ne.0) call out_%error("INFO.ne.0 in DSYEVD (Inertial momentum)")
      call mem_man%dealloc(work, "work_rotate_principal_axis")
      call mem_man%dealloc(ipiv, "ipiv_rotate_principal_axis")
      
      if(out_%ivrb.ge.3) then
         Write(out_%iunit,'(/a/)') '  EigenValues Inertial Momentum' 
         write(out_%iunit,'(3x,i1,4x,e13.6)') 1, eigenvalues(1)
         write(out_%iunit,'(3x,i1,4x,e13.6)') 2, eigenvalues(2)
         write(out_%iunit,'(3x,i1,4x,e13.6)') 3, eigenvalues(3)
         call out_%print_matrix('EigenVector Inertial Momentum', &
                                 inertial_moment,3,3)
      endif
  
      call mem_man%alloc(internal_coordinates, 3,target_%n_atoms, &
                         "internal_coordinates")
      !calculation of new coordinates stored temporarily 
      call dgemm ('t', 'n',              &
                   3,                    &
                   target_%n_atoms,      &
                   3,                    &
                   one,                  &
                   inertial_moment,      &
                   3,                    &
                   target_%coord,        &
                   3,                    &
                   zero,                 &
                   internal_coordinates, &
                   3)
      !$omp parallel do private(i)
      do i = 1, target_%n_atoms
         target_%coord(1,i) = internal_coordinates(1,i) 
         target_%coord(2,i) = internal_coordinates(2,i) 
         target_%coord(3,i) = internal_coordinates(3,i) 
      enddo
      !$omp end parallel do
  
      call mem_man%dealloc(internal_coordinates,"internal_coordinates")
     
      call target_%print_coord("Rotated") !printing rotated coordinates
  
      !save the information in the save_info file if requested
      if(save_info) then 
         if(out_%ivrb.ge.3) then 
            write(out_%iunit,'(1x,a)') "Rotated Geometry saved in : "//&
                                        out_%info_file
            write(out_%iunit,out_%sticks) 
         endif
         !put coordinates in the info file
         open(unit=out_%unit_info,file=out_%info_file,status="old", &
              iostat=iost,position='append')
            write(out_%unit_info,'(a)') 'Rotated Geometry (Angstrom)'
            do i = 1, target_%n_atoms
               write(out_%unit_info,'(f4.1,2x,3(1x,f25.16))')  &
                                target_%atomic_number(i), &
                                target_%coord(1,i),       &
                                target_%coord(2,i),       &
                                target_%coord(3,i)
            enddo 
         close(out_%unit_info)
      endif
  
   end subroutine rotate_principal_axis

!-------------------------------------------------------------------------------
!  SUBROUTINES FOR ALLOCATIONS
!-------------------------------------------------------------------------------

   !> Subroutine for allocating the GS matrices
   !!    In/Out  : target_             -- model 
   subroutine allocate_gs_matrices(target_)

      implicit none
      !input/output variables
      class(target_type), intent(inout) :: target_

      if(target_%name_.ne.'fq'        .and. &
         target_%name_.ne.'fqfmu'     .and. &
         target_%name_.ne.'wfq'       .and. &
         target_%name_.ne.'wfqfmu'    .and. &
         target_%name_.ne.'bem') &
         call out_%error("Target name: "//trim(target_%name_)// &
                         " not recognised")
     
   end subroutine allocate_gs_matrices


   !> Subroutine for allocating the matrices for static field
   !!    In/Out  : target_             -- model 
   subroutine allocate_static_field_matrices(target_)
  
      implicit none
       
      !input/output variables
      class(target_type), intent(inout) :: target_
  
      if(target_%name_.ne.'fq'        .and. &
         target_%name_.ne.'fqfmu'     .and. &
         target_%name_.ne.'wfq'       .and. &
         target_%name_.ne.'wfqfmu'    .and. &
         target_%name_.ne.'bem') &
         call out_%error("Target name: "//trim(target_%name_)//" not recognised")
       
   end subroutine allocate_static_field_matrices


   !> Subroutine for allocating the matrices for dynamic field
   !! for inversion, iterative, on the fly
   !!    In/Out  : target_             -- model 
   !!    In/Out  : variables_w         -- w-variables
   !!    In/Out  : rhs_w               -- dynamic RHS
   !!    In/Out  : polar_w             -- dynamic polar
   subroutine allocate_dynamic_field_general(target_,variables_w,rhs_w,polar_w)
  
      use field_module
  
      implicit none
  
      !input/output variables
      class(target_type)    :: target_
      complex(dp), dimension(:,:), allocatable :: variables_w
      complex(dp), dimension(:,:), allocatable :: rhs_w
      complex(dp), dimension(:,:), allocatable :: polar_w
  
      if(target_%name_.ne.'wfq'       .and. &
         target_%name_.ne.'wfqfmu'    .and. &
         target_%name_.ne.'bem') &
         call out_%error("Target name: "//trim(target_%name_)//&
                         " not recognised")
       
      target_%n_rhs = 3
      call mem_man%alloc(variables_w, target_%n_var,3, "variables_w")
      call mem_man%alloc(rhs_w, target_%n_var, 3, "rhs_w") 
      call mem_man%alloc(polar_w, 3,3, "polar_w")
      call mem_man%alloc(target_%results, field%n_freq, 20, "target_%results")
       
   end subroutine allocate_dynamic_field_general


   !> Subroutine for allocating dynamic field variables for memory inv/iter
   !!    In/Out  : target_             -- model 
   !!    In/Out  : matrix_constant     -- part of matrix_w that is constant
   !!    In/Out  : matrix_w            -- dynamic Matrix
   subroutine allocate_dynamic_field_memory(target_,matrix_constant,matrix_w)
  
      implicit none
       
      !input/output variables
      class(target_type)    :: target_
      real(dp), dimension(:,:), allocatable    :: matrix_constant
      complex(dp), dimension(:,:), allocatable :: matrix_w
  
      if(target_%name_.ne.'wfq'       .and. &
         target_%name_.ne.'wfqfmu'    .and. &
         target_%name_.ne.'bem') &
         call out_%error("Target name: "//trim(target_%name_)//&
                         " not recognised")
       
   end subroutine allocate_dynamic_field_memory


   !> Subroutine for allocating constant potential variables for BEM charge
   !! constraints
   !!    In/Out  : target_             -- model 
   !!    In/Out  : variables_1         -- new w-variables
   !!    In/Out  : rhs_1               -- rhs = 1
   subroutine allocate_constant_potential_variables(target_,variables_1,rhs_1)
  
      implicit none
  
      !input/output variables
      class(target_type)    :: target_
      complex(dp), dimension(:,:), allocatable :: rhs_1
      complex(dp), dimension(:,:), allocatable :: variables_1
  
      if(target_%name_.ne.'bem') &
         call out_%error("Target name: "//trim(target_%name_)// &
                         " not recognised")
       
      call mem_man%alloc(variables_1, target_%n_var,3, "variables_1")
  
      call mem_man%alloc(rhs_1, target_%n_var, 3, "rhs_1")
      rhs_1 = cmplx(one,zero,kind=dp)
      
   end subroutine allocate_constant_potential_variables


   !> Subroutine for deallocating the target_ variables
   !!    In/Out : target_      -- target_type
   subroutine deallocate_target(target_)

      implicit none
  
      !input/output variables
      class(target_type), intent(inout) :: target_

      if(target_%name_.ne.'fq'        .and. &
         target_%name_.ne.'fqfmu'     .and. &
         target_%name_.ne.'wfq'       .and. &
         target_%name_.ne.'wfqfmu'    .and. &
         target_%name_.ne.'bem') &
         call out_%error("Target name: "//trim(target_%name_)//&
                         " not recognised")
 
   end subroutine deallocate_target


   !> Subroutine for deallocating the atomtypes variables
   !!    In/Out : target_      -- target_type
   subroutine deallocate_atomtypes(target_)

      implicit none
  
      !input/output variables
      class(target_type), intent(inout) :: target_

      call mem_man%dealloc(target_%atom_type, 'target_%atom_type')
      call mem_man%dealloc(target_%chi, 'target_%chi')
      call mem_man%dealloc(target_%eta, 'target_%eta')
      call mem_man%dealloc(target_%alpha, 'target_%alpha')
      call mem_man%dealloc(target_%r_q, 'target_%r_q')
      call mem_man%dealloc(target_%r_mu, 'target_%r_mu')
 
   end subroutine deallocate_atomtypes


   !> Subroutine for deallocating the parameters variables
   !!    In/Out : target_      -- target_type
   subroutine deallocate_parameters(target_)

      implicit none
  
      !input/output variables
      class(target_type), intent(inout) :: target_

      call mem_man%dealloc(parameters%tau, "parameters%tau")
      call mem_man%dealloc(parameters%sigma_0, "parameters%sigma0")
      call mem_man%dealloc(parameters%scaling, "parameters%scaling")
      call mem_man%dealloc(parameters%A_ij, "parameters%A_ij")
      call mem_man%dealloc(parameters%fermi_d, "parameters%fermi_d")
      call mem_man%dealloc(parameters%fermi_s, "parameters%fermi_s")
      call mem_man%dealloc(parameters%fermi_energy, "parameters%fermi_energy")
      call mem_man%dealloc(parameters%permittivity_type, & 
                           "parameters%permittivity_type")
      call mem_man%dealloc(parameters%wfqfmu_file, "parameters%wfqfmu_file")
      if(allocated(parameters%density)) &
         call mem_man%dealloc(parameters%density, "parameters%density")
      if(target_%name_.eq.'wfqfmu') &
         call mem_man%dealloc(parameters%alpha_w, "parameters%alpha_w")
 
   end subroutine deallocate_parameters


   !> Subroutine for deallocating the atomic quantities
   !!    In/Out : target_      -- target_type
   subroutine deallocate_atomic_quantities(target_)

      implicit none
  
      !input/output variables
      class(target_type), intent(inout) :: target_

      call mem_man%dealloc(target_%i_mol, "target_%i_mol")
      call mem_man%dealloc(target_%coord, "target_%coord")
      call mem_man%dealloc(target_%atomic_number, "target_%atomic_numbers")
      call mem_man%dealloc(target_%map_atomtypes, "target_%map_atomtypes")
      call mem_man%dealloc(target_%n_atoms_per_molecule, &
                           "target_%n_atoms_per_molecule")
      call mem_man%dealloc(target_%neighbours,"target_%neighbours")
      call mem_man%dealloc(target_%n_neighbours,"target_%n_neighbours")
 
   end subroutine deallocate_atomic_quantities


   !> Subroutine for allocating the GS matrices
   !!    In/Out  : target_             -- model 
   subroutine deallocate_gs_matrices(target_)

      implicit none
      !input/output variables
      class(target_type), intent(inout) :: target_

      if(target_%name_.ne.'fq'        .and. &
         target_%name_.ne.'fqfmu'     .and. &
         target_%name_.ne.'wfq'       .and. &
         target_%name_.ne.'wfqfmu'    .and. &
         target_%name_.ne.'bem') &
         call out_%error("Target name: "//trim(target_%name_)// &
                         " not recognised")
     
   end subroutine deallocate_gs_matrices


   !> Subroutine for allocating the matrices for static field
   !!    In/Out  : target_             -- model 
   subroutine deallocate_static_field_matrices(target_)
  
      implicit none
       
      !input/output variables
      class(target_type), intent(inout) :: target_
  
      if(target_%name_.ne.'fq'        .and. &
         target_%name_.ne.'fqfmu'     .and. &
         target_%name_.ne.'wfq'       .and. &
         target_%name_.ne.'wfqfmu'    .and. &
         target_%name_.ne.'bem') &
         call out_%error("Target name: "//trim(target_%name_)//" not recognised")
       
   end subroutine deallocate_static_field_matrices


   !> Subroutine for allocating the matrices for dynamic field
   !! for inversion, iterative, on the fly
   !!    In/Out  : target_             -- model 
   !!    In/Out  : variables_w         -- w-variables
   !!    In/Out  : rhs_w               -- dynamic RHS
   !!    In/Out  : polar_w             -- dynamic polar
   subroutine deallocate_dynamic_field_general(target_,variables_w,rhs_w,polar_w)
  
      use field_module
  
      implicit none
  
      !input/output variables
      class(target_type)    :: target_
      complex(dp), dimension(:,:), allocatable :: variables_w
      complex(dp), dimension(:,:), allocatable :: rhs_w
      complex(dp), dimension(:,:), allocatable :: polar_w
  
      if(target_%name_.ne.'wfq'       .and. &
         target_%name_.ne.'wfqfmu'    .and. &
         target_%name_.ne.'bem') &
         call out_%error("Target name: "//trim(target_%name_)//&
                         " not recognised")
       
      call mem_man%dealloc(variables_w, "variables_w")
      call mem_man%dealloc(rhs_w, "rhs_w") 
      call mem_man%dealloc(polar_w, "polar_w")
      call mem_man%dealloc(target_%results, "target_%results")
       
   end subroutine deallocate_dynamic_field_general


   !> Subroutine for allocating dynamic field variables for memory inv/iter
   !!    In/Out  : target_             -- model 
   !!    In/Out  : matrix_constant     -- part of matrix_w that is constant
   !!    In/Out  : matrix_w            -- dynamic Matrix
   subroutine deallocate_dynamic_field_memory(target_,matrix_w,matrix_constant)
  
      implicit none
       
      !input/output variables
      class(target_type)    :: target_
      complex(dp), dimension(:,:), allocatable :: matrix_w
      real(dp), dimension(:,:), allocatable, optional :: matrix_constant
  
      if(target_%name_.ne.'wfq'       .and. &
         target_%name_.ne.'wfqfmu'    .and. &
         target_%name_.ne.'bem') &
         call out_%error("Target name: "//trim(target_%name_)//&
                         " not recognised")
       
   end subroutine deallocate_dynamic_field_memory


   !> Subroutine for allocating constant potential variables for BEM charge
   !! constraints
   !!    In/Out  : target_             -- model 
   !!    In/Out  : variables_1         -- new w-variables
   !!    In/Out  : rhs_1               -- rhs = 1
   subroutine deallocate_constant_potential_variables(target_,variables_1,rhs_1)
  
      implicit none
  
      !input/output variables
      class(target_type)    :: target_
      complex(dp), dimension(:,:), allocatable :: rhs_1
      complex(dp), dimension(:,:), allocatable :: variables_1
  
      if(target_%name_.ne.'bem') &
         call out_%error("Target name: "//trim(target_%name_)// &
                         " not recognised")
       
      call mem_man%dealloc(variables_1, "variables_1")
      call mem_man%dealloc(rhs_1, "rhs_1")
      
   end subroutine deallocate_constant_potential_variables

!-------------------------------------------------------------------------------
!  SUBROUTINES FOR CONSTRUCTING MATRICES/VECTORS
!-------------------------------------------------------------------------------

   !> Subroutine for constructing the static matrix for FQ/FQFMu models
   !!    In/Out  : target_      -- model 
   subroutine construct_static_matrix(target_)
  
      implicit none
  
      !input/output variables
      class(target_type), intent(inout) :: target_

      if(target_%name_.ne.'fq'        .and. &
         target_%name_.ne.'fqfmu'     .and. &
         target_%name_.ne.'wfq'       .and. &
         target_%name_.ne.'wfqfmu'    .and. &
         target_%name_.ne.'bem') &
         call out_%error("Target name: "//trim(target_%name_)//&
                         " not recognised")

   end subroutine construct_static_matrix


   !> Subroutine for constructing constant matrix
   !!    In/Out  : target_             -- model 
   !!    In/Out  : matrix_constant     -- part of matrix_w that is constant
   subroutine construct_constant_matrix(target_,matrix_constant)
  
      implicit none
       
      !input/output variables
      class(target_type), intent(inout) :: target_
      real(dp), dimension(target_%n_var,target_%n_var), intent(inout) :: &
                                                                matrix_constant
  
      if(target_%name_.ne.'wfq'       .and. &
         target_%name_.ne.'wfqfmu'    .and. &
         target_%name_.ne.'bem') &
         call out_%error("Target name: "//trim(target_%name_)//&
                         " not recognised")
       
   end subroutine construct_constant_matrix


   !> Subroutine for constructing dynamic matrix
   !!    In/Out  : target_          -- model 
   !!    Input   : i_freq           -- index of frequency
   !!    Input   : matrix_constant  -- constant part of the matrix
   !!    In/Out  : matrix_w         -- dynamic matrix
   subroutine construct_dynamic_matrix(target_,i_freq,matrix_constant,matrix_w)
  
      implicit none
  
      !input/output
      class(target_type), intent(inout)  :: target_
      integer, intent(in)                :: i_freq
      real(dp), dimension(target_%n_var,target_%n_var), intent(in) :: &
                                                                  matrix_constant
      complex(dp), dimension(target_%n_var,target_%n_var), intent(inout) :: &
                                                                         matrix_w
  
      if(target_%name_.ne.'wfq'       .and. &
         target_%name_.ne.'wfqfmu'    .and. &
         target_%name_.ne.'bem') &
         call out_%error("Target name: "//trim(target_%name_)//&
                         " not recognised")
  
   end subroutine construct_dynamic_matrix


   !> Subroutine for constructing dynamic matrix [GMRES algorithm]
   !! This is only needed if the constant matrix depends on the frequency
   !! i.e. right now only for heterogeneous structures
   !!    In/Out  : target_          -- model 
   !!    Input   : i_freq           -- index of frequency
   !!    In/Out  : matrix_iterative -- dynamic matrix iterative
   subroutine construct_dynamic_matrix_gmres(target_,i_freq,matrix_iterative)
  
      implicit none
       
      !input/output variables
      class(target_type), intent(inout) :: target_
      integer, intent(in)               :: i_freq
      complex(dp), dimension(target_%n_var,target_%n_var), intent(inout) ::&
                                                                matrix_iterative
  
      if(target_%name_.ne.'wfqfmu'    .and. &
         target_%name_.ne.'wfqfmu_bem') &
         call out_%error("Target name: "//trim(target_%name_)//&
                         " not recognised")
       
   end subroutine construct_dynamic_matrix_gmres


   !> Subroutine for constructing the GS RHS 
   !!    In/Out  : target_      -- model 
   subroutine construct_ground_state_rhs(target_)

      implicit none

      !input/output variables
      class(target_type), intent(inout) :: target_

      if(target_%name_.ne.'fq'        .and. &
         target_%name_.ne.'fqfmu'     .and. &
         target_%name_.ne.'wfq'       .and. &
         target_%name_.ne.'wfqfmu'    .and. &
         target_%name_.ne.'bem') &
         call out_%error("Target name: "//trim(target_%name_)//&
                         " not recognised")
     
   end subroutine construct_ground_state_rhs


   !> Subroutine for constructing the RHS when static field is applied
   !!    In/Out  : target_      -- model 
   subroutine construct_static_field_rhs(target_)

      implicit none

      !input/output variables
      class(target_type), intent(inout) :: target_

      if(target_%name_.ne.'fq'        .and. &
         target_%name_.ne.'fqfmu'     .and. &
         target_%name_.ne.'wfq'       .and. &
         target_%name_.ne.'wfqfmu'    .and. &
         target_%name_.ne.'bem') &
         call out_%error("Target name: "//trim(target_%name_)// &
                         " not recognised")

   end subroutine construct_static_field_rhs


   !> Subroutine for constructing dynamic field RHS
   !!    In/Out  : target_   -- model 
   !!    In/Out  : rhs_w     -- dynamic RHS
   subroutine construct_dynamic_field_rhs(target_,rhs_w)
   
      implicit none
       
      !input/output variables
      class(target_type), intent(inout) :: target_
      complex(dp), dimension(target_%n_var,3), intent(inout) :: rhs_w
   
      if(target_%name_.ne.'wfq'       .and. &
         target_%name_.ne.'wfqfmu'    .and. &
         target_%name_.ne.'bem') &
         call out_%error("Target name: "//trim(target_%name_)//&
                         " not recognised")
       
   end subroutine construct_dynamic_field_rhs


   !> Subroutine for constructing dynamic field RHS [on the fly algorithm]
   !!    In/Out  : target_   -- model 
   !!    In/Out  : rhs_w     -- dynamic RHS
   subroutine construct_dynamic_field_rhs_on_the_fly(target_,rhs_w)
   
      implicit none
      !input/output variables
      class(target_type), intent(inout) :: target_
      complex(dp), dimension(target_%n_var,3), intent(inout) :: rhs_w
   
      if(target_%name_.ne.'wfq'       .and. &
         target_%name_.ne.'wfqfmu'    .and. &
         target_%name_.ne.'bem') &
         call out_%error("Target name: "//trim(target_%name_)//" not recognised")
        
   end subroutine construct_dynamic_field_rhs_on_the_fly


   !> Subroutine for performing the matrix vector multiplication on the fly
   !!    In/Out  : target_   -- model 
   !!    Input   : i_freq    -- index of frequency
   !!    In/Out  : x         -- input vector 
   !!    In/Out  : y         -- output vector 
   subroutine product_matrix_vector(target_, i_freq, x, y)
  
      implicit none
  
      !input/output variables
      class(target_type), intent(inout) :: target_
      integer, intent(in)               :: i_freq
      complex(dp), dimension(target_%n_var), intent(inout) :: x
      complex(dp), dimension(target_%n_var), intent(inout) :: y
  
      if(target_%name_.ne.'wfq'       .and. &
         target_%name_.ne.'wfqfmu'    .and. &
         target_%name_.ne.'bem') &
         call out_%error("Target name: "//trim(target_%name_)//&
                         " not recognised")
       
   end subroutine product_matrix_vector


   !> Subroutine for applying the diagonal shift for GMRES algorithm
   !!    In/Out  : target_      -- model 
   !!    Input   : i_freq       -- index of the frequency
   !!    Input   : colx         -- index of the column x (GMRES)
   !!    Input   : colz         -- index of the column z (GMRES)
   !!    Input   : lwork        -- dimension of the work array
   !!    In/Out  : work         -- array work
   subroutine gmres_diagonal_shift(target_,i_freq,colx,colz,lwork,work)

      implicit none

      !input/output variables
      class(target_type)    :: target_
      integer  :: i_freq
      integer  :: colx
      integer  :: colz
      integer  :: lwork
      complex(dp), dimension(lwork), intent(inout) :: work

      if(target_%name_.ne.'wfq'       .and. &
         target_%name_.ne.'wfqfmu'    .and. &
         target_%name_.ne.'bem') &
         call out_%error("Target name: "//trim(target_%name_)//&
                         " not recognised")

   end subroutine gmres_diagonal_shift


   !> Subroutine for applying the diagonal shift for GMRES algorithm
   !!    In/Out  : target_      -- model 
   !!    Input   : i_freq       -- index of the frequency
   !!    Input   : r            -- rhs
   !!    In/Out  : w            -- solution
   subroutine new_gmres_diagonal_shift(target_,i_freq,r,w)

      implicit none

      !input/output variables
      class(target_type)    :: target_
      integer  :: i_freq
      complex(dp), dimension(target_%n_var), intent(in) :: r
      complex(dp), dimension(target_%n_var), intent(inout) :: w

      if(target_%name_.ne.'wfq'       .and. &
         target_%name_.ne.'wfqfmu'    .and. &
         target_%name_.ne.'bem') &
         call out_%error("Target name: "//trim(target_%name_)//&
                         " not recognised")

   end subroutine new_gmres_diagonal_shift

!-------------------------------------------------------------------------------
!  SUBROUTINES FOR PRINTING & SAVING FILES
!-------------------------------------------------------------------------------

   !> Subroutine for printing the atomtypes
   !!    Input  : target_       -- model 
   !!    Input  : used_defaults -- if defaults parameters have been used
   subroutine print_info_target(target_, used_defaults)
  
      use field_module

      implicit none
       
      !input/output variables
      class(target_type), intent(in) :: target_
      logical, optional, intent(in) :: used_defaults
  
      write(out_%iunit,'(1x,a)') "Target             : "// &
                                 trim(target_%name_)
      if(target_%name_.ne.'bem') then
         write(out_%iunit,'(1x,a)') "ForceField         : "// &
                                    trim(target_%forcefield)
         write(out_%iunit,'(1x,a)') "Kernel             : "// &
                                    trim(target_%kernel)
         write(out_%iunit,'(1x,a)') "RHS Type           : "// &
                                     trim(field%rhs_form)
         !atomtypes informations
         write(out_%iunit,'(1x,a,i3)') "Num. AtomTypes     : ", &
                                       target_%n_atomtypes
         write(out_%iunit,'(1x,a,l1)') "Heterogeneous      : ", &
                                       target_%heterogeneous
         if (present(used_defaults)) &
            write(out_%iunit,'(1x,a,l1)') "Used defaults      : ", &
                                          used_defaults
         call target_%print_atomtypes()
         write(out_%iunit,out_%sticks)
         flush(out_%iunit)
      endif
       
   end subroutine print_info_target


   !> Subroutine for printing the atomtypes
   !!    Input  : target_   -- model 
   subroutine print_atomtypes(target_)
  
      implicit none
       
      !input/output variables
      class(target_type), intent(in) :: target_
  
      if(target_%name_.ne.'fq'        .and. &
         target_%name_.ne.'fqfmu'     .and. &
         target_%name_.ne.'wfq'       .and. &
         target_%name_.ne.'wfqfmu'    .and. &
         target_%name_.ne.'bem') &
         call out_%error("Target name: "//trim(target_%name_)//&
                         " not recognised")
       
   end subroutine print_atomtypes


   !> Subroutine for printing the coordinates
   !!    In/Out  : target_   -- model 
   !!    Input   : what      -- what to calculate
   subroutine print_coord(target_,what)
  
      implicit none
       
      !input/output variables
      class(target_type)    :: target_
      character(len=*)      :: what
  
      !internal variables
      integer :: i 
      character(len=76) :: format_1 = "(13x,'Atom',15x,'X',19x,'Y',19x,'Z')"
      character(len=50) :: format_2 = "(3x,i9,1x,a4,3(f20.6))"
  
      if(target_%name_.ne.'fq'        .and. &
         target_%name_.ne.'fqfmu'     .and. &
         target_%name_.ne.'wfq'       .and. &
         target_%name_.ne.'wfqfmu'    .and. &
         target_%name_.ne.'bem') &
         call out_%error("Target name: "//trim(target_%name_)//&
                         " not recognised")

      write(out_%iunit,'(/,31x,a,/)') trim(what)//" Geometry (Ang)"
      write(out_%iunit,out_%sticks)
      write(out_%iunit,format_1)
      write(out_%iunit,out_%sticks)

      if(target_%n_atoms.le.20000.or.out_%ivrb.eq.3) then
         do i = 1, target_%n_atoms
            write(out_%iunit,format_2) i,                        &
              trim(target_%atom_type(target_%map_atomtypes(i))), &
                                       target_%coord(1,i)*ToAng, &
                                       target_%coord(2,i)*ToAng, &
                                       target_%coord(3,i)*ToAng
         enddo
      else
         do i = 1, 2000
            write(out_%iunit,format_2) i,                        &
              trim(target_%atom_type(target_%map_atomtypes(i))), &
                                       target_%coord(1,i)*ToAng, &
                                       target_%coord(2,i)*ToAng, &
                                       target_%coord(3,i)*ToAng
  
         enddo
         write(out_%iunit,'(a,/,/)')
         do i = target_%n_atoms-2000, target_%n_atoms
            write(out_%iunit,format_2) i,                        &
              trim(target_%atom_type(target_%map_atomtypes(i))), &
                                       target_%coord(1,i)*ToAng, &
                                       target_%coord(2,i)*ToAng, &
                                       target_%coord(3,i)*ToAng
         enddo
      endif
  
      write(out_%iunit,out_%sticks)
       
   end subroutine print_coord


   !> Subroutine for printing GS variables
   !!    In/Out  : target_      -- model 
   subroutine print_gs_variables(target_)

      implicit none

      !input/output variables
      class(target_type), intent(in) :: target_

      if(target_%name_.ne.'fq'        .and. &
         target_%name_.ne.'fqfmu'     .and. &
         target_%name_.ne.'wfq'       .and. &
         target_%name_.ne.'wfqfmu'    .and. &
         target_%name_.ne.'bem') &
         call out_%error("Target name: "//trim(target_%name_)// &
                         " not recognised")

   end subroutine print_gs_variables


   !> Subroutine for printing the variables when static field is applied
   !!    In/Out  : target_      -- model 
   subroutine print_static_field_variables(target_)

      implicit none

      !input/output variables
      class(target_type), intent(in) :: target_

      if(target_%name_.ne.'fq'        .and. &
         target_%name_.ne.'fqfmu'     .and. &
         target_%name_.ne.'wfq'       .and. &
         target_%name_.ne.'wfqfmu'    .and. &
         target_%name_.ne.'bem') &
         call out_%error("Target name: "//trim(target_%name_)// &
                         " not recognised")

   end subroutine print_static_field_variables


   !> Subroutine for printing the variables when dynamic field is applied
   !!    Input   : target_      -- model 
   !!    Input   : variables_w  -- w-variables (N_var, 3 [x,y,z])
   subroutine print_dynamic_field_variables(target_,variables_w)

      implicit none

      !input/output variables
      class(target_type),intent(in) :: target_
      complex(dp), dimension(target_%n_var,3), intent(in) :: variables_w

      if(target_%name_.ne.'wfq'       .and. &
         target_%name_.ne.'wfqfmu'    .and. &
         target_%name_.ne.'bem') &
         call out_%error("Target name: "//trim(target_%name_)//&
                         " not recognised")

   end subroutine print_dynamic_field_variables


   !> Subroutine for printing dynamic polarizability
   !!    In/Out  : target_      -- model 
   !!    Input   : i_freq       -- index of the frequency
   !!    Input   : polar_w      -- dynamic polarizability (3,3)
   subroutine print_dynamic_polar(target_,i_freq,polar_w)

      use field_module
      implicit none
       
      !input/output variables
      class(target_type)    :: target_
      integer, intent(in)                     :: i_freq
      complex(dp), dimension(3,3), intent(in) :: polar_w
 
      !internal variables
      real(dp)          :: freq_eV
      real(dp)          :: freq_au
      real(dp)          :: freq_nm
 
      character(len=70) :: format_1 = "(1x,'Results for w = ', E9.3, ' a.u.',&
                                        &3x, E9.3, ' nm',3x,E9.3,' eV',/)"
      character(len=28) :: format_2 = "(28x,'Real',19x,'Imaginary')"
      character(len=55) :: format_iso = "(1x,'Pol Isotr   = ',6x,E14.6,' a.u.',&
                                          &6x,E14.6,' a.u.')"
      character(len=55) :: format_x = "(1x,'Pol Long X  = ',6x,E14.6,' a.u.',&
                                        &6x,E14.6,' a.u.')"
      character(len=55) :: format_y = "(1x,'Pol Long Y  = ',6x,E14.6,' a.u.',&
                                        &6x,E14.6,' a.u.')"
      character(len=55) :: format_z = "(1x,'Pol Long Z  = ',6x,E14.6,' a.u.',&
                                        &6x,E14.6,' a.u.')"
 
      if(target_%name_.ne.'wfq'       .and. &
         target_%name_.ne.'wfqfmu'    .and. &
         target_%name_.ne.'bem') &
         call out_%error("Target name: "//trim(target_%name_)// &
                         " not recognised")
 
      call freqautonm(field%freq(i_freq),freq_nm)
      call freqautoev(field%freq(i_freq),freq_eV)
      freq_au = field%freq(i_freq)
 
      if (out_%ivrb.ge.2) &
          write(out_%iunit,format_1) freq_au,freq_nm, freq_ev
       
      If (out_%ivrb.ge.2) then 
         call out_%print_matrix('Real Polarizability Tensor',DBle(polar_w),3,3)
         call out_%print_matrix('Imaginary Polarizability Tensor', &
                                aimag(polar_w),3,3)
      endIf
       
      if(out_%ivrb.ge.2) then
         write(out_%iunit,format_2) 
         write(out_%iunit,format_iso) target_%results(i_freq,1), &
                                      target_%results(i_freq,2)
         write(out_%iunit,format_x)   target_%results(i_freq,3), &
                                      target_%results(i_freq,6)
         write(out_%iunit,format_y)   target_%results(i_freq,4), &
                                      target_%results(i_freq,7)
         write(out_%iunit,format_z)   target_%results(i_freq,5), &
                                      target_%results(i_freq,8)
         write(out_%iunit,out_%sticks) 
      endif
      
   end subroutine print_dynamic_polar


   !> Subroutine for printing dynamic results
   !!    In/Out  : target_      -- model 
   subroutine print_dynamic_results(target_)
  
      use field_module
  
      implicit none
       
      !input/output variables 
      class(target_type)    :: target_
   
      !internal variables
      integer           :: i
      real(dp)          :: freq_eV
      real(dp)          :: freq_au
      real(dp)          :: freq_nm
      character(len=70) :: format_1 = "(1x,'Results for w = ', E9.3, ' a.u.',&
                                      &3x, E9.3, ' nm',3x,E9.3,' eV',/)"
      character(len=23) :: format_2 = "(1x,a,6x,E14.6,' a.u.')"
  
      if(target_%name_.ne.'wfq'       .and. &
         target_%name_.ne.'wfqfmu'    .and. &
         target_%name_.ne.'bem') &
         call out_%error("Target name: "//trim(target_%name_)//&
                         " not recognised")
  
      do i = 1, field%n_freq
 
         call freqautonm(field%freq(i),freq_nm)
         call freqautoev(field%freq(i),freq_eV)
         freq_au = field%freq(i)
         write(out_%iunit,format_1) freq_au, freq_nm, freq_ev
         write(out_%iunit,format_2) 'Isotr. Real Polar.       = ', &
                                     target_%results(i,1)
         write(out_%iunit,format_2) 'Isotr. Imag Polar.       = ', &
                                     target_%results(i,2)
         write(out_%iunit,format_2) 'Long. Real Polar. X      = ', &
                                     target_%results(i,3)
         write(out_%iunit,format_2) 'Long. Real Polar. Y      = ', &
                                     target_%results(i,4)
         write(out_%iunit,format_2) 'Long. Real Polar. Z      = ', &
                                     target_%results(i,5)
         write(out_%iunit,format_2) 'Long. Imag Polar. X      = ', &
                                     target_%results(i,6)
         write(out_%iunit,format_2) 'Long. Imag Polar. Y      = ', &
                                     target_%results(i,7)
         write(out_%iunit,format_2) 'Long. Imag Polar. Z      = ', &
                                     target_%results(i,8)
         write(out_%iunit,format_2) 'Iso. Abs. Cross. Sect.   = ', &
                                     target_%results(i,9)
         write(out_%iunit,format_2) 'Long Abs. Cross. Sect. X = ', &
                                     target_%results(i,10)
         write(out_%iunit,format_2) 'Long Abs. Cross. Sect. Y = ', &
                                     target_%results(i,11)
         write(out_%iunit,format_2) 'Long Abs. Cross. Sect. Z = ', &
                                     target_%results(i,12)
         write(out_%iunit,format_2) 'Iso. Sca. Cross. Sect.   = ', &
                                     target_%results(i,13)
         write(out_%iunit,format_2) 'Long Sca. Cross. Sect. X = ', &
                                     target_%results(i,14)
         write(out_%iunit,format_2) 'Long Sca. Cross. Sect. Y = ', &
                                     target_%results(i,15)
         write(out_%iunit,format_2) 'Long Sca. Cross. Sect. Z = ', &
                                     target_%results(i,16)
         write(out_%iunit,format_2) 'Iso. Ext. Cross. Sect.   = ', &
                                     target_%results(i,17)
         write(out_%iunit,format_2) 'Long Ext. Cross. Sect. X = ', &
                                     target_%results(i,18)
         write(out_%iunit,format_2) 'Long Ext. Cross. Sect. Y = ', &
                                     target_%results(i,19)
         write(out_%iunit,format_2) 'Long Ext. Cross. Sect. Z = ', &
                                     target_%results(i,20)
         write(out_%iunit,out_%sticks) 
      enddo
       
   end subroutine print_dynamic_results


   !> Subroutine for saving the intermediate results for restarting
   !! by creating the Results array
   !!    Results(FreqI,1)       Isotr. Real Polar.         
   !!    Results(FreqI,2)       Isotr. Imag Polar.         
   !!    Results(FreqI,3)       Long. Real Polar. X        
   !!    Results(FreqI,4)       Long. Real Polar. Y        
   !!    Results(FreqI,5)       Long. Real Polar. Z        
   !!    Results(FreqI,6)       Long. Imag Polar. X        
   !!    Results(FreqI,7)       Long. Imag Polar. Y        
   !!    Results(FreqI,8)       Long. Imag Polar. Z        
   !!    Results(FreqI,9)       Iso. Abs. Cross. Sect.     
   !!    Results(FreqI,10)      Long Abs. Cross. Sect. X   
   !!    Results(FreqI,11)      Long Abs. Cross. Sect. Y   
   !!    Results(FreqI,12)      Long Abs. Cross. Sect. Z   
   !!    Results(FreqI,13)      Iso. Sca. Cross. Sect.     
   !!    Results(FreqI,14)      Long Sca. Cross. Sect. X   
   !!    Results(FreqI,15)      Long Sca. Cross. Sect. Y   
   !!    Results(FreqI,16)      Long Sca. Cross. Sect. Z   
   !!    Results(FreqI,17)      Iso. Ext. Cross. Sect.     
   !!    Results(FreqI,18)      Long Ext. Cross. Sect. X   
   !!    Results(FreqI,19)      Long Ext. Cross. Sect. Y   
   !!    Results(FreqI,20)      Long Ext. Cross. Sect. Z   
   !!    Results(FreqI,21)      Current through Origin X   
   !!    Results(FreqI,22)      Current through Origin Y   
   !!    Results(FreqI,23)      Current through Origin Z   
   !!
   !!    In/Out  : target_   -- model 
   !!    Input   : i_freq    -- index of frequency
   subroutine save_intermediate_results(target_,i_freq)
  
      use field_module
  
      implicit none
       
      !input/output variables
      class(target_type), intent(inout)  :: target_
      integer, intent(in)                :: i_freq
  
      !internal variables
      integer               :: iunit_bk
      integer               :: iost
      real(dp)              :: freq_ev
      real(dp)              :: freq_au
      real(dp)              :: freq_nm
      character(len=99)     :: output_bk
      character(len=7)      :: freq_char
      character(len=70)     :: format_1 = "(1x,'Results for w = ', E9.3, &
                                 &' a.u.',3x, E9.3, ' nm',3x,E9.3,' eV',/)"
      character(len=23)     :: format_2 = "(1x,a,6x,E14.6,' a.u.')"
  
      if(target_%name_.ne.'wfq'       .and. &
         target_%name_.ne.'wfqfmu'    .and. &
         target_%name_.ne.'bem') &
         call out_%error("Target name: "//trim(target_%name_)//&
                         " not recognised")
       
      call freqautonm(field%freq(i_freq),freq_nm)
      call freqautoev(field%freq(i_freq),freq_eV)
      freq_au = field%freq(i_freq)
      write(freq_char,'(f7.5)') freq_ev
      write(output_bk,'(a)') out_%filename(1:(len_trim(out_%filename)-4)) // &
                                                '-' // freq_char // '.plasmonX.bk'

      iunit_bk = 15
      open(unit=iunit_bk,file=output_bk,status="unknown",iostat=iost)
         write(iunit_bk,out_%sticks) 
         write(iunit_bk,format_1) freq_au, freq_nm, freq_ev
         write(iunit_bk,format_2) 'Isotr. Real Polar.       = ', &
                                  target_%results(i_freq,1)
         write(iunit_bk,format_2) 'Isotr. Imag Polar.       = ', &
                                  target_%results(i_freq,2)
         write(iunit_bk,format_2) 'Long. Real Polar. X      = ', &
                                  target_%results(i_freq,3)
         write(iunit_bk,format_2) 'Long. Real Polar. Y      = ', &
                                  target_%results(i_freq,4)
         write(iunit_bk,format_2) 'Long. Real Polar. Z      = ', &
                                  target_%results(i_freq,5)
         write(iunit_bk,format_2) 'Long. Imag Polar. X      = ', &
                                  target_%results(i_freq,6)
         write(iunit_bk,format_2) 'Long. Imag Polar. Y      = ', &
                                  target_%results(i_freq,7)
         write(iunit_bk,format_2) 'Long. Imag Polar. Z      = ', &
                                  target_%results(i_freq,8)
         write(iunit_bk,format_2) 'Iso. Abs. Cross. Sect.   = ', &
                                  target_%results(i_freq,9)
         write(iunit_bk,format_2) 'Long Abs. Cross. Sect. X = ', &
                                  target_%results(i_freq,10)
         write(iunit_bk,format_2) 'Long Abs. Cross. Sect. Y = ', &
                                  target_%results(i_freq,11)
         write(iunit_bk,format_2) 'Long Abs. Cross. Sect. Z = ', &
                                  target_%results(i_freq,12)
         write(iunit_bk,format_2) 'Iso. Sca. Cross. Sect.   = ', &
                                  target_%results(i_freq,13)
         write(iunit_bk,format_2) 'Long Sca. Cross. Sect. X = ', &
                                  target_%results(i_freq,14)
         write(iunit_bk,format_2) 'Long Sca. Cross. Sect. Y = ', &
                                  target_%results(i_freq,15)
         write(iunit_bk,format_2) 'Long Sca. Cross. Sect. Z = ', &
                                  target_%results(i_freq,16)
         write(iunit_bk,format_2) 'Iso. Ext. Cross. Sect.   = ', &
                                  target_%results(i_freq,17)
         write(iunit_bk,format_2) 'Long Ext. Cross. Sect. X = ', &
                                  target_%results(i_freq,18)
         write(iunit_bk,format_2) 'Long Ext. Cross. Sect. Y = ', &
                                  target_%results(i_freq,19)
         write(iunit_bk,format_2) 'Long Ext. Cross. Sect. Z = ', &
                                  target_%results(i_freq,20)
         write(iunit_bk,out_%sticks) 
         flush(iunit_bk)
      close(iunit_bk)
  
   end subroutine save_intermediate_results


   !> Subroutine for saving the csv file containing the results vector
   !!    In/Out  : target_   -- model 
   subroutine save_csv_file(target_)
  
      use field_module
  
      implicit none
       
      !input/output variables
      class(target_type)    :: target_
  
      !internal variables
      integer               :: iunit_csv
      integer               :: iost
      integer               :: i,j
      real(dp)              :: freq_ev
      character(len=25)     :: format_1 ="(1x,E12.5,20(3X,E14.6))"
      character(len=99)     :: output_csv
  
      if(target_%name_.ne.'wfq'       .and. &
         target_%name_.ne.'wfqfmu'    .and. &
         target_%name_.ne.'bem') &
         call out_%error("Target name: "//trim(target_%name_)//&
                         " not recognised")
       
      iunit_csv = 18
      Write(output_csv,'(a)') out_%filename(1:(len_trim(out_%filename)-4)) // &
                                                                          '.csv'
      open(unit=iunit_csv,file=output_csv,status="unknown",iostat=iost)
         write(iunit_csv,'(a)') '#------------------------------------------'
         write(iunit_csv,'(a)') '#    Column          Quantity'
         write(iunit_csv,'(a)') '#------------------------------------------'
         write(iunit_csv,'(a)') '#      1        Frequency (eV)             '
         write(iunit_csv,'(a)') '#      2        Isotr. Real Polar.         '
         write(iunit_csv,'(a)') '#      3        Isotr. Imag Polar.         '
         write(iunit_csv,'(a)') '#      4        Long. Real Polar. X        '
         write(iunit_csv,'(a)') '#      5        Long. Real Polar. Y        '
         write(iunit_csv,'(a)') '#      6        Long. Real Polar. Z        '
         write(iunit_csv,'(a)') '#      7        Long. Imag Polar. X        '
         write(iunit_csv,'(a)') '#      8        Long. Imag Polar. Y        '
         write(iunit_csv,'(a)') '#      9        Long. Imag Polar. Z        '
         write(iunit_csv,'(a)') '#      10       Iso. Abs. Cross. Sect.     '
         write(iunit_csv,'(a)') '#      11       Long Abs. Cross. Sect. X   '
         write(iunit_csv,'(a)') '#      12       Long Abs. Cross. Sect. Y   '
         write(iunit_csv,'(a)') '#      13       Long Abs. Cross. Sect. Z   '
         write(iunit_csv,'(a)') '#      14       Iso. Sca. Cross. Sect.     '
         write(iunit_csv,'(a)') '#      15       Long Sca. Cross. Sect. X   '
         write(iunit_csv,'(a)') '#      16       Long Sca. Cross. Sect. Y   '
         write(iunit_csv,'(a)') '#      17       Long Sca. Cross. Sect. Z   '
         write(iunit_csv,'(a)') '#      18       Iso. Ext. Cross. Sect.     '
         write(iunit_csv,'(a)') '#      19       Long Ext. Cross. Sect. X   '
         write(iunit_csv,'(a)') '#      20       Long Ext. Cross. Sect. Y   '
         write(iunit_csv,'(a)') '#      21       Long Ext. Cross. Sect. Z   '
         Write(iunit_csv,'(a)') '#------------------------------------------'
         do i = 1, field%n_freq        
            call freqautoev(field%freq(i), freq_ev)
            write(iunit_csv,format_1) freq_ev, (target_%results(i,j), j=1,20)
         enddo
      close(iunit_csv)
  
   end subroutine save_csv_file

!-------------------------------------------------------------------------------
!  SUBROUTINES FOR CALCULATING PROPERTIES
!-------------------------------------------------------------------------------

   !> Subroutine for calculating the energy of the GS 
   !!    In/Out  : target_      -- model 
   subroutine calculate_energy(target_)

      implicit none

      !input/output variables
      class(target_type), intent(inout) :: target_

      if(target_%name_.ne.'fq'        .and. &
         target_%name_.ne.'fqfmu'     .and. &
         target_%name_.ne.'wfq'       .and. &
         target_%name_.ne.'wfqfmu'    .and. &
         target_%name_.ne.'bem') &
         call out_%error("Target name: "//trim(target_%name_)//&
                         " not recognised")

   end subroutine calculate_energy


   !> Subroutine for calculating the static polar
   !!    In/Out  : target_      -- model 
   subroutine calculate_static_polar(target_)

      implicit none

      !input/output variables
      class(target_type), intent(inout) :: target_

      if(target_%name_.ne.'fq'        .and. &
         target_%name_.ne.'fqfmu'     .and. &
         target_%name_.ne.'wfq'       .and. &
         target_%name_.ne.'wfqfmu'    .and. &
         target_%name_.ne.'bem') &
         call out_%error("Target name: "//trim(target_%name_)//&
                         " not recognised")

   end subroutine calculate_static_polar


   !> Subroutine for calculating the dynamic polar
   !!    In/Out  : target_      -- model 
   !!    Input   : i_freq       -- index of the frequency
   !!    Input   : variables_w  -- w-variables
   !!    Output  : polar_w      -- dynamic polar
   subroutine calculate_dynamic_polar(target_,i_freq,variables_w,polar_w)
   
      use field_module
  
      implicit none
       
      !input/output variables
      class(target_type), intent(inout) :: target_
      integer, intent(in)               :: i_freq
      complex(dp), dimension(target_%n_var,3), intent(in)    :: variables_w
      complex(dp), dimension(3,3), intent(inout)             :: polar_w
  
      if(target_%name_.ne.'wfq'       .and. &
         target_%name_.ne.'wfqfmu'    .and. &
         target_%name_.ne.'bem') &
         call out_%error("Target name: "//trim(target_%name_)//&
                         " not recognised")
  
      if(i_freq.eq.0) call out_%error('Frequency index is equal to 0 in &
                                      &calculate_dynamic polar')
       
   end subroutine calculate_dynamic_polar


   !> Subroutine for calculating the cross sections
   !!    In/Out  : target_      -- model 
   !!    Input   : i_freq       -- index of frequency
   subroutine calculate_cross_section(target_,i_freq)
 
      use field_module
  
      implicit none
       
      !input/output variables
      class(target_type)    :: target_
      integer               :: i_freq
       
      !internal variables
      real(dp)              :: freq_au
      real(dp)              :: polar_module_iso
      real(dp)              :: polar_module_x
      real(dp)              :: polar_module_y
      real(dp)              :: polar_module_z
  
      if(target_%name_.ne.'wfq'       .and. &
         target_%name_.ne.'wfqfmu'    .and. &
         target_%name_.ne.'bem') &
         call out_%error("Target name: "//trim(target_%name_)//&
                         " not recognised")
          
      freq_au = field%freq(i_freq)
  
      polar_module_iso = target_%results(i_freq,1)**2 + &
                         target_%results(i_freq,2)**2
      polar_module_x   = target_%results(i_freq,3)**2 + &
                         target_%results(i_freq,6)**2
      polar_module_y   = target_%results(i_freq,4)**2 + &
                         target_%results(i_freq,7)**2
      polar_module_z   = target_%results(i_freq,5)**2 + &
                         target_%results(i_freq,8)**2
  
      !absorption cross section : 4piw/c * imag_polar 
      target_%results(i_freq,9)  = four*pi*freq_au*target_%results(i_freq,2)/& 
                                   light ! iso 
      target_%results(i_freq,10) = four*pi*freq_au*target_%results(i_freq,6)/&
                                   light ! X
      target_%results(i_freq,11) = four*pi*freq_au*target_%results(i_freq,7)/&
                                   light ! Y
      target_%results(i_freq,12) = four*pi*freq_au*target_%results(i_freq,8)/&
                                   light ! Z
  
      !scattering cross section 
      target_%results(i_freq,13) = (eight*pi*freq_au**4)*polar_module_iso/&
                                   (three*light**4)! iso 
      target_%results(i_freq,14) = (eight*pi*freq_au**4)*polar_module_x/&
                                   (three*light**4)  ! X
      target_%results(i_freq,15) = (eight*pi*freq_au**4)*polar_module_y/&
                                   (three*light**4)  ! Y
      target_%results(i_freq,16) = (eight*pi*freq_au**4)*polar_module_z/&
                                   (three*light**4)  ! Z
  
      !extintion cross section : abs + sca
      target_%results(i_freq,17) = target_%results(i_freq,9)  + &
                                   target_%results(i_freq,13)! iso 
      target_%results(i_freq,18) = target_%results(i_freq,10) + &
                                   target_%results(i_freq,14)! X
      target_%results(i_freq,19) = target_%results(i_freq,11) + &
                                   target_%results(i_freq,15)! Y
      target_%results(i_freq,20) = target_%results(i_freq,12) + &
                                   target_%results(i_freq,16)! Z
  
       
   end subroutine calculate_cross_section


   !> Function for calculating the induced field at a specific point
   !!    Input   : target_      -- model 
   !!    Input   : i_pol        -- index of the polarization
   !!    Input   : point_coord  -- coordinates of the point
   !!    Input   : variables_w  -- w-variables
   !!    Output  : EField       -- Electric Field
   function calculate_induced_field_at_point(target_,i_pol,point_coord, &
                                             variables_w) result(EField)
   
      implicit none
       
      !input/output variables
      class(target_type), intent(in)     :: target_
      integer, intent(in) :: i_pol
      real(dp), dimension(3), intent(in) :: point_coord
      complex(dp), dimension(target_%n_var), intent(in) :: variables_w
      complex(dp), dimension(3) :: EField
  
      if(target_%name_.ne.'wfq'.and. &
         target_%name_.ne.'wfqfmu') &
         call out_%error("Target name: "//trim(target_%name_)//&
                         " not recognised")
  
      EField = zero
  
   end function calculate_induced_field_at_point


   !> Function for calculating the plasmon density at a specific point
   !!    Input   : target_      -- model 
   !!    Input   : point_coord  -- coordinates of the point
   !!    Input   : variables_w  -- w-variables
   !!    Output  : density      -- Plasmon density
   function calculate_density_at_point(target_,point_coord,variables_w) &
            result(density)
   
      implicit none
       
      !input/output variables
      class(target_type), intent(in)     :: target_
      real(dp), dimension(3), intent(in) :: point_coord
      complex(dp), dimension(target_%n_var), intent(in) :: variables_w
      complex(dp) :: density
   
      if(target_%name_.ne.'wfq'.and. &
         target_%name_.ne.'wfqfmu') &
         call out_%error("Target name: "//trim(target_%name_)//&
                         " not recognised")
   
      density = zero
   
   end function calculate_density_at_point


   !> Function for calculating the plasmon density at a specific point 
   !! separating q and mu contributions
   !!    Input   : target_      -- model 
   !!    Input   : point_coord  -- coordinates of the point
   !!    Input   : variables_w  -- w-variables
   !!    Output  : densities    -- 3D array: (1) FQ (2) FMu (3) Total
   function calculate_density_at_point_q_mu_separated(target_,     &
                                                      point_coord, &
                                                      variables_w) &
            result(densities)
   
      implicit none
       
      !input/output variables
      class(target_type), intent(in)     :: target_
      real(dp), dimension(3), intent(in) :: point_coord
      complex(dp), dimension(target_%n_var), intent(in) :: variables_w
      complex(dp), dimension(3) :: densities
   
      if(target_%name_.ne.'wfqmu') &
         call out_%error("Target name: "//trim(target_%name_)//&
                         " not recognised")
   
      densities(:) = zero
   
   end function calculate_density_at_point_q_mu_separated

end module target_module
