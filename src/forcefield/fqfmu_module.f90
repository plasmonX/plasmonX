!> FQFMu module
!!
!! This module contains the subroutines for the fq type
!!
!! Date         : 2024 
!!
module fqfmu_module
      

   use output_module
   use parameters_module
   use string_manipulation_module
   use array_manipulation_module
   use field_module
   use target_module
   use fq_module

   Implicit None

   !public variables
   public fqfmu

   type, extends(fq_type) :: fqfmu_type
      integer   :: n_fq = 0 ! number of charges + molecules (if present)
 
      contains
       
      !assign/get variables
      procedure :: assign_model_dimensions        => &
                   assign_model_dimensions_fqfmu
      !construct matrices/vectors
      procedure :: construct_static_matrix        => &
                   construct_static_matrix_fqfmu
      procedure :: construct_ground_state_rhs     => &
                   construct_ground_state_rhs_fqfmu
      procedure :: construct_static_field_rhs     => &
                   construct_static_field_rhs_fqfmu
      !printing & saving files
      procedure :: print_atomtypes                => &
                   print_atomtypes_fqfmu
      procedure :: print_gs_variables             => &
                   print_gs_variables_fqfmu
      procedure :: print_static_field_variables   => &
                   print_static_field_variables_fqfmu
      !calculate properties
      procedure :: calculate_energy               => &
                   calculate_energy_fqfmu
      procedure :: calculate_static_polar         => &
                   calculate_static_polar_fqfmu
    end type fqfmu_type
     
   type (fqfmu_type), target, save :: fqfmu

contains

!-------------------------------------------------------------------------------
!  SUBROUTINES FOR ASSIGNING & GETTING VARIABLES
!-------------------------------------------------------------------------------

   !> Subroutine for assigning the dimensions of the FQFMu model
   !!    In/Out : target_      -- FQFMu
   subroutine assign_model_dimensions_fqfmu(target_)
  
      implicit none

      !input/output variables
      class(fqfmu_type), intent(inout) :: target_
       
      target_%n_var = 4*target_%n_atoms + target_%n_mol
      target_%n_q   = target_%n_atoms
      target_%n_mu  = 3*target_%n_atoms
  
   end subroutine assign_model_dimensions_fqfmu

!-------------------------------------------------------------------------------
!  SUBROUTINES FOR CONSTRUCTING MATRICES/VECTORS
!-------------------------------------------------------------------------------

   !> Subroutine for constructing the static matrix for FQFMu models
   !!    In/Out  : target_      -- FQ
   subroutine construct_static_matrix_fqfmu(target_)
  
      implicit none
       
      !input/output variables
      class(fqfmu_type), intent(inout) :: target_
  
      !internal variables
      real(dp), dimension(:,:), allocatable :: tmp_qq
      real(dp), dimension(:,:), allocatable :: tmp_qmu
      real(dp), dimension(:,:), allocatable :: tmp_mumu
      real(dp), dimension(:,:), allocatable :: tmp_lang

      !Tqq block
      call mem_man%alloc(tmp_qq, target_%n_q, target_%n_q, "tmp_qq")
      call construct_static_t_qq(target_, tmp_qq)
      target_%matrix(1:target_%n_q,1:target_%n_q) = tmp_qq
      call mem_man%dealloc(tmp_qq, "tmp_qq")

      !Lagrangian block
      if(target_%n_var.ne.target_%n_q.and.target_%name_.ne.'wfqfmu') then
         target_%n_fq = target_%n_q + target_%n_mol
         call mem_man%alloc(tmp_lang, target_%n_fq, target_%n_fq, "tmp_lang")
         call add_lagrangian_blocks(target_, tmp_lang)
         target_%matrix(target_%n_q+1:target_%n_fq,1:target_%n_q) = &
               tmp_lang(target_%n_q+1:target_%n_fq,1:target_%n_q)
         target_%matrix(1:target_%n_q,target_%n_q+1:target_%n_fq) = &
               tmp_lang(1:target_%n_q,target_%n_q+1:target_%n_fq)
         call mem_man%dealloc(tmp_lang, "tmp_lang")
      else
         target_%n_fq = target_%n_q
         if(out_%ivrb.ge.2) &
            write(out_%iunit,'(1x,a)') ' Charge Constraints not setted'
      endif
      !Tqmu block
      call mem_man%alloc(tmp_qmu, target_%n_q, target_%n_mu, "tmp_qmu")
      call construct_static_t_qmu(target_, tmp_qmu)
      target_%matrix(1:target_%n_q,target_%n_fq+1:) = tmp_qmu
      !Tmuq block = Tqmu ** T
      target_%matrix(target_%n_fq+1:,1:target_%n_q) = transpose(tmp_qmu)
      call mem_man%dealloc(tmp_qmu, "tmp_qmu")
      !Tmumu block
      call mem_man%alloc(tmp_mumu, target_%n_mu, target_%n_mu, "tmp_mumu")
      call construct_static_t_mumu(target_, tmp_mumu)
      target_%matrix(target_%n_fq+1:,target_%n_fq+1:) = tmp_mumu
      call mem_man%dealloc(tmp_mumu, "tmp_mumu")
      !printing
      if(out_%ivrb.ge.3) &
         call out_%print_matrix('FQFMu Matrix',target_%matrix, &
                                 target_%n_var,target_%n_var)
       
   end subroutine construct_static_matrix_fqfmu


   !> Subroutine for constructing the GS RHS 
   !!    In/Out  : target_      -- FQFMu
   subroutine construct_ground_state_rhs_fqfmu(target_)
  
      implicit none
       
      !input/output variables
      class(fqfmu_type), intent(inout) :: target_
  
      !internal variables
      integer :: i
      character(len=14) :: format_1 = "(14x,'RHS ',/)"
      character(len=15) :: format_2 = "(i4,(4x,e14.6))"
  
      !from n_atoms to n_var, rhs = 0 [Lagrangian block]
      !this should be changed if the MM molecules has charge != 0
      !$omp parallel do private(i)
      do i=1,target_%n_atoms
         target_%rhs(i,1) = -target_%chi(target_%map_atomtypes(i))
      enddo
      !$omp end parallel do
      
      !printing
      if(out_%ivrb.ge.3) then
         write(out_%iunit,format_1)
         do i = 1,target_%n_var
            write(out_%iunit,format_2) i,target_%rhs(i,1)
         enddo
         write(out_%iunit,out_%sticks)
         flush(out_%iunit)
      endif
       
   end subroutine construct_ground_state_rhs_fqfmu


   !> Subroutine for constructing the RHS when static field is applied
   !!    In/Out  : target_      -- FQFMu
   subroutine construct_static_field_rhs_fqfmu(target_)
  
      implicit none
       
      !input/output variables
      class(fqfmu_type), intent(inout) :: target_
  
      !Q block
      call construct_static_field_rhs_q(target_, &
                                        target_%rhs(1:target_%n_atoms, :))
      !Mu block
      call construct_static_field_rhs_mu(target_, &
                                         target_%rhs(target_%n_fq+1:, :))
      !printing
      if(out_%ivrb.ge.3) call print_static_field_rhs(target_,target_%rhs)
       
   end subroutine construct_static_field_rhs_fqfmu

!-------------------------------------------------------------------------------
!  SUBROUTINES FOR PRINTING & SAVING FILES
!-------------------------------------------------------------------------------

   !> Subroutine for printing the atomtypes
   !!    Input  : target_   -- FQFMu
   subroutine print_atomtypes_fqfmu(target_)
       
      implicit none
  
      !input/output variables
      class(fqfmu_type), intent(in) :: target_
  
      !internal variables
      integer :: i
       
      do i = 1, target_%n_atomtypes
         write(out_%iunit,out_%sticks) 
         write(out_%iunit,'(1x,a)') "AtomType           : "//&
                                    trim(target_%atom_type(i))
         write(out_%iunit,'(1x,a,f8.3,a)') "Chi                : ", &
                                           target_%chi(i), " a.u."
         write(out_%iunit,'(1x,a,f8.3,a)') "Eta                : ", &
                                           target_%eta(i), " a.u."
         write(out_%iunit,'(1x,a,f8.3,a)') "Alpha              : ", &
                                           target_%alpha(i), &
                                            " a.u."
         write(out_%iunit,'(1x,a,f8.3,a)') "R_q                : ", &
                                           target_%r_q(i), " a.u."
         write(out_%iunit,'(1x,a,f8.3,a)') "R_mu               : ", &
                                           target_%r_mu(i)," a.u."
      enddo
       
   end subroutine print_atomtypes_fqfmu


   !> Subroutine for printing GS variables
   !!    In/Out  : target_      -- FQFMu
   subroutine print_gs_variables_fqfmu(target_)
  
      implicit none
       
      !input/output variables
      class(fqfmu_type), intent(in) :: target_
  
      !internal variables
      integer :: i
      real(dp) :: sum_charges = zero
      character(len=18) :: format_1 = "(14x,'Charges ',/)"
      character(len=18) :: format_2 = "(1x,I4,(4x,E14.6))"
      character(len=28) :: format_3 = "(/,1x,'ChErr',3X,(E14.6,4x))"
      character(len=18) :: format_4 = "(14x,'Dipoles ',/)"
      character(len=18) :: format_5 = "(1x,i4,(4x,e14.6))"
  
       !charges
       write(out_%iunit,format_1) 
       do i = 1, target_%n_atoms
          write(out_%iunit,format_2) i,target_%variables(i,1)
          sum_charges = sum_charges + target_%variables(i,1)
       enddo
       write(out_%iunit,format_3) sum_charges
       write(out_%iunit,out_%sticks) 
  
       !dipoles
       write(out_%iunit,format_4)
       do i = target_%n_atoms+target_%n_mol+1,target_%n_var
          write(out_%iunit,format_5) i-target_%n_atoms-target_%n_mol, &
                                     target_%variables(i,1)
       enddo
       write(out_%iunit,out_%sticks) 
       
   end subroutine print_gs_variables_fqfmu


   !> Subroutine for printing the variables when static field is applied
   !!    In/Out  : target_      -- model 
   subroutine print_static_field_variables_fqfmu(target_)
  
      implicit none
       
      !input/output variables
      class(fqfmu_type), intent(in) :: target_
  
      !internal variables
      integer :: i
      real(dp) :: sum_X = zero
      real(dp) :: sum_Y = zero
      real(dp) :: sum_Z = zero
      character(len=18) :: format_1 = "(31x,'Charges ',/)"
      character(len=55) :: format_2 = "(11x,'X Component',7x,'Y Component',7x,&
                                      &'Z Component',/)"
      character(len=19) :: format_3 = "(1x,I4,3(4x,E14.6))"
      character(len=31) :: format_4 = "(/,1x,'ChErr',3X,3(D14.6,4x),/)"
      character(len=18) :: format_5 = "(31x,'Dipoles ',/)"
  
       !charges
       write(out_%iunit,format_1)
       write(out_%iunit,format_2)
       do i = 1, target_%n_atoms
          write(out_%iunit,format_3) i,target_%variables(i,:)
          sum_x = sum_x + target_%variables(i,1)
          sum_y = sum_y + target_%variables(i,2)
          sum_z = sum_z + target_%variables(i,3)
       enddo
       write(out_%iunit,format_4) sum_x, sum_y, sum_z
       write(out_%iunit,out_%sticks) 
  
       !dipoles
       write(out_%iunit,format_5)
       write(out_%iunit,format_2)
       do i = target_%n_atoms+target_%n_mol+1,target_%n_var
          write(out_%iunit,format_3) i-target_%n_atoms-target_%n_mol, &
                                     target_%variables(i,:)
       enddo
       write(out_%iunit,out_%sticks) 
       
   end subroutine print_static_field_variables_fqfmu

!-------------------------------------------------------------------------------
!  SUBROUTINES FOR CALCULATING PROPERTIES
!-------------------------------------------------------------------------------

   !> Subroutine for calculating the energy of the GS 
   !!    In/Out  : target_      -- FQFMu
   subroutine calculate_energy_fqfmu(target_)
  
      implicit none
       
      !input/output variables
      class(fqfmu_type), intent(inout) :: target_
  
      !internal variables
      real(dp)          :: ddot
      character(len=14) :: format_1 = "(1x,a,f14.8,a)"
  
      !Energy = -1/2 * LHS * RHS 
      target_%energy = -half * &
                        ddot(target_%n_var,target_%variables,1,target_%rhs,1)
      write(out_%iunit,format_1) 'Energy = ', target_%energy,' a.u.'
      write(out_%iunit,out_%sticks)
       
   end subroutine calculate_energy_fqfmu


   !> Subroutine for calculating the static polar
   !!    In/Out  : target_      -- FQFMu
   subroutine calculate_static_polar_fqfmu(target_)
  
      implicit none
       
      !input/output variables
      class(fqfmu_type), intent(inout) :: target_
  
      if (target_%name_.ne.'fqfmu') &
          call out_%error("Static polar fqfmu but not FQFMu forcefield")
      call out_%error("Static Polar for FQFMu not yet implemented")
  
   end subroutine calculate_static_polar_fqfmu

end module fqfmu_module
