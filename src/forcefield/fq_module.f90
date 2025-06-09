!> FQ module
!!
!! This module contains the subroutines for the fq type
!!
!! Date         : 2025
!!
module fq_module
      
   use output_module
   use parameters_module
   use string_manipulation_module
   use array_manipulation_module
   use field_module
   use target_module
   use matrix_module
   use rhs_module

   implicit none

   !public variables
   public fq

   type, extends(target_type) :: fq_type

      contains
      
      !assign/get variables
      procedure :: assign_model_dimensions        => &
                   assign_model_dimensions_fq
      !allocations
      procedure :: allocate_gs_matrices           => &
                   allocate_gs_matrices_fq
      procedure :: allocate_static_field_matrices => &
                   allocate_static_field_matrices_fq
      procedure :: dealloc                        => &
                   deallocate_fq
      procedure :: deallocate_gs_matrices           => &
                   deallocate_gs_matrices_fq
      procedure :: deallocate_static_field_matrices => &
                   deallocate_static_field_matrices_fq
      !construct matrices/vectors
      procedure :: construct_static_matrix        => &
                   construct_static_matrix_fq
      procedure :: construct_ground_state_rhs     => &
                   construct_ground_state_rhs_fq
      procedure :: construct_static_field_rhs     => &
                   construct_static_field_rhs_fq
      !printing & saving files
      procedure :: print_atomtypes                => &
                   print_atomtypes_fq
      procedure :: print_gs_variables             => &
                   print_gs_variables_fq
      procedure :: print_static_field_variables   => &
                   print_static_field_variables_fq
      !calculate properties
      procedure :: calculate_energy               => &
                   calculate_energy_fq
      procedure :: calculate_static_polar         => &
                   calculate_static_polar_fq
   end type fq_type
    
   type (fq_type), target, save :: fq

contains

!-------------------------------------------------------------------------------
!  SUBROUTINES FOR ASSIGNING & GETTING VARIABLES
!-------------------------------------------------------------------------------

   !> Subroutine for assigning the dimensions of the FQ model
   !!    In/Out : target_      -- FQ
   subroutine assign_model_dimensions_fq(target_)
  
      implicit none

      !input/output variables
      class(fq_type), intent(inout) :: target_
       
      target_%n_var = target_%n_atoms + target_%n_mol
      target_%n_q   = target_%n_atoms
  
   end subroutine assign_model_dimensions_fq

!-------------------------------------------------------------------------------
!  SUBROUTINES FOR ALLOCATIONS
!-------------------------------------------------------------------------------

   !> Subroutine for allocating the GS matrices
   !!    In/Out  : target_             -- FQ
   subroutine allocate_gs_matrices_fq(target_)
  
      implicit none
       
      !input/output variables
      class(fq_type), intent(inout) :: target_
  
      target_%n_rhsre = 1
      call mem_man%alloc(target_%matrix,target_%n_var,target_%n_var, &
                         "target_%matrix")
      call mem_man%alloc(target_%rhs, target_%n_var, 1, "target_%rhs")
      call mem_man%alloc(target_%variables, target_%n_var,1,"target_%variables")
       
   end subroutine allocate_gs_matrices_fq


   !> Subroutine for allocating the matrices for static field
   !!    In/Out  : target_             -- FQ
   subroutine allocate_static_field_matrices_fq(target_)
  
      implicit none
       
      !input/output variables
      class(fq_type), intent(inout) :: target_
  
      target_%n_rhsre = 3
      call mem_man%alloc(target_%matrix, target_%n_var, target_%n_var, &
                         "target_%matrix")
      call mem_man%alloc(target_%rhs, target_%n_var, 3, "target_%rhs")
      call mem_man%alloc(target_%variables, target_%n_var,3,"target_%variables")
      call mem_man%alloc(target_%polar, 3, 3, "target_%polar")
  
   end subroutine allocate_static_field_matrices_fq


   !> Subroutine for deallocating the target_ variables
   !!    In/Out : target_      -- target_type
   subroutine deallocate_fq(target_)

      implicit none
  
      !input/output variables
      class(fq_type), intent(inout) :: target_

      call target_%deallocate_atomtypes
      call target_%deallocate_parameters()
      call target_%deallocate_atomic_quantities()

   end subroutine deallocate_fq


   !> Subroutine for allocating the GS matrices
   !!    In/Out  : target_             -- FQ
   subroutine deallocate_gs_matrices_fq(target_)
  
      implicit none
       
      !input/output variables
      class(fq_type), intent(inout) :: target_
  
      call mem_man%dealloc(target_%matrix, "target_%matrix")
      call mem_man%dealloc(target_%rhs, "target_%rhs")
      call mem_man%dealloc(target_%variables, "target_%variables")
       
   end subroutine deallocate_gs_matrices_fq


   !> Subroutine for allocating the matrices for static field
   !!    In/Out  : target_             -- FQ
   subroutine deallocate_static_field_matrices_fq(target_)
  
      implicit none
       
      !input/output variables
      class(fq_type), intent(inout) :: target_
  
      call mem_man%dealloc(target_%matrix, "target_%matrix")
      call mem_man%dealloc(target_%rhs, "target_%rhs")
      call mem_man%dealloc(target_%variables, "target_%variables")
      call mem_man%dealloc(target_%polar, "target_%polar")
  
   end subroutine deallocate_static_field_matrices_fq

!-------------------------------------------------------------------------------
!  SUBROUTINES FOR CONSTRUCTING MATRICES/VECTORS
!-------------------------------------------------------------------------------

   !> Subroutine for constructing the static matrix for FQ models
   !!    In/Out  : target_      -- FQ
   subroutine construct_static_matrix_fq(target_)
  
      implicit none
       
      !input/output variables
      class(fq_type), intent(inout) :: target_

      !internal variables
      real(dp), dimension(:,:), allocatable :: tmp_qq
  
      !Tqq block
      call mem_man%alloc(tmp_qq, target_%n_q, target_%n_q, "tmp_qq")
      call construct_static_t_qq(target_, tmp_qq)
      target_%matrix(1:target_%n_q,1:target_%n_q) = tmp_qq
      call mem_man%dealloc(tmp_qq, "tmp_qq")
  
      !Lagrangian multipliers
      if(target_%n_var.ne.target_%n_atoms.and.target_%name_.ne.'wfqfmu') then
         call add_lagrangian_blocks(target_, target_%matrix)
      else
         if(out_%ivrb.ge.2) &
            write(out_%iunit,'(1x,a)') ' Charge Constraints not setted'
      endif
      !printing
      if(out_%ivrb.ge.3) &
         call out_%print_matrix('FQ Matrix',target_%matrix,&
                                 target_%n_var,target_%n_var)
       
   end subroutine construct_static_matrix_fq


   !> Subroutine for constructing the GS RHS 
   !!    In/Out  : target_      -- FQ
   subroutine construct_ground_state_rhs_fq(target_)
  
      implicit none
       
      !input/output variables
      class(fq_type), intent(inout) :: target_
  
      !internal variables
      integer :: i
      character(len=14) :: format_1 = "(14x,'RHS ',/)"
      character(len=15) :: format_2 = "(i4,(4x,e14.6))"
  
      !from n_atoms to n_var, rhs = 0 [Lagrangian block]
      !this should be changed if the MM molecules has charge != 0
      !$omp parallel do private(i)
      do i=1,target_%n_atoms
         target_%rhs(i,1) = -fq%chi(target_%map_atomtypes(i))
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
       
   end subroutine construct_ground_state_rhs_fq


   !> Subroutine for constructing the RHS when static field is applied
   !!    In/Out  : target_      -- FQ 
   subroutine construct_static_field_rhs_fq(target_)
  
      implicit none
       
      !input/output variables
      class(fq_type), intent(inout) :: target_
  
      call construct_static_field_rhs_q(target_, &
                                        target_%rhs(1:target_%n_atoms, :))
      !printing
      if(out_%ivrb.ge.3) call print_static_field_rhs(target_,target_%rhs)
       
   end subroutine construct_static_field_rhs_fq

!-------------------------------------------------------------------------------
!  SUBROUTINES FOR PRINTING & SAVING FILES
!-------------------------------------------------------------------------------

   !> Subroutine for printing the atomtypes
   !!    Input  : target_   -- FQ
   subroutine print_atomtypes_fq(target_)
       
      implicit none
  
      !input/output variables
      class(fq_type), intent(in) :: target_
  
      !internal variables
      integer :: i
       
      do i = 1, target_%n_atomtypes
         write(out_%iunit,out_%sticks) 
         write(out_%iunit,'(1x,a)') "AtomType           : "// &
                                      trim(target_%atom_type(i))
         write(out_%iunit,'(1x,a,f9.3,a)') "Chi                : ", &
                                            target_%chi(i), " a.u."
         write(out_%iunit,'(1x,a,f9.3,a)') "Eta                : ", &
                                            target_%eta(i), " a.u."
         write(out_%iunit,'(1x,a,f9.3,a)') "R_q                : ", &
                                            target_%r_q(i), " a.u."
      enddo
       
   end subroutine print_atomtypes_fq


   !> Subroutine for printing GS variables
   !!    In/Out  : target_      -- FQ
   subroutine print_gs_variables_fq(target_)
  
      implicit none
       
      !input/output variables
      class(fq_type), intent(in) :: target_
  
      !internal variables
      integer :: i
      real(dp) :: sum_charges = zero
      character(len=20) :: format_1 = "(/,14x,'Charges ',/)"
      character(len=18) :: format_2 = "(1x,I4,(4x,E14.6))"
      character(len=28) :: format_3 = "(/,1x,'ChErr',3X,(E14.6,4x))"
  
      write(out_%iunit,format_1) 
      do i = 1, target_%n_atoms
         write(out_%iunit,format_2) i,target_%variables(i,1)
         sum_charges = sum_charges + target_%variables(i,1)
      enddo
      write(out_%iunit,format_3) sum_charges
      write(out_%iunit,out_%sticks) 
  
   end subroutine print_gs_variables_fq


   !> Subroutine for printing the variables when static field is applied
   !!    In/Out  : target_      -- model 
   subroutine print_static_field_variables_fq(target_)
  
      implicit none
       
      !input/output variables
      class(fq_type), intent(in) :: target_
  
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
       
   end subroutine print_static_field_variables_fq

!-------------------------------------------------------------------------------
!  SUBROUTINES FOR CALCULATING PROPERTIES
!-------------------------------------------------------------------------------

   !> Subroutine for calculating the energy of the GS 
   !!    In/Out  : target_      -- FQ
   subroutine calculate_energy_fq(target_)
  
      implicit none
       
      !input/output variables
      class(fq_type), intent(inout) :: target_
  
      real(dp)          :: ddot
      character(len=14) :: format_1 = "(1x,a,f14.8,a)"
  
      !energy = -1/2 q * RHS
      target_%energy = - half * &
                         ddot(target_%n_var,target_%variables,1,target_%rhs,1)
      write(out_%iunit,format_1) 'Energy = ', target_%energy,' a.u.'
      write(out_%iunit,out_%sticks)
       
   end subroutine calculate_energy_fq


   !> Subroutine for calculating the static polar
   !!    In/Out  : target_      -- FQ
   subroutine calculate_static_polar_fq(target_)
  
      implicit none
       
      !input/output variables
      class(fq_type), intent(inout) :: target_
  
      !internal variables
      integer :: info
      integer :: lwork
      real(dp) :: PolarReIso
      integer, dimension(:), allocatable    :: ipiv
      real(dp), dimension(1)                :: work_dummy
      real(dp), dimension(:),   allocatable :: work
      real(dp), dimension(:,:), allocatable :: Jmatrix
      real(dp), dimension(:,:), allocatable :: intermediate_matrix
      real(dp), dimension(:,:), allocatable :: PolarRe_2
      real(dp), dimension(:,:), allocatable :: lang_mult_matrix
      real(dp), dimension(:,:), allocatable :: J_1
      real(dp), dimension(:,:), allocatable :: J_1_1
      real(dp), dimension(:,:), allocatable :: J_1_R
      real(dp), dimension(:,:), allocatable :: J_1_1_R
      !R_a (for FQ = coord, for FQFMu = coord - Tqmu_tmumu e_a
      real(dp), dimension(:,:), allocatable :: R_a
      character(len=35) :: format_1 = "(1x,'Polar Iso  = ',f19.10,' a.u.')"
      character(len=39) :: format_2 = "(1x,'Results for w = ', E9.3, ' a.u.')"
  
      !first calculation of J**-1
      call mem_man%alloc(jmatrix, target_%n_q, target_%n_q, "jmatrix")
      JMatrix = target_%matrix(1:target_%n_q,1:target_%n_q)

      !allocation of the IPIV matrix
      call mem_man%alloc(ipiv, target_%n_q, "ipiv")

      !LU factorization
      call dgetrf(target_%n_q,target_%n_q,jmatrix,target_%n_q,ipiv,info)
      if(info.ne.0) call out_%error('info.ne.0 in FQ dgetrf in &
                                    &calculate_static_polar')
      !get the dimension for inversion
      call dgetri(target_%n_q,jmatrix,target_%n_q,ipiv,work_dummy(1),-1,info)
      lwork = nint(work_dummy(1))
      call mem_man%alloc(work, lwork, "work")

      !inversion
      call dgetri(target_%n_q,jmatrix,target_%n_q,ipiv,work,lwork,info)
      if(info.ne.0) call out_%error('info.ne.0 in dgetri in FQ &
                                    &calculate_static_polar')
      
      !printing of the inverse matrix
      if(out_%ivrb.ge.3) &
         call out_%print_matrix('JMatrix-1',JMatrix,target_%n_q,target_%n_q)
      call mem_man%dealloc(work,"work")
      call mem_man%dealloc(ipiv,"ipiv")
  
      !calculation of -R_a**T J**-1 R_b (a,b = x,y,z)
      ! 1) allocation
      ! 2) definition of R_a
      ! 3) Matrix product with DGEMM
      call mem_man%alloc(intermediate_matrix, target_%n_q, 3, &
                         "intermediate_matrix")
      call mem_man%alloc(R_a, 3, target_%n_q, "R_a")
      !definition of R_a
      R_a = target_%coord
      !J**-1 R_b
      Call DGEMM('N','T',              &   
                  target_%n_q,         &
                  3,                   &
                  target_%n_q,         &
                  one,                 &
                  JMatrix,             &
                  target_%n_q,         &
                  R_a,                 &
                  3,                   &
                  zero,                &
                  Intermediate_matrix, &
                  target_%n_q)      
      !R_a* J**-1 R_b
      Call DGEMM('N','N',              &
                  3,                   &
                  3,                   &
                  target_%n_q,         &
                  one,                 &
                  R_a,                 &
                  3,                   &
                  Intermediate_matrix, &
                  target_%n_q,         &
                  zero,                &
                  target_%polar,       &
                  3)                   
  
      !Now, let's add the contribution by the Lagrangian blocks
      !calculation of 1**T J**-1 * R_a * 1**T J**-1 * R_b
      !new equation by T. Giovannini and M. Ambrosetti (2020) unpublished
      !1) 1**T J-1 1   : (target_%n_mol,NMol)
      !2) 1**T J-1 R_a : (target_%n_mol, 3)
      call mem_man%alloc(lang_mult_matrix, target_%n_q, target_%n_mol, &
                         "lang_mult_matrix")
      call mem_man%alloc(J_1_R, target_%n_mol, 3, "J_1_R")
      call mem_man%alloc(J_1, target_%n_q, target_%n_mol, "J_1")
      call mem_man%alloc(J_1_1, target_%n_mol, target_%n_mol, "J_1_1")
      call mem_man%alloc(PolarRe_2, 3, 3, "PolarRe_2")
      !Get the Lagrangian block (1:Natoms, Natoms+1:NVar)
      lang_mult_matrix = &
               target_%matrix(1:target_%n_q,target_%n_q+1:target_%n_var)
      !calculation of 1**T J**-1 * 1
      !J**-1 * 1
      Call DGEMM('N','N',            &
                  target_%n_q,   &
                  target_%n_mol,     &
                  target_%n_q,       &
                  one,               &
                  JMatrix,           &
                  target_%n_q,       &
                  lang_mult_matrix,  &
                  target_%n_q,       &
                  zero,              &
                  J_1,               &
                  target_%n_q) 
      !1**T J**-1 * 1
      Call DGEMM('T','N',            &
                   target_%n_mol,    &
                   target_%n_mol,    &
                   target_%n_q,      &
                   one,              &
                   lang_mult_matrix, &
                   target_%n_q,      &
                   J_1,              &
                   target_%n_q,      &
                   zero,             &
                   J_1_1,            &
                   target_%n_mol)    ! 1**T J**-1 * 1
  
      !2) 1**T J-1 R_a : (target_%n_mol, 3)
      call mem_man%alloc(ipiv, target_%n_mol, "ipiv")
      !LU decomposition of J_1_1
      call dgetrf(target_%n_mol,target_%n_mol,J_1_1,target_%n_mol,ipiv,info)
      if(info.ne.0) call out_%error('info.ne.0 in dgetrf @2 in FQ &
                                    &calculate_static_polar')
      !Get the dimension of inversion
      call dgetri(target_%n_mol,J_1_1,target_%n_mol,ipiv,work_dummy(1),-1,info)
      lwork = nint(work_dummy(1))
      call mem_man%alloc(work, lwork, "work")
      !inversion of J_1_1
      call dgetri(target_%n_mol,J_1_1,target_%n_mol,ipiv,work,lwork,info)
      if(info.ne.0) call out_%error('info.ne.0 in dgetri @2 in FQ &
                                    &calculate_static_polar')
      !printing
      If(out_%ivrb.ge.3) &
         call out_%print_matrix('J_1_1 in calculate_static_polar',J_1_1, &
                                target_%n_mol,target_%n_mol)
      call mem_man%dealloc(work, "work")
      call mem_man%dealloc(ipiv, "ipiv")
  
      !polarre_2: contribution due to Lagrangian blocks
      call mem_man%alloc(J_1_1_R, target_%n_mol, 3, "J_1_1_R")
      ! 1**T * J**-1 R_b
      Call DGEMM('T','N',              &
                  target_%n_mol,       &
                  3,                   &
                  target_%n_q,         &
                  one,                 &
                  lang_mult_matrix,    &
                  target_%n_q,         &
                  Intermediate_matrix, &
                  target_%n_q,         &
                  zero,                &
                  J_1_R,               &
                  target_%n_mol) 
      ! J_1_1**-1 J_1_R
      Call DGEMM('N','N',        &
                  target_%n_mol, &
                  3,             &
                  target_%n_mol, &
                  one,           &
                  J_1_1,         &
                  target_%n_mol, &
                  J_1_R,         &
                  target_%n_mol, &
                  zero,          &
                  J_1_1_R,       &
                  target_%n_mol)    
      ! 1**T * J**-1 R_b
      Call DGEMM('T','N',        &
                  3,             &
                  3,             &
                  target_%n_mol, &
                  one,           &
                  J_1_R,         &
                  target_%n_mol, &
                  J_1_1_R,       &
                  target_%n_mol, &
                  zero,          &
                  PolarRe_2,     &
                  3) 
      !update polar
      target_%polar = target_%polar - PolarRe_2
      !isotropic polar
      PolarReIso = (target_%polar(1,1)+target_%polar(2,2)+target_%polar(3,3))/3

      !printing in output
      write(out_%iunit,format_2) zero
      call out_%print_matrix('polarizability tensor (a.u.)',target_%polar,3,3)
      Write(out_%iunit,format_1) PolarReIso 
      Write(out_%iunit,out_%sticks) 
  
      !deallocations
      call mem_man%dealloc(Jmatrix, 'jmatrix')
      call mem_man%dealloc(intermediate_matrix, 'intermediate_matrix')
      call mem_man%dealloc(PolarRe_2, 'PolarRe_2')
      call mem_man%dealloc(lang_mult_matrix, 'lang_mult_matrix')
      call mem_man%dealloc(J_1, 'J_1')
      call mem_man%dealloc(J_1_1, 'J_1_1')
      call mem_man%dealloc(J_1_R, 'J_1_R')
      call mem_man%dealloc(J_1_1_R, 'J_1_1_R')
      call mem_man%dealloc(R_a, 'R_a')
  
   end subroutine calculate_static_polar_fq

end module fq_module
