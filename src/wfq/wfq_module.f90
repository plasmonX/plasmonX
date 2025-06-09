!> wFQ module
!!
!! This module contains the subroutines for the wfq type
!!
!! Date         : 2024 
!!
module wfq_module
 
   use output_module
   use parameters_module
   use string_manipulation_module
   use array_manipulation_module
   use field_module
   use target_module
   use fq_module
 
   implicit none
 
   !public variables
   public wfq
 
   type, extends(fq_type) :: wfq_type
      !K_minus_P matrix
      real(dp), dimension(:,:), allocatable :: K_minus_P  
      !RHS static field
      real(dp), dimension(:,:), allocatable :: rhs_static 
      ! H_minus_L matrix -- for alloys
      complex(dp), dimension(:,:), allocatable :: H_minus_L  
        
      contains 
    
      !assign/get variables
      procedure :: assign_model_dimensions                => &
                   assign_model_dimensions_wfq
      procedure :: assign_model_parameters                => &
                   assign_model_parameters_wfq
      !allocations
      procedure :: allocate_dynamic_field_memory          => &
                   allocate_dynamic_field_memory_wfq
      procedure :: deallocate_dynamic_field_memory        => &
                   deallocate_dynamic_field_memory_wfq
      !construct matrices/vectors
      procedure :: construct_constant_matrix              => &
                   construct_constant_matrix_wfq
      procedure :: construct_dynamic_matrix               => &
                   construct_dynamic_matrix_wfq
      procedure :: construct_dynamic_field_rhs            => &
                   construct_dynamic_field_rhs_wfq
      procedure :: construct_dynamic_field_rhs_on_the_fly => &
                   construct_dynamic_field_rhs_on_the_fly_wfq
      procedure :: product_matrix_vector                  => &
                   product_matrix_vector_wfq
      procedure :: gmres_diagonal_shift                   => &
                   gmres_diagonal_shift_wfq
      procedure :: new_gmres_diagonal_shift               => &
                   new_gmres_diagonal_shift_wfq
      !printing & saving files
      procedure :: print_atomtypes                        => &
                   print_atomtypes_wfq
      procedure :: print_dynamic_field_variables          => &
                   print_dynamic_field_variables_wfq
      !calculate properties
      procedure :: calculate_dynamic_polar                => &
                   calculate_dynamic_polar_wfq
      procedure :: calculate_induced_field_at_point       => &
                   calculate_induced_field_at_point_wfq
      procedure :: calculate_density_at_point             => &
                   calculate_density_at_point_wfq      
   end type wfq_type
    
   type (wfq_type), target, save :: wfq

contains

!-------------------------------------------------------------------------------
!  SUBROUTINES FOR ASSIGNING & GETTING VARIABLES
!-------------------------------------------------------------------------------

   !> Subroutine for assigning the dimensions of the wFQ model
   !!    In/Out : target_      -- wFQ
   subroutine assign_model_dimensions_wfq(target_)
  
      implicit none

      !input/output variables
      class(wfq_type), intent(inout) :: target_
       
      target_%n_var = target_%n_atoms 
      target_%n_q   = target_%n_atoms
  
   end subroutine assign_model_dimensions_wfq


   !> Subroutine for assigning the wfq parameters (sigma_0, tau, density)
   !!    In/Out : target_      -- wfq
   subroutine assign_model_parameters_wfq(target_)
  
      implicit none
       
      !input/output variables
      class(wfq_type), intent(inout) :: target_
      !internal variables
      integer :: i
  
      call mem_man%alloc(parameters%density, target_%n_atoms, &
                         "parameters%density")
      do i = 1, target_%n_atomtypes
         if(is_equal(target_%atomic_number(i),6.0d0)) then !this is graphene-based
            !the Fermi energy is in input
            !Considering atomic units 
            ! [See eqs. 1-3 in JPCC 2023, 127, 611 - 10.1021/acs.jpcc.3c01565]
            !
            ! E_F = v_F * sqrt(pi * n_2D) --> n_2D = E_F^2 / ( v_F^2 * pi )
            !
            ! n = sqrt(n_2D) * v_F / sqrt(pi) --> 
            ! [ E_F / (v_F * sqrt(pi) ) ] * v_F / sqrt(pi) --> E_F / pi
            parameters%density(i) = parameters%fermi_energy(i)/(pi*au_to_ev)
         else
            parameters%density(i) = parameters%sigma_0(i)/parameters%tau(i)
         endif
         if(.not.is_equal(parameters%scaling(i),one)) then !this is used for sodium 
            parameters%sigma_0(i) = parameters%sigma_0(i)/parameters%scaling(i)
            parameters%tau(i)     = parameters%tau(i)/parameters%scaling(i)
         endif
      enddo 
  
   end subroutine assign_model_parameters_wfq

!-------------------------------------------------------------------------------
!  SUBROUTINES FOR ALLOCATIONS
!-------------------------------------------------------------------------------

   !> Subroutine for allocating dynamic field variables for memory (iter/inv)
   !!    In/Out  : target_             -- wfq
   !!    In/Out  : matrix_constant     -- part of matrix_w that is constant
   !!    In/Out  : matrix_w            -- dynamic Matrix
   subroutine allocate_dynamic_field_memory_wfq(target_,matrix_constant, &
                                                matrix_w)
  
      implicit none
       
      !input/output variables
      class(wfq_type)    :: target_
      real(dp), dimension(:,:), allocatable    :: matrix_constant
      complex(dp), dimension(:,:), allocatable :: matrix_w
  
      target_%n_rhs = 3
      call mem_man%alloc(matrix_w, target_%n_var, target_%n_var, "matrix_w") 
      call mem_man%alloc(matrix_constant, target_%n_var, target_%n_var, &
                        "matrix_constant") 
      call mem_man%alloc(target_%K_minus_P, target_%n_q,target_%n_q, &
                         "target_%K_minus_P")
      if(target_%heterogeneous) &
         !allocation of H-L [eq. 16 - Front. Photon. 2023, 4, 1199598] 
         call mem_man%alloc(target_%H_minus_L, target_%n_q, target_%n_q, &
                            "target_%H_minus_L")
       
   end subroutine allocate_dynamic_field_memory_wfq


   !> Subroutine for allocating dynamic field variables for memory (iter/inv)
   !!    In/Out  : target_             -- wfq
   !!    In/Out  : matrix_constant     -- part of matrix_w that is constant
   !!    In/Out  : matrix_w            -- dynamic Matrix
   subroutine deallocate_dynamic_field_memory_wfq(target_,matrix_w, &
                                                matrix_constant)
  
      implicit none
       
      !input/output variables
      class(wfq_type)    :: target_
      complex(dp), dimension(:,:), allocatable :: matrix_w
      real(dp), dimension(:,:), allocatable, optional :: matrix_constant
  
      call mem_man%dealloc(matrix_w, "matrix_w") 
      if(present(matrix_constant)) &
         call mem_man%dealloc(matrix_constant,"matrix_constant") 
      call mem_man%dealloc(target_%K_minus_P,"target_%K_minus_P")
      if(target_%heterogeneous) &
         !allocation of H-L [eq. 16 - Front. Photon. 2023, 4, 1199598] 
         call mem_man%dealloc(target_%H_minus_L,"target_%H_minus_L")
       
   end subroutine deallocate_dynamic_field_memory_wfq

!-------------------------------------------------------------------------------
!  SUBROUTINES FOR CONSTRUCTING MATRICES/VECTORS
!-------------------------------------------------------------------------------

   !> Subroutine for constructing the constant matrix
   !!    In/Out : target_         -- wfq
   !!    In/Out : matrix_constant -- constant part of the dynamic matrix
   subroutine construct_constant_matrix_wfq(target_,matrix_constant)
  
      implicit none
       
      !input/output variables
      class(wfq_type), intent(inout) :: target_
      real(dp), dimension(target_%n_var,target_%n_var), intent(inout) :: &
                                                                 matrix_constant
  
      ! K-P matrix : see Eq. 17 in J. Phys. Chem. C 2021, 125, 23848-23863
      call construct_K_minus_P_matrix(target_, target_%K_minus_P)
      if(.not.target_%heterogeneous) &
         call construct_Aqq_matrix(target_, target_%K_minus_P, matrix_constant)
  
   end subroutine construct_constant_matrix_wfq


   !> Subroutine for constructing the wfq dynamic matrix
   !!    In/Out : target_         -- wfq
   !!    In/Out : i_freq          -- index of frequency
   !!    Input  : matrix_constant -- constant part of the dynamic matrix
   !!    In/Out : matrix_w        -- dynamic matrix
   subroutine construct_dynamic_matrix_wfq(target_,i_freq,matrix_constant, &
                                           matrix_w)
  
      implicit none
       
      !input/output variables
      class(wfq_type), intent(inout) :: target_
      integer, intent(in)            :: i_freq
      real(dp), dimension(target_%n_var,target_%n_var), intent(in)    :: &
                                                                 matrix_constant
      complex(dp), dimension(target_%n_var,target_%n_var), intent(inout) :: &
                                                                        matrix_w
  
      !internal variables
      integer           :: i,j
      real(dp)          :: freq_au
      complex(dp)       :: z_w 
  
      freq_au = field%freq(i_freq)

      do i = 1, target_%n_var
         do j = 1, target_%n_var
            matrix_w(i,j) = dcmplx(matrix_constant(i,j),zero)
         enddo
         !z(w) = - \frac{ w / (2 * n_0 * tau) * ( w * tau + i ) }
         z_w = - ( freq_au /(two*parameters%density(target_%map_atomtypes(i))* &
                                 parameters%tau(target_%map_atomtypes(i))) ) * &
                    dcmplx(freq_au*parameters%tau(target_%map_atomtypes(i)),one)
         !matrix_w(i,i) = matrix_w(i,i) - z(w)
         matrix_w(i,i) = matrix_w(i,i) - z_w
  
      enddo
      if(out_%ivrb.ge.3) then
         call out_%print_matrix("wFQ Matrix: Real Part",dble(matrix_w), &
                                target_%n_var,target_%n_var)
         call out_%print_matrix("wFQ Matrix: IMag Part",dimag(matrix_w),&
                                target_%n_var,target_%n_var)
      endif
  
   end subroutine construct_dynamic_matrix_wfq


   !> Subroutine for constructing the wfq RHS on memory
   !!    In/Out : target_         -- wfq
   !!    In/Out : rhs_w           -- dynamic RHS
   subroutine construct_dynamic_field_rhs_wfq(target_,rhs_w)
  
      implicit none
       
      !input/output variables
      class(wfq_type), intent(inout) :: target_
      complex(dp), dimension(target_%n_var,3), intent(inout) :: rhs_w
  
      target_%n_rhsre = 3
      call mem_man%alloc(target_%rhs_static,target_%n_var,3,"target_%rhs_static")
      call construct_static_field_rhs_q(target_, &
                                        target_%rhs_static(1:target_%n_atoms,:))
      call construct_rhs_w_q(target_,target_%rhs_static,target_%K_minus_P,rhs_w)
      
      if (out_%ivrb.ge.3) call print_dynamic_field_rhs(target_,rhs_w)
      call mem_man%dealloc(target_%rhs_static,"target_%rhs_static")
  
   end subroutine construct_dynamic_field_rhs_wfq


   !> Subroutine for constructing the wfq RHS on the fly
   !!    In/Out : target_         -- wfq
   !!    In/Out : rhs_w           -- dynamic RHS
   subroutine construct_dynamic_field_rhs_on_the_fly_wfq(target_,rhs_w)
  
      implicit none
       
      !input/output variables
      class(wfq_type), intent(inout) :: target_
      complex(dp), dimension(target_%n_var,3), intent(inout) :: rhs_w
  
      target_%n_rhsre = 3
      call mem_man%alloc(target_%rhs_static,target_%n_var,3, &
                         "target_%rhs_static")
      call construct_static_field_rhs_q(target_, &
                                        target_%rhs_static(1:target_%n_atoms,:))
      call construct_rhs_w_q_on_the_fly(target_, target_%rhs_static, rhs_w)

      if (out_%ivrb.ge.3) call print_dynamic_field_rhs(target_,rhs_w)
      call mem_man%dealloc(target_%rhs_static, "target_%rhs_static")
  
   end subroutine construct_dynamic_field_rhs_on_the_fly_wfq


   !> Subroutine for performing the matrix vector multiplication on the fly
   !!     
   !!    La matrice A ha la seguente forma
   !!       A_{ij} = \sum_k M_{ik} (-D_{ij+D_{kj})
   !!    e definendo la matrice diagonale
   !!       P_{ij} = \sum_k K_{ik}\delta_{ij}
   !!    può essere espressa come prodotto di due matrici simmetriche
   !!       A_{ij} = \sum_k (K-P)_{ik} D_{kj}
   !!    Quindi prima calcoliamo il prodotto 
   !!       t = D*x (product_q_block_x)
   !!    Il vettore risultante e' chiamato static_mat_dot_x 
   !!       y=(M-P)*t (update y_q)
   !!
   !!    In/Out  : target_   -- model 
   !!    Input   : i_freq    -- index of frequency
   !!    In/Out  : x         -- input vector 
   !!    In/Out  : y         -- output vector 
   subroutine product_matrix_vector_wfq(target_,i_freq, x, y)
  
      implicit none
       
      !input/output variables
      class(wfq_type), intent(inout) :: target_
      integer, intent(in)            :: i_freq
      complex(dp), dimension(target_%n_var), intent(inout) :: x
      complex(dp), dimension(target_%n_var), intent(inout) :: y
  
      !internal variables
      integer :: i
      real(dp) :: freq_au
      complex(dp), dimension(:), allocatable :: static_matrix_dot_x
  
      call mem_man%alloc(static_matrix_dot_x,target_%n_var, &
                         "static_matrix_dot_x")
      freq_au = field%freq(i_freq)
      !$omp parallel do
      do i=1,target_%n_var
         y(i) = dcmplx(zero,zero)
      enddo
      !$omp end parallel do
      call product_q_block_x(target_,x,static_matrix_dot_x)
      call update_y_q(target_, freq_au, static_matrix_dot_x, y) 
      call mem_man%dealloc(static_matrix_dot_x,"static_matrix_dot_x")
       
   end subroutine product_matrix_vector_wfq


   !> Subroutine for applying the diagonal shift for GMRES algorithm
   !!    In/Out  : target_      -- model 
   !!    Input   : i_freq       -- index of the frequency
   !!    Input   : colx         -- index of the column x (GMRES)
   !!    Input   : colz         -- index of the column z (GMRES)
   !!    Input   : lwork        -- dimension of the work array
   !!    In/Out  : work         -- array work
   subroutine gmres_diagonal_shift_wfq(target_,i_freq,colx,colz,lwork,work)
  
      implicit none
       
      !input/output variables
      class(wfq_type)    :: target_
      integer  :: i_freq
      integer  :: colx
      integer  :: colz
      integer  :: lwork
      complex(dp), dimension(lwork), intent(inout) :: work
  
      !internal variables
      integer  :: i
      integer  :: start_x
      integer  :: start_z
      real(dp) :: freq_au
      complex(dp) :: z_i
  
      freq_au = field%freq(i_freq)
  
      start_z = colz - 1
      start_x = colx - 1
      do i = 1, target_%n_var
         z_i = freq_au*(dcmplx( -freq_au* &
                                parameters%tau(target_%map_atomtypes(i)),  &
                                -one ))/ & 
                       (two*parameters%density(target_%map_atomtypes(i)) * &
                            parameters%tau(target_%map_atomtypes(i)))
         work(start_z+i) = work(start_z + i) - z_i * work(start_x + i)
      enddo
       
   end subroutine gmres_diagonal_shift_wfq


   !> Subroutine for applying the diagonal shift for GMRES algorithm
   !!    In/Out  : target_      -- model 
   !!    Input   : i_freq       -- index of the frequency
   !!    Input   : r            -- rhs
   !!    In/Out  : w            -- solution
   subroutine new_gmres_diagonal_shift_wfq(target_,i_freq,r,w)
  
      implicit none
       
      !input/output variables
      class(wfq_type)    :: target_
      integer  :: i_freq
      complex(dp), dimension(target_%n_var), intent(in) :: r
      complex(dp), dimension(target_%n_var), intent(inout) :: w
  
      !internal variables
      integer  :: i
      real(dp) :: freq_au
      complex(dp) :: z_i
  
      freq_au = field%freq(i_freq)
  
      do i = 1, target_%n_var
         z_i = freq_au*(dcmplx( -freq_au* &
                                parameters%tau(target_%map_atomtypes(i)),  &
                                -one ))/ & 
                       (two*parameters%density(target_%map_atomtypes(i)) * &
                            parameters%tau(target_%map_atomtypes(i)))
         w(i) = w(i) - z_i * r(i)
      enddo
       
   end subroutine new_gmres_diagonal_shift_wfq

!-------------------------------------------------------------------------------
!  SUBROUTINES FOR PRINTING & SAVING FILES
!-------------------------------------------------------------------------------

   !> Subroutine for printing the atomtypes
   !!    Input  : target_   -- model 
   subroutine print_atomtypes_wfq(target_)
      
      implicit none
  
      !input/output variables
      class(wfq_type), intent(in) :: target_
  
      !internal variables
      integer            :: i
       
      do i = 1, target_%n_atomtypes
         write(out_%iunit,out_%sticks) 
         write(out_%iunit,'(1x,a)') "AtomType           : "// &
                                    trim(target_%atom_type(i))
         write(out_%iunit,'(1x,a,f10.3,a)') "Chi                : ", &
                                            target_%chi(i), " a.u."
         write(out_%iunit,'(1x,a,f10.3,a)') "Eta                : ", &
                                            target_%eta(i), " a.u."
         write(out_%iunit,'(1x,a,f10.3,a)') "R_q                : ", &
                                            target_%r_q(i), " a.u."
  
         write(out_%iunit,'(/,1x,a,e10.3,a)') "tau                : ", &
                                               parameters%tau(i)," a.u."
         write(out_%iunit,'(1x,a,e10.3,a)')   "sigma_0            : ", &
                                               parameters%sigma_0(i)," a.u."
         if(allocated(parameters%scaling)) &
            write(out_%iunit,'(1x,a,e10.3,a)')   "scaling sigma0-tau : ", &
                                                  parameters%scaling(i)," a.u."
         write(out_%iunit,'(1x,a,e10.3,a)')   "A_ij               : ", &
                                               parameters%A_ij(i)," a.u."
         write(out_%iunit,'(1x,a,e10.3,a)')   "fermi_d            : ", &
                                               parameters%fermi_d(i,i)," a.u."
         write(out_%iunit,'(1x,a,e10.3,a)')   "fermi_s            : ", &
                                               parameters%fermi_s(i,i)," a.u."
         if(trim(target_%atom_type(i)).eq.'C') &
            write(out_%iunit,'(1x,a,e10.3,a)')   "Fermi energy       : ", &
                                               parameters%fermi_energy(i), " eV"
         if(allocated(parameters%density)) &
            write(out_%iunit,'(1x,a,e10.3,a)')   "density            : ", &
                                               parameters%density(i), " a.u."
      enddo
       
   end subroutine print_atomtypes_wfq


   !> Subroutine for printing the variables when dynamic field is applied
   !!    Input   : target_      -- model 
   !!    Input   : variables_w  -- w-variables (N_var, 3 [x,y,z])
   subroutine print_dynamic_field_variables_wfq(target_,variables_w)
  
      implicit none
       
      !input/output variables
      class(wfq_type),intent(in) :: target_
      complex(dp), dimension(target_%n_var,3), intent(in) :: variables_w
  
      !internal variables
      integer :: i
      complex(dp) :: sum_X = dcmplx(zero,zero)
      complex(dp) :: sum_Y = dcmplx(zero,zero)
      complex(dp) :: sum_Z = dcmplx(zero,zero)
      character(len=21) :: format_1 = "(38x,'wFQ Charges',/)"
      character(len=64) :: format_2 = "(13x,'X Component',14x,'Y Component',&
                                      &14x,'Z Component',/)"
      character(len=60) :: format_3 = "(1x,I4,2x,2(E10.3,2x,E10.3,1X,'|',1x),&
                                      &E10.3,2X,E10.3)"
      character(len=64) :: format_4 = "(/,'ChErr',2x,2(E10.3,2x,E10.3,1X,'|',&
                                      &1x),E10.3,2X,E10.3)"
  
      write(out_%iunit,format_1)
      write(out_%iunit,format_2)
      do i = 1, target_%n_atoms
         write(out_%iunit,format_3) i,variables_w(i,:)
         sum_x = sum_x + variables_w(i,1)
         sum_y = sum_y + variables_w(i,2)
         sum_z = sum_z + variables_w(i,3)
      enddo
      write(out_%iunit,format_4) sum_x, sum_y, sum_z
      write(out_%iunit,out_%sticks) 
       
   end subroutine print_dynamic_field_variables_wfq

!-------------------------------------------------------------------------------
!  SUBROUTINES FOR CALCULATING PROPERTIES
!-------------------------------------------------------------------------------

   !> Subroutine for calculating the dynamic polar
   !!    In/Out  : target_      -- model 
   !!    Input   : i_freq       -- index of the frequency
   !!    Input   : variables_w  -- w-variables
   !!    Output  : polar_w      -- dynamic polar
   subroutine calculate_dynamic_polar_wfq(target_,i_freq,variables_w,polar_w)
  
      implicit none
       
      !input/output variables
      class(wfq_type), intent(inout) :: target_
      integer, intent(in)            :: i_freq
      complex(dp), dimension(target_%n_var,3), intent(in) :: variables_w
      complex(dp), dimension(3,3), intent(inout)          :: polar_w
  
      !internal variables
      integer :: i
      complex(dp), dimension(3,3) :: internal_polar
  
      internal_polar = dcmplx(zero,zero)
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
  
   end subroutine calculate_dynamic_polar_wfq


   !> Function for calculating the induced field at a specific point
   !!    Input   : target_      -- model 
   !!    Input   : i_pol        -- index of the polarization
   !!    Input   : point_coord  -- coordinates of the point
   !!    Input   : variables_w  -- w-variables
   !!    Output  : EField       -- Electric Field
   function calculate_induced_field_at_point_wfq(target_,i_pol,point_coord, &
                                                 variables_w) result(EField)
  
      implicit none
       
      !input/output variables
      class(wfq_type), intent(in) :: target_
      integer, intent(in) :: i_pol
      real(dp), dimension(3), intent(in) :: point_coord
      complex(dp), dimension(target_%n_var), intent(in) :: variables_w
      complex(dp), dimension(3) :: EField
  
      !internal variables
      integer  :: i
      real(dp), dimension(3) :: d_IJ
       
      EField = dcmplx(zero, zero)
      EField(i_pol) = dcmplx(field%e_0, zero)
      do i = 1, target_%n_var
         d_IJ(1)   = (target_%coord(1,i)-point_coord(1))*tobohr
         d_IJ(2)   = (target_%coord(2,i)-point_coord(2))*tobohr
         d_IJ(3)   = (target_%coord(3,i)-point_coord(3))*tobohr
         call calculate_induced_field_q(target_,i,d_IJ,variables_w(i),EField)
      enddo                  
       
   end function calculate_induced_field_at_point_wfq


   !> Function for calculating the plasmon density at a specific point
   !!    Input   : target_      -- model 
   !!    Input   : point_coord  -- coordinates of the point
   !!    Input   : variables_w  -- w-variables
   !!    Output  : density      -- Plasmon density
   function calculate_density_at_point_wfq(target_,point_coord, &
                                           variables_w) result(density)
  
      implicit none
       
      !input/output variables
      class(wfq_type), intent(in)        :: target_
      real(dp), dimension(3), intent(in) :: point_coord
      complex(dp), dimension(target_%n_var), intent(in) :: variables_w
      complex(dp) :: density
  
      !internal variables
      integer  :: i
      real(dp), dimension(3) :: d_IJ
       
      density = zero
      do i = 1, target_%n_var
         d_IJ(1)   = (target_%coord(1,i)-point_coord(1))*tobohr
         d_IJ(2)   = (target_%coord(2,i)-point_coord(2))*tobohr
         d_IJ(3)   = (target_%coord(3,i)-point_coord(3))*tobohr
         call calculate_density_q(target_, i, d_IJ, variables_w(i), density)
      enddo                  
       
   end function calculate_density_at_point_wfq

end module wfq_module
