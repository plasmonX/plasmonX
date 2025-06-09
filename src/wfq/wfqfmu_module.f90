!> wFQFmu module
!!
!! This module contains the subroutines for the wfqfmu type
!!
!! Date         : 2025
!!
module wfqfmu_module
      
   use output_module
   use parameters_module
   use string_manipulation_module
   use array_manipulation_module
   use field_module
   use target_module
   use fqfmu_module
   use wfq_module

   implicit none

   !public variables
   public wfqfmu

   type, extends(wfq_type) :: wfqfmu_type
      integer :: n_dipoles
         
      contains 

      !assign/get variables
      procedure :: assign_model_dimensions                 => &
                   assign_model_dimensions_wfqfmu
      !assign/get variables
      procedure :: assign_variables_w                     => &
                   assign_variables_w_wfqfmu              
      procedure :: assign_variables_w_on_the_fly          => &
                   assign_variables_w_on_the_fly_wfqfmu             
      !construct matrices/vectors
      procedure :: construct_constant_matrix              => &
                   construct_constant_matrix_wfqfmu
      procedure :: construct_dynamic_field_rhs            => &
                   construct_dynamic_field_rhs_wfqfmu
      procedure :: construct_dynamic_field_rhs_on_the_fly => &
                   construct_dynamic_field_rhs_on_the_fly_wfqfmu
      procedure :: construct_dynamic_matrix               => &
                   construct_dynamic_matrix_wfqfmu
      procedure :: construct_dynamic_matrix_gmres         => &
                   construct_dynamic_matrix_gmres_wfqfmu
      procedure :: product_matrix_vector                  => &
                   product_matrix_vector_wfqfmu
      procedure :: gmres_diagonal_shift                   => &
                   gmres_diagonal_shift_wfqfmu
      procedure :: new_gmres_diagonal_shift               => &
                   new_gmres_diagonal_shift_wfqfmu
      !printing & saving files
      procedure :: print_dynamic_field_variables          => &
                   print_dynamic_field_variables_wfqfmu
      procedure :: print_atomtypes                        => &
                   print_atomtypes_wfqfmu
      !calculate properties
      procedure :: calculate_dynamic_polar                => &
                   calculate_dynamic_polar_wfqfmu
      procedure :: calculate_induced_field_at_point       => &
                   calculate_induced_field_at_point_wfqfmu
      procedure :: calculate_density_at_point             => &
                   calculate_density_at_point_wfqfmu      
   end type wfqfmu_type
    
   type (wfqfmu_type), target, save :: wfqfmu

contains

!-------------------------------------------------------------------------------
!  SUBROUTINES FOR ASSIGNING & GETTING VARIABLES
!-------------------------------------------------------------------------------

   !> Subroutine for assigning the dimensions of the wFQFMu model
   !!    In/Out : target_      -- wFQFMu
   subroutine assign_model_dimensions_wfqfmu(target_)
  
      implicit none

      !input/output variables
      class(wfqfmu_type), intent(inout) :: target_
       
      target_%n_var = 4*target_%n_atoms 
      target_%n_q   = target_%n_atoms
      target_%n_mu  = 3*target_%n_atoms
  
   end subroutine assign_model_dimensions_wfqfmu

!-------------------------------------------------------------------------------
!  SUBROUTINES FOR wFQFMu PARAMETERS
!-------------------------------------------------------------------------------

   !> Subroutine for reading alpha(w) from file or tabulated data
   subroutine read_alpha()
  
      use field_module
       
      implicit none
       
      !internal variables
      integer :: i
  
      call mem_man%alloc(parameters%alpha_w, target_%n_atomtypes, field%n_freq, &
                         "parameters%alpha_w")
  
      do i = 1, target_%n_atomtypes
         if(trim(parameters%permittivity_type(i)).ne.'none') then
            call set_fitted_permittivity(i)
         else if(trim(parameters%wfqfmu_file(i)).ne.'none') then
            call set_permittivity_from_file(i)
         endif
     enddo
  
   end subroutine read_alpha


   !> Subroutine for setting alpha(w) from fitted tabulated data
   !!    Input  : i           -- atomtype
   subroutine set_fitted_permittivity(i)
       
      implicit none
       
      !input/output variables
      integer, intent(in)               :: i
  
      if(parameters%permittivity_type(i).eq.'silver etchegoin') then
         call silver_etchegoin(out_%iunit,      &
                               'polar',         &
                               field%n_freq,   &
                               field%freq,     &
                               parameters%alpha_w(i,:))
      else if(parameters%permittivity_type(i).eq.'gold etchegoin') then
         call gold_etchegoin(out_%iunit,      &
                             'polar',         &
                             field%n_freq,   &
                             field%freq,     &
                             parameters%alpha_w(i,:))
      else 
         call out_%error("Fitted permittivity for material: "//&
                         trim(parameters%permittivity_type(i))// &
                         " not yet implemented")
      endif
  
   end subroutine set_fitted_permittivity


   !> Subroutine for setting alpha(w) from data reported in csv file
   !! specified in the input file
   !!    Input  : i           -- atomtype
   subroutine set_permittivity_from_file(i)
       
      implicit none
       
      !input variables
      integer, intent(in)               :: i 
  
      call read_permittivity_info_from_file(out_%iunit,                &
                                            parameters%wfqfmu_file(i), &
                                            field%n_freq,              &
                                            field%freq,                &
                                            parameters%alpha_w(i,:))
    
   end subroutine set_permittivity_from_file

!-------------------------------------------------------------------------------
!  SUBROUTINES FOR ASSIGNING & GETTING VARIABLES
!-------------------------------------------------------------------------------

   !> Subroutine for assigning the dynamic variables
   !!   Default : variables_w = rhs_w. This changes for heterogeneous
   !!    Input  : target_      -- model 
   !!    In/Out : rhs_w        -- dynamic RHS
   !!    Output : variables_w  -- dynamic variables
   subroutine assign_variables_w_wfqfmu(target_, rhs_w, variables_w)
       
      implicit none
  
      !input/output variables
      class(wfqfmu_type), intent(in)   :: target_
      complex(dp), dimension(target_%n_var,3), intent(inout) :: rhs_w
      complex(dp), dimension(target_%n_var,3), intent(out)   :: variables_w
  
      !internal variables
      complex(dp), dimension(:,:), allocatable :: rhs_0_complex
  
      if(target_%heterogeneous) then
         call mem_man%alloc(rhs_0_complex, target_%n_q, 3, "rhs_0_complex")
         !rhs_0_complex = rhs_w(1:target_%n_q,:)
         call zcopy(target_%n_q*3, rhs_w(1:target_%n_q,:), 1, rhs_0_complex, 1)
         call construct_rhs_w_q_complex(target_,half,rhs_0_complex, &
                                        target_%H_minus_L,          &
                                        rhs_w(1:target_%n_q,:))
         if (out_%ivrb.ge.3) call print_dynamic_field_rhs(target_,rhs_w)
         call mem_man%dealloc(rhs_0_complex, "rhs_0_complex")
      endif
  
      !variables_w = rhs_w
      call zcopy(target_%n_var*3, rhs_w, 1, variables_w, 1)
  
   end subroutine assign_variables_w_wfqfmu
   

   !> Subroutine for assigning the dynamic variables on the fly algorithm
   !!   Default : variables_w = rhs_w. This changes for heterogeneous
   !!    Input  : target_      -- model 
   !!    Input  : i_freq       -- index of the frequency
   !!    In/Out : rhs_w        -- dynamic RHS
   !!    Output : variables_w  -- dynamic variables
   subroutine assign_variables_w_on_the_fly_wfqfmu(target_, i_freq, rhs_w, &
                                                   variables_w)
       
      implicit none
  
      !input/output variables
      class(wfqfmu_type), intent(in)   :: target_
      integer, intent(in) :: i_freq
      complex(dp), dimension(target_%n_var,3), intent(inout) :: rhs_w
      complex(dp), dimension(target_%n_var,3), intent(out)   :: variables_w
  
      !internal variables
      real(dp) :: freq_au
  
      if(target_%heterogeneous) then
         freq_au = field%freq(i_freq)
         call construct_rhs_w_q_complex_on_the_fly(target_,half,freq_au, &
                                                   rhs_w(1:target_%n_q,:),&
                                                   rhs_w(1:target_%n_q,:))
         if (out_%ivrb.ge.3) call print_dynamic_field_rhs(target_,rhs_w)
      endif
  
      !variables_w = rhs_w !default
      call zcopy(target_%n_var*3, rhs_w, 1, variables_w, 1)
  
   end subroutine assign_variables_w_on_the_fly_wfqfmu

!-------------------------------------------------------------------------------
!  SUBROUTINES FOR CONSTRUCTING MATRICES/VECTORS
!-------------------------------------------------------------------------------

   !> Subroutine for constructing the constant matrix
   !!    In/Out : target_         -- wfqfmu
   !!    In/Out : matrix_constant -- constant part of the dynamic matrix
   subroutine construct_constant_matrix_wfqfmu(target_,matrix_constant)
  
      implicit none
       
      !input/output variables
      class(wfqfmu_type), intent(inout) :: target_
      real(dp), dimension(target_%n_var,target_%n_var), intent(inout) :: &
                                                                 matrix_constant
      !internal variables
      real(dp), dimension(:,:), allocatable :: tmp_qq
      real(dp), dimension(:,:), allocatable :: tmp_qmu
      real(dp), dimension(:,:), allocatable :: tmp_muq
      real(dp), dimension(:,:), allocatable :: tmp_mumu
  
      !wFQ block
      call construct_K_minus_P_matrix(target_, target_%K_minus_P)
      if(.not.target_%heterogeneous) then
         call mem_man%alloc(tmp_qq, target_%n_q, target_%n_q, "tmp_qq")
         call construct_Aqq_matrix(target_, target_%K_minus_P, tmp_qq)
         matrix_constant(1:target_%n_q,1:target_%n_q) = tmp_qq
         call mem_man%dealloc(tmp_qq, "tmp_qq")
      endif

      !wFMu - wFQ Block
      call mem_man%alloc(tmp_muq, target_%n_mu, target_%n_q, "tmp_muq")
      call construct_static_t_muq(target_, tmp_muq)
      matrix_constant(target_%n_q+1:,1:target_%n_q) = tmp_muq

      !wFQ - wFMu Block
      if(.not.target_%heterogeneous) then
         call mem_man%alloc(tmp_qmu, target_%n_q, target_%n_mu, "tmp_qmu")
         call construct_Aqmu_matrix(target_, target_%K_minus_P, &
                                    tmp_qmu, T_muq=tmp_muq)
         matrix_constant(1:target_%n_q,target_%n_q+1:) = tmp_qmu
         call mem_man%dealloc(tmp_qmu, "tmp_qmu")
      endif
      call mem_man%dealloc(tmp_muq, "tmp_muq")

      !Mu - Mu Block
      call mem_man%alloc(tmp_mumu, target_%n_mu, target_%n_mu, "tmp_mumu")
      call construct_static_t_mumu(target_, tmp_mumu)
      matrix_constant(target_%n_q+1:,target_%n_q+1:) = tmp_mumu
      call mem_man%dealloc(tmp_mumu, "tmp_mumu")
  
      if(out_%ivrb.ge.3) call out_%print_matrix('wfqfmu Constant Matrix',   &
                                                matrix_constant,            &
                                                target_%n_var, target_%n_var)
  
   end subroutine construct_constant_matrix_wfqfmu


   !> Subroutine for constructing the wfqfmu RHS on memory
   !!    In/Out : target_         -- wfqfmu
   !!    In/Out : rhs_w           -- dynamic RHS
   subroutine construct_dynamic_field_rhs_wfqfmu(target_,rhs_w)
  
      implicit none
       
      !input/output variables
      class(wfqfmu_type), intent(inout) :: target_ 
      complex(dp), dimension(target_%n_var,3), intent(inout) :: rhs_w

      !internal variables
      complex(dp), dimension(:,:), allocatable :: tmp_q
      complex(dp), dimension(:,:), allocatable :: tmp_mu
  
      target_%n_rhsre = 3
      call mem_man%alloc(target_%rhs_static,target_%n_q,3,"target_%rhs_static")
      !static RHS
      call construct_static_field_rhs_q(target_, &
                                        target_%rhs_static(1:target_%n_atoms,:))
      !wfq block
      if(.not.target_%heterogeneous) then
         call mem_man%alloc(tmp_q, target_%n_q, 3, "tmp_q")
         call construct_rhs_w_q(target_, target_%rhs_static, target_%K_minus_P,&
                                tmp_q)
         rhs_w(1:target_%n_q,:) = tmp_q
         call mem_man%dealloc(tmp_q, "tmp_q")
      else
         !rhs_w = (rhs_0,0) 
         !if heterogeneous depends on frequency. Thus, it is constructed later
         if (field%polarization.eq.'all') then
            rhs_w(1:target_%n_q,:) = &
                               dcmplx(target_%rhs_static(1:target_%n_q,:),zero)
         else
            rhs_w(1:target_%n_q,field%polarization_index) = &
            dcmplx(target_%rhs_static(1:target_%n_q,field%polarization_index),&
                                      zero)
         endif
      endif
      !wfmu block
      call mem_man%alloc(tmp_mu, target_%n_mu, 3, "tmp_mu")
      call construct_rhs_w_mu(target_, tmp_mu)
      rhs_w(target_%n_q + 1:, :) = tmp_mu
      call mem_man%dealloc(tmp_mu, "tmp_mu")
  
      !printing
      if (out_%ivrb.ge.3) call print_dynamic_field_rhs(target_,rhs_w)
      call mem_man%dealloc(target_%rhs_static,"target_%rhs_static")
  
   end subroutine construct_dynamic_field_rhs_wfqfmu


   !> Subroutine for constructing the wfq RHS on the fly
   !!    In/Out : target_         -- wfqfmu
   !!    In/Out : rhs_w           -- dynamic RHS
   subroutine construct_dynamic_field_rhs_on_the_fly_wfqfmu(target_,rhs_w)
  
      implicit none
       
      !input/output variables
      class(wfqfmu_type), intent(inout) :: target_
      complex(dp), dimension(target_%n_var,3), intent(inout) :: rhs_w

      !internal variables
      complex(dp), dimension(:,:), allocatable :: tmp_q
      complex(dp), dimension(:,:), allocatable :: tmp_mu
  
      target_%n_rhsre = 3
      call mem_man%alloc(target_%rhs_static,target_%n_q,3,"target_%rhs_static")
  
      !construct static field RHS
      call construct_static_field_rhs_q(target_, target_%rhs_static)
      !wfq block
      if(.not.target_%heterogeneous) then
         call mem_man%alloc(tmp_q, target_%n_q, 3, "tmp_q")
         call construct_rhs_w_q_on_the_fly(target_, target_%rhs_static, tmp_q)
         rhs_w(1:target_%n_q,:) = tmp_q
         call mem_man%dealloc(tmp_q, "tmp_q")
      else
         !rhs_w = (rhs_0,0) 
         !if heterogeneous depends on frequency. Thus, it is constructed later
         if (field%polarization.eq.'all') then
            rhs_w(1:target_%n_q,:) = &
                               dcmplx(target_%rhs_static(1:target_%n_q,:),zero)
         else
            rhs_w(1:target_%n_q,field%polarization_index) = &
            dcmplx(target_%rhs_static(1:target_%n_q,field%polarization_index),&
                                      zero)
         endif
      endif
      !wfmu block
      call mem_man%alloc(tmp_mu, target_%n_mu, 3, "tmp_mu")
      call construct_rhs_w_mu(target_, tmp_mu)
      rhs_w(target_%n_q + 1:, :) = tmp_mu
      call mem_man%dealloc(tmp_mu, "tmp_mu")
      
      !printing
      if (out_%ivrb.ge.3) call print_dynamic_field_rhs(target_,rhs_w)
      call mem_man%dealloc(target_%rhs_static,"target_%rhs_static")
  
   end subroutine construct_dynamic_field_rhs_on_the_fly_wfqfmu


   !> Subroutine for constructing the wfqfmu dynamic matrix
   !!    In/Out : target_         -- wfqmu
   !!    In/Out : i_freq          -- index of frequency
   !!    Input  : matrix_constant -- constant part of the dynamic matrix
   !!    In/Out : matrix_w        -- dynamic matrix
   subroutine construct_dynamic_matrix_wfqfmu(target_,i_freq,matrix_constant, &
                                              matrix_w)
  
      implicit none
       
      !input/output variables
      class(wfqfmu_type), intent(inout) :: target_
      integer, intent(in)               :: i_freq
      real(dp), dimension(target_%n_var,target_%n_var), intent(in) :: &
                                                                 matrix_constant
      complex(dp), dimension(target_%n_var,target_%n_var), intent(inout) :: &
                                                                        matrix_w
  
      !internal variables
      integer     :: index_1
      integer     :: i,j
      real(dp)    :: freq_au
      complex(dp) :: z_w_q
      complex(dp) :: z_w_mu
  
      freq_au = field%freq(i_freq)
       
      !wFQ Block if heterogeneous
      if(target_%heterogeneous) then
         !1) construct H_minus_L matrix (\bar{H})
         call construct_H_minus_L_matrix(target_, freq_au, target_%H_minus_L)
         !2) H_minus_L = \bar{K} + \bar{H}
         target_%H_minus_L = target_%H_minus_L + target_%K_minus_P 
         !3) construct Aqq block of the matrix
         call construct_Aqq_matrix_complex(target_, half, target_%H_minus_L, &
                                          matrix_w(1:target_%n_q,1:target_%n_q))
         !4)diagonal
         do i = 1, target_%n_atoms 
            !z(w) = - \frac{ w }{(2 * n_0 * tau)} * ( w * tau + i ) }
            z_w_q = -( freq_au / ( two* &
                                  parameters%density(target_%map_atomtypes(i))*&
                                  parameters%tau(target_%map_atomtypes(i)) ) )*&
                      dcmplx( freq_au*parameters%tau(target_%map_atomtypes(i)),&
                              one ) 
            !matrix_w(i,i) = matrix_w(i,i) - z_q(w)
            matrix_w(i,i) = matrix_w(i,i) - z_w_q
         enddo 
         !5) construct the Aqmu block
         call construct_Aqmu_matrix_complex(target_, half, target_%H_minus_L,  &
                                       matrix_w(1:target_%n_q,target_%n_q+1:), &
                          T_muq = matrix_constant(target_%n_q+1:,1:target_%n_q))
      !wFQ Block if homogeneous material
      else 
         do i = 1, target_%n_atoms
            do j = 1, target_%n_atoms
               !matrix_w = constant_matrix
               matrix_w(i,j) = dcmplx(matrix_constant(i,j),zero)
            enddo
            !diagonal
            !z(w) = - \frac{ w }{(2 * n_0 * tau)} * ( w * tau + i ) }
            z_w_q = - ( freq_au / ( two* &
                                  parameters%density(target_%map_atomtypes(i))*&
                                  parameters%tau(target_%map_atomtypes(i)) ) )*&
                      dcmplx( freq_au*parameters%tau(target_%map_atomtypes(i)),&
                              one ) 
            !matrix_w(i,i) = matrix_w(i,i) - z_q(w)
            matrix_w(i,i) = matrix_w(i,i) - z_w_q
         enddo
         !Aqmu block
         matrix_w(1:target_%n_atoms,target_%n_atoms+1:) = &
             dcmplx(matrix_constant(1:target_%n_atoms,target_%n_atoms+1:),zero)
      endif
  
      !wFMu-wFQ blocks
      matrix_w(target_%n_atoms+1:,1:target_%n_atoms) = &
             dcmplx(matrix_constant(target_%n_atoms+1:,1:target_%n_atoms),zero)
  
      !wFMus block
      matrix_w(target_%n_atoms+1:,target_%n_atoms+1:) = &
            dcmplx(matrix_constant(target_%n_atoms+1:,target_%n_atoms+1:),zero)
  
      !diagonal mu-mu block
      do i = 1, target_%n_atoms 
         z_w_mu = dcmplx(zero,zero)
         do j = 1, target_%n_atomtypes
            !this is written for the general case of heterogeneous
            !this reduces to 1/alpha(w) for homogeneous
            z_w_mu = z_w_mu - (real(target_%neighbours(j,i),kind=dp)/ &
                               real(target_%n_neighbours(i),kind=dp))/&
                               parameters%alpha_w(j,i_freq)
         enddo
         index_1 = 3*(i-1) + target_%n_atoms
         matrix_w(index_1+1,index_1+1) = - z_w_mu
         matrix_w(index_1+2,index_1+2) = - z_w_mu
         matrix_w(index_1+3,index_1+3) = - z_w_mu
      enddo
  
      If(out_%ivrb.ge.3) then
         call out_%print_matrix("wfqfmu Matrix: Real Part",dble(matrix_w), &
                                                   target_%n_var,target_%n_var)
         call out_%print_matrix("wfqfmu Matrix: IMag Part",dimag(matrix_w),&
                                                   target_%n_var,target_%n_var)
      endif
  
   end subroutine construct_dynamic_matrix_wfqfmu


   !> Subroutine for constructing dynamic matrix [GMRES algorithm]
   !! This is only needed if the constant matrix depends on the frequency
   !! i.e. right now only for heterogeneous structures
   !!    In/Out  : target_          -- model 
   !!    Input   : i_freq           -- index_ of frequencies
   !!    In/Out  : matrix_iterative -- dynamic matrix iterative
   subroutine construct_dynamic_matrix_gmres_wfqfmu(target_,i_freq, &
                                                    matrix_iterative)
  
      implicit none
       
      !input/output variables
      class(wfqfmu_type), intent(inout) :: target_
      integer, intent(in)               :: i_freq
      complex(dp), dimension(target_%n_var,target_%n_var), intent(inout) :: &
                                                                 matrix_iterative
  
      !internal variables
      real(dp)          :: freq_au
  
      freq_au = field%freq(i_freq)
      !1) construct H-L (\bar{H})
      call construct_H_minus_L_matrix(target_, freq_au, target_%H_minus_L)
      !2) H_minus_L = \bar{K} + \bar{H}
      target_%H_minus_L = target_%H_minus_L + target_%K_minus_P 
      !3) construct A_qq block
      call construct_Aqq_matrix_complex(target_, half, target_%H_minus_L, &
                                 matrix_iterative(1:target_%n_q, 1:target_%n_q))
      !3) construct A_qmu block
      call construct_Aqmu_matrix_complex(target_, half, target_%H_minus_L, &
                               matrix_iterative(1:target_%n_q,target_%n_q+1:), &
                   T_muq = dble(matrix_iterative(target_%n_q+1:,1:target_%n_q)))
  
   end subroutine construct_dynamic_matrix_gmres_wfqfmu


   !> Subroutine for performing the matrix vector multiplication on the fly
   !! See comment on wfq_module file
   !!
   !!    In/Out  : target_   -- model 
   !!    Input   : i_freq    -- index of frequency
   !!    In/Out  : x         -- input vector 
   !!    In/Out  : y         -- output vector 
   subroutine product_matrix_vector_wfqfmu(target_,i_freq, x, y)
  
      implicit none
       
      !input/output variables
      class(wfqfmu_type), intent(inout) :: target_
      integer, intent(in)               :: i_freq
      complex(dp), dimension(target_%n_var), intent(inout) :: x
      complex(dp), dimension(target_%n_var), intent(inout) :: y
  
      !internal variables
      integer :: i
      real(dp) :: freq_au
      complex(dp), dimension(:), allocatable :: static_matrix_dot_x
      complex(dp), dimension(:), allocatable :: q_block_x
       
      call mem_man%alloc(static_matrix_dot_x,target_%n_var, &
                         "static_matrix_dot_x")
      !$omp parallel do
      do i=1,target_%n_var
         y(i) = dcmplx(zero,zero)
      enddo
      !$omp end parallel do
  
      freq_au = field%freq(i_freq)
      !product q block * x
      call mem_man%alloc(q_block_x,target_%n_q, "q_block_x")
      !call product_q_block_x(target_, x(1:target_%n_q), &
      !                       static_matrix_dot_x(1:target_%n_q))
      call product_q_block_x(target_, x(1:target_%n_q), q_block_x)
      static_matrix_dot_x(1:target_%n_q) = q_block_x
      call mem_man%dealloc(q_block_x, "q_block_x")
      
      !product mu block * x
      call product_mu_block_x(target_, x, static_matrix_dot_x) 
      !update y by q block
      call update_y_q(target_, freq_au, static_matrix_dot_x, y) 
      !update y by mu block
      call update_y_mu(target_, static_matrix_dot_x, y)  
  
      call mem_man%dealloc(static_matrix_dot_x, "static_matrix_dot_x")
       
   end subroutine product_matrix_vector_wfqfmu


   !> Subroutine for applying the diagonal shift for GMRES algorithm
   !!    In/Out  : target_      -- model 
   !!    Input   : i_freq       -- index of the frequency
   !!    Input   : colx         -- index of the column x (GMRES)
   !!    Input   : colz         -- index of the column z (GMRES)
   !!    Input   : lwork        -- dimension of the work array
   !!    In/Out  : work         -- array work
   subroutine gmres_diagonal_shift_wfqfmu(target_,i_freq,colx,colz,lwork,work)
  
      implicit none
       
      !input/output variables
      class(wfqfmu_type)    :: target_
      integer  :: i_freq
      integer  :: colx
      integer  :: colz
      integer  :: lwork
      complex(dp), dimension(lwork), intent(inout) :: work
  
      !internal variables
      integer  :: i, j
      integer  :: index_1
      integer  :: start_x
      integer  :: start_z
      real(dp) :: freq_au
      complex(dp) :: z_i
  
      freq_au = field%freq(i_freq)
  
      !diagonal part wfq 
      start_z = colz - 1
      start_x = colx - 1
      do i = 1, target_%n_atoms
         z_i = freq_au * &
               (dcmplx(-freq_au*parameters%tau(target_%map_atomtypes(i)), &
                       -one)) / & 
               (two*parameters%density(target_%map_atomtypes(i)) * &
                    parameters%tau(target_%map_atomtypes(i)))
         work(start_z+i) = work(start_z + i) - z_i * work(start_x + i)
      enddo
  
      !diagonal part wfqfmu
      start_z = colz + target_%n_atoms - 1
      start_x = colx + target_%n_atoms - 1
      do i = 1, target_%n_atoms
         z_i = dcmplx(zero,zero)
         do j = 1, target_%n_atomtypes
            !this is written for the general case of alloys
            z_i = z_i + (real(target_%neighbours(j,i),kind=dp)/ &
                         real(target_%n_neighbours(i),kind=dp))/&
                         parameters%alpha_w(j,i_freq)
         enddo
         index_1 = 3*(i-1)
         work(start_z + index_1 + 1) = work(start_z + index_1 + 1) + & 
                                       z_i * work(start_x + index_1 + 1)
         work(start_z + index_1 + 2) = work(start_z + index_1 + 2) + &
                                       z_i * work(start_x + index_1 + 2)
         work(start_z + index_1 + 3) = work(start_z + index_1 + 3) + &
                                       z_i * work(start_x + index_1 + 3)
      enddo
       
   end subroutine gmres_diagonal_shift_wfqfmu


   !> Subroutine for applying the diagonal shift for GMRES algorithm
   !!    In/Out  : target_      -- model 
   !!    Input   : i_freq       -- index of the frequency
   !!    Input   : r            -- rhs
   !!    In/Out  : w            -- solution
   subroutine new_gmres_diagonal_shift_wfqfmu(target_,i_freq,r,w)
  
      implicit none
       
      !input/output variables
      class(wfqfmu_type)    :: target_
      integer  :: i_freq
      complex(dp), dimension(target_%n_var), intent(in) :: r
      complex(dp), dimension(target_%n_var), intent(inout) :: w
  
      !internal variables
      integer  :: i, j
      integer  :: index_1
      integer  :: start
      real(dp) :: freq_au
      complex(dp) :: z_i
  
      freq_au = field%freq(i_freq)
  
      !diagonal part wfq 
      do i = 1, target_%n_atoms
         z_i = freq_au * &
               (dcmplx(-freq_au*parameters%tau(target_%map_atomtypes(i)), &
                       -one)) / & 
               (two*parameters%density(target_%map_atomtypes(i)) * &
                    parameters%tau(target_%map_atomtypes(i)))
         w(i) = w(i) - z_i * r(i)
      enddo
  
      !diagonal part wfqfmu
      start = target_%n_atoms 
      do i = 1, target_%n_atoms
         z_i = dcmplx(zero,zero)
         do j = 1, target_%n_atomtypes
            !this is written for the general case of alloys
            z_i = z_i + (real(target_%neighbours(j,i),kind=dp)/ &
                         real(target_%n_neighbours(i),kind=dp))/&
                         parameters%alpha_w(j,i_freq)
         enddo
         index_1 = 3*(i-1)
         w(start + index_1 + 1) = w(start + index_1 + 1) + & 
                                  z_i * r(start + index_1 + 1)
         w(start + index_1 + 2) = w(start + index_1 + 2) + &
                                  z_i * r(start + index_1 + 2)
         w(start + index_1 + 3) = w(start + index_1 + 3) + &
                                  z_i * r(start + index_1 + 3)
      enddo
       
   end subroutine new_gmres_diagonal_shift_wfqfmu

!-------------------------------------------------------------------------------
!  SUBROUTINES FOR PRINTING & SAVING FILES
!-------------------------------------------------------------------------------

   !> Subroutine for printing the variables when dynamic field is applied
   !!    Input   : target_      -- model 
   !!    Input   : variables_w  -- w-variables (N_var, 3 [x,y,z]) 
   subroutine print_dynamic_field_variables_wfqfmu(target_,variables_w)
  
      implicit none
       
      !input/output variables
      class(wfqfmu_type), intent(in) :: target_
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
      character(len=22) :: format_5 = "(38x,'wFMu Dipoles',/)"
  
      !wFQs
      write(out_%iunit,format_1)
      write(out_%iunit,format_2)
      do i = 1, target_%n_atoms
         write(out_%iunit,format_3) i,variables_w(i,:)
         sum_x = sum_x + variables_w(i,1)
         sum_y = sum_y + variables_w(i,2)
         sum_z = sum_z + variables_w(i,3)
      enddo
      write(out_%iunit,format_4) sum_x, sum_y, sum_z
       
      !wFMus
      write(out_%iunit,out_%sticks) 
      write(out_%iunit,format_5)
      do i = target_%n_atoms+1, target_%n_var
         write(out_%iunit,format_3) i - target_%n_atoms, variables_w(i,:)
      enddo
      write(out_%iunit,out_%sticks) 
       
   end subroutine print_dynamic_field_variables_wfqfmu


   !> Subroutine for printing the atomtypes
   !!    Input  : target_   -- model 
   subroutine print_atomtypes_wfqfmu(target_)
       
      implicit none
  
      !input/output variables
      class(wfqfmu_type), intent(in) :: target_
  
      !internal variables
      integer            :: i
      integer            :: j
       
      do i = 1, target_%n_atomtypes
         write(out_%iunit,out_%sticks) 
         write(out_%iunit,'(1x,a)') "AtomType           : "// &
                                    trim(target_%atom_type(i))
         write(out_%iunit,'(1x,a,f10.3,a)') "Chi                : ", &
                                            target_%chi(i), " a.u."
         write(out_%iunit,'(1x,a,f10.3,a)') "Eta                : ", &
                                            target_%eta(i), " a.u."
         write(out_%iunit,'(1x,a,f10.3,a)') "Alpha              : ", &
                                            target_%alpha(i)," a.u."
         write(out_%iunit,'(1x,a,f10.3,a)') "R_q                : ", &
                                            target_%r_q(i), " a.u."
         write(out_%iunit,'(1x,a,f10.3,a)') "R_mu               : ", &
                                            target_%r_mu(i), " a.u."
      
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
         if(allocated(parameters%density)) &
            write(out_%iunit,'(1x,a,e10.3,a)')   "density            : ", &
                                               parameters%density(i), " a.u."
         if(allocated(parameters%permittivity_type)) & 
            write(out_%iunit,'(1x,a)') "wFQFMu IB alpha    :  "// &
                                        trim(parameters%permittivity_type(i))
         if(allocated(parameters%wfqfmu_file)) & 
            write(out_%iunit,'(1x,a)') "wFQFMu file        :  "// &
                                        trim(parameters%wfqfmu_file(i))
      enddo
      if (target_%heterogeneous) then
         do i = 1, target_%n_atomtypes
           do j = 1, target_%n_atomtypes
             if (i.ne.j) then
               write(out_%iunit,out_%sticks) 
               write(out_%iunit,'(1x,a)') "Interaction        : "// &
                   trim(target_%atom_type(i))//"->"//trim(target_%atom_type(j))
               write(out_%iunit,'(1x,a,e10.3,a)') "fermi_d            : ", &
                                                   parameters%fermi_d(i,j),&
                                                   " a.u."
               write(out_%iunit,'(1x,a,e10.3,a)') "fermi_s            : ", &
                                                   parameters%fermi_s(i,j),&
                                                   " a.u."
             endif
           enddo
         enddo
      end if

   end subroutine print_atomtypes_wfqfmu

!-------------------------------------------------------------------------------
!  SUBROUTINES FOR CALCULATING PROPERTIES
!-------------------------------------------------------------------------------

   !> Subroutine for calculating the dynamic polar
   !!    In/Out  : target_      -- model 
   !!    Input   : i_freq       -- index of the frequency
   !!    Input   : variables_w  -- w-variables
   !!    Output  : polar_w      -- dynamic polar
   subroutine calculate_dynamic_polar_wfqfmu(target_,i_freq,variables_w,polar_w)
  
      implicit none
       
      !input/output variables
      class(wfqfmu_type), intent(inout) :: target_
      integer, intent(in)               :: i_freq
      complex(dp), dimension(target_%n_var,3), intent(in) :: variables_w
      complex(dp), dimension(3,3), intent(inout)          :: polar_w
  
      !internal variables
      integer :: i
      integer :: index_1
      complex(dp), dimension(3,3) :: internal_polar
  
      internal_polar = dcmplx(zero,zero)
      !wfq contribution
      do i = 1, target_%n_atoms
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
       
      !wfmu contribution
      do i = 1, target_%n_atoms
        index_1 = 3*(I-1) + target_%n_atoms
        internal_polar(1,1) = internal_polar(1,1) + &
                              variables_w(index_1+1,1)/field%e_0
        internal_polar(1,2) = internal_polar(1,2) + &
                              variables_w(index_1+1,2)/field%e_0
        internal_polar(1,3) = internal_polar(1,3) + &
                              variables_w(index_1+1,3)/field%e_0
        internal_polar(2,1) = internal_polar(2,1) + &
                              variables_w(index_1+2,1)/field%e_0
        internal_polar(2,2) = internal_polar(2,2) + &
                              variables_w(index_1+2,2)/field%e_0
        internal_polar(2,3) = internal_polar(2,3) + &
                              variables_w(index_1+2,3)/field%e_0
        internal_polar(3,1) = internal_polar(3,1) + &
                              variables_w(index_1+3,1)/field%e_0
        internal_polar(3,2) = internal_polar(3,2) + &
                              variables_w(index_1+3,2)/field%e_0
        internal_polar(3,3) = internal_polar(3,3) + &
                              variables_w(index_1+3,3)/field%e_0
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
  
   end subroutine calculate_dynamic_polar_wfqfmu


   !> Function for calculating the induced field at a specific point
   !!    In/Out  : target_      -- model 
   !!    Input   : i_pol        -- index of the polarization
   !!    Input   : point_coord  -- coordinates of the point
   !!    Input   : variables_w  -- w-variables
   !!    Output  : EField       -- Electric Field
   function calculate_induced_field_at_point_wfqfmu(target_,i_pol,point_coord, &
                                                    variables_w) result(EField)
  
      implicit none
       
      !input/output variables
      class(wfqfmu_type), intent(in)     :: target_
      integer, intent(in) :: i_pol
      real(dp), dimension(3), intent(in) :: point_coord
      complex(dp), dimension(target_%n_var), intent(in) :: variables_w
      complex(dp), dimension(3) :: EField
  
      !internal variables
      integer  :: i
      integer  :: index_mu
      real(dp), dimension(3) :: d_IJ
       
      EField = dcmplx(zero, zero)
      EField(i_pol) = dcmplx(field%e_0, zero)
      do i = 1, target_%n_atoms 
         d_IJ(1)   = (target_%coord(1,i)-point_coord(1))*tobohr
         d_IJ(2)   = (target_%coord(2,i)-point_coord(2))*tobohr
         d_IJ(3)   = (target_%coord(3,i)-point_coord(3))*tobohr
         !wfq contribution
         call calculate_induced_field_q(target_,i,d_IJ,variables_w(i),EField)
         !wfmu contribution
         index_mu = 3*(i-1) + target_%n_atoms + 1
         call calculate_induced_field_mu(target_,i,d_IJ,variables_w(index_mu), &
                                         EField)
      enddo                  
 
   end function calculate_induced_field_at_point_wfqfmu


   !> Function for calculating the plasmon density at a specific point
   !!    Input   : target_      -- model 
   !!    Input   : point_coord  -- coordinates of the point
   !!    Input   : variables_w  -- w-variables
   !!    Output  : density      -- Plasmon density
   function calculate_density_at_point_wfqfmu(target_,point_coord, &
                                                    variables_w) result(density)
  
      implicit none
       
      !input/output variables
      class(wfqfmu_type), intent(in)     :: target_
      real(dp), dimension(3), intent(in) :: point_coord
      complex(dp), dimension(target_%n_var), intent(in) :: variables_w
      complex(dp) :: density
  
      !internal variables
      integer  :: i
      integer  :: index_mu
      real(dp), dimension(3) :: d_IJ
       
      density = zero
      do i = 1, target_%n_atoms  
         d_IJ(1)   = (target_%coord(1,i)-point_coord(1))*tobohr
         d_IJ(2)   = (target_%coord(2,i)-point_coord(2))*tobohr
         d_IJ(3)   = (target_%coord(3,i)-point_coord(3))*tobohr
         !wfq contribution
         call calculate_density_q(target_,i,d_IJ,variables_w(i),density)
         !wfmu contribution
         index_mu = 3*(i-1) + target_%n_atoms + 1
         call calculate_density_mu(target_,i,d_IJ,variables_w(index_mu),density)
      enddo                  
       
   end function calculate_density_at_point_wfqfmu

end module wfqfmu_module
