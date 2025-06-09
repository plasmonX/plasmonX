!> Matrix module
!!
!! This module contains the subroutines for calculating quantities
!! common to more than one model.
!!
!! Date         : 2024 
!!
module matrix_module

   use output_module
   use parameters_module
   use string_manipulation_module
   use array_manipulation_module
   use field_module
   use target_module
   use memory_manager_module

   Implicit None

   !Interaction kernel
   public construct_static_t_qq            ! q-q   kernel
   public construct_static_t_qmu           ! q-mu  kernel
   public construct_static_t_muq           ! mu-q  kernel 
   public construct_static_t_mumu          ! mu-mu kernel
   public add_lagrangian_blocks            ! 1 blocks for Lagrangian

   !Matrices for w-models
   public construct_K_minus_P_matrix       ! \bar{K} matrix
   public construct_H_minus_L_matrix       ! \bar{H} matrix
   public construct_Aqq_matrix             ! A_qq  matrix
   public construct_Aqmu_matrix            ! A_qmu matrix
   public construct_Aqq_matrix_complex     ! A_qq(w)  matrix for hetero
   public construct_Aqmu_matrix_complex    ! A_qmu(w) matrix for hetero 
   public fermi_function                   ! Fermi function

   !GMRES-related operations (on the fly)
   public product_q_block_x                ! q-block  dot x (GMRES otf)
   public product_mu_block_x               ! mu-block dot x (GMRES otf)
   public update_y_q                       ! update  q-block y (GMRES otf)
   public update_y_mu                      ! update mu-block y (GMRES otf)

   !Analysis (post-process)
   public calculate_induced_field_q        ! q-induced field (analysis)
   public calculate_induced_field_mu       ! mu-induced field (analysis)
   public calculate_field_gradient_num_q   ! q-field gradient (analysis)
   public calculate_field_gradient_num_mu  ! mu-field gradient (analysis)
   public calculate_induced_field_q_bem    ! q-induced field bem (analysis)

contains

!-------------------------------------------------------------------------------
!  INTERACTION KERNEL SUBROUTINES
!-------------------------------------------------------------------------------

   !> Subroutine for constructing Tqq block of the interaction kernel
   !!    Input  : target_      -- model 
   !!    Output : matrix_block -- charge-charge block
   subroutine construct_static_t_qq(target_, matrix_block)

      implicit none

      !input/output variables
      class(target_type), intent(in)                            :: target_
      real(dp), dimension(target_%n_q,target_%n_q), intent(out) :: matrix_block

      !internal variables
      integer :: i, j
      real(dp) :: distIJ

      !$omp parallel do collapse(2) private(j,i,distij)
      do j = 1, target_%n_q
        do i = 1, target_%n_q
          distIJ = dsqrt((target_%coord(1,i)-target_%coord(1,j))**2 + &
                         (target_%coord(2,i)-target_%coord(2,j))**2 + &
                         (target_%coord(3,i)-target_%coord(3,j))**2 )
          ! Coulomb kernel: 
          ! 1/r_ij (off diagonal)
          ! eta_i  (diagonal)
          If (target_%kernel.eq.'coulomb') then !coulomb
             if(i.ne.j) then
                matrix_block(i,j) = one/distij
             else
                matrix_block(i,j) = (target_%eta(target_%map_atomtypes(i)) + &
                                     target_%eta(target_%map_atomtypes(j)))/two
             endIf
          ! Ohno kernel (see JCTC 2011, XX, ZZZ)
          ! 1/r_ij (off diagonal)
          ! eta_i  (diagonal)
          else If(target_%kernel.eq.'ohno') then !ohno
             matrix_block(i,j) = half*(target_%eta(target_%map_atomtypes(i)) + &
                                       target_%eta(target_%map_atomtypes(j)))/ & 
                        (one + ((half*(target_%eta(target_%map_atomtypes(i)) + &
                               target_%eta(target_%map_atomtypes(j))))**two) * &
                               (distij**two))**(half)
          ! Gaussian kernel (see JCTC 2011, XX, ZZZ)
          ! 1/r_ij (off diagonal)
          ! eta_i  (diagonal)
          else if (target_%kernel.eq.'gaussian') then !gaus
             if(i.ne.j) then
                matrix_block(i,j) = erf(distIJ/ &
                         dsqrt(target_%r_q(target_%map_atomtypes(i))**2 + &
                               target_%r_q(target_%map_atomtypes(j))**2))/distIJ
             else
                matrix_block(i,j) = (target_%eta(target_%map_atomtypes(i)) + &
                                     target_%eta(target_%map_atomtypes(j)))/two
             endif
          endif
        enddo
      enddo
      !$omp end parallel do

      if(out_%ivrb.ge.4) &
         call out_%print_matrix('Tqq block',matrix_block,target_%n_q,target_%n_q)
     
   end subroutine construct_static_t_qq

   !> Subroutine for constructin Tqmu block of the interaction kernel
   !!    Input  : target_      -- model 
   !!    Output : matrix_block -- charge-dipole block
   subroutine construct_static_t_qmu(target_, matrix_block)
     implicit none
     
     !input/output variables
     class(target_type), intent(in)                             :: target_
     real(dp), dimension(target_%n_q,target_%n_mu), intent(out) :: matrix_block

     !internal variables
     integer :: i, j
     integer :: index_1
     real(dp) :: distIJ
     real(dp) :: distIJ_2
     real(dp) :: distIJ_3
     real(dp) :: f_erf_ij
     real(dp) :: f_erf_ij_2
     real(dp) :: x_ij, y_ij, z_ij

     do i = 1, target_%n_atoms
        do j = 1 , target_%n_atoms
           if(i.ne.j) then ! charge and dipole on the same atom do not interact
              x_ij   = target_%coord(1,i)-target_%coord(1,j)
              y_ij   = target_%coord(2,i)-target_%coord(2,j)
              z_ij   = target_%coord(3,i)-target_%coord(3,j)

              distij   = dsqrt(x_ij**2 + y_ij**2 + z_ij**2)
              distij_2 = distij**two
              distij_3 = distij**three

              f_erf_IJ   = dsqrt(target_%r_q(target_%map_atomtypes(i))**2 + &
                                 target_%r_mu(target_%map_atomtypes(j))**2)
              f_erf_IJ_2 = f_erf_IJ**two

              index_1 = 3*(J-1) 

              matrix_block(i,index_1+1) = x_ij/distij_3 *                      & 
                                          (derf(distij/f_erf_ij) -             &
                                          ((two*distij)/(dsqrt(pi)*f_erf_ij))* &
                                            dexp(-distij_2/f_erf_ij_2))
              matrix_block(i,index_1+2) = y_ij/distij_3 *                      &
                                          (derf(distij/f_erf_ij) -             &
                                          ((two*distij)/(dsqrt(pi)*f_erf_ij))* &
                                            dexp(-distij_2/f_erf_ij_2))
              matrix_block(i,index_1+3) = z_ij/distij_3 *                      &
                                          (derf(distij/f_erf_ij) -             &
                                          ((two*distij)/(dsqrt(pi)*f_erf_ij))* &
                                            dexp(-distij_2/f_erf_ij_2))
           endif
        enddo
     enddo

     if(out_%ivrb.ge.4) &
        call out_%print_matrix('Tqmu block',matrix_block,target_%n_q,target_%n_mu)

   end subroutine construct_static_t_qmu
   

   !> Subroutine for constructin Tmuq block of the interaction kernel
   !!    Input  : target_      -- model 
   !!    Output : matrix_block -- dipole-charge block
   subroutine construct_static_t_muq(target_, matrix_block)

      implicit none
      
      !input/output variables
      class(target_type), intent(in)                             :: target_
      real(dp), dimension(target_%n_mu,target_%n_q), intent(out) :: matrix_block

      !internal variables
      integer :: i, j
      integer :: index_1
      real(dp) :: distIJ
      real(dp) :: distIJ_2
      real(dp) :: distIJ_3
      real(dp) :: f_erf_ij
      real(dp) :: f_erf_ij_2
      real(dp) :: x_ij, y_ij, z_ij

      do i = 1, target_%n_atoms
         do j = 1 , target_%n_atoms
            if(i.ne.j) then ! exclude dipole-charge interaction in the same atom

               x_ij   = target_%coord(1,i)-target_%coord(1,j)
               y_ij   = target_%coord(2,i)-target_%coord(2,j)
               z_ij   = target_%coord(3,i)-target_%coord(3,j)

               distij   = dsqrt(x_ij**2 + y_ij**2 + z_ij**2)
               distij_2 = distij**two
               distij_3 = distij**three

               f_erf_IJ   = dsqrt(target_%r_q(target_%map_atomtypes(i))**2 + &
                                  target_%r_mu(target_%map_atomtypes(j))**2)
               f_erf_IJ_2 = f_erf_IJ**two

               index_1 = 3*(J-1) 

               matrix_block(index_1+1,i) = x_ij/distij_3 *                      & 
                                           (derf(distij/f_erf_ij) -             &
                                           ((two*distij)/(dsqrt(pi)*f_erf_ij))* &
                                             dexp(-distij_2/f_erf_ij_2))
               matrix_block(index_1+2,i) = y_ij/distij_3 *                      &
                                           (derf(distij/f_erf_ij) -             &
                                           ((two*distij)/(dsqrt(pi)*f_erf_ij))* &
                                             dexp(-distij_2/f_erf_ij_2))
               matrix_block(index_1+3,i) = z_ij/distij_3 *                      &
                                           (derf(distij/f_erf_ij) -             &
                                           ((two*distij)/(dsqrt(pi)*f_erf_ij))* &
                                             dexp(-distij_2/f_erf_ij_2))
            endif
         enddo
      enddo

      if(out_%ivrb.ge.4) &
         call out_%print_matrix('Tmuq block',matrix_block,target_%n_mu,target_%n_q)

   end subroutine construct_static_t_muq


   !> Subroutine for constructin Tmuq block of the interaction kernel
   !!    Input  : target_      -- model 
   !!    Output : matrix_block -- dipole-charge block
   subroutine construct_static_t_mumu(target_, matrix_block)

      implicit none

      !input/output variables
      class(target_type), intent(in)                              :: target_
      real(dp), dimension(target_%n_mu,target_%n_mu), intent(out) :: matrix_block

      !internal variables
      integer :: i, j
      integer :: index_1
      integer :: index_2
      real(dp) :: distIJ
      real(dp) :: distIJ_2
      real(dp) :: distIJ_3
      real(dp) :: distIJ_5
      real(dp) :: f_erf_ij
      real(dp) :: f_erf_ij_2
      real(dp) :: f_erf_ij_3
      real(dp) :: fD1
      real(dp) :: fD2
      real(dp) :: x_ij, y_ij, z_ij
      real(dp) :: Txx, Txy, Txz
      real(dp) :: Tyy, Tyx, Tyz
      real(dp) :: Tzz, Tzx, Tzy

      do i = 1, target_%n_atoms
         do j = 1, i
            x_ij   = target_%coord(1,i)-target_%coord(1,j)
            y_ij   = target_%coord(2,i)-target_%coord(2,j)
            z_ij   = target_%coord(3,i)-target_%coord(3,j)

            distij = dsqrt(x_ij**2 + y_ij**2 + z_ij**2)
            distij_2 = distij**two
            distij_3 = distij**three
            distij_5 = distij**five

            index_1 = 3*(i-1) 
            index_2 = 3*(j-1) 

            if(I.eq.J) then ! diagonal = 1/alpha if .not. w-models
               if(target_%name_(1:3) .ne. 'wfq') then 
                  matrix_block(index_1+1,index_2+1) = one/ &
                                         target_%alpha(target_%map_atomtypes(i))
                  matrix_block(index_1+2,index_2+2) = one/ &
                                         target_%alpha(target_%map_atomtypes(i))
                  matrix_block(index_1+3,index_2+3) = one/ &
                                         target_%alpha(target_%map_atomtypes(i))
               endif
            else
               f_erf_IJ = dsqrt(target_%r_mu(target_%map_atomtypes(i))**2 + &
                                target_%r_mu(target_%map_atomtypes(j))**2)
               f_erf_IJ_2 = f_erf_IJ**two
               f_erf_IJ_3 = f_erf_IJ**three

               fD1 = derf(distIJ/f_erf_IJ) - ( (two*distIJ) /             &
                                               (dsqrt(pi)*f_erf_IJ) ) *   &
                                                dexp(-distIJ_2/f_erf_IJ_2)

               fD2 = (four/(sqrt(pi)*f_erf_IJ_3)) * dexp(-distIJ_2/f_erf_IJ_2)

               Txx = ((-three*x_IJ*x_IJ)/distIJ_5 + one/distij_3)*fD1 + &
                      fD2*(x_IJ*x_IJ)/distIJ_2
               Tyx = ((-three*y_IJ*x_IJ)/distIJ_5               )*fD1 + &
                      fD2*(y_IJ*x_IJ)/distIJ_2
               Tzx = ((-three*z_IJ*x_IJ)/distIJ_5               )*fD1 + &
                      fD2*(z_IJ*x_IJ)/distIJ_2
               Txy = ((-three*x_IJ*y_IJ)/distIJ_5               )*fD1 + &
                      fD2*(x_IJ*y_IJ)/distIJ_2
               Tyy = ((-three*y_IJ*y_IJ)/distIJ_5 + one/distIJ_3)*fD1 + &
                      fD2*(y_IJ*y_IJ)/distIJ_2
               Tzy = ((-three*z_IJ*y_IJ)/distIJ_5               )*fD1 + &
                      fD2*(z_IJ*y_IJ)/distIJ_2
               Txz = ((-three*x_IJ*z_IJ)/distIJ_5               )*fD1 + &
                      fD2*(x_IJ*z_IJ)/distIJ_2
               Tyz = ((-three*y_IJ*z_IJ)/distIJ_5               )*fD1 + &
                      fD2*(y_IJ*z_IJ)/distIJ_2
               Tzz = ((-three*z_IJ*z_IJ)/distIJ_5 + one/distIJ_3)*fD1 + &
                      fD2*(z_IJ*z_IJ)/distIJ_2


               matrix_block(index_1+1,index_2+1) = Txx
               matrix_block(index_1+2,index_2+1) = Tyx
               matrix_block(index_1+3,index_2+1) = Tzx
               matrix_block(index_1+1,index_2+2) = Txy
               matrix_block(index_1+2,index_2+2) = Tyy
               matrix_block(index_1+3,index_2+2) = Tzy
               matrix_block(index_1+1,index_2+3) = Txz
               matrix_block(index_1+2,index_2+3) = Tyz
               matrix_block(index_1+3,index_2+3) = Tzz
                                                                         
               matrix_block(index_2+1,index_1+1) = matrix_block(index_1+1,index_2+1)
               matrix_block(index_2+2,index_1+1) = matrix_block(index_1+2,index_2+1)
               matrix_block(index_2+3,index_1+1) = matrix_block(index_1+3,index_2+1)
               matrix_block(index_2+1,index_1+2) = matrix_block(index_1+1,index_2+2)
               matrix_block(index_2+2,index_1+2) = matrix_block(index_1+2,index_2+2)
               matrix_block(index_2+3,index_1+2) = matrix_block(index_1+3,index_2+2)
               matrix_block(index_2+1,index_1+3) = matrix_block(index_1+1,index_2+3)
               matrix_block(index_2+2,index_1+3) = matrix_block(index_1+2,index_2+3)
               matrix_block(index_2+3,index_1+3) = matrix_block(index_1+3,index_2+3)
            endif
         enddo
      enddo
      
!     
   end subroutine construct_static_t_mumu


   !> Subroutine for adding the Lagrangian blocks to conserve the cahrge
   !!    Input  : target_      -- model
   !!    Output : matrix       -- This the total FQ matrix (n_q + n_mol)
   subroutine add_lagrangian_blocks(target_, matrix)

     implicit none
     
     !input/output variables
     class(target_type), intent(in)   :: target_
     real(dp), dimension(target_%n_q+target_%n_mol, &
                         target_%n_q+target_%n_mol) :: matrix
 
     !internal variables
     integer :: i, j, j2, k

     j2 = 0
     do k = 1, target_%n_mol
        j2 = j2 + 1
        j  = j2 + target_%n_q
        If(k.eq.1) then
           do i = 1, target_%n_atoms_per_molecule(k)
              matrix(i,j) = one
              matrix(j,i) = one
           enddo
        else
           do i = (target_%n_atoms_per_molecule(k-1)*(k-1)) + 1, &
                   target_%n_atoms_per_molecule(k-1)*(k-1) +   &
                   target_%n_atoms_per_molecule(k)
              matrix(i,j) = one
              matrix(j,i) = one
           enddo
        endif
     enddo

   end subroutine add_lagrangian_blocks

!-------------------------------------------------------------------------------
!  MATRICES FOR w-MODELS
!-------------------------------------------------------------------------------

   !> Subroutine for constructing K_bar matrix (K - P) 
   !!    See J. Phys. Chem. C, 2021, 125, 23848-23863
   !! 
   !!    Input  : target_      -- model 
   !!    Output : K-P matrix 
   subroutine construct_K_minus_P_matrix(target_,K_minus_P)

      implicit none
 
      !input/output variables
      class(target_type), intent(in)                            :: target_
      real(dp), dimension(target_%n_q,target_%n_q), intent(out) :: K_minus_P

      !internal variables
      integer  :: i,j 
      real(dp) :: distij

      !$omp parallel do private(j,i,distij)
      do i = 1, target_%n_atoms
         do j = i + 1, target_%n_atoms
            distij    = sqrt((target_%coord(1,i)-target_%coord(1,j))**2 + &
                             (target_%coord(2,i)-target_%coord(2,j))**2 + &
                             (target_%coord(3,i)-target_%coord(3,j))**2 )
            K_minus_P(i,j) = (one-fermi_function(i,j,distij))* & 
                              parameters%a_ij(target_%map_atomtypes(i))/distij
            K_minus_P(j,i) = (one-fermi_function(j,i,distij))* & 
                              parameters%a_ij(target_%map_atomtypes(j))/distij
         enddo
      enddo
      !$omp end parallel do

      !The diagonal contains - the sum of the row
      !$omp parallel do private(i)
      do i=1,target_%n_atoms
         K_minus_P(i,i) = - sum(K_minus_P(i,:))
      enddo
      !$omp end parallel do

   end subroutine construct_K_minus_P_matrix


   !> Subroutine for constructing H_bar matrix (H - L) 
   !!    See eq. 16 - Front. Photon. 2023, 4, 1199598
   !! 
   !!    Input  : target_      -- model 
   !!    Input  : freq         -- frequency
   !!    Output : H-L matrix 
   subroutine construct_H_minus_L_matrix(target_, freq, H_minus_L)

      implicit none
      
      !input/output variables
      class(target_type), intent(in)                               :: target_
      real(dp), intent(in)                                         :: freq
      complex(dp), dimension(target_%n_q,target_%n_q), intent(out) :: H_minus_L

      !internal variables
      integer  :: i,j 
      real(dp) :: distij
      complex(dp) :: w_i
      complex(dp) :: w_j

      do i = 1, target_%n_atoms
         ! w_i = 2 * n_i * tau_i / (1 - i w tau_i)
         w_i = two * parameters%density(target_%map_atomtypes(i)) /   &
                     dcmplx(one/parameters%tau(target_%map_atomtypes(i)), -freq)
         do j = i + 1, target_%n_atoms
            ! w_j = 2 * n_j * tau_j / (1 - i w tau_j)
            w_j = two * parameters%density(target_%map_atomtypes(j)) /   &
                     dcmplx(one/parameters%tau(target_%map_atomtypes(j)), -freq)

            distij    = sqrt((target_%coord(1,i)-target_%coord(1,j))**2 + &
                             (target_%coord(2,i)-target_%coord(2,j))**2 + &
                             (target_%coord(3,i)-target_%coord(3,j))**2 )
            ! i - j elements : K_ji * w_j/w_i
            H_minus_L(i,j) = (w_j/w_i) *                        &
                             (one-fermi_function(i,j,distij)) * & 
                             parameters%a_ij(target_%map_atomtypes(j))/distij

            ! j - i elements : K_ij * w_i/w_j
            H_minus_L(j,i) =  (w_i/w_j) *                        &  
                              (one-fermi_function(j,i,distij)) * & 
                               parameters%a_ij(target_%map_atomtypes(i))/distij
         enddo
      enddo

      do i=1,target_%n_atoms
         H_minus_L(i,i) = - sum(H_minus_L(i,:))
      enddo

   end subroutine construct_H_minus_L_matrix


   !> Subroutine for constructing A_qq matrix 
   !!    See J. Phys. Chem. C, 2021, 125, 23848-23863
   !! 
   !!    Input  : target_      -- mode
   !!    Input  : K-P matrix   
   !!    Output : A_qq
   subroutine construct_Aqq_matrix(target_, K_minus_P,A_qq)

      implicit none
   
      !input/output variables
      class(target_type), intent(in)                            :: target_
      real(dp), dimension(target_%n_q,target_%n_q), intent(in)  :: K_minus_P
      real(dp), dimension(target_%n_q,target_%n_q), intent(out) :: A_qq

      !internal variables
      real(dp), dimension(:,:), allocatable :: T_qq

      call mem_man%alloc(T_qq, target_%n_q, target_%n_q, "T_qq")

      call construct_static_t_qq(target_,T_qq)
      call dgemm('n','n',         &
                 target_%n_q,     &
                 target_%n_q,     &
                 target_%n_q,     &
                 one,             &
                 K_minus_P,       &
                 target_%n_q,     &
                 T_qq,            &
                 target_%n_q,     &
                 zero,            &
                 A_qq,            &
                 target_%n_q) ! internal_matrix*T_qq

      If(out_%ivrb.ge.4) &
         call out_%print_matrix("A_qmu Matrix",A_qq,target_%n_q,target_%n_q)

      call mem_man%dealloc(T_qq, "T_qq")

   end subroutine construct_Aqq_matrix


   !> Subroutine for constructing A_qmu matrix 
   !!    See J. Phys. Chem. C, 2021, 125, 23848-23863
   !! 
   !!    Input  : target_          -- mode
   !!    Input  : K-P matrix   
   !!    Input  : T_muq (optional) -- static T_muq matrix
   !!    Output : A_qmu
   subroutine construct_Aqmu_matrix(target_,K_minus_P,A_qmu,T_muq)

      implicit none
       
      !input/output variables
      class(target_type), intent(in)                                :: target_
      real(dp), dimension(target_%n_q,target_%n_q), intent(in)      :: K_minus_P
      real(dp), dimension(target_%n_q,target_%n_mu), intent(out)    :: A_qmu
      real(dp), dimension(target_%n_mu, target_%n_q), intent(in), &
                                                        optional    :: T_muq

      !internal variables
      real(dp), dimension(:,:), allocatable :: T_qmu

      if(present(T_muq)) then 
         call dgemm('n','t',      &
                    target_%n_q,  &
                    target_%n_mu, &
                    target_%n_q,  &
                    one,          &
                    K_minus_P,    &
                    target_%n_q,  &
                    T_muq,        &
                    target_%n_mu, &
                    zero,         &
                    A_qmu,        &
                    target_%n_q) ! K-P*T_qmu
      else ! not passe optional argument
         call mem_man%alloc(T_qmu,target_%n_q,target_%n_mu,"T_qmu")
         call construct_static_t_qmu(target_,T_qmu)
         call dgemm('n','n',      &
                    target_%n_q,  &
                    target_%n_mu, &
                    target_%n_q,  &
                    one,          &
                    K_minus_P,    &
                    target_%n_q,  &
                    T_qmu,        &
                    target_%n_q,  &
                    zero,         &
                    A_qmu,        &
                    target_%n_q) ! K-P*T_qmu
         call mem_man%dealloc(T_qmu, "T_qmu")
      endif
      
      If(out_%ivrb.ge.4) &
         call out_%print_matrix("A_qmu Matrix",A_qmu,target_%n_q,target_%n_mu)

   end subroutine construct_Aqmu_matrix


   !> Subroutine for constructing A_qq matrix [ complex version for hetero ]
   !!    See Front. Photonics.
   !! 
   !!    Input  : target_          -- mode
   !!    Input  : scale_           -- generally 1/2
   !!    Input  : K+H bar matrix   (this is complex)
   !!    Output : A_qq             (this is complex)
   subroutine construct_Aqq_matrix_complex(target_,scale_, K_plus_H, A_qq)

      implicit none
      
      !input/output variables
      class(target_type), intent(in)                               :: target_
      real(dp)                                                     :: scale_
      complex(dp), dimension(target_%n_q,target_%n_q), intent(in)  :: K_plus_H
      complex(dp), dimension(target_%n_q,target_%n_q), intent(out) :: A_qq

      !internal variables
      complex(dp) :: scale_cmp
      real(dp), dimension(:,:), allocatable :: T_qq
      complex(dp), dimension(:,:), allocatable :: T_qq_complex

      call mem_man%alloc(T_qq,target_%n_q,target_%n_q,"T_qq")
      call mem_man%alloc(T_qq_complex,target_%n_q,target_%n_q,"T_qq_complex")
      call construct_static_t_qq(target_,T_qq)

      T_qq_complex = dcmplx ( T_qq, zero)
      scale_cmp    = dcmplx(scale_,zero)

      call zgemm('n','n',          &
                 target_%n_q,      &
                 target_%n_q,      &
                 target_%n_q,      &
                 scale_cmp,        &
                 K_plus_H,         &
                 target_%n_q,      &
                 T_qq_complex,     &
                 target_%n_q,      &
                 dcmplx(zero,zero),&
                 A_qq,             &
                 target_%n_q) ! (K+P)*T_qq*scale_
      
      call mem_man%dealloc(T_qq, "T_qq")
      call mem_man%dealloc(T_qq_complex, "T_qq_complex")

   end subroutine construct_Aqq_matrix_complex


   !> Subroutine for constructing A_qmu matrix [ complex version for hetero]
   !!    See Front. Photonics 
   !! 
   !!    Input  : target_          -- model
   !!    Input  : scale_           -- generally 1/2
   !!    Input  : K+H bar          -- complex  
   !!    Input  : T_muq (optional) -- static T_muq matrix
   !!    Output : A_qmu            -- complex
   subroutine construct_Aqmu_matrix_complex(target_,scale_,K_plus_H,A_qmu,T_muq)

      implicit none
      
      !input/output variables
      class(target_type), intent(in)                                :: target_
      real(dp), intent(in)                                          :: scale_
      complex(dp), dimension(target_%n_q,target_%n_q), intent(in)   :: K_plus_H
      real(dp), dimension(target_%n_mu, target_%n_q), intent(in), &
                                                           optional :: T_muq
      complex(dp), dimension(target_%n_q,target_%n_mu), intent(out) :: A_qmu

      !internal variables
      complex(dp) :: scale_cmp
      real(dp), dimension(:,:), allocatable :: T_qmu
      complex(dp), dimension(:,:), allocatable :: T_qmu_complex

      call mem_man%alloc(T_qmu_complex,target_%n_q,target_%n_mu,"T_qmu_complex")

      if(present(T_muq)) then 
         T_qmu_complex = dcmplx(transpose(T_muq), zero)
      else 
         call mem_man%alloc(T_qmu,target_%n_q,target_%n_mu,"T_qmu")
         call construct_static_t_qmu(target_,T_qmu)
         T_qmu_complex = dcmplx(T_qmu, zero)
         call mem_man%dealloc(T_qmu, "T_qmu")
      endif

      scale_cmp = dcmplx(scale_,zero)
      call zgemm('n','n',           &
                 target_%n_q,       &
                 target_%n_mu,      &
                 target_%n_q,       &
                 scale_cmp,         &
                 K_plus_H,          &
                 target_%n_q,       &
                 T_qmu_complex,     &
                 target_%n_q,       &
                 dcmplx(zero,zero), &
                 A_qmu,             &
                 target_%n_q) ! scale_*(H+K)*T_qmu

      call mem_man%dealloc(T_qmu_complex, "T_qmu_complex")
     
   end subroutine construct_Aqmu_matrix_complex


   !> Function for calculating the Fermi function between i-j atoms
   !!    Damping Function recalling TS eq. 12
   !! 
   !!    Input  : i, j             -- Atom indeces
   !!    Input  : rij              -- i-j distance (in Bohr)
   !!    Output : Fermi_function   -- Value of the function
   real(dp) function fermi_function(i,j,rij)
   
      implicit none
      
      !input/output variables
      integer  :: i,j 
      real(dp) :: rij

      !internal variables
      real(dp) :: rij0
      real(dp) :: d_ij
      real(dp) :: s_ij

      d_ij = parameters%fermi_d(target_%map_atomtypes(i), &
                                target_%map_atomtypes(j))
      s_ij = parameters%fermi_s(target_%map_atomtypes(i), &
                                target_%map_atomtypes(j))
      rij0 = target_%nearest_distance(i,j)

      fermi_function = (one/(one+exp(-d_ij * ((rij/(s_ij*rij0))-one))))

   end function fermi_function

!-------------------------------------------------------------------------------
!  Subroutines for GMRES-related operations (on the fly)
!-------------------------------------------------------------------------------

   !> Subroutine for calculating Tqq dot x (where x comes from GMRES)
   !!    See J. Phys. Chem. C, 2021, 125, 23848-23863
   !! 
   !!    Input  : target_          -- model
   !!    Input  : x   
   !!    In/Out : M_dot_x          -- Matrix dot x
   subroutine product_q_block_x(target_,x,M_dot_x)

      implicit none
      
      !input/output variables
      class(target_type), intent(in)                     :: target_
      complex(dp), dimension(target_%n_q), intent(in)    :: x
      complex(dp), dimension(target_%n_q), intent(inout) :: M_dot_x

      !internal variables
      integer  :: i, j
      real(dp) :: distij
      real(dp) :: etaij, dij
      complex(dp), dimension(:), allocatable :: temp_array

      call mem_man%alloc(temp_array,target_%n_q,"temp_array_q_block_x")
 
      !Perform the matrix-vector multiplication on the fly
      !See construct_static_t_qq for the definition of the matrices 
      if(target_%kernel.eq.'coulomb') then
         !$omp parallel private(i,j,distij,dij) firstprivate(temp_array)
         !$omp do schedule(guided)
         do i = 1, target_%n_q
            temp_array(i) = temp_array(i) +  & 
                          dcmplx(target_%eta(target_%map_atomtypes(i)),zero)*x(i)
            do j = 1, i-1
               distij = dsqrt((target_%coord(1,i)-target_%coord(1,j))**2 + &
                              (target_%coord(2,i)-target_%coord(2,j))**2 + &
                              (target_%coord(3,i)-target_%coord(3,j))**2 )
               temp_array(i) = temp_array(i)+dcmplx(one/distij,zero)*x(j)
               temp_array(j) = temp_array(j)+dcmplx(one/distij,zero)*x(i)
            enddo
         enddo
         !$omp end do
         !$omp critical
         do i = 1, target_%n_q
            M_dot_x(i) = M_dot_x(i) + temp_array(i)
         enddo
         !$omp end critical
         !$omp end parallel
      else if(target_%kernel.eq.'ohno') then
         !$omp parallel private(i,j,distij,etaij,dij) firstprivate(temp_array)
         !$omp do schedule(guided)
         do i = 1, target_%n_q
            temp_array(i) = temp_array(i) +  &
                          dcmplx(target_%eta(target_%map_atomtypes(i)),zero)*x(i)
            do j = 1, i-1
               distij = dsqrt((target_%coord(1,i)-target_%coord(1,j))**2 + &
                              (target_%coord(2,i)-target_%coord(2,j))**2 + &
                              (target_%coord(3,i)-target_%coord(3,j))**2 )
               etaij = half*(target_%eta(target_%map_atomtypes(i))+ &
                             target_%eta(target_%map_atomtypes(j)))
               dij   = (one+(etaij**two)*(distij**two))**(half)
               temp_array(i) = temp_array(i)+dcmplx(etaij/dij,zero)*x(j)
               temp_array(j) = temp_array(j)+dcmplx(etaij/dij,zero)*x(i)
            enddo
         enddo
         !$omp end do
         !$omp critical
         do i = 1, target_%n_q
            M_dot_x(i) = M_dot_x(i) + temp_array(i)
         enddo
         !$omp end critical
         !$omp end parallel
      else if(target_%kernel.eq.'gaussian') then
         !$omp parallel private(i,j,distij,dij) firstprivate(temp_array)
         !$omp do schedule(guided)
         do i = 1, target_%n_q
            temp_array(i) = temp_array(i) +  &
                          dcmplx(target_%eta(target_%map_atomtypes(i)),zero)*x(i)
            do j = 1, i-1
               distij = dsqrt((target_%coord(1,i)-target_%coord(1,j))**2 + &
                              (target_%coord(2,i)-target_%coord(2,j))**2 + &
                              (target_%coord(3,i)-target_%coord(3,j))**2 )
               dij = erf(distIJ/ &
                           dsqrt(target_%r_q(target_%map_atomtypes(i))**2 + &
                                 target_%r_q(target_%map_atomtypes(j))**2))/ &
                                 distIJ
               temp_array(i) = temp_array(i)+dcmplx(dij,zero)*x(j)
               temp_array(j) = temp_array(j)+dcmplx(dij,zero)*x(i)
            enddo
         enddo
         !$omp end do
         !$omp critical
         do i = 1, target_%n_q
            M_dot_x(i) = M_dot_x(i) + temp_array(i)
         enddo
         !$omp end critical
         !$omp end parallel
      endif
      
      call mem_man%dealloc(temp_array,"temp_array_q_block_x")

   end subroutine product_q_block_x


   !> Subroutine for calculating Tqmu x, Tmuq x, Tmumu x (x comes from GMRES)
   !!    See J. Phys. Chem. C, 2021, 125, 23848-23863
   !! 
   !!    Input  : target_          -- model
   !!    Input  : x   
   !!    In/Out : M_dot_x          -- Matrix dot x
   subroutine product_mu_block_x(target_,x,M_dot_x)

      implicit none
      
      !input/output variables
      class(target_type), intent(in)                       :: target_
      complex(dp), dimension(target_%n_var), intent(in)    :: x
      complex(dp), dimension(target_%n_var), intent(inout) :: M_dot_x
 
      !internal variables
      integer  :: i, j
      integer  :: index_1
      real(dp) :: distij
      real(dp) :: distij_2
      real(dp) :: distij_3
      real(dp) :: distij_5
      real(dp) :: f_erf_ij
      real(dp) :: f_erf_ij_2
      real(dp) :: f_erf_ij_3
      real(dp) :: fD1
      real(dp) :: fD2
      real(dp) :: x_ij, y_ij, z_ij
      real(dp) :: tqmux, tqmuy, tqmuz
      real(dp) :: Txx, Tyy, Tzz, Txy, Txz, Tyz
      complex(dp), dimension(:), allocatable   :: temp_array

      ! Alloca array principale
      call mem_man%alloc(temp_array, target_%n_var, "temp_array_mu_block_x")
      
      !Perform the matrix-vector multiplication on the fly
      !See construct_static_t_qmu, t_muq, t_mumu for the definition 
      !of the matrices 
      !$omp parallel private(i,j,x_ij,y_ij,z_ij,distij,f_erf_ij,tqmux,tqmuy,tqmuz) &
      !$omp private(fd1,fd2,txx,txy,txz,tyy,tyz,tzz) firstprivate(temp_array)
      !$omp do schedule(guided)
      do i = 1, target_%n_atoms
         do j = 1, i-1

            x_ij   = target_%coord(1,i)-target_%coord(1,j)
            y_ij   = target_%coord(2,i)-target_%coord(2,j)
            z_ij   = target_%coord(3,i)-target_%coord(3,j)
            distij   = dsqrt(x_ij**2 + y_ij**2 + z_ij**2)
            distij_2 = distij**two
            distij_3 = distij**three
            distij_5 = distij**five
            f_erf_IJ   = dsqrt(target_%r_q(target_%map_atomtypes(i))**2 + &
                               target_%r_mu(target_%map_atomtypes(j))**2)
            f_erf_IJ_2 = f_erf_IJ**two
            index_1 = 3*(j-1) 
            TqmuX = x_ij/distij_3 * (derf(distij/f_erf_ij) - ((two*distij)/ & 
                               (dsqrt(pi)*f_erf_ij))*dexp(-distij_2/f_erf_ij_2))
            TqmuY = y_ij/distij_3 * (derf(distij/f_erf_ij) - ((two*distij)/ &
                               (dsqrt(pi)*f_erf_ij))*dexp(-distij_2/f_erf_ij_2))
            TqmuZ = z_ij/distij_3 * (derf(distij/f_erf_ij) - ((two*distij)/ &
                               (dsqrt(pi)*f_erf_ij))*dexp(-distij_2/f_erf_ij_2))
 
            !T^{q mu} x
            temp_array(i) = temp_array(i) + TqmuX * x(3*(j-1)+target_%n_q+1) + &
                                            TqmuY * x(3*(j-1)+target_%n_q+2) + &
                                            TqmuZ * x(3*(j-1)+target_%n_q+3)
 
            !T^{mu q} x
            temp_array(3*(j-1)+target_%n_q+1) = & 
            temp_array(3*(j-1)+target_%n_q+1) + tqmux * x(i)
            temp_array(3*(j-1)+target_%n_q+2) = &
            temp_array(3*(j-1)+target_%n_q+2) + tqmuy * x(i)
            temp_array(3*(j-1)+target_%n_q+3) = &
            temp_array(3*(j-1)+target_%n_q+3) + tqmuz * x(i)
 
            !The transpose is not the same for heterostructures
            f_erf_IJ   = dsqrt(target_%r_q(target_%map_atomtypes(j))**2 + &
                               target_%r_mu(target_%map_atomtypes(i))**2)
            f_erf_IJ_2 = f_erf_IJ**two

            TqmuX = x_ij/distij_3 * (derf(distij/f_erf_ij) - ((two*distij)/ &
                               (dsqrt(pi)*f_erf_ij))*dexp(-distij_2/f_erf_ij_2))
            TqmuY = y_ij/distij_3 * (derf(distij/f_erf_ij) - ((two*distij)/ &
                               (dsqrt(pi)*f_erf_ij))*dexp(-distij_2/f_erf_ij_2))
            TqmuZ = z_ij/distij_3 * (derf(distij/f_erf_ij) - ((two*distij)/ &
                               (dsqrt(pi)*f_erf_ij))*dexp(-distij_2/f_erf_ij_2))
 
            !T^{q mu} x
            temp_array(j) = temp_array(j) - TqmuX * x(3*(i-1)+target_%n_q+1) - &
                                            TqmuY * x(3*(i-1)+target_%n_q+2) - &
                                            TqmuZ * x(3*(i-1)+target_%n_q+3)
 
            !T^{mu q} x
            temp_array(3*(i-1)+target_%n_q+1) = &
            temp_array(3*(i-1)+target_%n_q+1) - tqmux * x(j)
            temp_array(3*(i-1)+target_%n_q+2) = &
            temp_array(3*(i-1)+target_%n_q+2) - tqmuy * x(j)
            temp_array(3*(i-1)+target_%n_q+3) = &
            temp_array(3*(i-1)+target_%n_q+3) - tqmuz * x(j)
 
            !T^{mu mu} x
            f_erf_IJ = dsqrt(target_%r_mu(target_%map_atomtypes(i))**2 + &
                             target_%r_mu(target_%map_atomtypes(j))**2)
            f_erf_IJ_2 = f_erf_IJ**two
            f_erf_IJ_3 = f_erf_IJ**three

            fD1 = derf(distIJ/f_erf_IJ) - ((two*distIJ)/(dsqrt(pi)*f_erf_IJ))* &
                  dexp(-distIJ_2/f_erf_IJ_2)
            fD2 = (four/(sqrt(pi)*f_erf_IJ_3)) * dexp(-distIJ_2/f_erf_IJ_2)

            Txx = ((-three*x_ij*x_ij)/distij_5+one/distij**3)*fD1 + &
                    fD2*(x_ij*x_ij)/distij_2
            Txy = ((-three*x_ij*y_ij)/distij_5              )*fD1 + &
                    fD2*(x_ij*y_ij)/distij_2
            Txz = ((-three*x_ij*z_ij)/distij_5              )*fD1 + &
                    fD2*(x_ij*z_ij)/distij_2
            Tyy = ((-three*y_ij*y_ij)/distij_5+one/distij**3)*fD1 + &
                    fD2*(y_ij*y_ij)/distij_2
            Tyz = ((-three*y_ij*z_ij)/distij_5              )*fD1 + &
                    fD2*(y_ij*z_ij)/distij_2
            Tzz = ((-three*z_ij*z_ij)/distij_5+one/distij**3)*fD1 + &
                    fD2*(z_ij*z_ij)/distij_2
 
            !Symmetrization
            temp_array(3*(i-1)+target_%n_q+1) = &
            temp_array(3*(i-1)+target_%n_q+1) + Txx * x(3*(j-1)+target_%n_q+1)
            temp_array(3*(i-1)+target_%n_q+1) = &
            temp_array(3*(i-1)+target_%n_q+1) + Txy * x(3*(j-1)+target_%n_q+2)
            temp_array(3*(i-1)+target_%n_q+1) = &
            temp_array(3*(i-1)+target_%n_q+1) + Txz * x(3*(j-1)+target_%n_q+3)
            temp_array(3*(i-1)+target_%n_q+2) = &
            temp_array(3*(i-1)+target_%n_q+2) + Txy * x(3*(j-1)+target_%n_q+1)
            temp_array(3*(i-1)+target_%n_q+2) = &
            temp_array(3*(i-1)+target_%n_q+2) + Tyy * x(3*(j-1)+target_%n_q+2)
            temp_array(3*(i-1)+target_%n_q+2) = &
            temp_array(3*(i-1)+target_%n_q+2) + Tyz * x(3*(j-1)+target_%n_q+3)
            temp_array(3*(i-1)+target_%n_q+3) = &
            temp_array(3*(i-1)+target_%n_q+3) + Txz * x(3*(j-1)+target_%n_q+1)
            temp_array(3*(i-1)+target_%n_q+3) = &
            temp_array(3*(i-1)+target_%n_q+3) + Tyz * x(3*(j-1)+target_%n_q+2)
            temp_array(3*(i-1)+target_%n_q+3) = &
            temp_array(3*(i-1)+target_%n_q+3) + Tzz * x(3*(j-1)+target_%n_q+3)

            temp_array(3*(j-1)+target_%n_q+1) = &
            temp_array(3*(j-1)+target_%n_q+1) + Txx * x(3*(i-1)+target_%n_q+1)
            temp_array(3*(j-1)+target_%n_q+1) = &
            temp_array(3*(j-1)+target_%n_q+1) + Txy * x(3*(i-1)+target_%n_q+2)
            temp_array(3*(j-1)+target_%n_q+1) = &
            temp_array(3*(j-1)+target_%n_q+1) + Txz * x(3*(i-1)+target_%n_q+3)
            temp_array(3*(j-1)+target_%n_q+2) = &
            temp_array(3*(j-1)+target_%n_q+2) + Txy * x(3*(i-1)+target_%n_q+1)
            temp_array(3*(j-1)+target_%n_q+2) = &
            temp_array(3*(j-1)+target_%n_q+2) + Tyy * x(3*(i-1)+target_%n_q+2)
            temp_array(3*(j-1)+target_%n_q+2) = &
            temp_array(3*(j-1)+target_%n_q+2) + Tyz * x(3*(i-1)+target_%n_q+3)
            temp_array(3*(j-1)+target_%n_q+3) = &
            temp_array(3*(j-1)+target_%n_q+3) + Txz * x(3*(i-1)+target_%n_q+1)
            temp_array(3*(j-1)+target_%n_q+3) = &
            temp_array(3*(j-1)+target_%n_q+3) + Tyz * x(3*(i-1)+target_%n_q+2)
            temp_array(3*(j-1)+target_%n_q+3) = &
            temp_array(3*(j-1)+target_%n_q+3) + Tzz * x(3*(i-1)+target_%n_q+3)      
         enddo
      enddo
      !$omp end do
      !$omp critical
      do i = 1, target_%n_var
         M_dot_x(i) = M_dot_x(i) + temp_array(i)
      enddo
      !$omp end critical
      !$omp end parallel

      call mem_man%dealloc(temp_array,"temp_array_mu_block_x")

   end subroutine product_mu_block_x


   !> Subroutine for the updating GMRES y vector by M_dot_x (charges)
   !!    See J. Phys. Chem. C, 2021, 125, 23848-23863
   !! 
   !!    Input  : target_          -- model
   !!    Input  : freq             -- for heterostructures    
   !!    Input  : M_dot_x          -- Matrix dot x 
   !!    Output : y                -- GMRES y solution
   subroutine update_y_q(target_, freq, M_dot_x, y)

      implicit none
      
      !input/output variables
      class(target_type), intent(in)                     :: target_
      real(dp), intent(in)                               :: freq
      complex(dp), dimension(target_%n_q), intent(in)    :: M_dot_x
      complex(dp), dimension(target_%n_q), intent(inout) :: y
 
      !internal variables
      integer  :: i, j
      real(dp) :: distij
      complex(dp) :: w_i
      complex(dp) :: w_j
      complex(dp) :: K_ij, K_ji
      complex(dp) :: H_ij, H_ji
      
      complex(dp), dimension(:), allocatable :: temp_array

      call mem_man%alloc(temp_array, target_%n_q, "temp_array_update_y_q")

      if(target_%heterogeneous) then ! I will have 1/2 [ H(w) + K ] * M_dot_x
         !$omp parallel private(i,j,distij,K_ij,K_ji,H_ij,H_ji,w_i,w_j) &
         !$omp firstprivate(temp_array) 
         !$omp do schedule(guided)
         do i = 1, target_%n_q
            ! w_i = 2 * n_i * tau_i / (1 - i w tau_i)
            w_i = two * parameters%density(target_%map_atomtypes(i)) /       &
                        dcmplx(one/parameters%tau(target_%map_atomtypes(i)), &
                               - freq)
            do j = 1, i-1
               ! w_j = 2 * n_j * tau_j / (1 - i w tau_j)
               distij    = sqrt((target_%coord(1,i)-target_%coord(1,j))**2 + &
                                (target_%coord(2,i)-target_%coord(2,j))**2 + &
                                (target_%coord(3,i)-target_%coord(3,j))**2 )
               K_ij = dcmplx( (one-fermi_function(i,j,distij))* &
                      (parameters%a_ij(target_%map_atomtypes(i)))/distij, zero)
               K_ji = dcmplx( (one-fermi_function(j,i,distij))* &
                      (parameters%a_ij(target_%map_atomtypes(j)))/distij, zero)
               w_j  = two * parameters%density(target_%map_atomtypes(j)) /   &
                         dcmplx(one/parameters%tau(target_%map_atomtypes(j)),&
                                -freq)
               H_ij = (w_j / w_i) * K_ji
               H_ji = (w_i / w_j) * K_ij
               temp_array(i) = temp_array(i) + half * (K_ij + H_ij) *  &
                               (M_dot_x(j) - M_dot_x(i))
               temp_array(j) = temp_array(j) + half * (K_ji + H_ji) *  &
                               (M_dot_x(i) - M_dot_x(j))
            enddo
         enddo
         !$omp end do
         !$omp critical
         do i = 1,target_%n_q
            y(i) = y(i) + temp_array(i)
         enddo
         !$omp end critical
         !$omp end parallel 
      else !Standard Case : (K-P)*M_dot_x
         !$omp parallel private(i,j,distij,k_ij) firstprivate(temp_array)
         !$omp do schedule(guided)
         do i = 1, target_%n_q
            do j = 1, i-1
               distij = dsqrt((target_%coord(1,i)-target_%coord(1,j))**2 + &
                              (target_%coord(2,i)-target_%coord(2,j))**2 + &
                              (target_%coord(3,i)-target_%coord(3,j))**2 )
               !K_ij is symmetric in this case
               K_ij = dcmplx((one-fermi_function(i,j,distij))* &
                      (parameters%a_ij(target_%map_atomtypes(i)))/distij, zero)
               temp_array(i) = temp_array(i) + K_ij * (M_dot_x(j) - M_dot_x(i))
               temp_array(j) = temp_array(j) + K_ij * (M_dot_x(i) - M_dot_x(j))
            enddo
         enddo
         !$omp end do
         !$omp critical
         do i = 1,target_%n_q
            y(i) = y(i) + temp_array(i)
         enddo
         !$omp end critical
         !$omp end parallel 
      endif

      call mem_man%dealloc(temp_array, "temp_array_update_y_q")

   end subroutine update_y_q


   !> Subroutine for the updating GMRES y vector by M_dot_x (dipoles)
   !!    See J. Phys. Chem. C, 2021, 125, 23848-23863
   !! 
   !!    Input  : target_          -- model
   !!    Input  : M_dot_x          -- Matrix dot x 
   !!    Output : y                -- GMRES y solution
   subroutine update_y_mu(target_,M_dot_x, y)

      implicit none
      
      !input/output variables
      class(target_type), intent(in)                      :: target_
      complex(dp), dimension(target_%n_var), intent(in)    :: M_dot_x
      complex(dp), dimension(target_%n_var), intent(inout) :: y

      !internal variables
      integer :: i
      
      !$omp parallel do schedule(static)
      do i = target_%n_q + 1, target_%n_var
         y(i) = M_dot_x(i)
      enddo
      !$omp end parallel do

   end subroutine update_y_mu

!-------------------------------------------------------------------------------
!  Subroutines for Analysis - Post-Process
!-------------------------------------------------------------------------------


   !> Subroutine for calculating the induced field due to wfq
   !!    See ACS Photonics, 2022, 9, 3064
   !! 
   !!    Input  : target_      -- model
   !!    Input  : index_       -- Index of the atom
   !!    Input  : d_IJ         -- distance between atom and grid point(in Bohr)
   !!    Input  : variable_w   -- The charge
   !!    In/Out : EField       -- Electric field
   subroutine calculate_induced_field_q(target_,index_,d_IJ, variable_w, EField)

      implicit none

      !input/output variables
      class(target_type), intent(in)         :: target_
      integer, intent(in)                    :: index_     ! index of the atom
      real(dp), dimension(3), intent(in)     :: d_IJ       ! i-point distance
      complex(dp), intent(in)                :: variable_w ! i-th charge
      complex(dp), dimension(3), intent(out) :: EField     ! Field
 
      !internal variables
      real(dp) :: R_Q_I, R_Q_I_2
      real(dp) :: distIJ, distIJ_2,distIJ_3
      real(dp) :: const
      real(dp) :: factor
        
      distIJ = dsqrt(d_IJ(1)**2 + d_IJ(2)**2 + d_IJ(3)**2)

      if(distIJ.gt.1.0d-1) then
         distIJ_2 = distIJ**two
         distIJ_3 = distIJ*distIJ_2

         R_Q_I  = target_%r_q(target_%map_atomtypes(index_))
         R_Q_I_2 = R_Q_I**two

         const = two/(sqrt(pi)*R_Q_I)
         factor = derf(distIJ/R_Q_I) - const *distIJ * dexp(-distIJ_2/R_Q_I_2)

         EField(:) = EField(:) - factor*d_IJ(:)*variable_w/distIJ_3
      endif

   end subroutine calculate_induced_field_q


   !> Subroutine for calculating the induced field due to wfmu
   !!    See ACS Photonics, 2022, 9, 3064
   !! 
   !!    Input  : target_      -- model
   !!    Input  : index_       -- Index of the atom
   !!    Input  : d_IJ         -- distance between atom and grid point(in Bohr)
   !!    Input  : variable_w   -- The dipole
   !!    In/Out : EField       -- Electric field
   subroutine calculate_induced_field_mu(target_,index_,d_IJ, variable_w, EField)

      implicit none
 
      !input/output variables
      class(target_type), intent(in)           :: target_
      integer, intent(in)                      :: index_     ! index
      real(dp), dimension(3), intent(in)       :: d_IJ       ! i-point distance
      complex(dp), dimension(3), intent(in)    :: variable_w ! i-th dipole
      complex(dp), dimension(3), intent(inout) :: EField     ! Field

      !internal variables
      real(dp) :: R_Mu_I, R_Mu_I_2, R_Mu_I_3
      real(dp) :: distIJ, distIJ_2,distIJ_3, distIJ_5
      real(dp) :: fD1, fD2
      real(dp) :: Txx, Txy, Txz
      real(dp) :: Tyy, Tyx, Tyz
      real(dp) :: Tzz, Tzx, Tzy
        
      distIJ = dsqrt(d_IJ(1)**2 + d_IJ(2)**2 + d_IJ(3)**2)
      
      if(distIJ.gt.1.0d-1) then
         distIJ_2 = distIJ**two
         distIJ_3 = distIJ*distIJ_2
         distIJ_5 = distIJ_2 * distIJ_3

         R_Mu_I    = target_%r_mu(target_%map_atomtypes(index_))
         R_Mu_I_2 = R_Mu_I**two
         R_Mu_I_3 = R_Mu_I**three

         fD1 = derf(distIJ/R_Mu_I) - ((two*distIJ)/(dsqrt(pi)*R_Mu_I))* &
                                       dexp(-distIJ_2/R_Mu_I_2)
         fD2 = (four/(sqrt(pi)*R_Mu_I_3)) * dexp(-distIJ_2/R_Mu_I_2)
      
         Txx = (((-three*d_IJ(1)*d_IJ(1))/distIJ_5 + one/distij_3)*fD1 + &
                fD2*(d_IJ(1)*d_IJ(1))/distIJ_2) 
         Txy = (((-three*d_IJ(1)*d_IJ(2))/distIJ_5               )*fD1 + & 
                fD2*(d_IJ(1)*d_IJ(2))/distIJ_2) 
         Txz = (((-three*d_IJ(1)*d_IJ(3))/distIJ_5               )*fD1 + & 
                fD2*(d_IJ(1)*d_IJ(3))/distIJ_2)
         Tyx = (((-three*d_IJ(2)*d_IJ(1))/distIJ_5               )*fD1 + & 
                fD2*(d_IJ(2)*d_IJ(1))/distIJ_2)
         Tyy = (((-three*d_IJ(2)*d_IJ(2))/distIJ_5 + one/distij_3)*fD1 + & 
                fD2*(d_IJ(2)*d_IJ(2))/distIJ_2) 
         Tyz = (((-three*d_IJ(2)*d_IJ(3))/distIJ_5               )*fD1 + & 
                fD2*(d_IJ(2)*d_IJ(3))/distIJ_2) 
         Tzx = (((-three*d_IJ(3)*d_IJ(1))/distIJ_5               )*fD1 + & 
                fD2*(d_IJ(3)*d_IJ(1))/distIJ_2)
         Tzy = (((-three*d_IJ(3)*d_IJ(2))/distIJ_5               )*fD1 + & 
                fD2*(d_IJ(3)*d_IJ(2))/distIJ_2) 
         Tzz = (((-three*d_IJ(3)*d_IJ(3))/distIJ_5 + one/distij_3)*fD1 + & 
                fD2*(d_IJ(3)*d_IJ(3))/distIJ_2)
         
         EFIeld(1) = EField(1) + variable_w(1) * Txx + &
                                 variable_w(2) * Txy + &
                                 variable_w(3) * Txz
         EFIeld(2) = EField(2) + variable_w(1) * Tyx + &
                                 variable_w(2) * Tyy + &
                                 variable_w(3) * Tyz 
         EFIeld(3) = EField(3) + variable_w(1) * Tzx + &
                                 variable_w(2) * Tzy + &
                                 variable_w(3) * Tzz
      endif
     
   end subroutine calculate_induced_field_mu


   !> Subroutine for calculating the plasmon density due to wfq
   !!    See ACS Photonics, 2022, 9, 3064
   !! 
   !!    Input  : target_      -- model
   !!    Input  : index_       -- Index of the atom
   !!    Input  : d_IJ         -- distance between atom and grid point(in Bohr)
   !!    Input  : variable_w   -- The charge
   !!    In/Out : Density      -- Plasmon Density
   subroutine calculate_density_q(target_,index_,d_IJ, variable_w, density)
 
      implicit none
       
      !input/output variables
      class(target_type), intent(in)     :: target_
      integer, intent(in)                :: index_              ! index atom
      real(dp), dimension(3), intent(in) :: d_IJ ! distance from point coord
      complex(dp), intent(in)            :: variable_w      ! i-th charge
      complex(dp), intent(inout)         :: density  ! density
 
      !internal variables
      real(dp) :: R_Q_I, R_Q_I_2
      real(dp) :: distIJ, distIJ_2
      real(dp) :: const
      real(dp) :: factor
         
      distIJ   = dsqrt(d_IJ(1)**2 + d_IJ(2)**2 + d_IJ(3)**2)
      distIJ_2 = distIJ**two
      R_Q_I    = target_%r_q(target_%map_atomtypes(index_))
      R_Q_I_2  = R_Q_I**two
      const    = one/((pi**(three/two))*R_Q_I*R_Q_I_2)
      factor   = const*dexp(-distIJ_2/R_Q_I_2)
 
      density = density + factor*variable_w
      
   end subroutine calculate_density_q


   !> Subroutine for calculating the plasmon density due to wfmu
   !!    See ACS Photonics, 2022, 9, 3064
   !! 
   !!    Input  : target_      -- model
   !!    Input  : index_       -- Index of the atom
   !!    Input  : d_IJ         -- distance between atom and grid point(in Bohr)
   !!    Input  : variable_w   -- The dipole
   !!    In/Out : Density      -- Plasmon Density
   subroutine calculate_density_mu(target_,index_,d_IJ, variable_w,density)
  
      implicit none
       
      !input/output variables
      class(target_type), intent(in)        :: target_
      integer, intent(in)                   :: index_     ! index
      real(dp), dimension(3), intent(in)    :: d_IJ       ! i-point distance
      complex(dp), dimension(3), intent(in) :: variable_w ! i-th dipole
      complex(dp), intent(inout)            :: density    ! Plasmon density
  
      !internal variables
      real(dp) :: R_Mu_I, R_Mu_I_2
      real(dp) :: distIJ, distIJ_2
      real(dp) :: const
      real(dp) :: factor
         
      distIJ   = dsqrt(d_IJ(1)**2 + d_IJ(2)**2 + d_IJ(3)**2)
      distIJ_2 = distIJ**two
      R_Mu_I   = target_%r_mu(target_%map_atomtypes(index_))
      R_Mu_I_2 = R_Mu_I**two
      const    = one/((pi**(three/two))*R_Mu_I*R_Mu_I_2)
      factor   = const*dexp(-distIJ_2/R_Mu_I_2)
       
      density = density - two*factor/R_Mu_I_2 * &  ! derivate wrt r
                          (variable_w(1) * d_IJ(1) + & 
                           variable_w(2) * d_IJ(2) + &
                           variable_w(3) * d_IJ(3)) 
       
   end subroutine calculate_density_mu


   !> Subroutine for calculating the field gradient due to wfq
   !!    This is performed by a two-point numerical differentiation of the 
   !!    electric field.
   !! 
   !!    Input  : target_      -- model
   !!    Input  : index_       -- Index of the atom
   !!    Input  : coord_atom   -- Coordinate of the atom
   !!    Input  : coord_point  -- Coordinate of the grid point
   !!    Input  : variable_w   -- The charge
   !!    In/Out : FieldGrad    -- Field Gradient
   subroutine calculate_field_gradient_num_q(target_,index_, coord_atom, &
                                             coord_point, variable_w, FieldGrad)
  
      implicit none
       
      !input/output variables
      class(target_type), intent(in)             :: target_
      integer, intent(in)                        :: index_      ! index
      real(dp), dimension(3), intent(in)         :: coord_atom  ! coord atom i
      real(dp), dimension(3), intent(in)         :: coord_point ! coord point 
      complex(dp), intent(in)                    :: variable_w  ! i-th charge
      complex(dp), dimension(3,3), intent(inout) :: FieldGrad   ! Field Gradient
  
      !internal variables
      complex(dp), dimension(3) :: EField_p, EField_m
      real(dp), dimension(3)    :: point_coord ! coord point (atoms index) copy
      real(dp) :: shift = 0.003d0 ! in Angstrom
      integer  :: i
      real(dp), dimension(3) :: d_IJ
  
      point_coord = coord_point
      do i = 1, 3 
         EField_p = zero
         point_coord(i) = point_coord(i) + shift
         d_IJ(1)   = (coord_atom(1)-point_coord(1))*tobohr
         d_IJ(2)   = (coord_atom(2)-point_coord(2))*tobohr
         d_IJ(3)   = (coord_atom(3)-point_coord(3))*tobohr
         call calculate_induced_field_q(target_, index_, d_IJ, variable_w, &
                                        EField_p)

         EField_m = zero
         point_coord(i) = point_coord(i) - two*shift
         d_IJ(1)   = (coord_atom(1)-point_coord(1))*tobohr
         d_IJ(2)   = (coord_atom(2)-point_coord(2))*tobohr
         d_IJ(3)   = (coord_atom(3)-point_coord(3))*tobohr
         call calculate_induced_field_q(target_, index_, d_IJ, variable_w, &
                                        EField_m)
  
         point_coord(i) = point_coord(i) + shift

         FieldGrad(:, i) =  FieldGrad(:, i) + (EField_p - EField_m)/ (two*shift)
      enddo
     
   end subroutine calculate_field_gradient_num_q


   !> Subroutine for calculating the field gradient due to wfq
   !!    This is performed by a two-point numerical differentiation of the 
   !!    electric field.
   !! 
   !!    Input  : target_      -- model
   !!    Input  : index_       -- Index of the atom
   !!    Input  : coord_atom   -- Coordinate of the atom
   !!    Input  : coord_point  -- Coordinate of the grid point
   !!    Input  : variable_w   -- The dipole
   !!    In/Out : FieldGrad    -- Field Gradient
   subroutine calculate_field_gradient_num_mu(target_,index_, coord_atom, &
                                              coord_point, variable_w, FieldGrad)
  
      implicit none
       
      !input/output variables
      class(target_type), intent(in)             :: target_
      integer, intent(in)                        :: index_      ! index
      real(dp), dimension(3), intent(in)         :: coord_atom  ! coord atom i 
      real(dp), dimension(3), intent(in)         :: coord_point ! coord point 
      complex(dp), dimension(3), intent(in)      :: variable_w  ! i-th dipole
      complex(dp), dimension(3,3), intent(inout) :: FieldGrad   ! Field Gradient 
  
      !internal variables
      complex(dp), dimension(3) :: EField_p, EField_m
      real(dp), dimension(3)    :: point_coord ! coord point (atoms index) copy
      real(dp) :: shift = 0.003d0 !in Angstrom
      integer  :: i
      real(dp), dimension(3) :: d_IJ
  
      point_coord = coord_point
      do i = 1, 3 
         EField_p = zero
         point_coord(i) = point_coord(i) + shift
         d_IJ(1)   = (coord_atom(1)-point_coord(1))*tobohr
         d_IJ(2)   = (coord_atom(2)-point_coord(2))*tobohr
         d_IJ(3)   = (coord_atom(3)-point_coord(3))*tobohr
         call calculate_induced_field_mu(target_, index_, d_IJ, variable_w, &
                                         EField_p)
  
         EField_m = zero
         point_coord(i) = point_coord(i) - two*shift
         d_IJ(1)   = (coord_atom(1)-point_coord(1))*tobohr
         d_IJ(2)   = (coord_atom(2)-point_coord(2))*tobohr
         d_IJ(3)   = (coord_atom(3)-point_coord(3))*tobohr
         call calculate_induced_field_mu(target_, index_, d_IJ, variable_w, &
                                         EField_m)
  
         point_coord(i) = point_coord(i) + shift
   
         FieldGrad(:, i) = FieldGrad(:, i) + (EField_p - EField_m)/ (two*shift)
      enddo
       
   end subroutine calculate_field_gradient_num_mu


   !> Subroutine for calculating the induced field due to BEM
   !! 
   !!    Input  : d_IJ         -- distance between atom and grid point(in Bohr)
   !!    Input  : variable_w   -- The charge
   !!    In/Out : EField       -- Electric field
   subroutine calculate_induced_field_q_bem(d_IJ, variable_w, EField)

      implicit none

      !input/output variables
      real(dp), dimension(3), intent(in)     :: d_IJ       ! i-point distance
      complex(dp), intent(in)                :: variable_w ! i-th charge
      complex(dp), dimension(3), intent(out) :: EField     ! Field
 
      !internal variables
      real(dp) :: distIJ, distIJ_2,distIJ_3
        
      distIJ = dsqrt(d_IJ(1)**2 + d_IJ(2)**2 + d_IJ(3)**2)

      if(distIJ.gt.1.0d-1) then 
         distIJ_2 = distIJ**two
         distIJ_3 = distIJ*distIJ_2
         EField(:) = EField(:) - d_IJ(:)*variable_w/distIJ_3
      endif

   end subroutine calculate_induced_field_q_bem

end module matrix_module
