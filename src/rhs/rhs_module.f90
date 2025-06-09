!> RHS module
!!
!! This module contains the subroutines for calculating quantities
!! and arrays related to RHS of the linear equations 
!!
!! Date         : 2025
!!
module rhs_module
   use output_module
   use parameters_module
   use string_manipulation_module
   use array_manipulation_module
   use field_module
   use target_module
   use memory_manager_module

   implicit none

   !static field
   public construct_static_field_rhs_q
   public construct_static_field_rhs_mu

   !dynamic field
   public construct_rhs_w_q
   public construct_rhs_w_q_complex
   public construct_rhs_w_q_on_the_fly
   public construct_rhs_w_q_complex_on_the_fly
   public construct_rhs_w_mu

   !printing
   public print_static_field_rhs
   public print_dynamic_field_rhs
    
contains

!-------------------------------------------------------------------------------
!  Static Field
!-------------------------------------------------------------------------------

   !> Subroutine for constructing static field RHS (charge)
   !!    Input  : target_      -- model 
   !!    Output : rhs_q        -- RHS-charge
   subroutine construct_static_field_rhs_q(target_, rhs_q)
  
      implicit none
       
      !input variables
      class(target_type), intent(in)                  :: target_
      real(dp), dimension(target_%n_q,3), intent(out) :: rhs_q
  
      !internal variables
      integer  :: i
      real(dp) :: x_min, y_min, z_min
      real(dp) :: dx, dy, dz
  
      if(field%rhs_form.eq.'potential') then
         x_min = zero
         y_min = zero
         z_min = zero
         do i = 1, target_%n_atoms
            if ( i .eq. 1)             x_min = target_%coord(1,1)
            if ( target_%coord(1,i).lt.x_min) x_min = target_%coord(1,i)
            if ( i .eq. 1)            y_min = target_%coord(2,1)
            if ( target_%coord(2,i).lt.y_min) y_min = target_%coord(2,i)
            if ( i .eq. 1)            z_min = target_%coord(3,1)
            if ( target_%coord(3,i).lt.z_min) z_min = target_%coord(3,i)
         enddo
         !there is no chi when static/dynamic response
         !$omp parallel do private(i,dx,dy,dz)
         do i=1,target_%n_atoms
            dx = abs(target_%coord(1,i) - x_min)
            dy = abs(target_%coord(2,i) - y_min)
            dz = abs(target_%coord(3,i) - z_min)
            rhs_q(i,1) =  - field%e_0*dx
            rhs_q(i,2) =  - field%e_0*dy
            rhs_q(i,3) =  - field%e_0*dz
         enddo
         !$omp end parallel do
      else if(field%rhs_form.eq.'field') then
        !there is no chi when static/dynamic response
        !$omp parallel do private(i)
         do i=1,target_%n_atoms
            rhs_q(i,1) = - field%e_0*target_%coord(1,i) 
            rhs_q(i,2) = - field%e_0*target_%coord(2,i) 
            rhs_q(i,3) = - field%e_0*target_%coord(3,i) 
         enddo
         !$omp end parallel do
      endif
  
   end subroutine construct_static_field_rhs_q


   !> Subroutine for constructing static field RHS (dipole)
   !!    Input  : target_      -- model 
   !!    Output : rhs_q        -- RHS-charge
   subroutine construct_static_field_rhs_mu(target_, rhs_mu)
  
      implicit none
       
      !input/output variables
      class(target_type), intent(in)                   :: target_
      real(dp), dimension(target_%n_mu,3), intent(out) :: rhs_mu
  
      !internal variables
      integer :: i,index_1
  
      if (is_equal(field%e_0,zero)) then
         !$omp parallel do private(i, index_1)
         do i=1,target_%n_q
            index_1 = 3*(i-1)
            rhs_mu(index_1+1,1) = field%e_0
            rhs_mu(index_1+2,2) = field%e_0
            rhs_mu(index_1+3,3) = field%e_0
         enddo
         !$omp end parallel do
      else
         !$omp parallel do private(i, index_1)
         do i=1,target_%n_q
            index_1 = 3*(i-1)
            rhs_mu(index_1+1,1) = one
            rhs_mu(index_1+2,2) = one
            rhs_mu(index_1+3,3) = one
         enddo
         !$omp end parallel do
      endif
       
   end subroutine construct_static_field_rhs_mu

!-------------------------------------------------------------------------------
!  Dynamic Field
!-------------------------------------------------------------------------------

   !> Subroutine for constructing dynamic field RHS (q)
   !!    Input  : target_      -- model 
   !!    Input  : rhs_0        -- static RHS
   !!    Input  : K_minus_P    -- K-P matrix
   !!    Output : rhs_w_q      -- RHS-charge(w)
   !!
   !!    RHS-q(w) = -(K-P)*V
   subroutine construct_rhs_w_q(target_, rhs_0, K_minus_P, rhs_w_q)
 
      implicit none
       
      !input/output variables
      class(target_type), intent(in)   :: target_
      real(dp), dimension(target_%n_q,target_%n_q), intent(in) :: K_minus_P
      real(dp), dimension(target_%n_q,3), intent(in)           :: rhs_0
      complex(dp), dimension(target_%n_q,3), intent(inout)     :: rhs_w_q
 
      !internal variables
      integer :: i
      real(dp), dimension(:,:), allocatable :: temporary_vector
 
      call mem_man%alloc(temporary_vector, target_%n_q, 3, &
                         "temporary_vector_rhs_w_q")
      call dgemm('N','N',           &
                  target_%n_q,      &
                  3,                &
                  target_%n_q,      &
                  -one,             &
                  K_minus_P,        &
                  target_%n_q,      &
                  rhs_0,            &
                  target_%n_q,      &
                  zero,             &
                  temporary_vector, &
                  target_%n_q)
 
      if (field%polarization .eq. 'all') then
         !$omp parallel do private(i)
         do i=1,target_%n_q
            rhs_w_q(i,:) = dcmplx(temporary_vector(i,:),zero)
         enddo
         !$omp end parallel do
      else 
         !$omp parallel do private(i)
         do i=1,target_%n_q
            rhs_w_q(i,field%polarization_index) = &
                    dcmplx(temporary_vector(i,field%polarization_index),zero)
         enddo
         !$omp end parallel do
      endif
      call mem_man%dealloc(temporary_vector, "temporary_vector_rhs_w_q")
 
   end subroutine construct_rhs_w_q


   !> Subroutine for constructing dynamic field RHS (q) on the fly
   !!    Input  : target_      -- model 
   !!    Input  : rhs_0        -- static RHS
   !!    Output : rhs_w_q      -- RHS-charge(w)
   !!
   !!    RHS-q(w) = -(K-P)*V [on the fly]
   subroutine construct_rhs_w_q_on_the_fly(target_, rhs_0, rhs_w_q)

      use matrix_module, only: fermi_function
       
      implicit none
       
      !input/output variables
      class(target_type), intent(in)   :: target_
      real(dp), dimension(target_%n_q,3), intent(in)           :: rhs_0
      complex(dp), dimension(target_%n_q,3), intent(inout)     :: rhs_w_q
  
      !internal variables
      integer  :: i,j 
      real(dp) :: distij
      real(dp), dimension(:,:), allocatable :: tmp_vector

      call mem_man%alloc(tmp_vector, target_%n_q, 3, "tmp_vector_rhs_w_q_otf")

      !$omp parallel private(i,j,distij) firstprivate(tmp_vector)
      !$omp do schedule(guided) 
      do i = 1, target_%n_atoms
         do j = 1, i-1
            distij    = sqrt((target_%coord(1,i)-target_%coord(1,j))**2 + &
                             (target_%coord(2,i)-target_%coord(2,j))**2 + &
                             (target_%coord(3,i)-target_%coord(3,j))**2 )
            tmp_vector(i,:) = tmp_vector(i,:) - &
                              (one-fermi_function(i,j,distij)) * &
                              ((parameters%a_ij(target_%map_atomtypes(i))) / &
                                distij)*(rhs_0(j,:)-rhs_0(i,:))
            tmp_vector(j,:) = tmp_vector(j,:) - &
                              (one-fermi_function(j,i,distij)) * &
                              ((parameters%a_ij(target_%map_atomtypes(j))) / &
                                distij)*(rhs_0(i,:)-rhs_0(j,:))
         enddo
      enddo
      !$omp end do
      if (field%polarization.eq.'all') then
         !$omp critical
         do i = 1, target_%n_q
            rhs_w_q(i,:) = rhs_w_q(i,:) + dcmplx(tmp_vector(i,:),zero)
         enddo
         !$omp end critical
      else
         !$omp critical
         do i = 1, target_%n_q
            rhs_w_q(i,field%polarization_index) = &
            rhs_w_q(i,field%polarization_index) + &
                           dcmplx(tmp_vector(i,field%polarization_index),zero)
         enddo
         !$omp end critical
      endif
      !$omp end parallel
      call mem_man%dealloc(tmp_vector, "tmp_vector_rhs_w_q_otf")

   end subroutine construct_rhs_w_q_on_the_fly


   !> Subroutine for constructing dynamic field RHS (q) for hetero
   !!    Input  : target_      -- model 
   !!    Input  : scale_       -- 1/2 generally
   !!    Input  : rhs_0        -- static RHS
   !!    Input  : K_plus_H     -- K+H bar matrix [this depends on the freq]
   !!    Output : rhs_w_q      -- RHS-charge(w)
   !!
   !!    RHS-q(w) = -(K+H)*V
   subroutine construct_rhs_w_q_complex(target_, scale_, rhs_0, K_plus_H, &
                                        rhs_w_q)

      implicit none
  
      !input/output variables
      class(target_type), intent(in)   :: target_
      real(dp), intent(in) :: scale_
      complex(dp), dimension(target_%n_q,target_%n_q), intent(in) :: K_plus_H
      complex(dp), dimension(target_%n_q,3), intent(in)           :: rhs_0
      complex(dp), dimension(target_%n_q,3), intent(inout)        :: rhs_w_q
  
      !internal variables
      complex(dp) :: scale_cmp 

      scale_cmp = dcmplx(scale_,zero)
      if (field%polarization .eq. 'all') then
         call zgemm('N','N',           &
                     target_%n_q,      &
                     3,                &
                     target_%n_q,      &
                     -scale_cmp,       &
                     K_plus_H,         &
                     target_%n_q,      &
                     rhs_0,            &
                     target_%n_q,      &
                     dcmplx(zero,zero),&
                     rhs_w_q,          &
                     target_%n_q)
      else
         call zgemm('N','N',                               &
                     target_%n_q,                          &
                     1,                                    &
                     target_%n_q,                          &
                     -scale_cmp,                           &
                     K_plus_H,                             &
                     target_%n_q,                          &
                     rhs_0(:,field%polarization_index),   &
                     target_%n_q,                          &
                     dcmplx(zero,zero),                    &
                     rhs_w_q(:,field%polarization_index), &
                     target_%n_q)
      endif

   end subroutine construct_rhs_w_q_complex


   !> Subroutine for constructing dynamic field RHS (q) for hetero
   !!    Input  : target_      -- model 
   !!    Input  : scale_       -- 1/2 generally
   !!    Input  : freq         -- frequency
   !!    Input  : rhs_0        -- static RHS
   !!    Output : rhs_w_q      -- RHS-charge(w)
   !!
   !!    RHS-q(w) = -(K+H)*V [on the fly]
   subroutine construct_rhs_w_q_complex_on_the_fly(target_, scale_, freq, &
                                                   rhs_0, rhs_w_q)
  
      use matrix_module, only: fermi_function
  
      implicit none
       
      !input/output variables
      class(target_type), intent(in)   :: target_
      real(dp), intent(in) :: scale_
      real(dp), intent(in) :: freq
      complex(dp), dimension(target_%n_q,3), intent(in)    :: rhs_0
      complex(dp), dimension(target_%n_q,3), intent(inout) :: rhs_w_q
  
      !internal variables
      integer  :: i,j
      real(dp) :: distij
      complex(dp) :: w_i
      complex(dp) :: w_j
      complex(dp) :: K_ij, K_ji
      complex(dp) :: H_ij, H_ji
      complex(dp), dimension(:,:), allocatable :: temporary_vector
  
      call mem_man%alloc(temporary_vector, target_%n_q, 3, &
                         "temporary_vector_rhs_w_q_cmplx_otf")
      !$omp parallel do private(i,j,distij,w_i,w_j,K_ij,K_ji,H_ij,H_ji) &
      !$omp reduction(+:temporary_vector)
      do i = 1, target_%n_atoms
         ! w_i = 2 * n_i * tau_i / (1 - i w tau_i)
         w_i = two * parameters%density(target_%map_atomtypes(i)) /   &
                     dcmplx(one/parameters%tau(target_%map_atomtypes(i)),-freq)
         do j = 1, i-1
            ! w_j = 2 * n_j * tau_j / (1 - i w tau_j)
            w_j = two * parameters%density(target_%map_atomtypes(j)) /      &
                        dcmplx(one/parameters%tau(target_%map_atomtypes(j)),&
                               -freq)
            distij    = sqrt((target_%coord(1,i)-target_%coord(1,j))**2 + &
                             (target_%coord(2,i)-target_%coord(2,j))**2 + &
                             (target_%coord(3,i)-target_%coord(3,j))**2 )
            K_ij = dcmplx( (one-fermi_function(i,j,distij)) * &
                   (parameters%a_ij(target_%map_atomtypes(i)))/distij, zero)
            K_ji = dcmplx( (one-fermi_function(j,i,distij)) * &
                   (parameters%a_ij(target_%map_atomtypes(j)))/distij, zero)
            H_ij = (w_j / w_i) * K_ji
            H_ji = (w_i / w_j) * K_ij
            temporary_vector(i,:) = temporary_vector(i,:) - scale_ * &
                                    (K_ij + H_ij) * (rhs_0(j,:) - rhs_0(i,:)) 
            temporary_vector(j,:) = temporary_vector(j,:) - scale_ * &
                                    (K_ji + H_ji) * (rhs_0(i,:) - rhs_0(j,:)) 
         enddo
      enddo
      !$omp end parallel do
  
      if (field%polarization .eq. 'all') then
         !$omp parallel do private(i)
         do i = 1, target_%n_q
            rhs_w_q(i,:) = temporary_vector(i,:)
         enddo
         !$omp end parallel do
      else
         !$omp parallel do private(i)
         do i = 1, target_%n_q
            rhs_w_q(i,field%polarization_index) = &
                    temporary_vector(i,field%polarization_index)
         enddo
         !$omp end parallel do
      endif

      call mem_man%dealloc(temporary_vector, &
                           "temporary_vector_rhs_w_q_cmplx_otf")
  
   end subroutine construct_rhs_w_q_complex_on_the_fly


   !> Subroutine for constructing dynamic field RHS (mu) 
   !!    Input  : target_      -- model 
   !!    Output : rhs_w_mu     -- RHS-dipole(w)
   !!
   !!    RHS-q(w) = E(w)
   subroutine construct_rhs_w_mu(target_, rhs_w_mu)
  
      implicit none
       
      !input/output variables
      class(target_type), intent(in)   :: target_
      complex(dp), dimension(target_%n_mu,3), intent(inout)     :: rhs_w_mu
  
      !internal variables
      integer :: i
      integer :: index_1
  
      if (field%polarization.eq.'all') then
         !$omp parallel do private(i, index_1)
         do i=1,target_%n_atoms
            index_1 = 3*(i-1)
            rhs_w_mu(index_1+1,1) = dcmplx(field%e_0,zero)
            rhs_w_mu(index_1+2,2) = dcmplx(field%e_0,zero)
            rhs_w_mu(index_1+3,3) = dcmplx(field%e_0,zero)
         enddo
         !$omp end parallel do
      else 
         !$omp parallel do private(i, index_1)
         do i=1,target_%n_atoms
            index_1 = 3*(i-1)
            rhs_w_mu(index_1+field%polarization_index,&
                             field%polarization_index) = &
                                                         dcmplx(field%e_0,zero)
         enddo
         !$omp end parallel do
      endif
       
   end subroutine construct_rhs_w_mu

!-------------------------------------------------------------------------------
!  Printing
!-------------------------------------------------------------------------------

   !> Subroutine for printing the RHS for a static field
   !!    Input  : target_      -- model 
   !!    Output : rhs          -- RHS-static
   subroutine print_static_field_rhs(target_, rhs)
  
      implicit none
       
      !input/output variables
      class(target_type), intent(in)   :: target_
      real(dp), dimension(target_%n_var,3) :: rhs
  
      !internal variables
      integer :: i,j
      character(len=14) :: format_1 = "(32x,'RHS ',/)"
      character(len=63) :: format_2 = "(10x,'X Component',7x,'Y Component'7x,&
                                      &'Z Component',/)"
      character(len=16) :: format_3 = "(I4,3(4x,E14.6))"
      
      write(out_%iunit,format_1)
      write(out_%iunit,format_2)
      do i = 1,target_%n_var
         write(out_%iunit,format_3) i, (rhs(i,j), j = 1,3)
      enddo
      write(out_%iunit,out_%sticks)
      flush(out_%iunit)
       
   end subroutine print_static_field_rhs


   !> Subroutine for printing the RHS for a dynamic field
   !!    Input  : target_      -- model 
   !!    Output : rhs_w        -- RHS(w)
   subroutine print_dynamic_field_rhs(target_, rhs_w)
  
      implicit none
       
      !input/output variables
      class(target_type), intent(in)   :: target_
      complex(dp), dimension(target_%n_var,3), intent(in) :: rhs_w
  
      !internal variables
      integer :: i,j
      character(len=43) :: format_1 = "(33x,'RHS Complex Response',/)"
      character(len=63) :: format_2 = "(13x,'X Component',14x,'Y Component'&
                                      &14x,'Z Component',/)"
      character(len=60) :: format_3 = "(1x,I4,2x,2(E10.3,2x,E10.3,1X,'|',1x),&
                                      &E10.3,2X,E10.3)"
      
      write(out_%iunit,format_1)
      write(out_%iunit,format_2)
      do i = 1,target_%n_var
         write(out_%iunit,format_3) i,(rhs_w(i,j),j=1,3)
      enddo
      write(out_%iunit,out_%sticks)
      flush(out_%iunit)
       
   end subroutine print_dynamic_field_rhs

end module rhs_module
