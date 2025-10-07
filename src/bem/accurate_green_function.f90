!> Submodule Accurate Green Function (BEM module)
!!
!! This module contains the subroutine that perform the following steps:
!!  1. flat_integration_diagonal ! diagonal elements
!!     2. quadrature rule
!!     3. tesserae_integration_polar_coordinates
!!        4. legendre_gauss_lobatto_integration(npol1, points1,weights1)
!!        5. legendre_gauss_lobatto_integration(npol2, points2,weights2)
!!  6. flat_integration_off_diagonal
!!  7. selection_of_refined_elements
!!  8. refinement_diagonal_elements
!!  9. refinement_off_diagonal_elements (function)
!!
!! Author       : Tommaso Giovannini
!! Date         : 2025
!!
submodule (bem_module) accurate_green_function
 
   implicit none
 
contains

   !> Subroutine for constructing the accurate Green function 
   !!    In/Out  : bem      -- bem_type
   module subroutine construct_accurate_green_function(bem)
  
      implicit none
       
      !input/output variables
      class(bem_type), intent(inout) :: bem
  
      call bem%flat_integration_diagonal()
      call bem%flat_integration_off_diagonal()
      call bem%selection_of_refined_elements()
      call bem%refinement_diagonal_elements()
  
   end subroutine construct_accurate_green_function


   !> Subroutine for flat integration diagonal part
   !! This is a subroutine valid for flat tesserae only in polar coordinates.
   !!    In/Out  : bem      -- bem_type
   module subroutine flat_integration_diagonal(bem)
  
      implicit none
       
      !input/output variables
      class(bem_type), intent(inout) :: bem
  
      !internal variables
      integer :: i, k, k_1 
      integer :: n_update
      real(dp), dimension(:), allocatable :: int_coord_x
      real(dp), dimension(:), allocatable :: int_coord_y
      real(dp), dimension(:,:), allocatable :: internal_matrix
      real(dp), dimension(:,:), allocatable :: tmp_coord
      real(dp), dimension(:,:), allocatable :: product_matrix
  
      call bem%quadrature_rule()
  
      call mem_man%alloc(int_coord_x,bem%npol1*bem%npol2*3,"int_coord_x")
      call mem_man%alloc(int_coord_y,bem%npol1*bem%npol2*3,"int_coord_y")
      call mem_man%alloc(bem%integration_weights, bem%npol1*bem%npol2*3, &
                         "bem%integration_weights")

      !do the tesserae integration polar coordinates
      call bem%tesserae_integration_polar_coordinates(int_coord_x,int_coord_y)
      !allocation matrix for int points and weight
      !internal_matrix corresponds to tri in quadpol_flat
      call mem_man%alloc(internal_matrix,bem%npol1*bem%npol2*3,3, &
                         "internal_matrix")
  
      do i = 1, bem%npol1*bem%npol2*3
         internal_matrix(i,1) = int_coord_x(i)
         internal_matrix(i,2) = int_coord_y(i)
         internal_matrix(i,3) = one - int_coord_x(i) - int_coord_y(i)
      enddo
       
      call mem_man%alloc(tmp_coord, 3,3, "tmp_coord")
      call mem_man%alloc(product_matrix, bem%npol1*bem%npol2*3,3, &
                        "product_matrix")
      call mem_man%alloc(bem%polar_coord, 3,bem%npol1*bem%npol2*3*bem%n_var, &
                         "bem%polar_coord")
      call mem_man%alloc(bem%distance_polar_coord, &
                         bem%npol1*bem%npol2*3*bem%n_var, &
                         "bem%distance_polar_coord")
  
      n_update = 0
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
         product_matrix = matmul(internal_matrix,tmp_coord)
         Do k = 1, bem%npol1*bem%npol2*3
            K_1 = k + n_update
            bem%polar_coord(1,k_1) = product_matrix(k,1)
            bem%polar_coord(2,k_1) = product_matrix(k,2)
            bem%polar_coord(3,k_1) = product_matrix(k,3)
            bem%distance_polar_coord(k_1) = sqrt((bem%coord(1,i) - &
                                                  bem%polar_coord(1,k_1))**2 + &
                                                 (bem%coord(2,i) - &
                                                  bem%polar_coord(2,k_1))**2 + &
                                                 (bem%coord(3,i) - &
                                                  bem%polar_coord(3,k_1))**2)
         enddo
         n_update = n_update + bem%npol1*bem%npol2*3
      enddo
      !printing
      if(out_%ivrb.ge.4) call out_%print_matrix("BEM Polar coord", &
                              bem%polar_coord,3,bem%npol1*bem%npol2*3*bem%n_var)
      !deallocation
      call mem_man%dealloc(int_coord_x, "int_coord_x")
      call mem_man%dealloc(int_coord_y, "int_coord_y")
      call mem_man%dealloc(internal_matrix, "internal_matrix")
      call mem_man%dealloc(tmp_coord, "tmp_coord")
      call mem_man%dealloc(product_matrix, "product_matrix")
  
   end subroutine flat_integration_diagonal


   !> Subroutine for quadrature rules setting for an unit triangular mesh
   !! Strang and Fix, "An Analysis of the Finite Element Method"
   !! Prentice Hall, 1973.
   !! Other quadrature rule can be implemented, with different precisions  
   !!    In/Out  : bem      -- bem_type
   module subroutine quadrature_rule(bem)
  
      implicit none
  
      !input/output variables
      class(bem_type), intent(inout) :: bem
  
      !allocation of bem%x_tab, y_tab, w_tab
      call mem_man%alloc(bem%x_tab, 28,"bem%x_tab")
      call mem_man%alloc(bem%y_tab, 28,"bem%y_tab")
      call mem_man%alloc(bem%w_tab, 28,"bem%w_tab")
  
      !Definition of parameters from Strang and Fix
      bem%x_tab(1)  = 0.3333333333333333d0 
      bem%x_tab(2)  = 0.9480217181434233d0 
      bem%x_tab(3)  = 0.02598914092828833d0
      bem%x_tab(4)  = 0.02598914092828833d0
      bem%x_tab(5)  = 0.8114249947041546d0
      bem%x_tab(6)  = 0.09428750264792270d0 
      bem%x_tab(7)  = 0.09428750264792270d0 
      bem%x_tab(8)  = 0.01072644996557060d0 
      bem%x_tab(9)  = 0.4946367750172147d0 
      bem%x_tab(10) = 0.4946367750172147d0 
      bem%x_tab(11) = 0.5853132347709715d0 
      bem%x_tab(12) = 0.2073433826145142d0 
      bem%x_tab(13) = 0.2073433826145142d0 
      bem%x_tab(14) = 0.1221843885990187d0 
      bem%x_tab(15) = 0.4389078057004907d0 
      bem%x_tab(16) = 0.4389078057004907d0 
      bem%x_tab(17) = 0.6779376548825902d0 
      bem%x_tab(18) = 0.6779376548825902d0 
      bem%x_tab(19) = 0.04484167758913055d0 
      bem%x_tab(20) = 0.04484167758913055d0 
      bem%x_tab(21) = 0.27722066752827925d0 
      bem%x_tab(22) = 0.27722066752827925d0 
      bem%x_tab(23) = 0.8588702812826364d0 
      bem%x_tab(24) = 0.8588702812826364d0 
      bem%x_tab(25) = zero
      bem%x_tab(26) = zero
      bem%x_tab(27) = 0.1411297187173636d0
      bem%x_tab(28) = 0.1411297187173636d0
  
      bem%y_tab(1)  = 0.3333333333333333d0
      bem%y_tab(2)  = 0.02598914092828833d0
      bem%y_tab(3)  = 0.9480217181434233d0
      bem%y_tab(4)  = 0.02598914092828833d0
      bem%y_tab(5)  = 0.09428750264792270d0
      bem%y_tab(6)  = 0.8114249947041546d0
      bem%y_tab(7)  = 0.09428750264792270d0
      bem%y_tab(8)  = 0.4946367750172147d0
      bem%y_tab(9)  = 0.01072644996557060d0
      bem%y_tab(10) = 0.4946367750172147d0
      bem%y_tab(11) = 0.2073433826145142d0
      bem%y_tab(12) = 0.5853132347709715d0
      bem%y_tab(13) = 0.2073433826145142d0
      bem%y_tab(14) = 0.4389078057004907d0
      bem%y_tab(15) = 0.1221843885990187d0
      bem%y_tab(16) = 0.4389078057004907d0
      bem%y_tab(17) = 0.04484167758913055d0
      bem%y_tab(18) = 0.27722066752827925d0
      bem%y_tab(19) = 0.6779376548825902d0
      bem%y_tab(20) = 0.27722066752827925d0
      bem%y_tab(21) = 0.6779376548825902d0
      bem%y_tab(22) = 0.04484167758913055d0
      bem%y_tab(23) = zero
      bem%y_tab(24) = 0.1411297187173636d0
      bem%y_tab(25) = 0.8588702812826364d0
      bem%y_tab(26) = 0.1411297187173636d0
      bem%y_tab(27) = 0.8588702812826364d0
      bem%y_tab(28) = zero
             
      bem%w_tab(1)  = 0.08797730116222190d0 
      bem%w_tab(2)  = 0.008744311553736190d0 
      bem%w_tab(3)  = 0.008744311553736190d0 
      bem%w_tab(4)  = 0.008744311553736190d0 
      bem%w_tab(5)  = 0.03808157199393533d0
      bem%w_tab(6)  = 0.03808157199393533d0
      bem%w_tab(7)  = 0.03808157199393533d0
      bem%w_tab(8)  = 0.01885544805613125d0
      bem%w_tab(9)  = 0.01885544805613125d0
      bem%w_tab(10) = 0.01885544805613125d0
      bem%w_tab(11) = 0.07215969754474100d0
      bem%w_tab(12) = 0.07215969754474100d0
      bem%w_tab(13) = 0.07215969754474100d0
      bem%w_tab(14) = 0.06932913870553720d0
      bem%w_tab(15) = 0.06932913870553720d0
      bem%w_tab(16) = 0.06932913870553720d0
      bem%w_tab(17) = 0.04105631542928860d0
      bem%w_tab(18) = 0.04105631542928860d0
      bem%w_tab(19) = 0.04105631542928860d0
      bem%w_tab(20) = 0.04105631542928860d0
      bem%w_tab(21) = 0.04105631542928860d0
      bem%w_tab(22) = 0.04105631542928860d0
      bem%w_tab(23) = 0.007362383783300573d0 
      bem%w_tab(24) = 0.007362383783300573d0 
      bem%w_tab(25) = 0.007362383783300573d0 
      bem%w_tab(26) = 0.007362383783300573d0 
      bem%w_tab(27) = 0.007362383783300573d0 
      bem%w_tab(28) = 0.007362383783300573d0
  
   end subroutine quadrature_rule


   !> Subroutine for the integration over tesserae in polar coordinates
   !!    In/Out  : bem      -- bem_type
   !!    In/Out  : int_coord_x -- integration coordinates x
   !!    In/Out  : int_coord_y -- integration coordinates y
   module subroutine tesserae_integration_polar_coordinates(bem,int_coord_x, &
                                                            int_coord_y)
  
      implicit none
       
      !input/output variables
      class(bem_type), intent(inout) :: bem
      real(dp), dimension(bem%npol1*bem%npol2*3), intent(inout) :: int_coord_x
      real(dp), dimension(bem%npol1*bem%npol2*3), intent(inout) :: int_coord_y
  
      !internal variables
      real(dp), dimension(:), allocatable :: points1
      real(dp), dimension(:), allocatable :: weight1
      real(dp), dimension(:), allocatable :: points2
      real(dp), dimension(:), allocatable :: weight2
      real(dp), dimension(:), allocatable :: rho
      real(dp), dimension(:), allocatable :: phi
      real(dp), dimension(:), allocatable :: tmp_phi
      real(dp), dimension(:), allocatable :: tmp_rho
      real(dp), dimension(:), allocatable :: tmp_x
      real(dp), dimension(:), allocatable :: tmp_y
      real(dp), dimension(:), allocatable :: tmp_w
      real(dp), parameter :: phi0 = (120.0d0/180.0d0)*pi
      integer :: i
      integer :: npol_tot
  
      call mem_man%alloc(points1, bem%npol1, "points1")
      call mem_man%alloc(weight1, bem%npol1, "weight1")
      call mem_man%alloc(points2, bem%npol2, "points2")
      call mem_man%alloc(weight2, bem%npol2, "weight2")
  
      call legendre_gauss_lobatto_integration(bem%npol1,points1,weight1)
      call legendre_gauss_lobatto_integration(bem%npol2,points2,weight2)
  
      !allocation of array for angles
      call mem_man%alloc(rho, bem%npol1, "rho")
      call mem_man%alloc(phi, bem%npol2, "phi")
  
      !Calculation of rho and phi
      rho = half*(points1 + one + 1.0d-06)
      phi = pi*(270.0d0 + 60.0d0*points2)/180.0d0
       
      !2D array for rho and phi
      npol_tot = bem%npol1*bem%npol2
  
      call mem_man%alloc(tmp_phi, npol_tot, "tmp_phi")
      do i = 1, npol_tot
         if (i.ge.1 .and. i.le.bem%npol1) then
            tmp_phi(i) = phi(1)
         else if (i .gt. 1*bem%npol1 .and. i .le. 2*bem%npol1) then
            tmp_phi(i) = phi(2)
         else if (i .gt. 2*bem%npol1 .and. i .le. 3*bem%npol1) then
            tmp_phi(i) = phi(3)
         else if (i .gt. 3*bem%npol1 .and. i .le. 4*bem%npol1) then
            tmp_phi(i) = phi(4)
         else if (i .gt. 4*bem%npol1 .and. i .le. 5*bem%npol1) then
            tmp_phi(i) = phi(5)
         else if (i .gt. 5*bem%npol1 .and. i .le. 6*bem%npol1) then
            tmp_phi(i) = phi(6)
         endif 
      enddo
  
      call mem_man%alloc(tmp_rho, npol_tot, "tmp_rho")
      tmp_rho(1:8)   = rho
      tmp_rho(9:16)  = rho
      tmp_rho(17:24) = rho
      tmp_rho(25:32) = rho
      tmp_rho(33:40) = rho
      tmp_rho(41:48) = rho
       
      !Creazione vettori finali per coordinate e pesi
      call mem_man%alloc(tmp_x, npol_tot*3, "tmp_x")
      call mem_man%alloc(tmp_y, npol_tot*3, "tmp_y")
  
      do i = 1, npol_tot
         tmp_x(i + 0*npol_tot) = cos(tmp_phi(i)             )*tmp_rho(i)*one/ &
                                 (abs(two*sin(tmp_phi(i))))   
         tmp_x(i + 1*npol_tot) = cos(tmp_phi(i) + phi0      )*tmp_rho(i)*one/ &
                                 (abs(two*sin(tmp_phi(i))))   
         tmp_x(i + 2*npol_tot) = cos(tmp_phi(i) + (two*phi0))*tmp_rho(i)*one/ &
                                 (abs(two*sin(tmp_phi(i))))   
         tmp_y(i + 0*npol_tot) = sin(tmp_phi(i)             )*tmp_rho(i)*one/ &
                                 (abs(two*sin(tmp_phi(i))))   
         tmp_y(i + 1*npol_tot) = sin(tmp_phi(i) + phi0      )*tmp_rho(i)*one/ &
                                 (abs(two*sin(tmp_phi(i))))   
         tmp_y(i + 2*npol_tot) = sin(tmp_phi(i) + (two*phi0))*tmp_rho(i)*one/ &
                                 (abs(two*sin(tmp_phi(i))))   
      enddo
      do i = 1, npol_tot * 3 
         int_coord_x(i) = (one - sqrt(three)*tmp_x(i) - tmp_y(i))/three   
         int_coord_y(i) = (one + sqrt(three)*tmp_x(i) - tmp_y(i))/three   
      enddo
       
      !weights
      call mem_man%alloc(tmp_w, npol_tot, "tmp_w")
      call dgemm('N', 'T',   &
                 bem%npol1,  &
                 bem%npol2,  &
                 1,          &
                 one,        &
                 weight1,    &
                 bem%npol1,  &
                 weight2,    &
                 bem%npol2,  &
                 zero,       &
                 tmp_w,      &
                 bem%npol1)
      do i = 1, npol_tot
         bem%integration_weights(i + 0*npol_tot) = tmp_w(i)*tmp_rho(i)* &
                                         ((abs(one/(two*sin(tmp_phi(i)))))**Two)
         bem%integration_weights(i + 1*npol_tot) = tmp_w(i)*tmp_rho(i)* &
                                         ((abs(one/(two*sin(tmp_phi(i)))))**Two)
         bem%integration_weights(i + 2*npol_tot) = tmp_w(i)*tmp_rho(i)* &
                                         ((abs(one/(two*sin(tmp_phi(i)))))**Two)
      enddo
      bem%integration_weights = bem%integration_weights/ &
                                sum(bem%integration_weights)
  
      !deallocation
      call mem_man%dealloc(rho, "rho")
      call mem_man%dealloc(phi, "phi")
      call mem_man%dealloc(tmp_rho, "tmp_rho")
      call mem_man%dealloc(tmp_phi, "tmp_phi")
      call mem_man%dealloc(tmp_x, "tmp_x")
      call mem_man%dealloc(tmp_y, "tmp_y")
      call mem_man%dealloc(tmp_w, "tmp_w")
      call mem_man%dealloc(points1, "points1")
      call mem_man%dealloc(weight1, "weight1")
      call mem_man%dealloc(points2, "points2")
      call mem_man%dealloc(weight2, "weight2")
  
   end subroutine tesserae_integration_polar_coordinates


   !> Subroutine for the integration based on Legendre_Gauss-Lobatto
   !!    In/Out  : npol        -- npoints
   !!    In/Out  : points      -- points
   !!    In/Out  : weights     -- associated weights
   subroutine legendre_gauss_lobatto_integration(npol,points,weights)
  
      implicit none
       
      !input/output variables
      integer, intent(in)                      :: npol
      real(dp), dimension(npol), intent(inout) :: points
      real(dp), dimension(npol), intent(inout) :: weights
   
      !internal variables
      integer :: i, j
      real(dp), dimension(:), allocatable   :: tmp
      real(dp), dimension(:), allocatable   :: first_array
      real(dp), dimension(:), allocatable   :: old_array
      real(dp), dimension(:,:), allocatable :: legendre_vandermonde_matrix
      real(dp) :: check_value 
      real(dp), parameter :: eps = 2.2204d-16
  
  
      call mem_man%alloc(tmp, npol, "tmp")
      if (npol.eq.8) then
         tmp(1) = zero
         tmp(2) = one
         tmp(3) = two
         tmp(4) = three
         tmp(5) = four
         tmp(6) = five
         tmp(7) = six
         tmp(8) = seven
      else if (npol.eq.6) then 
         tmp(1) = zero
         tmp(2) = one
         tmp(3) = two
         tmp(4) = three
         tmp(5) = four
         tmp(6) = five
      else         
         call out_%error('npol.ne.8 e npol.ne.6 in &
                         &legendre_gauss_lobatto_integration')
      endif
      !First guess: Chebyshev-Gauss-Lobatto nodes
      do i = 1, npol
         points(i) = cos((pi*tmp(i))/(npol - 1))
      enddo
  
      call mem_man%dealloc(tmp, "tmp")
      !allocation of the Legendre Vandermonde Matrix
      call mem_man%alloc(legendre_vandermonde_matrix, npol,npol, &
                         "legendre_vandermonde_matrix")
      !allocation of arrays for iterative method for construct the matrix
      call mem_man%alloc(first_array, npol, "first_array")
      call mem_man%alloc(old_array, npol, "old_array")
      old_array = two

      !iterative method
      first_array = abs(points - old_array)
      check_value = maxval(first_array)
      do while(check_value .gt. eps) 
         first_array = abs(points - old_array)
         check_value = maxval(first_array)
         old_array   = points
         do i = 1, npol
            do j = 1, npol
               if (j .eq. 1) then
                  legendre_vandermonde_matrix(i,j) = one
               elseif (j .eq. 2) then 
                  legendre_vandermonde_matrix(i,j) = points(i)
               else 
                  legendre_vandermonde_matrix(i,j) = &
                                                  ( (two*(j-1)-one)*points(i)* &
                                          legendre_vandermonde_matrix(i,j-1) - &
                        ( (j-1-one) * legendre_vandermonde_matrix(i,j-2) ) ) / &
                                                                           (j-1)
               endif  
            enddo
         enddo  
  
         !vector update
         do i = 1, npol
            points(i) = old_array(i)-( points(i) * &
                                         legendre_vandermonde_matrix(i,npol) - &
                                      legendre_vandermonde_matrix(i,npol-1) )/ &
                                      (npol*legendre_vandermonde_matrix(i,npol))
         enddo
      enddo
      !weights calculation
      weights = two/(npol*(npol-one)*(legendre_vandermonde_matrix(:,npol)**two))

      !dealloc
      call mem_man%dealloc(legendre_vandermonde_matrix, &
                           "legendre_vandermonde_matrix")
      call mem_man%dealloc(first_array, "first_array")
      call mem_man%dealloc(old_array, "old_array")
  
   end subroutine legendre_gauss_lobatto_integration


   !> Subroutine for flat integration off diagonal part
   !!    In/Out  : bem      -- bem_type
   module subroutine flat_integration_off_diagonal(bem)
  
      implicit none
       
      !input/output variables
      class(bem_type), intent(inout) :: bem
  
      !internal variables
      integer :: i, k, k_1 
      integer :: n_update
      real(dp), dimension(:,:), allocatable :: internal_matrix
      real(dp), dimension(:,:), allocatable :: tmp_coord
      real(dp), dimension(:,:), allocatable :: product_matrix
  
      !allocation of useless vectors 
      call mem_man%alloc(internal_matrix, 28,3, "internal_matrix")
      do i = 1, 28
         internal_matrix(i,1) = bem%x_tab(i)
         internal_matrix(i,2) = bem%y_tab(i)
         internal_matrix(i,3) = one - bem%x_tab(i) - bem%y_tab(i)
      enddo
      call mem_man%alloc(tmp_coord, 3,3, "tmp_coord")
      call mem_man%alloc(product_matrix, 28,3, "product_matrix")
      call mem_man%alloc(bem%integration_coord,3,bem%n_var*28,&
                        "bem%integration_coord")
  
      n_update = 0
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
         product_matrix = matmul(internal_matrix,tmp_coord)
         do k = 1, 28
            k_1 = k + n_update
            bem%integration_coord(1,k_1) = product_matrix(k,1)
            bem%integration_coord(2,k_1) = product_matrix(k,2)
            bem%integration_coord(3,k_1) = product_matrix(k,3)
         enddo
         n_update = n_update + 28
      enddo
      !printing
      if(out_%ivrb.ge.4) call out_%print_matrix("BEM Integration Coord", &
                                     bem%integration_coord,3,bem%n_var*28)
      !deallocation
      call mem_man%dealloc(internal_matrix, "internal_matrix")
      call mem_man%dealloc(tmp_coord, "tmp_coord")
      call mem_man%dealloc(product_matrix, "product_matrix")
  
   end subroutine flat_integration_off_diagonal


   !> Subroutine for selecting the refined elements
   !!    In/Out  : bem      -- bem_type
   module subroutine selection_of_refined_elements(bem)
  
      implicit none
       
      !input/output variables
      class(bem_type), intent(inout) :: bem
   
      !internal variables
      integer :: i,j,k 
      real(dp) :: radius_x
      real(dp) :: radius_y
      real(dp) :: radius_z
      real(dp), dimension(:), allocatable   :: radius
      real(dp), dimension(:), allocatable   :: norm_coord
      real(dp), dimension(:,:), allocatable :: tmp_coord
      real(dp), dimension(:,:), allocatable :: prod_coord
      real(dp), dimension(:,:), allocatable :: tmp_1
      real(dp), dimension(:,:), allocatable :: tmp_2
  
      call mem_man%alloc(tmp_coord,3,3,"tmp_coord")
      call mem_man%alloc(radius,bem%n_var,"radius")
       
      !vertices and distance between vertices and centroids coord.
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
         radius_x = sqrt((bem%coord(1,i) - tmp_coord(1,1))**2 + &
                         (bem%coord(2,i) - tmp_coord(1,2))**2 + &
                         (bem%coord(3,i) - tmp_coord(1,3))**2)
         radius_y = sqrt((bem%coord(1,i) - tmp_coord(2,1))**2 + &
                         (bem%coord(2,i) - tmp_coord(2,2))**2 + &
                         (bem%coord(3,i) - tmp_coord(2,3))**2)
         radius_z = sqrt((bem%coord(1,i) - tmp_coord(3,1))**2 + &
                         (bem%coord(2,i) - tmp_coord(3,2))**2 + &
                         (bem%coord(3,i) - tmp_coord(3,3))**2)
         !max value for the minimal radius containing the tesserae
         if(radius_x.ge.radius_y) then
            if(radius_x.ge.radius_z) then
               radius(i) = radius_x        
            else
               radius(i) = radius_z        
            endif
         else
            if(radius_y.ge.radius_z) then
               radius(i) = radius_y
            else
               radius(i) = radius_z
            endif
         endif
      enddo
      !distance array between positions p1 and p2 (centroids)
      call mem_man%alloc(norm_coord, bem%n_var, "norm_coord")
      call mem_man%alloc(prod_coord, bem%n_var,bem%n_var, "prod_coord")
      call mem_man%alloc(tmp_1, bem%n_var,bem%n_var, "tmp_1")
      call mem_man%alloc(tmp_2, bem%n_var,bem%n_var, "tmp_2")
      call mem_man%alloc(bem%selected_elements, bem%n_var,bem%n_var, &
                         "bem%selected_elements")

      !norm of the coordinates
      do i = 1, bem%n_var
         norm_coord(i) = bem%coord(1,i)**2 + & 
                         bem%coord(2,i)**2 + & 
                         bem%coord(3,i)**2 
      enddo
      call dgemm('T', 'N',   &
                 bem%n_var,  &
                 bem%n_var,  &
                 3,          &
                 two,        &
                 bem%coord,  &
                 3,          &
                 bem%coord,  &
                 3,          &
                 zero,       &
                 prod_coord, &
                 bem%n_var)
      do i = 1, bem%n_var
         do j = 1, bem%n_var
            tmp_1(i,j) = norm_coord(i) + norm_coord(j) - prod_coord(i,j)
            if (tmp_1(i,j).lt.1.0e-10) tmp_1(i,j) = zero
            tmp_2(i,j) = radius(j)
         enddo
      enddo
      !selected elements = sqrt(tmp_1) - tmp_2
      bem%selected_elements = sqrt(tmp_1) - tmp_2
      !printing
      if (out_%ivrb.eq.4) &
          call out_%print_matrix("selected elements",bem%selected_elements, &
                                                     bem%n_var,bem%n_var)
      !deallocate
      call mem_man%dealloc(radius, "radius")
      call mem_man%dealloc(tmp_coord, "tmp_coord")
      call mem_man%dealloc(prod_coord, "prod_coord")
      call mem_man%dealloc(norm_coord, "norm_coord")
      call mem_man%dealloc(tmp_1, "tmp_1")
      call mem_man%dealloc(tmp_2, "tmp_2")
  
   end subroutine selection_of_refined_elements


   !> Subroutine for refining the diagonal elements
   !!    In/Out  : bem      -- bem_type
   module subroutine refinement_diagonal_elements(bem)
  
      implicit none
       
      !input/output variables
      class(bem_type), intent(inout) :: bem
  
      !internal variables
      integer :: i, j, k
      real(dp) :: vec_x
      real(dp) :: vec_y
      real(dp) :: vec_z
      real(dp) :: green_sum_normal    
      real(dp) :: green_sum_tangent_1 
      real(dp) :: green_sum_tangent_2 
      real(dp) :: max_distance_polar_coord_1e_4
  
      call mem_man%alloc(bem%green_diagonal_accurate, 3,bem%n_var, &
                         "bem%green_diagonal_accurate")
      max_distance_polar_coord_1e_4 = maxval(bem%distance_polar_coord)*1.0d-4
  
      k = 0
      Do i = 1, bem%n_var
         green_sum_normal    = zero
         green_sum_tangent_1 = zero
         green_sum_tangent_2 = zero
         do j = 1, bem%npol1*bem%npol2*3
            k = k + 1
            vec_x = bem%coord(1,i) - bem%polar_coord(1,k)
            vec_y = bem%coord(2,i) - bem%polar_coord(2,k)
            vec_z = bem%coord(3,i) - bem%polar_coord(3,k)
            green_sum_normal = green_sum_normal - &
                               bem%integration_weights(j)*bem%area(i) * &
                               (vec_x*bem%normal(1,i) + &
                                vec_y*bem%normal(2,i) + &
                                vec_z*bem%normal(3,i)) / &
                                bem%distance_polar_coord(k)**3
            if(bem%distance_polar_coord(k).lt.1.0E-04) &
                     bem%distance_polar_coord(k) = max_distance_polar_coord_1e_4
            green_sum_tangent_1 = green_sum_tangent_1 - &
                                  bem%integration_weights(j)*bem%area(i)* &
                                  ( vec_x*bem%tangent_1(1,i)   + &
                                    vec_y*bem%tangent_1(2,i)   + &
                                    vec_z*bem%tangent_1(3,i) ) / &
                                   bem%distance_polar_coord(k)**3
            green_sum_tangent_2 = green_sum_tangent_2 -        &
                                  bem%integration_weights(j)*bem%area(i)* &
                                  ( vec_x*bem%tangent_2(1,i)  + &
                                    vec_y*bem%tangent_2(2,i)  + &
                                    vec_z*bem%tangent_2(3,i) ) / &
                                   bem%distance_polar_coord(k)**3
         enddo
         bem%green_diagonal_accurate(1,i) = bem%normal   (1,i)  * &
                                         green_sum_normal    + &
                                         bem%tangent_1(1,i)  * &
                                         green_sum_tangent_1 + & 
                                         bem%tangent_2(1,i)  * &
                                         green_sum_tangent_2  
         bem%green_diagonal_accurate(2,i) = bem%normal   (2,i)  * &
                                         green_sum_normal    + &
                                         bem%tangent_1(2,i)  * &
                                         green_sum_tangent_1 + & 
                                         bem%tangent_2(2,i)  * &
                                         green_sum_tangent_2 
         bem%green_diagonal_accurate(3,i) = bem%normal   (3,i)  * &
                                         green_sum_normal    + &
                                         bem%tangent_1(3,i)  * &
                                         green_sum_tangent_1 + & 
                                         bem%tangent_2(3,i)  * &
                                         green_sum_tangent_2  
      enddo
      !printing
      if(out_%ivrb.ge.4) call out_%print_matrix("green_diagonal_accurate", &
                              transpose(bem%green_diagonal_accurate),bem%n_var,3)
  
   end subroutine refinement_diagonal_elements


   !> Subroutine for refining the off diagonal elements
   !!    Input   : bem                  -- bem_type
   !!    Input   : i_row                -- index of the row
   !!    Input   : j_column             -- index of the column
   !!    Output  : off_diagonal_element -- i,j element
   module function refinement_off_diagonal_elements(bem,i_row,j_column) &
                   result(off_diagonal_element)
  
      implicit none
       
      !input/output variables
      class(bem_type), intent(in) :: bem
      integer, intent(in) :: i_row
      integer, intent(in) :: j_column
      real(dp) :: off_diagonal_element
  
      !internal variables
      integer  :: i, k
      real(dp) :: off_diag_x
      real(dp) :: off_diag_y
      real(dp) :: off_diag_z
      real(dp) :: tmp_x
      real(dp) :: tmp_y
      real(dp) :: tmp_z
      real(dp) :: tmp_norm
  
      off_diag_x = zero
      off_diag_y = zero
      off_diag_z = zero
  
      k = (j_column-1)*28
      do i = 1, 28
         k = k + 1
         tmp_x = bem%coord(1,i_row) - bem%integration_coord(1,k)
         tmp_y = bem%coord(2,i_row) - bem%integration_coord(2,k)
         tmp_z = bem%coord(3,i_row) - bem%integration_coord(3,k)
         tmp_norm = sqrt(tmp_x**2 + tmp_y**2 + tmp_z**2)
         off_diag_x = off_diag_x - bem%w_tab(i)*bem%area(j_column)*tmp_x/ &
                                   (tmp_norm**3)
         off_diag_y = off_diag_y - bem%w_tab(i)*bem%area(j_column)*tmp_y/ &
                                   (tmp_norm**3)
         off_diag_z = off_diag_z - bem%w_tab(i)*bem%area(j_column)*tmp_z/ &
                                   (tmp_norm**3)
       enddo
       off_diagonal_element = off_diag_x*bem%normal(1,i_row) + & 
                              off_diag_y*bem%normal(2,i_row) + & 
                              off_diag_z*bem%normal(3,i_row)               
  
   end function refinement_off_diagonal_elements


   !> Subroutine for the final refinement of the Green function
   !!    In/Out  : bem                  -- bem_type
   !!    Input   : i_row                -- index of the row
   !!    Input   : j_column             -- index of the column
   !!    Output  : off_diagonal_element -- i,j element
   module subroutine final_refinement_green_function(bem, matrix_constant)
  
      implicit none
       
      !input/output variables
      class(bem_type) :: bem
      real(dp), dimension(bem%n_var,bem%n_var), intent(inout) :: matrix_constant
  
      !internal variables
      integer :: i 
      real(dp), dimension(:,:), allocatable :: matrix_tmp
  
      call mem_man%alloc(matrix_tmp, bem%n_var, bem%n_var, "matrix_tmp")
  
      matrix_tmp = matrix_constant
      do i = 1, bem%n_var
         matrix_tmp(:,i) = matrix_tmp(:,i) / bem%area(i)
      enddo
      do i = 1, bem%n_var
         matrix_tmp(i,:) = matrix_tmp(i,:) * bem%area(i)
      enddo
      do i = 1, bem%n_var
          
         bem%green_diagonal_accurate(1,i) = bem%normal(1,i) * &
                                            (-2*pi + sum(matrix_tmp(:,i))) + &
                                            bem%green_diagonal_accurate(1,i)
         bem%green_diagonal_accurate(2,i) = bem%normal(2,i) * &
                                            (-2*pi + sum(matrix_tmp(:,i))) + &
                                            bem%green_diagonal_accurate(2,i)
         bem%green_diagonal_accurate(3,i) = bem%normal(3,i) * &
                                            (-2*pi + sum(matrix_tmp(:,i))) + &
                                            bem%green_diagonal_accurate(3,i)
         matrix_constant(i,i) = - ( bem%green_diagonal_accurate(1,i) * &
                                    bem%normal(1,i)                  + &
                                    bem%green_diagonal_accurate(2,i) * &
                                    bem%normal(2,i)                  + &
                                    bem%green_diagonal_accurate(3,i) * &
                                    bem%normal(3,i) )
      enddo
  
      call mem_man%dealloc(matrix_tmp, "matrix_tmp")
  
   end subroutine final_refinement_green_function

end submodule accurate_green_function
