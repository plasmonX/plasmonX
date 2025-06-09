!> Grid module   
!!
!! This module contains the subroutines for the grid
!!
!! Date         : 2024
!!
module grid_module

   use parameters_module
   use control_analysis_module
   use target_module
   use output_analysis_module
   use field_module

   implicit none 

   public grid

   type :: grid_type
      integer :: nx
      integer :: ny
      integer :: nz
      !number of points for 2D
      integer :: n_2d_points(2)
      !offset of the grid X_max + offset
      real(dp) :: offset
      !min_coord
      real(dp), dimension(3) :: min_coord
      !max_coord
      real(dp), dimension(3) :: max_coord  
      !step coord
      real(dp), dimension(3) :: step
      !quantities on grid
      real(dp), dimension(:,:), allocatable      :: field_2d
      real(dp), dimension(:,:,:), allocatable    :: field_3d
      complex(dp), dimension(:,:,:), allocatable :: density_3d
      complex(dp), dimension(:,:), allocatable   :: density_2d
      !grid defined by the used
      logical :: user_defined = .false.
      !format = plt or cube
      character(len=4) :: format_ 

      contains

      procedure :: create_cubic_grid
      procedure :: create_plane_grid
      procedure :: print_grid_3d_info
      procedure :: calculate_field_2d
      procedure :: calculate_field_3d
      procedure :: calculate_density_2d
      procedure :: calculate_density_3d
      procedure :: get_coordinate_index

   end type grid_type
   
   type (grid_type) :: grid

contains

   !> Subroutine for creating the cubic grid
   !!    In/Out : grid        -- grid
   subroutine create_cubic_grid(grid)
       
      implicit none
  
      !input/output variables
      class(grid_type), intent(inout) :: grid
  
      if(.not.grid%user_defined) then
         grid%min_coord(1) = minval(target_%coord(1,:)) - grid%offset
         grid%min_coord(2) = minval(target_%coord(2,:)) - grid%offset
         grid%min_coord(3) = minval(target_%coord(3,:)) - grid%offset
  
         grid%max_coord(1) = maxval(target_%coord(1,:)) + grid%offset
         grid%max_coord(2) = maxval(target_%coord(2,:)) + grid%offset
         grid%max_coord(3) = maxval(target_%coord(3,:)) + grid%offset
  
      endif
      grid%step(1) = (grid%max_coord(1) - grid%min_coord(1))/(grid%nx-1)
      grid%step(2) = (grid%max_coord(2) - grid%min_coord(2))/(grid%ny-1)
      grid%step(3) = (grid%max_coord(3) - grid%min_coord(3))/(grid%nz-1)
  
      !print info
      call grid%print_grid_3d_info()
  
   end subroutine create_cubic_grid


   !> Subroutine for creating the grid for the plane
   !!    In/Out : grid        -- grid
   subroutine create_plane_grid(grid)
      
      implicit none
 
      !input/output variables
      class(grid_type), intent(inout) :: grid
 
      !xy plane
      if(control_analysis%plane.eq.'xy') then
         if(.not.grid%user_defined) then
            grid%min_coord(1) = minval(target_%coord(1,:)) - grid%offset
            grid%min_coord(2) = minval(target_%coord(2,:)) - grid%offset
            grid%min_coord(3) = zero
 
            grid%max_coord(1) = maxval(target_%coord(1,:)) + grid%offset
            grid%max_coord(2) = maxval(target_%coord(2,:)) + grid%offset
            grid%max_coord(3) = zero
         endif
         grid%step(1) = (grid%max_coord(1) - grid%min_coord(1))/(grid%nx-1)
         grid%step(2) = (grid%max_coord(2) - grid%min_coord(2))/(grid%ny-1)
         grid%step(3) = zero
         grid%nz = 0
         grid%n_2d_points(1) = grid%nx
         grid%n_2d_points(2) = grid%ny
      !xz plane
      else if(control_analysis%plane.eq.'xz') then
         if(.not.grid%user_defined) then
            grid%min_coord(1) = minval(target_%coord(1,:)) - grid%offset
            grid%min_coord(2) = zero
            grid%min_coord(3) = minval(target_%coord(3,:)) - grid%offset
            grid%max_coord(1) = maxval(target_%coord(1,:)) + grid%offset
            grid%max_coord(2) = zero
            grid%max_coord(3) = maxval(target_%coord(3,:)) + grid%offset
         endif
         grid%step(1) = (grid%max_coord(1) - grid%min_coord(1))/(grid%nx-1)
         grid%step(2) = zero
         grid%step(3) = (grid%max_coord(3) - grid%min_coord(3))/(grid%nz-1)
         grid%ny = 0
         grid%n_2d_points(1) = grid%nx
         grid%n_2d_points(2) = grid%nz
      !yz plane
      else if(control_analysis%plane.eq.'yz') then
         if(.not.grid%user_defined) then
            grid%min_coord(1) = zero
            grid%min_coord(2) = minval(target_%coord(2,:)) - grid%offset
            grid%min_coord(3) = minval(target_%coord(3,:)) - grid%offset
            grid%max_coord(1) = zero
            grid%max_coord(2) = maxval(target_%coord(2,:)) + grid%offset
            grid%max_coord(3) = maxval(target_%coord(3,:)) + grid%offset
         endif
         grid%step(1) = zero
         grid%step(2) = (grid%max_coord(2) - grid%min_coord(2))/(grid%ny-1)
         grid%step(3) = (grid%max_coord(3) - grid%min_coord(3))/(grid%nz-1)
         grid%nx = 0
         grid%n_2d_points(1) = grid%ny
         grid%n_2d_points(2) = grid%nz
      endif
 
      !printing info
      call grid%print_grid_3d_info()
 
   end subroutine create_plane_grid


   !> Subroutine for printing the info of the grid
   !!    Input  : grid        -- grid
   subroutine print_grid_3d_info(grid)
       
      implicit none
  
      !input/output variables
      class(grid_type), intent(in) :: grid
       
      !internal variables
      character(len=19) :: format_1 = "(10x,3(a,f12.4,2x))"
      character(len=24) :: format_2 = "(10x,2(a,f10.4,2x),a,i7)"
  
      write(out_analysis%iunit,'(27x,a)') 'Grid Information (Angstrom)'
      write(out_analysis%iunit,out_analysis%sticks) 
  
      write(out_analysis%iunit,format_1) 'StepX: ', grid%step(1), &
                                         'StepY: ', grid%step(2), &
                                         'StepZ: ', grid%step(3)
      write(out_analysis%iunit,out_analysis%sticks) 
  
      write(out_analysis%iunit,format_2) 'XminGrid:', grid%min_coord(1), &
                                         'XmaxGrid:', grid%max_coord(1), &
                                         'X NPoints:', grid%nx
      write(out_analysis%iunit,out_analysis%sticks)           
      write(out_analysis%iunit,format_2) 'YminGrid:', grid%min_coord(2), &
                                         'YmaxGrid:', grid%max_coord(2), &
                                         'Y NPoints:', grid%ny
      write(out_analysis%iunit,out_analysis%sticks) 
      write(out_analysis%iunit,format_2) 'ZminGrid:', grid%min_coord(3), &
                                         'ZmaxGrid:', grid%max_coord(3), & 
                                         'Z NPoints:', grid%nz
      write(out_analysis%iunit,out_analysis%sticks) 
      flush(out_analysis%iunit)
       
   end subroutine print_grid_3d_info     


   !> Subroutine for calculating the induced field in a 3D grid
   !! This calculates |E|/|E_0|
   !!    In/Out : grid        -- grid
   !!    Input  : i_pol       -- polarization index (x,y,z)
   !!    Input  : variables_w -- dynamic variables producing the field
   subroutine calculate_field_3d(grid, i_pol, variables_w)
  
      implicit none
  
      !input/output variables
      class(grid_type), intent(inout) :: grid
      integer, intent(in)             :: i_pol
      complex(dp), dimension(target_%n_var), intent(in) :: variables_w
  
      !internal variables
      integer :: k, j, i
      real(dp) :: norm_field
      real(dp), dimension(3)    :: point_coord
      complex(dp), dimension(3) :: efield
  
      integer :: idx, n_total
  
      n_total = grid%nx * grid%ny * grid%nz
      !$omp parallel do private(idx,i,j,k,point_coord,efield,norm_field) &
      !$omp schedule(dynamic)
      do idx = 0, n_total - 1
         i = mod(idx, grid%nx) + 1
         j = mod((idx / grid%nx), grid%ny) + 1
         k = idx / (grid%nx * grid%ny) + 1
      
         point_coord(1) = grid%min_coord(1) + (i-1)*grid%step(1)
         point_coord(2) = grid%min_coord(2) + (j-1)*grid%step(2)
         point_coord(3) = grid%min_coord(3) + (k-1)*grid%step(3)
  
         efield = target_%calculate_induced_field_at_point(i_pol,       &
                                                           point_coord, &
                                                           variables_w)
  
         norm_field = dsqrt(dble(efield(1))**2 + dimag(efield(1))**2 + &
                            dble(efield(2))**2 + dimag(efield(2))**2 + &
                            dble(efield(3))**2 + dimag(efield(3))**2)
  
         grid%field_3d(i,j,k) = norm_field/field%e_0
  
      enddo
      !$omp end parallel do
       
   end subroutine calculate_field_3d


   !> Subroutine for calculating the induced density in a 3D grid
   !!    In/Out : grid        -- grid
   !!    Input  : variables_w -- dynamic variables producing the field
   subroutine calculate_density_3d(grid, variables_w)
  
      implicit none
  
      !input/output variables
      class(grid_type), intent(inout) :: grid
      complex(dp), dimension(target_%n_var), intent(in) :: variables_w
  
      !internal variables
      integer :: k, j, i
      real(dp), dimension(3) :: point_coord
      complex(dp) :: density
      integer :: idx, n_total
  
      n_total = grid%nx * grid%ny * grid%nz
      !$omp parallel do private(idx,i,j,k,point_coord,density) schedule(dynamic)
      do idx = 0, n_total - 1
         i = mod(idx, grid%nx) + 1
         j = mod((idx / grid%nx), grid%ny) + 1
         k = idx / (grid%nx * grid%ny) + 1
         point_coord(1) = grid%min_coord(1) + (i-1)*grid%step(1)
         point_coord(2) = grid%min_coord(2) + (j-1)*grid%step(2)
         point_coord(3) = grid%min_coord(3) + (k-1)*grid%step(3)
  
         density = target_%calculate_density_at_point(point_coord, variables_w)
         grid%density_3d(i,j,k) = density
      enddo
      !$omp end parallel do
       
   end subroutine calculate_density_3d


   !> Subroutine for calculating the induced field in a 2D grid
   !! This calculates |E|/|E_0|
   !!    In/Out : grid        -- grid
   !!    Input  : i_pol       -- polarization index (x,y,z)
   !!    Input  : variables_w -- dynamic variables producing the field
   !!    Input  : third_coord -- the third coordinate of the plane
   subroutine calculate_field_2d(grid, i_pol, variables_w, third_coord)
  
      implicit none
  
      !input/output variables
      class(grid_type), intent(inout) :: grid
      integer, intent(in)             :: i_pol
      real(dp), intent(in)            :: third_coord 
      complex(dp), dimension(target_%n_var), intent(in) :: variables_w
  
      !internal variables
      integer :: j, i
      real(dp) :: norm_field
      real(dp), dimension(3)    :: point_coord
      complex(dp), dimension(3) :: efield
  
      !xy plane
      if(control_analysis%plane.eq.'xy') then
         !$omp parallel do private(i,j,point_coord,efield,norm_field) &
         !$omp collapse(2) schedule(dynamic)
         do j = 1, grid%ny  !y
            do i = 1, grid%nx !x
               point_coord(1) = grid%min_coord(1) + (i-1)*grid%step(1)
               point_coord(2) = grid%min_coord(2) + (j-1)*grid%step(2)
               point_coord(3) = third_coord
               efield = target_%calculate_induced_field_at_point(i_pol, &
                                                                 point_coord, &
                                                                 variables_w)
               norm_field = dsqrt(dble(efield(1))**2 + dimag(efield(1))**2 + &
                                  dble(efield(2))**2 + dimag(efield(2))**2 + &
                                  dble(efield(3))**2 + dimag(efield(3))**2)
               grid%field_2d(i,j) = norm_field/field%e_0
  
            enddo
         enddo
         !$omp end parallel do
      !xz plane
      else if(control_analysis%plane.eq.'xz') then
         !$omp parallel do private(i,j,point_coord,efield, norm_field) &
         !$omp collapse(2) schedule(dynamic)
         do j = 1, grid%nz  !z
            do i = 1, grid%nx !x
               point_coord(1) = grid%min_coord(1) + (i-1)*grid%step(1)
               point_coord(2) = third_coord
               point_coord(3) = grid%min_coord(3) + (j-1)*grid%step(3)
               efield = target_%calculate_induced_field_at_point(i_pol, &
                                                                 point_coord, &
                                                                 variables_w)
               norm_field = dsqrt(dble(efield(1))**2 + dimag(efield(1))**2 + &
                                  dble(efield(2))**2 + dimag(efield(2))**2 + &
                                  dble(efield(3))**2 + dimag(efield(3))**2)
               grid%field_2d(i,j) = norm_field/field%e_0
            enddo
         enddo
         !$omp end parallel do
      !yz plane
      else if(control_analysis%plane.eq.'yz') then
         !$omp parallel do private(i,j,point_coord,efield,norm_field) &
         !$omp collapse(2) schedule(dynamic)
         do j = 1, grid%nz  !z
            do i = 1, grid%ny !y
  
               point_coord(1) = third_coord
               point_coord(2) = grid%min_coord(2) + (i-1)*grid%step(2)
               point_coord(3) = grid%min_coord(3) + (j-1)*grid%step(3)
               efield = target_%calculate_induced_field_at_point(i_pol, &
                                                                 point_coord, &
                                                                 variables_w)
               norm_field = dsqrt(dble(efield(1))**2 + dimag(efield(1))**2 + &
                                  dble(efield(2))**2 + dimag(efield(2))**2 + &
                                  dble(efield(3))**2 + dimag(efield(3))**2)
               grid%field_2d(i,j) = norm_field/field%e_0
            enddo
         enddo
         !$omp end parallel do
      endif
  
       
   end subroutine calculate_field_2d


   !> Subroutine for calculating the induced density in a 2D grid
   !!    In/Out : grid        -- grid
   !!    Input  : variables_w -- dynamic variables producing the field
   !!    Input  : third_coord -- the third coordinate of the plane
   subroutine calculate_density_2d(grid, variables_w, third_coord)
  
      implicit none
  
      !input/output variables
      class(grid_type), intent(inout) :: grid
      real(dp), intent(in) :: third_coord ! third coordinate -- plane
      complex(dp), dimension(target_%n_var), intent(in) :: variables_w
  
      !internal variables
      integer :: j, i
      real(dp), dimension(3)    :: point_coord
      complex(dp) :: density
  
      !xy plane
      if(control_analysis%plane.eq.'xy') then
         !$omp parallel do private(i,j,point_coord,density), collapse(2)
         do j = 1, grid%ny  !y
            do i = 1, grid%nx !x
               point_coord(1) = grid%min_coord(1) + (i-1)*grid%step(1)
               point_coord(2) = grid%min_coord(2) + (j-1)*grid%step(2)
               point_coord(3) = third_coord
               density = target_%calculate_density_at_point(point_coord, &
                                                            variables_w)
               grid%density_2d(i,j) = density
            enddo
         enddo
         !$omp end parallel do
      !xz plane
      else if(control_analysis%plane.eq.'xz') then
         !$omp parallel do private(i,j,point_coord,density), collapse(2)
         do j = 1, grid%nz  !z
            do i = 1, grid%nx !x
               point_coord(1) = grid%min_coord(1) + (i-1)*grid%step(1)
               point_coord(2) = third_coord
               point_coord(3) = grid%min_coord(3) + (j-1)*grid%step(3)
               density = target_%calculate_density_at_point(point_coord, &
                                                            variables_w)
               grid%density_2d(i,j) = density
            enddo
         enddo
         !$omp end parallel do
      !yz plane
      else if(control_analysis%plane.eq.'yz') then
         !$omp parallel do private(i,j,point_coord,density), collapse(2)
         do j = 1, grid%nz  !z
            do i = 1, grid%ny !y
               point_coord(1) = third_coord
               point_coord(2) = grid%min_coord(2) + (i-1)*grid%step(2)
               point_coord(3) = grid%min_coord(3) + (j-1)*grid%step(3)
               density = target_%calculate_density_at_point(point_coord, &
                                                            variables_w)
               grid%density_2d(i,j) = density
            enddo
         enddo
         !$omp end parallel do
      endif
  
   end subroutine calculate_density_2d


   !> Function to getting the coordinate given the index
   !!    Input  : grid        -- grid
   !!    Input  : index_      -- index of the coordinate
   !!    Output : coord       -- x,y,z of the requested index
   function get_coordinate_index(grid, index_) result (coord)
  
      implicit none
  
      !input/output variables
      class(grid_type), intent(in) :: grid
      integer, dimension(3), intent(in) :: index_
      real(dp), dimension(3) :: coord
       
      !internal variables
      real(dp), dimension(3) :: dummy_coord 
      integer :: i,j,k
  
      do k = 1, grid%nz  !z
         do j = 1, grid%ny  !y
            do i = 1, grid%nx !x
               dummy_coord(1) = grid%min_coord(1) + (i-1)*grid%step(1)
               dummy_coord(2) = grid%min_coord(2) + (j-1)*grid%step(2)
               dummy_coord(3) = grid%min_coord(3) + (k-1)*grid%step(3)
               if (index_(1).eq.i .and. index_(2).eq.j.and.index_(3).eq.k) then
                  coord = dummy_coord
                  return
               endif
            enddo
         enddo
      enddo
       
   end function get_coordinate_index


end module grid_module
