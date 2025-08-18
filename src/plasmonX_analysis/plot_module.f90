!> Plot module   
!!
!! This module contains the subroutines for plotting
!!
!! Date         : 2025
!!
module plot_module

   use parameters_module
   use input_analysis_module
   use control_analysis_module
   use target_module
   use grid_module

   implicit none 

contains

   !> Subroutine for creating xyz file
   subroutine plot_xyz()
  
      implicit none
  
      !internal variables
      character(len=1000) :: xyz_file
      integer :: unit_xyz
      integer :: iost
      integer :: i
      character(len=100) :: str_natoms
      character(len=20) :: format_1 = "(a2,2x,3(f25.16,1x))"
      
      xyz_file = trim(out_analysis%root_folder)// &
                 trim(out_analysis%root_filename)//".xyz"
      write(out_analysis%iunit,'(1x,a)')"The requested xyz file is located in: "
      write(out_analysis%iunit,'(1x,a)') trim(xyz_file)
      write(out_analysis%iunit,out_analysis%sticks)
       
      unit_xyz = 13
      open(unit=unit_xyz,file=trim(xyz_file),status="unknown",iostat=iost)
         if (iost.ne.0) call out_analysis%error('Error in opening(w) '//&
                                                trim(xyz_file))
         if(target_%name_.ne.'bem') then
            write(str_natoms,'(i10)') target_%n_atoms
            !n_atoms + blank line
            write(unit_xyz,'(a/)') adjustl(trim(str_natoms))
            !xyz atoms
            do i = 1, target_%n_atoms
               write(unit_xyz,format_1) adjustl(target_%atom_name(i)), &
                                        target_%coord(1,i), &
                                        target_%coord(2,i), &
                                        target_%coord(3,i)
            enddo
         else
            write(str_natoms,'(i10)') target_%n_var
            !n_var + blank line
            write(unit_xyz,'(a/)') adjustl(trim(str_natoms))
            !xyz atoms
            do i = 1, target_%n_var
               write(unit_xyz,format_1) "X ", &
                                        target_%coord(1,i), &
                                        target_%coord(2,i), &
                                        target_%coord(3,i)
            enddo
         endif
      close(unit_xyz)
  
   end subroutine plot_xyz


   !> Subroutine for creating pqr file
   subroutine plot_pqr()
       
      implicit none
       
      !internal variables
      complex(dp), dimension(:,:), allocatable :: variables_w
      character(len=1000) :: pqr_file
      integer :: i, j
      integer :: ierr
      character(len=12)  :: string_freq
  
      if(target_%forcefield.eq.'fqfmu') &
         call out_analysis%error("PQR not implemented for FQFMu")
       
      !if the field is not dynamic then error 
      if(field%static) &
         call out_analysis%error("PQR Analysis only for Field Dynamic")
       
      allocate(variables_w(target_%n_var,3),stat=ierr)
      if(ierr.gt.0) &
         call out_%error('Not enough space for dynamic field variables_w')
  
      do i = 1, field%n_freq
         write(string_freq,'(f12.10)') field%freq(i)
         string_freq = adjustl(string_freq)

         !read variables
         call read_variables_w_freq_file(string_freq, variables_w)

         do j = 1, field%n_polarization
            !Real charges
            pqr_file = trim(out_analysis%root_folder) //         &
                       trim(out_analysis%root_filename)// "-" // & 
                       field%polarization_name(j)// "-"//       &
                       trim(string_freq)//'-Re.pqr'
  
            call write_pqr_file(pqr_file, dble(variables_w), &
                                field%index_rhs_polarization(j))
  
            !Imaginary charges
            pqr_file = trim(out_analysis%root_folder) //         &
                       trim(out_analysis%root_filename)// "-" // & 
                       field%polarization_name(j)// "-"//       &
                       trim(string_freq)//'-Im.pqr'
            call write_pqr_file(pqr_file, aimag(variables_w), &
                                field%index_rhs_polarization(j))
         enddo
      enddo
  
   end subroutine plot_pqr


   !> Subroutine for creating field or density cube file
   subroutine plot_field_or_density()
       
      !plot field or density
      implicit none
       
      !internal variables
      complex(dp), dimension(:,:), allocatable :: variables_w
      integer :: i, j
      integer :: ierr
      character(len=12)  :: string_freq
      character(len=1000) :: plt_file
      character(len=1000) :: cube_file
      character(len=1000) :: root_3d
      real(dp) :: integral_field_volume, maximum_field_volume
      real(dp), dimension(2,3) :: max_quantity = zero ! maximum
      integer,  dimension(2,3,3) :: i_max ! index_of_the_max
       
      i_max = 0
      !if the field is not dynamic then error 
      if(field%static) & 
         call out_analysis%error("Field Plot only for Dynamic Fields")
       
      allocate(variables_w(target_%n_var,3),stat=ierr)
      if(ierr.gt.0) &
         call out_%error('Not enough space for dynamic field variables_w')
  
      !before starting the frequency cycle create the cubic grid or the plane
      if (control_analysis%plane_requested) then 
         call grid%create_plane_grid()
         call out_analysis%create_plane_folder()
      else !for volume or standard option
         call grid%create_cubic_grid()
      endif
  
      !cycle on external frequencies
      do i = 1, field%n_freq
         write(string_freq,'(f12.10)') field%freq(i)
         string_freq = adjustl(string_freq)
         !recover from file
         call read_variables_w_freq_file(string_freq, variables_w)

         !cycle on polarizations
         do j = 1, field%n_polarization
            !plane
            if(control_analysis%plane_requested) then 
               if(control_analysis%what.eq.'field') then
                  call plot_field_on_planes(j, &
                               variables_w(:,field%index_rhs_polarization(j)), &
                               string_freq)
               else if(control_analysis%what.eq.'density') then
                  if (target_%name_.eq.'bem') &
                     call out_analysis%error('Density plot not allowed for BEM')
                  call plot_density_on_planes(j, &
                               variables_w(:,field%index_rhs_polarization(j)), &
                               string_freq)
               endif
            !cube or plt plot
            else 
               if(control_analysis%what .eq. 'field') then
                  !calculation field on grid
                  allocate(grid%field_3d(grid%nx , grid%ny , grid%nz))
                  call grid%calculate_field_3d(field%index_rhs_polarization(j), &
                                               variables_w(:, &
                                               field%index_rhs_polarization(j)))
                  max_quantity(1,j) = maxval(grid%field_3d)
                  i_max(1,j,:)      = maxloc(grid%field_3d)
                  !plot based on plt or cube option
                  if(trim(grid%format_).eq.'plt') then
                     plt_file = trim(out_analysis%root_folder) //         &
                                trim(out_analysis%root_filename)// "-" // & 
                                field%polarization_name(j)// "-"//        &
                                trim(string_freq)//'.plt'
                     call write_plt_file(plt_file, grid%field_3d)
                  else if(trim(grid%format_).eq.'cube') then
                     cube_file = trim(out_analysis%root_folder) //         &
                                 trim(out_analysis%root_filename)// "-" // & 
                                 field%polarization_name(j)// "-"//        &
                                 trim(string_freq)//'.cube'
                     call write_cube_file(cube_file, grid%field_3d)
                  endif
                  !here volume calculation if requested
                  if(control_analysis%volume_requested) then 
                     integral_field_volume = sum(grid%field_3d**2)* &
                                                 grid%step(1) * &
                                                 grid%step(2) * &
                                                 grid%step(3) 
                     maximum_field_volume = max_quantity(1,j)**2
                     integral_field_volume = integral_field_volume
                     call print_effective_volume_area('volume',                &
                                                   field%polarization_name(j), &
                                                        integral_field_volume, &
                                                        maximum_field_volume)
                  endif
                  deallocate(grid%field_3d)
               else if(control_analysis%what.eq.'density') then
                  if (target_%name_.eq.'bem') &
                     call out_analysis%error('Density plot not allowed for BEM')
                  !calculation density on grid
                  allocate(grid%density_3d(grid%nx, grid%ny, grid%nz))
                  !if q and mu separated, allocate the proper arrays
                  if (control_analysis%separate_q_mu) then
                     allocate(grid%density_3d_q(grid%nx, grid%ny, grid%nz))
                     allocate(grid%density_3d_mu(grid%nx, grid%ny, grid%nz))
                  endif
                  call grid%calculate_density_3d(variables_w(:, &
                                               field%index_rhs_polarization(j)))
                  max_quantity(1,j) = maxval(dble(grid%density_3d))
                  i_max(1,j,:)      = maxloc(dble(grid%density_3d))
                  max_quantity(2,j) = maxval(aimag(grid%density_3d))
                  i_max(2,j,:)      = maxloc(aimag(grid%density_3d))
                  !common root for output files
                  root_3d = trim(out_analysis%root_folder) //         &
                            trim(out_analysis%root_filename)// "-" // & 
                            field%polarization_name(j)// "-"//        &
                            trim(string_freq)
                  !plot based on plt or cube option
                  if(trim(grid%format_).eq.'plt') then
                     !write total density plt files
                     plt_file = trim(root_3d) // '-densityRe.plt'
                     call write_plt_file(plt_file, dble(grid%density_3d))
                     plt_file = trim(root_3d) // '-densityIm.plt'
                     call write_plt_file(plt_file, aimag(grid%density_3d))
                     if (control_analysis%separate_q_mu) then
                        !write FQ density plt files
                        plt_file = trim(root_3d) // '-densityRe_wFQ.plt'
                        call write_plt_file(plt_file, dble(grid%density_3d_q))
                        plt_file = trim(root_3d) // '-densityIm_wFQ.plt'
                        call write_plt_file(plt_file, aimag(grid%density_3d_q))
                        !write FMu density plt files
                        plt_file = trim(root_3d) // '-densityRe_wFMu.plt'
                        call write_plt_file(plt_file, dble(grid%density_3d_mu))
                        plt_file = trim(root_3d) // '-densityIm_wFMu.plt'
                        call write_plt_file(plt_file, aimag(grid%density_3d_mu))
                     endif
                  else if(trim(grid%format_).eq.'cube') then
                     !write total density cube files
                     cube_file = trim(root_3d) // '-densityRe.cube'
                     call write_cube_file(cube_file, dble(grid%density_3d))
                     cube_file = trim(root_3d) // '-densityIm.cube'
                     call write_cube_file(cube_file, aimag(grid%density_3d))
                     if (control_analysis%separate_q_mu) then
                        !write FQ density cube files
                        cube_file = trim(root_3d) // '-densityRe_wFQ.cube'
                        call write_cube_file(cube_file, dble(grid%density_3d_q))
                        cube_file = trim(root_3d) // '-densityIm_wFQ.cube'
                        call write_cube_file(cube_file, aimag(grid%density_3d_q))
                        !write FMu density cube files
                        cube_file = trim(root_3d) // '-densityRe_wFMu.cube'
                        call write_cube_file(cube_file, dble(grid%density_3d_mu))
                        cube_file = trim(root_3d) // '-densityIm_wFMu.cube'
                        call write_cube_file(cube_file, aimag(grid%density_3d_mu))
                     endif
                  endif
                  deallocate(grid%density_3d)
                  if (control_analysis%separate_q_mu) then
                     deallocate(grid%density_3d_q)
                     deallocate(grid%density_3d_mu)
                  endif
               endif ! field or density for cube or plt
               call print_quantity_analysis_3d(control_analysis%what, &
                                               j,                     &
                                               max_quantity(:,j),     &
                                               i_max,                 &
                                               field%polarization_name(j))
            endif !plane, volume, cube or plt
         enddo ! polarization
      enddo !freq
  
   end subroutine plot_field_or_density


   !> Subroutine for plotting field on planes
   !!    Input  : i_pol        -- polarization index
   !!    Input  : variables_w  -- freq-dep variables
   !!    Input  : string_freq  -- frequency as a string
   subroutine plot_field_on_planes(i_pol,variables_w,string_freq)
       
      implicit none
  
      !input/output variables
      integer, intent(in)           :: i_pol
      character(len=12), intent(in) :: string_freq
      complex(dp), dimension(target_%n_var), intent(in) :: variables_w
       
      !internal variables
      integer :: i 
      real(dp) :: third_coord
      character(len=20)   :: string_coord
      character(len=1000) :: csv_file
      character(len=21) :: format_1 = "(1x,a,i9,2x,a,f20.2/)"
  
      allocate(grid%field_2d(grid%n_2d_points(1), grid%n_2d_points(2)))
  
      do i = 1, control_analysis%n_plane
         !third coordinate is updated
         third_coord = control_analysis%start_plane + &
                       (i-1)*control_analysis%step_plane
         write(string_coord,'(f20.2)') third_coord
         string_coord = adjustl(string_coord)
         write(out_analysis%iunit,format_1) 'Plane: ', i, 'Coord: ', third_coord
         flush(out_analysis%iunit)
  
         call grid%calculate_field_2d(field%index_rhs_polarization(i_pol), &
                                      variables_w,third_coord)
  
         !write csv
         csv_file = trim(out_analysis%plane_folder) //        &
                    trim(out_analysis%root_filename)// "-" // & 
                    field%polarization_name(i_pol)// "-"//    &
                    trim(string_freq)//"-p-"//                &
                    trim(string_coord)//'.csv'
         call write_csv_file(csv_file, grid%field_2d)
      enddo
      deallocate(grid%field_2d)
  
   end subroutine plot_field_on_planes


   !> Subroutine for plotting density on planes
   !!    Input  : i_pol        -- polarization index
   !!    Input  : variables_w  -- freq-dep variables
   !!    Input  : string_freq  -- frequency as a string
   subroutine plot_density_on_planes(i_pol,variables_w,string_freq)
       
      implicit none
  
      !input/output variables
      integer, intent(in)           :: i_pol
      character(len=12), intent(in) :: string_freq
      complex(dp), dimension(target_%n_var), intent(in) :: variables_w
       
      !internal variables
      integer :: i 
      real(dp) :: third_coord
  
      character(len=20)   :: string_coord
      character(len=1000) :: csv_file
  
      character(len=21) :: format_1 = "(1x,a,i9,2x,a,f20.2/)"
  
      allocate(grid%density_2d(grid%n_2d_points(1), grid%n_2d_points(2)))
  
      do i = 1, control_analysis%n_plane
         !third coordinate is updated for each plane
         third_coord = control_analysis%start_plane + &
                       (i-1)*control_analysis%step_plane
         write(string_coord,'(f20.2)') third_coord
         string_coord = adjustl(string_coord)
         write(out_analysis%iunit,format_1) 'Plane: ', i, 'Coord: ', third_coord
         flush(out_analysis%iunit)
          
         call grid%calculate_density_2d(variables_w,third_coord)
  
         !write csv real part
         csv_file = trim(out_analysis%plane_folder) //         &
                    trim(out_analysis%root_filename)// "-" //  & 
                    field%polarization_name(i_pol)// "-"//     &
                    trim(string_freq)//"-p-"//                 &
                    trim(string_coord)//'-densityRe.csv'
         call write_csv_file(csv_file, dble(grid%density_2d))
         !write csv imaginary part
         csv_file = trim(out_analysis%plane_folder) //         &
                    trim(out_analysis%root_filename)// "-" //  & 
                    field%polarization_name(i_pol)// "-"//     &
                    trim(string_freq)//"-p-"//                 &
                    trim(string_coord)//'-densityIm.csv'
         call write_csv_file(csv_file, aimag(grid%density_2d))
      enddo
      deallocate(grid%density_2d)
  
   end subroutine plot_density_on_planes


   !> Subroutine for writing PQR files
   !!    Input  : pqr_file     -- name of the file
   !!    Input  : variables    -- variables to be written (Re o Im)
   !!    Input  : index_rhs    -- which component of the RHS
   subroutine write_pqr_file(pqr_file, variables, index_rhs)
       
      implicit none
  
      !input/output variables
      integer,intent(in) :: index_rhs
      character(len=1000), intent(in) :: pqr_file
      real(dp), dimension(target_%n_var, 3), intent(in) :: variables
  
      !internal variables
      integer :: unit_pqr = 14
      integer :: iost
      integer :: k
      character(len=65) :: format_pqr = &
                   "('ATOM',3X,I4,2x,a2,2x,'  X',5x,'1',4x,3f8.3,E11.3,1x,f5.3)"
      character(len=65) :: format_pqr_2 = &
                   "('ATOM',3X,I4,2x,a2,2x,'DUM',5x,'1',4x,3f8.3,E11.3,1x,f5.3)"
  
      write(out_analysis%iunit,'(1x,a)') "Created PQR file: "// trim(pqr_file)
      write(out_analysis%iunit,out_analysis%sticks)
      
      open(unit=unit_pqr,file=trim(pqr_file),status="unknown",iostat=iost)
         if (iost.ne.0) call out_analysis%error('Error in opening(w) '// &
                                                trim(pqr_file))
         if (target_%name_.ne.'bem') then
            do k = 1, target_%n_var
               write(unit_pqr,format_pqr) k, adjustl(target_%atom_name(k)), &
                                          target_%coord(1,k),               &
                                          target_%coord(2,k),               &
                                          target_%coord(3,k),               &
                                          variables(k,index_rhs),           &
                                          assign_atomic_radius(target_%atom_name(k))
            enddo
         else
            do k = 1, target_%n_var
               write(unit_pqr,format_pqr_2) k, "X", &
                                          target_%coord(1,k),               &
                                          target_%coord(2,k),               &
                                          target_%coord(3,k),               &
                                          variables(k,index_rhs),           &
                                          1.0d0
            enddo
         endif
      close(unit_pqr)
  
   end subroutine write_pqr_file


   !> Subroutine for writing plt files
   !!    Input  : plt_file     -- name of the file
   !!    Input  : quantity     -- quantity to be written in the grid
   subroutine write_plt_file(plt_file, quantity)
       
      implicit none
  
      !input/output variables
      character(len=1000) :: plt_file
      real(dp), dimension(grid%nx*grid%ny*grid%nz) :: quantity
  
      !internal variables
      integer :: unit_plt = 14
      integer :: ntot, nlength
      integer(i6)  :: int1, int2
       
      write(out_analysis%iunit,'(1x,a)') "Created PLT file: "// trim(plt_file)
      write(out_analysis%iunit,out_analysis%sticks)
       
      ntot = grid%nx * grid%ny * grid%nz + 11
      nlength = ntot*4
      
      open(unit=unit_plt, file=trim(plt_file), access='direct', &
           action='write', recl=nlength, status='unknown', &
           form='unformatted')
         int1 = 3 
         int2 = 200 
         write(unit_plt, rec=1) int1, &                    
                                int2, &                    
                                int(grid%nz,i6), &
                                int(grid%ny,i6), &
                                int(grid%nx,i6), &  
                                real(grid%min_coord(3), 4), &
                                real(grid%max_coord(3), 4), &      
                                real(grid%min_coord(2), 4), &
                                real(grid%max_coord(2), 4), &      
                                real(grid%min_coord(1), 4), &
                                real(grid%max_coord(1), 4), &      
                                real(quantity, 4)
      close(unit_plt)
  
   end subroutine write_plt_file


   !> Subroutine for writing cube files
   !!    Input  : cube_file    -- name of the file
   !!    Input  : quantity     -- quantity to be written in the grid
   subroutine write_cube_file(cube_file, quantity)
       
      implicit none
  
      !input/output variables
      character(len=1000) :: cube_file
      real(dp), dimension(grid%nx,grid%ny,grid%nz) :: quantity
  
      !internal variables
      integer :: unit_cube = 14
      integer :: i, j, k
      integer :: counter
      character(len=11) :: format_1     = "(i5,3f12.6)"
      character(len=20) :: format_atoms = "(i5,1x,f10.6,3f12.6)"
  
      write(out_analysis%iunit,'(1x,a)') "Created CUBE file: "// trim(cube_file)
      write(out_analysis%iunit,out_analysis%sticks)
      
      open(unit=unit_cube,file=trim(cube_file),status="unknown") 
  
         ! Write the header information
         write(unit_cube, '(a)') 'Generated by plasmonX analysis'
         write(unit_cube, '(a)') ''
  
         ! Write the number of atoms and the origin of the grid
         if(target_%name_.ne.'bem') then
            write(unit_cube, format_1) target_%n_atoms, grid%min_coord(:)*tobohr
         else
            write(unit_cube, format_1) target_%n_var, grid%min_coord(:)*tobohr
         endif
  
         ! Write the grid information
         write(unit_cube, format_1) grid%nx, grid%step(1)*tobohr, 0.0, 0.0
         write(unit_cube, format_1) grid%ny, 0.0, grid%step(2)*tobohr, 0.0
         write(unit_cube, format_1) grid%nz, 0.0, 0.0, grid%step(3)*tobohr
  
         if (target_%name_.ne.'bem') then
            do i = 1, target_%n_atoms
               write(unit_cube, format_atoms) int(target_%atomic_number(i)), &
                                              0.0, target_%coord(:,i)*tobohr
            enddo
         else
            do i = 1, target_%n_var
               write(unit_cube, format_atoms) 0, &
                                              0.0, target_%coord(:,i)*tobohr
            enddo
         endif
  
         do i = 1, grid%nx
            do j= 1, grid%ny
               counter = 0
               do k = 1, grid%nz
                   if(abs(quantity(i,j,k)).gt.1.0d-99) then
                      write(unit_cube, '(1x, E14.6)', advance='no') &
                         quantity(i,j,k)
                   else
                      write(unit_cube, '(1x, E14.6)', advance='no') zero
                   endif
                   counter = counter + 1
                   if (mod(counter, 6) .eq. 0) write(unit_cube, *)
               end do
               if (mod(counter, 6) .ne. 0) write(unit_cube, *)
            end do
         end do
      close(unit_cube)
  
   end subroutine write_cube_file


   !> Subroutine for writing csv files
   !!    Input  : cube_file    -- name of the file
   !!    Input  : quantity     -- quantity to be written in the grid
   subroutine write_csv_file(csv_file, quantity)
       
      implicit none
  
      !input/output variables
      character(len=1000) :: csv_file
      real(dp), dimension(grid%n_2d_points(1), grid%n_2d_points(2)) :: quantity
  
      !internal variables
      integer :: unit_csv = 14
      integer :: i, j
      integer :: i_1 = 0 
      integer :: i_2 = 0
      real(dp) :: x, y
      character(len=14) :: format_csv = "(3(E25.16,3x))"
  
      write(out_analysis%iunit,'(1x,a)') "Created CSV file: "// trim(csv_file)
      write(out_analysis%iunit,out_analysis%sticks) 
      if(control_analysis%plane.eq.'xy') then 
         i_1 = 1
         i_2 = 2
      else if(control_analysis%plane.eq.'xz') then 
         i_1 = 1
         i_2 = 3
      else if(control_analysis%plane.eq.'yz') then 
         i_1 = 2
         i_2 = 3
      endif
      open(unit=unit_csv, file=trim(csv_file)  ,status="unknown")
         do i = 1, grid%n_2d_points(1)
            do j = 1, grid%n_2d_points(2)
                x = grid%min_coord(i_1) + (i-1)*grid%step(i_1)
                y = grid%min_coord(i_2) + (j-1)*grid%step(i_2)
                if(abs(quantity(i,j)).gt.1.0d-99) then
                   write(unit_csv,format_csv)  x, y, quantity(i,j)
                else
                   write(unit_csv,format_csv)  x, y, zero
                endif
                if(j.eq.grid%n_2d_points(2))  write(unit_csv,'(a)') ""
            enddo
         enddo
      close(unit_csv)
  
   end subroutine write_csv_file


   !> Subroutine for analyzing the input quantity (3D)
   !!    Input  : what         -- field or density
   !!    Input  : i_pol        -- index of polarization -- external cycle
   !!    Input  : maximum      -- maximum of the what quantity
   !!    Input  : i_max        -- indices of the maximum of the what quantity
   !!    Input  : pol_name     -- polarization name (Ex, Ey, Ez)
   subroutine print_quantity_analysis_3d(what, i_pol, maximum, i_max, pol_name)
       
      implicit none
  
      !input/output variables
      character(len=*), intent(in) :: what
      integer, intent(in) :: i_pol
      real(dp), dimension(2, 1), intent(in) :: maximum
      integer, dimension(2, 3, 3), intent(in) :: i_max
      character(len=2), intent(in) :: pol_name
  
      !internal variables
      real(dp), dimension(3) :: coord_max
      integer, dimension(3) :: index_max
      character(len=27) :: format_1 = '(2x,a,11x,a,2(16x,a),13x,a)'
      character(len=30) :: format_2 = '(3x,a,2x,3(f15.5,2x),3x,E13.6)'

      !printing in output
      if(trim(what).eq.'field') then
         index_max = i_max(1,i_pol,:)
         coord_max = grid%get_coordinate_index(index_max)
         write(out_analysis%iunit,format_1) "Field","X","Y","Z","Max."
         write(out_analysis%iunit,format_2) trim(pol_name), &
                                            coord_max(1), &
                                            coord_max(2), &
                                            coord_max(3), &
                                            maximum(1,1)
      else if (trim(what).eq.'density') then
         !real
         index_max = i_max(1,i_pol,:)
         coord_max = grid%get_coordinate_index(index_max)
         write(out_analysis%iunit,format_1) "RhoRe","X","Y","Z","Max."
         write(out_analysis%iunit,format_2) trim(pol_name), &
                                            coord_max(1), &
                                            coord_max(2), &
                                            coord_max(3), &
                                            maximum(1,1)
         !imaginary
         index_max = i_max(2,i_pol,:)
         coord_max = grid%get_coordinate_index(index_max)
         write(out_analysis%iunit,format_1) "RhoIm","X","Y","Z","Max."
         write(out_analysis%iunit,format_2) trim(pol_name), &
                                            coord_max(1), &
                                            coord_max(2), &
                                            coord_max(3), &
                                            maximum(2,1)
      endif
      write(out_analysis%iunit,out_analysis%sticks)
  
   end subroutine print_quantity_analysis_3d


   !> Subroutine for printing the effective volume or area (3D or 2D)
   !! This calculates:
   !!   int_V |E_ind(x,y,z)|^2 / |E^max_ind|^2 dV 
   !!   see Nanolett 2015, 15, 3410
   !!    Input  : what         -- field or density
   !!    Input  : pol_name     -- polarization name (Ex, Ey, Ez)
   !!    Input  : integral     -- maximum of the integral in the grid
   !!    Input  : maximum      -- maximum of the field in the grid
   subroutine print_effective_volume_area(what, pol_name, integral, maximum)
       
      implicit none
  
      !input/output variables
      character(len=*), intent(in) :: what
      character(len=2), intent(in) :: pol_name
      real(dp), intent(in) :: integral
      real(dp), intent(in) :: maximum
  
      !internal variables
      real(dp), dimension(3) :: diff_grid
      real(dp) :: effective_volume
      real(dp) :: effective_area
      real(dp) :: minimum
      character(len=1) :: direction

      !printing in output
      if(trim(what).eq.'volume') then
         diff_grid(1) = grid%max_coord(1) - grid%min_coord(1)
         diff_grid(2) = grid%max_coord(2) - grid%min_coord(2)
         diff_grid(3) = grid%max_coord(3) - grid%min_coord(3)
         if (diff_grid(1).lt.diff_grid(2).and.diff_grid(1).lt.diff_grid(3)) then
            minimum = diff_grid(1)
            direction = "X"
         else if (diff_grid(2).lt.diff_grid(1).and.&
                  diff_grid(2).lt.diff_grid(3)) then
            minimum = diff_grid(2)
            direction = "Y"
         else if (diff_grid(3).lt.diff_grid(1).and.&
                  diff_grid(3).lt.diff_grid(2)) then
            minimum = diff_grid(3)
            direction = "Z"
         else
            minimum = zero
            call out_analysis%error("The provided grid cannot be used for &
                                    &effective volume or area calculations")
         endif
         if (minimum .lt. 1.0d-14) &
            call out_analysis%error("The minimum distance in the grid is zero")
         effective_volume = integral/maximum
         effective_area   = effective_volume/minimum
         write(out_analysis%iunit,'(1x,a)') "Calculation of effective volume &
                                         &and area for polarization: "//&
                                         trim(pol_name)
         write(out_analysis%iunit,'(1x,a,f6.3,a)') "The slice has minimum distance&
                                         & in "//direction//" direction: ", &
                                         minimum, " Angstrom."
         write(out_analysis%iunit,'(/,1x,a,f15.5,a)') "Effective Volume : ", &
                                        (effective_volume)/10**3, " nm^3"
         write(out_analysis%iunit,'(1x,a,f15.5,a)') "Effective Area   : ", &
                                        (effective_area)/10**2, " nm^2"
      endif
      write(out_analysis%iunit,out_analysis%sticks)
  
   end subroutine print_effective_volume_area

end module plot_module
