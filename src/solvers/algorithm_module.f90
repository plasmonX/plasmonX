!> Algorithm module
!!
!! This module contains the subroutines for the algorithm
!!
!! Date         : 2025
!!
module algorithm_module
      
   use output_module
   use parameters_module
   use target_module
   use field_module
   use memory_manager_module
   use gmres_module
  
   implicit none
  
   public algorithm
  
   type :: algorithm_type
      integer  :: n_iter = 1000 ! max number of cycles
      integer  :: n_dim_gmres = 1000 ! number of dimension GMRES
      integer  :: norm_for_convergence = 1 !which norm to use
      real(dp) :: threshold  = 1.0d-6      !convergence threshold
      logical  :: inversion = .false.                !inversion
      logical  :: iterative = .false.                !iterative
      logical  :: on_the_fly = .false.               !on the fly
      logical  :: do_not_change_thresholds = .false. !do not change input thr.
      logical  :: adaptive_tuning = .false.          !adaptive tuning algorithm
      logical  :: matrix_in_parallel = .false.       !matrix construction openmp
      logical  :: freq_in_parallel   = .false.       !frequencies openmp
      logical  :: changed_thresholds = .false.       !changed_input_threshold
      logical  :: left_preconditioner = .false.      !GMRES left preconditioner
      logical  :: right_preconditioner = .false.     !GMRES right preconditioner
  
      contains
  
      procedure :: print_info => print_info_algorithm
      procedure :: solve_ground_state
      procedure :: solve_static_field
      procedure :: solve_dynamic_field
   end type algorithm_type
    
   type (algorithm_type), target, save :: algorithm
  
contains

   !> Subroutine for writing the info of algorithm section
   !!    Input  : algorithm   -- algorithm
   subroutine print_info_algorithm(algorithm)
  
      implicit none
       
      !input/output variables
      class(algorithm_type), intent(in) :: algorithm
       
      if(algorithm%matrix_in_parallel) then
         write(out_%iunit,'(1x,a,l1)') "Parallel execution : matrix"
      else if (algorithm%freq_in_parallel) then
         write(out_%iunit,'(1x,a,l1)') "Parallel execution : frequencies"
      endif
      if (algorithm%inversion) then
         write(out_%iunit,'(1x,a)') "Algorithm method   : Inversion"
      else if (algorithm%iterative) then 
         if (algorithm%on_the_fly) then
            write(out_%iunit,'(1x,a)') "Algorithm method   : Iterative on the fly"
         else
            write(out_%iunit,'(1x,a)') "Algorithm method   : Iterative"
         endif
      endif
      write(out_%iunit,'(1x,a,l1)') "Adaptive Tuning    : ", &
         algorithm%adaptive_tuning
      if(algorithm%iterative.or.algorithm%on_the_fly) then
         write(out_%iunit,'(1x,a,i5)') "Iteration number   :", algorithm%n_iter
         write(out_%iunit,'(1x,a,i5)') "Max GMRES cycles   :", &
                                       algorithm%n_dim_gmres
         write(out_%iunit,'(1x,a,e10.3)') "Tolerance          :", &
                                       algorithm%threshold/ &
                                       sqrt(dble(target_%n_var))
         write(out_%iunit,'(1x,a,e10.3)') "Tolerance (RMSE)   :", &
                                       algorithm%threshold
      endif
      write(out_%iunit,out_%sticks) 
      flush(out_%iunit)
  
   end subroutine print_info_algorithm


   !> Subroutine for solving ground state
   !!    In/Out : algorithm   -- algorithm
   subroutine solve_ground_state(algorithm)
  
      implicit none
       
      !input/output variables
      class(algorithm_type) :: algorithm
       
      if(algorithm%inversion) then 
        !allocation of needed variables
        call target_%allocate_gs_matrices()
        !construct matrix
        call target_%construct_static_matrix()
        !construct right hand side
        call target_%construct_ground_state_rhs()
        !solve by inversion
        call real_inversion()
      else
        call out_%error('Iterative not implemented for ground state')
      endif
  
      if(out_%ivrb.ge.2) call target_%print_gs_variables()
      call target_%calculate_energy()
      call target_%deallocate_gs_matrices()
       
   end subroutine solve_ground_state


   !> Subroutine for solving ground state
   !!    In/Out : algorithm   -- algorithm
   subroutine solve_static_field(algorithm)
  
      implicit none
       
      !input/output variables
      class(algorithm_type) :: algorithm
       
      !allocation of needed variables
      call target_%allocate_static_field_matrices()
      !construct matrix
      call target_%construct_static_matrix()
      !construct right hand side
      call target_%construct_static_field_rhs()
      !solve system
      if(algorithm%inversion) then 
         call real_inversion()
      else
         call out_%error('Iterative not implemented for static field')
      endif
  
      if(out_%ivrb.ge.2) call target_%print_static_field_variables()
      !calculate static polar
      call target_%calculate_static_polar()
      call target_%deallocate_static_field_matrices()
       
   end subroutine solve_static_field


   !> Subroutine for solving ground state
   !!    In/Out : algorithm   -- algorithm
   subroutine solve_dynamic_field(algorithm)
      
      implicit none
       
      !input/output variables
      class(algorithm_type) :: algorithm
  
      !internal variables
      complex(dp), dimension(:,:), allocatable :: variables_w
      complex(dp), dimension(:,:), allocatable :: polar_w
      complex(dp), dimension(:,:), allocatable :: rhs_w
       
      call target_%allocate_dynamic_field_general(variables_w,rhs_w,polar_w)
      !solve depending on algorithm
      if(algorithm%inversion) then
         call solve_w_inversion(variables_w,rhs_w,polar_w)
      else if(algorithm%iterative) then 
         if(algorithm%on_the_fly) then
            call solve_w_iter_on_the_fly(variables_w,rhs_w,polar_w)
         else
            call solve_w_iter_mem(variables_w,rhs_w,polar_w)
         endif
      endif
      call out_%clean_up_scratch()
      call target_%print_dynamic_results()
      if(field%n_freq.gt.1) call print_maxima()
      if(field%n_freq.gt.1) call target_%save_csv_file()
      call target_%deallocate_dynamic_field_general(variables_w,rhs_w,polar_w)
       
   end subroutine solve_dynamic_field

!-------------------------------------------------------------------------------
!  INTERNAL SUBROUTINES 
!-------------------------------------------------------------------------------

   !> Subroutine for solving w inversion
   !!    In/Out : variables_w -- w-variables
   !!    In/Out : rhs_w       -- w-RHS
   !!    In/Out : polar_w     -- dynamic polar
   subroutine solve_w_inversion(variables_w,rhs_w,polar_w)
  
      use control_module
      use bem_module
       
      implicit none
  
      !input/output variables
      complex(dp), dimension(target_%n_var,3), intent(inout) :: variables_w
      complex(dp), dimension(target_%n_var,3), intent(inout) :: rhs_w
      complex(dp), dimension(3,3), intent(inout) :: polar_w
  
      !internal variables
      integer :: i
      integer :: i_cycle
      real(dp), dimension(:,:), allocatable    :: matrix_constant
      complex(dp), dimension(:,:), allocatable :: matrix_w
      logical :: check_exist
      !for bem charge constraint
      complex(dp), dimension(:,:), allocatable :: rhs_1
      complex(dp), dimension(:,:), allocatable :: variables_1
  
      !allocation of dynamic field matrices
      call target_%allocate_dynamic_field_memory(matrix_constant, matrix_w)
      if(bem%charge_constraint) &
         call target_%allocate_constant_potential_variables(variables_1,rhs_1)
  
      !construct matrix
      call target_%construct_constant_matrix(matrix_constant)
      !construct rhs
      call target_%construct_dynamic_field_rhs(rhs_w)
  
      i_cycle = 0
      if(.not.bem%charge_constraint) then
         !cycles over frequencies
         !$omp parallel do firstprivate(matrix_constant,rhs_w)     &
         !$omp private(matrix_w,variables_w,polar_w) private(i) &
         !$omp schedule(static) if(algorithm%freq_in_parallel)
         do i = 1, field%n_freq
            !$omp critical
            if(control%restart) then 
               call check_if_restart(i,check_exist)
            else
               check_exist = .false.
            endif
            !$omp end critical
            call check_i_cycle(i_cycle)
            if(.not.check_exist) then
               !construct dynamic matrix
               call target_%construct_dynamic_matrix(i,matrix_constant,matrix_w)
               !assign the w-variables
               call target_%assign_variables_w(rhs_w,variables_w)
               !solve complex system by inversion
               call complex_inversion(matrix_w,variables_w)
               !$omp critical
               if(out_%ivrb.ge.2) & 
                  call target_%print_dynamic_field_variables(variables_w)
               !$omp end critical
               !calculate dynamic polar
               call target_%calculate_dynamic_polar(i,variables_w,polar_w)
               !print dynamic polar
               call target_%print_dynamic_polar(i,polar_w)
               !calculate cross sections
               call target_%calculate_cross_section(i)
               !$omp critical   
               if(control%save_info) call save_dynamic_variables(i,variables_w)
               !$omp end critical
               !$omp critical   
               call target_%save_intermediate_results(i)
               !$omp end critical
            else !check_exist: there is a backup
               !$omp critical   
               call recover_from_backup_file(i)
               !$omp end critical
            endif !check_exist
         enddo
         !$omp end parallel do 
      else if (bem%charge_constraint) then
         !cycles over frequencies
         !$omp parallel do firstprivate(matrix_constant,rhs_w,rhs_1)     &
         !$omp private(matrix_w,variables_w,polar_w,variables_1) private(i) &
         !$omp schedule(static) if(algorithm%freq_in_parallel)
         do i = 1, field%n_freq
            !$omp critical
            if(control%restart) then 
               call check_if_restart(i,check_exist)
            else
               check_exist = .false.
            endif
            !$omp end critical
            call check_i_cycle(i_cycle)
            if(.not.check_exist) then
               !construct dynamic matrix
               call target_%construct_dynamic_matrix(i,matrix_constant,matrix_w)
               !assign variables w
               call target_%assign_variables_w(rhs_w,variables_w)
               !solve complex system by inversion
               call complex_inversion(matrix_w,variables_w)
               !reassign variables_1 to rhs_1 (see above)
               variables_1 = rhs_1
               !solve complex system by inversion
               call complex_inversion(matrix_w,variables_1)
               variables_w(:,1) = variables_w(:,1) - variables_1(:,1)* &
                                     sum(variables_w(:,1))/sum(variables_1(:,1))
               variables_w(:,2) = variables_w(:,2) - variables_1(:,2)* &
                                     sum(variables_w(:,2))/sum(variables_1(:,2))
               variables_w(:,3) = variables_w(:,3) - variables_1(:,3)* &
                                     sum(variables_w(:,3))/sum(variables_1(:,3))
               !$omp critical
               if(out_%ivrb.ge.2) &
                  call target_%print_dynamic_field_variables(variables_w)
               !$omp end critical
               !calculate dynamic polar
               call target_%calculate_dynamic_polar(i,variables_w,polar_w)
               !print dynamic polar
               call target_%print_dynamic_polar(i,polar_w)
               !calculate cross section
               call target_%calculate_cross_section(i)
               !$omp critical   
               if(control%save_info) call save_dynamic_variables(i,variables_w)
               !$omp end critical
               !$omp critical   
               call target_%save_intermediate_results(i)
               !$omp end critical
            else !check_exist
               !$omp critical   
               call recover_from_backup_file(i)
               !$omp end critical
            endif !check_exist
         enddo
         !$omp end parallel do 
      endif

      !deallocation of dynamic field matrices
      call target_%deallocate_dynamic_field_memory(matrix_w, matrix_constant)
      if(bem%charge_constraint) &
         call target_%deallocate_constant_potential_variables(variables_1,rhs_1)
  
   end subroutine solve_w_inversion


   !> Subroutine for solving w with iterative algorithm in memory
   !!    In/Out : variables_w -- w-variables
   !!    In/Out : rhs_w       -- w-RHS
   !!    In/Out : polar_w     -- dynamic polar
   subroutine solve_w_iter_mem(variables_w,rhs_w,polar_w)
  
      use control_module
       
      implicit none
  
      !input/output variables
      complex(dp), dimension(target_%n_var,3), intent(inout) :: variables_w
      complex(dp), dimension(target_%n_var,3), intent(inout) :: rhs_w
      complex(dp), dimension(3,3), intent(inout)             :: polar_w
  
      !internal variables
      integer :: i
      integer :: i_cycle
      real(dp), dimension(:,:), allocatable    :: matrix_constant
      complex(dp), dimension(:,:), allocatable :: matrix_iterative
      logical :: check_exist
  
      !allocation of dynamic field matrices
      call target_%allocate_dynamic_field_memory(matrix_constant, &
                                                 matrix_iterative)
      !construct matrix
      call target_%construct_constant_matrix(matrix_constant)
      !construct matrix for gmres
      matrix_iterative = dcmplx(matrix_constant,zero)
      call mem_man%dealloc(matrix_constant, "matrix_constant")
      !construct rhs
      call target_%construct_dynamic_field_rhs(rhs_w)
  
      i_cycle = 0
      !cycles over frequencies
      !$omp parallel do firstprivate(matrix_iterative,rhs_w)  &
      !$omp private(variables_w,polar_w) private(i)           &
      !$omp schedule(static) if(algorithm%freq_in_parallel)
      do i = 1, field%n_freq
         !$omp critical
         if(control%restart) then 
            call check_if_restart(i,check_exist)
         else
          check_exist = .false.
         endif
         !$omp end critical
         call check_i_cycle(i_cycle)
         if(.not.check_exist) then
            if(target_%heterogeneous) & !here it is changed if it is hetero
               call target_%construct_dynamic_matrix_gmres(i,matrix_iterative) 
            !assign variables w
            call target_%assign_variables_w(rhs_w,variables_w) 
            !$omp critical
            !call GMRES algorithm iterative on memory
            call complex_gmres_iterative(i,matrix_iterative,variables_w)
            !$omp end critical
            !$omp critical
            if(out_%ivrb.ge.2) &
               call target_%print_dynamic_field_variables(variables_w)
            !$omp end critical
            call target_%calculate_dynamic_polar(i,variables_w,polar_w)
            call target_%print_dynamic_polar(i,polar_w)
            call target_%calculate_cross_section(i)
            !$omp critical   
            if(control%save_info) call save_dynamic_variables(i,variables_w)
            !$omp end critical
            !$omp critical   
            call target_%save_intermediate_results(i)
            !$omp end critical
         else !check_exist
            !$omp critical   
            call recover_from_backup_file(i)
            !$omp end critical
         endif !check_exist
      enddo
      !$omp end parallel do 

      !deallocation of dynamic field matrices
      call target_%deallocate_dynamic_field_memory(matrix_iterative)
  
   end subroutine solve_w_iter_mem


   !> Subroutine for solving w with iterative on the fly
   !!    In/Out : variables_w -- w-variables
   !!    In/Out : rhs_w       -- w-RHS
   !!    In/Out : polar_w     -- dynamic polar
   subroutine solve_w_iter_on_the_fly(variables_w,rhs_w,polar_w)
  
      use control_module
       
      implicit none
  
      !input/output variables
      complex(dp), dimension(target_%n_var,3), intent(inout) :: variables_w
      complex(dp), dimension(target_%n_var,3), intent(inout) :: rhs_w
      complex(dp), dimension(3,3), intent(inout)             :: polar_w
  
      !internal variables
      integer :: i
      integer :: i_cycle
      logical :: check_exist
  
      !construct rhs on the fly
      call target_%construct_dynamic_field_rhs_on_the_fly(rhs_w)
  
      i_cycle = 0
      !cycles over frequencies
      do i = 1, field%n_freq
         if(control%restart) then 
            call check_if_restart(i,check_exist)
         else
            check_exist = .false.
         endif
         call check_i_cycle(i_cycle)
         if(.not.check_exist) then
            !assign w variables 
            call target_%assign_variables_w_on_the_fly(i, rhs_w, variables_w)
            !solve complex system with GMRES on the fly
            call complex_gmres_iterative_on_the_fly(i,variables_w)
            if(out_%ivrb.ge.2) &
               call target_%print_dynamic_field_variables(variables_w)
            !calculate dynamic polar
            call target_%calculate_dynamic_polar(i,variables_w,polar_w)
            !print dynamic polar
            call target_%print_dynamic_polar(i,polar_w)
            !calculate cross section
            call target_%calculate_cross_section(i)
            if(control%save_info) call save_dynamic_variables(i,variables_w)
            call target_%save_intermediate_results(i)
         else !check_exist
            call recover_from_backup_file(i)
         endif !check_exist
      enddo

   end subroutine solve_w_iter_on_the_fly


   !> Subroutine for solving real system with inversion
   subroutine real_inversion()
  
      implicit none
  
      !internal variables
      integer :: info
      integer :: i, j
      integer, dimension(target_%n_var)     :: IPV
      real(dp), dimension(target_%n_var)    :: Work
      real(dp), dimension(:,:), allocatable :: matrix_cp

      call mem_man%alloc(matrix_cp, target_%n_var,target_%n_var, "matrix_cp")

      !$omp parallel do collapse(2)
      do i = 1, target_%n_var
         do j = 1, target_%n_var
            matrix_cp(i,j) = target_%matrix(i,j)
         enddo
      enddo
      !$omp end parallel do
      target_%variables = target_%rhs
      call dsysv('U',               &
                 target_%n_var,     &
                 target_%n_rhsre,   &
                 matrix_cp,         &
                 target_%n_var,     &
                 IPV,               &
                 target_%variables, &
                 target_%n_var,     &
                 work,              &
                 target_%n_var,     &
                 info)
      if(info.ne.0) call out_%error('INFO.ne.0 in DSYSV in real inversion')
      call mem_man%dealloc(matrix_cp, "matrix_cp")
       
   end subroutine real_inversion


   !> Subroutine for solving complex system with inversion 
   !!    In/Out : matrix_w        -- dynamic matrix
   !!    In/Out : variables_w     -- w-variables
   subroutine complex_inversion(matrix_w,variables_w)
  
      implicit none
  
      !input/output variables
      complex(dp), dimension(target_%n_var,3), intent(inout) :: variables_w
      complex(dp), dimension(target_%n_var,target_%n_var), intent(inout) :: &
                                                                        matrix_w
  
      !internal variables
      integer :: info
      integer, dimension(target_%n_var) :: ipiv
  
      call zgetrf(target_%n_var, &
                  target_%n_var, &
                  matrix_w,      &
                  target_%n_var, &
                  ipiv,          &
                  info)
      if(info.ne.0) call out_%error('INFO.ne.0 in ZGETRF in complex inversion')
      call zgetrs('n',                 &
                  target_%n_var,       &
                  3,                   &
                  matrix_w,            &
                  target_%n_var,       &
                  ipiv,                &
                  variables_w, &
                  target_%n_var,       &
                  info)
      if(info.ne.0) call out_%error('INFO.ne.0 in ZGETRS in complex inversion')
  
   end subroutine complex_inversion


   !> Subroutine for solving complex system with GMRES iterative (memory)
   !!    Input  : i_freq            -- index of frequency
   !!    In/Out : matrix_iterative  -- dynamic matrix
   !!    In/Out : variables_w       -- w-variables
   subroutine complex_gmres_iterative(i_freq,matrix_iterative,variables_w)
  
      implicit none
  
      !input/output variables
      integer, intent(in) :: i_freq
      complex(dp), dimension(target_%n_var,target_%n_var), intent(in) :: &
                                                               matrix_iterative
      complex(dp), dimension(target_%n_var,3), intent(inout) :: variables_w
  
      !internal variables
      real(dp) :: limit_ = 1.0d-15
      real(dp) :: freq_eV

      call freqautoev(field%freq(i_freq), freq_eV)
      write(out_%iunit,'(/1x,a,f5.3,a,e13.5)') "GMRES Iterative Algorithm. &
                                               &Freq: ", freq_eV, " eV. &
                                               &Threshold =", &
                                               algorithm%threshold
  
      if((field%polarization.eq.'x'   .or.  &
         field%polarization.eq.'all').and. &
         any(dble(variables_w(:,1)).gt.limit_)) then
         write(out_%iunit,"(/1x,a)") "- Polarization: X"
         call complex_gmres_worker(i_freq,matrix_iterative,variables_w(:,1))
      endif
      if((field%polarization.eq.'y'   .or.  &
         field%polarization.eq.'all').and. &
         any(dble(variables_w(:,2)).gt.limit_)) then
         write(out_%iunit,"(/1x,a)") "- Polarization: Y"
         call complex_gmres_worker(i_freq,matrix_iterative,variables_w(:,2))
      endif
      if((field%polarization.eq.'z'   .or.  &
         field%polarization.eq.'all').and. &
         any(dble(variables_w(:,3)).gt.limit_)) then
         write(out_%iunit,"(/1x,a)") "- Polarization: Z"
         call complex_gmres_worker(i_freq,matrix_iterative,variables_w(:,3))
      endif
       
  end subroutine complex_gmres_iterative


   !> Subroutine for solving complex system with GMRES iterative (memory)
   !!    Input  : i_freq            -- index of frequency
   !!    In/Out : variables_w       -- w-variables
   subroutine complex_gmres_iterative_on_the_fly(i_freq,variables_w)
  
      implicit none
  
      !input/output variables
      integer :: i_freq
      complex(dp), dimension(target_%n_var,3), intent(inout) :: variables_w
  
      !internal variables
      real(dp) :: limit_ = 1.0d-15
      real(dp) :: freq_eV
       
      call freqautoev(field%freq(i_freq), freq_eV)
      write(out_%iunit,'(/1x,a,f5.3,a,e13.5)') "GMRES Iterative Algorithm. &
                                               &Freq: ", freq_eV, " eV. &
                                               &Threshold =", &
                                               algorithm%threshold
  
      if((field%polarization.eq.'x'   .or.  &
         field%polarization.eq.'all').and. &
         any(dble(variables_w(:,1)).gt.limit_)) then
         write(out_%iunit,"(/1x,a)") "- Polarization: X"
         call complex_gmres_worker_on_the_fly(i_freq,variables_w(:,1))
      endif
      if((field%polarization.eq.'y'   .or.  &
         field%polarization.eq.'all').and. &
         any(dble(variables_w(:,2)).gt.limit_)) then
         write(out_%iunit,"(/1x,a)") "- Polarization: Y"
         call complex_gmres_worker_on_the_fly(i_freq,variables_w(:,2))
      endif
      if((field%polarization.eq.'z'   .or.  &
         field%polarization.eq.'all').and. &
         any(dble(variables_w(:,3)).gt.limit_)) then
         write(out_%iunit,"(/1x,a)") "- Polarization: Z"
         call complex_gmres_worker_on_the_fly(i_freq,variables_w(:,3))
      endif
       
   end subroutine complex_gmres_iterative_on_the_fly


    !> Subroutine for complex GMRES (worker)
    !!    Input  : i_freq            -- index of frequency
    !!    Input  : matrix_iterative  -- dynamic matrix
    !!    In/Out : variables_w       -- w-variables
    subroutine complex_gmres_worker(i_freq,matrix_iterative,variables_w)
   
       implicit none
   
       !input/output variables
       integer :: i_freq
       complex(dp), dimension(target_%n_var,target_%n_var), intent(in) :: &
                                                                 matrix_iterative
       complex(dp), dimension(target_%n_var), intent(inout) :: variables_w
 
       call gmres_solver(i_freq           = i_freq,              &
                         n                = target_%n_var,       &
                         x                = variables_w,         &
                         matrix_iterative = matrix_iterative,    &
                         onthefly         = .false.,             &
                         tol              = algorithm%threshold, &
                         max_iter         = algorithm%n_iter,    &
                         restart          = algorithm%n_dim_gmres)
 
    end subroutine complex_gmres_worker
 
 
    !> Subroutine for complex GMRES (worker) 
    !!    Input  : i_freq            -- index of frequency
    !!    In/Out : variables_w       -- w-variables
    subroutine complex_gmres_worker_on_the_fly(i_freq,variables_w)
   
       implicit none
   
       !input/output variables
       integer :: i_freq
       complex(dp), dimension(target_%n_var), intent(inout) :: variables_w
 
       call gmres_solver(i_freq           = i_freq,              &
                         n                = target_%n_var,       &
                         x                = variables_w,         &
                         onthefly         = .true.,              &
                         tol              = algorithm%threshold, &
                         max_iter         = algorithm%n_iter,    &
                         restart          = algorithm%n_dim_gmres)
 
    end subroutine complex_gmres_worker_on_the_fly


   !> Subroutine for checking at what cycle we are 
   !!    Input  : i_cycle   -- index of the cycle
   subroutine check_i_cycle(i_cycle)
     
      use field_module
      use output_module
      !$ use omp_lib
       
      implicit none
  
      !input/output variables
      integer :: i_cycle
  
      !internal variables
      integer :: i_thread
      integer :: n_threads
      integer :: n_cycles
  
      n_threads = 1
      !$ i_thread = OMP_GET_THREAD_NUM()
      !$ if(i_thread .eq. 0) then 
      !$ n_threads   = OMP_GET_NUM_THREADS() 
         n_cycles   = ceiling ((1.0d0*field%n_freq) / n_threads)
         i_cycle    = i_cycle + 1
         if(i_cycle.ne.1) write(out_%iunit,out_%sticks)
         write(out_%iunit,'(a6,i6,a7,i6)') ' Cycle ',i_cycle,' out of ',n_cycles
         flush(out_%iunit)
      !$ endif
  
   end subroutine check_i_cycle


   !> Subroutine for checking if some bk file is present
   !!    Input  : i_freq      -- index of the frequency
   !!    In/Out : check_exist -- logical for exist
   subroutine check_if_restart(i_freq,check_exist)
  
      use field_module
     
      implicit none
       
      !input/output variables
      integer           :: i_freq
      logical           :: check_exist
  
      !internal variables
      character(len=99) :: output_bk
      character(len=7)  :: freq_char
      real(dp)          :: freq_eV
  
      call freqautoev(field%freq(i_freq),freq_ev)
      write(freq_char,'(f7.5)') freq_ev
      write(output_bk,'(a)') out_%filename(1:(len_trim(out_%filename) - 4)) // &
                                                '-' // freq_char // '.plasmonX.bk'
      inquire (file=output_bk,exist=check_exist)
  
   end subroutine check_if_restart


   !> Subroutine for recovering info from the backup file
   !!    Input  : i_freq      -- index of the frequency
   subroutine recover_from_backup_file(i_freq)
     
      use field_module
      use control_module
       
      implicit none
       
      !input/output variables
      integer, intent(in)   :: i_freq
  
      !internal variables
      character(len=99)     :: output_bk
      character(len=7)      :: freq_char
      character(len=81)     :: line
      character(len=26)     :: format_1 = "(17x,e9.3,8x,e9.3,6x,e9.3)"
      character(len=11)     :: format_2 = "(34x,e14.6)"
      integer               :: iunit_bk
      integer               :: iost
      real(dp)              :: freq_ev
      real(dp)              :: freq_au
      real(dp)              :: freq_nm
       
      iunit_bk = 15
      freq_au = field%freq(i_freq)
      call freqautoev(field%freq(i_freq),freq_ev)
      call freqautonm(field%freq(i_freq),freq_nm)
      if (out_%ivrb.ge.1) then
         !$omp critical(write_backup)
         write(out_%iunit, '(a,f7.5,a)') " Recovering from backup frequency: ",&
                                         freq_ev, " eV."
         flush(out_%iunit)
         !$omp end critical(write_backup)
      endif
      write(freq_char,'(f7.5)') freq_ev
      write(output_bk,'(a)') out_%filename(1:(len_trim(out_%filename)-4)) // &
                                              '-' // freq_char // '.plasmonX.bk'
      open(unit=iunit_bk,file=output_bk,status="old",iostat=iost,err=05)
         read(iunit_bk,*) line
         read(iunit_bk,format_1) freq_au,freq_nm,freq_ev
         read(iunit_bk,*)
         read(iunit_bk,format_2) target_%results(i_freq,1)
         read(iunit_bk,format_2) target_%results(i_freq,2)
         read(iunit_bk,format_2) target_%results(i_freq,3)
         read(iunit_bk,format_2) target_%results(i_freq,4)
         read(iunit_bk,format_2) target_%results(i_freq,5)
         read(iunit_bk,format_2) target_%results(i_freq,6)
         read(iunit_bk,format_2) target_%results(i_freq,7)
         read(iunit_bk,format_2) target_%results(i_freq,8)
         read(iunit_bk,format_2) target_%results(i_freq,9)
         read(iunit_bk,format_2) target_%results(i_freq,10)
         read(iunit_bk,format_2) target_%results(i_freq,11)
         read(iunit_bk,format_2) target_%results(i_freq,12)
         read(iunit_bk,format_2) target_%results(i_freq,13)
         read(iunit_bk,format_2) target_%results(i_freq,14)
         read(iunit_bk,format_2) target_%results(i_freq,15)
         read(iunit_bk,format_2) target_%results(i_freq,16)
         read(iunit_bk,format_2) target_%results(i_freq,17)
         read(iunit_bk,format_2) target_%results(i_freq,18)
         read(iunit_bk,format_2) target_%results(i_freq,19)
         read(iunit_bk,format_2) target_%results(i_freq,20)
         read(iunit_bk,*) line
      05 continue      
      close(iunit_bk)
      
   end subroutine recover_from_backup_file


   !> Subroutine for saving the dynamic variables
   !!    Input  : i_freq      -- index of the frequency
   !!    In/Out : variables_w -- w-variables           
   subroutine save_dynamic_variables(i_freq,variables_w)
  
      use target_module
      use bem_module
  
      implicit none
       
      !input/output variables
      integer, intent(in):: i_freq
      complex(dp), dimension(target_%n_var,3), intent(inout) :: variables_w
  
      !internal variables
      integer            :: i
      integer            :: IOST
      integer            :: unit_freq
      real(dp)           :: freq_eV
      character(len=12)  :: freq_char
      character(len=100) :: variable
      character(len=100) :: save_file
      character(len=20) :: format_1 = "(i8,5x,2(f25.16,3x))"
     
      variable = "Charges"
      if(target_%name_.eq.'wfqfmu') variable = "Charges/Dipoles"
  
      call freqautoev(field%freq(i_freq),freq_ev)
  
      write(freq_char,'(f12.10)') freq_ev
      write(save_file,'(a)') out_%filename(1:len_trim(out_%filename)-4) // &
                                                  '-'//trim(freq_char)//'.freq'
  
      if(out_%ivrb.ge.3) then 
         write(out_%iunit,'(1x,a,e9.3,a)') "complex "//trim(variable)//&
                             " for freq = ",freq_ev, " saved in : "//save_file
         write(out_%iunit,out_%sticks)
      endif
      unit_freq = 13
      open(unit=unit_freq,file=save_file,status="unknown",iostat=iost)
         write(unit_freq,'(a,f25.16)') 'complex '//trim(variable)// &
                                                     ': x. freq = ', freq_ev
         ! X 
         do i = 1, target_%n_var
            if (trim(target_%name_).eq. 'bem') then
               if (trim(bem%variant).eq.'dpcm') then
                  write(unit_freq,format_1) i, &
                                       dble(variables_w(i,1))*bem%area(i), &
                                      dimag(variables_w(i,1))*bem%area(i)
               else
                  write(unit_freq,format_1) i, dble(variables_w(i,1)), &
                                               dimag(variables_w(i,1))
               endif
            else
               write(unit_freq,format_1) i, dble(variables_w(i,1)), &
                                            dimag(variables_w(i,1))
            endif
         enddo 
         ! Y 
         write(unit_freq,'(a,f25.16)') 'complex '//trim(variable)// &
                                                     ': y. freq = ', freq_ev
         do i = 1, target_%n_var
            if (trim(target_%name_).eq. 'bem') then
               if (trim(bem%variant).eq.'dpcm') then
                  write(unit_freq,format_1) i, &
                                       dble(variables_w(i,2))*bem%area(i), &
                                      dimag(variables_w(i,2))*bem%area(i)
               else
                  write(unit_freq,format_1) i, dble(variables_w(i,2)), &
                                               dimag(variables_w(i,2))
               endif
            else
               write(unit_freq,format_1) i, dble(variables_w(i,2)), &
                                            dimag(variables_w(i,2))
            endif
         enddo 
         write(unit_freq,'(a,f25.16)') 'complex '//trim(variable)// &
                                                     ': z. freq = ', freq_ev
         ! Z
         do i = 1, target_%n_var
            if (trim(target_%name_).eq. 'bem') then
               if (trim(bem%variant).eq.'dpcm') then
                  write(unit_freq,format_1) i, &
                                       dble(variables_w(i,3))*bem%area(i), &
                                      dimag(variables_w(i,3))*bem%area(i)
               else
                  write(unit_freq,format_1) i, dble(variables_w(i,3)), &
                                               dimag(variables_w(i,3))
               endif
            else
               write(unit_freq,format_1) i, dble(variables_w(i,3)), &
                                            dimag(variables_w(i,3))
            endif
         enddo 
      close(unit_freq)
  
   end subroutine save_dynamic_variables


   !> Subroutine for printing the maxima depending on the required solution
   subroutine print_maxima()
  
      use field_module
      use control_module
  
      implicit none
       
      !internal variables
      integer                            :: i, j 
      integer                            :: n_max
      real(dp)                           :: freq_eV
      real(dp), dimension(field%n_freq) :: freq_max
      real(dp), dimension(field%n_freq) :: results_max
      character(len=40) :: format_1 = "(1x,'Maxima Analysis: NumExFreq = ',i6/)"
      character(len=47) :: format_2 = "(10x,'NState      Freq(eV)     ',a24,&
                                      &' (a.u.)')"
      character(len=27) :: format_3 = "(10x,I3,5X,F10.3,13x,E12.5)"
  
      j    = 0
      n_max = 0
      do i=2,field%n_freq-1
         call freqautoev(field%freq(i), freq_eV)
         !absorption 
         if(trim(control%maxima_analysis).eq."absorption") then
            if(abs(target_%results(i-1,9)).lt.abs(target_%results(i,9)).and.&
               abs(target_%results(i+1,9)).lt.abs(target_%results(i,9))) then
               j = j + 1
               freq_max(j) = freq_eV
               results_max(j) = target_%results(i,9)
            endif
         !scattering
         else if(trim(control%maxima_analysis).eq.'scattering') then
            if(abs(target_%results(i-1,13)).lt.abs(target_%results(i,13)).and. &
               abs(target_%results(i+1,13)).lt.abs(target_%results(i,13))) then
               j = j + 1
               freq_max(j) = freq_eV
               results_max(j) = target_%results(i,13)
            endif
         !excition
         else if(trim(control%maxima_analysis).eq.'exctinction') then
            if(abs(target_%results(i-1,17)).lt.abs(target_%results(i,17)).and. &
               abs(target_%results(i+1,17)).lt.abs(target_%results(i,17))) then
               j = j + 1
               freq_max(j) = freq_eV
               results_max(j) = target_%results(i,17)
            endif
         endif
         n_max = j
      enddo
  
      if(n_max.gt.0) then
         write(out_%iunit,format_1) field%n_freq
         if(trim(control%maxima_analysis).eq."absorption") then 
            write(out_%iunit,format_2) 'Isotr. Abs. Cross. Sec. '
         else if(trim(control%maxima_analysis).eq.'scattering') then 
            write(out_%iunit,format_2) 'Isotr. Sca. Cross. Sec. '
         else if(trim(control%maxima_analysis).eq.'exctinction') then 
            write(out_%iunit,format_2) 'Isotr. Ext. Cross. Sec. '
         endif
         do i = 1,n_max 
            write(out_%iunit,format_3) i, freq_max(i),results_max(i)
         enddo
         write(out_%iunit,out_%sticks) 
      endif
  
   end subroutine print_maxima

end module algorithm_module
