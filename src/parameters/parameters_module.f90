!> Parameters module
!!
!! This module contains the parameters used in the code, the parameters of 
!! the models (chi, eta, etc...), and some subroutines
!!
!! Date         : 2024 
!!
module parameters_module
      
   implicit none

   public ! the public variables
   integer, parameter                       :: dp=kind(1.0d0)
   integer, parameter                       :: i6  = selected_int_kind(6)

   !numbers
   real(dp), parameter                      :: zero   = 0.0d0
   real(dp), parameter                      :: one    = 1.0d0
   real(dp), parameter                      :: two    = 2.0d0
   real(dp), parameter                      :: three  = 3.0d0
   real(dp), parameter                      :: four   = 4.0d0
   real(dp), parameter                      :: five   = 5.0d0
   real(dp), parameter                      :: six    = 6.0d0
   real(dp), parameter                      :: seven  = 7.0d0
   real(dp), parameter                      :: eight  = 8.0d0
   real(dp), parameter                      :: Half   = 0.5d0
   real(dp), parameter                      :: pi     = four*atan(one)

   !tolerance floats
   real(dp), parameter                      :: tol_float = 1.0d-14

   !OMP threads
   integer :: n_threads_omp      ! number of OMP threads
   integer :: n_threads_omp_used ! number of OMP threads used

   !conversion memory
   real(dp), parameter :: bytes_per_gb = 1024.0d0 * 1024.0d0 * 1024.0d0

   !conversion factors
   real(dp), parameter                      :: ToBohr = 1.8897261254578281d0
   real(dp), parameter                      :: ToAng  = 1.0d0/ToBohr
   !S/m to a.u.
   real(dp), parameter                      :: FSmAu  = 0.00000021739d0
   !codata    au to eV https://physics.nist.gov/cgi-bin/cuu/Value?hrev
   real(dp), parameter                      :: au_to_ev =  27.211386245981d0 
   !wikipedia au to cm-1 
   real(dp), parameter                      :: au_to_cm =  219474.63136320d0 

   !nearest neighbour distances
   real(dp), parameter                      :: RAg0 = 2.885d0 * ToBohr
   real(dp), parameter                      :: RAu0 = 2.885d0 * ToBohr
                                               !3.66328d0 * ToBohr
   real(dp), parameter                      :: RNa0 = 6.9226100144600462d0 
   real(dp), parameter                      :: RC0  = 1.418d0 * ToBohr
   real(dp), parameter                      :: ROH0 = 1.65d0
   real(dp), parameter                      :: RAl0 = 2.861d0 * ToBohr
   real(dp), parameter                      :: RCu0 = 2.543463d0 * ToBohr

   !Physical quantities in a.u.
   real(dp), parameter                      :: Light  = 137.035227d0
                                               !needed for Graphene 1.0d6 m/s
   real(dp), parameter                      :: fermi_velocity = 0.4573138778d0

   type :: parameters_type
      ! Graphene geometry
      character(len=200), dimension(:), allocatable :: graphene_geometry
      ! File for wfqfmu polarizability
      character(len=200), dimension(:), allocatable :: wfqfmu_file
      !dynamic epsilon of the material [wfqfmu]
      character(len=200), dimension(:), allocatable :: permittivity_type
      ! number of interacting structures (for bowtie)
      integer, dimension(:), allocatable    :: n_structures  
      ! tau
      real(dp), dimension(:), allocatable   :: tau          
      ! sigma_0 in a.u.
      real(dp), dimension(:), allocatable   :: sigma_0      
      ! scaling factor for tau and sigma_0
      real(dp), dimension(:), allocatable   :: scaling      
      ! Effective Area
      real(dp), dimension(:), allocatable   :: A_ij         
      ! Fermi Energy for Graphene
      real(dp), dimension(:), allocatable   :: fermi_energy 
      ! density for wFQ/wFQFMu
      real(dp), dimension(:), allocatable   :: density      
      ! Fermi d parameter
      real(dp), dimension(:,:), allocatable :: fermi_d      
      ! Fermi s parameter
      real(dp), dimension(:,:), allocatable :: fermi_s      
      !dynamic polar [wfqfmu]
      complex(dp), dimension(:,:), allocatable :: alpha_w
      ! check defined interactions
      logical, dimension(:,:), allocatable  :: interaction   
   end type parameters_type
    
   type (parameters_type), save :: parameters

contains

   !> Subroutine for constructing Ag permittivity or interband polarizability
   !! based on Etchegoin parametrization.
   !! P. B. Johnson and R.-W. Christy, Phys. Rev. B 6, 4370 (1972)
   !!    Input  : iunit        -- output unit
   !!    Input  : what         -- epsilon (bem) or polar(wfqfmu)
   !!    Input  : n_freq_in    -- number of frequencies given
   !!    Input  : freq_in      -- frequencies given
   !!    Output : vector       -- permittivity or interband polar
   subroutine silver_etchegoin(iunit,what,n_freq_in,freq_in,vector)

      implicit none
       
      !input variables
      integer, intent(in)                            :: iunit   !output unit
      character(len=*), intent(in)                   :: what    !eps or polar
      integer, intent(in)                            :: n_freq_in !num freq
      real(dp), dimension(n_freq_in), intent(in)     :: freq_in !frequencies
      complex(dp), dimension(n_freq_in), intent(out) :: vector  !eps (bem) or
                                                                !alpha(wfqfmu)
      !internal variables
      integer                    :: i
      real(dp)                   :: freqnm
      real(dp)                   :: w_p, n_0
      real(dp)                   :: lambda, gammap, freq3
      real(dp), dimension(3)     :: aj, lamj, gamj
      complex(dp)                :: bsline
      complex(dp), dimension (3) :: pij
      real(dp), parameter        :: Ag_density = 10.49d0    !in g/cm**3
      real(dp), parameter        :: N_avogadro = 6.022d23   !in 1/mol
      real(dp), parameter        :: MM_Ag      = 107.8682d0 !in g/mol
  
      !fitted parameters 
      bsline = dcmplx(1.469322d0,0.2450325d0)
      aj(1)  = 0.1632212d0
      aj(2)  = 0.9672775d0
      aj(3)  = 0.6150938d0
      pij(1) = (zero, -1.709221d0)
      pij(2) = (zero,  0.1546664d0)
      pij(3) = (zero, -1.782307d0)

      lambda  = 137.9643d0
      gammap  = one/74484.21d0
      lamj(1) = one/307.9550d0
      lamj(2) = one/260.7301d0
      lamj(3) = one/236.2128d0
      gamj(1) = one/4277.721d0
      gamj(2) = one/1051.474d0
      gamj(3) = one/1008.622d0
  
      n_0 = Ag_density/MM_Ag*N_avogadro/((1.0d8*ToBohr)**3) !in au
      w_p = dsqrt(four*pi*n_0) ! plasma frequency
  
      if(what.eq.'epsilon') then
         do i = 1, n_freq_in
            call freqautonm(freq_in(i),freqnm)
            !freq3 = one/freq_in(i) !1/nm
            freq3 = one/freqnm
            !calcolo epsilon
            vector(i) = bsline-1.0d0/(lambda**2*freq3*dcmplx(freq3,gammap)) + &
                        aj(1)*lamj(1)*(exp(pij(1)) /                          &
                        dcmplx(lamj(1)-freq3,-gamj(1)) +                      &
                        exp(-pij(1))/dcmplx(lamj(1)+freq3,gamj(1))) +         &
                        aj(2)*lamj(2)*(exp(pij(2)) /                          &
                        dcmplx(lamj(2)-freq3,-gamj(2)) +                      &
                        exp(-pij(2))/dcmplx(lamj(2)+freq3, gamj(2))) +        &
                        aj(3)*lamj(3)*(exp(pij(3)) /                          &
                        dcmplx( lamj(3)-freq3,-gamj(3)) +                     &
                        exp(-pij(3))/dcmplx(lamj(3)+freq3, gamj(3)))
         enddo
      else if(what.eq.'polar') then
         do i = 1, n_freq_in
            call freqautonm(freq_in(i),freqnm)
            !freq3 = one/freq_in(i) !1/nm
            freq3 = one/freqnm
            !calcolo polar
            vector(i) = ( bsline-dble(bsline) +                           &
                          aj(1)*lamj(1)*(exp(pij(1)) /                    &
                          dcmplx(lamj(1)-freq3,-gamj(1)) +                &
                          exp(-pij(1))/dcmplx(lamj(1)+freq3, gamj(1))) +  &
                          aj(2)*lamj(2)*(exp(pij(2)) /                    &
                          dcmplx( lamj(2)-freq3,-gamj(2))  +              &
                          exp(-pij(2))/dcmplx(lamj(2)+freq3, gamj(2))) +  &
                          aj(3)*lamj(3)*(exp(pij(3)) /                    &
                          dcmplx(lamj(3)-freq3,-gamj(3)) +                &
                          exp(-pij(3))/dcmplx(lamj(3)+freq3,gamj(3))) ) / &
                        w_p**2
         enddo
      else
         write(iunit,'(/1x,a)') "Error during the execution of plasmonX"
         write(iunit,'(1x,a/)') "What: "//trim(what)//" not recognised in &
                                &silver Etchegoin"
         stop
      endif
  
   end subroutine silver_etchegoin


   !> Subroutine for constructing Ag permittivity 
   !! based on Brendel Bormann parametrization.
   !! Appl. Opt. 37, 5271 (1998)
   !!    Input  : iunit        -- output unit
   !!    Input  : what         -- epsilon (bem) 
   !!    Input  : n_freq_in    -- number of frequencies given
   !!    Input  : freq_in      -- frequencies given
   !!    Output : vector       -- permittivity 
   subroutine silver_brendel_bormann(iunit, what, n_freq_in, freq_in, vector)
   
      implicit none
   
      !input variables
      integer, intent(in)                            :: iunit   !output unit
      character(len=*), intent(in)                   :: what    !eps or polar
      integer, intent(in)                            :: n_freq_in !num freq
      real(dp), dimension(n_freq_in), intent(in)     :: freq_in !frequencies
      complex(dp), dimension(n_freq_in), intent(out) :: vector  !eps (bem)
   
      !internal variables
      integer :: i, j
      real(dp) :: w !freq eV
      complex(dp) :: eps_drude, bb_corr, a_val(5), eps_final
      complex(dp) :: z_m(5), z_p(5), u_m(5), u_p(5), u(5), chi_j(5)
      real(dp) :: G(5), wj(5), sj(5), fj(5), G0, f0, wp, omega_p, sqrt2

      write(iunit,'(/1x,a)') "Error during the execution of plasmonX"
      write(iunit,'(1x,a/)') "Silver Brendel Bormann permittivity NYI"
      stop
   
      ! AG-specific Brendel-Bormann parameters
      f0 = 0.821d0
      G0 = 0.049d0
      fj = (/ 0.050d0, 0.133d0, 0.051d0, 0.467d0, 4.000d0 /)
      G  = (/ 0.189d0, 0.067d0, 0.019d0, 0.117d0, 0.052d0 /)
      wj = (/ 2.025d0, 5.185d0, 4.343d0, 9.809d0, 18.56d0 /)
      sj = (/ 1.894d0, 0.665d0, 0.189d0, 1.170d0, 0.516d0 /)
      wp = 9.01d0
      omega_p = sqrt(f0) * wp
      sqrt2 = sqrt(2.0d0)

      if (trim(what).eq.'epsilon') then
         do i = 1, n_freq_in
            call freqautoev(freq_in(i), w)
   
            !Drude term
            eps_drude = (omega_p**two) / (w * (w + dcmplx(0.0d0, G0)))
   
            !a coefficients (complex)
            do j = 1, 5
              a_val(j) = dcmplx( w/sqrt2*sqrt(sqrt(one+(G(j)/w)**two)+one), &
                                 w/sqrt2*sqrt(sqrt(one+(G(j)/w)**two)-one))
   
              z_m(j) = (a_val(j)-dcmplx(wj(j),zero))/(dcmplx(sqrt2*sj(j),zero))
              z_p(j) = (a_val(j)+dcmplx(wj(j),zero))/(dcmplx(sqrt2*sj(j),zero))
   
              ! U function via Faddeeva (approximation or external lib required)
              ! Here we fake with a simple lorentzian-like model as placeholder
              ! we do not have wofz ... not yet implemented
              !u_m(j) = sqrt(pi) / wofz(z_p(j))
              !u_p(j) = sqrt(pi) / wofz(z_m(j))
   
              u(j) = u_m(j) + u_p(j)
   
              chi_j(j) = dcmplx(zero,(fj(j)*wp**2)) / &
                         (two * sqrt2 * sj(j) * a_val(j)) * u(j)
            end do
   
            bb_corr = dcmplx(zero, zero)
            do j = 1, 5
              bb_corr = bb_corr + chi_j(j)
            end do
   
            eps_final = dcmplx(one, zero) - eps_drude + bb_corr
   
            vector(i) = eps_final
         end do
      else
         write(iunit,'(/1x,a)') "Error during the execution of plasmonX"
         write(iunit,'(1x,a/)') "What: "//trim(what)//" not recognised in silver &
                                &Brendel Bormann"
         stop
      end if

   end subroutine silver_brendel_bormann


   !> Subroutine for constructing Au permittivity or interband polarizability
   !! based on Etchegoin parametrization. 
   !! P. B. Johnson and R.-W. Christy, Phys. Rev. B 6, 4370 (1972)
   !!    Input  : iunit        -- output unit
   !!    Input  : what         -- epsilon (bem) or polar(wfqfmu)
   !!    Input  : n_freq_in    -- number of frequencies given
   !!    Input  : freq_in      -- frequencies given
   !!    Output : vector       -- permittivity or interband polar
   subroutine gold_etchegoin(iunit,what,n_freq_in,freq_in,vector)
      
      implicit none
       
      !input variables
      integer, intent(in)                            :: iunit   !output unit
      character(len=*), intent(in)                   :: what    !eps or polar
      integer, intent(in)                            :: n_freq_in !num freq
      real(dp), dimension(n_freq_in), intent(in)     :: freq_in !frequencies
      complex(dp), dimension(n_freq_in), intent(out) :: vector  !eps (bem) or
                                                                !alpha(wfqfmu)
      !internal variables
      integer                         :: i
      real(dp)                        :: freqnm
      real(dp)                        :: w_p, n_0
      real(dp)                        :: lambda, gammap, freq3
      real(dp), dimension(2)          :: aj, lamj, gamj
      complex(dp)                     :: bsline
      complex(dp), dimension (2)      :: pij
      real(dp), parameter             :: Au_density = 19.320d0    !in g/cm**3
      real(dp), parameter             :: N_avogadro = 6.022d23    !in 1/mol
      real(dp), parameter             :: MM_Au      = 196.96657d0 !in g/mol
  
      !fitted parameters
      bsline = dcmplx(1.54d0 ,zero )
      aj(1)  = 1.27d0  !1.27 JC
      aj(2)  = 1.10d0  !1.1
      pij(1) = dcmplx(zero, - pi / four)
      pij(2) = dcmplx(zero, - pi / four)
       
      !the following data are in nm
      lambda = 143.0d0
      gammap = one/14500.0d0
      lamj(1)= one/  470.0d0
      lamj(2)= one/  325.0d0
      gamj(1)= one/ 1900.0d0
      gamj(2)= one/ 1060.0d0
  
      n_0 = Au_density/MM_Au*N_avogadro/((1.0d8*ToBohr)**3) !in a.u.
      w_p = dsqrt(four*pi*n_0) !plasma frequency
  
      if(what.eq.'epsilon') then
         do i = 1, n_freq_in
            call freqautonm(freq_in(i),freqnm)
            !freq3 = one/freq_in(i) !1/nm
            freq3 = one/freqnm
            ! epsilon for BEM
            vector(i) = bsline-1.0d0/(lambda**2*freq3*dcmplx(freq3,gammap)) + &
                        aj(1)*lamj(1)*(exp(pij(1)) /                          &
                        dcmplx(lamj(1)-freq3,-gamj(1)) +                      &
                        exp(-pij(1))/dcmplx(lamj(1)+freq3,gamj(1))) +         &
                        aj(2)*lamj(2)*(exp(pij(2)) /                          &
                        dcmplx(lamj(2)-freq3,-gamj(2)) +                      &
                        exp(-pij(2))/dcmplx(lamj(2)+freq3,gamj(2)))                  
         enddo
      else if(what.eq.'polar') then
         do i = 1, n_freq_in
            call freqautonm(freq_in(i),freqnm)
            !freq3 = one/freq_in(i) !1/nm
            freq3 = one/freqnm
            ! calcolo polar for wfqfmu
            vector(i) = ( aj(1)*lamj(1)*(exp(pij(1)) /                    &
                          dcmplx(lamj(1)-freq3,-gamj(1)) +                &
                          exp(-pij(1))/dcmplx(lamj(1)+freq3,gamj(1))) +   &
                          aj(2)*lamj(2)*(exp(pij(2))/                     &
                          dcmplx(lamj(2)-freq3,-gamj(2)) +                &
                          exp(-pij(2))/dcmplx(lamj(2)+freq3,gamj(2))) ) / &
                          w_p**2
            ! correct if a negative polar is obtained
            if(dimag(vector(i)) .le. 0.0d0) &
               vector(i) = dcmplx(dble(vector(i)),1.0d-10)
         enddo
      else
         write(iunit,'(/1x,a)') "Error during the execution of plasmonX"
         write(iunit,'(1x,a/)') "What: "//trim(what)//" not recognised in &
                                &gold Etchegoin"
         stop
      endif
    
   end subroutine gold_etchegoin


   !> Subroutine for reading permittivity info (epsilon or polar) from file
   !!    Input  : iunit        -- output unit
   !!    Input  : file_        -- file to be read
   !!    Input  : n_freq_in    -- number of frequencies given
   !!    Input  : freq_in      -- frequencies given
   !!    Output : vector       -- permittivity or interband polar
   subroutine read_permittivity_info_from_file(iunit,file_,n_freq_in,freq_in, &
                                               vector)
      use string_manipulation_module

      implicit none

      !input variables
      integer, intent(in)                            :: iunit   !output unit
      character(len=*), intent(in)                   :: file_   !file to be read
      integer, intent(in)                            :: n_freq_in !num freq
      real(dp), dimension(n_freq_in), intent(in)     :: freq_in !frequencies
      complex(dp), dimension(n_freq_in), intent(out) :: vector  !eps (bem) or
                                                                !alpha(wfqfmu)
      !internal variables
      integer :: num_lines_commented
      integer :: num_freq
      integer :: i,j
      character(len=200) :: line_file
      logical :: exists
      real(dp) :: freq_eV
      real(dp), dimension(:), allocatable :: freq
      real(dp), dimension(:), allocatable :: tmp_re
      real(dp), dimension(:), allocatable :: tmp_im
  
      inquire(file=file_,exist=exists) 
      if(.not.exists) then 
         write(iunit,'(/1x,a)') "Error during the execution of plasmonX"
         write(iunit,'(1x,a/)') "File: "//trim(file_)//" does not exist. &
                                &Check the name"
         stop
      endif
          
      open(unit=11, file=file_, status="old") !open file 
         rewind(11) 
         read(11,"(a)") line_file
         num_lines_commented = 0 
         !skypping any comments in the file: they should be at the beginning
         do while (index(line_file,"#").gt.0)
            read(11,"(a)") line_file
            num_lines_commented = num_lines_commented + 1
         enddo
         call get_number_lines(11,num_freq)
         num_freq = num_freq + 1
          
         rewind(11) 
         do i = 1, num_lines_commented !skip commented lines
            read(11,"(a)") line_file
         enddo

         allocate(freq(num_freq)) 
         allocate(tmp_re(num_freq))
         allocate(tmp_im(num_freq))

         do i=1, num_freq
            ! format of file : freq (eV)    Re    Im (epsilon or polar)
            read(11,*) freq(i), tmp_re(i), tmp_im(i)
         enddo 
         !adjusting the frequencies in external file with requested frequencies
         do i=1, n_freq_in 
            call freqautoev(freq_in(i), freq_eV)
            do j = 1, num_freq ! num freq file
               if(abs(freq_eV-freq(j)).le.1.0d-08) then
                  vector(i) = dcmplx(tmp_re(j), tmp_im(j))
                  exit
               else if(j.eq.num_freq) then
                  write(iunit,'(/1x,a)') "Error during the execution of plasmonX"
                  write(iunit,'(1x,a,f6.2,/)') &
                    "I did not find data for freq = ", freq_eV
                  stop
               endif
            enddo !num frequency file
         enddo!num frequency in input
          
         deallocate(freq)
         deallocate(tmp_re)
         deallocate(tmp_im)
      close(unit=11) !close file
  
   end subroutine read_permittivity_info_from_file


   !> Subroutine for converting frequency from eV to nm
   !!    Input  : freqeV -- in eV
   !!    Output : freqnm -- in nm
   subroutine freqevtonm(freqeV,freqnm)
 
      implicit none
       
      !input variables
      real(dp), intent(in)  :: freqeV  ! in eV 
      real(dp), intent(out) :: freqnm  ! in nm
 
      !internal variables
      real(dp) :: freqau 
      real(dp) :: freqcm
      
      freqau = freqeV / au_to_ev
      freqcm = freqau * au_to_cm
      freqnm  = 10000000.0d0/ freqcm
    
   end subroutine freqevtonm


   !> Subroutine for converting frequency from ev to au
   !!    Input  : freqev -- in eV
   !!    Output : freqau -- in au
   subroutine freqevtoau(freqev,freqau)
  
      implicit none
       
      !input variables
      real(dp), intent(in)  :: freqev  !in ev
      real(dp), intent(out) :: freqau  !in au
  
      freqau = freqev/au_to_ev
    
   end subroutine freqevtoau


   !> Subroutine for converting frequency from ev to au
   !!    Input  : freqau -- in au
   !!    Output : freqev -- in eV
   subroutine freqautoev(freqau,freqev)
  
      implicit none
       
      !input variables
      real(dp), intent(in)  :: freqau  !in au
      real(dp), intent(out) :: freqev  !in ev
  
      freqev = freqau*au_to_ev
    
   end subroutine freqautoev


   !> Subroutine for converting frequency from nm to au
   !!    Input  : freqnm -- in nm
   !!    Output : freqau -- in au
   subroutine freqnmtoau(freqnm,freqau)
  
      implicit none
       
      !input variables
      real(dp), intent(in)  :: freqnm  !in nm
      real(dp), intent(out) :: freqau  !in au
  
      !internal variables
      real(dp) :: freqcm
       
      freqcm = 10000000.0d0/freqnm
      freqau = freqcm / au_to_cm
    
   end subroutine freqnmtoau


   !> Subroutine for converting frequency from nm to au
   !!    Input  : freqau -- in au
   !!    Output : freqnm -- in nm
   subroutine freqautonm(freqau,freqnm)
  
      implicit none
       
      !input variables
      real(dp), intent(in)  :: freqau  !in au
      real(dp), intent(out) :: freqnm  !in nm
  
      !internal variables
      real(dp) :: freqcm
       
      freqcm = freqau * au_to_cm
      freqnm = 10000000.0d0 / freqcm 
    
   end subroutine freqautonm


   !> Subroutine for converting frequency from nm to eV
   !!    Input  : freqnm -- in nm
   !!    Output : freqeV -- in eV
   subroutine freqnmtoev(freqnm,freqeV)
  
      implicit none
       
      !input variables
      real(dp) :: freqnm  ! in nm
      real(dp) :: freqeV  ! in ev
  
      !internal variables
      real(dp) :: freqau 
      
      call freqnmtoau(freqnm,freqau)
      freqeV = freqau * au_to_ev
     
   end subroutine freqnmtoev

end module parameters_module
