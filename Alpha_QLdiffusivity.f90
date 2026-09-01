!---------------------------------------------------------
! Alpha_QLdiffusivity.f90
!
! PURPOSE:
! 
!   compute the the diffusivity profile for two EP specied 
!   OUTPUT:   QLdiffEP(:,:)  [m**2/sec]  
!    from
!   INPUT:    den(:,:) transported density  [10**19 1/m**3] 
!   
!    calld by general time dependent transport 
!     Alpha_time_dep_transport.f90
!
!   notation in Alpha_transport.f90
!
!  like Alpha_diffusivity which is designed for steady state 
!       ITG/TEM  background EP diffusion with stiff CGM AE EP transport
!  except
!       Alpha_QLdiffusivity is designed for intermittent 
!       Quasi-Linear CGM tranport EP based on TGLFEP CG and QL rates
!
!       There is a low level ITG/TEM  background EP diffusion D_bkg
!        It is possble call Alpha_diffusivity with the CGM off
!        (set i_tot_TAE=-100) then call Alpha_QLdiffusivity
!        for the QL EP model

!---------------------------------------------------------

subroutine Alpha_QLdiffusivity(den,QLdiffEP,dt_update,i_test_print,i_sav_QLdiffEP,i_bkg)

    use Alpha_use_input

    use Alpha_use_output

    use Alpha_use_other

    use Alpha_use_sav_diffEP
    use Alpha_use_sav_QLdiffEP  !communication only in Alpha_QLdiffusivity.f90


  !--------------------------------------
  implicit none
  
     integer :: i  !i=1,n_rho_grid)   rho_hat(1) = 0.   rho_hat(n_rho_grid) =1.0
     integer :: i_s !species index   i_s=1,i_s_max
     integer :: i_s_max
!     i_s_max=1 NBI_flag=0 or 1,
!     i_s_max=2 NBI_flag=2 

     integer :: km  !mode label
     integer :: km_max ! max mode label .le. 15

     integer :: i_read_TGLFEP
     integer :: i_one_spec_test


     integer :: i_threshold
     integer :: i_threshold2

     real, intent(in) :: dt_update  !should be same as time_run steps
     integer, intent(in) :: i_test_print  !for test on first call
     integer, intent(in) :: i_sav_QLdiffEP  !=1 yes save  =0 already saved
     integer, intent(in) :: i_bkg   !=1 yes add D_bkg  =0 no D_bkg



     real,  dimension(2,n_rho_grid), intent(in) :: den  !n_alpha_tran_rho, n_alpha2_tran_rho
     real,  dimension(2,n_rho_grid), intent(out) :: QLdiffEP    !EP diffusivity  [m**2/sec]


!general working

!slowing down
     real, dimension(n_rho_grid) :: rg_n_alpha_rho
     real, dimension(n_rho_grid) :: rg_n_alpha2_rho
     real, dimension(n_rho_grid) :: rg_T_alpha_rho
     real, dimension(n_rho_grid) :: rg_T_alpha2_rho

     real, dimension(n_rho_grid) :: Ln_alpha_rho
     real, dimension(n_rho_grid) :: Ln_alpha2_rho
     real, dimension(n_rho_grid) :: LT_alpha_rho
     real, dimension(n_rho_grid) :: LT_alpha2_rho


!time dependent transported     
     real, dimension(n_rho_grid) :: n_alpha_tran_rho
     real, dimension(n_rho_grid) :: n_alpha2_tran_rho

     real, dimension(n_rho_grid) :: rg_n_alpha_tran_rho
     real, dimension(n_rho_grid) :: rg_n_alpha2_tran_rho


     real, dimension(n_rho_grid) :: rg_n_alpha_th_rho
     real, dimension(n_rho_grid) :: rg_n_alpha2_th_rho

     real, dimension(n_rho_grid) :: Ln_alpha_tran_rho
     real, dimension(n_rho_grid) :: Ln_alpha2_tran_rho

     real :: dr
     real :: pi

     real :: C_nl  !   C_nl*E_hat  nonlinear damping rate
     real :: CZ_nl !   CZ_nl*Z_hat  Z_hat nonlinear transfer rate

     real :: D_bkg
     real, dimension(n_rho_grid) :: D_alpha
     real, dimension(n_rho_grid) :: D_alpha2

!RQL 
    real, dimension(2,15,n_rho_grid)  :: gamma_star_hat
     !gB units
    real, dimension(2,15,n_rho_grid)  :: diff_star
     ! [m**2/sec]
    real, dimension(2,15,n_rho_grid)  :: rg_n_crit_star
     !same units as  rg_n_alpha_th_rho,  rg_n_alpha2_th_rho

    real, dimension(15,n_rho_grid)  :: gamma_AE_hat !total alpha & alpha2
    real, dimension(15,n_rho_grid)  :: gamma_Z_hat
    real, dimension(15,n_rho_grid)  :: Z_hat_prev

    real, dimension(2,15,n_rho_grid)  :: gamma_hat  !each sp acting alone
     !gB_units

    real :: F_frac
    real :: F_frac_one_spec_test
    real, dimension(2,n_rho_grid) ::  f_cor
    integer :: i_f_cor
    real :: C_R
    real, dimension(15) :: F_mode_drive

    integer :: i_del2
    integer :: i_del3
    integer :: i_del4
    integer :: i_del5
  
    real :: frac_crit1
    real :: frac_crit2
    real :: frac_crit3
    real :: frac_crit4
    real :: frac_crit5

    
    integer :: i_RBF   !1 = on  PRBF; 0 = off  resonance bradening
                       !2 = on  NRBF; 0 = off  resonance bradening
    real :: RBF(15,n_rho_grid)  !resonance broadening factor
    real :: omega_RBF

    integer :: i_RBF1  !1 = on  PRBF1; 0 = off  resonance bradening
    real :: RBF1(2,15,n_rho_grid)  !resonance broadening factor each species
    real :: omega_RBF1

!  set up grid
!  real, dimension(n_rho_grid) :: rho_hat

   dr = rmin/float(n_rho_grid-1)  ![m]

   pi = 3.141592565

  do i = 1,n_rho_grid
   rho_hat(i) = float(i-1)/float(n_rho_grid-1)
  enddo

  i_threshold  = 10
  i_threshold2 = 10

  i_one_spec_test = 0 !0 = 2 specs ;  1= 1 spec #1 alpha 
                      !2 = 1 spec #2 alpha2

  i_one_spec_test = 1
 
  if(NBI_flag .eq. 0 .or. NBI_flag .eq. 1) i_one_spec_test =1

  if (i_one_spec_test .eq. 0) F_frac_one_spec_test = 0.5 !#1 and #2
  if (i_one_spec_test .eq. 1) F_frac_one_spec_test = 1.0 !#1 only
  if (i_one_spec_test .eq. 2) F_frac_one_spec_test = 0.0 !#2 only

!   general setup
!   --------------------------------------------------------------

!  setup critical gradients and slowing down gradients

  if(i_sav_QLdiffEP .eq. 1) then

! slowing down alpha
  rg_n_alpha_rho(1) = 0.
  do i = 2,n_rho_grid-1
    rg_n_alpha_rho(i)= &
   -(n_alpha_rho(i+1)-n_alpha_rho(i-1))/rmin/(rho_hat(i+1)-rho_hat(i-1))
  enddo
   i=n_rho_grid
    rg_n_alpha_rho(i)= &
    -(n_alpha_rho(i)-n_alpha_rho(i-1))/rmin/(rho_hat(i)-rho_hat(i-1))

  rg_T_alpha_rho(1) = 0.
  do i = 2,n_rho_grid-1
    rg_T_alpha_rho(i)= &
   -(T_alpha_equiv_rho(i+1)-T_alpha_equiv_rho(i-1))/rmin/(rho_hat(i+1)-rho_hat(i-1))
  enddo
   i=n_rho_grid
    rg_T_alpha_rho(i)= &
    -(T_alpha_equiv_rho(i)-T_alpha_equiv_rho(i-1))/rmin/(rho_hat(i)-rho_hat(i-1))

! Ln_sd & LT_sd for alpha
   do i = 2,n_rho_grid
    Ln_alpha_rho(i) = n_alpha_rho(i)/rg_n_alpha_rho(i)
    LT_alpha_rho(i) = T_alpha_equiv_rho(i)/rg_T_alpha_rho(i)
   enddo
    Ln_alpha_rho(1) = Ln_alpha_rho(2)
    LT_alpha_rho(1) = LT_alpha_rho(2)
  

  if(i_threshold .eq. 10) then
   rg_n_alpha_th_rho(:) = crit_grad_n_alpha_rho(:)   !TGLF input
  endif


    sav_rg_n_alpha_th_rho(:) = rg_n_alpha_th_rho(:)

    sav_rg_n_alpha_rho(:) = rg_n_alpha_rho(:)
    sav_Ln_alpha_rho(:) = Ln_alpha_rho(:)
    sav_LT_alpha_rho(:) = LT_alpha_rho(:)
   endif ! i_sav_QLdiffEP .eq. 1


   if(NBI_flag .eq. 2) then

   if(i_sav_QLdiffEP .eq. 1) then

! slowing down for alpha2
  rg_n_alpha2_rho(1) = 0.
  do i = 2,n_rho_grid-1
    rg_n_alpha2_rho(i)= &
   -(n_alpha2_rho(i+1)-n_alpha2_rho(i-1))/rmin/(rho_hat(i+1)-rho_hat(i-1))
  enddo
   i=n_rho_grid
    rg_n_alpha2_rho(i)= &
    -(n_alpha2_rho(i)-n_alpha2_rho(i-1))/rmin/(rho_hat(i)-rho_hat(i-1))

  rg_T_alpha2_rho(1) = 0.
  do i = 2,n_rho_grid-1
    rg_T_alpha2_rho(i)= &
   -(T_alpha2_equiv_rho(i+1)-T_alpha2_equiv_rho(i-1))/rmin/(rho_hat(i+1)-rho_hat(i-1))
  enddo
   i=n_rho_grid
    rg_T_alpha2_rho(i)= &
    -(T_alpha2_equiv_rho(i)-T_alpha2_equiv_rho(i-1))/rmin/(rho_hat(i)-rho_hat(i-1))

! Ln_sd & LT_sd for alpha2
   do i = 2,n_rho_grid
    Ln_alpha2_rho(i) = n_alpha2_rho(i)/rg_n_alpha2_rho(i)
    LT_alpha2_rho(i) = T_alpha_equiv_rho(i)/rg_T_alpha2_rho(i)
   enddo
    Ln_alpha2_rho(1) = Ln_alpha2_rho(2)
    LT_alpha2_rho(1) = LT_alpha2_rho(2)

   if(i_threshold2 .eq. 10) then
    rg_n_alpha2_th_rho(:) = crit_grad_n_alpha2_rho(:)   !TGLF input
   endif

   sav_rg_n_alpha2_rho(:) = rg_n_alpha2_rho(:)
   sav_Ln_alpha2_rho(:) = Ln_alpha2_rho(:)
   sav_LT_alpha2_rho(:) = LT_alpha2_rho(:)
   
  endif !i_sav_QLdiffEP .eq. 1

  endif !NBI_flag = 2

!end critical gradient setup


!   time_run dependent
!   transported densities

    n_alpha_tran_rho(:) = den(1,:)    ![10**19/m**3]
    n_alpha2_tran_rho(:) = den(2,:)   ![10**19/m**3]

!   radial gradients  for transported density profile
!     rg_n_alpha_tran_rho(i), rg_n_alphai2_tran_rho(i)   -dn/dr [10**19/m**3]/m 

!alpha
  rg_n_alpha_tran_rho(1) = 0.
  do i = 2,n_rho_grid-1
    rg_n_alpha_tran_rho(i)= &
   -(n_alpha_tran_rho(i+1)-n_alpha_tran_rho(i-1))/rmin/(rho_hat(i+1)-rho_hat(i-1))
  enddo
   i=n_rho_grid
    rg_n_alpha_tran_rho(i)= &
    -(n_alpha_tran_rho(i)-n_alpha_tran_rho(i-1))/rmin/(rho_hat(i)-rho_hat(i-1))

   do i=2,n_rho_grid
    Ln_alpha_tran_rho(i) = n_alpha_tran_rho(i)/rg_n_alpha_tran_rho(i)
   enddo
    Ln_alpha_tran_rho(1) = Ln_alpha_tran_rho(2)


  if(NBI_flag .eq. 2) then
!alpha2
   rg_n_alpha2_tran_rho(1) = 0.
  do i = 2,n_rho_grid-1
    rg_n_alpha2_tran_rho(i)= &
     -(n_alpha2_tran_rho(i+1)-n_alpha2_tran_rho(i-1))/rmin/(rho_hat(i+1)-rho_hat(i-1))
  enddo
   i=n_rho_grid
    rg_n_alpha2_tran_rho(i)= &
     -(n_alpha2_tran_rho(i)-n_alpha2_tran_rho(i-1))/rmin/(rho_hat(i)-rho_hat(i-1))

   do i=2,n_rho_grid
    Ln_alpha2_tran_rho(i) = n_alpha2_tran_rho(i)/rg_n_alpha2_tran_rho(i)
   enddo
    Ln_alpha2_tran_rho(1) = Ln_alpha2_tran_rho(2)

    endif  !NBI_flag = 2

! default
!    QLdiffEP(:,:) = 1.0
!    return
  
!start models
   QLdiffEP(:,:) = 0.

! some HARDWIRES

   i_read_TGLFEP = 0 
!   i_f_cor = 0
   i_f_cor = 1! 1.22.18  1:39pm
   C_R=2.0

!   i_read_TGLFEP = 1 ! read "star's" from TGLFEP files
!   km_max = 1
    km_max = 5
!   km_max = 2


   D_bkg = 0.0
  if(i_bkg .eq. 1) then
   D_bkg = 0.001
!   D_bkg = 0.003
!   D_bkg = 0.01
!   D_bkg = 0.05
  endif

  D_alpha(:) = D_bkg
  D_alpha2(:) =  D_bkg        !progamming oly for alpha (1 species) QL

   C_nl = 0.0   !normal QL
!  C_nl =  0.5  !nonlinear damping  C_nl*E_hat 
  CZ_nl = 0.0  !normal QL
!   CZ_nl = 1.0 !nonlinear transfer CZ_nl*Z_hat transfer rate
!   CZ_nl = 0.5
!   CZ_nl = 3.
!   CZ_nl = 15.
!   CZ_nl = 50.
!   CZ_nl = 250.
!    CZ_nl = 100.
!  set up for i_one_spec_test =1 first species "alpha"  only

  i_RBF = 0 !normal  i.e. no RBF
!  i_RBF = 1 !PRBF
!   i_RBF = 1  !abs(gamma)/(gamma**2) + omega**2)
!   i_RBF = 2 !NRBF  !RBF1 = 0.0 if .lt. 0 
   omega_RBF = 0.1
   omega_RBF = 1.0  !1.19.18  11:03am

  i_RBF1 = 0 !normal
!  i_RBF1 = 1 ! replace  RBF1 with RBF  !NaN problem
  omega_RBF1 = 0.1
  omega_RBF1 = 1.0  !1.19.18  11:03am

    F_mode_drive(:) = 1.0

! test multi-mode drive 
!   F_mode_drive(1) = 0.8
!   F_mode_drive(2) = 0.9
!   F_mode_drive(3) = 1.0
!   F_mode_drive(4) = 1.1
!   F_mode_drive(5) = 1.2

  if(i_test_print .eq. 1) then
    print *, 'km,  F_mode_drive(km)'
   do km=1,km_max
    print *, km, ' ',  F_mode_drive(km)
   enddo
  endif
 

!NOTE only one species QLdiffusity
!  if(i_test_print .eq. 1)    print *, 'only one species QLdiffusity'
!    QLdiffEP(2,:) = D_alpha2(:)
!NOTE only one species QLdiffusity

!   D_W(is,i_mode,i) 
!   is=1,is_max<2, i_mode = 1,i_mode_max<5, i=1,n_rho_grid (nn=51) 
!   see Alpha_use_sav_QLdiffEP  HARDWIRES  
!
!   per mode k
!   D_W_k = (dFlux_W_k/dn_sf)/(-d n_sd/dr) is the QL diffusivity weight
!   Flux_W_k is the QL flux
!   n_sf is the scale factor on the profile of n_EP = n_sd at n_sf =1

!   (dFlux_W_k/dn_sf) is obtained by running TGLFEP with thermal gradients 
!      off, for say  n_sf = 0.5, 1.0, 1.5 or whatever needed
!       the a/Ln_EP profile = a/Ln_EP_sd  is unchanged 

!       must multiply TGLF Flux_W_k by the local gyroBohm unit of particle
!       flux so that D_W_k appears in m**2/sec

!       key assumption is that Flux_W_k is linear in n_sf 

!       

!GOING FORWARD as a test, a  D_W_k 

!---------------------------------------------------------------------
   if(i_sav_QLdiffEP .eq. 1) then !start up at time_run = 0.

     E_hat(:,:) = 0.01  !time_run = 0 star one time call
     Z_hat(:,:) = 0.01  !if Z_hat started at 0 it will stay zero   
!Note:  Could be a "posivity" problem if the "intensity" Z_hat goes negative
     time_called = 0.0

   print *, 'i_read_TGLFEP=',i_read_TGLFEP
   print *, 'i_one_spec_test=',i_one_spec_test
   print *, 'i_f_cor=',i_f_cor
   print *, 'C_R=',C_R
   print *, 'km_max=',km_max
   print *, 'i_bkg=',i_bkg
   print *, 'D_bkg=',D_bkg
   print *, 'E_hat(1,1)=',E_hat(1,1)
   print *, 'C_nl=',C_nl
   print *, 'CZ_nl=',CZ_nl
   print *, 'i_RBF=',i_RBF
   print *, 'i_RBF=1 has abs(gamma)/(gamma**2 + omega_RBF**2)'
   print *, 'omega_RBF=',omega_RBF
   print *, 'i_RBF1=',i_RBF1
   print *, 'omega_RBF1=',omega_RBF1

       !RQL values: gamma_star_hat, diff_star, rg_n_crit_star 
       gamma_star_hat(:,:,:) = 0.0
       diff_star(:,:,:) = 0.0
       rg_n_crit_star(:,:,:) = 0.0

! one km mode test model
    if(i_read_TGLFEP .eq. 0) then !make RQL values

!      gamma_star_hat (1,:,:) = 0.1*1.0 !most unstable
!      gamma_star_hat (2,:,:) = 0.1*1.0 !most unstable
     gamma_star_hat (:,:,:) = 0.1*1.0  !1.25.18pm

!--------------------------------------------
!      gamma_star_hat(1,:,:) = 1.0 ! [m**2/sec]
!      gamma_star_hat(2,:,:) = 1.0
!     if(i_test_print .eq. 1) then
!      print *, 'gamma_star_hat 0.1 --> 1.0 TEST'
!     endif
!------------------------------------------

!12.12.17pm  1/e, 1, 1/e
!   gamma_star profile factors
     if(i_test_print .eq. 1) then
      print *, 'gamma_star profile factor OFF'
     endif
!       do i=1,n_rho_grid
!        gamma_star_hat (:,1,i) = gamma_star_hat (:,1,i)*exp( &
!           -(rho_hat(i)-rho_hat(25))**2/(rho_hat(25)**2))
!       enddo

!      diff_star(1,:,:) = 1.0 ! [m**2/sec]
!      diff_star(2,:,:) = 1.0
      diff_star(:,:,:) = 1.0  !1.25.18

!--------------------------------------------
!      diff_star(1,:,:) = 0.1 ! [m**2/sec]
!      diff_star(2,:,:) = 0.1
!     if(i_test_print .eq. 1) then
!      print *, 'diff_star 1.0 --> 0.1 TEST'
!     endif
!------------------------------------------

      
     frac_crit1 =1.0 !normal
!      frac_crit1  =2.2
     do km=1,km_max
      rg_n_crit_star(1,km,:) = frac_crit1*rg_n_alpha_th_rho(:)
     enddo !km

!start shifted "island" chain game for 1 species alpha 
!   mode 2  

!       i_del2 = 5
!      i_del2 = 10
      i_del2 = 20
!       i_del2=0


!  frac_crit2 = 0.3
   frac_crit2 = 1.0

      rg_n_crit_star(1,2,:) = frac_crit2*rg_n_alpha_th_rho(1)
      do i=1,n_rho_grid
       if(i-i_del2 .ge. 1) then 
        rg_n_crit_star(1,2,i) = frac_crit2*rg_n_alpha_th_rho(i-i_del2)
       endif
      enddo

   if(i_test_print .eq. 1) then
    print *, 'island chain drive'
    print *, 'frac_crit1=',frac_crit1
    print *, 'i_del2=',i_del2,' ','frac_crit2=',frac_crit2
 
    do i=1,n_rho_grid
     print *, i, ' ', rg_n_crit_star(1,1,i)/rg_n_alpha_rho(i), &
        '  ', rg_n_crit_star(1,2,i)/rg_n_alpha_rho(i)
    enddo
   endif
!add more islands between 1 and 2
   i_del3 = 5
   i_del4 =10
   i_del5 =15
   frac_crit3 = 1.0
   frac_crit4 = 1.0
   frac_crit5 = 1.0
   if(i_test_print .eq. 1 .and. km_max .gt. 2) then
    print *, 'add mode islands between 1 and 2'
    print *, 'km_max=',km_max
    print *, 'i_del3=',i_del3
    print *, 'i_del4=',i_del4
    print *, 'i_del5=',i_del5
   endif

   if(km_max .gt. 2) then

      rg_n_crit_star(1,3,:) = frac_crit3*rg_n_alpha_th_rho(1)
      do i=1,n_rho_grid
       if(i-i_del2 .ge. 1) then
        rg_n_crit_star(1,3,i) = frac_crit3*rg_n_alpha_th_rho(i-i_del3)
       endif
      enddo

      rg_n_crit_star(1,4,:) = frac_crit4*rg_n_alpha_th_rho(1)
      do i=1,n_rho_grid
       if(i-i_del4.ge. 1) then
        rg_n_crit_star(1,4,i) = frac_crit4*rg_n_alpha_th_rho(i-i_del4)
       endif
      enddo

      rg_n_crit_star(1,5,:) = frac_crit5*rg_n_alpha_th_rho(1)
      do i=1,n_rho_grid
       if(i-i_del5 .ge. 1) then
        rg_n_crit_star(1,5,i) = frac_crit5*rg_n_alpha_th_rho(i-i_del5)
       endif
      enddo
   endif
!end shifed "island" chain game

    if(NBI_flag .eq. 2) then
!01.02.18PM ERROR      rg_n_crit_star(2,1,:) = 1.0*rg_n_alpha_th_rho(:)
    do km=1,km_max
     rg_n_crit_star(2,km,:) = 1.0*rg_n_alpha2_th_rho(:)
    enddo !km
    endif
!end one km mode test model
!    NOTE: It is really hard to make up relative strength of modes

    endif !i_read_TGLFEP .eq. 0
!---------------------------------------------------------------------

!---------------------------------------------------------------------
    if(i_read_TGLFEP .eq. 1) then !read  RQL values from TGLFEP files
!fill-in later
    endif !i_read_TGLFEP .eq. 1
!---------------------------------------------------------------------

    sav_gamma_star_hat(:,:,:)=gamma_star_hat(:,:,:)
    sav_diff_star(:,:,:) = diff_star(:,:,:)
    sav_rg_n_crit_star(:,:,:) = rg_n_crit_star(:,:,:)

    if(i_one_spec_test .eq. 1) then
     sav_gamma_star_hat(2,:,:) = 0.0
     sav_diff_star(2,:,:) = 0.0
    endif
    if(i_one_spec_test .eq. 2) then
     sav_gamma_star_hat(1,:,:) = 0.0
     sav_diff_star(1,:,:) = 0.0
    endif

   if(i_test_print .eq. 1) then
    print *, 'alpha slowing down drive factor'
    print *, '-------------------------'
   do i=2,n_rho_grid
   print *, &
        (sav_rg_n_alpha_rho(i) - sav_rg_n_crit_star(1,1,i))/&
          sav_rg_n_alpha_rho(i)
   enddo
    print *, '-------------------------'

    print *, 't=0  alpha drive factor'
    print *, '-------------------------'
   do i=2,n_rho_grid
   print *, &
        (rg_n_alpha_tran_rho(i) - sav_rg_n_crit_star(1,1,i))/&
          sav_rg_n_alpha_rho(i)
   enddo
    print *, '-------------------------'

 if(NBI_flag .eq. 2) then
   print *, 'alpha2 slowing down drive factor'
    print *, '-------------------------'
   do i=2,n_rho_grid
   print *, &
        (sav_rg_n_alpha2_rho(i) - sav_rg_n_crit_star(2,1,i))/&
          sav_rg_n_alpha2_rho(i)
   enddo
    print *, '-------------------------'

    print *, 't=0  alpha2 drive factor'
    print *, '-------------------------'
   do i=2,n_rho_grid
   print *, &
        (rg_n_alpha2_tran_rho(i) - sav_rg_n_crit_star(2,1,i))/&
          sav_rg_n_alpha2_rho(i)
   enddo
    print *, '-------------------------'
  endif

   endif !i_test_print .eq. 1

   endif !i_sav_QLdiffEP .eq. 1
!---------------------------------------------------------------------
   if(i_sav_QLdiffEP .eq. 2 .and.  km_max .eq. 2) then
    if(i_test_print.eq. 1) then
      print *, 'final time shifted "island" chain game'
     do i=2,n_rho_grid-1
      print *, i, rg_n_alpha_tran_rho(i)/sav_rg_n_crit_star(1,1,i), &
                  rg_n_alpha_tran_rho(i)/sav_rg_n_crit_star(1,2,i)
     enddo
    endif    
   endif !i_sav_QLdiffEP .eq. 2
  if(i_sav_QLdiffEP .eq. 2 .and.  km_max .eq. 1) then
    if(i_test_print.eq. 1) then
      print *, 'final time one "island" '
     do i=2,n_rho_grid-1
      print *, i, rg_n_alpha_tran_rho(i)/sav_rg_n_crit_star(1,1,i)
     enddo
    endif
   endif !i_sav_QLdiffEP .eq. 2

!---------------------------------------------------------------------
    if(i_sav_QLdiffEP .eq. 0) then !do time step update

    f_cor(:,:) = 1.0
    if(i_f_cor .eq. 1) then !set-up f_cor correction factor
     f_cor(1,i) =  &
    (1.0 + Ln_alpha_tran_rho(i)/sav_LT_alpha_rho(i) &
                  - C_R*Ln_alpha_tran_rho(i)/Rmaj_rho(i))
     f_cor(1,i) = f_cor(1,i)/ &
    (1.0 + sav_Ln_alpha_rho(i)/sav_LT_alpha_rho(i) &
                 - C_R*sav_Ln_alpha_rho(i)/Rmaj_rho(i))

     if(NBI_flag .eq. 2) then
     f_cor(2,i) =  &
    (1.0 + Ln_alpha2_tran_rho(i)/sav_LT_alpha2_rho(i) &
                  - C_R*Ln_alpha2_tran_rho(i)/Rmaj_rho(i))
     f_cor(1,i) = f_cor(1,i)/ &
    (1.0 + sav_Ln_alpha2_rho(i)/sav_LT_alpha2_rho(i) &
                 - C_R*sav_Ln_alpha2_rho(i)/Rmaj_rho(i))
     endif !NBI_flag .eq. 2
    endif !i_f_cor .eq. 1


!time_run advance of mode intensities
  if(NBI_flag .eq. 1 .or. NBI_flag .eq. 0) then
! only one species
   do km=1,km_max
    do i=1,n_rho_grid
      gamma_AE_hat(km,i) = sav_gamma_star_hat(1,km,i)*f_cor(1,i) &
          *(F_mode_drive(km)*rg_n_alpha_tran_rho(i)  &
              -sav_rg_n_crit_star(1,km,i)) &
             /sav_rg_n_alpha_rho(i)

       gamma_Z_hat(km,i) = sav_gamma_star_hat(1,km,i)*f_cor(2,i) &
             *(sav_rg_n_crit_star(1,km,i)) &
             /sav_rg_n_alpha_rho(i)

      RBF(km,i) = gamma_AE_hat(km,i)/ &
                      (gamma_AE_hat(km,i)**2 + omega_RBF**2)
!     RBF(km,i) = abs(gamma_AE_hat(km,i))/ &
!                      (gamma_AE_hat(km,i)**2 + omega_RBF**2)

     gamma_hat(1,km,i) = sav_gamma_star_hat(1,km,i)*f_cor(1,i) &
          *(F_mode_drive(km)*rg_n_alpha_tran_rho(i) &
              -sav_rg_n_crit_star(1,km,i)) &
             /sav_rg_n_alpha_rho(i)

     RBF1(1,km,i) =  gamma_hat(1,km,i)/ &
                (gamma_hat(1,km,i)**2 + omega_RBF1**2)
    enddo
   enddo !km  
  endif ! NBI_flag .eq. 1 .or. NBI_flag .eq. 0

  if(NBI_flag .eq. 2) then
   RBF(:,:) = 1.0
   do km=1,km_max
    do i=1,n_rho_grid

       F_frac = F_frac_one_spec_test

      gamma_AE_hat(km,i) = sav_gamma_star_hat(1,km,i)*f_cor(1,i) &
          *(F_mode_drive(km)*rg_n_alpha_tran_rho(i) &
              -F_frac*sav_rg_n_crit_star(1,km,i)) &
             /sav_rg_n_alpha_rho(i)


     gamma_Z_hat(km,i) = sav_gamma_star_hat(1,km,i)*f_cor(1,i) &
             *(F_frac*sav_rg_n_crit_star(1,km,i)) &
             /sav_rg_n_alpha_rho(i)

     gamma_hat(1,km,i) = sav_gamma_star_hat(1,km,i)*f_cor(1,i) &
          *(F_mode_drive(km)*rg_n_alpha_tran_rho(i) &
              -sav_rg_n_crit_star(1,km,i)) &
             /sav_rg_n_alpha_rho(i)

     RBF1(1,km,i) =  gamma_hat(1,km,i)/ &
                (gamma_hat(1,km,i)**2 + omega_RBF1**2)

!add second species

       F_frac = 1.0-F_frac

      gamma_AE_hat(km,i) = gamma_AE_hat(km,i) &
                      +sav_gamma_star_hat(2,km,i)*f_cor(2,i) &
             *(F_mode_drive(km)*rg_n_alpha2_tran_rho(i) &
              -F_frac*sav_rg_n_crit_star(2,km,i)) &
             /sav_rg_n_alpha2_rho(i)

!  Z_hat damping rated same as E_hat damping rate

      gamma_Z_hat(km,i) = gamma_Z_hat(km,i) + sav_gamma_star_hat(2,km,i)*f_cor(2,i) &
             *(F_frac*sav_rg_n_crit_star(2,km,i)) &
             /sav_rg_n_alpha2_rho(i)

      RBF(km,i) = gamma_AE_hat(km,i)/ &
                      (gamma_AE_hat(km,i)**2 + omega_RBF**2)
!     RBF(km,i) = abs(gamma_AE_hat(km,i))/ &
!                      (gamma_AE_hat(km,i)**2 + omega_RBF**2)

      gamma_hat(2,km,i) = sav_gamma_star_hat(2,km,i)*f_cor(2,i) &
          *(F_mode_drive(km)*rg_n_alpha2_tran_rho(i) &
              -sav_rg_n_crit_star(2,km,i)) &
             /sav_rg_n_alpha2_rho(i)

      RBF1(2,km,i) =  gamma_hat(2,km,i)/ &
                (gamma_hat(2,km,i)**2 + omega_RBF1**2)
    enddo
   enddo !km
  endif !NBI_flag.eq. 2

  if(CZ_nl .gt. 0.0) then
! update Z_hat
   do km=1,km_max
    do i=1,n_rho_grid
      Z_hat_prev(km,i) = Z_hat(km,i)

      Z_hat(km,i) = Z_hat(km,i) +dt_update*2.*CZ_nl*E_hat(km,i)*Z_hat(km,i) &
         -dt_update*2.*gamma_Z_hat(km,i)*Z_hat(km,i)
    enddo
   enddo
  endif

! update E_hat
  do km=1,km_max
    do i=1,n_rho_grid
   E_hat(km,i) = E_hat(km,i) + dt_update*2.*gamma_AE_hat(km,i)*E_hat(km,i) &
         -dt_update*2.*C_nl*E_hat(km,i)**2 &
         -dt_update*2.*CZ_nl*E_hat(km,i)*Z_hat_prev(km,i)
!Z_hat_prev  in place of Z_hat insures numerical conservation of CZ_nl transfer
    enddo !i
  enddo !km

  time_called = time_called + dt_update

  if(i_RBF1 .eq. 0) then
   RBF1(1,:,:) = RBF(:,:)
   RBF1(2,:,:) = RBF(:,:)
  endif

   if(i_RBF .eq. 2) then
  do km=1,km_max
   do i=1,n_rho_grid
    if(RBF1(1,km,i).lt. 0) RBF1(1,km,i) = 0.
    if(RBF1(2,km,i).lt. 0) RBF1(2,km,i) = 0.
   enddo
  enddo
   endif

  if(i_RBF .eq. 0)  RBF1(:,:,:) = 1.0  !no RBF factor
  
   
  QLdiffEP(:,:) = 0.0 
  do km=1,km_max
   do i=2,n_rho_grid-1
    QLdiffEP(:,i)=QLdiffEP(:,i) +  &
      E_hat(km,i)*sav_diff_star(:,km,i)*f_cor(:,i)*RBF1(:,km,i)
   enddo !i
  enddo
   QLdiffEP(:,1) = QLdiffEP(:,2)
   QLdiffEP(:,n_rho_grid) = QLdiffEP(:,n_rho_grid-1)
   
   QLdiffEP(:,:) = QLdiffEP(:,:)/float(km_max) + D_bkg
! be careful with /float(km_max)  norm
! really only appropriate for test model and not direct TGLF


  endif !i_sav_QLdiffEP .eq. 0 

  if(i_test_print .eq. 1) then
   print *,  'time_called'
   print *, 'E_hat(1,25)  gamma_AE_hat(1,25)  QLdiffEP(1,25)'
   print *, '------------------------'
  endif

  print *, time_called
  print *, E_hat(1,25), ' ', gamma_AE_hat(1,25) , ' ', QLdiffEP(1,25)


   if(i_test_print .eq. 1)    print *, 'Alpha_QLdiffusivity done'

end subroutine Alpha_QLdiffusivity 
