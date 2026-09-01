!!---------------------------------------------------------
! Alpha_diffusivity.f90
!
! PURPOSE:
! 
!   compute the the diffusivity profile for two EP specied 
!   OUTPUT:   diffEP(:,:)  [m**2/sec]  
!    from
!   INPUT:    den(:,:) transported density  [10**19 1/m**3] 
!   
!    calld by general time dependent transport 
!     Alpha_time_dep_transport.f90
!
!   notation in Alpha_transport.f90

!---------------------------------------------------------

subroutine Alpha_diffusivity(den,diffEP,i_test_print,i_sav_diffEP,i_cgm)

    use Alpha_use_input

    use Alpha_use_output

    use Alpha_use_other

    use Alpha_use_sav_diffEP


  !--------------------------------------
  implicit none
  
     integer :: i  !i=1,n_rho_grid)   rho_hat(1) = 0.   rho_hat(n_rho_grid) =1.0
     integer :: i_s !species index   i_s=1,i_s_max
     integer :: i_s_max
!     i_s_max=1 NBI_flag=0 or 1,
!     i_s_max=2 NBI_flag=2 

     integer :: i_bkg_model
     integer :: i_AE_model
     
     
     integer :: i_tot_TAE

     integer :: i_threshold
     integer :: i_threshold2

     integer, intent(in) :: i_test_print  !for test on first call
     integer, intent(in) :: i_sav_diffEP  !=1 yes save   =0 already saved
     integer, intent(in) :: i_cgm         !=1 yes CGM on =0 CGM off

     real :: pi
     real :: dr
     real, dimension(n_rho_grid) :: V_prime_rho


     real,  dimension(2,n_rho_grid), intent(in) :: den  !n_alpha_tran_rho, n_alpha2_tran_rho
     real,  dimension(2,n_rho_grid), intent(out) :: diffEP    !EP diffusivity  [m**2/sec]


!general working
     real, dimension(n_rho_grid) :: n_alpha_tran_rho
     real, dimension(n_rho_grid) :: n_alpha2_tran_rho
     real, dimension(n_rho_grid) :: p_alpha_rho
     real, dimension(n_rho_grid) :: p_alpha2_rho
     real, dimension(n_rho_grid) :: p_alpha_tot_rho

     real, dimension(n_rho_grid) :: rg_n_alpha_tran_rho
     real, dimension(n_rho_grid) :: rg_n_alpha2_tran_rho
     real, dimension(n_rho_grid) :: rg_p_alpha_tran_rho
     real, dimension(n_rho_grid) :: rg_p_alpha2_tran_rho

     real, dimension(n_rho_grid) :: rg_p_alpha_tot_tran_rho
     real, dimension(n_rho_grid) :: rg_p_alpha_fac_rho
     real, dimension(n_rho_grid) :: rg_p_alpha2_fac_rho

     real, dimension(n_rho_grid) :: rg_n_alpha_th_rho
     real, dimension(n_rho_grid) :: rg_n_alpha2_th_rho
     real, dimension(n_rho_grid) :: rg_p_alpha_th_rho
     real, dimension(n_rho_grid) :: rg_p_alpha2_th_rho
     real, dimension(n_rho_grid) :: rg_p_alpha_tot_th_rho

     real :: D_bkg
     real :: D_TAE
     real, dimension(n_rho_grid) :: D_alpha
     real, dimension(n_rho_grid) :: D_alpha2

! Angioni model controls and profiles
     integer :: i_bkg_Angioni
     real :: Angi_Pinch_fac
     real :: Angi_exp
     integer :: Angi_negative
     real :: pinch_check

     real :: Q_fus
     real, dimension(n_rho_grid) :: flux_source_rho
     real, dimension(n_rho_grid) :: chi_eff
     real, dimension(n_rho_grid) :: D_bkg_Angioni
     real, dimension(n_rho_grid) :: C_p_alpha
     real, dimension(n_rho_grid) :: D_bkg_Angioni2
     real, dimension(n_rho_grid) :: C_p_alpha2


! some HARDWIRES
   i_bkg_model = 1  !Angioni model
   D_bkg = 0.3 !m**2/sec

   i_threshold = 10
   i_threshold2 = 10
   D_TAE = 1.0

   i_tot_TAE = 1
! end some HARDWIRES

!  set up grid
!  real, dimension(n_rho_grid) :: rho_hat

   dr = rmin/float(n_rho_grid-1)  ![m]

   pi = 3.141592565

  do i = 1,n_rho_grid
   rho_hat(i) = float(i-1)/float(n_rho_grid-1)
  enddo

  do i= 1,n_rho_grid
   V_prime_rho(i) = 2.*pi*kappa_rho(i)*rmin_rho(i)*2.*pi*Rmaj_rho(i) !in m**2
  enddo

!duplicated in Alpha_transport.f90 & Alpha_comp_eq_plasma.90
!caution: normal n_rho_grid=51 input  and nn=51 in Alpha_use_input &
!Alpha_use_output


!   general setup
!   --------------------------------------------------------------

!  setup critical gradients

  if(i_sav_diffEP .eq. 1) then

  if(i_threshold .eq. 10) then
   rg_n_alpha_th_rho(:) = crit_grad_n_alpha_rho(:)   !TGLF input
  endif

! set-up alpha critical pressure gardient
  do i = 2,n_rho_grid-1
   rg_p_alpha_th_rho(i)= T_alpha_equiv_rho(i)*rg_n_alpha_th_rho(i)*0.16022* &
    (1. +(T_alpha_equiv_rho(i+1)-T_alpha_equiv_rho(i-1))/T_alpha_equiv_rho(i)/ &
           (n_alpha_rho(i+1)-n_alpha_rho(i-1))*n_alpha_rho(i))
!   rg_p_alpha_th_rho(i)= T_alpha_equiv_rho(i)*rg_n_alpha_th_rho(i)*0.16022  + &
!        ((T_alpha_equiv_rho(i+1)-T_alpha_equiv_rho(i-1))/(rmin_rho(i+1)-rmin_rho(i)))*n_alpha_rho(i)*0.16022
  enddo
   rg_p_alpha_th_rho(1) = rg_p_alpha_th_rho(2)
   rg_p_alpha_th_rho(n_rho_grid) = rg_p_alpha_th_rho(n_rho_grid-1)

    sav_rg_n_alpha_th_rho(:) = rg_n_alpha_th_rho(:)
    sav_rg_p_alpha_th_rho(:) = rg_p_alpha_th_rho(:)
   endif ! i_sav_diffEP .eq. 1

   if(i_sav_diffEP .eq. 0) then
    if(i_test_print .eq. 1) print *, 'reading rg_n_alpha_th_rho from sav'
    rg_n_alpha_th_rho(:) = sav_rg_n_alpha_th_rho(:)
    rg_p_alpha_th_rho(:) = sav_rg_p_alpha_th_rho(:)
   endif ! i_sav_diffEP .eq. 0

   if(NBI_flag .eq. 2) then

   if(i_sav_diffEP .eq. 1) then

   if(i_threshold2 .eq. 10) then
    rg_n_alpha2_th_rho(:) = crit_grad_n_alpha2_rho(:)   !TGLF input
   endif


! set-up alpha2 critical pressure gradient
  do i = 2,n_rho_grid-1
   rg_p_alpha2_th_rho(i)= T_alpha2_equiv_rho(i)*rg_n_alpha2_th_rho(i)*0.16022* &
 (1. +(T_alpha2_equiv_rho(i+1)-T_alpha2_equiv_rho(i-1))/T_alpha2_equiv_rho(i)/ &
           (n_alpha2_rho(i+1)-n_alpha2_rho(i-1))*n_alpha2_rho(i))
  enddo
   rg_p_alpha2_th_rho(1) = rg_p_alpha2_th_rho(2)
   rg_p_alpha2_th_rho(n_rho_grid) = rg_p_alpha2_th_rho(n_rho_grid-1)

! set-up alpha_tot critical pressure gradient
   do i=1,n_rho_grid
    rg_p_alpha_tot_th_rho(i) = sqrt(rg_p_alpha_th_rho(i)*rg_p_alpha2_th_rho(i))
   enddo
! set-up alpha and alpha2 rg_p factors
   do i=1,n_rho_grid
    rg_p_alpha_fac_rho(i) = sqrt(rg_p_alpha2_th_rho(i)/rg_p_alpha_th_rho(i))
    rg_p_alpha2_fac_rho(i) = sqrt(rg_p_alpha_th_rho(i)/rg_p_alpha2_th_rho(i))
   enddo

   sav_rg_n_alpha2_th_rho(:) = rg_n_alpha2_th_rho(:)
   sav_rg_p_alpha2_th_rho(:) = rg_p_alpha2_th_rho(:)

   sav_rg_p_alpha_tot_th_rho(:) = rg_p_alpha_tot_th_rho(:)
   sav_rg_p_alpha_fac_rho(:) = rg_p_alpha_fac_rho(:)
   sav_rg_p_alpha2_fac_rho(:) = rg_p_alpha2_fac_rho(:)
   
  endif !i_sav_diffEP .eq. 1

  if(i_sav_diffEP .eq. 0) then
   if(i_test_print .eq. 1) print *, 'reading rg_n_alpha2_th_rho from sav'
   rg_n_alpha2_th_rho(:) = sav_rg_n_alpha2_th_rho(:)
   rg_p_alpha2_th_rho(:) = sav_rg_p_alpha2_th_rho(:)

   rg_p_alpha_tot_th_rho(:) = sav_rg_p_alpha_tot_th_rho(:)
   rg_p_alpha_fac_rho(:) = sav_rg_p_alpha_fac_rho(:)
   rg_p_alpha2_fac_rho(:) = sav_rg_p_alpha2_fac_rho(:)
  endif !i_sav_diffEP .eq. 0

  endif !NBI_flag = 2

!end critical gradient setup


!   transported densities
    n_alpha_tran_rho(:) = den(1,:)    ![10**19/m**3]
    n_alpha2_tran_rho(:) = den(2,:)   ![10**19/m**3]

!   radial gradients  for transported density profile
!     rg_n_alpha_tran_rho(i), rg_n_alphai2_tran_rho(i)   -dn/dr [10**19/m**3]/m 
!     rg_p_alpha_tran_rho(i), rg_p_alpha_tran_rho(i)     -dp/dr keV*[10**19/m**3]/m 

!alpha
  rg_n_alpha_tran_rho(1) = 0.
  do i = 2,n_rho_grid-1
    rg_n_alpha_tran_rho(i)= &
   -(n_alpha_tran_rho(i+1)-n_alpha_tran_rho(i-1))/rmin/(rho_hat(i+1)-rho_hat(i-1))
  enddo
   i=n_rho_grid
    rg_n_alpha_tran_rho(i)= &
    -(n_alpha_tran_rho(i)-n_alpha_tran_rho(i-1))/rmin/(rho_hat(i)-rho_hat(i-1))

   do i =1, n_rho_grid
    p_alpha_rho(i) = n_alpha_rho(i)*T_alpha_equiv_rho(i)*0.16022
   enddo

    rg_p_alpha_tran_rho(1) = 0.
   do i = 2,n_rho_grid-1
    rg_p_alpha_tran_rho(i) = &
        -(n_alpha_tran_rho(i+1)*T_alpha_equiv_rho(i+1) &
           -n_alpha_tran_rho(i-1)*T_alpha_equiv_rho(i-1))*0.16022/rmin/(rho_hat(i+1)-rho_hat(i-1))
   enddo
   i=n_rho_grid
   rg_p_alpha_tran_rho(i) = -(n_alpha_tran_rho(i)*T_alpha_equiv_rho(i) &
           -n_alpha_tran_rho(i-1)*T_alpha_equiv_rho(i-1))*0.16022/rmin/(rho_hat(i)-rho_hat(i-1))


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

   do i =1, n_rho_grid
    p_alpha2_rho(i) = n_alpha2_rho(i)*T_alpha2_equiv_rho(i)*0.16022
   enddo

    rg_p_alpha2_tran_rho(1) = 0.
   do i = 2,n_rho_grid-1
    rg_p_alpha2_tran_rho(i) = &
        -(n_alpha2_tran_rho(i+1)*T_alpha2_equiv_rho(i+1) &
           -n_alpha2_tran_rho(i-1)*T_alpha2_equiv_rho(i-1))*0.16022/rmin/(rho_hat(i+1)-rho_hat(i-1))
   enddo
   i=n_rho_grid
   rg_p_alpha2_tran_rho(i) = -(n_alpha2_tran_rho(i)*T_alpha2_equiv_rho(i) &
           -n_alpha2_tran_rho(i-1)*T_alpha2_equiv_rho(i-1))*0.16022/rmin/(rho_hat(i)-rho_hat(i-1))
       
!  total transported EP pressure  gradient 

    rg_p_alpha_tot_tran_rho(1) = 0.
    do i = 2,n_rho_grid-1
     rg_p_alpha_tot_tran_rho(i) = &
        -(n_alpha_tran_rho(i+1)*T_alpha_equiv_rho(i+1) &
         -n_alpha_tran_rho(i-1)*T_alpha_equiv_rho(i-1))*0.16022/rmin/(rho_hat(i+1)-rho_hat(i-1)) &
               *rg_p_alpha_fac_rho(i) &
        -(n_alpha2_tran_rho(i+1)*T_alpha2_equiv_rho(i+1) &
         -n_alpha2_tran_rho(i-1)*T_alpha2_equiv_rho(i-1))*0.16022/rmin/(rho_hat(i+1)-rho_hat(i-1)) &
               *rg_p_alpha2_fac_rho(i)
     enddo
     i=n_rho_grid
     rg_p_alpha_tot_tran_rho(i) = &
        -(n_alpha_tran_rho(i)*T_alpha_equiv_rho(i) &
          -n_alpha_tran_rho(i-1)*T_alpha_equiv_rho(i-1))*0.16022/rmin/(rho_hat(i)-rho_hat(i-1)) &
               *rg_p_alpha_fac_rho(i) &
        -(n_alpha2_tran_rho(i)*T_alpha2_equiv_rho(i) &
          -n_alpha2_tran_rho(i-1)*T_alpha2_equiv_rho(i-1))*0.16022/rmin/(rho_hat(i)-rho_hat(i-1)) &
               *rg_p_alpha2_fac_rho(i)

!       'alpha_tot presure in 10**4 N/m**2'
    do i = 1,n_rho_grid
      p_alpha_tot_rho(i) = p_alpha_rho(i)+p_alpha2_rho(i)
    enddo

    endif  !NBI_flag = 2

! default
!    diffEP(:,:) = 1.0
!    return
  
!start models
   diffEP(:,:) = 0.

  D_alpha(:) = D_bkg
  if(i_bkg_model .eq. 1) then
!--------------Angioni model for bkg microturbulent ITG/TEM EP transport -----
!-----------------------------------------------------------------------------
   i_bkg_Angioni = 1
  

   Q_fus = 10.  !default
!  Q_fus = 20.  !for the 2x baseline case
    if(NBI_flag .eq. 1)      Q_fus = 1000000. 
      !for  DIIID NBI test with effective  E_alpha
  Angi_Pinch_fac = 1.0 ! default
!9.18.14  error fix
  Angi_exp = -1.0   !correct    +1.0  !incorrect used in NF  BW2014 ITER paper
  Angi_negative = 1

   chi_eff(:) = 0.

   if(i_sav_diffEP .eq. 1) then

!get chi_eff thermal plasma  normalization for Angioni model 

  !source flux without SDsink
  ! primary EP energy source S)_rho  NBI_flag = 0 or 2  fusion alpha 
  !                                  NBI_flag = 1   NBI like DIIID
    flux_source_rho(1)=0.
     do i=2,n_rho_grid
    flux_source_rho(i) = V_prime_rho(i-1)/V_prime_rho(i)*flux_source_rho(i-1) &
       + 0.5*dr/V_prime_rho(i)* &
         (V_prime_rho(i)*S0_rho(i)*(1.0 -  &
           0.0*n_alpha_tran_rho(i)/n_alpha_rho(i)) + &
            V_prime_rho(i-1)*S0_rho(i-1)*(1.0 - &
             0.0*n_alpha_tran_rho(i-1)/n_alpha_rho(i-1)))
     enddo


  do i = 2,n_rho_grid-1
    chi_eff(i) =  (1.+5./Q_fus)*flux_source_rho(i)*(E_alpha*1000.)/&
     (0.5*n_i_rho(i)*(-T_i_rho(i+1)+T_i_rho(i-1))/(2.*dr) + &
      0.5*n_e_rho(i)*(-T_e_rho(i+1)+T_e_rho(i-1))/(2.*dr))
  enddo
   i=n_rho_grid
    chi_eff(i) = (1.+5./Q_fus)*flux_source_rho(i)*(E_alpha*1000.)/ &
         (0.5*n_i_rho(i)*(-T_i_rho(i)+T_i_rho(i-1))/dr + &
          0.5*n_e_rho(i)*(-T_e_rho(i)+T_e_rho(i-1))/dr)

    chi_eff(1) = chi_eff(2)

   do i = 1,n_rho_grid
        D_bkg_Angioni(i)=chi_eff(i)*&
               (0.02+4.5*(T_e_rho(i)/(E_alpha*1000.))&
               +8.0*(T_e_rho(i)/(E_alpha*1000.))**2&
               +350.*(T_e_rho(i)/(E_alpha*1000.))**3)
  enddo
!DEBUG
!  if(i_test_print .eq. 1) then
!   print *, 'D_bkg_Angioni(i)'
!   do i=1,n_rho_grid
!    print *, D_bkg_Angioni(i)
!   enddo
!  endif
!DEBUG
  do i = 2,n_rho_grid-1
   C_p_alpha(i) = &
      (3./2.)*Rmaj_rho(i)*(-T_e_rho(i+1)+T_e_rho(i-1))/(2.*dr)/T_e_rho(i)*&
                (1./(1.+1./E_c_hat_rho(i)**(Angi_exp*1.5))/log(1.+1./E_c_hat_rho(i)**1.5)-1.)

  enddo
   i=n_rho_grid
   C_p_alpha(i) = (3./2.)*Rmaj_rho(i)*(-T_e_rho(i)+T_e_rho(i-1))/dr/T_e_rho(i)*&
                (1./(1.+1./E_c_hat_rho(i)**(Angi_exp*1.5))/log(1.+1./E_c_hat_rho(i)**1.5)-1.)
   C_p_alpha(1) = 0.

!DEBUG
!  if(i_test_print .eq. 1) then
!   print *, 'C_p_alpha(i)'
!   do i=1,n_rho_grid
!    print *, C_p_alpha(i)
!   enddo
!  endif
!DEBUG

!for speed up
!need only save D_bkg_Angioni(i) and C_p_alpha(i) for later calls 
    sav_D_bkg_Angioni(:) = D_bkg_Angioni(:)
    sav_C_p_alpha(:) = C_p_alpha(:)
    endif !i_sav_diffEP .eq. 1

    if(i_sav_diffEP .eq. 0) then
     if(i_test_print .eq. 1) print *, 'reading D_bkg_Angioni from sav'
     D_bkg_Angioni(:) = sav_D_bkg_Angioni(:)
     C_p_alpha(:) = sav_C_p_alpha(:)
    endif ! i_sav_diffEP eq. 0

    if(i_test_print .eq. 1) then

      print *, '---------------------------------------'
      print *, 'chi_eff(i)'
      do i = 1,n_rho_grid
       print *, chi_eff(i)
      enddo
      print *, '---------------------------------------'

      print *, '---------------------------------------'
      print *, 'pinched D_bkg_Angioni(i) with SD alpha'

   print *,  D_bkg_Angioni(1)
  do i = 2,n_rho_grid-1
   print *,     D_bkg_Angioni(i)*&
   (1.+Angi_Pinch_fac*C_p_alpha(i)*(2.*dr)*n_alpha_rho(i)/(-n_alpha_rho(i+1)+n_alpha_rho(i-1))/Rmaj_rho(i))
  enddo
  i=n_rho_grid
   print *, D_bkg_Angioni(i)*&
   (1.+Angi_Pinch_fac*C_p_alpha(i)*dr*n_alpha_rho(i)/(-n_alpha_rho(i)+n_alpha_rho(i-1))/Rmaj_rho(i))

     print *, '---------------------------------------'

    endif !i_test_print =1

   
   
!----------------------------------------
! set up finished for alpha
!---------------------------------------

!  for any  NBI_flag         !do alpha 
!----------------------------------------


!  start Angioni model
   do i=1,n_rho_grid
      if(i .ne. 1 .and. i .ne. n_rho_grid) then
      D_alpha(i) = D_bkg_Angioni(i)*&
       (1.+Angi_Pinch_fac*C_p_alpha(i)/Rmaj_rho(i)*n_alpha_tran_rho(i)/rg_n_alpha_tran_rho(i))
       if(Angi_negative .eq. 1) D_alpha(i) = abs(D_alpha(i))
      endif
   enddo
!  end Angioni model
  endif !i_bkg_model .eq. 1

   if(i_cgm .eq. 1) then
!  start stiff critical gradient  transport for alpha
    if(i_tot_TAE .eq. 0) then
    do i=1,n_rho_grid
     if (rg_n_alpha_tran_rho(i) .gt. rg_n_alpha_th_rho(i)) D_alpha(i) =  &
            D_alpha(i) + &
               D_TAE*(rg_n_alpha_tran_rho(i)-rg_n_alpha_th_rho(i))*rmin/n_alpha_rho(i)
    enddo
    endif
    if(i_tot_TAE .eq. 1) then
    do i=1,n_rho_grid
      if (rg_p_alpha_tot_tran_rho(i) .gt. rg_p_alpha_tot_th_rho(i)) D_alpha(i) = &
             D_alpha(i) + &
                D_TAE*(rg_p_alpha_tot_tran_rho(i)-rg_p_alpha_tot_th_rho(i))*rmin/p_alpha_tot_rho(i)
    enddo
    endif
    if(i_tot_TAE .eq. -1) then  !11.22.16
    do i=1,n_rho_grid
     if (rg_p_alpha_tran_rho(i) .gt. rg_p_alpha_th_rho(i)) D_alpha(i) = &
             D_alpha(i) + &
                D_TAE*(rg_p_alpha_tran_rho(i)-rg_p_alpha_th_rho(i))*rmin/p_alpha_rho(i)
    enddo
    endif
!   end stiff critical gradient transport for alpha
   endif !i_cgm .eq. 1  alpha

   if(i_bkg_Angioni .eq. 1) then
    D_alpha(1) = D_alpha(2)
    D_alpha(n_rho_grid) =D_alpha(n_rho_grid-1)
   endif
!  end transort for alpha


    if(NBI_flag .eq. 2) then  !do alpha2
!----------------------------------------
    D_alpha2(:) = D_bkg
  if(i_bkg_model .eq. 1) then

   if(i_sav_diffEP .eq. 1) then

  do i = 1,n_rho_grid
        D_bkg_Angioni2(i)=chi_eff(i)*&
               (0.02+4.5*(T_e_rho(i)/(E_alpha2*1000.))&
               +8.0*(T_e_rho(i)/(E_alpha2*1000.))**2&
               +350.*(T_e_rho(i)/(E_alpha2*1000.))**3)
  enddo
!DEBUG
!  if(i_test_print .eq. 1) then
!   print *, 'D_bkg_Angioni2(i)'
!   do i=1,n_rho_grid
!    print *, D_bkg_Angioni2(i)
!   enddo
!  endif
!DEBUG

  do i = 2,n_rho_grid-1
   C_p_alpha2(i) = &
     (3./2.)*Rmaj_rho(i)*(-T_e_rho(i+1)+T_e_rho(i-1))/(2.*dr)/T_e_rho(i)*&
                (1./(1.+1./E_c_hat2_rho(i)**(Angi_exp*1.5))/log(1.+1./E_c_hat2_rho(i)**1.5)-1.)

  enddo
   i=n_rho_grid
   C_p_alpha2(i) =  &
     (3./2.)*Rmaj_rho(i)*(-T_e_rho(i)+T_e_rho(i-1))/dr/T_e_rho(i)*&
                (1./(1.+1./E_c_hat2_rho(i)**(Angi_exp*1.5))/log(1.+1./E_c_hat2_rho(i)**1.5)-1.)
   C_p_alpha2(1) = 0.

!DEBUG
!  if(i_test_print .eq. 1) then
!   print *, 'C_p_alpha2(i)'
!   do i=1,n_rho_grid
!    print *, C_p_alpha2(i)
!   enddo
!  endif
!DEBUG

!for speed up
!need only save D_bkg_Angioni2(i) and C_p_alpha2(i) for later calls
    sav_D_bkg_Angioni2(:) = D_bkg_Angioni2(:)
    sav_C_p_alpha2(:) = C_p_alpha2(:)
   endif !i_sav_diffEP .eq. 1

   if(i_sav_diffEP .eq. 0) then
     if(i_test_print .eq. 1) print *, 'reading D_bkg_Angioni2 from sav'
    D_bkg_Angioni2(:) = sav_D_bkg_Angioni2(:)
    C_p_alpha2(:) = sav_C_p_alpha2(:)
   endif



   if(i_test_print .eq. 1) then
    print *, '---------------------------------------'
      print *, 'pinched D_bkg_Angioni2(i) with SD alpha2'
    print *,  D_bkg_Angioni2(1)
  do i = 2,n_rho_grid-1
     print *,  D_bkg_Angioni2(i)*&
           (1.+Angi_Pinch_fac*C_p_alpha2(i)*(2.*dr)*n_alpha2_rho(i)/(-n_alpha2_rho(i+1)+n_alpha2_rho(i-1))/Rmaj_rho(i))
  enddo
    print *, '---------------------------------------'
   endif

!----------------------------------------
! set up finished for alpha2
!---------------------------------------


!  start Angioni model
    do i=1,n_rho_grid
      if(i .ne. 1 .and. i .ne. n_rho_grid) then
!      D_alpha2(i) = D_bkg_Angioni2(i)*&
!       (1.+Angi_Pinch_fac*C_p_alpha2(i)/Rmaj_rho(i)*n_alpha2_tran_rho(i)/rg_n_alpha2_tran_rho(i))
!!!!!02.03.16       if(Angi_negative .eq. 1) D_alpha2(i) = abs(D_alpha2(i))
!02.04.16AM
      pinch_check = &
       Angi_Pinch_fac*C_p_alpha2(i)/Rmaj_rho(i)*n_alpha2_tran_rho(i)/rg_n_alpha2_tran_rho(i)
      if(pinch_check .lt. -0.9)  pinch_check = -0.9
      D_alpha2(i) = D_bkg_Angioni2(i)*(1.+pinch_check)
    endif !i ne 1 & i ne n_rho_grid
    enddo 
!  end Angioni model
   endif !i_bkg_model .eq. 1

  if(i_cgm .eq. 1) then
! start stiff critical gradient transport for alpha2
    if(i_tot_TAE .eq. 0) then
    do i=1,n_rho_grid
      if (rg_n_alpha2_tran_rho(i) .gt. rg_n_alpha2_th_rho(i)) D_alpha2(i) = &
         D_alpha2(i) + &
             D_TAE*(rg_n_alpha2_tran_rho(i)-rg_n_alpha2_th_rho(i))*rmin/n_alpha2_rho(i)
    enddo
    endif
    if(i_tot_TAE .eq. 1) then
    do i=1,n_rho_grid
      if (rg_p_alpha_tot_tran_rho(i) .gt. rg_p_alpha_tot_th_rho(i)) D_alpha2(i) = &
          D_alpha2(i) + &
             D_TAE*(rg_p_alpha_tot_tran_rho(i)-rg_p_alpha_tot_th_rho(i))*rmin/p_alpha_tot_rho(i)
    enddo
    endif
    if(i_tot_TAE .eq. -1) then  !11.22.16
    do i=1,n_rho_grid
     if (rg_p_alpha2_tran_rho(i) .gt. rg_p_alpha2_th_rho(i)) D_alpha2(i) = &
          D_alpha2(i) + &
              D_TAE*(rg_p_alpha2_tran_rho(i)-rg_p_alpha2_th_rho(i))*rmin/p_alpha2_rho(i)
    enddo
    endif
!  end stiff critical gradient transport for alpha2 
   endif !i_cgm .eq. 1 for alpha2

   if(i_bkg_Angioni .eq. 1) then
    D_alpha2(1) = D_alpha2(2)
    D_alpha2(n_rho_grid) =D_alpha2(n_rho_grid-1)
   endif

! end transport for alpha2

    endif !NBI_flag = 2

     diffEP(1,:) = diffEP(1,:) + D_alpha(:)
     if(NBI_flag .eq. 2) diffEP(2,:) = diffEP(2,:) + D_alpha2(:)




   if(i_test_print .eq. 1)    print *, 'Alpha_diffusivity done'





end subroutine Alpha_diffusivity 
