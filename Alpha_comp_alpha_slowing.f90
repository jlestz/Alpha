!---------------------------------------------------------
! Alpha_comp_alpha_slowing.f90
!
! PURPOSE:
!  compute alpha slowing down profiles
!
!---------------------------------------------------------

subroutine Alpha_comp_alpha_slowing

    use Alpha_use_input      !n_rho_grid,beta_N_ped,n14_ped,I_p,B_t,Rmaj0,delR0oa,rmin,
!                        !kappa_1,kappa_0,delta_1,delta_0,q_1,q_0,rho_q,eps_q,
!                        !Zeff,M_DT,E_alpha,Fpeak_T,Fpeak_n,Fpeak_ei,NBI_flag
!2A
!                        !E_alpha2


    use Alpha_use_output     !beta_ped_percent,T_ped,beta_N_glob,arho,
!                        !q_rho(i),Rmaj_rho(i),rmin_rho(i),
!                        !kappa_rho(i),delta_rho(i),T_i_rho(i),T_e_rho(i),
!                        !n_i_rho(i),n_e_rho(i),
!                        !n_alpha_rho(i),T_alpha_equiv_rho(i),E_c_hat_rho(i)
!2A
!                        !n_alpha2_rho(i),T_alpha2_equiv_rho(i),E_c_hat2_rho(i)
   
    use Alpha_use_other      !n_alpha_ave_rho(i),S0_rho(i),bannana_rho(i)
!2A                          !n_alpha2_ave_rho(i),S02_rho(i),bannana2_rho(i)

    use expro
  !--------------------------------------
  implicit none
     integer :: i

     integer :: ii

     integer :: i_en
     integer :: n_en
     integer :: n_en_max

     integer :: j_newt
     integer :: j_newt_max
 
     real :: Numer0
     real :: Denom0
     real :: G_D
     real :: En
     real :: EnoTe
     real :: Convec_Factor

     real :: Z1
     real :: S0
     real, dimension(n_rho_grid) :: S0_fusion_rho
     real, dimension(n_rho_grid) :: S02_fusion_rho
     real :: tau_ee
     real :: tau_s
     real, dimension(n_rho_grid) :: tau_s_i

     real :: ln_lambda
     
     real :: pi
     real :: a

     real :: numer
     real :: denom
     real :: delta
     real :: rho_t
     real :: rho_p

     real, dimension(n_rho_grid) :: I2
     real, dimension(n_rho_grid) :: I4
     real, dimension(n_rho_grid) :: dI2dE
     real, dimension(n_rho_grid) :: dI4dE
     real, dimension(n_rho_grid) :: E_alpha_rho_0
     real, dimension(n_rho_grid) :: E_alpha2_rho_0
     real, dimension(n_rho_grid) :: T_0
     real, dimension(n_rho_grid) :: dTdE

     real, dimension(n_rho_grid) :: rho_hat_p
!2A
     real, dimension(n_rho_grid) :: z01
     real, dimension(n_rho_grid) :: z02
     real, dimension(n_rho_grid) :: z03
     real, dimension(n_rho_grid) :: z04

     real, dimension(n_rho_grid) :: sigma_v_rho
!        z=0 at pensil beam entrance
!        z=z_NBI at pencil beam tanget to Rmaj_NBI
!        z=2.*z_NBI at pencil beam exit
!        Rmaj_rho(i) = Rmaj0*(1.0+delR0oa*(rmin/Rmaj0)*(1-rho_hat(i)))
!        z_NBI = sqrt(Rmaj0**2-Rmaj_NBI**2)
!
!        Rmaj_NBI  >  Rmaj0 -rmin  so as to miss inner wall minor radius
!             i.e.   pencil beam has only one entrance and exit of plasma
     real, dimension(n_rho_grid) :: V_prime_rho

     integer :: i_NBI
     integer :: i_NBI_max
     integer :: NBI_model
     real :: NBI_profile(7)
     real :: NBI_height(7)
     real :: NBI_center
     real :: NBI_del_height
     
     real :: z_NBI 
     real :: h_NBI
     real :: Rmaj_NBI
     real :: Pow_NBI
     real :: Icur_NBI
     real :: Icur_shine
     real :: lambda_NBI
     real :: z_path   !distance along beam path
     real :: I_path
     real :: lnI_path !ln(I_Path)

     real :: z1m
     real :: z2m
     real :: z3m
     real :: z4m
  !
  character(len=80) :: comment

  logical :: l_exp_T_alpha
  real :: xdum
  character :: dummy

  !--------------------------------------

  pi = 3.141592565

  E_alpha2=1.0  !Mev
  print *, 'E_alpha2 HARDWIRE =',E_alpha2

  l_exp_T_alpha = .false.
  if (E_alpha .lt. 0) then
    l_exp_T_alpha = .true.
    E_alpha = abs(E_alpha)
  endif


!  source_method = 0     ! Calculate source from given classical distribution
!  source_method = 1     ! Use expro source from input.gacode
!  source_method = 2     ! Use externally specified source
  if (trim(source_file) .eq. 'null') then
    source_method = 0
  else if (trim(source_file) .eq. 'null_new') then
    source_method = 0
    source_write_flag = 1
  else if (trim(source_file) .eq. 'input.gacode') then
    source_method = 1
  else if (trim(source_file) .eq. 'fusion') then
    source_method = 3
  else
    source_method = 2
  endif

  if (trim(source_file_2) .eq. 'null') then
    source_method_2 = 0
  else if (trim(source_file_2) .eq. 'null_new') then
    source_method_2 = 0
    source_write_flag_2 = 1
  else if (trim(source_file_2) .eq. 'input.gacode') then
    source_method_2 = 1
  else if (trim(source_file_2) .eq. 'fusion') then
    source_method_2 = 3
  else
    source_method_2 = 2
  endif

  call cub_spline_irregular(DT_Tgrid,DT_sigma_v,301,T_i_rho,sigma_v_rho,n_rho_grid)
  open (unit=1, file='out.Alpha.sigma_v_vs_rho', status='replace')
    write (1,*) '(i_r,rho,Ti,<sigma*v>)'
    do i=1, n_rho_grid
      write (1,*) i, EXPRO_rho(i), T_i_rho(i), sigma_v_rho(i)
    enddo
  close (1)


! rho_hat grid
  rho_hat(1) = 0.
  print *, 'n_rho_grid', n_rho_grid
  do i = 1,n_rho_grid
!   rho_hat(i) = float(i-1)/float(n_rho_grid-1)
   E_alpha_grid(i) = E_alpha*1.E3 - 0.01*i     !not used ?
!2A
   E_alpha2_grid(i) = E_alpha2*1.E3 - 0.01*i     !not used ?

   V_prime_rho(i) = 2.*pi*kappa_rho(i)*rmin_rho(i)*2.*pi*Rmaj_rho(i) !in m**2
!   print *, 'Test of V_prime_rho parts:', i, V_prime_rho(i), kappa_rho(i), rmin_rho(i), Rmaj_rho(i), pi
  enddo

!  print *, 'Midpoint V_prime_rho(25)', V_prime_rho(25)
  
!------------------------------------------------------------------------------------
!------------------------------------------------------------------------------------


 if ((l_read_exp_profile .ne. 1) .and. (NBI_flag .eq. 0 .or. NBI_flag .eq. 2)) then
!  Z1=5./3.
   Z1 = M_beam * (f_i1*Q_i1**2/m_i1 + (1.0-f_i1)*Q_i2**2/m_i2)

! alpha injection energy normed cross over energy E_c_hat = E_c/E_aplha 
   E_c_hat_rho(:) = 0.
  do i = 1,n_rho_grid
   E_c_hat_rho(i) = (T_e_rho(i)/E_alpha)*1.E-3*(4.*1836.)**(1./3.)*(3.*sqrt(pi)*Z1/4.)**(2./3.)

   a=sqrt(E_c_hat_rho(i))

   I2(i) = 1./3.*log((1.+a**3)/a**3)

   I4(i) = 1./2.-a**2*(1./6.*log((1.-a+a**2)/(1.+a)**2)+1./sqrt(3.)*(atan((2.-a)/a/sqrt(3.))+pi/6.)) 
  enddo

! alpha_slowing down density in 10**19 1/m**3

  n_alpha_rho(:) = 0.
  do i = 1,n_rho_grid
    S0=2.5E-6*(n_i_rho(i))**2*(T_i_rho(i))**2 !in [10**19 1/m**3]/sec
    S0_rho(i)=S0

    ln_lambda=17.
!error 9.24.13   tau_ee=1.088E-4*(T_e_rho(i))**1.5/n_e_rho(i)/ln_lambda ! in sec
    tau_ee=1.088E-3*(T_e_rho(i))**1.5/n_e_rho(i)/ln_lambda ! in sec
!error 9.24.13    tau_s = 9.23E4*tau_ee    
    tau_s = 1836.*tau_ee
    
    I2(i) = 1./3.*log((1.+E_c_hat_rho(i)**1.5)/E_c_hat_rho(i)**1.5)
    
    n_alpha_rho(i)=S0*tau_s*I2(i)  !in [10**19 1/m**3]
  enddo

!  print *, 'tau_s error fixed making n_alpha_rho 5.023x smaller'
  endif  !NBI_flag = 0  or 2
!------------------------------------------------------------------------------------

!------------------------------------------------------------------------------------
  if ((l_read_exp_profile .eq. 1) .or. (NBI_flag .eq. 1)) then
! put in alpha to NBI mass & charge corrections

!  Z1=(5./3.) &                ! This is appropriate for
!       *(M_DT/4.)/(M_DT/2.5)  ! DT into DT                ; EMB

!  Z1 = 1.0                    ! Appropriate for D into D or H into H ; EMB

!  Z1 = 0.5                    ! Appropriate for H into D ; EMB

   Z1 = M_beam * (f_i1*Q_i1**2/m_i1 + (1.0-f_i1)*Q_i2**2/m_i2)

   print *, 'Z1 = ', Z1

!   print *, 'Midpoint V_prime_rho(25)', V_prime_rho(25)

! alpha injection energy normed cross over energy E_c_hat = E_c/E_aplha 
   E_c_hat_rho(:) = 0.
   E_alpha_rho(:) = E_alpha
   j_newt_max = 300
   ! Newton iteration to find E_alpha_rho consistent with known n_sd and T_eff
   if ((l_read_exp_profile == 1) .and. (l_exp_T_alpha)) then
     do j_newt = 1, j_newt_max
       E_alpha_rho_0(:) = E_alpha_rho(:)
       print *, 'Newton iteration ', j_newt, ' E_alpha_rho(50)=', E_alpha_rho(50)
       do i = 1,n_rho_grid
!       E_c_hat_rho(i) = (T_e_rho(i)/E_alpha)*1.E-3*(4.*1836.)**(1./3.)*(3.*sqrt(pi)*Z1/4.)**(2./3.) &
!           *(M_DT/4.)**(1./3.)
         E_c_hat_rho(i) = (T_e_rho(i)/E_alpha_rho_0(i))*1.E-3*(M_beam*1836.)**(1./3.)*(3.*sqrt(pi)*Z1/4.)**(2./3.)
  
         a=sqrt(E_c_hat_rho(i))
  
         I2(i) = 1./3.*log((1.+a**3)/a**3)
!         dI2dE(i) = -(1./E_alpha_rho(i))/(1+a**3)
         dI2dE(i) = -0.5*(1/a**2)*(1./E_alpha_rho(i))/(1+a**3)

         I4(i) = 1./2.-a**2*(1./6.*log((1.-a+a**2)/(1.+a)**2)+1./sqrt(3.)*(atan((2.-a)/a/sqrt(3.))+pi/6.))
!         dI4dE(i) = (1.-2.*I4(i))/E_alpha_rho(i) - a**3/(E_alpha_rho(i)*(1+a)*(a**2-a+1))
         dI4dE(i) = (1./(2.*a*E_alpha_rho(i))) * ( -(1.-2.*I4(i))/a + a**2/((1+a)*(a**2-a+1)) )

         T_0(i)=2./3.*I4(i)/I2(i)*E_alpha_rho(i)*10.**3
         dTdE(i) = (2./3.)*E_alpha_rho(i)*10.**3*(dI4dE(i)/I2(i)-I4(i)*dI2dE(i)/I2(i)**2) &
                    + T_0(i)/E_alpha_rho(i)

                                      ! Damping factor included
                                      !    \/
         E_alpha_rho(i) = E_alpha_rho_0(i)+0.3*(T_alpha_equiv_rho(i)-T_0(i))/dTdE(i)

       enddo  ! Radial grid
     enddo ! Newton iteration j_newt

     else   ! No experimental temperature given, use E_alpha
            ! alpha equivalent Maxwellian temperature in keV
       T_alpha_equiv_rho(:) = 0.     
       do i = 1, n_rho_grid
         E_c_hat_rho(i) = (T_e_rho(i)/E_alpha_rho(i))*1.E-3*(M_beam*1836.)**(1./3.)*(3.*sqrt(pi)*Z1/4.)**(2./3.)

         a=sqrt(E_c_hat_rho(i))

         I2(i) = 1./3.*log((1.+a**3)/a**3)
         I4(i) = 1./2.-a**2*(1./6.*log((1.-a+a**2)/(1.+a)**2)+1./sqrt(3.)*(atan((2.-a)/a/sqrt(3.))+pi/6.))
         T_alpha_equiv_rho(i)=2./3.*I4(i)/I2(i)*E_alpha_rho(i)*10.**3
       enddo

     endif ! Experimental profiles or not


! alpha_slowing down density in 10**19 1/m**3  !ACTUALLY the NBI_slowing down density
! invert above to get Source profile S0_rho(i)=S0 in [10**19 1/m**3]/sec 
! from n_alpha_rho(i) in 10**19 1/m**3
   S0_rho(:)=0.
   do i = 1,n_rho_grid
       ln_lambda=17.
    tau_ee=1.088E-3*(T_e_rho(i))**1.5/n_e_rho(i)/ln_lambda ! in sec

!    tau_s = 1836.*tau_ee &            ! For a "DT" beam actually meaning
!       *(M_DT/4.0)/(1.0/2.0)**2       ! beam same as main ion.  EMB

!    tau_s = 1836.*tau_ee*(1.0/4.0)/(1.0/2.0)**2  ! For a hydrogen beam into D.  EMB

    tau_s = 1836.*tau_ee*(M_beam/4.0)/(1.0/2.0)**2
    tau_s_i(i) = tau_s

    I2(i) = 1./3.*log((1.+E_c_hat_rho(i)**1.5)/E_c_hat_rho(i)**1.5)

! Harcoded fusion source assume D and T are first two ion species
! in input.gacode file
!   The approximate fusion source formula is inadequate for performance projections.
!   Need to use exact formula using DT_sigma_v.dat instead.
!    S0_fusion_rho(i)=2.5E-6*(EXPRO_ni(1,i)*EXPRO_ni(2,i))*(EXPRO_ti(1,i)*EXPRO_ti(2,i)) !in [10**19 1/m**3]/sec
    S0_fusion_rho(i) = EXPRO_ni(1,i)*EXPRO_ni(2,i)*sigma_v_rho(i)*1.0E19 / 4.365 ! in [10**19 1/m^3]/s
                                                                     !      ^ kluge to agree with ONETWO/TGYRO
    S0_rho(i) = n_alpha_rho(i)/(tau_s*I2(i))     ! Source for source_method = 0

   enddo
   if (source_method .eq. 1) S0_rho(:) = expro_qpar_beam(:) / 10.**19   ! Overwrite to use input.gacode particle source
   if (source_method .eq. 2) then                                       ! or separte file source
      open(unit=1, file=trim(source_file), status='old')
      read (1,*) dummy
      do i=1, n_rho_grid
        read (1,*) xdum, S0_rho(i)
        S0_rho(i) = S0_rho(i) / 10.**19
      enddo
   endif
   if (source_method .eq. 3) then
     S0_rho(:) = S0_fusion_rho(:)
     n_alpha_rho(:) = S0_rho(:) * tau_s_i(:) * I2(:)
!     print *, 'Confirm source_method=3'
!     print *, 'n_alpha_rho(1)=', n_alpha_rho(1)
!     print *, 'S0,tau_s,I2(1,25)', S0_rho(1), tau_s_i(1), I2(1), S0_rho(25), tau_s_i(25), I2(25)
   endif
! For source_method/=0, given n_alpha_rho is not classical. Overwrite with classical slowing-down density.
   if (source_method .ne. 0) n_alpha_rho(:) = S0_rho(:) * tau_s_i(:) * I2(:)
!   print *, 'n_alpha_rho(1) check 1', n_alpha_rho(1)
  endif ! NBI_flag = 1 or experimental profile

  if ((l_read_exp_profile .eq. 1) .and. (NBI_flag .eq. 2)) then

     Z1 = M_beam2 * (f_i1*Q_i1**2/m_i1 + (1.0-f_i1)*Q_i2**2/m_i2)

! alpha injection energy normed cross over energy E_c_hat = E_c/E_aplha 
     E_c_hat2_rho(:) = 0.
     E_alpha2_rho(:) = E_alpha2
     j_newt_max = 300
   ! Newton iteration to find E_alpha_rho consistent with known n_sd and T_eff
     do j_newt = 1, j_newt_max
       E_alpha2_rho_0(:) = E_alpha2_rho(:)
       print *, 'Newton iteration ', j_newt, ' E_alpha2_rho(50)=', E_alpha2_rho(50)
       do i = 1,n_rho_grid
         E_c_hat2_rho(i) = (T_e_rho(i)/E_alpha2_rho_0(i))*1.E-3*(M_beam2*1836.)**(1./3.)*(3.*sqrt(pi)*Z1/4.)**(2./3.)

         a=sqrt(E_c_hat2_rho(i))
 
         I2(i) = 1./3.*log((1.+a**3)/a**3)
         dI2dE(i) = -(1./E_alpha2_rho(i))/(1+a**3)

         I4(i) = 1./2.-a**2*(1./6.*log((1.-a+a**2)/(1.+a)**2)+1./sqrt(3.)*(atan((2.-a)/a/sqrt(3.))+pi/6.))
         dI4dE(i) = (1.-2.*I4(i))/E_alpha2_rho(i) - a**3/(E_alpha2_rho(i)*(1+a)*(a**2-a+1))

         T_0(i)=2./3.*I4(i)/I2(i)*E_alpha2_rho(i)*10.**3
         dTdE(i) = (2./3.)*E_alpha2_rho(i)*10.**3*(dI4dE(i)/I2(i)-I4(i)*dI2dE(i)/I2(i)**2) &
                    + T_0(i)/E_alpha2_rho(i)

                                        ! Damping factor included
                                        !    \/
         E_alpha2_rho(i) = E_alpha2_rho_0(i)+0.1*(T_alpha2_equiv_rho(i)-T_0(i))/dTdE(i)

       enddo  ! Radial grid
     enddo ! Newton iteration j_newt

! alpha_slowing down density in 10**19 1/m**3  !ACTUALLY the NBI_slowing down density
! invert above to get Source profile S0_rho(i)=S0 in [10**19 1/m**3]/sec 
! from n_alpha_rho(i) in 10**19 1/m**3
   S02_rho(:)=0.
   do i = 1,n_rho_grid
       ln_lambda=17.
    tau_ee=1.088E-3*(T_e_rho(i))**1.5/n_e_rho(i)/ln_lambda ! in sec

!    tau_s = 1836.*tau_ee &            ! For a "DT" beam actually meaning
!       *(M_DT/4.0)/(1.0/2.0)**2       ! beam same as main ion.  EMB

!    tau_s = 1836.*tau_ee*(1.0/4.0)/(1.0/2.0)**2  ! For a hydrogen beam into D.  EMB

    tau_s = 1836.*tau_ee*(M_beam2/4.0)/(1.0/2.0)**2

    I2(i) = 1./3.*log((1.+E_c_hat2_rho(i)**1.5)/E_c_hat2_rho(i)**1.5)

    S02_rho(i) = n_alpha2_rho(i)/(tau_s*I2(i))
   enddo
   if (source_method_2 .eq. 1) S0_rho(:) = expro_qpar_beam(:) / 10.**19   ! Overwrite to use input.gacode particle source
   if (source_method_2 .eq. 2) then                                       ! or separte file source
      open(unit=1, file=trim(source_file_2), status='old')
      read (1,*) dummy
      do i=1, n_rho_grid
        read (1,*) xdum, S02_rho(i)
!        S02_rho(i) = S02_rho(i) / 10.**19
      enddo
   endif
   if (source_method_2 .eq. 3) then
     S02_rho(:) = S02_fusion_rho(:)
     n_alpha2_rho(:) = S02_rho(:) * tau_s_i(:) * I2(:)
!     print *, 'Confirm source_method=3'
!     print *, 'n_alpha_rho(1)=', n_alpha_rho(1)
!     print *, 'S0,tau_s,I2(1,25)', S0_rho(1), tau_s_i(1), I2(1), S0_rho(25), tau_s_i(25), I2(25)
   endif
! For source_method/=0, given n_alpha_rho is not classical. Overwrite with classical slowing-down density.
   if (source_method_2 .ne. 0) n_alpha2_rho(:) = S02_rho(:) * tau_s_i(:) * I2(:)


  endif  ! experimental profiles and NBI_flag = 2


!------------------------------------------------------------------------------------

!------------------------------------------------------------------------------------

!calculate smeared alpha density  n_alpha_ave_rho(i)

   print *, 'Midpoint V_prime_rho(25)', V_prime_rho(25)   
   rho_hat_p(:)=rho_hat(:)
   do i=2,n_rho_grid
    rho_t=1.02E2*(1./2.)*sqrt(4.)*sqrt(2.*T_alpha_equiv_rho(i)*1000.)/(B_t*10000.)/(rmin*100.)
    rho_p=rho_t*rmaj_rho(i)/rmin_rho(i)*q_rho(i)
!   rho_p always greater than rho_t
    delta=rho_p
    if(rho_p.gt.rho_hat(i)) delta=rho_hat(i)  !banana half width less than radius not allowed
    if(delta.lt.rho_t) delta=rho_t   ! there is aways a minimum smearing over rho_t
    bannana_rho(i)=delta
    numer=0.
    denom=0.
    do ii=1,n_rho_grid
     denom=denom+exp(-(rho_hat(i)-rho_hat_p(ii))**2/delta**2)
     numer=numer+exp(-(rho_hat(i)-rho_hat_p(ii))**2/delta**2)*n_alpha_rho(ii)
    enddo
    n_alpha_ave_rho(i)=numer/denom
   enddo
   n_alpha_ave_rho(1)=n_alpha_ave_rho(2)


! compute the increase in energy flux over simple 3/2 T_alpha convection

   n_en = 1000
!   n_en_max = 200
   n_en_max = 1000
   print *, 'n_en'
   print *, 'Midpoint V_prime_rho(25)', V_prime_rho(25) 
   print *, 'Convec_factor'
!   do i=1,n_rho_grid
      i=5
!     i=25
!     i=45
     print *, 'I4/I2=',I4(i)/I2(i)
     print *, 'T_e_rho=',T_e_rho(i)
     print *, 'E_c_hat_rho=',E_c_hat_rho(i)
     Numer0=0.0
     Denom0=0.0
     do i_en = 1,n_en_max
      En = (float(i_en)-0.5)/float(n_en)
      EnoTe = En*E_alpha*10**3/T_e_rho(i)
      G_D=exp(-8.14E-5*EnoTe**4+3.77E-3*EnoTe**3-0.0553*EnoTe**2+0.036*EnoTe+0.45)
!      if(EnoTe .le. 2.7) G_D =1.25  !Angioni-Peters Fig 2 Eq. 32  
!test   passed G_D = 1.0 test
!      G_D=1.0
       if(EnoTe .le. 2.7) then
        G_D=1.0
        if(EnoTe .gt. 0.3) G_D = 0.25*(EnoTe-0.3)/(2.7-0.3)+1.0
       endif
      Numer0=Numer0+1./float(n_en)*En*G_D*En**0.5/(E_c_hat_rho(i)**1.5+En**1.5)
      Denom0=Denom0+1./float(n_en)*G_D*En**0.5/(E_c_hat_rho(i)**1.5+En**1.5)
      print *, i_en,En,EnoTe,G_D,Numer0,Denom0,Numer0/Denom0
     enddo
     Convec_factor = Numer0/Denom0/(I4(i)/I2(i))
     print *, i,Numer0,Denom0,Convec_factor
!   enddo


!------------------------------------------------------------------------------------

!------------------------------------------------------------------------------------
!2A
  if(NBI_flag .eq. 2) then   !alpha2 

! put in alpha to NBI mass & charge corrections
  Z1=(5./3.) &
       *(M_DT/4.)/(M_DT/2.5)

! alpha2 injection energy normed cross over energy E_c_hat = E_c/E_aplha
  E_c_hat2_rho(:) = 0.
  print *, 'Midpoint V_prime_rho(25)', V_prime_rho(25)
  do i = 1,n_rho_grid
   E_c_hat2_rho(i) = (T_e_rho(i)/E_alpha2)*1.E-3*(4.*1836.)**(1./3.)*(3.*sqrt(pi)*Z1/4.)**(2./3.) &
       *(M_DT/4.)**(1./3.)

   a=sqrt(E_c_hat2_rho(i))

   I2(i) = 1./3.*log((1.+a**3)/a**3)

   I4(i) = 1./2.-a**2*(1./6.*log((1.-a+a**2)/(1.+a)**2)+1./sqrt(3.)*(atan((2.-a)/a/sqrt(3.))+pi/6.))
  enddo

! alpha2 equivalent Maxwellian temperature in keV
   T_alpha2_equiv_rho(:) = 0.
   do i = 1,n_rho_grid
    T_alpha2_equiv_rho(i)=2./3.*I4(i)/I2(i)*E_alpha2*10.**3
   enddo


!------------------------------------------------------------------------------------
!------------------------------------------------------------------------------------

  if (l_read_exp_profile .eq. 0) then
!get S02_rho(i) from pencil beam NBI power Pow2 in MW calculation.  Pencil beam tangent to R2 in meters 
!and decays by lambda normed by desnity variation

   S02_rho(:) = 0.0
    !2A
!add pensil beam path
  print *, '------------------------------------------------------------------'
  print *, 'pencil beam path setup'
  print *, '------------------------------------------------------------------'

!  Rmaj_NBI = 0.8*Rmaj0
  Rmaj_NBI =1.0*Rmaj0   ! center tangent
!  Rmaj_NBI =1.2*Rmaj0  !no good
!   Rmaj_NBI = 1.1*Rmaj0 !no good
  Pow_NBI = 30.0 !MW
  Icur_NBI = Pow_NBI/E_alpha2/1.6022E-19  ! 1/sec
  lambda_NBI = 2.5*rmin
  print *, 'Pow_NBI=',Pow_NBI,' Icur_NBI=',Icur_NBI,' lambda_NBI=',lambda_NBI

  print *, 'Midpoint V_prime_rho(25)', V_prime_rho(25)

!  NBI_model = 3  !beam profile and heights
  NBI_model = 2 ! 2*i_NBI_max -1 equal beamlets
!  NBI_model = 1 ! 2*i_NBI_max -1 beamlets with center beamlet having 2x more current
   
  print *, 'NBI_model=',NBI_model
  
  i_NBI_max = 7
!  i_NBI_max = 5
!   i_NBI_max = 3
!   i_NBI_max = 2
  print *, 'i_NBI_max=',i_NBI_max

!set NBI_model = 3

  NBI_center =0.0
!  NBI_center = 0.3
!  NBI_del_height = 0.05
   NBI_del_height = 0.1

  print *, 'NBI_center=',NBI_center,' NBI_del_height=',NBI_del_height

  NBI_profile(1)=0.2
  NBI_profile(2)=0.4
  NBI_profile(3)=0.8
  NBI_profile(4)=1.0
  NBI_profile(5)=0.8
  NBI_profile(6)=0.4
  NBI_profile(7)=0.2

  NBI_profile(:) = NBI_profile(:)/3.8

  NBI_height(1)=NBI_center -3.*NBI_del_height
  NBI_height(2)=NBI_center -2.*NBI_del_height
  NBI_height(3)=NBI_center -1.*NBI_del_height
  NBI_height(4)=NBI_center 
  NBI_height(5)=NBI_center +1.*NBI_del_height
  NBI_height(6)=NBI_center +2.*NBI_del_height
  NBI_height(7)=NBI_center +3.*NBI_del_height

  NBI_height(:) = abs(NBI_height(:))

  Icur_shine = 0.0

  print *, 'Midpoint V_prime_rho(25)', V_prime_rho(25)

  do i_NBI = 1,i_NBI_max  ! i_NBI = 1  is zero beam hight

  print *, '***************************'
  print *, 'i_NBI=',i_NBI
  print *, '***************************'


   if(NBI_model .le. 2)  h_NBI = rmin*rho_hat(i_NBI)
!    the actual height above midplane is more like kappa X rmin*rho_hat(i_NBI)
   if(NBI_model .eq. 3)  h_NBI = NBI_height(i_NBI)
!    NBI_model = 3 is wrong

  z_NBI = sqrt((Rmaj_rho(n_rho_grid)-h_NBI+rmin*rho_hat(n_rho_grid))**2-Rmaj_NBI**2)

  z1m =  sqrt((Rmaj_rho(n_rho_grid)-h_NBI+rmin*rho_hat(n_rho_grid))**2-Rmaj_NBI**2) &
          -sqrt(abs((Rmaj_rho(i_NBI)+rmin*(rho_hat(1)-rho_hat(i_NBI)))**2-Rmaj_NBI**2))

  z2m = z_NBI

  z3m = sqrt((Rmaj_rho(n_rho_grid)-h_NBI+rmin*rho_hat(n_rho_grid))**2-Rmaj_NBI**2) &
          +sqrt(abs((Rmaj_rho(i_NBI)-rmin*(rho_hat(1)-rho_hat(i_NBI)))**2-Rmaj_NBI**2))

  z4m = 2.*z_NBI

  print *, 'i_NBI=',i_NBI
  print *, 'z_NBI=',z_NBI, ' Rmaj_NBI=',Rmaj_NBI, ' h_NBI=',h_NBI
  print *, 'z1m=',z1m,' z2m=',z2m,' z3m=',z3m,' z4m=',z4m

  print *, 'Midpoint V_prime_rho(25)', V_prime_rho(25)

  do i = 1,n_rho_grid

   z01(i) = sqrt((Rmaj_rho(n_rho_grid)-h_NBI+rmin*rho_hat(n_rho_grid))**2-Rmaj_NBI**2) &
          -sqrt(abs((Rmaj_rho(i)-h_NBI+rmin*rho_hat(i))**2-Rmaj_NBI**2))

   z02(i) = sqrt((Rmaj_rho(n_rho_grid)-h_NBI+rmin*rho_hat(n_rho_grid))**2-Rmaj_NBI**2) &
          -sqrt(abs((Rmaj_rho(i)+h_NBI-rmin*rho_hat(i))**2-Rmaj_NBI**2))

   z03(i) = sqrt((Rmaj_rho(n_rho_grid)-h_NBI+rmin*rho_hat(n_rho_grid))**2-Rmaj_NBI**2) &
          +sqrt(abs((Rmaj_rho(i)+h_NBI-rmin*rho_hat(i))**2-Rmaj_NBI**2))

   z04(i) = sqrt((Rmaj_rho(n_rho_grid)-h_NBI+rmin*rho_hat(n_rho_grid))**2-Rmaj_NBI**2) &
          +sqrt(abs((Rmaj_rho(i)-h_NBI+rmin*rho_hat(i))**2-Rmaj_NBI**2))

   print *, 'rho_hat=',rho_hat(i)
   print *, ' ',z01(i),' ',z02(i),' ',z03(i),' ',z04(i)
  enddo

  print *, '------------------------------------------------------------------'

   z_path = 0.

   if(NBI_model .eq. 3) I_path = NBI_profile(i_NBI)*Icur_NBI

   if(NBI_model.eq. 2) then
    I_path = Icur_NBI/float(2*i_NBI_max-1)
    if(i_NBI .gt. 1) I_path = 2.*I_path
   endif

   if(NBI_model .eq. 1)   I_path = Icur_NBI/float(i_NBI_max)
  
   print *, 'Midpoint V_prime_rho(25)', V_prime_rho(25)

   lnI_path = alog(I_path)
!check
   print *, 'check', I_path, '=?', exp(lnI_path)
   print *, 'rmin     z_path      I_path'
   print *, 'start pathe 1'
   print *, rmin*rho_hat(n_rho_grid),' ',z_path, ' ',I_path
!start path 1
   do i=n_rho_grid-1,i_NBI+1,-1
    if(Rmaj_rho(i) .gt.  Rmaj_NBI) then
     lnI_path = lnI_path - abs(z01(i)-z01(i-1))/lambda_NBI*(n_e_rho(i)/n_e_rho(n_rho_grid))
     z_path = z01(i)
     I_path = exp(lnI_path)
     print *, rmin*rho_hat(i),' ',z_path, ' ',I_path
     S02_rho(i) = S02_rho(i) + I_path*abs(z01(i)-z01(i-1))/lambda_NBI* &
              (n_e_rho(i)/n_e_rho(n_rho_grid))/(abs(rho_hat(i)-rho_hat(i-1))*rmin)/V_prime_rho(i)
     print *, 'Test of S02_rho parts: ', i, I_path, z01(i), lambda_NBI, n_e_rho(i), V_prime_rho(i)
    endif
   enddo

   print *, 'start path 2'
!start path 2
   do i=1+i_NBI,n_rho_grid
    if(Rmaj_rho(i) .gt.  Rmaj_NBI) then
     lnI_path = lnI_path - abs(z02(i)-z02(i-1))/lambda_NBI*(n_e_rho(i)/n_e_rho(n_rho_grid))
     z_path = z02(i)
     I_path = exp(lnI_path)
     print *, rmin*rho_hat(i),' ',z_path, ' ',I_path
     S02_rho(i) = S02_rho(i) + I_path*abs(z02(i)-z02(i-1))/lambda_NBI* &
              (n_e_rho(i)/n_e_rho(n_rho_grid))/(abs(rho_hat(i)-rho_hat(i-1))*rmin)/V_prime_rho(i)
    endif
   enddo

   print *, 'start path 3'
!start path 3 
   do i=n_rho_grid-1,1+i_NBI,-1
    if(Rmaj_rho(i) .gt.  Rmaj_NBI) then
     lnI_path = lnI_path - abs(z03(i)-z03(i-1))/lambda_NBI*(n_e_rho(i)/n_e_rho(n_rho_grid))
     z_path = z03(i)
     I_path = exp(lnI_path)
     print *, rmin*rho_hat(i),' ',z_path, ' ',I_path
     S02_rho(i) = S02_rho(i) + I_path*abs(z03(i)-z03(i-1))/lambda_NBI* &
              (n_e_rho(i)/n_e_rho(n_rho_grid))/(abs(rho_hat(i)-rho_hat(i-1))*rmin)/V_prime_rho(i)
    endif
   enddo

   print *, 'start path 4'
!start path4 
   do i=1+i_NBI,n_rho_grid
    if(Rmaj_rho(i) .gt.  Rmaj_NBI) then
     lnI_path = lnI_path - abs(z04(i)-z04(i-1))/lambda_NBI*(n_e_rho(i)/n_e_rho(n_rho_grid))
     z_path = z04(i)
     I_path = exp(lnI_path)
     print *, rmin*rho_hat(i),' ',z_path, ' ',I_path
     S02_rho(i) = S02_rho(i) + I_path*abs(z04(i)-z04(i-1))/lambda_NBI* &
              (n_e_rho(i)/n_e_rho(n_rho_grid))/(abs(rho_hat(i)-rho_hat(i-1))*rmin)/V_prime_rho(i)
    endif
   enddo
   
   Icur_shine = Icur_shine+I_path

 enddo !i_NBI loop
   print *, 'Icur_shine=',Icur_shine,' Icur_NBI=',Icur_NBI

   S02_rho(1) = S02_rho(2)

!3.07.16
   if(S02_rho(n_rho_grid) .eq. 0.0) S02_rho(n_rho_grid) = S02_rho(n_rho_grid-1)

!renorm 1/m**3/sec to 10**19/m**2/sec
   S02_rho(:) = S02_rho(:)/(1.0E19)
!------------------------------------------------------------------------------------
!------------------------------------------------------------------------------------


! alpha2_slowing down density in 10**19 1/m**3  !ACTUALLY the NBI_slowing down density
! from Source profile S02_rho(i)=S0 in [10**19 1/m**3]/sec
   do i = 1,n_rho_grid
       ln_lambda=17.
    tau_ee=1.088E-3*(T_e_rho(i))**1.5/n_e_rho(i)/ln_lambda ! in sec
    tau_s = 1836.*tau_ee &            ! For a "DT" beam actually meaning
       *(M_DT/4.0)/(1.0/2.0)**2       ! beam same as main ion.

    I2(i) = 1./3.*log((1.+E_c_hat2_rho(i)**1.5)/E_c_hat2_rho(i)**1.5)

!    S0_rho(i) = n_alpha_rho(i)/(tau_s*I2(i))
    n_alpha2_rho(i) = (tau_s*I2(i))*S02_rho(i)  ![10**19/m**3]
!    print *, 'Test of n_alpha2 parts:', i, tau_s, I2(i), S02_rho(i), E_c_hat2_rho(i)
   enddo

   endif ! l_read_exp_profiles = 0

  endif !NBI_flag = 2

!  print *, 'n_alpha_rho(1) at end of Alpha_comp_alpha_slowing', n_alpha_rho(1)

end subroutine Alpha_comp_alpha_slowing

subroutine cub_spline_irregular(x,y,n,xi,yi,ni)

  !-------------------------------------------------------------
  implicit none
  !
  integer :: i,ii
  double precision :: x0
  !
  integer, intent(in) :: n
  double precision, intent(in), dimension(n) :: x,y
  !
  integer, intent(in) :: ni
  double precision, dimension(ni) :: xi,yi
  !
  ! LAPACK working variables
  !
  integer :: info
  integer, dimension(n) :: ipiv
  double precision, dimension(n) :: c,z
  double precision, dimension(n-1) :: zl,zu,h
  double precision, dimension(n-2) :: zu2
  double precision, dimension(n-1) :: b,d
  !-------------------------------------------------------------

  !-------------------------------------------------------------
  ! Check to see that interpolated point is inside data interval
  !
  if (xi(ni) > x(n)) then
     print *,'ERROR: (cub_spline) Data above upper bound'
     print *,'xi(ni) > x(n)',xi(ni),x(n)
     print *,'y(:)',y(:)
     stop
  endif
  if (xi(1) < x(1)) then
     print *,'ERROR: (cub_spline) Data below lower bound'
     print *,'xi(1) < x(1)',xi(1),x(1)
     stop
  endif
  !-------------------------------------------------------------

  !-------------------------------------------------------------
  ! Define coefficients of spline matrix 
  !
  do i=1,n-1
     h(i)  = x(i+1)-x(i)
     zl(i) = h(i)
     zu(i) = h(i)
  enddo
  zl(n-1) = 0d0
  zu(1)   = 0d0

  z(1) = 1d0
  c(1) = 0d0
  do i=2,n-1
     z(i) = 2d0*(h(i-1)+h(i))
     c(i) = 3d0*((y(i+1)-y(i))/h(i)-(y(i)-y(i-1))/h(i-1))
  enddo
  z(n) = 1d0
  c(n) = 0d0
  !-------------------------------------------------------------

  !-------------------------------------------------------------
  ! Solve the system using LAPACK
  !
  call DGTTRF(n,zl,z,zu,zu2,ipiv,info)
  call DGTTRS('N',n,1,zl,z,zu,zu2,ipiv,c,n,info)
  !-------------------------------------------------------------

  !-------------------------------------------------------------
  ! Find remaining polynomial coefficients:
  !
  c(n) = 0d0
  do i=1,n-1
     b(i) = (y(i+1)-y(i))/h(i)-h(i)*(2d0*c(i)+c(i+1))/3d0
     d(i) = (c(i+1)-c(i))/(3d0*h(i))
  enddo
  !-------------------------------------------------------------

  !-------------------------------------------------------------
  ! Using known polynomial coefficients, perform interpolation.
  !
  !  S(x) = y(i) + b(i) [x-x(i)] + c(i) [x-x(i)]^2
  !                                        + d(i) [x-x(i)]^3
  !
  do ii = 1, ni
    i = 1
    x0 = xi(ii)
    do while (x(i) <= x0)
      i = i+1
    enddo
!    print *, 'Spline', ii, x0, i, x(i)
    yi(ii) = y(i)+(x0-x(i))*(b(i)+(x0-x(i))*(c(i)+(x0-x(i))*d(i)))
  enddo

  !-------------------------------------------------------------

end subroutine cub_spline_irregular

