!!---------------------------------------------------------
! Alpha_time_dep_transport.f90
!
! PURPOSE:
! 
!   compute time trace of EP denity profiles
!    with a general time dependent transport equation
!
!   notation in Alpha_transport.f90
!
!   steady state results should agree with Alpha_transport.f90  for same
!   diffusivit model
!
!   profile grid:  rho_hat (r/a) distinction from rho/rho_a  not important 
!    V_prime_rho(i) = 2.*pi*kappa_rho(i)*rmin_rho(i)*2.*pi*Rmaj_rho(i) !in m**2
!
!   main inputs:
!    source rates:  S0_rho, S02_rho                              [10**19 (1/m**3)/sec]
!    slowing down density profiles:  n_alpha_rho, n_alpha2_rho   [10**19 (1/m**3)]
!
!   main external:  
!    total effective density diffusivies:  D_alpha, D_alpha2    [m**2/sec]
!
!   main outputs at each output time:
!    transported EP density:  n_alpha_tran_rho, n_alpha2_tran_rho [10**19 (1/m**3)]
!    transport EP flows:   V_prim_rho*flux_rho, V_prim_rho*flux2_rho [1/sec]
!
!   NBI_flag  =0 (alpha=fusion alpha)
!   NBI_flag  =1 (alpha=NBI)
!   NBI_flag  =2 (alpha=fusion alpha,  alpha2=NBI)
!---------------------------------------------------------

subroutine Alpha_time_dep_transport

    use Alpha_use_input

    use Alpha_use_output

    use Alpha_use_other


  !--------------------------------------
  implicit none


     integer :: i  !i=1,n_rho_grid)   rho_hat(1) = 0.   rho_hat(n_rho_grid) =1.0
     integer :: i_s !species index   i_s=1,i_s_max
     integer :: i_s_max
!     i_s_max=1 NBI_flag=0 or 1,
!     i_s_max=2 NBI_flag=2 

     real :: delta1  !r=a  BC fraction of EP slowing down density 

!    time step indices
     integer :: it       !time step index
     integer :: it_max   !max time steps
     integer :: it_write !write out every it_write it steps
     integer :: iRK
     real :: time_run
     real :: diffEP_test

!    time plots
     integer :: it_plot   ! e.g. it_plot=10 means plot it=1,11,21, ...<=it_max
     integer :: ir_plot   ! e.g. ir_plot=2  means plot  i=1,3,5,....<=51


     real :: pi
     real :: dr

! simplifies working variables
     real, dimension(n_rho_grid) :: rhat    !rhat_rho  
     real, dimension(n_rho_grid) :: V_p     !V_prime_rho   [m**2]
     real, dimension(n_rho_grid) :: net_of_sink

     real, dimension(2,n_rho_grid) :: den       !n_alpha_tran_rho, n_alpha2_tran_rho 
     real, dimension(2,n_rho_grid) :: den_begin !begin each time step
     real, dimension(2,n_rho_grid) :: den_sd    !n_alpha_rho,n_alpha2_rho
     real, dimension(2,n_rho_grid) :: flux      !flux_rho, flux2_rho
     real, dimension(2,n_rho_grid) :: flux_half !flux at cell interfaces i+1/2, i=1,...,n_rho_grid-1
     real, dimension(2,n_rho_grid) :: S0        !S0_rho, S02_rho
     real, dimension(2,n_rho_grid) :: RHS       ! d den /dt = RHS
     real, dimension(2,n_rho_grid) :: diffEP    !EP diffusivity  [m**2/sec]
     real, dimension(2,n_rho_grid) :: QLdiffEP  !EP diffusivity  [m**2/sec]

! useful gradient profiles
     real, dimension(n_rho_grid) :: n_alpha_tran_rho

     real, dimension(n_rho_grid) :: rg_n_alpha_th_rho
     real, dimension(n_rho_grid) :: rg_n_alpha_tran_rho
     real, dimension(n_rho_grid) :: rg_p_alpha_th_rho
     real, dimension(n_rho_grid) :: rg_p_alpha_tran_rho

     real, dimension(n_rho_grid) :: n_alpha2_tran_rho

     real, dimension(n_rho_grid) :: rg_n_alpha2_th_rho
     real, dimension(n_rho_grid) :: rg_n_alpha2_tran_rho
     real, dimension(n_rho_grid) :: rg_p_alpha2_th_rho
     real, dimension(n_rho_grid) :: rg_p_alpha2_tran_rho
! time advance variable
     real :: ts
     real :: ts2
     real :: dt
     real :: dt_frac    !HARDWIRE
     real :: start_frac !HARDWIRE
     real :: dtRK

     real :: SDsink !slowing down sink normally 1.0
     real :: DS     !drive strength normally 1.0

! number EP's
     i_s_max =1
     if(NBI_flag .eq. 2) i_s_max=2

!HARDWIRE controls
     dt_frac = .0005                     !HARDWIRE
     start_frac = .1  !start from below  !HARDWIRE  !OK with 30000/300 .0000001
!     pcase28 compare pcase39

!    start_frac = 2.0 !start from above  !HARDWIRE
!     start_frac = 1.0
!     start_frac = 1.5 !OK with 50000/500  .00000001

!     start_frac = 1.0 !pcase39
!     it_max=10000
!     it_write=100                        !HARDWIRE
!     it_max =100
!     it_write =1
     
     it_max=30000
     it_write=300

     it_max=60000
     it_write=600

    it_max=600000
    it_write=6000

    it_max=1200000
    it_write=12000
!1.15.18
    it_max=3600000
    it_write=12000
!    it_max=3600000
!    it_write=4000

!     it_max=50000
!     it_write=500
!     diffEP_test = 1.                   !HARDWIRE
     diffEP_test =  0.5                  !HARDWIRE
     delta1 = 0.0  !r=a BC fraction of sd!HARDWIRE

     SDsink = 1.0 !normal
!     SDsink = 0. !test
!      SDsink = 0.5
!      SDsink = 0.2

     DS =1.0 !normal
!      DS = 0.5
!      DS = 0.3
!       DS = 0.7
!        DS = 0.6
!        DS = 0.55

!keep  n_sd fixed while increasing tau_sd   DS/SDsink=1.  
! (-d n_sd/dr)/(-d n_sd/dr)_crit  is fixed
! tau_sd --> tau_sd/SDsink
!!!    SDsink=0.1
!!!    DS = 0.1
       

     ir_plot =1
     it_plot = 10

     print *, 'starting Alpha_time_dep_transport'
     print *, 'HARDWIRES'
     print *, 'dt_frac=',dt_frac 
     print *, 'start_frac=',start_frac
     print *, 'it_max=',it_max
     print *, 'it_write=',it_write
     print *, 'diffEP_test=',diffEP_test
     print *, 'delta1=',delta1
     print *, 'SDsink=',SDsink
     print *, 'DS=',DS
   
     print *, 'ir_plot=',ir_plot
     print *, 'it_plot=',it_plot

!  set up grid
!  real, dimension(n_rho_grid) :: rho_hat

  rho_hat(1) = 0.
  do i = 1,n_rho_grid
   rho_hat(i) = float(i-1)/float(n_rho_grid-1)
  enddo
!duplicated in Alpha_transport.f90 & Alpha_comp_eq_plasma.90
!caution: normal n_rho_grid=51 input  and nn=51 in Alpha_use_input & Alpha_use_output

   dr = rmin/float(n_rho_grid-1)  ![m]

   rhat(:) = rho_hat(:)

!set V_prime_rho(i)   "V_p"
!since  V' enters as 1/V' d[V' D dn/dr]/dr  size of V' or accuracy of V' not vey
!important

  pi = 3.141592565

  do i= 1,n_rho_grid
   V_p(i) = 2.*pi*kappa_rho(i)*rmin_rho(i)*2.*pi*Rmaj_rho(i) ![m**2]
  enddo

! get sources 
    S0(:,:) = 0.
    S0(1,:) = S0_rho(:)
    if(NBI_flag .eq. 2) S0(2,:)=S02_rho(:)

! get slowing down densities
    den_sd(:,:) = 0.1
    den_sd(1,:) = n_alpha_rho(:)
    if(NBI_flag .eq. 2) den_sd(2,:)=n_alpha2_rho(:)

! find typical slowing down time scale
    ts=den_sd(1,25)/S0(1,25)  ![sec]
    ts2=0.0
    if(NBI_flag .eq. 2) ts2=den_sd(2,25)/S0(2,25)
! time step
!    dt = dt_frac**ts !crazy error
     dt = dt_frac*ts
   print *, 'ts=',ts,' dt=',dt
   print *, 'ts2=',ts2
   print *, 'dr=',dr
! set inital density
   den(:,:) = start_frac*den_sd(:,:)

   time_run = 0.0
   print *, 'time_run=',time_run
   print *, 'den_sd(:,1)=',den_sd(:,1)
   print *, 'den(:,1)=',den(:,1)
   print *, '---------------------------------------------'

      call Alpha_diffusivity(den,diffEP,1,1,1)  !yes print, yes sav,  yes CGM
!!!      call Alpha_diffusivity(den,diffEP,1,1,0)  !yes print, yes sav,  no CGM
!      call Alpha_QLdiffusivity(den,QLdiffEP,dt,1,1,0)  !no D_bgk
      call Alpha_QLdiffusivity(den,QLdiffEP,dt,1,1,1)   ! yes D_bgk

!   return   !test before time steps

! start time step loop

! open and start time plot files
  open(unit=6,file='Alpha_time_run.plot',status='replace')
  write(6,*) it_plot  !it_plot=10 means plot it=1,11,21, ...<=it_max
  write(6,*) it_max       !it=1,2,3....it_max
  write(6,*) ir_plot  !ir_plot=2  means plot  i=1,3,5,....<=51
  write(6,*) n_rho_grid   !i=1,2,.....51 normally

  open(unit=5,file='Alpha_den1.plot',status='replace')
  write(5,*) it_plot  !it_plot=10 means plot it=1,11,21, ...<=it_max
  write(5,*) it_max       !it=1,2,3....it_max
  write(5,*) ir_plot  !ir_plot=2  means plot  i=1,3,5,....<=51
  write(5,*) n_rho_grid   !i=1,2,.....51 normally
  
 if(NBI_flag .eq. 2) then
  open(unit=4,file='Alpha_den2.plot',status='replace')
  write(4,*) it_plot  !it_plot=10 means plot it=1,11,21, ...<=it_max
  write(4,*) it_max       !it=1,2,3....it_max
  write(4,*) ir_plot  !ir_plot=2  means plot  i=1,3,5,....<=51
  write(4,*) n_rho_grid   !i=1,2,.....51 normally
 endif

    dtRK = 0.5*dt
    do it=1,it_max
     den_begin(:,:) = den(:,:)
     do iRK =1,2   !2nd RK
      if(iRK .eq. 2) then
        dtRK = dt

!      get EP diffusivity model:  diffEP
!       diffEP(:,:) = diffEP_test  !test diffEP

!         call Alpha_diffusivity(den,diffEP,0,0,1) ! call CGM
         call Alpha_diffusivity(den,diffEP,0,0,0) !do not call CGM 

           diffEP(:,:) =0.0   !no Angioni ITG/TEM  or CGM 

!          QLdiffEP(:,:) = 0.0

         call Alpha_QLdiffusivity(den,QLdiffEP,dtRK,0,0,1) ! call D_bgk
!         call Alpha_QLdiffusivity(den,QLdiffEP,dtRK,0,0,0) ! do not call D_bgk

         diffEP(:,:) = diffEP(:,:) + QLdiffEP(:,:)  !add QL
      endif

!      get EP transport flux at cell interfaces (i+1/2), i=1,...,n_rho_grid-1
!      flux_half(i_s,i) = flux between grid points i and i+1
!      this compact, flux-conservative form replaces the old (i+1)-(i-1)
!      "wide" stencil, which decoupled odd/even grid points and produced
!      grid-scale checkerboard oscillations
       do i_s=1,i_s_max
        do i=1,n_rho_grid-1
         flux_half(i_s,i) = -0.5*(diffEP(i_s,i)+diffEP(i_s,i+1)) &
                              *(den(i_s,i+1)-den(i_s,i))/dr
        enddo
        flux(i_s,1) = 0.
        i=n_rho_grid
        flux(i_s,i) = -diffEP(i_s,i-1)*(den(i_s,i)-den(i_s,i-1))/dr 
       enddo !i_s
      

!      get RHS
      RHS(:,:) = 0.0
      do i_s=1,i_s_max

!      get net source RHS
       do i=1,n_rho_grid
         RHS(i_s,i) = RHS(i_s,i) + &
                DS*S0(i_s,i)*(1.0 - SDsink*den(i_s,i)/(DS*den_sd(i_s,i)))
       enddo ! i
 
!     get transport loss RHS  (flux-conservative: compact 3-point stencil)
       
    !   i=1
    !   RHS(i_s,i) = RHS(i_s,i)  &
    !      -flux(i_s,i+1)/dr
       do i=2,n_rho_grid-1     
         RHS(i_s,i) = RHS(i_s,i)  &
           -( 0.5*(V_p(i)+V_p(i+1))*flux_half(i_s,i)   &
             -0.5*(V_p(i-1)+V_p(i))*flux_half(i_s,i-1) ) &
             /(V_p(i)*dr)
       enddo ! i
       i=n_rho_grid
       RHS(i_s,i) = RHS(i_s,i)  &
           -1./V_p(i)/dr*(V_p(i)*flux(i_s,i)-V_p(i-1)*flux_half(i_s,i-1))

      enddo ! i_s


! update den
     den(:,:) = den_begin(:,:) + dtRK*RHS(:,:)
     den(:,1) = den(:,2)  ! r=0 BC
     den(:,n_rho_grid) = delta1*den_sd(:,n_rho_grid)  !r=a BC

     enddo ! iRK
     time_run = time_run + dt

   if(modulo(it,it_write) .eq. 0) then
     print *,it, den(:,1)
   endif

! write time plot files


   if(modulo(it,it_plot) .eq. 0) then
     write(6,*) time_run

     do i=1,n_rho_grid,ir_plot
      write(5,*) den(1,i) 
     enddo
    if(NBI_flag .eq. 2) then
     do i=1,n_rho_grid,ir_plot
      write(4,*) den(2,i)
     enddo
    endif !NBI_flag .eq. 2
   endif !modulo
! end time plot file writes

    enddo ! it

   if(NBI_flag .eq. 2) then
    close(4)  !den2
   endif
    close(5)  !den1
    close(6)  !time_run
    
    print *, 'den final','time_run=',time_run
    do i=1,n_rho_grid
     print *, den(:,i)
    enddo

    call Alpha_QLdiffusivity(den,QLdiffEP,dt,1,2,0)

!  print DriveStrength  measure 
    print *, 'DriveStrength profiles based on density'
    print *, '--------------------------'
    print *, '--------------------------'

    print *, 'DriveStrength profile for alpha'
    print *, '--------------------------'
   do i=2,n_rho_grid-1
    print *, i, ' ', & 
    -(den(1,i+1)-den(1,i-1))/rmin/(rho_hat(i+1)-rho_hat(i-1))/ &
      crit_grad_n_alpha_rho(i)
   enddo
    print *, '--------------------------'

   if(NBI_flag .eq. 2) then

   print *, 'DriveStrength profile for alpha2'
    print *, '--------------------------'
   do i=2,n_rho_grid-1
    print *, i, ' ', &
    -(den(2,i+1)-den(2,i-1))/rmin/(rho_hat(i+1)-rho_hat(i-1))/ &
      crit_grad_n_alpha2_rho(i)
   enddo
    print *, '--------------------------'

  print *, 'Total DriveStrength profile for alpha1+alpha2'
    print *, '--------------------------'
   do i=2,n_rho_grid-1
    print *, i, ' ', &
    -(den(2,i+1)-den(2,i-1))/rmin/(rho_hat(i+1)-rho_hat(i-1))/ &
      crit_grad_n_alpha2_rho(i) &
    -(den(1,i+1)-den(1,i-1))/rmin/(rho_hat(i+1)-rho_hat(i-1))/ &
      crit_grad_n_alpha_rho(i)
   enddo
    print *, '--------------------------'

   endif


    print *, 'DriveStrength profiles based on pressure'
    print *, '--------------------------'
    print *, '--------------------------'

  rg_n_alpha_th_rho(:) = crit_grad_n_alpha_rho(:)

  n_alpha_tran_rho(:) = den(1,:)

  do i = 2,n_rho_grid-1
  rg_n_alpha_tran_rho(i) = &
      -(den(1,i+1)-den(1,i-1))/rmin/(rho_hat(i+1)-rho_hat(i-1))
  enddo
  rg_n_alpha_tran_rho(1) = rg_n_alpha_tran_rho(2)
  rg_n_alpha_tran_rho(n_rho_grid) = rg_n_alpha_tran_rho(n_rho_grid-1)
  
  do i = 2,n_rho_grid-1
   rg_p_alpha_th_rho(i)= T_alpha_equiv_rho(i)*rg_n_alpha_th_rho(i)*0.16022* &
    (1. +(T_alpha_equiv_rho(i+1)-T_alpha_equiv_rho(i-1))/T_alpha_equiv_rho(i)/ &
           (n_alpha_rho(i+1)-n_alpha_rho(i-1))*n_alpha_rho(i))
  enddo
   rg_p_alpha_th_rho(1) = rg_p_alpha_th_rho(2)
   rg_p_alpha_th_rho(n_rho_grid) = rg_p_alpha_th_rho(n_rho_grid-1)


    rg_p_alpha_tran_rho(1) = 0.
   do i = 2,n_rho_grid-1
    rg_p_alpha_tran_rho(i) = &
        -(n_alpha_tran_rho(i+1)*T_alpha_equiv_rho(i+1) &
           -n_alpha_tran_rho(i-1)*T_alpha_equiv_rho(i-1))*0.16022/rmin/(rho_hat(i+1)-rho_hat(i-1))
   enddo
   i=n_rho_grid
   rg_p_alpha_tran_rho(i) = -(n_alpha_tran_rho(i)*T_alpha_equiv_rho(i) &
           -n_alpha_tran_rho(i-1)*T_alpha_equiv_rho(i-1))*0.16022/rmin/(rho_hat(i)-rho_hat(i-1))


    print *, 'presure DriveStrength profile for alpha'
    print *, '--------------------------'
   do i=2,n_rho_grid-1
    print *, i, ' ', &
      rg_p_alpha_tran_rho(i)/rg_p_alpha_th_rho(i)
   enddo
    print *, '--------------------------'
   

  if(NBI_flag.eq.2) then
     rg_n_alpha2_th_rho(:) = crit_grad_n_alpha2_rho(:)

  n_alpha2_tran_rho(:) = den(2,:)

  do i = 2,n_rho_grid-1
  rg_n_alpha2_tran_rho(i) = &
      -(den(2,i+1)-den(2,i-1))/rmin/(rho_hat(i+1)-rho_hat(i-1))
  enddo
  rg_n_alpha2_tran_rho(1) = rg_n_alpha2_tran_rho(2)
  rg_n_alpha2_tran_rho(n_rho_grid) = rg_n_alpha2_tran_rho(n_rho_grid-1)

  do i = 2,n_rho_grid-1
   rg_p_alpha2_th_rho(i)= T_alpha2_equiv_rho(i)*rg_n_alpha2_th_rho(i)*0.16022* &
    (1. +(T_alpha2_equiv_rho(i+1)-T_alpha2_equiv_rho(i-1))/T_alpha2_equiv_rho(i)/ &
           (n_alpha2_rho(i+1)-n_alpha2_rho(i-1))*n_alpha2_rho(i))
  enddo
   rg_p_alpha2_th_rho(1) = rg_p_alpha2_th_rho(2)
   rg_p_alpha2_th_rho(n_rho_grid) = rg_p_alpha2_th_rho(n_rho_grid-1)


    rg_p_alpha2_tran_rho(1) = 0.
   do i = 2,n_rho_grid-1
    rg_p_alpha2_tran_rho(i) = &
        -(n_alpha2_tran_rho(i+1)*T_alpha2_equiv_rho(i+1) &
           -n_alpha2_tran_rho(i-1)*T_alpha2_equiv_rho(i-1))*0.16022/rmin/(rho_hat(i+1)-rho_hat(i-1))
   enddo
   i=n_rho_grid
   rg_p_alpha2_tran_rho(i) = -(n_alpha2_tran_rho(i)*T_alpha2_equiv_rho(i) &
           -n_alpha2_tran_rho(i-1)*T_alpha2_equiv_rho(i-1))*0.16022/rmin/(rho_hat(i)-rho_hat(i-1))


    print *, 'presure DriveStrength profile for alpha2'
    print *, '--------------------------'
   do i=2,n_rho_grid-1
    print *, i, ' ', &
      rg_p_alpha2_tran_rho(i)/rg_p_alpha2_th_rho(i)
   enddo
    print *, '--------------------------'

    print *, 'presure DriveStrength profile for alpha+alpha2'
    print *, '--------------------------'
   do i=2,n_rho_grid-1
    print *, i, ' ', &
      rg_p_alpha_tran_rho(i)/rg_p_alpha_th_rho(i) + &
      rg_p_alpha2_tran_rho(i)/rg_p_alpha2_th_rho(i)
   enddo
    print *, '--------------------------'
  endif !NBI_flag .eq. 2

!print gradients

  print *, 'n_grad_sd, n_grad_tran, n_grad_crit alpha'
  print *, '--------------------------------------------'
    do i=2,n_rho_grid-1
  print *, -(n_alpha_rho(i+1)-n_alpha_rho(i-1))/ &
         rmin/(rho_hat(i+1)-rho_hat(i-1)) , &
          -(den(1,i+1)-den(1,i-1))/ &
         rmin/(rho_hat(i+1)-rho_hat(i-1)), &
        crit_grad_n_alpha_rho(i)
   enddo
  print *, '--------------------------------------------'

  if(NBI_flag .eq. 2) then
  print *, 'n_grad_sd, n_grad_tran, n_grad_crit alpha2'
  print *, '--------------------------------------------'
    do i=2,n_rho_grid-1
  print *, -(n_alpha2_rho(i+1)-n_alpha2_rho(i-1))/ &
         rmin/(rho_hat(i+1)-rho_hat(i-1)) , &
          -(den(2,i+1)-den(2,i-1))/ &
         rmin/(rho_hat(i+1)-rho_hat(i-1)), &
        crit_grad_n_alpha2_rho(i)
   enddo
  print *, '--------------------------------------------'
  endif


!---------------------------------------------------
!write final results
!--------------------------------------------------
  open(unit=3,file='Alpha_time_dep_transport.out',status='replace')
    write(3,*) '--------------------------------------------------------------'
    write(3,*) 'final results from Alpha_time_dep_transport run'
    write(3,*) 'can be compared to results in Alpha_transport.out'
    write(3,*) 'can in printed in separate .out files'
    write(3,*) '--------------------------------------------------------------'
    write(3,*) 'time_run_final=',time_run,' dt=',dt
    write(3,*) 'it_max=',it_max,' it_write=', it_write
    write(3,*) '--------------------------------------------------------------'
    write(3,*) 'rho_hat'
    do i=1,n_rho_grid
     write(3,*)   rho_hat(i)
    enddo
    write(3,*) '--------------------------------------------------------------'
    write(3,*) 'alpha desults'
    write(3,*) '--------------------------------------------------------------'
    write(3,*) '--------------------------------------------------------------'

  write(3,*) '--------------------------------------------------------------'
  write(3,*) 'slowing down alpha density profile'
  write(3,*) '--------------------------------------------------------------'
  write(3,*) 'alpha density in 10**19 1/m**3'
  do i = 1,n_rho_grid
   write(3,*) n_alpha_rho(i)
  enddo
  write(3,*) '--------------------------------------------------------------'
  write(3,*) 'finished transported alpha density profile'
  write(3,*) '--------------------------------------------------------------'
  write(3,*) 'alpha density in 10**19 1/m**3'
  do i = 1,n_rho_grid
   ! write(3,*) n_alpha_tran_rho(i)
     write(3,*) den(1,i)
  enddo

  write(3,*) '--------------------------------------------------------------'
  write(3,*) 'finished D_alpha'
  write(3,*) '--------------------------------------------------------------'
  write(3,*) 'D_alpha in m**2/sec'
  do i = 1,n_rho_grid
  ! write(3,*) D_alpha(i)
    write(3,*) diffEP(1,i)
  enddo
  write(3,*) '--------------------------------------------------------------'


  write(3,*) '--------------------------------------------------------------'
  write(3,*) 'source flow'
  write(3,*) '--------------------------------------------------------------'
  write(3,*) 'in [10**19 1/m**3]*[m/sec]*m**2'

   flux(1,1)=0.
     do i=2,n_rho_grid
   flux(1,i) = V_p(i-1)/V_p(i)*flux(1,i-1) &
         + 0.5*dr/V_p(i)* &
         DS*(V_p(i)*S0_rho(i)*(1.0 - &
             0.0*den(1,i)/n_alpha_rho(i)) + &
           V_p(i-1)*S0_rho(i-1)*(1.0 - &
             0.0*den(1,i-1)/n_alpha_rho(i-1)))
     enddo
  do i = 1,n_rho_grid
   write(3,*) flux(1,i)*V_p(i)
  enddo
  write(3,*) '--------------------------------------------------------------'

  write(3,*) '--------------------------------------------------------------'
  write(3,*) 'source energy flow'
  write(3,*) '--------------------------------------------------------------'
  write(3,*) 'in MW'
! 1 X [10**19 1/m**3]*[m/sec]*m**2 
!    X [1.6022 10**(-19) Coul] X 0.001kAmp X 1000.MV = MW
  do i = 1,n_rho_grid
   write(3,*) flux(1,i)*V_p(i)*1.6022*E_alpha
  enddo
  write(3,*) '--------------------------------------------------------------'


  write(3,*) '--------------------------------------------------------------'
  write(3,*) 'transport flow'
  write(3,*) '--------------------------------------------------------------'
  write(3,*) 'in [10**19 1/m**3]*[m/sec]*m**2'

!   flux(1,:) = 0.
!   flux(1,1)=0.
!     do i=2,n_rho_grid
!   flux(1,i) = V_p(i-1)/V_p(i)*flux(1,i-1) &
!         + 0.5*dr/V_p(i)* &
!         (V_p(i)*S0_rho(i)*(1.0  &
!             - 1.0*den(1,i)/n_alpha_rho(i)) + &
!           V_p(i-1)*S0_rho(i-1)*(1.0  &
!             - 1.0*den(1,i-1)/n_alpha_rho(i-1)))
!     enddo

     do i=1,n_rho_grid
      net_of_sink(i) = 1.0 - SDsink*den(1,i)/(DS*n_alpha_rho(i))
     enddo
!NOTE:
!    He source sate is  S0_rho(i)*den(1,i)/n_alpha_rho(i) in [10**19 1/m**3]/sec

   flux(1,:) = 0.
   flux(1,1)=0.
     do i=2,n_rho_grid
      flux(1,i) = V_p(i-1)/V_p(i)*flux(1,i-1) + &
             0.5*dr/V_p(i)* &
         DS*(V_p(i)*S0_rho(i)*net_of_sink(i) + &
           V_p(i-1)*S0_rho(i-1)*net_of_sink(i-1))
     enddo

  do i = 1,n_rho_grid
   write(3,*) flux(1,i)*V_p(i)
  enddo
  write(3,*) '--------------------------------------------------------------'

  write(3,*) '--------------------------------------------------------------'
  write(3,*) 'transport energy flow'
  write(3,*) '--------------------------------------------------------------'
  write(3,*) 'in MW'
! 1 X [10**19 1/m**3]*[m/sec]*m**2 
!    X [1.6022 10**(-19) Coul] X 0.001kAmp X 1000.MV = MW
  do i = 1,n_rho_grid
   write(3,*) flux(1,i)*V_p(i)*1.6022*E_alpha
  enddo
  write(3,*) '--------------------------------------------------------------'

  if(NBI_flag .eq. 2) then
    write(3,*) '--------------------------------------------------------------'
    write(3,*) 'alpha2 desults'
    write(3,*) '--------------------------------------------------------------'
    write(3,*) '--------------------------------------------------------------'

  write(3,*) '--------------------------------------------------------------'
  write(3,*) 'slowing down alpha2 density profile'
  write(3,*) '--------------------------------------------------------------'
  write(3,*) 'alpha density in 10**19 1/m**3'
  do i = 1,n_rho_grid
   write(3,*) n_alpha2_rho(i)
  enddo
  write(3,*) '--------------------------------------------------------------'
  write(3,*) 'finished transported alpha2 density profile'
  write(3,*) '--------------------------------------------------------------'
  write(3,*) 'alpha2 density in 10**19 1/m**3'
  do i = 1,n_rho_grid
   ! write(3,*) n_alpha2_tran_rho(i)
     write(3,*) den(2,i)
  enddo

  write(3,*) '--------------------------------------------------------------'
  write(3,*) 'finished D_alpha2'
  write(3,*) '--------------------------------------------------------------'
  write(3,*) 'D_alpha in m**2/sec'
  do i = 1,n_rho_grid
  ! write(3,*) D_alpha2(i)
    write(3,*) diffEP(2,i)
  enddo
  write(3,*) '--------------------------------------------------------------'

  write(3,*) '--------------------------------------------------------------'
  write(3,*) 'source2 flow'
  write(3,*) '--------------------------------------------------------------'
  write(3,*) 'in [10**19 1/m**3]*[m/sec]*m**2'

   flux(2,1)=0.
     do i=2,n_rho_grid
   flux(2,i) = V_p(i-1)/V_p(i)*flux(2,i-1) &
         + 0.5*dr/V_p(i)* &
         DS*(V_p(i)*S02_rho(i)*(1.0 - &
             0.0*den(2,i)/n_alpha2_rho(i)) + &
           V_p(i-1)*S02_rho(i-1)*(1.0 - &
             0.0*den(2,i-1)/n_alpha2_rho(i-1)))
     enddo
  do i = 1,n_rho_grid
   write(3,*) flux(2,i)*V_p(i)
  enddo
  write(3,*) '--------------------------------------------------------------'

  write(3,*) '--------------------------------------------------------------'
  write(3,*) 'source2 energy flow'
  write(3,*) '--------------------------------------------------------------'
  write(3,*) 'in MW'
! 1 X [10**19 1/m**3]*[m/sec]*m**2 
!    X [1.6022 10**(-19) Coul] X 0.001kAmp X 1000.MV = MW
  do i = 1,n_rho_grid
   write(3,*) flux(2,i)*V_p(i)*1.6022*E_alpha2
  enddo
  write(3,*) '--------------------------------------------------------------'

  write(3,*) '--------------------------------------------------------------'
  write(3,*) 'transport flow2'
  write(3,*) '--------------------------------------------------------------'
  write(3,*) 'in [10**19 1/m**3]*[m/sec]*m**2'

!   flux(2,:) = 0.
!   flux(2,1)=0.
!     do i=2,n_rho_grid
!   flux(2,i) = V_p(i-1)/V_p(i)*flux(2,i-1) &
!         + 0.5*dr/V_p(i)* &
!         (V_p(i)*S02_rho(i)*(1.0  &
!            - 1.0*den(2,i)/n_alpha2_rho(i)) + &
!           V_p(i-1)*S02_rho(i-1)*(1.0  &
!            -  1.0*den(2,i-1)/n_alpha2_rho(i-1)))
!     enddo

     do i=1,n_rho_grid
      net_of_sink(i) = 1.0 - SDsink*den(2,i)/(DS*n_alpha2_rho(i))
     enddo

   flux(2,:) = 0.
   flux(2,1)=0.
     do i=2,n_rho_grid
      flux(2,i) = V_p(i-1)/V_p(i)*flux(2,i-1) + &
             0.5*dr/V_p(i)* &
         DS*(V_p(i)*S02_rho(i)*net_of_sink(i) + &
           V_p(i-1)*S02_rho(i-1)*net_of_sink(i-1))
     enddo
  do i = 1,n_rho_grid
   write(3,*) flux(2,i)*V_p(i)
  enddo
  write(3,*) '--------------------------------------------------------------'

  write(3,*) '--------------------------------------------------------------'
  write(3,*) 'transport energy flow2'
  write(3,*) '--------------------------------------------------------------'
  write(3,*) 'in MW'
! 1 X [10**19 1/m**3]*[m/sec]*m**2 
!    X [1.6022 10**(-19) Coul] X 0.001kAmp X 1000.MV = MW
  do i = 1,n_rho_grid
   write(3,*) flux(2,i)*V_p(i)*1.6022*E_alpha2
  enddo
  write(3,*) '--------------------------------------------------------------'

  endif !NBI_flag .eq. 2

  close(3)

    print *, 'Alpha_time_dep_transport done'
end subroutine Alpha_time_dep_transport
