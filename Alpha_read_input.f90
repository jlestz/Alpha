!---------------------------------------------------------
! Alpha_read_input.f90
!
! PURPOSE:
!  Reads all input parameters from input. input.Alpha 
!
!---------------------------------------------------------

subroutine Alpha_read_input

    use Alpha_use_input      !n_rho_grid,beta_N_ped,n14_ped,I_p,B_t,Rmaj0,delR0oa,rmin,
!                        !kappa_1,kappa_0,delta_1,delta_0,q_1,q_0,rho_q,eps_q,
!                        !Zeff,M_DT,E_alpha,Fpeak_T,Fpeak_n,Fpeak_ei,NBI_flag


    use Alpha_use_output     !beta_ped_percent,T_ped,beta_N_glob,arho,
!                        !q_rho(i),Rmaj_rho(i),rmin_rho(i),
!                        !kappa_rho(i),delta_rho(i),T_i_rho(i),T_e_rho(i),
!                        !n_i_rho(i),n_e_rho(i),n_alpha_rho(i),T_alpha_equiv_rho(i),
!                        !E_c_hat_rho(i)

    use EXPRO

  !--------------------------------------
  implicit none
  !
  integer :: i_read_default
  integer :: i
  integer :: j
  integer :: jstat
  !
  character(len=80) :: comment
  character(len=5) :: header
  character(len=150) :: ALPHA_DIR
  !--------------------------------------

  i_read_default = 0 !1=yes 0=no

  !defaults
  !Kinsey et al Standard case for ITER:  
  !J.E. Kinsey, G.M. Staebler, R.E. Waltz, and R. V. Budny, Nucl. Fusion 51 92012) 083001
  !see Figs 7,11,13
  !defaults

  n_rho_grid=51
  beta_N_ped=0.92
  n14_ped=0.9
  I_p=15.
  B_t=5.3
  Rmaj0=6.2
  delR0oa=.3
  rmin=2.0
  kappa_1=1.75
  kappa_0=1.0
  delta_1=0.2
  delta_0=0.0
  q_1=3.8
  q_0=1.0
  rho_q=0.5
  eps_q=0.1
  Zeff=1.0
  M_DT=2.5
  E_alpha=3.5  
  Fpeak_T=3. !20%  3.6  Q=1->20 
  Fpeak_n=1.1 !20%  1.3  Q=10->20  1.2**4 = 2.07
  Fpeak_ei=1.1
  use_t_equiv=0
  NBI_flag=0  !0=Alpha   1=NBI
  n_critgrad_input = 1  ! Number of critical gradient files to read
  critgrad_file(:) = 'alpha_ncrit.input'

  if(i_read_default .eq. 0) open(unit=1,file='Alpha_input',status='old')

  !----------------------------------------------------------
  ! Order and variable format in Alpha_input  must match here.
  !

  l_read_exp_profile = 0

  if(i_read_default .eq. 0) then

  read(1,*) n_rho_grid
!  print *, 'n_rho_grid', n_rho_grid
  read(1,*) beta_N_ped
  read(1,*) n14_ped
  read(1,*) I_p
  read(1,*) B_t
  read(1,*) Rmaj0
  read(1,*) delR0oa
  read(1,*) rmin
  read(1,*) kappa_1
  read(1,*) kappa_0
  read(1,*) delta_1
  read(1,*) delta_0
  read(1,*) q_1
  read(1,*) q_0
  read(1,*) rho_q
  read(1,*) eps_q
  read(1,*) Zeff
  read(1,*) M_i1
  read(1,*) Q_i1
  read(1,*) M_i2
  read(1,*) Q_i2
  read(1,*) f_i1
  read(1,*) M_beam
  read(1,*) Q_beam
  read(1,*) E_alpha
  read(1,*) Fpeak_T
  read(1,*) Fpeak_n
  read(1,*) Fpeak_ei
  read(1,'(I1)') use_t_equiv
  read(1,'(I1)') NBI_flag  !0=Alpha   1=NBI 2=Alpha+1MevNBI
                           !  alpha --> n_alpha & n_s ; 1MevNBI --> n_alpha2 & n_s2, S02 etc

  read(1,'(I1)') l_Angioni_flag

  M_DT = f_i1*M_i1 + (1.0-f_i1)*M_i2                         
  print *, NBI_flag
  if (NBI_flag .ne. 0) then
    print *, 'Reading critical gradient input files.'
!    read(1,'(I1)') n_critgrad_input      ! Number and names
    read(1,*) n_critgrad_input      ! Number and names
    print *, n_critgrad_input, ' files'
    do j = 1, n_critgrad_input           ! of files containing
      print *, j
      read(1,'(A)') critgrad_file(j)    ! critcal-gradient profiles
      print *, critgrad_file(j) 
    enddo                                ! (minimum used, single species) 
    read(1,*) i_EP
    read(1,*) nbeam_scale

    read(1,'(A)',iostat=jstat) D_bkg_file
    if (jstat<0) D_bkg_file = 'null'
    read(1,'(A)',iostat=jstat) source_file
    if (jstat<0) source_file = 'null'
  endif

  if ((n_rho_grid .eq. 0) .and. (NBI_flag .eq. 2)) then
    print *, 'Reading critical gradient input files of second EP.'
    read(1,*) n_critgrad_input_2      ! Number and names
    print *, n_critgrad_input_2, ' files for second species'
    do j = 1, n_critgrad_input_2           ! of files containing
      print *, j
      read(1,'(A)') critgrad_file_2(j)    ! critcal-gradient profiles
      print *, critgrad_file_2(j)
    enddo                                ! (minimum used, single species)
    read(1,*) E_alpha2
    read(1,*) M_beam2
    read(1,*) Q_beam2
    read(1,*) i_EP_2
    read(1,*) l_coupled
    read(1,*) nbeam_scale_2

    read(1,'(A)',iostat=jstat) D_bkg_file_2
    if (jstat<0) D_bkg_file_2 = 'null'
    read(1,'(A)',iostat=jstat) source_file_2
    if (jstat<0) source_file_2 = 'null'
  endif

  close(1)

  if (n_rho_grid .eq. 0) then
!    call EXPRO_alloc('./',1)
!    EXPRO_ctrl_numeq_flag = 0
!    EXPRO_ctrl_z(:) = 1.0
    EXPRO_ctrl_quasineutral_flag = 0
!    EXPRO_ctrl_n_ion = 5     ! Thermal+Carbon++mix+beam hardcoded for now
    call expro_read('input.gacode',0)
    print *, 'Finished EXPRO read.'
    l_read_exp_profile = 1
    n_rho_grid = EXPRO_n_exp
    call Alpha_allocate(1)
    print *, n_rho_grid
!    print *, EXPRO_ne
!    print *, EXPRO_b_ref
!    print *, EXPRO_q
!    print *, EXPRO_kappa
  else
    call Alpha_allocate(1)
  endif

!  goto 1000

  if ((l_read_exp_profile .ne. 1) .and. (NBI_flag .eq. 2)) then
!rew 11/22/16
  print *, 'reading TGLF_ITERalpha_ncrit.input'

  open(unit=4,file='TGLF_ITERalpha_ncrit.input',status='old')
   do i=1,n_rho_grid
    read(4,20) crit_grad_n_alpha_rho(i) 
   enddo
  close(4)

  print *, 'n_rho_grid', n_rho_grid
  print *, 'reading TGLF_ITER_NBI_ncrit.input'

  open(unit=5,file='TGLF_ITER_NBI_ncrit.input',status='old')
   do i=1,n_rho_grid
    read(5,20) crit_grad_n_alpha2_rho(i)
   enddo
  close(5)

  print *, 'n_rho_grid', n_rho_grid
  print *, 'crit_grad_n_alpha_rho(25)=',crit_grad_n_alpha_rho(25)
  print *, 'crit_grad_n_alpha2_rho(25)=',crit_grad_n_alpha2_rho(25)

 endif

 if (l_read_exp_profile .eq. 1) then
  crit_grad_n_alpha_rho(:) = 1000.0
  crit_grad_p_alpha_rho(:) = 1000.0
  do j = 1, n_critgrad_input
    open(unit=6,file=trim(critgrad_file(j)),status='old')

    read(6,'(A5)') header

    if (header .eq. 'Densi') then
      print *, 'Density critical gradient file detected:  ', critgrad_file(j)
      if (j .eq. 1) then
        l_critgrad_method = 1
      else if (l_critgrad_method .ne. 1) then
        print *, 'Error: Critical gradient files must all be of same type.'
        print *, 'Terminating run.'
        STOP
      endif
    else if (header .eq. 'Press') then
      print *, 'Pressure critical gradient file detected:  ', critgrad_file(j)
      if (j .eq. 1) then
        l_critgrad_method = 2
      else if (l_critgrad_method .ne. 2) then
        print *, 'Error: Critical gradient files must all be of same type.'
        print *, 'Terminating run.'
        STOP
      endif
    else
      print *, 'Error: ', critgrad_file(j), ' is of unspecified type.'
      print *, 'Terminating run.'
      STOP
    endif

    select case (l_critgrad_method)
      case (1)
        do i=1,n_rho_grid
          read(6,20) crit_grad_n_alpha_rho_0
          crit_grad_n_alpha_rho(i) = min(crit_grad_n_alpha_rho(i),crit_grad_n_alpha_rho_0)
        enddo
      case (2)
        do i=1,n_rho_grid
          read(6,20) crit_grad_p_alpha_rho_0
          crit_grad_p_alpha_rho(i) = min(crit_grad_p_alpha_rho(i),crit_grad_p_alpha_rho_0)
        enddo
    end select

    close(6)
  enddo

!  print *, 'crit_grad_n_alpha_rho(25)=',crit_grad_n_alpha_rho(25)

  select case (l_critgrad_method)
    case (1)
      open(unit=1,file='Alpha_dndr_crit_composite.out',status='replace')
      do i=1,n_rho_grid
        write(1,*) crit_grad_n_alpha_rho(i)
      enddo
      close(1)
    case (2)
      open(unit=1,file='Alpha_dpdr_crit_composite.out',status='replace')
      do i=1,n_rho_grid
        write(1,*) crit_grad_p_alpha_rho(i)
      enddo
      close(1)
  end select

  if (NBI_flag .eq. 2) then
    print *, 'reading crit_grad_n_alpha.input for 2nd EP'
    crit_grad_n_alpha2_rho(:) = 1000.0
    crit_grad_p_alpha2_rho(:) = 1000.0
  do j = 1, n_critgrad_input_2
      open(unit=6,file=trim(critgrad_file_2(j)),status='old')

      read(6,'(A5)') header

      if (header .eq. 'Densi') then
        print *, 'Density critical gradient file detected:  ', critgrad_file_2(j)
        if (j .eq. 1) then
          l_critgrad_method = 1
        else if (l_critgrad_method .ne. 1) then
          print *, 'Error: Critical gradient files must all be of same type.'
          print *, 'Terminating run.'
          STOP
        endif
      else if (header .eq. 'Press') then
        print *, 'Pressure critical gradient file detected:  ', critgrad_file_2(j)
        if (j .eq. 1) then
          l_critgrad_method = 2
        else if (l_critgrad_method .ne. 2) then
          print *, 'Error: Critical gradient files must all be of same type.'
          print *, 'Terminating run.'
          STOP
        endif
      else
        print *, 'Error: ', critgrad_file_2(j), ' is of unspecified type.'
        print *, 'Terminating run.'
        STOP
      endif

      select case (l_critgrad_method)
        case (1)
          do i=1,n_rho_grid
            read(6,20) crit_grad_n_alpha2_rho_0
            crit_grad_n_alpha2_rho(i) = min(crit_grad_n_alpha2_rho(i),crit_grad_n_alpha2_rho_0)
          enddo
        case (2)
          do i=1,n_rho_grid
            read(6,20) crit_grad_p_alpha2_rho_0
            crit_grad_p_alpha2_rho(i) = min(crit_grad_p_alpha2_rho(i),crit_grad_p_alpha2_rho_0)
          enddo
      end select

      close(6)
    enddo

    select case (l_critgrad_method)
      case (1)
        open(unit=1,file='Alpha2_dndr_crit_composite.out',status='replace')
        do i=1,n_rho_grid
          write(1,*) crit_grad_n_alpha2_rho(i)
        enddo
        close(1)
      case (2)
        open(unit=1,file='Alpha2_dpdr_crit_composite.out',status='replace')
        do i=1,n_rho_grid
          write(1,*) crit_grad_p_alpha2_rho(i)
        enddo
        close(1)
    end select
  
  endif ! NBI_flag = 2, two EP species

 endif ! Experimental profile critical gradient inputs

 print *, 'Reading DT fusion cross data.'
! open(unit=1,file='/global/homes/b/bassem/gacode_add_ebass_perlmutter/Alpha/DT_sigma_v.dat',status='old')
 call get_environment_variable('ALPHA_DIR',ALPHA_DIR)
 open(unit=1,file=trim(ALPHA_DIR)//'/DT_sigma_v.dat',status='old')
 read (1,*) header
 do i = 1, 301       ! Format of these data is hardcoded at 300 grid points.
   read (1,*) DT_tgrid(i), DT_sigma_v(i)
 enddo
 close (1)

20 format(ES14.7)

  endif

1000 continue

end subroutine Alpha_read_input
