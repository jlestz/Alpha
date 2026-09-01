!------------------------------------------------------------------
! Alpha_driver.f90
!
! PURPOSE:
!  driver calls Alpha_mainsub
!------------------------------------------------------------------

program Alpha_driver


    use Alpha_use_input      !n_rho_grid,beta_N_ped,n14_ped,I_p,B_t,Rmaj0,delR0oa,rmin,
!                        !kappa_1,kappa_0,delta_1,delta_0,q_1,q_0,rho_q,eps_q,
!                        !Zeff,M_DT,E_alpha,Fpeak_T,Fpeak_n,Fpeak_ei


    use Alpha_use_output     !beta_ped_percent,T_ped,beta_N_glob,arho,
!                        !q_rho(i),Rmaj_rho(i),rmin_rho(i),
!                        !kappa_rho(i),delta_rho(i),T_i_rho(i),T_e_rho(i),
!                        !n_i_rho(i),n_e_rho(i),n_alpha_rho(i),T_alpha_equiv_rho(i),
!                        !E_c_hat_rho(i)

    use Alpha_use_other      !n_alpha_ave_rho(i),S0_rho(i) 

    use expro

  !---------------------------------------------------------------
  implicit none
  !---------------------------------------------------------------

  !  Alpha_mainsub
    print *, '----------------------------------------'
    print *, 'Alpha_driver calling Alpha_mainsub'
    print *, '----------------------------------------'
   call Alpha_mainsub 

!   goto 1000

  ! Alpha_transport
    print *, '----------------------------------------'
    print *, 'Alpha_driver calling Alpha_transport '
    print *, '----------------------------------------'
   call Alpha_transport

!11.15.17 REW
  ! Alpha_time_dep_transport
    print *, '----------------------------------------'
    print *, 'Alpha_driver calling Alpha_time_dep_transport '
    print *, '----------------------------------------'
!   call Alpha_time_dep_transport
    print *, 'Alpha_time_dep_transport NOT called 7.20.18'

!    if (l_read_exp_profile .eq. 1) call EXPRO_alloc('.',0)

    call Alpha_allocate(0)

end program Alpha_driver
