!--------------------------------------------------
! Alpha_allocate.f90
!
! Allocate routine for global radial arrays.
!
!--------------------------------------------------     

subroutine Alpha_allocate(i_control_flag)

  use Alpha_use_input
  use Alpha_use_output
  use Alpha_use_other
  use Alpha_use_sav_QLdiffEP
  use Alpha_use_sav_diffEP

  integer, intent(in) :: i_control_flag
  integer :: nn

  nn = n_rho_grid

  if (i_control_flag .eq. 1) then

    ! ----- Inputs ---------      
    allocate(crit_grad_n_alpha_rho(nn))
    allocate(crit_grad_n_alpha2_rho(nn))
    allocate(crit_grad_p_alpha_rho(nn))
    allocate(crit_grad_p_alpha2_rho(nn))

    ! ----- Outputs --------
    allocate(E_alpha_rho(nn))
    allocate(E_alpha2_rho(nn))
    allocate(rho_hat(nn))
    allocate(q_rho(nn))
    allocate(Rmaj_rho(nn))
    allocate(rmin_rho(nn))
    allocate(kappa_rho(nn))
    allocate(delta_rho(nn))
    allocate(T_i_rho(nn))
    allocate(T_e_rho(nn))
    allocate(n_i_rho(nn))
    allocate(n_e_rho(nn))
    allocate(n_alpha_rho(nn))
    allocate(T_alpha_equiv_rho(nn))
    allocate(E_c_hat_rho(nn))
    allocate(E_alpha_grid(nn))
    allocate(n_alpha2_rho(nn))
    allocate(T_alpha2_equiv_rho(nn))
    allocate(E_c_hat2_rho(nn))
    allocate(E_alpha2_grid(nn))

    ! ----- Other ----------
    allocate(n_alpha_ave_rho(nn))
    allocate(S0_rho(nn))
    allocate(bannana_rho(nn))
    allocate(n_alpha2_ave_rho(nn))
    allocate(S02_rho(nn))
    allocate(bannana2_rho(nn))

    ! ----- Saves ----------
    allocate(sav_rg_n_alpha_th_rho(nn))
    allocate(sav_rg_n_alpha2_th_rho(nn))

    allocate(sav_rg_n_alpha_rho(nn))
    allocate(sav_Ln_alpha_rho(nn))
    allocate(sav_LT_alpha_rho(nn))

    allocate(sav_rg_n_alpha2_rho(nn))
    allocate(sav_Ln_alpha2_rho(nn))
    allocate(sav_LT_alpha2_rho(nn))

    allocate(E_hat(15,nn))
    allocate(Z_hat(15,nn))
    allocate(sav_gamma_star_hat(2,15,nn))
    allocate(sav_diff_star(2,15,nn))
    allocate(sav_rg_n_crit_star(2,15,nn))
!    allocate(sav_rg_n_alpha_th_rho(nn))
    allocate(sav_rg_p_alpha_th_rho(nn))
!    allocate(sav_rg_n_alpha2_th_rho(nn))
    allocate(sav_rg_p_alpha2_th_rho(nn))
    allocate(sav_rg_p_alpha_tot_th_rho(nn))
    allocate(sav_rg_p_alpha_fac_rho(nn))
    allocate(sav_rg_p_alpha2_fac_rho(nn))
    
    allocate(sav_D_bkg_Angioni(nn))
    allocate(sav_C_p_alpha(nn))
    allocate(sav_D_bkg_Angioni2(nn))
    allocate(sav_C_p_alpha2(nn))

  else

    ! ----- Inputs ---------
    deallocate(crit_grad_n_alpha_rho)
    deallocate(crit_grad_n_alpha2_rho)
    deallocate(crit_grad_p_alpha_rho)
    deallocate(crit_grad_p_alpha2_rho)

    ! ----- Outputs --------
    deallocate(E_alpha_rho)
    deallocate(E_alpha2_rho)
    deallocate(rho_hat)
    deallocate(q_rho)
    deallocate(Rmaj_rho)
    deallocate(rmin_rho)
    deallocate(kappa_rho)
    deallocate(delta_rho)
    deallocate(T_i_rho)
    deallocate(T_e_rho)
    deallocate(n_i_rho)
    deallocate(n_e_rho)
    deallocate(n_alpha_rho)
    deallocate(T_alpha_equiv_rho)
    deallocate(E_c_hat_rho)
    deallocate(E_alpha_grid)
    deallocate(n_alpha2_rho)
    deallocate(T_alpha2_equiv_rho)
    deallocate(E_c_hat2_rho)
    deallocate(E_alpha2_grid)
    deallocate(E_hat)
    deallocate(Z_hat)

    ! ----- Other ----------
    deallocate(n_alpha_ave_rho)
    deallocate(S0_rho)
    deallocate(bannana_rho)
    deallocate(n_alpha2_ave_rho)
    deallocate(S02_rho)
    deallocate(bannana2_rho)
          
    ! ----- Saves ----------
    deallocate(sav_gamma_star_hat)
    deallocate(sav_diff_star)
    deallocate(sav_rg_n_crit_star)
    deallocate(sav_rg_n_alpha_th_rho)
    deallocate(sav_rg_p_alpha_th_rho)
    deallocate(sav_rg_n_alpha2_th_rho)
    deallocate(sav_rg_p_alpha2_th_rho)
    deallocate(sav_rg_p_alpha_tot_th_rho)
    deallocate(sav_rg_p_alpha_fac_rho)
    deallocate(sav_rg_p_alpha2_fac_rho)

    deallocate(sav_D_bkg_Angioni)
    deallocate(sav_C_p_alpha)
    deallocate(sav_D_bkg_Angioni2)
    deallocate(sav_C_p_alpha2)

  endif

end subroutine
