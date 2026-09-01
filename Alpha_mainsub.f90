!---------------------------------------------------------
! Alpha_mainsub.f90
!
!
!--------------------------------------------------------
! PURPOSE:
!  compute ITER model plasma profiles from beta_N_ped of the EPED1 model
! 
!  compute the fusion alpha slowing down density and cross-over profiles 
!       as well as the alpha equivalent Maxwellian temperature
!  
!  The profiles are equally space "rho_hat" profiles (typically n_rho_grid =51 
!    points  rho_hat = [0.0, 1.0]) convenient for typical GYRO exp. profiles.
! 
!  GYRO inteprolates the rho_hat profiles onto a typically not equatlly spaced
!      r_hat = r/a profiles [r_hat_in, rho_hat_out]
!
!   INPUTS:
!     number of rho_hat grids: n_rho_grid
! 
!     ITER parameters:   
!                         beta_N_ped
!                         n14_ped pedestal density in 10**14 1/cm**3 
!                         I_p     current in MA
!                         B_t     vaccuum toroidal field in Tesla
!                         Rmaj0   vaccuum center major radius in m
!                         delR0oa R0/a  Shafranov shift parameter
!                         rmin    midplane minor radius in m
!                         kappa_1 elongation at rho_hat = 1
!                         kappa_0 elongation at rho_hat = 0
!                         delta_1 triangularity rho_hat = 1
!                         delta_0 triangularity at rho_hat = 0
!                         q_1     safety factor at rho_hat = 1
!                         q_0     safety factor at rho_hat = 0
!                         rho_q     flat central q control radius
!                         eps_q   minimum core s_hat control
!                         Zeff    Zeff for ei collisions
!                         M_DT    average DT proton mass atomic number  (2.5)
!                         E_alpha alpha injection energy 3.5Mev
!
!    Profile peaking parameters
!                         Fpeak_T   T(0)/T_ped
!                         Fpeak_n   n(0)/n_ped 
!                         Fpeak_ei  Te(0)/Ti(0)
!
!  OUTPUT: 
!                         beta_ped_percent  beta_ped in percent 
!                         T_ped      Ti=Te pedestal temperature in keV
!                         arho       rho(a) in m (needed by GYRO Miller geo)
!
!    equillibrium profiles
!
!                         q_rho      q profile      
!                         Rmaj_rho   Rmajor profile in m
!                         rmin_rho   midplane minor radius in m
!                         kappa_rho  elongation profile
!                         delta_rho  triangularity profile
!                         
!    Plasma profiles
!                         T_i_rho    ion temperature profile in keV
!                         T_e_rho    electron temperature profile in keV
!                         n_i_rho    ion density profile in 10**19 1/m**3
!                         n_e_rho    electron density profile in 10**19 1/m**3
!
!    alpha_slowing down profles
!                         n_alpha_rho slowing down density in10**19 1/m**3
!                         T_alpha_equiv_rho  equivalent Maxwellian temperatire in keV
!                         E_c_hat_rho alpha_cross_over E_c/E_alpha 
!
!   
!---------------------------------------------------------

subroutine Alpha_mainsub


    use Alpha_use_input      !n_rho_grid,beta_N_ped,n14_ped,I_p,B_t,Rmaj0,delR0oa,rmin,
!                        !kappa_1,kappa_0,delta_1,delta_0,q_1,q_0,rho_q,eps_q,
!                        !Zeff,M_DT,E_alpha,Fpeak_T,Fpeak_n,Fpeak_ei
 

    use Alpha_use_output     !beta_ped_percent,T_ped,beta_N_glob,arho,
!                        !q_rho(i),Rmaj_rho(i),rmin_rho(i),
!                        !kappa_rho(i),delta_rho(i),T_i_rho(i),T_e_rho(i),
!                        !n_i_rho(i),n_e_rho(i),n_alpha_rho(i),T_alpha_equiv_rho(i),
!                        !E_c_hat_rho(i)



!
  implicit none

    integer :: i
    integer :: i_chk_input
 
    i_chk_input = 1

   if(i_chk_input .eq. 1) then
    print *, '      '
    print *, '----------------------------------------------------------'
    print *, '----------------------------------------------------------'

   endif

   
    call Alpha_read_input

    print *, 'n_rho_grid', n_rho_grid
!  goto 1000

    call Alpha_comp_eq_plasma

    print *, 'n_rho_grid', n_rho_grid

    call Alpha_comp_alpha_slowing 

    print *, 'in Alpha_mainsub #1', rho_hat(1:10)

    call Alpha_write_output

!    print *, 'in Alpha_mainsub #2', rho_hat(1:10)

    print *, 'Alpha_mainsub done'

1000 continue
end subroutine Alpha_mainsub
