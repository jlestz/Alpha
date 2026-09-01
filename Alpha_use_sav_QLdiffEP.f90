!-------------------------------------------------------------------------
! Alpha_use_sav_QLdiffEP.f90
!
! PURPOSE:
! save stuff for  Alpha_QLdiffusivity
!
!-------------------------------------------------------------------------

module Alpha_use_sav_QLdiffEP


  implicit none

!    integer, parameter :: nn = 51  !must be same as n_rho_grid


!   minimum critical gardients
!    real, dimension(:), allocatable :: sav_rg_n_alpha_th_rho 
!    real, dimension(:), allocatable :: sav_rg_n_alpha2_th_rho

!   slowing down
    real, dimension(:), allocatable :: sav_rg_n_alpha_rho  ![10**19/m**3]/m
    real, dimension(:), allocatable :: sav_Ln_alpha_rho    ![m]
    real, dimension(:), allocatable :: sav_LT_alpha_rho    ![m]

    real, dimension(:), allocatable :: sav_rg_n_alpha2_rho ![10**19/m**3]/m
    real, dimension(:), allocatable :: sav_Ln_alpha2_rho   ![m]
    real, dimension(:), allocatable :: sav_LT_alpha2_rho   ![m]

!RQL runing mode intensites

    real, dimension(:,:), allocatable :: E_hat  !function of time_run
    real, dimension(:,:), allocatable :: Z_hat
    real :: time_called

!RQL  from TGLFEP

    real, dimension(:,:,:), allocatable :: sav_gamma_star_hat
     !gB units
    real, dimension(:,:,:), allocatable :: sav_diff_star
     ! [m**2/sec]
    real, dimension(:,:,:), allocatable :: sav_rg_n_crit_star
     !same units as  rg_n_alpha_th_rho,  rg_n_alpha2_th_rho
     ! [10**19/m**3]/m

end       module Alpha_use_sav_QLdiffEP
