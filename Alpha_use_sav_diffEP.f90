!-------------------------------------------------------------------------
! Alpha_use_sav_diffEP.f90
!
! PURPOSE:
! save stuff for diffEP
!
!-------------------------------------------------------------------------

module Alpha_use_sav_diffEP


  implicit none

!    integer, parameter :: nn = 51  !must be same as n_rho_grid


    real, dimension(:), allocatable :: sav_rg_n_alpha_th_rho 
    real, dimension(:), allocatable :: sav_rg_p_alpha_th_rho 

    real, dimension(:), allocatable :: sav_rg_n_alpha2_th_rho
    real, dimension(:), allocatable :: sav_rg_p_alpha2_th_rho
    real, dimension(:), allocatable :: sav_rg_p_alpha_tot_th_rho
    real, dimension(:), allocatable :: sav_rg_p_alpha_fac_rho
    real, dimension(:), allocatable :: sav_rg_p_alpha2_fac_rho

    real, dimension(:), allocatable :: sav_D_bkg_Angioni
    real, dimension(:), allocatable :: sav_C_p_alpha
 
    real, dimension(:), allocatable :: sav_D_bkg_Angioni2
    real, dimension(:), allocatable :: sav_C_p_alpha2

end       module Alpha_use_sav_diffEP
