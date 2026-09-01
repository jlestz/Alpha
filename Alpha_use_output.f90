!-------------------------------------------------------------------------
! Alpha_use_output.f90
!
! PURPOSE:
!    output profiled
!
!-------------------------------------------------------------------------

module Alpha_use_output


  implicit none

!    integer, parameter :: nn = 51  !must be same as n_rho_grid

    real :: beta_ped_percent
    real :: T_ped
    real :: beta_N_glob
    real :: arho

    real, dimension(:), allocatable :: E_alpha_rho
    real, dimension(:), allocatable :: E_alpha2_rho

    real, dimension(:), allocatable :: rho_hat
    real, dimension(:), allocatable :: q_rho
    real, dimension(:), allocatable :: Rmaj_rho
    real, dimension(:), allocatable :: rmin_rho
    real, dimension(:), allocatable :: kappa_rho
    real, dimension(:), allocatable :: delta_rho
    real, dimension(:), allocatable :: T_i_rho
    real, dimension(:), allocatable :: T_e_rho
    real, dimension(:), allocatable :: n_i_rho
    real, dimension(:), allocatable :: n_e_rho

    real, dimension(:), allocatable :: n_alpha_rho
    real, dimension(:), allocatable :: T_alpha_equiv_rho
    real, dimension(:), allocatable :: E_c_hat_rho
    real, dimension(:), allocatable :: E_alpha_grid
!2A
    real, dimension(:), allocatable :: n_alpha2_rho
    real, dimension(:), allocatable :: T_alpha2_equiv_rho
    real, dimension(:), allocatable :: E_c_hat2_rho
    real, dimension(:), allocatable :: E_alpha2_grid

!!    real :: n_alpha_ave_rho
 
end       module Alpha_use_output
