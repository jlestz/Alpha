!-------------------------------------------------------------------------
! Alpha_use_other.f90
!
! PURPOSE:
!   common for other than input & output variables 
!
!-------------------------------------------------------------------------

module Alpha_use_other


  implicit none

!    integer, parameter :: nn = 51  !must be same as n_rho_grid

    real, dimension(:), allocatable :: n_alpha_ave_rho
    real, dimension(:), allocatable :: S0_rho
    real, dimension(:), allocatable :: bannana_rho
!2A
    real, dimension(:), allocatable :: n_alpha2_ave_rho
    real, dimension(:), allocatable :: S02_rho
    real, dimension(:), allocatable :: bannana2_rho
 
end       module Alpha_use_other
