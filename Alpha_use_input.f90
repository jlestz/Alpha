!-------------------------------------------------------------------------
! Alpha_use_input.f90
!
! PURPOSE:
!   inputs
!
!-------------------------------------------------------------------------

module Alpha_use_input 

  implicit none
 
!    integer, parameter :: nn = 51  !must be same as n_rho_grid

    integer:: n_rho_grid


    real :: beta_N_ped
    real :: n14_ped
    real :: I_p
    real :: B_t
    real :: Rmaj0
    real :: delR0oa
    real :: rmin
    real :: kappa_1
    real :: kappa_0
    real :: delta_1
    real :: delta_0
    real :: q_1
    real :: q_0
    real :: rho_q
    real :: eps_q
    real :: Zeff
    real :: M_i1
    real :: Q_i1
    real :: M_i2
    real :: Q_i2
    real :: f_i1
    real :: M_beam
    real :: Q_beam
    real :: M_beam2
    real :: Q_beam2
    real :: M_DT
    real :: E_alpha
!2A
    real :: E_alpha2
    real :: Fpeak_T
    real :: Fpeak_n
    real :: Fpeak_ei
    integer :: use_t_equiv
    integer :: NBI_flag
    integer :: l_Angioni_flag
    integer :: i_EP
    integer :: i_EP_2
    integer :: l_coupled

    integer :: n_critgrad_input
    integer :: n_critgrad_input_2
    character(len=60), dimension(20) :: critgrad_file
    character(len=60), dimension(20) :: critgrad_file_2

    real, dimension(:), allocatable :: crit_grad_n_alpha_rho
    real :: crit_grad_n_alpha_rho_0
    real :: crit_grad_n_alpha2_rho_0
    real, dimension(:), allocatable :: crit_grad_n_alpha2_rho
    real, dimension(:), allocatable :: crit_grad_p_alpha_rho
    real, dimension(:), allocatable :: crit_grad_p_alpha2_rho
    real :: crit_grad_p_alpha_rho_0
    real :: crit_grad_p_alpha2_rho_0

    integer :: l_read_exp_profile
    integer :: l_critgrad_method

    real :: nbeam_scale
    real :: nbeam_scale_2

    character(len=40) :: D_bkg_file
    character(len=40) :: D_bkg_file_2

    character(len=80) :: source_file
    character(len=80) :: source_file_2
    integer :: source_method
    integer :: source_method_2
    integer :: source_write_flag=0
    integer :: source_write_flag_2=0

    real, dimension(301) :: DT_Tgrid
    real, dimension(301) :: DT_sigma_v

end module Alpha_use_input
