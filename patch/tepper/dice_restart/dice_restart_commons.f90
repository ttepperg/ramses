module dice_restart_commons
  use amr_commons
  use hydro_commons

  ! misc
  integer::nrestart2               = 0
  real(dp)::IG_rho_restart         = 1.0D-5
  real(dp)::IG_T2_restart          = 1.0D7
  real(dp)::IG_metal_restart       = 0.01
  real(dp),dimension(1:3)::ic_center_restart = (/ 0.0, 0.0, 0.0 /)
  integer,dimension(1:100)::restart_vars=0
  integer::nvar_min
  logical::dice_restart_init=.false.
  logical::add_dice_ic=.false. ! load additional DICE ICs at restart
  real(dp)::restart_boxlen
  real(dp)::restart_unit_t
  real(dp)::restart_unit_l
  real(dp)::restart_unit_d

end module dice_restart_commons
