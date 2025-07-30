subroutine units(scale_l,scale_t,scale_d,scale_v,scale_nH,scale_T2)
  use amr_commons
  use hydro_commons
  use constants, only: Mpc2cm, mH, kB, rhoc, factG_in_cgs, M_sun, Myr2sec&
  &, kpc2cm
  use cooling_module, only: X
  use hydro_parameters, only: hydro, mass_sph
  implicit none

  real(dp)::scale_nH,scale_T2,scale_t,scale_v,scale_d,scale_l
  real(dp)::scale_m
  real(dp),parameter::km2cm=1.0d+5
  logical,save::first_call=.true.
  !-----------------------------------------------------------------------
  ! Conversion factors from user units into cgs units
  ! For gravity runs, make sure that G=1 in user units.
  !-----------------------------------------------------------------------

  ! scale_d converts mass density from user units into g/cc
  scale_d = units_density
  if(cosmo) scale_d = omega_m * rhoc *(h0/100)**2 / aexp**3

  ! scale_t converts time from user units into seconds
  scale_t = units_time
  if(cosmo) scale_t = aexp**2 / (h0*1d5/Mpc2cm)
  ! Ensure G = 1 in user units for all gravity runs
  if(poisson.and..not.cosmo) scale_t = 1.0/sqrt(factG_in_cgs * scale_d)

  ! scale_l converts distance from user units into cm
  scale_l = units_length
  if(cosmo) scale_l = aexp * boxlen_ini * Mpc2cm / (h0/100)

  ! scale_v converts velocity in user units into cm/s
  scale_v = scale_l / scale_t

  ! mass conversion factor
  scale_m = scale_d * scale_l**3

  ! scale_T2 converts (P/rho) in user unit into (T/mu) in Kelvin
  scale_T2 = mH/kB * scale_v**2

  ! scale_nH converts rho in user units into nH in H/cc
  scale_nH = X/mH * scale_d

  ! output unit information; only main process
  if ((first_call).and.(myid==1)) then
    first_call = .false.
    write(*,*)
    write(*,*)
    write(*,*) "Physical Units:"
    write(*,*) "---------------"
    write(*,'(a32,1pe10.2)') "Unit density [g/cm^3]: ", scale_d
    write(*,'(a32,1pe10.2)') "Unit mass [Msun]: ", scale_m / M_sun
    write(*,'(a32,1pe10.2)') "Unit time [Myr]: ", scale_t / Myr2sec
    write(*,'(a32,1pe10.2)') "Unit length [kpc]: ", scale_l / kpc2cm
    write(*,'(a32,1pe10.2)') "Unit velocity [km/s]: ", scale_v / km2cm
    write(*,'(a32,1pe10.2)') "Unit temperature/mu [K]: ", scale_T2
    write(*,'(a32,1pe10.2)') "Unit hydrogen density [cm^-3]: ", scale_nH
    write(*,*)
    if (hydro) then
      write(*,'(a44,1pe10.2)') "Refinement > gas mass [Msun] at levelmin: ",&
      & m_refine(levelmin) * mass_sph * scale_m / M_sun
      write(*,'(a44,1pe10.2)') "Refinement > gas mass [Msun] at nlevelmax: ",&
      & m_refine(nlevelmax) * mass_sph * scale_m / M_sun
    end if
    write(*,*)
    write(*,*)
  end if

end subroutine units
