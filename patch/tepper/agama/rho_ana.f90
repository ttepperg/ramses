!#########################################################
!#########################################################
!#########################################################
!#########################################################
subroutine rho_ana(x,d,dx,ncell)
  use amr_parameters
  use amr_commons, only: myid,t
  use hydro_parameters
  use poisson_parameters
  use constants
  use agama_commons
  implicit none
  integer ::ncell                         ! Number of cells
  real(dp)::dx                            ! Cell size
  real(dp),dimension(1:nvector)     ::d   ! Density
  real(dp),dimension(1:nvector,1:ndim)::x ! Cell center position.
  !================================================================
  ! This routine generates analytical Poisson source term.
  ! Positions are in user units:
  ! x(i,1:3) are in [0,boxlen]**ndim.
  ! d(i) is the density field in user units.
  !================================================================
  integer::i
  real(dp)::scale_l,scale_t,scale_d,scale_v,scale_nH,scale_T2,scale_m

  ! ----------- AGAMA patch BELOW this line --------------------

  ! Functions provided by the AGAMA library
  real(dp)::agama_potential ! not used here, but to remind me of it
  real(dp)::agama_density

  ! The following is not an actual string, but a placeholder to keep the pointer to the C++ object
  character(len=8),save::c_obj5     ! <- VERY important to save

  character(len=89)::filename,infile
  logical::file_exists
  real(dp),dimension(1:3)::xyz
  real(dp)::dummy_dp, mass_factor
  real(dp),dimension(1:3)::x_c      ! potential's centre coordinates

  logical,save::read_file=.true.    ! affects input file read
  real(dp),save::t_prev=0.       ! affects output info

  ! Factor to transform density in AGAMA units (G=1, V=1km/s, L=1kpc, T~1Gyr) to Msun/kpc^3; potential is direclty given in (km/s)^2
  mass_factor = 2.33d5

  ! set I/O flags
  if(t_prev < t) then
     !read_file=.true.
  endif

  if(myid==1.and.read_file)then
     write(*,*)'Entering rho_ana at t=',t,' (previous: ',t_prev,')'
     t_prev = t             ! save current time
  endif

  ! get factors to scale between physical units and code units
  call units(scale_l,scale_t,scale_d,scale_v,scale_nH,scale_T2)

  ! coordinates of potential centre IN CODE UNITS
  x_c(1:3) = ( INT( 0.5d0 * boxlen / dx ) ) * dx

  ! File needs to be read only *once* per *main* (coarse) time step
  if(read_file)then

  ! Constructing an AGAMA potential from parameters stored in an INI file
      if(TRIM(initfile(levelmin)).NE.' ')then
         filename=TRIM(initfile(levelmin))//'/'//TRIM(agama_file)
         INQUIRE(FILE=filename,EXIST=file_exists)
         if(.not.file_exists) then
            if(myid==1) write(*,*) TRIM(filename)," not found"
            call clean_stop
         endif
      endif

      if(myid==1.and.agama_verbose)then
         write(*,*) "Reading AGAMA INI file: "
         write(*,*) TRIM(filename)
      end if

      ! VERY important to TRIM the file name when passing to routine:
      call agama_initfromfile(c_obj5, TRIM(filename))

      if(myid==1.and.agama_verbose) write(*,*) 'DONE'

  endif

  ! update I/O flag
  read_file=.false.

  ! Assign density to grid
  do i=1,ncell

    ! position in physical units (kpc)
    ! MUST take centre shift into account -> TO DO!
    xyz(1:3) = (x(i,1:3) - x_c(1:3)) * scale_l / agama_scale_l

    ! density in physical units (Msun/kpc^3)
    dummy_dp = agama_density(c_obj5, xyz) * mass_factor

    ! convert from physical units to code units
    d(i) = dummy_dp * (agama_scale_m / agama_scale_l**3) / scale_d

  end do

  ! ----------- AGAMA patch ABOVE this line --------------------

end subroutine rho_ana
