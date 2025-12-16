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

  ! AGAMA patch
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
  character(len=8),save::c_obj     ! <- VERY important to save

  character(len=80)::potfilename
  logical::file_exists
  real(dp),dimension(1:3)::xyz
  real(dp)::dummy_dp
  real(dp),dimension(1:3)::x_c      ! potential's centre coordinates
  logical::read_file=.true.         ! ensure file read at t=0
  real(dp),save::t_prev=0.          ! file read at this previous time
  real(dp)::t_input=0.,t_output=0.  ! convenience variables

  if(verbose)write(*,*)' Entering rho_ana'

  if(agama_debug)then
     agama_verbose = .true.
  endif

  ! get factors to scale between physical units and code units
  call units(scale_l,scale_t,scale_d,scale_v,scale_nH,scale_T2)


  ! the following determines the read file frequency and deletes the previously allocated potential
  t_input = t_prev+t_step
  if(t.gt.t_input) then
     read_file=.true.       ! update I/O flag
     t_output = t_prev      ! save previous output time
     t_prev = t             ! save current time
  endif

  ! coordinates of potential centre IN CODE UNITS
  x_c(1:3) = 0.5d0 * boxlen

  ! File needs to be read only *once* per *main* (coarse) time step
  if(read_file) then

     ! Constructing an AGAMA potential from parameters stored in an INI file
     if(TRIM(initfile(levelmin)).NE.' ')then
        potfilename=TRIM(initfile(levelmin))//'/'//TRIM(agama_pot_file)
        INQUIRE(FILE=potfilename,EXIST=file_exists)
        if(.not.file_exists) then
           if(myid==1) write(*,*) TRIM(potfilename)," not found"
           call clean_stop
        endif
     endif

     ! output info
     if(myid==1.and.agama_verbose)then
        write(*,*) '  Reading AGAMA INI file:'
        write(*,*) '  ', TRIM(potfilename)
        write(*,*) '  At t=', t, ' (previous t=', t_output, ')'
     end if

     ! VERY important to TRIM the file name when passing to routine:
     call agama_initfromfile(c_obj, TRIM(potfilename))

     if(myid==1.and.agama_verbose) write(*,*) '  DONE'

  endif

  ! update I/O flag
  read_file=.false.

  ! Assign density to grid
  do i=1,ncell

    ! position in physical units (kpc)
    ! MUST take centre shift into account -> TO DO!
    xyz(1:3) = (x(i,1:3) - x_c(1:3)) * scale_l / agama_scale_l

    ! density in physical units (Msun/kpc^3)
    dummy_dp = agama_density(c_obj, xyz) * dens_scale

    ! convert from physical units to code units
    d(i) = dummy_dp * (agama_scale_m / agama_scale_l**3) / scale_d

  end do

  ! ----------- AGAMA patch ABOVE this line --------------------

end subroutine rho_ana
