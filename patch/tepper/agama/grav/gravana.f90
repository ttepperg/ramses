module agama_utils

    use amr_parameters, only: dp
    implicit none

contains

    function find_closest_t(tarr,ti) result(iti)
    ! returns the index of the element in array tarr
    ! closest to ti
      real(dp),intent(in)::tarr(:),ti
      integer::iti
      iti = minloc(abs(tarr-ti),dim=1)
    end function find_closest_t

    !------------------------------------------------------------
    ! Helper function: check if all characters in a string are digits
    !------------------------------------------------------------
    pure function is_all_digits(s) result(res)
      character(len=*), intent(in) :: s
      logical :: res
      integer :: j
      res = .true.
      do j = 1, len(s)
        if (s(j:j) < '0' .or. s(j:j) > '9') then
          res = .false.
          return
        end if
      end do
    end function is_all_digits

    !------------------------------------------------------------
    ! set_pot_filename
    ! Replaces the first occurrence of a numeric substring in a filename
    ! immediately following a specified prefix.
    !
    ! Arguments:
    !   template : input filename (character)
    !   prefix   : prefix string preceding the numeric field
    !   width    : number of digits to replace (integer)
    !   newnum   : integer to insert
    !
    ! Returns:
    !   new filename (character(len=template))
    !------------------------------------------------------------
    function set_pot_filename(template, prefix, width, newnum) result(newfile)
      implicit none
      character(len=*), intent(in) :: template
      character(len=*), intent(in) :: prefix
      integer,          intent(in) :: width
      integer,          intent(in) :: newnum
      character(len=len(template))  :: newfile

      integer :: i, j, k, nzeros, prefix_len
      character(len=32) :: numstr
      logical :: found_pattern

      newfile = template
      found_pattern = .false.
      prefix_len = len_trim(prefix)

      ! Scan for prefix followed by numeric substring of given width
      do i = 1, len(template) - prefix_len - width + 1
        if (template(i:i+prefix_len-1) == prefix .and. &
            is_all_digits(template(i+prefix_len:i+prefix_len+width-1))) then
          j = i + prefix_len          ! start of numeric substring
          k = j + width - 1           ! end of numeric substring
          found_pattern = .true.
          exit
        end if
      end do

      if (.not. found_pattern) then
        print *, "Error: no pattern found for prefix '", prefix, "' with width ", width
        stop
      end if

      ! Write new number into buffer
      write(numstr, '(I0)') newnum

      ! Zero-pad to desired width
      nzeros = max(0, width - len_trim(numstr))
      numstr = repeat('0', nzeros) // trim(numstr)

      ! Copy into template slice
      newfile(j:k) = numstr(1:width)

    end function set_pot_filename

end module agama_utils

!#########################################################
!#########################################################
!#########################################################
!#########################################################
subroutine gravana(x,f,dx,ncell)
  use amr_parameters
  use amr_commons, only: myid,t
  use poisson_parameters
  use constants

  ! AGAMA patch
  use agama_utils
  use agama_commons

  implicit none
  integer ::ncell                         ! Size of input arrays
  real(dp)::dx                            ! Cell size
  real(dp),dimension(1:nvector,1:ndim)::f ! Gravitational acceleration
  real(dp),dimension(1:nvector,1:ndim)::x ! Cell center position.
  !================================================================
  ! This routine computes the acceleration using analytical models.
  ! x(i,1:ndim) are cell center position in [0,boxlen] (user units).
  ! f(i,1:ndim) is the gravitational acceleration in user units.
  ! Only if there is no self-gravity
  !================================================================
  integer::idim,i
  real(dp)::scale_l,scale_t,scale_d,scale_v,scale_nH,scale_T2
  real(dp)::scale_m
  real(dp),parameter::km2cm=1.0d5

  ! ----------- AGAMA patch BELOW this line --------------------

  ! Functions provided by the AGAMA library
  real(dp)::agama_potential ! not used here, but to remind me of it
  real(dp)::agama_potforce

  ! The following is not an actual string, but a placeholder to keep the pointer to the C++ object
  character(len=8),save::c_obj     ! <- VERY important to save

  character(len=80)::potfilename, stzfilename
  logical::pot_file_exists
  logical,save::stz_file_exists
  real(dp),dimension(1:3)::xyz
  real(dp)::dummy_dp
  real(dp),dimension(1:3)::force_vector    ! agama force vector
  real(dp)::agama_force_factor, ramses_force_factor ! unit conversion factors
  real(dp),dimension(1:3)::x_c      ! potential's centre coordinates
  logical::read_pot_file=.true.     ! ensure pot file read at t_step
  logical::read_stz_file=.true.     ! ensure stz file read at t_step
  real(dp),save::t_prev=0.          ! file read at this previous time
  real(dp)::t_input=0.,t_output=0.  ! convenience variables

  integer::ios                      ! control file read
  integer::count_snapshots          ! store the number of lines in stzfilename
  integer,allocatable,dimension(:),save::snapnum_agama            ! snapshot number corresponding to pot file
  real(dp),allocatable,dimension(:),save::t_agama,dt_agama        ! time [step] corresponding to snapshot
  real(dp):: t_Myr                  ! simulation time in Myr
  integer::t_index                  ! index of stz table time closest to t
  integer,save::t_index_prev=-1     ! previous index


  if(verbose)write(*,*)' Entering grav_ana'

  if(agama_debug)then
     agama_verbose = .true.
  end if

  ! get factors to scale between physical units and code units
  call units(scale_l,scale_t,scale_d,scale_v,scale_nH,scale_T2)
  ! mass conversion factor
  scale_m = scale_d * scale_l**3

  !----------------------------------------------------
  ! The following block refers to the snapshot-time-redshift table

  ! File needs to be read only *once* per [re]start (at t=0 or nrestart!=0)
  read_stz_file = read_stz_file .and. ( (t==0) .or. (nrestart.ne.0))
  if(read_stz_file) then

    ! update I/O flag
    read_stz_file = .false.

    if(TRIM(initfile(levelmin)).NE.' ')then
      stzfilename=TRIM(initfile(levelmin))//'/'//TRIM(agama_stz_file)
      INQUIRE(FILE=stzfilename,EXIST=stz_file_exists)
      if(.not.stz_file_exists.and.myid==1) then
        write(*,*) '  '//TRIM(stzfilename)//' not found'
        write(*,*) '  Will consider only one potential file (if existent):'
        write(*,*) '  '//TRIM(agama_pot_file)
      end if
    end if

    if (stz_file_exists) then

        ! output info
        if(myid==1.and.agama_verbose)then
          write(*,*) '  Reading AGAMA STZ file:'
          write(*,*) '  ', TRIM(agama_stz_file)
          write(*,*) '  At t=', t
        end if

        ! determine how many snapshots
        count_snapshots = 0
        open(1234,file=stzfilename,form='formatted')
          read(1234,*) ! skip column header
          do
            read(1234,*,iostat=ios) dummy_dp
            if (ios < 0) exit      ! EOF
            if (ios > 0) error stop "I/O error"
            count_snapshots = count_snapshots + 1
          end do
        close(1234)
        if(myid==1.and.agama_verbose) &
          & write(*,*)'  Found ',count_snapshots,' snapshots'

        ! Allocate arrays
        allocate(snapnum_agama(1:count_snapshots))
        allocate(t_agama(1:count_snapshots))
        allocate(dt_agama(1:count_snapshots))

        count_snapshots = 0
        open(1234,file=stzfilename,form='formatted')
          read(1234,*) ! skip column header
          do
            read(1234,*,iostat=ios) snapnum_agama(count_snapshots+1), t_agama(count_snapshots+1), dt_agama(count_snapshots+1)
            if(ios/=0) exit
            count_snapshots = count_snapshots + 1
          end do
        close(1234)

    end if ! file exists

  end if

  !----------------------------------------------------
  ! The following block refers to the AGAMA INI potential file

  ! Determine the read file frequency and delete the previously allocated potential
  if (.not.stz_file_exists) then ! do only if no snap|time|z table available

    t_input = t_prev+t_step
    if(t.gt.t_input) then
      read_pot_file = .true.   ! update I/O flag
      t_output = t_prev        ! save previous output time
      t_prev = t               ! save current time
    end if

    ! define potential filename
    potfilename=TRIM(initfile(levelmin))//'/'//TRIM(agama_pot_file)

  else

    ! Determine which potential file to read depending on t
    t_Myr = t * (scale_t / Myr2sec)

    ! IMPORTANT: May need to use dt_agama rather than t_agama -> TO CHECK
    t_index = find_closest_t(t_agama,t_Myr)

    if (ABS(t_Myr - t_agama(t_index)).le.t_diff) then

        ! Define potential filename
        ! NB
        ! - 'agama_pot_file' is used as as string template and is *always* required
        ! - adjust 'snap' and the width (5) if necessary; the settings below are appropriate for a template filename of the form "str1_snap00042_str2.str3"
        potfilename = &
        &set_pot_filename(agama_pot_file, "snap", 5, snapnum_agama(t_index))
        potfilename = TRIM(initfile(levelmin))//'/'//TRIM(potfilename)

        if (t_index.gt.t_index_prev) then

          if (myid==1) then
            write(*,*)'  Reading next potential file at'
            write(*,*)'  sim time | snap time | snap number'
            write(*,*) t_Myr, '| ', t_agama(t_index), '| ', snapnum_agama(t_index)
          end if

          t_index_prev = t_index   ! save current index
          read_pot_file = .true.   ! update I/O flag
          t_output = t_prev        ! save previous output time
          t_prev = t               ! save current time

        end if ! (t_index.gt.t_index_prev)

    end if ! (ABS(t_Myr - t_agama(t_index)).le.t_diff)

  end if ! (.not.stz_file_exists)

  ! File needs to be read only *once* per time step; the latter can be variously defined as 't_step' or determined from the 'agama_stz_file'
  !
  if(read_pot_file) then

    ! Constructing an AGAMA potential from parameters stored in an INI file
    if(TRIM(initfile(levelmin)).NE.' ')then
      INQUIRE(FILE=potfilename,EXIST=pot_file_exists)
      if(.not.pot_file_exists) then
        if(myid==1) write(*,*) TRIM(potfilename)," not found"
        call clean_stop
      end if
    end if

    ! output info
    if(myid==1.and.agama_verbose)then
      write(*,*) '  Reading AGAMA INI file:'
      write(*,*) '  ', TRIM(potfilename)
      write(*,*) '  At t=', t, ' (previous t=', t_output, ')'
    end if

    ! VERY important to TRIM the file name when passing to routine:
    call agama_initfromfile(c_obj, TRIM(potfilename))

    if(myid==1.and.agama_verbose) write(*,*) '  DONE'

  end if

  ! update I/O flag: file will be read *only* at each t_step
  read_pot_file=.false.

  !----------------------------------------------------
  ! The following block refers to mapping AGAMA forces onto the AMR grid

  ! coordinates of potential centre IN CODE UNITS
  x_c(1:3) = 0.5d0 * boxlen

  ! Assign density to grid
  do i=1,ncell

    ! position in physical units (kpc)
    ! taking centre shift into account
    xyz(1:3) = (x(i,1:3) - x_c(1:3)) * scale_l / agama_scale_l

    ! DEVELOPMENT
    ! force in physical units (Msun/kpc)(km/s)^2
    dummy_dp = agama_potforce(c_obj, xyz, force_vector) * force_scale
    !if(myid==1)write(*,*)'Force [c.u. agama] = ', force_vector
    force_vector = force_vector * force_scale
    !if(myid==1)write(*,*)'Force [p.u.] = ', force_vector

    ! DEVELOPMENT
    ! convert from physical units to code units
    agama_force_factor = (agama_scale_m / agama_scale_l) * ( (agama_scale_v * km2cm)**2 )
    ramses_force_factor = (scale_m / scale_l) * (scale_v**2)
    force_vector(1:3) = force_vector(1:3) * (agama_force_factor / ramses_force_factor)

    !if(myid==1)write(*,*)agama_scale_m, agama_scale_l, agama_scale_v
    !if(myid==1)write(*,*)scale_m, scale_l, scale_v
    !if(myid==1)write(*,*)agama_force_factor, ramses_force_factor

    !if(myid==1)write(*,*)'Force = [c.u. ramses] ', force_vector

    ! Assign force to grid
    f(i,1:3)=force_vector(1:3)

  end do

! CHANGE density for force
!  if (agama_debug) then
!    if (myid==1) write(*,'(A,2ES20.10)') 'grav_ana min/max = ', MINVAL(d), MAXVAL(d)
!    if (myid==1) write(*,*) (x(MAXLOC(d,dim=1),1:3) - x_c(1:3))
!  end if


end subroutine gravana
!#########################################################
!#########################################################
!#########################################################
!#########################################################
!
! TTG: DO NOT TOUCH THIS ROUTINE. IT IS NEEDED BY
!   poisson/phi_fine_cg.f90
!   poisson/boundary_potential.f90
!
subroutine phi_ana(rr,pp,ngrid)
  use amr_commons
  use poisson_commons
  use constants, only: twopi
  implicit none
  integer::ngrid
  real(dp),dimension(1:nvector)::rr,pp
  ! -------------------------------------------------------------------
  ! This routine set up boundary conditions for fine levels.
  ! -------------------------------------------------------------------

  integer :: i
  real(dp):: fourpi

  fourpi=2*twopi

  do i=1,ngrid
#if NDIM==1
     pp(i)=multipole(1)*fourpi/2*rr(i)
#elif NDIM==2
     pp(i)=multipole(1)*2*log(rr(i))
#elif NDIM==3
     pp(i)=-multipole(1)/rr(i)
#endif
  end do
end subroutine phi_ana
