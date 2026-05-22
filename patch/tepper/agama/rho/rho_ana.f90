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
subroutine rho_ana(x,d,dx,ncell)
  use amr_parameters
  use amr_commons, only: myid,t
  use hydro_parameters
  use poisson_parameters
  use constants

  ! AGAMA patch
  use agama_utils
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

  character(len=80)::potfilename, stzfilename
  logical::pot_file_exists
  logical,save::stz_file_exists
  real(dp),dimension(1:3)::xyz
  real(dp)::dummy_dp
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


  if(verbose)write(*,*)' Entering rho_ana'

  if(agama_debug)then
     agama_verbose = .true.
  end if

  ! get factors to scale between physical units and code units
  call units(scale_l,scale_t,scale_d,scale_v,scale_nH,scale_T2)

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
  ! The following block refers to mapping AGAMA densities onto the AMR grid

  ! coordinates of potential centre IN CODE UNITS
  x_c(1:3) = 0.5d0 * boxlen

  ! Assign density to grid
  do i=1,ncell

    ! position in physical units (kpc)
    ! taking centre shift into account
    xyz(1:3) = (x(i,1:3) - x_c(1:3)) * scale_l / agama_scale_l

    ! density in physical units (Msun/kpc^3)
    dummy_dp = agama_density(c_obj, xyz) * dens_scale

    ! ensure positive density throughout
    ! Explanation: because density is a spline, it may be negative at times
    ! dummy_dp = MAX(0.0d0, dummy_dp)
    ! A smoother approach that -> the above when smallr->0 (now: smallr=1e-10):
    dummy_dp = 0.5d0*(dummy_dp + sqrt(dummy_dp*dummy_dp + smallr*smallr))

    ! convert from physical units to code units
    d(i) = dummy_dp * (agama_scale_m / agama_scale_l**3) / scale_d

  end do

  if (agama_debug) then
    if (myid==1) write(*,'(A,2ES20.10)') 'rho_ana min/max = ', MINVAL(d), MAXVAL(d)
    if (myid==1) write(*,*) (x(MAXLOC(d,dim=1),1:3) - x_c(1:3))
  end if

end subroutine rho_ana
