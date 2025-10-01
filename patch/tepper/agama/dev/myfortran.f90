!  Demonstrates how to create and use potentials from AGAMA library in FORTRAN.
!  Due to the absence of a native pointer type in FORTRAN, the pointer to the !  C++ object should be stored in a placeholder variable of type CHAR*8, which !  is passed as the first argument to all functions in this module.

!  IMPORTANT: AGAMA adopts units where G=1, [V] = km/s, [L]=kpc, [T] = 0.978 Gyr
!  which fix [M] = 2.33e5 Msun
!  The potential is therefore automatically given in (km/s), while the density
!  must be multiplied by the mass unit to arrive at a density in Msun.kpc^3

    program example
    implicit none
!  This is not an actual string, but a placeholder to keep the pointer to the C++ object
    character(len=8)::c_obj5
    character(len=80)::filename
!  Functions provided by the AGAMA library
    double precision::agama_potential, agama_density
!  Local variables
    integer::i
    double precision::xyz(3),dummy_dp,scale_m
    logical::success
    success = .true.
    scale_m = 2.33d5

!  Constructing a potential from parameters stored in an INI file
!    call agama_initfromfile(c_obj5, 'BT08.ini')
    filename='hbd_26_comp1_snap00001_pot_skip20.dat'
    call agama_initfromfile(c_obj5, TRIM(filename))
    xyz(1)= 1d0
    xyz(2)= 1d0
    xyz(3)= 1d0
    do i=1,100000 ! evaluate a large number of times
       dummy_dp=agama_potential(c_obj5, xyz)
       dummy_dp=agama_density(c_obj5, xyz)
    enddo
    write(*,*)'Potential=-(', SQRT(ABS(agama_potential(c_obj5, xyz))),' km/s)^2'
    write(*,*)'Density=', agama_density(c_obj5, xyz)*scale_m,' Msun/kpc^3'

    if(success) then
        print*, 'ALL TESTS PASSED'
    else
        print*, 'SOME TESTS FAILED'
    endif
    end program example
