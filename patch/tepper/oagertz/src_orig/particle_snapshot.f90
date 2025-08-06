module particle_snapshot
  use amr_commons, only: myid, verbose
  use amr_parameters, only: dp, i8b, ndim, imov, nvector, poisson, particle_snapshot_file
  use pm_commons, only: xp, vp, partp, is_tracer, is_DM, is_gas_tracer, npartmax, typep, idp, part2int, itmpp, FAM_DM
  use hydro_commons, only: uold, smallr, nener, gamma
  use poisson_commons, only: phi
  use tracer_utils, only: find_grid_containing
  use dump_utils, only: generic_dump, dump_header_info, dim_keys

  implicit none

  integer(i8b), allocatable, dimension(:) :: particle_ids
  integer :: npart_subset
  integer :: selected_particle_type = -999  ! Dummy value to prevent any error

contains
  subroutine load_ids_subset()
    integer :: ptype

    ! This will only read the file once
    if (.not. allocated(particle_ids) .and. trim(particle_snapshot_file) /= "") then
      call read_particle_ids(&
        particle_snapshot_file, &
        npart_subset, &
        particle_ids, &
        selected_particle_type&
      )
    end if
  end subroutine load_ids_subset


  ! Perform selection, returns != 0 if a particle is selected
  subroutine select_particle(ipart, npart_selected)
    integer, intent(in) :: ipart
    integer, intent(inout) :: npart_selected

    itmpp(ipart) = 0
    if (is_tracer(typep(ipart))) then
      itmpp(ipart) = 1
    else if (typep(ipart)%family == selected_particle_type) then
      ! The binary search returns /= 0 iff the value is found
      itmpp(ipart) = binary_search(idp(ipart), particle_ids, npart_subset)
    end if

    if (itmpp(ipart) > 0) then
      npart_selected = npart_selected + 1
    end if
  end subroutine select_particle

  subroutine output_particle_snapshot(filename)

    character(len=1024), intent(in) :: filename

    integer :: idim, j, jcache, ipart, npart_selected
    real(dp), allocatable, dimension(:) :: xdp_out
    integer(i8b), allocatable, dimension(:) :: i8b_out
    integer, allocatable, dimension(:) :: int_out
    integer, allocatable, dimension(:) :: ind_part, ind_cell_part, ind_level_part

    real(dp), dimension(1:nvector, 1:ndim) :: xp_cache
    integer, dimension(1:nvector) :: ind_part_cache, ind_grid_cache, ind_cell_cache, ind_level_cache

    integer :: ind, ncache

    integer :: unit_out, unit_info, ivar

    logical :: ok, dump_info

#if NENER > 0
    integer :: irad
#endif

    if (verbose) write(*,*) 'Entering particle_snapshot'

    ! Set ivar to 1 for first variable
    ivar = 1

    npart_selected = 0
    do j = 1, npartmax
      call select_particle(j, npart_selected)
    end do

    allocate(ind_part(npart_selected))
    allocate(ind_cell_part(npart_selected))
    allocate(ind_level_part(npart_selected))
    allocate(xdp_out(npart_selected))
    allocate(i8b_out(npart_selected))
    allocate(int_out(npart_selected))

    ! Store particle indices and their host cell
    ipart = 1
    ncache = 0
    do j = 1, npartmax
      if (itmpp(j) > 0) then
        ind_part(ipart) = j

        ncache = ncache + 1
        xp_cache(ncache, :) = xp(j, :)
        ind_part_cache(ncache) = ipart
        if (ncache == nvector) then
          call find_grid_containing( &
            xp_cache,                &
            ind_grid_cache,          &
            ind_cell_cache,          &
            ind_level_cache,         &
            ncache                   &
          )
          do jcache = 1, ncache
            ind_cell_part(ind_part_cache(jcache)) = ind_cell_cache(jcache)
          end do
          ncache = 0
        end if
        ipart = ipart + 1
      end if
    end do

    ! Empty cache if necessary
    if (ncache > 0) then
      call find_grid_containing( &
        xp_cache,                &
        ind_grid_cache,          &
        ind_cell_cache,          &
        ind_level_cache,         &
        ncache                   &
      )
      do jcache = 1, ncache
        ind_cell_part(ind_part_cache(jcache)) = ind_cell_cache(jcache)
      end do
      ncache = 0
    end if

    ! If on main process
    if (myid == 1) then
      open(newunit=unit_info, file="movie1/particle_file_descriptor.txt", form="formatted")
      call dump_header_info(unit_info)
      dump_info = .true.
    else
      dump_info = .false.
    end if

    ! Now write file
    open(newunit=unit_out, file=filename, status="unknown", form="unformatted", position="append")

    ! Write headers
    write(unit_out) imov, npart_selected

    ! Write particle ids
    do ipart = 1, npart_selected
      i8b_out(ipart) = idp(ind_part(ipart))
    end do
    call generic_dump("identity", ivar, i8b_out, unit_out, dump_info, unit_info)

    ! Write particle type+tag
    do ipart = 1, npart_selected
      int_out(ipart) = part2int(typep(ind_part(ipart)))
    end do
    call generic_dump("familytag", ivar, int_out, unit_out, dump_info, unit_info)

    ! Write position
    do idim = 1, ndim
      do ipart = 1, npart_selected
        xdp_out(ipart) = xp(ind_part(ipart), idim)
      end do
      call generic_dump("position_"//dim_keys(idim), ivar, xdp_out, unit_out, dump_info, unit_info)
    end do

    ! Write gas/star velocity
    do idim = 1, ndim
      do ipart = 1, npart_selected
        if (is_gas_tracer(typep(ind_part(ipart)))) then
          ind = ind_cell_part(ipart)
          xdp_out(ipart) = uold(ind, 1 + idim) / max(uold(ind, 1), smallr)
        else if (is_tracer(typep(ind_part(ipart)))) then
          ind = partp(ind_part(ipart))
          xdp_out(ipart) = vp(ind, idim)
        else
          xdp_out(ipart) = vp(ind_part(ipart), idim)
        end if
      end do
      call generic_dump("velocity_"//dim_keys(idim), ivar, xdp_out, unit_out, dump_info, unit_info)
    end do

    ! Write *gas* density
    do ipart = 1, npart_selected
      xdp_out(ipart) = uold(ind_cell_part(ipart), 1)
    end do
    call generic_dump("density", ivar, xdp_out, unit_out, dump_info, unit_info)

    ! Write pressure
    do ipart = 1, npart_selected
      ind = ind_cell_part(ipart)
      xdp_out(ipart) = uold(ind, ndim + 2)
      xdp_out(ipart) = xdp_out(ipart) - 0.5d0 * uold(ind, 2) ** 2 / max(uold(ind, 1), smallr)
#if NDIM > 1
      xdp_out(ipart) = xdp_out(ipart) - 0.5d0 * uold(ind, 3) ** 2 / max(uold(ind, 1), smallr)
#endif
#if NDIM > 2
      xdp_out(ipart) = xdp_out(ipart) - 0.5d0 * uold(ind, 4) ** 2 / max(uold(ind, 1), smallr)
#endif
#if NENER > 0
      do irad = 1, nener
          xdp_out(ipart) = xdp_out(ipart) - uold(ind, ndim + 2 + irad)
      end do
#endif
      xdp_out(ipart) = (gamma - 1d0) * xdp_out(ipart)
    end do
    call generic_dump("pressure", ivar, xdp_out, unit_out, dump_info, unit_info)

    ! Write potential
    if (poisson) then
      do ipart = 1, npart_selected
        ind = ind_cell_part(ipart)
        xdp_out(ipart) = phi(ind)
      end do
      call generic_dump("potential", ivar, xdp_out, unit_out, dump_info, unit_info)
    end if

    close(unit_out)
    close(unit_info)

  end subroutine output_particle_snapshot

  pure function binary_search(id, ids_list, np)
    ! Uses a binary search to find `id` in `ids_list`, assuming the latter is sorted in increasing order
    integer(i8b), intent(in) :: id
    integer(i8b), dimension(np), intent(in) :: ids_list
    integer, intent(in) :: np

    integer :: binary_search

    integer(i8b) :: v
    integer :: L, M, R
    L = 1; R = np

    ! Early break if the id is not within min/max
    if ((id < ids_list(L)) .or. (id > ids_list(R))) then
      binary_search = 0
      return
    end if

    do while (L <= R)
       M = (L + R) / 2
       v = ids_list(M)
       if (v < id) then
          L = M + 1
       else if (v > id) then
          R = M - 1
       else
          binary_search = M
          return
       end if
    end do

    binary_search = 0

 end function binary_search

 recursive subroutine quicksort(a, first, last)
    integer(i8b), intent(inout), dimension(:) :: a
    integer(i8b) :: x, t
    integer, intent(in) :: first, last
    integer :: i, j

    x = a( (first+last) / 2 )
    i = first
    j = last
    do
       do while (a(i) < x)
          i=i+1
       end do
       do while (x < a(j))
          j=j-1
       end do
       if (i >= j) exit
       t = a(i);  a(i) = a(j);  a(j) = t
       i=i+1
       j=j-1
    end do
    if (first < i-1) call quicksort(a, first, i-1)
    if (j+1 < last)  call quicksort(a, j+1, last)
 end subroutine quicksort

  subroutine read_particle_ids(particle_source_file, n, ids, particle_type)
    character(len=*), intent(in) :: particle_source_file
    integer, intent(out) :: n
    integer(i8b), intent(out), allocatable, dimension(:) :: ids
    integer, intent(out) :: particle_type

    integer :: read_unit, ios
    integer(i8b) :: idpart

    open(newunit=read_unit, file=trim(particle_source_file), form='formatted', iostat=ios, status="old")
    if (ios /= 0) then
      print*, "ERROR: Could not read file `", trim(particle_source_file), "`"
      call clean_stop
    end if

    read(read_unit, *, iostat=ios) particle_type
    if (particle_type < -5 .or. particle_type > 5) then
      print*, "ERROR: the first line of", trim(particle_source_file), "should be a valid particle family. Got", particle_type
      call clean_stop
    end if
    n = 0
    do
      read(read_unit, *, iostat=ios) idpart
      if (ios /= 0) exit
      n = n + 1
    end do

    allocate(ids(1:n))
    rewind(read_unit)
    n = 0
    read(read_unit, *) particle_type ! Pass first line
    do
      read(read_unit, *, iostat=ios) idpart
      if (ios /= 0) exit
      n = n + 1
      ids(n) = idpart
    end do

    close(read_unit)

    call quicksort(ids, 1, n)

  end subroutine read_particle_ids


end module particle_snapshot
