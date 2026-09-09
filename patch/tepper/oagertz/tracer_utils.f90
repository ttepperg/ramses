! Hooks for RAMSES
module tracer_utils
   use amr_parameters      ! nx, ny, nz, dp, ngridmax, nvector, …
   use amr_commons         ! ncoarse, father, xg, son, myid, cpu_map, cpu_map2, ncpu, nbor
   use pm_commons          ! xp, tp, idp, levelp, headp, mp, localseed, numbp, nextp, move_flag
   use random, only        : ranf
   use hydro_commons, only : uold, if1, if2, jf1, jf2, kf1, kf2, nvar
#ifndef WITHOUTMPI
   use mpi_mod
#endif
   implicit none

   real(kind=dp), dimension(1:3) :: skip_loc
   real(kind=dp) :: scale
   integer :: nx_loc

   logical :: ddebug = .false.
   ! logical :: ddebug = .true.
   integer :: dbpart = -1 !105

   real(dp), allocatable, dimension(:) :: proba_yield  ! Working array

   contains

   subroutine initialize_skip_loc
      logical, save :: firstCall = .true.

      if (firstCall) then
         skip_loc=(/0.0d0, 0.0d0, 0.0d0/)
         if(ndim>0) skip_loc(1) = dble(icoarse_min)
         if(ndim>1) skip_loc(2) = dble(jcoarse_min)
         if(ndim>2) skip_loc(3) = dble(kcoarse_min)

         nx_loc=(icoarse_max-icoarse_min+1)
         scale = boxlen/dble(nx_loc)

         firstCall = .false.
      end if

   end subroutine initialize_skip_loc

   !---------------------------------------------------------------------------------
   ! Hooks before and after grid creation
   !---------------------------------------------------------------------------------
   subroutine pre_kill_grid_hook(ind_cell, ilevel, nn, ibound, boundary_region)

      integer, intent(in) :: nn, ilevel, ibound
      logical, intent(in) :: boundary_region
      integer, intent(in), dimension(1:nvector) :: ind_cell

      !######################
      ! Customize here
      !######################

      integer :: ipart, igrid, i, j, dim
      real(dp) :: dx, prevxp(1:ndim)

      real(dp), dimension(1:ndim) :: tmp_xp

      call initialize_skip_loc

      ! print*, 'in pre_kill_grid_hook'
      dx = 0.5D0**ilevel

      if (MC_tracer) then
         ! For all particles, recenter them in the center of the grid that's being deleted
         ! (becoming a cell)
         do j = 1, nn
            igrid = son(ind_cell(j))
            ipart = headp(igrid)

            do i = 1, numbp(igrid)
               if (is_gas_tracer(typep(ipart))) then
                  prevxp(:) = xp(ipart, 1:ndim)
                  do dim = 1, ndim
                     xp(ipart, dim) = (xg(igrid, dim) - skip_loc(dim)) * scale
                  end do
                  ! Attach to the new cell
                  partp(ipart) = ind_cell(j)
                  tmp_xp(:) = xp(ipart, :)
                  call safe_move(tmp_xp, prevxp(1:ndim))
                  xp(ipart, :) = tmp_xp(:)
               end if
               ipart = nextp(ipart)
            end do
         end do
      end if
   end subroutine pre_kill_grid_hook

   subroutine post_kill_grid_hook(ind_cell, ilevel, nn, ibound, boundary_region)
      use amr_commons

      integer, intent(in) :: nn, ilevel, ibound
      logical, intent(in) :: boundary_region
      integer, intent(in), dimension(:) :: ind_cell

      !######################
      ! Customize here
      !######################

      call initialize_skip_loc

   end subroutine post_kill_grid_hook

   !---------------------------------------------------------------------------------
   ! Hooks before and after grid creation
   !---------------------------------------------------------------------------------
   subroutine pre_make_grid_fine_hook(ind_grid, ind_cell, ind, &
      ilevel, nn, ibound, boundary_region)
      use amr_commons
      integer, intent(in) :: nn, ind, ilevel, ibound
      logical, intent(in) :: boundary_region
      integer, dimension(1:nvector), intent(in) :: ind_grid, ind_cell

      !######################
      ! Customize here
      !######################
      call initialize_skip_loc

   end subroutine pre_make_grid_fine_hook

   subroutine post_make_grid_fine_hook(ind_grid, ind_cell, ind, &
      ilevel, nn, ibound, boundary_region)
      use amr_commons
      use hydro_commons
      use pm_commons
      use random
      integer, intent(in) :: nn, ind, ilevel, ibound
      logical, intent(in) :: boundary_region
      integer, dimension(1:nvector), intent(in) :: ind_grid, ind_cell

      !######################
      ! Customize here
      !######################
      real(dp) :: dx, dxcoarse

      real(dp) :: x(1:ndim), rand, mass(0:twotondim), prevxp(1:ndim)
      integer :: j, i, fgrid, igrid, ipart, ison, iskip, icell
      integer :: loc(1:3)
      logical :: ok

      real(dp), dimension(1:ndim) :: tmp_xp

      call initialize_skip_loc

      dx = 0.5D0**ilevel           ! dx of the new level
      dxcoarse = 0.5D0**(ilevel-1) ! dx of the previous level

      ! print*, 'in post_make_grid_fine_hook'
      if (MC_tracer) then
         ! Compute the expected location of particles relative to xg in dx units
         loc(3) = (ind-1) / 4
         loc(2) = (ind-1-loc(3)*4) / 2
         loc(1) = (ind-1-loc(3)*4-loc(2)*2)

         do j = 1, nn
            fgrid = ind_grid(j)
            igrid = son(ind_cell(j))
            ipart = headp(fgrid)

            ! Load masses
            do ison = 1, twotondim
               iskip = ncoarse + (ison-1)*ngridmax
               icell = iskip + igrid
               mass(ison) = uold(icell, 1)
            end do

            mass(0) = sum(mass(1:twotondim))

            do i = 1, numbp(fgrid)
               if (is_gas_tracer(typep(ipart))) then

                  ! Check whether the particle was in the refined cell
                  x(1:ndim) = cellCenter(ind, fgrid, dxcoarse)

                  ok = all(xp(ipart, 1:ndim) == x(1:ndim))
                  ! If the particle is in refined cell, spread it accordingly
                  if (ok) then

                     ! Pick a random direction
                     call ranf(tracer_seed, rand)

                     do ison = 1, twotondim
                        if (rand < mass(ison) / mass(0)) then
                           ! Move particle to center of new cells
                           prevxp(:) = xp(ipart, :)
                           xp(ipart, :) = cellCenter(ison, igrid, dx)

                           tmp_xp(:) = xp(ipart, :)
                           call safe_move(tmp_xp, prevxp(1:ndim))
                           xp(ipart, :) = tmp_xp(:)
                           iskip = ncoarse + (ison-1)*ngridmax
                           partp(ipart) = igrid + iskip
                           exit
                        else
                           rand = rand - mass(ison) / mass(0)
                        end if
                     end do
                  end if
               end if
               ipart = nextp(ipart)
            end do
         end do
      end if
   end subroutine post_make_grid_fine_hook

   ! Local function that returns the *cell* center given a grid and the
   ! index of the cell within the grid.
   function cellCenter(ison, ind_grid, dx) result (x)
      implicit none

      integer, intent(in)  :: ison, ind_grid
      real(dp), intent(in) :: dx
      real(dp), dimension(1:ndim) :: x

      real(dp), dimension(1:3) :: xc

      integer :: ixn, iyn, izn, dim
      ! Get the location of neighbor cell in its grid
      izn = (ison-1)/4
      iyn = (ison-1-4*izn)/2
      ixn = (ison-1-4*izn-2*iyn)

      ! Compute the expected location of the particles
      xc(1) = (dble(ixn)-0.5d0)*dx
      xc(2) = (dble(iyn)-0.5d0)*dx
      xc(3) = (dble(izn)-0.5d0)*dx

      do dim = 1, ndim
         x(dim) = (xg(ind_grid, dim) + xc(dim) - skip_loc(dim)) * scale
      end do
   end function cellCenter

   ! Move particles taking into account periodic boundary conditions
   subroutine safe_move(new_xp, prev_xp)
      real(dp), dimension(1:ndim), intent(inout) :: new_xp
      real(dp), dimension(1:ndim), intent(in)    :: prev_xp

      integer :: dim

      do dim = 1, ndim
         if (new_xp(dim) - prev_xp(dim) > 0.5d0*scale) then
            new_xp(dim) = new_xp(dim) - scale
         else if (new_xp(dim) - prev_xp(dim) < -0.5d0*scale) then
            new_xp(dim) = new_xp(dim) + scale
         end if
      end do

   end subroutine safe_move

   ! Helper function to determine the relative level of neighbouring cells
   subroutine relative_level(ind_ngrid, ind_ncell, direction, level)
      integer, intent(in) :: direction
      integer, intent(in), dimension(0:twondim) :: ind_ncell ! the neighboring cells
      integer, intent(in), dimension(0:twondim) :: ind_ngrid ! the neighboring parent grids
      integer, intent(out) :: level

      integer :: pos, ind_ngrid_ncell

      ! If the neighbor cell is a grid
      if (son(ind_ncell(direction)) > 0 ) then
         level = +1
      else
         ! get the grid containing the neighbor cell
         pos = (ind_ncell(direction)-ncoarse-1)/ngridmax+1
         ind_ngrid_ncell = ind_ncell(direction)-ncoarse-(pos-1)*ngridmax

         ! if the father grid and the neighbor cell's grid are neighors / the same
         if (ind_ngrid_ncell == ind_ngrid(direction) .or. ind_ngrid_ncell == ind_ngrid(0)) then
            level = 0
         else
            level = -1
         end if

      end if

   end subroutine relative_level

   ! Compute all four cells (in 3D) on the cell pointing in
   ! the given direction
   subroutine get_cells_on_face(direction, locs)
      implicit none

      integer, intent(in) :: direction
      integer, dimension(1:twotondim/2), intent(out) :: locs

      integer, save, dimension(1:6, 1:4) :: mapping
      logical, save :: firstCall = .true.

      integer, save :: twotondimo2 = twotondim/2

      if (firstCall) then
         mapping(1, 1:4) = (/1, 3, 5, 7/) ! left cells
         mapping(2, 1:4) = (/2, 4, 6, 8/) ! right cells
         mapping(3, 1:4) = (/1, 2, 5, 6/) ! top cells
         mapping(4, 1:4) = (/3, 4, 7, 8/) ! bottom cells
         mapping(5, 1:4) = (/1, 2, 3, 4/) ! front cells
         mapping(6, 1:4) = (/5, 6, 7, 8/) ! back cells

         firstCall = .false.
      end if

      locs(1:twotondimo2) = mapping(direction, 1:twotondimo2)

   end subroutine get_cells_on_face

   !----------------------------------------------------
   ! Handle tracers being released from other particles
   ! This is used e.g. when doing feedback

   ! This should be called before doing particle yield
   subroutine pre_particle_yield()
      if (.not. allocated(proba_yield)) then
         allocate(proba_yield(npartmax))
      end if
      proba_yield = 0
   end subroutine

   ! This should be called whenever a particle yields mass (and hence tracers)
   subroutine mark_yielding_particle(ipart, yield_proba)
      integer, intent(in) :: ipart ! Index of the yielding particle
      real(dp), intent(in) :: yield_proba ! Probability for a tracer to be yielded

      proba_yield(ipart) = yield_proba
   end subroutine

   ! This detaches the tracer from the particle
   subroutine yield_tracers(icpu, ilevel, type_to_yield)
      integer, intent(in)      :: icpu, ilevel
      type(part_t), intent(in) :: type_to_yield

      integer  :: igrid, jgrid, ipart, jpart, npart1
      real(dp) :: rand

      igrid = headl(icpu,ilevel)

      ! Loop over grids
      do jgrid=1,numbl(icpu,ilevel)
         npart1=numbp(igrid)  ! Number of particles in the grid

         ! Loop over tracer particles
         ipart = headp(igrid)
         do jpart = 1, npart1
            if (typep(ipart)%family == type_to_yield%family) then
               call ranf(tracer_seed, rand)

               ! Detach particles
               if (rand < proba_yield(partp(ipart))) then
                  typep(ipart)%family = FAM_TRACER_GAS
                  ! Change here to tag particles that were on a star
                  move_flag(ipart) = 1
               end if
            end if
            ipart = nextp(ipart)
         end do
         ! End loop over tracer

         igrid = next(igrid) ! Go to next grid
      end do
   end subroutine

   subroutine yield_tracers_within_radius(icpu, ilevel, radius, type_to_yield)
      use constants, only: pi
      integer, intent(in)      :: icpu, ilevel
      real(dp)                 :: radius
      type(part_t), intent(in) :: type_to_yield

      integer  :: igrid, jgrid, ipart, jpart, npart1, next_part
      real(dp) :: rand1, rand2, rand3, rand_r, rand_theta, rand_phi

      real(dp), dimension(1:ndim) :: dxp

      igrid = headl(icpu,ilevel)

      do jgrid = 1, numbl(icpu,ilevel)
         ! Loop over particles
         npart1 = numbp(igrid)  ! Number of particles in the grid

         ipart=headp(igrid)
         ! Loop over particles
         do jpart=1,npart1
            next_part = nextp(ipart)
            if (typep(ipart)%family == type_to_yield%family) then
               call ranf(tracer_seed, rand1)

               ! Detach particles
               if (rand1 < proba_yield(partp(ipart))) then
                  move_flag(ipart) = 1
                  ! Generate a random position within the explosion radius.
                  ! See https://math.stackexchange.com/questions/87230/picking-random-points-in-the-volume-of-sphere-with-uniform-probability
                  call ranf(tracer_seed, rand1)
                  call ranf(tracer_seed, rand2)
                  call ranf(tracer_seed, rand3)
                  rand_r = rand1**(1._dp/3._dp) * radius
                  rand_theta = acos(2*rand2 - 1)
                  rand_phi = 2*pi * rand3
                  dxp(1) = rand_r * sin(rand_theta) * cos(rand_phi)
#if NDIM > 1
                  dxp(2) = rand_r * sin(rand_theta) * sin(rand_phi)
#if NDIM > 2
                  dxp(3) = rand_r * cos(rand_theta)
#endif
#endif
                  xp(ipart, :) = xp(ipart, :) + dxp(:)
                  typep(ipart)%family = FAM_TRACER_GAS
                  typep(ipart)%tag = typep(ipart)%tag + 1_1
                  ! Change here to tag particles that were on a star
               end if
            end if
            ipart = next_part  ! Go to next particle
            ! End loop over particles
         end do
         igrid = next(igrid) ! Go to next grid
      end do ! End loop over cached grids
   end subroutine

   ! This should be called after doing particle yield
   subroutine post_particle_yield()
      proba_yield = 0
   end subroutine

   ! This routines decrease by one the move_flag of the MC tracer at
   ! level ilevel
   subroutine reset_tracer_move_flag(ilevel)
      use pm_commons
      use amr_commons

      implicit none

      integer, intent(in) :: ilevel

      integer :: ipart, jpart, next_part, jgrid, npart1, igrid

      ! Loop over grids
      igrid = headl(myid, ilevel)
      do jgrid = 1, numbl(myid, ilevel)
         npart1 = numbp(igrid)  ! Number of particles in the grid
         if (npart1 > 0) then
            ipart = headp(igrid)
            ! Loop over particles
            do jpart = 1, npart1
               ! Save next particle  <---- Very important !!!
               next_part = nextp(ipart)

               if (is_tracer(typep(ipart))) then
                  move_flag(ipart) = max(move_flag(ipart) - 1, 0)
               end if
               ipart = next_part  ! Go to next particle
            end do
         end if
         igrid = next(igrid)   ! Go to next grid
      end do

    end subroutine reset_tracer_move_flag


   ! Attach tracer particles to other particles (e.g. stars)
   subroutine attach_tracer(ind_tracer, proba, x_target_part, ind_target_part, nattach)
      use amr_commons
      use random
      use pm_commons
      implicit none

      integer, intent(in) :: nattach
      integer, dimension(1:nvector), intent(in) :: ind_tracer, ind_target_part
      real(dp), dimension(1:nvector), intent(in) :: proba
      real(dp), dimension(1:nvector, 1:3), intent(in) :: x_target_part

      logical, dimension(1:nvector), save :: attach = .false.
      integer :: i, idim
      real(dp) :: r

      do i = 1, nattach
         call ranf(tracer_seed, r)
         attach(i) = r < proba(i)
      end do

      ! Change particle pointer and kind
      do i = 1, nattach
         if (attach(i)) then
            ! Tag particle with target particle index
            partp(ind_tracer(i)) = ind_target_part(i)

            ! The tracer family is the opposite of the family of the target particle
            typep(ind_tracer(i))%family = -typep(ind_target_part(i))%family
         end if
      end do

      ! Change particle location
      do idim = 1, ndim
         do i = 1, nattach
            if (attach(i)) then
               xp(ind_tracer(i), idim) = x_target_part(i, idim)
            end if
         end do
      end do

      ! Save state of the particle
      do i = 1, nattach
         if (attach(i)) then
            ! Set move_flag to 1 to prevent further move
            move_flag(ind_tracer(i)) = 1
         end if
      end do

    end subroutine attach_tracer

    subroutine find_grid_containing_max(pos, ind_grid, ind_cell, ind_level, npart, maxlevel)
      ! Find the finest grid and cell containing the positions.
      use amr_commons, only : levelmin, ncoarse, son, xg
      use amr_parameters, only : ndim, ngridmax, nvector
      real(dp), intent(in), dimension(1:nvector, 1:ndim) :: pos
      integer, intent(in) :: npart, maxlevel

      integer, intent(out), dimension(1:nvector) :: ind_grid, ind_cell, ind_level

      integer :: igrid, ipart, idim, ilevel, istart, ind

      real(dp) :: xgrid(3)

      call initialize_skip_loc()

      ind_grid(1:npart) = 0; ind_cell(1:npart) = 0; ind_level(1:npart) = 0

      ! Find coarse grid containing position
      istart = 1
      igrid = 1

      ! Initialize at level = 1
      xgrid = (xg(igrid, :) - skip_loc) * scale

      if (any(pos(:, :) < 0._dp)) then
         print*, 'Got a position < 0'
         stop
      else if (any(pos(:, :) > boxlen)) then ! assuming x_box == [1, 1, 1]
         print*, 'Got a position > 1'
         stop
      end if

      do ipart = 1, npart
         ind = 1
         ! Loop over cells to compute cell containing object
         do idim = 1, ndim
            if (pos(ipart, idim) > xgrid(idim)) then
               ind = ind + 2**(idim-1)
            end if
         end do

         ! Compute cell and grid indexes
         ind_grid(ipart)  = igrid
         ind_cell(ipart)  = igrid + ncoarse + (ind - 1) * ngridmax
         ind_level(ipart) = levelmin
      end do

      ! Loop over levels
      do ipart = 1, npart
         level_loop: do ilevel = istart+1, maxlevel+1
            igrid = son(ind_cell(ipart))

            ! True if cell is refined → go down the tree
            if (igrid > 0) then
               xgrid = (xg(igrid, :) - skip_loc) * scale
               ind = 1

               ! Loop over cells to compute cell containing object
               do idim = 1, ndim
                  if (pos(ipart, idim) > xgrid(idim))then
                     ind = ind + 2**(idim-1)
                  end if
               end do

               ! Compute cell and grid indexes
               ind_grid(ipart)  = igrid
               ind_cell(ipart)  = igrid + ncoarse + (ind - 1) * ngridmax
               ind_level(ipart) = ilevel
            else
               exit level_loop
            end if
         end do level_loop
      end do ! particle loop

    end subroutine find_grid_containing_max

    subroutine find_grid_containing(pos, ind_grid, ind_cell, ind_level, npart)
      use amr_commons, only: nlevelmax
      ! Find the finest grid and cell containing the positions.
      real(dp), intent(in), dimension(1:nvector, 1:ndim) :: pos
      integer, intent(in) :: npart

      integer, intent(out), dimension(1:nvector) :: ind_grid, ind_cell, ind_level
      call find_grid_containing_max(pos, ind_grid, ind_cell, ind_level, npart, nlevelmax)
    end subroutine find_grid_containing

end module tracer_utils
