subroutine load_dice_restart
!!! DICE
  dice_init=.true.
  ! Conversion factor from user units to cgs units
  call units(scale_l,scale_t,scale_d,scale_v,scale_nH,scale_T2)
  scale_m = scale_d*scale_l**3
  ! Reading header of the Gadget file
  error=.false.
  ipart    = 0
  do ifile=1,ic_nfile
	 write(ifile_str,*) ifile
	 if(ic_nfile.eq.1) then
		filename=TRIM(initfile(levelmin))//'/'//TRIM(ic_file)
	 else
		filename=TRIM(initfile(levelmin))//'/'//TRIM(ic_file)//'.'//ADJUSTL(ifile_str)
	 endif
	 INQUIRE(FILE=filename,EXIST=file_exists)
	 if(.not.file_exists) then
		if(myid==1) write(*,*) TRIM(filename)," not found"
		call clean_stop
	 endif
	 if(myid==1)then
		write(*,'(A12,A)') " Opening -> ",filename
		if((ic_format.ne.'Gadget1').and.(ic_format.ne.'Gadget2')) then
		   if(myid==1) write(*,*) 'Specify a valid IC file format [ic_format=Gadget1/Gadget2]'
		   error=.true.
		endif
		OPEN(unit=1,file=filename,status='old',action='read',form='unformatted',access="stream")
		! Init block address
		head_blck  = -1
		pos_blck   = -1
		vel_blck   = -1
		id_blck    = -1
		u_blck     = -1
		mass_blck  = -1
		metal_blck = -1
		age_blck   = -1

		if(ic_format .eq. 'Gadget1') then
		   ! Init block counter
		   jump_blck = 1
		   blck_cnt = 1
		   do while(.true.)
			  ! Reading data block header
			  read(1,POS=jump_blck,iostat=stat) blck_size
			  if(stat /= 0) exit
			  ! Saving data block positions
			  if(blck_cnt .eq. 1) then
				 head_blck  = jump_blck+sizeof(blck_size)
				 head_size  = blck_size
				 write(*,*)blck_cnt,blck_size
			  endif
			  if(blck_cnt .eq. 2) then
				 pos_blck   = jump_blck+sizeof(blck_size)
				 pos_size   = blck_size/(3*sizeof(dummy_real))
				 write(*,*)blck_cnt,blck_size
			  endif
			  if(blck_cnt .eq. 3) then
				 vel_blck   = jump_blck+sizeof(blck_size)
				 vel_size   = blck_size/(3*sizeof(dummy_real))
				 write(*,*)blck_cnt,blck_size
			  endif
			  if(blck_cnt .eq. 4) then
				 id_blck    = jump_blck+sizeof(blck_size)
				 id_size    = blck_size/sizeof(dummy_int)
				 write(*,*)blck_cnt,blck_size
			  endif
			  if(blck_cnt .eq. 5) then
				 u_blck     = jump_blck+sizeof(blck_size)
				 u_size     = blck_size/sizeof(dummy_real)
				 write(*,*)blck_cnt,blck_size
			  endif
			  if(blck_cnt .eq. 6) then
				 mass_blck  = jump_blck+sizeof(blck_size)
				 mass_size  = blck_size/sizeof(dummy_real)
				 write(*,*)blck_cnt,blck_size
			  endif
			  if(blck_cnt .eq. 7) then
				 metal_blck = jump_blck+sizeof(blck_size)
				 metal_size = blck_size/sizeof(dummy_real)
				 write(*,*)blck_cnt,blck_size
			  endif
			  if(blck_cnt .eq. 8) then
				 age_blck   = jump_blck+sizeof(blck_size)
				 age_size   = blck_size/sizeof(dummy_real)
				 write(*,*)blck_cnt,blck_size
			  endif
			  jump_blck = jump_blck+blck_size+2*sizeof(dummy_int)
			  blck_cnt = blck_cnt+1
		   enddo
		endif

		if(ic_format .eq. 'Gadget2') then
		   ! Init block counter
		   jump_blck = 1
		   write(*,'(A50)')"__________________________________________________"
		   do while(.true.)
			  ! Reading data block header
			  read(1,POS=jump_blck,iostat=stat) dummy_int
			  if(stat /= 0) exit
			  read(1,POS=jump_blck+sizeof(dummy_int),iostat=stat) blck_name
			  if(stat /= 0) exit
			  read(1,POS=jump_blck+sizeof(dummy_int)+sizeof(blck_name),iostat=stat) dummy_int
			  if(stat /= 0) exit
			  read(1,POS=jump_blck+2*sizeof(dummy_int)+sizeof(blck_name),iostat=stat) dummy_int
			  if(stat /= 0) exit
			  read(1,POS=jump_blck+3*sizeof(dummy_int)+sizeof(blck_name),iostat=stat) blck_size
			  if(stat /= 0) exit
			  ! Saving data block positions
			  if(blck_name .eq. ic_head_name) then
				 head_blck  = jump_blck+sizeof(blck_name)+4*sizeof(dummy_int)
				 head_size  = blck_size
				 write(*,*) '-> Found ',blck_name,' block'
			  endif
			  if(blck_name .eq. ic_pos_name) then
				 pos_blck   = jump_blck+sizeof(blck_name)+4*sizeof(dummy_int)
				 pos_size   = blck_size/(3*sizeof(dummy_real))
				 write(*,*) '-> Found ',blck_name,' block'
			  endif
			  if(blck_name .eq. ic_vel_name) then
				 vel_blck  = jump_blck+sizeof(blck_name)+4*sizeof(dummy_int)
				 vel_size  = blck_size/(3*sizeof(dummy_real))
				 write(*,*) '-> Found ',blck_name,' block'
			  endif
			  if(blck_name .eq. ic_id_name) then
				 id_blck    = jump_blck+sizeof(blck_name)+4*sizeof(dummy_int)
				 id_size    = blck_size/sizeof(dummy_int)
				 write(*,*) '-> Found ',blck_name,' block'
			  endif
			  if(blck_name .eq. ic_mass_name) then
				 mass_blck  = jump_blck+sizeof(blck_name)+4*sizeof(dummy_int)
				 mass_size  = blck_size/sizeof(dummy_real)
				 write(*,*) '-> Found ',blck_name,' block'
			  endif
			  if(blck_name .eq. ic_u_name) then
				 u_blck     = jump_blck+sizeof(blck_name)+4*sizeof(dummy_int)
				 u_size     = blck_size/sizeof(dummy_real)
				 write(*,*) '-> Found ',blck_name,' block'
			  endif
			  if(blck_name .eq. ic_metal_name) then
				 metal_blck = jump_blck+sizeof(blck_name)+4*sizeof(dummy_int)
				 metal_size = blck_size/sizeof(dummy_real)
				 write(*,*) '-> Found ',blck_name,' block'
			  endif
			  if(blck_name .eq. ic_age_name) then
				 age_blck   = jump_blck+sizeof(blck_name)+4*sizeof(dummy_int)
				 age_size   = blck_size/sizeof(dummy_real)
				 write(*,*) '-> Found ',blck_name,' block'
			  endif
			  jump_blck = jump_blck+blck_size+sizeof(blck_name)+5*sizeof(dummy_int)
		   enddo
		endif

		if((head_blck.eq.-1).or.(pos_blck.eq.-1).or.(vel_blck.eq.-1)) then
		   write(*,*) 'Gadget file does not contain handful data'
		   error=.true.
		endif
		if(head_size.ne.256) then
		   write(*,*) 'Gadget header is not 256 bytes'
		   error=.true.
		endif

		! Byte swapping doesn't appear to work if you just do READ(1)header
		READ(1,POS=head_blck) header%npart,header%mass,header%time,header%redshift, &
			 header%flag_sfr,header%flag_feedback,header%nparttotal, &
			 header%flag_cooling,header%numfiles,header%boxsize, &
			 header%omega0,header%omegalambda,header%hubbleparam, &
			 header%flag_stellarage,header%flag_metals,header%totalhighword, &
			 header%flag_entropy_instead_u, header%flag_doubleprecision, &
			 header%flag_ic_info, header%lpt_scalingfactor

		nstar_tot = sum(header%npart(3:5))
		npart     = sum(header%npart)
		ngas      = header%npart(1)
		nhalo     = header%npart(2)
		if(cosmo) T2_start = 1.356d-2/aexp**2

		write(*,'(A50)')"__________________________________________________"
		write(*,*)"Found ",npart," particles"
		skip=.false.
		do j=1,6
		   if(ic_skip_type(j).eq.0) skip=.true.
		enddo
		if(.not.skip) write(*,*)"----> ",header%npart(1)," type 0 particles with header mass ",header%mass(1)
		skip=.false.
		do j=1,6
		   if(ic_skip_type(j).eq.1) skip=.true.
		enddo
		if(.not.skip) write(*,*)"----> ",header%npart(2)," type 1 particles with header mass ",header%mass(2)
		skip=.false.
		do j=1,6
		   if(ic_skip_type(j).eq.2) skip=.true.
		enddo
		if(.not.skip) write(*,*)"----> ",header%npart(3)," type 2 particles with header mass ",header%mass(3)
		skip=.false.
		do j=1,6
		   if(ic_skip_type(j).eq.3) skip=.true.
		enddo
		if(.not.skip) write(*,*)"----> ",header%npart(4)," type 3 particles with header mass ",header%mass(4)
		skip=.false.
		do j=1,6
		   if(ic_skip_type(j).eq.4) skip=.true.
		enddo
		if(.not.skip) write(*,*)"----> ",header%npart(5)," type 4 particles with header mass ",header%mass(5)
		skip=.false.
		do j=1,6
		   if(ic_skip_type(j).eq.5) skip=.true.
		enddo
		if(.not.skip) write(*,*)"----> ",header%npart(6)," type 5 particles with header mass ",header%mass(6)

		write(*,'(A50)')"_____________________progress_____________________"
		if((pos_size.ne.npart).or.(vel_size.ne.npart)) then
		   write(*,*) 'POS =',pos_size
		   write(*,*) 'VEL =',vel_size
		   write(*,*) 'Number of particles does not correspond to block sizes'
		   error=.true.
		endif

	 endif
	 if(error) call clean_stop
#ifndef WITHOUTMPI
	 call MPI_BCAST(nstar_tot,1,MPI_INTEGER,0,MPI_COMM_WORLD,info)
#endif
	 eob      = .false.
	 kpart    = 0
	 lpart    = 0
	 mpart    = 0
	 gpart    = 0
	 opart    = 0
	 mgas_tot = 0.
	 ipbar    = 0.
	 do while(.not.eob)
		xx=0.
		vv=0.
		ii=0
		mm=0.
		tt=0.
		zz=0.
		uu=0.
		if(myid==1)then
		   jpart=0
		   do i=1,nvector
			  jpart=jpart+1

			  ! All particles counter
			  kpart=kpart+1
			  if(kpart.le.header%npart(1)) type_index = 1
			  do j=1,5
				 if(kpart.gt.sum(header%npart(1:j)).and.kpart.le.sum(header%npart(1:j+1))) type_index = j+1
			  enddo
			  if((sum(header%npart(3:5)).gt.0).and.(kpart.gt.(header%npart(1)+header%npart(2)))) mpart=mpart+1
			  if(type_index.ne.2) gpart=gpart+1

			  ! Reading Gadget1 or Gadget2 file line-by-line
			  ! Mandatory data
			  read(1,POS=pos_blck+3*sizeof(dummy_real)*(kpart-1)) xx_sp(i,1:3)
			  read(1,POS=vel_blck+3*sizeof(dummy_real)*(kpart-1)) vv_sp(i,1:3)
			  if(header%mass(type_index).gt.0) then
				 mm_sp(i) = header%mass(type_index)
			  else
				 opart=opart+1
				 read(1,POS=mass_blck+sizeof(dummy_real)*(opart-1)) mm_sp(i)
			  endif
			  ! Optional data
			  if(id_blck.ne.-1) then
				 read(1,POS=id_blck+sizeof(dummy_int)*(kpart-1)) ii(i)
			  else
				 ii(i) = kpart
			  endif
			  if(kpart.le.header%npart(1)) then
				 if((u_blck.ne.-1).and.(u_size.eq.header%npart(1))) then
					read(1,POS=u_blck+sizeof(dummy_real)*(kpart-1)) uu_sp(i)
				 endif
			  endif
			  if(metal) then
				 if((metal_blck.ne.-1).and.(metal_size.eq.npart)) then
					read(1,POS=metal_blck+sizeof(dummy_real)*(kpart-1)) zz_sp(i)
				 endif
				 if((metal_blck.ne.-1).and.(metal_size.eq.ngas+nstar_tot)) then
					read(1,POS=metal_blck+sizeof(dummy_real)*(gpart-1)) zz_sp(i)
				 endif
			  endif
			  if(star) then
				 if((age_blck.ne.-1).and.(age_size.eq.sum(header%npart(3:5)))) then
					if((sum(header%npart(3:5)).gt.0).and.(kpart.gt.(header%npart(1)+header%npart(2)))) then
					   read(1,POS=age_blck+sizeof(dummy_real)*(mpart-1)) tt_sp(i)
					endif
				 endif
			  endif
			  ! Scaling to ramses code units
			  if(cosmo) then
				 gadget_scale_l = scale_l/header%boxsize
				 gadget_scale_v = 1e3*SQRT(aexp)/header%boxsize*aexp/100.
			  endif
			  xx(i,:)   = xx_sp(i,:)*(gadget_scale_l/scale_l)*ic_scale_pos
			  vv(i,:)   = vv_sp(i,:)*(gadget_scale_v/scale_v)*ic_scale_vel
			  mm(i)     = mm_sp(i)*(gadget_scale_m/scale_m)*ic_scale_mass
			  if(cosmo) then
				 if(type_index .eq. 1) mass_sph = mm(i)
				 if(xx(i,1)<  0.0d0  )xx(i,1)=xx(i,1)+dble(nx)
				 if(xx(i,1)>=dble(nx))xx(i,1)=xx(i,1)-dble(nx)
				 if(xx(i,2)<  0.0d0  )xx(i,2)=xx(i,2)+dble(ny)
				 if(xx(i,2)>=dble(ny))xx(i,2)=xx(i,2)-dble(ny)
				 if(xx(i,3)<  0.0d0  )xx(i,3)=xx(i,3)+dble(nz)
				 if(xx(i,3)>=dble(nz))xx(i,3)=xx(i,3)-dble(nz)
			  endif

			  if(metal) then
				 if(metal_blck.ne.-1) then
					zz(i) = zz_sp(i)*ic_scale_metal
				 else
					zz(i) = 0.02*z_ave
				 endif
			  endif
			  if(kpart.gt.header%npart(1)+header%npart(2)) then
				 if(age_blck.ne.-1) then
					if(cosmo) then
					   tt(i) = tt_sp(i)
					else
					   tt(i) = tt_sp(i)*(gadget_scale_t/(scale_t/aexp**2))*ic_scale_age
					endif
				 else
					tt(i) = -13.8*1d9*3.15360d7/scale_t ! Age of the universe
				 endif
			  endif
			  if(kpart.le.header%npart(1)) then
				 if(cosmo) then
					uu(i) = T2_start/scale_T2
				 else
					! Temperature stored in units of K/mu
					uu(i) = uu_sp(i)*mu_mol*(gadget_scale_v/scale_v)**2*ic_scale_u
				 endif

			  endif
			  if(kpart.le.header%npart(1)) mgas_tot = mgas_tot+mm(i)
			  ! Check the End Of Block
			  if(kpart.ge.ipbar*(npart/49.0))then
				 write(*,'(A1)',advance='no') "_"
				 ipbar = ipbar+1.0
			  endif
			  if(kpart.ge.npart) then
				 write(*,'(A1)') " "
				 write(*,'(A,A7,A)') ' ',TRIM(ic_format),' file successfully loaded'
				 write(*,'(A50)')"__________________________________________________"
				 eob=.true.
				 exit
			  endif
		   enddo
		endif
#ifndef WITHOUTMPI
		call MPI_BCAST(eob,1         ,MPI_LOGICAL         ,0,MPI_COMM_WORLD,info)
		call MPI_BCAST(xx,nvector*3  ,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,info)
		call MPI_BCAST(vv,nvector*3  ,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,info)
		call MPI_BCAST(ii,nvector    ,MPI_INTEGER         ,0,MPI_COMM_WORLD,info)
		call MPI_BCAST(mm,nvector    ,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,info)
		call MPI_BCAST(zz,nvector    ,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,info)
		call MPI_BCAST(tt,nvector    ,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,info)
		call MPI_BCAST(uu,nvector    ,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,info)
		call MPI_BCAST(jpart,1       ,MPI_INTEGER         ,0,MPI_COMM_WORLD,info)
		call MPI_BCAST(header%npart,6,MPI_INTEGER         ,0,MPI_COMM_WORLD,info)
		call cmp_cpumap(xx,cc,jpart)
#endif
		do i=1,jpart
#ifndef WITHOUTMPI
		   ! Check the CPU map
		   if(cc(i)==myid)then
#endif
			  ! Determine current particle type
			  if((lpart+i).le.header%npart(1)) type_index = 1
			  do j=1,5
				 if((lpart+i).gt.sum(header%npart(1:j)).and.(lpart+i).le.sum(header%npart(1:j+1))) type_index = j+1
			  enddo
			  skip           = .false.
			  do j=1,6
				 if(ic_skip_type(j).eq.type_index-1) skip=.true.
			  enddo
			  if(.not.skip) then
				 if(abs(xx(i,1)-ic_center(1)).ge.boxlen/2d0) cycle
				 if(abs(xx(i,2)-ic_center(2)).ge.boxlen/2d0) cycle
				 if(abs(xx(i,3)-ic_center(3)).ge.boxlen/2d0) cycle
				 ipart          = ipart+1
				 if(ipart.gt.npartmax) then
					write(*,*) "Increase npartmax"
#ifndef WITHOUTMPI
					call MPI_ABORT(MPI_COMM_WORLD,1,info)
#else
					stop
#endif
				 endif
				 xp(ipart,1:3)  = xx(i,1:3)+boxlen/2.0D0-ic_center(1:3)
				 vp(ipart,1:3)  = vv(i,1:3)
				 ! Flag gas particles with idp=1
				 if(type_index.gt.1)then
					idp(ipart)   = ii(i)+1
				 else
					idp(ipart)   = 1
				 endif
				 mp(ipart)      = mm(i)
				 levelp(ipart)  = levelmin
				 if(star) then
					tp(ipart)    = tt(i)
					! Particle metallicity
					if(metal) then
					   zp(ipart)  = zz(i)
					endif
				 endif
				 if(type_index.gt.2)then
					if(star)then
					   typep(ipart)%family = FAM_STAR
					   typep(ipart)%tag    = 0
					end if
				 else if(type_index.eq.2)then
					typep(ipart)%family = FAM_DM
					typep(ipart)%tag    = 0
				 end if
				 up(ipart)      = uu(i)
				 if(ic_mask_ptype.gt.-1)then
					if(ic_mask_ptype.eq.type_index-1)then
					   maskp(ipart) = 1.0
					else
					   maskp(ipart) = 0.0
					endif
				 endif
				 ! Add a gas particle outside the zoom region
				 if(cosmo) then
					do j=1,6
					   if(type_index.eq.cosmo_add_gas_index(j)) then
						  ! Add a gas particle
						  xp(ipart+1,1:3) = xp(ipart,1:3)
						  vp(ipart+1,1:3) = vp(ipart,1:3)
						  idp(ipart+1)    = -1
						  mp(ipart+1)     = mp(ipart)*(omega_b/omega_m)
						  levelp(ipart+1) = levelmin
						  up(ipart+1)     = T2_start/scale_T2
						  if(metal) then
							 zp(ipart+1)  = z_ave*0.02
						  endif
						  ! Remove mass from the DM particle
						  mp(ipart) = mp(ipart)-mp(ipart+1)
						  ! Update index
						  ipart           = ipart+1
					   endif
					end do
				 endif
			  endif
#ifndef WITHOUTMPI
		   endif
#endif
		enddo
		lpart = lpart+jpart
	 enddo
	 if(myid==1)then
		write(*,'(A,E10.3,A)') ' Gas mass in AMR grid -> ',mgas_tot,' unit_m'
		write(*,'(A50)')"__________________________________________________"
		close(1)
	 endif
  enddo
  npart = ipart
  ! Compute total number of particle
  npart_cpu       = 0
  npart_all       = 0
  npart_cpu(myid) = npart
#ifndef WITHOUTMPI
  call MPI_ALLREDUCE(npart_cpu,npart_all,ncpu,MPI_INTEGER,MPI_SUM,MPI_COMM_WORLD,info)
  npart_cpu(1) = npart_all(1)
#else
  npart_all       = npart
#endif
  if(myid==1)then
	 write(*,*) ' npart_tot -> ',sum(npart_all)
	 write(*,'(A50)')"__________________________________________________"
	 close(1)
  endif
  do icpu=2,ncpu
	 npart_cpu(icpu)=npart_cpu(icpu-1)+npart_all(icpu)
  end do
  if(debug)write(*,*)'npart=',npart,'/',npart_cpu(ncpu)
  ifout = ic_ifout
  t = ic_t_restart
  ! DICE patch
end subroutine load_dice_restart
