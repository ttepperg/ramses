! Patch restart
 case ('ramses')

    ! Initialisation
    restart_init = .true.
    dice_init    = .true.			!-> VERY IMPORTANT: switch on DICE patch
!         amr_struct   = .false.		!-> does not conserve mass but works better; do NOT hardcode, though
    eocpu        = .false.
    error        = .false.
    icpu         = 1
    lpart_restart        = 0
    ipart        = 0
    mhalo_tot    = 0.
    mstar_tot    = 0.
    mgas_tot_restart     = 0.
    nhalo_tot    = 0
    nstar_tot    = 0
    if(myid==1) then
       write(*,'(A50)')"__________________________________________________"
       write(*,*)" RAMSES restart"
       write(*,'(A50)')"__________________________________________________"
    endif
    do while(.not.eocpu)
       nstar_loc = 0
       nhalo_loc = 0
       ngas_loc = 0
       nsink_loc = 0
       if(myid==1)then
          call title(abs(nrestart),nchar)
          if(icpu==1) write(*,'(a)') " Loading -> "//TRIM(output_dir)//'output_'//TRIM(nchar)
          fileloc=TRIM(output_dir)//'output_'//TRIM(nchar)//'/part_'//TRIM(nchar)//'.out'
          call title(icpu,nchar)
          fileloc=TRIM(fileloc)//TRIM(nchar)
          INQUIRE(file=fileloc,exist=ok)
          if(.not.ok)then
             write(*,*) TRIM(fileloc),' not found'
             call clean_stop
          endif
          ilun1 = 1
          OPEN(unit=ilun1,file=fileloc,status='old',action='read',form='unformatted',access="stream")

          mypos = 1

          read(ilun1,pos=mypos)size_blck_restart; mypos=mypos+sizeof(dummy_int)
          read(ilun1,pos=mypos)header_part%ncpu; mypos=mypos+sizeof(dummy_int)+size_blck_restart
          read(ilun1,pos=mypos)size_blck_restart; mypos=mypos+sizeof(dummy_int)
          read(ilun1,pos=mypos)header_part%ndim; mypos=mypos+sizeof(dummy_int)+size_blck_restart
          read(ilun1,pos=mypos)size_blck_restart; mypos=mypos+sizeof(dummy_int)
          read(ilun1,pos=mypos)header_part%npart; mypos=mypos+sizeof(dummy_int)+size_blck_restart
          read(ilun1,pos=mypos)size_blck_restart; mypos=mypos+sizeof(dummy_int)
          read(ilun1,pos=mypos)header_part%localseed; mypos=mypos+sizeof(dummy_int)+size_blck_restart
          read(ilun1,pos=mypos)size_blck_restart; mypos=mypos+sizeof(dummy_int)
          read(ilun1,pos=mypos)header_part%nstar_tot; mypos=mypos+sizeof(dummy_int)+size_blck_restart
          read(ilun1,pos=mypos)size_blck_restart; mypos=mypos+sizeof(dummy_int)
          read(ilun1,pos=mypos)header_part%mstar_tot; mypos=mypos+sizeof(dummy_int)+size_blck_restart
          read(ilun1,pos=mypos)size_blck_restart; mypos=mypos+sizeof(dummy_int)
          read(ilun1,pos=mypos)header_part%mstar_lost; mypos=mypos+sizeof(dummy_int)+size_blck_restart
          read(ilun1,pos=mypos)size_blck_restart; mypos=mypos+sizeof(dummy_int)
          read(ilun1,pos=mypos)header_part%nsink; mypos=mypos+sizeof(dummy_int)+size_blck_restart

          ! Init block address
          read(ilun1,pos=mypos) size_blck_restart; mypos = mypos+sizeof(dummy_int)
          pos_blck_restart(1) = mypos; mypos = mypos+size_blck_restart+sizeof(dummy_int)

          read(ilun1,pos=mypos) size_blck_restart; mypos = mypos+sizeof(dummy_int)
          pos_blck_restart(2) = mypos; mypos = mypos+size_blck_restart+sizeof(dummy_int)

          read(ilun1,pos=mypos) size_blck_restart; mypos = mypos+sizeof(dummy_int)
          pos_blck_restart(3) = mypos; mypos = mypos+size_blck_restart+sizeof(dummy_int)

          read(ilun1,pos=mypos) size_blck_restart; mypos = mypos+sizeof(dummy_int)
          vel_blck_restart(1) = mypos; mypos = mypos+size_blck_restart+sizeof(dummy_int)

          read(ilun1,pos=mypos) size_blck_restart; mypos = mypos+sizeof(dummy_int)
          vel_blck_restart(2) = mypos; mypos = mypos+size_blck_restart+sizeof(dummy_int)

          read(ilun1,pos=mypos) size_blck_restart; mypos = mypos+sizeof(dummy_int)
          vel_blck_restart(3) = mypos; mypos = mypos+size_blck_restart+sizeof(dummy_int)

          read(ilun1,pos=mypos) size_blck_restart; mypos = mypos+sizeof(dummy_int)
          mass_blck_restart = mypos; mypos = mypos+size_blck_restart+sizeof(dummy_int)

          read(ilun1,pos=mypos) size_blck_restart; mypos = mypos+sizeof(dummy_int)
          id_blck_restart = mypos; mypos = mypos+size_blck_restart+sizeof(dummy_int)

          read(ilun1,pos=mypos) size_blck_restart; mypos = mypos+sizeof(dummy_int)
          level_blck_restart = mypos; mypos = mypos+size_blck_restart+sizeof(dummy_int)

          if(star.or.sink)then
             read(ilun1,pos=mypos) size_blck_restart; mypos = mypos+sizeof(dummy_int)
             age_blck_restart = mypos; mypos = mypos+size_blck_restart+sizeof(dummy_int)
             if(metal) then
                read(ilun1,pos=mypos) size_blck_restart; mypos = mypos+sizeof(dummy_int)
                metal_blck_restart = mypos
             endif
          endif
       endif

!	TTG 2018: Loop to read particles
       eob_restart      = .false.
       kpart_restart    = 0
       do while(.not.eob_restart)
          xx=0.
          vv=0.
          ii=0.
          mm=0.
          tt=0.
          zz=0.
			  !ADDED BY TTG MAY 2017
			  uu    = 0.
          if(myid==1)then
             jpart=0
             do i=1,nvector
                jpart=jpart+1
                ! All particles counter
                kpart_restart=kpart_restart+1
                ! Reading ramses part file line-by-line
                do idim=1,ndim
                   read(ilun1,pos=pos_blck_restart(idim)+sizeof(dummy_real_restart)*(kpart_restart-1)) xx(jpart,idim)
                end do
                do idim=1,ndim
                   read(ilun1,pos=vel_blck_restart(idim)+sizeof(dummy_real_restart)*(kpart_restart-1)) vv(jpart,idim)
                end do
                read(ilun1,pos=mass_blck_restart+sizeof(dummy_real_restart)*(kpart_restart-1)) mm(jpart)
                read(ilun1,pos=id_blck_restart+sizeof(dummy_int)*(kpart_restart-1)) ii(jpart)
                if(star.or.sink) then
                   read(ilun1,pos=age_blck_restart+sizeof(dummy_real_restart)*(kpart_restart-1)) tt(jpart)
                   if(metal) then
                      read(ilun1,pos=metal_blck_restart+sizeof(dummy_real_restart)*(kpart_restart-1)) zz(jpart)
                   endif
                endif
                ! Updating total masses
                if(tt(jpart)==0d0) then
                   mhalo_tot = mhalo_tot+mm(jpart)
                   nhalo_tot = nhalo_tot+1
                   nhalo_loc = nhalo_loc+1
                else
                   mstar_tot = mstar_tot+mm(jpart)
                   nstar_tot = nstar_tot+1
                   nstar_loc = nstar_loc+1
                endif
                ! Check the End Of Block
                if(kpart_restart.ge.header_part%npart) then
                   eob_restart=.true.
                   exit
                endif
             enddo
          endif
#ifndef WITHOUTMPI
          call MPI_BCAST(eob_restart,1        ,MPI_LOGICAL         ,0,MPI_COMM_WORLD,info)
          call MPI_BCAST(xx,nvector*3         ,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,info)
          call MPI_BCAST(vv,nvector*3         ,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,info)
          call MPI_BCAST(ii,nvector           ,MPI_INTEGER         ,0,MPI_COMM_WORLD,info)
          call MPI_BCAST(mm,nvector           ,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,info)
          call MPI_BCAST(zz,nvector           ,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,info)
          call MPI_BCAST(tt,nvector           ,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,info)
			  !ADDED BY TTG MAY 2017
          call MPI_BCAST(uu,nvector           ,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,info)
          call MPI_BCAST(jpart,1              ,MPI_INTEGER         ,0,MPI_COMM_WORLD,info)
          call MPI_BARRIER(MPI_COMM_WORLD,info)
          call cmp_cpumap(xx,cc,jpart)
#endif
          do i=1,jpart
#ifndef WITHOUTMPI
              ! Check the CPU map
              if(cc(i)==myid)then
#endif
						!HERE: transform particle positions and velocities
                 xx(i,1:3) = xx(i,1:3)-restart_ic_center(1:3)

                 if(xx(i,1).ge.0d0.and.xx(i,1).le.boxlen.and. &
                        & xx(i,2).ge.0d0.and.xx(i,2).le.boxlen.and. &
                        & xx(i,3).ge.0d0.and.xx(i,3).le.boxlen) then
                    ipart          = ipart+1
                    if(ipart.gt.npartmax) then
                       write(*,*) "Increase npartmax"
                       error=.true.
#ifndef WITHOUTMPI
                       call MPI_BCAST(error,1,MPI_LOGICAL,0,MPI_COMM_WORLD,info)
#endif
                    endif
                    xp(ipart,1:3)  = xx(i,1:3)
                    vp(ipart,1:3)  = vv(i,1:3)
							!CHANGED BY TTG MAY 2017
!                         idp(ipart)     = ii(i)+1
                    idp(ipart)     = ii(i)
                    mp(ipart)      = mm(i)
							!ADDED BY TTG MAY 2017
							up(ipart)      = uu(i)
                    levelp(ipart)  = levelmin
                    if(star) then
                       tp(ipart)    = tt(i)
                       if(metal) then
                          zp(ipart) = zz(i)
                       endif
                    endif
                 endif
#ifndef WITOUTMPI
              endif
#endif
          enddo
#ifndef WITOUTMPI
          call MPI_BARRIER(MPI_COMM_WORLD,info)
#endif
          if(error) call clean_stop
       enddo

       if(myid==1)then
          ! Conversion factor from user units to cgs units
          call units(scale_l,scale_t,scale_d,scale_v,scale_nH,scale_T2)
          scale_m = scale_d*scale_l**3
          call title(abs(nrestart),nchar)
          fileloc=TRIM(output_dir)//'output_'//TRIM(nchar)//'/amr_'//TRIM(nchar)//'.out'
          call title(icpu,nchar)
          fileloc=TRIM(fileloc)//TRIM(nchar)
          INQUIRE(file=fileloc,exist=ok)
          if(.not.ok)then
             write(*,*) TRIM(fileloc),' not found'
             call clean_stop
          endif
          ilun2 = 2
          OPEN(unit=ilun2,file=fileloc,status='old',action='read',form='unformatted',access="stream")
          mypos=1
          read(ilun2,pos=mypos)size_blck_restart; mypos=mypos+sizeof(dummy_int)
          read(ilun2,pos=mypos)header_amr%ncpu; mypos=mypos+sizeof(dummy_int)+size_blck_restart
          read(ilun2,pos=mypos)size_blck_restart; mypos=mypos+sizeof(dummy_int)
          read(ilun2,pos=mypos)header_amr%ndim; mypos=mypos+sizeof(dummy_int)+size_blck_restart
          read(ilun2,pos=mypos)size_blck_restart; mypos=mypos+sizeof(dummy_int)
          read(ilun2,pos=mypos)header_amr%nx,header_amr%ny,header_amr%nz; mypos=mypos+sizeof(dummy_int)+size_blck_restart
          read(ilun2,pos=mypos)size_blck_restart; mypos=mypos+sizeof(dummy_int)
          read(ilun2,pos=mypos)header_amr%nlevelmax; mypos=mypos+sizeof(dummy_int)+size_blck_restart
          read(ilun2,pos=mypos)size_blck_restart; mypos=mypos+sizeof(dummy_int)
          read(ilun2,pos=mypos)header_amr%ngridmax; mypos=mypos+sizeof(dummy_int)+size_blck_restart
          read(ilun2,pos=mypos)size_blck_restart; mypos=mypos+sizeof(dummy_int)
          read(ilun2,pos=mypos)header_amr%nboundary; mypos=mypos+sizeof(dummy_int)+size_blck_restart
          read(ilun2,pos=mypos)size_blck_restart; mypos=mypos+sizeof(dummy_int)
          read(ilun2,pos=mypos)header_amr%ngrid_current; mypos=mypos+sizeof(dummy_int)+size_blck_restart
          read(ilun2,pos=mypos)size_blck_restart; mypos=mypos+sizeof(dummy_int)
          read(ilun2,pos=mypos)header_amr%boxlen; mypos=mypos+sizeof(dummy_int)+size_blck_restart
          ! Read time variables
          read(ilun2,pos=mypos)size_blck_restart; mypos=mypos+sizeof(dummy_int)
          read(ilun2,pos=mypos)header_amr%noutput,header_amr%iout,header_amr%ifout; mypos=mypos+sizeof(dummy_int)+size_blck_restart
          read(ilun2,pos=mypos)size_blck_restart; mypos=mypos+sizeof(dummy_int)
          read(ilun2,pos=mypos)header_amr%tout(1:header_amr%noutput); mypos=mypos+sizeof(dummy_int)+size_blck_restart
          read(ilun2,pos=mypos)size_blck_restart; mypos=mypos+sizeof(dummy_int)
          read(ilun2,pos=mypos)header_amr%aout(1:header_amr%noutput); mypos=mypos+sizeof(dummy_int)+size_blck_restart

          ! Old output times
          read(ilun2,pos=mypos)size_blck_restart; mypos=mypos+sizeof(dummy_int)
          read(ilun2,pos=mypos)header_amr%t; mypos=mypos+sizeof(dummy_int)+size_blck_restart
          read(ilun2,pos=mypos)size_blck_restart; mypos=mypos+sizeof(dummy_int)
          read(ilun2,pos=mypos)header_amr%dtold(1:header_amr%nlevelmax); mypos=mypos+sizeof(dummy_int)+size_blck_restart
          read(ilun2,pos=mypos)size_blck_restart; mypos=mypos+sizeof(dummy_int)
          read(ilun2,pos=mypos)header_amr%dtnew(1:header_amr%nlevelmax); mypos=mypos+sizeof(dummy_int)+size_blck_restart
          read(ilun2,pos=mypos)size_blck_restart; mypos=mypos+sizeof(dummy_int)
          read(ilun2,pos=mypos)header_amr%nstep,header_amr%nstep_coarse; mypos=mypos+sizeof(dummy_int)+size_blck_restart
          read(ilun2,pos=mypos)size_blck_restart; mypos=mypos+sizeof(dummy_int)

          read(ilun2,pos=mypos)header_amr%einit,header_amr%mass_tot_0,header_amr%rho_tot; mypos=mypos+sizeof(dummy_int)+size_blck_restart
          read(ilun2,pos=mypos)size_blck_restart; mypos=mypos+sizeof(dummy_int)

          read(ilun2,pos=mypos)header_amr%omega_m,header_amr%omega_l,header_amr%omega_k, &
             & header_amr%omega_b,header_amr%h0,header_amr%aexp_ini,header_amr%boxlen_ini; mypos=mypos+sizeof(dummy_int)+size_blck_restart
          read(ilun2,pos=mypos)size_blck_restart; mypos=mypos+sizeof(dummy_int)
          read(ilun2,pos=mypos)header_amr%aexp,header_amr%hexp,header_amr%aexp_old,header_amr%epot_tot_int,header_amr%epot_tot_old; mypos=mypos+sizeof(dummy_int)+size_blck_restart
          read(ilun2,pos=mypos)size_blck_restart; mypos=mypos+sizeof(dummy_int)
          read(ilun2,pos=mypos)header_amr%mass_sph; mypos=mypos+sizeof(dummy_int)+size_blck_restart

          header_amr%twotondim=2**header_amr%ndim
          header_amr%twondim=2*header_amr%ndim

          ! Read levels variables
          if(icpu==1) then
             allocate(header_amr%headl(1:header_amr%ncpu,1:header_amr%nlevelmax))
             allocate(header_amr%taill(1:header_amr%ncpu,1:header_amr%nlevelmax))
             allocate(header_amr%numbl(1:header_amr%ncpu,1:header_amr%nlevelmax))
             allocate(header_amr%numbtot(1:10,1:header_amr%nlevelmax))
          endif

          read(ilun2,pos=mypos)size_blck_restart; mypos=mypos+sizeof(dummy_int)
          read(ilun2,pos=mypos)header_amr%headl(1:header_amr%ncpu,1:header_amr%nlevelmax); mypos=mypos+sizeof(dummy_int)+size_blck_restart
          read(ilun2,pos=mypos)size_blck_restart; mypos=mypos+sizeof(dummy_int)
          read(ilun2,pos=mypos)header_amr%taill(1:header_amr%ncpu,1:header_amr%nlevelmax); mypos=mypos+sizeof(dummy_int)+size_blck_restart
          read(ilun2,pos=mypos)size_blck_restart; mypos=mypos+sizeof(dummy_int)
          read(ilun2,pos=mypos)header_amr%numbl(1:header_amr%ncpu,1:header_amr%nlevelmax); mypos=mypos+sizeof(dummy_int)+size_blck_restart
          read(ilun2,pos=mypos)size_blck_restart; mypos=mypos+sizeof(dummy_int)
          read(ilun2,pos=mypos)header_amr%numbtot(1:10,1:header_amr%nlevelmax); mypos=mypos+sizeof(dummy_int)+size_blck_restart

          ! Read boundary linked list
          if(icpu==1) then
             allocate(header_amr%headb   (1:MAXBOUND,1:header_amr%nlevelmax))
             allocate(header_amr%tailb   (1:MAXBOUND,1:header_amr%nlevelmax))
             allocate(header_amr%numbb   (1:MAXBOUND,1:header_amr%nlevelmax))
             allocate(header_amr%boundary(1:MAXBOUND,1:header_amr%nlevelmax))
          endif

          if(simple_boundary)then
             read(ilun2,pos=mypos)size_blck_restart; mypos=mypos+sizeof(dummy_int)
             read(ilun2,pos=mypos)header_amr%headb(1:header_amr%nboundary,1:header_amr%nlevelmax); mypos=mypos+sizeof(dummy_int)+size_blck_restart
             read(ilun2,pos=mypos)size_blck_restart; mypos=mypos+sizeof(dummy_int)
             read(ilun2,pos=mypos)header_amr%tailb(1:header_amr%nboundary,1:header_amr%nlevelmax); mypos=mypos+sizeof(dummy_int)+size_blck_restart
             read(ilun2,pos=mypos)size_blck_restart; mypos=mypos+sizeof(dummy_int)
             read(ilun2,pos=mypos)header_amr%numbb(1:header_amr%nboundary,1:header_amr%nlevelmax); mypos=mypos+sizeof(dummy_int)+size_blck_restart
          end if
          ! Read free memory
          read(ilun2,pos=mypos)size_blck_restart; mypos=mypos+sizeof(dummy_int)
          read(ilun2,pos=mypos)header_amr%headf,header_amr%tailf,header_amr%numbf,header_amr%used_mem,header_amr%used_mem_tot; mypos=mypos+sizeof(dummy_int)+size_blck_restart
          ! Read cpu boundaries
          read(ilun2,pos=mypos)size_blck_restart; mypos=mypos+sizeof(dummy_int)
          read(ilun2,pos=mypos)header_amr%ordering; mypos=mypos+sizeof(dummy_int)+size_blck_restart

          header_amr%nbilevelmax=ceiling(log(dble(header_amr%ncpu))/log(2.0))
          header_amr%nbinodes=2**(nbilevelmax+1)-1
          header_amr%ndomain=header_amr%ncpu*overload

          if(icpu==1) allocate(header_amr%bound_key (0:header_amr%ndomain))

          if(header_amr%ordering=='bisection') then
             read(ilun2,pos=mypos)size_blck_restart; mypos=mypos+sizeof(dummy_int)
             read(ilun2,pos=mypos)header_amr%bisec_wall(1:header_amr%nbinodes); mypos=mypos+sizeof(dummy_int)+size_blck_restart
             read(ilun2,pos=mypos)size_blck_restart; mypos=mypos+sizeof(dummy_int)
             read(ilun2,pos=mypos)header_amr%bisec_next(1:header_amr%nbinodes,1:2); mypos=mypos+sizeof(dummy_int)+size_blck_restart
             read(ilun2,pos=mypos)size_blck_restart; mypos=mypos+sizeof(dummy_int)
             read(ilun2,pos=mypos)header_amr%bisec_indx(1:header_amr%nbinodes); mypos=mypos+sizeof(dummy_int)+size_blck_restart
             read(ilun2,pos=mypos)size_blck_restart; mypos=mypos+sizeof(dummy_int)
             read(ilun2,pos=mypos)header_amr%bisec_cpubox_min(1:header_amr%ncpu,1:header_amr%ndim); mypos=mypos+sizeof(dummy_int)+size_blck_restart
             read(ilun2,pos=mypos)size_blck_restart; mypos=mypos+sizeof(dummy_int)
             read(ilun2,pos=mypos)header_amr%bisec_cpubox_max(1:header_amr%ncpu,1:header_amr%ndim); mypos=mypos+sizeof(dummy_int)+size_blck_restart
          else
             read(ilun2,pos=mypos)size_blck_restart; mypos=mypos+sizeof(dummy_int)
             read(ilun2,pos=mypos)header_amr%bound_key(0:header_amr%ndomain); mypos=mypos+sizeof(dummy_int)+size_blck_restart
          endif
          ! Read coarse level
          ! Son array
          read(ilun2,pos=mypos)size_blck_restart; mypos=mypos+2*sizeof(dummy_int)+size_blck_restart
          ! Flag array
          read(ilun2,pos=mypos)size_blck_restart; mypos=mypos+2*sizeof(dummy_int)+size_blck_restart
          ! cpu_map array
          read(ilun2,pos=mypos)size_blck_restart; mypos=mypos+2*sizeof(dummy_int)+size_blck_restart

          if(icpu==1) then
             allocate(amr_pos_blck_restart(1:header_amr%ndim,1:header_amr%nlevelmax,1:header_amr%nboundary+header_amr%ncpu))
             allocate(amr_son_blck_restart(1:header_amr%nlevelmax,1:header_amr%nboundary+header_amr%ncpu,1:header_amr%twotondim))
          endif

          ! Read fine levels
          do ilevel=1,header_amr%nlevelmax
             do ibound=1,header_amr%nboundary+header_amr%ncpu
                if(ibound<=header_amr%ncpu)then
                   ncache=header_amr%numbl(ibound,ilevel)
                else
                   ncache=header_amr%numbb(ibound-header_amr%ncpu,ilevel)
                end if
                if(ncache>0)then
                   ! Read grid index
                   read(ilun2,pos=mypos)size_blck_restart; mypos=mypos+2*sizeof(dummy_int)+size_blck_restart
                   ! Read next index
                   read(ilun2,pos=mypos)size_blck_restart; mypos=mypos+2*sizeof(dummy_int)+size_blck_restart
                   ! Read prev index
                   read(ilun2,pos=mypos)size_blck_restart; mypos=mypos+2*sizeof(dummy_int)+size_blck_restart
                   ! Read grid center
                   do idim=1,header_amr%ndim
                      read(ilun2,pos=mypos)size_blck_restart; mypos=mypos+sizeof(dummy_int)
                      amr_pos_blck_restart(idim,ilevel,ibound) = mypos; mypos=mypos+sizeof(dummy_int)+size_blck_restart
                      if(size_blck_restart.ne.ncache*sizeof(dummy_real_restart)) then
                         write(*,*) "AMR -> Grid center block size does not correspond to ncache"
                         call clean_stop
                      endif
                   end do
                   ! Read father index
                   read(ilun2,pos=mypos)size_blck_restart; mypos=mypos+2*sizeof(dummy_int)+size_blck_restart
                   ! Read nbor index
                   do ind=1,header_amr%twondim
                      read(ilun2,pos=mypos)size_blck_restart; mypos=mypos+2*sizeof(dummy_int)+size_blck_restart
                   end do
                   ! Read son index
                   do ind=1,header_amr%twotondim
                      read(ilun2,pos=mypos)size_blck_restart; mypos=mypos+sizeof(dummy_int)
                      amr_son_blck_restart(ilevel,ibound,ind) = mypos; mypos=mypos+sizeof(dummy_int)+size_blck_restart
                      if(size_blck_restart.ne.ncache*sizeof(dummy_int)) then
                         write(*,*) "AMR -> Son index block size does not correspond to ncache"
                         call clean_stop
                      endif
                   end do
                   ! Read cpu map
                   do ind=1,header_amr%twotondim
                      read(ilun2,pos=mypos)size_blck_restart; mypos=mypos+2*sizeof(dummy_int)+size_blck_restart
                   end do
                   ! Read refinement map
                   do ind=1,header_amr%twotondim
                      read(ilun2,pos=mypos)size_blck_restart; mypos=mypos+2*sizeof(dummy_int)+size_blck_restart
                   end do
                endif
             end do
          end do

          call title(abs(nrestart),nchar)
          fileloc=TRIM(output_dir)//'output_'//TRIM(nchar)//'/hydro_'//TRIM(nchar)//'.out'
          call title(icpu,nchar)
          fileloc=TRIM(fileloc)//TRIM(nchar)
          INQUIRE(file=fileloc,exist=ok)
          if(.not.ok)then
             write(*,*) TRIM(fileloc),' not found'
             call clean_stop
          endif
          ilun3 = 3
          OPEN(unit=ilun3,file=fileloc,status='old',action='read',form='unformatted',access="stream")
          mypos=1
          read(ilun3,pos=mypos)size_blck_restart; mypos=mypos+sizeof(dummy_int)
          read(ilun3,pos=mypos)header_hydro%ncpu; mypos=mypos+sizeof(dummy_int)+size_blck_restart
          read(ilun3,pos=mypos)size_blck_restart; mypos=mypos+sizeof(dummy_int)
          read(ilun3,pos=mypos)header_hydro%nvar; mypos=mypos+sizeof(dummy_int)+size_blck_restart
          read(ilun3,pos=mypos)size_blck_restart; mypos=mypos+sizeof(dummy_int)
          read(ilun3,pos=mypos)header_hydro%ndim; mypos=mypos+sizeof(dummy_int)+size_blck_restart
          read(ilun3,pos=mypos)size_blck_restart; mypos=mypos+sizeof(dummy_int)
          read(ilun3,pos=mypos)header_hydro%nlevelmax; mypos=mypos+sizeof(dummy_int)+size_blck_restart
          read(ilun3,pos=mypos)size_blck_restart; mypos=mypos+sizeof(dummy_int)
          read(ilun3,pos=mypos)header_hydro%nboundary; mypos=mypos+sizeof(dummy_int)+size_blck_restart
          read(ilun3,pos=mypos)size_blck_restart; mypos=mypos+sizeof(dummy_int)
          read(ilun3,pos=mypos)header_hydro%gamma; mypos=mypos+sizeof(dummy_int)+size_blck_restart

          if(icpu==1) &
          &allocate(hydro_var_blck_restart(1:header_amr%nlevelmax,1:header_amr%nboundary+header_amr%ncpu,1:header_amr%twotondim,1:header_hydro%nvar))

          do ilevel=1,header_hydro%nlevelmax
             do ibound=1,header_hydro%nboundary+header_hydro%ncpu
                if(ibound<=header_hydro%ncpu)then
                   ncache=header_amr%numbl(ibound,ilevel)
                else
                   ncache=header_amr%numbb(ibound-header_hydro%ncpu,ilevel)
                end if
                read(ilun3,pos=mypos)size_blck_restart; mypos=mypos+sizeof(dummy_int)
                read(ilun3,pos=mypos)ilevel2; mypos=mypos+sizeof(dummy_int)+size_blck_restart

                read(ilun3,pos=mypos)size_blck_restart; mypos=mypos+sizeof(dummy_int)
                read(ilun3,pos=mypos)numbl2; mypos=mypos+sizeof(dummy_int)+size_blck_restart
                if(numbl2.ne.ncache)then
                   write(*,*)'File hydro.tmp is not compatible'
                   write(*,*)'Found   =',numbl2,' for level ',ilevel2
                   write(*,*)'Expected=',ncache,' for level ',ilevel
                end if
                if(ncache>0)then
                   ! Loop over cells
                   do ind=1,header_amr%twotondim
                      ! Read all hydro var
                      do ivar=1,header_hydro%nvar
                         read(ilun3,pos=mypos)size_blck_restart; mypos=mypos+sizeof(dummy_int)
                         hydro_var_blck_restart(ilevel,ibound,ind,ivar)=mypos
                         mypos=mypos+sizeof(dummy_int)+size_blck_restart
                         if(size_blck_restart.ne.ncache*sizeof(dummy_real_restart)) then
                            write(*,*) "HYDRO -> Var block size does not correspond to ncache"
                            call clean_stop
                         endif
                      end do
                   end do
                end if
             end do
          end do
          if(nvector/header_amr%twotondim==0) then
             write(*,*) "Recompile with a greater value for NVECTOR"
             call clean_stop
          endif
       endif

       if(icpu==1) then
#ifndef WITHOUTMPI
          call MPI_BCAST(header_amr%nlevelmax,1   ,MPI_INTEGER,0,MPI_COMM_WORLD,info)
          call MPI_BCAST(header_amr%ndim,1        ,MPI_INTEGER,0,MPI_COMM_WORLD,info)
          call MPI_BCAST(header_amr%twotondim,1   ,MPI_INTEGER,0,MPI_COMM_WORLD,info)
          call MPI_BCAST(header_amr%ncpu,1        ,MPI_INTEGER,0,MPI_COMM_WORLD,info)
          call MPI_BCAST(header_amr%nboundary,1   ,MPI_INTEGER,0,MPI_COMM_WORLD,info)

          call MPI_BCAST(header_amr%ifout,1       ,MPI_INTEGER,0,MPI_COMM_WORLD,info)
          call MPI_BCAST(header_amr%t,1           ,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,info)
			  !ADDED BY TTG MAY 2017
          call MPI_BCAST(header_amr%iout,1        ,MPI_INTEGER,0,MPI_COMM_WORLD,info)
!               call MPI_BCAST(header_amr%tout,MAXOUT   ,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,info)
!               call MPI_BCAST(header_amr%aout,MAXOUT   ,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,info)
!               call MPI_BCAST(header_amr%noutput,1     ,MPI_INTEGER,0,MPI_COMM_WORLD,info)
          call MPI_BCAST(header_amr%nstep_coarse,1,MPI_INTEGER,0,MPI_COMM_WORLD,info)
          call MPI_BCAST(header_amr%nstep,1       ,MPI_INTEGER,0,MPI_COMM_WORLD,info)

          call MPI_BCAST(header_amr%boxlen,1      ,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,info)
			  !ADDED BY TTG MAY 2017
          call MPI_BCAST(header_amr%einit,1       ,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,info)
          call MPI_BCAST(header_amr%mass_tot_0,1  ,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,info)
          call MPI_BCAST(header_amr%rho_tot,1     ,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,info)

          call MPI_BCAST(header_hydro%nvar,1      ,MPI_INTEGER,0,MPI_COMM_WORLD,info)
          call MPI_BCAST(header_hydro%ndim,1      ,MPI_INTEGER,0,MPI_COMM_WORLD,info)
          if(cosmo) then
             call MPI_BCAST(header_amr%omega_m,1     ,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,info)
             call MPI_BCAST(header_amr%omega_l,1     ,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,info)
             call MPI_BCAST(header_amr%h0,1          ,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,info)
             call MPI_BCAST(header_amr%boxlen_ini,1  ,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,info)
             call MPI_BCAST(header_amr%aexp,1        ,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,info)
          endif
          call MPI_BARRIER(MPI_COMM_WORLD,info)
#endif

          if(myid==1) write(*,'(A50)')"__________________________________________________"
          ! Number of passive scalars to load (does NOT include internal energy)
          nvar_min       = min(header_hydro%nvar-header_hydro%ndim-1,nvar-ndim-1)-1
          if(myid==1) write(*,*) '[Info] Number of passive variables [restart_passive_var]: ', nvar_min
!              restart_boxlen = header_amr%boxlen
          if(cosmo) then
             omega_m = header_amr%omega_m
             omega_l = header_amr%omega_l
             h0 = header_amr%h0
             boxlen_ini = header_amr%boxlen_ini
             aexp = header_amr%aexp
             aexp_ini = aexp
          endif
          allocate(pasvar(1:nvector,1:nvar_min))
          ! Check passive hydro variables indices
!               if(myid==1) write(*,'(A50)')"__________________________________________________"
          do ivar=1,nvar_min
             if(myid==1) write(*,'(A,I2,A,I2)') ' [Info] Passive variable ',restart_passive_vars(ivar),' loaded in var',ndim+2+ivar
             if(restart_passive_vars(ivar).gt.header_hydro%nvar) then
                if(myid==1) write(*,*) '[Error] ivar=',restart_passive_vars(ivar),' is empty in the restart output'
                call clean_stop
             endif
!                  if(restart_passive_vars(ivar).lt.header_amr%ndim+2) then
             if(restart_passive_vars(ivar).le.header_amr%ndim+2) then  ! include internal energy / pressure
                if(myid==1) write(*,*) '[Error] ivar=',restart_passive_vars(ivar),' is an active variable'
                if(myid==1) write(*,*) '[Info] Additional variables have indices from ', ndim+2+nvar_min-1, 'to', header_hydro%nvar
                call clean_stop
             endif
          enddo
          if(myid==1) write(*,'(A50)')"__________________________________________________"


          ! (Re-)setting simulation time
          if(.not.reset_time) then
				  t              = header_amr%t
				  ifout          = header_amr%ifout
			     !ADDED BY TTG MAY 2017
! 				     noutput        = header_amr%noutput
				  iout           = header_amr%iout
! 					  tout(1:noutput)= header_amr%tout(1:header_amr%noutput)
! 					  aout(1:noutput)= header_amr%aout(1:header_amr%noutput)
! 					  if(ifout.gt.abs(nrestart)+1) ifout=abs(nrestart)+1
				  nstep_coarse   = header_amr%nstep_coarse
				  nstep          = header_amr%nstep
			  endif

			  ! (Re-)setting mass / energy / density conservation
			  !ADDED BY TTG MAY 2017
          if(.not.reset_cvars) then
				  einit     = header_amr%einit
				  mass_tot_0= header_amr%mass_tot_0
				  rho_tot   = header_amr%rho_tot
			  endif

          ! Compute movie frame number if applicable
          if(imovout>0) then
             do i=2,imovout
                if(aendmov>0)then
                   if(aexp>amovout(i-1).and.aexp<amovout(i)) then
                      imov=i
                   endif
                else
                   if(t>tmovout(i-1).and.t<tmovout(i)) then
                      imov=i
                   endif
                endif
             enddo
          endif

       endif

!	TTG 2018: Loop to map AMR grid cells onto pseudo particles
       do ilevel=1,header_amr%nlevelmax
          ! Mesh spacing in that level
          dx=0.5D0**ilevel
          nx_loc=(icoarse_max-icoarse_min+1)
          skip_loc=(/0.0d0,0.0d0,0.0d0/)
          if(ndim>0)skip_loc(1)=dble(icoarse_min)
          if(ndim>1)skip_loc(2)=dble(jcoarse_min)
          if(ndim>2)skip_loc(3)=dble(kcoarse_min)
          scale=header_amr%boxlen/dble(nx_loc)
          dx_loc=dx*scale
          vol_loc=dx_loc**header_amr%ndim
          do ind=1,header_amr%twotondim
             iz=(ind-1)/4
             iy=(ind-1-4*iz)/2
             ix=(ind-1-2*iy-4*iz)
             if(ndim>0)xc(ind,1)=(dble(ix)-0.5D0)*dx
             if(ndim>1)xc(ind,2)=(dble(iy)-0.5D0)*dx
             if(ndim>2)xc(ind,3)=(dble(iz)-0.5D0)*dx
          end do
          do ibound=1,header_amr%nboundary+header_amr%ncpu
             if(myid==1) then
                if(ibound<=header_amr%ncpu)then
                   ncache=header_amr%numbl(ibound,ilevel)
                else
                   ncache=header_amr%numbb(ibound-header_amr%ncpu,ilevel)
                end if
             endif
#ifndef WITHOUTMPI
             call MPI_BARRIER(MPI_COMM_WORLD,info)
             call MPI_BCAST(ncache,1,MPI_INTEGER,0,MPI_COMM_WORLD,info)
             call MPI_BARRIER(MPI_COMM_WORLD,info)
#endif
             if(ncache>0.and.ibound.eq.icpu) then
                kpart_restart = 0
                eob_restart = .false.
                do while(.not.eob_restart)
                   xx    = 0.
                   xxg   = 0.
                   vv    = 0.
                   mm    = 0.
                   tt    = 0.
                   pasvar    = 0.
                   uu    = 0.
                   jpart = 0
                   if(myid==1) then
                      do i=1,nvector/header_amr%twotondim
                         ! All particles counter
                         kpart_restart=kpart_restart+1
                         ! Read the grid center one time per twotondim cells
                         read_center=.true.
                         ! Loop over cells
                         do ind=1,header_amr%twotondim
                            read(ilun2,pos=amr_son_blck_restart(ilevel,ibound,ind)+sizeof(dummy_int)*(kpart_restart-1)) son1
                            ! Consider leaf cells only
                            if(son1==0) then
                               jpart = jpart+1
                               lpart_restart = lpart_restart+1
                               ! Reading grid center
										  ! TTG 2017: Converting to physical units is this taken care of when reading velocities below;
										  ! This needs to be done this way since reading will be executed only once per oct!
                               if(read_center) then
                                  do idim=1,ndim
                                     read(ilun2,pos=amr_pos_blck_restart(idim,ilevel,ibound)+sizeof(dummy_real_restart)*(kpart_restart-1)) xxg(idim)
                                   end do
                                  read_center=.false.
                               endif
                               ! Reading cell density
                               read(ilun3,pos=hydro_var_blck_restart(ilevel,ibound,ind,1)+sizeof(dummy_real_restart)*(kpart_restart-1)) mm(jpart)
                               if(mm(jpart).lt.restart_rho_min) restart_rho_min = mm(jpart)
                               if(mm(jpart).gt.restart_rho_max) restart_rho_max = mm(jpart)
                               ! Converting to mass
                               mm(jpart)=mm(jpart)*vol_loc
                               ! Updating total mass
                               mgas_tot_restart=mgas_tot_restart+mm(jpart)
                               ! Updating leaf cells counter
                               ngas_loc = ngas_loc+1
                               ! Reading velocities AND convert grid centre to physical units
                               do idim=1,ndim
                                  read(ilun3,pos=hydro_var_blck_restart(ilevel,ibound,ind,idim+1)+sizeof(dummy_real_restart)*(kpart_restart-1)) vv(jpart,idim)
                                  xx(jpart,idim)=(xxg(idim)+xc(ind,idim)-skip_loc(idim))*scale
                               end do
                               ! TTG MAY 2017: reading gas thermal pressure
                               read(ilun3,pos=hydro_var_blck_restart(ilevel,ibound,ind,header_amr%ndim+2)+sizeof(dummy_real_restart)*(kpart_restart-1)) uu(jpart)
                               ! TTG MAY 2017: transform to SPECIFIC internal energy
                               uu(jpart)=uu(jpart)/(header_hydro%gamma-1.d0)/(mm(jpart)/vol_loc)
                               ! Reading passive hydro variables
                               do ivar=1,nvar_min
                                 if(restart_passive_vars(ivar).gt.header_amr%ndim+2.and.restart_passive_vars(ivar).le.header_hydro%nvar) then
                                    read(ilun3,pos=hydro_var_blck_restart(ilevel,ibound,ind,restart_passive_vars(ivar))+sizeof(dummy_real_restart)*(kpart_restart-1))&
                                   &pasvar(jpart,ivar)
                                 	if(metal.and.((ndim+2+ivar)==imetal)) zz(jpart)   = pasvar(jpart,ivar)
                                   if((ndim+2+ivar)==itracer)            trc(jpart)  = pasvar(jpart,ivar)
                                 endif
                               enddo
                            endif
                         enddo
                         ! Check the End Of Block
                         if(kpart_restart.ge.ncache) then
                            eob_restart=.true.
                            exit
                         endif
                      enddo
                   endif
#ifndef WITHOUTMPI
                   call MPI_BARRIER(MPI_COMM_WORLD,info)
                   call MPI_BCAST(eob_restart,1                     ,MPI_LOGICAL         ,0,MPI_COMM_WORLD,info)
                   call MPI_BCAST(jpart,1                           ,MPI_INTEGER         ,0,MPI_COMM_WORLD,info)
                   call MPI_BCAST(kpart_restart,1                   ,MPI_INTEGER         ,0,MPI_COMM_WORLD,info)
                   call MPI_BCAST(xx,nvector*3                      ,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,info)
                   call MPI_BCAST(vv,nvector*3                      ,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,info)
                   call MPI_BCAST(mm,nvector                        ,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,info)
						  !ADDED BY TTG MAY 2017
                   call MPI_BCAST(uu,nvector                        ,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,info)
!                        call MPI_BCAST(pasvar,nvector*nvar_min           ,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,info)
                   call MPI_BCAST(zz,nvector                        ,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,info)
                   call MPI_BCAST(trc,nvector                       ,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,info)
                   call cmp_cpumap(xx,cc,jpart)
                   call MPI_BARRIER(MPI_COMM_WORLD,info)
#endif
                   do i=1,jpart
#ifndef WITHOUTMPI
                       ! Check the CPU map
                       if(cc(i)==myid)then
#endif
									!HERE: transform pseudo particle positions and velocities
                          xx(i,1:3) = xx(i,1:3)-restart_ic_center(1:3)

                          if(xx(i,1).ge.0d0.and.xx(i,1).le.boxlen.and. &
                           & xx(i,2).ge.0d0.and.xx(i,2).le.boxlen.and. &
                           & xx(i,3).ge.0d0.and.xx(i,3).le.boxlen) then
                             ipart          = ipart+1
                             if(ipart.gt.npartmax) then
                                write(*,*) "Increase npartmax"
                                error=.true.
#ifndef WITHOUTMPI
                                call MPI_BCAST(error,1,MPI_LOGICAL,0,MPI_COMM_WORLD,info)
#endif
                             endif
                             xp(ipart,1:3)  = xx(i,1:3)
                             vp(ipart,1:3)  = vv(i,1:3)
                             mp(ipart)      = mm(i)
						            !ADDED BY TTG MAY 2017
                             up(ipart)      = uu(i)
                             ! Flagged as gas particles
!                                  idp(ipart)     = 1
							         !CHANGED BY TTG MAY 2017
                             ! Flagged as gas particles according to myDICE (gas -> even id; collisionless -> odd id)
                             idp(ipart)     = 2
						            !CHANGED BY TTG MAY 2017
!                                  levelp(ipart)  = ilevel
                             levelp(ipart)  = levelmin
!                                  do ivar=1,nvar_min
!                                     if(metal.and.((ndim+2+ivar)==imetal)) zp(ipart)   = pasvar(i,ivar)
!                                     if((ndim+2+ivar)==itracer)            trcp(ipart) = pasvar(i,ivar)
!                                  end do
                             zp(ipart)   = zz(i)
                             trcp(ipart) = trc(i)
                          endif
#ifndef WITOUTMPI
                       endif
#endif
                   enddo
#ifndef WITHOUTMPI
                   call MPI_BARRIER(MPI_COMM_WORLD,info)
#endif
                   if(error) call clean_stop
                enddo
             endif
          enddo
       enddo

       if(myid==1) then
          call title(abs(nrestart),nchar)
          fileloc=TRIM(output_dir)//'output_'//TRIM(nchar)//'/sink_'//TRIM(nchar)//'.out'
          call title(icpu,nchar)
          fileloc=TRIM(fileloc)//TRIM(nchar)
          INQUIRE(file=fileloc,exist=ok)
          if(sink.and.ok)then
             ilun4 = 4
             OPEN(unit=ilun4,file=fileloc,status='old',action='read',form='unformatted',access="stream")
             mypos=1
             read(ilun4,pos=mypos)size_blck_restart; mypos=mypos+sizeof(dummy_int)
             read(ilun4,pos=mypos)header_sink%nsink; mypos=mypos+sizeof(dummy_int)+size_blck_restart
             read(ilun4,pos=mypos)size_blck_restart; mypos=mypos+sizeof(dummy_int)
             read(ilun4,pos=mypos)header_sink%nindsink; mypos=mypos+sizeof(dummy_int)+size_blck_restart

             if(header_sink%nsink>0)then
                ! Sink mass
                read(ilun4,pos=mypos) size_blck_restart; mypos = mypos+sizeof(dummy_int)
                sink_mass_blck_restart = mypos; mypos = mypos+size_blck_restart+sizeof(dummy_int)
                ! Sink birth epoch
                read(ilun4,pos=mypos) size_blck_restart; mypos = mypos+sizeof(dummy_int)
                sink_birth_blck_restart = mypos; mypos = mypos+size_blck_restart+sizeof(dummy_int)
                ! Sink position
                do idim=1,ndim
                   read(ilun4,pos=mypos) size_blck_restart; mypos = mypos+sizeof(dummy_int)
                   sink_pos_blck_restart(idim) = mypos; mypos = mypos+size_blck_restart+sizeof(dummy_int)
                end do
                ! Sink velocitiy
                do idim=1,ndim
                   read(ilun4,pos=mypos) size_blck_restart; mypos = mypos+sizeof(dummy_int)
                   sink_vel_blck_restart(idim) = mypos; mypos = mypos+size_blck_restart+sizeof(dummy_int)
                end do
                ! Sink angular momentum
                do idim=1,ndim
                   read(ilun4,pos=mypos) size_blck_restart; mypos = mypos+sizeof(dummy_int)
                   sink_am_blck_restart(idim) = mypos; mypos = mypos+size_blck_restart+sizeof(dummy_int)
                end do
                ! Sink accumulated rest mass energy
                read(ilun4,pos=mypos) size_blck_restart; mypos = mypos+sizeof(dummy_int)
                sink_dm_blck_restart = mypos; mypos = mypos+size_blck_restart+sizeof(dummy_int)
                ! Sink accretion rate
                read(ilun4,pos=mypos) size_blck_restart; mypos = mypos+sizeof(dummy_int)
                sink_ar_blck_restart = mypos; mypos = mypos+size_blck_restart+sizeof(dummy_int)
                ! Sink index
                read(ilun4,pos=mypos) size_blck_restart; mypos = mypos+sizeof(dummy_int)
                sink_index_blck_restart = mypos; mypos = mypos+size_blck_restart+sizeof(dummy_int)
                ! Sink new born boolean
                read(ilun4,pos=mypos) size_blck_restart; mypos = mypos+sizeof(dummy_int)
                sink_newb_blck_restart = mypos; mypos = mypos+size_blck_restart+sizeof(dummy_int)
                ! Sink int level
                read(ilun4,pos=mypos) sinkint_level; mypos=mypos+sizeof(dummy_int)+size_blck_restart
             endif
             eob_restart      = .false.
             kpart_restart    = 0
             do while(.not.eob_restart)
                xx=0.
                vv=0.
                ii=0.
                mm=0.
                tt=0.
                zz=0.
                aa=0.
                dd=0.
                ll=0.
					  !ADDED BY TTG MAY 2017
                uu=0.
                nn=.false.
                if(myid==1)then
                   jpart=0
                   do i=1,nvector
                      jpart=jpart+1
                      ! All particles counter
                      kpart_restart=kpart_restart+1
                      ! Reading ramses sink file line-by-line
                      do idim=1,ndim
                         read(ilun4,pos=sink_pos_blck_restart(idim)+sizeof(dummy_real_restart)*(kpart_restart-1)) xx(jpart,idim)
                         read(ilun4,pos=sink_vel_blck_restart(idim)+sizeof(dummy_real_restart)*(kpart_restart-1)) vv(jpart,idim)
                         read(ilun4,pos=sink_am_blck_restart(idim)+sizeof(dummy_real_restart)*(kpart_restart-1)) ll(jpart,idim)
                      end do
                      read(ilun4,pos=sink_mass_blck_restart+sizeof(dummy_real_restart)*(kpart_restart-1)) mm(jpart)
                      read(ilun4,pos=sink_index_blck_restart+sizeof(dummy_int)*(kpart_restart-1)) ii(jpart)
                      read(ilun4,pos=sink_dm_blck_restart+sizeof(dummy_real_restart)*(kpart_restart-1)) dd(jpart)
                      read(ilun4,pos=sink_ar_blck_restart+sizeof(dummy_real_restart)*(kpart_restart-1)) aa(jpart)
                      read(ilun4,pos=sink_newb_blck_restart+sizeof(dummy_logical)*(kpart_restart-1)) nn(jpart)
                      ! Updating total masses
                      if(tt(jpart)==0d0) then
                         msink_tot = msink_tot+mm(jpart)
                         nsink_tot = nsink_tot+1
                         nsink_loc = nsink_loc+1
                      endif
                      ! Check the End Of Block
                      if(kpart_restart.ge.header_part%npart) then
                         eob_restart=.true.
                         exit
                      endif
                   enddo
                endif
#ifndef WITHOUTMPI
                call MPI_BCAST(eob_restart,1        ,MPI_LOGICAL         ,0,MPI_COMM_WORLD,info)
                call MPI_BCAST(xx,nvector*3         ,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,info)
                call MPI_BCAST(vv,nvector*3         ,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,info)
                call MPI_BCAST(ll,nvector*3         ,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,info)
                call MPI_BCAST(ii,nvector           ,MPI_INTEGER         ,0,MPI_COMM_WORLD,info)
                call MPI_BCAST(mm,nvector           ,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,info)
                call MPI_BCAST(nn,nvector           ,MPI_LOGICAL         ,0,MPI_COMM_WORLD,info)
                call MPI_BCAST(tt,nvector           ,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,info)
					  !ADDED BY TTG MAY 2017
                call MPI_BCAST(uu,nvector           ,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,info)
                call MPI_BCAST(jpart,1              ,MPI_INTEGER         ,0,MPI_COMM_WORLD,info)
                call cmp_cpumap(xx,cc,jpart)
                call MPI_BARRIER(MPI_COMM_WORLD,info)
#endif
                do i=1,jpart
#ifndef WITHOUTMPI
                    ! Check the CPU map
                    if(cc(i)==myid)then
#endif
								!HERE: transform sink particle positions and velocities
                       xx(i,1:3) = xx(i,1:3)-restart_ic_center(1:3)

                       if(xx(i,1).ge.0d0.and.xx(i,1).le.boxlen.and. &
                              & xx(i,2).ge.0d0.and.xx(i,2).le.boxlen.and. &
                              & xx(i,3).ge.0d0.and.xx(i,3).le.boxlen) then
                          ipart          = ipart+1
                          if(ipart.gt.nsinkmax) then
                             write(*,*) "Increase nsinkmax"
                             error=.true.
#ifndef WITHOUTMPI
                             call MPI_BCAST(error,1,MPI_LOGICAL,0,MPI_COMM_WORLD,info)
#endif
                          endif
                          xsink(ipart,1:3)  = xx(i,1:3)
                          vsink(ipart,1:3)  = vv(i,1:3)
                          lsink(ipart,1:3)  = ll(i,1:3)
                          msink(ipart)      = mm(i)
                          tsink(ipart)      = tt(i)
                          delta_mass(ipart) = dd(i)
                          acc_rate(ipart)   = aa(i)
                          idsink(ipart)     = ii(i)
                          new_born(ipart)   = nn(i)
                       endif
#ifndef WITOUTMPI
                    endif
#endif
                enddo
#ifndef WITOUTMPI
                call MPI_BARRIER(MPI_COMM_WORLD,info)
                if(error) call clean_stop
#endif
             enddo
             call compute_ncloud_sink
             if(ir_feedback)then
                do i=1,nsink
                   acc_lum(i)=ir_eff*acc_rate(i)*msink(i)/(5*6.955d10/scale_l)
                end do
             end if

          endif
       endif

       if(myid==1) then
          close(ilun1)
          close(ilun2)
          close(ilun3)
          close(ilun4)
          write(*,*) 'CPU',icpu,' -> [',ngas_loc,'cells/',nstar_loc,'stars/',nhalo_loc,'dm]'
          icpu = icpu+1
          if(icpu.gt.header_part%ncpu) eocpu=.true.
       endif
#ifndef WITHOUTMPI
       call MPI_BCAST(eocpu,1      ,MPI_LOGICAL,0,MPI_COMM_WORLD,info)
       call MPI_BCAST(icpu,1       ,MPI_INTEGER,0,MPI_COMM_WORLD,info)
       call MPI_BCAST(nstar_tot,1  ,MPI_INTEGER,0,MPI_COMM_WORLD,info)
       call MPI_BCAST(mstar_tot,1  ,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,info)

		  !ADDED BY TTG MAY 2017
       call MPI_BCAST(restart_rho_min,1  ,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,info)
       call MPI_BCAST(restart_rho_max,1  ,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,info)
#endif
    enddo
    if(myid==1) then
       write(*,'(a50)')"__________________________________________________"
       write(*,*)" RAMSES restart summary"
       write(*,'(A50)')"__________________________________________________"
       write(*,'(a32,i10)') 'Restart snapshot [ifout]: ',ifout-1
       write(*,'(a32,1pe10.4)') 'Restart time [t]: ',t
       write(*,'(a32,1pe10.4)') 'Next snapshot time [tout(iout)]: ',tout(iout)
       write(*,'(a32,i10)') 'Total snapshots [noutput]: ',noutput
       write(*,'(a32,i10)') 'Main step [nstep_coarse]: ',nstep_coarse
       write(*,'(a32,i10)') 'Fine step [nstep]: ',nstep
       write(*,'(a32,i10)') 'CPUs: ',header_amr%ncpu
       write(*,'(a32,i10)') 'DM particles: ',nhalo_tot
       write(*,'(a32,1pe10.2)') 'M_DM  [Msun]: ',mhalo_tot*scale_m/msol
       write(*,'(a32,i10)') 'STAR particles: ',nstar_tot
       write(*,'(a32,1pe10.2)') 'M_stars  [Msun]: ',mstar_tot*scale_m/msol
       if(hydro) write(*,'(a32,i10)') 'GAS particles [leaf cells]: ',lpart_restart
       if(hydro) write(*,'(a32,1pe10.2)') 'M_gas  [Msun]: ',mgas_tot_restart*scale_m/msol
       if(hydro) write(*,'(a32,1pe10.2)') 'rho_min  [cm^-3]: ',restart_rho_min
       if(hydro) write(*,'(a32,1pe10.2)') 'rho_max  [cm^-3]: ',restart_rho_max
       if(sink)  write(*,'(a32,i10)') 'SINK particles: ',nsink_tot
       if(sink)  write(*,'(a32,1pe10.2)') 'M_sink  [Msun]: ',msink_tot*scale_m/msol
       write(*,'(a50)')"__________________________________________________"
    endif
    npart = ipart

    ! Compute total number of particles
    npart_cpu       = 0
    npart_all       = 0
    npart_cpu(myid) = npart
#ifndef WITHOUTMPI
    call MPI_ALLREDUCE(npart_cpu,npart_all,ncpu,MPI_INTEGER,MPI_SUM,MPI_COMM_WORLD,info)
    npart_cpu(1) = npart_all(1)
#endif
    if(myid==1)then
       write(*,*) ' npart_tot -> ',sum(npart_all)
       write(*,'(A50)')"__________________________________________________"
       close(1)
    endif
    do icpu=2,ncpu
       npart_cpu(icpu)=npart_cpu(icpu-1)+npart_all(icpu)
    end do
!         if(debug)write(*,*)'RAMSES file: npart(myid)=',npart,'(',myid,') / ',npart_cpu(ncpu)

  ! RESTART patch
 !----------------------------------------------------------------------
