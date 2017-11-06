Program driver

 ! --------------------------------------------------------------------------
 ! MPI Version of the Cartesian Acoustic Sun Simulator. 
 ! Copyright 2006, Shravan Hanasoge
 
 ! Hansen Experimental Physics Laboratory
 ! 455 Via Palou way, Stanford
 ! CA 94305, USA
 ! Email: shravan@stanford.edu 
 ! ---------------------------------------------------------------------------
  use initialize
  use all_modules
  use derivatives
  use PHYSICS_RHS
  use time_step

  implicit none

    integer ii,jj,kk, nn,init,ierr
    real(dp) start_time,end_time,t1,T00,nor1
    real(dp) start_mp_time
    character(len=140) message
    real(dp) CFL, CFLred
    logical :: fail, fail_all

 fail = .false.
 fail_all = .false.

 ! Initializing the computation
 ! --------------------------------------------

 call INIT_MPI()
 call INITANDALLOCATE()

  ! Excitation series computation -----------------------------

  num_steps = FLOOR(DBLE(maxtime)*timestep/cadforcing) + 2
  num_steps = 2*(num_steps/2)+2

 start_mp_time = MPI_WTIME()
 start_time = MPI_WTIME()

 !!!!!!! SOUND-SPEED PERTURBATIONS HERE !!!!!!

 call INITIALISE_RHS()

 call Initialize_step()

 if (TWOD) then
   npulses = 100
 else
   npulses = 500
 endif
 
 allocate(psr_t0(npulses),psr_amp(npulses), psr_sr(npulses),psr_sz(npulses),psr_px(npulses),psr_py(npulses),psr_pz(npulses), &
  psr_t1(npulses),pulseend(npulses))

 end_time = MPI_WTIME()

 if (rank == 0) print *, 'Initialization took: ', (end_time -start_time)

 !--------------------------------------------- 
  if (RESTART) then
    init = time0*fulloutputsteps+1
    call read_in_full_state(time0)
  else
  init = 0
  endif

  ! Excitation series computation -----------------------------

  if (.not. pulse) then

   allocate(vr(dim1(rank),dim2(rank),num_steps))

   call readfits_xy_distribute(forcingfunc,vr,num_steps)

   time_old = 0
   vr = (randomexciteampl/dimc)*vr
 endif

 if (magnetic) then

  CFL = 100.

  do ii = 1,dim1(rank)
   do jj = 1, dim2(rank)
    do kk = 1,dim3(rank)-1
     CFL = min(diml*(z(kk+1)-z(kk))/sqrt(dimb**(2.0)*reduction(ii,jj,kk)*(box(ii,jj,kk)**2 + boy(ii,jj,kk)**2.+boz(ii,jj,kk)**2.) &
        /(4.*pi*dimrho*rho0(ii,jj,kk)) + (dimc*c_speed(ii,jj,kk))**(2.0)),&
        diml*(z(kk+1)-z(kk))/sqrt(dimb**(2.0)*reduction(ii,jj,kk+1)*(box(ii,jj,kk+1)**2 + boy(ii,jj,kk+1)**2.+boz(ii,jj,kk+1)**2.) &
        /(4.*pi*rho0(ii,jj,kk+1)*dimrho) + (dimc*c_speed(ii,jj,kk+1))**(2.0)),&
        CFL)
    enddo
   enddo
  enddo

 endif
 if (.not. magnetic) then
  CFL = 100.

  do ii = 1,dim1(rank)
   do jj = 1, dim2(rank)
    do kk = 1,dim3(rank)-1
     CFL = min(diml/dimc*(z(kk+1)-z(kk))/(c_speed(ii,jj,kk)),CFL)
    enddo
   enddo
  enddo

 endif

  call MPI_REDUCE(CFL,CFLred,1,MPI_DOUBLE_PRECISION, MPI_MIN,0,MPI_COMM_WORLD,ierr)

  if (rank ==0) then
  print *, 'CFL estimate', timestep/CFLred , 'Integration scheme can probaby handle 1.5'
  endif

   !----------------Driver Kernel----------------------------
   ! RMS_HIST STORES THE RMS TIME EVOLUTION OF v_z
   ! A VERY USEFUL INDICATOR OF STABILITY

  do nn = init,maxtime

    time = nn

  ! code to be timed
    call cpu_time(t1)
    T00 = t1
    if (.not. pulse) then
      call lagrange_interp ()
      forcing = forcing*source_damp(:,:,1)
    else
      forcing_p = 0.
   
      if (.not. PSRPULSE) then 
        do ii = 1,dim1(rank)
          do jj = 1,dim2(rank)
            do kk = 1,dim3(rank)
              forcing_p(ii,jj,kk) = pulse_amp/0.56458272_dp/dimc*&
              sin(2.0_dp*pi*(time*timestep-pulse_t1+2.0_dp*pulse_t0)/pulse_t0)*&
              exp(-1.0_dp*(time*timestep-pulse_t1)**2/(2.0_dp*pulse_st**2))*&
              exp(-1.0_dp*((x(ii)*xlength/rsun-pulse_px/diml)**2+(y(jj)*ylength/rsun&
              -pulse_py/diml)**2)/(2.0_dp*(pulse_sr/diml)**2))*&
              exp(-1.0_dp*((z(kk)-(pulse_pz/diml+1.0_dp))**2)/(2.0_dp*(pulse_sz/diml)**2))
            enddo
          enddo
        enddo

      else
        call RUN_PSR_PULSE (nn,init)

        forcing_p = forcing_p*source_damp*src_scale_beta
      endif
    endif

    call step()
    end_time = MPI_WTIME() 
    call cpu_time(t1) 
   
    nor1 = norm((rho0+a(:,:,:,1))*(a(:,:,:,2)**2+a(:,:,:,3)**2+a(:,:,:,4)**2)/2.0)*dimc**2*dimrho !* dimc
    end_time = MPI_WTIME()

  ! Timing and 'convergence' statistics
  ! The criterion used to see if the results are 'convergent' is the
  ! L2 norm of the radial velocity. Quite silly, really.

   if (rank ==0) then
    
    write(message,"(A11,I9.2,A12, F8.4,A12, F8.4,A11, E12.5)") 'Iteration: ',time,'  CPU_time: ',end_time-start_time,'  Real_time: ',time*timestep/3600,'  norm_KE: ',nor1

    print *,trim(message)
   endif
   
   if (nor1 .GT. 10.0**(3.0)) then
     if (rank ==0) print *, 'Norm of VZ is getting very large, something odd is happening!'
     fail = .true.
   endif

   if((minval(rho+rho0) .LE. 0.) .OR. (minval(p+p0) .LE. 0.)) then
     print *, 'pressure or density negative in core'
     fail = .true.
   endif

   start_time = MPI_WTIME()

   call MPI_ALLREDUCE(fail, fail_all, 1, MPI_LOGICAL, MPI_LOR, MPI_COMM_WORLD, ierr)

   if (fail_all) then

     if (rank == 0) print *, 'A crash has occured, saving and dumping'
     call write_out_full_state(fail_all)
     call MPI_BARRIER(MPI_COMM_WORLD, ierr)
     call MPI_FINALIZE(ierr)
     stop
   endif

   if ((FULLOUT .AND. (mod(time,fulloutputsteps) == 0)) .OR. (time .EQ. 0)) call write_out_full_state(fail_all)
   if (CHUNKOUT .AND. (mod(time,chunkoutputsteps) == 0) .and. (time*timestep .GE. minsavtime*3600) .and. (time .GT. init)) call write_out_chunk()
   if (time .NE. 0) call filter_vars()  
 
 enddo

 call MPI_BARRIER(MPI_COMM_WORLD, ierr) 
 call MPI_FINALIZE(ierr)


end program driver
