module all_modules

! --------------------------------------------------------------------------
! MPI Version of the Spherical Acoustic Sun Simulator.
! Copyright 2006, Shravan Hanasoge
                                                                                                                                                        
! Hansen Experimental Physics Laboratory
! 455 Via Palou way, Stanford
! CA 94305, USA
! Email: shravan@stanford.edu
! --------------------------------------------------------------------------

  USE H5LT
  
  implicit none
  include 'mpif.h'
  include 'fftw3.f'
 
! ----------------------------
!! Parameters Variables

  integer, parameter :: dp=kind(0.d0)           ! double precision
  integer, parameter :: sp=kind(0.0 )           ! single precision

  integer nx,ny,nz, timestamp_size,time0,nxpmlbot,nxpmltop,nxpml,nypmlbot, &
  nypmltop,nypml,nzpmlbot,nzpmltop,nzpml,pulse_dir,TWOP5D

  !! PML COefficients, contains 1-reflection coefficient, 2- decay value n, 3-
  !! kappa value, 4 -interior sponge strength
  real(dp),dimension(4) :: xpmlb,xpmlt,ypmlb,ypmlt,zpmlb,zpmlt

  real(dp) xlength,ylength,timestep, wall_time, solartime, fulloutputcad, chunkoutputcad,minsavtime, &
  sp_str,sp_length,sp_decay,pulse_t0,pulse_t1,pulse_st,pulse_px,pulse_py, &
  pulse_pz,pulse_sr, pulse_sz,cadforcing,randomexcitedepth,randomexcitesigma,randomexciteampl,pulse_amp, &
  maxva

  logical DAMPING, STABILISE, PULSE,MAGNETIC,TWOD,&
  USE_PML,USE_HPML,PERIODIC,PSRPULSE, QUIET,DISPL,FLOWS,RESTART,CHUNKOUT,FULLOUT
  character(LEN=200) :: forcingfunc, file_data, dir_bkg, dir_out,simname
  real(dp) s_z,t_z, s_xy, t_xy

  !! Chunk output stuff
  
  integer xcb,xct,ycb,yct,zcb,zct,my_xcb,my_xct,my_ycb,my_yct,my_zcb,my_zct, &
          nxc,nyc,nzc,my_nxc,my_nyc,my_nzc,my_offsetx,my_offsety,my_offsetz
  logical IAM_CHUNKOUT

! ---------------------------
  integer time,step_rk,dimfive,option
  integer MAXTIME,e_rad,num_steps, fulloutputsteps, chunkoutputsteps, zfreq_filtering,xyfreq_filtering

  real(dp)  pi,deltat,rsun, dimc, diml, dimrho
  real(dp),dimension(:), allocatable :: z,x,y, spongex, spongey, spongez,&
  cq,invrho0q, stretch, unstretch,xaxis,yaxis,zaxis

  parameter(rsun = 69598946770.0_dp, pi = 3.14159265358979_dp)
  real(dp) stretchx, stretchy
 
!! RK4 stuff
  real(dp), dimension(5) :: betas, optimals

!! MPI Related definitions  

  integer numtasks, rank, xrank,xtasks,yrank,ytasks,zrank,ztasks

!! Domain distribution 
!!
  integer, dimension(:), allocatable :: dim1, dim2, dim3,dim_damp
  integer, dimension(:), allocatable :: ystart, xstart, zstart,yend,xend,zend
  integer, dimension(3) :: block_size,coords_xyz, starts, &
  parcel_xlds, parcel_xldr,parcel_xlfs,parcel_xlfr,parcel_xrds, parcel_xrdr, &
  parcel_xrfs,parcel_xrfr, parcel_ylds,parcel_yldr,parcel_ylfs,parcel_ylfr,parcel_yrds, parcel_yrdr, &
  parcel_yrfs,parcel_yrfr,parcel_zlds, parcel_zldr,parcel_zlfs,parcel_zlfr,parcel_zrds, parcel_zrdr, &
  parcel_zrfs,parcel_zrfr
  integer, allocatable, dimension(:,:) :: coords_xyz_all
  integer, allocatable, dimension(:,:,:) :: rank_xyz_all
  integer, allocatable, dimension(:) :: mysource3D
  integer, allocatable, dimension(:,:) :: damp_2
  integer :: chunkpmlbot,chunkpmltop,datatype_zcol,allpmltop,allpmlbot,damp_1

  logical, dimension(3) :: mpi_periodic
  integer comm_cart_3D, xsplit_comm, ysplit_comm, zsplit_comm,comm_cart_savechunk

  logical :: IAM_XBOT,IAM_XTOP,IAM_YBOT,IAM_YTOP,IAM_ZBOT,IAM_ZTOP,IAM_NOTB
! Common to all calculations
! div - divergence of the velocity vector
! vr - forcing applied on the RHS of the radial momentum equation
! p0, rho0, T0, c_speed, and g  - pressure, density, temperature, sound speed
! and gravity of the background solar state
! gradp, gradvr - radial velocity and pressure gradients in 3-space 
! omega - vorticities vector
! a - array that contains the 5 variables in the calculation: density, radial
! velocity, latitudinal velocity, longitudinal velocity and pressure in that order
! temp_step, scr - scratch space arrays used in the time evolution algorithm

  real(dp), allocatable, dimension(:,:,:) :: div, vr, c2, src_scale_beta,spongexyz

  real(dp), allocatable, dimension(:,:) :: forcing
  real(dp), allocatable, dimension(:,:,:) :: source_damp,forcing_p

  real(dp), allocatable, dimension(:,:,:) :: p0, gradp0_z, c_speed, rho0, gradrho0_z, rhoinv, c2rho0, c2div_v0 ,reduction
  real(dp), allocatable, dimension(:) :: source_dep, g
  real(dp), allocatable, target, dimension(:,:,:,:) ::  gradp
  real(dp), allocatable, target, dimension(:,:,:,:) ::  a,temp_step,scr

  real(dp), pointer, dimension(:,:,:) :: rho,RHSv_z,    &
  RHSv_y,RHScont,p,RHSp,gradp_z,gradp_y,gradp_x,RHSv_x

  real(dp), allocatable, target, dimension(:,:,:) :: dvxdx, dvydy,dvzdz

!-----Quadratic and Lagrange interpolation stuff------!
  integer time1, time_old
  real(dp) x0, x1, x2, x3, x4, x5, x6
  real(dp), dimension(:,:), allocatable :: LC0, LC1, LC2, LC3, LC4, LC5, LC6

!-----Magnetic field stuff-----!

  real(dp) dimb
  real(dp), allocatable, dimension(:,:,:) :: box, boy, boz, gradp0_x, gradp0_y, &
              curlbox, curlboy, curlboz, curlbx, curlby, curlbz, gradrho0_x, gradrho0_y 
  real(dp), pointer, dimension(:,:,:) :: bx, by, bz, RHSb_x, RHSb_y, RHSb_z, v_x, v_y, v_z

!-----DISPLACEMENT--------!
  real(dp), allocatable, target, dimension(:,:,:,:) :: scratch
  real(dp), pointer, dimension(:,:,:) :: dxixdx, dxiydy, dxizdz, RHSxi_x, RHSxi_y, RHSxi_z, &
                                   xi_x, xi_y, xi_z

!-----FLOWS-----!

  real(dp), allocatable, target, dimension(:,:,:,:) :: v0, dv0x,dv0y,dv0z
  real(dp), pointer, dimension(:,:,:) :: v0_x, v0_y, v0_z,dv0xdx,dv0xdy,dv0xdz,dv0ydx,dv0ydy, &
                                    dv0ydz, dv0zdx,dv0zdy,dv0zdz

!------PML STUFF------!
  real(dp), allocatable, dimension(:,:,:) :: az,bzpml,ax,bxpml,ay,bypml,kappax,kappay,kappaz
  real(dp), target, allocatable, dimension(:,:,:,:) :: scrzpml, psizpml, zpmlvars, &
  scrxpml, psixpml, xpmlvars,scrypml, psiypml, ypmlvars
  real(dp), pointer, dimension(:,:,:) :: psizvz, psizp, psizdzbx, psizdzby, RHSpsizinductionbx,&
                                   psizinductionbx, psizinductionby, RHSpsizvz, RHSpsizp, &
                                   RHSpsizdzbx, RHSpsizdzby, RHSpsizinductionby
  real(dp), pointer, dimension(:,:,:) :: psiyvy, psiyp, psiydybx, psiydybz, RHSpsiyinductionbx,&
                                   psiyinductionbz, psiyinductionbx, RHSpsiyvy, RHSpsiyp, &
                                   RHSpsiydybx, RHSpsiydybz, RHSpsiyinductionbz
  real(dp), pointer, dimension(:,:,:) :: psixvx, psixp, psixdxby, psixdxbz, RHSpsixinductionby,&
                                   psixinductionbz, psixinductionby, RHSpsixvx, RHSpsixp, &
                                   RHSpsixdxbz, RHSpsixdxby, RHSpsixinductionbz

!-----Damiens explicit derivatives stuff---!
  real(dp) :: dpz(6,11),dpxy(6,11) !dpxy(3,5)
  real(dp), allocatable, dimension(:) :: dfz, dfxy

  integer :: g_xy_filt, g_z_filt,g_xy_der, g_z_der
  PARAMETER(g_xy_der = 5, g_z_der = 5)
!-----------------------------
!! PSR_STUFF

  integer :: npulses
  real(dp) :: psr_st
  real(dp),allocatable, dimension(:) :: psr_t0,psr_amp, psr_sr,psr_sz,psr_px,psr_py,psr_pz, &
  psr_t1,pulseend

   !-----DAMPING STUFF---!
   integer nk, norder, filteralloc
   parameter(nk = 10)
   real(dp), allocatable, dimension(:,:,:) :: transfermatrix
   complex*16, allocatable, dimension(:,:,:) :: transfertemp
   real(dp), allocatable, dimension(:,:) :: damping_rates,kay2d
   real(dp), allocatable, dimension(:,:,:,:,:) :: filt
   real(dp), allocatable, dimension(:,:,:,:) :: temp1
   real(dp), allocatable, dimension(:,:) :: afilt
   real(dp), allocatable, dimension(:) :: bfilt, cfilt
   real(dp) dfilt,  dampcoe
   integer*8 fftw_plan_fwd_2d,fftw_plan_inv_2d


contains
!==================================================================================

SUBROUTINE READ_IN_HDF5()

  integer(HID_T) :: file_id, dset_id  , fspace_id, mspace_id, fapl_id      , dxpl_id

  integer(HSIZE_T), dimension(3) :: count,xyzdims
  integer(HSSIZE_T), dimension(3) :: offset 
  integer(HSIZE_T), dimension(1) :: zdim
  integer,dimension(3) :: dims

  real(dp), dimension(nz) :: grav
 
  integer :: ierr, stat(MPI_STATUS_SIZE)
  integer :: ii,jj,kk
  real(dp) :: temp, temp2

  !! Open file and setup access over all processors
  call h5open_f(ierr)
  call h5pcreate_f(H5P_FILE_ACCESS_F, fapl_id, ierr)
  call h5pset_fapl_mpio_f(fapl_id, MPI_COMM_WORLD, MPI_INFO_NULL, ierr)
  call h5fopen_f(dir_bkg, H5F_ACC_RDONLY_F, file_id, ierr, access_prp = fapl_id)
  call h5pcreate_f(H5P_DATASET_XFER_F, dxpl_id, ierr)
  call h5pset_dxpl_mpio_f(dxpl_id, H5FD_MPIO_COLLECTIVE_F, ierr)

  !! get the attributes, for now only dimensions needed. I should make a check
  !! That this matches that entered in the parameters file
  call h5ltget_attribute_int_f(file_id,"/","xyzdims", dims, ierr)

  !! Size and offset of the array to be read in
  count(1) = dim1(rank)
  count(2) = dim2(rank)
  count(3) = dim3(rank)
  offset(1) = xstart(rank)-1
  offset(2) = ystart(rank)-1
  offset(3) = zstart(rank)-1

  xyzdims = dims
 
  call h5screate_simple_f(3, count, mspace_id, ierr)
  call h5dopen_f(file_id,"rho0",dset_id,ierr)
  call h5dget_space_f(dset_id, fspace_id, ierr)
  call h5sselect_hyperslab_f(fspace_id, H5S_SELECT_SET_F, offset, count, ierr)

  ! Read in Thermodynamic variables
  call h5dread_f(dset_id,H5T_NATIVE_DOUBLE,rho0(1:dim1(rank),1:dim2(rank),1:dim3(rank)),xyzdims,ierr,file_space_id=fspace_id,&
       mem_space_id = mspace_id, xfer_prp=dxpl_id)
  call h5dclose_f(dset_id, ierr)
  call h5dopen_f(file_id,"pre0",dset_id,ierr)            
  call h5dread_f(dset_id,H5T_NATIVE_DOUBLE,p0(1:dim1(rank),1:dim2(rank),1:dim3(rank)),xyzdims,ierr,file_space_id=fspace_id,&
       mem_space_id = mspace_id, xfer_prp=dxpl_id)
  call h5dclose_f(dset_id, ierr)
  call h5dopen_f(file_id,"cs0",dset_id,ierr)            
  call h5dread_f(dset_id,H5T_NATIVE_DOUBLE,c_speed(1:dim1(rank),1:dim2(rank),1:dim3(rank)),xyzdims,ierr,file_space_id=fspace_id,&
       mem_space_id = mspace_id, xfer_prp=dxpl_id)
  call h5dclose_f(dset_id, ierr)

  !! if magnetic read in magnetic field

  if (MAGNETIC) then
      call h5dopen_f(file_id,"bx0",dset_id,ierr)
      call h5dread_f(dset_id,H5T_NATIVE_DOUBLE,box(1:dim1(rank),1:dim2(rank),1:dim3(rank)),xyzdims,ierr,file_space_id=fspace_id,&
      mem_space_id = mspace_id, xfer_prp=dxpl_id)
      call h5dclose_f(dset_id, ierr)
      call h5dopen_f(file_id,"by0",dset_id,ierr)
      call h5dread_f(dset_id,H5T_NATIVE_DOUBLE,boy(1:dim1(rank),1:dim2(rank),1:dim3(rank)),xyzdims,ierr,file_space_id=fspace_id,&
      mem_space_id = mspace_id, xfer_prp=dxpl_id)
      call h5dclose_f(dset_id, ierr)
      call h5dopen_f(file_id,"bz0",dset_id,ierr)
      call h5dread_f(dset_id,H5T_NATIVE_DOUBLE,boz(1:dim1(rank),1:dim2(rank),1:dim3(rank)),xyzdims,ierr,file_space_id=fspace_id,&
       mem_space_id = mspace_id, xfer_prp=dxpl_id)
      call h5dclose_f(dset_id, ierr)
  endif

  !! if flows read in background flows
  if (FLOWS) then

  call h5dopen_f(file_id,"vx0",dset_id,ierr)
  call h5dread_f(dset_id,H5T_NATIVE_DOUBLE,v0_x(1:dim1(rank),1:dim2(rank),1:dim3(rank)),xyzdims,ierr,file_space_id=fspace_id,&
       mem_space_id = mspace_id, xfer_prp=dxpl_id)
  call h5dclose_f(dset_id, ierr)
  call h5dopen_f(file_id,"vy0",dset_id,ierr)
  call h5dread_f(dset_id,H5T_NATIVE_DOUBLE,v0_y(1:dim1(rank),1:dim2(rank),1:dim3(rank)),xyzdims,ierr,file_space_id=fspace_id,&
           mem_space_id = mspace_id, xfer_prp=dxpl_id)
  call h5dclose_f(dset_id, ierr)
  call h5dopen_f(file_id,"vz0",dset_id,ierr)
  call h5dread_f(dset_id,H5T_NATIVE_DOUBLE,v0_z(1:dim1(rank),1:dim2(rank),1:dim3(rank)),xyzdims,ierr,file_space_id=fspace_id,&
       mem_space_id = mspace_id, xfer_prp=dxpl_id)
  call h5dclose_f(dset_id, ierr)

  endif

  call h5dopen_f(file_id,"zaxis",dset_id,ierr)
  call h5dread_f(dset_id,H5T_NATIVE_DOUBLE,zaxis(:),zdim,ierr)
  call h5dclose_f(dset_id,ierr)

  call h5dopen_f(file_id,"grav",dset_id,ierr)
  call h5dread_f(dset_id,H5T_NATIVE_DOUBLE,grav(:),zdim,ierr)
  call h5dclose_f(dset_id,ierr)

  call h5sclose_f(fspace_id, ierr)
  call h5sclose_f(mspace_id, ierr)
  call h5pclose_f(fapl_id, ierr)   
  call h5pclose_f(dxpl_id, ierr)   
  call h5fclose_f(file_id,ierr)
  call h5close_f(ierr)

  g(:) = grav(zstart(rank):zend(rank))

  zaxis = 1.0_dp + zaxis*1.0e8_dp/Rsun
  z(:) = zaxis(zstart(rank):zend(rank))

  temp2 = (1.0_dp + randomexcitedepth/Rsun)
  temp = 1.0_dp

  do kk = 1,nz
    if (abs(zaxis(kk) - temp2) .LE. temp) then
      temp = abs(zaxis(kk) - temp2)
      e_rad = kk
    endif
  enddo

  do kk = 1,dim3(rank)
    source_dep(kk) = exp(-0.5_dp*(z(kk) - zaxis(e_rad))**2/2_dp/(randomexcitesigma/rsun)**2)
  enddo

  if ((coords_xyz(1) .EQ. 0) .AND. (coords_xyz(2) .EQ. 0) .AND. (coords_xyz(3) .EQ. 0)) then
    dimc = c_speed(1,1,1)
    dimrho = rho0(1,1,1)
    diml = rsun
    dimb = (4_dp*pi*dimc**2*dimrho)**0.5_dp
  endif

  if ((coords_xyz(1) .EQ. 0) .AND. (coords_xyz(2) .eq. 0)) then
    cq = c_speed(1,1,:)

    do ii = 0,block_size(1)-1
      do jj = 0,block_size(2)-1
        if (ii .EQ. 0 .AND. jj .EQ. 0) then

        else
          call MPI_SEND(cq(1),dim3(rank),MPI_DOUBLE_PRECISION,&
          rank_xyz_all(ii,jj,coords_xyz(3)),(1+ii)*(1+jj),comm_cart_3D,ierr)
        endif
      enddo
    enddo
  else
    call MPI_RECV(cq(1),dim3(rank),MPI_DOUBLE_PRECISION,&
    rank_xyz_all(0,0,coords_xyz(3)),(1+coords_xyz(2))*(1+coords_xyz(1)),comm_cart_3D,stat,ierr)
  endif

  
  if ((coords_xyz(1) .EQ. 0) .AND. (coords_xyz(2) .eq. 0)) then
    invrho0q = 2.5440390e-7_dp/rho0(1,1,:)

    do ii = 0,block_size(1)-1
      do jj = 0,block_size(2)-1
        if (ii .EQ. 0 .AND. jj .EQ. 0) then

        else
          call MPI_SEND(invrho0q(1),dim3(rank),MPI_DOUBLE_PRECISION,&
          rank_xyz_all(ii,jj,coords_xyz(3)),(1+ii)*(1+jj),comm_cart_3D,ierr)
        endif
      enddo
    enddo
  else
    call MPI_RECV(invrho0q(1),dim3(rank),MPI_DOUBLE_PRECISION,&
    rank_xyz_all(0,0,coords_xyz(3)),(1+coords_xyz(2))*(1+coords_xyz(1)),comm_cart_3D,stat,ierr)
  endif


  call MPI_BCAST(dimc,1,MPI_DOUBLE_PRECISION,rank_xyz_all(0,0,0),comm_cart_3D,ierr)
  call MPI_BCAST(dimrho,1,MPI_DOUBLE_PRECISION,rank_xyz_all(0,0,0),comm_cart_3D,ierr)
  call MPI_BCAST(diml,1,MPI_DOUBLE_PRECISION,rank_xyz_all(0,0,0),comm_cart_3D,ierr)
  call MPI_BCAST(dimb,1,MPI_DOUBLE_PRECISION,rank_xyz_all(0,0,0),comm_cart_3D,ierr)

  if (QUIET) then
    if (coords_xyz(1) .EQ. 0 .AND. coords_xyz(2) .eq. 0) then
      do ii = 0,block_size(1)-1
        do jj = 0,block_size(2)-1
          if (ii .EQ. 0 .AND. jj .EQ. 0) then
          else
            call MPI_SEND(p0(1,1,1),1,datatype_zcol,&
            rank_xyz_all(ii,jj,coords_xyz(3)),1,comm_cart_3D,ierr)
            call MPI_SEND(rho0(1,1,1),1,datatype_zcol,&
            rank_xyz_all(ii,jj,coords_xyz(3)),2,comm_cart_3D,ierr)
            call MPI_SEND(c_speed(1,1,1),1,datatype_zcol,&
            rank_xyz_all(ii,jj,coords_xyz(3)),3,comm_cart_3D,ierr)
          endif
        enddo
      enddo
    else
      call MPI_RECV(p0(1,1,1),1,datatype_zcol,&
      rank_xyz_all(0,0,coords_xyz(3)),1,comm_cart_3D,stat,ierr)
      call MPI_RECV(rho0(1,1,1),1,datatype_zcol,&
      rank_xyz_all(0,0,coords_xyz(3)),2,comm_cart_3D,stat,ierr)
      call MPI_RECV(c_speed(1,1,1),1,datatype_zcol,&
      rank_xyz_all(0,0,coords_xyz(3)),3,comm_cart_3D,stat,ierr)
    endif
    call MPI_BARRIER(comm_cart_3D,ierr)
    do ii = 1,dim1(rank)
      do jj = 1,dim2(rank)
        p0(ii,jj,:) = p0(1,1,:)
        rho0(ii,jj,:) = rho0(1,1,:)
        c_speed(ii,jj,:) = c_speed(1,1,:)
      enddo
    enddo
  endif

  allocate(src_scale_beta(dim1(rank),dim2(rank),dim3(rank)))

  if (MAGNETIC) then
    src_scale_beta = 0.01*8.0*pi*p0/MAX(Box**2 + Boy**2 + Boz**2,1.0e-20_dp)
    do ii = 1,dim1(rank)
      do jj = 1,dim2(rank)
        do kk = 1,dim3(rank)
          src_scale_beta(ii,jj,kk) = min(src_scale_beta(ii,jj,kk),1.0)
        enddo
      enddo
    enddo
  else
    src_scale_beta = 1.0
  endif

  g = g*diml/dimc**2.0
  cq = cq/dimc
  p0 = p0/(dimrho * dimc**2.0)
  rho0 = rho0/dimrho
  c_speed = c_speed/dimc
  if (MAGNETIC) then
    Box = Box/dimb
    Boy = Boy/dimb
    Boz = Boz/dimb
  endif

  c2=c_speed*c_speed

  call MPI_BARRIER(comm_cart_3D,ierr)

end subroutine READ_IN_HDF5


!================================================================================================

function norm2_mpi(matrix)
  
   implicit none
   integer i,j,ierr
   real(dp) matrix(dim1(rank),dim2(rank))
   real(dp) norm2_mpi, sum

   sum = 0.0  
   norm2_mpi  = 0.0
   do j =1,dim2(rank)
    do i =1,dim1(rank)
      sum = sum + matrix(i,j)**2.0
    end do     
   end do     

   call MPI_REDUCE(sum, norm2_mpi, 1, MPI_DOUBLE_PRECISION, MPI_SUM, 0, comm_cart_3D, ierr)

   norm2_mpi = (norm2_mpi/(DBLE(nx)*DBLE(ny)))**0.5

end function norm2_mpi

!================================================================================

subroutine printerror(status)

! See cookbook.f on FITSIO for details

  integer status
  character errtext*30,errmessage*80

  if (status .le. 0)return

  call ftgerr(status,errtext)
  print *,'FITSIO Error Status =',status,': ',errtext

  call ftgmsg(errmessage)
  do while (errmessage .ne. ' ')
    print *,errmessage
     call ftgmsg(errmessage)
  end do

 end subroutine printerror

!================================================================================

!! rewrite this to send dimz * dim1*dim2 squares with a stride of nz
!! Legacy fits, leftover to read in a forcing function. WIll be replaced one
!! day....

  subroutine readfits_xy_distribute(filename,read,dimz)

  integer status,unit,readwrite,blocksize,naxes(3), kk,ierr
  integer group,firstpix, stat(MPI_STATUS_SIZE,numtasks-1)
  integer stat1(MPI_STATUS_SIZE), req1, req(numtasks-1)
  integer*8 nelements
  integer dimz
  real(dp) nullval
  real(dp), intent(out) :: read(dim1(rank),dim2(rank),dimz)
  real(dp) :: temp(nx,ny,dimz)
  logical anynull
  character*(*) filename

  if (rank == 0) then
    status=0
    call ftgiou(unit,status)
    readwrite=0
    call ftopen(unit,filename,readwrite,blocksize,status)

    naxes(1) = nx
    naxes(2) = ny
    naxes(3) = dimz
    nelements=naxes(1)*naxes(2)*naxes(3)
    group=1
    firstpix=1
    nullval=-999

    call ftgpvd(unit,group,firstpix,nelements,nullval, &
    &            temp,anynull,status)
    call ftclos(unit, status)
    call ftfiou(unit, status)

    if (status .gt. 0) call printerror(status)
    
    read = temp(xstart(0):xend(0),ystart(0):yend(0),:) 

    do kk = 1,numtasks-1
      call MPI_ISEND(temp(xstart(kk),ystart(kk),1), 1, mysource3D(kk), &
      kk, kk, comm_cart_3D,req(kk), ierr)
    enddo
      call MPI_WAITALL(numtasks-1,req,stat,ierr)
  else
    call MPI_RECV(read, dim1(rank)*dim2(rank)*dimz, MPI_DOUBLE_PRECISION,&
    0, rank, comm_cart_3D, req1, ierr)
    call MPI_WAIT(req1,stat1,ierr)
  endif
   
end subroutine readfits_xy_distribute
 
!================================================================================
  
  subroutine RUN_PSR_PULSE (ii,init)

   integer :: ii,init,nn,jj,kk 
   real(dp),allocatable,dimension(:) :: mu_amp,sigma_amp,mu_sz, sigma_sz
   real(dp) :: mu_sr, sigma_sr,mu_freq,sigma_freq
   real(dp) :: norm

   ! Temp Array for calculating Uniforms + Normals
   real(dp),allocatable, dimension(:,:) :: T1

   allocate(mu_amp(npulses),sigma_amp(npulses),mu_sz(npulses),sigma_sz(npulses),T1(9,npulses))
   !! width of horizontal pulses based on horizontal resolution
   !! to ensure well resolved, with slight randomness is size.
   !! For high-def runs probably reccomend more pulses, since the gaussian
   !! will be smaller. I feel like I should also change this to be dependant on
   !! depth, as size of pulse will be proportional to length-scale. But for now
   !! they are resolution dependant. 

   if (TWOD) then
     mu_sr = 1.5_dp*(x(2)-x(1))
   else
     mu_sr = 1.5_dp*max(x(2)-x(1),y(2)-y(1))
   endif
   sigma_sr = 0.2_dp*mu_sr

   psr_st = 120_dp
   mu_freq = 3.3e-3_dp
   sigma_freq = 0.5e-3_dp

   if (periodic) then
     norm = 1.0_dp
   else
     norm = 10.0_dp
   endif

     call init_random_seed()
     call RANDOM_NUMBER(T1)

   !! IF timestep = 0, then initialize sources. fill array of npulses.
   if (ii == init) then

   ! Uniforms
   ! R3 and R4 between 0 and 1
   ! R5 uniform between -0.2 Mm and -1. Mm
         
     psr_px = nxpmlbot*1.0_dp/nx + T1(1,:)*(nx-nxpmltop-nxpmlbot)*1.0_dp/nx
   if (TWOD) then
     psr_py = 0.0
   else
     psr_py = nypmlbot*1.0_dp/ny + T1(2,:)*(ny-nypmltop-nypmlbot)*1.0_dp/ny
   endif  

     psr_pz = 1.0_dp - 0.1e8_dp/diml - T1(3,:)*(1.9e8_dp/diml)
     mu_sz = 0.05e8_dp/diml + 0.05e8_dp/diml/2.0e8_dp*(1.0_dp-psr_pz(:))
     sigma_sz = 0.2_dp*mu_sz

     psr_sr = mu_sr + sigma_sr*sqrt(-2.0*LOG(T1(4,:)))*cos(2.*pi*T1(5,:))
     psr_sz = mu_sz + sigma_sz*sqrt(-2.0*LOG(T1(4,:)))*sin(2.*pi*T1(5,:))

     mu_amp = 1.0e-1_dp
     sigma_amp = 0.2_dp*mu_amp
     psr_amp = mu_amp + sigma_amp*sqrt(-2.0_dp*LOG(T1(6,:)))*cos(2.0_dp*pi*T1(7,:))
 
     psr_t0 = 1_dp/(mu_freq + sigma_freq*sqrt(-2.0_dp*LOG(T1(8,:)))*cos(2.0_dp*pi*T1(9,:)))
     ! Randomise starting time over first 400 seconds
     psr_t1 = init*timestep+(3_dp+1_dp*T1(7,:))*psr_st
     ! Pulseend = the first iteration after some multiple of the pulses decay
     pulseend=1.0d100 !pulseend = psr_t1 + 3_dp*psr_st

     !! which pulses do i need does (p_x +/- 3*s_r,p_y +/- 3*s_r, p_z +/- 3*s_z)
     !! lie in my domain
     do nn = 1,npulses

     if(((MODULO(psr_px(nn)+3.0_dp*psr_sr(nn),1.0_dp) .GE. x(1)) .OR. (psr_px(nn)+3.0_dp*psr_sr(nn) .GE. xaxis(nx))) .AND. &
        ((MODULO(psr_px(nn)-3.0_dp*psr_sr(nn),1.0_dp) .LE. x(dim1(rank))) .OR. (psr_px(nn)-3.0_dp*psr_sr(nn) .LE. xaxis(1))) .AND. &
        ((MODULO(psr_py(nn)+3.0_dp*psr_sr(nn),1.0_dp) .GE. y(1)) .OR. (psr_py(nn)+3.0_dp*psr_sr(nn) .GE. yaxis(ny))) .AND. &
        ((MODULO(psr_py(nn)-3.0_dp*psr_sr(nn),1.0_dp) .LE. y(dim2(rank))) .OR. (psr_py(nn)-3.0_dp*psr_sr(nn) .LE. yaxis(1))) .AND. &
        (psr_pz(nn)+3.0_dp*psr_sz(nn) .GE. z(1)) .AND. (psr_pz(nn)-3.0_dp*psr_sz(nn) .LE. z(dim3(rank)))) then

       do kk = 1,dim3(rank)
         do jj = 1,dim2(rank)
           forcing_p(:,jj,kk) = forcing_p(:,jj,kk)+invrho0q(kk)*psr_amp(nn)*&
           sin(2.d0*pi*(time*timestep-psr_t1(nn))/psr_t0(nn))*&
           exp(-(time*timestep-psr_t1(nn))**2/(2.d0*psr_st**2))*&
           exp(-(min((x(:)-psr_px(nn))**2,(x(:)-1.0-psr_px(nn))**2,(x(:)+1.0-psr_px(nn))**2)+&
           min((y(jj)-psr_py(nn))**2,(y(jj)-1.0-psr_py(nn))**2,(y(jj)+1.0-psr_py(nn))**2))/(2.0_dp*psr_sr(nn)**2))*&
           exp(-(z(kk)-psr_pz(nn))**2/(2.0_dp*psr_sz(nn)**2))
         enddo
       enddo
     endif
     enddo
   endif

   do nn = 1,npulses
     if (pulseend(nn) .LE. ii*timestep) then

       psr_px(nn) = nxpmlbot*1.0_dp/nx + T1(1,nn)*(nx-nxpmltop-nxpmlbot)*1.0_dp/nx
      
       if (TWOD) then
         psr_py = 0.0_dp 
       else
         psr_py(nn) = nypmlbot*1.0_dp/ny + T1(2,nn)*(ny-nypmltop-nypmlbot)*1.0_dp/ny
       endif         

       psr_pz(nn) = 1.0_dp - 0.1e8_dp/diml - T1(3,nn)*(1.9e8_dp/diml)
       mu_sz(nn) = 0.05e8_dp/diml + 0.05e8_dp/diml/2.0e8_dp*(1.0_dp-psr_pz(nn))   
       sigma_sz(nn) = 0.1_dp*mu_sz(nn)

       psr_sr(nn) = mu_sr + sigma_sr*sqrt(-2.0*LOG(T1(4,nn)))*cos(2.*pi*T1(5,nn))
       psr_sz(nn) = mu_sz(nn) + sigma_sz(nn)*sqrt(-2.0*LOG(T1(4,nn)))*sin(2.*pi*T1(5,nn))

       mu_amp(nn) = 1.0e-1_dp
       sigma_amp(nn) = 0.1_dp*mu_amp(nn)

       psr_amp(nn) = mu_amp(nn) + sigma_amp(nn)*sqrt(-2.0_dp*LOG(T1(6,nn)))*cos(2.0_dp*pi*T1(7,nn))
 
       psr_t0(nn) = 1_dp/(mu_freq + sigma_freq*sqrt(-2.0_dp*LOG(T1(8,nn)))*cos(2.0_dp*pi*T1(9,nn)))

       !! New pulseend
       psr_t1(nn) = time*timestep + (3_dp+3_dp*T1(7,nn))*psr_st
       pulseend=1.0d100 !pulseend(nn) = psr_t1(nn) + 3_dp*psr_st
       endif

     if(((MODULO(psr_px(nn)+3.0_dp*psr_sr(nn),1.0_dp) .GE. x(1)) .OR. (psr_px(nn)+3.0_dp*psr_sr(nn) .GE. xaxis(nx))) .AND. &
        ((MODULO(psr_px(nn)-3.0_dp*psr_sr(nn),1.0_dp) .LE. x(dim1(rank))) .OR. (psr_px(nn)-3.0_dp*psr_sr(nn) .LE. xaxis(1))) .AND. & 
        ((MODULO(psr_py(nn)+3.0_dp*psr_sr(nn),1.0_dp) .GE. y(1)) .OR. (psr_py(nn)+3.0_dp*psr_sr(nn) .GE. yaxis(ny))) .AND. &
        ((MODULO(psr_py(nn)-3.0_dp*psr_sr(nn),1.0_dp) .LE. y(dim2(rank))) .OR. (psr_py(nn)-3.0_dp*psr_sr(nn) .LE. yaxis(1))) .AND. &
        (psr_pz(nn)+3.0_dp*psr_sz(nn) .GE. z(1)) .AND. (psr_pz(nn)-3.0_dp*psr_sz(nn) .LE. z(dim3(rank)))) then

     do kk = 1,dim3(rank)
       do jj = 1,dim2(rank)
         forcing_p(:,jj,kk) = forcing_p(:,jj,kk)+invrho0q(kk)*psr_amp(nn)*&
         sin(2.d0*pi*(time*timestep-psr_t1(nn))/psr_t0(nn))*&
         exp(-(time*timestep-psr_t1(nn))**2/(2.d0*psr_st**2))*&
         exp(-(min((x(:)-psr_px(nn))**2,(x(:)-1.0-psr_px(nn))**2,(x(:)+1.0-psr_px(nn))**2)+&
         min((y(jj)-psr_py(nn))**2,(y(jj)-1.0-psr_py(nn))**2,(y(jj)+1.0-psr_py(nn))**2))/(2.0_dp*psr_sr(nn)**2))*&
         exp(-(z(kk)-psr_pz(nn))**2/(2.0_dp*psr_sz(nn)**2))
       enddo
     enddo
    endif 
    enddo

end subroutine RUN_PSR_PULSE
!================================================================================
subroutine init_random_seed()

  integer :: clock,n,i,ierr
  integer, DIMENSION(:), ALLOCATABLE :: seed

  call RANDOM_SEED(size = n)
  ALLOCATE(seed(n))

  call system_clock(count=clock)

  call  MPI_BCAST(clock,1,MPI_INTEGER,0,comm_cart_3D,ierr)
  seed = clock + 37*(/(i-1,i=1,n)/)
  call RANDOM_SEED(PUT=seed)

  DEALLOCATE(seed)
end subroutine init_random_seed

!================================================================================

function norm(matrix)

   integer i, j, k,ierr
   real(dp), dimension(dim1(rank),dim2(rank),dim3(rank)) :: matrix
   real(dp) norm, sum
 
   norm  = 0.0
   sum = 0.0

   do k = 1,dim3(rank)
    do j =1,dim2(rank)
     do i =1,dim1(rank)
       sum = sum + matrix(i,j,k)**2.0
     end do     
    end do     
   end do 

   call MPI_REDUCE(sum, norm, 1, MPI_DOUBLE_PRECISION, MPI_SUM, 0, comm_cart_3D, ierr) 

   norm = (norm/(DBLE(nx)*DBLE(ny)*DBLE(nz)))**0.5_dp

   call MPI_BARRIER(comm_cart_3D, ierr)
end function norm

!================================================================================

 subroutine convert_to_string(numbe,sting,length_string)

  integer i,length_string,numbe,n(1:length_string),number_temp
  character*(length_string) sting
  character*1 charc(10)

  charc(1)  = '0'
  charc(2)  = '1'
  charc(3)  = '2'
  charc(4)  = '3'
  charc(5)  = '4'
  charc(6)  = '5'
  charc(7)  = '6'
  charc(8)  = '7'
  charc(9)  = '8'
  charc(10) = '9'


  number_temp = numbe
  do i=length_string,1,-1

    n(length_string-i+1) = floor(number_temp*10.0**(-(i-1.0)))
    number_temp = number_temp - n(length_string-i+1)*10**(i-1)
    sting(length_string-i+1:length_string-i+1) = charc(n(length_string-i+1)+1)

  enddo

  end subroutine convert_to_string

!================================================================================
 
  subroutine READ_IN_FULL_STATE (readintime)
 
  integer(HID_T) :: file_id, dset_id, fspace_id, mspace_id, fapl_id, dxpl_id
  integer(HSIZE_T), dimension(3) :: xyzdims,count
  integer(HSSIZE_T), dimension(3) :: offset
  integer(HSSIZE_T), dimension(4) :: offset_pml
  integer(HSIZE_T), dimension(4) :: size_pml,size_mypml
  integer :: readintime
  integer :: pmlsize
  integer, dimension(3) :: dims
  character*(timestamp_size) tempc
  integer ierr

  call convert_to_string(readintime,tempc,timestamp_size)

  if (rank == 0) print *,'Reading in full snapshot from: ', trim(dir_out)//trim(simname)//'full_'//trim(tempc)//'.h5', ' iteration:',readintime,' and time', readintime*timestep/3600.0_dp, ' hours'

  call h5open_f(ierr)
  call h5pcreate_f(H5P_FILE_ACCESS_F, fapl_id, ierr)
  call h5pset_fapl_mpio_f(fapl_id, MPI_COMM_WORLD, MPI_INFO_NULL, ierr)
  call h5fopen_f(trim(dir_out)//trim(simname)//'full_'//trim(tempc)//'.h5', H5F_ACC_RDONLY_F, file_id, ierr, access_prp = fapl_id)
  call h5pcreate_f(H5P_DATASET_XFER_F, dxpl_id, ierr)
  call h5pset_dxpl_mpio_f(dxpl_id, H5FD_MPIO_COLLECTIVE_F, ierr)
  call h5ltget_attribute_int_f(file_id,"/","xyzdims", dims, ierr)

  count(1)  = dim1(rank) 
  count(2)  = dim2(rank) 
  count(3)  = dim3(rank) 
  offset(1) = xstart(rank)-1
  offset(2) = ystart(rank)-1
  offset(3) = zstart(rank)-1
  xyzdims = dims

  call h5screate_simple_f(3, count, mspace_id, ierr)
  call h5dopen_f(file_id,"vx",dset_id,ierr)
  call h5dget_space_f(dset_id, fspace_id, ierr)
  call h5sselect_hyperslab_f (fspace_id, H5S_SELECT_SET_F, offset, count, ierr)

  if (.not. DISPL) then

    call h5dread_f(dset_id, H5T_NATIVE_DOUBLE,a(:,:,:,2),xyzdims, ierr, &
file_space_id=fspace_id, mem_space_id = mspace_id, xfer_prp = dxpl_id)
    call h5dclose_f(dset_id, ierr)

    a(:,:,:,2) = a(:,:,:,2)/dimc

    call h5dopen_f(file_id, "rho", dset_id, ierr)
    call h5dread_f(dset_id, H5T_NATIVE_DOUBLE,a(:,:,:,1),xyzdims, ierr, &
file_space_id=fspace_id, mem_space_id = mspace_id, xfer_prp = dxpl_id)
    call h5dclose_f(dset_id, ierr)

    a(:,:,:,1) = a(:,:,:,1)/dimrho

    call h5dopen_f(file_id, "pre", dset_id, ierr)
    call h5dread_f(dset_id, H5T_NATIVE_DOUBLE,a(:,:,:,5),xyzdims, ierr, &
file_space_id=fspace_id, mem_space_id = mspace_id, xfer_prp = dxpl_id)
    call h5dclose_f(dset_id, ierr)

    a(:,:,:,5) = a(:,:,:,5)/dimrho/dimc**2

    call h5dopen_f(file_id, "vy", dset_id,ierr)
    call h5dread_f(dset_id, H5T_NATIVE_DOUBLE,a(:,:,:,3),xyzdims, ierr, &
file_space_id=fspace_id, mem_space_id = mspace_id, xfer_prp = dxpl_id)
    call h5dclose_f(dset_id, ierr)

    a(:,:,:,3) = a(:,:,:,3)/dimc

    call h5dopen_f(file_id, "vz", dset_id,ierr)
    call h5dread_f(dset_id, H5T_NATIVE_DOUBLE,a(:,:,:,4),xyzdims, ierr, &
file_space_id=fspace_id, mem_space_id = mspace_id, xfer_prp = dxpl_id)
    call h5dclose_f(dset_id, ierr)

    a(:,:,:,4) = a(:,:,:,4)/dimc


    if (MAGNETIC) then

      call h5dopen_f(file_id, "bx", dset_id,ierr)
      call h5dread_f(dset_id, H5T_NATIVE_DOUBLE,a(:,:,:,6),xyzdims, ierr, &
file_space_id=fspace_id, mem_space_id = mspace_id, xfer_prp = dxpl_id)
      call h5dclose_f(dset_id, ierr)

      a(:,:,:,6) = a(:,:,:,6)/dimb

      call h5dopen_f(file_id, "by", dset_id,ierr)
      call h5dread_f(dset_id, H5T_NATIVE_DOUBLE,a(:,:,:,7),xyzdims, ierr, &
file_space_id=fspace_id, mem_space_id = mspace_id, xfer_prp = dxpl_id)
      call h5dclose_f(dset_id, ierr)

      a(:,:,:,7) = a(:,:,:,7)/dimb

      call h5dopen_f(file_id, "bz", dset_id,ierr)
      call h5dread_f(dset_id, H5T_NATIVE_DOUBLE,a(:,:,:,8),xyzdims, ierr, &
file_space_id=fspace_id, mem_space_id = mspace_id, xfer_prp = dxpl_id)
      call h5dclose_f(dset_id, ierr)

      a(:,:,:,8) = a(:,:,:,8)/dimb

    endif
  else

    call h5dread_f(dset_id, H5T_NATIVE_DOUBLE,a(:,:,:,4),xyzdims, ierr, &
file_space_id=fspace_id, mem_space_id = mspace_id, xfer_prp = dxpl_id)
    call h5dclose_f(dset_id, ierr)

    a(:,:,:,4) = a(:,:,:,4)/dimc

    call h5dopen_f(file_id, "vy", dset_id,ierr)
    call h5dread_f(dset_id, H5T_NATIVE_DOUBLE,a(:,:,:,5),xyzdims, ierr, &
file_space_id=fspace_id, mem_space_id = mspace_id, xfer_prp = dxpl_id)
    call h5dclose_f(dset_id, ierr)

    a(:,:,:,5) = a(:,:,:,5)/dimc

    call h5dopen_f(file_id, "vz", dset_id,ierr)
    call h5dread_f(dset_id, H5T_NATIVE_DOUBLE,a(:,:,:,6),xyzdims, ierr, &
file_space_id=fspace_id, mem_space_id = mspace_id, xfer_prp = dxpl_id)
    call h5dclose_f(dset_id, ierr)

    a(:,:,:,6) = a(:,:,:,6)/dimc

    call h5dopen_f(file_id, "xix", dset_id,ierr)
    call h5dread_f(dset_id, H5T_NATIVE_DOUBLE,a(:,:,:,1),xyzdims, ierr, &
file_space_id=fspace_id, mem_space_id = mspace_id, xfer_prp = dxpl_id)
    call h5dclose_f(dset_id, ierr)

    a(:,:,:,1) = a(:,:,:,1)/diml

    call h5dopen_f(file_id, "xiy", dset_id,ierr)
    call h5dread_f(dset_id, H5T_NATIVE_DOUBLE,a(:,:,:,2),xyzdims, ierr, &
file_space_id=fspace_id, mem_space_id = mspace_id, xfer_prp = dxpl_id)
    call h5dclose_f(dset_id, ierr)

    a(:,:,:,2) = a(:,:,:,2)/diml

    call h5dopen_f(file_id, "xiz", dset_id,ierr)
    call h5dread_f(dset_id, H5T_NATIVE_DOUBLE,a(:,:,:,3),xyzdims, ierr, &
file_space_id=fspace_id, mem_space_id = mspace_id, xfer_prp = dxpl_id)
    call h5dclose_f(dset_id, ierr)

    a(:,:,:,3) = a(:,:,:,3)/diml

  endif

  call h5sclose_f(fspace_id, ierr)
  call h5sclose_f(mspace_id, ierr)
  call h5pclose_f(fapl_id, ierr)
  call h5pclose_f(dxpl_id, ierr)
  call h5fclose_f(file_id, ierr)


  if (MAGNETIC) then
    pmlsize = 6
  else 
    pmlsize = 2
  endif
  !! Top PML in Z

  if (USE_PML) then
    if (IAM_ZTOP) then
    offset_pml(1) = xstart(rank)-1
    offset_pml(2) = ystart(rank)-1
    offset_pml(3) = 0
    offset_pml(4) = 0
    size_mypml = (/dim1(rank),dim2(rank),nzpmltop,pmlsize/)
    size_pml = (/nx,ny,nzpmltop,pmlsize/)
    call h5pcreate_f(H5P_FILE_ACCESS_F, fapl_id, ierr)
    call h5pset_fapl_mpio_f(fapl_id, zsplit_comm, MPI_INFO_NULL, ierr)
    call h5fopen_f(trim(dir_out)//trim(simname)//'full_'//trim(tempc)//'.h5',H5F_ACC_RDONLY_F,file_id,ierr,access_prp = fapl_id)
    call h5screate_simple_f(4, size_mypml, mspace_id, ierr)
    call h5dopen_f(file_id, "PML_ZTOP",dset_id,ierr)
    call h5dget_space_f(dset_id, fspace_id, ierr)
    call h5sselect_hyperslab_f (fspace_id, H5S_SELECT_SET_F, offset_pml,size_mypml, ierr)
    call h5pcreate_f(H5P_DATASET_XFER_F, dxpl_id, ierr)
    call h5pset_dxpl_mpio_f(dxpl_id, H5FD_MPIO_COLLECTIVE_F, ierr)
    call h5dread_f(dset_id, H5T_NATIVE_DOUBLE,zpmlvars, size_pml, ierr, &
                 file_space_id = fspace_id, mem_space_id = mspace_id, xfer_prp =dxpl_id)
    call h5sclose_f(fspace_id,ierr)
    call h5sclose_f(mspace_id,ierr)
    call h5dclose_f(dset_id, ierr)
    call h5pclose_f(fapl_id, ierr)
    call h5pclose_f(dxpl_id, ierr)
    call h5fclose_f(file_id, ierr)
  endif

  call MPI_BARRIER(MPI_COMM_WORLD,ierr)

  if (IAM_ZBOT) then
    size_mypml = (/dim1(rank),dim2(rank),nzpmlbot,pmlsize/)
    size_pml = (/nx,ny,nzpmlbot,pmlsize/)
    offset_pml(1) = xstart(rank)-1
    offset_pml(2) = ystart(rank)-1
    offset_pml(3) = 0
    offset_pml(4) = 0
    call h5pcreate_f(H5P_FILE_ACCESS_F, fapl_id, ierr)
    call h5pset_fapl_mpio_f(fapl_id, zsplit_comm, MPI_INFO_NULL, ierr)
    call h5fopen_f(trim(dir_out)//trim(simname)//'full_'//trim(tempc)//'.h5',H5F_ACC_RDONLY_F,file_id,ierr,access_prp = fapl_id)
    call h5dopen_f(file_id, "PML_ZBOT",dset_id,ierr)
    call h5screate_simple_f(4, size_mypml, mspace_id, ierr)
    call h5dget_space_f(dset_id, fspace_id, ierr)
    call h5sselect_hyperslab_f (fspace_id, H5S_SELECT_SET_F, offset_pml,size_mypml, ierr)
    call h5pcreate_f(H5P_DATASET_XFER_F, dxpl_id, ierr)
    call h5pset_dxpl_mpio_f(dxpl_id, H5FD_MPIO_COLLECTIVE_F, ierr)
    call h5dread_f(dset_id, H5T_NATIVE_DOUBLE,zpmlvars, size_pml, ierr, &
                 file_space_id = fspace_id, mem_space_id = mspace_id, xfer_prp =dxpl_id)
    call h5sclose_f(fspace_id,ierr)
    call h5sclose_f(mspace_id,ierr)
    call h5dclose_f(dset_id, ierr)
    call h5pclose_f(fapl_id, ierr)
    call h5pclose_f(dxpl_id, ierr)
    call h5fclose_f(file_id, ierr)
    endif
    call MPI_BARRIER(MPI_COMM_WORLD,ierr)
  endif 

  if (USE_HPML) then
    if (IAM_XTOP) then
    size_mypml = (/nxpmltop,dim2(rank),dim3(rank),pmlsize/)
    size_pml = (/nxpmltop,ny,nz,pmlsize/)
    offset_pml(1) = 0
    offset_pml(2) = ystart(rank)-1
    offset_pml(3) = zstart(rank)-1
    offset_pml(4) = 0
    call h5pcreate_f(H5P_FILE_ACCESS_F, fapl_id, ierr)
    call h5pset_fapl_mpio_f(fapl_id, xsplit_comm, MPI_INFO_NULL, ierr)
    call h5fopen_f(trim(dir_out)//trim(simname)//'full_'//trim(tempc)//'.h5',H5F_ACC_RDONLY_F,file_id,ierr,access_prp = fapl_id)
    call h5dopen_f(file_id, "PML_XTOP",dset_id,ierr)
    call h5screate_simple_f(4, size_mypml, mspace_id, ierr)
    call h5dget_space_f(dset_id, fspace_id, ierr)
    call h5sselect_hyperslab_f (fspace_id, H5S_SELECT_SET_F, offset_pml,size_mypml, ierr)
    call h5pcreate_f(H5P_DATASET_XFER_F, dxpl_id, ierr)
    call h5pset_dxpl_mpio_f(dxpl_id, H5FD_MPIO_COLLECTIVE_F, ierr)
    call h5dread_f(dset_id, H5T_NATIVE_DOUBLE,xpmlvars, size_pml, ierr, &
                 file_space_id = fspace_id, mem_space_id = mspace_id, xfer_prp =dxpl_id)
    call h5sclose_f(fspace_id,ierr)
    call h5sclose_f(mspace_id,ierr)
    call h5dclose_f(dset_id, ierr)
    call h5pclose_f(fapl_id, ierr)
    call h5pclose_f(dxpl_id, ierr)
    call h5fclose_f(file_id, ierr)
    endif
    call MPI_BARRIER(MPI_COMM_WORLD,ierr)
    
    if (IAM_XBOT) then
    size_mypml = (/nxpmlbot,dim2(rank),dim3(rank),pmlsize/)
    size_pml = (/nxpmlbot,ny,nz,pmlsize/)
    offset_pml(1) = 0
    offset_pml(2) = ystart(rank)-1
    offset_pml(3) = zstart(rank)-1
    offset_pml(4) = 0 
    call h5pcreate_f(H5P_FILE_ACCESS_F, fapl_id, ierr)
    call h5pset_fapl_mpio_f(fapl_id, xsplit_comm, MPI_INFO_NULL, ierr)
    call h5fopen_f(trim(dir_out)//trim(simname)//'full_'//trim(tempc)//'.h5',H5F_ACC_RDONLY_F,file_id,ierr,access_prp = fapl_id)
    call h5dopen_f(file_id, "PML_XBOT",dset_id,ierr)
    call h5screate_simple_f(4, size_mypml, mspace_id, ierr)
    call h5dget_space_f(dset_id, fspace_id, ierr)
    call h5sselect_hyperslab_f (fspace_id, H5S_SELECT_SET_F, offset_pml,size_mypml, ierr)
    call h5pcreate_f(H5P_DATASET_XFER_F, dxpl_id, ierr)
    call h5pset_dxpl_mpio_f(dxpl_id, H5FD_MPIO_COLLECTIVE_F, ierr)
    call h5dread_f(dset_id, H5T_NATIVE_DOUBLE,xpmlvars, size_pml, ierr, &
                 file_space_id = fspace_id, mem_space_id = mspace_id, xfer_prp =dxpl_id)
    call h5sclose_f(fspace_id,ierr)
    call h5sclose_f(mspace_id,ierr)
    call h5dclose_f(dset_id, ierr)
    call h5pclose_f(fapl_id, ierr)
    call h5pclose_f(dxpl_id, ierr)
    call h5fclose_f(file_id, ierr)
   endif
    call MPI_BARRIER(MPI_COMM_WORLD,ierr) 

    if (IAM_YTOP) then
    size_mypml = (/dim1(rank),nypmlbot,dim3(rank),pmlsize/)
    size_pml = (/nx,nypmlbot,nz,pmlsize/)
    offset_pml(1) = xstart(rank)-1
    offset_pml(2) = 0
    offset_pml(3) = zstart(rank)-1
    offset_pml(4) = 0 
    
    call h5pcreate_f(H5P_FILE_ACCESS_F, fapl_id, ierr)
    call h5pset_fapl_mpio_f(fapl_id, ysplit_comm, MPI_INFO_NULL, ierr)
    call h5fopen_f(trim(dir_out)//trim(simname)//'full_'//trim(tempc)//'.h5',H5F_ACC_RDONLY_F,file_id,ierr,access_prp = fapl_id)
    call h5dopen_f(file_id, "PML_YTOP",dset_id,ierr)
    call h5screate_simple_f(4, size_mypml, mspace_id, ierr)
    call h5dget_space_f(dset_id, fspace_id, ierr)
    call h5sselect_hyperslab_f (fspace_id, H5S_SELECT_SET_F, offset_pml,size_mypml, ierr)
    call h5pcreate_f(H5P_DATASET_XFER_F, dxpl_id, ierr)
    call h5pset_dxpl_mpio_f(dxpl_id, H5FD_MPIO_COLLECTIVE_F, ierr)
    call h5dread_f(dset_id, H5T_NATIVE_DOUBLE,ypmlvars, size_pml, ierr, &
                 file_space_id = fspace_id, mem_space_id = mspace_id, xfer_prp =dxpl_id)
    call h5sclose_f(fspace_id,ierr)
    call h5sclose_f(mspace_id,ierr)
    call h5dclose_f(dset_id, ierr)
    call h5pclose_f(fapl_id, ierr)
    call h5pclose_f(dxpl_id, ierr)
    call h5fclose_f(file_id, ierr)
    endif
    call MPI_BARRIER(MPI_COMM_WORLD,ierr)
    
    if (IAM_YBOT) then
    size_mypml = (/dim1(rank),nypmlbot,dim3(rank),pmlsize/)
    size_pml = (/nx,nypmlbot,nz,pmlsize/)
    offset_pml(1) = xstart(rank)-1
    offset_pml(2) = 0
    offset_pml(3) = zstart(rank)-1
    offset_pml(4) = 0
    call h5pcreate_f(H5P_FILE_ACCESS_F, fapl_id, ierr)
    call h5pset_fapl_mpio_f(fapl_id, ysplit_comm, MPI_INFO_NULL, ierr)
    call h5fopen_f(trim(dir_out)//trim(simname)//'full_'//trim(tempc)//'.h5',H5F_ACC_RDONLY_F,file_id,ierr,access_prp = fapl_id)
    call h5dopen_f(file_id, "PML_YBOT",dset_id,ierr)
    call h5screate_simple_f(4, size_mypml, mspace_id, ierr)
    call h5dget_space_f(dset_id, fspace_id, ierr)
    call h5sselect_hyperslab_f (fspace_id, H5S_SELECT_SET_F, offset_pml,size_mypml, ierr)
    call h5pcreate_f(H5P_DATASET_XFER_F, dxpl_id, ierr)
    call h5pset_dxpl_mpio_f(dxpl_id, H5FD_MPIO_COLLECTIVE_F, ierr)
    call h5dread_f(dset_id, H5T_NATIVE_DOUBLE, ypmlvars, size_pml, ierr, &
                 file_space_id = fspace_id, mem_space_id = mspace_id, xfer_prp =dxpl_id)
    call h5sclose_f(fspace_id,ierr)
    call h5sclose_f(mspace_id,ierr)
    call h5dclose_f(dset_id, ierr)
    call h5pclose_f(fapl_id, ierr)
    call h5pclose_f(dxpl_id, ierr)
    call h5fclose_f(file_id, ierr)
    endif
    call MPI_BARRIER(MPI_COMM_WORLD,ierr)
  endif

  call MPI_BARRIER(comm_cart_3D, ierr)

  end subroutine READ_IN_FULL_STATE
!!================================================================================

subroutine WRITE_OUT_FULL_STATE(fail)

  integer(HID_T) :: file_id, dset_id, fspace_id, mspace_id, fapl_id, dxpl_id

  integer(HSIZE_T), dimension(1) :: zdims,zcount
  integer(HSIZE_T), dimension(3) :: xyzdims,count
  integer(HSIZE_T), dimension(4) :: size_pml,size_mypml
  integer(HSSIZE_T), dimension(3) :: offset
  integer(HSSIZE_T), dimension(4) :: offset_pml
  integer(SIZE_T) :: size = 1

  logical,intent(in) :: fail
  real(dp), dimension(nz) :: grav
  character*(timestamp_size) tempc
  integer ierr
  REAL(dp), dimension(1) :: timeout
  integer outtime,pmlsize

  if (time .EQ. 0) then
   outtime = 0
  else
   outtime = time/fulloutputsteps
  endif

  if (fail) then
    tempc = 'XXXX'
  else
   call convert_to_string(outtime,tempc,timestamp_size)
  endif

  if (rank == 0) print *,'Writing out full snapshot from: ', trim(dir_out)//trim(simname)//'full_'//trim(tempc)//'.h5', ' iteration:',time,' and time', time*timestep/3600.0_dp, ' hours'


  timeout = time*timestep/3600.0_dp

  call h5open_f(ierr)

  call h5pcreate_f(H5P_FILE_ACCESS_F, fapl_id, ierr)
  call h5pset_fapl_mpio_f(fapl_id,MPI_COMM_WORLD,MPI_INFO_NULL,ierr)

  call h5fcreate_f(trim(dir_out)//trim(simname)//'full_'//trim(tempc)//'.h5', &
H5F_ACC_TRUNC_F, file_id, ierr, access_prp = fapl_id)
  call h5pcreate_f(H5P_DATASET_XFER_F, dxpl_id, ierr)
  call h5pset_dxpl_mpio_f(dxpl_id, H5FD_MPIO_COLLECTIVE_F, ierr)

  call h5ltset_attribute_int_f(file_id,"/","xyzdims", &
(/nx,ny,nz/),3*size, ierr)
  call h5ltset_attribute_double_f(file_id,"/","time", &
timeout, size, ierr)
    
  count(1)  = dim1(rank) 
  count(2)  = dim2(rank) 
  count(3)  = dim3(rank) 
  offset(1) = xstart(rank)-1
  offset(2) = ystart(rank)-1
  offset(3) = zstart(rank)-1
  xyzdims = (/nx,ny,nz/)

  call h5screate_simple_f(3, xyzdims, fspace_id, ierr)
  call h5screate_simple_f(3, count, mspace_id, ierr)
  call h5sselect_hyperslab_f(fspace_id, H5S_SELECT_SET_F, offset, count, &
ierr)

  if (.not. DISPL) then

    call h5dcreate_f(file_id, "pre", H5T_NATIVE_DOUBLE, fspace_id, dset_id, &
ierr)
    call h5dwrite_f(dset_id, H5T_NATIVE_DOUBLE,a(:,:,:,5)*dimrho*dimc**2,xyzdims, ierr, &
           file_space_id=fspace_id, mem_space_id = mspace_id, xfer_prp = dxpl_id)
    call h5dclose_f(dset_id, ierr)

    call h5dcreate_f(file_id, "rho", H5T_NATIVE_DOUBLE, fspace_id, dset_id, &
ierr)
    call h5dwrite_f(dset_id, H5T_NATIVE_DOUBLE,a(:,:,:,1)*dimrho,xyzdims, ierr, &
           file_space_id=fspace_id, mem_space_id = mspace_id, xfer_prp = dxpl_id)
    call h5dclose_f(dset_id, ierr)

    call h5dcreate_f(file_id, "vx", H5T_NATIVE_DOUBLE, fspace_id, dset_id, &
ierr)
    call h5dwrite_f(dset_id, H5T_NATIVE_DOUBLE,a(:,:,:,2)*dimc,xyzdims, ierr, &
           file_space_id=fspace_id, mem_space_id = mspace_id, xfer_prp = dxpl_id)
    call h5dclose_f(dset_id, ierr)

    call h5dcreate_f(file_id, "vy", H5T_NATIVE_DOUBLE, fspace_id, dset_id, &
ierr)
    call h5dwrite_f(dset_id, H5T_NATIVE_DOUBLE,a(:,:,:,3)*dimc,xyzdims, ierr, &
           file_space_id=fspace_id, mem_space_id = mspace_id, xfer_prp = dxpl_id)
    call h5dclose_f(dset_id, ierr)

    call h5dcreate_f(file_id, "vz", H5T_NATIVE_DOUBLE, fspace_id, dset_id, &
ierr)
    call h5dwrite_f(dset_id, H5T_NATIVE_DOUBLE,a(:,:,:,4)*dimc,xyzdims, ierr, &
           file_space_id=fspace_id, mem_space_id = mspace_id, xfer_prp = dxpl_id)
    call h5dclose_f(dset_id, ierr)

    if (MAGNETIC) then

      call h5dcreate_f(file_id, "bx", H5T_NATIVE_DOUBLE, fspace_id, dset_id, &
ierr)
      call h5dwrite_f(dset_id, H5T_NATIVE_DOUBLE,a(:,:,:,6)*dimb,xyzdims, ierr, &
file_space_id=fspace_id, mem_space_id = mspace_id, xfer_prp = dxpl_id)
      call h5dclose_f(dset_id, ierr)

      call h5dcreate_f(file_id, "by", H5T_NATIVE_DOUBLE, fspace_id, dset_id, &
ierr)
      call h5dwrite_f(dset_id, H5T_NATIVE_DOUBLE,a(:,:,:,7)*dimb,xyzdims, ierr, &
file_space_id=fspace_id, mem_space_id = mspace_id, xfer_prp = dxpl_id)
      call h5dclose_f(dset_id, ierr)

      call h5dcreate_f(file_id, "bz", H5T_NATIVE_DOUBLE, fspace_id, dset_id, &
ierr)
      call h5dwrite_f(dset_id, H5T_NATIVE_DOUBLE,a(:,:,:,8)*dimb,xyzdims, ierr, &
file_space_id=fspace_id, mem_space_id = mspace_id, xfer_prp = dxpl_id)
      call h5dclose_f(dset_id, ierr)
      
    endif
  else

    call h5dcreate_f(file_id, "vx", H5T_NATIVE_DOUBLE, fspace_id, dset_id, &
ierr)
    call h5dwrite_f(dset_id, H5T_NATIVE_DOUBLE,a(:,:,:,4)*dimc,xyzdims, ierr, &
file_space_id=fspace_id, mem_space_id = mspace_id, xfer_prp = dxpl_id)
    call h5dclose_f(dset_id, ierr)

    call h5dcreate_f(file_id, "vy", H5T_NATIVE_DOUBLE, fspace_id, dset_id, &
ierr)
    call h5dwrite_f(dset_id, H5T_NATIVE_DOUBLE,a(:,:,:,5)*dimc,xyzdims, ierr, &
file_space_id=fspace_id, mem_space_id = mspace_id, xfer_prp = dxpl_id)
    call h5dclose_f(dset_id, ierr)

    call h5dcreate_f(file_id, "vz", H5T_NATIVE_DOUBLE, fspace_id, dset_id, &
ierr)
    call h5dwrite_f(dset_id, H5T_NATIVE_DOUBLE,a(:,:,:,6)*dimc,xyzdims, ierr, &
file_space_id=fspace_id, mem_space_id = mspace_id, xfer_prp = dxpl_id)
    call h5dclose_f(dset_id, ierr)

    call h5dcreate_f(file_id, "xix", H5T_NATIVE_DOUBLE, fspace_id, dset_id, &
ierr)
    call h5dwrite_f(dset_id, H5T_NATIVE_DOUBLE,a(:,:,:,1)*diml,xyzdims, ierr, &
file_space_id=fspace_id, mem_space_id = mspace_id, xfer_prp = dxpl_id)
    call h5dclose_f(dset_id, ierr)

    call h5dcreate_f(file_id, "xiy", H5T_NATIVE_DOUBLE, fspace_id, dset_id, &
ierr)
    call h5dwrite_f(dset_id, H5T_NATIVE_DOUBLE,a(:,:,:,2)*diml,xyzdims, ierr, &
file_space_id=fspace_id, mem_space_id = mspace_id, xfer_prp = dxpl_id)
    call h5dclose_f(dset_id, ierr)

    call h5dcreate_f(file_id, "xiz", H5T_NATIVE_DOUBLE, fspace_id, dset_id, &
ierr)
    call h5dwrite_f(dset_id, H5T_NATIVE_DOUBLE,a(:,:,:,3)*diml,xyzdims, ierr, &
file_space_id=fspace_id, mem_space_id = mspace_id, xfer_prp = dxpl_id)
    call h5dclose_f(dset_id, ierr)

  endif

      ! Close filespace and memspace from non-PML data
  call h5sclose_f(fspace_id, ierr)
  call h5sclose_f(mspace_id, ierr)

      ! Close property lists and file
  call h5pclose_f(fapl_id, ierr)   ! Collective access to the file
  call h5pclose_f(dxpl_id, ierr)   ! Collective access to the dataset
  call h5fclose_f(file_id, ierr)   ! The file


  if (MAGNETIC) then
    pmlsize = 6
  else 
    pmlsize = 2
  endif
  !! Top PML in Z

  if (USE_PML) then
    if (IAM_ZTOP) then
    offset_pml(1) = xstart(rank)-1
    offset_pml(2) = ystart(rank)-1
    offset_pml(3) = 0
    offset_pml(4) = 0
    size_mypml = (/dim1(rank),dim2(rank),nzpmltop,pmlsize/)
    size_pml = (/nx,ny,nzpmltop,pmlsize/)
    call h5pcreate_f(H5P_FILE_ACCESS_F, fapl_id, ierr)
    call h5pset_fapl_mpio_f(fapl_id, zsplit_comm, MPI_INFO_NULL, ierr)
    call h5fopen_f(trim(dir_out)//trim(simname)//'full_'//trim(tempc)//'.h5',H5F_ACC_RDWR_F,file_id,ierr,access_prp = fapl_id)
    call h5screate_simple_f(4, size_pml, fspace_id, ierr)
    call h5dcreate_f(file_id, "PML_ZTOP",   H5T_NATIVE_DOUBLE, fspace_id, dset_id,ierr)
    call h5sclose_f(fspace_id, ierr)
    call h5screate_simple_f(4, size_mypml, mspace_id, ierr)
    call h5dget_space_f(dset_id, fspace_id, ierr)
    call h5sselect_hyperslab_f (fspace_id, H5S_SELECT_SET_F, offset_pml,size_mypml, ierr)
    call h5pcreate_f(H5P_DATASET_XFER_F, dxpl_id, ierr)
    call h5pset_dxpl_mpio_f(dxpl_id, H5FD_MPIO_COLLECTIVE_F, ierr)
    call h5dwrite_f(dset_id, H5T_NATIVE_DOUBLE,zpmlvars, size_pml, ierr, &
                 file_space_id = fspace_id, mem_space_id = mspace_id, xfer_prp =dxpl_id)
    call h5sclose_f(fspace_id,ierr)
    call h5sclose_f(mspace_id,ierr)
    call h5dclose_f(dset_id, ierr)
    call h5pclose_f(fapl_id, ierr)
    call h5pclose_f(dxpl_id, ierr)
    call h5fclose_f(file_id, ierr)
  endif

  call MPI_BARRIER(MPI_COMM_WORLD,ierr)

  if (IAM_ZBOT) then
    size_mypml = (/dim1(rank),dim2(rank),nzpmlbot,pmlsize/)
    size_pml = (/nx,ny,nzpmlbot,pmlsize/)
    offset_pml(1) = xstart(rank)-1
    offset_pml(2) = ystart(rank)-1
    offset_pml(3) = 0
    offset_pml(4) = 0
    call h5pcreate_f(H5P_FILE_ACCESS_F, fapl_id, ierr)
    call h5pset_fapl_mpio_f(fapl_id, zsplit_comm, MPI_INFO_NULL, ierr)
    call h5fopen_f(trim(dir_out)//trim(simname)//'full_'//trim(tempc)//'.h5',H5F_ACC_RDWR_F,file_id,ierr,access_prp = fapl_id)
    call h5screate_simple_f(4, size_pml, fspace_id, ierr)
    call h5dcreate_f(file_id, "PML_ZBOT",   H5T_NATIVE_DOUBLE, fspace_id, dset_id,ierr)
    call h5sclose_f(fspace_id, ierr)
    call h5screate_simple_f(4, size_mypml, mspace_id, ierr)
    call h5dget_space_f(dset_id, fspace_id, ierr)
    call h5sselect_hyperslab_f (fspace_id, H5S_SELECT_SET_F, offset_pml,size_mypml, ierr)
    call h5pcreate_f(H5P_DATASET_XFER_F, dxpl_id, ierr)
    call h5pset_dxpl_mpio_f(dxpl_id, H5FD_MPIO_COLLECTIVE_F, ierr)
    call h5dwrite_f(dset_id, H5T_NATIVE_DOUBLE,zpmlvars, size_pml, ierr, &
                 file_space_id = fspace_id, mem_space_id = mspace_id, xfer_prp =dxpl_id)
    call h5sclose_f(fspace_id,ierr)
    call h5sclose_f(mspace_id,ierr)
    call h5dclose_f(dset_id, ierr)
    call h5pclose_f(fapl_id, ierr)
    call h5pclose_f(dxpl_id, ierr)
    call h5fclose_f(file_id, ierr)
    endif
    call MPI_BARRIER(MPI_COMM_WORLD,ierr)
  endif 

  if (USE_HPML) then
    if (IAM_XTOP) then
    size_mypml = (/nxpmltop,dim2(rank),dim3(rank),pmlsize/)
    size_pml = (/nxpmltop,ny,nz,pmlsize/)
    offset_pml(1) = 0
    offset_pml(2) = ystart(rank)-1
    offset_pml(3) = zstart(rank)-1
    offset_pml(4) = 0
    call h5pcreate_f(H5P_FILE_ACCESS_F, fapl_id, ierr)
    call h5pset_fapl_mpio_f(fapl_id, xsplit_comm, MPI_INFO_NULL, ierr)
    call h5fopen_f(trim(dir_out)//trim(simname)//'full_'//trim(tempc)//'.h5',H5F_ACC_RDWR_F,file_id,ierr,access_prp = fapl_id)
    call h5screate_simple_f(4, size_pml, fspace_id, ierr)
    call h5dcreate_f(file_id, "PML_XTOP",   H5T_NATIVE_DOUBLE, fspace_id,dset_id,ierr)
    call h5sclose_f(fspace_id, ierr)
    call h5screate_simple_f(4, size_mypml, mspace_id, ierr)
    call h5dget_space_f(dset_id, fspace_id, ierr)
    call h5sselect_hyperslab_f (fspace_id, H5S_SELECT_SET_F, offset_pml,size_mypml, ierr)
    call h5pcreate_f(H5P_DATASET_XFER_F, dxpl_id, ierr)
    call h5pset_dxpl_mpio_f(dxpl_id, H5FD_MPIO_COLLECTIVE_F, ierr)
    call h5dwrite_f(dset_id, H5T_NATIVE_DOUBLE,xpmlvars, size_pml, ierr, &
                 file_space_id = fspace_id, mem_space_id = mspace_id, xfer_prp =dxpl_id)
    call h5sclose_f(fspace_id,ierr)
    call h5sclose_f(mspace_id,ierr)
    call h5dclose_f(dset_id, ierr)
    call h5pclose_f(fapl_id, ierr)
    call h5pclose_f(dxpl_id, ierr)
    call h5fclose_f(file_id, ierr)
    endif
    call MPI_BARRIER(MPI_COMM_WORLD,ierr)
    
    if (IAM_XBOT) then
    size_mypml = (/nxpmlbot,dim2(rank),dim3(rank),pmlsize/)
    size_pml = (/nxpmlbot,ny,nz,pmlsize/)
    offset_pml(1) = 0
    offset_pml(2) = ystart(rank)-1
    offset_pml(3) = zstart(rank)-1
    offset_pml(4) = 0 
    call h5pcreate_f(H5P_FILE_ACCESS_F, fapl_id, ierr)
    call h5pset_fapl_mpio_f(fapl_id, xsplit_comm, MPI_INFO_NULL, ierr)
    call h5fopen_f(trim(dir_out)//trim(simname)//'full_'//trim(tempc)//'.h5',H5F_ACC_RDWR_F,file_id,ierr,access_prp = fapl_id)
    call h5screate_simple_f(4, size_pml, fspace_id, ierr)
    call h5dcreate_f(file_id, "PML_XBOT",   H5T_NATIVE_DOUBLE, fspace_id,dset_id,ierr)
    call h5sclose_f(fspace_id, ierr)
    call h5screate_simple_f(4, size_mypml, mspace_id, ierr)
    call h5dget_space_f(dset_id, fspace_id, ierr)
    call h5sselect_hyperslab_f (fspace_id, H5S_SELECT_SET_F, offset_pml,size_mypml, ierr)
    call h5pcreate_f(H5P_DATASET_XFER_F, dxpl_id, ierr)
    call h5pset_dxpl_mpio_f(dxpl_id, H5FD_MPIO_COLLECTIVE_F, ierr)
    call h5dwrite_f(dset_id, H5T_NATIVE_DOUBLE,xpmlvars, size_pml, ierr, &
                 file_space_id = fspace_id, mem_space_id = mspace_id, xfer_prp =dxpl_id)
    call h5sclose_f(fspace_id,ierr)
    call h5sclose_f(mspace_id,ierr)
    call h5dclose_f(dset_id, ierr)
    call h5pclose_f(fapl_id, ierr)
    call h5pclose_f(dxpl_id, ierr)
    call h5fclose_f(file_id, ierr)
   endif
    call MPI_BARRIER(MPI_COMM_WORLD,ierr) 

    if (IAM_YTOP) then
    size_mypml = (/dim1(rank),nypmlbot,dim3(rank),pmlsize/)
    size_pml = (/nx,nypmlbot,nz,pmlsize/)
    offset_pml(1) = xstart(rank)-1
    offset_pml(2) = 0
    offset_pml(3) = zstart(rank)-1
    offset_pml(4) = 0 
    
    call h5pcreate_f(H5P_FILE_ACCESS_F, fapl_id, ierr)
    call h5pset_fapl_mpio_f(fapl_id, ysplit_comm, MPI_INFO_NULL, ierr)
    call h5fopen_f(trim(dir_out)//trim(simname)//'full_'//trim(tempc)//'.h5',H5F_ACC_RDWR_F,file_id,ierr,access_prp = fapl_id)
    call h5screate_simple_f(4, size_pml, fspace_id, ierr)
    call h5dcreate_f(file_id, "PML_YTOP",   H5T_NATIVE_DOUBLE, fspace_id,dset_id,ierr)
    call h5sclose_f(fspace_id, ierr)
    call h5screate_simple_f(4, size_mypml, mspace_id, ierr)
    call h5dget_space_f(dset_id, fspace_id, ierr)
    call h5sselect_hyperslab_f (fspace_id, H5S_SELECT_SET_F, offset_pml,size_mypml, ierr)
    call h5pcreate_f(H5P_DATASET_XFER_F, dxpl_id, ierr)
    call h5pset_dxpl_mpio_f(dxpl_id, H5FD_MPIO_COLLECTIVE_F, ierr)
    call h5dwrite_f(dset_id, H5T_NATIVE_DOUBLE,ypmlvars, size_pml, ierr, &
                 file_space_id = fspace_id, mem_space_id = mspace_id, xfer_prp =dxpl_id)
    call h5sclose_f(fspace_id,ierr)
    call h5sclose_f(mspace_id,ierr)
    call h5dclose_f(dset_id, ierr)
    call h5pclose_f(fapl_id, ierr)
    call h5pclose_f(dxpl_id, ierr)
    call h5fclose_f(file_id, ierr)
    endif
    call MPI_BARRIER(MPI_COMM_WORLD,ierr)
    
    if (IAM_YBOT) then
    size_mypml = (/dim1(rank),nypmlbot,dim3(rank),pmlsize/)
    size_pml = (/nx,nypmlbot,nz,pmlsize/)
    offset_pml(1) = xstart(rank)-1
    offset_pml(2) = 0
    offset_pml(3) = zstart(rank)-1
    offset_pml(4) = 0
    call h5pcreate_f(H5P_FILE_ACCESS_F, fapl_id, ierr)
    call h5pset_fapl_mpio_f(fapl_id, ysplit_comm, MPI_INFO_NULL, ierr)
    call h5fopen_f(trim(dir_out)//trim(simname)//'full_'//trim(tempc)//'.h5',H5F_ACC_RDWR_F,file_id,ierr,access_prp = fapl_id)
    call h5screate_simple_f(4, size_pml, fspace_id, ierr)
    call h5dcreate_f(file_id, "PML_YBOT",   H5T_NATIVE_DOUBLE, fspace_id,dset_id,ierr)
    call h5sclose_f(fspace_id, ierr)
    call h5screate_simple_f(4, size_mypml, mspace_id, ierr)
    call h5dget_space_f(dset_id, fspace_id, ierr)
    call h5sselect_hyperslab_f (fspace_id, H5S_SELECT_SET_F, offset_pml,size_mypml, ierr)
    call h5pcreate_f(H5P_DATASET_XFER_F, dxpl_id, ierr)
    call h5pset_dxpl_mpio_f(dxpl_id, H5FD_MPIO_COLLECTIVE_F, ierr)
    call h5dwrite_f(dset_id, H5T_NATIVE_DOUBLE, ypmlvars, size_pml, ierr, &
                 file_space_id = fspace_id, mem_space_id = mspace_id, xfer_prp =dxpl_id)
    call h5sclose_f(fspace_id,ierr)
    call h5sclose_f(mspace_id,ierr)
    call h5dclose_f(dset_id, ierr)
    call h5pclose_f(fapl_id, ierr)
    call h5pclose_f(dxpl_id, ierr)
    call h5fclose_f(file_id, ierr)
    endif
    call MPI_BARRIER(MPI_COMM_WORLD,ierr)
  endif
  !!! If time == 0 then save backgeround

  if (time == 0) then

    call h5open_f(ierr)

    call h5pcreate_f(H5P_FILE_ACCESS_F, fapl_id, ierr)
    call h5pset_fapl_mpio_f(fapl_id,MPI_COMM_WORLD,MPI_INFO_NULL,ierr)

    call h5fcreate_f(trim(dir_out)//trim(simname)//'bkg.h5', &
H5F_ACC_TRUNC_F, file_id, ierr, access_prp = fapl_id)
    call h5pcreate_f(H5P_DATASET_XFER_F, dxpl_id, ierr)
    call h5pset_dxpl_mpio_f(dxpl_id, H5FD_MPIO_COLLECTIVE_F, ierr)

    call h5ltset_attribute_int_f(file_id,"/","xyzdims", &
(/nx,ny,nz/),3*size, ierr)
    call h5ltset_attribute_double_f(file_id,"/","time", &
timeout, size, ierr)

    count(1)  = dim1(rank) 
    count(2)  = dim2(rank) 
    count(3)  = dim3(rank) 
    offset(1) = xstart(rank)-1
    offset(2) = ystart(rank)-1
    offset(3) = zstart(rank)-1
    xyzdims = (/nx,ny,nz/)

    call h5screate_simple_f(3, xyzdims, fspace_id, ierr)
    call h5screate_simple_f(3, count, mspace_id, ierr)
    call h5sselect_hyperslab_f(fspace_id, H5S_SELECT_SET_F, offset, count, &
ierr)

    call h5dcreate_f(file_id, "pre0", H5T_NATIVE_DOUBLE, fspace_id, dset_id, &
ierr)
    call h5dwrite_f(dset_id, H5T_NATIVE_DOUBLE,p0(:,:,:)*dimrho*dimc**2,xyzdims, ierr, &
file_space_id=fspace_id, mem_space_id = mspace_id, xfer_prp = dxpl_id)
    call h5dclose_f(dset_id, ierr)

    call h5dcreate_f(file_id, "rho0", H5T_NATIVE_DOUBLE, fspace_id, dset_id, &
ierr)
    call h5dwrite_f(dset_id, H5T_NATIVE_DOUBLE,rho0(:,:,:)*dimrho,xyzdims, ierr, &
file_space_id=fspace_id, mem_space_id = mspace_id, xfer_prp = dxpl_id)
    call h5dclose_f(dset_id, ierr)

    call h5dcreate_f(file_id, "sponge", H5T_NATIVE_DOUBLE, fspace_id, dset_id, &
ierr)
    call h5dwrite_f(dset_id, H5T_NATIVE_DOUBLE,spongexyz(:,:,:),xyzdims, ierr, &
file_space_id=fspace_id, mem_space_id = mspace_id, xfer_prp = dxpl_id)
    call h5dclose_f(dset_id, ierr)

  if (MAGNETIC) then 
    call h5dcreate_f(file_id, "reduction", H5T_NATIVE_DOUBLE, fspace_id, dset_id, &
ierr)
    call h5dwrite_f(dset_id, H5T_NATIVE_DOUBLE,reduction(:,:,:),xyzdims, ierr, &
file_space_id=fspace_id, mem_space_id = mspace_id, xfer_prp = dxpl_id)
    call h5dclose_f(dset_id, ierr)

      call h5dcreate_f(file_id, "bx0", H5T_NATIVE_DOUBLE, fspace_id, dset_id, &
ierr)
      call h5dwrite_f(dset_id, H5T_NATIVE_DOUBLE,box*dimb,xyzdims, ierr,&
file_space_id=fspace_id, mem_space_id = mspace_id, xfer_prp = dxpl_id)
      call h5dclose_f(dset_id, ierr)

      call h5dcreate_f(file_id, "by0", H5T_NATIVE_DOUBLE, fspace_id, dset_id, &
ierr)
      call h5dwrite_f(dset_id, H5T_NATIVE_DOUBLE,boy*dimb,xyzdims, ierr,&
file_space_id=fspace_id, mem_space_id = mspace_id, xfer_prp = dxpl_id)
      call h5dclose_f(dset_id, ierr)

      call h5dcreate_f(file_id, "bz0", H5T_NATIVE_DOUBLE, fspace_id, dset_id, &
ierr)
      call h5dwrite_f(dset_id, H5T_NATIVE_DOUBLE,boz*dimb,xyzdims, ierr,&
file_space_id=fspace_id, mem_space_id = mspace_id, xfer_prp = dxpl_id)
      call h5dclose_f(dset_id, ierr)

    endif

    if (FLOWS) then

      call h5dcreate_f(file_id, "vx0", H5T_NATIVE_DOUBLE, fspace_id, dset_id, &
ierr)
      call h5dwrite_f(dset_id, H5T_NATIVE_DOUBLE,v0_x*dimc,xyzdims, ierr,&
file_space_id=fspace_id, mem_space_id = mspace_id, xfer_prp = dxpl_id)
      call h5dclose_f(dset_id, ierr)

      call h5dcreate_f(file_id, "vy0", H5T_NATIVE_DOUBLE, fspace_id, dset_id, &
ierr)
      call h5dwrite_f(dset_id, H5T_NATIVE_DOUBLE,v0_y*dimc,xyzdims, ierr,&
file_space_id=fspace_id, mem_space_id = mspace_id, xfer_prp = dxpl_id)
      call h5dclose_f(dset_id, ierr)

      call h5dcreate_f(file_id, "vz0", H5T_NATIVE_DOUBLE, fspace_id, dset_id, &
ierr)
      call h5dwrite_f(dset_id, H5T_NATIVE_DOUBLE,v0_z*dimc,xyzdims, ierr,&
file_space_id=fspace_id, mem_space_id = mspace_id, xfer_prp = dxpl_id)
      call h5dclose_f(dset_id, ierr)

    endif

    ! Close filespace and memspace from non-PML data
    call h5sclose_f(fspace_id, ierr)
    call h5sclose_f(mspace_id, ierr)

    ! Close property lists and file
    call h5pclose_f(fapl_id, ierr)   ! Collective access to the file
    call h5pclose_f(dxpl_id, ierr)   ! Collective access to the dataset
    call h5fclose_f(file_id, ierr)   ! The file

    call MPI_BARRIER(MPI_COMM_WORLD,ierr)

    zdims = (/nz/)
    zcount = (/nz/) 

    if (rank == 0) then

      call h5screate_simple_f(1, zdims, fspace_id, ierr)
      call h5screate_simple_f(1, zcount, mspace_id, ierr)

      call h5fopen_f(trim(dir_out)//trim(simname)//'bkg.h5', H5F_ACC_RDWR_F, file_id, ierr)
      
      call h5dcreate_f(file_id,"zaxis",H5T_NATIVE_DOUBLE,fspace_id,dset_id,ierr)
      call h5dwrite_f(dset_id, H5T_NATIVE_DOUBLE,(zaxis-1.0_dp)*rsun/1.0e8_dp,zdims, ierr)
      call h5dclose_f(dset_id,ierr)

      call h5dcreate_f(file_id,"grav",H5T_NATIVE_DOUBLE,fspace_id,dset_id,ierr)
      call h5dwrite_f(dset_id, H5T_NATIVE_DOUBLE,grav,zdims, ierr)
      call h5dclose_f(dset_id,ierr)

      call h5sclose_f(fspace_id, ierr)
      call h5sclose_f(mspace_id, ierr)

      call h5fclose_f(file_id,ierr)

    endif

  endif
  call MPI_BARRIER(comm_cart_3D, ierr)

end subroutine WRITE_OUT_FULL_STATE

!!================================================================================

subroutine WRITE_OUT_CHUNK()

  integer(HID_T) :: file_id, dset_id, fspace_id, mspace_id, fapl_id, dxpl_id

  integer(HSIZE_T), dimension(1) :: xdims,ydims,zdims
  integer(HSIZE_T), dimension(3) :: size_chunk,size_mychunk
  integer(HSSIZE_T), dimension(3) :: offset
  integer(SIZE_T) :: size = 1

  character*(timestamp_size) tempc
  integer ierr
  REAL(dp), dimension(1) :: timeout
  integer outtime

  if (minsavtime .GT. 0) then
   outtime = time/chunkoutputsteps - floor(minsavtime*3600/timestep+1.0d-5)/chunkoutputsteps
  else
   outtime = time/chunkoutputsteps
  endif

  call convert_to_string(outtime,tempc,timestamp_size)

  if (rank == 0) print *,'Writing out chunk snapshot from: ', trim(dir_out)//trim(simname)//'chunk_'//trim(tempc)//'.h5', ' iteration:',time,' and time', time*timestep/3600.0_dp, ' hours'

  timeout = time*timestep/3600.0_dp

  
  if (rank == 0) then

    xdims = (/nxc/)
    ydims=(/nyc/)
    zdims = (/nzc/)

    call h5open_f(ierr)

    call h5fcreate_f(trim(dir_out)//trim(simname)//'chunk_'//trim(tempc)//'.h5',H5F_ACC_TRUNC_F, file_id, ierr)

    call h5ltset_attribute_int_f(file_id,"/","xyzdims", &
(/nxc,nyc,nzc/),3*size, ierr)
    call h5ltset_attribute_double_f(file_id,"/","time", &
timeout, size, ierr)

    call h5screate_simple_f(1, xdims, fspace_id, ierr)
    call h5screate_simple_f(1, xdims, mspace_id, ierr)

    call h5dcreate_f(file_id,"xaxis",H5T_NATIVE_DOUBLE,fspace_id,dset_id,ierr)
    call h5dwrite_f(dset_id, H5T_NATIVE_DOUBLE,xaxis(xcb:xct)*xlength,xdims,ierr)
    call h5dclose_f(dset_id,ierr)

    call h5screate_simple_f(1, ydims, fspace_id, ierr)
    call h5screate_simple_f(1, ydims, mspace_id, ierr)

    call h5dcreate_f(file_id,"yaxis",H5T_NATIVE_DOUBLE,fspace_id,dset_id,ierr)
    call h5dwrite_f(dset_id,H5T_NATIVE_DOUBLE,yaxis(ycb:yct)*ylength,ydims,ierr)
    call h5dclose_f(dset_id,ierr)

    call h5screate_simple_f(1, zdims, fspace_id, ierr)
    call h5screate_simple_f(1, zdims, mspace_id, ierr)

    call h5dcreate_f(file_id,"zaxis",H5T_NATIVE_DOUBLE,fspace_id,dset_id,ierr)
    call h5dwrite_f(dset_id,H5T_NATIVE_DOUBLE,(zaxis(zcb:zct)-1.0_dp)*rsun/1.0e8_dp,zdims,ierr)

    call h5dclose_f(dset_id,ierr)
    call h5sclose_f(fspace_id, ierr)
    call h5sclose_f(mspace_id, ierr)
    call h5fclose_f(file_id,ierr)

  endif

  call MPI_BARRIER(MPI_COMM_WORLD,ierr)

  if (IAM_CHUNKOUT) then

    size_chunk(1) = nxc
    size_chunk(2) = nyc
    size_chunk(3) = nzc
   
    size_mychunk(1) = my_nxc
    size_mychunk(2) = my_nyc
    size_mychunk(3) = my_nzc

    offset(1) = my_offsetx
    offset(2) = my_offsety
    offset(3) = my_offsetz

    call h5pcreate_f(H5P_FILE_ACCESS_F, fapl_id, ierr)
    call h5pset_fapl_mpio_f(fapl_id, comm_cart_savechunk, MPI_INFO_NULL, ierr)
    call h5fopen_f(trim(dir_out)//trim(simname)//'chunk_'//trim(tempc)//'.h5',H5F_ACC_RDWR_F,file_id,ierr,access_prp= fapl_id)
    !! vx
   
    call h5screate_simple_f(3, size_chunk, fspace_id, ierr)
    call h5dcreate_f(file_id, "vx",   H5T_NATIVE_DOUBLE, fspace_id, dset_id,ierr)
    call h5sclose_f(fspace_id, ierr)
    call h5screate_simple_f(3, size_mychunk, mspace_id, ierr)
    call h5dget_space_f(dset_id, fspace_id, ierr)
    call h5sselect_hyperslab_f (fspace_id, H5S_SELECT_SET_F, offset,size_mychunk, ierr)
    call h5pcreate_f(H5P_DATASET_XFER_F, dxpl_id, ierr)
    call h5pset_dxpl_mpio_f(dxpl_id, H5FD_MPIO_COLLECTIVE_F, ierr)
    if (DISPL) then
    call h5dwrite_f(dset_id, H5T_NATIVE_DOUBLE,a(my_xcb:my_xct,my_ycb:my_yct,my_zcb:my_zct,4)*dimc, size_chunk, ierr, &
                 file_space_id = fspace_id, mem_space_id = mspace_id, xfer_prp =dxpl_id)
    else
    call h5dwrite_f(dset_id,H5T_NATIVE_DOUBLE,a(my_xcb:my_xct,my_ycb:my_yct,my_zcb:my_zct,2)*dimc, size_chunk,ierr, &
                 file_space_id = fspace_id, mem_space_id = mspace_id, xfer_prp=dxpl_id)
    endif

    call h5sclose_f(fspace_id,ierr)
    call h5dclose_f(dset_id, ierr)
 
    call h5screate_simple_f(3, size_chunk, fspace_id, ierr)
    call h5dcreate_f(file_id, "vy",   H5T_NATIVE_DOUBLE, fspace_id, dset_id,ierr)
    call h5sclose_f(fspace_id, ierr)
    call h5dget_space_f(dset_id, fspace_id, ierr)
    call h5sselect_hyperslab_f (fspace_id, H5S_SELECT_SET_F,offset,size_mychunk, ierr)
    if (DISPL) then
    call h5dwrite_f(dset_id,H5T_NATIVE_DOUBLE,a(my_xcb:my_xct,my_ycb:my_yct,my_zcb:my_zct,5)*dimc, size_chunk,ierr, &
                 file_space_id = fspace_id, mem_space_id = mspace_id, xfer_prp=dxpl_id)
    else
    call h5dwrite_f(dset_id,H5T_NATIVE_DOUBLE,a(my_xcb:my_xct,my_ycb:my_yct,my_zcb:my_zct,3)*dimc, size_chunk,ierr, &
                 file_space_id = fspace_id, mem_space_id = mspace_id, xfer_prp=dxpl_id)
    endif

    call h5sclose_f(fspace_id,ierr)
    call h5dclose_f(dset_id, ierr)

    call h5screate_simple_f(3, size_chunk, fspace_id, ierr)
    call h5dcreate_f(file_id, "vz",   H5T_NATIVE_DOUBLE, fspace_id, dset_id,ierr)
    call h5sclose_f(fspace_id, ierr)
    call h5dget_space_f(dset_id, fspace_id, ierr)
    call h5sselect_hyperslab_f (fspace_id, H5S_SELECT_SET_F,offset,size_mychunk, ierr)
    if (DISPL) then
    call h5dwrite_f(dset_id,H5T_NATIVE_DOUBLE,a(my_xcb:my_xct,my_ycb:my_yct,my_zcb:my_zct,6)*dimc, size_chunk,ierr, &
                 file_space_id = fspace_id, mem_space_id = mspace_id, xfer_prp=dxpl_id)
    else
    call h5dwrite_f(dset_id,H5T_NATIVE_DOUBLE,a(my_xcb:my_xct,my_ycb:my_yct,my_zcb:my_zct,4)*dimc, size_chunk,ierr, &
                 file_space_id = fspace_id, mem_space_id = mspace_id, xfer_prp=dxpl_id)
    endif
    call h5sclose_f(fspace_id,ierr)
    call h5dclose_f(dset_id, ierr)

    if (DISPL) then
    call h5screate_simple_f(3, size_chunk, fspace_id, ierr)
    call h5dcreate_f(file_id, "xix",   H5T_NATIVE_DOUBLE, fspace_id, dset_id,ierr)
    call h5sclose_f(fspace_id, ierr)
    call h5dget_space_f(dset_id, fspace_id, ierr)
    call h5sselect_hyperslab_f (fspace_id, H5S_SELECT_SET_F,offset,size_mychunk, ierr)
    call h5dwrite_f(dset_id,H5T_NATIVE_DOUBLE,a(my_xcb:my_xct,my_ycb:my_yct,my_zcb:my_zct,1)*diml, size_chunk,ierr, &
                 file_space_id = fspace_id, mem_space_id = mspace_id, xfer_prp=dxpl_id)
    call h5sclose_f(fspace_id,ierr)
    call h5dclose_f(dset_id, ierr)

    call h5screate_simple_f(3, size_chunk, fspace_id, ierr)
    call h5dcreate_f(file_id, "xiy",   H5T_NATIVE_DOUBLE, fspace_id, dset_id,ierr)
    call h5sclose_f(fspace_id, ierr)
    call h5dget_space_f(dset_id, fspace_id, ierr)
    call h5sselect_hyperslab_f (fspace_id, H5S_SELECT_SET_F,offset,size_mychunk, ierr)
    call h5dwrite_f(dset_id,H5T_NATIVE_DOUBLE,a(my_xcb:my_xct,my_ycb:my_yct,my_zcb:my_zct,2)*diml, size_chunk,ierr, &
                 file_space_id = fspace_id, mem_space_id = mspace_id, xfer_prp=dxpl_id)
    call h5sclose_f(fspace_id,ierr)
    call h5dclose_f(dset_id, ierr)

    call h5screate_simple_f(3, size_chunk, fspace_id, ierr)
    call h5dcreate_f(file_id, "xiz",   H5T_NATIVE_DOUBLE, fspace_id,dset_id,ierr)
    call h5sclose_f(fspace_id, ierr)
    call h5dget_space_f(dset_id, fspace_id, ierr)
    call h5sselect_hyperslab_f (fspace_id, H5S_SELECT_SET_F,offset,size_mychunk, ierr)
    call h5dwrite_f(dset_id,H5T_NATIVE_DOUBLE,a(my_xcb:my_xct,my_ycb:my_yct,my_zcb:my_zct,3)*diml, size_chunk,ierr, &
                 file_space_id = fspace_id, mem_space_id = mspace_id, xfer_prp=dxpl_id)
    call h5sclose_f(fspace_id,ierr)
    call h5dclose_f(dset_id, ierr)

    else
    call h5screate_simple_f(3, size_chunk, fspace_id, ierr)
    call h5dcreate_f(file_id, "rho",   H5T_NATIVE_DOUBLE, fspace_id, dset_id,ierr)
    call h5sclose_f(fspace_id, ierr)
    call h5dget_space_f(dset_id, fspace_id, ierr)
    call h5sselect_hyperslab_f (fspace_id, H5S_SELECT_SET_F,offset,size_mychunk, ierr)
    call h5dwrite_f(dset_id,H5T_NATIVE_DOUBLE,a(my_xcb:my_xct,my_ycb:my_yct,my_zcb:my_zct,1)*dimrho, size_chunk,ierr, &
                 file_space_id = fspace_id, mem_space_id = mspace_id, xfer_prp=dxpl_id)
    call h5sclose_f(fspace_id,ierr)
    call h5dclose_f(dset_id, ierr)

    call h5screate_simple_f(3, size_chunk, fspace_id, ierr)
    call h5dcreate_f(file_id, "pre",   H5T_NATIVE_DOUBLE, fspace_id, dset_id,ierr)
    call h5sclose_f(fspace_id, ierr)
    call h5dget_space_f(dset_id, fspace_id, ierr)
    call h5sselect_hyperslab_f (fspace_id, H5S_SELECT_SET_F,offset,size_mychunk, ierr)
    call h5dwrite_f(dset_id,H5T_NATIVE_DOUBLE,a(my_xcb:my_xct,my_ycb:my_yct,my_zcb:my_zct,5)*dimrho*dimc**2, size_chunk,ierr, &
                 file_space_id = fspace_id, mem_space_id = mspace_id, xfer_prp=dxpl_id)
    call h5sclose_f(fspace_id,ierr)
    call h5dclose_f(dset_id, ierr)
    
    if (MAGNETIC) then

    call h5screate_simple_f(3, size_chunk, fspace_id, ierr)
    call h5dcreate_f(file_id, "bx",   H5T_NATIVE_DOUBLE, fspace_id,dset_id,ierr)
    call h5sclose_f(fspace_id, ierr)
    call h5dget_space_f(dset_id, fspace_id, ierr)
    call h5sselect_hyperslab_f (fspace_id, H5S_SELECT_SET_F,offset,size_mychunk,ierr)
    call h5dwrite_f(dset_id,H5T_NATIVE_DOUBLE,a(my_xcb:my_xct,my_ycb:my_yct,my_zcb:my_zct,6)*dimb,size_chunk,ierr, &
                 file_space_id = fspace_id, mem_space_id = mspace_id,xfer_prp=dxpl_id)
    call h5sclose_f(fspace_id,ierr)
    call h5dclose_f(dset_id, ierr)

    call h5screate_simple_f(3, size_chunk, fspace_id, ierr)
    call h5dcreate_f(file_id, "by",   H5T_NATIVE_DOUBLE, fspace_id, dset_id,ierr)
    call h5sclose_f(fspace_id, ierr)
    call h5dget_space_f(dset_id, fspace_id, ierr)
    call h5sselect_hyperslab_f (fspace_id, H5S_SELECT_SET_F,offset,size_mychunk, ierr)
    call h5dwrite_f(dset_id,H5T_NATIVE_DOUBLE,a(my_xcb:my_xct,my_ycb:my_yct,my_zcb:my_zct,7)*dimb, size_chunk,ierr, &
                 file_space_id = fspace_id, mem_space_id = mspace_id, xfer_prp=dxpl_id)
    call h5sclose_f(fspace_id,ierr)
    call h5dclose_f(dset_id, ierr)

    call h5screate_simple_f(3, size_chunk, fspace_id, ierr)
    call h5dcreate_f(file_id, "bz",   H5T_NATIVE_DOUBLE, fspace_id, dset_id,ierr)
    call h5sclose_f(fspace_id, ierr)
    call h5dget_space_f(dset_id, fspace_id, ierr)
    call h5sselect_hyperslab_f (fspace_id, H5S_SELECT_SET_F,offset,size_mychunk, ierr)
    call  h5dwrite_f(dset_id,H5T_NATIVE_DOUBLE,a(my_xcb:my_xct,my_ycb:my_yct,my_zcb:my_zct,8)*dimb, size_chunk,ierr, &
                 file_space_id = fspace_id, mem_space_id = mspace_id, xfer_prp=dxpl_id)
    call h5sclose_f(fspace_id,ierr)
    call h5dclose_f(dset_id,ierr)
    endif
    endif
    call h5pclose_f(dxpl_id, ierr)
    call h5sclose_f(mspace_id,ierr)
    call h5pclose_f(fapl_id, ierr)
    call h5fclose_f(file_id, ierr)
  endif

  call MPI_BARRIER(MPI_COMM_WORLD,ierr)

end subroutine WRITE_OUT_CHUNK

!================================================================================================

  subroutine damp_velocity (arr,term)

    real(dp), intent(in), dimension(dim1(rank),dim2(rank),dim3(rank)) :: arr
    real(dp), intent(out), dimension(dim1(rank),dim2(rank),dim3(rank)) :: term
    integer ii,jj,kk
    integer ierr,recrank
    integer, allocatable, dimension(:) :: req_send,req_recv, req2_send,req2_recv
    integer, allocatable, dimension(:,:) :: stat_send,stat_recv,stat2_send,stat2_recv
    integer counter, count_1

    count_1 = min(dim3(rank),block_size(1)*block_size(2))

    allocate(req_recv(block_size(1)*block_size(2)),stat_recv(MPI_STATUS_SIZE,block_size(1)*block_size(2)), &
    req2_send(block_size(1)*block_size(2)),stat2_send(MPI_STATUS_SIZE,block_size(1)*block_size(2)), &
    req2_recv(count_1),stat2_recv(MPI_STATUS_SIZE,count_1), &
    req_send(count_1),stat_send(MPI_STATUS_SIZE,count_1))

    !! Send undamped

    kk = 1
    counter = 1
    do ii = 0,block_size(1)-1
      do jj = 0,block_size(2)-1
          recrank = rank_xyz_all(ii,jj,coords_xyz(3))

          if (dim_damp(rank) .GT. 0) then
            ! recv
            call MPI_IRECV(transfermatrix(xstart(recrank),ystart(recrank),1),1,damp_1, recrank,recrank,comm_cart_3D,&
            req_recv(ii+1+block_size(1)*(jj)),ierr)
          endif
          
          ! send  
       
         if (dim_damp(recrank) .GT. 0) then
           call MPI_ISEND(arr(1,1,kk),1, damp_2(ii,jj),recrank,rank,comm_cart_3D, req_send(counter),ierr)
           counter = counter + 1
           ! count
           kk = kk+dim_damp(recrank)
        endif
      enddo
    enddo

    if (dim_damp(rank) .GT. 0) then
      !! WAIT
      call MPI_WAITALL(block_size(1)*block_size(2),req_recv, stat_recv,ierr)

      !! do damping
      call dfftw_execute_dft_r2c(fftw_plan_fwd_2D, transfermatrix, transfertemp)
      do kk = 1,dim_damp(rank)
        transfertemp(:,:,kk) = transfertemp(:,:,kk) * damping_rates
      enddo
      call dfftw_execute_dft_c2r(fftw_plan_inv_2D, transfertemp, transfermatrix)
    endif

    call MPI_WAITALL(count_1,req_send,stat_send,ierr)
    !! Return Damped

    kk = 1
    counter = 1
    do ii = 0,block_size(1)-1
      do jj = 0,block_size(2)-1
          recrank = rank_xyz_all(ii,jj,coords_xyz(3))
 
          if (dim_damp(rank) .GT. 0) then
            ! send
            call MPI_ISEND(transfermatrix(xstart(recrank),ystart(recrank),1),1,damp_1, recrank,recrank,comm_cart_3D,&
            req2_send(ii+block_size(1)*(jj)+1),ierr)
          endif
          ! recv           
          if (dim_damp(recrank) .GT. 0) then
            call MPI_IRECV(term(1,1,kk),1,damp_2(ii,jj),recrank,rank,comm_cart_3D, req2_recv(counter),ierr)
          ! count
          counter = counter+1
          kk = kk+dim_damp(recrank)
        endif
      enddo
    enddo

    !! WAIT
    if (dim_damp(rank) .GT. 0) then
      call MPI_WAITALL(block_size(1)*block_size(2),req2_send, stat2_send,ierr)
    endif

    call MPI_WAITALL(count_1,req2_recv,stat2_recv,ierr)

  end subroutine damp_velocity

!================================================================================================

  subroutine LAGRANGE_INTERP()

  implicit none
  real(dp) xt

  xt = DBLE(time)*timestep/cadforcing
  time1 = FLOOR(xt)+1

  xt = time*timestep

  if (time_old .LT. time1) then
   x0 = (time1-4)*cadforcing
   x1 = (time1-3)*cadforcing
   x2 = (time1-2)*cadforcing
   x3 = (time1-1)*cadforcing
   x4 = (time1)*cadforcing
   x5 = (time1+1)*cadforcing
   x6 = (time1+2)*cadforcing

   LC0 = 0.0
   LC1 = 0.0
   LC2 = 0.0
   LC3 = 0.0
   LC4 = 0.0
   LC5 = 0.0
   LC6 = 0.0

   if (x0 .GE. 0) &
    LC0 = vr(:,:,(time1-3))/((x0-x1) * (x0-x2) * (x0-x3) * (x0-x4) * (x0-x5) *(x0-x6))

   if (x1 .GE. 0) &
    LC1 = vr(:,:,(time1-2))/((x1-x0) * (x1-x2) * (x1-x3) * (x1-x4) * (x1-x5) *(x1-x6))

   if (x2 .GE. 0) &
    LC2 = vr(:,:,(time1-1))/((x2-x0) * (x2-x1) * (x2-x3) * (x2-x4) * (x2-x5) *(x2-x6))

   if (x3 .GE. 0) &
    LC3 = vr(:,:,(time1))/((x3-x0) * (x3-x1) * (x3-x2) * (x3-x4) * (x3-x5) *(x3-x6))

   LC4 = vr(:,:,(time1+1))/((x4-x0) * (x4-x1) * (x4-x2) * (x4-x3) * (x4-x5) *(x4-x6))

   if (time1 .LT. num_steps) &
    LC5 = vr(:,:,(time1+2))/((x5-x0) * (x5-x1) * (x5-x2) * (x5-x3) * (x5-x4) *(x5-x6))

   if (time1 .LT. (num_steps-1)) &
    LC6 = vr(:,:,(time1+3))/((x6-x0) * (x6-x1) * (x6-x2) * (x6-x3) * (x6-x4) *(x6-x5))

   time_old = time1
  endif

  forcing = LC0 * (xt-x1) * (xt-x2) * (xt-x3) * (xt-x4) * (xt-x5) *(xt-x6)+ &
        LC1 * (xt-x0) * (xt-x2) * (xt-x3) * (xt-x4) * (xt-x5) *(xt-x6)+ &
      LC2 * (xt-x0) * (xt-x1) * (xt-x3) * (xt-x4) * (xt-x5) *(xt-x6)+ &
      LC3 * (xt-x0) * (xt-x1) * (xt-x2) * (xt-x4) * (xt-x5) *(xt-x6)+ &
      LC4 * (xt-x0) * (xt-x1) * (xt-x2) * (xt-x3) * (xt-x5) *(xt-x6)+ &
      LC5 * (xt-x0) * (xt-x1) * (xt-x2) * (xt-x3) * (xt-x4) *(xt-x6)+ &
      LC6 * (xt-x0) * (xt-x1) * (xt-x2) * (xt-x3) * (xt-x4) *(xt-x5)

end subroutine LAGRANGE_INTERP

!==================================================================================

end module all_modules

