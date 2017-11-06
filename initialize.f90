MODULE INITIALIZE

use ALL_MODULES
use derivatives

implicit none

CONTAINS
!==========================================================================


 subroutine INIT_MPI()

   integer chunk, chunk2D, xyz_send, rank_o,ierr, sourcechunk, damp_temp
   integer remtemp
   integer, dimension(1000) :: dimtemp
   integer (KIND=MPI_ADDRESS_KIND) dubdum,intdum,extent2D,extentofcol,sizeofdouble
   integer (KIND=MPI_ADDRESS_KIND) sizeofint
   integer ii,kk
   integer, dimension(3) :: blocksize3D,xdersize3D,&
   xfiltsize3D,ydersize3D,yfiltsize3D,zdersize3D,zfiltsize3D, parcel_size_x_d, &
   parcel_size_x_f, parcel_size_y_d, parcel_size_y_f, parcel_size_z_d, parcel_size_z_f

   call MPI_INIT(ierr)
   call MPI_COMM_RANK(MPI_COMM_WORLD, rank_o, ierr)
   call MPI_COMM_SIZE(MPI_COMM_WORLD, numtasks, ierr)
   call MPI_TYPE_GET_EXTENT(MPI_DOUBLE_PRECISION,dubdum,sizeofdouble,ierr)
   call MPI_TYPE_GET_EXTENT(MPI_INTEGER,intdum,sizeofint,ierr)
   rank = rank_o
   call READ_SPARC_PARAMETERS()

   if (PERIODIC) then
       mpi_periodic = (/.true.,.true.,.false./)
   else
       mpi_periodic =   (/.false.,.false.,.false./)
   endif

   !! Calculate block-size from available number of cores.
  
   if (numtasks .NE. block_size(1)*block_size(2)*block_size(3)) then
     if (rank .EQ. 0) print *, 'number of cores ',numtasks, ' not equal to B_x*B_y*B_z', block_size(1),block_size(2),block_size(3), block_size(1)*block_size(2)*block_size(3)
     call MPI_BARRIER(MPI_COMM_WORLD, ierr)
     call MPI_FINALIZE(ierr)
   endif

   call MPI_CART_CREATE(MPI_COMM_WORLD,3,block_size,mpi_periodic,.TRUE.,comm_cart_3D,ierr)
   call MPI_COMM_RANK(comm_cart_3D,rank,ierr)
   call MPI_CART_COORDS(comm_cart_3D,rank,3,coords_xyz,ierr)

   call MPI_CART_SUB(comm_cart_3D, (/.FALSE.,.TRUE.,.TRUE./), xsplit_comm,ierr)
   call MPI_CART_SUB(comm_cart_3D, (/.TRUE.,.FALSE.,.TRUE./), ysplit_comm,ierr)
   call MPI_CART_SUB(comm_cart_3D, (/.TRUE.,.TRUE.,.FALSE./), zsplit_comm,ierr)

   !! coords_xyz gives the x,y,z position of each core.

   IAM_XBOT = .FALSE.
   IAM_XTOP = .FALSE.
   IAM_YBOT = .FALSE.
   IAM_YTOP = .FALSE.
   IAM_ZBOT = .FALSE.
   IAM_ZTOP = .FALSE.
   IAM_NOTB = .TRUE.

   !!! If i am a boundary, then which am i, and set IAM_NOTB to false

   if (coords_xyz(1) .EQ. 0) then ; IAM_XBOT = .TRUE. ; IAM_NOTB = .FALSE. ; endif
   if (coords_xyz(1) .EQ. block_size(1)-1) then ; IAM_XTOP = .TRUE. ; IAM_NOTB = .FALSE. ; endif 
   if ((.not. TWOD) .AND. (coords_xyz(2) .EQ. 0)) then ; IAM_YBOT = .TRUE. ; IAM_NOTB = .FALSE. ; endif 
   if ((.not. TWOD) .AND. (coords_xyz(2) .EQ. block_size(2)-1)) then ; IAM_YTOP = .TRUE. ; IAM_NOTB = .FALSE. ; endif
   if (coords_xyz(3) .EQ. 0) then ; IAM_ZBOT = .TRUE. ; IAM_NOTB = .FALSE. ; endif
   if (coords_xyz(3) .EQ. block_size(3)-1) then ; IAM_ZTOP = .TRUE. ; IAM_NOTB = .FALSE. ; endif
 

   allocate(dim1(0:numtasks-1), dim2(0:numtasks-1), dim3(0:numtasks-1),xstart(0:numtasks-1),ystart(0:numtasks-1),zstart(0:numtasks-1),&
   xend(0:numtasks-1),yend(0:numtasks-1),zend(0:numtasks-1))

   dimtemp = nx/block_size(1)
   remtemp = nx-dimtemp(1)*block_size(1)

   do ii=1,remtemp
     dimtemp(ii) =  dimtemp(ii) +1
   enddo
   dim1(rank) = dimtemp(coords_xyz(1)+1)

   xstart(rank) = 1
   xend(rank) = dimtemp(1)   
   do ii=1,coords_xyz(1)
     xstart(rank) = xstart(rank)+dimtemp(ii)
     xend(rank) = xend(rank) + dimtemp(ii+1) 
   enddo

   dimtemp = ny/block_size(2)
   remtemp = ny-dimtemp(1)*block_size(2)

   do ii=1,remtemp
     dimtemp(ii) =  dimtemp(ii) +1
   enddo
   dim2(rank) = dimtemp(coords_xyz(2)+1)

   ystart(rank) = 1
   yend(rank) = dimtemp(1)      
   do ii=1,coords_xyz(2)
     ystart(rank) = ystart(rank)+dimtemp(ii)
     yend(rank) = yend(rank) + dimtemp(ii+1)  
   enddo
   
   dimtemp = nz/block_size(3)
   remtemp = nz-dimtemp(1)*block_size(3)

   do ii=1,remtemp                
     dimtemp(ii) =  dimtemp(ii) +1
   enddo                      
   dim3(rank) = dimtemp(coords_xyz(3)+1)

   zstart(rank) = 1
   zend(rank) = dimtemp(1)      
   do ii=1,coords_xyz(3)
     zstart(rank) = zstart(rank)+dimtemp(ii)
     zend(rank) = zend(rank) + dimtemp(ii+1)  
   enddo

   allocate(coords_xyz_all(3,0:numtasks-1), rank_xyz_all(0:block_size(1)-1,0:block_size(2)-1, 0:block_size(3)-1))
   
   !!set of dim1 columns
   call MPI_TYPE_VECTOR(3,1,1,MPI_INTEGER,xyz_send,ierr)
   !!set of dim2 columns with a stride of nx
   call MPI_TYPE_COMMIT(xyz_send,ierr)

   call MPI_ALLGATHER(coords_xyz, 3, MPI_INTEGER,coords_xyz_all,3,MPI_INTEGER, comm_cart_3D,ierr)

   do ii = 0,numtasks-1
     rank_xyz_all(coords_xyz_all(1,ii),coords_xyz_all(2,ii),coords_xyz_all(3,ii)) = ii
   enddo

   call MPI_Allgather(dim1(rank), 1, MPI_INTEGER,dim1,1,MPI_INTEGER,comm_cart_3D,ierr)
   call MPI_Allgather(dim2(rank), 1, MPI_INTEGER,dim2,1,MPI_INTEGER,comm_cart_3D,ierr)
   call MPI_Allgather(dim3(rank), 1, MPI_INTEGER,dim3,1,MPI_INTEGER,comm_cart_3D,ierr)
   call MPI_Allgather(xstart(rank), 1, MPI_INTEGER,xstart,1,MPI_INTEGER,comm_cart_3D,ierr)
   call MPI_Allgather(ystart(rank), 1, MPI_INTEGER,ystart,1,MPI_INTEGER,comm_cart_3D,ierr)
   call MPI_Allgather(zstart(rank), 1, MPI_INTEGER,zstart,1,MPI_INTEGER,comm_cart_3D,ierr)
   call MPI_Allgather(xend(rank), 1, MPI_INTEGER,xend,1,MPI_INTEGER,comm_cart_3D,ierr)
   call MPI_Allgather(yend(rank), 1, MPI_INTEGER,yend,1,MPI_INTEGER,comm_cart_3D,ierr)
   call MPI_Allgather(zend(rank), 1, MPI_INTEGER,zend,1,MPI_INTEGER,comm_cart_3D,ierr)

   if ((MINVAL(dim1) .LT. 20) .or. ((.not. TWOD) .AND. (MINVAL(dim2) .LT. 20)) .or. (MINVAL(dim3) .LT. 20)) then
     if (rank .EQ. 0) then
       print *, 'number of cells per block is too small for efficient computation'
       print *, minval(dim1), '  in x.  ',minval(dim2), ' in y.  ', minval(dim3), ' in z.'
       print *, 'reccomend decreasing b_x,b_y or b_z until it is greater than or equal to 20'
       print *, 'or using a smaller stencil derivative than the current one'
     endif
     call MPI_BARRIER(MPI_COMM_WORLD, ierr)          
     call MPI_FINALIZE(ierr)
   endif


   blocksize3D = (/dim1(rank),dim2(rank),dim3(rank)/)
   
   xdersize3D = (/dim1(rank)+2*g_xy_der,dim2(rank),dim3(rank)/)
   xfiltsize3D = (/dim1(rank) + 2*g_xy_filt,dim2(rank),dim3(rank)/)
   
   ydersize3D = (/dim1(rank),dim2(rank)+2*g_xy_der,dim3(rank)/)
   yfiltsize3D = (/dim1(rank),dim2(rank) + 2*g_xy_filt,dim3(rank)/)
   
   zdersize3D = (/dim1(rank),dim2(rank),dim3(rank)+2*g_z_der/)
   zfiltsize3D = (/dim1(rank),dim2(rank),dim3(rank) + 2*g_z_filt/)
   
   call MPI_BARRIER(comm_cart_3D,ierr)

  !! Create datatypes, we need: x,y,z 5 cell for derivatives, x,y, 3-cell for
  !! filtering, z 5-cell for filtering (so can use derivatives)
  !! Using g_xy_der,g_z_der,g_xy_filt, g_z_filt
  !! 3Dblocksize contains the array of data
  !! To be honest, this is uglier than just making the data-types the easy way with
  !! MPI_TYPE_HVECTOR etc. But I wanted to see how it works, so it is done now, it works & I am to lazy to do it again.
  
  parcel_size_x_d = (/g_xy_der,dim2(rank),dim3(rank)/)
  parcel_size_y_d = (/dim1(rank),g_xy_der,dim3(rank)/)
  parcel_size_z_d = (/dim1(rank),dim2(rank),g_z_der/)
 
  parcel_size_x_f = (/g_xy_filt,dim2(rank),dim3(rank)/)
  parcel_size_y_f = (/dim1(rank),g_xy_filt,dim3(rank)/)
  parcel_size_z_f = (/dim1(rank),dim2(rank),g_z_filt/)

  !! NB MPI datatypes start at 0

  !! START X DATATYPES

  

  !! X, left, send, derivative
  starts = (/0,0,0/)

  call MPI_TYPE_CREATE_SUBARRAY(3,blocksize3D,parcel_size_x_d,starts,MPI_ORDER_FORTRAN,MPI_DOUBLE_PRECISION,parcel_xlds,ierr)
  call MPI_TYPE_COMMIT(parcel_xlds,ierr)

  !! X, left, send, filter

  call MPI_TYPE_CREATE_SUBARRAY(3,blocksize3D,parcel_size_x_f,starts,MPI_ORDER_FORTRAN,MPI_DOUBLE_PRECISION,parcel_xlfs,ierr)
  call MPI_TYPE_COMMIT(parcel_xlfs,ierr)
 
  !! X, left, recv, derivative

  call MPI_TYPE_CREATE_SUBARRAY(3,xdersize3D,parcel_size_x_d,starts,MPI_ORDER_FORTRAN,MPI_DOUBLE_PRECISION,parcel_xldr,ierr)
  call MPI_TYPE_COMMIT(parcel_xldr,ierr)

  !! X, left, recv, filter

  call MPI_TYPE_CREATE_SUBARRAY(3,xfiltsize3D,parcel_size_x_f,starts,MPI_ORDER_FORTRAN,MPI_DOUBLE_PRECISION,parcel_xlfr,ierr)
  call MPI_TYPE_COMMIT(parcel_xlfr,ierr)

  !! X, right, send, derivative
  starts = (/dim1(rank)-g_xy_der,0,0/)
  
  call MPI_TYPE_CREATE_SUBARRAY(3,blocksize3D,parcel_size_x_d,starts,MPI_ORDER_FORTRAN,MPI_DOUBLE_PRECISION,parcel_xrds,ierr)
  call MPI_TYPE_COMMIT(parcel_xrds,ierr)

  !! X, right, send, filter
  starts = (/dim1(rank)-g_xy_filt,0,0/)
  
  call MPI_TYPE_CREATE_SUBARRAY(3,blocksize3D,parcel_size_x_f,starts,MPI_ORDER_FORTRAN,MPI_DOUBLE_PRECISION,parcel_xrfs,ierr)
  call MPI_TYPE_COMMIT(parcel_xrfs,ierr)

  !! X, right, recv, derivative
  starts = (/dim1(rank)+g_xy_der,0,0/)
  
  call MPI_TYPE_CREATE_SUBARRAY(3,xdersize3D,parcel_size_x_d,starts,MPI_ORDER_FORTRAN,MPI_DOUBLE_PRECISION,parcel_xrdr,ierr)
  call MPI_TYPE_COMMIT(parcel_xrdr,ierr)

  !! X, right, recv, filter
  starts = (/dim1(rank)+g_xy_filt,0,0/)

  call MPI_TYPE_CREATE_SUBARRAY(3,xfiltsize3D,parcel_size_x_f,starts,MPI_ORDER_FORTRAN,MPI_DOUBLE_PRECISION,parcel_xrfr,ierr)
  call MPI_TYPE_COMMIT(parcel_xrfr,ierr)

  if (.not. TWOD) then
  !!! START Y Datatypes

  starts = (/0,0,0/)

  !! Y, left, send, derivative

  call MPI_TYPE_CREATE_SUBARRAY(3,blocksize3D,parcel_size_y_d,starts,MPI_ORDER_FORTRAN,MPI_DOUBLE_PRECISION,parcel_ylds,ierr)
  call MPI_TYPE_COMMIT(parcel_ylds,ierr)

  !! Y, left, send, filter

  call MPI_TYPE_CREATE_SUBARRAY(3,blocksize3D,parcel_size_y_f,starts,MPI_ORDER_FORTRAN,MPI_DOUBLE_PRECISION,parcel_ylfs,ierr)
  call MPI_TYPE_COMMIT(parcel_ylfs,ierr)
 
  !! Y, left, recv, derivative
 
  call MPI_TYPE_CREATE_SUBARRAY(3,ydersize3D,parcel_size_y_d,starts,MPI_ORDER_FORTRAN,MPI_DOUBLE_PRECISION,parcel_yldr,ierr)
  call MPI_TYPE_COMMIT(parcel_yldr,ierr)

  !! Y, left, recv, filter

  call MPI_TYPE_CREATE_SUBARRAY(3,yfiltsize3D,parcel_size_y_f,starts,MPI_ORDER_FORTRAN,MPI_DOUBLE_PRECISION,parcel_ylfr,ierr)
  call MPI_TYPE_COMMIT(parcel_ylfr,ierr)

  starts = (/0,dim2(rank)-g_xy_der,0/)

  !! Y, right, send, derivative

  call MPI_TYPE_CREATE_SUBARRAY(3,blocksize3D,parcel_size_y_d,starts,MPI_ORDER_FORTRAN,MPI_DOUBLE_PRECISION,parcel_yrds,ierr)
  call MPI_TYPE_COMMIT(parcel_yrds,ierr)

  starts = (/0,dim2(rank)-g_xy_filt,0/)

  !! Y, right, send, filter

  call MPI_TYPE_CREATE_SUBARRAY(3,blocksize3D,parcel_size_y_f,starts,MPI_ORDER_FORTRAN,MPI_DOUBLE_PRECISION,parcel_yrfs,ierr)
  call MPI_TYPE_COMMIT(parcel_yrfs,ierr)

  !! Y, right, recv, derivative
  starts = (/0,dim2(rank)+g_xy_der,0/)
  
  call MPI_TYPE_CREATE_SUBARRAY(3,ydersize3D,parcel_size_y_d,starts,MPI_ORDER_FORTRAN,MPI_DOUBLE_PRECISION,parcel_yrdr,ierr)
  call MPI_TYPE_COMMIT(parcel_yrdr,ierr)

  !! Y, right, recv, filter

  starts = (/0,dim2(rank)+g_xy_filt,0/)

  call MPI_TYPE_CREATE_SUBARRAY(3,yfiltsize3D,parcel_size_y_f,starts,MPI_ORDER_FORTRAN,MPI_DOUBLE_PRECISION,parcel_yrfr,ierr)
  call MPI_TYPE_COMMIT(parcel_yrfr,ierr)
  
  endif

  !!! START Z Datatypes

  starts = (/0,0,0/)

  !! Z, left, send, derivative

  call MPI_TYPE_CREATE_SUBARRAY(3,blocksize3D,parcel_size_z_d,starts,MPI_ORDER_FORTRAN,MPI_DOUBLE_PRECISION,parcel_zlds,ierr)
  call MPI_TYPE_COMMIT(parcel_zlds,ierr)

  !! Z, left, send, filter

  call MPI_TYPE_CREATE_SUBARRAY(3,blocksize3D,parcel_size_z_f,starts,MPI_ORDER_FORTRAN,MPI_DOUBLE_PRECISION,parcel_zlfs,ierr)
  call MPI_TYPE_COMMIT(parcel_zlfs,ierr)
 
  !! Z, left, recv, derivative
 
  call MPI_TYPE_CREATE_SUBARRAY(3,zdersize3D,parcel_size_z_d,starts,MPI_ORDER_FORTRAN,MPI_DOUBLE_PRECISION,parcel_zldr,ierr)
  call MPI_TYPE_COMMIT(parcel_zldr,ierr)

  !! Z, left, recv, filter

  call MPI_TYPE_CREATE_SUBARRAY(3,zfiltsize3D,parcel_size_z_f,starts,MPI_ORDER_FORTRAN,MPI_DOUBLE_PRECISION,parcel_zlfr,ierr)
  call MPI_TYPE_COMMIT(parcel_zlfr,ierr)

  starts = (/0,0,dim3(rank)-g_z_der/)

  !! Z, right, send, derivative

  call MPI_TYPE_CREATE_SUBARRAY(3,blocksize3D,parcel_size_z_d,starts,MPI_ORDER_FORTRAN,MPI_DOUBLE_PRECISION,parcel_zrds,ierr)
  call MPI_TYPE_COMMIT(parcel_zrds,ierr)

  starts = (/0,0,dim3(rank)-g_z_filt/)

  !! Z, right, send, filter

  call MPI_TYPE_CREATE_SUBARRAY(3,blocksize3D,parcel_size_z_f,starts,MPI_ORDER_FORTRAN,MPI_DOUBLE_PRECISION,parcel_zrfs,ierr)
  call MPI_TYPE_COMMIT(parcel_zrfs,ierr)

  !! Z, right, recv, derivative
  starts = (/0,0,dim3(rank)+g_z_der/)
  
  call MPI_TYPE_CREATE_SUBARRAY(3,zdersize3D,parcel_size_z_d,starts,MPI_ORDER_FORTRAN,MPI_DOUBLE_PRECISION,parcel_zrdr,ierr)
  call MPI_TYPE_COMMIT(parcel_zrdr,ierr)

  !! Z, right, recv, filter

  starts = (/0,0,dim3(rank)+g_z_filt/)

  call MPI_TYPE_CREATE_SUBARRAY(3,zfiltsize3D,parcel_size_z_f,starts,MPI_ORDER_FORTRAN,MPI_DOUBLE_PRECISION,parcel_zrfr,ierr)
  call MPI_TYPE_COMMIT(parcel_zrfr,ierr)
 
  !! MPI datatype for receiving (or sending) a dim1xdim2xdim3 sub-cube of a
  !! nx,ny,nz array

  allocate(mysource3D(0:numtasks-1))

  do kk = 0,numtasks-1
    !!set of dim1 columns
    call MPI_TYPE_VECTOR(dim1(kk),1,1,MPI_DOUBLE_PRECISION,chunk,ierr)
    extentofcol = nx*sizeofdouble 
    call MPI_TYPE_CREATE_HVECTOR(dim2(kk),1,extentofcol,chunk,chunk2D,ierr)
    extent2D = nx*ny*sizeofdouble
    call MPI_TYPE_CREATE_HVECTOR(num_steps,1,extent2D,chunk2D,sourcechunk,ierr)
    call MPI_TYPE_COMMIT(sourcechunk,ierr)
    mysource3D(kk) = sourcechunk
  enddo

  call MPI_TYPE_VECTOR(dim1(rank),1,1,MPI_DOUBLE_PRECISION,chunk,ierr)
  extentofcol = dim1(rank)*sizeofdouble
  
  call MPI_TYPE_CREATE_HVECTOR(dim2(rank),1,extentofcol,chunk,chunk2D,ierr)
  extent2D = dim1(rank)*dim2(rank)*sizeofdouble

  call MPI_TYPE_CREATE_HVECTOR(nzpmltop,1,extent2D,chunk2D,chunkpmltop,ierr)
  call MPI_TYPE_COMMIT(chunkpmltop,ierr)

  call MPI_TYPE_CREATE_HVECTOR(nzpmlbot,1,extent2D,chunk2D,chunkpmlbot,ierr)
  call MPI_TYPE_COMMIT(chunkpmlbot,ierr)

  call MPI_TYPE_VECTOR(dim1(rank),1,1,MPI_DOUBLE_PRECISION,chunk,ierr)
  extentofcol = nx*sizeofdouble

  call MPI_TYPE_CREATE_HVECTOR(dim2(rank),1,extentofcol,chunk,chunk2D,ierr)
  extent2D = nx*ny*sizeofdouble

  call MPI_TYPE_CREATE_HVECTOR(nzpmltop,1,extent2D,chunk2D,allpmltop,ierr)
  call MPI_TYPE_COMMIT(allpmltop,ierr)

  call MPI_TYPE_CREATE_HVECTOR(nzpmlbot,1,extent2D,chunk2D,allpmlbot,ierr)
  call MPI_TYPE_COMMIT(allpmlbot,ierr)


  call MPI_TYPE_CREATE_HVECTOR(dim3(rank),1,dim1(rank)*dim2(rank)*sizeofdouble,MPI_DOUBLE_PRECISION,&
  datatype_zcol,ierr)
  call MPI_TYPE_COMMIT(datatype_zcol,ierr)

  !! We want to split all cores with the same vertical rank into numcores_x*numcores_y
  !! even-ishly distributed xy slices for FFTS.

  if (DAMPING) then

    allocate(dim_damp(0:numtasks-1),damp_2(0:block_size(1)-1,0:block_size(2)-1))
    if (dim3(rank) .LE. block_size(1)*block_size(2)) then

      dim_damp = 0
      dim_damp(0:dim3(rank)-1) = 1
    else

      dimtemp = dim3(rank)/block_size(1)/block_size(2)
      remtemp = dim3(rank)-dimtemp(1)*block_size(1)*block_size(2)

      ii = 1
      kk = 1

      do while (remtemp .GT. 0)
        
        if (ii .GT. block_size(1)) then
          ii = 1
          kk = kk+1
        endif

        dimtemp(ii+block_size(1)*(kk-1)) =  dimtemp(ii+block_size(1)*(kk-1)) + 1
        remtemp = remtemp-1
        ii = ii + 1
      enddo

      dim_damp(rank) = dimtemp(coords_xyz(1)+1 + block_size(1)*coords_xyz(2))

      !!! share dim_damp

      call MPI_Allgather(dim_damp(rank), 1, MPI_INTEGER,dim_damp,1,MPI_INTEGER,comm_cart_3D,ierr)

    endif

    !! Create datatypes

    !! First datatype, dim1(rank)*dim2(rank)*dim_damp(rank) of a nx*ny*dim_damp(rank) array
    
    call MPI_TYPE_VECTOR(dim1(rank),1,1,MPI_DOUBLE_PRECISION,chunk,ierr)
    extentofcol = nx*sizeofdouble
    call MPI_TYPE_CREATE_HVECTOR(dim2(rank),1,extentofcol,chunk,chunk2D,ierr)
    extent2D = nx*ny*sizeofdouble
    call MPI_TYPE_CREATE_HVECTOR(dim_damp(rank),1,extent2D,chunk2D,damp_1,ierr)
    call MPI_TYPE_COMMIT(damp_1,ierr)

    !! second datatype, a set of datatypes to send xyslices to different cores

    do ii = 0,block_size(1)-1
      do kk = 0,block_size(2)-1

          if (dim_damp(rank_xyz_all(ii,kk,coords_xyz(3))) .GT. 0) then
          call MPI_TYPE_VECTOR(dim1(rank),1,1,MPI_DOUBLE_PRECISION,chunk,ierr)
          extentofcol = dim1(rank)*sizeofdouble
          call MPI_TYPE_CREATE_HVECTOR(dim2(rank),1,extentofcol,chunk,chunk2D,ierr)
          extent2D = dim1(rank)*dim2(rank)*sizeofdouble
          call MPI_TYPE_CREATE_HVECTOR(dim_damp(rank_xyz_all(ii,kk,coords_xyz(3))),1,extent2D,chunk2D,damp_temp,ierr)
          call MPI_TYPE_COMMIT(damp_temp,ierr)

          damp_2(ii,kk) = damp_temp

          elseif (dim_damp(rank_xyz_all(ii,kk,coords_xyz(3))) .EQ. 0) then
           !! no more datatype needed
          endif
      enddo
    enddo
 endif

  call MPI_BARRIER(comm_cart_3D,ierr)

  end subroutine INIT_MPI

!==========================================================================

  subroutine INITANDALLOCATE ()
     
  integer :: ii,jj,kk
  integer ierr

  allocate(z(dim3(rank)),xaxis(nx),yaxis(ny),zaxis(nz),x(dim1(rank)),y(dim2(rank)))
  
  if (.NOT. (periodic)) then
    do ii=1,nx
      xaxis(ii) = DBLE(ii-1_dp)/DBLE(nx - 1_dp)
    enddo
  else
    do ii=1,nx
      xaxis(ii) = DBLE(ii-1_dp)/DBLE(nx)
    enddo
  endif

  if (TWOD) then
    yaxis = 0.0_dp
  else
    if (.NOT. (periodic)) then
      do jj=1,ny
        yaxis(jj) = DBLE(jj-1_dp)/DBLE(ny - 1_dp)
      enddo
    else
      do jj=1,ny
        yaxis(jj) =  DBLE(jj-1_dp)/DBLE(ny)
      enddo
    endif
  endif

  x = xaxis(xstart(rank):xend(rank))
  y = yaxis(ystart(rank):yend(rank))
 
  allocate(spongex(dim1(rank)), spongey(dim2(rank)), spongez(dim3(rank)),&
  cq(dim3(rank)), stretch(dim3(rank)), unstretch(dim3(rank)), &
  invrho0q(dim3(rank)), g(dim3(rank)), source_dep(dim3(rank)))

  call PARSE_OPTIONS ()

  if (periodic) then
    stretchx = rsun/xlength*nx
  else
    stretchx = rsun/xlength*(nx-1)
  endif

  if (.not. TWOD) then
    if (periodic) then
      stretchy = rsun/ylength*ny
    else
      stretchy = rsun/ylength*(ny-1)
    endif
  else 
    stretchy = 1_dp
  endif

  ! The size of the fifth dimension

  dimfive = 5
  if (magnetic) dimfive = 8
  if (displ) dimfive = 6

  if (magnetic) then 
    if (.not. displ) then 
      allocate(a(dim1(rank),dim2(rank),dim3(rank),8), temp_step(dim1(rank),dim2(rank),dim3(rank),8),&
      scr(dim1(rank),dim2(rank),dim3(rank),8))
    else
      allocate(a(dim1(rank),dim2(rank),dim3(rank),6), temp_step(dim1(rank),dim2(rank),dim3(rank),6),&
      scr(dim1(rank),dim2(rank),dim3(rank),6),scratch(dim1(rank),dim2(rank),dim3(rank),5))
    endif

    allocate(box(dim1(rank),dim2(rank),dim3(rank)), boy(dim1(rank), dim2(rank), dim3(rank)), boz(dim1(rank), dim2(rank), dim3(rank)),&
    curlbox(dim1(rank),dim2(rank),dim3(rank)), curlboy(dim1(rank),dim2(rank),dim3(rank)), curlboz(dim1(rank),dim2(rank),dim3(rank)),&
    curlbx(dim1(rank),dim2(rank),dim3(rank)), curlby(dim1(rank),dim2(rank),dim3(rank)), curlbz(dim1(rank),dim2(rank),dim3(rank)),&
    gradrho0_x(dim1(rank), dim2(rank), dim3(rank)), gradrho0_y(dim1(rank), dim2(rank), dim3(rank)), &
    gradp0_x(dim1(rank), dim2(rank), dim3(rank)),&
    gradp0_y(dim1(rank), dim2(rank), dim3(rank)),reduction(dim1(rank),dim2(rank),dim3(rank)))

    if (USE_PML .AND. IAM_ZBOT) allocate(scrzpml(dim1(rank),dim2(rank),nzpmlbot,6),zpmlvars(dim1(rank),dim2(rank),nzpmlbot,6), &
    psizpml(dim1(rank),dim2(rank),nzpmlbot,6))
    
    if (USE_PML .AND. IAM_ZTOP) allocate(scrzpml(dim1(rank),dim2(rank),nzpmltop,6),zpmlvars(dim1(rank),dim2(rank),nzpmltop,6), &
    psizpml(dim1(rank),dim2(rank),nzpmltop,6))

    if (USE_HPML .AND. IAM_XBOT) allocate(scrxpml(nxpmlbot,dim2(rank),dim3(rank),6),xpmlvars(nxpmlbot,dim2(rank),dim3(rank),6), &
    psixpml(nxpmlbot,dim2(rank),dim3(rank),6))   

    if (USE_HPML .AND. IAM_XTOP) allocate(scrxpml(nxpmltop,dim2(rank),dim3(rank),6),xpmlvars(nxpmltop,dim2(rank),dim3(rank),6), &
    psixpml(nxpmltop,dim2(rank),dim3(rank),6))

    if (USE_HPML .AND. IAM_YBOT) allocate(scrypml(dim1(rank),nypmlbot,dim3(rank),6),ypmlvars(dim1(rank),nypmlbot,dim3(rank),6), &
    psiypml(dim1(rank),nypmlbot,dim3(rank),6))

    if (USE_HPML .AND. IAM_YTOP) allocate(scrypml(dim1(rank),nypmltop,dim3(rank),6),ypmlvars(dim1(rank),nypmltop,dim3(rank),6), &
    psiypml(dim1(rank),nypmltop,dim3(rank),6))


  else
    if (.not. displ) then
      allocate(a(dim1(rank),dim2(rank),dim3(rank),5), temp_step(dim1(rank),dim2(rank),dim3(rank),5),&
      scr(dim1(rank),dim2(rank),dim3(rank),5))
    else
      allocate(a(dim1(rank),dim2(rank),dim3(rank),6), temp_step(dim1(rank),dim2(rank),dim3(rank),6),&
      scr(dim1(rank),dim2(rank),dim3(rank),6),scratch(dim1(rank),dim2(rank),dim3(rank),2))
    endif

    if (USE_PML .AND. IAM_ZBOT) allocate(scrzpml(dim1(rank),dim2(rank),nzpmlbot,2),zpmlvars(dim1(rank),dim2(rank),nzpmlbot,2), &
    psizpml(dim1(rank),dim2(rank),nzpmlbot,2))

    if (USE_PML .AND. IAM_ZTOP) allocate(scrzpml(dim1(rank),dim2(rank),nzpmltop,2),zpmlvars(dim1(rank),dim2(rank),nzpmltop,2), &
    psizpml(dim1(rank),dim2(rank),nzpmltop,2))

    if (USE_HPML .AND. IAM_XBOT) allocate(scrxpml(nxpmlbot,dim2(rank),dim3(rank),2),xpmlvars(nxpmlbot,dim2(rank),dim3(rank),2), &
    psixpml(nxpmlbot,dim2(rank),dim3(rank),2))

    if (USE_HPML .AND. IAM_XTOP) allocate(scrxpml(nxpmltop,dim2(rank),dim3(rank),2),xpmlvars(nxpmltop,dim2(rank),dim3(rank),2), &
    psixpml(nxpmltop,dim2(rank),dim3(rank),2))

    if (USE_HPML .AND. IAM_YBOT) allocate(scrypml(dim1(rank),nypmlbot,dim3(rank),2),ypmlvars(dim1(rank),nypmlbot,dim3(rank),2), &
    psiypml(dim1(rank),nypmlbot,dim3(rank),2))

    if (USE_HPML .AND. IAM_YTOP) allocate(scrypml(dim1(rank),nypmltop,dim3(rank),2),ypmlvars(dim1(rank),nypmltop,dim3(rank),2), &
    psiypml(dim1(rank),nypmltop,dim3(rank),2))


  endif
   
  allocate(gradp(dim1(rank),dim2(rank),dim3(rank),3),&
  dvzdz(dim1(rank),dim2(rank),dim3(rank)),c2(dim1(rank),dim2(rank),dim3(rank)), div(dim1(rank),dim2(rank),dim3(rank)),&
  dvxdx(dim1(rank),dim2(rank),dim3(rank)), dvydy(dim1(rank),dim2(rank),dim3(rank)),&
  spongexyz(dim1(rank), dim2(rank), dim3(rank)), gradrho0_z(dim1(rank), dim2(rank), dim3(rank)),&
  gradp0_z(dim1(rank), dim2(rank), dim3(rank)), rhoinv(dim1(rank), dim2(rank), dim3(rank)), rho0(dim1(rank), dim2(rank), dim3(rank)), &
  p0(dim1(rank), dim2(rank), dim3(rank)), c_speed(dim1(rank), dim2(rank), dim3(rank)), c2rho0(dim1(rank),dim2(rank),dim3(rank)),&
  LC0(dim1(rank),dim2(rank)), LC1(dim1(rank),dim2(rank)), LC2(dim1(rank),dim2(rank)),&
  LC3(dim1(rank),dim2(rank)), LC4(dim1(rank),dim2(rank)), LC5(dim1(rank),dim2(rank)), LC6(dim1(rank),dim2(rank)))

  if (FLOWS) then

    allocate(v0(dim1(rank),dim2(rank),dim3(rank),3))

    v0 = 0.0_dp
    v0_x => v0(:,:,:,1)
    v0_y => v0(:,:,:,2)
    v0_z => v0(:,:,:,3)

    if (.not. DISPL) then

    allocate(dv0x(dim1(rank),dim2(rank),dim3(rank),3),dv0y(dim1(rank),dim2(rank),dim3(rank),3),&
    dv0z(dim1(rank),dim2(rank),dim3(rank),3))

    dv0x = 0.0_dp    
    dv0y = 0.0_dp
    dv0z = 0.0_dp

    dv0xdx => dv0x(:,:,:,1)
    dv0xdy => dv0x(:,:,:,2)
    dv0xdz => dv0x(:,:,:,3)
    dv0ydx => dv0y(:,:,:,1)
    dv0ydy => dv0y(:,:,:,2)
    dv0ydz => dv0y(:,:,:,3)
    dv0zdx => dv0z(:,:,:,1)
    dv0zdy => dv0z(:,:,:,2)
    dv0zdz => dv0z(:,:,:,3)

    endif

  endif

  if (USE_PML .AND. IAM_ZBOT) allocate(az(dim1(rank),dim2(rank),nzpmlbot),bzpml(dim1(rank),dim2(rank),nzpmlbot), &
  kappaz(dim1(rank),dim2(rank),nzpmlbot))
  if (USE_PML .AND. IAM_ZTOP) allocate(az(dim1(rank),dim2(rank),nzpmltop),bzpml(dim1(rank),dim2(rank),nzpmltop), &
  kappaz(dim1(rank),dim2(rank),nzpmltop))

  if (USE_HPML .AND. IAM_XBOT) allocate(ax(nxpmlbot,dim2(rank),dim3(rank)),bxpml(nxpmlbot,dim2(rank),dim3(rank)), &
  kappax(nxpmlbot,dim2(rank),dim3(rank)))
  if (USE_HPML .AND. IAM_XTOP) allocate(ax(nxpmltop,dim2(rank),dim3(rank)),bxpml(nxpmltop,dim2(rank),dim3(rank)), &
  kappax(nxpmltop,dim2(rank),dim3(rank)))
  if (USE_HPML .AND. IAM_YBOT) allocate(ay(dim1(rank),nypmlbot,dim3(rank)),bypml(dim1(rank),nypmlbot,dim3(rank)), &
  kappay(dim1(rank),nypmlbot,dim3(rank)))
  if (USE_HPML .AND. IAM_YTOP) allocate(ay(dim1(rank),nypmltop,dim3(rank)),bypml(dim1(rank),nypmltop,dim3(rank)), &
  kappay(dim1(rank),nypmltop,dim3(rank)))

    rho0 = 0.0_dp
    p0 = 0.0_dp
  
    if (magnetic) then
      box = 0.0_dp
      boy = 0.0_dp
      boz = 0.0_dp
    endif

    call READ_IN_HDF5 ()

    if (pulse) then
    allocate(forcing_p(dim1(rank),dim2(rank),dim3(rank)),source_damp(dim1(rank),dim2(rank),dim3(rank)))
    source_damp = 0.0_dp
    else
    allocate(forcing(dim1(rank),dim2(rank)),source_damp(dim1(rank),dim2(rank),1))
    endif
 
  source_damp = 0.0_dp
  a = 0.0_dp
  scr = 0.0_dp
  temp_step = 0.0_dp

  call initexplicitders ()

  !! This needs to be fixed up when the derivative routines are fixed.

  !! send z's needed to zaxis

  if (IAM_ZBOT) then

    do ii = 1,g_z_der
        unstretch(ii) = sum(zaxis(1:2*g_z_der+1)*dpz(ii,:))*(nz-1)
    enddo

    do ii = g_z_der+1,dim3(rank)
 
        unstretch(ii) = ((-zaxis(zstart(rank)-1 + ii-5) + zaxis(zstart(rank)-1 + ii+5))*dpz(6,11) + &
        (-zaxis(zstart(rank)-1 + ii-4) + zaxis(zstart(rank)-1 + ii+4))*dpz(6,10) + &
        (-zaxis(zstart(rank)-1 + ii-3) + zaxis(zstart(rank)-1 + ii+3))*dpz(6,9) + &
        (-zaxis(zstart(rank)-1 + ii-2) + zaxis(zstart(rank)-1 + ii+2))*dpz(6,8) + &
        (-zaxis(zstart(rank)-1 + ii-1) + zaxis(zstart(rank)-1 + ii+1))*dpz(6,7))*(nz-1)
    enddo

  elseif (IAM_ZTOP) then

    do ii = 1,dim3(rank)-g_z_der
        unstretch(ii) = ((-zaxis(zstart(rank)-1 + ii-5) + zaxis(zstart(rank)-1 + ii+5))*dpz(6,11) + &
        (-zaxis(zstart(rank)-1 + ii-4) + zaxis(zstart(rank)-1 + ii+4))*dpz(6,10) + &
        (-zaxis(zstart(rank)-1 + ii-3) + zaxis(zstart(rank)-1 + ii+3))*dpz(6,9) + &
        (-zaxis(zstart(rank)-1 + ii-2) + zaxis(zstart(rank)-1 + ii+2))*dpz(6,8) + &
        (-zaxis(zstart(rank)-1 + ii-1) + zaxis(zstart(rank)-1 + ii+1))*dpz(6,7))*(nz-1)
    enddo

    do ii = 1,g_z_der
        jj = dim3(rank)+1-ii

        unstretch(jj) = -(zaxis(nz)*dpz(ii,1) + &
        zaxis(nz-1)*dpz(ii,2) + &
        zaxis(nz-2)*dpz(ii,3) + &
        zaxis(nz-3)*dpz(ii,4) + &
        zaxis(nz-4)*dpz(ii,5) + &
        zaxis(nz-5)*dpz(ii,6) + &
        zaxis(nz-6)*dpz(ii,7) + &
        zaxis(nz-7)*dpz(ii,8) + &
        zaxis(nz-8)*dpz(ii,9) + &
        zaxis(nz-9)*dpz(ii,10) + &
        zaxis(nz-10)*dpz(ii,11))*(nz-1)

    enddo

  else

    do ii = 1,dim3(rank)
        unstretch(ii) = ((-zaxis(zstart(rank)-1 + ii-5) + zaxis(zstart(rank)-1 + ii+5))*dpz(6,11) + &
        (-zaxis(zstart(rank)-1 + ii-4) + zaxis(zstart(rank)-1 + ii+4))*dpz(6,10) + &
        (-zaxis(zstart(rank)-1 + ii-3) + zaxis(zstart(rank)-1 + ii+3))*dpz(6,9) + &
        (-zaxis(zstart(rank)-1 + ii-2) + zaxis(zstart(rank)-1 + ii+2))*dpz(6,8) + &
        (-zaxis(zstart(rank)-1 + ii-1) + zaxis(zstart(rank)-1 + ii+1))*dpz(6,7))*(nz-1)
    enddo

  endif

  do kk=1,dim3(rank)
    stretch(kk) = (nz-1)*1_dp/unstretch(kk)
  enddo

   ! Non-dimensional form, for the calculations
   deltat = timestep*dimc/diml
 
  if (rank == 0) print *,'With cadence of 60 s, the timestep is', deltat*diml/dimc

  if (USE_PML .AND. (IAM_ZTOP .OR. IAM_ZBOT)) then
    psizpml = 0.0_dp
    scrzpml = 0.0_dp
    zpmlvars = 0.0_dp

    psizvz => psizpml(:,:,:,1)
    psizp => psizpml(:,:,:,2)

    RHSpsizvz => scrzpml(:,:,:,1)
    RHSpsizp => scrzpml(:,:,:,2)
  endif

  if (USE_HPML .AND. (IAM_XBOT .OR. IAM_XTOP)) then

    psixpml = 0.0_dp
    scrxpml = 0.0_dp
    xpmlvars = 0.0_dp

    psixvx => psixpml(:,:,:,1)
    psixp => psixpml(:,:,:,2)

    RHSpsixvx => scrxpml(:,:,:,1)
    RHSpsixp => scrxpml(:,:,:,2)

  endif
  if (USE_HPML .AND. (IAM_YBOT .OR. IAM_YTOP)) then

    psiypml = 0.0_dp
    scrypml = 0.0_dp
    ypmlvars = 0.0_dp

    psiyvy => psiypml(:,:,:,1)
    psiyp => psiypml(:,:,:,2)

    RHSpsiyvy => scrypml(:,:,:,1)
    RHSpsiyp => scrypml(:,:,:,2)

  endif

  gradp_x => gradp(:,:,:,1)
  gradp_y => gradp(:,:,:,2)
  gradp_z => gradp(:,:,:,3)

  if (.not. displ) then
    rho => temp_step(:,:,:,1)
    v_x => temp_step(:,:,:,2)
    v_y => temp_step(:,:,:,3)
    v_z => temp_step(:,:,:,4)
    p => temp_step(:,:,:,5)

    RHScont => scr(:,:,:,1) ! Continuity eqn.
    RHSv_x  => scr(:,:,:,2) ! Radial momentum
    RHSv_y  => scr(:,:,:,3) ! Latitudinal momentum
    RHSv_z  => scr(:,:,:,4) ! Longitudinal momentum
    RHSp  => scr(:,:,:,5) ! Pressure eqn.

  else
    xi_x => temp_step(:,:,:,1) ! Continuity eqn.
    xi_y  => temp_Step(:,:,:,2) ! Radial momentum
    xi_z  => temp_step(:,:,:,3) ! Latitudinal momentum
    v_x  => temp_step(:,:,:,4) ! Longitudinal momentum
    v_y  => temp_step(:,:,:,5) ! Pressure eqn.
    v_z  => temp_step(:,:,:,6) ! Pressure eqn.

    RHSxi_x => scr(:,:,:,1) ! Continuity eqn.
    RHSxi_y  => scr(:,:,:,2) ! Radial momentum
    RHSxi_z  => scr(:,:,:,3) ! Latitudinal momentum
    RHSv_x  => scr(:,:,:,4) ! Longitudinal momentum
    RHSv_y  => scr(:,:,:,5) ! Pressure eqn.
    RHSv_z  => scr(:,:,:,6) ! Pressure eqn.

    dxizdz => dvzdz
    dxixdx => dvxdx
    dxiydy => dvydy

    p => scratch(:,:,:,2)
    rho => scratch(:,:,:,1)

  endif

  if (magnetic) then 

    if (.not. displ) then
      bx => temp_step(:,:,:,6)
      by => temp_step(:,:,:,7)
      bz => temp_step(:,:,:,8)

      RHSb_x => scr(:,:,:,6) ! equation for b_x
      RHSb_y => scr(:,:,:,7) ! equation for b_y
      RHSb_z => scr(:,:,:,8) ! equation for b_z
    else
      bx => scratch(:,:,:,3)
      by => scratch(:,:,:,4)
      bz => scratch(:,:,:,5)
    endif

   if (USE_PML .AND. (IAM_ZTOP .OR. IAM_ZBOT)) then

    RHSpsizdzbx => scrzpml(:,:,:,3)
    RHSpsizdzby => scrzpml(:,:,:,4)
    RHSpsizinductionbx => scrzpml(:,:,:,5)
    RHSpsizinductionby => scrzpml(:,:,:,6)
 
    psizdzbx => psizpml(:,:,:,3)
    psizdzby => psizpml(:,:,:,4)
    psizinductionbx => psizpml(:,:,:,5) 
    psizinductionby => psizpml(:,:,:,6)
   
   endif

   if (USE_HPML .AND. (IAM_XBOT .OR. IAM_XTOP)) then

    RHSpsixdxby => scrxpml(:,:,:,3)
    RHSpsixdxbz => scrxpml(:,:,:,4)
    RHSpsixinductionby => scrxpml(:,:,:,5)
    RHSpsixinductionbz => scrxpml(:,:,:,6)

    psixdxby => psixpml(:,:,:,3)
    psixdxbz => psixpml(:,:,:,4)
    psixinductionby => psixpml(:,:,:,5)
    psixinductionbz => psixpml(:,:,:,6)
   endif

   if (USE_HPML .AND. (IAM_YBOT .OR. IAM_YTOP)) then

    RHSpsiydybx => scrypml(:,:,:,3)
    RHSpsiydybz => scrypml(:,:,:,4)
    RHSpsiyinductionbx => scrypml(:,:,:,5)
    RHSpsiyinductionbz => scrypml(:,:,:,6)

    psiydybx => psiypml(:,:,:,3)
    psiydybz => psiypml(:,:,:,4) 
    psiyinductionbx => psiypml(:,:,:,5)
    psiyinductionbz => psiypml(:,:,:,6)

   endif

  endif 

  if (DAMPING) call INIT_DAMPING()

  call MPI_BARRIER(comm_cart_3D, ierr)

end subroutine INITANDALLOCATE
!============================================================================

 subroutine INITIALISE_RHS()

  integer ii,jj,kk,bc,key
  real(dp),allocatable,dimension(:) :: alpha
  real(dp) start_sponge, finish_sponge,va,dtemp,temp,pml_l,pml_f0
  real(dp), allocatable, dimension(:,:,:) :: tempin, damp
  real(dp), allocatable, dimension(:) :: source_dampx, source_dampy
  real(dp), dimension(dim1(rank),dim2(rank),dim3(rank)) :: curltemp
  integer ierr
! Computing the absorbent boundary adjoining sponge. 
 
! First compute sponge_z, for the PML case this is a small absorbing sponge located over the PML CELLS
! For the non-PML case we have it as defined in the original sparc (nobody uses this anymore but)
! Realisticly the strength of the top sponge should be higher the higher in the chromosphere/corona we go,
! I am too lazy to do this now.

  if (periodic) then
    bc = 5
  else
    bc = 1
  endif

  if (.not. USE_PML) then
    do kk = 1,dim3(rank)
      spongez(kk) = 400.0_dp/(1_dp + exp((1.002_dp - z(kk))/0.0001_dp)) &
      + 40.0_dp/(1_dp + exp((z(kk) - (zaxis(1) + 0.005_dp))/0.0005_dp))
    enddo
  endif

! Compute x (and if 3D, y sponges)

! Damping length from the horizontal boundaries
  start_sponge = sp_length/xlength
  finish_sponge = 1_dp -  sp_length/xlength

  spongexyz =  0.0_dp
  spongex = 0.0_dp
  spongey = 0.0_dp

  if (PERIODIC) then
    do kk=1,dim3(rank)
      spongexyz(:,:,kk) =  spongez(kk)
    enddo
  else
    if (.not. USE_HPML) then
      do ii=1,dim1(rank)
        spongex(ii) = sp_str*(1_dp - 1_dp/( (1_dp + exp((start_sponge-x(ii))/sp_decay))*&
        (1_dp + exp((x(ii)-finish_sponge)/sp_decay)) ))
      enddo
      if (.not. TWOD) then
        do jj=1,dim2(rank)
          spongey(jj) = sp_str*(1_dp - 1_dp/( (1_dp + exp((start_sponge-y(jj))/sp_decay))&
          *(1_dp + exp((y(jj)-finish_sponge)/sp_decay)) ))
        enddo
      endif 
    endif

    do kk = 1,dim3(rank)
      do jj = 1,dim2(rank)
        do ii = 1,dim1(rank)
          spongexyz(ii,jj,kk) = spongex(ii) + spongey(jj)+spongez(kk)
        enddo
      enddo
    enddo
  endif

  call MPI_BARRIER(comm_cart_3D,ierr)

  if (rank == 0) print *,' Sponges initialized!'

!!! Reduction stuff for magnetic field runs, has been rewritten to use the max_VA in the 

  if (magnetic) then
   
    !! This little bit sets the reduction amount for the Alfven limiter, with
    !! the max VA set in the parameters file. There is probably a better way of
    !! doing it, but this works and it is in the initialisation routines so the
    !! speed probably doesn't matter
 
    do kk = 1,dim3(rank)
            reduction(:,:,kk) = (boz(:,:,kk)**2+ box(:,:,kk)**2+boy(:,:,kk)**2)/(rho0(:,:,kk)*cq(kk)**2)
    enddo

    va = 0.0_dp
    dtemp = 1.0e9_dp
    temp = -dtemp
 
    if (maxval(reduction) .GT. 0) then
      do while ((dtemp .GT. 1.0e-5_dp) .AND. (temp .LT. 1.0e10_dp-1.1_dp))
        do while ((va .LT. maxva) .AND. (temp .LT. 1.0e10_dp))
          temp = temp + dtemp
          va = maxval(sqrt(dimb**2*temp/(temp+reduction)*(box**2+boy**2+boz**2)/(4.0_dp*pi*dimrho*rho0)))
        enddo
        temp = temp-dtemp
        dtemp = dtemp/10.0_dp
        va = 0.0_dp
      enddo
    else
      temp = 1.0d10
    endif

    call MPI_ALLREDUCE(temp,dtemp,1,MPI_DOUBLE_PRECISION, MPI_MIN,comm_cart_3D,ierr)

    reduction = dtemp/(dtemp + reduction)

    call MPI_REDUCE(maxval(sqrt(dimb**2*reduction*(box**2+boy**2+boz**2)/(4.0_dp*pi*dimrho*rho0))),&
    va,1,MPI_DOUBLE_PRECISION, MPI_MAX,0,comm_cart_3D,ierr)

    if (rank ==0) print *, 'REDUCTION DONE: MaxVA =', va*1.0e-5_dp, 'km/s'

    curltemp = 0.0_dp
    call ddx(boy, curlboz, bc)
    if (.not. TWOD) call ddy(box, curltemp, bc)
    curlboz = curlboz - curltemp

    curltemp=0.0_dp
    if (.not. TWOD) call ddy(boz, curltemp, bc)
    call ddz(boy, curlbox, 1)
    curlbox=curltemp-curlbox

    curltemp=0.0_dp
    call ddz(box, curlboy, 1)
    call ddx(boz, curltemp, bc)
    curlboy=curlboy-curltemp

  endif

! Stabilize the background model dpdz, according to Schunker et al. (2011)

  rhoinv = 1.0/rho0
  c2rho0 = c2 * rho0

  call ddz(rho0, gradrho0_z, 1)
  call ddz(p0, gradp0_z, 1)

  if (STABILISE) then
    do kk=1,dim3(rank)
      do jj = 1,dim2(rank)
        do ii = 1,dim1(rank)
          ! FOR NEAR SURFACE LAYERS ONLY, (-0.15<z<0.1)
          if ((z(kk) .GT. 0.9998) .AND. (z(kk) .LT. 1.0000)) then           
            gradp0_z(ii,jj,kk)=max(c2(ii,jj,kk)*gradrho0_z(ii,jj,kk)+1.e-4_dp*Rsun/(dimrho*dimc**2),-rho0(ii,jj,kk)*g(kk))
          ! AND FOR THE REST OF THE BOX z<= -0.15 and z>=0.1 Mm
          elseif (z(kk) .LE. 0.9998) then
            gradp0_z(ii,jj,kk)=max(c2(ii,jj,kk)*gradrho0_z(ii,jj,kk),-rho0(ii,jj,kk)*g(kk))
          elseif (z(kk) .GE. 1.0001) then
            gradp0_z(ii,jj,kk)=max(c2(ii,jj,kk)*gradrho0_z(ii,jj,kk),-rho0(ii,jj,kk)*g(kk))
         endif
        enddo
      enddo
    enddo
  endif

  if (USE_PML) then
    allocate(tempin(dim1(rank),dim2(rank),dim3(rank)),source_dampx(dim1(rank)),source_dampy(dim2(rank)))
    source_dampx = 0.0_dp
    source_dampy = 0.0_dp
    if (magnetic) then
      tempin= sqrt(c2 + reduction*(box**2 + boy**2 + boz**2)/rho0)
    else
      tempin = sqrt(c2)
    endif
    
    pml_f0 = pi*0.005_dp*diml/dimc

    if (IAM_ZBOT) then

      allocate(damp(dim1(rank),dim2(rank),nzpmlbot),alpha(nzpmlbot))

      pml_L = (zaxis(nzpmlbot) - zaxis(1))

      do kk = 1,nzpmlbot
        damp(:,:,kk) = zpmlb(1)*tempin(:,:,nzpmlbot)*(zpmlb(2)+1_dp)/(2_dp*pml_L)* &
        ((zaxis(nzpmlbot)-z(kk))/pml_L)**zpmlb(2)
        alpha(kk) = pml_f0*(z(kk)-zaxis(1))/pml_L
        kappaz(:,:,kk) = 1_dp+(zpmlb(3)-1_dp)*((zaxis(nzpmlbot)-z(kk))/pml_L)**zpmlb(2)

        bzpml(:,:,kk) = - alpha(kk) - damp(:,:,kk)/kappaz(:,:,kk)
        az(:,:,kk) = -damp(:,:,kk)/kappaz(:,:,kk)**2

        gradrho0_z(:,:,kk) = gradrho0_z(:,:,kk)*&
        alpha(kk)/kappaz(:,:,kk)/(damp(:,:,kk)/kappaz(:,:,kk) + alpha(kk))
        gradp0_z(:,:,kk) = gradp0_z(:,:,kk)*&
        alpha(kk)/kappaz(:,:,kk)/(damp(:,:,kk)/kappaz(:,:,kk) + alpha(kk))
        g(kk)  =g(kk) * alpha(kk)/kappaz(1,1,kk)/(damp(1,1,kk)/kappaz(1,1,kk) + alpha(kk))
        spongexyz(:,:,kk) = spongexyz(:,:,kk) + zpmlb(4)*tempin(:,:,nzpmlbot)* &
        ((zaxis(nzpmlbot)-z(kk))/pml_L)**(zpmlb(2)+1)

      enddo
      deallocate(damp,alpha)

    endif
    if (IAM_ZTOP) then

      allocate(damp(dim1(rank),dim2(rank),nzpmltop),alpha(nzpmltop))

      pml_L = (zaxis(nz) - zaxis(nz-nzpmltop+1))

      do kk = 1,nzpmltop
        damp(:,:,kk) = zpmlt(1)*tempin(:,:,dim3(rank)-nzpmltop+1)*(zpmlt(2)+1_dp)/(2_dp*pml_L)* &
        (-(z(dim3(rank)-nzpmltop +1)-z(dim3(rank)-nzpmltop+kk))/pml_L)**zpmlt(2)
        alpha(kk) = pml_f0*(z(dim3(rank)) -z(dim3(rank)-nzpmltop+kk))/pml_L
        kappaz(:,:,kk) = 1_dp+(zpmlt(3)-1_dp)*(-(z(dim3(rank)-nzpmltop +1) &
-z(dim3(rank)-nzpmltop+kk))/pml_L)**zpmlt(2)

        bzpml(:,:,kk) = - alpha(kk) - damp(:,:,kk)/kappaz(:,:,kk)
        az(:,:,kk) = -damp(:,:,kk)/kappaz(:,:,kk)**2.0

        gradrho0_z(:,:,dim3(rank)-nzpmltop+kk) = gradrho0_z(:,:,dim3(rank)-nzpmltop+kk)*&
alpha(kk)/kappaz(:,:,kk)/(damp(:,:,kk)/kappaz(:,:,kk) +  alpha(kk))
        gradp0_z(:,:,dim3(rank)-nzpmltop+kk) = gradp0_z(:,:,dim3(rank)-nzpmltop+kk)*&
alpha(kk)/kappaz(:,:,kk)/(damp(:,:,kk)/kappaz(:,:,kk) +  alpha(kk))
        g(dim3(rank)-nzpmltop+kk)  =g(dim3(rank)-nzpmltop+kk) * &
alpha(kk)/kappaz(1,1,kk)/(damp(1,1,kk)/kappaz(1,1,kk) +  alpha(kk))

        spongexyz(:,:,dim3(rank)-nzpmltop+kk) = spongexyz(:,:,dim3(rank)-nzpmltop+kk) + zpmlt(4)*tempin(:,:,dim3(rank)-nzpmltop+1)* &
        (-(z(dim3(rank)-nzpmltop +1) -z(dim3(rank)-nzpmltop+kk))/pml_L)**(zpmlt(2)+1)

      enddo
      deallocate(damp,alpha)

    endif

    if (USE_HPML) then

    if (IAM_XBOT) then

      allocate(damp(nxpmlbot,dim2(rank),dim3(rank)),alpha(nxpmlbot))

      pml_L = (xaxis(nxpmlbot) - xaxis(1))

      do ii = 1,nxpmlbot
        damp(ii,:,:) = xpmlb(1)*tempin(nxpmlbot,:,:)*(xpmlb(2)+1_dp)/(2_dp*pml_L)* &
        ((xaxis(nxpmlbot)-x(ii))/pml_L)**xpmlb(2)
        alpha(ii) = pml_f0*(x(ii)-xaxis(1))/pml_L
        kappax(ii,:,:) = 1_dp+(xpmlb(3)-1_dp)*((xaxis(nxpmlbot)-x(ii))/pml_L)**xpmlb(2)

        bxpml(ii,:,:) = - alpha(ii) - damp(ii,:,:)/kappax(ii,:,:)
        ax(ii,:,:) = -damp(ii,:,:)/kappax(ii,:,:)**2

        source_dampx(ii) = ((xaxis(nxpmlbot)-x(ii))/pml_L)**xpmlb(2)
        spongexyz(ii,:,:) = spongexyz(ii,:,:) + xpmlb(4)*tempin(nxpmlbot,:,:)* &
        ((xaxis(nxpmlbot)-x(ii))/pml_L)**(xpmlb(2)+1)
      enddo
      deallocate(damp,alpha)

   endif
   if (IAM_XTOP) then

      allocate(damp(nxpmltop,dim2(rank),dim3(rank)),alpha(nxpmltop))

      pml_L = (xaxis(nx) - xaxis(nx-nxpmltop+1))

      do ii = 1,nxpmltop
        damp(ii,:,:) = xpmlt(1)*tempin(dim1(rank)-nxpmltop+1,:,:)*(xpmlt(2)+1_dp)/(2_dp*pml_L)* &
        (-(xaxis(nx-nxpmltop+1) -xaxis(nx-nxpmltop+ii))/pml_L)**xpmlt(2)
        alpha(ii) = pml_f0*(xaxis(nx) -xaxis(nx-nxpmltop+ii))/pml_L 
        kappax(ii,:,:) = 1_dp+(xpmlt(3)-1_dp)*(-(xaxis(nx-nxpmltop+1) &
-x(dim1(rank)-nxpmltop+ii))/pml_L)**xpmlt(2)

        source_dampx(dim1(rank)-nxpmltop+ii) = (-(xaxis(nx-nxpmltop+1)-xaxis(nx-nxpmltop+ii))/pml_L)**xpmlt(2)

        bxpml(ii,:,:) = - alpha(ii) - damp(ii,:,:)/kappax(ii,:,:)
        ax(ii,:,:) = -damp(ii,:,:)/kappax(ii,:,:)**2
        spongexyz(dim1(rank)-nxpmltop+ii,:,:) = spongexyz(dim1(rank)-nxpmltop+ii,:,:) + xpmlt(4)*tempin(dim1(rank)-nxpmltop+1,:,:)* &
        (-(xaxis(nx-nxpmltop+1) -xaxis(nx-nxpmltop+ii))/pml_L)**(xpmlt(2)+1)
      enddo
      deallocate(damp,alpha)

    endif

    if (IAM_YBOT) then

      allocate(damp(dim1(rank),nypmlbot,dim3(rank)), alpha(nypmlbot))

      pml_L = (yaxis(nypmlbot) - yaxis(1))

      do jj = 1,nypmlbot
        damp(:,jj,:) = ypmlb(1)*tempin(:,nypmlbot,:)*(ypmlb(2)+1_dp)/(2_dp*pml_L)* &
        ((yaxis(nypmlbot)-yaxis(jj))/pml_L)**ypmlb(2)
        alpha(jj) = pml_f0*(yaxis(jj)-yaxis(1))/pml_L
        kappay(:,jj,:) = 1_dp+(ypmlb(3)-1_dp)*((yaxis(nypmlbot)-yaxis(jj))/pml_L)**ypmlb(2)
 
        source_dampy(jj) = ((yaxis(nypmlbot)-yaxis(jj))/pml_L)**ypmlb(2) 
       
        bypml(:,jj,:) = - alpha(jj) - damp(:,jj,:)/kappay(:,jj,:)
        ay(:,jj,:) = -damp(:,jj,:)/kappay(:,jj,:)**2
        spongexyz(:,jj,:) = spongexyz(:,jj,:) + ypmlb(4)*tempin(:,nypmlbot,:)* &
        ((yaxis(nypmlbot)-yaxis(jj))/pml_L)**(ypmlb(2)+1)
      enddo

      deallocate(damp,alpha)

    endif
    if (IAM_YTOP) then

      allocate(damp(dim1(rank),nypmltop,dim3(rank)), alpha(nypmltop))

      pml_L = (yaxis(ny) - yaxis(ny-nypmltop+1))

      do jj = 1,nypmltop
        damp(:,jj,:) = ypmlt(1)*tempin(:,dim2(rank)-nypmltop+1,:)*(ypmlt(2)+1_dp)/(2_dp*pml_L)* &
        (-(yaxis(ny-nypmltop+1) -yaxis(ny-nypmltop+jj))/pml_L)**ypmlt(2)
        alpha(jj) = pml_f0*(yaxis(ny) -yaxis(ny-nypmltop+jj))/pml_L
        kappay(:,jj,:) = 1_dp+(ypmlt(3)-1_dp)*((-(yaxis(ny-nypmltop+1)-yaxis(ny-nypmltop+jj)))/pml_L)**ypmlt(2)

        source_dampy(dim2(rank)-nypmltop+jj) = (-(yaxis(ny-nypmltop+1) -yaxis(ny-nypmltop+jj))/pml_L)**ypmlt(2)

        bypml(:,jj,:) = - alpha(jj) - damp(:,jj,:)/kappay(:,jj,:)
        ay(:,jj,:) = -damp(:,jj,:)/kappay(:,jj,:)**2
        spongexyz(:,dim2(rank)-nypmltop+jj,:) = spongexyz(:,dim2(rank)-nypmltop+jj,:) + ypmlt(4)*tempin(:,dim2(rank)-nypmltop+1,:)* &
        (-(yaxis(ny-nypmltop+1) -yaxis(ny-nypmltop+jj))/pml_L)**(ypmlt(2)+1)

      enddo
      deallocate(damp,alpha)

    endif

    endif

    do ii = 1,dim1(rank)
      do jj = 1,dim2(rank)
        source_damp(ii,jj,:) = max(source_dampx(ii),source_dampy(jj))
      enddo
    enddo

    deallocate(tempin,source_dampx,source_dampy)
  endif

  if (PERIODIC) then
    source_damp = 1_dp
  else
    source_damp = 1_dp - source_damp
  endif


  if (magnetic) then
    
    call ddx(rho0, gradrho0_x, bc)
    call ddx(p0, gradp0_x, bc)

    if (.not. TWOD) then
      call ddy(rho0, gradrho0_y, bc)
      call ddy(p0, gradp0_y, bc)
    endif
  endif


  !! Chunk Output Stuff
  IAM_CHUNKOUT = .FALSE.

  do ii = 0,dim1(rank)-1
    do jj = 0,dim2(rank)-1
      do kk = 0,dim3(rank)-1
        if ((xstart(rank) + ii .GE. xcb) .AND. (xstart(rank)+ii .LE. xct) .AND. (ystart(rank)+jj .GE. ycb) &
        .AND.(ystart(rank)+jj .LE. yct) .AND. (zstart(rank)+kk .GE. zcb) .AND.(zstart(rank) +kk .LE. zct)) then
          IAM_CHUNKOUT = .TRUE.
        endif
      enddo
    enddo
  enddo
  key = 0
  if (IAM_CHUNKOUT) then
      key = 1
      my_xcb = max(xcb-xstart(rank)+1,1)
      my_xct = min(xct-xstart(rank)+1,dim1(rank))
      my_ycb = max(ycb-ystart(rank)+1,1)
      my_yct = min(yct-ystart(rank)+1,dim2(rank))
      my_zcb = max(zcb-zstart(rank)+1, 1)
      my_zct = min(zct- zstart(rank)+1,dim3(rank))

      my_offsetx = xstart(rank)+my_xcb-xcb-1
      my_offsety = ystart(rank)+my_ycb-ycb-1
      my_offsetz = zstart(rank)+my_zcb-zcb-1

      my_nxc = my_xct-my_xcb+1
      my_nyc = my_yct-my_ycb+1
      my_nzc = my_zct-my_zcb+1

   endif

  call MPI_COMM_SPLIT(comm_cart_3D,key,rank,comm_cart_savechunk,ierr)

  call MPI_BARRIER(comm_cart_3D,ierr)

   if (rank == 0) then
     print *,'Source excitation at radius', zaxis(e_rad)
     print *,'The corresponding radial gridpoint', e_rad
   endif
 
 end subroutine INITIALISE_RHS
!==========================================================================

subroutine INIT_DAMPING()

 integer ii,jj
 integer narr(2),inembed(2),onembed(2)
 real*8 constx, consty, kayref

 if (dim_damp(rank) .GT. 0) then

 allocate(transfermatrix(nx,ny,dim_damp(rank)),transfertemp(nx/2+1,ny,dim_damp(rank)), &
 damping_rates(nx/2+1,ny),kay2d(nx/2+1,ny))

 constx = 2_dp*pi*Rsun/xlength
 consty = 2_dp*pi*Rsun/ylength

 if (TWOD) then

 jj =1 
 do ii=1,nx/2+1
   kay2d(ii,jj) = (constx**2*(ii-1.0)**2 + consty**2*(jj-1.0)**2)**0.5_dp
  enddo

 else

 do jj=1,ny/2+1
  do ii=1,nx/2+1
   kay2d(ii,jj) = (constx**2*(ii-1.0)**2 + consty**2*(jj-1.0)**2)**0.5_dp
  enddo
 enddo
 
do jj=ny/2+2,ny
  do ii=1,nx/2+1
   kay2d(ii,jj) = (constx**2*(ii-1.0)**2 + consty**2*(ny-jj+1.0)**2)**0.5_dp
  enddo
 enddo
 endif

 kayref = 902_dp

 damping_rates = diml/dimc/nx/ny*1.0e-6_dp*(500.0_dp*(kay2d/kayref)**2.2_dp)*20.0_dp

  narr(1) = nx
  narr(2) = ny
  inembed(1) = nx
  inembed(2) = ny
  onembed(1) = nx/2+1
  onembed(2) = ny
  call dfftw_plan_many_dft_r2c(fftw_plan_fwd_2D,2,narr,dim_damp(rank),transfermatrix(1,1,1), &
inembed,1,nx*ny,transfertemp(1,1,1),onembed,1,(nx/2+1)*ny,FFTW_MEASURE)

   narr(1) = nx
   narr(2) = ny
   inembed(1) = nx/2+1
   inembed(2) = ny
   onembed(1) = nx
   onembed(2) = ny

   call dfftw_plan_many_dft_c2r(fftw_plan_inv_2D,2,narr,dim_damp(rank),transfertemp(1,1,1), &
inembed,1,(nx/2+1)*ny,transfermatrix(1,1,1),onembed,1,nx*ny,FFTW_MEASURE)

   transfermatrix = 0.0_dp
   transfertemp = 0.0_dp
 
   endif

end subroutine INIT_DAMPING
!==========================================================================

subroutine PARSE_OPTIONS()

 if (.not. magnetic) then
   if (.not. displ) then
     if (.not. TWOD)  option = 1
     if (TWOD) option = 5
   else
     if (.not. TWOD)  option = 2
     if (TWOD) option = 6
   endif
 else
   if (.not. displ) then
     if (.not. TWOD)  option = 3
     if (TWOD) option = 7
   else
     if (.not. TWOD)  option = 4
     if (TWOD) option = 8
   endif
 endif

end subroutine PARSE_OPTIONS

!==========================================================================
  subroutine INITIALIZE_STEP

     optimals(1) =  1_dp
     optimals(2) =  0.5_dp
     optimals(3) =  0.1665579725791151184_dp
     optimals(4) =  0.03950410250586729987_dp
     optimals(5) =  0.007810706393337838236_dp

     betas(5) = optimals(2)
     betas(4) = optimals(3)/betas(5)
     betas(3) = optimals(4)/(betas(5)*betas(4))
     betas(2) = optimals(5)/(betas(5)*betas(4)*betas(3))
     betas(1) = 0.0

     betas = betas * deltat
  end subroutine INITIALIZE_STEP
!=====================================================================================


subroutine READ_SPARC_PARAMETERS()

    !! Message contains a buffer to be used to repeat input settings in the log
    !! file. Filename reads the input filename from the first command arguement when executed.

    character(len=140) message,filename

    !! Line Buffer includes new line, if it contains an equals, linecase looks for the name to
    !! the left of the =, line data takes the value to the right of the
    !! equal to be read into the variable.

    CHARACTER(LEN=200) :: linebuffer, linevalue
    CHARACTER(LEN=30)  :: linename

    !! If line contains an = toread will be flagged, when EOF is reached
    !! linestate will become non-zero and the routine will end. Switch is used
    !! when a true/false or other is needed. Temp is just a temporary array
    !! which is sometimes needed.

    integer :: linestat, switch, temp,toread

    !_____________________________________________________________________

    linestat = 0

    call get_command_argument(1,filename)

    if (rank == 0) print *, "Reading input parameters from file:", trim(filename)
    
    open(0, file=trim(filename), action='READ')

    do while (linestat == 0)
       read(0, '(A)', iostat=linestat) linebuffer
       if (linestat == 0) then

          toread = SCAN(linebuffer, '=')

          if (toread .EQ. 0) then
             linename = ""
          else
             linename = trim(linebuffer(1:toread-1))
          endif
          linevalue = linebuffer(toread+1:)

          select case (linename)

          case('dir_out')
             read(linevalue, *, iostat=linestat) dir_out

          case('dir_bkg')
             read(linevalue, *, iostat=linestat) dir_bkg

          case('simname')
             read(linevalue, *, iostat=linestat) simname

          case('nx')
             read(linevalue, *, iostat=linestat) nx
          case('ny')
             read(linevalue, *, iostat=linestat) ny
          case('nz')
             read(linevalue, *, iostat=linestat) nz

          case('xlength')
             read(linevalue, *, iostat=linestat) xlength

          case('ylength')
             read(linevalue, *, iostat=linestat) ylength

          case('restart')
             read(linevalue, *,iostat=linestat) switch
             if (switch ==0) then
               restart = .false.
             else
               restart = .true.
             endif

          case('time0')
             read(linevalue, *, iostat=linestat) time0

          case('timestep')
             read(linevalue, *, iostat=linestat) timestep

          case('wall_time')
             read(linevalue, *, iostat=linestat) wall_time
             wall_time = wall_time*3600

          case('solartime')
             read(linevalue, *, iostat=linestat) solartime

          case('chunkoutputcad')
             read(linevalue, *, iostat=linestat) chunkoutputcad

          case('xcb')
             read(linevalue, *, iostat=linestat) xcb

          case('xct')
             read(linevalue, *, iostat=linestat) xct

          case('ycb')
             read(linevalue, *, iostat=linestat) ycb

          case('yct')
             read(linevalue, *, iostat=linestat) yct

          case('zcb')
             read(linevalue, *, iostat=linestat) zcb

          case('zct')
             read(linevalue, *, iostat=linestat) zct

          case('fulloutputcad')
             read(linevalue, *, iostat=linestat) fulloutputcad

          case('minsavtime')
             read(linevalue, *, iostat=linestat) minsavtime

          case('timestamp_size')
             read(linevalue, *, iostat=linestat) timestamp_size

          case('STABILISE')
             read(linevalue, *,iostat=linestat) switch
             if (switch ==0) then
               STABILISE = .false.
             else
               STABILISE = .true.
             endif 

          case('DAMPING')
             read(linevalue, *,iostat=linestat) switch
             if (switch ==0) then
               DAMPING = .false.
             else
               DAMPING = .true.
             endif


          case('MODE')
            read(linevalue, *,  iostat=linestat) switch
            if (switch == 0) then
              DISPL = .false.
            else 
              DISPL = .true.
            endif

            case('FLOWS')
              read(linevalue,*, iostat=linestat) switch
              if (switch == 0) then
                FLOWS = .false.
              else
                FLOWS = .true.
              endif

          case('PULSE')
             read(linevalue, *,iostat=linestat) switch
             if (switch ==0) then
               PULSE = .false.
             else  
               PULSE = .true.
               if (switch == 2) then
                 PSRPULSE = .true.
               endif
             endif

          case('filter_mode_xy')
             read(linevalue, *,iostat=linestat) switch
             if (switch ==1) then
               g_xy_filt = 5
             else if (switch == 2) then
               g_xy_filt = 3
             endif

          case('filter_mode_z')
             read(linevalue, *,iostat=linestat) switch
             if (switch == 1) then
               g_z_filt = 5
             else if (switch == 2) then
               g_z_filt = 3
             endif


          case('MAGNETIC')
             read(linevalue, *,iostat=linestat) switch
             if (switch ==0) then
               MAGNETIC = .false.
               QUIET = .true.
             elseif (switch == 1) then
               MAGNETIC = .false.
               QUIET = .false.
             elseif (switch == 2) then
               MAGNETIC = .true.
               QUIET = .false.
             endif

          case('maxva')
            read(linevalue, *,iostat=linestat) maxva

          case('TWOD')
             read(linevalue, *,iostat=linestat) switch
             if (switch ==0) then
               TWOD = .false.
             else
               TWOD = .true.
             endif

          case('TWOP5D')
             read(linevalue, *,iostat=linestat) TWOP5D

          case('USE_PML')
             read(linevalue, *,iostat=linestat) switch
             if (switch ==0) then 
               USE_PML = .false.
             else
               USE_PML = .true.
             endif

          case('USE_HPML')
            read(linevalue, *,iostat=linestat) switch
            if (switch ==0) then
              if (.not. USE_PML) print *,'USE_PML must be on for HPMLs, turning on USE_PML .true.'
              USE_PML = .true.
              USE_HPML = .false.
            else
              USE_HPML = .true.
            endif

          case('PERIODIC')
             read(linevalue, *,iostat=linestat) switch
             if (switch ==0) then
                 PERIODIC = .false.
             else
               if (USE_HPML) then
                 PERIODIC = .false.
                 if (rank ==0) then
                    print *, "conflict in parameters, both periodic and horizontal pmls are on"
                    print *, "turning periodic off"
                 endif
               else
                 PERIODIC = .true.
               endif
             endif

          case('xpmlt')
            read(linevalue, *, iostat=linestat) nxpmltop,xpmlt
            if (.not. USE_HPML) then
              nxpmltop = 0
            endif

          case('xpmlb')
            read(linevalue, *, iostat=linestat) nxpmlbot, xpmlb
            if (.not. USE_HPML) then
              nxpmlbot = 0
            endif

          case('ypmlt')
            read(linevalue, *, iostat=linestat) nypmltop,ypmlt
            if (.not. USE_HPML) then
              nypmltop = 0
            endif

          case('ypmlb')
            read(linevalue, *, iostat=linestat) nypmlbot,ypmlb
            if (.not. USE_HPML) then
              nypmlbot = 0
            endif

          case('zpmlt')
            read(linevalue, *, iostat=linestat) nzpmltop, zpmlt
            if (.not. USE_PML) then
              nzpmltop = 0
            endif

          case('zpmlb')
            read(linevalue, *, iostat=linestat) nzpmlbot,zpmlb
            if (.not. USE_PML) then
              nzpmlbot = 0
            endif

          case('sp_str')
             read(linevalue, *, iostat=linestat) sp_str

          case('sp_length')
             read(linevalue, *, iostat=linestat) sp_length

          case('sp_decay')
             read(linevalue, *, iostat=linestat) sp_decay

          case('forcingfunc')
             read(linevalue, *, iostat=linestat) forcingfunc

          case('cadforcing')
             read(linevalue, *, iostat=linestat) cadforcing

          case('randomexcitedepth')
             read(linevalue, *, iostat=linestat) randomexcitedepth

          case('randomexcitesigma')
             read(linevalue, *, iostat=linestat) randomexcitesigma

          case('randomexciteampl')
             read(linevalue, *, iostat=linestat) randomexciteampl

          case('pulse_t0')
             if (pulse) then 
             read(linevalue, *, iostat=linestat) pulse_t0
             else
             pulse_t0 = 1.0
             endif

          case('pulse_dir')
            if (pulse) then
              read(linevalue,*,iostat=linestat) pulse_dir
            else
              pulse_dir = 1
            endif

          case('pulse_t1')
             if (pulse) then 
             read(linevalue, *, iostat=linestat) pulse_t1
             else
             pulse_t1 = 1.0
             endif

          case('pulse_st')
             if (pulse) then 
             read(linevalue, *, iostat=linestat) pulse_st
             else
             pulse_st = 1.0
             endif

          case('pulse_px')
             if (pulse) then 
             read(linevalue, *, iostat=linestat) pulse_px
             else
             pulse_px = 1.0
             endif

          case('pulse_py')
             if (pulse) then 
             read(linevalue, *, iostat=linestat) pulse_py
             else
             pulse_py = 1.0
             endif

          case('pulse_pz')
            if (pulse) then
              read(linevalue,*,iostat=linestat) pulse_pz
            else
              pulse_pz = 1_dp
            endif

          case('pulse_sr')
             if (pulse) then 
             read(linevalue, *, iostat=linestat) pulse_sr
             else
             pulse_sr = 1.0
             endif

          case('pulse_sz')
             if (pulse) then 
             read(linevalue, *, iostat=linestat) pulse_sz
             else
             pulse_sz = 1.0
             endif

          case('pulse_amp')
             if (pulse) then 
             read(linevalue, *, iostat=linestat) pulse_amp
             else
             pulse_amp = 0.0
             endif

          case('s_xy')
             read(linevalue, *, iostat=linestat) s_xy

          case('s_z')
             read(linevalue, *, iostat=linestat) s_z

          case('t_xy')
             read(linevalue, *, iostat=linestat) t_xy
             if (t_xy .LT. timestep) t_xy = timestep

          case('t_z')
             read(linevalue, *, iostat=linestat) t_z
             if (t_z .LT. timestep) t_z = timestep

          case('block_x')
             read(linevalue,*, iostat=linestat) temp
             block_size(1) = temp

          case('block_y')
             read(linevalue,*, iostat=linestat) temp
             block_size(2) = temp

          case('block_z')
             read(linevalue,*, iostat=linestat) temp
             block_size(3) = temp

          endselect
       endif
    
    enddo

    if (TWOD) then
      ylength = 0.0_dp
      ny = 1
      nypmltop = 0
      nypmlbot = 0
      pulse_py = 0.0_dp
    endif

    if (TWOD) then
      simname = trim(simname)//'_2D'
    else
      simname = trim(simname)//'_3D'
    endif

    if (DISPL) then
      simname = trim(simname)//'_d'
    else
      simname = trim(simname)//'_v'
    endif
    if (MAGNETIC) then
      simname = trim(simname)//'m'
    elseif (QUIET) then
      simname = trim(simname)//'q'
    else
      simname = trim(simname)//'t'
    endif

    if (FLOWS) then
      simname = trim(simname)//'f'
    endif

    simname = trim(simname)//'_'
    nzpml = nzpmlbot+nzpmltop
    nxpml = nxpmlbot+nxpmltop
    nypml = nypmlbot+nypmltop

    if (TWOD) then
      ylength = 0.0_dp
      ny = 1
      nypmltop = 0
      nypmlbot = 0
      nypml = 0
    endif

    !! Do I have data in the chunk being output

    nxc = xct-xcb+1
    nyc = yct-ycb+1
    nzc = zct-zcb+1

    !! temporary

       if (rank == 0) then
          write(message,"(A9, I7.1,A5, I7.1,A5, I7.1)") 'Grid Nx: ',nx,' Ny: ',ny, ' Nz: ',nz
          print *,trim(message)
          PRINT *, "Box x-y length (Mm)", xlength/1.d8, ylength/1.d8
          print *, "3D Fits Cubes: ", dir_bkg
          print *, "Output: ", dir_out
          PRINT *, "Timestep: ", timestep," for a length of: ", solartime
          PRINT *, "Output every ", fulloutputcad, " starting at ", minsavtime, &
          " with a timestamp size of ", timestamp_size, " and a max walltime of", &
          wall_time/3600
          PRINT *, "CONFIG SETTINGS:"
          PRINT *, "STABILISE:", STABILISE
          PRINT *, "PULSE:", PULSE
          PRINT *, "MAGNETIC:", MAGNETIC, " and QUIET:", QUIET
          PRINT *, "2D?:", TWOD, " and extended to 2.5 D?", TWOP5D
          print *, "Boundary conditions:"
          if (USE_PML) then
            print *, "Z Pmls: ", USE_PML
            print *,  "Top PML, Size:", nzpmltop, " pml_RC:", zpmlt(1)," pml_N:", zpmlt(2)," pml_kappa:", zpmlt(3)," Sponge Sigma:", zpmlt(4)
            print *,  "Bot PML, Size:", nzpmlbot, " pml_RC:", zpmlb(1)," pml_N:", zpmlb(2)," pml_kappa:", zpmlb(3)," Sponge Sigma:", zpmlb(4) 
            if (USE_HPML) then
            print *, "X Pmls: ", USE_HPML 
            print *,  "Top PML, Size:", nxpmltop, " pml_RC:", xpmlt(1)," pml_N:", xpmlt(2)," pml_kappa:", xpmlt(3)," Sponge Sigma:", xpmlt(4) 
            print *,  "Bot PML, Size:", nxpmlbot, " pml_RC:", xpmlb(1)," pml_N:", xpmlb(2)," pml_kappa:", xpmlb(3)," Sponge Sigma:", xpmlb(4)    
            print *, "Y Pmls: ", USE_HPML 
            print *,  "Top PML, Size:", nypmltop, " pml_RC:", ypmlt(1)," pml_N:", ypmlt(2)," pml_kappa:", ypmlt(3)," Sponge Sigma:", ypmlt(4) 
            print *,  "Bot PML, Size:", nypmlbot, " pml_RC:", ypmlb(1)," pml_N:", ypmlb(2)," pml_kappa:", ypmlb(3)," Sponge Sigma:", ypmlb(4)    
            else
              print *, "XY Pmls: ", USE_HPML
            endif
          else
            print *, "Pmls: ", USE_PML
          endif
          if ((.NOT. PERIODIC) .AND. (.NOT. USE_HPML)) then
            print *, "Horizontal Sponges:"
            print *, "Strength: ", sp_str, " Length (Mm):", sp_length/1.d8, " Decay Sigma: ", sp_decay
          endif
          if (.not. PULSE) then
          print *, "Forcing Function:  ", forcingfunc
          print *, "at Cadence: ", cadforcing, " Amplitude: ", randomexciteampl
          print *, "Depth:", randomexcitedepth, " and Sigma: ", randomexcitesigma
          else
            print *, "pulse parameters:"
            print *, "t0: ", pulse_t0, " t1:", pulse_t1, "st: ", pulse_st
            print *, "px: ", pulse_px/1.d8, " py:", pulse_py/1.d8, "pz:",pulse_pz/1.d8
            print *, "sr: ", pulse_sr/1.d8, " sz:", pulse_sz/1.d8, "amplitude: ", pulse_amp
          endif
  
            print *, "Filtering, t_xy = ", t_xy, "s_xy", s_xy, "t_z:", t_z
       endif

    CLOSE(0)

   ! Length of simulation

   maxtime = floor(solartime*3600.0/timestep)
   maxtime = 2*(maxtime/2)

   fulloutputsteps = floor(fulloutputcad/timestep + 1.0e-5_dp)
   chunkoutputsteps = floor(chunkoutputcad/timestep + 1.0e-5_dp)

   CHUNKOUT = .TRUE.
   FULLOUT = .TRUE.

   if (chunkoutputsteps == 0) then
     CHUNKOUT=.FALSE.
     chunkoutputsteps=1
   endif
   if (fulloutputsteps == 0) then
     FULLOUT=.FALSE.
     fulloutputsteps=1
   endif
   xyfreq_filtering = floor(t_xy/timestep + 1.0e-5_dp)
   zfreq_filtering = floor(t_z/timestep + 1.0e-5_dp)
   if (rank ==0) then
     print *,'Chunk output:', CHUNKOUT, ' at a frequency of ', chunkoutputsteps, ' timesteps'
     print *,'Full output:', FULLOUT, ' at a frequency of ', fulloutputsteps, ' timesteps'
     print *, 'NUMBER OF TIMESTEPS =', maxtime
     print *, 'TIMESTEP =', timestep
       print *, 'FREQ OF FILTERING', xyfreq_filtering,zfreq_filtering
   endif

  num_steps = FLOOR(DBLE(maxtime)*timestep/cadforcing) + 2
  num_steps = 2*(num_steps/2)+1

  end subroutine READ_SPARC_PARAMETERS

  !==========================================================================
subroutine initexplicitders

 real(dp) :: coeffs(-3:3),coeffs2(-5:5)
 
 !! Coefficients for filter 1
 coeffs(-3) = -0.015625000_dp
 coeffs(-2) = 0.093750000_dp
 coeffs(-1) = -0.234375000_dp
 coeffs(0) = 0.312500000_dp
 coeffs(1) = -0.234375000_dp 
 coeffs(2) = 0.093750000_dp
 coeffs(3) = -0.015625000_dp

 !! Coefficients for filter 2
 coeffs2(-5) = 0.5*0.00195312500000_dp
 coeffs2(-4) = -0.5*0.01953125000000_dp
 coeffs2(-3) = 0.5*0.08789062500000_dp
 coeffs2(-2) = -0.5*0.23437500000000_dp
 coeffs2(-1) = 0.5*0.41015625000000_dp
 coeffs2(0) = 0.75390625000000_dp
 coeffs2(1) = 0.5*0.41015625000000_dp
 coeffs2(2) = -0.5*0.23437500000000_dp
 coeffs2(3) = 0.5*0.08789062500000_dp
 coeffs2(4) = -0.5*0.01953125000000_dp
 coeffs2(5) = 0.5*0.00195312500000_dp

 allocate(dfz(-g_z_filt:g_z_filt),dfxy(-g_xy_filt:g_xy_filt))
 if (g_z_filt == 5) dfz = coeffs2
 if (g_z_filt == 3) dfz = coeffs
 if (g_xy_filt == 5) dfxy = coeffs2
 if (g_xy_filt == 3) dfxy = coeffs

 dpz(6,1) = - 0.002484594688_dp
 dpz(6,2) = + 0.020779405824_dp
 dpz(6,3) = - 0.090320001280_dp
 dpz(6,4) = + 0.286511173973_dp
 dpz(6,5) = - 0.872756993962_dp
 dpz(6,6) =   0.0_dp
 dpz(6,7) = + 0.872756993962_dp
 dpz(6,8) = - 0.286511173973_dp
 dpz(6,9) = + 0.090320001280_dp
 dpz(6,10) = - 0.020779405824_dp
 dpz(6,11) = + 0.002484594688_dp

 dpz(1,1) = -2.391602219538_dp
 dpz(1,2) =  5.832490322294_dp
 dpz(1,3) = -7.650218001182_dp + 0.000000000001_dp
 dpz(1,4) =  7.907810563576_dp
 dpz(1,5) = -5.922599052629_dp
 dpz(1,6) =  3.071037015445_dp
 dpz(1,7) = -1.014956769726_dp
 dpz(1,8) =  0.170022256519_dp
 dpz(1,9) =  0.002819958377_dp
 dpz(1,10)= -0.004791009708_dp
 dpz(1,11)= -0.000013063429_dp

 dpz(2,1) = -0.180022054228_dp
 dpz(2,2) = -1.237550583044_dp
 dpz(2,3) =  2.484731692990_dp
 dpz(2,4) = -1.810320814061_dp
 dpz(2,5) =  1.112990048440_dp
 dpz(2,6) = -0.481086916514_dp
 dpz(2,7) =  0.126598690230_dp
 dpz(2,8) = -0.015510730165_dp
 dpz(2,9) =  0.000024609059_dp - 0.000003_dp
 dpz(2,10)=  0.000156447571_dp - 0.000000000001_dp
 dpz(2,11)= -0.000007390277_dp

 dpz(3,1) =  0.057982271137_dp
 dpz(3,2) = -0.536135360383_dp
 dpz(3,3) = -0.264089548967_dp + 0.000000000002_dp
 dpz(3,4) =  0.917445877606_dp - 0.000000000002_dp
 dpz(3,5) = -0.169688364841_dp
 dpz(3,6) = -0.029716326170_dp
 dpz(3,7) =  0.029681617641_dp
 dpz(3,8) = -0.005222483773_dp
 dpz(3,9) = -0.000118806260_dp
 dpz(3,10)= -0.000118806260_dp
 dpz(3,11)= -0.000020069730_dp

 dpz(4,1) = -0.013277273810_dp
 dpz(4,2) =  0.115976072920_dp
 dpz(4,3) = -0.617479187931_dp
 dpz(4,4) = -0.274113948206_dp + 0.000000000002_dp
 dpz(4,5) =  1.086208764655_dp - 0.000000000002_dp
 dpz(4,6) = -0.402951626982_dp
 dpz(4,7) =  0.131066986242_dp
 dpz(4,8) = -0.028154858354_dp
 dpz(4,9) =  0.002596328316_dp
 dpz(4,10)=  0.000128743150_dp
 dpz(4,11)=  0.0_dp

 dpz(5,1) =  0.016756572303_dp
 dpz(5,2) = -0.117478455239_dp
 dpz(5,3) =  0.411034935097_dp
 dpz(5,4) = -1.130286765151_dp
 dpz(5,5) =  0.341435872100_dp - 0.000000000001_dp
 dpz(5,6) =  0.556396830543_dp
 dpz(5,7) = -0.082525734207_dp
 dpz(5,8) =  0.003565834658_dp
 dpz(5,9) =  0.001173034777_dp
 dpz(5,10)= -0.000071772671_dp + 0.000000000064_dp
 dpz(5,11)= -0.000000352273_dp

 dpxy = dpz


end subroutine initexplicitders  

end MODULE INITIALIZE
