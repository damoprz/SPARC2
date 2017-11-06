MODULE DERIVATIVES

 use all_modules 
 implicit none

Contains

!================================================================================================ 

SUBROUTINE ddxyz(var1,dvar1,var2,dvar2,var3,dvar3,bcx,bcy,bcz)

real(dp), INTENT(IN), dimension(dim1(rank),dim2(rank),dim3(rank)) :: var1,var2,var3
real(dp), INTENT(OUT), dimension(dim1(rank),dim2(rank),dim3(rank)) :: dvar1,dvar2,dvar3
integer :: bcx,bcy,bcz

call ddx(var1,dvar1,bcx)
call ddy(var2,dvar2,bcy)
call ddz(var3,dvar3,bcz)

END SUBROUTINE ddxyz

!================================================================================================ 

SUBROUTINE ddx(input,output, bc)

  integer ii,jj,kk,ii2,bc
  real(dp), intent(in),dimension(dim1(rank), dim2(rank),dim3(rank)) :: input
  real(dp), intent(out),dimension(dim1(rank), dim2(rank),dim3(rank)) :: output
  real(dp), dimension(1-g_xy_der:dim1(rank)+g_xy_der,dim2(rank),dim3(rank)) :: temp

  integer d,s,ierr
  integer, dimension(2) :: req1,req2
  integer, dimension(MPI_STATUS_SIZE,2) :: stat1,stat2
  
  temp = 0.0_dp
  output = 0.0_dp

  !! BEGIN DATA TRANSFER

  CALL MPI_CART_SHIFT(comm_cart_3D,0,1,s,d,ierr)
  CALL MPI_ISSEND(input,1,parcel_xrds,d,0,comm_cart_3D,req1(1),ierr)
  CALL MPI_ISSEND(input,1,parcel_xlds,s,1,comm_cart_3D,req1(2),ierr)

  CALL MPI_IRECV(temp,1,parcel_xldr,s,0,comm_cart_3D,req2(1),ierr)
  CALL MPI_IRECV(temp,1,parcel_xrdr,d,1,comm_cart_3D,req2(2),ierr)

  !! DO STUFF WHILE WAITING :) I should calculate all interior derivatives.
  !! insides stay the same

  temp(1:dim1(rank),:,:) = input
  
  do kk = 1,dim3(rank)
    do jj = 1,dim2(rank)
      do ii = g_xy_der+1,dim1(rank)-g_xy_der
        output(ii,jj,kk) = stretchx*SUM(temp(ii-g_xy_der:ii+g_xy_der,jj,kk)*dpxy(g_xy_der+1,:))
      enddo
    enddo
  enddo
    
    call MPI_WAITALL(2,req2,stat2,ierr) !! After the centre is done, wait for the communication


  if (.not. (IAM_XBOT .OR. IAM_XTOP)) then 
  !! IF I AM IN THE MIDDLE

    do kk = 1,dim3(rank)
      do jj = 1,dim2(rank)
        do ii = 1,g_xy_der
          ii2 = dim1(rank)-g_xy_der+ii
          output(ii2,jj,kk) = stretchx*SUM(temp(ii2-g_xy_der:ii2+g_xy_der,jj,kk)*dpxy(g_xy_der+1,:))
          output(ii,jj,kk) = stretchx*SUM(temp(ii-g_xy_der:ii+g_xy_der,jj,kk)*dpxy(g_xy_der+1,:))
        enddo
      enddo
    enddo
   
  else if (IAM_XBOT) then !! I AM THE BOTTOM

    if (bc .EQ. 0) temp(1-g_xy_der:0,:,:) = temp(g_xy_der:1:-1,:,:)  !! Derivative = 0, symetric
    if (bc .EQ. 2) temp(1-g_xy_der:0,:,:) = -temp(g_xy_der:1:-1,:,:)   !! value = 0, antisymetric
    if (bc .EQ. 1) then
      do kk = 1,dim3(rank)
        do jj = 1,dim2(rank)
          do ii = 1,g_xy_der
          ii2 = dim1(rank)-g_xy_der+ii
          output(ii2,jj,kk) = stretchx*SUM(temp(ii2-g_xy_der:ii2+g_xy_der,jj,kk)*dpxy(g_xy_der+1,:))
          output(ii,jj,kk) = stretchx*SUM(temp(1:2*g_xy_der+1,jj,kk)*dpxy(ii,:))
          enddo
        enddo
      enddo
    else
      do kk = 1,dim3(rank)
        do jj = 1,dim2(rank)
          do ii = 1,g_xy_der
          ii2 = dim1(rank)-g_xy_der+ii
          output(ii2,jj,kk) = stretchx*SUM(temp(ii2-g_xy_der:ii2+g_xy_der,jj,kk)*dpxy(g_xy_der+1,:))
          output(ii,jj,kk) = stretchx*SUM(temp(ii-g_xy_der:ii+g_xy_der,jj,kk)*dpxy(g_xy_der+1,:))
          enddo
        enddo
      enddo
    endif

  else !! I AM THE TOP

    if (bc .EQ. 0) temp(dim1(rank)+1:dim1(rank)+1+g_xy_der,:,:) = temp(dim1(rank)-1:dim1(rank)-g_xy_der+1:-1,:,:)!! Derivative = 0, symetric
    if (bc .EQ. 2) temp(dim1(rank)+1:dim1(rank)+1+g_xy_der,:,:) = -temp(dim1(rank)-1:dim1(rank)-g_xy_der+1:-1,:,:)!! Value = 0, antisymetric
    if (bc .EQ. 1) then
      do kk = 1,dim3(rank)
        do jj = 1,dim2(rank)
          do ii = 1,g_xy_der
          ii2 = dim1(rank)-ii+1
          output(ii,jj,kk) = stretchx*SUM(temp(ii-g_xy_der:ii+g_xy_der,jj,kk)*dpxy(g_xy_der+1,:))
          output(ii2,jj,kk) = -stretchx*SUM(temp(dim1(rank):dim1(rank)-2*g_xy_der:-1,jj,kk)*dpxy(ii,:))
          enddo
        enddo
      enddo
    else
      do kk = 1,dim3(rank)
        do jj = 1,dim2(rank)
          do ii = 1,g_xy_der
            ii2 = dim1(rank)-ii+1
            output(ii,jj,kk) = stretchx*SUM(temp(ii-g_xy_der:ii+g_xy_der,jj,kk)*dpxy(g_xy_der+1,:))
            output(ii2,jj,kk) = stretchx*SUM(temp(ii2-g_xy_der:ii2+g_xy_der,jj,kk)*dpxy(g_xy_der+1,:))
          enddo
        enddo
      enddo
    endif
  endif
  call MPI_WAITALL(2,req1,stat1,ierr)

END SUBROUTINE ddx
!================================================================================================ 

SUBROUTINE filtx(input,output)

  integer ii,jj,kk,ii2
  real(dp), intent(in),dimension(dim1(rank), dim2(rank),dim3(rank)) :: input
  real(dp), intent(out),dimension(dim1(rank), dim2(rank),dim3(rank)) :: output
  real(dp), dimension(1-g_xy_filt:dim1(rank)+g_xy_filt,dim2(rank),dim3(rank)) :: temp

  integer d,s,ierr
  integer, dimension(2) :: req1,req2
  integer, dimension(MPI_STATUS_SIZE,2) :: stat1,stat2

  temp = 0.0_dp
  output = 0.0_dp

  !! BEGIN DATA TRANSFER

  CALL MPI_CART_SHIFT(comm_cart_3D,0,1,s,d,ierr)
  CALL MPI_ISSEND(input,1,parcel_xrfs,d,0,comm_cart_3D,req1(1),ierr)
  CALL MPI_ISSEND(input,1,parcel_xlfs,s,1,comm_cart_3D,req1(2),ierr)

  CALL MPI_IRECV(temp,1,parcel_xlfr,s,0,comm_cart_3D,req2(1),ierr)
  CALL MPI_IRECV(temp,1,parcel_xrfr,d,1,comm_cart_3D,req2(2),ierr)

  !! DO STUFF WHILE WAITING :) I should calculate all interior derivatives.
  !! insides stay the same

  temp(1:dim1(rank),:,:) = input

  do kk = 1,dim3(rank)
    do jj = 1,dim2(rank)
      do ii = g_xy_filt+1,dim1(rank)-g_xy_filt
        output(ii,jj,kk) = SUM(temp(ii-g_xy_filt:ii+g_xy_filt,jj,kk)*dfxy(:))
     enddo
   enddo
 enddo
  
 call MPI_WAITALL(2,req2,stat2,ierr) !! After the centre is done, wait for the communication

 if (.not. (IAM_XBOT .OR. IAM_XTOP)) then
   !! IF I AM IN THE MIDDLE

     do kk = 1,dim3(rank)
      do jj = 1,dim2(rank)
        do ii = 1,g_xy_filt
          ii2 = dim1(rank)-g_xy_filt+ii
          output(ii,jj,kk) = SUM(temp(ii-g_xy_filt:ii+g_xy_filt,jj,kk)*dfxy(:))
          output(ii2,jj,kk) = SUM(temp(ii2-g_xy_filt:ii2+g_xy_filt,jj,kk)*dfxy(:))
        enddo
      enddo
    enddo
 
  else if (IAM_XBOT) then !! I AM THE BOTTOM

    if (PERIODIC) then
      do kk = 1,dim3(rank)
        do jj = 1,dim2(rank)
          do ii = 1,g_xy_filt
            ii2 = dim1(rank)+1-ii
            output(ii,jj,kk) = SUM(temp(ii-g_xy_filt:ii+g_xy_filt,jj,kk)*dfxy(:))
            output(ii2,jj,kk) = SUM(temp(ii2-g_xy_filt:ii2+g_xy_filt,jj,kk)*dfxy(:))
          enddo
        enddo
      enddo
    else
      do kk = 1,dim3(rank)
        do jj = 1,dim2(rank)
          do ii = dim1(rank)-g_xy_filt+1,dim1(rank)
            output(ii,jj,kk) = SUM(temp(ii-g_xy_filt:ii+g_xy_filt,jj,kk)*dfxy(:))
          enddo
        enddo
      enddo
      if (g_xy_filt == 3) output(1:g_xy_filt,:,:) = 0.0_dp
      if (g_xy_filt == 5) output(1:g_xy_filt,:,:) = &
      temp(1:g_xy_filt,:,:)
    endif

  else !! I AM THE TOP

    if (PERIODIC) then
      do kk = 1,dim3(rank)
        do jj = 1,dim2(rank)
          do ii = 1,g_xy_filt
            ii2 = dim1(rank)-g_xy_filt+ii
            output(ii,jj,kk) = SUM(temp(ii-g_xy_filt:ii+g_xy_filt,jj,kk)*dfxy(:))
            output(ii2,jj,kk) = SUM(temp(ii2-g_xy_filt:ii2+g_xy_filt,jj,kk)*dfxy(:))
          enddo
        enddo
      enddo
    else
      do kk = 1,dim3(rank)
        do jj = 1,dim2(rank)
          do ii = 1,g_xy_filt
            output(ii,jj,kk) = SUM(temp(ii-g_xy_filt:ii+g_xy_filt,jj,kk)*dfxy(:))
          enddo
        enddo
      enddo
      if (g_xy_filt == 3) output(dim1(rank)-g_xy_filt+1:dim1(rank),:,:) = 0.0_dp
      if (g_xy_filt == 5) output(dim1(rank)-g_xy_filt+1:dim1(rank),:,:) = &
      temp(dim1(rank)-g_xy_filt+1:dim1(rank),:,:)
    endif
  endif
 
  call MPI_WAITALL(2,req1,stat1,ierr)

END SUBROUTINE filtx

!================================================================================================ 
SUBROUTINE ddy(input, output,bc)

  integer ii,jj,kk,ii2,bc
  real*8, intent(in), dimension(dim1(rank), dim2(rank),dim3(rank)) :: input
  real*8, intent(out),dimension(dim1(rank), dim2(rank),dim3(rank)) :: output
  real(dp), dimension(dim1(rank),1-g_xy_der:dim2(rank)+g_xy_der,dim3(rank)) :: temp

  integer d,s,ierr
  integer, dimension(2) :: req1,req2
  integer, dimension(MPI_STATUS_SIZE,2) :: stat1,stat2
 ! Implementing the optimized 11-point centered and non-centred differences of
 ! Bogey, Bailly
 ! 2004 and Berland, Bogey, Marsden and Bailly Marsden 2007, JCP.

  temp = 0.0_dp
  output = 0.0_dp

  !! BEGIN DATA TRANSFER

  CALL MPI_CART_SHIFT(comm_cart_3D,1,1,s,d,ierr)

  CALL MPI_ISSEND(input,1,parcel_yrds,d,1,comm_cart_3D,req1(1),ierr)
  CALL MPI_ISSEND(input,1,parcel_ylds,s,0,comm_cart_3D,req1(2),ierr)

  CALL MPI_IRECV(temp,1,parcel_yldr,s,1,comm_cart_3D,req2(1),ierr)
  CALL MPI_IRECV(temp,1,parcel_yrdr,d,0,comm_cart_3D,req2(2),ierr)

  !! DO STUFF WHILE WAITING :) I should calculate all interior derivatives.
  !! insides stay the same

  temp(:,1:dim2(rank),:) = input

  do kk = 1,dim3(rank)
    do jj = g_xy_der+1,dim2(rank)-g_xy_der
      do ii = 1,dim1(rank)
        output(ii,jj,kk) = stretchy*SUM(temp(ii,jj-g_xy_der:jj+g_xy_der,kk)*dpxy(g_xy_der+1,:))
      enddo
    enddo
  enddo

  call MPI_WAITALL(2,req2,stat2,ierr)

  if (.not. (IAM_YBOT .or. IAM_YTOP)) then
  !! IF I AM IN THE MIDDLE

    do kk = 1,dim3(rank)
      do jj = 1,g_xy_der
        ii2 = dim2(rank)-g_xy_der+jj
        do ii = 1,dim1(rank)

          output(ii,ii2,kk) = stretchy*SUM(temp(ii,ii2-g_xy_der:ii2+g_xy_der,kk)*dpxy(g_xy_der+1,:))
          output(ii,jj,kk) = stretchy*SUM(temp(ii,jj-g_xy_der:jj+g_xy_der,kk)*dpxy(g_xy_der+1,:))

        enddo
      enddo
    enddo


  else if (IAM_YBOT) then !! I AM THE BOTTOM
    
    if (bc .EQ. 0) temp(:,1-g_xy_der:0,:) = temp(:,g_xy_der:1:-1,:)  !! Derivative = 0, symetric
    if (bc .EQ. 2) temp(:,1-g_xy_der:0,:) = -temp(:,g_xy_der:1:-1,:)   !! value = 0, antisymetric
    
    if (bc .EQ. 1) then
      do kk = 1,dim3(rank)
        do jj = 1,g_xy_der
          ii2 = dim2(rank)-g_xy_der+jj
          do ii = 1,dim1(rank)
            output(ii,ii2,kk) = stretchy*SUM(temp(ii,ii2-g_xy_der:ii2+g_xy_der,kk)*dpxy(g_xy_der+1,:))
            output(ii,jj,kk) = stretchy*SUM(temp(ii,1:2*g_xy_der+1,kk)*dpxy(jj,:))
          enddo
        enddo
      enddo
    else
      do kk = 1,dim3(rank)
        do jj = 1,g_xy_der
          ii2 = dim2(rank)-g_xy_der+jj
          do ii = 1,dim1(rank)
            output(ii,ii2,kk) = stretchy*SUM(temp(ii,ii2-g_xy_der:ii2+g_xy_der,kk)*dpxy(g_xy_der+1,:))
            output(ii,jj,kk) = stretchy*SUM(temp(ii,jj-g_xy_der:jj+g_xy_der,kk)*dpxy(g_xy_der+1,:))
          enddo
        enddo
      enddo
    endif

  else if (IAM_YTOP) then !! I AM THE TOP

    if (bc .eq. 0) temp(:,dim2(rank)+1:dim2(rank)+1+g_xy_der,:) = temp(:,dim2(rank)-1:dim2(rank)-g_xy_der+1:-1,:)!! Derivative = 0, symetric
    if (bc .eq. 2) temp(:,dim2(rank)+1:dim2(rank)+1+g_xy_der,:) = -temp(:,dim2(rank)-1:dim2(rank)-g_xy_der+1:-1,:)!! Value = 0, antisymetric

    if (bc .EQ. 1) then
      do kk = 1,dim3(rank)
        do jj = 1,g_xy_der
          ii2 = dim2(rank)+1-jj
          do ii = 1,dim1(rank)
            output(ii,jj,kk) = stretchy*SUM(temp(ii,jj-g_xy_der:jj+g_xy_der,kk)*dpxy(g_xy_der+1,:))
            output(ii,ii2,kk) = -stretchy*SUM(temp(ii,dim2(rank):dim2(rank)-2*g_xy_der:-1,kk)*dpxy(jj,:))
          enddo
        enddo
      enddo
    else
      do kk = 1,dim3(rank)
        do jj = 1,g_xy_der
          ii2 = dim2(rank)-g_xy_der+jj 
          do ii = 1,dim1(rank)
            output(ii,ii2,kk) = stretchy*SUM(temp(ii,ii2-g_xy_der:ii2+g_xy_der,kk)*dpxy(g_xy_der+1,:))
            output(ii,jj,kk) = stretchy*SUM(temp(ii,jj-g_xy_der:jj+g_xy_der,kk)*dpxy(g_xy_der+1,:))
          enddo
        enddo
      enddo
    endif
  endif

  call MPI_WAITALL(2,req1,stat1,ierr)

END SUBROUTINE ddy

!================================================================================================ 

SUBROUTINE filty(input,output)

  integer ii,jj,kk,ii2
  real(dp), intent(in), dimension(dim1(rank), dim2(rank),dim3(rank)) :: input
  real(dp), intent(out), dimension(dim1(rank), dim2(rank),dim3(rank)) :: output
  real(dp), dimension(dim1(rank),1-g_xy_filt:dim2(rank)+g_xy_filt,dim3(rank)) :: temp

  integer d, s,ierr
  integer, dimension(2) :: req1,req2
  integer, dimension(MPI_STATUS_SIZE,2) :: stat1,stat2

  temp = 0.0_dp
  output = 0.0_dp

  !! BEGIN DATA TRANSFER

  CALL MPI_CART_SHIFT(comm_cart_3D,1,1,s,d,ierr)
  CALL MPI_ISSEND(input,1,parcel_yrfs,d,0,comm_cart_3D,req1(1),ierr)
  CALL MPI_ISSEND(input,1,parcel_ylfs,s,1,comm_cart_3D,req1(2),ierr)

  CALL MPI_IRECV(temp,1,parcel_ylfr,s,0,comm_cart_3D,req2(1),ierr)
  CALL MPI_IRECV(temp,1,parcel_yrfr,d,1,comm_cart_3D,req2(2),ierr)

  !! DO STUFF WHILE WAITING :) I should calculate all interior derivatives.
  !! insides stay the same

  temp(:,1:dim2(rank),:) = input

  do kk = 1,dim3(rank)
    do ii = 1,dim1(rank)
      do jj = g_xy_filt+1,dim2(rank)-g_xy_filt
        output(ii,jj,kk) = sum(temp(ii,jj-g_xy_filt:jj+g_xy_filt,kk)*dfxy(:))
      enddo
    enddo
  enddo

  call MPI_WAITALL(2,req2,stat2,ierr)

  if (.not. (IAM_YBOT .OR. IAM_YTOP)) then
    do kk = 1,dim3(rank)
        do ii = 1,dim1(rank)
            do jj = 1,g_xy_filt
                ii2=dim2(rank)-g_xy_filt+jj
                output(ii,jj,kk) = sum(temp(ii,jj-g_xy_filt:jj+g_xy_filt,kk)*dfxy(:))
                output(ii,ii2,kk) = sum(temp(ii,ii2-g_xy_filt:ii2+g_xy_filt,kk)*dfxy(:))
            enddo
        enddo
    enddo
  
  else if (IAM_YBOT) then !! I AM THE BOTTOM

    if (periodic) then
      do kk = 1,dim3(rank)
        do ii = 1,dim1(rank)
          do jj = 1,g_xy_filt !do top and bottom
            ii2 = dim2(rank)-g_xy_filt+jj
                output(ii,jj,kk) = sum(temp(ii,jj-g_xy_filt:jj+g_xy_filt,kk)*dfxy(:))
                output(ii,ii2,kk) = sum(temp(ii,ii2-g_xy_filt:ii2+g_xy_filt,kk)*dfxy(:))
          enddo          
        enddo
      enddo
    else
      do kk = 1,dim3(rank)
        do ii = 1,dim1(rank)
           do jj = dim2(rank)-g_xy_filt+1,dim2(rank) !do top
             output(ii,jj,kk) = sum(temp(ii,jj-g_xy_filt:jj+g_xy_filt,kk)*dfxy(:))
           enddo
         enddo
       enddo
      if (g_xy_filt == 3) output(:,1:g_xy_filt,:) = 0.0_dp
      if (g_xy_filt == 5) output(:,1:g_xy_filt,:) = &
      temp(:,1:g_xy_filt,:)
    endif

  else !! I AM THE TOP

    if (periodic) then
      do kk = 1,dim3(rank)
        do ii = 1,dim1(rank)
          do jj = 1,g_xy_filt !! Do Bottom + Top
            ii2=dim2(rank)-g_xy_filt+jj
                output(ii,jj,kk) = sum(temp(ii,jj-g_xy_filt:jj+g_xy_filt,kk)*dfxy(:))
                output(ii,ii2,kk) = sum(temp(ii,ii2-g_xy_filt:ii2+g_xy_filt,kk)*dfxy(:))
          enddo
        enddo
      enddo
    else
      do kk = 1,dim3(rank)
        do ii = 1,dim1(rank)
          do jj = 1,g_xy_filt !! Do Bottom
            output(ii,jj,kk) = sum(temp(ii,jj-g_xy_filt:jj+g_xy_filt,kk)*dfxy(:))
          enddo
        enddo
      enddo
      if (g_xy_filt == 3) output(:,dim2(rank)-g_xy_filt+1:dim2(rank),:) = 0.0_dp
      if (g_xy_filt == 5) output(:,dim2(rank)-g_xy_filt+1:dim2(rank),:) = &
      temp(:,dim2(rank)-g_xy_filt+1:dim2(rank),:)
    endif
  endif

  call MPI_WAITALL(2,req1,stat1,ierr)

END SUBROUTINE filty

!================================================================================================ 

SUBROUTINE ddz(input, output,bc)

  integer ii,jj,kk,ii2,bc
  real(dp), intent(in), dimension(dim1(rank), dim2(rank),dim3(rank)) :: input
  real(dp), intent(out), dimension(dim1(rank), dim2(rank),dim3(rank)) :: output
  real(dp), dimension(dim1(rank),dim2(rank),1-g_z_der:dim3(rank)+g_z_der) :: temp

  integer d, s,ierr
  integer, dimension(2) :: req1,req2
  integer, dimension(MPI_STATUS_SIZE,2) :: stat1,stat2
 ! Implementing the optimized 11-point centered and non-centred differences of
 ! Bogey, Bailly
 ! 2004 and Berland, Bogey, Marsden and Bailly 2007, JCP.

  temp = 0.0_dp
  output = 0.0_dp

  !! BEGIN DATA TRANSFER

  CALL MPI_CART_SHIFT(comm_cart_3D,2,1,s,d,ierr)
  CALL MPI_ISSEND(input,1,parcel_zrds,d,0,comm_cart_3D,req1(1),ierr)
  CALL MPI_ISSEND(input,1,parcel_zlds,s,1,comm_cart_3D,req1(2),ierr)

  CALL MPI_IRECV(temp,1,parcel_zldr,s,0,comm_cart_3D,req2(1),ierr)
  CALL MPI_IRECV(temp,1,parcel_zrdr,d,1,comm_cart_3D,req2(2),ierr)

  !! DO STUFF WHILE WAITING :) I should calculate all interior derivatives.
  !! insides stay the same
  temp(:,:,1:dim3(rank)) = input

  do kk = g_z_der+1,dim3(rank)-g_z_der
    do jj = 1,dim2(rank)
      do ii = 1,dim1(rank)
        output(ii,jj,kk) = stretch(kk)*SUM(temp(ii,jj,kk-g_z_der:kk+g_z_der)*dpz(g_z_der+1,:))
      enddo
    enddo
  enddo

  call MPI_WAITALL(2,req2,stat2,ierr)

    
  if (.NOT. (IAM_ZTOP .or. IAM_ZBOT)) then !! I AM IN THE MIDDLE

    do kk = 1,g_z_der
      ii2 = dim3(rank)-g_z_der+kk
      do jj = 1,dim2(rank)
        do ii = 1,dim1(rank)
          output(ii,jj,ii2) = stretch(ii2)*SUM(temp(ii,jj,ii2-g_z_der:ii2+g_z_der)*dpz(g_z_der+1,:))
          output(ii,jj,kk) = stretch(kk)*SUM(temp(ii,jj,kk-g_z_der:kk+g_z_der)*dpz(g_z_der+1,:))
        enddo
      enddo
    enddo

  else if (IAM_ZBOT) then !!!! Now do the edges, depending on location

    if (bc .EQ. 0) temp(:,:,1-g_z_der:0) = temp(:,:,g_z_der:1:-1)  !! Derivative = 0, symetric
    if (bc .EQ. 2) temp(:,:,1-g_z_der:0) = -temp(:,:,g_z_der:1:-1)   !! value = 0, antisymetric

    if (bc .EQ. 1) then
      do kk = 1,g_z_der
        ii2 = dim3(rank)-g_z_der+kk
        do jj = 1,dim2(rank)
          do ii = 1,dim1(rank)
            output(ii,jj,ii2) = stretch(ii2)*SUM(temp(ii,jj,ii2-g_z_der:ii2+g_z_der)*dpz(g_z_der+1,:))
            output(ii,jj,kk) = stretch(kk)*SUM(temp(ii,jj,kk-g_z_der:kk+g_z_der)*dpz(g_z_der+1,:))
          enddo
        enddo
      enddo
    else
      do kk = 1,g_z_der
        ii2 = dim3(rank)-g_z_der+kk
        do jj = 1,dim2(rank)
          do ii = 1,dim1(rank)
              output(ii,jj,ii2) = stretch(ii2)*SUM(temp(ii,jj,ii2-g_z_der:ii2+g_z_der)*dpz(g_z_der+1,:))
              output(ii,jj,kk) = stretch(kk)*SUM(temp(ii,jj,1:2*g_z_der+1)*dpz(kk,:))
          enddo
        enddo
      enddo
    endif

  else !! I AM THE TOP
    
    if (bc .eq. 0) temp(:,:,dim3(rank)+1:dim3(rank)+1+g_z_der) = temp(:,:,dim3(rank)-1:dim3(rank)-g_z_der+1:-1)!! Derivative = 0, symetric
    if (bc .eq. 2) temp(:,:,dim3(rank)+1:dim3(rank)+1+g_z_der) = -temp(:,:,dim3(rank)-1:dim3(rank)-g_z_der+1:-1)!! Value = 0, antisymetric

    if (BC .EQ. 1) then
      do kk = 1,g_z_der
        ii2 = dim3(rank)+1-kk
        do jj = 1,dim2(rank)
          do ii = 1,dim1(rank)
            output(ii,jj,kk) = stretch(kk)*SUM(temp(ii,jj,kk-g_z_der:kk+g_z_der)*dpz(g_z_der+1,:))
            output(ii,jj,ii2) = -stretch(ii2)*SUM(temp(ii,jj,dim3(rank):dim3(rank)-2*g_z_der:-1)*dpz(kk,:))
          enddo
        enddo
      enddo
    else
      do kk = 1,g_z_der
        ii2 = dim3(rank)+1-kk
        do jj = 1,dim2(rank)
          do ii = 1,dim1(rank)
            output(ii,jj,ii2) = stretch(ii2)*SUM(temp(ii,jj,ii2-g_z_der:ii2+g_z_der)*dpz(g_z_der+1,:))
            output(ii,jj,kk) = stretch(kk)*SUM(temp(ii,jj,kk-g_z_der:kk+g_z_der)*dpz(g_z_der+1,:))
          enddo
        enddo
      enddo
    endif

  endif
 
  call MPI_WAITALL(2,req1,stat1,ierr)

  END SUBROUTINE ddz

!================================================================================================
SUBROUTINE filtz(input,output)

  integer ii,jj,kk,ii2
  real(dp), INTENT(IN), dimension(dim1(rank), dim2(rank),dim3(rank)) :: input
  real(dp), INTENT(OUT), dimension(dim1(rank), dim2(rank),dim3(rank)) :: output
  real(dp), dimension(dim1(rank),dim2(rank),1-g_z_filt:dim3(rank)+g_z_filt) :: temp

  integer d, s,ierr
  integer, dimension(2) :: req1,req2
  integer, dimension(MPI_STATUS_SIZE,2) :: stat1,stat2

  temp = 0.0_dp
  output = 0.0_dp

  !! BEGIN DATA TRANSFER


  CALL MPI_CART_SHIFT(comm_cart_3D,2,1,s,d,ierr)
  CALL MPI_ISSEND(input,1,parcel_zrfs,d,0,comm_cart_3D,req1(1),ierr)
  CALL MPI_ISSEND(input,1,parcel_zlfs,s,1,comm_cart_3D,req1(2),ierr)

  CALL MPI_IRECV(temp,1,parcel_zlfr,s,0,comm_cart_3D,req2(1),ierr)
  CALL MPI_IRECV(temp,1,parcel_zrfr,d,1,comm_cart_3D,req2(2),ierr)

  !! DO STUFF WHILE WAITING :) I should calculate all interior derivatives.
  !! insides stay the same
  temp(:,:,1:dim3(rank)) = input

  do kk = g_z_filt+1,dim3(rank)-g_z_filt
    do jj = 1,dim2(rank)
      do ii = 1,dim1(rank)
        output(ii,jj,kk) = sum(temp(ii,jj,kk-g_z_filt:kk+g_z_filt)*dfz(:))
      enddo
    enddo
  enddo

  call MPI_WAITALL(2,req2,stat2,ierr)
  
  if (.NOT. (IAM_ZTOP .or. IAM_ZBOT)) then
   !! IF I AM IN THE MIDDLE


    do kk = 1,g_z_filt
    ii2 = dim3(rank)-g_z_filt+kk
      do jj = 1,dim2(rank)
        do ii = 1,dim1(rank)
          output(ii,jj,kk) = sum(temp(ii,jj,kk-g_z_filt:kk+g_z_filt)*dfz(:)) 
          output(ii,jj,ii2) = sum(temp(ii,jj,ii2-g_z_filt:ii2+g_z_filt)*dfz(:)) 
        enddo
      enddo
    enddo

   else if (IAM_ZBOT) then !! I AM THE BOTTOM

   do jj = 1,dim2(rank)
        do ii = 1,dim1(rank)
          do kk = dim3(rank)-g_z_filt+1,dim3(rank)
            output(ii,jj,kk) = sum(temp(ii,jj,kk-g_z_filt:kk+g_z_filt)*dfz(:))
          enddo
          output(ii,jj,1:g_z_filt) = temp(ii,jj,1:g_z_filt)
       enddo
    enddo

  else !! I AM THE TOP

   do jj = 1,dim2(rank)
        do ii = 1,dim1(rank)
            do kk = 1,g_z_filt
                output(ii,jj,kk) = sum(temp(ii,jj,kk-g_z_filt:kk+g_z_filt)*dfz(:))
            enddo
            output(ii,jj,dim3(rank)-g_z_filt+1:dim3(rank)) = temp(ii,jj,dim3(rank)-g_z_filt+1:dim3(rank))
      enddo
    enddo
  endif
  
  call MPI_WAITALL(2,req1,stat1,ierr)

  END SUBROUTINE filtz

!================================================================================================
                                                                                                                                                            
SUBROUTINE CROSS(aa_x, aa_y, aa_z, bb_x, bb_y, bb_z, cc_x, cc_y, cc_z)
  
  real(dp), intent(in), dimension(dim1(rank), dim2(rank), dim3(rank)) :: aa_x, aa_y, aa_z, bb_x, bb_y, bb_z
  real(dp), intent(out), dimension(dim1(rank), dim2(rank), dim3(rank)) :: cc_x, cc_y, cc_z

  cc_x = aa_y * bb_z - aa_z * bb_y
  cc_y = bb_x * aa_z - aa_x * bb_z
  cc_z = aa_x * bb_y - aa_y * bb_x

END SUBROUTINE CROSS

!================================================================================================

SUBROUTINE FILTER_VARS()
  
  integer ll
  real(dp), dimension(dim1(rank),dim2(rank),dim3(rank)) :: filt
  integer ierr

  !! For this routine: do l = 1,dimfive: if xy-filt timestep then filter (x,y) rotate order, If z-filt timestep then filter z.

  do ll=1,dimfive
    
    if (mod(time,xyfreq_filtering) .EQ. 0) then
      if (.not. TWOD) then
        if (.not. PERIODIC) then
          if (IAM_YTOP) a(:,dim2(rank)-g_xy_filt+1:dim2(rank)-1,:,ll) = &
          (a(:,dim2(rank)-g_xy_filt+2:dim2(rank),:,ll) + a(:,dim2(rank)-g_xy_filt:dim2(rank)-2,:,ll))/2.0_dp
          if (IAM_YBOT) a(:,2:g_xy_filt,:,ll) = (a(:,1:g_xy_filt-1,:,ll)+a(:,3:g_xy_filt+1,:,ll))/2.0_dp
        endif
        call filty(a(:,:,:,ll),filt)
        if (g_xy_filt == 3) a(:,:,:,ll) = a(:,:,:,ll)-s_xy*filt
        if (g_xy_filt == 5) a(:,:,:,ll) = filt
      endif
      if (.not. PERIODIC) then
        if (IAM_XTOP) a(dim1(rank)-g_xy_filt+1:dim1(rank)-1,:,:,ll) = &
        (a(dim1(rank)-g_xy_filt+2:dim1(rank),:,:,ll) + a(dim1(rank)-g_xy_filt:dim1(rank)-2,:,:,ll))/2.0_dp
        if (IAM_XBOT) a(2:g_xy_filt,:,:,ll) = (a(1:g_xy_filt-1,:,:,ll)+a(3:g_xy_filt+1,:,:,ll))/2.0_dp
      endif
      call filtx(a(:,:,:,ll),filt)
      if (g_xy_filt == 3) a(:,:,:,ll) = a(:,:,:,ll)-s_xy*filt
      if (g_xy_filt == 5) a(:,:,:,ll) = filt
    endif

    if (mod(time,zfreq_filtering) .EQ. 0) then
      if (IAM_ZTOP) a(:,:,dim3(rank)-g_z_filt+1:dim3(rank)-1,ll) = &
      (a(:,:,dim3(rank)-g_z_filt+2:dim3(rank),ll) + a(:,:,dim3(rank)-g_z_filt:dim3(rank)-2,ll))/2.0_dp
      if (IAM_ZBOT) a(:,:,2:g_z_filt,ll) = (a(:,:,1:g_z_filt-1,ll)+a(:,:,3:g_z_filt+1,ll))/2.0_dp
      call filtz(a(:,:,:,ll),filt)
      if (g_z_filt == 3) a(:,:,:,ll) = a(:,:,:,ll)-s_z*filt
      if (g_z_filt == 5) a(:,:,:,ll) = filt
    endif

  enddo

  call MPI_BARRIER(MPI_COMM_WORLD,ierr)

END SUBROUTINE FILTER_VARS

!================================================================================================ 

END MODULE DERIVATIVES
