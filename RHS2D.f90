Module PHYSICS_RHS_2D

! --------------------------------------------------------------------------
! MPI Version of the Cartesian Solar Wave Simulator.
! Copyright 2008, Shravan Hanasoge
                                                                                                                                                         
! W. W. Hansen Experimental Physics Laboratory
! Stanford University, Stanford, CA 94305, USA
! Email: shravan@solar.stanford.edu
! --------------------------------------------------------------------------

  use all_modules
  use physics_RHS
  use derivatives
  
  implicit none

Contains
!============================================================================

SUBROUTINE MP_MHD_2D


 ! NON-CONSERVATIVE FORM OF THE FLUXES.

 integer kk
 real(dp), dimension(dim1(rank),1,dim3(rank)) :: flux1, flux2, flux3
 real(dp), dimension(dim1(rank),1,dim3(rank)) :: dxflux2,dxflux3,dzflux2,dzflux1
 real(dp), dimension(dim1(rank),1,dim3(rank)) :: dzbx,dzby,dxby,dxbz
 real(dp), dimension(dim1(rank),dim2(rank),dim3(rank)) :: term,gradrho_x,gradrho_z,&
        div0,dvxdz,dvydx,dvydz,dvzdx
 
 integer bh_vel,bv_vel,bh_mag,bv_mag,bh_prerho,bv_prerho

 if (periodic) then
 bh_vel = 5
 bv_vel = 1
 bh_mag = 5
 bv_mag = 1
 bh_prerho = 5
 bv_prerho = 0
 else
 bh_vel = 1
 bv_vel = 1
 bh_mag = 1
 bv_mag = 1
 bh_prerho = 0
 bv_prerho = 0
 endif

 if (TWOP5D == 1) flux1 = - v_z * boy + v_y * boz
 flux2 = - v_x * boz + v_z * box
 if (TWOP5D == 1) flux3 = - v_y * box + v_x * boy

 !!! Call the flux derivatives

 call ddz(flux2,dzflux2,bv_mag)
 call ddz(flux1,dzflux1,bv_mag)
 call ddx(flux2,dxflux2,bh_mag)
 call ddx(flux3,dxflux3,bh_mag)

 call ddx(v_x,dvxdx,bh_vel)
 call ddz(v_z,dvzdz,bv_vel)

 call ddz(by,dzby,bv_mag)
 call ddz(bx,dzbx,bv_mag)
 call ddx(bz,dxbz,bh_mag)
 call ddx(by,dxby,bh_mag)

 call ddx(p,gradp_x,bh_prerho)
 call ddz(p,gradp_z,bv_prerho)

  if (USE_PML) then
    if (IAM_ZBOT) then
      RHSpsizp = az*gradp_z(:,:,1:nzpmlbot) + bzpml*psizp
      gradp_z(:,:,1:nzpmlbot) =  gradp_z(:,:,1:nzpmlbot)/kappaz + psizp
      RHSpsizvz = az*dvzdz(:,:,1:nzpmlbot) + bzpml*psizvz
      dvzdz(:,:,1:nzpmlbot) = dvzdz(:,:,1:nzpmlbot)/kappaz + psizvz
    endif
    if (IAM_ZTOP) then
      RHSpsizp = az*gradp_z(:,:,dim3(rank)-nzpmltop+1:dim3(rank)) + bzpml*psizp
      gradp_z(:,:,dim3(rank)-nzpmltop+1:dim3(rank)) = &
gradp_z(:,:,dim3(rank)-nzpmltop+1:dim3(rank))/kappaz + psizp
      RHSpsizvz = az*dvzdz(:,:,dim3(rank)-nzpmltop+1:dim3(rank)) + bzpml*psizvz
      dvzdz(:,:,dim3(rank)-nzpmltop+1:dim3(rank)) = &
dvzdz(:,:,dim3(rank)-nzpmltop+1:dim3(rank))/kappaz + psizvz
    endif
  endif


  if (USE_HPML) then
    if (IAM_XBOT) then
      RHSpsixp = ax*gradp_x(1:nxpmlbot,:,:) + bxpml*psixp
      gradp_x(1:nxpmlbot,:,:) =  gradp_x(1:nxpmlbot,:,:)/kappax + psixp
      RHSpsixvx = ax*dvxdx(1:nxpmlbot,:,:) + bxpml*psixvx
      dvxdx(1:nxpmlbot,:,:) = dvxdx(1:nxpmlbot,:,:)/kappax + psixvx
    endif
    if (IAM_XTOP) then
      RHSpsixp = ax*gradp_x(dim1(rank)-nxpmltop+1:dim1(rank),:,:) + bxpml*psixp
      gradp_x(dim1(rank)-nxpmltop+1:dim1(rank),:,:) = &
gradp_x(dim1(rank)-nxpmltop+1:dim1(rank),:,:)/kappax + psixp
      RHSpsixvx = ax*dvxdx(dim1(rank)-nxpmltop+1:dim1(rank),:,:) + bxpml*psixvx
      dvxdx(dim1(rank)-nxpmltop+1:dim1(rank),:,:) = &
dvxdx(dim1(rank)-nxpmltop+1:dim1(rank),:,:)/kappax + psixvx
    endif
  endif

    if (USE_PML) then
    if (IAM_ZBOT) then
    RHSpsizinductionbx = az*dzflux2(:,:,1:nzpmlbot) + bzpml*psizinductionbx 
    dzflux2(:,:,1:nzpmlbot) = dzflux2(:,:,1:nzpmlbot)/kappaz + psizinductionbx
    if (TWOP5D == 1) RHSpsizinductionby = az*dzflux1(:,:,1:nzpmlbot) + bzpml*psizinductionby
    if (TWOP5D == 1) dzflux1(:,:,1:nzpmlbot) = dzflux1(:,:,1:nzpmlbot)/kappaz + psizinductionby
    RHSpsizdzbx = az*dzbx(:,:,1:nzpmlbot) + bzpml*psizdzbx
    dzbx(:,:,1:nzpmlbot) = dzbx(:,:,1:nzpmlbot)/kappaz + psizdzbx
    if (TWOP5D == 1) RHSpsizdzby = az*dzby(:,:,1:nzpmlbot) + bzpml*psizdzby
    if (TWOP5D == 1) dzby(:,:,1:nzpmlbot) = dzby(:,:,1:nzpmlbot)/kappaz + psizdzby
    endif
    if (IAM_ZTOP) then
    RHSpsizinductionbx = az*dzflux2(:,:,dim3(rank)-nzpmltop+1:dim3(rank)) + bzpml*psizinductionbx
    dzflux2(:,:,dim3(rank)-nzpmltop+1:dim3(rank)) = &
dzflux2(:,:,dim3(rank)-nzpmltop+1:dim3(rank))/kappaz + psizinductionbx
    if (TWOP5D == 1) RHSpsizinductionby = az*dzflux1(:,:,dim3(rank)-nzpmltop+1:dim3(rank)) + bzpml*psizinductionby
    if (TWOP5D == 1) dzflux1(:,:,dim3(rank)-nzpmltop+1:dim3(rank)) = &
    dzflux1(:,:,dim3(rank)-nzpmltop+1:dim3(rank))/kappaz + psizinductionby
    RHSpsizdzbx = az*dzbx(:,:,dim3(rank)-nzpmltop+1:dim3(rank)) + bzpml*psizdzbx
    dzbx(:,:,dim3(rank)-nzpmltop+1:dim3(rank)) = &
    dzbx(:,:,dim3(rank)-nzpmltop+1:dim3(rank))/kappaz + psizdzbx
    if (TWOP5D == 1) RHSpsizdzby = az*dzby(:,:,dim3(rank)-nzpmltop+1:dim3(rank)) + bzpml*psizdzby
    if (TWOP5D == 1) dzby(:,:,dim3(rank)-nzpmltop+1:dim3(rank)) = &
    dzby(:,:,dim3(rank)-nzpmltop+1:dim3(rank))/kappaz + psizdzby
    endif
  endif
  if (USE_HPML) then
    if (IAM_XBOT) then
      if (TWOP5D == 1) RHSpsixinductionby = ax*dxflux3(1:nxpmlbot,:,:) + bxpml*psixinductionby
      if (TWOP5D == 1) dxflux3(1:nxpmlbot,:,:) = dxflux3(1:nxpmlbot,:,:)/kappax + psixinductionby
      RHSpsixinductionbz = ax*dxflux2(1:nxpmlbot,:,:) + bxpml*psixinductionbz
      dxflux2(1:nxpmlbot,:,:) = dxflux2(1:nxpmlbot,:,:)/kappax + psixinductionbz
      if (TWOP5D == 1) RHSpsixdxby = ax*dxby(1:nxpmlbot,:,:) + bxpml*psixdxby
      if (TWOP5D == 1) dxby(1:nxpmlbot,:,:) = dxby(1:nxpmlbot,:,:)/kappax + psixdxby
      RHSpsixdxbz = ax*dxbz(1:nxpmlbot,:,:) + bxpml*psixdxbz
      dxbz(1:nxpmlbot,:,:) = dxbz(1:nxpmlbot,:,:)/kappax + psixdxbz
    endif
    if (IAM_XTOP) then
      if (TWOP5D == 1) RHSpsixinductionby = ax*dxflux3(dim1(rank)-nxpmltop+1:dim1(rank),:,:) + &
      bxpml*psixinductionby
      if (TWOP5D == 1) dxflux3(dim1(rank)-nxpmltop+1:dim1(rank),:,:) = &
      dxflux3(dim1(rank)-nxpmltop+1:dim1(rank),:,:)/kappax + psixinductionby
      RHSpsixinductionbz = ax*dxflux2(dim1(rank)-nxpmltop+1:dim1(rank),:,:) + &
      bxpml*psixinductionbz
      dxflux2(dim1(rank)-nxpmltop+1:dim1(rank),:,:) = &
      dxflux2(dim1(rank)-nxpmltop+1:dim1(rank),:,:)/kappax + psixinductionbz
      if (TWOP5D == 1) RHSpsixdxby = ax*dxby(dim1(rank)-nxpmltop+1:dim1(rank),:,:) + bxpml*psixdxby
      if (TWOP5D == 1) dxby(dim1(rank)-nxpmltop+1:dim1(rank),:,:) = &
      dxby(dim1(rank)-nxpmltop+1:dim1(rank),:,:)/kappax + psixdxby
      RHSpsixdxbz = ax*dxbz(dim1(rank)-nxpmltop+1:dim1(rank),:,:) + bxpml*psixdxbz
      dxbz(dim1(rank)-nxpmltop+1:dim1(rank),:,:) = &
      dxbz(dim1(rank)-nxpmltop+1:dim1(rank),:,:)/kappax + psixdxbz
    endif
  endif

 RHSb_x = -reduction*dzflux2
 if (TWOP5D == 1) RHSb_y = reduction*(dzflux1 - dxflux3)
 RHSb_z = reduction*dxflux2 

 div = dvxdx + dvzdz

 curlbx =  -dzby*TWOP5D
 if (TWOP5D == 1) curlby = -dxbz + dzbx
 curlbz = dxby*TWOP5D


 ! CONTIUNUITY 
 ! --------------------------------------

  RHScont = - gradrho0_x * v_x &
           - gradrho0_z * v_z - rho0 * div - spongexyz*rho

 ! PRESSURE


 RHSp = - c2rho0 * div - v_x * gradp0_x &
        - v_z * gradp0_z - spongexyz*p

 call cross(curlbx, curlby, curlbz, box, boy, boz, flux1, flux2, flux3)

 RHSv_x = -rhoinv*(gradp_x - flux1*reduction)-spongexyz*v_x
 if (TWOP5D == 1) RHSv_y = -rhoinv*( - flux2*reduction)-spongexyz*v_y
 RHSv_z = -rhoinv*(gradp_z - flux3*reduction)-spongexyz*v_z

 call cross(curlbox, curlboy, curlboz, bx, by, bz, flux1, flux2, flux3)

 RHSv_x = RHSv_x + rhoinv*(flux1*reduction)
 if (TWOP5D == 1) RHSv_y = RHSv_y + TWOP5D*rhoinv*(flux2*reduction)
 RHSv_z = RHSv_z + rhoinv*(flux3*reduction)

 if (FLOWS) then

 if (TWOP5D == 1) flux1 = - v0_z * by + v0_y * bz
 flux2 = - v0_x * bz + v0_z * bx
 if (TWOP5D == 1) flux3 = - v0_y * bx + v0_x * by
 
 call ddz(flux2,dzflux2,bv_mag)
 if (TWOP5D == 1) call ddz(flux1,dzflux1,bv_mag)
 call ddx(flux2,dxflux2,bh_mag)
 if (TWOP5D == 1) call ddx(flux3,dxflux3,bh_mag)

 RHSb_x = RHSb_x -reduction*dzflux2
 if (TWOP5D == 1) RHSb_y = RHSb_y + TWOP5D*reduction*(dzflux1 - dxflux3)
 RHSb_z = RHSb_z + reduction*dxflux2
 
 RHSp = RHSp - gradp_x * v0_x &
           - gradp_z * v0_z

 call ddz(rho,gradrho_z,bv_prerho)
 call ddx(rho,gradrho_x,bh_prerho)

 RHScont = RHScont - gradrho_x * v0_x &
           - gradrho_z * v0_z - rho * div0

  call ddz(v_x,dvxdz,bv_vel)
  call ddx(v_z, dvzdx,bh_vel)

  RHSv_x = RHSv_x - 2.0_dp *(v0_x*dvxdx - v0_z*dvxdz - &
v_x*dv0xdx - v_z*dv0xdz)
  RHSv_z = RHSv_z - 2.0_dp *(v0_z*dvzdz - v0_x*dvzdx - &
v_z*dv0zdz - v_x*dv0zdx)
  if (TWOP5D == 1) then

   call ddz(v_y,dvydz,bv_vel)
   call ddx(v_y, dvydx,bh_vel)
   RHSv_y = RHSv_y - TWOP5D*2.0_dp *(v0_z*dvydz - v0_x*dvydx - &
v_z*dv0ydz - v_x*dv0ydx)
  endif
 endif

 do kk=1,dim3(rank)
  RHSv_z(:,:,kk) = RHSv_z(:,:,kk)  - rhoinv(:,:,kk) * rho(:,:,kk) * g(kk) 

  if (pulse) then
    if (pulse_dir == 1) then
      RHSv_x(:,:,kk) = RHSv_x(:,:,kk) + forcing_p(:,:,kk)
    else if (pulse_dir == 2) then
      RHSv_y(:,:,kk) = RHSv_y(:,:,kk) + TWOP5D*forcing_p(:,:,kk)
    else if (pulse_dir == 3) then
      RHSv_z(:,:,kk) = RHSv_z(:,:,kk) + forcing_p(:,:,kk)
    endif
  else
    RHSv_z(:,:,kk) = RHSv_z(:,:,kk) +  forcing*source_dep(kk)
  endif
 enddo

 if (IAM_NOTB) then
 else
   if (.not. PERIODIC) then
     if (IAM_XTOP) scr(dim1(rank),:,:,:) = 0.0_dp
     if (IAM_XBOT) scr(1,:,:,:) = 0.0_dp
   endif
   if (IAM_ZTOP) scr(:,:,dim3(rank),:) = 0.0_dp
   if (IAM_ZBOT) scr(:,:,1,:) = 0.0_dp
 endif

 if (DAMPING) then
  call damp_velocity(v_x,term)
  RHSv_x = RHSv_x - term
  if (TWOP5D == 1) call damp_velocity(v_y,term)
  if (TWOP5D == 1) RHSv_y = RHSv_y - TWOP5D*term
  call damp_velocity(v_z,term)
  RHSv_z = RHSv_z - term
 endif

 END SUBROUTINE MP_MHD_2D

!================================================================================================
SUBROUTINE MP_MHD_DISPL_2D

 ! CONSERVATIVE FORM OF THE FLUXES.

 ! I can't really do the conservative form - doesn't work for the incredibly small density 
 ! and forcing functions anyway. I don't know what the error is anyway
 implicit none
 integer kk 
 real(dp), dimension(dim1(rank), 1, dim3(rank)) :: flux1, flux2, flux3
 real(dp), dimension(dim1(rank),1,dim3(rank)) :: dxby,dxbz,dzbx,dzby
 real(dp), dimension(dim1(rank),dim3(rank),dim3(rank)) :: dzflux1,dzflux2,dxflux2,dxflux3
 real(dp), dimension(dim1(rank),1,dim3(rank)) :: term
integer bh_displ,bv_displ,bh_prerho,bv_prerho,bh_mag,bv_mag

 if (periodic) then
   bh_displ = 5
   bh_prerho = 5
   bh_mag = 5
 else
   bh_displ = 1
   bh_prerho = 0
   bh_mag = 1
 endif
 bv_displ = 1
 bv_prerho = 0
 bv_mag = 1

 if (TWOP5D == 1) flux1 = - (xi_z * boy - xi_y * boz)
 flux2 = - (xi_x * boz - xi_z * box)
 if (TWOP5D == 1) flux3 = - (xi_y * box - xi_x * boy)

 call ddx(xi_x,dxixdx,bh_displ)
 call ddz(xi_z,dxizdz,bv_displ)

  if (USE_PML) then
    if (IAM_ZBOT) then
      RHSpsizvz = az*dxizdz(:,:,1:nzpmlbot) + bzpml*psizvz
      dxizdz(:,:,1:nzpmlbot) = dxizdz(:,:,1:nzpmlbot)/kappaz + psizvz
    elseif (IAM_ZTOP) then
      RHSpsizvz = az*dxizdz(:,:,dim3(rank)-nzpmltop+1:dim3(rank)) + bzpml*psizvz
      dxizdz(:,:,dim3(rank)-nzpmltop+1:dim3(rank)) = &
dxizdz(:,:,dim3(rank)-nzpmltop+1:dim3(rank))/kappaz + psizvz
    endif
  endif

  if (USE_HPML) then
    if (IAM_XBOT) then
      RHSpsixvx = ax*dxixdx(1:nxpmlbot,:,:) + bxpml*psixvx
      dxixdx(1:nxpmlbot,:,:) = dxixdx(1:nxpmlbot,:,:)/kappax + psixvx
    elseif (IAM_XTOP) then
      RHSpsixvx = ax*dxixdx(dim1(rank)-nxpmltop+1:dim1(rank),:,:) + bxpml*psixvx
      dxixdx(dim1(rank)-nxpmltop+1:dim1(rank),:,:) = &
dxixdx(dim1(rank)-nxpmltop+1:dim1(rank),:,:)/kappax + psixvx
    endif
  endif

 div = dxixdx + dxizdz

 p = - c2rho0 * div - xi_z * gradp0_z - xi_x * gradp0_x
 rho = - rho0 * div - xi_z * gradrho0_z - xi_x * gradrho0_x

 call ddx(p,gradp_x,bh_prerho)
 call ddz(p,gradp_z,bv_prerho)

  if (USE_PML) then
    if (IAM_ZBOT) then
      RHSpsizp = az*gradp_z(:,:,1:nzpmlbot) + bzpml*psizp
      gradp_z(:,:,1:nzpmlbot) = gradp_z(:,:,1:nzpmlbot)/kappaz + psizp
    elseif (IAM_ZTOP) then
      RHSpsizp = az*gradp_z(:,:,dim3(rank)-nzpmltop+1:dim3(rank)) + bzpml*psizp
      gradp_z(:,:,dim3(rank)-nzpmltop+1:dim3(rank)) = &
gradp_z(:,:,dim3(rank)-nzpmltop+1:dim3(rank))/kappaz + psizp
    endif
  endif

  if (USE_HPML) then
    if (IAM_XBOT) then
      RHSpsixp = ax*gradp_x(1:nxpmlbot,:,:) + bxpml*psixp
      gradp_x(1:nxpmlbot,:,:) = gradp_x(1:nxpmlbot,:,:)/kappax + psixp
    elseif (IAM_XTOP) then
      RHSpsixp = ax*gradp_x(dim1(rank)-nxpmltop+1:dim1(rank),:,:) + bxpml*psixp
      gradp_x(dim1(rank)-nxpmltop+1:dim1(rank),:,:) = &
gradp_x(dim1(rank)-nxpmltop+1:dim1(rank),:,:)/kappax + psixp
    endif
  endif

  if (TWOP5D == 1) call ddx(flux3,dxflux3,bh_mag)
  call ddx(flux2,dxflux2,bh_mag)
  if (TWOP5D == 1) call ddz(flux1,dzflux1,bv_mag)
  call ddz(flux2,dzflux2,bv_mag)

  if (USE_PML) then
    if (IAM_ZBOT) then
    RHSpsizinductionbx = az*dzflux2(:,:,1:nzpmlbot) + bzpml*psizinductionbx
    dzflux2(:,:,1:nzpmlbot) = dzflux2(:,:,1:nzpmlbot)/kappaz + psizinductionbx
    if (TWOP5D == 1) RHSpsizinductionby = az*dzflux1(:,:,1:nzpmlbot) + bzpml*psizinductionby
    if (TWOP5D == 1) dzflux1(:,:,1:nzpmlbot) = dzflux1(:,:,1:nzpmlbot)/kappaz + psizinductionby
    elseif (IAM_ZTOP) then
    RHSpsizinductionbx = az*dzflux2(:,:,dim3(rank)-nzpmltop+1:dim3(rank)) + &
bzpml*psizinductionbx
    dzflux2(:,:,dim3(rank)-nzpmltop+1:dim3(rank)) = &
dzflux2(:,:,dim3(rank)-nzpmltop+1:dim3(rank))/kappaz + psizinductionbx
    if (TWOP5D == 1) RHSpsizinductionby = az*dzflux1(:,:,dim3(rank)-nzpmltop+1:dim3(rank)) + &
bzpml*psizinductionby
    if (TWOP5D == 1) dzflux1(:,:,dim3(rank)-nzpmltop+1:dim3(rank)) = &
    dzflux1(:,:,dim3(rank)-nzpmltop+1:dim3(rank))/kappaz + psizinductionby
    endif
  endif
  if (USE_HPML) then
    if (IAM_XBOT) then
      if (TWOP5D == 1) RHSpsixinductionby = ax*dxflux3(1:nxpmlbot,:,:) + bxpml*psixinductionby
      if (TWOP5D == 1) dxflux3(1:nxpmlbot,:,:) = dxflux3(1:nxpmlbot,:,:)/kappax + psixinductionby
      RHSpsixinductionbz = ax*dxflux2(1:nxpmlbot,:,:) + bxpml*psixinductionbz
      dxflux2(1:nxpmlbot,:,:) = dxflux2(1:nxpmlbot,:,:)/kappax + psixinductionbz
    elseif (IAM_XTOP) then
      if (TWOP5D == 1) RHSpsixinductionby = ax*dxflux3(dim1(rank)-nxpmltop+1:dim1(rank),:,:) + &
      bxpml*psixinductionby
      if (TWOP5D == 1) dxflux3(dim1(rank)-nxpmltop+1:dim1(rank),:,:) = &
      dxflux3(dim1(rank)-nxpmltop+1:dim1(rank),:,:)/kappax + psixinductionby
      RHSpsixinductionbz = ax*dxflux2(dim1(rank)-nxpmltop+1:dim1(rank),:,:) + &
      bxpml*psixinductionbz
      dxflux2(dim1(rank)-nxpmltop+1:dim1(rank),:,:) = &
      dxflux2(dim1(rank)-nxpmltop+1:dim1(rank),:,:)/kappax + psixinductionbz
    endif
  endif
 
 bx =   - dzflux2
 if (TWOP5D == 1) by = - dxflux3 + dzflux1
 bz =   dxflux2 

 if (TWOP5D == 1) call ddx(by,dxby, bh_mag)
 call ddx(bz,dxbz,bh_mag)
 call ddz(bx,dzbx,bv_mag)
 if (TWOP5D == 1) call ddz(by,dzby,bv_mag)

 ! NOTE THAT FLUX1 IS OVERWRITTEN IN THE FINAL CALL

    if (USE_PML) then
    if (IAM_ZBOT) then
    RHSpsizdzbx = az*dzbx(:,:,1:nzpmlbot) + bzpml*psizdzbx
    dzbx(:,:,1:nzpmlbot) = dzbx(:,:,1:nzpmlbot)/kappaz + psizdzbx
    if (TWOP5D == 1) RHSpsizdzby = az*dzby(:,:,1:nzpmlbot) + bzpml*psizdzby
    if (TWOP5D == 1) dzby(:,:,1:nzpmlbot) = dzby(:,:,1:nzpmlbot)/kappaz + psizdzby
    elseif (IAM_ZTOP) then
    RHSpsizdzbx = az*dzbx(:,:,dim3(rank)-nzpmltop+1:dim3(rank)) + bzpml*psizdzbx
    dzbx(:,:,dim3(rank)-nzpmltop+1:dim3(rank)) = &
    dzbx(:,:,dim3(rank)-nzpmltop+1:dim3(rank))/kappaz + psizdzbx
    if (TWOP5D == 1) RHSpsizdzby = az*dzby(:,:,dim3(rank)-nzpmltop+1:dim3(rank)) + bzpml*psizdzby
    if (TWOP5D == 1) dzby(:,:,dim3(rank)-nzpmltop+1:dim3(rank)) = &
    dzby(:,:,dim3(rank)-nzpmltop+1:dim3(rank))/kappaz + psizdzby
    endif
  endif
  if (USE_HPML) then
    if (IAM_XBOT) then
      if (TWOP5D == 1) RHSpsixdxby = ax*dxby(1:nxpmlbot,:,:) + bxpml*psixdxby
      if (TWOP5D == 1) dxby(1:nxpmlbot,:,:) = dxby(1:nxpmlbot,:,:)/kappax + psixdxby
      RHSpsixdxbz = ax*dxbz(1:nxpmlbot,:,:) + bxpml*psixdxbz
      dxbz(1:nxpmlbot,:,:) = dxbz(1:nxpmlbot,:,:)/kappax + psixdxbz
    elseif (IAM_XTOP) then
      if (TWOP5D == 1) RHSpsixdxby = ax*dxby(dim1(rank)-nxpmltop+1:dim1(rank),:,:) + &
bxpml*psixdxby
      if (TWOP5D == 1) dxby(dim1(rank)-nxpmltop+1:dim1(rank),:,:) = &
      dxby(dim1(rank)-nxpmltop+1:dim1(rank),:,:)/kappax + psixdxby
      RHSpsixdxbz = ax*dxbz(dim1(rank)-nxpmltop+1:dim1(rank),:,:) + &
bxpml*psixdxbz
      dxbz(dim1(rank)-nxpmltop+1:dim1(rank),:,:) = &
      dxbz(dim1(rank)-nxpmltop+1:dim1(rank),:,:)/kappax + psixdxbz
    endif
  endif

 curlbx =  - dzby
 if (TWOP5D == 1) curlby = -dxbz + dzbx
 curlbz = dxby 

 call cross(curlbx, curlby, curlbz, box, boy, boz, flux1, flux2, flux3)

 RHSv_x = rhoinv * (- gradp_x + flux1)
 if (TWOP5D == 1) RHSv_y = rhoinv * (- flux2)
 RHSv_z = rhoinv * (- gradp_z + flux3)

 call cross(curlbox, curlboy, curlboz, bx, by, bz, flux1, flux2, flux3)

 RHSv_x = RHSv_x + rhoinv*flux1 - spongexyz*v_x 
 if (TWOP5D == 1) RHSv_y = RHSv_y + rhoinv*flux2 - spongexyz*v_y
 RHSv_z = RHSv_z + rhoinv*flux3 - spongexyz*v_z

 do kk= 1,dim3(rank)
   RHSv_z(:,:,kk) = RHSv_z(:,:,kk) - rhoinv(:,:,kk) * rho(:,:,kk)*g(kk)
 enddo 

  if (FLOWS) then

  call ddx(v_x,dxixdx,bh_displ)
  RHSv_x = RHSv_x - 2.0_dp * v0_x * dxixdx

  call ddz(v_x,dxixdx,bh_displ)
  RHSv_x = RHSv_x - 2.0_dp * v0_z * dxixdx

  if (TWOP5D == 1) then
  call ddx(v_y,dxixdx,bh_displ)
  RHSv_x = RHSv_x - 2.0_dp * v0_y * dxixdx

  call ddz(v_y,dxixdx,bh_displ)
  RHSv_x = RHSv_x - 2.0_dp * v0_y * dxixdx

  endif

  call ddx(v_z,dxixdx,bv_displ)
  RHSv_z = RHSv_z - 2.0_dp * v0_x * dxixdx

  call ddz(v_z,dxixdx,bv_displ)
  RHSv_z = RHSv_z - 2.0_dp * v0_z * dxixdx

 endif


 RHSxi_x = v_x
 if (TWOP5D == 1) RHSxi_y = v_y
 RHSxi_z = v_z

 if (pulse) then
   if (pulse_dir == 1) then
     RHSv_x = RHSv_x + forcing_p
   elseif (pulse_dir == 2) then
     RHSv_y = RHSv_y + forcing_p
   elseif (pulse_dir == 3) then
     RHSv_z = RHSv_z + forcing_p
   endif
 else
   do kk = 1,dim3(rank)
     RHSv_z(:,:,kk) = RHSv_z(:,:,kk) +forcing*source_dep(kk)
   enddo
 endif
 if (IAM_NOTB) then
 else
   if (.not. PERIODIC) then
     if (IAM_XTOP) scr(dim1(rank),:,:,:) = 0.0_dp
     if (IAM_XBOT) scr(1,:,:,:) = 0.0_dp
   endif
   if (IAM_ZTOP) scr(:,:,dim3(rank),:) = 0.0_dp
   if (IAM_ZBOT) scr(:,:,1,:) = 0.0_dp
 endif

  if (DAMPING) then
  call damp_velocity(v_x,term)
  RHSv_x = RHSv_x - term
  if (TWOP5D == 1) call damp_velocity(v_y,term)
  if (TWOP5D == 1) RHSv_y = RHSv_y - term
  call damp_velocity(v_z,term)
  RHSv_z = RHSv_z - term
 endif

 END SUBROUTINE MP_MHD_DISPL_2D

!================================================================================================

SUBROUTINE MP_QUIET_2D

 ! CONSERVATIVE FORM OF THE FLUXES.

 ! I can't really do the conservative form - doesn't work for the incredibly small density 
 ! and forcing functions anyway. I don't know what the error is anyway
 integer kk
 integer bh_vel, bh_prerho, bv_vel,bv_prerho
 real(dp), dimension(dim1(rank),dim2(rank),dim3(rank)) :: term,gradrho_x,gradrho_z,&
        div0,dvxdz,dvzdx
 
 if (periodic) then
   bh_vel = 5
   bh_prerho = 5
 else
   bh_vel = 1
   bh_prerho = 0
 endif
 bv_vel = 1
 bv_prerho = 0

 call ddz(v_z, dvzdz, bv_vel)
 call ddz(p, gradp_z, bv_prerho)
 call ddx(v_x, dvxdx, bh_vel)
 call ddx(p, gradp_x, bh_prerho)

 if (USE_PML) then
    if (IAM_ZBOT) then
      RHSpsizp = az*gradp_z(:,:,1:nzpmlbot) + bzpml*psizp
      gradp_z(:,:,1:nzpmlbot) =  gradp_z(:,:,1:nzpmlbot) + psizp
      RHSpsizvz = az*dvzdz(:,:,1:nzpmlbot) + bzpml*psizvz
      dvzdz(:,:,1:nzpmlbot) = dvzdz(:,:,1:nzpmlbot) + psizvz
    endif
    if (IAM_ZTOP) then
      RHSpsizp = az*gradp_z(:,:,dim3(rank)-nzpmltop+1:dim3(rank)) + bzpml*psizp
      gradp_z(:,:,dim3(rank)-nzpmltop+1:dim3(rank)) =  gradp_z(:,:,dim3(rank)-nzpmltop+1:dim3(rank)) + psizp
      RHSpsizvz = az*dvzdz(:,:,dim3(rank)-nzpmltop+1:dim3(rank)) + bzpml*psizvz
      dvzdz(:,:,dim3(rank)-nzpmltop+1:dim3(rank)) = dvzdz(:,:,dim3(rank)-nzpmltop+1:dim3(rank)) + psizvz
    endif
  endif

  if (USE_HPML) then
    if (IAM_XBOT) then
      RHSpsixp = ax*gradp_x(1:nxpmlbot,:,:) + bxpml*psixp
      gradp_x(1:nxpmlbot,:,:) =  gradp_x(1:nxpmlbot,:,:) + psixp
      RHSpsixvx = ax*dvxdx(1:nxpmlbot,:,:) + bxpml*psixvx
      dvxdx(1:nxpmlbot,:,:) = dvxdx(1:nxpmlbot,:,:) + psixvx
    endif
    if (IAM_XTOP) then
      RHSpsixp = ax*gradp_x(dim1(rank)-nxpmltop+1:dim1(rank),:,:) + bxpml*psixp
      gradp_x(dim1(rank)-nxpmltop+1:dim1(rank),:,:) =  gradp_x(dim1(rank)-nxpmltop+1:dim1(rank),:,:) + psixp
      RHSpsixvx = ax*dvxdx(dim1(rank)-nxpmltop+1:dim1(rank),:,:) + bxpml*psixvx
      dvxdx(dim1(rank)-nxpmltop+1:dim1(rank),:,:) = dvxdx(dim1(rank)-nxpmltop+1:dim1(rank),:,:) + psixvx
    endif
  endif


 div = dvxdx + dvzdz

 RHScont =  - rho0 * div - gradrho0_z * v_z-spongexyz*rho
 RHSp = - c2rho0 * div - v_z * gradp0_z - spongexyz*p
 RHSv_x = - rhoinv * gradp_x -spongexyz*v_x

 do kk= 1,dim3(rank)
  RHSv_z(:,:,kk) =  - rhoinv(:,:,kk) * (gradp_z(:,:,kk) + rho(:,:,kk)*g(kk))-&
  spongexyz(:,:,kk)*v_z(:,:,kk)
 enddo

 if (FLOWS) then

 RHSp = RHSp - gradp_x * v0_x &
           - gradp_z * v0_z

 call ddz(rho,gradrho_z,bv_prerho)
 call ddx(rho,gradrho_x,bh_prerho)

 RHScont = RHScont - gradrho_x * v0_x &
           - gradrho_z * v0_z - rho * div0

  call ddz(v_x,dvxdz,bv_vel)
  call ddx(v_z, dvzdx,bh_vel)
  RHSv_x = RHSv_x - 2.0_dp * v0_x * dvxdx - 2.0_dp * v0_z * dvxdz - &
2.0_dp * v_x * dv0xdx - 2.0_dp * v_z * dv0xdz
  RHSv_z = RHSv_z - 2.0_dp * v0_z * dvzdz - 2.0_dp * v0_x * dvzdx - &
2.0_dp * v_z * dv0zdz - 2.0_dp* v_x * dv0zdx

 endif

 if (pulse) then
   if (pulse_dir == 1) then
     RHSv_x = RHSv_x + forcing_p
   elseif (pulse_dir == 3) then
     RHSv_z = RHSv_z + forcing_p
   endif
 else
   do kk = 1,dim3(rank)
     RHSv_z(:,:,kk) = RHSv_z(:,:,kk) + forcing*source_dep(kk)
   enddo
 endif

  if (IAM_NOTB) then
 else
   if (.not. PERIODIC) then
     if (IAM_XTOP) scr(dim1(rank),:,:,:) = 0.0_dp
     if (IAM_XBOT) scr(1,:,:,:) = 0.0_dp
   endif
   if (IAM_ZTOP) scr(:,:,dim3(rank),:) = 0.0_dp
   if (IAM_ZBOT) scr(:,:,1,:) = 0.0_dp
 endif

  if (DAMPING) then
  call damp_velocity(v_x,term)
  RHSv_x = RHSv_x - term
  call damp_velocity(v_z,term)
  RHSv_z = RHSv_z - term
 endif

END SUBROUTINE MP_QUIET_2D

!================================================================================================

SUBROUTINE MP_QUIET_DISPL_2D

 ! CONSERVATIVE FORM OF THE FLUXES.

 ! I can't really do the conservative form - doesn't work for the incredibly small density 
 ! and forcing functions anyway. I don't know what the error is anyway
 integer kk
 integer bh_displ,bv_displ,bh_prerho,bv_prerho
  real(dp), dimension(dim1(rank),1,dim3(rank)) :: term

 if (periodic) then
   bh_displ = 5
   bh_prerho = 5
 else
   bh_displ = 1
   bh_prerho = 0
 endif
 bv_displ = 1
 bv_prerho = 0

 call ddz(xi_z, dxizdz, bv_displ)
 call ddx(xi_x, dxixdx, bh_displ)

  if (USE_PML) then
    if (IAM_ZBOT) then
      RHSpsizvz = az*dxizdz(:,:,1:nzpmlbot) + bzpml*psizvz
      dxizdz(:,:,1:nzpmlbot) = dxizdz(:,:,1:nzpmlbot) + psizvz
    endif
    if (IAM_ZTOP) then
      RHSpsizvz = az*dxizdz(:,:,dim3(rank)-nzpmlbot+1:dim3(rank)) + bzpml*psizvz
      dxizdz(:,:,dim3(rank)-nzpmlbot+1:dim3(rank)) = &
dxizdz(:,:,dim3(rank)-nzpmlbot+1:dim3(rank)) + psizvz
    endif
  endif

  if (USE_HPML) then
    if (IAM_XBOT) then
      RHSpsixvx = ax*dxixdx(1:nxpmlbot,:,:) + bxpml*psixvx
      dxixdx(1:nxpmlbot,:,:) = dxixdx(1:nxpmlbot,:,:) + psixvx
    endif
    if (IAM_XTOP) then
      RHSpsixvx = ax*dxixdx(dim1(rank)-nxpmltop+1:dim1(rank),:,:) + bxpml*psixvx
      dxixdx(dim1(rank)-nxpmltop+1:dim1(rank),:,:) = &
dxixdx(dim1(rank)-nxpmltop+1:dim1(rank),:,:) + psixvx
    endif
  endif
 
 ! Compute divergence

 div = dxixdx + dxizdz

 ! RHS continuity and pressure

 p = - c2rho0 * div - gradp0_z *xi_z - gradp0_x*xi_x
 
 call ddx(p, gradp_x, bh_prerho)
 call ddz(p, gradp_z, bv_prerho)


  if (USE_PML) then
    if (IAM_ZBOT) then
      RHSpsizp = az*gradp_z(:,:,1:nzpmlbot) + bzpml*psizp
      gradp_z(:,:,1:nzpmlbot) = gradp_z(:,:,1:nzpmlbot) + psizp
    endif
    if (IAM_ZTOP) then
      RHSpsizp = az*gradp_z(:,:,dim3(rank)-nzpmlbot+1:dim3(rank)) + bzpml*psizp
      gradp_z(:,:,dim3(rank)-nzpmlbot+1:dim3(rank)) = &
gradp_z(:,:,dim3(rank)-nzpmlbot+1:dim3(rank)) + psizp
    endif
  endif

  if (USE_HPML) then
    if (IAM_XBOT) then
      RHSpsixp = ax*gradp_x(1:nxpmlbot,:,:) + bxpml*psixp
      gradp_x(1:nxpmlbot,:,:) = gradp_x(1:nxpmlbot,:,:) + psixp
    endif
    if (IAM_XTOP) then
      RHSpsixp = ax*gradp_x(dim1(rank)-nxpmltop+1:dim1(rank),:,:) + bxpml*psixp
      gradp_x(dim1(rank)-nxpmltop+1:dim1(rank),:,:) = &
gradp_x(dim1(rank)-nxpmltop+1:dim1(rank),:,:) + psixp
    endif
  endif

 RHSv_x = - rhoinv * gradp_x
 rho =  - rho0 * div - gradrho0_z * xi_z
 do kk= 1,dim3(rank)
  RHSv_z(:,:,kk) =  - rhoinv(:,:,kk) * (gradp_z(:,:,kk) + rho(:,:,kk)*g(kk))
 enddo

 RHSxi_x = v_x
 RHSxi_z = v_z

 if (pulse) then
   if (pulse_dir == 1) then
     RHSv_x = RHSv_x + forcing_p
   elseif (pulse_dir == 3) then
     RHSv_z = RHSv_z + forcing_p
   endif
 else
   do kk = 1,dim3(rank)
     RHSv_z(:,:,kk) = RHSv_z(:,:,kk) + forcing*source_dep(kk)
   enddo
 endif

 if (FLOWS) then

  call ddx(v_x,dxixdx,bh_displ)
  RHSv_x = RHSv_x - 2.0_dp * v0_x * dxixdx

  call ddz(v_x,dxixdx,bh_displ)
  RHSv_x = RHSv_x - 2.0_dp * v0_z * dxixdx
 
  call ddx(v_z,dxixdx,bv_displ)
  RHSv_z = RHSv_z - 2.0_dp * v0_x * dxixdx

  call ddz(v_z,dxixdx,bv_displ)
  RHSv_z = RHSv_z - 2.0_dp * v0_z * dxixdx

 endif

 if (IAM_NOTB) then
 else
   if (.not. PERIODIC) then
     if (IAM_XTOP) scr(dim1(rank),:,:,:) = 0.
     if (IAM_XBOT) scr(1,:,:,:) = 0.
   endif
   if (IAM_ZTOP) scr(:,:,dim3(rank),:) = 0.
   if (IAM_ZBOT) scr(:,:,1,:) = 0.
 endif

 if (DAMPING) then
  call damp_velocity(v_x,term)
  RHSv_x = RHSv_x - term
  call damp_velocity(v_z,term)
  RHSv_z = RHSv_z - term
 endif

 END SUBROUTINE MP_QUIET_DISPL_2D


!================================================================================================

END MODULE PHYSICS_RHS_2D
