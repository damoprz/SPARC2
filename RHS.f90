module PHYSICS_RHS

!! THE RHS for Quiet + PML routines
!! A simplified version of physics.f90 and pml.f90 from
!! Shravans SPARC. All initialization has been moved to 
!! init.f90

  use all_modules
  use derivatives

  implicit none
  
  contains

!==================================================================================


subroutine MP_QUIET_3D

! CONSERVATIVE FORM OF THE FLUXES.

! I can't really do the conservative form - doesn't work for the incredibly small density 
! and forcing functions anyway. I don't know what the error is anyway

  integer bh_vel,bv_vel,bh_prerho,bv_prerho
  integer kk
  real(dp), dimension(dim1(rank),dim2(rank),dim3(rank)) :: term,gradrho_x,gradrho_y,gradrho_z,&
        div0,dvxdy,dvxdz,dvydx,dvydz,dvzdx,dvzdy
  
  if (periodic) then
    bh_vel = 5
    bh_prerho = 5
  else
    bh_vel = 1 
    bh_prerho = 1
  endif
  bv_vel = 1
  bv_prerho = 1

! CALL THE OVERLAP-DERIVATIVE ROUTINE

  call ddxyz(v_x, dvxdx, v_y, dvydy, v_z, dvzdz,bh_vel,bh_vel,bv_vel)
  call ddxyz(p, gradp_x, p, gradp_y, p, gradp_z,bh_prerho,bh_prerho,bv_prerho)

! NOTE THAT flux3 is overwritten in the final ddxyz call

  if (USE_PML) then
    if (IAM_ZBOT) then
      RHSpsizp = az*gradp_z(:,:,1:nzpmlbot) + bzpml*psizp
      gradp_z(:,:,1:nzpmlbot) =  gradp_z(:,:,1:nzpmlbot)/kappaz + psizp
      RHSpsizvz = az*dvzdz(:,:,1:nzpmlbot) + bzpml*psizvz
      dvzdz(:,:,1:nzpmlbot) = dvzdz(:,:,1:nzpmlbot)/kappaz + psizvz
    elseif (IAM_ZTOP) then
      RHSpsizp = az*gradp_z(:,:,dim3(rank)-nzpmltop+1:dim3(rank)) + bzpml*psizp
      gradp_z(:,:,dim3(rank)-nzpmltop+1:dim3(rank)) =  gradp_z(:,:,dim3(rank)-nzpmltop+1:dim3(rank))/kappaz + psizp
      RHSpsizvz = az*dvzdz(:,:,dim3(rank)-nzpmltop+1:dim3(rank)) + bzpml*psizvz
      dvzdz(:,:,dim3(rank)-nzpmltop+1:dim3(rank)) = dvzdz(:,:,dim3(rank)-nzpmltop+1:dim3(rank))/kappaz + psizvz
    endif
  endif

  if (USE_HPML) then
    if (IAM_XBOT) then
      RHSpsixp = ax*gradp_x(1:nxpmlbot,:,:) + bxpml*psixp
      gradp_x(1:nxpmlbot,:,:) =  gradp_x(1:nxpmlbot,:,:)/kappax + psixp
      RHSpsixvx = ax*dvxdx(1:nxpmlbot,:,:) + bxpml*psixvx
      dvxdx(1:nxpmlbot,:,:) = dvxdx(1:nxpmlbot,:,:)/kappax + psixvx
    elseif (IAM_XTOP) then
      RHSpsixp = ax*gradp_x(dim1(rank)-nxpmltop+1:dim1(rank),:,:) + bxpml*psixp
      gradp_x(dim1(rank)-nxpmltop+1:dim1(rank),:,:) =  gradp_x(dim1(rank)-nxpmltop+1:dim1(rank),:,:)/kappax &
      + psixp
      RHSpsixvx = ax*dvxdx(dim1(rank)-nxpmltop+1:dim1(rank),:,:) + bxpml*psixvx
      dvxdx(dim1(rank)-nxpmltop+1:dim1(rank),:,:) = dvxdx(dim1(rank)-nxpmltop+1:dim1(rank),:,:)/kappax &
      + psixvx
    endif

    if (IAM_YBOT) then
      RHSpsiyp = ay*gradp_y(:,1:nypmlbot,:) + bypml*psiyp
      gradp_y(:,1:nypmlbot,:) =  gradp_y(:,1:nypmlbot,:)/kappay + psiyp
      RHSpsiyvy = ay*dvydy(:,1:nypmlbot,:) + bypml*psiyvy
      dvydy(:,1:nypmlbot,:) = dvydy(:,1:nypmlbot,:)/kappay + psiyvy
    elseif (IAM_YTOP) then
      RHSpsiyp = ay*gradp_y(:,dim2(rank)-nypmltop+1:dim2(rank),:) + bypml*psiyp
      gradp_y(:,dim2(rank)-nypmltop+1:dim2(rank),:) =  gradp_y(:,dim2(rank)-nypmltop+1:dim2(rank),:)/kappay + &
      psiyp
      RHSpsiyvy = ay*dvydy(:,dim2(rank)-nypmltop+1:dim2(rank),:) + bypml*psiyvy
      dvydy(:,dim2(rank)-nypmltop+1:dim2(rank),:) = dvydy(:,dim2(rank)-nypmltop+1:dim2(rank),:)/kappay + &
      psiyvy
    endif

  endif
  
  div = dvxdx + dvydy + dvzdz

  RHScont = - gradrho0_z * v_z - rho0 * div - spongexyz*rho

  RHSv_x = - rhoinv * gradp_x - spongexyz * v_x

  RHSv_y = - rhoinv * gradp_y - spongexyz * v_y

 if (FLOWS) then

 RHSp = RHSp - gradp_x * v0_x &
 -gradp_y*v0_y - gradp_z * v0_z

 call ddz(rho,gradrho_z,bv_prerho)
 call ddx(rho,gradrho_x,bh_prerho)
 call ddy(rho,gradrho_y,bh_prerho)

 RHScont = RHScont - gradrho_x * v0_x &
  - gradrho_y*v0_y - gradrho_z * v0_z - rho * div0

  call ddz(v_x,dvxdz,bv_vel)
  call ddx(v_z, dvzdx,bh_vel)
  call ddz(v_y,dvydz,bv_vel)
  call ddx(v_y, dvydx,bh_vel)
  call ddy(v_z, dvzdy,bh_vel)
  call ddy(v_x, dvxdy,bh_vel)

  RHSv_x = RHSv_x - 2.0_dp*(v0_x*dvxdx + v0_y*dvxdy + v0_z*dvxdz + &
  v_x * dv0xdx + v_y*dv0xdy + v_z * dv0xdz)

  RHSv_z = RHSv_z - 2.0_dp*(v0_x*dvzdx + v0_y*dvzdy + v0_z*dvzdz + &
  v_x * dv0zdx + v_y*dv0zdy + v_z * dv0zdz)

  RHSv_y = RHSv_z - 2.0_dp*(v0_x*dvydx + v0_y*dvydy + v0_z*dvydz + &
  v_x * dv0ydx + v_y*dv0ydy + v_z * dv0ydz)

 endif


  do kk= 1,dim3(rank)

    RHSv_z(:,:,kk) =  - rhoinv(:,:,kk) * (gradp_z(:,:,kk) + rho(:,:,kk)*g(kk)) - v_z(:,:,kk)*spongexyz(:,:,kk)

    if (pulse) then
      if (pulse_dir == 1) then
        RHSv_x(:,:,kk) = RHSv_x(:,:,kk) + forcing_p(:,:,kk)
      elseif (pulse_dir == 2) then
        RHSv_y(:,:,kk) = RHSv_y(:,:,kk) + forcing_p(:,:,kk)
      elseif (pulse_dir == 3) then
        RHSv_z(:,:,kk) = RHSv_z(:,:,kk) + forcing_p(:,:,kk)
      endif
    else
      RHSv_z(:,:,kk) = RHSv_z(:,:,kk) +  forcing*source_dep(kk)
    endif

 enddo

 RHSp = - c2rho0 * div - v_z * gradp0_z - spongexyz*p

 if (IAM_NOTB) then
 else 
   if (.not. PERIODIC) then
     if (IAM_XTOP) scr(dim1(rank),:,:,:) = 0.0_dp
     if (IAM_XBOT) scr(1,:,:,:) = 0.0_dp 
     if (IAM_YTOP) scr(:,dim2(rank),:,:) = 0.0_dp
     if (IAM_YBOT) scr(:,1,:,:) = 0.0_dp 
   endif
   if (IAM_ZTOP) scr(:,:,dim3(rank),:) = 0.0_dp 
   if (IAM_ZBOT) scr(:,:,1,:) = 0.0_dp
 endif
 
 if (DAMPING) then
  call damp_velocity(v_x,term)
  RHSv_x = RHSv_x - term
  call damp_velocity(v_y,term)
  RHSv_y = RHSv_y - term
  call damp_velocity(v_z,term)
  RHSv_z = RHSv_z - term
 endif
 
end subroutine MP_QUIET_3D
!================================================================================================

subroutine MP_QUIET_DISPL_3D

 ! CONSERVATIVE FORM OF THE FLUXES.

 ! I can't really do the conservative form - doesn't wo rk for the incredibly small density 
 ! and forcing functions anyway. I don't know what the error is anyway
 implicit none
 integer bh_disp,bv_disp,bh_prerho,bv_prerho
 integer kk
 real(dp), allocatable, dimension(:,:,:) ::temp_x,temp_y, temp_z
 real(dp), dimension(dim1(rank),dim2(rank),dim3(rank)) :: term

  if (periodic) then
   bh_disp = 5
   bh_prerho = 5
 else
   bh_disp = 1
   bh_prerho = 1
 endif
 bv_disp = 1
 bv_prerho = 1

 call ddxyz(xi_x, dxixdx, xi_y, dxiydy, xi_z, dxizdz, bh_disp,bh_disp,bv_disp)


  if (USE_PML) then
    if (IAM_ZBOT) then
      dxizdz(:,:,1:nzpmlbot) = dxizdz(:,:,1:nzpmlbot)/kappaz + psizvz
      RHSpsizvz = az*dxizdz(:,:,1:nzpmlbot) + bzpml*psizvz
    endif
    if (IAM_ZTOP) then
      dxizdz(:,:,dim3(rank)-nzpmltop+1:dim3(rank)) = &
dxizdz(:,:,dim3(rank)-nzpmltop+1:dim3(rank))/kappaz + psizvz
      RHSpsizvz = az*dxizdz(:,:,dim3(rank)-nzpmltop+1:dim3(rank)) + bzpml*psizvz
    endif
  endif

  if (USE_HPML) then
    if (IAM_XBOT) then
      dxixdx(1:nxpmlbot,:,:) = dxixdx(1:nxpmlbot,:,:)/kappax + psixvx
      RHSpsixvx = ax*dxixdx(1:nxpmlbot,:,:) + bxpml*psixvx
    endif
    if (IAM_XTOP) then
      dxixdx(dim1(rank)-nxpmltop+1:dim1(rank),:,:) = &
dxixdx(dim1(rank)-nxpmltop+1:dim1(rank),:,:)/kappax + psixvx
      RHSpsixvx = ax*dxixdx(dim1(rank)-nxpmltop+1:dim1(rank),:,:) + bxpml*psixvx
    endif

    if (IAM_YBOT) then
      dxiydy(:,1:nypmlbot,:) = dxiydy(:,1:nypmlbot,:)/kappay + psiyvy
      RHSpsiyvy = ay*dxiydy(:,1:nypmlbot,:) + bypml*psiyvy
    endif

    if (IAM_YTOP) then
      dxiydy(:,dim2(rank)-nypmltop+1:dim2(rank),:) = &
dxiydy(:,dim2(rank)-nypmltop+1:dim2(rank),:)/kappay + psiyvy
      RHSpsiyvy = ay*dxiydy(:,dim2(rank)-nypmltop+1:dim2(rank),:) + bypml*psiyvy
    endif

  endif

 div = dxixdx + dxiydy + dxizdz

 p = - c2rho0 * div - xi_z * gradp0_z
 rho = - rho0 * div - xi_z * gradrho0_z

 call ddxyz(p, gradp_x, p, gradp_y, p, gradp_z, bh_prerho,bh_prerho,bv_prerho)

  if (USE_PML) then
    if (IAM_ZBOT) then
      RHSpsizp = az*gradp_z(:,:,1:nzpmlbot) + bzpml*psizp
      gradp_z(:,:,1:nzpmlbot) = gradp_z(:,:,1:nzpmlbot)/kappaz + psizp
    endif
    if (IAM_ZTOP) then
      RHSpsizp = az*gradp_z(:,:,dim3(rank)-nzpmltop+1:dim3(rank)) + bzpml*psizp
      gradp_z(:,:,dim3(rank)-nzpmltop+1:dim3(rank)) = &
gradp_z(:,:,dim3(rank)-nzpmltop+1:dim3(rank))/kappaz + psizp
    endif
  endif

  if (USE_HPML) then
    if (IAM_XBOT) then
      RHSpsixp = ax*gradp_x(1:nxpmlbot,:,:) + bxpml*psixp
      gradp_x(1:nxpmlbot,:,:) = gradp_x(1:nxpmlbot,:,:)/kappax + psixp
    endif
    if (IAM_XTOP) then
      RHSpsixp = ax*gradp_x(dim1(rank)-nxpmltop+1:dim1(rank),:,:) + bxpml*psixp
      gradp_x(dim1(rank)-nxpmltop+1:dim1(rank),:,:) = &
gradp_x(dim1(rank)-nxpmltop+1:dim1(rank),:,:)/kappax + psixp
    endif
 
    if (IAM_YBOT) then
      RHSpsiyp = ay*gradp_y(:,1:nypmlbot,:) + bypml*psiyp
      gradp_y(:,1:nypmlbot,:) = gradp_y(:,1:nypmlbot,:)/kappay + psiyp
    endif

    if (IAM_YTOP) then
      RHSpsiyp = ay*gradp_y(:,dim2(rank)-nypmltop+1:dim2(rank),:) + bypml*psiyp
      gradp_y(:,dim2(rank)-nypmltop+1:dim2(rank),:) = &
gradp_y(:,dim2(rank)-nypmltop+1:dim2(rank),:)/kappay + psiyp
    endif

  endif

  RHSv_x = - rhoinv * gradp_x - spongexyz*v_x
  RHSv_y = - rhoinv * gradp_y - spongexyz*v_y
  RHSv_z = - rhoinv * gradp_z - spongexyz*v_z

if (FLOWS) then

  allocate(temp_x(dim1(rank),dim2(rank),dim3(rank)),temp_y(dim1(rank),dim2(rank),dim3(rank)),&
   temp_z(dim1(rank),dim2(rank),dim3(rank)))

   call ddxyz(v_x, temp_x, v_x, temp_y, v_x, temp_z, bh_disp,bh_disp,bv_disp)

   RHSv_x = RHSv_x- 2.0_dp * (v0_x * temp_x + v0_y * temp_y + v0_z * temp_z)

   call ddxyz(v_y, temp_x, v_y, temp_y, v_y, temp_z, bh_disp,bh_disp,bv_disp)

   RHSv_y = RHSv_y- 2.0_dp * (v0_x * temp_x + v0_y * temp_y + v0_z * temp_z)

   call ddxyz(v_z, temp_x, v_z, temp_y, v_z, temp_z, bh_disp,bh_disp,bv_disp)

   RHSv_z = RHSv_z- 2.0_dp * (v0_x * temp_x + v0_y * temp_y + v0_z * temp_z)

   deallocate(temp_x, temp_y, temp_z)
endif

 do kk= 1,dim3(rank)
   RHSv_z(:,:,kk) = RHSv_z(:,:,kk) - rho(:,:,kk)*g(kk)*rhoinv(:,:,kk)
  enddo

 RHSxi_x = v_x
 RHSxi_y = v_y
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
     RHSv_z(:,:,kk) = RHSv_z(:,:,kk) + forcing*source_dep(kk)
   enddo
 endif

 if (IAM_NOTB) then
 else
   if (.not. PERIODIC) then
     if (IAM_XTOP) scr(dim1(rank),:,:,:) = 0.0_dp
     if (IAM_XBOT) scr(1,:,:,:) = 0.0_dp
     if (IAM_YTOP) scr(:,dim2(rank),:,:) = 0.0_dp
     if (IAM_YBOT) scr(:,1,:,:) = 0.0_dp
   endif
   if (IAM_ZTOP) scr(:,:,dim3(rank),:) = 0.0_dp
   if (IAM_ZBOT) scr(:,:,1,:) = 0.0_dp
 endif

 if (DAMPING) then
  call damp_velocity(v_x,term)
  RHSv_x = RHSv_x - term
  call damp_velocity(v_y,term)
  RHSv_y = RHSv_y - term
  call damp_velocity(v_z,term)
  RHSv_z = RHSv_z - term
 endif

 end subroutine MP_QUIET_DISPL_3D

!================================================================================================

 subroutine MP_MHD_3D

 real(dp), dimension(dim1(rank),dim2(rank),dim3(rank)) :: flux1,flux2,flux3
 real(dp), dimension(dim1(rank),dim2(rank),dim3(rank)) :: dzflux1,dzflux2,dxflux2,dxflux3,dyflux1,dyflux3
 real(dp), dimension(dim1(rank),dim2(rank),dim3(rank)) :: dxbz,dxby,dybx,dybz,dzbx,dzby
 real(dp), dimension(dim1(rank),dim2(rank),dim3(rank)) :: term,gradrho_x,gradrho_y,gradrho_z,&
        div0,dvxdy,dvxdz,dvydx,dvydz,dvzdx,dvzdy
 integer kk
 integer bh_vel,bv_vel,bh_mag,bv_mag,bh_prerho,bv_prerho


 if (periodic) then
 bh_vel = 5
 bv_vel = 1
 bh_mag = 5
 bv_mag = 1
 bh_prerho = 5
 bv_prerho = 1
 else
 bh_vel = 1
 bv_vel = 1
 bh_mag = 1
 bv_mag = 1
 bh_prerho = 1
 bv_prerho = 1
 endif
 ! CONSERVATIVE FORM OF THE FLUXES.

 ! I can't really do the conservative form - doesn't work for the incredibly small density 
 ! and forcing functions anyway. I don't know what the error is anyway

 flux1 = - v_z * boy + v_y * boz
 flux2 = - v_x * boz + v_z * box
 flux3 = - v_y * box + v_x * boy

 ! CALLING OVERLAP DERIVATIVE ROUTINE

 call ddxyz(v_x, dvxdx, v_y, dvydy, v_z, dvzdz,bh_vel,bh_vel,bv_vel)
 call ddxyz(flux3, dxflux3, p, gradp_y, flux1, dzflux1,bh_mag,bh_prerho,bv_mag)
 call ddxyz(p, gradp_x, flux3, dyflux3, flux2, dzflux2,bh_prerho,bh_mag,bv_mag)
 call ddxyz(flux2, dxflux2, bz, dybz,bx, dzbx,bh_mag,bh_mag,bv_mag)
 call ddxyz(by, dxby, bx, dybx, p, gradp_z,bh_mag,bh_mag,bv_prerho)
 call ddxyz(bz, dxbz, flux1, dyflux1, by, dzby,bh_mag,bh_mag,bv_mag)

  if (USE_PML) then
    if (IAM_ZBOT) then
      RHSpsizp = az*gradp_z(:,:,1:nzpmlbot) + bzpml*psizp
      gradp_z(:,:,1:nzpmlbot) =  gradp_z(:,:,1:nzpmlbot)/kappaz + psizp
      RHSpsizvz = az*dvzdz(:,:,1:nzpmlbot) + bzpml*psizvz
      dvzdz(:,:,1:nzpmlbot) = dvzdz(:,:,1:nzpmlbot)/kappaz + psizvz
    elseif (IAM_ZTOP) then
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
    elseif (IAM_XTOP) then
      RHSpsixp = ax*gradp_x(dim1(rank)-nxpmltop+1:dim1(rank),:,:) + bxpml*psixp
      gradp_x(dim1(rank)-nxpmltop+1:dim1(rank),:,:) = &
gradp_x(dim1(rank)-nxpmltop+1:dim1(rank),:,:)/kappax + psixp
      RHSpsixvx = ax*dvxdx(dim1(rank)-nxpmltop+1:dim1(rank),:,:) + bxpml*psixvx
      dvxdx(dim1(rank)-nxpmltop+1:dim1(rank),:,:) = &
dvxdx(dim1(rank)-nxpmltop+1:dim1(rank),:,:)/kappax + psixvx
    endif
    if (IAM_YBOT) then
      RHSpsiyp = ay*gradp_y(:,1:nypmlbot,:) + bypml*psiyp
      gradp_y(:,1:nypmlbot,:) =  gradp_y(:,1:nypmlbot,:)/kappay + psiyp
      RHSpsiyvy = ay*dvydy(:,1:nypmlbot,:) + bypml*psiyvy
      dvydy(:,1:nypmlbot,:) = dvydy(:,1:nypmlbot,:)/kappay + psiyvy
    elseif (IAM_YTOP) then
      RHSpsiyp = ay*gradp_y(:,dim2(rank)-nypmltop+1:dim2(rank),:) + bypml*psiyp
      gradp_y(:,dim2(rank)-nypmltop+1:dim2(rank),:) = &
gradp_y(:,dim2(rank)-nypmltop+1:dim2(rank),:)/kappay + psiyp
      RHSpsiyvy = ay*dvydy(:,dim2(rank)-nypmltop+1:dim2(rank),:) + bypml*psiyvy
      dvydy(:,dim2(rank)-nypmltop+1:dim2(rank),:) = &
dvydy(:,dim2(rank)-nypmltop+1:dim2(rank),:)/kappay + psiyvy
    endif

  endif

  if (USE_PML) then
    if (IAM_ZBOT) then
    RHSpsizinductionbx = az*dzflux2(:,:,1:nzpmlbot) + bzpml*psizinductionbx
    dzflux2(:,:,1:nzpmlbot) = dzflux2(:,:,1:nzpmlbot)/kappaz + psizinductionbx
    RHSpsizinductionby = az*dzflux1(:,:,1:nzpmlbot) + bzpml*psizinductionby
    dzflux1(:,:,1:nzpmlbot) = dzflux1(:,:,1:nzpmlbot)/kappaz + psizinductionby
    RHSpsizdzbx = az*dzbx(:,:,1:nzpmlbot) + bzpml*psizdzbx
    dzbx(:,:,1:nzpmlbot) = dzbx(:,:,1:nzpmlbot)/kappaz + psizdzbx
    RHSpsizdzby = az*dzby(:,:,1:nzpmlbot) + bzpml*psizdzby
    dzby(:,:,1:nzpmlbot) = dzby(:,:,1:nzpmlbot)/kappaz + psizdzby
    elseif (IAM_ZTOP) then
    RHSpsizinductionbx = az*dzflux2(:,:,dim3(rank)-nzpmltop+1:dim3(rank)) + bzpml*psizinductionbx
    dzflux2(:,:,dim3(rank)-nzpmltop+1:dim3(rank)) = &
dzflux2(:,:,dim3(rank)-nzpmltop+1:dim3(rank))/kappaz + psizinductionbx
    RHSpsizinductionby = az*dzflux1(:,:,dim3(rank)-nzpmltop+1:dim3(rank)) + bzpml*psizinductionby
    dzflux1(:,:,dim3(rank)-nzpmltop+1:dim3(rank)) = &
    dzflux1(:,:,dim3(rank)-nzpmltop+1:dim3(rank))/kappaz + psizinductionby
    RHSpsizdzbx = az*dzbx(:,:,dim3(rank)-nzpmltop+1:dim3(rank)) + bzpml*psizdzbx
    dzbx(:,:,dim3(rank)-nzpmltop+1:dim3(rank)) = &
    dzbx(:,:,dim3(rank)-nzpmltop+1:dim3(rank))/kappaz + psizdzbx
    RHSpsizdzby =az*dzby(:,:,dim3(rank)-nzpmltop+1:dim3(rank)) + bzpml*psizdzby
    dzby(:,:,dim3(rank)-nzpmltop+1:dim3(rank)) = &
    dzby(:,:,dim3(rank)-nzpmltop+1:dim3(rank))/kappaz + psizdzby
    endif
  endif
  if (USE_HPML) then
    if (IAM_XBOT) then
      RHSpsixinductionby = ax*dxflux3(1:nxpmlbot,:,:) + bxpml*psixinductionby
      dxflux3(1:nxpmlbot,:,:) = dxflux3(1:nxpmlbot,:,:)/kappax + psixinductionby
      RHSpsixinductionbz = ax*dxflux2(1:nxpmlbot,:,:) + bxpml*psixinductionbz
      dxflux2(1:nxpmlbot,:,:) = dxflux2(1:nxpmlbot,:,:)/kappax + psixinductionbz
      RHSpsixdxby = ax*dxby(1:nxpmlbot,:,:) + bxpml*psixdxby
      dxby(1:nxpmlbot,:,:) = dxby(1:nxpmlbot,:,:)/kappax + psixdxby
      RHSpsixdxbz = ax*dxbz(1:nxpmlbot,:,:) + bxpml*psixdxbz
      dxbz(1:nxpmlbot,:,:) = dxbz(1:nxpmlbot,:,:)/kappax + psixdxbz
    elseif (IAM_XTOP) then
      RHSpsixinductionby = ax*dxflux3(dim1(rank)-nxpmltop+1:dim1(rank),:,:) + &
      bxpml*psixinductionby
      dxflux3(dim1(rank)-nxpmltop+1:dim1(rank),:,:) = &
      dxflux3(dim1(rank)-nxpmltop+1:dim1(rank),:,:)/kappax + psixinductionby
      RHSpsixinductionbz = ax*dxflux2(dim1(rank)-nxpmltop+1:dim1(rank),:,:) + &
      bxpml*psixinductionbz
      dxflux2(dim1(rank)-nxpmltop+1:dim1(rank),:,:) = &
      dxflux2(dim1(rank)-nxpmltop+1:dim1(rank),:,:)/kappax + psixinductionbz
      RHSpsixdxby = ax*dxby(dim1(rank)-nxpmltop+1:dim1(rank),:,:) + bxpml*psixdxby
      dxby(dim1(rank)-nxpmltop+1:dim1(rank),:,:) = &
      dxby(dim1(rank)-nxpmltop+1:dim1(rank),:,:)/kappax + psixdxby
      RHSpsixdxbz = ax*dxbz(dim1(rank)-nxpmltop+1:dim1(rank),:,:) + bxpml*psixdxbz
      dxbz(dim1(rank)-nxpmltop+1:dim1(rank),:,:) = &
      dxbz(dim1(rank)-nxpmltop+1:dim1(rank),:,:)/kappax + psixdxbz
    endif    
    if (IAM_YBOT) then
      RHSpsiyinductionbx = ay*dyflux3(:,1:nypmlbot,:) + bypml*psiyinductionbx
      dyflux3(:,1:nypmlbot,:) = dyflux3(:,1:nypmlbot,:)/kappay + psiyinductionbx
      RHSpsiyinductionbz = ay*dyflux1(:,1:nypmlbot,:) + bypml*psiyinductionbz
      dyflux1(:,1:nypmlbot,:) = dyflux1(:,1:nypmlbot,:)/kappay + psiyinductionbz
      RHSpsiydybx = ay*dybx(:,1:nypmlbot,:) + bypml*psiydybx
      dybx(:,1:nypmlbot,:) = dybx(:,1:nypmlbot,:)/kappay + psiydybx
      RHSpsiydybz = ay*dybz(:,1:nypmlbot,:) + bypml*psiydybz
      dybz(:,1:nypmlbot,:) = dybz(:,1:nypmlbot,:)/kappay + psiydybz
    elseif (IAM_YTOP) then
      RHSpsiyinductionbx = ay*dyflux3(:,dim2(rank)-nypmltop+1:dim2(rank),:) + &
      bypml*psiyinductionbx
      dyflux3(:,dim2(rank)-nypmltop+1:dim2(rank),:) = &
      dyflux3(:,dim2(rank)-nypmltop+1:dim2(rank),:)/kappay + psiyinductionbx
      RHSpsiyinductionbz = ay*dyflux1(:,dim2(rank)-nypmltop+1:dim2(rank),:) + &
      bypml*psiyinductionbz
      dyflux1(:,dim2(rank)-nypmltop+1:dim2(rank),:) = &
      dyflux1(:,dim2(rank)-nypmltop+1:dim2(rank),:)/kappay + psiyinductionbz
      RHSpsiydybx = ay*dybx(:,dim2(rank)-nypmltop+1:dim2(rank),:) + bypml*psiydybx
      dybx(:,dim2(rank)-nypmltop+1:dim2(rank),:) = &
      dybx(:,dim2(rank)-nypmltop+1:dim2(rank),:)/kappay + psiydybx
      RHSpsiydybz = ay*dybz(:,dim2(rank)-nypmltop+1:dim2(rank),:) + bypml*psiydybz
      dybz(:,dim2(rank)-nypmltop+1:dim2(rank),:) = &
      dybz(:,dim2(rank)-nypmltop+1:dim2(rank),:)/kappay + psiydybz
    endif
  endif

 RHSb_x = reduction*(dyflux3 - dzflux2)
 RHSb_y = reduction*(dzflux1 - dxflux3)
 RHSb_z = reduction*(dxflux2 - dyflux1)

 div = dvxdx + dvydy + dvzdz

 curlbx = dybz - dzby
 curlby = -dxbz + dzbx
 curlbz = dxby - dybx


 ! CONTIUNUITY 
 ! --------------------------------------

  RHScont = - gradrho0_x * v_x - gradrho0_y * v_y &
           - gradrho0_z * v_z - rho0 * div -spongexyz*rho

 ! PRESSURE

 RHSp = - c2rho0 * div - v_x * gradp0_x - v_y * gradp0_y &
        - v_z * gradp0_z - spongexyz*p

 call cross(curlbx, curlby, curlbz, box, boy, boz, flux1, flux2, flux3)

 RHSv_x = -rhoinv*(gradp_x - (flux1*reduction))-spongexyz*v_x
 RHSv_y = -rhoinv*(gradp_y - (flux2*reduction))-spongexyz*v_y
 RHSv_z = -rhoinv*(gradp_z - (flux3*reduction))-spongexyz*v_z

 call cross(curlbox, curlboy, curlboz, bx, by, bz, flux1, flux2, flux3)

 RHSv_x = RHSv_x + rhoinv*flux1*reduction
 RHSv_y = RHSv_y + rhoinv*flux2*reduction
 RHSv_z = RHSv_z + rhoinv*flux3*reduction

 if (FLOWS) then

 flux1 = - v0_z * by + v0_y * bz
 flux2 = - v0_x * bz + v0_z * bx
 flux3 = - v0_y * bx + v0_x * by

 call ddz(flux2,dzflux2,bv_mag)
 call ddz(flux1,dzflux1,bv_mag)
 call ddx(flux2,dxflux2,bh_mag)
 call ddx(flux3,dxflux3,bh_mag)
 call ddy(flux3,dyflux3,bh_mag)
 call ddy(flux1,dyflux1,bh_mag)

 RHSb_x = RHSb_x + (dyflux3 -dzflux2)
 RHSb_y = RHSb_y + (dzflux1 - dxflux3)
 RHSb_z = RHSb_z + (dxflux2 - dyflux1)

 RHSp = RHSp - gradp_x * v0_x &
 -gradp_y*v0_y - gradp_z * v0_z

 call ddz(rho,gradrho_z,bv_prerho)
 call ddx(rho,gradrho_x,bh_prerho)
 call ddy(rho,gradrho_y,bh_prerho)

 RHScont = RHScont - gradrho_x * v0_x &
  - gradrho_y*v0_y - gradrho_z * v0_z - rho * div0

  call ddz(v_x,dvxdz,bv_vel)
  call ddx(v_z, dvzdx,bh_vel)
  call ddz(v_y,dvydz,bv_vel)
  call ddx(v_y, dvydx,bh_vel)
  call ddy(v_z, dvzdy,bh_vel)
  call ddy(v_x, dvxdy,bh_vel)

  RHSv_x = RHSv_x - 2.0_dp*(v0_x*dvxdx + v0_y*dvxdy + v0_z*dvxdz + &
  v_x * dv0xdx + v_y*dv0xdy + v_z * dv0xdz)

  RHSv_z = RHSv_z - 2.0_dp*(v0_x*dvzdx + v0_y*dvzdy + v0_z*dvzdz + &
  v_x * dv0zdx + v_y*dv0zdy + v_z * dv0zdz)

  RHSv_y = RHSv_z - 2.0_dp*(v0_x*dvydx + v0_y*dvydy + v0_z*dvydz + &
  v_x * dv0ydx + v_y*dv0ydy + v_z * dv0ydz)

 endif


 do kk= 1,dim3(rank)

   RHSv_z(:,:,kk) = RHSv_z(:,:,kk) - rhoinv(:,:,kk) * rho(:,:,kk)*g(kk)

 if (pulse) then
   if (pulse_dir == 1) then
     RHSv_x(:,:,kk) = RHSv_x(:,:,kk) + forcing_p(:,:,kk)
   elseif (pulse_dir == 2) then
     RHSv_y(:,:,kk) = RHSv_y(:,:,kk) + forcing_p(:,:,kk)
   elseif (pulse_dir == 3) then
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
     if (IAM_YTOP) scr(:,dim2(rank),:,:) = 0.0_dp
     if (IAM_YBOT) scr(:,1,:,:) = 0.0_dp
   endif
   if (IAM_ZTOP) scr(:,:,dim3(rank),:) = 0.0_dp
   if (IAM_ZBOT) scr(:,:,1,:) = 0.0_dp
 endif

 if (DAMPING) then
  call damp_velocity(v_x,term)
  RHSv_x = RHSv_x - term
  call damp_velocity(v_y,term)
  RHSv_y = RHSv_y - term
  call damp_velocity(v_z,term)
  RHSv_z = RHSv_z - term
 endif

 end subroutine MP_MHD_3D
!==================================================================================

subroutine MP_MHD_DISPL_3D

 ! I can't really do the conservative form - doesn't work for the incredibly small density 
 ! and forcing functions anyway. I don't know what the error is anyway
 implicit none
 integer bh_disp,bv_disp,bh_mag,bv_mag,bh_prerho,bv_prerho
 real(dp), dimension(dim1(rank),dim2(rank),dim3(rank)) :: flux1,flux2,flux3
 real(dp), dimension(dim1(rank),dim2(rank),dim3(rank)) :: dzflux1,dzflux2,dxflux2,dxflux3,dyflux1,dyflux3
 real(dp), dimension(dim1(rank),dim2(rank),dim3(rank)) :: dxbz,dxby,dybx,dybz,dzbx,dzby
 real(dp), dimension(dim1(rank),dim2(rank),dim3(rank)) :: term
 real(dp), allocatable, dimension(:,:,:) ::temp_x,temp_y, temp_z
 integer kk

 if (periodic) then
 bh_disp = 5
 bv_disp = 1
 bh_mag = 5
 bv_mag = 1
 bh_prerho = 5
 bv_prerho = 1
 else
 bh_disp = 1
 bv_disp = 1
 bh_mag = 1
 bv_mag = 1
 bh_prerho = 1
 bv_prerho = 1
 endif

 flux1 = - (xi_z * boy - xi_y * boz)
 flux2 = - (xi_x * boz - xi_z * box)
 flux3 = - (xi_y * box - xi_x * boy)

 call ddxyz(xi_x, dxixdx, xi_y, dxiydy, xi_z, dxizdz, bh_disp,bh_disp,bv_disp)

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

    if (IAM_YBOT) then
      RHSpsiyvy = ay*dxiydy(:,1:nypmlbot,:) + bypml*psiyvy
      dxiydy(:,1:nypmlbot,:) = dxiydy(:,1:nypmlbot,:)/kappay + psiyvy
    elseif (IAM_YTOP) then
      RHSpsiyvy = ay*dxiydy(:,dim2(rank)-nypmltop+1:dim2(rank),:) + bypml*psiyvy
      dxiydy(:,dim2(rank)-nypmltop+1:dim2(rank),:) = &
dxiydy(:,dim2(rank)-nypmltop+1:dim2(rank),:)/kappay + psiyvy
    endif

  endif

 div = dxixdx + dxiydy + dxizdz

 p = - c2rho0 * div - xi_z * gradp0_z - xi_x * gradp0_x - xi_y * gradp0_y
 rho = - rho0 * div - xi_z * gradrho0_z - xi_x * gradrho0_x - xi_y * gradrho0_y

 call ddxyz(p, gradp_x, p, gradp_y, p, gradp_z, bh_prerho,bh_prerho,bv_prerho)

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

    if (IAM_YBOT) then
      RHSpsiyp = ay*gradp_y(:,1:nypmlbot,:) + bypml*psiyp
      gradp_y(:,1:nypmlbot,:) = gradp_y(:,1:nypmlbot,:)/kappay + psiyp
    elseif (IAM_YTOP) then
      RHSpsiyp = ay*gradp_y(:,dim2(rank)-nypmltop+1:dim2(rank),:) + bypml*psiyp
      gradp_y(:,dim2(rank)-nypmltop+1:dim2(rank),:) = &
gradp_y(:,dim2(rank)-nypmltop+1:dim2(rank),:)/kappay + psiyp
    endif

  endif

 call ddxyz(flux3, dxflux3, flux3, dyflux3, flux2, dzflux2,bh_mag,bh_mag,bv_mag )
 call ddxyz(flux2, dxflux2, flux1, dyflux1, flux1, dzflux1, bh_mag,bh_mag,bv_mag)


    if (USE_PML) then
    if (IAM_ZBOT) then
    RHSpsizinductionbx = az*dzflux2(:,:,1:nzpmlbot) + bzpml*psizinductionbx
    dzflux2(:,:,1:nzpmlbot) = dzflux2(:,:,1:nzpmlbot)/kappaz + psizinductionbx
    RHSpsizinductionby = az*dzflux1(:,:,1:nzpmlbot) + bzpml*psizinductionby
    dzflux1(:,:,1:nzpmlbot) = dzflux1(:,:,1:nzpmlbot)/kappaz + psizinductionby
    elseif (IAM_ZTOP) then
    RHSpsizinductionbx = az*dzflux2(:,:,dim3(rank)-nzpmltop+1:dim3(rank)) + &
bzpml*psizinductionbx
    dzflux2(:,:,dim3(rank)-nzpmltop+1:dim3(rank)) = &
dzflux2(:,:,dim3(rank)-nzpmltop+1:dim3(rank))/kappaz + psizinductionbx
    RHSpsizinductionby = az*dzflux1(:,:,dim3(rank)-nzpmltop+1:dim3(rank)) + &
bzpml*psizinductionby
    dzflux1(:,:,dim3(rank)-nzpmltop+1:dim3(rank)) = &
    dzflux1(:,:,dim3(rank)-nzpmltop+1:dim3(rank))/kappaz + psizinductionby
    endif
  endif
  if (USE_HPML) then
    if (IAM_XBOT) then
      RHSpsixinductionby = ax*dxflux3(1:nxpmlbot,:,:) + bxpml*psixinductionby
      dxflux3(1:nxpmlbot,:,:) = dxflux3(1:nxpmlbot,:,:)/kappax + psixinductionby
      RHSpsixinductionbz = ax*dxflux2(1:nxpmlbot,:,:) + bxpml*psixinductionbz
      dxflux2(1:nxpmlbot,:,:) = dxflux2(1:nxpmlbot,:,:)/kappax + psixinductionbz
    elseif (IAM_XTOP) then
      RHSpsixinductionby = ax*dxflux3(dim1(rank)-nxpmltop+1:dim1(rank),:,:) + &
      bxpml*psixinductionby
      dxflux3(dim1(rank)-nxpmltop+1:dim1(rank),:,:) = &
      dxflux3(dim1(rank)-nxpmltop+1:dim1(rank),:,:)/kappax + psixinductionby
      RHSpsixinductionbz = ax*dxflux2(dim1(rank)-nxpmltop+1:dim1(rank),:,:) + &
      bxpml*psixinductionbz
      dxflux2(dim1(rank)-nxpmltop+1:dim1(rank),:,:) = &
      dxflux2(dim1(rank)-nxpmltop+1:dim1(rank),:,:)/kappax + psixinductionbz
    endif
    if (IAM_YBOT) then
      RHSpsiyinductionbx = ay*dyflux3(:,1:nypmlbot,:) + bypml*psiyinductionbx
      dyflux3(:,1:nypmlbot,:) = dyflux3(:,1:nypmlbot,:)/kappay + psiyinductionbx
      RHSpsiyinductionbz = ay*dyflux1(:,1:nypmlbot,:) + bypml*psiyinductionbz
      dyflux1(:,1:nypmlbot,:) = dyflux1(:,1:nypmlbot,:)/kappay + psiyinductionbz
    elseif (IAM_YTOP) then
      RHSpsiyinductionbx = ay*dyflux3(:,dim2(rank)-nypmltop+1:dim2(rank),:) + &
      bypml*psiyinductionbx
      dyflux3(:,dim2(rank)-nypmltop+1:dim2(rank),:) = &
      dyflux3(:,dim2(rank)-nypmltop+1:dim2(rank),:)/kappay + psiyinductionbx
      RHSpsiyinductionbz = ay*dyflux1(:,dim2(rank)-nypmltop+1:dim2(rank),:) + &
      bypml*psiyinductionbz
      dyflux1(:,dim2(rank)-nypmltop+1:dim2(rank),:) = &
      dyflux1(:,dim2(rank)-nypmltop+1:dim2(rank),:)/kappay + psiyinductionbz
    endif
  endif

 
 bx =   dyflux3 - dzflux2
 by = - dxflux3 + dzflux1
 bz =   dxflux2 - dyflux1

 call ddxyz(by, dxby, bz, dybz,  bx, dzbx, bh_mag,bh_mag,bv_mag)
 call ddxyz(bz, dxbz, bx, dybx, by, dzby, bh_mag,bh_mag,bv_mag)

 ! NOTE THAT FLUX1 IS OVERWRITTEN IN THE FINAL CALL

    if (USE_PML) then
    if (IAM_ZBOT) then
    RHSpsizdzbx = az*dzbx(:,:,1:nzpmlbot) + bzpml*psizdzbx
    dzbx(:,:,1:nzpmlbot) = dzbx(:,:,1:nzpmlbot)/kappaz + psizdzbx
    RHSpsizdzby = az*dzby(:,:,1:nzpmlbot) + bzpml*psizdzby
    dzby(:,:,1:nzpmlbot) = dzby(:,:,1:nzpmlbot)/kappaz + psizdzby
    elseif (IAM_ZTOP) then
    RHSpsizdzbx = az*dzbx(:,:,dim3(rank)-nzpmltop+1:dim3(rank)) + bzpml*psizdzbx
    dzbx(:,:,dim3(rank)-nzpmltop+1:dim3(rank)) = &
    dzbx(:,:,dim3(rank)-nzpmltop+1:dim3(rank))/kappaz + psizdzbx
    RHSpsizdzby = az*dzby(:,:,dim3(rank)-nzpmltop+1:dim3(rank)) + bzpml*psizdzby
    dzby(:,:,dim3(rank)-nzpmltop+1:dim3(rank)) = &
    dzby(:,:,dim3(rank)-nzpmltop+1:dim3(rank))/kappaz + psizdzby
    endif
  endif
  if (USE_HPML) then
    if (IAM_XBOT) then
      RHSpsixdxby = ax*dxby(1:nxpmlbot,:,:) + bxpml*psixdxby
      dxby(1:nxpmlbot,:,:) = dxby(1:nxpmlbot,:,:)/kappax + psixdxby
      RHSpsixdxbz = ax*dxbz(1:nxpmlbot,:,:) + bxpml*psixdxbz
      dxbz(1:nxpmlbot,:,:) = dxbz(1:nxpmlbot,:,:)/kappax + psixdxbz
    elseif (IAM_XTOP) then
      RHSpsixdxby = ax*dxby(dim1(rank)-nxpmltop+1:dim1(rank),:,:) + &
bxpml*psixdxby
      dxby(dim1(rank)-nxpmltop+1:dim1(rank),:,:) = &
      dxby(dim1(rank)-nxpmltop+1:dim1(rank),:,:)/kappax + psixdxby
      RHSpsixdxbz = ax*dxbz(dim1(rank)-nxpmltop+1:dim1(rank),:,:) + &
bxpml*psixdxbz
      dxbz(dim1(rank)-nxpmltop+1:dim1(rank),:,:) = &
      dxbz(dim1(rank)-nxpmltop+1:dim1(rank),:,:)/kappax + psixdxbz
    endif
    if (IAM_YBOT) then
      RHSpsiydybx = ay*dybx(:,1:nypmlbot,:) + bypml*psiydybx
      dybx(:,1:nypmlbot,:) = dybx(:,1:nypmlbot,:)/kappay + psiydybx
      RHSpsiydybz = ay*dybz(:,1:nypmlbot,:) + bypml*psiydybz
      dybz(:,1:nypmlbot,:) = dybz(:,1:nypmlbot,:)/kappay + psiydybz
    elseif (IAM_YTOP) then
      RHSpsiydybx = ay*dybx(:,dim2(rank)-nypmltop+1:dim2(rank),:) + &
bypml*psiydybx
      dybx(:,dim2(rank)-nypmltop+1:dim2(rank),:) = &
      dybx(:,dim2(rank)-nypmltop+1:dim2(rank),:)/kappay + psiydybx
      RHSpsiydybz = ay*dybz(:,dim2(rank)-nypmltop+1:dim2(rank),:) + &
bypml*psiydybz
      dybz(:,dim2(rank)-nypmltop+1:dim2(rank),:) = &
      dybz(:,dim2(rank)-nypmltop+1:dim2(rank),:)/kappay + psiydybz
    endif
  endif

 curlbx = dybz - dzby
 curlby = -dxbz + dzbx
 curlbz = dxby - dybx

 call cross(curlbx, curlby, curlbz, box, boy, boz, flux1, flux2, flux3)

 RHSv_x = rhoinv * (- gradp_x + flux1)
 RHSv_y = rhoinv * (- gradp_y + flux2)
 RHSv_z = rhoinv * (- gradp_z + flux3)

 call cross(curlbox, curlboy, curlboz, bx, by, bz, flux1, flux2, flux3)

 RHSv_x = RHSv_x + rhoinv*flux1 - spongexyz*v_x 
 RHSv_y = RHSv_y + rhoinv*flux2 - spongexyz*v_y
 RHSv_z = RHSv_z + rhoinv*flux3 - spongexyz*v_z

if (FLOWS) then

  allocate(temp_x(dim1(rank),dim2(rank),dim3(rank)),temp_y(dim1(rank),dim2(rank),dim3(rank)),&
   temp_z(dim1(rank),dim2(rank),dim3(rank)))

   call ddxyz(v_x, temp_x, v_x, temp_y, v_x, temp_z, bh_disp,bh_disp,bv_disp)

   RHSv_x = RHSv_x- 2.0_dp * (v0_x * temp_x + v0_y * temp_y + v0_z * temp_z)

   call ddxyz(v_y, temp_x, v_y, temp_y, v_y, temp_z, bh_disp,bh_disp,bv_disp)

   RHSv_y = RHSv_y- 2.0_dp * (v0_x * temp_x + v0_y * temp_y + v0_z * temp_z)

   call ddxyz(v_z, temp_x, v_z, temp_y, v_z, temp_z, bh_disp,bh_disp,bv_disp)

   RHSv_z = RHSv_z- 2.0_dp * (v0_x * temp_x + v0_y * temp_y + v0_z * temp_z)

   deallocate(temp_x, temp_y, temp_z)
endif


 do kk= 1,dim3(rank)
 RHSv_z(:,:,kk) = RHSv_z(:,:,kk) - rhoinv(:,:,kk) * rho(:,:,kk)*g(kk)
 enddo 

 RHSxi_x = v_x
 RHSxi_y = v_y
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
     RHSv_z(:,:,kk) = RHSv_z(:,:,kk) + forcing*source_dep(kk)
   enddo
 endif
 if (IAM_NOTB) then
 else
   if (.not. PERIODIC) then
     if (IAM_XTOP) scr(dim1(rank),:,:,:) = 0.0_dp
     if (IAM_XBOT) scr(1,:,:,:) = 0.0_dp
     if (IAM_YTOP) scr(:,dim2(rank),:,:) = 0.0_dp
     if (IAM_YBOT) scr(:,1,:,:) = 0.0_dp
   endif
   if (IAM_ZTOP) scr(:,:,dim3(rank),:) = 0.0_dp
   if (IAM_ZBOT) scr(:,:,1,:) = 0.0_dp
 endif

  if (DAMPING) then
  call damp_velocity(v_x,term)
  RHSv_x = RHSv_x - term
  call damp_velocity(v_y,term)
  RHSv_y = RHSv_y - term
  call damp_velocity(v_z,term)
  RHSv_z = RHSv_z - term
 endif

 end subroutine MP_MHD_DISPL_3D

 
end module PHYSICS_RHS
