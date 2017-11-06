Module Time_Step

! --------------------------------------------------------------------------
! MPI Version of the Spherical Acoustic Sun Simulator.
! Copyright 2006, Shravan Hanasoge
                                                                                                                                                         
! Hansen Experimental Physics Laboratory
! 455 Via Palou way, Stanford
! CA 94305, USA
! Email: shravan@stanford.edu
! --------------------------------------------------------------------------

  use all_modules
  use PHYSICS_RHS
  use PHYSICS_RHS_2D
  Implicit None

Contains
!=====================================================================================
  SUBROUTINE STEP

   implicit none

  ! All new and improved optimized RK 5, 4th order accurate integrator. 50% increase in time step.

   do step_rk = 1,5
     if (step_rk == 1) then
       temp_step = a
       if (USE_PML .AND. (IAM_ZTOP .OR. IAM_ZBOT)) psizpml = zpmlvars
       if (USE_HPML .AND. (IAM_XTOP .OR. IAM_XBOT)) psixpml = xpmlvars
       if (USE_HPML .AND. (IAM_YTOP .OR. IAM_YBOT))  psiypml = ypmlvars
     else
       temp_step = a + betas(step_rk)*scr
       if (USE_PML .AND. (IAM_ZTOP .OR. IAM_ZBOT))   psizpml = zpmlvars + betas(step_rk) * scrzpml
       if (USE_HPML .AND. (IAM_XTOP .OR. IAM_XBOT))  psixpml = xpmlvars + betas(step_rk)*scrxpml
       if (USE_HPML .AND. (IAM_YTOP .OR. IAM_YBOT))  psiypml = ypmlvars + betas(step_rk)*scrypml
     endif
 
     if (option==1) then
       call MP_QUIET_3D()
     else if (option == 2) then
       call MP_QUIET_DISPL_3D()
     else if (option == 3) then
       call MP_MHD_3D()
     else if (option == 4) then
       call MP_MHD_DISPL_3D()
     else if (option == 5) then
       call MP_QUIET_2D()
     else if (option == 6) then
       call MP_QUIET_DISPL_2D()
     else if (option == 7) then
       call MP_MHD_2D()
     else if (option == 8) then
       call MP_MHD_DISPL_2D
     endif  
   enddo
    
   a = a + deltat*scr
   if (USE_PML .AND. (IAM_ZTOP .OR. IAM_ZBOT)) zpmlvars = zpmlvars + deltat * scrzpml
   if (USE_HPML .AND. (IAM_XTOP .OR. IAM_XBOT)) xpmlvars = xpmlvars + deltat*scrxpml
   if (USE_HPML .AND. (IAM_YTOP .OR. IAM_YBOT)) ypmlvars = ypmlvars + deltat*scrypml

  End SUBROUTINE step
!=====================================================================================
End Module Time_Step
