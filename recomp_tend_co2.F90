!
! $Id $
!
SUBROUTINE recomp_tend_co2(tendencies_co2_ice_phys,vmr_co2_gcm,ps_GCM_2,global_ave_press_GCM,global_ave_press_new,timelen,ngrid)

      IMPLICIT NONE

!=======================================================================
!
!  Routine that compute the evolution of the tendencie for co2 ice
!
!=======================================================================

!   arguments:
!   ----------

!   INPUT
  INTEGER, intent(in) :: timelen,ngrid
  REAL, INTENT(in) ::  vmr_co2_gcm(ngrid,timelen)                ! physical point field : Volume mixing ratio of co2 in the first layer
  REAL, intent(in) :: ps_GCM_2(ngrid,timelen)                 ! physical point field : Surface pressure in the GCM
!  REAL, INTENT(in) ::  q_co2_GCM(ngrid)                ! physical point field : Density of co2 in the first layer
!  REAL, intent(in) :: ps_GCM(ngrid)                 ! physical point field : Density of co2 in the first layer
  REAL, intent(in) :: global_ave_press_GCM
  REAL, intent(in) :: global_ave_press_new


!   OUTPUT
  REAL, intent(inout) ::  tendencies_co2_ice_phys(ngrid) ! physical point field : Evolution of perenial ice over one year


!   local:
!   ----

  INTEGER :: i,t
  REAL :: eps, sigma, L, beta, alpha, coef, ave

  eps=0.95
  sigma=5.678E-8
  L=5.71*10**5
  beta=3182.48
  alpha=23.3494

  coef=669*24*3600*eps*sigma/L

  print *, "coef", coef
  print *, "global_ave_press_GCM", global_ave_press_GCM
  print *, "global_ave_press_new", global_ave_press_new

! Evolution of the water ice for each physical point
  do i=1,ngrid
       ave=0.
       if (tendencies_co2_ice_phys(i).LT. 1E-5 .and. tendencies_co2_ice_phys(i).GT.-1E-5) then
       else
    do t=1,timelen
       ave=ave+(beta/(alpha-log(vmr_co2_gcm(i,t)*ps_GCM_2(i,t)/100)))**4  &
              -(beta/(alpha-log(vmr_co2_gcm(i,t)*ps_GCM_2(i,t)*global_ave_press_GCM/global_ave_press_new/100)))**4
    enddo
    print *, "i", i
    print *, "tendencies bef", tendencies_co2_ice_phys(i)
    tendencies_co2_ice_phys(i)=tendencies_co2_ice_phys(i)+coef*ave/timelen
       endif
    print *, "tendencies after", tendencies_co2_ice_phys(i)
    print *, "ave", ave
    print *, "timelen", timelen
    print *, "vmr_co2_gcm(i,t)*ps_GCM_2(i,t)", vmr_co2_gcm(i,t)*ps_GCM_2(i,t)
    print *, "-------------------"
  enddo


END SUBROUTINE recomp_tend_co2
