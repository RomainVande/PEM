!
! $Id $
!
SUBROUTINE recomp_tend_co2_slope(tendencies_co2_ice_phys,tendencies_co2_ice_phys_ini,vmr_co2_gcm,vmr_co2_pem,ps_GCM_2,global_ave_press_GCM,global_ave_press_new,timelen,ngrid,nslope)

      IMPLICIT NONE

!=======================================================================
!
!  Routine that compute the evolution of the tendencie for co2 ice
!
!=======================================================================

!   arguments:
!   ----------

!   INPUT
  INTEGER, intent(in) :: timelen,ngrid,nslope
  REAL, INTENT(in) ::  vmr_co2_gcm(ngrid,timelen)                ! physical point field : Volume mixing ratio of co2 in the first layer
  REAL, INTENT(in) ::  vmr_co2_pem(ngrid,timelen)                ! physical point field : Volume mixing ratio of co2 in the first layer
  REAL, intent(in) :: ps_GCM_2(ngrid,timelen)                 ! physical point field : Surface pressure in the GCM
!  REAL, INTENT(in) ::  q_co2_GCM(ngrid)                ! physical point field : Density of co2 in the first layer
!  REAL, intent(in) :: ps_GCM(ngrid)                 ! physical point field : Density of co2 in the first layer
  REAL, intent(in) :: global_ave_press_GCM
  REAL, intent(in) :: global_ave_press_new
  REAL, intent(in) ::  tendencies_co2_ice_phys_ini(ngrid,nslope) ! physical point field : Evolution of perenial ice over one year


!   OUTPUT
  REAL, intent(inout) ::  tendencies_co2_ice_phys(ngrid,nslope) ! physical point field : Evolution of perenial ice over one year


!   local:
!   ----

  INTEGER :: i,t,islope
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
    do islope=1,nslope
       ave=0.
    do t=1,timelen

!       write(*,*)'i,t=',i,t,islope, alpha,beta,ave,vmr_co2_gcm(i,t),vmr_co2_pem(i,t),ps_GCM_2(i,t),global_ave_press_GCM,global_ave_press_new
       ave=ave+(beta/(alpha-log(vmr_co2_gcm(i,t)*ps_GCM_2(i,t)/100)))**4  &
              -(beta/(alpha-log(vmr_co2_pem(i,t)*ps_GCM_2(i,t)*global_ave_press_GCM/global_ave_press_new/100)))**4
    enddo
!    print *, "i", i
!    print *, "tendencies_co2_ice_phys_ini bef", tendencies_co2_ice_phys_ini(i,islope)
!    tendencies_co2_ice_phys(i,islope)=tendencies_co2_ice_phys_ini(i,islope)+coef*ave/timelen
    tendencies_co2_ice_phys(i,islope)=tendencies_co2_ice_phys_ini(i,islope)-coef*ave/timelen

!    print *, "tendencies after", tendencies_co2_ice_phys(i,islope)
!    print *, "ave", ave
!    print *, "timelen", timelen
!    print *, "vmr_co2_pem(i,t)*ps_GCM_2(i,t)", vmr_co2_pem(i,t)*ps_GCM_2(i,t)
!    print *, "-------------------"
  enddo
  enddo


END SUBROUTINE recomp_tend_co2_slope
