!
! $Id $
!
SUBROUTINE compute_tendencies_slope(tendencies_h2o_ice,min_h2o_ice_Y1,&
     min_h2o_ice_Y2,iim_input,jjm_input,ngrid,tendencies_h2o_ice_phys,nslope)

      IMPLICIT NONE


!=======================================================================
!
!  Compute the tendencies of the evolution of water ice over the years
!
!=======================================================================

!   arguments:
!   ----------

!   INPUT

     INTEGER, intent(in) :: iim_input,jjm_input,ngrid  ,nslope                           ! # of grid points along longitude/latitude/ total
     REAL, intent(in) , dimension(iim_input+1,jjm_input+1,nslope):: min_h2o_ice_Y1       ! LON x LAT field : minimum of water ice at each point for the first year
     REAL, intent(in) , dimension(iim_input+1,jjm_input+1,nslope):: min_h2o_ice_Y2       ! LON x LAT field : minimum of water ice at each point for the second year

!   OUTPUT
     REAL, intent(out) , dimension(iim_input+1,jjm_input+1,nslope) :: tendencies_h2o_ice ! LON x LAT field : difference between the minima = evolution of perenial ice
     REAL, intent(out) , dimension(ngrid,nslope)   :: tendencies_h2o_ice_phys            ! physical point field : difference between the minima = evolution of perenial ice

!   local:
!   ------

     INTEGER :: i,j,ig0,islope                                                           ! loop variable

!=======================================================================


!  We compute the difference
!  tendencies_h2o_ice(:,:,:)=min_h2o_ice_Y2(:,:,:)-min_h2o_ice_Y1(:,:,:)

  DO j=1,jjm_input+1
    DO i = 1, iim_input
       DO islope = 1, nslope
         tendencies_h2o_ice(i,j,islope)=min_h2o_ice_Y2(i,j,islope)-min_h2o_ice_Y1(i,j,islope)
       enddo
    ENDDO
  ENDDO

     print *, "jjm_input+1", jjm_input+1
     print *, "iim_input+1", iim_input+1
     print *, "nslope+1", nslope+1

!  If the difference is too small; there is no evolution
  DO j=1,jjm_input+1
    DO i = 1, iim_input
       DO islope = 1, nslope
         print *, "tendencies_h2o_ice(i,j,islope)LAAA", tendencies_h2o_ice(i,j,islope)
         if(abs(tendencies_h2o_ice(i,j,islope)).LT.1.0E-10) then
            tendencies_h2o_ice(i,j,islope)=0.
         endif
         print *, "tendencies_h2o_ice(i,j,islope)HERE", tendencies_h2o_ice(i,j,islope)
       enddo
    ENDDO
  ENDDO


!  We reorganise the difference on the physical grid
  tendencies_h2o_ice_phys(1,:)=tendencies_h2o_ice(1,1,:)

  ig0 = 2
  DO j=2,jjm_input
    DO i = 1, iim_input
       tendencies_h2o_ice_phys(ig0,:)  =tendencies_h2o_ice(i,j,:)
       ig0= ig0 + 1
    ENDDO
  ENDDO

  tendencies_h2o_ice_phys(ig0,:) = tendencies_h2o_ice(1,jjm_input+1,:)

  print *, "tendencies_h2o_ice_physze", tendencies_h2o_ice_phys(:,:)


END SUBROUTINE compute_tendencies_slope





