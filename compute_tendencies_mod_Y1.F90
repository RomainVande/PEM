!
! $Id $
!
SUBROUTINE compute_tendencies_Y1(tendencies_h2o_ice,min_h2o_ice_Y1,&
     h2o_ice_first_last_day,iim_input,jjm_input,ngrid,tendencies_h2o_ice_phys)

      IMPLICIT NONE


!=======================================================================
!
!  Compute the tendencies of the evolution of water ice over the years
!
!=======================================================================

!   arguments:
!   ----------

!   INPUT

     INTEGER, intent(in) :: iim_input,jjm_input,ngrid                             ! # of grid points along longitude/latitude/ total
     REAL, intent(in) , dimension(iim_input+1,jjm_input+1):: min_h2o_ice_Y1       ! LON x LAT field : minimum of water ice at each point for the first year
     REAL, intent(in) , dimension(iim_input+1,jjm_input+1,2):: h2o_ice_first_last_day       ! LON x LAT field : minimum of water ice at each point for the second year

!   OUTPUT
     REAL, intent(out) , dimension(iim_input+1,jjm_input+1) :: tendencies_h2o_ice ! LON x LAT field : difference between the minima = evolution of perenial ice
     REAL, intent(out) , dimension(ngrid)   :: tendencies_h2o_ice_phys            ! physical point field : difference between the minima = evolution of perenial ice

!   local:
!   ------

     INTEGER :: i,j,ig0                                                           ! loop variable

!=======================================================================


!  We compute the difference
  tendencies_h2o_ice(:,:)=h2o_ice_first_last_day(:,:,2)-h2o_ice_first_last_day(:,:,1)

!  If the difference is too small; there is no evolution
  DO j=1,jjm_input+1
    DO i = 1, iim_input
       if(abs(tendencies_h2o_ice(i,j)).LT.1.0E-10) then
          tendencies_h2o_ice(i,j)=0.
       endif
       !if(min_h2o_ice_Y1(i,j).LE.1.0E-10 .and. tendencies_h2o_ice(i,j).LT.0) then
       if(min_h2o_ice_Y1(i,j).LE.1.0E-10) then
         tendencies_h2o_ice(i,j)=0.
       endif
    ENDDO
  ENDDO


!  We reorganise the difference on the physical grid
  tendencies_h2o_ice_phys(1)=tendencies_h2o_ice(1,1)

  ig0 = 2
  DO j=2,jjm_input
    DO i = 1, iim_input
       tendencies_h2o_ice_phys(ig0)  =tendencies_h2o_ice(i,j)
       ig0= ig0 + 1
    ENDDO
  ENDDO

  tendencies_h2o_ice_phys(ig0) = tendencies_h2o_ice(1,jjm_input+1)


END SUBROUTINE compute_tendencies_Y1





