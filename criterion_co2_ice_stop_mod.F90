!
! $Id $
!
SUBROUTINE criterion_co2_ice_stop(cell_area,initial_co2_ice,co2ice,STOPPING,ngrid,latitude,n_band_lat)

  USE temps_mod_evol, ONLY: alpha_criterion

      IMPLICIT NONE

!=======================================================================
!
!  Routine that checks if the criterion to stop the PEM is reached
!
!=======================================================================

!   arguments:
!   ----------

!   INPUT
  INTEGER, intent(in) :: ngrid                  ! # of grid physical grid points 
  REAL,    intent(in) :: cell_area(ngrid)       ! physical point field : Area of the cells
  REAL,    intent(in) ::  co2ice(ngrid)          ! physical point field : Actual density of water ice
  REAL,    intent(in) ::  latitude(ngrid)          ! physical point field : Latitude
  REAL,    intent(in) ::  initial_co2_ice(n_band_lat)  ! Initial/Actual surface of water ice



!   OUTPUT
  LOGICAL, intent(out) :: STOPPING              ! Logical : is the criterion reached?

!   local:
!   -----
  INTEGER :: i,j,n_band_lat                    ! Loop
  REAL :: present_co2(n_band_lat)  ! Initial/Actual surface of water ice
  REAL :: pi

!=======================================================================

      pi=4.D0*DATAN(1.D0)

!   initialisation to false
    STOPPING=.FALSE.

     do j=1,n_band_lat
        present_co2(j)=0.
     enddo

  do i=1,ngrid
            j=floor((latitude(i)+(pi/2))/(pi)*n_band_lat)+1
      if(j.GT.n_band_lat) then
          j=n_band_lat
      endif
      present_co2(j)=present_co2(j)+co2ice(i)*cell_area(i)
  enddo
  
!   check of the criterion
  do j=1,n_band_lat
    if(present_co2(j).LT.initial_co2_ice(j)*(1-alpha_criterion) .OR. &
       present_co2(j).GT.initial_co2_ice(j)*(1+alpha_criterion)) then
         STOPPING=.TRUE. 
         print *, "j", j
         print *, "present_co2(j)", present_co2(j)
         print *, "initial_co2_ice(j)", initial_co2_ice(j)
    endif
  enddo

END SUBROUTINE criterion_co2_ice_stop





