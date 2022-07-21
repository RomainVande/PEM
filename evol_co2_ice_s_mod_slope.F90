!
! $Id $
!
SUBROUTINE evol_co2_ice_s_slope(qsurf,tendencies_co2_ice_phys,&
                             iim_input,jjm_input,ngrid,cell_area,STOPPING,nslope)

  USE temps_mod_evol, ONLY: dt_pem

      IMPLICIT NONE

!=======================================================================
!
!  Routine that compute the evolution of the water ice
!
!=======================================================================

!   arguments:
!   ----------

!   INPUT

  INTEGER, intent(in) :: iim_input,jjm_input, ngrid,nslope   ! # of grid points along longitude/latitude/ total
!  REAL, intent(in) ::  tendencies_h2o_ice_phys(ngrid) ! physical point field : Evolution of perenial ice over one year
  REAL, intent(in) ::  cell_area(ngrid)

!   OUTPUT
  REAL, INTENT(INOUT) ::  qsurf(ngrid,nslope)                ! physical point field : Previous and actual density of water ice
  LOGICAL :: STOPPING
  REAL, intent(inout) ::  tendencies_co2_ice_phys(ngrid,nslope) ! physical point field : Evolution of perenial ice over one year


!   local:
!   ----

  INTEGER :: i,j,ig0,islope                                  ! loop variable
!  REAL :: not_budget, budget
  REAL :: pos_tend, neg_tend, real_coefficient,negative_part
  REAL ::  new_tendencies(ngrid)

  STOPPING=.false.


! Evolution of the water ice for each physical point
  do i=1,ngrid
    do islope=1,nslope
      qsurf(i,islope)=qsurf(i,islope)+tendencies_co2_ice_phys(i,islope)*dt_pem
      if (qsurf(i,islope).lt.0) then
        qsurf(i,islope)=0.
        tendencies_co2_ice_phys(i,islope)=0.
      endif
    enddo
  enddo


END SUBROUTINE evol_co2_ice_s_slope
