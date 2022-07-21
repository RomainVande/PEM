!
! $Id $
!
SUBROUTINE evol_co2_ice_s(qsurf,tendencies_co2_ice_phys,&
                             iim_input,jjm_input,ngrid,cell_area,STOPPING)

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

  INTEGER, intent(in) :: iim_input,jjm_input, ngrid   ! # of grid points along longitude/latitude/ total
!  REAL, intent(in) ::  tendencies_h2o_ice_phys(ngrid) ! physical point field : Evolution of perenial ice over one year
  REAL, intent(in) ::  cell_area(ngrid)

!   OUTPUT
  REAL, INTENT(INOUT) ::  qsurf(ngrid)                ! physical point field : Previous and actual density of water ice
  LOGICAL :: STOPPING
  REAL, intent(inout) ::  tendencies_co2_ice_phys(ngrid) ! physical point field : Evolution of perenial ice over one year


!   local:
!   ----

  INTEGER :: i,j,ig0                                  ! loop variable
!  REAL :: not_budget, budget
  REAL :: pos_tend, neg_tend, real_coefficient,negative_part
  REAL ::  new_tendencies(ngrid)

  STOPPING=.false.


! Evolution of the water ice for each physical point
  do i=1,ngrid
      qsurf(i)=qsurf(i)+tendencies_co2_ice_phys(i)*dt_pem
      if (qsurf(i).lt.0) then
        qsurf(i)=0.
        tendencies_co2_ice_phys(i)=0.
      endif
  enddo


END SUBROUTINE evol_co2_ice_s
