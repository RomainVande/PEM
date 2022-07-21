!
! $Id $
!
SUBROUTINE evol_h2o_ice_s_slope(qsurf,tendencies_h2o_ice_phys,&
                             iim_input,jjm_input,ngrid,cell_area,STOPPING,nslope)

  USE temps_mod_evol, ONLY: dt_pem
  	  use comslope_mod, ONLY: subslope_dist

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
  REAL, intent(inout) ::  tendencies_h2o_ice_phys(ngrid,nslope) ! physical point field : Evolution of perenial ice over one year


!   local:
!   ----

  INTEGER :: i,j,ig0,islope                                  ! loop variable
!  REAL :: not_budget, budget
  REAL :: pos_tend, neg_tend, real_coefficient,negative_part
  REAL ::  new_tendencies(ngrid,nslope)


!=======================================================================

!  budget=sum(qsurf(:))

  pos_tend=0.
  neg_tend=0.

  do i=1,ngrid
     do islope=1,nslope
     if (qsurf(i,islope).GT.0) then
         if (tendencies_h2o_ice_phys(i,islope).GT.0) then
            pos_tend=pos_tend+tendencies_h2o_ice_phys(i,islope)*cell_area(i)*subslope_dist(i,islope)
         else
            neg_tend=neg_tend-tendencies_h2o_ice_phys(i,islope)*cell_area(i)*subslope_dist(i,islope)
         endif
     endif
     enddo
  enddo

  print *, "pos_tend", pos_tend
  print *, "neg_tend", neg_tend

  if(neg_tend.GT.pos_tend .and. pos_tend.GT.0) then
     do i=1,ngrid
       do islope=1,nslope
       if(tendencies_h2o_ice_phys(i,islope).LT.0) then
          print *, "pos_tend/neg_tend", pos_tend/neg_tend
          new_tendencies(i,islope)=tendencies_h2o_ice_phys(i,islope)*(pos_tend/neg_tend)
       else
          new_tendencies(i,islope)=tendencies_h2o_ice_phys(i,islope)
       endif
       enddo
     enddo
  elseif(neg_tend.LT.pos_tend .and. neg_tend.GT.0) then
          print *, "neg_tend/pos_tend", neg_tend/pos_tend
     do i=1,ngrid
       do islope=1,nslope
       if(tendencies_h2o_ice_phys(i,islope).LT.0) then
          new_tendencies(i,islope)=tendencies_h2o_ice_phys(i,islope)
       else
          new_tendencies(i,islope)=tendencies_h2o_ice_phys(i,islope)*(neg_tend/pos_tend)
       endif
       enddo
     enddo
  elseif(pos_tend.EQ.0 .OR. neg_tend.EQ.0) then 

!      call criterion_ice_stop(cell_area,1,qsurf*0.,STOPPING,ngrid,cell_area)
      call criterion_ice_stop_water_slope(cell_area,1,qsurf(:,:)*0.,STOPPING,ngrid,cell_area)
      do i=1,ngrid
         do islope=1,nslope
          new_tendencies(i,islope)=0
         enddo
      enddo
  endif

 
  



! Evolution of the water ice for each physical point
  do i=1,ngrid
    do islope=1, nslope
!      qsurf(i)=qsurf(i)+tendencies_h2o_ice_phys(i)*dt_pem
      qsurf(i,islope)=qsurf(i,islope)+new_tendencies(i,islope)*dt_pem
!      budget=budget+tendencies_h2o_ice_phys(i)*dt_pem
      if (qsurf(i,islope).lt.0) then
!        not_budget=not_budget+qsurf(i)
        print *, "NNqsurf(i,islope)", qsurf(i,islope)
        print *, "NNnew_tendencies(i,islope)", new_tendencies(i,islope)
        print *, "NNtendencies_h2o_ice_phys(i,islope)", tendencies_h2o_ice_phys(i,islope)
        negative_part=negative_part-qsurf(i,islope)*cell_area(i)*subslope_dist(i,islope)
        qsurf(i,islope)=0.
        tendencies_h2o_ice_phys(i,islope)=0.
        print *, "NNineg", i
      endif
      if(qsurf(i,islope).NE.qsurf(i,islope)) then
          print *, "qsurf(i,islope)",qsurf(i,islope)
          print *, "new_tendencies",new_tendencies(i,islope)
          print *, "tendencies_h2o_ice_phys",tendencies_h2o_ice_phys(i,islope)
          print *, "i", i
          print *,"islope",islope
      endif
    enddo
  enddo

  print *, "negative_part", negative_part
  real_coefficient=negative_part/pos_tend
  print *, "real_coefficient", real_coefficient

  do i=1,ngrid
    do islope=1, nslope
     if(new_tendencies(i,islope).GT.0) then
         qsurf(i,islope)=qsurf(i,islope)-new_tendencies(i,islope)*real_coefficient*dt_pem
     endif
    enddo
  enddo



! Conservation of water ice
!  qsurf(:)=qsurf(:)*budget/sum(qsurf(:))


END SUBROUTINE evol_h2o_ice_s_slope
