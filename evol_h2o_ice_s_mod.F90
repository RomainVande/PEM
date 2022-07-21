!
! $Id $
!
SUBROUTINE evol_h2o_ice_s(qsurf,tendencies_h2o_ice_phys,&
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
  REAL, intent(inout) ::  tendencies_h2o_ice_phys(ngrid) ! physical point field : Evolution of perenial ice over one year


!   local:
!   ----

  INTEGER :: i,j,ig0                                  ! loop variable
!  REAL :: not_budget, budget
  REAL :: pos_tend, neg_tend, real_coefficient,negative_part
  REAL ::  new_tendencies(ngrid)


!=======================================================================

!  budget=sum(qsurf(:))

  pos_tend=0.
  neg_tend=0.

  do i=1,ngrid
     if (qsurf(i).GT.0) then
         if (tendencies_h2o_ice_phys(i).GT.0) then
            pos_tend=pos_tend+tendencies_h2o_ice_phys(i)*cell_area(i)
         else
            neg_tend=neg_tend-tendencies_h2o_ice_phys(i)*cell_area(i)
         endif
     endif
  enddo

  print *, "pos_tend", pos_tend
  print *, "neg_tend", neg_tend

  if(neg_tend.GT.pos_tend .and. pos_tend.GT.0) then
     do i=1,ngrid
       if(tendencies_h2o_ice_phys(i).LT.0) then
          print *, "pos_tend/neg_tend", pos_tend/neg_tend
          new_tendencies(i)=tendencies_h2o_ice_phys(i)*(pos_tend/neg_tend)
       else
          new_tendencies(i)=tendencies_h2o_ice_phys(i)
       endif
     enddo
  elseif(neg_tend.LT.pos_tend .and. neg_tend.GT.0) then
          print *, "neg_tend/pos_tend", neg_tend/pos_tend
     do i=1,ngrid
       if(tendencies_h2o_ice_phys(i).LT.0) then
          new_tendencies(i)=tendencies_h2o_ice_phys(i)
       else
          new_tendencies(i)=tendencies_h2o_ice_phys(i)*(neg_tend/pos_tend)
       endif
     enddo
  elseif(pos_tend.EQ.0 .OR. neg_tend.EQ.0) then 
      call criterion_ice_stop(cell_area,1,qsurf*0.,STOPPING,ngrid,cell_area)
      do i=1,ngrid
          new_tendencies(i)=0
      enddo
  endif

  print *, qsurf(1), qsurf(10), qsurf(100)

  



! Evolution of the water ice for each physical point
  do i=1,ngrid
!      qsurf(i)=qsurf(i)+tendencies_h2o_ice_phys(i)*dt_pem
      qsurf(i)=qsurf(i)+new_tendencies(i)*dt_pem
!      budget=budget+tendencies_h2o_ice_phys(i)*dt_pem
      if (qsurf(i).lt.0) then
!        not_budget=not_budget+qsurf(i)
        print *, "NNqsurf(i)", qsurf(i)
        print *, "NNnew_tendencies(i)", new_tendencies(i)
        print *, "NNtendencies_h2o_ice_phys(i)", tendencies_h2o_ice_phys(i)
        negative_part=negative_part-qsurf(i)*cell_area(i)
        qsurf(i)=0.
        tendencies_h2o_ice_phys(i)=0.
        print *, "NNineg", i
      endif
      if(qsurf(i).NE.qsurf(i)) then
          print *, "qsurf(i)",qsurf(i)
          print *, "new_tendencies",new_tendencies(i)
          print *, "tendencies_h2o_ice_phys",tendencies_h2o_ice_phys(i)
          print *, "i", i
      endif
  enddo

  print *, "negative_part", negative_part
  real_coefficient=negative_part/pos_tend
  print *, "real_coefficient", real_coefficient

  do i=1,ngrid
     if(new_tendencies(i).GT.0) then
         qsurf(i)=qsurf(i)-new_tendencies(i)*real_coefficient*dt_pem
     endif
  enddo



! Conservation of water ice
!  qsurf(:)=qsurf(:)*budget/sum(qsurf(:))


END SUBROUTINE evol_h2o_ice_s
