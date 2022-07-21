!
! $Id $
!
SUBROUTINE criterion_ice_stop_water(cell_area,ini_surf,qsurf,STOPPING,ngrid,initial_h2o_ice)

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
  REAL,    intent(in) ::  qsurf(ngrid)          ! physical point field : Actual density of water ice
  REAL,    intent(in) :: ini_surf
  REAL,    intent(in) :: initial_h2o_ice(ngrid)


!   OUTPUT
  LOGICAL, intent(out) :: STOPPING              ! Logical : is the criterion reached?

!   local:
!   -----
  INTEGER :: i                    ! Loop
  REAL :: present_surf  ! Initial/Actual surface of water ice

!=======================================================================

!   initialisation to false
    STOPPING=.FALSE.

!   computation of the actual surface
  present_surf=0.
  do i=1,ngrid
      if (initial_h2o_ice(i).GT.0.5 .and. qsurf(i).GT.0.) then
         print *, "i", i
         print *, "initial_h2o_ice(i)", initial_h2o_ice(i)
         print *, "qsurf(i)", qsurf(i)
         print *, "cell_area(i)", cell_area(i)
         print *, "present_surf",present_surf
         present_surf=present_surf+cell_area(i)
      endif
  enddo

!  print *, "initial_h2o_ice", initial_h2o_ice
!  print *, "qsurf", qsurf

  print *, "present_surf", present_surf
  print *, "ini_surf", ini_surf
  print *, "ini_surf*0.8", ini_surf*(1-alpha_criterion)
  
!   check of the criterion
  if(present_surf.LT.ini_surf*(1-alpha_criterion) .OR. &
     present_surf.GT.ini_surf*(1+alpha_criterion)) then
  STOPPING=.TRUE. 
  endif

  if (ini_surf.LT. 1E-5 .and. ini_surf.GT. -1E-5) then
       STOPPING=.FALSE.
  endif

END SUBROUTINE criterion_ice_stop_water





