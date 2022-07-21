!
! $Id: conf_evol.F $

SUBROUTINE conf_evol( tapedef, etatinit )

#ifdef CPP_IOIPSL
  use IOIPSL
#else
  ! if not using IOIPSL, we still need to use (a local version of) getin
  use ioipsl_getincom
#endif

  USE temps_mod_evol, ONLY: nyear, dt_pem, alpha_criterion

  IMPLICIT NONE
!-----------------------------------------------------------------------
! Read the run_pem.def file
!
!     Arguments :
!
!     nyear   : Maximun number of year over which the PEM can interpolate
!     dt_pem  : in years, the time step used by the PEM    
!     alpha_criterion : percentage of change before stopping the PEM 
!
  LOGICAL,INTENT(IN) :: etatinit
  INTEGER,INTENT(IN) :: tapedef

!   Declarations :
!   --------------
  include "dimensions.h"
  include "paramet.h"
  include "comdissnew.h"
  include "iniprint.h"

!   local:
!   ------

  REAL clonn,clatt,grossismxx,grossismyy
  REAL dzoomxx,dzoomyy, tauxx,tauyy
  LOGICAL  fxyhypbb, ysinuss
  LOGICAL use_filtre_fft


  nyear=100
  CALL getin('nyear', nyear)

  dt_pem=1
  CALL getin('dt_pem', dt_pem)

  alpha_criterion=0.2
  CALL getin('alpha_criterion', alpha_criterion)


END SUBROUTINE conf_evol
