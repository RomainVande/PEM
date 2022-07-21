MODULE temps_mod_evol

IMPLICIT NONE  

  INTEGER   nyear           !     nyear   : Maximun number of year over which the PEM can interpolate
  INTEGER   dt_pem          !     dt_pem  : in years, the time step used by the PEM    
  REAL      alpha_criterion !     alpha_criterion : percentage of change before stopping the PEM 

!$OMP THREADPRIVATE(nyear)        

!WARNING: when adding a threadprivate variable in this module
!        do not forget to add it to the copyin clause when opening an OpenMP
!        parallel section. e.g. in gcm before call leapfrog_loc and/or
!        possibly in iniphysiq

END MODULE temps_mod_evol
