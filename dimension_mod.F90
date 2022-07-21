MODULE dimension_mod

    IMPLICIT NONE
    INTEGER,SAVE    :: x 
!$OMP THREADPRIVATE(x)
    INTEGER,SAVE    :: y
!$OMP THREADPRIVATE(y)
    INTEGER,SAVE    :: z     
!$OMP THREADPRIVATE(z)

CONTAINS

  SUBROUTINE init_dim(x_, y_, z_)
    USE ioipsl_getin_p_mod, ONLY : getin_p
    IMPLICIT NONE
    INTEGER,INTENT(IN) :: x_
    INTEGER,INTENT(IN) :: y_
    INTEGER,INTENT(IN) :: z_
    
    x=x_
    y=y_
    z=z_
    
  END SUBROUTINE init_dim

END MODULE dimension_mod      
