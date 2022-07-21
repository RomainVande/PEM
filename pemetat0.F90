!
! $Id $
!
SUBROUTINE pemetat0(fichnom,min_h2o_ice_s,min_co2_ice_s,iim_input,jjm_input,nlayer,vmr_co2_gcm,ps_GCM,timelen, & 
             min_co2_ice_slope,min_h2o_ice_slope,nslope,tsurf_ave,tsoil_ave,tsurf_gcm,tsoil_gcm,TI_ave,q_co2_GCM,q_h2o_GCM,co2_ice_slope)

      use netcdf, only: nf90_open,NF90_NOWRITE,nf90_noerr,nf90_strerror, &
                        nf90_get_var, nf90_inq_varid, nf90_inq_dimid, &
                        nf90_inquire_dimension,nf90_close
      use comsoil_h, only: nsoilmx

      IMPLICIT NONE

!=======================================================================
!
! Read initial confitions file
!
!=======================================================================

  include "dimensions.h"

!===============================================================================
! Arguments:
  CHARACTER(LEN=*), INTENT(IN) :: fichnom          !--- FILE NAME

  INTEGER :: iim_input,jjm_input,nlayer,nslope
  REAL, ALLOCATABLE ::  h2o_ice_s(:,:,:)                       ! h2o_ice_s of the concatenated file
  REAL, ALLOCATABLE ::  co2_ice_s(:,:,:)                       ! co2_ice_s of the concatenated file

  REAL, ALLOCATABLE ::  h2o_ice_s_slope(:,:,:,:)                       ! co2_ice_s of the concatenated file

  REAL, INTENT(OUT) ::  min_h2o_ice_s(iim_input+1,jjm_input+1) ! Minimum of h2o_ice_s of the year
  REAL, INTENT(OUT) ::  min_co2_ice_s(iim_input+1,jjm_input+1) ! Minimum of co2_ice_s of the year
  REAL, INTENT(OUT) ::  min_co2_ice_slope(iim_input+1,jjm_input+1,nslope) ! Minimum of co2_ice slope of the year
  REAL, INTENT(OUT) ::  min_h2o_ice_slope(iim_input+1,jjm_input+1,nslope) ! Minimum of co2_ice slope of the year
!  REAL, ALLOCATABLE ::  vmr_co2_gcm(:,:,:)                     !!!!vmr_co2_phys_gcm(iim_input+1,jjm_input+1,timelen)
  REAL, INTENT(OUT) ::  vmr_co2_gcm(iim_input+1,jjm_input+1,2676)                     !!!!vmr_co2_phys_gcm(iim_input+1,jjm_input+1,timelen)
!  REAL, ALLOCATABLE ::  q_h2o_GCM(:,:,:)
  REAL, INTENT(OUT) ::  q_h2o_GCM(iim_input+1,jjm_input+1,2676)
  REAL, INTENT(OUT) ::  q_co2_GCM(iim_input+1,jjm_input+1,2676)
!  REAL, ALLOCATABLE ::  q_co2_GCM(:,:,:)
  REAL, ALLOCATABLE ::  q1_co2_GCM(:,:,:)
!  real, INTENT(OUT) ::  vmr_co2_phys_gcm(:,:)                  !!!!vmr_co2_gcm(ngrid,timelen)
!  REAL, ALLOCATABLE ::  ps_GCM(:,:,:)
  REAL,  INTENT(OUT) ::  ps_GCM(iim_input+1,jjm_input+1,2676)


!SOIL
  REAL, INTENT(OUT) ::  tsurf_ave(iim_input+1,jjm_input+1,nslope) ! Average surface temperature of the concatenated file
  REAL, INTENT(OUT) ::  tsoil_ave(iim_input+1,jjm_input+1,nsoilmx,nslope) ! Average soil temperature of the concatenated file

  REAL ,INTENT(OUT) ::  tsurf_gcm(iim_input+1,jjm_input+1,nslope,2676) ! Surface temperature of the concatenated file
  REAL , INTENT(OUT) ::  tsoil_gcm(iim_input+1,jjm_input+1,nsoilmx,nslope,2676) ! Soil temperature of the concatenated file

  REAL ::  TI_gcm(iim_input+1,jjm_input+1,nsoilmx,nslope,2676) ! Thermal Inertia  of the concatenated file
  REAL, INTENT(OUT) ::  TI_ave(iim_input+1,jjm_input+1,nsoilmx,nslope) ! Average Thermal Inertia  of the concatenated file
  REAL, INTENT(OUT) ::  co2_ice_slope(iim_input+1,jjm_input+1,nslope,2676) ! Minimum of co2_ice slope of the year
!===============================================================================
!   Local Variables 
  CHARACTER(LEN=256) :: msg, var, modname
  INTEGER,PARAMETER :: length=100
  INTEGER :: iq, fID, vID, idecal
  INTEGER :: ierr
  CHARACTER(len=12) :: start_file_type="earth" ! default start file type

  REAL,ALLOCATABLE :: time(:) ! times stored in start
  INTEGER :: timelen ! number of times stored in the file
  INTEGER :: indextime ! index of selected time

  INTEGER :: edges(4),corner(4)
  INTEGER :: i,j,t
  real,save :: m_co2, m_noco2, A , B, mmean

  INTEGER :: islope
  CHARACTER*2 :: num


!-----------------------------------------------------------------------
  modname="pemetat0"

      m_co2 = 44.01E-3  ! CO2 molecular mass (kg/mol)   
      m_noco2 = 33.37E-3  ! Non condensible mol mass (kg/mol)   
      A =(1/m_co2 - 1/m_noco2)
      B=1/m_noco2

!  Open initial state NetCDF file
  var=fichnom
  CALL err(NF90_OPEN(var,NF90_NOWRITE,fID),"open",var)

      ierr = nf90_inq_varid (fID, "temps", vID)
      IF (ierr .NE. nf90_noerr) THEN
        write(*,*)"pemetat0: Le champ <temps> est absent"
        write(*,*)"pemetat0: J essaie <Time>"
        ierr = nf90_inq_varid (fID, "Time", vID)
        IF (ierr .NE. nf90_noerr) THEN
           write(*,*)"pemetat0: Le champ <Time> est absent"
           write(*,*)trim(nf90_strerror(ierr))
     print *, "ICIIII0"
           CALL ABORT_gcm("pemetat0", "", 1)
        ENDIF
        ! Get the length of the "Time" dimension
     print *, "ICIIIITIME"
        ierr = nf90_inq_dimid(fID,"Time",vID)
        ierr = nf90_inquire_dimension(fID,vID,len=timelen)
      ELSE   
        ! Get the length of the "temps" dimension
     print *, "ICIIIITEMPS"
        ierr = nf90_inq_dimid(fID,"temps",vID)
        ierr = nf90_inquire_dimension(fID,vID,len=timelen)
      ENDIF

     print *, "ICIIII"

      allocate(co2_ice_s(iim+1,jjm+1,timelen))

     print *, "ICIIIIAAAA"

!      allocate(q_co2_GCM(iim+1,jjm+1,timelen))

     print *, "ICIIIIBBBB"

!      allocate(q_h2o_GCM(iim+1,jjm+1,timelen))

     print *, "ICIIIICCCC"

      allocate(q1_co2_GCM(iim+1,jjm+1,timelen))

     print *, "ICIIII2"



      allocate(h2o_ice_s_slope(iim+1,jjm+1,nslope,timelen))

          print *, "ICIIII3"

! Get h2o_ice_s of the concatenated file
  CALL get_var3("h2o_ice_s"   ,h2o_ice_s)

  print *, "A"

  CALL get_var3("co2ice"   ,co2_ice_s)
  CALL get_var3("co2_cropped"   ,q_co2_GCM)
  CALL get_var3("h2o_cropped"   ,q_h2o_GCM)

  print *, "B"

  CALL get_var3("ps"   ,ps_GCM)

  print *, "C"

  print *, "nslope", nslope

DO islope=1,nslope
  write(num,fmt='(i2.2)') islope
  call get_var3("co2ice_slope"//num,co2_ice_slope(:,:,islope,:))
ENDDO

  print *, "co2ice_slope"

DO islope=1,nslope
  write(num,fmt='(i2.2)') islope
  call get_var3("h2o_ice_s_slope"//num,h2o_ice_s_slope(:,:,islope,:))
ENDDO

  print *, "h2o_ice_s_slope"

DO islope=1,nslope
  write(num,fmt='(i2.2)') islope
  call get_var3("tsurf_slope"//num,tsurf_gcm(:,:,islope,:))
ENDDO

  print *, "tsurf_slope"

DO islope=1,nslope
  write(num,fmt='(i2.2)') islope
  call get_var4("tsoil_slope"//num,tsoil_gcm(:,:,:,islope,:))
ENDDO

  print *, "tsoil_slope"

DO islope=1,nslope
  write(num,fmt='(i2.2)') islope
  call get_var4("inertiesoil_slope"//num,TI_gcm(:,:,:,islope,:))
ENDDO

  print *, "inertiesoil_slope"





! Compute the minimum over the year for each point
  min_h2o_ice_s(:,:)=minval(h2o_ice_s,3)
  min_co2_ice_s(:,:)=minval(co2_ice_s,3)

  min_co2_ice_slope(:,:,:)=minval(co2_ice_slope,4)
  min_h2o_ice_slope(:,:,:)=minval(h2o_ice_s_slope,4)

!Compute averages

!  DO i=1,timelen
    tsurf_ave(:,:,:)=SUM(tsurf_gcm(:,:,:,:),4)/timelen
    tsoil_ave(:,:,:,:)=SUM(tsoil_gcm(:,:,:,:,:),5)/timelen
    TI_ave(:,:,:,:)=SUM(TI_gcm(:,:,:,:,:),5)/timelen
!  ENDDO




! By definition, a density is positive, we get rid of the negative values
  DO i=1,iim+1
    DO j = 1, jjm+1
       if (min_co2_ice_s(i,j).LT.0) then
          min_h2o_ice_s(i,j)  = 0.
          min_co2_ice_s(i,j)  = 0.
       endif
       DO islope=1,nslope
          if (min_co2_ice_slope(i,j,islope).LT.0) then
            min_co2_ice_slope(i,j,islope)  = 0.
          endif
          if (min_h2o_ice_slope(i,j,islope).LT.0) then
            min_h2o_ice_slope(i,j,islope)  = 0.
          endif
       ENDDO
    ENDDO
  ENDDO

  DO i=1,iim+1
    DO j = 1, jjm+1
      DO t = 1, timelen
!         q1_co2_GCM(i,j,t)=q_co2_GCM(i,j,t)
!         ps_GCM(i,j,t)=ps(i)
!c             Mean air molecular mass = 1/(q(ico2)/m_co2 + (1-q(ico2))/m_noco2)
         if (q_co2_GCM(i,j,t).LT.0) then
              q_co2_GCM(i,j,t)=1E-10
         elseif (q_co2_GCM(i,j,t).GT.1) then
              q_co2_GCM(i,j,t)=1.
         endif
         if (q_h2o_GCM(i,j,t).LT.0) then
              q_h2o_GCM(i,j,t)=1E-30
         elseif (q_h2o_GCM(i,j,t).GT.1) then
              q_h2o_GCM(i,j,t)=1.
         endif
         mmean=1/(A*q_co2_GCM(i,j,t) +B)
         vmr_co2_gcm(i,j,t) = q_co2_GCM(i,j,t)*mmean/m_co2
      ENDDO
    ENDDO
  ENDDO






  CONTAINS

SUBROUTINE check_dim(n1,n2,str1,str2)
  INTEGER,          INTENT(IN) :: n1, n2
  CHARACTER(LEN=*), INTENT(IN) :: str1, str2
  CHARACTER(LEN=256) :: s1, s2
  IF(n1/=n2) THEN
    s1='value of '//TRIM(str1)//' ='
    s2=' read in starting file differs from parametrized '//TRIM(str2)//' ='
    WRITE(msg,'(10x,a,i4,2x,a,i4)')TRIM(s1),n1,TRIM(s2),n2
    CALL ABORT_gcm(TRIM(modname),TRIM(msg),1)
  END IF
END SUBROUTINE check_dim


SUBROUTINE get_var1(var,v)
  CHARACTER(LEN=*), INTENT(IN)  :: var
  REAL,             INTENT(OUT) :: v(:)
  CALL err(NF90_INQ_VARID(fID,var,vID),"inq",var)
  CALL err(NF90_GET_VAR(fID,vID,v),"get",var)
END SUBROUTINE get_var1


SUBROUTINE get_var3(var,v) ! on U grid
  CHARACTER(LEN=*), INTENT(IN)  :: var
  REAL,             INTENT(OUT) :: v(:,:,:)
  CALL err(NF90_INQ_VARID(fID,var,vID),"inq",var)
  CALL err(NF90_GET_VAR(fID,vID,v),"get",var)

END SUBROUTINE get_var3

SUBROUTINE get_var4(var,v) 
  CHARACTER(LEN=*), INTENT(IN)  :: var
  REAL,             INTENT(OUT) :: v(:,:,:,:)
  CALL err(NF90_INQ_VARID(fID,var,vID),"inq",var)
  CALL err(NF90_GET_VAR(fID,vID,v),"get",var)
END SUBROUTINE get_var4

SUBROUTINE err(ierr,typ,nam)
  INTEGER,          INTENT(IN) :: ierr   !--- NetCDF ERROR CODE
  CHARACTER(LEN=*), INTENT(IN) :: typ    !--- TYPE OF OPERATION
  CHARACTER(LEN=*), INTENT(IN) :: nam    !--- FIELD/FILE NAME
  IF(ierr==NF90_NoERR) RETURN
  SELECT CASE(typ)
    CASE('inq');   msg="Field <"//TRIM(nam)//"> is missing"
    CASE('get');   msg="Reading failed for <"//TRIM(nam)//">"
    CASE('open');  msg="File opening failed for <"//TRIM(nam)//">"
    CASE('close'); msg="File closing failed for <"//TRIM(nam)//">"
  END SELECT
  CALL ABORT_gcm(TRIM(modname),TRIM(msg),ierr)
END SUBROUTINE err

END SUBROUTINE pemetat0
