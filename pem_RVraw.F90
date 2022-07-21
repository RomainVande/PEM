
!------------------------

! I   Initialisation
!    I_a READ run.def , run_pem.def
!    I_b READ starfi_0.nc

! II  Run
!    II_a Compute tendencies
!    II_b Save initial situation
!    II_c Time loop

! III Output
!    III_a Write startfi.nc

!------------------------


PROGRAM pem

!module needed for INITIALISATION
      use phyetat0_mod, only: phyetat0
      use comsoil_h, only: tsoil, nsoilmx, ini_comsoil_h,inertiedat, mlayer,volcapa
      use surfdat_h, only: tsurf, co2ice, emis,&
      &                    qsurf,watercap, ini_surfdat_h, &
                           albedodat, zmea, zstd, zsig, zgam, zthe, &
      &                    hmons, summit, base
      use dimradmars_mod, only: totcloudfrac, albedo, & 
                                ini_dimradmars_mod
      use turb_mod, only: q2, wstar, ini_turb_mod
      use dust_param_mod, only: tauscaling, ini_dust_param_mod
!      use co2cloud_mod, only: mem_Mccn_co2, mem_Mh2o_co2,&
!      &                        mem_Nccn_co2,ini_co2cloud
      use netcdf, only: nf90_open,NF90_NOWRITE,nf90_noerr,nf90_strerror, &
                        nf90_get_var, nf90_inq_varid, nf90_inq_dimid, &
                        nf90_inquire_dimension,nf90_close
      use phyredem, only: physdem0, physdem1
      use tracer_mod, only: noms,nqmx ! tracer names

! For phyredem :
      USE control_mod, ONLY: iphysiq, day_step,nsplit_phys
      use time_phylmdz_mod, only: daysec
      use mod_phys_lmdz_para, only: is_parallel, is_sequential, &
                                   is_mpi_root, is_omp_root,    &
                                   is_master
      use time_phylmdz_mod, only: dtphys
      USE mod_const_mpi, ONLY: COMM_LMDZ
      USE comconst_mod, ONLY: rad,g,r,cpp
      USE logic_mod, ONLY: iflag_phys
      USE iniphysiq_mod, ONLY: iniphysiq
      USE infotrac
      USE temps_mod_evol, ONLY: nyear, dt_pem
!     USE vertical_layers_mod, ONLY: ap,bp

      USE comslope_mod, ONLY: nslope,def_slope,def_slope_mean, &
                           subslope_dist,co2ice_slope, &
                           tsurf_slope,tsoil_slope,fluxgrd_slope,&
                           fluxrad_sky_slope,sky_slope,callsubslope,&
                           co2iceflow, beta_slope, capcal_slope,&
                           albedo_slope,emiss_slope,qsurf_slope,&
                           iflat

      USE geometry_mod, only: longitude_deg,latitude_deg 


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! SOIL

      use comsoil_h_PEM, only: ini_comsoil_h_PEM,end_comsoil_h_PEM,nsoilmx_PEM, &
                              TI_PEM,alph_PEM, beta_PEM, &        ! soil thermal inertia          
                              tsoil_PEM ,       &        !number of subsurface layers
                              mlayer_PEM,layer_PEM, &       ! soil mid layer depths
                              fluxgeo ! geothermal flux

  IMPLICIT NONE

  include "dimensions.h"
  include "paramet.h"
  include "comdissnew.h"
  include "comgeom.h"
  include "iniprint.h"

! Same variable's name as in the GCM

      INTEGER :: ngrid      !Number of physical grid points
      INTEGER :: nlayer     !Number of vertical layer
      INTEGER :: nq         !Number of tracer
      INTEGER :: day_ini    !First day of the simulation
      REAL :: pday          !Physical day
      REAL :: time_phys     !Same as GCM
      REAL :: ptimestep     !Same as GCM
      REAL :: ztime_fin     !Same as GCM

! Variable for reading start.nc
      character (len = *), parameter :: FILE_NAME_start = "start_0.nc" !Name of the file used for initialsing the PEM
  !   variables dynamiques
  REAL vcov(ip1jm,llm),ucov(ip1jmp1,llm) ! vents covariants
  REAL teta(ip1jmp1,llm)                 ! temperature potentielle 
  REAL, ALLOCATABLE, DIMENSION(:,:,:):: q! champs advectes
  REAL ps(ip1jmp1)                       ! pression  au sol
  REAL, dimension(:),allocatable :: ps_phys !(ngrid)                       ! pression  au sol
!  REAL p (ip1jmp1,llmp1  )               ! pression aux interfac.des couches
  REAL masse(ip1jmp1,llm)                ! masse d'air
  REAL phis(ip1jmp1)                     ! geopotentiel au sol
!  REAL phi(ip1jmp1,llm)                  ! geopotentiel
!  REAL w(ip1jmp1,llm)                    ! vitesse verticale
  REAL time_0

! Variable for reading starfi.nc

      character (len = *), parameter :: FILE_NAME = "startfi_0.nc" !Name of the file used for initialsing the PEM
      integer :: ncid, varid,status                                !Variable for handling opening of files
      integer :: phydimid, subdimid, nlayerdimid, nqdimid          !Variable ID for Netcdf files
      integer :: lonvarid, latvarid, areavarid,sdvarid             !Variable ID for Netcdf files
      integer :: apvarid,bpvarid				   !Variable ID for Netcdf files

! Variable for reading starfi.nc and writting restartfi.nc

      REAL, dimension(:),allocatable :: longitude                  !Longitude read in FILE_NAME and written in restartfi 
      REAL, dimension(:),allocatable  :: latitude                  !Latitude read in FILE_NAME and written in restartfi 
      REAL, dimension(:),allocatable :: ap                         !Coefficient ap read in FILE_NAME_start and written in restart 
      REAL, dimension(:),allocatable :: bp                         !Coefficient bp read in FILE_NAME_start and written in restart
      REAL, dimension(:),allocatable  :: cell_area                 !Cell_area read in FILE_NAME and written in restartfi 
      REAL :: Total_surface                                        !Total surface of the planet

! Variable for h2o_ice evolution

      REAL , dimension(:,:), allocatable :: tendencies_h2o_ice     ! LON x LAT field : Tendency of evolution of perenial ice over a year
      REAL, dimension(:),allocatable  :: tendencies_h2o_ice_phys   ! physical point field : Tendency of evolution of perenial ice over a year

      REAL , dimension(:,:), allocatable :: tendencies_co2_ice     ! LON x LAT field : Tendency of evolution of perenial co2 ice over a year
      REAL, dimension(:),allocatable  :: tendencies_co2_ice_phys   ! physical point field : Tendency of evolution of perenial co2 ice over a year

      REAL :: ini_surf                                             ! Initial surface of sublimating water ice
      REAL :: ini_surf_h2o                                             ! Initial surface of sublimating water ice
      REAL, dimension(:),allocatable  :: initial_h2o_ice           ! physical point field : Logical array indicating sublimating point

      REAL :: ini_surf_co2                                         ! Initial surface of sublimating co2 ice
      REAL, dimension(:),allocatable  :: initial_co2_ice           ! physical point field : Logical array indicating sublimating point of co2 ice

      REAL , dimension(:,:), allocatable :: min_h2o_ice_s_1        ! LON x LAT field : minimum of water ice at each point for the first year
      REAL , dimension(:,:), allocatable :: min_h2o_ice_s_2        ! LON x LAT field : minimum of water ice at each point for the second year
      REAL , dimension(:,:,:), allocatable :: h2o_ice_first_last_day     ! LON x LAT x 2 field : First and Last value of a GCM year simulation of water ice

      REAL , dimension(:,:), allocatable :: min_co2_ice_s_1        ! LON x LAT field : minimum of water ice at each point for the first year
      REAL , dimension(:,:), allocatable :: min_co2_ice_s_2        ! LON x LAT field : minimum of water ice at each point for the second year
      REAL , dimension(:,:,:), allocatable :: co2_ice_first_last_day     ! LON x LAT x 2 field : First and Last value of a GCM year simulation of water ice

      REAL, dimension(:),allocatable  :: local_old_press           ! physical point field : Local pressure of initial/previous time step
      REAL, dimension(:),allocatable  :: local_new_press           ! physical point field : Local pressure of current time step

      REAL :: global_ave_press_GCM
      REAL :: global_ave_press_old           ! physical point field : Global average pressure of initial/previous time step
      REAL :: global_ave_press_new           ! physical point field : Global average pressure of current time step

      REAL , dimension(:,:), allocatable ::  zplev_new
      REAL , dimension(:,:), allocatable :: zplev_old

      REAL :: tot_co2_atm,tot_var_co2_atm

!      INTEGER :: n_band_lat                    ! Number of latitude band to estimate the stopping criterion for co2_ice
!      REAL, dimension(:),allocatable  :: initial_co2_ice ! Initial amount of co2_ice per lat band

      INTEGER :: year_iter  !Counter for the number of PEM iteration
      LOGICAL :: STOPPING_water   ! Logical : is the criterion (% of change in the surface of sublimating water ice) reached?
      LOGICAL :: STOPPING_1_water ! Logical : is there still water ice to sublimate?
      LOGICAL :: STOPPING_co2   ! Logical : is the criterion (% of change in the surface of sublimating water ice) reached?
      LOGICAL :: STOPPING_1_co2 ! Logical : is there still water ice to sublimate?

      REAL, dimension(:,:,:),allocatable  :: q_co2_GCM ! Initial amount of co2 in the first layer
      REAL ps_GCM(ip1jmp1)                       ! pression  au sol donné par le GCM
      real,save :: m_co2, m_noco2, A , B, mmean
      real ,allocatable :: vmr_co2_gcm_phys(:,:) !(ngrid) ! co2 volume mixing ratio
      real ,allocatable :: vmr_co2_pem_phys(:,:) !(ngrid) ! co2 volume mixing ratio
      real ,allocatable :: q_h2o_GCM_phys(:,:)
      real ,allocatable :: q_co2_GCM_phys(:,:)
      real ,allocatable :: q_co2_PEM_phys(:,:)
      real ,allocatable :: q_co2_PEM_phys_ave(:)
      REAL, ALLOCATABLE ::  ps_GCM_1(:,:,:)
      REAL, ALLOCATABLE ::  ps_GCM_2(:,:,:)
      REAL, ALLOCATABLE ::  T_cond_GCM(:,:,:)
      REAL, ALLOCATABLE ::  vmr_co2_gcm(:,:,:)
      REAL, ALLOCATABLE ::  q_h2o_GCM(:,:,:)
      REAL ,allocatable ::  q_h2o_PEM_phys(:,:)
      integer :: timelen
      REAL :: ave

      REAL, ALLOCATABLE :: p(:,:)  !(ngrid,llmp1)






!!!!!!!!!!!!!!!!!!!!!!!! SLOPE
      character*2 :: str2
      REAL ,allocatable :: watercap_slope(:,:)    !(ngrid,nslope)
      REAL , dimension(:,:,:), allocatable :: min_co2_ice_slope_1        ! LON x LAT field : minimum of water ice at each point for the first year
      REAL , dimension(:,:,:), allocatable :: min_co2_ice_slope_2        ! LON x LAT field : minimum of water ice at each point for the second year
      REAL , dimension(:,:,:), allocatable :: min_h2o_ice_slope_1        ! LON x LAT field : minimum of water ice at each point for the first year
      REAL , dimension(:,:,:), allocatable :: min_h2o_ice_slope_2        ! LON x LAT field : minimum of water ice at each point for the second year

      REAL, dimension(:,:),allocatable  :: initial_co2_ice_sublim_slope           ! physical point field : Logical array indicating sublimating point of co2 ice
      REAL, dimension(:,:),allocatable  :: initial_h2o_ice_slope           ! physical point field : Logical array indicating sublimating point of h2o ice
      REAL, dimension(:,:),allocatable  :: initial_co2_ice_slope           ! physical point field : Logical array indicating sublimating point of co2 ice

      REAL , dimension(:,:,:), allocatable :: tendencies_co2_ice_slope     ! LON x LAT field : Tendency of evolution of perenial co2 ice over a year
      REAL , dimension(:,:,:), allocatable :: tendencies_h2o_ice_slope     ! LON x LAT field : Tendency of evolution of perenial co2 ice over a year
      REAL, dimension(:,:),allocatable  :: tendencies_co2_ice_phys_slope   ! physical point field : Tendency of evolution of perenial co2 ice over a year
      REAL, dimension(:,:),allocatable  :: tendencies_co2_ice_phys_slope_ini ! physical point field x nslope: Tendency of evolution of perenial co2 ice over a year in the GCM
      REAL, dimension(:,:),allocatable  :: tendencies_h2o_ice_phys_slope   ! physical point field : Tendency of evolution of perenial co2 ice over a year

      REAL, PARAMETER :: co2_hmax = 10 			  ! Maximum height  for CO2 deposit on slopes (m)
      REAL, PARAMETER :: rho_co2 = 1600           ! CO2 ice density (kg/m^3)
      INTEGER :: iaval                            ! Index of the neighboord slope ()
      REAL , dimension(:,:), allocatable :: flag_co2flow(:,:)   !(ngrid,nslope)          ! To flag where there is a glacier flow
      REAL , dimension(:), allocatable :: flag_co2flow_mesh(:)  !(ngrid)          ! To flag where there is a glacier flow


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! SURFACE/SOIL 


     REAL, ALLOCATABLE :: tsurf_ave(:,:,:) ! LON x LAT x SLOPE field : Averaged Surface Temperature [K]
     REAL, ALLOCATABLE  :: tsurf_ave_phys(:,:) ! IG x LAT x SLOPE field : Averaged Surface Temperature [K]
     REAL, ALLOCATABLE :: tsoil_ave(:,:,:,:) ! LON x LAT x SLOPE field : Averaged Soil Temperature [K]
     REAL, ALLOCATABLE :: tsoil_ave_phys(:,:,:,:) !IG x SLOPE field : Averaged Soil Temperature [K]
     REAL, ALLOCATABLE :: TI_GCM_ave(:,:,:,:) ! LON x LAT x SLOPE field : Averaged Thermal Inertia  [SI]
     REAL, ALLOCATABLE :: tsurf_ave_phys_inst(:,:)
     REAL, ALLOCATABLE :: tsoil_ave_phys_inst(:,:,:)

     LOGICAL :: firstcall
!     REAL, PARAMETER :: daysec=88775.              ! duree du sol (s)  ~88775 s
     REAL, PARAMETER :: year_day = 669
     REAL, PARAMETER :: year_step = 1
     REAL :: timestep

     REAL, ALLOCATABLE :: TI_GCM_phys(:,:,:) ! Averaged GCM Thermal Inertia  [SI]
     REAL, ALLOCATABLE :: TI_GCM_start(:,:,:) ! Averaged GCM Thermal Inertia  [SI]

     REAL,ALLOCATABLE  :: q_h2o_PEM_phys_ave(:) ! averaged water vapor content
     REAL,ALLOCATABLE  :: interp_coef(:)

     REAL,ALLOCATABLE  :: ice_depth(:,:)
     REAL,ALLOCATABLE  :: TI_locslope(:,:)
     REAL,ALLOCATABLE  :: Tsoil_locslope(:,:)
     REAL,ALLOCATABLE  :: Tsurf_locslope(:)
     REAL,ALLOCATABLE  :: alph_locslope(:,:)
     REAL,ALLOCATABLE  :: beta_locslope(:,:)   
     REAL :: kcond

!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! Loop variable
     LOGICAL :: bool_sublim
     INTEGER :: i,j,ig0,l,ig,nnq,t,islope,ig_loop,islope_loop,iloop
     REAL :: pi,beta,alpha

! Parallel variables
      is_sequential=.true.
      is_parallel=.false.
      is_mpi_root=.true.
      is_omp_root=.true.
      is_master=.true.

      day_ini=0    !test
      time_phys=0. !test

!      n_band_lat=18
!      pi=4.D0*DATAN(1.D0)
      m_co2 = 44.01E-3  ! CO2 molecular mass (kg/mol)   
      m_noco2 = 33.37E-3  ! Non condensible mol mass (kg/mol)   
      A =(1/m_co2 - 1/m_noco2)
      B=1/m_noco2



!------------------------

! I   Initialisation
!    I_a READ run.def 


!------------------------

!----------------------------READ run.def ---------------------
      CALL conf_gcm( 99, .TRUE. )
      CALL conf_evol( 99, .TRUE. )


!------------------------

! I   Initialisation
!    I_a READ run.def 
!    I_b READ starfi_0.nc

!----------------------------READ startfi.nc --------------------- 

!----------------------------initialisation --------------------- 

      status =nf90_open(FILE_NAME, NF90_NOWRITE, ncid)
      status =nf90_inq_dimid(ncid, "physical_points", phydimid)
      status =nf90_inquire_dimension(ncid, phydimid, len = ngrid)

      status =nf90_inq_dimid(ncid, "nlayer", nlayerdimid)
      status =nf90_inquire_dimension(ncid, nlayerdimid, len = nlayer)

      status =nf90_inq_dimid(ncid,"number_of_advected_fields",nqdimid)
      status =nf90_inquire_dimension(ncid, nqdimid, len = nq)

      allocate(longitude(ngrid))
      allocate(latitude(ngrid))
      allocate(cell_area(ngrid))

      status = nf90_inq_varid(ncid, "longitude", lonvarid)
      status = nf90_get_var(ncid, lonvarid, longitude)

      status = nf90_inq_varid(ncid, "latitude", latvarid)
      status = nf90_get_var(ncid, latvarid, latitude)

      status = nf90_inq_varid(ncid, "area", areavarid)
      status = nf90_get_var(ncid, areavarid, cell_area)

      call ini_comsoil_h(ngrid)

      status = nf90_inq_varid(ncid, "soildepth", sdvarid)
      status = nf90_get_var(ncid, sdvarid, mlayer)

      status =nf90_close(ncid)



!----------------------------READ start.nc --------------------- 

    call infotrac_init

     allocate(q(ip1jmp1,llm,nqtot))

!     CALL iniacademic(vcov,ucov,teta,q,masse,ps,phis,time_0)

     CALL dynetat0(FILE_NAME_start,vcov,ucov, &
                    teta,q,masse,ps,phis, time_0)
  

     CALL iniconst !new
     CALL inigeom
!     CALL inifilr !new

     allocate(ap(nlayer+1))
     allocate(bp(nlayer+1))

      status =nf90_open(FILE_NAME_start, NF90_NOWRITE, ncid)

      status = nf90_inq_varid(ncid, "ap", apvarid)
      status = nf90_get_var(ncid, apvarid, ap)
      status = nf90_inq_varid(ncid, "bp", bpvarid)
      status = nf90_get_var(ncid, bpvarid, bp)

      status =nf90_close(ncid)
     
   
!    daysec=1 !test
!    dtphys=iphysiq*daysec/REAL(day_step) !test
!    rad=1. !test
!    g=1. !test
!    r=1. !test?
!    cpp=1. !test

!    daysec=10 !test
!    dtphys=iphysiq*daysec/REAL(day_step) !test
!    rad=3397200. !test
!    g=3.72000002861023 !test
!    r=10. !test
!    cpp=744.499 !test




    CALL iniphysiq(iim,jjm,llm, &
          (jjm-1)*iim+2,comm_lmdz, &
          daysec,day_ini,dtphys/nsplit_phys, &
          rlatu,rlatv,rlonu,rlonv,aire,cu,cv,rad,g,r,cpp, &
          iflag_phys)


!----------------------------reading ---------------------


! First we read the initial state (starfi.nc)

!    CALL phyetat0 (FILE_NAME,0,0, &
!              nsoilmx,ngrid,nlayer,nq, &
!              day_ini,time_phys, &
!              tsurf,tsoil,albedo,emis, &
!              q2,qsurf,co2ice,tauscaling,totcloudfrac,wstar, &
!              mem_Mccn_co2,mem_Nccn_co2, &
!              mem_Mh2o_co2,watercap) 

       allocate(watercap_slope(ngrid,nslope))
       allocate(TI_GCM_start(ngrid,nsoilmx,nslope))
       allocate(q_h2o_PEM_phys_ave(ngrid))



         CALL phyetat0 (FILE_NAME,0,0, &
              nsoilmx,ngrid,nlayer,nq,   &
              day_ini,time_phys,         &
              tsurf,tsoil,albedo,emis,   &
              q2,qsurf,co2ice,tauscaling,totcloudfrac,wstar,     &
              watercap,nslope,tsurf_slope,                       &
              tsoil_slope,co2ice_slope,def_slope,def_slope_mean, &
              subslope_dist,albedo_slope,emiss_slope, TI_GCM_start,     &
              qsurf_slope,watercap_slope)

       iflat=1
       DO islope=2,nslope
         IF(abs(def_slope_mean(islope)).lt. &
           abs(def_slope_mean(iflat)))THEN
           iflat = islope
         ENDIF     

       ENDDO
       PRINT*,'Flat slope for islope = ',iflat
       PRINT*,'corresponding criterium = ',def_slope_mean(iflat)

       flag_co2flow(:,:) = 0.     
       flag_co2flow_mesh(:) = 0.




!!!!!!!!!!!!!!!!!!!!!

! Then we read the evolutaion of water and co2 ice over the first year of the GCM run, saving only the minimum value

     allocate(min_h2o_ice_s_1(iim+1,jjm+1))
     allocate(min_co2_ice_s_1(iim+1,jjm+1))
     allocate(vmr_co2_gcm(iim+1,jjm+1,2676))
     allocate(q_h2o_GCM(iim+1,jjm+1,2676))
     allocate(q_co2_GCM(iim+1,jjm+1,2676))
     allocate(ps_GCM_1(iim+1,jjm+1,2676))
     allocate(ps_GCM_2(iim+1,jjm+1,2676))

     allocate(min_co2_ice_slope_1(iim+1,jjm+1,nslope))
     allocate(min_h2o_ice_slope_1(iim+1,jjm+1,nslope))

     allocate(tsurf_ave(iim+1,jjm+1,nslope))
     allocate(tsoil_ave(iim+1,jjm+1,nsoilmx,nslope))
     allocate(tsurf_ave_phys(ngrid,nslope))
     allocate(tsoil_ave_phys_inst(ngrid,nsoilmx,nslope))
     allocate(tsurf_ave_phys_inst(ngrid,nslope))
     allocate(TI_GCM_ave(iim+1,jjm+1,nsoilmx,nslope))

!  call pemetat0 ("concat_year_one.nc", min_h2o_ice_s_1,min_co2_ice_s_1,iim,jjm,nlayer,vmr_co2_gcm,ps_GCM_1,timelen,min_co2_ice_slope_1,nslope)
     call pemetat0 ("concat_year_one.nc", min_h2o_ice_s_1,min_co2_ice_s_1,iim,jjm,nlayer,vmr_co2_gcm,ps_GCM_1,timelen,min_co2_ice_slope_1,min_h2o_ice_slope_1,&   
                       nslope,tsurf_ave,tsoil_ave,TI_GCM_ave,q_co2_GCM,q_h2o_GCM)

! Then we read the evolutaion of water and co2 ice over the second year of the GCM run, saving only the minimum value

     allocate(min_h2o_ice_s_2(iim+1,jjm+1))
     allocate(min_co2_ice_s_2(iim+1,jjm+1))
     allocate(min_co2_ice_slope_2(iim+1,jjm+1,nslope))
     allocate(min_h2o_ice_slope_2(iim+1,jjm+1,nslope))

!  call pemetat0 ("concat_year_two.nc", min_h2o_ice_s_2,min_co2_ice_s_2,iim,jjm,nlayer,vmr_co2_gcm,ps_GCM_2,timelen,min_co2_ice_slope_2,nslope)
  call pemetat0 ("concat_year_two.nc", min_h2o_ice_s_2,min_co2_ice_s_2,iim,jjm,nlayer,vmr_co2_gcm,ps_GCM_2,timelen,min_co2_ice_slope_2,min_h2o_ice_slope_2, &
                  nslope,tsurf_ave,tsoil_ave,TI_GCM_ave,q_co2_GCM,q_h2o_GCM)


!!!!!!!!!!!!!!!!!!!!! Initialisation for soil_PEM

      call end_comsoil_h_PEM
      call ini_comsoil_h_PEM(ngrid,nslope)

      timestep=year_day*daysec/year_step
      
     


     allocate(ice_depth(ngrid,nslope))






! We read the evolutaion of water ice over the first year of the GCM run, saving the first and last state

     allocate(vmr_co2_gcm_phys(ngrid,timelen))
     allocate(vmr_co2_pem_phys(ngrid,timelen))
     allocate(q_h2o_GCM_phys(ngrid,timelen))
     allocate(q_h2o_PEM_phys(ngrid,timelen))
     allocate(q_co2_GCM_phys(ngrid,timelen))
     allocate(q_co2_PEM_phys(ngrid,timelen))
     allocate(q_co2_PEM_phys_ave(ngrid))

  vmr_co2_gcm_phys(1,:)=vmr_co2_gcm(1,1,:)
  q_h2o_GCM_phys(1,:)=q_h2o_GCM(1,1,:)
  q_co2_GCM_phys(1,:)=q_co2_GCM(1,1,:)

  ig0 = 2
  DO j=2,jjm
    DO i = 1, iim
       vmr_co2_gcm_phys(ig0,:)  =vmr_co2_gcm(i,j,:)
       q_h2o_GCM_phys(ig0,:)  =q_h2o_GCM(i,j,:)
       q_co2_GCM_phys(ig0,:)  =q_co2_GCM(i,j,:)
       ig0= ig0 + 1
    ENDDO
  ENDDO

  vmr_co2_gcm_phys(ig0,:) = vmr_co2_gcm(1,jjm+1,:)
  q_h2o_GCM_phys(ig0,:) = q_h2o_GCM(1,jjm+1,:)
  q_co2_GCM_phys(ig0,:) = q_co2_GCM(1,jjm+1,:)
 
  q_co2_PEM_phys(:,:)=  q_co2_GCM_phys(:,:)
  q_h2o_PEM_phys(:,:)=  q_h2o_GCM_phys(:,:)
!!!!!!! Soil!!!!!
! We initialize the TI
      allocate(TI_GCM_phys(ngrid,nsoilmx,nslope))

DO l=1,nsoilmx
 TI_GCM_phys(1,l,:)=TI_GCM_ave(1,1,l,:)
ENDDO
  ig0 = 2
  DO j=2,jjm
    DO i = 1, iim
       DO l=1,nsoilmx
       TI_GCM_phys(ig0,l,:)=TI_GCM_ave(i,j,l,:)
!Tsoil_PEM=tsoil_GCM_average
       ENDDO
    ig0= ig0 + 1
    ENDDO
  ENDDO

DO l=1,nsoilmx
  TI_GCM_phys(ig0,l,:) = TI_GCM_ave(1,jjm+1,l,:)
ENDDO

      call soil_settings_PEM(ngrid,nslope,nsoilmx_PEM,nsoilmx,TI_GCM_phys,TI_PEM)

!! now the temperature

DO l=1,nsoilmx
  tsoil_PEM(1,l,:)=tsoil_ave(1,1,l,:)
ENDDO
DO l=nsoilmx,nsoilmx_PEM
  tsoil_PEM(1,l,:)=tsoil_ave(1,1,nsoilmx,:)
ENDDO

     allocate(interp_coef(nslope))

  ig0 = 2
  DO j=2,jjm
    DO i = 1, iim
       DO l=1,nsoilmx
       tsoil_PEM(ig0,l,:)=tsoil_ave(i,j,l,:)             !Tsoil_PEM=tsoil_GCM_average
       ENDDO

       interp_coef(:)=( tsoil_ave(i,j,nsoilmx-1,:)-tsoil_ave(i,j,nsoilmx,:) ) / (mlayer_PEM(nsoilmx-1)-mlayer_PEM(nsoilmx))
!        write(*,*) 'interpcoeff=',interp_coef(:),mlayer_PEM
       DO l=nsoilmx+1,nsoilmx_PEM
       tsoil_PEM(ig0,l,:)=tsoil_ave(i,j,nsoilmx,:)+interp_coef(:)*(mlayer_PEM(l)-mlayer_PEM(nsoilmx))             !Tsoil_PEM=tsoil_GCM_average     
        ENDDO
       ig0= ig0 + 1
    ENDDO
  ENDDO


DO l=1,nsoilmx
  tsoil_PEM(ig0,l,:) = tsoil_ave(1,jjm+1,l,:)
ENDDO

interp_coef(:)=( tsoil_ave(i,j,nsoilmx-1,:)-tsoil_ave(i,j,nsoilmx,:) ) / (mlayer_PEM(nsoilmx-1)-mlayer_PEM(nsoilmx))

DO l=nsoilmx,nsoilmx_PEM
  tsoil_PEM(ig0,l,:) = tsoil_ave(1,jjm+1,nsoilmx,:)+interp_coef(:)*(mlayer_PEM(l)-mlayer_PEM(nsoilmx))
ENDDO

DO ig = 1,ngrid
DO iloop = 1, nsoilmx_PEM
DO islope = 1,nslope
   if (isnan(tsoil_PEM(ig,iloop,islope))) write(*,*) ig,iloop,islope
ENDDO
ENDDO
ENDDO



DO ig = 1,ngrid
  DO iloop = nsoilmx+1,nsoilmx_PEM
    DO islope = 1,nslope
      kcond = (TI_PEM(ig,iloop,islope)*TI_PEM(ig,iloop,islope))/volcapa
      tsoil_PEM(ig,iloop,islope) = tsoil_PEM(ig,nsoilmx,islope) + fluxgeo/kcond*(mlayer_PEM(iloop)-mlayer_PEM(nsoilmx))
    ENDDO
  ENDDO
ENDDO

do iloop = 1,nsoilmx_PEM
write(*,*) 'Tsoil1 = ',tsoil_PEM(10,iloop,1)
enddo

do iloop = 1,nsoilmx_PEM
write(*,*) 'Tsoil2 = ',tsoil_PEM(10,iloop,2)
enddo

do iloop = 1,nsoilmx_PEM
write(*,*) 'Tsoil3 = ',tsoil_PEM(10,iloop,3)
enddo

do iloop = 1,nsoilmx_PEM
write(*,*) 'Tsoil4 = ',tsoil_PEM(10,iloop,4)
enddo

do iloop = 1,nsoilmx_PEM
write(*,*) 'Tsoil5 = ',tsoil_PEM(10,iloop,5)
enddo

do iloop = 1,nsoilmx_PEM
write(*,*) 'Tsoil6 = ',tsoil_PEM(10,iloop,6)
enddo

do iloop = 1,nsoilmx_PEM
write(*,*) 'Tsoil7 = ',tsoil_PEM(10,iloop,7)
enddo


do iloop = 1,nsoilmx_PEM
write(*,*) 'Tsoil1 = ',tsoil_PEM(50,iloop,1)
enddo

do iloop = 1,nsoilmx_PEM
write(*,*) 'Tsoil2 = ',tsoil_PEM(50,iloop,2)
enddo

do iloop = 1,nsoilmx_PEM
write(*,*) 'Tsoil3 = ',tsoil_PEM(50,iloop,3)
enddo

do iloop = 1,nsoilmx_PEM
write(*,*) 'Tsoil4 = ',tsoil_PEM(50,iloop,4)
enddo

do iloop = 1,nsoilmx_PEM
write(*,*) 'Tsoil5 = ',tsoil_PEM(50,iloop,5)
enddo

do iloop = 1,nsoilmx_PEM
write(*,*) 'Tsoil6 = ',tsoil_PEM(50,iloop,6)
enddo

do iloop = 1,nsoilmx_PEM
write(*,*) 'Tsoil7 = ',tsoil_PEM(50,iloop,7)
enddo



do iloop = 1,nsoilmx_PEM
write(*,*) 'TI_PEM1 = ',TI_PEM(10,iloop,1)
enddo

do iloop = 1,nsoilmx_PEM
write(*,*) 'TI_PEM2 = ',TI_PEM(10,iloop,2)
enddo

do iloop = 1,nsoilmx_PEM
write(*,*) 'TI_PEM3 = ',TI_PEM(10,iloop,3)
enddo

do iloop = 1,nsoilmx_PEM
write(*,*) 'TI_PEM4 = ',TI_PEM(10,iloop,4)
enddo

do iloop = 1,nsoilmx_PEM
write(*,*) 'TI_PEM5 = ',TI_PEM(10,iloop,5)
enddo

do iloop = 1,nsoilmx_PEM
write(*,*) 'TI_PEM6 = ',TI_PEM(10,iloop,6)
enddo

do iloop = 1,nsoilmx_PEM
write(*,*) 'TI_PEM7 = ',TI_PEM(10,iloop,7)
enddo


do iloop = 1,nsoilmx_PEM
write(*,*) 'TI_PEM1 = ',TI_PEM(50,iloop,1)
enddo

do iloop = 1,nsoilmx_PEM
write(*,*) 'TI_PEM2 = ',TI_PEM(50,iloop,2)
enddo

do iloop = 1,nsoilmx_PEM
write(*,*) 'TI_PEM3 = ',TI_PEM(50,iloop,3)
enddo

do iloop = 1,nsoilmx_PEM
write(*,*) 'TI_PEM4 = ',TI_PEM(50,iloop,4)
enddo

do iloop = 1,nsoilmx_PEM
write(*,*) 'TI_PEM5 = ',TI_PEM(50,iloop,5)
enddo

do iloop = 1,nsoilmx_PEM
write(*,*) 'TI_PEM6 = ',TI_PEM(50,iloop,6)
enddo

do iloop = 1,nsoilmx_PEM
write(*,*) 'TI_PEM7 = ',TI_PEM(50,iloop,7)
enddo




tsoil_ave_phys_inst(:,:,:) = tsoil_PEM(:,1:nsoilmx,:)


  tsurf_ave_phys(1,:)=tsurf_ave(1,1,:)

  ig0 = 2
  DO j=2,jjm
    DO i = 1, iim
       tsurf_ave_phys(ig0,:)  =tsurf_ave(i,j,:)
       ig0= ig0 + 1
    ENDDO
  ENDDO

  tsurf_ave_phys(ig0,:) = tsurf_ave(1,jjm+1,:)
  tsurf_ave_phys_inst(:,:) = tsurf_ave_phys(:,:)


DO ig = 1,ngrid

DO islope = 1,nslope
   if (isnan(tsurf_ave_phys(ig,islope))) write(*,*) ig,iloop,islope
ENDDO

ENDDO

write(*,*) "ok test nan"



!     allocate(h2o_ice_first_last_day(iim+1,jjm+1,2))
!  call pemetat0_Y1 ("concat_year.nc", h2o_ice_first_last_day,iim,jjm)

!------------------------

! II  Run
!    II_a Compute tendencies

!------------------------

!---------------------------- RUN ---------------------




!----- Compute tendencies


     allocate(tendencies_h2o_ice(iim+1,jjm+1))
     allocate(tendencies_h2o_ice_phys(ngrid))

     allocate(tendencies_co2_ice(iim+1,jjm+1))
     allocate(tendencies_co2_ice_phys(ngrid))

     allocate(tendencies_co2_ice_slope(iim+1,jjm+1,nslope))
     allocate(tendencies_co2_ice_phys_slope(ngrid,nslope))
     allocate(tendencies_co2_ice_phys_slope_ini(ngrid,nslope))


     allocate(tendencies_h2o_ice_slope(iim+1,jjm+1,nslope))
     allocate(tendencies_h2o_ice_phys_slope(ngrid,nslope))

!  Compute the tendencies of the evolution of water ice over the years

      call compute_tendencies(tendencies_h2o_ice,min_h2o_ice_s_1,&
             min_h2o_ice_s_2,iim,jjm,ngrid,tendencies_h2o_ice_phys)

      call compute_tendencies(tendencies_co2_ice,min_co2_ice_s_1,&
             min_co2_ice_s_2,iim,jjm,ngrid,tendencies_co2_ice_phys)

      call compute_tendencies_slope(tendencies_co2_ice_slope,min_co2_ice_slope_1,&
             min_co2_ice_slope_2,iim,jjm,ngrid,tendencies_co2_ice_phys_slope,nslope)


      tendencies_co2_ice_phys_slope_ini(:,:)=tendencies_co2_ice_phys_slope(:,:)

      call compute_tendencies_slope(tendencies_h2o_ice_slope,min_h2o_ice_slope_1,&
             min_h2o_ice_slope_2,iim,jjm,ngrid,tendencies_h2o_ice_phys_slope,nslope)

!  Compute the tendencies of the evolution of water ice over one year

!      call compute_tendencies_Y1(tendencies_h2o_ice,min_h2o_ice_s_1,&
!             h2o_ice_first_last_day,iim,jjm,ngrid,tendencies_h2o_ice_phys)

     

!------------------------

! II  Run
!    II_a Compute tendencies
!    II_b Save initial situation

!------------------------

!----- Save initial situation

     allocate(initial_h2o_ice(ngrid))
!     allocate(initial_co2_ice(n_band_lat))
!     do j=1,n_band_lat
!        initial_co2_ice(j)=0.
!     enddo

     allocate(initial_co2_ice(ngrid))
!     allocate(q_co2_GCM(ngrid))
!     allocate(vmr_co2_gcm_phys(ngrid,timelen))

     allocate(initial_co2_ice_sublim_slope(ngrid,nslope))
     allocate(initial_co2_ice_slope(ngrid,nslope))
     allocate(initial_h2o_ice_slope(ngrid,nslope))
     year_iter=0

! We save the places where water ice is sublimating
  do i=1,ngrid
      if (tendencies_h2o_ice_phys(i).LT.0) then
         initial_h2o_ice(i)=1.
      else
         initial_h2o_ice(i)=0.         
      endif
    do islope=1,nslope

      if (tendencies_co2_ice_phys_slope(i,islope).LT.0) then
         initial_co2_ice_sublim_slope(i,islope)=1.
      else
         initial_co2_ice_sublim_slope(i,islope)=0.         
      endif
      if (co2ice_slope(i,islope).GT.0) then
         initial_co2_ice_slope(i,islope)=1.
      else
         initial_co2_ice_slope(i,islope)=0.         
      endif
      if (tendencies_h2o_ice_phys_slope(i,islope).LT.0) then
         initial_h2o_ice_slope(i,islope)=1.
      else
         initial_h2o_ice_slope(i,islope)=0.         
      endif
!      j=floor((latitude(i)+(pi/2))/(pi)*n_band_lat)+1
!      print *, "latitudelatitude",latitude(i)
!      print *, "jjjj", j
!      if(j.GT.n_band_lat) then
!          j=n_band_lat
!      endif
!      initial_co2_ice(j)=initial_co2_ice(j)+co2ice(i)*cell_area(i)
!!      q_co2_GCM(i)=q(i,1,1)
!!      ps_GCM(i)=ps(i)
!c             Mean air molecular mass = 1/(q(ico2)/m_co2 + (1-q(ico2))/m_noco2)
!!      mmean=1/(A*q_co2_GCM(i) +B)
!!      vmr_co2_GCM(i) = q_co2_GCM(i)*mmean/m_co2
    enddo
  enddo


     allocate(T_cond_GCM(iim+1,jjm+1,timelen))

  DO i=1,iim+1
    DO j = 1, jjm+1
      DO t = 1, timelen
        beta=3182.48
        alpha=23.3494
        T_cond_GCM(i,j,t)=beta/(alpha-log(vmr_co2_gcm(i,j,t)*ps_GCM_2(i,j,t)/100))
      ENDDO
    ENDDO
  ENDDO

! We compute the surface of water ice sublimating
  ini_surf=0.
  ini_surf_co2=0.
  ini_surf_h2o=0.
  Total_surface=0.
  do i=1,ngrid
    if (initial_h2o_ice(i).GT.0.5) then
       ini_surf=ini_surf+cell_area(i)
    endif
    do islope=1,nslope
      if (initial_co2_ice_sublim_slope(i,islope).GT.0.5) then
         ini_surf_co2=ini_surf_co2+cell_area(i)*subslope_dist(i,islope)
      endif
      if (initial_h2o_ice_slope(i,islope).GT.0.5) then
         ini_surf_h2o=ini_surf_h2o+cell_area(i)*subslope_dist(i,islope)
      endif
    enddo
    Total_surface=Total_surface+cell_area(i)
  enddo

     print *, "ini_surf_co2=", ini_surf_co2
     print *, "ini_surf=", ini_surf
     print *, "ini_surf_h2o=", ini_surf_h2o

!------------------------

! II  Run
!    II_a Compute tendencies
!    II_b Save initial situation
!    II_c Time loop

!------------------------

     allocate(local_old_press(ngrid))
     allocate(local_new_press(ngrid))
     allocate(ps_phys(ngrid))
     allocate(zplev_new(ngrid,nlayer+1))
     allocate(zplev_old(ngrid,nlayer+1))

  ig0 = iim+1
  ps_phys(1)=ps(ig0)
  ig0=ig0+1

  DO i = 2, ngrid-1
     if(modulo(ig0,iim+1).eq.0 .and. i.gt.3) then
       ig0=ig0+1
     endif
     ps_phys(i) = ps(ig0)
     ig0= ig0 + 1
  ENDDO

    ps_phys(ngrid)=ps(ig0)

     global_ave_press_old=0.
     do i=1,ngrid
       global_ave_press_old=global_ave_press_old+cell_area(i)*ps_phys(i)/Total_surface
     enddo

     global_ave_press_GCM=global_ave_press_old
     print *, "global_ave_press_old", global_ave_press_old
        DO l=1,nlayer+1
         DO ig=1,ngrid
          zplev_old(ig,l) = ap(l)  + bp(l)*ps(ig)
         ENDDO
        ENDDO
!        zplev_old(:,nlayer+1) = 0.


!----- Time loop

     allocate(flag_co2flow(ngrid,nslope))
     allocate(flag_co2flow_mesh(ngrid))

     firstcall=.TRUE.



!     do while (year_iter.LT.nyear)
     do while (.true.)

     global_ave_press_new=global_ave_press_old
     do i=1,ngrid
       do islope=1,nslope
           global_ave_press_new=global_ave_press_new-g*cell_area(i)*tendencies_co2_ice_phys_slope(i,islope)*subslope_dist(ig,islope)/cos(pi*def_slope_mean(islope)/180.)/Total_surface
       enddo
     enddo


     print *, "global_ave_press_new", global_ave_press_new
     print *, "pb here ?"
     do i=1,ngrid
       ps_phys(i)=ps_phys(i)*global_ave_press_new/global_ave_press_old
     enddo
      print *, "0.1 is okay"
     do i=1,ip1jmp1
!       local_new_press(i)=local_new_press(i)*global_ave_press_old/global_ave_press_new
       ps(i)=ps(i)*global_ave_press_new/global_ave_press_old
     enddo
      print *, "0.2 is okay"
        DO l=1,nlayer+1
         DO ig=1,ngrid
          zplev_new(ig,l) = ap(l)  + bp(l)*ps(ig)
         ENDDO
        ENDDO
!        zplev_new(:,nlayer+1) = 0.

      print *, "1 is okay"

        DO nnq=1,nqtot
        if (nnq.NE.1) then
          DO l=1,llm-1
            DO ig=1,ngrid
               q(ig,l,nnq)=q(ig,l,nnq)*(zplev_old(ig,l)-zplev_old(ig,l+1))/(zplev_new(ig,l)-zplev_new(ig,l+1))
            ENDDO
            q(:,llm,nnq)=q(:,llm-1,nnq)
          ENDDO
        else
!          tot_co2_atm=0.
!          tot_var_co2_atm=0.
!          DO l=1,llm
!            DO i=1,ngrid
!               tot_co2_atm=tot_co2_atm+q(i,l,nnq)*cell_area(i)
!               tot_var_co2_atm=tot_var_co2_atm+cell_area(i)*tendencies_co2_ice_phys(i)
!            ENDDO
!          ENDDO
!          DO l=1,llm
!            DO i=1,ngrid
!               q(i,l,nnq)=q(i,l,nnq)*(1-tot_var_co2_atm/tot_co2_atm)
!            ENDDO
!          ENDDO
          DO l=1,llm-1
            DO ig=1,ngrid
!                q(i,l,nnq)=q(i,l,nnq)*(global_ave_press_old/global_ave_press_new)+(global_ave_press_new-global_ave_press_old)/global_ave_press_new
                q(ig,l,nnq)=q(ig,l,nnq)*(zplev_old(ig,l)-zplev_old(ig,l+1))/(zplev_new(ig,l)-zplev_new(ig,l+1))  &
                                + (   (zplev_new(ig,l)-zplev_new(ig,l+1))  -       &
                                      (zplev_old(ig,l)-zplev_old(ig,l+1))     )  / &
                                      (zplev_new(ig,l)-zplev_new(ig,l+1))
            ENDDO
           q(:,llm,nnq)=q(:,llm-1,nnq)
          ENDDO
        endif
        ENDDO


                print *, "2 is okay"  
  l=1
  DO ig=1,ngrid
      DO t = 1, timelen
         q_h2o_PEM_phys(ig,t)=q_h2o_PEM_phys(ig,t)*(zplev_old(ig,l)-zplev_old(ig,l+1))/(zplev_new(ig,l)-zplev_new(ig,l+1))
     enddo
  enddo



  DO ig=1,ngrid
      DO t = 1, timelen
         q_co2_PEM_phys(ig,t)=q_co2_PEM_phys(ig,t)*(zplev_old(ig,l)-zplev_old(ig,l+1))/(zplev_new(ig,l)-zplev_new(ig,l+1))  &
                                + (   (zplev_new(ig,l)-zplev_new(ig,l+1))  -       &
                                      (zplev_old(ig,l)-zplev_old(ig,l+1))     )  / &
                                      (zplev_new(ig,l)-zplev_new(ig,l+1))
         if (q_co2_PEM_phys(ig,t).LT.0) then
              q_co2_PEM_phys(ig,t)=1E-10
         elseif (q_co2_PEM_phys(ig,t).GT.1) then
              q_co2_PEM_phys(ig,t)=1.
         endif
         if (q_co2_PEM_phys(ig,t).LT.0) then
              q_co2_PEM_phys(ig,t)=1E-30
         elseif (q_co2_PEM_phys(ig,t).GT.1) then
              q_co2_PEM_phys(ig,t)=1.
         endif
         mmean=1/(A*q_co2_PEM_phys(ig,t) +B)
         vmr_co2_pem_phys(ig,t) = q_co2_PEM_phys(ig,t)*mmean/m_co2
      ENDDO
  ENDDO
   

  DO ig = 1,ngrid
    ave=0.
    DO t=1,timelen
               ave=ave+q_co2_PEM_phys(ig,t)
    ENDDO
    q_co2_PEM_phys_ave=ave/timelen
    ave = 0.
  ENDDO

      print *, "3 is okay"
!We apply the tendency

!     call evol_h2o_ice_s(qsurf(:,7),tendencies_h2o_ice_phys,iim,jjm,ngrid,cell_area,STOPPING_1_water)
!TO DOOOOOOOOO
     call evol_h2o_ice_s_slope(qsurf_slope(:,6,:),tendencies_h2o_ice_phys_slope,iim,jjm,ngrid,cell_area,STOPPING_1_water,nslope)
!     call evol_co2_ice_s(qsurf(:,2),tendencies_co2_ice_phys,iim,jjm,ngrid,cell_area,STOPPING_1)
!     call evol_co2_ice_s(co2ice,tendencies_co2_ice_phys,iim,jjm,ngrid,cell_area,STOPPING_1)
      print *, "4 is okay"

     call evol_co2_ice_s_slope(co2ice_slope,tendencies_co2_ice_phys_slope,iim,jjm,ngrid,cell_area,STOPPING_1_co2,nslope)
!     call evol_h2o_ice_s(qsurf(:,5),tendencies_h2o_ice_phys,iim,jjm,ngrid)
      print *, "5 is okay"
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! FLOW     

!      print *, "2"
!      print *, "ngrid", ngrid
!      print *, "nslope", nslope
!      print *, "iflat", iflat

!       DO ig = 1,ngrid
!        DO islope = 1,nslope
!          IF(islope.ne.iflat) THEN ! ice can be infinite on flat ground
!! First: check that CO2 ice must flow (excess of ice on the slope), ice can accumulate on flat ground
!            IF(co2ice_slope(ig,islope).ge.rho_co2*co2_hmax * &
!                  cos(pi*def_slope_mean(islope)/180.)) THEN 

!! Second: determine the flatest sloptsoil_slopees possible:
       
!                IF(islope.gt.iflat) THEN
!                  iaval=islope-1
!                ELSE 
!                 iaval=islope+1
!                ENDIF
!                do while ((iaval.ne.iflat).or.      &
!                    (subslope_dist(ig,iaval).ne.0))
!                    print *, "ig", ig
!                    print *, "islope", islope
!                    print *, "iaval", iaval
!                    print *, "iflat", iflat
!                    print *, "subslope_dist(ig,iaval)", subslope_dist(ig,iaval)
!                    print *, "----------------------------"
!                  IF(iaval.gt.iflat) THEN
!                     iaval=iaval-1
!                  ELSE 
!                     iaval=iaval+1
!                  ENDIF
!                enddo


!                co2ice_slope(ig,iaval) = co2ice_slope(ig,iaval) + &
!               (co2ice_slope(ig,islope) - rho_co2* co2_hmax *     &
!               cos(pi*def_slope_mean(islope)/180.)) *            &
!               subslope_dist(ig,islope)/subslope_dist(ig,iaval) * &
!               cos(pi*def_slope_mean(iaval)/180.) /              &
!               cos(pi*def_slope_mean(islope)/180.)

!                co2ice_slope(ig,islope)=rho_co2*co2_hmax *        &
!                  cos(pi*def_slope_mean(islope)/180.) 

!             flag_co2flow(ig,islope) = 1.
!             flag_co2flow_mesh(ig) = 1.
!            ENDIF ! co2ice > hmax
!          ENDIF ! iflat
!        ENDDO !islope 
!       ENDDO !ig

       DO ig = 1,ngrid
        DO islope = 1,nslope
          IF(islope.ne.iflat) THEN ! ice can be infinite on flat ground
! First: check that CO2 ice must flow (excess of ice on the slope), ice can accumulate on flat ground
            IF(co2ice_slope(ig,islope).ge.rho_co2*co2_hmax * &
                  cos(pi*def_slope_mean(islope)/180.)) THEN 

! Second: determine the flatest slopes possible:
       
                IF(islope.gt.iflat) THEN
                  iaval=islope-1
                ELSE 
                 iaval=islope+1
                ENDIF
                do while ((iaval.ne.iflat).and.  &
                    (subslope_dist(ig,iaval).eq.0))

                  IF(iaval.gt.iflat) THEN
                     iaval=iaval-1
                  ELSE 
                     iaval=iaval+1
                  ENDIF
             
                enddo


                co2ice_slope(ig,iaval) = co2ice_slope(ig,iaval) + &
               (co2ice_slope(ig,islope) - rho_co2* co2_hmax *     &
               cos(pi*def_slope_mean(islope)/180.)) *             &
               subslope_dist(ig,islope)/subslope_dist(ig,iaval) * &
               cos(pi*def_slope_mean(iaval)/180.) /               &
               cos(pi*def_slope_mean(islope)/180.)                

                co2ice_slope(ig,islope)=rho_co2*co2_hmax *        &
                  cos(pi*def_slope_mean(islope)/180.) 

             flag_co2flow(ig,islope) = 1.
             flag_co2flow_mesh(ig) = 1.
            ENDIF ! co2ice > hmax
          ENDIF ! iflattsoil_lope
        ENDDO !islope 
       ENDDO !ig

      print *, "6 is okay"
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! END OF FLOW

       DO ig = 1,ngrid
        DO islope = 1,nslope
           if(initial_co2_ice_slope(ig,islope).gt.0.5 .and. co2ice_slope(ig,islope).LT. 1E-5) THEN !co2ice disappeared, look for closest point without co2ice
              if(latitude_deg(ig).gt.0) then
                do ig_loop=ig,ngrid
                  DO islope_loop=islope,iflat,-1
                     if(initial_co2_ice_slope(ig_loop,islope_loop).lt.0.5 .and. co2ice_slope(ig_loop,islope_loop).LT. 1E-5) then
                            tsurf_ave_phys(ig,islope)=tsurf_ave_phys(ig_loop,islope_loop)
                            bool_sublim=1
                            exit
                     endif
                  enddo
                  if (bool_sublim.eq.1) then
                   exit
                  endif
                enddo
              else
               do ig_loop=ig,1,-1
                  DO islope_loop=islope,iflat
                     if(initial_co2_ice_slope(ig_loop,islope_loop).lt.0.5 .and. co2ice_slope(ig_loop,islope_loop).LT. 1E-5) then
                            tsurf_ave_phys(ig,islope)=tsurf_ave_phys(ig_loop,islope_loop)
                            bool_sublim=1
                            exit
                     endif
                  enddo
                  if (bool_sublim.eq.1) then
                   exit
                  endif
                enddo
              endif
              initial_co2_ice_slope(ig,islope)=0
           elseif(initial_co2_ice_slope(ig,islope).gt.0.5 .and. co2ice_slope(ig,islope).GT. 1E-5) THEN !Put tsurf as tcond co2
             ave=0.
             do t=1,timelen
               ave=ave+vmr_co2_pem_phys(ig,t)*ps_phys(ig)/100
             enddo
             tsurf_ave_phys(ig,islope)=beta/(alpha-log(ave/timelen))
           endif
        enddo
       enddo

      print *, "7 is okay"


!Let's update the soil_temperature & soil thermal inertia

allocate(TI_locslope(ngrid,nsoilmx_PEM))
allocate(Tsoil_locslope(ngrid,nsoilmx_PEM))
allocate(alph_locslope(ngrid,nsoilmx_PEM))
allocate(beta_locslope(ngrid,nsoilmx_PEM))
allocate(Tsurf_locslope(ngrid))

      print *, "8 is okay"
if (firstcall) then
do islope = 1,nslope
     TI_locslope(:,:) = TI_PEM(:,:,islope)
     Tsoil_locslope(:,:) = tsoil_PEM(:,:,islope)
     Tsurf_locslope(:) = tsurf_ave_phys(:,islope) 
     alph_locslope(:,:) = alph_PEM(:,:,islope)
     beta_locslope(:,:) = beta_PEM(:,:,islope)
!     call soil_pem(ngrid,nsoilmx_PEM,.true.,TI_locslope,timestep,Tsurf_locslope,Tsoil_locslope, mlayer_PEM,layer_PEM,alph_locslope,beta_locslope)  
     TI_PEM(:,:,islope) = TI_locslope 
     tsoil_PEM(:,:,islope) = Tsoil_locslope
     tsurf_ave_phys(:,islope)  =  Tsurf_locslope
     alph_PEM(:,:,islope) = alph_locslope(:,:) 
     beta_PEM(:,:,islope) = beta_locslope(:,:) 
 enddo
firstcall = .false.
endif



do islope = 1,nslope
     TI_locslope(:,:) = TI_PEM(:,:,islope)
     Tsoil_locslope(:,:) = tsoil_PEM(:,:,islope)
     Tsurf_locslope(:) = tsurf_ave_phys(:,islope) 
     alph_locslope(:,:) = alph_PEM(:,:,islope)
     beta_locslope(:,:) = beta_PEM(:,:,islope)
!     call soil_pem(ngrid,nsoilmx_PEM,.false.,TI_locslope,timestep,Tsurf_locslope,Tsoil_locslope, mlayer_PEM,layer_PEM,alph_locslope,beta_locslope)  
     TI_PEM(:,:,islope) = TI_locslope 
     tsoil_PEM(:,:,islope) = Tsoil_locslope
     tsurf_ave_phys(:,islope)  =  Tsurf_locslope
     alph_PEM(:,:,islope) = alph_locslope(:,:) 
     beta_PEM(:,:,islope) = beta_locslope(:,:) 
 enddo

      print *, "9 is okay"



!! Let's compute the average of qh2o_vap
   do  ig=1,ngrid
      q_h2o_PEM_phys_ave(ig) = 0.
      do  t = 1, timelen
       q_h2o_PEM_phys_ave(ig) = q_h2o_PEM_phys_ave(ig) + q_h2o_PEM_phys(ig,t)/timelen
     enddo
  enddo

      print *, "10 is okay"

 call update_soil(ngrid,nslope,nsoilmx_PEM,tendencies_h2o_ice_phys_slope,tendencies_co2_ice_phys_slope,co2ice_slope,qsurf_slope(:,6,:),ps_phys, &
          q_h2o_PEM_phys_ave,tsoil_PEM,cell_area,subslope_dist,ice_depth,TI_PEM)

      print *, "11 is okay"
   
!    call recomp_tend_co2(tendencies_co2_ice_phys,q(:,1,1),ps,q_co2_GCM,ps_GCM,ngrid)
!    call recomp_tend_co2(tendencies_co2_ice_phys,vmr_co2_gcm_phys,ps_GCM_2,global_ave_press_GCM,global_ave_press_new,timelen,ngrid)
    call recomp_tend_co2_slope(tendencies_co2_ice_phys_slope,tendencies_co2_ice_phys_slope_ini,vmr_co2_gcm_phys,vmr_co2_pem_phys,ps_GCM_2,&     
                               global_ave_press_GCM,global_ave_press_new,timelen,ngrid,nslope)

          print *, "12 is okay"
      year_iter=year_iter+dt_pem

!We check with we should stop



!      call criterion_ice_stop(cell_area,initial_h2o_ice,qsurf(:,5),STOPPING,ngrid)
!      call criterion_ice_stop(cell_area,initial_h2o_ice,qsurf(:,3),STOPPING,ngrid)physdem
!!!!!!! CORRECT      call criterion_ice_stop_water(cell_area,ini_surf,qsurf(:,7),STOPPING_water,ngrid,initial_h2o_ice)
!TO DOOOOOOOOO_slope
      call criterion_ice_stop_water_slope(cell_area,ini_surf_h2o,qsurf_slope(:,6,:),STOPPING_water,ngrid,initial_h2o_ice)
!      call criterion_ice_stop(cell_area,ini_surf_co2,qsurf(:,2),STOPPING,ngrid,initial_co2_ice)
!      call criterion_ice_stop(cell_area,ini_surf_co2,co2ice,STOPPING,ngrid,initial_co2_ice,global_ave_press_GCM,global_ave_press_new)
          print *, "13 is okay"
      call criterion_ice_stop_slope(cell_area,ini_surf_co2,co2ice_slope,STOPPING_co2,ngrid,initial_co2_ice_sublim_slope,global_ave_press_GCM,global_ave_press_new,nslope)

 print *, "criterion ok"
!      call criterion_co2_ice_stop(cell_area,initial_co2_ice,co2ice,STOPPING,ngrid,latitude,n_band_lat)

      print *, "STOPPING_water", STOPPING_water, "STOPPING1_water", STOPPING_1_water
      print *, "STOPPING_co2", STOPPING_co2, "STOPPING1_co2", STOPPING_1_co2
      print *, "year_iter=", year_iter
      if (STOPPING_water .or. STOPPING_1_water .or. STOPPING_co2 .or. STOPPING_1_co2)  then
        exit
      endif

      global_ave_press_old=global_ave_press_new
      zplev_old=zplev_new


      enddo


      DO ig = 1,ngrid
          co2ice(ig) = 0.
             DO islope = 1,nslope
                 co2ice(ig) = co2ice(ig) + co2ice_slope(ig,islope) &
                            * subslope_dist(ig,islope) /          &
                           cos(pi*def_slope_mean(islope)/180.)
             ENDDO
      ENDDO ! of DO ig=1,ngrid

      DO ig = 1,ngrid !!!!! TO DOC GHANGE IF QSURF ?ü!?!?!?
          qsurf(ig,6) = 0.
             DO islope = 1,nslope
                 qsurf(ig,6) = qsurf(ig,6) + qsurf_slope(ig,6,islope) &
                            * subslope_dist(ig,islope) /          &
                           cos(pi*def_slope_mean(islope)/180.)
             ENDDO
      ENDDO ! of DO ig=1,ngrid




!! Last soil update
 call update_soil(ngrid,nslope,nsoilmx_PEM,tendencies_h2o_ice_phys_slope,tendencies_co2_ice_phys_slope,co2ice_slope,qsurf_slope(:,6,:),ps_phys, &
          q_h2o_PEM_phys_ave,tsoil_PEM,tsurf_ave_phys,cell_area,subslope_dist,ice_depth,TI_PEM)


   call interpolate_TIPEM_TIGCM(ngrid,nslope,nsoilmx_PEM,nsoilmx,TI_PEM,TI_GCM_phys)


!! NOW that temperature averaged have been fixed, let's update the instantanous
!temperature

do ig = 1,ngrid
 do islope = 1,nslope
!   tsurf_slope(ig,islope) =  tsurf_slope(ig,islope) - (tsurf_ave_phys_inst(ig,islope)-tsurf_ave_phys(ig,islope))
tsurf_slope(ig,islope) = tsurf_ave_phys(ig,islope) 
  do iloop = 1,nsoilmx
!     tsoil_slope(ig,iloop,islope) = tsoil_slope(ig,iloop,islope) - (tsoil_ave_phys_inst(ig,iloop,islope) - tsoil_PEM(ig,iloop,islope))
tsoil_slope(ig,iloop,islope) = tsoil_PEM(ig,iloop,islope) 

  enddo  
 enddo
enddo



!------------------------

! III Output
!    III_a Write startfi.nc

!------------------------

!----------------------------WRITE restart.nc ---------------------

      ptimestep=iphysiq*daysec/REAL(day_step)/nsplit_phys
      pday=day_ini
      ztime_fin=0.

         print *, "RVip1jmp1", ip1jmp1
     allocate(p(ip1jmp1,nlayer+1))
         print *, "ngrid", ngrid
         print *, "nlayer", nlayer
          CALL pression (ip1jmp1,ap,bp,ps,p)
     print *, "masse 1", masse(1,:)
          CALL massdair(p,masse)
     print *, "masse 2", masse(1,:)

      CALL dynredem0("restart_evol.nc", day_ini, phis)

      CALL dynredem1("restart_evol.nc",   &
      		time_0,vcov,ucov,teta,q,masse,ps)


!----------------------------WRITE restartfi.nc ---------------------

!      call physdem0("restartfi_evol.nc",longitude,latitude,     &
!            nsoilmx,ngrid,nlayer,nq,                       &
!            ptimestep,pday,time_phys,cell_area,            &
!            albedodat,inertiedat,zmea,zstd,zsig,zgam,zthe, &
!            hmons,summit,base)

             call physdem0("restartfi_evol.nc",longitude,latitude, &
                        nsoilmx,ngrid,nlayer,nq,              &
                        ptimestep,pday,0.,cell_area,          &
                        albedodat,inertiedat,zmea,zstd,zsig,zgam,zthe, &
                        hmons,summit,base,nslope,def_slope,   &
                        subslope_dist)

!      call physdem1("restartfi_evol.nc",nsoilmx,ngrid,nlayer,nq,&
!     	    ptimestep,ztime_fin,                           &
!            tsurf,tsoil,co2ice,albedo,emis,                &
!            q2,qsurf,tauscaling,totcloudfrac,wstar,        &
!            mem_Mccn_co2,mem_Nccn_co2,mem_Mh2o_co2,watercap)

          call physdem1("restartfi_evol.nc",nsoilmx,ngrid,nlayer,nq,  &
                     ptimestep,ztime_fin,                        &
                     tsurf,tsoil,co2ice,albedo,emis,             &
                     q2,qsurf,tauscaling,totcloudfrac,wstar,     &
                     watercap,nslope,co2ice_slope,               &
                     tsurf_slope,tsoil_slope, albedo_slope,      &
                     emiss_slope,qsurf_slope,watercap_slope, TI_GCM_phys)

      do islope = 1,nslope
       write(str2(1:2),'(i2.2)') islope
       call WRITEDIAGFI(ngrid,'ice_depth'//str2,'Ice Depth', &
                  'm',2,ice_depth(:,islope))
      enddo


      print *, "RV : So far so good"




!      deallocate(longitude)
!      deallocate(latitude)
!      deallocate(cell_area)


END PROGRAM pem

