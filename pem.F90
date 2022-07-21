
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
  REAL, dimension(:,:),allocatable :: ps_phys_timeseries !(ngrid x timelen) ! pression  au sol instantannées
  REAL masse(ip1jmp1,llm)                ! masse d'air
  REAL phis(ip1jmp1)                     ! geopotentiel au sol
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
      REAL , dimension(:,:,:), allocatable ::  zplev_new_timeseries  ! same but with the time series
      REAL , dimension(:,:,:), allocatable :: zplev_old_timeseries! same but with the time series




      REAL :: tot_co2_atm,tot_var_co2_atm


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


      REAL , dimension(:,:,:,:), allocatable :: co2_ice_GCM_slope        ! LON x LATX NSLOPE x Times field : co2 ice given by the GCM
      REAL , dimension(:,:,:), allocatable :: co2_ice_GCM_phys_slope        ! LON x LATX NSLOPE x Times field : co2 ice given by the GCM


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
     REAL, ALLOCATABLE :: tsurf_GCM_timeseries(:,:,:,:) !LONX LAT x SLOPE XTULES field : NOn averaged Surf Temperature [K]
     REAL, ALLOCATABLE :: tsurf_phys_GCM_timeseries_update(:,:,:) !IG x SLOPE XTULES field : NOn averaged Surf Temperature [K]
     REAL, ALLOCATABLE :: tsurf_phys_GCM_timeseries(:,:,:) !IG x SLOPE XTULES field : NOn averaged Surf Temperature [K]
     REAL, ALLOCATABLE :: tsoil_ave_phys_inst(:,:,:)
     REAL, ALLOCATABLE :: tsoil_phys_PEM_timeseries(:,:,:,:) !IG x SLOPE XTULES field : NOn averaged Soil Temperature [K]
     REAL, ALLOCATABLE :: tsoil_GCM_timeseries(:,:,:,:,:) !IG x SLOPE XTULES field : NOn averaged Soil Temperature [K]
     REAL, ALLOCATABLE :: tsoil_phys_PEM_timeseries_update(:,:,:,:) !IG x SLOPE XTULES field : NOn averaged Soil Temperature [K]

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
     REAL,ALLOCATABLE  :: Tsoilave_before_saved(:,:,:)
     REAL,ALLOCATABLE  :: Tsurfave_before_saved(:,:)
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
     CALL dynetat0(FILE_NAME_start,vcov,ucov, &
                    teta,q,masse,ps,phis, time_0)
     CALL iniconst !new
     CALL inigeom
     allocate(ap(nlayer+1))
     allocate(bp(nlayer+1))
     status =nf90_open(FILE_NAME_start, NF90_NOWRITE, ncid)
     status = nf90_inq_varid(ncid, "ap", apvarid)
     status = nf90_get_var(ncid, apvarid, ap)
     status = nf90_inq_varid(ncid, "bp", bpvarid)
     status = nf90_get_var(ncid, bpvarid, bp)
     status =nf90_close(ncid)
     
   



    CALL iniphysiq(iim,jjm,llm, &
          (jjm-1)*iim+2,comm_lmdz, &
          daysec,day_ini,dtphys/nsplit_phys, &
          rlatu,rlatv,rlonu,rlonv,aire,cu,cv,rad,g,r,cpp, &
          iflag_phys)


!----------------------------reading ---------------------


! First we read the initial state (starfi.nc)


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
!     Define some slope statistics

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

!------------------------
! I   Initialisation
!    I_a READ run.def 
!    I_b READ starfi_0.nc
!    I_c COMPUTE tendencies based on the GCM
!------------------------

! First we read the evolution of water and co2 ice (and the mass mixing ratio) over the first year of the GCM run, saving only the minimum value

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
     allocate(tsoil_GCM_timeseries(iim+1,jjm+1,nsoilmx,nslope,2676))
     allocate(tsoil_phys_PEM_timeseries(ngrid,nsoilmx_PEM,nslope,2676))
     allocate(Tsoilave_before_saved(ngrid,nsoilmx_PEM,nslope))
     allocate(tsoil_phys_PEM_timeseries_update(ngrid,nsoilmx_PEM,nslope,2676))
     allocate(tsurf_GCM_timeseries(iim+1,jjm+1,nslope,2676))
     allocate(tsurf_phys_GCM_timeseries(ngrid,nslope,2676))
     allocate(tsurf_phys_GCM_timeseries_update(ngrid,nslope,2676))
     allocate(co2_ice_GCM_phys_slope(ngrid,nslope,2676))
     allocate(co2_ice_GCM_slope(iim+1,jjm+1,nslope,2676))
     allocate(Tsurfave_before_saved(ngrid,nslope))


     call pemetat0 ("concat_year_one.nc", min_h2o_ice_s_1,min_co2_ice_s_1,iim,jjm,nlayer,vmr_co2_gcm,ps_GCM_1,timelen,min_co2_ice_slope_1,min_h2o_ice_slope_1,&   
                       nslope,tsurf_ave,tsoil_ave, tsurf_GCM_timeseries,tsoil_GCM_timeseries,TI_GCM_ave,q_co2_GCM,q_h2o_GCM,co2_ice_GCM_slope)

! Then we read the evolution of water and co2 ice (and the mass mixing ratio) over the second year of the GCM run, saving only the minimum value

     allocate(min_h2o_ice_s_2(iim+1,jjm+1))
     allocate(min_co2_ice_s_2(iim+1,jjm+1))
     allocate(min_co2_ice_slope_2(iim+1,jjm+1,nslope))
     allocate(min_h2o_ice_slope_2(iim+1,jjm+1,nslope))

     call pemetat0 ("concat_year_two.nc", min_h2o_ice_s_2,min_co2_ice_s_2,iim,jjm,nlayer,vmr_co2_gcm,ps_GCM_2,timelen,min_co2_ice_slope_2,min_h2o_ice_slope_2, &
                  nslope,tsurf_ave,tsoil_ave, tsurf_GCM_timeseries,tsoil_GCM_timeseries,TI_GCM_ave,q_co2_GCM,q_h2o_GCM,co2_ice_GCM_slope)




! The variables in the dynamic grid are transfered to the physical grid


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
       vmr_co2_gcm_phys(ig0,:)  = vmr_co2_gcm(i,j,:)
       q_h2o_GCM_phys(ig0,:)    = q_h2o_GCM(i,j,:)
       q_co2_GCM_phys(ig0,:)    = q_co2_GCM(i,j,:)
       ig0= ig0 + 1
      ENDDO
     ENDDO

    vmr_co2_gcm_phys(ig0,:) = vmr_co2_gcm(1,jjm+1,:)
    q_h2o_GCM_phys(ig0,:) = q_h2o_GCM(1,jjm+1,:)
    q_co2_GCM_phys(ig0,:) = q_co2_GCM(1,jjm+1,:)
    q_co2_PEM_phys(:,:)=  q_co2_GCM_phys(:,:)
    q_h2o_PEM_phys(:,:)=  q_h2o_GCM_phys(:,:)




!------------------------
! I   Initialisation
!    I_a READ run.def 
!    I_b READ starfi_0.nc
!    I_c COMPUTE tendencies based on the GCM
!    I_d Initialize the soil for the PEM (will be read in startfiPEM later)
!------------------------

!!!!!!!!!!!!!!!!!!!!! Initialisation for soil_PEM

      call end_comsoil_h_PEM
      call ini_comsoil_h_PEM(ngrid,nslope)
      timestep=year_day*daysec/year_step
      allocate(ice_depth(ngrid,nslope))
      allocate(TI_GCM_phys(ngrid,nsoilmx,nslope))

! The variables in the dynamic grid are transfered to the physical grid
      DO l=1,nsoilmx
         TI_GCM_phys(1,l,:)=TI_GCM_ave(1,1,l,:)
      ENDDO
      ig0 = 2
      DO j=2,jjm
        DO i = 1, iim
          DO l=1,nsoilmx
            TI_GCM_phys(ig0,l,:)=TI_GCM_ave(i,j,l,:)
          ENDDO
        ig0= ig0 + 1
        ENDDO
       ENDDO
      DO l=1,nsoilmx
       TI_GCM_phys(ig0,l,:) = TI_GCM_ave(1,jjm+1,l,:)
      ENDDO
! Thermal inertia are transfered from the GCM depth grid to the PEM depth grid
      call soil_settings_PEM(ngrid,nslope,nsoilmx_PEM,nsoilmx,TI_GCM_phys,TI_PEM)

!!!! now same but for the soil temperature. Soil temperature is initialized before reading in a start pem
!! First the average  
      DO l=1,nsoilmx
        tsoil_PEM(1,l,:)=tsoil_ave(1,1,l,:)
      ENDDO
      DO l=nsoilmx,nsoilmx_PEM
         tsoil_PEM(1,l,:)=tsoil_ave(1,1,nsoilmx,:)
      ENDDO

      ig0 = 2
      DO j=2,jjm
         DO i = 1, iim
           DO l=1,nsoilmx
             tsoil_PEM(ig0,l,:)=tsoil_ave(i,j,l,:)             !Tsoil_PEM=tsoil_GCM_average
            ENDDO
           
            ig0= ig0 + 1
         ENDDO
       ENDDO


      DO l=1,nsoilmx
        tsoil_PEM(ig0,l,:) = tsoil_ave(1,jjm+1,l,:)
      ENDDO


    DO ig = 1,ngrid
      DO iloop=nsoilmx+1,nsoilmx_PEM
         DO islope = 1,nslope
            kcond = (TI_PEM(ig,iloop,islope)*TI_PEM(ig,iloop,islope))/volcapa
            tsoil_PEM(ig,iloop,islope) = tsoil_PEM(ig,nsoilmx,islope) + fluxgeo/kcond*(mlayer_PEM(iloop)-mlayer_PEM(nsoilmx))
         ENDDO      
    ENDDO
    ENDDO

      Tsoilave_before_saved(:,:,:) = tsoil_PEM(:,:,:)

   
!! And now the time series
      DO l=1,nsoilmx
        tsoil_phys_PEM_timeseries(1,l,:,:)=tsoil_GCM_timeseries(1,1,l,:,:)
      ENDDO
      DO l=nsoilmx,nsoilmx_PEM
         tsoil_phys_PEM_timeseries(1,l,:,:)=tsoil_GCM_timeseries(1,1,nsoilmx,:,:)
      ENDDO
      
      ig0 = 2
      DO j=2,jjm
         DO i = 1, iim
           DO l=1,nsoilmx
             tsoil_phys_PEM_timeseries(ig0,l,:,:)=tsoil_GCM_timeseries(i,j,l,:,:)             !Tsoil_PEM=tsoil_GCM_average
            ENDDO
            ig0= ig0 + 1
         ENDDO
       ENDDO


      DO l=1,nsoilmx
       tsoil_phys_PEM_timeseries(ig0,l,:,:) = tsoil_GCM_timeseries(1,jjm+1,l,:,:)
      ENDDO


    DO ig = 1,ngrid
     DO islope = 1,nslope
      DO iloop=nsoilmx+1,nsoilmx_PEM
            kcond = (TI_PEM(ig,iloop,islope)*TI_PEM(ig,iloop,islope))/volcapa
            tsoil_phys_PEM_timeseries(ig,iloop,islope,:) = tsoil_phys_PEM_timeseries(ig,nsoilmx,islope,:) + fluxgeo/kcond*(mlayer_PEM(iloop)-mlayer_PEM(nsoilmx))
      ENDDO
     ENDDO
    ENDDO

    DO ig = 1,ngrid
      DO islope = 1,nslope
        DO iloop = 1,nsoilmx+1
           if(isnan(tsoil_PEM(ig,iloop,islope))) then 
         write(*,*) "failed nan construction", ig, iloop, islope
           stop
           endif
           do t = 1,timelen
              if(isnan(tsoil_phys_PEM_timeseries(ig,iloop,islope,t))) then


                write(*,*) "failed construction",ig,iloop,islope,t
                stop
              endif
            
           enddo

        ENDDO
      ENDDO
     ENDDO

      write(*,*) "construction ok, no nan" 






!! now the surface temperature

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




     Tsurfave_before_saved(:,:) = tsurf_ave_phys(:,:)

      tsurf_phys_GCM_timeseries(1,:,:)=tsurf_GCM_timeseries(1,1,:,:)
      ig0 = 2
      DO j=2,jjm
        DO i = 1, iim
          tsurf_phys_GCM_timeseries(ig0,:,:)  =tsurf_GCM_timeseries(i,j,:,:)
          ig0= ig0 + 1
        ENDDO
      ENDDO
      tsurf_phys_GCM_timeseries(ig0,:,:) = tsurf_GCM_timeseries(1,jjm+1,:,:)



! Now the CO2 ice 


      co2_ice_GCM_phys_slope(1,:,:)=co2_ice_GCM_slope(1,1,:,:)
      ig0 = 2
      DO j=2,jjm
        DO i = 1, iim
          co2_ice_GCM_phys_slope(ig0,:,:)  = co2_ice_GCM_slope(i,j,:,:)
          ig0= ig0 + 1
        ENDDO
      ENDDO
      co2_ice_GCM_phys_slope(ig0,:,:) = co2_ice_GCM_slope(1,jjm+1,:,:)
!---------------------------- END INITIALISATION ---------------------




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

     

!------------------------

! II  Run
!    II_a Compute tendencies
!    II_b Save initial situation

!------------------------

!----- Save initial situation

     allocate(initial_h2o_ice(ngrid))
     allocate(initial_co2_ice(ngrid))
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
     allocate(ps_phys_timeseries(ngrid,timelen))
     allocate(zplev_new(ngrid,nlayer+1))
     allocate(zplev_old(ngrid,nlayer+1))
     allocate(zplev_new_timeseries(ngrid,nlayer+1,timelen))
     allocate(zplev_old_timeseries(ngrid,nlayer+1,timelen))

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


! for the time series:





  ps_phys_timeseries(1,:)= ps_GCM_2(1,1,:)
  ig0=2
  DO j=2,jjm
    DO i = 1, iim
      ps_phys_timeseries(ig0,:) = ps_GCM_2(i,j,:)
     ig0= ig0 + 1
    ENDDO
  ENDDO
  ps_phys_timeseries(ig0,:) = ps_GCM_2(1,jjm+1,:)



     global_ave_press_old=0.
     do i=1,ngrid
       global_ave_press_old=global_ave_press_old+cell_area(i)*ps_phys(i)/Total_surface
     enddo

     global_ave_press_GCM=global_ave_press_old
     print *, "global_ave_press_old", global_ave_press_old
        DO l=1,nlayer+1
         DO ig=1,ngrid
          zplev_old(ig,l) = ap(l)  + bp(l)*ps_phys(ig)
          zplev_old_timeseries(ig,l,:) = ap(l)  + bp(l)*ps_phys_timeseries(ig,:)
         ENDDO
        ENDDO



!----- Time loop

     allocate(flag_co2flow(ngrid,nslope))
     allocate(flag_co2flow_mesh(ngrid))

     firstcall=.TRUE.
     do while (.true.)

     global_ave_press_new=global_ave_press_old
     do i=1,ngrid
       do islope=1,nslope
           global_ave_press_new=global_ave_press_new-g*cell_area(i)*tendencies_co2_ice_phys_slope(i,islope)*subslope_dist(ig,islope)/cos(pi*def_slope_mean(islope)/180.)/Total_surface
       enddo
     enddo


     do i=1,ip1jmp1

       ps(i)=ps(i)*global_ave_press_new/global_ave_press_old

     enddo
     do i = 1,ngrid
       ps_phys(i)=ps_phys(i)*global_ave_press_new/global_ave_press_old
       ps_phys_timeseries(i,:) = ps_phys_timeseries(i,:)*global_ave_press_new/global_ave_press_old
     enddo


        do l=1,nlayer+1
         do ig=1,ngrid
          zplev_new(ig,l) = ap(l)  + bp(l)*ps_phys(ig)
          zplev_new_timeseries(ig,l,:)  = ap(l)  + bp(l)*ps_phys_timeseries(ig,:)
         enddo
        enddo


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
          DO l=1,llm-1
            DO ig=1,ngrid
                q(ig,l,nnq)=q(ig,l,nnq)*(zplev_old(ig,l)-zplev_old(ig,l+1))/(zplev_new(ig,l)-zplev_new(ig,l+1))  &
                                + (   (zplev_new(ig,l)-zplev_new(ig,l+1))  -       &
                                      (zplev_old(ig,l)-zplev_old(ig,l+1))     )  / &
                                      (zplev_new(ig,l)-zplev_new(ig,l+1))
            ENDDO
           q(:,llm,nnq)=q(:,llm-1,nnq)
          ENDDO
        endif
        ENDDO

  l=1
  DO ig=1,ngrid
      DO t = 1, timelen
         q_h2o_PEM_phys(ig,t)=q_h2o_PEM_phys(ig,t)*(zplev_old_timeseries(ig,l,t)-zplev_old_timeseries(ig,l+1,t))/(zplev_new_timeseries(ig,l,t)-zplev_new_timeseries(ig,l+1,t))
     enddo
  enddo



  DO ig=1,ngrid
      DO t = 1, timelen
         q_co2_PEM_phys(ig,t)=q_co2_PEM_phys(ig,t)*(zplev_old_timeseries(ig,l,t)-zplev_old_timeseries(ig,l+1,t))/(zplev_new_timeseries(ig,l,t)-zplev_new_timeseries(ig,l+1,t))  &
                                + (   (zplev_new_timeseries(ig,l,t)-zplev_new_timeseries(ig,l+1,t))  -       &
                                      (zplev_old_timeseries(ig,l,t)-zplev_old_timeseries(ig,l+1,t))     )  / &
                                      (zplev_new_timeseries(ig,l,t)-zplev_new_timeseries(ig,l+1,t))
 
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


!We apply the tendency


     call evol_h2o_ice_s_slope(qsurf_slope(:,6,:),tendencies_h2o_ice_phys_slope,iim,jjm,ngrid,cell_area,STOPPING_1_water,nslope)
      print *, "4 is okay"

     call evol_co2_ice_s_slope(co2ice_slope,tendencies_co2_ice_phys_slope,iim,jjm,ngrid,cell_area,STOPPING_1_co2,nslope)
      print *, "5 is okay"
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! FLOW     

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


!-------------- Update Tsurf -------------

  
       DO ig = 1,ngrid
        DO islope = 1,nslope
           if(initial_co2_ice_slope(ig,islope).gt.0.5 .and. co2ice_slope(ig,islope).LT. 1E-5) THEN !co2ice disappeared, look for closest point without co2ice

              write(*,*) 'in if',ig,islope,initial_co2_ice_slope(ig,islope), co2ice_slope(ig,islope),tsurf_ave_phys(ig,islope)
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
               write(*,*) 'out if',ig,islope,ig_loop,islope_loop,tsurf_ave_phys(ig_loop,islope_loop),initial_co2_ice_slope(ig_loop,islope_loop),co2ice_slope(ig_loop,islope_loop)
              endif
              initial_co2_ice_slope(ig,islope)=0
           elseif( co2ice_slope(ig,islope).GT. 1E-5) THEN !Put tsurf as tcond co2
  
             ave=0.
             do t=1,timelen
              if(co2_ice_GCM_phys_slope(ig,islope,t).gt.1e-5) then
                 ave = ave + beta/(alpha-log(vmr_co2_pem_phys(ig,t)*ps_phys_timeseries(ig,t)/100.))
              else
                 ave = ave + tsurf_phys_GCM_timeseries(ig,islope,t)
              endif
             enddo
	     write(*,*) 'in elseif',ig,islope,initial_co2_ice_slope(ig,islope), co2ice_slope(ig,islope), tsurf_ave_phys(ig,islope),ave/timelen
             tsurf_ave_phys(ig,islope)=ave/timelen
           endif
        enddo
       enddo

      print *, "7 is okay"
      write(*,*) 'diff Tsurf after =', tsurf_ave_phys(:,:) -Tsurfave_before_saved(:,:)

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
     call soil_pem(ngrid,nsoilmx_PEM,.true.,TI_locslope,timestep,Tsurf_locslope,Tsoil_locslope, mlayer_PEM,layer_PEM,alph_locslope,beta_locslope)  
     TI_PEM(:,:,islope) = TI_locslope(:,:) 
     tsoil_PEM(:,:,islope) = Tsoil_locslope(:,:)
     tsurf_ave_phys(:,islope)  =  Tsurf_locslope(:)
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
     call soil_pem(ngrid,nsoilmx_PEM,.false.,TI_locslope,timestep,Tsurf_locslope,Tsoil_locslope, mlayer_PEM,layer_PEM,alph_locslope,beta_locslope)  
     TI_PEM(:,:,islope) = TI_locslope(:,:) 
     tsoil_PEM(:,:,islope) = Tsoil_locslope(:,:)
     tsurf_ave_phys(:,islope)  =  Tsurf_locslope(:)
     alph_PEM(:,:,islope) = alph_locslope(:,:) 
     beta_PEM(:,:,islope) = beta_locslope(:,:) 
 enddo

      print *, "9 is okay"


!!!!!! Now we recompute 
     do t = 1,timelen
     tsoil_phys_PEM_timeseries_update(:,:,:,t) = tsoil_phys_PEM_timeseries(:,:,:,t) +( tsoil_PEM(:,:,:) -Tsoilave_before_saved(:,:,:)) 
     tsurf_phys_GCM_timeseries_update(:,:,t) = tsurf_phys_GCM_timeseries(:,:,t) +( tsurf_ave_phys(:,:) -Tsurfave_before_saved(:,:)) 
     enddo

!     write(*,*) 'tsoil_update', tsoil_phys_PEM_timeseries_update(:,:,:,timelen)


!     write(*,*) 'tsoil_avephys_update',tsoil_PEM(:,:,:)

!     write(*,*) 'tssoil saved',Tsoilave_before_saved(:,:,:)
!

!     write(*,*) 'diff avg', tsoil_PEM(:,:,:) -Tsoilave_before_saved(:,:,:)
   

     call computeice_table(timelen,ngrid,nslope,nsoilmx_PEM,tsoil_phys_PEM_timeseries_update,tsurf_phys_GCM_timeseries_update,q_co2_PEM_phys,q_h2o_PEM_phys, & 
                           ps_phys_timeseries, ice_depth)
     DO ig = 1,ngrid
         DO islope = 1,nslope
          write(*,*)'zice =',ig,islope, ice_depth(ig,islope)
         ENDDO
      ENDDO
      
      print *, "10 is okay"

! call update_soil(timelen,ngrid,nslope,nsoilmx_PEM,tendencies_h2o_ice_phys_slope,tendencies_co2_ice_phys_slope,co2ice_slope,qsurf_slope(:,6,:),ps_phys, &
!        cell_area,ice_depth,TI_PEM)

      print *, "11 is okay"
   

    call recomp_tend_co2_slope(tendencies_co2_ice_phys_slope,tendencies_co2_ice_phys_slope_ini,vmr_co2_gcm_phys,vmr_co2_pem_phys,ps_GCM_2,&     
                               global_ave_press_GCM,global_ave_press_new,timelen,ngrid,nslope)

          print *, "12 is okay"
      year_iter=year_iter+dt_pem

!We check with we should stop


      call criterion_ice_stop_water_slope(cell_area,ini_surf_h2o,qsurf_slope(:,6,:),STOPPING_water,ngrid,initial_h2o_ice)


          print *, "13 is okay"
      call criterion_ice_stop_slope(cell_area,ini_surf_co2,co2ice_slope,STOPPING_co2,ngrid,initial_co2_ice_sublim_slope,global_ave_press_GCM,global_ave_press_new,nslope)

 print *, "criterion ok"


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


!   call interpolate_TIPEM_TIGCM(ngrid,nslope,nsoilmx_PEM,nsoilmx,TI_PEM,TI_GCM_phys)


!! NOW that temperature averaged have been fixed, let's update the instantanous
!temperature

do ig = 1,ngrid
 do islope = 1,nslope

   tsurf_slope(ig,islope) =  tsurf_slope(ig,islope) - (tsurf_ave_phys_inst(ig,islope)-tsurf_ave_phys(ig,islope))
!tsurf_slope(ig,islope) = tsurf_ave_phys(ig,islope) 
  do iloop = 1,nsoilmx
     tsoil_slope(ig,iloop,islope) = tsoil_slope(ig,iloop,islope) - (tsoil_ave_phys_inst(ig,iloop,islope) - tsoil_PEM(ig,iloop,islope))
!tsoil_slope(ig,iloop,islope) = tsoil_PEM(ig,iloop,islope) 
  enddo
       
 enddo
enddo



      DO ig = 1,ngrid
          tsurf(ig) = 0.
             DO islope = 1,nslope
                 tsurf(ig) = tsurf(ig) + tsurf_slope(ig,islope) &
                            * subslope_dist(ig,islope)      
             ENDDO
      ENDDO ! of DO ig=1,ngrid

      DO ig = 1,ngrid
         DO iloop = 1,nsoilmx
          tsoil(ig,iloop) = 0.
             DO islope = 1,nslope
                  tsoil(ig,iloop) =  tsoil(ig,iloop) + tsoil_slope(ig,iloop,islope) &
                            * subslope_dist(ig,islope)      
             ENDDO
          ENDDO
      ENDDO ! of DO ig=1,ngrid


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


             call physdem0("restartfi_evol.nc",longitude,latitude, &
                        nsoilmx,ngrid,nlayer,nq,              &
                        ptimestep,pday,0.,cell_area,          &
                        albedodat,inertiedat,zmea,zstd,zsig,zgam,zthe, &
                        hmons,summit,base,nslope,def_slope,   &
                        subslope_dist)


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

