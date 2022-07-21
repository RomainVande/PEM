module comsoil_h_PEM

implicit none
! nsoilmx : number of subterranean layers
!EM: old soil routine:      integer, parameter :: nsoilmx = 10
  integer, parameter :: nsoilmx_PEM = 25 

  real,save,allocatable,dimension(:) :: layer_PEM      ! soil layer depths
  real,save,allocatable,dimension(:) :: mlayer_PEM     ! soil mid-layer depths
  real,save,allocatable,dimension(:,:,:) :: TI_PEM ! soil thermal inertia
  real,save :: volcapa    ! soil volumetric heat capacity
       ! NB: volcapa is read fromn control(35) from physicq start file
       !     in physdem (or set via tabfi, or initialized in
       !                 soil_settings.F)

  ! variables (FC: built in firstcall in soil.F)
  REAL,SAVE,ALLOCATABLE :: tsoil_PEM(:,:,:)       ! sub-surface temperatures (K)
  real,save,allocatable :: mthermdiff_PEM(:,:)  ! (FC) mid-layer thermal diffusivity
  real,save,allocatable :: thermdiff_PEM(:,:)   ! (FC) inter-layer thermal diffusivity
  real,save,allocatable :: coefq_PEM(:)         ! (FC) q_{k+1/2} coefficients
  real,save,allocatable :: coefd_PEM(:,:)       ! (FC) d_k coefficients
  real,save,allocatable :: alph_PEM(:,:,:)        ! (FC) alpha_k coefficients
  real,save,allocatable :: beta_PEM(:,:,:)        ! beta_k coefficients
  real,save :: mu_PEM
  real,parameter :: fluxgeo = 30e-3 !W/m^2

contains

  subroutine ini_comsoil_h_PEM(ngrid,nslope)
  
  implicit none
  integer,intent(in) :: ngrid ! number of atmospheric columns
  integer,intent(in) :: nslope ! number of slope within a mesh 

    allocate(layer_PEM(nsoilmx_PEM)) !soil layer depths
    allocate(mlayer_PEM(0:nsoilmx_PEM)) ! soil mid-layer depths
    allocate(TI_PEM(ngrid,nsoilmx_PEM,nslope)) ! soil thermal inertia
    allocate(tsoil_PEM(ngrid,nsoilmx_PEM,nslope)) ! soil temperatures
    allocate(mthermdiff_PEM(ngrid,0:nsoilmx_PEM-1))
    allocate(thermdiff_PEM(ngrid,nsoilmx_PEM-1))
    allocate(coefq_PEM(0:nsoilmx_PEM-1))
    allocate(coefd_PEM(ngrid,nsoilmx_PEM-1))
    allocate(alph_PEM(ngrid,nsoilmx_PEM-1,nslope))
    allocate(beta_PEM(ngrid,nsoilmx_PEM-1,nslope))

 
  end subroutine ini_comsoil_h_PEM


  subroutine end_comsoil_h_PEM

  implicit none

    if (allocated(layer_PEM)) deallocate(layer_PEM)
    if (allocated(mlayer_PEM)) deallocate(mlayer_PEM)
    if (allocated(TI_PEM)) deallocate(TI_PEM)
    if (allocated(tsoil_PEM)) deallocate(tsoil_PEM)
    if (allocated(mthermdiff_PEM)) deallocate(mthermdiff_PEM)
    if (allocated(thermdiff_PEM)) deallocate(thermdiff_PEM)
    if (allocated(coefq_PEM)) deallocate(coefq_PEM) 
    if (allocated(coefd_PEM)) deallocate(coefd_PEM)
    if (allocated(alph_PEM)) deallocate(alph_PEM)
    if (allocated(beta_PEM)) deallocate(beta_PEM)

  end subroutine end_comsoil_h_PEM

end module comsoil_h_PEM
