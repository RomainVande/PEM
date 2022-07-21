      subroutine soil_settings_PEM(ngrid,nslope,nsoil_PEM,nsoil_GCM,
     &   TI_GCM,TI_PEM)


!      use netcdf
      use comsoil_h_PEM, only: layer_PEM, mlayer_PEM
      use comsoil_h, only: inertiedat,layer,mlayer, volcapa
      use iostart, only: inquire_field_ndims, get_var, get_field,
     &                   inquire_field, inquire_dimension_length

      implicit none

!======================================================================
!  Author: Ehouarn Millour (07/2006)
!
!  Purpose: Read and/or initialise soil depths and properties
!
! Modifications: Aug.2010 EM : use NetCDF90 to load variables (enables using
!                      r4 or r8 restarts independently of having compiled
!                      the GCM in r4 or r8)
!                June 2013 TN : Possibility to read files with a time axis
!
!
!  This subroutine reads from a NetCDF file (opened by the caller)
!  of "startfi.nc" format.
!  The various actions and variable read/initialized are:
!  1. Check out the number of soil layers (in datafile); if it isn't equal
!     to nsoil, then some interpolation will be required
!     Also check if data in file "startfi.nc" is in older format (ie:
!     thermal inertia was depth-independent; and there was no "depth"
!     coordinate.
!     Read/build layer (and midlayer) depths
!  2. Read volumetric specific heat (or initialise it to default value)
!  3. Read Thermal inertia
!  4. Read soil temperatures
!  5. Interpolate thermal inertia and temperature on the new grid, if
!     necessary
!======================================================================

!======================================================================
!  arguments
!  ---------
!  inputs:
      integer,intent(in) :: ngrid	! # of horizontal grid points
      integer,intent(in) :: nslope	! # of subslope wihtin the mesh
      integer,intent(in) :: nsoil_PEM	! # of soil layers in the PEM
      integer,intent(in) :: nsoil_GCM	! # of soil layers in the GCM
      real,intent(in) :: TI_GCM(ngrid,nsoil_GCM,nslope)	! # of soil layers in the GCM
      real,intent(inout) :: TI_PEM(ngrid,nsoil_PEM,nslope)	! # of soil layers in the PEM

!======================================================================
! local variables:
      integer ig,iloop,islope	! loop counters
      logical found

      real alpha,lay1 ! coefficients for building layers
      real xmin,xmax ! to display min and max of a field

      real,parameter :: default_volcapa=1.e6
      
!======================================================================

! 1. Depth coordinate
! -------------------

! 1.4 Build mlayer(), if necessary
!      if (interpol) then
      ! default mlayer distribution, following a power law:
      !  mlayer(k)=lay1*alpha**(k-1/2)
        lay1=2.e-4
	alpha=2
        do iloop=0,nsoil_PEM
	  mlayer_PEM(iloop)=lay1*(alpha**(iloop-0.5))
	enddo
!      endif

! 1.5 Build layer(); following the same law as mlayer()
      ! Assuming layer distribution follows mid-layer law:
      ! layer(k)=lay1*alpha**(k-1)
      lay1=sqrt(mlayer_PEM(0)*mlayer_PEM(1))
      alpha=mlayer_PEM(1)/mlayer_PEM(0)
      do iloop=1,nsoil_PEM
        layer_PEM(iloop)=lay1*(alpha**(iloop-1))
      enddo


! 2. Volumetric heat capacity (note: it is declared in comsoil_h)
! ---------------------------
! "volcapa" is (so far) 0D and written in "controle" table of startfi file 
! volcapa is read or set when "controle" is read (see tabfi.F)
! Just in case, we check here that it is not zero. If it is, we
! set it to "default_volcapa"

      if (volcapa.le.0.0) then
        write(*,*)'soil_settings: Warning, volcapa = ',volcapa
	write(*,*)'               That doesn t seem right'
        write(*,*)'        Initializing Volumetric heat capacity to ',
     &             default_volcapa
	volcapa=default_volcapa
      endif


! 3. Thermal inertia (note: it is declared in comsoil_h)
! ------------------


      do ig = 1,ngrid
        do islope = 1,nslope
          do iloop = 1,nsoil_GCM
            TI_PEM(ig,iloop,islope) = TI_GCM(ig,iloop,islope)
          enddo
          if(nsoil_PEM.gt.nsoil_GCM) then
           do iloop = nsoil_GCM+1,nsoil_PEM
             TI_PEM(ig,iloop,islope) = TI_GCM(ig,nsoil_GCM,islope)
           enddo
          endif
        enddo
      enddo
      end subroutine soil_settings_PEM
