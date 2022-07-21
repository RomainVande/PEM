      subroutine soil_pem(ngrid,nsoil,firstcall, &
               therm_i,                          &
               timestep,tsurf,tsoil,mlayer_PEM,layer_PEM,alph_PEM,beta_PEM)


      use comsoil_h_PEM, only:   &
                          mthermdiff_PEM, thermdiff_PEM, coefq_PEM, &
                          coefd_PEM, mu_PEM,fluxgeo
      use comsoil_h,only: volcapa
      implicit none


!      subroutine soil(ngrid,nsoil,nslope,firstcall,
!     &          therm_i,
!     &          timestep,tsurf,tsoil,
!     &          capcal,fluxgrd)

!-----------------------------------------------------------------------
!  Author: Ehouarn Millour
!
!  Purpose: Compute soil temperature using an implict 1st order scheme
!  
!  Note: depths of layers and mid-layers, soil thermal inertia and 
!        heat capacity are commons in comsoil.h
!-----------------------------------------------------------------------

#include "dimensions.h"
!#include "dimphys.h"

!#include"comsoil.h"


!-----------------------------------------------------------------------
!  arguments
!  ---------
!  inputs:
      integer,intent(in) :: ngrid	! number of (horizontal) grid-points 
      integer,intent(in) :: nsoil	! number of soil layers  
      logical,intent(in) :: firstcall ! identifier for initialization call 
      real,intent(in) :: therm_i(ngrid,nsoil) ! thermal inertia
      real,intent(in) :: timestep	    ! time step
      real,intent(in) :: tsurf(ngrid)   ! surface temperature
      real,intent(in) :: mlayer_PEM(0:nsoil-1) ! layer's mid depth 
      real,intent(in) :: layer_PEM(nsoil) ! layer's  depth 
! outputs:
      real,intent(inout) :: tsoil(ngrid,nsoil) ! soil (mid-layer) temperature
      real,intent(inout) :: alph_PEM(ngrid,nsoil-1)
      real,intent(inout) :: beta_PEM(ngrid,nsoil-1)


!      real capcal(ngrid) ! surface specific heat
!      real fluxgrd(ngrid) ! surface diffusive heat flux

! local saved variables:
!      real,save :: layer(ngridmx,nsoilmx)	! layer depth
!!      real,save :: mthermdiff_PEM(ngridmx,0:nsoilmx-1) ! mid-layer thermal diffusivity
!      real,save :: thermdiff(ngridmx,nsoilmx-1) ! inter-layer thermal diffusivity
!      real,save :: coefq_PEM(0:nsoilmx-1)		! q_{k+1/2} coefficients
!      real,save :: coefd_PEM(ngridmx,nsoilmx-1)	! d_k coefficients
!      real,save :: alph_PEM(ngridmx,nsoilmx-1)	! alph_PEMa_k coefficients
!      real,save :: beta_PEM(ngridmx,nsoilmx-1)	! beta_PEM_k coefficients
!      real,save :: mu_PEM

      
! local variables:
      integer ig,ik

! 0. Initialisations and preprocessing step
 if (firstcall) then
    
! 0.1 Build mthermdiff_PEM(:), the mid-layer thermal diffusivities
      do ig=1,ngrid
        do ik=0,nsoil-1
	  mthermdiff_PEM(ig,ik)=therm_i(ig,ik+1)*therm_i(ig,ik+1)/volcapa
    
	enddo
      enddo


! 0.2 Build thermdiff(:), the "interlayer" thermal diffusivities
      do ig=1,ngrid
        do ik=1,nsoil-1
      thermdiff_PEM(ig,ik)=((layer_PEM(ik)-mlayer_PEM(ik-1))*mthermdiff_PEM(ig,ik) &
                     +(mlayer_PEM(ik)-layer_PEM(ik))*mthermdiff_PEM(ig,ik-1))  &
                         /(mlayer_PEM(ik)-mlayer_PEM(ik-1))

	enddo
      enddo

! 0.3 Build coefficients mu_PEM, q_{k+1/2}, d_k, alph_PEMa_k and capcal
      ! mu_PEM
      mu_PEM=mlayer_PEM(0)/(mlayer_PEM(1)-mlayer_PEM(0))

      ! q_{1/2}
      coefq_PEM(0)=volcapa*layer_PEM(1)/timestep
	! q_{k+1/2}
        do ik=1,nsoil-1
          coefq_PEM(ik)=volcapa*(layer_PEM(ik+1)-layer_PEM(ik)) 		 &
                      /timestep
	enddo
      print *,'0.3 ok'
      do ig=1,ngrid
	! d_k
	do ik=1,nsoil-1
	  coefd_PEM(ig,ik)=thermdiff_PEM(ig,ik)/(mlayer_PEM(ik)-mlayer_PEM(ik-1))
	enddo
	
	! alph_PEM_{N-1}
	alph_PEM(ig,nsoil-1)=coefd_PEM(ig,nsoil-1)/ 			 &
                       (coefq_PEM(nsoil-1)+coefd_PEM(ig,nsoil-1))
        ! alph_PEM_k
        do ik=nsoil-2,1,-1
	  alph_PEM(ig,ik)=coefd_PEM(ig,ik)/(coefq_PEM(ik)+coefd_PEM(ig,ik+1)*    &
                                   (1.-alph_PEM(ig,ik+1))+coefd_PEM(ig,ik))
	enddo

        ! capcal
! Cstar
!        capcal(ig)=volcapa*layer(1)+
!     &              (thermdiff(ig,1)/(mlayer(1)-mlayer(0)))*
!     &              (timestep*(1.-alph_PEM(ig,1)))
! Cs
!        capcal(ig)=capcal(ig)/(1.+mu_PEM*(1.0-alph_PEM(ig,1))*
!     &                         thermdiff(ig,1)/mthermdiff_PEM(ig,0))
!      write(*,*)'soil: ig=',ig,' capcal(ig)=',capcal(ig)
      enddo ! of do ig=1,ngrid

      
  
endif ! of if (firstcall



      IF (.not.firstcall) THEN
!  2. Compute soil temperatures
! First layer:
      do ig=1,ngrid

        tsoil(ig,1)=(tsurf(ig)+mu_PEM*beta_PEM(ig,1)* &  
                                  thermdiff_PEM(ig,1)/mthermdiff_PEM(ig,0))/ &
                   (1.+mu_PEM*(1.0-alph_PEM(ig,1))*&
                    thermdiff_PEM(ig,1)/mthermdiff_PEM(ig,0))
      enddo

! Other layers:
      do ik=1,nsoil-1
        do ig=1,ngrid
          tsoil(ig,ik+1)=alph_PEM(ig,ik)*tsoil(ig,ik)+beta_PEM(ig,ik)
        enddo

!  ********* Tentative flux at the bottom: temporaire

        do ig=1,ngrid

          tsoil(ig,nsoil) = tsoil(ig,nsoil)                      &
        + timestep*fluxgeo/(volcapa*(layer_PEM(nsoil)-layer_PEM(nsoil-1)))


      enddo
      enddo
     ENDIF




!  2. Compute beta_PEM coefficients (preprocessing for next time step)
! Bottom layer, beta_PEM_{N-1}
      do ig=1,ngrid
        beta_PEM(ig,nsoil-1)=coefq_PEM(nsoil-1)*tsoil(ig,nsoil)		 &
                        /(coefq_PEM(nsoil-1)+coefd_PEM(ig,nsoil-1))
      enddo
! Other layers
      do ik=nsoil-2,1,-1
        do ig=1,ngrid
	  beta_PEM(ig,ik)=(coefq_PEM(ik)*tsoil(ig,ik+1)+		 &
                      coefd_PEM(ig,ik+1)*beta_PEM(ig,ik+1))/		 &
                      (coefq_PEM(ik)+coefd_PEM(ig,ik+1)*(1.0-alph_PEM(ig,ik+1)) &
                       +coefd_PEM(ig,ik))
	enddo
      enddo








    
! ******************************************




!  3. Compute surface diffusive flux & calorific capacity
!      do ig=1,ngrid
! Cstar
!        capcal(ig)=volcapa(ig,1)*layer(ig,1)+
!     &              (thermdiff(ig,1)/(mlayer(ig,1)-mlayer(ig,0)))*
!     &              (timestep*(1.-alph_PEM(ig,1)))
! Fstar
!        fluxgrd(ig)=(thermdiff(ig,1)/(mlayer(1)-mlayer(0)))*
!     &              (beta_PEM(ig,1)+(alph_PEM(ig,1)-1.0)*tsoil(ig,1))

!        mu_PEM=mlayer(ig,0)/(mlayer(ig,1)-mlayer(ig,0))
!        capcal(ig)=capcal(ig)/(1.+mu_PEM*(1.0-alph_PEM(ig,1))*
!     &                         thermdiff(ig,1)/mthermdiff_PEM(ig,0))
! Fs
!        fluxgrd(ig)=fluxgrd(ig)+(capcal(ig)/timestep)*
!     &              (tsoil(ig,1)*(1.+mu_PEM*(1.0-alph_PEM(ig,1))*
!     &                         thermdiff(ig,1)/mthermdiff_PEM(ig,0))
!     &               -tsurf(ig)-mu_PEM*beta_PEM(ig,1)*
!     &                          thermdiff(ig,1)/mthermd
      end

