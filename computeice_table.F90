   SUBROUTINE computeice_table(timelen,ngrid,nslope,nsoil_PEM,tsoil,tsurf,q_co2,q_h2o,ps,ice_table)
    USE comsoil_h, only:  inertiedat, volcapa
    USE vertical_layers_mod, ONLY: ap,bp
    USE comsoil_h_PEM, only: layer_PEM 

    integer,intent(in) :: timelen,ngrid,nslope,nsoil_PEM
    real,intent(in) :: tsoil(ngrid,nsoil_PEM,nslope,timelen)    ! soil temperature [K]
    real,intent(in) :: tsurf(ngrid,nslope,timelen)              ! surface temperature [K]
    real,intent(in) :: q_co2(ngrid,timelen)                     ! MMR tracer co2 [kg/kg]
    real,intent(in) :: q_h2o(ngrid,timelen)                     ! MMR tracer h2o [kg/kg]
    real,intent(in) :: ps(ngrid,timelen)                        ! surface pressure [Pa]
    real,intent(out) :: ice_table(ngrid,nslope)                 ! ice table [m]


    real :: m_h2o = 18.01528E-3
    real :: m_co2 = 44.01E-3  
    real :: m_noco2 = 33.37E-3  
    real :: A,B,z1,z2
    real :: alpha = -6143.7 
    real :: beta = 29.9074


    integer ig, islope,isoil,it
    real :: mass_mean(ngrid,timelen)                            ! mean mass above the surface
    real :: zplev_mean(ngrid,timelen)                           ! pressure above the surface
    real :: pvapor(ngrid,timelen)                               ! partial pressure above the surface
    real :: pvapor_slope(ngrid,nslope,timelen)                  ! partial pressure above the slopped surface
    real :: rhovapor(ngrid,nslope,timelen)
    real :: rhovapor_avg(ngrid,nslope)                          ! mean vapor_density at the surface yearly averaged


    real :: rho_soil(ngrid,nslope,nsoil_PEM,timelen)            ! water vapor in the soil
    real :: rho_soil_avg(ngrid,nslope,nsoil_PEM)                ! water vapor in the soil yearly averaged

    real :: diff_rho(ngrid,nslope,nsoil_PEM)                    ! difference of vapor content
! 0. Some initializations

      A =(1/m_co2 - 1/m_noco2)
      B=1/m_noco2
! 1. Compute rho surface yearly averaged

!   1.1 Compute the partial pressure of vapor
!a. the molecular mass into the column
     do ig = 1,ngrid
       mass_mean(ig,:) = 1/(A*q_co2(ig,:) +B)
     enddo

     write(*,*) '1.a ok'
! b. pressure level
     do it = 1,timelen
       do ig = 1,ngrid
         zplev_mean(ig,it) = ap(1) + bp(1)*ps(ig,it)
       enddo
     enddo
     write(*,*) '1.b ok'
! c. Vapor pressure
     pvapor(:,:) = mass_mean(:,:)/m_h2o*q_h2o(:,:)*zplev_mean(:,:)
!    1.2   Check if there is frost at the surface and then compute the density
!    at the surface
     do ig = 1,ngrid
       do islope = 1,nslope
         do it = 1,timelen
           psv_surf = exp(alpha/tsurf(ig,islope,it) +beta)
           pvapor_slope(ig,islope,it) = min(psv_surf,pvapor(ig,it))
          if (isnan(pvapor_slope(ig,islope,it))) then
          write(*,*) 'pb vapor',ig,islope,it
         stop
         endif       
        rhovapor(ig,islope,it) = pvapor_slope(ig,islope,it)/tsurf(ig,islope,it)
         enddo
       enddo
     enddo

     write(*,*) '1.c ok'



!    1.3 Density at the surface yearly averaged
     rhovapor_avg(:,:) = SUM(rhovapor(:,:,:),3)/timelen
     write(*,*) '1.3 ok'
! 2. Compute rho soil vapor
   
     do ig = 1,ngrid
       do islope = 1,nslope
         do isoil = 1,nsoil_PEM
            do it = 1,timelen 
              rho_soil(ig,islope,isoil,it) = exp(alpha/tsoil(ig,isoil,islope,it) +beta)/tsoil(ig,isoil,islope,it) 
         
!             write(*,*) 'rho soil', ig,islope,isoil,it,exp(alpha/tsoil(ig,isoil,islope,it) +beta),rho_soil(ig,islope,isoil,it) 
            
            if(isnan(rho_soil(ig,islope,isoil,it))) then
              write(*,*) 'pb rho_soil',ig,islope,isoil,it
             stop
            endif
          enddo
           enddo
       enddo
     enddo

    rho_soil_avg(:,:,:) = SUM( rho_soil(:,:,:,:),4)/timelen
     write(*,*) '2 ok'
! 3. Computing ice table
   
    ice_table (:,:) = 1.e8

    
 do ig = 1,ngrid
       do islope = 1,nslope
         do isoil = 1,nsoil_PEM
            do it = 1,timelen 


             diff_rho(ig,islope,isoil) = rhovapor_avg(ig,islope) - rho_soil_avg(ig,islope,isoil)
!             write(*,*) 'diff =',ig,islope,isoil,diff_rho(ig,islope,isoil),rhovapor_avg(ig,islope) ,rho_soil_avg(ig,islope,isoil)
            enddo
         enddo
       enddo
     enddo








     write(*,*) 'diff_rho ok'


     do ig = 1,ngrid
       do islope = 1,nslope
         if(diff_rho(ig,islope,1) > 0) then
           ice_table(ig,islope) = 0.
         else
           do isoil = 1,nsoil_PEM -1
             if((diff_rho(ig,islope,isoil).lt.0).and.(diff_rho(ig,islope,isoil+1).gt.0.)) then
               z1 = (diff_rho(ig,islope,isoil) - diff_rho(ig,islope,isoil+1))/(layer_PEM(isoil) - layer_PEM(isoil+1))
               z2 = -layer_PEM(isoil+1)*z1 +  diff_rho(ig,islope,isoil+1)
               ice_table(ig,islope) = -z2/z1
               exit
             endif
           enddo
          endif
        enddo
      enddo

     write(*,*) 'end'





!=======================================================================
      RETURN





      END
