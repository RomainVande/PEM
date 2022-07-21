   SUBROUTINE update_soil(timelen,ngrid,nslope,nsoil_PEM,tend_h2oglaciers,tend_co2glaciers,co2ice,waterice,ps_new,&
                          cellarea,ice_depth,TI_new)


 USE comsoil_h, only:  inertiedat, volcapa
 USE comsoil_h_PEM, only: layer_PEM 
 USE vertical_layers_mod, ONLY: ap,bp
! Input: 
 INTEGER,INTENT(IN) :: ngrid, nslope, nsoil_PEM
 REAL,INTENT(IN) :: tend_h2oglaciers(ngrid,nslope),tend_co2glaciers(ngrid,nslope)
 REAL,INTENT(IN) :: ps_new(ngrid)
 REAL,INTENT(IN) :: cellarea(ngrid)
 REAL,INTENT(IN) :: co2ice(ngrid,nslope)
 REAL,INTENT(IN) :: waterice(ngrid,nslope)
 REAL,INTENT(INOUT) :: ice_depth(ngrid,nslope)
 

! Outputs:

 REAL,INTENT(OUT) :: TI_new(ngrid,nsoil_PEM,nslope)
  
! Constants:

 REAL ::  alpha = 0.2
 REAL ::  beta = 1.08e7
 REAL ::  To = 273.15
 REAL ::  R =  8.314
 REAL ::  L =  51058.
 REAL ::  inertie_thresold = 800 ! look for ice
 REAL ::  inertie_co2glaciers = 2000 ! Mellon et al. 2000
 REAL ::  inertie_averaged = 250 ! Mellon et al. 2000
 REAL ::  ice_inertia = 2000  ! Inertia of ice
 REAL ::  P610 = 610.0
 REAL ::  m_h2o = 18.01528E-3
 REAL ::  m_co2 = 44.01E-3  ! CO2 molecular mass (kg/mol)   
 REAL ::  m_noco2 = 33.37E-3  ! Non condensible mol mass (kg/mol)   
 REAL ::  A,B,mmean             

! Local variables:

 INTEGER :: ig,islope,iloop,iref,k
 REAL :: regolith_inertia(ngrid,nslope) ! TI of the regolith 
 REAL :: d(ngrid,nslope)
 REAL :: p_avg_new
 REAL :: zplev_new(ngrid,2)
 REAL :: pp_water_new(ngrid)
 REAL :: T_h2o(ngrid)
 REAL :: Total_surface
 INTEGER :: ispermanent_co2glaciers(ngrid,nslope)
 INTEGER :: ispermanent_h2oglaciers(ngrid,nslope)
 REAL :: es,qsat,x

! 0. Initialisation


  do islope = 1,nslope
    regolith_inertia(:,islope) = inertiedat(:,1)
  enddo

write(*,*) 'in update'

  Total_surface = sum(cellarea)
  p_avg_new = 0.
  do ig = 1,ngrid
    p_avg_new = p_avg_new + ps_new(ig)*cellarea(ig)/Total_surface
  enddo
  write(*,*) 'pnew=',p_avg_new
  do l=1,2
    do ig=1,ngrid
      zplev_new(ig,l) = ap(l)  + bp(l)*ps_new(ig)
    enddo
  enddo

  do ig = 1,ngrid
    do islope = 1,nslope
     if((abs(tend_h2oglaciers(ig,islope)).lt.1e-5).and.(abs(waterice(ig,islope)).gt.0)) then
        ispermanent_h2oglaciers(ig,islope) = 1
     else
        ispermanent_h2oglaciers(ig,islope) = 0
     endif 

     if((abs(tend_co2glaciers(ig,islope)).lt.1e-5).and.(abs(co2ice(ig,islope)).gt.0)) then
        ispermanent_co2glaciers(ig,islope) = 1
     else
        ispermanent_co2glaciers(ig,islope) = 0
     endif
    if(ig.eq.12) write(*,*) 'permanent 12',islope, ispermanent_co2glaciers(12,islope), ispermanent_h2oglaciers(12,islope)
    enddo

  enddo


 ispermanent_co2glaciers(:,:) = 0
 ispermanent_h2oglaciers(:,:) = 0

! 1. Modification of the regolith thermal inertia.

! a) First check if the surface is not anymore covered by ice

  do ig = 1,ngrid
   do islope = 1,nslope
     if((regolith_inertia(ig,islope).gt.inertie_thresold).and.(ispermanent_co2glaciers(ig,islope).eq.0).and.(ispermanent_h2oglaciers(ig,islope).eq.0)) then
        regolith_inertia(ig,islope) = inertie_averaged
        write(*,*) 'in there',ig,islope 
     endif
    enddo
  enddo


 

! b) Then we update the thermal inertia

  do ig=1,ngrid
    do islope=1,nslope
	      d(ig,islope) = ((regolith_inertia(ig,islope)*regolith_inertia(ig,islope))/(volcapa*alpha*P610**0.6))**(-1/(0.11*log(P610/beta)))
 
    write(*,*) 'regolith_inertia before',regolith_inertia(ig,islope)

   if(regolith_inertia(ig,islope).lt.inertie_thresold+1) then !           we are modifying the regolith properties, not ice
          regolith_inertia(ig,islope) = sqrt(volcapa*alpha*(p_avg_new**0.6)*d(ig,islope)**(-0.11*log(p_avg_new/beta)))
    write(*,*) 'regolith_inertia after',regolith_inertia(ig,islope)
endif 
     enddo
   enddo


! b) Third: put infinite ice below permanent CO2/H2O glaciers


!  3. Looking for the icedepth

  do ig=1,ngrid
    do islope=1,nslope
      if ((ispermanent_co2glaciers(ig,islope).eq.1).or.(ispermanent_h2oglaciers(ig,islope).eq.1)) then
           ice_depth(ig,islope) = 0.
          write(*,*), 'ig,is,ip',ig,islope,ispermanent_co2glaciers(ig,islope),ispermanent_h2oglaciers(ig,islope)
   
      endif
    enddo
   enddo


!  4. Build new TI for the PEM
! a) For the regolith
  do ig=1,ngrid
    do islope=1,nslope
! 4.0 FIrst if permanent ice
   if (ice_depth(ig,islope).lt.1e-10) then
      do iloop = 1,nsoil_PEM
       TI_new(ig,iloop,islope)=ice_inertia
      enddo
     else 
      ! 4.1 find the index of the mixed layer
      iref=0 ! initialize iref
      do k=1,nsoil_PEM ! loop on layers
        if (ice_depth(ig,islope).ge.layer_PEM(k)) then
          iref=k ! pure regolith layer up to here
        else
         ! correct iref was obtained in previous cycle
         exit
        endif
       
      enddo

! put infinite ice below permanent CO2/H2O glaciers


      if (ig.eq.12) write(*,*) 'iref=',iref,ice_depth(ig,islope)
      ! 4.2 Build the new ti
      do iloop=1,iref
      if (ig.eq.12) write(*,*) 'in=',iloop,islope
         TI_new(ig,iloop,islope) =regolith_inertia(ig,islope) 
      enddo
      if (iref.lt.nsoil_PEM) then
        if (iref.ne.0) then
          ! mixed layer
           TI_new(ig,iref+1,islope)=sqrt((layer_PEM(iref+1)-layer_PEM(iref))/ &
            (((ice_depth(ig,islope)-layer_PEM(iref))/(regolith_inertia(ig,islope)**2))+ &
                      ((layer_PEM(iref+1)-ice_depth(ig,islope))/(ice_inertia**2))))
        else ! first layer is already a mixed layer
          ! (ie: take layer(iref=0)=0)
          TI_new(ig,1,islope)=sqrt((layer_PEM(1))/ &
                          (((ice_depth(ig,islope))/(regolith_inertia(ig,islope)**2))+ &
                           ((layer_PEM(1)-ice_depth(ig,islope))/(ice_inertia**2))))
        endif ! of if (iref.ne.0)        
        ! lower layers of pure ice
        do iloop=iref+2,nsoil_PEM
          TI_new(ig,iloop,islope)=ice_inertia
        enddo
      endif ! of if (iref.lt.(nsoilmx))
      endif ! permanent glaciers

     enddo !islope
    enddo !ig







!=======================================================================
      RETURN
      END
