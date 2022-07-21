         ! Update surface ice distribution to iterate to steady state if requested
         if(ice_update)then

            do ig=1,ngrid

               delta_ice = (qsurf(ig,igcm_h2o_ice)-ice_initial(ig))
               
               ! add multiple years of evolution
               qsurf_hist(ig,igcm_h2o_ice) = qsurf_hist(ig,igcm_h2o_ice) + delta_ice*icetstep 

               ! if ice has gone -ve, set to zero
               if(qsurf_hist(ig,igcm_h2o_ice).lt.0.0)then
                  qsurf_hist(ig,igcm_h2o_ice) = 0.0 
               endif

               ! if ice is seasonal, set to zero (NEW)
               if(ice_min(ig).lt.0.01)then
                  qsurf_hist(ig,igcm_h2o_ice) = 0.0 
               endif

            enddo

            ! enforce ice conservation
            ice_tot= SUM(qsurf_hist(:,igcm_h2o_ice)*cell_area(:) )/SUM(cell_area(:))
            qsurf_hist(:,igcm_h2o_ice) = qsurf_hist(:,igcm_h2o_ice)*(icesrf/ice_tot)

         endif
