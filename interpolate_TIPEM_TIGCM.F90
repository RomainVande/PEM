   subroutine interpolate_TIPEM_TIGCM(ngrid,nslope,nsoil_PEM,nsoil_GCM,TI_PEM,TI_GCM)



      implicit none


!======================================================================
!  arguments
!  ---------
!  inputs:
      integer,intent(in) :: ngrid	! # of horizontal grid points
      integer,intent(in) :: nslope	! # of subslope wihtin the mesh
      integer,intent(in) :: nsoil_PEM	! # of soil layers in the PEM
      integer,intent(in) :: nsoil_GCM	! # of soil layers in the GCM
      real,intent(in) :: TI_PEM(ngrid,nsoil_PEM,nslope)	! # of soil layers in the PEM
      real,intent(inout) :: TI_GCM(ngrid,nsoil_GCM,nslope)	! # of soil layers in the PEM

!local variable
      integer :: ig,islope,iloop



     do ig = 1,ngrid
       do islope = 1,nslope
         do iloop = 1,nsoil_GCM
           TI_GCM(ig,iloop,islope) = TI_PEM(ig,iloop,islope)
         enddo
       enddo
     enddo


  end subroutine interpolate_TIPEM_TIGCM
