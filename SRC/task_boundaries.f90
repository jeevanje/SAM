
        subroutine task_boundaries(flag)
        	
!  These routines exchanges overlaping information for various subdomains to
!  be able to calculate various quantities in common /com3d/
!  near the subdomain boundaries.

use vars
use microphysics
use params, only: dosmagor, doscalar, dotracers, dosgs
use tracers
implicit none

integer flag,i

if(flag.eq.0) then

 call task_exchange(u,dimx1_u,dimx2_u,dimy1_u,dimy2_u,nzm,1,1,1,1,1)
 call task_exchange(v,dimx1_v,dimx2_v,dimy1_v,dimy2_v,nzm,1,1,1,1,2)
 ! use w at the top level  - 0s anyway - to exchange the sst boundaries (for
 ! surface fluxes call
 w(1:nx,1:ny,nz) = sstxy(1:nx,1:ny)
 call task_exchange(w,dimx1_w,dimx2_w,dimy1_w,dimy2_w,nz,1,1,1,1,3)	
 sstxy(0:nx,1-YES3D:ny) = w(0:nx,1-YES3D:ny,nz)
 w(0:nx+1,1-YES3D:ny+YES3D,nz) = 0. ! fill it back with 0s

endif

if(flag.eq.1) then

 call task_exchange(u,dimx1_u,dimx2_u,dimy1_u,dimy2_u,nzm,2,3,2+NADV,2+NADV,1)
 call task_exchange(v,dimx1_v,dimx2_v,dimy1_v,dimy2_v,nzm,2+NADV,2+NADV,2,3,2)
 call task_exchange(w,dimx1_w,dimx2_w,dimy1_w,dimy2_w,nz,2+NADV,2+NADV,2+NADV,2+NADV,3)	

endif


if(flag.eq.2) then

 call task_exchange(t,dimx1_s,dimx2_s,dimy1_s,dimy2_s,nzm,3+NADVS,3+NADVS,3+NADVS,3+NADVS,4)
 if(dosgs.and..not.dosmagor.or.doscalar) &
     call task_exchange(tke,dimx1_s,dimx2_s,dimy1_s,dimy2_s,nzm,3,3,3,3,5)
 call task_exchange(tk,0,nxp1,1-YES3D,nyp1,nzm,1,1,1,1,6)
 call task_exchange(tkh,0,nxp1,1-YES3D,nyp1,nzm,1,1,1,1,7)        
 do i = 1,nmicro_fields
    if(   i.eq.index_water_vapor             &
     .or. docloud.and.flag_precip(i).ne.1    &
     .or. doprecip.and.flag_precip(i).eq.1 ) &
     call task_exchange(micro_field(:,:,:,i),dimx1_s,dimx2_s,dimy1_s,dimy2_s,nzm,3+NADVS,3+NADVS,3+NADVS,3+NADVS,7+i)
 end do
 if(dotracers) then
   do i=1,ntracers
     call task_exchange(tracer(:,:,:,i),dimx1_s,dimx2_s,dimy1_s,dimy2_s,nzm,3+NADVS,3+NADVS,3+NADVS,3+NADVS,7+nmicro_fields+i)
   end do
 end if

endif



if(flag.eq.3) then

 call task_exchange(t,dimx1_s,dimx2_s,dimy1_s,dimy2_s,nzm,1,1,1,1,4)
 if(dosgs.and..not.dosmagor.or.doscalar) &
     call task_exchange(tke,dimx1_s,dimx2_s,dimy1_s,dimy2_s,nzm,1,1,1,1,5)
 do i = 1,nmicro_fields
    if(   i.eq.index_water_vapor             &
     .or. docloud.and.flag_precip(i).ne.1    &
     .or. doprecip.and.flag_precip(i).eq.1 ) &
     call task_exchange(micro_field(:,:,:,i),dimx1_s,dimx2_s,dimy1_s,dimy2_s,nzm,1,1,1,1,7+i)
 end do
 if(dotracers) then
   do i=1,ntracers
     call task_exchange(tracer(:,:,:,i),dimx1_s,dimx2_s,dimy1_s,dimy2_s,nzm,1,1,1,1,7+nmicro_fields+i)
   end do
 end if

end if

if(dosgs.and.flag.eq.4) then


endif

end
	
	
