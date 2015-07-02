
subroutine task_dispatch(buff,tag)
	
! dispathes the messages according to the field sent.

use vars
use microphysics
use tracers
implicit none
	
real buff(*)	! buff for sending data
integer tag	! tag of the message
integer field   ! id of field
	
field = tag/100000
	
if(field.eq.1) then        
  call task_assign_bnd(u,dimx1_u,dimx2_u,dimy1_u,dimy2_u,nzm,buff,tag)
elseif(field.eq.2) then        
  call task_assign_bnd(v,dimx1_v,dimx2_v,dimy1_v,dimy2_v,nzm,buff,tag)
elseif(field.eq.3) then        
  call task_assign_bnd(w,dimx1_w,dimx2_w,dimy1_w,dimy2_w,nz,buff,tag)	
elseif(field.eq.4) then        
  call task_assign_bnd(t,dimx1_s,dimx2_s,dimy1_s,dimy2_s,nzm,buff,tag)
elseif(field.eq.5) then        
  call task_assign_bnd(tke,dimx1_s,dimx2_s,dimy1_s,dimy2_s,nzm,buff,tag)
elseif(field.eq.6) then        
  call task_assign_bnd(tk,0,nxp1,1-YES3D,nyp1,nzm,buff,tag)
elseif(field.eq.7) then        
  call task_assign_bnd(tkh,0,nxp1,1-YES3D,nyp1,nzm,buff,tag)        
elseif(field.gt.7.and.field.le.7+nmicro_fields) then
  call task_assign_bnd(micro_field(:,:,:,field-7),dimx1_s,dimx2_s,dimy1_s,dimy2_s,nzm,buff,tag)	 
elseif(field.gt.7+nmicro_fields.and.field.le.7+nmicro_fields+ntracers) then
  call task_assign_bnd(tracer(:,:,:,field-7-nmicro_fields),dimx1_s,dimx2_s,dimy1_s,dimy2_s,nzm,buff,tag)	 
end if
	
end
	     
	     
