
subroutine setperturb

!  Random noise

use vars
use params
use microphysics, only: micro_field, index_water_vapor
use grid

implicit none

integer i,j,k,ptype,it,jt
real rrr,ranf_, L, x0, y0, r, noise
real xxx,yyy,zzz

call ranset_(3*rank)

ptype = perturb_type

select case (ptype)

  case(0)

     do k=1,nzm
      do j=1,ny
       do i=1,nx
         rrr=1.-2.*ranf_()
         if(k.le.5) then
            t(i,j,k)=t(i,j,k)+0.02*rrr*(6-k)
         endif
         if(k.le.4.and..not.dosmagor) then
            tke(i,j,k)=0.04*(5-k)
         endif 
       end do
      end do
     end do

  case(1)

     do k=1,nzm
      do j=1,ny
       do i=1,nx
         rrr=1.-2.*ranf_()
         if(q0(k).gt.6.e-3.and..not.dosmagor) then
            t(i,j,k)=t(i,j,k)+0.1*rrr
            tke(i,j,k)=1.
         endif
       end do
      end do
     end do

  case(2) ! warm bubble

     if(masterproc) then
       print*, 'initialize with warm bubble:'
       print*, 'bubble_x0=',bubble_x0
       print*, 'bubble_y0=',bubble_y0
       print*, 'bubble_z0=',bubble_z0
       print*, 'bubble_radius_hor=',bubble_radius_hor
       print*, 'bubble_radius_ver=',bubble_radius_ver
       print*, 'bubble_dtemp=',bubble_dtemp
       print*, 'bubble_dq=',bubble_dq
     end if

     call task_rank_to_index(rank,it,jt)
     do k=1,nzm
       zzz = z(k)
       do j=1,ny
         yyy = dy*(j+jt)
         do i=1,nx
          xxx = dx*(i+it)
           if((xxx-bubble_x0)**2+YES3D*(yyy-bubble_y0)**2.lt.bubble_radius_hor**2 &
            .and.(zzz-bubble_z0)**2.lt.bubble_radius_ver**2) then
              rrr = cos(pi/2.*(xxx-bubble_x0)/bubble_radius_hor)**2 &
               *cos(pi/2.*(yyy-bubble_y0)/bubble_radius_hor)**2 &
               *cos(pi/2.*(zzz-bubble_z0)/bubble_radius_ver)**2
              t(i,j,k) = t(i,j,k) + bubble_dtemp*rrr
              micro_field(i,j,k,index_water_vapor) = &
                  micro_field(i,j,k,index_water_vapor) + bubble_dq*rrr
           end if
         end do
       end do
     end do

  case(3)   ! gcss wg1 smoke-cloud case

     do k=1,nzm
      do j=1,ny
       do i=1,nx
         rrr=1.-2.*ranf_()
         if(q0(k).gt.0.5e-3.and..not.dosmagor) then
            t(i,j,k)=t(i,j,k)+0.1*rrr
            tke(i,j,k)=1.
         endif
       end do
      end do
     end do

  case(4)  ! gcss wg1 arm case

     do k=1,nzm
      do j=1,ny
       do i=1,nx
         rrr=1.-2.*ranf_()
         if(z(k).le.200.) then
            t(i,j,k)=t(i,j,k)+0.1*rrr*(1.-z(k)/200.)
         endif
         if(z(k).le.150..and..not.dosmagor) then
            tke(i,j,k)=0.15*(1.-z(k)/150.)
         endif
       end do
      end do
     end do

  case(5)  ! gcss wg1 BOMEX case

     do k=1,nzm
      do j=1,ny
       do i=1,nx
         rrr=1.-2.*ranf_()
         if(z(k).le.1600.) then
            t(i,j,k)=t(i,j,k)+0.1*rrr
            micro_field(i,j,k,index_water_vapor)= &
                      micro_field(i,j,k,index_water_vapor)+0.025e-3*rrr
         endif
         if(z(k).le.3000..and..not.dosmagor) then
            tke(i,j,k)=1.-z(k)/3000.
         endif
       end do
      end do
     end do

  case(6)  ! GCSS Lagragngian ASTEX


     do k=1,nzm
      do j=1,ny
       do i=1,nx
         rrr=1.-2.*ranf_()
         if(q0(k).gt.6.e-3) then
            t(i,j,k)=t(i,j,k)+0.1*rrr
            micro_field(i,j,k,index_water_vapor)= &
                      micro_field(i,j,k,index_water_vapor)+2.5e-5*rrr
            tke(i,j,k)=1.
         endif
       end do
      end do
     end do

  case(7) ! Aggregated intial conditions for water vapor

     if(masterproc) then
       print*, 'initialize aggregated initial conditions'
     end if

     L = float(nx_gl)*dx 
     x0 = L/2.
     y0 = L/2.

     call task_rank_to_index(rank,it,jt)
     do k=1,nzm
       zzz = z(k)
       do j=1,ny
         yyy = dy*(j+jt)   ! Global y coordinate
         do i=1,nx
          xxx = dx*(i+it)  ! Global x coordinate
          r = sqrt((xxx-x0)**2 + YES3D*(yyy-y0)**2)  
          rrr = exp(-zzz/3000.)
          micro_field(i,j,k,index_water_vapor) = rrr*(0.016 - (2./L)*0.008*r) ! Aggregated  water vapor field 
         end do
       end do
     end do 

 case(8) ! Strip-Aggregated IC 

     if(masterproc) then
       print*, 'initialize strip-aggregated initial conditions'
     end if

     L = float(nx_gl)*dx 
     x0 = L/2.
     y0 = L/2.

     call task_rank_to_index(rank,it,jt)
     do k=1,nzm
       zzz = z(k)
       do j=1,ny
         yyy = dy*(j+jt)   ! Global y coordinate
         do i=1,nx
          xxx = dx*(i+it)  ! Global x coordinate
          r = sqrt((xxx-x0)**2 + YES3D*(yyy-y0)**2)  
          rrr = exp(-zzz/3000.)
          micro_field(i,j,k,index_water_vapor) = rrr*(0.016 - (2./L)*0.008*ABS(xxx-x0)) ! Strip-aggregated water vapor field
          noise=1.-2.*ranf_()
          if(k.le.5) then
            t(i,j,k)=t(i,j,k)+0.02*noise*(6-k) ! Add white noise to temp field
          endif
          if(k.le.4.and..not.dosmagor) then
            tke(i,j,k)=0.04*(5-k)
          endif
         end do
       end do
     end do 



  case default

       if(masterproc) print*,'perturb_type is not defined in setperturb(). Exitting...'
       call task_abort()

end select


if(doscalar) then
   tke(1:nx,1:ny,1:nzm) = qv(1:nx,1:ny,1:nzm)
end if

end

