program crm

!       Main module.

use vars
use hbuffer
use microphysics
use tracers
implicit none

integer k, icyc, nn, nstatsteps
real dummy(nz), tmp
real fluxbtmp(nx,ny), fluxttmp(nx,ny) !bloss
!-------------------------------------------------------------------
! determine the rank of the current task and of the neighbour's ranks

call task_init() 
!------------------------------------------------------------------
! print time, version, etc

if(masterproc) call header()	
!------------------------------------------------------------------
! Initialize timing library.  2nd arg 0 means disable, 1 means enable

   call t_setoptionf (1, 0)
   call t_initializef ()

   call t_startf ('total')
   call t_startf ('initialize')
!------------------------------------------------------------------

call init()     ! initialize some statistics arrays
call setparm()	! set all parameters and constants

!------------------------------------------------------------------
! Initialize or restart from the save-dataset:

if(nrestart.eq.0) then
   day=day0 
   call setgrid() ! initialize vertical grid structure
   call setdata() ! initialize all variables
elseif(nrestart.eq.1) then
   call read_all()
   call setgrid() ! initialize vertical grid structure
   call diagnose()
   call micro_init()  !initialize microphysics
elseif(nrestart.eq.2) then  ! branch run
   call read_all()
   call setgrid() ! initialize vertical grid structure
   call diagnose()
   call setparm() ! overwrite the parameters
   call micro_init()  !initialize microphysics
   nstep = 0
   day0 = day
else
   print *,'Error: confused by value of NRESTART'
   call task_abort() 
endif

if(mod(nsubdomains_x,4).ne.0) then
  print *,'Error: nsubdomains_x not divisible by 4. Aborting'
  call task_abort()
endif 

call stat_2Dinit()
call tracers_init() ! initialize tracers
call setforcing()
if(masterproc) call printout()
!------------------------------------------------------------------
!  Initialize statistics buffer:

call hbuf_init()
	
!------------------------------------------------------------------
nstatis = nstat/nstatfrq
nstat = nstatis * nstatfrq
nstatsteps = 0
call t_stopf ('initialize')
!------------------------------------------------------------------
!   Main time loop    
!------------------------------------------------------------------

do while(nstep.lt.nstop.and.nelapse.gt.0) 
        
  nstep = nstep + 1
  time = time + dt
  day = day0 + nstep*dt/86400.
  nelapse = nelapse - 1
!------------------------------------------------------------------
!  Check if the dynamical time step should be decreased 
!  to handle the cases when the flow being locally linearly unstable
!------------------------------------------------------------------

  ncycle = 1

  call kurant()

  total_water_before = 0.
  total_water_after = 0.
  total_water_evap = 0.
  total_water_prec = 0.
  total_water_ls = 0.

  do icyc=1,ncycle

     icycle = icyc
     dtn = dt/ncycle
     dt3(na) = dtn
     dtfactor = dtn/dt

     if(mod(nstep,nstatis).eq.0.and.icycle.eq.ncycle) then
        nstatsteps = nstatsteps + 1
        dostatis = .true.
        if(masterproc) print *,'Collecting statistics...'
     else
        dostatis = .false.
     endif

     !bloss:make special statistics flag for radiation,since it's only updated at icycle==1.
     dostatisrad = .false.
     if(mod(nstep,nstatis).eq.0.and.icycle.eq.1) dostatisrad = .true.

!---------------------------------------------
!  	the Adams-Bashforth scheme in time

     call abcoefs()
 
!---------------------------------------------
!  	initialize stuff: 
	
     call zero()

     total_water_before = total_water_before + total_water()

!-----------------------------------------------------------
!       Buoyancy term:
	     
     call buoyancy()

!------------------------------------------------------------
!       Large-scale and surface forcing:

     total_water_ls =  total_water_ls - total_water()
    
     call forcing()

!----------------------------------------------------------
!       Nadging:

     call nudging()

!----------------------------------------------------------
!   	suppress turbulence near the upper boundary (spange):

     if(dodamping) call damping()

     total_water_ls =  total_water_ls + total_water()

!----------------------------------------------------------
!      Update the subdomain's boundaries for velocity

     call boundaries(0)

!---------------------------------------------------------
!	SGS TKE equation:     	
	   
     if(dosgs) call tke_full()
!---------------------------------------------------------
!   Ice fall-out
   
      if(docloud) then
          call ice_fall()
      end if

!---------------------------------------------------------
!        Update boundaries for scalars, sst,  SGS exchange coefficients 

     call boundaries(2)

!-----------------------------------------------
!       advection of momentum:

     call advect_mom()

!-----------------------------------------------
!   	surface fluxes:

     if(dosurface) then

       call surface()

     end if
!----------------------------------------------------------
!	SGS diffusion of momentum:

     if(dosgs) call diffuse_mom()

!-----------------------------------------------------------
!       Coriolis force:
	     
     if(docoriolis) call coriolis()
	 
!---------------------------------------------------------
!       compute rhs of the Poisson equation and solve it for pressure. 

     call pressure()


!---------------------------------------------------------
!       find velocity field at n+1/2 timestep needed for advection of scalars:
	 
     call adams()

!----------------------------------------------------------
!     Update boundaries for velocity fields to use for advection of scalars:

     call boundaries(1)

!---------------------------------------------------------
!      advection of scalars :

     call advect_scalar(t,tadv,twle,t2leadv,t2legrad,twleadv,.true.)
     
     if(dosgs.and..not.dosmagor) then
      call advect_scalar(tke,dummy,tkewle,dummy,dummy,dummy,.false.)
     else if(doscalar) then
      call advect_scalar(tke,dummy,tkewle,s2leadv,s2legrad,swleadv,.true.)
     end if

!
!    Advection of microphysics prognostics:
!

     do k = 1,nmicro_fields
        if(  k.eq.index_water_vapor   &! transport water-vapor variable no metter what
        .or. docloud.and.flag_precip(k).ne.1    & ! transport non-precipitation vars
        .or. doprecip.and.flag_precip(k).eq.1 ) then
           call advect_scalar(micro_field(:,:,:,k),mkadv(:,k),mkwle(:,k),dummy,dummy,dummy,.false.)
        end if
     end do


!
!   Precipitation fallout:
!
    if(doprecip) then

       total_water_prec = total_water_prec + total_water()
 
       call micro_precip_fall()

       total_water_prec = total_water_prec - total_water()


    end if

 ! advection of tracers:

     if(dotracers) then

        do k = 1,ntracers
         call advect_scalar(tracer(:,:,:,k),tradv(:,k),trwle(:,k),dummy,dummy,dummy,.false.)
        end do

     end if

!---------------------------------------------------------
!      diffusion of scalars :

!        Update boundaries for scalars:

      if(dosgs) call boundaries(3)

      call diffuse_scalar(t,fluxbt,fluxtt,tdiff,twsb, &
                           t2lediff,t2lediss,twlediff,.true.)
     
      if(.not.dosmagor) then
          call diffuse_scalar(tke,fzero,fzero,dummy,tkewsb, &
                                    dummy,dummy,dummy,.false.)
      else if(doscalar) then
          call diffuse_scalar(tke,fluxbq,fluxtq,dummy,tkewsb, &
                           s2lediff,s2lediss,swlediff,.true.)
      end if


!
!    diffusion of microphysics prognostics:
!
      call micro_flux()

      total_water_evap = total_water_evap - total_water()

      do k = 1,nmicro_fields
        if(   k.eq.index_water_vapor             &! transport water-vapor variable no metter what
         .or. docloud.and.flag_precip(k).ne.1    & ! transport non-precipitation vars
         .or. doprecip.and.flag_precip(k).eq.1 ) then
           fluxbtmp(1:nx,1:ny) = fluxbmk(1:nx,1:ny,k)
           fluxttmp(1:nx,1:ny) = fluxtmk(1:nx,1:ny,k)
           call diffuse_scalar(micro_field(:,:,:,k),fluxbtmp,fluxttmp, &
                mkdiff(:,k),mkwsb(:,k), dummy,dummy,dummy,.false.)
       end if
      end do

      total_water_evap = total_water_evap + total_water()

 ! diffusion of tracers:

      if(dotracers) then

        call tracers_flux()

        do k = 1,ntracers

          fluxbtmp = fluxbtr(:,:,k)
          fluxttmp = fluxttr(:,:,k)
          call diffuse_scalar(tracer(:,:,:,k),fluxbtmp,fluxttmp, &
               trdiff(:,k),trwsb(:,k), &
               dummy,dummy,dummy,.false.)
!!$          call diffuse_scalar(tracer(:,:,:,k),fluxbtr(:,:,k),fluxttr(:,:,k),trdiff(:,k),trwsb(:,k), &
!!$                           dummy,dummy,dummy,.false.)
 
        end do

      end if

!-----------------------------------------------------------
!    Convert back from Courant numbers and Updatee the velocity field:

      call uvw()

!-----------------------------------------------------------
!       Handle upper boundary for scalars

     if(doupperbound) call upperbound()

!-----------------------------------------------------------
!       Cloud condensation/evaporation and precipitation processes:

      if(docloud.or.dosmoke) call micro_proc()

!----------------------------------------------------------
!  Tracers' physics:

      call tracers_physics()

!-----------------------------------------------------------
!	Radiation

      if(dolongwave.or.doshortwave) then 
	call radiation()     
      end if

!-----------------------------------------------------------
!    Compute diagnostic fields:

      call diagnose()

      total_water_after = total_water_after + total_water()

!----------------------------------------------------------
! Rotate the dynamic tendency arrays for Adams-bashforth scheme:

      nn=na
      na=nc
      nc=nb
      nb=nn

   end do ! icycle	
          
!----------------------------------------------------------
!  collect statistics, write save-file, etc.

   call stepout(nstatsteps)
  
!----------------------------------------------------------
!----------------------------------------------------------

end do ! main loop

!----------------------------------------------------------
!----------------------------------------------------------

   call t_stopf('total')
   if(masterproc) call t_prf(rank)

call task_stop()

end program crm
