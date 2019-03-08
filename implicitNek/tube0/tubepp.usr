#include "implicit/nek_implicit.f"
c-----------------------------------------------------------------------
c 
c
c-----------------------------------------------------------------------
      subroutine uservp (ix,iy,iz,ieg)
      include 'SIZE'
      include 'TOTAL'
      include 'NEKUSE'

      common /cdsmag/ ediff(lx1,ly1,lz1,lelv)

      ie     = gllel(ieg)
      udiff  = 1
      utrans = 1.

      return
      end
c-----------------------------------------------------------------------
      subroutine userf  (ix,iy,iz,ieg)
      include 'SIZE'
      include 'TOTAL'
      include 'NEKUSE'
      common /cforce/ ffx_new,ffy_new,ffz_new

      ffx = 0.0 ! This value determined from fixed U_b=1 case.
      ffy = 0.0
      ffz = ffz_new

      return
      end
c-----------------------------------------------------------------------
      subroutine userq  (ix,iy,iz,ieg)
      include 'SIZE'
      include 'TOTAL'
      include 'NEKUSE'

      qvol   = 0.0
      source = 0.0

      return
      end
c-----------------------------------------------------------------------
      subroutine userchk
      include 'SIZE'
      include 'TOTAL'

      common /cforce/ ffx_new,ffy_new,ffz_new
      common /vr/ vr(lx1,ly1,lz1,lelt),vphi(lx1,ly1,lz1,lelt),
     & vz2(lx1,ly1,lz1,lelt),vr2(lx1,ly1,lz1,lelt),
     & vphi2(lx1,ly1,lz1,lelt)
      real*8 xx,yy,rr,cos,sin
	  
      ntot = nx1*ny1*nz1*nelv 
	  
      do i=1,ntot
      xx = xm1(i,1,1,1)
      yy = ym1(i,1,1,1)
      rr = sqrt(xx**2+yy**2)
      cos = xx/rr
      sin = yy/rr
      vr(i,1,1,1) = vx(i,1,1,1)*cos + vy(i,1,1,1)*sin
      vphi(i,1,1,1) = vy(i,1,1,1)*cos - vx(i,1,1,1)*sin

      vz2(i,1,1,1) =  vz(i,1,1,1)**2
      vr2(i,1,1,1) =  vr(i,1,1,1)**2
      vphi2(i,1,1,1) = vphi(i,1,1,1)**2

      enddo
	  
      call getvelocityprofile()
	  
      ! driving force for Ubar = 1
      call set_forcing(ffz_new,vz,1)


      iostep = param(15)
      if((istep.gt.0).and.(mod(istep,iostep).eq.0)) then
      call outpost(vz,vr,vphi,pr,pr,'zrp')
      call outpost(vz2,vr2,vphi2,pr,pr,'rm2')
      endif


      return
      end
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
      subroutine getvelocityprofile 
      include 'SIZE'
      include 'TOTAL'
	  
      integer icalld1
      save    icalld1
      data    icalld1 /0/
      real atime,timel
      save atime,timel

      common /vr/ vr(lx1,ly1,lz1,lelt),vphi(lx1,ly1,lz1,lelt),
     & vz2(lx1,ly1,lz1,lelt),vr2(lx1,ly1,lz1,lelt),
     & vphi2(lx1,ly1,lz1,lelt)

      parameter(nphi=18,nr=300,nz=10)
      parameter(nint=nphi*nz*nr)
      integer in,ir,iz,iphi
      real*8 xyz(3,nphi*nz*nr)
      real*8 uvw(3,nphi*nz*nr)
      real*8 uvw2(3,nphi*nz*nr)
      real*8 vr_int(nr),vz_int(nr)
      real*8 vz2_int(nr),vr2_int(nr),vphi2_int(nr)
      real*8 vr_avg(nr),vz_avg(nr),yr(nr),yr_plus(nr)
      real*8 vz2_avg(nr),vr2_avg(nr),vphi2_avg(nr)
      real*8 radius,length,alpha1,rr,pi
      integer n_interp
      real*8 nu,rho,tau_w,u_tau,delta_v,Re_tau

      n_interp =0

      radius = 0.5
      length = 5.0
      pi = 3.1415926
cc set up xyz

      call rzero(vr_int,nr)
      call rzero(vz_int,nr)
      call rzero(vr2_int,nr)
      call rzero(vz2_int,nr)
      call rzero(vphi2_int,nr)

      if(nid.eq.0) then
      n_interp = nint
      do ir=1,nr
      do iz=1,nz
      do iphi=1,nphi
      in =(ir-1)*nphi*nz + (iz-1)*nphi + iphi
      alpha1 = 2*pi*dble(iphi)/dble(nphi)
      rr =(1.0-dble(ir)/dble(nr+1))*radius
      yr(ir) = radius - rr
      xyz(1,in) = rr*cos(alpha1)
      xyz(2,in) = rr*sin(alpha1)
      xyz(3,in) = length*dble(iz)/dble(nz+1)
      enddo
      enddo
      enddo
      endif

      call interp_v(uvw,xyz,n_interp,1)
      call interp_v(uvw2,xyz,n_interp,2)

      if(nid.eq.0) then
cc return vr_avg,vz_avg,	  
      do ir=1,nr
      do iz=1,nz
      do iphi=1,nphi
      in =(ir-1)*nphi*nz + (iz-1)*nphi + iphi
      vz_int(ir) = vz_int(ir) + uvw(1,in)
      vr_int(ir) = vr_int(ir) + uvw(2,in)
      vz2_int(ir) = vz2_int(ir) + uvw2(1,in)
      vr2_int(ir) = vr2_int(ir) + uvw2(2,in)
      vphi2_int(ir) = vphi2_int(ir) + uvw2(3,in)
      enddo
      enddo
      vz_int(ir) = vz_int(ir) /(dble(nz)*dble(nphi))
      vr_int(ir) = vr_int(ir) /(dble(nz)*dble(nphi))
      vz2_int(ir) = vz2_int(ir) /(dble(nz)*dble(nphi))
      vr2_int(ir) = vr2_int(ir) /(dble(nz)*dble(nphi))
      vphi2_int(ir) = vphi2_int(ir) /(dble(nz)*dble(nphi))
      enddo

      write(6,*) "vr2_int(250):",vr2_int(250)
      write(6,*) "vr2_int(275):",vr2_int(275)
      write(6,*) "vr2_int(300):",vr2_int(300)

      write(6,*) "vphi2_int(250):",vphi2_int(250)
      write(6,*) "vphi2_int(275):",vphi2_int(275)
      write(6,*) "vphi2_int(300):",vphi2_int(300)

      if(icalld1.eq.0) then
        call rzero(vr_avg,nr)
        call rzero(vz_avg,nr)
        call rzero(vz2_avg,nr)
        call rzero(vr2_avg,nr)
        call rzero(vphi2_avg,nr)
        atime = 0.
        timel = time
        icalld1 = 1
      endif

        dtime = time - timel
        atime = atime + dtime
      if (atime.ne.0. .and. dtime.ne.0.) then
        beta      = dtime/atime
        alpha     = 1.-beta 
        do ir=1,nr
        vz_avg(ir) = vz_avg(ir)*alpha + vz_int(ir)*beta
        vr_avg(ir) = vr_avg(ir)*alpha + vr_int(ir)*beta
        vz2_avg(ir) = vz2_avg(ir)*alpha + vz2_int(ir)*beta
        vr2_avg(ir) = vr2_avg(ir)*alpha + vr2_int(ir)*beta
        vphi2_avg(ir) = vphi2_avg(ir)*alpha + vphi2_int(ir)*beta
        enddo
      endif
      timel = time

      rho = param(1)
      nu = param(2)
	  
      if(nid.eq.0) write(6,*) "rho :",rho 
      if(nid.eq.0) write(6,*) "nu :",nu 
      if(nid.eq.0) write(6,*) "vz_avg(1) :",vz_avg(1) 
      if(nid.eq.0) write(6,*) "yr(1) :",yr(1) 
	  
      tau_w = rho*nu*(vz_avg(2)-vz_avg(1))/(yr(2)-yr(1))
      u_tau = sqrt(tau_w/rho)
      delta_v = nu/u_tau
      Re_tau = radius/delta_v
	  
      if(nid.eq.0) write(6,*) "tau_w :",tau_w 
      if(nid.eq.0) write(6,*) "u_tau :",u_tau 
      if(nid.eq.0) write(6,*) "delta_v :",delta_v 
      if(nid.eq.0) write(6,*) "Re_tau :",Re_tau    
	  
	  
      iostep = param(15)
      
c dump averaged force to file at iostep
      if (mod(istep,iostep).eq.0.and.istep.gt.0.and.nid.eq.0) then
      open(unit=57,file='velocity.dat')
      write(57,'(A)')
     & '   i    yr   yplus    vz_avg    vz_rms  vr_rms  vphi_rms'
	  
   10 format(I5,',',F14.10,',',F14.10,',',F14.10,',',
     & F14.10,',',F14.10,',',F14.10)
   
      do ir = 1,nr
      write(57,10) ir,
     & yr(ir),
     & yr(ir)/delta_v,
     & vz_avg(ir)/u_tau,
     & (vz2_avg(ir)-vz_avg(ir)**2)/u_tau**2,
     & vr2_avg(ir)/u_tau**2,
     & vphi2_avg(ir)/u_tau**2
      enddo

      close(57)      

      endif
		
      endif
      return
      end
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
      subroutine interp_v(uvw,xyz,n,iflag)
c
c     evaluate velocity for list of points xyz
c
      include 'SIZE'
      include 'TOTAL'
      common /vr/ vr(lx1,ly1,lz1,lelt),vphi(lx1,ly1,lz1,lelt),
     & vz2(lx1,ly1,lz1,lelt),vr2(lx1,ly1,lz1,lelt),
     & vphi2(lx1,ly1,lz1,lelt)

      real uvw(ldim,n),xyz(ldim,n)

      parameter(nmax=60000,nfldmx=ldim) 

      common /rwk_intp/ 
     $       fwrk(lx1*ly1*lz1*lelt,nfldmx),
     $       rwk(nmax,ldim+1),
     $       fpts(ldim*nmax),
     $       pts(ldim*nmax)
      common /iwk_intp/ ih_intp,
     $       iwk(nmax,3)

      integer iflag
	 
      integer icalld,e
      save    icalld
      data    icalld /0/

      nxyz  = nx1*ny1*nz1
      ntot  = nxyz*nelt

      if (n.gt.nmax) call exitti ('n > nmax in interp_v!$',n)
      
      if (nelgt.ne.nelgv) call exitti
     $   ('nelgt.ne.nelgv not yet supported in interp_v!$',nelgv)

      do i=1,n				! ? not moving -> save?
         pts(i)     = xyz(1,i)
         pts(i + n) = xyz(2,i)
         if (if3d) pts(i + n*2) = xyz(3,i)
      enddo

      if (icalld.eq.0) then		! interpolation setup
        icalld = 1
        tolin  = 1.e-8
        !call intp_setup(tolin)
        call interp_setup(ih_intp,tolin,0,nelt)
      endif

      ! pack working array
      if(iflag.eq.1)then
      call opcopy(fwrk(1,1),fwrk(1,2),fwrk(1,3),vz,vr,vphi)
      elseif(iflag.eq.2)then
      call opcopy(fwrk(1,1),fwrk(1,2),fwrk(1,3),vz2,vr2,vphi2)	
      endif
	  
      ! interpolate
      call interp_nfld(fpts,fwrk,ndim,
     $             pts(1),pts(1+n),pts(2*n+1),n,
     $             iwk,rwk,nmax,.true.,ih_intp)
c      call intp_do(fpts,fwrk,ndim,
c     $             pts(1),pts(1+n),pts(2*n+1),n,
c     $             iwk,rwk,nmax,.true.)

      do i=1,n
         uvw(1,i) = fpts(i)
         uvw(2,i) = fpts(i + n)
         if(if3d) uvw(3,i) = fpts(i + n*2)
      enddo

      return
      end
c-----------------------------------------------------------------------
      subroutine userbc (ix,iy,iz,iside,ieg)
      include 'SIZE'
      include 'TOTAL'
      include 'NEKUSE'

      ux=0.0
      uy=0.0
      uz=0.0

      temp=0.0

      return
      end
c-----------------------------------------------------------------------
      subroutine useric (ix,iy,iz,ieg)
      include 'SIZE'
      include 'TOTAL'
      include 'NEKUSE'

      integer idum
      save    idum 
      data    idum / 0 /
      real*8 u_targ,diameter,radius

      if (idum.eq.0) idum = 99 + nid
   
      u_targ = 1.0
      diameter = 1.0
      radius = diameter/2.0	  
      r = sqrt(x*x + y*y)/radius

c     blunt profile w/ random perturbations
      eps = 0.1*u_targ
      ux  = eps*(ran1(idum)-.5)
      uy  = eps*(ran1(idum)-.5)
      uz  = 1.25*(1-r**8.0)*u_targ + eps*(ran1(idum)-.5)

      temp=0

      return
      end
c-----------------------------------------------------------------------
      subroutine usrdat   ! This routine to modify element vertices
      include 'SIZE'      ! _before_ mesh is generated, which 
      include 'TOTAL'     ! guarantees GLL mapping of mesh.

      return
      end
c-----------------------------------------------------------------------
      subroutine usrdat2   ! This routine to modify mesh coordinates
      include 'SIZE'
      include 'TOTAL'

      ntot = nx1*ny1*nz1*nelt
 
      do i=1,ntot
         xm1(i,1,1,1) = xm1(i,1,1,1)*0.5
         ym1(i,1,1,1) = ym1(i,1,1,1)*0.5
      enddo
      param(59) = 1     ! force nek5 to recognize element deformation.
   
      if (nid.eq.0) then 
      write(6,*) 'geometry scale completed'	  
      endif	  
	  
      return
      end
c-----------------------------------------------------------------------
      subroutine usrdat3
      include 'SIZE'
      include 'TOTAL'

      return
      end
c-----------------------------------------------------------------------
      subroutine set_forcing(f_new,u,idir)  ! driving force for Ubar = 1

      include 'SIZE'
      include 'TOTAL'

      common /cforce/ ffx_new,ffy_new,ffz_new
      common /uforce/ utldo(ldim)

      n=nx1*ny1*nz1*nelv

      if (istep.eq.0) f_new=param(71)  ! Update forcing for Ubar = 1

      u_targ  = 1.
      ubar    = glsc2(u,bm1,n)/volvm1
      utilde  = ubar/u_targ

      if (istep.gt.0) f_new = 0.5 * (f_new + f_new/utilde)

      if (istep.gt.5) then
         alpha = abs(utilde-utldo(idir))/dt
         alpha = min(alpha,0.05)
         if     (utilde.gt.1.and.utilde.gt.utldo(idir)) then
            f_new = (1-alpha)*f_new
         elseif (utilde.lt.1.and.utilde.lt.utldo(idir)) then
            f_new = (1+alpha)*f_new
         endif
      endif
      utldo(idir)   = utilde

      f_min = 0.00001                 ! ad hoc limits
      f_new = max(f_new,f_min)
      
      f_max = 0.10000                 ! ad hoc limits
      f_new = min(f_new,f_max)

      if(nid.eq.0) write(6,*) "new force :",f_new     
 
      return
      end
c-----------------------------------------------------------------------
c automatically added by makenek
      subroutine usrsetvert(glo_num,nel,nx,ny,nz) ! to modify glo_num
      integer*8 glo_num(1)
      return
      end
