c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
      subroutine heat_im (igeom)
C
C     Driver for temperature or passive scalar.
C
C     Current version:
C     (1) Varaiable properties.
C     (2) Implicit time stepping.
C     (3) User specified tolerance for the Helmholtz solver
C         (not based on eigenvalues).
C     (4) A passive scalar can be defined on either the 
C         temperatur or the velocity mesh.
C     (5) A passive scalar has its own multiplicity (B.C.).  
C
      include 'SIZE'
      include 'INPUT'
      include 'TSTEP'
      include 'DEALIAS'

      real*8 ts, dnekclock

      ts = dnekclock()

      if (nio.eq.0 .and. igeom.ge.2) 
     &    write(*,'(13x,a)') 'Solving for Hmholtz scalars'

      do ifield = 2,nfield
         if (idpss(ifield-1).eq.0) then      ! helmholtz
            intype        = -1
            if (.not.iftmsh(ifield)) imesh = 1
            if (     iftmsh(ifield)) imesh = 2
            call unorm
            call settolt
            call cdscal_im(igeom)
         endif
      enddo

      if (nio.eq.0 .and. igeom.ge.2)
     &   write(*,'(4x,i7,a,1p2e12.4)') 
     &   istep,'  Scalars done',time,dnekclock()-ts

      return
      end
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
      subroutine cdscal_im (igeom)
C
C     Solve the convection-diffusion equation for passive scalar IPSCAL
C
      include 'SIZE'
      include 'INPUT'
      include 'GEOM'
      include 'MVGEOM'
      include 'SOLN'
      include 'MASS'
      include 'TSTEP'
      COMMON  /CPRINT/ IFPRINT
      LOGICAL          IFPRINT
      LOGICAL          IFCONV

      COMMON /SCRNS/ TA(LX1,LY1,LZ1,LELT)
     $              ,TB(LX1,LY1,LZ1,LELT)
      COMMON /SCRVH/ H1(LX1,LY1,LZ1,LELT)
     $              ,H2(LX1,LY1,LZ1,LELT)




      include 'ORTHOT'

      if (ifdgfld(ifield)) then
         call cdscal_dg(igeom)
         return
      endif


      ifld1 = ifield-1
      napproxt(1,ifld1) = laxtt

      nel    = nelfld(ifield)
      n   = lx1*ly1*lz1*nel

      if (igeom.eq.1) then   ! geometry at t^{n-1}

         call makeq_im(igeom)
         call lagscal

      else                   ! geometry at t^n
         if((ngeom.gt.2).and.(igeom.gt.2)) call makeq_im(igeom)

         IF (IFPRINT) THEN
            IF (IFIELD.EQ.2.AND.NID.EQ.0)
     $          WRITE (6,*) ' Temperature/Passive scalar solution'
         ENDIF

         if1=ifield-1
         write(name4t,1) if1-1
    1    format('PS',i2)
         if(ifield.eq.2) write(name4t,'(A4)') 'TEMP'

C
C        New geometry
C
         isd = 1
         if (ifaxis.and.ifaziv.and.ifield.eq.2) isd = 2
c        if (ifaxis.and.ifmhd) isd = 2 !This is a problem if T is to be T!

         do 1000 iter=1,nmxnl ! iterate for nonlin. prob. (e.g. radiation b.c.)

         intype = 0
         if (iftran) intype = -1
         call sethlm  (h1,h2,intype)
         call bcneusc (ta,-1)
         call add2    (h2,ta,n)
         call bcdirsc (t(1,1,1,1,ifield-1))
         call axhelm  (ta,t(1,1,1,1,ifield-1),h1,h2,imesh,ISD)
         call sub3    (tb,bq(1,1,1,1,ifield-1),ta,n)
         call bcneusc (ta,1)
         call add2    (tb,ta,n)

         if(iftmsh(ifield)) then
           call hsolve  (name4t,TA,TB,H1,H2 
     $                   ,tmask(1,1,1,1,ifield-1)
     $                   ,tmult(1,1,1,1,ifield-1)
     $                   ,imesh,tolht(ifield),nmxt(ifield-1),1
     $                   ,approxt(1,0,ifld1),napproxt(1,ifld1),bintm1)
         else
           call hsolve  (name4t,TA,TB,H1,H2 
     $                   ,tmask(1,1,1,1,ifield-1)
     $                   ,tmult(1,1,1,1,ifield-1)
     $                   ,imesh,tolht(ifield),nmxt(ifield-1),1
     $                   ,approxt(1,0,ifld1),napproxt(1,ifld1),binvm1)
         endif 

         call add2    (t(1,1,1,1,ifield-1),ta,n)

         call cvgnlps (ifconv) ! Check convergence for nonlinear problem 
         if (ifconv) goto 2000

C        Radiation case, smooth convergence, avoid flip-flop (ER).
         call cmult (ta,0.5,n)
         call sub2  (t(1,1,1,1,ifield-1),ta,n)

 1000    continue
 2000    continue

      endif

      return
      end

c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
      subroutine makeq_im(igeom)
c construct explicit term for energy eqn
c convective term + partial time marching term. 
c
C     Generate forcing function for the solution of a passive scalar.
C     !! NOTE: Do not change the content of the array BQ until the current

      include 'SIZE'
      include 'TOTAL'
      include 'CTIMER'
#     include "IMPLICIT"

      logical  if_conv_std
      common /SCRUZ/ w1(lx1,ly1,lz1,lelt)

      nxyz = lx1*ly1*lz1
      ntot = nxyz*nelv

      etime0 = dnekclock()      

      if (nio.eq.0.and.loglevel.gt.2)
     $   write(6,*) 'makeq', ifield

      if_conv_std = .true.

      call makeq_aux                             ! add user q, bq = usrq or zero.
      call convab_im                             ! calculate convective term , add to bq
      call makeabq_im(igeom)                      ! use 3rd order Adams-Bashforth scheme to obtain stable convective term. bq is overwrite here.
      if(igeom.eq.1) call makebdq_im             ! partial time marching term (bq_tlag, is only calculated at igeom =1, after that no change)
      call ADD2(bq(1,1,1,1,ifield-1),bq_tlag(1,1,1,1,ifield-1),ntot)
	  
      if(igeom.eq.1) tmakq=tmakq+(dnekclock()-etime0)

      var_max = glmax(bq_tlag(1,1,1,1,ifield-1),ntot)
      var_min = glmin(bq_tlag(1,1,1,1,ifield-1),ntot)
      if(nio.eq.0) write(6,*) "bq_tlag ",ifield," min/max",var_min,"-",var_max 
	  
      var_max = glmax(bq(1,1,1,1,ifield-1),ntot)
      var_min = glmin(bq(1,1,1,1,ifield-1),ntot)
      if(nio.eq.0) write(6,*) "bq ",ifield," min/max",var_min,"-",var_max 
	  
      return
      end
c-----------------------------------------------------------------------
      subroutine convab_im
C---------------------------------------------------------------
C
C     Eulerian scheme, add convection term to forcing function 
C     at current time step.
C
C---------------------------------------------------------------
      include 'SIZE'
      include 'SOLN'
      include 'MASS'
      include 'TSTEP'

      common /scruz/ ta (lx1*ly1*lz1*lelt)

      nel = nelfld(ifield)
      n   = lx1*ly1*lz1*nel

      call convop  (ta,t(1,1,1,1,ifield-1))
      do i=1,n
        bq(i,1,1,1,ifield-1) = bq (i,1,1,1,ifield-1)
     $                       - bm1(i,1,1,1)*ta(i)*vtrans(i,1,1,1,ifield)
      enddo

      return
      end
c-----------------------------------------------------------------------
      subroutine makeabq_im(igeom) 
c convective term.
c fully explicit term now, extrapolated from previous time-steps
c should be implicit or partially implicit 
c
C
C     Sum up contributions to 3rd order Adams-Bashforth scheme.
C
      include 'SIZE'
      include 'SOLN'
      include 'TSTEP'

      ab0   = ab(1)
      ab1   = ab(2)
      ab2   = ab(3)
      nel   = nelfld(ifield)
      n     = lx1*ly1*lz1*nel

c use extrapolation
c only do extrapolation at igeom=1
      if(igeom.eq.1) then
      do i=1,n
         ta=ab1*vgradt1(i,1,1,1,ifield-1)+ab2*vgradt2(i,1,1,1,ifield-1)
         vgradt2(i,1,1,1,ifield-1)=vgradt1(i,1,1,1,ifield-1)
         vgradt1(i,1,1,1,ifield-1)=bq     (i,1,1,1,ifield-1)
         bq     (i,1,1,1,ifield-1)=bq     (i,1,1,1,ifield-1)*ab0+ta
      enddo
      endif

      return
      end
c-----------------------------------------------------------------------
      subroutine makebdq_im
C-----------------------------------------------------------------------
C
C     Add contributions to F from lagged BD terms.
C
C-----------------------------------------------------------------------
      include 'SIZE'
      include 'TOTAL'
#     include "IMPLICIT"

      parameter (lt=lx1*ly1*lz1*lelt)
      common /scrns/ tb(lt),h2(lt)

      nel   = nelfld(ifield)
      n     = lx1*ly1*lz1*nel
   
      const = 1./dt
      do i=1,n
         h2(i)=const*vtrans(i,1,1,1,ifield)
         tb(i)=bd(2)*bm1(i,1,1,1)*t(i,1,1,1,ifield-1)
      enddo

      do ilag=2,nbd
         if (ifgeom) then
            do i=1,n
               ta=bm1lag(i,1,1,1,ilag-1)*tlag(i,1,1,1,ilag-1,ifield-1)
               tb(i)=tb(i)+ta*bd(ilag+1)
            enddo
         else
            do i=1,n
               ta=bm1(i,1,1,1)*tlag(i,1,1,1,ilag-1,ifield-1)
               tb(i)=tb(i)+ta*bd(ilag+1)
            enddo
         endif
      enddo

c      call addcol3 (bq(1,1,1,1,ifield-1),tb,h2,n)
      call rzero(bq_tlag(1,1,1,1,ifield-1),n)
      call addcol3(bq_tlag(1,1,1,1,ifield-1),tb,h2,n)

      return
      end
c-----------------------------------------------------------------------