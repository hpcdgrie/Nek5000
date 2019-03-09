C-----------------------------------------------------------------------
      subroutine fluid_im (igeom)
C
C     Driver for solving the incompressible Navier-Stokes equations.
C
C     Current version:
C     (1) Velocity/stress formulation.
C     (2) Constant/variable properties.
C     (3) Implicit/explicit time stepping.
C     (4) Automatic setting of tolerances .
C     (5) Lagrangian/"Eulerian"(operator splitting) modes
C
C-----------------------------------------------------------------------
      include 'SIZE'
      include 'DEALIAS'
      include 'INPUT'
      include 'SOLN'
      include 'TSTEP'

      real*8 ts, dnekclock
 
      ifield = 1
      imesh  = 1
      call unorm
      call settolv

      ts = dnekclock() 

      if(nio.eq.0 .and. igeom.ge.2)
     &   write(*,'(13x,a)') 'Solving for fluid'

      ! if (ifsplit) then

c        PLAN 4: TOMBO SPLITTING
c                - Time-dependent Navier-Stokes calculation (Re>>1).
c                - Same approximation spaces for pressure and velocity.
c                - Incompressibe or Weakly compressible (div u .ne. 0).

         call plan4_im (igeom)                                           
         if (igeom.ge.2) call chkptol         ! check pressure tolerance 
         if (igeom.ge.2) call printdiverr

      if(nio.eq.0 .and. igeom.ge.2) 
     &   write(*,'(4x,i7,a,1p2e12.4)') 
     &   istep,'  Fluid done',time,dnekclock()-ts

      return
      end
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
      subroutine plan4_im (igeom)

C     Splitting scheme A.G. Tomboulides et al.
c     Journal of Sci.Comp.,Vol. 12, No. 2, 1998
c
C     NOTE: QTL denotes the so called thermal
c           divergence and has to be provided
c           by userqtl.
c
      INCLUDE 'SIZE'
      INCLUDE 'INPUT'
      INCLUDE 'GEOM'
      INCLUDE 'MASS'
      INCLUDE 'SOLN'
      INCLUDE 'MVGEOM'
      INCLUDE 'TSTEP'
      INCLUDE 'ORTHOP'
      INCLUDE 'CTIMER'
#     include "IMPLICIT"
C
      COMMON /SCRNS/ RES1  (LX1,LY1,LZ1,LELV)
     $ ,             RES2  (LX1,LY1,LZ1,LELV)
     $ ,             RES3  (LX1,LY1,LZ1,LELV)
     $ ,             DV1   (LX1,LY1,LZ1,LELV)
     $ ,             DV2   (LX1,LY1,LZ1,LELV)
     $ ,             DV3   (LX1,LY1,LZ1,LELV)
     $ ,             RESPR (LX2,LY2,LZ2,LELV)
      common /scrvh/ h1    (lx1,ly1,lz1,lelv)
     $ ,             h2    (lx1,ly1,lz1,lelv)
 
      REAL           DPR   (LX2,LY2,LZ2,LELV)
      EQUIVALENCE   (DPR,DV1)
      LOGICAL        IFSTSP

      REAL DVC (LX1,LY1,LZ1,LELV), DFC(LX1,LY1,LZ1,LELV)
      REAL DIV1, DIV2, DIF1, DIF2, QTL1, QTL2
c
      INTYPE = -1
      NTOT1  = lx1*ly1*lz1*NELV

c      call copy(vx_e,vx,ntot1)
c      call copy(vy_e,vy,ntot1)
c      call copy(vz_e,vz,ntot1)
	  
      if (igeom.eq.1) then

         ! compute explicit contributions bfx,bfy,bfz 
         call makef_im(igeom)
		 ! call makef()

         call sumab(vx_e,vx,vxlag,ntot1,ab,nab)     ! vx_e(extrapolated velocities are used in pressure eqn)
         call sumab(vy_e,vy,vylag,ntot1,ab,nab)     ! (in crespsp, Compute startresidual/right-hand-side in the pressure)
         if (if3d) call sumab(vz_e,vz,vzlag,ntot1,ab,nab)

      else
         if((ngeom.gt.2).and.(igeom.gt.2)) call makef_im(igeom)
		 if((ngeom.gt.2).and.(igeom.gt.2)) then     ! use update vx_e with current vx
            call copy(vx_e,vx,ntot1)
            call copy(vy_e,vy,ntot1)
            call copy(vz_e,vz,ntot1)
         endif
         if(iflomach) call opcolv(bfx,bfy,bfz,vtrans)

         ! add user defined divergence to qtl 
         call add2 (qtl,usrdiv,ntot1)
         if (igeom.eq.2) call lagvel
         ! mask Dirichlet boundaries
         call bcdirvc  (vx,vy,vz,v1mask,v2mask,v3mask) 

         ! compute pressure
         call copy(prlag,pr,ntot1)
         if (icalld.eq.0) tpres=0.0
         icalld=icalld+1
         npres=icalld
         etime1=dnekclock()

         call crespsp  (respr)
         call invers2  (h1,vtrans,ntot1)
         call rzero    (h2,ntot1)
         call ctolspl  (tolspl,respr)
         napproxp(1) = laxtp
         call hsolve   ('PRES',dpr,respr,h1,h2 
     $                        ,pmask,vmult
     $                        ,imesh,tolspl,nmxp,1
     $                        ,approxp,napproxp,binvm1)
         call add2    (pr,dpr,ntot1)
         call ortho   (pr)

         tpres=tpres+(dnekclock()-etime1)

         ! compute velocity
         if(ifstrs .and. .not.ifaxis) then
            call bcneutr
            call cresvsp_weak(res1,res2,res3,h1,h2)
         else
            call cresvsp     (res1,res2,res3,h1,h2)
         endif
         call ophinv       (dv1,dv2,dv3,res1,res2,res3,h1,h2,tolhv,nmxv)
         call opadd2       (vx,vy,vz,dv1,dv2,dv3)

      endif

      return
      END

c-----------------------------------------------------------------------
      subroutine makef_im(igeom)
C---------------------------------------------------------------------
C
C     Compute and add: (1) user specified forcing function (FX,FY,FZ)
C                      (2) driving force due to natural convection
C                      (3) convection term
C
C     !! NOTE: Do not change the arrays BFX, BFY, BFZ until the
C              current time step is completed.
C
C----------------------------------------------------------------------
      include 'SIZE'
      include 'SOLN'
      include 'MASS'
      include 'INPUT'
      include 'TSTEP'
      include 'CTIMER'
      include 'MVGEOM'
#     include "IMPLICIT"

      etime1 = dnekclock()
	  
      if(ifchar.and.(nid.eq.0)) then
      write(6,*) 'ERROR: cannot use characteristics'
      endif
	  
      call makeuf            ! define bfx/y/z = usr define force
      ! if (filterType.eq.2)   call make_hpf  ! add high pass filter to bfx/y/z
      call makevis           ! add bfx/y/z with explicit part of stress term
      call advab_im          ! add bfx/y/z with convection term
      call makeabf_im(igeom)        ! use extrapolation to obtain stable bfx/y/z
                             ! however, this extrapolation is turned off now.
      if(igeom.eq.1) call makebdf_im
      CALL OPADD2(bfx,bfy,bfz,bfx_tlag,bfy_tlag,bfz_tlag)

      if(igeom.eq.1) tmakf=tmakf+(dnekclock()-etime1)

c      ntot = lx1*ly1*lz1*NELV
c      var_max = glmax(bfx_tlag,ntot)
c      var_min = glmin(bfx_tlag,ntot)
c      if(nio.eq.0) write(6,*) "bfx_tlag min/max",var_min,"-",var_max 
	  
c      var_max = glmax(bfx,ntot)
c      var_min = glmin(bfx,ntot)
c      if(nio.eq.0) write(6,*) "bfx min/max",var_min,"-",var_max 
	  
	  
      return
      end
c---------------------------------------------------------------
      subroutine advab_im
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
C
      COMMON /SCRUZ/ TA1 (LX1,LY1,LZ1,LELV)
     $ ,             TA2 (LX1,LY1,LZ1,LELV)
     $ ,             TA3 (LX1,LY1,LZ1,LELV)
C
      NTOT1 = lx1*ly1*lz1*NELV
      CALL CONVOP  (TA1,VX) 
      CALL CONVOP  (TA2,VY)
      CALL SUBCOL3 (BFX,BM1,TA1,NTOT1)  ! BFX = BFX  - BM1*TA1 
      CALL SUBCOL3 (BFY,BM1,TA2,NTOT1)
      IF (ldim.EQ.2) THEN
         CALL RZERO (TA3,NTOT1)
      ELSE
         CALL CONVOP  (TA3,VZ)
         CALL SUBCOL3 (BFZ,BM1,TA3,NTOT1)
      ENDIF
C

      return
      END
c-----------------------------------------------------------------------
      subroutine makeabf_im(igeom) 
C-----------------------------------------------------------------------
C
C     Sum up contributions to kth order extrapolation scheme.
C
C-----------------------------------------------------------------------
      include 'SIZE'
      include 'INPUT'
      include 'SOLN'
      include 'TSTEP'
C
      COMMON /SCRUZ/ TA1 (LX1,LY1,LZ1,LELV)
     $ ,             TA2 (LX1,LY1,LZ1,LELV)
     $ ,             TA3 (LX1,LY1,LZ1,LELV)
C
      NTOT1 = lx1*ly1*lz1*NELV
C
      AB0 = AB(1)
      AB1 = AB(2)
      AB2 = AB(3)
	  
c use extrapolation
c only do extrapolation at igeom=1
      if(igeom.eq.1) then
      CALL ADD3S2 (TA1,ABX1,ABX2,AB1,AB2,NTOT1)  ! TA1 = ABX1*AB1 + ABX2*AB2
      CALL ADD3S2 (TA2,ABY1,ABY2,AB1,AB2,NTOT1)
      CALL COPY   (ABX2,ABX1,NTOT1)
      CALL COPY   (ABY2,ABY1,NTOT1)
      CALL COPY   (ABX1,BFX,NTOT1)
      CALL COPY   (ABY1,BFY,NTOT1)
      CALL ADD2S1 (BFX,TA1,AB0,NTOT1)            ! BFX = BFX*AB0 + TA1 => BFX = BFX*AB0 +ABX1*AB1 + ABX2*AB2
      CALL ADD2S1 (BFY,TA2,AB0,NTOT1)
      IF (ldim.EQ.3) THEN
         CALL ADD3S2 (TA3,ABZ1,ABZ2,AB1,AB2,NTOT1)
         CALL COPY   (ABZ2,ABZ1,NTOT1)
         CALL COPY   (ABZ1,BFZ,NTOT1)
         CALL ADD2S1 (BFZ,TA3,AB0,NTOT1)
      ENDIF
      endif
C
c turn off extrapolation now. 
c

      if(.not.iflomach) CALL COL2   (BFX,VTRANS,NTOT1)          ! multiply by density
      if(.not.iflomach) CALL COL2   (BFY,VTRANS,NTOT1)
      IF (ldim.EQ.3) THEN
         if(.not.iflomach) CALL COL2   (BFZ,VTRANS,NTOT1)
      ENDIF

      return
      END
C
c-----------------------------------------------------------------------
      subroutine makebdf_im
C
C     Add contributions to F from lagged BD terms.
C
      include 'SIZE'
      include 'SOLN'
      include 'MASS'
      include 'GEOM'
      include 'INPUT'
      include 'TSTEP'
#     include "IMPLICIT"

      COMMON /SCRNS/ TA1(LX1,LY1,LZ1,LELV)
     $ ,             TA2(LX1,LY1,LZ1,LELV)
     $ ,             TA3(LX1,LY1,LZ1,LELV)
     $ ,             TB1(LX1,LY1,LZ1,LELV)
     $ ,             TB2(LX1,LY1,LZ1,LELV)
     $ ,             TB3(LX1,LY1,LZ1,LELV)
     $ ,             H2 (LX1,LY1,LZ1,LELV)
C
      NTOT1 = lx1*ly1*lz1*NELV
      CONST = 1./DT

      if(iflomach) then
        call cfill(h2,CONST,ntot1)
      else
        call cmult2(h2,vtrans(1,1,1,1,ifield),const,ntot1)
      endif

      CALL OPCOLV3c (TB1,TB2,TB3,vx,vy,vz,BM1,bd(2))
C
      DO 100 ILAG=2,NBD
         IF (IFGEOM) THEN
            CALL OPCOLV3c(TA1,TA2,TA3,VXLAG (1,1,1,1,ILAG-1),
     $                                VYLAG (1,1,1,1,ILAG-1),
     $                                VZLAG (1,1,1,1,ILAG-1),
     $                                BM1LAG(1,1,1,1,ILAG-1),bd(ilag+1))
         ELSE
            CALL OPCOLV3c(TA1,TA2,TA3,VXLAG (1,1,1,1,ILAG-1),
     $                                VYLAG (1,1,1,1,ILAG-1),
     $                                VZLAG (1,1,1,1,ILAG-1),
     $                                BM1                   ,bd(ilag+1))
         ENDIF
         CALL OPADD2  (TB1,TB2,TB3,TA1,TA2,TA3)
 100  CONTINUE
c       CALL OPADD2col (BFX,BFY,BFZ,TB1,TB2,TB3,h2)
	  
      call rzero(bfx_tlag,NTOT1)
      call rzero(bfy_tlag,NTOT1)
      call rzero(bfz_tlag,NTOT1)
      CALL OPADD2col (bfx_tlag,bfy_tlag,bfz_tlag,TB1,TB2,TB3,h2)
	  
C
      return
      END
c-----------------------------------------------------------------------