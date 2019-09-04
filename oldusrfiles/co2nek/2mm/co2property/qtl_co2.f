c-----------------------------------------------------------------------
      subroutine qtl_co2()
c calculate thermal divergence based on real gas property.
c    Compute the thermal divergence QTL 
C
C     QTL := div(v) = -1/rho * Drho/Dt
c
c     If we use the ideal gas law and assume
c     that p,R is const we end up with
c     QTL = 1/(rho*cp) rho*cp*DT/Dt
C
C     where rho*cp*DT/Dt represents the RHS of the
C     energy equation expressed in terms of temperature.
c 
c 
c==============================================================================
c Haomin Yuan:
c    modified thermal divergence term for non-ideal gas assumption.
c    using enthalpy for energy equation.
c    
c    assuming pressure is constant
c    div(v) = -1/rho * Drho/Dt = -1/rho *(drho/dh)* Dh/Dt = -1/rho**2 *(drho/dh)*(rho * Dh/Dt)
c
c    assume -(drho/dh) = beta
c
c     div(v) = beta/rho**2 *(rho * Dh/Dt)
c     div(v) = beta/rho**2 * h_RHS
c    
c    (rho * Dh/Dt) represents the RHS of the energy equation expressed in terms of enthalpy.
c
c    (rho * Dh/Dt) = delta ((k/cp)delta(h))
c  
c===============================================================================

      include 'SIZE'
      include 'TOTAL'
#     include "CO2PROP"

      COMMON /SCRNS/ w1(LX1,LY1,LZ1,LELT)
     $              ,w2(LX1,LY1,LZ1,LELT)
     $              ,tx(LX1,LY1,LZ1,LELT)
     $              ,ty(LX1,LY1,LZ1,LELT)
     $              ,tz(LX1,LY1,LZ1,LELT)
      real varmin,varmax
      ntot = nx1*ny1*nz1*nelv

      if (.not.iflomach) then
         call rzero(qtl,ntot)
         return
      endif

      ifld_save = ifield

c - - Assemble RHS of T-eqn
      ifield=2
      call setqvol (QTL) ! volumetric heating source
      call col2    (QTL,BM1,ntot)

      ifield=1     !set right gs handle (QTL is only defined on the velocity mesh)
      call opgrad  (tx,ty,tz,t(1,1,1,1,1))
      call opdssum (tx,ty,tz)
      call opcolv  (tx,ty,tz,binvm1)
      call opcolv  (tx,ty,tz,vdiff(1,1,1,1,2))    
      ! vdiff2 is thermal conductivity (kappa) in temperature equation, but is kappa/Cp in enthalpy equation
      call opdiv   (w2,tx,ty,tz)

      call add2    (QTL,w2,ntot)
      call col2    (QTL,beta_co2,ntot)
      call col3    (w1,rho_co2,rho_co2,ntot) 
      call invcol2 (QTL,w1,ntot)
      
      call dssum   (QTL,nx1,ny1,nz1)
      call col2    (QTL,binvm1,ntot)

      varmax=glmax(QTL,ntot)  
      varmin=glmin(QTL,ntot)
      if (nid.eq.0) write(6,*) "QTL: ",varmin," - ",varmax

      ifield = ifld_save

      return
      end
c-----------------------------------------------------------------------