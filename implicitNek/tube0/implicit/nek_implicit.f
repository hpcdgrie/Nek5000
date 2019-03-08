c Created by Haomin Yuan, Argonne National Lab, 2/21/2019
c 
C  The idea is to update convection term to make it pseudo implicit
c  so we sacrfise speed to improve accuracy.
C 
c   nek_advance_implicit should replace nek_advance
c   1. ngeom shuold be greate than 2 to apply mutiple iterations within
c      one time-step
c   2. explicit time term is only updated at igeom =1
c   3. in old solver, convection term is only updated at igeom =1,
c      but now convection term is updated at every igeom.
c   4. at igeom =1,2, convection term is obtained from extrapolation. 
c      but at igeom =3 or greater, convection term is updated using current velocity
c
c
c-----------------------------------------------------------------------
#include "fluid_im.f"
#include "heat_im.f"
c-----------------------------------------------------------------------
      subroutine nek_advance_implicit

      include 'SIZE'
      include 'TOTAL'
      include 'CTIMER'
#     include "IMPLICIT"
	 
      common /cgeom/ igeom

      ntot = lx1*ly1*lz1*nelv

      call nekgsync

      call setup_convect(2) ! Save conv vel

      if (iftran) call settime
      if (ifmhd ) call cfl_check
      call setsolv
      call comment

      ngeom = 4 ! or usr define in usrdat
	  
      if(ngeom.eq.2) then	  
        if(nid.eq.0) then
      write(6,*) "NOTE: ngeom should greater than 2 "
      write(6,*) "for implicit solver"
      write(6,*) "reset ngeom to 3"
        endif
        ngeom = 3
      endif
	  
      do igeom=1,ngeom

         if (ifgeom) then
            call gengeom (igeom)
            call geneig  (igeom)
         endif

         if (igeom.ge.2) then  
            call setprop
            call rzero(qtl,ntot)
            if (iflomach) call qthermal
         endif

         call setup_convect(2)  ! reset convection velocity

         if (ifheat)          call heat_im (igeom)
         if (ifflow)          call fluid_im (igeom)
         !if (ifflow)          call fluid(igeom)
         if (igeom.eq.ngeom.and.filterType.eq.1)
     $                        call q_filter(param(103))

      enddo

      return
      end


