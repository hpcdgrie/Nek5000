c gmsh2nek
c convert gmsh .msh file (version 2, ascii format) to re2 file
c 
c trying to add periodicity setting when doing mesh converting.
c and set boundary condition direct in 
c casename.bc 
c Haomin Yuan. 9/24/18
c
c update: 
c 1.
c modified the boundary set up part. 
c old scheme is very slow when element number gets large. 
c old scheme is O(N^2)
c new scheme is much faster.
c new scheme store the information of node - > hex.
c thus, for each quad, only need to search nearby hex to setup boundaries.
c new scheme is O(N)
c
c 2.
c optimize memory usage of arrays. in msh file, node number is greatly 
c larger than quad and hex numbers. its array size for nodes is now
c 10 times of array size of quad and hex.
c
c Haomin Yuan. 3/14/2019
c
c-----------------------------------------------------------------------
      program gmsh2nek
#     include "SIZE"

      call read_input_name  ! read input name
      call gmsh_read        ! read gmsh .msh file (nodes,quads,hexes)
      call convert          ! convert gmsh hex20 to nek hex20
      call setbc            ! set boundary condition
      call gen_re2          ! write nek mesh re2 file 

      end 
c-----------------------------------------------------------------------
      subroutine read_input_name

#     include "SIZE"

      character*1 re2nam1(80)
      character*1 mshnam1(32)

      equivalence(re2name,re2nam1)
      equivalence(mshname,mshnam1) 

      call blank (mshname, 32)
      call blank (re2name, 80)

      write(6,*) 'Input (.msh) file name:'
      read (5,'(a32)') mshname
      len = ltrunc(mshname,32)

      call chcopy(re2name        ,mshname,32)
      call chcopy(mshnam1(len+1) ,'.msh' , 4)
      call chcopy(re2nam1(len+1) ,'.re2' , 4)

      return 
      end
c-----------------------------------------------------------------------
      subroutine gmsh_read
c read .msh file (version 2, ascii format)
#     include "SIZE"

      character*32  mshnam2
      character*1   mshnam3(32)
      character*80 charline
      character*1  charlin1(80)
      integer A,B,C,elemType

      equivalence(mshnam2,mshnam3)
      equivalence(charline,charlin1) 	  

      call chcopy(mshnam2,mshname,32)	  
      len = ltrunc(mshnam2,32)
      call chcopy(mshnam3(len+1) ,'_1' , 2)


      call blank (charline,80)
      call chcopy(charlin1(1),'cp ',3)
	
      len = ltrunc(mshname,32)
      call chcopy(charlin1(4),mshname,len)
	
      len = ltrunc(charline,80)
      call chcopy(charlin1(len+1),' ',1)

      len2 = ltrunc(mshnam2,32)
      call chcopy(charlin1(len+2),mshnam2,len2)

c      write(6,*) charline
	  
      call system(charline)	  

      call blank (charline,80)

      open(299,file=mshname)
      open(300,file=mshnam2)

      write(6,*) 'Starting reading ',mshname

c loop to find Nodes section
      do while (.true.) 
        read(299,*) charline
        read(300,*)
        charline = trim(charline)
        if (charline.eq."$Nodes") goto 1010
      enddo
c end loop to "$Nodes"

c read all nodes xyz
1010  read(299,*) totalNode
      read(300,*)

      if(totalNode.gt.(max_num_elem*10)) then
       write(6,*) 'ERROR, increase MAXNEL to ',totalNode/10+1
       write(6,*) 'and recompile gmsh2nek'
       call exit(1)
      endif

      do inode = 1,totalNode
      read(299,*)A,node_xyz(1,inode),node_xyz(2,inode)
     &,node_xyz(3,inode)
      read(300,*)
      enddo
c end read all nodes xyz

      call rzero(node_quad,4*max_num_elem*10) ! node -> quad
      call rzero(node_hex,8*max_num_elem*10)  ! node -> hex

      read(299,*)charline ! "$EndNodes"
      read(300,*)
      read(299,*)charline ! "$Elements"
      read(300,*)
	  
      read(299,*)totalElem
      read(300,*)
      totalQuad = 0
      totalHex = 0

      do iElem= 1,totalElem
      read(299,*) A,elemType

	  ! detemine element type
      if (elemType.eq.16) then ! if quad8
      totalQuad = totalQuad + 1
      read(300,*) A,B,C,
     &quad_array(1,totalQuad),quad_array(2,totalQuad),
     &quad_array(3,totalQuad),quad_array(4,totalQuad),
     &quad_array(5,totalQuad),quad_array(6,totalQuad),
     &quad_array(7,totalQuad),quad_array(8,totalQuad),
     &quad_array(9,totalQuad),quad_array(10,totalQuad)

      do inode = 1,4
      call addTo_node_quad(quad_array(2+inode,totalQuad),totalQuad)
      enddo
 
      elseif (elemType.eq.10) then ! if quad9
      totalQuad = totalQuad + 1
      read(300,*) A,B,C,
     &quad_array(1,totalQuad),quad_array(2,totalQuad),
     &quad_array(3,totalQuad),quad_array(4,totalQuad),
     &quad_array(5,totalQuad),quad_array(6,totalQuad),
     &quad_array(7,totalQuad),quad_array(8,totalQuad),
     &quad_array(9,totalQuad),quad_array(10,totalQuad),
     &quad_array(11,totalQuad)

      do inode = 1,4
      call addTo_node_quad(quad_array(2+inode,totalQuad),totalQuad)
      enddo 
	 
      elseif (elemType.eq.17) then ! if hex20
      totalHex = totalHex + 1
      read(300,*) A,B,C,
     &hex_array(1,totalHex),hex_array(2,totalHex),
     &hex_array(3,totalHex),hex_array(4,totalHex),
     &hex_array(5,totalHex),hex_array(6,totalHex),
     &hex_array(7,totalHex),hex_array(8,totalHex),
     &hex_array(9,totalHex),hex_array(10,totalHex),
     &hex_array(11,totalHex),hex_array(12,totalHex),
     &hex_array(13,totalHex),hex_array(14,totalHex),
     &hex_array(15,totalHex),hex_array(16,totalHex),
     &hex_array(17,totalHex),hex_array(18,totalHex),
     &hex_array(19,totalHex),hex_array(20,totalHex),
     &hex_array(21,totalHex),hex_array(22,totalHex)
	 
      do inode = 1,8
      call addTo_node_hex(hex_array(2+inode,totalHex),totalHex)
      enddo

      elseif (elemType.eq.12) then ! if hex27
      totalHex = totalHex + 1
      read(300,*) A,B,C,
     &hex_array(1,totalHex),hex_array(2,totalHex),
     &hex_array(3,totalHex),hex_array(4,totalHex),
     &hex_array(5,totalHex),hex_array(6,totalHex),
     &hex_array(7,totalHex),hex_array(8,totalHex),
     &hex_array(9,totalHex),hex_array(10,totalHex),
     &hex_array(11,totalHex),hex_array(12,totalHex),
     &hex_array(13,totalHex),hex_array(14,totalHex),
     &hex_array(15,totalHex),hex_array(16,totalHex),
     &hex_array(17,totalHex),hex_array(18,totalHex),
     &hex_array(19,totalHex),hex_array(20,totalHex),
     &hex_array(21,totalHex),hex_array(22,totalHex),
     &hex_array(23,totalHex),hex_array(24,totalHex),
     &hex_array(25,totalHex),hex_array(26,totalHex),
     &hex_array(27,totalHex),hex_array(28,totalHex),
     &hex_array(29,totalHex)

      do inode = 1,8
      call addTo_node_hex(hex_array(2+inode,totalHex),totalHex)
      enddo

      else 
c proceed to next line in file handle 300.
      read(300,*)
      endif

      enddo

      if(totalQuad.gt.max_num_elem) then
       write(6,*) 'ERROR, increase MAXNEL to ',totalQuad
       write(6,*) 'and recompile gmsh2nek'
       call exit(1)
      endif

      if(totalHex.gt.max_num_elem) then
       write(6,*) 'ERROR, increase MAXNEL to ',totalHex
       write(6,*) 'and recompile gmsh2nek'
       call exit(1)
      endif

      write (6,*) 'total node number is ', totalNode
      write (6,*) 'total quad element number is ', totalQuad
      write (6,*) 'total hex element number is ', totalHex

      close(299)
      close(300)

      write(6,*) 'Done:: reading ',mshname  

      num_dim = 3
      num_elem = totalHex
 
      call blank (charline,80)
      call chcopy(charlin1(1),'rm ',3)
	
      len = ltrunc(mshnam2,32)
      call chcopy(charlin1(4),mshnam2,len)
	
      call system(charline)

      return
      end
C-----------------------------------------------------------------
      subroutine addTo_node_quad(inode,iquad)
#     include "SIZE"
      integer inode,iquad
      integer iquadstart
      iquadstart = 1

      do i = 1,4
        if(node_quad(i,inode).eq.0) then
        iquadstart = i
        node_quad(iquadstart,inode) = iquad
        return
        endif
      enddo

      return 
      end
C-----------------------------------------------------------------
      subroutine addTo_node_hex(inode,ihex)
#     include "SIZE"
      integer inode,ihex
      integer ihexstart
      ihexstart = 1

      do i = 1,8
        if(node_hex(i,inode).eq.0) then
        ihexstart = i
        node_hex(ihexstart,inode) = ihex
        return
        endif
      enddo

      return 
      end
C-----------------------------------------------------------------
      subroutine convert
c
c  Subroutine to convert gmsh hex20/hex27 to nek hex20/hex27 elements.

#     include "SIZE"

      integer msh_to_nek(27) ! gmsh hex27 node order map to nek hex27 node order
      data    msh_to_nek
     &      /1,3,9,7,19,21,27,25,2,4,10,6,12,
     &       8,18,16,20,22,24,26,5,11,13,15,17,23,14/
 
      integer hex_face_node(4,6)
      data hex_face_node
     &      /1,2,6,5,2,3,7,6,3,4,8,7,1,4,8,5,
     &       1,2,3,4,5,6,7,8/
 
      integer fnode(4),quadnode(4),fhex(32),fhex_nd(32)
      integer ihex,imshvert,inekvert
      integer ifoundquad,imatch
      integer physicalTag
      integer iaddhex,addhex
      logical ifnew

	  
      write(6,'(A)') ' '
      write(6,'(A)') 'Converting elements ... '
      do ihex = 1, totalHex
        do imshvert = 1,27 ! 20
          inekvert = msh_to_nek(imshvert)
          xm1(inekvert,1,1,ihex)= node_xyz(1,hex_array(imshvert+2,ihex))
          ym1(inekvert,1,1,ihex)= node_xyz(2,hex_array(imshvert+2,ihex))
          zm1(inekvert,1,1,ihex)= node_xyz(3,hex_array(imshvert+2,ihex))
        enddo
      enddo
      write(6,'(A)') 'Done :: Converting elements '  

      write(6,*) 'Converting boundary conditions...' 
c zero-out bc and curve sides arrays
      call blank   (cbc,3*2*ldim*max_num_elem)
      call rzero   (bc,5*2*ldim*max_num_elem)
      call blank   (ccurve,(4+8*(ldim-2))*max_num_elem)
      call rzero   (curve,2*ldim*12*max_num_elem)

c currently, does consider converting boundary condition now.
c if need to assigne boundary condition,
c for each hex face, should search which quad is associate with it.
c set bc's
c only need to associate quad4 to hex8 faces.

      do ihex = 1, totalHex
        do iface = 1,6
          hex_face_array(iface,ihex) = 0
        enddo
      enddo

C---- new scheme to search for boundaries.
c---- because node_quad, and node_hex now contains information of this node.
c node_quad((1-4),inode) contains the quad ids associate with this inode
c node_hex((1-8),inode) contains the hex ids associate with this inode
c this new scheme O(N) should be much faster than old scheme O(N^2)

c search boundary using loop in all quads.
      do iquad = 1,totalQuad
         physicalTag = quad_array(1,iquad)
	                !quad_array(2,iquad) ! or  geometrical tag
	     do inode = 1,4
          quadnode(inode) = quad_array(inode+2,iquad) ! first 4 nodes of quad
         enddo

C find all hexes that share the same nodes of this quad
         call rzero(fhex,32)
         call rzero(fhex_nd,32)
         iaddhex = 0
         do inode = 1,4
            do ihex = 1,8
               ! find all hex id related to the nodes in this quad.
			   ! there will be duplicated hex id 
               if(node_hex(ihex,quadnode(inode)).gt.0) then
               iaddhex = iaddhex + 1
               fhex(iaddhex) = node_hex(ihex,quadnode(inode))
               endif
            enddo
         enddo

c fhex(32) contains all hexes that share the same nodes of this quad.
c now, only need to loop over this fhex(32) to find which hex face correspond to this quad.
c however, there are duplicated hex id in fhex(32)
C eliminate duplicated hex id in fhex, and store to fhex_nd

       addhex = iaddhex
       iaddhex_nd = 1
       fhex_nd(1) = fhex(1)
       do iaddhex1 = 2,addhex  ! loop in fhex 
            ifnew = .TRUE.
            do iaddhex2 =1,iaddhex_nd ! loop in fhex_nd
			! if duplicate hex id, ifnew = false
            if(fhex(iaddhex1).eq.fhex_nd(iaddhex2)) ifnew=.FALSE.		
            enddo
            if(ifnew) then
              iaddhex_nd = iaddhex_nd + 1
              fhex_nd(iaddhex_nd) = fhex(iaddhex1)
            endif
       enddo

c look over in fhex_nd to find 
c which hex which face is corresponding to this quad

       do iaddhex = 1,iaddhex_nd
         ihex = fhex_nd(iaddhex) ! hex id that are related to the nodes of quad
         do iface = 1,6       ! loop all faces
         do ifnode = 1,4
         fnode(ifnode)=hex_array(hex_face_node(ifnode,iface)+2,ihex)
         enddo
			  ! fnode is the node ids of hex face
			  ! quadnode is the node ids of quad
         call ifquadmatch(imatch,fnode,quadnode)
         ! write(6,*) imatch 
         if(imatch.eq.1) then ! if match. jump out.
c         write(6,*) 'found hex face for quad id ',iquad 
         goto 1040
         endif
         enddo
       enddo

c if no hex face is found for this quad, report error
      write(6,*) 'ERROR: cannot find hex face for quad id ',iquad
      write(6,*) 'ERROR: this should not happen, please check your mesh'
      write(6,*) 'ERROR: or your mesh exporting process in gmsh'

c assign physical tag to hex face
1040  hex_face_array(iface,ihex) = physicalTag
      enddo

C---- old scheme to search for boundary setup
C---- this is very slow O(N^2) scheme. 	  
c      do ihex = 1, totalHex
c        do iface = 1,6
c            ! obtain node id for this face.
c            do ifnode = 1,4
c            fnode(ifnode)=hex_array(hex_face_node(ifnode,iface)+2,ihex)
c            enddo
c            !input fnode, return quad elements number iquad
c            call findquad(fnode,ifoundquad)   ! ifoundquad is the matching quad number
c                                              ! ifoundquad is 0 if no quad element find
c            if(ifoundquad.ne.0) then
c            hex_face_array(iface,ihex) = quad_array(1,ifoundquad) ! physical tag
c                                         !quad_array(2,ifoundquad) ! geometrical tag
c            endif
c        enddo
c      enddo

cc assign dummy boundary condition and id to bc array	  
      do ihex = 1, totalHex
        do iface = 1,6
         if(hex_face_array(iface,ihex).ne.0) then       ! if on boundary, with physical tag
          cbc(iface,ihex) = 'MSH'                       ! dummy boundary condition
          bc(5,iface,ihex) = hex_face_array(iface,ihex) ! assign tag 
         endif
        enddo
      enddo

      write(6,*) 'Done::Converting boundary conditions' 
 
      return
      end
C--------------------------------------------------------------------
      subroutine findquad(fnode,ifoundquad)
#     include "SIZE"
      integer fnode(4),ifoundquad,iquad
      integer quadnode(4),imatch

      ifoundquad = 0
      iquad = 0
      imatch = 0
c loop over all quad to find the quad has the fnode numbers.
c
      do iquad = 1,totalQuad

         do inode = 1,4
          quadnode(inode) = quad_array(inode+2,iquad)
         enddo

         call ifquadmatch(imatch,fnode,quadnode)

         if(imatch.eq.1) then
          ifoundquad = iquad
          return
         endif
      enddo

      return
      end
C--------------------------------------------------------------------
      subroutine ifquadmatch(imatch,fnode,quadnode)
      integer fnode(4),quadnode(4)
      integer imatch

      imatch = 0
  
      do ifnode = 1,4
         do iquadnode = 1,4
         if(fnode(ifnode).eq.quadnode(iquadnode)) imatch=imatch+1
         enddo
      enddo
 
c imatch should equal 4 if fnode and quadnode is matching.
      if(imatch.eq.4) then
       imatch = 1
      else
       imatch =0
      endif
 
      return
      end
C-----------------------------------------------------------------
      subroutine setbc
#     include "SIZE"
c  set boundary condition 
c  read from a file casename.bc

      parameter (npbc_max=10) ! maximum pairs of periodic boundary condition

      integer hex_face_node(4,6)
      data hex_face_node
     &      /1,2,6,5,2,3,7,6,3,4,8,7,1,4,8,5,
     &       1,2,3,4,5,6,7,8/
	 
      integer parray(2,2,max_num_elem)

      character*3 ubc
      integer tags(2),ibc,nbc,io
      integer ip,np,ipe,ipe2,nipe(2)
      integer ptags(2,npbc_max),fnode(4)
      real pvec(3,npbc_max),fpxyz(3,2)
      real AB_v(3),AD_v(3),farea,product_v(3)
      real dist,distMax,ptol
	  
      character*32 bcname
      character*1 bcnam1(32)
      equivalence(bcname,bcnam1)
	  
      call blank (bcname,32)

      len = ltrunc(mshname,32)
      call chcopy(bcnam1,mshname,(len-3))
      call chcopy(bcnam1(len-3),'.bc' , 3)
      len = ltrunc(bcnam1,32)
	  
      open(301,file=bcname,err=1010)
      write(6,*) 'Setting boundary condition from ',bcname(:len),' file' 

      read(301,*,iostat=io) nbc
      if(io.ne.0) then
         write(6,*) bcname(:len),' file is empty ',
     &'please set boundary condition in .usr file'
         return
      endif		   

        do ibc = 1,nbc
        call blank (ubc,3)
        read(301,*) tags(1),ubc
c not periodic boundary condition, direct set it up
          write(6,*) 'setting ',ubc,' to surface ',tags(1)
          do ihex = 1, totalHex
            do iface = 1,6
               if(bc(5,iface,ihex).eq.tags(1)) then
                 cbc(iface,ihex) = ubc
               endif
            enddo
          enddo
      enddo
      ip = 0
      read(301,*,iostat=io) nbc,ptol ! nbc is the number of pairs of periodic boundary condition
                           ! ptol is the tolerence used to search periodic boundary conditin
      if(io.ne.0) then
         write(6,*) 'No periodic boundary condition set.'
         return
      endif					
						  
      if(nbc.gt.npbc_max) then
         write(6,*) 'ERROR: increase npbc_max to ',nbc
         return
      endif
	  
      do ibc = 1,nbc
	    ip = ip + 1
        read(301,*) tags(1),tags(2),pvec(1,ip),pvec(2,ip),pvec(3,ip)
c set periodic boundary condition.
c read periodic mapping vector
        ptags(1,ip) = tags(1)
        ptags(2,ip) = tags(2)
      enddo
      close(301)
 
      write(6,*) 'Setting periodic boundary condition' 
      np = ip
      do ip = 1,np
          ipe = 0
          do ihex = 1, totalHex
            do iface = 1,6
               if(bc(5,iface,ihex).eq.ptags(1,ip)) then
                ipe = ipe + 1
                parray(1,1,ipe) = ihex
                parray(2,1,ipe) = iface
               endif
            enddo
          enddo
          nipe(1) = ipe
	  
          ipe = 0
	      do ihex = 1, totalHex
            do iface = 1,6
               if(bc(5,iface,ihex).eq.ptags(2,ip)) then
                ipe = ipe + 1
                parray(1,2,ipe) = ihex
                parray(2,2,ipe) = iface
               endif
            enddo
          enddo
          nipe(2) = ipe
		  
          write(6,*)'maping surface',ptags(1,ip),'with',nipe(1),
     & 'faces'
          write(6,*)'to surface',ptags(2,ip),'with',nipe(2),
     & 'faces'

          if(nipe(1).ne.nipe(2))  then
            write(6,*) 'EORROR, face numbers are not matching'
          endif

          do ipe = 1,nipe(1)
             ihex = parray(1,1,ipe)
             iface = parray(2,1,ipe)
c get face center xyz
             call rzero(fpxyz(1,1),3)
             do ifnode = 1,4
             fnode(ifnode)=hex_array(hex_face_node(ifnode,iface)+2,ihex)
             fpxyz(1,1) = fpxyz(1,1) + node_xyz(1,fnode(ifnode))/4.0
             fpxyz(2,1) = fpxyz(2,1) + node_xyz(2,fnode(ifnode))/4.0
             fpxyz(3,1) = fpxyz(3,1) + node_xyz(3,fnode(ifnode))/4.0
             enddo

c lopp over mapped faces to find its mapping face
             distMax = 1000.0
             do ipe2 = 1,nipe(2)
                ihex2 = parray(1,2,ipe2)
                iface2 = parray(2,2,ipe2)
c get face center xyz
               call rzero(fpxyz(1,2),3)
               do ifnode = 1,4
           fnode(ifnode)=hex_array(hex_face_node(ifnode,iface2)+2,ihex2)
               fpxyz(1,2) = fpxyz(1,2) + node_xyz(1,fnode(ifnode))/4.0
               fpxyz(2,2) = fpxyz(2,2) + node_xyz(2,fnode(ifnode))/4.0
               fpxyz(3,2) = fpxyz(3,2) + node_xyz(3,fnode(ifnode))/4.0
               enddo
               
              dist = sqrt((fpxyz(1,2) - fpxyz(1,1) - pvec(1,ip))**2
     & + (fpxyz(2,2) - fpxyz(2,1) - pvec(2,ip))**2
     & + (fpxyz(3,2) - fpxyz(3,1) - pvec(3,ip))**2)
 
               if(dist.lt.distMax) then 
                  distMax = dist
                  !write(6,*) distMax
                  if(distMax.le.ptol) then
                  bc(1,iface,ihex) = ihex2*1.0
                  bc(2,iface,ihex) = iface2*1.0
                  bc(1,iface2,ihex2) = ihex*1.0
                  bc(2,iface2,ihex2) = iface*1.0
                  cbc(iface,ihex) = 'P  '
                  cbc(iface2,ihex2) = 'P  '
c               write(6,*) ihex,iface,bc(1,iface,ihex),bc(2,iface,ihex)
c          write(6,*) ihex2,iface2,bc(1,iface2,ihex2),bc(2,iface2,ihex2)
c          write(6,*) fpxyz(1,1),fpxyz(2,1),fpxyz(3,1)
c          write(6,*) fpxyz(1,2),fpxyz(2,2),fpxyz(3,2)
                  endif
               endif
             enddo
          enddo

          nperror = 0
		  
          write(6,*)'doing periodic check for surface',ptags(1,ip)

          do ipe = 1,nipe(1)
             ihex = parray(1,1,ipe)
             iface = parray(2,1,ipe)
             if (cbc(iface,ihex).ne.'P  ') then
                  nperror = nperror +1 
             endif
          enddo
          if (nperror.gt.0) write(6,*) 'ERROR,',nperror,
     & 'faces did not map'

          nperror = 0

          do ipe = 1,nipe(1)
             ihex = parray(1,1,ipe)
             iface = parray(2,1,ipe)
             ihex2 = int(bc(1,iface,ihex))
             iface2 = int(bc(2,iface,ihex))
             ihex3 = int(bc(1,iface2,ihex2))
             iface3 = int(bc(2,iface2,ihex2))
             if ((ihex.ne.ihex3).or.(ihex.ne.ihex3)) then
			 
cc for debug use only
cc
c                write(6,*) 'ERROR,',ihex,iface,' map to ',ihex2,iface2
c                write(6,*) 'but,',ihex2,iface2,' map to ',ihex3,iface3
c
c		     do ifnode = 1,4
c         fnode(ifnode)=hex_face_node(ifnode,iface)
c      write(6,*)xm1(fnode(ifnode),1,1,ihex),ym1(fnode(ifnode),1,1,ihex),
c     & zm1(fnode(ifnode),1,1,ihex) 
c             enddo
c			 
c		     do ifnode = 1,4
c         fnode(ifnode)=hex_face_node(ifnode,iface2)
c      write(6,*)xm1(fnode(ifnode),1,1,ihex2),
c     & ym1(fnode(ifnode),1,1,ihex2),zm1(fnode(ifnode),1,1,ihex2) 
c             enddo
c
c		     do ifnode = 1,4
c         fnode(ifnode)=hex_face_node(ifnode,iface3)
c      write(6,*)xm1(fnode(ifnode),1,1,ihex3),
c     & ym1(fnode(ifnode),1,1,ihex3),zm1(fnode(ifnode),1,1,ihex3) 
c             enddo

                nperror = nperror + 1
             endif
          enddo		  

          if (nperror.gt.0) then
          write(6,*) 'ERROR,',nperror,'faces are wrong',
     & 'out of total ',nipe(1),' faces'
          endif
		  
      enddo

1010  return
      end
c--------------------------------------------------------------------
c-----------------------------------------------------------------------
      subroutine cross_product(AB_v,AC_v,prod_v,area)
cc calculate cross product of two vectors

      real*8 AB_v(3),AC_v(3),prod_v(3)
      real*8 area

      prod_v(1) = AB_v(2)*AC_v(3)-AB_v(3)*AC_v(2)
      prod_v(2) = AB_v(3)*AC_v(1)-AB_v(1)*AC_v(3)
      prod_v(3) = AB_v(1)*AC_v(2)-AB_v(2)*AC_v(1)
	  
      area = prod_v(1)**2 + prod_v(2)**2+prod_v(3)**2
      area = sqrt(area)
c      area = area/2.0
	  
      return 
      end
c------------------------------------------------------------------------------------------
      subroutine gen_re2
#     include "SIZE"

      write(6,*)
      write(6,'(A,A)') 'writing ', re2name

      call open_re2
      call write_xyz
      call write_curve
      call write_bc
      call close_re2

      return
      end
C--------------------------------------------------------------------
      subroutine open_re2

#     include "SIZE"

      character*80  hdr


      real*4 test
      data   test  / 6.54321 /

      call byte_open(re2name,ierr)
            
c  Write the header
      call blank     (hdr,80)    
      write(hdr,1) num_elem, num_dim, num_elem
    1 format('#v003',i9,i3,i9,' this is the hdr')
      call byte_write(hdr,20,ierr)         
      call byte_write(test,1,ierr)     ! write the endian discriminator

      return
      end
C--------------------------------------------------------------------
      subroutine write_xyz

#     include "SIZE"

      real     xx(8), yy(8), zz(8)
      real*8   rgroup, buf2(30)

      integer isym2pre(8)   ! Symmetric-to-prenek vertex ordering
      data    isym2pre    / 1 , 2 , 4 , 3 , 5 , 6 , 8 , 7 /

      nxs = 3-1  ! always nx1=ny1=nz1=3 here
      nys = 3-1
      nzs = 3-1

      igr    = 0
      rgroup = igr

      do iel=1,num_elem
        l = 0
        if (num_dim.eq.3) then
          do k=0,1
          do j=0,1
          do i=0,1
            l=l+1
            li=isym2pre(l)
            xx(li) = xm1(1+i*nxs,1+j*nys,1+k*nzs,iel)
            yy(li) = ym1(1+i*nxs,1+j*nys,1+k*nzs,iel)
            zz(li) = zm1(1+i*nxs,1+j*nys,1+k*nzs,iel)
          enddo
          enddo
          enddo
        else
          do j=0,1
          do i=0,1
            l=l+1
            li=isym2pre(l)
            xx(li) = xm1(1+i*nxs,1+j*nys,1,iel)
            yy(li) = ym1(1+i*nxs,1+j*nys,1,iel)
          enddo
          enddo
        endif

        call byte_write(rgroup, 2,ierr)

        if (num_dim.eq.3) then
          call copy        (buf2(1) ,xx,8)
          call copy        (buf2(9) ,yy,8)
          call copy        (buf2(17),zz,8)
          call byte_write  (buf2(1) ,16, ierr)
          call byte_write  (buf2(9) ,16, ierr)
          call byte_write  (buf2(17),16, ierr)
        else
          call copy        (buf2(1),xx,4)
          call copy        (buf2(5),yy,4)
          call byte_write  (buf2(1),8, ierr)
          call byte_write  (buf2(5),8, ierr)
        endif
      enddo

      return
      end
C-----------------------------------------------------------------------
      subroutine write_curve

#     include "SIZE"

      real*8     buf2(30)
      real*8     rcurve

      character*1 cc

      do iel=1,num_elem
         call gen_rea_midside_e(iel)
      enddo

      nedge = 4 + 8*(num_dim-2)
      ncurv = 0
      do iel=1,num_elem
        do iedge=1,nedge
          if (ccurve(iedge,iel).ne.' ') ncurv = ncurv+1
        enddo
      enddo
 
      rcurve = ncurv
      call byte_write(rcurve,2, ierr)

      do iel=1,num_elem
        do iedge=1,nedge
          if (ccurve(iedge,iel).ne.' ') then
            if (ccurve(iedge,iel).eq.'C') cc='C'
            if (ccurve(iedge,iel).eq.'s') cc='s'
            if (ccurve(iedge,iel).eq.'m') cc='m'
            buf2(1) = iel
            buf2(2) = iedge
            call copy       (buf2(3),curve(1,iedge,iel),5)
            call blank      (buf2(8),8)
            call chcopy     (buf2(8),cc,1)
            call byte_write (buf2,16, ierr)
          endif
        enddo
      enddo

      return
      end
C-----------------------------------------------------------------------
      subroutine write_bc
      
#     include "SIZE"

      real*8  rbc, buf2(30)

      character*3 ch3
      character*1 chdum
      data        chdum /' '/

      nface = 2*num_dim
      nbc   = 0

      do iel=1,num_elem
        do ifc=1,nface
          if (cbc(ifc,iel).ne.'   ')  nbc = nbc + 1
        enddo
      enddo

      rbc = nbc
      call byte_write (rbc,2, ierr)

      do iel = 1,num_elem
        do ifc = 1,nface
          ch3 = cbc(ifc,iel)
c          if (ch3.eq.'MSH') then
          if (ch3.ne.'   ') then ! setting boundary condition.
            buf2(1)=iel
            buf2(2)=ifc
            call copy   (buf2(3),bc(1,ifc,iel),5)
            call blank  (buf2(8),8)
            call chcopy (buf2(8),ch3,3)
            if (num_elem.ge.1 000 000) then
              ibc     = bc(1,ifc,iel)
              buf2(3) = ibc
            endif
            call byte_write (buf2,16, ierr)
          endif
        enddo
      enddo

      return
      end 
C-----------------------------------------------------------------------
      subroutine close_re2

      call byte_close (ierr)

      return
      end
c-----------------------------------------------------------------------
      subroutine gen_rea_midside_e(e)

#     include "SIZE"

      real        len
      real        x3(27),y3(27),z3(27),xyz(3,3)
      character*1 ccrve(12)
      integer     e,edge

      integer e3(3,12)
      save    e3
      data    e3 /  1, 2, 3,    3, 6, 9,    9, 8, 7,    7, 4, 1
     $           , 19,20,21,   21,24,27,   27,26,25,   25,22,19
     $           ,  1,10,19,    3,12,21,    9,18,27,    7,16,25 /

      call chcopy(ccrve,ccurve(1,e),12)

      call map2reg     (x3,3,xm1(1,1,1,e),1)  ! Map to 3x3x3 array
      call map2reg     (y3,3,ym1(1,1,1,e),1)
      if (num_dim.eq.3) call map2reg    (z3,3,zm1(1,1,1,e),1)

c     Take care of spherical curved face defn
      if (ccurve(5,e).eq.'s') then
         call chcopy (ccrve(1),'ssss',4) ! face 5
         call chcopy (ccrve(5),' ',1)    ! face 5
      endif
      if (ccurve(6,e).eq.'s') then
         call chcopy (ccrve(5),'ssss',4) ! face 6
      endif

      tol   = 1.e-4
      tol2  = tol**2
      nedge = 4 + 8*(num_dim-2)

      do i=1,nedge
         if (ccrve(i).eq.' ') then
            do j=1,3
               xyz(1,j) = x3(e3(j,i))
               xyz(2,j) = y3(e3(j,i))
               xyz(3,j) = z3(e3(j,i))
            enddo
            len = 0.
            h   = 0.
            do j=1,num_dim
               xmid = .5*(xyz(j,1)+xyz(j,3))
               h    = h   + (xyz(j,2)-xmid)**2
               len  = len + (xyz(j,3)-xyz(j,1))**2
            enddo
            if (h.gt.tol2*len) ccurve(i,e) = 'm'
            if (h.gt.tol2*len) call copy (curve(1,i,e),xyz(1,2),num_dim)
         endif
      enddo

      return
      end
c-----------------------------------------------------------------------
      subroutine map2reg(ur,n,u,nel)
c
c     Map scalar field u() to regular n x n x n array ur

#     include "SIZE"

      real    ur(1), u(3*3*3,1)
      integer e

      ldr = n**num_dim

      k=1
      do e=1,nel
         if (num_dim.eq.2) call map2reg_2di_e(ur(k),n,u(1,e),3)
         if (num_dim.eq.3) call map2reg_3di_e(ur(k),n,u(1,e),3)
         k = k + ldr
      enddo

      return
      end
c-----------------------------------------------------------------------
      subroutine map2reg_2di_e(uf,n,uc,m) ! Fine, uniform pt

      real uf(n,n),uc(m,m)

      parameter (l=50)
      common /cmap2d/ j(l*l),jt(l*l),w(l*l),z(l)

      integer mo,no
      save    mo,no
      data    mo,no / 0,0 /

      if (m.ne.mo .or. n.ne.no ) then

          call zwgll (z,w,m)
          call zuni  (w,n)

          call gen_int_gz(j,jt,w,n,z,m)

      endif

      call mxmf2(j,n,uc,m,w ,m)
      call mxmf2(w,n,jt,m,uf,n)

      return
      end
c-----------------------------------------------------------------------
      subroutine map2reg_3di_e(uf,n,uc,m) ! Fine, uniform pt

      real uf(n,n,n),uc(m,m,m)

      parameter (l=50)
      common /cmap3d/ j(l*l),jt(l*l),v(l*l*l),w(l*l*l),z(l)

      integer mo,no
      save    mo,no
      data    mo,no / 0,0 /

      if (m.ne.mo .or. n.ne.no ) then

          call zwgll (z,w,m)
          call zuni  (w,n)

          call gen_int_gz(j,jt,w,n,z,m)

      endif

      mm = m*m
      mn = m*n
      nn = n*n

      call mxmf2(j,n,uc,m,v ,mm)
      iv=1
      iw=1
      do k=1,m
         call mxmf2(v(iv),n,jt,m,w(iw),n)
         iv = iv+mn
         iw = iw+nn
      enddo
      call mxmf2(w,nn,jt,m,uf,n)

      return
      end
c-----------------------------------------------------------------------
      subroutine gen_int_gz(j,jt,g,n,z,m)

c     Generate interpolater from m z points to n g points

c        j   = interpolation matrix, mapping from z to g
c        jt  = transpose of interpolation matrix
c        m   = number of points on z grid
c        n   = number of points on g grid

      real j(n,m),jt(m,n),g(n),z(m)

      mpoly  = m-1
      do i=1,n
         call fd_weights_full(g(i),z,mpoly,0,jt(1,i))
      enddo

      call transpose(j,n,jt,m)

      return
      end
c-----------------------------------------------------------------------
      SUBROUTINE BLANK(A,N)
      CHARACTER*1 A(1)
      CHARACTER*1 BLNK
      SAVE        BLNK
      DATA        BLNK /' '/
C
      DO 10 I=1,N
         A(I)=BLNK
   10 CONTINUE
      RETURN
      END
c-----------------------------------------------------------------------
      subroutine exitt

      stop
      return
      end
c-----------------------------------------------------------------------
      subroutine copy(a,b,n)
      real a(1),b(1)

      do i=1,n
         a(i)=b(i)
      enddo

      return
      end
c-----------------------------------------------------------------------
      subroutine chcopy(a,b,n)
      CHARACTER*1 A(1), B(1)
C
      DO 100 I = 1, N
 100     A(I) = B(I)
      return
      END
c-----------------------------------------------------------------------
      subroutine icopy(a,b,n)
      INTEGER A(1), B(1)
C
      DO 100 I = 1, N
 100     A(I) = B(I)
      return
      END
c-----------------------------------------------------------------------
      subroutine fd_weights_full(xx,x,n,m,c)
c
c     This routine evaluates the derivative based on all points
c     in the stencils.  It is more memory efficient than "fd_weights"
c
c     This set of routines comes from the appendix of
c     A Practical Guide to Pseudospectral Methods, B. Fornberg
c     Cambridge Univ. Press, 1996.   (pff)
c
c     Input parameters:
c       xx -- point at wich the approximations are to be accurate
c       x  -- array of x-ordinates:   x(0:n)
c       n  -- polynomial degree of interpolant (# of points := n+1)
c       m  -- highest order of derivative to be approxxmated at xi
c
c     Output:
c       c  -- set of coefficients c(0:n,0:m).
c             c(j,k) is to be applied at x(j) when
c             the kth derivative is approxxmated by a
c             stencil extending over x(0),x(1),...x(n).
c
c
      real x(0:n),c(0:n,0:m)
c
      c1       = 1.
      c4       = x(0) - xx
c
      do k=0,m
      do j=0,n
         c(j,k) = 0.
      enddo
      enddo
      c(0,0) = 1.
c
      do i=1,n
         mn = min(i,m)
         c2 = 1.
         c5 = c4
         c4 = x(i)-xx
         do j=0,i-1
            c3 = x(i)-x(j)
            c2 = c2*c3
            do k=mn,1,-1
               c(i,k) = c1*(k*c(i-1,k-1)-c5*c(i-1,k))/c2
            enddo
            c(i,0) = -c1*c5*c(i-1,0)/c2
            do k=mn,1,-1
               c(j,k) = (c4*c(j,k)-k*c(j,k-1))/c3
            enddo
            c(j,0) = c4*c(j,0)/c3
         enddo
         c1 = c2
      enddo

      return
      end
c-----------------------------------------------------------------------
      subroutine transpose(a,lda,b,ldb)
      real a(lda,1),b(ldb,1)
c
      do j=1,ldb
         do i=1,lda
            a(i,j) = b(j,i)
         enddo
      enddo

      return
      end
c-----------------------------------------------------------------------
      function ltrunc(string,l)
      CHARACTER*1 STRING(L)
      CHARACTER*1   BLNK
      DATA BLNK/' '/
C
      DO 100 I=L,1,-1
         L1=I
         IF (STRING(I).NE.BLNK) GOTO 200
  100 CONTINUE
      L1=0
  200 CONTINUE
      LTRUNC=L1
      return
      END
c-----------------------------------------------------------------------
      subroutine zuni(z,np)
c
c     Generate equaly spaced np points on the interval [-1:1]
c
      real z(1)

      dz = 2./(np-1)
      z(1) = -1.
      do i = 2,np-1
         z(i) = z(i-1) + dz
      enddo
      z(np) = 1.

      return
      end
c-----------------------------------------------------------------------
      subroutine rzero(a,n)
      real A(1)
      DO 100 I = 1, N
 100     A(I ) = 0.0
      return
      END

