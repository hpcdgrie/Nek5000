c----------------------------------------------------------------------
c Noted by Haomin Yuan, 8/17/2018
c
c the original exo2nek convert pure hex20 exodus file to re2 file.
c side set information is stored passed to re2 file.
c
c now we want to add some feature to exo2nek.
c 1. convert 1 tet element to 4 hex element
c 2. convert 1 wedge element to 3 hex element, for boundary layer.
c 3. maintain all connectivity and  sideset information.
c
c use the exodus_read to read exodus file.
c use a new subroutine, split, to split 1 tet to 4 hex, 1 wedge to 3 hex.
c Haomin Yuan, 8/21/2018
c----------------------------------------------------------------------
c Noted by  Haomin Yuan, 8/29/2018
c in the fianl version of the code, there are 3 options. please see code
c output.
c 
c in option 3, sideset infomation is successfully passed to nek mesh.
c option 2 has not been tested.
c 
c now, we slightly improve the sideset setting in option 1, which is 
c very unefficient right now. 
c
c----------------------------------------------------------------------
      program exo2nek

#     include "SIZE"

      integer option

      write(6,*) 'Please input option:'
      write(6,*) '1: assume pure hex20 mesh, original exo2nek'
      write(6,*) '2: assume hybrid tetra4 and wedge6 mesh'
      write(6,*) '    all tetra4 elements in block 1 (or first block)'
      write(6,*) '    all wedge6 elements in block 2 (or second block)'
      write(6,*) '3: assume hybrid hex20, tetra10 and wedge15 mesh'
      write(6,*) '    all tetra4 elements in block 1 (or first block)'
      write(6,*) '    all hex8   elements in block 2 (or second block)'
      write(6,*) '    all wedge6 elements in block 3 (or third block)'

      read (5,'(I1)') option
	  
      if(option.EQ.1) then
cc assumes all element to be HEX20.
cc this part is 
      call read_input_name
      call exodus_read
      call convert

      elseif(option.EQ.2) then
cc assume all element to be TRTRA10 and WEDGE15
cc TRTRA10 elements in block 1
cc WEDGE15 elements in block 2, one WEDGE15 exo element is converted to 3 hex20 nek elements.
      call read_input_name
      call exodus_read_new ! exodus_read_new does not check element type in exodus file.
      call split_convert1

      elseif(option.EQ.3) then
cc assume all element to be HEX20 TRTRA10 and WEDGE15
cc TRTRA10 elements in block 1
cc HEX20 elements in block 2, one HEX20 exo element is converted into 8 hex20 nek elements
cc WEDGE15 elements in block 3, one WEDGE15 exo element is converted to 6 hex20 nek elements.
c
c this feature is added to handle the urban case, in which far field mesh is extruded.
c some other cases may also be benefited.
      call read_input_name
      call exodus_read_new ! exodus_read_new does not check element type in exodus file.
      call split_convert2

      else
       write(6,*) 'Unknown option for exo2nek, ABORT.'
      endif
      
      call setbc
      call gen_re2
	  
      end 
c-----------------------------------------------------------------------
      subroutine read_input_name

#     include "SIZE"

      character*1 re2nam1(80)
      character*1 exonam1(32)

      equivalence(re2name,re2nam1)
      equivalence(exoname,exonam1) 

      call blank (exoname, 32)
      call blank (re2name, 80)

      write(6,*) 'Input (.exo) file name:'
      read (5,'(a32)') exoname
      len = ltrunc(exoname,32)

      call chcopy(re2name        ,exoname,32)
      call chcopy(exonam1(len+1) ,'.exo' , 4)
      call chcopy(re2nam1(len+1) ,'.re2' , 4)

      return 
      end
c-----------------------------------------------------------------------
      subroutine exodus_read_new
cc assume tet10 and wedge15 elements.
cc return tet10 element number and wedge 15 element number
cc or a flag showing wether this element is a tet 10 or wedge 15 elements
c
c  Subroutine to read an exodusII binary file containing a mesh.
c  It uses exodus fortran binding subroutines, which depend on
c  the netcdf library for low level data access.
c
      include 'exodusII.inc'
#     include "SIZE"

      integer exoid, cpu_ws, io_ws

      character*(MXSTLN) typ, qa_record(4,10)
      character*(MXLNLN) titl
      character*1        cdum

      integer idblk              (max_num_elem_blk)
      integer num_attr           (max_num_elem_blk)  ! not used
      integer num_nodes_per_elem (max_num_elem_blk)
c
c open EXODUS II file
c
      cpu_ws = 8 ! use real*8 to communicate with exodus
      io_ws  = 0
      exoid  = exopen (exoname, EXREAD, cpu_ws, io_ws, vers, ierr)
      if (ierr.lt.0) then
        write(6,'(2a)') "ERROR: cannot open file ", exoname 
        STOP
      endif
      write(6,*)
      write(6,'(a32,a,f4.2)') 
     &      exoname," is an EXODUSII file; version ",vers
      write(6,'(a,i2)') "I/O word size", io_ws
c
c read database parameters
c
      call exgini (exoid, titl, num_dim, num_nodes, num_elem,
     &             num_elem_blk, num_node_sets, num_side_sets, ierr)
      if (ierr.lt.0) then
        write(6,'(a)') "ERROR: cannot read exodusII parameters (exgini)"
        STOP
      endif
      write (6,    '(/"database parameters:"/ /
     &               "title         = ", a81 / /
     &               "num_dim       = ", i8 /
     &               "num_nodes     = ", i8 /
     &               "num_elem      = ", i8 /
     &               "num_elem_blk  = ", i8 /
     &               "num_side_sets = ", i8)')
     &               titl,num_dim, num_nodes, num_elem,
     &               num_elem_blk, num_side_sets
      write (6,*)
c
c perform some checks
c
      if (num_elem.gt.max_num_elem) then
        write(6,*) 'Abort: number of elements too large',num_elem
        write(6,*) 'change MAXNEL and recompile'
        STOP
      endif
      if (num_side_sets.gt.max_num_sidesets) then
        write(6,'(a)')
     &    "ERROR: number of sidesets > max_num_sidesets. "
        write(6,'(a,i2,a)') "Set max_num_sidesets >= ",num_side_sets,
     &                      " and recompile exo2nek. "
        STOP
      endif
c
c read element block parameters
c
      call exgebi (exoid, idblk, ierr)
      if (ierr.lt.0) then
        write(6,'(a)') "ERROR: cannot read block ids (exgebi)"
        STOP
      endif

      do i = 1, num_elem_blk
        call exgelb (exoid, idblk(i), typ, num_elem_in_block(i),
     &               num_nodes_per_elem(i), num_attr(i), ierr)
        if (ierr.lt.0) then
          write(6,'(a,i3,a)')
     &    "ERROR: cannot read parameters for block ",i," (exgelb)"
          STOP
        endif
        write (6, '("element block id   = ", i8,/
     &              "element type       = ", a8,/
     &              "num_elem_in_block  = ", i8,/
     &              "num_nodes_per_elem = ", i8)')
     &              idblk(i), typ, num_elem_in_block(i),
     &              num_nodes_per_elem(i)

      enddo
c
c read nodal coordinates values from database
c
      call exgcor (exoid, x_exo, y_exo, z_exo, ierr)
      if (ierr.lt.0) then
        write(6,'(a)') "ERROR: cannot read nodal coordinates (exgcor)"
        STOP
      endif
c
c read element connectivity
c
      iend = 0
      do 60 i = 1, num_elem_blk
        istart = iend + 1
        call exgelc (exoid, idblk(i), connect(istart), ierr)
        iend = iend + num_nodes_per_elem(i)*num_elem_in_block(i) ! kirk modifed!! 8/23/2018
        if (ierr.lt.0) then
          write(6,'(a)') "ERROR: cannot read elm. connectivity (exgelc)"
          STOP
        endif
60    continue
c
c read individual side sets
c
      if (num_side_sets .gt. 0) then
        call exgssi (exoid, idss, ierr)
        if (ierr.lt.0) then
          write(6,'(a)') "ERROR: cannot read SideSet ids (exgssi)"
          STOP
        endif

        do i = 1, num_side_sets
          call exgsp (exoid,idss(i),num_sides_in_set(i),idum,ierr)
          if (ierr.lt.0) then
            write(6,'(a,i3,a)')
     &      "ERROR: cannot read parameters for SideSet No.",i," (exgsp)"
            STOP
          endif
          write (6, '("side set ", i2, " num_sides = ", i8)')
     &           idss(i), num_sides_in_set(i)
          call exgss (exoid,idss(i),elem_list(1,i),side_list(1,i),ierr)
          if (ierr.lt.0) then
            write(6,'(a,i3,a)')
     &      "ERROR: cannot read parameters for SideSet No.",i," (exgss)"
            STOP
          endif
        enddo
        write (6,*)
      else
        write(6,'(a)') "WARNING: No SideSets in exodus file!"
      endif
c
c read QA records
c
      call exinq (exoid, EXQA, num_qa_rec, fdum, cdum, ierr)
      if (ierr.lt.0) then
        write(6,'(a)') "ERROR: cannot read QA records (exinq QA) "
        STOP
      elseif (ierr.gt.0) then
        write(6,'(a)') "INFO: file does not contain any QA records"
      else
        call exgqa (exoid, qa_record, ierr)
        if (ierr.lt.0) then
          write(6,'(a,i3)') "WARNING: cannot read QA records (exgqa)"
        else
          write (6, '("QA records = ")')
          if (num_qa_rec.gt.10) then
            write(6,'(2a)')
     &        'WARNING: Cannot handle more than 10 QA records',
     &        'Printing only the first 10...'
          else
            do i = 1, num_qa_rec
              do j = 1, 4
                write (6,'(a)') qa_record(j,i)
              enddo
            enddo
          endif
        endif
      endif

      return
      end
C-----------------------------------------------------------------
      subroutine split_convert1
c split 1 TETRA10 to 4 HEX20
C split 1 WEDGE15 to 3 HEX20
C
c  Subroutine to convert an already read exodusII mesh to a nek
c  mesh. The idea is to fill each element's node coordinates of
c  size lx1**3 (3D) or lx1**2 (2D) (lx1=3) with the hex27/quad9
c  coordinates.
c
      include 'exodusII.inc'
#     include "SIZE"

cc in fortran arraysm higher order is at back.
      real*8 hexver(3,27,8) ! hex coordinates
      real*8 tetver(3,10),wedgever(3,15)
      real*8 ehexver(3,20)

cc exoss is used to store sideset information for all exo elements.
cc tet,hex,wedge
      integer exoss(6,max_num_elem) 

      integer tetss(4),wedgess(5),ehexss(6),hexss(6,8)
      integer ehexnumber,tetnumber,wedgenumber
      integer vert_index_exo
      integer iel_nek,iel_exo,ifc_exo


cc Haomin Yuan, I do not know why, rzero cannot used for exoss
      do i=1,6*max_num_elem
      exoss(i,1)=0
      enddo

c store sideset information
      if (num_side_sets.ne.0) then
      write(6,'(a)') ''
      write(6,'(a)') 'Store SideSet information from EXO file'
        do iss=1,num_side_sets   ! loop over ss 
        write(6,'(a,i2,a)') 'Sideset ',idss(iss), ' ...'
           do i=1,num_sides_in_set(iss) ! loop over sides in ss
             iel_exo = elem_list(i,iss) ! exo element number
             ifc_exo = side_list(i,iss) ! exo element face number
             exoss(ifc_exo,iel_exo) = idss(iss)
           enddo
        enddo
      endif

      write(6,'(A)') ' '
      write(6,'(A)') 'Converting elements ... '
	  
      call blank   (cbc,3*2*ldim*max_num_elem)
      call rzero   (bc,5*2*ldim*max_num_elem)
      call blank   (ccurve,(4+8*(ldim-2))*max_num_elem)
      call rzero   (curve,2*ldim*12*max_num_elem)

cc assume block 1 contains all tet elements
cc assume block 2 contants all wedge elements

      tetnumber  = num_elem_in_block(1)

      vert_index_exo = 0
      iel_nek = 0
      do iel_exo = 1, num_elem  
	  
!! number of elements == 4*number of tet + 3*number of wedge

cc if tet.
       if (iel_exo.le.tetnumber) then
cc now convert 1 tet10 element to 4 hex20 elements.

c       do ivert = 1, 10
c       vert_index_exo = vert_index_exo + 1  
c       tetver(1,ivert) = x_exo(connect(vert_index_exo))
c       tetver(2,ivert) = y_exo(connect(vert_index_exo))
c       tetver(3,ivert) = z_exo(connect(vert_index_exo))
c       enddo

cc read tet 4 . 
cc linear interpolate to tet10
       do ivert = 1, 4
       vert_index_exo = vert_index_exo + 1  
       tetver(1,ivert) = x_exo(connect(vert_index_exo))
       tetver(2,ivert) = y_exo(connect(vert_index_exo))
       tetver(3,ivert) = z_exo(connect(vert_index_exo))
       enddo

       call average2vec(tetver(1,5),tetver(1,1),tetver(1,2))
       call average2vec(tetver(1,6),tetver(1,2),tetver(1,3))
       call average2vec(tetver(1,7),tetver(1,1),tetver(1,3))
       call average2vec(tetver(1,8),tetver(1,1),tetver(1,4))
       call average2vec(tetver(1,9),tetver(1,2),tetver(1,4))
       call average2vec(tetver(1,10),tetver(1,3),tetver(1,4))

cc assign sideset to tet elements
       call rzero(tetss,4)
       if (num_side_sets.ne.0) then
        do ifc_exo=1,4
         tetss(ifc_exo) = exoss(ifc_exo,iel_exo)
        enddo
       endif
cc given tet10 vertices, and return you four hex coords.
       call tettohex(hexver,tetver,hexss,tetss)	   
       
       do ihex = 1, 4
          iel_nek = iel_nek + 1
          do inekvert = 1,27
             xm1(inekvert,1,1,iel_nek) = hexver(1,inekvert,ihex)
             ym1(inekvert,1,1,iel_nek) = hexver(2,inekvert,ihex)
             zm1(inekvert,1,1,iel_nek) = hexver(3,inekvert,ihex)
          enddo
          do ifc=1,6
             if(hexss(ifc,ihex).gt.0) then
              bc(5,ifc,iel_nek) = hexss(ifc,ihex)
              cbc(ifc,iel_nek)   = 'EXO' ! dummy boundary condition
             endif
          enddo
       enddo

      else ! if wedge elements 
      
c       do ivert = 1,15
c       vert_index_exo = vert_index_exo + 1  
c       wedgever(1,ivert) = x_exo(connect(vert_index_exo))
c       wedgever(2,ivert) = y_exo(connect(vert_index_exo))
c       wedgever(3,ivert) = z_exo(connect(vert_index_exo))
c       enddo

cc read wedge 6,
cc linear interpolate to wedge 15
       do ivert = 1,6
       vert_index_exo = vert_index_exo + 1  
       wedgever(1,ivert) = x_exo(connect(vert_index_exo))
       wedgever(2,ivert) = y_exo(connect(vert_index_exo))
       wedgever(3,ivert) = z_exo(connect(vert_index_exo))
       enddo
	   
       call average2vec(wedgever(1,7),wedgever(1,1),
     & wedgever(1,2))
       call average2vec(wedgever(1,8),wedgever(1,2),
     & wedgever(1,3))
       call average2vec(wedgever(1,9),wedgever(1,1),
     & wedgever(1,3))
       call average2vec(wedgever(1,10),wedgever(1,1),
     & wedgever(1,4))
       call average2vec(wedgever(1,11),wedgever(1,2),
     & wedgever(1,5))
       call average2vec(wedgever(1,12),wedgever(1,3),
     & wedgever(1,6))
       call average2vec(wedgever(1,13),wedgever(1,4),
     & wedgever(1,5))
       call average2vec(wedgever(1,14),wedgever(1,5),
     & wedgever(1,6))
       call average2vec(wedgever(1,15),wedgever(1,4),
     & wedgever(1,6))

cc assign sideset to wedge elements
       call rzero(wedgess,5)
       if (num_side_sets.ne.0) then
        do ifc_exo=1,5
          wedgess(ifc_exo) = exoss(ifc_exo,iel_exo)
        enddo
       endif
	   
cc given wedge15 vertices, and return you four hex coords.
       call wedgetohex(hexver,wedgever,hexss,wedgess)	   
       
       do ihex = 1,3
          iel_nek = iel_nek + 1
          do inekvert = 1,27
             xm1(inekvert,1,1,iel_nek) = hexver(1,inekvert,ihex)
             ym1(inekvert,1,1,iel_nek) = hexver(2,inekvert,ihex)
             zm1(inekvert,1,1,iel_nek) = hexver(3,inekvert,ihex)
          enddo
          do ifc=1,6
             if(hexss(ifc,ihex).gt.0) then
              bc(5,ifc,iel_nek) = hexss(ifc,ihex)
              cbc(ifc,iel_nek)   = 'EXO' ! dummy boundary condition
             endif
          enddo
       enddo

      endif ! end if if (iel_exo.le.tetnumber) then
	   
      enddo

      if(iel_nek.gt.max_num_elem) then
       write(6,*) 'ERROR, increase MAXNEL to ',iel_nek
       write(6,*) 'and recompile exo2nek'
       call exit(1)
      endif
	  
      write(6,*) 'Converted elements in nek:',iel_nek
      write(6,'(A)') 'Done :: Converting elements '
	  
      num_elem = iel_nek

      if(num_elem.gt.max_num_elem) then
       write(6,*) 'ERROR, increase MAXNEL to ',num_elem
       write(6,*) 'and recompile exo2nek'
       call exit(1)
      endif

      return
      end
C--------------------------------------------------------------------
      subroutine split_convert2
c split 1 TETRA10 to 4 HEX20, block 1
c split 1 exo HEX20 to 8 hex20, block 2
C split 1 WEDGE15 to 6 HEX20, block 3
C
c  Subroutine to convert an already read exodusII mesh to a nek
c  mesh. The idea is to fill each element's node coordinates of
c  size lx1**3 (3D) or lx1**2 (2D) (lx1=3) with the hex27/quad9
c  coordinates.
c
      include 'exodusII.inc'
#     include "SIZE"

cc in fortran arraysm higher order is at back.
      real*8 hexver(3,27,8) ! hex coordinates
      real*8 tetver(3,10),wedgever(3,15)
      real*8 ehexver(3,20)

cc exoss is used to store sideset information for all exo elements.
cc tet,hex,wedge
      integer exoss(6,max_num_elem) 

      integer tetss(4),wedgess(5),ehexss(6),hexss(6,8)
      integer ehexnumber,tetnumber,wedgenumber
      integer vert_index_exo
      integer iel_nek,iel_exo,ifc_exo

cc Haomin Yuan, I do not know why, rzero cannot used for exoss
      do i=1,6*max_num_elem
      exoss(i,1)=0
      enddo

c store sideset information

      if (num_side_sets.ne.0) then
      write(6,'(a)') ''
      write(6,'(a)') 'Store SideSet information from EXO file'
        do iss=1,num_side_sets   ! loop over ss 
        write(6,'(a,i2,a)') 'Sideset ',idss(iss), ' ...'
           do i=1,num_sides_in_set(iss) ! loop over sides in ss
             iel_exo = elem_list(i,iss) ! exo element number
             ifc_exo = side_list(i,iss) ! exo element face number
             exoss(ifc_exo,iel_exo) = idss(iss)
           enddo
        enddo
      endif
c
      write(6,'(A)') ' '
      write(6,'(A)') 'Converting elements ... '
	  
      call blank   (cbc,3*2*ldim*max_num_elem)
      call rzero   (bc,5*2*ldim*max_num_elem)
      call blank   (ccurve,(4+8*(ldim-2))*max_num_elem)
      call rzero   (curve,2*ldim*12*max_num_elem)

      call rzero(hexss,6*8) 
 
      tetnumber  = num_elem_in_block(1) ! all exo tet elements in block 1 (or first block)
      ehexnumber = num_elem_in_block(2) ! all exo hex elements in block 2 (or second block)
      wedgenumber = num_elem_in_block(3) ! all exo wedge elements in block 3 (or third block)

      vert_index_exo = 0
      iel_nek = 0
	  
      do iel_exo = 1, num_elem  
	  
!! number of elements == 4*number of tet + 6*number of wedge + 8*number of hex

cc block 1, if tet.
       if (iel_exo.le.tetnumber) then
cc now convert 1 tet10 element to 4 hex20 elements.
c       do ivert = 1, 10
c       vert_index_exo = vert_index_exo + 1  
c       tetver(1,ivert) = x_exo(connect(vert_index_exo))
c       tetver(2,ivert) = y_exo(connect(vert_index_exo))
c       tetver(3,ivert) = z_exo(connect(vert_index_exo))
c       enddo

cc read tet 4 . 
cc linear interpolate to tet10
       do ivert = 1, 4
       vert_index_exo = vert_index_exo + 1  
       tetver(1,ivert) = x_exo(connect(vert_index_exo))
       tetver(2,ivert) = y_exo(connect(vert_index_exo))
       tetver(3,ivert) = z_exo(connect(vert_index_exo))
       enddo

       call average2vec(tetver(1,5),tetver(1,1),tetver(1,2))
       call average2vec(tetver(1,6),tetver(1,2),tetver(1,3))
       call average2vec(tetver(1,7),tetver(1,1),tetver(1,3))
       call average2vec(tetver(1,8),tetver(1,1),tetver(1,4))
       call average2vec(tetver(1,9),tetver(1,2),tetver(1,4))
       call average2vec(tetver(1,10),tetver(1,3),tetver(1,4))
	   
	   
cc assign sideset to tet elements
       call rzero(tetss,4)
       if (num_side_sets.ne.0) then
        do ifc_exo=1,4
          tetss(ifc_exo) = exoss(ifc_exo,iel_exo)
        enddo
       endif

cc given tet10 vertices, and return you four hex coords.
       call tettohex(hexver,tetver,hexss,tetss) 

       do ihex = 1,4
          iel_nek = iel_nek + 1
          do inekvert = 1,27
             xm1(inekvert,1,1,iel_nek) = hexver(1,inekvert,ihex)
             ym1(inekvert,1,1,iel_nek) = hexver(2,inekvert,ihex)
             zm1(inekvert,1,1,iel_nek) = hexver(3,inekvert,ihex)
          enddo
          do ifc=1,6
             if(hexss(ifc,ihex).gt.0) then
              bc(5,ifc,iel_nek) = hexss(ifc,ihex)
              cbc(ifc,iel_nek)   = 'EXO' ! dummy boundary condition
             endif
          enddo
       enddo

cc block 2, if exo hex element, convert to 8 nek hex elements
       elseif (iel_exo.gt.tetnumber
     &.and.(iel_exo.le.(tetnumber+ehexnumber))) then

c       do ivert = 1,20
c       vert_index_exo = vert_index_exo + 1  
c       ehexver(1,ivert) = x_exo(connect(vert_index_exo))
c       ehexver(2,ivert) = y_exo(connect(vert_index_exo))
c       ehexver(3,ivert) = z_exo(connect(vert_index_exo))
c       enddo

cc read hex8
cc linear interpolate to hex20	   
       do ivert = 1,8
       vert_index_exo = vert_index_exo + 1  
       ehexver(1,ivert) = x_exo(connect(vert_index_exo))
       ehexver(2,ivert) = y_exo(connect(vert_index_exo))
       ehexver(3,ivert) = z_exo(connect(vert_index_exo))
       enddo

       call average2vec(ehexver(1,9),ehexver(1,1),ehexver(1,2))
       call average2vec(ehexver(1,10),ehexver(1,2),ehexver(1,3))
       call average2vec(ehexver(1,11),ehexver(1,3),ehexver(1,4))
       call average2vec(ehexver(1,12),ehexver(1,1),ehexver(1,4))
       call average2vec(ehexver(1,13),ehexver(1,1),ehexver(1,5))
       call average2vec(ehexver(1,14),ehexver(1,2),ehexver(1,6))
       call average2vec(ehexver(1,15),ehexver(1,3),ehexver(1,7))
       call average2vec(ehexver(1,16),ehexver(1,4),ehexver(1,8))
       call average2vec(ehexver(1,17),ehexver(1,5),ehexver(1,6))
       call average2vec(ehexver(1,18),ehexver(1,6),ehexver(1,7))
       call average2vec(ehexver(1,19),ehexver(1,7),ehexver(1,8))
       call average2vec(ehexver(1,20),ehexver(1,5),ehexver(1,8))

cc assign sideset
       call rzero(ehexss,6)
       if (num_side_sets.ne.0) then
        do ifc_exo=1,6
          ehexss(ifc_exo) = exoss(ifc_exo,iel_exo)
        enddo
       endif

       call ehexto8hex(hexver,ehexver,hexss,ehexss) 

       do ihex = 1,8
          iel_nek = iel_nek + 1
          do inekvert = 1,27
             xm1(inekvert,1,1,iel_nek) = hexver(1,inekvert,ihex)
             ym1(inekvert,1,1,iel_nek) = hexver(2,inekvert,ihex)
             zm1(inekvert,1,1,iel_nek) = hexver(3,inekvert,ihex)
          enddo
          do ifc=1,6
             if(hexss(ifc,ihex).gt.0) then
              bc(5,ifc,iel_nek) = hexss(ifc,ihex)
              cbc(ifc,iel_nek)   = 'EXO' ! dummy boundary condition
             endif
          enddo
       enddo

! block 3, if WEDGE15 elements
      elseif (iel_exo.gt.(tetnumber+ehexnumber)) then
	  
c       do ivert = 1,15
c       vert_index_exo = vert_index_exo + 1
c       wedgever(1,ivert) = x_exo(connect(vert_index_exo))
c       wedgever(2,ivert) = y_exo(connect(vert_index_exo))
c       wedgever(3,ivert) = z_exo(connect(vert_index_exo))
c       enddo

cc read wedge 6,
cc linear interpolate to wedge 15
       do ivert = 1,6
       vert_index_exo = vert_index_exo + 1  
       wedgever(1,ivert) = x_exo(connect(vert_index_exo))
       wedgever(2,ivert) = y_exo(connect(vert_index_exo))
       wedgever(3,ivert) = z_exo(connect(vert_index_exo))
       enddo
	   
       call average2vec(wedgever(1,7),wedgever(1,1),
     & wedgever(1,2))
       call average2vec(wedgever(1,8),wedgever(1,2),
     & wedgever(1,3))
       call average2vec(wedgever(1,9),wedgever(1,1),
     & wedgever(1,3))
       call average2vec(wedgever(1,10),wedgever(1,1),
     & wedgever(1,4))
       call average2vec(wedgever(1,11),wedgever(1,2),
     & wedgever(1,5))
       call average2vec(wedgever(1,12),wedgever(1,3),
     & wedgever(1,6))
       call average2vec(wedgever(1,13),wedgever(1,4),
     & wedgever(1,5))
       call average2vec(wedgever(1,14),wedgever(1,5),
     & wedgever(1,6))
       call average2vec(wedgever(1,15),wedgever(1,4),
     & wedgever(1,6))


cc assign sideset to wedge elements
       call rzero(wedgess,5)
       if (num_side_sets.ne.0) then
        do ifc_exo=1,5
          wedgess(ifc_exo) = exoss(ifc_exo,iel_exo)
        enddo
       endif

cc given wedge15 vertices, and return you 6 hex coords.
       call wedgetohex2(hexver,wedgever,hexss,wedgess)	   

       do ihex = 1,6
          iel_nek = iel_nek + 1
          do inekvert = 1,27
             xm1(inekvert,1,1,iel_nek) = hexver(1,inekvert,ihex)
             ym1(inekvert,1,1,iel_nek) = hexver(2,inekvert,ihex)
             zm1(inekvert,1,1,iel_nek) = hexver(3,inekvert,ihex)
          enddo
          do ifc=1,6
             if(hexss(ifc,ihex).gt.0) then
              bc(5,ifc,iel_nek) = hexss(ifc,ihex)
              cbc(ifc,iel_nek)   = 'EXO' ! dummy boundary condition
             endif
          enddo
       enddo

      endif

      enddo ! do iel_exo = 1, num_elem

      if(iel_nek.gt.max_num_elem) then
       write(6,*) 'ERROR, increase MAXNEL to ',iel_nek
       write(6,*) 'and recompile exo2nek'
       call exit(1)
      endif
	  
      write(6,*) 'Converted elements in nek:',iel_nek
      write(6,'(A)') 'Done :: Converting elements '

      num_elem = iel_nek
	  
      if(num_elem.gt.max_num_elem) then
       write(6,*) 'ERROR, increase MAXNEL to ',num_elem
       write(6,*) 'and recompile exo2nek'
       call exit(1)
      endif
	  
      return
      end
C--------------------------------------------------------------------
      subroutine tettohex(hexver,tetver,hexss,tetss)
      real*8 tetver(3,10) ! tet vertices
      real*8 tetface(3,4) ! tet face center
      real*8 tetcen(3,1)  ! tet vol center
      real*8 hex8(3,8)
      real*8 hexver(3,27,4) ! four hex vertices in nek format.
      real*8 tempvec(3,3) ! temperary vector variable

      integer tetss(4),hexss(6,4)
	  
      do i=1,6*4
      hexss(i,1)=0
      enddo

cc get face center coords.
      call average3vec(tetface(1,1),tetver(1,1),tetver(1,2),tetver(1,4))
      call average3vec(tetface(1,2),tetver(1,2),tetver(1,3),tetver(1,4))
      call average3vec(tetface(1,3),tetver(1,1),tetver(1,3),tetver(1,4))
      call average3vec(tetface(1,4),tetver(1,1),tetver(1,2),tetver(1,3))

cc get tet vol center
      call average4vec(tetcen(1,1),tetver(1,1),tetver(1,2),
     &tetver(1,3),tetver(1,4))

cc assign coordinates to four hex.
cc hex 1
      hexss(1,1) = tetss(1)
      hexss(4,1) = tetss(3)
      hexss(5,1) = tetss(4)

      call assignvec(hex8(1,1),tetver(1,1))
      call assignvec(hex8(1,2),tetver(1,5))
      call assignvec(hex8(1,3),tetface(1,4))
      call assignvec(hex8(1,4),tetver(1,7))
      call assignvec(hex8(1,5),tetver(1,8))
      call assignvec(hex8(1,6),tetface(1,1))
      call assignvec(hex8(1,7),tetcen(1,1))
      call assignvec(hex8(1,8),tetface(1,3))
      call hex8tohex27(hexver(1,1,1),hex8(1,1))	 
cc hex 2
      hexss(1,2) = tetss(1)
      hexss(2,2) = tetss(2)
      hexss(5,2) = tetss(4)

      call assignvec(hex8(1,1),tetver(1,5))
      call assignvec(hex8(1,2),tetver(1,2))
      call assignvec(hex8(1,3),tetver(1,6))
      call assignvec(hex8(1,4),tetface(1,4))
      call assignvec(hex8(1,5),tetface(1,1))
      call assignvec(hex8(1,6),tetver(1,9))
      call assignvec(hex8(1,7),tetface(1,2))
      call assignvec(hex8(1,8),tetcen(1,1))
      call hex8tohex27(hexver(1,1,2),hex8(1,1))
cc hex 3
      hexss(2,3) = tetss(2)
      hexss(3,3) = tetss(3)
      hexss(5,3) = tetss(4)

      call assignvec(hex8(1,1),tetface(1,4))
      call assignvec(hex8(1,2),tetver(1,6))
      call assignvec(hex8(1,3),tetver(1,3))
      call assignvec(hex8(1,4),tetver(1,7))
      call assignvec(hex8(1,5),tetcen(1,1))
      call assignvec(hex8(1,6),tetface(1,2))
      call assignvec(hex8(1,7),tetver(1,10))
      call assignvec(hex8(1,8),tetface(1,3))
      call hex8tohex27(hexver(1,1,3),hex8(1,1))

cc hex 4
      hexss(2,4) = tetss(2)
      hexss(3,4) = tetss(3)
      hexss(6,4) = tetss(1)

      call assignvec(hex8(1,1),tetcen(1,1))
      call assignvec(hex8(1,2),tetface(1,2))
      call assignvec(hex8(1,3),tetver(1,10))
      call assignvec(hex8(1,4),tetface(1,3))
      call assignvec(hex8(1,5),tetface(1,1))
      call assignvec(hex8(1,6),tetver(1,9))
      call assignvec(hex8(1,7),tetver(1,4))
      call assignvec(hex8(1,8),tetver(1,8))
      call hex8tohex27(hexver(1,1,4),hex8(1,1))

      return
      end	  
C--------------------------------------------------------------------
      subroutine wedgetohex(hexver,wedgever,hexss,wedgess)
cc convert 1 wedge to 3 hex
      real*8 wedgever(3,15) ! tet vertices
      real*8 wedgeface(3,5) ! tet face center
      real*8 wedgecen(3,1) ! tet vol center
      real*8 hex8(3,8)
      real*8 hexver(3,27,4) ! four hex vertices in nek format.
      real*8 tempvec(3,3) ! temperary vector variable

      integer wedgess(5),hexss(6,4)
	  
      do i=1,6*4
      hexss(i,1)=0
      enddo

      call average3vec(wedgeface(1,4),wedgever(1,1),wedgever(1,2),
     &wedgever(1,3))
      call average3vec(wedgeface(1,5),wedgever(1,4),wedgever(1,5),
     &wedgever(1,6))

cc assign coordinates to 3 hex.
cc hex 1
      hexss(1,1) = wedgess(1)
      hexss(4,1) = wedgess(3)
      hexss(5,1) = wedgess(4)
      hexss(6,1) = wedgess(5)

      call assignvec(hex8(1,1),wedgever(1,1))
      call assignvec(hex8(1,2),wedgever(1,7))
      call assignvec(hex8(1,3),wedgeface(1,4))
      call assignvec(hex8(1,4),wedgever(1,9))
      call assignvec(hex8(1,5),wedgever(1,4))
      call assignvec(hex8(1,6),wedgever(1,13))
      call assignvec(hex8(1,7),wedgeface(1,5))
      call assignvec(hex8(1,8),wedgever(1,15))
      call hex8tohex27(hexver(1,1,1),hex8(1,1))	 
cc hex 2
      hexss(1,2) = wedgess(1)
      hexss(2,2) = wedgess(2)
      hexss(5,2) = wedgess(4)
      hexss(6,2) = wedgess(5)

      call assignvec(hex8(1,1),wedgever(1,7))
      call assignvec(hex8(1,2),wedgever(1,2))
      call assignvec(hex8(1,3),wedgever(1,8))
      call assignvec(hex8(1,4),wedgeface(1,4))
      call assignvec(hex8(1,5),wedgever(1,13))
      call assignvec(hex8(1,6),wedgever(1,5))
      call assignvec(hex8(1,7),wedgever(1,14))
      call assignvec(hex8(1,8),wedgeface(1,5))
      call hex8tohex27(hexver(1,1,2),hex8(1,1))
cc hex 3
      hexss(2,3) = wedgess(2)
      hexss(3,3) = wedgess(3)
      hexss(5,3) = wedgess(4)
      hexss(6,3) = wedgess(5)

      call assignvec(hex8(1,1),wedgeface(1,4))
      call assignvec(hex8(1,2),wedgever(1,8))
      call assignvec(hex8(1,3),wedgever(1,3))
      call assignvec(hex8(1,4),wedgever(1,9))
      call assignvec(hex8(1,5),wedgeface(1,5))
      call assignvec(hex8(1,6),wedgever(1,14))
      call assignvec(hex8(1,7),wedgever(1,6))
      call assignvec(hex8(1,8),wedgever(1,15))
      call hex8tohex27(hexver(1,1,3),hex8(1,1))

      return
      end	  
c--------------------------------------------------------------------
C--------------------------------------------------------------------
      subroutine wedgetohex2(hexver,wedgever,hexss,wedgess)
cc convert 1 wedge to 6 nek hex20 elements
      real*8 wedgever(3,15) ! tet vertices
      real*8 wedgeface(3,5) ! tet face center
      real*8 wedgecen(3,1) ! tet vol center
      real*8 hex8(3,8)
      real*8 hexver(3,27,8) ! four hex vertices in nek format.
      real*8 tempvec(3,3) ! temperary vector variable

      integer wedgess(5),hexss(6,8)

      do i=1,6*8
      hexss(i,1)=0
      enddo

      call average4vec(wedgeface(1,1),wedgever(1,1),wedgever(1,2),
     &wedgever(1,4),wedgever(1,5))	 
      call average4vec(wedgeface(1,2),wedgever(1,2),wedgever(1,3),
     &wedgever(1,5),wedgever(1,6))
      call average4vec(wedgeface(1,3),wedgever(1,1),wedgever(1,3),
     &wedgever(1,4),wedgever(1,6))
      call average3vec(wedgeface(1,4),wedgever(1,1),wedgever(1,2),
     &wedgever(1,3))
      call average3vec(wedgeface(1,5),wedgever(1,4),wedgever(1,5),
     &wedgever(1,6))
      call average2vec(wedgecen(1,1),wedgeface(1,4),wedgeface(1,5))

cc assign coordinates to 6 hex.
cc hex 1
      hexss(1,1) = wedgess(1)
      hexss(4,1) = wedgess(3)
      hexss(5,1) = wedgess(4)

      call assignvec(hex8(1,1),wedgever(1,1))
      call assignvec(hex8(1,2),wedgever(1,7))
      call assignvec(hex8(1,3),wedgeface(1,4))
      call assignvec(hex8(1,4),wedgever(1,9))
      call assignvec(hex8(1,5),wedgever(1,10))
      call assignvec(hex8(1,6),wedgeface(1,1))
      call assignvec(hex8(1,7),wedgecen(1,1))
      call assignvec(hex8(1,8),wedgeface(1,3))
      call hex8tohex27(hexver(1,1,1),hex8(1,1))

cc hex 2
      hexss(1,2) = wedgess(1)
      hexss(2,2) = wedgess(2)
      hexss(5,2) = wedgess(4)

      call assignvec(hex8(1,1),wedgever(1,7))
      call assignvec(hex8(1,2),wedgever(1,2))
      call assignvec(hex8(1,3),wedgever(1,8))
      call assignvec(hex8(1,4),wedgeface(1,4))
      call assignvec(hex8(1,5),wedgeface(1,1))
      call assignvec(hex8(1,6),wedgever(1,11))
      call assignvec(hex8(1,7),wedgeface(1,2))
      call assignvec(hex8(1,8),wedgecen(1,1))
      call hex8tohex27(hexver(1,1,2),hex8(1,1))
cc hex 3
      hexss(2,3) = wedgess(2)
      hexss(3,3) = wedgess(3)
      hexss(5,3) = wedgess(4)

      call assignvec(hex8(1,1),wedgeface(1,4))
      call assignvec(hex8(1,2),wedgever(1,8))
      call assignvec(hex8(1,3),wedgever(1,3))
      call assignvec(hex8(1,4),wedgever(1,9))
      call assignvec(hex8(1,5),wedgecen(1,1))
      call assignvec(hex8(1,6),wedgeface(1,2))
      call assignvec(hex8(1,7),wedgever(1,12))
      call assignvec(hex8(1,8),wedgeface(1,3))
      call hex8tohex27(hexver(1,1,3),hex8(1,1))

cc hex 4
      hexss(1,4) = wedgess(1)
      hexss(4,4) = wedgess(3)
      hexss(6,4) = wedgess(5)

      call assignvec(hex8(1,1),wedgever(1,10))
      call assignvec(hex8(1,2),wedgeface(1,1))
      call assignvec(hex8(1,3),wedgecen(1,1))
      call assignvec(hex8(1,4),wedgeface(1,3))
      call assignvec(hex8(1,5),wedgever(1,4))
      call assignvec(hex8(1,6),wedgever(1,13))
      call assignvec(hex8(1,7),wedgeface(1,5))
      call assignvec(hex8(1,8),wedgever(1,15))
      call hex8tohex27(hexver(1,1,4),hex8(1,1))	 
cc hex 5
      hexss(1,5) = wedgess(1)
      hexss(2,5) = wedgess(2)
      hexss(6,5) = wedgess(5)

      call assignvec(hex8(1,1),wedgeface(1,1))
      call assignvec(hex8(1,2),wedgever(1,11))
      call assignvec(hex8(1,3),wedgeface(1,2))
      call assignvec(hex8(1,4),wedgecen(1,1))
      call assignvec(hex8(1,5),wedgever(1,13))
      call assignvec(hex8(1,6),wedgever(1,5))
      call assignvec(hex8(1,7),wedgever(1,14))
      call assignvec(hex8(1,8),wedgeface(1,5))
      call hex8tohex27(hexver(1,1,5),hex8(1,1))
cc hex 6
      hexss(2,6) = wedgess(2)
      hexss(3,6) = wedgess(3)
      hexss(6,6) = wedgess(5)

      call assignvec(hex8(1,1),wedgecen(1,1))
      call assignvec(hex8(1,2),wedgeface(1,2))
      call assignvec(hex8(1,3),wedgever(1,12))
      call assignvec(hex8(1,4),wedgeface(1,3))
      call assignvec(hex8(1,5),wedgeface(1,5))
      call assignvec(hex8(1,6),wedgever(1,14))
      call assignvec(hex8(1,7),wedgever(1,6))
      call assignvec(hex8(1,8),wedgever(1,15))
      call hex8tohex27(hexver(1,1,6),hex8(1,1))

      return
      end	  
c--------------------------------------------------------------------
      subroutine ehexto8hex(hexver,ehexver,hexss,ehexss)
cc convert exodus hex20 elements 
cc to 8 nek hex20 elements
      real*8 ehexver(3,20)
      real*8 ehexface(3,6)
      real*8 ehexcen(3,1) ! tet vol center
      real*8 hex8(3,8)
      real*8 hexver(3,27,8)
      integer ehexss(6),hexss(6,8)

      do i=1,6*8
      hexss(i,1)=0
      enddo

      call average4vec(ehexface(1,1),ehexver(1,1),ehexver(1,2),
     &ehexver(1,5),ehexver(1,6))
      call average4vec(ehexface(1,2),ehexver(1,2),ehexver(1,3),
     &ehexver(1,6),ehexver(1,7))
      call average4vec(ehexface(1,3),ehexver(1,3),ehexver(1,4),
     &ehexver(1,7),ehexver(1,8))
      call average4vec(ehexface(1,4),ehexver(1,1),ehexver(1,4),
     &ehexver(1,5),ehexver(1,8))
      call average4vec(ehexface(1,5),ehexver(1,1),ehexver(1,2),
     &ehexver(1,3),ehexver(1,4))
      call average4vec(ehexface(1,6),ehexver(1,5),ehexver(1,6),
     &ehexver(1,7),ehexver(1,8))
      call average2vec(ehexcen(1,1),ehexface(1,5),ehexface(1,6))

cc assign coordinates to 8 hex.
cc hex 1
      hexss(1,1) = ehexss(1)
      hexss(4,1) = ehexss(4)
      hexss(5,1) = ehexss(5)

      call assignvec(hex8(1,1),ehexver(1,1))
      call assignvec(hex8(1,2),ehexver(1,9))
      call assignvec(hex8(1,3),ehexface(1,5))
      call assignvec(hex8(1,4),ehexver(1,12))
      call assignvec(hex8(1,5),ehexver(1,13))
      call assignvec(hex8(1,6),ehexface(1,1))
      call assignvec(hex8(1,7),ehexcen(1,1))
      call assignvec(hex8(1,8),ehexface(1,4))
      call hex8tohex27(hexver(1,1,1),hex8(1,1))	 

cc hex 2
      hexss(1,2) = ehexss(1)
      hexss(2,2) = ehexss(2)
      hexss(5,2) = ehexss(5)

      call assignvec(hex8(1,1),ehexver(1,9))
      call assignvec(hex8(1,2),ehexver(1,2))
      call assignvec(hex8(1,3),ehexver(1,10))
      call assignvec(hex8(1,4),ehexface(1,5))
      call assignvec(hex8(1,5),ehexface(1,1))
      call assignvec(hex8(1,6),ehexver(1,14))
      call assignvec(hex8(1,7),ehexface(1,2))
      call assignvec(hex8(1,8),ehexcen(1,1))
      call hex8tohex27(hexver(1,1,2),hex8(1,1))
cc hex 3
      hexss(2,3) = ehexss(2)
      hexss(3,3) = ehexss(3)
      hexss(5,3) = ehexss(5)

      call assignvec(hex8(1,1),ehexface(1,5))
      call assignvec(hex8(1,2),ehexver(1,10))
      call assignvec(hex8(1,3),ehexver(1,3))
      call assignvec(hex8(1,4),ehexver(1,11))
      call assignvec(hex8(1,5),ehexcen(1,1))
      call assignvec(hex8(1,6),ehexface(1,2))
      call assignvec(hex8(1,7),ehexver(1,15))
      call assignvec(hex8(1,8),ehexface(1,3))
      call hex8tohex27(hexver(1,1,3),hex8(1,1))

cc hex 4
      hexss(3,4) = ehexss(3)
      hexss(4,4) = ehexss(4)
      hexss(5,4) = ehexss(5)

      call assignvec(hex8(1,1),ehexver(1,12))
      call assignvec(hex8(1,2),ehexface(1,5))
      call assignvec(hex8(1,3),ehexver(1,11))
      call assignvec(hex8(1,4),ehexver(1,4))
      call assignvec(hex8(1,5),ehexface(1,4))
      call assignvec(hex8(1,6),ehexcen(1,1))
      call assignvec(hex8(1,7),ehexface(1,3))
      call assignvec(hex8(1,8),ehexver(1,16))
      call hex8tohex27(hexver(1,1,4),hex8(1,1))

cc hex 5
      hexss(1,5) = ehexss(1)
      hexss(4,5) = ehexss(4)
      hexss(6,5) = ehexss(6)

      call assignvec(hex8(1,1),ehexver(1,13))
      call assignvec(hex8(1,2),ehexface(1,1))
      call assignvec(hex8(1,3),ehexcen(1,1))
      call assignvec(hex8(1,4),ehexface(1,4))
      call assignvec(hex8(1,5),ehexver(1,5))
      call assignvec(hex8(1,6),ehexver(1,17))
      call assignvec(hex8(1,7),ehexface(1,6))
      call assignvec(hex8(1,8),ehexver(1,20))
      call hex8tohex27(hexver(1,1,5),hex8(1,1))	 
cc hex 6
      hexss(1,6) = ehexss(1)
      hexss(2,6) = ehexss(2)
      hexss(6,6) = ehexss(6)

      call assignvec(hex8(1,1),ehexface(1,1))
      call assignvec(hex8(1,2),ehexver(1,14))
      call assignvec(hex8(1,3),ehexface(1,2))
      call assignvec(hex8(1,4),ehexcen(1,1))
      call assignvec(hex8(1,5),ehexver(1,17))
      call assignvec(hex8(1,6),ehexver(1,6))
      call assignvec(hex8(1,7),ehexver(1,18))
      call assignvec(hex8(1,8),ehexface(1,6))
      call hex8tohex27(hexver(1,1,6),hex8(1,1))
cc hex 7
      hexss(2,7) = ehexss(2)
      hexss(3,7) = ehexss(3)
      hexss(6,7) = ehexss(6)

      call assignvec(hex8(1,1),ehexcen(1,1))
      call assignvec(hex8(1,2),ehexface(1,2))
      call assignvec(hex8(1,3),ehexver(1,15))
      call assignvec(hex8(1,4),ehexface(1,3))
      call assignvec(hex8(1,5),ehexface(1,6))
      call assignvec(hex8(1,6),ehexver(1,18))
      call assignvec(hex8(1,7),ehexver(1,7))
      call assignvec(hex8(1,8),ehexver(1,19))
      call hex8tohex27(hexver(1,1,7),hex8(1,1))

cc hex 8
      hexss(3,8) = ehexss(3)
      hexss(4,8) = ehexss(4)
      hexss(6,8) = ehexss(6)

      call assignvec(hex8(1,1),ehexface(1,4))
      call assignvec(hex8(1,2),ehexcen(1,1))
      call assignvec(hex8(1,3),ehexface(1,3))
      call assignvec(hex8(1,4),ehexver(1,16))
      call assignvec(hex8(1,5),ehexver(1,20))
      call assignvec(hex8(1,6),ehexface(1,6))
      call assignvec(hex8(1,7),ehexver(1,19))
      call assignvec(hex8(1,8),ehexver(1,8))
      call hex8tohex27(hexver(1,1,8),hex8(1,1))	 

      return
      end
c--------------------------------------------------------------------
      subroutine hex8tohex27(hex27,hex8)
cc convert hex8 coordinates to hex27 coordinates in nek.
cc 
      real*8 hex8(3,8)
      real*8 hex27(3,27)
      real*8 tempvec(3,3)

cc  hex27 vert 1
      call assignvec(hex27(1,1),hex8(1,1))
cc  hex27 vert 2
      call average2vec(tempvec(1,1),hex8(1,1),hex8(1,2))
      call assignvec(hex27(1,2),tempvec(1,1))
cc  hex27 vert 3
      call assignvec(hex27(1,3),hex8(1,2))
cc  hex27 vert 4
      call average2vec(tempvec(1,1),hex8(1,1),hex8(1,4))
      call assignvec(hex27(1,4),tempvec(1,1))
cc  hex27 vert 5
      call average4vec(tempvec(1,1),hex8(1,1),
     &hex8(1,2),hex8(1,3),hex8(1,4))
      call assignvec(hex27(1,5),tempvec(1,1))
cc  hex27 vert 6
      call average2vec(tempvec(1,1),hex8(1,2),hex8(1,3))
      call assignvec(hex27(1,6),tempvec(1,1))
cc  hex27 vert 7
      call assignvec(hex27(1,7),hex8(1,4))
cc  hex27 vert 8
      call average2vec(tempvec(1,1),hex8(1,3),hex8(1,4))
      call assignvec(hex27(1,8),tempvec(1,1))
cc  hex27 vert 9
      call assignvec(hex27(1,9),hex8(1,3))

cc  hex27 vert 10
      call average2vec(tempvec(1,1),hex8(1,1),hex8(1,5))
      call assignvec(hex27(1,10),tempvec(1,1))  
cc  hex27 vert 11
      call average4vec(tempvec(1,1),hex8(1,1),
     &hex8(1,2),hex8(1,5),hex8(1,6))
      call assignvec(hex27(1,11),tempvec(1,1))
cc  hex27 vert 12
      call average2vec(tempvec(1,1),hex8(1,2),hex8(1,6))
      call assignvec(hex27(1,12),tempvec(1,1))

cc  hex27 vert 13
      call average4vec(tempvec(1,1),hex8(1,1),
     &hex8(1,4),hex8(1,5),hex8(1,8))
      call assignvec(hex27(1,13),tempvec(1,1))
cc  hex27 vert 14
      call average4vec(tempvec(1,1),hex8(1,1),
     &hex8(1,4),hex8(1,5),hex8(1,8))
      call average4vec(tempvec(1,2),hex8(1,2),
     &hex8(1,3),hex8(1,6),hex8(1,7))
      call average2vec(tempvec(1,3),tempvec(1,1),tempvec(1,2))
      call assignvec(hex27(1,14),tempvec(1,3))
cc  hex27 vert 15
      call average4vec(tempvec(1,1),hex8(1,2),
     &hex8(1,3),hex8(1,6),hex8(1,7))
      call assignvec(hex27(1,15),tempvec(1,1))

cc  hex27 vert 16
      call average2vec(tempvec(1,1),hex8(1,4),hex8(1,8))
      call assignvec(hex27(1,16),tempvec(1,1))
cc  hex27 vert 17
      call average4vec(tempvec(1,1),hex8(1,3),
     &hex8(1,4),hex8(1,7),hex8(1,8))
      call assignvec(hex27(1,17),tempvec(1,1))
cc  hex27 vert 18
      call average2vec(tempvec(1,1),hex8(1,3),hex8(1,7))
      call assignvec(hex27(1,18),tempvec(1,1))
	
cc  hex27 vert 19
      call assignvec(hex27(1,19),hex8(1,5))
cc  hex27 vert 20
      call average2vec(tempvec(1,1),hex8(1,5),hex8(1,6))
      call assignvec(hex27(1,20),tempvec(1,1))	  
cc  hex27 vert 21
      call assignvec(hex27(1,21),hex8(1,6))

cc  hex27 vert 22
      call average2vec(tempvec(1,1),hex8(1,5),hex8(1,8))
      call assignvec(hex27(1,22),tempvec(1,1))	  
cc  hex27 vert 23
      call average4vec(tempvec(1,1),hex8(1,5),
     &hex8(1,6),hex8(1,7),hex8(1,8))
      call assignvec(hex27(1,23),tempvec(1,1))
cc  hex27 vert 24
      call average2vec(tempvec(1,1),hex8(1,6),hex8(1,7))
      call assignvec(hex27(1,24),tempvec(1,1))
	  
cc  hex27 vert 25
      call assignvec(hex27(1,25),hex8(1,8))
cc  hex27 vert 26
      call average2vec(tempvec(1,1),hex8(1,7),hex8(1,8))
      call assignvec(hex27(1,26),tempvec(1,1))	  
cc  hex27 vert 27
      call assignvec(hex27(1,27),hex8(1,7))
	 
      return
      end
c--------------------------------------------------------------------
      subroutine assignvec(a,b)
      real*8 a(3),b(3)
      a(1) = b(1)
      a(2) = b(2)
      a(3) = b(3)
      return
      end
c--------------------------------------------------------------------
      subroutine average2vec(avg,a,b)
      real*8 a(3),b(3),avg(3)      
      avg(1) = (a(1)+b(1))*0.5
      avg(2) = (a(2)+b(2))*0.5
      avg(3) = (a(3)+b(3))*0.5
      return
      end
C--------------------------------------------------------------------
      subroutine average3vec(avg,a,b,c)
      real*8 a(3),b(3),c(3),avg(3)      
      avg(1) = (a(1)+b(1)+c(1))*0.33333333
      avg(2) = (a(2)+b(2)+c(2))*0.33333333
      avg(3) = (a(3)+b(3)+c(3))*0.33333333

      return
      end
c--------------------------------------------------------------------
      subroutine average4vec(avg,a,b,c,d)
      real*8  a(3),b(3),c(3),d(3),avg(3)       
      avg(1) = (a(1)+b(1)+c(1)+d(1))*0.25
      avg(2) = (a(2)+b(2)+c(2)+d(2))*0.25
      avg(3) = (a(3)+b(3)+c(3)+d(3))*0.25
      return
      end
c--------------------------------------------------------------------
C-----------------------------------------------------------------
      subroutine setbc
#     include "SIZE"
c  set boundary condition 
c  read from a file casename.bc

      parameter (npbc_max=10) ! maximum pairs of periodic boundary condition

      character*3 ubc
      integer tags(2),ibc,nbc,io,er
      integer ip,np
      integer ptags(2,npbc_max)
      real  pvecs(3,npbc_max)
	  
      character*32 bcname
      character*1 bcnam1(32)
      equivalence(bcname,bcnam1)
	  
      call blank (bcname,32)

      len = ltrunc(exoname,32)
      call chcopy(bcnam1,exoname,(len-3))
      call chcopy(bcnam1(len-3),'.bc' , 3)
      len = ltrunc(bcnam1,32)

      open(301,file=bcname,err=1010) ! if error, direct go to label 1010, return
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
          do ihex = 1, num_elem
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
        read(301,*) ptags(1,ip),ptags(2,ip),
     & pvecs(1,ip),pvecs(2,ip),pvecs(3,ip)
cc ptags(1,ip) is surface 1 number
cc ptags(2,ip) is surface 2 number
cc pvecs(1,ip),pvecs(2,ip),pvecs(3,ip) is user defined projection vector
      enddo
      close(301)

      write(6,*) 'Setting periodic boundary condition' 
      np = ip
c np is the pairs of periodic bc.
      do ip = 1,np
c mapping surface ptags(1,ip) to surface ptags(2,ip)
         call setPeriodic(ptags(1,ip),pvecs(1,ip),ptol)
      enddo

1010  return
      end
c--------------------------------------------------------------------
      subroutine setPeriodic(ptags,pvec,ptol)
#     include "SIZE"

      integer hex_face_node(4,6)
      data hex_face_node
     &      /1,3,21,19,3,9,27,21,7,9,27,25,1,7,25,19,
     &       1,7,9,3,19,21,27,25/

      integer parray(2,2,max_num_elem)
      real parea(2,max_num_elem)

      integer fnode(4),ifnode,ptags(2),ipe,nipe(2),nperror
      real pvec(3),ptol,fpxyz(3,2)
      real dist,distMax

      ipe = 0

c collect ihex,iface for surface ptags(1)
      do ihex = 1,num_elem
         do iface = 1,6
            if(bc(5,iface,ihex).eq.ptags(1)) then 
              ipe = ipe + 1
              parray(1,1,ipe) = ihex
              parray(2,1,ipe) = iface
            endif
         enddo
      enddo
      nipe(1) = ipe

c collect ihex,iface for surface ptags(2)
      ipe = 0
      do ihex = 1, num_elem
         do iface = 1,6
            if(bc(5,iface,ihex).eq.ptags(2)) then
                ipe = ipe + 1
                parray(1,2,ipe) = ihex
                parray(2,2,ipe) = iface
             endif
         enddo
      enddo
      nipe(2) = ipe

      write(6,*)'maping surface',ptags(1),'with',nipe(1),'faces'
      write(6,*)'to surface',ptags(2),'with',nipe(2),'faces' 

      if(nipe(1).ne.nipe(2))  then
         write(6,*) 'EORROR, face numbers are not matching'
         return
      endif

c 1st loop, loop faces on surface 1
      do ipe = 1,nipe(1)
         ihex = parray(1,1,ipe)
         iface = parray(2,1,ipe)
c get face center xyz
         call rzero(fpxyz(1,1),3)				
	     do ifnode = 1,4
             fnode(ifnode)=hex_face_node(ifnode,iface)
             fpxyz(1,1) = fpxyz(1,1)+xm1(fnode(ifnode),1,1,ihex)*0.25
             fpxyz(2,1) = fpxyz(2,1)+ym1(fnode(ifnode),1,1,ihex)*0.25
             fpxyz(3,1) = fpxyz(3,1)+zm1(fnode(ifnode),1,1,ihex)*0.25			 
         enddo

c 2nd,loop over surface 2
         distMax = 1000.0
         do ipe2 = 1,nipe(2)
                ihex2 = parray(1,2,ipe2)
                iface2 = parray(2,2,ipe2)
c get face center xyz
               call rzero(fpxyz(1,2),3)
               do ifnode = 1,4
               fnode(ifnode)=hex_face_node(ifnode,iface2)
               fpxyz(1,2) = fpxyz(1,2)+xm1(fnode(ifnode),1,1,ihex2)*0.25
               fpxyz(2,2) = fpxyz(2,2)+ym1(fnode(ifnode),1,1,ihex2)*0.25
               fpxyz(3,2) = fpxyz(3,2)+zm1(fnode(ifnode),1,1,ihex2)*0.25
              enddo
 
               dist = sqrt((fpxyz(1,2) - fpxyz(1,1) - pvec(1))**2
     & + (fpxyz(2,2) - fpxyz(2,1) - pvec(2))**2
     & + (fpxyz(3,2) - fpxyz(3,1) - pvec(3))**2)

               if (dist.lt.distMax) then 
                  distMax = dist
                  !write(6,*) distMax
                  if(distMax.le.ptol) then
                  bc(1,iface,ihex) = ihex2*1.0
                  bc(2,iface,ihex) = iface2*1.0
                  bc(1,iface2,ihex2) = ihex*1.0
                  bc(2,iface2,ihex2) = iface*1.0
                  cbc(iface,ihex) = 'P  '
                  cbc(iface2,ihex2) = 'P  '
cc for debug use only
cc          write(6,*) ihex,iface,bc(1,iface,ihex),bc(2,iface,ihex)
cc          write(6,*) ihex2,iface2,bc(1,iface2,ihex2),bc(2,iface2,ihex2)
cc          write(6,*) fpxyz(1,1),fpxyz(2,1),fpxyz(3,1)
cc          write(6,*) fpxyz(1,2),fpxyz(2,2),fpxyz(3,2)
cc          write(6,*) dist,areadiff,parea(1,ipe),parea(2,ipe2)
                  endif
               endif
         enddo
      enddo

          nperror = 0
		  
          write(6,*)'doing periodic check for surface',ptags(1)

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

      return
      end
c--------------------------------------------------------------------
C all blow part is direct takend from origin exo2nek.f, totally unchanged
C
c-----------------------------------------------------------------------
      subroutine exodus_read
c
c  Subroutine to read an exodusII binary file containing a mesh.
c  It uses exodus fortran binding subroutines, which depend on
c  the netcdf library for low level data access.
c
      include 'exodusII.inc'
#     include "SIZE"

      integer exoid, cpu_ws, io_ws

      character*(MXSTLN) typ, qa_record(4,10)
      character*(MXLNLN) titl
      character*1        cdum

      integer idblk              (max_num_elem_blk)
      integer num_attr           (max_num_elem_blk)  ! not used
      integer num_nodes_per_elem (max_num_elem_blk)
c
c open EXODUS II file
c
      cpu_ws = 8 ! use real*8 to communicate with exodus
      io_ws  = 0
      exoid  = exopen (exoname, EXREAD, cpu_ws, io_ws, vers, ierr)
      if (ierr.lt.0) then
        write(6,'(2a)') "ERROR: cannot open file ", exoname 
        STOP
      endif
      write(6,*)
      write(6,'(a32,a,f4.2)') 
     &      exoname," is an EXODUSII file; version ",vers
      write(6,'(a,i2)') "I/O word size", io_ws
c
c read database parameters
c
      call exgini (exoid, titl, num_dim, num_nodes, num_elem,
     &             num_elem_blk, num_node_sets, num_side_sets, ierr)
      if (ierr.lt.0) then
        write(6,'(a)') "ERROR: cannot read exodusII parameters (exgini)"
        STOP
      endif
      write (6,    '(/"database parameters:"/ /
     &               "title         = ", a81 / /
     &               "num_dim       = ", i8 /
     &               "num_nodes     = ", i8 /
     &               "num_elem      = ", i8 /
     &               "num_elem_blk  = ", i8 /
     &               "num_side_sets = ", i8)')
     &               titl,num_dim, num_nodes, num_elem,
     &               num_elem_blk, num_side_sets
      write (6,*)
c
c perform some checks
c
      if (num_elem.gt.max_num_elem) then
        write(6,*) 'Abort: number of elements too large',num_elem
        write(6,*) 'change MAXNEL and recompile'
        STOP
      endif
      if (num_side_sets.gt.max_num_sidesets) then
        write(6,'(a)')
     &    "ERROR: number of sidesets > max_num_sidesets. "
        write(6,'(a,i2,a)') "Set max_num_sidesets >= ",num_side_sets,
     &                      " and recompile exo2nek. "
        STOP
      endif
c
c read element block parameters
c
      call exgebi (exoid, idblk, ierr)
      if (ierr.lt.0) then
        write(6,'(a)') "ERROR: cannot read block ids (exgebi)"
        STOP
      endif

      do i = 1, num_elem_blk
        call exgelb (exoid, idblk(i), typ, num_elem_in_block(i),
     &               num_nodes_per_elem(i), num_attr(i), ierr)
        if (ierr.lt.0) then
          write(6,'(a,i3,a)')
     &    "ERROR: cannot read parameters for block ",i," (exgelb)"
          STOP
        endif
        write (6, '("element block id   = ", i8,/
     &              "element type       = ", a8,/
     &              "num_elem_in_block  = ", i8,/
     &              "num_nodes_per_elem = ", i8)')
     &              idblk(i), typ, num_elem_in_block(i),
     &              num_nodes_per_elem(i)
        if ( num_dim.eq.3.and.(num_nodes_per_elem(i).ne.27 .and.
     &                         num_nodes_per_elem(i).ne.20) ) then
          write(6,'(a)')
     &     "ERROR: Only HEX20/HEX27 elements are allowed in a 3D mesh!"
          write(6,'(a,i3)') "num_nodes_per_elem= ",num_nodes_per_elem(i)
          STOP
        elseif (num_dim.eq.2.and.(num_nodes_per_elem(i).ne.9 .and.
     &                            num_nodes_per_elem(i).ne.8) ) then
          write(6,'(a)')
     &      "ERROR: Only QUAD8/QUAD9 elements are allowed in a 2D mesh!"
          write(6,'(a,i3)') "num_nodes_per_elem= ",num_nodes_per_elem(i)
          STOP
        endif
        if (i.eq.1) nvert=num_nodes_per_elem(i)
        if (num_nodes_per_elem(i).ne.nvert) then
          write(6,'(a)') 
     &     "ERROR: All blocks should contain elements of the same type!"
          write(6,'(a,i3)') "num_nodes_per_elem= ",num_nodes_per_elem(i)
          write(6,'(a,i3)') "num_nodes_per_elem of block 1 ",nvert
          STOP
        endif
        write (6,*)
      enddo
c
c read nodal coordinates values from database
c
      call exgcor (exoid, x_exo, y_exo, z_exo, ierr)
      if (ierr.lt.0) then
        write(6,'(a)') "ERROR: cannot read nodal coordinates (exgcor)"
        STOP
      endif
c
c read element connectivity
c
      iend = 0
      do 60 i = 1, num_elem_blk
        istart = iend + 1
        call exgelc (exoid, idblk(i), connect(istart), ierr)
        !iend = num_nodes_per_elem(i)*num_elem_in_block(i)
        iend = iend + num_nodes_per_elem(i)*num_elem_in_block(i) ! kirk modifed!! 8/23/2018
        if (ierr.lt.0) then
          write(6,'(a)') "ERROR: cannot read elm. connectivity (exgelc)"
          STOP
        endif
60    continue
c
c read individual side sets
c
      if (num_side_sets .gt. 0) then
        call exgssi (exoid, idss, ierr)
        if (ierr.lt.0) then
          write(6,'(a)') "ERROR: cannot read SideSet ids (exgssi)"
          STOP
        endif

        do i = 1, num_side_sets
          call exgsp (exoid,idss(i),num_sides_in_set(i),idum,ierr)
          if (ierr.lt.0) then
            write(6,'(a,i3,a)')
     &      "ERROR: cannot read parameters for SideSet No.",i," (exgsp)"
            STOP
          endif
          write (6, '("side set ", i2, " num_sides = ", i8)')
     &           idss(i), num_sides_in_set(i)
          call exgss (exoid,idss(i),elem_list(1,i),side_list(1,i),ierr)
          if (ierr.lt.0) then
            write(6,'(a,i3,a)')
     &      "ERROR: cannot read parameters for SideSet No.",i," (exgss)"
            STOP
          endif
        enddo
        write (6,*)
      else
        write(6,'(a)') "WARNING: No SideSets in exodus file!"
      endif
c
c read QA records
c
      call exinq (exoid, EXQA, num_qa_rec, fdum, cdum, ierr)
      if (ierr.lt.0) then
        write(6,'(a)') "ERROR: cannot read QA records (exinq QA) "
        STOP
      elseif (ierr.gt.0) then
        write(6,'(a)') "INFO: file does not contain any QA records"
      else
        call exgqa (exoid, qa_record, ierr)
        if (ierr.lt.0) then
          write(6,'(a,i3)') "WARNING: cannot read QA records (exgqa)"
        else
          write (6, '("QA records = ")')
          if (num_qa_rec.gt.10) then
            write(6,'(2a)')
     &        'WARNING: Cannot handle more than 10 QA records',
     &        'Printing only the first 10...'
          else
            do i = 1, num_qa_rec
              do j = 1, 4
                write (6,'(a)') qa_record(j,i)
              enddo
            enddo
          endif
        endif
      endif

      return
      end
c--------------------------------------------------------------------
      subroutine convert
c
c  Subroutine to convert an already read exodusII mesh to a nek
c  mesh. The idea is to fill each element's node coordinates of
c  size lx1**3 (3D) or lx1**2 (2D) (lx1=3) with the hex27/quad9
c  coordinates.
c
      include 'exodusII.inc'
#     include "SIZE"

c node and face conversion (it works at least for cubit):
      integer exo_to_nek_vert3D(27)
      data    exo_to_nek_vert3D
     &      / 19,  1,  7, 25, 21,  3,  9, 27, 10                  ! hex27 to nek numbering
     &      ,  4, 16, 22, 20,  2,  8, 26, 12,  6
     &      , 18, 24, 14, 13, 15, 23,  5, 11, 17 /

      integer exo_to_nek_vert2D(9)
      data    exo_to_nek_vert2D   / 1, 3, 9, 7, 2, 6, 8, 4, 5  /  ! quad9 to nek numbering

      integer exo_to_nek_face3D(6)
      data    exo_to_nek_face3D  / 1, 5, 3, 6, 4, 2 /    ! symmetric face numbering

      integer exo_to_nek_face2D(4)
      data    exo_to_nek_face2D  / 1, 2, 3, 4 /          ! symmetric face numbering


c zero-out bc and curve sides arrays
      call blank   (cbc,3*2*ldim*max_num_elem)
      call rzero   (bc,5*2*ldim*max_num_elem)
      call blank   (ccurve,(4+8*(ldim-2))*max_num_elem)
      call rzero   (curve,2*ldim*12*max_num_elem)


      write(6,'(A)') ' '
      write(6,'(A)') 'Converting elements ... '
      do iel = 1, num_elem
        do ivert = 1, nvert
          if (num_dim.eq.2) then
            jvert = exo_to_nek_vert2D(ivert)
            xm1(jvert,1,1,iel)=x_exo(connect(nvert*(iel-1)+ivert))
            ym1(jvert,1,1,iel)=y_exo(connect(nvert*(iel-1)+ivert))
          else
            jvert = exo_to_nek_vert3D(ivert)
            xm1(jvert,1,1,iel)=x_exo(connect(nvert*(iel-1)+ivert))
            ym1(jvert,1,1,iel)=y_exo(connect(nvert*(iel-1)+ivert))
            zm1(jvert,1,1,iel)=z_exo(connect(nvert*(iel-1)+ivert))
          endif
        enddo
      enddo
      write(6,'(A)') 'done :: Converting elements '

      if (num_side_sets.ne.0) then
      write(6,'(a)') ''
      write(6,'(a)') 'Converting SideSets ...'
        do iss=1,num_side_sets   ! loop over ss 
        write(6,'(a)') ''
        write(6,'(a,i2,a)') 'Sideset ',idss(iss), ' ...'
           do i=1,num_sides_in_set(iss) ! loop over sides in ss
             iel = elem_list(i,iss)
             ifc = side_list(i,iss)
             jfc = exo_to_nek_face3D(ifc)
             cbc(jfc,iel)   = 'EXO' ! dummy exodus bc 
             bc (5,jfc,iel) = idss(iss)
           enddo
        enddo
      endif

c
c zero-out bc and curve sides arrays
c      call blank   (cbc,3*2*ldim*max_num_elem)
c      call rzero   (bc,5*2*ldim*max_num_elem)
c      call blank   (ccurve,(4+8*(ldim-2))*max_num_elem)
c      call rzero   (curve,2*ldim*12*max_num_elem)

c currently, this part is very unefficient.
c for large mesh, this part could take hours.
c 
c the old method do a full loop for each elements for each sideset, which is very slow.
c 
c now we do element loop first.
c
c set bc's
c
c      if (num_side_sets.eq.0) return   ! no sidesets
c
c      write(6,'(a)') ''
c      write(6,'(a)') 'Converting SideSets ...'
c the expensive part, improve it...
c      do iss=1,num_side_sets   ! loop over ss 
c        write(6,'(a)') ''
c        write(6,'(a,i2,a)') 'Sideset ',idss(iss), ' ...'
c        do iel=1,num_elem
c          do ifc=1,2*num_dim             ! loop over faces
c            do i=1,num_sides_in_set(iss) ! loop over sides in ss
c              if    ( (iel.eq.elem_list(i,iss))
c     &        .and. (ifc.eq.side_list(i,iss)) ) then
c                if (num_dim.eq.2) then
c                  jfc = exo_to_nek_face2D(ifc)
c                else
c                  jfc = exo_to_nek_face3D(ifc)
c                endif
c                cbc(jfc,iel)   = 'EXO' ! dummy exodus bc 
c                bc (5,jfc,iel) = idss(iss)
c              endif
c            enddo
c          enddo
c        enddo
c        write(6,'(A,I2)') 'done :: Sideset ',idss(iss)
c      enddo
c      write(6,'(a)') ''
c      write(6,'(a)') 'done :: Converting SideSets '

      return
      end
C--------------------------------------------------------------------
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

      if (num_side_sets.eq.0) return

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
          if (ch3.ne.'   ') then
c          if (ch3.eq.'EXO') then
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

