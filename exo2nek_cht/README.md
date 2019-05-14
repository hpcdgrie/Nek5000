This is a special version of exo2nek to for conjugate heat transfer mesh (cht)
All elements must be hex20.

note about usage:
1. fluid mesh blocks must be placed before solid mesh blocks, in terms of block ID.
(CUBIT can shift block order, not sure about ICEM)

2. sideset must be given to all boundary surfaces including fluid-solid interfaces.
However, NO sideset given to fluid-fluid interface and solid-solid interface.

3. example about sideset setup

block 1 (fluid): 
sideset 1 -> inlet
sideset 2 -> outlet
sideset 3 -> fluid-solid interface

block 2 (fluid):  
sideset 4 -> inlet
sideset 5 -> outlet
sideset 6 -> fluid-solid interface

block 3 (solid):  
sideset 7 -> insulate boundary
sideset 8 -> fluid-solid interface

bloci 4 (solid):  
sideset 9 -> insulate boundary
sideset 10 -> fluid-solid interface

to set actual boundary conditions
in usrdat2()

       do iel=1,nelt                       ! must use nelt
        do ifc=1,2*ndim 
          
          ! set boundary for velocity
          id_face1 = bc(5,ifc,iel,1)       
          if (id_face1.eq.1) then 
             cbc(ifc,iel,1) = 'v  '
          elseif (id_face1.eq.2) then
             cbc(ifc,iel,1) = 'O  '
          elseif (id_face1.eq.3) then
             cbc(ifc,iel,1) = 'W  '
	      elseif (id_face1.eq.4) then
             cbc(ifc,iel,1) = 'v  '
          elseif (id_face1.eq.5) then
             cbc(ifc,iel,1) = 'O  '
          elseif (id_face1.eq.6) then
             cbc(ifc,iel,1) = 'W  '
          endif
          
          ! set boundary for thermal field, should skip fluid-solid interface
          id_face2 = bc(5,ifc,iel,2)
          if (id_face2.eq.1) then 
             cbc(ifc,iel,1) = 't  '
          elseif (id_face2.eq.2) then
             cbc(ifc,iel,1) = 'O  '
	      elseif (id_face2.eq.4) then
             cbc(ifc,iel,1) = 't  '
          elseif (id_face2.eq.5) then
             cbc(ifc,iel,1) = 'O  '
          elseif (id_face2.eq.7) then
             cbc(ifc,iel,1) = 'I  '
          elseif (id_face2.eq.9) then
             cbc(ifc,iel,1) = 'I  '
          endif
        enddo
        enddo

		
4. to identify which blocks are elements in
now:
nelgv = num_elements_block1 + num_elements_block2
nelgt =  num_elements_block1 + num_elements_block2 +  num_elements_block3 + num_elements_block4

but we need to distinguish all four blocks
eg : global element order

eg <= num_elements_block1 : in block 1
num_elements_block1 < eg <= (num_elements_block1 + num_elements_block2 ) : block 2
etc....

now, we can set different properties for different blocks.

