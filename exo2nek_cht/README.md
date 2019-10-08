This is a special version of exo2nek to for conjugate heat transfer mesh (cht)
All elements must be hex20.

Workflow:
assuming CUBIT or Trelis (The commercial version of CUBIT) is not available.
1. Mesh generation in ANSYS meshing or ICEM
2. Convert all mesh to hex20 in ICEM. Having fluid mesh and solid mesh in different parts.
3. In ICEM, exporting fluid mesh to one exo file. This can be done be removing mesh in other parts.
   Remember to set sideset before exporting to exo file.
4. Exporting solid mesh in the same way.
5. Running exo2nek, input fluid exo files first, then solid exo files. 
   Currently only accept fluid-solid interface. 
   That is to say, one integrated fluid domain should not be exported to different exo files.
   This is also true for solid domains.

Usage in Nek:

1. You need to identify which exo files are elements in
now, assuming exo files that are feed in are named exo1,exo2,exo3,etc
assume exo1 is fluid, exo2,exo3 are solid

now:
nelgv = num_elements_exo1
nelgt = num_elements_exo1+num_elements_exo2+num_elements_exo3

eg : global element order

elements in exo1: eg <= num_elements_exo1 : in exo1
elements in exo2: num_elements_exo1 < eg <= (num_elements_exo1 + num_elements_exo2 ) 
elements in exo3: (num_elements_exo1 + num_elements_exo2 )  < eg 

now, we can set different properties for different blocks.

it is suggested to make a subroutine in your usr file to to this conversion

c----------------------------------------------
      subroutine get_exo_id(eglobal)
      parameter (n_exos=3)
      integer num_elements_exo(n_exos)
      integer exo_id,eglobal
      integer exo_ground,exo_ceiling
      
      num_elements_exo(1) = ? 
      num_elements_exo(2) = ? 
      num_elements_exo(3) = ? 
      
      do i_exo =1,n_exos
      
      if(i_exo.eq.1) then
      exo_ground = 0
      exo_ceiling = num_elements_exo(1)
         if((eglobal.gt.exo_ground).and.(eglobal.le.exo_ceiling)) then
           exo_id = i_exo
           exit
         endif 
      else 
      exo_ground = exo_ground + num_elements_exo(i_exo-1)
      exo_ceiling = exo_ceiling + num_elements_exo(i_exo)
         if((eglobal.gt.exo_ground).and.(eglobal.le.exo_ceiling)) then
           exo_id = i_exo
           exit
         endif 
      endif

      enddo
      
      return exo_id
      end
c----------------------------------------------

2. setup bc
assuming you are using get_exo_id to obtain exo id for a given element

to set actual boundary conditions
in usrdat2()

       do iel=1,nelt                       ! must use nelt
        eglobal = lglel(iel) ! convert local element number to global element number
        iexo = get_exo_id(eglobal)
        if(iexo.eq.1) then                 ! set up boundary condition for elements in exo1

          do ifc=1,2*ndim 
          ! set boundary for velocity
          id_face1 = bc(5,ifc,iel,1)       
          if (id_face1.eq.2) then 
             cbc(ifc,iel,1) = 'v  '
          elseif (id_face1.eq.3) then
             cbc(ifc,iel,1) = 'O  '
          elseif (id_face1.eq.4) then
             cbc(ifc,iel,1) = 'W  '
	      elseif (id_face1.eq.5) then      ! 
             cbc(ifc,iel,1) = 'W  '
          endif
          
          ! set boundary for thermal field, should skip fluid-solid interface
          id_face2 = bc(5,ifc,iel,2)
          if (id_face2.eq.2) then 
             cbc(ifc,iel,2) = 't  '
          elseif (id_face2.eq.3) then
             cbc(ifc,iel,2) = 'O  '
	      elseif (id_face2.eq.4) then
             cbc(ifc,iel,2) = 'I  '
          endif
          enddo

        else if (iexo.eq.2) then 
         ! follow the set up for iexo =1 
        else if (iexo.eq.3) then
         ! follow the set up for iexo =1 
        endif
 
        enddo
