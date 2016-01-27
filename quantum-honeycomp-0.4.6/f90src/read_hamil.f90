! get onsite matrix, onsiteite energies and hoppingping from a file 
subroutine read_hamil 
use inputs
use system_variables
use sintax_hamiltonian
implicit none 
integer :: ats2,el,i1 
real (kind=8) :: holdr,holdi 
real (kind=8), allocatable :: tmp_pos(:) ! temporal position
integer :: i,j,ii,jj,kk,ihop 
character (len=200) :: line='  '
character (len=40) :: namehop
! list of names of the hoppings
integer :: io ! status of the read


! get the number of orbitals
open(78,file=hamiltonian_file) 
if (show_read_hamiltonian) write(*,*) "Reading ",hamiltonian_file
do while (.true.)
  read(78,*) line
  if (line=='DIMENSIONALITY_OF_THE_SYSTEM') exit
enddo
read(78,*) system_dimension  ! dimension of the hamiltonian

if (show_read_hamiltonian) then
  write(*,*) 'Dimensionality of the system',system_dimension
endif

close(78)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!! read dimension of the hamiltonian !!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

open(78,file=hamiltonian_file) 
do while (.true.)
  read(78,*) line
  if (line=='DIMENSION_OF_THE_HAMILTONIAN') exit
enddo
read(78,*) norbitals  ! dimension of the hamiltonian
close(78)


call initialize_hamiltonian  ! allocate hamiltonian
call setup_hop_names ! get the names of the hoppings


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!! read the real lattice !!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
open(78,file=hamiltonian_file) 
do while (.true.)
  if (hamiltonian%dimensionality==0) exit  ! if zero dimensional skip
  read(78,*,IOSTAT=io) line
  if (line=='LATTICE_VECTORS') then  ! if lattice vectors have been found
    hamiltonian%has_lattice = .true. ! put that information in the hamiltonian
    ! allocate the space
    allocate(hamiltonian%lattice_vectors(system_dimension,system_dimension))
    allocate(hamiltonian%reciprocal_vectors(system_dimension,system_dimension))
    do i=1,system_dimension
      read(78,*) hamiltonian%lattice_vectors(i,:) ! read each vector
      if (show_read_hamiltonian) then
        write(*,*) hamiltonian%lattice_vectors(i,:) ! read each vector
      endif
    enddo
    ! get the reciprocal vectors
    call create_reciprocal_lattice
    exit  ! exit the loop which searchs LATTICE_VECTORS
  endif
  if (io.ne.0) exit
enddo
close(78)


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! read the coordinates of the system !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
open(78,file=hamiltonian_file) 
do while (.true.)
  read(78,*,IOSTAT=io) line
  if (io.ne.0) exit ! if end is reached skip
  if (line=='POSITIONS') then  ! if lattice vectors have been found
    hamiltonian%has_positions = .true. ! put that information in the hamiltonian
    ! allocate the space
    if (allocated(hamiltonian%positions))  deallocate(hamiltonian%positions)
    allocate(hamiltonian%positions(hamiltonian%norbitals,3))
    if (allocated(tmp_pos)) deallocate(tmp_pos)
    allocate(tmp_pos(3)) ! allocate temporal position
    ! read a position
    do ii=1,hamiltonian%norbitals  ! loop over the number of orbitals
      read(78,*,IOSTAT=io) (tmp_pos(j),j=1,3)
      hamiltonian%positions(ii,:) = tmp_pos(:) ! save position
    enddo
    exit  ! end if the search has ended 
  endif
enddo
close(78)






!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! check if it has spin polarization !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

open(78,file=hamiltonian_file) 
do while (.true.)
  read(78,*,IOSTAT=io) line
  if (line=='WITHOUT_SPIN_POLARIZATION') then
    hamiltonian%has_spinpol=.false.  ! put hamiltonian variable to false
    if (show_read_hamiltonian)   write(*,*) 'Without spinpol on input'
    exit
  endif
  if (io.lt.0) exit  ! end of file
enddo
close(78)








!!!!!!!!!!!!!!!!!!!!!!!!!
!! get onsite elements !!
!!!!!!!!!!!!!!!!!!!!!!!!!

open(78,file=hamiltonian_file) 
do while (.true.)
  read(78,*,IOSTAT=io) line
  if (line=='ONSITE_MATRIX') exit
enddo

! if found read elements
if (io==0) then
  do while (.true.)
    read(78,*,IOSTAT=io) i,j,holdr,holdi 
    if (io.ne.0) exit
    onsite(i,j)=holdr+im*holdi 
  enddo 
endif
close(78)

! save in the class
hamiltonian%onsite = onsite


! loop over the different hoppings
do ihop=1, hamiltonian%num_hoppings
! get hopping elements 
  namehop = list_hop_names(ihop)
  if (show_read_hamiltonian) then
    write(*,*) 'Searching for  ',namehop
    write(*,*) list_hop_dirs(ihop,:)
    write(*,*) 
  endif
  open(78,file=hamiltonian_file) 
  hopping = (0.d00,0.d00)  ! initialize hopping to zero
  io = 0 ! file is ok
  do while (io==0)   ! search for the file the string
    read(78,*,IOSTAT=io) line
    if (line==namehop) exit
  enddo
  if (io==0) then ! do if the string has been found
    if (show_read_hamiltonian)  write(*,*) 'Matrix  ',namehop,'found'
    do while (.true.)
      read(78,*,IOSTAT=io) i,j,holdr,holdi 
      if (io.ne.0) exit
      hopping(i,j)=holdr+im*holdi
    enddo
    write(*,*)
  endif
  ! save name 
  hamiltonian%directional_hoppings(ihop)%hopping = hopping
  ! save direction 
  hamiltonian%directional_hoppings(ihop)%direction(:) = &
                             list_hop_dirs(ihop,:)
  hamiltonian%directional_hoppings(ihop)%label = namehop
  close(78)
enddo

! save the geometry
! call save_geometry

 
return 
end subroutine 



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! create the reciprocal lattice !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine create_reciprocal_lattice
use inputs
use system_variables
implicit none
real (kind=8) :: rvec(system_dimension,system_dimension) ! direct
complex (kind=8) :: kvec(system_dimension,system_dimension) ! reciprocal


kvec(:,:) = hamiltonian%lattice_vectors(:,:)
rvec(:,:) = hamiltonian%lattice_vectors(:,:)
call inverse(kvec,system_dimension)
hamiltonian%reciprocal_vectors(:,:) = real(kvec(:,:))




return
end subroutine create_reciprocal_lattice


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! save geometry as an independent file !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine save_geometry
use inputs
use system_variables
implicit none
integer :: i1
!!!!! write POSITIONS !!
!if (hamiltonian%has_positions) then ! if has positions
!  open(41,file="POSITIONS.OUT")
!  write(41,*) '# x    y     z   (without spin degree)'
!  do i1=1,hamiltonian%norbitals/2 ! go only for up channel
!    write(41,*) hamiltonian%positions(2*i1,:)
!  enddo
!  close(41) ! close the file
!endif ! close positions if

!!!!!!! write lattice !!
!if (hamiltonian%has_lattice) then
!  open(41,file="LATTICE.OUT")
!  do i1=1,system_dimension
!    write(41,*) hamiltonian%lattice_vectors(i1,:) ! read each vector
!  enddo
!  close(41) ! close the file
!endif  ! close lattice if



return
end subroutine


