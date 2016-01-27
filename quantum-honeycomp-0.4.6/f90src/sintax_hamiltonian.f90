
!!!! module to manage the sintax of the file hamiltonian.in !!!
module sintax_hamiltonian
use inputs
use system_variables
implicit none

! list of names of hte hoppings
character (len=40), dimension(:), allocatable :: list_hop_names
real (kind=8), allocatable :: list_hop_dirs(:,:)

contains
  subroutine setup_hop_names
  implicit none
  
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! Create the list of names of the hoppings !
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  if (allocated(list_hop_names))  deallocate(list_hop_names)   ! deallocate
  if (allocated(list_hop_dirs))  deallocate(list_hop_dirs)   ! deallocate
  allocate(list_hop_names(hamiltonian%num_hoppings))  ! and allocate
  allocate(list_hop_dirs(hamiltonian%num_hoppings,hamiltonian%dimensionality))  ! and allocate

  list_hop_dirs = 0.d00  
  ! each dimension has its own names
  if (hamiltonian%dimensionality==1) then
    list_hop_names(1) = 'HOPPING_MATRIX_1'
    list_hop_dirs(1,:) = (1.d00)
    list_hop_names(2) = 'HOPPING_MATRIX_-1'
    list_hop_dirs(2,:) = (-1.d00)
  endif
  
  if (hamiltonian%dimensionality==2) then
    list_hop_names(1) = 'HOPPING_MATRIX_1'
    list_hop_dirs(1,1) = 1.d00
    list_hop_names(2) = 'HOPPING_MATRIX_-1'
    list_hop_dirs(2,1) = -1.d00
    list_hop_names(3) = 'HOPPING_MATRIX_0_1'
    list_hop_dirs(3,2) = 1.d00
    list_hop_names(4) = 'HOPPING_MATRIX_0_-1'
    list_hop_dirs(4,2) = -1.d00
    list_hop_names(5) = 'HOPPING_MATRIX_1_1'
    list_hop_dirs(5,1) = 1.d00
    list_hop_dirs(5,2) = 1.d00
    list_hop_names(6) = 'HOPPING_MATRIX_-1_1'
    list_hop_dirs(6,1) = -1.d00
    list_hop_dirs(6,2) = 1.d00
    list_hop_names(7) = 'HOPPING_MATRIX_1_-1'
    list_hop_dirs(7,1) = 1.d00
    list_hop_dirs(7,2) = -1.d00
    list_hop_names(8) = 'HOPPING_MATRIX_-1_-1'
    list_hop_dirs(8,1) = -1.d00
    list_hop_dirs(8,2) = -1.d00
  endif
  
  
  if (hamiltonian%dimensionality==3) then
    write(*,*) '3d not implemented in read_hamil.f90'
    stop
  endif
  
  
  end subroutine setup_hop_names

end module sintax_hamiltonian
