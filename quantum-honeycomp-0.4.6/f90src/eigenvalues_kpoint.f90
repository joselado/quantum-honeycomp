! calculate all the eigenfunctions in all the k-points
subroutine eigenvalues_kpoint(kpoint) 
use inputs
use system_variables
use offdiagonal
implicit none
complex (kind=8) :: tk(hamiltonian%norbitals,hamiltonian%norbitals) 
complex (kind=8) :: hk(hamiltonian%norbitals,hamiltonian%norbitals) 
real (kind=8) :: evals(hamiltonian%norbitals) 
complex (kind=8) :: phik  ! complex bloch phase
real (kind=8) :: kpoint(hamiltonian%dimensionality)  ! kpoint
integer :: ii,jj,kk,ik,ihop,num_offdiagonals


! print type of diagonalization
!if (.not.use_ewindow) write(*,*) "***Calculating all eigenvectors***"
!if (use_ewindow) write(*,*) "***Calculating eigenvectors in window***", &
!                 emin,emax


! loop over kvectors 
kk = 1
num_wf = 0 ! initialize to zero the total wavefunctions
onsite = hamiltonian%onsite ! store onsite part
hk = onsite ! add onsite
do ihop = 1, hamiltonian%num_hoppings
    hopping = hamiltonian%directional_hoppings(ihop)%hopping
    ! calculate the complex phase
    phik = sum(hamiltonian%directional_hoppings(ihop)%direction*kpoint)
    phik = im*pi*2.d00*phik
    tk = hopping*exp(phik) ! calculate bloch hopping
    hk = hk + tk   ! add to the hamiltonian
enddo ! end loop over hoppings


  
  !! full diagonalization
if (.not.use_ewindow) then
  write (*,*) norbitals,hamiltonian%norbitals
  call eigenvalues(hamiltonian%norbitals,hk,energies_k)   ! diagonalize
  num_wf_k = norbitals
  nev_found = norbitals
  ! store eigenfunctions and eigenvalues
endif

! energy window diagonalization
if (use_ewindow) then
  write (*,*) norbitals,hamiltonian%norbitals
  call get_number_offdiagonals(hk,hamiltonian%norbitals, &
            num_offdiagonals)
  call eigenvalues_ewindow(hamiltonian%norbitals,hk,nev_found, &
                            energies_k,emin,emax,num_offdiagonals)
  write(*,*) "found",nev_found
  num_wf_k = nev_found  ! add this number of waves
endif



return
end subroutine eigenvalues_kpoint

!
!subroutine offdiagonals(hamiltonian)::
!  use system_variables, only hamiltonian_type
!  implicit none
!  type(hamiltonian_type) :: hamiltonian ! hamiltonian
!  integer num_off
!
!  num_offdiagonals = 0 ! initial number
!  call get_number_offdiagonals(hamiltonian%onsite,hamiltonian%norbitals, &
!            num_offdiagonals)  ! first the onsite
!  ! now the hoppings
!  do ihop = 1, hamiltonian%num_hoppings
!      hopping(:,:) = hamiltonian%directional_hoppings(ihop)%hopping(:,:)
!      call get_number_offdiagonals(hopping,hamiltonian%norbitals, &
!            tempoff)  ! first the onsite
!  ! if bigger number, update
!  if (tempoff.gt.num_offdiagonals) num_offdiagonals = tempoff
!  write(*,*) "Offdiagonals",num_offdiagonals,tempoff
!  enddo
!end subroutine offdiagonals
!
!
!
