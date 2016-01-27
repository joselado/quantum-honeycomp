! calculate all the eigenfunctions in all the k-points
subroutine eigenvectors_kpoint(kpoint) 
use inputs
use system_variables
use offdiagonal
implicit none
complex (kind=8) :: tk(hamiltonian%norbitals,hamiltonian%norbitals) 
complex (kind=8) :: hk(hamiltonian%norbitals,hamiltonian%norbitals) 
complex (kind=8) :: phik  ! complex bloch phase
real (kind=8) :: kpoint(hamiltonian%dimensionality)  ! kpoint
integer :: ii,jj,kk,ik,ihop
integer :: num_offdiagonals

! print type of diagonalization
if (show_diagonalization) then
  if (.not.use_ewindow) write(*,*) "***Calculating all eigenvectors***"
  if (use_ewindow) write(*,*) "***Calculating eigenvectors in window***", &
                 emin,emax
endif

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
  call diagonalize(norbitals,hk,energies_k,wf_k)   ! diagonalize
  num_wf_k = norbitals
  ! store eigenfunctions and eigenvalues
endif

! energy window diagonalization
if (use_ewindow) then
  call get_number_offdiagonals(hk,hamiltonian%norbitals, &
            num_offdiagonals)
  call eigenvectors_ewindow(norbitals,hk,nev_found,energies_k, &
         wf_k,emin,emax, num_offdiagonals)
  num_wf_k = nev_found  ! add this number of waves
endif



return
end subroutine eigenvectors_kpoint
