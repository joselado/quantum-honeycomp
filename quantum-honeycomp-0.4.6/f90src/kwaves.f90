! calculate all the eigenfunctions in one kpoint
subroutine kwaves(kpoint,priv_waves,priv_energies) 
use inputs
use system_variables
implicit none
real (kind=8), intent(in) :: kpoint(hamiltonian%dimensionality)  ! kpoint
complex (kind=8), intent(out) :: priv_energies(hamiltonian%norbitals) 
complex (kind=8), intent(out) :: priv_waves(hamiltonian%norbitals,hamiltonian%norbitals) 

complex (kind=8) :: tk(hamiltonian%norbitals,hamiltonian%norbitals) 
complex (kind=8) :: hop(hamiltonian%norbitals,hamiltonian%norbitals) 
complex (kind=8) :: hk(hamiltonian%norbitals,hamiltonian%norbitals) 
! wavefunctions
! eigenenergies
complex (kind=8) :: phik  ! complex bloch phase
integer :: ii,jj,kk,ik,ihop



hk = hamiltonian%onsite ! store onsite part
do ihop = 1, hamiltonian%num_hoppings
    hop = hamiltonian%directional_hoppings(ihop)%hopping
    ! calculate the complex phase
    phik = sum(hamiltonian%directional_hoppings(ihop)%direction*kpoint)
    phik = im*pi*2.d00*phik
    tk = hop*exp(phik) ! calculate bloch hopping
    hk = hk + tk   ! add to the hamiltonian
enddo ! end loop over hoppings

! diagonalize
call diagonalize(hamiltonian%norbitals,hk,priv_energies,priv_waves)   



return
end subroutine kwaves
