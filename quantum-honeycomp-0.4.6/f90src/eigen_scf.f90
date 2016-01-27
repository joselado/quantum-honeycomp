! calculate all the eigenfunctions in all the k-points
subroutine eigen_scf
use inputs
use system_variables
use offdiagonal ! for ewindow diagonalization
implicit none 
complex (kind=8) :: hk(norbitals,norbitals)  ! k-dependent hamiltonian
complex (kind=8) :: e(norbitals,norbitals)  ! full onsite
complex (kind=8) :: tk(norbitals,norbitals)  ! k-dependent hopping
complex (kind=8) :: wfk(norbitals,norbitals)  ! wavefunctions for a k
complex (kind=8) :: phik  ! complex bloch phase
real (kind=8) :: enek(norbitals)  ! eigenvalues for a k
real (kind=8) :: kpoint(hamiltonian%dimensionality)  ! kpoint
real (kind=8) :: ndir(hamiltonian%dimensionality)  ! kpoint
integer :: ii,jj,kk,ik,ihop
integer :: num_offdiagonals,tempoff ! for the number of offdiagonals






! print type of diagonalization
if (show_diagonalization) then
  if (.not.use_ewindow) write(*,*) "***Calculating all eigenvectors***"
  if (use_ewindow) write(*,*) "***Calculating eigenvectors in window***", &
                 emin,emax
endif





! if (.not.use_ewindow) then   ! full diagoalization
  ! loop over kvectors 
  kk = 1
  num_wf = 0 ! initialize to zero the total wavefunctions
  onsite = hamiltonian%onsite ! store onsite part
  if ((hamiltonian%dimensionality.eq.0)) then
    call diagonalize(norbitals,onsite,energies,wf)   ! diagonalize
    klist%nev_found(1) = norbitals  ! saves this number found
  else if (hamiltonian%dimensionality.gt.0) then
!$omp parallel default (shared) private (kpoint,hk,hopping,ndir,phik,tk,enek,wfk,ik,ihop,ii,kk)
!$omp do
  do ik=1,klist%nkpoints
    kpoint(:)=klist%kpoint(ik,:) ! select the kpoint
    hk = onsite ! add onsite
    do ihop = 1, hamiltonian%num_hoppings
        hopping(:,:) = hamiltonian%directional_hoppings(ihop)%hopping(:,:)
        ! calculate the complex phase
        ndir(:) = hamiltonian%directional_hoppings(ihop)%direction(:)
        phik = dot_product(ndir,kpoint)
        phik = im*pi*2.d00*phik
        tk = hopping*exp(phik) ! calculate bloch hopping
        hk = hk + tk   ! add to the hamiltonian
    enddo ! end loop over hoppings
      !! full diagonalization
    if (is_collinear) then
      call diagonalize_collinear(norbitals,hk,enek,wfk)   ! diagonalize
    else
      call diagonalize(norbitals,hk,enek,wfk)   ! diagonalize
    endif
    klist%nev_found(ik) = norbitals  ! saves this number found
    ! store eigenfunctions and eigenvalues
    do ii=1,hamiltonian%norbitals
      kk = norbitals*(ik-1) + ii ! index of the state
      wf(kk,:) = wfk(ii,:)  ! saving wavefunctions
      energies(kk) = enek(ii) ! saving eigenvalues
    enddo 
  enddo  ! close k-point loop
!$omp end do
!$omp end parallel
  else
   write(*,*) "ERROR in all_eigen.f90"
   stop
  
  endif ! close parallel if
  num_wf = klist%nkpoints*hamiltonian%norbitals
! endif  ! close full diagonalization


return
end subroutine 
