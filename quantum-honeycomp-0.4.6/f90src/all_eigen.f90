! calculate all the eigenfunctions in all the k-points
subroutine all_eigen 
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





! if use ewindow, calculate the offdigonals
if (use_ewindow) then
  num_offdiagonals = 0 ! initial number
  call get_number_offdiagonals(hamiltonian%onsite,hamiltonian%norbitals, &
            num_offdiagonals)  ! first the onsite
  ! now the hoppings
  do ihop = 1, hamiltonian%num_hoppings
      hopping(:,:) = hamiltonian%directional_hoppings(ihop)%hopping(:,:)
      call get_number_offdiagonals(hopping,hamiltonian%norbitals, &
            tempoff)  ! first the onsite
  ! if bigger number, update
  if (tempoff.gt.num_offdiagonals) num_offdiagonals = tempoff
  write(*,*) "Offdiagonals",num_offdiagonals,tempoff
  enddo

  ! Now diagonlize
  ! loop over kvectors 
  kk = 1
  num_wf = 0 ! initialize to zero the total wavefunctions
  onsite = hamiltonian%onsite ! store onsite part
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
    call eigenvectors_ewindow(norbitals,hk,nev_found,enek,wfk, &
         emin,emax,num_offdiagonals)
    klist%nev_found(ik) = nev_found  ! saves this number found
    num_wf = num_wf + nev_found  ! add this number of waves
    ! store eigenfunctions and eigenvalues
    do ii=1,nev_found
      wf(kk,:) = wfk(ii,:)  ! saving wavefunctions
      energies(kk) = enek(ii) ! saving eigenvalues
      kk = kk + 1 ! update the index
    enddo
  enddo  ! close k-point loop

endif




if (.not.use_ewindow) then   ! full diagoalization
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
    call diagonalize(norbitals,hk,enek,wfk)   ! diagonalize
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
endif  ! close full diagonalization


return
end subroutine 
