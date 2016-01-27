! new mean field hamiltonian 
subroutine new_mf_ham
use inputs
use system_variables 
use mean_field_variables
implicit none 
integer :: info,ind,i1,i2,i3,i4 
integer :: ie 
real (kind=4) :: energies2(num_wf) 
complex (kind=8) ::  evmfa(num_AB),wfi(norbitals) 
complex (kind=8) ::  evmfb(num_AB) 
complex (kind=8) ::  cholda,choldb,holdr 
integer :: nfermi 
integer :: n_occ_states
! smearing 
real (kind=8) :: csmearing,sumcsmearing 

sumcsmearing=0.d00 

energies2=real(energies) 
! sort the eigenvalues 
call slasrt('I', num_wf, energies2, info ) 
if (info /= 0) then 
  write(*,*) 'Error sorting the eigenvalues' 
  stop
  return 
endif 


! obtain fermi energy 
  n_occ_states = dble(num_wf)*filling   ! number of occupied states


! add extra electrons!!!
! for periodix systems multiply by the number of k-points
  if (system_dimension==0) then
    n_occ_states = n_occ_states + extra_electrons 
  else
    n_occ_states = n_occ_states + extra_electrons * nkpoints
  endif

! check if too many electrons
  if (n_occ_states.gt.norbitals*nkpoints) then 
    write(*,*) 'Error, too few eigenvalues' ,n_occ_states
    stop  ! pause the full program
    return 
  endif 

  ! get the fermi energy from the filling
  fermi=dble(energies2(n_occ_states)+energies2(n_occ_states+1))/2.d00 
  ! and shift it according to input
  fermi = fermi + shift_fermi

!write(*,*) 'Fermi energy = ',fermi

! calculate new mean field hamiltonian 
! ==================================== 
mf_ham=(0.d00,0.d00) 
evmfa=(0.d00,0.d00) 
evmfb=(0.d00,0.d00) 

! calculate expectation value of the operators 

if (show_debug) then
  write(*,*) num_AB,'numab'
endif

!$omp parallel default (shared) private (i1,cholda,choldb,ie,wfi,csmearing,sumcsmearing,i2,i3,i4)
!$omp do
do i1=1,num_AB 
  cholda=(0.d00,0.d00) 
  choldb=(0.d00,0.d00) 
  do ie=1,num_wf
    wfi = wf(ie,:) 
    if (smearing_type_scf =="linear" ) then ! if the semaring is of linear type
      csmearing=0.d00 
      ! calculate coefficient 
      if (energies(ie) < (fermi-smearing)) then 
      csmearing=1.00 
      else if (energies(ie) > (fermi+smearing)) then 
      csmearing=0.00 
      else 
      csmearing=1.00-(energies(ie)-fermi+smearing)/(2.d00*smearing) 
      if (dabs(csmearing).gt.11.d-01) then
         write(*,*) "Error in linear smaering"
         stop
      endif
      endif 
      sumcsmearing=sumcsmearing+csmearing 
    ! if the semaring is of fermi-dirac
    else if (smearing_type_scf =="fermi-dirac" ) then
      call fermi_dirac(energies(ie),smearing,fermi,csmearing)  ! fermi dirac
    else
      write(*,*) "Unrecognised"
      stop
    endif  ! close smearing linear 
      ! matrix A 
      do i2=1,ab_matrices(i1)%a%nv   ! number of non vanishing elements
        i3=ab_matrices(i1)%a%i(i2) 
        i4=ab_matrices(i1)%a%j(i2) 
        cholda=cholda+conjg(wfi(i3))*ab_matrices(i1)%a%mij(i2) & 
        *wfi(i4)*csmearing 
      enddo 
      ! matrix B 
      do i2=1,ab_matrices(i1)%b%nv   ! number of non vanishing elements
        i3=ab_matrices(i1)%b%i(i2) 
        i4=ab_matrices(i1)%b%j(i2) 
        choldb=choldb+conjg(wfi(i3))*ab_matrices(i1)%b%mij(i2) & 
        *wfi(i4)*csmearing 
      enddo 
  enddo 
  ! store in the mean field expectation values array 
  evmfa(i1)=cholda
  evmfb(i1)=choldb 
  
enddo 
!$omp end do
!$omp end parallel

!!! if there are kpoints normalize !!!!
evmfa = evmfa/dble(klist%nkpoints)
evmfb = evmfb/dble(klist%nkpoints)
  
  ! now calculate the mean field hamiltonian 
do i1=1,num_AB 
! <A>B contribution 
  do i2=1,ab_matrices(i1)%b%nv 
    i3=ab_matrices(i1)%b%i(i2) 
    i4=ab_matrices(i1)%b%j(i2) 
    mf_ham(i3,i4)=mf_ham(i3,i4)+evmfa(i1)*ab_matrices(i1)%b%mij(i2) &
      *ab_matrices(i1)%lambda 
  enddo 
  ! <B>A contribution 
  do i2=1,ab_matrices(i1)%a%nv
    i3=ab_matrices(i1)%a%i(i2) 
    i4=ab_matrices(i1)%a%j(i2) 
    mf_ham(i3,i4)=mf_ham(i3,i4)+evmfb(i1)*ab_matrices(i1)%a%mij(i2) &
      *ab_matrices(i1)%lambda 
  enddo 
enddo 

! calculate <A><B> contribution (only used when calculating energy) 

holdr=0.d00 
do i1=1,num_AB
  holdr=holdr+dble(evmfa(i1)*evmfb(i1)*ab_matrices(i1)%lambda) 
enddo 




! guarantee that it is hermitic 
mf_ham=(mf_ham+transpose(conjg(mf_ham)))/2.d00 

! energy of the system 
ene_tot=0.d00 
! input is filling 
!if (filloref == 0) then 
  do i1=1,n_occ_states 
    ene_tot=ene_tot+dble(energies2(i1)) 
  enddo 


! normalize for the kpoints
  if (system_dimension.gt.0) then
    ene_tot = ene_tot/dble(nkpoints)
  endif


! substract DC term 
ene_tot=ene_tot-dble(holdr) 


! write energy in file
open(19,file="ENERGY.OUT",ACCESS = 'APPEND')
write(19,*) ene_tot,ene_tot/dble(norbitals)
write(*,*) 'Total energy = ',ene_tot
close(19)


return 
end subroutine 
