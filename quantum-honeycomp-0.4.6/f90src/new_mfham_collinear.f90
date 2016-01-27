! new mean field hamiltonian 
subroutine new_mf_ham_collinear
use inputs
use system_variables 
use mean_field_variables
implicit none 
integer :: info,ind,i1,i2,i3,i4 
integer :: ie 
real (kind=4) :: energies2(num_wf) 
complex (kind=8) ::  wfi(norbitals) 
real (kind=8) ::  deni(norbitals) 
real (kind=8) ::  enedc ! double counting energy
real (kind=8) ::  enedc1(norbitals/2) ! double counting energy
real (kind=8) ::  enedc2(norbitals/2) ! double counting energy
integer :: norbspinless
integer :: nfermi 
integer :: n_occ_states
! smearing 
real (kind=8) :: csmearing,sumcsmearing 

sumcsmearing=0.d00 
norbspinless = norbitals/2

deni =0.d00
enedc1 = 0.d00
enedc2 = 0.d00

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

! calculate expectation value of the operators 

if (show_debug) then
  write(*,*) num_AB,'numab'
endif

do ie=1,num_wf
  wfi = wf(ie,:) 
  if (smearing_type_scf =="linear" ) then ! if the semaring is of linear type
    csmearing=0.d00 
    ! calculate coefficient 
    if (energies(ie) < (fermi-smearing)) then 
    csmearing=1.00 
    else if (energies(ie) > (fermi+smearing)) then 
    csmearing=0.00 
    cycle  ! next iteration
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

  deni = deni + wfi*conjg(wfi)*csmearing ! density of this wave
enddo  ! close loop over wf
  ! store in the mean field expectation values array 

do i1 = 1,norbspinless
  mf_ham(2*i1,2*i1) = mf_ham(2*i1,2*i1) + deni(2*i1-1) ! down density
  mf_ham(2*i1-1,2*i1-1) = mf_ham(2*i1-1,2*i1-1) +  deni(2*i1) ! up density
  enedc1(i1) = enedc1(i1) + deni(2*i1-1) ! double counting
  enedc2(i1) = enedc2(i1) + deni(2*i1) ! double counting
enddo ! end loop over sites




  

!!! if there are kpoints normalize !!!!
mf_ham = mf_ham/dble(klist%nkpoints)*hubbard_scf
enedc1 = enedc1/dble(klist%nkpoints)
enedc2 = enedc2/dble(klist%nkpoints)

enedc = sum(enedc1*enedc2)*hubbard_scf ! double counting correction
write(*,*) enedc  
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
ene_tot=ene_tot-enedc


! write energy in file
open(19,file="ENERGY.OUT",ACCESS = 'APPEND')
write(19,*) ene_tot,ene_tot/dble(norbitals)
write(*,*) 'Total energy = ',ene_tot
close(19)


return 
end subroutine 
