! Calculates the dielectric response function
subroutine dielectric_response
use inputs
use system_variables
use dielectric_variables
use dielectric_operators
implicit none
integer :: ik,ii


write(*,*) "============================================"
write(*,*) "Calculating the dielectric response function"
write(*,*) "============================================"
write(*,*)

call read_inputs  ! read the input file
nkpoints = nkpoints*refine_kmesh_dielectric  ! refine mesh for dielectric
call read_hamil   ! obtain hamiltonian
call initialize_eigenvectors_k  ! initialize eigenvectors




do_rpa = .true.  ! activate RPA routine

! energy window options
use_ewindow = .true.
use_ewindow = .false.
emin = emin_chi - ewindow_chi 
emax = emax_chi + ewindow_chi 


! generate all the ab matrices
if (dielectric_type=='global_all') then
  call setup_all_ab_matrices
else if (dielectric_type=='local_density') then
  call setup_local_density
!else if (dielectric_type=='global_density_rpa') then
!  call setup_ab_global_density
else if (dielectric_type=='local_xy_spin') then
  call setup_local_xy_spin
else if (dielectric_type=='local_zz_spin') then
  call setup_local_zz_spin
else if (dielectric_type=='from_file') then
  call setup_response_operators_from_file
else
  write(*,*) 'Unrecognised option in dielectric_type, stopping'
  stop
endif




!!!! check dimension of the system and initialize klist if neccesary !!!!
call initialize_klist  




call initialize_chi_ab  ! initialize chi_tensor_w variables

! initialize the qvector
! for 0d, no vetor is needed
if (hamiltonian%dimensionality.eq.1)  qvector(1) = q_chi  ! for 1d
if (hamiltonian%dimensionality.eq.2)  then
   qvector(1) = q_chi*cos(3.141592*2.d00*phi_chi)  ! for 2d
   qvector(2) = q_chi*sin(3.141592*2.d00*phi_chi)  ! for 2d
endif





!!! 0d systems !!!
if (system_dimension==0) then  ! open 0d
  write(*,*) "*** Performing Chi calculation for 0d system"
  write(*,*)
  call initialize_eigenstates ! initialize variables
  call all_eigen  ! diagonalize
  ene_chi_i = energies     ! first energies
  ene_chi_j = energies     ! first energies
  wf_chi_i = wf            ! first waves
  wf_chi_j = wf            ! first waves
  !!! retain only the energies which lie in an interval !!!!
  call calculate_chi_ab ! calculate chi_tensor_w
  
endif ! close 0d

chi_total = (0.d00,0.d00)



!!!!! 1d systems !!!!
if (system_dimension>0) then  ! open 1d
  write(*,*) "*** Q-vector for dielectric is",qvector
  write(*,*)
  do ik=1,klist%nkpoints
    write (*,*) "Completed",ik,"of",nkpoints
    call eigenvectors_kpoint(klist%kpoint(ik,:))   ! eigenvalues in k
    ene_chi_i = energies_k     ! first energies
    wf_chi_i = wf_k            ! first waves
    num_wf_i = num_wf_k
    call eigenvectors_kpoint(klist%kpoint(ik,:)+qvector(:))   ! eigenvalues in k
    ene_chi_j = energies_k      ! second energies
    wf_chi_j = wf_k             ! second waves
    num_wf_j = num_wf_k
    !!! retain only the energies which lie in an interval !!!!
    call calculate_chi_ab ! calculate chi_tensor_w
  enddo
  chi_total = chi_total/dble(klist%nkpoints)  ! normalize for the kpoints
endif ! close 1d



if (do_rpa) then ! do RPA
  write(*,*) '==============================='
  write(*,*) 'Calculating RPA linear response'
  write(*,*) '==============================='
  call rpa_response(int(sqrt(real(num_chi_ab))))
endif


! write trace of the response function
!call write_trace_chi_ab


! write the full response
!if (write_full_chi) then
  call write_full_chi_ab ! write chi_w
!endif


call clean_chi_ab ! clean chi_w variables

return
end subroutine

