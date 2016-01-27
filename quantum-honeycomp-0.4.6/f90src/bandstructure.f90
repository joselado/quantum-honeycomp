subroutine bandstructure
use system_variables  ! variables from this particular system 
use bands_variables ! routines to write outputs
use expectation_values_variables ! for the calculation of expectation values

implicit none 
integer :: ik ! counter for waves


write(*,*) "========================="
write(*,*) "Calculating bandstructure"
write(*,*) "========================="
write(*,*)



! setup variables if ewindow
use_ewindow = use_ewindow_bands
emin = emin_bands
emax = emax_bands




call read_inputs ! read the input file to reinitialice the parameters
nkpoints = nkpoints*refine_kmesh_bands  ! refine the mesh for bands calculation
call read_hamil   ! obtain hamiltonian


!!!! klist for diferent dimensions !!!!!
! for 0 is not needed, nothing to do
if ((klist_bands=='default').or.(hamiltonian%dimensionality==0)) then
  call initialize_klist  ! obtain klist
!  write(*,*) "Generated klist"
else if (klist_bands=='from_file') then
  call read_klist
else
  write(*,*) 'Invalid option in klist_bands',klist_bands
  stop
endif
! calculate eigenvectors if neccesary, otherwise only eigenvalues


num_operators = 0  ! assuming no operators

if (bands_operators_option=='from_file') then
  num_operators = 1  ! just one operator
  call read_expectation_operator  ! read the operator
else if (bands_operators_option=='default') then
  call default_operators  ! read the operator
else if (bands_operators_option=='Sx') then
  call operator_sx  ! read the operator
else if (bands_operators_option=='Sy') then
  call operator_sy  ! read the operator
else if (bands_operators_option=='Sz') then
  call operator_sz  ! read the operator
else if (bands_operators_option=='index') then
  call operator_index  ! read the operator
else if (bands_operators_option=='None') then
  num_operators=0  ! no operators
else
  write(*,*) 'Unrecognised option in bands_operator option'
  write(*,*) "Valid options are"
  write(*,*) "   - from_file"
  write(*,*) "   - default"
  write(*,*) "   - Sx"
  write(*,*) "   - Sy"
  write(*,*) "   - Sz"
  write(*,*) "   - Index"
  write(*,*) "   - None"
  write(*,*) "STOPING CALCULATION"
  stop
endif



if (.not.economic_memory_bands) then  
  write(*,*) "Memory demanding bands"
  ! for only band structure
  if (num_operators.eq.0) then 
    call initialize_eigenvalues ! calculate only eigenvalues
    call all_eigenvalues ! calculate only eigenvalues
  endif
  ! for band structure and expectation value 
  if (num_operators.gt.0) then
    call initialize_eigenstates   ! calculate eigenvectors as well
    call all_eigen   ! calculate eigenvectors
    call calculate_exp_val_values ! calculate the values
  endif
  call write_bands ! write bandstructure in BANDS.OUT
endif


if (economic_memory_bands) then ! diagonalice point by point
  write(*,*) "Memory economic bands"
!  write(*,*) "BUGGED"
!  stop
  ! clean file
  call write_bands_kpoint(1,.true.) 
  
  ! initialice vectors
  if (num_operators.gt.0) then
    call initialize_eigenvectors_k ! initialice keigenvector
  endif
  if (num_operators.eq.0) then
    call initialize_eigenvalues_k ! initialice keigenvector
  endif

    do ik=1,klist%nkpoints ! go point by point
      write(*,*) "Done",ik,"from",klist%nkpoints
      ! eigenvalues and eigenvectors
      if (num_operators.gt.0) then
        call eigenvectors_kpoint(klist%kpoint(ik,:)) ! get eigenvectors
        num_wf = num_wf_k ! number of wavefunctions
        call calculate_eev_kpoint(wf_k,energies_k,num_wf, &
              hamiltonian%norbitals) ! calculate EV for kp
      endif
      ! only eigenvalues
      if (num_operators.eq.0) then
        call eigenvalues_kpoint(klist%kpoint(ik,:)) ! get eigenvectors
        num_wf = num_wf_k ! number of wavefunctions
      endif
      klist%nev_found(ik)  = num_wf ! number of eigenvalues found
      call write_bands_kpoint(ik,.false.) 
! write bandstructure in BANDS.OUT
    enddo
endif









! write klist in KLIST.OUT
if (hamiltonian%dimensionality.gt.0)  call write_klist_bands 

write(*,*) "Written bandstructure in BANDS.OUT"


return 
end subroutine 
