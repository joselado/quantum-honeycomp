! Calculates the density of states of the system
subroutine density_of_states
use inputs
use system_variables
use dos_variables
use expectation_values_variables
implicit none


write(*,*) "============================="
write(*,*) "Calculating density of states"
write(*,*) "============================="
write(*,*)

call read_inputs  ! read the input file
nkpoints = nkpoints*refine_kmesh_dos  ! refine the mesh for DOS calculation
call read_hamil   ! obtain hamiltonian



! setup the energy window
emin = emin_dos
emax = emax_dos 
use_ewindow = use_ewindow_dos ! use energy window


call initialize_klist  ! obtain klist


! setup the matrices for expectation value if needed
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

num_operators = 0  ! assuming no operators

if (dos_operators_option=='from_file') then
  num_operators = 1  ! just one operator
  call read_expectation_operator  ! read the operator
  write(*,*) 'Reading operator from file'
else if (dos_operators_option=='default') then
  write(*,*) 'Using default operators'
  call default_operators  ! read the operator
else if (dos_operators_option=='None') then
  num_operators = 0
endif


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!



! calculate eigenvalues
if (num_operators==0) then
  call initialize_eigenstates  ! initialize eigenstates
  call all_eigenvalues ! get all eigenvalues
endif

! calculate eigenvalues
if (num_operators.gt.0) then
  call initialize_eigenstates   ! calculate eigenvectors as well
  call all_eigen   ! calculate eigenvectors
endif

call calculate_dos ! calculate density of states



! calculate LDOS as well


if (system_dimension.gt.0) then  ! if periodic normalize by the kmesh
  value_dos = value_dos/dble(nkpoints) ! normalize by the number of kpoints
endif


call write_dos ! write density of states

call clean_dos ! clean DOS variables
call clean_waves

return
end subroutine

