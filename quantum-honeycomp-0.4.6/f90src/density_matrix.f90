subroutine ground_density_matrix
use inputs  ! input variables
use system_variables  ! variables from this particular system 
use density_matrix_variables  ! variables and subroutines for denmat


implicit none 

character (len=8) :: aa

write(*,*) "======================================="
write(*,*) "Calculating ground state density matrix"
write(*,*) "======================================="
write(*,*)


use_ewindow = .false.   ! calculate all the eigenvectors
!emin = emin_bands
!emax = emax_bands




!!!! klist for diferent dimensions !!!!!
!nkpoints = nkpoints*refine_kmesh_denmat  ! refine the mesh for bands calculation
call read_inputs ! read the input file to reinitialice the parameters
hamiltonian_file = "hamiltonian.in"
call read_hamil   ! obtain hamiltonian
call initialize_klist ! This has to be initialized after reading hamiltonian
call initialize_eigenstates   ! calculate eigenvectors 

write(*,*) "nkpoints",klist%nkpoints


call initialize_density_matrix  ! allocate space for the density matrix
write(*,*) "nkpoints",klist%nkpoints
! calculate eigenvectors 
call all_eigen   ! calculate all eigenvectors
call calculate_density_matrix  ! calculate the density matrix

if (system_dimension.gt.0) then ! if not dot, normalize
  den_mat = den_mat/dble(nkpoints)
  magnetism = magnetism / dble(nkpoints)
  density = density / dble(nkpoints)
endif



! write the different outputs!!!!!
if (write_denmat) then
  call write_density_matrix  ! write the density matrix
endif
call write_density  ! write the density
call write_magnetism  ! write the magnetism


call clean_density_matrix  ! clean the density matrix




return 
end subroutine 
