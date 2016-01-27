
! program to solve a mean field hamiltonian in a ribbon 
program tb90 
use inputs  ! global variables taken from input
implicit none 
integer :: stat

!!!!!!!! remove the DONE file !!!!!!!!
open(unit=1234, iostat=stat, file='DONE', status='old')
if (stat.eq.0) close(1234, status='delete')
!!!!!! create a RUNNING file !!!!!!
open(unit=1234,file='RUNNING')
close(1234)

call read_inputs ! read the input file


! calculate everything
if (do_scf) call selfconsistency  ! run SCF mean field calculation
if (do_bands) call bandstructure  ! calculate bandstructure
if (do_dos) call density_of_states  ! calculate density of states
if (do_dielectric) call dielectric_response  ! calculate the dielectric response function
if (do_denmat) call ground_density_matrix  ! calculate the dielectric response function
if (do_berry) call berry ! calculate berry stuff


!!!!!!!!!!! remove the RUNNING file  !!!!!!!!!!!
open(unit=1234, iostat=stat, file='RUNNING', status='old')
if (stat.eq.0) close(1234, status='delete')
!!!!!! create a DONE file !!!!!!
open(unit=1234,file='DONE')
close(1234)

end program 

