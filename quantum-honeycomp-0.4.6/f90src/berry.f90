
! main routine for the calculation of topological properties

subroutine berry
use inputs
use system_variables
use berry_variables
use bands_variables
implicit none


call read_inputs ! reread the input file
call read_hamil ! reread hamiltonian


! For two dimensional, select the path
if (hamiltonian%dimensionality==2)  then
  if (.not.has_tr) then ! if the system does not hve TR do nothing
    if (klist_berry == "default") then
      call initialize_klist ! initialize klist
    else if (klist_berry=="from_file") then
      call read_klist
    else
      write(*,*) "Unrecognised option in klist_berry, STOPPING"
      stop
    endif
  else
    write(*,*) "Z2 invariant calculation"
  endif
endif


! For one dimensional, default path
if (hamiltonian%dimensionality==1)  then
  call initialize_klist ! initialize klist
endif  
call initialize_berry ! initialize berry variables


if (hamiltonian%dimensionality==1) then
  write(*,*) '======================='
  write(*,*) 'Berry phase calculation'
  write(*,*) '======================='
  ! calculate the Berry phase for 1d systems
  call calculate_berry_phase(kpoints_berry,nkpoints_berry,berry_phase)
  call write_berry_phase   ! write it in file
else if (hamiltonian%dimensionality==2) then
  write(*,*) '==========================='
  write(*,*) 'Berry curvature calculation'
  write(*,*) '==========================='
  if (has_tr) then  ! if it has time reversal symmetry
    call calculate_z2_invariant ! calculate the Z2 for a TRS hamiltonian
  else  ! Berry curvature calculation
    call calculate_berry_curvature
    call write_berry_curvature   ! write it in file
  endif
else
  write(*,*) 'No Berry stuff implemented for dimension', &
                 hamiltonian%dimensionality
endif 

write(*,*) "============"
write(*,*) "End of Berry"
write(*,*) "============"

call clean_waves  ! clean all the wavefunctions

return
end subroutine berry








