subroutine selfconsistency
use inputs  ! input variables
use system_variables  ! variables from this particular system 
use mean_field_variables  ! variables for mean field
use save_mf_routines ! routines to write scf results

implicit none 
integer :: j,i1,i2,i3,i4,i5 
logical :: dostop=.false.  ! skip the loop by the presence of a STOP file
integer :: stat ! status of the STOP file
integer :: iite

call read_inputs

hamiltonian_file='hamiltonian_0.in'
call read_hamil   ! obtain hamiltonian
hamiltonian_file='hamiltonian.in'
! onsite is modified along the SCF cycle, onsite_0 is the same always


! initialize mean_field is inside read_mf_op
if (mean_field_operators=='from_file') then
  call read_mf_op  ! read mean field matrices
else if (mean_field_operators=='hubbard') then
  call mean_field_operators_hubbard  ! read mean field matrices
else if (mean_field_operators=='hubbard-collinear') then
  call mean_field_operators_collinear_hubbard  ! read mean field matrices
else 
  write(*,*) 'Invalid option in mean_field_operators,\n Valid  &
              options are &
              \n    from_file  (read from mean_field_operators.in)&
              \n    hubbard   (generate assuming a Hubbard interaction)&
              \n',mean_field_operators
    stop
endif

! onsite_0 has been allocated in initialize mean_field
onsite_0 = onsite ! save the onsite term in the free one


! do not use energy window (at least for now...)
use_ewindow = .false.


call initialize_klist  ! obtain 1d klist

!!! This iniitalization before klist due to nkpoints redefinition
call initialize_eigenstates  ! initialize eigenstates
write(*,*) "Generated klist"


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! Begin guess of the hamiltonian !!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

if (mean_field_matrix=='from_file') then
  call read_mean_field  ! get the initial guess for the hamiltonian
else if (mean_field_matrix=='ferro') then
  call generate_ferro_mean_field  ! ferromagnetic initial configuration 
else if (mean_field_matrix=='antiferro') then
  call generate_antiferro_mean_field  ! antiferromagnetic initial configuration 
else if (mean_field_matrix=='random') then
  call generate_random_mean_field  ! random initial configuration 
else if (mean_field_matrix=='random_xz') then
  call generate_random_xz_mean_field  ! random initial configuration, XZ plane 
else
  write(*,*) 'Invalid option in mean_field_matrix,\n Valid  &
              options are &
              \n    from_file (read from mean_field.in) &
              \n    ferro (assume a ferromagnetic configuration) &
              \n    antiferro (assume a ferromagnetic configuration) &
              \n    random (generate a random initial guess) &
              \n',mean_field_matrix
  stop
endif

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! End guess of the hamiltonian !!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


scf_err = 1.d05   ! ensure that the cycle begins



call setup_mixing  ! initialice the mixing scheme


! remove the STOP file
open(unit=1234, iostat=stat, file='STOP', status='old')
if (stat.eq.0) close(1234, status='delete')


! open file to write energies
open(19,file="ENERGY.OUT")
write(19,*) "# TOTAL ENERGY,  ENERGY PER ORBITAL"
close(19) ! close the file


! and to write the error
open(23,file="ERROR.OUT")
write(23,*) "# ERROR"
close(23) ! close the file



iite = 1

do while (.true.)  ! do infinite loop 
  write(*,*) 'Iteration number',iite
  hamiltonian%onsite = onsite_0 + mf_ham ! add the mean field to the onsite matrix
  call eigen_scf   ! diagonalize in all the kpoints 
  if (is_collinear.and.(mean_field_operators=='hubbard-collinear').and. &
         ultrafast_mf) then
    call new_mf_ham_collinear   ! new MF hamiltonian 
    write(*,*) "Ultrafast calculation of mean field" 
  else ! conventional mean field calculation
    call new_mf_ham   ! new MF hamiltonian  
  endif
  call mixing    ! mix and get error
  iite = iite + 1 ! increase the iteration counter

!!!!!!!!!!!!!!!!
!! stop part !!!
!!!!!!!!!!!!!!!!
  if (scf_err < max_scf_err) exit   ! if convergence has been reached
  inquire(file="STOP", exist=dostop)  ! check if a stop file exists
  if (dostop) exit ! if stop file exist, end the selfconsistent loop
enddo 

call write_hamiltonian  ! save the final MF hamiltonian

if (save_mean_field)  then
  write (*,*) "Saved mean field"
  call write_mean_field  ! save the final SCF field
endif

call clean_mean_field ! deallocate mean field variables


return 
end subroutine 
