module bands_variables
use io_files
use expectation_values_variables
implicit none
integer :: nels_bands_op ! number of nonvanishing elements of the operator 
integer,save, allocatable :: ind_bands_op(:,:)  ! indices of bands operators
complex (kind=8),save, allocatable :: value_bands_op(:,:)  ! value of BOP
real (kind=8),allocatable :: vals_bands_ops(:,:)


contains
! subroutine to write the 1d bandstructure
  subroutine write_bands
  use inputs
  use system_variables
  use expectation_values_variables
  implicit none
  integer :: i,j,ii,k
  open(31,file=bands_file)   !!! open the file
  write(31,*) "#  TB90 band-structure"
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!! print d dimensional bands !!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    write(31,*) "#  TB90 system_dimension = ",hamiltonian%dimensionality
    ii = 0
    if (num_operators.eq.0) then
      ! zero dimensional case
      if (hamiltonian%dimensionality==0) then
        write(31,*) "#  index, eigenvalue"
        do j=1,klist%nev_found(1)  
          write(31,*) j,energies(j)
        enddo
      return ! skip the rest
      endif

      write(31,*) "#  k-index, eigenvalue"
      do i=1,klist%nkpoints ! write the kpointund in it
        do j=1,klist%nev_found(i)  
! write the eigenvalue sfound in the kpoint        
          ii = ii+1 ! increase the counter
          write(31,*) i,energies(ii)
        enddo
      enddo
    endif

    if (num_operators.gt.0) then

    ! 0 dimensional system
      if (hamiltonian%dimensionality==0) then
        write(31,*) "#  index, eigenvalue"
        do j=1,klist%nev_found(1)  
          write(31,*) j,energies(j),(exp_val_values(k,j) &
                 ,k=1,num_operators)
        enddo
      return ! skip the rest
      endif

      write(31,*) "#  k-index     eigenvalue    operator"
      do i=1,klist%nkpoints ! write the kpoint the number of eigenvalues found in it
        do j=1,klist%nev_found(i)  
! write the eigenvalue sfound in the kpoint        
          ii = ii+1 ! increase the counter
          write(31,*) i,energies(ii),(exp_val_values(k,ii) &
                 ,k=1,num_operators)
        enddo
      enddo
    endif
  close(31)
  return
  end subroutine






  subroutine write_klist_bands ! writes the klist in a file
  implicit none
  integer :: i,j
  real (kind=8) :: shift(system_dimension)
  real (kind=8) :: kp(system_dimension)
  shift(1) = expand_bz
  open(17,file='KPOINTS_BANDS.OUT')  ! open output file
  write(17,*) '# index, k in reciprocal, k in cartesian' 
  do i = 1, klist%nkpoints
      do j=1,klist%nev_found(i)
         ! vector in cartesian coordinates
           if (hamiltonian%has_lattice) then
             kp = matmul(hamiltonian%reciprocal_vectors(:,:),&
                         klist%kpoint(i,:)) 
        ! write in file
             write(17,*) i,(klist%kpoint(i,:)),(kp(:)) 
           else
             write(17,*) i,(klist%kpoint(i,:))
           endif 
      enddo
  enddo

  close(17)  ! close output file

  end subroutine write_klist_bands



! subroutine to setup the expectation values
  subroutine setup_bands_vals
  use inputs
  use system_variables
  implicit none
  ! deallocate if allocated
  if (allocated(vals_bands_ops)) deallocate(vals_bands_ops)
  ! if at least one operator, allocate
  if (num_operators.gt.0)  allocate(vals_bands_ops(num_operators,num_wf))


  return
  end subroutine setup_bands_vals


! read klist from input file
  subroutine read_klist
  use inputs
  implicit none
  integer ik,jk,id,kk,nk,ak

  open(11,file='klist.in')
  write(*,*) 'Reading klist from klist.in' 
  read(11,*) nkpoints  ! number of kpoints
  write(*,*) nkpoints

  klist%nkpoints = nkpoints  ! store number of kpoints in the class

  ! if already allocated deallocate
  if (allocated(klist%kpoint)) deallocate(klist%kpoint)
  if (allocated(klist%weight)) deallocate(klist%weight)
  if (allocated(klist%kname)) deallocate(klist%kname)
  if (allocated(klist%nev_found)) deallocate(klist%nev_found)
  ! allocate the kpoints
  allocate(klist%kpoint(klist%nkpoints,hamiltonian%dimensionality))
  ! allocate the name of the kpoints
  allocate(klist%kname(klist%nkpoints))
  allocate(klist%weight(klist%nkpoints))
  allocate(klist%nev_found(klist%nkpoints))
  ! create kpoint list

  ! read each of the kpoints
  do ik=1,nkpoints  ! loop over kpoints
    read(11,*) (klist%kpoint(ik,jk),jk=1,hamiltonian%dimensionality)
  enddo

  klist%nev_found = 0   ! initialice the number of eigenvalues per k
  klist%weight(:) = 1.d00  ! initialize the weight of the kpoint
  klist%weight = klist%weight/sum(klist%weight)

  close(11) ! close klist.in file!!!!

  return
  end subroutine read_klist


  subroutine write_klist ! writes the klist in a file
  implicit none
  integer :: i,j
  real (kind=8) :: shift(system_dimension)
  real (kind=8) :: kp(system_dimension)
  shift(1) = expand_bz
  open(17,file='KPOINTS.OUT')  ! open output file
  write(17,*) '# index, k in reciprocal, k in cartesian'
  do i = 1, klist%nkpoints
         write(17,*) klist%kpoint(i,:)
  enddo

  close(17)  ! close output file
  return
  end subroutine write_klist


  subroutine write_bands_kpoint(ik,clean)
  use inputs
  use system_variables
  use expectation_values_variables
  implicit none
  integer :: i,j,ii,k
  logical :: clean
  integer :: ik ! index of the kpoint
 
  if (clean) then ! first iteration, just clean the file 
    open(31,file=bands_file)   !!! open the file
    write(31,*) "#  TB90 system_dimension = ",hamiltonian%dimensionality
    close(31)
    return ! end the function
  endif


  open(31,file=bands_file,ACCESS = 'APPEND')   !!! open the file
  write(31,*) "#  TB90 band-structure"
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!! print d dimensional bands !!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ii = 0
    if (num_operators.eq.0) then
      ! zero dimensional case


      if (hamiltonian%dimensionality==0) then
        do j=1,klist%nev_found(1)  
          write(31,*) j,energies_k(j)
        enddo
      return ! skip the rest
      endif



        do j=1,num_wf
! write the eigenvalue sfound in the kpoint        
          write(31,*) ik,energies_k(j)
        enddo
    endif

    if (num_operators.gt.0) then

    ! 0 dimensional system
      if (hamiltonian%dimensionality==0) then
        do j=1,num_wf 
          write(31,*) j,energies_k(j),(exp_val_values(k,j) &
                 ,k=1,num_operators)
        enddo
      return ! skip the rest
      endif





      do j=1,num_wf  
! write the eigenvalue sfound in the kpoint        
          write(31,*) ik,energies_k(j),(exp_val_values(k,j) &
                 ,k=1,num_operators)
      enddo
    endif
  close(31)
  return
  end subroutine







end module
