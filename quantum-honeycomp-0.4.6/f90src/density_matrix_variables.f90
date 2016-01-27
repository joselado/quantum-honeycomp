module density_matrix_variables
use system_variables ! for know norbitals
use inputs ! for know norbitals
use system_variables ! for know norbitals
implicit none
complex (kind=8),allocatable :: den_mat(:,:)  ! density matrix
real (kind=8),allocatable :: density(:)  ! density of the system
real (kind=8),allocatable :: magnetism(:,:)  ! magnetism of the system

contains

  !!!! for initializing the den_mat !!!
  subroutine initialize_density_matrix
  if (allocated(den_mat)) deallocate(den_mat)  ! deallocate if allocated
  allocate(den_mat(norbitals,norbitals))  ! allocate density matrix
  if (allocated(density)) deallocate(density)  ! deallocate if allocated
  allocate(density(norbitals/2))  ! allocate density
  if (allocated(magnetism)) deallocate(magnetism)  ! deallocate if allocated
  allocate(magnetism(norbitals/2,3))  ! allocate magnetism
  return
  end subroutine initialize_density_matrix

  !!! for cleaning the den_mat !!!
  subroutine clean_density_matrix
  if (allocated(den_mat)) deallocate(den_mat)  ! deallocate if allocated
  if (allocated(density)) deallocate(density)  ! deallocate if allocated
  if (allocated(magnetism)) deallocate(magnetism)  ! deallocate if allocated
  return
  end subroutine clean_density_matrix

  !!! for calculating the den_mat !!!
  subroutine calculate_density_matrix
  integer :: i,j,iw,n
  complex (kind = 8) :: wf_tmp(norbitals) ! temporal wf storage
  complex (kind=8) :: occ ! ocupation of the eigenvalue
  den_mat(:,:) = (0.d00,0.d00)  ! initialize to 0
  do iw=1,num_wf ! loop over wavefunction
!    if (energies(iw).lt.0.d00) then  ! if occupied calculate
    if (energies(iw).lt.dabs(temperature_denmat)) then  ! if occupied calculate
!      write(*,*) energies(iw),iw
      occ = 1.d00 ! default is fully occupied
      if (energies(iw).gt.(-dabs(temperature_denmat))) then
        occ = energies(iw)/temperature_denmat ! to interval (-1,1) 
        occ = 1.d00 - occ ! to interval (2,0)
        occ = occ*5.d-1 ! to interval (1,0)
      endif
      wf_tmp(:) = wf(iw,:)  ! store in handy variable
      if (write_denmat) then ! calculate full density matrix
        do i =1,hamiltonian%norbitals  ! loop over i
          do j=1,hamiltonian%norbitals  ! loop over j
          ! add contribution
            den_mat(i,j) = den_mat(i,j) + wf_tmp(j)*conjg(wf_tmp(i)) 
          enddo
        enddo
      else  ! if only magnetism is needed
        do n =1,hamiltonian%norbitals/2  ! loop over atoms
          i = 2*n-1
          j = 2*n-1
          den_mat(i,j) = den_mat(i,j) + wf_tmp(j)*conjg(wf_tmp(i)) 
          i = 2*n
          j = 2*n-1
          den_mat(i,j) = den_mat(i,j) + wf_tmp(j)*conjg(wf_tmp(i)) 
          i = 2*n-1
          j = 2*n
          den_mat(i,j) = den_mat(i,j) + wf_tmp(j)*conjg(wf_tmp(i)) 
          i = 2*n
          j = 2*n
          den_mat(i,j) = den_mat(i,j) + wf_tmp(j)*conjg(wf_tmp(i)) 
        enddo ! close loop over atoms
      endif  ! close reduced matrix if
    endif ! close occupied if
  enddo ! close loop loop over wavefunctions

  ! calculate density and magnetism
  do i=1,hamiltonian%norbitals/2
    j = 2*i
    magnetism(i,1) = den_mat(j-1,j) + den_mat(j,j-1) ! Mx
    magnetism(i,2) = -aimag(den_mat(j-1,j) - den_mat(j,j-1)) !My
    magnetism(i,3) = den_mat(j-1,j-1) - den_mat(j,j) !Mz
    density(i) = den_mat(j-1,j-1) + den_mat(j,j)
  enddo
  return
  end subroutine calculate_density_matrix

  !!! write density matrix !!!
  subroutine write_density_matrix
  integer :: i,j
  open(49,file="DENSITY_MATRIX.OUT")  ! file where to write
  do i=1,norbitals
    do j = 1,norbitals
      write(49,*) i,j,real(den_mat(i,j)),aimag(den_mat(i,j))
    enddo
  enddo
  
  close(49)

  return
  end subroutine write_density_matrix


  ! write file with the density
  subroutine write_density
  integer :: i,j
  open(49,file=density_file)  ! file where to write
  write(49,*) "# Index (without spin)         Density" 
  do i=1,norbitals/2
    write(49,*) i,density(i)
  enddo  
  close(49)
  write(*,*) "Saved density in  ",density_file 
  return
  end subroutine write_density



  ! write file with the magnetism
  subroutine write_magnetism
  integer :: i,j
  open(49,file=magnetism_file)  ! file where to write
  write(49,*) "# Index (without spin)         Mx        My      Mz" 
  do i=1,norbitals/2
    write(49,*) i,magnetism(i,1),magnetism(i,2),magnetism(i,3)
  enddo  
  close(49)
  write(*,*) "Saved magnetism in  ",magnetism_file 
  return
  end subroutine write_magnetism


  ! rotate the hamiltonian to a basis with z magnetization






end module density_matrix_variables
