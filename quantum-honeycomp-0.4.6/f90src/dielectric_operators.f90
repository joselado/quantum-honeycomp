! variables for calculation of dielectric response
module dielectric_operators
use inputs 
use system_variables
use sparse ! sparse library
use dielectric_variables
  

contains

  ! sets up all the matrices for linear response
      ! dielectric response
      ! SzSz response
      ! SxSy response
  subroutine setup_all_ab_matrices
  implicit none

  if (allocated(ab_matrices))  deallocate(ab_matrices)  ! deallocate 
  num_chi_ab = 3 ! number of pairs of operators
  allocate(ab_matrices(num_chi_ab))  ! allocate array with pairs 

  ! Generate matrices for dielectric response
  call sparse_id(norbitals,ab_matrices(1)%a)
  call sparse_id(norbitals,ab_matrices(1)%b)
  ab_matrices(1)%label = '$\rho\rho$  '
  ! Generate matrices for SzSz response
  call sparse_sz(norbitals,ab_matrices(2)%a)
  call sparse_sz(norbitals,ab_matrices(2)%b)
  ab_matrices(2)%label = '$S_zS_z$  '
  ! Generate matrices for SpSm response
  call sparse_sx(norbitals,ab_matrices(3)%a)
  call sparse_sy(norbitals,ab_matrices(3)%b)
  ab_matrices(3)%label = '$S_xS_y$  '

  return
  end subroutine setup_all_ab_matrices






  subroutine setup_ab_global_density
  implicit none
  if (allocated(ab_matrices))  deallocate(ab_matrices)  ! deallocate 
  num_chi_ab = 1 ! number of pairs of operators
  allocate(ab_matrices(num_chi_ab))  ! allocate array with pairs 
  ! Generate matrices for dielectric response
  call sparse_id(norbitals,ab_matrices(1)%a)
  call sparse_id(norbitals,ab_matrices(1)%b)
  ab_matrices(1)%label = '$\rho\rho$  '
  return
  end subroutine setup_ab_global_density












  subroutine setup_local_density ! sets up all the pair of matrices rhoirhoj
  implicit none
  integer :: isite,jsite,ii ! counter for site


  if (allocated(ab_matrices))  deallocate(ab_matrices)  ! deallocate 
  ii = 1
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !! with spin polarization !!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!

  if (hamiltonian%has_spinpol) then
    num_chi_ab = norbitals/2 ! number of pairs of operators
    num_chi_ab = num_chi_ab*num_chi_ab ! number of pairs of operators
    allocate(ab_matrices(num_chi_ab))  ! allocate array with pairs 
    do isite = 1,norbitals/2   ! loop over the sites
      do jsite = 1,norbitals/2   ! loop over the sites
        call sparse_local_rho(isite,ab_matrices(ii)%a)
        call sparse_local_rho(jsite,ab_matrices(ii)%b)
        ab_matrices(ii)%label = '  ii  '
        ii = ii + 1
      enddo
    enddo
  endif


  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !! without spin polarization !!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  if (.not.hamiltonian%has_spinpol) then
    num_chi_ab = norbitals ! number of pairs of operators
    num_chi_ab = num_chi_ab*num_chi_ab ! number of pairs of operators
    allocate(ab_matrices(num_chi_ab))  ! allocate array with pairs 
    do isite = 1,norbitals   ! loop over the sites
      do jsite = 1,norbitals   ! loop over the sites
        call sparse_local_rho_nonsp(norbitals,isite,ab_matrices(ii)%a)
        call sparse_local_rho_nonsp(norbitals,jsite,ab_matrices(ii)%b)
        ab_matrices(ii)%label = '  ii  '
        ii = ii + 1
      enddo
    enddo
  endif
  return
  end subroutine setup_local_density




  subroutine setup_local_xy_spin ! sets up all the pair of matrices rhoirhoj
  implicit none
  integer :: isite,jsite,ii ! counter for site
  if (allocated(ab_matrices))  deallocate(ab_matrices)  ! deallocate 
  ii = 1
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !! with spin polarization !!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!
  if (hamiltonian%has_spinpol) then
    num_chi_ab = norbitals/2 ! number of pairs of operators
    num_chi_ab = num_chi_ab*num_chi_ab ! number of pairs of operators
    allocate(ab_matrices(num_chi_ab))  ! allocate array with pairs 
    do isite = 1,norbitals/2   ! loop over the sites
      do jsite = 1,norbitals/2   ! loop over the sites
        call sparse_local_sp(isite,ab_matrices(ii)%a)
        call sparse_local_sm(jsite,ab_matrices(ii)%b)
        ab_matrices(ii)%label = '  ii  '
        ii = ii + 1
      enddo
    enddo
  endif
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !! without spin polarization !!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  if (.not.hamiltonian%has_spinpol) then
    write(*,*) 'No spinpol in magnetic response, stoping'
    stop  
  endif
  return
  end subroutine setup_local_xy_spin




  subroutine setup_local_zz_spin ! sets up all the pair of matrices rhoirhoj
  implicit none
  integer :: isite,jsite,ii ! counter for site
  if (allocated(ab_matrices))  deallocate(ab_matrices)  ! deallocate 
  ii = 1
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !! with spin polarization !!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!
  if (hamiltonian%has_spinpol) then
    num_chi_ab = norbitals/2 ! number of pairs of operators
    num_chi_ab = num_chi_ab*num_chi_ab ! number of pairs of operators
    allocate(ab_matrices(num_chi_ab))  ! allocate array with pairs 
    do isite = 1,norbitals/2   ! loop over the sites
      do jsite = 1,norbitals/2   ! loop over the sites
        call sparse_local_sz(isite,ab_matrices(ii)%a)
        call sparse_local_sz(jsite,ab_matrices(ii)%b)
        ab_matrices(ii)%label = '  ii  '
        ii = ii + 1
      enddo
    enddo
  endif
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !! without spin polarization !!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  if (.not.hamiltonian%has_spinpol) then
    write(*,*) 'No spinpol in magnetic response, stoping'
    stop  
  endif
  return
  end subroutine setup_local_zz_spin



  ! read response operators from file
  subroutine setup_response_operators_from_file
  implicit none
  integer :: i,iatom,jatom  ! counter for pairs
  character (len=10) :: itype,jtype



  open(71,file='response_operators.in')
  ! inicitalice pairs of operators
  read(71,*)  ! number of pairs of operators
  read(71,*)  num_chi_ab
  read(71,*)  ! atom index, type, atom index, type
  if (allocated(ab_matrices))  deallocate(ab_matrices)  ! deallocate 
  allocate(ab_matrices(num_chi_ab))  ! allocate array with pairs 
  do i=1,num_chi_ab
    read(71,*)  iatom,itype,jatom,jtype
    ! Sz correlation
    if (itype=='Sz') call sparse_local_sz(iatom,ab_matrices(i)%a)    
    if (jtype=='Sz') call sparse_local_sz(jatom,ab_matrices(i)%b)    
    ! Sp correlation
    if (itype=='S+') call sparse_local_sp(iatom,ab_matrices(i)%a)    
    if (jtype=='S+') call sparse_local_sp(jatom,ab_matrices(i)%b)    
    ! Sm correlation
    if (itype=='S-') call sparse_local_sm(iatom,ab_matrices(i)%a)    
    if (jtype=='S-') call sparse_local_sm(jatom,ab_matrices(i)%b)    
        ab_matrices(i)%label = itype//jtype
  enddo
  close(71)


  end subroutine setup_response_operators_from_file



end module
