

! variables and routines for berry calculations


module berry_variables
use system_variables
use inputs
use io_files
use bands_variables, only : write_klist, read_klist
implicit none
integer :: num_occ_bands ! number of occupied bands
real (kind=8) :: berry_phase ! berry phase for 1d systems
! for 2d systems
real (kind=8), allocatable :: berry_curvature(:) ! berry curvature
real (kind=8), allocatable :: kpoints_berry(:,:) ! path for berry calculation
integer :: nkpoints_berry ! number of kpoints for the path
complex (kind=8), allocatable :: wf_k0(:,:) ! wavefunction at k=0
complex (kind=8), allocatable :: wf_k_old(:,:) ! wavefunction at the old kpoint
! berry curvature for 2d systems



contains

  subroutine initialize_berry ! initialize the berry variables
  implicit none
  
  ! for 1d no allocatation is needed

  ! initialize berry curvature for 2d
  if (hamiltonian%dimensionality==2) then
    if (allocated(berry_curvature)) deallocate(berry_curvature)
    allocate(berry_curvature(klist%nkpoints)) ! allocate berry curvature array
    berry_curvature(:) = 0.d00 ! initialice to 0
! allocate the small kpath for Berry curv calculation
    if (allocated(kpoints_berry)) deallocate(kpoints_berry)
    nkpoints_berry = 4 ! number of kpoints (small path!!!)
    allocate(kpoints_berry(nkpoints_berry,hamiltonian%dimensionality))
  endif

  ! initialize a set of wf for the first kpoint in 1d
  if (hamiltonian%dimensionality==1) then
! allocate the kpoint path
    if (allocated(kpoints_berry)) deallocate(kpoints_berry)
    allocate(kpoints_berry(klist%nkpoints,hamiltonian%dimensionality))
    kpoints_berry(:,:) = klist%kpoint(:,:)  ! use as values the klist
    nkpoints_berry = klist%nkpoints ! number of kpoints
  endif

! allocate first wavefunction
  if (allocated(wf_k0)) deallocate(wf_k0)
  allocate(wf_k0(hamiltonian%norbitals,hamiltonian%norbitals)) 
! allocate old wavefunction
  if (allocated(wf_k_old)) deallocate(wf_k_old)
  allocate(wf_k_old(hamiltonian%norbitals,hamiltonian%norbitals)) 



  call initialize_eigenvectors_k ! initialize wf_k

  return
  end subroutine initialize_berry

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine write_berry_phase
  implicit none
  open(61,file=berry_phase_file) ! open the file
  write(61,*) "# Berry phase for a 1d system"
  write(61,*) berry_phase
  close(61) ! close the file
  end subroutine  write_berry_phase


  subroutine write_berry_curvature
  implicit none
  integer :: ik
  real (kind=8) :: kp(2)
  open(61,file=berry_curv_file) ! open the file
  write(61,*) "# Berry curvature for a 2d system"
  
  if (hamiltonian%has_lattice) then
    write(61,*) "# kvectors in real coordinates"
  else 
    write(61,*) "# kvectors in normal coordinates"
  endif  

  do ik=1,klist%nkpoints
    ! get the kpoint in desired coordinates
    if (hamiltonian%has_lattice) then
       kp(:) = matmul(hamiltonian%reciprocal_vectors(:,:),&
                           klist%kpoint(ik,:))
    else
       kp(:) = klist%kpoint(ik,:)
    endif
    ! now write the result
    write(61,*) kp(:),berry_curvature(ik)
  enddo
  close(61) ! close the file
  end subroutine  write_berry_curvature


 

 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! calculates the berry phase for a 1d system, summing over all
  ! the occupied bands by using the determinant algorithm
 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine calculate_berry_phase(kpberry,nk,phiberry) 
  implicit none
  integer, intent(in) :: nk  ! number of kpoints
  ! path for berry calculation
  real (kind=8), intent(in) :: kpberry(nk,hamiltonian%dimensionality) 
  ! berry phase 
  real (kind=8), intent(out) :: phiberry 
  ! storage of wavefunctions
  complex (kind=8), dimension(norbitals,norbitals) :: wf0,wf1,wf2
  ! storage of energies
  real (kind=8), dimension(norbitals) :: e0,e1,e2
  integer :: nbnd ! number of occupied bands
! matrix for berry calculations
  complex (kind=8), allocatable :: berrymat(:,:)  
  complex (kind=8), allocatable :: tmpmat(:,:)  
  complex (kind=8) :: zdet ! determinant of the matrix multiplication 
  integer :: ik,ib,i,j ! counter for kpoints
  ! first kpoint
  call kwaves(kpberry(1,:),wf0,e0) ! get eigenvectors
  nbnd = 0

  ! get the number of occupied bands
  do ib =1, hamiltonian%norbitals 
    if (e0(ib).gt.0.d00) exit ! break if eigenvalue from cond band
    nbnd = nbnd + 1 ! increase the number of occ bands
  enddo
  if (nbnd.eq.0) then ! if no bands occupied, return
    phiberry = 0.d00 ! zero berry phase
    return 
  endif
  ! initialize berry matrix
  if (allocated(berrymat)) deallocate(berrymat)
  if (allocated(tmpmat)) deallocate(tmpmat)
  allocate(berrymat(nbnd,nbnd))
  allocate(tmpmat(nbnd,nbnd))
  berrymat = (0.d00,0.d00)
  do ib=1,nbnd
    berrymat(ib,ib) = (1.d00,0.d00)
  enddo
  wf1 = wf0 ! first point is also the last
  ! loop over kpoints_berry
  do ik = 2, nk
    tmpmat = (0.d00,0.d00)
    call kwaves(kpberry(ik,:),wf2,e0) ! get eigenvectors
    ! build berry matrix by performing brakets
!$omp parallel default (shared) private (i,j)
!$omp do
    do i=1,nbnd
      do j=1,nbnd
         tmpmat(i,j) = sum(conjg(wf1(i,:))*wf2(j,:)) ! <old,i|new,j>
      enddo
    enddo
!$omp end do
!$omp end parallel
    wf1 = wf2 ! the old in the next iteration is the new in this one
    ! perform matrix multiplication
    berrymat = matmul(berrymat,tmpmat) ! rho_old*rho_new
  enddo
  ! perform the last multiplication
  tmpmat = (0.d00,0.d00)

!$omp parallel default (shared) private (i,j)
!$omp do
  do i=1,nbnd
    do j=1,nbnd
       tmpmat(i,j) = sum(conjg(wf2(i,:))*wf0(j,:)) ! <old,i|new,j>
    enddo
  enddo
!$omp end do
!$omp end parallel

  berrymat = matmul(berrymat,tmpmat) ! rho_old*rho_new
  call determinant(nbnd,berrymat,zdet)
  ! now calculate the phase of the complex number
  phiberry = atan2(aimag(zdet),real(zdet)) ! calculate the berry phase
  phiberry = phiberry/(3.141592*2.d00) ! normalize by 2pi
  ! deallocate matrices
  deallocate(berrymat)
  deallocate(tmpmat)
  return
  end subroutine calculate_berry_phase

 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine calculate_berry_curvature ! initialize the berry variables
  implicit none
  integer :: ik,i,j ! counter for kpoints
  real (kind=8) :: phiberry

! allocate berry curvature
  if (allocated(berry_curvature)) deallocate(berry_curvature)
  allocate(berry_curvature(klist%nkpoints)) ! allocate berry curvature array

! loop over kpoints
  nkpoints_berry = 4 ! number of kpoints (small path!!!)
  if (allocated(kpoints_berry)) deallocate(kpoints_berry)
  allocate(kpoints_berry(nkpoints_berry,hamiltonian%dimensionality))
!$omp parallel default (shared) private (ik,kpoints_berry,phiberry)
!$omp do
  do ik=1,klist%nkpoints  
    ! now create a small path enclosing that kpoint
    ! upper right
    kpoints_berry(1,1) = klist%kpoint(ik,1) + dk_becurv
    kpoints_berry(1,2) = klist%kpoint(ik,2) + dk_becurv
    ! upper left
    kpoints_berry(2,1) = klist%kpoint(ik,1) - dk_becurv
    kpoints_berry(2,2) = klist%kpoint(ik,2) + dk_becurv
    ! lower left
    kpoints_berry(3,1) = klist%kpoint(ik,1) - dk_becurv
    kpoints_berry(3,2) = klist%kpoint(ik,2) - dk_becurv
    ! lower right
    kpoints_berry(4,1) = klist%kpoint(ik,1) + dk_becurv
    kpoints_berry(4,2) = klist%kpoint(ik,2) - dk_becurv
    call calculate_berry_phase(kpoints_berry,4,phiberry) 
    berry_curvature(ik) = phiberry ! store each berry phase
  enddo
!$omp end do
!$omp end parallel

  ! normalize the berry curvature
  berry_curvature = berry_curvature/((dk_becurv*dk_becurv)*4.d00)

  return
  end subroutine calculate_berry_curvature



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! calculates the Z2 invariant for a system with TR symmetry
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine calculate_z2_invariant
  implicit none
  real (kind=8) :: phiberry
  real (kind=8) :: hchern
  integer :: ik
  character (len=30) :: kmesh
  integer :: z2 

  ! create the klist
! creates a klist which goes in half of the BZ 
  kmesh = "bulk"
  call create_half_bz_klist(kmesh) 
! calculate the berry curvature over 1/2 of the BZ
  call calculate_berry_curvature   
  hchern = sum(berry_curvature)/klist%nkpoints
! creates a klist which goes in half of the BZ, in the boundary 
  kmesh = "boundary"
!  call create_half_bz_klist(kmesh) 
  call read_klist
  if (allocated(kpoints_berry)) deallocate(kpoints_berry)
  allocate(kpoints_berry(klist%nkpoints,hamiltonian%dimensionality))
  do ik=1,klist%nkpoints  ! copy points for using as input in subroutine
    kpoints_berry(ik,:) = klist%kpoint(ik,:)  ! copy the point
  enddo
! calculate the berry phase over 1/2 of the BZ
  call integrate_berry_connection(kpoints_berry,klist%nkpoints,phiberry)

  z2 = abs(mod(nint(hchern-phiberry),2))  ! calculte z2 invariant
  write(*,*) "Z2 is", z2

  !!!!!!!!!!!!!!!!!!!
  !! write in file !!
  !!!!!!!!!!!!!!!!!!!
  call write_klist
  open(79,file="Z2.OUT")
    write(79,*) "#  Chern contribution = ",hchern
    write(79,*) "#  Berry connection contribution = ",phiberry
    if (z2==1) then
      write(79,*) "#  The system is TOPOLOGICAL"
    else
      write(79,*) "#  The system is TRIVIAL"
    endif
    write(79,*) z2
  close(79)
  return
  end subroutine calculate_z2_invariant




! integrate the berry connection along a path
  subroutine integrate_berry_connection(kpberry,nk,phiberry) 
  implicit none
  integer, intent(in) :: nk  ! number of kpoints
  ! path for berry calculation
  real (kind=8), intent(in) :: kpberry(nk,hamiltonian%dimensionality) 
  ! berry phase 
  real (kind=8), intent(out) :: phiberry 
  ! storage of wavefunctions
  complex (kind=8), dimension(norbitals,norbitals) :: wf0,wf1,wf2
  ! storage of energies
  real (kind=8), dimension(norbitals) :: e0,e1,e2
  integer :: nbnd ! number of occupied bands
! matrix for berry calculations
  complex (kind=8), allocatable :: tmpmat(:,:)  

  real (kind=8) :: sumberry 

  complex (kind=8) :: zdet ! determinant of the matrix multiplication 
  integer :: ik,ib,i,j ! counter for kpoints
  ! first kpoint
  call kwaves(kpberry(1,:),wf0,e0) ! get eigenvectors
  nbnd = 0

  sumberry = (0.d00,0.d00) 
  ! get the number of occupied bands
  do ib =1, hamiltonian%norbitals 
    if (e0(ib).gt.0.d00) exit ! break if eigenvalue from cond band
    nbnd = nbnd + 1 ! increase the number of occ bands
  enddo
  ! initialize berry matrix
  if (allocated(tmpmat)) deallocate(tmpmat)
  allocate(tmpmat(nbnd,nbnd))
  wf1 = wf0 ! first point is also the last
  ! loop over kpoints_berry
  do ik = 2, nk
    tmpmat = (0.d00,0.d00)
    call kwaves(kpberry(ik,:),wf2,e0) ! get eigenvectors
    ! build berry matrix by performing brakets
!$omp parallel default (shared) private (i,j)
!$omp do
    do i=1,nbnd
      do j=1,nbnd
         tmpmat(i,j) = sum(conjg(wf1(i,:))*wf2(j,:)) ! <old,i|new,j>
      enddo
    enddo
!$omp end do
!$omp end parallel
    call determinant(nbnd,tmpmat,zdet)
    sumberry = sumberry+atan2(aimag(zdet),real(zdet))
    wf1 = wf2 ! the old in the next iteration is the new in this one
    ! perform matrix multiplication
  enddo
  ! perform the last multiplication
  tmpmat = (0.d00,0.d00)

!$omp parallel default (shared) private (i,j)
!$omp do
  do i=1,nbnd
    do j=1,nbnd
       tmpmat(i,j) = sum(conjg(wf2(i,:))*wf0(j,:)) ! <old,i|new,j>
    enddo
  enddo
!$omp end do
!$omp end parallel
  call determinant(nbnd,tmpmat,zdet)
  sumberry = sumberry+atan2(aimag(zdet),real(zdet))
  ! now calculate the phase of the complex number
  sumberry = sumberry/(2.d00*3.141592)
!  phiberry = mod(real(sumberry),2.d00)
  phiberry = real(sumberry)
!  write(*,*) "Berry phase =",sumberry,mod(real(sumberry),2.d00)

  ! deallocate matrices
  deallocate(tmpmat)
  return
  end subroutine integrate_berry_connection


end module



