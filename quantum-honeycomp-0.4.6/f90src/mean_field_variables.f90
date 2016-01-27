! this module stores the system variables as hamiltonian, mean field operators
! wavefunctions,etc

module mean_field_variables
use sparse ! use sparse library
implicit none

  type ab_mf
    type(spmatrix) :: a,b  ! operators for the mean field
    complex (kind=8) :: lambda ! coupling betweeen operators
  endtype ab_mf

type(ab_mf), allocatable :: ab_matrices(:) ! ab matrices for mean field

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! Mean field operator variables !!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
integer, save :: num_AB ! number of mean field operators
! elements of mean field operators
  ! first index is index of element, second is index of operator
complex (kind=8), allocatable, save :: mf_A(:,:),mf_B(:,:)
! first are ij, second index is index of element, third is index of operator
! index which overlap
integer, allocatable, save :: nv_ind_A(:,:,:),nv_ind_B(:,:,:)
! indexes
integer, allocatable, save :: ind_mfmata(:,:,:),nind_mfmata(:)
integer, allocatable, save :: ind_mfmatb(:,:,:),nind_mfmatb(:)
integer, save :: max_mf_el=4 ! maximun non vanishing elements
!# of elements in A,B
integer, allocatable, save :: num_mf_el_A(:),num_mf_el_B(:)
complex (kind=8), allocatable, save :: lambda_AB(:) ! AB coupling 

!!!!!!!!!!!!!!!!!!!
!! SCF variables !!
!!!!!!!!!!!!!!!!!!!
complex (kind=8), allocatable, save :: mf_ham(:,:) ! mean field hamiltonian
complex (kind=8), allocatable, save :: onsite_0(:,:) ! free onsite hamiltonian
complex (kind=8), allocatable, save :: old_mf_ham(:,:,:) ! old MF hamiltonians
real (kind=8), allocatable, save :: coef(:) ! coefficients for mixing old MF
real (kind=8), save :: scf_err  ! SCF error
real (kind=8), save :: ene_tot ! total energy


contains
! initialize mean field stuff
  subroutine initialize_mean_field
  use inputs
  use system_variables
  implicit none

  ! clean everything
!  call clean_mean_field

!  allocate(ab_matrices(num_AB)) ! number of non AB pairs

  ! old mean fields
  allocate(mf_ham(norbitals,norbitals))  ! allocate mean field
  allocate(onsite_0(norbitals,norbitals))  ! allocate free onsite
  allocate(old_mf_ham(num_old_ham,norbitals,norbitals))  ! old mean fields
  allocate(coef(num_old_ham))
  ! AB matrices
  allocate(num_mf_el_A(num_AB)) ! number of non vanishing elements in A
  ! A matrix
  allocate(mf_A(max_mf_el,num_AB))  ! values on the entries of A matrices
  allocate(nv_ind_A(2,max_mf_el,num_AB)) ! indexes of the previous value
  ! B matrix
  allocate(num_mf_el_B(num_AB)) ! number of non vanishing elements in B
  allocate(mf_B(max_mf_el,num_AB))  ! values on the entries of B matrices
  allocate(nv_ind_B(2,max_mf_el,num_AB)) ! indexes of the previous value
  ! coupling
  allocate(lambda_AB(num_AB))  ! mean field couplings
 
  ! initialize everything to zero
  mf_ham = (0.d00,0.d00)
  old_mf_ham = (0.d00,0.d00)
  coef = 0.d00
  num_mf_el_A = 0
  num_mf_el_B = 0
  nv_ind_A = 0
  nv_ind_B = 0
  mf_A = (0.d00,0.d00)
  mf_B = (0.d00,0.d00)
  write(*,*) "Allocated mean field variables"
  write(*,*)
  return
  end subroutine initialize_mean_field


  ! clear mean field stuff
  subroutine clean_mean_field
  implicit none
  if (allocated(mf_ham)) deallocate(mf_ham)  ! mean field
  if (allocated(onsite_0)) deallocate(onsite_0)  ! free onsite
  if (allocated(old_mf_ham)) deallocate(old_mf_ham) ! old mean field
  if (allocated(coef)) deallocate(coef)  ! mixing coefficients
  if (allocated(num_mf_el_A)) deallocate(num_mf_el_A)
  if (allocated(num_mf_el_B)) deallocate(num_mf_el_B)
  if (allocated(mf_A)) deallocate(mf_A)
  if (allocated(mf_B)) deallocate(mf_B)
  if (allocated(nv_ind_A)) deallocate(nv_ind_A)
  if (allocated(nv_ind_B)) deallocate(nv_ind_B)
  if (allocated(lambda_AB)) deallocate(lambda_AB)
  write(*,*) "Deallocated mean field variables"
  write(*,*)
  return
  end subroutine clean_mean_field

! READ mena field
! initial guess for the mean field hamiltonian 
  subroutine read_mean_field
  use inputs
  use system_variables
  implicit none
  
  real (kind=8) :: hr,hi 
  integer :: i,j,nel,i1,i2
  
  ! initialize 
  old_mf_ham=(0.d00,0.d00) 
  mf_ham=(0.d00,0.d00) 
  
  ! open file and read first guess 
  open(37,file='mean_field.in') 
  ! read only nonvanishing elements 
  read(37,*) 
  read(37,*) nel 
  read(37,*) 
  do i1=1,nel 
    read(37,*) i,j,hr,hi 
    if ((i.gt.norbitals).or.(j.gt.norbitals)) then
    write(*,*) " Error at reading mean_field.in, requesting element",i,j
    stop
    endif
    mf_ham(i,j)=hr+im*hi 
    do i2=1,num_old_ham 
    old_mf_ham(i2,i,j)=hr+im*hi 
    enddo 
  enddo 
  
  close(37) 
  
  
  return 
  end subroutine read_mean_field
 


  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !! generate a ferromagnetic matrix as initial guess !!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine generate_ferro_mean_field
  use inputs
  use system_variables
  implicit none  
  real (kind=8) :: hr,hi 
  integer :: i,j,nel,i1,i2
  
  ! initialize 
  old_mf_ham=(0.d00,0.d00) 
  mf_ham=(0.d00,0.d00) 
  
  do i=1,norbitals
    mf_ham(i,i) = mean_field_amplitude * (-1.d00)**i
  enddo

  write(*,*) '<<Generated ferromagnetic initial guess>>'
  
  return 
  end subroutine generate_ferro_mean_field




  subroutine generate_antiferro_mean_field
  use inputs
  use system_variables
  implicit none  
  real (kind=8) :: hr,hi 
  integer :: i,j,nel,i1,i2,iat
  ! initialize 
  old_mf_ham=(0.d00,0.d00) 
  mf_ham=(0.d00,0.d00) 
  
  do i=1,norbitals
    iat = (i+1)/2 ! gives the atom index
    mf_ham(i,i) = mean_field_amplitude * (-1.d00)**(i+iat) ! spin + atom
  enddo

  write(*,*) '<<Generated collinear antiferromagnetic initial guess>>'
  
  return 
  end subroutine generate_antiferro_mean_field










  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !! generate a random matrix as initial guess !!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine generate_random_mean_field
  use inputs
  use system_variables
  implicit none  
  real (kind=8) :: hr,hi,harvest_r,harvest_i
  integer :: i,j,nel,i1,i2
  ! initialize 
  old_mf_ham=(0.d00,0.d00) 
  mf_ham=(0.d00,0.d00) 
  call random_seed ! initialize
! initialice upper triangular part 
  do i=1,norbitals
    do j=i,norbitals
      call random_number(harvest_r)
      call random_number(harvest_i)
      mf_ham(i,j) = mean_field_amplitude * (harvest_r &
             *(-1.d00)**i+im*harvest_i)
    enddo
  enddo
! and the other part of the matrix
  mf_ham = mf_ham + transpose(conjg(mf_ham))  
  write(*,*) '<<Generated random initial guess>>' 
  return 
  end subroutine generate_random_mean_field





  subroutine generate_random_xz_mean_field
  use inputs
  use system_variables
  implicit none  
  real (kind=8) :: hr,hi,harvest_r,harvest_i
  integer :: i,j,nel,i1,i2
  ! initialize 
  old_mf_ham=(0.d00,0.d00) 
  mf_ham=(0.d00,0.d00) 
  call random_seed ! initialize
! initialice upper triangular part 
  do i=1,norbitals
    do j=i,norbitals
      call random_number(harvest_r)
      call random_number(harvest_i)
      mf_ham(i,j) = mean_field_amplitude * (harvest_r*(-1.d00)**i)
    enddo
  enddo
! and the other part of the matrix
  mf_ham = mf_ham + transpose(conjg(mf_ham))  
  write(*,*) '<<Generated random initial guess>>' 
  return 
  end subroutine generate_random_xz_mean_field










 


end module
