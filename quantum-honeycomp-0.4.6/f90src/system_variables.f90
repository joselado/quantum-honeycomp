! this module stores the system variables as hamiltonian, mean field operators
! wavefunctions,etc

module system_variables
implicit none

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! Hamiltonian variables !!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
integer :: system_dimension=0 ! dimensionality of the system
integer :: norbitals ! number of orbitals (dimension)
complex (kind=8), allocatable ::  onsite(:,:),hopping(:,:)  ! hamiltonian
complex (kind=8), allocatable :: wf(:,:)  ! all eigenfunctions
complex (kind=8), allocatable :: wf_k(:,:)  ! eigenfunctions at k
real (kind=8), allocatable :: energies(:) ! eigenvalues of the hamiltonian
real (kind=8), allocatable :: energies_k(:) ! eigenvalues at k
integer :: num_wf ! number of wavefunctions
integer :: num_wf_k
real (kind=8) :: fermi=0.d00



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! Diagonalization variables !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
logical :: use_ewindow=.false.   !! Use ewindow in diagonalization 
real (kind=8) :: emin=-1.d00,emax=1.d00 ! energy window
integer :: nev_found  ! eigenvalues found in diagonalization



!!!!!!!!!!!!!!!!!!!!!!!
!! Global parameters !!
!!!!!!!!!!!!!!!!!!!!!!!
complex (kind=8) :: im=(0.d00,1.d00)
real (kind=8) :: pi=acos(-1.d00)


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! define class for the hamiltonian !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
type hopping_type
  real (kind=8), allocatable :: direction(:)  ! direction of the hopping
  complex (kind=8), allocatable :: hopping(:,:)  ! matrix of the hopping
  character (len=40):: label ! name of the hopping
endtype hopping_type

type hamiltonian_type
  ! type for the directional hoppings
  logical :: has_spinpol=.true.   ! know if it has spin polarization
  logical :: has_lattice=.false.  ! know if it has lattice vectors
  logical :: has_positions=.false.  ! know if the orbitals have positions
  real (kind=8), allocatable :: positions(:,:)  ! postions of the orbitals
  real (kind=8), allocatable :: lattice_vectors(:,:)  ! lattice vectors
  real (kind=8), allocatable :: reciprocal_vectors(:,:)  ! lattice vectors
  integer :: dimensionality=0   ! dimensionality of the hamiltonian
  integer :: norbitals=2   ! number of elements
  integer :: num_hoppings=0   ! number of hoppings
  !!! diferent matrices !!!!
  complex (kind=8), allocatable :: onsite(:,:)
  type(hopping_type), allocatable :: directional_hoppings(:)  
endtype hamiltonian_type

type(hamiltonian_type),save :: hamiltonian  ! define tha hamiltonian variable



!!!!!!!!!!!!!!!!!!!
!!!! klist type !!!
!!!!!!!!!!!!!!!!!!!

type klist_type
  integer :: nkpoints = 0
  integer, allocatable :: nev_found(:) ! number of eigenfunctions found
  real (kind=8), allocatable :: kpoint(:,:) ! array of kpoints
  real (kind=8), allocatable :: weight(:) ! weight of each kpoint
  character (len=10), allocatable :: kname(:) ! array of kpoints
endtype klist_type

type(klist_type),save :: klist    ! type for the klist



contains
! initialize the hamiltonian
  subroutine initialize_hamiltonian
  use inputs  ! needed for the number of kpoints
  implicit none
  integer :: ihop
  ! deallocate if neccesary
  if (allocated(onsite)) deallocate(onsite)
  if (allocated(hopping)) deallocate(hopping)
  ! allocate onsite and hopping
  allocate(onsite(norbitals,norbitals))  ! allocate onsite matrix
!  write(*,*) "Allocated onsite matrix"
  onsite =(0.d00,0.d00)  ! initialize to 0

! allocate hopping matrix for 1d system
  if (system_dimension>0)  then 
    allocate(hopping(norbitals,norbitals)) 
!    write(*,*) "Allocated hopping matrix"
    hopping =(0.d00,0.d00)  ! initialize to 0
  endif


  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !! object version !!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  ! onsite matrix
  hamiltonian%dimensionality = system_dimension ! store dimensionality
  hamiltonian%norbitals = norbitals ! store size of the matrix
  hamiltonian%num_hoppings = 3**(system_dimension) - 1
  if (allocated(hamiltonian%onsite)) deallocate(hamiltonian%onsite)
  allocate(hamiltonian%onsite(norbitals,norbitals))  ! allocate onsite matrix
  if (allocated(hamiltonian%directional_hoppings)) &
           deallocate(hamiltonian%directional_hoppings)
  allocate(hamiltonian%directional_hoppings(hamiltonian%num_hoppings))
  ! allocate each hopping
  do ihop = 1, hamiltonian%num_hoppings
    allocate(hamiltonian%directional_hoppings(ihop) &
                 %hopping(norbitals,norbitals))
    allocate(hamiltonian%directional_hoppings(ihop) &
                 %direction(hamiltonian%dimensionality))
  enddo
  ! deallocate the real and reciprocal lattice
  if (allocated(hamiltonian%lattice_vectors)) & 
              deallocate(hamiltonian%lattice_vectors)
  if (allocated(hamiltonian%reciprocal_vectors)) & 
              deallocate(hamiltonian%reciprocal_vectors)


  return
  end subroutine initialize_hamiltonian



! initialize eigenfunctions and eigenvalues
  subroutine initialize_eigenstates
  use inputs  ! needed for the number of kpoints
  implicit none
  if (system_dimension.gt.0) then  ! if periodic
    num_wf = norbitals*klist%nkpoints  ! number of wavefunctions
!    write(*,*) " Allocating ",norbitals,"norbitals x",nkpoints, &
!              "nkpoints =", num_wf,"wavefunctions"
  endif ! end periodic
  if (system_dimension==0) then  ! if finite
    num_wf = norbitals  ! number of wavefunctions
!    write(*,*) " Allocating", num_wf,"wavefunctions"
  endif ! end finite



  ! if something is allocated deallocate
  if (allocated(wf))  deallocate(wf)   
  if (allocated(energies))  deallocate(energies)   
  ! allocate
  write(*,*) "Allocated",num_wf,norbitals,klist%nkpoints
  allocate(wf(num_wf,norbitals))  ! allocate wavefunctions
  allocate(energies(num_wf))  ! allocate eigenvalues
  ! initialice everything to zero
  wf =(0.d00,0.d00)
  energies = 0.d00
  return
  end subroutine initialize_eigenstates


! initialize eigenfunctions and eigenvalues
  subroutine initialize_eigenvalues
  use inputs  ! needed for the number of kpoints
  implicit none
  num_wf = norbitals*klist%nkpoints  ! number of wavefunctions
!  write(*,*) " Allocated ",norbitals,"norbitals x",nkpoints, &
!              "nkpoints =", num_wf,"eigenvalues"
  ! if something is allocated deallocate
  if (allocated(wf))  deallocate(wf)   
  if (allocated(energies))  deallocate(energies)   
  ! allocate
  allocate(energies(num_wf))  ! allocate eigenvalues
  ! initialice everything to zero
  energies = 0.d00
  return
  end subroutine initialize_eigenvalues



! initialize eigenvalues
  subroutine initialize_eigenvalues_k
  use inputs  ! needed for the number of kpoints
  implicit none
  ! if something is allocated deallocate
  if (allocated(energies_k))  deallocate(energies_k)   
  ! allocate
  write(*,*) "alloca",hamiltonian%norbitals 
  allocate(energies_k(hamiltonian%norbitals))  ! allocate eigenvalues
  ! initialice everything to zero
  energies_k(:) = 0.d00
  return
  end subroutine initialize_eigenvalues_k



! initialize eigenfunctions and eigenvalues
  subroutine initialize_eigenvectors_k
  use inputs  ! needed for the number of kpoints
  implicit none
  ! if something is allocated deallocate
  if (allocated(energies_k))  deallocate(energies_k)   
  if (allocated(wf_k))  deallocate(wf_k)   
  ! allocate
  allocate(energies_k(hamiltonian%norbitals))  ! allocate eigenvalues
  allocate(wf_k(hamiltonian%norbitals,hamiltonian%norbitals))  ! allocate eigenvalues
  ! initialice everything to zero
  energies_k = 0.d00
  return
  end subroutine initialize_eigenvectors_k








  ! initialice a general klist
  subroutine initialize_klist
  use inputs
  implicit none
  integer ik,jk,id,kk,nk,ak
  ! recalculate the number of kpoints according to the dimensionality
  if (hamiltonian%dimensionality==0) klist%nkpoints = 1 ! one kpoint
  if (hamiltonian%dimensionality==1) klist%nkpoints = nkpoints ! 1d
  if (hamiltonian%dimensionality==2) then
        klist%nkpoints = (nint(sqrt(float(nkpoints))))**2 ! 2d
  endif
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
!  write(*,*) 'nkpoints'
!  write(*,*) nkpoints
  ! create kpoint list
  ! 1d
  if (hamiltonian%dimensionality==1) then
    do ik=1,klist%nkpoints 
      klist%kpoint(ik,1) = dble(ik)/dble(klist%nkpoints)
    enddo
  endif


  ! 2d
  if (hamiltonian%dimensionality==2) then
    ak = 1
    nk = nint(sqrt(float(klist%nkpoints)))
    write(*,*) nk
    do ik=1,nk
      do jk=1,nk
        klist%kpoint(ak,1) = dble(jk)/dble(nk) - 5.d-01
        klist%kpoint(ak,2) = dble(ik)/dble(nk) - 5.d-01
        ak = ak+1
      enddo
    enddo
  endif


  ! 3d
  if (hamiltonian%dimensionality==3) then
    ak = 1
    do ik=1,nk
      do jk=1,nk
        do kk=1,nk
          klist%kpoint(ak,1) = dble(ik)/dble(nk)
          klist%kpoint(ak,2) = dble(jk)/dble(nk)
          klist%kpoint(ak,3) = dble(kk)/dble(nk)
          ak = ak + 1
        enddo
      enddo
    enddo
  endif


  klist%nev_found = 0   ! initialice the number of eigenvalues per k
  klist%weight(:) = 1.d00  ! initialize the weight of the kpoint
  klist%weight = klist%weight/sum(klist%weight) 

  !!!! expand the klist !!!
  klist%kpoint(:,:) = klist%kpoint(:,:)*expand_bz  



  return
  end subroutine initialize_klist


  subroutine clean_waves
  if (allocated(wf))  deallocate(wf)
  if (allocated(wf_k))  deallocate(wf_k)
  if (allocated(energies))  deallocate(energies)
  if (allocated(energies_k))  deallocate(energies_k)

  return
  end subroutine clean_waves





  ! initialice a half bz klist
  subroutine create_half_bz_klist(mesh)
  use inputs
  implicit none
  integer ik,jk,id,kk,nk,ak
  character (len=30) :: mesh ! type of mesh, surface or contour
  real (kind=8) :: dk_path

  ! recalculate the number of kpoints according to the dimensionality
  if (mesh=="bulk") then
    if (hamiltonian%dimensionality==1) then
      nk = nkpoints  ! do nothing
      klist%nkpoints = nk
    else
      nk = 2*(int(nkpoints**(0.5)))
      klist%nkpoints = nk*nk/2
    endif
  elseif (mesh=="boundary") then
    if (hamiltonian%dimensionality==1) then
      klist%nkpoints = 2  ! only two points
    else
      klist%nkpoints = 4*nkpoints ! path four times larger  
      klist%nkpoints = 4*nkpoints_connection ! manual kpoints, highy unstable...
    endif
  else
    write(*,*) "Unrecognised type of kmesh"
  endif

  ! if already allocated deallocate
  if (allocated(klist%kpoint)) deallocate(klist%kpoint)
  ! allocate the kpoints
  allocate(klist%kpoint(klist%nkpoints,hamiltonian%dimensionality))
  ! allocate the name of the kpoints
  klist%kpoint(:,:) = 0.d00

  !!!!!!!!!!!!!!!!!!!!!!!!!
  ! create kpoint list
  ! 1d
  if (hamiltonian%dimensionality==1) then
    if (mesh=="bulk") then  ! bulk contribution
      do ik=1,nk 
! extra two to do only half BZ
        klist%kpoint(ik,1) = dble(ik)/dble(2*nk) 
      enddo
    elseif (mesh=="boundary") then  ! boundary contribution
        klist%kpoint(1,1) = 0.d00   ! Gamma
        klist%kpoint(2,1) = 1.d00   ! X
    endif  ! close bulk/boundary if
  endif


  ! 2d
  if (hamiltonian%dimensionality==2) then
    ak = 1
    if (mesh=="bulk") then  ! bulk contribution
      do ik=1,nk/2
        do jk=1,nk
          klist%kpoint(ak,1) = dble(jk)/dble(nk) ! full axis
          klist%kpoint(ak,2) = dble(ik)/dble(2*nk)   ! half the axis
          ak = ak+1
        enddo
      enddo
    elseif (mesh=="boundary") then  ! boundary contribution
      nk = klist%nkpoints ! number of kpoints
      do ik=1,nk/4 ! loop over points of each part
          ! first path
          klist%kpoint(ik,1) = dble(ik)/dble(nk/4) - 5.d-01 
          klist%kpoint(ik,2) = 0.d00 
          ! second path
          klist%kpoint(ik+nk/4,1) = 5.d-01 
          klist%kpoint(ik+nk/4,2) = dble(ik)/dble(2*nk/4) 
          ! third path
          klist%kpoint(ik+2*nk/4,1) = 5.d-01 - dble(ik)/dble(nk/4) 
          klist%kpoint(ik+2*nk/4,2) = 5.d-01 
          ! fourth path
          klist%kpoint(ik+3*nk/4,1) = -5.d-01
          klist%kpoint(ik+3*nk/4,2) = 5.d-01 - dble(ik)/dble(2*nk/4) 
      enddo
      dk_path = 1.d00/(nk/4)
      dk_path = 1.d-01
      ! now compress the path
      klist%kpoint(:,2) = klist%kpoint(:,2)*(1.d00-dk_path) + &
                                           dk_path/4.d00
      klist%kpoint(:,1) = klist%kpoint(:,1)*(1.d00+dk_path)
      ! deform the path according to the black magic recipe
      do ik=1,nk ! loop over points of each part
        if ((klist%kpoint(ik,1).gt.dk_path).and. &
               (klist%kpoint(ik,1).lt.(5.d-01-dk_path)).and. &
               (klist%kpoint(ik,2).lt.25.d-02)) then
          klist%kpoint(ik,2) = klist%kpoint(ik,2) - dk_path/2 
        endif
        if ((-klist%kpoint(ik,1).gt.dk_path).and. &
               (-klist%kpoint(ik,1).lt.(5.d-01-dk_path)).and. &
               (klist%kpoint(ik,2).gt.(25.d-02))) then
          klist%kpoint(ik,2) = klist%kpoint(ik,2) + dk_path/2
        endif
      enddo           
    endif  ! close bulk/boundary if
  endif


  ! 3d
  if (hamiltonian%dimensionality==3) then
    ak = 1
    do ik=1,nk
      do jk=1,nk
        do kk=1,nk
          klist%kpoint(ak,1) = dble(ik)/dble(nk) - 0.5 ! full axis
          klist%kpoint(ak,2) = dble(jk)/dble(nk) - 0.5 ! full axis
          klist%kpoint(ak,3) = dble(kk)/dble(2*nk)  ! half axis
          ak = ak + 1
        enddo
      enddo
    enddo
  endif

  return
  end subroutine create_half_bz_klist


end module
