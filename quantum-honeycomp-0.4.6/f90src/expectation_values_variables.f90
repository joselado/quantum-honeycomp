module expectation_values_variables

use sparse  ! use sparse routines !!!
use inputs
use system_variables
implicit none
type(spmatrix), allocatable ::  operators(:)  ! operators
integer :: num_operators ! number of operators
real (kind=8), allocatable :: exp_val_values(:,:) ! expectation value

contains

  ! default operators
  subroutine default_operators
  implicit none
  num_operators = 0  ! initial number of operators
  if (allocated(operators)) deallocate(operators)  ! deallocate operators
  !!!! spin polarization !!!!
  if (hamiltonian%has_spinpol) then   ! spin polarized case
    allocate(operators(3))
    call sparse_sx(norbitals,operators(1))
    call sparse_sy(norbitals,operators(2))
    call sparse_sz(norbitals,operators(3))
    num_operators = num_operators + 3
  endif 
  return
  end subroutine default_operators



!!!!!!!!!!!!!!!!!
!! Sx operator !!
!!!!!!!!!!!!!!!!!

  subroutine operator_sx
  implicit none
  num_operators = 0  ! initial number of operators
  if (allocated(operators)) deallocate(operators)  ! deallocate operators
  !!!! spin polarization !!!!
  if (hamiltonian%has_spinpol) then   ! spin polarized case
    allocate(operators(1))
    call sparse_sx(norbitals,operators(1))
    num_operators = num_operators + 1
  endif
  end subroutine operator_sx



!!!!!!!!!!!!!!!!!
!! Sy operator !!
!!!!!!!!!!!!!!!!!

  subroutine operator_sy
  implicit none
  num_operators = 0  ! initial number of operators
  if (allocated(operators)) deallocate(operators)  ! deallocate operators
  !!!! spin polarization !!!!
  if (hamiltonian%has_spinpol) then   ! spin polarized case
    allocate(operators(1))
    call sparse_sy(norbitals,operators(1))
    num_operators = num_operators + 1
  endif
  end subroutine operator_sy

!!!!!!!!!!!!!!!!!
!! Sz operator !!
!!!!!!!!!!!!!!!!!

  subroutine operator_sz
  implicit none
  num_operators = 0  ! initial number of operators
  if (allocated(operators)) deallocate(operators)  ! deallocate operators
  !!!! spin polarization !!!!
  if (hamiltonian%has_spinpol) then   ! spin polarized case
    allocate(operators(1))
    call sparse_sz(norbitals,operators(1))
    num_operators = num_operators + 1
  endif
  end subroutine operator_sz



!!!! index operator !!!!
  subroutine operator_index
  implicit none
  num_operators = 0  ! initial number of operators
  if (allocated(operators)) deallocate(operators)  ! deallocate operators
    allocate(operators(1))
  !!!! spin polarization !!!!
  call sparse_index(norbitals,operators(1))
  num_operators = num_operators + 1
  end subroutine operator_index




  ! read a operator from a file
  subroutine read_expectation_operator
  implicit none
  integer i,j,iel,exp_val_num_entries,nv
  real (kind=8) :: rp,ip

  ! allocate just one operator
  if (allocated(operators)) deallocate(operators)
  allocate(operators(1))
  ! deallocate if allocated
  open(71,file="operator.in")
  read(71,*)   ! name
  read(71,*)   ! number of element string
  read(71,*) nv
  operators(1)%nv = nv
  read(71,*) ! notation string
  ! allocate arrays
  allocate(operators(1)%mij(nv))   ! values
  allocate(operators(1)%i(nv))   ! index i
  allocate(operators(1)%j(nv))   ! index j
  do iel=1,nv
    read(71,*) i,j,rp,ip
    operators(1)%i(iel) = i
    operators(1)%j(iel) = j
    operators(1)%mij(iel) = rp + (0.d00,1.d00)*ip
  enddo
  close(71)

  operators(1)%name = '  from_file  '

  return
  end subroutine



  ! routine to calculate the expectation value of an array of matrices 
  subroutine calculate_exp_val_values 
  implicit none 
  integer :: i,j,k,l,iop
  complex (kind=8) :: hold
  if (allocated(exp_val_values)) deallocate(exp_val_values)
  allocate(exp_val_values(num_operators,num_wf))  ! allocate the values
  ! calculate expectation value 
  do iop=1,num_operators
    do i=1,num_wf 
      hold = (0.d00,0.d00) 
      do l=1,operators(iop)%nv
        j=operators(iop)%i(l) 
        k=operators(iop)%j(l) 
        hold=hold+conjg(wf(i,k))*operators(iop)%mij(l)*wf(i,j) 
      enddo 
      exp_val_values(iop,i)=real(hold) 
    enddo 
  enddo
  return 
  end subroutine 

  ! routine to calculate the expectation value of an array of matrices 
  subroutine calculate_eev_kpoint(wfk,enek,num_wf,n) 
  implicit none 
  integer :: i,j,k,l,iop
  integer :: n  ! dimension
  integer :: num_wf ! number of wavefunctions
  complex (kind=8) :: hold
  complex (kind=8), intent(in) :: wfk(n,n) ! wavefunctions
  real (kind=8), intent(in) :: enek(n,n) ! wavefunctions
  if (allocated(exp_val_values)) deallocate(exp_val_values)
  allocate(exp_val_values(num_operators,num_wf))  ! allocate the values
  ! calculate expectation value 
  do iop=1,num_operators
    do i=1,num_wf 
      hold = (0.d00,0.d00) 
      do l=1,operators(iop)%nv
        j=operators(iop)%i(l) 
        k=operators(iop)%j(l) 
        hold=hold+conjg(wfk(i,k))*operators(iop)%mij(l)*wfk(i,j) 
      enddo 
      exp_val_values(iop,i)=real(hold) 
    enddo 
  enddo
  return 
  end subroutine 



  



end module
