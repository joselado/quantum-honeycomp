! variables for calculation of dielectric response
module dielectric_variables
use inputs 
use system_variables
use sparse ! sparse library
! class for a pair of operators in linear response
  type ab_pair
    type(spmatrix) :: a,b
    character(20) :: label
  endtype ab_pair


! chi at different energies   
type(ab_pair), allocatable :: ab_matrices(:)  ! pairs of matrices AB for LR
complex (kind=8), allocatable :: chi_total(:,:) ! total LR 
real (kind=8), allocatable :: eigen_chi_r(:,:) ! normal states of Real(chi)   
real (kind=8), allocatable :: chi_r(:,:) ! normal states of Real(chi)   
real (kind=8), allocatable :: eigen_chi_i(:,:) ! normal states of Imag(chi)   
real (kind=8), allocatable :: chi_i(:,:) ! normal states of Imag(chi)   
! this is not a true density matrix is something like
! n*_i m_i m*_j n_j
! in the future this should be revised....
complex (kind=8), allocatable :: den2_wf_ij(:,:) ! pseudo density matrix ij   
real (kind=8), allocatable :: ene_chi(:)  ! energies of chi 
real (kind=8), allocatable :: ene_chi_initial(:)  ! energies of chi 

real (kind=8), allocatable :: ene_chi_i(:)  ! i energies of the system 
real (kind=8), allocatable :: ene_chi_j(:)  ! j energies of the system 
complex (kind=8), allocatable :: wf_chi_i(:,:)  ! i wavefunctions of the system 
complex (kind=8), allocatable :: wf_chi_j(:,:)  ! j wavefunctions of the system 
integer :: num_wf_i, num_wf_j ! number of wavefunctions inside the interval 
integer :: num_chi_ab ! number of pair of operators for linear response
logical :: do_rpa = .false. ! perform rpa calculation
! rpa calculation, useful for looking for excitations 


!!!!!!!!!!!!!!!!!!!!!!!!
!! variables for RPA !!!
!!!!!!!!!!!!!!!!!!!!!!!!
complex (kind=8), allocatable :: chi_total_rpa(:,:)  ! total RPA chi

! qvector of the response
real (kind=8), allocatable :: qvector(:)




contains
  subroutine initialize_chi_ab  ! initialice chi_w
  implicit none
  integer :: i
  real (kind=8) :: ene_dif
  ! deallocate if already allocated
  if (allocated(chi_total))    deallocate(chi_total)   !total Chi
  if (allocated(chi_total_rpa))    deallocate(chi_total_rpa)   !total Chi
  if (allocated(chi_r))    deallocate(chi_r) !Eigenstates of Chi
  if (allocated(chi_i))    deallocate(chi_i) !Eigenstates of Chi
  if (allocated(ene_chi))      deallocate(ene_chi)     ! E    
  if (allocated(ene_chi_i))    deallocate(ene_chi_i)   !E_i
  if (allocated(ene_chi_j))    deallocate(ene_chi_j)   !E_j
  if (allocated(wf_chi_i))     deallocate(wf_chi_i)   !Wf_i
  if (allocated(wf_chi_j))     deallocate(wf_chi_j)   !Wf_j
  ! allocate new ones
  allocate(chi_total(num_ene_chi,num_chi_ab))  !total Chi
  if (do_rpa)  allocate(chi_total_rpa(num_ene_chi,num_chi_ab))  !total Chi
  allocate(chi_r(num_ene_chi,num_chi_ab))            !Eigenstates of chi
  allocate(chi_i(num_ene_chi,num_chi_ab))            !Eigenstates of chi
  allocate(ene_chi(num_ene_chi))                        ! E
  allocate(ene_chi_i(norbitals))                       !E_i
  allocate(ene_chi_j(norbitals))                       !E_j
  allocate(wf_chi_i(norbitals,norbitals))               !Wf_i
  allocate(wf_chi_j(norbitals,norbitals))               !Wf_j

!  write(*,*) "Allocated dielectric Chi(w) response"
!  write(*,*)
  ! initialize arrays
  chi_total = (0.d00,0.d00)
  ene_dif = (emax_chi-emin_chi)/dble(num_ene_chi)
  do i=1,num_ene_chi
    ene_chi(i) = emin_chi + ene_dif*dble(i)
  enddo

  ! initialize number of wavefunctions
  num_wf_i = norbitals
  num_wf_j = norbitals



  ! initialize the qvector
  if (allocated(qvector)) deallocate(qvector)
  allocate(qvector(hamiltonian%dimensionality))
  qvector(:) = 0.d00  ! initialize the q vector



  return
  end subroutine initialize_chi_ab


  subroutine clean_chi_ab  ! clean chi_w
  implicit none
  ! deallocate if already allocated
  if (allocated(chi_total))    deallocate(chi_total)   !total Chi
  if (allocated(chi_total_rpa))    deallocate(chi_total_rpa)   !total Chi
  if (allocated(chi_r))    deallocate(chi_r) !Eigenstates of Chi
  if (allocated(chi_i))    deallocate(chi_i) !Eigenstates of Chi
  if (allocated(ene_chi))      deallocate(ene_chi)     ! E    
  if (allocated(ene_chi_i))    deallocate(ene_chi_i)   !E_i
  if (allocated(ene_chi_j))    deallocate(ene_chi_j)   !E_j
  if (allocated(wf_chi_i))     deallocate(wf_chi_i)   !Wf_i
  if (allocated(wf_chi_j))     deallocate(wf_chi_j)   !Wf_j
!  write(*,*) "Deallocated dielectric Chi(w) response"
!  write(*,*)
  return
  end subroutine clean_chi_ab



  subroutine write_full_chi_ab  ! write the output
  implicit none
  integer :: i
  real (kind=8) :: chi2sigma ! factor to transform chi2sigma
  real (kind=8) :: eps(num_ene_chi,num_chi_ab) ! reflectivity
  real (kind=8) :: ref(num_ene_chi,num_chi_ab) ! reflectivity
  real (kind=8) :: nij(num_ene_chi,num_chi_ab) ! reflective index
  real (kind=8) :: kij(num_ene_chi,num_chi_ab) ! k alpha alphay
  real (kind=8) :: loss(num_ene_chi,num_chi_ab) ! k alpha alphay

  chi2sigma = 1.d00/sum(qvector*qvector)

  ! check if max_die_branches has a legal value
  chi_r = real(chi_total)
  chi_i = aimag(chi_total)
  !!!! if charge susceptibility !!!!
!  if ((dielectric_type=='local_density').or. &
!            (dielectric_type=='local_density_tpa')) then
    open(71,file="CHI_R.OUT")
    open(72,file="CHI_I.OUT")
!    open(73,file="SIGMA_R.OUT")
!    open(74,file="SIGMA_I.OUT")
!    open(75,file="ABSORPTION.OUT")
    write(71,*) "# energy  ",ab_matrices(:)%label
    write(72,*) "# energy  ",ab_matrices(:)%label
    do i=1,num_ene_chi
      write(71,*) ene_chi(i),chi_r(i,:)
      write(72,*) ene_chi(i),chi_i(i,:)
!      write(73,*) ene_chi(i),chi_i(i,:)*ene_chi(i)*chi2sigma
!      write(74,*) ene_chi(i),chi_r(i,:)*ene_chi(i)*chi2sigma
!      write(75,*) ene_chi(i),chi_i(i,:)*ene_chi(i)*chi2sigma*ene_chi(i)
    enddo
    close(71)
    close(72)
!    close(73)
!    close(74)
!    close(75)
!  endif
  ! write now the RPA result
  if (do_rpa) then
    chi_r = real(chi_total_rpa)
    chi_i = aimag(chi_total_rpa)
    open(71,file="CHI_RPA_R.OUT")
    open(72,file="CHI_RPA_I.OUT")
    write(71,*) "# energy  ",ab_matrices(:)%label
    write(72,*) "# energy  ",ab_matrices(:)%label
    do i=1,num_ene_chi
      write(71,*) ene_chi(i),chi_r(i,:)
      write(72,*) ene_chi(i),chi_i(i,:)
    enddo
    close(71)
    close(72)
  endif

  return
  end subroutine write_full_chi_ab








  ! give the occupation at a given temperature
  subroutine occupation(energy,temperature,occ)
  implicit none
  real (kind=8), intent(out) :: occ
  real (kind=8), intent(in) :: temperature,energy
  ! if below threhold
  if (energy.lt.(-temperature)) then
  occ = 1.d00
  ! if above threhold
  else if (energy.gt.(temperature)) then
  occ = 0.d00
  else 
  ! if in interval
  occ = 5.d-01 - (energy)/(2.d00*temperature)
  end if
  end subroutine occupation







  !!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!! AB linear response !!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine calculate_chi_ab
  implicit none
  ! index for energy wave i and wave j
  integer :: ie,iwf1,iwf2
  real (kind=8) :: enetmp,etmp1,etmp2 ! temporal energies
  ! temporal wavefunctions
  complex (kind=8) :: wftmp1(norbitals),wftmp2(norbitals)
  complex (kind=8) :: chitmp  ! temporal chi
  complex (kind=8) :: ieps  ! smearing of Chi
  complex (kind=8) :: vawwbv,holdab  ! matrix element of the response
  real (kind=8) :: occ_fac,occ1,occ2 ! occupation factor 
  real (kind=8) :: den_res ! denominator of the response function
  integer :: iab ! counter for pair of operators

  ieps = im*smearing_chi  ! complex infinitesimal


    do iwf1 = 1, num_wf_i ! first loop over states
      etmp1 = ene_chi_i(iwf1)
      wftmp1(:) = wf_chi_i(iwf1,:)
      do iwf2 = 1, num_wf_j  ! second loop over states
        etmp2 = ene_chi_j(iwf2)
        ! fermi energy has been put in 0
        call occupation(etmp1,temperature_chi,occ1)
        call occupation(etmp2,temperature_chi,occ2)
        occ_fac = occ2-occ1  ! occupation factor
        if (dabs(occ_fac).lt.1.d-04)  cycle  ! next iteration if too far
        ! if contribution is different from zero continue
        wftmp2(:) = wf_chi_j(iwf2,:)
!$omp parallel default (shared) private (iab,ie,enetmp,den_res,chitmp,holdab,vawwbv)
!$omp do
! loop over the different matrices of linear response
        do iab=1,num_chi_ab   
    ! calculate only at zero temperature
        ! calculate the matrix elementes <A><B>
          call sparse_vmw(wftmp1,wftmp2,ab_matrices(iab)%a,norbitals,holdab)  
          vawwbv = holdab  ! first matrix element
          call sparse_vmw(wftmp2,wftmp1,ab_matrices(iab)%b,norbitals,holdab) 
          vawwbv = vawwbv * holdab    ! second matrix element
  ! save in more accesible variables
          do ie = 1, num_ene_chi
            enetmp = ene_chi(ie)
            den_res = ((etmp1-etmp2) - enetmp)
            chitmp = occ_fac*vawwbv/(den_res+ieps)  ! add contribution
            ! ADD to the chi matrix array
            chi_total(ie,iab) = chi_total(ie,iab) + chitmp 
          enddo ! close loop over energies
        enddo ! close loop over AB
!$omp end do
!$omp end parallel
    enddo  ! close loop over energies
  enddo ! close loop over the different matrices of linear response
  return
  end subroutine calculate_chi_ab







  !!!! RPA dielectric function !!!!
  subroutine rpa_response(chi_dim)
  implicit none  
  integer :: ie ! counter for energies
  integer :: chi_dim ! dimension of the response
  complex (kind=8) :: chi0_rpa(chi_dim,chi_dim)   ! initial chi of rpa
  complex (kind=8) :: chi_rpa(chi_dim,chi_dim)   ! final chi of rpa
  complex (kind=8) :: chi_rpa2(chi_dim,chi_dim)
  complex (kind=8) :: mat_tmp(chi_dim,chi_dim)   ! holder for rpa
  complex (kind=8) :: mat_lr(chi_dim,chi_dim)   ! long range matrix
  complex (kind=8) :: det_rpa  ! determinant of the
  integer :: i,j,k
  real (kind=8) :: interaction_factor ! factor from interaction
  complex (kind=8) det
  external det

  chi_dim = int(sqrt(real(num_chi_ab))) ! dimension of the chi matrix
  if (.not.(chi_dim*chi_dim == num_chi_ab)) then
    write(*,*) "Something wrong with the matrices in RPA",chi_dim,&
                num_chi_ab
    stop
  endif 


  mat_lr = (0.d00,0.d00)


  ! create the interaction matrix 
  do i=1,chi_dim
!   mat_lr(i,i) = mat_lr(i,i) + hubbard_scf  
    mat_lr(i,i) = hubbard_scf  
  enddo

  open(57,file='DETERMINANT_RPA.OUT')
    write(57,*) '# energy   RPA determinant'

  do ie = 1, num_ene_chi
  ! convert to matrix
    mat_tmp = (0.d00,0.d00)  
    k = 1
    do i=1,chi_dim ! loop over rows of chi
      do j=1,chi_dim  ! loop over columns of chi
        chi0_rpa(i,j) = chi_total(ie,k) ! get the element
        k = k+1  ! increase the counter
      enddo
      mat_tmp(i,i) = (1.d00,0.d00)  ! initialice identity
    enddo
    ! apply dyson equation
    mat_tmp = mat_tmp - matmul(mat_lr,chi0_rpa)   ! I - U \Chi_0
    call inverse(mat_tmp,chi_dim)  ! 1/(I -U \Chi_0)
    call determinant(chi_dim,mat_tmp,det_rpa) ! determinant for looking for exitations
    det_rpa = sqrt(det_rpa*conjg(det_rpa))
    write(57,*) ene_chi(ie),dble(det_rpa)
    chi_rpa2 = matmul(mat_tmp,chi0_rpa)  ! oposite order
    chi_rpa = matmul(chi0_rpa,mat_tmp)  ! correct order

    ! check if there is diffference in the order
    chi_rpa2 = chi_rpa2-chi_rpa
    chi_rpa2 = chi_rpa2*conjg(chi_rpa2)
    if (real(sum(chi_rpa2)).gt.1.d-03) then
       write(*,*) "Warning!!!, RPA matrices do not commute",sum(chi_rpa2)
    endif

    ! apply the wrapper to the energy array
    k = 1
    do i=1,chi_dim
      do j=1,chi_dim
        chi_total_rpa(ie,k) = chi_rpa(i,j)  ! get the element
        k = k+1  ! increase the counter
      enddo
    enddo
  enddo 

  close(57)

  ! deallocate temporal arrays

  return
  end subroutine rpa_response


  subroutine write_trace_chi_ab  ! write the output
  implicit none
  integer :: i


  ! check if max_die_branches has a legal value
  chi_r = real(chi_total)
  chi_i = aimag(chi_total)
  !!!! if charge susceptibility !!!!
!  if ((dielectric_type=='local_density').or. &
!            (dielectric_type=='local_density_tpa')) then
    open(71,file="CHI_TRACE_R.OUT")
    open(72,file="CHI_TRACE_I.OUT")
    write(71,*) "# energy,      Trace(Chi_R) "
    write(72,*) "# energy,      Trace(Chi_I) "
    do i=1,num_ene_chi
      write(71,*) ene_chi(i),chi_r(i,:)
      write(72,*) ene_chi(i),chi_i(i,:)
    enddo
    close(71)
    close(72)
    close(73)
    close(74)
    close(75)
!  endif
  ! write now the RPA result
  if (do_rpa) then
    chi_r = real(chi_total_rpa)
    chi_i = aimag(chi_total_rpa)
    open(71,file="CHI_RPA_R.OUT")
    open(72,file="CHI_RPA_I.OUT")
    write(71,*) "# energy  ",ab_matrices(:)%label
    write(72,*) "# energy  ",ab_matrices(:)%label
    do i=1,num_ene_chi
      write(71,*) ene_chi(i),chi_r(i,:)
      write(72,*) ene_chi(i),chi_i(i,:)
    enddo
    close(71)
    close(72)
  endif

  return
  end subroutine write_trace_chi_ab







end module
