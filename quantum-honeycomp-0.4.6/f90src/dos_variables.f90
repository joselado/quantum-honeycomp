! varibles for dos calculation
module dos_variables
use inputs
use io_files
use system_variables
use expectation_values_variables
implicit none
real (kind=8), allocatable, save :: value_dos(:,:)  ! array for the value of DOS
real (kind=8), allocatable, save :: energies_dos(:)  ! array for the energies of DOS

contains
  ! subroutine to initialize the previous arrays
  subroutine initialize_dos
  implicit none
  integer :: i1
  call setup_dos_parameters  ! setup the energy spacing and that stuff
  ! deallocate if allocated
  if (allocated(value_dos)) deallocate(value_dos)
  if (allocated(energies_dos)) deallocate(energies_dos)
  ! and allocate space
  allocate(value_dos(num_operators+1,num_dos))
  allocate(energies_dos(num_dos))
  
  write(*,*) "Allocated DOS",num_dos,estep_dos
  write(*,*) 

  ! and initialize arrays
  value_dos = 0.d00
  do i1=1,num_dos
    energies_dos(i1) = emin_dos + dble(i1)*estep_dos
  enddo


  return
  end subroutine initialize_dos


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


! subroutine to clean the dos
  subroutine clean_dos
  implicit none
  integer :: i1
  ! deallocate DOS
  if (allocated(value_dos)) deallocate(value_dos)
  if (allocated(energies_dos)) deallocate(energies_dos)
  return
  end subroutine clean_dos

! subroutine to decide energy interval of DOS
  subroutine setup_dos_parameters
  implicit none

  !!!! First choose the energy window !!!!
  if  (.not. use_ewindow_dos) then ! read the energy window from the eigenvalues
    if (.not. allocated(energies)) then
      write(*,*) "Error in DOS, eigenvalues not calculated"
      stop
    endif
    emin_dos = minval(energies,num_wf)
    emax_dos = maxval(energies,num_wf)
    write(*,*) "Energy window in DOS taken from the whole spectra"
  endif
  if (use_ewindow_dos) then  ! take input ewindow
    write(*,*) "Energy window in DOS taken from input"
  endif
  ! print the energy window 
  write(*,*) "Emin_DOS=",emin_dos,"Emax_DOS=",emax_dos
  write(*,*) 

  !!! now get the number of intervals and energy step
  if (.not. define_num_dos) then ! if not defined calculate
    write(*,*) "Number of energy steps generated from estep_dos"
    num_dos = int((emax_dos-emin_dos)/estep_dos)
  endif
    write(*,*) num_dos, "energy steps in DOS"

  if (define_num_dos) then ! if defined calculate energy increasing
    write(*,*) "Energy step calculated from num_dos"
    estep_dos = (emax_dos -emin_dos)/dble(num_dos)
  endif
    write(*,*) estep_dos, "energy increasing in DOS"



  return
  end subroutine setup_dos_parameters

!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! subroutine to write the DOS
  subroutine write_dos
  implicit none
  integer :: i,k
  ! check that everithing has been calculated
  if (.not. (allocated(energies_dos).and.allocated(value_dos))) then
  write(*,*) "Error in write_dos, arrays not allocated"
  stop
  endif

  open(51,file=dos_file)  ! open file to write DOS
    write(51,*) '#  Energy   DOS  ',(operators(k)%name,k=1,num_operators)
  do i=1,num_dos  ! loop over energies
    write(51,*) energies_dos(i),(value_dos(k,i),k=1,num_operators+1)
  enddo
  write(*,*) "Written density of states in DOS.OUT"
  write(*,*)
  close(51)
  end subroutine write_dos

! Calculate DOS (and LDOS as well)
  subroutine calculate_dos
  implicit none
  
  real (kind=8) :: tmp,tmpl,tmpe
  real (kind=8) :: tmp_arrdos(num_dos) ! temporal array for dos
  integer :: i,j,iop
  real (kind=8) :: smearing_dos2, tmp_dos, max_edif
  complex (kind=8) :: exp_val
  
  
  
  call initialize_dos  ! initialize dos
  
   
  ! square of the smering for lorentzian scheme
  smearing_dos2 = smearing_dos**2
  
  ! setup a maximun energy difference
  max_edif = smearing_dos*2.d00
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! calculation without operators!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!$omp parallel default (shared) private (j,tmpe,tmp_dos,i,tmpl,tmp)
!$omp do
    do j=1,num_dos ! loop over energies
      tmpe=energies_dos(j)
      tmp_dos = 0.d00
      do i=1,num_wf ! loop over wavefunctions
        tmpl=energies(i)
        ! lorentzian (Green function) scheme
        tmp=tmpl-tmpe
        if (dabs(tmp).gt.max_edif) cycle  ! avoid this energy
        tmp=tmp*tmp+smearing_dos2
        tmp=1.d00/tmp
        tmp_dos=tmp_dos + tmp ! store in temporal variable
      enddo
      tmp_arrdos(j) = tmp_dos  ! save in global array
    enddo
!$omp end do
!$omp end parallel  
  
  value_dos(1,:) = tmp_arrdos(:)*smearing_dos  ! multiply for a common factor
  
  
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! calculation with operators!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  do iop=1,num_operators
    do j=1,num_dos ! loop over energies
      tmpe=energies_dos(j)
      tmp_dos = 0.d00
      do i=1,num_wf ! loop over wavefunctions
        tmpl=energies(i)
        ! lorentzian (Green function) scheme
        tmp=tmpl-tmpe
        if (dabs(tmp).gt.max_edif) cycle  ! avoid this energy
        tmp=tmp*tmp+smearing_dos2
        tmp=1.d00/tmp
        !!!! multiply by the expectation value !!!
        call sparse_vmw(wf(i,:),wf(i,:),operators(iop),norbitals,exp_val)
        tmp = tmp * real(exp_val)
        tmp_dos=tmp_dos + tmp ! store in temporal variable
      enddo
      tmp_arrdos(j) = tmp_dos  ! save in global array
    enddo
  value_dos(iop+1,:) = tmp_arrdos(:)*smearing_dos  ! multiply for a common factor
  enddo
  return
  end subroutine  calculate_dos



end module




