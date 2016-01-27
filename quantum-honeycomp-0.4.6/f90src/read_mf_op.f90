! routine to read the mean field operators wanted to evaluate 
subroutine read_mf_op
use inputs
use system_variables 
use mean_field_variables 
implicit none 
! counters 
integer :: i1,i2,iii,jjj,nv
real (kind=8) :: holdr,holdi 



! open file to read 
open(47,file='mean_field_operators.in') 
read(47,*) 
read(47,*) num_AB 
write(*,*) "Read num_AB"

call initialize_mean_field ! allocate and initialize mean field variables


allocate(ab_matrices(num_AB))  ! allocate pairs of matrices


! loop over mean field terms 
do i1=1,num_AB
  ! read matrix A 
  read(47,*) 
  read(47,*) 
  read(47,*) nv
  call allocate_ab_mf(nv,ab_matrices(i1)%a)
  do i2=1,ab_matrices(i1)%a%nv
    read(47,*) iii,jjj,holdr,holdi 
    write(*,*) iii,jjj,holdr,holdi,i2 
    ab_matrices(i1)%a%i(i2) = iii
    ab_matrices(i1)%a%j(i2) = jjj
    ab_matrices(i1)%a%mij(i2) = holdr+im*holdi
  enddo 
  ! read matrix B 
  read(47,*) 
  read(47,*) 
  read(47,*) nv
  call allocate_ab_mf(nv,ab_matrices(i1)%b)
  do i2=1,ab_matrices(i1)%b%nv 
    read(47,*) iii,jjj,holdr,holdi 
    ab_matrices(i1)%b%i(i2) = iii
    ab_matrices(i1)%b%j(i2) = jjj
    ab_matrices(i1)%b%mij(i2) = holdr+im*holdi
  enddo 
  ! read couplings between matrices 
  read(47,*) 
  read(47,*) holdr,holdi 
  ab_matrices(i1)%lambda = holdr+im*holdi
  read(47,*) 
  
enddo 


close(47) 


return 
end subroutine 


subroutine allocate_ab_mf(n,mat)
  use sparse
  implicit none
  integer, intent(in) :: n ! dimension
  type(spmatrix), intent(out) :: mat ! sparse matrix
  integer :: ii
  if (allocated(mat%i)) deallocate(mat%i)
  if (allocated(mat%j)) deallocate(mat%j)
  if (allocated(mat%mij)) deallocate(mat%mij)

  allocate(mat%i(n))
  allocate(mat%j(n))
  allocate(mat%mij(n))

  mat%nv = n
  return
end subroutine 



  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!! generate mean field operators for hubbard !!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine mean_field_operators_hubbard
use inputs
use system_variables 
use mean_field_variables 
implicit none 
! counters 
integer :: i1,i2,iii,jjj,nv,iat,ia,ja,ib,jb
real (kind=8) :: holdr,holdi,lambda 


! number of operators is equal to the number of orbitals
num_AB = norbitals

call initialize_mean_field ! allocate and initialize mean field variables

allocate(ab_matrices(num_AB))  ! allocate pairs of matrices

write(*,*) 'Generating',num_AB,'mean field Hubbard matrices'

nv = 1  ! number of non vanishing elements

! non vanishin element of the matrix
holdr = 1.d00
holdi = 0.d00

do i1=1,num_AB
  iat = mod(i1,(num_AB/2))    ! atom index
  if (iat==0) iat = num_AB/2  ! for the highest value
  call allocate_ab_mf(nv,ab_matrices(i1)%a)
  call allocate_ab_mf(nv,ab_matrices(i1)%b)
  do i2=1,nv
    if (i1.le.num_AB/2) then   !!! density density pair !!!
      ia = 2*iat - 1   ! element   ia=1
      ja = 2*iat - 1   ! element   ja=1
      ib = 2*iat   ! element   ib=2
      jb = 2*iat   ! element   jb=2
      lambda = hubbard_scf ! assume positive hubbard
    else    !!! exchange exchange pair !!
      ia = 2*iat - 1   ! element   ia=1
      ja = 2*iat       ! element   ja=2
      ib = 2*iat    ! element   ib=2
      jb = 2*iat - 1      ! element   jb=1
      lambda =  - hubbard_scf ! in the exchange there is sign change
    endif
    !! A matrix !!
    ab_matrices(i1)%a%i(i2) = ia
    ab_matrices(i1)%a%j(i2) = ja
    ab_matrices(i1)%a%mij(i2) = holdr+im*holdi
    !! B matrix !!
    ab_matrices(i1)%b%i(i2) = ib
    ab_matrices(i1)%b%j(i2) = jb
    ab_matrices(i1)%b%mij(i2) = holdr+im*holdi
    !! coupling !!
    ab_matrices(i1)%lambda = lambda
  enddo   ! close loop over elements
enddo ! close loop over matrices
return 
end subroutine mean_field_operators_hubbard






!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!! generate mean field operators for hubbard !!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine mean_field_operators_collinear_hubbard
use inputs
use system_variables 
use mean_field_variables 
implicit none 
! counters 
integer :: i1,i2,iii,jjj,nv,iat,ia,ja,ib,jb
real (kind=8) :: holdr,holdi,lambda 


! number of operators is equal to the number of orbitals
num_AB = norbitals/2

call initialize_mean_field ! allocate and initialize mean field variables

allocate(ab_matrices(num_AB))  ! allocate pairs of matrices

write(*,*) 'Generating',num_AB,'mean field collineal Hubbard matrices'

nv = 1  ! number of non vanishing elements

! non vanishin element of the matrix
holdr = 1.d00
holdi = 0.d00

do i1=1,num_AB
  iat = i1    ! atom index
  call allocate_ab_mf(nv,ab_matrices(i1)%a)
  call allocate_ab_mf(nv,ab_matrices(i1)%b)
  do i2=1,nv
    ia = 2*iat - 1   ! element   ia=1
    ja = 2*iat - 1   ! element   ja=1
    ib = 2*iat   ! element   ib=2
    jb = 2*iat   ! element   jb=2
    lambda = hubbard_scf ! assume positive hubbard
    !! A matrix !!
    ab_matrices(i1)%a%i(i2) = ia
    ab_matrices(i1)%a%j(i2) = ja
    ab_matrices(i1)%a%mij(i2) = holdr+im*holdi
    !! B matrix !!
    ab_matrices(i1)%b%i(i2) = ib
    ab_matrices(i1)%b%j(i2) = jb
    ab_matrices(i1)%b%mij(i2) = holdr+im*holdi
    !! coupling !!
    ab_matrices(i1)%lambda = lambda
  enddo   ! close loop over elements
enddo ! close loop over matrices
return 
end subroutine mean_field_operators_collinear_hubbard








