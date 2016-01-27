! this functions calculates the surface green function by
!using the dyson equation
subroutine dyson(onsite,hopping,energy,ite,eps,mixing,&
green_guess,max_error,green,norbitals,num_redo)
implicit none
integer, intent(in) :: norbitals ! dimension of the hamiltonian
complex (kind=8), intent(in) :: onsite(norbitals,norbitals)
complex (kind=8), intent(in) :: hopping(norbitals,norbitals)
  ! onsite and hopping
real (kind=8), intent(in) :: energy ! energy at which compute the green function
real (kind=8), intent(in), optional :: eps ! infinitesimal part
real (kind=8), intent(in), optional :: mixing ! infinitesimal part
real (kind=8) :: mix ! infinitesimal part
complex (kind=8) :: cenergy ! energy at which compute the green function
complex (kind=8) :: store_green(norbitals,norbitals)
! output green function
complex (kind=8), intent(out) :: green(norbitals,norbitals)
! initial guess for the green function
complex (kind=8), intent(in), optional :: green_guess(norbitals,norbitals)
! maximun error
real (kind=8), intent(in), optional :: max_error
logical :: redo_ite
complex (kind=8)  :: id(norbitals,norbitals)
! temporal matrix
complex (kind=8) :: tmpm(norbitals,norbitals)
! selfenergy
complex (kind=8) :: selfenergy(norbitals,norbitals)
! dagger of hopping
complex (kind=8) :: hoppingH(norbitals,norbitals)
!error in iterative procedure
real (kind=8) :: error
complex (kind=8) :: cerrori
integer :: ite  ! number of iterations
integer :: i,j
integer, intent(out) :: num_redo

! initial guess for the green function
green = green_guess

! setup complex energy
cenergy = energy + (0.d00,1.d00)*eps

! put mixing
mix = mixing

! initialice identity
id = -onsite ! add onsite
do i=1,norbitals
  id(i,i) = id(i,i) + cenergy 
enddo

! calculate dagger of hopping
hoppingH = transpose(conjg(hopping))

redo_ite = .true.  ! do iterations no matter what
num_redo = 0
do while (redo_ite)  !! redo the iteration cycle 
  ! loop over iterations
  do i = 1, ite
    ! denominator
    call multiply(hopping,green,tmpm,norbitals)
    call multiply(tmpm,hoppingH,selfenergy,norbitals)
!    store_green = id - &
!    matmul(matmul(hopping,green),transpose(conjg(hopping)))
    store_green = id - selfenergy 
    call inverse(store_green,norbitals)
    green = store_green*mix +(1.d00-mix)*green
  enddo


  store_green = id -  &
  matmul(matmul(hopping,green),transpose(conjg(hopping)))
  call inverse(store_green,norbitals)


  
  
  error = 0.d00  ! initialize error accumulator
  
  do i=1,norbitals
  !  do j=1,norbitals
      cerrori = green(i,i)-store_green(i,i)
      error = error + real(cerrori*conjg(cerrori)) 
  !  enddo
  enddo
  error = sqrt(error/norbitals)

!  write(*,*) error,num_redo
    if (error.lt.max_error) then  ! max_error check
    redo_ite = .false.
    endif
  num_redo = num_redo + 1  ! increasde the number of redos
enddo  


return
end subroutine


! multiply matrices, assume same dimension
subroutine multiply(A,B,C,d)
implicit none
integer :: n,k,m
integer :: d
complex (kind=8) :: A(d,d),B(d,d),C(d,d)
complex (kind=8) :: alpha,beta


m = d
k = d
n = d

alpha=(1.d00,0.d00)
beta=(0.d00,0.d00)

call ZGEMM('n','n',m,n,k,alpha,A,m,B,k,beta,C,m)

! C is the output

return
end subroutine multiply









! computes the inverse of the matrix
subroutine inverse(A,n)
implicit none
! dimension
integer :: n
! matrix on imput inverse on output
complex (kind=8) :: A(n,n)
integer :: ipiv(n),info,lwork
complex (kind=8) :: work(1)

call ZGETRF( n, n, A, n, ipiv, info )

lwork=-1

call ZGETRI( n, A, n, ipiv, work, lwork, info )

lwork=n

call inv_aux(A,n,lwork,ipiv)


return
end subroutine

! inverse
subroutine inv_aux(A,n,lwork,ipiv)
implicit none
! dimension
integer :: n
!c matrix on imput inverse on output
complex (kind=8) A(n,n)
integer :: ipiv(n),info,lwork
complex (kind=8) :: work(lwork)

call ZGETRI( n, A, n, ipiv, work, lwork, info )
if (info.ne.0) then
write(*,*) 'Error in inverse',info,n
endif

return
end subroutine

