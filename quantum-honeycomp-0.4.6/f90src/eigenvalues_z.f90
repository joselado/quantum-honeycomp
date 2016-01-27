subroutine eigenvalues_z(n,matrix,eigen)
implicit none
! input variables
integer :: n ! dimension of the matrix
complex (kind=8) :: matrix(n,n) ! matrix
complex (kind=8) :: eigen(n) ! eigenvalues
! internal variables
complex (kind=8) :: m(n,n) ! matrix
complex (kind=8) :: vl(1,n) ! matrix
complex (kind=8) :: vr(1,n) ! matrix
complex (kind=8) :: work(4*n)
complex (kind=8) :: vs(1,n)
real (kind=8) :: rwork(2*n)
integer :: info,ldvr,ldvl
m = matrix ! copy the matrix
ldvr = 1
ldvl = 1

write(*,*) m

call CGEEV( 'N', 'N', n, m, n, eigen, vl, ldvl, vr, ldvr, &
         work, 4*n, rwork, info )

write(*,*) info

return
end subroutine





