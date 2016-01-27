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

