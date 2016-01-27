subroutine determinant(N, matin,det)
implicit none
integer :: N 
complex (kind=8):: matin(N,N)
!complex (kind=8):: matin(:,:)
complex (kind=8) :: mat(N,N)
integer :: i,j, info
integer :: ipiv(N)
complex (kind=8) ::  sgn
complex (kind=8) :: ONE,det

ONE=(1.d00,0.d00)


do i=1,N
  do j=1,N
    mat(i,j) = matin(i,j)
  enddo
enddo


ipiv = 0

call zgetrf(N, N, mat, N, ipiv, info)

det = ONE

do i = 1, N
   det = det*mat(i, i)
end do
sgn = ONE
do i = 1, N
   if(ipiv(i) /= i) then
       sgn = -sgn
   end if
end do
det = sgn*det   

end subroutine determinant
