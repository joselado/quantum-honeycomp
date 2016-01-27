subroutine diagonalize_collinear(ndim,matrix,evals,evec)
implicit none
! number of atoms
integer :: ndim
! eigenvalues
real (kind=8) :: evals(ndim)
real (kind=8) :: e1(ndim/2)
real (kind=8) :: e2(ndim/2)
complex (kind=8) matrix(ndim,ndim)
complex (kind=8) m1(ndim/2,ndim/2)
complex (kind=8) m2(ndim/2,ndim/2)
complex (kind=8) :: evec(ndim,ndim)
complex (kind=8) :: v1(ndim/2,ndim/2)
complex (kind=8) :: v2(ndim/2,ndim/2)
! complex matrix to diagonalize
integer i,j,n

m1 = (0.d00,0.d00)
m2 = (0.d00,0.d00)
evec = (0.d00,0.d00)

n = ndim/2 ! spinless orbitals

! assign spin up down sectors
do i=1,n
  do j=1,n
    m1(i,j) = matrix(2*i-1,2*j-1)
    m2(i,j) = matrix(2*i,2*j)
  enddo
enddo 


! diagonalize both sectors
call diagonalize(n,m1,e1,v1)
call diagonalize(n,m2,e2,v2)

! assign eigenvalues
do i=1,n
  evals(i) = e1(i)
  evals(i+n) = e2(i)
  do j=1,n
    evec(i,2*j-1) = v1(i,j)  ! up channel
    evec(i+n,2*j) = v2(i,j)  ! down channel
  enddo
enddo



return
end subroutine diagonalize_collinear
