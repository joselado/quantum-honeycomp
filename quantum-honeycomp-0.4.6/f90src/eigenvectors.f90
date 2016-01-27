subroutine eigenvectors(ndim,matrix,evals,evec)
implicit none
! number of atoms
integer :: ndim
! eigenvalues
real (kind=8) :: evals(ndim)
! complex matrix to diagonalize
complex (kind=8) matrix(ndim,ndim)
complex (kind=8) work(ndim*ndim)
! more variables...
integer :: info,ldq,ldz
complex (kind=8) :: evec(ndim,ndim),q(ndim,ndim)
real (kind=8) :: rwork(3*ndim-2)
integer :: iwork(5*ndim)

! all the eigenvalues
    call ZHEEV('V', 'U', ndim, matrix, ndim, evals, &
    work, ndim*ndim, rwork,info )

 
! output matrix is eigenvectors
    evec = transpose(matrix)  ! transpose it in order to be wf(1,:) an eigen

    return
    end subroutine
