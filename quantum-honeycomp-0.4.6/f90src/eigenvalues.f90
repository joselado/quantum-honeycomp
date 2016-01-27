subroutine eigenvalues(ndim,matrix,evals)
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
complex (kind=8) :: q(ndim,ndim)
real (kind=8) :: rwork(3*ndim-2)
integer :: iwork(5*ndim)

! all the eigenvalues
    call ZHEEV('N', 'U', ndim, matrix, ndim, evals, &
    work, ndim*ndim, rwork,info )

 

    return
    end subroutine
