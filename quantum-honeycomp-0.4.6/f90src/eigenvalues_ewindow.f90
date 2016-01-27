subroutine eigenvalues_ewindow(ndim,matrix,nev_found,evals,&
             emin,emax,num_offdiagonals)
implicit none
!dimension of the matrixrix
integer ndim
!onsite energy and hoppings        
!energies
real (kind=8) evals(ndim)
!number of energies found
integer nev_found
!minimun and maximun values of energy
real (kind=8) emin,emax
!leading dimension of ab
integer ldab
!band matrixrix
complex (kind=8) ab(num_offdiagonals+1,ndim),work(ndim)
!matrixrix to diagonalize
complex (kind=8) matrix(ndim,ndim)
!absolute error tolerance
real (kind=8) abstol
!more variables...
integer info,ldq,ldz
integer ifail(ndim) 
complex (kind=8) evecs(ndim,ndim),q(ndim,ndim)
real (kind=8) rwork(7*ndim)
integer iwork(5*ndim)          
integer il,iu       
!counters
integer c1,c2,c3
integer :: num_offdiagonals

ldq=ndim
ldz=ndim
abstol=1.d-06


!calculate ab
ab=0.d00
do c1=1,num_offdiagonals+1
  do c2=1,ndim-c1+1
    ab(c1,c2)=matrix(c1-1+c2,c2)
  enddo
enddo



!get eigenvalues

ldab = num_offdiagonals+1

call zhbevx('N','V','L',ndim,ldab-1,ab,ldab,q, &
ldq,emin,emax,il,iu,abstol,nev_found,evals,evecs,ldz,work, &
rwork,iwork,ifail,info)


return
end subroutine
