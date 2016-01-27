! calculate the expansion in chevychev polynomials
subroutine  get_momentsf90(is,js,ms,v,numc,numij,numm,mus)
implicit none
! inputs
integer, intent(in) :: numc ! number of coefficients
integer, intent(in) :: numij ! number of entries of the matrix
integer, intent(in) :: numm ! dimension of the vector
integer, intent(in) :: is(numij),js(numij) ! indexes of the matrix
complex (kind=8), intent(in) :: ms(numij) ! entries of the matrix
complex (kind=8), intent(in) :: v(numm) ! input vector

! output
complex (kind=8), intent(out) :: mus(numc*2) ! coefficients of the expansion

! working varibles
complex (kind=8) :: a(numm)
complex (kind=8) :: am(numm)
complex (kind=8) :: ap(numm)
complex (kind=8) :: bk,bk1
complex (kind=8) :: mu0,mu1
integer i,j,k,l ! counter
! now go to work

! this routine is just the transcription in fortran of the
! routine in kpm.py


! first iteration
a(:) = 0.d00 ! initialize
am(:) = v(:) ! initialize
do k=1,numij  ! loop over non vanishing elements
  i = is(k) ! get index
  j = js(k) ! get index
  a(i) = a(i) + ms(k)*v(j)  ! perform product
enddo

bk = sum(conjg(v(:))*v(:)) ! scalar product
bk1 = sum(conjg(a(:))*v(:)) ! scalar product
! store first elements
mus(1) = bk
mus(2) = bk1

! now do the rest
do l=2,numc
  ap(:) = 0.d00 ! initialize
  do k=1,numij  ! loop over non vanishing elements
    i = is(k) ! get index
    j = js(k) ! get index
    ap(i) = ap(i) + 2.d00*ms(k)*a(j)  ! perform product
  enddo
  ap(:) = ap(:) - am(:) ! substract am, recursion relation
  bk = sum(conjg(a(:))*a(:)) ! scalar product
  bk1 = sum(conjg(ap(:))*a(:)) ! scalar product
  mus(2*l-1) = 2.d00*bk
  mus(2*l) = 2.d00*bk1
  am(:) = a(:) ! next iteration
  a(:) = ap(:) ! next iteration
enddo

! now substract first term
mu0 = mus(1)
mu1 = mus(2)
do l=2,numc
  mus(2*l-1) = mus(2*l-1) - mu0
  mus(2*l) = mus(2*l) - mu1
enddo



return
end subroutine get_momentsf90
