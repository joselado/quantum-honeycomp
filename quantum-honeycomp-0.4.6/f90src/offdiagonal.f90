! this module only has one variable, used for determining the
! number of of offdiagonals of the hamiltonian
module offdiagonal
implicit none
!integer :: num_offdiagonals


contains

  ! calculate the number of offdiagonals of a matrix
  subroutine get_number_offdiagonals(m,mdim,numoff)
  implicit none
  integer, intent(in) :: mdim ! dimension of the input matrix
  integer, intent(out) :: numoff ! number of offdiagonals
  complex (kind=8), intent(in) :: m(mdim,mdim) ! matrix 
  integer :: i,j,tmpoff 
  real (kind=8) :: cutoff ! cutoff of the elements 
  real (kind=8) :: rr,ii ! real and imaginary parts

  cutoff = 1.d-05 ! default cutoff

  numoff = 0

  ! loop over the matrix
  do i=1,mdim
    do j=1,mdim
      rr = real(m(i,j)) ! real part
      ii = aimag(m(i,j)) ! imaginary part
      if ((rr.gt.cutoff).or.(ii.gt.cutoff)) then
        tmpoff = abs(i-j) ! off diagonal number...
        if (tmpoff.gt.numoff) numoff=tmpoff ! maximum offdiagonal 
      endif
    enddo
  enddo

  return

  end subroutine get_number_offdiagonals

end module offdiagonal
