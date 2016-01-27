subroutine setup_mixing
use inputs
use system_variables
use mean_field_variables
implicit none

integer :: i ! counter
! equal mixing of old hamiltonians

coef = (1.d00-mix_coef)/dble(num_old_ham)


! the newest the better
do i = 1, num_old_ham
  coef(i) = dble(i)
enddo
coef = coef/sum(coef)  ! normalize the sum of the coefficients
coef = coef - mix_coef/dble(num_old_ham) ! substract a part to each coeffcient

! write the coefficients used in a file
open(45,file="mixing_coefficients")
  write(45,*)  mix_coef, "! new hamiltonian"
  do i = 1, num_old_ham
    write(45,*) coef(i)
  enddo
close(45)

end subroutine
