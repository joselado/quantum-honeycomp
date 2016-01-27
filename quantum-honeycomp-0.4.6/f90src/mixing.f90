! subroutine to mix hmf and calculate error 
subroutine mixing
use inputs 
use system_variables
use mean_field_variables
implicit none 
complex (kind=8) :: mfold(norbitals,norbitals) 
complex (kind=8) :: hold
real (kind=8) :: cnew,tra1,tra2 
integer :: i1,i2,i3,ierror 

i1=1 



!
!! create new mean field hamiltonian by mixing old solutions 
mfold=(0.d00,0.d00) 
do i3=1,num_old_ham 
      mfold(:,:)=mfold(:,:)+old_mf_ham(i3,:,:)*coef(i3) 
enddo
mf_ham=mf_ham*mix_coef+mfold 
!

! calculate maximun desviation in the mean field hamiltonian 
scf_err=0.d00 
do i1=1,norbitals 
  do i2=1,norbitals 
  hold=old_mf_ham(num_old_ham,i1,i2)-mf_ham(i1,i2) 
!  hold=old_mf_ham(num_old_ham,i1,i2)
  hold=hold*conjg(hold) 
    if (dble(hold) > scf_err) then 
    scf_err=dble(hold) 
    endif 
  enddo 
enddo 

scf_err=sqrt(scf_err) 

open(21,file="ERROR.OUT",ACCESS = 'APPEND')
write(21,*) scf_err
close(21)

write(*,*) "Selfconsistency error = ", scf_err

! update old hamiltonians 
do i3=2,num_old_ham 
  old_mf_ham(i3-1,:,:)=old_mf_ham(i3,:,:) 
enddo 
old_mf_ham(num_old_ham,:,:)=mf_ham(:,:) ! update the most recent one


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! update the mixing coefficient !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


!mix_coef = 1 - scf_err/(1.d-02+scf_err)
!if (mix_coef.lt.0.7) mix_coef=7.d-01
!write(*,*) 'mix_coef',mix_coef




return 
end subroutine 
