! routine to read the kpoints 
subroutine get_klist(klist,nik) 
implicit none 
integer :: nik,i 
real*8 :: klist(nik) 
logical :: ex 

! check if the file exists 
inquire(file="klist_scf.in", exist=ex ) 




! if exists read kpoints 
if (ex) then 
open(71,file='klist_scf.in') 
read(71,*) 
read(71,*) 
read(71,*) 
do i=1,nik 
read(71,*) klist(i) 
enddo 
endif 
! if does not exist, create klist 
if ( .NOT. ex) then 
do i=1,nik 
klist(i)=dble(i)/dble(nik) 
enddo 
endif 


return 
end subroutine 
