! module to write results of mean field calculation
module save_mf_routines


contains
  ! routine to write the full final hamiltonian
   subroutine write_hamiltonian
   use inputs
   use system_variables
   use mean_field_variables
   complex (kind=8) full_onsite(norbitals,norbitals)
   integer :: i1,i2,nel,ihop
   real (kind=8) :: hold
   character (len=80) :: namehop

   full_onsite=onsite_0+mf_ham

   ! shift the hamiltonian to have the fermi level in zero
   if (shift_to_zero) then
     do i1=1,norbitals
       full_onsite(i1,i1) = full_onsite(i1,i1)-fermi
     enddo
     write(*,*) "Shifted fermi energy to zero"
     write(*,*) 
   endif
   
   open(41,file=hamiltonian_file)    ! file where to write the final hamiltonian
   
   write(41,*) 'DIMENSIONALITY_OF_THE_SYSTEM'
   write(41,*) system_dimension
   write(41,*)
   write(41,*)
   write(41,*) 'DIMENSION_OF_THE_HAMILTONIAN'
   write(41,*) norbitals
   write(41,*)
   write(41,*)
   ! write onsite part
   write(41,*) '     ONSITE_MATRIX'
   do i1=1,norbitals
       do i2=1,norbitals
           hold=dble(full_onsite(i1,i2)*conjg(full_onsite(i1,i2)))
           if (hold > 1.d-08) then
               write(41,*) i1,i2,dble(full_onsite(i1,i2)), &
               aimag(full_onsite(i1,i2))
           endif
       enddo
   enddo
   write(41,*) 
   write(41,*) 
   do ihop = 1, hamiltonian%num_hoppings 
!     write(*,*) hamiltonian%num_hoppings
     write(41,*)
!     write(*,*) hamiltonian%directional_hoppings(ihop)%label
     write(41,*)
     write(41,*)   
     ! write hopping part
     write(41,*) hamiltonian%directional_hoppings(ihop)%label ! write the name
     ! store in temporal array
     hopping = hamiltonian%directional_hoppings(ihop)%hopping
     do i1=1,norbitals
         do i2=1,norbitals
             hold=dble(hopping(i1,i2)*conjg(hopping(i1,i2)))
             if (hold > 1.d-08) then
                 write(41,*) i1,i2,dble(hopping(i1,i2)),aimag(hopping(i1,i2))
             endif
         enddo
     enddo
     write(41,*)
   enddo


   !!!!! write POSITIONS !!
   if (hamiltonian%has_positions) then ! if has positions
     write(41,*) 'POSITIONS'
     do i1=1,hamiltonian%norbitals ! go only for up channel
       write(41,*) hamiltonian%positions(i1,:)  
     enddo
   endif ! close positions if

   !!!!!!! write lattice !!
   if (hamiltonian%has_lattice) then
     write(41,*) 
     write(41,*) 
     write(41,*) 
     write(41,*) 'LATTICE_VECTORS'
     do i1=1,system_dimension
       write(41,*) hamiltonian%lattice_vectors(i1,:) ! read each vector
     enddo
   endif  ! close lattice if





   close(41)  ! close hamiltonian.in
   end subroutine write_hamiltonian 
   

   ! write mean field part
   subroutine write_mean_field
   use inputs
   use system_variables
   use mean_field_variables
   integer :: i1,i2,nel
   real (kind=8) :: hold
   open(91,file='mean_field.in')
   ! number of elements
   nel=0
   do i1=1,norbitals
       do i2=1,norbitals
           hold=dble(mf_ham(i1,i2)*conjg(mf_ham(i1,i2)))
           if (hold > 1.d-08) then
               nel=nel+1
           endif
       enddo
   enddo
   
   write(91,*) 'number of nonvanishing elements of H_mf'
   write(91,*) nel
   write(91,*) 'i   j   real imaginary'
   do i1=1,norbitals
       do i2=1,norbitals
           hold=dble(mf_ham(i1,i2)*conjg(mf_ham(i1,i2)))
           if (hold > 1.d-08) then
               write(91,*) i1,i2,dble(mf_ham(i1,i2)),aimag(mf_ham(i1,i2))
           endif
       enddo
   enddo
   
   close(91)  ! close hamil_mf.in
   
   return
   end subroutine write_mean_field
end module

