! sorting function, modified version taken from
! http://www.cs.mtu.edu/~shene/COURSES/cs201/NOTES/chap08/sorting.f90


! --------------------------------------------------------------------
! INTEGER FUNCTION  FindMinimum():
!    This function returns the location of the minimum in the section
! between Start and End.
! --------------------------------------------------------------------

   subroutine  FindMinimum(x, Start, Siz,location)
      IMPLICIT  NONE
      integer :: Siz
      real (kind=8) :: x(Siz)
      INTEGER                :: Start, End
      real (kind=8)                            :: Minimum
      INTEGER                            :: Location
      INTEGER                            :: i

      Minimum  = abs(x(Start))		! assume the first is the min
      Location = Start			! record its position
      DO i = Start+1, Siz		! start with next elements
      ! modified for absolute values
         IF (abs(x(i)) < Minimum) THEN	!   if abs(x(i)) less than the min?
            Minimum  = abs(x(i))		!      Yes, a new minimum found
            Location = i                !      record its position
         END IF
      END DO
   END subroutine  FindMinimum

! --------------------------------------------------------------------
! SUBROUTINE  Swap():
!    This subroutine swaps the values of its two formal arguments.
! --------------------------------------------------------------------

   SUBROUTINE  Swap(a, b)
      IMPLICIT  NONE
      real (kind=8), INTENT(INOUT) :: a, b
      real (kind=8)                :: Temp

      Temp = a
      a    = b
      b    = Temp
   END SUBROUTINE  Swap


subroutine reverse_array(array,n)
! reverse the order of an array
integer :: n,i
real (kind=8) :: array(n),tarr(n)

tarr = array

do i=1,n
  array(i) = tarr(n+1-i)
enddo


return
end subroutine



! --------------------------------------------------------------------
!    This subroutine receives an array x() and sorts it into ascending
! order.
! --------------------------------------------------------------------

SUBROUTINE  sort_absolute(x, Siz)
IMPLICIT  NONE
INTEGER                   :: Siz
real (kind=8) :: x(Siz)
INTEGER                               :: i
INTEGER                               :: Location
DO i = 1, Siz-1			! except for the last
   call FindMinimum(x, i, Siz,Location)	! find min from this to last
   CALL  Swap(x(i), x(Location))	! swap this and the minimum
END DO
call reverse_array(x,Siz)
END SUBROUTINE  sort_absolute
