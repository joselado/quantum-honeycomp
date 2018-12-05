subroutine density_matrix(es,vs,ne,dimv,dm)
implicit none
integer, intent(in) :: dimv,ne ! dimensions
real (kind=8), intent(in) :: es(ne) ! vectors
complex (kind=8), intent(in) :: vs(ne,dimv) ! vectors
complex (kind=8), intent(out) :: dm(dimv,dimv) ! output matrix
complex (kind=8) :: v(dimv) ! vector
integer :: i,j,k

dm(:,:) = (0.d00,0.d00) ! initialize


do k=1,ne
  if (es(k).lt.0.d00) then
    v(:) = vs(k,:) ! store
    do i=1,dimv
      do j=1,dimv
        dm(i,j) = dm(i,j) + conjg(v(i))*v(j)
      enddo
    enddo
  endif
enddo


return
end subroutine
