!!! suboutines to calculate ocuupations


subroutine fermi_dirac(e,t,fermi,occ)
implicit none
real (kind=8),intent(in) :: e,fermi,t  ! energy and fermi level and temperature
real (kind=8),intent(out) :: occ ! occupation

if (t.gt.0.d00) then  ! for T>0
  occ = 1.d00/(1.d00 + exp((e-fermi)/t)) ! fermi-dirac distribution
else ! for T=0
  if (e.lt.fermi) then  
    occ=1.d00
  else 
    occ = 0.d00
  endif ! close if
endif  ! close T=0

return
end subroutine fermi_dirac

