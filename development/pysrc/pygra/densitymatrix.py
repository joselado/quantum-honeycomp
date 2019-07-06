from __future__ import print_function, division
import numpy as np

try:
  from . import density_matrixf90
  use_fortran = True
except:
  use_fortran = False
  print("Fortran routines not working in densitymatrix.py")
 

def full_dm(h,use_fortran=True,nk=10,delta=1e-2):
  if h.dimensionality == 0: fac = 1.
  elif h.dimensionality == 1: fac = 1./nk
  elif h.dimensionality == 2: fac = 1./nk**2
  elif h.dimensionality == 3: fac = 1./nk**3
  else: raise
  es,vs = h.eigenvectors() # get eigenvectors
  if use_fortran:
    dm = density_matrixf90.density_matrix(np.array(es),np.array(vs),delta)
    return dm*fac
  else:
    return np.matrix(full_dm_python(h.intra.shape[0],es,np.array(vs)))*fac # call hte function



def full_dm_python(n,es,vs):
  """Calculate the density matrix"""
#  dm = [[0. for i in range(n)] for j in range(n)] # zero matrix
  dm = np.zeros((n,n)) +0j
  for ie in range(len(es)): # loop
    if es[ie]<0.: # if below Fermi energy
      for i in range(n):
        for j in range(n): 
          dm[i,j] += vs[ie][i].conjugate()*vs[ie][j] # add contribution
  return dm


def restricted_dm(h,use_fortran=True,mode="KPM",pairs=[],
                   scale=10.0,npol=400,ne=None):
  """Calculate certain elements of the density matrix"""
  if h.dimensionality != 0 : raise
  if mode=="full": # full inversion and then select
    dm = full_dm(h,use_fortran=use_fortran) # Full DM
    outm = np.array([dm[j,i] for (i,j) in pairs]) # get the desired ones
    return outm # return elements
  elif mode=="KPM": # use Kernel polynomial method
    if ne is None: ne = npol*4
    from . import kpm
    xin = np.linspace(-.99*scale,0.0,ne) # input x array
    out = np.zeros(len(pairs),dtype=np.complex)
    ii = 0
    for (i,j) in pairs: # loop over inputs
      (x,y) = kpm.dm_ij_energy(h.intra,i=i,j=j,scale=scale,npol=npol,
                      ne=ne,x=xin)
      out[ii] = np.trapz(y,x=x)/np.pi # pi is here so it normalizes to 0.5
      ii += 1
    return out
  else: raise
       



