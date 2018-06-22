from scipy.linalg import eigvalsh
import numpy as np


def gsenergy(h,nk=10):
  """Calculate the ground state energy"""
  if h.dimensionality!=2: raise  # continue if two dimensional
  hk_gen = h.get_hk_gen() # gets the function to generate h(k)
  kxs = np.linspace(0.,1.,nk)  # generate kx
  kys = np.linspace(0.,1.,nk)  # generate ky
  e0 = 0.0 # initalize
  for x in kxs:
    for y in kys:
      hk = hk_gen([x,y]) # get hamiltonian
      es = eigvalsh(hk) # eigenvalues
      e0 += np.sum(es[es<0.]) # add contribution
  return e0/nk**2
