from __future__ import print_function
import numpy as np
from scipy.sparse import csc_matrix

def equal(m1,m2):
  """Check if two matrices are the same"""
  if np.max(np.abs(m1-m2))>0.000001:
    print(csc_matrix(m1-m2))
#    print("\n")
#    print(csc_matrix(m2))
    return False
  else: return True

def check_hamiltonian(h):
  """Do various checks in Hamiltonian, to ensure that nothing weird happens"""
  hk = h.get_hk_gen() # get generator
  m = hk(np.random.random(3)) # random k-point
  if not equal(m,m.H):
    print("CHECK FAILED, Hamiltonian is not Hermitian")
    raise # not hermitian
  if h.has_eh: # if it has electron hole degree of freedom
    v = np.random.random(3) # random kpoint
    m1 = hk(v) # Hamiltonian
    m2 = hk(-v) # Hamiltonian in time reversal point
    from superconductivity import eh_operator
    eh = eh_operator(m1) # get the function
    if not equal(m1,-eh(m2)): 
      print("CHECK FAILED, Hamiltonian does not have electron-hole symmetry")
      raise
    print("CHECKED that the Hamiltonian has electron-hole symmetry")

