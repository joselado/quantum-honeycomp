import numpy as np

def hk_gen(h):
  """Generate a k dependent hamiltonian"""
  if h.is_multicell==False: raise
  if h.dimensionality == 0: return None
  if h.dimensionality == 1: # one dimensional
    def hk(k):
      """k dependent hamiltonian, k goes from 0 to 1"""
      mout = h.intra # intracell term
      for t in h.hopping: # loop over matrices
        phi = t.dir[0]*k # phase
        tk = t.m * np.exp(1j*np.pi*2.*phi) # k hopping
        mout = mout + tk 
      return mout
  if h.dimensionality == 2: # two dimensional
    def hk(k):
      """k dependent hamiltonian, k goes from 0 to 1"""
      k = np.array([k[0],k[1]]) # convert to array
      mout = h.intra # intracell term
      for t in h.hopping: # loop over matrices
        d = t.dir
        d = np.array([d[0],d[1]]) # vector director of hopping
        phi = d.dot(k) # phase
        tk = t.m * np.exp(1j*np.pi*2.*phi) # k hopping
        mout = mout + tk 
      return mout

    return hk  # return the function

def turn_spinful(h):
  """Turn a hamiltonian spinful"""
  from increase_hilbert import spinful
  if h.has_eh: raise
  if h.has_spin: return # return if already has spin
  h.has_spin = True # put spin
  h.intra = spinful(h.intra) # spinful intra
  for i in range(len(h.hopping)): 
    h.hopping[i].m = spinful(h.hopping[i].m) # spinful hopping

