from __future__ import print_function
import numpy as np


def anderson(h,w=0.0,p=1.0):
  """Return a Hamiltonian with Anderson disorder"""
  if h.has_eh: raise # not yet
  if h.has_spin: raise # not yet
  ho = h.copy() # copy the Hamiltonian
  n = h.intra.shape[0] # dimension
  cs = (np.random.random(n) - .5)*2*w # disorder
  # decide if an impurity should be added
  ps = []
  for ip in np.random.random(n): # loop over random numbers
    if p<ip: ps.append(0.)
    else: ps.append(1.)
  ps = np.array(ps)
  ho.intra = ho.intra + np.diag(cs*ps)
  return ho # return Hamiltonian




def phase(h,w=0.0):
  """Random phase disorder"""
  if h.has_eh: raise # not yet
  if h.has_spin: raise # not yet
  ho = h.copy() # copy the Hamiltonian
  n = h.intra.shape[0] # dimension
  cs = np.array(np.random.random((n,n)) - .5)*2*w # disorder
  cs = (cs - cs.T)/2. # Hermitian
  cs = np.exp(1j*np.pi*cs) # phases
  ho.intra = np.matrix(np.array(ho.intra)*np.array(cs))
  return ho


