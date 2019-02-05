from scipy.sparse import issparse
from scipy.sparse import csc_matrix as csc
import scipy.linalg as dlg
import numpy as np


def braket_wAw(w,A,wi=None):
  """
  Compute the braket of a wavefunction
  """
  if wi is None: wi = w
  if issparse(A): # sparse matrices
    return (np.conjugate(wi)@A@w) # modern way
  else: # matrices and arrays
    return (np.conjugate(wi)@A@w)[0,0] # modern way


def disentangle_manifold(wfs,A):
  """
  Disentangles the wavefunctions of a degenerate manifold
  by expressing them in terms of eigenvalues of an input operator
  """
  ma = get_representation(wfs,A) # get the matrix form of the operator
  wfsout = [] # empty list
  evals,evecs = dlg.eigh(ma) # diagonalize

  evecs = evecs.transpose() # transpose eigenvectors
  for v in evecs: # loop over eigenvectors
    wf = wfs[0]*0.0j
    for (i,iv) in zip(range(len(v)),v): # loop over components
      wf += iv*wfs[i] # add contribution
    wfsout.append(wf.copy()) # store wavefunction
  return wfsout



def get_representation(wfs,A):
  """
  Gets the matrix representation of a certain operator
  """
  n = len(wfs) # number of eigenfunctions
  ma = np.zeros((n,n),dtype=np.complex) # representation of A
  sa = csc(A) # sparse matrix
  for i in range(n):
    vi = csc(np.conjugate(wfs[i])) # first wavefunction
    for j in range(n):
      vj = csc(wfs[j]).transpose() # second wavefunction
      data = (vi*sa*vj).todense()[0,0]
      ma[i,j] = data
  return ma






