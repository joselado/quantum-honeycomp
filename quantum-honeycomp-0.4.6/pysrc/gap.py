import numpy.linalg as lg
from scipy.optimize import minimize_scalar
import numpy as np
import scipy.sparse.linalg as lgs
from scipy.sparse import csc_matrix

def minimize_gap(f,tol=0.001):
  """Miimizes the gap of the system, the argument is between 0 and 1"""
  return f(minimize_scalar(f,method="Bounded",bounds=(0.,1.),tol=tol).x)





def gap_line(h,kpgen,assume_eh = False,sparse=True):
  """Return a function with argument between 0,1, which returns the gap"""
  hk_gen = h.get_hk_gen() # get hamiltonian generator
  def f(k):
    kp = kpgen(k) # get kpoint
    hk = hk_gen(kp) # generate hamiltonian
    if sparse: 
      es,ew = lgs.eigsh(csc_matrix(hk),k=4,which="LM",sigma=0.0)
    else:
      es = lg.eigvalsh(hk) # get eigenvalues
    if assume_eh: g = np.min(es[es>0.])
    else:  g = np.min(es[es>0.]) - np.max(es[es<0.])
    return g  # return gap
  return f  # return gap


def raw_gap(h,kpgen,sparse=True,nk=100):
  hk_gen = h.get_hk_gen() # get hamiltonian generator
  ks = np.linspace(0.,1.,nk)
  etot = [] # total eigenvalues
  for k in ks:
    kp = kpgen(k)
    hk = hk_gen(kp) # generate hamiltonian
    if sparse: 
      es,ew = lgs.eigsh(csc_matrix(hk),k=4,which="LM",sigma=0.0)
    else:
      es = lg.eigvalsh(hk) # get eigenvalues
    etot.append(es)
  etot = np.array(etot)
  return min(etot[etot>0.])


def gap2d(h,nk=40,k0=np.array([0.,0.]),rmap=1.0,recursive=False,
           iterations=5):
  """Calculates the gap for a 2d Hamiltonian by doing
  a kmesh sampling"""
  if h.dimensionality != 2: raise
  hk_gen = h.get_hk_gen() # get hamiltonian generator
  emin = 1000. # initial values
  for ix in np.linspace(-.5,.5,nk):  
    for iy in np.linspace(-.5,.5,nk):  
      k = np.array([ix,iy]) # generate kvector
      if recursive: k = k0 + k*rmap # scale vector
      hk = hk_gen(k) # generate hamiltonian
      es = lg.eigvalsh(hk) # get eigenvalues
      es = es[es>0.] # retain positive
      if min(es)<emin:
        emin = min(es) # store new minimum 
        kbest = k.copy() # store the best k
  if recursive: # if it has been chosen recursive
    if iterations>0: # if still iterations left
      emin = gap2d(h,nk=nk,k0=kbest,rmap=rmap/4,recursive=recursive,
                      iterations=iterations-1)
  return emin # gap



