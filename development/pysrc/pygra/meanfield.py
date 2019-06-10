from __future__ import print_function
from scipy.sparse import csc_matrix,bmat
import numpy as np
#from .scftypes import selfconsistency

dup = csc_matrix(([1.0],[[0],[0]]),shape=(2,2),dtype=np.complex)
ddn = csc_matrix(([1.0],[[1],[1]]),shape=(2,2),dtype=np.complex)
sp = csc_matrix(([1.0],[[1],[0]]),shape=(2,2),dtype=np.complex)
sm = csc_matrix(([1.0],[[0],[1]]),shape=(2,2),dtype=np.complex)
def zero(d): return csc_matrix(([],[[],[]]),shape=(d,d),dtype=np.complex)


def element(i,n,p,d=2,j=None):
  if j is None: j=i
  o = csc_matrix(([1.0],[[d*i+p[0]],[d*j+p[1]]]),
                   shape=(n*d,n*d),dtype=np.complex)  
  return o
  o = csc_matrix(([1.0],[[p[0]],[p[1]]]),shape=(d,d),dtype=np.complex)
  m = [[None for i1 in range(n)] for i2 in range(n)]
  for i1 in range(n): m[i1][i1] = zero(d)
  m[i][j] = o.copy()
  return bmat(m) # return matrix



class interaction(): 
  def __init__(self):
    self.contribution = "AB" # assume that both operators contribute
    self.dir = [0,0,0] # unit cell to which to hop
    self.dhop = [0,0,0] # matrix that is affected

def hubbard_density(i,n,g=1.0):
  """Return pair of operators for a Hubbard mean field"""
  v = interaction()
  v.a = element(i,n,[0,0]) # value of the coupling
  v.b = element(i,n,[1,1]) # value of the coupling
  v.dir = [0,0,0] # direction of the interaction 
  v.g = g
  v.i = i
  v.j = i
  return v



def spinless_hubbard_density(i,n,g=1.0):
  """Return pair of operators for a Hubbard mean field"""
  v = interaction()
  v.a = element(i,n,[0,0],d=1) # value of the coupling
  v.b = element(i,n,[0,0],d=1) # value of the coupling
  v.dir = [0,0,0] # direction of the interaction 
  v.g = g
  v.i = i
  v.j = i
  return v







def hubbard_exchange(i,n,g=1.0):
  """Return pair of operators for a Hubbard mean field"""
  v = interaction()
  v.a = element(i,n,[0,1]) # value of the coupling
  v.b = element(i,n,[1,0]) # value of the coupling
  v.dir = [0,0,0] # direction of the interaction 
  v.g = -g # minus from fermion
  v.i = i
  v.j = i
  return v




def hubbard_pairing_ud(i,n,g=1.0):
  """Return pair of operators for a Hubbard mean field"""
  v = interaction()
  v.a = element(i,n,[0,2],d=4) # cc
  v.b = element(i,n,[2,0],d=4) # cdcd
  v.dir = [0,0,0] # direction of the interaction 
  v.g = g
  v.i = i
  v.j = i
  return v


def hubbard_pairing_du(i,n,g=1.0):
  """Return pair of operators for a Hubbard mean field"""
  v = interaction()
  v.a = element(i,n,[1,3],d=4) # cc
  v.b = element(i,n,[3,1],d=4) # cdcd
  v.dir = [0,0,0] # direction of the interaction 
  v.g = g
  v.i = i
  v.j = i
  return v




def v_pairing_uu(i,j,n,g=1.0,d=[0,0,0],channel="ee"):
  """Return pair of operators for a Hubbard mean field"""
  v = interaction()
  if channel=="ee": # ee channel
    v.a = element(i,n,[0,3],d=4,j=j) # cc
    v.b = element(j,n,[3,0],d=4,j=i) # cdcd
  elif channel=="hh":
    v.a = element(i,n,[3,0],d=4,j=j) # cc
    v.b = element(j,n,[0,3],d=4,j=i) # cdcd
  else: raise
  v.dir = d # direction of the interaction 
  v.g = g
  v.contribution = "A"
  v.i = i
  v.j = j
  return v


def v_pairing_dd(i,j,n,g=1.0,d=[0,0,0],channel="ee"):
  """Return pair of operators for a V mean field"""
  v = interaction()
  if channel=="ee": # ee channel
    v.a = element(i,n,[1,2],d=4,j=j) # cc
    v.b = element(j,n,[2,1],d=4,j=i) # cdcd
  elif channel=="hh":
    v.a = element(i,n,[2,1],d=4,j=j) # cc
    v.b = element(j,n,[1,2],d=4,j=i) # cdcd
  else: raise
  v.dir = d # direction of the interaction 
  v.g = g
  v.contribution = "A"
  v.i = i
  v.j = j
  return v


def v_pairing_du(i,j,n,g=1.0,d=[0,0,0]):
  """Return pair of operators for a V mean field"""
  v = interaction()
  v.a = element(i,n,[1,3],d=4,j=j) # cc
  v.b = element(j,n,[3,1],d=4,j=i) # cdcd
  v.dir = d # direction of the interaction 
  v.g = g
  v.i = i
  v.j = j
  return v




def v_ij(i,j,n,g=1.0,d=[0,0,0],spini=0,spinj=0):
  """Return pair of operators for a V mean field"""
  v = interaction()
  v.a = element(i,n,[spini,spini],d=2,j=j) # cc
  v.b = element(j,n,[spinj,spinj],d=2,j=i) # cdc
  v.dir = d # direction of the interaction 
  v.g = -g # this minus comes from commutation relations
  v.contribution = "A"
  v.i = i
  v.j = j
  return v


def v_ij_spinless(i,j,n,g=1.0,d=[0,0,0]):
  """Return pair of operators for a V mean field"""
  v = interaction()
  v.a = csc_matrix(([1.0],[[i],[j]]),shape=(n,n),dtype=np.complex) # cc
  v.b = csc_matrix(([1.0],[[j],[i]]),shape=(n,n),dtype=np.complex) # cdc
  v.dir = d # direction of the interaction 
  v.dhop = d # direction of the interaction 
  v.g = -g # this minus comes from commutation relations
  v.contribution = "A"
  v.i = i
  v.j = j
  return v



def v_ij_density_spinless(i,j,n,g=1.0,d=[0,0,0]):
  """Return pair of operators for a V mean field"""
  v = interaction()
  v.a = csc_matrix(([1.0],[[i],[i]]),shape=(n,n),dtype=np.complex) # cc
  v.b = csc_matrix(([1.0],[[j],[j]]),shape=(n,n),dtype=np.complex) # cdc
  v.dir = d # direction of the neighbor
  v.dhop = [0,0,0] # hopping that is affected
  v.g = g 
  v.contribution = "AB"
  v.i = i
  v.j = j
  return v




def guess(h,mode="ferro",fun=0.01):
  """Return a mean field matrix guess given a certain Hamiltonian"""
  h0 = h.copy() # copy Hamiltonian
  h0 = h0.get_multicell() # multicell
#  h0.intra *= 0. # initialize
  h0.clean() # clean the Hamiltonian
  if mode=="ferro":
    h0.add_zeeman(fun)
  elif mode=="random":
    h0.add_magnetism([np.random.random(3)-0.5 for i in h.geometry.r])
  elif mode=="potential":
    h0.add_onsite(fun)
  elif mode=="antiferro":
    h0.add_antiferromagnetism(fun)
  elif mode=="imbalance":
    h0.add_sublattice_imbalance(fun)
  elif mode=="swave":
    h0.add_swave(fun)
  elif mode=="pwave":
    for t in h0.hopping: t.m *= 0. # clean
    h0.add_pwave(fun)
    hop = dict()
    hop[(0,0,0)] = h0.intra
    for t in h0.hopping: hop[tuple(t.dir)] = t.m
    return hop
  else: raise
  return h0.intra # return matrix

from .algebra import braket_wAw
#from numba import jit

#@jit
def expectation_value(wfs,A,phis):
  """Return the expectation value of a set of wavevectors"""
  out = 0.0j
  for (p,w) in zip(phis,wfs):
    out += braket_wAw(w,A)*p
#    w = np.matrix(w) # convert to matrix
#    out += ((w.T).H*A*w.T)[0,0]*p # expectation value
  return np.conjugate(out) # return value




def enforce_pwave(mf):
  """Enforce pwave symmetry in a mean field Hamiltonian"""
  for key in mf: mf[key] = np.matrix(mf[key]) # dense matrix
  for key in mf: 
#    print(mf[key],type(mf[key]))
    n = mf[key].shape[0]//4 # number of sites 
    dm = tuple([-di for di in key]) # the opposite matrix
    for i in range(n): # loop over positions
      for j in range(n): # loop over positions
        for (ii,jj) in [(1,2),(0,3),(2,1),(3,0)]:
          mf[key][4*i+ii,4*j+jj] = (mf[key][4*i+ii,4*j+jj] - mf[dm][4*j+ii,4*i+jj])/2.
  return mf


def enforce_eh(h,mf):
  """Enforce eh symmetry in a mean field Hamiltonian"""
  from .superconductivity import eh_operator
  eh = eh_operator(h.intra) # get the function
  for key in mf: mf[key] = np.matrix(mf[key]) # dense matrix
  mfout = dict()
  for key in mf:
    mkey = (-key[0],-key[1],-key[2]) 
    mfout[key] = (mf[key] - eh(mf[mkey].H))/2.
  return mfout



def broyden_solver(scf):
    """Broyden solver for selfconsistency"""
    scf.mixing = 1.0 # perfect mixing
    scf.iterate() # do one iteration
    x0 = scf.cij # get the array with expectation values
    def fun(x): # function to solve
        scf.cij = x # impose those expectation values
        scf.cij2v() # update mean field
        scf.iterate() # perform an iteration
        return scf.cij - x # difference with respect to input
    from scipy.optimize import broyden1,broyden2,anderson
#    x = anderson(fun,x0)
    x = broyden1(fun,x0,f_tol=1e-7) # broyden solver
    scf.cij = x # impose those expectation values
    scf.cij2v() # update mean field
    scf.iterate() # perform an iteration
    return scf # return scf




def coulomb_interaction(g,**kwargs):
    """Return a list with the Coulomb interaction terms"""
    from .selfconsistency.coulomb import coulomb_density_matrix
    m = coulomb_density_matrix(g,**kwargs)
    interactions = [] # empty list
    nat = len(g.r)
    for i in range(nat):
      for j in range(nat):
        if abs(m[i,j])>1e-5:
          interactions.append(v_ij_density_spinless(i,j,nat,
            g=m[i,j],d=[0,0,0]))
    return interactions





