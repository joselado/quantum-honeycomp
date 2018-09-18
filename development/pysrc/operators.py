# library to create operators
from __future__ import division
import numpy as np
from scipy.sparse import csc_matrix as csc
from scipy.sparse import bmat
from superconductivity import build_eh
import superconductivity
import scipy.linalg as lg
#from bandstructure import braket_wAw
import current


def braket_wAw(w,A):
  w = np.matrix(w) # convert to matrix
  return ((w.T).H*A*w.T)[0,0] # expectation value




def index(h,n=[0]):
  """Return a projector onto a site"""
  num = len(h.geometry.r)
  val = [1. for i in n]
  m = csc((val,(n,n)),shape=(num,num),dtype=np.complex)
  return h.spinless2full(m) # return matrix






def operator2list(operator):
  """Convert an input operator in a list of operators"""
  if operator is None: # no operator given on input
    operator = [] # empty list
  elif not isinstance(operator,list): # if it is not a list
    operator = [operator] # convert to list
  return operator



def surface(h,cut = 2.,which="both"):
  """Return an operator which is non-zero in the upper surface"""
  zmax = np.max(h.geometry.z) # maximum z
  zmin = np.min(h.geometry.z) # maximum z
  dind = 1 # index to which divide the positions
  if h.has_spin:  dind *= 2 # duplicate for spin
  if h.has_eh:  dind *= 2  # duplicate for eh
  n = h.intra.shape[0] # number of elments of the hamiltonian
  data = [] # epmty list
  for i in range(n): # loop over elements
    z = h.geometry.z[i//dind]
    if which=="upper": # only the upper surface
      if np.abs(z-zmax) < cut:  data.append(1.)
      else: data.append(0.)
    elif which=="lower": # only the upper surface
      if np.abs(z-zmin) < cut:  data.append(1.)
      else: data.append(0.)
    elif which=="both": # only the upper surface
      if np.abs(z-zmax) < cut:  data.append(1.)
      elif np.abs(z-zmin) < cut:  data.append(-1.)
      else: data.append(0.)
    else: raise
  row, col = range(n),range(n)
  m = csc((data,(row,col)),shape=(n,n),dtype=np.complex)
  return m # return the operator



def interface1d(h,cut = 3.):
  dind = 1 # index to which divide the positions
  if h.has_spin:  dind *= 2 # duplicate for spin
  if h.has_eh:  dind *= 2  # duplicate for eh
  n = h.intra.shape[0] # number of elments of the hamiltonian
  data = [] # epmty list
  for i in range(n): # loop over elements
    y = h.geometry.y[i//dind]
    if np.abs(y)<cut: data.append(1.) # if it belongs to the interface
    else:  data.append(0.)  # otherwise
  row, col = range(n),range(n)
  m = csc((data,(row,col)),shape=(n,n),dtype=np.complex)
  return m # return the operator



def get_interface(h,fun=None):
  """Return an operator that projects onte the interface"""
  dind = 1 # index to which divide the positions
  if h.has_spin:  dind *= 2 # duplicate for spin
  if h.has_eh:  dind *= 2  # duplicate for eh
  iden = csc(np.matrix(np.identity(dind,dtype=np.complex))) # identity matrix
  r = h.geometry.r # positions
  out = [[None for ri in r] for rj in r] # initialize
  if fun is None: # no input function
    cut = 2.0 # cutoff
    if h.dimensionality==1: index = 1
    elif h.dimensionality==2: index = 2
    else: raise
    def fun(ri): # define the function
      if np.abs(ri[index])<cut: return 1.0
      else: return 0.0
  for i in range(len(r)): # loop over positions
    out[i][i] = fun(r[i])*iden 
  return bmat(out) # return matrix



def get_pairing(h,ptype="s"):
  """Return an operator that calculates the expectation value of the
  s-wave pairing"""
  if not h.has_eh: raise # only for e-h systems
  if ptype=="s": op = superconductivity.spair
  elif ptype=="deltax": op = superconductivity.deltax
  elif ptype=="deltay": op = superconductivity.deltay
  elif ptype=="deltaz": op = superconductivity.deltaz
  else: raise
  r = h.geometry.r
  out = [[None for ri in r] for rj in r]
  for i in range(len(r)): # loop over positions
    out[i][i] = op
  return bmat(out) # return matrix



def get_electron(h):
  """Operator to project on the electron sector"""
  if not h.has_eh: raise # only for e-h systems
  op = superconductivity.proje
  r = h.geometry.r
  out = [[None for ri in r] for rj in r]
  for i in range(len(r)): # loop over positions
    out[i][i] = op
  return bmat(out)


def get_hole(h):
  """Operator to project on the electron sector"""
  if not h.has_eh: raise # only for e-h systems
  op = superconductivity.projh
  r = h.geometry.r
  out = [[None for ri in r] for rj in r]
  for i in range(len(r)): # loop over positions
    out[i][i] = op
  return bmat(out)





def bulk1d(h,p = 0.5):
  dind = 1 # index to which divide the positions
  if h.has_spin:  dind *= 2 # duplicate for spin
  if h.has_eh:  dind *= 2  # duplicate for eh
  n = h.intra.shape[0] # number of elments of the hamiltonian
  data = [] # epmty list
  cut = np.max(h.geometry.y)*p
  for i in range(n): # loop over elements
    y = h.geometry.y[i//dind]
    if y < -cut:  data.append(-1.)
    elif y > cut: data.append(1.)
    else: data.append(0.)
  row, col = range(n),range(n)
  m = csc((data,(row,col)),shape=(n,n),dtype=np.complex)
  return m # return the operator


def get_xposition(h):  return get_position(h,mode="x")
def get_yposition(h):  return get_position(h,mode="y")
def get_zposition(h):  return get_position(h,mode="z")




def get_position(h,mode="z"):
  dind = 1
  if h.has_spin:  dind *= 2 # duplicate for spin
  if h.has_eh:  dind *= 2  # duplicate for eh
  n = h.intra.shape[0] # number of elments of the hamiltonian
  if len(h.geometry.z)!=n//dind: raise # dimensions do not match
  data = [] # epmty list
  if mode=="x": pos = h.geometry.x
  elif mode=="y": pos = h.geometry.y
  elif mode=="z":  pos = h.geometry.z
  else: raise
  for i in range(n): # loop over elements
    z = pos[i//dind]
    data.append(z)
  row, col = range(n),range(n)
  m = csc((data,(row,col)),shape=(n,n),dtype=np.complex)
  return m # return the operator



from rotate_spin import sx,sy,sz # import pauli matrices
 

def get_si(h,i=1):
  """Return a certain Pauli matrix for the full Hamiltonian"""
  if i==1: si = sx # sx matrix
  elif i==2: si = sy # sy matrix
  elif i==3: si = sz # sz matrix
  else: raise # unknown pauli matrix
  if h.has_eh: ndim = h.intra.shape[0]//4 # half the dimension
  else: ndim = h.intra.shape[0]//2 # dimension
  if not h.has_spin: raise # it does not have spin
  if h.has_spin: # spinful system
    op = [[None for i in range(ndim)] for j in range(ndim)] # initialize
    for i in range(ndim): op[i][i] = si # store matrix
    op = bmat(op) # create matrix
  if h.has_eh: op = build_eh(op,is_sparse=True) # add electron and hole parts 
  return op

# define the functions for the three spin components
get_sx = lambda h: get_si(h,i=1) # sx matrix
get_sy = lambda h: get_si(h,i=2) # sy matrix
get_sz = lambda h: get_si(h,i=3) # sz matrix






def get_z(h):
  """Operator for the calculation of z expectation value"""
  if h.intra.shape[0]==len(h.geometry.z): # if as many positions as entries
    op = np.zeros(h.intra.shape,dtype=np.complex) # initialize matrix
    for i in range(len(h.geometry.z)):
      op[i,i] = h.geometry.z[i]
    return op
  raise
  if h.has_eh: raise
  if not h.has_spin: raise
  if h.has_spin:
    op = np.zeros(h.intra.shape,dtype=np.complex) # initialize matrix
    for i in range(len(op)//2):   
      op[2*i,2*i+1] = -1j
      op[2*i+1,2*i] = 1j
  



def get_rop(h,fun):
  """Operator for the calculation of a position expectation value"""
  rep = 1 # repetitions 
  if h.has_spin: rep *= 2
  if h.has_eh: rep *= 2
  data = []
  for ri in h.geometry.r: 
    for i in range(rep): data.append(fun(ri)) # store
  n = h.intra.shape[0]
  row = range(n)
  col = range(n)
  m = csc((data,(row,col)),shape=(n,n),dtype=np.complex)
  return m




def get_sublattice(h,mode="both"):
  """Sublattice operator"""
  if not h.geometry.has_sublattice: raise
  rep = 1 # repetitions 
  if h.has_spin: rep *= 2
  if h.has_eh: rep *= 2
  data = []
  for s in h.geometry.sublattice: 
    for i in range(rep): 
      if mode=="both": data.append(s) # store
      elif mode=="A": data.append((s+1.)/2.) # store
      elif mode=="B": data.append((-s+1.)/2.) # store
      else: raise
  n = h.intra.shape[0]
  row = range(n)
  col = range(n)
  m = csc((data,(row,col)),shape=(n,n),dtype=np.complex)
  return m


def get_velocity(h):
  """Return the velocity operator"""
  vk = current.current_operator(h)
  def f(w,k=[0.,0.,0.]):
    return braket_wAw(w,vk(k)).real
  return f


get_current = get_velocity

def get_spin_current(h):
  vk = current.current_operator(h)
  sz = get_sz(h)
  def f(w,k=[0.,0.,0.]):
    return braket_wAw(w,vk(k)).real*braket_wAw(w,sz).real
  return f





def get_valley(h,projector=False,delta=None):
  """Return a callable that calculates the valley expectation value
  using the modified Haldane coupling"""
  ho = h.copy() # copy Hamiltonian
  ho.clean() # set to zero
  ho.add_modified_haldane(1.0/4.5) # add modified Haldane coupling
  hkgen = ho.get_hk_gen() # get generator for the hk Hamiltonian
  from scipy.sparse import issparse
  def sharpen(m):
    """Sharpen the eigenvalues of a matrix"""
#    return m
    if delta is None: return m # do nothing
    if issparse(m): return m # temporal workaround
    if issparse(m): m = m.todense() # I should fix this
    (es,vs) = lg.eigh(m) # diagonalize
    es = es/(np.abs(es)+delta) # renormalize the valley eigenvalues
    vs = np.matrix(vs) # convert
    m0 = np.matrix(np.diag(es)) # build new hamiltonian
    return vs*m0*vs.H # return renormalized operator
  if projector: # function returns a matrix
    def fun(m,k=None):
      if h.dimensionality>0 and k is None: raise # requires a kpoint
      hk = hkgen(k) # evaluate Hamiltonian
      hk = sharpen(hk) # sharpen the valley
      return m*hk # return the projector
  else: # conventional way
    def fun(w,k=None):
      if h.dimensionality>0 and k is None: raise # requires a kpoint
      hk = hkgen(k) # evaluate Hamiltonian
      hk = sharpen(hk) # sharpen the valley
      return braket_wAw(w,hk).real # return the braket
  return fun # return function




def get_inplane_valley(h):
  """Returns an operator that computes the absolute value
  of the intervalley mixing"""
  ho = h.copy() # copy Hamiltonian
  ho.clean() # set to zero
  ho.add_modified_haldane(1.0/4.5) # add modified Haldane coupling
  hkgen = ho.get_hk_gen() # get generator for the hk Hamiltonian
  hkgen0 = h.get_hk_gen() # get generator for the hk Hamiltonian
  def fun(w,k=None):
#    return abs(np.sum(w*w))
    if h.dimensionality>0 and k is None: raise # requires a kpoint
    hk = hkgen(k) # evaluate Hamiltonian
    hk0 = hkgen0(k) # evaluate Hamiltonian
    A = hk*hk0 - hk0*hk # commutator
    A = -A*A
    return abs(braket_wAw(w,A)) # return the braket
  return fun # return function










