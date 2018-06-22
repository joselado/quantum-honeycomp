from __future__ import print_function
import numpy as np
from scipy.sparse import coo_matrix, bmat, csc_matrix
import scipy.sparse as sp



def get_eh_sector_odd_even(m,i=0,j=0):
  """ Return the electron hole sector of a matrix, 
  assuming the form is even index for electrons and odd for holes"""
  if i>1 or j>1: raise # meaningless
  n = m.shape[0]//2 # dimension of the output
  nat = m.shape[0]//4 # number of sites
  mout = np.matrix(np.zeros((n,n)),dtype=np.complex) # define matrix
  for ii in range(nat): # loop over index
    for jj in range(nat): # loop over index
      if i==0: indi = 4*ii # electron sector
      else: indi = 4*ii + 2 # electron sector
      if j==0: indj = 4*jj # electron sector
      else: indj = 4*jj + 2 # electron sector
      mout[2*ii,2*jj] = m[indi,indj] # assign
      if i==0: indi = 4*ii+1 # electron sector
      else: indi = 4*ii + 3 # electron sector
      if j==0: indj = 4*jj+1 # electron sector
      else: indj = 4*jj + 3 # electron sector
      mout[2*ii+1,2*jj+1] = m[indi,indj] # assign
  return mout # return matrix



def get_nambu_tauz(m,has_eh=False):
  """Return the nambu matrix tauz for electron-hole"""
  raise # this function is not consistent with the Nambu notation (see below)
  n = m.shape[0] # number of sites 
  if has_eh: n = n//2 # half
  mout = np.matrix(np.zeros((n*2,n*2)),dtype=np.complex) # define matrix
  for ii in range(n): # loop over index
    mout[2*ii,2*ii] = 1. # assign
    mout[2*ii+1,2*ii+1] = -1. # assign
  return mout # return tauz



def project_electrons(m):
  """Return the nambu matrix tauz for electron-hole"""
  n = m.shape[0] # number of sites 
  mout = m*0.0 # define matrix
  for ii in range(n): # loop over index
    for jj in range(n): # loop over index
      if ii%4<2 and jj%4<2:  mout[ii,jj] = m[ii,jj] # assign
      else: continue
  return mout # return tauz



def project_holes(m):
  """Return the nambu matrix tauz for electron-hole"""
  n = m.shape[0] # number of sites 
  mout = m*0.0 # define matrix
  for ii in range(n): # loop over index
    for jj in range(n): # loop over index
      if ii%4>1 and jj%4>1:  mout[ii,jj] = m[ii,jj] # assign
      else: continue
  return mout # return tauz




def eh_operator(m):
  """Return the electron hole symmetry operator, as a function"""
  n = m.shape[0]//4 # number of sites 
  out = [[None for i in range(n)] for j in range(n)] # output matrix
  tau = csc_matrix([[0.,0.,0.,-1.],[0,0,1,0],[0,1,0,0],[-1,0,0,0]])
  for i in range(n):
    out[i][i] = tau
  out = bmat(out) # sparse matrix
  def ehop(inm):
    """Function that applies electron-hole symmetry to a matrix"""
    return out*np.conjugate(inm)*out # return matrix
  return ehop # return operator
    





def non_unitary_pairing(v):
  """Calculate the vector that defines the non-unitary pairing,
  the input is the three components of the pairing matrix"""
  vc = [np.conjugate(v[i]) for i in range(3)]
  ux = v[1]*vc[2] - v[2]*vc[1]
  uy = v[2]*vc[0] - v[0]*vc[2]
  uz = v[0]*vc[1] - v[1]*vc[0]
  return (1j*np.array([ux,uy,uz])).real



# from now on, the spinors will be
# Psi =
# Psi_up
# Psi_dn
# Psi^dag_dn
# -Psi^dag_up


spair = csc_matrix([[0.,0.,1.,0.],[0.,0.,0.,1.],[0.,0.,0.,0.],[0.,0.,0.,0.]])
proje = csc_matrix([[1.,0.,0.,0.],[0.,1.,0.,0.],[0.,0.,0.,0.],[0.,0.,0.,0.]])
projh = csc_matrix([[0.,0.,0.,0.],[0.,0.,0.,0.],[0.,0.,1.,0.],[0.,0.,0.,1.]])
#pzpair = csc_matrix([[0.,0.,0.,1.],[0.,0.,-1.,0.],[0.,0.,0.,0.],[0.,0.,0.,0.]])
#pxpair = csc_matrix([[0.,0.,1.,0.],[0.,0.,0.,-1.],[0.,0.,0.,0.],[-0.,0.,0.,0.]])
#pypair = csc_matrix([[0.,0.,-1j,0.],[0.,0.,0.,1j],[0.,0,0.,0.],[0,0.,0.,0.]])
# Definition using Manfred's notes
# This are the operators that give the d-vectors
# Minus sign due to the nambu spinor!
deltauu = csc_matrix([[0.,0.,0.,0.],[0.,0.,0.,0.],[0.,0.,0.,0.],[-1.,0.,0.,0.]])
# and the rest are simple
deltadd = csc_matrix([[0.,0.,0.,0.],[0.,0.,0.,0.],[0.,1.,0.,0.],[0.,0.,0.,0.]])
deltauu = deltauu.H
deltadd = deltadd.H
deltax = (deltadd - deltauu)/2.
deltay = (deltadd + deltauu)/2j
# this one is tricky, we only want the antisymmetric part
deltaz = csc_matrix([[0.,0.,0.,0.],[0.,0.,0.,0.],[1.,0.,0.,0.],[0.,-1.,0.,0.]])/2.

deltaz = deltaz.H



# redefine the functions according to the previous notation

def time_reversal(m):
  """Do the spinful time reversal of this matrix"""
  from hamiltonians import sy
  n = m.shape[0]//2 # number of spinful blocks
  msy = [[None for ri in range(n)] for j in range(n)]
  for i in range(n): msy[i][i] = sy # add sy
  msy = bmat(msy) # as block matrix
  return msy*np.conjugate(m)*msy # return time reversal




def build_eh(hin,coupling=None,is_sparse=True):
  if coupling is not None:
    c12 = coupling
    c21 = coupling.H
  else:
    c12 = None
    c21 = None
  return build_nambu_matrix(hin,is_sparse=is_sparse,c12=c12,c21=c21)

#  n = hin.shape[0]  # dimension of input
#  bout = [[None,None],[None,None]] # initialize None matrix
#  bout[0][0] = csc_matrix(hin) # electron part
#  bout[1][1] = -time_reversal(csc_matrix(hin)) # hole part
#  if coupling is not None:
#    bout[0][1] = csc_matrix(coupling) # pairing part
#    bout[1][0] = csc_matrix(coupling).H # pairing part
#  bout = bmat(bout) # return matrix
#  out = reorder(bout) # reorder matrix
#  if is_sparse: return out
#  else: return out.todense()



def build_nambu_matrix(hin,c12=None,c21=None,is_sparse=True):
  n = hin.shape[0]  # dimension of input
  bout = [[None,None],[None,None]] # initialize None matrix
  bout[0][0] = csc_matrix(hin) # electron part
  bout[1][1] = -time_reversal(csc_matrix(hin)) # hole part
  if c12 is not None:
    bout[0][1] = csc_matrix(c12) # pairing part
  if c21 is not None:
    bout[1][0] = csc_matrix(c21) # pairing part
  bout = bmat(bout) # return matrix
  out = reorder(bout) # reorder matrix
  if is_sparse: return out
  else: return out.todense()






def add_swave(delta=0.0,is_sparse=False,rs=None):
  """ Adds swave pairing """
  if rs is None: raise # raise error to signal that this is temporal
  n = len(rs) # number of sites
  if callable(delta): # delta is a function
    datar = [delta(ri) for ri in rs] # generate data for the different positions
    data = []
    # the coupling involves up and down
    for dr in datar: # loop over positions
      data.append(dr)  # up e dn h
      data.append(dr)  # dne up h
    iis = range(n*2) # indexes
    coupling = csc_matrix((data,(iis,iis)),dtype=np.complex)  # generate matrix
  else:
    coupling = sp.identity(n*2)*delta # delta matrix
  zero = coupling*0.
  return build_eh(zero,coupling=coupling,is_sparse=is_sparse) # return matrix


def add_pxipy(delta=0.0,is_sparse=False,r1=None,r2=None):
  """Add px x + py y pairing"""
  def deltafun(r1i,r2i):
    """Function to calculate the pairing"""
    dr = r2i-r1i
    dr2 = dr.dot(dr)
    if 0.9<dr2<1.1: # first neighbor
      dr = delta*dr/np.sqrt(dr2) # unit vector
#      dr = np.cross(dr,np.array([0.,0.,1.]))
      dr = [dr[0],1j*dr[1],0.]      
#      dr = [0.,0.,dr[0]+1j*dr[1]]
#      return dr
      return dvector2deltas(dr) # return delta
    else: return [0.,0.,0.] # zero vector
  return add_pairing(deltas=deltafun,r1=r1,r2=r2)


add_pwave = add_pxipy


def dvector2deltas(ds):
  """Transform a certain dvector into deltauu, deltadd and deltaud"""
  deltas = [0.,0.,0.]
  deltas[0] = ds[0]+ds[1]
  deltas[1] = -1j*(ds[0]-ds[1]) # this sign might not be ok
  deltas[2] = ds[2]
  return deltas


def add_pairing(deltas=[0.,0.,0.],is_sparse=True,r1=[],r2=[]):
  """ Adds a general pairing in real space"""
  def get_pmatrix(r1i,r2j): # return the different pairings
    if callable(deltas): dv = deltas(r1i,r2j) # get the components
    else: dv = deltas
    duu = dv[0] # delta up up
    ddd = dv[1] # delta dn dn
    dud = dv[2] # delta up dn
    # be aware of the minus signs coming from the definition of the
    # Nambu spinor!!!!!!!!!!!!!!!!!
    # c_up d_dn d^\dagger_dn -d^\dagger_up
    D = csc_matrix([[dud,-duu],[ddd,-dud]]) # SC matrix
    D = bmat([[None,D],[-D.H,None]]) # the minus sign comes from triplet
    return D
  # superconducting coupling
  n = len(r1)  # number of sites
  bout = [[None for i in range(n)] for j in range(n)] # initialize None matrix
  # zeros in the diagonal
  for i in range(n): bout[i][i] = csc_matrix(np.zeros((4,4),dtype=np.complex))
  for i in range(n): # loop over sites
    for j in range(n): # loop over sites
      bout[i][j] = get_pmatrix(r1[i],r2[j]) # get this pairing
  return bmat(bout) # return pairing matrix







def reorder(m):
  '''Reorder a matrix that has electrons and holes, so that they are
  ordered in an analogous way as the initial matrix'''
  R = np.matrix(np.zeros(m.shape)) # zero matrix
  nr = m.shape[0]//4 # number of positions
  for i in range(nr): # electrons
    R[2*i,4*i] = 1.0 # up electron
    R[2*i+1,4*i+1] = 1.0 # down electron
    R[2*i+2*nr,4*i+2] = 1.0 # down holes
    R[2*i+1+2*nr,4*i+3] = 1.0 # up holes
  R = csc_matrix(R) # to sparse
  return R.H*m*R


