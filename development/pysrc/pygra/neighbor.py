from __future__ import print_function
import numpy as np
from scipy.sparse import csc_matrix,bmat


minimum_hopping = 1e-3





try: 
  from . import first_neighborsf90
  def find_first_neighbor(r1,r2):
      """Calls the fortran routine"""
      from . import first_neighborsf90 as fn
      r1t = np.matrix(r1).T
      r2t = np.matrix(r2).T
      nn = fn.number_neighborsf90(r1t,r2t)
      if nn==0: return []  # if no neighbors found
      pairs = np.array(fn.first_neighborsf90(r1t,r2t,nn))
      return pairs.T # return the pairs


except:
  print("ERROR, neighbor fortran routine is not well compiled")
  def find_first_neighbor(r1,r2):
     """Calls the fortran routine"""
     print("Using ultraslow function!!!!!!!!!!")
     pairs = []
     for i in range(len(r1)):
       for j in range(len(r2)):
         ri = r1[i]
         rj = r2[j]
         dr = ri-rj
         dr = dr.dot(dr)
         if 0.8<dr<1.2: pairs.append([i,j])
     return np.array(pairs)



def connections(r1,r2):
  """Return a list with the connections of each atom"""
  pairs = find_first_neighbor(r1,r2) # get the pairs of first neighbors
  out = [[] for i in range(len(r1))] # output list
  for p in pairs:
    out[p[0]].append(p[1]) 
  return out # return list





def parametric_hopping(r1,r2,fc,is_sparse=False):
  """ Generates a parametric hopping based on a function"""
  if is_sparse: # sparse matrix
    print("Sparse parametric hopping")
    m = np.matrix([[0.0j for i in range(len(r2))] for j in range(len(r1))])
    rows,cols,data = [],[],[]
    for i in range(len(r1)):
      for j in range(len(r2)):
        val = fc(r1[i],r2[j]) # add hopping based on function
        if abs(val) > minimum_hopping: # retain this hopping
         data.append(val)
         rows.append(i)
         cols.append(j)
    n = len(r2)
    m = csc_matrix((data,(rows,cols)),shape=(n,n))
  #  if not is_sparse: m = m.todense() # dense matrix
    return m
  else:
    n = len(r2)
    m = np.matrix(np.zeros((n,n),dtype=np.complex)) # complex matrix
    for i in range(len(r1)):
      for j in range(len(r2)):
        m[i,j] = fc(r1[i],r2[j])
    return m

 





def parametric_hopping_spinful(r1,r2,fc,is_sparse=False):
  """ Generates a parametric hopping based on a function, that returns
  a 2x2 matrix"""
  m = [[None for i in range(len(r2))] for j in range(len(r1))]
  for i in range(len(r1)):
    for j in range(len(r2)):
      val = fc(r1[i],r2[j]) # add hopping based on function
      m[i][j] = val # store this result
  m = bmat(m) # convert to matrix
  if not is_sparse: m = m.todense() # dense matrix
  return m







def generate_parametric_hopping(h,f=None,mgenerator=None,
             spinful_generator=False):
  """ Adds a parametric hopping to the hamiltonian based on an input function"""
  rs = h.geometry.r # positions
  g = h.geometry # geometry
  has_spin = h.has_spin # check if it has spin
  is_sparse = h.is_sparse
  if mgenerator is None: # no matrix generator given on input
    if f is None: raise # no function given on input
    if spinful_generator:
      raise
      print("WARNING, I am not sure why I programmed this")
      h.has_spin = True
      generator = parametric_hopping_spinful
    else:
      h.has_spin = False
      generator = parametric_hopping
    def mgenerator(r1,r2):
      return generator(r1,r2,f,is_sparse=is_sparse)
  else:
    if h.dimensionality==3: raise
  h.intra = mgenerator(rs,rs)
  if h.dimensionality == 0: pass
  elif h.dimensionality == 1:
    dr = g.a1
    h.inter = mgenerator(rs,rs+dr)
  elif h.dimensionality == 2:
    h.tx = mgenerator(rs,rs+g.a1)
    h.ty = mgenerator(rs,rs+g.a2)
    h.txy = mgenerator(rs,rs+g.a1+g.a2)
    h.txmy = mgenerator(rs,rs+g.a1-g.a2)
  elif h.dimensionality == 3:
    if spinful_generator: raise # not implemented
    h.is_multicell = True # multicell Hamiltonian
    from . import multicell
    multicell.parametric_hopping_hamiltonian(h,fc=f)
  else: raise
  # check that the sparse mde is set ok
  if is_sparse and type(h.intra)==type(np.matrix([[]])):
    print("Matrices should be sparse, fixing")
    h.is_sparse = False
    h.turn_sparse() # turn the matrix sparse
  if not is_sparse and type(h.intra)!=type(np.matrix([[]])):
    print("Matrices should be dense, fixing",type(h.intra),type(np.matrix))
    h.is_sparse = True
    h.turn_dense() # turn the matrix sparse
  if has_spin: # Hamiltonian should be spinful
    print("Adding spin degree of freedom")
    h.has_spin = False
    h.turn_spinful()
  return h

