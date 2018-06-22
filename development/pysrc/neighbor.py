from __future__ import print_function
import numpy as np
from scipy.sparse import csc_matrix,bmat


minimum_hopping = 1e-3


def find_first_neighbor(r1,r2):
   """Calls the fortran routine"""
   import first_neighborsf90 as fn
   r1t = np.matrix(r1).T
   r2t = np.matrix(r2).T
   nn = fn.number_neighborsf90(r1t,r2t)
#   print nn
   if nn==0: return []  # if no neighbors found
   pairs = np.array(fn.first_neighborsf90(r1t,r2t,nn))
   return pairs.T # return the pairs



def connections(r1,r2):
  """Return a list with the connections of each atom"""
  pairs = find_first_neighbor(r1,r2) # get the pairs of first neighbors
  out = [[] for i in range(len(r1))] # output list
  for p in pairs:
    out[p[0]].append(p[1]) 
  return out # return list






try: 
  import first_neighborsf90
except:
  print("ERROR, neighbor fortran routine is not well compiled")
  def find_first_neighbor(r1,r2):
     """Calls the fortran routine"""
     print("Using ultraslow function!!!!!!!!!!")
     pairs = []
     for i in range(len(r1)):
       for j in range(len(r1)):
         ri = r1[i]
         rj = r2[j]
         dr = ri-rj
         dr = dr.dot(dr)
         if 0.8<dr<1.2: pairs.append([i,j])
     return np.array(pairs)



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

