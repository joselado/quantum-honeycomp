import numpy as np

def find_first_neighbor(r1,r2):
   """Calls the fortran routine"""
   import first_neighborsf90 as fn
   r1t = r1.T
   r2t = r2.T
   nn = fn.number_neighborsf90(r1t,r2t)
   print nn
   if nn==0: return []  # if no neighbors found
   pairs = np.array(fn.first_neighborsf90(r1t,r2t,nn))
   return pairs.T # return the pairs


