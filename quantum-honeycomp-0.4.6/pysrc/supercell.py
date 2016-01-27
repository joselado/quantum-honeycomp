

def intra_super2d(h,n=1,central=None):
  """ Calculates the intra of a 2d system"""
  from scipy.sparse import bmat
  from scipy.sparse import csc_matrix as csc
  tx = csc(h.tx)
  ty = csc(h.ty)
  txy = csc(h.txy)
  txmy = csc(h.txmy)
  intra = csc(h.intra)
  for i in range(n):
    intrasuper[i][i] = intra # intracell
    (x1,y1) = inds[i]
    for j in range(n):
      (x2,y2) = inds[j]
      dx = x2-x1
      dy = y2-y1
      if dx==1 and  dy==0: intrasuper[i][j] = tx
      if dx==-1 and dy==0: intrasuper[i][j] = tx.H
      if dx==0 and  dy==1: intrasuper[i][j] = ty
      if dx==0 and  dy==-1: intrasuper[i][j] = ty.H
      if dx==1 and  dy==1: intrasuper[i][j] = txy
      if dx==-1 and dy==-1: intrasuper[i][j] = txy.H
      if dx==1 and  dy==-1: intrasuper[i][j] = txmy
      if dx==-1 and dy==1: intrasuper[i][j] = txmy.H
  # substitute the central cell if it is the case
  if central!=None:
    ic = (n-1)/2
    intrasuper[ic][ic] = central
  # now convert to matrix
  intrasuper = bmat(intrasuper).todense() # supercell
  return intrasuper


def supercell1d(h,nsuper,sparse=False):
  """ Get a supercell of the system,
      h is the hamiltonian
      nsuper is the number of replicas """
  from scipy.sparse import csc_matrix as csc
  from scipy.sparse import bmat
  from copy import deepcopy
  hout = deepcopy(h) # output hamiltonian
  intra = csc(h.intra)  # intracell matrix
  inter = csc(h.inter)  # intercell matrix
  # crease supercells block matrices
  intra_super = [[None for i in range(nsuper)] for j in range(nsuper)]
  inter_super = [[None for i in range(nsuper)] for j in range(nsuper)]
  zero = csc(0.0*intra)
  for i in range(nsuper):  # onsite part
    intra_super[i][i] = intra
    inter_super[i][i] = zero
  for i in range(nsuper-1):  # inter minicell part
    intra_super[i][i+1] = inter
    intra_super[i+1][i] = inter.H
  inter_super[nsuper-1][0] = inter # inter supercell part
  inter_super[0][nsuper-1] = inter.H*0. # inter supercell part
  # create dense matrices
  if sparse:
    intra_super = bmat(intra_super)
    inter_super = csc(bmat(inter_super))
  else:
    intra_super = bmat(intra_super).todense()
    inter_super = bmat(inter_super).todense()
  # add to the output hamiltonian
  del intra
  del inter
  hout.intra = intra_super
  hout.inter = inter_super
  try: hout.geometry = h.geometry.supercell(nsuper)
  except: print "No geometry given"
  return hout



