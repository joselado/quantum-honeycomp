import numpy as np



from scipy.sparse import csc_matrix,bmat
sx = csc_matrix([[0.,1.],[1.,0.]])
sy = csc_matrix([[0.,-1j],[1j,0.]])
sz = csc_matrix([[1.,0.],[0.,-1.]])



def rotation_matrix(m,vectors):
  """ Rotates a matrix, to align its components with the direction
  of the magnetism """
  if not len(m)==2*len(vectors): # stop if they don't have
                                  # compatible dimensions
     raise
  # pauli matrices
  n = len(m)/2 # number of sites
  R = [[None for i in range(n)] for j in range(n)] # rotation matrix
  from scipy.linalg import expm  # exponenciate matrix
  for (i,v) in zip(range(n),vectors): # loop over sites
    vv = np.sqrt(v.dot(v)) # norm of v
    if vv>0.000001: # if nonzero scale
      u = v/vv
    else: # if zero put to zero
      u = np.array([0.,0.,0.,])
#    rot = u[0]*sx + u[1]*sy + u[2]*sz 
    uxy = np.sqrt(u[0]**2 + u[1]**2) # component in xy plane
    phi = np.arctan2(u[1],u[0])
    theta = np.arctan2(uxy,u[2])
    r1 =  phi*sz/2.0 # rotate along z
    r2 =  theta*sy/2.0 # rotate along y
    # a factor 2 is taken out due to 1/2 of S
    rot = expm(1j*r2) * expm(1j*r1)   
    R[i][i] = rot  # save term
  R = bmat(R)  # convert to full sparse matrix
  return R.todense()



def align_magnetism(m,vectors):
  """ Align matrix with the magnetic moments"""
  R = rotation_matrix(m,vectors) # return the rotation matrix
  mout = R * csc_matrix(m) * R.H  # rotate matrix
  return mout.todense() # return dense matrix





def global_spin_rotation(m,vector = np.array([0.,0.,1.]),angle = 0.0,
                             spiral = False,atoms = None):
  """ Rotates a matrix along a certain qvector """
  # pauli matrices
  from scipy.sparse import csc_matrix,bmat
  iden = csc_matrix([[1.,0.],[0.,1.]])
  n = m.shape[0]//2 # number of sites
  if atoms==None: 
    atoms = range(n) # all the atoms
  else: 
    raise
  R = [[None for i in range(n)] for j in range(n)] # rotation matrix
  from scipy.linalg import expm  # exponenciate matrix
  for i in range(n): # loop over sites
    u = np.array(vector)
    u = u/np.sqrt(u.dot(u)) # unit vector
    rot = u[0]*sx + u[1]*sy + u[2]*sz 
    # a factor 2 is taken out due to 1/2 of S
    rot = expm(np.pi*1j*rot*angle/2.0)
#    if i in atoms:
    R[i][i] = rot  # save term
#    else:
#      R[i][i] = iden  # save term
  R = bmat(R)  # convert to full sparse matrix
  if spiral:  # for spin spiral
    mout = R * csc_matrix(m)  # rotate matrix
  else:  # normal global roration
    mout = R * csc_matrix(m) * R.H  # rotate matrix
  return mout.todense() # return dense matrix



