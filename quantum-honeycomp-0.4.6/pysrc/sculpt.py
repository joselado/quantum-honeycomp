import geometry
from copy import deepcopy
import numpy as np


def remove(g,l):
  """ Remove certain atoms from the geometry"""
  nr = len(l) # number of removed atoms
  go = g.copy() # copy the geometry
  xo = [] # copy the list
  yo = [] # copy the list
  zo = [] # copy the list
  for i in range(len(g.x)):
    if not i in l:
      xo.append(g.x[i])
      yo.append(g.y[i])
      zo.append(g.z[i])
  go.x = np.array(xo)
  go.y = np.array(yo)
  go.z = np.array(zo)
  go.xyz2r() # update the revectors
  ##### if has sublattice ####
  if g.has_sublattice: # if has sublattice, keep the indexes
    ab = [] # initialize
    for i in range(len(g.x)):
      if not i in l:
        ab.append(g.sublattice[i]) # keep the index
    go.sublattice = ab # store the keeped atoms
  return go

def intersec(g,f):
  """ Intersec coordinates with a certain function which yields True or False,
  output is resultant geometry """
  gout = g.copy() # copy the geometry
  x = [] # out x
  y = [] # out y
  z = [] # out y
  for (ix,iy,iz) in zip(g.x,g.y,g.z): # loop over positions
    if f(ix,iy): # if the function yields true
      x.append(ix)
      y.append(iy)
      z.append(iz)
  gout.x = np.array(x)  # copy x
  gout.y = np.array(y)  # copy y
  gout.z = np.array(z)  # copy y
  gout.xyz2r() # update r
  return gout


def circle(r=1.0,out=False):
  """ Returns a function which encondes a circle"""
  if out:  # if true is inside
    def f(x,y):
      if r*r > x*x + y*y:
        return True
      else:
        return False  
  else: # if true is outside
    def f(x,y):
      if r*r < x*x + y*y:
        return True
      else:
        return False  
  return f # return the function



def rotate(g,angle):
  """ Rotates a geometry"""
#  phi = np.pi*angle # angle in radians
  phi = angle
  go = g.copy()
  # modify x and y, z is the same
  x,y,z = [],[],[]  # initialize lists
  c,s = np.cos(phi), np.sin(phi)  # sin and cos of the anggle
  for (ix,iy,iz) in zip(g.x,g.y,g.z):
    x.append(c*ix + s*iy)    # x coordinate 
    y.append(-s*ix + c*iy)    # y coordinate
    z.append(iz)
  go.x = np.array(x)
  go.y = np.array(y)
  go.z = np.array(z)
  go.xyz2r() # update r
  if go.dimensionality==2:  # two dimensional
    x,y,z = go.a1
    go.a1 = np.array([c*x + s*y,-s*x + c*y,z])
    x,y,z = go.a2
    go.a2 = np.array([c*x + s*y,-s*x + c*y,z])
  elif go.dimensionality==1: raise
  elif go.dimensionality==0: pass
  return go


def center(g,angle):
  """Center a geometry"""
  g.x = g.x -sum(g.x)/len(g.x)
  g.y = g.y -sum(g.y)/len(g.y)
 


def remove_unibonded(g,d=1.0,tol=0.01):
  """Removes from the geometry atoms with only one bond"""
  sb = []
  for i in range(len(g.r)): 
    r1 = g.r[i] # first position
    nb = 0 # initialize
    for r2 in g.r:
      dr = r1-r2
      if d-tol < dr.dot(dr) < d+tol: # if sirdt neighbor
        nb += 1 # increase counter
    if nb<2:
      sb.append(i+0) # add to the list
  return remove(g,sb) # remove those atoms


def remove_central(g,n):
  """Removes n atoms from the center of the crystal"""
  rr = [r.dot(r) for r in g.r] # norm of the distances
  inds = range(len(rr)) # indexes
  sort_inds = [x for (y,x) in sorted(zip(rr,inds))] # indexes sorted by distance
  rind = [sort_inds[i] for i in range(n)]
#  print "remove",rind
  return remove(g,rind)  # return geometry with removed 


def get_central(g,n=1):
  """Gets n atoms from the center of the crystal"""
  rr = [r.dot(r) for r in g.r] # norm of the distances
  inds = range(len(rr)) # indexes
  sort_inds = [x for (y,x) in sorted(zip(rr,inds))] # indexes sorted by distance
  rind = [sort_inds[i] for i in range(n)]
  return rind # return the indexes



def get_furthest(g,n=1,angle=0.,tol=5.):
  """Gets n atoms in a certain direction"""
  rs = [] # norm of the distances
  inds = [] # norm of the distances
  for ir in range(len(g.r)): # store only vectors with a certain angle
    r = g.r[ir]
    a = np.arctan2(r[1],r[0])/np.pi*180.
    if (np.abs(angle-a)%360)<tol: # if direction is right
      rs.append(r) # store vector
      inds.append(ir) # indexes
  rr = [-r.dot(r) for r in rs] # norm of the distances
  sort_inds = [x for (y,x) in sorted(zip(rr,inds))] # indexes sorted by distance
  rind = [sort_inds[i] for i in range(n)]
  return rind # return the indexes





def rotate_a2b(g,a,b):
  """ Rotates the geometry making a original vector a pointing along b"""
  da = a.dot(a)
  da = a/np.sqrt(da) # unit vector
  db = b.dot(b)
  db = b/np.sqrt(db) # unit vector
  angle = da.dot(db) # angle to rotate
  angle = np.arccos(angle)
  return rotate(g,angle)



