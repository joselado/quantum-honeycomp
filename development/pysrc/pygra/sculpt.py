from __future__ import print_function
from . import geometry
from copy import deepcopy
import numpy as np
try: 
  import clean_geometryf90
  use_fortran = True
except: use_fortran = False


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
  iis = []
  for (ii,ix,iy,iz) in zip(range(len(g.x)),g.x,g.y,g.z): # loop over positions
    if f([ix,iy,iz]): # if the function yields true
      x.append(ix)
      y.append(iy)
      z.append(iz)
      iis.append(ii)
  gout.x = np.array(x)  # copy x
  gout.y = np.array(y)  # copy y
  gout.z = np.array(z)  # copy y
  gout.xyz2r() # update r
  if gout.has_sublattice: # if has sublattice, keep the indexes
    gout.sublattice = [g.sublattice[i] for i in iis]
  return gout






def intersected_indexes(g,f):
  """Return the indexes of the atoms located in this function"""
  iis = []
  for (ii,ix,iy,iz) in zip(range(len(g.x)),g.x,g.y,g.z): # loop over positions
    if f([ix,iy,iz]): # if the function yields true
      iis.append(ii)
  return iis











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
  if np.abs(angle)<0.0000001: 
    print("No rotation performed")
    return g
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
  elif go.dimensionality==1: 
    x,y,z = go.a1
    go.a1 = np.array([c*x + s*y,-s*x + c*y,z])
  elif go.dimensionality==0: pass
  else: raise # 
  go.get_fractional() # get fractional coordinates 
  return go


def center(g,angle):
  """Center a geometry"""
  g.x = g.x -sum(g.x)/len(g.x)
  g.y = g.y -sum(g.y)/len(g.y)
 


def remove_unibonded(g,d=1.0,tol=0.01,use_fortran=use_fortran,iterative=False):
  """Removes from the geometry atoms with only one bond"""
  sb = []
  if g.dimensionality>0: use_fortran = False # only for 0d
  if use_fortran: # use fortran routine
# array with true in right atoms
    retain = clean_geometryf90.clean_geometry(g.r.transpose(),2) 
    sb = [] # remove this atoms
    for i in range(len(retain)):
      if not retain[i]: sb.append(i) # remove
    gout = remove(g,sb) # remove those atoms
  else: # use python routine
    for i in range(len(g.r)): 
      r1 = g.r[i] # first position
      nb = 0 # initialize
      for direc in g.neighbor_directions(): # loop over directions
        for r2 in g.replicas(d=direc): # loop over replicas
          dr = r1-r2
          if d-tol < dr.dot(dr) < d+tol: # if sirdt neighbor
            nb += 1 # increase counter
      if nb<2:
        sb.append(i+0) # add to the list
    gout = remove(g,sb) # remove those atoms
  if iterative: # ensure that it hs the same number of atoms by calling again
    if len(g.x) != len(gout.x): # call again
      print("Iterative cleaning")
      return remove_unibonded(gout,d=d,tol=tol,
         use_fortran=use_fortran,iterative=iterative)
  return gout # return the geometry



def remove_central(g,n):
  """Removes n atoms from the center of the crystal"""
  rr = [r.dot(r) for r in g.r] # norm of the distances
  inds = range(len(rr)) # indexes
  sort_inds = [x for (y,x) in sorted(zip(rr,inds))] # indexes sorted by distance
  rind = [sort_inds[i] for i in range(n)]
  return remove(g,rind)  # return geometry with removed 


def get_closest(g,n=1,r0=[0.,0.,0.]):
  """Gets n atoms from the center of the crystal"""
  r0 = np.array(r0)
  rr = [(r-r0).dot(r-r0) for r in g.r] # norm of the distances
  inds = range(len(rr)) # indexes
  sort_inds = [x for (y,x) in sorted(zip(rr,inds))] # indexes sorted by distance
  rind = [sort_inds[i] for i in range(n)]
  return rind # return the indexes


def get_central(g,n=1):
  """Get the index of the central atom"""
  return get_closest(g,n=n,r0=[0.,0.,0.])

def shift(g,r=np.array([0.,0.,0.])):
  """Shift the geometry by a certain vector"""
  g.r = np.array([ri-r for ri in g.r])
  g.r2xyz() # update


def get_angle(v1,v2):
  """Get the angle between two vectors"""
  v3 = v1/np.sqrt(v1.dot(v1)) # normalize
  v4 = v2/np.sqrt(v2.dot(v2)) # normalize
  alpha = np.arccos(v3.dot(v4))
  return alpha





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
  # define two complex numbers
  za = da[0] + da[1]*1j
  zb = db[0] + db[1]*1j
  angle = np.angle(za/zb) # angle in the complex plane
  return rotate(g,angle)


def build_island(gin,n=5,angle=20,nedges=6,clear=True):
  """ Build an island starting from a 2d geometry"""
  nf = float(n)   # get the desired size, in float
  if gin.dimensionality!=2: raise 
  g = gin.copy()
  g = g.supercell(8*n)   # create supercell
  g.set_finite() # set as finite system
  g.center() # center the geometry
  # now scuplt the geometry
  g = rotate(g,angle*2.*np.pi/360) # initial rotation
  def f(x,y): return x>-nf*(np.cos(np.pi/3)+1.)  # function to use as cut
  for i in range(nedges): # loop over rotations, 60 degrees
    g = intersec(g,f) # retain certain atoms
    g = rotate(g,2.*np.pi/nedges) # rotate 60 degrees
  if clear:  g = remove_unibonded(g)  # remove single bonded atoms
  g.center() # center the geometry
  return g # return the new geometry





def reciprocal(v1,v2,v3=np.array([0.,0.,1.])):
  """Return the reciprocal vectors"""
  vol = v1.dot(np.cross(v2,v3)) # volume
  w1 = np.cross(v2,v3)/vol
  w2 = np.cross(v3,v1)/vol
  w3 = np.cross(v1,v2)/vol
  return (w1,w2,w3)



def build_ribbon(g,n):
  """ Return a geometry of a ribbon based on this cell"""
  if g.dimensionality!=2: raise # if it is not two dimensional
  angle = sculpt.get_angle(g.a1,g.a2)/np.pi*180 # get the angle
  if np.abs(angle-90)<1.: # if it is square
    gout = g.copy() # copy geometry
    gout.dimensionality = 1
    rs = []
    for ir in g.r:
      for i in range(n):
        rs.append(ir+g.a1*i) # append position
    gout.r = rs
    gout.r2xyz() # update
    raise
    return gout


def image2island(impath,g,nsuper=4,size=10,color="black"):
  """Build an island using a certain image"""
  from PIL import Image
  im = Image.open(impath)
  im = im.convert('RGBA')
  data = np.array(im) # convert to array
  red, green, blue, alpha = data.T # store data for readbility
  if color=="black": #retain the black color
    retain = (red < 20) & (blue < 20) & (green < 20)
  elif color=="red": #retain the black color
    retain = (red > 200) & (blue < 20) & (green < 20)
  elif color=="blue": #retain the black color
    retain = (red < 20) & (blue > 200) & (green < 20)
  elif color=="green": #retain the black color
    retain = (red < 20) & (blue < 20) & (green > 200)
  else: raise # unrecognized
  data[..., :-1][retain.T] = (0, 0, 0) # set as black
  data[..., :-1][np.logical_not(retain.T)] = (255, 255, 255) # set as white
#  data[..., :-1][not retain.T] = (255, 255, 255) # set as black
  im2 = Image.fromarray(data) # convert to image
  im2 = im2.convert("L") # to black and white
  bw = np.asarray(im2).copy() # convert to array
  bw[bw < 128] = 0  # Black
  bw[bw >= 128] = 1 # White
  bw = bw.transpose() # transpose image
  # now create a supercell
  nx,ny = bw.shape # size of the image
  go = g.supercell(nsuper*size) # build supercell
  go.center()
  minx = -size
  maxx = size
  miny = -size*bw.shape[1]/bw.shape[0]
  maxy = size*bw.shape[1]/bw.shape[0]
  def finter(rtmp):
    x = rtmp[0]
    y = rtmp[1]
    x = (x - minx)/(maxx-minx)
    y = (y - miny)/(maxy-miny)
    xi = (nx-1)*x # normalized
    yi = (ny-1)*y # normalized
    xi,yi = int(round(xi)),int(round(yi)) # integer
    if not 0<xi<bw.shape[0]: return False
    if not 0<yi<bw.shape[1]: return False
    if bw[xi,yi]==0: return True
    else: return False
  go = intersec(go,finter)
  go.dimensionality = 0 # zero dimensional
  go.celldis = None
  return go 


def common(g1,g2,tol=0.1):
  """Return the indexes of atoms common in both structures"""
  indexes = []
  for i in range(len(g1.r)):
    r1 = g1.r[i] # position
    for r2 in g2.r:
      dr = r1-r2
      if dr.dot(dr)<tol:
        indexes.append(i)
        break # next iteration
  return indexes

def add(g1,g2):
  g = g1.copy() # copy geometry
  g.x = np.concatenate([g1.x,g2.x])
  g.y = np.concatenate([g1.y,g2.y])
  g.z = np.concatenate([g1.z,g2.z])
  g.xyz2r()
  g.has_sublattice = False
  return g


def set_xy_plane(g):
  """Modify a geometry so the lattice vectors lie in the xy plane"""
#  if g.dimensionality != 2: raise # only for 2d
  go = g.copy() # copy geometry
  nv = np.cross(g.a1,g.a2)
  nv /= np.sqrt(nv.dot(nv)) # unitary normal vector
  rho = np.sqrt(nv[0]*nv[0] + nv[1]*nv[1]) # planar component
  theta = np.arctan2(rho,nv[2]) # theta angle
  phi = np.arctan2(nv[1],nv[0]) # phi angle
  # matrix that transforms (0,0,1) to that vector
  ct,st = np.cos(theta),np.sin(theta)
  cp,sp = np.cos(phi),np.sin(phi)
  Rt = np.matrix([[ct,0,st],[0.,1.,0],[-st,0,ct]]) # rotate along y
  Rp = np.matrix([[cp,-sp,0.],[sp,cp,0],[0.,0.,1.]]) # rotate along z
  R = Rp*Rt # transforms (0,0,1) to nv
  U = R.I # inverse transformation
  # now transform everything
  def transform(r): return np.array((U*np.matrix(r).T).T)[0] 
  go.a1 = transform(g.a1)
  go.a2 = transform(g.a2)
  go.a3 = transform(g.a3)
  go.r = np.array([transform(ri) for ri in g.r])
  go.r2xyz()
  return go




def retain_unit_cell(r,a1,a2,a3,dim=3):
  """Retain position located in the unit cell defined by a1,a2,a3"""
  R = np.matrix([a1,a2,a3]).T # transformation matrix
  L = R.I # inverse matrix
  rs = [] # output
#  d0 = -np.random.random()*0.001 - .5 # accuracy
  d0 = 0.00234231421 - 0.5 # random number
  d1 = 1.0 + d0 # accuracy

  for ri in r: # loop over positions
    rn = L*np.matrix(ri).T  # transform
    rn = np.array([rn.T[0,i] for i in range(3)]) # convert to array
    n1,n2,n3 = rn[0],rn[1],rn[2]
    if dim==3:
      if d0<n1<d1 and d0<n2<d1 and d0<n3<d1:
         rs.append(ri)
    if dim==2:
      if d0<n1<d1 and d0<n2<d1:
         rs.append(ri)
  return np.array(rs) # return positions








