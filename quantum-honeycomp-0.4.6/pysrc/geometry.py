import numpy as np
from copy import deepcopy
from scipy.sparse import bmat
from scipy.sparse import csc_matrix as csc

class geometry:
  """ Class for a geometry in a system """
  has_sublattice = False # has sublattice index
  dimensionality = 1 # dimension of the hamiltonian
  x = [] # positions in x
  y = [] # positions in y
  z = [] # positions in z
  r = [] # full positions 
  celldis = 1.0 # distance to the nearest cell (for 1d)
  a1 = np.array([1.0,0.0,0.])  # first vector to the nearest cell
  a2 = np.array([0.0,1.0,0.])  # first vector to the nearest cell
  shift_kspace = False # shift the klist when plotting
  name = "None"
  def plot_geometry(self):
    """Plots the system"""
    return plot_geometry(self)
  def set_finite(self):
    """ Transfrom the geometry into a finite system"""
    self.celldis = None  # create a finite system
    self.dimensionality = 0
  def supercell(self,nsuper):
    """Creates a supercell"""
    if self.dimensionality==1:
      s = supercell1d(self,nsuper)
    elif self.dimensionality==2:
      try: # two number given
        nsuper1 = nsuper[0]
        nsuper2 = nsuper[1]
      except: # one number given
        nsuper1 = nsuper
        nsuper2 = nsuper
      s = supercell2d(self,n1=nsuper1,n2=nsuper2)
    else: raise
    return s
  def xyz2r(self):
    """Updates r atributte according to xyz"""
    self.r = np.array([self.x,self.y,self.z]).transpose()
  def r2xyz(self):
    """Updates x,y,z atributtes according to r"""
    r = self.r.transpose()
    self.x = r[0]
    self.y = r[1]
    self.z = r[2]
  def get_hamiltonian(self):
    """ Create the hamiltonian for this geometry"""
    from hamiltonians import hamiltonian
    h = hamiltonian(self)  # create the object
    h.first_neighbors()  # create first neighbor hopping
    return h # return the object
  def write(self):
    """ Writes the geometry in file"""
    from input_tb90 import write_geometry
    write_geometry(self)
    from input_tb90 import write_lattice
    write_lattice(self)
  def copy(self):
    """Copy the geometry"""
    return deepcopy(self)
  def center(self):
    """ Centers the geometry in (0,0,0)"""
    self.x = self.x - sum(self.x)/len(self.x)
    self.y = self.y - sum(self.y)/len(self.y)
    self.z = self.z - sum(self.z)/len(self.z)
    self.xyz2r() # update r


def add(g1,g2):
  """ Adds two geometries """
  raise
  gs = geometry()
  gs.x = np.array(g1.x.tolist() + g2.x.tolist())  
  gs.y = np.array(g1.y.tolist() + g2.y.tolist())  
  gs.z = np.array(g1.z.tolist() + g2.z.tolist())  
  gs.celldis = max([g1.celldis,g2.celldis]) 
  return gs






def squid_square(width=4,inner_radius=6,arm_length=8,arm_width=2,fill=False):
  nt = width + inner_radius # half side of the big square
#  import geometry
  g = geometry() # create the geometry of the system
  xc = [] # empty list
  yc = [] # empty list
  shift_y = float(arm_width-1)/2.0
  for i in range(-nt,nt+1): 
    for j in range(-nt,nt+arm_width):
      # if in the ring
      yy = float(j)-shift_y  # temporal y
      if abs(i)>inner_radius or abs(yy)>(inner_radius+shift_y) or fill: 
        xc.append(float(i))  # add x coordinate
        yc.append(yy)  # add y coordinate
  # now add right and left parts
  xr = [] # empty list
  yr = [] # empty list
  xl = [] # empty list
  yl = [] # empty list
  shift_y = float(arm_width-1)/2.0
  min_x = min(xc) - 1.0
  max_x = max(xc) + 1.0
  for i in range(arm_length):
    for j in range(arm_width): # double width of the arms
      xr.append(float(i)+max_x)
      xl.append(-float(i)+min_x)
      yr.append(float(j)-shift_y)
      yl.append(float(j)-shift_y)
  x = np.array(xr+xc+xl)  # all the x positions
  y = np.array(yr+yc+yl)  # all the y positions
  g.x = x  # add positions in x
  g.y = y  # add positions in y
  g.celldis = max(x) - min(x) +1.0 # distance to neighbor +1.0
  g.dimensionality = 1 # 0 dimensional system
  return g










def honeycomb_armchair_ribbon(ntetramers):
  """ Creates the positions of an armchair ribbon
  of width ntetramers, return a geometry class """
  from numpy import array, sqrt
  n = ntetramers
  x=array([0.0 for i in range(4*n)])
  y=array([0.0 for i in range(4*n)])
  s3=sqrt(3.0)/2.0 
  for ii in range(n):
    fi=float(ii)*s3*2.0
    i=4*ii
    x[i]=0.0
    x[i+1]=1.0
    x[i+2]=1.5
    x[i+3]=2.5 
    y[i]=fi
    y[i+1]=fi
    y[i+2]=fi+s3
    y[i+3]=fi+s3
  x=x-sum(x)/float(4*n)
  y=y-sum(y)/float(4*n)
  g = geometry() # create geometry class
  g.x = x  # add to the x atribute
  g.y = y  # add to the y atribute
  g.z = y*0.0  # add to the y atribute
  g.celldis = 3.0 # add distance to the nearest cell
  g.a1 = np.array([3.0,0.,0.]) # add distance to the nearest cell
  g.shift_kspace = True # shift kpoint when plotting
  g.has_sublattice = True # has sublattice index
  g.sublattice = [(-1.)**i for i in range(len(x))] # subattice number
  g.name = "honeycomb_armchair_ribbon"  # name of the geometry
  g.xyz2r() # create r coordinates
  g.dimensionality = 1
  return g


def square_ribbon(natoms):
  """ Creates the hamiltonian of a square ribbon lattice"""
  from numpy import array
  x=array([0.0 for i in range(natoms)]) # create x coordinates
  y=array([float(i) for i in range(natoms)])  # create y coordinates
  y=y-sum(y)/float(natoms) # shift to the center
  g = geometry() # create geometry class
  g.x = x  # add to the x atribute
  g.y = y  # add to the y atribute
  g.z = y*0.0  # add to the y atribute
  g.celldis = 1.0 # add distance to the nearest cell
  g.a1 = np.array([1.0,0.,0.]) # add distance to the nearest cell
  g.xyz2r() # create r coordinates
  g.has_sublattice = False # does not have sublattice
  return g

def chain():
  """ Create a chain """
  g = square_ribbon(1) 
  g = g.supercell(2)
  g.has_sublattice = True
  g.sublattice = [1.,-1.]
  return g



def square_tetramer_ribbon(ntetramers):
  """ Creates the hamiltonian of a square tetramer ribbon lattice"""
  from numpy import array
  natoms = ntetramers*4
  x=array([0.0 for i in range(natoms)]) # create x coordinates
  y=array([0.0 for i in range(natoms)])  # create y coordinates
  for i in range(ntetramers):
    x[4*i] = 0.0
    x[4*i+1] = 1.0
    x[4*i+2] = 1.0
    x[4*i+3] = 0.0
    y[4*i] = 2.*i
    y[4*i+1] = 2.*i
    y[4*i+2] = 2.*i +1.0
    y[4*i+3] = 2.*i +1.0
  y=y-sum(y)/float(natoms) # shift to the center
  x=x-sum(x)/float(natoms) # shift to the center
  g = geometry() # create geometry class
  g.x = x  # add to the x atribute
  g.y = y  # add to the y atribute
  g.z = y*0.  # add to the z atribute
  g.celldis = 2.0 # add distance to the nearest cell
  g.a1 = np.array([2.0,0.,0.]) # add distance to the nearest cell
  g.shift_kspace = True # add distance to the nearest cell
  g.xyz2r() # create r coordinates
  g.has_sublattice = True # has sublattice index
  g.sublattice = [(-1.)**i for i in range(len(x))] # subattice number
  g.dimensionality = 1
  return g







def square_zigzag_ribbon(npairs):
  """ Creates the hamiltonian of a square zigzag (11) lattice"""
  from numpy import array,sqrt
  s2 = sqrt(2.) # square root of 2
  natoms = 2*npairs
  x=array([s2/4.*(-1)**i for i in range(natoms)]) # create x coordinates
  y=array([0.0 for i in range(natoms)])  # create y coordinates of pairs
  yp=array([s2*float(i) for i in range(npairs)])  # create y coordinates of pairs
  for i in range(npairs): # y position in each pair
    y[2*i] = yp[i]
    y[2*i+1] = yp[i] + s2/2.
  y=y-sum(y)/float(natoms) # shift to the center
  g = geometry() # create geometry class
  g.x = x  # add to the x atribute
  g.y = y  # add to the y atribute
  g.celldis = s2 # add distance to the nearest cell
  g.xyz2r() # create r coordinates
  g.dimensionality = 1
  return g




def honeycomb_zigzag_ribbon(ntetramers):
  from numpy import array, sqrt
  n = ntetramers
  x=array([0.0 for i in range(4*n)])
  y=array([0.0 for i in range(4*n)])
  s3=sqrt(3.0)/2.0
  for ii in range(n):
    fi=-float(ii)*3.0
    i=4*ii
    x[i]=0.0
    x[i+1]=s3
    x[i+2]=s3
    x[i+3]=0.0
    y[i]=fi
    y[i+1]=fi-0.5
    y[i+2]=fi-1.5
    y[i+3]=fi-2.0
  x=x-sum(x)/float(4*n)
  y=y-sum(y)/float(4*n)
  g = geometry() # create geometry class
  g.x = x  # add to the x atribute
  g.y = y  # add to the y atribute
  g.z = y*0.0  # add to the z atribute
  g.celldis = sqrt(3.0) # add distance to the neares cell
  g.a1 = np.array([sqrt(3.0),0.,0.]) # add distance to the nearest cell
  g.has_sublattice = True # has sublattice index
  g.sublattice = [(-1.)**i for i in range(len(x))] # subattice number
  g.name = "honeycomb_zigzag_ribbon"
  g.xyz2r() # create r coordinates
  g.dimensionality = 1
  return g


def plot_geometry(g):
   """Shows a 2d plot of the current geometry,
      returns a figure"""
   import pylab
   fig = pylab.figure() # create the figure
   sp = fig.add_subplot(111)
   x = g.x # x coordinates
   y = g.y # y coordinates
   sp.scatter(x,y,marker = "o",s=80,color="red") # create central cell
   celldis = g.celldis # distance to neighboring cell
   if not celldis== None: # if there is neighbor
     sp.scatter(x+celldis,y,marker = "o",s=80,color="black") # create right cell
     sp.scatter(x-celldis,y,marker = "o",s=80,color="black") # create left cell
   sp.set_xlabel("X")
   sp.set_xlabel("Y")
   sp.axis("equal") # same scale in axes
   fig.set_facecolor("white") # white figure  
   return fig


def supercell1d(g,nsuper):
  """Creates a supercell of the system"""
  # get the old geometry 
  y = g.y
  x = g.x
  z = g.z
  celldis = g.celldis
  # position of the supercell
  yout = []
  xout = []
  for i in range(nsuper):
    yout += y.tolist()
    xout += (x+i*celldis).tolist()
  # now modify the geometry
  go = deepcopy(g)
  go.x = np.array(xout)
  go.y = np.array(yout)
  # and shift to zero
  go.z = np.array(z.tolist()*nsuper)
  go.center() # center the unit cell
  go.celldis = celldis*nsuper
  go.a1 = g.a1*nsuper # supercell
  go.xyz2r() # update r
  if g.has_sublattice: # if has sublattice, keep the indexes
    go.sublattice = (np.array(go.sublattice)).tolist() # convert in list
    go.sublattice = g.sublattice*nsuper  # store the keeped atoms

  return go






################################################
########### begin 2d geometries ################
################################################

def honeycomb_lattice():
  """ Creates a honeycomb lattice """
  g = geometry() # create geometry
  g.x = np.array([-0.5,0.5])
  g.y = np.array([0.0,0.0])
  g.z = np.array([0.0,0.0])
  g.a1 = np.array([3./2.,np.sqrt(3.)/2,0.]) # first lattice vector
  g.a2 = np.array([-3./2.,np.sqrt(3.)/2,0.]) # second lattice vector
  g.dimensionality = 2 # two dimensional system
  g.xyz2r() # create r coordinates
  g.has_sublattice = True # has sublattice index
  g.sublattice = [(-1.)**i for i in range(len(g.x))] # subattice number
  return g


def square_lattice():
  """ Creates a square lattice """
  g = geometry() # create geometry
  g.x = np.array([-0.5,0.5,0.5,-0.5])
  g.y = np.array([-0.5,-0.5,0.5,0.5])
  g.z = g.x*0.
  g.a1 = np.array([2.,0.,0.]) # first lattice vector
  g.a2 = np.array([0.,2.,0.]) # second lattice vector
  g.dimensionality = 2 # two dimensional system
  g.xyz2r() # create r coordinates
  g.has_sublattice = True # has sublattice index
  g.sublattice = [(-1.)**i for i in range(len(g.x))] # subattice number
  return g



def lieb_lattice():
  """ Create a Lieb lattice"""
  g = geometry() # create geometry
  g.x = np.array([-0.5,0.5,0.5])
  g.y = np.array([-0.5,-0.5,0.5])
  g.z = g.x*0.
  g.a1 = np.array([2.,0.,0.]) # first lattice vector
  g.a2 = np.array([0.,2.,0.]) # second lattice vector
  g.dimensionality = 2 # two dimensional system
  g.xyz2r() # create r coordinates
  g.has_sublattice = True # has sublattice index
  g.sublattice = [(-1.)**i for i in range(len(g.x))] # subattice number
  return g




def kagome_lattice():
  """ Creates a honeycomb lattice """
  g = geometry() # create geometry
  dx = 1./2.
  dy = np.sqrt(3)/2.
  g.x = np.array([-dx,dx,0.])
  g.y = np.array([-dy,-dy,0.0])
  g.z = np.array([0.0,0.0,0.])
  g.a1 = np.array([2.,0.,0.]) # first lattice vector
  g.a2 = np.array([1.,np.sqrt(3),0.]) # second lattice vector
  g.dimensionality = 2 # two dimensional system
  g.xyz2r() # create r coordinates
  g.has_sublattice = False # does not have sublattice index
  return g


def rectangular_kagome_lattice():
  """ Creates a square kagome lattice"""
  g = kagome_lattice()
  g = g.supercell(2) # create a supercell
  go = g.copy()
  go.a2 = np.array([0.0,np.sqrt(3)*2.0,0.0]) 
  return go







def honeycomb_lattice_square_cell():
  """ Creates a honeycomb lattice """
  g = honeycomb_lattice() # create geometry
  go = deepcopy(g)
  go.a1 =  g.a1 + g.a2
  go.a2 = g.a1 - g.a2
  go.x = np.concatenate([g.x,g.x-g.a1[0]])  # new x coordinates
  go.y = np.concatenate([g.y,g.y-g.a1[1]])  # new y coordinates
  go.z = np.concatenate([g.z,g.z])  # new z coordinates
  go.has_sublattice = True # does not have sublattice index
  go.sublattice = [1.,-1.,1.,-1.] # sublattice index
  go.xyz2r() # create r coordinates
  return go

def honeycomb_lattice_square_cell_v2():
  """ Creates a honeycomb lattice """
  g = honeycomb_lattice() # create geometry
  go = deepcopy(g)
  go.a1 =  g.a1
  go.a2 = - g.a1 + 2*g.a2
  go.x = np.concatenate([g.x,g.x+g.a2[0]])  # new x coordinates
  go.y = np.concatenate([g.y,g.y+g.a2[1]])  # new y coordinates
  go.z = np.concatenate([g.z,g.z])  # new z coordinates
  go.xyz2r() # create r coordinates
  return go



def kagome_ribbon(n=5):
  """Create a Kagome ribbon"""
  g = rectangular_kagome_lattice() # 2d geometry
  import ribbonizate
  g = ribbonizate.bulk2ribbon(g,n=n) # create ribbon from 2d
  return g










def supercell2d(g,n1=1,n2=1):
  """ Creates a supercell for a 2d system"""
  nc = len(g.x) # number of atoms in a cell
  n = nc*n1*n2 # total number of positions
  xo = np.array([0. for i in range(n)])
  yo = np.array([0. for i in range(n)])
  ik = 0 # index of the atom
  a1 = g.a1 # first vector
  a2 = g.a2 # second vector
  for i in range(n1): 
    for j in range(n2): 
      for k in range(nc): 
        xo[ik] = i*a1[0] + j*a2[0] + g.x[k]
        yo[ik] = i*a1[1] + j*a2[1] + g.y[k]
        ik += 1
  go = deepcopy(g) # copy geometry
  xo = xo - sum(xo)/len(xo) # shift x
  yo = yo - sum(yo)/len(yo) # shift y
  go.x = xo  
  go.y = yo  
  go.z = np.array((n1*n2)*((go.z).tolist()))
  go.a1 = a1*n1
  go.a2 = a2*n2
  # shift to zero
  go.center() 
  if g.has_sublattice: # supercell sublattice
    go.sublattice = g.sublattice*(n1*n2)
  return go



def read(input_file="POSITIONS.OUT"):
  """ Reads a geometry """
  m = np.genfromtxt(input_file).transpose()
  g = geometry() # cretae geometry
  g.dimensionality = 0
  g.x = m[0]
  g.y = m[1]
  g.x = g.x - sum(g.x)/len(g.x) # normalize
  g.y = g.y - sum(g.y)/len(g.y) # normalize
  g.z = m[2]
  g.xyz2r() # create r coordinates
  try:   lat = np.genfromtxt("LATTICE.OUT")   # read lattice
  except: lat = None
  try: # two dimensional
    g.a1 = np.array([lat[0,0],lat[0,1],0.0])
    g.a2 = np.array([lat[1,0],lat[1,1],0.0])
    g.dimensionality=2
    return g
  except: pass
  try:
    g.celldis = lat
    g.dimensionality = 1
    return g
  except: pass
  g.dimensionality = 0
  return g




def bulk2ribbon(h,n=10):
  """Converts a hexagonal bulk hamiltonian into a ribbon hamiltonian"""
  from hexagonal import bulk2ribbon as br
  return br(h,n=n)


def get_reciprocal2d(a1,a2):
  """Get reciprocal vectors"""
  b1 = np.array([a2[1],-a2[0],0.])
  b2 = np.array([-a1[1],a1[0],0.])
  b1 = b1 / np.sqrt(b1.dot(b1))
  b2 = b2 / np.sqrt(b2.dot(b2))
  return (b1,b2)


def apilate(g,drs=[np.array([0.,0.,0.])]):
  """ generate a geometry by appilating another geometry, displacements
  are given on input """
  nl = len(drs) # number of layers
  ro = np.concatenate([g.r + dr for dr in drs])
  go = g.copy()
  go.r = ro
  go.r2xyz()
  return go





