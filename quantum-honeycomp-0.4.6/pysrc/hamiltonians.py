import numpy as np
from scipy.sparse import csc_matrix,bmat
import multicell
import scipy.linalg as lg

optimal = False

class hamiltonian():
  """ Class for a hamiltonian """
  has_spin = True # has spin degree of freedom
  has_eh = False # has electron hole pairs
  fermi_energy = 0.0 # fermi energy at zero
  is_sparse = False
  is_multicell = False # for hamiltonians with hoppings to several neighbors
  def eigenvectors(self,nk=10):
    return eigenvectors(self,nk=nk)
  def __init__(self,geometry=None):
    if not geometry == None:
# dimensionality of the system
      self.dimensionality = geometry.dimensionality 
      self.geometry = geometry # add geometry object
      self.num_orbitals = len(geometry.x)
  def get_hk_gen(self):
    """ Generate kdependent hamiltonian"""
    if self.is_multicell: return multicell.hk_gen(self) # for multicell
    else: return hk_gen(self) # for normal cells
  def print_hamiltonian(self):
    """Print hamiltonian on screen"""
    print_hamiltonian(self)
  def diagonalize(self,nkpoints=100):
    """Return eigenvalues"""
    return diagonalize(self,nkpoints=nkpoints)
  def get_bands(self,nkpoints=100,use_lines=False,kpath=None,operator=None):
    """ Returns a figure with teh bandstructure"""
    if self.dimensionality==0: # for 0d and 1d system
      get_bands(self,nkpoints=nkpoints)
    if self.dimensionality==1: # for 0d and 1d system
      get_bands1d(self,nkpoints=nkpoints,operator=operator)
    if self.dimensionality==2:  # for 2d system
      get_bands_2d(self,kpath=kpath,operator=operator)
  def plot_bands(self,nkpoints=100,use_lines=False):
    """ Returns a figure with teh bandstructure"""
    return plot_bands(self,nkpoints=nkpoints,use_lines=use_lines)
  def add_sublattice_imbalance(self,mass):
    """ Adds a sublattice imbalance """
    add_sublattice_imbalance(self,mass)
  def add_antiferromagnetism(self,mass,d="z"):
    """ Adds antiferromagnetic imbalanc """
    add_antiferromagnetism(self,mass,d=d)
  def add_pwave_electron_hole_pairing(self,delta=0.0,mu=0.0,phi=0.0):
    """ Adds spin conserving first neighbor electron hole pairing"""
    (self.intra,self.inter) = add_pwave_electron_hole_pairing(self,
                                delta=delta,mu=mu,phi=phi)
  def add_swave_electron_hole_pairing(self,delta=0.0,mu=0.0,phi=0.0):
    """ Adds spin mixing insite electron hole pairing"""
    (self.intra,self.inter) = add_swave_electron_hole_pairing(self,
                                delta=delta,mu=mu,phi=phi)
  def add_swave(self,delta=0.0,mu=0.0,phi=0.0):
    """ Adds spin mixing insite electron hole pairing"""
   # self.intra = add_swave(self.intra,delta=delta,mu=mu,phi=phi)
    (self.intra,a) = add_swave_electron_hole_pairing(self,
                                delta=delta,mu=mu,phi=phi)
    if self.dimensionality==1: # one dimensional systems
      self.inter = nambu(self.inter)
    if self.dimensionality==2: # two dimensional systems
      self.tx = nambu(self.tx)
      self.ty = nambu(self.ty)
      self.txy = nambu(self.txy)
      self.txmy = nambu(self.txmy)
    self.has_eh = True
  def supercell(self,nsuper):
    """ Creates a supercell of a one dimensional system"""
    import supercell
    if self.dimensionality==1:
      hout = supercell.supercell1d(self,nsuper,sparse=self.is_sparse) # modify hoppings
    return hout
  def set_finite_system(self):
    """ Transforms the system into a finite system"""
    set_finite_system(self) 
  def write(self,output_file="hamiltonian.in"):
    """ Write the hamiltonian in hamiltonian_0.in"""
    from input_tb90 import write_hamiltonian
    write_hamiltonian(self,output_file=output_file)
    from input_tb90 import write_geometry
    write_geometry(self.geometry)
    from input_tb90 import write_lattice
    write_lattice(self.geometry)
  def read(self,input_file="hamiltonian.in"):
    """ Write the hamiltonian in hamiltonian_0.in"""
    from input_tb90 import read_hamiltonian
    read_hamiltonian(self,input_file=input_file)
  def total_energy(self,nkpoints=100):
    """ Get total energy of the system"""
    return total_energy(self,nkpoints=nkpoints)
  def add_zeeman(self,zeeman):
    """Adds zeeman to the matrix """
    if self.has_spin:  # if it has spin degree of freedom
      add_zeeman(self,zeeman=zeeman)
  def turn_spinful(self):
    """Turn the hamiltonian spinful""" 
    if self.is_multicell: # if multicell
      from multicell import turn_spinful as ts
      ts(self) # turn spinful
    else:
      from increase_hilbert import spinful
      if not self.has_spin:
        self.has_spin = True
        self.intra = spinful(self.intra) 
        if self.dimensionality==1:
          self.inter = spinful(self.inter) 
        if self.dimensionality==2:
          self.tx = spinful(self.tx) 
          self.ty = spinful(self.ty) 
          self.txy = spinful(self.txy) 
          self.txmy = spinful(self.txmy) 
  def remove_spin(self):
    """Removes spin degree of freedom"""
    if self.has_spin: # if has spin remove the component
      self.intra = des_spin(self.intra,component=0)
      if self.dimensionality==1: # if one dimensional
        self.inter = des_spin(self.inter,component=0)
      if self.dimensionality==2: # if one dimensional
        self.tx = des_spin(self.tx,component=0)
        self.txy = des_spin(self.txy,component=0)
        self.ty = des_spin(self.ty,component=0)
        self.txmy = des_spin(self.txmy,component=0)
    self.has_spin = False  # flag for nonspin calculation
  def shift_fermi(self,fermi):
    """ Move the Fermi energy of the system"""
    shift_fermi(self,fermi)
  def first_neighbors(self):
    """ Create first neighbor hopping"""
    if self.dimensionality == 0:
      first_neighbors0d(self)
    if self.dimensionality == 1:
      first_neighbors1d(self)
    if self.dimensionality == 2:
      first_neighbors2d(self)
  def dos_semiinfinite(self,energies=[0.0],num_rep=1000,
                      mixing=0.7,eps=0.0001,green_guess=None,max_error=0.0001):
    """ Calculates the surface density of states by using a 
    green function approach"""
    if self.dimensionality==1:  # if one dimensional
      from green import dos_semiinfinite
      dos = dos_semiinfinite(self.intra,self.inter,energies=energies,
             num_rep=num_rep,mixing=mixing,
             eps=eps,green_guess=green_guess,max_error=max_error)
      return dos # return the density of states
  def dos_infinite(self,energies=[0.0],num_rep=1000,
                      mixing=0.7,eps=0.0001,green_guess=None,max_error=0.0001):
    """ Calculates the density of states by using a 
    green function approach"""
    if self.dimensionality==1:  # if one dimensional
      from green import dos_infinite
      dos = dos_infinite(self.intra,self.inter,energies=energies,
             num_rep=num_rep,mixing=mixing,
             eps=eps,green_guess=green_guess,max_error=max_error)
      return dos # return the density of states

  def plot_dos_semiinfinite(self,energies=[0.0],num_rep=1000,
                      mixing=0.7,eps=0.0001,green_guess=None,max_error=0.0001):
    """ Plots the surface density of states by using a 
    green function approach"""
    if self.dimensionality==1:  # if one dimensional
      from green import plot_dos_semiinfinite
      fig = plot_dos_semiinfinite(self.intra,self.inter,energies=energies,
             num_rep=num_rep,mixing=mixing,
             eps=eps,green_guess=green_guess,max_error=max_error)
      return fig


  def plot_dos_infinite(self,energies=[0.0],num_rep=1000,
                      mixing=0.7,eps=0.0001,green_guess=None,max_error=0.0001):
    """ Plots the surface density of states by using a 
    green function approach"""
    if self.dimensionality==1:  # if one dimensional
      from green import plot_dos_infinite
      fig = plot_dos_infinite(self.intra,self.inter,energies=energies,
             num_rep=num_rep,mixing=mixing,
             eps=eps,green_guess=green_guess,max_error=max_error)
      return fig
  def copy(self):
    """ Returns a copy of the hamiltonian"""
    from copy import deepcopy
    return deepcopy(self)
  def check(self):
    """Checks if the Hamiltonian is hermitic"""
    hh = self.intra -self.intra.H
    if np.max(np.abs(hh)) > 0.000001:
      print "Warning, non hermitic hamiltonian"
      self.intra = (self.intra + self.intra.H)/2.0      
  def to_sparse(self):
    """ Transforms the hamiltonian into a sparse hamiltonian"""
    from scipy.sparse import csc_matrix
    self.intra = csc_matrix(self.intra)
    if self.dimensionality == 1: # for one dimensional
      self.inter = csc_matrix(self.inter)
  def add_rashba(self,c):
    """Adds Rashba coupling"""
    g = self.geometry
    self.intra += rashba(g.r,c=c)
    if g.dimensionality==1:  # one dimensional
      r2 = [ir + g.a1 for ir in g.r]
      self.inter += rashba(g.r,r2=r2,c=c)
    if g.dimensionality==2:  # two dimensional
      r2 = [ir + g.a1 for ir in g.r]
      self.tx += rashba(g.r,r2=r2,c=c)
      r2 = [ir + g.a2 for ir in g.r]
      self.ty += rashba(g.r,r2=r2,c=c)
      r2 = [ir + g.a1+g.a2 for ir in g.r]
      self.txy += rashba(g.r,r2=r2,c=c)
      r2 = [ir + g.a1-g.a2 for ir in g.r]
      self.txmy += rashba(g.r,r2=r2,c=c)
  def add_kane_mele(self,t):
    """ Adds a Kane-Mele SOC term"""  
    g = self.geometry
    if not self.has_spin: raise  # only for spinfull
#    x , y = g.x ,g.y  
#    if self.dimensionality==0:  # zero dimensional
#      rs = [[x,y]]
#    if h.dimensionality==1:  # one dimensional
#      rs = [[x,y],[x+g.celldis,y]]
    if self.dimensionality==0:  # two dimensional
      self.intra += t*kane_mele(g.x,g.y,x2=g.x,y2=g.y,x3=g.x,y3=g.y)
    elif self.dimensionality==2:  # two dimensional
      rs = [] # all the cells
      for i in (-1,0,1): # loop over neighbouring cells
        for j in (-1,0,1):
          x = g.x + i*g.a1[0] + j*g.a2[0]
          y = g.y + i*g.a1[1] + j*g.a2[1]
          rs.append([x,y]) # add to the array
      tlist = [(self.intra,(0,0))]
      tlist += [(self.tx,(1,0))]
      tlist += [(self.ty,(0,1))]
      tlist += [(self.txy,(1,1))]
      tlist += [(self.txmy,(1,-1))]
      for (hop,(i,j)) in tlist: # loop over hoppings
        dx = i*g.a1[0] + j*g.a2[0]
        dy = i*g.a1[1] + j*g.a2[1]
        for (xm,ym) in rs:
          hop += t*kane_mele(g.x,g.y,x2=g.x+dx,y2=g.y+dy,x3=xm,y3=ym)

#    self.intra += t*kane_mele(x,y)
#    if h.dimensionality < 
    elif self.dimensionality==1:  # one dimensional
      (intra_soc, inter_soc) = tmp_kane_mele(g.x,g.y,celldis = g.celldis,
                                          lambda_soc = t)
      self.intra += intra_soc
      self.inter += inter_soc
    else: raise
  def add_peierls(self,mag_field):
    """Add magnetic field"""
    add_peierls(self,mag_field=mag_field)
  def align_magnetism(self,vectors):
    """ Rotate the Hamiltonian to have magnetism in the z direction"""
    from rotate_spin import align_magnetism as align
    self.intra = align(self.intra,vectors)
    self.inter = align(self.inter,vectors)
  def global_spin_rotation(self,vector=np.array([0.,0.,1.]),angle=0.):
    from rotate_spin import global_spin_rotation as gsr
    self.intra = gsr(self.intra,vector=vector,angle=angle)
    self.inter = gsr(self.inter,vector=vector,angle=angle)
  def generate_spin_spiral(self,vector=np.array([1.,0.,0.]),
                            angle=0.,atoms=None):
    from rotate_spin import global_spin_rotation as gsr
    self.inter = gsr(self.inter,vector=vector,angle=angle,
                     spiral=True,atoms=atoms)
  def add_transverse_efield(self,efield=0.0):
    """Adds a transverse electric field"""
    def f(x,y):  # define function which shifts fermi energy
      return y*efield
    self.shift_fermi(f) # shift fermi energy locally







from scipy.sparse import coo_matrix,bmat
sx = coo_matrix([[0.,1.],[1.,0.]])
sy = coo_matrix([[0.,-1j],[1j,0.]])
sz = coo_matrix([[1.,0.],[0.,-1.]])





def rashba(r1,r2=None,c=0.,d=[0.,0.,1.]):
  """Add Rashba coupling, returns a spin polarized matrix"""
  zero = coo_matrix([[0.,0.],[0.,0.]])
  if r2==None:
    r2 = r1
  nat = len(r1) # number of atoms
  m = [[zero for i in range(nat)] for j in range(nat)] # cretae amtrix
  for i in range(nat): # loop over first atoms
    for j in range(nat):  # loop over second atoms
      rij = r2[j] - r1[i]   # x component
      dx,dy,dz = rij[0],rij[1],rij[2]  # different components
      if 0.9<(rij.dot(rij))<1.1: # if nearest neighbor
        rxs = [dy*sz-dz*sy,dz*sx-dx*sz,dx*sy-dy*sx]  # cross product
        ms = 1j*(d[0]*rxs[0] + d[1]*rxs[1] + d[2]*rxs[2]) # E dot r times s
        if callable(c): # call the coupling strength
          s = ms*c(r1,r2) 
        else:
          s = ms*c # multiply
        m[i][j] = s # rashba term
  m = bmat(m).todense()  # to normal matrix
  return m


# routine to obtain haldane hopping to second neighbors
def haldane(x1,y1,x2=None,y2=None,xm=None,ym=None):
  from scipy.sparse import coo_matrix, bmat
  from numpy import array, matrix, sqrt
  # use x1,y1 as default if the others not provided
  if x2==None:
    x2=x1
  if y2==None:
    y2=y1
  if xm==None:
    xm=x1
  if ym==None:
    ym=y1
  # do all the stuff
  s3=sqrt(3.0)
  n=len(x1)
  mat=[[0.0j for i in range(n)] for j in range(n)]
  for i in range(n):
    r1=array([x1[i],y1[i]])
    for j in range(n):
      r2=array([x2[j],y2[j]])
# if it is second neighbor...
      if is_neigh(r1,r2,s3,0.1):
# find the intermediate atom
        for k in range(n):
          rm=array([xm[k],ym[k]])
          if is_neigh(r1,rm,1.0,0.1):
            if is_neigh(r2,rm,1.0,0.1):
# and write the chiral dependent hopping
              mat[i][j]=1j*vec_chi(r1-rm,r2-rm)
  mat=matrix(mat)
  return mat











def tmp_kane_mele(x,y,first = 1.0, celldis = 3.0,lambda_soc=0.0,mag_field=0.0):
  """ Computes the Kane-Mele matrix,
      input are x and y coordinates,
      distance to first neighbors,
      distance to nearest cell,
      lambda_soc is SOC strength"""
  from increase_hilbert import m2spin
  from numpy import array
  def tmp_haldane(value,minus=False): # calculate Haldane hopping
    onsite_hal = get_hal_hop_sec_neigh(x,y,x,y,x,y,
                 mag_field=mag_field,value=value)
    onsite_hal += get_hal_hop_sec_neigh(x,y,x,y,x+celldis,
                 y,mag_field=mag_field,value=value)
    onsite_hal += get_hal_hop_sec_neigh(x,y,x,y,x-celldis,
                 y,mag_field=mag_field,value=value)
# check hermiticity if intra_hal
    if not is_hermitic(onsite_hal):
      print "onsite_hal is not hermitic"
      raise
# hopping matrix
    hopping_hal = get_hal_hop_sec_neigh(x,y,x+celldis,y,x,y,
               mag_field=mag_field,value=value)
    hopping_hal += get_hal_hop_sec_neigh(x,y,x+celldis,y,
               x+celldis,y,mag_field=mag_field,value=value)
    hopping_hal += get_hal_hop_sec_neigh(x,y,x+celldis,y,
               x-celldis,y,mag_field=mag_field,value=value)
# total hopping and onsite
    if not minus:
      intra_hal = onsite_hal
      inter_hal = hopping_hal
    if minus:
      intra_hal = -onsite_hal
      inter_hal = -hopping_hal
    return (intra_hal,inter_hal)  # return haldane hoppings
  # if there is SOC
  if lambda_soc == 0.0:
    return (0.0,0.0)  # don't do anything
  (intra_up,inter_up) = tmp_haldane(lambda_soc,minus=False)  # up component
  (intra_dn,inter_dn) = tmp_haldane(lambda_soc,minus=True)  # up component
# add spin degree of freedom
  intra = m2spin(intra_up,intra_dn) # intracell term 
  inter = m2spin(inter_up,inter_dn) # intercell term
  return (intra,inter)




def kane_mele(x1,y1,x2=None,y2=None,x3=None,y3=None,has_spin=True):
  """ Generate Kane Mele matrices"""
  xs,ys = [x2,x3],[y2,y3]
  for ix in xs: 
    if ix==None: 
      ix=x1
  for iy in ys: 
    if iy==None: 
      iy=y1
  m = get_hal_hop_sec_neigh(x1,y1,x2,y2,x3,y3)
  from increase_hilbert import m2spin
  if has_spin:
    m = m2spin(m,np.conjugate(m))
  return m



def find_first_neighbor(x1,y1,x2,y2,optimal=False):
  """ Gets the first neighbors, returns a list of tuples """
  if optimal:
    import neighbor
    pairs = neighbor.find_first_neighbor(x1,y1,x2,y2)
  else:
    from numpy import array
    n=len(x1)
    pairs = [] # pairs of neighbors
    for i in range(n):
      r1=array([x1[i],y1[i]])
      for j in range(n):
        r2=array([x2[j],y2[j]])
        if is_neigh(r1,r2,1.0,0.1): # check if distance is 1
          pairs.append([i,j])  # add to the list
  return pairs # return pairs of first neighbors




def get_first_neighbors(r1,r2,optimal=optimal):
  """Gets the fist neighbors, input are arrays"""
  if optimal:
    import neighbor
    pairs = neighbor.find_first_neighbor(r1,r2)
    return pairs
  else:
    from numpy import array
    n=len(r1)
    pairs = [] # pairs of neighbors
    for i in range(n):
      ri=r1[i]
      for j in range(n):
        rj=r2[j]
        dr = ri - rj ; dr = dr.dot(dr)
        if .9<dr<1.1 : # check if distance is 1
          pairs.append([i,j])  # add to the list
    return pairs # return pairs of first neighbors









def find_first_neighbor_ordered(x1,y1,x2,y2):
  """ Gets the ordered first neighbors, returns a list of tuples """
  from numpy import array
  n=len(x1)
  pairs = [] # pairs of neighbors
  for i in range(n):
    r1=array([x1[i],y1[i]])
    for j in range(n):
      r2=array([x2[j],y2[j]])
      if is_neigh(r1,r2,1.0,0.1): # check if distance is 1
        pairs.append([i,j])  # add to the list
  return pairs # return pairs of first neighbors


def check_second_neighbor(x1,y1,x2,y2,xm,ym):
  """ Check if three positions are second neighbor,
  and returns 0.0,+1 or -1"""
  r1 = np.array([x1,y1]) 
  r2 = np.array([x1,y1]) 
  rm = np.array([xm,ym]) 





# routine to obtain haldane hopping to second neighbors
def get_hal_hop_sec_neigh(x1,y1,x2,y2,xm,ym,mag_field=0.0,value=1.0):
  from scipy.sparse import coo_matrix, bmat
  from numpy import array, matrix, sqrt
  s3=sqrt(3.0)
  n=len(x1)
  mat=np.matrix([[0.0j for i in range(n)] for j in range(n)]) # create matrix
# get pairs of first neighbors ordered
  pairs_1 = find_first_neighbor_ordered(x1,y1,xm,ym) 
# get pairs of first neighbors ordered
  pairs_2 = find_first_neighbor_ordered(xm,ym,x2,y2) 
#  print pairs_1
#  print pairs_2
  for (i1,im1) in pairs_1: # loop over first neighbors of 1
    for (im2,i2) in pairs_2:  # loop over fist neighbors of two
      if im1==im2: # if connected by a common atom
         rm=array([xm[im1],ym[im1]]) # position of the intermediate
         r1=array([x1[i1],y1[i1]])  # position of the first one
         r2=array([x2[i2],y2[i2]])  # position of the second one
         phi_mag = peierls(x1[i1],y1[i1],x2[i2],y2[i2],mag_field) # mag phase
         if is_number(value):  # if value is number assig its value
           fac = value
         if callable(value):  # if value is a function call it
           fac1 = value(x1[i1],y1[i1],x2[i2],y2[i2],xm[i2],ym[i2])
           fac2 = value(x1[i1],y1[i1],x2[i2],y2[i2],xm[i2],ym[i2])
           if not fac1==fac2:
             print "Function in Haldane is not symmetric"
             print fac1,fac2,i1,i2
             raise
           fac = fac1
# add element to the matrix
         mat[i1,i2] += 1j*vec_chi(r1-rm,r2-rm)*phi_mag*fac 
  return mat

# routine to obtain first hoppings
def get_hop_first_neigh(x1,y1,x2,y2,mag_field=0.0):
  n=len(x1)
  mat=np.matrix([[0.0j for i in range(n)] for j in range(n)])
  pairs = find_first_neighbor(x1,y1,x2,y2) # get pairs of first neighbors
  for p in pairs: # loop over pairs
    phi_mag = peierls(x1[p[0]],y1[p[0]],x2[p[1]],y2[p[1]],mag_field) # mag phase
    mat[p[0],p[1]] = phi_mag # element of the matrix
  return mat



def create_fn_hopping(r1,r2):
  n=len(r1)
  mat=np.matrix([[0.0j for i in range(n)] for j in range(n)])
  pairs = get_first_neighbors(r1,r2) # get pairs of first neighbors
  for p in pairs: # loop over pairs
    mat[p[0],p[1]] = 1.0 
  return mat



# function to calculate the chirality between two vectors (vectorial productc)
def vec_chi(r1,r2):
  """Return clockwise or anticlockwise"""
  z=r1[0]*r2[1]-r1[1]*r2[0]
  zz = r1-r2
  zz = sum(zz*zz)
  if zz>0.01:
    if z>0.01: # clockwise
      return 1.0
    if z<-0.01: # anticlockwise
      return -1.0
  return 0.0




# routine to check if two atoms arein adistance d
def is_neigh(r1,r2,d,tol):
  r=r1-r2
  x=r[0]
  y=r[1]
  dt=abs(d*d-x*x-y*y)
  if dt<tol:
    return True
  return False




#################################3

def diagonalize(h,nkpoints=100):
  """ Diagonalice a hamiltonian """
  import scipy.linalg as lg
  # for one dimensional systems
  if h.dimensionality==1:  # one simensional system
    klist = np.arange(0.0,1.0,1.0/nkpoints)  # create array with the klist
    if h.geometry.shift_kspace:
      klist = np.arange(-0.5,0.5,1.0/nkpoints)  # create array with the klist
    intra = h.intra  # assign intraterm
    inter = h.inter  # assign interterm
    energies = [] # list with the energies
    for k in klist: # loop over kpoints
      bf = np.exp(1j*np.pi*2.*k)  # bloch factor for the intercell terms
      inter_k = inter*bf  # bloch hopping
      hk = intra + inter_k + inter_k.H # k dependent hamiltonian
      energies += [lg.eigvalsh(hk)] # get eigenvalues of the current hamiltonian
    energies = np.array(energies).transpose() # each kpoint in a line
    return (klist,energies) # return the klist and the energies
# for zero dimensional systems system
  if h.dimensionality==0:  
    intra = h.intra  # assign intraterm
    energies = lg.eigvalsh(intra) # get eigenvalues of the current hamiltonian
    return (range(len(intra)),energies) # return indexes and energies


def eigenvectors(h,nk=10):
  import scipy.linalg as lg
  from scipy.sparse import csc_matrix as csc
  if h.dimensionality==0:
    vv = lg.eigh(h.intra)
    vecs = [csc(np.matrix(v).T) for v in vv[1].transpose()]
    return vv[0],vecs
  elif h.dimensionality==1:
    kp = np.linspace(0.,1.0,nk)
    eigvecs = [] # empty list with eigenvectors
    eigvals = [] # empty list with eigenvectors
    f = h.get_hk_gen()
    for k in kp: # loop over kpoints
      hk = f(k)  # kdependent hamiltonian
      vv = lg.eigh(hk) # diagonalize k hamiltonian
      eigvals += vv[0].tolist() # store eigenvalues
      vecs = [csc(np.matrix(v).T) for v in vv[1].transpose()]
      eigvecs += vecs # store eigenvectors in sparse form
    return eigvals,eigvecs
  else:
    raise



def plot_bands(h,nkpoints=100,use_lines=False):
  """ Returns a figure with the bandstructure of the system"""
  import pylab
  fig = pylab.figure() # create the figure
  sp = fig.add_subplot(111)
  sp.set_ylabel("Energy")
  sp.set_xlabel("k-vector")
  fig.set_facecolor("white") # white figure  
  if h.dimensionality==1:  # one simensional system
    (klist,energies) = diagonalize(h,nkpoints)
    for band in energies:
      sp.scatter(klist,band,marker = "o",color="black") # print the bands
      sp.set_xlim([min(klist),max(klist)])

  if h.dimensionality==0:  # one simensional system
    (index,energies) = diagonalize(h,nkpoints)
    sp.scatter(index,energies,marker = "o",color="black") # print the bands
  return fig # and return the figure 


def get_bands1d(h,nkpoints=100,operator=None):
  hkgen = h.get_hk_gen() # generator hamiltonian
  ks = np.linspace(0.,1.,nkpoints)
  f = open("BANDS.OUT","w") # open bands file
  if operator!=None: operator=np.matrix(operator) # convert to matrix
  for k in ks: # loop over kpoints
    hk = hkgen(k) # get hamiltonian
    if operator==None:      
      es = lg.eigvalsh(hk)  
      for e in es:  # loop over energies
        f.write(str(k)+"   "+str(e)+"\n") # write in file
    else:
      es,ws = lg.eigh(hk)  
      ws = ws.transpose() # transpose eigenvectors
      for (e,w) in zip(es,ws):  # loop over waves
        w = np.matrix(w) # convert to matrix
        waw = (w.T).H*operator*w.T # expectation value
        waw = waw[0,0].real # real part
        f.write(str(k)+"   "+str(e)+"  "+str(waw)+"\n") # write in file
  f.close()




def get_bands(h,nkpoints=100,use_lines=False,kpath=None):
  """ Returns a figure with the bandstructure of the system"""
  if h.dimensionality==1:  # one simensional system
    (klist,energies) = diagonalize(h,nkpoints)
    energies = energies.transpose()
    fb = open("BANDS.OUT","w")
    for (k,b) in zip(klist,energies): # write in file
      for e in b:
        fb.write(str(k)+"    ")
        fb.write(str(e)+"\n")
    fb.close() # close file 
  elif h.dimensionality==0:  # one simensional system
    (klist,energies) = diagonalize(h,nkpoints)
    fb = open("BANDS.OUT","w")
    for (k,e) in zip(klist,energies): # write in file
        fb.write(str(k)+"    ")
        fb.write(str(e)+"\n")
    fb.close() # close file 


def get_bands_2d(h,kpath=None,operator=None):
  """Get a 2d bandstructure"""
  fo = open("BANDS.OUT","w")
  (klist,energies) = diagonalize_kpath(h,kpath) # get the energies
  for (k,e) in zip(klist,energies): # loop over eigenvalues
    fo.write(str(k)+"      ")
    fo.write(str(e)+"\n")
  fo.close()


def diagonalize_kpath(h,kpath):
  """Diagonalice in a certain path"""
  energies = [] # empty list with energies
  import scipy.linalg as lg
  ik = 0.
  iks = [] # empty list
  for k in kpath:
    f = h.get_hk_gen() # get Hk generator
    hk = f(k)  # k dependent hamiltonian
    es = (lg.eigvalsh(hk)).tolist() # get eigenvalues for current hamiltonian
    energies += es # append energies 
    iks += [ik for i in es]
    ik += 1.
  iks = np.array(iks)
  iks = iks/max(iks) # normalize path
  return (iks,energies)





def print_hamiltonian(h):
  """ Print the hamilotnian on screen """
  from scipy.sparse import coo_matrix as coo # import sparse matrix
  intra = coo(h.intra) # intracell
  inter = coo(h.inter) # intracell
  print "Intracell matrix"
  print intra
  print "Intercell matrix"
  print inter
  return



def add_sublattice_imbalance(h,mass):
  """ Adds to the intracell matrix a sublattice imbalance """
  intra = h.intra # intracell hopping
  if h.geometry.has_sublattice:  # if has sublattice
    def ab(i): 
      return h.geometry.sublattice[i]
  else: # if does not have sublattice
    def ab(i): 
      print "WARNING, system does not have sublattice"
      return (-1)**i

  if h.has_spin:
    natoms = len(intra)/2 # assume spinpolarized calculation 
    if is_number(mass): # if mass is a float
      for i in range(natoms):
        intra[2*i,2*i] += mass*ab(i)
        intra[2*i+1,2*i+1] += mass*ab(i)
    elif callable(mass):  # if mass is a function
      x = h.geometry.x
      y = h.geometry.y
      for i in range(natoms):
        intra[2*i,2*i] += mass(x[i],y[i])*ab(i)
        intra[2*i+1,2*i+1] += mass(x[i],y[i])*ab(i)
    else:  # something else error
      print "Unrecognized lattice imbaance"
      raise
  if not h.has_spin:
    natoms = len(intra) # assume spinpolarized calculation 
    if type(mass)==float: # if mass is a float
      for i in range(natoms):
        intra[i,i] = mass*ab(i)
    elif callable(mass):  # if mass is a function
      x = h.geometry.x
      y = h.geometry.y
      for i in range(natoms):
        intra[i,i] = mass(x[i],y[i])*ab(i)
  h.intra = intra # update onsite







def add_antiferromagnetism(h,mass,d="z"):
  """ Adds to the intracell matrix an antiferromagnetic imbalance """
  intra = h.intra # intracell hopping
  if h.geometry.has_sublattice:  # if has sublattice
    def ab(k): 
      return h.geometry.sublattice[k]
  else: # if does not have sublattice
    def ab(k): 
      print "WARNING, system does not have sublattice"
      return (-1)**k
  if h.has_spin:
    natoms = len(h.geometry.x) # assume spinpolarized calculation 
    if is_number(mass): # if mass is a float
      for i in range(natoms):
        if d=="z":
          intra[2*i,2*i] += mass*ab(i)
          intra[2*i+1,2*i+1] += -mass*ab(i)
        elif d=="x":
          intra[2*i,2*i+1] += mass*ab(i)
          intra[2*i+1,2*i] += mass*ab(i)
        else: raise
    elif callable(mass):  # if mass is a function
      x = h.geometry.x
      y = h.geometry.y
      for i in range(natoms):
        if d=="z":
          intra[2*i,2*i] += mass(x[i],y[i])*ab(i)
          intra[2*i+1,2*i+1] += -mass(x[i],y[i])*ab(i)
        elif d=="x":
          intra[2*i,2*i+1] += mass(x[i],y[i])*ab(i)
          intra[2*i+1,2*i] += mass(x[i],y[i])*ab(i)
        else: raise
    else:  # something else error
      print "unrecognized AF option"
      raise
  else:
    print "no AF for unpolarized hamiltonian"
    raise
















def add_pwave_electron_hole_pairing(h,delta=0.0,mu=0.0,phi=0.0):
  """ Adds an pwave (spin conserving) electron-hole pairing
   potential to first neighbors,
   delta is pairing potential
   mu is chemical potential
   phi is superconducting phase"""
  h.has_eh = True # has electron hole pairs
  intra = h.intra  # intra cell potential
  intra_sc = 0.0*intra  # intracell pairing potential
  x = h.geometry.x # x position
  y = h.geometry.y # y position
  from scipy.sparse import coo_matrix as coo
  from scipy.sparse import bmat
  deltac = delta*np.exp(1j*phi*np.pi) # complex SC order parameter
  # intracell contribution
  pairs = find_first_neighbor(x,y,x,y) # get intra pairs of first neighbors
  for p in pairs: # loop over intracell pairs
 # geometric phase, dependent on bonding orientation
    gp = np.arctan2(y[p[0]]-y[p[1]],x[p[0]]-x[p[1]]) 
    gp = np.exp(1j*gp)  # get complex phase
    if h.has_spin:
      intra_sc[2*p[0],2*p[1]] = deltac*gp    # up up, times phase
      intra_sc[2*p[0]+1,2*p[1]+1] = deltac*gp  # down down, times phase
    else:
      intra_sc[p[0],p[1]] = deltac*gp    # up up, times phase
  # intracell eh matrix
  intra_eh = nambu(intra,coupling=intra_sc)
  if h.dimensionality == 0:   # if 0 dimensional return
    return (intra_eh,None)
  # intercell contribution
  if h.dimensionality == 1: # if one dimensional calculate neighbor contribution
    celldis = h.geometry.celldis # distance to neighboring cell
    inter = h.inter  # inter cell potential
    inter_sc = 0.0*inter  # intercell pairing potential
    pairs = find_first_neighbor(x,y,x+celldis,y) # get inter pairs of first neighbors
    for p in pairs: # loop over intercell pairs
      gp = np.arctan2(x[p[0]]-x[p[1]-celldis],y[p[0]]-y[p[1]]) 
      gp = np.exp(1j*gp)  # get complex phase
      if h.has_spin:
        inter_sc[2*p[0],2*p[1]] = deltac*gp   # up up, times phase
        inter_sc[2*p[0]+1,2*p[1]+1] = deltac*gp  # down down, times phase
      else:
        inter_sc[p[0],p[1]] = deltac*gp   # up up, times phase
    # intercell eh matrix
    inter_eh = nambu(inter,coupling=inter_sc)
#    inter_eh = [[coo(inter),coo(inter_sc)],[-coo(inter_sc.H),-coo(inter)]]
#    inter_eh = bmat(inter_eh).todense()
    return (intra_eh,inter_eh)



def add_swave_electron_hole_pairing(h,delta=0.0,mu=0.0,phi=0.0):
  """ Adds an swave electron-hole pairing
   potential to first neighbors,
   delta is pairing potential
   mu is chemical potential
   phi is superconducting phase,
   delta and phi can be functions!!!!!"""
  h.has_eh = True # has electron hole pairs
  intra = h.intra  # intra cell potential
  intra_sc = 0.0*intra  # intracell pairing potential
  x = h.geometry.x # x position
  y = h.geometry.y # y position
  from scipy.sparse import coo_matrix as coo
  from scipy.sparse import bmat
  # intracell contribution
  pairs = find_first_neighbor(x,y,x,y) # get intra pairs of first neighbors
  for i in range(len(x)): # loop over intracell pairs

    if is_number(delta) and is_number(phi):  # if both are numbers
      deltac = delta*np.exp(1j*phi*np.pi) # complex SC order parameter
    elif callable(delta) and callable(phi): # if both are functions
      deltac = delta(x[i],y[i])*np.exp(1j*phi(x[i],y[i])*np.pi) 
    elif callable(phi) and is_number(delta): # if phi is function
      deltac = delta*np.exp(1j*phi(x[i],y[i])*np.pi) 
    elif callable(delta) and is_number(phi): # if delta is function
      deltac = delta(x[i],y[i])*np.exp(1j*phi*np.pi) 
    else:  # error if unrecognized
      print type(delta),type(phi)  
      raise
    # notation is
        # psi_up
        # psi_dn
        # psi_dn_dag
        # psi_up_dag
    if h.has_spin:
      # couples electron in even orbital with hole in odd orbital
            # in the curren notation even corresponds to up
            # and odd corresponds to down, so this term couples
            # at the same atom, up electrons with down holes
                     # index 0 is first atom up
                     # index 1 is first atom down
      intra_sc[2*i,2*i+1] = deltac    # up down electron hole
      intra_sc[2*i+1,2*i] = -deltac    # down up
  #    intra_sc[2*i,2*i] = deltac    # up down electron hole
  #    intra_sc[2*i+1,2*i+1] = -deltac    # down up
    else:
      raise
  # intracell eh matrix
#  intra_eh = [[coo(intra),coo(intra_sc)],[coo(intra_sc.H)
#             ,-coo(np.conjugate(intra))]]
#             ,-coo(intra)]]
#             ,-coo(intra.H)]]
#  intra_eh = bmat(intra_eh).todense()
  intra_eh = build_eh(intra,coupling=intra_sc)
  # intercell contribution
  if h.dimensionality == 1: # if one dimensional calculate neighbor contribution
    inter = h.inter  # inter cell potential
    # intercell eh matrix, no mixing terms
#    inter_eh = [[coo(inter),None],[None,-coo(np.conjugate(inter))]]
#    inter_eh = bmat(inter_eh).todense()
    inter_eh = build_eh(inter)  # create electron hole matrix
    return (intra_eh,inter_eh)
  else:   # if 0 dimensional return
    return (intra_eh,None)





def build_eh(hin,coupling=None):
  """Creates a electron hole matrix, from an input matrix, coupling couples
     electrons and holes
      - hin is the hamiltonian for electrons, which has the usual common form
      - coupling is the matrix which tells the coupling between electron
        on state i woth holes on state j, for exmaple, with swave pairing
        the non vanishing elments are (0,1),(2,3),(4,5) and so on..."""
  n = len(hin)  # dimension of input
  nn = 2*n  # dimension of output
  hout = np.matrix(np.zeros((nn,nn),dtype=complex))  # output hamiltonian
  for i in range(n):
    for j in range(n):
      hout[2*i,2*j] = hin[i,j]  # electron term
      hout[2*i+1,2*j+1] = -np.conjugate(hin[i,j])  # hole term
  if not coupling == None: # if there is coupling
    for i in range(n):
      for j in range(n):
        # couples electron in i with hole in j
        hout[2*i,2*j+1] = coupling[i,j]  # electron hole term
        # couples hole in i with electron in j
        hout[2*j+1,2*i] = np.conjugate(coupling[i,j])  # hole electron term
  return hout 











def set_finite_system(h):
  """ Transforms the hamiltonian into a finite system,
  removing the hoppings """
  from copy import deepcopy
  if h.dimensionality == 1:
    h.intra = h.intra + h.inter + h.inter.H 
  if h.dimensionality == 2:
    h.intra = h.intra +  h.tx + h.tx.H 
    h.intra = h.intra +  h.ty + h.ty.H 
    h.intra = h.intra +  h.txy + h.txy.H 
    h.intra = h.intra +  h.txmy + h.txmy.H 
  h.dimensionality = 0 # put dimensionality = 0
  

def peierls(x1,y1,x2,y2,mag_field):
  """ Returns the complex phase with magnetic field """
  if is_number(mag_field): # if it is a number assumen Landau gauge
    phase = mag_field*(x1-x2)*(y1+y2)/2.0  
  elif callable(mag_field): # if it is callable 
    phase = mag_field(x1,y1,x2,y2)  # specific call to the value
  else:  # if anything else error
    raise
  return np.exp(1j*phase)


def total_energy(h,nkpoints=100):
  """ Returns total energy of the system (assuming Efermi = 0)"""
  if h.dimensionality == 0: # 0 dimensional systems
    (k,eigen) = diagonalize(h,nkpoints=nkpoints)
  etot = 0.0
  for e in eigen:
    if e<0.0: # return if reached positive energy
      etot += e
  return etot # return if reached positive energy

def add_zeeman(h,zeeman=[0.0,0.0,0.0]):
  """ Add Zeeman to the hamiltonian """
  # convert the input into a list
  def convert(z):
    try: # try to get the first element
      a = z[0]
    except: # convert to list
      z = [0.,0.,z]
    return z
  from scipy.sparse import coo_matrix as coo
  from scipy.sparse import bmat
  if h.has_spin: # only if the system has spin
    sx = coo([[0.0,1.0],[1.0,0.0]])   # sigma x
    sy = coo([[0.0,-1j],[1j,0.0]])   # sigma y
    sz = coo([[1.0,0.0],[0.0,-1.0]])   # sigma z
   # no = h.num_orbitals # number of orbitals (without spin)
    no = len(h.geometry.x) # number of orbitals (without spin)
    # create matrix to add to the hamiltonian
    bzee = [[None for i in range(no)] for j in range(no)]
    # assign diagonal terms
    if not callable(zeeman): # if it is a number
      for i in range(no):
        zeeman = convert(zeeman) # convert to list
        bzee[i][i] = zeeman[0]*sx+zeeman[1]*sy+zeeman[2]*sz
    elif callable(zeeman): # if it is a function
      x = h.geometry.x  # x position
      y = h.geometry.y  # y position
      for i in range(no):
        z = zeeman(x[i],y[i])  # get the value of the zeeman
        z = convert(z) # convert to list
        bzee[i][i] = z[0]*sx+z[1]*sy+z[2]*sz
    else:
      raise
    bzee = bmat(bzee) # create matrix
    if h.has_eh: # if electron hole, double the dimension
      bzee = bmat([[bzee,None],[None,-bzee]])
    bzee = bzee.todense() # create dense matrix
    h.intra = h.intra + bzee
  if not h.has_spin:  # still have to implement this...
    raise


def des_spin(m,component=0):
  """ Removes the spin degree of freedom"""
  d = len(m) # dimension of the matrix
  if d%2==1: # if the hamiltonian doesn't have the correct dimension
    print "Hamiltonian dimension is odd"
    raise
  mout = np.matrix([[0.0j for i in range(d/2)] for j in range(d/2)])
  for i in range(d/2):
    for j in range(d/2):
      mout[i,j] = m[2*i+component,2*j+component]  # assign spin up part
  return mout


def shift_fermi(h,fermi):
  """ Moves the fermi energy of the system, the new value is at zero"""
  intra = h.intra # get intraterm
  n = len(intra) # number of elements in the matrix
  if is_number(fermi):  # if fermi is a number
#    print "Shifting FERMI energy, input is NUMBER"
    if h.has_eh:  # change signs for electrons and holes
      for i in range(n/2):      
        intra[2*i,2*i] += fermi
        intra[2*i+1,2*i+1] += -fermi
    else:  # same sign for everyone
      for i in range(n):
        h.intra[i,i] += fermi
  if callable(fermi):  # if fermi is a function
    print "Shifting FERMI energy, input is FUNCTION"
    x = h.geometry.x
    y = h.geometry.y
    z = h.geometry.z
    if h.has_spin: indd = 2  # with spin
    else: indd = 1   # without spin
    if h.has_eh:  # change signs for electrons and holes
      for i in range(n/2):      
        intra[2*i,2*i] += fermi(x[i/indd],y[i/indd],z[i/indd])
        intra[2*i+1,2*i+1] += -fermi(x[i/indd],y[i/indd],z[i/indd])
    else:  # same sign for everyone
      for i in range(n):
        intra[i,i] += fermi(x[i/indd],y[i/indd],z[i/indd])

  h.intra = intra # assign new matrix



def is_number(s):
    try:
        float(s)
        return True
    except:
        return False

def is_hermitic(m):
  mh = m.H
  hh = m - m.H
  for i in range(len(hh)):
    for j in range(len(hh)):
      if np.abs(hh[i,j]) > 0.000001:
        print "No hermitic element", i,j,m[i,j],m[j,i]
        return False
  return True
  






#########################
# first neighbor hopping
##########################
def first_neighbors2d_old(h):
  """ Gets a first neighbor hamiltonian"""
  x = h.geometry.x    # x coordinate 
  y = h.geometry.y    # x coordinate 
  g = h.geometry
  celldis = h.geometry.celldis    # distance to neighboring cell
# first neighbors hopping
  h.intra = get_hop_first_neigh(x,y,x,y,mag_field=0.0)
#  if h.dimensionality==2: # if one dimensional
  if True: # if one dimensional
    dx = g.a1[0] # vector to the cell
    dy = g.a1[1] # vector to the cell
    h.tx = get_hop_first_neigh(x,y,x+dx,y+dy,mag_field=0.0)
    dx = g.a2[0] # vector to the cell
    dy = g.a2[1] # vector to the cell
    h.ty = get_hop_first_neigh(x,y,x+dx,y+dy,mag_field=0.0)
    dx = g.a1[0] + g.a2[0] # vector to the cell
    dy = g.a1[1] + g.a2[1] # vector to the cell
    h.txy = get_hop_first_neigh(x,y,x+dx,y+dy,mag_field=0.0)
    dx = g.a1[0] - g.a2[0] # vector to the cell
    dy = g.a1[1] - g.a2[1] # vector to the cell
    h.txmy = get_hop_first_neigh(x,y,x+dx,y+dy,mag_field=0.0)
  from increase_hilbert import m2spin
  if h.has_spin: # if it has spin degree of freedom
    h.intra = m2spin(h.intra,h.intra) # intracell term 
    h.tx = m2spin(h.tx,h.tx) # hopping
    h.ty = m2spin(h.ty,h.ty) # hopping
    h.txy = m2spin(h.txy,h.txy) # hopping
    h.txmy = m2spin(h.txmy,h.txmy) # hopping



def first_neighbors2d(h):
  """ Gets a first neighbor hamiltonian"""
  r = h.geometry.r    # x coordinate 
  g = h.geometry
# first neighbors hopping, all the matrices
  a1, a2 = g.a1, g.a2
  h.intra = create_fn_hopping(r,r)
  h.tx = create_fn_hopping(r,r+a1)
  h.ty = create_fn_hopping(r,r+a2)
  h.txy = create_fn_hopping(r,r+a1+a2)
  h.txmy = create_fn_hopping(r,r+a1-a2)
  from increase_hilbert import m2spin
  if h.has_spin: # if it has spin degree of freedom
    h.intra = m2spin(h.intra,h.intra) # intracell term 
    h.tx = m2spin(h.tx,h.tx) # hopping
    h.ty = m2spin(h.ty,h.ty) # hopping
    h.txy = m2spin(h.txy,h.txy) # hopping
    h.txmy = m2spin(h.txmy,h.txmy) # hopping



def first_neighbors1d(h):
  """ Gets a first neighbor hamiltonian"""
  r = h.geometry.r    # x coordinate 
  a1 = np.array([h.geometry.celldis,0.,0.])
# first neighbors hopping, all the matrices
  h.intra = create_fn_hopping(r,r)
  h.inter = create_fn_hopping(r,r+a1)
  if h.has_spin: # if it has spin degree of freedom
    from increase_hilbert import m2spin
    h.intra = m2spin(h.intra,h.intra) # intracell term 
    h.inter = m2spin(h.inter,h.inter) # hopping


def first_neighbors0d(h):
  """ Gets a first neighbor hamiltonian"""
  r = h.geometry.r    # x coordinate 
  if h.is_sparse: # for sparse matrix
    import neighbor
    pairs = neighbor.find_first_neighbor(r,r)
    rows,cols = pairs.T # transpose
    data = np.array([1. for c in cols])
    n = len(r)
    h.intra = csc_matrix((data,(rows,cols)),shape=(n,n))
  else:
    h.intra = create_fn_hopping(r,r)  # intracell
    if h.has_spin: # if it has spin degree of freedom
      from increase_hilbert import m2spin
      h.intra = m2spin(h.intra,h.intra) # intracell term 




def nambu(m,coupling=None):
  """ Creates a nambu matrix"""
  mo = [[None,None],[None,None]] # create output matrix
  mo[0][0] = csc_matrix(m) # electron 
  mo[1][1] = -np.conjugate(csc_matrix(m)) # hole
  if not coupling == None:
    mo[0][1] = csc_matrix(coupling) # electron-hole 
    mo[1][0] = csc_matrix(coupling.H) # hole-electron 
  return bmat(mo).todense() # return matrix

nambu = build_eh


def add_swave(intra,delta=0.0,mu=0.0,phi=0.0):
  """ Adds swave pairing """
  intra_sc = 0.0*intra  # intracell pairing potential
  from scipy.sparse import coo_matrix as coo
  from scipy.sparse import bmat
  # intracell contribution
  deltac = delta*np.exp(1j*phi*np.pi) # complex SC order parameter
  for i in range(len(intra)/2): # loop over intracell pairs
    intra_sc[2*i,2*i+1] = deltac    # up down electron hole
    intra_sc[2*i+1,2*i] = -deltac    # down up
  intra_eh = nambu(intra,coupling=intra_sc)
  return intra_eh





def lowest_bands(h,nkpoints=100,nbands=10,operator = None,info = False):
  """ Returns a figure with the bandstructure of the system"""
  from scipy.sparse import csc_matrix
  intra = csc_matrix(h.intra)
  t = csc_matrix(h.inter)
  k = np.arange(0.,1.,1./nkpoints) # list of kpoints
  import scipy.sparse.linalg as lg
  fo = open("BANDS.OUT","w")
  if operator==None: # if there is not an operator
    if h.dimensionality==0:  # dot
      eig,eigvec = lg.eigsh(csc_matrix(h.intra),k=nbands,which="LM",sigma=0.0)
      for i in range(len(eig)):
        fo.write(str(i)+"     "+str(eig[i])+"\n")


    if h.dimensionality==1:
      for ik in k:  # ribbon
        tk = t*np.exp(1j*2.*np.pi*ik)
        hk = intra +tk + tk.H
        eig,eigvec = lg.eigsh(hk,k=nbands,which="LM",sigma=0.0)
        for e in eig:
          fo.write(str(ik)+"     "+str(e)+"\n")
  else:  # if there is an operator
    if h.dimensionality==1:
      for ik in k:
        tk = t*np.exp(1j*2.*np.pi*ik)
        hk = intra +tk + tk.H
        eig,eigvec = lg.eigsh(hk,k=nbands,which="LM",sigma=0.0)
        eigvec = eigvec.transpose() # tranpose the matrix
        if info:  print "Done",ik
        for (e,v) in zip(eig,eigvec): # loop over eigenvectors
          v = csc_matrix(v)
        #  print v.shape,type(v),operator.shape
          a = np.conjugate(v) * operator * v.T
          a = a.data[0]
          a = a.real # real part
        #  print a 
          fo.write(str(ik)+"     "+str(e)+"     "+str(a)+"\n")
  fo.close()


def add_peierls(h,mag_field=0.0):
  """ Adds Peierls phase to the Hamiltonian"""
  x = h.geometry.x    # x coordinate 
  y = h.geometry.y    # x coordinate 
  celldis = h.geometry.celldis    # distance to neighboring cell
  from numpy import array
  norb = h.intra.shape[0]  # number of orbitals
  if h.is_sparse: # sparse hamiltonian
    from scipy.sparse import coo_matrix,csc_matrix
    if h.has_spin: raise # raise if spinful
    if h.dimensionality==0: # zero dimensional
      m = coo_matrix(h.intra) # convert to sparse matrix
      row,col = m.row,m.col
      data = m.data +0j
      for k in range(len(m.data)): # loop over non vanishing elements
        i = m.row[k]
        j = m.col[k]
        p = peierls(x[i],y[i],x[j],y[j],mag_field) # peierls phase
        data[k] *= p # add phase
      h.intra = csc_matrix((data,(row,col)),shape=(norb,norb)) # convert to csc
    if h.dimensionality==1: # one dimensional
      m = coo_matrix(h.inter) # convert to sparse matrix
      row,col = m.row,m.col
      data = m.data +0j
      for k in range(len(m.data)): # loop over non vanishing elements
        i = m.row[k]
        j = m.col[k]
        p = peierls(x[i],y[i],x[j]+celldis,y[j],mag_field) # peierls phase
        data[k] *= p # add phase
      h.inter = csc_matrix((data,(row,col)),shape=(norb,norb)) # convert to csc
    if h.dimensionality>1: raise # error if greater than 1
  else: # not sparse
    def gaugeize(m,d=0.0):
      """Add gauge phase to a matrix"""
      for i in range(len(x)):
        for j in range(len(x)):
          p = peierls(x[i],y[i],x[j]+d,y[j],mag_field) # peierls phase
          if h.has_spin:
            m[2*i,2*j] *= p
            m[2*i,2*j+1] *= p
            m[2*i+1,2*j] *= p
            m[2*i+1,2*j+1] *= p
          else:
            m[i,j] *= p
    gaugeize(h.intra,d=0.0)  # gaugeize intraterm
    if h.dimensionality==0: pass # if zero dimensional
    elif h.dimensionality==1: # if one dimensional
      gaugeize(h.inter,d=celldis) # gaugeize interterm
    elif h.dimensionality==2: # if bigger dimensional
      print "WARNING, is your gauge periodic?"
      gaugeize(h.tx,d=h.geometry.a1[0]) # gaugeize interterm
    else:
      raise


def hk_gen(h):
  """ Returns a function that generates a k dependent hamiltonian"""
  if h.dimensionality == 0: return None
  if h.dimensionality == 1: 
    def hk(k):
      """k dependent hamiltonian, k goes from 0 to 1"""
      tk = h.inter * np.exp(1j*np.pi*2.*k)
      ho = h.intra + tk + tk.H
      return ho
    return hk  # return the function
  if h.dimensionality == 2: 
    def hk(k):
      """k dependent hamiltonian, k goes from (0,0) to (1,1)"""
      if len(k)==3:
        k = np.array([k[0],k[1]]) # redefine for 2d
      k = np.array(k)
      ux = np.array([1.,0.])
      uy = np.array([0.,1.])
      ptk = [[h.tx,ux],[h.ty,uy],[h.txy,ux+uy],[h.txmy,ux-uy]] 
      ho = (h.intra).copy() # intraterm
      for p in ptk: # loop over hoppings
        tk = p[0]*np.exp(1j*np.pi*2.*(p[1].dot(k)))  # add bloch hopping
        ho += tk + tk.H  # add bloch hopping
      return ho
    return hk




def ldos(h,e=0.0,delta=0.01):
  """Calculates the local density of states of a hamiltonian and
     writes it in file"""
  if h.dimensionality==0:
    iden = np.identity(h.intra.shape[0],dtype=np.complex) # create identity
    g = ( (e+1j*delta)*iden -h.intra ).I # calculate green function
  else: raise # not implemented...
  if h.has_spin: # calculate LDOS for spinpol
    d = [ -(g[2*i,2*i] + g[2*i+1,2*i+1]).imag for i in range(len(g)/2)]
  else: # calculate LDOS for spin unpol
    d = [ -(g[i,i]).imag for i in range(len(g))]
  fd = open("LDOS.OUT","w")   # open file
  g = h.geometry  # store geometry
  for (ix,iy,idos) in zip(g.x,g.y,d): # write everything
    fd.write(str(ix) +"   "+ str(iy) + "   "+ str(idos)+"\n")
  fd.close() # close file





def generate_parametric_hopping(h,f):
  """ Adds a parametric hopping to the hamiltonian based on an input function"""
  rs = h.geometry.r # positions
  g = h.geometry # geometry
  h.has_spin = False
  if h.is_sparse:
    data = []
    rows = []
    cols = []
    if h.dimensionality == 0:
      for i in range(len(rs)):
        for j in range(len(rs)):
          c = f(rs[i],rs[j]) # get coupling
          if np.abs(c)>0.001: #cutoff
             rows.append(i)
             cols.append(j)
             data.append(c)
      n = len(rs) # dimension of the matrix
      h.intra = csc_matrix((data,(rows,cols)),shape=(n,n)) # store in hamil
    else: raise # error if not 0d
    return h
  else: # not sparse 
    h.intra = parametric_hopping(rs,rs,f)
    if h.dimensionality == 0: pass
    elif h.dimensionality == 1:
      dr = np.array([g.celldis,0.,0.])
      h.inter = parametric_hopping(rs,rs+dr,f)
    else: raise
    return h


def parametric_hopping(r1,r2,f):
  """ Generates a parametric hopping based on a function"""
  m = np.matrix([[0.0j for i in range(len(r1))] for j in range(len(r2))])
  for i in range(len(r1)):
    for j in range(len(r2)):
      m[i,j] = f(r1[i],r2[j]) # add hopping based on function
  return m





