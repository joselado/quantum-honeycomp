from __future__ import print_function
from __future__ import division
import numpy as np
from scipy.sparse import csc_matrix,bmat,csr_matrix
import scipy.linalg as lg
import scipy.sparse.linalg as slg
import numpy as np
from . import operators
from . import inout
from . import superconductivity
from . import kanemele 
from . import magnetism
from . import checkclass
from . import extract
from . import multicell
from . import spectrum
from . import kekule
from . import algebra
from . import groundstate
from . import rotate_spin
from .bandstructure import get_bands_nd

from scipy.sparse import coo_matrix,bmat
from .rotate_spin import sx,sy,sz
from .increase_hilbert import get_spinless2full,get_spinful2full
from . import tails
from scipy.sparse import diags as sparse_diag
import pickle

#import data

from .limits import densedimension as maxmatrix
#maxmatrix = 4000 # maximum dimension
optimal = False

class hamiltonian():
  """ Class for a hamiltonian """
  def get_tails(self,discard=None):
    """Write the tails of the wavefunctions"""
    if self.dimensionality!=0: raise
    else: return tails.matrix_tails(self.intra,discard=discard)
  def spinless2full(self,m,time_reversal=False):
    """Transform a spinless matrix in its full form"""
    return get_spinless2full(self,time_reversal=time_reversal)(m) # return
  def spinful2full(self,m):
    """Transform a spinless matrix in its full form"""
    return get_spinful2full(self)(m) # return
  def kchain(self,k=0.):
    return kchain(self,k)
  def eigenvectors(self,nk=10,kpoints=False,k=None,sparse=False,numw=None):
    return eigenvectors(self,nk=nk,kpoints=kpoints,k=k,
                                 sparse=sparse,numw=numw)
  def get_filling(self,energy=0.5,nk=10):
    """Get the filling of a Hamiltonian at this energy"""
    es = spectrum.eigenvalues(self,nk=nk) # eigenvalues
    es = np.array(es)
    esf = es[es<energy]
    return len(esf)/len(es) # return filling
  def set_filling(self,filling=0.5,nk=10,extrae=0.):
    """Set the filling of the Hamiltonian"""
    es = spectrum.eigenvalues(self,nk=nk)
    from .scftypes import get_fermi_energy
    fill = filling + extrae/self.intra.shape[0] # filling
    self.shift_fermi(-get_fermi_energy(es,fill)) # shift the fermi energy
  def __init__(self,geometry=None):
    self.data = dict() # empty dictionary with various data
    self.has_spin = True # has spin degree of freedom
    self.prefix = "" # a string used a prefix for different files
    self.path = "" # a path used for different files
    self.has_eh = False # has electron hole pairs
    self.get_eh_sector = None # no function for getting electrons
    self.fermi_energy = 0.0 # fermi energy at zero
    self.dimensionality = 0 # dimensionality of the Hamiltonian
    self.is_sparse = False
    self.is_multicell = False # for hamiltonians with hoppings to several neighbors
    self.hopping_dict = {} # hopping dictonary
    self.non_hermitian = False # non hermitian Hamiltonian
    if not geometry is None:
# dimensionality of the system
      self.dimensionality = geometry.dimensionality 
      self.geometry = geometry # add geometry object
      self.num_orbitals = len(geometry.x)
  def get_hk_gen(self):
    """ Generate kdependent hamiltonian"""
    if self.is_multicell: return multicell.hk_gen(self) # for multicell
    else: return hk_gen(self) # for normal cells
  def get_ldos(self,nk=4,e=0.0,mode="arpack",delta=0.05,nrep=3):
      from . import ldos
      ldos.ldos(self,e=e,delta=delta,nk=nk,mode=mode,nrep=nrep)
  def get_gk_gen(self,delta=0.05,operator=None,canonical_phase=False):
    """Return the Green function generator"""
    hkgen = self.get_hk_gen() # Hamiltonian generator
    def f(k=[0.,0.,0.],e=0.0):
      hk = hkgen(k) # get matrix
      if canonical_phase: # use a Bloch phase in all the sites
        frac_r = self.geometry.frac_r # fractional coordinates
        # start in zero
        U = np.diag([self.geometry.bloch_phase(k,r) for r in frac_r])
        U = np.matrix(U) # this is without .H
        U = self.spinless2full(U) # increase the space if necessary
        hk = U.H@hk@U
#        print(csc_matrix(np.angle(hk)))
#        exit()
      if operator is not None: hk = operator.H@hk@operator # project
      out = (np.identity(hk.shape[0])*(e+1j*delta) - hk).I 
#      print(self.geometry.frac_r) 
#      exit()
      return out
    return f
  def print_hamiltonian(self):
    """Print hamiltonian on screen"""
    print_hamiltonian(self)
  def diagonalize(self,nkpoints=100):
    """Return eigenvalues"""
    return diagonalize(self,nkpoints=nkpoints)
  def get_dos(self,**kwargs):
      from . import dos
      return dos.dos(self,**kwargs)
  def get_bands(self,**kwargs):
    """ Returns a figure with teh bandstructure"""
    return get_bands_nd(self,**kwargs)
  def plot_bands(self,nkpoints=100,use_lines=False):
    """Dummy function"""
    self.get_bands()
  def add_sublattice_imbalance(self,mass):
    """ Adds a sublattice imbalance """
    if self.geometry.has_sublattice and self.geometry.sublattice_number==2:
      add_sublattice_imbalance(self,mass)
  def add_antiferromagnetism(self,mass):
    """ Adds antiferromagnetic imbalanc """
    if self.geometry.has_sublattice:
      if self.geometry.sublattice_number==2:
        magnetism.add_antiferromagnetism(self,mass)
      elif self.geometry.sublattice_number>2:
        magnetism.add_frustrated_antiferromagnetism(self,mass)
      else: raise
  def turn_nambu(self):
    """Add electron hole degree of freedom"""
    self.get_eh_sector = get_eh_sector_odd_even # assign function
    turn_nambu(self)
  def add_swave(self,delta=0.0,phi=None):
    """ Adds spin mixing insite electron hole pairing"""
    self.turn_nambu() # add electron hole
    if phi is not None: delta = delta*np.exp(1j*phi*np.pi)
    self.intra = self.intra + add_swave(delta=delta,rs=self.geometry.r,is_sparse=self.is_sparse)
  def add_pairing(self,delta=0.0,**kwargs):
    """ Add a general pairing matrix, uu,dd,ud"""
    superconductivity.add_pairing_to_hamiltonian(self,delta=delta,**kwargs)
  def same_hamiltonian(self,h,ntries=10):
      """Check if two hamiltonians are the same"""
      hk1 = self.get_hk_gen()
      hk2 = h.get_hk_gen()
      for i in range(ntries):
        k = np.random.random(3)
        m = hk1(k) - hk2(k)
        if np.max(np.abs(m))>0.000001: return False
      return True
      
  def supercell(self,nsuper):
    """ Creates a supercell of a one dimensional system"""
    if nsuper==1: return self
    if self.dimensionality==0: return self
    elif self.dimensionality==1: ns = [nsuper,1,1]
    elif self.dimensionality==2: ns = [nsuper,nsuper,1]
    elif self.dimensionality==3: ns = [nsuper,nsuper,nsuper]
    else: raise
    return multicell.supercell_hamiltonian(self,nsuper=ns)
  def set_finite_system(self,periodic=True):
    """ Transforms the system into a finite system"""
    return set_finite_system(self,periodic=periodic) 
  def get_gap(self):
    """Returns the gap of the Hamiltonian"""
    from . import gap
    return gap.indirect_gap(self) # return the gap
  def save(self,output_file="hamiltonian.pkl"):
    """ Write the hamiltonian in a file"""
    inout.save(self,output_file) # write in a file
  write = save # just in case
  def read(self,output_file="hamiltonian.pkl"):
    """ Read the Hamiltonian"""
    return load(output_file) # read Hamiltonian
  def load(self,**kwargs): self.read(**kwargs)
  def get_total_energy(self,**kwargs):
    """ Get total energy of the system"""
    from .spectrum import total_energy
    return total_energy(self,**kwargs)
  def total_energy(self,**kwargs): return self.get_total_energy(**kwargs)
  def add_zeeman(self,zeeman):
    """Adds zeeman to the matrix """
    if self.has_spin:  # if it has spin degree of freedom
      from .magnetism import add_zeeman
      add_zeeman(self,zeeman=zeeman)
  def add_magnetism(self,m):
    """Adds magnetism, new version of zeeman"""
    if self.has_spin:  # if it has spin degree of freedom
      from .magnetism import add_magnetism
      add_magnetism(self,m)
  def turn_spinful(self,enforce_tr=False):
    """Turn the hamiltonian spinful""" 
    if self.is_sparse: # sparse Hamiltonian
      self.turn_dense() # dense Hamiltonian
      self.turn_spinful(enforce_tr=enforce_tr) # spinful
      self.turn_sparse()
    else: # dense Hamiltonian
      if self.has_spin: return # already spinful
      if self.is_multicell: # if multicell
        from .multicell import turn_spinful as ts
        ts(self) # turn spinful
      else:
        from .increase_hilbert import spinful
        from .superconductivity import time_reversal
        def fun(m):
            if enforce_tr: return spinful(m,np.conjugate(m))
            else: return spinful(m)
        if not self.has_spin:
          self.has_spin = True
          self.intra = fun(self.intra) 
          if self.dimensionality==0: pass
          elif self.dimensionality==1:
            self.inter = fun(self.inter) 
          elif self.dimensionality==2:
            self.tx = fun(self.tx) 
            self.ty = fun(self.ty) 
            self.txy = fun(self.txy) 
            self.txmy = fun(self.txmy) 
          else: raise
  def remove_spin(self):
    """Removes spin degree of freedom"""
    if self.has_spin: # if has spin remove the component
      self.intra = des_spin(self.intra,component=0)
      if self.is_multicell: # conventional
        for i in range(len(self.hopping)):
          self.hopping[i].m = des_spin(self.hopping[i].m)
      else:
        if self.dimensionality==1: # if one dimensional
          self.inter = des_spin(self.inter,component=0)
        elif self.dimensionality==2: # if one dimensional
          self.tx = des_spin(self.tx,component=0)
          self.txy = des_spin(self.txy,component=0)
          self.ty = des_spin(self.ty,component=0)
          self.txmy = des_spin(self.txmy,component=0)
        else: raise
      self.has_spin = False  # flag for nonspin calculation
  def add_onsite(self,fermi):
    """ Move the Fermi energy of the system"""
    shift_fermi(self,fermi)
  def shift_fermi(self,fermi): self.add_onsite(fermi)  
  def first_neighbors(self):
    """ Create first neighbor hopping"""
    if 0<=self.dimensionality<3:
      first_neighborsnd(self)
    elif self.dimensionality == 3:
      from .multicell import first_neighbors as fnm
      fnm(self)
    else: raise
  def add_hopping_matrix(self,fm):
      """
      Add a certain hopping matrix to the Hamiltonian
      """
      if not self.is_multicell: raise # this may not work for multicell
      h = self.geometry.get_hamiltonian(has_spin=self.has_spin,
              is_multicell=self.is_multicell,
              mgenerator=fm) # generate a new Hamiltonian
      self.add_hamiltonian(h) # add this contribution
  def add_hamiltonian(self,h):
      """
      Add the hoppings of another Hamiltonian
      """
      if not self.is_multicell: raise # not implemented
      hd = h.get_dict() # get the dictionary
      self.intra = self.intra + hd[(0,0,0)] # add the matrix
      for i in range(len(self.hopping)):
          d = tuple(self.hopping[i].dir)
          if d in hd:
            self.hopping[i].m = self.hopping[i].m + hd[d]
  def get_dict(self):
      """
      Return the dictionary that yields the hoppings
      """
      if not self.is_multicell: raise # not implemented
      hop = dict()
      hop[(0,0,0)] = self.intra
      for t in self.hopping: hop[tuple(t.dir)] = t.m
      return hop # return dictionary
  def copy(self):
    """
    Return a copy of the hamiltonian
    """
    from copy import deepcopy
    return deepcopy(self)
  def check(self):
    """
    Check if the Hamiltonian is OK
    """
    from . import check
    check.check_hamiltonian(self) # check the Hamiltonian
  def turn_sparse(self):
    """
    Transforms the hamiltonian into a sparse hamiltonian
    """
    from scipy.sparse import csc_matrix
#    if self.is_sparse: return # if it is sparse return
    self.is_sparse = True # sparse flag to true
    self.intra = csc_matrix(self.intra)
    if self.is_multicell: # multicell Hamiltonian 
      for i in range(len(self.hopping)): # loop
        self.hopping[i].m = csc_matrix(self.hopping[i].m)
      return
    else:  # no multicell
      if self.dimensionality == 1: # for one dimensional
        self.inter = csc_matrix(self.inter)
      if self.dimensionality == 2: # for one dimensional
        self.tx = csc_matrix(self.tx)
        self.ty = csc_matrix(self.ty)
        self.txy = csc_matrix(self.txy)
        self.txmy = csc_matrix(self.txmy)
      if self.dimensionality > 2: raise # for one dimensional
  def turn_dense(self):
    """ Transforms the hamiltonian into a sparse hamiltonian"""
#    if not self.is_sparse: return
    from scipy.sparse import csc_matrix
#    if not self.is_sparse: return # if it is sparse return
    if self.intra.shape[0]>maxmatrix: raise
    self.is_sparse = False # sparse flag to false
    from scipy.sparse import issparse
    def densify(m):
        if issparse(m): return m.todense()
        else: return m

    self.intra = densify(self.intra)
    if self.is_multicell: # multicell Hamiltonian 
      for i in range(len(self.hopping)): # loop
        self.hopping[i].m = densify(self.hopping[i].m)
      return
    else:  # no multicell
      if self.dimensionality == 0: pass # for one dimensional
      elif self.dimensionality == 1: # for one dimensional
        self.inter = densify(self.inter)
      elif self.dimensionality == 2: # for one dimensional
        self.tx = densify(self.tx)
        self.ty = densify(self.ty)
        self.txy = densify(self.txy)
        self.txmy = densify(self.txmy)
      else: raise

  def add_rashba(self,c):
    """Adds Rashba coupling"""
    from . import rashba
    rashba.add_rashba(self,c)
  def add_kane_mele(self,t):
    """ Adds a Kane-Mele SOC term"""  
    kanemele.add_kane_mele(self,t) # return kane-mele SOC
  def add_haldane(self,t):
    """ Adds a Haldane term"""  
    kanemele.add_haldane(self,t) # return Haldane SOC
  def add_kekule(self,t):
      """
      Add Kekule coupling
      """
      if self.dimensionality==0: # zero dimensional
        m = kekule.kekule_matrix(self.geometry.r,t=t)
        self.intra = self.intra + self.spinless2full(m)
      else: # workaround for higher dimensionality
        r = self.geometry.multireplicas(2) # get many replicas
        fm = kekule.kekule_function(r,t=t)
        self.add_hopping_matrix(fm) # add the Kekule hopping
  def add_chiral_kekule(self,**kwargs):
      """
      Add a chiral kekule hopping
      """
      fun = kekule.chiral_kekule(self.geometry,**kwargs)
      self.add_kekule(fun)

  def add_modified_haldane(self,t):
    """
    Adds a Haldane term
    """  
    kanemele.add_modified_haldane(self,t) # return Haldane SOC
  def add_anti_kane_mele(self,t):
    """
    Adds an anti Kane-Mele term
    """  
    kanemele.add_anti_kane_mele(self,t) # return anti kane mele SOC
  def add_antihaldane(self,t): self.add_modified_haldane(t) # second name
  def add_crystal_field(self,v):
      from . import crystalfield
      crystalfield.hartree(self,v=v) # store
  def add_peierls(self,mag_field,new=False):
    """
    Add magnetic field
    """
    from .peierls import add_peierls
    add_peierls(self,mag_field=mag_field,new=new)
  def add_inplane_bfield(self,**kwargs):
      """Add in-plane magnetic field"""
      from .peierls import add_inplane_bfield
      add_inplane_bfield(self,**kwargs)
  def align_magnetism(self,vectors):
    """ Rotate the Hamiltonian to have magnetism in the z direction"""
    if self.has_eh: raise
    from .rotate_spin import align_magnetism as align
    self.intra = align(self.intra,vectors)
    self.inter = align(self.inter,vectors)
  def global_spin_rotation(self,**kwargs):
      """ Perform a global spin rotation """
      return rotate_spin.hamiltonian_spin_rotation(self,**kwargs)
  def generate_spin_spiral(self,**kwargs):
      """ Generate a spin spiral antsaz in the Hamiltonian """
      return rotate_spin.generate_spin_spiral(self,**kwargs)
  def get_magnetization(self,nkp=10):
      """ Return the magnetization """
      mx = self.extract(name="mx")
      my = self.extract(name="my")
      mz = self.extract(name="mz")
      return np.array([mx,my,mz]).T # return array
  def compute_vev(self,name="sz",**kwargs):
      """
      Compute a VEV of a spatially resolved operator
      """
      n = len(self.geometry.r) # number of sites
      ops = [operators.index(self,n=[i]) for i in range(n)]
      if name=="sx": op = operators.get_sx(self)
      elif name=="sy": op = operators.get_sy(self)
      elif name=="sz": op = operators.get_sz(self)
      elif name=="density": op = operators.index(self,n=range(n))
      else: raise
      ops = [o@op for o in ops] # define operators
      return spectrum.ev(self,operator=ops,**kwargs).real
      
#    from .magnetism import get_magnetization
#    return get_magnetization(self,nkp=nkp)
  def get_1dh(self,k=0.0):
      """Return a 1d Hamiltonian"""
      if self.is_multicell: raise # not implemented
      if not self.dimensionality==2: raise # not implemented
      intra,inter = kchain(self,k) # generate intra and inter
      hout = self.copy() # copy the Hamiltonian
      hout.intra = intra # store
      hout.inter = inter # store
      hout.dimensionality = 1 # one dimensional
      hout.geometry.dimensionality = 1 # one dimensional
      return hout
  def get_multicell(self):
    """Return a multicell Hamiltonian"""
    return multicell.turn_multicell(self)
  def turn_multicell(self):
    """Conver to multicell Hamiltonian"""
    h = multicell.turn_multicell(self)
    self.is_multicell = True
    self.hopping = h.hopping # assign hopping
  def get_no_multicell(self):
    """Return a multicell Hamiltonian"""
    return multicell.turn_no_multicell(self)
  def clean(self):
    """Clean a Hamiltonian"""
    from .clean import clean_hamiltonian
    clean_hamiltonian(self)
  def get_operator(self,name,projector=False):
    """Return a certain operator"""
    if name=="sx": return operators.get_sx(self)
    elif name=="sy": return operators.get_sy(self)
    elif name=="sz": return operators.get_sz(self)
    elif name=="current": 
        if self.dimensionality==1: return operators.get_current(self)
        else: raise
    elif name=="sublattice": return operators.get_sublattice(self,mode="both")
    elif name=="sublatticeA": return operators.get_sublattice(self,mode="A")
    elif name=="sublatticeB": return operators.get_sublattice(self,mode="B")
    elif name=="interface": return operators.get_interface(self)
    elif name=="spair": return operators.get_pairing(self,ptype="s")
    elif name=="deltax": return operators.get_pairing(self,ptype="deltax")
    elif name=="deltay": return operators.get_pairing(self,ptype="deltay")
    elif name=="deltaz": return operators.get_pairing(self,ptype="deltaz")
    elif name=="electron": return operators.get_electron(self)
    elif name=="hole": return operators.get_hole(self)
    elif name=="zposition": return operators.get_zposition(self)
    elif name=="yposition": return operators.get_yposition(self)
    elif name=="xposition": return operators.get_xposition(self)
    elif name=="velocity": return operators.get_velocity(self)
    # total magnetizations
    elif name=="mx": 
      return self.get_operator("sx")*self.get_operator("electron")
    elif name=="my": 
      return self.get_operator("sy")*self.get_operator("electron")
    elif name=="mz": 
      return self.get_operator("sz")*self.get_operator("electron")
    elif name=="valley": return operators.get_valley(self,projector=projector)
    elif name=="inplane_valley": return operators.get_inplane_valley(self)
    elif name=="valley_upper": 
      print("This operator only makes sense for TBG")
      ht = self.copy()
      ht.geometry.sublattice = self.geometry.sublattice * (np.sign(self.geometry.z)+1.0)/2.0
      return operators.get_valley(ht)
    elif name=="inplane_valley_upper": 
      print("This operator only makes sense for TBG")
      ht = self.copy()
      ht.geometry.sublattice = self.geometry.sublattice * (np.sign(self.geometry.z)+1.0)/2.0
      return operators.get_inplane_valley(ht)
    elif name=="valley_lower": 
      print("This operator only makes sense for TBG")
      ht = self.copy()
      ht.geometry.sublattice = self.geometry.sublattice * (-np.sign(self.geometry.z)+1.0)/2.0
      return operators.get_valley(ht)
    elif name=="ipr": return operators.ipr 
    else: raise
  def extract(self,name="mz"): 
    """Extract somethign from the Hamiltonian"""
    if self.has_eh: raise # not implemented
    if name=="density":
      return extract.onsite(self.intra,has_spin=self.has_spin)
    elif name=="mx" and self.has_spin:
      return extract.mx(self.intra)
    elif name=="my" and self.has_spin:
      return extract.my(self.intra)
    elif name=="mz" and self.has_spin:
      return extract.mz(self.intra)
    else: raise
  def write_magnetization(self,nrep=5):
    """Extract the magnetization and write it in a file"""
    if not self.has_eh: # without electron hole
      if self.has_spin: # for spinful
        mx = extract.mx(self.intra)
        my = extract.my(self.intra)
        mz = extract.mz(self.intra)
        g = self.geometry
        g.write_profile(mx,name="MX.OUT",normal_order=True,nrep=nrep)
        g.write_profile(my,name="MY.OUT",normal_order=True,nrep=nrep)
        g.write_profile(mz,name="MZ.OUT",normal_order=True,nrep=nrep)
        # this is just a workaround
        m = np.genfromtxt("MX.OUT").transpose()
        (x,y,z,mx) = m[0],m[1],m[2],m[3]
        my = np.genfromtxt("MY.OUT").transpose()[3]
        mz = np.genfromtxt("MZ.OUT").transpose()[3]
        np.savetxt("MAGNETISM.OUT",np.array([x,y,z,mx,my,mz]).T)
        return np.array([mx,my,mz])
#    return np.array([mx,my,mz]).transpose()
  def write_onsite(self,nrep=5,normal_order=False):
      """Extract onsite energy"""
      if self.has_eh: raise
      d = extract.onsite(self.intra,has_spin=self.has_spin)
      d = d - np.mean(d)
      self.geometry.write_profile(d,name="ONSITE.OUT",
              normal_order=normal_order,nrep=nrep)
  def write_hopping(self,**kwargs):
      groundstate.hopping(self,**kwargs)
  def get_ipr(self,**kwargs):
      """Return the IPR"""
      from . import ipr
      if self.dimensionality==0:
          return ipr.ipr(self.intra,**kwargs) 
      else: raise # not implemented








def get_first_neighbors(r1,r2,optimal=optimal):
  """Gets the fist neighbors, input are arrays"""
  if optimal:
    from . import neighbor
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
  elif h.dimensionality==0:  
    intra = h.intra  # assign intraterm
    energies = lg.eigvalsh(intra) # get eigenvalues of the current hamiltonian
    return (range(len(intra)),energies) # return indexes and energies
  else: raise




def diagonalize_hk(k):
  return lg.eigh(hk(k))





def eigenvectors(h,nk=10,kpoints=False,k=None,sparse=False,numw=None):
  import scipy.linalg as lg
  from scipy.sparse import csc_matrix as csc
  shape = h.intra.shape
  if h.dimensionality==0:
    vv = algebra.eigh(h.intra)
    vecs = [v for v in vv[1].transpose()]
    if kpoints: return vv[0],vecs,[[0.,0.,0.] for e in vv[0]]
    else: return vv[0],vecs
  elif h.dimensionality>0:
    f = h.get_hk_gen()
    if k is None: 
      from .klist import kmesh
      kp = kmesh(h.dimensionality,nk=nk) # generate a mesh
    else:  kp = np.array([k]) # kpoint given on input
#    vvs = [lg.eigh(f(k)) for k in kp] # diagonalize k hamiltonian
    nkp = len(kp) # total number of k-points
    if sparse: # sparse Hamiltonians
      vvs = [slg.eigsh(csc(f(k)),k=numw,which="LM",sigma=0.0,tol=1e-10) for k in kp] # 

    else: # dense Hamiltonians
      from . import parallel
      if parallel.cores>1: # in parallel
#        vvs = parallel.multieigh([f(k) for k in kp]) # multidiagonalization
        vvs = parallel.pcall(lambda k: algebra.eigh(f(k)),kp)
      else: vvs = [algebra.eigh(f(k)) for k in kp] # 
    nume = sum([len(v[0]) for v in vvs]) # number of eigenvalues calculated
    eigvecs = np.zeros((nume,h.intra.shape[0]),dtype=np.complex) # eigenvectors
    eigvals = np.zeros(nume) # eigenvalues

    #### New way ####
#    eigvals = np.array([iv[0] for iv in vvs]).reshape(nkp*shape[0],order="F")
#    eigvecs = np.array([iv[1].transpose() for iv in vvs]).reshape((nkp*shape[0],shape[1]),order="F")
#    if kpoints: # return also the kpoints
#      kvectors = [] # empty list
#      for ik in kp: 
#        for i in range(h.intra.shape[0]): kvectors.append(ik) # store
#      return eigvals,eigvecs,kvectors
#    else:
#      return eigvals,eigvecs

    #### Old way, slightly slower but clearer ####
    iv = 0
    kvectors = [] # empty list
    for ik in range(len(kp)): # loop over kpoints
      vv = vvs[ik] # get eigenvalues and eigenvectors
      for (e,v) in zip(vv[0],vv[1].transpose()):
        eigvecs[iv] = v.copy()
        eigvals[iv] = e.copy()
        kvectors.append(kp[ik])
        iv += 1
    if kpoints: # return also the kpoints
#      for iik in range(len(kp)):
#        ik = kp[iik] # store kpoint 
#        for e in vvs[iik][0]: kvectors.append(ik) # store
      return eigvals,eigvecs,kvectors
    else:
      return eigvals,eigvecs
  else:
    raise






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
  print("Intracell matrix")
  print(intra)
  print("Intercell matrix")
  print(inter)
  return



def add_sublattice_imbalance(h,mass):
  """ Adds to the intracell matrix a sublattice imbalance """
  if h.geometry.has_sublattice:  # if has sublattice
    def ab(i): 
      return h.geometry.sublattice[i]
  else: 
    print("WARNING, no sublattice present")
    return 0. # if does not have sublattice
  natoms = len(h.geometry.r) # assume spinpolarized calculation 
  rows = range(natoms)
  if callable(mass):  # if mass is a function
    r = h.geometry.r
    data = [mass(r[i])*ab(i) for i in range(natoms)]
  else: data = [mass*ab(i) for i in range(natoms)]
  massterm = csc_matrix((data,(rows,rows)),shape=(natoms,natoms)) # matrix
  h.intra = h.intra + h.spinless2full(massterm)








def build_eh_nonh(hin,c1=None,c2=None):
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
  if not c1 is None: # if there is coupling
    for i in range(n):
      for j in range(n):
        # couples electron in i with hole in j
        hout[2*i,2*j+1] = c1[i,j]  # electron hole term
  if not c2 is None: # if there is coupling
    for i in range(n):
      for j in range(n):
        # couples hole in i with electron in j
        hout[2*j+1,2*i] = np.conjugate(c2[i,j])  # hole electron term
  return hout 









def set_finite_system(hin,periodic=True):
  """ Transforms the hamiltonian into a finite system,
  removing the hoppings """
  from copy import deepcopy
  h = hin.copy() # copy Hamiltonian
  h.dimensionality = 0 # put dimensionality = 0
  if periodic: # periodic boundary conditions
    if h.dimensionality == 1:
      h.intra = h.intra + h.inter + h.inter.H 
    if h.dimensionality == 2:
      h.intra = h.intra +  h.tx + h.tx.H 
      h.intra = h.intra +  h.ty + h.ty.H 
      h.intra = h.intra +  h.txy + h.txy.H 
      h.intra = h.intra +  h.txmy + h.txmy.H 
  else: pass
  return h
  



def des_spin(m,component=0):
  """ Removes the spin degree of freedom"""
  d = len(m) # dimension of the matrix
  if d%2==1: # if the hamiltonian doesn't have the correct dimension
    print("Hamiltonian dimension is odd")
    raise
  mout = np.matrix([[0.0j for i in range(d//2)] for j in range(d//2)])
  for i in range(d//2):
    for j in range(d//2):
      mout[i,j] = m[2*i+component,2*j+component]  # assign spin up part
  return mout


def shift_fermi(h,fermi):
  """ Moves the fermi energy of the system, the new value is at zero"""
  r = h.geometry.r # positions
  n = len(r) # number of sites
  if checkclass.is_iterable(fermi): # iterable
    if len(fermi)==n: # same number of sites
      h.intra = h.intra + h.spinless2full(sparse_diag([fermi],[0]))
    else: raise
  else:
    rc = [i for i in range(n)]  # index
    datatmp = [] # data
    for i in range(n): # loop over positions 
      if callable(fermi): fshift = fermi(r[i]) # fermi shift
      else: fshift = fermi # assume it is a number
      datatmp.append(fshift) # append value
    m = csc_matrix((datatmp,(rc,rc)),shape=(n,n)) # matrix with the shift
    h.intra = h.intra + h.spinless2full(m) # Add matrix 
    return



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
        print("No hermitic element", i,j,m[i,j],m[j,i])
        return False
  return True
  





def first_neighborsnd(h):
  """ Gets a first neighbor hamiltonian"""
  r = h.geometry.r    # x coordinate 
  g = h.geometry
# first neighbors hopping, all the matrices
  a1, a2 = g.a1, g.a2
  def gett(r1,r2):
    """Return hopping given two sets of positions"""
    from . import neighbor
    pairs = neighbor.find_first_neighbor(r1,r2)
    if len(pairs)==0: rows,cols = [],[]
    else: rows,cols = np.array(pairs).T # transpose
    data = np.array([1. for c in cols])
    n = len(r1)
    m = csc_matrix((data,(rows,cols)),shape=(n,n),dtype=np.complex)
    m = h.spinless2full(m) # add spin degree of freedom if necessary
    if h.is_sparse: return m
    else: return m.todense()
  if h.dimensionality==0:
    h.intra = gett(r,r)
  elif h.dimensionality==1:
    h.intra = gett(r,r)
    h.inter = gett(r,r+a1)
  elif h.dimensionality==2:
    h.intra = gett(r,r)
    h.tx = gett(r,r+a1)
    h.ty = gett(r,r+a2)
    h.txy = gett(r,r+a1+a2)
    h.txmy = gett(r,r+a1-a2)



from .superconductivity import add_swave
from .superconductivity import build_eh
nambu_nonh = build_eh_nonh



from .bandstructure import lowest_bands



def hk_gen(h):
  """ Returns a function that generates a k dependent hamiltonian"""
  if h.dimensionality == 0: return lambda x: h.intra
  elif h.dimensionality == 1: 
    def hk(k):
      """k dependent hamiltonian, k goes from 0 to 1"""
      try: kp = k[0]
      except: kp = k
      tk = h.inter * h.geometry.bloch_phase([1.],kp) # get the bloch phase
      ho = h.intra + tk + tk.H # hamiltonian
      return ho
    return hk  # return the function
  elif h.dimensionality == 2: 
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
#        tk = p[0]*np.exp(1j*np.pi*2.*(p[1].dot(k)))  # add bloch hopping
        tk = p[0]*h.geometry.bloch_phase(p[1],k)  # add bloch hopping
        ho = ho + tk + tk.H  # add bloch hopping
      return ho
    return hk
  else: raise



from .neighbor import parametric_hopping
from .neighbor import parametric_hopping_spinful
from .neighbor import generate_parametric_hopping




from .superconductivity import get_nambu_tauz
from .superconductivity import project_electrons
from .superconductivity import project_holes
from .superconductivity import get_eh_sector_odd_even




def kchain(h,k):
  """ Return the kchain Hamiltonian """
  if h.dimensionality != 2: raise
  tky = h.ty*np.exp(1j*np.pi*2.*k)
  tkxy = h.txy*np.exp(1j*np.pi*2.*k)
  tkxmy = h.txmy*np.exp(-1j*np.pi*2.*k)  # notice the minus sign !!!!
  # chain in the x direction
  ons = h.intra + tky + tky.H  # intra of k dependent chain
  hop = h.tx + tkxy + tkxmy  # hopping of k-dependent chain
  return (ons,hop)



# import the function written in the library
from .kanemele import generalized_kane_mele




def turn_nambu(self):
  """Turn a Hamiltonian an Nambu Hamiltonian"""
  from .superconductivity import build_eh as nambu # redefine
  if self.has_eh: return # do nothing if already has eh
#  self.get_eh_sector = get_eh_sector_odd_even # assign function
#  if not self.has_eh: # if has not been assigned yet
#    self.nambu_tauz = get_nambu_tauz(self.intra) # assign matrix
  # add pairing
  self.intra = nambu(self.intra,is_sparse=self.is_sparse)
  if self.is_multicell: # for multicell hamiltonians
    for i in range(len(self.hopping)): # loop over hoppings
      self.hopping[i].m = nambu(self.hopping[i].m,is_sparse=self.is_sparse) # put in nambu form
  else: # conventional way
    if self.dimensionality==0: pass # one dimensional systems
    elif self.dimensionality==1: # one dimensional systems
      self.inter = nambu(self.inter,is_sparse=self.is_sparse)
    elif self.dimensionality==2: # two dimensional systems
      self.tx = nambu(self.tx,is_sparse=self.is_sparse)
      self.ty = nambu(self.ty,is_sparse=self.is_sparse)
      self.txy = nambu(self.txy,is_sparse=self.is_sparse)
      self.txmy = nambu(self.txmy,is_sparse=self.is_sparse)
    else: raise
  self.has_eh = True


from . import inout



def load(input_file="hamiltonian.pkl"):  return inout.load(input_file)


def print_hopping(h):
    """Print all the hoppings in a user friendly way"""
    from pandas import DataFrame
    def pprint(m): 
        if np.max(np.abs(m.real))>0.00001: print(DataFrame(m.real))
        if np.max(np.abs(m.imag))>0.00001: print(DataFrame(m.imag*1j))
        print("\n")
    print("Onsite")
    pprint(h.intra)
    if h.dimensionality==0: return
    h = h.get_multicell()
    for t in h.hopping:
        print("Hopping",t.dir)
        pprint(t.m)







