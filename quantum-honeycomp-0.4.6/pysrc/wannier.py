# routines to work with wannier hamiltonians

import numpy as np
import hamiltonians
import geometry
from scipy.sparse import csc_matrix as csc
from scipy.sparse import bmat


def read_geometry(input_file="wannier.win"):
  """Reads the geometry of the wannier calculation"""
  ll = read_between("begin unit_cell_cart","end unit_cell_cart",input_file)
  a1 = ll[1].split()
  a2 = ll[2].split()
  a3 = ll[3].split()
  # read the unit vectors
  a1 = np.array([float(a1[0]),float(a1[1]),float(a1[2])]) # first vector
  a2 = np.array([float(a2[0]),float(a2[1]),float(a2[2])]) # second vector
  a3 = np.array([float(a3[0]),float(a3[1]),float(a3[2])]) # second vector
  g = geometry.geometry()
  g.dimensionality = 2
  g.has_spin = False
  g.has_sublattice = False
  g.a1 = a1  # store vector
  g.a2 = a2  # store vector
  # read the coordinates 
  ll = read_between("begin projections","end projections",input_file)
  rs = [] # empty list for positions
  for l in ll:
    name = l.split(":")[0] # get name of hte atom
    r = get_positions(name,input_file) # get positins of the atoms
    for i in range(len(r)): # to real coordinates
      r[i] = r[i][0]*a1 + r[i][1]*a2 + r[i][2]*a3
    rs += r # store positions
  g.r = np.array(rs) # store in the class
  g.r2xyz() # store also in xyz atributes
  return g


def get_positions(atom,input_file):
  """Get positions of certain orbitals"""
  ll = read_between("begin atoms_frac","end atoms_frac",input_file)
  rs = [] # empty list
  for l in ll: # loop over lines
    l = l.split()
    name = l[0]
    if atom==name: # found atom
      r = np.array([float(l[1]),float(l[2]),float(l[3])]) # position
      rs.append(r) # add to the list
  return rs # return positions


def get_projected_atom_names(input_file):
  """Get the names of all the atoms who have projections"""
  lorb = read_between("begin projections","end projections",input_file)
  names = []
  for l in lorb:
    name = l.split(":")[0]
    if name not in names:
      names.append(name) # add this name
  return names


def get_orbitals(specie,input_file="wannier.win"):
  """ Get the names of the orbitals of a cetain atom"""
  lorb = read_between("begin projections","end projections",input_file)
  orbs = []
  for l in lorb: # loop over orbitals
    sname = l.split(":")[0] # name of the atom
    oname = l.split(":")[1] # name of the orbital
    oname = oname.split()[0] # remove the \n
    if specie==sname: orbs.append(oname)
  return orbs



def get_atoms_specie(specie,input_file="wannier.win"):
  """ Get number of atoms o a certain specie"""
  ll = read_between("begin atoms_frac","end atoms_frac",input_file)
  nat = 0
  for l in ll:
    name = l.split()[0] # name of the atom
    if name==specie: nat +=1 
  return nat # return umber of atoms






def get_index_orbital(specie,atom,orbital,input_file="wannier.win"):
  """Get the index of a certain orbital"""
  ll = read_between("begin atoms_frac","end atoms_frac",input_file)
  lorb = read_between("begin projections","end projections",input_file)
  iorb = 0 # index in the matrix
  for l in lorb: # loop over orbitals
    sname = l.split(":")[0] # name of the atom
    oname = l.split(":")[1] # name of the orbital
    oname = oname.split()[0] # remove the \n
    ifound = 0 # number of atoms found
    for la in ll: # loop over atoms in the structure
      rname = la.split()[0] # name of the atom in the structure
      if rname==sname: 
        iorb += 1 # increase de counter of the hamiltonian index 
        ifound += 1 # increase the counter in the atoms
      # if the desired orbital and atom has been reached
      if specie==sname and atom==ifound and orbital==oname:
        return iorb
  raise # error if this point is reached



def get_indexes(input_file="wannier.win"):
  """Returns a list with the indexes of the different orbitals for
  each atom """
  ll = read_between("begin atoms_frac","end atoms_frac",input_file)
  lorb = read_between("begin projections","end projections",input_file)
  # first get all the different names
  names = get_projected_atom_names(input_file) # get names of the atoms
  dat = dict() # create a diccionary
  for name in names: # loop over different atoms
    out = [] # empty list with the different indexes and names
    ind = 0 # initialize index of the hamiltonian
    for l in ll: # loop over atoms in the structure
      out_atom = [] # empty list for the next atom
      l = l.split() # split the line
      label = l[0]  # get label of the atom
      if label==name: # if same name as input
        for lo in lorb: # go to the orbitals list and look for this atom
          name_o = lo.split(":")[0] # get name of the atom in the orbital list
          if name_o==name: # if an orbital of this atom has been found 
            orb = lo.split(":")[1] # name of the orbital
            orb = orb.split()[0] # name of the orbital
            out_atom.append((ind,orb)) # append tuple with index and norm
        out.append([out_atom]) # append this atom to the final 
    dat[name] = out # dd this atom to the diccionary
  return dat  # return the diccionary




def read_between(a,b,input_file):
  ll = open(input_file).readlines()
  start = False # found the klist
  out = []
  for (i,l) in zip(range(len(ll)),ll):
    if b in l: break # end of klist
    if start: # sotre line
      out.append(l)
    if a in l: start = True # found beginning 
  return out # return output lines





def read_hamiltonian(input_file="hr_truncated.dat",is_real=False):
  """Reads an output hamiltonian from wannier"""
  mt = np.genfromtxt(input_file) # get file
  m = mt.transpose() # transpose matrix
  # read the hamiltonian matrices
  class Hopping: pass # create empty class
  tlist = []
  def get_t(i,j,k):
    norb = np.max([np.max(np.abs(m[3])),np.max(np.abs(m[4]))])
    mo = np.matrix(np.zeros((norb,norb),dtype=np.complex))  
    for l in mt: # look into the file
      if i==int(l[0]) and j==int(l[1]) and k==int(l[2]):
        if is_real:
          mo[int(l[3])-1,int(l[4])-1] = l[5] # store element
        else:
          mo[int(l[3])-1,int(l[4])-1] = l[5] + 1j*l[6] # store element
    return mo # return the matrix
#  for i in range(-nmax,nmax):
#    for j in range(-nmax,nmax):
#      for k in range(-nmax,nmax):
#        t = Hopping() # create hopping
#        t.dir = [i,j,k] # direction
#        t.m = get_t(i,j,k) # read the matrix
#        tlist.append(t) # store hopping
  # the previous is not used yet...
  g = geometry.kagome_lattice() # create geometry
  h = g.get_hamiltonian() # build hamiltonian
  h.intra = get_t(0,0,0)
  h.tx = get_t(1,0,0)
  h.ty = get_t(0,1,0)
  h.txy = get_t(1,1,0)
  h.txmy = get_t(1,-1,0)
  h.has_spin = False  # if not spin polarized
  h.geometry = read_geometry() # read the geometry of the system
  if len(h.geometry.r)!=len(h.intra): 
    print "Dimensions do not match",len(g.r),len(h.intra)
    print h.geometry.r
    raise # error if dimensions dont match
  return h



def read_multicell_hamiltonian(input_file="hr_truncated.dat",ncells=[3,3,0]):
  """Reads an output hamiltonian from wannier"""
  mt = np.genfromtxt(input_file) # get file
  m = mt.transpose() # transpose matrix
  # read the hamiltonian matrices
  class Hopping: pass # create empty class
  tlist = []
  def get_t(i,j,k):
    norb = np.max([np.max(np.abs(m[3])),np.max(np.abs(m[4]))])
    mo = np.matrix(np.zeros((norb,norb),dtype=np.complex))  
    for l in mt: # look into the file
      if i==int(l[0]) and j==int(l[1]) and k==int(l[2]):
        mo[int(l[3])-1,int(l[4])-1] = l[5] + 1j*l[6] # store element
    return mo # return the matrix
  for i in range(-ncells[0],ncells[0]+1):
    for j in range(-ncells[1],ncells[1]+1):
      for k in range(-ncells[2],ncells[2]+1):
        if (i,j,k)==(0,0,0): continue # skip intracell
        print i,j,k
        t = Hopping() # create hopping
        t.dir = [i,j,k] # direction
        t.m = get_t(i,j,k) # read the matrix
        tlist.append(t) # store hopping
  # the previous is not used yet...
  g = geometry.kagome_lattice() # create geometry
  h = g.get_hamiltonian() # build hamiltonian
  h.is_multicell = True
  h.hopping = tlist # list of hoppings
  h.has_spin = False  # if not spin polarized
  h.geometry = read_geometry() # read the geometry of the system
  h.intra = get_t(0,0,0)
  if len(h.geometry.r)!=len(h.intra): 
    print "Dimensions do not match",len(g.r),len(h.intra)
    print h.geometry.r
    raise # error if dimensions dont match
  return h







def get_klist(input_file="wannier.win"):
  """ Get the klist for bands calculation"""
  ll = read_between("begin kpoint_path","end kpoint_path",input_file)
  kp = [] # empty vertex
  for l in ll: 
    l2 = l.split() # split the numbers
    kp.append([[float(l2[1]),float(l2[2]),float(l2[3])],[float(l2[5]),float(l2[6]),float(l2[7])]])
  klist = [] # empty klist
  nk = 500/len(kp) # number of kpoints
  for (k1,k2) in kp: # loop over pairs
    k1 = np.array(k1)
    k2 = np.array(k2)
    dk = (k2-k1)/nk # create dk
    for i in range(nk): # loop over ks
      klist.append(k1+i*dk)
  return klist


def get_num_wannier(input_file):
  """Get the number of wannier orbitals"""
  ll = read_between("begin atoms_frac","end atoms_frac",input_file)
  lorb = read_between("begin projections","end projections",input_file)
  norb = 0
  for o in lorb: # loop over orbitals
    name = o.split(":")[0] # name of the orbital
    for l in ll:
      rname = l.split()[0] # name of the atom
      if name==rname: norb += 1
  return norb # return number of orbitals



def get_num_atoms(specie,input_file):
  """Get the number of wannier orbitals"""
  ll = read_between("begin atoms_frac","end atoms_frac",input_file)
  nat = 0 # initialize
  for l in ll:
    name = l.split()[0] # name of the atom
    if name==specie: nat += 1
  return nat # return number of orbitals






def generate_soc(specie,value,input_file="wannier.win",nsuper=1):
  """Add SOC to a hamiltonian based on wannier.win"""
  fo = open(".soc.status","w")
  iat = 1 # atom counter
  orbnames = names_soc_orbitals(specie) # get which are the orbitals
  ls = soc_l((len(orbnames)-1)/2) # get the SOC matrix
  norb = get_num_wannier(input_file) # number of wannier orbitals
  m = np.matrix([[0.0j for i in range(norb*2)] for j in range(norb*2)]) # matrix
  nat = get_num_atoms(specie,input_file) # number of atoms of this specie
  for iat in range(nat):
   # try:
      fo.write("Attempting "+specie+"  "+str(iat+1)+"\n")
      for i in range(len(orbnames)): # loop over one index
        orbi = orbnames[i] 
        for j in range(len(orbnames)): # loop over other index
          orbj = orbnames[j]
          ii = get_index_orbital(specie,iat+1,orbi)  # index in wannier
          jj = get_index_orbital(specie,iat+1,orbj) # index in wannier
          fo.write(str(ii)+"   "+str(jj)+"\n")
          ii += -1 # python starts in 0
          jj += -1 # python starts in 0
          m[2*ii,2*jj] = ls[2*i,2*j] # store the soc coupling
          m[2*ii+1,2*jj] = ls[2*i+1,2*j] # store the soc coupling
          m[2*ii,2*jj+1] = ls[2*i,2*j+1] # store the soc coupling
          m[2*ii+1,2*jj+1] = ls[2*i+1,2*j+1] # store the soc coupling
   #   return
   # except: break
  fo.close()
  n = nsuper**2 # supercell
  mo = [[None for i in range(n)] for j in range(n)]
  for i in range(n): mo[i][i] = csc(m) # diagonal elements
  mo = bmat(mo).todense() # dense matrix
  return mo*value # return matrix
#  for name in atoms: # loop over atoms

def get_element(name):
  """Gets the true name of an element, assumes that there might be a number"""
  out = "" # initialize
  if name[-1] in [str(i) for i in range(10)]:
    for i in range(len(name)-1): # except the last one
      out += name[i] # add to the string
    return get_element(out) # recall the function, just in case two digits
  return name



def names_soc_orbitals(specie):
  name = get_element(specie) # get the true name of the atom
  dorbs = ["dz2","dxz","dyz","dxy","dx2-y2"]
  datoms = ["Sc","Ti","V","Cr","Mn","Fe","Co","Ni","Cu","Zn"]
  datoms += ["Y","Zr","Nb","Mo","Tc","Ru","Rh","Pd","Ag","Cd"]
  datoms += ["Hf","Ta","W","Re","Os","Ir","Pt","Au","Hg"]
  if name in datoms: return dorbs # for dorbitals
  porbs = ["px","py","pz"]
  patoms = ["B","C","N","O","F","P"]
  patoms += ["Al","Si","P","S","Cl"]
  patoms += ["Ga","Ge","As","Se","Br"]
  patoms += ["In","Sn","An","Te","I"]
  if name in patoms: return porbs # for dorbitals




def ylm2xyz_l2():
  """Return the matrix that converts the cartesian into spherical harmonics"""
  m = np.matrix([[0.0j for i in range(5)] for j in range(5)]) 
  s2 = np.sqrt(2.)
  m[2,0] = 1. # dz2
  m[1,1] = 1./s2 # dxz
  m[3,1] = -1./s2 # dxz
  m[1,2] = 1j/s2 # dyz
  m[3,2] = 1j/s2 # dyz
  m[0,3] = 1j/s2 # dxy
  m[4,3] = -1j/s2 # dxy
  m[0,4] = 1./s2 # dx2y2
  m[4,4] = 1./s2 # dx2y2
#  m = m.H # no fucking idea if there is a transpose...
#  m = m.T # inverse
  from increase_hilbert import spinful
  m = spinful(m) # with spin degree of freedom
  return m # return change of bassi matrix



def ylm2xyz_l1():
  """Return the matrix that converts the cartesian into spherical harmonics"""
  m = np.matrix([[0.0j for i in range(3)] for j in range(3)]) 
  s2 = np.sqrt(2.)
  m[1,2] = 1. # pz
  m[0,0] = 1./s2 # dxz
  m[2,0] = -1./s2 # dxz
  m[0,1] = 1j/s2 # dyz
  m[2,1] = 1j/s2 # dyz
  from increase_hilbert import spinful
  m = spinful(m) # with spin degree of freedom
  return m








def soc_l(l):
  """Calculate the spin orbit coupling in a basis of spherical harmonics"""
  nm = 2*l + 1 # number of components
  zero = np.matrix([[0.0j for i in range(2*nm)] for j in range(2*nm)]) 
  # initialize matrices
  lz = zero.copy()
  lx = zero.copy()
  ly = zero.copy()
  lm = zero.copy()
  lp = zero.copy()
  # create l+ and l- and lz
  for m in range(-l,l): # loop over m components
    val = np.sqrt((l-m)*(l+m+1)) # value of the cupling
    im = m + l
    lp[2*(im+1),2*im] = val # up channel
    lp[2*(im+1)+1,2*im+1] = val # down channel
  for m in range(-l,l+1): # loop over m components
    im = m + l
    lz[2*im,2*im] = m # value of lz, up channel
    lz[2*im+1,2*im+1] = m # value of lz, down channel
  lm = lp.H # adjoint
  lx = (lp + lm) /2.
  ly = -1j*(lp - lm) /2.
  # create spin matrices
  sz = zero.copy()
  sx = zero.copy()
  sy = zero.copy()
  for m in range(-l,l+1): # loop over m components
    im = m + l
    sx[2*im,2*im+1] = 1.0 
    sx[2*im+1,2*im] = 1.0 
    sy[2*im,2*im+1] = -1j 
    sy[2*im+1,2*im] = 1j 
    sz[2*im,2*im] = 1.
    sz[2*im+1,2*im+1] = -1.
  # check that the matrix is fine
  sx = sx/2.
  sy = sy/2.
  sz = sz/2.
  if True:
    comm_zero(sx,lx)
    comm_zero(sx,ly)
    comm_zero(sx,lz)
    comm_zero(sy,ly)
    comm_zero(sz,lz)
    comm_angular(sx,sy,sz)
    comm_angular(sy,sz,sx)
    comm_angular(sz,sx,sy)
    comm_angular(lx,ly,lz)
    comm_angular(ly,lz,lx)
    comm_angular(lz,lx,ly)
  ls = lx*sx + ly*sy + lz*sz  # SOC matrix
  import scipy.linalg as lg
  from scipy.sparse import csc_matrix as csc
  if l==2: R = ylm2xyz_l2() # get change of basis matrix
  if l==1: R = ylm2xyz_l1() # get change of basis matrix
  ls = R.H * ls * R # change to cartesian orbitals
#  print csc(ls-ls.H)
#  print lg.eigvalsh(ls)
  return ls # return the matrix



def comm_angular(x,y,z):
  from scipy.sparse import csc_matrix as csc
  xy = x*y - y*x
  xy = xy - 1j*z
  if np.abs(np.max(xy))>0.001:
    print csc(xy)
    raise



def comm_zero(x,y):
  xy = x*y - y*x
  if np.abs(np.max(xy))>0.01:
    print x*y - y*x
    raise



def symmetrize_atoms(h,specie,input_file="wannier.win"):
  """Symmetrizes a certain atom"""
  orbs = get_orbitals(specie,input_file=input_file) # read the orbitals
  nat = get_atoms_specie(specie,input_file=input_file) # number of atoms
  if h.has_spin: raise
  for iorb in orbs: # loop over orbitals
    avg = 0.
    for iat in range(nat): # loop over atoms 
      i = get_index_orbital(specie,iat+1,iorb) - 1 
      ons = h.intra[i,i] # add to the average
      avg += ons # add to the average
    avg = avg/nat # average value
    for iat in range(nat): # loop over atoms
      i = get_index_orbital(specie,iat+1,iorb) - 1
      h.intra[i,i] = avg # substitute by average


def get_hoppings(h,specie1,specie2,input_file="wannier.win"):
  """Get hoppings between two atoms"""
  orbs1 = get_orbitals(specie1,input_file=input_file) # read the orbitals
  orbs2 = get_orbitals(specie2,input_file=input_file) # read the orbitals
  nat1 = get_atoms_specie(specie1,input_file=input_file) # number of atoms
  nat2 = get_atoms_specie(specie2,input_file=input_file) # number of atoms
  mats = [h.intra,h.tx,h.ty,h.txy,h.txmy]
  for iorb1 in orbs1: # loop over orbitals
    for iorb2 in orbs2: # loop over orbitals
      for iat1 in range(nat1): # loop over atoms 
        for iat2 in range(nat2): # loop over atoms 
          i = get_index_orbital(specie1,iat1+1,iorb1) - 1 
          j = get_index_orbital(specie2,iat2+1,iorb1) - 1 
          for m in mats:
            ons = m[i,j].real # add to the average
            print iat1,iorb1,iat2,iorb2,ons



def get_atomic_projection(specie,input_file="wannier.win"):
  """Get the matrix that projects onto a certain atom"""
  orbs = get_orbitals(specie,input_file=input_file) # read the orbitals
  norb = get_num_wannier(input_file) # number of wannier orbitals
  proj = np.matrix([[0.0j for i in range(len(norb))] for j in range(len(norb))])
  for iorb in orbs: # loop over orbitals
    i = get_index_orbital(specie,iat+1,iorb) - 1 
    proj[i,i] = 1.0 # non vanishing
  return proj





def read_supercell_hamiltonian(input_file="hr_truncated.dat",is_real=False,nsuper=1):
  """Reads an output hamiltonian for a supercell from wannier"""
  mt = np.genfromtxt(input_file) # get file
  m = mt.transpose() # transpose matrix
  # read the hamiltonian matrices
  class Hopping: pass # create empty class
  tlist = []
  def get_t(i,j,k):
    print i,j,k
    norb = np.max([np.max(np.abs(m[3])),np.max(np.abs(m[4]))])
    mo = np.matrix(np.zeros((norb,norb),dtype=np.complex))  
    for l in mt: # look into the file
      if i==int(l[0]) and j==int(l[1]) and k==int(l[2]):
        if is_real:
          mo[int(l[3])-1,int(l[4])-1] = l[5] # store element
        else:
          mo[int(l[3])-1,int(l[4])-1] = l[5] + 1j*l[6] # store element
    return mo # return the matrix
  # this function will be called in a loop
  g = geometry.kagome_lattice() # create geometry
  h = g.get_hamiltonian() # build hamiltonian
  h.has_spin = False
  nstot = nsuper**2
  intra = [[None for i in range(nstot)] for j in range(nstot)]
  tx = [[None for i in range(nstot)] for j in range(nstot)]
  ty = [[None for i in range(nstot)] for j in range(nstot)]
  txy = [[None for i in range(nstot)] for j in range(nstot)]
  txmy = [[None for i in range(nstot)] for j in range(nstot)]
  from scipy.sparse import csc_matrix as csc
  vecs = []
  # create the identifacion vectors
  inds = []
  acu = 0
  for i in range(nsuper): # loop over first replica
    for j in range(nsuper): # loop over second replica
      vecs.append(np.array([i,j])) # append vector
      inds.append(acu)
      acu += 1 # add one to the accumulator
  for i in inds: # loop over first vector
    for j in inds:  # loop over second vector
      v1 = vecs[i] # position of i esim cell
      v2 = vecs[j] # position of j esim cell
      dv = v2 - v1 # difference in vector
      print dv
      # get the different block elements
      intra[i][j] = csc(get_t(dv[0],dv[1],0))
      tx[i][j] = csc(get_t(dv[0]+nsuper,dv[1],0))
      ty[i][j] = csc(get_t(dv[0],dv[1]+nsuper,0))
      txy[i][j] = csc(get_t(dv[0]+nsuper,dv[1]+nsuper,0))
      txmy[i][j] = csc(get_t(dv[0]+nsuper,dv[1]-nsuper,0))
  h.intra = bmat(intra).todense()
  h.tx = bmat(tx).todense()
  h.ty = bmat(ty).todense()
  h.txy = bmat(txy).todense()
  h.txmy = bmat(txmy).todense()
  h.geometry = read_geometry() # read the geometry of the system
  h.geometry = h.geometry.supercell(nsuper) # create supercell
  if len(h.geometry.r)!=len(h.intra): 
    print "Dimensions do not match",len(g.r),len(h.intra)
    print h.geometry.r
    raise # error if dimensions dont match
  return h




def rotate90(h):
  """ Rotate by 90 degrees th unit cell"""
  h.geometry.x, h.geometry.y = h.geometry.y, h.geometry.x
  h.geometry.xyz2r()
  h.geometry.a1, h.geometry.a2 = h.geometry.a2,  h.geometry.a1  # geometry
  h.tx , h.ty = h.ty, h.tx
  h.txy , h.tmy = h.txy, h.txmy.H
  return h


