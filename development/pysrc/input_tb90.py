from __future__ import print_function
from scipy.sparse import csc_matrix as csc
from scipy.sparse import coo_matrix as coo
import numpy as np
import os
# library to write inputs for tb90.x

class tb90in:
  """Class with all the inputs of tb90.in"""
  class electrons:
    filling = 0.5    # filling of the system
    shift_fermi = 0.0   # shift fermi energy on input
    extra_electrons = 0.0 # number of extra electrons
  class kpoints:
    nkpoints = 10   # number of kpoints
    expand_bz = 1.0 # expand the sampling
    refine_kmesh_bands = 10 # kmesh refination for bandstructure calculation
    refine_kmesh_dos = 10 # kmesh refination for DOS calculation
    refine_kmesh_dielectric = 1 # kmesh refination for dielectric calculation
    refine_kmesh_denmat = 1 # kmesh refination for bandstructure calculation
  class bands:
    bands_operators_option = 'Sz' # operator
    use_ewindow_bands = True # use energy window in bands calculation
    economic_memory_bands = False
    emin_bands = -1.0
    emax_bands = 1.0
    do_bands = True  # perform bands calculation
    klist_bands = "default"
  class scf_convergence:
    do_scf = False   # perform mean field calculation
    is_collinear = False   # perform mean field calculation
    mean_field_operators = 'hubbard'  # mean field operators
    mean_field_matrix = 'from_file'  # initial mean field guess
    save_mean_field = True
    mean_field_amplitude = 0.3  # mean field amplitude if generated
    hubbard_scf = 1.0 # hubbard term in SCF
    ultrafast_mf = False
    num_old_ham = 1   # number of old hamiltonians stored for SCF
    max_ite = 1000  # maximun number of iterations
    smearing = 0.00001  # smearing in SCF
    smearing_type_scf = "fermi-dirac"
    max_scf_err = 0.0001  # maximun error in SCF
    mix_coef = 1.1  # mixing coefficient in the SCF damping
    shift_to_zero = True # shift fermi energy to zero
  class dos:
    do_dos = False # perform dos calculation
    dos_operators_option = 'None' # default operators
    num_dos = 300   # number of intervals for the dos calculation
    define_num_dos = True # define the number of energy steps
    estep_dos = 0.02 # energy step in DOS calculation
    use_ewindow_dos = False # use energy window for the DOS calculation
    emin_dos = -1.0  # minimun energy
    emax_dos = 1.0 # maximun energy
    autosmearing_dos = False # autocalculate smearing
    smearing_dos = 0.04 # smearing in the DOS calculation
  class dielectric:
    do_dielectric = False # perform dielectric calculation
    dielectric_type = 'local_xy_spin' # type of dielectric calculation
    num_ene_chi = 100   # number of energies for dielectric
    ewindow_chi = 10.0  # energy window to accept any contribution to chi
    emin_chi = -0.5   # minimun energy for dielectric
    emax_chi = 0.5  # minimun energy for dielectric
    smearing_chi = 0.01 # smearing in dielectric response function
    temperature_chi = 0.001 # temperature in dielectric response function
    q_chi = 0.0001  # q-vector for the response
    phi_chi = 0.0  # phi-angle of q_chi for the response (2D)
  class density_matrix:
    do_denmat = False # perform density matrix ground state calculation
    write_denmat = False # write the full density matrix
  class development:
    show_debug = False # show debug information
    show_read_hamiltonian = False
  class topology:
    do_berry = False  # do berry stuff calculation
    has_tr = False  # time reversal symmetry
    nkpoints_connection = 10 # points in the path, highly unstable
    dk_becurv= 0.01  # delta k for berry curvature
    klist_berry = "default" # default klist



  def write(self): # write in file
    """ Write the tb90.in """
    write_tb90in(self)
  def run(self):
    """ Run tb90.x"""
    os.system("tb90.x")

  def set_hubbard(self,value):
    self.scf_convergence.hubbard_scf = value

  def mode(self,name):
    """Different default modes for tb90.x"""
    if name=="all_eigenvalues":
      self.bands.bands_operators_option = "None"
      self.bands.use_ewindow_bands = False
    elif name=="eigenvalues_small_window":
      self.bands.bands_operators_option = "None"
      self.bands.use_ewindow_bands = True
      self.bands.emin_bands = -0.1
      self.bands.emax_bands = 0.1
    elif name=="eigenvalues_Sz_window":
      self.bands.bands_operators_option = "Sz"
      self.bands.use_ewindow_bands = True
      self.bands.emin_bands = -1.0
      self.bands.emax_bands = 1.0
    elif name=="DOS":
      self.dos.do_dos = True
    elif name=="SCF":
      self.scf_convergence.do_scf = True
    elif name=="total_energy":
      self.scf_convergence.do_scf = True
      self.scf_convergence.hubbard_scf = 0.0
      self.scf_convergence.mean_field_amplitude = 0.0
    elif name=="nothing":
      self.dos.do_dos = False 
      self.scf_convergence.do_scf = False 
      self.dielectric.do_dielectric = False 
      self.bands.do_bands = False 
      self.density_matrix.do_denmat = False 
      self.topology.do_berry = False 
    elif name=="spin-RPA":
      self.dielectric.do_dielectric = True
      self.dielectric.dielectric_type = "local_xy_spin"
    elif name=="bands":
      self.bands.do_bands = True
    elif name=="berry":
      self.topology.do_berry = True
    elif name == "magnetism":
      self.density_matrix.do_denmat = True
    elif name == "total_energy":
      self.dos.do_dos = False 
      self.scf_convergence.do_scf = False 
      self.dielectric.do_dielectric = False 
      self.bands.do_bands = False 
      self.density_matrix.do_denmat = False 
      self.scf_convergence.do_scf = True   # activate SCF
      self.scf_convergence.max_scferr = 0.000001  # high precission
    else:
      print("Unimplemeted mode in tb90in")
      raise


def write_tb90in(p):
  """ Writes the file tb90.in """
  fin = open("tb90.in","w")
  def w(l,indent=True):
    if not indent:
      fin.write("\n") # write the indentation
    if indent:
      fin.write("   ") # write the indentation
    fin.write(l+"\n") # write the line
  s = stringtb90 # function to write parameters
  # write electrons
  nl = p.electrons # get the subclass
  w("&electrons",indent=False)
  w("filling = "+s(nl.filling))    # filling of the system
  w("shift_fermi = "+s(nl.shift_fermi))   # fermi energy on input
  w("extra_electrons = "+s(nl.extra_electrons)) # number of extra electrons
  w("/")

  # write kpoints
  nl = p.kpoints # get the subclass
  w("&kpoints",indent=False)
  w("nkpoints = "+s(nl.nkpoints))   # number of kpoints
  w("expand_bz = "+s(nl.expand_bz)) # expand the sampling
  w("refine_kmesh_bands = "+s(nl.refine_kmesh_bands))
  w("refine_kmesh_dos = "+s(nl.refine_kmesh_dos))
  w("refine_kmesh_dielectric = "+s(nl.refine_kmesh_dielectric))
  w("refine_kmesh_denmat = "+s(nl.refine_kmesh_denmat)) 
  w("/")

  # write bands
  nl = p.bands # get the subclass
  w("&bands",indent=False)
  w("bands_operators_option = "+s(nl.bands_operators_option)) # operator
  w("use_ewindow_bands = "+s(nl.use_ewindow_bands)) # use energy window in bands calculation
  w("emin_bands = "+s(nl.emin_bands))
  w("emax_bands = "+s(nl.emax_bands))
  w("do_bands = "+s(nl.do_bands))  # perform bands calculation
  w("klist_bands = "+s(nl.klist_bands))  # perform bands calculation
  w("economic_memory_bands = "+s(nl.economic_memory_bands))  # perform bands calculation
  w("/")

  # scf calculations
  nl = p.scf_convergence # get the subclass
  w("&scf_convergence",indent=False)
  w("do_scf = "+s(nl.do_scf))   
  w("ultrafast_mf = "+s(nl.ultrafast_mf))   
  w("is_collinear = "+s(nl.is_collinear))   
  w("mean_field_operators = "+s(nl.mean_field_operators))
  w("mean_field_matrix = "+s(nl.mean_field_matrix))
  w("save_mean_field = "+s(nl.save_mean_field))
  w("mean_field_amplitude = "+s(nl.mean_field_amplitude))
  w("hubbard_scf = "+s(nl.hubbard_scf))
  w("num_old_ham = "+s(nl.num_old_ham))
  w("max_ite =  "+s(nl.max_ite))
  w("smearing =  "+s(nl.smearing))
  w("smearing_type_scf =  "+s(nl.smearing_type_scf))
  w("max_scf_err = "+s(nl.max_scf_err))
  w("mix_coef = "+s(nl.mix_coef))
  w("shift_to_zero = "+s(nl.shift_to_zero))
  w("/")


  # dos calculations
  nl = p.dos # get the subclass
  w("&dos",indent=False)
  w("do_dos = "+s(nl.do_dos))
  w("dos_operators_option = "+s(nl.dos_operators_option))
  w("num_dos = "+s(nl.num_dos))
  w("define_num_dos = "+s(nl.define_num_dos))
  w("estep_dos = "+s(nl.estep_dos))
  w("use_ewindow_dos = "+s(nl.use_ewindow_dos))
  w("emin_dos = "+s(nl.emin_dos))
  w("emax_dos = "+s(nl.emax_dos))
  w("autosmearing_dos = "+s(nl.autosmearing_dos))
  w("smearing_dos = "+s(nl.smearing_dos))
  w("/")


  # dielectric calculation
  nl = p.dielectric # get the subclass
  w("&dielectric",indent=False)
  w("do_dielectric = " +s(nl.do_dielectric)) 
  w("dielectric_type = "+s(nl.dielectric_type))
  w("num_ene_chi = "+s(nl.num_ene_chi))
  w("ewindow_chi = "+s(nl.ewindow_chi)) 
  w("emin_chi = "+s(nl.emin_chi)) 
  w("emax_chi = "+s(nl.emax_chi)) 
  w("smearing_chi = "+s(nl.smearing_chi)) 
  w("temperature_chi = "+s(nl.temperature_chi)) 
  w("q_chi = "+s(nl.q_chi))
  w("phi_chi = "+s(nl.phi_chi))
  w("/")


  # density matrix calculation
  nl = p.density_matrix # get the subclass
  w("&density_matrix",indent=False)
  w("write_denmat = "+s(nl.write_denmat)) 
  w("do_denmat = "+s(nl.do_denmat))
  w("/")

  nl = p.topology # get the subclass
  w("&topology",indent=False)
  w("do_berry = "+s(nl.do_berry))
  w("has_tr = "+s(nl.has_tr)+"   ! z2 calculation for systes with TR")
  w("nkpoints_connection = "+s(nl.nkpoints_connection)+"  ! highly unstable")
  w("dk_becurv = "+s(nl.dk_becurv))
  w("klist_berry = "+s(nl.klist_berry))
  w("/")
  # development
  nl = p.development # get the subclass
  w("&development",indent=False)
  w("show_debug = "+s(nl.show_debug))
  w("/")


  fin.close()


def stringtb90(a):
  """Returns a string with the format for tb90"""
  if type(a) is bool:
    if a: return ".true."
    if not a: return ".false."
  if type(a) is str:
    return "'"+a+"'"
  if type(a) is float:
    o = str(a)
    o = o.replace("e","d")
    return o
  if type(a) is int:
    return str(a)





# write hamiltonian
def write_hamiltonian(h,output_file="hamiltonian.in"):
  ons = h.intra
  x = h.geometry.x
  y = h.geometry.y
  n=ons.shape[0]
  f=open(output_file,'w')
  f.write("DIMENSIONALITY_OF_THE_SYSTEM\n")
  f.write(str(h.dimensionality))
  f.write("\n\n\n")
  f.write('\nDIMENSION_OF_THE_HAMILTONIAN\n')
  f.write(str(n)+'\n')
  f.write('\n\n       ONSITE_MATRIX\n')
# write non vanishing elements
  nv=nv_el(ons)  
  for iv in range(len(nv)):
    f.write(str(int(nv[iv][0]))+'   ')
    f.write(str(int(nv[iv][1]))+'   ')
    f.write('{0:.8f}'.format(float(nv[iv][2]))+'   ')
    f.write('{0:.8f}'.format(float(nv[iv][3]))+'   ')
    f.write('  !!!  i  j   Real   Imag\n')
  f.write('\n')
  def write_hopping(hop,name):
    """Writes a particular hopping"""
    f.write('\n\n       '+name+"\n")
    nv=nv_el(hop)
    for iv in range(len(nv)):
      f.write(str(int(nv[iv][0]))+'   ')
      f.write(str(int(nv[iv][1]))+'   ')
      f.write('{0:.8f}'.format(float(nv[iv][2]))+'   ')
      f.write('{0:.8f}'.format(float(nv[iv][3]))+'   ')
      f.write('  !!!  i  j   Real   Imag\n')
  #################################
  # now write all the hoppings!!!!
  #################################
  if h.dimensionality == 1: # if one dimensional system
    hop = h.inter # read the hopping
    write_hopping(hop,"HOPPING_MATRIX_1") # write the first hopping
    write_hopping(hop.H,"HOPPING_MATRIX_-1") # write the second hopping
  if h.dimensionality == 2: # if two dimensional system
    write_hopping(h.tx,"HOPPING_MATRIX_1") # write the first hopping
    write_hopping(h.tx.H,"HOPPING_MATRIX_-1") # write the second hopping
    write_hopping(h.ty,"HOPPING_MATRIX_0_1") # write the second hopping
    write_hopping(h.ty.H,"HOPPING_MATRIX_0_-1") # write the second hopping
    write_hopping(h.txy,"HOPPING_MATRIX_1_1") # write the second hopping
    write_hopping(h.txy.H,"HOPPING_MATRIX_-1_-1") # write the second hopping
    write_hopping(h.txmy,"HOPPING_MATRIX_1_-1") # write the second hopping
    write_hopping(h.txmy.H,"HOPPING_MATRIX_-1_1") # write the second hopping
# print the coordinates
  f.write('\n\nPOSITIONS\n')
  for i in range(len(x)):
    for j in range(2):
      f.write('{0:.8f}'.format(x[i])+'  ')
      f.write('{0:.8f}'.format(y[i])+'  ')
      f.write('{0:.8f}'.format(0.)+'        ')
      f.write('\n')

  f.write('\n\n')
  # write lattice vectors
  if h.dimensionality == 1:  # if  1d system
    f.write("LATTICE_VECTORS\n"+str(h.geometry.celldis)+"\n\n")
  if h.dimensionality == 2:  # if  2d system
    f.write("LATTICE_VECTORS\n")
    f.write(str(h.geometry.a1[0])+"   "+str(h.geometry.a1[1])+"\n")
    f.write(str(h.geometry.a2[0])+"   "+str(h.geometry.a2[1])+"\n")

  f.close()


# detect non vanishing elements of a matrix
def nv_el(m):
  """ get the non vanishing elments of a matrix"""
  from scipy.sparse import csc_matrix as csc
  mc = csc(m) # to coo_matrixi
  mc.eliminate_zeros()
  mc = mc.tocoo()
  data = mc.data # get data
  col = mc.col # get column index
  row = mc.row # get row index
  nv = []
  nt=len(data)
  for i in range(nt):
   # save the nonvanishing values
   nv.append([row[i]+1,col[i]+1,data[i].real,data[i].imag])
  return nv


def write_operator(h,name):
  """ Writes an operator in operator.in"""
  f = open("operator.in","w")
  f.write("# "+name+"\n") # write the name of the operator
# edges of a hybrid ribbon (2 last atoms)
  if name=="edges":  
    nind = 2 # two first atoms
    if h.has_spin:
      nind = nind*2 # times spin
    if h.has_eh:
      nind = nind*2 # times electron hole
    f.write("#Number of elements\n") # number of elements
    f.write(str(nind*2)+"\n") # number of elements
    f.write("#i  j real imag\n") # indexes
    n = len(h.intra) # number of elments of the hamiltonian
    for i in range(nind): 
      f.write(str(i+1)+"  "+str(i+1)+"  1.0  0.0\n")  # lower part
      f.write(str(n-i)+"  "+str(n-i)+"  1.0  0.0\n")  # upper part
# center of a hybrid ribbon (4 last atoms)
  if name=="interface":  
    nind = 2 # two first atoms
    dind = 1 # index to which divide the positions
    if h.has_spin:
      nind = nind*2 # times spin
      dind *= 2
    if h.has_eh:
      nind = nind*2 # times electron hole
      dind *= 2
    n = len(h.intra) # number of elments of the hamiltonian
    f.write("#Number of elements\n") # number of elements
    f.write(str(n)+"\n") # number of elements
    f.write("#i  j real imag\n") # indexes
    for i in range(n):
      y = h.geometry.y[i/dind]
      if y < -3.:  
        f.write(str(i+1)+"  "+str(i+1)+"  -1.0  0.0\n")  # middle part
      elif y > 3.:  
        f.write(str(i+1)+"  "+str(i+1)+"  1.0  0.0\n")  # middle part
      else:  
        f.write(str(i+1)+"  "+str(i+1)+"  0.0  0.0\n")  # middle part
# upper edge
  if name=="upper_edge":  
    nind = 2 # two first atoms
    if h.has_spin:
      nind = nind*2 # times spin
    if h.has_eh:
      nind = nind*2 # times electron hole
    f.write("#Number of elements\n") # number of elements
    f.write(str(nind)+"\n") # number of elements
    f.write("#i  j real imag\n") # indexes
    n = len(h.intra) # number of elments of the hamiltonian
    for i in range(nind): 
      f.write(str(n-i)+"  "+str(n-i)+"  1.0  0.0\n")  # upper part
# lower edge
  if name=="lower_edge":  
    nind = 2 # two first atoms
    if h.has_spin:
      nind = nind*2 # times spin
    if h.has_eh:
      nind = nind*2 # times electron hole
    f.write("#Number of elements\n") # number of elements
    f.write(str(nind)+"\n") # number of elements
    f.write("#i  j real imag\n") # indexes
    n = len(h.intra) # number of elments of the hamiltonian
    for i in range(nind): 
      f.write(str(i+1)+"  "+str(i+1)+"  1.0  0.0\n")  # lower part
  if name=="sublattice": 
    if h.has_eh: raise 
    f.write("#Number of elements\n") # number of elements
    if h.has_spin:   # dimension of the operator 
      f.write(str(len(h.geometry.x))+"\n") # number of elements
    else:    
      f.write(str(len(h.geometry.x)*2)+"\n") # number of elements
    f.write("#i  j real imag\n") # indexes
    if h.geometry.has_sublattice:
      sl = h.geometry.sublattice
      for i in range(len(sl)):
        if h.has_spin:
          f.write(str(2*i+1)+"  "+str(2*i+1)+"   "+str(sl[i])+"  0.0\n")
          f.write(str(2*i+2)+"  "+str(2*i+2)+"   "+str(sl[i])+"  0.0\n")
        else:
          f.write(str(i+1)+"  "+str(i+1)+"   "+str(sl[i])+"  0.0\n")  # lower part
    else: raise # error if there is no sublattice
  f.close()







def read_hamiltonian(h,input_file="hamiltonian.in"):
  """ Read a hamiltonian from hamiltonian.in """
  print("reading hamiltonian from",input_file)
  lines = open(input_file,"r").readlines() # read the file
  nl = range(len(lines))
  # get the dimension
  for (l,i) in zip(lines,nl):
    if "DIMENSION_OF_THE_HAMILTONIAN" in l: 
      dd = int(lines[i+1].split()[0]) # dimension of the hamiltonian
      break
  for (l,i) in zip(lines,nl):
    if "DIMENSIONALITY_OF_THE_SYSTEM" in l: 
      dim = int(lines[i+1].split()[0]) # dimensionality of the hamiltonian
      h.dimensionality = dim
      break
  # get the onsite
  def read_name(mm):
    mat = np.matrix([[0.0j for i in range(dd)] for j in range(dd)])
    for (l,i) in zip(lines,nl):
      if mm in l: 
        first = i+1 # first line with onsite element
        break
    for il in range(first,len(lines)):  # loop over lines
      l = lines[il].split()
      if len(l)<4:  # stop if blanck line
        break
      i = int(l[0])-1
      j = int(l[1])-1
      re = float(l[2])
      im = float(l[3])
      mat[i,j] = re + 1j*im # add element
    return mat
   # now read the matrices
  h.intra = read_name("ONSITE_MATRIX")  # read onsite
  if h.dimensionality == 1:
    h.inter = read_name("HOPPING_MATRIX")  # read hopping
  if h.dimensionality == 2:
    h.tx = read_name("HOPPING_MATRIX_1")  # read hopping
    h.ty = read_name("HOPPING_MATRIX_0_1")  # read hopping
    h.txy = read_name("HOPPING_MATRIX_1_1")  # read hopping
    h.txmy = read_name("HOPPING_MATRIX_1_-1")  # read hopping
  h.has_eh = False
  from geometry import read as read_geometry
  try:
    h.geometry = read_geometry(input_file="POSITIONS.OUT")
  except:
    print("geometry hasn't been read")




def write_mean_field_operators(h,hubbard = 1.0,
             output_file = "mean_field_operators.in"):
  """ Writes the mean field operators to the file,
   mandatory input is hamiltonian"""
  if not h.has_spin: # check if it has spin
    print("Hamiltonian doesn't have spin degree")
    raise
  fo = open(output_file,"w")  # open the file
  # part without electron hole sector
  if not h.has_eh: # check that it doesn't have electron-hole sector
    norb = len(h.intra)/2 # number of orbitals (without spin)
    fo.write("Number of matrices, noncollinear calculation, no electron-hole\n")  
    fo.write(str(norb*2)+"\n")  # number of matrices
    for iorb in range(norb): # loop over orbitals
      i = 2*iorb + 1 # index i
      j = 2*iorb + 2 # index j
      # density density term
      fo.write("Matrix A "+str(i)+"\n")  # matrix A
      fo.write("Number of non vanishing elements\n1\n")  # nv elements
      fo.write(str(i)+"  "+str(i)+"  1.d00  0.d00\n") # up density 
      fo.write("Matrix B "+str(i)+"\n")  # matrix A
      fo.write("Number of non vanishing elements\n1\n")  # nv elements
      fo.write(str(j)+"  "+str(j)+"  1.d00  0.d00\n") # down density 
      fo.write("Coupling\n"+str(hubbard)+"  0.d00\n\n")
      # exchange exchange term
      fo.write("Matrix A "+str(j)+"\n")  # matrix A
      fo.write("Number of non vanishing elements\n1\n")  # nv elements
      fo.write(str(i)+"  "+str(j)+"  1.d00  0.d00\n") # up density 
      fo.write("Matrix B "+str(j)+"\n")  # matrix A
      fo.write("Number of non vanishing elements\n1\n")  # nv elements
      fo.write(str(j)+"  "+str(i)+"  1.d00  0.d00\n") # down density 
      fo.write("Coupling\n-"+str(hubbard)+"  0.d00\n\n")

  # with electron hole
  if h.has_eh: # check that it has electron-hole sector
    print("wrong way!!!!")
    raise
    print("=== writting mean field matrices with electron hole sector ===")
    norb = len(h.intra)/4 # number of orbitals (without spin)
    fo.write("Number of matrices, noncollinear calculation, electron-hole\n")  
    fo.write(str(norb*2)+"\n")  # number of matrices
    for iorb in range(norb): # loop over orbitals
      ii = 2*iorb + 1 # index i
      jj = 2*iorb + 3 # index j
      for pap in range(2):  # loop for particle-antiparticle
        i = ii + pap
        j = jj + pap
        ehc = (-1.0)**pap # sign for electron hole
        # density density term
        fo.write("Matrix A "+str(i)+"\n")  # matrix A
        fo.write("Number of non vanishing elements\n1\n")  # nv elements
        fo.write(str(i)+"  "+str(i)+"  1.d00  0.d00\n") # up density 
        fo.write("Matrix B "+str(i)+"\n")  # matrix A
        fo.write("Number of non vanishing elements\n1\n")  # nv elements
        fo.write(str(j)+"  "+str(j)+"  1.d00  0.d00\n") # down density 
        fo.write("Coupling\n"+str(ehc*hubbard)+"  0.d00\n\n")
        # exchange exchange term
        fo.write("Matrix A "+str(j)+"\n")  # matrix A
        fo.write("Number of non vanishing elements\n1\n")  # nv elements
        fo.write(str(i)+"  "+str(j)+"  1.d00  0.d00\n") # up density 
        fo.write("Matrix B "+str(j)+"\n")  # matrix A
        fo.write("Number of non vanishing elements\n1\n")  # nv elements
        fo.write(str(j)+"  "+str(i)+"  1.d00  0.d00\n") # down density 
        fo.write("Coupling\n"+str(-ehc*hubbard)+"  0.d00\n\n")

  fo.close() # close the file




class mean_field():
  """ Class for mean field hamiltonians """
  matrix = None
  def read(self,input_file = "mean_field.in"):
    """ Read the matrix """
    self.matrix = read_mean_field_hamiltonian(input_file=input_file,
                                               traceless = False)
  def write(self,output_file = "mean_field.in"):
    write_mean_field_hamiltonian(self.matrix,output_file=output_file)
    







def read_mean_field(input_file="mean_field.in",traceless = True):
  """Reads the file input_file and returns a traceless matrix"""
  lines = open(input_file,"r").readlines()
  # read the file as a matrix
  if len(lines)==3: return 0. # return 0 if no elements
  i = []
  j = []
  a = []
  for ii in range(3,len(lines)):
    l = lines[ii].split()
    i += [int(l[0])] # index i
    j += [int(l[1])] # index j
    a += [float(l[2])+1j*float(l[3])]   # a_ij
  # create matrix
  n = max([max(i),max(j)]) # dimension of the matrix
  m = np.matrix([[0.0j for ii in range(n)] for jj in range(n)])
  for (ii,jj,aa) in zip(i,j,a):
    m[ii-1,jj-1] = aa
  if traceless:  # if the desired matrix is traceless
    tr = m.trace()[0,0]/n
    for ii in range(n):
      m[ii,ii] += -tr
  return m



def write_mean_field(m,output_file="mean_field.in"):
  """Reads the file input_file and returns a traceless matrix"""
  f = open(output_file,"w") # open file
  mf = coo(m) # sparse matrix
  f.write("# number of non vanishing elements\n"+str(len(mf.col))+"\n")
  f.write("# i   j     real      imaginary\n")
  for (i,j,d) in zip(mf.row,mf.col,mf.data):
    f.write(str(i+1)+"    "+str(j+1)+"   "+str(d.real)+"   "+str(d.imag)+"\n")
  f.close() # close file










def write_magnetization(f,x,y):
  """Writes the magnetization associated with that function"""
  fm = open("MAGNETIZATION.OUT","w")
  for i in range(len(x)):
    m = f(x[i],y[i]) # get magnetization
    fm.write(str(i+1)+"  ")
    fm.write(str(m[0])+"  ")
    fm.write(str(m[1])+"  ")
    fm.write(str(m[2])+"  ")
    fm.write("\n")
  fm.close()


def read_magnetization(input_file = "MAGNETIZATION.OUT"):
  m = np.genfromtxt("MAGNETIZATION.OUT")
  mag = np.array([[m[i,1],m[i,2],m[i,3]] for i in range(len(m))])
  return mag



def spin_stiffness(h,vector=np.array([0.,0.,1.]),
                     qrange = None,atoms = None):
  """ Calculates the spin stiffness by performing a set
  of selfconsistent calculations"""
  from copy import deepcopy
  if qrange == None:
    qrange = np.arange(0.01,2.,0.02)
  feq = open("EnergyVSq.OUT","w")
  feq.write("#    qx       qy       qz      Energy\n")
  os.system("rm -r scfs")
  os.system("mkdir scfs")
  for q in qrange:  # loop over qvectors
    hi = deepcopy(h)  # copy the hamiltonian
    hi.generate_spin_spiral(vector=vector,angle=q,atoms=atoms) 
# generate spin spiral in the hamiltonian
    os.system("mkdir scfs/"+str(q))
    os.system("cp tb90.in scfs/"+str(q)+"/")
    os.chdir("scfs/"+str(q))
    hi.write(output_file="hamiltonian_0.in") # write te hamiltonian
    os.system("tb90.x") # run the calculation
    try:
      energy = np.genfromtxt("ENERGY.OUT").transpose()[0][-1] # last energy
    except:  # only one iteration
      energy = np.genfromtxt("ENERGY.OUT")[0] # last energy
    os.chdir("../..")
    iq = vector*q
    feq.write(str(iq[0])+"       ")
    feq.write(str(iq[1])+"       ")
    feq.write(str(iq[2])+"       ")
    feq.write(str(energy)+"\n") # write in file
  feq.close()




def get_energy():
  try:
      energy = np.genfromtxt("ENERGY.OUT").transpose()[0][-1] # last energy
  except:  # only one iteration
      energy = np.genfromtxt("ENERGY.OUT")[0] # last energy
  return energy




def script_spin_stiffness(h,vector=np.array([0.,0.,1.]),
                     qrange = None,atoms = None):
  """ Writes a script to calculates the spin stiffness by performing a set
  of selfconsistent calculations"""
  from copy import deepcopy
  if qrange == None:
    qrange = np.arange(0.01,2.,0.02)
  # move to a new folder
  os.system("rm -r scfs")
  os.system("mkdir scfs")
  os.chdir("scfs")
  pwd = os.getcwd()
  # write script
  fpy = open("multiple_scf.py","w")
  fq = open("qvectors.dat","w")
  ffol = open("folders.dat","w")
  fpy.write("import os\n\n")
  fq.write("#    qx       qy       qz    name\n")
  ifol = 0
  fpy.write("pwd = os.getcwd()\n")
  for q in qrange:  # loop over qvectors
    name = "folder"+str(ifol)
    hi = deepcopy(h)  # copy the hamiltonian
    hi.generate_spin_spiral(vector=vector,angle=q,atoms=atoms) 
# generate spin spiral in the hamiltonian
    os.system("mkdir "+name)
    os.chdir(name)
    fpy.write("\nos.chdir('"+name+"')\n") # run the calculation
    hi.write(output_file="hamiltonian_0.in") # write te hamiltonian
    fpy.write("os.system('cp ../../tb90.pbs ./')\n") # run the calculation
    fpy.write("os.system('cp ../../tb90.in ./')\n") # run the calculation
    fpy.write("os.system('qsub tb90.pbs')\n") # run the calculation
    fpy.write("os.chdir(pwd)\n") # run the calculation
    iq = vector*q
    fq.write(str(iq[0])+"       ")
    fq.write(str(iq[1])+"       ")
    fq.write(str(iq[2])+"\n")
    ffol.write(name+"\n") # write in file
    ifol += 1  # increase counter
    os.chdir(pwd)
  fq.close()
  ffol.close()
  fpy.close()

def get_spiral_energies():
  """Get energies in a spin spiral calculation"""
  fo = open("E_VS_Q.OUT","w")
  os.chdir("scfs/")
  folders = open("folders.dat").readlines()
  qs = open("qvectors.dat").readlines()
  for (f,q) in zip(folders,qs):
    f = f.split()[0]
    os.chdir(f) # move to that file
    try: 
      energy = get_energy() # get the energy
      q = q.split()[0]
      fo.write(q+"    ") 
      fo.write(str(energy)+"\n") 
    except: pass 
    os.chdir("..") # move to that file
  os.chdir("..") # move to that file
  fo.close()




def replace(name,value,input_file="tb90.in",output_file="tb90.in"):
  """ Modifies a parameter in the tb90.in file"""
  from copy import deepcopy
  lines = open(input_file,"r").readlines()
  lout = []
  for l in lines:
    if " "+name+" " in l: # if name in line
      lo = "   "+name+" = "+value+" \n" # replace line
    else:
      lo = l
    lout.append(lo) # add to the list
  fo = open(output_file,"w") # open output file
  for l in lout:
    fo.write(l) 
  fo.close()









