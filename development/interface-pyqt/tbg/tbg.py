#!/usr/bin/python3

import sys
import os

# Root path
qhroot = os.environ["QHROOT"] # root path

# Add path of the wrapper
sys.path.append(qhroot+"/interface-pyqt/qtwrap")


sys.path.append(qhroot+"/pysrc/") # python libraries

import qtwrap # import the library with simple wrappaers to qt4
get = qtwrap.get  # get the value of a certain variable
getbox = qtwrap.getbox  # get the value of a certain variable
window = qtwrap.main() # this is the main interface


from qh_interface import * # import all the libraries needed




def get_geometry2d():
  """ Create a 2d honeycomb lattice"""
  n = int(get("cell_size")) # size of the unit cell
  import specialgeometry
  g = specialgeometry.twisted_bilayer(n)
#  g = geometry.honeycomb_lattice()
#  g = g.supercell(n)
  return g



def initialize():
  """ Initialize the calculation"""
  g = get_geometry2d() # get the geometry
  from specialhopping import twisted,twisted_matrix
  has_spin = False
  h = g.get_hamiltonian(is_sparse=True,has_spin=has_spin,is_multicell=False,
     mgenerator=twisted_matrix(ti=get("tinter"),lambi=7.0))
#  return h
#  h.turn_dense()
  h.add_sublattice_imbalance(get("mAB"))  # sublattice imbalance
  efield = get("interlayer_bias")
  def bias(r):
    if r[2]<0.0: return efield
    else: return -efield
  h.shift_fermi(bias)
  if h.has_spin:
    h.add_zeeman([get("Bx"),get("By"),get("Bz")]) # Zeeman fields
    h.add_rashba(get("rashba"))  # Rashba field
    h.add_antiferromagnetism(get("mAF"))  # AF order
    h.add_kane_mele(get("kanemele")) # intrinsic SOC
  h.shift_fermi(get("fermi")) # shift fermi energy
  if False:
    h.add_swave(get("swave"))
  if False:
    if h.has_eh:
        print("SCF not implemented with Nambu")
        raise
  if False:
    custom_scf(h) # create the tb90.in
  else:
    h.write("hamiltonian.in")
  klist.default(g,nk=int(get("nkpoints")))  # write klist
#  klist.tr_path(nk=int(get("nkpoints")))  # write klist
  return h



def show_bands(self):
  h = pickup_hamiltonian()  # get the hamiltonian
  if h.intra.shape[0]<4000:
    num_bands = int(get("nbands"))
    if num_bands> h.intra.shape[0]-3: 
      h.turn_dense()
      num_bands = None
    if num_bands<0: 
      h.turn_dense()
      num_bands = None
  else: num_bands = max(20,int(get("nbands")))
  opname = getbox("bands_operator")
  kpath = klist.default(h.geometry,nk=int(get("nk_bands")))  # write klist
  if opname=="None": op = None # no operators
  elif opname=="Valley": op = h.get_operator("valley_upper") # no operators
  h.get_bands(kpath=kpath,num_bands=num_bands,operator=op) 
  execute_script("qh-bands2d ")
  




  

def show_dos(self):
  h = pickup_hamiltonian()  # get the hamiltonian
  nk = int(round(np.sqrt(get("nk_dos"))))
  ndos = int(get("nume_dos"))
  npol = int(get("numpol_dos"))
  ndos = npol*10
  delta = get("delta_dos")
  scale = 10.0 # scale for KPM
  ewindow = get("ewindow_dos")/scale
  dos.dos2d(h,use_kpm=True,nk=nk,ntries=1,delta=delta,random=True,
              ndos=ndos,kpm_window=ewindow,scale=scale)
  execute_script("tb90-dos  ")
  return



def show_magnetism(self):
  h = pickup_hamiltonian() # get hamiltonian
  h.get_magnetization(nkp=int(get("nkpoints"))) # get the magnetization
  execute_script("qh-magnetism 2  ")
#  execute_script("qh-structure 2  ")





def show_chern(self):
  h = pickup_hamiltonian()  # get the hamiltonian
  nk = get("nkpoints")
  nk = int(round(np.sqrt(nk)))
  topology.chern(h,nk=nk) # calculate chern number
  execute_script("tb90-chern") 
#  execute_script("qh-plotchern") 



def show_berry_2d(self):
  h = pickup_hamiltonian()  # get the hamiltonian
  nk = get("nkpoints")
  nk = int(round(np.sqrt(nk)))
  topology.berry_map(h,nk=nk) 
  execute_script("qh-berry2d BERRY_MAP.OUT") 

 

def show_berry_1d(self):
  h = pickup_hamiltonian()  # get the hamiltonian
  ks = klist.default(h.geometry,nk=int(get("nkpoints")))  # write klist
  topology.write_berry(h,ks) 
  execute_script("tb90-berry1d  label  ")


def show_structure():
  g = get_geometry2d() # get the geometry
  nsuper = int(get("nsuper_struct")) 
  g = g.supercell(nsuper) # build a supercell
  g.write()
  execute_script("qh-potential POSITIONS.OUT ")



def pickup_hamiltonian():
  if False:
    return read_hamiltonian()
  else: # generate from scratch
    return initialize()



def read_hamiltonian():
  g = get_geometry2d() # get the geometry
  h = g.get_hamiltonian() # get the hamiltonian
  h.read("hamiltonian.in") # read hamiltonian
  h.has_eh = builder.get_object("has_eh").get_active()
  h.has_spin = builder.get_object("has_spin").get_active()
  return h









def show_ldos():
  h = pickup_hamiltonian()  # get the hamiltonian
#  if h.intra.shape[0]<2000: h.turn_dense()
  e = get("energy_ldos")
  delta = get("delta_ldos")
  nk = get("nk_ldos")
  nk = int(round(np.sqrt(nk)))
  nsuper = int(get("nsuper_ldos"))
  ldos.ldos2d(h,e=e,delta=delta,nk=nk,mode="arpack",nrep=nsuper)
  execute_script("qh-fast-ldos LDOS.OUT  ")
#  execute_script("qh-multildos ")
#  execute_script("tb90-calculate-ldos "+str(get("stm_bias")))
#  execute_script("qh-interpolate LDOS.OUT ")
#  execute_script("tb90-cmap LDOS.OUT-interpolated ")


def show_fermi_surface(self):
  h = pickup_hamiltonian()  # get the hamiltonian
#  spectrum.fermi_surface(h,reciprocal=True,nk=get("nkpoints")/2) # calculate Fs
  spectrum.boolean_fermi_surface(h,reciprocal=True,nk=get("nkpoints")/2) # calculate Fs
#  execute_script("tb90-cmap FERMI_MAP.OUT  ") # plot the result
  execute_script("tb90-cmap BOOL_FERMI_MAP.OUT  ") # plot the result

def show_z2_invariant(self):
  h = pickup_hamiltonian()  # get the hamiltonian
  nk = get("nkpoints")/4
  topology.z2_vanderbilt(h,nk=nk,nt=nk/2) # calculate z2 invariant
  execute_script("qh-wannier-center  ") # plot the result




def show_kdos(self):
  h = pickup_hamiltonian()  # get the hamiltonian
  ew = get("e_kdos")
  new = int(get("nkpoints")/10) # scale as kpoints
  energies = np.linspace(-ew,ew,new) # number of ene
  klist = np.linspace(0.,1.,new)
  kdos.write_surface_2d(h,energies=energies,delta=ew/new,klist=klist)
  execute_script("qh-kdos KDOS.OUT  ")


def show_dosbands(self=0):
  h = pickup_hamiltonian() # get hamiltonian
  kpath = klist.default(h.geometry,nk=int(get("nk_kbands")))
  kdos.kdos_bands(h,scale=get("scale_kbands"),ewindow=get("window_kbands"),
                   ne=int(get("ne_kbands")),delta=get("delta_kbands"),
                   ntries=int(get("nv_kbands")),kpath=kpath)
  execute_script("qh-dosbands  KDOS_BANDS.OUT ")




def show_2dband(self):
  h = pickup_hamiltonian()  # get the hamiltonian
  nk = get("nkpoints")/4
  ns = get_text("num_2dband") # get indexes of the bands
  if "," in ns: ns = [int(n) for n in ns.split(",")] # get the different numbers
  else: ns=[int(ns)] # single one
  if 0 in ns: # in case all eigenvalues wanted
    ns = [i+1 for i in range(h.intra.shape[0]//2)]
    ns += [-i for i in ns]
  spectrum.get_bands(h,nindex=ns,nk=nk,reciprocal=True)
  string = ""
  for n in ns: string += "BANDS2D__"+str(n)+".OUT "
  execute_script("qh-plot3d "+string +"  ")










def show_structure_3d(self):
  """Show the lattice of the system"""
  g = get_geometry2d() # get the geometry
  nsuper = int(get("nsuper_struct"))
  g = g.supercell(nsuper)
  g.write()
#  execute_script("qh-structure3d POSITIONS.OUT")
#  execute_script("qh-magnetism nomag nobonds POSITIONS.OUT")
  execute_script("qh-pick ")








save_results = lambda x: save_outputs(inipath,tmppath) # function to save

# create signals
signals = dict()
#signals["on_window_destroy"] = gtk.main_quit  # close the window
signals["show_bands"] = show_bands  # show bandstructure
signals["show_dos"] = show_dos  # show DOS
#signals["show_chern"] = show_chern  # show Chern number 
#signals["show_berry_1d"] = show_berry_1d  # show Berry curvature
#signals["show_berry_2d"] = show_berry_2d  # show Berry curvature
signals["show_ldos"] = show_ldos  # show Berry curvature
#signals["show_fermi_surface"] = show_fermi_surface  # show Berry curvature
#signals["show_z2_invariant"] = show_z2_invariant  # show Berry curvature
#signals["show_magnetism"] = show_magnetism  # show magnetism
signals["show_structure"] = show_structure  # show magnetism
signals["show_structure_3d"] = show_structure_3d
signals["show_dosbands"] = show_dosbands  # show magnetism
#signals["show_2dband"] = show_2dband  # show magnetism
#signals["show_kdos"] = show_kdos  # show kdos
#signals["save_results"] = save_results  # save the results


window.connect_clicks(signals)
folder = create_folder()
tmppath = os.getcwd() # get the initial directory
window.run()

