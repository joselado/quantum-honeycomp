#!/usr/bin/python3

import sys
import os

# Root path
qhroot = os.environ["QHROOT"] # root path

# Add path of the wrapper


sys.path.append(qhroot+"/pysrc/") # python libraries

import qtwrap # import the library with simple wrappaers to qt4
get = qtwrap.get  # get the value of a certain variable
is_checked = qtwrap.is_checked  # get the value of a certain variable
getbox = qtwrap.getbox  # get the value of a certain variable
window = qtwrap.main() # this is the main interface


from qh_interface import * # import all the libraries needed


from interfacetk import interfacetk
modify_geometry = lambda x: interfacetk.modify_geometry(x,qtwrap)


def get_geometry(modify=True):
  """ Create a 2d honeycomb lattice"""
  n = int(get("cell_size")) # size of the unit cell
  g = specialgeometry.twisted_bilayer(n)
#  g = geometry.honeycomb_lattice()
#  g = g.supercell(n)
  if modify: g = modify_geometry(g) # remove atoms if necessary
  return g



def initialize():
  """ Initialize the calculation"""
  g = get_geometry() # get the geometry
  twisted_matrix = specialhopping.twisted_matrix
  has_spin = False
  h = g.get_hamiltonian(is_sparse=True,has_spin=has_spin,is_multicell=True,
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
  if is_checked("set_half_filling"): h.set_filling(nk=2)
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
  comp = computing() # create the computing window
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
  comp.kill()
  execute_script("qh-bands2d ")
  
  




  

def show_dos(self):
  comp = computing() # create the computing window
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
  comp.kill()
  return




def show_fermi_surface(silent=False):
  h = pickup_hamiltonian() # get hamiltonian
  ndos = int(get("ne_dos"))
  if h.dimensionality==2:
    spectrum.fermi_surface(h,e=get("energy_fs"),nk=int(get("nk_fs")),
            nsuper = 2,reciprocal=True,delta=get("delta_fs"),
            mode="lowest",num_waves=10)
  else: raise
  if not silent:
      execute_script("qh-fermi-surface FERMI_MAP.OUT") # show the result


def show_structure():
  g = get_geometry() # get the geometry
  nsuper = int(get("nsuper_struct")) 
  g = g.supercell(nsuper) # build a supercell
  g.write()
  execute_script("qh-potential POSITIONS.OUT ")



def pickup_hamiltonian():
    return initialize()







def show_ldos():
  h = pickup_hamiltonian()  # get the hamiltonian
#  if h.intra.shape[0]<2000: h.turn_dense()
  e = get("energy_ldos_single")
  delta = get("delta_ldos_single")
  nk = get("nk_ldos_single")
  nk = int(round(np.sqrt(nk)))
  nsuper = int(get("nsuper_ldos_single"))
  ldos.ldos2d(h,e=e,delta=delta,nk=nk,mode="arpack",nrep=nsuper)
  execute_script("qh-fast-ldos LDOS.OUT  ")
#  execute_script("qh-multildos ")
#  execute_script("tb90-calculate-ldos "+str(get("stm_bias")))
#  execute_script("qh-interpolate LDOS.OUT ")
#  execute_script("tb90-cmap LDOS.OUT-interpolated ")



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
  g = get_geometry() # get the geometry
  nsuper = int(get("nsuper_struct"))
  g = g.supercell(nsuper)
  g.write()
#  execute_script("qh-structure3d POSITIONS.OUT")
#  execute_script("qh-magnetism nomag nobonds POSITIONS.OUT")
  execute_script("qh-structure-tbg ")




def show_interactive_ldos():
  comp = computing() # create the computing window
  h = pickup_hamiltonian()  # get the hamiltonian
  h.turn_dense()
  ewin = get("window_ldos")
  nrep = int(get("nsuper_ldos"))
  nk = int(np.sqrt(get("nk_ldos")))
  ne = int(get("ne_ldos"))
  delta = get("delta_ldos")
  ldos.multi_ldos(h,es=np.linspace(-ewin,ewin,ne),nk=nk,delta=delta,
          nrep=nrep)
  comp.kill()
  execute_script("qh-multildos ")


def select_atoms_removal(self):
  g = get_geometry(modify=False) # get the unmodified geometry
  g.write() # write geometry
  execute_script("qh-remove-atoms-geometry-3d") # remove the file


save_results = lambda x: save_outputs(inipath,tmppath) # function to save

# create signals
signals = dict()
#signals["on_window_destroy"] = gtk.main_quit  # close the window
signals["show_bands"] = show_bands  # show bandstructure
signals["show_dos"] = show_dos  # show DOS
#signals["show_chern"] = show_chern  # show Chern number 
#signals["show_berry_1d"] = show_berry_1d  # show Berry curvature
#signals["show_berry_2d"] = show_berry_2d  # show Berry curvature
signals["show_ldos_single"] = show_ldos  # show Berry curvature
#signals["show_z2_invariant"] = show_z2_invariant  # show Berry curvature
#signals["show_magnetism"] = show_magnetism  # show magnetism
signals["show_structure"] = show_structure  # show magnetism
signals["show_structure_3d"] = show_structure_3d
signals["show_dosbands"] = show_dosbands  # show magnetism
signals["show_fermi_surface"] = show_fermi_surface  # show magnetism
signals["show_interactive_ldos"] = show_interactive_ldos  # show magnetism
signals["select_atoms_removal"] = select_atoms_removal
#signals["show_2dband"] = show_2dband  # show magnetism
#signals["show_kdos"] = show_kdos  # show kdos
#signals["save_results"] = save_results  # save the results


window.connect_clicks(signals)
folder = create_folder()
tmppath = os.getcwd() # get the initial directory
window.run()

