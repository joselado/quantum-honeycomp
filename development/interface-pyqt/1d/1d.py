#!/usr/bin/python

from __future__ import print_function

import sys
import os

qhroot = os.environ["QHROOT"] # root path
sys.path.append(qhroot+"/pysrc/") # python libraries


from interfacetk import qtwrap # import the library with simple wrappaers to qt4
get = qtwrap.get  # get the value of a certain variable
getbox = qtwrap.getbox  # get the value of a certain variable
window = qtwrap.main() # this is the main interface



from interfacetk.qh_interface import * # import all the libraries needed
from interfacetk import common # common routines for all the geometries
common.initialize(qtwrap) # do several common initializations






def get_geometry(modify=True):
  """ Create a 0d island"""
  lattice_name = getbox("lattice") # get the option
  n = int(get("width")) # thickness of the system
#  lattice_name = builder.get_object("lattice").get_active_text()
  if lattice_name=="Chain":
    g = geometry.chain()
  if lattice_name=="Honeycomb":
    g = geometry.honeycomb_lattice()
  elif lattice_name=="Square":
    g = geometry.square_lattice()
  elif lattice_name=="Kagome":
    g = geometry.kagome_lattice()
  elif lattice_name=="Lieb":
    g = geometry.lieb_lattice()
  elif lattice_name=="Triangular":
    g = geometry.triangular_lattice_tripartite()
  elif lattice_name=="Honeycomb zigzag":
    g = geometry.honeycomb_zigzag_ribbon(n)
  elif lattice_name=="Honeycomb armchair":
    g = geometry.honeycomb_armchair_ribbon(n)
  if g.dimensionality==2: # original is a 2d geometry
    g = ribbon.bulk2ribbon(g,n=n)
  nsuper = int(get("nsuper"))
  g = g.supercell(nsuper)
  if modify: g = modify_geometry(g) # modify the geometry
  return g




def select_atoms_removal(self):
  g = get_geometry(modify=False) # get the unmodified geometry
  g.write() # write geometry
  execute_script("qh-remove-atoms-geometry") # remove the file


def modify_geometry(g):
  """Modify the geometry according to the interface"""
  if qtwrap.is_checked("remove_selected"): # remove some atoms
      try:
        inds = np.array(np.genfromtxt("REMOVE_ATOMS.INFO",dtype=np.int))
        if inds.shape==(): inds = [inds]
      except: inds = [] # Nothing
      print(inds)
      g = sculpt.remove(g,inds) # remove those atoms
  if qtwrap.is_checked("remove_single_bonded"): # remove single bonds
      g = sculpt.remove_unibonded(g,iterative=True)
  return g # return geometry








def initialize():
  """ Initialize the calculation"""
  g = get_geometry() # get the geometry
  h = g.get_hamiltonian(has_spin=True)
  h.add_peierls(get("peierls")) # magnetic field
  h.add_zeeman([get("Bx"),get("By"),get("Bz")]) # Zeeman fields
  h.add_sublattice_imbalance(get("mAB"))  # sublattice imbalance
  if abs(get("rashba")) > 0.0: h.add_rashba(get("rashba"))  # Rashba field
  h.add_antiferromagnetism(get("mAF"))  # AF order
  h.shift_fermi(get("fermi")) # shift fermi energy
  if abs(get("kanemele"))>0.0:  h.add_kane_mele(get("kanemele")) # intrinsic SOC
  if abs(get("haldane"))>0.0:  h.add_haldane(get("haldane")) # intrinsic SOC
  if abs(get("antihaldane"))>0.0:  h.add_antihaldane(get("antihaldane")) 
  if abs(get("antikanemele"))>0.0:  h.add_anti_kane_mele(get("antikanemele")) 
  if abs(get("swave"))>0.0:  h.add_swave(get("swave")) 
#  h.add_peierls(get("peierls")) # shift fermi energy

  return h


def show_bands(self=0):
  h = pickup_hamiltonian() # get hamiltonian
  common.get_bands(h,qtwrap) # wrapper


def show_dosbands(self=0):
  comp = computing() # create the computing window
  h = pickup_hamiltonian() # get hamiltonian
  kdos.kdos_bands(h,scale=get("scale_kbands"),ewindow=get("window_kbands"),
                   ne=int(get("ne_kbands")),delta=get("delta_kbands"),
                   ntries=int(get("nv_kbands")))
  comp.kill()
  execute_script("qh-dosbands1d  KDOS_BANDS.OUT ")




def show_interactive_ldos():
  comp = computing() # create the computing window
  h = pickup_hamiltonian()  # get the hamiltonian
  ewin = get("window_ldos")
  nrep = int(get("nsuper_ldos"))
  nk = int(get("nk_ldos"))
  ne = int(get("ne_ldos"))
  delta = get("delta_ldos")
  opname = getbox("operator_ldos") # get the operator
  if opname=="None": op = None
  elif opname=="Current": op = operators.get_current(h)
  elif opname=="Spin current": op = operators.get_spin_current(h)
  else: raise
  ldos.multi_ldos(h,es=np.linspace(-ewin,ewin,ne),nk=nk,
          delta=delta,nrep=nrep,op=op)
  comp.kill()
  if op is None: execute_script("qh-multildos ")
  else: execute_script("qh-multildos-signed ")






def show_dos(self):
  comp = computing() # create the computing window
  h = pickup_hamiltonian() # get hamiltonian
#  mode = getbox("mode_dos") # mode for the DOS
  if h.dimensionality==0:
    dos.dos0d(h,es=np.linspace(-3.1,3.1,500),delta=get("DOS_smearing"))
  elif h.dimensionality==1:
#    dos.dos1d(h,ndos=400,delta=get("DOS_smearing"))
    dos.dos1d(h,ndos=400)
  elif h.dimensionality==2:
    dos.dos2d(h,ndos=500,delta=get("DOS_smearing"))
  else: raise
  comp.kill()
  execute_script("qh-dos  DOS.OUT")



def pickup_hamiltonian():
  if qtwrap.is_checked("do_scf"):
    return hamiltonians.load() # load the Hamiltonian
  else: # generate from scratch
    return initialize()



def show_magnetism():
  """Show the magnetism of the system"""
  h = pickup_hamiltonian() # get the Hamiltonian
  common.show_exchange(h,qtwrap)



  


def show_structure(self):
  """Show the lattice of the system"""
  g = get_geometry() # get the geometry
  nsuper = int(get("nsuper_struct"))
  g = g.supercell(nsuper)
  g.write()
#  execute_script("qh-light-structure POSITIONS.OUT")
  execute_script("qh-structure-bond POSITIONS.OUT")
#  execute_script("qh-structure  ")



def show_structure_3d(self):
  """Show the lattice of the system"""
  g = get_geometry() # get the geometry
  nsuper = int(get("nsuper_struct"))
  g = g.supercell(nsuper)
  g.write()
  execute_script("qh-structure3d POSITIONS.OUT")





def solve_scf():
  """Perform a selfconsistent calculation"""
  comp = computing() # create the computing window
  scfin = getbox("scf_initialization")
  h = initialize() # initialize the Hamiltonian
  mf = scftypes.guess(h,mode=scfin)
  nk = int(get("nk_scf"))
  U = get("hubbard")
  filling = get("filling_scf")
  filling = filling%1.
  scf = scftypes.selfconsistency(h,nkp=nk,filling=filling,g=U,
                mf=mf,mode="U",smearing=get("smearing_scf"),
                mix = get("mix_scf"))
  scf.hamiltonian.save() # save in a file
  comp.kill()






save_results = lambda x: save_outputs(inipath,tmppath) # function to save


# create signals
signals = dict()
#signals["initialize"] = initialize  # initialize and run
signals["show_bands"] = show_bands  # show bandstructure
signals["show_structure"] = show_structure  # show bandstructure
signals["show_dos"] = show_dos  # show DOS
signals["show_dosbands"] = show_dosbands  # show DOS
signals["show_interactive_ldos"] = show_interactive_ldos  # show DOS
signals["show_structure_3d"] = show_structure_3d
signals["show_magnetism"] = show_magnetism 
signals["solve_scf"] = solve_scf
signals["select_atoms_removal"] = select_atoms_removal




window.connect_clicks(signals)
folder = create_folder()
tmppath = os.getcwd() # get the initial directory
window.run()

