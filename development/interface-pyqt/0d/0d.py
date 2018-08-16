#!/usr/bin/python

from __future__ import print_function

import sys
import os

qhroot = os.environ["QHROOT"] # root path
sys.path.append(qhroot+"/interface-pyqt/qtwrap")
sys.path.append(qhroot+"/pysrc/") # python libraries


import qtwrap # import the library with simple wrappaers to qt4
get = qtwrap.get  # get the value of a certain variable
getbox = qtwrap.getbox  # get the value of a certain variable
window = qtwrap.main() # this is the main interface



from qh_interface import * # import all the libraries needed




def get_geometry():
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
  import islands
  rot = get("rotation")*np.pi/180.
  g = islands.get_geometry(n=n,nedges=int(get("nsides")),rot=rot,geo=g)
  return g











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
  if abs(get("swave"))>0.0:  h.add_swave(get("swave")) 
#  h.add_peierls(get("peierls")) # shift fermi energy

  return h


def show_bands(self=0):
  comp = computing() # create the computing window
  h = pickup_hamiltonian() # get hamiltonian
  opname = getbox("bands_color")
  if opname=="None": op = None # no operators
  elif opname=="Sx": op = h.get_operator("sx") # off plane case
  elif opname=="Sy": op = h.get_operator("sy")# off plane case
  elif opname=="Sz": op = h.get_operator("sz")# off plane case
  elif opname=="Valley": op = h.get_operator("valley")
  elif opname=="Vy-position": op = h.get_operator("yposition")
  else: op =None
  h.get_bands(operator=op)
  comp.kill()
  execute_script("qh-bands1d  ")



def show_dosbands(self=0):
  h = pickup_hamiltonian() # get hamiltonian
  kdos.kdos_bands(h,scale=get("scale_kbands"),ewindow=get("window_kbands"),
                   ne=int(get("ne_kbands")),delta=get("delta_kbands"),
                   ntries=int(get("nv_kbands")))
  execute_script("qh-dosbands1d  KDOS_BANDS.OUT ")




def show_interactive_ldos():
  comp = computing() # create the computing window
  h = pickup_hamiltonian()  # get the hamiltonian
  ewin = get("window_ldos")
  nrep = 1
  nk = 1
  ne = int(get("ne_ldos"))
  delta = get("delta_ldos")
  ldos.multi_ldos(h,es=np.linspace(-ewin,ewin,ne),nk=nk,delta=delta,nrep=nrep)
  comp.kill()
  execute_script("qh-multildos ")






def show_dos(self):
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
  execute_script("qh-dos  ")


def pickup_hamiltonian():
  return initialize()
  if builder.get_object("activate_scf").get_active():
    return read_hamiltonian()
  else: # generate from scratch
    return initialize()




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
  nk = 1
  U = get("hubbard")
  filling = get("filling_scf")
  filling = filling%1.
  scf = scftypes.selfconsistency(h,nkp=nk,filling=filling,g=U,
                mf=mf,mode="U",smearing=get("smearing_scf"),
                mix = get("mix_scf"))
  scf.hamiltonian.save() # save in a file
  comp.kill()


def pickup_hamiltonian():
  if qtwrap.is_checked("do_scf"):
    return hamiltonians.load() # load the Hamiltonian
  else: # generate from scratch
    return initialize()





def show_magnetism():
  """Show the magnetism of the system"""
  h = pickup_hamiltonian() # get the Hamiltonian
  h.write_magnetization() # write the magnetism
  execute_script("qh-moments",mayavi=True)








save_results = lambda x: save_outputs(inipath,tmppath) # function to save


# create signals
signals = dict()
#signals["initialize"] = initialize  # initialize and run
signals["solve_scf"] = solve_scf
signals["show_bands"] = show_bands  # show bandstructure
signals["show_structure"] = show_structure  # show bandstructure
#signals["show_dos"] = show_dos  # show DOS
signals["show_structure_3d"] = show_structure_3d
signals["show_interactive_ldos"] = show_interactive_ldos  # show DOS
signals["show_magnetism"] = show_magnetism





#from qh_interface import create_folder # import all the libraries needed

window.connect_clicks(signals)
folder = create_folder()
tmppath = os.getcwd() # get the initial directory
window.run()

