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



def custom_scf(h):
  """Modifies the properties of tb90.in input"""
  # begin SCF options
  scfin = builder.get_object("scf_initialization").get_active_text()
  import random
  rand = random.random
  def rv():
    return .5 - np.random.random(3) # 3 component random vector
  def rvxy():
    v = np.random.random(3) -.5 # 3 component random vector
    v[0] = 0.0
    return v
  ###################################
  if scfin == "Reconverge": # no new mean field
    mf = np.load("MEAN_FIELD.npy") # load the mean field from file
  else: # if create new mean field matrix
    mf = interactions.initialization(h.geometry,scfin=scfin)  # create the mean field mat
  U = get("hubbard") # value of Hubbard
  mix = 0.5
  if U > 1.5: mix = 0.1 # slow mixing
  # perform the selfconsistent calculation
  scf = interactions.hubbardscf(h,U=U,nkp=1,mf=mf,silent=False,mix=mix,
                                  maxerror=1e-5,filling=get("filling"))
  scf.hamiltonian.write() # write in a file
  np.save("MEAN_FIELD.npy",np.array(scf.mean_field)) # save in a file




def get_geometry():
  """ Create a 0d island"""
  lattice_name = getbox("lattice") # get the option
#  lattice_name = builder.get_object("lattice").get_active_text()
  if lattice_name=="Honeycomb":
    geometry_builder = geometry.honeycomb_lattice
  elif lattice_name=="Square":
    geometry_builder = geometry.square_lattice
  elif lattice_name=="Kagome":
    geometry_builder = geometry.kagome_lattice
  elif lattice_name=="Lieb":
    geometry_builder = geometry.lieb_lattice
  elif lattice_name=="Triangular":
    geometry_builder = geometry.triangular_lattice
  g = geometry_builder() # call the geometry
  nsuper = int(get("nsuper"))
  g = g.supercell(nsuper)
  return g











def initialize():
  """ Initialize the calculation"""
  g = get_geometry() # get the geometry
  h = g.get_hamiltonian(has_spin=True)
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
  h = pickup_hamiltonian() # get hamiltonian
  opname = getbox("bands_color")
  if opname=="None": op = None # no operators
  elif opname=="Sx": op = h.get_operator("sx") # off plane case
  elif opname=="Sy": op = h.get_operator("sy")# off plane case
  elif opname=="Sz": op = h.get_operator("sz")# off plane case
  elif opname=="Valley": op = h.get_operator("valley")
  else: op =None
  h.get_bands(operator=op)
  execute_script("qh-bands2d  ")



def show_dosbands(self=0):
  h = pickup_hamiltonian() # get hamiltonian
  kdos.kdos_bands(h,scale=get("scale_kbands"),ewindow=get("window_kbands"),
                   ne=int(get("ne_kbands")),delta=get("delta_kbands"),
                   ntries=int(get("nv_kbands")))
  execute_script("qh-dosbands  KDOS_BANDS.OUT ")



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
  execute_script("tb90-dos  ")


def pickup_hamiltonian():
  return initialize()
  if builder.get_object("activate_scf").get_active():
    return read_hamiltonian()
  else: # generate from scratch
    return initialize()










def show_stm(self):
  h = pickup_hamiltonian() # get hamiltonian
#  ldos.multi_ldos()
  ewin = abs(get("window_ldos")) # energy window
  ne = int(get("num_ldos")) # number of LDOS
  delta = ewin/ne # delta
  ldos.multi_ldos(h,es=np.linspace(-ewin,ewin,ne),nk=1,delta=delta)
  execute_script("qh-multildos ")
#  hamiltonians.ldos(h,e=get("stm_bias"),delta=get("DOS_smearing")) # calculate the stm spectra
#  print("Using semaring",get("DOS_smearing"))
#  execute_script("qh-ldos  LDOS.OUT")
  return


def show_berry2d():
  h = pickup_hamiltonian() # get hamiltonian
  nk = int(get("nk_topology"))
  topology.berry_map(h,nk=nk)
  execute_script("qh-berry2d BERRY_MAP.OUT")

  

def show_magnetism(self):
  h = pickup_hamiltonian() # get hamiltonian
  h.get_magnetization() # get the magnetization
  execute_script("tb90-magnetism  ")
#  execute_script("qh-magnetism  ")


def show_structure(self):
  """Show the lattice of the system"""
  g = get_geometry() # get the geometry
  nsuper = int(get("nsuper_struct"))
  g = g.supercell(nsuper)
  g.write()
#  execute_script("qh-light-structure POSITIONS.OUT")
  execute_script("qh-structure-bond POSITIONS.OUT")
#  execute_script("qh-structure  ")



def show_kdos(self):
  h = pickup_hamiltonian()  # get the hamiltonian
  ew = get("ewindow_kdos")
  new = int(get("mesh_kdos")) # scale as kpoints
  energies = np.linspace(-ew,ew,new) # number of ene
  klist = np.linspace(0.,1.,new)
  kdos.write_surface_2d(h,energies=energies,delta=ew/new,klist=klist)
  execute_script("qh-kdos-both KDOS.OUT  ")



def show_berry1d(self):
  h = pickup_hamiltonian()  # get the hamiltonian
  ks = klist.default(h.geometry,nk=int(get("nk_topology")))  # write klist
  topology.write_berry(h,ks)
  execute_script("qh-berry1d  label  ")


def show_z2(self):
  h = pickup_hamiltonian()  # get the hamiltonian
  nk = get("nk_topology")
  topology.z2_vanderbilt(h,nk=nk,nt=nk/2) # calculate z2 invariant
  execute_script("qh-wannier-center  ") # plot the result





save_results = lambda x: save_outputs(inipath,tmppath) # function to save


# create signals
signals = dict()
signals["initialize"] = initialize  # initialize and run
signals["show_bands"] = show_bands  # show bandstructure
signals["show_structure"] = show_structure  # show bandstructure
signals["show_dos"] = show_dos  # show DOS
signals["show_berry2d"] = show_berry2d  # show DOS
signals["show_berry1d"] = show_berry1d  # show DOS
signals["show_kdos"] = show_kdos  # show DOS
signals["show_dosbands"] = show_dosbands  # show DOS
signals["show_z2"] = show_z2  # show DOS
#signals["show_stm"] = show_stm  # show STM
#signals["show_magnetism"] = show_magnetism  # show magnetism
#signals["show_lattice"] = show_lattice  # show magnetism
#signals["save_results"] = save_results  # save the results





#from qh_interface import create_folder # import all the libraries needed

window.connect_clicks(signals)
folder = create_folder()
tmppath = os.getcwd() # get the initial directory
window.run()

