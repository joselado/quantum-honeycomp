#!/usr/bin/python

from __future__ import print_function

import sys
import os

sys.path.append(os.getcwd()+"/../../pysrc") # append library
os.environ["QHROOT"] = os.getcwd()+"/../.." # fine
#print(os.environ["QHROOT"])

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
  dim = getbox("geometry_dimensionality") # get the option
  if dim=="Island (0D)": raise
  elif dim=="Ribbon (1D)":
    g = geometry_builder() # call the geometry
    g = g.supercell((1,10))
    g.dimensionality = 1
  elif dim=="Sheet (2D)": 
    g = geometry_builder() # call the geometry
  return g











def initialize(self=0):
  """ Initialize the calculation"""
  g = get_geometry() # get the geometry
  h = g.get_hamiltonian(has_spin=True)
  h.add_zeeman([get("Bx"),get("By"),get("Bz")]) # Zeeman fields
  print("Bz",get("Bz"))
  h.add_sublattice_imbalance(get("mAB"))  # sublattice imbalance
  if abs(get("rashba")) > 0.0: h.add_rashba(get("rashba"))  # Rashba field
  h.add_antiferromagnetism(get("mAF"))  # AF order
  if abs(get("kanemele"))>0.0:  h.add_kane_mele(get("kanemele")) # intrinsic SOC
  h.shift_fermi(get("fermi")) # shift fermi energy
  h.add_peierls(get("peierls")) # shift fermi energy
  h.write("hamiltonian.in") # write the Hamiltonian in file

#  if builder.get_object("activate_scf").get_active():
#    custom_scf(h) # perform selfconsistent calculation
#  else: pass
  return h


def show_bands(self=0):
  h = pickup_hamiltonian() # get hamiltonian
#  opname = builder.get_object("bands_color").get_active_text()
  opname = getbox("bands_color")
  if opname=="None": op = None # no operators
  elif opname=="Sx": op = h.get_operator("sx") # off plane case
  elif opname=="Sy": op = h.get_operator("sy")# off plane case
  elif opname=="Sz": op = h.get_operator("sz")# off plane case
  elif opname=="Valley": op = h.get_operator("valley")# off plane case
  else: op =None
  h.get_bands(operator=op)

  execute_script("qh-bands0d  ")


  

def show_dos(self):
  h = pickup_hamiltonian() # get hamiltonian
  if h.dimensionality==0:
    dos.dos0d(h,es=np.linspace(-3.1,3.1,500),delta=get("DOS_smearing"))
  elif h.dimensionality==1:
#    dos.dos1d(h,ndos=400,delta=get("DOS_smearing"))
    dos.dos1d(h,ndos=400)
  elif h.dimensionality==2:
    dos.dos2d(h,es=np.linspace(-3.1,3.1,500),delta=get("DOS_smearing"))
  else: raise
  execute_script("tb90-dos  ")


def pickup_hamiltonian():
  return initialize()
  if builder.get_object("activate_scf").get_active():
    return read_hamiltonian()
  else: # generate from scratch
    return initialize()



def read_hamiltonian():
  g = get_geometry0d() # get the geometry
  h = g.get_hamiltonian() # get the hamiltonian
  h.read("hamiltonian.in") # read hamiltonian
#  h.has_eh = builder.get_object("has_eh").get_active()
#  h.has_spin = builder.get_object("has_spin").get_active()
  return h







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
  nk = 100
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
#  h = g.get_hamiltonian() 
#  h.write()
  g.write()
  execute_script("qh-light-structure POSITIONS.OUT")
#  execute_script("qh-structure  ")




save_results = lambda x: save_outputs(inipath,tmppath) # function to save


# create signals
signals = dict()
signals["initialize"] = initialize  # initialize and run
signals["show_bands"] = show_bands  # show bandstructure
signals["show_structure"] = show_structure  # show bandstructure
signals["show_dos"] = show_dos  # show DOS
signals["show_berry2d"] = show_berry2d  # show DOS
#signals["show_stm"] = show_stm  # show STM
#signals["show_magnetism"] = show_magnetism  # show magnetism
#signals["show_lattice"] = show_lattice  # show magnetism
#signals["save_results"] = save_results  # save the results





#from qh_interface import create_folder # import all the libraries needed

window.connect_clicks(signals)
folder = create_folder()
tmppath = os.getcwd() # get the initial directory
window.run()

