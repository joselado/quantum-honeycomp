#!/usr/bin/python

fname = "interface1d"  # name of the folder

import sys
if len(sys.argv)>1: # if input provided
  sys.path.append(sys.argv[1]+ "pysrc")
  xmlpath = sys.argv[1] + "systems/"+fname+"/"  # path of the xml file
else: # default path
  import os
  sys.path.append(os.getcwd()+"/../../pysrc")
  xmlpath = ""  # path of the xml file




from gi.repository import Gtk as gtk
builder = gtk.Builder()

from qh_interface import * # import all the libraries needed
## These are the libraries imported
#######################################
##import subprocess
##import os
##import sys
##import hamiltonians
##import klist
##import geometry
##import input_tb90 
#######################################



def get(name):
  """Get the value of a certain variable"""
  return float(builder.get_object(name).get_text())



def custom_tb90_scf():
  """Modifies the properties of tb90.in input"""
  intb90 = input_tb90.tb90in() # create input
  intb90.mode("nothing")  # no calculation
  intb90.mode("SCF")  # no calculation
  # begin SCF options
  intb90.scf_convergence.mean_field_matrix = "from_file"
  scfin = builder.get_object("scf_initialization").get_active_text()
  import random
  rand = random.random
  ###################################
  print scfin
  if scfin == "Reconverge": # no new mean field
    pass
  else:  # if create new mean field matrix
    mf = interactions.initialization(g,scfin=scfin)  # create the mean field mat
    input_tb90.write_mean_field(mf,output_file="mean_field.in") # write in file
  # end SCF options
  intb90.kpoints.nkpoints = int(get("nkpoints"))
  intb90.electrons.filling = get("filling")  # get the filling
  intb90.scf_convergence.hubbard_scf = get("hubbard") # Hubbard
  intb90.scf_convergence.mix_coef = 0.8 # Hubbard
  intb90.write() # write input



def custom_tb90_bands():
  """Modifies the properties of tb90.in input"""
  intb90 = input_tb90.tb90in() # create input
  intb90.kpoints.nkpoints = int(get("nkpoints"))
  intb90.mode("nothing")  # no calculation
  intb90.mode("bands")   # bands calculation
  # type of operator for bands
  op = builder.get_object("bands_color").get_active_text()
  if op=="None": intb90.bands.bands_operators_option = "None"
  elif op=="Off-plane S": intb90.bands.bands_operators_option = "Sz"
  elif op=="In-plane S": intb90.bands.bands_operators_option = "Sx"
  elif op=="Position": intb90.bands.bands_operators_option = "index"
  else: intb90.bands.bands_operators_option = "from_file"
  intb90.kpoints.nkpoints = int(get("nkpoints"))
  intb90.bands.use_ewindow_bands = False
  # end type of operator
  intb90.write() # write input


def get_geometry():
  if builder.get_object("ribbon_type").get_active_text()=="Honeycomb zigzag":
    geometry_builder = geometry.honeycomb_zigzag_ribbon
  if builder.get_object("ribbon_type").get_active_text()=="Honeycomb armchair":
    geometry_builder = geometry.honeycomb_armchair_ribbon
  if builder.get_object("ribbon_type").get_active_text()=="Square (10)":
    geometry_builder = geometry.square_tetramer_ribbon
  g = geometry_builder(int(get("width")))
  return g



def initialize(self):
  """ Initialize the calculation"""
  global g
  os.system("rm *.OUT")   # clean the outputs
  g = get_geometry() # get the geometry
  h = hamiltonians.hamiltonian(g) # get the hamiltonian
  h.first_neighbors()  # first neighbor hoppin
  def ipol(x,y,z): # interpola
    if y<0.0: return -1.0
    else: return 1.0
  #### Zeeman ####
  def zee(x=0.0,y=0.0,z=0.0):  # zeeman field
    if y<0.0:
      return np.array([get("Bx_1"),get("By_1"),get("Bz_1")]) # b field
    else: 
      return np.array([get("Bx_2"),get("By_2"),get("Bz_2")]) # b field
  h.add_zeeman(zee) # Zeeman fields
  ### Sublatice imbalance ###
  def mab(x=0.0,y=0.0,z=0.0):  # sublattice imbalance
    if y<0.0: return get("mAB_1") # AB
    else: return get("mAB_2") # AB
  h.add_sublattice_imbalance(mab)  # sublattice imbalance
  ### rashba ###
#  h.add_rashba(get("rashba_soc"))  # Rashba field
#  h.add_antiferromagnetism(get("af_imbalance"))  # AF order
#  h.add_kane_mele(get("km_soc")) # intrinsic SOC

  ### Magnetic field ###
  def peierls_hack(x1=0.0,y1=0.0,x2=0.0,y2=0.0):  # sublattice imbalance
    phi = (x1-x2)*(y1+y2) # phase
    if (y1+y2)>0.0: return get("peierls_1")*phi # AB
    else: return get("peierls_2")*phi # AB
  h.add_peierls(peierls_hack) # gauge field




  if builder.get_object("activate_scf").get_active():
    h.write("hamiltonian_0.in")
    custom_tb90_scf() # create the tb90.in
    os.system("tb90.x") # run the mean field calculation
  else:
    h.write("hamiltonian.in")


def show_bands(self):
  intb90 = input_tb90.tb90in() # create input
  custom_tb90_bands() # create the tb90.in
  os.system("tb90.x") # run the calculation
  os.system("tb90-bands &")


def show_bdg_bands(self):
  h = hamiltonians.hamiltonian(get_geometry())
  h.read("hamiltonian.in")  # read the hamiltonian
# move to a temporal file
  os.system("mv hamiltonian.in hamiltonian-noBdG.in")  
  h.add_swave_electron_hole_pairing(get("pairing"))
  h.shift_fermi(get("fermi_bdg"))
  h.write()
# create tb90 input
  intb90 = input_tb90.tb90in() # create input
  intb90.kpoints.nkpoints = int(get("nkpoints"))
  intb90.mode("nothing")  # no calculation
  intb90.mode("bands")   # bands calculation
  intb90.bands.bands_operators_option = "None"
  intb90.write()
  os.system("tb90.x ") 
  os.system("tb90-bands &")
# recover from temporal file
  os.system("mv hamiltonian.in hamiltonian-BdG.in")  
  os.system("mv hamiltonian-noBdG.in hamiltonian.in")

  

def show_dos(self):
  intb90 = input_tb90.tb90in() # create input
  intb90.mode("nothing")  # no calculation
  intb90.mode("DOS")
  intb90.kpoints.nkpoints = 10*int(get("nkpoints"))
  intb90.write() # write input
  os.system("tb90.x")
  os.system("tb90-dos &")


def show_stm(self):
  os.system("tb90-calculate-ldos "+str(get("stm_bias")))
  os.system("tb90-stmmap &")


def show_spin_response(self):
  intb90 = input_tb90.tb90in() # create input
  intb90.mode("nothing") # do nothing
  intb90.mode("spin-RPA") # do RPA
  intb90.dielectric.q_chi = get("qvector") # qvector in RPA
  intb90.scf_convergence.hubbard_scf = get("hubbard") # Hubbard
  intb90.kpoints.nkpoints = int(get("nkpoints"))
  intb90.write() # write input
  os.system("tb90.x")
  os.system("tb90-imchi &")
  

def show_magnetism(self):
  intb90 = input_tb90.tb90in() # create input
  intb90.mode("nothing") # do nothing
  intb90.mode("magnetism") # do RPA
  intb90.write() # write input
  os.system("tb90.x")
  os.system("tb90-magnetism 3 &")





# create signals
signals = dict()
signals["on_window_destroy"] = gtk.main_quit  # close the window
signals["initialize"] = initialize  # initialize and run
signals["show_bands"] = show_bands  # show bandstructure
signals["show_dos"] = show_dos  # show DOS
signals["show_stm"] = show_stm  # show STM
#signals["show_spin_response"] = show_spin_response  # show STM
#signals["show_bdg_bands"] = show_bdg_bands  # show BdG bands
signals["show_magnetism"] = show_magnetism  # show magnetism

class RibbonApp(object):       
	def __init__(self):
	    builder.add_from_file(xmlpath+fname+".xml")
	    builder.connect_signals(signals)
	    self.window = builder.get_object("hybrid_window")
# creates a temporal folder in /tmp
            folder = create_folder()  
	    self.window.show()



#if __name__ == "__main__":
app = RibbonApp()  # launch the app
gtk.main()  # infinite loop

