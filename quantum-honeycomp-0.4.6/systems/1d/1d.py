#!/usr/bin/python


import sys
if len(sys.argv)>1: # if input provided
  sys.path.append(sys.argv[1]+ "pysrc")
  xmlpath = sys.argv[1] + "systems/1d/"  # path of the xml file
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
  lattice_name = builder.get_object("ribbon_type").get_active_text()
  if lattice_name=="Honeycomb zigzag":
    geometry_builder = geometry.honeycomb_zigzag_ribbon
  if lattice_name=="Honeycomb armchair":
    geometry_builder = geometry.honeycomb_armchair_ribbon
  if lattice_name=="Square (10)":
    geometry_builder = geometry.square_tetramer_ribbon
  g = geometry_builder(int(get("tetramers")))
  return g



def initialize(self):
  """ Initialize the calculation"""
  global g
  os.system("rm *.OUT")   # clean the outputs
  g = get_geometry() # get the geometry
  g = g.supercell(int(get("nsuper"))) # number of supercells
  g = sculpt.remove_central(g,int(get("nvac"))) # add the vacancies 
  h = g.get_hamiltonian() # get the hamiltonian
#  h.first_neighbors()  # first neighbor hoppin
  h.add_zeeman([get("Bx"),get("By"),get("Bz")]) # Zeeman fields
  h.add_sublattice_imbalance(get("AB_imbalance"))  # sublattice imbalance
  h.add_rashba(get("rashba_soc"))  # Rashba field
  h.add_antiferromagnetism(get("af_imbalance"))  # AF order
  h.add_kane_mele(get("km_soc")) # intrinsic SOC
  h.add_peierls(get("peierls")) # gauge field
#  h.add_transverse_efield(get("efield")) # electric field field
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
  g = get_geometry() # get the geometry
  h = g.get_hamiltonian() # get the hamiltonian
  h.read("hamiltonian.in") # read hamiltonian
  ldos.ldos1d(h,e=get("stm_bias"),delta=0.01,nrep=10) 
#  os.system("tb90-calculate-ldos "+str(get("stm_bias")))
  os.system("qh-ldos &")


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
#  os.system("qh-magnetism 3 &")





# create signals
signals = dict()
signals["on_window_destroy"] = gtk.main_quit  # close the window
signals["initialize"] = initialize  # initialize and run
signals["show_bands"] = show_bands  # show bandstructure
signals["show_dos"] = show_dos  # show DOS
signals["show_stm"] = show_stm  # show STM
signals["show_spin_response"] = show_spin_response  # show STM
signals["show_bdg_bands"] = show_bdg_bands  # show BdG bands
signals["show_magnetism"] = show_magnetism  # show magnetism

class RibbonApp(object):       
	def __init__(self):
	    builder.add_from_file(xmlpath+"1d.xml")
	    builder.connect_signals(signals)
	    self.window = builder.get_object("ribbon_window")
# creates a temporal folder in /tmp
            folder = create_folder()  
	    self.window.show()



#if __name__ == "__main__":
app = RibbonApp()  # launch the app
gtk.main()  # infinite loop

