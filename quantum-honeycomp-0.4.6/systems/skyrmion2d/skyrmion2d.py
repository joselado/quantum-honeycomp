#!/usr/bin/python


import sys
if len(sys.argv)>1: # if input provided
  sys.path.append(sys.argv[1]+ "pysrc")
  xmlpath = sys.argv[1] + "systems/skyrmion2d/"  # path of the xml file
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


def build_skyrmion(h):
  """ Creates the skyrmion and adds the zeeman field to the hamiltonian"""
  # write magnetization in file
  g = h.geometry # get geometry of the system
  rc = max(g.x)*get("rcut")  # cutoff radius
  def step(r):  # step function
    if r<rc:  return 0.  # no z magnetism
    else:  return 1.  # maximun z magnetism

  def linear(r):  # step function
    if r<rc:  return (-.5 + r/rc)*2.  # no z magnetism
    else:  return 1.  # maximun z magnetism


  #### define skyrmion function ###
  nsky = int(get("winding"))
  def skyrmion(x,y):
    """Function which generates the magnetism"""
    r = (x**2 + y**2)
    if builder.get_object("interpolation").get_active_text()=="Step":
      mz = step(np.sqrt(r)) # magnetism in z
    if builder.get_object("interpolation").get_active_text()=="Linear":
      mz = linear(np.sqrt(r)) # magnetism in z
    mxy = 1. - mz**2 # magnetism in-plane
    phi = np.arctan2(np.array(y),np.array(x))
    mr = mxy
    mx = mr * np.cos(nsky*phi) # magnetism in x
    my = mr * np.sin(nsky*phi) # magnetism in y
    m = np.array([mx,my,mz])*get("Jsky")  
    return m # return magnetization
  input_tb90.write_magnetization(skyrmion,g.x,g.y) # write magnetization in file
  h.add_zeeman(skyrmion)   # add zeeman field to the hamiltonian 



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
  intb90.kpoints.nkpoints = int(get("nkpoints")**2)
  intb90.bands.klist_bands = "from_file"
  intb90.bands.use_ewindow_bands = False
  # end type of operator
  intb90.write() # write input




def get_geometry2d():
  """ Create a 2d honeycomb lattice"""
  if builder.get_object("lattice").get_active_text()=="Honeycomb":
    geometry_builder = geometry.honeycomb_lattice
  if builder.get_object("lattice").get_active_text()=="Square":
    geometry_builder = geometry.square_lattice
  if builder.get_object("lattice").get_active_text()=="Kagome":
    geometry_builder = geometry.kagome_lattice
  if builder.get_object("lattice").get_active_text()=="Honeycomb 4 atoms":
    geometry_builder = geometry.honeycomb_lattice_square_cell
  g = geometry_builder()
  # create supercell
  g = g.supercell(int(get("nsuper")))
  return g








def initialize_and_show(self):
  """ Initialize the calculation"""
  global g
  g = get_geometry2d() # get the geometry
  h = hamiltonians.hamiltonian(g) # get the hamiltonian
  h.first_neighbors()  # first neighbor hoppin
  h.add_zeeman([get("Bx"),get("By"),get("Bz")]) # Zeeman fields
  h.add_sublattice_imbalance(get("mAB"))  # sublattice imbalance
  h.add_rashba(get("rashba"))  # Rashba field
  h.add_antiferromagnetism(get("mAF"))  # AF order
  if abs(get("kanemele")>0.0): h.add_kane_mele(get("kanemele")) # intrinsic SOC
  h.shift_fermi(get("fermi")) # shift fermi energy
  klist.default(g,nk=int(get("nkpoints")))  # write klist
  build_skyrmion(h) # add the zeeman field and write magnetization
  h.write("hamiltonian.in")  # write in file
  os.system("tb90-magnetism 1 nobonds &") # show unit cell with magnetism


def show_bands(self):
  intb90 = input_tb90.tb90in() # create input
  custom_tb90_bands() # create the tb90.in
  klist.default(g,nk=int(get("nkpoints")))  # write klist
  os.system("tb90.x") # run the calculation
  os.system("tb90-bands &")



def show_dos(self):
  intb90 = input_tb90.tb90in() # create input
  intb90.mode("nothing")  # no calculation
  intb90.mode("DOS")
  intb90.kpoints.nkpoints = int(get("nkpoints")**2)
  intb90.write() # write input
  os.system("tb90.x")
  os.system("tb90-dos &")


def show_stm(self):
  os.system("tb90-calculate-ldos "+str(get("stm_bias")))
  os.system("tb90-stmmap &")



def show_magnetism(self):
  intb90 = input_tb90.tb90in() # create input
  intb90.mode("nothing") # do nothing
  intb90.mode("magnetism") # do RPA
  intb90.write() # write input
  os.system("tb90.x")
  os.system("tb90-magnetism 2 &")
#  os.system("tb90-magnetism &")


def show_chern(self):
  intb90 = input_tb90.tb90in() # create input
  intb90.mode("nothing") # do nothing
  intb90.mode("berry") # do Berry calculation
  intb90.topology.klist_berry = "default"  
  intb90.kpoints.nkpoints = int(get("nkpoints")) # number of kpoints
  intb90.write() # write input
  os.system("tb90.x")
  os.system("tb90-chern") 
  os.system("tb90-cmap BERRY_CURVATURE.OUT") 
 

def show_berry(self):
  intb90 = input_tb90.tb90in() # create input
  klist.default(g,nk=int(get("nkpoints")))  # write klist
  intb90.mode("nothing") # do nothing
  intb90.mode("berry") # do Berry calculation
  intb90.topology.klist_berry = "from_file"  
  intb90.write() # write input
  os.system("tb90.x")
  os.system("tb90-berry1d  label &")





# create signals
signals = dict()
signals["on_window_destroy"] = gtk.main_quit  # close the window
signals["initialize"] = initialize_and_show  # initialize and show
signals["show_bands"] = show_bands  # show bandstructure
signals["show_dos"] = show_dos  # show DOS
signals["show_chern"] = show_chern  # show BdG bands
signals["show_berry"] = show_berry  # show Berry curvature
signals["show_magnetism"] = show_magnetism  # show magnetism

class BulkApp(object):       
	def __init__(self):
	    builder.add_from_file(xmlpath+"skyrmion2d.xml")
	    builder.connect_signals(signals)
	    self.window = builder.get_object("bulk_window")
            folder = create_folder()
	    self.window.show()



if __name__ == "__main__":
	app = BulkApp()
	gtk.main()

