#!/usr/bin/python2


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
  def skyrmion(rin):
    """Function which generates the magnetism"""
    x,y,z = rin[0],rin[1],rin[2]
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
  mag = [skyrmion(r) for r in g.r]
  scftypes.write_magnetization(mag) # write magnetization in file
  h.add_zeeman(skyrmion)   # add zeeman field to the hamiltonian 



def get_geometry2d():
  """ Create a 2d honeycomb lattice"""
  lattice_name = builder.get_object("lattice").get_active_text()
  if lattice_name=="Honeycomb":
    geometry_builder = geometry.honeycomb_lattice
  if lattice_name=="Square":
    geometry_builder = geometry.square_lattice
  if lattice_name=="Kagome":
    geometry_builder = geometry.kagome_lattice
  if lattice_name=="Rectangular Kagome":
    geometry_builder = geometry.rectangular_kagome_lattice
  if lattice_name=="Honeycomb 4 atoms":
    geometry_builder = geometry.honeycomb_lattice_square_cell
  if lattice_name=="Lieb":
    geometry_builder = geometry.lieb_lattice
  if lattice_name=="Triangular":
    geometry_builder = geometry.triangular_lattice
  g = geometry_builder()
  g = g.supercell(int(get("nsuper")))
  return g










def initialize(self):
  """ Initialize the calculation"""
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
  return h


def show_skyrmion(self):
  initialize(1)
  execute_script("qh-magnetism 1 nobonds  ") # show unit cell with magnetism


def show_bands(self):
  h = pickup_hamiltonian()  # get the hamiltonian
  opname = builder.get_object("bands_color").get_active_text()
  kpath = klist.default(h.geometry,nk=int(get("nkpoints")))  # write klist
  if opname=="None": # no operators
    h.get_bands(kpath=kpath)
  elif opname=="Sx": # off plane case
    h.get_bands(operator=operators.get_sx(h))
  elif opname=="Sy": # off plane case
    h.get_bands(operator=operators.get_sy(h))
  elif opname=="Sz": # off plane case
    h.get_bands(operator=operators.get_sz(h))
  elif opname=="Berry curvature":
    bandstructure.berry_bands(h,klist=kpath) # Berry curvature 
  elif opname=="Sz Berry curvature": # Sz Berry curvature 
    bandstructure.berry_bands(h,klist=kpath,operator=operators.get_sz(h))
  execute_script("qh-bands2d  ")




def show_dos(self):
  h = pickup_hamiltonian()  # get the hamiltonian
  dos.dos2d(h,use_kpm=False,nk=int(get("nkpoints")))
  execute_script("tb90-dos  ")
  return



def show_stm(self):
  execute_script("tb90-calculate-ldos "+str(get("stm_bias")))
  execute_script("tb90-stmmap  ")


def show_magnetism(self):
  h = pickup_hamiltonian() # get hamiltonian
  h.get_magnetization(nkp=int(get("nkpoints"))) # get the magnetization
  execute_script("qh-magnetism 2  ")




def show_chern(self):
  h = pickup_hamiltonian()  # get the hamiltonian
  topology.chern(h,nk=get("nkpoints")/4) # calculate chern number
  execute_script("tb90-chern")
  execute_script("qh-plotchern")


def show_ldos(self):
  h = pickup_hamiltonian()  # get the hamiltonian
  nk = int(2000/h.intra.shape[0]) + 1 # number of kpoints
  ldos.ldos2d(h,e=get("energy_ldos"),delta=0.01,nrep=3,nk=nk)
  execute_script("qh-ldos LDOS.OUT ")



def show_berry_2d(self):
  h = pickup_hamiltonian()  # get the hamiltonian
  topology.berry_map(h,nk=int(round(np.sqrt(get("nkpoints")))))
  execute_script("qh-berry2d BERRY_MAP.OUT")



def show_berry_1d(self):
  h = pickup_hamiltonian()  # get the hamiltonian
  ks = klist.default(h.geometry,nk=int(get("nkpoints")))  # write klist
  topology.write_berry(h,ks)
  execute_script("tb90-berry1d  label  ")



def pickup_hamiltonian():
    return initialize(1)





save_results = lambda x: save_outputs(inipath,tmppath) # function to save




# create signals
signals = dict()
signals["on_window_destroy"] = gtk.main_quit  # close the window
signals["initialize"] = initialize  # initialize and show
signals["show_skyrmion"] = show_skyrmion  # initialize and show
signals["show_bands"] = show_bands  # show bandstructure
signals["show_dos"] = show_dos  # show DOS
signals["show_chern"] = show_chern  # show BdG bands
signals["show_berry_1d"] = show_berry_1d  # show Berry curvature
signals["show_berry_2d"] = show_berry_2d  # show Berry curvature
signals["show_magnetism"] = show_magnetism  # show magnetism
signals["show_ldos"] = show_ldos  # show magnetism
signals["save_results"] = save_results  # save the results

class BulkApp(object):       
      def __init__(self):
            global tmppath # temporal path
            builder.add_from_file(xmlpath+"skyrmion2d.xml")
            builder.connect_signals(signals)
            self.window = builder.get_object("bulk_window")
            folder = create_folder()
            tmppath = os.getcwd() # get the initial directory
            print("Temporal path is",tmppath)
            self.window.show()



if __name__ == "__main__":
  inipath = os.getcwd() # get the initial directory
  print("Initial path is",inipath)
  app = BulkApp()
  gtk.main()

