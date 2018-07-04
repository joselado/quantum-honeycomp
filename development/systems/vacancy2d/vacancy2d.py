#!/usr/bin/python


import sys
if len(sys.argv)>1: # if input provided
  sys.path.append(sys.argv[1]+ "pysrc")
  xmlpath = sys.argv[1] + "systems/vacancy2d/"  # path of the xml file
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
#######################################


def get(name):
  """Get the value of a certain variable"""
  return float(builder.get_object(name).get_text())




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
  nk = int(np.sqrt(get("nkpoints"))) # number of kpoints
  scf = interactions.hubbardscf(h,U=U,nkp=nk,mf=mf,
            silent=False,mix=mix,maxerror=1e-5,filling=get("filling"))
  scf.hamiltonian.write() # write in a file
  np.save("MEAN_FIELD.npy",np.array(scf.mean_field)) # save in a file







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
  g = g.supercell(int(get("nsuper"))) # create supercell
  g = sculpt.remove_central(g,int(get("numvac"))) # removes several cen atoms
  return g











def initialize(self):
  """ Initialize the calculation"""
  g = get_geometry2d() # get the geometry
  h = hamiltonians.hamiltonian(g) # get the hamiltonian
  h.first_neighbors()  # first neighbor hoppin
  h.add_zeeman([get("Bx"),get("By"),get("Bz")]) # Zeeman fields
  h.add_sublattice_imbalance(get("mab"))  # sublattice imbalance
  h.add_rashba(get("rashba"))  # Rashba field
  h.add_antiferromagnetism(get("maf"))  # AF order
  h.add_kane_mele(get("kanemele")) # intrinsic SOC
  h.shift_fermi(get("fermi")) # shift fermi energy
  if builder.get_object("activate_scf").get_active():
    custom_scf(h) # create the tb90.in
  else:
    h.write()
  klist.default(g,nk=int(get("nkpoints")))  # write klist
#  klist.tr_path(nk=int(get("nkpoints")))  # write klist
  ## show the cell
  return h



def show_bands(self):
  h = pickup_hamiltonian()  # get the hamiltonian
  opname = builder.get_object("bands_color").get_active_text()
  kpath = klist.default(h.geometry,nk=int(get("nkpoints")))  # write klist
  if opname=="None": # no operators
    h.get_bands(kpath=kpath)
  elif opname=="Sz": # off plane case
    h.get_bands(kpath=kpath,operator=operators.get_sz(h))
  elif opname=="Sx":
    h.get_bands(kpath=kpath,operator=operators.get_sx(h))
  elif opname=="Sy":
    h.get_bands(kpath=kpath,operator=operators.get_sy(h))
  elif opname=="Sublattice":
    h.get_bands(kpath=kpath,operator=operators.get_sublattice(h))
  elif opname=="Berry curvature":
    bandstructure.berry_bands(h,klist=kpath) # Berry curvature 
  elif opname=="Sz Berry curvature": # Sz Berry curvature 
    bandstructure.berry_bands(h,klist=kpath,operator=operators.get_sz(h))
  else: raise
  execute_script("qh-bands2d  ")





def show_dos(self):
  h = pickup_hamiltonian()  # get the hamiltonian
  dos.dos2d(h,use_kpm=False,nk=int(get("nkpoints")))
  execute_script("tb90-dos  ")
  return



def show_stm(self):
  h = pickup_hamiltonian()  # get the hamiltonian
  ldos.ldos2d(h,e=get("e_ldos"),delta=0.01)
#  execute_script("tb90-calculate-ldos "+str(get("stm_bias")))
  execute_script("qh-interpolate LDOS.OUT ")
  execute_script("tb90-cmap LDOS.OUT-interpolated ")



  

def show_magnetism(self):
  h = pickup_hamiltonian() # get hamiltonian
  h.get_magnetization(nkp=int(round(np.sqrt(get("nkpoints"))))) # get the magnetization
  execute_script("qh-magnetism ")





def show_crystal(self):
  """Show the crystal structure"""
  g = get_geometry2d() # get the geometry
  g.write()
  execute_script("qh-magnetism 1 nomag nobonds  ")






def show_chern(self):
  h = pickup_hamiltonian()  # get the hamiltonian
  topology.chern(h,nk=get("nkpoints")/4) # calculate chern number
  execute_script("tb90-chern")
  execute_script("qh-plotchern")


def show_berry_2d(self):
  h = pickup_hamiltonian()  # get the hamiltonian
  topology.berry_map(h,nk=int(np.sqrt(get("nkpoints"))))
  execute_script("qh-berry2d BERRY_MAP.OUT")



def show_berry_1d(self):
  h = pickup_hamiltonian()  # get the hamiltonian
  ks = klist.default(h.geometry,nk=int(get("nkpoints")))  # write klist
  topology.write_berry(h,ks)
  execute_script("tb90-berry1d  label  ")



def pickup_hamiltonian():
  if builder.get_object("activate_scf").get_active():
    return read_hamiltonian()
  else: # generate from scratch
    return initialize(1)



def read_hamiltonian():
  return hamiltonians.load()


 


save_results = lambda x: save_outputs(inipath,tmppath) # function to save



# create signals
signals = dict()
signals["on_window_destroy"] = gtk.main_quit  # close the window
signals["initialize"] = initialize  # initialize and run
signals["show_bands"] = show_bands  # show bandstructure
signals["show_dos"] = show_dos  # show DOS
signals["show_chern"] = show_chern  # show BdG bands
signals["show_berry_1d"] = show_berry_1d  # show Berry curvature
signals["show_berry_2d"] = show_berry_2d  # show Berry curvature
signals["show_magnetism"] = show_magnetism  # show magnetism
signals["show_crystal"] = show_crystal  # show magnetism
signals["save_results"] = save_results

class BulkApp(object):       
     def __init__(self):
            global tmppath # temporal path
            builder.add_from_file(xmlpath+"vacancy2d.xml")
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

