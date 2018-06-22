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
#######################################



def get(name):
  """Get the value of a certain variable"""
  return float(builder.get_object(name).get_text())


def custom_scf(h):
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
  noncol = builder.get_object("noncollinear").get_active()
  scf = interactions.hubbardscf(h,U=U,nkp=int(get("nkpoints")//10),mf=mf,
            silent=False,mix=mix,maxerror=1e-5,filling=get("filling"),
            collinear= not noncol)
  scf.hamiltonian.write() # write in a file
  np.save("MEAN_FIELD.npy",np.array(scf.mean_field)) # save in a file
  return scf.hamiltonian




def get_geometry():
  lattice_name = builder.get_object("ribbon_type").get_active_text()
  import ribbon
  def chiral_generator(gfun):
    def geometry_builder(n):
      g = gfun()
      cx,cy = builder.get_object("chirality").get_text().split(",")
      cx,cy = int(cx),int(cy)
      m = [[cx,cy,0],[0,1,0],[0,0,1]]
      ##################################
      # this is a dirty workaround
      g.has_sublattice = False
      g.z = 0.1*np.array(g.sublattice) # workaround
      g = supercell.non_orthogonal_supercell(g,m)
      g.has_sublattice = True
      g.sublattice = np.sign(g.z)
      g.z *= 0.0 # set to zero
      ################################
      g = g.supercell((1,n))
      g.dimensionality = 1
      return g
    return geometry_builder
  if lattice_name=="Honeycomb zigzag":
    geometry_builder = geometry.honeycomb_zigzag_ribbon
  elif lattice_name=="Honeycomb armchair":
    geometry_builder = geometry.honeycomb_armchair_ribbon
  elif lattice_name=="Square (10)":
    geometry_builder = geometry.square_tetramer_ribbon
  elif lattice_name=="Square chiral":
    geometry_builder = chiral_generator(geometry.square_lattice)
  elif lattice_name=="Honeycomb chiral":
    geometry_builder = chiral_generator(geometry.honeycomb_lattice)
  elif lattice_name=="Triangular":
    def geometry_builder(n): 
      g = geometry.triangular_lattice()
      return ribbon.bulk2ribbon(g,n=n)
  elif lattice_name=="Kagome":
    def geometry_builder(n): 
      g = geometry.kagome_lattice()
      return ribbon.bulk2ribbon(g,n=n)
  g = geometry_builder(int(get("tetramers")))
  g = g.supercell(int(get("nsuper"))) # number of supercells
  g = sculpt.remove_central(g,int(get("nvac"))) # add the vacancies 
  n_messed = int(get("messed_edge")) # number of messed edge atoms
  mess = sculpt.get_furthest(g,n=2*n_messed,tol=400) # get edge atoms
  g = sculpt.remove(g,mess) # remove messed atoms
  g.write()
  return g



def initialize(self):
  """ Initialize the calculation"""
#  global g
  global hscf
  os.system("rm *.OUT")   # clean the outputs
  g = get_geometry() # get the geometry
  h = g.get_hamiltonian(is_multicell=False) # get the hamiltonian
#  h.first_neighbors()  # first neighbor hoppin
  h.add_zeeman([get("Bx"),get("By"),get("Bz")]) # Zeeman fields
  h.add_sublattice_imbalance(get("AB_imbalance"))  # sublattice imbalance
  h.add_rashba(get("rashba_soc"))  # Rashba field
  h.add_antiferromagnetism(get("af_imbalance"))  # AF order
  h.add_kane_mele(get("km_soc")) # intrinsic SOC
  h.add_peierls(get("peierls")) # gauge field
  if builder.get_object("has_eh").get_active():
    h.add_swave(get("swave"))
  else: # check if spin should be removed
    if not builder.get_object("has_spin").get_active():
      h.remove_spin()
#  h.add_transverse_efield(get("efield")) # electric field field
  if builder.get_object("activate_scf").get_active():
    if h.has_eh: 
        print("SCF not implemented with Nambu")
        raise
    hscf = custom_scf(h) 
  else:
    h.write("hamiltonian.in")
  return h # return the Hamiltonian


def pickup_hamiltonian():
  global hscf
  if builder.get_object("activate_scf").get_active():
    return hscf.copy()
  else: # generate from scratch
    return initialize(1)
  


def read_hamiltonian():
  g = get_geometry() # get the geometry
  h = g.get_hamiltonian() # get the hamiltonian
  h.read("hamiltonian.in") # read hamiltonian
  h.has_eh = builder.get_object("has_eh").get_active()
  h.has_spin = builder.get_object("has_spin").get_active()
  return h


def show_bands(self):
  h = pickup_hamiltonian()  # get the hamiltonian
  opname = builder.get_object("bands_color").get_active_text()
  kpath = klist.default(h.geometry,nk=int(get("nkpoints")))  # write klist
  nkpoints = get("nkpoints")
  if opname=="None": # no operators
    h.get_bands(kpath=kpath,nkpoints=nkpoints)
  elif opname=="Sx": # off plane case
    h.get_bands(operator=operators.get_sx(h))
  elif opname=="Sy": # off plane case
    h.get_bands(operator=operators.get_sy(h))
  elif opname=="Sz": # off plane case
    h.get_bands(operator=operators.get_sz(h))
  elif opname=="Position":
    h.get_bands(operator=operators.bulk1d(h).todense(),nkpoints=nkpoints)
  elif opname=="Sublattice":
    h.get_bands(operator=operators.get_sublattice(h).todense(),nkpoints=nkpoints)
  elif opname=="Current": # current density
    bandstructure.current_bands(h,klist=kpath)
  else: raise
  execute_script("qh-bands1d  ")


  

def show_dos(self):
  h = pickup_hamiltonian()  # get the hamiltonian
  dos.dos1d(h,use_kpm=False,nk=int(get("nkpoints")))
  execute_script("tb90-dos  ")
  return



def show_stm(self):
  h = pickup_hamiltonian()  # get the hamiltonian
  ewin = get("stm_bias")
  nrep = int(len(h.geometry.r)//5) + 1
  nk = get("nkpoints")
  ne = nk # number of energies
  delta = 5./nk
  ldos.multi_ldos(h,es=np.linspace(-ewin,ewin,ne),nk=nk,delta=delta,nrep=nrep)
  execute_script("qh-multildos ")
#  ldos.ldos1d(h,e=get("stm_bias"),delta=0.01,nrep=10) 
#  execute_script("tb90-calculate-ldos "+str(get("stm_bias")))
#  execute_script("qh-ldos  LDOS.OUT")


  

def show_magnetism(self):
  h = pickup_hamiltonian() # get hamiltonian
  h.get_magnetization() # get the magnetization
  execute_script("qh-magnetism 2  ")

def show_structure(self):
  g = get_geometry()
#  execute_script("qh-structure 2  ")
  execute_script("qh-magnetism 1 nomag nobonds  ")



save_results = lambda x: save_outputs(inipath,tmppath) # function to save



# create signals
signals = dict()
signals["on_window_destroy"] = gtk.main_quit  # close the window
signals["initialize"] = initialize  # initialize and run
signals["show_bands"] = show_bands  # show bandstructure
signals["show_dos"] = show_dos  # show DOS
signals["show_stm"] = show_stm  # show STM
signals["show_magnetism"] = show_magnetism  # show magnetism
signals["show_structure"] = show_structure  # show structure
signals["save_results"] = save_results  # save the results

class RibbonApp(object):       
  def __init__(self):
    global tmppath # temporal path
    builder.add_from_file(xmlpath+"1d.xml")
    builder.connect_signals(signals)
    self.window = builder.get_object("ribbon_window")
# creates a temporal folder in /tmp
    folder = create_folder()  
    tmppath = os.getcwd() # get the initial directory
    print("Temporal path is",tmppath)
    self.window.show()



if __name__ == "__main__":
  inipath = os.getcwd() # get the initial directory
  print("Initial path is",inipath)
  app = RibbonApp()  # launch the app
  gtk.main()  # infinite loop

