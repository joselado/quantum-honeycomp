#!/usr/bin/python

from __future__ import print_function

import sys
if len(sys.argv)>1: # if input provided
  sys.path.append(sys.argv[1]+ "pysrc")
  xmlpath = sys.argv[1] + "systems/0d/"  # path of the xml file
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
  return scf.hamiltonian




def get_geometry0d():
  """ Create a 0d island"""
  lattice_name = builder.get_object("lattice").get_active_text()
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
  # first create a raw unit cell
  g = geometry_builder()  # build a 2d unit cell
  n = 1+int(round(get("size")))   # get the desired size
  nf = get("size")   # get the desired size, in float
  g = g.supercell(6*n)   # create supercell
  g.set_finite() # set as finite system
  g.center() # center the geometry
  # now scuplt the geometry
  nedges = int(get("nedges")) # number of edges
  if nedges>0: # if edges are provided
    g = sculpt.rotate(g,get("rotation")*2.*np.pi/360) # initial rotation
    def f(r): return r[0]>-nf*(np.cos(np.pi/3)+1.)  # function to use as cut
    for i in range(nedges): # loop over rotations, 60 degrees
      g = sculpt.intersec(g,f) # retain certain atoms
      g = sculpt.rotate(g,2.*np.pi/nedges) # rotate 60 degrees
    g = sculpt.remove_unibonded(g)  # remove single bonded atoms
    g.center() # center the geometry
    g = sculpt.remove_central(g,int(get("numvac"))) # removes several cen atoms
  else: # if number of edges nor defined, user the formula defined
    print("Using function")
    expression = builder.get_object("geometry_formula").get_text()
    fo = open(os.getcwd()+"/geometry_formula.py","w") # write in file
    fo.write("import geometry\n\n\n") # define function
    fo.write("import sculpt\n\n\n") # define function
    fo.write("def sculpt_fun(x,y):\n") # define function
    fo.write("  if "+expression+":  return True\n") # define function
    fo.write("  else:  return False\n\n") # define function
    fo.write("g = geometry.read()\n") # define function
    fo.write("g = sculpt.intersec(g,sculpt_fun)\n") # define function
    fo.write("g.write()\n") # define function
    fo.close()
    g.write()
#    execute_script("pygra-update")
    execute_script("python geometry_formula.py")
    g = geometry.read()
    g.dimensionality = 0
  return g











def initialize(self):
  """ Initialize the calculation"""
  global g,hscf
  g = get_geometry0d() # get the geometry
  h = hamiltonians.hamiltonian(g) # get the hamiltonian
  h.first_neighbors()  # first neighbor hoppin
  h.add_zeeman([get("Bx"),get("By"),get("Bz")]) # Zeeman fields
  h.add_sublattice_imbalance(get("mAB"))  # sublattice imbalance
  if abs(get("rashba")) > 0.0: h.add_rashba(get("rashba"))  # Rashba field
  h.add_antiferromagnetism(get("mAF"))  # AF order
  if abs(get("kanemele"))>0.0:  h.add_kane_mele(get("kanemele")) # intrinsic SOC
  h.shift_fermi(get("fermi")) # shift fermi energy
  h.add_peierls(get("peierls")) # shift fermi energy
  h.write() # write the Hamiltonian in file
  if builder.get_object("activate_scf").get_active():
    hscf = custom_scf(h) # perform selfconsistent calculation
  else: pass
  h.save() # save hamiltonian
  return h


def show_bands(self):
  h = pickup_hamiltonian() # get hamiltonian
  h.get_bands()
  opname = builder.get_object("bands_color").get_active_text()
  if opname=="None": # no operators
    h.get_bands()
  elif opname=="Sx": # off plane case
    h.get_bands(operator=operators.get_sx(h))
  elif opname=="Sy": # off plane case
    h.get_bands(operator=operators.get_sy(h))
  elif opname=="Sz": # off plane case
    h.get_bands(operator=operators.get_sz(h))

  execute_script("qh-bands0d  ")


  

def show_dos(self):
  h = pickup_hamiltonian() # get hamiltonian
  dos.dos0d(h,es=np.linspace(-3.1,3.1,500),delta=get("DOS_smearing"))
  execute_script("tb90-dos  ")


def pickup_hamiltonian():
  global hscf
  if builder.get_object("activate_scf").get_active():
    return hscf.copy()
  else: # generate from scratch
    return initialize(1)



def read_hamiltonian():
  return hamiltonians.load() # load Hamiltonian







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



def show_magnetism(self):
  h = pickup_hamiltonian() # get hamiltonian
  h.get_magnetization() # get the magnetization
  execute_script("tb90-magnetism  ")
#  execute_script("qh-magnetism  ")


def show_lattice(self):
  """Show the lattice of the system"""
  g = get_geometry0d() # get the geometry
  h = g.get_hamiltonian() 
  h.write()
  execute_script("qh-structure")
#  execute_script("qh-structure  ")




save_results = lambda x: save_outputs(inipath,tmppath) # function to save


# create signals
signals = dict()
signals["on_window_destroy"] = gtk.main_quit  # close the window
signals["initialize"] = initialize  # initialize and run
signals["show_bands"] = show_bands  # show bandstructure
signals["show_dos"] = show_dos  # show DOS
signals["show_stm"] = show_stm  # show STM
signals["show_magnetism"] = show_magnetism  # show magnetism
signals["show_lattice"] = show_lattice  # show magnetism
signals["save_results"] = save_results  # save the results

class BulkApp(object):       
       def __init__(self):
            global tmppath
            builder.add_from_file(xmlpath+"0d.xml")
            builder.connect_signals(signals)
            self.window = builder.get_object("bulk_window")
            folder = create_folder()
            tmppath = os.getcwd() # get the initial directory
            self.window.show()



if __name__ == "__main__":
        inipath = os.getcwd() # get the initial directory
        print("Initial path is",inipath)
        app = BulkApp()
        gtk.main()

