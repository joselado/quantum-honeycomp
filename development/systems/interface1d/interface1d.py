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



def get(name):
  """Get the value of a certain variable"""
  return float(builder.get_object(name).get_text())








def get_geometry():
  name = builder.get_object("ribbon_type").get_active_text()
  if name=="Honeycomb zigzag":
    geometry_builder = geometry.honeycomb_zigzag_ribbon
  if name=="Honeycomb armchair":
    geometry_builder = geometry.honeycomb_armchair_ribbon
  if name=="Square (10)":
    geometry_builder = geometry.square_tetramer_ribbon
  if name=="Square (11)":
    gb = geometry.square_lattice() # bulk square lattice
    import supercell
    gb = supercell.non_orthogonal_supercell(gb,
                     [[1,1,0],[-1,1,0],[0,0,1]],mode="brute")
    def geometry_builder(n):
      g = gb.supercell((1,n))
      g.dimensionality = 1
      g = sculpt.rotate_a2b(g,g.a1,np.array([1.,0.,0.]))
      return g
  g = geometry_builder(int(get("width")))
  return g





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
  scf = interactions.hubbardscf(h,U=U,nkp=int(get("nkpoints")),mf=mf,
            silent=False,mix=mix,maxerror=1e-5,filling=get("filling"))
  scf.hamiltonian.write() # write in a file
  np.save("MEAN_FIELD.npy",np.array(scf.mean_field)) # save in a file









def initialize(self):
  """ Initialize the calculation"""
  g = get_geometry() # get the geometry
  os.system("rm -f *.OUT")   # clean the outputs
  def create_h(name):
    g = get_geometry() # get the geometry
    h = g.get_hamiltonian(is_multicell=False) # get the hamiltonian
    h.add_zeeman([get("Bx"+name),get("By"+name),get("Bz"+name)]) # Zeeman fields
    h.add_sublattice_imbalance(get("mAB"+name))  # sublattice imbalance
    h.add_rashba(get("rashba"+name))  # Rashba field
    h.add_antiferromagnetism(get("mAF"+name))  # AF order
    h.add_kane_mele(get("km_soc"+name)) # intrinsic SOC
    h.add_peierls(get("peierls"+name)) # gauge field
    h.shift_fermi(get("fermi"+name)) # Fermi energy
    if builder.get_object("has_eh").get_active():
      h.add_swave(get("swave"+name))
    if builder.get_object("activate_scf").get_active():
      if h.has_eh:
          print("SCF not implemented with Nambu")
          raise
  
    return h
  # get the two Hamiltonians
  h1 = create_h("_1")
  h2 = create_h("_2")
  tlen = get("tran_length") # transition length
  def fun(r):
    if r[1]<-tlen: return 0.0
    elif r[1]>tlen: return 1.0
    else: return (r[1]+tlen)/(2*tlen) # interpolate
  h = hybrid.half_and_half(h1,h2,fun=fun) # create hybrid Hamiltonian
  h.write() # write Hamiltonian
  return h



  if builder.get_object("activate_scf").get_active():
    custom_scf(h) # cr
  else:
    h.write()




def pickup_hamiltonian():
  if builder.get_object("activate_scf").get_active():
    return read_hamiltonian()
  else: # generate from scratch
    return initialize(1)



def read_hamiltonian():
  return hamiltonians.load()








def show_bands(self):
  h = pickup_hamiltonian()  # get the hamiltonian
  opname = builder.get_object("bands_color").get_active_text()
  kpath = klist.default(h.geometry,nk=int(get("nkpoints")))  # write klist
  nkpoints = get("nkpoints")
  if opname=="None": # no operators
    h.get_bands(kpath=kpath)
  elif opname=="Sx": # off plane case
    h.get_bands(operator=operators.get_sx(h),nkpoints=nkpoints)
  elif opname=="Sy": # off plane case
    h.get_bands(operator=operators.get_sy(h),nkpoints=nkpoints)
  elif opname=="Sz": # off plane case
    h.get_bands(operator=operators.get_sz(h),nkpoints=nkpoints)
  elif opname=="Position":
    h.get_bands(operator=operators.bulk1d(h),nkpoints=nkpoints)
  elif opname=="Sublattice":
    h.get_bands(operator=operators.get_sublattice(h),nkpoints=nkpoints)
  elif opname=="Interface":
    h.get_bands(operator=operators.interface1d(h,cut=6.).todense(),nkpoints=nkpoints)
  execute_script("qh-bands1d  ")



  
def show_dos(self):
  h = pickup_hamiltonian()  # get the hamiltonian
  dos.dos1d(h,use_kpm=False,nk=int(get("nkpoints")))
  execute_script("tb90-dos  ")
  return




def show_stm(self):
  h = pickup_hamiltonian()  # get the hamiltonian
  g = h.geometry
  ldos.ldos1d(h,e=get("stm_bias"),delta=0.01,nrep=int(len(g.y)/2))
  execute_script("qh-ldos  LDOS.OUT")

def show_structure(self):
  h = pickup_hamiltonian()  # get the hamiltonian
  execute_script("tb90-magnetism 2")


  


def show_magnetism(self):
  h = pickup_hamiltonian() # get hamiltonian
  h.get_magnetization(nkp=int(get("nkpoints"))) # get the magnetization
  execute_script("tb90-magnetism  ")



save_results = lambda x: save_outputs(inipath,tmppath) # function to save




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
signals["show_structure"] = show_structure  # show magnetism
signals["save_results"] = save_results  # save the results

class RibbonApp(object):       
     def __init__(self):
            global tmppath # temporal path
            builder.add_from_file(xmlpath+fname+".xml")
            builder.connect_signals(signals)
            self.window = builder.get_object("hybrid_window")
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

