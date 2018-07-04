#!/usr/bin/python


import sys
if len(sys.argv)>1: # if input provided
  sys.path.append(sys.argv[1]+ "pysrc")
  xmlpath = sys.argv[1] + "systems/grain_boundary/"  # path of the xml file
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
  def magf(r):
    if abs(r[2])<2: return [0.,0.,0.1]
    else: return [0.0,0.,0.]

  scfin = builder.get_object("scf_initialization").get_active_text()
  mf = scftypes.guess(h,"ferro",fun=magf)
  filling = (len(h.geometry.r)-1.0+get("extra_electrons"))/(2.0*len(h.geometry.r))
  nk = int(get("nkpoints")/10)
  scf = scftypes.selfconsistency(h,nkp=nk,filling=filling,mf=mf,
               mode="Hubbard collinear",g=get("hubbard"))
  return scf.hamiltonian
  






def get_geometry():
  n = int(get("width"))
  g1 = geometry.honeycomb_zigzag_ribbon(n)
  g2 = geometry.honeycomb_armchair_ribbon(n*2)
  g1 = g1.supercell(2)
  g2 = g2.supercell(1)
  fac = g2.a1[0]/g1.a1[0]
  
  print("Compression factor",fac)
  
  g = g1.copy()
  g.has_sublattice = False
  
  g1.y -= np.min(g1.y) - 0.9
  g2.x /= fac
  g1.x += -0.44
  g1.xyz2r()
  g2.y -= np.max(g2.y)
  g2.xyz2r()
  
  g.r = np.concatenate([g1.r,g2.r])
  g.r2xyz()
  g.write()
  return g





def is_neigh(r1,r2):
  """Check if this is a neighbor"""
  dr = r1-r2
  dr = dr.dot(dr)
  if 0.7<dr<1.4: return True
  else: return False




def initialize(self):
  """ Initialize the calculation"""
#  global g
  global hscf # global hamiltonian
  os.system("rm *.OUT")   # clean the outputs
  
  def fun(r1,r2):
    if is_neigh(r1,r2): return -1.0
    else: return 0.0

  g = get_geometry() # get the geometry
  h = g.get_hamiltonian(is_multicell=False,fun=fun,has_spin=False) # get the hamiltonian
  # remove the zigzag state
  def get_switch(name): return builder.get_object(name).get_active()
  if get_switch("remove_zz"):
    def ferf(r):
      maxy = np.max(h.geometry.y)
      if r[1]>(maxy-0.5): return +1.0
      else: return 0.0
  
    h.shift_fermi(ferf)


  if get_switch("has_spin"):  
    h.turn_spinful() # turn spinful
  h.add_zeeman([get("Bx"),get("By"),get("Bz")]) # Zeeman fields
#  h.add_sublattice_imbalance(get("AB_imbalance"))  # sublattice imbalance
#  h.add_rashba(get("rashba_soc"))  # Rashba field
#  h.add_antiferromagnetism(get("af_imbalance"))  # AF order
#  h.add_kane_mele(get("km_soc")) # intrinsic SOC
  h.add_peierls(get("peierls")) # gauge field
  if get_switch("has_spin"):  
    if get_switch("has_eh"):
      h.add_swave(get("swave"))
#  else: # check if spin should be removed
#    if not builder.get_object("has_spin").get_active():
#      h.remove_spin()
#  h.add_transverse_efield(get("efield")) # electric field field
  if get_switch("activate_scf"):
    h.turn_spinful()
    if h.has_eh: 
        print("SCF not implemented with Nambu")
        raise
    h = custom_scf(h) 
#  else:

  hscf = h.copy()
  return h # return the Hamiltonian


def pickup_hamiltonian():
  if builder.get_object("activate_scf").get_active():
    return hscf.copy()
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
    h.get_bands(kpath=kpath,nkpoints=nkpoints)
  elif opname=="Sx": # off plane case
    h.get_bands(operator=operators.get_sx(h))
  elif opname=="Sy": # off plane case
    h.get_bands(operator=operators.get_sy(h))
  elif opname=="Sz": # off plane case
    h.get_bands(operator=operators.get_sz(h))
  elif opname=="Position":
    h.get_bands(operator=operators.bulk1d(h).todense(),nkpoints=nkpoints)
  elif opname=="Interface":
    h.get_bands(operator=h.get_operator("interface"),nkpoints=nkpoints)
  elif opname=="Sublattice":
    h.get_bands(operator=operators.get_sublattice(h).todense(),nkpoints=nkpoints)
  elif opname=="Current": # current density
    bandstructure.current_bands(h,klist=kpath)
  else: raise
  execute_script("qh-bands1d  ")


  

def show_dos(self):
  h = pickup_hamiltonian()  # get the hamiltonian
  nk = int(get("nkpoints"))
  dos.dos1d(h,use_kpm=False,nk=nk,ndos=50*nk)
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
  execute_script("tb90-magnetism 2 nobonds  ")

def show_structure(self):
  g = get_geometry()
#  execute_script("qh-structure 2  ")
  execute_script("tb90-magnetism 3 nomag nobonds  ")



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
    builder.add_from_file(xmlpath+"grain_boundary.xml")
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

