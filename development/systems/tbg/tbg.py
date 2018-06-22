#!/usr/bin/python2


import sys
if len(sys.argv)>1: # if input provided
  sys.path.append(sys.argv[1]+ "pysrc")
  xmlpath = sys.argv[1] + "systems/tbg/"  # path of the xml file
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

def get_text(name):
  """Get the value of a certain variable"""
  return builder.get_object(name).get_text()



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
  n = int(get("size")) # size of the unit cell
  import specialgeometry
  g = specialgeometry.twisted_bilayer(n)
#  g = geometry.honeycomb_lattice()
#  g = g.supercell(n)
  return g



def initialize(self):
  """ Initialize the calculation"""
  g = get_geometry2d() # get the geometry
  from specialhopping import twisted,twisted_matrix
  has_spin = builder.get_object("has_spin").get_active()
  h = g.get_hamiltonian(is_sparse=True,has_spin=has_spin,is_multicell=False,
     mgenerator=twisted_matrix(ti=get("tinter"),lambi=7.0))
#  return h
#  h.turn_dense()
  h.add_sublattice_imbalance(get("mAB"))  # sublattice imbalance
  efield = get("interlayer_bias")
  def bias(r):
    if r[2]<0.0: return efield
    else: return -efield
  h.shift_fermi(bias)
  if h.has_spin:
    h.add_zeeman([get("Bx"),get("By"),get("Bz")]) # Zeeman fields
    h.add_rashba(get("rashba"))  # Rashba field
    h.add_antiferromagnetism(get("mAF"))  # AF order
    h.add_kane_mele(get("kanemele")) # intrinsic SOC
  h.shift_fermi(get("fermi")) # shift fermi energy
  if builder.get_object("has_eh").get_active():
    h.add_swave(get("swave"))
  if builder.get_object("activate_scf").get_active():
    if h.has_eh:
        print("SCF not implemented with Nambu")
        raise
  if builder.get_object("activate_scf").get_active():
    custom_scf(h) # create the tb90.in
  else:
    h.write("hamiltonian.in")
  klist.default(g,nk=int(get("nkpoints")))  # write klist
#  klist.tr_path(nk=int(get("nkpoints")))  # write klist
  return h



def show_bands(self):
  h = pickup_hamiltonian()  # get the hamiltonian
  if h.intra.shape[0]<4000:
    num_bands = int(get("num_bands"))
    if num_bands> h.intra.shape[0]-3: 
      h.turn_dense()
      num_bands = None
    if num_bands<0: 
      h.turn_dense()
      num_bands = None
  else: num_bands = max(20,int(get("num_bands")))
  opname = builder.get_object("bands_color").get_active_text()
  kpath = klist.default(h.geometry,nk=int(get("nkpoints")))  # write klist
  if opname=="None": # no operators
    h.get_bands(kpath=kpath,num_bands=num_bands) 
  elif opname=="Sx": # off plane case
    h.get_bands(operator=operators.get_sx(h))
  elif opname=="Sy": # off plane case
    h.get_bands(operator=operators.get_sy(h))
  elif opname=="Sz": # off plane case
    h.get_bands(operator=operators.get_sz(h))
  elif opname=="Sublattice": # off plane case
    h.get_bands(operator=operators.get_sublattice(h))
  elif opname=="Berry curvature":  
    bandstructure.berry_bands(h,klist=kpath) # Berry curvature 
  elif opname=="Sz Berry curvature": # Sz Berry curvature 
    bandstructure.berry_bands(h,klist=kpath,operator=operators.get_sz(h)) 
  else: raise
  execute_script("qh-bands2d ")
  




  

def show_dos(self):
  h = pickup_hamiltonian()  # get the hamiltonian
  nk = int(round(np.sqrt(get("nkpoints"))))
  dos.dos2d(h,use_kpm=False,nk=nk)
  execute_script("tb90-dos  ")
  return



def show_magnetism(self):
  h = pickup_hamiltonian() # get hamiltonian
  h.get_magnetization(nkp=int(get("nkpoints"))) # get the magnetization
  execute_script("qh-magnetism 2  ")
#  execute_script("qh-structure 2  ")





def show_chern(self):
  h = pickup_hamiltonian()  # get the hamiltonian
  nk = get("nkpoints")
  nk = int(round(np.sqrt(nk)))
  topology.chern(h,nk=nk) # calculate chern number
  execute_script("tb90-chern") 
#  execute_script("qh-plotchern") 



def show_berry_2d(self):
  h = pickup_hamiltonian()  # get the hamiltonian
  nk = get("nkpoints")
  nk = int(round(np.sqrt(nk)))
  topology.berry_map(h,nk=nk) 
  execute_script("qh-berry2d BERRY_MAP.OUT") 

 

def show_berry_1d(self):
  h = pickup_hamiltonian()  # get the hamiltonian
  ks = klist.default(h.geometry,nk=int(get("nkpoints")))  # write klist
  topology.write_berry(h,ks) 
  execute_script("tb90-berry1d  label  ")


def show_structure(self=None):
  g = get_geometry2d() # get the geometry
  g.write()
  execute_script("qh-potential POSITIONS.OUT ")



def pickup_hamiltonian():
  if builder.get_object("activate_scf").get_active():
    return read_hamiltonian()
  else: # generate from scratch
    return initialize(1)



def read_hamiltonian():
  g = get_geometry2d() # get the geometry
  h = g.get_hamiltonian() # get the hamiltonian
  h.read("hamiltonian.in") # read hamiltonian
  h.has_eh = builder.get_object("has_eh").get_active()
  h.has_spin = builder.get_object("has_spin").get_active()
  return h









def show_ldos(self):
  h = pickup_hamiltonian()  # get the hamiltonian
#  if h.intra.shape[0]<2000: h.turn_dense()
  e = get("energy_ldos")
  delta = get("delta_ldos")
  nk = get("nkpoints")/10
  nk = int(round(np.sqrt(nk)))
  ldos.ldos2d(h,e=e,delta=delta,nk=nk,mode="arpack")
  execute_script("qh-fast-ldos LDOS.OUT  ")
#  execute_script("qh-multildos ")
#  execute_script("tb90-calculate-ldos "+str(get("stm_bias")))
#  execute_script("qh-interpolate LDOS.OUT ")
#  execute_script("tb90-cmap LDOS.OUT-interpolated ")


def show_fermi_surface(self):
  h = pickup_hamiltonian()  # get the hamiltonian
#  spectrum.fermi_surface(h,reciprocal=True,nk=get("nkpoints")/2) # calculate Fs
  spectrum.boolean_fermi_surface(h,reciprocal=True,nk=get("nkpoints")/2) # calculate Fs
#  execute_script("tb90-cmap FERMI_MAP.OUT  ") # plot the result
  execute_script("tb90-cmap BOOL_FERMI_MAP.OUT  ") # plot the result

def show_z2_invariant(self):
  h = pickup_hamiltonian()  # get the hamiltonian
  nk = get("nkpoints")/4
  topology.z2_vanderbilt(h,nk=nk,nt=nk/2) # calculate z2 invariant
  execute_script("qh-wannier-center  ") # plot the result




def show_kdos(self):
  h = pickup_hamiltonian()  # get the hamiltonian
  ew = get("e_kdos")
  new = int(get("nkpoints")/10) # scale as kpoints
  energies = np.linspace(-ew,ew,new) # number of ene
  klist = np.linspace(0.,1.,new)
  kdos.write_surface_2d(h,energies=energies,delta=ew/new,klist=klist)
  execute_script("qh-kdos KDOS.OUT  ")





def show_2dband(self):
  h = pickup_hamiltonian()  # get the hamiltonian
  nk = get("nkpoints")/4
  ns = get_text("num_2dband") # get indexes of the bands
  if "," in ns: ns = [int(n) for n in ns.split(",")] # get the different numbers
  else: ns=[int(ns)] # single one
  if 0 in ns: # in case all eigenvalues wanted
    ns = [i+1 for i in range(h.intra.shape[0]//2)]
    ns += [-i for i in ns]
  spectrum.get_bands(h,nindex=ns,nk=nk,reciprocal=True)
  string = ""
  for n in ns: string += "BANDS2D__"+str(n)+".OUT "
  execute_script("qh-plot3d "+string +"  ")

save_results = lambda x: save_outputs(inipath,tmppath) # function to save

# create signals
signals = dict()
signals["on_window_destroy"] = gtk.main_quit  # close the window
signals["initialize"] = initialize  # initialize and run
signals["show_bands"] = show_bands  # show bandstructure
signals["show_dos"] = show_dos  # show DOS
signals["show_chern"] = show_chern  # show Chern number 
signals["show_berry_1d"] = show_berry_1d  # show Berry curvature
signals["show_berry_2d"] = show_berry_2d  # show Berry curvature
signals["show_ldos"] = show_ldos  # show Berry curvature
signals["show_fermi_surface"] = show_fermi_surface  # show Berry curvature
signals["show_z2_invariant"] = show_z2_invariant  # show Berry curvature
signals["show_magnetism"] = show_magnetism  # show magnetism
signals["show_structure"] = show_structure  # show magnetism
signals["show_2dband"] = show_2dband  # show magnetism
signals["show_kdos"] = show_kdos  # show kdos
signals["save_results"] = save_results  # save the results

class BulkApp(object):       
        def __init__(self):
            global tmppath # temporal path
            print(xmlpath)
            builder.add_from_file(xmlpath+"tbg.xml")
            builder.connect_signals(signals)
            self.window = builder.get_object("bulk_window")
            folder = create_folder()
            tmppath = os.getcwd() # get the initial directory
            self.window.show()



if __name__ == "__main__":
        inipath = os.getcwd() # get the initial directory
        app = BulkApp()
        gtk.main()

