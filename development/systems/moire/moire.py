#!/usr/bin/python2


import sys
if len(sys.argv)>1: # if input provided
  sys.path.append(sys.argv[1]+ "pysrc")
  xmlpath = sys.argv[1] + "systems/moire/"  # path of the xml file
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
  g = sculpt.rotate_a2b(g,g.a1,np.array([1.0,0.0,0.0])) # 
  g = g.supercell(int(get("nsuper"))) # supercell
  return g





def get_function(name,g):
  """Get the function to create the potential for a certain parameter"""
  a = np.sqrt(g.a1.dot(g.a1)) # get the lattice constant
  val = get(name) # value of the potential
  n = int(round(get("Cn_"+name))) # get the overall symmetry
  mode = get("mode_"+name) # get the mode
  mfun = potentials.cnpot(n=n,k=2*2*np.pi/a*mode,v=val)
 
  def dummyfun(r1,r2=None): # dummy function in case there are two arguments
    if r2 is None: return mfun(r1)
    else: return mfun((r1+r2)/2.)
  ms = [mfun(ir) for ir in g.r]
  np.savetxt("POTENTIAL_"+name+".OUT",np.matrix([g.x,g.y,ms]).T)
  return dummyfun
 








def initialize(self):
  """ Initialize the calculation"""
  g = get_geometry2d() # get the geometry
  h = hamiltonians.hamiltonian(g) # get the hamiltonian
  h.first_neighbors()  # first neighbor hoppin

  fbx = get_function("Bx",g)
  fby = get_function("By",g)
  fbz = get_function("Bz",g)
  h.add_zeeman(lambda r: [fbx(r),fby(r),fbz(r)]) # Zeeman fields
  h.add_sublattice_imbalance(get_function("mAB",g))  # sublattice imbalance
  h.add_rashba(get_function("rashba",g))  # Rashba field
  h.add_antiferromagnetism(get_function("mAF",g))  # AF order
  if abs(get("kanemele"))>0.0:
    h.add_kane_mele(get_function("kanemele",g)) # intrinsic SOC
  h.shift_fermi(get_function("fermi",g)) # shift fermi energy
  if builder.get_object("has_eh").get_active():
    h.add_swave(get_function("swave",g))
    h.write("hamiltonian.in")
  else:
    h.write("hamiltonian.in")
  klist.default(g,nk=int(get("nkpoints")))  # write klist
#  klist.tr_path(nk=int(get("nkpoints")))  # write klist
  return h



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
  dos.dos2d(h,use_kpm=False,nk=int(get("nkpoints")))
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
  execute_script("tb90-cmap BERRY_MAP.OUT") 

 

def show_berry_1d(self):
  h = pickup_hamiltonian()  # get the hamiltonian
  ks = klist.default(g,nk=int(get("nkpoints")))  # write klist
  topology.write_berry(h,ks) 
  execute_script("tb90-berry1d  label  ")






def pickup_hamiltonian():
    return initialize(1)








def show_ldos(self):
  h = pickup_hamiltonian()  # get the hamiltonian
  ldos.ldos2d(h,e=get("e_ldos"),delta=0.01,nk=4)
#  execute_script("tb90-calculate-ldos "+str(get("stm_bias")))
  execute_script("qh-ldos LDOS.OUT ")
#  execute_script("qh-ldos-bond LDOS.OUT ")


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


def plot_potential(self,name=""):
  """Plot a certain poential"""
  print("plotting",name)
  g = get_geometry2d() # get the geometry
  get_function(name,g) # get this function
  filename = "POTENTIAL_"+name+".OUT"
  execute_script("qh-potential "+filename) # plot this potential
  





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
signals["show_2dband"] = show_2dband  # show magnetism
signals["show_kdos"] = show_kdos  # show kdos
signals["save_results"] = save_results  # save the results

# plot the potentials
signals["plot_potential_mAB"] = lambda x: plot_potential(x,"mAB")
signals["plot_potential_mAF"] = lambda x: plot_potential(x,"mAF")
signals["plot_potential_fermi"] = lambda x: plot_potential(x,"fermi")
signals["plot_potential_Bx"] = lambda x: plot_potential(x,"Bx")
signals["plot_potential_By"] = lambda x: plot_potential(x,"By")
signals["plot_potential_Bz"] = lambda x: plot_potential(x,"Bz")
signals["plot_potential_rashba"] = lambda x: plot_potential(x,"rashba")
signals["plot_potential_kanemele"] = lambda x: plot_potential(x,"kanwmele")
signals["plot_potential_swave"] = lambda x: plot_potential(x,"swave")



class BulkApp(object):       
      def __init__(self):
            global tmppath # temporal path
            builder.add_from_file(xmlpath+"moire.xml")
            builder.connect_signals(signals)
            self.window = builder.get_object("bulk_window")
            folder = create_folder()
            tmppath = os.getcwd() # get the initial directory
            self.window.show()



if __name__ == "__main__":
  inipath = os.getcwd() # get the initial directory
  app = BulkApp()
  gtk.main()

