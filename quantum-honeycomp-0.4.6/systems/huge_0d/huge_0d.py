#!/usr/bin/python


import sys
if len(sys.argv)>1: # if input provided
  sys.path.append(sys.argv[1]+ "pysrc")
  xmlpath = sys.argv[1] + "systems/huge_0d/"  # path of the xml file
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


def get_vacancies():
  """Get the value of a certain variable"""
  name = "vacancies" # name of the object
  ats = builder.get_object(name).get_text()
  ats = ats.replace(","," ") # substitute comma by space
  ats = ats.split() # separte bu commas
  print "Remove atoms = ",ats
  ats = [int(float(a)) for a in ats] # convert to int
  return ats # return list



def get_numbers(name):
  """Get the values of a certain variable"""
  ats = builder.get_object(name).get_text()
  ats = ats.replace(","," ") # substitute comma by space
  ats = ats.split() # separte bu commas
  print name,"selected = ",ats
  ats = [int(float(a)) for a in ats] # convert to int
  return ats




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
  def rv():
    return .5 - np.random.random(3) # 3 component random vector
  def rvxy():
    v = np.random.random(3) -.5 # 3 component random vector
    v[0] = 0.0
    return v
  ###################################
  if scfin == "Reconverge": # no new mean field
    pass
  else: # if create new mean field matrix
    mf = interactions.initialization(g,scfin=scfin)  # create the mean field mat
    input_tb90.write_mean_field(mf,output_file="mean_field.in") # write in file
  # end SCF options
  ##################
  # end SCF options
  intb90.electrons.filling = get("filling")  # get the filling
  u = get("hubbard")
  intb90.scf_convergence.hubbard_scf = u # Hubbard
  if u > 1.5:
    intb90.scf_convergence.mix_coef = 0.1 
  intb90.write() # write input



def custom_tb90_bands():
  """Modifies the properties of tb90.in input"""
  intb90 = input_tb90.tb90in() # create input
  intb90.mode("nothing")  # no calculation
  intb90.mode("bands")   # bands calculation
  # type of operator for bands
  op = builder.get_object("bands_color").get_active_text()
  if op=="None": intb90.bands.bands_operators_option = "None"
  elif op=="Off-plane S": intb90.bands.bands_operators_option = "Sz"
  elif op=="In-plane S": intb90.bands.bands_operators_option = "Sx"
  elif op=="Position": intb90.bands.bands_operators_option = "index"
  else: intb90.bands.bands_operators_option = "from_file"
  intb90.bands.klist_bands = "from_file"
  intb90.bands.use_ewindow_bands = False
  # end type of operator
  intb90.write() # write input




def get_geometry0d():
  """ Create a 0d island"""
  lattice_name = builder.get_object("lattice").get_active_text()
  if lattice_name=="Honeycomb":
    geometry_builder = geometry.honeycomb_lattice
  if lattice_name=="Square":
    geometry_builder = geometry.square_lattice
  if lattice_name=="Kagome":
    geometry_builder = geometry.kagome_lattice
  if lattice_name=="Lieb":
    geometry_builder = geometry.lieb_lattice
  # first create a raw unit cell
  g = geometry_builder()  # build a 2d unit cell
  n = 1+int(round(get("size")))   # get the desired size
  nf = get("size")   # get the desired size, in float
  g = g.supercell(4*n)   # create supercell
  g.set_finite() # set as finite system
  g.center() # center the geometry
  # now scuplt the geometry
  nedges = int(get("nedges")) # number of edges
  g = sculpt.rotate(g,get("rotation")*2.*np.pi/360) # initial rotation
  def f(x,y): return x>-nf*(np.cos(np.pi/3)+1.)  # function to use as cut
  for i in range(nedges): # loop over rotations, 60 degrees
    g = sculpt.intersec(g,f) # retain certain atoms
    g = sculpt.rotate(g,2.*np.pi/nedges) # rotate 60 degrees
  g = sculpt.remove_unibonded(g)  # remove single bonded atoms
  g.center() # center the geometry
  #g = sculpt.remove_central(g,int(get("numvac"))) # removes several cen atoms
  g = sculpt.remove(g,get_vacancies()) # removes several atoms
  return g











def initialize(self):
  """ Initialize the calculation"""
  global g
  g = get_geometry0d() # get the geometry
  h = hamiltonians.hamiltonian(g) # get the hamiltonian
  h.has_spin = builder.get_object("is_spinful").get_active() # spin treatmen
  if h.has_spin: # spinful hamiltonian
    h.first_neighbors()  # first neighbor hoppin
    h.add_zeeman([get("Bx"),get("By"),get("Bz")]) # Zeeman fields
    h.add_sublattice_imbalance(get("mAB"))  # sublattice imbalance
    if abs(get("rashba")) > 0.0: h.add_rashba(get("rashba"))  # Rashba field
    h.add_antiferromagnetism(get("mAF"))  # AF order
    if abs(get("kanemele"))>0.0:  h.add_kane_mele(get("kanemele")) # intrinsic SOC
    h.shift_fermi(get("fermi")) # shift fermi energy
  else: # spinless treatment
    print "Spinless hamiltonian"
    h.is_sparse = True
    h.first_neighbors()  # first neighbor hopping
  h.add_peierls(get("peierls")) # add magnetic field
  h.write("hamiltonian.in")


  

def show_ldos(self):
  h = hamiltonians.hamiltonian() # create hamiltonian object
  h.read() # read hamiltonian from file
  h.has_spin = builder.get_object("is_spinful").get_active()
  print "Hamiltonian has been read"
#  points = 600 # number of polynomials
  points = int(get("LDOS_polynomials")) # number of polynomials
  x = np.linspace(-.9,.9,points*10) # energies
  h.intra = h.intra/6.0 # normalize
  atoms = get_numbers("LDOS_num_atom")
  os.system("rm LDOS_*")
  for iatom in atoms: # loop over atoms
#  iatom = get("LDOS_num_atom")
    if h.has_spin: iatom = iatom*2 # if spinful
    mus = kpm.local_dos(h.intra,n=points,i=iatom) # calculate moments
    y = kpm.generate_profile(mus,x) # calculate DOS 
    dos.write_dos(x*6.,y) # write dos in file
    fname = "LDOS_"+str(iatom)+".OUT" # name of the file
    os.system("cp DOS.OUT " + fname) # save the file
  os.system("qh-several-ldos &")


def show_dos(self):
  h = hamiltonians.hamiltonian() # create hamiltonian object
  h.read() # read hamiltonian from file
  h.has_spin = builder.get_object("is_spinful").get_active()
  print "Hamiltonian has been read"
  points = int(get("DOS_polynomials")) # number of polynomials
  x = np.linspace(-.9,.9,points*10) # energies
  h.intra = h.intra/6.0 # normalize
  ntries = get("DOS_iterations")
  mus = kpm.random_trace(h.intra,n=points,ntries=100) # calculate moments
  y = kpm.generate_profile(mus,x) # calculate DOS 
  dos.write_dos(x*6.,y) # write dos in file
  os.system("tb90-dos &")


def show_spatial_dos(self):
  h = hamiltonians.hamiltonian() # create hamiltonian object
  h.read() # read hamiltonian from file
  h.has_spin = builder.get_object("is_spinful").get_active()
  print "Hamiltonian has been read"
#  ldos.ldos0d_wf(h,e = get("energy_spatial_DOS"),
#                   delta = get("smearing_spatial_DOS"),
#                   num_wf = get("waves_spatial_DOS"))
  ldos.ldos0d(h,e = get("energy_spatial_DOS"),
                   delta = get("smearing_spatial_DOS"))
  os.system("qh-ldos &")



def show_lattice(self):
  """Show the lattice of the system"""
  g = get_geometry0d() # get the geometry
  g.write()
  print "Structure has been created"
#  os.system("qh-light-structure &")
  os.system("qh-pick &")



# create signals
signals = dict()
signals["on_window_destroy"] = gtk.main_quit  # close the window
signals["initialize"] = initialize  # initialize and run
signals["show_ldos"] = show_ldos  # show DOS
signals["show_dos"] = show_dos  # show DOS
signals["show_spatial_dos"] = show_spatial_dos  # show DOS
signals["show_lattice"] = show_lattice  # show magnetism

class IslandApp(object):       
	def __init__(self):
	    builder.add_from_file(xmlpath+"huge_0d.xml")
	    builder.connect_signals(signals)
	    self.window = builder.get_object("bulk_window")
            folder = create_folder()
	    self.window.show()



if __name__ == "__main__":
	app = IslandApp()
	gtk.main()

