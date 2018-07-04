#!/usr/bin/python2
from __future__ import print_function
import time


import sys
if len(sys.argv)>1: # if input provided
  sys.path.append(sys.argv[1]+ "pysrc")
  xmlpath = sys.argv[1] + "systems/huge_0d/"  # path of the xml file
else: # default path
  import os
  sys.path.append(os.getcwd()+"/../../pysrc")
  xmlpath = ""  # path of the xml file


import gi 
gi.require_version('Gtk', '3.0')
from gi.repository import Gtk as gtk
builder = gtk.Builder()

from qh_interface import * # import all the libraries needed


def getactive(name):
  return builder.get_object(name).get_active() # spin treatment



def modify(name,value,active=False):
  """Modify a certain text entry"""
  entry = builder.get_object(name)
  entry.set_text(str(value))
  if active: entry.activate() # set active



def get(name):
  """Get the value of a certain variable"""
  return float(builder.get_object(name).get_text())


def getfile(name):
  """Get the name of the file"""
  return builder.get_object(name).get_filename()



def get_vacancies():
  """Get the value of a certain variable"""
  name = "vacancies" # name of the object
  ats = builder.get_object(name).get_text()
  ats = ats.replace(","," ") # substitute comma by space
  ats = ats.split() # separte bu commas
#  print "Remove atoms = ",ats
  ats = [int(float(a)) for a in ats] # convert to int
  return ats # return list



def get_numbers(name):
  """Get the values of a certain variable"""
  ats = builder.get_object(name).get_text()
  ats = ats.replace(","," ") # substitute comma by space
  ats = ats.split() # separte bu commas
  ats = [int(float(a)) for a in ats] # convert to int
  return ats



def getbox(name):
  return builder.get_object(name).get_active_text()


def get_geometry0d(second_call=False):
  """ Create a 0d island"""
  t0 = time.clock() # initial time
  lattice_name = getbox("lattice")
  if lattice_name=="Honeycomb":
    geometry_builder = geometry.honeycomb_lattice
  if lattice_name=="Square":
    geometry_builder = geometry.square_lattice
  if lattice_name=="Kagome":
    geometry_builder = geometry.kagome_lattice
  if lattice_name=="Lieb":
    geometry_builder = geometry.lieb_lattice
  if lattice_name=="Triangular":
    geometry_builder = geometry.triangular_lattice
  # first create a raw unit cell
  gbulk = geometry_builder()  # build a 2d unit cell
  if builder.get_object("activate_bilayer").get_active(): 
    gbulk = multilayers.bilayer_geometry(gbulk) 
  # now scuplt the geometry
  nf = 1+get("size")   # get the desired size, in float
  ####################################
  ####################################
  if getbox("geometry_mode") == "Positions": # generate a perfect island
    os.system("cp "+getfile("positions_file")+" POSITIONS.OUT")
    g = geometry.read()
    g.center()
    return g
  if getbox("geometry_mode") == "Recipe": # generate a perfect island
    nedges = int(get("nedges")) # number of edges
    angle = get("rotation")*2.*np.pi/360 # angle to rotate
    g = islands.get_geometry(geo=gbulk,n=nf,nedges=nedges,
                               rot=angle,clean=False)
  elif getbox("geometry_mode") == "Image": # generate from an image
    print("Direction",getfile("image_path"))
    g = sculpt.image2island(getfile("image_path"),gbulk,size=int(nf),color="black")
  ####################################
  #############################################
# if a precise diameter is wanted, only for the first call
  if getactive("target_diameter") and not second_call: 
    diameter = get("desired_diameter")
    ratio = diameter/g.get_diameter() # ratio between wanted and obtained
    print("\nChecking that it has the desired size",ratio)
    if not 0.99<ratio<1.01: # if outside the tolerance
      newsize = round(ratio*float(get("size"))) # new size
      modify("size",newsize,active=True) # modify the value
      print("Recalling the geometry with size",newsize)
      return get_geometry0d(second_call = True)
  else: pass
  ####################################
  # clean the island
  g.center() # center the geometry
  ############################################3
  g = modify_geometry(g) # modify the geometry in several ways
  print("Total number of atoms =",len(g.r))
  print("Time spent in creating the geometry =",time.clock() - t0)
  if builder.get_object("clean_island").get_active(): # if it is cleaned
    g = sculpt.remove_unibonded(g,iterative=True)  # remove single bonded atoms
  return g



def modify_geometry(g):
  """Modify the geometry according to the interface"""
  mtype = getbox("modify_geometry")
  print("Modifying geometry according to",mtype)
  if mtype == "None": return g # do nothing
  elif mtype == "Index": 
    return sculpt.remove(g,get_vacancies()) # removes several atoms
  elif mtype=="Choose atoms": # special case
    print("Removing as chosen\n")
    try:
      inds = np.genfromtxt("REMOVE_ATOMS.INFO") # selected atoms
      print("Removed indexes",inds)
    except: return g
    try:  
      inds = [int(i) for i in inds] # as integer
    except: inds = [int(inds)]
    try: return sculpt.remove(g,inds) # removes several atoms
    except: return g









def initialize(self):
  """ Initialize the calculation"""
  t0 = time.clock()
  os.system("rm SELECTED_ATOMS.INFO") # remove the file for the DOS
  g = get_geometry0d() # get the geometry
  h = hamiltonians.hamiltonian(g) # get the hamiltonian
  h.has_spin = builder.get_object("is_spinful").get_active() # spin treatment
  if h.has_spin: # spinful hamiltonian
    print("Spinful hamiltonian, DO NOT USE VERY LARGE ISLANDS!!!")
    h.is_sparse = False
    h.first_neighbors()  # first neighbor hoppin
    h.add_zeeman([get("Bx"),get("By"),get("Bz")]) # Zeeman fields
    if abs(get("rashba")) > 0.0: h.add_rashba(get("rashba"))  # Rashba field
    h.add_antiferromagnetism(get("mAF"))  # AF order
    if abs(get("kanemele"))>0.0:  h.add_kane_mele(get("kanemele")) # intrinsic SOC
    h.shift_fermi(get("fermi")) # shift fermi energy
    h.turn_sparse() # turn it sparse
  else: # spinless treatment
    h.is_sparse = True
    print("Spinless hamiltonian")
    if builder.get_object("activate_bilayer").get_active(): 
      tlayer = get("interlayer_coupling") # get the interlayer coupling
      efield = get("efield") # get the interlayer coupling
      def f(r1,r2):
        dr=r1-r2  # distance
        dr = dr.dot(dr) # absolute distance
        if dr<0.0001: return efield*np.sign(r1[2])# electric field
        if np.abs(r1[2]-r2[2])<0.0001: # same layer
          if 0.9<dr<1.1: return 1.0 # first neighbor
        else:
          if 3.9<dr<4.1: return tlayer # first neighbor
        
        return 0.0
      hamiltonians.generate_parametric_hopping(h,f) # hamiltonian
    else: # conventional first neighbors 
      h.first_neighbors()  # first neighbor hopping
  h.add_sublattice_imbalance(get("mAB"))  # sublattice imbalance
  h.add_peierls(get("peierls")) # add magnetic field
  modh = getbox("modify_hamiltonian") # way of modifying the Hamiltonian
  edgesites = edge_atoms(h.geometry) # get the edge atoms
  if modh == "None":
    print("No modifications to the Hamiltonian")
  elif modh == "Edge potential": 
    print("Added an edge potential")
    if not h.has_spin:
      from scipy.sparse import diags
      h.intra = h.intra + get("potval")*diags([edgesites],[0]) # Add
    else: raise
  elif modh == "Disorder": 
    print("Added Anderson disorder")
    adis = get("anderson_disorder") # value of the disorder
    def rfun(x,y,z=0): return (np.random.random()-0.5)*adis
    h.shift_fermi(rfun) # random onsites
  # part for bilayer systems
  #####
  print("Time spent in creating the Hamiltonian =",time.clock() - t0)
  h.geometry.write()
  h.save() # save the Hamiltonian


  

def show_ldos(self):
  h = load_hamiltonian() # get the hmiltonian
  points = int(get("LDOS_polynomials")) # number of polynomials
  x = np.linspace(-.9,.9,int(get("num_ene_ldos"))) # energies
  h.intra = h.intra/6.0 # normalize
  # make the update of the atoms number if SELECTED_ATOMS file exist
  if os.path.isfile("SELECTED_ATOMS.INFO"): 
    ind_atoms = open("SELECTED_ATOMS.INFO").read().replace("\n",", ")
    modify("LDOS_num_atom",ind_atoms,active=True)

  # now continue in the usual way
  atoms = get_numbers("LDOS_num_atom")
  ecut = get("energy_cutoff_local_dos")/6.
  os.system("rm LDOS_*") # remove
  for iatom in atoms: # loop over atoms
    if h.has_spin: iatom = iatom*2 # if spinful
    x = np.linspace(-ecut,ecut,int(get("num_ene_ldos")),endpoint=True) # energies
    mus = kpm.local_dos(h.intra,n=points,i=iatom) # calculate moments
    y = kpm.generate_profile(mus,x) # calculate DOS 
    x,y = x*6.,y/6. # renormalize
    y = dos.convolve(x,y,delta=get("smearing_local_dos")) # add broadening
    dos.write_dos(x,y) # write dos in file
    fname = "LDOS_"+str(iatom)+".OUT" # name of the file
    os.system("mv DOS.OUT " + fname) # save the file
  execute_script("qh-several-ldos ")



def show_full_spectrum(self):
  h = load_hamiltonian() # get the hmiltonian
  nmax = 10000
  if h.intra.shape[0]<nmax:
    h.get_bands()    
    execute_script("qh-bands0d ")
  else:
    print("Too large Hamiltonian ",nmax)



def load_hamiltonian():
  h = hamiltonians.load() # load the Hamiltonian
  return h



def show_dos(self):
  h = load_hamiltonian() # get the hmiltonian
  h.has_spin = builder.get_object("is_spinful").get_active()
  points = int(get("DOS_polynomials")) # number of polynomials
  x = np.linspace(-.9,.9,int(get("num_ene_dos"))) # energies
  h.intra = h.intra/6.0 # normalize
  ntries = get("DOS_iterations")
  t0 = time.clock() # initial time
  mus = kpm.random_trace(h.intra,n=points,ntries=100) # calculate moments
  y = kpm.generate_profile(mus,x) # calculate DOS 
  x,y = x*6.,y/6. # renormalize
  y = dos.convolve(x,y,delta=get("smearing_dos")) # add broadening
  dos.write_dos(x,y) # write dos in file
  execute_script("tb90-dos ")
  print("Time spent in Kernel PM DOS calculation =",time.clock() - t0)


def show_spatial_dos(self):
  t0 = time.clock()
  h = load_hamiltonian() # get the hamiltonian
  mode_stm = getbox("mode_stm") # get the way the images will be calculated
  delta = get("smearing_spatial_DOS")
  def shot_dos(energy):
    if mode_stm == "Full": # matrix inversion
      ldos.ldos0d(h,e = energy,
                     delta = delta)
    elif mode_stm == "Eigen": # Using Arnoldi
      ldos.ldos0d_wf(h,e = energy,
                     delta = delta,
                     num_wf = int(get("nwaves_dos")),
                     robust=True,tol=delta/1000)
    if builder.get_object("show_upper").get_active(): 
      execute_script("qh-upper-layer LDOS.OUT")
  mode_dosmap = getbox("mode_dosmap") # one shot or video
  if mode_dosmap=="Single shot":
    energy = get("energy_spatial_DOS") # one shot
    shot_dos(energy) # calculate
    print("Time spent in STM calculation =",time.clock() - t0)
    if builder.get_object("3d_stm").get_active(): # how to plot
      execute_script("qh-ldos ") # usig mayavi
    else: 
      execute_script("qh-fast-ldos LDOS.OUT  ") # using matplotlib
  if mode_dosmap=="Movie": # do a sweep
    energies = np.linspace(get("mine_movie"),get("maxe_movie"),
                               int(get("stepse_movie")))
    fof = open("FRAMES.OUT","w") # file with frames
    for energy in energies:
      print("Calculating",energy)
      shot_dos(energy) # calculate
      name = "ENERGY_"+'{0:.8f}'.format(energy)+"_LDOS.OUT"
      os.system("mv LDOS.OUT "+name) # save the data
      execute_script("qh-silent-ldos "+name) # save the image
      namepng = name+".png" # name of the image
      fof.write(namepng+"\n") # save the name
    fof.close() # close file
    os.system("xdg-open "+namepng) # open the last file

def edge_atoms(g,nn=3):
  """Get the edge potential"""
  cs = g.get_connections() # get the connections
  v1 = np.array([int(len(c)<nn) for c in cs]) # check if the atom is on the edge or not
  # and the first neighbors to the edge
  v2 = np.zeros(len(v1)) # initialize
  for i in range(len(v1)): # loop
    if v1[i]==0: # not in the edge yet
      for ic in cs[i]: # loop over neighbors
        if v1[ic]==1: # edge atom
          v2[i] = 1 # assign
          break
  v = v1 + v2*2./3. # sum
  np.savetxt("EDGE.OUT",np.matrix([g.x,g.y,v]).T) # save
  return v # return the array


def show_potential(self):
  execute_script("qh-absolute-potential EDGE.OUT  ")




def show_lattice(self):
  """Show the lattice of the system"""
  g = get_geometry0d() # get the geometry
  g.write()
  print("Structure has been created")
  if getactive("show3d"): execute_script("qh-pick ")
  else: execute_script("qh-fast-pick  ")


def show_path_dos(self):
  """Show the the DOS in the path"""
  calculate_path_dos() # calculate the path
  h = load_hamiltonian() # get the hmiltonian
  pols = int(get("pols_path")) # number of polynomials
  h.intra = h.intra/6.0 # normalize
  atoms = np.genfromtxt("PATH.OUT").transpose()[0] # indexes
  atoms = [int(a) for a in atoms] # indexes of the atoms
  os.system("rm LDOS_*") # clean all the data
  ecut = np.abs(get("ecut_path")) # maximum energy in path
  if ecut>5.5: ecut = 5.9
  ecut /= 6.0 # to interval 0,1
  for iatom in atoms: # loop over atoms
    print("Calculating DOS in ",iatom)
    if h.has_spin: iatom = iatom*2 # if spinful
    mus = kpm.local_dos(h.intra,n=pols,i=iatom) # calculate moments
    x = np.linspace(-ecut,ecut,int(get("num_ene_path"))) # energies
    y = kpm.generate_profile(mus,x) # calculate DOS 
    xout,yout = x*6.,y/6. #renormalize
    yout = dos.convolve(xout,yout,delta=get("smearing_path_dos")) # add broadening
    dos.write_dos(xout,yout) # write dos in file
    fname = "LDOS_"+str(iatom)+".OUT" # name of the file
    os.system("cp DOS.OUT " + fname) # save the file
  execute_script("qh-dos-path  ")











def show_path(self):
  """Show the path followed in the DOS, when calculating several
  atoms"""
  calculate_path_dos() # calculate the path
  execute_script("qh-path  ")




def calculate_path_dos():
  """Calculate all the DOS in a path"""
  i0 = int(get("initial_atom"))
  i1 = int(get("final_atom"))
  h = load_hamiltonian() # get the hamiltonian
  r0 = h.geometry.r[i0] # initial point
  r1 = h.geometry.r[i1] # final point
  # now do something dangerous, rotate the geometry to check which
  # atoms to accept
  print("Initial position",r0)
  print("Final position",r1)
  print("Created new rotated geometry")
  gtmp = h.geometry.copy() # copy geometry
  gtmp.shift(r0) # shift the geometry
  gtmp = sculpt.rotate_a2b(gtmp,r1-r0,np.array([1.,0.,0.]))
  # new positions
  r0 = gtmp.r[i0]
  r1 = gtmp.r[i1]
  dr = r1 - r0 # vector between both points
  dy = np.abs(get("width_path")) # width accepted
  dx = np.sqrt(dr.dot(dr))+0.1
  print("Initial ROTATED position",r0)
  print("Final ROTATED position",r1)
  # and every coordinate
  x1,y1 = r0[0],r0[1]
  x2,y2 = r1[0],r1[1]
  ym = (y1+y2)/2. # average y
  def fun(r):
    """Function that decides which atoms to calculate"""
    x0 = r[0]
    y0 = r[1]
    if (x1-dy)<x0<(x2+dy) and np.abs(y0-ym)<dy: return True
    else: return False
  inds = sculpt.intersected_indexes(gtmp,fun) # indexes of the atoms
  fo = open("PATH.OUT","w") # output file
  fo.write("# index of the atom, step in the path, x, y\n")
  ur = dr/np.sqrt(dr.dot(dr)) # unitary vector
  steps = [(gtmp.r[i] - r0).dot(ur) for i in inds] # proyect along that line
  inds = [i for (s,i) in sorted(zip(steps,inds))] # sort by step
  steps = sorted(steps) # sort the steps
  # only select those between the initial and final points
  inds0,steps0 = [],[]
#  for (i,s) in (inds,steps):
#    if 0<s<1.
  ####################
  g = h.geometry
  for (i,s) in zip(inds,steps):
    fo.write(str(i)+"    "+str(s)+"   "+str(g.x[i])+"    "+str(g.y[i])+"\n")
  fo.close() # close file



save_results = lambda x: save_outputs(inipath,tmppath) # function to save


def clear_removal(self):
  os.system("rm SELECTED.INFO") # remove the file


def select_atoms(self):
  execute_script("qh-fast-pick write remove  ") # remove the file


def select_atoms_dos(self):
  execute_script("qh-fast-pick write  ") # remove the file


# create signals
signals = dict()
signals["on_window_destroy"] = gtk.main_quit  # close the window
signals["initialize"] = initialize  # initialize and run
signals["show_ldos"] = show_ldos  # show LDOS
signals["show_dos"] = show_dos  # show DOS
signals["show_spatial_dos"] = show_spatial_dos  # show DOS
signals["show_lattice"] = show_lattice  # show magnetism
signals["show_full_spectrum"] = show_full_spectrum  # show all the eigenvalues
signals["show_path"] = show_path  # show the path
signals["show_path_dos"] = show_path_dos  # show the path
signals["show_potential"] = show_potential  # show the potential added
signals["save_results"] = save_results  # save the results
signals["clear_removal"] = clear_removal  # clear the file
signals["select_atoms"] = select_atoms  # select_atoms
signals["select_atoms_dos"] = select_atoms_dos  # select_atoms

class IslandApp(object):       
     def __init__(self):
            global tmppath # temporal path
            builder.add_from_file(xmlpath+"huge_0d.xml")
            builder.connect_signals(signals)
            self.window = builder.get_object("bulk_window")
            folder = create_folder()
            tmppath = os.getcwd() # get the initial directory
            print("Temporal path is",tmppath)
            self.window.show()



if __name__ == "__main__":
        inipath = os.getcwd() # get the initial directory
        print("Initial path is",inipath)
        app = IslandApp()
        gtk.main()

