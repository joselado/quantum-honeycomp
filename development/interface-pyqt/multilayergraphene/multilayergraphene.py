#!/usr/bin/python

from __future__ import print_function

import sys
import os

qhroot = os.environ["QHROOT"] # root path
sys.path.append(qhroot+"/pysrc/") # python libraries


import qtwrap # import the library with simple wrappaers to qt4
get = qtwrap.get  # get the value of a certain variable
getbox = qtwrap.getbox  # get the value of a certain variable
window = qtwrap.main() # this is the main interface



from qh_interface import * # import all the libraries needed








def get_geometry(modify=True):
  """ Create a 0d island"""
  lattice = getbox("lattice") # get the option
  ns = [] # empt list
  for il in lattice:
      if il=="A": ns += [-1]
      elif il=="B": ns += [0]
      elif il=="C": ns += [1]
      else: print("Error",il)
  print(ns)
  g = specialgeometry.multilayer_graphene(l=ns)
  nsuper = int(get("nsuper"))
  g = g.supercell(nsuper)
  if modify: g = modify_geometry(g) # modify the geometry
  return g





def select_atoms_removal(self):
  g = get_geometry(modify=False) # get the unmodified geometry
  g.write() # write geometry
  execute_script("qh-remove-atoms-geometry-3d") # remove the file


def modify_geometry(g):
  """Modify the geometry according to the interface"""
  if qtwrap.is_checked("remove_selected"): # remove some atoms
      try:
        inds = np.array(np.genfromtxt("REMOVE_ATOMS.INFO",dtype=np.int))
        if inds.shape==(): inds = [inds]
      except: inds = [] # Nothing
      print(inds)
      g = sculpt.remove(g,inds) # remove those atoms
  if qtwrap.is_checked("remove_single_bonded"): # remove single bonds
      g = sculpt.remove_unibonded(g,iterative=True)
  return g # return geometry











def initialize():
  """ Initialize the calculation"""
  g = get_geometry() # get the geometry
  ti = get("interlayer")
  h = g.get_hamiltonian(has_spin=True,fun=specialhopping.multilayer(ti=ti))
  h.add_zeeman([get("Bx"),get("By"),get("Bz")]) # Zeeman fields
  h.add_sublattice_imbalance(get("mAB"))  # sublattice imbalance
  if abs(get("rashba")) > 0.0: h.add_rashba(get("rashba"))  # Rashba field
  h.add_antiferromagnetism(get("mAF"))  # AF order
  h.shift_fermi(get("fermi")) # shift fermi energy
  h.shift_fermi(lambda r: get("bias")*r[2]) # interlayer bias
  if abs(get("kanemele"))>0.0:  h.add_kane_mele(get("kanemele")) # intrinsic SOC
  if abs(get("haldane"))>0.0:  h.add_haldane(get("haldane")) # intrinsic SOC
  if abs(get("antihaldane"))>0.0:  h.add_antihaldane(get("antihaldane")) 
  if abs(get("antikanemele"))>0.0:  h.add_anti_kane_mele(get("antikanemele")) 
  if abs(get("swave"))>0.0:  h.add_swave(get("swave")) 
#  h.add_peierls(get("peierls")) # shift fermi energy

  return h


def show_bands(self=0):
  comp = computing() # create the computing window
  h = pickup_hamiltonian() # get hamiltonian
  opname = getbox("bands_color")
  if opname=="None": op = None # no operators
  elif opname=="Sx": op = h.get_operator("sx") # off plane case
  elif opname=="Sy": op = h.get_operator("sy")# off plane case
  elif opname=="Sz": op = h.get_operator("sz")# off plane case
  elif opname=="Valley": op = h.get_operator("valley")
  else: op =None
  kpath = klist.default(h.geometry,nk=int(get("nk_bands")))
  h.get_bands(operator=op,kpath=kpath)
  execute_script("qh-bands2d  ")
  comp.kill()



def show_dosbands(self=0):
  h = pickup_hamiltonian() # get hamiltonian
  kdos.kdos_bands(h,scale=get("scale_kbands"),ewindow=get("window_kbands"),
                   ne=int(get("ne_kbands")),delta=get("delta_kbands"),
                   ntries=int(get("nv_kbands")))
  execute_script("qh-dosbands  KDOS_BANDS.OUT ")


def show_fermi_surface(silent=False):
  h = pickup_hamiltonian() # get hamiltonian
  ndos = int(get("ne_dos"))
  if h.dimensionality==2:
    spectrum.fermi_surface(h,e=get("energy_fs"),nk=int(get("nk_fs")),
            nsuper = 2,reciprocal=True,delta=get("delta_fs"))
#    dos.dos2d(h,ndos=500,delta=get("delta_dos"),nk=int(get("nk_dos")),
#            window=get("window_dos"))
  else: raise
  if not silent: 
      execute_script("qh-fermi-surface FERMI_MAP.OUT") # show the result




def show_dos(silent=False):
  h = pickup_hamiltonian() # get hamiltonian
  ndos = int(get("ne_dos"))
  if h.dimensionality==2:
    dos.dos2d(h,ndos=500,delta=get("delta_dos"),nk=int(get("nk_dos")),
            window=get("window_dos"))
  else: raise
  if not silent: execute_script("qh-dos DOS.OUT") # show the result


def pickup_hamiltonian():
  if qtwrap.is_checked("do_scf"):
    return hamiltonians.load() # load the Hamiltonian
  else: # generate from scratch
    return initialize()












def show_berry2d():
  h = pickup_hamiltonian() # get hamiltonian
  nk = int(np.sqrt(get("nk_topology")))
  opname = getbox("operator_topology")
  if opname=="None": op=None
  elif opname=="Valley": op = operators.get_valley(h,projector=True)
  else: raise 
  topology.berry_map(h,nk=nk,operator=op)
  execute_script("qh-berry2d BERRY_MAP.OUT")

  


def show_structure(self):
  """Show the lattice of the system"""
  g = get_geometry() # get the geometry
  nsuper = int(get("nsuper_struct"))
  g = g.supercell(nsuper)
  g.write()
#  execute_script("qh-light-structure POSITIONS.OUT")
  execute_script("qh-structure-bond POSITIONS.OUT")
#  execute_script("qh-structure  ")



def show_kdos(self):
  h = pickup_hamiltonian()  # get the hamiltonian
  ew = get("ewindow_kdos")
  new = int(get("mesh_kdos")) # scale as kpoints
  energies = np.linspace(-ew,ew,new) # number of ene
  kpath = [[i,0.,0.] for i in np.linspace(0.,1.,new)]
  kdos.surface(h,energies=energies,delta=ew/new,kpath=kpath)
  execute_script("qh-kdos-both KDOS.OUT  ")



def show_berry1d(self):
  h = pickup_hamiltonian()  # get the hamiltonian
  ks = klist.default(h.geometry,nk=int(get("nk_topology")))  # write klist
  opname = getbox("operator_topology")
  if opname=="None": op=None
  elif opname=="Valley": op = operators.get_valley(h,projector=True)
  else: raise 
  topology.write_berry(h,ks,operator=op)
  execute_script("qh-berry1d  label  ")


def show_z2(self):
  h = pickup_hamiltonian()  # get the hamiltonian
  nk = get("nk_topology")
  topology.z2_vanderbilt(h,nk=nk,nt=nk/2) # calculate z2 invariant
  execute_script("qh-wannier-center  ") # plot the result



def solve_scf():
  """Perform a selfconsistent calculation"""
  comp = computing() # create the computing window
  scfin = getbox("scf_initialization")
  h = initialize() # initialize the Hamiltonian
  mf = scftypes.guess(h,mode=scfin)
  nk = int(get("nk_scf"))
  U = get("hubbard")
  filling = get("filling_scf")
  filling = filling%1.
  scf = scftypes.selfconsistency(h,nkp=nk,filling=filling,g=U,
                mf=mf,mode="U",smearing=get("smearing_scf"),
                mix = get("mix_scf"))
  scf.hamiltonian.save() # save in a file
  comp.kill()



def show_magnetism():
  """Show the magnetism of the system"""
  h = pickup_hamiltonian() # get the Hamiltonian
  h.write_magnetization() # write the magnetism
  execute_script("qh-moments",mayavi=True)




def show_structure_3d(self):
  """Show the lattice of the system"""
  g = get_geometry() # get the geometry
  nsuper = int(get("nsuper_struct"))
  g = g.supercell(nsuper)
  g.write()
  execute_script("qh-structure3d POSITIONS.OUT")



def show_interactive_ldos():
  comp = computing() # create the computing window
  h = pickup_hamiltonian()  # get the hamiltonian
  ewin = get("window_ldos")
  nrep = int(get("nsuper_ldos"))
  nk = int(get("nk_ldos"))
  ne = int(get("ne_ldos"))
  delta = get("delta_ldos")
  ldos.multi_ldos(h,es=np.linspace(-ewin,ewin,ne),nk=nk,delta=delta,
          nrep=nrep)
  comp.kill()
  execute_script("qh-multildos ")



def sweep_parameter():
    """Perform a sweep in a parameter"""
    pname = getbox("sweep_parameter") # get the parameter
    ps = np.linspace(get("sweep_initial"),get("sweep_final"),
               int(get("sweep_steps"))) # parameters
    def modify(p): # function to change the parameter
        if pname=="Sublattice imbalance": qtwrap.modify("mAB",p)
        elif pname=="Kane-Mele": qtwrap.modify("kanemele",p)
        elif pname=="Jx": qtwrap.modify("Bx",p)
        elif pname=="Jy": qtwrap.modify("By",p)
        elif pname=="Jz": qtwrap.modify("Bz",p)
        elif pname=="Rashba": qtwrap.modify("rashba",p)
        elif pname=="Haldane": qtwrap.modify("haldane",p)
        elif pname=="Anti-Haldane": qtwrap.modify("antihaldane",p)
        elif pname=="s-wave pairing": qtwrap.modify("swave",p)
        elif pname=="Fermi": qtwrap.modify("fermi",p)
        else: raise # not implemented
    cname = getbox("sweep_task") # type of computation
    out = [] # empty list
    for p in ps: # loop over the values
        modify(p) # modify the Hamiltonian
        h = pickup_hamiltonian() # get hamiltonian
        if cname=="DOS": 
            show_dos(silent=True) # compute DOS
            m = np.genfromtxt("DOS.OUT").transpose() # get dos
            es,ds = m[0],m[1]
            for (e,d) in zip(es,ds): # loop over energies
                out.append([p,e,d])
        elif cname=="Indirect gap": # compute the gap
            g = h.get_gap() # compute Gap
            out.append([p,g]) # store result
        elif cname=="Chern number": # compute the gap
            c = topology.chern(h,nk=int(np.sqrt(get("nk_topology"))))
            out.append([p,c]) # store result
        elif cname=="Eigenvalues": # compute the gap
            kpath = [np.random.random(3) for i in range(int(get("nk_bands")))]
            (ks,es) = h.get_bands() # compute eigenvalues
            for e in es: out.append([p,e]) # store
        else: raise
    np.savetxt("SWEEP.OUT",np.matrix(out)) # store result
    if cname=="DOS":
        execute_script("qh-sweep-dos SWEEP.OUT") # remove the file
    elif cname=="Indirect gap":
        execute_script("qh-indirect-gap SWEEP.OUT") # remove the file
    elif cname=="Chern number":
        execute_script("qh-chern-evolution SWEEP.OUT") # remove the file
    else:
        execute_script("qh-indirect-gap SWEEP.OUT") # remove the file
    






save_results = lambda x: save_outputs(inipath,tmppath) # function to save


# create signals
signals = dict()
signals["solve_scf"] = solve_scf  # initialize and run
signals["show_bands"] = show_bands  # show bandstructure
signals["show_structure"] = show_structure  # show bandstructure
signals["show_dos"] = show_dos  # show DOS
signals["show_berry2d"] = show_berry2d  # show DOS
signals["show_berry1d"] = show_berry1d  # show DOS
signals["show_kdos"] = show_kdos  # show DOS
signals["show_dosbands"] = show_dosbands  # show DOS
signals["show_z2"] = show_z2  # show DOS
signals["show_magnetism"] = show_magnetism  # show magnetism
signals["compute_sweep"] = sweep_parameter  
signals["show_structure_3d"] = show_structure_3d
signals["select_atoms_removal"] = select_atoms_removal
signals["show_interactive_ldos"] = show_interactive_ldos
signals["show_fermi_surface"] = show_fermi_surface






window.connect_clicks(signals)
folder = create_folder()
tmppath = os.getcwd() # get the initial directory
window.run()

