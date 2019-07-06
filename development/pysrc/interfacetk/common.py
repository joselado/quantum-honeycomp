from .qh_interface import execute_script
import numpy as np
from pygra import klist
from qh_interface import *
from pygra import parallel

def get_operator(h,opname,projector=False):
    """Return an operator"""
    if opname=="None": op = None # no operators
    elif opname=="Sx": op = h.get_operator("sx") # off plane case
    elif opname=="Sy": op = h.get_operator("sy")# off plane case
    elif opname=="Sz": op = h.get_operator("sz")# off plane case
    elif opname=="Valley": op = h.get_operator("valley",projector=projector)
    elif opname=="IPR": op = h.get_operator("ipr")
    elif opname=="y-position": op = h.get_operator("yposition")
    elif opname=="x-position": op = h.get_operator("xposition")
    elif opname=="z-position": op = h.get_operator("zposition")
    elif opname=="Interface": op = h.get_operator("interface")
    elif opname=="Surface": op = h.get_operator("surface")
    elif opname=="Layer": op = h.get_operator("zposition")
    else: op = None
    return op



def get_bands(h,window):
    """Compute the bandstructure of the system"""
    opname = window.getbox("bands_color")
    op = get_operator(h,opname) # get operator
    kpath = klist.default(h.geometry,nk=int(window.get("nk_bands")))
    num_bands = int(window.get("nbands"))
    if num_bands<1: num_bands = None # all the eigenvalues
    check_parallel(window) # check if use parallelization
    h.get_bands(operator=op,kpath=kpath,num_bands=num_bands)
    command = "qh-bands --dim "+str(h.dimensionality) 
    if op is not None: command += " --cblabel "+opname
#    if window.getbox("bands_colormap") is not None: 
#        command += " --cmap "+window.getbox("bands_colormap")
    execute_script(command) # execute the command




def get_kdos(h,window):
    """Show the KDOS"""
    ew = window.get("kdos_ewindow")
    new = int(window.get("kdos_mesh")) # scale as kpoints
    energies = np.linspace(-ew,ew,new) # number of ene
    kpath = [[i,0.,0.] for i in np.linspace(0.,1.,new)]
    kdos.surface(h,energies=energies,delta=ew/new,kpath=kpath)
    command = "qh-kdos-both --input KDOS.OUT"
    execute_script(command) # execute the script



def show_exchange(h,window):
    """Show the exchange field"""
    nrep = max([int(window.get("magnetization_nrep")),1]) # replicas
    h.write_magnetization(nrep=nrep) # write the magnetism
    execute_script("qh-moments",mayavi=True)


def get_dos(h,window,silent=False):
  nk = max([int(window.get("dos_nk")),1])
  delta = window.get("dos_delta")
  ewindow = abs(window.get("dos_ewindow"))
  energies = np.linspace(-ewindow,ewindow,ewindow/delta*5) # get the energies
  if window.getbox("dos_mode")=="Green":
    dos.dos(h,delta=delta,nk=nk,energies=energies,mode="Green") # compute DOS
  else:
    dos.dos(h,delta=delta,nk=nk,energies=energies) # compute DOS
  if not silent: execute_script("qh-dos --input DOS.OUT")



def get_berry1d(h,window):
    """Get the one dimensional Berry curvature"""
    ks = klist.default(h.geometry,
            nk=int(window.get("topology_nk")))  # write klist
    opname = window.getbox("topology_operator")
    op = get_operator(h,opname,projector=True) # get operator
    topology.write_berry(h,ks,operator=op)
    command = "qh-berry1d  --label True " 
    if opname!="None": command += " --mode "+opname
    execute_script(command)





def get_berry2d(h,window):
    """Get the Berry curvature"""
    nk = int(np.sqrt(window.get("topology_nk")))
    opname = window.getbox("topology_operator")
    op = get_operator(h,opname,projector=True) # get operator
    topology.berry_map(h,nk=nk,operator=op)
    execute_script("qh-berry2d BERRY_MAP.OUT")


def get_chern(h,window):
    """Get the Chern number"""
    nk = int(np.sqrt(window.get("topology_nk")))
    opname = window.getbox("topology_operator")
    op = get_operator(h,opname,projector=True) # get operator
    topology.chern(h,nk=nk,operator=op)
    execute_script("qh-chern BERRY_CURVATURE.OUT")

def get_fermi_surface(h,window):
    check_parallel(window) # check if use parallelization
    e = window.get("fs_ewindow")
    energies = np.linspace(-e,e,100)
    nk = int(window.get("fs_nk")) # number of kpoints
    numw = int(window.get("fs_numw")) # number of waves for sparse
    delta = window.get("fs_delta")
    spectrum.multi_fermi_surface(h,nk=nk,energies=energies,
        delta=delta,nsuper=1,numw=numw)
    execute_script("qh-multifermisurface")





def get_z2(h,window):
    nk = int(np.sqrt(window.get("topology_nk")))
    topology.z2_vanderbilt(h,nk=nk,nt=nk/2) # calculate z2 invariant
    execute_script("qh-wannier-center  ") # plot the result



def get_multildos(h,window):
    check_parallel(window) # check if use parallelization
    ewin = window.get("multildos_ewindow")
    nrep = int(max([1,window.get("multildos_nrep")]))
    nk = int(max([1,window.get("multildos_nk")]))
    numw = int(window.get("multildos_numw"))
    ne = 100 # 100 points
    delta = window.get("multildos_delta")
    ldos.multi_ldos(h,es=np.linspace(-ewin,ewin,ne),
            nk=nk,delta=delta,nrep=nrep,numw=numw)
    execute_script("qh-multildos ")







def check_parallel(qtwrap):
  """Check if there is parallelization"""
  if qtwrap.getbox("use_parallelization") =="Yes":
      parallel.cores = parallel.maxcpu
  else: parallel.cores = 1 # single core



def set_colormaps(form,name,cs=[]):
    """Add the different colormaps to a combox"""
    try: cb = getattr(form,name)
    except:
        print("Combobox",name,"not found")
        return
    cb.clear() # clear the items
    cb.addItems(cs)



def initialize(window):
    """Do various initializations"""
    cs = ["RGB","hot","inferno","plasma","bwr","rainbow","gnuplot"]
    set_colormaps(window.form,"bands_colormap",cs=cs) # set the bands



