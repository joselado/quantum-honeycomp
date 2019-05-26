from .qh_interface import execute_script
import numpy as np
from pygra import klist
from qh_interface import *
from pygra import parallel

def get_bands(h,window):
    """Compute the bandstructure of the system"""
    opname = window.getbox("bands_color")
    if opname=="None": op = None # no operators
    elif opname=="Sx": op = h.get_operator("sx") # off plane case
    elif opname=="Sy": op = h.get_operator("sy")# off plane case
    elif opname=="Sz": op = h.get_operator("sz")# off plane case
    elif opname=="Valley": op = h.get_operator("valley")
    elif opname=="IPR": op = h.get_operator("ipr")
    elif opname=="y-position": op = h.get_operator("yposition")
    elif opname=="x-position": op = h.get_operator("xposition")
    elif opname=="z-position": op = h.get_operator("zposition")
    elif opname=="Interface": op = h.get_operator("interface")
    elif opname=="Layer": op = h.get_operator("zposition")
    else: op = None
    kpath = klist.default(h.geometry,nk=int(window.get("nk_bands")))
    num_bands = int(window.get("nbands"))
    if num_bands<1: num_bands = None # all the eigenvalues
    check_parallel(window) # check if use parallelization
    h.get_bands(operator=op,kpath=kpath,num_bands=num_bands)
    command = "qh-bands --dim "+str(h.dimensionality) 
    if op is not None: command += " --cblabel "+opname
    if window.getbox("bands_colormap") is not None: 
        command += " --cmap "+window.getbox("bands_colormap")
    execute_script(command) # execute the command




def get_kdos(h,window):
    """Show the KDOS"""
    ew = window.get("kdos_ewindow")
    new = int(window.get("kdos_mesh")) # scale as kpoints
    energies = np.linspace(-ew,ew,new) # number of ene
    kpath = [[i,0.,0.] for i in np.linspace(0.,1.,new)]
    kdos.surface(h,energies=energies,delta=ew/new,kpath=kpath)
    command = "qh-kdos-both --input KDOS.OUT"
    try: command += " --cmap "+ window.getbox("kdos_colormap")
    except: pass
    execute_script(command) # execute the script







def check_parallel(qtwrap):
  """Check if there is parallelization"""
  if qtwrap.getbox("use_parallelization") =="Yes":
      parallel.cores = parallel.maxcpu
  else: parallel.cores = 1 # single core



def set_colormaps(form,name,cs=[]):
    """Add the different colormaps to a combox"""
    cb = getattr(form,name)
    try: cb = getattr(form,name)
    except:
        print("Combox",name,"not found")
        return
    cb.clear() # clear the items
    cb.addItems(cs)



def initialize(window):
    """Do various initializations"""
    cs = ["RGB","hot","inferno","plasma","bwr","rainbow","gnuplot"]
    set_colormaps(window.form,"bands_colormap",cs=cs) # set the bands
    cs = ["hot","inferno","plasma","rainbow","gnuplot","cool"]
    set_colormaps(window.form,"kdos_colormap",cs=cs) # set the bands



