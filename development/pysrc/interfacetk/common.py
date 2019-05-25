from .qh_interface import execute_script
from pygra import klist
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



def check_parallel(qtwrap):
  """Check if there is parallelization"""
  if qtwrap.getbox("use_parallelization") =="Yes":
      parallel.cores = parallel.maxcpu
  else: parallel.cores = 1 # single core



def set_colormaps(form,name):
    """Add the different colormaps to a combox"""
    cb = getattr(form,name)
    try: cb = getattr(form,name)
    except:
        print("Combox",name,"not found")
        return
    cb.clear() # clear the items
    cs = ["RGB","hot","inferno","plasma","bwr","rainbow","gnuplot"]
    cb.addItems(cs)


def initialize(window):
    """Do various initializations"""
    set_colormaps(window.form,"bands_colormap")



