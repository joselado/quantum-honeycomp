from __future__ import print_function
import subprocess
import os
import sys
import numpy as np
# import the different libraries for quantum honeycomp
from pygra import hamiltonians
from pygra import klist
from pygra import geometry
from pygra import sculpt
from pygra import multilayers
from pygra import dos
from pygra import ldos
from pygra import films
from pygra import kpm
from pygra import current
from pygra import spectrum
from pygra import topology
#from pygra import heterostructures
from pygra import inout
from pygra import operators
from pygra import bandstructure
from pygra import islands
from pygra import ribbon
from pygra import hybrid
from pygra import kdos
from pygra import potentials
from pygra import supercell
from pygra import scftypes
from pygra import indexing
from pygra import meanfield
from pygra import specialgeometry
from pygra import specialhopping
from pygra import timeevolution

import platform


dirname = os.path.dirname(os.path.realpath(__file__))
sys.path.append(dirname+"/../interpreter") # add this path
from interpreter import pycommand


def get_python():
  return pycommand.get_python()

get_anaconda_command = get_python







def get_qhroot():
  """Gets the root path of quantum honeycomp"""
  return os.environ["QHROOT"]+"/"



def create_folder():
  """Creates a temporal folder and goes to that one"""
  os.chdir("/tmp")
  # get the name of the folder
  i = 0
  forig = "qh-tmp-"
  folders = os.listdir(os.getcwd()) # list all the folders 
  while True:
    folder = forig + str(i)
    if not folder in folders:
      break # stop if folder doesn't exist
    i += 1 # increase the number
  os.system("mkdir "+folder)  # create the temporal folder
  os.chdir(folder)  # go to the temporal folder
  return folder  # return the name of the folder






def save_outputs(inipath,tmppath):
  """Save all the results in the original folder"""
  savepath = inipath+"/QH_save" # name of the fodler where ot save
  print("Saving results in",savepath)
  os.system("rm -rf "+savepath) # remove the folder
  os.system("cp -r "+tmppath+"  "+savepath) # remove the folder



def execute_script(name,background=True,mayavi=False):
  """Executes a certain script from the folder utilities"""
  try: qhpath = get_qhroot() # get the main path
  except: qhpath = "" 
  print("Root path",qhpath)
  scriptpath = qhpath+"utilities/"+name # name of the script
  try:
    python = get_anaconda_command("python") # get the anaconda python
  except:
    python = get_python() # get the correct interpreter
  python = pycommand.get_python()
  if background: os.system(python+" "+scriptpath+" &") # execute the script
  else: os.system(python+" "+scriptpath) # execute the script



def computing():
  """Return an object that shows up a window saying computing"""
  qhpath = get_qhroot()
  python = get_python()
  name = qhpath + "interface-pyqt/timer/timer.py"
  import subprocess
#  import psutil
#  print(name)
  subp = subprocess.Popen([python,name]) # execute the command
#  p = psutil.Process(subp.pid) # return process
  return subp


