from __future__ import print_function
import subprocess
import os
import sys
# import the different libraries for the quantum honeycomp
import hamiltonians
import klist
import geometry
import input_tb90 
import sculpt
import interactions
import multilayers
import dos
import ldos
import kpm
import current
import spectrum
import topology
import heterostructures
import numpy as np
import inout
import operators
import bandstructure
import islands
import hybrid
import kdos
import potentials
import supercell
import scftypes
import indexing
import specialgeometry
import specialhopping


import platform

def get_python():
  try:
    return get_anaconda_command("python") # return anaconda
  except:
    if platform.system()=="Linux":
      python = "/usr/bin/python3" # Python 3
    else:
      python ="python" # Python for mac
    return python





def get_anaconda_command(name="python"):
  """Return the path for Anaconda Python, which has pyqt by default"""
  os.system("rm -f /tmp/qh_commands.txt") # remove
  os.system("which -a "+name+"  > /tmp/qh_commands.txt") # run the command
  lines = open("/tmp/qh_commands.txt").read() # read the lines
  lines = lines.split("\n") # split the lines
  del lines[-1] # remove the last one
  print("Found ",len(lines),"python paths\n")
  for l in lines: print(l)
  for l in lines: # loop over pythons
    l = l.split(" ")[-1] # get last line 
    if "anaconda" in l:
      print("\nFound Anaconda ",name,"in",l)
      return l
  print("Anaconda",name,"not found")
  raise

















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
  fends = ["*.OUT","*.png","*.in"]
  name = "QH_save" # name of the folder
  savepath = inipath + "/"+name
  print("Saving results in",savepath)
  os.system("mkdir "+savepath) # create folder
  for e in fends:
    os.system("rm -f "+savepath+"/"+e) # clean
    os.system("cp "+tmppath+"/"+e+"  "+savepath) # copy 
  return


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
#    if mayavi: python = "/usr/bin/python3" # use the native
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


