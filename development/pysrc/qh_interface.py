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


import platform

def get_python():
  
  if platform.system()=="Linux":
    python = "/usr/bin/python3" # Python 3
  else:
    python ="python" # Python for mac
  return python





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
  if mayavi: python = "/usr/bin/python2"
  else: python = get_python() # get the correct interpreter
  if background: os.system(python+" "+scriptpath+" &") # execute the script
  else: os.system(python+" "+scriptpath) # execute the script


