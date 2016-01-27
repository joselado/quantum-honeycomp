import subprocess
import os
import sys
# import the different libraries for the quantum honeycomp
import hamiltonians
import klist
import geometry
import input_tb90 
import sculpt
import numpy as np
import interactions
import multilayers
import dos
import ldos
import kpm

def get_qhroot():
  """Gets the root path of quantum honeycomp"""
  # get the path
  m = subprocess.Popen(["which","quantum-honeycomp"],stdout=subprocess.PIPE)
  m = m.communicate()
  p = m[0] # get the element of the tuple
  fs = p.split()[0].split("/") # list with folder steps
  p = ""
  for i in range(len(fs)-2):  # create path to the main folder
    f = fs[i] # folder name
    p += f+"/" # add to path
  return p # return the path


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




# setup the different paths to needed libraries
qh_root = get_qhroot() # get the path to the main folder
pyroot = qh_root + "pyroot" # path to python scripts
f90root = qh_root + "f90root" # path to tb90 program


