
# this script installs the quantum honeycomp package
import os
import sys


path = os.path.dirname(os.path.realpath(__file__)) # current path
sys.path.append(path+"/pysrc/") # add the interpreter
from interpreter import pycommand

#pycommand.install_python() # install the correct python dist
#pycommand.install_dependencies() # install the dependencies
pycommand.add_to_path() # install the correct python dist

