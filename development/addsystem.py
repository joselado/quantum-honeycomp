#!/usr/bin/python
from __future__ import print_function
import os


pwd = os.path.dirname(os.path.realpath(__file__))
#pwd = os.getcwd() # get the current location

qhpath = pwd + "/bin/quantum-honeycomp"


# different files for Linux and Mac
import platform
if platform.system()=="Linux":
  if os.path.isdir(os.environ["HOME"]+"/.local/bin"):
    binpath = os.environ["HOME"]+"/.local/bin/quantum-honeycomp" # path to .bashrc
    os.system("rm -f "+binpath)
    os.system("ln -s "+qhpath+"  "+binpath)
    print("Added link to ",binpath)
    exit()



# different files for Linux and Mac
import platform
if platform.system()=="Linux":
  bashrc = os.environ["HOME"]+"/.bashrc" # path to .bashrc
  print("Detected Linux system")
else:
  bashrc = os.environ["HOME"]+"/.bash_profile" # path to .bashrc
  print("Detected Mac system")
print("Adding Quantum Honeycomp to your ",bashrc)

route = "\n###############################\n"
route += "# Added by Quantum Honeycomp\n"
route += "###############################\n"
route += "  export PATH=\""+pwd+"/bin\":$PATH\n"
#route += "  export QHROOT=\""+pwd+"\"\n"
route += "###############################\n"


f = open(bashrc).read() # read bashrc
open(bashrc,"w").write(f+route) # write in the bashrc

print("Added \n"+route+"\n route to ",bashrc)




