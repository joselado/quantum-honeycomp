### list with the different modes of the main interface

modes = []
modes.append("1d")
modes.append("0d")
modes.append("2d")
modes.append("skyrmion2d")
modes.append("vacancy2d")
modes.append("multilayers2d")
modes.append("huge_0d")
modes.append("interface1d")
modes.append("moire")
modes.append("tbg")

import os
def call_mode(path,name):
  """ Returns a function which calls a particular quantum honeycomp mode"""
  def launch(self):
    """Function to launch a particular calculation, passes qhroot as input"""
    os.system("/usr/bin/python3 "+path+"systems/"+name+"/"+name+".py   "+path +"  &")  
#    os.system("python "+path+"systems/"+name+"/"+name+".py   "+path +"  &")  
  return launch



