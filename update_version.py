#!/usr/bin/python
from __future__ import print_function
import os
import numpy as np
import sys


if len(sys.argv)>1: # if input provided
  if "update" in sys.argv[1]:
    l = np.genfromtxt("VERSION.TXT") # read the version
    n = l[2] # get the version number
    a,b,c = str(int(l[0])),str(int(l[1])),str(int(n)) # create the three numbers
    c= str(int(n + 1)) # increase number
    open("VERSION.TXT","w").write(a+" "+b+" "+c)
else:
  l = np.genfromtxt("VERSION.TXT") # read the version
  n = l[2] # get the version number
  a,b,c = str(int(l[0])),str(int(l[1])),str(int(n)) # create the three numbers
  print("Version number not updated, current is",a,b,c)



l = np.genfromtxt("VERSION.TXT") # read the version
n = l[2] # get the version number
a,b,c = str(int(l[0])),str(int(l[1])),str(int(n)) # create the three numbers


name = "quantum-honeycomp-"+a+"."+b+"."+c # create the new name

os.system("cp VERSION.TXT development") # copy to the development folder
os.system("cp -r development "+name) # copy to a new folder
os.chdir(name) # go to that folder
os.system("./clean") # clean files
os.chdir("..") # go back to the main folder

os.system("mv *.tar old_versions/") # move the old tar files
os.system("tar -cvf "+name+".tar "+name) # create new tar file
os.system("rm -r "+name) # remove the folder
os.system("cp "+name+".tar quantum-honeycomp-latest.tar") # create new tar file

#exit()



# update the version number
#open("VERSION.TXT","w").write(a+"  "+b+"  "+c)

os.system("git add .")
os.system("git commit -m 'Update'")
os.system("git push")


