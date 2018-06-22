#!/usr/bin/python

### This script will download and install the latest quantum honeycomp ###

import os

name = "quantum-honeycomp-latest.tar" # name of the file
path = "https://github.com/joselado/quantum-honeycomp/raw/master/quantum-honeycomp-latest.tar"
fname = "qh_temporal" # temporal folder for quantum honeycomp 
os.system("rm -rf "+fname) # remove temporal folder
os.system("mkdir "+fname) # create temporal folder
os.chdir(fname) # go to the temporal folder

import platform
if platform.system()=="Linux":
  os.system("wget "+path) # download latest version from git
else: # for Mac
  os.system("curl -LO "+path) # download latest version from git



os.system("tar -xvf "+name) # untar file
fols = next(os.walk('.'))[1] # get the directory
if len(fols)>1:
  print("Something wrong")
  exit()

folder = fols[0] # quantum honeycomp folder
os.chdir("..") # go back
os.system("rm -rf "+folder) # remove folder in case it exists
os.system("mv "+fname+"/"+folder+" ./") # move the folder
os.system("rm -rf "+fname) # remove temporal folder
os.chdir(folder) # go to the quantum honeycomp folder
os.system("python install") # install quantum honeycomp



