import os

os.chdir("fortran/first_neigh") # move to folder
os.system("sh compile.sh") # compile
os.chdir("../..") # return to parent

os.chdir("fortran/kpm") # move to folder
os.system("sh compile.sh") # compile
os.chdir("../..") # return to parent
