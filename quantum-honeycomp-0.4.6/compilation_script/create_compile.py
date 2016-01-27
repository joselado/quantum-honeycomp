
# PART TO EDIT
#################################################
#################################################
# options for the creation of the compilation script 
compiler = "gfortran"  # compiler used
#options = "-O3 -fopenmp -ffast-math -fbounds-check -fbacktrace"  # full optimization
options = "-O3 -ffast-math -fbounds-check -fbacktrace -fopenmp"  # full optimization
#options = "-O3"  # full optimization
#options = "-O3 "  # full optimization
#options = "-fbounds-check  -fopenmp"  # full optimization
#lapack = "-llapack"  # lapack library flag
lapack = "-L/usr/lib/ -llapack"  # lapack library flag
#lapack = "-L/usr/lib/lapack/ -l lapack  -L/usr/lib/libblas -l blas"  # lapack library flag
lapack = "-L/home/jose/apps/lapack-3.6.0 -l lapack -l blas "  # lapack library flag
#blas =" -L/usr/lib/libblas.a"
#################################################
#################################################





# The following doesn't have to be edited!!!
#################################################
#################################################
#################################################
#################################################
#################################################
#################################################
#################################################
#################################################
#################################################
#################################################
#################################################
#################################################


import os 
from os.path import isfile, join

def get_lapack():
  """Get the location of the lapack library"""
  os.system("locate liblapack.a > .lapack.txt")  


direction = "../f90src"  # direction where the files are


os.chdir(direction)

# get all the files
pwd = os.getcwd()

files = os.listdir(pwd)

# retain the .f90 files
allf90 = []
for f in files:
  end = f.split(".")
  if end[-1]=="f90":
      allf90.append(f)

# f90 files, ordered according to compilation

# remove from the whole list
taskf90 = []
taskf90 += ["density_of_states.f90"]
taskf90 += ["berry.f90"]
taskf90 += ["dielectric_response.f90"]


for f in taskf90:
  allf90.remove(f)




# first files
f90 = ["io_files.f90"]
f90 += ["inputs.f90"]
f90 += ["sparse.f90"]
f90 += ["system_variables.f90"]
f90 += ["expectation_values_variables.f90"]
f90 += ["sintax_hamiltonian.f90"]
f90 += ["bands_variables.f90"]
f90 += ["mean_field_variables.f90"]
f90 += ["berry_variables.f90"]
f90 += ["save_mf_routines.f90"]
f90 += ["offdiagonal.f90"]
#f90 += ["dielectric_operators.f90"]

# f90 variables
varf90 = []
for f in allf90:
  if "_variables" in f:
    varf90.append(f)
# remove from the whole list
for f in varf90:
  allf90.remove(f)

allf90 = varf90 + allf90 + taskf90

allf90.remove("tb90.f90")

f90 += allf90

fc = open("compilation.sh","w")

# compile each file
for f in f90:
  fc.write(compiler+" "+options+"  -c "+f+" "+lapack+"\n")


fc.write(compiler+" "+options+" -o tb90.x *.f90 "+lapack+"\n")


fc.close()
os.system("chmod +x compilation.sh")
os.system("./compilation.sh")
os.system("rm *.mod")



