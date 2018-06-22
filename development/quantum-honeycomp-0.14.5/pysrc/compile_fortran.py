import os

# names is a lists with pairs of name of folder, f90 file and .so file

names = [("first_neigh","first_neighborsf90.f90","first_neighborsf90")] 
names += [("kpm","kpm.f90","kpmf90")] 
names += [("dos","dos.f90","dosf90")] 
names += [("berry","berry_curvature.f90","berry_curvaturef90")] 
names += [("gauss_inv","gauss_inv.f90","gauss_invf90")] 
names += [("clean_geometry","clean_geometryf90.f90","clean_geometryf90")] 
names += [("correlators","correlatorsf90.f90","correlatorsf90")]
names += [("supercell","supercellf90.f90","supercellf90")] 
names += [("classicalspin","classicalspinf90.f90","classicalspinf90")] 
names += [("density_matrix","density_matrixf90.f90","density_matrixf90")] 
names += [("kanemele","kanemelef90.f90","kanemelef90")] 
names += [("green","greenf90.f90","greenf90")] 
names += [("specialhopping","specialhoppingf90.f90","specialhoppingf90")] 
names += [("tails","tailsf90.f90","tailsf90")] 

compiler = "/usr/bin/f2py3" # name of the compiler
compiler = "f2py" # name of the compiler
os.system("rm *.so") # clean all the .so files

for name in names:
  folder,f90,so = name[0],name[1],name[2] # different names
  os.chdir("fortran/"+folder) # go to the folder
#  open("compile_custom.sh","w").write(open("compile.sh").read().replace("XXX",compiler))
#  script = open("compile.sh").read().replace("f2py2.7","XXX")
#  open("compile.sh","w").write(script)
  os.system("rm -f *.so") # remove old libraries
#  os.system("bash compile_custom.sh") # compile
#  os.system("cp compile_custom.sh compile.sh")
  os.system(compiler+" -llapack -c -m "+so+"  "+f90) # compile
  os.system("cp *.so ../../"+so+".so") # copy library
  os.chdir("../..") # return to parent

