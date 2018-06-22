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



def get_anaconda_command(name="python"):
  """Return the path for Anaconda Python, which has pyqt by default"""
  os.system("rm -f /tmp/qh_commands.txt") # remove
  os.system("which -a "+name+"  > /tmp/qh_commands.txt") # run the command
  lines = open("/tmp/qh_commands.txt").read() # read the lines
  lines = lines.split("\n") # split the lines
  del lines[-1] # remove the last one
  print("Found ",len(lines),"python paths\n")
  for l in lines: print(l)
  for l in lines: # loop over pythons
    l = l.split(" ")[-1] # get last line 
    if "anaconda" in l:
      print("\nFound Anaconda ",name,"in",l)
      return l
  print("Anaconda",name,"not found")
  raise




import platform
if platform.system()=="Linux":
  compiler = "/usr/bin/f2py3" # name of the compiler
else: 
  compiler = get_anaconda_command("f2py") # name of the compiler


for name in names:
  folder,f90,so = name[0],name[1],name[2] # different names
  os.chdir("pysrc/fortran/"+folder) # go to the folder
  os.system("rm -f *.so") # remove old libraries
  os.system(compiler+" -c -m "+so+"  "+f90) # compile
  os.system("cp *.so ../../"+so+".so") # copy library
  os.chdir("../../..") # return to parent

