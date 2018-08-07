from __future__ import print_function
import numpy as np

import extract

def swave(h0,name="SWAVE.OUT",rep=1):
  """Write the swave pairing of a Hamiltonian"""
  h = h0.supercell(rep)
  f1 = open("AMPLITUDE_"+name,"w")
  f2 = open("PHASE_"+name,"w")
  f3 = open(name,"w")
  if not h.has_spin: raise
  if not h.has_eh: raise
  ds = extract.swave(h.intra) # get the pairing
  r = h.geometry.r
  f1.write("# x y |Delta|\n")
  f2.write("# x y |phi|\n")
  f3.write("# x y ReD  ImD\n")
  for i in range(len(r)):
    ri = r[i] # position
    di = ds[i]
    f1.write(str(ri[0])+"    ")
    f2.write(str(ri[0])+"    ")
    f3.write(str(ri[0])+"    ")
    f1.write(str(ri[1])+"    ")
    f2.write(str(ri[1])+"    ")
    f3.write(str(ri[1])+"    ")
    f1.write(str(np.abs(di))+"    ")
    f2.write(str(np.angle(di))+"    ")
    f3.write(str(di.real)+"    ")
    f3.write(str(di.imag)+"    ")
    f1.write("\n")
    f2.write("\n")
    f3.write("\n")
  f1.close()
  f2.close()
  f3.close()


def hopping(h,name="HOPPING.OUT",reps=0):
  """Write the magnitude of the hopping in a file"""
  if h.has_eh: raise
  if h.has_spin: (ii,jj,ts) = extract.hopping_spinful(h.intra)
  else: (ii,jj,ts) = extract.hopping_spinless(h.intra)
  f = open(name,"w") # write file
  for (i,j,t) in zip(ii,jj,ts):
    f.write(str(h.geometry.r[i][0])+"  ")
    f.write(str(h.geometry.r[i][1])+"  ")
    f.write(str(h.geometry.r[j][0])+"  ")
    f.write(str(h.geometry.r[j][1])+"  ")
    f.write(str(t)+"\n")
  f.close()



def mz(h,name="MZ.OUT"):
  if h.has_eh: raise
  if h.has_spin: ms = extract.mz(h.intra)
  else: raise
  np.savetxt(name,np.matrix([range(len(ms)),ms]).T)



def magnetization(h):
  """Write all the magnetizations"""
  if h.has_eh: raise
  if h.has_spin: 
    mx = extract.mx(h.intra)
    my = extract.my(h.intra)
    mz = extract.mz(h.intra)
  else: raise
  np.savetxt("MAGNETIZATION_X.OUT",np.matrix([h.geometry.x,h.geometry.y,mx]).T)
  np.savetxt("MAGNETIZATION_Y.OUT",np.matrix([h.geometry.x,h.geometry.y,my]).T)
  np.savetxt("MAGNETIZATION_Z.OUT",np.matrix([h.geometry.x,h.geometry.y,mz]).T)




