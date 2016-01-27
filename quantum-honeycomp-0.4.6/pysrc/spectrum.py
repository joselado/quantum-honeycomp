# library to deal with the spectral properties of the hamiltonian
import numpy as np
import scipy.linalg as lg

def fermi_surface(h,write=True,output_file="FERMI_SURFACE.OUT",
                    e=0.0,nk=50,nsuper=1,reciprocal=False,
                    delta=0.1):
  """Calculates the Fermi surface of a 2d system"""
  if h.dimensionality!=2: raise  # continue if two dimensional
  hk_gen = h.get_hk_gen() # gets the function to generate h(k)
  kxs = np.linspace(0.,nsuper,nk)  # generate kx
  kys = np.linspace(0.,nsuper,nk)  # generate ky
  iden = np.identity(h.intra.shape[0],dtype=np.complex)
  kdos = [] # empty list
  kxout = []
  kyout = []
  if reciprocal: 
    from geometry import get_reciprocal2d
    (ux,uy) = get_reciprocal2d(h.geometry.a1,h.geometry.a2)
  else:
    ux = np.array([1.,0.])
    uy = np.array([0.,1.])
  for x in kxs:
    for y in kxs:
      k = x*ux + y*uy # reciprocal vectors
      hk = hk_gen(k) # get hamiltonian
      gf = ((e+1j*delta)*iden - hk).I # get green function
      kdos.append(-gf.trace()[0,0].imag) # add to the list
      kxout.append(x)
      kyout.append(y)
  if write:  # optionally, write in file
    f = open(output_file,"w") 
    for (x,y,d) in zip(kxout,kyout,kdos):
      f.write(str(x)+ "   "+str(y)+"   "+str(d)+"\n")
    f.close() # close the file
  return (kxout,kyout,d) # return result

