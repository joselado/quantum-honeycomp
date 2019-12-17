# library to perform embedding calculations
from __future__ import print_function

from . import green
from . import parallel
from . import algebra
from . import timing
import numpy as np
import scipy.linalg as lg


def dos_impurity(h,vc=None,energies=np.linspace(-.5,.5,20),
                   mode="adaptive",delta=0.01,nk=50,silent=True,
                   use_generator=False):
  """ Calculates the green function using the embedding technique"""
  if vc is None: vc = h.intra  # assume perfect
  iden = np.identity(h.intra.shape[0],dtype=np.complex)
  if use_generator:
    getgreen = green.green_generator(h,nk=nk) # get the function
    def get_green(energy): return getgreen(energy,delta=delta)
  else: # use the conventional method
    def get_green(energy):
      return green.bloch_selfenergy(h,energy=energy,delta=delta,nk=nk,
                                       mode=mode)
  def pfun(energy): # function to parallelize
    g,selfe = get_green(energy) # get Green and selfenergy
    emat = iden*(energy + delta*1j)  # E +i\delta 
    gv = lg.inv(emat - vc -selfe)   # Green function of a vacancy, with selfener
    d = -np.trace(g).imag  # save DOS of the pristine
    dv = -np.trace(gv).imag  # save DOS of the defected
    if not silent: print("Done",energy)
    return [d,dv]
  out = np.array(parallel.pcall(pfun,energies)) # compute
  ds,dsv = out[:,0],out[:,1] # get the different data
  np.savetxt("DOS_PRISTINE.OUT",np.array([energies,ds]).T)
  np.savetxt("DOS_DEFECTIVE.OUT",np.array([energies,dsv]).T)
  return ds,dsv # return object





def bulk_and_surface(h1,nk=100,energies=np.linspace(-1.,1.,100),**kwargs):
  """Get the surface DOS of an interface"""
  from scipy.sparse import csc_matrix,bmat
  if h1.dimensionality==2:
      kpath = [[k,0.,0.] for k in np.linspace(0.,1.,nk)]
  else: raise
  ik = 0
  h1 = h1.get_multicell() # multicell Hamiltonian
  tr = timing.Testimator("DOS") # generate object
  dos_bulk = energies*0.0
  dos_sur = energies*0.0
  for k in kpath:
    tr.remaining(ik,len(kpath)) # generate object
    ik += 1
    outs = green.surface_multienergy(h1,k=k,energies=energies,**kwargs)
    dos_bulk += np.array([-algebra.trace(g[1]).imag for g in outs])
    dos_sur += np.array([-algebra.trace(g[0]).imag for g in outs])
  dos_bulk /= len(kpath)
  dos_sur /= len(kpath)
  np.savetxt("DOS.OUT",np.array([energies,dos_bulk,dos_sur]).T)
  return energies,dos_bulk,dos_sur




