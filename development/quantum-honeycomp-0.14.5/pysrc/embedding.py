# library to perform embedding calculations
from __future__ import print_function

import green
import numpy as np


def dos_impurity(h,vc=None,energies=np.linspace(-.5,.5,20),
                   mode="adaptive",delta=0.01,nk=50,silent=True,
                   use_generator=False):
  """ Calculates the green function using the embedding technique"""
  if vc is None: vc = h.intra  # assume perfect
  ds,dsv = [],[] # empty lists
  iden = np.matrix(np.identity(h.intra.shape[0],dtype=np.complex))
  if use_generator:
    getgreen = green.green_generator(h,nk=nk) # get the function
    def get_green(energy): return getgreen(energy,delta=delta)
  else: # use the conventional method
    def get_green(energy):
      return green.bloch_selfenergy(h,energy=energy,delta=delta,nk=nk,
                                       mode=mode)
  for energy in energies:
    g,selfe = get_green(energy) # get Green and selfenergy
    emat = iden*(energy + delta*1j)  # E +i\delta 
    gv = (emat - vc -selfe).I   # Green function of a vacancy, with selfener
    ds.append(-g.trace()[0,0].imag)  # save DOS of the pristine
    dsv.append(-gv.trace()[0,0].imag)  # save DOS of the defected
    if not silent: print("Done",energy)
  ds = np.array(ds)
  dsv = np.array(dsv)
  np.savetxt("DOS_PRISTINE.OUT",np.matrix([energies,ds]).T)
  np.savetxt("DOS_DEFECTIVE.OUT",np.matrix([energies,dsv]).T)
  return ds,dsv # return object





