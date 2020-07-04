# library to perform embedding calculations
from __future__ import print_function

from . import green
from . import parallel
from . import algebra
from . import timing
import numpy as np
import scipy.linalg as lg
import os
from .increase_hilbert import full2profile



class Embedding():
    """Define an embedding object"""
    def __init__(self,h,m=None):
        self.h0 = h.copy() # Pristine Hamiltonian
        if m is not None: self.m = m # provided matrix
        else: self.m = h.intra.copy() # pristine one
    def ldos(self,e=0.0,delta=1e-2,nsuper=3,nk=100,**kwargs):
        """Compute the local density of states"""
        h = self.h0
        g,selfe = green.supercell_selfenergy(h,e=e,delta=delta,nk=nk,
                nsuper=nsuper)
        ms = onsite_supercell(h,nsuper)
        n = self.m.shape[0]
        ms = onsite_defective_central(h,self.m,nsuper)
        ns = ms.shape[0] # dimension of the supercell
        iden = np.identity(ns,dtype=np.complex) # identity
        emat = iden*(e + delta*1j) # energy matrix
        gv = algebra.inv(emat - ms -selfe)   # Defective Green function 
        ds = [-gv[i,i].imag for i in range(ns)] # LDOS
        ds = full2profile(h,ds,check=False) # resum if necessary
        ds = np.array(ds) # convert to array
        gs = self.h0.geometry.supercell(nsuper)
        x,y = gs.x,gs.y
        return x,y,ds
    def multildos(self,es=np.linspace(-2.,2.,20),**kwargs):
        """Compute the ldos at different energies"""
        os.system("rm -rf MULTILDOS")
        os.system("mkdir MULTILDOS")
        ds = [] # total DOS
        fo = open("MULTILDOS/MULTILDOS.TXT","w")
        for e in es:
            (x,y,d) = self.ldos(e=e,**kwargs) # compute LDOS
            ds.append(np.mean(d)) # total DOS
            name0 = "LDOS_"+str(e)+"_.OUT" # name
            name = "MULTILDOS/"+name0
            fo.write(name0+"\n") # name of the file
            np.savetxt(name,np.array([x,y,d]).T) # save data
        np.savetxt("MULTILDOS/DOS.OUT",np.array([es,ds]).T)








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




def onsite_defective_central(h,m,nsuper):
    return onsite_supercell(h,nsuper,mc=m)


def onsite_supercell(h,nsuper,mc=None):
    if h.dimensionality!=2: return NotImplemented
    inds = []
    k = 0
    n = nsuper*nsuper # number of cells
    intrasuper = [[None for j in range(n)] for i in range(n)]
    for i in range(nsuper):
      for j in range(nsuper):
        inds += [(i,j)]
        k += 1
    from scipy.sparse import bmat
    from scipy.sparse import csc_matrix as csc
    tx = csc(h.tx)
    ty = csc(h.ty)
    txy = csc(h.txy)
    txmy = csc(h.txmy)
    intra = csc(h.intra)
    if mc is None: mc = intra
    else: mc = csc(mc)
    for i in range(n):
      intrasuper[i][i] = intra # intracell
      (x1,y1) = inds[i]
      for j in range(n):
        (x2,y2) = inds[j]
        dx = x2-x1
        dy = y2-y1
        if dx==1 and  dy==0: intrasuper[i][j] = tx
        if dx==-1 and dy==0: intrasuper[i][j] = tx.H
        if dx==0 and  dy==1: intrasuper[i][j] = ty
        if dx==0 and  dy==-1: intrasuper[i][j] = ty.H
        if dx==1 and  dy==1: intrasuper[i][j] = txy
        if dx==-1 and dy==-1: intrasuper[i][j] = txy.H
        if dx==1 and  dy==-1: intrasuper[i][j] = txmy
        if dx==-1 and dy==1: intrasuper[i][j] = txmy.H
    if nsuper%2==1: # central cell
        ii=int(n/2)
        intrasuper[ii][ii] = mc # central onsite
    else:
        ii=int(n/2)
        ii = ii - int(nsuper/2)
        intrasuper[ii][ii] = mc # central onsite
    intrasuper = bmat(intrasuper).todense() # supercell
    return intrasuper






