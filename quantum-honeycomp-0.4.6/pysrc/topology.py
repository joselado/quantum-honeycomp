# library to calculate topological properties

import numpy as np
from scipy.sparse import bmat, csc_matrix
import scipy.linalg as lg

def majorana_invariant(intra,inter,kp=1000):
  """ Calculates the topological invariant for a 1d topological
  superconductor, returns a list of determinants of the upper diagonal"""
  ks = np.arange(0.5,1.5,1.0/kp) # create kpoints
  dets = [0. for k in ks]  # create list with determinants
  # rotate the matrices into a non diagonal block for
  # assume that h = sz*h_0 + i sy * h_delta
  # and perfor a rotation e^-i*pi/4 sy
  rot = intra * 0.0
  n = len(intra)/2 # number of orbitals including spin
  # create rotation matrix
  intra_a = np.matrix([[0.0j for i in range(n)] for j in range(n)])
  inter_a = np.matrix([[0.0j for i in range(n)] for j in range(n)])
  print csc_matrix(intra-intra.H)
  for i in range(n):
    for j in range(n):
      # couples electron in i with hole in j
      s = 1.0
      if i<j:
        s = -1.0
      intra_a[i,j] = 1j*s*intra[2*i,2*j+1] + intra[2*i,2*j]
#      intra_a[i,j] = intra[2*i,2*j+1]
      inter_a[i,j] = 1j*s*inter[2*i,2*j+1] + inter[2*i,2*j]
#  print csc_matrix(intra_a)
  fm = open("WINDING_MAJORANA.OUT","w")
  fm.write("# kpoint,     phase")
  for k in ks:
    tk = inter_a * np.exp(1j*2.*np.pi*k)
    hk = intra_a + tk + tk.H # kdependent k
    det = lg.det(hk)
    phi = np.arctan2(det.imag,det.real)
    fm.write(str(k)+"    "+str(phi)+"\n")
  fm.close()  # close the file



def write_berry(h,kpath,dk=0.01):
  """Calculate and write in file the Berry curvature"""
  be = berry_curvature(h,kpath,dk=dk)
  fo = open("BERRY_CURVATURE.OUT","w") # open file
  for (k,b) in zip(kpath,be):
    fo.write(str(k[0])+"   ")
    fo.write(str(k[1])+"   ")
    fo.write(str(b)+"\n")
  fo.close() # close file





def berry_curvature(h,k,dk=0.01):
  """ Calculates the Berry curvature of a 2d hamiltonian"""
  try: 
    a = k[0][0] # check if array of k's on input
    return [berry_curvature(h,ik,dk=dk) for ik in k] # recursive call
  except: pass # if no list, continue
  if h.dimensionality != 2: raise # only for 2d
  k = np.array([k[0],k[1]]) 
  dx = np.array([dk,0.])
  dy = np.array([0.,dk])
  occ = occ_states2d  # functio to get the occupied states
  # get the waves
  wf1 = occ(h,k-dx-dy) 
  wf2 = occ(h,k+dx-dy) 
  wf3 = occ(h,k+dx+dy) 
  wf4 = occ(h,k-dx+dy) 
  # get the uij  
  m = uij(wf1,wf2)*uij(wf2,wf3)*uij(wf3,wf4)*uij(wf4,wf1)
  d = lg.det(m) # calculate determinant
  phi = np.arctan2(d.imag,d.real)/(4.*dk*dk)
  return phi

def occ_states2d(h,k):
  """ Returns the WF of the occupied states in a 2d hamiltonian"""
  hk_gen = h.get_hk_gen() # get hamiltonian generator
  hk = hk_gen(k) # get hamiltonian
  es,wfs = lg.eigh(hk) # diagonalize
  wfs = wfs.transpose()
  occwf = []
  for (ie,iw) in zip(es,wfs):  # loop over states
    if ie < 0:  # if below fermi
      if np.abs(iw[0])>0.0001:
        iw = iw*np.conjugate(iw[0])/np.abs(iw[0])
      occwf.append(iw)  # add to the list
  return occwf 



def uij(wf1,wf2):
  """ Calcultes the matrix produxt of two sets of input wavefunctions"""
  m = np.matrix(np.zeros((len(wf1),len(wf2)),dtype=np.complex))
  for i in range(len(wf1)):
    for j in range(len(wf2)):
      m[i,j] = sum(np.conjugate(wf1[i])*wf2[j])
  return m


def precise_chern(h,dk=0.01):
  """ Calculates the chern number of a 2d system """
  from scipy import integrate
  err = {"epsabs" : 0.01, "epsrel": 0.01,"limit" : 20}
#  err = [err,err]
  def f(x,y): # function to integrate
    return berry_curvature(h,np.array([x,y]),dk=dk)
  c = integrate.dblquad(f,0.,1.,lambda x : 0., lambda x: 1.,epsabs=1.,
                          epsrel=1.)
  return c[0]/(2.*np.pi)


def chern(h,dk=0.01,nk=10):
  """ Calculates the chern number of a 2d system """
  c = 0.0
  for x in np.linspace(0.,1.,nk):
    for y in np.linspace(0.,1.,nk):
      c += berry_curvature(h,np.array([x,y]),dk=dk)
  c = c/(2.*np.pi*nk*nk)
  return c


def z2_invariant(h,nk=100):
  """Calculates the Z2 invariant by the path evolution algorithm"""
  import klist
  kst = klist.tr_klist(nk=nk)  # class with the klist
  wfs1 = [occ_states2d(h,k) for k in kst.path1]  # lower path
  wfs2 = [occ_states2d(h,k) for k in kst.path2]  # lower path
  wfsc1 = [occ_states2d(h,k) for k in kst.common]  # lower path
  # related with the other by symmetry
  wfsc2 = [wfsc1[-i] for i in range(len(wfsc1))] # invert the order
  wfs = wfs1 + wfsc1 + wfs2 + wfsc2  # all the paths
  phi = 0.0 # initialize phase
  # connection contribution
  for i in range(len(wfs)):
    m = uij(wfs[i-1],wfs[i]) # matrix of wavefunctions
    d = lg.det(m) # calculate determinant
    phi += np.arctan2(d.imag,d.real) # add contribution
  phi = phi/(2.*np.pi)
  # curvature contribution
  halfc = 0.0
  for ik in np.linspace(-.5,.5,nk):
    for jk in np.linspace(-.0,.5,nk/2):
      halfc += berry_curvature(h,np.array([ik,jk]))
  halfc = halfc/(2.*np.pi*nk*nk)
  print halfc,phi
  return phi-halfc
