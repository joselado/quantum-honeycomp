import numpy as np
from scipy.sparse import csc_matrix
import scipy.sparse.linalg as lg

def curvature(h,k=np.array([0.,0.]),dk=0.01,n=2):
  """Calculates the Berry curvature of a kpoint
  of the n lowest bands below E=0"""
  hkgen = h.get_hk_gen()  # function which generates hk
  def get_wf(hk,n2=n):
    """Get lowest bands"""
    e,wf = lg.eigsh(csc_matrix(hk),k=n2,which="LM",sigma=0.0) # get wf
    wfsout = []
    ni = 0 # number of waves found
    wf = wf.transpose() # transpose waves
    for ie,iwf in zip(e,wf):
      if ie<0.0:
        wfsout.append(iwf)
        ni += 1
      if ni==n:
        n = ni + 0 # change global value
        return wfsout # return waves
    return get_wf(hk,n2=n2+2)  # if not eough waves, try with two more...
  # now calculate berry 
  dx = np.array([dk/2.,0.])
  dy = np.array([0.,dk/2.])
  k1 = k - dx - dy  # kpoints
  k2 = k + dx - dy
  k3 = k + dx + dy
  k4 = k - dx + dy
  wk1 = get_wf(k1)  # waves
  wk2 = get_wf(k2)
  wk3 = get_wf(k3)
  wk4 = get_wf(k4)
  def mij(wis,wjs):
    """Calculate matrix"""
    nw = len(wis)
    m = np.matrix(np.zeros((nw,nw))) # create matrix
    for i in range(nw):
      for j in range(nw):
        m[i,j] = wis[i].H * wjs  # matrix element
  phi = np.det(mij(wk1,wk2)*mij(wk2,wk3)*mij(wk3,wk4)*mij(wk4*wk1) ) # deter
  phi = np.arctan2(phi.imag,phi.real)/(dk*dk) # phase
  return phi
   
