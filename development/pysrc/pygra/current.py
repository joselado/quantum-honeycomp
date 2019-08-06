import numpy as np
import scipy.sparse.linalg as lg
import scipy.linalg as lg

def braket_wAw(w,A):
  w = np.matrix(w) # convert to matrix
  return ((w.T).H*A*w.T)[0,0] # expectation value

def current_operator(h0):
  """Get the current operator"""
  h = h0.copy()
  h = h0.get_multicell() # multicell Hamiltonian
  if h.dimensionality == 0: return None
  elif h.dimensionality == 1:
    if not h.is_multicell: # no multicell
      def fj(k0):
        k = k0[0]
        phik = np.exp(1j*2.*np.pi*k) # complex phase
        jk = 1j*(h.inter*phik - h.inter.H*np.conjugate(phik)) 
        return jk
    else: # multicell Hamiltonian
      def fj(k0):
        k = k0[0]
        jk = h.intra*0. # initialize
        for t in h.hopping:
          phik = np.exp(1j*2.*np.pi*k*t.dir[0]) # complex phase
          jk = jk + 1j*t.m*phik*t.dir[0]
        return jk
    return fj
  else: raise



def gs_current(h,nk=400):
  weighted_current(h,nk=nk)




def fermi_current(h,nk=400,delta=0.5):
  def fun(e):
    return delta/(delta**2+e**2)*2/np.pi
  weighted_current(h,nk=nk,fun=fun)



def weighted_current(h,nk=400,fun=None):
  """Calculate the Ground state current"""
  if fun is None:
    delta = 0.01
    def fun(e): return (-np.tanh(e/delta) + 1.0)/2.0
  jgs = np.zeros(h.intra.shape[0]) # current array
  hkgen = h.get_hk_gen() # generator
  fj = current_operator(h) # current operator
  ks = np.linspace(0.0,1.0,nk,endpoint=False) # k-points
  for k in ks: # loop
    hk = hkgen(k) # Hamiltonian
    (es,ws) = lg.eigh(hk) # diagonalize
    ws = ws.transpose() # transpose
    jk = fj(k) # get the generator
    for (e,w) in zip(es,ws): # loop
      weight = fun(e) # weight
      print(weight)
      d = np.conjugate(w)*ket_Aw(jk,w) # current density
      jgs += d.real*weight # add contribution
#      jgs += (np.abs(w)**2*weight).real # add contribution
  jgs /= nk # normalize
  print("Total current",np.sum(jgs))
  np.savetxt("CURRENT1D.OUT",np.matrix([range(len(jgs)),jgs]).T)

