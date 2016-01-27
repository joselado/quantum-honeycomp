# kernel polynomial method libraries

import scipy.sparse.linalg as lg
from scipy.sparse import csc_matrix as csc
import numpy.random as rand
import numpy as np

def get_moments(v,m,n=100,use_fortran=True):
  """ Get the first n moments of a certain vector
  using the Chebychev recursion relations"""
  if use_fortran:
    from kpmf90 import get_momentsf90 # fortran routine
    from scipy.sparse import coo_matrix
    mo = coo_matrix(m) # convert to coo matrix
    vo = v.todense() # convert to conventional vector
    vo = np.array([vo[i,0] for i in range(len(vo))])
# call the fortran routine
    mus = get_momentsf90(mo.row+1,mo.col+1,mo.data,vo,n) 
    return mus # return fortran result
  else:
    mus = np.array([0.0j for i in range(2*n)]) # empty arrray for the moments
    a = v.copy() # first vector
    am = v.copy() # zero vector
    a = m*v  # vector number 1
    bk = (np.transpose(np.conjugate(v))*v)[0,0] # scalar product
    bk1 = (np.transpose(np.conjugate(v))*a)[0,0] # scalar product
    mus[0] = bk  # mu0
    mus[1] = bk1 # mu1
    for i in range(1,n): 
      ap = 2*m*a - am # recursion relation
      bk = (np.transpose(np.conjugate(a))*a)[0,0] # scalar product
      bk1 = (np.transpose(np.conjugate(ap))*a)[0,0] # scalar product
      mus[2*i] = 2.*bk
      mus[2*i+1] = 2.*bk1
      am = a +0. # new variables
      a = ap+0. # new variables
    mu0 = mus[0] # first
    mu1 = mus[1] # second
    for i in range(1,n): 
      mus[2*i] +=  - mu0
      mus[2*i+1] += -mu1 
    return mus




def get_momentsA(v,m,n=100,A=None):
  """ Get the first n moments of a certain vector
  using the Chebychev recursion relations"""
  mus = np.array([0.0j for i in range(2*n)]) # empty arrray for the moments
  am = v.copy() # zero vector
  a = m*v  # vector number 1
  bk = (np.transpose(np.conjugate(v))*A*v)[0,0] # scalar product
  bk1 = (np.transpose(np.conjugate(v))*A*a)[0,0] # scalar product
  mus[0] = bk  # mu0
  mus[1] = bk1 # mu1
  for i in range(1,n): 
    ap = 2*m*a - am # recursion relation
    bk = (np.transpose(np.conjugate(a))*A*a)[0,0] # scalar product
    bk1 = (np.transpose(np.conjugate(ap))*A*a)[0,0] # scalar product
    mus[2*i] = 2.*bk
    mus[2*i+1] = 2.*bk1
    am = a +0. # new variables
    a = ap+0. # new variables
  mu0 = mus[0] # first
  mu1 = mus[1] # second
  for i in range(1,n): 
    mus[2*i] +=  - mu0
    mus[2*i+1] += -mu1 
  return mus












def full_trace(m_in,n=200,use_fortran=True):
  """ Get full trace of the matrix"""
  m = csc(m_in) # saprse matrix
  nd = m.shape[0] # length of the matrix
  mus = np.array([0.0j for i in range(2*n)])
#  for i in range(ntries):
  for i in range(nd):
    mus += local_dos(m_in,i=i,n=n,use_fortran=use_fortran)
  return mus/nd









def local_dos(m_in,i=0,n=200,use_fortran=True):
  """ Calculates local DOS using the KPM"""
  m = csc(m_in) # saprse matrix
  nd = m.shape[0] # length of the matrix
  mus = np.array([0.0j for j in range(2*n)])
  v = rand.random(nd)*0.
  v[i] = 1.0 # vector only in site i 
  v = csc(v).transpose()
# get the chebychev moments
  mus += get_moments(v,m,n=n,use_fortran=use_fortran) 
  return mus



def random_trace(m_in,ntries=20,n=200):
  """ Calculates local DOS using the KPM"""
  m = csc(m_in) # saprse matrix
  nd = m.shape[0] # length of the matrix
  mus = np.array([0.0j for j in range(2*n)])
  for i in range(ntries): # loop over tries
    #v = rand.random(nd) - .5
    v = rand.random(nd) -.5 + 1j*rand.random(nd) -.5j
    v = v/np.sqrt(v.dot(np.conjugate(v))) # normalize the vector
    v = csc(v).transpose()
    mus += get_moments(v,m,n=n) # get the chebychev moments
  return mus/ntries



def random_trace_AB(m_in,ntries=20,n=200,A=None,B=None):
  """ Calculates local DOS using the KPM"""
  m = csc(m_in) # saprse matrix
  nd = m.shape[0] # length of the matrix
  mus = np.array([0.0j for j in range(2*n)])
  for i in range(ntries): # loop over tries
    #v = rand.random(nd) - .5
    v = rand.random(nd) -.5 + 1j*rand.random(nd) -.5j
    v = v/np.sqrt(v.dot(v)) # normalize the vector
    v = csc(v).transpose()
    mus += get_momentsA(v,m,n=n,A=A*B) # get the chebychev moments
    mus += get_momentsA(v,-m,n=n,A=B*A) # get the chebychev moments
  return mus/ntries



def full_trace_AB(m_in,ntries=20,n=200,A=None,B=None):
  """ Calculates local DOS using the KPM"""
  m = csc(m_in) # saprse matrix
  nd = m.shape[0] # length of the matrix
  mus = np.array([0.0j for j in range(2*n)])
  for i in range(nd): # loop over tries
    #v = rand.random(nd) - .5
    v = rand.random(nd)*0.
    v[i] = 1.0 # vector only in site i 
    v = csc(v).transpose()
    mus += get_momentsA(v,m,n=n,A=A*B) # get the chebychev moments
    mus += get_momentsA(v,-m,n=n,A=B*A) # get the chebychev moments
  return mus/nd











def generate_profile(mus,xs):
  """ Uses the Chebychev expansion to create a certain profile"""












def generate_profile(mus,xs):
  """ Uses the Chebychev expansion to create a certain profile"""
  # initialize polynomials
#  xs = np.array([0.])
  tm = xs.copy()*0. +1.
  t = xs.copy()
  ys = mus[0] # first term
  mus = jackson_kernel(mus)
#  mus = fejer_kernel(mus)
  # loop over all contributions
  for i in range(1,len(mus)):
    mu = mus[i]
    ys += 2.*mu*t # add contribution
#    print tm
    tp = 2.*xs*t - tm # chebychev recursion relation
    tm = t + 0.
    t = 0. + tp # next iteration
  ys = ys/np.sqrt(1.-xs*xs) # prefactor
  ys = ys.real
  return ys



def jackson_kernel(mus):
  """ Modify coeficient using the Jackson Kernel"""
  mo = mus.copy() # copy array
  n = len(mo)
  pn = np.pi/(n+1.) # factor
  for i in range(n):
    fac = ((n-i+1)*np.cos(pn*i)+np.sin(pn*i)/np.tan(pn))/(n+1)
    mo[i] *= fac
#    print fac,pn
  return mo


def fejer_kernel(mus):
  """Default kernel"""
  n = len(mus)
  mo = mus.copy()
  for i in range(len(mus)):
    mo[i] *= (1.-float(i)/n) 
  return mo














