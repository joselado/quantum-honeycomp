
from __future__ import print_function
from scipy.interpolate import interp1d
import numpy as np

# perform interpolation a matrix 
def intermatrix(fin,xs=np.linspace(-5.0,5.0,20)):
  """Return a function capable to interpolate the
  input function between values in the interval. Values
  outside the intervale will be returned as zero"""
  m = fin(0) # call once
  ar = np.zeros((len(xs),m.shape[0],m.shape[1])) # empty array, real part
  ai = np.zeros((len(xs),m.shape[0],m.shape[1])) # empty array, imag part
  for i in range(len(xs)): # loop
    m = fin(xs[i]) # call the function
    ar[i,:,:] = m.real # real part
    ai[i,:,:] = m.imag # imaginary part
  zero = np.matrix(np.zeros(m.shape,dtype=np.complex)) # zero matrix
  fr = interp1d(xs, ar, axis=0,kind=3,fill_value=zero)
  fi = interp1d(xs, ai, axis=0,kind=3,fill_value=zero)
  def fout(e): # output function
    return np.matrix(fr(e) + 1j*fi(e)) # return 
  return fout # return function


