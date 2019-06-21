from __future__ import print_function
import numpy as np

def bulk1d(g,fac=0.3):
  """Return the indexes of the bulk sites"""
  ymin = np.min(g.y)
  ymax = np.max(g.y)
  out = np.zeros(len(g.r)) # output array
  for i in range(len(g.r)):
    if fac<(g.y[i]-ymin)/(ymax-ymin)<(1.0-fac): out[i] = 1.0
  return out


def bulk0d(g0,fac=0.3):
  """Return the indexes of the bulk sites"""
  g = g0.copy()
  g.center()
  r2 = np.sqrt(g.x**2 + g.y**2 + g.z**2)
  rmax = np.max(r2)
  out = np.zeros(r2.shape) # initialize
  out[r2<(fac*rmax)] = 1.0
  return out


def bulk(g,**kwargs):
    if g.dimensionality==0:  return bulk0d(g,**kwargs)
    elif g.dimensionality==1:  return bulk1d(g,**kwargs)



