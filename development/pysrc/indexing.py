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



