from __future__ import print_function, division
import numpy as np

def cnpot(n=4,k=0.0,v=1.0,angle=0.):
  """Returns a function that generates a potential
  with C_n symmetry"""
  if n==0: return lambda r: v
  if n%2==0: f = np.cos # even 
  if n%2==1: f = np.sin # even 
  def fun(r):
    """Function with the potential"""
    x0,y0 = r[0],r[1]
    acu = 0. # result
    for i in range(n):
      x = np.cos(angle)*x0 + np.sin(angle)*y0
      y = np.cos(angle)*y0 - np.sin(angle)*x0
      arg = np.cos(i*np.pi*2/n)*x+np.sin(i*np.pi*2/n)*y
      acu += f(k*arg) 
    return v*acu/n
  return fun




def aahf1d(n0=0,beta=0.0000001,k=None,b=None,v=1.0,normalize=False):
  """Return the generalized AAHF potential"""
  tau = (1.+np.sqrt(5))/2.
  if b is None: b = 1/tau # default field
  if k is None: k = 3*np.pi*b # default phase
  if beta==0.0: beta=0.000001 # just in case
  def fun(r):
    """Function"""
    ns = r[0] # first coordinate
    ys = np.tanh(beta*(np.cos(2.*np.pi*b*ns+k)-np.cos(np.pi*b)))
    ys /= np.tanh(beta)
    return v*ys
  return fun # return function



def commensurate_potential(g):
    """Return a potential that is commensurate with
    the lattice"""
    a12 = g.a2.dot(g.a1)/(np.sqrt(g.a1.dot(g.a1))*np.sqrt(g.a1.dot(g.a1)))
    print(a12)
    if 0.49<abs(a12)<0.51: # angle is 60 degrees
      angle = np.pi/3.
      return cnpot(n=6,k=2.*np.pi/np.sqrt(g.a1.dot(g.a1)),angle=angle)
