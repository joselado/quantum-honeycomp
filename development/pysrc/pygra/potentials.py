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
    else: raise






def fibonacci(n,n0=0):
  """Generate the Fibonacci sequence"""
  a1 = [0] # first
  a2 = [0,1] # first
  while True:
    a3 = a1 + a2
    a1 = a2
    a2 = a3
    if len(a3)>n+n0: break
  out = [a3[i] for i in range(n0,n0+n)] # output
  return np.array(out)




def thue_morse(n,n0=0):
  """Generate the Thue-Morse sequence"""
  out = []
  for i in range(n0,n0+n): # loop over integers
    ib = bin(i).split("b")[1] # in binary
    acu = 0 # start
    for iib in ib: # loop
      if iib=="1": acu += 1
    acu = acu%2 # modulus
    out.append(acu) # output
  return np.array(out)



def tbgAA(g):
    """Return a function that yields +1 for AA and -1 otherwise"""
    from . import geometry
    from .specialhopping import twisted_matrix
    h = g.get_hamiltonian(mgenerator=twisted_matrix(ti=0.3),
            has_spin=False,is_sparse=True)
    h.set_filling(0.5,nk=1)
    ks = [[.5,0.,0.],[0.,.5,0.],[.5,.5,0.]]
    (x,y,d) = h.get_ldos(e=0.0,delta=0.01,ks=ks,mode="arpack",
            write=False,silent=True)
    d = d - np.mean(d) # average in zero
    d = d - np.min(d)
    d = d/np.max(d) # set maximum in 1
    d = 2*(d-0.5) # between -1 and +1
    funr = geometry.get_fractional_function(g,center=True)
    rf = np.array([funr(ri) for ri in g.r])
    fint = interpolate2d(rf[:,0:2],d) # interpolation
    return lambda ri: fint(funr(ri))[0] 
#    def f(ri):
#        dr = rf[:,0:3] - funr(ri)[0:3]
#        dr = np.sum(np.abs(dr),axis=1)
#        i = np.argmin(dr)
#        return d[i] # this one
#    return f



def interpolate2d(r,v):
    """Return a function that does 2d interpolation"""
    from scipy.interpolate import interp2d
    x,y = r[:,0],r[:,1] # data
    grid_x, grid_y = np.mgrid[np.min(x):np.max(x):100j,np.min(y):np.max(y):100j]
    from scipy.interpolate import griddata
    grid_z = griddata(r,v, (grid_x, grid_y), method='nearest')
    f = interp2d(grid_x, grid_y, grid_z, kind='linear')
    return lambda ri: f(ri[0],ri[1])

