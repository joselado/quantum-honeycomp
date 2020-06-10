import numpy as np
from scipy.sparse import csc_matrix
from numba import jit

try:
    from . import specialhoppingf90
    use_fortran=True
except:
    use_fortran=False
    print("FORTRAN not working in specialhopping")



def twisted(cutoff=5.0,ti=0.3,lambi=8.0,
        lamb=12.0,dl=3.0,lambz=10.0,b=0.0,phi=0.0):
  """Hopping for twisted bilayer graphene"""
  cutoff2 = cutoff**2 # cutoff in distance
  def fun(r1,r2):
    rr = (r1-r2) # distance
    rr = rr.dot(rr) # distance
    if rr>cutoff2: return 0.0 # too far
    if rr<0.001: return 0.0 # same atom
    dx = r1[0]-r2[0]
    dy = r1[1]-r2[1]
    dz = r1[2]-r2[2]
    r = np.sqrt(rr)
#    if r2>100.0: return 0.0 # too far
    if (r-1.0)<-0.1: 
      print(r)
      raise
    out = -(dx*dx + dy*dy)/rr*np.exp(-lamb*(r-1.0))*np.exp(-lambz*dz*dz)
    out += -ti*(dz*dz)/rr*np.exp(-lambi*(r-dl))
    #### fix for magnetic field
    cphi = np.cos(phi*np.pi)
    sphi = np.sin(phi*np.pi)
    r = (r1+r2)/2.
    dr = r1-r2
    p = 2*r[2]*(dr[0]*sphi - dr[1]*cphi)
    out *= np.exp(1j*b*p)
    #####
    return out
  return fun


def twisted_matrix(cutoff=5.0,ti=0.3,lambi=8.0,
        lamb=12.0,dl=3.0,lambz=10.0,**kwargs):
  """Function capable of returning the hopping matrix
  for twisted bilayer graphene"""
  if use_fortran:
    from . import specialhoppingf90
    def funhop(r1,r2):
      """Function that returns a hopping matrix"""
      nr = len(r1) # 
      nmax = len(r1)*int(10*cutoff**2) # maximum number of hoppings
      (ii,jj,ts,nout) = specialhoppingf90.twistedhopping(r1,r2,nmax,
                                  cutoff,ti,lamb,lambi,lambz,1e-5,dl)
      if nout>nmax: raise # sanity check
      ts = ts[0:nout]
      ii = ii[0:nout]
      jj = jj[0:nout]
      out = csc_matrix((ts,(ii-1,jj-1)),shape=(nr,nr),dtype=np.complex) # matrix
      return out
  else:
    print("Using Python function in twisted")
    def funhop(r1,r2):
      fh = twisted(cutoff=cutoff,ti=ti,lambi=lambi,lamb=lamb,dl=dl,**kwargs)
      m = np.array([[fh(r1i,r2j) for r1i in r1] for r2j in r2],dtype=np.complex)
      m = csc_matrix(m,dtype=np.complex).T
      m.eliminate_zeros()
      return m
#      raise
  return funhop # function




def multilayer(ti=0.3,dz=3.0):
    """Return hopping for a multilayer"""
    def fhop(ri,rj):
      """Function to compute the hopping"""
      dr = ri-rj ; dr2 = dr.dot(dr) # distance
      if abs(1.0-dr2)<0.01: return 1.0 # first neighbors
      # interlayer hopping (distance between the layers is 3)
      if abs(dz**2-dr2)<0.01 and abs(dz-abs(dr[2]))<0.01: return ti
      return 0.0 # else
    return fhop






def phase_C3(g,phi=0.0):
    """Create a fucntion that computes hoppings that alternate
    between +\phi and -\phi every 60 degrees"""
    if len(g.r)==1:
      g = g.supercell(2) # create a supercell
    ds = g.get_connections()
    i = 0 # first site
    j = ds[i][0] # connected site
    dr = g.r[i] - g.r[j] # distance between sites
    z = dr[0] + 1j*dr[1] # comple vector
    z0 = np.exp(1j*np.pi*2./3.) # rotation
    zs = [z,z*z0,z*z0**2] # the three vectors
    def fun(r1,r2):
        """Function to compute hoppings"""
        dr = r1-r2
        if 0.99<dr.dot(dr)<1.01: # first neighbors
            zi = dr[0]+1j*dr[1]
            for zj in zs: # one of the three directions
                d = np.abs(zi/zj-1.0)
                print(d)
                if d<1e-2: 
                    return np.exp(1j*phi)
            else: return np.exp(-1j*phi)
        return 0.0
    return fun






def distance_hopping_matrix(vs,ds):
    """Return a hopping that to the 1-th neighbor is vs"""
    vs = np.array(vs)
    ds = np.array(ds)
    def mgenerator(r1,r2):
        r1 = np.array(r1)
        r2 = np.array(r2)
        n = len(r1)
        out = np.zeros((n,n),dtype=np.complex) # output
        return distance_hopping_matrix_jit(r1,r2,vs,ds**2,out) 
    return mgenerator

@jit(nopython=True)
def distance_hopping_matrix_jit(r1,r2,vs,ds2,out):
    """Return a hopping that to the 1-th neighbor is vs"""
    n = len(r1) # number of sites
    nn = len(ds2) # number of neighbors
    for i in range(n):
      for j in range(n):
          dr = r1[i] - r2[j] # difference
          dr2 = dr[0]*dr[0] + dr[1]*dr[1] + dr[2]*dr[2]
          for k in range(nn):
              if np.abs(ds2[k]-dr2)<1e-6: out[i,j] = vs[k]
    return out






