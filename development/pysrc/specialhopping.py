import numpy as np
from scipy.sparse import csc_matrix


def twisted(cutoff=5.0,ti=0.3,lambi=8.0,lamb=12.0):
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
    out = -(dx*dx + dy*dy)/rr*np.exp(-lamb*(r-1.0))
    out += -ti*(dz*dz)/rr*np.exp(-lambi*(r-3.0))
    return out
  return fun


def twisted_matrix(cutoff=5.0,ti=0.3,lambi=8.0,lamb=12.0):
  """Function capable of returning the hopping matrix
  for twisted bilayer graphene"""
  try: # use fortran routine
    import specialhoppingf90
    def funhop(r1,r2):
      """Function that returns a hopping matrix"""
      nr = len(r1) # 
      nmax = len(r1)*100 # maximum number of hoppings
      (ii,jj,ts,nout) = specialhoppingf90.twistedhopping(r1,r2,nmax,
                                  cutoff,ti,lamb,lambi,1e-3)
      ts = -ts[0:nout]
      ii = ii[0:nout]
      jj = jj[0:nout]
      out = csc_matrix((ts,(ii-1,jj-1)),shape=(nr,nr),dtype=np.complex) # matrix
      return out
  except:
    print("FORTRAN not working in specialhopping")
    def funhop(r1,r2):
      fh = twisted(cutoff=cutoff,ti=ti,lambi=lambi,lamb=lamb)
      m = np.matrix([[fh(r1i,r2j) for r1i in r1] for r2j in r2])
      return csc_matrix(m,dtype=np.complex)
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

