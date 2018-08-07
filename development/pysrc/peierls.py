import numpy as np

from scipy.sparse import coo_matrix,csc_matrix


def add_phase(m1,r1,r2,phasefun,has_spin=False):
  if m1.shape[0] != len(r1): raise
  m = coo_matrix(m1) # convert to sparse matrix
  row,col = m.row,m.col
  data = m.data +0j
  for k in range(len(m.data)): # loop over non vanishing elements
    i = m.row[k]
    j = m.col[k]
    if has_spin: i,j = i/2,j/2 # if spinful
    # peierls phase
    p = phasefun(r1[i],r2[j]) # function yielding the phase
    data[k] *= p # add phase
  out = csc_matrix((data,(row,col)),shape=m1.shape) # convert to csc
  if type(m1) is type(np.matrix): out = out.todense() # dense matrix
  out = out.todense()
  return out






def add_peierls(h,mag_field=0.0,new=False):
  """ Adds Peierls phase to the Hamiltonian"""
  if h.has_eh: raise
  x = h.geometry.x    # x coordinate 
  y = h.geometry.y    # x coordinate 
  a1 = h.geometry.a1    # distance to neighboring cell
  celldis = np.sqrt(a1.dot(a1))
  from numpy import array
  norb = h.intra.shape[0]  # number of orbitals

  if new:
    print("New method to calculate peierls")
    g = h.geometry # geometry
    has_spin = h.has_spin # has spin degree of freedom
    def phasefun(ri,rj): 
      if callable(mag_field): return mag_field(ri[0],ri[1],rj[0],rj[1])
      else: raise
    h.intra = add_phase(h.intra,g.r,g.r,phasefun,has_spin) 
    if h.dimensionality==2:
      h.tx = add_phase(h.tx,g.r,g.replicas([1.,0.,0.]),phasefun,has_spin) 
      h.ty = add_phase(h.ty,g.r,g.replicas([0.,1.,0.]),phasefun,has_spin) 
      h.txy = add_phase(h.txy,g.r,g.replicas([1.,1.,0.]),phasefun,has_spin) 
      h.txmy = add_phase(h.txmy,g.r,g.replicas([1.,-1.,0.]),phasefun,has_spin) 
      return
  else: # old method


    if h.is_sparse: # sparse hamiltonian
      if True: # zero dimensional
        m = coo_matrix(h.intra) # convert to sparse matrix
        row,col = m.row,m.col
        data = m.data +0j
        for k in range(len(m.data)): # loop over non vanishing elements
          i = m.row[k]
          j = m.col[k]
          if h.has_spin: i,j = i/2,j/2 # raise if spinful
          p = peierls(x[i],y[i],x[j],y[j],mag_field) # peierls phase
          data[k] *= p # add phase
        h.intra = csc_matrix((data,(row,col)),shape=(norb,norb)) # convert to csc
      if h.dimensionality==1: # one dimensional
        # check that celldis is right
        if np.abs(celldis - h.geometry.a1[0])>0.001: raise
        def phaseize(inter,numn=1):
          m = coo_matrix(inter) # convert to sparse matrix
          row,col = m.row,m.col
          data = m.data +0j
          for k in range(len(m.data)): # loop over non vanishing elements
            i = m.row[k]
            j = m.col[k]
            if h.has_spin: i,j = i//2,j//2 # raise if spinful
            # peierls phase
            p = peierls(x[i],y[i],x[j]+numn*celldis,y[j],mag_field) 
            data[k] *= p # add phase
          return csc_matrix((data,(row,col)),shape=(norb,norb)) # convert to csc
        # for normal hamiltonians
        if not h.is_multicell: h.inter = phaseize(h.inter,numn=1) 
        # for multicell hamiltonians
        if h.is_multicell: 
          hopping = [] # empty list
          for i in range(len(h.hopping)):
            h.hopping[i].m = phaseize(h.hopping[i].m,numn=h.hopping[i].dir[0]) 
    else: # not sparse
      def gaugeize(m,d=0.0):
        """Add gauge phase to a matrix"""
        for i in range(len(x)):
          for j in range(len(x)):
            p = peierls(x[i],y[i],x[j]+d,y[j],mag_field) # peierls phase
            if h.has_spin:
              m[2*i,2*j] *= p
              m[2*i,2*j+1] *= p
              m[2*i+1,2*j] *= p
              m[2*i+1,2*j+1] *= p
            else:
              m[i,j] *= p
      gaugeize(h.intra,d=0.0)  # gaugeize intraterm
      if h.dimensionality==0: pass # if zero dimensional
      elif h.dimensionality==1: # if one dimensional
        if h.is_multicell: raise
        gaugeize(h.inter,d=celldis) # gaugeize interterm
      elif h.dimensionality==2: # if bigger dimensional
        print("WARNING, is your gauge periodic?")
        gaugeize(h.tx,d=h.geometry.a1[0]) # gaugeize interterm
      else:
        raise




def peierls(x1,y1,x2,y2,mag_field):
  """ Returns the complex phase with magnetic field """
  if is_number(mag_field): b = mag_field 
  elif callable(mag_field): 
      b = mag_field(np.array([x1,y1,0.0]),np.array([x2,y2,0.0]))
  else: raise
  phase = b*(x1-x2)*(y1+y2)/2.0
  return np.exp(1j*phase)



def is_number(s):
    try:
        float(s)
        return True
    except:
        return False

