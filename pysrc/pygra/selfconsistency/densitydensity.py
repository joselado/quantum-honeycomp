# specialized routine to perform an SCF, taking as starting point an
# attractive local interaction in a spinless Hamiltonian

from .. import inout
import numpy as np
import time
import os
from .. import densitymatrix
from copy import deepcopy
from numba import jit
from .. import utilities

class Interaction():
    def __init__(self,h=None):
        self.dimensionality = 0
        if h is not None: self.dimensionality = h.dimensionality
        self.v_dict = dict() # store dictionary
    def __mult__(self,a):
        """Function to multiply"""
        out = 0
        for key in self: out = out + self[key]*a[key]
        return out



def normal_term(v,dm):
    """Return the normal term of the mean field"""
    out = dm*0.0 # initialize
    return normal_term_jit(v,dm,out) # return the normal term



def normal_term_ii(v,dm):
    """Return the normal term of the mean field"""
    out = dm*0.0 # initialize
    return normal_term_ii_jit(v,dm,out) # return the normal term


def normal_term_jj(v,dm):
    """Return the normal term of the mean field"""
    out = dm*0.0 # initialize
    return normal_term_jj_jit(v,dm,out) # return the normal term


def normal_term_ij(v,dm):
    """Return the normal term of the mean field"""
    out = dm*0.0 # initialize
    return normal_term_ij_jit(v,dm,out) # return the normal term


def normal_term_ji(v,dm):
    """Return the normal term of the mean field"""
    out = dm*0.0 # initialize
    return normal_term_ji_jit(v,dm,out) # return the normal term

@jit(nopython=True)
def normal_term_jit(v,dm,out):
    """Return the normal terms, jit function"""
    n = len(v[0])
    for i in range(n): # loop
      for j in range(n): # loop
        out[i,j] = out[i,j] - v[i,j]*dm[j,i]
        out[j,i] = out[j,i] - v[i,j]*dm[i,j]
        out[i,i] = out[i,i] + v[i,j]*dm[j,j]
        out[j,j] = out[j,j] + v[i,j]*dm[i,i]
    return out


@jit(nopython=True)
def normal_term_ii_jit(v,dm,out):
    """Return the normal terms, jit function"""
    n = len(v[0])
    for i in range(n): # loop
      for j in range(n): # loop
        out[i,i] = out[i,i] + v[i,j]*dm[j,j]
    return out

@jit(nopython=True)
def normal_term_jj_jit(v,dm,out):
    """Return the normal terms, jit function"""
    n = len(v[0])
    for i in range(n): # loop
      for j in range(n): # loop
        out[j,j] = out[j,j] + v[i,j]*dm[i,i]
    return out


@jit(nopython=True)
def normal_term_ij_jit(v,dm,out):
    """Return the normal terms, jit function"""
    n = len(v[0])
    for i in range(n): # loop
      for j in range(n): # loop
        out[i,j] = out[i,j] - v[i,j]*dm[j,i]
    return out


@jit(nopython=True)
def normal_term_ji_jit(v,dm,out):
    """Return the normal terms, jit function"""
    n = len(v[0])
    for i in range(n): # loop
      for j in range(n): # loop
        out[j,i] = out[j,i] - v[i,j]*dm[i,j]
    return out



def update_hamiltonian(tdict,mf):
    """Update the hoppings with the mean field"""
    out = deepcopy(tdict) # copy
    for key in mf:
        out[key] = tdict[key] + mf[key] # add contribution
    return out # return dictionary


def mix_mf(mf,mf0,mix=0.8):
    """Mix mean fields"""
    out = dict() # initialize
    for key in mf: # loop
        if key not in mf0: out[key] = mf[key]
        else:
            #v0 = np.tanh(mf0[key]/mix)
            #v1 = np.tanh(mf[key]/mix)
            #out[key] = np.arctanh((v0+v1)/2.)*mix
            out[key] = mf0[key]*(1.-mix) + mf[key]*mix # add contribution
        #out[key] = mf0[key]*(1.-mix) + mf[key]*mix # add contribution
    return out



def diff_mf(mf0,mf):
    """Difference mean fields"""
    out = 0.0 # initialize
    for key in mf: # loop
        if key not in mf0: out += np.mean(np.abs(mf[key]))
        else: out += np.mean(np.abs(mf0[key] - mf[key])) # add contribution
        #out += np.mean(np.abs(mf0[key] - mf[key])) # add contribution
    return out # return


def hamiltonian2dict(h):
    out = dict() # create dictionary
    if not h.is_multicell: raise
    out[(0,0,0)] = h.intra
    for t in h.hopping: out[tuple(t.dir)] = t.m # store
    return out


def set_hoppings(h,hop):
    """Add the hoppings to the Hamiltonian"""
    h.is_multicell = True
    from ..multicell import Hopping as Hop
#    h.intra = h.intra*0.0 # set to zero
    h.intra = hop[(0,0,0)]
    hopping = [] # empty list
    for key in hop: # loop
        if key==(0,0,0): continue
        t = Hop() # generate
        t.dir = np.array(key) # transform to array
        t.m = hop[key] # matrix
        hopping.append(t) # store
    h.hopping = hopping # store

def get_dm(h,nk=1):
    """Get the density matrix"""
    ds = [(0,0,0)] # directions
    if h.dimensionality>0:
      for t in h.hopping: ds.append(tuple(t.dir)) # store
    dms = densitymatrix.full_dm(h,ds=ds,nk=nk) # get all the density matrices
    dm = dict()
    for i in range(len(ds)): 
        dm[ds[i]] = dms[i] # store
    return dm # return dictionary with the density matrix



def get_mf(v,dm,has_eh=False,**kwargs):
    """Get the mean field matrix"""
    if has_eh:
        # let us assume that it is a full Nambu spinor
        # (this may not be general, but good enough in the meantime)
        from .. import superconductivity
        dme = dict() # dictionary
        dma01 = dict() # dictionary
        dma10 = dict() # dictionary
        op = superconductivity.get_nambu2signless(dm[(0,0,0)]) # transform
        for key in dm: # extract the electron part 
            m = op.T@dm[key]@op # transform to the new basis
            dme[key] = superconductivity.get_eh_sector(m,i=0,j=0)
            dma01[key] = superconductivity.get_eh_sector(m,i=0,j=1)
            dma10[key] = superconductivity.get_eh_sector(m,i=1,j=0)
        mfe = get_mf_normal(v,dme,**kwargs) # electron part of the mean field
        # anomalous part
        mfa01 = get_mf_normal(v,dma01,compute_dd=False,add_dagger=False) 
        mfa10 = get_mf_normal(v,dma10,compute_dd=False,add_dagger=False) 
        # now rebuild the Hamiltonian
        mf = dict()
        for key in v:
            m = superconductivity.build_nambu_matrix(mfe[key],
                    c12 = -mfa10[key],c21=-mfa01[key]
                    )
            m = op.T@m@op # undo the transformation
            mf[key] = m # store this matrix
        return mf # return mean field matrix
    else: return get_mf_normal(v,dm,**kwargs) # no BdG Hamiltonian




def get_mf_normal(v,dm,compute_dd=True,add_dagger=True,
        compute_cross=True):
    """Get the mean field"""
    zero = dm[(0,0,0)]*0. # zero
    mf = dict()
    for d in v: mf[d] = zero.copy()  # initialize
    # compute the contribution to the mean field
    # onsite term
#    mf[(0,0,0)] = normal_term(v[(0,0,0)],dm[(0,0,0)]) 
    def dag(m): return m.T.conjugate()
    for d in v: # loop over directions
        d2 = (-d[0],-d[1],-d[2]) # minus this direction
        # add the normal terms
        if compute_dd: # only density density terms
            m = normal_term_ij(v[d],dm[d2]) # get matrix
            mf[d] = mf[d] + m # add normal term
            if add_dagger:
                mf[d2] = mf[d2] + dag(m) # add normal term
        if compute_cross: # density density terms
            m = normal_term_ii(v[d],dm[(0,0,0)]) # get matrix
            mf[(0,0,0)] = mf[(0,0,0)] + m # add normal term
            m = normal_term_jj(v[d2],dm[(0,0,0)]) # get matrix
            mf[(0,0,0)] = mf[(0,0,0)] + m # add normal term
    return mf

def get_dc_energy(v,dm):
    """Compute double counting energy"""
    out = 0.0
    for d in v: # loop over interactions
        d2 = (-d[0],-d[1],-d[2]) # minus this direction
        n = v[d].shape[0] # shape
        for i in range(n): # loop
          for j in range(n): # loop
              out -= v[d][i,j]*dm[(0,0,0)][i,i]*dm[(0,0,0)][j,j]
              c = dm[d][i,j] # cross term
              out += v[d][i,j]*c*np.conjugate(c) # add contribution
    print("DC energy",out.real)
    return out.real


def obj2mf(a):
    if type(a)==np.ndarray or type(a)==np.matrix:
        return {(0,0,0):a}
    else: return a


mf_file = "MF.pkl" 

def generic_densitydensity(h0,mf=None,mix=0.9,v=None,nk=8,solver="plain",
        maxerror=1e-5,filling=None,callback_mf=None,callback_dm=None,
        load_mf=True,compute_cross=True,
        callback_h=None,**kwargs):
    """Perform the SCF mean field"""
#    if not h0.check_mode("spinless"): raise # sanity check
    mf = obj2mf(mf)
    h1 = h0.copy() # initial Hamiltonian
    h1.turn_dense()
    h1.nk = nk # store the number of kpoints
    if mf is None:
      try: 
          if load_mf: mf = inout.load(mf_file) # load the file
          else: raise
      except: 
          mf = dict()
          for d in v: mf[d] = np.exp(1j*np.random.random(h1.intra.shape))
          mf[(0,0,0)] = mf[(0,0,0)] + mf[(0,0,0)].T.conjugate()
    else: pass # initial guess
    ii = 0
    os.system("rm -f STOP") # remove stop file
    hop0 = hamiltonian2dict(h1) # create dictionary
    def f(mf,h=h1):
      """Function to minimize"""
#      print("Iteration #",ii) # Iteration
      mf0 = deepcopy(mf) # copy
      h = h1.copy()
      hop = update_hamiltonian(hop0,mf) # add the mean field to the Hamiltonian
      set_hoppings(h,hop) # set the new hoppings in the Hamiltonian
      if callback_h is not None:
          h = callback_h(h) # callback for the Hamiltonian
      t0 = time.perf_counter() # time
      dm = get_dm(h,nk=nk) # get the density matrix
      if callback_dm is not None:
          dm = callback_dm(dm) # callback for the density matrix
      t1 = time.perf_counter() # time
      mf = get_mf(v,dm,compute_cross=compute_cross,
              has_eh=h0.has_eh) # return the mean field
      if callback_mf is not None:
          mf = callback_mf(mf) # callback for the mean field
      t2 = time.perf_counter() # time
      print("Time in density matrix = ",t1-t0) # Difference
      print("Time in the normal term = ",t2-t1) # Difference
      scf = SCF() # create object
      scf.hamiltonian = h # store
      scf.mf = mf # store mean field
      if os.path.exists("STOP"): scf.mf = mf0 # use the guess
      scf.dm = dm # store density matrix
      scf.v = v # store interaction
      return scf
    if solver=="plain":
      do_scf = True
      while do_scf:
        scf = f(mf) # new vector
        mfnew = scf.mf # new vector
        t0 = time.clock() # time
        diff = diff_mf(mfnew,mf) # mix mean field
        mf = mix_mf(mfnew,mf,mix=mix) # mix mean field
        t1 = time.clock() # time
        print("Time in mixing",t1-t0)
        print("ERROR",diff)
        #print("Mixing",dmix)
        print()
        if diff<maxerror: 
            inout.save(scf.mf,mf_file) # save the mean field
            return scf
    else: # use different solvers
        scf = f(mf) # perform one iteration
        fmf2a = get_mf2array(scf) # convert MF to array
        fa2mf = get_array2mf(scf) # convert array to MF
        def fsol(x): # define the function to solve
            mf1 = fa2mf(x) # convert to a MF
            scf1 = f(mf1) # compute function
            xn = fmf2a(scf1.mf) # new vector
            diff = x - xn # difference vector
            print("ERROR",np.max(np.abs(diff)))
            print()
            return x - xn # return vector
        x0 = fmf2a(scf.mf) # initial guess
        # these methods do seem too efficient, but lets have them anyway
        if solver=="krylov":
            from scipy.optimize import newton_krylov
            x = newton_krylov(fsol,x0,rdiff=1e-3) # use the solver
        elif solver=="anderson":
            from scipy.optimize import anderson
            x = anderson(fsol,x0) # use the solver
        elif solver=="broyden1":
            from scipy.optimize import broyden1
            x = broyden1(fsol,x0,f_tol=maxerror*100) # use the solver
        elif solver=="linear":
            from scipy.optimize import linearmixing
            x = linearmixing(fsol,x0,f_tol=maxerror*100) # use the solver
        else: raise # unrecognised solver
        mf = fa2mf(x) # transform to MF
        scf = f(mf) # compute the SCF with the solution
        inout.save(scf.mf,mf_file) # save the mean field
        return scf # return the mean field


def get_mf2array(scf):
    """Function to transform the mean field in an array"""
    nt = len(scf.mf) # number of terms in the dictionary
    n = scf.mf[(0,0,0)].shape[0]
    def fmf2a(mf):
        print(mf[(0,0,0)].real)
        out = [mf[key].real for key in mf] # to plain array
        out += [mf[key].imag for key in mf] # to plain array
        out = np.array(out)
        print(out.shape)
        out = out.reshape(nt*n*n*2) # reshape
        return out
    return fmf2a # return function

def get_array2mf(scf):
    """Function to transform an array into a mean field"""
    ds = [key for key in scf.mf] # store keys
    nt = len(scf.mf) # number of terms in the dictionary
    n = scf.mf[(0,0,0)].shape[0] # size
    def fa2mf(a):
        a = a.copy().reshape((2*nt,n*n)) # reshape array
        mf =  dict()
        print(a.shape)
        for i in range(len(ds)):
            d = ds[i]
            m = a[i,:] + 1j*a[i+nt,:] # get matrix
            mf[d] = m.reshape((n,n)) # store
        return mf
    return fa2mf # return function


def densitydensity(h,filling=0.5,**kwargs):
    """Function for density-density interactions"""
    if h.has_eh: 
        if not h.has_spin: return NotImplemented # only for spinful
    h = h.get_multicell()
    h.turn_dense()
    def callback_h(h):
        """Set the filling"""
        fermi = h.get_fermi4filling(filling,nk=h.nk) # get the filling
        print("Fermi energy",fermi)
        h.fermi = fermi
        h.shift_fermi(-fermi) # shift by the fermi energy
        return h
#    callback_h = None
    scf = generic_densitydensity(h,callback_h=callback_h,**kwargs)
    # Now compute the total energy
    h = scf.hamiltonian
    etot = h.get_total_energy(nk=h.nk)
    etot += h.fermi*h.intra.shape[0]*filling # add the Fermi energy
    print("Occupied energies",etot)
    etot += get_dc_energy(scf.v,scf.dm) # add the double counting energy
    etot = etot.real
    scf.total_energy = etot
    print("##################")
    print("Total energy",etot)
    print("##################")
    return scf




def hubbard(h,U=1.0,**kwargs):
    """Wrapper to perform a Hubbard model calculation"""
    h = h.copy() # copy Hamiltonian
    h.turn_multicell() # multicell Hamiltonian
    U = utilities.obj2fun(U) # redefine as a function
    if h.has_spin:
      n = len(h.geometry.r) # number of spinless sites
      zero = np.zeros((2*n,2*n),dtype=np.complex)
      for i in range(n): zero[2*i,2*i+1] = U(i) # Hubbard interaction
    else: 
      zero = np.zeros((n,n),dtype=np.complex)
      n = len(h.geometry.r) # number of spinless sites
      for i in range(n): zero[i,i] = U(i) # Hubbard interaction
    v = dict() # dictionary
    v[(0,0,0)] = zero 
    if h.has_spin:
      return densitydensity(h,v=v,**kwargs)
    else:
      return densitydensity(h,v=v,compute_cross=False,**kwargs)


def Vinteraction(h,V1=0.0,V2=0.0,U=0.0,**kwargs):
    """Wrapper to perform a Hubbard model calculation"""
    h = h.get_multicell() # multicell Hamiltonian
    h.turn_dense()
    # define the function
    nd = h.geometry.neighbor_distances() # distance to first neighbors
 #   def fun(r1,r2):
 #       dr = r1-r2
 #       dr = np.sqrt(dr.dot(dr)) # distance
 #       if abs(dr-nd[0])<1e-6: return V1/2.
 #       if abs(dr-nd[1])<1e-6: return V2/2.
 #       return 0.0
    from .. import specialhopping
    mgenerator = specialhopping.distance_hopping_matrix([V1/2.,V2/2.],nd[0:2])
    hv = h.geometry.get_hamiltonian(has_spin=False,is_multicell=True,
            mgenerator=mgenerator) 
 #   hv = h.geometry.get_hamiltonian(has_spin=False,is_multicell=True,
 #           fun=fun) 
    v = hv.get_hopping_dict() # hopping dictionary
    if h.has_spin: #raise # not implemented
        for d in v: # loop
            m = v[d] ; n = m.shape[0]
            m1 = np.zeros((2*n,2*n),dtype=np.complex)
            for i in range(n):
              for j in range(n): 
                  m1[2*i,2*j] = m[i,j]
                  m1[2*i+1,2*j] = m[i,j]
                  m1[2*i,2*j+1] = m[i,j]
                  m1[2*i+1,2*j+1] = m[i,j]
            v[d] = m1 # store
        for i in range(n):
            v[(0,0,0)][2*i,2*i+1] += U # add
        #    v[(0,0,0)][2*i+1,2*i] += U/2. # add
    return densitydensity(h,v=v,**kwargs)





class SCF(): pass

