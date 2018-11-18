
# special band structures

from __future__ import print_function
import topology
import operators
import scipy.linalg as lg
from scipy.sparse import csc_matrix
import scipy.sparse.linalg as slg
import numpy as np
import timing
import klist

from limits import densedimension as maxdim
arpack_tol = 1e-6
arpack_maxiter = 10000


def berry_bands(h,klist=None,mode=None,operator=None):
  """Calculate band structure resolved with Berry curvature"""
  ks = [] # list of kpoints
  if mode is not None: # get the mode
    if mode=="sz": operator = operators.get_sz(h)
    else: raise

  fo = open("BANDS.OUT","w")
  for ik in range(len(klist)): # loop over kpoints
    (es,bs) = topology.operator_berry_bands(h,k=klist[ik],operator=operator)
    for (e,b) in zip(es,bs):
      fo.write(str(ik)+"    "+str(e)+"    "+str(b)+"\n")
  fo.close()



def current_bands(h,klist=None):
  """Calcualte the band structure, with the bands"""
  if h.dimensionality != 1: raise # only 1 dimensional
  # go for the rest
  hkgen = h.get_hk_gen() # get generator of the hamiltonian
  if klist is None:  klist = np.linspace(0,1.,100) # generate k points
  fo = open("BANDS.OUT","w") # output file
  import current
  fj = current.current_operator(h) # function that generates the operator
  for k in klist: # loop over kpoints
    hk = hkgen(k) # get k-hamiltonian
    jk = fj(k) # get current operator
    evals,evecs = lg.eigh(hk) # eigenvectors and eigenvalues
    evecs = np.transpose(evecs) # transpose eigenvectors
    for (e,w) in zip(evals,evecs): # do the loop
        waw = braket_wAw(w,jk).real # product
        fo.write(str(k)+"    "+str(e)+"   "+str(waw)+"\n")
  fo.close()








def braket_wAw(w,A):
  w = np.matrix(w) # convert to matrix
  return ((w.T).H*A*w.T)[0,0] # expectation value


def ket_Aw(A,w):
  w = np.matrix([w]).T # transpose()
  wo = A*w # get a vector
  wo = np.array(wo)[:,0] # wo
  return wo


#def get_bands_0d(h,operator=None):
#  """ Returns a figure with the bandstructure of the system"""
#  if h.is_sparse:
#    energies = lg.eigvalsh(h.intra.todense()) # eigenvalues
#  else:
#    if operator is None: energies = lg.eigvalsh(h.intra) # eigenvalues
#    else: # matrix
#      energies,ws = lg.eigh(h.intra)
#      vals = []
#      ws = np.transpose(ws) # transpose
#      for w in ws: # loop over waves
#        if callable(operator):
#          waw = operator(w)
#        else:
#          waw = braket_wAw(w,operator).real # calcualte expectation value
#        vals.append(waw) # store
#  klist = np.linspace(0,1,len(energies))
#  if operator is None: msave = np.matrix([klist,energies]).T
#  else: msave = np.matrix([klist,energies,vals]).T
#  np.savetxt("BANDS.OUT",msave) # save all
#  return np.genfromtxt("BANDS.OUT").transpose() # return data







#def get_bands_1d(h,nkpoints=100,operator=None,num_bands=None,callback=None):
#  if num_bands is None: # all the bands
#    if operator is not None: diagf = lg.eigh # all eigenvals and eigenfuncs
#    else: diagf = lg.eigvalsh # all eigenvals and eigenfuncs
#  else: # using arpack
#    h = h.copy()
#    h.turn_sparse() # sparse Hamiltonian
#    def diagf(m):
#      eig,eigvec = slg.eigsh(m,k=num_bands,which="LM",sigma=0.0,
#                               tol=arpack_tol,maxiter=arpack_maxiter)
#      if operator is None: return eig
#      else: return (eig,eigvec)
#  hkgen = h.get_hk_gen() # generator hamiltonian
#  ks = np.linspace(0.,1.,nkpoints,endpoint=True)
#  f = open("BANDS.OUT","w") # open bands file
#  f.write("# system_dimension = 1\n")
##  if operator is not None: operator=np.matrix(operator) # convert to matrix
#  tr = timing.Testimator("BANDSTRUCTURE") # generate object
#  ik = 0 
#  for k in ks: # loop over kpoints
#    ik += 1
#    tr.remaining(ik,ks.shape[0])
#    hk = hkgen(k) # get hamiltonian
##    if h.is_sparse: hk = hk.todense() # turn the matrix dense
#    if operator is None:
#      es = diagf(hk)
#      for e in es:  # loop over energies
#        f.write(str(k)+"   "+str(e)+"\n") # write in file
#      if callback is not None: callback(k,es) # call the function
#    else:
#      es,ws = diagf(hk)
#      ws = ws.transpose() # transpose eigenvectors
#      for (e,w) in zip(es,ws):  # loop over waves
#        if callable(operator):
#          waw = operator(w,k=k)
#        else:
#          waw = braket_wAw(w,operator).real # calcualte expectation value
#        f.write(str(k)+"   "+str(e)+"  "+str(waw)+"\n") # write in file
#      if callback is not None: callback(k,es,ws) # call the function
#  f.close()
#  return np.genfromtxt("BANDS.OUT").transpose() # return data



def get_bands_nd(h,kpath=None,operator=None,num_bands=None,
                    callback=None,central_energy=0.0):
  """Get a 2d bandstructure"""
  if type(operator)==str: operator = self.get_operator(operator)
  if num_bands is None: # all the bands
    if operator is not None: 
      def diagf(m): # diagonalization routine
        if h.is_sparse and h.intra.shape[0]<maxdim: 
          return lg.eigh(m.todense()) # all eigenvals and eigenfuncs
        else:
          return lg.eigh(m) # all eigenvals and eigenfuncs
    else: 
      def diagf(m): # diagonalization routine
        if h.is_sparse and h.intra.shape[0]<maxdim: 
          return lg.eigvalsh(m.todense()) # all eigenvals and eigenfuncs
        else:
          return lg.eigvalsh(m) # all eigenvals and eigenfuncs
  else: # using arpack
    h = h.copy()
    h.turn_sparse() # sparse Hamiltonian
    def diagf(m):
      eig,eigvec = slg.eigsh(m,k=num_bands,which="LM",sigma=central_energy,
                                  tol=arpack_tol,maxiter=arpack_maxiter)
      if operator is None: return eig
      else: return (eig,eigvec)
  # open file and get generator
  f = open("BANDS.OUT","w") # open bands file
  hkgen = h.get_hk_gen() # generator hamiltonian
  if kpath is None:
    import klist
    kpath = klist.default(h.geometry) # generate default klist
  tr = timing.Testimator("BANDSTRUCTURE") # generate object
  for k in range(len(kpath)): # loop over kpoints
#    print("Bands in kpoint",k,"of",len(kpath),end="\r")
    tr.remaining(k,len(kpath))
    hk = hkgen(kpath[k]) # get hamiltonian
    if operator is None:
      es = diagf(hk)
      for e in es:  # loop over energies
        f.write(str(k)+"   "+str(e)+"\n") # write in file
      if callback is not None: callback(k,es) # call the function
    else:
      es,ws = diagf(hk)
      ws = ws.transpose() # transpose eigenvectors
      for (e,w) in zip(es,ws):  # loop over waves
        if callable(operator):  
          try: waw = operator(w,k=kpath[k]) # call the operator
          except: 
            print("Check out the k optional argument in operator")
            waw = operator(w) # call the operator
        else: waw = braket_wAw(w,operator).real # calculate expectation value
        f.write(str(k)+"   "+str(e)+"  "+str(waw)+"\n") # write in file
      # callback function in each iteration
      if callback is not None: callback(k,es,ws) # call the function
    f.flush()
  f.close()
  print("\nBANDS finished")
  return np.genfromtxt("BANDS.OUT").transpose() # return data



def smalleig(m,numw=10,evecs=False):
  """Return the smallest eigenvalues using arpack"""
  tol = arpack_tol
  eig,eigvec = slg.eigsh(m,k=numw,which="LM",sigma=0.0,
                                  tol=tol,maxiter=arpack_maxiter)
  if evecs:  return eig,eigvec.transpose()  # return eigenvectors
  else:  return eig  # return eigenvalues


def lowest_bands(h,nkpoints=100,nbands=10,operator = None,
                   info = False,kpath = None,discard=None):
  """ Returns a figure with the bandstructure of the system"""
  from scipy.sparse import csc_matrix
  if kpath is None: 
    k = klist.default(h.geometry) # default path
  else: k = kpath
  import gc # garbage collector
  fo = open("BANDS.OUT","w")
  if operator is None: # if there is not an operator
    if h.dimensionality==0:  # dot
      eig,eigvec = slg.eigsh(csc_matrix(h.intra),k=nbands,which="LM",sigma=0.0,
                                  tol=arpack_tol,maxiter=arpack_maxiter)
      eigvec = eigvec.transpose() # transpose
      iw = 0
      for i in range(len(eig)):
        if discard is not None: # use the function
          v = eigvec[i] # eigenfunction
          if discard(v): continue
        fo.write(str(iw)+"     "+str(eig[i])+"\n")
        iw += 1 # increase counter
    elif h.dimensionality>0: 
      hkgen = h.get_hk_gen() # get generator
      for ik in k:  # ribbon
        hk = hkgen(ik) # get hamiltonians
        gc.collect() # clean memory
        eig,eigvec = slg.eigsh(hk,k=nbands,which="LM",sigma=0.0)
        del eigvec # clean eigenvectors
        del hk # clean hk
        for e in eig:
          fo.write(str(ik)+"     "+str(e)+"\n")
        if info:  print("Done",ik,end="\r")
    else: raise # ups
  else:  # if there is an operator
    if h.dimensionality==1:
      hkgen = h.get_hk_gen() # get generator
      for ik in k:
        hk = hkgen(ik) # get hamiltonians
        eig,eigvec = slg.eigsh(hk,k=nbands,which="LM",sigma=0.0)
        eigvec = eigvec.transpose() # tranpose the matrix
        if info:  print("Done",ik,end="\r")
        for (e,v) in zip(eig,eigvec): # loop over eigenvectors
          a = braket_wAw(v,operator)
          fo.write(str(ik)+"     "+str(e)+"     "+str(a)+"\n")
  fo.close()



