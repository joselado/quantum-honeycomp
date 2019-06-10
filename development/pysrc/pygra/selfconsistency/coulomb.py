import os
import numpy as np
from ..scftypes import scfclass
from ..scftypes import get_fermi_energy
from ..scftypes import get_occupied_states
import time
from .. import limits
from .. import inout



mf_file = "MF.pkl" # mean field file

def coulombscf(h,g=1.0,nkp = 100,filling=0.5,mix=0.9,
                  maxerror=1e-05,silent=False,mf=None,
                  smearing=None,fermi_shift=0.0,
                  maxite=1000,save=False,**kwargs):
  """ Solve a selfconsistent Hubbard mean field"""
  vc = coulomb_density_matrix(h.geometry,vc=2*g) # get the matrix
  mix = 1. - mix
  U = g # redefine
  if h.has_spin: raise # not implemented
  os.system("rm -f STOP") # remove stop file
  from scipy.linalg import eigh
  nat = h.intra.shape[0] # number of atoms
  htmp = h.copy()  # copy hamiltonian
  htmp.turn_dense() # turn into a dense Hamiltonian
  # generate the necessary list of correlators
  if mf is None: # generate initial mean field
    try:  
        old_mf = inout.load(mf_file) # load the file
        print("Mean field read from file")
    except: # generate mean field
        old_mf = np.random.random((nat,nat))
        old_mf += old_mf.T
  else: old_mf = mf # use guess
  # get the pairs for the correlators
  ndim = nat # dimension
  totkp = nkp**(h.dimensionality) # total number of kpoints
  file_etot = open("SCF_ENERGY.OUT","w")
  file_error = open("SCF_ERROR.OUT","w")
  ite = 0 # iteration counter
  scf = scfclass(h) # create scf class
  while True: # infinite loop
    ite += 1 # increase counter
    htmp.intra = h.intra + old_mf # add mean field 
    t1 = time.time()
# get eigenvectors
    eigvals,eigvecs,kvectors = htmp.eigenvectors(nkp,kpoints=True)
    eigvecs = np.conjugate(eigvecs)
# fermi energy
    t2 = time.time()
    fermi = get_fermi_energy(eigvals,filling,fermi_shift=fermi_shift)
# occupied states
    eoccs,voccs,koccs = get_occupied_states(eigvals,eigvecs,kvectors,fermi)
# mean field
    mf,edc = charge_mean_field(voccs,vc) # get the new mean field
    mf = mf/totkp # normalize
    edc = edc/totkp # normalize
    t3 = time.time()
    error = np.max(np.abs(old_mf-mf)) # absolute difference
    # total energy
    etot = np.sum(eoccs)/totkp + edc  # eigenvalues and double counting
    file_etot.write(str(ite)+"    "+str(etot)+"\n") # write energy in file
    file_error.write(str(ite)+"    "+str(error)+"\n") # write energy in file
    file_etot.flush()
    file_error.flush()
#    totcharge = np.sum(charge).real # total charge
#    avcharge = totcharge/nat # average charge
    if save: inout.save(mf,mf_file) # save the mean field
    ######
    if not silent:
      print("Times in diagonalization",t2-t1)
      print("Times in new mean field",t3-t2)
      print("\n")
      print("Iteration number =",ite)
      print("Error in SCF =",error)
      print("Fermi energy =",fermi)
      print("Total energy =",etot)
#      print("Total charge =",totcharge)
    old_mf = old_mf*mix + mf*(1.-mix) # mixing
    if error<maxerror or os.path.exists("STOP"): # if converged break
      break
    if ite>=maxite: break # if too many iterations
  file_etot.close() # close file
  file_error.close() # close file
  # output result
  class scftmp(): pass # create an empty class
  scf = scftmp() # empty class
  scf.hamiltonian = htmp.copy() # copy Hamiltonian
  scf.hamiltonian.intra -= fermi*np.identity(ndim) # shift Fermi energy
  scf.total_energy = etot # store total energy
  scf.mf = mf # store mean field matrix
  return scf # return mean field





def charge_mean_field(voccs,vc):
    """Return the charge mean field"""
    n = voccs[0].shape[0]
    charge = np.zeros(n) # initialize
    mf = np.zeros(n) # initialize
    for w in voccs:
        charge += np.abs(w)**2 # add contribution
    for i in range(len(vc)): # loop over rows of the matrix
        mf[i] = charge.dot(vc[i]) # get contribution
    # This formula is only ok if there is a single pair per ij
    edc = -mf@(vc@mf)
    return np.diag(mf),edc # return diagonal matrix



def coulomb_density_matrix(g,**kwargs):
    """Return a list with the Coulomb interaction terms"""
    from ..crystalfield import hartree_onsite
    m = hartree_onsite(g,**kwargs) # get array
    m = m - np.mean(m) # remove average
    nat = len(g.r) # number of atoms
    ind = range(nat) # indexes
    mat = csc_matrix((m,(ind,ind)),shape=(nat,nat),dtype=np.complex)
    return mat # return

