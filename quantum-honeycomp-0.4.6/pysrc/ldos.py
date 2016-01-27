import scipy.sparse.linalg as slg
from scipy.sparse import csc_matrix as csc
from scipy.sparse import csc_matrix 
from scipy.sparse import bmat
import numpy as np

def ldos0d(h,e=0.0,delta=0.01):
  """Calculates the local density of states of a hamiltonian and
     writes it in file"""
  if h.dimensionality==0:  # only for 0d
    iden = np.identity(h.intra.shape[0],dtype=np.complex) # create identity
    g = ( (e+1j*delta)*iden -h.intra ).I # calculate green function
  else: raise # not implemented...
  d = [ -(g[i,i]).imag for i in range(len(g))] # get imaginary part
  d = spatial_dos(h,d) # convert to spatial resolved DOS
  g = h.geometry  # store geometry
  write_ldos(g.x,g.y,d) # write in file
  return d





def ldos0d_wf(h,e=0.0,delta=0.01,num_wf = 10):
  """Calculates the local density of states of a hamiltonian and
     writes it in file, using arpack"""
  if h.dimensionality==0:  # only for 0d
    intra = csc_matrix(h.intra) # matrix
  else: raise # not implemented...
  eig,eigvec = slg.eigsh(intra,k=int(num_wf),which="LM",sigma=e) # number of wfs
  d = np.array([0.0 for i in range(intra.shape[0])]) # initialize
  for (v,ie) in zip(eigvec.transpose(),eig): # loop over wavefunctions
    v2 = (np.conjugate(v)*v).real # square of wavefunction
    fac = delta/((e-ie)**2 + delta**2) # factor to create a delta
    d += fac*v2 # add contribution
  d /= num_wf # normalize
  d = spatial_dos(h,d) # resum if necessary
  g = h.geometry  # store geometry
  write_ldos(g.x,g.y,d) # write in file










def ldos1d(h,e=0.0,delta=0.001,nrep=3):
  """ Calculate DOS for a 1d system"""
  import green
  if h.dimensionality!=1: raise # only for 1d
  gb,gs = green.green_renormalization(h.intra,h.inter,energy=e,delta=delta)
  d = [ -(gb[i,i]).imag for i in range(len(gb))] # get imaginary part
  d = spatial_dos(h,d) # convert to spatial resolved DOS
  g = h.geometry  # store geometry
  x,y = g.x,g.y # get the coordinates
  go = h.geometry.copy() # copy geometry
  go = go.supercell(nrep) # create supercell
  write_ldos(go.x,go.y,d.tolist()*nrep) # write in file
  return d




def ldos2d(h,e=0.0,delta=0.001,nrep=3):
  """ Calculate DOS for a 1d system"""
  import green
  if h.dimensionality!=2: raise # only for 1d
  gb,gs = green.bloch_selfenergy(h,energy=e,delta=delta,mode="adaptive")
  d = [ -(gb[i,i]).imag for i in range(len(gb))] # get imaginary part
  d = spatial_dos(h,d) # convert to spatial resolved DOS
  g = h.geometry  # store geometry
  x,y = g.x,g.y # get the coordinates
  go = h.geometry.copy() # copy geometry
  go = go.supercell(nrep) # create supercell
  write_ldos(go.x,go.y,d.tolist()*(nrep**2)) # write in file












def spatial_dos(h,dos):
  """Resums a certain DOS to show only the spatial dependence"""
  if h.has_spin == False and h.has_eh==False: return np.array(dos)
  elif h.has_spin == True and h.has_eh==False: 
    return np.array([dos[2*i]+dos[2*i+1] for i in range(len(dos)/2)])
  elif h.has_spin == False and h.has_eh==True: 
    return np.array([dos[2*i]+dos[2*i+1] for i in range(len(dos)/2)])
  elif h.has_spin == True and h.has_eh==True: 
    return np.array([dos[4*i]+dos[4*i+1]+dos[4*i+2]+dos[4*i+3] for i in range(len(dos)/4)])
  else: raise


def write_ldos(x,y,dos):
  """ Write LDOS in a file"""
  fd = open("LDOS.OUT","w")   # open file
  fd.write("# x,  y, local density of states\n")
  for (ix,iy,idos) in zip(x,y,dos): # write everything
    fd.write(str(ix) +"   "+ str(iy) + "   "+ str(idos)+"\n")
  fd.close() # close file




def ldos_finite(h,e=0.0,n=10,nwf=4,delta=0.0001):
  """Calculate the density of states for a finite system"""
  if h.dimensionality!=1: raise # if it is not one dimensional
  intra = csc(h.intra) # convert to sparse
  inter = csc(h.inter) # convert to sparse
  interH = inter.H # hermitian
  m = [[None for i in range(n)] for j in range(n)] # full matrix
  for i in range(n): # add intracell
    m[i][i] = intra
  for i in range(n-1): # add intercell
    m[i][i+1] = inter
    m[i+1][i] = interH
  m = bmat(m) # convert to matrix
  (ene,wfs) = slg.eigsh(m,k=nwf,which="LM",sigma=0.0) # diagonalize
  wfs = wfs.transpose() # transpose wavefunctions
  dos = (wfs[0].real)*0.0 # calculate dos
  for (ie,f) in zip(ene,wfs): # loop over waves
    c = 1./(1.+((ie-e)/delta)**2) # calculate coefficient
    dos += np.abs(f)*c # add contribution
  odos = spatial_dos(h,dos) # get the spatial distribution
  go = h.geometry.supercell(n) # get the supercell
  write_ldos(go.x,go.y,odos) # write in a file
  return dos # return the dos





def ldos_defect(h,v,e=0.0,delta=0.001,n=1):
  """Calculates the LDOS of a cell with a defect, writting the n
  neighring cells"""
  raise # still not finished
  import green
  # number of repetitions
  rep = 2*n +1
  # calculate pristine green function
  g,selfe = green.supercell_selfenergy(h,e=e,delta=delta,nk=100,nsuper=rep)
  # now calculate defected green function 
  ez = e + 1j*delta # complex energy
  emat = np.matrix(np.identity(len(g)))*ez  # E +i\delta 
  import supercell
  pintra = supercell.intra_super2d(h,n=rep) # pristine
  vintra = supercell.intra_super2d(h,n=rep,central=v) # defective
  selfe = emat - pintra - g.I # dyson euqation, get selfenergy
  gv = (emat - vintra -selfe).I   # Green function of a vacancy, with selfener
  return



