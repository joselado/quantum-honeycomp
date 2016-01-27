import scipy.sparse.linalg as lg 
from scipy.sparse import csc_matrix as csc
from scipy.sparse import bmat
import numpy as np


class pair_mf:
  """ Classs for mean field matrices"""
  def __init__(self,a,b,g=0.0):
    self.a = a # B matrix
    self.b = b # A matrix
    self.g = g # coupling


def hubbard(h,U=1.0,pairs = []):
  """ Creates mean field matrices for the input hamiltonian,
      returns list with mean field amtrix pairs in csc form"""
  hubm = [] # list with the hubbard mean field matrices
  if len(pairs)==0:
    de = len(h.intra) # dimension of the hamiltonian
    pairs = range(de/2) # create list with oribtal indices
  for i in pairs: # loop over up down pairs
    # density density term    
    data = [1.0+0.j] # value of the term in the matrix
    ij = [[2*i],[2*i]] # indexes in the matrix
    a = csc((data,ij),shape=(de,de))
    ij = [[2*i+1],[2*i+1]] # indexes in the matrix
    b = csc((data,ij),shape=(de,de))
    hubm += [pair_mf(a,b,g=U)] # add mean field pair
    # exchange exchange term    
    data = [1.0+0.j] # value of the term in the matrix
    ij = [[2*i],[2*i+1]] # indexes in the matrix
    a = csc((data,ij),shape=(de,de))
    ij = [[2*i+1],[2*i]] # indexes in the matrix
    b = csc((data,ij),shape=(de,de))
    hubm += [pair_mf(a,b,g= -U)] # add mean field pair
  return hubm # return list with the pairs


def selfconsistency(h,ab_list=None,nkp = 100,filling=0.5,old_mf=None):
  """ Solve a selfocnsistent mean field one dimensional system"""
  from scipy.linalg import eigh
  from operator import itemgetter
  if ab_list==None: ab_list=hubbard(h)
  intra = h.intra  # intramatrix
  inter = h.inter  # intermatrix
  if old_mf==None:  old_mf = intra*0.
  while True: # infinite loop
    htmp = h.copy()  # copy hamiltonian
    htmp.intra += old_mf # add mean field 
    eigvals,eigvecs = htmp.eigenvectors(nkp) # get eigenvectors
    # get the fermi energy
    ne = len(eigvals) ; ifermi = int(round(ne*filling))
    fermi = sorted(eigvals)[ifermi]
    mf = ab_list[0].a*0. # initialize mean field
    for (e,v) in zip(eigvals,eigvecs): # loop over eigenvals,eigenvecs
      for iab in ab_list: # loop over mean field matrices
        if e<fermi:  # if level is filled, add contribution
          a = iab.a # A matrix
          b = iab.b # B matrix
          g = iab.g # coupling
          vav = (v.H * a * v).todense()[0,0] # <vAv>
          vbv = (v.H * b * v).todense()[0,0] # <vBv>
          mf = mf + (vav*b + vbv*a)*g # mean field hamiltonian 
    mf = mf.todense() # new intramatrix
    if h.dimensionality>0: mf = mf/float(nkp) # normalize by nkvectors
    error = np.max(np.abs(old_mf-mf)) # absolute difference
    print "Error in SCF =",error
    old_mf = old_mf*0.1 + mf*0.9 # mixing
    if error<0.001: # if converged break
      break
  return mf # return mean field

def ferro_initialization(h,value = 0.1):
  mf = h.intra*0.0
  for i in range(len(mf)/2):
    mf[2*i,2*i] = value
    mf[2*i+1,2*i+1] = -value
  return mf



def antiferro_initialization(h,value = 0.1):
  mf = h.intra*0.0
  for i in range(len(mf)/2):
    v = value*(-1)**i
    mf[2*i,2*i] = v
    mf[2*i+1,2*i+1] = -v
  return mf



def directional_abg(v = np.array([0.,0.,1.]),i=0,d=1,g=1.0):
  """Calculates a pair of mean field matrices for a particular site"""
  r = np.sqrt(v[0]**2 + v[1]**2)  # radial vector in xy
  theta = np.arctan2(r,v[2])  # calculate theta
  phi = np.arctan2(v[1],v[0])  # calculate phi
  t,p = theta, phi  # store in easiest variables
  print t,p
  ct,st = np.cos(t/2.0), np.sin(t/2.0)  # sin and cosine
  e, ce = np.exp(1j*p), np.exp(-1j*p)  # complex phase
  ij=[[2*i,2*i],[2*i,2*i+1],[2*i+1,2*i],[2*i+1,2*i+1]] # index
  row=[2*i,2*i,2*i+1,2*i+1] # index
  col=[2*i,2*i+1,2*i,2*i+1] # index
  data1=[ct*ct,ct*st*ce,ct*st*e,st*st] # values of the matrix A
  data2=[st*st,-ct*st*ce,-ct*st*e,ct*ct] # values of B
  a = csc((data1,(row,col)),shape=(d,d))   # create A
  b = csc((data2,(row,col)),shape=(d,d))   # create B
  abg = pair_mf(a,b,g=g) # create the A, B, coupling object
  return abg


def directional_hubbard(vecs,g=1.):
  """Calculates the mena field matrices for a collinear hubbard""" 
  abgs = []
  d = len(vecs) # number of sites
  i = 0
  print vecs
  for v in vecs: # loop over vectors
    abgs.append(directional_abg(v=v,i=i,d=2*d,g=g))
    i += 1
  return abgs






def write_ab_mf(abgs):
  """Write mean field operators in mean_field_operators.in,
     input is list of ab_pairs"""
  f=open('mean_field_operators.in','w')
  f.write('Number of matrices\n')
  f.write(str(len(abgs))+'\n')  # number of pairs
  def write_a(a,f):
      """Function to write a particular matrix"""
      f.write('    number of non vanishing elements\n') 
      f.write(str(len(a.data))+'\n') # number of nonvanishing
      rc,data = a.nonzero(), a.data # row, columns and data
      row,col = rc[0],rc[1] # store row and column
      for (i,j,d) in zip(row,col,data): # loop over element
        f.write(str(i+1)+"   "+str(j+1)+"   ")
        f.write('{0:.8f}'.format(d.real)+'   ')
        f.write('{0:.8f}'.format(d.imag)+'   ')
        f.write('\n')
  for (abg,i) in zip(abgs,range(len(abgs))): # loop over pairs
    a,b,g = abg.a , abg.b ,abg.g  # store in simplest variables
    f.write('Matrix A  '+str(i+1)+'\n')  # title
    write_a(a,f) # write matrix a
    f.write('Matrix B  '+str(i+1)+'\n')  # title
    write_a(b,f) # write matrix b
    f.write('Coupling  '+str(i+1)+'\n') # now coupling
    f.write(str(g.real) +"   "+str(g.imag)+"\n\n")  # write the value
  f.close() # close the file



def directional_mean_field(vecs):
  """ Creates an initial mean field accoring to certain vectors"""
  from hamiltonians import sx,sy,sz
  mf = [[None for i in vecs] for j in vecs] # create
  for i in range(len(vecs)):
    v = vecs[i]
    mf[i][i] = v[0]*sx + v[1]*sy + v[2]*sz # add contribution
  return bmat(mf).todense()



def initialization(g,scfin="Random"):
  import random
  def rv():
    return .5 - np.random.random(3) # 3 component random vector
  def rvxy():
    v = np.random.random(3) -.5 # 3 component random vector
    v[0] = 0.0
    return v
  ###################################
  if scfin == "Reconverge": # no new mean field
    pass
  else:
    if scfin == "Ferro":
      v = rv() # one random vector
      vecs = [v for i in range(len(g.x))] # create vectors
    elif scfin == "Ferro Z":
      v = np.array([0.,0.,1.]) # one random vector
      vecs = [v for i in range(len(g.x))] # create vectors
    elif scfin == "Ferro X":
      v = np.array([1.,0.,0.]) # one random vector
      vecs = [v for i in range(len(g.x))] # create vectors
    elif scfin == "Ferro Y":
      v = np.array([0.,1.,0.]) # one random vector
      vecs = [v for i in range(len(g.x))] # create vectors
    elif scfin == "Ferro XY":
      vecs = [rvxy() for i in range(len(g.x))] # create vectors
    elif scfin == "Antiferro":
      v = rv() # one random vector
      vecs = [v*g.sublattice[i] for i in range(len(g.x))] # create vectors
    elif scfin == "Antiferro Y":
      v = np.array([0.,1.,0.]) # one random vector
      vecs = [v*g.sublattice[i] for i in range(len(g.x))] # create vectors
    elif scfin == "Antiferro X":
      v = np.array([1.,0.,0.]) # one random vector
      vecs = [v*g.sublattice[i] for i in range(len(g.x))] # create vectors
    elif scfin == "Antiferro Z":
      v = np.array([0.,0.,1.]) # one random vector
      vecs = [v*g.sublattice[i] for i in range(len(g.x))] # create vectors
    elif scfin == "Random":
      vecs = [rv() for i in range(len(g.x))] # create vectors
    else:
      raise
    mf = directional_mean_field(vecs) # create mean field matrix
    return mf # return matrix


def tb90_scf(h,ab,old_mf,nkpoints=100,error=0.0000001,filling=.5,do_dos=False):
  """Run SCF cycle using the tb90 code"""
  h.write("hamiltonian_0.in") # write in file
  import input_tb90
  in90 = input_tb90.tb90in() # create class
  in90.mode("nothing") # no mode
  in90.mode("SCF") # no mode
  in90.kpoints.nkpoints = nkpoints
  if do_dos: in90.dos.do_dos = True # do DOS
  in90.electrons.filling = filling
  in90.scf_convergence.mean_field_operators = "from_file"
  in90.scf_convergence.mean_field_matrix = "from_file"
  in90.scf_convergence.max_scf_err = error # error in SCF
  in90.write() # write in file
  write_ab_mf(ab) # write men field operators 
  input_tb90.write_mean_field(old_mf,output_file="mean_field.in") # write old mean field
#  exit()
  import os
  os.system("tb90.x") # run SCF
  mf = input_tb90.read_mean_field(input_file="mean_field.in")
  return mf

