# functions to wrap around fortran code


import pylab as py
import numpy as np


class gf_convergence():
   """ Class to manage the convergence  options
   of the green functions """
   optimal = False
   refinement = False
   guess = True # use old green function
   def __init__(self,mode):
     if mode=="fast":   # fast mode,used for coule to finite systems 
       self.eps = 0.001
       self.max_error = 1.0
       self.num_rep = 10
       self.mixing = 1.0
     if mode=="lead":
       self.eps = 0.001
       self.max_error = 0.00001
       self.num_rep = 3
       self.mixing = 0.8
     if mode=="hundred":  
       self.eps = 0.001
       self.max_error = 1.0
       self.num_rep = 100
       self.mixing = 1.0



def dyson(intra,inter,energy=0.0,gf=None,is_sparse=False,initial = None):
  """ Solves the dyson equation for a one dimensional
  system with intra matrix 'intra' and inter to the nerest cell
  'inter'"""
  # get parameters
  if gf == None: gf = gf_convergence("lead")
  mixing = gf.mixing
  eps = gf.eps
  max_error = gf.max_error
  num_rep = gf.num_rep
  optimal = gf.optimal
  try:
    intra = intra.todense()
    inter = inter.todense()
  except:
    a = 1
  if initial==None:  # if green not provided. initialize at zero
    from numpy import zeros
   
    g_guess = intra*0.0j
  else:
    g_guess = initial
  # calculate using fortran
  if optimal:
    print "Fortran dyson calculation"
    from green_fortran import dyson  # import fortran subroutine
    (g,num_redo) = dyson(intra,inter,energy,num_rep,mixing=mixing,
               eps=eps,green_guess=g_guess,max_error=max_error)
    print "      Converged in ",num_redo,"iterations\n"
    from numpy import matrix
    g = matrix(g)
  # calculate using python
  if not optimal:
    g_old = g_guess # first iteration
    iden = np.matrix(np.identity(len(intra),dtype=complex)) # create identity
    e = iden*(energy+1j*eps) # complex energy
    while True: # loop over iterations
      self = inter*g_old*inter.H # selfenergy
      g = (e - intra - self).I # dyson equation
      if np.max(np.abs(g-g_old))<gf.max_error: break
#      print np.max(np.abs(g-g_old))
      g_old = mixing*g + (1.-mixing)*g_old # new green function
  if is_sparse: 
    from scipy.sparse import csc_matrix
    g = csc_matrix(g)
  return g











def dos_infinite(intra,inter,energies=[0.0],num_rep=100,
                      mixing=0.7,eps=0.0001,green_guess=None,max_error=0.0001):
   """ Calculates the surface density of states by using a 
    green function approach"""
   dos = [] # list with the density of states
   iden = np.matrix(np.identity(len(intra),dtype=complex)) # create idntity
   for energy in energies: # loop over energies
     # right green function
     gr = dyson(intra,inter,energy=energy,num_rep=num_rep,mixing=mixing,
          eps=eps,green_guess=green_guess,max_error=max_error)
     # left green function
     gl = dyson(intra,inter.H,energy=energy,num_rep=num_rep,mixing=mixing,
          eps=eps,green_guess=green_guess,max_error=max_error)
     # central green function
     selfl = inter.H*gl*inter # left selfenergy
     selfr = inter*gr*inter.H # right selfenergy
     gc = energy*iden -intra -selfl -selfr # dyson equation for the center
     gc = gc.I # calculate inverse
     dos.append(-gc.trace()[0,0].imag)  # calculate the trace of the Green function
   return dos




def dos_semiinfinite(intra,inter,energies=[0.0],num_rep=100,
                      mixing=0.7,eps=0.0001,green_guess=None,max_error=0.0001):
   """ Calculates the surface density of states by using a 
    green function approach"""
   dos = [] # list with the density of states
   for energy in energies: # loop over energies
     gf = dyson(intra,inter,energy=energy,num_rep=num_rep,mixing=mixing,
          eps=eps,green_guess=green_guess,max_error=max_error)
     dos.append(-gf.trace()[0,0].imag)  # calculate the trace of the Green function
   return dos







def plot_dos_semiinfinite(intra,inter,energies=[0.0],num_rep=100,
                      mixing=0.7,eps=0.0001,green_guess=None,max_error=0.0001):
   """ Plots the surface density of states by using a 
    green function approach"""
   # get the dos
   dos = dos_semiinfinite(intra,inter,energies=energies,
             num_rep=num_rep,mixing=mixing,
             eps=eps,green_guess=green_guess,max_error=max_error)
   # plot the figure
   fig = py.figure() # create figure
   fig.subplots_adjust(0.2,0.2)
   fig.set_facecolor("white") # face in white
   sdos = fig.add_subplot(111) # create subplot
   sdos.set_xlabel("Energy",size=20)
   sdos.set_ylabel("Surface DOS",size=20)
   sdos.plot(energies,dos,color="red") # create the plot
   sdos.plot(energies,dos,color="red") # create the plot
   sdos.fill_between(energies,0,dos,color="red")
   sdos.tick_params(labelsize=20)
   sdos.set_ylim([0.0,max(dos)])
   return fig



def plot_dos_infinite(intra,inter,energies=[0.0],num_rep=100,
                      mixing=0.7,eps=0.0001,green_guess=None,max_error=0.0001):
   """ Plots the density of states by using a 
    green function approach"""
   # get the dos
   dos = dos_infinite(intra,inter,energies=energies,
             num_rep=num_rep,mixing=mixing,
             eps=eps,green_guess=green_guess,max_error=max_error)
   # plot the figure
   fig = py.figure() # create figure
   fig.set_facecolor("white") # face in white
   sp = fig.add_subplot(111) # create subplot
   sp.set_xlabel("Energy",size=20)
   sp.set_ylabel("Surface DOS",size=20)
   sp.plot(energies,dos) # create the plot
   return fig






def dos_heterostructure(hetero,energies=[0.0],num_rep=100,
                      mixing=0.7,eps=0.0001,green_guess=None,max_error=0.0001):
   """ Calculates the density of states 
       of a heterostructure by a  
    green function approach, input is a heterostructure class"""
   dos = [] # list with the density of states
   iden = np.matrix(np.identity(len(intra),dtype=complex)) # create idntity
   for energy in energies: # loop over energies
     # right green function
     intra = hetero.right_intra
     inter = hetero.right_inter
     gr = dyson(intra,inter,energy=energy,num_rep=num_rep,mixing=mixing,
          eps=eps,green_guess=green_guess,max_error=max_error)
     # left green function
     intra = hetero.right_intra
     inter = hetero.right_inter
     gl = dyson(intra,inter,energy=energy,num_rep=num_rep,mixing=mixing,
          eps=eps,green_guess=green_guess,max_error=max_error)
     # central green function
     selfl = inter.H*gl*inter # left selfenergy
     selfr = inter*gr*inter.H # right selfenergy
     gc = energy*iden -intra -selfl -selfr # dyson equation for the center
     gc = gc.I # calculate inverse
     dos.append(-gc.trace()[0,0].imag)  # calculate the trace of the Green function
   return dos



def read_matrix(f):
  """Read green function from a file"""
  m = np.genfromtxt(f)
  d = int(max(m.transpose()[0]))+1 # dimension of the green functions
  print d
  g = np.matrix([[0.0j for i in range(d)] for j in range(d)]) # create matrix
  for r in m:
    i = int(r[0])
    j = int(r[1])
    ar = r[2]
    ai = r[3]
    g[i,j] = ar +1j*ai # store element
  return g # return green function



def write_matrix(f,g):
  """Write green function from a file"""
  fw = open(f,"w") # open file to write
  n = len(g) # dimension of the matrix
  for i in range(n):
    for j in range(n):
      fw.write(str(i)+"  ")
      fw.write(str(j)+"  ")
      fw.write(str(g[i,j].real)+"  ")
      fw.write(str(g[i,j].imag)+"\n")
  fw.close()   # close file




def write_sparse(f,g):
  """ Write a sparse matrix in a file"""
  from input_tb90 import nv_el
  fw = open(f,"w") # open the file
  fw.write("# dimension = "+str(g.shape[0])+"\n")
  nv=nv_el(g)
  for iv in range(len(nv)):
    fw.write(str(int(nv[iv][0]))+'   ')
    fw.write(str(int(nv[iv][1]))+'   ')
    fw.write('{0:.8f}'.format(float(nv[iv][2]))+'   ')
    fw.write('{0:.8f}'.format(float(nv[iv][3]))+'   ')
    fw.write('  !!!  i  j   Real   Imag\n')
  fw.close()




def read_sparse(f,sparse=True):
  """Read green function from a file"""
  l = open(f,"r").readlines()[0] # first line
  d = int(l.split("=")[1])
  m = np.genfromtxt(f)
  if not sparse:
# create matrix  
    g = np.matrix([[0.0j for i in range(d)] for j in range(d)]) 
    for r in m:
      i = int(r[0])-1
      j = int(r[1])-1
      ar = r[2]
      ai = r[3]
      g[i,j] = ar +1j*ai # store element
  if sparse:
    from scipy.sparse import coo_matrix
    g = coo_matrix([[0.0j for i in range(d)] for j in range(d)]) 
    row = np.array([0 for i in range(len(m))])
    col = np.array([0 for i in range(len(m))])
    data = np.array([0j for i in range(len(m))])
    for i in range(len(m)):
      r = m[i]
      row[i] = int(r[0])-1
      col[i] = int(r[1])-1
      ar = r[2]
      ai = r[3]
      data[i] = ar +1j*ai # store element
    g.col = col
    g.row = row
    g.data = data
  return g # return green function





def gauss_inverse(m,i=0,j=0):
  from gauss_inv import gauss_inv as ginv
  """ Calculates the inverso of a block diagonal
      matrix """
  nb = len(m) # number of blocks
  ca = [None for ii in range(nb)]
  ua = [None for ii in range(nb-1)]
  da = [None for ii in range(nb-1)]
  for ii in range(nb): # diagonal part
    ca[ii] = m[ii][ii]
  for ii in range(nb-1):
    ua[ii] = m[ii][ii+1]
    da[ii] = m[ii+1][ii]
  mout = ginv(ca,da,ua,i+1,j+1)
  return np.matrix(mout)






def green_renormalization(intra,inter,energy=0.0,nite=None,
                            error=0.0001,info=False,delta=0.001):
  """ Calculates bulk and surface Green function by a renormalization
  algorithm, as described in I. Phys. F: Met. Phys. 15 (1985) 851-858 """

  e = np.matrix(np.identity(len(intra))) * (energy + 1j*delta)
  ite = 0
  alpha = inter
  beta = inter.H
  epsilon = intra
  epsilon_s = intra + inter * (e-intra).I * inter.H
  epsilon_s = intra
  while True: # implementation of Eq 11
    einv = (e - epsilon).I # inverse
    epsilon_s = epsilon_s + alpha * einv * beta
    epsilon = epsilon + alpha * einv * beta + beta * einv * alpha
    alpha = alpha * einv * alpha  # new alpha
    beta = beta * einv * beta  # new beta
    ite += 1
    # stop conditions
    if not nite == None:
      if ite > nite:  break 
    else:
      if np.max(np.abs(alpha))<error and np.max(np.abs(beta))<error: break
  if info:
    print "Converged in ",ite,"iterations"
  g_surf = (e - epsilon_s).I # surface green function
  g_bulk = (e - epsilon).I  # bulk green function 
  return g_bulk,g_surf



def bloch_selfenergy(h,nk=100,energy = 0.0, delta = 0.01,mode="full",
                         error=0.00001):
  """ Calculates the selfenergy of a cell defect,
      input is a hamiltonian class"""
  if mode=="adaptative": mode = "adaptive"
  def gr(ons,hop):
    """ Calculates G by renormalization"""
    gf,sf = green_renormalization(ons,hop,energy=energy,nite=None,
                            error=error,info=False,delta=delta)
    return gf,sf
  hk_gen = h.get_hk_gen()  # generator of k dependent hamiltonian
  d = h.dimensionality # dimensionality of the system
  g = h.intra *0.0j # initialize green function
  e = np.matrix(np.identity(len(g)))*(energy + delta*1j) # complex energy
  if mode=="full":  # full integration
    if d==1: # one dimensional
      ks = np.linspace(0.,1.,nk,endpoint=False) 
    elif d==2: # two dimensional
      ks = []
      kk = np.linspace(0.,1.,nk,endpoint=False)  # interval 0,1
      for ikx in kk:
        for iky in kk:
          ks.append([ikx,iky])
      ks = np.array(ks)  # all the kpoints
    else: raise # raise error
    for k in ks:  # loop in BZ
      g += (e - hk_gen(k)).I  # add green function  
    g = g/len(ks)  # normalize
  #####################################################
  #####################################################
  if mode=="renormalization":
    if d==1: # full renormalization
      g,s = gr(h.intra,h.inter)  # perform renormalization
    elif d==2: # two dimensional, loop over k's
      ks = np.linspace(0.,1.,nk,endpoint=False) 
      for k in ks:  # loop over k in y direction
 # add contribution to green function
        g += green_kchain(h,k=k,energy=energy,delta=delta,error=error) 
      g = g/len(ks)
  #####################################################
  #####################################################
  if mode=="adaptive":
    if d==1: # full renormalization
      g,s = gr(h.intra,h.inter)  # perform renormalization
    elif d==2: # two dimensional, loop over k's
      ks = np.linspace(0.,1.,nk,endpoint=False) 
      import integration
      def fint(k):
        """ Function to integrate """
        return green_kchain(h,k=k,energy=energy,delta=delta,error=error) 
      # eps is error, might work....
      g = integration.integrate_matrix(fint,xlim=[0.,1.],eps=error) 
        # chain in the y direction
    else: raise
  # now calculate selfenergy
  selfenergy = e - h.intra - g.I
  return g,selfenergy



def green_kchain(h,k=0.,energy=0.,delta=0.01,only_bulk=True,error=0.0001):
  """ Calculates the green function of a kdependent chain for a 2d system """
  def gr(ons,hop):
    """ Calculates G by renormalization"""
    gf,sf = green_renormalization(ons,hop,energy=energy,nite=None,
                            error=error,info=False,delta=delta)
    if only_bulk:  return gf
    else:  return gf,sf
  tky = h.ty*np.exp(1j*np.pi*2.*k)
  tkx = h.ty*np.exp(1j*np.pi*2.*k)
  tkxy = h.txy*np.exp(1j*np.pi*2.*k)
  tkxmy = h.txmy*np.exp(-1j*np.pi*2.*k)  # notice the minus sign !!!!
  # chain in the x direction
  ons = h.intra + tky + tky.H  # intra of k dependent chain
  hop = h.tx + tkxy + tkxmy  # hopping of k-dependent chain
  return gr(ons,hop)  # return green function



def supercell_selfenergy(h,e=0.0,delta=0.001,nk=100,nsuper=[1,1]):
  """alculates the selfenergy of a certain supercell """
  try:   # if two number given
    nsuper1 = nsuper[0]
    nsuper2 = nsuper[1]
  except: # if only one number given
    nsuper1 = nsuper
    nsuper2 = nsuper
  print "Supercell",nsuper1,"x",nsuper2
  ez = e + 1j*delta # create complex energy
  import dyson2d
  g = dyson2d.dyson2d(h.intra,h.tx,h.ty,h.txy,h.txmy,nsuper1,nsuper2,300,ez)
  g = np.matrix(g) # convert to matrix
  n = nsuper1*nsuper2 # number of cells
  intrasuper = [[None for j in range(n)] for i in range(n)]
  # create indexes (same order as in fortran routine)
  k = 0
  inds = []
  for i in range(nsuper1):
    for j in range(nsuper2):
      inds += [(i,j)]
      k += 1 
  # create hamiltonian of the supercell
  from scipy.sparse import bmat
  from scipy.sparse import csc_matrix as csc
  tx = csc(h.tx)
  ty = csc(h.ty)
  txy = csc(h.txy)
  txmy = csc(h.txmy)
  intra = csc(h.intra)
  for i in range(n):
    intrasuper[i][i] = intra # intracell
    (x1,y1) = inds[i]
    for j in range(n):
      (x2,y2) = inds[j]
      dx = x2-x1
      dy = y2-y1
      if dx==1 and  dy==0: intrasuper[i][j] = tx  
      if dx==-1 and dy==0: intrasuper[i][j] = tx.H  
      if dx==0 and  dy==1: intrasuper[i][j] = ty  
      if dx==0 and  dy==-1: intrasuper[i][j] = ty.H  
      if dx==1 and  dy==1: intrasuper[i][j] = txy 
      if dx==-1 and dy==-1: intrasuper[i][j] = txy.H 
      if dx==1 and  dy==-1: intrasuper[i][j] = txmy  
      if dx==-1 and dy==1: intrasuper[i][j] = txmy.H  
  intrasuper = bmat(intrasuper).todense() # supercell
  eop = np.matrix(np.identity(len(g),dtype=np.complex))*(ez)
  selfe = eop - intrasuper - g.I
  return g,selfe
