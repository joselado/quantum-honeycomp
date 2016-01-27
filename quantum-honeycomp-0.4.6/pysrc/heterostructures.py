
import numpy as np
import pylab as py
from copy import deepcopy as dc


class heterostructure():
  """Class for a heterostructure"""
  file_right_green = "green_right.dat"  # in/out file for right green
  file_left_green = "green_left.dat"   # in/out file for left green
  file_heff = "heff.dat"   # out for the effective hamiltonian
  is_sparse = False
  block_diagonal = False
  def __init__(self,h):  # initialization using a hamiltonian
    self.heff = None  # effective hamiltonian
    self.right_intra = h.intra  # intraterm in the right lead
    self.right_inter = h.inter  # interterm in the right lead (to the right)
    self.right_green = None  # right green function
    self.left_intra = h.intra  # intraterm in the left lead
    self.left_inter = h.inter.H  # interterm in the left lead (to the left)
    self.left_green = None  # left green function
    self.central_intra = h.intra  # intraterm in the center
    self.right_coupling = h.inter # coupling from the center to the right lead
    self.left_coupling = h.inter.H # coupling from the center to the left lead
    # geometry of the central part
    gc = dc(h.geometry)
#    gc.set_finite()
    self.central_geometry = gc # geometry of the central part
    # additional degrees of freedom
    self.has_spin = h.has_spin   # spin degree of freedom
    self.has_eh = h.has_eh   # electron hole pairs
  def plot_central_dos(self,energies=[0.0],num_rep=100,
                      mixing=0.7,eps=0.0001,green_guess=None,max_error=0.0001):
    """ Plots the density of states by using a 
    green function approach"""
    return plot_central_dos(self,energies=energies,num_rep=num_rep,
                      mixing=mixing,eps=eps,green_guess=green_guess,
                      max_error=max_error)
  def plot_landauer(self,energy=[0.0],delta=0.001,has_eh=False):
    """ Plots the density of states by using a 
    green function approach"""
    return plot_landauer(self,energy=energy,delta=delta,has_eh=has_eh)
  def landauer(self,energy=[0.],delta=0.0001,do_leads=True):
    """ Return the Landauer transmission"""
    return landauer(self,energy=energy,delta=delta,do_leads=do_leads)
  def plot_local_central_dos(self,energies=0.0,gf=None):
    """ Plots the local densityt of states of the central part
        at a given energy"""
    return plot_local_central_dos(self,energies=energies,gf=gf)
  def write_green(self):
    """Writes the green functions in a file"""
    from green import write_matrix
    write_matrix(self.file_right_green,self.right_green)
    write_matrix(self.file_left_green,self.left_green)
  def read_green(self):
    """Reads the green functions from a file"""
    from green import read_matrix
    self.right_green = read_matrix(self.file_right_green)
    self.left_green = read_matrix(self.file_left_green)
  def write_heff(self):
    """ Writes effective hamiltonian in a file"""
    from green import write_sparse
    write_sparse(self.file_heff,self.heff)
  def eigenvalues(self,numeig = 10,effective=False,full=False):
    """ Calculates eigenvalues of the central part """
    return eigenvalues(self,numeig=numeig,effective=effective,
                        full=full)
  def replace_center(self,ht_replacement):
    """ Replaces the central part by the second argument"""
    self.central_intra = ht_replacement.central_intra  # make the change





def create_leads_and_central(h_right,h_left,h_central,num_central=1,
      interpolation="None",block_diagonal=True):
  """ Creates an heterojunction by giving the hamiltonian
     of the leads and the center, with a certain number of cells
     in the center  """
  # check the hamiltonians
  h_right.check()
  h_left.check()
  h_central.check()
  ht = heterostructure(h_central) # create heterostructure
  # assign matrices of the leads
  ht.right_intra = h_right.intra 
  ht.right_inter = h_right.inter 
  ht.left_intra = h_left.intra 
  ht.left_inter = h_left.inter.H
  # create matrix of the central part and couplings to the leads
  from scipy.sparse import csc_matrix,bmat
  z = csc_matrix(h_central.intra*0.0j) # zero matrix
  hc = [[None for i in range(num_central)] for j in range(num_central)]
  tcr = [[None] for i in range(num_central)]
  tcl = [[None] for i in range(num_central)]  
  ##############################################################
  # create the centrla hamiltonian according to different schemes
  ##############################################################
  if not block_diagonal: 
    ec = csc_matrix(h_central.intra) 
    er = csc_matrix(h_right.intra) 
    el = csc_matrix(h_left.intra) 
    tr = csc_matrix(h_right.inter) 
    tl = csc_matrix(h_right.inter) 
    tc = csc_matrix(h_central.inter) 
  if block_diagonal: 
    ec = h_central.intra 
    er = h_right.intra 
    el = h_left.intra 
    tr = h_right.inter 
    tl = h_right.inter
    tc = h_central.inter
  # central part is pure central input hamilotnian
  if interpolation=="None": # without central interpolation
    for i in range(num_central):  # intra term of the central blocks
      hc[i][i] = ec
    for i in range(num_central-1): # interterm of the central blocks
      hc[i][i+1] = tc
      hc[i+1][i] = tc.H

  # central part is a step of the right and left hamiltonians
  elif interpolation=="step": # with step interpolation
    # intraterm
    for i in range(num_central):  # intra term of the central blocks
      if i<num_central/2:
          hc[i][i] = el
      elif i>num_central/2:
          hc[i][i] = er
      elif i==num_central/2:
          hc[i][i] = ec
      else:
        raise
    # interterm
    for i in range(num_central-1): # interterm of the central blocks
      if i<num_central/2:
        hc[i][i+1] = tl
        hc[i+1][i] = tl.H
      elif i>num_central/2:
        hc[i][i+1] = tr
        hc[i+1][i] = tr.H
      elif i==num_central/2:
        hc[i][i+1] = tc
        hc[i+1][i] = tc.H
      else:
        raise


  # central part is a linear interpolation of right and left
  elif interpolation=="linear": 
    for i in range(num_central):  # intra term of the central blocks
      r = float(i+1)/float(num_central+1)   # 
      hc[i][i] = er*r + el*(1.-r)  # mean interpolation
    for i in range(num_central-1): # interterm of the central blocks
      r = float(i+1)/float(num_central+1)   # 
      hc[i][i+1] = tr*r + tl*(1.-r)
      hc[i+1][i] = (tr*r + tl*(1.-r)).H


  else:
    raise

  for i in range(num_central):  # intra term of the central blocks
    tcr[i][0] = csc_matrix(z) 
    tcl[i][0] = csc_matrix(z) 

  tcr[-1][0] = csc_matrix(h_right.inter) # hopping to the right lead
  tcl[0][0] = csc_matrix(h_left.inter.H) # hopping to the left lead
  # create dense matrices
  if not block_diagonal:
    hc = bmat(hc).todense()
    tcr = bmat(tcr).todense() 
    tcl = bmat(tcl).todense()
    ht.block_diagonal = False
  # do not create dense matrices
  if block_diagonal:
    hc = hc   # this is a list !!!!!
    tcr = h_right.inter   # this is a matrix
    tcl = h_left.inter.H  # this is a matrix
    ht.block_diagonal = True
  # assign to the heterostructure
  ht.right_coupling = tcr
  ht.left_coupling = tcl
  ht.central_intra = hc
  # and modify the geometry
  ht.central_geometry.supercell(num_central) 
  return ht



def central_dos(hetero,energies=[0.0],num_rep=100,
                      mixing=0.7,eps=0.0001,green_guess=None,max_error=0.0001):
   """ Calculates the density of states 
       of a heterostructure by a  
    green function approach, input is a heterostructure class"""
   from green import dyson
   dos = [] # list with the density of states
   intra = hetero.central_intra # central intraterm
   iden = np.matrix(np.identity(len(intra),dtype=complex)) # create idntity
   for energy in energies: # loop over energies
     # right green function
     intra = hetero.right_intra
     inter = hetero.right_inter
     gr = dyson(intra,inter,energy=energy,num_rep=num_rep,mixing=mixing,
          eps=eps,green_guess=green_guess,max_error=max_error)
     # left green function
     intra = hetero.left_intra
     inter = hetero.left_inter
     gl = dyson(intra,inter,energy=energy,num_rep=num_rep,mixing=mixing,
          eps=eps,green_guess=green_guess,max_error=max_error)
     # left selfenergy
     inter = hetero.left_coupling
     selfl = inter*gl*inter.H # left selfenergy
     # right selfenergy
     inter = hetero.right_coupling
     selfr = inter*gr*inter.H # right selfenergy
     # central green function
     intra = hetero.central_intra
     gc = energy*iden -intra -selfl -selfr # dyson equation for the center
     gc = gc.I # calculate inverse
     dos.append(-gc.trace()[0,0].imag)  # calculate the trace of the Green function
   return dos





def plot_central_dos(ht,energies=[0.0],num_rep=100,
                      mixing=0.7,eps=0.0001,green_guess=None,max_error=0.0001):
   """ Plots the density of states by using a 
    green function approach"""
   # get the dos
   dos = central_dos(ht,energies=energies,
             num_rep=num_rep,mixing=mixing,
             eps=eps,green_guess=green_guess,max_error=max_error)
   # plot the figure
   fig = py.figure() # create figure
   fig.set_facecolor("white") # face in white
   sp = fig.add_subplot(111) # create subplot
   sp.set_xlabel("Energy",size=20)
   sp.set_ylabel("Central DOS",size=20)
   sp.plot(energies,dos) # create the plot
   return fig





def landauer(hetero,energy=0.0,delta = 0.001,error=0.00001,do_leads=True,
             gr=None,gl=None,has_eh=False):
   """ Calculates transmission using Landauer formula"""
   try: # if it is a list
     return [landauer(hetero,energy=e,delta=delta,error=error,
                      do_leads=do_leads,gr=gr,gl=gl,has_eh=has_eh) for e in energy]
   except: # contnue if it is not a list
     pass
   import green
   from hamiltonians import is_number
   if not hetero.block_diagonal:
     intra = hetero.central_intra # central intraterm   
     dimhc = len(intra) # dimension of the central part
   if hetero.block_diagonal:
     intra = hetero.central_intra[0][0] # when it is diagonal
# dimension of the central part
     dimhc = len(hetero.central_intra)*len(intra)  
   iden = np.matrix(np.identity(len(intra),dtype=complex)) # create idntity
   intra = hetero.right_intra
   inter = hetero.right_inter
   if do_leads:
     grb,gr = green.green_renormalization(intra,inter,error=error,
                                          energy=energy,delta=delta)
     hetero.right_green = gr # save green function
   else:
     gr = hetero.right_green # get the saved green function
   # left green function
   intra = hetero.left_intra
   inter = hetero.left_inter
   if do_leads:
     glb,gl = green.green_renormalization(intra,inter,error=error,
                                          energy=energy,delta=delta)
     hetero.left_green = gl # save green function
   else:
     gl = hetero.left_green # get the saved green function
   # left selfenergy
   inter = hetero.left_coupling
   selfl = inter*gl*inter.H # left selfenergy
   # right selfenergy
   inter = hetero.right_coupling
   selfr = inter*gr*inter.H # right selfenergy

   #################################
   # calculate spectral functions
   #################################
   gammar = selfr-selfr.H
   gammal = selfl-selfl.H

   #################################
   # dyson equation for the center     
   #################################
   # central green function
   intra = hetero.central_intra
   # full matrix
   if not hetero.block_diagonal:
     heff = intra + selfl + selfr
     hetero.heff = heff
     gc = (energy+1j*delta)*iden - heff
     gc = gc.I # calculate inverse
     print has_eh
     if has_eh: # if it has electron-hole, trace over electrons
       raise
       G = (gammar*gc*gammal.H*gc.H)
       G = np.sum([G[2*i,2*i] for i in range(len(G)/2)]).real
     else:
       G = (gammar*gc*gammal.H*gc.H).trace()[0,0].real
   # reduced matrix
   if hetero.block_diagonal: 
     print has_eh
     from copy import deepcopy
     heff = deepcopy(intra)
     heff[0][0] = intra[0][0] + selfl
     heff[-1][-1] = intra[-1][-1] + selfr
     dd = (energy+1j*delta)*iden
     for i in range(len(intra)):  # add the diagonal energy part
       heff[i][i] = heff[i][i] - dd  # this has the wrong sign!!
    # now change the sign
     for i in range(len(intra)):  
       for j in range(len(intra)):  
         try:
           heff[i][j] = -heff[i][j]
         except:
           heff[i][j] = heff[i][j]
     # calculate green function
     from green import gauss_inverse  # routine to invert the matrix
     # calculate only some elements of the central green function
     gc1n = gauss_inverse(heff,0,len(heff)-1) # calculate element 1,n      
     # and apply Landauer formula
     if has_eh: # if it has electron-hole, trace over electrons
       raise
       G = (gammal*gc1n*gammar.H*gc1n.H)
       G = np.sum([G[2*i,2*i] for i in range(len(G)/2)]).real
     else:
       G = (gammal*gc1n*gammar.H*gc1n.H).trace()[0,0].real
   print "Landauer transmission E=",energy,"G=",G
   return G # return transmission










def plot_landauer(ht,energy=[0.0],delta=0.001,has_eh=False):
   """ Plots the density of states and Landauer transmission
    by using a 
    green function approach"""
   # get the transmission
   trans = landauer(ht,energy=energy,delta=delta,has_eh=has_eh)
   # plot the figure
   fig = py.figure() # create figure
   fig.subplots_adjust(0.2,0.2)
   fig.set_facecolor("white") # face in white
#   sdos = fig.add_subplot(121) # create subplot for the DOS
#   strans = fig.add_subplot(122) # create subplot for the transmission
   strans = fig.add_subplot(111) # create subplot for the transmission

   # DOS graph
#   sdos.set_xlabel("Energy",size=20)
#   sdos.set_ylabel("Central DOS",size=20)
#   sdos.plot(energies,dos,color="red") # create the plot
#   sdos.fill_between(energies,0,dos,color="red")
#   sdos.tick_params(labelsize=20)
#   sdos.set_ylim([0.0,max(dos)])
#   sdos.set_xlim([min(energies),max(energies)])

   # transport graph
   strans.set_xlabel("Energy",size=20)
   strans.set_ylabel("Transmission",size=20)
   strans.plot(energy,trans,color="green") # create the plot
   strans.fill_between(energy,0,trans,color="green")
   strans.set_ylim([0.0,max(trans)+0.5])
   strans.set_xlim([min(energy),max(energy)])
   strans.tick_params(labelsize=20)
   return fig




def plot_local_central_dos(hetero,energies=0.0,gf=None):
   """ Plots the local density of states in the central part"""
   from green import dyson
   from green import gf_convergence
   from hamiltonians import is_number
   if gf==None: gf = gf_convergence("lead")
   if not hetero.block_diagonal:
     intra = hetero.central_intra # central intraterm   
     dimhc = len(intra) # dimension of the central part
   if hetero.block_diagonal:
     intra = hetero.central_intra[0][0] # when it is diagonal
# dimension of the central part
     dimhc = len(hetero.central_intra)*len(intra)  
   iden = np.matrix(np.identity(len(intra),dtype=complex)) # create idntity
   ldos = np.array([0.0 for i in range(dimhc)]) # initialice ldos
   # initialize ldos
   for energy in energies: # loop over energies
     # right green function
     gr = None
     gl = None
     # perform dyson calculation

     intra = hetero.right_intra
     inter = hetero.right_inter
     gr = dyson(intra,inter,energy=energy,gf=gf)
     hetero.right_green = gr # save green function
     # left green function
     intra = hetero.left_intra
     inter = hetero.left_inter
     gl = dyson(intra,inter,energy=energy,gf=gf)
     hetero.left_green = gl # save green function
     # save green functions
#     hetero.write_green()
#     print "Saved green functions"
     # left selfenergy
     inter = hetero.left_coupling
     selfl = inter*gl*inter.H # left selfenergy
     # right selfenergy
     inter = hetero.right_coupling
     selfr = inter*gr*inter.H # right selfenergy
     # central green function
     intra = hetero.central_intra
# dyson equation for the center     
     # full matrix
     if not hetero.block_diagonal:
       heff = intra + selfl + selfr
       hetero.heff = heff
       gc = (energy+1j*eps)*iden - heff
       gc = gc.I # calculate inverse
       # get the local density of states
       ldos += np.array([-gc[i,i].imag for i in range(len(gc))])
#       if save_heff:
#         print "Saving effective hamiltonian in ",hetero.file_heff
#         hetero.write_heff() 
     # reduced matrix
     if hetero.block_diagonal: 
       from copy import deepcopy
       heff = deepcopy(intra)
       heff[0][0] = intra[0][0] + selfl
       heff[-1][-1] = intra[-1][-1] + selfr
       dd = (energy+1j*gf.eps)*iden
       for i in range(len(intra)):  # add the diagonal energy part
         heff[i][i] = heff[i][i] - dd  # this has the wrong sign!!
      # now change the sign
       for i in range(len(intra)):  
         for j in range(len(intra)):  
           try:
             heff[i][j] = -heff[i][j]
           except:
             heff[i][j] = heff[i][j]
      # save the green function
       hetero.heff = heff
#       if save_heff:
#         print "Saving effective hamiltonian in ",hetero.file_heff
#         from scipy.sparse import bmat
#         hetero.heff = bmat(heff)
#         hetero.write_heff() 
      # calculate the inverse
       from green import gauss_inverse  # routine to invert the matrix
       # list with the diagonal matrices
       ldos_e = ldos*0.0 # initialice ldos at this energy
       ii = 0 # counter for the element
       for i in range(len(heff)): # loop over blocks
         gci = gauss_inverse(heff,i,i) # calculate each block element      
         for j in range(len(heff[0][0])): # loop over each block
           ldos_e[ii] = -gci[j,j].imag 
           ii += 1 # increase counter
       if not ii==dimhc:
         print "Wrong dimensions",ii,dimhc
         raise
       ldos += ldos_e # add to the total ldos
# save the effective hamiltonian

   if hetero.has_spin: # resum ldos if there is spin degree of freedom
     ldos = [ldos[2*i]+ldos[2*i+1] for i in range(len(ldos)/2)]
   if hetero.has_eh: # resum ldos if there is eh 
     ldos = [ldos[2*i]+ldos[2*i+1] for i in range(len(ldos)/2)]
  #   ne = len(ldos)/2
  #   ldos = [ldos[i]+ldos[ne+i] for i in range(ne)]


   ldos = np.array(ldos) # transform into an array

   if min(ldos)<0.0:
     print "Negative density of states"
     print ldos
     raise

   g = hetero.central_geometry # geometry of the central part
   fldos = open("LDOS.OUT","w") # open file for ldos
   fldos.write("# X   Y    LDOS\n")
   for (ix,iy,il) in zip(g.x,g.y,ldos):
     fldos.write(str(ix)+"  "+str(iy)+"  "+str(il)+"\n")
   fldos.close()
   # scale the ldos
   # save the LDOS in a file
   if True:
#   if True>0.001:
#     ldos = np.sqrt(ldos)
#     ldos = np.arctan(7.*ldos/max(ldos))
     print "Sum of the LDOS =",sum(ldos)
     ldos = ldos*300/max(ldos)
   else:
     ldos = ldos*0.0





   # now create the figure
   fig = py.figure() # create figure
   fig.subplots_adjust(0.2,0.2)
   fig.set_facecolor("white") # face in white
   sldos = fig.add_subplot(111) # create subplot for the DOS

   # plot the lattice
   if not len(g.x)==len(ldos):
     raise 
   sldos.scatter(g.x,g.y,color="red",s=ldos) # plot the lattice
   sldos.scatter(g.x,g.y,color="black",s=4) # plot the lattice
   sldos.set_xlabel("X")
   sldos.set_xlabel("Y")
   sldos.axis("equal") # same scale in axes



   return fig







def create_leads_and_central_list(h_right,h_left,list_h_central):
  """ Creates an heterojunction by giving the hamiltonian
     of the leads and the list of the center """
  # check the hamiltonians
#  h_right.check()
#  h_left.check()
  ht = heterostructure(h_right) # create heterostructure
  # assign matrices of the leads
  ht.right_intra = h_right.intra 
  ht.right_inter = h_right.inter 
  ht.left_intra = h_left.intra 
  ht.left_inter = h_left.inter.H
  # create matrix of the central part and couplings to the leads
  from scipy.sparse import csc_matrix,bmat
  z = csc_matrix(h_right.intra*0.0j) # zero matrix
  num_central = len(list_h_central) # length of the central hamiltonian
  # create list of lists for the central part
  hc = [[None for i in range(num_central)] for j in range(num_central)]
  # create elements of the central hamiltonian
  for i in range(num_central):  # intra term of the central blocks
    hc[i][i] = list_h_central[i].intra # intra term of the iesim
  for i in range(num_central-1):  # intra term of the central blocks
    tr = list_h_central[i].inter + list_h_central[i+1].inter
    tr = tr/2.   # mean value of the hoppings
    hc[i][i+1] = tr # inter term of the iesim
    hc[i+1][i] = tr.H # inter term of the iesim

  # hoppings to the leads
  tcr = h_right.inter   # this is a matrix
  tcl = h_left.inter.H  # this is a matrix
  ht.block_diagonal = True
  # hoppings from the center to the leads
  ht.right_coupling = h_right.inter
  ht.left_coupling = h_left.inter.H
  # assign central hamiltonian
  ht.central_intra = hc
  # and modify the geometry of the central part
  ht.central_geometry.supercell(num_central) 
  # put if it si sparse
  ht.is_sparse = list_h_central[0].is_sparse
  return ht


def eigenvalues(hetero,numeig=10,effective=False,gf=None,full=False):
  """ Gets the lowest eigenvalues of the central part of the hamiltonian"""
  if not hetero.block_diagonal:
    print """ heterounction in eigenvalues must be block diagonal"""
    raise
  # if effective hamiltonian, just calculate the eigenvalues
  if effective: # effective hamiltonian
    print "Calculating eigenvalues of effective hamiltonian..."
    if hetero.heff == None:  #if hasn't been calculated so far
      effective_central_hamiltonian(hetero,write=False)
    heff = hetero.heff # store the list of list with central EFF ham
    import scipy.sparse.linalg as lg
    eig,eigvec = lg.eigs(heff,k=numeig,which="LM",sigma=0.0)
    return eig

  from scipy.sparse import csc_matrix,bmat
  import scipy.sparse.linalg as lg
  if not effective: # do not use the effective hamiltonian, only the central
    intra = hetero.central_intra # store the list of list with central ham
  numb = len(intra) # number of central blocks
  intrasp = [[None for i in range(numb)] for j in range(numb)]
  # assign onsite and couplings in sparse form
  if not hetero.is_sparse:
    for i in range(numb):
      intrasp[i][i] = csc_matrix(intra[i][i])
    for i in range(numb-1):
      intrasp[i][i+1] = csc_matrix(intra[i][i+1])
      intrasp[i+1][i] = csc_matrix(intra[i+1][i])
  if hetero.is_sparse:
    for i in range(numb):
      intrasp[i][i] = intra[i][i]
    for i in range(numb-1):
      intrasp[i][i+1] = intra[i][i+1]
      intrasp[i+1][i] = intra[i+1][i]
  intrasp = bmat(intrasp) # create sparse matrix
  if effective:
    evals,evecs = lg.eigs(intrasp,k=numeig,which="LM",sigma=0.0)
  if not effective:
    if full:  # full diagonalization
      from scipy.linalg import eigvalsh
      evals = eigvalsh(intrasp.todense())
    else:
      evals,evecs = lg.eigsh(intrasp,k=numeig,which="LM",sigma=0.0)
  return evals




def effective_central_hamiltonian(hetero,energy=0.0,delta=0.001,write=False):
   """ Plots the local density of states in the central part"""
   from green import green_renormalization
   from green import dyson
   from hamiltonians import is_number
   if not hetero.block_diagonal:
     intra = hetero.central_intra # central intraterm   
   if hetero.block_diagonal:
     intra = hetero.central_intra[0][0] # when it is diagonal
   # perform dyson calculation
   intra = hetero.right_intra
   inter = hetero.right_inter
#   gr = dyson(intra,inter,is_sparse=hetero.is_sparse)
   ggg,gr = green_renormalization(intra,inter,energy=energy,delta=delta)
   hetero.right_green = gr # save green function
   # left green function
   intra = hetero.left_intra
   inter = hetero.left_inter
#   gl = dyson(intra,inter,is_sparse=hetero.is_sparse)
   ggg,gl = green_renormalization(intra,inter,energy=energy,delta=delta)
   hetero.left_green = gl # save green function
   # save green functions
#   hetero.write_green()
#   print "Saved green functions"
   # left selfenergy
   inter = hetero.left_coupling
   selfl = inter*gl*inter.H # left selfenergy
   # right selfenergy
   inter = hetero.right_coupling
   selfr = inter*gr*inter.H # right selfenergy
   # central green function
   intra = hetero.central_intra
# dyson equation for the center     
   # full matrix
   if not hetero.block_diagonal:
     heff = intra + selfl + selfr
     hetero.heff = heff
     if save_heff:
       print "Saving effective hamiltonian in ",hetero.file_heff
       if write: # save hamiltonian if desired
         hetero.write_heff() 
   # reduced matrix
   if hetero.block_diagonal: 
     from copy import deepcopy
     heff = deepcopy(intra)
     heff[0][0] = intra[0][0] + selfl
     heff[-1][-1] = intra[-1][-1] + selfr
    # save the green function
     from scipy.sparse import bmat
     hetero.heff = bmat(heff)
     if write: # save hamiltonian
       print "Saving effective hamiltonian in ",hetero.file_heff
       hetero.write_heff() 

#
#def save_wave(hetero,energy=0.0,output_file="wave.out"):
#  """Saves in a file the closest eigenfunction to a certain energy"""
#  if not hetero.block_diagonal:
#    print """ heterounction in save_swave must be block diagonal"""
#    raise
#  # if effective hamiltonian, just calculate the eigenvalues
#  if effective: # do not use the effective hamiltonian, only the central
#    print "Calculating eigenvalues of effective hamiltonian..."
#    if hetero.heff == None:  #if hasn't been calculated so far
#      effective_central_hamiltonian(hetero,write=False,gf=gf)
#    heff = hetero.heff # store the list of list with central EFF ham
#    import scipy.sparse.linalg as lg
#    eig,eigvec = lg.eigs(heff,k=numeig,which="LM",sigma=0.0)
#
#  # continue if the central hamiltonian is right
#  from scipy.sparse import csc_matrix,bmat
#  import scipy.sparse.linalg as lg
#  if not effective: # do not use the effective hamiltonian, only the central
#    intra = hetero.central_intra # store the list of list with central ham
#  numb = len(intra) # number of central blocks
#  intrasp = [[None for i in range(numb)] for j in range(numb)]
#  # assign onsite and couplings in sparse form
#  if not hetero.is_sparse:
#    for i in range(numb):
#      intrasp[i][i] = csc_matrix(intra[i][i])
#    for i in range(numb-1):
#      intrasp[i][i+1] = csc_matrix(intra[i][i+1])
#      intrasp[i+1][i] = csc_matrix(intra[i+1][i])
#  if hetero.is_sparse:
#    for i in range(numb):
#      intrasp[i][i] = intra[i][i]
#    for i in range(numb-1):
#      intrasp[i][i+1] = intra[i][i+1]
#      intrasp[i+1][i] = intra[i+1][i]
#  intrasp = bmat(intrasp) # create sparse matrix
#  eigenval,eigenvec = lg.eigs(intrasp,k=1,which="LM",sigma=energy)
#  fo = open(output_file,"w"):
#    g = hetero.central_geometry
#    if not 4*len(g.x)==len(eigenvec):
#      print "Wrong dimensions in save_wave"
#      raise
#    print eigenvec
#    fo.write()
#  fo.close()
#  return evals
#
