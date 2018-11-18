from __future__ import print_function
import green
import numpy as np
import scipy.linalg as lg
from bandstructure import smalleig # arpack diagonalization
import time
import timing
import kpm
import checkclass

try:
#  raise
  import dosf90 
  use_fortran = True
except:
  print("Something wrong with FORTRAN in DOS")
  use_fortran = False


def calculate_dos(es,xs,d,use_fortran=use_fortran,w=None):
  if w is None: w = np.zeros(len(es)) + 1.0 # initialize
  if use_fortran: # use fortran routine
    import dosf90 
    return dosf90.calculate_dos(es,xs,d,w) # use the Fortran routine
  else:
#    if w is None: # no weights
#      ndos = len(xs)
#      delta = d/10.
#      ys,xs = np.histogram(es,bins=ndos) # create the histogram
#      lorentz = np.linspace(-1.,1.,len(ys)) # number of energies
#      lorentz = delta/(delta*delta + lorentz*lorentz) # smoothing function
#      ys = np.convolve(lorentz,ys,mode="same") # convolve lorentz and histogram
#      return ys
#    else: # conventional mode
      ys = np.zeros(xs.shape[0]) # initialize
      for (e,iw) in zip(es,w): # loop over energies
         de = xs - e # E - Ei
         de = d/(d*d + de*de) # 1/(delta^2 + (E-Ei)^2)
         ys += de*iw # times the weight
      return ys




def dos_surface(h,output_file="DOS.OUT",
                 energies=np.linspace(-1.,1.,20),delta=0.001):
  """Calculates the DOS of a surface, and writes in file"""
  if h.dimensionality!=1: raise # only for 1d
  fo = open(output_file,"w")
  fo.write("# energy, DOS surface, DOS bulk\n")
  for e in energies: # loop over energies
    print("Done",e)
    gb,gs = green.green_renormalization(h.intra,h.inter,energy=e,delta=delta)
    gb = -gb.trace()[0,0].imag
    gs = -gs.trace()[0,0].imag
    fo.write(str(e)+"     "+str(gs)+"    "+str(gb)+"\n")
  fo.close()




def dos0d(h,es=None,delta=0.01,i=None):
  """Calculate density of states of a 0d system"""
  if es is None: es = np.linspace(-4,4,500)
  ds = [] # empty list
  if h.dimensionality==0:  # only for 0d
    iden = np.identity(h.intra.shape[0],dtype=np.complex) # create identity
    for e in es: # loop over energies
      g = ( (e+1j*delta)*iden -h.intra ).I # calculate green function
      if i is None: d = -g.trace()[0,0].imag
      elif checkclass.is_iterable(i): # iterable list 
          d = sum([-g[ii,ii].imag for ii in i])
      else: d = -g[i,i].imag # assume an integer
      ds.append(d)  # add this dos
  else: raise # not implemented...
  write_dos(es,ds)
  return ds





def dos0d_kpm(h,use_kpm=True,scale=10,npol=100,ntries=100,fun=None):
  """ Calculate density of states of a 1d system"""
  if h.dimensionality!=0: raise # only for 0d
  if not use_kpm: raise # only using KPM
  h.turn_sparse() # turn the hamiltonian sparse
  mus = np.array([0.0j for i in range(2*npol)]) # initialize polynomials
  import kpm
  mus = kpm.random_trace(h.intra/scale,ntries=ntries,n=npol,fun=fun)
  xs = np.linspace(-0.9,0.9,4*npol) # x points
  ys = kpm.generate_profile(mus,xs) # generate the profile
  write_dos(xs*scale,ys) # write in file


def dos0d_sites(h,sites=[0],scale=10.,npol=500,ewindow=None,refine_e=1.0):
  """ Calculate density of states of a 1d system for a certain orbitals"""
  if h.dimensionality!=0: raise # only for 1d
  h.turn_sparse() # turn the hamiltonian sparse
  mus = np.array([0.0j for i in range(2*npol)]) # initialize polynomials
  import kpm
  hk = h.intra # hamiltonian
  for isite in sites:
    mus += kpm.local_dos(hk/scale,i=isite,n=npol)
  if ewindow is None:  xs = np.linspace(-0.9,0.9,int(npol*refine_e)) # x points
  else:  xs = np.linspace(-ewindow/scale,ewindow/scale,npol) # x points
  ys = kpm.generate_profile(mus,xs) # generate the profile
  write_dos(xs*scale,ys) # write in file











def write_dos(es,ds,output_file="DOS.OUT"):
  """ Write DOS in a file"""
  f = open(output_file,"w")
  for (e,d) in zip(es,ds):
    f.write(str(e)+"     ")
    f.write(str(d.real)+"\n")
  f.close()




def dos1d(h,use_kpm=False,scale=10.,nk=100,npol=100,ntries=2,
          ndos=1000,delta=0.01,ewindow=None,frand=None):
  """ Calculate density of states of a 1d system"""
  if h.dimensionality!=1: raise # only for 1d
  ks = np.linspace(0.,1.,nk,endpoint=False) # number of kpoints
  if not use_kpm: # conventional method
    hkgen = h.get_hk_gen() # get generator
#    delta = 16./(nk*h.intra.shape[0]) # smoothing
    calculate_dos_hkgen(hkgen,ks,ndos=ndos,delta=delta) # conventiona algorithm
  else:
    h.turn_sparse() # turn the hamiltonian sparse
    hkgen = h.get_hk_gen() # get generator
    yt = np.zeros(ndos) # number of dos
    import kpm
    ts = timing.Testimator("DOS") 
    for i in range(nk): # loop over kpoints
      k = np.random.random(3) # random k-point
      hk = hkgen(k) # hamiltonian
      (xs,yi) = kpm.tdos(hk,scale=scale,npol=npol,frand=frand,ewindow=ewindow,
                  ntries=ntries,ne=ndos)
      yt += yi # Add contribution
      ts.remaining(i,nk)
    yt /= nk # normalize
    write_dos(xs,yt) # write in file
    print()
    return xs,yt




def dos1d_sites(h,sites=[0],scale=10.,nk=100,npol=100,info=False,ewindow=None):
  """ Calculate density of states of a 1d system for a certain orbitals"""
  if h.dimensionality!=1: raise # only for 1d
  ks = np.linspace(0.,1.,nk,endpoint=False) # number of kpoints
  h.turn_sparse() # turn the hamiltonian sparse
  hkgen = h.get_hk_gen() # get generator
  mus = np.array([0.0j for i in range(2*npol)]) # initialize polynomials
  import kpm
  for k in ks: # loop over kpoints
    hk = hkgen(k) # hamiltonian
    for isite in sites:
      mus += kpm.local_dos(hk/scale,i=isite,n=npol)
    if info: print("Done",k)
  mus /= nk # normalize by the number of kpoints
  if ewindow is None:  xs = np.linspace(-0.9,0.9,npol) # x points
  else:  xs = np.linspace(-ewindow/scale,ewindow/scale,npol) # x points
  ys = kpm.generate_profile(mus,xs) # generate the profile
  write_dos(xs*scale,ys) # write in file


def calculate_dos_hkgen(hkgen,ks,ndos=100,delta=None,
         is_sparse=False,numw=10,window=None):
  """Calculate density of states using the ks given on input"""
  es = np.zeros((len(ks),hkgen(ks[0]).shape[0])) # empty list
  tr = timing.Testimator("DOS",maxite=len(ks))
  import parallel
  def fun(k): # function to execute
    if parallel.cores==1: tr.iterate() # print the info
    hk = hkgen(k) # Hamiltonian
    t0 = time.clock() # time
    if is_sparse: # sparse Hamiltonian 
      return smalleig(hk,numw=numw).tolist() # eigenvalues
    else: # dense Hamiltonian
      return lg.eigvalsh(hk).tolist() # get eigenvalues
#  for ik in range(len(ks)):  
  out = parallel.pcall(fun,ks) # launch all the processes
  es = [] # empty list
  for o in out: es += o # concatenate
#    tr.remaining(ik,len(ks))
#  es = es.reshape(len(es)*len(es[0])) # 1d array
  es = np.array(es) # convert to array
  nk = len(ks) # number of kpoints
  if delta is None: delta = 5./nk # automatic delta
  if window is None:
    xs = np.linspace(np.min(es)-.5,np.max(es)+.5,ndos) # create x
  else:
    xs = np.linspace(-window,window,ndos) # create x
  ys = calculate_dos(es,xs,delta) # use the Fortran routine
  ys /= nk # normalize 
  write_dos(xs,ys) # write in file
  print("\nDOS finished")











def dos2d(h,use_kpm=False,scale=10.,nk=100,ntries=1,delta=None,
          ndos=500,numw=20,random=True,kpm_window=1.0,window=None):
  """ Calculate density of states of a 2d system"""
  if h.dimensionality!=2: raise # only for 2d
  ks = []
  from klist import kmesh
  ks = kmesh(h.dimensionality,nk=nk)
  if random:
    ks = [np.random.random(2) for ik in ks]
    print("Random k-mesh")
  if not use_kpm: # conventional method
    hkgen = h.get_hk_gen() # get generator
    if delta is None: delta = 6./nk
# conventiona algorithm
    calculate_dos_hkgen(hkgen,ks,ndos=ndos,delta=delta,
                          is_sparse=h.is_sparse,numw=numw,window=window) 
  else: # use the kpm
    npol = ndos//10
    h.turn_sparse() # turn the hamiltonian sparse
    hkgen = h.get_hk_gen() # get generator
    mus = np.array([0.0j for i in range(2*npol)]) # initialize polynomials
    import kpm
    tr = timing.Testimator("DOS")
    ik = 0
    for k in ks: # loop over kpoints
#      print("KPM DOS at k-point",k)
      ik += 1
      tr.remaining(ik,len(ks))      
      if random: 
        kr = np.random.random(2)
        print("Random sampling in DOS")
        hk = hkgen(kr) # hamiltonian
      else: hk = hkgen(k) # hamiltonian
      mus += kpm.random_trace(hk/scale,ntries=ntries,n=npol)
    mus /= len(ks) # normalize by the number of kpoints
    xs = np.linspace(-0.9,0.9,ndos)*kpm_window # x points
    ys = kpm.generate_profile(mus,xs) # generate the profile
    write_dos(xs*scale,ys) # write in file
    return (xs,ys)




def dos3d(h,scale=10.,nk=20,delta=None,ndos=100,random=False):
  """ Calculate density of states of a 2d system"""
  if h.dimensionality!=3: raise # only for 2d
  ks = [np.random.random(3) for i in range(nk)] # number of kpoints
  hkgen = h.get_hk_gen() # get generator
  if delta is None: delta = 10./ndos # smoothing
  calculate_dos_hkgen(hkgen,ks,ndos=ndos,delta=delta) # conventional algorithm











def dos2d_ewindow(h,energies=np.linspace(-1.,1.,30),delta=None,info=False,
                    use_green=True,nk=300,mode="adaptive"):
  """Calculate the density of states in certain eenrgy window"""
  ys = [] # density of states
  if delta is None: # pick a good delta value
    delta = 0.1*(max(energies) - min(energies))/len(energies)
  if use_green:
    from green import bloch_selfenergy
    for energy in energies:
      (g,selfe) = bloch_selfenergy(h,nk=nk,energy=energy, delta=delta,
                   mode=mode)
      ys.append(-g.trace()[0,0].imag)
      if info: print("Done",energy)
    write_dos(energies,ys) # write in file
    return
  else: # do not use green function    
    import scipy.linalg as lg
    kxs = np.linspace(0.,1.,nk)
    kys = np.linspace(0.,1.,nk)
    hkgen= h.get_hk_gen() # get hamiltonian generator
    ys = energies*0.
    weight = 1./(nk*nk)
    for ix in kxs:
      for iy in kys:
        k = np.array([ix,iy,0.]) # create kpoint
        hk = hkgen(k) # get hk hamiltonian
        evals = lg.eigvalsh(hk) # get eigenvalues
        ys += weight*calculate_dos(evals,energies,delta) # add this contribution
      if info: print("Done",ix)
    write_dos(energies,ys) # write in file
    return






def dos1d_ewindow(h,energies=np.linspace(-1.,1.,30),delta=None,info=False,
                    use_green=True,nk=300):
  """Calculate the density of states in certain energy window"""
  ys = [] # density of states
  if delta is None: # pick a good delta value
    delta = 0.1*(max(energies) - min(energies))/len(energies)
  if True: # do not use green function    
    import scipy.linalg as lg
    kxs = np.linspace(0.,1.,nk)
    hkgen= h.get_hk_gen() # get hamiltonian generator
    ys = energies*0.
    weight = 1./(nk)
    for ix in kxs:
      hk = hkgen(ix) # get hk hamiltonian
      evals = lg.eigvalsh(hk) # get eigenvalues
      ys += weight*calculate_dos(evals,energies,delta) # add this contribution
    if info: print("Done",ix)
    write_dos(energies,ys) # write in file
    return











def dos_ewindow(h,energies=np.linspace(-1.,1.,30),delta=None,info=False,
                    use_green=True,nk=300):
  """ Calculate density of states in an energy window"""
  if h.dimensionality==2: # two dimensional
    dos2d_ewindow(h,energies=energies,delta=delta,info=info,
                    use_green=use_green,nk=nk)
  elif h.dimensionality==1: # one dimensional
    dos1d_ewindow(h,energies=energies,delta=delta,info=info,
                    use_green=use_green,nk=nk)
  else: raise # not implemented






def convolve(x,y,delta=None):
  """Add a broadening to a DOS"""
  if delta is None: return y # do nothing
  delta = np.abs(delta) # absolute value
  xnew = np.linspace(-1.,1.,len(y)) # array
  d2 = delta/(np.max(x) - np.min(x)) # effective broadening
  fconv = d2/(xnew**2 + d2**2) # convolving function
  yout = np.convolve(y,fconv,mode="same") # same size
  # ensure the normaliation is the same
  ratio = np.sum(np.abs(y))/np.sum(np.abs(yout))
  yout *= ratio
  return yout # return new array





def dos_kpm(h,scale=10.0,ewindow=4.0,ne=1000,delta=0.01,ntries=10,nk=100):
  """Calculate the KDOS bands using the KPM"""
  hkgen = h.get_hk_gen() # get generator
  numk = nk**h.dimensionality
  tr = timing.Testimator("DOS",maxite=numk) # generate object
  ytot = np.zeros(ne) # initialize
  for ik in range(numk): # loop over kpoints
    tr.iterate()
    hk = hkgen(np.random.random(3)) # get Hamiltonian
    npol = int(scale/delta) # number of polynomials
    (x,y) = kpm.tdos(hk,scale=scale,npol=npol,ne=ne,
                   ewindow=ewindow,ntries=ntries) # compute
    ytot += y # add contribution
  ytot /= nk # normalize
  np.savetxt("DOS.OUT",np.matrix([x,ytot]).T) # save in file



def dos(h,energies=np.linspace(-4.0,4.0,400),delta=0.01,nk=10,
            use_kpm=False,scale=10.,ntries=10,mode="ED"):
  """Calculate the density of states"""
  if use_kpm: # KPM
    ewindow = max([abs(min(energies)),abs(min(energies))]) # window
    dos_kpm(h,scale=scale,ewindow=ewindow,ne=len(energies),delta=delta,
                   ntries=ntries,nk=nk)
  else: # conventional methods
    if mode=="ED": # exact diagonalization
      if h.dimensionality==0:
        return dos0d(h,es=energies,delta=delta)
      elif h.dimensionality==1:
        return dos1d(h,ndos=len(energies),delta=delta,nk=nk)
      elif h.dimensionality==2:
        return dos2d(h,use_kpm=False,nk=100,ntries=1,delta=delta,
            ndos=len(energies),random=True,window=np.max(np.abs(energies)))
      else: raise
    elif mode=="Green": # Green function formalism
      if h.dimensionality==0:
        return dos0d(h,es=energies,delta=delta) # same as before
      elif h.dimensionality>0: # Bigger dimensionality
        from green import bloch_selfenergy
        tr = timing.Testimator("KDOS") # generate object
        ie = 0
        out = [] # storage
        for e in energies: # loop
          tr.remaining(ie,len(energies)) # print status
          ie += 1 # increase
          g = bloch_selfenergy(h,energy=e,delta=delta,mode="adaptive")[0]
          out.append(-g.trace()[0,0].imag) # store dos
        np.savetxt("DOS.OUT",np.matrix([energies,out]).T) # write in a file
        return energies,np.array(out) # return


def bulkandsurface(h1,energies=np.linspace(-1.,1.,100),operator=None,
                    delta=0.01,hs=None,nk=30):
  """Compute the DOS of the bulk and the surface"""
  tr = timing.Testimator("KDOS") # generate object
  ik = 0
  h1 = h1.get_multicell() # multicell Hamiltonian
  kpath = [np.random.random(3) for i in range(nk)] # loop
  dosout = np.zeros((2,len(energies))) # output DOS
  for k in kpath:
    tr.remaining(ik,len(kpath)) # generate object
    ik += 1
    outs = green.surface_multienergy(h1,k=k,energies=energies,delta=delta,hs=hs)
    ie = 0
    for (energy,out) in zip(energies,outs): # loop over energies
      # compute dos
      ig = 0
      for g in out: # loop
        if operator is None: d = -g.trace()[0,0].imag # only the trace 
        elif callable(operator): d = operator(g,k=k) # call the operator
        else:  d = -(g*operator).trace()[0,0].imag # assume it is a matrix
        dosout[ig,ie] += d # store
        ig += 1 # increase site
      ie += 1 # increase energy
  # write in file
  dosout/=nk # normalize
  np.savetxt("DOS_BULK_SURFACE.OUT",np.matrix([energies,dosout[0],dosout[1]]).T)






