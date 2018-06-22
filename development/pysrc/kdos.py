from __future__ import print_function,division
import numpy as np
import green
import dos
import klist
import kpm
import timing

def write_kdos(k=0.,es=[],ds=[],new=True):
  """ Write KDOS in a file"""
  if new: f = open("KDOS.OUT","w") # open new file
  else: f = open("KDOS.OUT","a") # continue writting
  for (e,d) in zip(es,ds): # loop over e and dos
    f.write(str(k)+"     ")
    f.write(str(e)+"     ")
    f.write(str(d)+"\n")
  f.close()





def kdos1d_sites(h,sites=[0],scale=10.,nk=100,npol=100,kshift=0.,
                  ewindow=None,info=False):
  """ Calculate kresolved density of states of
  a 1d system for a certain orbitals"""
  if h.dimensionality!=1: raise # only for 1d
  ks = np.linspace(0.,1.,nk) # number of kpoints
  h.turn_sparse() # turn the hamiltonian sparse
  hkgen = h.get_hk_gen() # get generator
  if ewindow is None:  xs = np.linspace(-0.9,0.9,nk) # x points
  else:  xs = np.linspace(-ewindow/scale,ewindow/scale,nk) # x points
  import kpm
  write_kdos() # initialize file
  for k in ks: # loop over kpoints
    mus = np.array([0.0j for i in range(2*npol)]) # initialize polynomials
    hk = hkgen(k+kshift) # hamiltonian
    for isite in sites:
      mus += kpm.local_dos(hk/scale,i=isite,n=npol)
    ys = kpm.generate_profile(mus,xs) # generate the profile
    write_kdos(k,xs*scale,ys,new=False) # write in file (append)
    if info: print("Done",k)


def surface(h,energies=None,klist=None,delta=0.01):
  bout = [] # empty list, bulk
  sout = [] # empty list, surface
  for k in klist:
    for energy in energies:
      gs,sf = green.green_kchain(h,k=k,energy=energy,delta=delta,only_bulk=False) 
      bout.append(gs.trace()[0,0].imag) # bulk
      sout.append(sf.trace()[0,0].imag) # surface
  bout = np.array(bout).reshape((len(energies),len(klist))) # convert to array
  sout = np.array(sout).reshape((len(energies),len(klist))) # convert to array
  return (bout.transpose(),sout.transpose())


def write_surface(h,energies=np.linspace(-.5,.5,100),klist=np.linspace(0.,1.,100),delta=None):
  if delta is None: delta = (np.max(energies)-np.min(energies))/len(energies)
  if h.dimensionality==2:
    write_surface_2d(h,energies=energies,klist=klist,delta=delta)
  elif h.dimensionality==3:
    write_surface_3d(h,energies=energies,klist=klist,delta=delta)
  else: raise


import multicell

def write_surface_2d(h,energies=None,klist=None,delta=0.01):
  if h.dimensionality != 2: raise # only for 2d
  if h.is_multicell: # multicell Hamiltonian
    h = multicell.turn_no_multicell(h) # convert into a normal Hamiltonian
  bout = [] # empty list, bulk
  sout = [] # empty list, surface
  if klist is None: klist = np.linspace(-.5,.5,50)
  if energies is None: energies = np.linspace(-.5,.5,50)
  fo  = open("KDOS.OUT","w") # open file
  for k in klist:
    print("Doing k-point",k)
    for energy in energies:
      gs,sf = green.green_kchain(h,k=k,energy=energy,delta=delta,only_bulk=False) 
      db = -gs.trace()[0,0].imag # bulk
      ds = -sf.trace()[0,0].imag # surface
      fo.write(str(k)+"   "+str(energy)+"   "+str(ds)+"   "+str(db)+"\n")
  fo.close()


def write_surface_3d(h,energies=None,klist=None,delta=0.01):
  raise # not implemented
  if h.dimensionality != 3: raise # only for 3d
  ho = h.copy() # copy Hamiltonian
  ho = ho.turn_multicell() # multicell Hamiltonian
  bout = [] # empty list, bulk
  sout = [] # empty list, surface
  if klist is None: raise
  if energies is None: energies = np.linspace(-.5,.5,50)
  fo  = open("KDOS.OUT","w") # open file
  for k in klist:
    for energy in energies:
      gs,sf = green.green_kchain(h,k=k,energy=energy,delta=delta,only_bulk=False) 
      db = -gs.trace()[0,0].imag # bulk
      ds = -sf.trace()[0,0].imag # surface
      fo.write(str(k)+"   "+str(energy)+"   "+str(ds)+"   "+str(db)+"\n")



def kdos_bands(h,use_kpm=True,kpath=None,scale=10.0,
                 ewindow=4.0,ne=1000,delta=0.01,ntries=10):
  """Calculate the KDOS bands using the KPM"""
  if not use_kpm: raise # nope
  hkgen = h.get_hk_gen() # get generator
  if kpath is None: kpath = klist.default(h.geometry) # default
  fo = open("KDOS_BANDS.OUT","w") # open file
  ik = 0
  tr = timing.Testimator("KDOS") # generate object
  for k in kpath: # loop over kpoints
    tr.remaining(ik,len(kpath))
    hk = hkgen(k) # get Hamiltonian
    npol = int(scale/delta) # number of polynomials
    (x,y) = kpm.tdos(hk,scale=scale,npol=npol,ne=ne,
                   ewindow=ewindow,ntries=ntries) # compute
    for (ix,iy) in zip(x,y): # loop
      fo.write(str(ik/len(kpath))+"   ")
      fo.write(str(ix)+"   ")
      fo.write(str(iy)+"\n")
    fo.flush()
    ik += 1
  fo.close()









def write_surface_kpm(h,ne=400,klist=None,scale=4.,npol=200,w=20,ntries=20):
  """Write the surface DOS using the KPM"""
  if klist is None: klist = np.linspace(-.5,.5,50)
  import kpm
  fo  = open("KDOS.OUT","w") # open file
  for k in klist:
    print("Doing kpoint",k)
    if h.dimensionality==2: 
      (intra,inter) = h.kchain(k) # k hamiltonian
      (es,ds,dsb) = kpm.edge_dos(intra,inter,scale=scale,w=w,npol=npol,
                            ne=ne,bulk=True)
    # if the Hamiltonian is 1d from the beginning
    elif h.dimensionality==1: 
      intra,inter = h.intra,h.inter # 1d hamiltonian
      dd = h.intra.shape[0] # dimension
      inde = np.zeros(dd) # array with zeros
      indb = np.zeros(dd) # array with zeros
      for i in range(dd//10): # one tenth
        inde[i] = 1. # use this one
        indb[4*dd//10 + i] = 1. # use this one
      def gedge(): return (np.random.random(len(inde))-0.5)*inde
      def gbulk(): return (np.random.random(len(indb))-0.5)*(indb)
      # hamiltonian
      h0 = intra + inter*np.exp(1j*np.pi*2.*k) + (inter*np.exp(1j*np.pi*2.*k)).H
      xs = np.linspace(-0.9,0.9,4*npol) # x points
      es = xs*scale
      # calculate the bulk
      mus = kpm.random_trace(h0/scale,ntries=ntries,n=npol,fun=gbulk)
      dsb = kpm.generate_profile(mus,xs) # generate the profile
      # calculate the edge
      mus = kpm.random_trace(h0/scale,ntries=ntries,n=npol,fun=gedge)
      ds = kpm.generate_profile(mus,xs) # generate the profile
    else: raise
    for (e,d1,d2) in zip(es,ds,dsb):
      fo.write(str(k)+"   "+str(e)+"   "+str(d1)+"    "+str(d2)+"\n")
  fo.close()




