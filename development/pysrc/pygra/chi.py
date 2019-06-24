import numpy as np
import scipy.linalg as lg
from . import parallel
try:
    from . import chif90
except:
    print("Error, chif90 not found")


def chargechi(h,i=0,j=0,es=np.linspace(-3.0,3.0,100),delta=0.01,temp=1e-7):
    """Compute charge response function"""
    if h.dimensionality!=0: raise
    hk = h.get_hk_gen() # get generator
    m = hk(0) # get Hamiltonian
    esh,ws = lg.eigh(m)
    ws = np.transpose(ws)
    if i<0: raise
    if j<0: raise
    return es,chif90.elementchi(ws,esh,ws,esh,es,i+1,j+1,temp,delta)



def chargechi_nowf(h,i=0,j=0,es=np.linspace(-3.0,3.0,100),delta=0.01,temp=1e-7):
    """Compute charge response function"""
    if h.dimensionality!=0: raise
    hk = h.get_hk_gen() # get generator
    m = hk(0) # get Hamiltonian
    esh,ws = lg.eigh(m)
    ws = np.transpose(ws)*0. + 1.
    if i<0: raise
    if j<0: raise
    return es,chif90.elementchi(ws,esh,ws,esh,es,i+1,j+1,temp,delta)







def chargechi_row(h,i=0,es=np.linspace(-3.0,3.0,100),delta=1e-6,temp=1e-7):
    """Compute charge response function"""
    if h.dimensionality!=0: raise
    hk = h.get_hk_gen() # get generator
    m = hk(0) # get Hamiltonian
    esh,ws = lg.eigh(m)
    ws = np.transpose(ws)
    out = [] # store
    if i<0: raise
    f = lambda j: chif90.elementchi(ws,esh,ws,esh,es,i+1,j+1,temp,delta)
    out = parallel.pcall(f,range(m.shape[0])) # parallel call
    return np.array(out)



def chargechi_nowf(h,i=0,j=0,es=np.linspace(-3.0,3.0,100),delta=0.01,temp=1e-7):
    """Compute charge response function"""
    if h.dimensionality!=0: raise
    hk = h.get_hk_gen() # get generator
    m = hk(0) # get Hamiltonian
    esh,ws = lg.eigh(m)
    ws = np.transpose(ws)*0. + 1.
    if i<0: raise
    if j<0: raise
    return es,chif90.elementchi(ws,esh,ws,esh,es,i+1,j+1,temp,delta)



def chargechi_reciprocal(h,i=None,
        es=np.linspace(-4.,4.,200),delta=1e-3,**kwargs):
    """Return the charge susceptibility in reciprocal space"""
    if h.dimensionality!=0: raise
    if i is None: i = h.geometry.get_central(1)[0]
    # compute correlators
    cs = chargechi_row(h,es=es,delta=delta,**kwargs) 
    # get the distances
    dx = np.abs(h.geometry.r[:,0] - h.geometry.r[i,0])
    # now do the interpolation
    from scipy.interpolate import interp1d
    xs = np.linspace(0.,np.max(dx),len(dx))
    cs = np.array([interp1d(xs,ic)(xs) for ic in cs.T]).T
#    if len(cs)!=len(h.geometry.r[:,0]): raise
    freq = np.fft.fftfreq(len(cs[:,0])) # get the frequencies
    # perform the fourier transform
    fcs = np.array([np.fft.fft(ic) for ic in cs.T])
    fo = open("CHIK.OUT","w")
    print(fcs.shape,es.shape,freq.shape)
    for i in range(len(freq)):
      for j in range(len(es)):
#          fo.write(str(h.geometry.r[i,0])+"  ")
          fo.write(str(freq[i])+"  ")
          fo.write(str(es[j])+"  ")
          fo.write(str(np.abs(fcs[j][i]))+"\n")
    fo.close()



