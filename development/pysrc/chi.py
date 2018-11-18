import numpy as np
import scipy.linalg as lg
try:
    import chif90
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



