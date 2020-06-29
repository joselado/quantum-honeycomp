import numpy as np
from numba import jit

import scipy.linalg as lg

#use_fortran = False
try:
    from . import dyson2df90
    use_fortran = True
except:
    use_fortran = False

#use_fortran = False

def dyson2d(intra,tx,ty,txy,txmy,nx,ny,nk,ez):
    if use_fortran:
        from . import dyson2df90
        return dyson2df90.dyson2d(intra,tx,ty,txy,txmy,nx,ny,nk,ez)
    else:
        ns = intra.shape[0]*nx*ny
        g = np.zeros((ns,ns),dtype=np.complex)
        return dyson2d_jit(intra,tx,ty,txy,txmy,nx,ny,nk,ez,g)


@jit(nopython=True)
def dyson2d_jit(intra,tx,ty,txy,txmy,nx,ny,nk,ez,g):
    """jit version of the function"""
    n = intra.shape[0] # size of the matrix
    gs = np.zeros((nk*nk,n,n),dtype=np.complex_) # GF in k-points
    gm = np.zeros((nx*2-1,ny*2-1,n,n),dtype=np.complex_) # GF in minicells
    ks = np.zeros((nk*nk,2)) # kpoints
    em = np.identity(n) # identity times energy
    em = em*ez # identity times energy
    ik = 0
    nk2 = nk*nk # total number of kpoints
    # Compute all the GF
    for i in range(nk):
        kx = 1./nk*np.pi*2*i # kpoint
        for j in range(nk):
            ky = 1./nk*np.pi*2*j # kpoint
            m = np.exp(1j*kx)*tx + np.exp(1j*ky)*ty
            m = m + np.exp(1j*(kx+ky))*txy + np.exp(1j*(kx-ky))*txmy
            m = m + np.conjugate(m.T) + intra # Bloch Hamiltonian
            gs[ik,:,:] = np.linalg.inv(em - m) # store GF
            ks[ik,0] = kx
            ks[ik,1] = ky
            ik += 1 # increase counter
    # GF in the different "minicells" of the supercell
    for i in range(-nx+1,nx):
        for j in range(-ny+1,ny):
            m = np.zeros((n,n),dtype=np.complex_) # initialize
            for ik in range(nk2): # loop over kpoints
                m = m + gs[ik,:,:]*np.exp(1j*(ks[ik,0]*i+ks[ik,1]*j))
            gm[nx+i-1,ny+j-1,:,:] = m[:,:]/nk2 # store
    # now create indexes for the minicell
    inds = np.zeros((nx*ny,2),dtype=np.int_)
    ic = 0
    for i in range(nx):
        for j in range(ny):
            inds[ic,0] = i
            inds[ic,1] = j
            ic += 1
    # now compute the supercell GF
    for i in range(nx*ny): # loop over rows
        i1 = inds[i,0] # minicell index
        j1 = inds[i,1] # minicell index
        for j in range(nx*ny):
            i2 = inds[j,0] # minicell index
            j2 = inds[j,1] # minicell index
            ii = i1-i2 + nx -1 
            jj = j1-j2 + ny -1
            m[:,:] = gm[ii,jj,:,:] # get this matrix
            ii0 = n*i
            ii1 = n*(i+1)
            jj0 = n*j
            jj1 = n*(j+1)
            g[ii0:ii1,jj0:jj1] = m[:,:] # store all this data
        #    g[jj0:jj1,ii0:ii1] = m[:,:] # store all this data
    return g






