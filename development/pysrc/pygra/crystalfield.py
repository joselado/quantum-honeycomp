import numpy as np
from scipy.sparse import csc_matrix


def hartree(h,v=0.0):
    """Add a crystal field to a Hamiltonian"""
    if v==0.0: return
    m = hartree_onsite(h.geometry,vc=v) # get array
    m = m - np.mean(m) # remove average
    nat = len(h.geometry.r) # number of atoms
    ind = range(nat) # indexes
    mat = csc_matrix((m,(ind,ind)),shape=(nat,nat),dtype=np.complex)
    h.intra = h.intra + h.spinless2full(mat) # add contribution






def hartree_onsite(g,rcut=6.0,vc=0.0):
    """Return an array with the Hartree terms"""
    g = g.copy() # copy geometry
    interactions = [] # empty list
    nat = len(g.r) # number of atoms
    mout = np.zeros(nat) # initialize matrix
    lat = np.sqrt(g.a1.dot(g.a1)) # size of the unit cell
    g.ncells = int(2*rcut/lat)+1 # number of unit cells to consider
    ri = g.r # positions
    for d in g.neighbor_directions(): # loop
        rj = np.array(g.replicas(d)) # positions
        for i in range(nat):
            dx = rj[:,0] - ri[i,0]
            dy = rj[:,1] - ri[i,1]
            dz = rj[:,2] - ri[i,2]
            dr = dx*dx + dy*dy + dz*dz
            dr = np.sqrt(dr) # square root
            dr[dr<1e-4] = 1e10
            v = vc/dr # Coulomb interaction
            v[dr<1e-4] = 0.0
            v[dr>rcut] = 0.0
            v = v*np.exp(-dr/rcut) # quench interaction
            mout[i] += np.sum(v) # store contribution
    return mout # return








