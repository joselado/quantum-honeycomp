import numpy as np
from numba import jit
import os
from .. import filesystem as fs

def multi_ldos(h,es=np.linspace(-2.0,2.0,100),delta=0.05,nrep=1,nk=20,
        ratomic=1.5,dr=0.2,**kwargs):
    """Compute the LDOS at different eenrgies, and add an envelop atomic
    orbital"""
    h = h.copy() # copy the Hamiltonian
    h.turn_dense() # dense hamiltonian
    def get_orbital(r0): # get a localized atomic orbital centered in r0
        def f(r):
            dr = r-r0
            dr2 = np.sum(dr*dr,axis=1) # sum
            return np.exp(-np.sqrt(dr2)/ratomic)
        return f # return the wavefunction
    evals,vs,ks = h.get_eigenvectors(nk=nk,kpoints=True) # compute wavefunctions
    dl = h.geometry.neighbor_directions(nrep+int(ratomic)*10) # directions of the neighbors
    # generate a dictionary with all the real space local orbitals
    ##########################################################
    lodict = dict() # dictionary for the local orbitals
    # get the grids
    x,y = get_grids(h.geometry,nrep=nrep,dr=dr,
            deltax=ratomic*3,deltay=ratomic*3)
    r = np.zeros((len(x),3)) ; r[:,0] = x ; r[:,1] = y
    for d in dl: # loop over directions
          rrep = h.geometry.replicas(d) # replicas in this direction
          for i in range(len(rrep)): # loop over the atoms
              r0 = rrep[i] # get this center
              if h.has_eh:
                if h.has_spin: # spinful
                  lodict[(tuple(d),4*i)] = get_orbital(r0)(r) # store 
                  lodict[(tuple(d),4*i+1)] = get_orbital(r0)(r) # store 
                  lodict[(tuple(d),4*i+2)] = 0. # store 
                  lodict[(tuple(d),4*i+3)] = 0. # store 
                else: raise
              else:
                if h.has_spin: # spinful
                  lodict[(tuple(d),2*i)] = get_orbital(r0)(r) # store 
                  lodict[(tuple(d),2*i+1)] = get_orbital(r0)(r) # store 
                else: # spinless
                  lodict[(tuple(d),i)] = get_orbital(r0)(r) # store 
    ##########################################################
    # now compute the real-space wavefunctions including the Bloch phase
    ds = np.zeros((len(vs),len(x))) # zero array
    for i in range(len(vs)): # loop over wavefunctions
        w = vs[i] # get the current Bloch wavefunction
        k = ks[i] # get the current bloch wavevector
        d = get_real_space_density(w,k,dl,lodict,h.geometry)
        ds[i] = d # store in the list
    # now compute all the LDOS
    fs.rmdir("MULTILDOS")
    fs.mkdir("MULTILDOS")
    fo = open("MULTILDOS/MULTILDOS.TXT","w") # files with the names
    for e in es: # loop over energies
        name0 = "LDOS_"+str(e)+"_.OUT" # name of the output
        name = "MULTILDOS/" + name0
        out = ldos_at_energy(evals,ds,e,delta) # compute the LDOS
        np.savetxt(name,np.array([x,y,out]).T) # save
        fo.write(name0+"\n") # name of the file
    fo.close()
    from ..dos import calculate_dos,write_dos
    es2 = np.linspace(min(es),max(es),len(es)*10)
    ys = calculate_dos(evals,es2,delta,w=None) # compute DOS
    write_dos(es2,ys,output_file="MULTILDOS/DOS.OUT")



def get_real_space_density(w,k,dl,lodict,g):
    """Compute the orbital in real space"""
    nc = len(w) # number of components of the Bloch wavefunction
    nd = len(dl) # number of unit cells to consider
    out = 0. # wavefunction in real space
    for d in dl:
        phi = g.bloch_phase(d,k) # get the Bloch phase
        for i in range(nc):
            out = out + w[i]*phi*lodict[(tuple(d),i)] # add contribution
    return (out*np.conjugate(out)).real # return 



def get_grids(g,nrep=0,dr=0.1,deltax=1.0,deltay=1.0):
    """Return the grids to plot the real space wavefunctions"""
    r = g.multireplicas(nrep) # get all the position
    xmin = np.min(r[:,0])
    xmax = np.max(r[:,0])
    ymin = np.min(r[:,1])
    ymax = np.max(r[:,1])
    nx = int((xmax-xmin+2*deltax)/dr) # number of x points
    ny = int((ymax-ymin+2*deltay)/dr) # number of y points
    xp = np.linspace(xmin-deltax,xmax+deltax,nx) # generate the points
    yp = np.linspace(ymin-deltay,ymax+deltay,ny) # generate the points
    gridx = np.zeros(nx*ny)
    gridy = np.zeros(nx*ny)
    gridx,gridy = get_grids_jit(xp,yp,gridx,gridy)
    return gridx,gridy # return the grids

#@jit
def get_grids_jit(x,y,gridx,gridy):
    nx = len(x)
    ny = len(y)
    k = 0
    for i in range(nx):
      for j in range(ny):
          gridx[k] = x[i]
          gridy[k] = y[j]
          k += 1
    return gridx,gridy

def ldos_at_energy(evals,ds,e,delta):
    """Compute the different local density of states at each energy"""
    de2 = (evals-e)**2 # difference in energy
    out = np.sum(ds.T*delta/(de2+delta**2),axis=1)
    return out # return that density



