# library to deal with the spectral properties of the hamiltonian
import numpy as np
import scipy.linalg as lg
import scipy.sparse.linalg as slg
import os
from numba import jit
from .. import filesystem as fs
from .. import parallel
from .. import interpolation


def get_qpi(self,mode="pm",**kwargs):
    if mode=="pm":
        return poor_man_qpi(self,**kwargs)



def poor_man_qpi(h,reciprocal=True,nk=20,energies=np.linspace(-4.0,4.0,80),
        output_folder="MULTIFERMISURFACE",nsuper=3,**kwargs):
    """Compute the QPI using a poor-mans convolution of the k-DOS"""
    if h.dimensionality!=2: raise
    from ..fermisurface import fermi_surface_generator
    if reciprocal: fR = h.geometry.get_k2K_generator() # get matrix
    else:  fR = lambda x: x # get identity
    qs0 = h.geometry.get_kmesh(nk=nk,nsuper=nsuper)
    qs0 = qs0 - np.mean(qs0,axis=0)
    qs = np.array([fR(q) for q in qs0]) # convert
    es,ks,ds = fermi_surface_generator(h,reciprocal=False,info=False,
            energies=energies,
            nsuper=1,nk=2*nk,**kwargs)
    # we now have the energies, k-points and DOS, lets do a convolution
    fs.rmdir(output_folder) # remove folder
    fs.mkdir(output_folder) # create folder
    fp = lambda i: poor_man_qpi_single_energy(ks,ds[:,i],qs) # parallel function
    kdos = parallel.pcall(fp,range(len(es))) # compute in parallel
    kdos = np.array(kdos).T # convert to array
    fo = open(output_folder+"/"+output_folder+".TXT","w")
    for i in range(len(es)): # loop over energies
        filename = output_folder+"_"+str(es[i])+"_.OUT" # name
        name = output_folder+"/"+filename
        np.savetxt(name,np.array([qs0[:,0],qs0[:,1],kdos[:,i]]).T)
        fo.write(filename+"\n")
        name = output_folder+"/DOS.OUT"
    name = "DOS.OUT"
    np.savetxt(name,np.array([es,np.sum(ds,axis=0)]).T)
    fo.close()



def poor_man_qpi_single_energy(ks,ds,qs):
#    return poor_man_qpi_single_energy_brute_force(ks,ds,qs)
    return poor_man_qpi_convolve(ks,ds,qs)



def poor_man_qpi_convolve(ks,ds,qs):
    """Convolve the DOS to simmulate the QPI"""
    nq = int(np.sqrt(len(qs))) # number of qpoints
    nk = int(np.sqrt(len(ks))) # number of qpoints
    grid_kx, grid_ky = np.mgrid[0:1:nq*1j, 0:1:nq*1j] # kx and ky
    # DOS in a grid
    f = interpolation.interpolator2d(ks[:,0],ks[:,1],ds,mode="periodic")
    ksg = np.array([grid_kx, grid_ky]).reshape((2,nq*nq)).T # create points
    dsg = f(ksg).reshape((nq,nq)) # interpolate on a grid
    from ..convolution import selfconvolve
    out = selfconvolve(dsg).reshape((nq*nq)) # do the convolution
    out = interpolation.interpolator2d(ksg[:,0],ksg[:,1],out,mode="periodic")(qs[:,0:2])
    print("Done QPI")
    return out


def poor_man_qpi_single_energy_brute_force(ks,ds,qs):
    """Do the convolution of the Fermi surfaces"""
    f0 = interpolation.interpolator2d(ks[:,0],ks[:,1],ds) # interpolated k-DOS
    def f(k):
        """Define a periodic function"""
        k = k[:,0:2]%1.
        o = f0(k)
        return o
#    return f(qs)
    print("Doing")
    out = [np.mean(f(ks)*f(ks+q)) for q in qs]
    return np.array(out) # return output


