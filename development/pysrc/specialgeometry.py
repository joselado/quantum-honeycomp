from __future__ import print_function,division
import numpy as np
import geometry
import sculpt
import geometry


def twisted_bilayer(m0=3,rotate=True,shift=[0.,0.],center="AB/BA",
  sublattice=True,r=1):
  """Return the geometry for twisted bilayer graphene"""
  g = geometry.honeycomb_lattice()
  g.has_sublattice = False
  if sublattice: # trace the sublattice using a dirty trick
    g.z[0] += 0.001
    g.z[1] -= 0.001
    g.xyz2r()
  else: pass
  g = geometry.non_orthogonal_supercell(g,m=[[-1,0,0],[0,1,0],[0,0,1]])
#  m0 = 3
#  r = 1
  theta = np.arccos((3.*m0**2+3*m0*r+r**2/2.)/(3.*m0**2+3*m0*r+r**2))
  print("Theta",theta*180.0/np.pi)
  nsuper = [[m0,m0+r,0],[-m0-r,2*m0+r,0],[0,0,1]]
  g = geometry.non_orthogonal_supercell(g,m=nsuper,
           reducef=lambda x: 3*np.sqrt(x))
  g1 = g.copy()
  g1.shift([1.,0.,0.]) 
  g.z -= 1.5
  g.xyz2r() # update
  if rotate: # rotate one of the layers
    g1 = g1.rotate(theta*180/np.pi)
    g1s = g1.supercell(2) # supercell
    g1s.z += 1.5
    g1s.x += shift[0] # shift upper layer
    g1s.y += shift[1] # shift upper layer
    g1s.xyz2r() # update
    g1.a1 = g.a1
    g1.a2 = g.a2
    rs = sculpt.retain_unit_cell(g1s.r,g.a1,g.a2,g.a3,dim=2) # get positions
  else: # do not rotate
    rs = np.array(g1.r) # same coordinates
    rs[:,2] = 1.0 # set as one
  g1.r = np.concatenate([rs,g.r])
 
  g1.r2xyz() # update
  g1.real2fractional() # update fractional coordinates 
  if center=="AB/BA": pass # do nothing 
  elif center=="AA": 
    g1.frac_x = (g1.frac_x)%1 # to the unit cell
    g1.fractional2real()
  elif center=="AB": 
    g1.frac_x = (g1.frac_x)%1 # to the unit cell
    g1.frac_y = (g1.frac_y)%1 # to the unit cell
    g1.frac_y += 0.5  # to the unit cell
    g1.frac_x = (g1.frac_x)%1 # to the unit cell
    g1.frac_y = (g1.frac_y)%1 # to the unit cell
    g1.frac_x -= 1./3.
    g1.frac_y -= 1./3.
    g1.frac_x = (g1.frac_x)%1 # to the unit cell
    g1.frac_y = (g1.frac_y)%1 # to the unit cell
    g1.fractional2real()
  else: raise
  if sublattice: # recover the sublattice
    g1.has_sublattice = True # say that it has
    sl = []
    for r in g1.r: # loop over positions
      if np.abs(r[2]-1.5)<0.01: # upper layer
        if r[2]-1.5>0.0: sl.append(-1.) # A sublattice
        else: sl.append(1.) # B sublattice
      elif np.abs(r[2]+1.5)<0.01: # lower layer
        if r[2]+1.5>0.0: sl.append(-1.) # A sublattice
        else: sl.append(1.) # B sublattice
    g1.z = np.round(g1.z,2) # delete the small shift
    g1.xyz2r() # update coordinates
    g1.sublattice = np.array(sl) # store sublattice
  g1 = sculpt.rotate_a2b(g1,g1.a1,np.array([1.,0.,0.])) # rotate
  g1.get_fractional() # get fractional coordinates 
  return g1



def twisted_multilayer(m0=3,rotate=True,shift=[0.,0.],
  sublattice=True,r=1,rot=[1,1,0,0]):
  """Return the geometry for twisted bilayer graphene"""
  g = geometry.honeycomb_lattice()
  g.has_sublattice = False # no sublattice
  g = geometry.non_orthogonal_supercell(g,m=[[-1,0,0],[0,1,0],[0,0,1]])
  theta = np.arccos((3.*m0**2+3*m0*r+r**2/2.)/(3.*m0**2+3*m0*r+r**2))
  print("Theta",theta*180.0/np.pi)
  nsuper = [[m0,m0+r,0],[-m0-r,2*m0+r,0],[0,0,1]]
  g = geometry.non_orthogonal_supercell(g,m=nsuper,
           reducef=lambda x: 3*np.sqrt(x))
  if rotate: # rotate one of the layers
    gs = [] # empty list with geometries
    ii = 0
    for i in rot: # loop
        print(i)
        if i!=0 and i!=1: raise # nope
        gs.append(rotate_layer(g,i*theta,dr=[0.,0.,3.0*(ii-len(rot)/2.+.5)]))
        ii += 1 # increase counter
  else: raise
#  g.r = np.concatenate([g1.r,g.r,g2.r]).copy()
  g.r = np.concatenate([gi.r for gi in gs]).copy() # all the positions
#  g.r = np.concatenate([g2.r,g.r]).copy()
 # g.r = g1.r
  g.r2xyz() # update
  g.real2fractional() # update fractional coordinates 
  g = sculpt.rotate_a2b(g,g.a1,np.array([1.,0.,0.])) # rotate
  g.get_fractional() # get fractional coordinates 
  return g


def twisted_trilayer(m0=3):
  return twisted_multilayer(m0=m0,shift=[-1.,0.,1.])


def rotate_layer(g,theta,dr=[0.,0.,0.]):
    """Rotate one of the layers, reatinign the same unit cell"""
    g1 = g.copy() # copy geometry
    g1 = g1.rotate(theta*180/np.pi)
    g1s = g1.supercell(2) # supercell
#    g1.a1 = g.a1.copy()
#    g1.a2 = g.a2.copy()
    g1.r = sculpt.retain_unit_cell(g1s.r,g.a1,g.a2,g.a3,dim=2) # get positions
    g1.r2xyz() # update
    g1.x += dr[0] # shift upper layer
    g1.y += dr[1] # shift upper layer
    g1.z += dr[2] # shift upper layer
    g1.xyz2r() # update
    return g1 # return the geometry


def multilayer_graphene(l=[0]):
    """Return a multilayer graphene system"""
    g = geometry.honeycomb_lattice() # get a honeycomb lattice
    dr = g.r[0] - g.r[1] # shift between two sites
    ss = [] # list for the sublattice
    rs = [] # list for the positions
    ii = 0 # start
    for il in l: # loop over layers
        for (ri,si) in zip(g.r,g.sublattice): # loop over positions
            rs.append(ri + dr*il + ii*np.array([0.,0.,3.])) # add position
            ss.append(si*(-1)**(il+ii)) # add sublattice
        ii += 1
    go = g.copy()
    go.r = np.array(rs)
    go.sublattice = np.array(ss)
    go.r2xyz()
    go.center()
    go = sculpt.rotate_a2b(go,go.a1,np.array([1.,0.,0.]))
    return go




