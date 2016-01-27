# library with specific function for hybrid ribbons

import numpy as np


def project_bands(h):
  """ Calculates the bands for different projections
    in a hybrid ribbon """
  from os import system
  import input_tb90 as in90
  # upper edge
  in90.write_operator(h,"upper_edge")   
  system("tb90.x")
  m_upper = np.genfromtxt("BANDS.OUT")
  # lower edge
  in90.write_operator(h,"lower_edge")   
  system("tb90.x")
  m_lower = np.genfromtxt("BANDS.OUT")
  # interface
  in90.write_operator(h,"interface")   
  system("tb90.x")
  m_inter = np.genfromtxt("BANDS.OUT")
  # now create the new file with everything
  mout = []
  mout += []




def create_hybrid(h1,h2,coupling=1.0):
  """ Creates a hybrid hamiltonian from the two inputs"""
  from copy import deepcopy
  hy = deepcopy(h1) # copy hamiltonian
  # check that thw two hamiltonians are mixable
  if not h1.has_spin == h2.has_spin:  # spinpol
    raise
  if not h1.has_eh == h2.has_eh:  # electron hole
    raise
  if not len(h1.intra) == len(h2.intra):  # same dimension
    print "Wrong dimensions", len(h1.intra),len(h2.intra)
    raise
  if not h1.dimensionality == h2.dimensionality:  # dimension
    raise
  dd = len(h1.intra) # dimension of the hamiltonian
  # we wnt a matrix like
  # ( h1_1   t12  )
  # ( t21   h2_2  )
  # substitute the block h2_2 by the second block in input h2
  for i in range(dd/2,dd):
    for j in range(dd/2,dd):
      hy.intra[i,j] = h2.intra[i,j] # substitute element
      hy.inter[i,j] = h2.inter[i,j] # substitute element
  # and now create t12 by mixing 
  for i in range(dd/2):
    for j in range(dd/2,dd):
      hy.intra[i,j] = coupling*(h1.intra[i,j]+h2.intra[i,j])/2.0 # substitute element
      hy.inter[i,j] = coupling*(h1.inter[i,j]+h2.inter[i,j])/2.0 # substitute element
  for i in range(dd/2,dd):
    for j in range(dd/2):
      hy.intra[i,j] = coupling*(h1.intra[i,j]+h2.intra[i,j])/2.0 # substitute element
      hy.inter[i,j] = coupling*(h1.inter[i,j]+h2.inter[i,j])/2.0 # substitute element
  hy.check()  # check that everything is fine
  return hy  # return the hamiltonian



