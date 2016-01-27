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







