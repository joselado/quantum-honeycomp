# library to calculate transpot in multiterminal devices

import numpy as np

class device():
  """ Device with leads and scattering part"""
  leads = [] # empty list of leads
  pass



class lead():
  """ Class for a lead"""
  intra = None  # intraterm
  inter = None  # interterm
  coupling = None  # coupling to the center
  def get_green(self,energy,error=0.0001,delta=0.0001):
    """ Get surface green function"""
    import green 
    grb,gr = green.green_renormalization(self.intra,self.inter,error=error,
                                          energy=energy,delta=delta)
    return gr
  def get_selfenergy(self,energy,error=0.0001,delta=0.0001):
    """ Get selfenergy"""
    gr = self.get_green(energy,error=error,delta=delta) # get greenfunction
    t = self.coupling # coupling
    selfenergy = t * gr * t.H 
    return selfenergy


def landauer(d,energy,ij=[(0,1)],error=0.0001,delta=0.0001):
  """ Calculate landauer tranmission between leads i,j """
  # get the selfenergies of the leads
  ss = [l.get_selfenergy(energy,error=error,delta=delta) for l in d.leads]
  # sum of selfenergies
  ssum = ss[0]*0.  
  for s in ss:  ssum += s
  # identity matrix
  iden = np.identity(len(d.center))
  # calculate central green function
  gc = ((energy + delta*1j)*iden -  d.center - ssum). I
  # calculate transmission
  ts = [] # empty list
  for (i,j) in ij: # loop over pairs
    # calculate spectral functions of the leads
    gammai = ss[i] - ss[i].H  
    gammaj = ss[j] - ss[j].H  
    # calculate transmission
    t = (gammai*gc*gammaj.H*gc.H).trace()[0,0].real 
    ts.append(t) # add to the list
  return ts # return list
 





