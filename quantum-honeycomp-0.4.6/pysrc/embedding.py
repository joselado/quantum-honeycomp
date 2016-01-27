# library to perform embedding calculations

import green


def dos_impurity(h,vc=None,energy=0.,mode="adaptive",delta=0.001,nk=200):
  """ Calculates the green function using the embedding technique"""
  if vc==None: vc = h.intra  # assume perfect
  if mode=="adaptive": mode = "adaptative" # stupid lexic mistake
  g,selfe = green.bloch_selfenergy(h,energy=energy,delta=delta,nk=nk,
                                     mode=mode)
  emat = np.matrix(np.identity(len(g)))*(energy + delta*1j)  # E +i\delta 
  gv = (emat - vc -selfe).I   # Green function of a vacancy, with selfener
  ds = (-g.trace()[0,0].imag)  # save DOS of the pristine
  dsv = (-gv.trace()[0,0].imag)  # save DOS of the defected
  class emb_dos: pass
  edos = emb_dos() # create object
  edos.dos_perfect = ds  # pristine dos
  edos.dos_defected = dsv  # defected dos
  return edos # return object





