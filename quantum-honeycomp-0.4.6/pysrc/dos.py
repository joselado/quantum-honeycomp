import green
import numpy as np

def dos_surface(h,output_file="DOS.OUT",
                 energies=np.linspace(-1.,1.,20),delta=0.001):
  """Calculates the DOS of a surface, and writes in file"""
  if h.dimensionality!=1: raise # only for 1d
  fo = open(output_file,"w")
  fo.write("# energy, DOS surface, DOS bulk\n")
  for e in energies: # loop over energies
    print "Done",e
    gb,gs = green.green_renormalization(h.intra,h.inter,energy=e,delta=delta)
    gb = -gb.trace()[0,0].imag
    gs = -gs.trace()[0,0].imag
    fo.write(str(e)+"     "+str(gs)+"    "+str(gb)+"\n")
  fo.close()




def dos0d(h,es=[],delta=0.001):
  """Calculate density of states of a 0d system"""
  ds = [] # empty list
  if h.dimensionality==0:  # only for 0d
    iden = np.identity(h.intra.shape[0],dtype=np.complex) # create identity
    for e in es: # loop over energies
      g = ( (e+1j*delta)*iden -h.intra ).I # calculate green function
      ds.append(-g.trace()[0,0].imag)  # add this dos
  else: raise # not implemented...
  write_dos(es,ds)
  return ds


def write_dos(es,ds):
  """ Write DOS in a file"""
  f = open("DOS.OUT","w")
  for (e,d) in zip(es,ds):
    f.write(str(e)+"     ")
    f.write(str(d)+"\n")
  f.close()

