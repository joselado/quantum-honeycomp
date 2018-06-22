# routines to call a function in parallel
from __future__ import print_function
import scipy.linalg as lg

cores = 1 # call in a single by default

def initialize():
  if cores>1:
    from multiprocessing import Pool
    mainpool = Pool(cores) # create pool
  return mainpool



def multieigh(ms):
  """Diagonalize a bunch of Hamiltonians at once"""
  mainpool = initialize()
  return mainpool.map(lg.eigh,ms)





def pcall_serial(fun,args):
  """Function to call in serial"""
  return [fun(a) for a in args]


def pcall_mp(fun,args,cores=1): return pcall_serial(fun,args)

try: # try to use the multiprocessing library
#  from pathos.multiprocessing import Pool
  def pcall_mp(fun,args,cores=cores):
    """Calls a function for every input in args"""
    print("Using",cores,"cores")
    return mainpool.map(fun,args) # return list

except:
  print("Multiprocessing not found, running in a single core")
  def pcall_mp(fun,args,cores=1): return pcall_serial(fun,args)
  





def pcall(fun,args): # define the function
  if cores==1: return pcall_serial(fun,args) # one core, simply iterate
  else: return pcall_mp(fun,args) # call in parallel

