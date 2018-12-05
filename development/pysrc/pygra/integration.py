
import numpy as np

def integrate_matrix(f,xlim=[0.,1.],eps=0.1,only_imag=False):
  """ Integrates a matrix, the measure is the maximun value of the matrix"""
  return adaptive_simpsons_rule(f,xlim[0],xlim[1],eps,only_imag=only_imag)




def simpsons_rule(f,a,b):
    c = (a+b) / 2.0
    h3 = abs(b-a) / 6.0
    return h3*(f(a) + 4.0*f(c) + f(b))
 
def recursive_asr(f,a,b,eps,whole):
    "Recursive implementation of adaptive Simpson's rule."
    c = (a+b) / 2.0
    left = simpsons_rule(f,a,c)
    right = simpsons_rule(f,c,b)
    if np.max(np.abs(left + right - whole)) <= 15*eps:
        return left + right + (left + right - whole)/15.0
    return recursive_asr(f,a,c,eps/2.0,left) + recursive_asr(f,c,b,eps/2.0,right)


def recursive_asr_imag(f,a,b,eps,whole):
    "Recursive implementation of adaptive Simpson's rule."
    c = (a+b) / 2.0
    left = simpsons_rule(f,a,c)
    right = simpsons_rule(f,c,b)
    if np.max(np.abs(np.imag(left + right - whole))) <= 15*eps:
        return left + right + (left + right - whole)/15.0
    return recursive_asr(f,a,c,eps/2.0,left) + recursive_asr(f,c,b,eps/2.0,right)



def simpson(f,xlim=[0.,1.],eps=0.1):
  """ Integrates a matrix, the measure is the maximun value of the matrix"""
  a = xlim[0]
  b = xlim[1]
  return scalar_asr(f,a,b,eps,simpsons_rule(f,a,b))



def scalar_asr(f,a,b,eps,whole):
    "Recursive implementation of adaptive Simpson's rule."
    c = (a+b) / 2.0
    left = simpsons_rule(f,a,c)
    right = simpsons_rule(f,c,b)
    if np.abs(left + right - whole) <= 15*eps:
        return left + right + (left + right - whole)/15.0
    return scalar_asr(f,a,c,eps/2.0,left) + scalar_asr(f,c,b,eps/2.0,right)



 
def adaptive_simpsons_rule(f,a,b,eps,only_imag=False):
    "Calculate integral of f from a to b with max error of eps."
    if only_imag:
      return recursive_asr_imag(f,a,b,eps,simpsons_rule(f,a,b))
    else:
      return recursive_asr(f,a,b,eps,simpsons_rule(f,a,b))








  



