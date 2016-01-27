


# puts the matrix in spinor form
def m2spin(matin,matin2=[]):
  n=len(matin)
  from numpy import matrix
  if len(matin2)==0:
    matin2=matrix([[0.0j for i in range(n)] for j in range(n)])
  matout=matrix([[0.0j for i in range(2*n)] for j in range(2*n)])
  for i in range(n):
    for j in range(n):
      matout[2*i,2*j]=matin[i,j].copy()
      matout[2*i+1,2*j+1]=matin2[i,j].copy()
  return matout


def spinful(m):
  """ Return a spinful hamiltonian"""
  return m2spin(m,matin2=m)

