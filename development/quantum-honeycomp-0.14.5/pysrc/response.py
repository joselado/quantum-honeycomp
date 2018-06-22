import numpy as np



def magnetic_response_map(h,nk=20,nq=20,j=[2.0,0.,0.],r=[0,0,1],
          kp=None):
  """Generate a magnetic susceptibility map"""
  h0 = h.copy() # copy Hamiltonian
# function with the energy
  energy_qvector = energy_qvector_generator(h,j=j,r=r,kp=kp,nk=nk) 
  f = open("SUSCEPTIBILITY.OUT","w")

  k2K = h.geometry.get_k2K_generator() # get the function
  xs = np.linspace(-2.0,2.0,nq) # kx
  ys = np.linspace(-2.0,2.0,nq) # ky

  # loop over qvectors
  for x in xs:
    for y in ys:
      q = np.array([x,y,0.])
      print("Doing",q)
      q = k2K(q)
      q2 = np.sqrt(q.dot(q))
      e = energy_qvector(angle=q2,q=q)  # energy
      f.write(str(x)+"  "+str(y)+"  "+str(e)+"\n")
      f.flush()
  
  f.close()
  print("writen SUSCEPTIBILITY.OUT")




def energy_qvector_generator(h,j=[1.0,0.0,0.0],r=[0,0,1],kp=None,nk=10):
  """Return a function that calculates the total energy
  in a specific qvector"""
  h0 = h.copy() # copy Hamiltonian
  def energy_qvector(random=False,angle=0.0,q=[0.,0.,0.]):
    h = h0.copy()
    h.generate_spin_spiral(vector=r,qspiral=q,angle=angle)
    h.add_zeeman(j)
    return h.total_energy(nkpoints=nk,kp=kp)
  return energy_qvector


