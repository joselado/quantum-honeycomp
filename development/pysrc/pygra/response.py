import numpy as np
from . import timing



def magnetic_response_map(h,nk=20,nq=20,j=[0.1,0.,0.],r=[0,0,1],
          kp=None,qs=None):
  """Generate a magnetic susceptibility map"""
  h0 = h.copy() # copy Hamiltonian
# function with the energy
  energy_qvector = energy_qvector_generator(h,j=j,r=r,kp=kp,nk=nk) 
  f = open("SUSCEPTIBILITY.OUT","w")

  k2K = h.geometry.get_k2K_generator() # get the function
  # loop over qvectors
  if qs is None: # not provided
    xs = np.linspace(-4.0,4.0,nq) # kx
    ys = np.linspace(-4.0,4.0,nq) # ky
    qs = [] # initialize qvectors
    for x in xs:
      for y in ys:
        qs.append(np.array([x,y,0.])) # store
  else: pass # qs provided on input
  est = timing.Testimator(maxite=len(qs)) # initialize estimator
  for q0 in qs: # loop
    est.iterate()
    q = k2K(q0) # put in reciprocal coordinates
    q2 = np.sqrt(q.dot(q))
    e = energy_qvector(angle=q2,q=q)  # energy
    f.write(str(q0[0])+"  "+str(q0[1])+"  "+str(e)+"\n")
    f.flush()
  
  f.close()
  print("writen SUSCEPTIBILITY.OUT")




def energy_qvector_generator(h,j=[1.0,0.0,0.0],r=[0,0,1],kp=None,nk=10):
  """Return a function that calculates the total energy
  in a specific qvector"""
  h0 = h.copy() # copy Hamiltonian
  def energy_qvector(random=False,angle=0.0,q=[0.,0.,0.]):
    h = h0.copy()
    h.generate_spin_spiral(vector=r,qspiral=2.*q,angle=angle)
    h.add_zeeman(j)
    return h.total_energy(nkpoints=nk,kp=kp)
  return energy_qvector


