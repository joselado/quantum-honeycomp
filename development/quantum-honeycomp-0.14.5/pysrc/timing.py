from __future__ import print_function,division
import time

class Testimator:
  def __init__(self,title=""):
    self.t0 = time.clock() # starting time
    self.title = title
    print(title)
  def remaining(self,i,tot):
    """Print the ramining time in this task"""
    t = time.clock() # current time
    dt = t - self.t0 # difference in time
    out = self.title + " " # empty line
    for j in range(10):
      if j<(i/tot*10): out += "#"
      else: out += " "
    out += str(round(i/tot*100))+"% completed,"
    trem = dt/(i+1)*(tot-i) # remaining time
    out += " remaining time "+str(round(trem,1))+"s"
    out += ", total time "+str(round(dt,1))+"s"
    print(out,end="\r")
