import numpy as np

from .algebra import isnumber

is_number = isnumber

class MultiHopping():
    """Class for a multihopping"""
    def __init__(self,a):
        if type(a)==dict: # dictionary type
            self.dict = a # dictionary
        else: return NotImplemented
    def __add__(self,a):
        if type(a)!=MultiHopping: return NotImplemented
        out = add_hopping_dict(self.dict,a.dict)
        out = MultiHopping(out) # create a new object
        return out
    def __mul__(self,a):
        if type(a)==MultiHopping:
            out = multiply_hopping_dict(self.dict,a.dict)
            out = MultiHopping(out) # create a new object
            return out
        elif isnumber(a):
            out = dict()
            for key in self.dict:
                out[key] = a*self.dict[key]
            out = MultiHopping(out) # create a new object
            return out
        else: return NotImplemented
    def __rmul__(self,a): 
        if isnumber(a): return self*a
        else: return NotImplemented
    def __neg__(self): return (-1)*self
    def __sub__(self,a): return self + (-a)
    def get_dict(self):
        return self.dict # dictionary
    def get_dagger(self):
        out = dict()
        for key in self.dict:
            out[key] = np.conjugate(self.dict[key]).T
        out = MultiHopping(out) # create a new object
        return out





def add_hopping_dict(hop1,hop2):
    """Multiply two hopping dictionaries"""
    out = dict() # create dictionary
    keys = [key1 for key1 in hop1]
    for key2 in hop2:
        if key2 not in keys: keys.append(key2) # store
    for key in keys:
        m = 0 # initialize
        if key in hop1: m = m + hop1[key]
        if key in hop2: m = m + hop2[key]
        out[key] = m # store
    return out



def multiply_hopping_dict(hop1,hop2):
    """Multiply two hopping dictionaries"""
    out = dict() # create dictionary
    for key1 in hop1:
        for key2 in hop2:
            key = tuple(np.array(key1) + np.array(key2))
            m = hop1[key1]@hop2[key2] # multiply
            if np.max(np.abs(m))>1e-6: # if non-zero
                if key in out: out[key] = out[key] + m # add
                else: out[key] = m # store
    return out # return output

