# library to verify the type of Hamiltonian

def check_mode(h,n):
    if n=="spinless_nambu":
        if (not h.has_spin) and h.has_eh: return True
        else: return False
    elif n=="spinful_nambu":
        if h.has_spin and h.has_eh: return True
        else: return False
    elif n=="spinful":
        if h.has_spin and (not h.has_eh): return True
        else: return False
    elif n=="spinless":
        if (not h.has_spin) and (not h.has_eh): return True
        else: return False
    else: raise
