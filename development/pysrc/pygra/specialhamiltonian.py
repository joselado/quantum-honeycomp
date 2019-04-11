from . import specialhopping
from . import specialgeometry


def tbg(n=7,ti=0.4,lambi=7.0,is_sparse=True,has_spin=False):
    """
    Return the Hamiltonian of twisted bilayer graphene
    """
    g = specialgeometry.twisted_bilayer(n)
    mgenerator = specialhopping.twisted_matrix(ti=ti,lambi=lambi)
    h = g.get_hamiltonian(is_sparse=is_sparse,has_spin=has_spin,
            is_multicell=True,mgenerator=mgenerator)
    return h




