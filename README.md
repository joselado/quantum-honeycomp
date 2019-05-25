## QUANTUM HONEYCOMP ##

# Aim #

This program allows to perform tight binding calculations with a user friendly interface


# How to install #

The program runs in Linux and Mac machines. 

Download the file "quantum-honeycomp-latest.tar", uncompress it, and execute the script "install". Afterwards, you can execute the program by writting in a terminal "quantum-honeycomp"

This program uses several Python libraries. The simplest way of getting all the dependencies is by installing Python Anaconda from https://www.anaconda.com/distribution/#download-section 

For using this program in Windows, the easiest solution is to create a virtual machine using Virtual Box, installing
a version of ubuntu in that virtual machine, and following the previous instructions.


# Examples
This program allows to study a variety of electronic states by means of tight binding models, in particular

## Quantum anomalous Hall state
Honeycomb lattice with Rashba spin-orbit coupling and exchange field, giving rise to a net Chern number and chiral edge states
https://journals.aps.org/prb/abstract/10.1103/PhysRevB.82.161414
![Alt text](screenshots/qah.png?raw=true "QAH state")



# Capabilities #
- Tight binding models in different lattices (triangular, square, honeycomb, Kagome, Lieb, diamond, pyrochlore)
- Tunable parameters in the Hamiltonian (Fermi energy, magnetic order, sublattice imbalance, magnetic field,  Rashba spin-orbit coupling, intrinsic spin-orbit coupling, Haldane coupling, anti-Haldane coupling, s-wave superconductivity)
- Different results are automatically plotted from the interface
- Band structure of 0d,1d,2d systems
- Density of states of 0d,1d,2d systems
- Selfconsistent mean field Hubbard calculations of 0d,1d,2d systems
- Berry curvature, Chern number and Z2 invariant in 2d systems
- Special module to deal with systems with 100000 atoms using the Kernel polynomial method
- Special modules to study 1d and 2d study interfaces between different systems

