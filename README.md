## QUANTUM HONEYCOMP ##

# Aim #

This program allows to perform tight binding calculations with a user friendly interface


# How to install #

The program has been developed to run in Linux machines. Execute the script "install" to compile all the necesary libraries,
afterwards add to your .bashrc the folder ../bin, where the script quantum-honeycomp is present. 

This program uses several python libraries. For Debian/Ubuntu machines, a script called "autoinstall" will install all the
necessary libraries.

For using this program in Windows, the easiest solution is to create a virtual machine using Virtual Box, installing
a version of ubuntu in that virtual machine, and following the previous instructions.


# Capabilities #
- Tight binding models in different lattices (triangular, square, honeycomb, Kagome, Lieb)
- Tunable parameters in the Hamiltonian (Fermi energy, magnetic order, sublattice imbalance, magnetic field,
  Rashba spin-orbit coupling)
- Different results are ploted from the interface
- Band structure of 0d,1d,2d systems
- Density of states of 0d,1d,2d systems
- Selfconsistent mean field Hubbard calculations of 0d,1d,2d systems
- Berry curvature, Chern number and Z2 invariant in 2d systems
- Special module to deal with systems with 100000 atoms using the Kernel polynomial method
- Special module to study the effect of skyrmion proximity effect in 2d systems

