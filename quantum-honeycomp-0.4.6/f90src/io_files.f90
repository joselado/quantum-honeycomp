! this file defines the name of the input and output files
module io_files
! file for the hamiltonian
character (len=80) :: hamiltonian_file='hamiltonian.in'   
character (len=80) :: magnetism_file='MAGNETIZATION.OUT'   
character (len=80) :: density_file='DENSITY.OUT'   
character (len=80) :: bands_file='BANDS.OUT'   
character (len=80) :: dos_file='DOS.OUT'   
character (len=80) :: berry_phase_file='BERRY_PHASE.OUT'   
character (len=80) :: berry_curv_file='BERRY_CURVATURE.OUT'   
character (len=80) :: denmat_file='DENSITY_MATRIX.OUT'   


end module io_files
