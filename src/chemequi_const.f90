module chemequi_const
  use iso_fortran_env, only: dp => real64
  implicit none
  public

  ! Constants
  real(dp), parameter :: R = 8.3144598d0
  real(dp), parameter :: amu = 1.660538921d-24
  real(dp), parameter :: kB = 1.3806488d-16
  real(dp), parameter :: mol = 6.02214129d23

  !> Atoms & masses, from http://www.science.co.il/PTelements.asp
  integer, parameter      :: N_atoms_save = 104
  character(2), parameter  :: names_atoms_save(N_atoms_save) = &
  (/ 'E ','H ','He','Li','Be','B ','C ','N ','O ','F ','Ne','Na', &
  'Mg','Al','Si','P ','S ','Cl','Ar','K ','Ca','Sc','Ti','V ','Cr','Mn', &
  'Fe','Co','Ni','Cu','Zn','Ga','Ge','As','Se','Br','Kr','Rb','Sr', &
  'Y ','Zr','Nb','Mo','Tc','Ru','Rh','Pd','Ag','Cd','In','Sn','Sb', &
  'Te','I ','Xe','Cs','Ba','La','Ce','Pr','Nd','Pm','Sm','Eu','Gd', &
  'Tb','Dy','Ho','Er','Tm','Yb','Lu','Hf','Ta','W ','Re','Os','Ir', &
  'Pt','Au','Hg','Tl','Pb','Bi','Po','At','Rn','Fr','Ra','Ac','Th', &
  'Pa','U ','Np','Pu','Am','Cm','Bk','Cf','Es','Fm','Md','No','Lr' /)
  real(dp), parameter :: masses_atoms_save(N_atoms_save) = &
  amu*(/ 0.000548579909d0, 1.0079d0,4.0026d0,6.941d0,9.0122d0,10.811d0,12.0107d0,14.0067d0 &
  ,15.9994d0,18.9984d0,20.1797d0,22.9897d0,24.305d0,26.9815d0,28.0855d0,30.9738d0,32.065d0 &
  ,35.453d0,39.948d0,39.0983d0,40.078d0,44.9559d0,47.867d0,50.9415d0,51.9961d0,54.938d0,55.845d0 &
  ,58.9332d0,58.6934d0,63.546d0,65.39d0,69.723d0,72.64d0,74.9216d0,78.96d0,79.904d0,83.8d0,85.4678d0 &
  ,87.62d0,88.9059d0,91.224d0,92.9064d0,95.94d0,98d0,101.07d0,102.9055d0,106.42d0,107.8682d0 &
  ,112.411d0,114.818d0,118.71d0,121.76d0,127.6d0,126.9045d0,131.293d0,132.9055d0,137.327d0,138.9055d0 &
  ,140.116d0,140.9077d0,144.24d0,145d0,150.36d0,151.964d0,157.25d0,158.9253d0,162.5d0,164.9303d0 &
  ,167.259d0,168.9342d0,173.04d0,174.967d0,178.49d0,180.9479d0,183.84d0,186.207d0,190.23d0,192.217d0 &
  ,195.078d0,196.9665d0,200.59d0,204.3833d0,207.2d0,208.9804d0,209d0,210d0,222d0,223d0,226d0,227d0,232.0381d0 &
  ,231.0359d0,238.0289d0,237d0,244d0,243d0,247d0,247d0,251d0,252d0,257d0,258d0,259d0,262d0/)

  ! For thermo arrays
  integer, parameter      :: N_coeffs = 10, N_temps = 10

end module