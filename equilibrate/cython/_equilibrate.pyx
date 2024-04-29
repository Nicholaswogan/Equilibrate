from numpy cimport ndarray, uint8_t
from libc.stdint cimport uintptr_t
from libcpp cimport bool as cbool
cimport ChemEquiAnalysis_pxd as cea_pxd
import numpy as np
import ctypes as ct
import os

DEF S_STR_LEN = 20;
DEF ERR_LEN = 1024;
  
cdef class ChemEquiAnalysis:

  cdef void *_ptr

  def __init__(self, thermofile = None):           
    """Initializes `ChemEquiAnalysis`

    Parameters
    ----------
    thermofile : str
        The input thermodynamic .yaml file
    """
    # Allocate memory
    cea_pxd.allocate_chemequianalysis(&self._ptr)

    # convert strings to char
    cdef bytes thermofile_b = pystring2cstring(thermofile)
    cdef char *thermofile_c = thermofile_b
    cdef char err[ERR_LEN+1]
    
    # Initialize
    cea_pxd.chemequianalysis_create_wrapper(&self._ptr, thermofile_c,
                                            err)
    if len(err.strip()) > 0:
      raise EquilibrateException(err.decode("utf-8").strip())

  def solve(self, double P, double T, molfracs_atoms = None, molfracs_species = None):
    """Computes chemical equilibrium given input atom or species mole fractions.
    If successful, then the equilibrium composition will be stored in a number
    of attributes (e.g., self%molfracs_species).

    Parameters
    ----------
    P : double
      Pressure in bars
    T : double
      Temperature in Kelvin
    molfracs_atoms : ndarray[double,ndim=1], optional
      Atom mole fractions in the same order and length as self.atoms_names.
    molfracs_species : ndarray[double,ndim=1], optional
      Species mole fractions in the same order and length as self.species_names.
    """

    cdef ndarray[double, ndim=1] molfracs_atoms_ = np.empty(1,dtype=np.double)
    cdef cbool molfracs_atoms_present = False
    if molfracs_atoms is not None:
      molfracs_atoms_present = True
      molfracs_atoms_ = molfracs_atoms
    cdef int molfracs_atoms_dim = molfracs_atoms_.shape[0]

    cdef ndarray[double, ndim=1] molfracs_species_ = np.empty(1,dtype=np.double)
    cdef cbool molfracs_species_present = False
    if molfracs_species is not None:
      molfracs_species_present = True
      molfracs_species_ = molfracs_species
    cdef int molfracs_species_dim = molfracs_species_.shape[0]

    cdef char err[ERR_LEN+1]

    cea_pxd.chemequianalysis_solve_wrapper(&self._ptr, &P, &T,
                         &molfracs_atoms_present, &molfracs_atoms_dim, <double *>molfracs_atoms_.data, 
                         &molfracs_species_present, &molfracs_species_dim, <double *>molfracs_species_.data, 
                         err)

    if len(err.strip()) > 0:
      raise EquilibrateException(err.decode("utf-8").strip())

  property atoms_names:
    "List. Names of atoms"
    def __get__(self):
      cdef int dim1
      cea_pxd.chemequianalysis_atoms_names_get_size(&self._ptr, &dim1)
      cdef ndarray arr_c = np.empty(dim1*S_STR_LEN + 1, 'S1')
      cea_pxd.chemequianalysis_atoms_names_get(&self._ptr, &dim1, <char *>arr_c.data)
      return c2stringarr(arr_c, S_STR_LEN, dim1)

  property species_names:
    "List. Names of species"
    def __get__(self):
      cdef int dim1
      cea_pxd.chemequianalysis_species_names_get_size(&self._ptr, &dim1)
      cdef ndarray arr_c = np.empty(dim1*S_STR_LEN + 1, 'S1')
      cea_pxd.chemequianalysis_species_names_get(&self._ptr, &dim1, <char *>arr_c.data)
      return c2stringarr(arr_c, S_STR_LEN, dim1)

  property gas_names:
    "List. Names of gases"
    def __get__(self):
      cdef int dim1
      cea_pxd.chemequianalysis_gas_names_get_size(&self._ptr, &dim1)
      cdef ndarray arr_c = np.empty(dim1*S_STR_LEN + 1, 'S1')
      cea_pxd.chemequianalysis_gas_names_get(&self._ptr, &dim1, <char *>arr_c.data)
      return c2stringarr(arr_c, S_STR_LEN, dim1)

  property condensate_names:
    "List. Names of condensates"
    def __get__(self):
      cdef int dim1
      cea_pxd.chemequianalysis_condensate_names_get_size(&self._ptr, &dim1)
      cdef ndarray arr_c = np.empty(dim1*S_STR_LEN + 1, 'S1')
      cea_pxd.chemequianalysis_condensate_names_get(&self._ptr, &dim1, <char *>arr_c.data)
      return c2stringarr(arr_c, S_STR_LEN, dim1)

  property molfracs_atoms:
    "ndarray[double,ndim=1]. Mole fractions of each atom."
    def __get__(self):
      cdef int dim1
      cea_pxd.chemequianalysis_molfracs_atoms_get_size(&self._ptr, &dim1)
      cdef ndarray arr = np.empty(dim1, np.double)
      cea_pxd.chemequianalysis_molfracs_atoms_get(&self._ptr, &dim1, <double *>arr.data)
      return arr

  property molfracs_species:
    "ndarray[double,ndim=1]. Mole fractions of each species."
    def __get__(self):
      cdef int dim1
      cea_pxd.chemequianalysis_molfracs_species_get_size(&self._ptr, &dim1)
      cdef ndarray arr = np.empty(dim1, np.double)
      cea_pxd.chemequianalysis_molfracs_species_get(&self._ptr, &dim1, <double *>arr.data)
      return arr

  property massfracs_species:
    "ndarray[double,ndim=1]. Mass fractions of each species."
    def __get__(self):
      cdef int dim1
      cea_pxd.chemequianalysis_massfracs_species_get_size(&self._ptr, &dim1)
      cdef ndarray arr = np.empty(dim1, np.double)
      cea_pxd.chemequianalysis_massfracs_species_get(&self._ptr, &dim1, <double *>arr.data)
      return arr

  property molfracs_atoms_gas:
    "ndarray[double,ndim=1]. Mole fractions of atoms in gas phase."
    def __get__(self):
      cdef int dim1
      cea_pxd.chemequianalysis_molfracs_atoms_gas_get_size(&self._ptr, &dim1)
      cdef ndarray arr = np.empty(dim1, np.double)
      cea_pxd.chemequianalysis_molfracs_atoms_gas_get(&self._ptr, &dim1, <double *>arr.data)
      return arr

  property molfracs_species_gas:
    "ndarray[double,ndim=1]. Mole fractions of species in gas phase."
    def __get__(self):
      cdef int dim1
      cea_pxd.chemequianalysis_molfracs_species_gas_get_size(&self._ptr, &dim1)
      cdef ndarray arr = np.empty(dim1, np.double)
      cea_pxd.chemequianalysis_molfracs_species_gas_get(&self._ptr, &dim1, <double *>arr.data)
      return arr

  property molfracs_atoms_condensate:
    "ndarray[double,ndim=1]. Mole fractions of atoms in condensed phase."
    def __get__(self):
      cdef int dim1
      cea_pxd.chemequianalysis_molfracs_atoms_condensate_get_size(&self._ptr, &dim1)
      cdef ndarray arr = np.empty(dim1, np.double)
      cea_pxd.chemequianalysis_molfracs_atoms_condensate_get(&self._ptr, &dim1, <double *>arr.data)
      return arr

  property molfracs_species_condensate:
    "ndarray[double,ndim=1]. Mole fractions of species in condensed phase."
    def __get__(self):
      cdef int dim1
      cea_pxd.chemequianalysis_molfracs_species_condensate_get_size(&self._ptr, &dim1)
      cdef ndarray arr = np.empty(dim1, np.double)
      cea_pxd.chemequianalysis_molfracs_species_condensate_get(&self._ptr, &dim1, <double *>arr.data)
      return arr

  property mass_tol:
    """float. Degree to which mass will be balanced. Gordon & McBride's default 
    is 1.0e-6, but it seems like 1.0e-2 is OK.
    """
    def __get__(self):
      cdef double val
      cea_pxd.chemequianalysis_mass_tol_get(&self._ptr, &val)
      return val
    def __set__(self, double val):
      cea_pxd.chemequianalysis_mass_tol_set(&self._ptr, &val)

# version
cdef extern void equilibrate_version_get(char *version_c)
def _equilibrate_version():
  cdef char version_c[100+1]
  equilibrate_version_get(version_c)
  return version_c.decode("utf-8").strip()
__version__ = _equilibrate_version()

# utils
cdef pystring2cstring(str pystring):
  # add a null c char, and convert to byes
  cdef bytes cstring = (pystring+'\0').encode('utf-8')
  return cstring

cdef c2stringarr(ndarray c_str_arr, int str_len, int arr_len):  
  bs = c_str_arr[:-1].tobytes()
  return [bs[i:i+str_len].decode().strip() for i in range(0, str_len*arr_len, str_len)]

class EquilibrateException(Exception):
    pass