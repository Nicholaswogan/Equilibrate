from libcpp cimport bool as cbool
cdef extern from "<stdbool.h>":
  pass

# Allocate and destroy
cdef extern void allocate_chemequianalysis(void *ptr);
cdef extern void deallocate_chemequianalysis(void *ptr);

# Wrappers for functions

cdef extern void chemequianalysis_create_wrapper(void *ptr, char *thermofile, char *err);
cdef extern void chemequianalysis_solve_wrapper(void *ptr, double *P, double *T,
                         cbool *molfracs_atoms_present, int *molfracs_atoms_dim, double *molfracs_atoms, 
                         cbool *molfracs_species_present, int *molfracs_species_dim, double *molfracs_species, 
                         char *err)

# Getters and setters

cdef extern void chemequianalysis_atoms_names_get_size(void *ptr, int *dim1)
cdef extern void chemequianalysis_atoms_names_get(void *ptr, int *dim1, char* arr)

cdef extern void chemequianalysis_species_names_get_size(void *ptr, int *dim1)
cdef extern void chemequianalysis_species_names_get(void *ptr, int *dim1, char* arr)

cdef extern void chemequianalysis_gas_names_get_size(void *ptr, int *dim1)
cdef extern void chemequianalysis_gas_names_get(void *ptr, int *dim1, char* arr)

cdef extern void chemequianalysis_condensate_names_get_size(void *ptr, int *dim1)
cdef extern void chemequianalysis_condensate_names_get(void *ptr, int *dim1, char* arr)

cdef extern void chemequianalysis_molfracs_atoms_get_size(void *ptr, int *dim1)
cdef extern void chemequianalysis_molfracs_atoms_get(void *ptr, int *dim1, double *arr)

cdef extern void chemequianalysis_molfracs_species_get_size(void *ptr, int *dim1)
cdef extern void chemequianalysis_molfracs_species_get(void *ptr, int *dim1, double *arr)

cdef extern void chemequianalysis_massfracs_species_get_size(void *ptr, int *dim1)
cdef extern void chemequianalysis_massfracs_species_get(void *ptr, int *dim1, double *arr)

cdef extern void chemequianalysis_molfracs_atoms_gas_get_size(void *ptr, int *dim1)
cdef extern void chemequianalysis_molfracs_atoms_gas_get(void *ptr, int *dim1, double *arr)

cdef extern void chemequianalysis_molfracs_species_gas_get_size(void *ptr, int *dim1)
cdef extern void chemequianalysis_molfracs_species_gas_get(void *ptr, int *dim1, double *arr)

cdef extern void chemequianalysis_molfracs_atoms_condensate_get_size(void *ptr, int *dim1)
cdef extern void chemequianalysis_molfracs_atoms_condensate_get(void *ptr, int *dim1, double *arr)

cdef extern void chemequianalysis_molfracs_species_condensate_get_size(void *ptr, int *dim1)
cdef extern void chemequianalysis_molfracs_species_condensate_get(void *ptr, int *dim1, double *arr)

cdef extern void chemequianalysis_mass_tol_get(void *ptr, double *val)
cdef extern void chemequianalysis_mass_tol_set(void *ptr, double *val)