module equilibrate
  use equilibrate_const, only: dp, atom_str_len, reac_str_len
  use equilibrate_cea, only: CEAData
  implicit none
  private

  public :: ChemEquiAnalysis, dp, atom_str_len, reac_str_len

  type :: ChemEquiAnalysis
    character(atom_str_len), allocatable :: atoms_names(:) !! Names of atoms
    character(reac_str_len), allocatable :: species_names(:) !! Names of species
    character(reac_str_len), allocatable :: gas_names(:) !! Names of gases
    character(reac_str_len), allocatable :: condensate_names(:) !! Names of condensates

    !> Describes the composition of each species (i.e., how many atoms)
    real(dp), allocatable, private :: species_composition(:,:)

    ! All below are results of an equilibrium solve.

    real(dp), allocatable :: molfracs_atoms(:) !! Mole fractions of each atom (size(atoms_names))
    real(dp), allocatable :: molfracs_species(:) !! Mole fractions of each species (size(species_names))
    real(dp), allocatable :: massfracs_species(:) !! Mass fractions of each species (size(species_names))

    real(dp), allocatable :: molfracs_atoms_gas(:) !! Mole fractions of atoms in gas phase (size(atoms_names))
    real(dp), allocatable :: molfracs_species_gas(:) !! Mole fractions of species in gas phase (size(gas_names))

    real(dp), allocatable :: molfracs_atoms_condensate(:) !! Mole fractions of atoms in condensed phase (size(atoms_names))
    real(dp), allocatable :: molfracs_species_condensate(:) !! Mole fractions of species in condensed 
                                                            !! phase (size(condensate_names))

    ! A few undocumented outputs.
    real(dp) :: nabla_ad, gamma2, MMW, rho, c_pe

    logical :: verbose = .false. !! Determines amount of printing.

    !> Driver class
    type(CEAData), allocatable :: dat

  contains
    procedure :: solve
  end type
  interface ChemEquiAnalysis
    module procedure :: create_ChemEquiAnalysis
  end interface

contains
    
  !> Initializes the chemical equilibrium solver given an input thermodynamic file.
  !> The file can have ".inp" or ".yaml" formats. If the file is ".inp" format, then
  !> both `atoms` and `species` must be inputs, describing the atoms and species to 
  !> consider in equilibrium chemistry. If the file is ".yaml" format, then
  !> `atoms` and `species` become optional inputs with the following effects:
  !> - If `atoms` are specified but not `species`, then the code will consider all 
  !>   species with the input `atoms`.
  !> - If `species` are specified but not `atoms`, then the code will consider all
  !>   atoms corresponding to the input `species`.
  !> - If neither `atoms` or `species` are specified, then the code will consider all
  !>   atoms and species in the input file.
  function create_ChemEquiAnalysis(thermopath, atoms, species, err) result(cea)
    character(*), intent(in) :: thermopath !! Path to file describing the thermodynamic 
                                           !! data of each species
    character(*), optional, intent(in) :: atoms(:) !! Names of atoms to include.
    character(*), optional, intent(in) :: species(:) !! Names of species to include.
    character(:), allocatable, intent(out) :: err
    type(ChemEquiAnalysis) :: cea

    integer :: i, j, jj, k

    if (present(atoms)) then
      if (size(atoms) < 1) then
        err = 'atoms and species must have size larger than 0'
        return
      endif
      if (len(atoms(1)) > atom_str_len) then
        err = 'atoms character array must have larger len'
        return
      endif
    endif
    if (present(species)) then
      if (size(species) < 1) then
        err = 'atoms and species must have size larger than 0'
        return
      endif 
      if (len(species(1)) > reac_str_len) then
        err = 'species character array must have larger len'
        return
      endif
    endif

    if (allocated(cea%dat)) deallocate(cea%dat)
    allocate(cea%dat)

    call cea%dat%set_data(thermopath, atoms, species)
    if (cea%dat%error) then
      err = trim(cea%dat%err_msg)
      return
    endif

    cea%atoms_names = cea%dat%names_atoms(1:size(cea%dat%names_atoms)-1) ! skip E
    cea%species_names = cea%dat%names_reactants_orig

    allocate(cea%condensate_names(cea%dat%N_cond))
    allocate(cea%gas_names(cea%dat%N_gas))
    j = 1
    jj = 1
    do i = 1,size(cea%species_names)
      if (cea%dat%reac_condensed(i)) then
        cea%condensate_names(j) = cea%species_names(i)
        j = j + 1
      else
        cea%gas_names(jj) = cea%species_names(i)
        jj = jj + 1
      endif
    enddo

    if (allocated(cea%molfracs_species)) then
      deallocate(cea%species_composition)
      deallocate(cea%molfracs_atoms)
      deallocate(cea%molfracs_atoms_gas)
      deallocate(cea%molfracs_atoms_condensate)
      deallocate(cea%molfracs_species)
      deallocate(cea%molfracs_species_gas)
      deallocate(cea%molfracs_species_condensate)
      deallocate(cea%massfracs_species)
    endif
    allocate(cea%species_composition(size(cea%atoms_names),size(cea%species_names)))
    allocate(cea%molfracs_atoms(size(cea%atoms_names)))
    allocate(cea%molfracs_atoms_gas(size(cea%atoms_names)))
    allocate(cea%molfracs_atoms_condensate(size(cea%atoms_names)))
    allocate(cea%molfracs_species(size(cea%species_names)))
    allocate(cea%molfracs_species_gas(size(cea%gas_names)))
    allocate(cea%molfracs_species_condensate(size(cea%condensate_names)))
    allocate(cea%massfracs_species(size(cea%species_names)))

    ! Get composition of each species
    cea%species_composition = 0.0_dp
    do j = 1,size(cea%species_names)
      do k = 1,size(cea%dat%reac_atoms_id,1)
        if (cea%dat%reac_atoms_id(k,j) > 0) then
          i = findloc(cea%dat%id_atoms, cea%dat%reac_atoms_id(k, j), 1)
          if (i == 0 .or. i > size(cea%atoms_names) + 1) then
            err = 'Indexing error during initialization.'
            return
          endif
          if (i == size(cea%atoms_names) + 1) then
            ! Electron so we skip
            cycle
          endif
          cea%species_composition(i,j) = cea%dat%reac_stoich(k,j)
        endif
      enddo
    enddo

  end function

  !> Computes chemical equilibrium given input atom or species mole fractions.
  !> If successful, then the equilibrium composition will be stored in a number
  !> of attributes (e.g., self%molfracs_species). If unsuccessful, then the variable
  !> `err` is allocated with an error message.
  subroutine solve(self, P, T, molfracs_atoms, molfracs_species, err)
    class(ChemEquiAnalysis), intent(inout) :: self
    real(dp), intent(in) :: P !! Pressure in bars
    real(dp), intent(in) :: T !! Temperature in Kelvin
    !> Atom mole fractions in the same order and length as self%atoms_names.
    real(dp), optional, intent(in) :: molfracs_atoms(:)
    !> Species mole fractions in the same order and length as self%species_names.
    real(dp), optional, intent(in) :: molfracs_species(:)
    character(:), allocatable, intent(out) :: err

    real(dp), allocatable :: molfracs_atoms_(:)
    integer :: i, j, jj

    if (present(molfracs_atoms) .and. present(molfracs_species)) then
      err = 'Both "molfracs_atoms" and "molfracs_species" are inputs, but only one is allowed.'
      return
    endif
    if (.not.present(molfracs_atoms) .and. .not.present(molfracs_species)) then
      err = 'Neither "molfracs_atoms" or "molfracs_species" are inputs, but one is needed.'
      return
    endif

    ! Input is an array of atom mole fractions
    if (present(molfracs_atoms)) then
      if (size(molfracs_atoms) /= size(self%atoms_names)) then
        err = 'Input "molfracs_atoms" has the wrong size'
        return
      endif
      if (any(molfracs_atoms < 0)) then
        err = 'Input "molfracs_atoms" can not be less than zero'
        return
      endif
      molfracs_atoms_ = molfracs_atoms/max(sum(molfracs_atoms),tiny(1.0_dp))
    endif

    ! Input is an array of species mole fractions
    if (present(molfracs_species)) then
      if (size(molfracs_species) /= size(self%species_names)) then
        err = 'Input "molfracs_species" has the wrong size'
        return
      endif
      if (any(molfracs_species < 0)) then
        err = 'Input "molfracs_species" can not be less than zero'
        return
      endif

      allocate(molfracs_atoms_(size(self%atoms_names)))
      molfracs_atoms_ = 0.0_dp
      do j = 1,size(self%species_names)
        molfracs_atoms_(:) = molfracs_atoms_(:) + self%species_composition(:,j)*molfracs_species(j)
      enddo
      molfracs_atoms_ = molfracs_atoms_/max(sum(molfracs_atoms_),tiny(1.0_dp))
    endif

    ! Compute chemical equilibrium
    call self%dat%solve(mode='q', verbo='  ', verbose2=self%verbose, N_atoms_in=size(self%atoms_names), &
                           N_reactants_in=size(self%species_names), molfracs_atoms=molfracs_atoms_, &
                           molfracs_reactants=self%molfracs_species, &
                           massfracs_reactants=self%massfracs_species, &
                           temp=T, press=P, &
                           nabla_ad=self%nabla_ad, gamma2=self%gamma2, MMW=self%MMW, rho=self%rho, c_pe=self%c_pe)
    if (self%dat%error) then
      err = trim(self%dat%err_msg)
      return
    endif 

    ! Compute mole fractions of gases and condensates in the solution
    j = 1
    jj = 1
    do i = 1,size(self%species_names)
      if (self%dat%reac_condensed(i)) then
        self%molfracs_species_condensate(j) = self%molfracs_species(i)
        j = j + 1
      else
        self%molfracs_species_gas(jj) = self%molfracs_species(i)
        jj = jj + 1
      endif
    enddo
    if (size(self%condensate_names) > 0) then
      self%molfracs_species_condensate = self%molfracs_species_condensate/max(sum(self%molfracs_species_condensate),tiny(1.0_dp))
    endif
    if (size(self%gas_names) > 0) then
      self%molfracs_species_gas = self%molfracs_species_gas/max(sum(self%molfracs_species_gas),tiny(1.0_dp))
    endif

    ! Compute the mole fractions of the atoms in the solution
    self%molfracs_atoms = 0.0_dp
    self%molfracs_atoms_gas = 0.0_dp
    self%molfracs_atoms_condensate = 0.0_dp
    do j = 1,size(self%species_names)
      self%molfracs_atoms(:) = self%molfracs_atoms(:) + self%species_composition(:,j)*self%molfracs_species(j)
      if (self%dat%reac_condensed(j)) then
        self%molfracs_atoms_condensate(:) = self%molfracs_atoms_condensate(:) &
                                            + self%species_composition(:,j)*self%molfracs_species(j)
      else
        self%molfracs_atoms_gas(:) = self%molfracs_atoms_gas(:) &
                                     + self%species_composition(:,j)*self%molfracs_species(j)
      endif
    enddo
    self%molfracs_atoms = self%molfracs_atoms/max(sum(self%molfracs_atoms),tiny(1.0_dp))
    self%molfracs_atoms_condensate = self%molfracs_atoms_condensate/max(sum(self%molfracs_atoms_condensate),tiny(1.0_dp))
    self%molfracs_atoms_gas = self%molfracs_atoms_gas/max(sum(self%molfracs_atoms_gas),tiny(1.0_dp))

  end subroutine

end module