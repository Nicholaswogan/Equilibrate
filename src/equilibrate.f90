module equilibrate
  use equilibrate_const, only: dp, atom_str_len, reac_str_len
  use equilibrate_cea, only: CEAData
  implicit none
  private

  public :: ChemEquiAnalysis, dp, atom_str_len, reac_str_len

  type :: ChemEquiAnalysis
    character(atom_str_len), allocatable :: atoms_names(:)
    character(reac_str_len), allocatable :: species_names(:)
    character(reac_str_len), allocatable :: gas_names(:)
    character(reac_str_len), allocatable :: condensate_names(:)

    real(dp), allocatable, private :: species_composition(:,:)

    real(dp), allocatable :: atom_fractions(:)
    real(dp), allocatable :: mole_fractions(:)
    real(dp), allocatable :: mass_fractions(:)

    real(dp), allocatable :: atom_fractions_gas(:)
    real(dp), allocatable :: mole_fractions_gas(:)
    real(dp), allocatable :: mass_fractions_gas(:)

    real(dp), allocatable :: atom_fractions_condensate(:)
    real(dp), allocatable :: mole_fractions_condensate(:)
    real(dp), allocatable :: mass_fractions_condensate(:)

    real(dp) :: nabla_ad, gamma2, MMW, rho, c_pe

    logical :: verbose = .false.

    !> driver routines
    type(CEAData), allocatable :: dat

  contains
    procedure :: solve
  end type
  interface ChemEquiAnalysis
    module procedure :: create_ChemEquiAnalysis
  end interface

contains
    
  function create_ChemEquiAnalysis(thermopath, atoms, species, err) result(cea)
    character(*), intent(in) :: thermopath
    character(*), optional, intent(in) :: atoms(:)
    character(*), optional, intent(in) :: species(:)
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

    if (allocated(cea%mole_fractions)) then
      deallocate(cea%species_composition)
      deallocate(cea%atom_fractions)
      deallocate(cea%atom_fractions_gas)
      deallocate(cea%atom_fractions_condensate)
      deallocate(cea%mole_fractions)
      deallocate(cea%mole_fractions_gas)
      deallocate(cea%mole_fractions_condensate)
      deallocate(cea%mass_fractions)
      deallocate(cea%mass_fractions_gas)
      deallocate(cea%mass_fractions_condensate)
    endif
    allocate(cea%species_composition(size(cea%atoms_names),size(cea%species_names)))
    allocate(cea%atom_fractions(size(cea%atoms_names)))
    allocate(cea%atom_fractions_gas(size(cea%atoms_names)))
    allocate(cea%atom_fractions_condensate(size(cea%atoms_names)))
    allocate(cea%mole_fractions(size(cea%species_names)))
    allocate(cea%mole_fractions_gas(size(cea%gas_names)))
    allocate(cea%mole_fractions_condensate(size(cea%condensate_names)))
    allocate(cea%mass_fractions(size(cea%species_names)))
    allocate(cea%mass_fractions_gas(size(cea%gas_names)))
    allocate(cea%mass_fractions_condensate(size(cea%condensate_names)))

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

  subroutine solve(self, P, T, X, err)
    class(ChemEquiAnalysis), intent(inout) :: self
    real(dp), intent(in) :: P
    real(dp), intent(in) :: T
    real(dp), intent(in) :: X(:)
    character(:), allocatable, intent(out) :: err

    real(dp), allocatable :: X_(:)
    integer :: i, j, jj

    if (size(X) /= size(self%atoms_names)) then
      err = 'Input "X" has the wrong size'
      return
    endif
    if (any(X < 0)) then
      err = 'Input "X" can not be less than zero'
      return
    endif

    ! normalize
    X_ = X/sum(X)

    call self%dat%easychem(mode='q', verbo='  ', verbose2=self%verbose, N_atoms_in=size(self%atoms_names), &
                           N_reactants_in=size(self%species_names), molfracs_atoms=X_, &
                           molfracs_reactants=self%mole_fractions, &
                           massfracs_reactants=self%mass_fractions, &
                           temp=T, press=P, &
                           nabla_ad=self%nabla_ad, gamma2=self%gamma2, MMW=self%MMW, rho=self%rho, c_pe=self%c_pe)
    if (self%dat%error) then
      err = trim(self%dat%err_msg)
      return
    endif 

    j = 1
    jj = 1
    do i = 1,size(self%species_names)
      if (self%dat%reac_condensed(i)) then
        self%mole_fractions_condensate(j) = self%mole_fractions(i)
        self%mass_fractions_condensate(j) = self%mass_fractions(i)
        j = j + 1
      else
        self%mole_fractions_gas(jj) = self%mole_fractions(i)
        self%mass_fractions_gas(jj) = self%mass_fractions(i)
        jj = jj + 1
      endif
    enddo

    if (size(self%condensate_names) > 0) then
      self%mole_fractions_condensate = self%mole_fractions_condensate/max(sum(self%mole_fractions_condensate),tiny(1.0_dp))
      self%mass_fractions_condensate = self%mass_fractions_condensate/max(sum(self%mass_fractions_condensate),tiny(1.0_dp))
    endif
    if (size(self%gas_names) > 0) then
      self%mole_fractions_gas = self%mole_fractions_gas/max(sum(self%mole_fractions_gas),tiny(1.0_dp))
      self%mass_fractions_gas = self%mass_fractions_gas/max(sum(self%mass_fractions_gas),tiny(1.0_dp))
    endif

    self%atom_fractions = 0.0_dp
    self%atom_fractions_gas = 0.0_dp
    self%atom_fractions_condensate = 0.0_dp
    do j = 1,size(self%species_names)
      self%atom_fractions(:) = self%atom_fractions(:) + self%species_composition(:,j)*self%mole_fractions(j)
      if (self%dat%reac_condensed(j)) then
        self%atom_fractions_condensate(:) = self%atom_fractions_condensate(:) &
                                            + self%species_composition(:,j)*self%mole_fractions(j)
      else
        self%atom_fractions_gas(:) = self%atom_fractions_gas(:) &
                                     + self%species_composition(:,j)*self%mole_fractions(j)
      endif
    enddo
    self%atom_fractions = self%atom_fractions/max(sum(self%atom_fractions),tiny(1.0_dp))
    self%atom_fractions_condensate = self%atom_fractions_condensate/max(sum(self%atom_fractions_condensate),tiny(1.0_dp))
    self%atom_fractions_gas = self%atom_fractions_gas/max(sum(self%atom_fractions_gas),tiny(1.0_dp))

  end subroutine

end module