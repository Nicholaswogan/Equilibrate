module chemequi
  use chemequi_const, only: dp, atom_str_len, reac_str_len
  use chemequi_cea, only: CEAData
  implicit none
  private

  public :: ChemEquiAnalysis, dp, atom_str_len, reac_str_len

  type :: ChemEquiAnalysis
    character(atom_str_len), allocatable :: atoms(:)
    character(reac_str_len), allocatable :: reactants(:)

    real(dp), allocatable :: molfracs_reactants(:)
    real(dp), allocatable :: massfracs_reactants(:)
    real(dp) :: nabla_ad, gamma2, MMW, rho, c_pe

    !> driver routines
    type(CEAData), allocatable :: dat

  contains
    procedure :: solve
  end type
  interface ChemEquiAnalysis
    module procedure :: create_ChemEquiAnalysis
  end interface

contains
    
  function create_ChemEquiAnalysis(thermopath, atoms, reactants, err) result(cea)
    character(*), intent(in) :: thermopath
    character(*), optional, intent(in) :: atoms(:)
    character(*), optional, intent(in) :: reactants(:)
    character(:), allocatable, intent(out) :: err
    type(ChemEquiAnalysis) :: cea

    if (present(atoms)) then
      if (size(atoms) < 1) then
        err = 'atoms and reactants must have size larger than 0'
        return
      endif
      if (len(atoms(1)) > atom_str_len) then
        err = 'atoms character array must have larger len'
        return
      endif
    endif
    if (present(reactants)) then
      if (size(reactants) < 1) then
        err = 'atoms and reactants must have size larger than 0'
        return
      endif 
      if (len(reactants(1)) > reac_str_len) then
        err = 'reactants character array must have larger len'
        return
      endif
    endif

    if (allocated(cea%dat)) deallocate(cea%dat)
    allocate(cea%dat)

    call cea%dat%set_data(thermopath, atoms, reactants)
    if (cea%dat%error) then
      err = trim(cea%dat%err_msg)
      return
    endif 

    cea%atoms = cea%dat%names_atoms(1:size(cea%dat%names_atoms)-1) ! skip E
    cea%reactants = cea%dat%names_reactants_orig
    if (allocated(cea%molfracs_reactants)) then
      deallocate(cea%molfracs_reactants)
      deallocate(cea%massfracs_reactants)
    endif
    allocate(cea%molfracs_reactants(size(cea%reactants)))
    allocate(cea%massfracs_reactants(size(cea%reactants)))

  end function

  subroutine solve(self, P, T, X, err)
    class(ChemEquiAnalysis), intent(inout) :: self
    real(dp), intent(in) :: P
    real(dp), intent(in) :: T
    real(dp), intent(in) :: X(:)
    character(:), allocatable, intent(out) :: err

    real(dp), allocatable :: X_(:)

    if (size(X) /= size(self%atoms)) then
      err = 'Input "X" has the wrong size'
      return
    endif
    if (any(X < 0)) then
      err = 'Input "X" can not be less than zero'
      return
    endif

    ! normalize
    X_ = X/sum(X)

    call self%dat%easychem(mode='q', verbo='  ', N_atoms_in=size(self%atoms), &
                           N_reactants_in=size(self%reactants), molfracs_atoms=X_, &
                           molfracs_reactants=self%molfracs_reactants, &
                           massfracs_reactants=self%massfracs_reactants, &
                           temp=T, press=P, &
                           nabla_ad=self%nabla_ad, gamma2=self%gamma2, MMW=self%MMW, rho=self%rho, c_pe=self%c_pe)
    if (self%dat%error) then
      err = trim(self%dat%err_msg)
      return
    endif 

  end subroutine

end module