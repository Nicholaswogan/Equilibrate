module chemequi
  use chemequi_const, only: dp, atom_str_len, reac_str_len
  use chemequi_cea, only: CEAData
  implicit none
  private

  type :: ChemEquiAnalysis
    character(atom_str_len), allocatable :: atoms(:)
    character(reac_str_len), allocatable :: reactants(:)

    !> driver routines
    type(CEAData) :: dat
  end type
  interface ChemEquiAnalysis
    module procedure :: create_ChemEquiAnalysis
  end interface

contains
    
  function create_ChemEquiAnalysis(atoms, reactants, thermopath, err) result(cea)
    character(*), intent(in) :: atoms(:)
    character(*), intent(in) :: reactants(:)
    character(*), intent(in) :: thermopath
    character(:), allocatable, intent(out) :: err
    type(ChemEquiAnalysis) :: cea

    if (size(atoms) < 1 .or. size(reactants) < 1) then
      err = 'atoms and reactants must have size larger than 0'
      return
    endif 
    if (len(atoms(1)) < atom_str_len) then
      err = 'atoms character array must have larger len'
      return
    endif
    if (len(reactants(1)) < reac_str_len) then
      err = 'reactants character array must have larger len'
      return
    endif
    
    cea%atoms = atoms
    cea%reactants = reactants

    call cea%dat%set_data(size(atoms), size(reactants), atoms, reactants, thermopath)
    if (cea%dat%error) then
      err = trim(cea%dat%err_msg)
      return
    endif 

  end function






end module