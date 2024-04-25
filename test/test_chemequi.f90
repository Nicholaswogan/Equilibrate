program main
  use chemequi, only: ChemEquiAnalysis, dp, reac_str_len, atom_str_len
  implicit none

  call test()

contains

  subroutine test()
    type(ChemEquiAnalysis) :: cea
    character(:), allocatable :: err
    character(reac_str_len), allocatable :: reactants(:)
    character(atom_str_len), allocatable :: atoms(:)
    real(dp), allocatable :: X(:), correct_answer(:)

    integer :: i

    reactants = [ &
      'H              ', &
      'H2             ', &
      'He             ', &
      'O              ', &
      'C              ', &
      'N              ', &
      'Mg             ', &
      'Si             ', &
      'Fe             ', &
      'S              ', &
      'Al             ', &
      'Ca             ', &
      'Na             ', &
      'Ni             ', &
      'P              ', &
      'K              ', &
      'Ti             ', &
      'CO             ', &
      'OH             ', &
      'SH             ', &
      'N2             ', &
      'O2             ', &
      'SiO            ', &
      'TiO            ', &
      'SiS            ', &
      'H2O            ', &
      'C2             ', &
      'CH             ', &
      'CN             ', &
      'CS             ', &
      'SiC            ', &
      'NH             ', &
      'SiH            ', &
      'NO             ', &
      'SN             ', &
      'SiN            ', &
      'SO             ', &
      'S2             ', &
      'C2H            ', &
      'HCN            ', &
      'C2H2,acetylene ', &
      'CH4            ', &
      'AlH            ', &
      'AlOH           ', &
      'Al2O           ', &
      'CaOH           ', &
      'MgH            ', &
      'MgOH           ', &
      'PH3            ', &
      'CO2            ', &
      'TiO2           ', &
      'Si2C           ', &
      'SiO2           ', &
      'FeO            ', &
      'NH2            ', &
      'NH3            ', &
      'CH2            ', &
      'CH3            ', &
      'H2S            ', &
      'VO             ', &
      'VO2            ', &
      'NaCl           ', &
      'KCl            ', &
      'e-             ', &
      'H+             ', &
      'H-             ', &
      'Na+            ', &
      'K+             ', &
      'PH2            ', &
      'P2             ', &
      'PS             ', &
      'PO             ', &
      'P4O6           ', &
      'PH             ', &
      'V              ', &
      'VO(c)          ', &
      'VO(L)          ', &
      'MgSiO3(c)      ', &
      'SiC(c)         ', &
      'Fe(c)          ', &
      'Al2O3(c)       ', &
      'Na2S(c)        ', &
      'KCl(c)         ', &
      'Fe(L)          ', &
      'SiC(L)         ', &
      'MgSiO3(L)      ', &
      'H2O(L)         ', &
      'H2O(c)         ', &
      'TiO(c)         ', &
      'TiO(L)         ', &
      'TiO2(c)        ', &
      'TiO2(L)        ', &
      'H3PO4(c)       ', &
      'H3PO4(L)       ' &
    ]

    atoms = [ &
      'H ', &
      'He', &
      'C ', &
      'N ', &
      'O ', &
      'Na', &
      'Mg', &
      'Al', &
      'Si', &
      'P ', &
      'S ', &
      'Cl', &
      'K ', &
      'Ca', &
      'Ti', &
      'V ', &
      'Fe', &
      'Ni' &
    ]

    X = [ &
      9.207539305000000e-01_dp, &
      7.836886940000000e-02_dp, &
      2.478241000000000e-04_dp, &
      6.225060569498810e-05_dp, &
      4.509658000000000e-04_dp, &
      1.600086943532050e-06_dp, &
      3.665587420553620e-05_dp, &
      2.595000000000000e-06_dp, &
      2.979500000000000e-05_dp, &
      2.366702019976680e-07_dp, &
      1.213790073460400e-05_dp, &
      2.911679584995890e-07_dp, &
      9.866056119256769e-08_dp, &
      2.014390114292550e-06_dp, &
      8.206228043663590e-08_dp, &
      7.836886940899920e-09_dp, &
      2.911679584995890e-05_dp, &
      1.528071168062810e-06_dp &
    ]

    correct_answer = [ &
      2.097572751458766e-09_dp, &
      8.531731613887943e-01_dp, &
      1.455015322134050e-01_dp, &
      1.042780525932260e-23_dp, &
      3.324518992641168e-32_dp, &
      2.215657306846010e-24_dp, &
      1.272590434869170e-05_dp, &
      1.407156949719381e-39_dp, &
      1.501096396346585e-14_dp, &
      1.279732322882177e-15_dp, &
      8.000386979842814e-27_dp, &
      5.942520535535404e-07_dp, &
      2.536655340883834e-06_dp, &
      2.837053768754091e-06_dp, &
      6.063950848321747e-15_dp, &
      7.669076761338264e-08_dp, &
      3.779282843519088e-31_dp, &
      1.230756622872010e-05_dp, &
      4.517417332085708e-15_dp, &
      1.018993933154473e-10_dp, &
      5.606388106792407e-05_dp, &
      4.395887817992532e-27_dp, &
      2.665884349002063e-27_dp, &
      6.406957289164702e-25_dp, &
      2.916480301529460e-28_dp, &
      6.482812282847260e-04_dp, &
      1.933633849489186e-38_dp, &
      2.835207940323098e-28_dp, &
      1.228293692087729e-22_dp, &
      4.463306559593289e-16_dp, &
      0.000000000000000e+00_dp, &
      1.729011724690401e-20_dp, &
      2.376605339725764e-37_dp, &
      3.817089469303068e-20_dp, &
      1.580758530397593e-19_dp, &
      8.630873323290326e-42_dp, &
      2.933670533266407e-17_dp, &
      3.477277095869979e-14_dp, &
      1.672217504051772e-27_dp, &
      2.901964760389726e-10_dp, &
      4.009508311969189e-14_dp, &
      4.477949124882182e-04_dp, &
      5.802425487208889e-25_dp, &
      1.507070706230124e-18_dp, &
      3.117988834775885e-33_dp, &
      3.145713280112587e-06_dp, &
      3.640291080449329e-11_dp, &
      1.212369761144567e-08_dp, &
      4.330783027573506e-07_dp, &
      1.340999161150335e-08_dp, &
      1.637887633641227e-23_dp, &
      0.000000000000000e+00_dp, &
      3.272522560297460e-33_dp, &
      1.822570028256660e-21_dp, &
      1.163861617512337e-14_dp, &
      3.447922902061464e-06_dp, &
      1.750125291174481e-20_dp, &
      1.952232981780049e-11_dp, &
      2.253541113131996e-05_dp, &
      3.869277213853278e-21_dp, &
      3.627438028270097e-18_dp, &
      4.341046527015537e-07_dp, &
      1.064847984094941e-07_dp, &
      1.038902490444878e-14_dp, &
      0.000000000000000e+00_dp, &
      3.271205495198492e-21_dp, &
      3.242755947534449e-17_dp, &
      1.035660061617894e-14_dp, &
      2.700056103682119e-09_dp, &
      1.809450315782121e-09_dp, &
      6.017508591220193e-12_dp, &
      1.809876338439503e-12_dp, &
      4.411013966372366e-34_dp, &
      2.509637985492699e-12_dp, &
      3.301416573163102e-24_dp, &
      1.455015321938003e-08_dp, &
      0.000000000000000e+00_dp, &
      5.531811528594553e-05_dp, &
      0.000000000000000e+00_dp, &
      5.405894509609382e-05_dp, &
      2.408969779610401e-06_dp, &
      0.000000000000000e+00_dp, &
      0.000000000000000e+00_dp, &
      0.000000000000000e+00_dp, &
      0.000000000000000e+00_dp, &
      0.000000000000000e+00_dp, &
      0.000000000000000e+00_dp, &
      0.000000000000000e+00_dp, &
      0.000000000000000e+00_dp, &
      0.000000000000000e+00_dp, &
      1.523588081832998e-07_dp, &
      0.000000000000000e+00_dp, &
      0.000000000000000e+00_dp, &
      0.000000000000000e+00_dp &
    ]

    cea = ChemEquiAnalysis(atoms, reactants, '../test/thermo_easy_chem_simp_own.inp', err)
    if (allocated(err)) then
      print*,err
      stop 1
    endif

    call cea%solve(1.0_dp, 1000.0_dp, X, err)
    if (allocated(err)) then
      print*,err
      stop 1
    endif

    if (.not.all(is_close(cea%molfracs_reactants,correct_answer))) then
      print*,'ChemEquiAnalysis failed to compute the right equilibrium.'
      stop 1
    endif

  end subroutine

  !> coppied from fortran stdlib v0.2.0
  elemental function is_close(a, b, tol, abs_tol, equal_nan) result(close)
    use, intrinsic :: ieee_arithmetic, only: ieee_is_nan
    real(dp), intent(in) :: a, b
    real(dp), intent(in), optional :: tol, abs_tol
    logical, intent(in), optional :: equal_nan
    logical :: close

    real(dp) :: rel_tol_, abs_tol_
    logical :: equal_nan_

    if (present(tol)) then
      rel_tol_ = tol
    else
      rel_tol_ = 1.0e-5_dp
    endif

    if (present(abs_tol)) then
      abs_tol_ = abs_tol
    else
      abs_tol_ = 0.0_dp
    endif

    if (present(equal_nan)) then
      equal_nan_ = equal_nan
    else
      equal_nan_ = .false.
    endif

    if (ieee_is_nan(a) .or. ieee_is_nan(b)) then
        close = merge(.true., .false., equal_nan_ .and. ieee_is_nan(a) .and. ieee_is_nan(b))
    else
        close = abs(a - b) <= max(abs(rel_tol_*max(abs(a), abs(b))), abs(abs_tol_))
    end if     

  end function


end program