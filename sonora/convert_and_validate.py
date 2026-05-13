import numpy as np
import pathlib
import convert_thermo_inp_to_yaml
import easychem.easychem as ec
from photochem.equilibrate import ChemEquiAnalysis

def convert_sonora():
    convert_thermo_inp_to_yaml.convert_file(
        input_path=pathlib.Path('./thermo-sonora-component.inp'),
        output_path=pathlib.Path('./thermo-sonora-component.yaml'),
        normalize_species_case=True,
        dedupe_names=True,
    )

# This cell inputs values from published solar abundance databases (or custom abundances)
# These include (in listed order) Lodders 2010, Lodders 2020, Asplund+ 2009, Asplund+ 2020
ABUNDS_DATA = ([
    #atom,Lo10,Lo20,Asp09,Asp20,custom,
    ('H',9.106254E-01,9.082387E-01,9.206789E-01,9.232609E-01,9.232609E-01),
    ('He',8.824980E-02,9.046346E-02,7.836248E-02,7.573985E-02,7.573985E-02),
    ('Li',1.954856E-09,2.050745E-09,1.033019E-11,8.420239E-12,8.420239E-12),
    #('Be',2.151748E-11,2.295826E-11,2.208555E-11,2.214749E-11,2.214749E-11),
    #('B',6.609945E-10,6.487419E-10,4.614325E-10,4.627266E-10,4.627266E-10),
    ('C',2.527952E-04,3.286959E-04,2.478039E-04,2.662713E-04,2.662713E-04),
    ('N',7.453768E-05,7.893027E-05,6.224553E-05,6.242010E-05,6.242010E-05),
    ('O',5.520007E-04,5.982842E-04,4.509290E-04,4.521936E-04,4.521936E-04),
    ('F',2.826806E-08,4.577235E-08,3.342783E-08,2.319126E-08,2.319126E-08),
    #('Ne',1.156740E-04,1.575001E-04,7.836248E-05,1.060045E-04,1.060045E-04),
    ('Na',2.028691E-06,2.083182E-06,1.599957E-06,1.532232E-06,1.532232E-06),
    ('Mg',3.621406E-05,3.712245E-05,3.665289E-05,3.275853E-05,3.275853E-05),
    #('Al',2.974475E-06,2.948892E-06,2.594826E-06,2.484989E-06,2.484989E-06),
    ('Si',3.515928E-05,3.604122E-05,2.979258E-05,2.987614E-05,2.987614E-05),
    ('P',2.918220E-07,2.977005E-07,2.366509E-07,2.373146E-07,2.373146E-07),
    ('S',1.480206E-05,1.575001E-05,1.213691E-05,1.217095E-05,1.217095E-05),
    ('Cl',1.817735E-07,1.906580E-07,2.911442E-07,1.885057E-07,1.885057E-07),
    #('Ar',3.259265E-06,3.521227E-06,2.312641E-06,2.214749E-06,2.214749E-06),
    ('K',1.321989E-07,1.301448E-07,9.865252E-08,1.084737E-07,1.084737E-07),
    #('Ca',2.123621E-06,2.062783E-06,2.014226E-06,1.842148E-06,1.842148E-06),
    #('Sc',1.209479E-09,1.214589E-09,1.300493E-09,1.274455E-09,1.274455E-09),
    ('Ti',8.684343E-08,8.862536E-08,8.205559E-08,8.616372E-08,8.616372E-08),
    ('V',1.005555E-08,9.911335E-09,7.836248E-09,7.333722E-09,7.333722E-09),
    ('Cr',4.605866E-07,4.732212E-07,4.018909E-07,3.848792E-07,3.848792E-07),
    #('Mn',3.241686E-07,3.276147E-07,2.478039E-07,2.428424E-07,2.428424E-07),
    ('Fe',2.981507E-05,3.142794E-05,2.911442E-05,2.662713E-05,2.662713E-05),
    #('Co',8.262431E-08,8.145315E-08,8.997217E-08,8.041266E-08,8.041266E-08),
    #('Ni',1.722805E-06,1.754126E-06,1.527947E-06,1.463270E-06,1.463270E-06),
    #('Cu',1.902117E-08,1.928205E-08,1.425963E-08,1.397412E-08,1.397412E-08),
    #('Zn',4.570707E-08,4.541193E-08,3.342783E-08,3.352158E-08,3.352158E-08),
    #('Ga',1.286830E-09,1.304692E-09,1.009504E-09,9.667728E-10,9.667728E-10),
    #('Ge',4.043317E-09,4.324946E-09,4.112522E-09,3.848792E-09,3.848792E-09),
    #('As',2.144716E-10,2.187702E-10,1.836996E-10,1.842148E-10,1.842148E-10),
    #('Se',2.373252E-09,2.436386E-09,2.014226E-09,2.019875E-09,2.019875E-09),
    #('Br',3.762043E-10,4.433070E-10,9.206789E-53,9.232609E-53,9.232609E-53),
    #('Kr',1.961888E-09,1.848914E-09,9.206789E-53,9.232609E-53,9.232609E-53),
    ('Rb',2.542016E-10,2.584155E-10,3.048654E-10,1.928965E-10,1.928965E-10),
    #('Sr',8.192113E-10,8.397604E-10,6.825087E-10,6.242010E-10,6.242010E-10),
    #('Y',1.627875E-10,1.567793E-10,1.493166E-10,1.497354E-10,1.497354E-10),
    #('Zr',3.797202E-10,3.928493E-10,3.500323E-10,3.591902E-10,3.591902E-10),
    #('Nb',2.742424E-11,2.811215E-11,2.655267E-11,2.724736E-11,2.724736E-11),
    #('Mo',8.965617E-11,9.370717E-11,6.984064E-11,7.003650E-11,7.003650E-11),
    #('Ru',6.258352E-11,6.523460E-11,5.177358E-11,5.191877E-11,5.191877E-11),
    #('Rh',1.300893E-11,1.218193E-11,7.483559E-12,5.563197E-12,5.563197E-12),
    #('Pd',4.781662E-11,4.973688E-11,3.420646E-11,3.430239E-11,3.430239E-11),
    #('Ag',1.719289E-11,1.791249E-11,8.018778E-12,8.420239E-12,8.420239E-12),
    #('Cd',5.520007E-11,5.694512E-11,4.721806E-11,4.735048E-11,4.735048E-11),
    #('In',6.258352E-12,6.451378E-12,5.809091E-12,5.825382E-12,5.825382E-12),
    #('Sn',1.265734E-10,1.293880E-10,1.009504E-10,9.667728E-11,9.667728E-11),
    #('Sb',1.100486E-11,1.293880E-11,9.421242E-12,9.447664E-12,9.447664E-12),
    #('Te',1.648970E-10,1.701145E-10,1.393504E-10,1.397412E-10,1.397412E-10),
    #('I',3.867521E-11,5.730554E-11,9.206789E-53,9.232609E-53,9.232609E-53),
    #('Xe',1.919697E-10,1.982267E-10,9.206789E-53,9.232609E-53,9.232609E-53),
    ('Cs',1.304409E-11,1.326317E-11,1.106899E-11,1.110004E-11,1.110004E-11),
    #('Ba',1.571620E-10,1.639875E-10,1.393504E-10,1.719192E-10,1.719192E-10),
    #('La',1.606779E-11,1.654292E-11,1.159066E-11,1.189390E-11,1.189390E-11),
    #('Ce',4.148795E-11,4.180781E-11,3.500323E-11,3.510140E-11,3.510140E-11),
    #('Pr',6.047397E-12,6.307213E-12,4.831791E-12,5.191877E-12,5.191877E-12),
    #('Nd',3.009635E-11,3.113961E-11,2.421632E-11,2.428424E-11,2.428424E-11),
    #('Sm',9.387528E-12,9.839253E-12,8.396691E-12,8.228571E-12,8.228571E-12),
    #('Eu',3.515928E-12,3.604122E-12,3.048654E-12,3.057204E-12,3.057204E-12),
    #('Gd',1.265734E-11,1.247026E-11,1.081703E-11,1.110004E-11,1.110004E-11),
    #('Tb',2.109557E-12,2.252576E-12,1.836996E-12,1.885057E-12,1.885057E-12),
    #('Dy',1.420435E-11,1.466878E-11,1.159066E-11,1.162317E-11,1.162317E-11),
    #('Ho',3.164335E-12,3.211273E-12,2.780406E-12,2.788203E-12,2.788203E-12),
    #('Er',9.211732E-12,9.226552E-12,7.657873E-12,7.858224E-12,7.858224E-12),
    #('Tm',1.406371E-12,1.452461E-12,1.159066E-12,1.189390E-12,1.189390E-12),
    #('Yb',9.000776E-12,9.082387E-12,6.369542E-12,6.536186E-12,6.536186E-12),
    #('Lu',1.336053E-12,1.373170E-12,1.159066E-12,1.162317E-12,1.162317E-12),
    #('Hf',5.484848E-12,5.586389E-12,6.517907E-12,6.536186E-12,6.536186E-12),
    #('Ta',7.383449E-13,7.748862E-13,6.984064E-13,6.536186E-13,6.536186E-13),
    #('W',4.816822E-12,5.189935E-12,6.517907E-12,5.692780E-12,5.692780E-12),
    #('Re',2.042754E-12,1.971455E-12,1.675360E-12,1.680059E-12,1.680059E-12),
    #('Os',2.383799E-11,2.349887E-11,2.312641E-11,2.066924E-11,2.066924E-11),
    #('Ir',2.362704E-11,2.281409E-11,2.208555E-11,1.928965E-11,1.928965E-11),
    #('Pt',4.465229E-11,4.469111E-11,3.838028E-11,3.761183E-11,3.761183E-11),
    #('Au',6.856060E-12,7.028038E-12,7.657873E-12,7.504546E-12,7.504546E-12),
    #('Hg',1.610295E-11,1.355150E-11,1.361784E-11,1.365603E-11,1.365603E-11),
    #('Tl',6.398989E-12,6.451378E-12,7.313212E-12,7.679349E-12,7.679349E-12),
    #('Pb',1.163772E-10,1.192964E-10,5.177358E-11,8.228571E-11,8.228571E-11),
    #('Bi',4.851981E-12,5.081812E-12,4.112522E-12,4.124055E-12,4.124055E-12),
    #('Th',1.547008E-12,1.517335E-12,9.640691E-13,9.892919E-13,9.892919E-13),
    #('U',6.680264E-14,8.610247E-13,2.655267E-13,2.662713E-13,2.662713E-13)
])

REACTANTS = np.array([

# Initial 31 species (up to COS) from *reported output* of original Sonora grids (Summer 2015).
# Li-bearing species and OH, C-gr added with updates following Gharib-Nezhad et al (2021).
# Additional atomics and ions added to output list (summer 2024).
# First 50 species (up to O+) are currently included in output of Sonora chemistry grids.
#
#
# NOTE: For stability and flexibility, this version uses component metal oxides as stoichiometric proxies
# for oxide and silicate condensates in substellar atmopsheres. More precise condensation curves
# may be found by replacing these oxides with the expected species in the condensate sequence.


    'e-','H2', 'H', 'H+', 'H-','H2-', 'H2+','H3+',
    'He',
    'H2O',
    'CH4','CO',
    'NH3','N2',
    'PH3','H2S',
    'TiO','VO','Fe','FeH','CrH',
    'Na','K','Rb','atCs',
    'CO2','HCN','C2H2','C2H4','C2H6','COS',
    'SiO','MgH',
    'Li','LiOH','LiH','LiCl','Li+','LiF',
    'OH','C-gr',
    'Mg','Mg+','Si','Fe+','Ti','Ti+','C','O','C+','O+',
    # Additional species included in calculation
    # MWE list that still approximates Sonora chemistry grids
    'He+',
    'C2','CH','CN',
    'CS','C2H','CH2','CH3','C3H8','HCHO','CH2OH','CH3OH','CH3O',
    'N','NH','NH2','NO','N2H2','N2H4',
    'O2','H2O2',
    'P','PH2', 'P2','PO','PH','P4O6(Gurvich)',
    #'P4O6(Gurvich)','P4O6(JANAF)','HPO2','H3PO4','PN','PS',
    'S','SH','SN',
    'SO', 'S2','SO2','S-','SH-',
    'Cr','Cr+','CrO','CrO2',
    'FeO','FeOH','FeS','Fe(OH)2','FeCl',
    'MgO','MgOH','MgS','Mg(OH)2',
    'Si+','SiS','SiH','SiO2','SiH2','SiH3','SiH4',
    'Na+', 'NaCl','NaOH','NaH',
    'K+','KCl','KH','KOH',
    'V', 'V+','VO2',
    'TiO2',
    'Cl-','Cl','HCl',#'Cl2',
    'RbCl','Rb+','RbH','RbO','RbOH','RbF',
    'CsCl','Cs+','CsH',
    'F','F-','HF','NaF',#'F2',
    # Al- and Ca-bearing species (optional)
    #'Al','AlH','AlO','AlOH','Al2O','AlCl','AlCl2','AlCl3',
    #'Ca','Ca+','CaH','CaO','CaOH',
# Condensates included in the calculation
    'NH4H2PO4(c)',
    'VO(c)','VO(L)',     #proxy
    'TiO2(c)','TiO2(L)', #proxy
    'MgO(c)','MgO(L)',   #proxy
    'SiO2(c)','SiO2(L)', #proxy
    #'MgSiO3(c)','Mg2SiO4(c)',
    'Cr(c)','Cr(L)',
    'Fe(c)','Fe(L)',
    'H2O(L)','H2O(c)',
    'Na2S(c)',
    'KCl(c)',
    'RbCl(c)',
    'CsCl(c)',
    'Li2S(c)',
    'LiF(cr)'
    #'Gehlenite(c)','Grossite(c)',
    # 'Al2O3(c)','CaO(c)'
    # Additional species as needed:
])

def process_abunds_data(data):
    return [list(col) if i == 0 else list(map(float, col))
            for i, col in enumerate(zip(*data))]

def set_easychem_atoms(exo, feh, co_factor):

    # feh = 0.00  # metallicity value in dex (i.e. +1.0 = 10x solar)
    # co_factor = 1.0 # C/O ratio **relative to solar** (i.e. 1 = solar ratio)

    # note - the filename output will include the *actual* C/O ratio of the run
    # for example co_factor = 1.0 will give co0.55 in filename for Lo20 abunds

    #******************************************************************

    natoms = len(exo.atomAbunds)

    # Metallicity adjustment for elements heavier than helium
    for i in range(2, natoms):
        exo.atomAbunds[i] *= 10.**feh

    # Calculate the new C/O ratio
    # Note: this keeps C + O = constant to keep total metallicity constant
    co_solar = exo.atomAbunds[3]/exo.atomAbunds[5]
    co_sum = exo.atomAbunds[3] + exo.atomAbunds[5]
    exo.atomAbunds[5] = co_sum/(1+co_factor*co_solar)
    exo.atomAbunds[3] = co_solar*co_factor*exo.atomAbunds[5]

def get_gridvals():
    "Values to test"
    P = np.array([1e3, 1.0, 1e-6])*1e6
    T = np.array([100, 1000, ])
    feh = np.array([])
    co_factor = np.array([])
    gridvals = (P, T, feh, co_factor)
    return gridvals

def benchmark():

    # Set up easychem
    exo = ec.ExoAtmos()
    exo._thermofpath = 'thermo-sonora-component.inp'

    # Set up equilibrate
    cea = ChemEquiAnalysis('thermo-sonora-component.yaml', species=list(REACTANTS))
    cea.mass_tol = 1e-2 # to match easyCHEM internals

    # Values to test
    test_values = [
        (1e9, 1000, 0.0, 1.0),
        (1e6, 1000, 0.0, 1.0),
        (1.0, 1000, 0.0, 1.0),
        (1e6, 200, 0.0, 1.0),
        (1.0, 200, 0.0, 1.0),
        (1e9, 1000, 3.0, 1.0),
        (1e6, 1000, 3.0, 1.0),
        (1.0, 1000, 3.0, 1.0),
        (1e6, 200, 3.0, 1.0),
        (1.0, 200, 3.0, 1.0),
        (1e9, 1000, 0.0, 5.0),
        (1e6, 1000, 0.0, 5.0),
        (1.0, 1000, 0.0, 5.0),
        (1e6, 200, 0.0, 5.0),
        (1.0, 200, 0.0, 5.0),
        (1e9, 1000, 3.0, 5.0),
        (1e6, 1000, 3.0, 5.0),
        (1.0, 1000, 3.0, 5.0),
        (1e6, 200, 3.0, 5.0),
        (1.0, 200, 3.0, 5.0),
    ] 

    for val in test_values:
        P, T, feh, co_factor = val
        test(exo, cea, P=P, T=T, feh=feh, co_factor=co_factor)

def test(exo, cea, P, T, feh, co_factor):

    # Set easychem atomic composition
    selected_atoms, Lo10_abunds, Lo20_abunds, Asp09_abunds, Asp20_abunds, custom_abunds = process_abunds_data(ABUNDS_DATA)
    exo.atoms = selected_atoms
    natoms = len(exo.atoms)
    exo.atomAbunds = np.array(Lo20_abunds)
    set_easychem_atoms(exo, feh=feh, co_factor=co_factor)
    exo.reactants = REACTANTS.copy()

    # Set equilibrate atomic composition to match
    molfracs_atoms_sun = cea.molfracs_atoms_sun*0.0
    for i,atom in enumerate(exo.atoms):
        ind = cea.atoms_names.index(atom)
        molfracs_atoms_sun[ind] = exo.atomAbunds[i]
    cea.molfracs_atoms_sun = molfracs_atoms_sun # Set solar composition

    # Solve
    P_bar = P/1e6
    assert cea.solve_metallicity(P, T, metallicity=1, CtoO=1)
    exo.solve(P_bar, T)

    # Test for consistenency
    assert cea.species_names == exo.reactants
    assert np.allclose(cea.molfracs_species, exo.reacMols, atol=1e-10)

def main():
    convert_sonora()
    benchmark()

if __name__ == '__main__':
    main()