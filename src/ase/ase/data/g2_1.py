"""
The following contains a database of small molecules

Data for the G2/97 database are from
Raghavachari, Redfern, and Pople, J. Chem. Phys. Vol. 106, 1063 (1997).
See http://www.cse.anl.gov/Catalysis_and_Energy_Conversion/Computational_Thermochemistry.shtml for the original files.

All numbers are experimental values, except for coordinates, which are
MP2(full)/6-31G(d) optimized geometries (from http://www.cse.anl.gov/OldCHMwebsiteContent/compmat/G2-97.htm)

Atomic species:
ref: Curtiss et al. JCP 106, 1063 (1997).
'Enthalpy' is the experimental enthalpies of formation at 0K
'thermal correction' is the thermal corrections H(298)-H(0)

Molecular species:
ref: Staroverov et al. JCP 119, 12129 (2003)
'Enthalpy' is the experimental enthalpies of formation at 298K
'ZPE' is the zero-point energies
'thermal correction' is the thermal enthalpy corrections H(298K) - H_exp(0K)
ZPE and thermal corrections are estimated from B3LYP geometries and vibrations.

For details about G2-1 and G2-2 sets see doi:10.1063/1.477422.

Experimental ionization potentials are from http://srdata.nist.gov/cccbdb/
Information presented on these pages is considered public information
and may be distributed or copied http://www.nist.gov/public_affairs/disclaimer.cfm

"""

from ase.atoms import string2symbols

atom_names = ['H','Li','Be','C','N','O','F','Na','Si','P','S','Cl']

molecule_names = ['LiH','BeH','CH','CH2_s3B1d','CH2_s1A1d','CH3','CH4','NH','NH2','NH3','OH','H2O','HF','SiH2_s1A1d','SiH2_s3B1d','SiH3','SiH4','PH2','PH3','SH2','HCl','Li2','LiF','C2H2','C2H4','C2H6','CN','HCN','CO','HCO','H2CO','CH3OH','N2','N2H4','NO','O2','H2O2','F2','CO2','Na2','Si2','P2','S2','Cl2','NaCl','SiO','CS','SO','ClO','ClF','Si2H6','CH3Cl','CH3SH','HOCl','SO2']

data = {
'H': {
    'name': 'Hydrogen',
    'database': 'G2-1',
    'symbols': 'H',
    'magmoms': [1.],
    'charges': None,
    'enthalpy': 51.63,
    'thermal correction': 1.01,
    'ionization energy': 13.60,
    'positions': [[ 0.  ,  0.  ,  0.]],
    },
'Li': {
    'name': 'Lithium',
    'database': 'G2-1',
    'symbols': 'Li',
    'magmoms': [1.],
    'charges': None,
    'enthalpy': 37.69,
    'thermal correction': 1.10,
    'ionization energy': 5.39,
    'positions': [[ 0.  ,  0.  ,  0.]],
    },
'Be': {
    'name': 'Beryllium',
    'database': 'G2-1',
    'symbols': 'Be',
    'magmoms': None,
    'charges': None,
    'enthalpy': 76.48,
    'thermal correction': 0.46,
    'ionization energy': 9.32,
    'positions': [[ 0.  ,  0.  ,  0.]],
    },
'C': {
    'name': 'Carbon',
    'database': 'G2-1',
    'symbols': 'C',
    'magmoms': [2.],
    'charges': None,
    'enthalpy': 169.98,
    'thermal correction': 0.25,
    'ionization energy': 11.26,
    'positions': [[ 0.  ,  0.  ,  0.]],
    },
'N': {
    'name': 'Nitrogen',
    'database': 'G2-1',
    'symbols': 'N',
    'magmoms': [3.],
    'charges': None,
    'enthalpy': 112.53,
    'thermal correction': 1.04,
    'ionization energy': 14.53,
    'positions': [[ 0.  ,  0.  ,  0.]],
    },
'O': {
    'name': 'Oxygen',
    'database': 'G2-1',
    'symbols': 'O',
    'magmoms': [2.],
    'charges': None,
    'enthalpy': 58.99,
    'thermal correction': 1.04,
    'ionization energy': 13.62,
    'positions': [[ 0.  ,  0.  ,  0.]],
    },
'F': {
    'name': 'Fluorine',
    'database': 'G2-1',
    'symbols': 'F',
    'magmoms': [1.],
    'charges': None,
    'enthalpy': 18.47,
    'thermal correction': 1.05,
    'ionization energy': 17.42,
    'positions': [[ 0.  ,  0.  ,  0.]],
    },
'Na': {
    'name': 'Sodium',
    'database': 'G2-1',
    'symbols': 'Na',
    'magmoms': [1.],
    'charges': None,
    'enthalpy': 25.69,
    'thermal correction': 1.54,
    'ionization energy': 5.14,
    'positions': [[ 0.  ,  0.  ,  0.]],
    },
'Si': {
    'name': 'Silicon',
    'database': 'G2-1',
    'symbols': 'Si',
    'magmoms': [2.],
    'charges': None,
    'enthalpy': 106.60,
    'thermal correction': 0.76,
    'ionization energy': 8.15,
    'positions': [[ 0.  ,  0.  ,  0.]],
    },
'P': {
    'name': 'Phosphorus',
    'database': 'G2-1',
    'symbols': 'P',
    'magmoms': [3.],
    'charges': None,
    'enthalpy': 75.42,
    'thermal correction': 1.28,
    'ionization energy': 10.49,
    'positions': [[ 0.  ,  0.  ,  0.]],
    },
'S': {
    'name': 'Sulfur',
    'database': 'G2-1',
    'symbols': 'S',
    'magmoms': [2.],
    'charges': None,
    'enthalpy': 65.66,
    'thermal correction': 1.05,
    'ionization energy': 10.36,
    'positions': [[ 0.  ,  0.  ,  0.]],
    },
'Cl': {
    'name': 'Chlorine',
    'database': 'G2-1',
    'symbols': 'Cl',
    'magmoms': [1.],
    'charges': None,
    'enthalpy': 28.59,
    'thermal correction': 1.10,
    'ionization energy': 12.97,
    'positions': [[ 0.  ,  0.  ,  0.]],
    },
'LiH': {
    'description': "Lithium hydride (LiH), C*v symm.",
    'name': "LiH",
    'database': 'G2-1',
    'enthalpy': 33.3,
    'ZPE': 2.0149,
    'thermal correction': 2.0783,
    'ionization energy': 7.90,
    'symbols': 'LiH',
    'magmoms': None,
    'charges': None,
    'positions': [[ 0.  ,  0.  ,  0.41],
                  [ 0.  ,  0.  , -1.23]]},
'BeH': {
    'description': "Beryllium hydride (BeH), D*h symm.",
    'name': "BeH",
    'database': 'G2-1',
    'enthalpy': 81.7,
    'ZPE': 2.9073,
    'thermal correction': 2.0739,
    'ionization energy': 8.21,
    'symbols': 'BeH',
    'magmoms': [ 0.8,  0.2],
    'charges': None,
    'positions': [[ 0.      ,  0.      ,  0.269654],
                  [ 0.      ,  0.      , -1.078616]]},
'CH': {
    'description': "CH radical. Doublet, C*v symm.",
    'name': "CH",
    'database': 'G2-1',
    'enthalpy': 142.5,
    'ZPE': 3.9659,
    'thermal correction': 2.0739,
    'ionization energy': 10.64,
    'symbols': 'CH',
    'magmoms': [ 1.,  0.],
    'charges': None,
    'positions': [[ 0.      ,  0.      ,  0.160074],
                  [ 0.      ,  0.      , -0.960446]]},
'CH2_s3B1d': {
    'description': "Triplet methylene (CH2), C2v symm, 3-B1.",
    'name': "CH_2 (^3B_1)",
    'database': 'G2-1',
    'enthalpy': 93.7,
    'ZPE': 10.6953,
    'thermal correction': 2.3877,
    'ionization energy': 10.40,
    'symbols': 'CHH',
    'magmoms': [ 2.,  0.,  0.],
    'charges': None,
    'positions': [[ 0.      ,  0.      ,  0.110381],
                  [ 0.      ,  0.982622, -0.331142],
                  [ 0.      , -0.982622, -0.331142]]},
'CH2_s1A1d': {
    'description': "Singlet methylene (CH2), C2v symm, 1-A1.",
    'name': "CH_2 (^1A_1)",
    'database': 'G2-1',
    'enthalpy': 102.8,
    'ZPE': 10.2422,
    'thermal correction': 2.3745,
    'symbols': 'CHH',
    'magmoms': None,
    'charges': None,
    'positions': [[ 0.      ,  0.      ,  0.174343],
                  [ 0.      ,  0.862232, -0.523029],
                  [ 0.      , -0.862232, -0.523029]]},
'CH3': {
    'description': "Methyl radical (CH3), D3h symm.",
    'name': "CH_3",
    'database': 'G2-1',
    'enthalpy': 35.0,
    'ZPE': 18.3383,
    'thermal correction': 2.5383,
    'ionization energy': 9.84,
    'symbols': 'CHHH',
    'magmoms': [ 1.,  0.,  0.,  0.],
    'charges': None,
    'positions': [[ 0.      ,  0.      ,  0.      ],
                  [ 0.      ,  1.07841 ,  0.      ],
                  [ 0.93393 , -0.539205,  0.      ],
                  [-0.93393 , -0.539205,  0.      ]]},
'CH4': {
    'description': "Methane (CH4), Td symm.",
    'name': "CH_4",
    'database': 'G2-1',
    'enthalpy': -17.9,
    'ZPE': 27.6744,
    'thermal correction': 2.3939,
    'ionization energy': 12.64,
    'vertical ionization energy': 13.60,
    'symbols': 'CHHHH',
    'magmoms': None,
    'charges': None,
    'positions': [[ 0.      ,  0.      ,  0.      ],
                  [ 0.629118,  0.629118,  0.629118],
                  [-0.629118, -0.629118,  0.629118],
                  [ 0.629118, -0.629118, -0.629118],
                  [-0.629118,  0.629118, -0.629118]]},
'NH': {
    'description': "NH, triplet, C*v symm.",
    'name': "NH",
    'database': 'G2-1',
    'enthalpy': 85.2,
    'ZPE': 4.5739,
    'thermal correction': 2.0739,
    'ionization energy': 13.10,
    'vertical ionization energy': 13.49,
    'symbols': 'NH',
    'magmoms': [ 2.,  0.],
    'charges': None,
    'positions': [[ 0.      ,  0.      ,  0.129929],
                  [ 0.      ,  0.      , -0.909501]]},
'NH2': {
    'description': "NH2 radical, C2v symm, 2-B1.",
    'name': "NH_2",
    'database': 'G2-1',
    'enthalpy': 45.1,
    'ZPE': 11.7420,
    'thermal correction': 2.3726,
    'ionization energy': 10.78,
    'vertical ionization energy': 12.00,
    'symbols': 'NHH',
    'magmoms': [ 1.,  0.,  0.],
    'charges': None,
    'positions': [[ 0.      ,  0.      ,  0.14169 ],
                  [ 0.      ,  0.806442, -0.495913],
                  [ 0.      , -0.806442, -0.495913]]},
'NH3': {
    'description': "Ammonia (NH3), C3v symm.",
    'name': "NH_3",
    'database': 'G2-1',
    'enthalpy': -11.0,
    'ZPE': 21.2462,
    'thermal correction': 2.3896,
    'ionization energy': 10.07,
    'vertical ionization energy': 10.82,
    'symbols': 'NHHH',
    'magmoms': None,
    'charges': None,
    'positions': [[ 0.      ,  0.      ,  0.116489],
                  [ 0.      ,  0.939731, -0.271808],
                  [ 0.813831, -0.469865, -0.271808],
                  [-0.813831, -0.469865, -0.271808]]},
'OH': {
    'description': "OH radical, C*v symm.",
    'name': "OH",
    'database': 'G2-1',
    'enthalpy': 9.4,
    'ZPE': 5.2039,
    'thermal correction': 2.0739,
    'ionization energy': 13.02,
    'symbols': 'OH',
    'magmoms': [ 0.5,  0.5],
    'charges': None,
    'positions': [[ 0.      ,  0.      ,  0.108786],
                  [ 0.      ,  0.      , -0.870284]]},
'H2O': {
    'description': "Water (H2O), C2v symm.",
    'name': "H_2O",
    'database': 'G2-1',
    'enthalpy': -57.8,
    'ZPE': 13.2179,
    'thermal correction': 2.3720,
    'ionization energy': 12.62,
    'symbols': 'OHH',
    'magmoms': None,
    'charges': None,
    'positions': [[ 0.      ,  0.      ,  0.119262],
                  [ 0.      ,  0.763239, -0.477047],
                  [ 0.      , -0.763239, -0.477047]]},
'HF': {
    'description': "Hydrogen fluoride (HF), C*v symm.",
    'name': "HF",
    'database': 'G2-1',
    'enthalpy': -65.1,
    'ZPE': 5.7994,
    'thermal correction': 2.0733,
    'ionization energy': 16.03,
    'vertical ionization energy': 16.12,
    'symbols': 'FH',
    'magmoms': None,
    'charges': None,
    'positions': [[ 0.      ,  0.      ,  0.093389],
                  [ 0.      ,  0.      , -0.840502]]},
'SiH2_s1A1d': {
    'description': "Singlet silylene (SiH2), C2v symm, 1-A1.",
    'name': "SiH_2 (^1A_1)",
    'database': 'G2-1',
    'enthalpy': 65.2,
    'ZPE': 7.1875,
    'thermal correction': 2.3927,
    'ionization energy': 8.92,
    'symbols': 'SiHH',
    'magmoms': None,
    'charges': None,
    'positions': [[ 0.      ,  0.      ,  0.131272],
                  [ 0.      ,  1.096938, -0.918905],
                  [ 0.      , -1.096938, -0.918905]]},
'SiH2_s3B1d': {
    'description': "Triplet silylene (SiH2), C2v symm, 3-B1.",
    'name': "SiH_2 (^3B_1)",
    'database': 'G2-1',
    'enthalpy': 86.2,
    'ZPE': 7.4203,
    'thermal correction': 2.4078,
    'symbols': 'SiHH',
    'magmoms': [ 2.,  0.,  0.],
    'charges': None,
    'positions': [[ 0.      ,  0.      ,  0.094869],
                  [ 0.      ,  1.271862, -0.664083],
                  [ 0.      , -1.271862, -0.664083]]},
'SiH3': {
    'description': "Silyl radical (SiH3), C3v symm.",
    'name': "SiH_3",
    'database': 'G2-1',
    'enthalpy': 47.9,
    'ZPE': 13.0898,
    'thermal correction': 2.4912,
    'ionization energy': 8.14,
    'vertical ionization energy': 8.74,
    'symbols': 'SiHHH',
    'magmoms': [ 1.,  0.,  0.,  0.],
    'charges': None,
    'positions': [[ 0.      ,  0.      ,  0.079299],
                  [ 0.      ,  1.41328 , -0.370061],
                  [ 1.223937, -0.70664 , -0.370061],
                  [-1.223937, -0.70664 , -0.370061]]},
'SiH4': {
    'description': "Silane (SiH4), Td symm.",
    'name': "SiH_4",
    'database': 'G2-1',
    'enthalpy': 8.2,
    'ZPE': 19.2664,
    'thermal correction': 2.5232,
    'ionization energy': 11.00,
    'vertical ionization energy': 12.30,
    'symbols': 'SiHHHH',
    'magmoms': None,
    'charges': None,
    'positions': [[ 0.      ,  0.      ,  0.      ],
                  [ 0.856135,  0.856135,  0.856135],
                  [-0.856135, -0.856135,  0.856135],
                  [-0.856135,  0.856135, -0.856135],
                  [ 0.856135, -0.856135, -0.856135]]},
'PH2': {
    'description': "PH2 radical, C2v symm.",
    'name': "PH_2",
    'database': 'G2-1',
    'enthalpy': 33.1,
    'ZPE': 8.2725,
    'thermal correction': 2.3845,
    'ionization energy': 9.82,
    'symbols': 'PHH',
    'magmoms': [ 1.,  0.,  0.],
    'charges': None,
    'positions': [[ 0.      ,  0.      ,  0.115396],
                  [ 0.      ,  1.025642, -0.865468],
                  [ 0.      , -1.025642, -0.865468]]},
'PH3': {
    'description': "Phosphine (PH3), C3v symm.",
    'name': "PH_3",
    'database': 'G2-1',
    'enthalpy': 1.3,
    'ZPE': 14.7885,
    'thermal correction': 2.4203,
    'ionization energy': 9.87,
    'vertical ionization energy': 10.95,
    'symbols': 'PHHH',
    'magmoms': None,
    'charges': None,
    'positions': [[ 0.      ,  0.      ,  0.124619],
                  [ 0.      ,  1.200647, -0.623095],
                  [ 1.039791, -0.600323, -0.623095],
                  [-1.039791, -0.600323, -0.623095]]},
'SH2': {
    'description': "Hydrogen sulfide (H2S), C2v symm.",
    'name': "SH_2",
    'database': 'G2-1',
    'enthalpy': -4.9,
    'ZPE': 9.3129,
    'thermal correction': 2.3808,
    'ionization energy': 10.46,
    'vertical ionization energy': 10.50,
    'symbols': 'SHH',
    'magmoms': None,
    'charges': None,
    'positions': [[ 0.      ,  0.      ,  0.102135],
                  [ 0.      ,  0.974269, -0.817083],
                  [ 0.      , -0.974269, -0.817083]]},
'HCl': {
    'description': "Hydrogen chloride (HCl), C*v symm.",
    'name': "HCl",
    'database': 'G2-1',
    'enthalpy': -22.1,
    'ZPE': 4.1673,
    'thermal correction': 2.0739,
    'ionization energy': 12.74,
    'symbols': 'ClH',
    'magmoms': None,
    'charges': None,
    'positions': [[ 0.      ,  0.      ,  0.07111 ],
                  [ 0.      ,  0.      , -1.208868]]},
'Li2': {
    'description': "Dilithium (Li2), D*h symm.",
    'name': "Li_2",
    'database': 'G2-1',
    'enthalpy': 51.6,
    'ZPE': 0.4838,
    'thermal correction': 2.3086,
    'ionization energy': 5.11,
    'symbols': 'LiLi',
    'magmoms': None,
    'charges': None,
    'positions': [[ 0.     ,  0.     ,  1.38653],
                  [ 0.     ,  0.     , -1.38653]]},
'LiF': {
    'description': "Lithium Fluoride (LiF), C*v symm.",
    'name': "LiF",
    'database': 'G2-1',
    'enthalpy': -80.1,
    'ZPE': 1.4019,
    'thermal correction': 2.0990,
    'ionization energy': 11.30,
    'symbols': 'LiF',
    'magmoms': None,
    'charges': None,
    'positions': [[ 0.      ,  0.      , -1.174965],
                  [ 0.      ,  0.      ,  0.391655]]},
'C2H2': {
    'description': "Acetylene (C2H2), D*h symm.",
    'name': "C_2H_2",
    'database': 'G2-1',
    'enthalpy': 54.2,
    'ZPE': 16.6001,
    'thermal correction': 2.4228,
    'ionization energy': 11.40,
    'vertical ionization energy': 11.49,
    'symbols': 'CCHH',
    'magmoms': None,
    'charges': None,
    'positions': [[ 0.     ,  0.     ,  0.60808],
                  [ 0.     ,  0.     , -0.60808],
                  [ 0.     ,  0.     , -1.67399],
                  [ 0.     ,  0.     ,  1.67399]]},
'C2H4': {
    'description': "Ethylene (H2C=CH2), D2h symm.",
    'name': "C_2H_4",
    'database': 'G2-1',
    'enthalpy': 12.5,
    'ZPE': 31.5267,
    'thermal correction': 2.5100,
    'ionization energy': 11.40,
    'vertical ionization energy': 11.49,
    'symbols': 'CCHHHH',
    'magmoms': None,
    'charges': None,
    'positions': [[ 0.      ,  0.      ,  0.66748 ],
                  [ 0.      ,  0.      , -0.66748 ],
                  [ 0.      ,  0.922832,  1.237695],
                  [ 0.      , -0.922832,  1.237695],
                  [ 0.      ,  0.922832, -1.237695],
                  [ 0.      , -0.922832, -1.237695]]},
'C2H6': {
    'description': "Ethane (H3C-CH3), D3d symm.",
    'name': "C_2H_6",
    'database': 'G2-1',
    'enthalpy': -20.1,
    'ZPE': 46.0950,
    'thermal correction': 2.7912,
    'symbols': 'CCHHHHHH',
    'magmoms': None,
    'charges': None,
    'positions': [[ 0.      ,  0.      ,  0.762209],
                  [ 0.      ,  0.      , -0.762209],
                  [ 0.      ,  1.018957,  1.157229],
                  [-0.882443, -0.509479,  1.157229],
                  [ 0.882443, -0.509479,  1.157229],
                  [ 0.      , -1.018957, -1.157229],
                  [-0.882443,  0.509479, -1.157229],
                  [ 0.882443,  0.509479, -1.157229]]},
'CN': {
    'description': "Cyano radical (CN), C*v symm, 2-Sigma+.",
    'name': "CN",
    'database': 'G2-1',
    'enthalpy': 104.9,
    'ZPE': 3.0183,
    'thermal correction': 2.0739,
    'ionization energy': 13.60,
    'symbols': 'CN',
    'magmoms': [ 1.,  0.],
    'charges': None,
    'positions': [[ 0.      ,  0.      , -0.611046],
                  [ 0.      ,  0.      ,  0.523753]]},
'HCN': {
    'description': "Hydrogen cyanide (HCN), C*v symm.",
    'name': "HCN",
    'database': 'G2-1',
    'enthalpy': 31.5,
    'ZPE': 10.2654,
    'thermal correction': 2.1768,
    'ionization energy': 13.60,
    'vertical ionization energy': 13.61,
    'symbols': 'CNH',
    'magmoms': None,
    'charges': None,
    'positions': [[ 0.      ,  0.      , -0.511747],
                  [ 0.      ,  0.      ,  0.664461],
                  [ 0.      ,  0.      , -1.580746]]},
'CO': {
    'description': "Carbon monoxide (CO), C*v symm.",
    'name': "CO",
    'database': 'G2-1',
    'enthalpy': -26.4,
    'ZPE': 3.1062,
    'thermal correction': 2.0739,
    'ionization energy': 14.01,
    'vertical ionization energy': 14.01,
    'symbols': 'OC',
    'magmoms': None,
    'charges': None,
    'positions': [[ 0.      ,  0.      ,  0.493003],
                  [ 0.      ,  0.      , -0.657337]]},
'HCO': {
    'description': "HCO radical, Bent Cs symm.",
    'name': "HCO",
    'database': 'G2-1',
    'enthalpy': 10.0,
    'ZPE': 8.0290,
    'thermal correction': 2.3864,
    'ionization energy': 8.12,
    'vertical ionization energy': 9.31,
    'symbols': 'COH',
    'magmoms': [ 1.,  0.,  0.],
    'charges': None,
    'positions': [[ 0.06256 ,  0.593926,  0.      ],
                  [ 0.06256 , -0.596914,  0.      ],
                  [-0.875835,  1.211755,  0.      ]]},
'H2CO': {
    'description': "Formaldehyde (H2C=O), C2v symm.",
    'name': "H_2CO",
    'database': 'G2-1',
    'enthalpy': -26.0,
    'ZPE': 16.4502,
    'thermal correction': 2.3927,
    'ionization energy': 10.88,
    'vertical ionization energy': 10.88,
    'symbols': 'OCHH',
    'magmoms': None,
    'charges': None,
    'positions': [[ 0.      ,  0.      ,  0.683501],
                  [ 0.      ,  0.      , -0.536614],
                  [ 0.      ,  0.93439 , -1.124164],
                  [ 0.      , -0.93439 , -1.124164]]},
'CH3OH': {
    'description': "Methanol (CH3-OH), Cs symm.",
    'name': "H_3COH",
    'database': 'G2-1',
    'enthalpy': -48.0,
    'ZPE': 31.6635,
    'thermal correction': 2.6832,
    'ionization energy': 10.84,
    'vertical ionization energy': 10.96,
    'symbols': 'COHHHH',
    'magmoms': None,
    'charges': None,
    'positions': [[-0.047131,  0.664389,  0.      ],
                  [-0.047131, -0.758551,  0.      ],
                  [-1.092995,  0.969785,  0.      ],
                  [ 0.878534, -1.048458,  0.      ],
                  [ 0.437145,  1.080376,  0.891772],
                  [ 0.437145,  1.080376, -0.891772]]},
'N2': {
    'description': "N2 molecule, D*h symm.",
    'name': "N_2",
    'database': 'G2-1',
    'enthalpy': 0.0,
    'ZPE': 3.4243,
    'thermal correction': 2.0733,
    'ionization energy': 15.58,
    'vertical ionization energy': 15.58,
    'symbols': 'NN',
    'magmoms': None,
    'charges': None,
    'positions': [[ 0.     ,  0.     ,  0.56499],
                  [ 0.     ,  0.     , -0.56499]]},
'N2H4': {
    'description': "Hydrazine (H2N-NH2), C2 symm.",
    'name': "H_2NNH_2",
    'database': 'G2-1',
    'enthalpy': 22.8,
    'ZPE': 32.9706,
    'thermal correction': 2.6531,
    'ionization energy': 8.10,
    'vertical ionization energy': 8.98,
    'symbols': 'NNHHHH',
    'magmoms': None,
    'charges': None,
    'positions': [[ 0.      ,  0.718959, -0.077687],
                  [ 0.      , -0.718959, -0.077687],
                  [ 0.211082,  1.092752,  0.847887],
                  [-0.948214,  1.005026, -0.304078],
                  [-0.211082, -1.092752,  0.847887],
                  [ 0.948214, -1.005026, -0.304078]]},
'NO': {
    'description': "NO radical, C*v symm, 2-Pi.",
    'name': "NO",
    'database': 'G2-1',
    'enthalpy': 21.6,
    'ZPE': 2.7974,
    'thermal correction': 2.0745,
    'ionization energy': 9.26,
    'vertical ionization energy': 9.26,
    'symbols': 'NO',
    'magmoms': [ 0.6,  0.4],
    'charges': None,
    'positions': [[ 0.      ,  0.      , -0.609442],
                  [ 0.      ,  0.      ,  0.533261]]},
'O2': {
    'description': "O2 molecule, D*h symm, Triplet.",
    'name': "O_2",
    'database': 'G2-1',
    'enthalpy': 0.0,
    'ZPE': 2.3444,
    'thermal correction': 2.0752,
    'ionization energy': 12.07,
    'vertical ionization energy': 12.30,
    'symbols': 'OO',
    'magmoms': [ 1.,  1.],
    'charges': None,
    'positions': [[ 0.      ,  0.      ,  0.622978],
                  [ 0.      ,  0.      , -0.622978]]},
'H2O2': {
    'description': "Hydrogen peroxide (HO-OH), C2 symm.",
    'name': "HOOH",
    'database': 'G2-1',
    'enthalpy': -32.5,
    'ZPE': 16.4081,
    'thermal correction': 2.6230,
    'ionization energy': 10.58,
    'vertical ionization energy': 11.70,
    'symbols': 'OOHH',
    'magmoms': None,
    'charges': None,
    'positions': [[ 0.      ,  0.734058, -0.05275 ],
                  [ 0.      , -0.734058, -0.05275 ],
                  [ 0.839547,  0.880752,  0.422001],
                  [-0.839547, -0.880752,  0.422001]]},
'F2': {
    'description': "F2 molecule, D*h symm.",
    'name': "F_2",
    'database': 'G2-1',
    'enthalpy': 0.0,
    'ZPE': 1.5179,
    'thermal correction': 2.0915,
    'ionization energy': 15.70,
    'vertical ionization energy': 15.70,
    'symbols': 'FF',
    'magmoms': None,
    'charges': None,
    'positions': [[ 0.      ,  0.      ,  0.710304],
                  [ 0.      ,  0.      , -0.710304]]},
'CO2': {
    'description': "Carbon dioxide (CO2), D*h symm.",
    'name': "CO_2",
    'database': 'G2-1',
    'enthalpy': -94.1,
    'ZPE': 7.3130,
    'thermal correction': 2.2321,
    'ionization energy': 13.78,
    'vertical ionization energy': 13.78,
    'symbols': 'COO',
    'magmoms': None,
    'charges': None,
    'positions': [[ 0.      ,  0.      ,  0.      ],
                  [ 0.      ,  0.      ,  1.178658],
                  [ 0.      ,  0.      , -1.178658]]},
'Na2': {
    'description': "Disodium (Na2), D*h symm.",
    'name': "Na_2",
    'database': 'G2-1',
    'enthalpy': 34.0,
    'ZPE': 0.2246,
    'thermal correction': 2.4699,
    'ionization energy': 4.89,
    'symbols': 'NaNa',
    'magmoms': None,
    'charges': None,
    'positions': [[ 0.      ,  0.      ,  1.576262],
                  [ 0.      ,  0.      , -1.576262]]},
'Si2': {
    'description': "Si2 molecule, D*h symm, Triplet (3-Sigma-G-).",
    'name': "Si_2",
    'database': 'G2-1',
    'enthalpy': 139.9,
    'ZPE': 0.7028,
    'thermal correction': 2.2182,
    'ionization energy': 7.90,
    'symbols': 'SiSi',
    'magmoms': [ 1.,  1.],
    'charges': None,
    'positions': [[ 0.      ,  0.      ,  1.130054],
                  [ 0.      ,  0.      , -1.130054]]},
'P2': {
    'description': "P2 molecule, D*h symm.",
    'name': "P_2",
    'database': 'G2-1',
    'enthalpy': 34.3,
    'ZPE': 1.1358,
    'thermal correction': 2.1235,
    'ionization energy': 10.53,
    'vertical ionization energy': 10.62,
    'symbols': 'PP',
    'magmoms': None,
    'charges': None,
    'positions': [[ 0.      ,  0.      ,  0.966144],
                  [ 0.      ,  0.      , -0.966144]]},
'S2': {
    'description': "S2 molecule, D*h symm, triplet.",
    'name': "S_2",
    'database': 'G2-1',
    'enthalpy': 30.7,
    'ZPE': 1.0078,
    'thermal correction': 2.1436,
    'ionization energy': 9.36,
    'vertical ionization energy': 9.55,
    'symbols': 'SS',
    'magmoms': [ 1.,  1.],
    'charges': None,
    'positions': [[ 0.      ,  0.      ,  0.960113],
                  [ 0.      ,  0.      , -0.960113]]},
'Cl2': {
    'description': "Cl2 molecule, D*h symm.",
    'name': "Cl_2",
    'database': 'G2-1',
    'enthalpy': 0.0,
    'ZPE': 0.7737,
    'thermal correction': 2.1963,
    'ionization energy': 11.48,
    'vertical ionization energy': 11.49,
    'symbols': 'ClCl',
    'magmoms': None,
    'charges': None,
    'positions': [[ 0.      ,  0.      ,  1.007541],
                  [ 0.      ,  0.      , -1.007541]]},
'NaCl': {
    'description': "Sodium Chloride (NaCl), C*v symm.",
    'name': "NaCl",
    'database': 'G2-1',
    'enthalpy': -43.6,
    'ZPE': 0.5152,
    'thermal correction': 2.2935,
    'ionization energy': 9.20,
    'vertical ionization energy': 9.80,
    'symbols': 'NaCl',
    'magmoms': None,
    'charges': None,
    'positions': [[ 0.     ,  0.     , -1.45166],
                  [ 0.     ,  0.     ,  0.93931]]},
'SiO': {
    'description': "Silicon monoxide (SiO), C*v symm.",
    'name': "SiO",
    'database': 'G2-1',
    'enthalpy': -24.6,
    'ZPE': 1.7859,
    'thermal correction': 2.0821,
    'ionization energy': 11.49,
    'symbols': 'SiO',
    'magmoms': None,
    'charges': None,
    'positions': [[ 0.      ,  0.      ,  0.560846],
                  [ 0.      ,  0.      , -0.98148 ]]},
'CS': {
    'description': "Carbon monosulfide (CS), C*v symm.",
    'name': "SC",
    'database': 'G2-1',
    'enthalpy': 66.9,
    'ZPE': 1.8242,
    'thermal correction': 2.0814,
    'ionization energy': 11.33,
    'symbols': 'CS',
    'magmoms': None,
    'charges': None,
    'positions': [[ 0.      ,  0.      , -1.123382],
                  [ 0.      ,  0.      ,  0.421268]]},
'SO': {
    'description': "Sulfur monoxide (SO), C*v symm, triplet.",
    'name': "SO",
    'database': 'G2-1',
    'enthalpy': 1.2,
    'ZPE': 1.6158,
    'thermal correction': 2.0877,
    'ionization energy': 11.29,
    'symbols': 'OS',
    'magmoms': [ 1.,  1.],
    'charges': None,
    'positions': [[ 0.      ,  0.      , -1.015992],
                  [ 0.      ,  0.      ,  0.507996]]},
'ClO': {
    'description': "ClO radical, C*v symm, 2-PI.",
    'name': "ClO",
    'database': 'G2-1',
    'enthalpy': 24.2,
    'ZPE': 1.1923,
    'thermal correction': 2.1172,
    'ionization energy': 10.89,
    'vertical ionization energy': 11.01,
    'symbols': 'ClO',
    'magmoms': [ 1.,  0.],
    'charges': None,
    'positions': [[ 0.      ,  0.      ,  0.514172],
                  [ 0.      ,  0.      , -1.092615]]},
'ClF': {
    'description': "ClF molecule, C*v symm, 1-SG.",
    'name': "FCl",
    'database': 'G2-1',
    'enthalpy': -13.2,
    'ZPE': 1.1113,
    'thermal correction': 2.1273,
    'ionization energy': 12.66,
    'vertical ionization energy': 12.77,
    'symbols': 'FCl',
    'magmoms': None,
    'charges': None,
    'positions': [[ 0.      ,  0.      , -1.084794],
                  [ 0.      ,  0.      ,  0.574302]]},
'Si2H6': {
    'description': "Disilane (H3Si-SiH3), D3d symm.",
    'name': "Si_2H_6",
    'database': 'G2-1',
    'enthalpy': 19.1,
    'ZPE': 30.2265,
    'thermal correction': 3.7927,
    'ionization energy': 9.74,
    'vertical ionization energy': 10.53,
    'symbols': 'SiSiHHHHHH',
    'magmoms': None,
    'charges': None,
    'positions': [[ 0.      ,  0.      ,  1.167683],
                  [ 0.      ,  0.      , -1.167683],
                  [ 0.      ,  1.393286,  1.68602 ],
                  [-1.206621, -0.696643,  1.68602 ],
                  [ 1.206621, -0.696643,  1.68602 ],
                  [ 0.      , -1.393286, -1.68602 ],
                  [-1.206621,  0.696643, -1.68602 ],
                  [ 1.206621,  0.696643, -1.68602 ]]},
'CH3Cl': {
    'description': "Methyl chloride (CH3Cl), C3v symm.",
    'name': "CH_3Cl",
    'database': 'G2-1',
    'enthalpy': -19.6,
    'ZPE': 23.3013,
    'thermal correction': 2.4956,
    'symbols': 'CClHHH',
    'ionization energy': 11.26,
    'vertical ionization energy': 11.29,
    'magmoms': None,
    'charges': None,
    'positions': [[ 0.      ,  0.      , -1.121389],
                  [ 0.      ,  0.      ,  0.655951],
                  [ 0.      ,  1.029318, -1.47428 ],
                  [ 0.891415, -0.514659, -1.47428 ],
                  [-0.891415, -0.514659, -1.47428 ]]},
'CH3SH': {
    'description': "Methanethiol (H3C-SH), Staggered, Cs symm.",
    'name': "H_3CSH",
    'database': 'G2-1',
    'enthalpy': -5.5,
    'ZPE': 28.3973,
    'thermal correction': 2.8690,
    'ionization energy': 9.44,
    'vertical ionization energy': 9.44,
    'symbols': 'CSHHHH',
    'magmoms': None,
    'charges': None,
    'positions': [[-0.047953,  1.149519,  0.      ],
                  [-0.047953, -0.664856,  0.      ],
                  [ 1.283076, -0.823249,  0.      ],
                  [-1.092601,  1.461428,  0.      ],
                  [ 0.432249,  1.551207,  0.892259],
                  [ 0.432249,  1.551207, -0.892259]]},
'HOCl': {
    'description': "HOCl molecule, Cs symm.",
    'name': "HOCl",
    'database': 'G2-1',
    'enthalpy': -17.8,
    'ZPE': 8.1539,
    'thermal correction': 2.4416,
    'ionization energy': 11.12,
    'symbols': 'OHCl',
    'magmoms': None,
    'charges': None,
    'positions': [[ 0.036702,  1.113517,  0.      ],
                  [-0.917548,  1.328879,  0.      ],
                  [ 0.036702, -0.602177,  0.      ]]},
'SO2': {
    'description': "Sulfur dioxide (SO2), C2v symm.",
    'name': "SO_2",
    'database': 'G2-1',
    'enthalpy': -71.0,
    'ZPE': 4.3242,
    'thermal correction': 2.5245,
    'ionization energy': 12.35,
    'vertical ionization energy': 12.50,
    'symbols': 'SOO',
    'magmoms': None,
    'charges': None,
    'positions': [[ 0.      ,  0.      ,  0.370268],
                  [ 0.      ,  1.277617, -0.370268],
                  [ 0.      , -1.277617, -0.370268]]},
}

def get_ionization_energy(name, vertical=True):
    """Return the experimental ionization energy from the database.

    If vertical is True, the vertical ionization energy is returned if
    available.
    """
    if name not in data:
        raise KeyError('System %s not in database.' % name)
    elif 'ionization energy' not in data[name]:
        raise KeyError('No data on ionization energy for system %s.' % name)
    else:
        if vertical and 'vertical ionization energy' in data[name]:
            return data[name]['vertical ionization energy']
        else:
            return data[name]['ionization energy']


def get_atomization_energy(name):
    """Determine extrapolated experimental atomization energy from the database.

    The atomization energy is extrapolated from experimental heats of
    formation at room temperature, using calculated zero-point energies
    and thermal corrections.

    The atomization energy is returned in kcal/mol = 43.36 meV:

    >>> from ase.units import *; print kcal / mol
    0.0433641146392

    """
    d = data[name]
    e = d['enthalpy']
    z = d['ZPE']
    dh = d['thermal correction']
    ae = -e + z + dh
    for a in string2symbols(d['symbols']):
        h = data[a]['enthalpy']
        dh = data[a]['thermal correction']
        ae += h - dh
    return ae