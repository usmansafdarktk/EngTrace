import math

# Speed of light in vacuum (m/s)
# @kind: defined
# @units: m/s
# @domain: none (exact at every condition: the 2019 SI fixes c by definition)
# [ON-DISK] codata_2022/allascii.txt @ quantity="speed of light in vacuum" precision=exact
#   CODATA 2022 prints 299 792 458 "(exact)": since the 2019 SI the metre is
#   defined through c, so no measurement can move it. The warrant is still an
#   artefact rather than a bare [BY-DEFINITION], because a mistyped digit is a
#   real failure that an artefact can catch and a definition cannot (D-071).
C0 = 299792458 

# Dictionary of common media and their approximate phase velocities for EM waves
# @kind: property
# @units: m/s
# @domain: wavelength=0.589 um, band=optical
#   Every row but the tissue row is an OPTICAL index, read at the sodium-line
#   wavelength the group comments assume; each tag states its dataset's conditions.
#   The one consumer, template_wave_parameters_basic, draws 50-500 MHz or a
#   0.1-2.0 m wavelength: radio, not optical. A correct index used outside its
#   band is the D-032 class - registered and referred to the repo owner (C3.5).
MEDIA_VELOCITIES = {
    # Gases (at 0°C and 1 atm, for visible light ~589 nm)
    # [BY-DEFINITION] n = 1 in vacuum, so the phase velocity is C0 itself
    "Vacuum": C0,                     
    # [ON-DISK] refractiveindex_info/refractiveindex.info-database-main.zip @ member="database/data/other/mixed gases/air/nk/Borzsonyi.yml" wavelength=589nm via="C0/x" tol=0.0005%
    #   n = 1.0002882 at 273 K, 100000 Pa - of the datasets covering 589 nm, the one nearest the
    #   group comment's 0 degC and 1 atm (its pressure is 100000 Pa, not 101325).
    #   The row's 1.000293 is +4.8e-06 on n but +1.7% on the
    #   refractivity n-1. Ciddor at 288.15 K, 101325 Pa gives 1.0002772.
    "Air (at sea level)": C0 / 1.000293, 
    # [KNOWN-DEFECTIVE] n = 1.000036 against 1.0000349 (Ermolov) and 1.0000349 (Mansfield), both at
    #   273.15 K and 101325 Pa - the group comment's own conditions - and the row is not
    #   a rounding of either at its own precision: +1.1e-06 on n, below every display in the
    #   corpus, but +3.2% on the refractivity n-1. Not corrected: a P6 event,
    #   referred (C3.5). First tagged tol=0.0002% (08d3c49); retagged under the C3 rule
    #   that tol= is only for an artefact whose conditions differ from the table's.
    "Helium": C0 / 1.000036,          
    # [ON-DISK] refractiveindex_info/refractiveindex.info-database-main.zip @ member="database/data/main/CO2/nk/Bideau-Mehu.yml" wavelength=589nm via="C0/x" precision=6sf
    #   n = 1.0004489 at 273.15 K, 101325 Pa; Old (273.15 K, 101325 Pa) gives 1.0004491, also 1.00045.
    "Carbon Dioxide": C0 / 1.00045,   
    
    # Liquids (for visible light ~589 nm)
    # [ON-DISK] refractiveindex_info/refractiveindex.info-database-main.zip @ member="database/data/main/H2O/nk/Daimon-20.0C.yml" wavelength=589nm via="C0/x" precision=4sf
    #   n = 1.3333584 at 293.15 K, the row's own 20 degC.
    "Water (distilled, 20°C)": C0 / 1.333,  
    # [ON-DISK] refractiveindex_info/refractiveindex.info-database-main.zip @ member="database/data/organic/C2H6O - ethanol/nk/Kedenburg.yml" wavelength=589nm via="C0/x" precision=3sf
    #   n = 1.3615137 at 293 K; Chang (295 K) 1.3613567 and Kozma (293 K)
    #   1.3607418 also round to 1.36.
    "Ethanol": C0 / 1.36,             
    # [ON-DISK] refractiveindex_info/refractiveindex.info-database-main.zip @ member="database/data/organic/C3H8O3 - glycerol/nk/Gupta.yml" wavelength=589nm via="C0/x" tol=0.12%
    #   n = 1.4712766 (Gupta; no temperature stated in the file). The row's 1.473 is +0.117%.
    #   Birkhoff's tabulated 1.4714963 interpolates between rows 1.48 and 1.47,
    #   so its 3-s.f. rounding depends on the interpolation, and it is not cited.
    #   Rheims' formula starts at 0.5893 um and is refused at 0.589.
    "Glycerine": C0 / 1.473,          
    # [ON-DISK] refractiveindex_info/refractiveindex.info-database-main.zip @ member="database/data/organic/C6H6 - benzene/nk/Chang.yml" wavelength=589nm via="C0/x" tol=0.1%
    #   n = 1.4995125 at 295 K (Chang). The row's 1.501 is +0.099%;
    #   Moutzouris at 300 K gives 1.4956337 (+0.36%).
    "Benzene": C0 / 1.501,            
    # [ON-DISK] refractiveindex_info/refractiveindex.info-database-main.zip @ member="database/data/main/CS2/nk/Chemnitz.yml" wavelength=589nm via="C0/x" precision=4sf
    #   n = 1.6281390 at 293 K (Chemnitz); Kedenburg (no temperature stated in the file) 1.6275063 also
    #   rounds to 1.628, Chang at 295 K, 1.6233939, does not.
    "Carbon Disulfide": C0 / 1.628,   # notable for high dispersion
    
    # Solids (for visible light ~589 nm)
    # [ON-DISK] refractiveindex_info/refractiveindex.info-database-main.zip @ member="database/data/main/H2O/nk/Warren-2008.yml" wavelength=589nm via="C0/x" precision=3sf
    #   n = 1.30973 interpolated at 266.15 K (Warren 2008); both bracketing rows,
    #   1.31 at 0.58 um and 1.3097 at 0.59 um, are 1.31 at 3 s.f.
    "Ice": C0 / 1.31,                 
    # [UNVERIFIED] no PTFE (polytetrafluoroethylene) dataset in the on-disk
    #   refractiveindex.info archive (member names searched: PTFE, teflon,
    #   tetrafluoroethylene, (C2F4)n). The nearest, (C2ClF3)n polychlorotrifluoroethylene,
    #   is a different polymer and has no row bracketing 589 nm. Residual register, C3.5.
    "Teflon (PTFE)": C0 / 1.35,       
    # [ON-DISK] refractiveindex_info/refractiveindex.info-database-main.zip @ member="database/data/main/SiO2/nk/Malitson.yml" wavelength=589nm via="C0/x" precision=4sf
    #   n = 1.4584132 at 293 K (Malitson); Arosa (no temperature stated in the file) 1.4584726.
    "Fused Silica (Glass)": C0 / 1.458, 
    # [ON-DISK] refractiveindex_info/refractiveindex.info-database-main.zip @ member="database/data/specs/schott/optical/N-BK7.yml" wavelength=589nm via="C0/x" precision=3sf
    #   n = 1.5167401 at 293 K for Schott N-BK7. The row names a FAMILY ("typical"),
    #   not a material; N-BK7 is the dataset this citation resolves against, chosen by
    #   the implementer. Whether the row should name one material is referred to the
    #   repo owner (C3.5).
    "Crown Glass (typical)": C0 / 1.52, 
    # [UNVERIFIED] the archive's two polyethylene datasets, (C2H4)n David and Smith,
    #   have no rows bracketing 589 nm; residual register, C3.5.
    "Polyethylene": C0 / 1.54,        
    # [ON-DISK] refractiveindex_info/refractiveindex.info-database-main.zip @ member="database/data/organic/(C8H8)n - polystyrene/nk/Sultanova.yml" wavelength=589nm via="C0/x" precision=3sf
    #   n = 1.5914841 at 293 K (Sultanova); Zhang (no temperature stated in the file) 1.5890780 also 1.59.
    "Polystyrene": C0 / 1.59,         
    # [ON-DISK] refractiveindex_info/refractiveindex.info-database-main.zip @ member="database/data/specs/schott/optical/N-SF2.yml" wavelength=589nm via="C0/x" precision=3sf
    #   n = 1.6475512 at 293 K for Schott N-SF2. The row names a FAMILY ("dense"),
    #   not a material; N-SF2 is the dataset this citation resolves against, chosen by
    #   the implementer. Referred to the repo owner (C3.5).
    "Flint Glass (dense)": C0 / 1.65, 
    # [ON-DISK] refractiveindex_info/refractiveindex.info-database-main.zip @ member="database/data/main/Al2O3/nk/Malitson-o.yml" wavelength=589nm via="C0/x" precision=3sf
    #   n = 1.7680928 at 293 K for the ORDINARY ray. Sapphire is birefringent and the
    #   row names no ray: the extraordinary ray gives 1.7600182, which is 1.76.
    "Sapphire": C0 / 1.77,            
    # [UNVERIFIED] names a class of materials, not one material with one index, and no
    #   dataset in the archive carries that name. Its amorphous-material datasets span
    #   n = 1.74 (other/amorphous/HAC/nk/Smith-Td250C.yml) to 3.88
    #   (main/Ge2Sb2Te5/nk/Frantz-amorphous.yml) at 589 nm, over 18 datasets - no one
    #   value is the row. Referred to the repo owner (C3.5).
    "Glass (amorphous semiconductor)": C0 / 1.8, 
    # [ON-DISK] refractiveindex_info/refractiveindex.info-database-main.zip @ member="database/data/main/C/nk/Peter.yml" wavelength=589nm via="C0/x" precision=3sf
    #   n = 2.4172982 (Peter; no temperature stated in the file); Phillip's tabulated 2.4166020, between rows
    #   2.4209 and 2.4114, also 2.42.
    "Diamond": C0 / 2.42,             
    # [KNOWN-DEFECTIVE] 3.5, against 3.36-3.41 from all five on-disk datasets that
    #   cover 589 nm (Adachi 3.4068, Aspnes 3.3774, Bond 3.3616, Jellison 3.3768, Khmelevskaia 3.3916):
    #   the row is +2.7% to +4.1% high, and every dataset is 3.4 at
    #   the row's own 1-dp precision. A correction moves emitted items (P6), so it is
    #   referred to the repo owner (C3.5), not made here.
    "Gallium Phosphide (GaP)": C0 / 3.5, 
    
    # Special Cases (Important for RF/Microwave Engineering)
    # [UNVERIFIED] a microwave value, not an optical index: 7.14 is the square root of
    #   the relative permittivity ~51 the row's own comment states (7.1414). No
    #   text, JSON or HTML source under docs/references names tissue or muscle (PDF text
    #   layers not searched). The archive's human-body datasets
    #   are optical, none is muscle, and they give n = 1.34-1.58 at 589 nm
    #   - a different quantity from the row's. Residual register, C3.5.
    "Human Body Tissue (muscle, ~3 GHz)": C0 / 7.14, # Relative permittivity ε_r ~51, n=√ε_r
}

# Permittivity of free space in Farads per meter (F/m)
# @kind: measured-constant
# @units: F/m
# @domain: none (a fundamental constant; CODATA's relative uncertainty, 1.6e-10, is below every display precision in the corpus)
# [ON-DISK] codata_2022/allascii.txt @ quantity="vacuum electric permittivity" precision=4sf
#   A ROUNDING of CODATA 2022's 8.854 187 8188e-12, not a transcription of it,
#   and the precision is part of the citation: the resolver checks the 4-s.f.
#   rounding exactly (D-071). It sits 2.12e-5 relative below CODATA. No
#   consuming question states epsilon_0 - the four electrostatics templates
#   print it only in the solution, at .4e or .3e, both lossless for a 4-s.f.
#   value - so a solver who uses the full CODATA value reaches an answer
#   2.1e-5 relative from gold (phaseC1_summary.md, C1.7).
EPSILON_0 = 8.854e-12

# Value ranges for random parameter generation to ensure diverse problems.
# Frequencies are kept as integers for clarity in the problem statement.
# @kind: range
# @units: Hz
# @given: stated (tests/constants_integrity/given_evidence.py: every drawn frequency is printed in the question, 40/40 seeds)
# @domain: none (a sampling window for generated problems, not a property measured at conditions)
# [POLICY: sampling-only] a window the C1.3 census classifies PLAUSIBILITY; census --check holds the tag to that class
FREQUENCY_RANGE_HZ = (50, 2000)
# @kind: range
# @units: 1
# @domain: none (a sampling window for generated problems, not a property measured at conditions)
# [POLICY: sampling-only] a window the C1.3 census classifies PLAUSIBILITY; census --check holds the tag to that class
AMPLITUDE_RANGE = (1.0, 50.0)
# @kind: range
# @units: deg
# @given: stated (tests/constants_integrity/given_evidence.py: each drawn phase is printed inside the signal expression the question states - 40/40 in the Nyquist template, and 24/24 of the seeds that draw a degree phase in the conversion template)
# @domain: none (a sampling window for generated problems, not a property measured at conditions)
# [POLICY: sampling-only] a window the C1.3 census classifies PLAUSIBILITY; census --check holds the tag to that class
PHASE_RANGE_DEG = (-180, 180)
# @kind: range
# @units: rad
# @domain: none (a sampling window for generated problems, not a property measured at conditions)
# [POLICY: sampling-only] a window the C1.3 census classifies PLAUSIBILITY; census --check holds the tag to that class
PHASE_RANGE_RAD = (-math.pi, math.pi)

# The continuous frequency Omega will be a multiple of pi. This range defines the multiplier.
# @kind: range
# @units: 1
# @domain: none (a sampling window for generated problems, not a property measured at conditions)
# [POLICY: sampling-only] a window the C1.3 census classifies PLAUSIBILITY; census --check holds the tag to that class
OMEGA_MULTIPLIER_RANGE = (100, 1000)

# @kind: range
# @units: Hz
# @domain: none (a sampling window for generated problems, not a property measured at conditions)
# [POLICY: sampling-only] a window the C1.3 census classifies PLAUSIBILITY; census --check holds the tag to that class
SAMPLING_FREQ_RANGE_HZ = (1000, 8000)

# @kind: range
# @units: Hz
# @domain: none (a sampling window for generated problems, not a property measured at conditions)
# [POLICY: sampling-only] a window the C1.3 census classifies PLAUSIBILITY; census --check holds the tag to that class
F0_RANGE_HZ = (500, 3000) 

# The gain of the discrete-time system
# @kind: range
# @units: 1
# @domain: none (a sampling window for generated problems, not a property measured at conditions)
# [POLICY: sampling-only] a window the C1.3 census classifies PLAUSIBILITY; census --check holds the tag to that class
GAIN_K_RANGE = (0.5, 5.0)

# The delay (in samples) of the discrete-time system
# @kind: range
# @units: sample
# @domain: none (a sampling window for generated problems, not a property measured at conditions)
# [POLICY: sampling-only] a window the C1.3 census classifies PLAUSIBILITY; census --check holds the tag to that class
DELAY_N0_RANGE = (1, 10)


# The integer factor by which the signal is downsampled.
# @kind: range
# @units: 1
# @given: stated (tests/constants_integrity/given_evidence.py: M is printed in the question, 40/40)
# @domain: none (a sampling window for generated problems, not a property measured at conditions)
# [POLICY: sampling-only] a window the C1.3 census classifies PLAUSIBILITY; census --check holds the tag to that class
DECIMATION_FACTOR_M_RANGE = (2, 5)

# Define the pool for denominators of the omega_0 fraction.
# Using larger numbers allows for more granularity in creating frequencies.
# @kind: range
# @units: 1
# @given: stated (tests/constants_integrity/given_evidence.py: the draw shapes omega_0, which the question states exactly as a reduced fraction, 40/40; the denominator itself is not printed)
# @domain: none (a sampling window for generated problems, not a property measured at conditions)
# [POLICY: sampling-only] a window the C1.3 census classifies PLAUSIBILITY; census --check holds the tag to that class
OMEGA_DENOMINATOR_RANGE = (8, 20)
