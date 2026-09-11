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
MEDIA_VELOCITIES = {
    # Gases (at 0°C and 1 atm, for visible light ~589 nm)
    "Vacuum": C0,                     
    "Air (at sea level)": C0 / 1.000293, 
    "Helium": C0 / 1.000036,          
    "Carbon Dioxide": C0 / 1.00045,   
    
    # Liquids (for visible light ~589 nm)
    "Water (distilled, 20°C)": C0 / 1.333,  
    "Ethanol": C0 / 1.36,             
    "Glycerine": C0 / 1.473,          
    "Benzene": C0 / 1.501,            
    "Carbon Disulfide": C0 / 1.628,   # notable for high dispersion
    
    # Solids (for visible light ~589 nm)
    "Ice": C0 / 1.31,                 
    "Teflon (PTFE)": C0 / 1.35,       
    "Fused Silica (Glass)": C0 / 1.458, 
    "Crown Glass (typical)": C0 / 1.52, 
    "Polyethylene": C0 / 1.54,        
    "Polystyrene": C0 / 1.59,         
    "Flint Glass (dense)": C0 / 1.65, 
    "Sapphire": C0 / 1.77,            
    "Glass (amorphous semiconductor)": C0 / 1.8, 
    "Diamond": C0 / 2.42,             
    "Gallium Phosphide (GaP)": C0 / 3.5, 
    
    # Special Cases (Important for RF/Microwave Engineering)
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
FREQUENCY_RANGE_HZ = (50, 2000)
# @kind: range
# @units: 1
AMPLITUDE_RANGE = (1.0, 50.0)
# @kind: range
# @units: deg
# @given: stated (tests/constants_integrity/given_evidence.py: each drawn phase is printed inside the signal expression the question states - 40/40 in the Nyquist template, and 24/24 of the seeds that draw a degree phase in the conversion template)
PHASE_RANGE_DEG = (-180, 180)
# @kind: range
# @units: rad
PHASE_RANGE_RAD = (-math.pi, math.pi)

# The continuous frequency Omega will be a multiple of pi. This range defines the multiplier.
# @kind: range
# @units: 1
OMEGA_MULTIPLIER_RANGE = (100, 1000)

# @kind: range
# @units: Hz
SAMPLING_FREQ_RANGE_HZ = (1000, 8000)

# @kind: range
# @units: Hz
F0_RANGE_HZ = (500, 3000) 

# The gain of the discrete-time system
# @kind: range
# @units: 1
GAIN_K_RANGE = (0.5, 5.0)

# The delay (in samples) of the discrete-time system
# @kind: range
# @units: sample
DELAY_N0_RANGE = (1, 10)


# The integer factor by which the signal is downsampled.
# @kind: range
# @units: 1
# @given: stated (tests/constants_integrity/given_evidence.py: M is printed in the question, 40/40)
DECIMATION_FACTOR_M_RANGE = (2, 5)

# Define the pool for denominators of the omega_0 fraction.
# Using larger numbers allows for more granularity in creating frequencies.
# @kind: range
# @units: 1
# @given: stated (tests/constants_integrity/given_evidence.py: the draw shapes omega_0, which the question states exactly as a reduced fraction, 40/40; the denominator itself is not printed)
OMEGA_DENOMINATOR_RANGE = (8, 20)
