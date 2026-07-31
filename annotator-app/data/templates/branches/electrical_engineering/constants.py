import math

# Speed of light in vacuum (m/s)
C0 = 299792458 

# Dictionary of common media and their approximate phase velocities for EM waves
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
EPSILON_0 = 8.854e-12

# Value ranges for random parameter generation to ensure diverse problems.
# Frequencies are kept as integers for clarity in the problem statement.
FREQUENCY_RANGE_HZ = (50, 2000)
AMPLITUDE_RANGE = (1.0, 50.0)
PHASE_RANGE_DEG = (-180, 180)
PHASE_RANGE_RAD = (-math.pi, math.pi)

# The continuous frequency Omega will be a multiple of pi. This range defines the multiplier.
OMEGA_MULTIPLIER_RANGE = (100, 1000)

SAMPLING_FREQ_RANGE_HZ = (1000, 8000)

F0_RANGE_HZ = (500, 3000) 

# The gain of the discrete-time system
GAIN_K_RANGE = (0.5, 5.0)

# The delay (in samples) of the discrete-time system
DELAY_N0_RANGE = (1, 10)


# The integer factor by which the signal is downsampled.
DECIMATION_FACTOR_M_RANGE = (2, 5)

# Define the pool for denominators of the omega_0 fraction.
# Using larger numbers allows for more granularity in creating frequencies.
OMEGA_DENOMINATOR_RANGE = (8, 20)
