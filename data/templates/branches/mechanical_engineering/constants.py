# E is the Modulus of Elasticity.
MATERIAL_PROPERTIES = {
    # Metals (Common)
    'Steel': {'E_GPa': 200, 'E_ksi': 29000, 'nu': 0.30},
    'Stainless Steel': {'E_GPa': 190, 'E_ksi': 27500, 'nu': 0.30}, # 304 Stainless
    'Aluminum': {'E_GPa': 69, 'E_ksi': 10000, 'nu': 0.33}, # General purpose 1100, 3003 Al
    'Aluminum 6061-T6': {'E_GPa': 68.9, 'E_ksi': 10000, 'nu': 0.33}, # A common specific alloy
    'Copper': {'E_GPa': 117, 'E_ksi': 17000, 'nu': 0.34},
    'Brass': {'E_GPa': 102, 'E_ksi': 14800, 'nu': 0.34}, # Yellow Brass (CuZn37)
    'Bronze': {'E_GPa': 110, 'E_ksi': 16000, 'nu': 0.34}, # Phosphor Bronze
    'Titanium': {'E_GPa': 116, 'E_ksi': 16800, 'nu': 0.34}, # Commercially Pure
    'Titanium Alloy (6Al-4V)': {'E_GPa': 114, 'E_ksi': 16500, 'nu': 0.34},
    'Magnesium': {'E_GPa': 45, 'E_ksi': 6500, 'nu': 0.29}, # AZ31B alloy
    'Tungsten': {'E_GPa': 411, 'E_ksi': 59600, 'nu': 0.28},
    'Cast Iron': {'E_GPa': 170, 'E_ksi': 24600, 'nu': 0.26}, # Gray Cast Iron
    'Nickel': {'E_GPa': 207, 'E_ksi': 30000, 'nu': 0.31},
    'Lead': {'E_GPa': 16, 'E_ksi': 2300, 'nu': 0.44}, # Added a common soft metal
    
    # Polymers/Plastics
    'Nylon': {'E_GPa': 2.1, 'E_ksi': 300, 'nu': 0.39}, # Nylon 6/6
    'Polycarbonate': {'E_GPa': 2.3, 'E_ksi': 334, 'nu': 0.38},
    'ABS': {'E_GPa': 2.0, 'E_ksi': 290, 'nu': 0.35},
    'PVC (rigid)': {'E_GPa': 3.0, 'E_ksi': 435, 'nu': 0.38},
    'PTFE (Teflon)': {'E_GPa': 0.5, 'E_ksi': 72, 'nu': 0.46},
    'Polyethylene (HDPE)': {'E_GPa': 1.1, 'E_ksi': 160, 'nu': 0.42},
    'Epoxy': {'E_GPa': 3.0, 'E_ksi': 435, 'nu': 0.38}, # Unreinforced epoxy resin
    'Natural Rubber': {'E_GPa': 0.0015, 'E_ksi': 0.22, 'nu': 0.4999}, # ~0.5 (nearly incompressible)
    
    # Composites & Other
    'Carbon Fiber Reinforced Polymer (CFRP)': {'E_GPa': 150, 'E_ksi': 21750, 'nu': 0.30}, # Direction-dependent, this is a typical in-plane value
    'Fiberglass (GFRP)': {'E_GPa': 45, 'E_ksi': 6500, 'nu': 0.25}, # Direction-dependent, typical in-plane value
    'Concrete': {'E_GPa': 30, 'E_ksi': 4350, 'nu': 0.15}, # Highly variable, common approx.
    'Glass (Borosilicate)': {'E_GPa': 70, 'E_ksi': 10150, 'nu': 0.20},
    'Ceramic (Alumina Al2O3)': {'E_GPa': 370, 'E_ksi': 53700, 'nu': 0.22},
}


# Typical Shear Modulus (G) Values for Engineering Materials (GPa)
# Note: Values are representative and can vary with alloy composition and heat treatment.
SHEAR_MODULUS_VALUES = {
    # Metals (Common Alloys)
    "Steel (A36)": 77.2,
    "Stainless Steel (304)": 77.0,
    "Aluminum 6061-T6": 26.0,
    "Aluminum 2024-T4": 28.0,
    "Brass (C36000)": 39.0,
    "Copper (CDA 110)": 44.7,
    "Bronze (Phosphor 510)": 41.4,
    "Titanium Alloy (Ti-6Al-4V)": 41.4,
    "Magnesium Alloy (AZ31B)": 16.5,
    "Tungsten": 160.0,
    "Molybdenum": 118.0,
    "Lead": 5.9,
    
    # Other Materials
    "Gray Cast Iron": 44.0, # Note: Anisotropic property, value is an approximation.
    "Nylon 6/6": 1.1,      # Polymers have much lower moduli
    "Polycarbonate": 0.9,
    "Concrete": 22.0,      # Highly variable, this is a common estimate for calculation
    "Glass": 26.2,
    
    # Composites (Highly variable based on fiber orientation and volume fraction)
    # Value given is an approximate in-plane shear modulus.
    "Carbon Fiber Epoxy (Unidirectional, in-plane)": 5.0,
}

# Standard atmospheric pressure in kPa
ATMOSPHERIC_PRESSURE_KPA = 101.325

# Standard acceleration due to gravity in m/s^2
GRAVITY = 9.81


# Densities of common fluids in kg/m^3 at 20°C and 1 atm, unless specified.
FLUID_DENSITIES = {
    # Water-based
    "Fresh Water": 998,         # More precise value for 20°C
    "Sea Water": 1025,
    "Brine (20% NaCl)": 1150,   # Common industrial solution
    
    # Oils and Hydrocarbons
    "SAE 10W Oil": 870,
    "SAE 30 Oil": 917,
    "Crude Oil": 870,           # Typical average value
    "Kerosene": 810,
    "Diesel Fuel": 850,
    "Gasoline": 726,
    "Engine Oil": 888,
    "Hydraulic Oil": 850,       # Typical average value
    
    # Alcohols and Solvents
    "Ethyl Alcohol (Ethanol)": 789,
    "Methyl Alcohol (Methanol)": 791,
    "Isopropyl Alcohol (IPA)": 786,
    "Acetone": 791,
    "Turpentine": 870,
    "Benzene": 876,
    "Toluene": 867,
    "Xylene": 870,
    
    # Cryogens & Liquefied Gases
    "Liquid Nitrogen (at -196°C)": 807,
    "Liquid Oxygen (at -183°C)": 1141,
    "Liquid Hydrogen (at -253°C)": 71,
    "Liquid Propane": 493,      # At 25°C, under pressure
    
    # High-Density Liquids
    "Mercury": 13550,
    "Glycerin": 1260,
    "Chloroform": 1489,
    "Carbon Tetrachloride": 1594,
    "Bromine": 3120,
    "Sulfuric Acid (98%)": 1830,
    "Hydrochloric Acid (37%)": 1190, # Also known as Muriatic Acid
    "Nitric Acid (68%)": 1350,
    
    # Common Liquids
    "Milk": 1035,
    "Whole Blood": 1060,        # Typical value for human blood
    "Ethylene Glycol (Antifreeze)": 1113,
    "Corn Syrup": 1380,
    "Olive Oil": 910,
    "Vegetable Oil": 920,
    "R-134a Refrigerant (Saturated Liquid at 25°C)": 1206,
    
    # Gases (at 1 atm, for comparison - though not for hydrostatic problems!)
    "Air (at 20°C)": 1.2,
    "Helium (at 20°C)": 0.166,
    "Carbon Dioxide (at 20°C)": 1.84
}

# Densities of various solid materials in kg/m^3
MATERIAL_DENSITIES = {
    # Woods & Natural Materials (generally float in water)
    "Pine Wood": 500,
    "Oak Wood": 750,
    "Cork": 240,
    "Teak Wood": 630,
    "Maple Wood": 740,
    "Ebony Wood": 1200,  # Sinks in water
    "Bamboo": 300,
    "Rubber (Natural)": 950,
    
    # Plastics & Synthetics (range from floating to sinking)
    "Polypropylene (PP)": 900,
    "Polyethylene (LDPE)": 920,
    "Polyethylene (HDPE)": 950,
    "Nylon": 1150,
    "Polystyrene (PS)": 1050,
    "PVC (Rigid)": 1380,
    "PTFE (Teflon)": 2200,
    "Acrylic (Plexiglas)": 1180,
    
    # Biological Materials
    "Human Body (avg.)": 985,  # Slightly less than water; most people float
    "Apple": 800,
    "Bone": 1850,
    
    # Lightweight & Porous Materials
    "Aerogel": 3,
    "Pumice": 700,  # The only rock that floats!
    "Brick": 1700,
    "Concrete (avg.)": 2400,
    
    # Common Metals & Alloys (will sink in water, float in mercury)
    "Aluminum": 2710,
    "Titanium": 4500,
    "Zinc": 7140,
    "Tin": 7280,
    "Iron (Wrought)": 7750,
    "Steel (Carbon)": 7850,
    "Brass": 8600,
    "Copper": 8940,
    "Nickel": 8900,
    "Silver": 10500,
    "Lead": 11340,
    "Uranium": 19100,
    "Gold": 19300,
    "Tungsten": 19600,
    "Platinum": 21450,
    "Osmium": 22590,  # The densest naturally occurring element
    
    # Other Materials
    "Ice (0°C)": 917,
    "Wax (Paraffin)": 900,
    "Glass (Window)": 2500,
    "Glass (Crystal)": 2900,
    "Graphite": 2100,
    "Diamond": 3500,
    "Quartz": 2650,
    "Granite": 2700,
    "Marble": 2700,
    "Sandstone": 2300,
    "Limestone": 2500,
    "Cork Board": 240,
    "Paper": 800,
    "Leather (Dry)": 860,
    "Chalk": 2500,
    "Asphalt": 1100,
    "Coal (Anthracite)": 1500,
}

OBJECT_SHAPES = [
    "sphere", "cube", "irregular block", "cylinder", "anchor",
    "rectangular prism", "cone", "pyramid", "torpedo", "submarine hull",
    "ball", "ingot", "pipe section", "bead", "marble",
    "statue", "metal ball", "wooden block", "concrete piling", "boat hull",
    "metal rod", "canister", "storage tank", "buoy", "weight",
    "metal box", "stone", "boulder", "metal cylinder", "plastic container"
]

OBJECT_MATERIALS = [
    "steel", "aluminum", "copper", "concrete", "plastic",
    "wood (oak)", "wood (pine)", "glass", "gold", "lead",
    "brass", "bronze", "iron", "titanium", "rubber",
    "cork", "foam", "ice", "wax", "ceramic",
    "polyethylene", "PVC", "fiberglass", "granite", "marble",
    "ebony", "balsa wood", "uranium", "magnesium", "tungsten"
]

# Densities of common fluids in kg/m^3
# Lighter fluids suitable for the pipe
PIPE_FLUIDS = {
    # Gases (Very Light)
    "Air": 1.225,
    "Helium": 0.1786,
    "Hydrogen": 0.0899,
    "Natural Gas (Methane)": 0.717,
    "Carbon Dioxide": 1.98,
    
    # Light Hydrocarbons & Fuels
    "Gasoline": 726,
    "Kerosene": 810,
    "Diesel Fuel": 850,
    "Jet Fuel (JP-4)": 770,
    "Ethanol": 789,
    "Methanol": 791,
    "Isopropyl Alcohol": 786,
    
    # Oils & Lubricants
    "SAE 10 Oil": 870,
    "SAE 20 Oil": 880,
    "SAE 30 Oil": 917,
    "SAE 40 Oil": 945,
    "SAE 50 Oil": 960,
    "Crude Oil (Light)": 825,
    "Crude Oil (Heavy)": 950,
    "Engine Oil": 888,
    "Hydraulic Oil": 860,
    
    # Common Liquids
    "Water": 1000,
    "Sea Water": 1025,
    "Milk": 1030,
    "Ethylene Glycol": 1110,
    "Antifreeze (50/50)": 1065,
    
    # Cryogenic Fluids
    "Liquid Nitrogen": 810,
    "Liquid Oxygen": 1141,
    
    # Chemical Solvents
    "Acetone": 784,
    "Toluene": 867,
    "Benzene": 876,
    "Hexane": 655
}

# Denser fluids suitable for the manometer
MANOMETER_FLUIDS = {
    # Standard Manometer Fluids
    "Mercury": 13550,
    "Water": 1000,  # For gas measurements
    "Sea Water": 1025,
    
    # Heavy Oils & Lubricants
    "SAE 50 Oil": 960,
    "SAE 90 Gear Oil": 920,
    
    # Organic Liquids
    "Carbon Tetrachloride": 1590,
    "Chloroform": 1480,
    "Bromine": 3120,
    "Tetrabromoethane": 2960,
    
    # Inorganic Solutions
    "Calcium Chloride Solution (40%)": 1390,
    "Zinc Chloride Solution (50%)": 1520,
    "Sodium Polysulfide": 1650,
    
    # Specialty Fluids
    "Glycerin": 1260,
    "Diiodomethane": 3325,
    "Acetylene Tetrabromide": 2960,
    
    # Molten Metals (for high-temperature applications)
    "Gallium": 6095,
    "Tin": 6980,
    "Zinc": 6570,
    
    # Very Heavy Fluids
    "Tellurium Mercury": 8100,
    "Tungsten Hexafluoride": 12900
}
