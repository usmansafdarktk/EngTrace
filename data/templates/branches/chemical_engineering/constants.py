# Standard: phases referenced at ~25 °C and 1 atm

LIQUID_PHASE_REACTANTS = [
    "Ethyl Acetate",
    "Propylene Glycol",
    "Benzene",
    "Toluene",
    "Acetone",
    "Methanol",
    "Ethanol",
    "Isopropanol",     # specify "isopropanol (2-propanol)" if needed
    "n-Butanol",       # specify isomer as needed
    "Methyl Ethyl Ketone (MEK) / 2-butanone",
    # Formaldehyde is typically a gas at 25°C; aqueous solution is "Formalin (37% aqueous)"
    "Acetic Acid",
    # Phenol is a solid at 25°C (mp ~40.5°C). Include only if working above that temperature.
    "Glycerol",
    "Xylene",
    "Styrene",
    "Aniline",
    "Cyclohexane",
    "Formic Acid",
    "Nitric Acid",
    "Sulfuric Acid",
    "Ethylene Glycol",
    "Diethyl Ether",
    "Tetrahydrofuran (THF)",
    "Chloroform",
    "Acrylonitrile",
    "Dimethylformamide (DMF)",
    "Methyl Methacrylate",
    "Acetic anhydride"   # corrected name (a.k.a. ethanoic anhydride)
]

GAS_PHASE_REACTANTS = [
    "Methane",
    "Ethane",
    "Propane",
    "Ammonia",
    "Ethylene",
    "Sulfur Dioxide",
    "Hydrogen Sulfide",
    "Vinyl Chloride",
    "Butadiene",
    "Hydrogen",
    "Oxygen",
    "Nitrogen",
    "Chlorine",
    "Carbon Monoxide",
    "Carbon Dioxide",
    "Propylene",
    "Butane",
    "Acetylene",
    "Nitric Oxide (NO)",
    "Nitrogen Dioxide (NO2)",
    "Hydrogen Chloride (HCl)",   # gas at 1 atm
    "Phosgene",
    "Ethylene Oxide",
    "Isobutane",
    "Formaldehyde",     # moved here: gas at 25°C
    "Acetaldehyde"      # volatile; bp ~20.2°C → effectively gas at 25°C
]

BIOCHEMICAL_SUBSTRATES = [
    "Glucose",
    "Sucrose",
    "Lactose",
    "Fructose",
    "Pyruvate",
    "Maltose",
    "Galactose",
    "Starch",        # polymeric solid
    "Cellulose",     # polymeric solid
    "Xylose",
    "Glutamate",     # often used as salts (e.g., sodium glutamate)
    "Alanine",
    "Lactate",       # usually as salt
    "Citrate",
    "Acetyl-CoA",    # coenzyme (large, charged)
    "Palmitic Acid",
    "Oleic Acid",
    "Triglycerides",
    "Urea",
    "Aspartate"
]

GENERAL_REACTANTS = LIQUID_PHASE_REACTANTS + GAS_PHASE_REACTANTS + BIOCHEMICAL_SUBSTRATES


PRODUCTS = [
    "Product Alpha",
    "Product Beta",
    "Product Gamma",
    "Product Delta",
    "Product Sigma",
    "Product Omega",
    "Product Theta",
    "Product Lambda",
    "Product Zeta",
    "Product Kappa",

    "Compound P",
    "Compound Q",
    "Compound R",
    "Compound S",
    "Compound T",
    "Compound V",
    "Compound W",
    "Compound X",
    "Compound Y",
    "Compound Z",

    "Species I",
    "Species II",
    "Species III",
    "Species IV",
    "Species V",
    "Species VI",
    "Species VII",
    "Species VIII",

    "Material A",
    "Material B",
    "Material C",
    "Material D",
    "Material E",
    "Material F",

    "Substance One",
    "Substance Two",
    "Substance Three",
    "Substance Four",
    "Substance Five",
    "Substance Six",

    "Molecule M",
    "Molecule N",
    "Molecule O",
    "Molecule P",
    "Molecule Q",
    "Molecule R",

    "Entity 1",
    "Entity 2",
    "Entity 3",
    "Entity 4",
    "Entity 5",
    "Entity 6"
]


# A list of common substances used in thermodynamics problems involving phase change.
THERMO_SUBSTANCES = [
    # Classic working fluids
    "Water",
    "Ammonia",
    "Carbon Dioxide",
    "Sulfur Dioxide",

    # Hydrocarbons (fuels and refrigerants)
    "Methane",
    "Ethane",
    "Propane",
    "Butane",
    "Isobutane",
    "Pentane",
    "Iso-pentane",

    # Refrigerants (with ASHRAE designations)
    "Refrigerant-11 (R-11, Trichlorofluoromethane)",
    "Refrigerant-12 (R-12, Dichlorodifluoromethane)",
    "Refrigerant-22 (R-22, Chlorodifluoromethane)",
    "Refrigerant-134a (R-134a, 1,1,1,2-Tetrafluoroethane)",
    "Refrigerant-123 (R-123, Dichlorotrifluoroethane)",
    "Refrigerant-410A (R-410A, blend of difluoromethane and pentafluoroethane)",

    # Common industrial/organic fluids used in Rankine/Organic Rankine cycles
    "Toluene",
    "Benzene",
    "Ethanol",
    "Methanol",
    "Acetone",
    "n-Hexane",
    "n-Octane",
    "Cyclohexane"
]


# A dictionary with comprehensive critical properties for various substances.
# Tc: Kelvin (K), Pc: bar, Vc: cm³/mol, Zc: dimensionless, omega: dimensionless.
CRITICAL_PROPERTIES = {
    "Methane": {"Tc": 190.6, "Pc": 45.99, "Vc": 99.0, "Zc": 0.286, "omega": 0.012},
    "Ethane": {"Tc": 305.3, "Pc": 48.72, "Vc": 146.0, "Zc": 0.279, "omega": 0.100},
    "Propane": {"Tc": 369.8, "Pc": 42.48, "Vc": 200.0, "Zc": 0.276, "omega": 0.152},
    "n-Butane": {"Tc": 425.1, "Pc": 37.96, "Vc": 255.0, "Zc": 0.274, "omega": 0.200},
    "n-Pentane": {"Tc": 469.7, "Pc": 33.7, "Vc": 311.0, "Zc": 0.269, "omega": 0.251},
    "n-Hexane": {"Tc": 507.6, "Pc": 30.12, "Vc": 368.0, "Zc": 0.264, "omega": 0.301},
    "n-Heptane": {"Tc": 540.2, "Pc": 27.36, "Vc": 426.0, "Zc": 0.263, "omega": 0.350},
    "n-Octane": {"Tc": 568.8, "Pc": 24.86, "Vc": 492.0, "Zc": 0.259, "omega": 0.400},
    "Ethylene": {"Tc": 282.4, "Pc": 50.42, "Vc": 131.0, "Zc": 0.281, "omega": 0.087},
    "Propylene": {"Tc": 365.0, "Pc": 46.0, "Vc": 181.0, "Zc": 0.275, "omega": 0.140},
    "Benzene": {"Tc": 562.2, "Pc": 48.95, "Vc": 259.0, "Zc": 0.271, "omega": 0.210},
    "Toluene": {"Tc": 591.8, "Pc": 41.09, "Vc": 316.0, "Zc": 0.264, "omega": 0.263},
    "p-Xylene": {"Tc": 616.2, "Pc": 35.12, "Vc": 379.0, "Zc": 0.260, "omega": 0.321},
    "Methanol": {"Tc": 512.6, "Pc": 80.84, "Vc": 118.0, "Zc": 0.224, "omega": 0.564},
    "Ethanol": {"Tc": 513.9, "Pc": 61.37, "Vc": 167.0, "Zc": 0.248, "omega": 0.645},
    "Acetone": {"Tc": 508.2, "Pc": 46.99, "Vc": 209.0, "Zc": 0.232, "omega": 0.304},
    "Water": {"Tc": 647.1, "Pc": 220.64, "Vc": 55.9, "Zc": 0.229, "omega": 0.345},
    "Ammonia": {"Tc": 405.5, "Pc": 113.53, "Vc": 72.5, "Zc": 0.242, "omega": 0.250},
    "Carbon dioxide": {"Tc": 304.2, "Pc": 73.83, "Vc": 94.0, "Zc": 0.274, "omega": 0.224},
    "Carbon monoxide": {"Tc": 132.9, "Pc": 34.99, "Vc": 93.1, "Zc": 0.295, "omega": 0.045},
    "Oxygen": {"Tc": 154.6, "Pc": 50.43, "Vc": 73.4, "Zc": 0.288, "omega": 0.022},
    "Nitrogen": {"Tc": 126.2, "Pc": 33.98, "Vc": 89.8, "Zc": 0.290, "omega": 0.039},
    "Hydrogen": {"Tc": 33.2, "Pc": 12.97, "Vc": 65.0, "Zc": 0.305, "omega": -0.216},
    "Helium": {"Tc": 5.2, "Pc": 2.27, "Vc": 57.8, "Zc": 0.301, "omega": -0.365},
    "Chlorine": {"Tc": 417.0, "Pc": 77.02, "Vc": 124.0, "Zc": 0.275, "omega": 0.090},
    "Sulfur dioxide": {"Tc": 430.8, "Pc": 78.84, "Vc": 122.0, "Zc": 0.269, "omega": 0.245},
    "Hydrogen sulfide": {"Tc": 373.5, "Pc": 89.63, "Vc": 98.6, "Zc": 0.284, "omega": 0.091}
}


# Common materials which undergo heating with their specific heat capacities in J/g·K
SUBSTANCES_FOR_HEATING = [
    # Metals & Solids
    {"name": "Iron", "state": "solid", "Cp": 0.449},
    {"name": "Copper", "state": "solid", "Cp": 0.385},
    {"name": "Aluminum", "state": "solid", "Cp": 0.897},
    {"name": "Gold", "state": "solid", "Cp": 0.129},
    {"name": "Lead", "state": "solid", "Cp": 0.16},
    {"name": "Silver", "state": "solid", "Cp": 0.235},
    {"name": "Tungsten", "state": "solid", "Cp": 0.134},
    {"name": "Silicon", "state": "solid", "Cp": 0.705},
    {"name": "Graphite (Carbon)", "state": "solid", "Cp": 0.709},
    {"name": "Glass (typical)", "state": "solid", "Cp": 0.84},
    {"name": "Ice (at 0°C)", "state": "solid", "Cp": 2.09},
    {"name": "Concrete", "state": "solid", "Cp": 0.88},
    {"name": "Wood (typical)", "state": "solid", "Cp": 1.7},
    {"name": "Polyethylene (plastic)", "state": "solid", "Cp": 2.3},

    # Liquids
    {"name": "Water", "state": "liquid", "Cp": 4.18},
    {"name": "Ethanol", "state": "liquid", "Cp": 2.44},
    {"name": "Methanol", "state": "liquid", "Cp": 2.53},
    {"name": "Acetone", "state": "liquid", "Cp": 2.17},
    {"name": "Mercury", "state": "liquid", "Cp": 0.14},
    {"name": "Glycerol", "state": "liquid", "Cp": 2.43},
    {"name": "Ethylene Glycol (Antifreeze)", "state": "liquid", "Cp": 2.36},
    {"name": "Olive Oil", "state": "liquid", "Cp": 1.97},
    {"name": "Engine Oil (typical)", "state": "liquid", "Cp": 1.9},
    {"name": "Sulfuric Acid", "state": "liquid", "Cp": 1.42},

    # Gases (at constant pressure, 25°C)
    {"name": "Air (dry)", "state": "gas", "Cp": 1.005},
    {"name": "Nitrogen", "state": "gas", "Cp": 1.04},
    {"name": "Oxygen", "state": "gas", "Cp": 0.918},
    {"name": "Hydrogen", "state": "gas", "Cp": 14.31},
    {"name": "Helium", "state": "gas", "Cp": 5.193},
    {"name": "Argon", "state": "gas", "Cp": 0.520},
    {"name": "Carbon Dioxide", "state": "gas", "Cp": 0.839},
    {"name": "Methane", "state": "gas", "Cp": 2.22},
    {"name": "Ammonia", "state": "gas", "Cp": 2.06},
    {"name": "Water Vapor (Steam, 100°C)", "state": "gas", "Cp": 2.01},
    {"name": "Chlorine", "state": "gas", "Cp": 0.48}
]


# A list of common substances with their molar heats of vaporization (delta_H_vap)
# at their normal boiling points. All values are in kJ/mol.
SUBSTANCES_FOR_VAPORIZATION = [
    # Alcohols & Water
    {"name": "Water", "delta_H_vap": 40.66},
    {"name": "Methanol", "delta_H_vap": 35.3},
    {"name": "Ethanol", "delta_H_vap": 38.6},
    {"name": "Isopropanol", "delta_H_vap": 39.85},

    # Alkanes
    {"name": "Propane", "delta_H_vap": 19.04},
    {"name": "n-Butane", "delta_H_vap": 22.44},
    {"name": "n-Hexane", "delta_H_vap": 28.85},

    # Organic Solvents
    {"name": "Acetone", "delta_H_vap": 29.1},
    {"name": "Benzene", "delta_H_vap": 30.8},
    {"name": "Toluene", "delta_H_vap": 33.48},
    {"name": "Carbon Tetrachloride", "delta_H_vap": 29.82},

    # Inorganic & Elemental Substances
    {"name": "Ammonia", "delta_H_vap": 23.3},
    {"name": "Mercury", "delta_H_vap": 59.11},
    {"name": "Nitrogen", "delta_H_vap": 5.57},
    {"name": "Oxygen", "delta_H_vap": 6.82},
    {"name": "Argon", "delta_H_vap": 6.43},
]


# Database of standard heats of formation (ΔH_f°) at 298.15 K in kJ/mol.
# A value of 0 indicates an element in its standard state.
HEATS_OF_FORMATION = {
    # Hydrocarbons
    "CH4(g)": -74.8,      # Methane
    "C2H6(g)": -84.7,     # Ethane
    "C3H8(g)": -103.8,    # Propane
    "C6H6(l)": 49.0,      # Benzene
    # Alcohols
    "CH3OH(l)": -238.6,   # Methanol
    "C2H5OH(l)": -277.7,  # Ethanol
    # Common Gases & Products
    "O2(g)": 0,
    "H2(g)": 0,
    "N2(g)": 0,
    "CO(g)": -110.5,
    "CO2(g)": -393.5,
    "H2O(g)": -241.8,
    "H2O(l)": -285.8,
    "NH3(g)": -46.1,      # Ammonia
    "NO(g)": 90.3,
    "NO2(g)": 33.2,
}

# A list of predefined, balanced chemical reactions.
REACTIONS = [
    {
        "name": "Combustion of Methane",
        "equation": "CH4(g) + 2O2(g) → CO2(g) + 2H2O(l)",
        "reactants": {"CH4(g)": 1, "O2(g)": 2},
        "products": {"CO2(g)": 1, "H2O(l)": 2}
    },
    {
        "name": "Combustion of Propane",
        "equation": "C3H8(g) + 5O2(g) → 3CO2(g) + 4H2O(l)",
        "reactants": {"C3H8(g)": 1, "O2(g)": 5},
        "products": {"CO2(g)": 3, "H2O(l)": 4}
    },
    {
        "name": "Oxidation of Ammonia",
        "equation": "4NH3(g) + 5O2(g) → 4NO(g) + 6H2O(g)",
        "reactants": {"NH3(g)": 4, "O2(g)": 5},
        "products": {"NO(g)": 4, "H2O(g)": 6}
    },
    {
        "name": "Steam Reforming of Methane",
        "equation": "CH4(g) + H2O(g) → CO(g) + 3H2(g)",
        "reactants": {"CH4(g)": 1, "H2O(g)": 1},
        "products": {"CO(g)": 1, "H2(g)": 3}
    }
]


# A dictionary of substances with their heat capacity parameters for the
# equation: Cp/R = A + B*T + C*T² + D*T⁻² where T is in Kelvin.
# Parameters are typically valid for temperature ranges around 298-1200K
CP_PARAMS = {
    # Original entries
    "CH4(g)": {"A": 1.702, "B": 9.081E-2, "C": -2.164E-5, "D": 0},
    "CO2(g)": {"A": 5.457, "B": 1.045E-2, "C": 0, "D": -1.157E5},
    "N2(g)": {"A": 3.280, "B": 0.593E-2, "C": 0, "D": 0.040E5},
    "H2O(l)": {"A": 8.712, "B": 1.25E-2, "C": -0.18E-5, "D": 0},
    "C2H5OH(g)": {"A": 3.518, "B": 20.001E-2, "C": -6.002E-5, "D": 0},
    
    # Common gases
    "O2(g)": {"A": 3.630, "B": 1.794E-2, "C": -0.658E-5, "D": 0.061E5},
    "H2(g)": {"A": 3.249, "B": 0.422E-2, "C": 0, "D": 0.083E5},
    "NH3(g)": {"A": 3.578, "B": 3.020E-2, "C": 0, "D": -0.186E5},
    "CO(g)": {"A": 3.376, "B": 0.557E-2, "C": 0, "D": -0.031E5},
    "H2S(g)": {"A": 3.931, "B": 1.490E-2, "C": 0, "D": -0.232E5},
    
    # Hydrocarbons
    "C2H6(g)": {"A": 1.131, "B": 19.225E-2, "C": -5.561E-5, "D": 0},
    "C3H8(g)": {"A": 1.213, "B": 28.785E-2, "C": -8.824E-5, "D": 0},
    "C4H10(g)": {"A": 1.935, "B": 36.915E-2, "C": -11.402E-5, "D": 0},
    "C2H4(g)": {"A": 1.424, "B": 14.394E-2, "C": -4.392E-5, "D": 0},
    "C2H2(g)": {"A": 6.132, "B": 8.914E-2, "C": -6.347E-5, "D": 0},
    
    # Common liquids
    "C6H6(l)": {"A": -0.206, "B": 39.064E-2, "C": -13.301E-5, "D": 0},
    "C7H8(l)": {"A": 0.290, "B": 47.052E-2, "C": -15.716E-5, "D": 0},
    "C3H6O(l)": {"A": 1.506, "B": 30.476E-2, "C": -9.127E-5, "D": 0},
    "CH3OH(l)": {"A": 5.052, "B": 16.561E-2, "C": -3.761E-5, "D": 0},
    "C6H14(l)": {"A": 2.738, "B": 45.854E-2, "C": -14.518E-5, "D": 0},
    
    # Inorganic compounds
    "SO2(g)": {"A": 5.699, "B": 0.801E-2, "C": 0, "D": -1.015E5},
    "NO(g)": {"A": 3.387, "B": 0.669E-2, "C": 0, "D": 0.095E5},
    "NO2(g)": {"A": 4.982, "B": 1.195E-2, "C": -0.792E-5, "D": -0.377E5},
    "Cl2(g)": {"A": 4.442, "B": 0.089E-2, "C": 0, "D": -0.344E5},
    "HCl(g)": {"A": 3.156, "B": 0.623E-2, "C": 0, "D": 0.151E5},
    
    # Noble gases
    "He(g)": {"A": 2.500, "B": 0, "C": 0, "D": 0},
    "Ar(g)": {"A": 2.500, "B": 0, "C": 0, "D": 0},
    "Ne(g)": {"A": 2.500, "B": 0, "C": 0, "D": 0},
    
    # Additional common substances
    "Air(g)": {"A": 3.355, "B": 0.575E-2, "C": 0, "D": -0.016E5},
    "H2O(g)": {"A": 3.470, "B": 1.450E-2, "C": 0, "D": 0.121E5},
    "H2SO4(l)": {"A": 2.850, "B": 13.400E-2, "C": 0, "D": 0},
    "NaCl(s)": {"A": 5.526, "B": 1.963E-2, "C": 0, "D": -0.333E5},
    "CaCO3(s)": {"A": 12.572, "B": 2.637E-2, "C": -3.120E-5, "D": -3.642E5}
}


# A list of predefined, balanced combustion reactions with theoretical air.
COMBUSTION_REACTIONS = [
    {
        "name": "Combustion of Methane",
        "fuel": "Methane",
        "equation": "CH4(g) + 2O2(g) + 7.52N2(g) → CO2(g) + 2H2O(g) + 7.52N2(g)",
        "reactants": {"CH4(g)": 1, "O2(g)": 2, "N2(g)": 7.52},
        "products": {"CO2(g)": 1, "H2O(g)": 2, "N2(g)": 7.52}
    },
    {
        "name": "Combustion of Propane",
        "fuel": "Propane",
        "equation": "C3H8(g) + 5O2(g) + 18.8N2(g) → 3CO2(g) + 4H2O(g) + 18.8N2(g)",
        "reactants": {"C3H8(g)": 1, "O2(g)": 5, "N2(g)": 18.8},
        "products": {"CO2(g)": 3, "H2O(g)": 4, "N2(g)": 18.8}
    },
    {
        "name": "Combustion of Hydrogen",
        "fuel": "Hydrogen",
        "equation": "H2(g) + 0.5O2(g) + 1.88N2(g) → H2O(g) + 1.88N2(g)",
        "reactants": {"H2(g)": 1, "O2(g)": 0.5, "N2(g)": 1.88},
        "products": {"H2O(g)": 1, "N2(g)": 1.88}
    },
    {
        "name": "Combustion of Carbon Monoxide",
        "fuel": "Carbon Monoxide",
        "equation": "CO(g) + 0.5O2(g) + 1.88N2(g) → CO2(g) + 1.88N2(g)",
        "reactants": {"CO(g)": 1, "O2(g)": 0.5, "N2(g)": 1.88},
        "products": {"CO2(g)": 1, "N2(g)": 1.88}
    },
    {
        "name": "Combustion of Acetylene",
        "fuel": "Acetylene",
        "equation": "C2H2(g) + 2.5O2(g) + 9.4N2(g) → 2CO2(g) + H2O(g) + 9.4N2(g)",
        "reactants": {"C2H2(g)": 1, "O2(g)": 2.5, "N2(g)": 9.4},
        "products": {"CO2(g)": 2, "H2O(g)": 1, "N2(g)": 9.4}
    },
    {
        "name": "Combustion of Ethane",
        "fuel": "Ethane",
        "equation": "C2H6(g) + 3.5O2(g) + 13.16N2(g) → 2CO2(g) + 3H2O(g) + 13.16N2(g)",
        "reactants": {"C2H6(g)": 1, "O2(g)": 3.5, "N2(g)": 13.16},
        "products": {"CO2(g)": 2, "H2O(g)": 3, "N2(g)": 13.16}
    },
    {
        "name": "Combustion of Butane",
        "fuel": "Butane",
        "equation": "C4H10(g) + 6.5O2(g) + 24.44N2(g) → 4CO2(g) + 5H2O(g) + 24.44N2(g)",
        "reactants": {"C4H10(g)": 1, "O2(g)": 6.5, "N2(g)": 24.44},
        "products": {"CO2(g)": 4, "H2O(g)": 5, "N2(g)": 24.44}
    },
    {
        "name": "Combustion of Octane",
        "fuel": "Octane",
        "equation": "C8H18(g) + 12.5O2(g) + 47.0N2(g) → 8CO2(g) + 9H2O(g) + 47.0N2(g)",
        "reactants": {"C8H18(g)": 1, "O2(g)": 12.5, "N2(g)": 47.0},
        "products": {"CO2(g)": 8, "H2O(g)": 9, "N2(g)": 47.0}
    },
    {
        "name": "Combustion of Ammonia",
        "fuel": "Ammonia", 
        "equation": "4NH3(g) + 3O2(g) + 11.28N2(g) → 2N2(g) + 6H2O(g) + 11.28N2(g)",
        "reactants": {"NH3(g)": 4, "O2(g)": 3, "N2(g)": 11.28},
        "products": {"N2(g)": 13.28, "H2O(g)": 6}  # Total N2 = 2 + 11.28
    },
    {
        "name": "Combustion of Methanol",
        "fuel": "Methanol",
        "equation": "CH3OH(g) + 1.5O2(g) + 5.64N2(g) → CO2(g) + 2H2O(g) + 5.64N2(g)",
        "reactants": {"CH3OH(g)": 1, "O2(g)": 1.5, "N2(g)": 5.64},
        "products": {"CO2(g)": 1, "H2O(g)": 2, "N2(g)": 5.64}
    },
    {
        "name": "Combustion of Ethanol",
        "fuel": "Ethanol",
        "equation": "C2H5OH(g) + 3O2(g) + 11.28N2(g) → 2CO2(g) + 3H2O(g) + 11.28N2(g)",
        "reactants": {"C2H5OH(g)": 1, "O2(g)": 3, "N2(g)": 11.28},
        "products": {"CO2(g)": 2, "H2O(g)": 3, "N2(g)": 11.28}
    }
]
