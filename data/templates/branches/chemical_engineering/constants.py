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
