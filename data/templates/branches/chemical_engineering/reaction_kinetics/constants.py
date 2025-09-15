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
