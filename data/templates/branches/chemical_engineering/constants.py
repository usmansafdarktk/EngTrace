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


# A dictionary of real-world saturation properties for common thermodynamic substances.
# Each entry contains: temperature (°C), specific volume of saturated liquid (m³/kg),
# and specific volume of saturated vapor (m³/kg)
# Data sources: NIST WebBook, Engineering Toolbox, and standard thermodynamic tables.
# @kind: property
# @units: temp_C=degC, v_f=m^3/kg, v_g=m^3/kg
# @domain: none (each row states its own saturation temperature in temp_C)
# C3: v_f and v_g re-derived from NIST saturation rows requested ON each row's temp_C
#   (fetch_references.py NIST_SATURATION_POINTS); each field's tag names its file and row.
REAL_FLUID_DATA = {
    #  Classic Working Fluids 
    # [ON-DISK] nist_fluid_properties/water_C7732185_saturation_373.15K.tsv @ T=373.15 col="Volume (l, m3/kg)" field=['v_f'] precision=4sf
    #   NIST 0.0010435 at 373.15 K (0.10142 MPa). Corrected at C3 from 0.001043, which was -0.05% from it.
    # [ON-DISK] nist_fluid_properties/water_C7732185_saturation_373.15K.tsv @ T=373.15 col="Volume (v, m3/kg)" field=['v_g'] precision=5sf
    #   NIST 1.6718 at 373.15 K (0.10142 MPa). Corrected at C3 from 1.6729, which was +0.07% from it.
    "Water": {"temp_C": 100, "v_f": 0.001044, "v_g": 1.6718},
    # [ON-DISK] nist_fluid_properties/ammonia_C7664417_saturation_298.15K.tsv @ T=298.15 col="Volume (l, m3/kg)" field=['v_f'] precision=4sf
    #   NIST 0.0016585 at 298.15 K (1.0027 MPa). Corrected at C3 from 0.001658, which was -0.03% from it.
    # [ON-DISK] nist_fluid_properties/ammonia_C7664417_saturation_298.15K.tsv @ T=298.15 col="Volume (v, m3/kg)" field=['v_g'] precision=4sf
    #   NIST 0.12819 at 298.15 K (1.0027 MPa). Corrected at C3 from 0.1284, which was +0.16% from it.
    "Ammonia": {"temp_C": 25, "v_f": 0.001659, "v_g": 0.1282},
    # [ON-DISK] nist_fluid_properties/carbon_dioxide_C124389_saturation_293.15K.tsv @ T=293.15 col="Volume (l, m3/kg)" field=['v_f'] precision=4sf
    #   NIST 0.0012930 at 293.15 K (5.7291 MPa). Corrected at C3 from 0.001286, which was -0.54% from it.
    # [ON-DISK] nist_fluid_properties/carbon_dioxide_C124389_saturation_293.15K.tsv @ T=293.15 col="Volume (v, m3/kg)" field=['v_g'] precision=3sf
    #   NIST 0.0051493 at 293.15 K (5.7291 MPa). Corrected at C3 from 0.0188, which was +265.10% from it.
    "Carbon Dioxide": {"temp_C": 20, "v_f": 0.001293, "v_g": 0.00515},
    # [ON-DISK] nist_fluid_properties/sulfur_dioxide_C7446095_saturation_298.15K.tsv @ T=298.15 col="Volume (l, m3/kg)" field=['v_f'] precision=3sf
    #   NIST 0.00073162 at 298.15 K (0.39215 MPa). Corrected at C3 from 0.000728, which was -0.49% from it.
    # [ON-DISK] nist_fluid_properties/sulfur_dioxide_C7446095_saturation_298.15K.tsv @ T=298.15 col="Volume (v, m3/kg)" field=['v_g'] precision=4sf
    #   NIST 0.091860 at 298.15 K (0.39215 MPa). Corrected at C3 from 0.1165, which was +26.82% from it.
    "Sulfur Dioxide": {"temp_C": 25, "v_f": 0.000732, "v_g": 0.09186},
    
    #  Hydrocarbons 
    # [ON-DISK] nist_fluid_properties/methane_C74828_saturation_112.15K.tsv @ T=112.15 col="Volume (l, m3/kg)" field=['v_f'] precision=4sf
    #   NIST 0.0023716 at 112.15 K (0.10542 MPa). Corrected at C3 from 0.002367, which was -0.19% from it.
    # [ON-DISK] nist_fluid_properties/methane_C74828_saturation_112.15K.tsv @ T=112.15 col="Volume (v, m3/kg)" field=['v_g'] precision=5sf
    #   NIST 0.53089 at 112.15 K (0.10542 MPa). Corrected at C3 from 1.8160, which was +242.07% from it.
    "Methane": {"temp_C": -161, "v_f": 0.002372, "v_g": 0.53089},  # Boiling Point
    # [ON-DISK] nist_fluid_properties/ethane_C74840_saturation_184.15K.tsv @ T=184.15 col="Volume (l, m3/kg)" field=['v_f'] precision=4sf
    #   NIST 0.0018370 at 184.15 K (0.099056 MPa). Corrected at C3 from 0.001830, which was -0.38% from it.
    # [ON-DISK] nist_fluid_properties/ethane_C74840_saturation_184.15K.tsv @ T=184.15 col="Volume (v, m3/kg)" field=['v_g'] precision=4sf
    #   NIST 0.49706 at 184.15 K (0.099056 MPa). Corrected at C3 from 0.7090, which was +42.64% from it.
    "Ethane": {"temp_C": -89, "v_f": 0.001837, "v_g": 0.4971},  # Boiling Point
    # [ON-DISK] nist_fluid_properties/propane_C74986_saturation_298.15K.tsv @ T=298.15 col="Volume (l, m3/kg)" field=['v_f'] precision=4sf
    #   NIST 0.0020310 at 298.15 K (0.95207 MPa). Corrected at C3 from 0.002028, which was -0.15% from it.
    # [ON-DISK] nist_fluid_properties/propane_C74986_saturation_298.15K.tsv @ T=298.15 col="Volume (v, m3/kg)" field=['v_g'] precision=3sf
    #   NIST 0.048502 at 298.15 K (0.95207 MPa). Corrected at C3 from 0.0461, which was -4.95% from it.
    "Propane": {"temp_C": 25, "v_f": 0.002031, "v_g": 0.0485},
    # [ON-DISK] nist_fluid_properties/butane_C106978_saturation_298.15K.tsv @ T=298.15 col="Volume (l, m3/kg)" field=['v_f'] precision=4sf
    #   NIST 0.0017457 at 298.15 K (0.24329 MPa). Corrected at C3 from 0.001726, which was -1.13% from it.
    # [ON-DISK] nist_fluid_properties/butane_C106978_saturation_298.15K.tsv @ T=298.15 col="Volume (v, m3/kg)" field=['v_g'] precision=4sf
    #   NIST 0.16202 at 298.15 K (0.24329 MPa). Corrected at C3 from 0.1632, which was +0.73% from it.
    "Butane": {"temp_C": 25, "v_f": 0.001746, "v_g": 0.1620},
    # [ON-DISK] nist_fluid_properties/isobutane_C75285_saturation_298.15K.tsv @ T=298.15 col="Volume (l, m3/kg)" field=['v_f'] precision=4sf
    #   NIST 0.0018160 at 298.15 K (0.35067 MPa). Corrected at C3 from 0.001815, which was -0.06% from it.
    # [ON-DISK] nist_fluid_properties/isobutane_C75285_saturation_298.15K.tsv @ T=298.15 col="Volume (v, m3/kg)" field=['v_g'] precision=4sf
    #   NIST 0.10958 at 298.15 K (0.35067 MPa). Corrected at C3 from 0.1118, which was +2.03% from it.
    "Isobutane": {"temp_C": 25, "v_f": 0.001816, "v_g": 0.1096},
    # [ON-DISK] nist_fluid_properties/pentane_C109660_saturation_298.15K.tsv @ T=298.15 col="Volume (l, m3/kg)" field=['v_f'] precision=4sf
    #   NIST 0.0016099 at 298.15 K (0.068355 MPa). Corrected at C3 from 0.001598, which was -0.74% from it.
    # [ON-DISK] nist_fluid_properties/pentane_C109660_saturation_298.15K.tsv @ T=298.15 col="Volume (v, m3/kg)" field=['v_g'] precision=4sf
    #   NIST 0.48595 at 298.15 K (0.068355 MPa). Corrected at C3 from 0.5050, which was +3.92% from it.
    "Pentane": {"temp_C": 25, "v_f": 0.001610, "v_g": 0.4860},
    # [ON-DISK] nist_fluid_properties/2_methylbutane_C78784_saturation_298.15K.tsv @ T=298.15 col="Volume (l, m3/kg)" field=['v_f'] precision=4sf
    #   NIST 0.0016261 at 298.15 K (0.091795 MPa). Corrected at C3 from 0.001614, which was -0.74% from it.
    # [ON-DISK] nist_fluid_properties/2_methylbutane_C78784_saturation_298.15K.tsv @ T=298.15 col="Volume (v, m3/kg)" field=['v_g'] precision=4sf
    #   NIST 0.35890 at 298.15 K (0.091795 MPa). Corrected at C3 from 0.4910, which was +36.81% from it.
    #   Iso-pentane is read as NIST's 2-methylbutane: the WebBook name search for "isopentane"
    #   returns "Butane, 2-methyl-", CAS 78-78-4, which is this file's CAS (C78784).
    "Iso-pentane": {"temp_C": 25, "v_f": 0.001626, "v_g": 0.3589},
    
    #  Refrigerants 
    # [ON-DISK] nist_fluid_properties/r11_C75694_saturation_298.15K.tsv @ T=298.15 col="Volume (l, m3/kg)" field=['v_f'] precision=3sf
    #   NIST 0.00067738 at 298.15 K (0.10605 MPa). Corrected at C3 from 0.000675, which was -0.35% from it.
    # [ON-DISK] nist_fluid_properties/r11_C75694_saturation_298.15K.tsv @ T=298.15 col="Volume (v, m3/kg)" field=['v_g'] precision=4sf
    #   NIST 0.16377 at 298.15 K (0.10605 MPa). Corrected at C3 from 0.1580, which was -3.52% from it.
    "Refrigerant-11 (R-11, Trichlorofluoromethane)": {"temp_C": 25, "v_f": 0.000677, "v_g": 0.1638},
    # [ON-DISK] nist_fluid_properties/r12_C75718_saturation_298.15K.tsv @ T=298.15 col="Volume (l, m3/kg)" field=['v_f'] precision=3sf
    #   NIST 0.00076276 at 298.15 K (0.65057 MPa).
    # [ON-DISK] nist_fluid_properties/r12_C75718_saturation_298.15K.tsv @ T=298.15 col="Volume (v, m3/kg)" field=['v_g'] precision=3sf
    #   NIST 0.027153 at 298.15 K (0.65057 MPa). Corrected at C3 from 0.0268, which was -1.30% from it.
    "Refrigerant-12 (R-12, Dichlorodifluoromethane)": {"temp_C": 25, "v_f": 0.000763, "v_g": 0.0272},
    # [ON-DISK] nist_fluid_properties/r22_C75456_saturation_298.15K.tsv @ T=298.15 col="Volume (l, m3/kg)" field=['v_f'] precision=3sf
    #   NIST 0.00083987 at 298.15 K (1.0439 MPa). Corrected at C3 from 0.000845, which was +0.61% from it.
    # [ON-DISK] nist_fluid_properties/r22_C75456_saturation_298.15K.tsv @ T=298.15 col="Volume (v, m3/kg)" field=['v_g'] precision=3sf
    #   NIST 0.022608 at 298.15 K (1.0439 MPa).
    "Refrigerant-22 (R-22, Chlorodifluoromethane)": {"temp_C": 25, "v_f": 0.000840, "v_g": 0.0226},
    # [ON-DISK] nist_fluid_properties/r134a_C811972_saturation_298.15K.tsv @ T=298.15 col="Volume (l, m3/kg)" field=['v_f'] precision=3sf
    #   NIST 0.00082870 at 298.15 K (0.66538 MPa).
    # [ON-DISK] nist_fluid_properties/r134a_C811972_saturation_298.15K.tsv @ T=298.15 col="Volume (v, m3/kg)" field=['v_g'] precision=3sf
    #   NIST 0.030912 at 298.15 K (0.66538 MPa). Corrected at C3 from 0.0323, which was +4.49% from it.
    "Refrigerant-134a (R-134a, 1,1,1,2-Tetrafluoroethane)": {"temp_C": 25, "v_f": 0.000829, "v_g": 0.0309},
    # [ON-DISK] nist_fluid_properties/r123_C306832_saturation_298.15K.tsv @ T=298.15 col="Volume (l, m3/kg)" field=['v_f'] precision=3sf
    #   NIST 0.00068311 at 298.15 K (0.091354 MPa). Corrected at C3 from 0.000680, which was -0.46% from it.
    # [ON-DISK] nist_fluid_properties/r123_C306832_saturation_298.15K.tsv @ T=298.15 col="Volume (v, m3/kg)" field=['v_g'] precision=4sf
    #   NIST 0.17032 at 298.15 K (0.091354 MPa). Corrected at C3 from 0.1810, which was +6.27% from it.
    "Refrigerant-123 (R-123, Dichlorotrifluoroethane)": {"temp_C": 25, "v_f": 0.000683, "v_g": 0.1703},
    # [UNVERIFIED] a BLEND, and neither it nor its components are on disk. The NIST fluid
    #   database here holds 39 PURE fluids (it held 36 before this phase added nonane,
    #   decane and dodecane), and R-410A is not among them - nor is R-32 (difluoromethane)
    #   or R-125 (pentafluoroethane), so the blend cannot even be reconstructed from its
    #   components. The only refrigerants on disk are R-11, R-12, R-22, R-123 and R-134a.
    #   NOT SUBSTITUTED, for two independent reasons. R-22 is the nearest on-disk
    #   refrigerant - at 298.15 K it is P=1.0439 MPa, v_f=0.00083987, v_g=0.022608
    #   (nist_fluid_properties/r22_C75456_saturation_298.15K.tsv) - and this row is
    #   +9.78% from its v_f and -33.21% from its v_g. Substituting would replace the
    #   row's VALUES, not relabel it, which is changing what the row says rather than
    #   sourcing it. And REAL_FLUID_DATA already carries an R-22 row, so it would put one
    #   fluid in the table twice under two keys - the defect just recorded for
    #   Tetrabromoethane / Acetylene Tetrabromide.
    #   The committed saturation volumes are therefore unsourced: no artefact on disk
    #   states them or contradicts them. C3.5 substance defect (D-033).
    "Refrigerant-410A (R-410A, blend of difluoromethane and pentafluoroethane)": {"temp_C": 25, "v_f": 0.000922, "v_g": 0.0151},
    
    #  Organic Solvents (Rankine Fluids) 
    # [ON-DISK] nist_fluid_properties/toluene_C108883_saturation_384.15K.tsv @ T=384.15 col="Volume (l, m3/kg)" field=['v_f'] precision=4sf
    #   NIST 0.0012841 at 384.15 K (0.10249 MPa). Corrected at C3 from 0.001308, which was +1.86% from it.
    # [ON-DISK] nist_fluid_properties/toluene_C108883_saturation_384.15K.tsv @ T=384.15 col="Volume (v, m3/kg)" field=['v_g'] precision=4sf
    #   NIST 0.32441 at 384.15 K (0.10249 MPa). Corrected at C3 from 0.3200, which was -1.36% from it.
    "Toluene": {"temp_C": 111, "v_f": 0.001284, "v_g": 0.3244},  # Boiling Point
    # [ON-DISK] nist_fluid_properties/benzene_C71432_saturation_353.15K.tsv @ T=353.15 col="Volume (l, m3/kg)" field=['v_f'] precision=4sf
    #   NIST 0.0012285 at 353.15 K (0.10111 MPa). Corrected at C3 from 0.001205, which was -1.91% from it.
    # [ON-DISK] nist_fluid_properties/benzene_C71432_saturation_353.15K.tsv @ T=353.15 col="Volume (v, m3/kg)" field=['v_g'] precision=4sf
    #   NIST 0.35918 at 353.15 K (0.10111 MPa). Corrected at C3 from 0.3660, which was +1.90% from it.
    "Benzene": {"temp_C": 80, "v_f": 0.001229, "v_g": 0.3592},  # Boiling Point
    # [UNVERIFIED] not in the NIST fluid database on disk, and its WebBook species pages on disk carry no saturated volumes. Residual register (C3.5).
    "Ethanol": {"temp_C": 78, "v_f": 0.001290, "v_g": 0.5770},  # Boiling Point
    # [ON-DISK] nist_fluid_properties/methanol_C67561_saturation_338.15K.tsv @ T=338.15 col="Volume (l, m3/kg)" field=['v_f'] precision=4sf
    #   NIST 0.0013372 at 338.15 K (0.10341 MPa). Corrected at C3 from 0.001375, which was +2.83% from it.
    # [ON-DISK] nist_fluid_properties/methanol_C67561_saturation_338.15K.tsv @ T=338.15 col="Volume (v, m3/kg)" field=['v_g'] precision=4sf
    #   NIST 0.80347 at 338.15 K (0.10341 MPa). Corrected at C3 from 0.8800, which was +9.52% from it.
    "Methanol": {"temp_C": 65, "v_f": 0.001337, "v_g": 0.8035},  # Boiling Point
    # [UNVERIFIED] not in the NIST fluid database on disk, and its WebBook species pages on disk carry no saturated volumes. Residual register (C3.5).
    "Acetone": {"temp_C": 56, "v_f": 0.001360, "v_g": 0.5370},  # Boiling Point
    # [ON-DISK] nist_fluid_properties/hexane_C110543_saturation_342.15K.tsv @ T=342.15 col="Volume (l, m3/kg)" field=['v_f'] precision=4sf
    #   NIST 0.0016311 at 342.15 K (0.10223 MPa). Corrected at C3 from 0.001660, which was +1.77% from it.
    # [ON-DISK] nist_fluid_properties/hexane_C110543_saturation_342.15K.tsv @ T=342.15 col="Volume (v, m3/kg)" field=['v_g'] precision=4sf
    #   NIST 0.30704 at 342.15 K (0.10223 MPa). Corrected at C3 from 0.3550, which was +15.62% from it.
    "n-Hexane": {"temp_C": 69, "v_f": 0.001631, "v_g": 0.3070},  # Boiling Point
    # [ON-DISK] nist_fluid_properties/octane_C111659_saturation_399.15K.tsv @ T=399.15 col="Volume (l, m3/kg)" field=['v_f'] precision=4sf
    #   NIST 0.0016370 at 399.15 K (0.10233 MPa). Corrected at C3 from 0.001690, which was +3.24% from it.
    # [ON-DISK] nist_fluid_properties/octane_C111659_saturation_399.15K.tsv @ T=399.15 col="Volume (v, m3/kg)" field=['v_g'] precision=4sf
    #   NIST 0.26783 at 399.15 K (0.10233 MPa). Corrected at C3 from 0.2500, which was -6.66% from it.
    "n-Octane": {"temp_C": 126, "v_f": 0.001637, "v_g": 0.2678},  # Boiling Point
    # [ON-DISK] nist_fluid_properties/cyclohexane_C110827_saturation_354.15K.tsv @ T=354.15 col="Volume (l, m3/kg)" field=['v_f'] precision=4sf
    #   NIST 0.0013905 at 354.15 K (0.10220 MPa). Corrected at C3 from 0.001420, which was +2.12% from it.
    # [ON-DISK] nist_fluid_properties/cyclohexane_C110827_saturation_354.15K.tsv @ T=354.15 col="Volume (v, m3/kg)" field=['v_g'] precision=4sf
    #   NIST 0.32946 at 354.15 K (0.10220 MPa). Corrected at C3 from 0.3600, which was +9.27% from it.
    "Cyclohexane": {"temp_C": 81, "v_f": 0.001391, "v_g": 0.3295}  # Boiling Point
}


# A dictionary with comprehensive critical properties for various substances.
# Tc: Kelvin (K), Pc: bar, Vc: cm³/mol, Zc: dimensionless, omega: dimensionless.
# @kind: property
# @units: Tc=K, Pc=bar, Vc=cm^3/mol, Zc=1, omega=1
# @domain: none (the critical point is itself the condition)
# [UNVERIFIED] no on-disk artefact covers this 27-row set. The WebBook species pages in
#   docs/references/nist_webbook_species are the only critical-point source on disk and
#   they cover 3 of the 27 rows, tagged individually below; the rest are asserted.
#   Registered as a C3.5 residual.
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
    # [ON-DISK] nist_webbook_species/p_xylene_C106423_phase_change.html @ text="T c (K)"
    #   The page's phase-change table gives Tc = 616.2 K (Majer and Svoboda, 1985),
    #   which is this row's Tc. Locator-only: the value sits in a positional HTML
    #   table the resolver does not parse, so it pins the page, not the number.
    "p-Xylene": {"Tc": 616.2, "Pc": 35.12, "Vc": 379.0, "Zc": 0.260, "omega": 0.321},
    "Methanol": {"Tc": 512.6, "Pc": 80.84, "Vc": 118.0, "Zc": 0.224, "omega": 0.564},
    # [ON-DISK] nist_webbook_species/ethanol_C64175_phase_change.html @ text="T c (K)"
    #   The page's phase-change table gives Tc = 513.9 K (Majer and Svoboda, 1985),
    #   which is this row's Tc. Locator-only: the value sits in a positional HTML
    #   table the resolver does not parse, so it pins the page, not the number.
    "Ethanol": {"Tc": 513.9, "Pc": 61.37, "Vc": 167.0, "Zc": 0.248, "omega": 0.645},
    # [ON-DISK] nist_webbook_species/acetone_C67641_phase_change.html @ text="T c (K)"
    #   The page's phase-change table gives Tc = 508.2 K (Majer and Svoboda, 1985),
    #   which is this row's Tc. Locator-only: the value sits in a positional HTML
    #   table the resolver does not parse, so it pins the page, not the number.
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
# Added 'min_temp' and 'max_temp' (in °C) to ensure phase stability.
# @kind: property
# @units: Cp=J/(g*K), min_temp=degC, max_temp=degC
# @domain: none (each row states its own min_temp/max_temp window)
# [UNVERIFIED] specific heats in J/(g*K) with no temperature stated per value. Nothing on
#   disk carries this table: the WebBook pages hold molar Cp for 5 of the 35 rows, in
#   different units and at different conditions. C3.5 residual.
SUBSTANCES_FOR_HEATING = [
    # Metals & Solids (Generally safe 20°C - 500°C)
    {"name": "Iron", "state": "solid", "Cp": 0.449, "min_temp": 20, "max_temp": 500},
    {"name": "Copper", "state": "solid", "Cp": 0.385, "min_temp": 20, "max_temp": 500},
    {"name": "Aluminum", "state": "solid", "Cp": 0.897, "min_temp": 20, "max_temp": 500},
    {"name": "Gold", "state": "solid", "Cp": 0.129, "min_temp": 20, "max_temp": 500},
    {"name": "Lead", "state": "solid", "Cp": 0.16, "min_temp": 20, "max_temp": 300},
    {"name": "Silver", "state": "solid", "Cp": 0.235, "min_temp": 20, "max_temp": 500},
    {"name": "Tungsten", "state": "solid", "Cp": 0.134, "min_temp": 20, "max_temp": 1000},
    {"name": "Silicon", "state": "solid", "Cp": 0.705, "min_temp": 20, "max_temp": 500},
    {"name": "Graphite (Carbon)", "state": "solid", "Cp": 0.709, "min_temp": 20, "max_temp": 1000},
    {"name": "Glass (typical)", "state": "solid", "Cp": 0.84, "min_temp": 20, "max_temp": 500},
    {"name": "Concrete", "state": "solid", "Cp": 0.88, "min_temp": 20, "max_temp": 500},
    {"name": "Wood (typical)", "state": "solid", "Cp": 1.7, "min_temp": 20, "max_temp": 150},
    {"name": "Polyethylene (plastic)", "state": "solid", "Cp": 2.3, "min_temp": 20, "max_temp": 100}, # Melts ~115°C

    # Phase-Sensitive Solids
    {"name": "Ice (at 0°C)", "state": "solid", "Cp": 2.09, "min_temp": -50, "max_temp": -2}, # Must be < 0°C

    # Liquids (Must stay between Freezing and Boiling Points)
    {"name": "Water", "state": "liquid", "Cp": 4.18, "min_temp": 5, "max_temp": 95},
    {"name": "Ethanol", "state": "liquid", "Cp": 2.44, "min_temp": -50, "max_temp": 75}, # Boils at 78°C
    {"name": "Methanol", "state": "liquid", "Cp": 2.53, "min_temp": -50, "max_temp": 60}, # Boils at 65°C
    {"name": "Acetone", "state": "liquid", "Cp": 2.17, "min_temp": -50, "max_temp": 50}, # Boils at 56°C
    {"name": "Mercury", "state": "liquid", "Cp": 0.14, "min_temp": -30, "max_temp": 300},
    {"name": "Glycerol", "state": "liquid", "Cp": 2.43, "min_temp": 20, "max_temp": 250},
    {"name": "Ethylene Glycol (Antifreeze)", "state": "liquid", "Cp": 2.36, "min_temp": -10, "max_temp": 190},
    {"name": "Olive Oil", "state": "liquid", "Cp": 1.97, "min_temp": 20, "max_temp": 200},
    {"name": "Engine Oil (typical)", "state": "liquid", "Cp": 1.9, "min_temp": 20, "max_temp": 250},
    {"name": "Sulfuric Acid", "state": "liquid", "Cp": 1.42, "min_temp": 20, "max_temp": 300},

    # Gases (Wide range, but avoid condensation for vapors)
    {"name": "Air (dry)", "state": "gas", "Cp": 1.005, "min_temp": -50, "max_temp": 500},
    {"name": "Nitrogen", "state": "gas", "Cp": 1.04, "min_temp": -100, "max_temp": 500},
    {"name": "Oxygen", "state": "gas", "Cp": 0.918, "min_temp": -100, "max_temp": 500},
    {"name": "Hydrogen", "state": "gas", "Cp": 14.31, "min_temp": -100, "max_temp": 500},
    {"name": "Helium", "state": "gas", "Cp": 5.193, "min_temp": -200, "max_temp": 500},
    {"name": "Argon", "state": "gas", "Cp": 0.520, "min_temp": -100, "max_temp": 500},
    {"name": "Carbon Dioxide", "state": "gas", "Cp": 0.839, "min_temp": -50, "max_temp": 500},
    {"name": "Methane", "state": "gas", "Cp": 2.22, "min_temp": -100, "max_temp": 500},
    {"name": "Ammonia", "state": "gas", "Cp": 2.06, "min_temp": -20, "max_temp": 400},
    {"name": "Water Vapor (Steam)", "state": "gas", "Cp": 2.01, "min_temp": 105, "max_temp": 400}, # Must be > 100°C
    {"name": "Chlorine", "state": "gas", "Cp": 0.48, "min_temp": 20, "max_temp": 200}
]


# A list of common substances with their molar heats of vaporization (delta_H_vap)
# at their normal boiling points. All values are in kJ/mol.
# @kind: property
# @units: delta_H_vap=kJ/mol
# @domain: none (each value is at its own substance's normal boiling point)
# [UNVERIFIED] as a table: 4 of the 16 rows resolve to a WebBook phase-change page and are
#   tagged below. The other 12 do not, and they are NOT [DERIVED] either, though they
#   could nearly be: dHvap = Enthalpy(v) - Enthalpy(l) at 1 atm from the NIST saturation
#   tables reproduces EVERY row that has a saturation file - 11 of 11 - to within
#   1.454%, eight of them within 0.2% (water -0.014, methanol +0.048, propane +1.454,
#   n-butane +0.095, n-hexane -0.108, benzene +0.153, toluene +0.737, ammonia -0.115,
#   nitrogen -0.174, oxygen +0.036, argon -0.112 %). Four are roundings of it.
#   CORRECTED at the C3 review (Reviewer H, chemical, F-2): this note previously said
#   "7 of them" and justified the class by "the rows that miss". NO ROW MISSES.
#   Mercury is the twelfth row and has no saturation file at all, so it is not a row
#   that fails the derivation - it is one the derivation cannot reach. The class stays
#   unverified because a ~1% agreement is evidence about provenance and not a
#   derivation (D-077), but the reason that was written for it was false.
SUBSTANCES_FOR_VAPORIZATION = [
    # Alcohols & Water
    {"name": "Water", "delta_H_vap": 40.66},
    {"name": "Methanol", "delta_H_vap": 35.3},
    # [ON-DISK] nist_webbook_species/ethanol_C64175_phase_change.html @ text="Enthalpy of vaporization"
    #   The page's first entry is 38.56 kJ/mol at 351.5 K (Majer and Svoboda, 1985).
    {"name": "Ethanol", "delta_H_vap": 38.6},
    # [ON-DISK] nist_webbook_species/isopropanol_C67630_phase_change.html @ text="Enthalpy of vaporization"
    #   The page's first entry is 39.85 kJ/mol at 355.4 K (Majer and Svoboda, 1985).
    {"name": "Isopropanol", "delta_H_vap": 39.85},

    # Alkanes
    {"name": "Propane", "delta_H_vap": 19.04},
    {"name": "n-Butane", "delta_H_vap": 22.44},
    {"name": "n-Hexane", "delta_H_vap": 28.85},

    # Organic Solvents
    # [ON-DISK] nist_webbook_species/acetone_C67641_phase_change.html @ text="Enthalpy of vaporization"
    #   The page's first entry is 29.1 kJ/mol at 329.3 K (Majer and Svoboda, 1985).
    {"name": "Acetone", "delta_H_vap": 29.1},
    {"name": "Benzene", "delta_H_vap": 30.8},
    {"name": "Toluene", "delta_H_vap": 33.48},
    # [ON-DISK] nist_webbook_species/carbon_tetrachloride_C56235_phase_change.html @ text="Enthalpy of vaporization"
    #   The page's first entry is 29.82 kJ/mol at 349.9 K (Majer and Svoboda, 1985).
    {"name": "Carbon Tetrachloride", "delta_H_vap": 29.82},

    # Inorganic & Elemental Substances
    {"name": "Ammonia", "delta_H_vap": 23.3},
    # [UNVERIFIED] no artefact on disk states a dHvap for mercury, and none can be derived.
    #   nist_webbook_species/mercury_C7439976_phase_change.html carries T boil (629.81 K,
    #   Marsh 1987), T fus, T triple and one Antoine row - and the strings "Enthalpy",
    #   "vaporization", "vapH" and "Hvap" do not occur anywhere in the file. Verified by
    #   reading it, not inferred from the section headings.
    #   The dHvap = H(v) - H(l) derivation that reaches the other eleven rows of this table
    #   cannot reach this one either: there is no mercury file in nist_fluid_properties at
    #   all, among 39 pure fluids. So 59.11 kJ/mol is neither supported nor
    #   contradicted - it is unsourced, and no further WebBook fetch would change that.
    #   C3.5 residual.
    {"name": "Mercury", "delta_H_vap": 59.11},
    {"name": "Nitrogen", "delta_H_vap": 5.57},
    {"name": "Oxygen", "delta_H_vap": 6.82},
    {"name": "Argon", "delta_H_vap": 6.43},
]


# Database of standard heats of formation (ΔH_f°) at 298.15 K in kJ/mol.
#
# UNITS. Standard enthalpy of formation at 298.15 K and 1 bar, in kJ/mol, for
# the species in the phase its key names. Elements in their standard state are
# exactly 0 by definition, not by measurement.
#
# PROVENANCE (Phase C2). Every value resolves to a specific measurement in
# docs/references/nist_webbook/shomate_coefficients.json, retrieved from the
# NIST Chemistry WebBook. Each citation names the CAS number, the value, its
# uncertainty and the MEASUREMENT the row uses.
#
# WHY THE MEASUREMENT AND NOT JUST THE SPECIES. NIST lists several independent
# determinations for most hydrocarbons, differing by more than any one of them
# claims. An earlier version of this header said all 21 were "inside their
# stated uncertainty", while citing NIST's recommended value for each row and
# using a different one - so three rows (C2H6, C3H8, CH3OH(l)) sat outside the
# uncertainty they themselves quoted, and the suite passed only because its
# tolerance was floored at 2 kJ/mol. The values were right; the citations named
# the wrong measurement (Phase C2 Reviewer G, findings G-2 and G-3).
#
# This table follows PROSEN AND ROSSINI (1945) for the hydrocarbons wherever
# NIST lists it - a self-consistent set from one laboratory - which is why
# C2H6, C3H8 and C8H18 do not carry NIST's recommended value. That is a source
# choice, and it is now written down instead of inferred.
#
# Reviewer H finding F-3: the evidence existed only in the test suite and the
# JSON, so the table did not carry its own provenance. It does now.
# @kind: property
# @units: kJ/mol
# @domain: T=298.15 K, p=1 bar
HEATS_OF_FORMATION = {
    # Hydrocarbons (Gases)
    # [ON-DISK] NIST 74-82-8, -74.6+/-0.3 (Manion 2002, adopting Gurvich 1991)
    #   C3 review, Reviewer H (chemical) F-5: the row uses -74.8, which is neither the
    #   cited -74.6 nor the JSON's alternative -74.87 (Chase 1998). It lies inside the
    #   cited +/-0.3, so this is not the C2 defect class, but the tag names a
    #   measurement the value is not. Referred (C3.5).
    "CH4(g)": -74.8,      # Methane
    # [ON-DISK] NIST 74-86-2, 226.73 (Chase 1998)
    "C2H2(g)": 226.7,     # Acetylene (Added)
    # [ON-DISK] NIST 74-84-0, -84.67+/-0.49 (Prosen and Rossini 1945).
    #   NIST's recommended value is -84.0+/-0.4 (Manion 2002); this row uses
    #   Prosen, which is where -84.7 comes from. Citing the recommendation
    #   while using Prosen is what made this row look 1.75 sigma out (G-2).
    "C2H6(g)": -84.7,     # Ethane
    # [ON-DISK] NIST 74-98-6, -103.8+/-0.59 (Prosen and Rossini 1945).
    #   NIST lists TWO measured values: -104.7+/-0.50 (Pittam and Pilcher 1972,
    #   recommended) and -103.8 (Prosen). This row uses Prosen. A source
    #   choice, not an error - Reviewer H retired the interim doc's "propane
    #   discrepancy" on exactly this point. Both are now on disk (G-3).
    "C3H8(g)": -103.8,    # Propane
    # [ON-DISK] NIST 106-97-8, -125.6+/-0.67
    "C4H10(g)": -125.7,   # Butane (Added)
    # [ON-DISK] NIST 111-65-9, -208.4+/-0.67 (Prosen and Rossini 1945).
    #   NIST also lists -208.7 (Good 1972, computed via the heat of
    #   vaporisation). No Cp fit exists for octane; the entry backs this
    #   citation only.
    "C8H18(g)": -208.4,   # Octane (Gas phase) (Added)
    # [ON-DISK] NIST 71-43-2 condensed phase, 49.0+/-0.9
    "C6H6(l)": 49.0,      # Benzene

    # Alcohols (Liquids & Gases)
    # [ON-DISK] NIST 67-56-1 condensed phase. This row is the one value in the
    #   table that transcribes NO single NIST measurement: -238.6 sits between
    #   Baroody and Carpenter 1972 (-238.4, no uncertainty stated) and Green
    #   1960 (-238.9+/-3.6), and is inside Green's uncertainty. It is NOT
    #   inside Chao and Rossini's -239.5+/-0.2, which this row used to cite
    #   (G-2). Kept, because two NIST-listed measurements bracket it.
    "CH3OH(l)": -238.6,   # Methanol (Liquid)
    # [ON-DISK] NIST 67-56-1 gas, -205+/-10 (NIST average of 9 values).
    #   NIST also lists -201.49+/-0.20 (Rossini 1932, flame calorimetry), which
    #   is the nearer of the two to this row.
    "CH3OH(g)": -200.7,   # Methanol (Gas) (Added - required for adiabatic flame temp)
    # [ON-DISK] NIST 64-17-5 condensed phase, -276+/-2 (NIST average of 6).
    #   The nearest itemised determination is Green 1960 at -277.6.
    "C2H5OH(l)": -277.7,  # Ethanol (Liquid)
    # [ON-DISK] NIST 64-17-5 gas, -234+/-2
    "C2H5OH(g)": -235.1,  # Ethanol (Gas) (Added - required for adiabatic flame temp)

    # Common Gases & Products
    # [BY-DEFINITION] element in its standard state, exactly 0
    "O2(g)": 0,
    # [BY-DEFINITION] element in its standard state, exactly 0
    "H2(g)": 0,
    # [BY-DEFINITION] element in its standard state, exactly 0
    "N2(g)": 0,
    # [ON-DISK] NIST 630-08-0, -110.53
    "CO(g)": -110.5,
    # [ON-DISK] NIST 124-38-9, -393.51+/-0.13 (CODATA)
    "CO2(g)": -393.5,
    # [ON-DISK] NIST 7732-18-5 gas, -241.826+/-0.040 (CODATA)
    "H2O(g)": -241.8,
    # [ON-DISK] NIST 7732-18-5 condensed, -285.83 (CODATA)
    "H2O(l)": -285.8,
    # [ON-DISK] NIST 7664-41-7, -45.94+/-0.35
    "NH3(g)": -46.1,
    # [ON-DISK] NIST 10102-43-9, 90.29
    "NO(g)": 90.3,
    # [ON-DISK] NIST 10102-44-0, 33.10 (Chase 1998) - the ONLY value NIST
    #   lists. This row read 33.2; corrected to 33.1. A 0.1 kJ/mol
    #   transcription slip, found by requiring the row to match a measurement
    #   on disk rather than merely to name a species (G-2).
    "NO2(g)": 33.1,
}

# A list of predefined, balanced chemical reactions.
# @kind: mathematical
# @units: reactants=1, products=1
# @domain: none (stoichiometry; no temperature or pressure condition)
# [DERIVED] every reaction balances. The patch script counts atoms on both sides of all recompute=atom_balance
#   4 reactions from the species keys themselves and refuses to write if any differs.
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
# Heat capacity parameters for the ideal-gas / condensed-phase polynomial
#
#     Cp/R = A + B*T + C*T^2 + D*T^-2       (T in Kelvin)
#
# UNITS. A, B, C and D are DIMENSIONLESS as written: the polynomial returns
# Cp/R, so a caller multiplies by R = 8.314 J/(mol K) to get Cp in J/(mol K).
# B, C and D carry implied units of K^-1, K^-2 and K^2 respectively, which is
# why each exponent is attached to its value rather than to a column heading -
# see the paragraph below.
#
# Valid roughly 298-1500 K for the gases. The coefficients are written with the
# exponent attached to each value - B as E-3, C as E-6, D as E5 - because the
# source tabulates them as the COLUMN HEADINGS 10^3 B, 10^6 C and 10^-5 D, and
# transcribing a column heading as if it were part of the value is exactly the
# defect Phase C2 found: every B was 10x too large and every C was 10x too
# large, which put Cp(N2, 300 K) at 42.4 J/mol/K against a true 29.1.
#
# PROVENANCE. Each row carries an [ON-DISK] citation resolving to
# docs/references/nist_webbook/shomate_coefficients.json, which holds the NIST
# Chemistry WebBook (SRD 69) data the row was checked against, retrieved
# 2026-09-06. NIST publishes the Shomate form - a different fit to the same
# thermochemistry - so agreement between the two is evidence, not a tautology.
# The "Cp298 x vs NIST y" note on each row is that check.
#
# Every row now carries an [ON-DISK] citation. Two species that could not be
# sourced at all - H2SO4(l) and CaCO3(s), both measurably wrong and both absent
# from NIST - were REPLACED by species NIST does cover, rather than shipped
# behind a tag (D-033).
#
# VALIDITY. CP_VALID_T_MAX below is the upper temperature this polynomial is
# fitted for. It is not decoration: the polynomial is a good fit inside it and
# degrades fast outside. Against NIST, CO2 is within 0.7% at 1500 K, +3.6% at
# 2000 K, +8.9% at 2500 K and +13.5% at 2900 K. Any template integrating Cp
# beyond CP_VALID_T_MAX is extrapolating and must say so (Phase C2 Reviewer H,
# finding F-2; DECISIONS D-032).
# @kind: property
# @units: A=1, B=1/K, C=1/K^2, D=K^2
# @domain: T_lo=298..1500 K, T_hi=298..1500 K
#   The per-row ceiling is CP_VALID_T_MAX below; 1500 K is the common case and the
#   widest this polynomial is fitted for.
CP_PARAMS = {
    # Key order is deliberately the original one.
    # template_sensible_heat_temp_dependent_cp draws its substance with
    # random.choice(list(CP_PARAMS.keys())), so the ORDER of this dict is part
    # of the item pool's identity: regrouping it by chemical family changed
    # which substance 274 of 300 seeds drew, for no correctness reason.

    # Hydrocarbons (gases)
    # [ON-DISK] NIST 74-82-8  Cp298 35.06 vs NIST 35.65
    "CH4(g)": {"A": 1.702, "B": 9.081E-3, "C": -2.164E-6, "D": 0.0E5},
    # Combustion products and common gases
    # [ON-DISK] NIST 124-38-9  Cp298 37.14 vs NIST 37.13
    "CO2(g)": {"A": 5.457, "B": 1.045E-3, "C": 0.0E-6, "D": -1.157E5},
    # [ON-DISK] NIST 7727-37-9  Cp298 29.12 vs NIST 29.12
    "N2(g)": {"A": 3.28, "B": 0.593E-3, "C": 0.0E-6, "D": 0.04E5},
    # Liquids and solids
    # [ON-DISK] NIST 7732-18-5  Cp298 75.40 vs NIST 75.3
    "H2O(l)": {"A": 8.712, "B": 1.25E-3, "C": -0.18E-6, "D": 0.0E5},
    # Oxygenates (gases)
    # [ON-DISK] NIST 64-17-5 (refit)  REFIT to the TRC 1997 Cp table; the previous row gave Cp298 74.40 against 65.21
    "C2H5OH(g)": {"A": 1.9799, "B": 22.1928E-3, "C": -6.8999E-6, "D": 0.0E5},
    # [ON-DISK] NIST 7782-44-7  REPLACED: the previous row gave Cp298 34.71 against NIST 29.38
    "O2(g)": {"A": 3.639, "B": 0.506E-3, "C": 0.0E-6, "D": -0.227E5},
    # [ON-DISK] NIST 1333-74-0  Cp298 28.84 vs NIST 28.84
    "H2(g)": {"A": 3.249, "B": 0.422E-3, "C": 0.0E-6, "D": 0.083E5},
    # Inorganic gases
    # [ON-DISK] NIST 7664-41-7  Cp298 35.50 vs NIST 35.65
    "NH3(g)": {"A": 3.578, "B": 3.02E-3, "C": 0.0E-6, "D": -0.186E5},
    # [ON-DISK] NIST 630-08-0  Cp298 29.16 vs NIST 29.15
    "CO(g)": {"A": 3.376, "B": 0.557E-3, "C": 0.0E-6, "D": -0.031E5},
    # [ON-DISK] NIST 7783-06-4  Cp298 34.15 vs NIST 34.19
    "H2S(g)": {"A": 3.931, "B": 1.49E-3, "C": 0.0E-6, "D": -0.232E5},
    # [ON-DISK] NIST 74-84-0  5-point Gurvich 1989 table, 298-1500 K, worst
    #   -2.6% at 1500 K. This row used to record only the Smith-Van Ness
    #   Cp298/R self-check - a check against the source's own column, which is
    #   the tautology this table's header disclaims. The multi-point NIST data
    #   was on disk and unconsulted (G-8).
    "C2H6(g)": {"A": 1.131, "B": 19.225E-3, "C": -5.561E-6, "D": 0.0E5},
    # [ON-DISK] NIST 74-98-6  Cp298 74.92 vs NIST 73.60
    "C3H8(g)": {"A": 1.213, "B": 28.785E-3, "C": -8.824E-6, "D": 0.0E5},
    # [ON-DISK] NIST 106-97-8  Cp298 99.17 vs NIST 98.49
    "C4H10(g)": {"A": 1.935, "B": 36.915E-3, "C": -11.402E-6, "D": 0.0E5},
    # [ON-DISK] NIST 74-85-1  Cp298 44.28 vs NIST 42.90
    "C2H4(g)": {"A": 1.424, "B": 14.394E-3, "C": -4.392E-6, "D": 0.0E5},
    # [ON-DISK] NIST 74-86-2 (refit)  REFIT: the previous row gave Cp298 68.39 against NIST 44.04
    "C2H2(g)": {"A": 3.4669, "B": 7.5416E-3, "C": -2.799E-6, "D": 0.0E5},
    # [ON-DISK] NIST 71-43-2  8-point TRC 1997 GAS table, 298-1500 K, worst
    #   -4.5% at 500 K. Was keyed C6H6(l) but holds gas coefficients; the
    #   re-key was justified by the liquid value (135.69) alone, which is
    #   evidence the row is misfiled, not evidence the row is right (G-8).
    "C6H6(g)": {"A": -0.206, "B": 39.064E-3, "C": -13.301E-6, "D": 0.0E5},
    # [ON-DISK] NIST 108-88-3  8-point Scott 1962 / Draeger 1985 GAS table,
    #   298-1500 K, worst +3.6% at 298 K. Was keyed C7H8(l); as for benzene,
    #   the liquid value (157.09) justified the re-key but verified nothing
    #   about the coefficients themselves (G-8).
    "C7H8(g)": {"A": 0.29, "B": 47.052E-3, "C": -15.716E-6, "D": 0.0E5},
    # [ON-DISK] NIST 67-64-1 (refit)  was keyed C3H6O(l) and held gas values that were
    # themselves 8.4% high; refitted to the Chao 1986 table, worst error 1.11%
    "C3H6O(g)": {"A": 2.4096, "B": 24.7630E-3, "C": -7.5941E-6, "D": 0.0E5},
    # [ON-DISK] NIST 67-56-1  Cp298 80.28 vs NIST 79.5
    "CH3OH(l)": {"A": 5.052, "B": 16.561E-3, "C": -3.761E-6, "D": 0.0E5},
    # [ON-DISK] NIST 110-54-3 (refit)  was keyed C6H14(l) and held values 11.8% below
    # the NIST GAS value; refitted to the Scott 1974 table
    "C6H14(g)": {"A": 3.1704, "B": 53.3808E-3, "C": -16.3758E-6, "D": 0.0E5},
    # [ON-DISK] NIST 7446-09-5  Cp298 39.88 vs NIST 39.87
    "SO2(g)": {"A": 5.699, "B": 0.801E-3, "C": 0.0E-6, "D": -1.015E5},
    # [ON-DISK] NIST 10102-43-9  Cp298 29.81 vs NIST 29.86
    "NO(g)": {"A": 3.387, "B": 0.669E-3, "C": 0.0E-6, "D": 0.095E5},
    # [ON-DISK] NIST 10102-44-0  REPLACED: D's mantissa -0.792 had been written into the C column
    "NO2(g)": {"A": 4.982, "B": 1.195E-3, "C": 0.0E-6, "D": -0.792E5},
    # [ON-DISK] NIST 7782-50-5  Cp298 33.85 vs NIST 33.95
    "Cl2(g)": {"A": 4.442, "B": 0.089E-3, "C": 0.0E-6, "D": -0.344E5},
    # [ON-DISK] NIST 7647-01-0  Cp298 29.14 vs NIST 29.14
    "HCl(g)": {"A": 3.156, "B": 0.623E-3, "C": 0.0E-6, "D": 0.151E5},
    # Monatomic gases - exact, Cp = 5R/2, no source required
    # [DERIVED] exact (monatomic ideal gas) recompute=monatomic_cp
    "He(g)": {"A": 2.5, "B": 0.0E-3, "C": 0.0E-6, "D": 0.0E5},
    # [DERIVED] exact (monatomic ideal gas) recompute=monatomic_cp
    "Ar(g)": {"A": 2.5, "B": 0.0E-3, "C": 0.0E-6, "D": 0.0E5},
    # [DERIVED] exact (monatomic ideal gas) recompute=monatomic_cp
    "Ne(g)": {"A": 2.5, "B": 0.0E-3, "C": 0.0E-6, "D": 0.0E5},
    # [DERIVED] mixture, no NIST entry  verified against the mole-weighted N2/O2/Ar average to 0.36%
    "Air(g)": {"A": 3.355, "B": 0.575E-3, "C": 0.0E-6, "D": -0.016E5},
    # [ON-DISK] NIST 7732-18-5  Cp500 35.28 vs NIST 35.22
    "H2O(g)": {"A": 3.47, "B": 1.45E-3, "C": 0.0E-6, "D": 0.121E5},
    # [ON-DISK] NIST 71-43-2 condensed phase (fit)
    # REPLACES H2SO4(l), which was ~59% low and which NIST carries no
    # condensed-phase Cp for, in either its free or its paid-linked sections.
    # Rather than ship a known-wrong row behind a tag, the species is swapped
    # for one that IS citable (D-033). Fitted to four NIST liquid measurements
    # spanning 293-322 K (Kalali 1987, Grolier 1993, Reddy 1986, Naziev 1986),
    # worst error 0.49%. Cp298 = 135.83 against NIST 135.69. Benzene boils at
    # 353 K, so the template's 280-350 K liquid window stays in phase.
    "C6H6(l)": {"A": 9.9167, "B": 21.5325E-3, "C": 0.0E-6, "D": 0.0E5},
    # [ON-DISK] NIST 7647-14-5 (refit)  the transcribed row read 47.70 against NIST
    # 50.50 (-5.5%). That was waved through as "a source disagreement" with no
    # argument; NIST does publish a solid-phase Shomate, so it is resolvable.
    # Refitted over 298-1073 K, worst error 0.19%.
    "NaCl(s)": {"A": 6.8261, "B": -1.6794E-3, "C": 2.6971E-6, "D": -0.4470E5},
    # [ON-DISK] NIST 1344-28-1 solid phase, corundum (fit)
    # REPLACES CaCO3(s), which was ~9% low and for which NIST's free tier
    # carries no Cp (verified against the calcite and calcium-carbonate pages
    # both). Swapped for a solid NIST covers with a full Shomate over
    # 298-2327 K (D-033). Fitted over 298-1200 K, the template's solid window,
    # worst error 0.06%. Cp298 = 78.76 against NIST 78.80.
    "Al2O3(s)": {"A": 12.5528, "B": 3.9106E-3, "C": -1.0715E-6, "D": -3.6899E5},
}


# Upper validity temperature of the CP_PARAMS polynomial, in Kelvin.
#
# 1500 K is the Smith-Van Ness Tmax for the gases and is the figure the source
# table states. The refitted rows carry the range they were actually fitted
# over. Condensed phases are capped at a temperature below their normal boiling
# or decomposition point.
# VALIDITY CEILING, in Kelvin, per species. ADVISORY, NOT ENFORCED: the only
# consumer is tests/constants_integrity/test_chemical_thermochemistry.py, which
# reports how far a template runs past it. No template reads this table, so
# nothing stops an item generating outside the range - that is the open half of
# Reviewer H's F-2, and the refit-vs-restrict decision is D-032 (Phase 2).
# Values are this repo's own fit ranges, not NIST's; each row says which.
# @kind: validity
# @units: K
# @copied-in: template_sensible_heat_temp_dependent_cp 350.0; template_sensible_heat_temp_dependent_cp 1200.0
#   C3.8: heat_effects.py draws T2 up to 350.0 and 1200.0 - this table's C6H6(l) and Al2O3(s)
#   limits - as literals; a refit that moved them would not reach the template
#   (phaseC3_literal_copies.md).
# @domain: none (this table IS a domain declaration)
# [DERIVED] each row is the top of the temperature interval this repo fitted that
#   species over - stated per row above, and not NIST's own range of validity.
CP_VALID_T_MAX = {
    **{k: 1500.0 for k in CP_PARAMS},
    "C2H5OH(g)": 1500.0,     # refit over 298-1500 K
    "C3H6O(g)": 1500.0,      # refit over 298-1500 K
    "C6H14(g)": 1500.0,      # refit over 298-1500 K
    "C2H2(g)": 1100.0,       # refit over 298-1100 K, NIST's lowest range top
    "NaCl(s)": 1073.0,       # refit over 298-1073 K
    "H2O(l)": 373.0,
    "CH3OH(l)": 337.0,
    "C6H6(l)": 350.0,        # fitted 293-322 K; benzene boils at 353 K
    "Al2O3(s)": 1200.0,      # fitted 298-1200 K; NIST Shomate runs to 2327 K
    "He(g)": 6000.0, "Ar(g)": 6000.0, "Ne(g)": 6000.0,   # exact at any T
}


# Dry-air composition used to verify the Air(g) row against its components.
# Recorded because Reviewer H noted the check was internal and the assumed
# composition was nowhere stated.
# Dry-air mole fractions, dimensionless. Standard dry-atmosphere values (US
# Standard Atmosphere 1976); CO2 and the remaining trace gases (~0.00036 mole
# fraction) are dropped, so these three sum to 0.99964 and a consumer must
# normalise if it needs unity. Used to verify the Air(g) Cp row, which is a
# mixture with no NIST entry of its own. ADVISORY like the table above: the
# test suite is its only consumer.
# @kind: property
# @units: 1
# @domain: none (dry air; the composition is the definition, not a measurement at a T)
# [UNVERIFIED] the US Standard Atmosphere 1976 is not on disk. The three fractions sum to
#   0.99964 by construction (trace gases dropped), as the note above says. Advisory:
#   the test suite is the only consumer, so nothing graded depends on it.
AIR_COMPOSITION = {"N2(g)": 0.78084, "O2(g)": 0.20946, "Ar(g)": 0.00934}


# HIGH-TEMPERATURE HEAT CAPACITIES FOR COMBUSTION PRODUCTS
#
#     Cp/R = A + B*T + C*T^2 + D*T^-2       (T in Kelvin)
#
# UNITS. As CP_PARAMS above: A, B, C, D are dimensionless and the polynomial
# returns Cp/R. Valid 298-3000 K.
#
# WHY A SECOND TABLE (DECISIONS D-032, D-036). CP_PARAMS is fitted to 1500 K.
# template_adiabatic_flame_temperature reaches 2908 K, so it was
# extrapolating ~1300 K past validity and every flame temperature came out low
# - methane 2311.2 K against a reference of 2326.35 K, acetylene 64.7 K below
# the value NIST gives. Phase C2 measured this (Reviewer H, F-2) and handed the
# refit-vs-restrict trade to the template owner rather than deciding it inside
# the constants track.
#
# Restricting the range is not available: the flame temperature is an OUTPUT,
# not a sampled input, so there is no knob to turn. Refitting CP_PARAMS itself
# would fix the flame template and damage template_sensible_heat_temp_dependent_cp,
# which lives at 298-1000 K where a wide-range fit is worse - for CO2 the
# wide-range residual reaches -5.4% at 298 K against -0.4% for the 1500 K fit.
#
# So the table is SPLIT BY CONSUMER. Two templates want different things from
# the same four species and one polynomial cannot serve both. CP_PARAMS is
# unchanged and keeps its low-temperature accuracy; this table is read ONLY by
# the flame template.
#
# RESULT: the flame temperatures this table produces agree with a direct NIST
# Shomate solve to within 1.93 K on all eleven reactions - worst carbon
# monoxide, 2662 against 2663.93 - where the 1500 K table was 5-65 K low.
# Methane computes 2327 K against the 2326.35 K complete-combustion
# reference: 0.03%.
#
# An earlier version of this note said 1.7 K. That was the refit-vs-NIST
# residual, not the figure that matters, which is what the TEMPLATE finally
# prints against NIST - a different and slightly larger number once the
# iteration and the rounding to whole kelvin are included (Phase 2 Reviewer C,
# finding C-1).
# @kind: property
# @units: A=1, B=1/K, C=1/K^2, D=K^2
# @domain: T=298..3000 K
CP_PARAMS_COMBUSTION = {
    # [DERIVED] least-squares refit of NIST 124-38-9 Shomate Cp over
    #   298-3000 K, 400 points. Worst residual -5.43% (at the
    #   298 K end), 1.30% above 1000 K where the flame integral
    #   has its mass.
    "CO2(g)": {"A": 5.2259, "B": 1.7079E-3, "C": -0.3275E-6, "D": -1.3180E5},   # Carbon dioxide
    # [DERIVED] least-squares refit of NIST 7732-18-5 Shomate Cp over
    #   298-3000 K, 400 points. Worst residual +0.87% (at the
    #   298 K end), 0.51% above 1000 K where the flame integral
    #   has its mass.
    #   CAVEAT: NIST's lowest gas-phase Shomate range for water starts at
    #   500 K, so the 298-500 K part of that grid is Shomate EXTRAPOLATED,
    #   not NIST data, and the +0.87% is measured against the extrapolation
    #   (Phase 2 Reviewer C, C-7). Checked and immaterial: the extrapolated
    #   Cp(298.15) is 33.590 against JANAF's 33.58 J/(mol K), 0.03%.
    "H2O(g)": {"A": 3.0520, "B": 2.2354E-3, "C": -0.3435E-6, "D": 0.3444E5},   # Water vapour
    # [DERIVED] least-squares refit of NIST 7727-37-9 Shomate Cp over
    #   298-3000 K, 400 points. Worst residual -1.05% (at the
    #   298 K end), 1.00% above 1000 K where the flame integral
    #   has its mass.
    "N2(g)": {"A": 3.1129, "B": 0.9776E-3, "C": -0.1819E-6, "D": 0.0693E5},   # Nitrogen
    # [DERIVED] least-squares refit of NIST 7782-44-7 Shomate Cp over
    #   298-3000 K, 400 points. Worst residual -3.60% (at the
    #   298 K end), 0.79% above 1000 K where the flame integral
    #   has its mass.
    "O2(g)": {"A": 3.6608, "B": 0.6232E-3, "C": -0.0847E-6, "D": -0.3842E5},   # Oxygen
}

# Validity ceiling for the table above, in Kelvin. Advisory, as CP_VALID_T_MAX
# is: the C2.2 suite checks the flame template stays inside it.
# @kind: validity
# @units: K
# @domain: none (this constant IS a domain declaration)
# [DERIVED] the top of the 298-3000 K interval the four rows above were refitted over.
CP_COMBUSTION_VALID_T_MAX = 3000.0


# A list of predefined, balanced combustion reactions with theoretical air.
# @kind: mathematical
# @units: reactants=1, products=1
# @domain: none (stoichiometry with theoretical air)
# [DERIVED] all 11 balance, checked atom by atom from the species keys, and every one recompute=atom_balance
#   carries theoretical air at N2/O2 = 3.76.
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


# Properties are given at 20°C (293.15 K) and standard pressure (101.325 kPa) unless otherwise noted.
# Viscosity can vary between sources. Values chosen for typical textbook accuracy.
# Format: { "Name": (Density [kg/m³], Dynamic Viscosity [Pa·s]) }
# @kind: property
# @units: [0]=kg/m^3, [1]=Pa*s
# @domain: T=293.15 K, p=101.325 kPa
# [UNVERIFIED] as a table - "values chosen for typical textbook accuracy", as the note above
#   says, and most of the 30 rows (honey, blood, engine oils) have no primary source on
#   disk at all. Three (row, field) pairs DO equal the NIST isobar at the declared
#   293.15 K and are tagged individually below.
#   FINDING (C3.5), DIAGNOSIS WITHDRAWN at the C3 review (Reviewer H, chemical, F-1):
#   for every organic liquid that has a NIST isobar on disk, the
#   VISCOSITY here is low by 4.6-7.2% at the very conditions this table declares - Methanol -7.06%, Benzene -7.16%, Toluene -4.62%, n-Hexane -6.12%.
#   The densities agree to 0.001-0.7%. SETTLED (C3.5 item 4): a 298.15 K grid was
#   acquired for these five species, and the answer is NEITHER reading offered.
#   The four organic VISCOSITIES fit 25 degC, and water - the control - does not:
#     methanol -7.06% at 20 degC -> +0.05% at 25 degC
#     benzene  -7.16% -> -0.33%
#     toluene  -4.62% -> +1.41%
#     n-hexane -6.12% -> -1.33%
#     WATER    +0.04% -> +12.58%   (the control, and it is decisive)
#   The DENSITY column fits 20 degC throughout. So the column is MIXED: the organic
#   viscosities are at 25 degC while water and every density are at 20 degC.
#   The review was right about the organics and wrong about the table. My own reading
#   - "a source difference, not scatter" - stays WITHDRAWN: nothing here tests it.
#   A mixed column cannot be repaired by editing @domain, because no single temperature
#   is true of the table. @domain therefore keeps 293.15 K, which is true of the
#   densities and of water. Re-sourcing the four organic viscosities at 293.15 K would
#   fix it and is a VALUE change; registered (C3.5), not made here.
COMMON_LIQUIDS = {
    # Water and Common Solvents
    # [ON-DISK] nist_fluid_properties/water_C7732185_isobar_1atm.tsv @ T=293.15 col="Density (kg/m3)" field=[0] precision=4sf
    # [ON-DISK] nist_fluid_properties/water_C7732185_isobar_1atm.tsv @ T=293.15 col="Viscosity (Pa*s)" field=[1] precision=4sf
    "Water": (998.2, 1.002e-3),
    "Seawater (3.5% salinity)": (1025, 1.07e-3), # Viscosity approx. 7% higher than pure water
    "Ethanol": (789.4, 1.074e-3),
    # [ON-DISK] nist_fluid_properties/methanol_C67561_isobar_1atm.tsv @ T=293.15 col="Viscosity (Pa*s)" field=[1] precision=3sf
    #   0.0005853 at the 293.15 K this table declares. Corrected from 0.000544
    #   (-7.06%), which is the 25 degC value: the 298.15 K grid acquired
    #   for C3.5 item 4 gives -7.06% at 20 degC against a few tenths at 25 degC.
    #   Water and the whole density column are at 20 degC, so the table could not be
    #   redeclared - these four rows had to come to the declared temperature. P6 event.
    "Methanol": (791.3, 0.000585),
    "Isopropyl Alcohol (IPA)": (781.8, 2.04e-3),
    "Acetone": (784.5, 0.306e-3),
    # [ON-DISK] nist_fluid_properties/benzene_C71432_isobar_1atm.tsv @ T=293.15 col="Viscosity (Pa*s)" field=[1] precision=3sf
    #   0.00064738 at the 293.15 K this table declares. Corrected from 0.000601
    #   (-7.16%), which is the 25 degC value: the 298.15 K grid acquired
    #   for C3.5 item 4 gives -7.16% at 20 degC against a few tenths at 25 degC.
    #   Water and the whole density column are at 20 degC, so the table could not be
    #   redeclared - these four rows had to come to the declared temperature. P6 event.
    "Benzene": (876.5, 0.000647),
    # [ON-DISK] nist_fluid_properties/toluene_C108883_isobar_1atm.tsv @ T=293.15 col="Density (kg/m3)" field=[0] precision=4sf
    # [ON-DISK] nist_fluid_properties/toluene_C108883_isobar_1atm.tsv @ T=293.15 col="Viscosity (Pa*s)" field=[1] precision=3sf
    #   0.00058714 at the 293.15 K this table declares. Corrected from 0.00056
    #   (-4.62%), which is the 25 degC value: the 298.15 K grid acquired
    #   for C3.5 item 4 gives -4.62% at 20 degC against a few tenths at 25 degC.
    #   Water and the whole density column are at 20 degC, so the table could not be
    #   redeclared - these four rows had to come to the declared temperature. P6 event.
    "Toluene": (866.9, 0.000587),
    "Diethyl Ether": (713.4, 0.223e-3),
    # [ON-DISK] nist_fluid_properties/hexane_C110543_isobar_1atm.tsv @ T=293.15 col="Viscosity (Pa*s)" field=[1] precision=3sf
    #   0.00031317 at the 293.15 K this table declares. Corrected from 0.000294
    #   (-6.12%), which is the 25 degC value: the 298.15 K grid acquired
    #   for C3.5 item 4 gives -6.12% at 20 degC against a few tenths at 25 degC.
    #   Water and the whole density column are at 20 degC, so the table could not be
    #   redeclared - these four rows had to come to the declared temperature. P6 event.
    "n-Hexane": (654.8, 0.000313),
    
    # Oils and Hydrocarbons
    "Engine Oil (SAE 10W)": (870, 0.065),
    "Engine Oil (SAE 30)": (891.0, 0.290), 
    "Engine Oil (SAE 50)": (902, 0.860),
    "Gear Oil (SAE 90)": (915, 0.700),
    "Crude Oil (light)": (850, 7.5e-3),
    "Kerosene": (810, 1.6e-3),
    "Diesel Fuel": (850, 3.5e-3),
    "Gasoline": (745, 0.29e-3), # ~0.4-0.8 cP, often approximated to water's order of magnitude
    
    # Organic & Food Grade Liquids
    "Glycerol (100%)": (1261.3, 1.490),
    "Olive Oil": (910, 0.081),
    "Corn Syrup": (1380, 5.0),     # Highly variable with concentration and temp
    "Honey": (1420, 10.0),         # Highly variable with type and temp
    "Milk (whole)": (1035, 2.0e-3),
    "Blood Plasma (human)": (1025, 1.5e-3),
    "Blood (whole, human)": (1060, 4.0e-3), # Shear-thinning, this is an approx. value
    
    # Cryogens & Liquefied Gases (at their boiling point @ 1 atm)
    "Liquid Nitrogen (77 K)": (804, 0.158e-3),   # Note: Temperature is 77 K (-196°C)
    "Liquid Oxygen (90 K)": (1141, 0.189e-3),    # Temperature is 90 K (-183°C)
    
    # Metals and Inorganics
    "Mercury": (13593, 1.526e-3),  
    "Sulfuric Acid (98%)": (1831, 25.4e-3),
    "Ethylene Glycol": (1113.4, 16.1e-3), # Common antifreeze
}


# Properties are given at 20°C (293.15 K) and 1 atm (101.325 kPa) unless otherwise noted.
# Format: { "Name": (Density [kg/m³], Dynamic Viscosity [Pa·s]) }
# @kind: property
# @units: [0]=kg/m^3, [1]=Pa*s
# @domain: T=293.15 K, p=101.325 kPa
# [UNVERIFIED] as a table; 10 (row, field) pairs equal the NIST isobar at 293.15 K and are
#   tagged below. One row is worse than that band and was missed because an automated
#   293.15 K sweep skips it (water is liquid there): "Steam (Water Vapor)" 0.747
#   kg/m^3 is +25.00% against NIST 0.59761 at the 373.15 K its own comment names, and
#   is within 0.25% of the ideal gas at 293.15 K (0.7489) - it is 20 degC vapour
#   mislabelled as 100 degC steam (Reviewer H, chemical, F-4; C3.5).
#   The rest sit 0.1-6% off it (worst: SF6 viscosity +6.06%, butane
#   viscosity +3.45%) with no on-disk source of their own. C3.5 residual.
COMMON_GASES = {
    # Common Gases
    "Air": (1.204, 1.81e-5),
    # [ON-DISK] nist_fluid_properties/nitrogen_C7727379_isobar_1atm.tsv @ T=293.15 col="Density (kg/m3)" field=[0] precision=4sf
    # [ON-DISK] nist_fluid_properties/nitrogen_C7727379_isobar_1atm.tsv @ T=293.15 col="Viscosity (Pa*s)" field=[1] precision=3sf
    "Nitrogen (N₂)": (1.165, 1.76e-5),
    # [ON-DISK] nist_fluid_properties/oxygen_C7782447_isobar_1atm.tsv @ T=293.15 col="Density (kg/m3)" field=[0] precision=4sf
    "Oxygen (O₂)": (1.331, 2.00e-5),
    # [ON-DISK] nist_fluid_properties/carbon_dioxide_C124389_isobar_1atm.tsv @ T=293.15 col="Viscosity (Pa*s)" field=[1] precision=3sf
    "Carbon Dioxide (CO₂)": (1.842, 1.47e-5), # Viscosity is temperature-dependent and increases for CO2
    # [ON-DISK] nist_fluid_properties/argon_C7440371_isobar_1atm.tsv @ T=293.15 col="Viscosity (Pa*s)" field=[1] precision=3sf
    "Argon": (1.661, 2.23e-5),
    # [ON-DISK] nist_fluid_properties/helium_C7440597_isobar_1atm.tsv @ T=293.15 col="Density (kg/m3)" field=[0] precision=3sf
    # [ON-DISK] nist_fluid_properties/helium_C7440597_isobar_1atm.tsv @ T=293.15 col="Viscosity (Pa*s)" field=[1] precision=3sf
    "Helium": (0.166, 1.96e-5),
    "Neon": (0.840, 3.18e-5),
    "Krypton": (3.479, 2.55e-5),
    "Xenon": (5.495, 2.28e-5), # Density high, but viscosity is similar to air
    
    # Hydrocarbons
    # [ON-DISK] nist_fluid_properties/methane_C74828_isobar_1atm.tsv @ T=293.15 col="Density (kg/m3)" field=[0] precision=3sf
    # [ON-DISK] nist_fluid_properties/methane_C74828_isobar_1atm.tsv @ T=293.15 col="Viscosity (Pa*s)" field=[1] precision=3sf
    "Methane (CH₄)": (0.668, 1.09e-5),
    "Ethane (C₂H₆)": (1.264, 9.15e-6),
    "Propane (C₃H₈)": (1.880, 8.00e-6), # Note: Viscosity decreases slightly with molecular weight in this series
    "Butane (C₄H₁₀)": (2.489, 7.50e-6),
    "Natural Gas (approx.)": (0.700, 1.10e-5), # Modeled after methane
    "Acetylene (C₂H₂)": (1.092, 9.80e-6),
    
    # Other Common Gases
    # [ON-DISK] nist_fluid_properties/hydrogen_C1333740_isobar_1atm.tsv @ T=293.15 col="Density (kg/m3)" field=[0] precision=3sf
    "Hydrogen (H₂)": (0.0838, 8.90e-6), # Lowest density, very low viscosity
    "Steam (Water Vapor)": (0.747, 1.02e-5), # At 100°C (373 K), 1 atm
    "Ammonia (NH₃)": (0.718, 1.01e-5),
    "Chlorine (Cl₂)": (2.994, 1.33e-5),
    "Sulfur Hexafluoride (SF₆)": (6.17, 1.59e-5), # High-density gas used in industry
    
    # Noble Gases
    "Radon": (9.23, 2.30e-5), # Theoretical value at 20°C; highly radioactive
}


# Molecular parameters for the Chapman-Enskog equation (Kinetic Theory)
# Sources: NIST, CRC Handbook, and standard chemical engineering texts.
# Format: { "Name": (Molar Mass [g/mol], Sigma σ [Å], Epsilon ε / k [K]) }
# Note: Epsilon ε / k (the Lennard-Jones energy parameter) is included for calculating the collision integral.
# @kind: property
# @units: [0]=g/mol, [1]=angstrom, [2]=K
# @domain: none (Lennard-Jones parameters; no fit temperature range is stated)
# [ON-DISK] nasa_tr_r132/svehla_1962_nasa_tr_r132.pdf @ page=22 text="FORCE CONSTANTS AND METHODS OF DETERMINATION"
#   Svehla, NASA TR R-132, Table I(a), which runs PDF pp.22-26. This tag pins the table's
#   FIRST page only - a page= locator takes one page, and each row's own page is named on
#   the row below. 16 of the 17 rows were re-sourced to it;
#   P6 event. The table previously mixed at least three compilations: sigma agreed with
#   Svehla for only 3 of 17 rows (Carbon Dioxide, Argon, Helium), and Air, N2 and H2
#   differed in BOTH columns.
#   This was classed known-defective rather than a citation nit because the sigma column is
#   wrong at the table's own purpose: Chapman-Enskog with the Neufeld collision integral
#   at 293.15 K, against the NIST isobars, gave n-butane -26.1%, SF6 -12.7%, propane
#   +3.5%, methane -0.7%, ethane +2.9%; Svehla's sigma brings all five inside 2%
#   (+0.6, +0.8, +1.2, +0.5, +1.9). Since mu ~ sigma^-2 a wrong sigma is a wrong
#   viscosity, and the template prints sigma into the question stem.
#   Literals carry Svehla's own printed digits (2.900, 3.711), not repr() of the float.
#   Three rows are not ordinary reads:
#     Cl2  the scan renders it "C_" (p.24). Accepted only because that token occurs
#          exactly once on the page, sits in alphabetical order after ClO, and its
#          eps/k 316.0 already equalled the committed value.
#     Xe   is NOT re-sourced and keeps its committed sigma. It has no token, and the
#          positional read the C3 review recorded is not reproducible: indexing p.26's
#          20th number tuple returns 2.608 / 10.22 / page 114 - helium's second
#          determination - because regex indices do not track the molecule list when
#          OCR-broken entries fail to match a clean tuple. Anchoring instead on the
#          committed eps/k 231.0 did not yield a unique tuple either. Left alone.
#     NH3  eps/k is NOT re-sourced and keeps its committed 558.3: the scan renders it
#          "55& 3", which the parser turns into 55.0. By eye it is 558.3, but by eye is
#          not a parse, so only NH3's sigma (2.900, cleanly printed) is taken. Its eps/k
#          is the one element of this table still carrying an unparsed value.
GAS_MOLECULAR_PARAMS = {
    # [ON-DISK] nasa_tr_r132/svehla_1962_nasa_tr_r132.pdf @ page=22 token="Air" col="sigma" field=[1] precision=exact
    # [ON-DISK] nasa_tr_r132/svehla_1962_nasa_tr_r132.pdf @ page=22 token="Air" col="eps_k" field=[2] precision=exact
    "Air": (28.97, 3.711, 78.6),
    # [ON-DISK] nasa_tr_r132/svehla_1962_nasa_tr_r132.pdf @ page=25 token="N_" col="sigma" field=[1] precision=exact
    # [ON-DISK] nasa_tr_r132/svehla_1962_nasa_tr_r132.pdf @ page=25 token="N_" col="eps_k" field=[2] precision=exact
    "Nitrogen (N₂)": (28.01, 3.798, 71.4),
    # [ON-DISK] nasa_tr_r132/svehla_1962_nasa_tr_r132.pdf @ page=25 token="02" col="sigma" field=[1] precision=exact
    # [ON-DISK] nasa_tr_r132/svehla_1962_nasa_tr_r132.pdf @ page=25 token="02" col="eps_k" field=[2] precision=exact
    "Oxygen (O₂)": (32.00, 3.467, 106.7),
    # [ON-DISK] nasa_tr_r132/svehla_1962_nasa_tr_r132.pdf @ page=23 token="CO2" col="sigma" field=[1] precision=exact
    # [ON-DISK] nasa_tr_r132/svehla_1962_nasa_tr_r132.pdf @ page=23 token="CO2" col="eps_k" field=[2] precision=exact
    "Carbon Dioxide (CO₂)": (44.01, 3.941, 195.2),
    # [ON-DISK] nasa_tr_r132/svehla_1962_nasa_tr_r132.pdf @ page=22 token="Ar" col="sigma" field=[1] precision=exact
    # [ON-DISK] nasa_tr_r132/svehla_1962_nasa_tr_r132.pdf @ page=22 token="Ar" col="eps_k" field=[2] precision=exact
    "Argon": (39.95, 3.542, 93.3),
    # [ON-DISK] nasa_tr_r132/svehla_1962_nasa_tr_r132.pdf @ page=24 token="He" col="sigma" field=[1] precision=exact
    # [ON-DISK] nasa_tr_r132/svehla_1962_nasa_tr_r132.pdf @ page=24 token="He" col="eps_k" field=[2] precision=exact
    "Helium": (4.003, 2.551, 10.22),
    # [ON-DISK] nasa_tr_r132/svehla_1962_nasa_tr_r132.pdf @ page=25 token="Ne" col="sigma" field=[1] precision=exact
    # [ON-DISK] nasa_tr_r132/svehla_1962_nasa_tr_r132.pdf @ page=25 token="Ne" col="eps_k" field=[2] precision=exact
    "Neon": (20.18, 2.820, 32.8),
    # [ON-DISK] nasa_tr_r132/svehla_1962_nasa_tr_r132.pdf @ page=24 token="Kr" col="sigma" field=[1] precision=exact
    # [ON-DISK] nasa_tr_r132/svehla_1962_nasa_tr_r132.pdf @ page=24 token="Kr" col="eps_k" field=[2] precision=exact
    "Krypton": (83.80, 3.655, 178.9),
    # [UNVERIFIED] sigma 4.10 is the one value in this table not re-sourced to
    #   Svehla. Xenon has no token on p.26, and the positional read the C3 review
    #   recorded is not reproducible: indexing the 20th number tuple returns
    #   2.608 / 10.22 / page 114, which is helium's second determination. Anchoring
    #   instead on the committed eps/k 231.0 found no tuple at all - the
    #   regex straddles OCR-broken entries, so p.26 cannot be read tuple-wise.
    #   Its eps/k DOES match Svehla. C3.5 residual.
    "Xenon": (131.29, 4.10, 231.0),
    # [ON-DISK] nasa_tr_r132/svehla_1962_nasa_tr_r132.pdf @ page=23 token="CH4" col="sigma" field=[1] precision=exact
    # [ON-DISK] nasa_tr_r132/svehla_1962_nasa_tr_r132.pdf @ page=23 token="CH4" col="eps_k" field=[2] precision=exact
    "Methane (CH₄)": (16.04, 3.758, 148.6),
    # [ON-DISK] nasa_tr_r132/svehla_1962_nasa_tr_r132.pdf @ page=23 token="C2H6" col="sigma" field=[1] precision=exact
    # [ON-DISK] nasa_tr_r132/svehla_1962_nasa_tr_r132.pdf @ page=23 token="C2H6" col="eps_k" field=[2] precision=exact
    "Ethane (C₂H₆)": (30.07, 4.443, 215.7),
    # [ON-DISK] nasa_tr_r132/svehla_1962_nasa_tr_r132.pdf @ page=23 token="C3Hs" col="sigma" field=[1] precision=exact
    # [ON-DISK] nasa_tr_r132/svehla_1962_nasa_tr_r132.pdf @ page=23 token="C3Hs" col="eps_k" field=[2] precision=exact
    "Propane (C₃H₈)": (44.10, 5.118, 237.1),
    # [ON-DISK] nasa_tr_r132/svehla_1962_nasa_tr_r132.pdf @ page=23 token="n-C4Hlo" col="sigma" field=[1] precision=exact
    # [ON-DISK] nasa_tr_r132/svehla_1962_nasa_tr_r132.pdf @ page=23 token="n-C4Hlo" col="eps_k" field=[2] precision=exact
    "Butane (C₄H₁₀)": (58.12, 4.687, 531.4), # n-butane
    # [ON-DISK] nasa_tr_r132/svehla_1962_nasa_tr_r132.pdf @ page=24 token="H2" col="sigma" field=[1] precision=exact
    # [ON-DISK] nasa_tr_r132/svehla_1962_nasa_tr_r132.pdf @ page=24 token="H2" col="eps_k" field=[2] precision=exact
    "Hydrogen (H₂)": (2.016, 2.827, 59.7),
    # [ON-DISK] nasa_tr_r132/svehla_1962_nasa_tr_r132.pdf @ page=25 token="NH3" col="sigma" field=[1] precision=exact
    # [UNVERIFIED] eps/k is not re-sourced: the scan prints it "55& 3" and the
    #   reader returns 55.0, so a tag claiming 558.3 would fail - and should.
    #   By eye it is 558.3; by eye is not a parse. C3.5 residual.
    "Ammonia (NH₃)": (17.03, 2.900, 558.3), # Polar molecule, value is an effective fit.
    # [ON-DISK] nasa_tr_r132/svehla_1962_nasa_tr_r132.pdf @ page=24 token="C_" col="sigma" field=[1] precision=exact
    # [ON-DISK] nasa_tr_r132/svehla_1962_nasa_tr_r132.pdf @ page=24 token="C_" col="eps_k" field=[2] precision=exact
    "Chlorine (Cl₂)": (70.90, 4.217, 316.0),
    # [ON-DISK] nasa_tr_r132/svehla_1962_nasa_tr_r132.pdf @ page=25 token="SF_" col="sigma" field=[1] precision=exact
    # [ON-DISK] nasa_tr_r132/svehla_1962_nasa_tr_r132.pdf @ page=25 token="SF_" col="eps_k" field=[2] precision=exact
    "Sulfur Hexafluoride (SF₆)": (146.06, 5.128, 222.1),
}


# Properties for non-Newtonian power-law fluids
# Format: { "Name": (Consistency Index K [Pa·s^n], Power-Law Index n [dimensionless]) }
# @kind: property
# @units: [0]=Pa*s^n, [1]=1
# @domain: none (K and n are shear-rate and temperature dependent; none is stated)
# [UNVERIFIED] ketchup, mayonnaise and toothpaste have no primary source on disk, and K for
#   a shear-thinning fluid is not a constant of the substance - it depends on the shear
#   range and temperature the fit was made over, neither of which this table records.
#   C3.5 residual.
POWER_LAW_FLUIDS = {
    # Common Household & Food (Shear-Thinning)
    "Ketchup": (32.5, 0.22),
    "Applesauce": (15.0, 0.3),
    "Mustard": (50.0, 0.28),
    "Mayonnaise": (85.0, 0.6),
    "Tomato Puree": (25.0, 0.5),
    "Yogurt": (12.5, 0.6),
    "Toothpaste": (120.0, 0.4),
    "Shampoo": (25.0, 0.6),
    "Hand Lotion": (80.0, 0.5),

    # Paints & Inks (Shear-Thinning)
    "Latex Paint": (45.0, 0.45),
    "Printing Ink": (10.0, 0.7),

    # Biological Fluids (Mostly Shear-Thinning)
    "Blood (Plasma)": (0.012, 0.95), # Very low K, nearly Newtonian
    "Mucus": (10.0, 0.5),

    # Polymer Solutions & Melts (Shear-Thinning)
    "0.5% Carboxymethylcellulose (CMC) in Water": (1.5, 0.6),
    "1.5% Polyacrylamide in Water": (5.0, 0.3),
    "Molten Polyethylene": (5000.0, 0.6), # K is very temperature-dependent

    # Newtonian Baseline (n = 1)
    "Water": (0.001, 1.0),           # K is the dynamic viscosity
    "Glycerol": (1.0, 1.0),
    "Air": (1.8e-5, 1.0),            # K is the dynamic viscosity

    # Shear-Thickening (Dilatant)
    "Corn Starch Suspension (40%)": (2.0, 1.9),
    "Corn Starch Suspension (50%)": (10.0, 2.5), # Higher concentration -> stronger effect
    "Silica Sand Suspension (60%)": (0.5, 1.6),
}


# Standard gravitational acceleration in m/s²
# @kind: defined
# @units: m/s^2
# @domain: none (a defined standard value, not a measurement at a condition)
# [ON-DISK] codata_2022/allascii.txt @ quantity="standard acceleration of gravity" precision=3sf
GRAVITATIONAL_ACCELERATION = 9.81
