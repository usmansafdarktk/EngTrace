# E is the Modulus of Elasticity.
# @kind: property
# @units: E_GPa=GPa, E_ksi=ksi, nu=1
# @domain: none (the table states no temperature, form or temper; each cited row names what its source covers)
MATERIAL_PROPERTIES = {
    # Metals (Common)
    # [UNVERIFIED] names no grade. MIL-HDBK-5J's carbon-steel table (p.62, AISI 1025) prints E
    #   29.0, G 11.0 x10^3 ksi, mu 0.32; choosing a grade is the repo owner's call (C3.5).
    'Steel': {'E_GPa': 200, 'E_ksi': 29000, 'nu': 0.30},
    # [KNOWN-DEFECTIVE] ['E_ksi'] 27500 against MIL-HDBK-5J p.277, Table 2.7.1.0(b) (AISI 301and Relateda,b,c Stainless Steels): E:L 29.0 x10^3 ksi, -5.2%.
    #   The row names what the table covers, so the literal should be a rounding of it and is not. Annealed column; the table's footnote b applies it to AISI 304 (AMS 5513), the row comment's alloy.
    #   Not corrected in this commit: a value change moves emitted items (P6); outcome in the C3.5 register.
    # [KNOWN-DEFECTIVE] ['E_GPa'] 190 against MIL-HDBK-5J p.277, Table 2.7.1.0(b) (AISI 301and Relateda,b,c Stainless Steels): E:L 29.0 x10^3 ksi = 199.9 GPa
    #   (SP 811 p.63: 1 ksi = 6.894757e+06 Pa, so 6.894757 GPa per 10^3 ksi); -5.0%, outside half a unit in the handbook's last digit (0.18%).
    #   Not corrected in this commit: a value change moves emitted items (P6); outcome in the C3.5 register.
    # [KNOWN-DEFECTIVE] ['nu'] 0.30 against MIL-HDBK-5J p.277, Table 2.7.1.0(b) (AISI 301and Relateda,b,c Stainless Steels): mu 0.27, +11.1%.
    #   The row names what the table covers, so the literal should be a rounding of it and is not.
    #   Not corrected in this commit: a value change moves emitted items (P6); outcome in the C3.5 register.
    'Stainless Steel': {'E_GPa': 190, 'E_ksi': 27500, 'nu': 0.30}, # 304 Stainless
    # [UNVERIFIED] the row comment names 1100 and 3003; MIL-HDBK-5J has no design table for either
    #   (its wrought aluminum tables are the 2xxx, 5xxx, 6xxx, 7xxx series). Referred (C3.5).
    'Aluminum': {'E_GPa': 69, 'E_ksi': 10000, 'nu': 0.33}, # General purpose 1100, 3003 Al
    # [KNOWN-DEFECTIVE] ['E_ksi'] 10000 against MIL-HDBK-5J p.566, Table 3.6.2.0(b1) (6061 Aluminum Alloy Sheet): E 9.9 x10^3 ksi, +1.0%.
    #   The row names what the table covers, so the literal should be a rounding of it and is not.
    #   Not corrected in this commit: a value change moves emitted items (P6); outcome in the C3.5 register.
    # [KNOWN-DEFECTIVE] ['E_GPa'] 68.9 against MIL-HDBK-5J p.566, Table 3.6.2.0(b1) (6061 Aluminum Alloy Sheet): E 9.9 x10^3 ksi = 68.26 GPa
    #   (SP 811 p.63: 1 ksi = 6.894757e+06 Pa, so 6.894757 GPa per 10^3 ksi); +0.9%, outside half a unit in the handbook's last digit (0.51%).
    #   Not corrected in this commit: a value change moves emitted items (P6); outcome in the C3.5 register.
    # [ON-DISK] mil_hdbk_5j/MIL-HDBK-5J_2003-01-31.pdf @ page=566 text="Table 3.6.2.0(b1)" mil="mu" field=['nu'] precision=exact
    #   MIL-HDBK-5J p.566, Table 3.6.2.0(b1) (6061 Aluminum Alloy Sheet): mu 0.33.
    'Aluminum 6061-T6': {'E_GPa': 68.9, 'E_ksi': 10000, 'nu': 0.33}, # A common specific alloy
    # [UNVERIFIED] no MIL-HDBK-5J design table is this material (its copper-base tables: C86500
    #   Manganese Bronze; C86300 Manganese Bronze; Copper Beryllium Strip; C17200 Copper Beryllium
    #   Rod and Bar; C17200 Copper Beryllium Mechanical Tubing), and no other source on disk.
    #   Residual register (C3.5).
    'Copper': {'E_GPa': 117, 'E_ksi': 17000, 'nu': 0.34},
    # [UNVERIFIED] no MIL-HDBK-5J design table is this material (its copper-base tables: C86500
    #   Manganese Bronze; C86300 Manganese Bronze; Copper Beryllium Strip; C17200 Copper Beryllium
    #   Rod and Bar; C17200 Copper Beryllium Mechanical Tubing), and no other source on disk.
    #   Residual register (C3.5).
    'Brass': {'E_GPa': 102, 'E_ksi': 14800, 'nu': 0.34}, # Yellow Brass (CuZn37)
    # [UNVERIFIED] no MIL-HDBK-5J design table is this material (its copper-base tables: C86500
    #   Manganese Bronze; C86300 Manganese Bronze; Copper Beryllium Strip; C17200 Copper Beryllium
    #   Rod and Bar; C17200 Copper Beryllium Mechanical Tubing), and no other source on disk.
    #   Residual register (C3.5).
    'Bronze': {'E_GPa': 110, 'E_ksi': 16000, 'nu': 0.34}, # Phosphor Bronze
    # [KNOWN-DEFECTIVE] ['E_ksi'] 16800 against MIL-HDBK-5J p.899, Table 5.2.1.0(b) (Commercially Pure Titanium): E 15.5 x10^3 ksi, +8.4%.
    #   The row names what the table covers, so the literal should be a rounding of it and is not. The row comment's commercially pure titanium.
    #   Not corrected in this commit: a value change moves emitted items (P6); outcome in the C3.5 register.
    # [KNOWN-DEFECTIVE] ['E_GPa'] 116 against MIL-HDBK-5J p.899, Table 5.2.1.0(b) (Commercially Pure Titanium): E 15.5 x10^3 ksi = 106.9 GPa
    #   (SP 811 p.63: 1 ksi = 6.894757e+06 Pa, so 6.894757 GPa per 10^3 ksi); +8.5%, outside half a unit in the handbook's last digit (0.33%).
    #   Not corrected in this commit: a value change moves emitted items (P6); outcome in the C3.5 register.
    # [UNVERIFIED] MIL-HDBK-5J p.899, Table 5.2.1.0(b) (Commercially Pure Titanium) prints no mu ("..."); no other source on disk
    'Titanium': {'E_GPa': 116, 'E_ksi': 16800, 'nu': 0.34}, # Commercially Pure
    # [KNOWN-DEFECTIVE] ['E_ksi'] 16500 against MIL-HDBK-5J p.945, Table 5.4.1.0(b) (Ti-6Al-4V Sheet, Strip, and Plate): E 16.0 x10^3 ksi, +3.1%.
    #   The row names what the table covers, so the literal should be a rounding of it and is not. Sheet, strip and plate; the bar table (p.946) prints E 16.9.
    #   Not corrected in this commit: a value change moves emitted items (P6); outcome in the C3.5 register.
    # [KNOWN-DEFECTIVE] ['E_GPa'] 114 against MIL-HDBK-5J p.945, Table 5.4.1.0(b) (Ti-6Al-4V Sheet, Strip, and Plate): E 16.0 x10^3 ksi = 110.3 GPa
    #   (SP 811 p.63: 1 ksi = 6.894757e+06 Pa, so 6.894757 GPa per 10^3 ksi); +3.3%, outside half a unit in the handbook's last digit (0.32%).
    #   Not corrected in this commit: a value change moves emitted items (P6); outcome in the C3.5 register.
    # [KNOWN-DEFECTIVE] ['nu'] 0.34 against MIL-HDBK-5J p.945, Table 5.4.1.0(b) (Ti-6Al-4V Sheet, Strip, and Plate): mu 0.31, +9.7%.
    #   The row names what the table covers, so the literal should be a rounding of it and is not.
    #   Not corrected in this commit: a value change moves emitted items (P6); outcome in the C3.5 register.
    'Titanium Alloy (6Al-4V)': {'E_GPa': 114, 'E_ksi': 16500, 'nu': 0.34},
    # [ON-DISK] mil_hdbk_5j/MIL-HDBK-5J_2003-01-31.pdf @ page=841 text="Table 4.2.1.0(b)" mil="E" field=['E_ksi'] scale=1e3 precision=exact
    #   MIL-HDBK-5J p.841, Table 4.2.1.0(b) (AZ31B Magnesium Alloy Sheet and Plate): E 6.5 x10^3 ksi. The row comment's AZ31B, sheet and plate.
    # [ON-DISK] mil_hdbk_5j/MIL-HDBK-5J_2003-01-31.pdf @ page=841 text="Table 4.2.1.0(b)" mil="E" field=['E_GPa'] scale=6.894757 tol=0.77%
    #   MIL-HDBK-5J p.841, Table 4.2.1.0(b) (AZ31B Magnesium Alloy Sheet and Plate): E 6.5 x10^3 ksi = 44.82 GPa. The row is +0.41% from it, inside
    #   half a unit in the handbook's last digit (0.77%) - the bound a unit conversion of a
    #   2-s.f. value supports (SP 811 p.63: 1 ksi = 6.894757e+06 Pa, so 6.894757 GPa per 10^3 ksi).
    # [KNOWN-DEFECTIVE] ['nu'] 0.29 against MIL-HDBK-5J p.841, Table 4.2.1.0(b) (AZ31B Magnesium Alloy Sheet and Plate): mu 0.35, -17.1%.
    #   The row names what the table covers, so the literal should be a rounding of it and is not.
    #   Not corrected in this commit: a value change moves emitted items (P6); outcome in the C3.5 register.
    'Magnesium': {'E_GPa': 45, 'E_ksi': 6500, 'nu': 0.29}, # AZ31B alloy
    # [UNVERIFIED] no MIL-HDBK-5J design table is this material (its copper-base tables: C86500
    #   Manganese Bronze; C86300 Manganese Bronze; Copper Beryllium Strip; C17200 Copper Beryllium
    #   Rod and Bar; C17200 Copper Beryllium Mechanical Tubing), and no other source on disk.
    #   Residual register (C3.5).
    'Tungsten': {'E_GPa': 411, 'E_ksi': 59600, 'nu': 0.28},
    # [UNVERIFIED] no MIL-HDBK-5J design table is this material (its copper-base tables: C86500
    #   Manganese Bronze; C86300 Manganese Bronze; Copper Beryllium Strip; C17200 Copper Beryllium
    #   Rod and Bar; C17200 Copper Beryllium Mechanical Tubing), and no other source on disk.
    #   Residual register (C3.5).
    'Cast Iron': {'E_GPa': 170, 'E_ksi': 24600, 'nu': 0.26}, # Gray Cast Iron
    # [UNVERIFIED] no MIL-HDBK-5J design table is this material (its copper-base tables: C86500
    #   Manganese Bronze; C86300 Manganese Bronze; Copper Beryllium Strip; C17200 Copper Beryllium
    #   Rod and Bar; C17200 Copper Beryllium Mechanical Tubing), and no other source on disk.
    #   Residual register (C3.5).
    'Nickel': {'E_GPa': 207, 'E_ksi': 30000, 'nu': 0.31},
    # [UNVERIFIED] no MIL-HDBK-5J design table is this material (its copper-base tables: C86500
    #   Manganese Bronze; C86300 Manganese Bronze; Copper Beryllium Strip; C17200 Copper Beryllium
    #   Rod and Bar; C17200 Copper Beryllium Mechanical Tubing), and no other source on disk.
    #   Residual register (C3.5).
    'Lead': {'E_GPa': 16, 'E_ksi': 2300, 'nu': 0.44}, # Added a common soft metal
    
    # Polymers/Plastics
    # [UNVERIFIED] no source on disk for this material. Residual register (C3.5).
    'Nylon': {'E_GPa': 2.1, 'E_ksi': 300, 'nu': 0.39}, # Nylon 6/6
    # [UNVERIFIED] no source on disk for this material. Residual register (C3.5).
    'Polycarbonate': {'E_GPa': 2.3, 'E_ksi': 334, 'nu': 0.38},
    # [UNVERIFIED] no source on disk for this material. Residual register (C3.5).
    'ABS': {'E_GPa': 2.0, 'E_ksi': 290, 'nu': 0.35},
    # [UNVERIFIED] no source on disk for this material. Residual register (C3.5).
    'PVC (rigid)': {'E_GPa': 3.0, 'E_ksi': 435, 'nu': 0.38},
    # [UNVERIFIED] no source on disk for this material. Residual register (C3.5).
    'PTFE (Teflon)': {'E_GPa': 0.5, 'E_ksi': 72, 'nu': 0.46},
    # [UNVERIFIED] no source on disk for this material. Residual register (C3.5).
    'Polyethylene (HDPE)': {'E_GPa': 1.1, 'E_ksi': 160, 'nu': 0.42},
    # [UNVERIFIED] no source on disk for this material. Residual register (C3.5).
    'Epoxy': {'E_GPa': 3.0, 'E_ksi': 435, 'nu': 0.38}, # Unreinforced epoxy resin
    # [UNVERIFIED] no source on disk for this material. Residual register (C3.5).
    'Natural Rubber': {'E_GPa': 0.0015, 'E_ksi': 0.22, 'nu': 0.4999}, # ~0.5 (nearly incompressible)
    
    # Composites & Other
    # [UNVERIFIED] no source on disk for this material. Residual register (C3.5).
    'Carbon Fiber Reinforced Polymer (CFRP)': {'E_GPa': 150, 'E_ksi': 21750, 'nu': 0.30}, # Direction-dependent, this is a typical in-plane value
    # [UNVERIFIED] no source on disk for this material. Residual register (C3.5).
    'Fiberglass (GFRP)': {'E_GPa': 45, 'E_ksi': 6500, 'nu': 0.25}, # Direction-dependent, typical in-plane value
    # [UNVERIFIED] no source on disk for this material. Residual register (C3.5).
    'Concrete': {'E_GPa': 30, 'E_ksi': 4350, 'nu': 0.15}, # Highly variable, common approx.
    # [UNVERIFIED] no source on disk for this material. Residual register (C3.5).
    'Glass (Borosilicate)': {'E_GPa': 70, 'E_ksi': 10150, 'nu': 0.20},
    # [UNVERIFIED] no source on disk for this material. Residual register (C3.5).
    'Ceramic (Alumina Al2O3)': {'E_GPa': 370, 'E_ksi': 53700, 'nu': 0.22},
}


# Typical Shear Modulus (G) Values for Engineering Materials (GPa)
# Note: Values are representative and can vary with alloy composition and heat treatment.
# @kind: property
# @units: GPa
# @domain: none (the table states no temperature, form or temper; each cited row names what its source covers)
SHEAR_MODULUS_VALUES = {
    # Metals (Common Alloys)
    # [UNVERIFIED] ASTM A36 has no MIL-HDBK-5J table; its carbon-steel table (p.62, AISI 1025)
    #   prints G 11.0 x10^3 ksi = 75.84 GPa. Residual register (C3.5).
    "Steel (A36)": 77.2,
    # [ON-DISK] mil_hdbk_5j/MIL-HDBK-5J_2003-01-31.pdf @ page=277 text="Table 2.7.1.0(b)" mil="G" scale=6.894757 tol=0.45%
    #   MIL-HDBK-5J p.277, Table 2.7.1.0(b) (AISI 301and Relateda,b,c Stainless Steels): G 11.2 x10^3 ksi = 77.22 GPa. The row is -0.29% from it, inside
    #   half a unit in the handbook's last digit (0.45%) - the bound a unit conversion of a
    #   3-s.f. value supports (SP 811 p.63: 1 ksi = 6.894757e+06 Pa, so 6.894757 GPa per 10^3 ksi). Annealed column; footnote b applies the table to AISI 304 (AMS 5513).
    "Stainless Steel (304)": 77.0,
    # [ON-DISK] mil_hdbk_5j/MIL-HDBK-5J_2003-01-31.pdf @ page=566 text="Table 3.6.2.0(b1)" mil="G" scale=6.894757 tol=1.32%
    #   MIL-HDBK-5J p.566, Table 3.6.2.0(b1) (6061 Aluminum Alloy Sheet): G 3.8 x10^3 ksi = 26.2 GPa. The row is -0.76% from it, inside
    #   half a unit in the handbook's last digit (1.32%) - the bound a unit conversion of a
    #   2-s.f. value supports (SP 811 p.63: 1 ksi = 6.894757e+06 Pa, so 6.894757 GPa per 10^3 ksi).
    "Aluminum 6061-T6": 26.0,
    # [UNVERIFIED] MIL-HDBK-5J's 2024 sheet table carrying T4 (p.374) does not parse from its text
    #   layer (4 labels and values [0, 0, 0, 0] fit neither layout). The 2024 tables that do parse
    #   print G 4.0 (p.373), 4.0 (p.382), 4.0 (p.384), 4.0 (p.385), 4.0 (p.386), 4.1 (p.387), 4.1
    #   (p.388) x10^3 ksi, i.e. 27.58 GPa for 4.0; whether T4's G is the same is not read here.
    #   Residual register (C3.5).
    "Aluminum 2024-T4": 28.0,
    # [UNVERIFIED] no MIL-HDBK-5J design table is this material, and no other source on disk.
    #   Residual register (C3.5).
    "Brass (C36000)": 39.0,
    # [UNVERIFIED] no MIL-HDBK-5J design table is this material, and no other source on disk.
    #   Residual register (C3.5).
    "Copper (CDA 110)": 44.7,
    # [UNVERIFIED] no MIL-HDBK-5J design table is this material, and no other source on disk.
    #   Residual register (C3.5).
    "Bronze (Phosphor 510)": 41.4,
    # [KNOWN-DEFECTIVE] value 41.4 against MIL-HDBK-5J p.945, Table 5.4.1.0(b) (Ti-6Al-4V Sheet, Strip, and Plate): G 6.2 x10^3 ksi = 42.75 GPa
    #   (SP 811 p.63: 1 ksi = 6.894757e+06 Pa, so 6.894757 GPa per 10^3 ksi); -3.2%, outside half a unit in the handbook's last digit (0.81%). Sheet, strip and plate; the bar table (p.946) prints G 6.2 too.
    #   Not corrected in this commit: a value change moves emitted items (P6); outcome in the C3.5 register.
    "Titanium Alloy (Ti-6Al-4V)": 41.4,
    # [ON-DISK] mil_hdbk_5j/MIL-HDBK-5J_2003-01-31.pdf @ page=841 text="Table 4.2.1.0(b)" mil="G" scale=6.894757 tol=2.09%
    #   MIL-HDBK-5J p.841, Table 4.2.1.0(b) (AZ31B Magnesium Alloy Sheet and Plate): G 2.4 x10^3 ksi = 16.55 GPa. The row is -0.29% from it, inside
    #   half a unit in the handbook's last digit (2.09%) - the bound a unit conversion of a
    #   2-s.f. value supports (SP 811 p.63: 1 ksi = 6.894757e+06 Pa, so 6.894757 GPa per 10^3 ksi). Sheet and plate.
    "Magnesium Alloy (AZ31B)": 16.5,
    # [UNVERIFIED] no MIL-HDBK-5J design table is this material, and no other source on disk.
    #   Residual register (C3.5).
    "Tungsten": 160.0,
    # [UNVERIFIED] no MIL-HDBK-5J design table is this material, and no other source on disk.
    #   Residual register (C3.5).
    "Molybdenum": 118.0,
    # [UNVERIFIED] no MIL-HDBK-5J design table is this material, and no other source on disk.
    #   Residual register (C3.5).
    "Lead": 5.9,
    
    # Other Materials
    # [UNVERIFIED] no MIL-HDBK-5J design table is this material, and no other source on disk.
    #   Residual register (C3.5).
    "Gray Cast Iron": 44.0, # Note: Anisotropic property, value is an approximation.
    # [UNVERIFIED] no MIL-HDBK-5J design table is this material, and no other source on disk.
    #   Residual register (C3.5).
    "Nylon 6/6": 1.1,      # Polymers have much lower moduli
    # [UNVERIFIED] no MIL-HDBK-5J design table is this material, and no other source on disk.
    #   Residual register (C3.5).
    "Polycarbonate": 0.9,
    # [UNVERIFIED] no MIL-HDBK-5J design table is this material, and no other source on disk.
    #   Residual register (C3.5).
    "Concrete": 22.0,      # Highly variable, this is a common estimate for calculation
    # [UNVERIFIED] no MIL-HDBK-5J design table is this material, and no other source on disk.
    #   Residual register (C3.5).
    "Glass": 26.2,
    
    # Composites (Highly variable based on fiber orientation and volume fraction)
    # Value given is an approximate in-plane shear modulus.
    # [UNVERIFIED] no MIL-HDBK-5J design table is this material, and no other source on disk.
    #   Residual register (C3.5).
    "Carbon Fiber Epoxy (Unidirectional, in-plane)": 5.0,
}

# Standard atmospheric pressure in kPa
# @kind: defined
# @units: kPa
# @domain: none (a defined standard value)
# [ON-DISK] codata_2022/allascii.txt @ quantity="standard atmosphere" scale=1e-3 precision=exact
#   CODATA 2022 prints 101325 Pa (exact).
# @copied-in: template_hydrostatic_pressure_at_depth 101.325
#   C3.8: fluid_statics.py writes P_atm_kpa = 101.325 instead of reading this table;
#   declared so the census counts the copy (phaseC3_literal_copies.md).
ATMOSPHERIC_PRESSURE_KPA = 101.325

# Standard acceleration due to gravity in m/s^2
# @kind: defined
# @units: m/s^2
# @domain: none (the standard value, fixed by definition; the corpus does not model local gravity)
# [ON-DISK] codata_2022/allascii.txt @ quantity="standard acceleration of gravity" precision=3sf
#   CODATA 2022 prints 9.80665 m s^-2 (exact); 9.81 is its 3-s.f. rounding.
GRAVITY = 9.81


# Densities of common fluids in kg/m^3 at 20°C and 1 atm, unless specified.
# @kind: property
# @units: kg/m^3
# @domain: T=293.15 K, P=101.325 kPa
#   the header's "20°C and 1 atm"; rows that name their own conditions say so in their tags
FLUID_DENSITIES = {
    # Water-based
    # [ON-DISK] nist_fluid_properties/water_C7732185_isobar_1atm.tsv @ T=293.15 col="Density (kg/m3)" precision=3sf
    #   998.21 at 293.15 K.
    "Fresh Water": 998,         # More precise value for 20°C
    # [UNVERIFIED] a mixture, solution or commercial product with no single composition in the row;
    #   no on-disk source, and none was searched for in C3. Residual register (C3.5).
    "Sea Water": 1025,
    # [UNVERIFIED] a mixture, solution or commercial product with no single composition in the row;
    #   no on-disk source, and none was searched for in C3. Residual register (C3.5).
    "Brine (20% NaCl)": 1150,   # Common industrial solution
    
    # Oils and Hydrocarbons
    # [UNVERIFIED] a mixture, solution or commercial product with no single composition in the row;
    #   no on-disk source, and none was searched for in C3. Residual register (C3.5).
    "SAE 10W Oil": 870,
    # [UNVERIFIED] a mixture, solution or commercial product with no single composition in the row;
    #   no on-disk source, and none was searched for in C3. Residual register (C3.5).
    "SAE 30 Oil": 917,
    # [UNVERIFIED] a mixture, solution or commercial product with no single composition in the row;
    #   no on-disk source, and none was searched for in C3. Residual register (C3.5).
    "Crude Oil": 870,           # Typical average value
    # [UNVERIFIED] a mixture, solution or commercial product with no single composition in the row;
    #   no on-disk source, and none was searched for in C3. Residual register (C3.5).
    "Kerosene": 810,
    # [UNVERIFIED] a mixture, solution or commercial product with no single composition in the row;
    #   no on-disk source, and none was searched for in C3. Residual register (C3.5).
    "Diesel Fuel": 850,
    # [UNVERIFIED] a mixture, solution or commercial product with no single composition in the row;
    #   no on-disk source, and none was searched for in C3. Residual register (C3.5).
    "Gasoline": 726,
    # [UNVERIFIED] a mixture, solution or commercial product with no single composition in the row;
    #   no on-disk source, and none was searched for in C3. Residual register (C3.5).
    "Engine Oil": 888,
    # [UNVERIFIED] a mixture, solution or commercial product with no single composition in the row;
    #   no on-disk source, and none was searched for in C3. Residual register (C3.5).
    "Hydraulic Oil": 850,       # Typical average value
    
    # Alcohols and Solvents
    # [UNVERIFIED] Ethanol is not in the NIST fluid database on disk; its WebBook species pages on
    #   disk (ethanol_C64175_condensed_phase.html, ethanol_C64175_phase_change.html) say density
    #   only as ['Critical density']. Residual register (C3.5).
    "Ethyl Alcohol (Ethanol)": 789,
    # [ON-DISK] nist_fluid_properties/methanol_C67561_isobar_1atm.tsv @ T=293.15 col="Density (kg/m3)" precision=3sf
    #   791.01 at 293.15 K.
    "Methyl Alcohol (Methanol)": 791,
    # [UNVERIFIED] Isopropanol is not in the NIST fluid database on disk; its WebBook species pages
    #   on disk (isopropanol_C67630_condensed_phase.html, isopropanol_C67630_phase_change.html) say
    #   density only as ['Critical density']. Residual register (C3.5).
    "Isopropyl Alcohol (IPA)": 786,
    # [UNVERIFIED] Acetone is not in the NIST fluid database on disk; its WebBook species pages on
    #   disk (acetone_C67641_condensed_phase.html, acetone_C67641_phase_change.html) say density
    #   only as ['Critical density']. Residual register (C3.5).
    "Acetone": 791,
    # [UNVERIFIED] a mixture, solution or commercial product with no single composition in the row;
    #   no on-disk source, and none was searched for in C3. Residual register (C3.5).
    "Turpentine": 870,
    # [KNOWN-DEFECTIVE] 876 against NIST 878.92 at 293.15 K (nist_fluid_properties/benzene_C71432_isobar_1atm.tsv), -0.33%: the table states
    #   these conditions, and 876 is not a 3-s.f. rounding of 878.92.
    #   Not corrected in this commit: a value change moves emitted items (P6); outcome in the C3.5 register.
    "Benzene": 876,
    # [ON-DISK] nist_fluid_properties/toluene_C108883_isobar_1atm.tsv @ T=293.15 col="Density (kg/m3)" precision=3sf
    #   866.89 at 293.15 K.
    "Toluene": 867,
    # [UNVERIFIED] Xylene is not in the NIST fluid database on disk; its WebBook species pages on
    #   disk (p_xylene_C106423_condensed_phase.html, p_xylene_C106423_phase_change.html) say
    #   density only as ['Critical density']. Residual register (C3.5).
    "Xylene": 870,
    
    # Cryogens & Liquefied Gases
    # [ON-DISK] nist_fluid_properties/nitrogen_C7727379_saturation_77.15K.tsv @ T=77.15 col="Density (l, kg/m3)" precision=3sf
    #   807.01 at 77.15 K. The row states -196 C under the table's 1 atm; the file is the saturated liquid at 77.15 K. Pressure in that row: 0.098899 MPa.
    "Liquid Nitrogen (at -196°C)": 807,
    # [ON-DISK] nist_fluid_properties/oxygen_C7782447_saturation_90.15K.tsv @ T=90.15 col="Density (l, kg/m3)" precision=4sf
    #   1141.4 at 90.15 K. The row states -183 C under the table's 1 atm; the file is the saturated liquid at 90.15 K. Pressure in that row: 0.10093 MPa.
    "Liquid Oxygen (at -183°C)": 1141,
    # [ON-DISK] nist_fluid_properties/hydrogen_C1333740_saturation_20.15K.tsv @ T=20.15 col="Density (l, kg/m3)" precision=2sf
    #   71.096 at 20.15 K. The row states -253 C under the table's 1 atm; the file is the saturated liquid at 20.15 K. Pressure in that row: 0.094929 MPa.
    "Liquid Hydrogen (at -253°C)": 71,
    # [ON-DISK] nist_fluid_properties/propane_C74986_saturation_298.15K.tsv @ T=298.15 col="Density (l, kg/m3)" tol=0.13%
    #   492.36 at 298.15 K; the row is +0.13%. The row says "At 25 C, under pressure" and names no pressure; the file is the saturated liquid. Pressure in that row: 0.95207 MPa.
    "Liquid Propane": 493,      # At 25°C, under pressure
    
    # High-Density Liquids
    # [UNVERIFIED] Mercury is not in the NIST fluid database on disk; its WebBook species pages on
    #   disk (mercury_C7439976_condensed_phase.html, mercury_C7439976_phase_change.html) carry no
    #   density. Residual register (C3.5).
    "Mercury": 13550,
    # [UNVERIFIED] Glycerol is not in the NIST fluid database on disk; its WebBook species pages on
    #   disk (glycerol_C56815_condensed_phase.html, glycerol_C56815_phase_change.html) carry no
    #   density. Residual register (C3.5).
    "Glycerin": 1260,
    # [UNVERIFIED] Chloroform is not in the NIST fluid database on disk (MANIFEST: 36 fluids), and
    #   no WebBook species page for it is on disk. Residual register (C3.5).
    "Chloroform": 1489,
    # [UNVERIFIED] Carbon tetrachloride is not in the NIST fluid database on disk; its WebBook
    #   species pages on disk (carbon_tetrachloride_C56235_condensed_phase.html,
    #   carbon_tetrachloride_C56235_phase_change.html) say density only as ['Critical density', 'of
    #   density', 'the density']. Residual register (C3.5).
    "Carbon Tetrachloride": 1594,
    # [UNVERIFIED] Bromine is not in the NIST fluid database on disk (MANIFEST: 36 fluids), and no
    #   WebBook species page for it is on disk. Residual register (C3.5).
    "Bromine": 3120,
    # [UNVERIFIED] a mixture, solution or commercial product with no single composition in the row;
    #   no on-disk source, and none was searched for in C3. Residual register (C3.5).
    "Sulfuric Acid (98%)": 1830,
    # [UNVERIFIED] a mixture, solution or commercial product with no single composition in the row;
    #   no on-disk source, and none was searched for in C3. Residual register (C3.5).
    "Hydrochloric Acid (37%)": 1190, # Also known as Muriatic Acid
    # [UNVERIFIED] a mixture, solution or commercial product with no single composition in the row;
    #   no on-disk source, and none was searched for in C3. Residual register (C3.5).
    "Nitric Acid (68%)": 1350,
    
    # Common Liquids
    # [UNVERIFIED] a mixture, solution or commercial product with no single composition in the row;
    #   no on-disk source, and none was searched for in C3. Residual register (C3.5).
    "Milk": 1035,
    # [UNVERIFIED] a mixture, solution or commercial product with no single composition in the row;
    #   no on-disk source, and none was searched for in C3. Residual register (C3.5).
    "Whole Blood": 1060,        # Typical value for human blood
    # [UNVERIFIED] Ethylene glycol is not in the NIST fluid database on disk; its WebBook species
    #   pages on disk (ethylene_glycol_C107211_condensed_phase.html,
    #   ethylene_glycol_C107211_phase_change.html) carry no density. Residual register (C3.5).
    "Ethylene Glycol (Antifreeze)": 1113,
    # [UNVERIFIED] a mixture, solution or commercial product with no single composition in the row;
    #   no on-disk source, and none was searched for in C3. Residual register (C3.5).
    "Corn Syrup": 1380,
    # [UNVERIFIED] a mixture, solution or commercial product with no single composition in the row;
    #   no on-disk source, and none was searched for in C3. Residual register (C3.5).
    "Olive Oil": 910,
    # [UNVERIFIED] a mixture, solution or commercial product with no single composition in the row;
    #   no on-disk source, and none was searched for in C3. Residual register (C3.5).
    "Vegetable Oil": 920,
    # [KNOWN-DEFECTIVE] 1206 against NIST 1206.7 at 298.15 K (nist_fluid_properties/r134a_C811972_saturation_298.15K.tsv), -0.06%: the table states
    #   these conditions, and 1206 is not a 4-s.f. rounding of 1206.7. The row states saturated liquid at 25 C, which is what the file is. Pressure in that row: 0.66538 MPa.
    #   Not corrected in this commit: a value change moves emitted items (P6); outcome in the C3.5 register.
    "R-134a Refrigerant (Saturated Liquid at 25°C)": 1206,
    
    # Gases (at 1 atm, for comparison - though not for hydrostatic problems!)
    # [UNVERIFIED] Air is not in the NIST fluid database on disk (MANIFEST: 36 fluids), and no
    #   WebBook species page for it is on disk. Residual register (C3.5).
    "Air (at 20°C)": 1.2,
    # [ON-DISK] nist_fluid_properties/helium_C7440597_isobar_1atm.tsv @ T=293.15 col="Density (kg/m3)" precision=3sf
    #   0.16631 at 293.15 K.
    "Helium (at 20°C)": 0.166,
    # [ON-DISK] nist_fluid_properties/carbon_dioxide_C124389_isobar_1atm.tsv @ T=293.15 col="Density (kg/m3)" precision=3sf
    #   1.8393 at 293.15 K.
    "Carbon Dioxide (at 20°C)": 1.84
}

# Densities of various solid materials in kg/m^3
# @kind: property
# @units: kg/m^3
# @domain: none (the table states no temperature)
MATERIAL_DENSITIES = {
    # Woods & Natural Materials (generally float in water)
    # [UNVERIFIED] names no species and no moisture content, and a wood density is not one number
    #   without both. The USDA Wood Handbook is on disk but was not read to row level in C3.
    #   Residual register (C3.5).
    "Pine Wood": 500,
    # [UNVERIFIED] names no species and no moisture content, and a wood density is not one number
    #   without both. The USDA Wood Handbook is on disk but was not read to row level in C3.
    #   Residual register (C3.5).
    "Oak Wood": 750,
    # [UNVERIFIED] names no species and no moisture content, and a wood density is not one number
    #   without both. The USDA Wood Handbook is on disk but was not read to row level in C3.
    #   Residual register (C3.5).
    "Cork": 240,
    # [UNVERIFIED] names no species and no moisture content, and a wood density is not one number
    #   without both. The USDA Wood Handbook is on disk but was not read to row level in C3.
    #   Residual register (C3.5).
    "Teak Wood": 630,
    # [UNVERIFIED] names no species and no moisture content, and a wood density is not one number
    #   without both. The USDA Wood Handbook is on disk but was not read to row level in C3.
    #   Residual register (C3.5).
    "Maple Wood": 740,
    # [UNVERIFIED] names no species and no moisture content, and a wood density is not one number
    #   without both. The USDA Wood Handbook is on disk but was not read to row level in C3.
    #   Residual register (C3.5).
    "Ebony Wood": 1200,  # Sinks in water
    # [UNVERIFIED] names no species and no moisture content, and a wood density is not one number
    #   without both. The USDA Wood Handbook is on disk but was not read to row level in C3.
    #   Residual register (C3.5).
    "Bamboo": 300,
    # [UNVERIFIED] no source on disk for this material. Residual register (C3.5).
    "Rubber (Natural)": 950,
    
    # Plastics & Synthetics (range from floating to sinking)
    # [UNVERIFIED] no source on disk for this material. Residual register (C3.5).
    "Polypropylene (PP)": 900,
    # [UNVERIFIED] no source on disk for this material. Residual register (C3.5).
    "Polyethylene (LDPE)": 920,
    # [UNVERIFIED] no source on disk for this material. Residual register (C3.5).
    "Polyethylene (HDPE)": 950,
    # [UNVERIFIED] no source on disk for this material. Residual register (C3.5).
    "Nylon": 1150,
    # [UNVERIFIED] no source on disk for this material. Residual register (C3.5).
    "Polystyrene (PS)": 1050,
    # [UNVERIFIED] no source on disk for this material. Residual register (C3.5).
    "PVC (Rigid)": 1380,
    # [UNVERIFIED] no source on disk for this material. Residual register (C3.5).
    "PTFE (Teflon)": 2200,
    # [UNVERIFIED] no source on disk for this material. Residual register (C3.5).
    "Acrylic (Plexiglas)": 1180,
    
    # Biological Materials
    # [UNVERIFIED] no source on disk for this material. Residual register (C3.5).
    "Human Body (avg.)": 985,  # Slightly less than water; most people float
    # [UNVERIFIED] no source on disk for this material. Residual register (C3.5).
    "Apple": 800,
    # [UNVERIFIED] no source on disk for this material. Residual register (C3.5).
    "Bone": 1850,
    
    # Lightweight & Porous Materials
    # [UNVERIFIED] no source on disk for this material. Residual register (C3.5).
    "Aerogel": 3,
    # [UNVERIFIED] no source on disk for this material. Residual register (C3.5).
    "Pumice": 700,  # The only rock that floats!
    # [UNVERIFIED] no source on disk for this material. Residual register (C3.5).
    "Brick": 1700,
    # [UNVERIFIED] no source on disk for this material. Residual register (C3.5).
    "Concrete (avg.)": 2400,
    
    # Common Metals & Alloys (will sink in water, float in mercury)
    # [UNVERIFIED] names no alloy: MIL-HDBK-5J prints 0.098 lb/in^3 for 6061 (p.566) = 2713 kg/m^3
    #   and 0.100 for 2024 (p.373) = 2768. Choosing one is the repo owner's call (C3.5).
    "Aluminum": 2710,
    # [ON-DISK] mil_hdbk_5j/MIL-HDBK-5J_2003-01-31.pdf @ page=899 text="Table 5.2.1.0(b)" mil="density" scale=27679.9 tol=0.31%
    #   MIL-HDBK-5J p.899, Table 5.2.1.0(b) (Commercially Pure Titanium): density 0.163 lb/in^3 = 4512 kg/m3. The row is -0.26% from it, inside
    #   half a unit in the handbook's last digit (0.31%) - the bound a unit conversion of a
    #   3-s.f. value supports (SP 811 p.66: 1 lb/in^3 = 27679.9 kg/m^3). Commercially pure; the row names no grade - the implementer's choice, referred (C3.5).
    "Titanium": 4500,
    # [UNVERIFIED] no source on disk for this material. Residual register (C3.5).
    "Zinc": 7140,
    # [UNVERIFIED] no source on disk for this material. Residual register (C3.5).
    "Tin": 7280,
    # [UNVERIFIED] no source on disk for this material. Residual register (C3.5).
    "Iron (Wrought)": 7750,
    # [ON-DISK] mil_hdbk_5j/MIL-HDBK-5J_2003-01-31.pdf @ page=62 text="Table 2.2.1.0(b)" mil="density" scale=27679.9 tol=0.18%
    #   MIL-HDBK-5J p.62, Table 2.2.1.0(b) (AISI 1025 Carbon Steel): density 0.284 lb/in^3 = 7861 kg/m3. The row is -0.14% from it, inside
    #   half a unit in the handbook's last digit (0.18%) - the bound a unit conversion of a
    #   3-s.f. value supports (SP 811 p.66: 1 lb/in^3 = 27679.9 kg/m^3). AISI 1025, the handbook's carbon-steel table; the row names no grade, so the grade is the implementer's choice - referred (C3.5).
    "Steel (Carbon)": 7850,
    # [UNVERIFIED] no source on disk for this material. Residual register (C3.5).
    "Brass": 8600,
    # [UNVERIFIED] no source on disk for this material. Residual register (C3.5).
    "Copper": 8940,
    # [UNVERIFIED] no source on disk for this material. Residual register (C3.5).
    "Nickel": 8900,
    # [UNVERIFIED] no source on disk for this material. Residual register (C3.5).
    "Silver": 10500,
    # [UNVERIFIED] no source on disk for this material. Residual register (C3.5).
    "Lead": 11340,
    # [UNVERIFIED] no source on disk for this material. Residual register (C3.5).
    "Uranium": 19100,
    # [UNVERIFIED] no source on disk for this material. Residual register (C3.5).
    "Gold": 19300,
    # [UNVERIFIED] no source on disk for this material. Residual register (C3.5).
    "Tungsten": 19600,
    # [UNVERIFIED] no source on disk for this material. Residual register (C3.5).
    "Platinum": 21450,
    # [UNVERIFIED] no source on disk for this material. Residual register (C3.5).
    "Osmium": 22590,  # The densest naturally occurring element
    
    # Other Materials
    # [UNVERIFIED] no source on disk for this material. Residual register (C3.5).
    "Ice (0°C)": 917,
    # [UNVERIFIED] no source on disk for this material. Residual register (C3.5).
    "Wax (Paraffin)": 900,
    # [UNVERIFIED] no source on disk for this material. Residual register (C3.5).
    "Glass (Window)": 2500,
    # [UNVERIFIED] no source on disk for this material. Residual register (C3.5).
    "Glass (Crystal)": 2900,
    # [UNVERIFIED] no source on disk for this material. Residual register (C3.5).
    "Graphite": 2100,
    # [UNVERIFIED] no source on disk for this material. Residual register (C3.5).
    "Diamond": 3500,
    # [UNVERIFIED] no source on disk for this material. Residual register (C3.5).
    "Quartz": 2650,
    # [UNVERIFIED] no source on disk for this material. Residual register (C3.5).
    "Granite": 2700,
    # [UNVERIFIED] no source on disk for this material. Residual register (C3.5).
    "Marble": 2700,
    # [UNVERIFIED] no source on disk for this material. Residual register (C3.5).
    "Sandstone": 2300,
    # [UNVERIFIED] no source on disk for this material. Residual register (C3.5).
    "Limestone": 2500,
    # [UNVERIFIED] names no species and no moisture content, and a wood density is not one number
    #   without both. The USDA Wood Handbook is on disk but was not read to row level in C3.
    #   Residual register (C3.5).
    "Cork Board": 240,
    # [UNVERIFIED] no source on disk for this material. Residual register (C3.5).
    "Paper": 800,
    # [UNVERIFIED] no source on disk for this material. Residual register (C3.5).
    "Leather (Dry)": 860,
    # [UNVERIFIED] no source on disk for this material. Residual register (C3.5).
    "Chalk": 2500,
    # [UNVERIFIED] no source on disk for this material. Residual register (C3.5).
    "Asphalt": 1100,
    # [UNVERIFIED] no source on disk for this material. Residual register (C3.5).
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
# @kind: property
# @units: kg/m^3
# @domain: none (the table states no conditions)
#   Measured against NIST at 1 atm, its rows agree at their own precision at DIFFERENT temperatures: Hydrogen at 273.15 K; Natural Gas (Methane) at 273.15 K; Carbon Dioxide at 273.15 K; Methanol at 293.15 K; Toluene at 293.15 K.
#   One table, at least two implied states: a D-032 finding, referred to the repo owner (C3.5).
PIPE_FLUIDS = {
    # Gases (Very Light)
    # [UNVERIFIED] no source on disk. Residual register (C3.5).
    "Air": 1.225,
    # [UNVERIFIED] the table states no temperature or pressure, so no row of an artefact is THE
    #   value. NIST at 1 atm: 0.17848 (vapor) at 273.15 K; 0.16631 (vapor) at 293.15 K
    #   (nist_fluid_properties/helium_C7440597_isobar_1atm.tsv). Residual register (C3.5).
    "Helium": 0.1786,
    # [UNVERIFIED] the table states no temperature or pressure, so no row of an artefact is THE
    #   value. NIST at 1 atm: 0.089883 (vapor) at 273.15 K; 0.083752 (vapor) at 293.15 K
    #   (nist_fluid_properties/hydrogen_C1333740_isobar_1atm.tsv). Residual register (C3.5).
    "Hydrogen": 0.0899,
    # [UNVERIFIED] the table states no temperature or pressure, so no row of an artefact is THE
    #   value. NIST at 1 atm: 0.71746 (vapor) at 273.15 K; 0.66816 (vapor) at 293.15 K
    #   (nist_fluid_properties/methane_C74828_isobar_1atm.tsv). Residual register (C3.5).
    "Natural Gas (Methane)": 0.717,
    # [UNVERIFIED] the table states no temperature or pressure, so no row of an artefact is THE
    #   value. NIST at 1 atm: 1.9768 (vapor) at 273.15 K; 1.8393 (vapor) at 293.15 K
    #   (nist_fluid_properties/carbon_dioxide_C124389_isobar_1atm.tsv). Residual register (C3.5).
    "Carbon Dioxide": 1.98,
    
    # Light Hydrocarbons & Fuels
    # [UNVERIFIED] no source on disk. Residual register (C3.5).
    "Gasoline": 726,
    # [UNVERIFIED] no source on disk. Residual register (C3.5).
    "Kerosene": 810,
    # [UNVERIFIED] no source on disk. Residual register (C3.5).
    "Diesel Fuel": 850,
    # [UNVERIFIED] no source on disk. Residual register (C3.5).
    "Jet Fuel (JP-4)": 770,
    # [UNVERIFIED] no source on disk. Residual register (C3.5).
    "Ethanol": 789,
    # [UNVERIFIED] the table states no temperature or pressure, so no row of an artefact is THE
    #   value. NIST at 1 atm: 809.73 (liquid) at 273.15 K; 791.01 (liquid) at 293.15 K
    #   (nist_fluid_properties/methanol_C67561_isobar_1atm.tsv). Residual register (C3.5).
    "Methanol": 791,
    # [UNVERIFIED] no source on disk. Residual register (C3.5).
    "Isopropyl Alcohol": 786,
    
    # Oils & Lubricants
    # [UNVERIFIED] no source on disk. Residual register (C3.5).
    "SAE 10 Oil": 870,
    # [UNVERIFIED] no source on disk. Residual register (C3.5).
    "SAE 20 Oil": 880,
    # [UNVERIFIED] no source on disk. Residual register (C3.5).
    "SAE 30 Oil": 917,
    # [UNVERIFIED] no source on disk. Residual register (C3.5).
    "SAE 40 Oil": 945,
    # [UNVERIFIED] no source on disk. Residual register (C3.5).
    "SAE 50 Oil": 960,
    # [UNVERIFIED] no source on disk. Residual register (C3.5).
    "Crude Oil (Light)": 825,
    # [UNVERIFIED] no source on disk. Residual register (C3.5).
    "Crude Oil (Heavy)": 950,
    # [UNVERIFIED] no source on disk. Residual register (C3.5).
    "Engine Oil": 888,
    # [UNVERIFIED] no source on disk. Residual register (C3.5).
    "Hydraulic Oil": 860,
    
    # Common Liquids
    # [UNVERIFIED] the table states no temperature or pressure, so no row of an artefact is THE
    #   value. NIST at 1 atm: no row at 273.15 K; 998.21 (liquid) at 293.15 K
    #   (nist_fluid_properties/water_C7732185_isobar_1atm.tsv). Residual register (C3.5).
    "Water": 1000,
    # [UNVERIFIED] no source on disk. Residual register (C3.5).
    "Sea Water": 1025,
    # [UNVERIFIED] no source on disk. Residual register (C3.5).
    "Milk": 1030,
    # [UNVERIFIED] no source on disk. Residual register (C3.5).
    "Ethylene Glycol": 1110,
    # [UNVERIFIED] no source on disk. Residual register (C3.5).
    "Antifreeze (50/50)": 1065,
    
    # Cryogenic Fluids
    # [UNVERIFIED] the table states no temperature or pressure for a cryogen. Residual register
    #   (C3.5).
    "Liquid Nitrogen": 810,
    # [UNVERIFIED] the table states no temperature or pressure for a cryogen. Residual register
    #   (C3.5).
    "Liquid Oxygen": 1141,
    
    # Chemical Solvents
    # [UNVERIFIED] no source on disk. Residual register (C3.5).
    "Acetone": 784,
    # [UNVERIFIED] the table states no temperature or pressure, so no row of an artefact is THE
    #   value. NIST at 1 atm: 885.42 (liquid) at 273.15 K; 866.89 (liquid) at 293.15 K
    #   (nist_fluid_properties/toluene_C108883_isobar_1atm.tsv). Residual register (C3.5).
    "Toluene": 867,
    # [UNVERIFIED] the table states no temperature or pressure, so no row of an artefact is THE
    #   value. NIST at 1 atm: no row at 273.15 K; 878.92 (liquid) at 293.15 K
    #   (nist_fluid_properties/benzene_C71432_isobar_1atm.tsv). Residual register (C3.5).
    "Benzene": 876,
    # [UNVERIFIED] the table states no temperature or pressure, so no row of an artefact is THE
    #   value. NIST at 1 atm: 677.19 (liquid) at 273.15 K; 659.38 (liquid) at 293.15 K
    #   (nist_fluid_properties/hexane_C110543_isobar_1atm.tsv). Residual register (C3.5).
    "Hexane": 655
}

# Denser fluids suitable for the manometer
# @kind: property
# @units: kg/m^3
# @domain: none (the table states no conditions)
MANOMETER_FLUIDS = {
    # Standard Manometer Fluids
    # [UNVERIFIED] Mercury is not in the NIST fluid database on disk; its WebBook species pages on
    #   disk (mercury_C7439976_condensed_phase.html, mercury_C7439976_phase_change.html) carry no
    #   density. Residual register (C3.5).
    "Mercury": 13550,
    # [UNVERIFIED] the table states no temperature or pressure, so no row of an artefact is THE
    #   value. NIST at 1 atm: no row at 273.15 K; 998.21 (liquid) at 293.15 K
    #   (nist_fluid_properties/water_C7732185_isobar_1atm.tsv). Residual register (C3.5).
    "Water": 1000,  # For gas measurements
    # [UNVERIFIED] no source on disk. Residual register (C3.5).
    "Sea Water": 1025,
    
    # Heavy Oils & Lubricants
    # [UNVERIFIED] no source on disk. Residual register (C3.5).
    "SAE 50 Oil": 960,
    # [UNVERIFIED] no source on disk. Residual register (C3.5).
    "SAE 90 Gear Oil": 920,
    
    # Organic Liquids
    # [UNVERIFIED] no source on disk. Residual register (C3.5).
    "Carbon Tetrachloride": 1590,
    # [UNVERIFIED] no source on disk. Residual register (C3.5).
    "Chloroform": 1480,
    # [UNVERIFIED] no source on disk. Residual register (C3.5).
    "Bromine": 3120,
    # [UNVERIFIED] no source on disk. Residual register (C3.5).
    "Tetrabromoethane": 2960,
    
    # Inorganic Solutions
    # [UNVERIFIED] no source on disk. Residual register (C3.5).
    "Calcium Chloride Solution (40%)": 1390,
    # [UNVERIFIED] no source on disk. Residual register (C3.5).
    "Zinc Chloride Solution (50%)": 1520,
    # [UNVERIFIED] no source on disk. Residual register (C3.5).
    "Sodium Polysulfide": 1650,
    
    # Specialty Fluids
    # [UNVERIFIED] no source on disk. Residual register (C3.5).
    "Glycerin": 1260,
    # [UNVERIFIED] no source on disk. Residual register (C3.5).
    "Diiodomethane": 3325,
    # [UNVERIFIED] no source on disk. Residual register (C3.5).
    "Acetylene Tetrabromide": 2960,
    
    # Molten Metals (for high-temperature applications)
    # [UNVERIFIED] no source on disk. Residual register (C3.5).
    "Gallium": 6095,
    # [UNVERIFIED] no source on disk. Residual register (C3.5).
    "Tin": 6980,
    # [UNVERIFIED] no source on disk. Residual register (C3.5).
    "Zinc": 6570,
    
    # Very Heavy Fluids
    # [UNVERIFIED] the WebBook name search returns "Name Not Found" for "tellurium mercury";
    #   whether the row means mercury telluride (a WebBook species) is not established. Referred to
    #   the repo owner (C3.5).
    "Tellurium Mercury": 8100,
    # [UNVERIFIED] its phase at manometer conditions, which the brief asks, from a source: the
    #   WebBook page on disk (tungsten_hexafluoride_C7783826_phase_change.html) gives Antoine
    #   parameters for 201.5-290.4 K under its own form log10(P/bar) = A - B/(T + C), A = 4.55569,
    #   B = 1021.208, C = -64.7. At 290.4 K, the top of that range, P = 1.074 bar, already above 1
    #   atm = 1.01325 bar (CODATA 101325 Pa; SP 811: 1 bar = 1.0E+05 Pa). A vapour pressure above 1
    #   atm at 290.4 K, and rising with T, makes it a GAS at 1 atm and room temperature, not a
    #   manometer liquid - [DERIVED] from those two artefacts, not read as a stated boiling point
    #   (the page states none). The row is a substance defect, not a value defect; replacing it is
    #   the repo owner's call (D-033, C3.5).
    "Tungsten Hexafluoride": 12900
}
