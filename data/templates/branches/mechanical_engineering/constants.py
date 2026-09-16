# E is the Modulus of Elasticity.
# @kind: property
# @units: E_GPa=GPa, E_ksi=ksi, nu=1
# @domain: none (the table states no temperature, form or temper; each cited row names what its source covers)
MATERIAL_PROPERTIES = {
    # Metals (Common)
    # C3 review, Reviewer H (mechanical) F2: this row was tagged unverified, "names no
    # grade",
    # while its own tag quoted the page that grades it. The branch accepts an implementer
    # grade choice when it RESOLVES the steel density row against this same table and when
    # it CONDEMNS the stainless rows; declining only here is the inconsistency. Steel is now
    # treated as Magnesium already is. p.62 prints E 29.0, G 11.0 x10^3 ksi, mu 0.32 (its
    # text layer spaces every digit: "0 . 3 2").
    # [ON-DISK] mil_hdbk_5j/MIL-HDBK-5J_2003-01-31.pdf @ page=62 text="Table 2.2.1.0(b)" mil="E" field=['E_ksi'] scale=1e3 precision=exact
    # [ON-DISK] mil_hdbk_5j/MIL-HDBK-5J_2003-01-31.pdf @ page=62 text="Table 2.2.1.0(b)" mil="E" field=['E_GPa'] scale=6.894757 tol=0.17% basis=half-unit
    #   29.0 x10^3 ksi = 199.95 GPa; the row's 200 is +0.03%, inside half a unit in
    #   the handbook's last printed digit (0.17%).
    # [ON-DISK] mil_hdbk_5j/MIL-HDBK-5J_2003-01-31.pdf @ page=62 text="Table 2.2.1.0(b)" mil="mu" field=['nu'] precision=exact
    #   AISI 1025 Carbon Steel: mu 0.32. Its text layer spaces every digit ("0 . 3 2").
    #   Corrected from 0.30 (-6.25%) on the repo owner's instruction; P6 event, instance
    #   dump recorded with the C3.5 outcome.
    'Steel': {'E_GPa': 200, 'E_ksi': 29000, 'nu': 0.32},
    # [ON-DISK] mil_hdbk_5j/MIL-HDBK-5J_2003-01-31.pdf @ page=277 text="Table 2.7.1.0(b)" mil="E:L" field=['E_ksi'] scale=1e3 precision=exact
    #   AISI 301 and Related Stainless Steels, annealed: E:L 29.0 x10^3 ksi. The table's
    #   footnote b applies it to AISI 304 (AMS 5513), the row comment's alloy.
    #   Corrected from 27500 (-5.2%); P6 event.
    # [ON-DISK] mil_hdbk_5j/MIL-HDBK-5J_2003-01-31.pdf @ page=277 text="Table 2.7.1.0(b)" mil="E:L" field=['E_GPa'] scale=6.894757 tol=0.17% basis=half-unit
    #   29.0 x10^3 ksi = 199.95 GPa (SP 811 p.63: 1 ksi = 6.894757e+06 Pa, so 6.894757 GPa per 10^3 ksi); the row is a 1-decimal
    #   rounding of it, inside half a unit in the handbook's last printed digit.
    #   Corrected from 190 (-5.0%); P6 event.
    # [ON-DISK] mil_hdbk_5j/MIL-HDBK-5J_2003-01-31.pdf @ page=277 text="Table 2.7.1.0(b)" mil="mu" field=['nu'] precision=exact
    #   AISI 301 and Related Stainless Steels, annealed: mu 0.27.
    #   Corrected from 0.30 (+11.1%); P6 event.
    'Stainless Steel': {'E_GPa': 199.9, 'E_ksi': 29000, 'nu': 0.27}, # 304 Stainless
    # [UNVERIFIED] the row comment names 1100 and 3003; MIL-HDBK-5J has no design table for either
    #   (its wrought aluminum tables are the 2xxx, 5xxx, 6xxx, 7xxx series). Referred (C3.5).
    'Aluminum': {'E_GPa': 69, 'E_ksi': 10000, 'nu': 0.33}, # General purpose 1100, 3003 Al
    # [ON-DISK] mil_hdbk_5j/MIL-HDBK-5J_2003-01-31.pdf @ page=566 text="Table 3.6.2.0(b1)" mil="E" field=['E_ksi'] scale=1e3 precision=exact
    #   6061 Aluminum Alloy Sheet: E 9.9 x10^3 ksi. Corrected from 10000 (+1.0%); P6 event.
    # [ON-DISK] mil_hdbk_5j/MIL-HDBK-5J_2003-01-31.pdf @ page=566 text="Table 3.6.2.0(b1)" mil="E" field=['E_GPa'] scale=6.894757 tol=0.51% basis=half-unit
    #   9.9 x10^3 ksi = 68.26 GPa (SP 811 p.63: 1 ksi = 6.894757e+06 Pa, so 6.894757 GPa per 10^3 ksi); the row is a 1-decimal
    #   rounding of it. Corrected from 68.9 (+0.9%); P6 event.
    # [ON-DISK] mil_hdbk_5j/MIL-HDBK-5J_2003-01-31.pdf @ page=566 text="Table 3.6.2.0(b1)" mil="mu" field=['nu'] precision=exact
    #   MIL-HDBK-5J p.566, Table 3.6.2.0(b1) (6061 Aluminum Alloy Sheet): mu 0.33.
    'Aluminum 6061-T6': {'E_GPa': 68.3, 'E_ksi': 9900, 'nu': 0.33}, # A common specific alloy
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
    # [ON-DISK] mil_hdbk_5j/MIL-HDBK-5J_2003-01-31.pdf @ page=899 text="Table 5.2.1.0(b)" mil="E" field=['E_ksi'] scale=1e3 precision=exact
    #   Commercially Pure Titanium: E 15.5 x10^3 ksi, the row comment's material.
    #   Corrected from 16800 (+8.4%); P6 event.
    # [ON-DISK] mil_hdbk_5j/MIL-HDBK-5J_2003-01-31.pdf @ page=899 text="Table 5.2.1.0(b)" mil="E" field=['E_GPa'] scale=6.894757 tol=0.32% basis=half-unit
    #   15.5 x10^3 ksi = 106.87 GPa (SP 811 p.63: 1 ksi = 6.894757e+06 Pa, so 6.894757 GPa per 10^3 ksi); the row is a 1-decimal
    #   rounding of it. Corrected from 116 (+8.5%); P6 event.
    # [UNVERIFIED] MIL-HDBK-5J p.899, Table 5.2.1.0(b) (Commercially Pure Titanium) prints no mu ("..."); no other source on disk
    'Titanium': {'E_GPa': 106.9, 'E_ksi': 15500, 'nu': 0.34}, # Commercially Pure
    # [ON-DISK] mil_hdbk_5j/MIL-HDBK-5J_2003-01-31.pdf @ page=945 text="Table 5.4.1.0(b)" mil="E" field=['E_ksi'] scale=1e3 precision=exact
    #   Ti-6Al-4V Sheet, Strip, and Plate: E 16.0 x10^3 ksi. Sheet, strip and plate; the
    #   bar table (p.946) prints E 16.9. Corrected from 16500 (+3.1%); P6 event.
    # [ON-DISK] mil_hdbk_5j/MIL-HDBK-5J_2003-01-31.pdf @ page=945 text="Table 5.4.1.0(b)" mil="E" field=['E_GPa'] scale=6.894757 tol=0.31% basis=half-unit
    #   16.0 x10^3 ksi = 110.32 GPa (SP 811 p.63: 1 ksi = 6.894757e+06 Pa, so 6.894757 GPa per 10^3 ksi); the row is a 1-decimal
    #   rounding of it. Corrected from 114 (+3.3%); P6 event.
    # [ON-DISK] mil_hdbk_5j/MIL-HDBK-5J_2003-01-31.pdf @ page=945 text="Table 5.4.1.0(b)" mil="mu" field=['nu'] precision=exact
    #   Ti-6Al-4V Sheet, Strip, and Plate: mu 0.31. Corrected from 0.34 (+9.7%); P6 event.
    'Titanium Alloy (6Al-4V)': {'E_GPa': 110.3, 'E_ksi': 16000, 'nu': 0.31},
    # [ON-DISK] mil_hdbk_5j/MIL-HDBK-5J_2003-01-31.pdf @ page=841 text="Table 4.2.1.0(b)" mil="E" field=['E_ksi'] scale=1e3 precision=exact
    #   MIL-HDBK-5J p.841, Table 4.2.1.0(b) (AZ31B Magnesium Alloy Sheet and Plate): E 6.5 x10^3 ksi. The row comment's AZ31B, sheet and plate.
    # [ON-DISK] mil_hdbk_5j/MIL-HDBK-5J_2003-01-31.pdf @ page=841 text="Table 4.2.1.0(b)" mil="E" field=['E_GPa'] scale=6.894757 tol=0.77% basis=half-unit
    #   MIL-HDBK-5J p.841, Table 4.2.1.0(b) (AZ31B Magnesium Alloy Sheet and Plate): E 6.5 x10^3 ksi = 44.82 GPa. The row is +0.41% from it, inside
    #   half a unit in the handbook's last digit (0.77%) - the bound a unit conversion of a
    #   2-s.f. value supports (SP 811 p.63: 1 ksi = 6.894757e+06 Pa, so 6.894757 GPa per 10^3 ksi).
    # [ON-DISK] mil_hdbk_5j/MIL-HDBK-5J_2003-01-31.pdf @ page=841 text="Table 4.2.1.0(b)" mil="mu" field=['nu'] precision=exact
    #   AZ31B Magnesium Alloy Sheet and Plate: mu 0.35, the row comment's alloy.
    #   Corrected from 0.29 (-17.1%); P6 event.
    'Magnesium': {'E_GPa': 45, 'E_ksi': 6500, 'nu': 0.35}, # AZ31B alloy
    # [UNVERIFIED] no MIL-HDBK-5J design table is this material under the names searched.
    #   C3 review, H (mechanical) F3: the reason here previously cited the handbook's
    #   COPPER-BASE tables as the search set. That is apt for Copper, Brass and Bronze
    #   and is no evidence at all about this material; the conclusion may hold, the
    #   evidence offered did not. Residual register (C3.5).
    'Tungsten': {'E_GPa': 411, 'E_ksi': 59600, 'nu': 0.28},
    # [UNVERIFIED] no MIL-HDBK-5J design table is this material under the names searched.
    #   C3 review, H (mechanical) F3: the reason here previously cited the handbook's
    #   COPPER-BASE tables as the search set. That is apt for Copper, Brass and Bronze
    #   and is no evidence at all about this material; the conclusion may hold, the
    #   evidence offered did not. Residual register (C3.5).
    'Cast Iron': {'E_GPa': 170, 'E_ksi': 24600, 'nu': 0.26}, # Gray Cast Iron
    # [UNVERIFIED] no MIL-HDBK-5J design table is this material under the names searched.
    #   C3 review, H (mechanical) F3: the reason here previously cited the handbook's
    #   COPPER-BASE tables as the search set. That is apt for Copper, Brass and Bronze
    #   and is no evidence at all about this material; the conclusion may hold, the
    #   evidence offered did not. Residual register (C3.5).
    #   C3 review, H (mechanical) F3: for NICKEL the conclusion is OPEN, not settled -
    #   the handbook carries 22 nickel-base design captions (A-286 p.1045, Hastelloy X
    #   p.1061, Inconel 600 pp.1067-1070), none of them commercially pure nickel.
    'Nickel': {'E_GPa': 207, 'E_ksi': 30000, 'nu': 0.31},
    # [UNVERIFIED] no MIL-HDBK-5J design table is this material under the names searched.
    #   C3 review, H (mechanical) F3: the reason here previously cited the handbook's
    #   COPPER-BASE tables as the search set. That is apt for Copper, Brass and Bronze
    #   and is no evidence at all about this material; the conclusion may hold, the
    #   evidence offered did not. Residual register (C3.5).
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
    # [ON-DISK] mil_hdbk_5j/MIL-HDBK-5J_2003-01-31.pdf @ page=277 text="Table 2.7.1.0(b)" mil="G" scale=6.894757 tol=0.45% basis=half-unit
    #   MIL-HDBK-5J p.277, Table 2.7.1.0(b) (AISI 301and Relateda,b,c Stainless Steels): G 11.2 x10^3 ksi = 77.22 GPa. The row is -0.29% from it, inside
    #   half a unit in the handbook's last digit (0.45%) - the bound a unit conversion of a
    #   3-s.f. value supports (SP 811 p.63: 1 ksi = 6.894757e+06 Pa, so 6.894757 GPa per 10^3 ksi). Annealed column; footnote b applies the table to AISI 304 (AMS 5513).
    "Stainless Steel (304)": 77.0,
    # [ON-DISK] mil_hdbk_5j/MIL-HDBK-5J_2003-01-31.pdf @ page=566 text="Table 3.6.2.0(b1)" mil="G" scale=6.894757 tol=1.32% basis=half-unit
    #   MIL-HDBK-5J p.566, Table 3.6.2.0(b1) (6061 Aluminum Alloy Sheet): G 3.8 x10^3 ksi = 26.2 GPa. The row is -0.76% from it, inside
    #   half a unit in the handbook's last digit (1.32%) - the bound a unit conversion of a
    #   2-s.f. value supports (SP 811 p.63: 1 ksi = 6.894757e+06 Pa, so 6.894757 GPa per 10^3 ksi).
    "Aluminum 6061-T6": 26.0,
    # [ON-DISK] mil_hdbk_5j/MIL-HDBK-5J_2003-01-31.pdf @ page=373 text="Table 3.2.3.0(b1)" mil="G" scale=6.894757 tol=1.25% basis=half-unit
    #   Bare 2024 Aluminum Alloy Sheet and Plate: G 4.0 x10^3 ksi = 27.58 GPa
    #   (SP 811 p.63: 1 ksi = 6.894757e+06 Pa), the same value in all 4 columns of p.373.
    #   WHY THIS APPLIES TO T4. This row's previous tag said the 2024 table "carrying T4
    #   (p.374) does not parse from its text layer". It does not parse because there is
    #   nothing on it to parse: p.374 and p.375 print "See Table 3.2.3.0(d)" in all four
    #   elastic cells. Following that deferral, p.376 carries "Table 3.2.3.0(d). Modulus
    #   Values and Poisson's Ratio for Bare 2024 Aluminum Alloy Sheet and Plate, ALL
    #   TEMPERS", printing E 10.5/10.7, Ec 10.7/10.9, G 4.0/4.0, mu 0.33/0.33 for the two
    #   thickness bands. So G does not vary by temper here, which is what the old tag left
    #   open ("whether T4's G is the same is not read here").
    #   p.376 is cited in prose and not as the locator because it is a Property/Thickness
    #   matrix with no "E, 10^3 ksi" leader-dot segment: milhdbk.design_values cannot read
    #   it, and a precision= tag against it would be unverifiable. p.373 carries the same
    #   G and IS machine-readable.
    #   Corrected from 28.0 (+1.5%, outside the 1.25% half-unit bound); P6 event.
    "Aluminum 2024-T4": 27.6,
    # [UNVERIFIED] no MIL-HDBK-5J design table is this material, and no other source on disk.
    #   Residual register (C3.5).
    "Brass (C36000)": 39.0,
    # [UNVERIFIED] no MIL-HDBK-5J design table is this material, and no other source on disk.
    #   Residual register (C3.5).
    "Copper (CDA 110)": 44.7,
    # [UNVERIFIED] no MIL-HDBK-5J design table is this material, and no other source on disk.
    #   Residual register (C3.5).
    "Bronze (Phosphor 510)": 41.4,
    # [ON-DISK] mil_hdbk_5j/MIL-HDBK-5J_2003-01-31.pdf @ page=945 text="Table 5.4.1.0(b)" mil="G" scale=6.894757 tol=0.81% basis=half-unit
    #   Ti-6Al-4V Sheet, Strip, and Plate: G 6.2 x10^3 ksi = 42.75 GPa
    #   (SP 811 p.63: 1 ksi = 6.894757e+06 Pa, so 6.894757 GPa per 10^3 ksi); the row is a 1-decimal rounding of it. Sheet, strip and plate; the bar
    #   table (p.946) prints G 6.2 too. Corrected from 41.4 (-3.2%); P6 event.
    "Titanium Alloy (Ti-6Al-4V)": 42.7,
    # [ON-DISK] mil_hdbk_5j/MIL-HDBK-5J_2003-01-31.pdf @ page=841 text="Table 4.2.1.0(b)" mil="G" scale=6.894757 tol=2.09% basis=half-unit
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
    # [UNVERIFIED] a mixture or commercial grade with no single composition, so no
    #   pure substance on disk IS this row. The nearest pure n-alkanes are now on
    #   disk and bracket it: nonane 718.03, decane 730.41, dodecane 749.44 kg/m3 at 293.15 K, 1 atm
    #   (nist_fluid_properties/<name>_isobar_1atm.tsv). This row holds 870, which is
    #   above that C9-C12 span.
    #   They were fetched to re-point rows like this one, and reading them settled
    #   that plan against itself: pointing kerosene at decane moves it 810 -> 730,
    #   an 11% error, because a real cut carries aromatics and cycloalkanes a pure
    #   n-alkane does not. Sourcing a row must not change what the row names.
    #   C3.5 residual - searched for, and the search is recorded.
    "SAE 10W Oil": 870,
    # [UNVERIFIED] a mixture or commercial grade with no single composition, so no
    #   pure substance on disk IS this row. The nearest pure n-alkanes are now on
    #   disk and bracket it: nonane 718.03, decane 730.41, dodecane 749.44 kg/m3 at 293.15 K, 1 atm
    #   (nist_fluid_properties/<name>_isobar_1atm.tsv). This row holds 917, which is
    #   above that C9-C12 span.
    #   They were fetched to re-point rows like this one, and reading them settled
    #   that plan against itself: pointing kerosene at decane moves it 810 -> 730,
    #   an 11% error, because a real cut carries aromatics and cycloalkanes a pure
    #   n-alkane does not. Sourcing a row must not change what the row names.
    #   C3.5 residual - searched for, and the search is recorded.
    "SAE 30 Oil": 917,
    # [UNVERIFIED] a mixture or commercial grade with no single composition, so no
    #   pure substance on disk IS this row. The nearest pure n-alkanes are now on
    #   disk and bracket it: nonane 718.03, decane 730.41, dodecane 749.44 kg/m3 at 293.15 K, 1 atm
    #   (nist_fluid_properties/<name>_isobar_1atm.tsv). This row holds 870, which is
    #   above that C9-C12 span.
    #   They were fetched to re-point rows like this one, and reading them settled
    #   that plan against itself: pointing kerosene at decane moves it 810 -> 730,
    #   an 11% error, because a real cut carries aromatics and cycloalkanes a pure
    #   n-alkane does not. Sourcing a row must not change what the row names.
    #   C3.5 residual - searched for, and the search is recorded.
    "Crude Oil": 870,           # Typical average value
    # [UNVERIFIED] a mixture or commercial grade with no single composition, so no
    #   pure substance on disk IS this row. The nearest pure n-alkanes are now on
    #   disk and bracket it: nonane 718.03, decane 730.41, dodecane 749.44 kg/m3 at 293.15 K, 1 atm
    #   (nist_fluid_properties/<name>_isobar_1atm.tsv). This row holds 810, which is
    #   above that C9-C12 span.
    #   They were fetched to re-point rows like this one, and reading them settled
    #   that plan against itself: pointing kerosene at decane moves it 810 -> 730,
    #   an 11% error, because a real cut carries aromatics and cycloalkanes a pure
    #   n-alkane does not. Sourcing a row must not change what the row names.
    #   C3.5 residual - searched for, and the search is recorded.
    "Kerosene": 810,
    # [UNVERIFIED] a mixture or commercial grade with no single composition, so no
    #   pure substance on disk IS this row. The nearest pure n-alkanes are now on
    #   disk and bracket it: nonane 718.03, decane 730.41, dodecane 749.44 kg/m3 at 293.15 K, 1 atm
    #   (nist_fluid_properties/<name>_isobar_1atm.tsv). This row holds 850, which is
    #   above that C9-C12 span.
    #   They were fetched to re-point rows like this one, and reading them settled
    #   that plan against itself: pointing kerosene at decane moves it 810 -> 730,
    #   an 11% error, because a real cut carries aromatics and cycloalkanes a pure
    #   n-alkane does not. Sourcing a row must not change what the row names.
    #   C3.5 residual - searched for, and the search is recorded.
    "Diesel Fuel": 850,
    # [UNVERIFIED] a mixture or commercial grade with no single composition, so no
    #   pure substance on disk IS this row. The nearest pure n-alkanes are now on
    #   disk and bracket it: nonane 718.03, decane 730.41, dodecane 749.44 kg/m3 at 293.15 K, 1 atm
    #   (nist_fluid_properties/<name>_isobar_1atm.tsv). This row holds 726, which is
    #   inside that C9-C12 span - COINCIDENTALLY, and that is NOT support.
    #   A petrol cut is C4-C12 with most of its mass in C5-C8, LIGHTER than every
    #   alkane in this bracket. Falling between nonane and decane is one number
    #   agreeing with a span of heavier pure substances; it says nothing about a
    #   mixture whose composition sits below them.
    #   They were fetched to re-point rows like this one, and reading them settled
    #   that plan against itself: pointing kerosene at decane moves it 810 -> 730,
    #   an 11% error, because a real cut carries aromatics and cycloalkanes a pure
    #   n-alkane does not. Sourcing a row must not change what the row names.
    #   C3.5 residual - searched for, and the search is recorded.
    "Gasoline": 726,
    # [UNVERIFIED] a mixture or commercial grade with no single composition, so no
    #   pure substance on disk IS this row. The nearest pure n-alkanes are now on
    #   disk and bracket it: nonane 718.03, decane 730.41, dodecane 749.44 kg/m3 at 293.15 K, 1 atm
    #   (nist_fluid_properties/<name>_isobar_1atm.tsv). This row holds 888, which is
    #   above that C9-C12 span.
    #   They were fetched to re-point rows like this one, and reading them settled
    #   that plan against itself: pointing kerosene at decane moves it 810 -> 730,
    #   an 11% error, because a real cut carries aromatics and cycloalkanes a pure
    #   n-alkane does not. Sourcing a row must not change what the row names.
    #   C3.5 residual - searched for, and the search is recorded.
    "Engine Oil": 888,
    # [UNVERIFIED] a mixture, solution or commercial product with no single composition in the row;
    #   no on-disk source, and none was searched for in C3. Residual register (C3.5).
    "Hydraulic Oil": 850,       # Typical average value
    
    # Alcohols and Solvents
    # [ON-DISK] pubchem/ethanol_density.json @ cid=702 heading="Density" ref=79 scale=1000 precision=3sf
    #   0.7893 g/cu cm at 20 °C
    #   Source: Hazardous Substances Data Bank (HSDB); conditions: 20 C.
    #   That entry cites: Haynes, W.M. (ed.). CRC Handbook of Chemistry and Physics. 95th Edition. CRC Press LLC, Boca Raton: 
    #   PubChem is an AGGREGATOR carrying several disagreeing densities, so the
    #   citation names ONE ReferenceNumber rather than the CID alone.
    "Ethyl Alcohol (Ethanol)": 789,
    # [ON-DISK] nist_fluid_properties/methanol_C67561_isobar_1atm.tsv @ T=293.15 col="Density (kg/m3)" precision=3sf
    #   791.01 at 293.15 K.
    "Methyl Alcohol (Methanol)": 791,
    # [UNVERIFIED] Isopropanol is not in the NIST fluid database on disk; its WebBook species pages
    #   on disk (isopropanol_C67630_condensed_phase.html, isopropanol_C67630_phase_change.html) say
    #   density only as ['Critical density']. Residual register (C3.5).
    "Isopropyl Alcohol (IPA)": 786,
    # [ON-DISK] pubchem/acetone_density.json @ cid=180 heading="Density" ref=50 scale=1000 precision=exact
    #   0.791 at 68 °F (USCG, 1999) - Less dense than water; will float
    #   Source: CAMEO Chemicals; conditions: 68 F.
    #   That entry cites: U.S. Coast Guard. 1999. Chemical Hazard Response Information System (CHRIS) - Hazardous Chemical Dat
    #   PubChem is an AGGREGATOR carrying several disagreeing densities, so the
    #   citation names ONE ReferenceNumber rather than the CID alone.
    "Acetone": 791,
    # [UNVERIFIED] a mixture, solution or commercial product with no single composition in the row;
    #   no on-disk source, and none was searched for in C3. Residual register (C3.5).
    "Turpentine": 870,
    # [ON-DISK] nist_fluid_properties/benzene_C71432_isobar_1atm.tsv @ T=293.15 col="Density (kg/m3)" precision=3sf
    #   878.92 at 293.15 K, 1 atm - the conditions this table states.
    #   Corrected from 876 (-0.33%, not a 3-s.f. rounding); P6 event.
    "Benzene": 879,
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
    # [ON-DISK] nist_fluid_properties/propane_C74986_saturation_298.15K.tsv @ T=298.15 col="Density (l, kg/m3)" tol=0.13% basis=condition
    #   492.36 at 298.15 K; the row is +0.13%. The row says "At 25 C, under pressure" and names no pressure; the file is the saturated liquid. Pressure in that row: 0.95207 MPa.
    #   C3 review, H (electrical) F-2 generalised: this tol= is the residual rounded up
    #   (0.12999% -> 0.13%) and is the only one outside the electrical branch. The
    #   artefact carries five significant figures, so half a unit in its last digit is
    #   0.001% and cannot be the basis. Half a unit in the ROW's own last digit is
    #   0.5/493 = 0.101%, which would NOT admit the value - so the warrant has to be
    #   the unstated pressure, and nothing here bounds it. Referred (C3.5).
    "Liquid Propane": 493,      # At 25°C, under pressure
    
    # High-Density Liquids
    # [ON-DISK] pubchem/mercury_density.json @ cid=23931 heading="Density" ref=8 scale=1000 precision=exact
    #   13.55 at 68 °F (USCG, 1999) - Denser than water; will sink
    #   Source: CAMEO Chemicals; conditions: 68 F.
    #   That entry cites: U.S. Coast Guard. 1999. Chemical Hazard Response Information System (CHRIS) - Hazardous Chemical Dat
    #   PubChem is an AGGREGATOR and this record carries several densities that
    #   disagree, so the citation names ONE ReferenceNumber rather than the CID alone.
    "Mercury": 13550,
    # [ON-DISK] pubchem/glycerol_density.json @ cid=753 heading="Density" ref=11 scale=1000 precision=3sf
    #   1.261 at 68 °F (USCG, 1999) - Denser than water; will sink
    #   Source: CAMEO Chemicals; conditions: 68 F.
    #   That entry cites: U.S. Coast Guard. 1999. Chemical Hazard Response Information System (CHRIS) - Hazardous Chemical Dat
    #   PubChem is an AGGREGATOR and this record carries several densities that
    #   disagree, so the citation names ONE ReferenceNumber rather than the CID alone.
    "Glycerin": 1260,
    # [UNVERIFIED] Chloroform is not in the NIST fluid database on disk (MANIFEST: 36 fluids), and
    #   no WebBook species page for it is on disk. Residual register (C3.5).
    "Chloroform": 1489,
    # [ON-DISK] pubchem/carbon_tetrachloride_density.json @ cid=5943 heading="Density" ref=76 scale=1000 precision=exact
    #   1.5940 @ 20°C
    #   Source: PAC Chemical Database, U.S. Department of Energy; conditions: 20 C.
    #   PubChem is an AGGREGATOR and this record carries several densities that
    #   disagree, so the citation names ONE ReferenceNumber rather than the CID alone.
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
    # [ON-DISK] nist_fluid_properties/r134a_C811972_saturation_298.15K.tsv @ T=298.15 col="Density (l, kg/m3)" precision=4sf
    #   1206.7 at 298.15 K, saturated liquid - what the row states and what the file is.
    #   Pressure in that row: 0.66538 MPa. Corrected from 1206 (-0.06%); P6 event.
    "R-134a Refrigerant (Saturated Liquid at 25°C)": 1207,
    
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
    # [UNVERIFIED] the row names a material FPL indexes by SPECIES and states no
    #   moisture condition, and both are needed before a density can be cited.
    #   FPL-GTR-282 pp.123-136 index several species under "Pine".
    #   Specific gravities READ from the metric tables (5-3a / 5-5a; 5-3b and 5-5b
    #   repeat the same species in inch-pound and are not counted twice):
    #     green 0.34-0.72 over 14 row(s)
    #     12%   0.35-0.56 over 4 row(s)
    #   Table 4-6a (p.107) converts those at 12% MC to 381-627 kg/m3
    #   (Gx 0.34-0.56); the row holds 500.
    #   THAT SPAN IS A LOWER BOUND, not FPL's range. FPL prints the genus heading
    #   once and then bare sub-names ("Bigleaf", "Black", "Red"), so a scan anchored
    #   on the genus truncates the block - for Maple it stops before Silver and
    #   Sugar. No claim is made here that any committed value lies outside FPL.
    #   C3.5 residual.
    "Pine Wood": 500,
    # [UNVERIFIED] the row names a material FPL indexes by SPECIES and states no
    #   moisture condition, and both are needed before a density can be cited.
    #   FPL-GTR-282 pp.121-135 index several species under "Oak".
    #   Specific gravities READ from the metric tables (5-3a / 5-5a; 5-3b and 5-5b
    #   repeat the same species in inch-pound and are not counted twice):
    #     green 0.56-0.80 over 9 row(s)
    #     12%   0.61-0.69 over 6 row(s)
    #   Table 4-6a (p.107) converts those at 12% MC to 672-762 kg/m3
    #   (Gx 0.60-0.68); the row holds 750.
    #   THAT SPAN IS A LOWER BOUND, not FPL's range. FPL prints the genus heading
    #   once and then bare sub-names ("Bigleaf", "Black", "Red"), so a scan anchored
    #   on the genus truncates the block - for Maple it stops before Silver and
    #   Sugar. No claim is made here that any committed value lies outside FPL.
    #   C3.5 residual.
    "Oak Wood": 750,
    # [UNVERIFIED] the row names a material FPL indexes by SPECIES and states no
    #   moisture condition, and both are needed before a density can be cited.
    #   FPL-GTR-282 pp.137-137 index several species under "Teak".
    #   Specific gravities READ from the metric tables (5-3a / 5-5a; 5-3b and 5-5b
    #   repeat the same species in inch-pound and are not counted twice):
    #     green 0.45-0.78 over 3 row(s)
    #   Table 5-5a prints no 12% specific gravity for it, only green.
    #   THAT SPAN IS A LOWER BOUND, not FPL's range. FPL prints the genus heading
    #   once and then bare sub-names ("Bigleaf", "Black", "Red"), so a scan anchored
    #   on the genus truncates the block - for Maple it stops before Silver and
    #   Sugar. No claim is made here that any committed value lies outside FPL.
    #   C3.5 residual.
    "Teak Wood": 630,
    # [UNVERIFIED] the row names a material FPL indexes by SPECIES and states no
    #   moisture condition, and both are needed before a density can be cited.
    #   FPL-GTR-282 pp.121-121 index several species under "Maple".
    #   Specific gravities READ from the metric tables (5-3a / 5-5a; 5-3b and 5-5b
    #   repeat the same species in inch-pound and are not counted twice):
    #     green 0.44-0.52 over 3 row(s)
    #     12%   0.48-0.57 over 2 row(s)
    #   Table 4-6a (p.107) converts those at 12% MC to 538-627 kg/m3
    #   (Gx 0.48-0.56); the row holds 740.
    #   THAT SPAN IS A LOWER BOUND, not FPL's range. FPL prints the genus heading
    #   once and then bare sub-names ("Bigleaf", "Black", "Red"), so a scan anchored
    #   on the genus truncates the block - for Maple it stops before Silver and
    #   Sugar. No claim is made here that any committed value lies outside FPL.
    #   C3.5 residual.
    "Maple Wood": 740,
    # [UNVERIFIED] the row names a material FPL indexes by SPECIES and states no
    #   moisture condition, and both are needed before a density can be cited.
    #   FPL-GTR-282 prints NO specific gravity under the name "Ebony": it is not
    #   indexed in Table 5-3a (US species, PDF pp.120-124) or Table 5-5a (imports,
    #   pp.134-137). The committed 1200 kg/m3 has no support on disk and no
    #   substitute was found - a different source is needed, not a better reading
    #   of this one. C3.5 residual.
    "Ebony Wood": 1200,  # Sinks in water
    # [UNVERIFIED] no source on disk for this material. Residual register (C3.5).
    "Rubber (Natural)": 950,
    
    # Plastics & Synthetics (range from floating to sinking)
    # [UNVERIFIED] SEARCHED, and found unreachable by a chemical database. PubChem
    #   returned PUGREST.NotFound for every polymer looked up (polyethylene,
    #   polypropylene, polystyrene, PVC, PTFE, PMMA): a polymer is not a single
    #   compound with a CID. That is the same reason its density depends on
    #   crystallinity and grade rather than on a formula, so no synonym would fix it.
    #   This row needs a materials handbook, not another fetch. C3.5.
    "Polypropylene (PP)": 900,
    # [UNVERIFIED] SEARCHED, and found unreachable by a chemical database. PubChem
    #   returned PUGREST.NotFound for every polymer looked up (polyethylene,
    #   polypropylene, polystyrene, PVC, PTFE, PMMA): a polymer is not a single
    #   compound with a CID. That is the same reason its density depends on
    #   crystallinity and grade rather than on a formula, so no synonym would fix it.
    #   This row needs a materials handbook, not another fetch. C3.5.
    "Polyethylene (LDPE)": 920,
    # [UNVERIFIED] SEARCHED, and found unreachable by a chemical database. PubChem
    #   returned PUGREST.NotFound for every polymer looked up (polyethylene,
    #   polypropylene, polystyrene, PVC, PTFE, PMMA): a polymer is not a single
    #   compound with a CID. That is the same reason its density depends on
    #   crystallinity and grade rather than on a formula, so no synonym would fix it.
    #   This row needs a materials handbook, not another fetch. C3.5.
    "Polyethylene (HDPE)": 950,
    # [UNVERIFIED] SEARCHED, and found unreachable by a chemical database. PubChem
    #   returned PUGREST.NotFound for every polymer looked up (polyethylene,
    #   polypropylene, polystyrene, PVC, PTFE, PMMA): a polymer is not a single
    #   compound with a CID. That is the same reason its density depends on
    #   crystallinity and grade rather than on a formula, so no synonym would fix it.
    #   This row needs a materials handbook, not another fetch. C3.5.
    "Nylon": 1150,
    # [UNVERIFIED] SEARCHED, and found unreachable by a chemical database. PubChem
    #   returned PUGREST.NotFound for every polymer looked up (polyethylene,
    #   polypropylene, polystyrene, PVC, PTFE, PMMA): a polymer is not a single
    #   compound with a CID. That is the same reason its density depends on
    #   crystallinity and grade rather than on a formula, so no synonym would fix it.
    #   This row needs a materials handbook, not another fetch. C3.5.
    "Polystyrene (PS)": 1050,
    # [UNVERIFIED] SEARCHED, and found unreachable by a chemical database. PubChem
    #   returned PUGREST.NotFound for every polymer looked up (polyethylene,
    #   polypropylene, polystyrene, PVC, PTFE, PMMA): a polymer is not a single
    #   compound with a CID. That is the same reason its density depends on
    #   crystallinity and grade rather than on a formula, so no synonym would fix it.
    #   This row needs a materials handbook, not another fetch. C3.5.
    "PVC (Rigid)": 1380,
    # [UNVERIFIED] SEARCHED, and found unreachable by a chemical database. PubChem
    #   returned PUGREST.NotFound for every polymer looked up (polyethylene,
    #   polypropylene, polystyrene, PVC, PTFE, PMMA): a polymer is not a single
    #   compound with a CID. That is the same reason its density depends on
    #   crystallinity and grade rather than on a formula, so no synonym would fix it.
    #   This row needs a materials handbook, not another fetch. C3.5.
    "PTFE (Teflon)": 2200,
    # [UNVERIFIED] SEARCHED, and found unreachable by a chemical database. PubChem
    #   returned PUGREST.NotFound for every polymer looked up (polyethylene,
    #   polypropylene, polystyrene, PVC, PTFE, PMMA): a polymer is not a single
    #   compound with a CID. That is the same reason its density depends on
    #   crystallinity and grade rather than on a formula, so no synonym would fix it.
    #   This row needs a materials handbook, not another fetch. C3.5.
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
    # [UNVERIFIED] the row names NO ALLOY, and the three readings available disagree by
    #   more than any of them is precise. PubChem cid=5359268 prints 2.7 g/cm3
    #   (pure aluminium, 2700); MIL-HDBK-5J prints 0.098 lb/in^3 for 6061 (p.566)
    #   = 2713 kg/m^3 and 0.100 for 2024 (p.373) = 2768. The row holds 2710, between the
    #   alloys and nearest 6061.
    #   NOT corrected: replacing it with the pure-metal figure would be CHOOSING which
    #   aluminium this row means and calling it sourcing. Naming the alloy in the key is
    #   the fix, and that moves the item pool. C3.5 residual.
    "Aluminum": 2710,
    # [ON-DISK] mil_hdbk_5j/MIL-HDBK-5J_2003-01-31.pdf @ page=899 text="Table 5.2.1.0(b)" mil="density" scale=27679.9 tol=0.31% basis=half-unit
    #   MIL-HDBK-5J p.899, Table 5.2.1.0(b) (Commercially Pure Titanium): density 0.163 lb/in^3 = 4512 kg/m3. The row is -0.26% from it, inside
    #   half a unit in the handbook's last digit (0.31%) - the bound a unit conversion of a
    #   3-s.f. value supports (SP 811 p.66: 1 lb/in^3 = 27679.9 kg/m^3). Commercially pure; the row names no grade - the implementer's choice, referred (C3.5).
    "Titanium": 4500,
    # [UNVERIFIED] no source on disk for this material. Residual register (C3.5).
    "Zinc": 7140,
    # [ON-DISK] pubchem/tin_density.json @ cid=5352426 heading="Density" ref=6 scale=1000 precision=exact
    #   7.28 (NIOSH, 2024) - Denser than water; will sink
    #   Source: CAMEO Chemicals; conditions: no temperature stated in the entry.
    #   PubChem is an AGGREGATOR and this record carries several densities that
    #   disagree, so the citation names ONE ReferenceNumber rather than the CID alone.
    "Tin": 7280,
    # [UNVERIFIED] no source on disk for this material. Residual register (C3.5).
    "Iron (Wrought)": 7750,
    # [ON-DISK] mil_hdbk_5j/MIL-HDBK-5J_2003-01-31.pdf @ page=62 text="Table 2.2.1.0(b)" mil="density" scale=27679.9 tol=0.18% basis=half-unit
    #   MIL-HDBK-5J p.62, Table 2.2.1.0(b) (AISI 1025 Carbon Steel): density 0.284 lb/in^3 = 7861 kg/m3. The row is -0.14% from it, inside
    #   half a unit in the handbook's last digit (0.18%) - the bound a unit conversion of a
    #   3-s.f. value supports (SP 811 p.66: 1 lb/in^3 = 27679.9 kg/m^3). AISI 1025, the handbook's carbon-steel table; the row names no grade, so the grade is the implementer's choice - referred (C3.5).
    "Steel (Carbon)": 7850,
    # [UNVERIFIED] no source on disk for this material. Residual register (C3.5).
    "Brass": 8600,
    # [ON-DISK] pubchem/copper_density.json @ cid=23978 heading="Density" ref=46 scale=1000 precision=exact
    #   8.94
    #   Source: Hazardous Substances Data Bank (HSDB); conditions: no temperature stated in the entry.
    #   That entry cites: Budavari, S. (ed.). The Merck Index - An Encyclopedia of Chemicals, Drugs, and Biologicals. Whitehou
    #   PubChem is an AGGREGATOR and this record carries several densities that
    #   disagree, so the citation names ONE ReferenceNumber rather than the CID alone.
    "Copper": 8940,
    # [ON-DISK] pubchem/nickel_density.json @ cid=935 heading="Density" ref=50 scale=1000 precision=2sf
    #   8.908
    #   Source: Hazardous Substances Data Bank (HSDB); conditions: no temperature stated in the entry.
    #   That entry cites: Lewis, R.J., Sr (Ed.). Hawley's Condensed Chemical  Dictionary. 13th ed. New York, NY: John Wiley & 
    #   PubChem is an AGGREGATOR and this record carries several densities that
    #   disagree, so the citation names ONE ReferenceNumber rather than the CID alone.
    "Nickel": 8900,
    # [ON-DISK] pubchem/silver_density.json @ cid=23954 heading="Density" ref=47 scale=1000 precision=3sf
    #   10.49 @ 15 °C
    #   Source: Hazardous Substances Data Bank (HSDB); conditions: 15 C.
    #   That entry cites: O'Neil, M.J. (ed.). The Merck Index - An Encyclopedia of  Chemicals, Drugs, and Biologicals. 13th Ed
    #   PubChem is an AGGREGATOR and this record carries several densities that
    #   disagree, so the citation names ONE ReferenceNumber rather than the CID alone.
    "Silver": 10500,
    # [ON-DISK] pubchem/lead_density.json @ cid=5352425 heading="Density" ref=10 scale=1000 precision=4sf
    #   11.3437 at 61 °F (NTP, 1992) - Denser than water; will sink
    #   Source: CAMEO Chemicals; conditions: 61 F.
    #   That entry cites: National Toxicology Program, Institute of Environmental Health Sciences, National Institutes of Heal
    #   PubChem is an AGGREGATOR and this record carries several densities that
    #   disagree, so the citation names ONE ReferenceNumber rather than the CID alone.
    "Lead": 11340,
    # [ON-DISK] pubchem/uranium_density.json @ cid=23989 heading="Density" ref=41 scale=1000 precision=exact
    #   19.1 @25 °C
    #   Source: PAC Chemical Database, U.S. Department of Energy; conditions: 25 C.
    #   PubChem is an AGGREGATOR and this record carries several densities that
    #   disagree, so the citation names ONE ReferenceNumber rather than the CID alone.
    "Uranium": 19100,
    # [ON-DISK] pubchem/gold_density.json @ cid=23985 heading="Density" ref=34 scale=1000 precision=exact
    #   19.3
    #   Source: Hazardous Substances Data Bank (HSDB); conditions: no temperature stated in the entry.
    #   That entry cites: O'Neil, M.J. (ed.). The Merck Index - An Encyclopedia of  Chemicals, Drugs, and Biologicals. 13th Ed
    #   PubChem is an AGGREGATOR and this record carries several densities that
    #   disagree, so the citation names ONE ReferenceNumber rather than the CID alone.
    "Gold": 19300,
    # [ON-DISK] pubchem/tungsten_density.json @ cid=23964 heading="Density" ref=46 scale=1000 precision=exact
    #   19.3 @25 °C
    #   Source: PAC Chemical Database, U.S. Department of Energy; conditions: 25 C.
    #   Corrected from 19600 (+1.55%). The previous tag said "no source on
    #   disk for this material"; the record is on disk and 5 independent entries in it
    #   print 19.3 - CAMEO Chemicals, ILO-WHO International , Occupational Safety an, PAC Chemical Database,, The National Institute.
    #   HSDB's entry is "18.7-19.3 @ 20 C/4 C; depends on extent of working" - a RANGE,
    #   which pubchem.py refuses as a value. Reading its leading number instead would
    #   make tungsten 18.7, which is how a range becomes a wrong value. P6 event.
    "Tungsten": 19300,
    # [ON-DISK] pubchem/platinum_density.json @ cid=23939 heading="Density" ref=28 scale=1000 precision=4sf
    #   21.447 (calc)
    #   Source: Hazardous Substances Data Bank (HSDB); conditions: no temperature stated in the entry.
    #   That entry cites: Budavari, S. (ed.). The Merck Index - An Encyclopedia of Chemicals, Drugs, and Biologicals. Whitehou
    #   PubChem is an AGGREGATOR and this record carries several densities that
    #   disagree, so the citation names ONE ReferenceNumber rather than the CID alone.
    "Platinum": 21450,
    # [ON-DISK] pubchem/osmium_density.json @ cid=23937 heading="Density" ref=22 scale=1000 precision=exact
    #   22.59 g/cu cm
    #   Source: Hazardous Substances Data Bank (HSDB); conditions: no temperature stated in the entry.
    #   That entry cites: Lide, DR (ed.). CRC Handbook of Chemistry and Physics. 81st Edition. CRC Press LLC, Boca Raton: FL 2
    #   PubChem is an AGGREGATOR and this record carries several densities that
    #   disagree, so the citation names ONE ReferenceNumber rather than the CID alone.
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
    # [UNVERIFIED] SEARCHED, and found unreachable by a chemical database.
    #   pubchem/silicon_dioxide_density.json (cid=24261) is on disk and gives
    #   2200, 2300, 2334, 2600 kg/m3 - amorphous silica, fumed silica and silica gel.
    #   None of those is crystalline alpha-quartz, which is what this row names and
    #   what 2650 is. The record is the wrong FORM of the right substance, so it
    #   neither supports nor challenges the committed value.
    #   A mineral density needs a materials or mineralogy source, not PubChem. C3.5.
    "Quartz": 2650,
    # [UNVERIFIED] no source on disk for this material. Residual register (C3.5).
    "Granite": 2700,
    # [UNVERIFIED] no source on disk for this material. Residual register (C3.5).
    "Marble": 2700,
    # [UNVERIFIED] no source on disk for this material. Residual register (C3.5).
    "Sandstone": 2300,
    # [UNVERIFIED] no source on disk for this material. Residual register (C3.5).
    "Limestone": 2500,
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
    # [UNVERIFIED] a mixture or commercial grade with no single composition, so no
    #   pure substance on disk IS this row. The nearest pure n-alkanes are now on
    #   disk and bracket it: nonane 718.03, decane 730.41, dodecane 749.44 kg/m3 at 293.15 K, 1 atm
    #   (nist_fluid_properties/<name>_isobar_1atm.tsv). This row holds 726, which is
    #   inside that C9-C12 span - COINCIDENTALLY, and that is NOT support.
    #   A petrol cut is C4-C12 with most of its mass in C5-C8, LIGHTER than every
    #   alkane in this bracket. Falling between nonane and decane is one number
    #   agreeing with a span of heavier pure substances; it says nothing about a
    #   mixture whose composition sits below them.
    #   They were fetched to re-point rows like this one, and reading them settled
    #   that plan against itself: pointing kerosene at decane moves it 810 -> 730,
    #   an 11% error, because a real cut carries aromatics and cycloalkanes a pure
    #   n-alkane does not. Sourcing a row must not change what the row names.
    #   C3.5 residual - searched for, and the search is recorded.
    "Gasoline": 726,
    # [UNVERIFIED] a mixture or commercial grade with no single composition, so no
    #   pure substance on disk IS this row. The nearest pure n-alkanes are now on
    #   disk and bracket it: nonane 718.03, decane 730.41, dodecane 749.44 kg/m3 at 293.15 K, 1 atm
    #   (nist_fluid_properties/<name>_isobar_1atm.tsv). This row holds 810, which is
    #   above that C9-C12 span.
    #   They were fetched to re-point rows like this one, and reading them settled
    #   that plan against itself: pointing kerosene at decane moves it 810 -> 730,
    #   an 11% error, because a real cut carries aromatics and cycloalkanes a pure
    #   n-alkane does not. Sourcing a row must not change what the row names.
    #   C3.5 residual - searched for, and the search is recorded.
    "Kerosene": 810,
    # [UNVERIFIED] no source on disk. Residual register (C3.5).
    "Diesel Fuel": 850,
    # [UNVERIFIED] no source on disk. Residual register (C3.5).
    "Jet Fuel (JP-4)": 770,
    # [ON-DISK] pubchem/ethanol_density.json @ cid=702 heading="Density" ref=79 scale=1000 precision=3sf
    #   0.7893 g/cu cm at 20 °C
    #   Source: Hazardous Substances Data Bank (HSDB); conditions: 20 C.
    #   That entry cites: Haynes, W.M. (ed.). CRC Handbook of Chemistry and Physics. 95th Edition. CRC Press LLC, Boca Raton: 
    #   PubChem is an AGGREGATOR carrying several disagreeing densities, so the
    #   citation names ONE ReferenceNumber rather than the CID alone.
    "Ethanol": 789,
    # [UNVERIFIED] the table states no temperature or pressure, so no row of an artefact is THE
    #   value. NIST at 1 atm: 809.73 (liquid) at 273.15 K; 791.01 (liquid) at 293.15 K
    #   (nist_fluid_properties/methanol_C67561_isobar_1atm.tsv). Residual register (C3.5).
    "Methanol": 791,
    # [UNVERIFIED] no source on disk. Residual register (C3.5).
    "Isopropyl Alcohol": 786,
    
    # Oils & Lubricants
    # [UNVERIFIED] a mixture or commercial grade with no single composition, so no
    #   pure substance on disk IS this row. The nearest pure n-alkanes are now on
    #   disk and bracket it: nonane 718.03, decane 730.41, dodecane 749.44 kg/m3 at 293.15 K, 1 atm
    #   (nist_fluid_properties/<name>_isobar_1atm.tsv). This row holds 870, which is
    #   above that C9-C12 span.
    #   They were fetched to re-point rows like this one, and reading them settled
    #   that plan against itself: pointing kerosene at decane moves it 810 -> 730,
    #   an 11% error, because a real cut carries aromatics and cycloalkanes a pure
    #   n-alkane does not. Sourcing a row must not change what the row names.
    #   C3.5 residual - searched for, and the search is recorded.
    "SAE 10 Oil": 870,
    # [UNVERIFIED] a mixture or commercial grade with no single composition, so no
    #   pure substance on disk IS this row. The nearest pure n-alkanes are now on
    #   disk and bracket it: nonane 718.03, decane 730.41, dodecane 749.44 kg/m3 at 293.15 K, 1 atm
    #   (nist_fluid_properties/<name>_isobar_1atm.tsv). This row holds 880, which is
    #   above that C9-C12 span.
    #   They were fetched to re-point rows like this one, and reading them settled
    #   that plan against itself: pointing kerosene at decane moves it 810 -> 730,
    #   an 11% error, because a real cut carries aromatics and cycloalkanes a pure
    #   n-alkane does not. Sourcing a row must not change what the row names.
    #   C3.5 residual - searched for, and the search is recorded.
    "SAE 20 Oil": 880,
    # [UNVERIFIED] a mixture or commercial grade with no single composition, so no
    #   pure substance on disk IS this row. The nearest pure n-alkanes are now on
    #   disk and bracket it: nonane 718.03, decane 730.41, dodecane 749.44 kg/m3 at 293.15 K, 1 atm
    #   (nist_fluid_properties/<name>_isobar_1atm.tsv). This row holds 917, which is
    #   above that C9-C12 span.
    #   They were fetched to re-point rows like this one, and reading them settled
    #   that plan against itself: pointing kerosene at decane moves it 810 -> 730,
    #   an 11% error, because a real cut carries aromatics and cycloalkanes a pure
    #   n-alkane does not. Sourcing a row must not change what the row names.
    #   C3.5 residual - searched for, and the search is recorded.
    "SAE 30 Oil": 917,
    # [UNVERIFIED] a mixture or commercial grade with no single composition, so no
    #   pure substance on disk IS this row. The nearest pure n-alkanes are now on
    #   disk and bracket it: nonane 718.03, decane 730.41, dodecane 749.44 kg/m3 at 293.15 K, 1 atm
    #   (nist_fluid_properties/<name>_isobar_1atm.tsv). This row holds 945, which is
    #   above that C9-C12 span.
    #   They were fetched to re-point rows like this one, and reading them settled
    #   that plan against itself: pointing kerosene at decane moves it 810 -> 730,
    #   an 11% error, because a real cut carries aromatics and cycloalkanes a pure
    #   n-alkane does not. Sourcing a row must not change what the row names.
    #   C3.5 residual - searched for, and the search is recorded.
    "SAE 40 Oil": 945,
    # [UNVERIFIED] a mixture or commercial grade with no single composition, so no
    #   pure substance on disk IS this row. The nearest pure n-alkanes are now on
    #   disk and bracket it: nonane 718.03, decane 730.41, dodecane 749.44 kg/m3 at 293.15 K, 1 atm
    #   (nist_fluid_properties/<name>_isobar_1atm.tsv). This row holds 960, which is
    #   above that C9-C12 span.
    #   They were fetched to re-point rows like this one, and reading them settled
    #   that plan against itself: pointing kerosene at decane moves it 810 -> 730,
    #   an 11% error, because a real cut carries aromatics and cycloalkanes a pure
    #   n-alkane does not. Sourcing a row must not change what the row names.
    #   C3.5 residual - searched for, and the search is recorded.
    "SAE 50 Oil": 960,
    # [UNVERIFIED] a mixture or commercial grade with no single composition, so no
    #   pure substance on disk IS this row. The nearest pure n-alkanes are now on
    #   disk and bracket it: nonane 718.03, decane 730.41, dodecane 749.44 kg/m3 at 293.15 K, 1 atm
    #   (nist_fluid_properties/<name>_isobar_1atm.tsv). This row holds 825, which is
    #   above that C9-C12 span.
    #   They were fetched to re-point rows like this one, and reading them settled
    #   that plan against itself: pointing kerosene at decane moves it 810 -> 730,
    #   an 11% error, because a real cut carries aromatics and cycloalkanes a pure
    #   n-alkane does not. Sourcing a row must not change what the row names.
    #   C3.5 residual - searched for, and the search is recorded.
    "Crude Oil (Light)": 825,
    # [UNVERIFIED] a mixture or commercial grade with no single composition, so no
    #   pure substance on disk IS this row. The nearest pure n-alkanes are now on
    #   disk and bracket it: nonane 718.03, decane 730.41, dodecane 749.44 kg/m3 at 293.15 K, 1 atm
    #   (nist_fluid_properties/<name>_isobar_1atm.tsv). This row holds 950, which is
    #   above that C9-C12 span.
    #   They were fetched to re-point rows like this one, and reading them settled
    #   that plan against itself: pointing kerosene at decane moves it 810 -> 730,
    #   an 11% error, because a real cut carries aromatics and cycloalkanes a pure
    #   n-alkane does not. Sourcing a row must not change what the row names.
    #   C3.5 residual - searched for, and the search is recorded.
    "Crude Oil (Heavy)": 950,
    # [UNVERIFIED] a mixture or commercial grade with no single composition, so no
    #   pure substance on disk IS this row. The nearest pure n-alkanes are now on
    #   disk and bracket it: nonane 718.03, decane 730.41, dodecane 749.44 kg/m3 at 293.15 K, 1 atm
    #   (nist_fluid_properties/<name>_isobar_1atm.tsv). This row holds 888, which is
    #   above that C9-C12 span.
    #   They were fetched to re-point rows like this one, and reading them settled
    #   that plan against itself: pointing kerosene at decane moves it 810 -> 730,
    #   an 11% error, because a real cut carries aromatics and cycloalkanes a pure
    #   n-alkane does not. Sourcing a row must not change what the row names.
    #   C3.5 residual - searched for, and the search is recorded.
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
    # [ON-DISK] pubchem/ethylene_glycol_density.json @ cid=174 heading="Density" ref=124 scale=1000 precision=3sf
    #   1.1135 @ 20°C
    #   Source: PAC Chemical Database, U.S. Department of Energy; conditions: 20 C.
    #   PubChem is an AGGREGATOR carrying several disagreeing densities, so the
    #   citation names ONE ReferenceNumber rather than the CID alone.
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
    # [ON-DISK] pubchem/mercury_density.json @ cid=23931 heading="Density" ref=8 scale=1000 precision=exact
    #   13.55 at 68 °F (USCG, 1999) - Denser than water; will sink
    #   Source: CAMEO Chemicals; conditions: 68 F.
    #   That entry cites: U.S. Coast Guard. 1999. Chemical Hazard Response Information System (CHRIS) - Hazardous Chemical Dat
    #   PubChem is an AGGREGATOR and this record carries several densities that
    #   disagree, so the citation names ONE ReferenceNumber rather than the CID alone.
    "Mercury": 13550,
    # [UNVERIFIED] the table states no temperature or pressure, so no row of an artefact is THE
    #   value. NIST at 1 atm: no row at 273.15 K; 998.21 (liquid) at 293.15 K
    #   (nist_fluid_properties/water_C7732185_isobar_1atm.tsv). Residual register (C3.5).
    "Water": 1000,  # For gas measurements
    # [UNVERIFIED] no source on disk. Residual register (C3.5).
    "Sea Water": 1025,
    
    # Heavy Oils & Lubricants
    # [UNVERIFIED] a mixture or commercial grade with no single composition, so no
    #   pure substance on disk IS this row. The nearest pure n-alkanes are now on
    #   disk and bracket it: nonane 718.03, decane 730.41, dodecane 749.44 kg/m3 at 293.15 K, 1 atm
    #   (nist_fluid_properties/<name>_isobar_1atm.tsv). This row holds 960, which is
    #   above that C9-C12 span.
    #   They were fetched to re-point rows like this one, and reading them settled
    #   that plan against itself: pointing kerosene at decane moves it 810 -> 730,
    #   an 11% error, because a real cut carries aromatics and cycloalkanes a pure
    #   n-alkane does not. Sourcing a row must not change what the row names.
    #   C3.5 residual - searched for, and the search is recorded.
    "SAE 50 Oil": 960,
    # [UNVERIFIED] a mixture or commercial grade with no single composition, so no
    #   pure substance on disk IS this row. The nearest pure n-alkanes are now on
    #   disk and bracket it: nonane 718.03, decane 730.41, dodecane 749.44 kg/m3 at 293.15 K, 1 atm
    #   (nist_fluid_properties/<name>_isobar_1atm.tsv). This row holds 920, which is
    #   above that C9-C12 span.
    #   They were fetched to re-point rows like this one, and reading them settled
    #   that plan against itself: pointing kerosene at decane moves it 810 -> 730,
    #   an 11% error, because a real cut carries aromatics and cycloalkanes a pure
    #   n-alkane does not. Sourcing a row must not change what the row names.
    #   C3.5 residual - searched for, and the search is recorded.
    "SAE 90 Gear Oil": 920,
    
    # Organic Liquids
    # [ON-DISK] pubchem/carbon_tetrachloride_density.json @ cid=5943 heading="Density" ref=6 scale=1000 precision=exact
    #   1.59 at 68 °F (USCG, 1999) - Denser than water; will sink
    #   Source: CAMEO Chemicals; conditions: 68 F.
    #   That entry cites: U.S. Coast Guard. 1999. Chemical Hazard Response Information System (CHRIS) - Hazardous Chemical Dat
    #   PubChem is an AGGREGATOR and this record carries several densities that
    #   disagree, so the citation names ONE ReferenceNumber rather than the CID alone.
    "Carbon Tetrachloride": 1590,
    # [ON-DISK] pubchem/chloroform_density.json @ cid=6212 heading="Density" ref=7 scale=1000 precision=3sf
    #   1.4832 at 68 °F (EPA, 1998) - Denser than water; will sink
    #   Source: CAMEO Chemicals; conditions: 68 F.
    #   That entry cites: U.S. Environmental Protection Agency. 1998. Extremely Hazardous Substances (EHS) Chemical Profiles a
    #   PubChem is an AGGREGATOR and this record carries several densities that
    #   disagree, so the citation names ONE ReferenceNumber rather than the CID alone.
    "Chloroform": 1480,
    # [ON-DISK] pubchem/bromine_density.json @ cid=24408 heading="Density" ref=55 scale=1000 precision=exact
    #   3.12
    #   Source: The National Institute for Occupational Safety and Health (NIOSH); conditions: no temperature stated in the entry.
    #   PubChem is an AGGREGATOR and this record carries several densities that
    #   disagree, so the citation names ONE ReferenceNumber rather than the CID alone.
    "Bromine": 3120,
    
    # Inorganic Solutions
    # [UNVERIFIED] no source on disk. Residual register (C3.5).
    "Calcium Chloride Solution (40%)": 1390,
    # [UNVERIFIED] no source on disk. Residual register (C3.5).
    "Zinc Chloride Solution (50%)": 1520,
    # [UNVERIFIED] no source on disk. Residual register (C3.5).
    "Sodium Polysulfide": 1650,
    
    # Specialty Fluids
    # [ON-DISK] pubchem/glycerol_density.json @ cid=753 heading="Density" ref=11 scale=1000 precision=3sf
    #   1.261 at 68 °F (USCG, 1999) - Denser than water; will sink
    #   Source: CAMEO Chemicals; conditions: 68 F.
    #   That entry cites: U.S. Coast Guard. 1999. Chemical Hazard Response Information System (CHRIS) - Hazardous Chemical Dat
    #   PubChem is an AGGREGATOR and this record carries several densities that
    #   disagree, so the citation names ONE ReferenceNumber rather than the CID alone.
    "Glycerin": 1260,
    # [UNVERIFIED] a source IS on disk, and it disagrees - the previous tag said "no
    #   source on disk", which was false once the PubChem records were acquired.
    #   pubchem/diiodomethane_density.json (cid=6346) gives 3.33 g/cm3
    #   = 3330 kg/m3 at 15 C, from PAC Chemical Database, U.S. Department of Energy,
    #   against this row's 3325. That entry names no primary, and it is the only
    #   unambiguous value in the record, so it is reported rather than cited: a single
    #   aggregator entry with no primary is not enough to move a value on.
    #   C3.5 residual - searched, and the search is recorded.
    "Diiodomethane": 3325,
    # [UNVERIFIED] the value is unsourced: no artefact on disk states or
    #   contradicts 2960.
    #   DUPLICATE RESOLVED (Phase 6, owner-directed). This table also carried
    #   "Tetrabromoethane" at the same 2960 - one substance under two names, which
    #   over-weighted it in every random.choice over this table. That row is now
    #   DELETED and this one kept.
    #   THIS name was kept because it is the unambiguous one: acetylene tetrabromide
    #   IS 1,1,2,2-tetrabromoethane, whereas a bare "tetrabromoethane" could in
    #   principle name the 1,1,1,2 isomer. No artefact was read to choose between
    #   them - the choice is about which NAME is less ambiguous, not about which
    #   VALUE is right. The value itself remains unsourced. C3.5 residual.
    "Acetylene Tetrabromide": 2960,
    
    # Molten Metals (for high-temperature applications)
    # [ON-DISK] pubchem/gallium_density.json @ cid=5360835 heading="Density" ref=23 scale=1000 precision=4sf
    #   6.0947 @ 29.8 °C (liquid); 5.9037 @ 29.65 °C (solid)
    #   Source: Hazardous Substances Data Bank (HSDB); conditions: 29.8 C.
    #   That entry cites: Budavari, S. (ed.). The Merck Index - An Encyclopedia of Chemicals, Drugs, and Biologicals. Whitehou
    #   PubChem is an AGGREGATOR and this record carries several densities that
    #   disagree, so the citation names ONE ReferenceNumber rather than the CID alone.
    "Gallium": 6095,
    # [UNVERIFIED] no source on disk. Residual register (C3.5).
    "Tin": 6980,
    # [UNVERIFIED] no source on disk. Residual register (C3.5).
    "Zinc": 6570,
    
}
