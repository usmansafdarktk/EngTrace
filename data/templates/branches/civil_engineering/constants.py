"""
Civil Engineering — authoritative constants for pilot template parameterization.

Stage B artifact (docs/pilot_template_authoring_spec.md). Every entry cites its
source inline. Sources marked [ON-DISK] were transcribed directly from the
documents in pilot/references/public/ (see MANIFEST.md); entries marked
[VERIFY:...] were curated from the named standard reference and MUST be
confirmed by the independent Data Reviewer before Stage C uses them.

Unit conventions:
  - SI values are primary; US customary provided where the dual-unit
    convention of the benchmark requires it.
  - AISC section properties follow the database's own units, noted per field.
"""

# ============================================================================
# UNIVERSAL CONSTANTS
# ============================================================================

GRAVITY_M_S2 = 9.81                  # standard gravity, m/s^2
GRAVITY_FT_S2 = 32.2                 # standard gravity, ft/s^2
UNIT_WEIGHT_WATER_KN_M3 = 9.81       # gamma_w at ~15-20 C, kN/m^3
UNIT_WEIGHT_WATER_PCF = 62.4         # gamma_w, lb/ft^3
WATER_DENSITY_KG_M3 = 998.2          # rho at 20 C  [VERIFY: CRC Handbook]
WATER_KINEMATIC_VISCOSITY_M2_S = 1.004e-6   # nu at 20 C  [VERIFY: CRC Handbook]

# ============================================================================
# DOMAIN 1 — STRUCTURAL ANALYSIS
# ============================================================================

# Steel properties — AISC Steel Construction Manual, 16th ed.
STEEL_E_KSI = 29000        # modulus of elasticity  [VERIFY: AISC Manual]
STEEL_E_GPA = 200          # SI companion value      [VERIFY: AISC Manual]
STEEL_FY_KSI = {           # yield stress by common grade  [VERIFY: AISC/ASTM]
    "ASTM A992": 50,       # wide-flange standard grade
    "ASTM A36": 36,        # plates, angles, legacy shapes
}
STEEL_FY_MPA = {"ASTM A992": 345, "ASTM A36": 250}

# Normal-weight concrete — ACI 318-19 Section 19.2.2.1
# Ec = 57000*sqrt(f'c) psi  |  Ec = 4700*sqrt(f'c) MPa   [VERIFY: ACI 318-19]
CONCRETE_FC_PSI = [3000, 4000, 5000]
CONCRETE_FC_MPA = [21, 28, 35]
CONCRETE_EC_COEFF_PSI = 57000
CONCRETE_EC_COEFF_MPA = 4700

# Typical uniform floor live loads — ASCE 7-22 Table 4.3-1  [VERIFY: ASCE 7]
LIVE_LOADS_PSF = {
    "office": 50,
    "residential dwelling": 40,
    "school classroom": 40,
    "corridor (first floor)": 100,
    "ordinary flat roof": 20,
}
LIVE_LOADS_KPA = {
    "office": 2.40,
    "residential dwelling": 1.92,
    "school classroom": 1.92,
    "corridor (first floor)": 4.79,
    "ordinary flat roof": 0.96,
}

# W-shape section properties — AISC Shapes Database v16.0  [ON-DISK: xlsx]
# 'us': W lb/ft, A in^2, d in, Ix in^4, Sx in^3, Zx in^3, rx in
# 'si': W kg/m,  A mm^2, d mm, Ix 10^6 mm^4, Sx 10^3 mm^3, Zx 10^3 mm^3, rx mm
AISC_W_SHAPES = {
    "W8X24":   {"us": {"W": 24,  "A": 7.08, "d": 7.93, "Ix": 82.7, "Sx": 20.9, "Zx": 23.1, "rx": 3.42},
                "si_label": "W200X35.9", "si": {"W": 35.9, "A": 4570, "d": 201, "Ix": 34.4, "Sx": 342, "Zx": 379, "rx": 86.9}},
    "W10X30":  {"us": {"W": 30,  "A": 8.84, "d": 10.5, "Ix": 170,  "Sx": 32.4, "Zx": 36.6, "rx": 4.38},
                "si_label": "W250X44.8", "si": {"W": 44.8, "A": 5700, "d": 267, "Ix": 70.8, "Sx": 531, "Zx": 600, "rx": 111}},
    "W12X26":  {"us": {"W": 26,  "A": 7.65, "d": 12.2, "Ix": 204,  "Sx": 33.4, "Zx": 37.2, "rx": 5.17},
                "si_label": "W310X38.7", "si": {"W": 38.7, "A": 4940, "d": 310, "Ix": 84.9, "Sx": 547, "Zx": 610, "rx": 131}},
    "W12X40":  {"us": {"W": 40,  "A": 11.7, "d": 11.9, "Ix": 307,  "Sx": 51.5, "Zx": 57,   "rx": 5.13},
                "si_label": "W310X60",   "si": {"W": 60,  "A": 7550, "d": 302, "Ix": 128,  "Sx": 844, "Zx": 934, "rx": 130}},
    "W14X30":  {"us": {"W": 30,  "A": 8.85, "d": 13.8, "Ix": 291,  "Sx": 42,   "Zx": 47.3, "rx": 5.73},
                "si_label": "W360X44",   "si": {"W": 44,  "A": 5710, "d": 351, "Ix": 121,  "Sx": 688, "Zx": 775, "rx": 146}},
    "W16X31":  {"us": {"W": 31,  "A": 9.13, "d": 15.9, "Ix": 375,  "Sx": 47.2, "Zx": 54,   "rx": 6.41},
                "si_label": "W410X46.1", "si": {"W": 46.1, "A": 5890, "d": 404, "Ix": 156, "Sx": 773, "Zx": 885, "rx": 163}},
    "W16X40":  {"us": {"W": 40,  "A": 11.8, "d": 16,   "Ix": 518,  "Sx": 64.7, "Zx": 73,   "rx": 6.63},
                "si_label": "W410X60",   "si": {"W": 60,  "A": 7610, "d": 406, "Ix": 216,  "Sx": 1060, "Zx": 1200, "rx": 168}},
    "W18X50":  {"us": {"W": 50,  "A": 14.7, "d": 18,   "Ix": 800,  "Sx": 88.9, "Zx": 101,  "rx": 7.38},
                "si_label": "W460X74",   "si": {"W": 74,  "A": 9480, "d": 457, "Ix": 333,  "Sx": 1460, "Zx": 1660, "rx": 187}},
    "W21X62":  {"us": {"W": 62,  "A": 18.3, "d": 21,   "Ix": 1330, "Sx": 127,  "Zx": 144,  "rx": 8.54},
                "si_label": "W530X92",   "si": {"W": 92,  "A": 11800, "d": 533, "Ix": 554, "Sx": 2080, "Zx": 2360, "rx": 217}},
    "W24X76":  {"us": {"W": 76,  "A": 22.4, "d": 23.9, "Ix": 2100, "Sx": 176,  "Zx": 200,  "rx": 9.69},
                "si_label": "W610X113",  "si": {"W": 113, "A": 14500, "d": 607, "Ix": 874, "Sx": 2880, "Zx": 3280, "rx": 246}},
    "W27X94":  {"us": {"W": 94,  "A": 27.6, "d": 26.9, "Ix": 3270, "Sx": 243,  "Zx": 278,  "rx": 10.9},
                "si_label": "W690X140",  "si": {"W": 140, "A": 17800, "d": 683, "Ix": 1360, "Sx": 3980, "Zx": 4560, "rx": 277}},
    "W30X108": {"us": {"W": 108, "A": 31.7, "d": 29.8, "Ix": 4470, "Sx": 299,  "Zx": 346,  "rx": 11.9},
                "si_label": "W760X161",  "si": {"W": 161, "A": 20500, "d": 757, "Ix": 1860, "Sx": 4900, "Zx": 5670, "rx": 302}},
    "W33X118": {"us": {"W": 118, "A": 34.7, "d": 32.9, "Ix": 5900, "Sx": 359,  "Zx": 415,  "rx": 13.0},
                "si_label": "W840X176",  "si": {"W": 176, "A": 22400, "d": 836, "Ix": 2460, "Sx": 5880, "Zx": 6800, "rx": 330}},
    "W36X135": {"us": {"W": 135, "A": 39.9, "d": 35.6, "Ix": 7800, "Sx": 439,  "Zx": 509,  "rx": 14.0},
                "si_label": "W920X201",  "si": {"W": 201, "A": 25700, "d": 904, "Ix": 3250, "Sx": 7190, "Zx": 8340, "rx": 356}},
}

# ============================================================================
# DOMAIN 2 — GEOTECHNICAL ENGINEERING
# ============================================================================

# Specific gravity of solids — Das PGE 10th ed., typical values  [VERIFY: Das]
SPECIFIC_GRAVITY_RANGES = {
    "sand": (2.65, 2.67),
    "silt": (2.67, 2.73),
    "inorganic clay": (2.70, 2.80),
}

# Coefficient of permeability by soil type, cm/s — Das PGE Table 7.1;
# cross-ref NAVFAC DM-7.01 Ch. 3 [ON-DISK]  [VERIFY: Das]
PERMEABILITY_RANGES_CM_S = {
    "clean gravel": (1.0, 100.0),
    "coarse sand": (0.01, 1.0),
    "fine sand": (0.001, 0.01),
    "silty clay": (1e-5, 0.001),
    "clay": (1e-8, 1e-6),
}

# Natural-state soil properties — Das & Sobhan, PGE 9th ed., Table 3.1
# ("Void Ratio, Moisture Content, and Dry Unit Weight for Some Typical Soils
# in a Natural State", text p. 68 / PDF p. 93)  [ON-DISK: full_books]
# Fields: void ratio e (-), saturated-state moisture content w (%), dry unit
# weight gamma_d (kN/m^3 and lb/ft^3). Tuples denote the table's own ranges.
# Derived quantities are computed, never assumed:
#   gamma_sat = gamma_d + (e/(1+e)) * gamma_w
#   gamma_moist = gamma_d * (1 + w/100) for any assumed w <= w_sat
DAS_NATURAL_STATE_SOILS = {
    "loose uniform sand":          {"e": 0.8,        "w_sat_pct": 30,        "gamma_d_kn_m3": 14.5,         "gamma_d_pcf": 92},
    "dense uniform sand":          {"e": 0.45,       "w_sat_pct": 16,        "gamma_d_kn_m3": 18.0,         "gamma_d_pcf": 115},
    "loose angular silty sand":    {"e": 0.65,       "w_sat_pct": 25,        "gamma_d_kn_m3": 16.0,         "gamma_d_pcf": 102},
    "dense angular silty sand":    {"e": 0.4,        "w_sat_pct": 15,        "gamma_d_kn_m3": 19.0,         "gamma_d_pcf": 121},
    "stiff clay":                  {"e": 0.6,        "w_sat_pct": 21,        "gamma_d_kn_m3": 17.0,         "gamma_d_pcf": 108},
    "soft clay":                   {"e": (0.9, 1.4), "w_sat_pct": (30, 50),  "gamma_d_kn_m3": (11.5, 14.5), "gamma_d_pcf": (73, 93)},
    "loess":                       {"e": 0.9,        "w_sat_pct": 25,        "gamma_d_kn_m3": 13.5,         "gamma_d_pcf": 86},
    "soft organic clay":           {"e": (2.5, 3.2), "w_sat_pct": (90, 120), "gamma_d_kn_m3": (6.0, 8.0),   "gamma_d_pcf": (38, 51)},
    "glacial till":                {"e": 0.3,        "w_sat_pct": 10,        "gamma_d_kn_m3": 21.0,         "gamma_d_pcf": 134},
}

# Drained friction angle ranges, degrees — Das PGE Ch. 12 typical values
# [VERIFY: Das]
FRICTION_ANGLE_RANGES_DEG = {
    "sand, rounded, loose": (27, 30),
    "sand, rounded, dense": (35, 38),
    "sand, angular, loose": (30, 35),
    "sand, angular, dense": (40, 45),
    "gravel with some sand": (34, 48),
    "silt": (26, 35),
}

# Terzaghi bearing-capacity factors (Nc, Nq, Ngamma) vs. friction angle —
# transcribed verbatim from Das & Sobhan PGE 9th ed., Table 16.1 (text
# p. 717 / PDF p. 742); Ngamma per Kumbhojkar (1993).  [ON-DISK: full_books]
# NOTE: an earlier memory-curated version of this table carried wrong
# Ngamma values at phi <= 25 and phi = 40 (mixed factor families) that
# web-based Stage B review failed to catch; corrected 2026-08-02 against
# the book text after an R2 finding (see data_review_log.md, Cycle 3).
TERZAGHI_BEARING_FACTORS = {
    0:  (5.70, 1.00, 0.00),
    5:  (7.34, 1.64, 0.14),
    10: (9.61, 2.69, 0.56),
    15: (12.86, 4.45, 1.52),
    20: (17.69, 7.44, 3.64),
    25: (25.13, 12.72, 8.34),
    30: (37.16, 22.46, 19.13),
    35: (57.75, 41.44, 45.41),
    40: (95.66, 81.27, 116.31),
}

# Terzaghi's MODIFIED bearing-capacity factors (N'c, N'q, N'gamma) for the
# LOCAL shear failure mode (phi' replaced by atan(2/3 tan phi')) —
# transcribed verbatim from Das & Sivakugan, PFE 9th ed., Table 3.2 (text
# p. 140 / PDF p. 160); local-shear strip equation per Eq. (3.9):
# qu = (2/3)c'N'c + qN'q + 0.5*gamma*B*N'gamma.  [ON-DISK: full_books]
TERZAGHI_MODIFIED_FACTORS = {
    20: (11.85, 3.88, 1.12),
    25: (14.80, 5.60, 2.25),
    30: (18.99, 8.31, 4.39),
    35: (25.18, 12.75, 8.35),
}

# Compression index correlation: Cc = 0.009 * (LL - 10)  — Skempton (1944),
# as given in Das PGE consolidation chapter  [VERIFY: Das]
SKEMPTON_CC_COEFF = 0.009
SKEMPTON_CC_OFFSET = 10

# Typical coefficient of consolidation, m^2/yr — order-of-magnitude ranges;
# reviewer evidence shows published low-plasticity values as low as ~2 m^2/yr,
# so the ranges deliberately overlap.  [POLICY: sampling-only. Das PGE has no
# general cv table and the DM-7.01 cv-LL chart is image-only; any sampled cv
# MUST appear verbatim as a given value in the question text, so template
# correctness never depends on these endpoints.]
CV_RANGES_M2_YR = {
    "high-plasticity clay": (0.3, 5.0),
    "low-plasticity clay": (2.0, 30.0),
}

# ============================================================================
# DOMAIN 3 — WATER RESOURCES & HYDRAULICS
# ============================================================================

# Manning's n, open channels — FHWA HDS-4 Table B.2, p. B-2 (PDF p.199)
# [ON-DISK]; canonical compilation: Chow (1959) Table 5-6, USGS WSP-2339
MANNINGS_N_CHANNELS = {
    "very smooth concrete": 0.011,
    "smooth concrete": 0.012,
    "ordinary concrete lining": 0.013,
    "wood": 0.014,
    "vitrified clay": 0.015,
    "shot concrete / earth channel, best condition": 0.017,
    "straight unlined earth canal, good condition": 0.020,
    "mountain stream, rocky bed": (0.040, 0.050),
    "minor stream, clean and straight": (0.025, 0.033),
    "minor stream, winding, some pools": (0.033, 0.045),
    "minor stream, sluggish and weedy": (0.050, 0.080),
    "floodplain, short grass pasture": (0.025, 0.035),
}

# Manning's n, closed conduits — FHWA HDS-4 Table B.3, p. B-4 (PDF p.201)
# [ON-DISK]
MANNINGS_N_CONDUITS = {
    "concrete pipe": (0.011, 0.013),
    "CMP, 2-2/3 x 1/2 in corrugations": (0.022, 0.027),
    "vitrified clay pipe": (0.012, 0.014),
    "steel pipe": (0.009, 0.013),
    "brick": (0.014, 0.017),
}

# Rational-method runoff coefficients — FHWA HEC-22 3rd ed., Table 3-1,
# p. 3-6 (PDF p.56)  [ON-DISK]
RATIONAL_C = {
    "business: downtown": (0.70, 0.95),
    "business: neighborhood": (0.50, 0.70),
    "residential: single-family": (0.30, 0.50),
    "residential: multi-unit detached": (0.40, 0.60),
    "residential: multi-unit attached": (0.60, 0.75),
    "residential: suburban": (0.25, 0.40),
    "residential: apartments": (0.50, 0.70),
    "industrial: light": (0.50, 0.80),
    "industrial: heavy": (0.60, 0.90),
    "parks and cemeteries": (0.10, 0.25),
    "playgrounds": (0.20, 0.40),
    "lawns: sandy soil, flat (2%)": (0.05, 0.10),
    "lawns: heavy soil, steep (7%)": (0.25, 0.35),
    "streets: asphaltic": (0.70, 0.95),
    "streets: concrete": (0.80, 0.95),
    "roofs": (0.75, 0.95),
}

# SCS runoff curve numbers by hydrologic soil group (A, B, C, D) —
# NRCS TR-55 (June 1986), Table 2-2a (urban, doc p.2-5) and Table 2-2b/c
# (agricultural, doc p.2-6/2-7)  [ON-DISK]
SCS_CURVE_NUMBERS = {
    # Table 2-2a — urban
    "open space, poor condition (<50% grass)": (68, 79, 86, 89),
    "open space, fair condition (50-75% grass)": (49, 69, 79, 84),
    "open space, good condition (>75% grass)": (39, 61, 74, 80),
    "paved parking, roofs, driveways": (98, 98, 98, 98),
    "streets: paved with curbs and storm sewers": (98, 98, 98, 98),
    "streets: paved with open ditches": (83, 89, 92, 93),
    "streets: gravel": (76, 85, 89, 91),
    "streets: dirt": (72, 82, 87, 89),
    "commercial and business (85% impervious)": (89, 92, 94, 95),
    "industrial (72% impervious)": (81, 88, 91, 93),
    "residential: 1/8 acre lots (65% impervious)": (77, 85, 90, 92),
    "residential: 1/4 acre lots (38% impervious)": (61, 75, 83, 87),
    "residential: 1/3 acre lots (30% impervious)": (57, 72, 81, 86),
    "residential: 1/2 acre lots (25% impervious)": (54, 70, 80, 85),
    "residential: 1 acre lots (20% impervious)": (51, 68, 79, 84),
    "residential: 2 acre lots (12% impervious)": (46, 65, 77, 82),
    "newly graded area (no vegetation)": (77, 86, 91, 94),
    # Table 2-2b — cultivated agricultural
    "fallow, bare soil": (77, 86, 91, 94),
    "row crops, straight row, poor condition": (72, 81, 88, 91),
    "row crops, straight row, good condition": (67, 78, 85, 89),
    "small grain, straight row, good condition": (63, 75, 83, 87),
    # Table 2-2c — other agricultural  [VERIFY: TR-55 p.2-7]
    "pasture, good condition": (39, 61, 74, 80),
    "meadow, protected from grazing": (30, 58, 71, 78),
    "woods, good condition": (30, 55, 70, 77),
    "farmsteads": (59, 74, 82, 86),
}

# SCS runoff equation constants — TR-55 Ch. 2: Q = (P - 0.2S)^2 / (P + 0.8S),
# S = 1000/CN - 10 (inches)  [ON-DISK]
SCS_IA_RATIO = 0.2

HYDROLOGIC_SOIL_GROUPS = ("A", "B", "C", "D")
