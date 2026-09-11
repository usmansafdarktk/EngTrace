"""
Industrial Engineering / Operations Research — authoritative constants for
pilot template parameterization.

Stage B artifact (docs/pilot_template_authoring_spec.md). Every entry cites its
source inline. Tags:
  [ON-DISK]  — transcribed verbatim from documents in pilot/references/public/
               (see MANIFEST.md), with the page noted. For image-only sources
               the transcription was done from rendered page images.
  [DERIVABLE] — mathematical constants; must be verified by independent
               derivation (not table lookup) by the Data Reviewer.
  [REALISM]  — plausibility screens for sampling, anchored to the named book's
               worked-example conventions. NOT empirical data. The GIVEN-VALUES
               RULE applies: every sampled value is stated verbatim in the
               question text, so these ranges affect scenario realism only,
               never numerical correctness (policy carried over from the Civil
               branch data review, Cycle 2).

Unit conventions:
  - This branch is dimensionally light (BOOKS.md §6): rates per hour/year,
    counts, probabilities, currency (USD). No SI/US dual-unit duality.
  - All rates are stated with their time unit; queueing templates that mix
    time units must state both explicitly in the question.
  - Normal distribution: Phi and Phi^-1 are computed with
    statistics.NormalDist() (stdlib, deterministic); solutions display z and
    Phi(z) at the declared precision so a student can follow with a z-table.
"""

# ============================================================================
# UNIVERSAL / MATHEMATICAL CONSTANTS
# ============================================================================

# Standard-normal quantiles z_p for common service levels / confidence levels.
# Cross-reference: Montgomery ISQC 7e Appendix Table II (cumulative normal);
# H&L 7e Appendix 5 statistical tables.
# @kind: mathematical
# @units: 1
# @domain: none (quantiles of the standard normal; no physical condition)
# [DERIVED] z_p = statistics.NormalDist().inv_cdf(p), rounded to 4 dp. Re-derived, not
#   looked up: all 6 entries are reproduced exactly by the stdlib at 4 dp (the patch
#   script asserts it, and tests/constants_integrity does the same on every run). The
#   Montgomery and H&L tables named below are a cross-reference, not the source.
Z_QUANTILES = {
    0.90:  1.2816,
    0.95:  1.6449,
    0.975: 1.9600,
    0.98:  2.0537,
    0.99:  2.3263,
    0.995: 2.5758,
}

# Three-sigma control-limit convention (Shewhart charts): +/- 3 standard
# deviations of the plotted statistic — Montgomery ISQC 7e §5.3.2.
# @kind: defined
# @units: 1
# @domain: none (the Shewhart limit convention; no condition applies)
# [ON-DISK] industrial/nist_sematech_pmc32.htm @ text="we speak of 3-sigma control charts"
#   The NIST/SEMATECH e-Handbook §6.3.2 states the convention in words; k is a choice,
#   not a measurement, so the locator pins the statement rather than a number.
SHEWHART_K_SIGMA = 3

# ============================================================================
# DOMAIN 1 — STOCHASTIC OPERATIONS
# ============================================================================

# Named queueing scenarios: arrival/service-rate windows chosen so templates
# can enforce the steady-state gate rho = lambda/(c*mu) < 1 by construction
# (BOOKS.md §4). Anchored to the worked-example conventions of
# H&L 7e Ch. 17 (county hospital ER, ~1-3 customers/hr) and Taha 10e Ch. 18
# (bank / tool crib / repair examples). Rates are per HOUR.
# @kind: range
# @units: lam_hr=1/h, mu_hr=1/h, servers=1
# @given: stated (C1 Reviewer G, reviews/phaseC1_reviewer_g_provenance.md §2 row 10 and reviews/phaseC1_reviewer_g_scratch/g_stated.py: rates stated 40/40 for M/M/1, lambda and c stated 40/40 for M/M/c)
# @domain: none (arrival/service-rate windows for scenario generation)
# [POLICY: sampling-only] Sampled values are stated verbatim in the question (the given-values rule), so this window changes which scenario is posed, never whether the answer is right. The windows are chosen so rho < 1 holds by
#   construction; H&L 7e Ch. 17 and Taha 10e Ch. 18 ground the TYPOLOGY, not the numbers.
QUEUE_SCENARIOS = {
    "bank teller line":        {"lam_hr": (8, 40),  "mu_hr": (12, 50),  "servers": (1, 4)},
    "call center":             {"lam_hr": (20, 90), "mu_hr": (10, 30),  "servers": (2, 6)},
    "hospital emergency room": {"lam_hr": (1, 6),   "mu_hr": (2, 8),    "servers": (1, 3)},
    "tool crib counter":       {"lam_hr": (4, 18),  "mu_hr": (6, 24),   "servers": (1, 3)},
    "drive-through window":    {"lam_hr": (15, 55), "mu_hr": (20, 70),  "servers": (1, 2)},
    "machine repair shop":     {"lam_hr": (0.5, 4), "mu_hr": (1, 6),    "servers": (1, 3)},
}

# Waiting-cost / service-cost windows (USD per hour) for queueing economic-
# comparison templates (S1-#3): per the cost-analysis conventions of
# H&L 7e Ch. 18 "The Application of Queueing Theory" (waiting cost >> server
# cost ratios in its examples) and Taha 10e §18.9 "Queuing Decision Models"
# (§18.9.1 Cost Models). Sampled values are always stated in the question.
# @kind: range
# @units: USD/h
# @domain: none (cost windows for scenario generation)
# [POLICY: sampling-only] Sampled values are stated verbatim in the question (the given-values rule), so this window changes which scenario is posed, never whether the answer is right.
QUEUE_COSTS_USD_HR = {
    "waiting_cost_per_customer": (10, 120),
    "server_cost_per_server":    (8, 60),
}

# Finite-capacity (M/M/1/K) system sizes: waiting-room capacities K
# (including the one in service) used by S1-#4. Small-buffer
# service systems per H&L 7e §17.6 finite-queue variation.
# @kind: range
# @units: 1
# @domain: none (a buffer-size window for scenario generation)
# [POLICY: sampling-only] Sampled values are stated verbatim in the question (the given-values rule), so this window changes which scenario is posed, never whether the answer is right.
FINITE_CAPACITY_K = (3, 10)

# Component reliability classes: per-component reliability windows for
# mission-time system-reliability templates (S3-#8).
# [POLICY: sampling-only] Ross 11e Ch. 9 grounds the TYPOLOGY only
# (its examples work symbolically with p_i — no numeric reliability bands to
# anchor to; data review 2026-08-05). Windows reflect common engineering
# grading practice; the given-values rule makes them correctness-neutral.
# Logged for the human-expert escalation list.
# @kind: range
# @units: 1
# @given: stated (C1 Reviewer G, reviews/phaseC1_reviewer_g_provenance.md §2 row 12 and reviews/phaseC1_reviewer_g_scratch/g_stated.py: every component reliability stated, 40/40)
# @domain: none (per-component reliability windows for scenario generation)
COMPONENT_RELIABILITY_CLASSES = {
    "commercial grade":  (0.90, 0.97),
    "industrial grade":  (0.95, 0.99),
    "high reliability":  (0.99, 0.999),
}

# Exponential failure-rate windows (failures per hour) by component family,
# used for MTTF templates (S3-#9).
# [POLICY: sampling-only] The NIST/SEMATECH e-Handbook Ch. 8 (apr/)
# defines failure/hazard rates (§8.1.2) but tabulates no magnitude bands per
# component family (data review 2026-08-05 — an earlier citation of
# §8.1.2/8.1.10 as a magnitude anchor was wrong and is withdrawn). Windows
# are order-of-magnitude engineering practice; given-values rule applies.
# Logged for the human-expert escalation list.
# @kind: range
# @units: 1/h
# @domain: none (order-of-magnitude failure-rate windows by component family)
FAILURE_RATE_PER_HR = {
    "electronic module":  (1e-6, 5e-5),
    "power supply":       (5e-6, 1e-4),
    "pump / mechanical":  (5e-5, 1e-3),
}

# Mission times for reliability evaluation, hours.
# @kind: range
# @units: h
# @domain: none (a mission-time window for scenario generation)
# [POLICY: sampling-only] Sampled values are stated verbatim in the question (the given-values rule), so this window changes which scenario is posed, never whether the answer is right.
MISSION_TIME_HR = (100, 5000)

# Markov-chain scenario families for S2 (weather, machine up/down, brand
# switching, credit-rating migration) carry no numeric constants — transition
# probabilities are sampled and fully stated in the question (given-values
# rule). Typology: Ross 11e Ch. 4 examples 4.1/4.3; H&L 7e Ch. 16.

# ============================================================================
# DOMAIN 2 — PRODUCTION & INVENTORY
# ============================================================================

# Annual holding-cost rate i (fraction of unit value per year): interest +
# storage + obsolescence components. Nahmias 7e §4.4 "Holding Cost"
# (text pp. 204-205, h = Ic; image-only copy — anchored by visual read). The
# section's own illustration builds I = 0.37 (28% capital + 2% taxes/ins. +
# 6% storage + 1% breakage); Ch. 4 problems use values near 0.22 — window
# spans both. Sampled i is always stated.
# @kind: range
# @units: 1/yr
# @given: stated (tests/constants_integrity/given_evidence.py: i is printed at the 2 dp it is drawn at, 40/40, all three consumers)
# @domain: none (an annual holding-rate window for scenario generation)
# [POLICY: sampling-only] Sampled values are stated verbatim in the question (the given-values rule), so this window changes which scenario is posed, never whether the answer is right.
HOLDING_RATE_PER_YR = (0.15, 0.40)

# Named inventory item classes for EOQ-family templates (P1): unit cost c,
# annual demand D, setup/order cost K, lead time. Magnitude conventions per
# Nahmias 7e Ch. 4 worked examples and H&L 7e Ch. 19 prototype examples
# (speaker/bicycle-class items). Never verbatim numbers.
# @kind: range
# @units: unit_cost_usd=USD, annual_demand=1/yr, order_cost_usd=USD
# @given: stated (tests/constants_integrity/given_evidence.py: c, D and K - or the price schedule - are printed, 40/40, all three consumers)
# @domain: none (item-class windows for scenario generation)
# [POLICY: sampling-only] Sampled values are stated verbatim in the question (the given-values rule), so this window changes which scenario is posed, never whether the answer is right.
INVENTORY_ITEMS = {
    "electronic component": {"unit_cost_usd": (2, 40),    "annual_demand": (2000, 60000), "order_cost_usd": (40, 300)},
    "machine spare part":   {"unit_cost_usd": (20, 400),  "annual_demand": (100, 5000),   "order_cost_usd": (60, 500)},
    "packaged food case":   {"unit_cost_usd": (8, 60),    "annual_demand": (5000, 80000), "order_cost_usd": (30, 200)},
    "industrial fastener":  {"unit_cost_usd": (0.5, 8),   "annual_demand": (10000, 200000), "order_cost_usd": (25, 150)},
    "retail appliance":     {"unit_cost_usd": (80, 900),  "annual_demand": (300, 8000),   "order_cost_usd": (50, 400)},
}

# Production rates for EPQ (finite production rate) templates: expressed as a
# multiple of the demand rate, P = m * D with m in the window below (m > 1
# guarantees feasibility). Nahmias 7e §4.6 convention.
# @kind: range
# @units: 1
# @domain: none (a production-to-demand ratio window; m > 1 keeps EPQ feasible)
# [POLICY: sampling-only] Sampled values are stated verbatim in the question (the given-values rule), so this window changes which scenario is posed, never whether the answer is right.
EPQ_PRODUCTION_MULTIPLE = (1.5, 6.0)

# All-units quantity-discount structures for P1-#3: 2-3 price breaks with
# successive discounts of 2-12% per break. Structure per Nahmias 7e
# §4.7 (all-units schedules) and H&L 7e Ch. 19 discount discussion. The full
# schedule is always printed in the question.
# @kind: range
# @units: 1
# @domain: none (a breakpoint-quantity window for scenario generation)
# [POLICY: sampling-only] Sampled values are stated verbatim in the question (the given-values rule), so this window changes which scenario is posed, never whether the answer is right.
DISCOUNT_BREAK_QTY = (100, 5000)          # candidate breakpoint window
# @kind: range
# @units: 1
# @domain: none (a per-break price-reduction window)
# [POLICY: sampling-only] Sampled values are stated verbatim in the question (the given-values rule), so this window changes which scenario is posed, never whether the answer is right.
DISCOUNT_STEP_FRACTION = (0.02, 0.12)     # per-break price reduction

# Lead times for reorder-point templates (P1-#4), weeks; the tau > T branch
# (lead time exceeding a cycle) is exercised deliberately per Nahmias 7e
# §4.5 "Inclusion of Order Lead Time" (text p. 213).
# @kind: range
# @units: week
# @domain: none (a lead-time window for scenario generation)
# [POLICY: sampling-only] Sampled values are stated verbatim in the question (the given-values rule), so this window changes which scenario is posed, never whether the answer is right. The earlier [ON-DISK visual] claimed a page read
#   off an image; no locator can check that, and the window is policy either way.
LEAD_TIME_WEEKS = (1, 10)

# Newsvendor item families for P2-#6: unit cost c, selling price p > c,
# salvage s < c. Typology per Nahmias 7e §5.3 (perishables) and
# H&L 7e Ch. 19 stochastic single-period model.
# @kind: range
# @units: cost_usd=USD, price_usd=USD, salvage_usd=USD
# @given: stated (C1 Reviewer G, reviews/phaseC1_reviewer_g_provenance.md §2 row 11 and reviews/phaseC1_reviewer_g_scratch/g_stated.py: c, p and s stated, 40/40)
# @domain: none (newsvendor item-class windows for scenario generation)
# [POLICY: sampling-only] Sampled values are stated verbatim in the question (the given-values rule), so this window changes which scenario is posed, never whether the answer is right.
NEWSVENDOR_ITEMS = {
    "daily newspaper stack": {"cost_usd": (0.2, 1.0),  "price_usd": (0.5, 2.5),  "salvage_usd": (0.0, 0.3)},
    "bakery batch":          {"cost_usd": (1.0, 6.0),  "price_usd": (3.0, 15.0), "salvage_usd": (0.2, 2.0)},
    "seasonal apparel lot":  {"cost_usd": (8, 40),     "price_usd": (20, 120),   "salvage_usd": (2, 15)},
    "fresh produce crate":   {"cost_usd": (4, 20),     "price_usd": (10, 45),    "salvage_usd": (0, 5)},
}

# Cycle-service levels offered to safety-stock / (Q,R) templates; z comes
# from Z_QUANTILES (Type 1 service per Nahmias 7e §5.5).
# @kind: range
# @units: 1
# @given: stated (C1 Reviewer G, phaseC1_reviewer_g_provenance.md §2 row 4: the drawn level is stated in the question 40/40; the probe raises only because the level keys Z_QUANTILES)
# @domain: none (the service levels a question may be posed at)
# [POLICY: sampling-only] Sampled values are stated verbatim in the question (the given-values rule), so this window changes which scenario is posed, never whether the answer is right. Each level keys Z_QUANTILES, which is [DERIVED];
#   the level itself is a menu, not a measurement.
SERVICE_LEVELS = [0.90, 0.95, 0.98, 0.99]

# Aggregate-planning cost windows (P3-#10), USD: hiring/firing per worker,
# inventory holding per unit-month; monthly demand windows. Magnitude
# conventions per Nahmias 7e §3.4-3.5 worked examples.
# @kind: range
# @units: USD
# @domain: none (aggregate-planning cost windows for scenario generation)
# [POLICY: sampling-only] Sampled values are stated verbatim in the question (the given-values rule), so this window changes which scenario is posed, never whether the answer is right.
AGGREGATE_COSTS_USD = {
    "hire_per_worker":       (300, 1500),
    "fire_per_worker":       (500, 2500),
    "hold_per_unit_month":   (2, 20),
}
# @kind: range
# @units: 1/month
# @domain: none (a monthly-demand window per product family)
# [POLICY: sampling-only] Sampled values are stated verbatim in the question (the given-values rule), so this window changes which scenario is posed, never whether the answer is right.
AGGREGATE_MONTHLY_DEMAND = (200, 5000)     # units/month, per product family
# @kind: range
# @units: 1/month
# @domain: none (a per-worker monthly output window)
# [POLICY: sampling-only] Sampled values are stated verbatim in the question (the given-values rule), so this window changes which scenario is posed, never whether the answer is right.
WORKER_MONTHLY_OUTPUT = (20, 200)          # units/worker-month

# Assembly-line balancing (P3-#8/#9): task-time and cycle-time windows,
# seconds; task counts. Per Nahmias 7e §9.10 (text p. 528) — its examples use
# second/minute-scale station tasks.
# @kind: range
# @units: s
# @domain: none (a task-time window for line-balancing scenarios)
# [POLICY: sampling-only] Sampled values are stated verbatim in the question (the given-values rule), so this window changes which scenario is posed, never whether the answer is right.
LINE_TASK_TIME_S = (10, 120)
# @kind: range
# @units: 1
# @domain: none (a task-count window for line-balancing scenarios)
# [POLICY: sampling-only] Sampled values are stated verbatim in the question (the given-values rule), so this window changes which scenario is posed, never whether the answer is right.
LINE_TASK_COUNT = (5, 9)
# @kind: range
# @units: 1/shift
# @domain: none (a per-shift demand window)
# [POLICY: sampling-only] Sampled values are stated verbatim in the question (the given-values rule), so this window changes which scenario is posed, never whether the answer is right.
LINE_DEMAND_PER_SHIFT = (100, 800)         # units per 8-hour shift

# ============================================================================
# DOMAIN 3 — QUALITY & RELIABILITY CONTROL
# ============================================================================

# ----------------------------------------------------------------------------
# Factors for Constructing Variables Control Charts
# [DERIVED] every column from its definition by tests/constants_integrity/control_chart_factors.py:
#   c4 in closed form and the A2/A3/B3/B4 limit formulas from docs/references/industrial/nist_sematech_pmc32.htm and pmc321.htm; d2, d3 as the mean and sd of the relative range.
#   C3.7 corrected 29 last digits the transcription below carried (the table it was
#   copied from computes from rounded intermediates); the module compares all 384 cells.
# Cross-check, not the source: Montgomery, Introduction to Statistical Quality Control, 7th ed.,
# Appendix Table VI (text p. 720; PDF p. 738 of the supplied copy, clean text
# layer). Transcribed VERBATIM 2026-08-05. Columns:
#   A    = 3/sqrt(n)          (X-bar chart, sigma known)
#   A2   = 3/(d2*sqrt(n))     (X-bar chart from R-bar)
#   A3   = 3/(c4*sqrt(n))     (X-bar chart from s-bar)
#   c4, 1/c4                  (s-chart center-line factors)
#   B3, B4                    (s-chart limits from s-bar)
#   B5, B6                    (s-chart limits, sigma known)
#   d2, 1/d2, d3              (relative-range mean / sd factors)
#   D1, D2                    (R-chart limits, sigma known)
#   D3, D4                    (R-chart limits from R-bar)
# The Data Reviewer MUST verify these by independent derivation (c4, d2, d3
# from their Gamma-function / integral definitions; remaining columns from
# the identities above) — BOOKS.md §5 provenance-check requirement.
# ----------------------------------------------------------------------------
# @kind: mathematical
# @units: 1
# @domain: n=2..25
#   The factors are functions of the subgroup size alone; the table carries n = 2..25,
#   and a consumer outside that range has no row to read.
CONTROL_CHART_FACTORS = {
    #  n:  (A,     A2,    A3,    c4,     1/c4,   B3,    B4,    B5,    B6,    d2,    1/d2,   d3,    D1,    D2,    D3,    D4)
    2:  (2.121, 1.880, 2.659, 0.7979, 1.2533, 0.0,   3.267, 0.0,   2.606, 1.128, 0.8862, 0.853, 0.0,   3.686, 0.0,   3.267),
    3:  (1.732, 1.023, 1.954, 0.8862, 1.1284, 0.0,   2.568, 0.0,   2.276, 1.693, 0.5908, 0.888, 0.0,   4.358, 0.0,   2.575),
    4:  (1.500, 0.729, 1.628, 0.9213, 1.0854, 0.0,   2.266, 0.0,   2.088, 2.059, 0.4857, 0.880, 0.0,   4.698, 0.0,   2.282),
    5:  (1.342, 0.577, 1.427, 0.9400, 1.0638, 0.0,   2.089, 0.0,   1.964, 2.326, 0.4299, 0.864, 0.0,   4.918, 0.0,   2.114),
    6:  (1.225, 0.483, 1.287, 0.9515, 1.0509, 0.030, 1.970, 0.029, 1.874, 2.534, 0.3946, 0.848, 0.0,   5.079, 0.0,   2.004),
    7:  (1.134, 0.419, 1.182, 0.9594, 1.0424, 0.118, 1.882, 0.113, 1.806, 2.704, 0.3698, 0.833, 0.205, 5.204, 0.076, 1.924),
    8:  (1.061, 0.373, 1.099, 0.9650, 1.0362, 0.185, 1.815, 0.179, 1.751, 2.847, 0.3512, 0.820, 0.388, 5.307, 0.136, 1.864),
    9:  (1.000, 0.337, 1.032, 0.9693, 1.0317, 0.239, 1.761, 0.232, 1.707, 2.970, 0.3367, 0.808, 0.547, 5.394, 0.184, 1.816),
    10: (0.949, 0.308, 0.975, 0.9727, 1.0281, 0.284, 1.716, 0.276, 1.669, 3.078, 0.3249, 0.797, 0.686, 5.469, 0.223, 1.777),
    11: (0.905, 0.285, 0.927, 0.9754, 1.0253, 0.321, 1.679, 0.313, 1.637, 3.173, 0.3152, 0.787, 0.811, 5.535, 0.256, 1.744),
    12: (0.866, 0.266, 0.886, 0.9776, 1.0230, 0.354, 1.646, 0.346, 1.610, 3.258, 0.3069, 0.778, 0.923, 5.594, 0.283, 1.717),
    13: (0.832, 0.249, 0.850, 0.9794, 1.0210, 0.382, 1.618, 0.374, 1.585, 3.336, 0.2998, 0.770, 1.025, 5.647, 0.307, 1.693),
    14: (0.802, 0.235, 0.817, 0.9810, 1.0194, 0.406, 1.594, 0.399, 1.563, 3.407, 0.2935, 0.763, 1.118, 5.696, 0.328, 1.672),
    15: (0.775, 0.223, 0.789, 0.9823, 1.0180, 0.428, 1.572, 0.421, 1.544, 3.472, 0.2880, 0.756, 1.203, 5.740, 0.347, 1.653),
    16: (0.750, 0.212, 0.763, 0.9835, 1.0168, 0.448, 1.552, 0.440, 1.526, 3.532, 0.2831, 0.750, 1.282, 5.782, 0.363, 1.637),
    17: (0.728, 0.203, 0.739, 0.9845, 1.0157, 0.466, 1.534, 0.458, 1.511, 3.588, 0.2787, 0.744, 1.356, 5.820, 0.378, 1.622),
    18: (0.707, 0.194, 0.718, 0.9854, 1.0148, 0.482, 1.518, 0.475, 1.496, 3.640, 0.2747, 0.739, 1.424, 5.856, 0.391, 1.609),
    19: (0.688, 0.187, 0.698, 0.9862, 1.0140, 0.497, 1.503, 0.490, 1.483, 3.689, 0.2711, 0.733, 1.489, 5.889, 0.404, 1.596),
    20: (0.671, 0.180, 0.680, 0.9869, 1.0132, 0.510, 1.490, 0.504, 1.470, 3.735, 0.2677, 0.729, 1.549, 5.921, 0.415, 1.585),
    21: (0.655, 0.173, 0.663, 0.9876, 1.0126, 0.523, 1.477, 0.516, 1.459, 3.778, 0.2647, 0.724, 1.606, 5.951, 0.425, 1.575),
    22: (0.640, 0.167, 0.647, 0.9882, 1.0120, 0.534, 1.466, 0.528, 1.448, 3.819, 0.2618, 0.720, 1.660, 5.979, 0.435, 1.565),
    23: (0.626, 0.162, 0.633, 0.9887, 1.0114, 0.545, 1.455, 0.539, 1.438, 3.858, 0.2592, 0.716, 1.711, 6.006, 0.443, 1.557),
    24: (0.612, 0.157, 0.619, 0.9892, 1.0109, 0.555, 1.445, 0.549, 1.429, 3.895, 0.2567, 0.712, 1.759, 6.032, 0.452, 1.548),
    25: (0.600, 0.153, 0.606, 0.9896, 1.0105, 0.565, 1.435, 0.559, 1.420, 3.931, 0.2544, 0.708, 1.805, 6.056, 0.459, 1.541),
}
CONTROL_CHART_COLUMNS = ("A", "A2", "A3", "c4", "inv_c4", "B3", "B4", "B5",
                         "B6", "d2", "inv_d2", "d3", "D1", "D2", "D3", "D4")


def chart_factor(n, name):
    """Convenience accessor: chart_factor(5, 'A2') -> 0.577."""
    return CONTROL_CHART_FACTORS[n][CONTROL_CHART_COLUMNS.index(name)]


# ----------------------------------------------------------------------------
# MIL-STD-105E acceptance sampling (Q4)
# MIL-STD-105E (10 May 1989) — transcribed from rendered page images 2026-08-05.
# The PDF is a PUBLIC on-disk artefact (docs/references/industrial/mil_std_105e_sampling.pdf), so the three
# tables below cite it directly. Its OCR text layer is unreadable (p.18 extracts as
# ":: .... m n mO ::,i:,o Vtm Loi o, batcb ah:e ."), so each locator pins the PAGE
# and no text= anchor: a text anchor that cannot be trusted is worse than none.
# ----------------------------------------------------------------------------

# Table I — Sample size code letters, GENERAL INSPECTION LEVEL II (the default
# level per §4.9.1 "Inspection Level": "Normally, Inspection Level II is
# used."; Table I itself is captioned "see 4.9.1 and 4.9.2").
# (lot_min, lot_max, code_letter); None = unbounded.
# @kind: standard
# @units: 1
# @domain: inspection_level=II
# [ON-DISK] industrial/mil_std_105e_sampling.pdf @ page=18
#   Table I, sample size code letters, general inspection level II (PDF p.18 = doc p.13).
MIL_STD_105E_CODE_LETTERS_GII = [
    (2,      8,      "A"),
    (9,      15,     "B"),
    (16,     25,     "C"),
    (26,     50,     "D"),
    (51,     90,     "E"),
    (91,     150,    "F"),
    (151,    280,    "G"),
    (281,    500,    "H"),
    (501,    1200,   "J"),
    (1201,   3200,   "K"),
    (3201,   10000,  "L"),
    (10001,  35000,  "M"),
    (35001,  150000, "N"),
    (150001, 500000, "P"),
    (500001, None,   "Q"),
]

# Table II-A — Single sampling plans for NORMAL inspection (master table),
# sample size by code letter (PDF p.19 = doc p.14).
# @kind: standard
# @units: 1
# @domain: inspection_level=II, plan=single-normal
# [ON-DISK] industrial/mil_std_105e_sampling.pdf @ page=19
#   Table II-A, single sampling plans for normal inspection (master table).
MIL_STD_105E_SAMPLE_SIZE = {
    "A": 2, "B": 3, "C": 5, "D": 8, "E": 13, "F": 20, "G": 32, "H": 50,
    "J": 80, "K": 125, "L": 200, "M": 315, "N": 500, "P": 800, "Q": 1250,
    "R": 2000,
}

# Table II-A acceptance numbers Ac for common AQLs (percent nonconforming);
# Re = Ac + 1 for every direct entry in this range. None = master-table
# ARROW cell ("use first sampling plan below/above the arrow") — templates
# MUST sample only (letter, AQL) pairs with a direct entry (non-None),
# keeping the transcription verbatim and the template logic table-faithful.
# @kind: standard
# @units: 1
# @domain: inspection_level=II, plan=single-normal
# [ON-DISK] industrial/mil_std_105e_sampling.pdf @ page=19
#   Table II-A acceptance numbers Ac (PDF p.19 = doc p.14).
MIL_STD_105E_SINGLE_NORMAL_AC = {
    #        AQL:  0.65   1.0    2.5    4.0    6.5   (percent)
    "A": {0.65: None, 1.0: None, 2.5: None, 4.0: None, 6.5: 0},
    "B": {0.65: None, 1.0: None, 2.5: None, 4.0: 0,   6.5: None},
    "C": {0.65: None, 1.0: None, 2.5: 0,   4.0: None, 6.5: None},
    "D": {0.65: None, 1.0: None, 2.5: None, 4.0: None, 6.5: 1},
    "E": {0.65: None, 1.0: 0,   2.5: None, 4.0: 1,   6.5: 2},
    "F": {0.65: 0,   1.0: None, 2.5: 1,   4.0: 2,   6.5: 3},
    "G": {0.65: None, 1.0: None, 2.5: 2,   4.0: 3,   6.5: 5},
    "H": {0.65: None, 1.0: 1,   2.5: 3,   4.0: 5,   6.5: 7},
    "J": {0.65: 1,   1.0: 2,   2.5: 5,   4.0: 7,   6.5: 10},
    "K": {0.65: 2,   1.0: 3,   2.5: 7,   4.0: 10,  6.5: 14},
    "L": {0.65: 3,   1.0: 5,   2.5: 10,  4.0: 14,  6.5: 21},
    "M": {0.65: 5,   1.0: 7,   2.5: 14,  4.0: 21,  6.5: None},
    "N": {0.65: 7,   1.0: 10,  2.5: 21,  4.0: None, 6.5: None},
    "P": {0.65: 10,  1.0: 14,  2.5: None, 4.0: None, 6.5: None},
    "Q": {0.65: 14,  1.0: 21,  2.5: None, 4.0: None, 6.5: None},
}

# Named measurable quality characteristics for SPC scenarios (Q1-Q2):
# target windows and within-process sigma expressed as a fraction of target.
# Typology per Montgomery ISQC 7e Ch. 6 examples (machined
# dimensions, fill volumes, electrical parameters); parameters sampled, never
# verbatim from worked examples (spec §3 copyright rule).
# @kind: range
# @units: target=per-row(key), sigma_frac=1
# @given: stated (C1 Reviewer G, reviews/phaseC1_reviewer_g_provenance.md §2 row 13 and reviews/phaseC1_reviewer_g_scratch/g_more.py: the grand mean or mu0 stated 40/40 in each of the 4 consumers)
# @domain: none (characteristic windows for scenario generation)
# [POLICY: sampling-only] Sampled values are stated verbatim in the question (the given-values rule), so this window changes which scenario is posed, never whether the answer is right.
SPC_CHARACTERISTICS = {
    "shaft diameter (mm)":        {"target": (10, 80),    "sigma_frac": (0.001, 0.01)},
    "bottle fill volume (mL)":    {"target": (250, 1000), "sigma_frac": (0.002, 0.015)},
    "coating thickness (micron)": {"target": (20, 200),   "sigma_frac": (0.01, 0.05)},
    "resistor resistance (ohm)":  {"target": (100, 5000), "sigma_frac": (0.005, 0.03)},
    "breaking strength (N)":      {"target": (200, 3000), "sigma_frac": (0.01, 0.06)},
}

# Attribute-chart baseline windows (Q3): historical fraction nonconforming
# for p-charts; mean defects per unit for c-charts. Magnitude
# conventions per Montgomery ISQC 7e Ch. 7 examples.
# @kind: range
# @units: 1
# @given: stated (C1 Reviewer G, phaseC1_reviewer_g_provenance.md §2 row 5: D, m and n are stated 40/40 and p-bar = D/(mn) is derived from them)
# @domain: none (a baseline fraction-nonconforming window)
# [POLICY: sampling-only] Sampled values are stated verbatim in the question (the given-values rule), so this window changes which scenario is posed, never whether the answer is right.
P_CHART_PBAR = (0.01, 0.15)
# @kind: range
# @units: 1
# @domain: none (a p-chart subgroup-size window)
# [POLICY: sampling-only] Sampled values are stated verbatim in the question (the given-values rule), so this window changes which scenario is posed, never whether the answer is right.
P_CHART_SUBGROUP_N = (50, 400)
# @kind: range
# @units: 1
# @domain: none (a baseline defects-per-unit window)
# [POLICY: sampling-only] Sampled values are stated verbatim in the question (the given-values rule), so this window changes which scenario is posed, never whether the answer is right.
C_CHART_CBAR = (2.0, 25.0)

# --- Q2 process-capability scenario data (spec R7, Stage D v2) ------
# Moved out of process_capability.py, which inlined all of it.
#
# PROVENANCE CAVEAT — READ BEFORE RELYING ON THESE. Unlike the Montgomery
# and MIL-STD-105E entries elsewhere in this file, the two standards
# below are NOT among the on-disk primary sources for this branch. These
# values are transcribed from the t26R c3 authoring round and were
# reviewed there for internal consistency, but under the branch's
# primary-source-only rule for tabulated data they count as UNVERIFIED.
# They are flagged for the Stage B Data Reviewer and the human
# checkpoint. Nothing here is a graded quantity: the series only supply
# plausible nominal resistances and a plausible thickness envelope.
#
# IEC 60063 preferred values. The E24 series is the 5%-tolerance series
# and E12 the 10% series; pairing a tolerance class with the wrong
# series was a c2 blocking defect, so the pairing is encoded here rather
# than left to the caller. E24 is filtered to multiples of 20 so that
# nominal*tol/100 is an integer and the stem's "nominal +/- tolerance"
# identity is literally true at the stated precision.
# @kind: standard
# @units: ohm
# @domain: none (preferred nominal resistances; no condition applies)
# [UNVERIFIED] IEC 60063 is not on disk. The prose above already says these count as
#   UNVERIFIED under the branch's primary-source-only rule; this makes the claim
#   machine-visible instead of leaving it to a reader. Nothing graded depends on it:
#   the series only supplies plausible nominal resistances, always stated in the item.
IEC60063_E24_5PCT = [n for n in
                     [100, 110, 120, 130, 150, 160, 180, 200, 220, 240,
                      270, 300, 330, 360, 390, 430, 470, 510, 560, 620,
                      680, 750, 820, 910, 1000, 1100, 1200, 1300, 1500,
                      1600, 1800, 2000] if n % 20 == 0]
# @kind: standard
# @units: ohm
# @domain: none (preferred nominal resistances; no condition applies)
# [UNVERIFIED] IEC 60063 is not on disk; see IEC60063_E24_5PCT.
IEC60063_E12_10PCT = [100, 120, 150, 180, 220, 270, 330, 390, 470, 560,
                      680, 820, 1000, 1200, 1500, 1800]
# @kind: standard
# @units: ohm
# @domain: none (the series a tolerance class selects)
# [UNVERIFIED] the pairing of tolerance class to series is IEC 60063's, which is not on
#   disk. Encoded here rather than left to the caller because mispairing them was a c2
#   blocking defect.
RESISTOR_SERIES_BY_TOLERANCE = {5: IEC60063_E24_5PCT,
                                10: IEC60063_E12_10PCT}

# MIL-A-8625 Type III (hard anodize) coating-thickness envelope, microns.
# Screened as a specification envelope, not as a graded value.
# @domain: none (a specification envelope for Type III coatings)
# [UNVERIFIED] MIL-A-8625 is not on disk. Unlike MIL-STD-105E, which is, this one
#   has no artefact to point at, so the envelope is asserted, not read.
# @kind: standard
# @units: um
HARD_ANODIZE_THICKNESS = (20, 110)

# Coating-thickness metrology resolution floor, microns: sigma below this
# is not measurable in practice, so draws under it are rejected.
# @kind: range
# @units: um
# @given: guard (C1 Reviewer G, phaseC1_reviewer_g_provenance.md §2 row 8: a screen, if sigma < FLOOR then redraw; the drawn sigma is stated in the question)
# @domain: none (a metrology resolution floor used as a screen)
# [POLICY: sampling-only] a redraw screen, not a measurement: C3.10's guard detector
#   confirmed statically that its only consumer reads it in `if sigma < FLOOR`, never as
#   a value, so no answer depends on where the floor sits.
COATING_METROLOGY_FLOOR = 1.5

# Subgroup-size windows for variables charts (Q1): the chart-type-selection
# branching template (Q1-#3) uses n <= 10 -> X-bar/R, n > 10 -> X-bar/s,
# per Montgomery ISQC 7e §6.3 guidance ("n moderately large—say, n > 10 or
# 12": the range method loses efficiency for moderate-to-large n).
# @kind: standard
# @units: 1
# @domain: none (subgroup-size windows selecting the chart type)
# [ON-DISK:LOCAL-ONLY] pilot/references/public/full_books_industrial_engineering/Montgomery, Introduction to Statistical Quality Control.pdf @ page=265 text="10 or 12"
#   The page carries "for large n—say, n> 10 or 12—it is probably best to use a control
#   chart for s", which is the n <= 10 / n > 10 split these two windows encode.
XBAR_R_SUBGROUP_N = (2, 10)
# @kind: standard
# @units: 1
# @domain: none (subgroup-size windows selecting the chart type)
# [ON-DISK:LOCAL-ONLY] pilot/references/public/full_books_industrial_engineering/Montgomery, Introduction to Statistical Quality Control.pdf @ page=265 text="10 or 12"
#   The upper half of the same split.
XBAR_S_SUBGROUP_N = (11, 25)
# @kind: range
# @units: 1
# @given: stated (C1 Reviewer G, phaseC1_reviewer_g_provenance.md §2 row 9: m is stated 40/40 in the p chart and in the c-chart question; the c chart reads the table only in an assert)
# @domain: none (a preliminary-sample-count window for trial limits)
# [POLICY: sampling-only] Sampled values are stated verbatim in the question (the given-values rule), so this window changes which scenario is posed, never whether the answer is right.
SPC_NUM_SUBGROUPS = (20, 30)   # preliminary samples m for trial limits
