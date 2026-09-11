"""Reviewer G, C1: independent check that each sampled PLAUSIBILITY table's
drawn value is STATED in the question (or the table is read only as a guard).
Runs templates from the frozen worktree; patches only in-memory module globals.
Usage: python g_stated.py   (cwd = worktree)
"""
import sys, os, random, re, inspect, importlib
WT = os.getcwd()
sys.path.insert(0, WT)
from tests.template_integrity.core import discover

REFS = {r.template_id: r for r in discover()}
SEEDS = range(40)


def gen(fn, seed):
    random.seed(seed)
    try:
        return fn()
    except Exception as exc:  # noqa
        return None, f"ERROR {type(exc).__name__}: {exc}"


def fn_of(tid):
    return REFS[tid].load()


def mod_of(tid):
    return sys.modules[fn_of(tid).__module__]


def consumers_in_source(name):
    """template ids whose FUNCTION source mentions name (my own, not census)."""
    out = []
    for tid, r in REFS.items():
        try:
            src = inspect.getsource(r.load())
        except Exception:
            continue
        if re.search(r"\b%s\b" % name, src):
            out.append(tid)
    return sorted(out)


def check(label, tid, extract, ok):
    fn = fn_of(tid)
    stated = inwin = err = 0
    bad = []
    for s in SEEDS:
        q, sol = gen(fn, s)
        if q is None:
            err += 1
            continue
        vals = extract(q)
        if vals is None:
            bad.append((s, "NOT-STATED"))
            continue
        stated += 1
        if ok(vals):
            inwin += 1
        else:
            bad.append((s, vals))
    print(f"[{label}] {tid}: stated {stated}/40, in-window {inwin}/40, errors {err}; bad={bad[:3]}")


C = importlib.import_module("data.templates.branches.civil_engineering.constants")
I = importlib.import_module("data.templates.branches.industrial_engineering.constants")
E = importlib.import_module("data.templates.branches.electrical_engineering.constants")

print("== my own function-source consumer scan (call-time reads only) ==")
for n in ["SPECIFIC_GRAVITY_RANGES", "FRICTION_ANGLE_RANGES_DEG", "PERMEABILITY_RANGES_CM_S",
          "SERVICE_LEVELS", "P_CHART_PBAR", "P_CHART_SUBGROUP_N", "C_CHART_CBAR",
          "COATING_METROLOGY_FLOOR", "SPC_NUM_SUBGROUPS", "QUEUE_SCENARIOS",
          "NEWSVENDOR_ITEMS", "COMPONENT_RELIABILITY_CLASSES", "SPC_CHARACTERISTICS",
          "AMPLITUDE_RANGE", "PHASE_RANGE_RAD"]:
    print(f"  {n}: {consumers_in_source(n)}")

# ---- 1. SPECIFIC_GRAVITY_RANGES: Gs stated and inside its soil window
GS = re.compile(r"specific gravity of (?:the )?(?:soil )?solids (?:is|of) (\d+\.\d+)")
SOIL = re.compile(r"\b(inorganic clay|silt|sand)\b")
def gs_ex(q):
    m = GS.search(q); s = SOIL.search(q)
    return (float(m.group(1)), s.group(1)) if m and s else None
def gs_ok(v):
    lo, hi = C.SPECIFIC_GRAVITY_RANGES[v[1]]
    return lo - 1e-9 <= v[0] <= hi + 1e-9
for tid in ["template_phase_relations_degree_of_saturation", "template_relative_density_of_sand",
            "template_borrow_pit_fill_volume", "template_effective_stress_profile",
            "template_upward_seepage_quick_condition"]:
    check("SG", tid, gs_ex, gs_ok)

# ---- 2. FRICTION_ANGLE_RANGES_DEG
FA = re.compile(r"cohesionless (\w+) \(([^)]+)\).*friction angle of the soil is (\d+) degrees", re.S)
def fa_ex(q):
    m = FA.search(q)
    return (f"{m.group(1)}, {m.group(2)}", int(m.group(3))) if m else None
def fa_ok(v):
    lo, hi = C.FRICTION_ANGLE_RANGES_DEG[v[0]]
    return lo <= v[1] <= hi
check("FA", "template_infinite_slope_factor_of_safety", fa_ex, fa_ok)

# ---- 3. PERMEABILITY_RANGES_CM_S: guard-only test (patch the template module's binding)
tid = "template_constant_head_permeability"; fn = fn_of(tid); M = mod_of(tid)
base = [gen(fn, s) for s in SEEDS]
orig = M.PERMEABILITY_RANGES_CM_S
M.PERMEABILITY_RANGES_CM_S = {k: (0.0, 1e9) for k in orig}
wide = [gen(fn, s) for s in SEEDS]
M.PERMEABILITY_RANGES_CM_S = {k: (1e-9, 1e-8) for k in orig}
narrow = [gen(fn, s) for s in SEEDS]
M.PERMEABILITY_RANGES_CM_S = orig
print(f"[PERM] {tid}: widened window -> identical items {sum(a == b for a, b in zip(base, wide))}/40; "
      f"window moved off the draws -> AssertionError {sum(1 for x in narrow if x[0] is None)}/40")

# ---- 4. SERVICE_LEVELS: stated? and what a numeric nudge does
tid = "template_safety_stock_reorder_point"; fn = fn_of(tid); M = mod_of(tid)
SL = re.compile(r"targets a (\d+(?:\.\d+)?)% cycle-service level")
check("SL", tid, lambda q: (float(SL.search(q).group(1)) / 100) if SL.search(q) else None,
      lambda v: any(abs(v - a) < 1e-9 for a in I.SERVICE_LEVELS))
orig = M.SERVICE_LEVELS
M.SERVICE_LEVELS = [round(a * 1.001, 6) for a in orig]
nud = [gen(fn, s) for s in SEEDS]
M.SERVICE_LEVELS = orig
print(f"[SL] nudged levels -> errors {sum(1 for x in nud if x[0] is None)}/40; first: {next((x[1] for x in nud if x[0] is None), '-')}")

# ---- 5. p-chart: m, n, D stated; p-bar = D/(mn) inside P_CHART_PBAR
PC = re.compile(r"D = (\d+) nonconforming \w+ in m = (\d+) samples of n = (\d+)")
def pc_ex(q):
    m = PC.search(q)
    return tuple(int(x) for x in m.groups()) if m else None
def pc_ok(v):
    d, m, n = v
    return (I.SPC_NUM_SUBGROUPS[0] <= m <= I.SPC_NUM_SUBGROUPS[1]
            and I.P_CHART_SUBGROUP_N[0] <= n <= I.P_CHART_SUBGROUP_N[1]
            and I.P_CHART_PBAR[0] <= d / (m * n) <= I.P_CHART_PBAR[1])
check("PCHART", "template_p_chart_limits_floor", pc_ex, pc_ok)

# ---- 6. reliability: every component reliability printed
RL = re.compile(r"reliabilit(?:y|ies) (?:of )?((?:0\.\d{3}[^.]*?)+)\.")
def rl_ex(q):
    v = [float(x) for x in re.findall(r"\b0\.\d{3}\b", q)]
    return v if len(v) >= 2 else None
check("REL", "template_system_reliability" if "template_system_reliability" in REFS else
      consumers_in_source("COMPONENT_RELIABILITY_CLASSES")[0], rl_ex,
      lambda v: all(0.90 <= x <= 0.99 for x in v))

# ---- 7. queue M/M/1 and M/M/c: rates / servers stated
Q1 = re.compile(r"(\d+) customers per hour")
check("QUEUE", "template_mm1_time_in_system",
      lambda q: [int(x) for x in Q1.findall(q)] if len(Q1.findall(q)) >= 2 else None,
      lambda v: all(0 < x <= 90 for x in v))
QC = re.compile(r"average rate of (\d+) customers per hour\. The facility has (\d+) identical servers")
check("QUEUE", "template_mmc_waiting_time",
      lambda q: tuple(int(x) for x in QC.search(q).groups()) if QC.search(q) else None,
      lambda v: 1 <= v[1] <= 6)

# ---- 8. newsvendor: c, p, s stated
NV = re.compile(r"costs \$(\d+\.\d\d), sells for \$(\d+\.\d\d).*salvaged for \$(\d+\.\d\d)", re.S)
check("NEWS", "template_newsvendor_normal_demand",
      lambda q: tuple(float(x) for x in NV.search(q).groups()) if NV.search(q) else None,
      lambda v: v[2] < v[0] < v[1])

# ---- 9. print one item from the consumers I have not read, for eyeballing
def show(tid, seed=3, want=None):
    fn = fn_of(tid)
    for s in range(seed, seed + 60):
        q, sol = gen(fn, s)
        if q is None:
            continue
        if want and want not in q:
            continue
        print(f"\n--- {tid} seed={s}\nQ: {q}\nS[:900]: {sol[:900]}")
        return
    print(f"\n--- {tid}: no matching item")

for tid in consumers_in_source("SPC_CHARACTERISTICS"):
    show(tid)
for tid in consumers_in_source("C_CHART_CBAR"):
    show(tid)
for tid in consumers_in_source("COATING_METROLOGY_FLOOR"):
    show(tid, want="microns")
for tid in sorted(set(consumers_in_source("AMPLITUDE_RANGE") + consumers_in_source("PHASE_RANGE_RAD"))):
    show(tid)
