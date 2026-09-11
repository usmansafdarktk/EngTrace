"""Reviewer G, C1: follow-ups. cwd = worktree. In-memory patches only."""
import sys, os, random, re, importlib, inspect
sys.path.insert(0, os.getcwd())
from tests.template_integrity.core import discover

REFS = {r.template_id: r for r in discover()}


def gen(fn, s):
    random.seed(s)
    try:
        return fn()
    except Exception as exc:  # noqa
        return None, f"ERROR {type(exc).__name__}: {exc}"


def cons(name):
    out = []
    for tid, r in REFS.items():
        try:
            if re.search(r"\b%s\b" % name, inspect.getsource(r.load())):
                out.append(tid)
        except Exception:
            pass
    return sorted(out)


C = importlib.import_module("data.templates.branches.civil_engineering.constants")
I = importlib.import_module("data.templates.branches.industrial_engineering.constants")

# ---- A. SCS_IA_RATIO: consumed by literal copy?
print("SCS_CURVE_NUMBERS readers:", cons("SCS_CURVE_NUMBERS"), "| SCS_IA_RATIO readers:", cons("SCS_IA_RATIO"))
tid = cons("SCS_CURVE_NUMBERS")[0]
fn = REFS[tid].load(); modname = fn.__module__
base = [gen(fn, s) for s in range(40)]
print(f"[SCS] {tid}: 'Ia = 0.2S' in question {sum('Ia = 0.2S' in q for q, _ in base if q)}/40; "
      f"'Ia = 0.2*S = 0.2 *' in solution {sum('Ia = 0.2*S = 0.2 *' in s for _, s in base)}/40")
print("[SCS] seed 0 Q:", base[0][0])
orig = C.SCS_IA_RATIO
C.SCS_IA_RATIO = 0.3
mod = importlib.reload(sys.modules[modname])
pat = [gen(getattr(mod, tid), s) for s in range(40)]
C.SCS_IA_RATIO = orig
importlib.reload(sys.modules[modname])
print(f"[SCS] SCS_IA_RATIO set to 0.3 and module reloaded -> identical items {sum(a == b for a, b in zip(base, pat))}/40")

# ---- B. census probe on SERVICE_LEVELS: what do the counts say?
from tests.constants_integrity import census as K
tid = "template_safety_stock_reorder_point"; fn = REFS[tid].load()
bl = {tid: [K._gen(fn, s) for s in range(10)]}
for fld in K.fields(I.SERVICE_LEVELS):
    res = K.probe(I, "SERVICE_LEVELS", fld, {tid: fn.__module__}, 10, bl)
    print(f"[PROBE] SERVICE_LEVELS field={fld!r}: {res[tid]} -> field_verdict {K.field_verdict(res[tid])}")
print("[PROBE] classify('range','NO-EFFECT') ->", K.classify("range", "NO-EFFECT"))
print("[PROBE] classify('range','SOME-HIDDEN') ->", K.classify("range", "SOME-HIDDEN"))

# ---- C. 40-seed stated checks
def stated(tid, rx, lo, hi):
    fn = REFS[tid].load(); st = inw = err = 0
    for s in range(40):
        q, sol = gen(fn, s)
        if q is None:
            err += 1; continue
        vals = [float(v) for v in re.findall(rx, q)]
        if vals:
            st += 1
            inw += all(lo <= v <= hi for v in vals)
    print(f"[STATED] {tid}: stated {st}/40, all in [{lo},{hi}] {inw}/40, errors {err}")

allt = [w["target"] for w in I.SPC_CHARACTERISTICS.values()]
tlo, thi = min(a for a, b in allt), max(b for a, b in allt)
for tid in cons("SPC_CHARACTERISTICS"):
    stated(tid, r"(?:x-double-bar|mu0|mean of mu) = (\d+(?:\.\d+)?)", tlo, thi)
E = importlib.import_module("data.templates.branches.electrical_engineering.constants")
for tid in cons("AMPLITUDE_RANGE"):
    stated(tid, r"(\d+\.\d+) \* (?:sin|cos)", *E.AMPLITUDE_RANGE)

tid = "template_continuous_to_discrete_conversion"; fn = REFS[tid].load()
rad = deg = other = 0
for s in range(40):
    q, _ = gen(fn, s)
    if q is None:
        continue
    if re.search(r"\+ -?\d+ deg\)", q):
        deg += 1
    elif re.search(r"\+ -?\d+\.\d+\)", q) or re.search(r"- \d+\.\d+\)", q):
        rad += 1
    else:
        other += 1
print(f"[STATED] {tid}: phase stated in deg {deg}, in rad {rad}, no phase/other {other} (of 40)")
