"""D5.6 -- the candidate extraction rules, scored, with the losers kept.

    python -m tests.comparators.n1_candidates            # rescore everything
    python -m tests.comparators.n1_candidates --check    # fail if the winner moved

N1 is that ``parse_number`` read the **first** number in an answer span, which is
often a fluid grade (``Engine Oil (SAE 50)``), a temperature (``at 541 K``) or a
chemical-formula subscript (the 4 of ``C4H10``). Picking a replacement from one
example is the error shape ``phase4_summary.md`` |S|8 records **six times**, so the
replacement was chosen by measurement over four corpora, and this module is what
makes that claim checkable rather than assertable.

**The held-out split is frozen before any candidate is scored.** ``held_out()``
is a pure function of the template id and a fixed salt. It does not look at any
result and it is defined above the candidates, which is the only ordering that
means anything.

## The two axes, and why one alone is worthless

**(a) false accepts** on gold x gold and archive x gold. Alone, the winner is
``unique_or_unresolved``: it achieves zero false accepts by deciding nothing.

**(b) decided rate on answers that ought to be decided** -- D4.4's 24 ``numeric``
cases, and every archive row paired with its own gold. The second needs no truth
label at all, which matters: labelling it from the deployed pipeline's
``final_answer_acc`` would be building the replacement to agree with the parser
it replaces (D-003, D-034).

## Three findings that no single example would have produced

1. **"Last number" is far worse than "first", not better** -- 424 held-out gold
   x gold false accepts against the incumbent's 16.
2. **The tokeniser matters more than the choice of number.** Excluding digits
   inside ``m^2`` and ``C4H10`` is orthogonal to which of the remaining numbers a
   rule picks, and it is what turns "last number" from 424 false accepts into 0.
3. **The rule must be per-kind.** ``check`` answers state the quantity first and
   the threshold second -- *"the deflection is 18.4 mm, less than the 25 mm
   limit"* -- so every last-ward rule picks the limit. On D4.4's 21 ``check``
   cases, first-number scores 13 decided / 13 correct; every last-ward rule
   scores 13 decided / **5** correct / 8 false rejects.

## A confound that had to be separated before the table meant anything

Two false accepts in the dev slice were ``400.0 MHz`` against ``400.0 kHz`` --
the **same number, a different unit**. No extraction rule can fix that; only a
declared unit can (D-052). Counted together with extraction errors it would have
credited or blamed every rule for something outside its reach, so the two are
reported in separate columns.

## And one where the incumbent's advantage was the defect

The incumbent's archive decided rate looked higher because it decided on a
*non-exponent incidental* number whose precision was computable, while the true
answer was often in exponent form and ``_decimals`` returned ``None`` for those.
So the incumbent decided **by reading the wrong number**. Implementing |S|7.1's
display tolerance for scientific notation (``extract.displayed_decimals``) lifted
the archive decided rate under *every* rule including the incumbent, which is
what shows it to be an independent fix rather than a way of paying for this one.
"""

from __future__ import annotations

import argparse
import csv
import glob
import hashlib
import json
import os
import re
import sys
from collections import Counter, defaultdict
from decimal import Decimal, InvalidOperation

REPO = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", ".."))
if REPO not in sys.path:
    sys.path.insert(0, REPO)

from tests.comparators import extract as EX, kinds  # noqa: E402
from tests.comparators.answer import compare_kind  # noqa: E402
from tests.comparators.bindings import DECLARED_UNITS  # noqa: E402
from tests.comparators.normalize import answer_span  # noqa: E402
from tests.template_integrity.core import discover, generate  # noqa: E402

INVENTORY = os.path.join(REPO, "docs", "re-implementation-sep", "audit", "template_inventory.csv")
ADVERSARIAL = os.path.join(os.path.dirname(__file__), "adversarial.json")
RESULTS = os.path.join(os.path.dirname(__file__), "n1_candidates.json")

# --------------------------------------------------------------------------
# THE SPLIT.  Frozen here, above every candidate, before any result exists.
# --------------------------------------------------------------------------
SALT = "phase5-d5.6-holdout-v1"


def held_out(template_id: str) -> bool:
    h = hashlib.sha256((SALT + template_id).encode()).hexdigest()
    return int(h[:8], 16) % 5 == 0          # 1 in 5 == 20%


# --------------------------------------------------------------------------
# Candidates.  Each returns the chosen match, or None meaning "I refuse".
# --------------------------------------------------------------------------
NUM = EX.NUM_RE


def c1_first(t, u=None):
    return NUM.search(t)


def c2_last(t, u=None):
    ms = list(NUM.finditer(t))
    return ms[-1] if ms else None


def c3_unit_adjacent(t, u=None):
    return EX._unit_adjacent(t, u)


def c4_unique_or_unresolved(t, u=None):
    ms = list(NUM.finditer(t))
    return ms[0] if len(ms) == 1 else None


def c5_anchored(t, u=None):
    return EX._after_anchor(t)


def c6_unit_then_last(t, u=None):
    return EX._unit_adjacent(t, u) or c2_last(t)


def c7_anchored_unit_last(t, u=None):
    return EX._after_anchor(t) or EX._unit_adjacent(t, u) or c2_last(t)


def c8_unit_then_anchor(t, u=None):
    return EX._unit_adjacent(t, u) or EX._after_anchor(t)


def c9_filtered_last(t, u=None):
    return EX._last_filtered(t)


def c10_ADOPTED(t, u=None):
    """unit -> anchor -> filtered last.  What `extract.answer_match` does."""
    return EX.answer_match(t, "numeric", u)


def c11_unit_or_single(t, u=None):
    hit = EX._unit_adjacent(t, u)
    if hit:
        return hit
    ms = EX.numbers(t)
    return ms[0] if len(ms) == 1 else None


CANDIDATES = {
    "C1_first_INCUMBENT": c1_first,
    "C2_last": c2_last,
    "C3_unit_adjacent": c3_unit_adjacent,
    "C4_unique_or_unresolved": c4_unique_or_unresolved,
    "C5_anchored": c5_anchored,
    "C6_unit_then_last": c6_unit_then_last,
    "C7_anchored_unit_last": c7_anchored_unit_last,
    "C8_unit_then_anchor": c8_unit_then_anchor,
    "C9_filtered_last": c9_filtered_last,
    "C10_ADOPTED": c10_ADOPTED,
    "C11_unit_or_single": c11_unit_or_single,
    # The configuration that ships: C10 for `numeric`, first-number for `check`.
    "SHIPPED_per_kind": c10_ADOPTED,
}
#: The one entry scored with per-kind routing rather than a single rule.
PER_KIND = {"SHIPPED_per_kind"}
ADOPTED = "SHIPPED_per_kind"

# --------------------------------------------------------------------------
# Patching.  `parse_number` and `_decimals` move TOGETHER, always: deriving the
# value from one number and the tolerance from another would judge an answer
# against a precision it never had.
# --------------------------------------------------------------------------
_CUR = {"pick": c1_first, "unit": None}


def _pick(text, kind):
    """The candidate rule under test, applied to EVERY kind.

    Deliberately NOT routed per-kind here.  The point of the `check` column is
    to show what each rule does to `check` answers when it is the only rule --
    which is the measurement that justified routing them separately at all.
    `SHIPPED_per_kind` is the one entry scored with the routing that ships.
    """
    if _CUR.get("per_kind"):
        return EX.answer_match(text, kind, _CUR["unit"])
    return _CUR["pick"](text, _CUR["unit"])


def _patched_parse(text, kind="numeric", unit=None):
    m = _pick(text, kind)
    if m is None:
        return None
    try:
        return Decimal(m.group(0).replace(",", ""))
    except InvalidOperation:
        return None


def _patched_decimals(text, kind="numeric", unit=None):
    m = _pick(text, kind)
    return None if m is None else EX.displayed_decimals(m.group(0))


def _strip_num(s):
    return re.sub(r"[-+]?[\d,]*\.?\d+(?:[eE][-+]?\d+)?", "#", s)


def _same_number_different_unit(sa, sb):
    na = [m.group(0).replace(",", "") for m in NUM.finditer(sa)]
    nb = [m.group(0).replace(",", "") for m in NUM.finditer(sb)]
    return na == nb and _strip_num(sa) != _strip_num(sb)


def nspan(t):
    return re.sub(r"\s+", " ", answer_span(t)[0]).strip()


# --------------------------------------------------------------------------

def build_corpora(n=12):
    inv = {r["template_id"]: r for r in csv.DictReader(open(INVENTORY, encoding="utf-8"))}
    gg = {}
    for ref in discover():
        if inv[ref.template_id]["answer_type"] != "scalar":
            continue
        sols = []
        for s in range(n):
            i = generate(ref, s, capture=False)
            if not i.ok:
                sols = []
                break
            sols.append(i.solution)
        if sols:
            gg[ref.template_id] = {
                "sols": sols, "spans": [nspan(x) for x in sols],
                "unit": DECLARED_UNITS.get(ref.template_id)}
    ids = {r.template_id for r in discover()}
    ag = defaultdict(list)
    for p in sorted(glob.glob(os.path.join(REPO, "error_analysis_annotation",
                                           "samples", "*.jsonl"))):
        for line in open(p, encoding="utf-8"):
            rec = json.loads(line)
            stem, inst = rec["question_id"].rsplit("__", 1)
            tid = "template_" + stem
            if tid in ids and inv[tid]["answer_type"] == "scalar":
                ag[tid].append((inst, rec["model_reasoning"], rec["gold_answer"]))
    d44 = [c for c in json.load(open(ADVERSARIAL, encoding="utf-8"))]
    return gg, ag, d44


def score(pick, gg, ag, d44):
    r = Counter()
    for tid, d in gg.items():
        _CUR["pick"], _CUR["unit"] = pick, d["unit"]
        sols, spans = d["sols"], d["spans"]
        b = "HELD" if held_out(tid) else "DEV"
        for a in range(len(sols)):
            for x in range(len(sols)):
                if a == x:
                    continue
                try:
                    v = compare_kind("numeric", sols[x], sols[a])
                except Exception:                              # noqa: BLE001
                    r[f"gg_{b}_ERROR"] += 1
                    continue
                r[f"gg_{b}_{v.outcome}"] += 1
                r[f"gg_{b}_decided"] += (v.outcome != "UNRESOLVED")
                if v.outcome == "MATCH" and spans[a] != spans[x]:
                    if _same_number_different_unit(spans[a], spans[x]):
                        r[f"gg_{b}_FA_unit"] += 1
                    else:
                        r[f"gg_{b}_FA"] += 1
    for tid, rs in ag.items():
        _CUR["pick"], _CUR["unit"] = pick, DECLARED_UNITS.get(tid)
        b = "HELD" if held_out(tid) else "DEV"
        golds = {}
        for inst, _c, g in rs:
            golds.setdefault(inst, g)
        for inst, cand, g in rs:
            try:
                v = compare_kind("numeric", g, cand)
                r[f"ag_{b}_pos"] += 1
                r[f"ag_{b}_pos_decided"] += (v.outcome != "UNRESOLVED")
            except Exception:                                  # noqa: BLE001
                r[f"ag_{b}_pos_ERROR"] += 1
            for inst2, g2 in golds.items():
                if inst2 == inst or nspan(g2) == nspan(g):
                    continue
                r[f"ag_{b}_neg"] += 1
                try:
                    if compare_kind("numeric", g2, cand).outcome == "MATCH":
                        r[f"ag_{b}_XI"] += 1
                except Exception:                              # noqa: BLE001
                    r[f"ag_{b}_ERROR"] += 1
    for c in d44:
        if c["kind"] not in ("numeric", "check"):
            continue
        opts = eval(c["options"]) if isinstance(c["options"], str) else c["options"]
        _CUR["pick"], _CUR["unit"] = pick, opts.get("unit")
        tag = c["kind"]
        try:
            v = compare_kind(c["kind"], c["gold"], c["candidate"], **opts)
        except Exception:                                      # noqa: BLE001
            r[f"d44_{tag}_ERROR"] += 1
            continue
        r[f"d44_{tag}_n"] += 1
        if v.outcome == "UNRESOLVED":
            continue
        r[f"d44_{tag}_decided"] += 1
        want, got = str(c["correct"]) == "True", v.outcome == "MATCH"
        if got == want:
            r[f"d44_{tag}_correct"] += 1
        elif got:
            r[f"d44_{tag}_FA"] += 1
        else:
            r[f"d44_{tag}_FR"] += 1
    return r


def _row(name, r, b):
    ggn = sum(r[f"gg_{b}_{k}"] for k in ("MATCH", "MISMATCH", "UNRESOLVED")) or 1
    pos = r[f"ag_{b}_pos"] or 1
    return (f"{name:26s} {r[f'gg_{b}_FA']:6d} {r[f'gg_{b}_FA_unit']:6d} "
            f"{100 * r[f'gg_{b}_decided'] / ggn:7.1f} {r[f'ag_{b}_XI']:5d} "
            f"{100 * r[f'ag_{b}_pos_decided'] / pos:7.1f}")


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--check", action="store_true",
                    help="fail if the adopted rule is no longer a best scorer")
    args = ap.parse_args(argv)

    orig = (kinds.parse_number, kinds._decimals)
    kinds.parse_number, kinds._decimals = _patched_parse, _patched_decimals
    try:
        gg, ag, d44 = build_corpora()
        n_h = sum(1 for t in gg if held_out(t))
        print(f"gold x gold   {len(gg)} scalar templates "
              f"({len(gg) - n_h} dev / {n_h} held out)")
        n_ah = sum(1 for t in ag if held_out(t))
        print(f"archive x gold {len(ag)} covered ({len(ag) - n_ah} dev / {n_ah} held out), "
              f"{sum(len(v) for v in ag.values())} rows")
        print(f"split salt {SALT!r}, frozen before any candidate was written")
        print()
        res = {}
        for name, fn in CANDIDATES.items():
            _CUR["per_kind"] = name in PER_KIND
            res[name] = score(fn, gg, ag, d44)
        _CUR["per_kind"] = False
        hdr = (f"{'rule':26s} {'ggFA':>6s} {'ggUNIT':>6s} {'ggDEC%':>7s} "
               f"{'agXI':>5s} {'agDEC%':>7s}")
        for b, title in (("DEV", "DEV SET"), ("HELD", "HELD-OUT SET (the one that counts)")):
            print(title)
            print(hdr)
            for name in CANDIDATES:
                print(_row(name, res[name], b))
            print()
        print(f"{'rule':26s} {'num dec':>8s} {'num FA':>7s} "
              f"{'chk dec':>8s} {'chk ok':>7s} {'chk FR':>7s}   D4.4")
        for name in CANDIDATES:
            r = res[name]
            print(f"{name:26s} {r['d44_numeric_decided']:3d}/{r['d44_numeric_n']:<4d} "
                  f"{r['d44_numeric_FA']:7d} {r['d44_check_decided']:3d}/{r['d44_check_n']:<4d} "
                  f"{r['d44_check_correct']:7d} {r['d44_check_FR']:7d}")
        json.dump({"salt": SALT,
                   "held_out_gg": sorted(t for t in gg if held_out(t)),
                   "held_out_ag": sorted(t for t in ag if held_out(t)),
                   "adopted": ADOPTED,
                   "results": {k: dict(v) for k, v in res.items()}},
                  open(RESULTS, "w", encoding="utf-8"), indent=1)
    finally:
        kinds.parse_number, kinds._decimals = orig

    a = res[ADOPTED]
    zero_fa = [n for n in CANDIDATES if res[n]["gg_HELD_FA"] == 0]
    best_dec = max(res[n]["ag_HELD_pos_decided"] for n in zero_fa)
    print()
    print(f"ADOPTED {ADOPTED}: held-out gold x gold false accepts "
          f"{a['gg_HELD_FA']}, archive decided {a['ag_HELD_pos_decided']}")
    print(f"  candidates with zero held-out false accepts: {len(zero_fa)}")
    print(f"  best archive decided among them: {best_dec}")
    ok = a["gg_HELD_FA"] == 0 and a["ag_HELD_pos_decided"] == best_dec
    print("  adopted rule is a best scorer among the zero-false-accept set:",
          "YES" if ok else "NO")
    if args.check and not ok:
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
