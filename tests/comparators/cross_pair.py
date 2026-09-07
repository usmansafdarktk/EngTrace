"""D5.4 / D5.5 -- cross-pairing as a standing check, not a review artefact.

    python -m tests.comparators.cross_pair              # both corpora
    python -m tests.comparators.cross_pair --gold-gold  # D5.4 only
    python -m tests.comparators.cross_pair --archive    # D5.5 only
    python -m tests.comparators.cross_pair --json out.json

Phase 4 closed with two gate items **met and uninformative**: the archive holds
no wrong answer at all for ``categorical`` or ``categorical[tuple]``, and no
trace of any kind for ``numeric``, ``check`` or ``narrative``. "No false accept
on real archived traces" therefore rested on 15 negative instances in two of six
kinds. Cross-pairing is what fills that in, and Reviewer E said of its own that
it *"should be run as a standing check rather than as a one-off review
artefact"*. This is that check.

## D5.4 -- gold x gold

Gold is by definition the correct answer to its own question, so pairing gold
*A* as a **candidate** against gold *B* as the **gold** has a truth known with no
label at all. That manufactures negatives for every kind, at arbitrary scale,
from the templates alone -- including for the three kinds the archive cannot
reach.

**Construction, stated so it can be checked rather than trusted:**

* for each of the 150 templates, generate ``n`` instances from seeds ``0..n-1``;
* take every **ordered** pair ``(a, b)`` with ``a != b`` -- ``n*(n-1)`` per
  template;
* score ``compare_template(t, gold=b, candidate=a)`` under ``t``'s declared
  binding;
* **truth**: the pair is expected to MATCH iff the two gold answer spans are the
  same text, after whitespace normalisation. A MATCH on a pair whose spans
  differ is a **false accept**.

The textual-identity proxy can err in one direction -- two instances could carry
the same answer written differently -- so **every false accept is reported with
both spans**, and the proxy is audited rather than believed.

## D5.5 -- archive x gold

Gold x gold cannot test *normalisation*, because gold text carries none of the
surface variation real model output does. So the same test runs again with a
real archived model answer as the candidate.

**Two constructions, and they answer different questions:**

* **positive** -- candidate = a model's answer for instance *i*, gold = the gold
  of instance *i*. Truth-free by design: the question is the **decided rate**,
  not correctness, so no label from the deployed pipeline is needed. Using
  ``final_answer_acc`` here would be building a binding to agree with the
  parser being replaced (D-003, D-034).
* **negative** -- candidate = a model's answer for instance *i*, gold = the gold
  of instance *j != i* whose answer differs. A MATCH is a **cross-instance
  accept**. It is *suspicious rather than certainly wrong* -- a model that
  answered *i* incorrectly could have produced *j*'s answer -- so these are
  reported separately and never folded into the false-accept count.

## What every run reports, and why each column is there

``pairs / templates-scored / templates-skipped``, **with every skip named.**
``19,668`` reads as complete and is not: it is 149 of 150, and the missing one
crashed. A pair count without its denominator in *templates* hides a crash as a
pass, which is exactly how N4 survived the sweep that found N1-N3.

**Errors are counted separately from verdicts** (D5.7b). A template that raises
has not passed. ``continuous_to_discrete_conversion`` contributed zero verdicts
and satisfied "zero false accepts" by crashing on all 132 of its pairs.

**Per-kind negative-instance counts sit beside every rate** (D5.9). Reviewer E:
*"precision without it is unreadable."* A kind with two negatives and a kind
with twelve thousand both report 100%.
"""

from __future__ import annotations

import argparse
import glob
import json
import os
import re
import sys
from collections import Counter, defaultdict

REPO = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", ".."))
if REPO not in sys.path:
    sys.path.insert(0, REPO)

from tests.comparators.bindings import BINDINGS, UNBOUND, compare_template  # noqa: E402
from tests.comparators.normalize import answer_span  # noqa: E402
from tests.template_integrity.core import discover, generate  # noqa: E402

ARCHIVE = os.path.join(REPO, "error_analysis_annotation", "samples", "*.jsonl")
DEFAULT_N = 12


def nspan(text: str) -> str:
    return re.sub(r"\s+", " ", answer_span(text)[0]).strip()


# --------------------------------------------------------------------------
# D5.4
# --------------------------------------------------------------------------

def gold_gold(n: int = DEFAULT_N) -> dict:
    scored, skipped = {}, {}
    tally = Counter()
    per_kind = defaultdict(Counter)
    false_accepts = []
    for ref in discover():
        tid = ref.template_id
        if tid not in BINDINGS:
            skipped[tid] = UNBOUND.get(tid, "no declared binding")
            continue
        sols = []
        gen_err = None
        for s in range(n):
            inst = generate(ref, s, capture=False)
            if not inst.ok:
                gen_err = f"seed {s}: {inst.error}"
                break
            sols.append(inst.solution)
        if gen_err:
            skipped[tid] = f"generation failed: {gen_err}"
            continue
        spans = [nspan(x) for x in sols]
        kind = BINDINGS[tid]["kind"]
        r = Counter()
        for a in range(n):
            for b in range(n):
                if a == b:
                    continue
                r["pairs"] += 1
                same = spans[a] == spans[b]
                r["negatives"] += (not same)
                try:
                    v = compare_template(tid, sols[b], sols[a])
                except Exception as exc:                      # noqa: BLE001
                    r["errors"] += 1
                    r.setdefault("error_example", None)
                    if r["errors"] == 1:
                        r["_err"] = f"{type(exc).__name__}: {exc}"
                    continue
                r[v.outcome] += 1
                if v.outcome == "MATCH" and not same:
                    r["false_accepts"] += 1
                    if len(false_accepts) < 25:
                        false_accepts.append(
                            {"template": tid, "kind": kind,
                             "cand": spans[a][:110], "gold": spans[b][:110]})
                elif v.outcome != "MATCH" and same:
                    r["false_rejects"] += 1
        scored[tid] = {"kind": kind, **{k: v for k, v in r.items()
                                        if not str(k).startswith("_")}}
        if "_err" in r:
            scored[tid]["error_example"] = r["_err"]
        for k in ("pairs", "negatives", "MATCH", "MISMATCH", "UNRESOLVED",
                  "errors", "false_accepts", "false_rejects"):
            tally[k] += r[k]
            per_kind[kind][k] += r[k]
    return {"n": n, "scored": scored, "skipped": skipped, "tally": dict(tally),
            "per_kind": {k: dict(v) for k, v in per_kind.items()},
            "false_accept_examples": false_accepts}


# --------------------------------------------------------------------------
# D5.5
# --------------------------------------------------------------------------

def load_archive() -> tuple[dict, dict]:
    ids = {r.template_id for r in discover()}
    rows = defaultdict(list)
    unmapped = Counter()
    for path in sorted(glob.glob(ARCHIVE)):
        with open(path, encoding="utf-8") as fh:
            for line in fh:
                rec = json.loads(line)
                stem, inst = rec["question_id"].rsplit("__", 1)
                tid = "template_" + stem
                if tid not in ids:
                    unmapped[stem] += 1
                    continue
                rows[tid].append((inst, rec["model_reasoning"], rec["gold_answer"]))
    return rows, dict(unmapped)


def archive_gold() -> dict:
    rows, unmapped = load_archive()
    scored, skipped = {}, {}
    tally = Counter()
    per_kind = defaultdict(Counter)
    xi_examples = []
    for tid, rs in sorted(rows.items()):
        if tid not in BINDINGS:
            skipped[tid] = UNBOUND.get(tid, "no declared binding")
            continue
        kind = BINDINGS[tid]["kind"]
        golds = {}
        for inst, _c, g in rs:
            golds.setdefault(inst, g)
        r = Counter()
        for inst, cand, g in rs:
            r["rows"] += 1
            try:                                      # positive: its own gold
                v = compare_template(tid, g, cand)
                r[f"pos_{v.outcome}"] += 1
                r["pos_decided"] += (v.outcome != "UNRESOLVED")
            except Exception:                         # noqa: BLE001
                r["pos_errors"] += 1
            for inst2, g2 in golds.items():           # negative: another instance
                if inst2 == inst or nspan(g2) == nspan(g):
                    continue
                r["neg_pairs"] += 1
                try:
                    v = compare_template(tid, g2, cand)
                except Exception:                     # noqa: BLE001
                    r["neg_errors"] += 1
                    continue
                if v.outcome == "MATCH":
                    r["cross_instance_accepts"] += 1
                    if len(xi_examples) < 25:
                        xi_examples.append(
                            {"template": tid, "own_gold": nspan(g)[:90],
                             "paired_gold": nspan(g2)[:90],
                             "model": nspan(cand)[-90:]})
        scored[tid] = {"kind": kind, **dict(r)}
        for k, v in r.items():
            tally[k] += v
            per_kind[kind][k] += v
    return {"scored": scored, "skipped": skipped, "unmapped_stems": unmapped,
            "tally": dict(tally),
            "per_kind": {k: dict(v) for k, v in per_kind.items()},
            "cross_instance_examples": xi_examples}


# --------------------------------------------------------------------------

def _report_gg(res):
    t, sc, sk = res["tally"], res["scored"], res["skipped"]
    print("D5.4  GOLD x GOLD")
    print(f"  pairs {t['pairs']} / templates-scored {len(sc)} / "
          f"templates-skipped {len(sk)}   [{len(sc)} x {res['n'] * (res['n'] - 1)}]")
    print(f"  negatives (spans differ)  {t['negatives']}")
    print(f"  MATCH {t['MATCH']}   MISMATCH {t['MISMATCH']}   "
          f"UNRESOLVED {t['UNRESOLVED']}")
    print(f"  ERRORS {t['errors']}  <- counted separately from every verdict (D5.7b)")
    print(f"  FALSE ACCEPTS {t['false_accepts']}    false rejects {t['false_rejects']}")
    print()
    print("  per kind, with the negative count beside the rate (D5.9):")
    print(f"    {'kind':14s} {'pairs':>8s} {'negatives':>10s} {'FA':>5s} "
          f"{'errors':>7s} {'decided':>8s}")
    for k, v in sorted(res["per_kind"].items()):
        dec = v["pairs"] - v["UNRESOLVED"] - v["errors"]
        print(f"    {k:14s} {v['pairs']:8d} {v['negatives']:10d} "
              f"{v['false_accepts']:5d} {v['errors']:7d} "
              f"{100 * dec / v['pairs'] if v['pairs'] else 0:7.1f}%")
    if sk:
        print()
        print(f"  SKIPPED, every one named ({len(sk)}):")
        for tid in sorted(sk):
            print(f"    {tid:46s} {sk[tid][:96]}")
    if res["false_accept_examples"]:
        print()
        print("  false accepts:")
        for e in res["false_accept_examples"][:8]:
            print(f"    {e['template']} [{e['kind']}]")
            print(f"        cand: {e['cand']}")
            print(f"        gold: {e['gold']}")


def _report_ag(res):
    t, sc, sk = res["tally"], res["scored"], res["skipped"]
    print("D5.5  ARCHIVE x GOLD")
    pos = t.get("rows", 0)
    neg = t.get("neg_pairs", 0)
    print(f"  pairs {pos + neg} ({pos} positive + {neg} negative) / "
          f"templates-scored {len(sc)} / templates-skipped {len(sk)}")
    print(f"  archive stems with no template: {len(res['unmapped_stems'])} "
          f"({sum(res['unmapped_stems'].values())} traces) "
          f"{sorted(res['unmapped_stems'])}")
    dec = t.get("pos_decided", 0)
    print(f"  decided rate on the positives  {dec}/{pos} "
          f"({100 * dec / pos if pos else 0:.1f}%)")
    print(f"  cross-instance accepts         {t.get('cross_instance_accepts', 0)}"
          f" / {neg}   <- suspicious, not certainly wrong; audited individually")
    print(f"  ERRORS  positive {t.get('pos_errors', 0)}  negative "
          f"{t.get('neg_errors', 0)}")
    print()
    print("  per kind, with the negative count beside the rate (D5.9):")
    print(f"    {'kind':14s} {'rows':>7s} {'negatives':>10s} {'XI':>4s} {'decided':>8s}")
    for k, v in sorted(res["per_kind"].items()):
        r_ = v.get("rows", 0)
        print(f"    {k:14s} {r_:7d} {v.get('neg_pairs', 0):10d} "
              f"{v.get('cross_instance_accepts', 0):4d} "
              f"{100 * v.get('pos_decided', 0) / r_ if r_ else 0:7.1f}%")
    if sk:
        print()
        print(f"  SKIPPED, every one named ({len(sk)}):")
        for tid in sorted(sk):
            print(f"    {tid:46s} {sk[tid][:96]}")


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--gold-gold", action="store_true")
    ap.add_argument("--archive", action="store_true")
    ap.add_argument("-n", type=int, default=DEFAULT_N)
    ap.add_argument("--json", default="")
    args = ap.parse_args(argv)
    both = not (args.gold_gold or args.archive)

    out, bad = {}, 0
    if both or args.gold_gold:
        gg = gold_gold(args.n)
        out["gold_gold"] = gg
        _report_gg(gg)
        bad += gg["tally"]["false_accepts"] + gg["tally"]["errors"]
        print()
    if both or args.archive:
        ag = archive_gold()
        out["archive_gold"] = ag
        _report_ag(ag)
        bad += ag["tally"].get("pos_errors", 0) + ag["tally"].get("neg_errors", 0)
    if args.json:
        with open(args.json, "w", encoding="utf-8") as fh:
            json.dump(out, fh, indent=1, default=str)
    print()
    print("GATE: zero false accepts and zero errors on gold x gold ->",
          "PASS" if bad == 0 else f"FAIL ({bad})")
    return 1 if bad else 0


if __name__ == "__main__":
    raise SystemExit(main())
