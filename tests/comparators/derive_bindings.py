"""Generate and validate ``bindings.py``'s tables (D5.8, D5.10).

    python -m tests.comparators.derive_bindings --regenerate   # rewrite the tables
    python -m tests.comparators.derive_bindings                # check them, exit 1 on drift

**Why the tables are generated rather than typed.** D-034: *a number written
beside the rule it describes is not evidence; a number a script recomputes from
the artefact is.* A binding table is exactly the artefact that rule is about --
150 declarations, each of which could be wrong and none of which a green test
suite would question. So it is derived from the corpus, validated against the
corpus, and written out; and the committed copy is checked against a fresh
derivation on every run.

That is not the same as *inferring* a binding at scoring time. The declaration
is fixed at generation time and reviewed; the comparator reads a static table.
A comparator that re-derived the unit from gold's string at scoring time would
agree with gold by construction (D-052).

## D5.10 -- "the gold answer carries a unit", defined

The Phase 5 brief records 68-87 templates as unit-carrying and says the spread
exists because the phrase has no operational definition: *at seed 0*, *on all 12
instances*, or *on any*. Here it is:

    unit_token(t, seed) = the token that ends the clause containing the LAST
                          number of the gold answer span, with trailing
                          punctuation stripped, if it is not in STOPWORDS

    carries_unit[seed0]     unit_token(t, 0) is not None
    carries_unit[any]       not None for at least one of 12 seeds
    carries_unit[always]    not None for ALL 12 seeds          <- ADOPTED
    carries_unit[invariant] always, AND the same token every time

**``always`` is adopted, not ``invariant``.** A unit is a property of the item,
but its *value* legitimately varies on the eight ``dual`` templates, which state
an answer in SI on some instances and imperial on others -- ``basic_stress_strain``
alternates ``MPa`` and ``ksi``, ``axial_deformation`` ranges over
``in``/``kN``/``kips``/``mm``. Requiring invariance would declare those
unit-less, which is the opposite of true. A unit is *declared* only where it is
also invariant, because a comparator can only check a unit it can name; the
varying ones are recorded as unit-carrying and left undeclared, which is a
smaller and honest residual.

Measured, and none of the four readings reproduces the brief's soft "87":

    seed0 115 | any 116 | always 113 | invariant 103      (of 150)

The stop-list is the only judgement in the definition and it is published in
this module rather than described.
"""

from __future__ import annotations

import argparse
import csv
import io
import os
import re
import sys
from collections import Counter

REPO = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", ".."))
if REPO not in sys.path:
    sys.path.insert(0, REPO)

from tests.comparators import bindings as B  # noqa: E402
from tests.comparators.extract import NUM_RE, numbers  # noqa: E402
from tests.comparators.normalize import answer_span, prepare  # noqa: E402
from tests.template_integrity.core import discover, generate  # noqa: E402

INVENTORY = os.path.join(REPO, "docs", "re-implementation-sep", "template_inventory.csv")
BINDINGS_PY = os.path.join(os.path.dirname(__file__), "bindings.py")

#: Instances per template for validation. 2,450 ordered pairs; resolves ~0.12%.
N_VALIDATE = 50
#: Instances used to derive the declaration itself.
N_DERIVE = 12

#: Words that can follow a number in an answer sentence and are not units.
#: Every one was observed in the corpus sweep that produced this list.
STOPWORDS = {
    "and", "or", "the", "a", "an", "is", "are", "was", "were", "to", "of", "at",
    "on", "for", "with", "from", "by", "as", "than", "then", "that", "this",
    "these", "which", "so", "if", "when", "while", "but", "x", "times",
    "approx", "approximately", "about", "over", "under", "less", "greater",
    "more", "each", "total", "value", "values", "point", "points", "step",
    "steps", "no", "yes", "e",
}

#: number, optional closing bracket, then a token that ENDS its clause.
#: Requiring the clause end is what keeps the preposition "in" out while keeping
#: the unit "in" (inches): `0.046 in.` matches, `in the pipe` does not.
_UNIT_TOKEN_RE = re.compile(
    r"\d(?:[\d,]*\.?\d*)?(?:[eE][-+]?\d+)?\s*[>\)\]]?\s*"
    r"([A-Za-zµΩ°][A-Za-z0-9µΩ°/^*_.-]{0,13})\s*(?=[.,;)\]]|$)")

#: **Declared, not derived.**  A label set is not recoverable from gold's
#: structure: `linear` and `not linear` are the same shape, and so are
#: `Underdamped` and `Overdamped`.  These six are declared by hand and then put
#: through exactly the same N>=50 cross-pairing as every derived binding --
#: declaring a binding is not exempting it from validation.
#:
#: The four Phase 4 templates are carried forward from `answer.PHASE4_BINDINGS`
#: rather than re-declared, so there is one definition of them and not two.
def _predeclared():
    from tests.comparators.answer import PHASE4_BINDINGS
    out = {t: dict(b) for t, b in PHASE4_BINDINGS.items()}
    return out


#: **A measurement that changed a decision, recorded because it did.**
#:
#: `damping_classification` and `reynolds_number_flow_regime` were declared here
#: with hand-written label sets -- `Underdamped`/`Overdamped`, `laminar`/
#: `turbulent` -- and both FAILED validation, at 544 and 1,650 false accepts in
#: 2,450 pairs.  The bindings are not wrong; the *pairing* is.  Both items state
#: a quantity AND a class, two instances routinely share the class while
#: differing in the quantity, and a label-only binding correctly says MATCH on a
#: pair whose spans differ.  A whole-span truth predicate reads that as a false
#: accept and cannot do otherwise.
#:
#: So a label-only binding on a mixed answer is **not validatable by this
#: instrument**, and the honest response is to say so rather than to weaken the
#: predicate until the binding passes.  `reynolds_number_flow_regime` binds as
#: `numeric` instead and is flagged `partial`; `damping_classification` stays
#: unbound.  A mixed numeric+categorical `AnswerSpec` is what both actually
#: need, and that is Phase 6's (D6.9).
_LABEL_ONLY_NOT_VALIDATABLE = {
    "template_damping_classification":
        "a label-only binding on a mixed answer is not validatable by whole-span "
        "cross-pairing: it states a damping ratio, a regime and a damped natural "
        "frequency, and two instances routinely share the regime while differing "
        "in both quantities. Measured: 544 false accepts in 2,450 pairs, none of "
        "them a comparator defect. Needs a mixed numeric+categorical AnswerSpec "
        "(D6.9)",
}


_BRACE_RE = re.compile(r"[{\[]\s*[-+*\d][^}\]]*[}\]]")
_FUNC_RE = re.compile(r"\b(?:sinc|cos|sin|exp|log|Q|tanh|sqrt)\s*[\^(]")


def unit_token(sol: str) -> str | None:
    t = prepare(answer_span(sol)[0])
    best = None
    for line in t.split("\n"):
        for m in _UNIT_TOKEN_RE.finditer(line):
            w = m.group(1).rstrip(".").rstrip("*/^_-")
            if w and not w[0].isdigit() and w.lower() not in STOPWORDS:
                best = w
    return best


def spans(ref, n):
    out = []
    for s in range(n):
        i = generate(ref, s, capture=False)
        if not i.ok:
            return None, f"generation error at seed {s}: {i.error}"
        out.append(i.solution)
    return out, None


def derive_units(ref):
    sols, err = spans(ref, N_DERIVE)
    if err:
        return {"error": err}
    toks = [unit_token(s) for s in sols]
    present = [t for t in toks if t]
    distinct = sorted(set(present))
    return {
        "seed0": toks[0] is not None,
        "any": bool(present),
        "always": len(present) == N_DERIVE,
        "invariant": len(present) == N_DERIVE and len(distinct) == 1,
        "distinct": distinct,
        "declared": distinct[0] if (len(present) == N_DERIVE and len(distinct) == 1) else None,
    }


def derive_kind(sols, unit):
    """Propose (kind, options, note) from what gold actually emits."""
    sp = [prepare(answer_span(s)[0]) for s in sols]
    counts = {len(numbers(x)) for x in sp}
    n_q = next(iter(counts)) if len(counts) == 1 else None
    if all(_BRACE_RE.search(x) for x in sp):
        return "sequence", {}, "a braced sequence in every gold span"
    if all(B.vector_components(s) is not None for s in sols):
        return "vector", {}, "an angle-bracket or unit-vector form in every gold span"
    if all(_FUNC_RE.search(x) for x in sp):
        return "symbolic", {"symbols": ("t", "f", "x", "n", "tau")}, \
            "a function call in every gold span"
    if n_q == 1:
        return "numeric", ({"unit": unit} if unit else {}), "exactly one asserted number"
    if n_q and n_q > 1:
        return "multipart", {"n": n_q, **({"unit": unit} if unit else {})}, \
            f"{n_q} asserted numbers in every gold span"
    return None, {}, f"the asserted-number count varies across instances: {sorted(counts)}"


def validate(tid, kind, opts, sols):
    """Cross-pair within the template. Zero false accepts AND zero errors."""
    saved = B.BINDINGS.get(tid)
    B.BINDINGS[tid] = {"kind": kind, "options": opts}
    sp = [re.sub(r"\s+", " ", answer_span(s)[0]).strip() for s in sols]
    r = Counter()
    try:
        for a in range(len(sols)):
            for b in range(len(sols)):
                if a == b:
                    continue
                r["pairs"] += 1
                try:
                    v = B.compare_template(tid, sols[b], sols[a])
                except Exception as exc:                       # noqa: BLE001
                    r["errors"] += 1
                    r["_err"] = r.get("_err") or f"{type(exc).__name__}: {exc}"
                    continue
                same = sp[a] == sp[b]
                if v.outcome == "UNRESOLVED":
                    r["unresolved"] += 1
                elif v.outcome == "MATCH" and not same:
                    r["false_accepts"] += 1
                elif v.outcome != "MATCH" and same:
                    r["false_rejects"] += 1
    finally:
        if saved is None:
            B.BINDINGS.pop(tid, None)
        else:
            B.BINDINGS[tid] = saved
    r["decided"] = r["pairs"] - r["unresolved"] - r["errors"]
    return r


def derive_all():
    inv = {r["template_id"]: r for r in csv.DictReader(open(INVENTORY, encoding="utf-8"))}
    units, binds, unbound, validation = {}, {}, {}, {}
    unit_readings = Counter()
    for ref in discover():
        tid = ref.template_id
        u = derive_units(ref)
        if "error" in u:
            unbound[tid] = u["error"]
            continue
        for k in ("seed0", "any", "always", "invariant"):
            unit_readings[k] += bool(u[k])
        if u["declared"]:
            units[tid] = u["declared"]
        sols, err = spans(ref, N_VALIDATE)
        if err:
            unbound[tid] = err
            continue
        pre = _predeclared().get(tid)
        if pre is not None:
            kind, opts = pre["kind"], pre.get("options", {})
            note = "declared: a label set is not derivable from gold's structure"
        else:
            kind, opts, note = derive_kind(sols[:N_DERIVE], u["declared"])
        if tid in _LABEL_ONLY_NOT_VALIDATABLE:
            unbound[tid] = _LABEL_ONLY_NOT_VALIDATABLE[tid]
            continue
        if kind is None:
            unbound[tid] = note
            continue
        r = validate(tid, kind, opts, sols)
        validation[tid] = {k: int(v) for k, v in r.items() if not k.startswith("_")}
        if r["false_accepts"] or r["errors"]:
            reason = []
            if r["false_accepts"]:
                reason.append(f"{r['false_accepts']} false accepts in {r['pairs']} pairs")
            if r["errors"]:
                reason.append(f"{r['errors']} errors ({r.get('_err', '')[:60]})")
            unbound[tid] = "; ".join(reason)
            continue
        entry = {"kind": kind, "options": opts, "note": note,
                 "answer_type": inv[tid]["answer_type"]}
        # **A binding can validate while checking only part of the answer**, and
        # that is worth recording rather than discovering later.  A `categorical`
        # binding on a gold span that also asserts numbers scores the label and
        # ignores the quantities: `reynolds_number_flow_regime` states a Reynolds
        # number and a regime, and the label binding checks the regime alone.  It
        # is a real binding -- it discriminates, and it never credits a wrong
        # label -- but it is narrower than the answer, so it is flagged.
        # **A binding can validate while checking only part of the answer.**
        # Recorded rather than left to be discovered: an item typed
        # `classification` whose binding is numeric scores the quantity and
        # ignores the class, so a candidate with the right number and the wrong
        # regime is credited.  That is a false accept the binding cannot see, and
        # it is declared here instead of being invisible.
        if kind in ("categorical", "categorical[tuple]"):
            n_num = max(len(numbers(prepare(answer_span(x)[0]))) for x in sols[:N_DERIVE])
            if n_num:
                entry["partial"] = (f"checks the label only; gold also asserts "
                                    f"{n_num} quantit{'y' if n_num == 1 else 'ies'}")
        elif inv[tid]["answer_type"] == "classification":
            entry["partial"] = ("checks the quantity only; gold also states a class "
                                "label, so a right number with a wrong class is "
                                "credited (D6.9)")
        binds[tid] = entry
    return units, binds, unbound, validation, unit_readings


def _fmt(name, d, annot):
    out = [f"{name}: {annot} = {{"]
    for k in sorted(d):
        out.append(f"    {k!r}: {d[k]!r},")
    out.append("}")
    return "\n".join(out)


def write_tables(units, binds, unbound, validation):
    src = io.open(BINDINGS_PY, encoding="utf-8").read()
    nl = "\r\n" if "\r\n" in src else "\n"
    body = "\n\n".join([
        _fmt("DECLARED_UNITS", units, "dict[str, str]"),
        _fmt("BINDINGS", binds, "dict[str, dict[str, Any]]"),
        _fmt("UNBOUND", unbound, "dict[str, str]"),
        _fmt("VALIDATION", validation, "dict[str, dict[str, int]]"),
    ])
    start = src.index("# BEGIN GENERATED")
    end = src.index("# END GENERATED")
    new = (src[:start] + "# BEGIN GENERATED" + nl
           + body.replace("\n", nl) + nl + src[end:])
    io.open(BINDINGS_PY, "wb").write(new.encode("utf-8"))


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--regenerate", action="store_true")
    args = ap.parse_args(argv)

    units, binds, unbound, validation, readings = derive_all()
    print(f"D5.10  carries_unit readings over 150 templates x {N_DERIVE} seeds:")
    for k in ("seed0", "any", "always", "invariant"):
        mark = "  <- ADOPTED" if k == "always" else ""
        print(f"         {k:10s} {readings[k]:3d}{mark}")
    print(f"       units DECLARED (always-present AND invariant): {len(units)}")
    print()
    pairs = sum(v["pairs"] for v in validation.values())
    print(f"D5.8   validated at N={N_VALIDATE} ({N_VALIDATE * (N_VALIDATE - 1)} ordered "
          f"pairs per template; smallest resolvable rate "
          f"~{3 / (N_VALIDATE * (N_VALIDATE - 1)):.3%})")
    print(f"       BOUND   {len(binds):3d} / 150   over {pairs} validated pairs")
    print(f"       UNBOUND {len(unbound):3d} / 150   every one named, with a measured reason")
    inv = {r["template_id"]: r["answer_type"]
           for r in csv.DictReader(open(INVENTORY, encoding="utf-8"))}
    bt, tt = Counter(), Counter()
    for t in inv:
        tt[inv[t]] += 1
        if t in binds:
            bt[inv[t]] += 1
    print()
    print("       by inventory answer_type (reporting only -- never used as a kind):")
    for k in sorted(tt):
        print(f"         {k:16s} {bt[k]:3d} / {tt[k]:3d}")
    print()
    print("       UNBOUND, with the measured reason:")
    for t in sorted(unbound):
        print(f"         {t:46s} [{inv.get(t, '?'):14s}] {unbound[t]}")

    if args.regenerate:
        write_tables(units, binds, unbound, validation)
        print("\nwrote the tables into bindings.py")
        return 0

    drift = []
    if units != B.DECLARED_UNITS:
        drift.append("DECLARED_UNITS")
    if {k: v for k, v in binds.items()} != B.BINDINGS:
        drift.append("BINDINGS")
    if unbound != B.UNBOUND:
        drift.append("UNBOUND")
    if drift:
        print(f"\n!! the committed tables disagree with a fresh derivation: {drift}")
        print("   run --regenerate and review the diff")
        return 1
    print("\ncommitted tables agree with a fresh derivation")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
