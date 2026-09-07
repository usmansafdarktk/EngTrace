"""D4.2 -- derive the normalisation vocabulary from the archive, with counts.

The spec requires the vocabulary to be "derived from archived model outputs,
not invented".  This recomputes every count that ``normalize.py``'s annotations
claim, from ``error_analysis_annotation/samples/`` on disk, so the annotations
cannot drift away from the evidence the way a transcribed reference value can
(D-034).

It also answers **Decision 2**, which the phase brief poses and does not settle:
the archive holds 16 / 24 / 15 / **6** traces for the four Phase 4 templates, and
six traces cannot establish a vocabulary -- by the standing rule (D-024, D-026)
N samples cannot resolve a variant rarer than ~3/N, so six are blind to anything
occurring in under half of outputs.  The brief offers three routes.  **Route 1
is taken**: draw the vocabulary from the whole 2,200 on the argument that these
are *model habits* rather than template-specific behaviour -- and then TEST that
argument rather than assert it, which is what ``habit_test`` below does.

**The test does not support the argument as stated**, and the result is reported
as measured: only 2 of 27 rules are cleanly model habits and 15 are item-driven.
What survives is narrower than route 1 claims and is what the rules actually
need -- the archive fixes that a form EXISTS and what it MEANS, while every
coverage number stays per template at its own resolution limit.  See
``habit_test``'s closing note for what that leaves unresolved.

Run: ``python -m tests.comparators.derive_vocabulary``
"""

from __future__ import annotations

import glob
import json
import re
from collections import Counter, defaultdict

from .normalize import (
    CHAR_MAP,
    HEDGES,
    NEITHER_NOR_RE,
    SUPERSCRIPT_MAP,
    answer_span,
)

ARCHIVE = "error_analysis_annotation/samples/*.jsonl"

PHASE4 = (
    "signal_operations",
    "system_properties_memory_causality",
    "system_property_linearity",
    "incompressible_continuity",
)

#: Every rule the normaliser applies, as a named detector over an answer span.
#: These are what the counts below are counts OF.
PROBES: dict[str, object] = {
    "unicode-minus": re.compile("−"),
    "unicode-superscript": re.compile("[" + "".join(SUPERSCRIPT_MAP) + "]"),
    "unicode-operator": re.compile("[·×≤≥≠→]"),
    "latex-inline-math": re.compile(r"\\\(|\\\)"),
    "latex-display-math": re.compile(r"\\\[|\\\]"),
    "latex-escaped-brace": re.compile(r"\\[{}]"),
    "latex-frac": re.compile(r"\\d?frac"),
    "latex-boxed": re.compile(r"\\boxed"),
    "latex-text": re.compile(r"\\(?:text|mathrm|mathbf)"),
    "latex-exponent-brace": re.compile(r"\^\s*\{"),
    "markdown-bold": re.compile(r"\*\*"),
    "trailing-annotation": re.compile(
        r"[\[(]\s*(?:units?\s*[:=]?\s*)?(?:N/?A|none|unitless|dimensionless"
        r"|system\s+properties)", re.I),
    "marker-final-answer-heading": re.compile(r"##\s*Final\s+Answer", re.I),
    "marker-bold-answer": re.compile(r"\*\*Answer:?\*\*", re.I),
    "marker-the-answer-is": re.compile(r"The\s+answer\s+is", re.I),
    "negator-not": re.compile(r"\bnot\b", re.I),
    "negator-fused-non": re.compile(r"\bnon-?(?:linear|causal)", re.I),
    "negator-neither-nor": NEITHER_NOR_RE,
    "negator-fails": re.compile(r"\bfails?\b", re.I),
    "hedge": re.compile(r"\b(?:" + "|".join(re.escape(h) for h in HEDGES) + r")\b", re.I),
    "brace-sequence": re.compile(r"\{[^{}]*\}"),
    "bracket-sequence": re.compile(r"\[[^\[\]]*,[^\[\]]*\]"),
    "origin-asterisk": re.compile(r"\*\s*-?\d+\s*\*"),
    "explicit-index-list": re.compile(r"for\s+n\s*[=:]", re.I),
    "per-element-assignment": re.compile(r"\b[a-z]\s*\[\s*-?\d+\s*\]\s*="),
    "arbitrary-function": re.compile(r"\b[A-Za-z]\s*\(\s*[a-z]\s*\)"),
    "thousands-separator": re.compile(r"\d,\d{3}\b"),
}


def load_all() -> list[dict]:
    rows = []
    for path in sorted(glob.glob(ARCHIVE)):
        with open(path, encoding="utf-8") as fh:
            for line in fh:
                rows.append(json.loads(line))
    return rows


def stem(row: dict) -> str:
    return row["question_id"].rsplit("__", 1)[0]


#: The answer-marker probes must run on the WHOLE trace, because `answer_span`
#: consumes the marker it matched -- probing the span for a marker counts zero
#: every time.  The first draft did exactly that and reported 0/61 for the
#: markers it also documents as appearing on 44/61 and 47/61, which is a
#: self-contradiction the annotation check below is there to catch
#: (phase4_summary.md S8, error 5).
_WHOLE_TRACE_PROBES = {
    "marker-final-answer-heading",
    "marker-bold-answer",
    "marker-the-answer-is",
}


def probe_counts(rows: list[dict]) -> Counter:
    c = Counter()
    for r in rows:
        whole = r["model_reasoning"]
        span, _ = answer_span(whole)
        for name, rx in PROBES.items():
            if rx.search(whole if name in _WHOLE_TRACE_PROBES else span):
                c[name] += 1
    return c


def report_phase4(rows: list[dict]) -> None:
    print("D4.2 -- vocabulary support, PER TEMPLATE, never averaged")
    print("=" * 96)
    subs = {t: [r for r in rows if stem(r) == t] for t in PHASE4}
    counts = {t: probe_counts(v) for t, v in subs.items()}
    ns = {t: len(v) for t, v in subs.items()}
    print(f"{'rule':30s} " + " ".join(f"{t[:13]:>14s}" for t in PHASE4) + f" {'all 2200':>10s}")
    print(f"{'(n)':30s} " + " ".join(f"{ns[t]:>14d}" for t in PHASE4)
          + f" {len(rows):>10d}")
    print("-" * 96)
    allc = probe_counts(rows)
    for name in PROBES:
        row = " ".join(f"{counts[t][name]:>14d}" for t in PHASE4)
        print(f"{name:30s} {row} {allc[name]:>10d}")
    print()
    print("RESOLUTION LIMIT, stated per template (D-024/D-026: N samples cannot")
    print("resolve a variant rarer than about 3/N):")
    for t in PHASE4:
        n = ns[t]
        print(f"  {t:36s} n={n:3d}  resolves down to ~{3 / n:.0%}")
    print(f"  {'whole archive':36s} n={len(rows):3d}  resolves down to ~{3 / len(rows):.2%}")
    print()
    print("**incompressible_continuity has six traces and resolves nothing rarer")
    print("than one in two.** That is the constraint Decision 2 exists to address.")
    print("The whole-archive column does NOT dissolve it -- see the habit test")
    print("below, which finds the rates item-driven, so a large number in the last")
    print("column is evidence that a form exists and is not evidence about how")
    print("often it occurs in these four templates.")
    print()


def habit_test(rows: list[dict]) -> None:
    """Decision 2's argument, tested rather than asserted.

    Route 1 draws the vocabulary from all 2,200 traces on the argument that
    these surface forms are **model habits** rather than template-specific
    behaviour.  If that is true, a rule's rate should be far more strongly
    determined by *which model* wrote the trace than by *which template* it
    answers.  If it is false -- if templates drive the rates -- then borrowing
    a vocabulary from other templates imports the wrong distribution and the
    route is unsound.

    The statistic: for each rule, the spread of its per-model rates versus the
    spread of its per-template rates, over the whole archive.  A rule whose
    model-spread dominates its template-spread is a habit.
    """
    print("DECISION 2 -- is the vocabulary a model habit or a template property?")
    print("=" * 96)
    by_model: dict[str, list[dict]] = defaultdict(list)
    by_template: dict[str, list[dict]] = defaultdict(list)
    for r in rows:
        by_model[r["model"]].append(r)
        by_template[stem(r)].append(r)
    # Keep only templates with enough traces for a rate to mean anything.
    by_template = {k: v for k, v in by_template.items() if len(v) >= 15}

    mc = {m: probe_counts(v) for m, v in by_model.items()}
    tc = {t: probe_counts(v) for t, v in by_template.items()}

    def spread(counts, groups):
        out = {}
        for name in PROBES:
            rates = [counts[g][name] / len(groups[g]) for g in groups]
            out[name] = max(rates) - min(rates)
        return out

    ms, ts = spread(mc, by_model), spread(tc, by_template)
    print(f"{'rule':30s} {'model spread':>13s} {'template spread':>16s} {'verdict':>22s}")
    print("-" * 96)
    habit = shared = template_driven = 0
    for name in PROBES:
        m, t = ms[name], ts[name]
        if m < 0.05 and t < 0.05:
            verdict, shared = "rare everywhere", shared
        elif m >= 2 * t:
            verdict = "MODEL HABIT"
            habit += 1
        elif t >= 2 * m:
            verdict = "template-driven"
            template_driven += 1
        else:
            verdict = "both"
            shared += 1
        print(f"{name:30s} {m:13.2f} {t:16.2f} {verdict:>22s}")
    print("-" * 96)
    print(f"model habit: {habit}   template-driven: {template_driven}   both/rare: "
          f"{len(PROBES) - habit - template_driven}")
    print()
    print("**THE ARGUMENT ROUTE 1 RESTS ON DOES NOT SURVIVE THIS TEST, AND THE")
    print("RESULT IS REPORTED AS MEASURED RATHER THAN AS EXPECTED.**")
    print()
    print(f"I expected the presentation rules to come out model habits. Only")
    print(f"{habit} of {len(PROBES)} do. {template_driven} are template-driven, "
          f"including every")
    print("Unicode and most LaTeX rule. A model does not write a superscript")
    print("because it is that model; it writes one because the item has an")
    print("exponent in it.")
    print()
    print("What that does and does not invalidate:")
    print()
    print("  It invalidates any COVERAGE claim borrowed across templates. 'The")
    print("  vocabulary handles 95% of observed forms' cannot be argued from the")
    print("  2,200 and applied to the 6, because the rates are item-driven. Every")
    print("  coverage number in D4.2 is therefore stated per template, at the")
    print("  resolution limit printed above.")
    print()
    print("  It does NOT invalidate the RULES. A normalisation rule is a claim")
    print("  about meaning, not about frequency: `\\frac{3x^2}{2}` denotes 1.5x^2")
    print("  whoever writes it and however often. The archive's job for a rule is")
    print("  to show the form exists and to fix what it means, and one occurrence")
    print("  anywhere in 2,200 traces does both. That is the part of route 1 that")
    print("  survives, and it is the part the rules actually use.")
    print()
    print("  A CONFOUND, named because it weakens the test itself: this cannot")
    print("  separate 'the model does not use this form' from 'the item gave it no")
    print("  opportunity'. Superscripts look template-driven partly because")
    print("  incompressible_continuity is the only Phase 4 template with an")
    print("  exponent at all. The better statistic is the rate GIVEN the")
    print("  opportunity, which needs an opportunity detector per rule and is not")
    print("  built here. So the 15 template-driven verdicts are an upper bound on")
    print("  how template-driven the vocabulary really is.")
    print()
    print("RESIDUAL RISK, unresolved and carried: the archive cannot tell me which")
    print("forms are MISSING for incompressible_continuity. Six traces resolve")
    print("nothing rarer than one in two, the forms are item-driven so the other")
    print("2,194 do not fill the gap, and route 2 (generating more) needs D-003.")
    print("For the symbolic comparator specifically the honest position is close")
    print("to the brief's **route 3**: the rules are stated, they are exercised")
    print("against 22 hand-written adversarial cases, and they are NOT validated")
    print("against a representative sample of real output. That is a risk on the")
    print("register, not a solved problem.")
    print()


def check_annotations(rows: list[dict]) -> int:
    """Every entry of ``normalize.SUPPORT``, recomputed from the archive.

    This is the D-034 guard on this phase's own vocabulary: the counts live in
    the module the rules live in, and they are verified against the artefact on
    disk rather than against a number typed beside them.  A disagreement is a
    failure, not a note.
    """
    from .normalize import NO_PHASE4_SUPPORT, SUPPORT

    p4 = [r for r in rows if stem(r) in PHASE4]
    c4, ca = probe_counts(p4), probe_counts(rows)

    print("SUPPORT CHECK -- normalize.SUPPORT recomputed from the archive")
    print("=" * 96)
    bad = 0
    for name, (want4, wantall) in SUPPORT.items():
        got4, gotall = c4[name], ca[name]
        if (got4, gotall) == (want4, wantall):
            continue
        bad += 1
        print(f"  MISMATCH {name:28s} declared ({want4}, {wantall}) "
              f"measured ({got4}, {gotall})")
    missing = sorted(set(PROBES) - set(SUPPORT))
    for name in missing:
        bad += 1
        print(f"  MISSING  {name:28s} probed but not declared in SUPPORT")
    print(f"  {len(SUPPORT) - bad}/{len(SUPPORT)} entries agree with the archive")
    print()
    print(f"Rules with ZERO Phase 4 support ({len(NO_PHASE4_SUPPORT)}), each carried")
    print("on the wider archive alone and each a place the comparator is specified")
    print("against forms these four templates have never been observed to produce:")
    for name in NO_PHASE4_SUPPORT:
        print(f"    {name:28s} 0 in 61, {SUPPORT[name][1]} in 2,200")
    print()
    return bad


def main() -> int:
    rows = load_all()
    report_phase4(rows)
    habit_test(rows)
    return 1 if check_annotations(rows) else 0


if __name__ == "__main__":
    raise SystemExit(main())
