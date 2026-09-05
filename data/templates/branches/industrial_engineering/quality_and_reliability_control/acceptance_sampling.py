import math
import random
from decimal import Decimal, ROUND_HALF_UP

from pilot.templates.branches.industrial_engineering.constants import (
    MIL_STD_105E_CODE_LETTERS_GII,
    MIL_STD_105E_SAMPLE_SIZE,
    MIL_STD_105E_SINGLE_NORMAL_AC,
)


# Usable single-sampling plans: direct (non-arrow) master-table entries
# with Ac <= 5 so the OC binomial sum stays hand-computable (<= 6 terms),
# letters D..L so lots and samples are realistically sized. Derived from
# the transcribed tables at import time — never re-keyed by hand.
_T29_PLANS = []
for _lo, _hi, _lt in MIL_STD_105E_CODE_LETTERS_GII:
    if _lt in "DEFGHJKL":
        for _aql, _ac in MIL_STD_105E_SINGLE_NORMAL_AC[_lt].items():
            if _ac is not None and _ac <= 5:
                _T29_PLANS.append(
                    (_lt, _lo, _hi, _aql, MIL_STD_105E_SAMPLE_SIZE[_lt], _ac))

_T29_PRODUCTS = [
    "molded electrical connectors",
    "stamped mounting brackets",
    "packaged pressure sensors",
]


# Template 29 (Intermediate) — Area Q4: Acceptance Sampling by Attributes
def template_single_sampling_oc_point():
    """
    Single-Sampling Plan from MIL-STD-105E: One Point on the OC Curve

    Scenario:
        Lots of N units are submitted under MIL-STD-105E, general
        inspection level II, normal inspection, at a stated AQL. The
        (quoted) tables assign a sample-size code letter and the
        single-sampling plan (n, Ac, Re). The SUPPLYING PROCESS runs
        at fraction nonconforming p (v2, R1 c1: process framing makes
        the type-B binomial EXACT by definition — Montgomery Sec.
        15.2.2; the c1 finite-lot wording implied a type-A
        hypergeometric situation in which the binomial answer was off
        at the 2nd decimal for the reachable n/N ratios). The number
        of nonconforming items in the sample is then binomial, and a
        lot is accepted when X <= Ac:

            Pa = sum_{x=0}^{Ac} C(n, x) p^x (1-p)^(n-x)

        Requested: the probability of acceptance Pa at the stated p,
        with each binomial term reported to 4 decimals (round half
        up) and Pa as the sum of the rounded terms (stated as the
        single term when Ac = 0). (The plan identification from the
        quoted table facts is a required intermediate.)

    Difficulty: Intermediate
    Grounding: MIL-STD-105E Tables I and II-A via constants.py
        (transcribed subsets, visual page reads logged in
        data_review_log.md); OC-curve computation per Montgomery,
        Introduction to Statistical Quality Control, 7th ed., Ch. 15
        (Sec. 15.2, type-B OC curve, binomial model — verbatim
        development verified in the on-disk copy, lesson 47 single
        source). The needed table rows are QUOTED IN-QUESTION, so the
        problem is self-contained (same convention as quoting chart
        constants).
    Physical bounds: plans drawn from the derived _T29_PLANS list
        (letters D..L, direct master-table entries only, Ac in
        [0, 5], n in [8, 200]); lot size N sampled inside the drawn
        letter's Table-I range; process quality p = k/1000 with k in
        [5, 120] (0.5% to 12.0%, one decimal in percent — an exact
        rational). EXACTNESS (lessons 51/65/76): every binomial term
        C(n,x) p^x (1-p)^(n-x) is computed as an EXACT Fraction;
        each displayed term is a single half-up rounding at 4
        decimals, screened at least 0.03*10^-4 from its boundary; Pa
        is the EXACT integer sum of the rounded terms (a 4-decimal
        multiple, no further rounding). PATH AGREEMENT screen: draws
        where sum-of-rounded-terms differs from the single-rounding
        of the exact sum are rejected, so a solver who carries full
        precision and rounds once at the end gets the IDENTICAL
        4-decimal answer (kills the t26-class path-dependence
        hazard); the exact sum is likewise screened off its 4-dp
        boundary. Pa screened into [0.05, 0.98] (an informative OC
        point). Asserts: n matches Table II-A for the letter; N
        inside the letter's range; term count = Ac + 1; Pa band.

    Returns:
        tuple(str, str): (question, solution)
    """
    from fractions import Fraction

    def _r4(fr):
        s = fr * 10000
        fpart = s - math.floor(s)
        if abs(fpart - Fraction(1, 2)) < Fraction(3, 100):
            return None                      # too close to a boundary
        return math.floor(s) + (1 if fpart >= Fraction(1, 2) else 0)

    for _ in range(300):
        letter, lo, hi, aql, n, ac = random.choice(_T29_PLANS)
        prod = random.choice(_T29_PRODUCTS)
        N = random.randint(lo, hi)
        k = random.randint(5, 120)           # p = k/1000
        p = Fraction(k, 1000)

        terms4 = []
        exact_sum = Fraction(0)
        ok = True
        for x in range(ac + 1):
            t = (math.comb(n, x) * p ** x * (1 - p) ** (n - x))
            r = _r4(t)
            if r is None or r < 1:           # no displayed 0.0000 terms
                ok = False
                break
            terms4.append(r)
            exact_sum += t
        if not ok:
            continue
        r_exact = _r4(exact_sum)
        if r_exact is None:
            continue
        pa4 = sum(terms4)
        if pa4 != r_exact:                   # path agreement
            continue
        if not 500 <= pa4 <= 9800:           # Pa in [0.05, 0.98]
            continue
        break
    else:
        raise AssertionError("resample loop exhausted")

    assert MIL_STD_105E_SAMPLE_SIZE[letter] == n, "n/table mismatch"
    assert lo <= N <= hi, "lot outside code-letter range"
    assert len(terms4) == ac + 1, "term count"
    assert 500 <= pa4 <= 9800, f"Pa out of band: {pa4}"

    p_pct = f"{k / 10:.1f}"
    pa_s = f"{pa4 / 10000:.4f}"
    aql_s = f"{aql:g}"

    if ac == 0:
        eval_clause = (
            f"evaluate the single term P(X = 0) to 4 decimals (round "
            f"half up), which is the acceptance probability"
        )
    else:
        eval_clause = (
            f"evaluate each term P(X = x) for x = 0 to {ac} to 4 "
            f"decimals (round half up), and report the acceptance "
            f"probability as their sum"
        )
    question = (
        f"Lots of N = {N} {prod} are submitted for acceptance sampling "
        f"under MIL-STD-105E, general inspection level II, normal "
        f"inspection, with an AQL of {aql_s}% nonconforming. Table I "
        f"assigns code letter {letter} to lot sizes {lo} to {hi} at "
        f"level II, and the single-sampling master table gives sample "
        f"size n = {n} with acceptance number Ac = {ac} (Re = {ac + 1}) "
        f"for code {letter} at this AQL. The supplying process runs at "
        f"a fraction nonconforming of p = {p_pct}%. Find the "
        f"probability that a submitted lot is accepted (the type-B OC "
        f"point): the count of nonconforming items in the sample is "
        f"binomial with parameters n and p; {eval_clause}."
    )

    term_lines = []
    for x, r in enumerate(terms4):
        term_lines.append(
            f"P(X = {x}) = C({n}, {x}) * {k / 1000:.3f}^{x} * "
            f"{(1000 - k) / 1000:.3f}^{n - x} = {r / 10000:.4f}")
    terms_block = ";\n".join(term_lines)
    sum_block = " + ".join(f"{r / 10000:.4f}" for r in terms4)

    item_clause = ("at most 1 nonconforming item is found" if ac == 1
                   else f"at most {ac} nonconforming items are found")
    step1 = (
        f"**Step 1:** Identify the sampling plan. At general "
        f"inspection level II, a lot of N = {N} falls in the "
        f"{lo}-to-{hi} band, so the code letter is {letter}; the "
        f"master table then gives n = {n}, Ac = {ac}, Re = {ac + 1} "
        f"at AQL {aql_s}%: inspect {n} units, accept the lot if "
        f"{item_clause}."
    )
    step2 = (
        f"**Step 2:** Type-B OC model at the stated process quality. "
        f"With the process running at p = {k / 1000:.3f}, the "
        f"nonconforming count in the sample is X ~ Binomial({n}, "
        f"{k / 1000:.3f}) (exact for process quality regardless of "
        f"lot size), and a lot is accepted when X <= {ac}."
    )
    step3 = (
        f"**Step 3:** Evaluate "
        f"{'the single binomial term' if ac == 0 else 'each binomial term'} "
        f"(4 decimals, round half up).\n{terms_block}"
    )
    if ac == 0:
        step4 = (
            f"**Step 4:** With Ac = 0 the acceptance probability is "
            f"that single term.\n"
            f"Pa = P(X = 0) = {pa_s}"
        )
    else:
        step4 = (
            f"**Step 4:** Acceptance probability is the sum of the "
            f"terms.\n"
            f"Pa = {sum_block} = {pa_s}"
        )

    solution = (
        f"**Given:**\n"
        f"N = {N}; general inspection level II, normal inspection; "
        f"AQL = {aql_s}%; quoted table facts: code letter {letter} "
        f"for lots {lo}-{hi}, plan n = {n}, Ac = {ac}, Re = {ac + 1}; "
        f"process quality p = {p_pct}%.\n\n"
        f"{step1}\n\n{step2}\n\n{step3}\n\n{step4}\n\n"
        f"**Answer:** The probability of accepting the lot is {pa_s}"
    )

    return question, solution


# Rectifying-inspection plans: subset of the t29 list with Ac <= 3 and
# letters D..H (n <= 50, lots <= 500), so the two OC evaluations stay
# hand-computable and the ATI path-agreement screen (below) is cheap
_T30_PLANS = [p for p in _T29_PLANS if p[5] <= 3 and p[0] in "DEFGH"]

_T30_PRODUCTS = [
    "die-cast pump housings",
    "crimped wiring harnesses",
    "anodized handle brackets",
]


# Template 30 (Advanced) — Area Q4: Acceptance Sampling by Attributes
def template_aoq_ati_rectifying():
    """
    Rectifying Inspection: AOQ at Two Quality Levels and the ATI

    Scenario:
        A single-sampling plan (n, Ac) from MIL-STD-105E (table facts
        quoted in-question) operates under RECTIFYING inspection:
        every nonconforming item found is replaced with a conforming
        one — in a rejected lot's 100% screen AND in the sample of a
        lot that is accepted (c4: this blanket-replacement assumption
        is exactly what produces the (N-n)/N factor; c3 R1 blocked
        because the stem had scoped replacement to rejected lots,
        leaving the quoted AOQ underivable). The AOQ and ATI FORMULAS
        are withheld from the stem so the Advanced label is earned by
        derivation, not plug-in (lesson 41; c3 R3 blind-labelled the
        formula-supplied version Intermediate). For process quality p (v2: process framing per the
        t29 type-B fix):

            Pa(p)  = sum_{x=0}^{Ac} C(n,x) p^x (1-p)^(n-x)
            AOQ(p) = Pa * p * (N - n) / N
            ATI(p) = n + (1 - Pa) * (N - n)

        The AOQ curve is NON-MONOTONE in p (it rises, peaks at the
        AOQL, then falls as rectification takes over). Requested:
        evaluate Pa and AOQ at BOTH stated quality levels p1 < p2,
        state which level yields the HIGHER average outgoing quality
        (the counterintuitive comparison is the design point), and
        report the average total inspection ATI at the current
        quality p1, to 1 decimal. (Construction-earned Advanced,
        lesson 41: two full OC evaluations, the AOQ non-monotonicity
        comparison, and the ATI chain.)

    Difficulty: Advanced
    Grounding: MIL-STD-105E Tables I and II-A via constants.py
        (transcribed subsets); rectifying-inspection formulas per
        Montgomery, Introduction to Statistical Quality Control, 7th
        ed., Ch. 15 (Sec. 15.2.4 Rectifying Inspection, AOQ at Eq. 15.4
        and ATI at Eq. 15.6:  AOQ = Pa*p*(N-n)/N,  ATI =
        n + (1-Pa)(N-n) — verbatim development verified in the
        on-disk copy, lesson 47 single source). Table rows QUOTED
        IN-QUESTION (t29 convention).
    Physical bounds: plans from _T30_PLANS (letters D..H, direct
        entries, Ac in [0, 3], n in [8, 50], lots N <= 500 sampled
        inside the letter's Table-I band); p1 = k1/1000 with k1 in
        [8, 60], p2 = k2/1000 with k2 - k1 in [15, 60], k2 <= 120.
        EXACTNESS (lessons 51/65/76/80): each Pa is the exact-integer
        sum of half-up-rounded 4-decimal binomial terms (every term
        AND the exact sum boundary-screened at 0.03*10^-4, and the
        two solver paths — round-each-term-then-sum vs
        sum-then-round-once — forced to agree, t29 pattern); AOQ is
        the exact rational pa4*k*(N-n)/(10*N) in ppm, half-up to a
        whole ppm, boundary-screened (>= 0.02 from the half), AND
        (v2, all three c1 reviewers) a solver carrying the EXACT Pa
        into AOQ must land on the SAME whole ppm — both AOQ paths
        screened per level; ATI*10 = 10*n + (10^4 - pa4)*(N-n)/10^3
        is an exact rational, half-up at 1 decimal, boundary-screened
        (>= 0.05), AND a solver carrying the EXACT Pa into the ATI
        formula must round to the same 1-decimal value (NOTE: the
        exact-Pa drift, up to 0.5e-4*(N-n)*10 ~ 0.25 in 10*ATI
        units, can EXCEED the 0.05 boundary margin, so the code never
        relies on the margin — it screens the exact-Pa path's
        rounding directly for equality; v2 corrects the c1
        docstring's self-contradictory bound claim).
        DECISIVENESS (lesson 77): |AOQ(p1) - AOQ(p2)| >= max(150 ppm,
        5% of the larger), so the comparison verdict is decisive, and
        BOTH verdict directions are reachable (the p2 window
        straddles typical AOQL locations). VERDICT PROSE (c3, R1/R3
        c2 blocking): each branch's explanation claims ONLY facts
        entailed by the screened inequality — the c2 p2-branch
        asserted 'still on its rising limb', which is FALSE whenever
        the AOQL peak lies between p1 and p2 (~32% of p2-verdict
        draws); the c3 text states the defectives-per-accepted-lot
        mechanism without locating p2 relative to the peak. The
        p1-branch text ('beyond its peak, rectification dominates')
        was verified true in every reachable p1-verdict draw by two
        independent c2 reviewers. AOQ additionally SCREENED
        to at most 60000 ppm (an outgoing quality worse than 6% is
        rejected as an implausible operating point; the analytic
        ceiling is AOQ < p2 <= 120000 ppm — derived, lesson 44, so
        the cap is a screen, never a bare assert). Asserts: n/Ac
        match the quoted table row; N in the letter band; AOQ in
        (0, 60000] ppm; ATI in (n, N); verdict margin.

    Returns:
        tuple(str, str): (question, solution)
    """
    from fractions import Fraction

    def _r4(fr):
        s = fr * 10000
        fpart = s - math.floor(s)
        if abs(fpart - Fraction(1, 2)) < Fraction(3, 100):
            return None
        return math.floor(s) + (1 if fpart >= Fraction(1, 2) else 0)

    def _pa4(n, ac, k):
        p = Fraction(k, 1000)
        terms4 = []
        exact = Fraction(0)
        for x in range(ac + 1):
            t = math.comb(n, x) * p ** x * (1 - p) ** (n - x)
            r = _r4(t)
            if r is None or r < 1:
                return None
            terms4.append(r)
            exact += t
        r_ex = _r4(exact)
        if r_ex is None or sum(terms4) != r_ex:
            return None
        return terms4, sum(terms4), exact

    for _ in range(500):
        letter, lo, hi, aql, n, ac = random.choice(_T30_PLANS)
        prod = random.choice(_T30_PRODUCTS)
        N = random.randint(max(lo, n * 2), hi)
        k1 = random.randint(8, 60)
        k2 = k1 + random.randint(15, 60)
        if k2 > 120:
            continue

        res1 = _pa4(n, ac, k1)
        res2 = _pa4(n, ac, k2)
        if res1 is None or res2 is None:
            continue
        terms1, pa1, exact1 = res1
        terms2, pa2, exact2 = res2
        if not (500 <= pa1 <= 9800 and 300 <= pa2 <= 9700):
            continue

        # AOQ in ppm: exact rational, half-up to whole ppm, screened
        def _aoq(pa4v, kv):
            a = Fraction(pa4v * kv * (N - n), 10 * N)
            fpart = a - math.floor(a)
            if abs(fpart - Fraction(1, 2)) < Fraction(1, 50):
                return None
            return math.floor(a) + (1 if fpart >= Fraction(1, 2) else 0)

        aoq1 = _aoq(pa1, k1)
        aoq2 = _aoq(pa2, k2)
        if aoq1 is None or aoq2 is None:
            continue

        # exact-Pa AOQ path agreement (v2; lesson 80): a solver
        # carrying the exact Pa must land on the same whole ppm
        def _aoqe(pa_exact, kv):
            a = pa_exact * Fraction(1000 * kv * (N - n), N)
            fpart = a - math.floor(a)
            if abs(fpart - Fraction(1, 2)) < Fraction(1, 50):
                return None
            return math.floor(a) + (1 if fpart >= Fraction(1, 2) else 0)

        if _aoqe(exact1, k1) != aoq1 or _aoqe(exact2, k2) != aoq2:
            continue
        # plausibility screen (NOT an assert; the analytic ceiling is
        # AOQ < p2 <= 120000 ppm — derived, lesson 44): outgoing
        # quality worse than 6% is rejected as an operating point
        if aoq1 > 60000 or aoq2 > 60000:
            continue
        if abs(aoq1 - aoq2) < max(150, 0.05 * max(aoq1, aoq2)):
            continue
        # the p1-verdict prose states that most lots at p2 are
        # rejected and screened: require it (c3 R1 minor)
        if aoq1 > aoq2 and pa2 >= 5000:
            continue

        # ATI at p1: exact rational *10, half-up 1 dp, screened; the
        # exact-Pa path must round identically
        ati10f = Fraction(10 * n) + Fraction((10000 - pa1) * (N - n), 1000)
        fpart = ati10f - math.floor(ati10f)
        if abs(fpart - Fraction(1, 2)) < Fraction(1, 20):
            continue
        ati10 = math.floor(ati10f) + (1 if fpart >= Fraction(1, 2) else 0)
        ati10_ex = 10 * n + (1 - exact1) * (N - n) * 10
        fe = ati10_ex - math.floor(ati10_ex)
        if abs(fe - Fraction(1, 2)) < Fraction(1, 50):
            continue
        if math.floor(ati10_ex) + (1 if fe >= Fraction(1, 2) else 0) != ati10:
            continue
        break
    else:
        raise AssertionError("resample loop exhausted")

    assert MIL_STD_105E_SAMPLE_SIZE[letter] == n, "n/table mismatch"
    assert lo <= N <= hi, "lot outside code-letter band"
    assert 0 < aoq1 <= 60000 and 0 < aoq2 <= 60000, "AOQ band"
    assert 10 * n < ati10 < 10 * N, "ATI must lie between n and N"
    assert abs(aoq1 - aoq2) >= 150, "AOQ verdict margin"

    p1_s = f"{k1 / 1000:.3f}"
    p2_s = f"{k2 / 1000:.3f}"
    q1_s = f"{(1000 - k1) / 1000:.3f}"
    q2_s = f"{(1000 - k2) / 1000:.3f}"
    pa1_s = f"{pa1 / 10000:.4f}"
    pa2_s = f"{pa2 / 10000:.4f}"
    ati_s = f"{ati10 / 10:.1f}"
    aql_s = f"{aql:g}"
    higher = "p1" if aoq1 > aoq2 else "p2"

    if ac == 0:
        acc_ask = ("by modelling the count of nonconforming items in "
                   "the sample as binomial with parameters n and p and "
                   "evaluating the single term P(X = 0) to 4 decimals "
                   "(round half up)")
    else:
        acc_ask = (f"by modelling the count of nonconforming items in "
                   f"the sample as binomial with parameters n and p, "
                   f"evaluating each term P(X = x) for x = 0 to {ac} to "
                   f"4 decimals (round half up), and summing them")
    question = (
        f"Lots of N = {N} {prod} are bought under MIL-STD-105E, "
        f"general inspection level II, normal inspection, AQL = "
        f"{aql_s}% (Table I: code letter {letter} for lot sizes {lo} "
        f"to {hi}; master table: n = {n}, Ac = {ac}, Re = {ac + 1}), "
        f"operated as RECTIFYING inspection: every rejected lot is "
        f"100% inspected, and every nonconforming item found — "
        f"whether in a rejected lot's full screen or in the sample "
        f"drawn from a lot that is accepted — is replaced with a "
        f"conforming one. The supplier's process currently runs at "
        f"p1 = {k1 / 10:.1f}% nonconforming and could deteriorate to "
        f"p2 = {k2 / 10:.1f}% (assume the plan stays on normal "
        f"inspection for this analysis). For each process quality "
        f"level, first find the probability that a lot is accepted "
        f"{acc_ask}"
        f". Then determine, at each level, the average outgoing "
        f"quality — the long-run fraction nonconforming in lots "
        f"leaving this station once rectification is accounted for "
        f"— to the nearest whole part per million (round half up), "
        f"and state which process quality level leaves the WORSE "
        f"(higher) average outgoing quality. Finally, at the current "
        f"quality, report the average total inspection per lot — the "
        f"long-run mean number of units inspected per lot submitted "
        f"— to 1 decimal (round half up)."
    )

    def _term_block(n_, k_, q_, terms):
        lines = []
        for x, r in enumerate(terms):
            lines.append(
                f"P(X = {x}) = C({n_}, {x}) * {k_ / 1000:.3f}^{x} * "
                f"{q_}^{n_ - x} = {r / 10000:.4f}")
        return ";\n".join(lines)

    acc_rule = (
        "accept only if no nonconforming item is found" if ac == 0
        else "accept if at most 1 nonconforming item is found"
        if ac == 1
        else f"accept if at most {ac} nonconforming items are found")
    step1 = (
        f"**Step 1:** The quoted plan: code {letter}, n = {n}, "
        f"Ac = {ac}, Re = {ac + 1} — inspect {n} units per lot, "
        f"{acc_rule}; rejected lots are "
        f"fully screened and made conforming (rectification)."
    )
    def _pa_line(label, terms, pa_s_):
        if ac == 0:
            return f"Pa({label}) = P(X = 0) = {pa_s_}"
        return (f"Pa({label}) = "
                + " + ".join(f"{r / 10000:.4f}" for r in terms)
                + f" = {pa_s_}")

    step2 = (
        f"**Step 2:** Acceptance probability at the current process "
        f"quality p1 = {p1_s} (X1 ~ Binomial({n}, {p1_s})).\n"
        f"{_term_block(n, k1, q1_s, terms1)};\n"
        f"{_pa_line('p1', terms1, pa1_s)}"
    )
    step3 = (
        f"**Step 3:** Acceptance probability at the deteriorated "
        f"process quality p2 = {p2_s} (X2 ~ Binomial({n}, {p2_s})).\n"
        f"{_term_block(n, k2, q2_s, terms2)};\n"
        f"{_pa_line('p2', terms2, pa2_s)}"
    )
    if higher == "p1":
        verdict_txt = (
            f"AOQ(p1) is the larger: outgoing quality is WORSE at the "
            f"BETTER incoming level — the AOQ curve is not monotone in "
            f"process quality; beyond its peak, rectification "
            f"dominates (most lots at p2 are rejected and screened "
            f"clean), so the deteriorated process actually leaves a "
            f"cleaner outgoing stream."
        )
    else:
        verdict_txt = (
            f"AOQ(p2) is the larger: although fewer lots pass at p2 "
            f"(Pa falls from {pa1_s} to {pa2_s}), the lots that do "
            f"pass carry a higher fraction nonconforming, and that "
            f"effect outweighs the extra screening of rejected lots — "
            f"the deteriorated process leaves the worse outgoing "
            f"quality."
        )
    step4 = (
        f"**Step 4:** Average outgoing quality. Under the stated "
        f"policy a rejected lot is screened 100% and leaves with no "
        f"nonconforming units at all, while an accepted lot leaves "
        f"with them only among the {N} - {n} = {N - n} units never "
        f"inspected (those found in the sample were replaced). A "
        f"fraction Pa of lots therefore carries about p*({N} - {n}) "
        f"nonconforming units in {N}, giving AOQ = Pa*p*(N - n)/N.\n"
        f"AOQ(p1) = {pa1_s} * {p1_s} * ({N} - {n})/{N} * 10^6 = "
        f"{aoq1} ppm;  "
        f"AOQ(p2) = {pa2_s} * {p2_s} * ({N} - {n})/{N} * 10^6 = "
        f"{aoq2} ppm (each to the nearest whole ppm).\n"
        f"{verdict_txt}"
    )
    step5 = (
        f"**Step 5:** Average total inspection at the current "
        f"quality. Every submitted lot costs the {n} sampled units; "
        f"the fraction 1 - Pa that is rejected costs the other "
        f"{N - n} as well, so ATI = n + (1 - Pa)*(N - n).\n"
        f"ATI = n + (1 - Pa(p1))*(N - n) = {n} + (1 - {pa1_s}) * "
        f"({N} - {n}) = {ati_s} units per lot"
    )

    solution = (
        f"**Given:**\n"
        f"N = {N}; plan n = {n}, Ac = {ac} (code {letter}, AQL "
        f"{aql_s}%); rectifying inspection; process quality currently "
        f"p1 = {k1 / 10:.1f}%, deteriorated case p2 = {k2 / 10:.1f}%."
        f"\n\n"
        f"{step1}\n\n{step2}\n\n{step3}\n\n{step4}\n\n{step5}\n\n"
        f"**Answer:** The average total inspection per lot at the "
        f"current quality is {ati_s} units"
    )

    return question, solution
