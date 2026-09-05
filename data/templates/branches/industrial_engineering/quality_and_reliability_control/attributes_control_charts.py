import math
import random
from decimal import Decimal, ROUND_HALF_UP

from pilot.templates.branches.industrial_engineering.constants import (
    C_CHART_CBAR,
    P_CHART_PBAR,
    P_CHART_SUBGROUP_N,
    SHEWHART_K_SIGMA,
    SPC_NUM_SUBGROUPS,
)


# p-chart contexts (attribute inspection of discrete items). The
# fraction nonconforming is screened into the branch's curated window
# P_CHART_PBAR = (0.01, 0.15); Montgomery's own flagship p-chart
# (Example 7.1, frozen orange juice cans) runs ABOVE that window at
# p-bar = 0.2313, so the curated ceiling is the conservative bound, not
# a typology limit (c5 R1 corrected an earlier comment here that had it
# backwards).
_T27_SETTINGS = [
    {"phrase": "wave-soldered printed circuit boards", "item": "boards"},
    {"phrase": "injection-molded bottle caps", "item": "caps"},
    {"phrase": "machine-stitched garment seams", "item": "seams"},
]

# Sampling plans (m, n), DERIVED from the curated windows rather than
# hand-listed (c6 R1/R2 major: the previous hand-written table claimed to
# be "exactly" the admissible set and silently omitted (25, 50), which is
# the second time a hand enumeration in this file was asserted without
# being run). Three constraints:
#   * m*n divides 10^4, so p-bar = D/(m*n) is EXACT at 4 decimals
#   * m in SPC_NUM_SUBGROUPS, n in P_CHART_SUBGROUP_N
#   * Montgomery's design floor n*p-bar > 3 must leave the plan able to
#     depict the LOWER part of the curated p-bar window. The floor forces
#     p-bar > 3/n, so a plan is kept only if 3/n stays under _T27_FLOOR_
#     FRAC of the window width above its base. This drops n = 50, where
#     the floor alone would force p-bar > 0.06 — the top 40% of the
#     window — so the item could only ever depict poor processes
#     (c6 R1: mean p-bar was 8.9%, "systematically depicts poor
#     processes").
_T27_FLOOR_FRAC = 0.30


def _t27_admissible_plans():
    lo, hi = P_CHART_PBAR
    plans = []
    for m in range(SPC_NUM_SUBGROUPS[0], SPC_NUM_SUBGROUPS[1] + 1):
        for n in range(P_CHART_SUBGROUP_N[0], P_CHART_SUBGROUP_N[1] + 1):
            if 10000 % (m * n):
                continue
            if 3.0 / n > lo + _T27_FLOOR_FRAC * (hi - lo):
                continue
            plans.append((m, n))
    return plans


_T27_PLANS = _t27_admissible_plans()


# Template 27 (Easy) — Area Q3: Attributes Control Charts
# NOT branching: see the docstring's "Branching withdrawn" note.
def template_p_chart_limits_floor():
    """
    p-Chart Trial Limits from Historical Data

    Scenario:
        Historical inspection of m samples of n discrete items each
        found D nonconforming in total. The p chart uses:

            p-bar = D / (m*n)
            se = sqrt(p-bar*(1 - p-bar)/n)
            UCL = p-bar + 3*se,  LCL = max(0, p-bar - 3*se)

        Requested: the upper control limit to 4 decimals. The center
        line and the standard error are required intermediates, and the
        lower limit is reported with Montgomery's max(0, .) convention
        applied — a fraction nonconforming cannot be negative — which
        binds in some draws and not others.

    Difficulty: Easy
    LABEL: one governing principle (three-sigma limits on the binomial
        proportion), applied by direct substitution; every formula and
        convention is supplied in the stem; no equation is selected by
        domain judgment. Blind reviewers returned Easy four times on
        earlier versions of this item. (The one blind Intermediate, at
        c5, was formed on the withdrawn branching design, where the
        graded quantity took two different closing formulas; that item
        no longer exists.)
    BRANCHING WITHDRAWN (c6, 2026-08-10): a Stage D remediation had made
        the max(0, .) floor decision answer-load-bearing by grading the
        distance between the limits in force. Three review cycles then
        chased a leak through three different surface features — the
        regime was readable off n (c3), off D (c4, 100% accurate), and
        finally off D/m (c5). The last one is STRUCTURAL: D/m is
        n*p-bar, whose crossover 9n/(9 + n) tends to 9 — it stays inside
        [7.63, 8.80] over the subgroup sizes THIS template admits, and
        inside [7.35, 8.84] over every (m, n) the exactness constraint
        alone admits — so no plan selection moves it, and the fix that
        would (subgroup sizes of order 10-40) is outside the curated
        P_CHART_SUBGROUP_N window. Measured accuracies are recorded in
        the review log rather than here, so they cannot drift against
        the code (c6 R1: the figure quoted here understated the shipped
        code, because the regimes are no longer balanced 50/50 once the
        branch stopped being sampled for). Rather than force-accept a
        branch a solver can shortcut nine times in ten, the claim is
        withdrawn: the floor still appears in the trace because it is
        part of the procedure, but it does not select the graded answer
        and this template is NOT counted toward the domain's branching
        quota (met by t23, t28 and t24). Recorded as a design finding
        for the branch report, not hidden.
    Grounding: Montgomery, Introduction to Statistical Quality Control,
        7th ed., Sec. 7.2.1 "Development and Operation of the Control
        Chart" (pp. 299-310): trial limits from a historical p-bar, the
        LCL convention verbatim ("sometimes the lower control limit
        LCL < 0. In these cases, we customarily set LCL = 0"), and the
        design floor "lambda = np must exceed 3.00" — all verified in
        the on-disk copy (lesson 47 single source). Typology only; the
        copyright screen below keeps Example 7.1's own data out.
    Physical bounds: all parameter windows come from the branch's
        curated constants.py rather than being inlined (c5 R1/R2 major,
        spec R7): m from SPC_NUM_SUBGROUPS, n from P_CHART_SUBGROUP_N,
        p-bar from P_CHART_PBAR, and the three-sigma multiplier from
        SHEWHART_K_SIGMA. The (m, n) plans are ENUMERATED IN CODE from
        those windows (see _t27_admissible_plans) rather than
        hand-listed, so the table cannot drift from the constants; the
        rule keeps every pair with m*n | 10^4 whose design floor still
        leaves the lower 30% of the p-bar window reachable — which drops
        the n = 50 pairs, where the floor alone would force p-bar > 6%.
        NOTE (c6 R1) that the mean p-bar remains about 8.7%, because D is
        drawn uniformly across a curated window whose CEILING is 15%;
        that is a property of P_CHART_PBAR, not of this sampler, and it
        is left to the curated constant rather than biased here. This
        makes
        p-bar = D/(m*n) EXACT at 4 decimals (stated "(exact)" in the
        trace) without confining the item to poor processes. Every draw also satisfies Montgomery's n*p-bar >= 3
        design floor, enforced strictly as D >= 3m + 1.
        EXACTNESS (lessons 51/65): the ONLY rounding in the chain is the
        standard error's single half-up at 4 decimals — screened at
        least 0.03*10^-4 from its boundary, far outside the ~1e-16 float
        sqrt error; the limits are then exact integer combinations in
        10^-4 units. PATH AGREEMENT (R3 c1 major): a solver who carries
        the EXACT standard error and rounds only the final limits must
        land on the same 4-dp UCL and displayed LCL as the prescribed
        rounded-se chain; draws where the two paths diverge are
        rejected, and the question pins the chain with "using this
        rounded standard error". The displayed raw LCL is screened
        0.0005 clear of zero so the floor decision is unambiguous at
        display precision. The asserts after the loop are REGRESSION
        GUARDS that restate the screens — they cannot fire on a
        reachable draw (c5 R1/R2 correctly called the previous
        docstring's presentation of them as live filters misleading).
        They restate the CURATED WINDOWS, not the attained ranges: the
        p-bar ceiling, both n edges and the design floor are live and
        attained, while the SPC_NUM_SUBGROUPS ceiling of 30 is not
        reachable at all, since no m above 25 divides 10^4 against an
        admissible n (c6 R1 minor — the previous sentence claimed every
        named bound was attained).

    Returns:
        tuple(str, str): (question, solution)
    """
    pb_lo, pb_hi = P_CHART_PBAR
    n_lo, n_hi = P_CHART_SUBGROUP_N
    m_lo, m_hi = SPC_NUM_SUBGROUPS
    k = SHEWHART_K_SIGMA
    for _ in range(400):
        cfg = random.choice(_T27_SETTINGS)
        m, n = random.choice(_T27_PLANS)
        mn = m * n
        mult = 10000 // mn
        assert mult * mn == 10000, "m*n must divide 10^4"
        # p-bar over the whole curated window, subject to Montgomery's
        # design floor. n*p-bar = D/m EXACTLY, and the source says the
        # quantity "must exceed 3.00" (strict), so D >= 3m + 1 (c6 R1
        # hairline: the previous >= 3 attained equality)
        d_lo = max(math.ceil(pb_lo * mn), 3 * m + 1)
        d_hi = math.floor(pb_hi * mn)
        if d_hi < d_lo:
            continue
        d = random.randint(d_lo, d_hi)
        pb4 = d * mult                      # p-bar in 10^-4 units, EXACT
        pbf = pb4 / 10000
        se = math.sqrt(pbf * (1 - pbf) / n)
        scaled = se * 10000.0
        if abs(scaled - math.floor(scaled) - 0.5) < 0.03:
            continue
        se4 = int(Decimal(repr(se)).quantize(Decimal("0.0001"),
                                             rounding=ROUND_HALF_UP)
                  .scaleb(4))
        ucl4 = pb4 + k * se4                # exact integer arithmetic
        lcl_raw = pb4 - k * se4             # may be negative
        # PATH AGREEMENT: a solver carrying the exact standard error and
        # rounding only the final limits must land on the SAME 4-dp
        # values as the prescribed rounded-se chain
        u_full = pb4 + k * scaled
        l_full = pb4 - k * scaled
        if abs(u_full - math.floor(u_full) - 0.5) < 0.03:
            continue
        if round(u_full) != ucl4:
            continue
        lf_round = math.copysign(math.floor(abs(l_full) + 0.5), l_full)
        if abs(abs(l_full) - math.floor(abs(l_full)) - 0.5) < 0.03:
            continue
        if int(lf_round) != lcl_raw:
            continue
        if abs(lcl_raw) < 5:                # floor decision unambiguous
            continue
        # spec §3 copyright rule: Montgomery Example 7.1 is m = 30
        # samples of n = 50 with D = 347 (p-bar = 0.2313), reused in
        # Example 7.2. m = 30 is not an admissible plan here and that
        # p-bar is above the curated ceiling, so the data are already
        # unreachable; screened explicitly so the guarantee survives any
        # later widening of the windows.
        if (m, n, d) == (30, 50, 347):
            continue
        break
    else:
        raise AssertionError("resample loop exhausted")

    lcl4 = max(0, lcl_raw)

    # post-loop regression guards (they restate the screens above and
    # cannot fire on a reachable draw; the bounds are the attained ones)
    assert m_lo <= m <= m_hi, f"m outside SPC_NUM_SUBGROUPS: {m}"
    assert n_lo <= n <= n_hi, f"n outside P_CHART_SUBGROUP_N: {n}"
    assert pb_lo * 10000 <= pb4 <= pb_hi * 10000, f"p-bar window: {pb4}"
    assert n * pb4 >= 30000, f"n*p-bar below Montgomery's floor: {n * pb4}"
    assert pb4 < ucl4, f"UCL must exceed the center line: {ucl4}"
    assert 0 <= lcl4 < pb4, f"LCL must lie in [0, p-bar): {lcl4}"

    pb_s = f"{pb4 / 10000:.4f}"
    se_s = f"{se4 / 10000:.4f}"
    ucl_s = f"{ucl4 / 10000:.4f}"

    question = (
        f"Historical inspection records for {cfg['phrase']} show a "
        f"total of D = {d} nonconforming {cfg['item']} in m = {m} "
        f"samples of n = {n} {cfg['item']} each. Set up the trial "
        f"p chart for the fraction nonconforming: compute the center "
        f"line p-bar = D/(m*n), the standard error "
        f"sqrt(p-bar*(1 - p-bar)/n) to 4 decimals (round half up), "
        f"and then, using this rounded standard error, the three-sigma "
        f"control limits to 4 decimals, recalling that a computed "
        f"lower limit below zero is set to zero. Report the upper "
        f"control limit to 4 decimals."
    )

    step1 = (
        f"**Step 1:** Center line.\n"
        f"p-bar = D/(m*n) = {d}/({m} * {n}) = {d}/{mn} = {pb_s} (exact)"
    )
    step2 = (
        f"**Step 2:** Standard error of the sample fraction.\n"
        f"se = sqrt(p-bar*(1 - p-bar)/n) = sqrt({pb_s} * "
        f"{(10000 - pb4) / 10000:.4f} / {n}) = {se_s} (4 decimals)"
    )
    if lcl_raw < 0:
        step3 = (
            f"**Step 3:** Three-sigma limits.\n"
            f"UCL = p-bar + 3*se = {pb_s} + 3 * {se_s} = {ucl_s};  "
            f"computed LCL = p-bar - 3*se = {pb_s} - 3 * {se_s} = "
            f"{lcl_raw / 10000:.4f}, which is negative — a fraction "
            f"nonconforming cannot fall below zero, so the lower limit "
            f"in force is LCL = max(0, {lcl_raw / 10000:.4f}) = 0.0000"
        )
    else:
        step3 = (
            f"**Step 3:** Three-sigma limits.\n"
            f"UCL = p-bar + 3*se = {pb_s} + 3 * {se_s} = {ucl_s};  "
            f"computed LCL = p-bar - 3*se = {pb_s} - 3 * {se_s} = "
            f"{lcl_raw / 10000:.4f}, which is positive, so the zero "
            f"floor does not bind and this is the lower limit in force"
        )

    solution = (
        f"**Given:**\n"
        f"D = {d} nonconforming in m = {m} samples of n = {n} "
        f"{cfg['item']}; three-sigma p chart.\n\n"
        f"{step1}\n\n{step2}\n\n{step3}\n\n"
        f"**Answer:** The upper control limit of the p chart is {ucl_s}"
    )

    return question, solution


# c-chart contexts: counts of nonconformities per fixed inspection unit
# (Montgomery Ch. 7 typology). Each is a FIXED-size unit with many defect
# opportunities of small probability, which is the Poisson c-chart model;
# the fixed 50-m roll deliberately avoids the variable-area case that
# Sec. 7.3.2 says requires a u chart.
_T28_SETTINGS = [
    {"phrase": "nonconformities found in final inspection of pallets "
               "of 25 dishwasher door assemblies", "unit_desc": "pallet"},
    {"phrase": "solder nonconformities in inspection units of 80 "
               "printed circuit boards", "unit_desc": "inspection unit"},
    {"phrase": "weaving flaws in 50-meter rolls of upholstery fabric",
     "unit_desc": "roll"},
]


# Template 28 (Intermediate) — Area Q3: Attributes Control Charts
# [BRANCHING: the solver's own out-of-control assessment decides whether
#  the limits are revised at all, and so which UCL is adopted]
def template_c_chart_revision():
    """
    c-Chart Trial Limits and the Decision Whether to Revise

    Scenario:
        m inspection units yield T total nonconformities. Three order
        statistics are given — the largest, second largest and smallest
        single-unit count — so every unnamed count lies between the
        smallest and the second largest. Trial limits:

            c-bar = T/m,   UCL = c-bar + 3*sqrt(c-bar),
                           LCL = c-bar - 3*sqrt(c-bar)

        The BRANCH, which is the solver's to resolve:
          * If the largest count plots above the trial UCL, that unit is
            out of control. The second largest is inside, so it is the
            ONLY one — entailed by the data, not asserted. Its cause is
            assignable and confirmed, so it is discarded and the chart
            revised from the remaining m-1 units,
                c-bar_rev = (T - c_max)/(m - 1) to 3 decimals,
                UCL_rev = c-bar_rev + 3*sqrt(c-bar_rev),
            and the REVISED UCL is adopted.
          * If the largest count is inside the trial band, then since no
            other count exceeds it, NO unit plots outside; no removal is
            warranted and the TRIAL UCL is adopted.
        Requested: the UCL actually adopted, to 3 decimals. The trial
        center line, the trial limits and the assessment are all
        load-bearing — a solver who skips them cannot know which limit
        to report.

    Difficulty: Intermediate
    EARNED: t28 was declared Intermediate for four cycles and
        blind-labeled Easy five times, because the stem pre-decided
        every judgment in it. Stage D required de-scaffolding (its §6(b)
        fallback). With the pre-decisions removed, two independent blind
        labels (c5 R3, c6 R3) returned Intermediate, both noting the
        early steps are now load-bearing. Later cycles changed the
        scaffolding no further; they fixed correctness and provenance.
    Grounding: Montgomery, Introduction to Statistical Quality Control,
        7th ed., Sec. 7.3 / 7.3.1 (c chart, three-sigma limits) and
        Sec. 6.2.2 (the preliminary-data workflow: plot the trial
        points, look for assignable causes, discard only confirmed ones,
        revise). Montgomery iterates that loop until every retained
        point plots in control; here ONE pass provably suffices, because
        the second-largest count is screened inside the revised limits
        too. NOTE (c6 R1): Montgomery's own Example 7.3 removes TWO
        points, so single-removal is a deliberate simplification bought
        by construction, not the general case. Typology only.
    Physical bounds: the regime and m are drawn OUTSIDE the resample
        loop (lesson 77 mix preservation) and T, c_max, c_2nd and c_min
        inside it. m is drawn as random.choice([20, 25, 25, 25]) — a 3:1
        weighting toward m = 25, measured 74.9%, disclosed here because
        c7 R2 found it undocumented. m in {20, 25} because those are the
        only values in SPC_NUM_SUBGROUPS that also divide 100, which is
        what makes c-bar = T/m exact at 2 decimals for EVERY integer T
        (proven, not screened — see the dead-screen list below).
        c-bar is confined to [10.0, 25.0]. The ceiling is C_CHART_CBAR's.
        The floor is a HARDCODED DESIGN BOUND, not a curated one: trial-
        LCL positivity needs only c-bar > 9, and 10.0 is a round number
        chosen above it, so the C_CHART_CBAR guard (floor 2.0) would not
        catch a regression in it (c7 R2 caught the previous wording
        claiming c-bar > 9 "forced" 10.0 — it does not).
        WHAT IS ROUNDED: c-bar = T/m is EXACT at 2 decimals. THREE
        quantities are rounded half-up to 3 decimals — both square roots
        and the revised center line (T - c_max)/(m - 1) — and the stem
        pins all three. Each is screened at least 0.03e-3 from its
        half-up boundary, INCLUDING both lower limits, which c7 R1 found
        had a path-agreement check but no boundary screen.
        PATH AGREEMENT is screened from the EXACT values, not from the
        rounded ones: the trial limits against exact sqrt(c-bar), and the
        revised limits against the exact quotient q and exact sqrt(q),
        plus a screen that round3(sqrt(q_exact)) equals the pinned root.
        c7 opened a hole here by building the revised screen from the
        already-rounded center line, which left a solver carrying the
        exact quotient diverging on 22.45% of revise draws (c7 R2 proved
        it over the whole reachable support); measured after the fix,
        0 divergences in 10,014 revise draws.
        THE ORDER STATISTICS ARE JOINTLY REALIZABLE. The m-3 unnamed
        counts each lie in [c_min, c_2nd], so the stated total is
        attainable only if
            (m-3)*c_min <= T - c_max - c_2nd - c_min <= (m-3)*c_2nd,
        enforced exactly as closed bounds on c_min. This condition is
        necessary AND sufficient for integer counts; c7 R1 verified it by
        constructing an explicit witness multiset for all 30,000 draws.
        (c6 shipped without it and 4.16% of draws — including a shipped
        instance — described an inspection record no data set could
        produce.) Both c_2nd and c_min are screened by their order-
        statistic probabilities under the chart's own model, Poisson at
        the mean of the m-1 units other than the largest, over m-1 units:
        each needs P(statistic >= value) >= 0.02 AND P(<= value) >= 0.02.
        The same lambda is used in BOTH regimes, which is the correct
        conditioning and keeps the screen from becoming a regime tell
        (c6 R2 found the earlier per-regime conditioning leaking). Under
        a per-regime-exact model c7 R1 measured 0.30% of in-control draws
        marginally below the 0.02 tail, none below 0.005 — a disclosed
        residual of that deliberate symmetry. The largest count is NOT
        order-statistic screened; in the revise regime it is an
        assignable-cause outlier the in-control model does not govern,
        and in the in-control regime it is empirically plausible anyway
        (c7 R1: worst two-sided p = 0.0864, none below 0.02).
        c_2nd sits at or above the revised center line and c_min strictly
        below it (c7 R1 found these placements undisclosed), and the
        integer-root rejection inside _root3 discards any root landing
        exactly on a whole number at 3 decimals.
        The assessment margin |c_max - trial UCL| >= 0.3 counts makes the
        verdict decisive at display precision — c7 R2 measured the worst
        case at 0.325 counts against a maximum path-induced wobble of
        0.0015, a 217x factor — and c_2nd and c_min are screened inside
        BOTH limit pairs by >= 0.3, which is what makes a single revision
        pass provably sufficient.
        NO CLAIM IS MADE THAT THE REGIME IS HARD TO PREDICT, and the
        earlier claim that "no single stem integer reveals the regime" is
        WITHDRAWN as false. c7 R1 predicted the regime well above chance
        from T alone, and from the gap c_max - c_2nd, neither of which
        reconstructs the tested comparison. The MECHANISM is worth
        recording because it is not obvious: the regime and m are fixed
        outside the resample loop while T is redrawn inside, so the
        regime-dependent acceptance rate prints itself onto T's marginal
        — rejection sampling leaks the branch into the marginals of
        whatever it resamples. (Accuracies are deliberately NOT quoted
        here. Every screen added since has moved them, and this docstring
        has already carried stale leak figures twice; the measurements
        belong in the cycle records, which are dated.) Nor could
        any such claim survive: the boundary IS the smooth curve
        c-bar + 3*sqrt(c-bar), and c5 R1 fitted a two-parameter curve to
        it with zero errors on 30,000 held-out draws. What is true, and
        is the reason this template is counted as branching, is narrower:
        the branch changes the graded answer, and every shortcut above
        still requires computing BOTH candidate limits correctly to be
        worth anything — a T-lookup solver that skips Steps 1-2 has no
        number to report. Figures live in the review log, not here, so
        they cannot drift against the code.
        DEAD SCREENS, retained as regression guards and listed rather
        than presented as live (c7 R2 proved each fires zero times over
        the entire reachable support): the c-bar exactness check, both
        LCL > 0.2 checks, the revised-tightens check, the regime
        cross-check and the assessment-margin check. They are kept
        because a later widening of any window could make them live
        again, which is exactly when a silently-removed guard costs
        something. The Montgomery Example 7.3 center-line screen is dead
        for the same reason and was already disclosed as such.
        ANSWER-SPACE RESIDUAL (restored after the c8 rewrite dropped
        it; Stage D v2 required action). Measured over 30,000 draws:
        152 distinct graded answers, modal share 3.16%, top-10
        concentration 20.3%, and P(two of a random 5-instance pack share
        an answer) about 9.5%. The ceiling is structural — the answer is
        a function of one centre line, and the path-agreement and
        order-statistic screens thin the support hard. It is far wider
        than the 69-value ceiling this template carried before
        de-scaffolding and the 78 it fell to at c6, and comfortable at
        five instances per template, but it must be re-examined before
        any larger instantiation. Escalated to the branch report rather
        than hidden.
        The post-loop asserts are likewise regression guards restating
        the screens; they cannot fire on a reachable draw.

    Returns:
        tuple(str, str): (question, solution)
    """

    def _round3(v):
        return int(Decimal(repr(v)).quantize(Decimal("0.001"),
                                             rounding=ROUND_HALF_UP)
                   .scaleb(3))

    def _boundary_ok(v):
        s = v * 1000.0
        return abs(s - math.floor(s) - 0.5) >= 0.03

    def _root3(v):
        r = math.sqrt(v)
        if not _boundary_ok(r):
            return None
        out = _round3(r)
        return None if out % 1000 == 0 else out

    def _pois_cdf(lam, kk):
        if kk < 0:
            return 0.0
        term = math.exp(-lam)
        tot = term
        for i in range(1, kk + 1):
            term *= lam / i
            tot += term
        return min(tot, 1.0)

    cb_lo, cb_hi = C_CHART_CBAR
    m_lo, m_hi = SPC_NUM_SUBGROUPS
    k = SHEWHART_K_SIGMA
    branch = random.choice(["revise", "in_control"])
    m = random.choice([20, 25, 25, 25])
    for _ in range(4000):
        cfg = random.choice(_T28_SETTINGS)
        T = random.randint(10 * m, 25 * m)
        cb100 = T * 100 // m
        if cb100 * m != T * 100:            # c-bar exact at 2 dp
            continue
        cbar = cb100 / 100
        if not (10.0 <= cbar <= cb_hi):
            continue
        r3 = _root3(cbar)
        if r3 is None:
            continue
        ucl3 = cb100 * 10 + k * r3
        lcl3 = cb100 * 10 - k * r3
        if lcl3 <= 200:
            continue
        u_full = cbar * 1000 + k * math.sqrt(cbar) * 1000
        l_full = cbar * 1000 - k * math.sqrt(cbar) * 1000
        if abs(u_full - math.floor(u_full) - 0.5) < 0.03:
            continue
        if abs(l_full - math.floor(l_full) - 0.5) < 0.03:
            continue                    # LCL boundary (c7 R1)
        if round(u_full) != ucl3 or round(l_full) != lcl3:
            continue

        # c_max: disjoint half-windows whose UNION straddles the trial
        # UCL, so the regime is where the draw lands
        uf = ucl3 / 1000.0
        half = max(3.0, 0.8 * math.sqrt(cbar))
        if branch == "revise":
            lo, hi = math.ceil(uf + 0.3), math.floor(uf + half)
        else:
            lo, hi = math.ceil(uf - half), math.floor(uf - 0.3)
        if hi < lo:
            continue
        c_max = random.randint(lo, hi)
        if (c_max * 1000 > ucl3) != (branch == "revise"):
            continue
        if abs(c_max * 1000 - ucl3) < 300:
            continue

        # revised center line: ROUNDED, with the rounding pinned in the
        # stem. This is what frees T and c_max from the exact-division
        # grid that collapsed the reachable support at c5/c6 — see
        # "WHAT IS ROUNDED" and the answer-space paragraph above.
        q_exact = (T - c_max) / (m - 1)
        if not _boundary_ok(q_exact):
            continue
        q1000 = _round3(q_exact)
        qf = q1000 / 1000
        rr3 = _root3(qf)
        if rr3 is None:
            continue
        uclr3 = q1000 + k * rr3
        lclr3 = q1000 - k * rr3
        if lclr3 <= 200 or uclr3 >= ucl3:
            continue
        # PATH AGREEMENT for the revised chain, screened from the EXACT
        # quotient rather than from the rounded center line. c7 built
        # these from qf, which left a solver who carries
        # (T - c_max)/(m - 1) unrounded to the end diverging on 22.45% of
        # revise draws (c7 R2 proved it over the whole support). The hole
        # was opened by pinning the rounding: while the quotient divided
        # exactly, qf and q_exact were the same number.
        if _round3(math.sqrt(q_exact)) != rr3:
            continue                    # exact q, pinned root (c7 R2)
        ur_full = q_exact * 1000 + k * math.sqrt(q_exact) * 1000
        lr_full = q_exact * 1000 - k * math.sqrt(q_exact) * 1000
        if abs(ur_full - math.floor(ur_full) - 0.5) < 0.03:
            continue
        if abs(lr_full - math.floor(lr_full) - 0.5) < 0.03:
            continue                    # LCL boundary (c7 R1)
        if round(ur_full) != uclr3 or round(lr_full) != lclr3:
            continue
        if cb100 in (1985, 1967) or q1000 in (19850, 19670):
            continue

        # order statistics of the m-1 units other than the largest:
        # Poisson at their own mean, in BOTH regimes (no regime tell)
        lam, kk = qf, m - 1
        hi2 = min(math.floor(min(ucl3, uclr3) / 1000 - 0.3), c_max - 1)
        lo2 = math.ceil(qf)
        if hi2 < lo2:
            continue
        c2_ok = []
        for c2 in range(lo2, hi2 + 1):
            p_le = _pois_cdf(lam, c2) ** kk          # P(max <= c2)
            p_ge = 1.0 - _pois_cdf(lam, c2 - 1) ** kk  # P(max >= c2)
            if p_le >= 0.02 and p_ge >= 0.02:
                c2_ok.append(c2)
        if not c2_ok:
            continue
        c_2nd = random.choice(c2_ok)

        # JOINT FEASIBILITY (c6 R1 blocking): the m-3 unnamed counts lie
        # in [c_min, c_2nd] and must absorb the remaining total exactly
        lo_b = max(math.ceil(lcl3 / 1000 + 0.3),
                   math.ceil(lclr3 / 1000 + 0.3),
                   T - c_max - (m - 2) * c_2nd)
        hi_b = min(int(qf) - 1, c_2nd - 1,
                   (T - c_max - c_2nd) // (m - 2))
        if hi_b < lo_b:
            continue
        cands, wts = [], []
        for cm in range(lo_b, hi_b + 1):
            p_ge = (1.0 - _pois_cdf(lam, cm - 1)) ** kk   # P(min >= cm)
            p_le = 1.0 - (1.0 - _pois_cdf(lam, cm)) ** kk  # P(min <= cm)
            if p_ge >= 0.02 and p_le >= 0.02:
                cands.append(cm)
                wts.append(max(p_ge - (1.0 - _pois_cdf(lam, cm)) ** kk,
                               1e-12))
        if not cands:
            continue
        c_min = random.choices(cands, weights=wts)[0]
        break
    else:
        raise AssertionError("resample loop exhausted")

    ans3 = uclr3 if branch == "revise" else ucl3
    rest = T - c_max - c_2nd - c_min

    # post-loop regression guards (they restate the screens above)
    assert m_lo <= m <= m_hi, f"m outside SPC_NUM_SUBGROUPS: {m}"
    assert cb_lo <= cbar <= cb_hi, f"c-bar outside C_CHART_CBAR: {cbar}"
    assert 10.0 <= cbar, f"c-bar below the LCL-positivity floor: {cbar}"
    assert c_min < c_2nd < c_max, "order statistics out of order"
    assert (m - 3) * c_min <= rest <= (m - 3) * c_2nd, \
        f"order statistics not jointly realizable: {rest}"
    assert (c_max * 1000 > ucl3) == (branch == "revise"), "regime mismatch"
    assert abs(c_max * 1000 - ucl3) >= 300, "assessment not decisive"
    assert c_2nd * 1000 <= min(ucl3, uclr3) - 300, "2nd must be inside both"
    assert c_min * 1000 >= max(lcl3, lclr3) + 300, "min must be inside both"
    assert uclr3 < ucl3, "revised UCL must tighten"
    assert (ans3 == ucl3) == (branch == "in_control"), "branch must matter"

    cb_s = f"{cbar:.2f}"
    r_s = f"{r3 / 1000:.3f}"
    ucl_s = f"{ucl3 / 1000:.3f}"
    lcl_s = f"{lcl3 / 1000:.3f}"
    cbr_s = f"{qf:.3f}"
    rr_s = f"{rr3 / 1000:.3f}"
    uclr_s = f"{uclr3 / 1000:.3f}"
    lclr_s = f"{lclr3 / 1000:.3f}"
    ans_s = f"{ans3 / 1000:.3f}"

    question = (
        f"A c chart is being established for {cfg['phrase']}. "
        f"Inspection of m = {m} consecutive units found T = {T} "
        f"nonconformities in total. Across those units the largest "
        f"count on any single {cfg['unit_desc']} was {c_max}, the "
        f"second largest was {c_2nd}, and the smallest was {c_min}. "
        f"Compute the trial center line c-bar = T/m, take sqrt(c-bar) "
        f"to 3 decimals (round half up), and form the three-sigma "
        f"trial limits from the rounded root. Then determine from "
        f"those counts whether any unit plots outside the trial "
        f"limits. Standard practice applies: a unit that plots outside "
        f"is traced to a confirmed assignable cause and discarded, and "
        f"the chart is revised from the units that remain — take the "
        f"revised center line (T - L)/(m - 1), where L is that largest "
        f"count, to 3 decimals (round "
        f"half up), its square root to 3 decimals (round half up), and "
        f"form the revised limits from those rounded values; if no "
        f"unit plots outside, the trial limits stand as they are. "
        f"Report the upper control limit that is adopted for ongoing "
        f"monitoring, to 3 decimals."
    )

    step1 = (
        f"**Step 1:** Trial center line.\n"
        f"c-bar = T/m = {T}/{m} = {cb_s} (exact)"
    )
    step2 = (
        f"**Step 2:** Trial three-sigma limits (Poisson model: "
        f"variance equals the mean).\n"
        f"sqrt(c-bar) = sqrt({cb_s}) = {r_s} (3 decimals);  "
        f"UCL = {cb_s} + 3 * {r_s} = {ucl_s};  "
        f"LCL = {cb_s} - 3 * {r_s} = {lcl_s}"
    )
    if branch == "revise":
        step3 = (
            f"**Step 3:** Assess the counts against the trial band "
            f"{lcl_s} to {ucl_s}. The largest, {c_max}, is above "
            f"UCL = {ucl_s}, so this {cfg['unit_desc']} is out of "
            f"control. The second largest, {c_2nd}, is below that UCL "
            f"and the smallest, {c_min}, is above LCL = {lcl_s}, and "
            f"every unnamed count lies between {c_min} and {c_2nd} — so "
            f"exactly one {cfg['unit_desc']} plots outside. Its cause is "
            f"assignable and confirmed, so it is discarded and the "
            f"chart revised."
        )
        step4 = (
            f"**Step 4:** Revised center line from the remaining "
            f"{m - 1} units.\n"
            f"c-bar_rev = (T - {c_max})/({m} - 1) = {T - c_max}/{m - 1} "
            f"= {cbr_s} (3 decimals)"
        )
        step5 = (
            f"**Step 5:** Revised three-sigma limits.\n"
            f"sqrt(c-bar_rev) = sqrt({cbr_s}) = {rr_s} (3 decimals);  "
            f"UCL_rev = {cbr_s} + 3 * {rr_s} = {uclr_s};  "
            f"LCL_rev = {cbr_s} - 3 * {rr_s} = {lclr_s}. Every retained "
            f"count lies between {c_min} and {c_2nd}, and both of those "
            f"are inside the revised limits, so no retained unit is out "
            f"of control and no second revision is needed.\n"
            f"adopted UCL = {uclr_s}"
        )
        body = f"{step1}\n\n{step2}\n\n{step3}\n\n{step4}\n\n{step5}"
    else:
        step3 = (
            f"**Step 3:** Assess the counts against the trial band "
            f"{lcl_s} to {ucl_s}. The largest, {c_max}, is below "
            f"UCL = {ucl_s}; the second largest, {c_2nd}, is smaller "
            f"still; and the smallest, {c_min}, is above "
            f"LCL = {lcl_s}. Every unnamed count lies between {c_min} "
            f"and {c_2nd}, so no {cfg['unit_desc']} plots outside the "
            f"trial band."
        )
        step4 = (
            f"**Step 4:** Decide on revision. No point is out of "
            f"control, so there is no assignable-cause unit to discard "
            f"and no basis for revising the center line; the trial "
            f"limits are adopted as they stand.\n"
            f"adopted UCL = {ucl_s}"
        )
        body = f"{step1}\n\n{step2}\n\n{step3}\n\n{step4}"

    solution = (
        f"**Given:**\n"
        f"m = {m} inspection units; T = {T} total nonconformities; "
        f"largest count {c_max}, second largest {c_2nd}, smallest "
        f"{c_min}; three-sigma c chart.\n\n"
        f"{body}\n\n"
        f"**Answer:** The adopted upper control limit of the c chart "
        f"is {ans_s}"
    )

    return question, solution
