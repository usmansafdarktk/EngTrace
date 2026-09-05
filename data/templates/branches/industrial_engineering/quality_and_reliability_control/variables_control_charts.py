import math
import random
from decimal import Decimal, ROUND_HALF_UP

from pilot.templates.branches.industrial_engineering.constants import (
    CONTROL_CHART_FACTORS,
    SPC_CHARACTERISTICS,
    chart_factor,
)


def _hu(x, places):
    """Half-up rounding via Decimal from the shortest float repr."""
    q = Decimal("1") if places == 0 else Decimal("0." + "0" * (places - 1) + "1")
    v = Decimal(repr(x)).quantize(q, rounding=ROUND_HALF_UP)
    return int(v) if places == 0 else float(v)


# Measurable characteristics for X-bar/R charting with display units and
# sensible decimal places for the sampled grand mean / average range
# (prose anchored to the SPC_CHARACTERISTICS [REALISM] classes).
_T21_SETTINGS = {
    # sf = sampled sigma_frac sub-window (inside the class window; the
    # shaft class is restricted to its lower half because micron-quoted
    # shaft measurements with sigma near 1% of a 10-80 mm target are
    # gauging-inconsistent — R1 cycle 1)
    "shaft diameter (mm)": {"phrase": "the diameter of machined shafts",
                            "unit": "mm", "dp": 3, "sf": (0.001, 0.004)},
    "bottle fill volume (mL)": {"phrase": "the fill volume of bottled juice",
                                "unit": "mL", "dp": 2, "sf": (0.002, 0.015)},
    "coating thickness (micron)": {"phrase": "the thickness of an anodized coating",
                                   "unit": "microns", "dp": 2, "sf": (0.01, 0.05)},
}


# Template 21 (Easy) — Area Q1: Variables Control Charts
def template_xbar_r_control_limits():
    """
    X-bar and R Charts: Trial Control Limits and Process Sigma

    Scenario:
        m preliminary subgroups of size n yield a grand mean x-double-bar
        and an average range R-bar. The trial three-sigma limits use the
        tabulated control-chart constants:

            UCL_x = xbb + A2 * Rbar        LCL_x = xbb - A2 * Rbar
            UCL_R = D4 * Rbar              LCL_R = D3 * Rbar
            sigma-hat = Rbar / d2

        Requested: the X-bar chart's upper control limit (the other
        limits and sigma-hat are required intermediate steps).

    Difficulty: Easy
    Grounding: Montgomery, Introduction to Statistical Quality Control,
        7th ed., Ch. 6 (X-bar and R charts, Sec. 6.2; factors from
        Appendix Table VI — transcribed and derivation-verified in
        constants.py [ON-DISK]). Typology only, never Montgomery's data.
    Physical bounds: subgroup size n in [4, 6] (A2/D3/D4/d2 pulled from
        CONTROL_CHART_FACTORS; D3 = 0 in this range, which the trace
        states explicitly); m in [20, 30] subgroups; characteristic
        class from SPC_CHARACTERISTICS with target window; grand mean
        sampled inside the class target window at the class's display
        precision. COHERENT SPREAD (v2, R1 c1): the true sigma is drawn
        from a per-class sigma-fraction sub-window INSIDE the class's
        SPC_CHARACTERISTICS window (shaft restricted to [0.001, 0.004]
        for metrological realism at micron resolution) and R-bar is
        DERIVED as round(d2*sigma, dp), so the implied sigma-hat over
        the target always lands inside the class window — asserted per
        draw. EXACTNESS: A2*Rbar and
        D4*Rbar are products of a 3-4-dp factor and a fixed-dp Rbar —
        computed in Decimal from the displayed strings and rounded
        half-up at the class dp (the precision the question prescribes
        for the limits); sigma-hat = Rbar/d2 rounded half-up at dp + 1;
        a full-precision solve matches every displayed value by
        construction (true decimal ties handled by Decimal half-up).
        Asserts: LCL_x > 0;
        UCL_R > Rbar; sigma-hat > 0; UCL_x within 10% of the grand
        mean.

    Returns:
        tuple(str, str): (question, solution)
    """
    key = random.choice(sorted(_T21_SETTINGS))
    cfg = _T21_SETTINGS[key]
    lo_t, hi_t = SPC_CHARACTERISTICS[key]["target"]
    dp = cfg["dp"]

    n = random.randint(4, 6)
    m = random.randint(20, 30)
    A2 = chart_factor(n, "A2")
    D3 = chart_factor(n, "D3")
    D4 = chart_factor(n, "D4")
    d2 = chart_factor(n, "d2")

    xbb = round(random.uniform(lo_t, hi_t), dp)
    # COHERENT spread sampling (lesson 53; R1 c1 major): draw the true
    # process sigma from the class sigma-fraction sub-window, then derive
    # R-bar = d2*sigma — so sigma-hat = R-bar/d2 always lands back inside
    # the class window (up to display rounding).
    sf_lo, sf_hi = cfg["sf"]
    sigma_true = xbb * random.uniform(sf_lo, sf_hi)
    rbar = round(d2 * sigma_true, dp)

    # Decimal-exact display chain (lessons 51/65/67): factors and the
    # sampled values are exact decimals.
    dA2, dD4, dD3 = (Decimal(str(A2)), Decimal(str(D4)), Decimal(str(D3)))
    dx, dr = Decimal(f"{xbb:.{dp}f}"), Decimal(f"{rbar:.{dp}f}")
    qq = Decimal(1).scaleb(-dp)          # exactly dp decimal places
    ucl_x = float((dx + dA2 * dr).quantize(qq, rounding=ROUND_HALF_UP))
    lcl_x = float((dx - dA2 * dr).quantize(qq, rounding=ROUND_HALF_UP))
    ucl_r = float((dD4 * dr).quantize(qq, rounding=ROUND_HALF_UP))
    lcl_r = float((dD3 * dr).quantize(qq, rounding=ROUND_HALF_UP))
    qq2 = Decimal(1).scaleb(-(dp + 1))   # sigma-hat carries dp + 1
    sig = float((dr / Decimal(str(d2))).quantize(qq2, rounding=ROUND_HALF_UP))

    implied_sf = (rbar / d2) / xbb
    assert 0.85 * sf_lo <= implied_sf <= 1.15 * sf_hi, \
        f"sigma coherence violated: {implied_sf} vs {cfg['sf']}"
    assert lcl_x > 0, f"LCL_x not positive: {lcl_x}"
    assert ucl_r > rbar and sig > 0, "R-chart/sigma sanity"

    question = (
        f"A quality engineer is setting up X-bar and R charts for "
        f"{cfg['phrase']}. From {m} preliminary subgroups of n = {n} "
        f"measurements each, the grand mean is x-double-bar = "
        f"{xbb:.{dp}f} {cfg['unit']} and the average range is R-bar = "
        f"{rbar:.{dp}f} {cfg['unit']}. Using the standard three-sigma "
        f"control-chart constants for n = {n} (A2 = {A2}, D3 = {D3}, "
        f"D4 = {D4}, d2 = {d2}), determine the upper control limit of "
        f"the X-bar chart, in {cfg['unit']} to {dp} decimals (round "
        f"half up). In your solution, also give both R-chart limits to "
        f"{dp} decimals and the estimated process standard deviation to "
        f"{dp + 1} decimals."
    )

    solution = (
        f"**Given:**\n"
        f"m = {m} subgroups of n = {n}; x-double-bar = {xbb:.{dp}f} "
        f"{cfg['unit']}; R-bar = {rbar:.{dp}f} {cfg['unit']}; "
        f"A2 = {A2}, D3 = {D3}, D4 = {D4}, d2 = {d2}.\n\n"
        f"**Step 1:** R-chart limits — the range chart is checked first "
        f"because the X-bar limits depend on R-bar being in control.\n"
        f"UCL_R = D4 * R-bar = {D4} * {rbar:.{dp}f} = {ucl_r:.{dp}f} "
        f"{cfg['unit']};  LCL_R = D3 * R-bar = {D3} * {rbar:.{dp}f} = "
        f"{lcl_r:.{dp}f} {cfg['unit']} (D3 = 0 for n = {n}, so the R "
        f"chart has no positive lower limit)\n\n"
        f"**Step 2:** X-bar chart limits.\n"
        f"UCL_x = x-double-bar + A2*R-bar = {xbb:.{dp}f} + {A2} * "
        f"{rbar:.{dp}f} = {ucl_x:.{dp}f} {cfg['unit']};  "
        f"LCL_x = {xbb:.{dp}f} - {A2} * {rbar:.{dp}f} = "
        f"{lcl_x:.{dp}f} {cfg['unit']}\n\n"
        f"**Step 3:** Estimate the process standard deviation from the "
        f"average range (requested deliverable).\n"
        f"sigma-hat = R-bar / d2 = {rbar:.{dp}f} / {d2} = "
        f"{sig:.{dp + 1}f} {cfg['unit']}\n\n"
        f"**Answer:** The upper control limit of the X-bar chart is "
        f"{ucl_x:.{dp}f} {cfg['unit']}"
    )

    return question, solution


# Known-sigma X-bar chart classes (lesson 75: display resolution matches
# credible gauging — coating carries dp=1; "machined pins" keeps the shaft
# sf window process-consistent, R1 c1). sf sub-windows as in t21.
_T22_SETTINGS = {
    "shaft diameter (mm)": {"phrase": "the diameter of machined pins",
                            "unit": "mm", "dp": 3, "sf": (0.001, 0.004)},
    "bottle fill volume (mL)": {"phrase": "the fill volume of a canning line",
                                "unit": "mL", "dp": 2, "sf": (0.002, 0.015)},
    "coating thickness (micron)": {"phrase": "the thickness of a paint coat",
                                   "unit": "microns", "dp": 1, "sf": (0.01, 0.05)},
}

# sigma_xbar = sigma/rootn is an exact decimal only when rootn divides a
# power of 10; the required extra digits beyond sigma's dp+1 are:
#   rootn 2 -> +1 digit, rootn 4 -> +2 digits, rootn 5 -> +1 digit.
# (n = 9 is EXCLUDED: /3 is non-terminating — the c1 blocking defect.)
_T22_SX_EXTRA = {2: 1, 4: 2, 5: 1}


# Template 22 (Easy) — Area Q1: Variables Control Charts
def template_xbar_known_sigma_classification():
    """
    X-bar Chart with Known Sigma: Classifying Plotted Subgroup Means

    Scenario:
        A process with KNOWN in-control mean mu0 and standard deviation
        sigma is monitored with an X-bar chart at subgroup size n:

            sigma_xbar = sigma / sqrt(n)
            UCL = mu0 + 3*sigma_xbar,  LCL = mu0 - 3*sigma_xbar

        Eight plotted subgroup means are classified against the limits,
        and the MOST EXTREME mean's standardized value
        z = (xbar - mu0) / sigma_xbar (negative below the center line)
        is requested, to two decimals. The out-of-control count is a
        required intermediate.

    Difficulty: Easy
    Grounding: Montgomery, Introduction to Statistical Quality Control,
        7th ed., Ch. 6 (Sec. 6.2.3, charts based on standard values —
        known mu and sigma). Typology only.
    Physical bounds: n in {4, 16, 25} — sqrt(n) in {2, 4, 5} divides a
        power of 10, so sigma_xbar is displayed EXACTLY (sigma carries
        dp+1 digits; sigma_xbar carries dp+1 plus 1/2/1 extra digits for
        rootn 2/4/5; n = 9 was excluded in cycle 2 because /3 is
        non-terminating — the c1 blocking defect). Limits = mu0 +/-
        3*sigma_xbar are exact at the same precision, so every printed
        equation reproduces from its displayed operands. Characteristic
        class with coherent sigma = sf*mu0 (lesson 74; coating dp = 1
        per lesson 75), with a post-rounding coherence assert
        (implied sf in [0.85*sf_lo, 1.15*sf_hi], as in t21). Plotted
        means at dp+1, placed by standardized offsets: all OOC
        excursions share ONE sampled direction (a single assignable
        cause, R1 c1), in-control |u| <= 2.6, OOC 3.3 <= |u| <= 4.5,
        k = 1..3 among 8; post-display margins re-verified (in-control
        |u| <= 2.75, OOC |u| >= 3.2); extreme unique by >= 0.15. The
        final z is computed in EXACT Fraction arithmetic and screened
        >= 0.002 away from every 2-dp rounding boundary (kills both
        true ties and the near-tie context hazard found by R2 c1), then
        rounded half-up. Asserts: coherence; k in [1, 3];
        |z_extreme| in [3.2, 4.6].

    Returns:
        tuple(str, str): (question, solution)
    """
    from fractions import Fraction

    for _ in range(300):
        key = random.choice(sorted(_T22_SETTINGS))
        cfg = _T22_SETTINGS[key]
        lo_t, hi_t = SPC_CHARACTERISTICS[key]["target"]
        dp = cfg["dp"]
        n = random.choice([4, 16, 25])
        rootn = int(math.isqrt(n))
        sxdp = dp + 1 + _T22_SX_EXTRA[rootn]

        mu0 = round(random.uniform(lo_t, hi_t), dp)
        sf = random.uniform(*cfg["sf"])
        sigma = round(mu0 * sf, dp + 1)
        if sigma <= 0:
            continue
        sf_lo, sf_hi = cfg["sf"]
        if not (0.85 * sf_lo <= sigma / mu0 <= 1.15 * sf_hi):
            continue
        d_mu = Decimal(f"{mu0:.{dp}f}")
        d_sig = Decimal(f"{sigma:.{dp + 1}f}")
        sx = (d_sig / rootn).quantize(Decimal(1).scaleb(-sxdp))
        assert sx * rootn == d_sig            # exact by construction
        ucl = d_mu + 3 * sx
        lcl = d_mu - 3 * sx
        if lcl <= 0:
            continue

        k = random.randint(1, 3)
        direction = random.choice([-1, 1])
        us = []
        for i in range(8):
            if i < k:
                us.append(random.uniform(3.3, 4.5) * direction)
            else:
                us.append(random.uniform(-2.6, 2.6))
        random.shuffle(us)
        means = [float((d_mu + Decimal(repr(u)) * sx)
                       .quantize(Decimal(1).scaleb(-(dp + 1)),
                                 rounding=ROUND_HALF_UP)) for u in us]
        d_means = [Decimal(f"{v:.{dp + 1}f}") for v in means]
        # exact standardized offsets as Fractions of the displayed values
        f_sx = Fraction(str(sx))
        f_mu = Fraction(str(d_mu))
        u_frac = [(Fraction(str(dm)) - f_mu) / f_sx for dm in d_means]
        if any(Fraction(11, 4) < abs(u) < Fraction(16, 5) for u in u_frac):
            continue
        ooc_idx = [i for i, u in enumerate(u_frac) if abs(u) > 3]
        if len(ooc_idx) != k:
            continue
        absu = sorted((abs(u) for u in u_frac), reverse=True)
        if absu[0] - absu[1] < Fraction(3, 20):
            continue
        # 2-dp boundary screen on the exact extreme z: distance of
        # 100*|z| from the nearest half-integer boundary >= 0.2
        zc = absu[0] * 100
        frac_part = zc - int(zc)
        if abs(frac_part - Fraction(1, 2)) < Fraction(1, 5):
            continue
        break
    else:
        raise AssertionError("resample loop exhausted")

    i_ext = max(range(8), key=lambda i: abs(u_frac[i]))
    zf = u_frac[i_ext]
    # half-up (away from zero) rounding of the exact Fraction at 2 dp;
    # the boundary screen guarantees no tie is ever hit
    sign = 1 if zf >= 0 else -1
    cents = abs(zf) * 100
    z_ext = sign * float((int(cents) + (1 if cents - int(cents) >= Fraction(1, 2) else 0))) / 100

    assert 1 <= k <= 3, f"k out of bounds: {k}"
    assert 3.2 <= abs(z_ext) <= 4.6, f"z_ext out of bounds: {z_ext}"

    sxf = float(sx)
    uclf, lclf = float(ucl), float(lcl)
    means_str = ", ".join(f"{v:.{dp + 1}f}" for v in means)

    question = (
        f"A process producing {cfg['phrase']} runs with a known "
        f"in-control mean of mu0 = {mu0:.{dp}f} {cfg['unit']} and a "
        f"known process standard deviation of sigma = "
        f"{sigma:.{dp + 1}f} {cfg['unit']}. Subgroup means of n = {n} "
        f"measurements are plotted on a three-sigma X-bar chart. The "
        f"last eight plotted subgroup means were: {means_str} "
        f"{cfg['unit']}. Compute the control limits, state how many of "
        f"the eight means plot outside them, and report the "
        f"standardized value z = (xbar - mu0)/(sigma/sqrt(n)) of the "
        f"MOST extreme mean (largest |z|; negative if below the center "
        f"line), to two decimals (round half up, away from zero)."
    )

    above = [i for i in ooc_idx if u_frac[i] > 0]
    below = [i for i in ooc_idx if u_frac[i] < 0]
    side_bits = []
    if above:
        side_bits.append(", ".join(f"{means[i]:.{dp + 1}f}" for i in above)
                         + f" {cfg['unit']} above UCL")
    if below:
        side_bits.append(", ".join(f"{means[i]:.{dp + 1}f}" for i in below)
                         + f" {cfg['unit']} below LCL")
    side_str = "; ".join(side_bits)
    verb = "falls" if k == 1 else "fall"

    solution = (
        f"**Given:**\n"
        f"mu0 = {mu0:.{dp}f} {cfg['unit']}; sigma = {sigma:.{dp + 1}f} "
        f"{cfg['unit']}; n = {n}; eight plotted means (listed in the "
        f"question).\n\n"
        f"**Step 1:** Standard error of the subgroup mean. Since "
        f"sqrt({n}) = {rootn},\n"
        f"sigma_xbar = sigma / sqrt(n) = {sigma:.{dp + 1}f} / {rootn} "
        f"= {sxf:.{sxdp}f} {cfg['unit']} (exact)\n\n"
        f"**Step 2:** Three-sigma control limits.\n"
        f"UCL = mu0 + 3*sigma_xbar = {mu0:.{dp}f} + 3 * {sxf:.{sxdp}f} "
        f"= {uclf:.{sxdp}f} {cfg['unit']};  "
        f"LCL = {mu0:.{dp}f} - 3 * {sxf:.{sxdp}f} = {lclf:.{sxdp}f} "
        f"{cfg['unit']}\n\n"
        f"**Step 3:** Classify the plotted means. Comparing each mean "
        f"with the limits, {k} of the eight {verb} outside: {side_str}. "
        f"The remaining {8 - k} lie between LCL and UCL.\n\n"
        f"**Step 4:** Standardize the most extreme mean.\n"
        f"The mean farthest from the center line is "
        f"{means[i_ext]:.{dp + 1}f} {cfg['unit']}:  "
        f"z = ({means[i_ext]:.{dp + 1}f} - {mu0:.{dp}f}) / {sxf:.{sxdp}f} "
        f"= {z_ext:.2f}\n\n"
        f"**Answer:** The standardized value of the most extreme "
        f"plotted mean is {z_ext:.2f}"
    )

    return question, solution


# Chart-selection template settings (same coherent classes; both R-bar
# and s-bar are reported in every instance — derived from ONE class-window
# sigma via d2 and c4, so they are mutually consistent, lesson 74).
_T23_SETTINGS = {
    "shaft diameter (mm)": {"phrase": "the diameter of turned bushings",
                            "unit": "mm", "dp": 3, "sf": (0.001, 0.004)},
    "bottle fill volume (mL)": {"phrase": "the fill volume of a sauce line",
                                "unit": "mL", "dp": 2, "sf": (0.002, 0.015)},
    "coating thickness (micron)": {"phrase": "the thickness of a powder coat",
                                   "unit": "microns", "dp": 1, "sf": (0.01, 0.05)},
}


# Template 23 (Intermediate) — Area Q1: Variables Control Charts
# [BRANCHING: subgroup size drives the X-bar/R vs X-bar/s chart choice]
def template_chart_pair_selection():
    """
    Choosing and Building the Right Variables Chart Pair: X-bar/R vs
    X-bar/s

    Scenario:
        Preliminary data from m subgroups of size n report the grand
        mean AND BOTH spread statistics (average range R-bar and average
        standard deviation s-bar, as SPC software does). The subgroup
        size BRANCHES the correct chart pair (Montgomery's guidance:
        the range loses statistical efficiency for moderate-to-large
        subgroups, so use X-bar/R for small n and X-bar/s once n
        exceeds the 10-to-12 threshold):

            small n:  UCL_x = xbb + A2*Rbar;  R limits D3*Rbar, D4*Rbar;
                      sigma-hat = Rbar/d2
            large n:  UCL_x = xbb + A3*sbar;  s limits B3*sbar, B4*sbar;
                      sigma-hat = sbar/c4

        Requested: the X-bar chart's upper control limit under the
        APPROPRIATE pair (the choice, the spread-chart limits, and
        sigma-hat are required intermediates).

    Difficulty: Intermediate
    Grounding: Montgomery, Introduction to Statistical Quality Control,
        7th ed., Sec. 6.3 (X-bar and s charts preferable when "n is
        moderately large — say, n > 10 or 12"; verbatim guidance
        verified in the on-disk copy — lesson 47 single source). Factors
        from Appendix Table VI via constants.py.
    Physical bounds: branch drawn 50/50 OUTSIDE the resample loop (the
        decisiveness screen below rejects at branch-dependent rates;
        fixing the branch first preserves the 50/50 mix) with n in
        [4, 8] (small) or [13, 20] (large) — the large window starts
        at 13 so every draw clears BOTH readings of Montgomery's
        "n > 10 or 12" guidance (n = 12 satisfies one threshold but
        not the other and is excluded as textbook-ambiguous, R3 c1);
        m in [20, 30]; class target/sigma coherence as in t21/t22
        (sigma from the class sf sub-window; R-bar = round(d2*sigma,
        dp+1) and s-bar = round(c4*sigma, dp+1) BOTH displayed at
        SPC-software precision dp+1 (R1 c1) and derived from the SAME
        sigma, so the two statistics are mutually consistent and
        neither leaks the decision — both are always given).
        DECISIVENESS (R2 c1): draws where the WRONG pair's X-bar UCL
        (the other constants family applied to its own statistic)
        quantizes to the SAME dp+1 value as the correct one are
        rejected and resampled, so the chart choice is load-bearing
        for the final answer in EVERY instance. B3 > 0 for n >= 13
        (no zero-clamp case in the large branch; D3 = 0 for n <= 6 is
        stated when it occurs). EXACTNESS: every displayed limit is a
        Decimal product/sum of exact-decimal operands quantized
        half-up at dp+1 (t21 pattern; the products are exact before
        the single final quantize); sigma-hat is a prescribed half-up
        division at dp+1. Both constants FAMILIES for the drawn n are
        quoted in the question (A2/D3/D4/d2 and A3/B3/B4/c4). Asserts:
        the branch's UCL_x within (xbb, xbb*1.08] (analytic corner:
        3*sigma/sqrt(n) <= 3*sf_hi/sqrt(4) = 7.5% of target at the
        coating ceiling); spread-chart UCL > the given statistic;
        sigma-hat coherence with the class window (0.85-1.15
        tolerance); wrong-pair UCL distinct from the answer.

    Returns:
        tuple(str, str): (question, solution)
    """
    branch = random.choice(["small", "large"])
    for _ in range(300):
        key = random.choice(sorted(_T23_SETTINGS))
        cfg = _T23_SETTINGS[key]
        lo_t, hi_t = SPC_CHARACTERISTICS[key]["target"]
        dp = cfg["dp"]
        m = random.randint(20, 30)
        n = (random.randint(4, 8) if branch == "small"
             else random.randint(13, 20))

        A2 = chart_factor(n, "A2")
        D3 = chart_factor(n, "D3")
        D4 = chart_factor(n, "D4")
        d2 = chart_factor(n, "d2")
        A3 = chart_factor(n, "A3")
        B3 = chart_factor(n, "B3")
        B4 = chart_factor(n, "B4")
        c4 = chart_factor(n, "c4")

        xbb = round(random.uniform(lo_t, hi_t), dp)
        sf_lo, sf_hi = cfg["sf"]
        sigma_true = xbb * random.uniform(sf_lo, sf_hi)
        rbar = round(d2 * sigma_true, dp + 1)
        sbar = round(c4 * sigma_true, dp + 1)

        d_x = Decimal(f"{xbb:.{dp}f}")
        d_r = Decimal(f"{rbar:.{dp + 1}f}")
        d_s = Decimal(f"{sbar:.{dp + 1}f}")
        qq = Decimal(1).scaleb(-(dp + 1))

        # decisiveness screen (R2 c1): compute BOTH candidate X-bar
        # UCLs; the wrong pair's value must differ at the answer
        # precision, else the chart choice is not load-bearing
        d_ucl_r = ((d_x + Decimal(str(A2)) * d_r)
                   .quantize(qq, rounding=ROUND_HALF_UP))
        d_ucl_s = ((d_x + Decimal(str(A3)) * d_s)
                   .quantize(qq, rounding=ROUND_HALF_UP))
        if d_ucl_r == d_ucl_s:
            continue

        if branch == "small":
            ucl_x, ucl_alt = float(d_ucl_r), float(d_ucl_s)
            lcl_x = float((d_x - Decimal(str(A2)) * d_r)
                          .quantize(qq, rounding=ROUND_HALF_UP))
            ucl_sp = float((Decimal(str(D4)) * d_r)
                           .quantize(qq, rounding=ROUND_HALF_UP))
            lcl_sp = float((Decimal(str(D3)) * d_r)
                           .quantize(qq, rounding=ROUND_HALF_UP))
            sig = float((d_r / Decimal(str(d2)))
                        .quantize(qq, rounding=ROUND_HALF_UP))
            spread_stat = rbar
        else:
            ucl_x, ucl_alt = float(d_ucl_s), float(d_ucl_r)
            lcl_x = float((d_x - Decimal(str(A3)) * d_s)
                          .quantize(qq, rounding=ROUND_HALF_UP))
            ucl_sp = float((Decimal(str(B4)) * d_s)
                           .quantize(qq, rounding=ROUND_HALF_UP))
            lcl_sp = float((Decimal(str(B3)) * d_s)
                           .quantize(qq, rounding=ROUND_HALF_UP))
            sig = float((d_s / Decimal(str(c4)))
                        .quantize(qq, rounding=ROUND_HALF_UP))
            spread_stat = sbar
        break
    else:
        raise AssertionError("resample loop exhausted")

    # analytic corner: A2*Rbar = 3*sigma/sqrt(n) <= 3*sf_hi*xbb/sqrt(4)
    # = 7.5% of xbb at the coating ceiling (derived, lesson 56)
    assert xbb < ucl_x <= 1.08 * xbb, f"UCL_x implausible: {ucl_x}"
    assert ucl_sp > spread_stat, "spread UCL sanity"
    assert 0.85 * sf_lo <= sig / xbb <= 1.15 * sf_hi, \
        f"sigma coherence violated: {sig / xbb}"
    assert ucl_x != ucl_alt, "wrong-pair UCL must be distinct"

    sdp = dp + 1
    question = (
        f"A quality engineer is designing control charts for "
        f"{cfg['phrase']}. SPC software summarizing {m} preliminary "
        f"subgroups of n = {n} measurements reports a grand mean of "
        f"x-double-bar = {xbb:.{dp}f} {cfg['unit']}, an average range "
        f"of R-bar = {rbar:.{sdp}f} {cfg['unit']}, and an average "
        f"subgroup standard deviation of s-bar = {sbar:.{dp + 1}f} "
        f"{cfg['unit']}. The tabulated three-sigma constants for "
        f"n = {n} are: A2 = {A2}, D3 = {D3}, D4 = {D4}, d2 = {d2}; "
        f"A3 = {A3}, B3 = {B3}, B4 = {B4}, c4 = {c4}. Select the "
        f"appropriate Shewhart variables chart pair for this subgroup "
        f"size, justify the choice, and determine the X-bar chart's "
        f"upper control limit under that pair, in {cfg['unit']} to "
        f"{sdp} decimals (round half up). In your solution, also give "
        f"the spread chart's limits to {sdp} decimals and the estimated "
        f"process standard deviation to {sdp} decimals."
    )

    if branch == "small":
        choice_step = (
            f"**Step 1:** Select the chart pair. With n = {n}, well "
            f"below the n = 10-to-12 threshold, the range remains a "
            f"statistically efficient spread estimator, so the standard "
            f"choice is the X-bar/R pair (the s chart is reserved for "
            f"larger subgroups, where the range loses efficiency)."
        )
        d3note = (" (D3 = 0 at this n: no positive lower limit)"
                  if D3 == 0 else "")
        spread_step = (
            f"**Step 2:** R-chart limits from the average range.\n"
            f"UCL_R = D4 * R-bar = {D4} * {rbar:.{sdp}f} = "
            f"{ucl_sp:.{sdp}f} {cfg['unit']};  "
            f"LCL_R = D3 * R-bar = {D3} * {rbar:.{sdp}f} = "
            f"{lcl_sp:.{sdp}f} {cfg['unit']}{d3note}"
        )
        x_step = (
            f"**Step 3:** X-bar chart limits from R-bar.\n"
            f"UCL_x = x-double-bar + A2*R-bar = {xbb:.{dp}f} + {A2} * "
            f"{rbar:.{sdp}f} = {ucl_x:.{sdp}f} {cfg['unit']};  "
            f"LCL_x = {xbb:.{dp}f} - {A2} * {rbar:.{sdp}f} = "
            f"{lcl_x:.{sdp}f} {cfg['unit']}"
        )
        sig_step = (
            f"**Step 4:** Estimated process standard deviation.\n"
            f"sigma-hat = R-bar / d2 = {rbar:.{sdp}f} / {d2} = "
            f"{sig:.{sdp}f} {cfg['unit']}"
        )
    else:
        choice_step = (
            f"**Step 1:** Select the chart pair. With n = {n}, above "
            f"the n = 10-to-12 threshold, the range loses statistical "
            f"efficiency as a spread estimator, so the appropriate "
            f"choice is the X-bar/s pair based on the average subgroup "
            f"standard deviation."
        )
        spread_step = (
            f"**Step 2:** s-chart limits from the average standard "
            f"deviation.\n"
            f"UCL_s = B4 * s-bar = {B4} * {sbar:.{dp + 1}f} = "
            f"{ucl_sp:.{sdp}f} {cfg['unit']};  "
            f"LCL_s = B3 * s-bar = {B3} * {sbar:.{dp + 1}f} = "
            f"{lcl_sp:.{sdp}f} {cfg['unit']}"
        )
        x_step = (
            f"**Step 3:** X-bar chart limits from s-bar.\n"
            f"UCL_x = x-double-bar + A3*s-bar = {xbb:.{dp}f} + {A3} * "
            f"{sbar:.{dp + 1}f} = {ucl_x:.{sdp}f} {cfg['unit']};  "
            f"LCL_x = {xbb:.{dp}f} - {A3} * {sbar:.{dp + 1}f} = "
            f"{lcl_x:.{sdp}f} {cfg['unit']}"
        )
        sig_step = (
            f"**Step 4:** Estimated process standard deviation.\n"
            f"sigma-hat = s-bar / c4 = {sbar:.{dp + 1}f} / {c4} = "
            f"{sig:.{sdp}f} {cfg['unit']}"
        )

    solution = (
        f"**Given:**\n"
        f"m = {m} subgroups of n = {n}; x-double-bar = {xbb:.{dp}f} "
        f"{cfg['unit']}; R-bar = {rbar:.{sdp}f} {cfg['unit']}; s-bar = "
        f"{sbar:.{dp + 1}f} {cfg['unit']}; constants for n = {n} as "
        f"listed in the question.\n\n"
        f"{choice_step}\n\n{spread_step}\n\n{x_step}\n\n{sig_step}\n\n"
        f"**Answer:** The X-bar chart's upper control limit under the "
        f"appropriate chart pair is {ucl_x:.{sdp}f} {cfg['unit']}"
    )

    return question, solution


# ARL-design template flavor classes (verbal context only — the
# computation is dimensionless in k and sqrt(n); characteristics are
# non-destructively gauged so large candidate subgroups are credible,
# R1 c1).
_T24_SETTINGS = [
    "the bore diameter of hydraulic valve bodies",
    "the net weight of filled detergent cartons",
    "the flange thickness of stamped steel brackets",
]

# candidate subgroup sizes (integer roots 3/4/5) and the k window in
# 0.01 units: k*sqrt(n) spans [1.74, 3.40] across candidates, so beta
# stays <= 0.8962 and the far-tail |z| = 3 + k*sqrt(n) >= 4.74
_T24_CANDIDATES = [9, 16, 25]
_T24_K_WINDOW = (58, 68)


# Template 24 (Advanced) — Area Q1: Variables Control Charts
def template_arl_beta_mean_shift():
    """
    Designing the Subgroup Size: Beta-Risk and Out-of-Control ARL
    Against a Stated Detection Target

    Scenario:
        An X-bar chart with conventional three-sigma limits will
        monitor a process; quality engineering must choose among
        candidate subgroup sizes n in {9, 16, 25}. The assignable
        cause of concern sustains a mean shift of k process standard
        deviations (sigma unchanged, limits not recomputed). For each
        candidate:

            beta(n) = Phi(3 - k*sqrt(n)) - Phi(-3 - k*sqrt(n))
            ARL1(n) = 1/(1 - beta(n))

        Management states a detection target ARL1 <= T. Requested:
        evaluate beta and ARL1 for EVERY candidate, select the
        SMALLEST subgroup size meeting the target, give the geometric
        probability of detection by the second post-shift subgroup at
        the selected size, and report the selected size's ARL1 to one
        decimal as the final answer. (Construction-earned Advanced,
        lesson 41: three full OC evaluations, a negligibility
        argument, a margined design selection, and geometric
        run-length reasoning — the c1 blind relabel found the
        single-candidate version one notch light.)

    Difficulty: Advanced
    Grounding: Montgomery, Introduction to Statistical Quality Control,
        7th ed., Sec. 6.2.6 (the OC curve of the X-bar chart) and
        Sec. 6.2.7 (average run length), whose
        operating-characteristic function for the
        X-bar chart: beta = Phi(L - k*sqrt(n)) - Phi(-L - k*sqrt(n))
        with L = 3; ARL1 = 1/(1 - beta); larger n increases the
        chart's ability to detect a given shift — verbatim development
        verified in the on-disk copy, lesson 47 single source). Phi
        values are the standard normal CDF rounded to 4 decimals as
        printed in Appendix Table II; the beyond-table convention is
        stated IN-QUESTION (R3 c1).
    Physical bounds: candidates fixed at {9, 16, 25}; k at 2 decimals
        SAMPLED from [0.58, 0.68], but the per-candidate Phi/ARL
        boundary screens permanently reject k in {0.58, 0.59, 0.60,
        0.62}, so the SHIPPING support is {0.61, 0.63..0.68} (7
        values; verified by full enumeration, c2 panel — screens
        narrow the sampled window and the narrowed support is the
        documented truth). Per candidate k*sqrt(n) is an EXACT
        2-decimal value (3k/4k/5k); informative z = 3 - k*sqrt(n)
        spans [0.96, 1.17] / [0.28, 0.56] / [-0.40, -0.05] on the
        shipping support — all inside any student table; far-tail
        argument -(3 + k*sqrt(n)) <= -4.83 with true Phi < 1.1e-6,
        far below the 5e-5 needed to round to 0.0000 at 4 decimals in
        EVERY candidate (no table-truncation ambiguity; convention
        also stated in-question).
        ANSWER-SPACE CEILING (Stage D, previously undisclosed): the
        graded answer is the selected candidate's 1-dp ARL, which the
        7-value shipping k support and 3 candidates confine to 19
        distinct values (modal share ~10%, P(duplicate answer in a
        5-instance set) ~50%). Structural: the Phi/ARL tie screens fix
        the k support, and widening k re-opens the table-rounding
        hazards those screens exist to close. Escalated to the branch
        report. Realized 1-dp ARL bands, RE-MEASURED over 40,000 draws
        after Stage D v2 falsified the previous figures: n=9 takes the 7
        values [5.9, 8.3], n=16 the 7 values [2.6, 3.5], n=25 the 5
        values [1.5, 1.9] (7 x 7 x 5 is not 19 because only the SELECTED
        candidate is graded). The superseded claim of
        [5.9, 9.6] / [2.6, 4.0] / [1.5, 2.2] quoted the k = 0.58 values
        that the paragraph above documents as permanently screened out —
        an upper bound that no shipped instance can reach. The bands are
        pairwise disjoint, and
        the REALIZED within-instance separations are >= 3.3 (9 vs 16)
        and >= 1.1 (16 vs 25), so selecting the wrong candidate always
        yields a decisively different final answer (lesson 77). The
        assert bands below are deliberate slight supersets. Target T on
        the 0.5 grid in [2.0, 11.0], screened: ARL1(25) <= T - 0.3
        (selection always exists) and |ARL1(n) - T| >= 0.3 for every
        candidate (unique, decisively margined selection); all three
        selections reachable. Beta <= 0.8962 removes the c1 latent
        regime where an exact-Phi solver's 1-dp ARL could diverge; in
        addition the ARL screen (below) makes the two paths provably
        agree. EXACTNESS (lessons 51/65/76): each beta equals the
        displayed 4-decimal Phi (far term exactly 0.0000; downward
        instances are reduced to the upward standardization by the
        symmetry sentence in Step 1); each ARL1 is the exact Fraction
        1/(1 - beta) of the DISPLAYED beta with a SINGLE half-up
        rounding at 1 decimal; the two-sample detection probability
        is exact integer arithmetic on the displayed beta (the
        8-decimal square is never printed - H5). Screens
        (reject-and-resample): each candidate's Phi at least
        0.03*10^-4 from its half-up boundary (printed table and erf
        agree); each candidate's 10*ARL1 at least 0.05 from its
        rounding boundary — this covers BOTH the 4-dp path's own tie
        AND a solver carrying unrounded Phi (max divergence
        10*ARL^2*5e-5 <= 0.047 at ARL <= 9.7); the selected
        candidate's beta^2 8-decimal tail at least 20e-8 from its
        4-decimal tie. Asserts: far-tail < 4.9e-5; ARL bands and
        strict ordering ARL(9) > ARL(16) > ARL(25); selection unique
        with the stated margins.

    Returns:
        tuple(str, str): (question, solution)
    """
    from fractions import Fraction

    for _ in range(300):
        phrase = random.choice(_T24_SETTINGS)
        k100 = random.randint(*_T24_K_WINDOW)
        k = k100 / 100
        direction = random.choice(["upward", "downward"])

        cand = []
        ok = True
        for n in _T24_CANDIDATES:
            rootn = int(math.isqrt(n))
            d100 = k100 * rootn            # 100*k*sqrt(n), exact integer
            z_info = (300 - d100) / 100
            z_far = -(300 + d100) / 100
            far_true = 0.5 * (1.0 + math.erf(z_far / math.sqrt(2.0)))
            assert far_true < 4.9e-5, "far tail must vanish at 4 dp"
            phi_true = 0.5 * (1.0 + math.erf(z_info / math.sqrt(2.0)))
            # screen: Phi at least 0.03*10^-4 from the half-up boundary
            scaled = phi_true * 10000.0
            if abs(scaled - math.floor(scaled) - 0.5) < 0.03:
                ok = False
                break
            b = int(Decimal(repr(phi_true)).quantize(Decimal("0.0001"),
                                                     rounding=ROUND_HALF_UP)
                    .scaleb(4))            # beta as a count of 10^-4
            # screen: 10*ARL1 at least 0.05 from its rounding boundary
            # (also forces exact-Phi and table-Phi paths to agree)
            q10 = Fraction(100000, 10000 - b)
            fpart = q10 - math.floor(q10)
            if abs(fpart - Fraction(1, 2)) < Fraction(1, 20):
                ok = False
                break
            arl = (math.floor(q10)
                   + (1 if fpart >= Fraction(1, 2) else 0)) / 10
            cand.append({"n": n, "rootn": rootn, "d100": d100,
                         "z_info": z_info, "z_far": z_far,
                         "b": b, "arl": arl})
        if not ok:
            continue

        arls = [c["arl"] for c in cand]
        t2 = random.randint(4, 22)         # T = t2/2 on [2.0, 11.0]
        T = t2 / 2
        if arls[2] > T - 0.3:              # n = 25 must meet the target
            continue
        if any(abs(a - T) < 0.3 for a in arls):
            continue
        sel = next(c for c in cand if c["arl"] <= T)
        b2 = sel["b"] * sel["b"]           # beta^2 * 10^8, exact
        if abs(b2 % 10000 - 5000) < 20:
            continue
        break
    else:
        raise AssertionError("resample loop exhausted")

    p2_num = 10 ** 8 - b2
    p2_4 = p2_num // 10000 + (1 if p2_num % 10000 >= 5000 else 0)

    assert cand[0]["arl"] > cand[1]["arl"] > cand[2]["arl"], \
        "ARL must strictly decrease in n"
    for c, (lo, hi) in zip(cand, [(5.9, 9.7), (2.5, 4.1), (1.4, 2.3)]):
        assert lo <= c["arl"] <= hi, f"ARL band: n={c['n']}: {c['arl']}"
    assert all(c["arl"] > T for c in cand if c["n"] < sel["n"]), \
        "every smaller candidate must miss the target"
    assert sel["arl"] < T, "selected candidate must meet the target"

    art = "an" if direction == "upward" else "a"
    question = (
        f"Quality engineering is choosing among candidate subgroup "
        f"sizes n = 9, 16, and 25 for an X-bar chart with conventional "
        f"three-sigma limits that will monitor {phrase}. The assignable "
        f"cause of concern sustains {art} {direction} shift of the "
        f"process mean by k = {k:.2f} process standard deviations "
        f"(sigma unchanged; limits not recomputed). Management requires "
        f"that such a shift be detected within an average of at most "
        f"{T:.1f} subgroups after it occurs, i.e. the out-of-control "
        f"average run length ARL1 = 1/(1 - beta) may not exceed "
        f"{T:.1f}. Using standard normal CDF values rounded to 4 "
        f"decimals (for arguments at or beyond -3.99 or 3.99, take Phi "
        f"as 0.0000 or 1.0000), compute beta (to 4 decimals) and ARL1 "
        f"(to 1 decimal, round half up) for each candidate subgroup "
        f"size, select the smallest subgroup size that meets the "
        f"requirement, give the probability (to 4 decimals) that the "
        f"shift is detected no later than the second post-shift "
        f"subgroup at the selected size, and report the ARL1 at the "
        f"selected subgroup size to 1 decimal."
    )

    sym = ("" if direction == "upward" else
           " (by normal symmetry, a downward shift of the same "
           "magnitude is detected with identical probabilities, so the "
           "upward standardization applies)")
    lines1 = ";  ".join(
        f"n = {c['n']}: z = 3 - {k:.2f}*{c['rootn']} = 3 - "
        f"{c['d100'] / 100:.2f} = {c['z_info']:.2f}" for c in cand)
    step1 = (
        f"**Step 1:** Standardized distance from the shifted mean to "
        f"the NEAR control limit for each candidate, z(n) = 3 - "
        f"k*sqrt(n){sym}.\n{lines1}.\n"
        f"The far limit sits at -3 - k*sqrt(n), at or beyond "
        f"{cand[0]['z_far']:.2f} for every candidate, so its Phi term "
        f"is 0.0000 at 4 decimals."
    )
    lines2 = ";  ".join(
        f"beta({c['n']}) = Phi({c['z_info']:.2f}) - 0.0000 = "
        f"{c['b'] / 10000:.4f}" for c in cand)
    step2 = (
        f"**Step 2:** Beta-risk of each candidate — the probability "
        f"the first post-shift subgroup mean still plots inside the "
        f"limits.\n{lines2}"
    )
    lines3 = ";  ".join(
        f"ARL1({c['n']}) = 1/(1 - {c['b'] / 10000:.4f}) = "
        f"1/{(10000 - c['b']) / 10000:.4f} = {c['arl']:.1f}"
        for c in cand)
    step3 = (
        f"**Step 3:** Out-of-control average run length of each "
        f"candidate.\n{lines3}"
    )
    fails = [c for c in cand if c["arl"] > T]
    if fails:
        fail_s = ", ".join(f"ARL1({c['n']}) = {c['arl']:.1f} > {T:.1f}"
                           for c in fails)
        step4 = (
            f"**Step 4:** Select the smallest adequate subgroup size. "
            f"{fail_s}, so "
            f"{'that candidate misses' if len(fails) == 1 else 'those candidates miss'} "
            f"the target; n = {sel['n']} is the smallest candidate "
            f"with ARL1 = {sel['arl']:.1f} <= {T:.1f}."
        )
    else:
        step4 = (
            f"**Step 4:** Select the smallest adequate subgroup size. "
            f"Already ARL1({sel['n']}) = {sel['arl']:.1f} <= {T:.1f}, "
            f"so the smallest candidate n = {sel['n']} meets the "
            f"target."
        )
    step5 = (
        f"**Step 5:** Detection no later than the second post-shift "
        f"subgroup at the selected n = {sel['n']} (geometric: the "
        f"shift survives two independent samples with probability "
        f"beta^2).\n"
        f"P(detected by sample 2) = 1 - beta^2 = 1 - "
        f"{sel['b'] / 10000:.4f}^2 = {p2_4 / 10000:.4f} (4 decimals, "
        f"round half up)"
    )

    solution = (
        f"**Given:**\n"
        f"Candidate subgroup sizes 9, 16, 25; three-sigma X-bar "
        f"limits; sustained {direction} mean shift of k = {k:.2f} "
        f"sigma; detection target: ARL1 at most {T:.1f} subgroups.\n\n"
        f"{step1}\n\n{step2}\n\n{step3}\n\n{step4}\n\n{step5}\n\n"
        f"**Answer:** At the selected subgroup size, the out-of-control "
        f"average run length is {sel['arl']:.1f} subgroups"
    )

    return question, solution
