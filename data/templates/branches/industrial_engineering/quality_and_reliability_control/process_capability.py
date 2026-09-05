import math
import random
from decimal import Decimal, ROUND_HALF_UP

from pilot.templates.branches.industrial_engineering.constants import (
    COATING_METROLOGY_FLOOR,
    HARD_ANODIZE_THICKNESS,
    RESISTOR_SERIES_BY_TOLERANCE,
    SPC_CHARACTERISTICS,
    chart_factor,
)


# Process-capability classes (coherent sigma sub-windows and display
# precisions follow the Q1 conventions: lessons 74/75; sigma-hat is a
# division result and carries dp+2 — prescribed in-question).
_T25_SETTINGS = {
    "shaft diameter (mm)": {"phrase": "the diameter of precision-ground spindles",
                            "unit": "mm", "dp": 3, "sf": (0.001, 0.004)},
    "bottle fill volume (mL)": {"phrase": "the fill volume of a beverage line",
                                "unit": "mL", "dp": 2, "sf": (0.002, 0.015)},
    "coating thickness (micron)": {"phrase": "the thickness of a galvanized zinc coat",
                                   "unit": "microns", "dp": 1, "sf": (0.01, 0.05)},
}


# Template 25 (Easy) — Area Q2: Process Capability
def template_cp_cpk_from_specs():
    """
    Cp and Cpk from Two-Sided Specification Limits with sigma-hat =
    R-bar/d2

    Scenario:
        A stable process (X-bar/R charts in control) with grand mean
        x-double-bar and average range R-bar is compared against
        two-sided specification limits LSL/USL:

            sigma-hat = R-bar / d2
            Cp  = (USL - LSL) / (6*sigma-hat)
            Cpu = (USL - x-double-bar) / (3*sigma-hat)
            Cpl = (x-double-bar - LSL) / (3*sigma-hat)
            Cpk = min(Cpu, Cpl)

        Requested: Cpk to 2 decimals (sigma-hat, Cp, both one-sided
        ratios, and the verdict against the stated 1.33 benchmark are
        required intermediates).

    Difficulty: Easy
    Grounding: Montgomery, Introduction to Statistical Quality Control,
        7th ed., Ch. 8 (Sec. 8.3 Process Capability Ratios; 8.3.1 Cp, 8.3.2 the off-centre process, eq. 8.9, process capability ratios Cp and Cpk;
        one-sided ratios and the min form verified in the on-disk copy
        — lesson 47 single source). d2 from Appendix Table VI via
        constants.py. The 1.33 benchmark is stated IN-QUESTION (no
        hidden assumption).
    Physical bounds: n in [4, 6], m in [20, 30]; class target/sigma
        coherence as in t21 (sigma from the class sf sub-window; R-bar
        = round(d2*sigma, dp+1)); sigma-hat prescribed at dp+2 (a
        division result at software precision — granularity relative
        to sigma stays below ~0.5%). CONSTRUCTIVE SPECS: the spec
        half-distances are built FROM target one-sided ratios sampled
        in [0.75, 2.05] (LSL = xbb - round(3*sigma-hat*Cpl_t, dp),
        USL = xbb + round(3*sigma-hat*Cpu_t, dp)), so specs land at
        gauge resolution dp and the realized ratios stay in a derived
        band (rounding shifts each ratio by at most
        0.5*10^-dp/(3*sigma-hat), largest in the coating class ~0.09).
        EXACTNESS (lessons 51/65/76): sigma-hat is a single half-up
        Decimal division of displayed R-bar by d2; Cp/Cpu/Cpl are exact
        Fractions of DISPLAYED operands (spec limits, x-double-bar,
        3*sigma-hat) rounded half-up ONCE at 2 decimals. Screens
        (reject-and-resample): each of Cp/Cpu/Cpl at least 0.002 from
        its 2-dp rounding boundary; |Cpu - Cpl| >= 0.03 so the min is
        DECISIVE at display precision (lesson 77 grading-model note:
        exact-match); |Cpk - 1.33| >= 0.03 so the verdict is decisive;
        LSL > 0. Asserts: Cp in (0.6, 2.3); Cpk in (0.55, 2.2);
        Cpk <= Cp; verdict margin and min margin as screened.

    Returns:
        tuple(str, str): (question, solution)
    """
    from fractions import Fraction

    for _ in range(300):
        key = random.choice(sorted(_T25_SETTINGS))
        cfg = _T25_SETTINGS[key]
        lo_t, hi_t = SPC_CHARACTERISTICS[key]["target"]
        dp = cfg["dp"]
        n = random.randint(4, 6)
        m = random.randint(20, 30)
        d2 = chart_factor(n, "d2")

        xbb = round(random.uniform(lo_t, hi_t), dp)
        sf_lo, sf_hi = cfg["sf"]
        sigma_true = xbb * random.uniform(sf_lo, sf_hi)
        rbar = round(d2 * sigma_true, dp + 1)
        if rbar <= 0:
            continue

        d_r = Decimal(f"{rbar:.{dp + 1}f}")
        sig = float((d_r / Decimal(str(d2)))
                    .quantize(Decimal(1).scaleb(-(dp + 2)),
                              rounding=ROUND_HALF_UP))
        if sig <= 0:
            continue

        cpl_t = random.uniform(0.75, 2.05)
        cpu_t = random.uniform(0.75, 2.05)
        off_l = round(3 * sig * cpl_t, dp)
        off_u = round(3 * sig * cpu_t, dp)
        lsl = round(xbb - off_l, dp)
        usl = round(xbb + off_u, dp)
        if lsl <= 0 or usl <= lsl:
            continue

        # exact ratios from DISPLAYED values only
        f_x = Fraction(f"{xbb:.{dp}f}")
        f_l = Fraction(f"{lsl:.{dp}f}")
        f_u = Fraction(f"{usl:.{dp}f}")
        f_3s = 3 * Fraction(f"{sig:.{dp + 2}f}")
        if f_3s == 0:
            continue
        cpu_f = (f_u - f_x) / f_3s
        cpl_f = (f_x - f_l) / f_3s
        cp_f = (f_u - f_l) / (2 * f_3s)

        # decisive min and decisive verdict at 2 dp
        if abs(cpu_f - cpl_f) < Fraction(3, 100):
            continue
        cpk_f = min(cpu_f, cpl_f)
        if abs(cpk_f - Fraction(133, 100)) < Fraction(3, 100):
            continue
        # 2-dp boundary screens (>= 0.002 from every half-up boundary)
        ok = True
        for v in (cp_f, cpu_f, cpl_f):
            cents = v * 100
            fpart = cents - math.floor(cents)
            if abs(fpart - Fraction(1, 2)) < Fraction(1, 5):
                ok = False
                break
        if not ok:
            continue
        break
    else:
        raise AssertionError("resample loop exhausted")

    def _r2(fr):
        cents = fr * 100
        return (math.floor(cents)
                + (1 if cents - math.floor(cents) >= Fraction(1, 2) else 0)) / 100

    cp = _r2(cp_f)
    cpu = _r2(cpu_f)
    cpl = _r2(cpl_f)
    cpk = _r2(cpk_f)
    capable = cpk_f >= Fraction(133, 100)
    side = "upper" if cpu_f < cpl_f else "lower"

    assert 0.6 < cp < 2.3, f"Cp out of band: {cp}"
    assert 0.55 < cpk < 2.2, f"Cpk out of band: {cpk}"
    assert cpk <= cp + 1e-9, "Cpk must not exceed Cp"
    assert abs(cpk - 1.33) >= 0.025, "verdict margin"

    sdp = dp + 2
    question = (
        f"An X-bar/R chart study of {cfg['phrase']} shows the process "
        f"in statistical control, with a grand mean of x-double-bar = "
        f"{xbb:.{dp}f} {cfg['unit']} and an average range of R-bar = "
        f"{rbar:.{dp + 1}f} {cfg['unit']} from {m} subgroups of n = {n} "
        f"measurements (d2 = {d2} for n = {n}). The specification "
        f"limits are LSL = {lsl:.{dp}f} {cfg['unit']} and USL = "
        f"{usl:.{dp}f} {cfg['unit']}. Estimate the process standard "
        f"deviation as sigma-hat = R-bar/d2 to {sdp} decimals, then "
        f"compute Cp, both one-sided ratios Cpu and Cpl, and the "
        f"process capability index Cpk, each to 2 decimals (round half "
        f"up), and state whether the process meets the commonly used "
        f"capability benchmark of 1.33."
    )

    verdict = (
        f"Since Cpk = {cpk:.2f} {'>=' if capable else '<'} 1.33, the "
        f"process {'meets' if capable else 'does not meet'} the stated "
        f"capability benchmark"
        + ("" if capable else
           f" (the {side} specification side is the limiting one)")
        + "."
    )

    solution = (
        f"**Given:**\n"
        f"x-double-bar = {xbb:.{dp}f} {cfg['unit']}; R-bar = "
        f"{rbar:.{dp + 1}f} {cfg['unit']}; n = {n}, d2 = {d2}; LSL = "
        f"{lsl:.{dp}f} {cfg['unit']}; USL = {usl:.{dp}f} {cfg['unit']}; "
        f"benchmark 1.33.\n\n"
        f"**Step 1:** Estimate the process standard deviation from the "
        f"average range.\n"
        f"sigma-hat = R-bar / d2 = {rbar:.{dp + 1}f} / {d2} = "
        f"{sig:.{sdp}f} {cfg['unit']}\n\n"
        f"**Step 2:** Potential capability (centering ignored).\n"
        f"Cp = (USL - LSL) / (6*sigma-hat) = ({usl:.{dp}f} - "
        f"{lsl:.{dp}f}) / (6 * {sig:.{sdp}f}) = {cp:.2f}\n\n"
        f"**Step 3:** One-sided capability ratios.\n"
        f"Cpu = (USL - x-double-bar) / (3*sigma-hat) = ({usl:.{dp}f} - "
        f"{xbb:.{dp}f}) / (3 * {sig:.{sdp}f}) = {cpu:.2f};  "
        f"Cpl = (x-double-bar - LSL) / (3*sigma-hat) = ({xbb:.{dp}f} - "
        f"{lsl:.{dp}f}) / (3 * {sig:.{sdp}f}) = {cpl:.2f}\n\n"
        f"**Step 4:** Actual capability is the smaller one-sided "
        f"ratio.\n"
        f"Cpk = min(Cpu, Cpl) = min({cpu:.2f}, {cpl:.2f}) = {cpk:.2f}. "
        f"{verdict}\n\n"
        f"**Answer:** The process capability index Cpk is {cpk:.2f}"
    )

    return question, solution


# Replacement for the discarded template_ppm_nonconforming_spec (cycle
# cap, see review log). Area Q2, Intermediate. No normal-table
# dependency anywhere (lesson 82).
#
# Standards-anchored classes (lesson 81):
#   anodize  — MIL-A-8625 Type III hardcoat; mu AND the whole spec band
#              inside 20-110 microns, band >= 10 microns and sigma >=
#              1.5 microns so the data stay above coating-metrology
#              repeatability (ASTM B244/B487) — c2 R1;
#   resistor — catalogue nominal +/- tolerance with EXACT bands: 5%
#              pairs with E24 nominals divisible by 20, 10% with E12
#              nominals (IEC 60063 assigns E24 to 5%, E12 to 10%), so
#              nominal*tol/100 is always a whole ohm and the quoted
#              "nominal +/- tol, i.e. LSL/USL" identity is literally
#              true (c2 R1+R2 blocking: rounded bands contradicted the
#              stated tolerance in ~15% of resistor draws and the
#              alternate reading changed the answer every time).
_T26R_ANODIZE = {
    "phrase": "the thickness of a hard-anodized (Type III) layer on "
              "hydraulic cylinder bores",
    "unit": "microns", "dp": 1, "sf": (0.01, 0.05), "tsub": (48, 88),
}
# IEC 60063 series and their tolerance pairing now live in constants.py
# (spec R7, Stage D v2 required action 3), with the provenance caveat
# recorded there.
_T26R_E24_5PCT = RESISTOR_SERIES_BY_TOLERANCE[5]
_T26R_E12_10PCT = RESISTOR_SERIES_BY_TOLERANCE[10]
_T26R_RESISTOR = {
    "phrase": "the resistance of thick-film power resistors",
    "unit": "ohms", "dp": 0, "sf": (0.005, 0.03),
}


# Template 26R (Intermediate) — Area Q2: Process Capability
def template_sigma_reduction_for_cpk():
    """
    Reducing Process Spread to Reach a Target Cpk

    Scenario:
        A stable process with known mean mu and standard deviation
        sigma sits inside a two-sided specification, off-centre and
        not capable. Crucially the process is not even POTENTIALLY
        capable — Cp = (USL - LSL)/(6*sigma) is itself below the
        target — so re-centring alone could not reach it and the
        spread must come down. The solver computes the current
        ratios, then INVERTS the capability definition (with the mean
        fixed, Cpk = d_min/(3*sigma), so the target holds exactly
        when sigma <= d_min/(3*Cpk*)) and expresses the shortfall as
        a percentage reduction in sigma.

        Requested: the required percentage reduction in sigma to 1
        decimal (Cp, Cpu, Cpl, Cpk, d_min and sigma_max are required
        intermediates).

    Difficulty: Intermediate
    EARNED (lesson 41; c2 R3 blind-labelled the c2 draft Easy because
        the stem handed over sigma_max = d_min/(3*Cpk*) and the
        percentage formula ready-made): the question now states only
        the capability DEFINITIONS and the target. The solver must
        see that a fixed mean makes Cpk inversely proportional to
        sigma, invert it, and convert to a percentage. No result
        formula is supplied.
    Grounding: Montgomery, Introduction to Statistical Quality Control,
        7th ed., Ch. 8 (Sec. 8.3 Process Capability Ratios; 8.3.1 Cp, 8.3.2 the off-centre process, eq. 8.9: Cp, the one-sided ratios and the min
        form for an off-centre process — verbatim development
        verified in the on-disk copy, lesson 47 single source).
        Distribution-free ratio algebra throughout: no normality claim
        appears in the stem (c2 R3: the c2 "normally distributed"
        opener was never used by the chain).
    Physical bounds: PREMISE SCREEN (c2 R1 major — 32.6% of c2 draws
        had Cp >= Cpk*, i.e. re-centring alone would have met the
        target, making the whole exercise indefensible): Cp is
        computed exactly and required to sit at least 0.03 BELOW the
        target in EXACT Fraction arithmetic, so variance reduction is genuinely
        unavoidable; the solution states this explicitly. Since
        Cp = Cpk/(1 - f) for a symmetric band with off-centre
        fraction f, and Cp = Cpk*(1 + r)/2 for a constructed band
        with distance ratio r, the drift/ratio windows are derived
        per draw from the target rather than sampled blindly.
        ANODIZE — mu in [48, 88] microns at dp 1, nested well inside the
        MIL-A-8625F Type III producible envelope (12.7-114 microns);
        the realized mean averages ~70 microns, i.e. above the 50.8
        micron default but well inside the envelope (c4 R1 corrected an
        earlier "straddles the default" claim). Band screened inside
        [20, 110] microns
        (round 5-micron drawing callouts were tried and reverted:
        snapping the limits moved the realized ratios enough to
        starve the Cp premise screen); sigma =
        round(sf*mu, dp+1) with sigma/mu STRICTLY inside (0.01, 0.05)
        in exact Fractions and sigma >= 1.5 microns; offsets built
        from a limiting ratio in [0.72, Cpk* - 0.07] and a second
        ratio at least 0.08 higher but capped so Cp stays under the
        target; band screened into [20, 110] microns and at least 10
        microns wide. RESISTOR — exact nominal +/- tolerance band as
        above; mean drifted off nominal by a fraction f of the
        half-band with f in [0.03, min(0.45, 0.9*(1 - lim/Cpk*))] so the Cp
        premise holds; sigma DERIVED as
        round(d_min/(3*limiting ratio), dp+1); sigma/mu strictly
        inside (0.005, 0.03). Target Cpk* from {1.33, 1.50, 1.67}.
        DISPLAY EXACTNESS (c1 blocking): d_min prints at
        max(dp, spec_dp) decimals — its EXACT value, asserted to
        round-trip — and Step 4 divides that printed value
        (lesson 83).
        DIRECTIONAL ROUNDING IS NOW CONSISTENT END TO END (c3 R3):
        sigma_max is floored AND the reported percentage is rounded
        UP, so applying the stated reduction genuinely reaches the
        target (asserted per draw); the c3 draft floored the ceiling
        but then rounded the percentage half-up, discarding that
        conservatism in ~35% of draws. The stem also states that the
        process is in statistical control (Montgomery Ch. 8 calls
        this precondition critical) and motivates the fixed mean.
        SIGMA_MAX PRECISION: the bound carries dp+2, one digit finer
        than sigma itself. c3 R1 read that as metrologically odd; it
        is deliberate and was re-tried at dp+1 in c4, which starved
        the Cpk* = 1.67 branch (1.3% of draws vs ~33% at dp+2)
        because a coarser bound collides with the reduction band and
        the sufficiency assert. A DERIVED bound quoted one digit
        beyond its datum is ordinary engineering practice, and the
        extra digit is load-bearing for the graded percentage.
        SIGMA_MAX DIRECTION (c2 R1+R3): sigma_max is rounded DOWN at
        dp+2, never half-up, so the printed bound genuinely meets the
        target (a half-up ceiling sat on the infeasible side in 43%
        of c2 draws while being described as the largest admissible
        sigma); the question states the downward rounding and a
        per-draw assert re-checks d_min/(3*sigma_max) >= Cpk*.
        WORDING: the reported quantity is a STANDARD-DEVIATION
        reduction and is named as such everywhere (c2 R1: the c2
        prose called it a variance reduction, a different number).
        Screens: Cp/Cpu/Cpl each >= 0.002 from their 2-dp boundary;
        |Cpu - Cpl| >= 0.03 exactly (decisive min, lesson 77);
        Cpk* - Cpk >= 0.05 and Cpk* - Cp >= 0.03 exactly (decisive
        shortfall and decisive premise);         reduction at least 0.002
        percentage points from its 1-dp boundary ON BOTH SIDES, with no
        exemption for values sitting exactly on the grid (c4 R1: the
        on-grid case is exactly where a float64 solver following a
        ceiling instruction flips), and sigma_max likewise screened off
        its floor boundary on both sides; AND path-agreement
        screened (lesson 80): a solver carrying the UNROUNDED
        sigma_max lands on the same 1-decimal percentage; reduction
        in [5.0, 60.0] percent. (The c2 "min-side consistency" screen
        was deleted as tautological — both ratios share the positive
        denominator 3*sigma — rather than left implying protection it
        never gave, c2 R3.)

    Returns:
        tuple(str, str): (question, solution)
    """
    from fractions import Fraction

    def _r2f(fr):
        cents = fr * 100
        fpart = cents - math.floor(cents)
        if abs(fpart - Fraction(1, 2)) < Fraction(1, 5):
            return None
        return (math.floor(cents)
                + (1 if fpart >= Fraction(1, 2) else 0))

    for _ in range(1200):
        cstar = random.choice(["1.33", "1.50", "1.67"])
        f_c = Fraction(cstar)
        c_val = float(f_c)
        lim_t = random.uniform(0.72, c_val - 0.07)
        is_res = random.choice([True, False])

        if is_res:
            cfg = _T26R_RESISTOR
            dp, sdp_spec = 0, 0
            tol = random.choice([5, 10])
            nominal = random.choice(_T26R_E24_5PCT if tol == 5
                                    else _T26R_E12_10PCT)
            step = nominal * tol // 100          # exact whole ohms
            assert step * 100 == nominal * tol, "band must be exact"
            lsl, usl = nominal - step, nominal + step
            # Cp = Cpk/(1-f) for a symmetric band: keep f under the
            # relative shortfall so Cp stays below the target
            f_hi = min(0.45, 0.9 * (1 - lim_t / c_val))
            if f_hi <= 0.03:
                continue
            drift = random.uniform(0.03, f_hi) * step
            mu = round(nominal + random.choice([-1, 1]) * drift)
            if not (lsl < mu < usl):
                continue
            d_lim = min(usl - mu, mu - lsl)
            if d_lim <= 0:
                continue
            sigma = round(d_lim / (3 * lim_t), dp + 1)
            band_note = (f"specified as a nominal {nominal} ohms with "
                         f"a {tol}% tolerance, i.e. LSL = {lsl} ohms "
                         f"and USL = {usl} ohms")
        else:
            cfg = _T26R_ANODIZE
            dp, sdp_spec = 1, 0
            lo_t, hi_t = cfg["tsub"]
            mu = round(random.uniform(lo_t, hi_t), dp)
            sf_lo, sf_hi = cfg["sf"]
            sigma = round(mu * random.uniform(sf_lo, sf_hi), dp + 1)
            if sigma < COATING_METROLOGY_FLOOR:
                continue
            # Cp = Cpk*(1+r)/2 with r = oth/lim: cap the second ratio
            # so Cp stays below the target
            d_hi = min(0.90, 0.9 * 2 * (c_val - lim_t))
            if d_hi <= 0.08:
                continue
            oth_t = lim_t + random.uniform(0.08, d_hi)
            off_lim = round(3 * sigma * lim_t, sdp_spec)
            off_oth = round(3 * sigma * oth_t, sdp_spec)
            if random.choice([True, False]):
                lsl, usl = round(mu - off_lim), round(mu + off_oth)
            else:
                lsl, usl = round(mu - off_oth), round(mu + off_lim)
            _an_lo, _an_hi = HARD_ANODIZE_THICKNESS
            if not (_an_lo <= lsl and usl <= _an_hi
                    and usl - lsl >= 10):
                continue
            band_note = (f"held to a specification of LSL = {lsl} "
                         f"microns and USL = {usl} microns")

        if sigma <= 0 or lsl <= 0 or usl <= lsl:
            continue
        mu_s = f"{mu:.{dp}f}"
        sig_s = f"{sigma:.{dp + 1}f}"
        sf_lo, sf_hi = cfg["sf"]
        ratio = Fraction(sig_s) / Fraction(mu_s)
        if not (Fraction(str(sf_lo)) < ratio < Fraction(str(sf_hi))):
            continue

        lsl_s = f"{lsl:.{sdp_spec}f}"
        usl_s = f"{usl:.{sdp_spec}f}"
        f_mu, f_sig = Fraction(mu_s), Fraction(sig_s)
        d_u = Fraction(usl_s) - f_mu
        d_l = f_mu - Fraction(lsl_s)
        if d_u <= 0 or d_l <= 0:
            continue
        cpu_f = d_u / (3 * f_sig)
        cpl_f = d_l / (3 * f_sig)
        cp_f = (d_u + d_l) / (6 * f_sig)
        if abs(cpu_f - cpl_f) < Fraction(3, 100):
            continue
        cpk_f = min(cpu_f, cpl_f)
        d_min = min(d_u, d_l)
        if not (Fraction(70, 100) <= cpk_f <= f_c - Fraction(5, 100)):
            continue
        # PREMISE: not even potentially capable (c2 R1)
        if cp_f > f_c - Fraction(3, 100):
            continue
        cpu2, cpl2, cp2 = _r2f(cpu_f), _r2f(cpl_f), _r2f(cp_f)
        if cpu2 is None or cpl2 is None or cp2 is None:
            continue
        cpk2 = min(cpu2, cpl2)

        dd = max(dp, sdp_spec)
        if d_min * (10 ** dd) != int(d_min * (10 ** dd)):
            continue
        dmin_s = f"{float(d_min):.{dd}f}"
        assert Fraction(dmin_s) == d_min, "d_min must display exactly"

        # sigma_max FLOORED at dp+2 so the printed ceiling is feasible
        smax_f = Fraction(dmin_s) / (3 * f_c)
        scale = 10 ** (dp + 2)
        sm_scaled = smax_f * scale
        sm_frac = sm_scaled - math.floor(sm_scaled)
        if sm_frac < Fraction(1, 50) or sm_frac > Fraction(49, 50):
            continue                      # c4 R1: floor boundary hole
        smax_i = math.floor(sm_scaled)
        smax_disp = Fraction(smax_i, scale)
        if smax_disp <= 0 or smax_disp >= f_sig:
            continue
        if Fraction(dmin_s) / (3 * smax_disp) < f_c:
            continue

        red_f = (1 - smax_disp / f_sig) * 100
        tenths = red_f * 10
        fpart = tenths - math.floor(tenths)
        if fpart < Fraction(1, 50) or fpart > Fraction(49, 50):
            continue                      # c4 R1: on-grid exemption
        red1 = math.ceil(tenths)
        red_ex = (1 - smax_f / f_sig) * 100
        t_ex = red_ex * 10
        fe = t_ex - math.floor(t_ex)
        if fe < Fraction(1, 50) or fe > Fraction(49, 50):
            continue
        if math.ceil(t_ex) != red1:
            continue
        if not 50 <= red1 <= 600:
            continue
        break
    else:
        raise AssertionError("resample loop exhausted")

    assert 0.70 <= float(cpk_f) <= float(f_c) - 0.049, "Cpk_now band"
    assert float(cp_f) <= float(f_c) - 0.029, "premise: Cp below target"
    assert 50 <= red1 <= 600, f"reduction out of band: {red1}"
    assert 0 < smax_disp < f_sig, "sigma_max must be a reduction"
    assert Fraction(dmin_s) / (3 * smax_disp) >= f_c, \
        "displayed sigma_max must meet the target"
    _sig_after = f_sig * (1 - Fraction(red1, 1000))
    assert Fraction(dmin_s) / (3 * _sig_after) >= f_c, \
        "the reported reduction must actually reach the target"

    recentre_note = (
        "Moving the mean would mean re-qualifying the bath schedule, "
        "which is deferred to a later process revision," if not is_res
        else "Moving the mean would require a new resistor paste lot "
        "and re-trim setup, which is deferred to a later process "
        "revision,")
    side = "upper" if d_u < d_l else "lower"
    sdp2 = dp + 2
    smax_s = f"{smax_i / scale:.{sdp2}f}"
    red_s = f"{red1 / 10:.1f}"
    tc_s = f"{3 * float(f_c):.2f}"

    question = (
        f"{cfg['phrase'][:1].upper()}{cfg['phrase'][1:]} is in "
        f"statistical control at a mean of mu = {mu_s} "
        f"{cfg['unit']}, with a process standard "
        f"deviation of sigma = {sig_s} {cfg['unit']}, {band_note}. "
        f"Capability is measured by Cp = (USL - LSL)/(6*sigma), "
        f"Cpu = (USL - mu)/(3*sigma), Cpl = (mu - LSL)/(3*sigma) and "
        f"Cpk = min(Cpu, Cpl); quality planning has set a target of "
        f"Cpk = {cstar}. {recentre_note} so for this "
        f"analysis the mean stays where it is and the target must be "
        f"met by reducing variation alone. Report Cp, Cpu, Cpl and "
        f"Cpk, each to 2 decimals (round half up). Then explain from "
        f"Cp whether re-centring alone could have reached the target. "
        f"Then report the largest process standard deviation that "
        f"would meet the target with the mean unchanged, to {sdp2} "
        f"decimals, rounded DOWN so that the value you state does "
        f"meet the target. Finally, report the percentage by which "
        f"the current standard deviation must be reduced to reach "
        f"that value, to 1 decimal, rounded UP so that the reduction "
        f"you state is sufficient."
    )

    step1 = (
        f"**Step 1:** Current capability ratios.\n"
        f"Cp = (USL - LSL)/(6*sigma) = ({usl_s} - {lsl_s}) / (6 * "
        f"{sig_s}) = {cp2 / 100:.2f};  "
        f"Cpu = (USL - mu)/(3*sigma) = ({usl_s} - {mu_s}) / (3 * "
        f"{sig_s}) = {cpu2 / 100:.2f};  "
        f"Cpl = (mu - LSL)/(3*sigma) = ({mu_s} - {lsl_s}) / (3 * "
        f"{sig_s}) = {cpl2 / 100:.2f};  "
        f"Cpk = min({cpu2 / 100:.2f}, {cpl2 / 100:.2f}) = "
        f"{cpk2 / 100:.2f}, short of the target {cstar}."
    )
    step2 = (
        f"**Step 2:** Check whether re-centring alone could fix it. "
        f"Cp is the "
        f"capability a perfectly centred process would achieve at "
        f"this spread, and Cp = {cp2 / 100:.2f} is itself below "
        f"{cstar}. So no amount of re-centring reaches the target at "
        f"the current spread — the standard deviation must come down."
    )
    step3 = (
        f"**Step 3:** With the mean fixed, Cpk is set by the nearer "
        f"specification limit, at distance\n"
        f"d_min = min(USL - mu, mu - LSL) = min({usl_s} - {mu_s}, "
        f"{mu_s} - {lsl_s}) = {dmin_s} {cfg['unit']} (the {side} "
        f"side), so Cpk = d_min/(3*sigma) throughout."
    )
    step4 = (
        f"**Step 4:** Invert that relation. Cpk = d_min/(3*sigma) "
        f"decreases as sigma grows, so Cpk >= {cstar} exactly when "
        f"sigma <= d_min/(3 * {cstar}).\n"
        f"sigma_max = {dmin_s} / {tc_s} = {smax_s} {cfg['unit']} "
        f"(rounded down, so this value still meets the target)"
    )
    step5 = (
        f"**Step 5:** Required reduction in the standard deviation "
        f"relative to its current value.\n"
        f"reduction = (1 - sigma_max/sigma) * 100 = "
        f"(1 - {smax_s}/{sig_s}) * 100 = {red_s} percent "
        f"(rounded up, so a reduction of this size does reach the "
        f"target)"
    )

    solution = (
        f"**Given:**\n"
        f"mu = {mu_s} {cfg['unit']}; sigma = {sig_s} {cfg['unit']}; "
        f"LSL = {lsl_s} {cfg['unit']}; USL = {usl_s} {cfg['unit']}; "
        f"target Cpk = {cstar}; mean held at its current value.\n\n"
        f"{step1}\n\n{step2}\n\n{step3}\n\n{step4}\n\n{step5}\n\n"
        f"**Answer:** The required reduction in the process standard "
        f"deviation is {red_s} percent"
    )

    return question, solution
