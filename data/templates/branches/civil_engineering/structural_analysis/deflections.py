import math
import random
from decimal import Decimal, ROUND_HALF_UP

from data.templates.branches.civil_engineering.constants import (
    AISC_W_SHAPES,
    STEEL_E_GPA,
    STEEL_E_KSI,
)


# Verifier tolerances (D1.3), keyed by template id.
#
# Each entry is the RELATIVE tolerance a deterministic verifier should apply
# when comparing a re-derived answer against this template's gold answer, on
# the understanding that the re-derived value is first quantised to the
# precision the trace quotes the answer to. That quantisation is the point:
# without it the tolerance has to absorb a full display step and stops testing
# anything (reviews/phase0_oracle_findings.md Section 6).
#
# After the Phase 1 fixes both items round-trip EXACTLY - measured worst
# relative error 0.0 over 1,000 seeds - so these values are headroom, not
# observed slack.
TEMPLATE_TOLERANCES = {
    # Answer quoted to 1 dp in mm over a 5-60 mm window (SI) or 3 dp in inches
    # over 0.2-2.5 in (US). One display step at the largest permitted
    # deflection is ~1.7e-3 relative, so 1e-3 still catches any genuine
    # single-step disagreement, while the real failure modes (wrong load-case
    # formula, dropped 1e6/1e-12 in the I conversion, missing ft->in span
    # conversion) are tens of percent.
    "template_beam_deflection_formula": 1e-3,
    # Same reasoning; the docstring caps tip deflection at [3, 45] mm, so one
    # 0.1 mm display step at the top of the window is 2.2e-3 relative.
    "template_cantilever_double_integration": 1e-3,
}


def _hu(x, places):
    """Round half-up to `places` dp, resolving the tie in DECIMAL.

    `round()` resolves a half-way tie on the binary value, so
    `round(0.01185 * 1000, 1)` gives 11.8 where a reader doing decimal
    arithmetic gets 11.9 (spec P2 as amended, DECISIONS D-012). Accepts a
    Decimal so an exact decimal product can be quantised without a detour
    through a float. Matches the industrial branch's `_hu` convention.
    """
    q = Decimal(1).scaleb(-places)
    d = x if isinstance(x, Decimal) else Decimal(repr(x))
    v = d.quantize(q, rounding=ROUND_HALF_UP)
    return int(v) if places == 0 else float(v)


def _shown(x, places):
    """The exact decimal value a reader sees when `x` is printed at `places` dp."""
    return Decimal(f"{x:.{places}f}")


def _as_printed(x, spec):
    """The value a reader recovers from `x` when it is printed with `spec`.

    P2 asks that the stored value and the printed value be the SAME value.
    Rounding alone does not achieve that: `34.4 * 1e-6` is
    3.4399999999999996e-05, prints as "3.4400e-05", and reparses to
    3.44e-05 - a different float, one ulp away. One ulp is enough to put the
    template's value and the reader's value on opposite sides of a display
    tie, which is how a non-closing Step-5 line survived the first round of
    this fix (Phase 1 review A, finding F-2; seeds 12639, 18910, 32721,
    47978). Binding through the display makes the two identical by
    construction, so the tie guard below tests the reader's value too.
    """
    return float(format(x, spec))


def _is_display_tie(x, places, rel_band=1e-12):
    """Is `x` at, or within a hair of, a half-way tie at `places` dp?

    A tie is the one case where NO rounding convention is defensible. A reader
    doing decimal arithmetic and applying half-up reads 0.02325 m as 23.3 mm; a
    reader using binary floats and `round()` reads it as 23.2 mm; and the
    printed line closes for exactly one of them whichever the template picks.
    Breaking the tie in decimal (P2 as amended) does not remove the ambiguity,
    it only moves it to the other reader - which is why such instances are
    RESAMPLED rather than resolved (D-016).

    Testing for an EXACT tie is not enough. Two independent evaluations of the
    same exact quantity - this template's, in kN and kN/m^2, and a solver's, in
    N and Pa - differ by a few ulps, so a rational tie such as 29.25 mm lands
    as 29.249999999999996 on one side and 29.250000000000004 on the other.
    Neither is exactly a tie, and the two then round in opposite directions:
    that is how eight non-closing instances per 60,000 seeds survived the first
    version of this guard (Phase 1 review A, finding F-2).

    So a narrow BAND is quarantined rather than a point. The band is a few
    thousand ulps wide where the tie lattice is 10^-places apart, so it removes
    nothing that is not genuinely ambiguous.
    """
    scaled = abs(x) * 10.0 ** places
    band = max(scaled * rel_band, 1e-9)
    return abs((scaled - math.floor(scaled)) - 0.5) <= band


def _fmt_g(x):
    """A value printed at its own precision - `82.7` stays 82.7, `170` stays 170.

    `f"{Ix:.0f}"` printed the AISC W8X24 second moment (82.7 in^4) as `83`
    while the arithmetic used 82.7, so the substitution line could not be
    reproduced from its own printed operands.
    """
    return f"{x:g}"


# Template 16 (Easy) — Area A3: Deflections
def template_beam_deflection_formula():
    """
    Midspan Deflection of a Simply Supported Steel Beam

    Scenario:
        A simply supported W-shape beam carries either a full-span uniform
        load or a central point load — the load case selects the standard
        elastic-deflection formula:

            uniform:  delta = 5*w*L^4 / (384*E*I)
            point:    delta = P*L^3 / (48*E*I)

        The instance is posed in either SI or US customary units (the
        AISC database carries both), so the unit-conversion chain also
        branches.

    Difficulty: Easy
    Grounding: Hibbeler, Structural Analysis, 10th ed. (SI), Ch. 8
        (Deflections — elastic-beam formulas); section properties from
        the AISC Shapes Database v16.0 (constants.AISC_W_SHAPES); E per
        AISC (200 GPa / 29,000 ksi).
    Physical bounds: span set from L/d in [15, 24]; the load is DERIVED
        from a target deflection ratio delta/L sampled inside a window
        jointly limited by serviceability (>= L/500) and a bending-stress
        screen (sigma <= ~180 MPa / 26 ksi, so the sampled load never
        overstresses the section it deflects) with a 1/245 ceiling for
        customary L/240 serviceability; recomputed delta/L asserted in
        [1/520, 1/238].

    Returns:
        tuple: (question, solution)
    """
    shape = random.choice(list(AISC_W_SHAPES.keys()))
    props = AISC_W_SHAPES[shape]
    use_si = random.choice([True, False])
    load_case = random.choice(["uniform", "point"])
    ld_ratio = random.uniform(15.0, 24.0)

    # Deflection-ratio window jointly bounded by serviceability and the
    # bending-stress screen (AUTHOR_NOTES lesson 14):
    #   point:   sigma ~ 6*E*(d/L)*r  -> r <= 1.5e-4 * (L/d)
    #   uniform: sigma ~ 4.8*E*(d/L)*r -> r <= 1.875e-4 * (L/d)
    # Margin 0.96 absorbs the whole-foot span rounding in the US branch
    # (R1, cycle 1: 0.98 let sigma reach 26.75 ksi); ceiling 1/245 keeps
    # the recomputed ratio inside customary L/240 serviceability.
    r_cap = (1.5e-4 if load_case == "point" else 1.875e-4) * ld_ratio * 0.96
    r_hi = min(r_cap, 1.0 / 245.0)

    if use_si:
        d_mm = props["si"]["d"]
        # Every operand that is displayed and then consumed is bound through
        # its own display, so the stored value IS the value a reader recovers
        # from the printed text (see _as_printed).
        Ix = _as_printed(props["si"]["Ix"] * 1e-6, ".4e")   # 10^6 mm^4 -> m^4
        E = _hu(STEEL_E_GPA * 1e6, 0)           # GPa -> kN/m^2
        L = _hu(ld_ratio * d_mm / 1000.0, 1)    # m
        shape_label = f"{shape} (SI designation {props['si_label']})"
        Ix_text = f"{props['si']['Ix']:.1f} x 10^6 mm^4"
        E_text = f"{STEEL_E_GPA} GPa"
        # Resample the load if the resulting deflection lands exactly on a
        # display tie (see _is_display_tie). Measured hit rate is well under
        # 1%, so this changes the draw sequence for almost no instance.
        for _attempt in range(200):
            r = random.uniform(1.0 / 500.0, r_hi)
            if load_case == "uniform":
                w = _hu(r * L * 384 * E * Ix / (5 * L ** 4), 1)
                delta = 5 * w * L ** 4 / (384 * E * Ix)            # m
            else:
                P = _hu(r * L * 48 * E * Ix / (L ** 3), 1)
                delta = P * L ** 3 / (48 * E * Ix)                 # m
            if not _is_display_tie(delta * 1000, 1):
                break
        else:
            raise RuntimeError("no tie-free load found in 200 draws")
        if load_case == "uniform":
            load_text = (f"a uniformly distributed load of {w:.1f} kN/m "
                         f"over the full span")
            assert 5.0 <= w <= 105.0, f"UDL implausible: {w}"
        else:
            load_text = f"a concentrated load of {P:.1f} kN at midspan"
            assert 15.0 <= P <= 900.0, f"point load implausible: {P}"
        # ONE rounding, at the precision the answer is quoted to, with the
        # metre display matched to it so the x1000 conversion is EXACT for
        # every reader (P2 as amended; D-016).
        #
        # This chain previously printed delta at 5 dp in metres and 1 dp in
        # mm. A 5-dp decimal times 1000 lands exactly on a 1-dp half-way tie
        # for ~10% of instances, and `round()` then resolved that tie on the
        # binary value: `round(0.01185 * 1000, 1)` gives 11.8 where decimal
        # half-up gives 11.9. Resolving the tie in decimal instead only moves
        # the disagreement to the other reader. Matching the precisions
        # removes the tie outright: 4 dp in metres IS 1 dp in mm, exactly.
        delta_out = _hu(delta * 1000, 1)                          # mm
        delta = _hu(delta_out / 1000.0, 4)                        # exact at 4 dp
        unit_out = "mm"
        span_text = f"L = {L:.1f} m"
    else:
        d_in = props["us"]["d"]
        Ix = _as_printed(props["us"]["Ix"], "g")   # in^4
        E = STEEL_E_KSI                          # ksi (an exact integer)
        L_ft = int(round(ld_ratio * d_in / 12.0))
        L_in = L_ft * 12
        shape_label = shape
        Ix_text = f"{_fmt_g(Ix)} in^4"
        E_text = f"{STEEL_E_KSI} ksi"
        L = L_in
        for _attempt in range(200):
            r = random.uniform(1.0 / 500.0, r_hi)
            if load_case == "uniform":
                w_try = _hu(
                    r * L_in * 384 * E * Ix / (5 * L_in ** 4) * 12, 2)
                d_try = 5 * (w_try / 12) * L_in ** 4 / (384 * E * Ix)
            else:
                P_try = _hu(r * L_in * 48 * E * Ix / (L_in ** 3), 1)
                d_try = P_try * L_in ** 3 / (48 * E * Ix)
            if not _is_display_tie(d_try, 3):
                break
        else:
            raise RuntimeError("no tie-free load found in 200 draws")
        if load_case == "uniform":
            w = _hu(
                r * L_in * 384 * E * Ix / (5 * L_in ** 4) * 12, 2)  # kip/ft
            # kip/ft -> kip/in is a division by 12, which does not terminate
            # in decimal, so ANY finite intermediate display is a rounding.
            # Computing from a 5-dp intermediate put the answer one unit in
            # the last place away from what the question implies on ~1.6% of
            # US uniform instances (T2 oracle, seeds 92/314/825/869).
            #
            # So Step 2 performs the conversion and states its value to 6 dp -
            # the step exists to demonstrate exactly this, and dropping it made
            # the step a no-op (Phase 1 review B, F-1) - while the SUBSTITUTION
            # carries the exact ratio, so nothing rounded sits between the
            # stated load and the answer (D-016).
            w_kip_in = _hu(w / 12, 6)          # DISPLAY ONLY - never consumed
            delta = 5 * (w / 12) * L_in ** 4 / (384 * E * Ix)      # in
            load_text = (f"a uniformly distributed load of {w:.2f} kip/ft "
                         f"over the full span")
            assert 0.3 <= w <= 7.5, f"UDL implausible: {w}"
        else:
            P = _hu(r * L_in * 48 * E * Ix / (L_in ** 3), 1)  # kips
            delta = P * L_in ** 3 / (48 * E * Ix)                 # in
            load_text = f"a concentrated load of {P:.1f} kips at midspan"
            assert 4.0 <= P <= 200.0, f"point load implausible: {P}"
        delta_out = _hu(delta, 3)
        delta = delta_out                    # the value the trace states
        unit_out = "in"
        span_text = f"L = {L_ft} ft = {L_in} in"

    # Recomputed serviceability ratio (from the presented, rounded load) —
    # dimensionless in either system.
    ratio_check = delta / (L if use_si else L_in)
    assert 1.0 / 520.0 <= ratio_check <= 1.0 / 238.0, (
        f"deflection ratio out of window: {ratio_check}")

    if load_case == "uniform":
        formula_step = (
            f"**Step 1:** Select the deflection formula.\n"
            f"For a simply supported beam under a full-span uniform load, "
            f"the maximum (midspan) deflection is:\n"
            f"delta = 5*w*L^4 / (384*E*I)\n\n"
        )
    else:
        formula_step = (
            f"**Step 1:** Select the deflection formula.\n"
            f"For a simply supported beam with a concentrated load at "
            f"midspan, the maximum (midspan) deflection is:\n"
            f"delta = P*L^3 / (48*E*I)\n\n"
        )

    if use_si:
        conv_step = (
            f"**Step 2:** Assemble consistent SI units.\n"
            f"E = {E_text} = {E:.0f} kN/m^2\n"
            f"I = {Ix_text} = {Ix:.4e} m^4\n"
            f"Span L = {L:.1f} m\n\n"
        )
        if load_case == "uniform":
            sub_step = (
                f"**Step 3:** Substitute and evaluate.\n"
                f"delta = 5 * {w:.1f} * ({L:.1f})^4 / (384 * {E:.0f} * "
                f"{Ix:.4e}) = {delta:.4f} m\n\n"
                f"**Step 4:** Express in millimeters.\n"
                f"delta = {delta:.4f} * 1000 = {delta_out:.1f} mm\n\n"
            )
        else:
            sub_step = (
                f"**Step 3:** Substitute and evaluate.\n"
                f"delta = {P:.1f} * ({L:.1f})^3 / (48 * {E:.0f} * "
                f"{Ix:.4e}) = {delta:.4f} m\n\n"
                f"**Step 4:** Express in millimeters.\n"
                f"delta = {delta:.4f} * 1000 = {delta_out:.1f} mm\n\n"
            )
    else:
        conv_step = (
            f"**Step 2:** Assemble consistent US customary units "
            f"(kips and inches).\n"
            f"E = {E_text}\n"
            f"I = {Ix_text}\n"
            f"Span L = {L_ft} ft = {L_in} in\n"
            + (f"w = {w:.2f} kip/ft = {w:.2f}/12 = {w_kip_in:.6f} kip/in\n\n"
               if load_case == "uniform" else "\n")
        )
        if load_case == "uniform":
            sub_step = (
                f"**Step 3:** Substitute and evaluate.\n"
                f"delta = 5 * ({w:.2f}/12) * ({L_in})^4 / (384 * {E} * "
                f"{_fmt_g(Ix)}) = {delta_out:.3f} in\n\n"
            )
        else:
            sub_step = (
                f"**Step 3:** Substitute and evaluate.\n"
                f"delta = {P:.1f} * ({L_in})^3 / (48 * {E} * {_fmt_g(Ix)}) "
                f"= {delta_out:.3f} in\n\n"
            )

    # A bare `{delta_out}` printed a repr, so the answer's precision drifted
    # with the value (0.38 vs 0.382 in). Fix it to the step's precision.
    delta_text = f"{delta_out:.1f}" if use_si else f"{delta_out:.3f}"

    question = (
        f"A simply supported steel beam consists of a {shape_label} "
        f"section (moment of inertia I = {Ix_text}, modulus of "
        f"elasticity E = {E_text}) spanning {span_text}. The beam "
        f"carries {load_text}. Determine the maximum deflection of the "
        f"beam in {'millimeters' if use_si else 'inches'}."
    )

    solution = (
        f"**Given:**\n"
        f"Section: {shape_label}, I = {Ix_text}\n"
        f"E = {E_text}; span {span_text}\n"
        f"Load: {load_text}\n\n"
        f"{formula_step}"
        f"{conv_step}"
        f"{sub_step}"
        f"**Answer:** The maximum deflection is {delta_text} {unit_out}"
    )

    return question, solution


# Template 17 (Intermediate) — Area A3: Deflections
def template_virtual_work_truss_deflection():
    """
    Truss Joint Deflection by the Method of Virtual Work

    Scenario:
        A two-member bracket truss (horizontal member AC, diagonal BC)
        carries a vertical load at joint C. The vertical deflection of C
        follows from the unit-load (virtual work) method:

            delta = sum( n_i * N_i * L_i / (A_i * E) )

        requiring one real-force analysis, one virtual-force analysis
        under a unit vertical load at C, and the tabulated sum.

    Difficulty: Intermediate
    Grounding: Hibbeler, Structural Analysis, 10th ed. (SI), Ch. 9
        (Deflections Using Energy Methods — method of virtual work,
        trusses); E per AISC steel (200 GPa).
    Physical bounds: geometry b in [2.0, 4.0] m, h in [1.5, 3.0] m;
        load P in [20, 60] kN; member areas (up to 3600 mm^2) sampled
        inside per-sample windows that jointly enforce the axial-stress
        screen (<= ~180 MPa, lesson 14) AND close the deflection band by
        construction; deflection asserted in [0.5, 15] mm.

    Returns:
        tuple: (question, solution)
    """
    E_kn_mm2 = STEEL_E_GPA                      # 200 GPa = 200 kN/mm^2
    b = round(random.uniform(2.0, 4.0), 1)
    h = round(random.uniform(1.5, 3.0), 1)
    P = random.randint(20, 60)

    L_ac = b
    L_bc = round(math.sqrt(b ** 2 + h ** 2), 3)

    # Real forces (tension positive) at joint C.
    N_bc = round(P * L_bc / h, 2)               # tension
    N_ac = round(-P * b / h, 2)                 # compression

    # Virtual forces under a unit vertical downward load at C.
    n_bc = round(L_bc / h, 4)
    n_ac = round(-b / h, 4)

    # Per-sample area windows close BOTH the stress screen and the
    # deflection band by construction (R1, cycle 1: stress floors alone
    # left a ~0.15% crash corner at delta ~ 19.6 mm). With
    # c_i = n_i*N_i*(L_i in mm)/E, any A_i >= c_sum/14.0 guarantees
    # delta <= ~14.5 mm and any A_i <= c_sum/0.55 guarantees
    # delta >= ~0.52 mm; both windows are provably non-empty over the
    # whole parameter space (worst case, member AC at b=2, h=3:
    # c_sum/0.55 >= ~12.1*|N_AC| vs the 5.83*|N| stress floor; verified
    # exhaustively by R1 over 13,776 parameter combos).
    c_bc = n_bc * N_bc * (L_bc * 1000) / E_kn_mm2
    c_ac = n_ac * N_ac * (L_ac * 1000) / E_kn_mm2
    c_sum = c_bc + c_ac

    def _sample_area(force):
        lo = max(500.0, 1.05 * abs(force) / 0.18, c_sum / 14.0)
        hi = min(lo + 600.0, c_sum / 0.55, 3600.0)
        assert lo < hi, f"empty area window: {lo}, {hi}"
        return int(round(random.uniform(lo, hi) / 10) * 10)

    A_bc = _sample_area(N_bc)
    A_ac = _sample_area(N_ac)

    # Terms of the virtual-work sum in consistent mm units: L in mm, A in
    # mm^2, E in kN/mm^2, N in kN -> delta in mm.
    term_bc = round(n_bc * N_bc * (L_bc * 1000) / (A_bc * E_kn_mm2), 3)
    term_ac = round(n_ac * N_ac * (L_ac * 1000) / (A_ac * E_kn_mm2), 3)
    delta = round(term_bc + term_ac, 3)

    assert abs(N_bc) / A_bc <= 0.185 and abs(N_ac) / A_ac <= 0.185, (
        "axial stress screen violated")
    assert 0.5 <= delta <= 15.0, f"deflection out of bounds: {delta}"

    question = (
        f"A two-member steel bracket truss lies in a vertical plane. "
        f"Joint A is pinned to a wall, joint B is pinned "
        f"to the wall directly above A at height h = {h:.1f} m, and "
        f"joint C is located a horizontal distance b = {b:.1f} m from A "
        f"at the same level as A. Member AC is horizontal "
        f"(cross-sectional area {A_ac:.0f} mm^2) and member BC is the "
        f"diagonal (area {A_bc:.0f} mm^2). A vertical downward load of "
        f"P = {P} kN acts at joint C. The members are adequately braced, "
        f"so buckling is not a concern. Taking E = {STEEL_E_GPA} GPa for "
        f"both members, use the method of virtual work to determine the "
        f"vertical deflection of joint C in mm."
    )

    solution = (
        f"**Given:**\n"
        f"Geometry: b = {b:.1f} m, h = {h:.1f} m; L_AC = {L_ac:.1f} m, "
        f"L_BC = sqrt(b^2 + h^2) = {L_bc:.3f} m\n"
        f"Areas: A_AC = {A_ac:.0f} mm^2, A_BC = {A_bc:.0f} mm^2; "
        f"E = {STEEL_E_GPA} GPa = {E_kn_mm2} kN/mm^2\n"
        f"Load (P): {P} kN downward at C\n\n"
        f"**Step 1:** Real member forces (tension positive), from "
        f"equilibrium of joint C.\n"
        f"Vertical: N_BC * (h / L_BC) = P  ->  N_BC = P * L_BC / h "
        f"= {P} * {L_bc:.3f} / {h:.1f} = {N_bc:.2f} kN (tension)\n"
        f"Horizontal: N_AC = -N_BC * (b / L_BC) = -P * b / h "
        f"= -{P} * {b:.1f} / {h:.1f} = {N_ac:.2f} kN (compression)\n\n"
        f"**Step 2:** Virtual member forces under a unit vertical "
        f"downward load at C (so a positive result means a downward "
        f"deflection); the joint equilibrium repeats with P = 1:\n"
        f"n_BC = L_BC / h = {L_bc:.3f} / {h:.1f} = {n_bc:.4f}\n"
        f"n_AC = -b / h = -{b:.1f} / {h:.1f} = {n_ac:.4f}\n\n"
        f"**Step 3:** Tabulate n*N*L/(A*E) for each member (L in mm).\n"
        f"BC: {n_bc:.4f} * {N_bc:.2f} * {L_bc * 1000:.0f} / "
        f"({A_bc:.0f} * {E_kn_mm2}) = {term_bc:.3f} mm\n"
        f"AC: ({n_ac:.4f}) * ({N_ac:.2f}) * {L_ac * 1000:.0f} / "
        f"({A_ac:.0f} * {E_kn_mm2}) = {term_ac:.3f} mm\n\n"
        f"**Step 4:** Sum the contributions.\n"
        f"delta_C = {term_bc:.3f} + {term_ac:.3f} = {delta:.3f} mm "
        f"(downward)\n\n"
        f"**Answer:** The vertical deflection of joint C is "
        f"{delta:.3f} mm"
    )

    return question, solution


# Template 18 (Intermediate) — Area A3: Deflections
def template_cantilever_double_integration():
    """
    Cantilever Tip Deflection by Double Integration

    Scenario:
        A cantilever W-shape beam, fixed at A, carries either a tip point
        load or a full-length uniform load — the load case changes the
        moment function and therefore the entire integration:

            point:   M(x) = -P*(L - x)   ->  delta_tip = P*L^3 / (3*E*I)
            uniform: M(x) = -w*(L - x)^2 / 2  ->  delta_tip = w*L^4 / (8*E*I)

        The trace performs the double integration explicitly with the
        fixed-end boundary conditions v(0) = v'(0) = 0.

    Difficulty: Intermediate
    Grounding: Hibbeler, Structural Analysis, 10th ed. (SI), Ch. 8
        (Deflections — double-integration method); section properties per
        AISC Shapes Database (SI values); E = 200 GPa.
    Physical bounds: cantilever length from L/d in [8, 14]; load DERIVED
        from a target tip-deflection ratio inside a window jointly limited
        by serviceability and a bending-stress screen at the fixed end
        (sigma <= ~180 MPa; point: r <= 6.0e-4*(L/d), uniform:
        r <= 4.5e-4*(L/d)); recomputed delta in [3, 45] mm.

    Returns:
        tuple: (question, solution)
    """
    shape = random.choice(list(AISC_W_SHAPES.keys()))
    props = AISC_W_SHAPES[shape]
    load_case = random.choice(["point", "uniform"])
    ld_ratio = random.uniform(8.0, 14.0)

    d_mm = props["si"]["d"]
    # Bound through their own displays - see _as_printed.
    Ix = _as_printed(props["si"]["Ix"] * 1e-6, ".4e")     # m^4
    E = _hu(STEEL_E_GPA * 1e6, 0)             # kN/m^2
    L = _hu(ld_ratio * d_mm / 1000.0, 1)

    # Per-sample feasibility (lesson 1): the ratio window is jointly
    # bounded by serviceability (1/350..1/90), the stress screen (r_cap),
    # AND the absolute tip-deflection band [3, 45] mm, which depends on L.
    r_cap = (6.0e-4 if load_case == "point" else 4.5e-4) * ld_ratio * 0.98
    r_lo = max(1.0 / 350.0, 1.02 * 0.003 / L)
    r_hi = min(1.0 / 90.0, r_cap, 0.98 * 0.045 / L)
    assert r_lo < r_hi, f"empty ratio window: {r_lo}, {r_hi}"

    # Resample a load whose deflection lands exactly on the 1-dp mm display
    # tie; such an instance has no defensible gold answer (D-016).
    for _attempt in range(200):
        r = random.uniform(r_lo, r_hi)
        if load_case == "point":
            P = _hu(r * L * 3 * E * Ix / (L ** 3), 1)
            delta = P * L ** 3 / (3 * E * Ix)
        else:
            w = _hu(r * L * 8 * E * Ix / (L ** 4), 1)
            delta = w * L ** 4 / (8 * E * Ix)
        if not _is_display_tie(delta * 1000, 1):
            break
    else:
        raise RuntimeError("no tie-free load found in 200 draws")

    if load_case == "point":
        assert 3.0 <= P <= 600.0, f"tip load implausible: {P}"
        load_text = f"a concentrated load of P = {P:.1f} kN at the free end"
        m_step = (
            f"**Step 1:** Write the internal moment function (x measured "
            f"from the fixed end A).\n"
            f"M(x) = -P * (L - x) = -{P:.1f} * ({L:.1f} - x)\n\n"
        )
        int_steps = (
            f"**Step 2:** Integrate E*I*v'' = M(x) once for the slope.\n"
            f"E*I*v' = -P*(L*x - x^2/2) + C1; the fixed end gives "
            f"v'(0) = 0, so C1 = 0.\n\n"
            f"**Step 3:** Integrate again for the deflection.\n"
            f"E*I*v = -P*(L*x^2/2 - x^3/6) + C2; v(0) = 0 gives C2 = 0.\n\n"
            f"**Step 4:** Evaluate at the tip, x = L.\n"
            f"E*I*v(L) = -P*(L^3/2 - L^3/6) = -P*L^3/3, so the downward "
            f"tip deflection is delta = P*L^3 / (3*E*I)\n\n"
        )
        sub_expr = (f"{P:.1f} * ({L:.1f})^3 / "
                    f"(3 * {E:.0f} * {Ix:.4e})")
    else:
        assert 1.5 <= w <= 120.0, f"UDL implausible: {w}"
        load_text = (f"a uniformly distributed load of w = {w:.1f} kN/m "
                     f"over its entire length")
        m_step = (
            f"**Step 1:** Write the internal moment function (x measured "
            f"from the fixed end A).\n"
            f"M(x) = -w * (L - x)^2 / 2 = -{w:.1f} * ({L:.1f} - x)^2 / "
            f"2\n\n"
        )
        int_steps = (
            f"**Step 2:** Integrate E*I*v'' = M(x) once for the slope.\n"
            f"E*I*v' = w*(L - x)^3/6 + C1; v'(0) = 0 gives "
            f"C1 = -w*L^3/6.\n\n"
            f"**Step 3:** Integrate again for the deflection.\n"
            f"E*I*v = -w*(L - x)^4/24 - w*L^3*x/6 + C2; v(0) = 0 gives "
            f"C2 = w*L^4/24.\n\n"
            f"**Step 4:** Evaluate at the tip, x = L.\n"
            f"E*I*v(L) = -w*L^4/6 + w*L^4/24 = -w*L^4/8, so the downward "
            f"tip deflection is delta = w*L^4 / (8*E*I)\n\n"
        )
        sub_expr = (f"{w:.1f} * ({L:.1f})^4 / "
                    f"(8 * {E:.0f} * {Ix:.4e})")

    delta_mm = _hu(delta * 1000, 1)
    delta = _hu(delta_mm / 1000.0, 4)            # exact at 4 dp; what Step 5 states
    sub_step = (
        f"**Step 5:** Substitute numerical values.\n"
        f"delta = {sub_expr} = {delta:.4f} m\n\n"
    )
    assert 3.0 <= delta_mm <= 45.0, f"tip deflection out of bounds: {delta_mm}"
    # The x1000 conversion Step 6 prints must be exact in decimal, not merely
    # within tolerance - this is the property the phase exists to restore.
    assert Decimal(f"{delta:.4f}") * 1000 == Decimal(f"{delta_mm:.1f}"), (
        f"Step 6 conversion does not close: {delta} -> {delta_mm}")

    question = (
        f"A cantilever steel beam is built from a {shape} section (SI "
        f"designation {props['si_label']}; moment of inertia I = "
        f"{props['si']['Ix']:.1f} x 10^6 mm^4) and extends L = {L:.1f} m "
        f"from its fixed support at A. It carries {load_text}. Taking "
        f"E = {STEEL_E_GPA} GPa, use the double-integration method to "
        f"determine the deflection at the free end in mm."
    )

    solution = (
        f"**Given:**\n"
        f"Section: {shape}, I = {props['si']['Ix']:.1f} x 10^6 mm^4 = "
        f"{Ix:.4e} m^4\n"
        f"E = {STEEL_E_GPA} GPa = {E:.0f} kN/m^2; length L = {L:.1f} m\n"
        f"Load: {load_text}\n\n"
        f"{m_step}"
        f"{int_steps}"
        f"{sub_step}"
        f"**Step 6:** Express in millimeters.\n"
        f"delta = {delta:.4f} * 1000 = {delta_mm:.1f} mm\n\n"
        f"**Answer:** The deflection at the free end is {delta_mm:.1f} mm"
    )

    return question, solution
