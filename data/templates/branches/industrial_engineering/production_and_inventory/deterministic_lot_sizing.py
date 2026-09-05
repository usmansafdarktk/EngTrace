import math
import random
from decimal import Decimal, ROUND_HALF_UP

from pilot.templates.branches.industrial_engineering.constants import (
    HOLDING_RATE_PER_YR,
    INVENTORY_ITEMS,
)


def _hu(x, places):
    """Half-up rounding of a float via its shortest decimal repr."""
    q = Decimal("1") if places == 0 else Decimal("0." + "0" * (places - 1) + "1")
    v = Decimal(repr(x)).quantize(q, rounding=ROUND_HALF_UP)
    return int(v) if places == 0 else float(v)


# Item classes whose constants.py windows keep the per-sample K window
# non-empty under the order-frequency screen (lesson 1 joint feasibility).
# "machine spare part" is excluded because its low corner (c=20, i=0.15,
# D=100) EMPTIES the window (k_lo=60 > floor(D*h/8)=37); "retail
# appliance" is excluded conservatively to keep lot-size scales within
# the three corner-verified classes (its windows are feasible but were
# not corner-audited for realism).
_T11_CLASSES = {
    "electronic component": "a purchased electronic component",
    "packaged food case": "a packaged food case SKU",
    "industrial fastener": "an industrial fastener",
}


# Template 11 (Easy) — Area P1: Deterministic Lot Sizing
def template_basic_eoq():
    """
    Economic Order Quantity: Optimal Lot Size

    Scenario:
        A stockroom replenishes an item with constant annual demand D,
        fixed cost K per order, unit cost c, and an annual holding-cost
        rate i (fraction of unit value per year):

            h  = i * c                     (holding cost, $/unit/year)
            Q* = sqrt(2 K D / h)           (EOQ)

        Requested: the optimal order quantity, rounded to the nearest
        whole unit. The trace also states the implied ordering frequency
        as a sanity check.

    Difficulty: Easy
    Grounding: Nahmias, Production and Operations Analysis, 7th ed.,
        Sec. 4.5 (the EOQ model; visually verified in the on-disk copy,
        text pp. 210-213) with holding cost h = I*c per Sec. 4.4
        (pp. 204-205). Cross-ref H&L 7e Ch. 19 (basic EOQ model).
    Physical bounds: item class from constants.py INVENTORY_ITEMS
        (classes listed in _T11_CLASSES); c (2 dp) and D (integer) from
        the class windows; i (2 dp) in HOLDING_RATE_PER_YR [0.15, 0.40];
        K integer in the per-sample window
        [max(K_lo, ceil(D*h/21632)), min(K_hi, floor(D*h/8))], which
        pins the ordering frequency D/Q* = sqrt(D*h/(2K)) into
        [2, 104] orders/year (window non-emptiness verified at all class
        corners). Draws whose exact root lies within 0.006 of a .5
        integer boundary are resampled (bounded deterministic loop,
        ~1.2% rejection), so the 2-dp displayed root and a full-precision
        solve always round to the SAME whole unit (lesson 51). ANALYTIC
        corner envelope (Q* monotone in K, D up and h down): min
        Q* = sqrt(2*40*2000/16.00) = 100.0 (electronic component at
        c = 40, i = 0.40, D = 2000, K = 40); max Q* =
        sqrt(2*150*200000/0.075) = 28284.3 (fastener at c = 0.50,
        i = 0.15, D = 200000, K = 150) — bulk fastener lots legitimately
        reach tens of thousands of units. Asserts: Q* in [95, 28500];
        frequency in [1.9, 104.5].

    Returns:
        tuple(str, str): (question, solution)
    """
    # Bounded deterministic resample loop: reject draws whose exact root
    # sits within 0.006 of a .5 integer boundary (keeps the 2-dp display
    # and full precision rounding to the same whole unit; measured ~1.2%
    # rejection, so the 50-attempt cap is never realistically reached).
    for _ in range(50):
        key = random.choice(sorted(_T11_CLASSES))
        cw = INVENTORY_ITEMS[key]
        c = round(random.uniform(*cw["unit_cost_usd"]), 2)
        D = random.randint(*cw["annual_demand"])
        i = round(random.uniform(*HOLDING_RATE_PER_YR), 2)
        # h is exact at 4 dp (2-dp * 2-dp); per-sample K window per
        # docstring.
        h = round(i * c, 4)
        k_lo = max(cw["order_cost_usd"][0], math.ceil(D * h / 21632))
        k_hi = min(cw["order_cost_usd"][1], math.floor(D * h / 8))
        K = random.randint(k_lo, k_hi)
        s_exact = math.sqrt(2 * K * D / h)
        if abs(s_exact - math.floor(s_exact) - 0.5) > 0.006:
            break
    setting = _T11_CLASSES[key]

    # Round-then-recompute: Q from the displayed exact h; the raw root is
    # shown at 2 dp, then rounded half-up to whole units (lesson 51).
    Q_raw = round(s_exact, 2)
    Q = _hu(Q_raw, 0)
    freq = round(D / Q, 1)

    assert 0.15 <= i <= 0.40, f"i out of bounds: {i}"
    assert 95 <= Q <= 28500, f"Q out of bounds: {Q}"
    assert 1.9 <= freq <= 104.5, f"frequency out of bounds: {freq}"
    assert Q == _hu(s_exact, 0), f"display-chain flip: {Q_raw} vs {s_exact}"

    question = (
        f"A stockroom manages {setting} with a steady annual demand of "
        f"{D} units. Each replenishment order costs ${K} to place and "
        f"process, the item costs ${c:.2f} per unit, and the annual "
        f"holding-cost rate is {i:.2f} (a fraction of unit value per "
        f"year). Using the EOQ model, determine the optimal order "
        f"quantity, rounded to the nearest whole unit. In your solution, "
        f"first compute the holding cost per unit per year, and finish by "
        f"stating the implied number of orders per year."
    )

    solution = (
        f"**Given:**\n"
        f"Annual demand D = {D} units/year; order cost K = ${K}; unit "
        f"cost c = ${c:.2f}; holding rate i = {i:.2f} per year.\n\n"
        f"**Step 1:** Compute the holding cost per unit per year.\n"
        f"h = i * c = {i:.2f} * {c:.2f} = {h:.4f} $/unit/year\n\n"
        f"**Step 2:** Apply the EOQ formula.\n"
        f"Q* = sqrt(2*K*D / h) = sqrt(2 * {K} * {D} / {h:.4f}) "
        f"= {Q_raw:.2f} units\n\n"
        f"**Step 3:** Round to a practical lot size and check the "
        f"ordering frequency.\n"
        f"Q* = {Q} units (nearest whole unit); the stockroom would place "
        f"about D / Q* = {D} / {Q} = {freq:.1f} orders per year, a "
        f"reasonable frequency.\n\n"
        f"**Answer:** The optimal order quantity is {Q} units"
    )

    return question, solution


# Production-appropriate framings for the EPQ template (cycle-2 fix: the
# EOQ display string "a purchased electronic component" contradicted the
# in-house production premise — R1/R3 c1).
_T12_SETTINGS = {
    "electronic component": "an electronic control module",
    "packaged food case": "a packaged food case SKU",
    "industrial fastener": "an industrial fastener",
}


# Template 12 (Easy) — Area P1: Deterministic Lot Sizing
def template_epq_finite_production():
    """
    Economic Production Quantity: Optimal Run Size at a Finite
    Production Rate

    Scenario:
        An item is produced in-house at a finite rate P > D rather than
        delivered all at once, so stock builds up only at rate P - D
        while producing. The EOQ logic applies with the reduced
        effective holding cost:

            h  = i * c
            h' = h * (1 - D/P)
            Q* = sqrt(2 K D / h')          (EPQ)

        Requested: the optimal production run size, rounded to the
        nearest TEN units (a practical run-size grid). The trace also
        reports the maximum on-hand inventory H = Q* * (1 - D/P),
        closing the loop on where the (1 - D/P) factor comes from.

    Difficulty: Easy
    Grounding: Nahmias, Production and Operations Analysis, 7th ed.,
        Sec. 4.6 "Extension to a Finite Production Rate" (visually
        verified in the on-disk copy, text p. 218); h = I*c per
        Sec. 4.4. Cross-ref H&L 7e Ch. 19 (EOQ with continuous
        replenishment typology).
    Physical bounds: item classes as in template_basic_eoq, with two
        cycle-2 tightenings that bound the display-chain wobble
        SCALE-AWARELY (R1/R2 c1): (a) c is floored per-sample at
        0.12/i + 0.01 so h = i*c >= ~0.12 (binds only for fasteners),
        and (b) the production-rate multiple m (1 dp, <= 6.0) is floored
        so h' = h*(1 - 1/m) >= 0.1 as well as satisfying the K-window
        feasibility floor(D*h'/8) >= K_lo. With r = D/P displayed at
        6 dp and h' at 6 dp, the chain shifts the raw root by at most
        Q*(dr/(2(1-r)) + dh'/(2h')) <= Q*(7.5e-7 + 2.5e-6) <= 0.09
        units at the 24,495-unit envelope max — so the 0.8-unit
        ten-unit-grid midpoint screen dominates the wobble with ~9x
        margin and the displayed chain and a full-precision solve
        always round to the SAME ten-unit answer. K window as in
        template_basic_eoq (frequency screen [2, 104]/yr; non-empty
        since frac_needed = max(8*K_lo/(D*h), 0.1/h) <= 5/6 at every
        class corner). ANALYTIC corner envelope (author QA 2026-08-06,
        corner search over class/i/c-floor/D/m/K endpoints): raw Q* in
        [109.5, 24494.9] — min at (electronic, c=40, i=0.40, D=2000,
        m=6.0, K=40), max at (fastener, c=0.50, i=0.40, D=200000,
        m=2.0, K=150). Asserts: Q in [100, 24700]; r = D/P in
        [0.16, 0.67]; h' >= 0.0999.

    Returns:
        tuple(str, str): (question, solution)
    """
    # Bounded deterministic resample loop; every draw is fully assembled
    # inside its own iteration and stored only if it passes ALL screens
    # (cycle-2 hygiene fix: no stale-state fallthrough is possible).
    draw = None
    for _ in range(50):
        key = random.choice(sorted(_T12_SETTINGS))
        cw = INVENTORY_ITEMS[key]
        i = round(random.uniform(*HOLDING_RATE_PER_YR), 2)
        c_lo = max(cw["unit_cost_usd"][0], 0.12 / i + 0.01)
        c = round(random.uniform(c_lo, cw["unit_cost_usd"][1]), 2)
        D = random.randint(*cw["annual_demand"])
        h = round(i * c, 4)
        frac_needed = max(8 * cw["order_cost_usd"][0] / (D * h), 0.1 / h)
        if frac_needed >= (1 - 1 / 6):
            continue
        m_lo = max(1.5, math.ceil((1 / (1 - frac_needed)) * 10) / 10)
        m = round(random.uniform(m_lo, 6.0), 1)
        P = int(round(m * D))
        r_exact = D / P
        hp_exact = h * (1 - r_exact)
        if hp_exact < 0.0999:
            continue
        k_lo = max(cw["order_cost_usd"][0], math.ceil(D * h / 21632))
        k_hi = min(cw["order_cost_usd"][1], math.floor(D * hp_exact / 8))
        if k_lo > k_hi:
            continue
        K = random.randint(k_lo, k_hi)
        s_exact = math.sqrt(2 * K * D / hp_exact)
        if abs((s_exact % 10) - 5.0) > 0.8:
            draw = (key, c, D, i, h, P, K, s_exact)
            break
    assert draw is not None, "resample loop exhausted"
    key, c, D, i, h, P, K, s_exact = draw
    setting = _T12_SETTINGS[key]

    # Round-then-recompute: r at 6 dp, h' at 6 dp, root at 2 dp, final on
    # the ten-unit grid half-up (wobble bound derived in the docstring).
    r = round(D / P, 6)
    hp = round(h * (1 - r), 6)
    Q_raw = round(math.sqrt(2 * K * D / hp), 2)
    Q = _hu(Q_raw / 10, 0) * 10
    H_max = _hu(Q * (1 - r), 0)

    assert 0.16 <= r <= 0.67, f"r out of bounds: {r}"
    assert hp >= 0.0999, f"hp out of bounds: {hp}"
    assert 100 <= Q <= 24700, f"Q out of bounds: {Q}"
    assert Q == _hu(s_exact / 10, 0) * 10, f"grid flip: {Q_raw} vs {s_exact}"

    question = (
        f"A plant produces {setting} in-house at a rate of {P} units per "
        f"year, while steady demand consumes it at {D} units per year. "
        f"Each production run costs ${K} to set up, the item is valued "
        f"at ${c:.2f} per unit, and the annual holding-cost rate is "
        f"{i:.2f} (a fraction of unit value per year). Because the item "
        f"is produced gradually rather than delivered all at once, use "
        f"the finite-production-rate (EPQ) model to determine the "
        f"optimal production run size, rounded to the nearest ten "
        f"units. In your solution, compute the holding cost per unit per "
        f"year, the effective holding cost reflecting the finite rate, "
        f"and finish by stating the maximum on-hand inventory."
    )

    solution = (
        f"**Given:**\n"
        f"Production rate P = {P} units/year; demand D = {D} units/year; "
        f"setup cost K = ${K}; unit cost c = ${c:.2f}; holding rate "
        f"i = {i:.2f} per year.\n\n"
        f"**Step 1:** Compute the holding cost per unit per year.\n"
        f"h = i * c = {i:.2f} * {c:.2f} = {h:.4f} $/unit/year\n\n"
        f"**Step 2:** Compute the effective holding cost. During a run, "
        f"inventory builds at rate P - D, so a run of size Q peaks at "
        f"only H = Q*(1 - D/P) on hand; the average inventory is H/2 "
        f"rather than Q/2, which is equivalent to charging the reduced "
        f"holding cost h' = h*(1 - D/P) against Q/2:\n"
        f"D/P = {D} / {P} = {r:.6f};  "
        f"h' = h * (1 - D/P) = {h:.4f} * (1 - {r:.6f}) = {hp:.6f} "
        f"$/unit/year\n\n"
        f"**Step 3:** Apply the EPQ formula and round to a practical "
        f"run size.\n"
        f"Q* = sqrt(2*K*D / h') = sqrt(2 * {K} * {D} / {hp:.6f}) "
        f"= {Q_raw:.2f}, so Q* = {Q} units (nearest ten units)\n\n"
        f"**Step 4:** State the maximum on-hand inventory.\n"
        f"H = Q* * (1 - D/P) = {Q} * (1 - {r:.6f}) = {H_max} units — the "
        f"gap between H and Q* is exactly the finite-production-rate "
        f"effect from Step 2.\n\n"
        f"**Answer:** The optimal production run size is {Q} units"
    )

    return question, solution


def _round50(x):
    """Round to the nearest multiple of 50 (breakpoint cosmetics)."""
    return int(round(x / 50.0)) * 50


# Template 13 (Intermediate) — Area P1: Deterministic Lot Sizing
# [BRANCHING: all-units discount-regime decision]
def template_quantity_discount_all_units():
    """
    All-Units Quantity Discount: Optimal Order Quantity

    Scenario:
        A supplier offers an all-units discount schedule with three
        price tiers (c0 below b1 units, c1 from b1 up to b2 - 1, c2 at
        b2 and above; c2 < c1 < c0). The standard all-units procedure:
        compute the EOQ at each tier price, Q_j = sqrt(2KD/(i*c_j));
        check which EOQs are REALIZABLE (fall inside their own tier);
        evaluate the total annual cost

            G(Q) = c_j(Q)*D + K*D/Q + i*c_j(Q)*Q/2

        at the highest-tier realizable EOQ and at every higher
        breakpoint; choose the cheapest. The sampled schedule BRANCHES
        the reasoning path into three classes:
          (a) Q2 realizable -> optimal outright (no comparison needed:
              lowest price AND interior optimum);
          (b) only Q1 realizable -> compare G(Q1) vs G(b2);
          (c) only Q0 realizable -> compare G(Q0), G(b1), G(b2).
        Requested: the optimal order quantity in whole units.

    Difficulty: Intermediate
    Grounding: Nahmias, Production and Operations Analysis, 7th ed.,
        Sec. 4.7 "Quantity Discount Models" (all-units schedules;
        visually verified in the on-disk copy, text pp. 220-221; Example
        4.4 typology only, never its numbers). Cross-ref H&L 7e Ch. 19.
    Physical bounds: item class electronic component or packaged food
        case; c0 (2 dp) sampled >= $8 within the class window; tier
        prices c1 = f1*c0, c2 = f2*c0 (2 dp) with f1 in [0.90, 0.96],
        f2 in [0.80, 0.88] (within the cumulative DISCOUNT_STEP_FRACTION
        span) and gaps c0 - c1 >= 0.02, c1 - c2 >= 0.01; D, i, K as in
        template_basic_eoq (frequency screen on the base-tier EOQ).
        Breakpoints are CONSTRUCTED from the tier EOQs per the sampled
        branch class (multiples of 50, 100 <= b1 < b2 <= 5000, inside
        DISCOUNT_BREAK_QTY), making every class reachable by design.
        Discount steps are SMALL (2-3.5% per break, inside
        DISCOUNT_STEP_FRACTION) and breakpoints in the comparison
        branches are sampled around the EXACT break-even quantities
        (the quadratic solving G_tier(b) = G_reference), so both/all
        comparison winners are genuinely reachable — with 12-20% steps
        the deepest discount won every comparison (degenerate decision;
        caught by the cycle-1 authoring QA sweep). A per-sample cap on
        D keeps the break-even scale below the 5000-unit breakpoint
        ceiling. Decisiveness screens (lessons 23/51/57): every tier
        EOQ sits >= 4% away from every breakpoint (1-dp EOQ displays
        cannot flip realizability); the top two candidate costs differ
        by >= 0.3%; an interior winning EOQ sits >= 0.06 from a .5
        integer boundary, so the displayed chain and a full-precision
        solve give the SAME whole-unit answer. Analytic dominance
        ceiling: b2 <= 5000 caps breakpoint winners at 5000 and
        branch-(a) interior winners at Q2 <= 5000/0.45 ~ 11111.
        Analytic floor: the base-tier EOQ shares template_basic_eoq's
        min corner (~100), so interior winners can reach ~100. Author
        QA (40,000-seed sweep, final design): branch mix and winner mix
        recorded in the review log with every comparison winner holding
        a >= 20% share within its branch; attained Q_opt spanned
        [121, 5000]. Asserts (analytic dominance, lesson 9): Q_opt in
        [95, 11150].

    Returns:
        tuple(str, str): (question, solution)
    """
    for _ in range(300):
        key = random.choice(["electronic component", "packaged food case"])
        cw = INVENTORY_ITEMS[key]
        i = round(random.uniform(*HOLDING_RATE_PER_YR), 2)
        c0 = round(random.uniform(max(cw["unit_cost_usd"][0], 8.0),
                                  cw["unit_cost_usd"][1]), 2)
        # Small commercial discount steps (2-3.5% per break, inside
        # DISCOUNT_STEP_FRACTION): the price saving Delta_c*D and the
        # extra holding h*b/2 are then the SAME order near realistic
        # breakpoints, so the cost comparison is a genuine decision
        # (with the 12-20% steps of the first design, the deepest
        # discount won every comparison - QA sweep, cycle 1 authoring).
        s1 = random.uniform(0.020, 0.030)
        s2 = random.uniform(0.020, 0.035)
        c1 = round(c0 * (1 - s1), 2)
        c2 = round(c0 * (1 - s1 - s2), 2)
        if not (c2 <= c1 - 0.01 and c1 <= c0 - 0.02):
            continue
        h0_e, h1_e, h2_e = i * c0, i * c1, i * c2
        # Cap D per-sample so the top-tier break-even quantity
        # BE2 = 2*(c1-c2)*D/h2 can sit BELOW the 5000-unit breakpoint
        # ceiling (b2 windows must straddle BE2 for both winners to be
        # reachable): D <= 4200*h2/(2*(c1-c2)).
        d_cap = int(4200 * h2_e / (2 * (c1 - c2)))
        d_lo, d_hi = cw["annual_demand"]
        if d_cap < d_lo:
            continue
        D = random.randint(d_lo, min(d_hi, d_cap))
        h0 = round(h0_e, 4)
        h1 = round(h1_e, 4)
        h2 = round(h2_e, 4)
        k_lo = max(cw["order_cost_usd"][0], math.ceil(D * h0 / 21632))
        k_hi = min(cw["order_cost_usd"][1], math.floor(D * h0 / 8))
        if k_lo > k_hi:
            continue
        K = random.randint(k_lo, k_hi)

        Q0e = math.sqrt(2 * K * D / h0)
        Q1e = math.sqrt(2 * K * D / h1)
        Q2e = math.sqrt(2 * K * D / h2)

        def _b_star(G_ref, cj, hj):
            """Upper break-even b solving cj*D + K*D/b + hj*b/2 = G_ref
            (beyond it the discounted breakpoint LOSES to G_ref); None if
            the tier's cost curve never rises back to G_ref."""
            dG = G_ref - cj * D
            disc = dG * dG - 2 * hj * K * D
            if disc <= 0:
                return None
            return (dG + math.sqrt(disc)) / hj

        branch = random.choice(["a", "b", "c", "c"])
        if branch == "a":
            b2 = _round50(Q2e * random.uniform(0.45, 0.80))
            b1 = _round50(b2 * random.uniform(0.35, 0.65))
        elif branch == "b":
            # sample b2 around the EXACT break-even vs the realizable
            # mid-tier EOQ cost, so both comparison winners are reachable
            G1 = c1 * D + math.sqrt(2 * K * D * h1)
            bs2 = _b_star(G1, c2, h2)
            if bs2 is None:
                continue
            b2 = _round50(bs2 * random.uniform(0.60, 1.45))
            b1 = _round50(Q1e * random.uniform(0.40, 0.78))
            if b2 < Q2e * 1.05:
                continue
        else:
            # 3-way comparison vs the base-tier EOQ cost: sample each
            # breakpoint around its own exact break-even
            G0 = c0 * D + math.sqrt(2 * K * D * h0)
            bs1 = _b_star(G0, c1, h1)
            bs2 = _b_star(G0, c2, h2)
            if bs1 is None or bs2 is None:
                continue
            b1 = _round50(max(Q1e * 1.05, bs1 * random.uniform(0.60, 1.45)))
            b2 = _round50(max(b1 * 1.25, bs2 * random.uniform(0.60, 1.45)))
        if not (100 <= b1 < b2 <= 5000):
            continue
        if any(abs(Qe - b) < 0.04 * b for Qe in (Q0e, Q1e, Q2e)
               for b in (b1, b2)):
            continue
        # branch class must hold after breakpoint rounding
        if branch == "a" and not (Q2e >= b2):
            continue
        if branch == "b" and not (Q2e < b2 and b1 <= Q1e):
            continue
        if branch == "c" and not (Q2e < b2 and Q1e < b1 and Q0e < b1):
            continue

        # interior candidate: whole-unit display-chain safety
        Qje = {"a": Q2e, "b": Q1e, "c": Q0e}[branch]
        if abs((Qje % 1) - 0.5) < 0.06:
            continue

        def _G(q, cj, hj):
            return cj * D + K * D / q + hj * q / 2.0

        q_int = _hu(round(Qje, 1), 0)
        if branch == "a":
            cands = [(_G(q_int, c2, h2), q_int, "the realizable lowest-price EOQ")]
        elif branch == "b":
            cands = [(_G(q_int, c1, h1), q_int, "the realizable mid-price EOQ"),
                     (_G(b2, c2, h2), b2, "the upper breakpoint")]
        else:
            cands = [(_G(q_int, c0, h0), q_int, "the realizable base-price EOQ"),
                     (_G(b1, c1, h1), b1, "the lower breakpoint"),
                     (_G(b2, c2, h2), b2, "the upper breakpoint")]
        cands_sorted = sorted(cands, key=lambda t: t[0])
        if len(cands_sorted) > 1 and \
                (cands_sorted[1][0] - cands_sorted[0][0]) < 0.003 * cands_sorted[0][0]:
            continue
        break
    else:
        raise AssertionError("resample loop exhausted")

    Q_opt = cands_sorted[0][1]
    assert 95 <= Q_opt <= 11150, f"Q_opt out of bounds: {Q_opt}"

    item_phrase = ("an electronic component" if key == "electronic component"
                   else "a packaged food case SKU")
    schedule = (
        f"1 to {b1 - 1} units at ${c0:.2f} per unit; {b1} to {b2 - 1} "
        f"units at ${c1:.2f} per unit; {b2} units or more at ${c2:.2f} "
        f"per unit"
    )

    question = (
        f"A purchasing manager buys {item_phrase} with steady annual "
        f"demand of {D} units. Each order costs ${K} to place, and the "
        f"annual holding-cost rate is {i:.2f} (a fraction of unit value "
        f"per year). The supplier offers an all-units quantity discount "
        f"— the quoted price applies to every unit in the order: "
        f"{schedule}. Using the standard all-units procedure (compute "
        f"the EOQ at each price, check which are realizable within their "
        f"own tier, then compare total annual costs — purchase plus "
        f"ordering plus holding — at the realizable EOQ and at any "
        f"higher breakpoints), determine the optimal order quantity, in "
        f"whole units."
    )

    q0d, q1d, q2d = round(Q0e, 1), round(Q1e, 1), round(Q2e, 1)
    steps = [
        f"**Given:**\n"
        f"Annual demand D = {D} units/year; order cost K = ${K}; holding "
        f"rate i = {i:.2f} per year; all-units prices c0 = ${c0:.2f} "
        f"(1 to {b1 - 1}), c1 = ${c1:.2f} ({b1} to {b2 - 1}), "
        f"c2 = ${c2:.2f} ({b2}+).",
        f"**Step 1:** Compute the holding cost per unit per year at each "
        f"tier price.\n"
        f"h0 = i*c0 = {i:.2f}*{c0:.2f} = {h0:.4f}; "
        f"h1 = i*c1 = {i:.2f}*{c1:.2f} = {h1:.4f}; "
        f"h2 = i*c2 = {i:.2f}*{c2:.2f} = {h2:.4f} $/unit/year",
        f"**Step 2:** Compute the EOQ at each tier price.\n"
        f"Q0 = sqrt(2*{K}*{D}/{h0:.4f}) = {q0d:.1f}; "
        f"Q1 = sqrt(2*{K}*{D}/{h1:.4f}) = {q1d:.1f}; "
        f"Q2 = sqrt(2*{K}*{D}/{h2:.4f}) = {q2d:.1f} units",
    ]

    if branch == "a":
        steps.append(
            f"**Step 3:** Check realizability from the lowest price "
            f"upward. Q2 = {q2d:.1f} >= {b2}, so the lowest-price EOQ "
            f"lies inside its own tier: it is REALIZABLE."
        )
        steps.append(
            f"**Step 4:** Conclude. A realizable EOQ at the lowest price "
            f"is optimal outright: no other candidate can beat the "
            f"cheapest price at its own cost-minimizing quantity. "
            f"Rounding to whole units, Q* = {Q_opt} units."
        )
    elif branch == "b":
        g_q1 = round(cands[0][0], 2)
        g_b2 = round(cands[1][0], 2)
        winner = cands_sorted[0]
        steps.append(
            f"**Step 3:** Check realizability from the lowest price "
            f"upward. Q2 = {q2d:.1f} < {b2}, so the lowest-price EOQ is "
            f"NOT realizable. Q1 = {q1d:.1f} lies in [{b1}, {b2 - 1}], "
            f"so the mid-price EOQ IS realizable (rounded to "
            f"{cands[0][1]} units). The candidates are the realizable "
            f"EOQ and the upper breakpoint {b2}."
        )
        steps.append(
            f"**Step 4:** Total annual cost at the realizable EOQ.\n"
            f"G({cands[0][1]}) = {c1:.2f}*{D} + {K}*{D}/{cands[0][1]} + "
            f"{h1:.4f}*{cands[0][1]}/2 = ${g_q1:,.2f}"
        )
        steps.append(
            f"**Step 5:** Total annual cost at the upper breakpoint.\n"
            f"G({b2}) = {c2:.2f}*{D} + {K}*{D}/{b2} + {h2:.4f}*{b2}/2 "
            f"= ${g_b2:,.2f}"
        )
        steps.append(
            f"**Step 6:** Choose the cheaper candidate. "
            f"${min(g_q1, g_b2):,.2f} < ${max(g_q1, g_b2):,.2f}, so "
            f"{winner[2]} wins: Q* = {Q_opt} units."
        )
    else:
        g_q0 = round(cands[0][0], 2)
        g_b1 = round(cands[1][0], 2)
        g_b2 = round(cands[2][0], 2)
        winner = cands_sorted[0]
        steps.append(
            f"**Step 3:** Check realizability from the lowest price "
            f"upward. Q2 = {q2d:.1f} < {b2} and Q1 = {q1d:.1f} < {b1}: "
            f"neither discounted EOQ reaches its own tier. Only the "
            f"base-price EOQ Q0 = {q0d:.1f} (rounded to {cands[0][1]} "
            f"units) is realizable. The candidates are Q0 and both "
            f"breakpoints."
        )
        steps.append(
            f"**Step 4:** Total annual cost at the base-price EOQ.\n"
            f"G({cands[0][1]}) = {c0:.2f}*{D} + {K}*{D}/{cands[0][1]} + "
            f"{h0:.4f}*{cands[0][1]}/2 = ${g_q0:,.2f}"
        )
        steps.append(
            f"**Step 5:** Total annual cost at the lower breakpoint.\n"
            f"G({b1}) = {c1:.2f}*{D} + {K}*{D}/{b1} + {h1:.4f}*{b1}/2 "
            f"= ${g_b1:,.2f}"
        )
        steps.append(
            f"**Step 6:** Total annual cost at the upper breakpoint.\n"
            f"G({b2}) = {c2:.2f}*{D} + {K}*{D}/{b2} + {h2:.4f}*{b2}/2 "
            f"= ${g_b2:,.2f}"
        )
        steps.append(
            f"**Step 7:** Choose the cheapest candidate. "
            f"${cands_sorted[0][0]:,.2f} beats "
            f"${cands_sorted[1][0]:,.2f} and "
            f"${cands_sorted[2][0]:,.2f}, so {winner[2]} wins: "
            f"Q* = {Q_opt} units."
        )

    solution = "\n\n".join(steps) + (
        f"\n\n**Answer:** The optimal order quantity is {Q_opt} units"
    )

    return question, solution


# Template 14 (Intermediate) — Area P1: Deterministic Lot Sizing
# [BRANCHING: lead time vs. order cycle regime]
def template_reorder_point_lead_time():
    """
    Reorder Point with Deterministic Lead Time: tau < T vs. tau > T

    Scenario:
        A warehouse SKU sells at a steady lambda units per week and is
        replenished in fixed lots of Q units, so the order cycle is
        T = Q / lambda weeks. Orders arrive after a deterministic lead
        time tau. The reorder point (on-hand inventory level that
        triggers an order) BRANCHES on the regime:

            tau < T:  R = lambda * tau
            tau > T:  orders are outstanding at all times; reduce the
                      lead time by whole cycles, tau_eff = tau - k*T
                      with k = floor(tau/T), and R = lambda * tau_eff
                      (equivalently R = lambda*tau - k*Q, exactly).

        Requested: the reorder point in whole units.

    Difficulty: Intermediate
    Grounding: Nahmias, Production and Operations Analysis, 7th ed.,
        Sec. 4.5 "Inclusion of Order Lead Time" (text p. 213; visually
        verified in the on-disk copy, including the tau > T treatment).
        Cross-ref H&L 7e Ch. 19 (reorder-point discussion).
    Physical bounds: lambda integer in [30, 300] units/week; lead time
        tau integer in [1, 10] weeks (inside LEAD_TIME_WEEKS); regime
        sampled 50/50. Under-regime: tau in [1, 6] with tau/T screened
        into [0.30, 0.80]; over-regime: k = floor(tau/T) in {1, 2} with
        the fractional part of tau/T screened into [0.12, 0.88] after Q
        rounding (lesson 23 boundary safety). Q is a positive multiple
        of 10. EXACTNESS: R equals lambda*tau (under) or
        lambda*tau - k*Q (over) — pure integer arithmetic, so the
        displayed chain (T at 3 dp) and a full-precision solve give the
        same whole-unit answer: the true R is an integer, and the
        display wobble (<= lambda*k*5e-4 <= 0.3 units) cannot move it
        across a rounding boundary. Analytic envelope: R_under =
        lambda*tau in [30, 1800]; R_over = lambda*(tau - k*T) with
        fractional part in [0.12, 0.88] of T <= 4.6 weeks, so R_over
        <= 300*0.88*4.6 ~ 1214 and R_over >= ~7. Asserts: R in
        [5, 1850]; 0 < R < Q in every instance.

    Returns:
        tuple(str, str): (question, solution)
    """
    for _ in range(200):
        lam = random.randint(30, 300)
        regime = random.choice(["under", "over"])
        if regime == "under":
            tau = random.randint(1, 6)
            u = random.uniform(1.3, 3.0)
            Q = int(round(lam * tau * u / 10.0)) * 10
            if Q <= 0:
                continue
            ratio = tau * lam / Q                 # = tau / T, exact
            if not (0.30 <= ratio <= 0.80):
                continue
            k = 0
        else:
            k = random.choice([1, 2])
            T_target = random.uniform(2.0, 4.5)
            f = random.uniform(0.20, 0.80)
            tau = int(round((k + f) * T_target))
            if not (2 <= tau <= 10):
                continue
            Q = int(round(lam * T_target / 10.0)) * 10
            if Q <= 0:
                continue
            ratio = tau * lam / Q                 # = tau / T, exact
            if math.floor(ratio) != k:
                continue
            if not (0.12 <= ratio - k <= 0.88):
                continue
        break
    else:
        raise AssertionError("resample loop exhausted")

    T = round(Q / lam, 3)
    if regime == "under":
        R = lam * tau                             # exact integer
        assert 0 < R < Q, f"R out of (0, Q): {R} vs {Q}"
    else:
        R = lam * tau - k * Q                     # exact integer
        tau_eff = round(tau - k * T, 3)
        assert 0 < R < Q, f"R out of (0, Q): {R} vs {Q}"
    assert 5 <= R <= 1850, f"R out of bounds: {R}"

    question = (
        f"A distribution warehouse sells a SKU at a steady rate of {lam} "
        f"units per week and replenishes it in fixed lots of {Q} units. "
        f"The supplier's delivery lead time is a constant {tau} weeks. "
        f"Determine the reorder point — the on-hand inventory level at "
        f"which a new order should be placed — in whole units. Note: "
        f"first compare the lead time with the order cycle T = Q/lambda; "
        f"if the lead time exceeds a full cycle, orders are outstanding "
        f"at all times, and the reorder point is computed from the lead "
        f"time reduced by as many WHOLE cycles as fit into it "
        f"(tau_eff = tau - k*T for the largest integer k with "
        f"k*T < tau)."
    )

    if regime == "under":
        solution = (
            f"**Given:**\n"
            f"Demand rate lambda = {lam} units/week; lot size Q = {Q} "
            f"units; lead time tau = {tau} weeks.\n\n"
            f"**Step 1:** Compute the order cycle.\n"
            f"T = Q / lambda = {Q} / {lam} = {T:.3f} weeks\n\n"
            f"**Step 2:** Compare the lead time with the cycle.\n"
            f"tau = {tau} weeks < T = {T:.3f} weeks, so the order placed "
            f"at the reorder point arrives within the same cycle: the "
            f"reorder point is simply the demand during the lead time.\n\n"
            f"**Step 3:** Compute the reorder point.\n"
            f"R = lambda * tau = {lam} * {tau} = {R} units\n\n"
            f"**Step 4:** Consistency check: R = {R} < Q = {Q}, "
            f"consistent with the lead time fitting inside one cycle.\n\n"
            f"**Answer:** The reorder point is {R} units"
        )
    else:
        cyc_word = "cycles fit" if k > 1 else "cycle fits"
        r_disp = round(lam * tau_eff, 1)
        solution = (
            f"**Given:**\n"
            f"Demand rate lambda = {lam} units/week; lot size Q = {Q} "
            f"units; lead time tau = {tau} weeks.\n\n"
            f"**Step 1:** Compute the order cycle.\n"
            f"T = Q / lambda = {Q} / {lam} = {T:.3f} weeks\n\n"
            f"**Step 2:** Compare the lead time with the cycle.\n"
            f"tau = {tau} weeks > T = {T:.3f} weeks: the lead time "
            f"spans more than one full cycle, so orders are outstanding "
            f"at all times and the lead time must be reduced by whole "
            f"cycles.\n\n"
            f"**Step 3:** Reduce the lead time by whole cycles.\n"
            f"tau / T = {tau} / {T:.3f} = {round(tau / T, 3):.3f}, so "
            f"k = {k} whole {cyc_word} into the lead time.\n"
            f"tau_eff = tau - k*T = {tau} - {k} * {T:.3f} = "
            f"{tau_eff:.3f} weeks\n\n"
            f"**Step 4:** Compute the reorder point from the effective "
            f"lead time.\n"
            f"R = lambda * tau_eff = {lam} * {tau_eff:.3f} = "
            f"{r_disp:.1f}, i.e. R = {R} units — exactly, since "
            f"R = lambda*tau - k*Q = {lam}*{tau} - {k}*{Q} = {R}.\n\n"
            f"**Step 5:** Consistency check: 0 < R = {R} < Q = {Q}, as "
            f"required for a reorder point applied to on-hand "
            f"inventory.\n\n"
            f"**Answer:** The reorder point is {R} units"
        )

    return question, solution
