import math
import random
from decimal import Decimal, ROUND_HALF_UP

from pilot.templates.branches.industrial_engineering.constants import (
    NEWSVENDOR_ITEMS,
    SERVICE_LEVELS,
    Z_QUANTILES,
)


def _hu(x, places):
    """Half-up rounding of a float via its shortest decimal repr."""
    q = Decimal("1") if places == 0 else Decimal("0." + "0" * (places - 1) + "1")
    v = Decimal(repr(x)).quantize(q, rounding=ROUND_HALF_UP)
    return int(v) if places == 0 else float(v)


# Template 15 (Easy) — Area P2: Stochastic Inventory
def template_safety_stock_reorder_point():
    """
    Safety Stock and Reorder Point at a Cycle-Service Level

    Scenario:
        Demand during the replenishment lead time is approximately
        normal with mean mu_L and standard deviation sigma_L. To hit a
        Type-1 (cycle) service level alpha — the probability of not
        stocking out during a lead time — the reorder point is

            SS = z_alpha * sigma_L        (safety stock)
            R  = mu_L + SS

        The z-quantile for the target level is given in the question
        (four decimals, from the branch constants), so the chain is
        reproducible without table interpolation. Requested: R in whole
        units.

    Difficulty: Easy
    Grounding: Nahmias, Production and Operations Analysis, 7th ed.,
        Sec. 5.5 "Service Levels in (Q,R) Systems" (Type-1 service,
        F(R) = alpha; visually verified in the on-disk copy, text
        p. 274). Cross-ref NIST/SEMATECH e-Handbook (normal quantiles).
    Physical bounds: mu_L integer in [50, 2000] units; sigma_L integer
        in the per-sample window [max(5, ceil(0.05*mu_L)),
        floor(0.25*mu_L)] (coefficient of variation 5-25%, always
        non-empty since floor(0.25*50) = 12 >= max(5, 3)); alpha from
        constants.py SERVICE_LEVELS with z from Z_QUANTILES (4 dp,
        stated verbatim in the question). EXACTNESS: z*sigma_L is exact
        at 4 dp (4-dp z times integer sigma); SS is displayed at 2 dp
        (rounding <= 0.005) and draws whose exact R = mu_L + z*sigma_L
        sits within 0.06 of a .5 integer boundary are resampled, so the
        displayed chain and a full-precision solve round to the SAME
        whole unit. Analytic corners: min R = 50 + 1.2816*5 = 56.4 ->
        56; max R = 2000 + 2.3263*500 = 3163.2 -> 3163. Asserts:
        R in [50, 3200]; SS in [6, 1170] (min 1.2816*5 = 6.4, max
        2.3263*500 = 1163.2).

    Returns:
        tuple(str, str): (question, solution)
    """
    for _ in range(100):
        mu = random.randint(50, 2000)
        s_lo = max(5, math.ceil(0.05 * mu))
        s_hi = math.floor(0.25 * mu)
        sigma = random.randint(s_lo, s_hi)
        alpha = random.choice(SERVICE_LEVELS)
        z = Z_QUANTILES[alpha]
        r_exact = mu + z * sigma
        if abs((r_exact % 1) - 0.5) > 0.06:
            break
    else:
        raise AssertionError("resample loop exhausted")

    SS = round(z * sigma, 2)               # display; exact value has 4 dp
    R = _hu(mu + SS, 0)

    assert 6 <= SS <= 1170, f"SS out of bounds: {SS}"
    assert 50 <= R <= 3200, f"R out of bounds: {R}"
    assert R == _hu(r_exact, 0), f"display-chain flip: {SS} vs {r_exact}"

    pct = f"{alpha * 100:g}%"
    question = (
        f"A regional warehouse reorders a product whenever its inventory "
        f"position falls to the reorder point. Demand during the "
        f"replenishment lead time is approximately normally distributed "
        f"with mean {mu} units and standard deviation {sigma} units. "
        f"Management targets a {pct} cycle-service level — the "
        f"probability of not stocking out during a replenishment lead "
        f"time. Using the standard normal quantile z = {z:.4f} for the "
        f"{pct} level, determine the reorder point, in whole units "
        f"(round half up). In your solution, state the safety stock."
    )

    solution = (
        f"**Given:**\n"
        f"Lead-time demand ~ Normal(mu_L = {mu}, sigma_L = {sigma}) "
        f"units; target cycle-service level alpha = {pct}; "
        f"z_alpha = {z:.4f}.\n\n"
        f"**Step 1:** Identify the required quantile. A {pct} "
        f"cycle-service level requires P(demand <= R) = {alpha:g}, so R "
        f"sits z = {z:.4f} standard deviations above the mean.\n\n"
        f"**Step 2:** Compute the safety stock.\n"
        f"SS = z * sigma_L = {z:.4f} * {sigma} = {SS:.2f} units\n\n"
        f"**Step 3:** Compute the reorder point.\n"
        f"R = mu_L + SS = {mu} + {SS:.2f} = {mu + SS:.2f}, which rounds "
        f"half-up to R = {R} units.\n\n"
        f"**Answer:** The reorder point is {R} units"
    )

    return question, solution


# Newsvendor framings: perishable/seasonal classes from constants.py
# NEWSVENDOR_ITEMS with per-class demand scales (units per selling period).
_T16_SETTINGS = {
    "bakery batch": {
        "phrase": "a morning pastry batch at a bakery",
        "period": "day", "mu": (40, 400),
    },
    "seasonal apparel lot": {
        "phrase": "a one-season fashion jacket",
        "period": "season", "mu": (200, 2000),
    },
    "fresh produce crate": {
        "phrase": "crates of fresh berries at a grocer",
        "period": "day", "mu": (30, 300),
    },
}


# Template 16 (Intermediate) — Area P2: Stochastic Inventory
def template_newsvendor_normal_demand():
    """
    Newsvendor with Normal Demand: Optimal Stocking Quantity

    Scenario:
        A perishable/seasonal item is stocked once for a selling period
        with normally distributed demand. Unsold units are salvaged;
        unmet demand is lost margin. With unit cost c, selling price p,
        and salvage value s:

            cu = p - c   (underage cost: margin lost per unsold demand)
            co = c - s   (overage cost: loss per leftover unit)
            CR = cu / (cu + co)          (critical ratio)
            Q* = mu + z_CR * sigma       (normal-demand newsvendor)

        The question PRESCRIBES the numerical scheme (lessons 29/62):
        critical ratio rounded to two decimals; the standard normal
        quantile for that rounded ratio to four decimals; Q* to the
        nearest whole unit (half up). Requested: Q*.

    Difficulty: Intermediate
    Grounding: Nahmias, Production and Operations Analysis, 7th ed.,
        Sec. 5.3 "The Newsvendor Model" (critical ratio; visually
        verified in the on-disk copy, text p. 258 — typology only,
        never its numbers). Cross-ref H&L 7e Ch. 19 Sec. 19.6
        (stochastic single-period model with the same optimality
        condition).
    Physical bounds: item class from constants.py NEWSVENDOR_ITEMS
        (three classes with per-class demand scales); prices sampled
        2 dp inside the class windows with per-sample feasibility
        p >= c + 0.50 and s <= c - 0.30 (margins meaningful; windows
        verified non-empty at all class corners since p_hi > c_hi + 0.5
        and s_lo < c_lo - 0.3 fail nowhere — bakery p_hi 15 vs c_hi 6;
        apparel 120 vs 40; produce 45 vs 20). mu integer from the class
        scale; sigma integer with CV in [0.15, 0.35]. Screens: CR_raw
        in [0.55, 0.95] (stocking above the mean but not extreme);
        |CR_raw - nearest .xx5| >= 0.002 (the 2-dp rounding is
        decisive); exact Q = mu + z*sigma at least 0.06 from a .5
        boundary (display chain and full-precision solve agree on the
        whole unit — z*sigma is exact at 4 dp so the only wobble is the
        2-dp sum display, <= 0.005). Author QA 2026-08-07 (20,000-seed
        sweep): Q in [31, 3077], CR in [0.55, 0.95], class mix roughly
        uniform; analytic floor mu_min + z_min*sigma_min ~ 30 + 0.1257*5
        ~ 30.6. Asserts: Q in [28, 3200]; CR in [0.55, 0.95]; mu < Q
        (CR > 0.5 stocks above the mean).

    Returns:
        tuple(str, str): (question, solution)
    """
    from statistics import NormalDist

    for _ in range(200):
        key = random.choice(sorted(_T16_SETTINGS))
        cfg = _T16_SETTINGS[key]
        cw = NEWSVENDOR_ITEMS[key]
        c = round(random.uniform(*cw["cost_usd"]), 2)
        p_lo = max(cw["price_usd"][0], c + 0.50)
        if p_lo > cw["price_usd"][1]:
            continue
        p = round(random.uniform(p_lo, cw["price_usd"][1]), 2)
        s_hi = min(cw["salvage_usd"][1], c - 0.30)
        if s_hi < cw["salvage_usd"][0]:
            continue
        s = round(random.uniform(cw["salvage_usd"][0], s_hi), 2)
        cu = round(p - c, 2)
        co = round(c - s, 2)
        if cu < 0.50 or co < 0.30:
            continue
        cr_raw = cu / (cu + co)
        if not (0.55 <= cr_raw <= 0.95):
            continue
        # 2-dp rounding decisiveness: stay clear of .xx5 boundaries
        if abs(cr_raw * 100 - math.floor(cr_raw * 100) - 0.5) < 0.2:
            continue
        CR = round(cr_raw, 2)
        z = round(NormalDist().inv_cdf(CR), 4)
        mu = random.randint(*cfg["mu"])
        sigma = random.randint(max(5, math.ceil(0.15 * mu)),
                               math.floor(0.35 * mu))
        q_exact = mu + z * sigma
        if abs((q_exact % 1) - 0.5) > 0.06:
            break
    else:
        raise AssertionError("resample loop exhausted")

    zs = round(z * sigma, 2)               # display; exact value has 4 dp
    Q = _hu(mu + zs, 0)

    assert 0.55 <= CR <= 0.95, f"CR out of bounds: {CR}"
    assert 28 <= Q <= 3200, f"Q out of bounds: {Q}"
    assert mu < Q, f"Q not above the mean: {Q} vs {mu}"
    assert Q == _hu(q_exact, 0), f"display-chain flip: {zs} vs {q_exact}"

    question = (
        f"A retailer stocks {cfg['phrase']} once per {cfg['period']}. "
        f"Each unit costs ${c:.2f}, sells for ${p:.2f} during the "
        f"{cfg['period']}, and any unsold unit is salvaged for "
        f"${s:.2f}. Demand for the {cfg['period']} is approximately "
        f"normally distributed with mean {mu} units and standard "
        f"deviation {sigma} units. Using the newsvendor model, determine "
        f"the stocking quantity that maximizes expected profit. Compute "
        f"the underage and overage costs, the critical ratio rounded to "
        f"two decimals, and the corresponding standard normal quantile "
        f"to four decimals; report the stocking quantity to the nearest "
        f"whole unit (round half up)."
    )

    solution = (
        f"**Given:**\n"
        f"Unit cost c = ${c:.2f}; price p = ${p:.2f}; salvage "
        f"s = ${s:.2f}; demand ~ Normal(mu = {mu}, sigma = {sigma}) "
        f"units per {cfg['period']}.\n\n"
        f"**Step 1:** Compute the underage and overage costs.\n"
        f"cu = p - c = {p:.2f} - {c:.2f} = {cu:.2f} $/unit (margin lost "
        f"per unit of unmet demand);  co = c - s = {c:.2f} - {s:.2f} = "
        f"{co:.2f} $/unit (loss per unsold unit)\n\n"
        f"**Step 2:** Compute the critical ratio (two decimals, as "
        f"prescribed).\n"
        f"CR = cu / (cu + co) = {cu:.2f} / ({cu:.2f} + {co:.2f}) = "
        f"{cr_raw:.4f}, rounded to {CR:.2f}\n\n"
        f"**Step 3:** Find the standard normal quantile for the critical "
        f"ratio (four decimals).\n"
        f"z = Phi^-1({CR:.2f}) = {z:.4f}\n\n"
        f"**Step 4:** Compute the optimal stocking quantity.\n"
        f"Q* = mu + z * sigma = {mu} + {z:.4f} * {sigma} = {mu} + "
        f"{zs:.2f} = {mu + zs:.2f}, which rounds half-up to Q* = {Q} "
        f"units.\n\n"
        f"**Answer:** The profit-maximizing stocking quantity is {Q} "
        f"units"
    )

    return question, solution


# Template 17 (Advanced) — Area P2: Stochastic Inventory
def template_qr_policy_one_iteration():
    """
    (Q, R) Policy with Stockout Penalty: One Iteration of the Q-R
    Procedure

    Scenario:
        A continuous-review (Q, R) system faces annual demand lambda,
        order cost K, holding cost h = i*c, a per-unit backorder
        penalty p, and lead-time demand ~ Normal(mu_L, sigma_L). The
        question states the COUPLED OPTIMALITY SYSTEM

            Q = sqrt(2*lambda*(K + p*n(R))/h)
            1 - F(R) = Q*h/(p*lambda)
            n(R) = sigma_L * L(z),  L(z) = phi(z) - z*(1 - Phi(z))

        with every symbol bound and an itemized precision list, and the
        solver carries out ONE iteration starting from the EOQ (n = 0):
        R0 from the second condition, then the updated Q1 from the
        first. Requested: Q1 in whole units.

    Difficulty: Intermediate
        (LABEL RECONCILIATION, cycle 4: blind relabels across review
        cycles ran Intermediate / Advanced / Intermediate; the majority
        judgment — a guided coupled-system substitution chain at the
        top of the Intermediate band (7 steps, regime-free but
        multi-concept) — is adopted, and the domain's second Advanced
        slot moves to the line-balancing template, whose station
        construction is authentic synthesis.)
    Grounding: Nahmias, Production and Operations Analysis, 7th ed.,
        Sec. 5.4 (the (Q,R) model with stockout cost and its iterative
        solution; visually verified in the on-disk copy — typology
        only). The unit normal loss function is given verbatim in the
        question for self-containment. Cross-ref H&L 7e Ch. 19.
    Physical bounds: lambda integer in [2000, 40000]/yr; c (2 dp) in
        [2, 40]; i (2 dp) in [0.15, 0.40] (h = i*c exact at 4 dp);
        K integer via the template_basic_eoq frequency window; penalty
        p (2 dp) in the per-sample window pinning the stockout
        probability r0 = Q0*h/(p*lambda) into [0.02, 0.25],
        intersected with the realism band [0.5*c, 5*c] (lesson 63).
        Lead time sampled 1-4 weeks (the accepted mix skews to ~69%
        lt=4 / ~30% lt=3; lt=1 is unreachable under the materiality
        intersection): mu_L = round(lambda*lt/52); sigma_L integer in
        the CV band [0.10, 0.30] INTERSECTED with the materiality band
        [0.3K/(p*L0), 2K/(p*L0)] so the stockout term is 30-200% of K
        and Q1/Q0 = sqrt(1 + p*n/K) >= 1.14 (the iteration visibly
        matters). EXACTNESS (lessons 51/65/67): every intermediate
        whose operands are exact decimals — r0, L, n, and R0 — is
        rounded half-up in DECIMAL arithmetic (float repr misrounds
        true ties); z, phi, Phi are transcendental (repr-safe); the
        final Q1 is a SINGLE direct half-up rounding of the raw root,
        screened 0.01 clear of the .5 boundary (dominates the 2-dp
        display wobble). Reviewer sweeps (c2/c3, 100k draws each):
        zero half-up violations on all channels after the fixes; gold
        Q1 equals direct half-up of the exact root everywhere.
        Attained behavior: Q1/Q0 in [1.138, 1.289]; r0 in
        [0.020, 0.051] (implied Type-1 service ~94.9-98.0% — the
        p >= 0.5c floor keeps most of the admissible r0 window
        unreached). Asserts (analytic dominance; the multi-screen
        coupling makes the tight envelope intractable in closed form):
        Q1 in [95, 23000]; Q1 >= 1.1*Q0; r0 in [0.019, 0.251].

    Returns:
        tuple(str, str): (question, solution)
    """
    from statistics import NormalDist

    nd = NormalDist()
    for _ in range(300):
        lam = random.randint(2000, 40000)
        c = round(random.uniform(2.0, 40.0), 2)
        i = round(random.uniform(0.15, 0.40), 2)
        h = round(i * c, 4)
        k_lo = max(40, math.ceil(lam * h / 21632))
        k_hi = min(300, math.floor(lam * h / 8))
        if k_lo > k_hi:
            continue
        K = random.randint(k_lo, k_hi)
        Q0 = _hu(math.sqrt(2 * K * lam / h), 1)

        # per-sample p window: r0 = Q0*h/(p*lam) in [0.02, 0.25]
        # => p in [Q0*h/(0.25*lam), Q0*h/(0.02*lam)], cut by [0.5c, 5c]
        p_lo = max(Q0 * h / (0.25 * lam), 0.5 * c)
        p_hi = min(Q0 * h / (0.02 * lam), 5.0 * c)
        if p_lo + 0.02 > p_hi:
            continue
        p = round(random.uniform(p_lo + 0.01, p_hi - 0.01), 2)
        # Decimal-exact: Q0 (1 dp), h (4 dp), p (2 dp) are exact decimals,
        # so the quotient's half-up 5-dp rounding must use Decimal ties
        # (cycle-3 fix: float repr rounds true ties DOWN - R1/R2 c2)
        r0 = float((Decimal(f"{Q0:.1f}") * Decimal(f"{h:.4f}")
                    / (Decimal(f"{p:.2f}") * lam))
                   .quantize(Decimal("0.00001"), rounding=ROUND_HALF_UP))
        if not (0.019 <= r0 <= 0.251):
            continue

        z0 = _hu(nd.inv_cdf(1 - r0), 4)      # transcendental: repr-safe
        phi0 = _hu(math.exp(-z0 * z0 / 2) / math.sqrt(2 * math.pi), 4)
        Phi0 = _hu(nd.cdf(z0), 4)
        # L0's operands are exact 4-dp decimals -> exact 8-dp result with
        # true .xxxx5 ties: compute in Decimal (cycle-3 fix)
        L0 = float((Decimal(f"{phi0:.4f}") - Decimal(f"{z0:.4f}")
                    * (1 - Decimal(f"{Phi0:.4f}")))
                   .quantize(Decimal("0.0001"), rounding=ROUND_HALF_UP))
        if L0 < 0.005:
            continue
        lt = random.randint(1, 4)
        mu = int(round(lam * lt / 52))
        # sigma window = CV window [0.10*mu, 0.30*mu] INTERSECTED with the
        # materiality window [0.3K/(p*L0), 2K/(p*L0)] (the stockout term
        # p*n = p*sigma*L0 must land between 30% and 200% of K, so the
        # iteration visibly changes the lot: Q1/Q0 = sqrt(1 + p*n/K) in
        # ~[1.14, 1.73]); empty intersections resample.
        s_lo = max(10, math.ceil(0.10 * mu), math.ceil(0.30 * K / (p * L0)))
        s_hi = min(math.floor(0.30 * mu), math.floor(2.0 * K / (p * L0)))
        if s_lo > s_hi:
            continue
        sigma = random.randint(s_lo, s_hi)
        # sigma (int) * L0 (exact 4 dp) is exact at 4 dp with true .xx5
        # ties: Decimal half-up (cycle-3 fix; float repr rounded 798/100k
        # true ties DOWN against the prescription in the c2 sweeps)
        n0 = float((sigma * Decimal(f"{L0:.4f}"))
                   .quantize(Decimal("0.01"), rounding=ROUND_HALF_UP))
        if not (0.29 * K <= p * n0 <= 2.05 * K):
            continue
        q1_raw = math.sqrt(2 * lam * (K + p * n0) / h)
        if abs((q1_raw % 1) - 0.5) < 0.01:
            continue
        break
    else:
        raise AssertionError("resample loop exhausted")

    # mu + z0*sigma is an exact 4-dp decimal: Decimal half-up (c3 fix -
    # the same tie class was fixed for r0/L0/n0 in c2 but R0 was missed)
    R0 = float((mu + Decimal(f"{z0:.4f}") * sigma)
               .quantize(Decimal("0.1"), rounding=ROUND_HALF_UP))
    Q1 = _hu(q1_raw, 0)          # SINGLE direct half-up rounding (R1/R2 c1)

    assert 0.019 <= r0 <= 0.251, f"r0 out of bounds: {r0}"
    assert 95 <= Q1 <= 23000, f"Q1 out of bounds: {Q1}"
    assert Q1 >= 1.1 * Q0, f"Q1 not materially above Q0: {Q1} vs {Q0}"

    question = (
        f"A distributor manages an item on a continuous-review (Q, R) "
        f"policy. Annual demand is lambda = {lam} units, each order "
        f"costs K = ${K} to place, the item costs c = ${c:.2f} per "
        f"unit, and the annual holding-cost rate is i = {i:.2f} (a "
        f"fraction of unit value per year), so h = i*c. Each unit of "
        f"unmet demand is backordered at a stockout penalty of "
        f"p = ${p:.2f}. Demand during the {lt}-week resupply lead time "
        f"is approximately normal with mean mu = {mu} units and "
        f"standard deviation sigma = {sigma} units. The optimal (Q, R) "
        f"pair satisfies the coupled optimality conditions "
        f"Q = sqrt(2*lambda*(K + p*n(R))/h) and "
        f"1 - F(R) = Q*h/(p*lambda), where n(R) = sigma*L(z) is the "
        f"expected shortage per cycle, z is the standard normal score "
        f"of R, and L(z) = phi(z) - z*(1 - Phi(z)) is the unit normal "
        f"loss function, F being the CDF of lead-time demand and phi and "
        f"Phi the standard normal density and CDF. Starting "
        f"from the EOQ as Q0 (which sets n = 0), carry out ONE "
        f"iteration of the standard procedure: determine R0 from the "
        f"second condition, then the updated lot size Q1 from the "
        f"first. Use these precisions, rounding every value half up: "
        f"(a) h to four decimals; (b) Q0 to one decimal; (c) the "
        f"stockout probability to five decimals; (d) z, phi(z), "
        f"Phi(z), and L(z) each to four decimals (use exact "
        f"calculator/software values of phi and Phi, then round); "
        f"(e) R0 to one decimal; (f) n to two decimals. Carry each "
        f"rounded value into all subsequent steps. Report Q1 to the "
        f"nearest whole unit (round half up)."
    )

    solution = (
        f"**Given:**\n"
        f"lambda = {lam} units/yr; K = ${K}; c = ${c:.2f}; i = {i:.2f}; "
        f"stockout penalty p = ${p:.2f}/unit; lead-time demand ~ "
        f"Normal(mu = {mu}, sigma = {sigma}).\n\n"
        f"**Step 1:** Holding cost.\n"
        f"h = i*c = {i:.2f} * {c:.2f} = {h:.4f} $/unit/year\n\n"
        f"**Step 2:** EOQ seed.\n"
        f"Q0 = sqrt(2*K*lambda/h) = sqrt(2 * {K} * {lam} / {h:.4f}) = "
        f"{Q0:.1f} units\n\n"
        f"**Step 3:** Allowed stockout probability. Balancing holding "
        f"against penalty cost, the optimality condition sets\n"
        f"1 - F(R0) = Q0*h/(p*lambda) = {Q0:.1f} * {h:.4f} / "
        f"({p:.2f} * {lam}) = {r0:.5f}\n\n"
        f"**Step 4:** Quantile and reorder point.\n"
        f"z0 = Phi^-1(1 - {r0:.5f}) = {z0:.4f};  "
        f"R0 = mu + z0*sigma = {mu} + {z0:.4f} * {sigma} = {R0:.1f} "
        f"units\n\n"
        f"**Step 5:** Unit normal loss function at z0.\n"
        f"phi({z0:.4f}) = {phi0:.4f}; Phi({z0:.4f}) = {Phi0:.4f};  "
        f"L(z0) = {phi0:.4f} - {z0:.4f} * (1 - {Phi0:.4f}) = {L0:.4f}\n\n"
        f"**Step 6:** Expected shortage per cycle.\n"
        f"n(R0) = sigma * L(z0) = {sigma} * {L0:.4f} = {n0:.2f} units\n\n"
        f"**Step 7:** Updated lot size. Stockouts act like an extra "
        f"ordering cost p*n per cycle, so\n"
        f"Q1 = sqrt(2*lambda*(K + p*n)/h) = sqrt(2 * {lam} * ({K} + "
        f"{p:.2f} * {n0:.2f}) / {h:.4f}) = {q1_raw:.2f}, which rounds "
        f"half-up to Q1 = {Q1} units. As expected, Q1 > Q0 = "
        f"{Q0:.1f}: the stockout penalty inflates the effective fixed "
        f"cost per cycle, so the system orders in larger lots.\n\n"
        f"**Answer:** The updated lot size after one iteration is {Q1} "
        f"units"
    )

    return question, solution
