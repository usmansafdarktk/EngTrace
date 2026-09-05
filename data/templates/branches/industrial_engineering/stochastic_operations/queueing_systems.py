import math
import random

from pilot.templates.branches.industrial_engineering.constants import (
    QUEUE_SCENARIOS,
)


def _mmc_p0(a, rho, c):
    """Erlang P0 for an M/M/c queue from offered load a = lam/mu and
    rho = a/c: P0 = 1 / [sum_{n=0}^{c-1} a^n/n! + a^c/(c!(1-rho))]."""
    s = sum(a ** n / math.factorial(n) for n in range(c))
    return 1.0 / (s + a ** c / (math.factorial(c) * (1.0 - rho)))

# Single-server framings for the M/M/1 template: display phrase per
# constants.py QUEUE_SCENARIOS key (Civil lesson 19: proofread rendered
# prose, not just f-strings). Only scenarios whose windows admit integer
# (lam, mu) with rho in [0.55, 0.92] for EVERY integer mu in the window
# (verified analytically; see sampling note below).
_MM1_SETTINGS = {
    "bank teller line": "a bank branch with a single open teller window",
    "drive-through window": "a coffee shop drive-through with one service window",
    "tool crib counter": "a factory tool crib with one attendant",
}


# Template 1 (Easy) — Area S1: Queueing Systems
def template_mm1_time_in_system():
    """
    M/M/1 Queue: Average Time in System

    Scenario:
        Customers arrive at a single-server service point according to a
        Poisson process (rate lambda per hour) and service times are
        exponential (rate mu per hour) — an M/M/1 queue. The steady-state
        results

            rho = lambda / mu          (must satisfy rho < 1)
            L   = lambda / (mu - lambda)
            W   = L / lambda           (Little's formula)

        give the long-run average number in the system and the average
        time a customer spends in the system (waiting plus service),
        which is requested in minutes.

    Difficulty: Easy
    Grounding: Hillier & Lieberman, Introduction to Operations Research,
        7th ed., Ch. 17 — Sec. 17.2 (Little's formula L = lambda*W) and
        Sec. 17.6 (birth-and-death queueing models; M/M/1 results).
        Cross-ref Ross, Introduction to Probability Models, 11th ed.,
        Ch. 8.
    Physical bounds: integer rates with lambda, mu inside the scenario's
        constants.py windows; utilization rho = lambda/mu in [0.55, 0.92];
        L in [1.0, 12.0] customers; answer W in [1.5, 61.0] minutes
        (analytic corners: min 60/(70-39) = 1.9 at the drive-through
        window extreme; max 60/(6-5) = 60.0 at the tool-crib extreme).

    Returns:
        tuple(str, str): (question, solution)
    """
    scenario_key = random.choice(sorted(_MM1_SETTINGS))
    setting = _MM1_SETTINGS[scenario_key]
    lam_lo, lam_hi = QUEUE_SCENARIOS[scenario_key]["lam_hr"]
    mu_lo, mu_hi = QUEUE_SCENARIOS[scenario_key]["mu_hr"]

    # Sample the service rate first, then the arrival rate inside the
    # per-sample window that guarantees 0.55 <= rho <= 0.92 exactly
    # (Civil lesson 1: per-sample joint feasibility, no rejection loops).
    # Non-emptiness holds for every integer mu in each scenario window:
    #   bank (12..50), drive-through (20..70), tool crib (6..24) all give
    #   ceil(0.55*mu) <= min(lam_hi, floor(0.92*mu)) — checked by hand.
    mu = random.randint(mu_lo, mu_hi)
    lam_min = max(lam_lo, math.ceil(0.55 * mu))
    lam_max = min(lam_hi, math.floor(0.92 * mu))
    lam = random.randint(lam_min, lam_max)

    # Gold trace derives only from the presented integer rates
    # (round-then-recompute).
    rho = round(lam / mu, 3)
    L = round(lam / (mu - lam), 3)
    W_hr = round(L / lam, 4)
    W_min = round(W_hr * 60, 1)

    # Physical bounds (docstring, verbatim)
    assert 0.55 <= lam / mu <= 0.92, f"utilization out of bounds: {lam}/{mu}"
    assert 1.0 <= L <= 12.0, f"L out of bounds: {L}"
    assert 1.5 <= W_min <= 61.0, f"W out of bounds: {W_min} min"

    question = (
        f"Customers arrive at {setting} according to a Poisson process at "
        f"an average rate of {lam} customers per hour. Service times are "
        f"exponentially distributed, and the single server completes "
        f"services at an average rate of {mu} customers per hour. Treating "
        f"the operation as an M/M/1 queue in steady state, determine the "
        f"average total time a customer spends in the system (waiting plus "
        f"service), in minutes. In your solution, verify that the system "
        f"reaches steady state and compute the average number of customers "
        f"in the system."
    )

    solution = (
        f"**Given:**\n"
        f"Arrival rate (lambda): {lam} customers/hour; service rate (mu): "
        f"{mu} customers/hour; single server (M/M/1).\n\n"
        f"**Step 1:** Verify the steady-state condition via the "
        f"utilization factor.\n"
        f"rho = lambda / mu = {lam} / {mu} = {rho:.3f}\n"
        f"Since rho = {rho:.3f} < 1, the queue is stable and steady-state "
        f"results apply.\n\n"
        f"**Step 2:** Compute the average number of customers in the "
        f"system.\n"
        f"L = lambda / (mu - lambda) = {lam} / ({mu} - {lam}) "
        f"= {L:.3f} customers\n\n"
        f"**Step 3:** Apply Little's formula to get the average time in "
        f"the system.\n"
        f"W = L / lambda = {L:.3f} / {lam} = {W_hr:.4f} hours\n\n"
        f"**Step 4:** Convert the time in system to minutes.\n"
        f"W = {W_hr:.4f} hours * 60 minutes/hour = {W_min:.1f} minutes\n\n"
        f"**Answer:** The average time a customer spends in the system is "
        f"{W_min:.1f} minutes"
    )

    return question, solution


# Multi-server framings for the M/M/c template (time-unit-mixing template of
# this area per AUTHOR_NOTES standing conventions: lambda per hour, mean
# service TIME in minutes).
_MMC_SETTINGS = {
    "call center": "a customer-support call center",
    "bank teller line": "a bank branch during the midday peak",
    "hospital emergency room": "a walk-in urgent-care clinic",
}

# Mean service times (minutes) per scenario, chosen from divisors of 60 so
# mu = 60/ts is an exact integer per hour inside the scenario's mu window.
_MMC_TS_MIN = {
    "call center": (2, 3, 4, 5, 6),
    "bank teller line": (2, 3, 4, 5),
    "hospital emergency room": (10, 12, 15, 20, 30),
}

# Feasible (scenario, c, ts) combos: c servers in the scenario window
# (capped at 3 so the P0 sum stays printable), and a non-empty integer
# lambda window enforcing rho = lam/(c*mu) in [0.60, 0.82] (AUTHOR_NOTES
# lesson 33 realism/precision cap). The lambda floor is then advanced to
# the first integer whose EXACT Wq is >= 2.0 minutes (Wq is monotone
# increasing in lambda at fixed c, mu), so the 2-dp minute display never
# quantizes worse than 0.005/2.0 = 0.25% (lessons 5/30). Built
# deterministically at import; per-sample windows are non-empty by
# construction.
_MMC_COMBOS = []
for _key in sorted(_MMC_SETTINGS):
    _lam_lo, _lam_hi = QUEUE_SCENARIOS[_key]["lam_hr"]
    _c_lo, _c_hi = QUEUE_SCENARIOS[_key]["servers"]
    for _c in range(max(2, _c_lo), min(3, _c_hi) + 1):
        for _ts in _MMC_TS_MIN[_key]:
            _mu = 60 // _ts
            _lo = max(_lam_lo, math.ceil(0.60 * _c * _mu))
            _hi = min(_lam_hi, math.floor(0.82 * _c * _mu))
            while _lo <= _hi:
                _ae, _rhoe = _lo / _mu, _lo / (_c * _mu)
                _Lqe = (_mmc_p0(_ae, _rhoe, _c) * _ae ** _c * _rhoe
                        / (math.factorial(_c) * (1.0 - _rhoe) ** 2))
                if _Lqe / _lo * 60 >= 2.0:
                    break
                _lo += 1
            if _lo <= _hi:
                _MMC_COMBOS.append((_key, _c, _ts, _lo, _hi))


# Template 2 (Intermediate) — Area S1: Queueing Systems
def template_mmc_waiting_time():
    """
    M/M/c Queue: Average Waiting Time in Queue

    Scenario:
        Customers arrive at a service facility with c identical parallel
        servers according to a Poisson process (lambda per hour); the mean
        service time is given in MINUTES, so the service rate must first
        be converted to per-hour units. The M/M/c steady-state chain is

            mu  = 60 / ts                 (per hour, ts in minutes)
            a   = lambda / mu             (offered load, Erlangs)
            rho = lambda / (c * mu) = a/c (must satisfy rho < 1)
            P0  = 1 / [ sum_{n=0}^{c-1} a^n/n! + a^c / (c! (1 - rho)) ]
            Lq  = P0 * a^c * rho / (c! * (1 - rho)^2)
            Wq  = Lq / lambda             (Little's formula on the queue)

        and the requested quantity is Wq in minutes.

    Difficulty: Intermediate
    Grounding: Hillier & Lieberman, Introduction to Operations Research,
        7th ed., Ch. 17, Sec. 17.6 (M/M/s model: P0, Lq formulas) and
        Sec. 17.2 (Little's formula). Cross-ref Ross 11e Ch. 8.
    Physical bounds: (scenario, c, ts) drawn from the precomputed feasible
        set; c in {2, 3}; rho in [0.60, 0.82] with the lambda floor
        advanced so exact Wq >= 2.0 min. Exhaustive enumeration of all 82
        reachable (combo, lambda) instances (author QA, 2026-08-05) gives
        rho in [0.6000, 0.8167], P0 in [0.05045, 0.25000], Lq in
        [0.6750, 3.2717] customers, answer Wq in [2.03, 38.57] minutes;
        asserts use P0 [0.045, 0.26], Lq [0.60, 3.40], Wq [1.9, 40.0]
        with margin for the rounding chain.

    Returns:
        tuple(str, str): (question, solution)
    """
    key, c, ts, lam_lo, lam_hi = random.choice(_MMC_COMBOS)
    setting = _MMC_SETTINGS[key]
    lam = random.randint(lam_lo, lam_hi)

    # Round-then-recompute: the gold chain derives only from the presented
    # (lam, ts, c); display precisions are sized to the (1-rho)^-2
    # amplification (lessons 5/30/33): a, rho at 4 dp; P0 at 5 dp; Lq 4 dp.
    mu = 60 // ts                      # exact integer by construction
    a = round(lam / mu, 4)
    rho = round(lam / (c * mu), 4)
    P0 = round(_mmc_p0(a, rho, c), 5)
    Lq = round(P0 * a ** c * rho / (math.factorial(c) * (1.0 - rho) ** 2), 4)
    Wq_hr = round(Lq / lam, 5)
    Wq_min = round(Wq_hr * 60, 2)

    # Physical bounds (docstring, verbatim; margins cover displayed rounding)
    assert 0.60 <= lam / (c * mu) <= 0.82, f"rho out of bounds: {lam}/({c}*{mu})"
    assert 0.045 <= P0 <= 0.26, f"P0 out of bounds: {P0}"
    assert 0.60 <= Lq <= 3.40, f"Lq out of bounds: {Lq}"
    assert 1.9 <= Wq_min <= 40.0, f"Wq out of bounds: {Wq_min} min"

    if c == 2:
        p0_eq = (
            f"P0 = 1 / (1 + a + a^2/(2*(1-rho))) "
            f"= 1 / (1 + {a:.4f} + ({a:.4f})^2/(2*(1 - {rho:.4f}))) "
            f"= {P0:.5f}"
        )
        lq_eq = (
            f"Lq = P0 * a^2 * rho / (2! * (1-rho)^2) "
            f"= {P0:.5f} * ({a:.4f})^2 * {rho:.4f} / (2 * (1 - {rho:.4f})^2) "
            f"= {Lq:.4f} customers"
        )
    else:
        p0_eq = (
            f"P0 = 1 / (1 + a + a^2/2 + a^3/(6*(1-rho))) "
            f"= 1 / (1 + {a:.4f} + ({a:.4f})^2/2 + "
            f"({a:.4f})^3/(6*(1 - {rho:.4f}))) = {P0:.5f}"
        )
        lq_eq = (
            f"Lq = P0 * a^3 * rho / (3! * (1-rho)^2) "
            f"= {P0:.5f} * ({a:.4f})^3 * {rho:.4f} / (6 * (1 - {rho:.4f})^2) "
            f"= {Lq:.4f} customers"
        )

    question = (
        f"Customers arrive at {setting} according to a Poisson process at "
        f"an average rate of {lam} customers per hour. The facility has "
        f"{c} identical servers working in parallel, and each service takes "
        f"an exponentially distributed time averaging {ts} minutes. "
        f"Treating the operation as an M/M/{c} queue, determine the average "
        f"time a customer waits in the queue before service begins, in "
        f"minutes. In your solution, first express the service rate in "
        f"customers per hour, verify that a steady state exists, and "
        f"compute the probability that the system is empty and the average "
        f"queue length."
    )

    solution = (
        f"**Given:**\n"
        f"Arrival rate (lambda): {lam} customers/hour; mean service time: "
        f"{ts} minutes per customer; parallel servers (c): {c}.\n\n"
        f"**Step 1:** Convert the mean service time to a service rate in "
        f"per-hour units.\n"
        f"mu = 60 / {ts} = {mu} customers/hour per server\n\n"
        f"**Step 2:** Compute the offered load and the utilization, and "
        f"verify that a steady state exists.\n"
        f"a = lambda / mu = {lam} / {mu} = {a:.4f} (Erlangs)\n"
        f"rho = a / c = {a:.4f} / {c} = {rho:.4f}\n"
        f"Since rho = {rho:.4f} < 1, a steady state exists.\n\n"
        f"**Step 3:** Compute the probability that the system is empty.\n"
        f"{p0_eq}\n\n"
        f"**Step 4:** Compute the average number waiting in the queue.\n"
        f"{lq_eq}\n\n"
        f"**Step 5:** Apply Little's formula to the queue.\n"
        f"Wq = Lq / lambda = {Lq:.4f} / {lam} = {Wq_hr:.5f} hours\n\n"
        f"**Step 6:** Convert the waiting time to minutes.\n"
        f"Wq = {Wq_hr:.5f} hours * 60 minutes/hour = {Wq_min:.2f} minutes\n\n"
        f"**Answer:** The average time a customer waits in the queue is "
        f"{Wq_min:.2f} minutes"
    )

    return question, solution


# Configuration-selection template (BRANCHING): one experienced (fast)
# server vs two standard servers. Anchored to the counter-service
# [REALISM] windows; both options' service rates and the arrival rate
# stay inside them.
#
# REMEDIATION (Stage E escalation #1, 2026-08-15): the grid was
# previously parameterised by mean service TIMES in whole minutes, which
# had to divide 60 to keep the hourly rates integral. Only eight values
# do, and with rho1 and rho2 both confined to [0.55, 0.85] that left just
# 21 reachable triples and 17 distinct answers — so the Stage E pack
# (seeds 201-205) drew only 3 distinct questions and shipped two
# duplicate question PAIRS into the 150-record testset. Parameterising by
# integer hourly RATES removes the divisibility constraint while keeping
# every exactness property (muA - lam is still a positive integer, so
# W1 = 60/(muA - lam) is still exact), and the scenario now carries four
# service contexts. Reachable triples: 21 -> 1845.
#
# Decisiveness rule (lesson 23): only (muA, muB, lam) triples whose EXACT
# mean times in system differ by >= 10% are reachable, so display
# rounding can never flip the winner; both winners occur in the reachable
# set (asserted in the builder tally below).
_SEL_SETTINGS = [
    {"place": "A bank branch", "server": "teller", "arrival": "customers",
     "arrive_at": "a single counter area"},
    {"place": "A pharmacy", "server": "pharmacist", "arrival": "customers",
     "arrive_at": "the prescription counter"},
    {"place": "A hospital admissions desk", "server": "admissions clerk",
     "arrival": "patients", "arrive_at": "the desk"},
    {"place": "A passport office", "server": "processing officer",
     "arrival": "applicants", "arrive_at": "the service window"},
]
_SEL_COMBOS = []
_SEL_WINNER_TALLY = {"single": 0, "pair": 0}
for _muA in range(10, 46):               # fast server, per hour
    for _muB in range(6, 31):            # standard server, per hour
        if _muB >= _muA:                 # "fast" must actually be faster
            continue
        for _lam in range(4, 51):
            _rhoA, _rhoB = _lam / _muA, _lam / (2 * _muB)
            if not (0.55 <= _rhoA <= 0.85 and 0.55 <= _rhoB <= 0.85):
                continue
            _WA = 60.0 / (_muA - _lam)                       # minutes, exact
            _LB = 2 * _rhoB / (1 - _rhoB ** 2)
            _WB = _LB / _lam * 60                            # minutes, exact
            if abs(_WA - _WB) / min(_WA, _WB) < 0.10:
                continue
            if not (1.9 <= _WA <= 41.0 and 1.9 <= _WB <= 41.0):
                continue
            _SEL_COMBOS.append((_muA, _muB, _lam))
            _SEL_WINNER_TALLY["single" if _WA < _WB else "pair"] += 1
assert _SEL_WINNER_TALLY["single"] > 0 and _SEL_WINNER_TALLY["pair"] > 0

# Template 3 (Intermediate) — Area S1: Queueing Systems  [BRANCHING]
def template_server_configuration_selection():
    """
    Service-Configuration Selection: One Fast Server vs. Two Slow Servers

    Scenario:
        A service counter can be staffed either with one experienced
        server (mean rate muA per hour) or with two standard servers in
        parallel (mean rate muB per hour each, muB < muA). Arrivals are
        Poisson at lambda per hour; service times are exponential. The
        two options are compared on the average time in the system:

            Option 1 (M/M/1):
                W1 = 1 / (muA - lambda)
            Option 2 (M/M/2):
                rho2 = lambda / (2*muB);  L2 = 2*rho2 / (1 - rho2^2);
                W2 = L2 / lambda   (Little's formula)

        Which option wins is parameter-dependent (BRANCHING): the
        reachable set contains both winners in near-equal proportion, and
        every reachable triple keeps the exact W1, W2 at least 10% apart
        so the decision is never a rounding artifact.

    Difficulty: Intermediate
    Grounding: Hillier & Lieberman, Introduction to Operations Research,
        7th ed., Ch. 17, Sec. 17.6 (M/M/s results; the one-fast-vs-
        several-slow comparison is a classic Ch. 17/18 decision typology).
        Cross-ref Taha 10e Sec. 18.9 (queueing decision models).
    Physical bounds: service rates are drawn as INTEGERS per hour —
        muA in [10, 45], muB in [6, 30] with muB < muA — and lambda as an
        integer in [4, 50]; the reachable set spans muA giving 1.3-6.0
        min per service and muB giving 2.0-10.0 min, with lambda 7-38/hr.
        Utilizations rho1 = lambda/muA and rho2 = lambda/(2*muB) are both
        held in [0.55, 0.85]; the EXACT |W1 - W2| / min(W1, W2) >= 0.10;
        both W in [1.9, 41.0] minutes, asserted after the draw.
        EXACTNESS: muA - lambda is a positive integer, so
        W1 = 60/(muA - lambda) minutes is exact before its single 2-dp
        display rounding; the Option-2 chain is round-then-recompute from
        the displayed rho2 and L2 (lessons 2/5/24).
        DIVERSITY (Stage E remediation, 2026-08-15): parameterising by
        integer RATES rather than by whole-minute service times — which
        had to divide 60 — took the reachable set from 21 triples and 17
        distinct answers to 1845 triples and 201, and four service
        contexts multiply the distinct question surface. The old ceiling
        put only 3 distinct questions in the 5-seed Stage E pack and
        shipped two duplicate question pairs; see the review log.

    Returns:
        tuple(str, str): (question, solution)
    """
    muA, muB, lam = random.choice(_SEL_COMBOS)
    cfg = random.choice(_SEL_SETTINGS)

    # Round-then-recompute from the presented (muA, muB, lam).
    rho1 = round(lam / muA, 4)
    rho2 = round(lam / (2 * muB), 4)
    W1_min = round(60.0 / (muA - lam), 2)          # exact integer denominator
    L2 = round(2 * rho2 / (1 - rho2 ** 2), 4)
    W2_hr = round(L2 / lam, 5)
    W2_min = round(W2_hr * 60, 2)

    winner_is_single = W1_min < W2_min
    W_best = W1_min if winner_is_single else W2_min

    assert muB < muA, f"fast server must be faster: {muA} vs {muB}"
    assert 0.55 <= lam / muA <= 0.85, f"rho1 out of bounds: {lam}/{muA}"
    assert 0.55 <= lam / (2 * muB) <= 0.85, f"rho2 out of bounds: {lam}/(2*{muB})"
    assert 1.9 <= W1_min <= 41.0, f"W1 out of bounds: {W1_min}"
    assert 1.9 <= W2_min <= 41.0, f"W2 out of bounds: {W2_min}"
    assert W1_min != W2_min, "winner must be decisive at display precision"

    if winner_is_single:
        conclusion = (
            f"Since W1 = {W1_min:.2f} min < W2 = {W2_min:.2f} min, the "
            f"single experienced {cfg['server']} gives the smaller average "
            f"time in the system."
        )
    else:
        conclusion = (
            f"Since W2 = {W2_min:.2f} min < W1 = {W1_min:.2f} min, the two "
            f"standard {cfg['server']}s give the smaller average time in "
            f"the system."
        )

    question = (
        f"{cfg['place']} expects {cfg['arrival']} to arrive at "
        f"{cfg['arrive_at']} according to a Poisson process at {lam} "
        f"{cfg['arrival']} per hour. Management can staff it either with "
        f"one experienced {cfg['server']}, who serves an average of "
        f"{muA} {cfg['arrival']} per hour, or with two standard "
        f"{cfg['server']}s working in parallel, each serving an average "
        f"of {muB} {cfg['arrival']} per hour. Service times are "
        f"exponentially distributed in both cases. Model the first option "
        f"as an M/M/1 queue and the second as an M/M/2 queue, verify that "
        f"a steady state exists for both, and compare the average time a "
        f"{cfg['arrival'][:-1]} spends in the system (waiting plus "
        f"service) under each option. Report, in minutes, the average "
        f"time in the system achieved by the better option."
    )

    solution = (
        f"**Given:**\n"
        f"Arrival rate (lambda): {lam} {cfg['arrival']}/hour; Option 1: "
        f"one server at muA = {muA}/hour; Option 2: two servers at "
        f"muB = {muB}/hour each.\n\n"
        f"**Step 1:** Verify that a steady state exists for both options.\n"
        f"rho1 = lambda / muA = {lam} / {muA} = {rho1:.4f} < 1;  "
        f"rho2 = lambda / (2*muB) = {lam} / {2 * muB} = {rho2:.4f} < 1\n"
        f"Both utilizations are below 1, so both options are stable.\n\n"
        f"**Step 2:** Average time in system for Option 1 (M/M/1).\n"
        f"W1 = 1 / (muA - lambda) = 1 / ({muA} - {lam}) hours "
        f"= 60 / {muA - lam} = {W1_min:.2f} minutes\n\n"
        f"**Step 3:** Average number in system for Option 2 (M/M/2), using "
        f"the standard result L = 2*rho / (1 - rho^2).\n"
        f"L2 = 2 * {rho2:.4f} / (1 - ({rho2:.4f})^2) = {L2:.4f} "
        f"{cfg['arrival']}\n\n"
        f"**Step 4:** Average time in system for Option 2 via Little's "
        f"formula.\n"
        f"W2 = L2 / lambda = {L2:.4f} / {lam} = {W2_hr:.5f} hours "
        f"= {W2_min:.2f} minutes\n\n"
        f"**Step 5:** Select the better configuration.\n"
        f"{conclusion}\n\n"
        f"**Answer:** The better option achieves an average time in the "
        f"system of {W_best:.2f} minutes"
    )

    return question, solution


# M/M/1/K finite-capacity template: the deliberate exception where rho >= 1
# is admissible (BOOKS.md branching plan) — the finite waiting line keeps the
# system stable and the steady-state distribution is DERIVED in the trace
# from the birth-death balance equations (cycle-2 revision: R3 relabeled the
# closed-form-substitution version Intermediate; Advanced is earned by
# construction per lesson 15/41). Anchored to the drive-through [REALISM]
# windows and FINITE_CAPACITY_K.
# Combos built at import: integer (lam, mu) with rho = lam/mu in
# [0.70, 0.92] or [1.08, 1.30] (the gap keeps the two regimes cleanly
# separated and bounds 4-dp-rho error amplification in rho^K; lessons
# 5/17/30), K in [4, 6] (inside FINITE_CAPACITY_K; capped at 6 so the
# derived normalization and expected-value sums stay printable), and only
# combos whose EXACT W is >= 2.0 minutes (2-dp quantization <= 0.25%).
_MM1K_COMBOS = []
for _mu in range(20, 71):                # drive-through mu window
    for _lam in range(15, 56):           # drive-through lam window
        _r = _lam / _mu
        if 0.70 <= _r <= 0.92 or 1.08 <= _r <= 1.30:
            for _K in range(4, 7):
                _S = sum(_r ** _n for _n in range(_K + 1))
                _P0e = 1.0 / _S
                _PKe = _P0e * _r ** _K
                _Le = _P0e * sum(_n * _r ** _n for _n in range(1, _K + 1))
                _We = _Le / (_lam * (1 - _PKe)) * 60
                if _We >= 2.0:
                    _MM1K_COMBOS.append((_lam, _mu, _K))
# Split by regime and sample the regime first (50/50) so overloaded
# (rho >= 1) instances appear with equal frequency — the regime commentary
# branch is this template's pedagogical point (Stage D branch balance).
_MM1K_OVER = [c for c in _MM1K_COMBOS if c[0] > c[1]]
_MM1K_UNDER = [c for c in _MM1K_COMBOS if c[0] < c[1]]
assert _MM1K_OVER and _MM1K_UNDER


# Template 4 (Advanced) — Area S1: Queueing Systems  [rho >= 1 admissible]
def template_mm1k_finite_capacity():
    """
    M/M/1/K Finite-Capacity Queue: Derived Distribution and Time in System

    Scenario:
        A single-server drive-through lane holds at most K cars (including
        the one being served); an arriving car that finds the lane full is
        lost. Arrivals are Poisson (lambda per hour), service exponential
        (rate mu per hour). The solver must CONSTRUCT the steady-state
        distribution from the birth-death balance equations before any
        performance measure can be computed:

            rate up = rate down across each cut:  lambda*p_(n-1) = mu*p_n
            =>  p_n = rho^n * p_0,  n = 0..K,  rho = lambda/mu
            normalization:  p_0 * (1 + rho + ... + rho^K) = 1
            PK = p_0 * rho^K            (blocking probability)
            lam_e = lambda * (1 - PK)   (effective arrival rate)
            L = p_0 * sum_{n=1}^{K} n * rho^n   (expected number, direct)
            W = L / lam_e               (Little's formula, effective rate)

        Because the state space is finite (at most K cars), the chain is
        stable for ANY rho — including rho >= 1 — which the trace reasons
        about explicitly in both regimes. Requested: W in minutes for the
        cars that actually join.

    Difficulty: Advanced
    Grounding: Hillier & Lieberman, Introduction to Operations Research,
        7th ed., Ch. 17, Sec. 17.5 (birth-and-death balance equations) and
        Sec. 17.6 ("The Finite Queue Variation of the M/M/s Model", s = 1);
        Little's formula with the effective arrival rate per Sec. 17.2.
        Cross-ref Ross 11e Ch. 8.
    Physical bounds: integer lam in [15, 55], mu in [20, 70] (drive-
        through windows) with rho = lam/mu in [0.70, 0.92] or
        [1.08, 1.30]; the regime (rho above vs. below 1) is sampled 50/50
        first; K in [4, 6] (inside FINITE_CAPACITY_K); builder keeps only
        combos with exact W >= 2.0 min. Exhaustive enumeration of all
        1949 reachable combos (author QA, 2026-08-06): PK in
        [0.0385, 0.3158], lam_e in [13.4, 51.6]/hr, L in
        [1.3232, 3.9935], answer W in [2.00, 12.70] min; asserts use
        PK [0.037, 0.32], lam_e [13.0, 52.0], L [1.30, 4.05],
        W [1.95, 13.0] with rounding-chain margin.

    Returns:
        tuple(str, str): (question, solution)
    """
    regime = random.choice(["under", "over"])
    lam, mu, K = random.choice(_MM1K_OVER if regime == "over" else _MM1K_UNDER)

    # Round-then-recompute from the presented (lam, mu, K): rho displayed
    # at 4 dp anchors the chain; the normalization sum S and weighted sum
    # SL are computed from displayed rho and themselves displayed at 4 dp.
    rho = round(lam / mu, 4)
    S = round(sum(rho ** n for n in range(K + 1)), 4)
    P0 = round(1.0 / S, 5)
    PK = round(P0 * rho ** K, 5)
    lam_e = round(lam * (1 - PK), 3)
    SL = round(sum(n * rho ** n for n in range(1, K + 1)), 4)
    L = round(P0 * SL, 4)
    W_hr = round(L / lam_e, 5)
    W_min = round(W_hr * 60, 2)

    r_exact = lam / mu
    assert (0.70 <= r_exact <= 0.92) or (1.08 <= r_exact <= 1.30), \
        f"rho out of bounds: {lam}/{mu}"
    assert 4 <= K <= 6, f"K out of bounds: {K}"
    assert 0.037 <= PK <= 0.32, f"PK out of bounds: {PK}"
    assert 13.0 <= lam_e <= 52.0, f"lam_e out of bounds: {lam_e}"
    assert 1.30 <= L <= 4.05, f"L out of bounds: {L}"
    assert 1.95 <= W_min <= 13.0, f"W out of bounds: {W_min} min"

    if r_exact > 1:
        regime_note = (
            f"The state space is finite — the lane never holds more than "
            f"K = {K} cars — so the chain has a proper steady state for "
            f"ANY value of rho. Here rho = {rho:.4f} >= 1, which would "
            f"make an unlimited queue grow without bound, but the finite "
            f"capacity keeps the system stable."
        )
    else:
        regime_note = (
            f"The state space is finite — the lane never holds more than "
            f"K = {K} cars — so the chain has a proper steady state for "
            f"ANY value of rho, even rho >= 1. Here rho = {rho:.4f} < 1, "
            f"but note that stability comes from the finite capacity, not "
            f"from rho being below 1."
        )

    s_terms = " + ".join(["1"] + [f"({rho:.4f})^{n}" if n > 1 else f"{rho:.4f}"
                                  for n in range(1, K + 1)])
    sl_terms = " + ".join([f"{n}*({rho:.4f})^{n}" if n > 1 else f"1*{rho:.4f}"
                           for n in range(1, K + 1)])

    question = (
        f"Cars arrive at a single-window drive-through according to a "
        f"Poisson process at {lam} cars per hour. Service times at the "
        f"window are exponentially distributed, completed at a rate of "
        f"{mu} cars per hour. The lane holds at most {K} cars in total, "
        f"including the car in service; a car that arrives to find the "
        f"lane full drives away and is lost. Starting from the "
        f"birth-and-death balance equations, derive the steady-state "
        f"probability distribution of the number of cars in the lane, "
        f"explain why the system is stable regardless of the ratio "
        f"lambda/mu, and determine the average time in the system "
        f"(waiting plus service), in minutes, experienced by the cars "
        f"that actually join the lane. In your solution, compute the "
        f"blocking probability, the effective arrival rate, and the "
        f"average number of cars in the lane from the derived "
        f"distribution."
    )

    solution = (
        f"**Given:**\n"
        f"Arrival rate (lambda): {lam} cars/hour; service rate (mu): {mu} "
        f"cars/hour; capacity: K = {K} cars (including the one in "
        f"service).\n\n"
        f"**Step 1:** Compute the traffic intensity and assess stability.\n"
        f"rho = lambda / mu = {lam} / {mu} = {rho:.4f}\n"
        f"{regime_note}\n\n"
        f"**Step 2:** Construct the steady-state distribution from the "
        f"balance equations. Across the cut between states n-1 and n "
        f"(for n = 1..{K}), rate up = rate down:\n"
        f"lambda * p_(n-1) = mu * p_n  =>  p_n = rho * p_(n-1)  =>  "
        f"p_n = rho^n * p_0 for n = 0..{K}; states above {K} do not exist "
        f"because arrivals to a full lane are lost.\n\n"
        f"**Step 3:** Normalize the distribution to find p_0.\n"
        f"p_0 * ({s_terms}) = 1\n"
        f"The bracket sums to S = {S:.4f}, so p_0 = 1 / {S:.4f} "
        f"= {P0:.5f}\n\n"
        f"**Step 4:** Compute the blocking probability (lane full).\n"
        f"PK = p_0 * rho^{K} = {P0:.5f} * ({rho:.4f})^{K} = {PK:.5f}\n\n"
        f"**Step 5:** Compute the effective arrival rate of cars that "
        f"join.\n"
        f"lam_e = lambda * (1 - PK) = {lam} * (1 - {PK:.5f}) "
        f"= {lam_e:.3f} cars/hour\n\n"
        f"**Step 6:** Compute the expected number of cars in the lane "
        f"directly from the derived distribution.\n"
        f"L = p_0 * ({sl_terms})\n"
        f"The weighted sum is SL = {SL:.4f}, so "
        f"L = {P0:.5f} * {SL:.4f} = {L:.4f} cars\n\n"
        f"**Step 7:** Apply Little's formula with the effective arrival "
        f"rate.\n"
        f"W = L / lam_e = {L:.4f} / {lam_e:.3f} = {W_hr:.5f} hours\n\n"
        f"**Step 8:** Convert the time in system to minutes.\n"
        f"W = {W_hr:.5f} hours * 60 minutes/hour = {W_min:.2f} minutes\n\n"
        f"**Answer:** The average time in the system for cars that join "
        f"the lane is {W_min:.2f} minutes"
    )

    return question, solution
