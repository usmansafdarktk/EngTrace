import math
import random
from decimal import Decimal, ROUND_HALF_UP

from data.templates.branches.industrial_engineering.constants import (
    QUEUE_SCENARIOS,
)


def _hu(x, places):
    """Half-up rounding of a float via its shortest decimal repr."""
    q = Decimal("1") if places == 0 else Decimal("0." + "0" * (places - 1) + "1")
    v = Decimal(repr(x)).quantize(q, rounding=ROUND_HALF_UP)
    return int(v) if places == 0 else float(v)


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

    Trace integrity (Layer 0, 2026-09-23):
        rho = lambda/mu (3 dp), L = lambda/(mu - lambda) (3 dp) and
        W = L/lambda (4 dp) are quotients by arbitrary integers, exact at
        no fixed display (13/16 = 0.8125 terminates, 13/21 does not), and
        the minute answer W*60 (a 4-dp value times 60, exact at 3 dp) is
        quoted at 1 dp, so lengthening it would change what the item asks
        for (D-044). A draw on which any of the four sits on a half-way tie
        at its display is therefore redrawn rather than rounded either way
        (D-016). Enumerated over all 1,317 reachable (scenario, mu, lambda)
        triples: 89 tie on at least one line, a sampling-weighted rejection
        of 5.9% (rho 2.1%, W in hours 1.6%, W in minutes 2.2%, L 1.0%).
        The minute-answer tie is invisible to T1, whose parser skips the
        "* 60 minutes/hour" line.

    Returns:
        tuple(str, str): (question, solution)
    """
    for _attempt in range(200):
        scenario_key = random.choice(sorted(_MM1_SETTINGS))
        setting = _MM1_SETTINGS[scenario_key]
        lam_lo, lam_hi = QUEUE_SCENARIOS[scenario_key]["lam_hr"]
        mu_lo, mu_hi = QUEUE_SCENARIOS[scenario_key]["mu_hr"]

        # Sample the service rate first, then the arrival rate inside the
        # per-sample window that guarantees 0.55 <= rho <= 0.92 exactly
        # (Civil lesson 1: per-sample joint feasibility, no feasibility
        # rejection; the loop around this block only removes display ties).
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
        # D-016: rho, L and W in hours are quotients exact at no fixed
        # display, and the 1-dp minute answer is the item's quoted precision
        # (D-044); a draw on which any of them sits on a half-way tie at its
        # display is redrawn rather than rounded either way. Screened last,
        # after every draw, so only a tie draw is redrawn.
        if (_is_display_tie(lam / mu, 3)
                or _is_display_tie(lam / (mu - lam), 3)
                or _is_display_tie(L / lam, 4)
                or _is_display_tie(W_hr * 60, 1)):
            continue
        break
    else:
        raise AssertionError("resample loop exhausted")

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

    Trace integrity (Layer 0, 2026-09-23):
        Wq in hours, Lq/lambda, is a 4-dp value over a small integer and
        was printed at 5 dp; on the five reachable instances with
        lambda = 4 or 20 the quotient is exact at 6 dp and sat on a
        half-way tie (0.8891 / 4 = 0.222275). No 6- or 7-dp display is
        tie-free either (the lambda = 32 and 16 quotients terminate at 7
        and 8 dp), so Wq is bound half-up at 8 dp, the length of the
        longest terminating quotient in the reachable set, and printed at
        8 dp in Steps 5 and 6. Enumerated over all 82 (combo, lambda)
        instances: no line ties at 8 dp, and the 2-dp minute answer is
        unchanged everywhere (D-016/D-037).
        rho was bound from lambda/(c*mu) and printed at 4 dp, while the
        Step 2 line shows it as a / c with a already at 4 dp: a 4-dp a over
        c = 2 is exact at 5 dp and sat on a 4-dp half-way tie on 16 of the
        82 instances (1.4667 / 2 = 0.73335). rho is now bound half-up at
        5 dp to the DISPLAYED a / c (D-016 part 2) and printed at 5 dp in
        Steps 2-4, so the chain consumes what the reader sees; over c = 3
        the quotient can never tie at any display (2*a*10^4 = 3*(2k+1)
        has no integer solution). Enumerated over all 82 instances: no
        rho tie at 5 dp, no P0/Lq/Wq tie, and the 2-dp minute answer
        moves on 6 of the 82 instances (accepted, Layer 0 round 3).

    Returns:
        tuple(str, str): (question, solution)
    """
    key, c, ts, lam_lo, lam_hi = random.choice(_MMC_COMBOS)
    setting = _MMC_SETTINGS[key]
    lam = random.randint(lam_lo, lam_hi)

    # Round-then-recompute: the gold chain derives only from the presented
    # (lam, ts, c); display precisions are sized to the (1-rho)^-2
    # amplification (lessons 5/30/33): a at 4 dp; rho at 5 dp; P0 at 5 dp;
    # Lq 4 dp.
    mu = 60 // ts                      # exact integer by construction
    a = round(lam / mu, 4)
    # rho is bound to the DISPLAYED a / c (D-016 part 2): a 4-dp a over
    # c = 2 is exact at 5 dp, and over c = 3 it can never sit on a half-way
    # tie at any display, so rho is half-up at 5 dp and tie-free on all 82
    # reachable instances (at 4 dp the a / c line tied on 16 of them).
    rho = _hu(a / c, 5)
    P0 = round(_mmc_p0(a, rho, c), 5)
    Lq = round(P0 * a ** c * rho / (math.factorial(c) * (1.0 - rho) ** 2), 4)
    # Lq/lambda is displayed at 8 dp: over the 82 reachable (combo, lambda)
    # instances every terminating quotient terminates within 8 dp (the
    # longest are /16 and /32), so nothing rounds and no half-way tie can
    # arise. At 5 dp the lambda = 4 and 20 instances (exact at 6 dp) tied;
    # 6 dp would move the tie to lambda = 32, 7 dp to lambda = 16 and 32.
    Wq_hr = _hu(Lq / lam, 8)
    Wq_min = round(Wq_hr * 60, 2)

    # Physical bounds (docstring, verbatim; margins cover displayed rounding)
    assert 0.60 <= lam / (c * mu) <= 0.82, f"rho out of bounds: {lam}/({c}*{mu})"
    assert 0.045 <= P0 <= 0.26, f"P0 out of bounds: {P0}"
    assert 0.60 <= Lq <= 3.40, f"Lq out of bounds: {Lq}"
    assert 1.9 <= Wq_min <= 40.0, f"Wq out of bounds: {Wq_min} min"

    if c == 2:
        p0_eq = (
            f"P0 = 1 / (1 + a + a^2/(2*(1-rho))) "
            f"= 1 / (1 + {a:.4f} + ({a:.4f})^2/(2*(1 - {rho:.5f}))) "
            f"= {P0:.5f}"
        )
        lq_eq = (
            f"Lq = P0 * a^2 * rho / (2! * (1-rho)^2) "
            f"= {P0:.5f} * ({a:.4f})^2 * {rho:.5f} / (2 * (1 - {rho:.5f})^2) "
            f"= {Lq:.4f} customers"
        )
    else:
        p0_eq = (
            f"P0 = 1 / (1 + a + a^2/2 + a^3/(6*(1-rho))) "
            f"= 1 / (1 + {a:.4f} + ({a:.4f})^2/2 + "
            f"({a:.4f})^3/(6*(1 - {rho:.5f}))) = {P0:.5f}"
        )
        lq_eq = (
            f"Lq = P0 * a^3 * rho / (3! * (1-rho)^2) "
            f"= {P0:.5f} * ({a:.4f})^3 * {rho:.5f} / (6 * (1 - {rho:.5f})^2) "
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
        f"rho = a / c = {a:.4f} / {c} = {rho:.5f}\n"
        f"Since rho = {rho:.5f} < 1, a steady state exists.\n\n"
        f"**Step 3:** Compute the probability that the system is empty.\n"
        f"{p0_eq}\n\n"
        f"**Step 4:** Compute the average number waiting in the queue.\n"
        f"{lq_eq}\n\n"
        f"**Step 5:** Apply Little's formula to the queue.\n"
        f"Wq = Lq / lambda = {Lq:.4f} / {lam} = {Wq_hr:.8f} hours\n\n"
        f"**Step 6:** Convert the waiting time to minutes.\n"
        f"Wq = {Wq_hr:.8f} hours * 60 minutes/hour = {Wq_min:.2f} minutes\n\n"
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

    Trace integrity (Layer 0, 2026-09-23):
        rho1 = lambda/muA and rho2 = lambda/(2*muB) are quotients of small
        integers exact at no fixed display, and at their 4-dp display
        they sit on a half-way tie on 92 of the 1,844 reachable triples
        (33 on rho1, 63 on rho2; 27 / 32 = 0.84375). Such a draw is
        REDRAWN rather than rounded either way (D-016); the rejection is
        5.0% of draws and leaves the winner split at 889 single / 863
        pair. L2 at 4 dp never ties. The W1 and W2 hours-to-minutes
        lines are registered parser limits and are untouched.

    Screen pass 1 (2026-09-23):
        All three judges flagged "a applicant" (passport-office setting):
        the question composed "a" + the singular of cfg['arrival'] whatever
        its initial sound. The article is now chosen per noun ("an" before
        a vowel-initial noun; the reachable nouns are customer, patient and
        applicant). Only that setting's question text changes; numbers,
        sampling and the registered W1/W2 lines are untouched.

    Returns:
        tuple(str, str): (question, solution)
    """
    for _attempt in range(200):
        muA, muB, lam = random.choice(_SEL_COMBOS)
        # D-016: a utilization quotient on a 4-dp half-way tie is redrawn
        # (92 of 1,844 triples); no other value in this chain can tie at
        # its display except the registered W1/W2 lines.
        if _is_display_tie(lam / muA, 4) or _is_display_tie(lam / (2 * muB), 4):
            continue
        break
    else:
        raise AssertionError("resample loop exhausted")
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

    # Screen pass 1 (2026-09-23): the article follows the noun's initial
    # sound -- "an applicant", "a customer", "a patient".
    unit = cfg["arrival"][:-1]
    article = "an" if unit[:1].lower() in "aeiou" else "a"

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
        f"a steady state exists for both, and compare the average time "
        f"{article} {unit} spends in the system (waiting plus "
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

    Trace integrity (Layer 0, 2026-09-23):
        lam_e = lambda*(1 - PK), an integer times a 5-dp probability, is
        exact at 5 dp and was printed at 3 dp, where it sat on a half-way
        tie on ~4% of draws (50 * (1 - 0.06047) = 46.9765); it is now bound
        half-up and printed at 5 dp in Steps 5 and 7, so nothing rounds
        (D-037); this moves the 2-dp minute answer by 0.01 on 8 of the
        1,917 surviving combos. rho = lambda/mu is a quotient exact at no
        fixed display: 550 of the 664 reachable (lambda, mu) pairs never
        terminate, and although the terminating ones end within 6 dp (so a
        6-dp display would carry no tie), rho anchors every later sum and a
        6-dp rho moves the answer on 42 combos, which D-044 rules out. A
        draw whose rho sits on a 4-dp half-way tie (32 of 1,949 combos, all
        with mu = 32 or 64) is therefore redrawn. The weighted sum SL (4 dp;
        a reader recomputes it from the printed rho) and the minute answer
        W*60 (2 dp, the item's quoted precision) can tie as well (30 and 41
        combos) and are screened the same way; both are invisible to T1.
        Enumerated over all combos, no other printed value can tie; the
        rejection is 5.9% (over) / 4.9% (under) of draws.

    Returns:
        tuple(str, str): (question, solution)
    """
    regime = random.choice(["under", "over"])
    for _attempt in range(200):
        lam, mu, K = random.choice(_MM1K_OVER if regime == "over" else _MM1K_UNDER)

        # Round-then-recompute from the presented (lam, mu, K): rho displayed
        # at 4 dp anchors the chain; the normalization sum S and weighted sum
        # SL are computed from displayed rho and themselves displayed at 4 dp.
        rho = round(lam / mu, 4)
        S = round(sum(rho ** n for n in range(K + 1)), 4)
        P0 = round(1.0 / S, 5)
        PK = round(P0 * rho ** K, 5)
        # lambda*(1 - PK) is exact at 5 dp (an integer times a 5-dp value);
        # its former 3-dp display sat on a half-way tie on ~4% of draws, so
        # it is bound and printed at 5 dp and nothing rounds (D-037).
        lam_e = _hu(lam * (1 - PK), 5)
        SL_exact = sum(n * rho ** n for n in range(1, K + 1))
        SL = round(SL_exact, 4)
        L = round(P0 * SL, 4)
        W_hr = round(L / lam_e, 5)
        W_min = round(W_hr * 60, 2)
        # D-016: rho (a quotient), SL (which a reader recomputes from the
        # printed rho) and the 2-dp minute answer may not sit on a half-way
        # tie at their display; the draw is repeated rather than rounded
        # either way. Screened last, after the draw, so only a tie draw is
        # redrawn. Enumerated over all 1,949 combos, no other printed value
        # can tie (S, p_0, PK, lam_e at 5 dp, L, W in hours).
        if (_is_display_tie(lam / mu, 4)
                or _is_display_tie(SL_exact, 4)
                or _is_display_tie(W_hr * 60, 2)):
            continue
        break
    else:
        raise AssertionError("resample loop exhausted")

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
        f"= {lam_e:.5f} cars/hour\n\n"
        f"**Step 6:** Compute the expected number of cars in the lane "
        f"directly from the derived distribution.\n"
        f"L = p_0 * ({sl_terms})\n"
        f"The weighted sum is SL = {SL:.4f}, so "
        f"L = {P0:.5f} * {SL:.4f} = {L:.4f} cars\n\n"
        f"**Step 7:** Apply Little's formula with the effective arrival "
        f"rate.\n"
        f"W = L / lam_e = {L:.4f} / {lam_e:.5f} = {W_hr:.5f} hours\n\n"
        f"**Step 8:** Convert the time in system to minutes.\n"
        f"W = {W_hr:.5f} hours * 60 minutes/hour = {W_min:.2f} minutes\n\n"
        f"**Answer:** The average time in the system for cars that join "
        f"the lane is {W_min:.2f} minutes"
    )

    return question, solution
