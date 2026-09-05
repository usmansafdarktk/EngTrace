import random
from decimal import Decimal, ROUND_HALF_UP

from pilot.templates.branches.industrial_engineering.constants import (
    COMPONENT_RELIABILITY_CLASSES,
)


def _d3(x):
    """The exact Decimal value of a 3-dp displayed reliability."""
    return Decimal(f"{x:.3f}")


def _q4(x):
    """Round an exact Decimal half-up to 4 dp (the convention a student
    applies); returns float for display/assert use."""
    return float(x.quantize(Decimal("0.0001"), rounding=ROUND_HALF_UP))

# Topology framings for the structure-reliability template (BRANCHING: the
# sampled topology changes the governing formula and the step structure).
# Component mission reliabilities are sampled at 3 dp from the named
# constants.py class windows ([POLICY: sampling-only]; given-values rule).
_T8_TOPOLOGIES = ("series", "parallel", "mixed")


# Template 8 (Easy) — Area S3: System Reliability  [BRANCHING: topology]
def template_system_reliability_topology():
    """
    System Reliability from Component Reliabilities: Series, Parallel, or
    Mixed Structure

    Scenario:
        An industrial control system is built from independent components
        whose mission reliabilities (probability of operating without
        failure over the stated production run) are given. The sampled
        topology BRANCHES the governing formula:

            series (3 units):    Rs = R1 * R2 * R3
            parallel (2 units):  Rs = 1 - (1 - R1)(1 - R2)
            mixed:               Rp = 1 - (1 - R1)(1 - R2);  Rs = Rp * R3

        Requested: the system mission reliability.

    Difficulty: Easy
    Grounding: Ross, Introduction to Probability Models, 11th ed., Ch. 9
        (Sec. 9.2 structure functions; series and parallel systems and
        their reliability functions). Cross-ref NIST/SEMATECH e-Handbook
        Ch. 8 (apr/, system reliability practice).
    Physical bounds: component reliabilities sampled at 3 dp from
        constants.py classes — series: commercial [0.90, 0.97] x
        industrial [0.95, 0.99] x industrial [0.95, 0.99]; parallel: two
        commercial units; mixed: commercial pair + industrial unit. All
        intermediates are carried EXACTLY in decimal arithmetic (products
        of 3-dp values shown at 6 dp), and the final display rounds
        HALF-UP to 4 dp — so a full-precision solve rounded by the
        standard half-up convention matches the gold answer
        digit-for-digit, including at true decimal ties. Rs is monotone
        increasing in every R_i, so the analytic corners are the window
        endpoints (author QA 2026-08-06, half-up): series Rs in
        [0.8123, 0.9507], parallel Rs in [0.9900, 0.9991], mixed Rs in
        [0.9405, 0.9891]; asserts use those per-branch windows with
        0.0002 margin.

    Returns:
        tuple(str, str): (question, solution)
    """
    topology = random.choice(_T8_TOPOLOGIES)
    lo_c, hi_c = COMPONENT_RELIABILITY_CLASSES["commercial grade"]
    lo_i, hi_i = COMPONENT_RELIABILITY_CLASSES["industrial grade"]

    if topology == "series":
        R1 = round(random.uniform(lo_c, hi_c), 3)   # flow sensor
        R2 = round(random.uniform(lo_i, hi_i), 3)   # controller
        R3 = round(random.uniform(lo_i, hi_i), 3)   # actuator
        P12 = _d3(R1) * _d3(R2)          # exact: product of two 3-dp values
        Rs = _q4(P12 * _d3(R3))          # only rounding: final 4-dp half-up
        assert 0.8121 <= Rs <= 0.9509, f"series Rs out of bounds: {Rs}"

        question = (
            f"An industrial control loop consists of a flow sensor, a "
            f"controller, and a valve actuator connected in series: the "
            f"loop functions only if all three components function. Over "
            f"the planned production run, the components "
            f"operate independently with mission reliabilities "
            f"{R1:.3f} (sensor), {R2:.3f} (controller), and {R3:.3f} "
            f"(actuator). Determine the mission reliability of the "
            f"control loop, to four decimal places."
        )
        solution = (
            f"**Given:**\n"
            f"Series system of three independent components; mission "
            f"reliabilities R1 = {R1:.3f} (sensor), R2 = {R2:.3f} "
            f"(controller), R3 = {R3:.3f} (actuator).\n\n"
            f"**Step 1:** Identify the structure formula. A series system "
            f"functions only if every component functions, so by "
            f"independence:\n"
            f"Rs = R1 * R2 * R3\n\n"
            f"**Step 2:** Multiply the first two reliabilities.\n"
            f"R1 * R2 = {R1:.3f} * {R2:.3f} = {P12:.6f}\n\n"
            f"**Step 3:** Multiply by the third reliability.\n"
            f"Rs = {P12:.6f} * {R3:.3f} = {Rs:.4f} (rounded to four "
            f"decimals)\n\n"
            f"**Answer:** The mission reliability of the control loop is "
            f"{Rs:.4f}"
        )

    elif topology == "parallel":
        R1 = round(random.uniform(lo_c, hi_c), 3)   # pump A
        R2 = round(random.uniform(lo_c, hi_c), 3)   # pump B
        Q1 = Decimal("1") - _d3(R1)                 # exact 3-dp complement
        Q2 = Decimal("1") - _d3(R2)
        QQ = Q1 * Q2                                # exact product of 3-dp
        Rs = _q4(Decimal("1") - QQ)                 # final 4-dp half-up
        assert 0.9898 <= Rs <= 0.9993, f"parallel Rs out of bounds: {Rs}"

        question = (
            f"A cooling station is served by two redundant pumps "
            f"installed in parallel: the station functions as long as at "
            f"least one pump functions. Over the planned production run, "
            f"the pumps operate independently with "
            f"mission reliabilities {R1:.3f} and {R2:.3f}. Determine the "
            f"mission reliability of the cooling station, to four "
            f"decimal places."
        )
        solution = (
            f"**Given:**\n"
            f"Parallel system of two independent pumps; mission "
            f"reliabilities R1 = {R1:.3f}, R2 = {R2:.3f}.\n\n"
            f"**Step 1:** Identify the structure formula. A parallel "
            f"system fails only if BOTH components fail, so by "
            f"independence:\n"
            f"Rs = 1 - (1 - R1) * (1 - R2)\n\n"
            f"**Step 2:** Compute the individual failure probabilities.\n"
            f"1 - R1 = 1 - {R1:.3f} = {Q1:.3f};  "
            f"1 - R2 = 1 - {R2:.3f} = {Q2:.3f}\n\n"
            f"**Step 3:** Compute the system reliability.\n"
            f"Rs = 1 - {Q1:.3f} * {Q2:.3f} = 1 - {QQ:.6f} = {Rs:.4f}\n\n"
            f"**Answer:** The mission reliability of the cooling station "
            f"is {Rs:.4f}"
        )

    else:  # mixed: redundant pair feeding one series unit
        R1 = round(random.uniform(lo_c, hi_c), 3)   # power supply A
        R2 = round(random.uniform(lo_c, hi_c), 3)   # power supply B
        R3 = round(random.uniform(lo_i, hi_i), 3)   # controller
        Q1 = Decimal("1") - _d3(R1)
        Q2 = Decimal("1") - _d3(R2)
        QQ = Q1 * Q2                     # exact: product of two 3-dp values
        Rp = Decimal("1") - QQ           # exact 6-dp complement
        Rs = _q4(Rp * _d3(R3))           # only rounding: final 4-dp half-up
        assert 0.9403 <= Rs <= 0.9893, f"mixed Rs out of bounds: {Rs}"

        question = (
            f"A process controller is fed by two redundant power "
            f"supplies wired in parallel (the controller has power as "
            f"long as at least one supply functions), and the system "
            f"functions only if the powered controller itself also "
            f"functions. Over the planned production run, "
            f"the components operate independently with mission "
            f"reliabilities {R1:.3f} and {R2:.3f} (power supplies) and "
            f"{R3:.3f} (controller). Determine the mission reliability "
            f"of the system, to four decimal places."
        )
        solution = (
            f"**Given:**\n"
            f"Parallel pair (power supplies, R1 = {R1:.3f}, "
            f"R2 = {R2:.3f}) in series with a controller "
            f"(R3 = {R3:.3f}); all components independent.\n\n"
            f"**Step 1:** Identify the structure. The power stage is a "
            f"parallel block — it fails only if both supplies fail — and "
            f"the system is that block in series with the controller:\n"
            f"Rs = [1 - (1 - R1)(1 - R2)] * R3\n\n"
            f"**Step 2:** Reduce the parallel power stage.\n"
            f"Rp = 1 - (1 - {R1:.3f}) * (1 - {R2:.3f}) "
            f"= 1 - {Q1:.3f} * {Q2:.3f} = 1 - {QQ:.6f} = {Rp:.6f}\n\n"
            f"**Step 3:** Combine the power stage in series with the "
            f"controller.\n"
            f"Rs = Rp * R3 = {Rp:.6f} * {R3:.3f} = {Rs:.4f} (rounded to "
            f"four decimals)\n\n"
            f"**Answer:** The mission reliability of the system is "
            f"{Rs:.4f}"
        )

    return question, solution


# MTTF template (BRANCHING: series vs parallel changes the formula family
# AND the given-data style). Coherence with constants.py (lesson 53):
# series rates are given per 1000 h in [0.20, 1.00] <=> [2e-4, 1e-3]/hr,
# and the parallel mean life theta in [1000, 5000] h <=> rates
# [2e-4, 1e-3]/hr — both inside FAILURE_RATE_PER_HR["pump / mechanical"]
# = (5e-5, 1e-3).
_T9_SERIES_RATE_CENTS = (20, 100)    # per 1000 h, sampled at 2 dp
_T9_THETA_RANGE = (1000, 5000)       # integer mean life, hours


# Template 9 (Intermediate) — Area S3: System Reliability  [BRANCHING]
def template_exponential_mttf_topology():
    """
    Mean Time to Failure of an Exponential-Component System: Series vs.
    Parallel

    Scenario:
        Pumps have independent, exponentially distributed lifetimes. The
        sampled configuration BRANCHES both the reasoning chain and the
        governing result:

            series (3 pumps, failure rates lam_i per 1000 h given): the
            system fails at the FIRST failure; the minimum of independent
            exponentials is exponential with the summed rate, so
                MTTF = 1000 / (lam1 + lam2 + lam3)   [hours]

            parallel (2 identical pumps, mean life theta given): the
            system fails when BOTH have failed; time to the first
            failure is exponential with doubled rate (mean theta/2), and
            by memorylessness the survivor is good-as-new, so
                MTTF = theta/2 + theta = 1.5 * theta   [hours]

        MTTF (mean time to failure) is glossed in the question. The
        requested quantity is the system MTTF in hours.

    Difficulty: Intermediate
    Grounding: Ross, Introduction to Probability Models, 11th ed., Ch. 5
        (Sec. 5.2.3 minimum of independent exponentials; memorylessness)
        and Ch. 9 (series/parallel system lifetimes). Cross-ref
        NIST/SEMATECH e-Handbook Ch. 8 (apr/).
    Physical bounds: series rates sampled at 2 dp in [0.20, 1.00] per
        1000 h (equivalently [2e-4, 1e-3]/hr, inside the constants.py
        pump/mechanical window); parallel theta integer in [1000, 5000] h
        (rates [2e-4, 1e-3]/hr, same window). Analytic corners (chains
        exact, single final half-up rounding): series MTTF = 1000/lams
        with lams in [0.60, 3.00] -> [333.3, 1666.7] h; parallel
        MTTF = 1.5*theta -> [1500.0, 7500.0] h. Asserts use series
        [333.0, 1667.0], parallel [1500.0, 7500.0].

    Returns:
        tuple(str, str): (question, solution)
    """
    configuration = random.choice(["series", "parallel"])

    if configuration == "series":
        lam1 = random.randint(*_T9_SERIES_RATE_CENTS) / 100
        lam2 = random.randint(*_T9_SERIES_RATE_CENTS) / 100
        lam3 = random.randint(*_T9_SERIES_RATE_CENTS) / 100
        lams = round(lam1 + lam2 + lam3, 2)          # exact 2-dp sum
        mttf = float((Decimal("1000") / Decimal(f"{lams:.2f}"))
                     .quantize(Decimal("0.1"), rounding=ROUND_HALF_UP))
        assert 333.0 <= mttf <= 1667.0, f"series MTTF out of bounds: {mttf}"

        question = (
            f"A booster station moves slurry through three pumps "
            f"connected in series along one line, so the station fails "
            f"as soon as ANY one pump fails. The pumps fail "
            f"independently, with exponentially distributed lifetimes "
            f"and constant failure rates of {lam1:.2f}, {lam2:.2f}, and "
            f"{lam3:.2f} failures per 1000 hours. Determine the mean "
            f"time to failure (MTTF — the expected operating time until "
            f"the station first fails) of the station, in hours, to one "
            f"decimal place (round half up)."
        )
        solution = (
            f"**Given:**\n"
            f"Three pumps in series; independent exponential lifetimes "
            f"with failure rates lam1 = {lam1:.2f}, lam2 = {lam2:.2f}, "
            f"lam3 = {lam3:.2f} per 1000 hours.\n\n"
            f"**Step 1:** Characterize the station lifetime. The station "
            f"fails at the FIRST pump failure, i.e., at the minimum of "
            f"the three lifetimes; the minimum of independent "
            f"exponential lifetimes is exponential with the SUM of the "
            f"rates.\n\n"
            f"**Step 2:** Sum the failure rates.\n"
            f"lams = {lam1:.2f} + {lam2:.2f} + {lam3:.2f} = {lams:.2f} "
            f"failures per 1000 hours\n\n"
            f"**Step 3:** Invert the rate and convert to hours. A rate "
            f"of {lams:.2f} per 1000 hours means a mean lifetime of\n"
            f"MTTF = 1000 / {lams:.2f} = {mttf:.1f} hours (rounded half "
            f"up to one decimal)\n\n"
            f"**Answer:** The mean time to failure of the booster "
            f"station is {mttf:.1f} hours"
        )

    else:
        theta = random.randint(*_T9_THETA_RANGE)
        t_first = theta / 2                          # exact (x.0 or x.5)
        mttf = theta + t_first                       # 1.5*theta, exact
        assert 1500.0 <= mttf <= 7500.0, f"parallel MTTF out of bounds: {mttf}"

        question = (
            f"A cooling loop is served by two identical pumps in "
            f"parallel; the loop keeps operating as long as at least one "
            f"pump operates, and fails only when BOTH pumps have "
            f"failed. The pumps fail independently, with exponentially "
            f"distributed lifetimes averaging {theta} hours each. Using "
            f"the memoryless property, determine the mean time to "
            f"failure (MTTF — the expected operating time until the "
            f"loop fails) of the cooling loop, in hours, to one decimal "
            f"place."
        )
        solution = (
            f"**Given:**\n"
            f"Two identical pumps in parallel; independent exponential "
            f"lifetimes with mean theta = {theta} hours (failure rate "
            f"1/{theta} per hour each).\n\n"
            f"**Step 1:** Split the loop lifetime into two phases: the "
            f"time until the FIRST pump failure, then the remaining "
            f"lifetime of the surviving pump.\n\n"
            f"**Step 2:** Expected time to the first failure. While both "
            f"pumps run, failures compete: the minimum of two "
            f"independent exponentials is exponential with doubled rate "
            f"2/{theta}, so\n"
            f"E[first failure] = theta / 2 = {theta} / 2 = {t_first:.1f} "
            f"hours\n\n"
            f"**Step 3:** Expected remaining lifetime of the survivor. "
            f"By the memoryless property, the surviving pump is "
            f"statistically good as new at the moment the other fails, "
            f"so its expected remaining lifetime is the full mean:\n"
            f"E[remaining] = theta = {theta} hours\n\n"
            f"**Step 4:** Add the phases.\n"
            f"MTTF = {t_first:.1f} + {theta} = {mttf:.1f} hours\n\n"
            f"**Answer:** The mean time to failure of the cooling loop "
            f"is {mttf:.1f} hours"
        )

    return question, solution
