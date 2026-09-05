import math
import random
from decimal import Decimal, ROUND_HALF_UP

from pilot.templates.branches.industrial_engineering.constants import (
    AGGREGATE_COSTS_USD,
    AGGREGATE_MONTHLY_DEMAND,
    LINE_DEMAND_PER_SHIFT,
    WORKER_MONTHLY_OUTPUT,
)


def _hu(x, places):
    """Half-up rounding via Decimal from the shortest float repr."""
    q = Decimal("1") if places == 0 else Decimal("0." + "0" * (places - 1) + "1")
    v = Decimal(repr(x)).quantize(q, rounding=ROUND_HALF_UP)
    return int(v) if places == 0 else float(v)


# Template 18 (Easy) — Area P3: Production Planning
def template_takt_time_line_efficiency():
    """
    Takt Time, Theoretical Minimum Stations, and Line Efficiency

    Scenario:
        An assembly line must meet a per-shift demand D with available
        working time A minutes per shift. The takt time, the theoretical
        minimum number of workstations for total work content W, and the
        line efficiency at that station count are

            takt = A*60 / D              (seconds per unit)
            N_min = ceil(W / takt)
            efficiency = W / (N_min * takt) * 100%

        Requested: the line efficiency at the theoretical minimum
        station count, in percent to one decimal.

    Difficulty: Easy
    Grounding: Nahmias, Production and Operations Analysis, 7th ed.,
        Sec. 9.10 "Assembly Line Balancing" (cycle time, minimum
        stations, balance efficiency; visually verified in the on-disk
        copy, text p. 528 ff. — typology only). Cross-ref NIST/SEMATECH
        e-Handbook (none needed; pure integer arithmetic).
    Physical bounds: shift 480 minutes minus a break allowance b in
        {0, 30, 60}; per-shift demand D sampled from the divisors of
        A*60 inside LINE_DEMAND_PER_SHIFT [100, 800], so takt = A*60/D
        is an EXACT integer number of seconds (chain is exact integer
        arithmetic; the only rounding is the final 1-dp half-up on the
        efficiency percentage). Work content W integer seconds sampled
        so that N_min lands in [3, 8] with the ceiling decision
        boundary-safe: W in [takt*(N-1) + m, takt*N - m] with margin
        m = ceil(0.1*takt) (lesson 23). Efficiency = W/(N*takt) then
        lies in [(N-1)/N + 0.1/N, 1 - 0.1/N] — analytic envelope
        [70.0%, 98.8%] over N in [3, 8]; asserts: takt in [32, 288] s
        (divisor extremes: 28800/800 = 36 down to A*60/D floors; margin
        low 32), N_min in [3, 8], efficiency in [69.5, 99.0].

    Returns:
        tuple(str, str): (question, solution)
    """
    b = random.choice([0, 30, 60])
    A = 480 - b
    total_s = A * 60
    divisors = [d for d in range(LINE_DEMAND_PER_SHIFT[0],
                                 LINE_DEMAND_PER_SHIFT[1] + 1)
                if total_s % d == 0]
    D = random.choice(divisors)
    takt = total_s // D                    # exact integer seconds

    N = random.randint(3, 8)
    m = math.ceil(0.1 * takt)
    W = random.randint(takt * (N - 1) + m, takt * N - m)

    ratio = W / takt
    N_min = math.ceil(ratio)
    eff = float((Decimal(W) / (Decimal(N_min) * takt) * 100)
                .quantize(Decimal("0.1"), rounding=ROUND_HALF_UP))

    assert 32 <= takt <= 288, f"takt out of bounds: {takt}"
    assert N_min == N and 3 <= N_min <= 8, f"N_min mismatch: {N_min} vs {N}"
    assert 69.5 <= eff <= 99.0, f"efficiency out of bounds: {eff}"

    break_phrase = ("with no scheduled breaks" if b == 0 else
                    f"with {b} minutes of scheduled breaks")
    question = (
        f"An assembly line runs one 480-minute shift per day {break_phrase}, "
        f"leaving {A} minutes of working time. The line must produce {D} "
        f"units per shift, and the total work content of one unit is {W} "
        f"seconds. Determine the line efficiency (total work content "
        f"divided by the time capacity of the stations) that results from "
        f"staffing the line with the THEORETICAL MINIMUM number of "
        f"workstations, in percent to one decimal (round half up). In "
        f"your solution, compute the takt time in seconds and the "
        f"theoretical minimum number of stations."
    )

    solution = (
        f"**Given:**\n"
        f"Working time A = {A} minutes/shift; demand D = {D} units/shift; "
        f"work content W = {W} seconds/unit.\n\n"
        f"**Step 1:** Compute the takt time — the pace at which units "
        f"must leave the line.\n"
        f"takt = A*60 / D = {A} * 60 / {D} = {total_s} / {D} = {takt} "
        f"seconds per unit\n\n"
        f"**Step 2:** Compute the theoretical minimum number of "
        f"workstations. Each station contributes at most one takt of "
        f"work per unit, so\n"
        f"N_min = ceil(W / takt) = ceil({W} / {takt}) = "
        f"ceil({ratio:.4f}) = {N_min} stations\n\n"
        f"**Step 3:** Compute the line efficiency at N_min stations.\n"
        f"efficiency = W / (N_min * takt) * 100 = {W} / ({N_min} * "
        f"{takt}) * 100 = {eff:.1f}%\n\n"
        f"**Answer:** The line efficiency at the theoretical minimum "
        f"number of stations is {eff:.1f} percent"
    )

    return question, solution


# Fixed 5-task precedence network for the line-balancing template:
#   a -> b, a -> c;  b and c -> d;  d -> e.
_T19_PRED = {"a": [], "b": ["a"], "c": ["a"], "d": ["b", "c"], "e": ["d"]}
_T19_ORDER = ("a", "b", "c", "d", "e")


def _t19_assign(times, CT):
    """Execute the longest-eligible-task-that-fits rule; return the
    station list [[task, ...], ...] and a per-station decision log
    [(station_no, [(eligible_snapshot, chosen_or_None, rem_after)])]."""
    assigned = []
    stations = []
    log = []
    while len(assigned) < 5:
        rem = CT
        station = []
        decisions = []
        while True:
            eligible = [x for x in _T19_ORDER
                        if x not in assigned
                        and all(p in assigned for p in _T19_PRED[x])]
            fitting = [x for x in eligible if times[x] <= rem]
            snapshot = [(x, times[x], times[x] <= rem) for x in eligible]
            if not fitting:
                decisions.append((snapshot, None, rem))
                break
            pick = max(fitting, key=lambda x: times[x])
            rem -= times[pick]
            assigned.append(pick)
            station.append(pick)
            decisions.append((snapshot, pick, rem))
        stations.append(station)
        log.append(decisions)
    return stations, log


# Template 19 (Advanced) — Area P3: Production Planning
# [carries the domain's second Advanced slot per the t17 reconciliation]
def template_line_balancing_heuristic():
    """
    Assembly-Line Balancing by the Longest-Eligible-Task Rule

    Scenario:
        Five assembly tasks with distinct integer durations and the
        precedence network a -> {b, c}; {b, c} -> d; d -> e must be
        grouped into workstations under a given cycle time CT. The
        solver CONSTRUCTS the assignment by executing the stated rule
        station by station over the network — at each point, among
        unassigned tasks whose predecessors are all assigned and that
        fit in the station's remaining time, assign the longest; when
        none fits, open the next station — then evaluates

            balance delay = (n*CT - sum(t_i)) / (n*CT) * 100%

        Requested: the balance delay in percent to one decimal.

    Difficulty: Advanced
        (Construction: the station-by-station assembly of the solution
        from the precedence network IS the derivation — nothing is
        substituted until the structure has been built; >= 6 steps.)
    Grounding: Nahmias, Production and Operations Analysis, 7th ed.,
        Sec. 9.10 "Assembly Line Balancing" (ranked-heuristic
        assignment under precedence, balance delay; visually verified
        in the on-disk copy, text pp. 528-533 — typology only, never
        its task networks or numbers).
    Physical bounds: task durations 5 DISTINCT integers in [15, 110] s
        (inside LINE_TASK_TIME_S; distinctness makes every "longest"
        choice unique); CT integer in [max task + 5, 150] s. Screens
        (bounded resample loop, cap 300): the executed rule yields
        n in {3, 4} stations; balance delay in [2, 40]%; and
        N_min = ceil(sum/CT) in {n, n - 1} is reported against the
        achieved n in the trace. All arithmetic is EXACT integer
        comparison and subtraction — fit decisions cannot wobble — and
        the single final rounding is a Decimal half-up of the exact
        rational percentage (lessons 51/65/71). Author QA 2026-08-07
        (20,000-seed sweep): n mix and delay range recorded in the
        review log; asserts: delay in [1.9, 40.5]; n in {3, 4};
        sum(t) <= n*CT.

    Returns:
        tuple(str, str): (question, solution)
    """
    for _ in range(300):
        ts = random.sample(range(15, 111), 5)
        times = dict(zip(_T19_ORDER, ts))
        total = sum(ts)
        ct_lo = max(ts) + 5
        ct_hi = 150
        if ct_lo > ct_hi:
            continue
        CT = random.randint(ct_lo, ct_hi)
        stations, log = _t19_assign(times, CT)
        n = len(stations)
        if n not in (3, 4):
            continue
        delay_exact = (Decimal(n * CT - total) / (n * CT)) * 100
        delay = float(delay_exact.quantize(Decimal("0.1"),
                                           rounding=ROUND_HALF_UP))
        if not (2.0 <= delay <= 40.0):
            continue
        break
    else:
        raise AssertionError("resample loop exhausted")

    N_min = math.ceil(total / CT)
    assert 1.9 <= delay <= 40.5, f"delay out of bounds: {delay}"
    assert n in (3, 4), f"n out of bounds: {n}"
    assert total <= n * CT, "capacity violated"
    assert len(set(ts)) == 5, "task times not distinct"

    task_list = "; ".join(f"task {x}: {times[x]} s" for x in _T19_ORDER)
    question = (
        f"An assembly line performs five tasks with these durations: "
        f"{task_list}. Precedence: task a must be completed before b and "
        f"before c; task d requires both b and c; task e requires d. The "
        f"line operates at a cycle time of {CT} seconds. Group the tasks "
        f"into workstations using the LONGEST-ELIGIBLE-TASK rule: fill "
        f"one station at a time; at each point, among the unassigned "
        f"tasks whose predecessors have all been assigned (at an earlier "
        f"station or earlier at the same station) AND whose duration "
        f"fits in the station's remaining time, assign the one with the "
        f"longest duration; when no eligible task fits, open the next "
        f"station. Then determine the balance delay of your line — the "
        f"idle fraction (n*CT - total work) / (n*CT) over the n stations "
        f"you used — in percent to one decimal (round half up). In your "
        f"solution, state the theoretical minimum number of stations and "
        f"the full assignment."
    )

    # Render the station-by-station construction.
    step_lines = []
    step_no = 2
    for st_idx, decisions in enumerate(log, start=1):
        parts = [f"**Step {step_no}:** Station {st_idx} (remaining time "
                 f"{CT} s)."]
        for snapshot, pick, rem in decisions:
            if not snapshot:
                parts.append("No tasks remain.")
                break
            desc = ", ".join(f"{x} ({tt} s{'' if fits else ', does not fit'})"
                             for x, tt, fits in snapshot)
            if pick is None:
                parts.append(f"Eligible: {desc} — none fits, so close "
                             f"the station.")
            else:
                parts.append(f"Eligible: {desc} — assign {pick} "
                             f"(remaining {rem} s).")
        station_str = ", ".join(stations[st_idx - 1])
        parts.append(f"Station {st_idx} contents: {{{station_str}}}.")
        step_lines.append("\n".join(parts))
        step_no += 1

    concord = ("matches" if n == N_min else "is one above")
    solution = (
        f"**Given:**\n"
        f"Task durations: {task_list}; precedence a -> b, a -> c; "
        f"b, c -> d; d -> e; cycle time CT = {CT} s.\n\n"
        f"**Step 1:** Total work content and the theoretical minimum "
        f"number of stations.\n"
        f"sum(t) = {' + '.join(str(times[x]) for x in _T19_ORDER)} = "
        f"{total} s;  N_min = ceil({total} / {CT}) = {N_min} stations "
        f"(stations must be whole, so round up)\n\n"
        + "\n\n".join(step_lines) + "\n\n"
        f"**Step {step_no}:** Compare with the bound and compute the "
        f"balance delay. The heuristic used n = {n} stations, which "
        f"{concord} the theoretical minimum of {N_min}.\n"
        f"balance delay = (n*CT - sum(t)) / (n*CT) * 100 = ({n}*{CT} - "
        f"{total}) / ({n}*{CT}) * 100 = {n * CT - total} / {n * CT} * "
        f"100 = {delay:.1f}%\n\n"
        f"**Answer:** The balance delay of the assembled line is "
        f"{delay:.1f} percent"
    )

    return question, solution


# Template 20 (Advanced) — Area P3: Production Planning
# [BRANCHING: chase-vs-level winner is cost-regime dependent]
def template_chase_vs_level_aggregate():
    """
    Aggregate Planning: Chase vs. Level Strategy over Three Months

    Scenario:
        A plant faces three months of demand d1, d2, d3 with worker
        productivity q units per worker-month, initial workforce W0,
        hiring cost cH and firing cost cF per worker, and holding cost
        h per unit of ending inventory per month. Two prescribed plans
        are costed and compared:

          CHASE: size month t's workforce as W_t = ceil(d_t / q)
            (ignoring carried inventory, as prescribed); pay hiring/
            firing on each change from the previous month (from W0 for
            month 1); ending inventory accumulates the small ceiling
            surpluses I_t = I_(t-1) + q*W_t - d_t.
          LEVEL: use the smallest CONSTANT workforce W_L that keeps
            cumulative production >= cumulative demand every month,
            W_L = max_t ceil((d1+...+dt) / (t*q)) — a max-over-prefixes
            construction; pay one workforce adjustment from W0; ending
            inventories I_t = t*q*W_L - (d1+...+dt).

        Both plans' total costs (workforce changes + holding) are
        computed and the cheaper plan identified — which one wins is
        cost-regime dependent (BRANCHING). Requested: the total cost of
        the cheaper plan, in whole dollars.

    Difficulty: Advanced
        (Construction: two complete multi-period plans are BUILT —
        workforce paths, inventory trajectories, cost accumulations —
        including the level plan's max-over-prefixes feasibility
        argument, before the final comparison; >= 6 steps.)
    Grounding: Nahmias, Production and Operations Analysis, 7th ed.,
        Sec. 3.4-3.5 (aggregate planning; zero-inventory/chase vs
        constant-workforce plans; visually verified in the on-disk
        copy, text p. 138 ff. — typology only, never its numbers).
    Physical bounds: monthly demands integer in [200, 5000]
        (AGGREGATE_MONTHLY_DEMAND) with month-to-month swings of
        10-45%; q integer in [20, 200] (WORKER_MONTHLY_OUTPUT); cH in
        [300, 1500], cF in [500, 2500], h in [2, 20] (AGGREGATE_COSTS
        windows, integers); W0 = chase month-1 need +/- up to 3
        workers. Screens (bounded loop, cap 400): every ceil operand's
        fractional part in [0.1, 0.9] (lesson 23 — both d_t/q and all
        prefix averages); a 50/50-sampled cost regime steers the
        winner (churn-heavy: high cH/cF with low h; holding-heavy: the
        reverse) and the achieved winner mix is verified by sweep;
        decisive margin |C_chase - C_level| >= 3% of the smaller. All
        arithmetic is EXACT integers; the answer needs no rounding at
        all. Analytic dominance floor: the prefix-fraction screen forces
        at least ~0.1*q units of ending inventory somewhere, so both
        totals exceed ~$10 even with zero workforce changes; the ~$800k
        ceiling dominates worst-case churn (~250 workers x $2500) plus
        holding. Author QA (20,000-seed sweep) winner mix and attained
        cost range recorded in the review log; asserts: both totals in
        [10, 800000]; margin >= 2.9%; all inventories >= 0.

    Returns:
        tuple(str, str): (question, solution)
    """
    for _ in range(400):
        q = random.randint(*WORKER_MONTHLY_OUTPUT)
        d1 = random.randint(*AGGREGATE_MONTHLY_DEMAND)
        swing2 = random.uniform(0.10, 0.45) * random.choice([-1, 1])
        swing3 = random.uniform(0.10, 0.45) * random.choice([-1, 1])
        d2 = int(round(d1 * (1 + swing2)))
        d3 = int(round(d2 * (1 + swing3)))
        if not (200 <= d2 <= 5000 and 200 <= d3 <= 5000):
            continue
        ds = (d1, d2, d3)
        # ceil decisiveness for chase needs and level prefixes
        ok = True
        for t in (1, 2, 3):
            fr = (ds[t - 1] / q) % 1
            if not (0.1 <= fr <= 0.9):
                ok = False
            cum = sum(ds[:t])
            fr2 = (cum / (t * q)) % 1
            if not (0.1 <= fr2 <= 0.9):
                ok = False
        if not ok:
            continue

        # cost regime steering (50/50): churn-heavy favors level,
        # holding-heavy favors chase
        regime = random.choice(["churn", "holding"])
        hire_w = AGGREGATE_COSTS_USD["hire_per_worker"]
        fire_w = AGGREGATE_COSTS_USD["fire_per_worker"]
        hold_w = AGGREGATE_COSTS_USD["hold_per_unit_month"]
        if regime == "churn":
            cH = random.randint(900, hire_w[1])
            cF = random.randint(1500, fire_w[1])
            h = random.randint(hold_w[0], 6)
        else:
            cH = random.randint(hire_w[0], 700)
            cF = random.randint(fire_w[0], 1000)
            h = random.randint(10, hold_w[1])

        # CHASE plan
        Wc = [math.ceil(d / q) for d in ds]
        W0 = Wc[0] + random.randint(-3, 3)
        if W0 < 1:
            continue
        prev = W0
        churn_cost = 0
        changes = []
        for w in Wc:
            delta = w - prev
            changes.append(delta)
            churn_cost += cH * delta if delta > 0 else cF * (-delta)
            prev = w
        inv = 0
        chase_invs = []
        for t in range(3):
            inv = inv + q * Wc[t] - ds[t]
            chase_invs.append(inv)
        if min(chase_invs) < 0:
            continue
        chase_hold = h * sum(chase_invs)
        C_chase = churn_cost + chase_hold

        # LEVEL plan
        prefixes = [math.ceil(sum(ds[:t]) / (t * q)) for t in (1, 2, 3)]
        W_L = max(prefixes)
        delta_L = W_L - W0
        level_adjust = cH * delta_L if delta_L > 0 else cF * (-delta_L)
        level_invs = [t * q * W_L - sum(ds[:t]) for t in (1, 2, 3)]
        if min(level_invs) < 0:
            continue
        level_hold = h * sum(level_invs)
        C_level = level_adjust + level_hold

        lo, hi = min(C_chase, C_level), max(C_chase, C_level)
        if lo <= 0 or (hi - lo) < 0.03 * lo:
            continue
        break
    else:
        raise AssertionError("resample loop exhausted")

    winner = "chase" if C_chase < C_level else "level"
    C_best = min(C_chase, C_level)

    assert 10 <= C_chase <= 800000 and 10 <= C_level <= 800000, \
        f"costs out of bounds: {C_chase}, {C_level}"
    assert (max(C_chase, C_level) - C_best) >= 0.029 * C_best, "margin"
    assert min(chase_invs) >= 0 and min(level_invs) >= 0, "backorder"

    question = (
        f"A plant must plan aggregate production for three months with "
        f"demands of {d1}, {d2}, and {d3} units. Each worker produces "
        f"q = {q} units per month. The plant starts with W0 = {W0} "
        f"workers; hiring costs cH = ${cH} per worker, firing costs "
        f"cF = ${cF} per worker, and ending inventory costs h = ${h} "
        f"per unit per month. Evaluate two plans. CHASE: each month t, "
        f"set the workforce to the smallest integer covering that "
        f"month's demand alone, W_t = ceil(d_t/q) (ignore carried "
        f"inventory when sizing), paying hiring/firing on every change "
        f"from the previous month's level (from W0 in month 1). LEVEL: "
        f"use the smallest CONSTANT workforce W_L such that cumulative "
        f"production covers cumulative demand in every month, adjusted "
        f"once from W0 at the start of month 1. For both plans, charge "
        f"h on each month's ending inventory. Compute both total "
        f"three-month costs (workforce changes plus holding), identify "
        f"the cheaper plan, and report ITS total cost in whole dollars."
    )

    ch_lines = []
    for t in range(3):
        delta = changes[t]
        verb = (f"hire {delta}" if delta > 0 else
                (f"fire {-delta}" if delta < 0 else "no change"))
        ch_lines.append(f"month {t + 1}: W{t + 1} = ceil({ds[t]}/{q}) = "
                        f"{Wc[t]} ({verb})")
    churn_terms = " + ".join(
        (f"{cH}*{d}" if d > 0 else (f"{cF}*{-d}" if d < 0 else "0"))
        for d in changes)
    pref_terms = "; ".join(
        f"ceil({sum(ds[:t])}/{t * q}) = {prefixes[t - 1]}" for t in (1, 2, 3))

    solution = (
        f"**Given:**\n"
        f"Demands ({d1}, {d2}, {d3}) units; q = {q} units/worker-month; "
        f"W0 = {W0}; cH = ${cH}; cF = ${cF}; h = ${h}/unit-month.\n\n"
        f"**Step 1:** CHASE workforce path — size each month to its own "
        f"demand (workers must be whole, so round up).\n"
        f"{'; '.join(ch_lines)}\n\n"
        f"**Step 2:** CHASE workforce-change cost.\n"
        f"changes cost = {churn_terms} = ${churn_cost}\n\n"
        f"**Step 3:** CHASE inventories and holding cost. Ending "
        f"inventory I_t = I_(t-1) + q*W_t - d_t:\n"
        f"I1 = {q}*{Wc[0]} - {d1} = {chase_invs[0]}; "
        f"I2 = {chase_invs[0]} + {q}*{Wc[1]} - {d2} = {chase_invs[1]}; "
        f"I3 = {chase_invs[1]} + {q}*{Wc[2]} - {d3} = {chase_invs[2]}\n"
        f"holding = {h} * ({chase_invs[0]} + {chase_invs[1]} + "
        f"{chase_invs[2]}) = ${chase_hold};  TOTAL chase = "
        f"${churn_cost} + ${chase_hold} = ${C_chase}\n\n"
        f"**Step 4:** LEVEL workforce — the smallest constant W_L "
        f"covering every cumulative-demand prefix:\n"
        f"{pref_terms};  W_L = max = {W_L} workers "
        f"(one adjustment from W0 = {W0}: "
        f"{'hire ' + str(delta_L) if delta_L > 0 else ('fire ' + str(-delta_L) if delta_L < 0 else 'no change')}, "
        f"cost ${level_adjust})\n\n"
        f"**Step 5:** LEVEL inventories and holding cost. "
        f"I_t = t*q*W_L - cumulative demand:\n"
        f"I1 = {q}*{W_L} - {d1} = {level_invs[0]}; "
        f"I2 = {2 * q}*{W_L} - {d1 + d2} = {level_invs[1]}; "
        f"I3 = {3 * q}*{W_L} - {d1 + d2 + d3} = {level_invs[2]}\n"
        f"holding = {h} * ({level_invs[0]} + {level_invs[1]} + "
        f"{level_invs[2]}) = ${level_hold};  TOTAL level = "
        f"${level_adjust} + ${level_hold} = ${C_level}\n\n"
        f"**Step 6:** Compare the plans.\n"
        f"chase ${C_chase} vs level ${C_level}: the {winner} plan is "
        f"cheaper.\n\n"
        f"**Answer:** The cheaper plan's total three-month cost is "
        f"{C_best} dollars"
    )

    return question, solution
