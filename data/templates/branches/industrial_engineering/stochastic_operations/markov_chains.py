import math
import random

# Markov-chain templates carry no numeric constants from constants.py by
# Stage B design: transition probabilities are sampled and fully stated in
# the question (given-values rule; see the Domain 1 note in constants.py).
# Typology anchors: Ross 11e Ch. 4 (Examples 4.1/4.3); H&L 7e Ch. 16.

# Two-state framings: (state1, state2, transition phrasing) — three distinct
# operational settings for surface variety (lesson 44). p = P(state1 ->
# state2), q = P(state2 -> state1); the requested quantity is the long-run
# fraction of periods in state1.
_TWO_STATE_SETTINGS = {
    "machine": {
        "s1": "operational", "s2": "down for repair", "period": "day",
        "intro": ("A packaging machine is inspected at the start of each "
                  "day and classified as operational or down for repair."),
        "p_phrase": ("If the machine is operational today, it breaks down "
                     "by tomorrow with probability {p:.2f}."),
        "q_phrase": ("If it is down today, the repair crew restores it by "
                     "tomorrow with probability {q:.2f}."),
        "ask": ("In the long run, what fraction of days is the machine "
                "operational?"),
    },
    "supplier": {
        "s1": "on time", "s2": "late", "period": "week",
        "intro": ("A parts supplier's weekly delivery is recorded as on "
                  "time or late, and next week's status depends only on "
                  "this week's."),
        "p_phrase": ("An on-time week is followed by a late week with "
                     "probability {p:.2f}."),
        "q_phrase": ("A late week is followed by an on-time week with "
                     "probability {q:.2f}."),
        "ask": ("In the long run, what fraction of weeks is the delivery "
                "on time?"),
    },
    "workstation": {
        "s1": "within specification", "s2": "out of specification",
        "period": "shift",
        "intro": ("A workstation's output each shift is judged within "
                  "specification or out of specification, and each "
                  "shift's status depends only on the previous shift's."),
        "p_phrase": ("A within-specification shift is followed by an "
                     "out-of-specification shift with probability {p:.2f}."),
        "q_phrase": ("An out-of-specification shift is followed by a "
                     "within-specification shift with probability {q:.2f}."),
        "ask": ("In the long run, what fraction of shifts is the output "
                "within specification?"),
    },
}


# Template 5 (Easy) — Area S2: Discrete-Time Markov Chains
def template_two_state_steady_state():
    """
    Two-State Markov Chain: Steady-State Probability

    Scenario:
        A system alternates between two states (e.g., a machine that is
        operational or down) and is observed once per period; the next
        state depends only on the current one — a two-state discrete-time
        Markov chain with

            p = P(state 1 -> state 2),   q = P(state 2 -> state 1).

        The steady-state balance equation pi1 * p = pi2 * q together with
        pi1 + pi2 = 1 gives

            pi1 = q / (p + q).

        Requested: the long-run fraction of periods spent in state 1.

    Difficulty: Easy
    Grounding: Ross, Introduction to Probability Models, 11th ed., Ch. 4
        (Sec. 4.4 limiting probabilities; two-state chains as in Examples
        4.1/4.8 typology). Cross-ref Hillier & Lieberman 7e, Ch. 16
        (steady-state probabilities).
    Physical bounds: p sampled in [0.05, 0.40], q in [0.30, 0.90] (both
        2 dp): failures/lapses are occasional, recoveries likelier than
        not. Analytic corners of pi1 = q/(p+q): min 0.30/(0.40+0.30)
        = 0.4286, max 0.90/(0.05+0.90) = 0.9474; assert pi1 in
        [0.42, 0.95].

    Returns:
        tuple(str, str): (question, solution)
    """
    key = random.choice(sorted(_TWO_STATE_SETTINGS))
    cfg = _TWO_STATE_SETTINGS[key]
    p = round(random.uniform(0.05, 0.40), 2)
    q = round(random.uniform(0.30, 0.90), 2)

    # Round-then-recompute: pi1 derives from the presented 2-dp (p, q).
    denom = round(p + q, 2)            # exact at 2 dp
    pi1 = round(q / denom, 4)
    pi2 = round(1 - pi1, 4)            # exact complement of displayed pi1

    assert 0.05 <= p <= 0.40 and 0.30 <= q <= 0.90, f"(p,q) out of bounds: {p},{q}"
    assert 0.42 <= pi1 <= 0.95, f"pi1 out of bounds: {pi1}"

    question = (
        f"{cfg['intro']} {cfg['p_phrase'].format(p=p)} "
        f"{cfg['q_phrase'].format(q=q)} Model the situation as a "
        f"two-state Markov chain and use the steady-state balance "
        f"equation to answer: {cfg['ask']} Give the steady-state "
        f"probability to four decimal places."
    )

    solution = (
        f"**Given:**\n"
        f"Two-state Markov chain observed each {cfg['period']}; state 1 = "
        f"{cfg['s1']}, state 2 = {cfg['s2']}; "
        f"P(1 -> 2) = p = {p:.2f}; P(2 -> 1) = q = {q:.2f}.\n\n"
        f"**Step 1:** Write the steady-state balance equation. In steady "
        f"state, the probability flow from state 1 to state 2 equals the "
        f"flow back:\n"
        f"pi1 * p = pi2 * q, together with pi1 + pi2 = 1\n\n"
        f"**Step 2:** Solve for pi1. Substituting pi2 = 1 - pi1:\n"
        f"pi1 * {p:.2f} = (1 - pi1) * {q:.2f}  =>  "
        f"pi1 * ({p:.2f} + {q:.2f}) = {q:.2f}  =>  "
        f"pi1 = {q:.2f} / {denom:.2f} = {pi1:.4f}\n\n"
        f"**Step 3:** Check normalization with the companion "
        f"probability.\n"
        f"pi2 = 1 - {pi1:.4f} = {pi2:.4f}; the two probabilities sum to "
        f"1, as required.\n\n"
        f"**Answer:** The long-run fraction of {cfg['period']}s in the "
        f"{cfg['s1']} state is {pi1:.4f}"
    )

    return question, solution


_BRANDS = ("A", "B", "C")


# Template 6 (Intermediate) — Area S2: Discrete-Time Markov Chains
def template_two_step_transition_probability():
    """
    Chapman-Kolmogorov: Two-Step Brand-Switching Probability

    Scenario:
        Customers buy one of three coffee brands (A, B, C) each week, and
        next week's choice depends only on this week's brand — a
        three-state Markov chain with a given one-step transition matrix.
        The two-step transition probability follows from the
        Chapman-Kolmogorov decomposition over the intermediate week:

            P2(i, j) = sum_k P(i, k) * P(k, j),   k in {A, B, C}

        The trace computes the three intermediate-brand path terms,
        sums them, and then verifies the result by computing the other
        two entries of the same P2 row and checking the row sums to 1
        (a non-circular check: each entry is computed independently).

    Difficulty: Intermediate
    Grounding: Ross, Introduction to Probability Models, 11th ed., Ch. 4,
        Sec. 4.2 (Chapman-Kolmogorov equations; n-step transition
        probabilities). Cross-ref Hillier & Lieberman 7e, Ch. 16
        (n-step transition matrices; brand-switching typology).
    Physical bounds: each row of the one-step matrix has a "loyalty"
        diagonal entry in [0.40, 0.70] (2 dp) and off-diagonal entries
        >= 0.05 (2 dp) summing the row to exactly 1. All two-step
        entries are exact 4-dp sums of products of 2-dp values. Analytic
        corners (author QA 2026-08-06): min P2 = 0.70*0.05 + 0.05*0.40
        + 0.25*0.05 = 0.0675 (i != j, everything anti-aligned); max P2 =
        0.70^2 + 0.25*0.55 + 0.05*0.55 = 0.6550 (i = j, competitors'
        return flows maximal); a 500,000-matrix Monte Carlo scan reached
        [0.0675, 0.6395], inside the analytic envelope. Asserts use
        [0.06, 0.66].

    Returns:
        tuple(str, str): (question, solution)
    """
    # Sample a loyalty-structured transition matrix: diagonal (staying
    # with the current brand) dominates, remainder split between the two
    # competitors with every entry >= 0.05 (2-dp cents arithmetic; rows
    # sum to exactly 1 by construction).
    P = {}
    for row in _BRANDS:
        loyal_c = random.randint(40, 70)                  # cents
        rest = 100 - loyal_c
        first_c = random.randint(5, rest - 5)
        others = [b for b in _BRANDS if b != row]
        entries = {row: loyal_c, others[0]: first_c, others[1]: rest - first_c}
        P[row] = {b: entries[b] / 100 for b in _BRANDS}

    i = random.choice(_BRANDS)
    j = random.choice(_BRANDS)

    # Round-then-recompute: all products/sums derive from the presented
    # 2-dp entries; products of two 2-dp values are EXACT at 4 dp, so the
    # only rounding is the final display itself.
    terms = {k: round(P[i][k] * P[k][j], 4) for k in _BRANDS}
    p2 = round(sum(terms.values()), 4)
    row_others = {}
    for jj in _BRANDS:
        if jj != j:
            row_others[jj] = round(sum(round(P[i][k] * P[k][jj], 4)
                                       for k in _BRANDS), 4)
    row_total = round(p2 + sum(row_others.values()), 4)

    for row in _BRANDS:
        assert 0.40 <= P[row][row] <= 0.70, f"loyalty out of bounds: {P[row]}"
        assert all(v >= 0.05 for v in P[row].values()), f"entry < 0.05: {P[row]}"
        assert abs(sum(P[row].values()) - 1.0) < 1e-9, f"row sum != 1: {P[row]}"
    assert 0.06 <= p2 <= 0.66, f"P2 out of bounds: {p2}"
    assert abs(row_total - 1.0) < 5e-4, f"P2 row check failed: {row_total}"

    matrix_text = "; ".join(
        f"P({r} -> {c}) = {P[r][c]:.2f}" for r in _BRANDS for c in _BRANDS
    )
    other_js = [jj for jj in _BRANDS if jj != j]
    term_lines = "\n".join(
        f"via brand {k}: P({i} -> {k}) * P({k} -> {j}) "
        f"= {P[i][k]:.2f} * {P[k][j]:.2f} = {terms[k]:.4f}"
        for k in _BRANDS
    )

    question = (
        f"A market survey tracks which of three coffee brands (A, B, C) "
        f"each customer buys every week; next week's choice depends only "
        f"on this week's brand, so the behavior is a three-state Markov "
        f"chain. The weekly transition probabilities are: {matrix_text}. "
        f"A customer buys brand {i} this week. Using the "
        f"Chapman-Kolmogorov decomposition over next week's brand, "
        f"determine the probability that this customer buys brand {j} "
        f"two weeks from now. Give the probability to four decimal "
        f"places, and verify your value by checking that the three "
        f"two-step probabilities out of brand {i} sum to 1."
    )

    solution = (
        f"**Given:**\n"
        f"Three-state weekly Markov chain over brands A, B, C with "
        f"one-step transition probabilities {matrix_text}. Current brand: "
        f"{i}; target: brand {j} two weeks from now.\n\n"
        f"**Step 1:** Write the two-step probability as a sum over the "
        f"intermediate week's brand (Chapman-Kolmogorov):\n"
        f"P2({i} -> {j}) = sum over k of P({i} -> k) * P(k -> {j}), "
        f"k in {{A, B, C}}\n\n"
        f"**Step 2:** Compute the three intermediate-brand path terms:\n"
        f"{term_lines}\n\n"
        f"**Step 3:** Sum the path terms.\n"
        f"P2({i} -> {j}) = {terms['A']:.4f} + {terms['B']:.4f} + "
        f"{terms['C']:.4f} = {p2:.4f}\n\n"
        f"**Step 4:** Verify with the row check. Computing the other two "
        f"two-step entries the same way:\n"
        f"P2({i} -> {other_js[0]}) = {row_others[other_js[0]]:.4f}; "
        f"P2({i} -> {other_js[1]}) = {row_others[other_js[1]]:.4f}\n"
        f"Row sum: {p2:.4f} + {row_others[other_js[0]]:.4f} + "
        f"{row_others[other_js[1]]:.4f} = {row_total:.4f} = 1, as "
        f"required for a probability distribution.\n\n"
        f"**Answer:** The probability that the customer buys brand {j} "
        f"two weeks from now is {p2:.4f}"
    )

    return question, solution


# Absorbing-chain template: expected time to absorption via first-step
# analysis — the trace CONSTRUCTS and solves the linear system (lesson 41:
# Advanced earned by construction). Sampling is windowed per-sample so the
# solve's denominator d1 = 1 - g - w*m/(1-s) stays >= 0.15 (precision
# lesson 30: d1 divides the answer, so its rounding error is amplified by
# 1/d1; 0.15 caps the drift at ~0.075%).
_ABS_G_RANGE = (55, 70)      # cents: P(stay Good)
_ABS_M_RANGE = (10, 30)      # cents: P(Worn -> Good) maintenance return
_ABS_D_RANGE = (10, 35)      # cents: P(Worn -> Failed)
_ABS_F_MAX = 10              # cents: P(Good -> Failed) direct failure cap


# Template 7 (Advanced) — Area S2: Discrete-Time Markov Chains
def template_absorbing_chain_time_to_failure():
    """
    Absorbing Markov Chain: Expected Weeks to Failure by First-Step
    Analysis

    Scenario:
        A CNC machine is inspected weekly and classified Good, Worn, or
        Failed; next week's condition depends only on this week's. Failed
        is absorbing (the machine is withdrawn for replacement). From
        Good the machine stays Good (g), degrades to Worn (w), or fails
        outright (f = 1 - g - w); from Worn, maintenance returns it to
        Good (m), it stays Worn (s), or it fails (d = 1 - m - s). The
        expected numbers of weeks to absorption satisfy the first-step
        system

            muG = 1 + g*muG + w*muW
            muW = 1 + m*muG + s*muW

        which the trace constructs and solves by substitution:

            muW = (1 + m*muG) / (1 - s)
            muG = (1 + w*c1) / (1 - g - w*m*c1),  c1 = 1/(1 - s)

        Requested: muG, the expected number of weeks until failure
        starting from Good.

    Difficulty: Advanced
    Grounding: Ross, Introduction to Probability Models, 11th ed., Ch. 4
        (mean time in transient states via first-step / conditioning
        analysis). Cross-ref Hillier & Lieberman 7e, Ch. 16 (absorbing
        states, first passage times).
    Physical bounds: 2-dp probabilities with g in [0.55, 0.70],
        m in [0.10, 0.30], d in [0.10, 0.35] (s = 1 - m - d in
        [0.35, 0.80]), f in a per-sample window [f_lo, 0.10] chosen so
        the exact solve denominator 1 - g - w*m/(1-s) >= 0.15 (window
        non-empty for every reachable (g, m, d) since 1 - g >= 0.30).
        Exhaustive enumeration of all 75,117 reachable parameter combos
        (author QA, 2026-08-06): muG in [4.78, 15.33] weeks, muW in
        [3.28, 13.71] weeks, max final-display deviation from the exact
        solution 0.111%, and the Step 6 back-substitution check equals
        muG exactly at 2 dp for every combo; asserts use muG
        [4.7, 15.4], muW [3.2, 13.8]; 2-dp display quantization cap
        documented at 0.15%.

    Returns:
        tuple(str, str): (question, solution)
    """
    g_c = random.randint(*_ABS_G_RANGE)
    m_c = random.randint(*_ABS_M_RANGE)
    d_c = random.randint(*_ABS_D_RANGE)
    g, m, d = g_c / 100, m_c / 100, d_c / 100
    s = round(1 - m - d, 2)              # exact cents arithmetic

    # Per-sample f window keeping the exact denominator >= 0.15:
    # w = 1 - g - f must satisfy w * m / (1 - s) <= (1 - g) - 0.15.
    w_cap = ((1 - g) - 0.15) * (1 - s) / m
    f_lo = max(2, math.ceil(round((1 - g - w_cap) * 100, 6)))
    f_c = random.randint(f_lo, _ABS_F_MAX)
    f = f_c / 100
    w = round(1 - g - f, 2)              # exact cents arithmetic

    # Round-then-recompute chain from the presented 2-dp probabilities:
    # c1 (4 dp) anchors the substitution; n1, d1 (4 dp); muG, muW (2 dp).
    c1 = round(1 / (1 - s), 4)
    n1 = round(1 + w * c1, 4)
    d1 = round(1 - g - w * m * c1, 4)
    muG = round(n1 / d1, 2)
    muW = round((1 + m * muG) * c1, 2)
    check = round(1 + g * muG + w * muW, 2)

    assert 0.55 <= g <= 0.70 and 0.02 <= f <= 0.10, f"(g,f) out of bounds: {g},{f}"
    assert 0.10 <= m <= 0.30 and 0.10 <= d <= 0.35, f"(m,d) out of bounds: {m},{d}"
    assert (1 - g - w * m / (1 - s)) >= 0.1499, "denominator screen failed"
    assert 4.7 <= muG <= 15.4, f"muG out of bounds: {muG}"
    assert 3.2 <= muW <= 13.8, f"muW out of bounds: {muW}"
    assert abs(check - muG) <= 0.06, f"back-substitution check failed: {check} vs {muG}"

    question = (
        f"A CNC machine is inspected at the start of each week and "
        f"classified as Good, Worn, or Failed; next week's condition "
        f"depends only on this week's classification. From Good, the "
        f"machine stays Good with probability {g:.2f}, degrades to Worn "
        f"with probability {w:.2f}, and fails outright with probability "
        f"{f:.2f}. From Worn, preventive maintenance returns it to Good "
        f"with probability {m:.2f}, it stays Worn with probability "
        f"{s:.2f}, and it fails with probability {d:.2f}. Failed is "
        f"absorbing: a failed machine is withdrawn from service. Set up "
        f"the first-step equations for the expected number of weeks "
        f"until failure from each working condition, solve the resulting "
        f"system, and report the expected number of weeks until failure "
        f"for a machine that is currently Good, to two decimal places."
    )

    solution = (
        f"**Given:**\n"
        f"Weekly chain over Good (G), Worn (W), Failed (F, absorbing); "
        f"P(G->G) = {g:.2f}, P(G->W) = {w:.2f}, P(G->F) = {f:.2f}; "
        f"P(W->G) = {m:.2f}, P(W->W) = {s:.2f}, P(W->F) = {d:.2f}.\n\n"
        f"**Step 1:** Set up the first-step equations. Let muG and muW be "
        f"the expected numbers of weeks until absorption (failure) "
        f"starting from Good and Worn. Conditioning on the first week's "
        f"transition (one week elapses, then the process restarts from "
        f"the new state; the F branch contributes no further time):\n"
        f"muG = 1 + {g:.2f}*muG + {w:.2f}*muW\n"
        f"muW = 1 + {m:.2f}*muG + {s:.2f}*muW\n\n"
        f"**Step 2:** Solve the second equation for muW in terms of muG.\n"
        f"muW * (1 - {s:.2f}) = 1 + {m:.2f}*muG  =>  "
        f"muW = (1 + {m:.2f}*muG) / {1 - s:.2f}\n"
        f"With 1/{1 - s:.2f} = {c1:.4f}:  "
        f"muW = (1 + {m:.2f}*muG) * {c1:.4f}\n\n"
        f"**Step 3:** Substitute into the first equation and isolate "
        f"muG.\n"
        f"muG = 1 + {g:.2f}*muG + {w:.2f}*(1 + {m:.2f}*muG)*{c1:.4f}\n"
        f"muG * (1 - {g:.2f} - {w:.2f}*{m:.2f}*{c1:.4f}) "
        f"= 1 + {w:.2f}*{c1:.4f}\n"
        f"Numerator: 1 + {w:.2f}*{c1:.4f} = {n1:.4f}; denominator: "
        f"1 - {g:.2f} - {w:.2f}*{m:.2f}*{c1:.4f} = {d1:.4f}\n\n"
        f"**Step 4:** Solve for muG.\n"
        f"muG = {n1:.4f} / {d1:.4f} = {muG:.2f} weeks (rounded to two "
        f"decimals)\n\n"
        f"**Step 5:** Back-substitute to obtain muW.\n"
        f"muW = (1 + {m:.2f}*{muG:.2f}) * {c1:.4f} = {muW:.2f} weeks\n\n"
        f"**Step 6:** Verify with the first equation (non-circular "
        f"check).\n"
        f"1 + {g:.2f}*{muG:.2f} + {w:.2f}*{muW:.2f} = {check:.2f}, which "
        f"matches muG = {muG:.2f} within display rounding.\n\n"
        f"**Answer:** The expected time until failure for a machine "
        f"currently in Good condition is {muG:.2f} weeks"
    )

    return question, solution
