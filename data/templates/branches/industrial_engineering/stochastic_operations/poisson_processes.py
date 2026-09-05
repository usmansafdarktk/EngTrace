import math
import random
from decimal import Decimal, ROUND_HALF_UP

# Poisson-process template. No numeric constants from constants.py by
# Stage B design (rates sampled and fully stated; given-values rule).
# Typology anchor: Ross 11e Ch. 5 (Sec. 5.3, the Poisson process and the
# Poisson distribution of N(t)); cross-ref H&L 7e Sec. 17.4.


def _hu4(x):
    """Half-up 4-dp rounding of a float via its shortest decimal repr."""
    return float(Decimal(repr(x)).quantize(Decimal("0.0001"),
                                           rounding=ROUND_HALF_UP))


# Reachable (rate_tenths, weeks, k) combos, built at import:
#   mu = lam * t <= 4.0 (keeps e^-mu >= 0.018 so the displayed 5-dp e^-mu
#   drives the chain within 0.03% of full precision; lessons 5/30);
#   exact P(N(t) = k) >= 0.06 (keeps the 4-dp display informative); and
#   the displayed chain (P from the 5-dp e^-mu) rounds HALF-UP to the SAME
#   4-dp value as the full-precision pmf, so a full-precision solve matches
#   the gold answer digit-for-digit (lesson 51).
_T10_COMBOS = []
for _lc in range(5, 21):                 # lam = 0.5 .. 2.0 per week, 1 dp
    for _t in (1, 2, 3, 4):
        _lam = _lc / 10
        _mu = round(_lam * _t, 1)        # exact 1-dp product
        if _mu > 4.0:
            continue
        _e5 = round(math.exp(-_mu), 5)
        for _k in range(0, 7):
            _p_exact = math.exp(-_mu) * _mu ** _k / math.factorial(_k)
            if _p_exact < 0.06:
                continue
            _p_chain = _e5 * _mu ** _k / math.factorial(_k)
            if _hu4(_p_chain) == _hu4(_p_exact):
                _T10_COMBOS.append((_lc, _t, _k))
assert len(_T10_COMBOS) >= 150


# Template 10 (Easy) — Area S4: Poisson Process & Exponential Distribution
def template_poisson_event_count():
    """
    Poisson Process: Probability of an Exact Event Count

    Scenario:
        Machine breakdowns occur according to a Poisson process at a
        known average rate per week. Over an observation window of t
        weeks the number of breakdowns N(t) is Poisson with mean
        mu = lambda * t:

            P(N(t) = k) = e^(-mu) * mu^k / k!

        Requested: the probability of exactly k breakdowns in the
        window.

    Difficulty: Easy
    Grounding: Ross, Introduction to Probability Models, 11th ed., Ch. 5,
        Sec. 5.3 (the Poisson process; N(t) ~ Poisson(lambda*t)).
        Cross-ref Hillier & Lieberman 7e, Sec. 17.4 (role of the
        exponential/Poisson in queueing inputs).
    Physical bounds: lambda in [0.5, 2.0] per week (1 dp), t in
        [1, 4] weeks, with mu = lambda*t <= 4.0; k in [0, 6] restricted
        to combos with exact P >= 0.06 AND displayed-chain/full-precision
        4-dp half-up agreement (builder-enforced; 212 reachable combos at
        import, 19 candidates rejected by the agreement filter; count
        asserted >= 150). Reachable displayed P(N=k) lies in
        [0.0602, 0.6065] (author QA 2026-08-06, exhaustive); assert P in
        [0.055, 0.61].

    Returns:
        tuple(str, str): (question, solution)
    """
    lc, t, k = random.choice(_T10_COMBOS)
    lam = lc / 10
    mu = round(lam * t, 1)               # exact 1-dp product

    e5 = round(math.exp(-mu), 5)
    P = _hu4(e5 * mu ** k / math.factorial(k))

    assert 0.5 <= lam <= 2.0 and 1 <= t <= 4, f"(lam,t) out of bounds: {lam},{t}"
    assert mu <= 4.0, f"mu out of bounds: {mu}"
    assert 0.055 <= P <= 0.61, f"P out of bounds: {P}"

    week_word = "week" if t == 1 else f"{t} weeks"
    question = (
        f"Breakdowns of a stamping press occur according to a Poisson "
        f"process at an average rate of {lam:.1f} breakdowns per week. "
        f"The plant schedules maintenance reviews every {week_word}. "
        f"Determine the probability that exactly {k} breakdown"
        f"{'s' if k != 1 else ''} occur{'s' if k == 1 else ''} during a "
        f"single review window, to four decimal places. In your "
        f"solution, first give the expected number of breakdowns in the "
        f"window."
    )

    solution = (
        f"**Given:**\n"
        f"Poisson breakdown process with rate lambda = {lam:.1f} per "
        f"week; window length t = {t} week{'s' if t != 1 else ''}; "
        f"requested count k = {k}.\n\n"
        f"**Step 1:** Compute the expected number of breakdowns in the "
        f"window. For a Poisson process, N(t) is Poisson with mean\n"
        f"mu = lambda * t = {lam:.1f} * {t} = {mu:.1f} breakdowns\n\n"
        f"**Step 2:** Evaluate the exponential factor of the Poisson "
        f"probability mass function P(N = k) = e^(-mu) * mu^k / k!.\n"
        f"e^(-mu) = e^(-{mu:.1f}) = {e5:.5f}\n\n"
        f"**Step 3:** Evaluate the probability of exactly {k} "
        f"breakdown{'s' if k != 1 else ''}.\n"
        f"P(N = {k}) = {e5:.5f} * ({mu:.1f})^{k} / {k}! = {P:.4f} "
        f"(rounded to four decimals)\n\n"
        f"**Answer:** The probability of that exact breakdown count in "
        f"the window is {P:.4f}"
    )

    return question, solution
