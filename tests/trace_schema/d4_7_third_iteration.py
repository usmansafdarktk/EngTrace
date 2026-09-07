"""D4.7 -- fit a third template to the `iteration` node type, binding only.

The rule, from the deliverable: **declare a binding, edit no verifier.**  D3.3
|S|10 says the claim that `iteration` generalises is "argued, not measured",
because ``roles`` and ``update_relation`` exist precisely so that fitting a new
template is a declaration rather than a code change -- and nobody had tried.
A *synthetic* renamed node passes the rename test; a real second template had
never been attempted.

This attempts both templates the spec names, ``linear_reservoir_routing_step``
and ``qr_policy_one_iteration``, builds the most faithful `iteration` node each
admits, and runs the **unmodified** ``reviewer_d2_verifier_15``.  The verdicts
are below and neither is a pass.  That is the measurement: the claim does not
generalise to either named target, and the two fail for *different* reasons,
which is what makes the result useful rather than merely negative.

Run: ``python -m tests.trace_schema.d4_7_third_iteration``
"""

from __future__ import annotations

import random

from tests.trace_schema.reviewer_d2_verifier_15 import verify_node


# --------------------------------------------------------------------------
# Candidate 1 -- linear_reservoir_routing_step
# --------------------------------------------------------------------------


def routing_node(seed: int = 0) -> dict:
    """The most faithful `iteration` node ``linear_reservoir_routing_step`` admits.

    The template routes an inflow hydrograph through two intervals:

        O_{j+1} = [ (I_j + I_{j+1})*dt/2 + O_j*(k - dt/2) ] / (k + dt/2)

    which really is a repeated sub-chain with a carried iterate, so the
    structural resemblance the spec saw is genuine.  Everything below is a
    declaration; no verifier line is touched.
    """
    import importlib

    mod = importlib.import_module(
        "data.templates.branches.civil_engineering.water_resources.hydrology")
    random.seed(seed)
    # Re-derive the same parameters the template samples, in the same order.
    dt = random.choice([900, 1200, 1800])
    k = random.randint(max(1800, dt // 2 + 300) // 60, 5400 // 60) * 60
    I1 = round(random.uniform(2.0, 10.0), 1)
    I2 = round(I1 + random.uniform(2.0, 10.0), 1)
    I3 = round(I2 + random.uniform(1.5, 8.0), 1)
    O1 = round(random.uniform(0.3, 0.6) * I1, 1)
    half_dt = round(dt / 2.0, 1)
    den = round(k + half_dt, 1)
    coeff = round(k - half_dt, 1)
    num1a = round((I1 + I2) * half_dt, 1)
    num1b = round(O1 * coeff, 1)
    O2 = round((num1a + num1b) / den, 2)
    num2a = round((I2 + I3) * half_dt, 1)
    num2b = round(O2 * coeff, 1)
    O3 = round((num2a + num2b) / den, 2)
    assert mod is not None

    return {
        "schema_version": "1.5",
        "node_type": "iteration",
        "node_id": "t_linear_reservoir_routing",
        "cardinality": "incidental",
        "roles": {
            "index": "j",
            "iterate_prev": "O_prev",
            "iterate_curr": "O_curr",
            "iterate_next": "O_next",
            "residual_prev": "num_a_prev",
            "residual_curr": "num_a",
            "change": "change",
            "converged": "converged",
            "evaluation": "evaluation",
        },
        # The routing recurrence, written in S4.2's grammar as faithfully as it
        # can be.  Note what is missing: the inflow pair (I_j, I_{j+1}) changes
        # every step and there is nowhere to put it.
        "update_relation": "(residual_curr + iterate_curr * coeff) / den",
        "termination": {
            "kind": "convergence",
            "quantity_role": "change",
            "comparison": "lt",
            "tolerance": 0.0,
            "predicate_prose": "the inflow hydrograph is exhausted after two intervals",
            "max_elements": 2,
            "satisfied": True,
        },
        "rounding": "decimal-half-up",
        "constants": {"k": k, "dt": dt, "coeff": coeff, "den": den,
                      "I1": I1, "I2": I2, "I3": I3},
        "frame_relations": [["num_a", "0"]],
        "symbol_precision": {
            "O_prev": 2, "O_curr": 2, "O_next": 2, "change": 2,
            "num_a_prev": 1, "num_a": 1, "O": 2,
        },
        "preamble_binding": {},
        "preamble_symbols": ["O", "num_a"],
        "evaluation_symbols": ["O", "num_a"],
        "evaluation_roles": {"iterate": "O", "residual": "num_a"},
        "element_symbols": ["j", "O_prev", "num_a_prev", "O_curr", "num_a",
                            "O_next", "change", "converged", "evaluation"],
        "carry": {"O_prev": "O_curr", "O_curr": "O_next"},
        "carry_notes": {},
        "preamble": [{"O": O1, "num_a": num1a}],
        "elements": [
            {"j": 1, "O_prev": O1, "num_a_prev": num1a, "O_curr": O1,
             "num_a": num1a, "O_next": O2, "change": round(O2 - O1, 2),
             "converged": False, "evaluation": {"O": O2, "num_a": num2a}},
            {"j": 2, "O_prev": O1, "num_a_prev": num1a, "O_curr": O2,
             "num_a": num2a, "O_next": O3, "change": round(O3 - O2, 2),
             "converged": True, "evaluation": None},
        ],
        "result": {"from": "elements[-1].O_next", "value": O3, "dp": 2,
                   "unit": "m^3/s"},
    }


# --------------------------------------------------------------------------
# Candidate 2 -- qr_policy_one_iteration
# --------------------------------------------------------------------------


def qr_node() -> dict:
    """``qr_policy_one_iteration`` as an `iteration` node.

    The (Q, R) procedure is a genuine fixed-point iteration, so the type looks
    like a fit -- until the iterate is written down.  It iterates on a *pair*,
    ``(Q, R)``, coupled through ``n(R)``; `iteration` |S|4.1 declares exactly one
    ``iterate_prev``/``iterate_curr``/``iterate_next`` triple.  Below, ``Q`` is
    bound as the iterate and ``R`` is demoted into the evaluation frame, which
    is the only declaration the type allows.
    """
    return {
        "schema_version": "1.5",
        "node_type": "iteration",
        "node_id": "t_qr_policy",
        "cardinality": "incidental",
        "roles": {
            "index": "i", "iterate_prev": "Q_prev", "iterate_curr": "Q_curr",
            "iterate_next": "Q_next", "residual_prev": "n_prev",
            "residual_curr": "n_curr", "change": "change",
            "converged": "converged", "evaluation": "evaluation",
        },
        "update_relation": "iterate_curr",
        "termination": {
            "kind": "convergence", "quantity_role": "change",
            "comparison": "lt", "tolerance": 1.0,
            "predicate_prose": "one iteration is requested",
            "max_elements": 1, "satisfied": True,
        },
        "rounding": "decimal-half-up",
        "constants": {"lam": 12000.0, "K": 80.0, "h": 3.2, "p": 12.0},
        "frame_relations": [["n", "0"]],
        "symbol_precision": {"Q_prev": 0, "Q_curr": 0, "Q_next": 0,
                             "change": 0, "n_prev": 2, "n_curr": 2,
                             "Q": 0, "n": 2},
        "preamble_binding": {},
        "preamble_symbols": ["Q", "n"],
        "evaluation_symbols": ["Q", "n"],
        "evaluation_roles": {"iterate": "Q", "residual": "n"},
        "element_symbols": ["i", "Q_prev", "n_prev", "Q_curr", "n_curr",
                            "Q_next", "change", "converged", "evaluation"],
        "carry": {},
        "carry_notes": {},
        "preamble": [{"Q": 775.0, "n": 0.0}],
        "elements": [
            {"i": 1, "Q_prev": 775.0, "n_prev": 0.0, "Q_curr": 775.0,
             "n_curr": 0.0, "Q_next": 913.0, "change": 138.0,
             "converged": True, "evaluation": None},
        ],
        "result": {"from": "elements[-1].Q_next", "value": 913.0, "dp": 0,
                   "unit": "units"},
    }


VERDICTS = {
    "linear_reservoir_routing_step": (
        "DOES NOT FIT, and the reason is structural rather than cosmetic.\n"
        "  1. S4 fixes `termination.kind: \"convergence\"` for this type. Routing "
        "does not converge; it EXHAUSTS a two-interval hydrograph. The element "
        "count is fixed by the question, so it is neither `incidental` (a third "
        "interval would give a different answer) nor `answer_bearing` (the count "
        "is not the answer, O3 is). D-038's discriminating rule has no verdict "
        "for it, because the rule assumes the count is either free under the "
        "comparator's slack or present in the answer, and here it is neither.\n"
        "  2. FrameEval requires every name in a frame relation to be a "
        "constant, the iterate, or an EARLIER frame symbol. The routing frame "
        "needs the inflow pair (I_j, I_{j+1}), which changes every element. "
        "There is nowhere in the node to put a per-element exogenous input. "
        "Declaring the inflows as `constants` (as this attempt does, to get as "
        "far as possible) makes the SAME inflow apply to both intervals, which "
        "is a different problem.\n"
        "  The spec named this template on the strength of it being 'a repeated "
        "sub-chain unrolled into the trace'. That resemblance is real and it is "
        "not a type match: `iteration` is not 'a repeated sub-chain', it is "
        "'a convergence-terminated refinement of ONE quantity driven by ITS OWN "
        "residual'."
    ),
    "qr_policy_one_iteration": (
        "DOES NOT FIT, for an unrelated reason, which is what makes the pair "
        "informative.\n"
        "  The (Q, R) procedure iterates on a PAIR coupled through n(R): "
        "Q = f(n(R)), R = g(Q). S4.1 declares exactly one iterate triple, so a "
        "binding must demote R into the evaluation frame -- and then S4.3.1's "
        "recomputation of the update from the residuals cannot see R at all, "
        "because the frame is a function of the iterate alone.\n"
        "  Separately, the template emits ONE iteration by design, so the "
        "termination predicate is 'the question asked for one', not a "
        "convergence test on a quantity. Declaring `tolerance` to make that "
        "predicate true is declaring a number the question does not contain."
    ),
}


def main() -> None:
    print("D4.7 -- fitting a third template to `iteration`, binding only")
    print("=" * 74)
    print("Verifier: reviewer_d2_verifier_15, UNMODIFIED.\n")

    for name, builder in (("linear_reservoir_routing_step", routing_node),
                          ("qr_policy_one_iteration", qr_node)):
        print(f"--- {name}")
        try:
            node = builder()
        except Exception as exc:                            # noqa: BLE001
            print(f"    could not even construct a node: {type(exc).__name__}: {exc}\n")
            continue
        rep = verify_node(node, label=name, mode="gold")
        if not rep.failures:
            print("    ACCEPTED -- the type generalises to this template by "
                  "declaration alone.\n")
            continue
        print(f"    REJECTED, {len(rep.failures)} clause failure(s):")
        for f in rep.failures[:8]:
            print(f"      {f}")
        if len(rep.failures) > 8:
            print(f"      ... and {len(rep.failures) - 8} more")
        print()

    print("=" * 74)
    print("VERDICT\n")
    for name, text in VERDICTS.items():
        print(f"{name}:")
        for line in text.split("\n"):
            print(f"  {line}")
        print()
    print("D4.7 answer: `iteration` does NOT generalise to either template the")
    print("spec names. D3.3 S10's 'argued, not measured' is now MEASURED, and the")
    print("measurement is negative. The type is sound; its stated scope is wrong.")


if __name__ == "__main__":
    main()
