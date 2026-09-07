import random
from data.templates.branches.chemical_engineering.constants import GENERAL_REACTANTS


def _as_printed(x, spec):
    """The value a reader recovers from `x` when it is printed with `spec`.

    P2 asks that the stored value and the printed value be the SAME value, so
    that every downstream line is computed from the number the reader sees
    (Phase 1, DECISIONS D-016).
    """
    return float(format(x, spec))


# Template 1 (Advanced)
def template_levenspiel_plot_interpretation():
    """
    Levenspiel Plot Data Interpretation for Reactor Volume Calculation

    Scenario:
        Experimental data of 1/(-r_A) vs. conversion (X) is provided in tabulated form.
        The goal is to determine the reactor volumes required to achieve a target
        conversion without explicitly knowing the rate law. Two cases are considered:

            - Continuous Stirred-Tank Reactor (CSTR):
            Volume is calculated using the rectangular approximation:
            V_CSTR = F_A0 * X * [1/(-r_A)]_exit

            - Plug Flow Reactor (PFR):
            Volume is calculated using trapezoidal integration of the curve:
            V_PFR = F_A0 * ∫₀ˣ (1/(-r_A)) dX

    Phase 2 notes:

        DETERMINISM. This template sampled its measurement noise from
        `np.random`, which has a generator of its own that `random.seed()` does
        not touch. It was the ONLY T3 failure in the 150-template corpus: a
        recorded seed did not regenerate the item, so a published item could
        not be reproduced, audited or corrected. All sampling is now stdlib
        `random`, and numpy is gone from this file entirely - seeding two
        generators correctly is a thing to get wrong every time, while having
        one is not.

        NO SILENT FALLBACKS. Six blanket `except Exception` handlers each
        substituted a different formula while the trace went on printing the
        original one. Measured over 5000 seeds, every one of them fired on
        0.00% of instances (resolvable to 0.06%), so removing them changes no
        answer - but had one ever fired it would have produced a confidently
        wrong gold trace, silently (P5).

        TRAPEZOIDAL SUM. The answer came from `np.trapezoid` while the trace
        showed a hand-accumulated sum. The two agreed to 9e-16, so nothing was
        wrong - but the answer is now taken from the sum the trace actually
        displays, so the printed derivation IS the calculation (P1).

    Returns:
        tuple: A tuple containing:
            - str: A question asking to compute both CSTR and PFR volumes from tabulated data.
            - str: A step-by-step solution showing the calculations, including trapezoidal rule details.
    """

    # 1. Sample the scenario. stdlib random only.
    reactant_name = random.choice(GENERAL_REACTANTS)
    F_A0 = round(random.uniform(1.0, 5.0), 2)  # mol/s

    target_conversion_options = [0.5, 0.6, 0.65, 0.7, 0.75, 0.8, 0.85, 0.9]
    target_conversion = random.choice(target_conversion_options)

    # 2. Conversion grid, evenly spaced and quoted to 2 dp - the reader sees
    #    these values, so they are what the arithmetic uses.
    if target_conversion <= 0.6:
        num_points = 7
    elif target_conversion <= 0.8:
        num_points = 9
    else:
        num_points = 10

    step = target_conversion / (num_points - 1)
    X_values = [_as_printed(i * step, '.2f') for i in range(num_points)]

    # 3. Tabulated 1/(-r_A), rising with conversion as the kinetics require.
    random_factor = random.uniform(0.8, 1.4)
    inv_rate_values = []
    for x in X_values:
        base = 1.5 + 2.0 * x + 3.0 * x * x
        noise = max(1 + 0.15 * random.uniform(-1, 1), 0.5)
        inv_rate_values.append(_as_printed(base * random_factor * noise, '.2f'))

    # The value at X = 0 anchors the table; keep it readable.
    if inv_rate_values[0] < 1.0:
        inv_rate_values[0] = round(random.uniform(1.0, 2.5), 2)
    inv_rate_values = [max(v, 0.5) for v in inv_rate_values]

    # Enforce the physical expectation that 1/(-r_A) rises with conversion.
    # This is data shaping, not a fallback: it fires on ~74% of draws and is
    # part of how the item is defined.
    for i in range(1, len(inv_rate_values)):
        if inv_rate_values[i] < inv_rate_values[i - 1]:
            inv_rate_values[i] = _as_printed(
                inv_rate_values[i - 1] + round(random.uniform(0.1, 0.3), 2), '.2f')
    inv_rate_values = [min(v, 50.0) for v in inv_rate_values]

    # 4. CSTR volume - the rectangle at the exit condition.
    inv_rate_at_final = inv_rate_values[-1]
    X_final = X_values[-1]
    V_CSTR = F_A0 * inv_rate_at_final * X_final

    # 5. PFR volume - the trapezoidal sum the trace prints, term by term, each
    #    term computed from the operands as displayed.
    # DISPLAY PRECISION IS CHOSEN SO THAT NOTHING ROUNDS. The tabulated rates
    # carry 2 dp, so a half-sum of two of them is exact in 3 dp, and that times
    # a 2-dp interval width is exact in 5 dp. Printing them any shorter would
    # put half of all intervals on a half-way display tie, where no rounding
    # convention closes in both directions (D-016). Measured: at 2 dp the
    # average height ties on 49.97% of intervals and 99.53% of instances carry
    # at least one, so resampling - the Phase 1 remedy - cannot work here.
    trapezoids = []
    for i in range(len(X_values) - 1):
        dx = _as_printed(X_values[i + 1] - X_values[i], '.2f')
        avg_height = _as_printed(
            (inv_rate_values[i] + inv_rate_values[i + 1]) / 2, '.3f')
        trapezoids.append((dx, avg_height, _as_printed(dx * avg_height, '.5f')))
    total_area = _as_printed(sum(t[2] for t in trapezoids), '.5f')
    V_PFR = F_A0 * total_area

    # --- invariants (T7) ---------------------------------------------------
    assert len(X_values) == len(inv_rate_values) == num_points, (
        "conversion grid and rate table are different lengths")
    assert all(b > a for a, b in zip(X_values, X_values[1:])), (
        f"conversion grid is not strictly increasing: {X_values}")
    assert all(v > 0 for v in inv_rate_values), (
        f"non-physical 1/(-r_A) in the table: {inv_rate_values}")
    assert all(b >= a for a, b in zip(inv_rate_values, inv_rate_values[1:])), (
        f"1/(-r_A) must not fall with conversion: {inv_rate_values}")
    assert V_CSTR > 0.0 and V_PFR > 0.0, (
        f"non-physical volumes: CSTR {V_CSTR}, PFR {V_PFR}")
    # Because 1/(-r_A) is non-decreasing, the exit rectangle contains the area
    # under the curve, so the PFR is never the larger reactor. This is a
    # theorem about the sampled data, not a hope - it is why the "unusually
    # large PFR" branch this template used to carry was unreachable.
    assert V_PFR <= V_CSTR + 1e-9, (
        f"PFR volume {V_PFR} exceeds CSTR volume {V_CSTR}, which a "
        f"non-decreasing Levenspiel curve makes impossible")

    # 6. Data table for display
    data_table = "X\t1/(-r_A) [L·s/mol]\n"
    data_table += "-" * 25 + "\n"
    for x, inv_r in zip(X_values, inv_rate_values):
        data_table += f"{x:.2f}\t{inv_r:.2f}\n"

    question = (
        f"A reaction involving {reactant_name} (A → products) is being studied for reactor design. "
        f"Experimental data has been collected and plotted as a Levenspiel plot (1/(-r_A) vs X). "
        f"The inlet molar flow rate is {F_A0} mol/s. Using the data below, calculate:\n\n"
        f"a) The volume of a CSTR required to achieve {target_conversion*100:.0f}% conversion\n"
        f"b) The volume of a PFR required to achieve {target_conversion*100:.0f}% conversion\n\n"
        f"Levenspiel Plot Data:\n"
        f"{data_table}"
    )

    # Steps are numbered contiguously across both parts. They used to restart
    # at Part (b), giving 1,2,3,1,2,3,4,5,6 - which breaks the output contract
    # (T4) and makes "Step 3" ambiguous in a milestone reference.
    solution = (
        f"**Given:**\n"
        f"- Reactant: {reactant_name}\n"
        f"- Inlet molar flow rate: F_A0 = {F_A0} mol/s\n"
        f"- Target conversion: X = {target_conversion} ({target_conversion*100:.0f}%)\n"
        f"- Data points: {len(X_values)} experimental values\n\n"

        f"**Part (a): CSTR Volume Calculation**\n\n"
        f"**Step 1:** For a CSTR, the design equation is:\n"
        f"V_CSTR = F_A0 × X × [1/(-r_A)]_exit\n\n"

        f"**Step 2:** From the data table, at X = {target_conversion}:\n"
        f"[1/(-r_A)]_at_X={target_conversion} = {inv_rate_at_final:.2f} L·s/mol\n\n"

        f"**Step 3:** Calculate CSTR volume:\n"
        f"V_CSTR = {F_A0} mol/s × {target_conversion} × {inv_rate_at_final:.2f} L·s/mol\n"
        f"V_CSTR = {round(V_CSTR, 2)} L\n\n"

        f"**Part (b): PFR Volume Calculation**\n\n"
        f"**Step 4:** For a PFR, the design equation is:\n"
        f"V_PFR = F_A0 × ∫[0 to X] (1/(-r_A)) dX\n\n"

        f"**Step 5:** The integral is the area under the Levenspiel plot curve.\n"
        f"Using the trapezoidal rule on the tabulated points:\n"
        f"Area = Σ[(X_i+1 - X_i) × (y_i + y_i+1)/2]\n"
        f"where y_i = [1/(-r_A)]_i\n\n"

        f"**Step 6:** Calculate individual trapezoid areas:\n"
    )

    for i, (dx, avg_height, trap_area) in enumerate(trapezoids):
        solution += (
            f"Interval [{X_values[i]:.2f} to {X_values[i+1]:.2f}]: "
            f"ΔX = {dx:.2f}, Avg height = ({inv_rate_values[i]:.2f} + "
            f"{inv_rate_values[i+1]:.2f})/2 = {avg_height:.3f}\n"
            f"Area = {dx:.2f} × {avg_height:.3f} = {trap_area:.5f}\n"
        )

    efficiency = _as_printed(((V_CSTR - V_PFR) / V_CSTR) * 100, '.1f')

    solution += (
        f"\n**Step 7:** Total area under curve:\n"
        f"Total area = {total_area:.5f}\n\n"

        f"**Step 8:** Calculate PFR volume:\n"
        f"V_PFR = F_A0 × Area = {F_A0} mol/s × {total_area:.5f}\n"
        f"V_PFR = {round(V_PFR, 2)} L\n\n"

        # The Note sits BEFORE the answer marker, not after it.  It is the only
        # trailing commentary in the corpus and it put a third number (the
        # efficiency percentage) inside the answer span on 100% of instances --
        # a number that is not an answer, in the region a comparator reads as
        # the answer.  Same text, same solution, same pedagogy; it just stops
        # being part of what the answer span means.
        f"**Note:** The PFR requires less volume than the CSTR "
        f"({round(V_PFR, 2)} L vs {round(V_CSTR, 2)} L), a {efficiency:.1f}% "
        f"volume reduction. For a rate that falls with conversion - the usual "
        f"case, and the one this data shows - the plug-flow reactor spends most "
        f"of its length at a higher rate than the CSTR, which operates entirely "
        f"at the exit condition.\n\n"

        f"**Answer:**\n"
        f"a) CSTR Volume = {round(V_CSTR, 2)} L\n"
        f"b) PFR Volume = {round(V_PFR, 2)} L"
    )

    return question, solution


def main():
    """
    Generate numerous instances of the Levenspiel plot template with different random seeds
    and write the results to a JSONL file.
    """

    import json
    import os

    # Define the output path (Modify this path according to where you are running the code from)
    output_file = "testset/chemical_engineering/reaction_kinetics/conversion_and_reactor_sizing.jsonl"

    # Create the directory if it doesn't exist
    os.makedirs(os.path.dirname(output_file), exist_ok=True)
    
    # List of template functions with their ID and level
    templates = [
        (template_levenspiel_plot_interpretation, "levenspiel_plot_interpretation", "Advanced")
    ]

    # List to store all generated problems
    all_problems = []

    # Generate problems for each template in the list
    for template_func, id_name, level in templates:
        for _ in range(50):
            # Generate a unique seed for reproducibility
            seed = random.randint(1_000_000_000, 4_000_000_000)
            random.seed(seed)

            # Generate the question and solution by calling the function
            question, solution = template_func()

            # Create a dictionary entry for the problem
            problem_entry = {
                "seed": seed,
                "branch": "chemical_engineering",
                "domain": "reaction_kinetics",
                "area": "conversion_and_reactor_sizing",
                "id": id_name,
                "level": level,
                "question": question,
                "solution": solution
            }

            # Add the new problem to our list
            all_problems.append(problem_entry)

    # Write all generated problems to a .jsonl file (JSON Lines format)
    with open(output_file, "w") as file:
        for problem in all_problems:
            file.write(json.dumps(problem))
            file.write("\n")

    print(f"\nSuccess! Generated {len(all_problems)} problems and saved them to '{output_file}'")


if __name__ == "__main__":
    main()
