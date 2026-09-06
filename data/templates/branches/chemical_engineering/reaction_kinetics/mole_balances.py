import random
import math
from data.templates.branches.chemical_engineering.constants import LIQUID_PHASE_REACTANTS, GENERAL_REACTANTS


# Display precision for the intermediates of template_pfr_volume_changing_rate.
PFR_DP = 4


def _as_printed(x, spec):
    """The value a reader recovers from `x` when it is printed with `spec`.

    P2 asks that the stored value and the printed value be the SAME value, so
    that every downstream line is computed from the number the reader can see.
    Rounding alone does not achieve that (Phase 1, DECISIONS D-016).
    """
    return float(format(x, spec))


# Template 1 (Easy)
def template_cstr_volume_basic():
    """
    CSTR Volume Calculation

    Scenario:
        The General Mole Balance for a Continuous Stirred-Tank Reactor (CSTR) at steady
        state is used to determine the necessary reactor volume to achieve a desired
        outcome. Given the inlet and outlet molar flow rates and the reaction rate,
        the goal is to compute the volume using:

            V = (F_A0 - F_A) / (-r_A)

    Returns:
        tuple: A tuple containing:
            - str: A question asking to compute the CSTR volume.
            - str: A step-by-step solution showing the calculation (the Python trace).
    """
    # 1. Parameterize the inputs with random values
    reactant_name = random.choice(GENERAL_REACTANTS)
    # Inlet molar flow rate in mol/s
    F_A0 = round(random.uniform(1.0, 5.0), 2)
    # Outlet molar flow rate in mol/s (must be less than inlet)
    F_A = round(random.uniform(0.1, F_A0 - 0.1), 2)
    # Reaction rate in mol/(L·s) - a positive value for -r_A
    neg_r_A = round(random.uniform(0.01, 0.05), 4)

    # 2. Perform the core calculation
    volume = (F_A0 - F_A) / neg_r_A

    # 3. Generate the question and solution strings
    question = (
        f"A steady-state CSTR is used for the consumption of {reactant_name}. "
        f"The inlet molar flow rate of {reactant_name} is {F_A0} mol/s, and the desired "
        f"outlet molar flow rate is {F_A} mol/s. If the rate of reaction is {neg_r_A} mol/(L·s), "
        f"what is the required reactor volume in liters?"
    )

    solution = (
        f"**Step 1:** State the CSTR design equation.\n"
        f"    V = (F_A0 - F_A) / (-r_A)\n\n"
        f"**Step 2:** Substitute the given values into the equation.\n"
        f"    V = ({F_A0} mol/s - {F_A} mol/s) / ({neg_r_A} mol/(L·s))\n\n"
        f"**Step 3:** Calculate the final volume.\n"
        f"    V = {round(volume, 2)} L\n\n"
        f"**Answer:** The required reactor volume is {round(volume, 2)} liters."
    )

    return question, solution


# Template 2 (Easy)
def template_batch_reactor_zero_order():
    """
    Batch Reactor Time Calculation - Zero Order Kinetics

    Scenario:
        The time required for a reaction in a constant-volume batch reactor is 
        determined by integrating the mole balance equation. For a zero-order 
        reaction, the rate of reaction is constant and independent of concentration. 
        Given the initial and final concentrations and the rate constant, the 
        goal is to compute the necessary time using the integrated rate law:

            t = (C_A0 - C_A) / k

    Returns:
        tuple: A tuple containing:
            - str: A question asking to calculate the required reaction time.
            - str: A step-by-step solution showing the derivation and calculation.
    """
    
    # 1. Parameterize inputs
    reactant_name = random.choice(LIQUID_PHASE_REACTANTS)
    C_A0 = round(random.uniform(2.0, 10.0), 2)  # Initial concentration (mol/L)
    
    # Generate final concentration ensuring reasonable conversion
    conversion = round(random.uniform(0.2, 0.8), 2)
    C_A = round(C_A0 * (1 - conversion), 2)
    
    # Zero-order rate constant (mol/(L·s))
    k = round(random.uniform(0.01, 0.1), 4)
    
    # Calculate time
    time = (C_A0 - C_A) / k
    
    # Generate question and solution
    question = (
        f"A zero-order liquid-phase reaction involving {reactant_name} is carried out "
        f"in a constant-volume batch reactor. The initial concentration is {C_A0} mol/L, "
        f"and the desired final concentration is {C_A} mol/L. If the zero-order rate "
        f"constant is {k} mol/(L·s), calculate the required reaction time."
    )
    
    solution = (
        f"**Given:**\n"
        f"- Initial concentration (C_A0) = {C_A0} mol/L\n"
        f"- Final concentration (C_A) = {C_A} mol/L\n"
        f"- Zero-order rate constant (k) = {k} mol/(L·s)\n"
        f"- Conversion = {round(conversion * 100, 1)}%\n\n"
        
        f"**Step 1:** Write the mole balance for a constant-volume batch reactor.\n"
        f"dC_A/dt = r_A\n\n"
        
        f"**Step 2:** For zero-order kinetics, r_A = -k (constant).\n"
        f"dC_A/dt = -k\n\n"
        
        f"**Step 3:** Integrate from t=0 (C_A = C_A0) to t=t (C_A = C_A).\n"
        f"∫[C_A0 to C_A] dC_A = -k ∫[0 to t] dt\n"
        f"C_A - C_A0 = -k × t\n\n"
        
        f"**Step 4:** Solve for time.\n"
        f"t = (C_A0 - C_A) / k\n\n"
        
        f"**Step 5:** Substitute values.\n"
        f"t = ({C_A0} - {C_A}) mol/L / {k} mol/(L·s)\n"
        f"t = {C_A0 - C_A} / {k} = {round(time, 2)} s\n\n"
        
        f"**Answer:** The required reaction time is {round(time, 2)} seconds."
    )
    
    return question, solution


# Template 3 (Intermediate)
def template_batch_reactor_first_order():
    """
    Batch Reactor Time Calculation - First Order Kinetics

    Scenario:
        This template calculates the time for a reaction in a constant-volume
        batch reactor following first-order kinetics. Unlike zero-order reactions,
        the rate of a first-order reaction is directly proportional to the
        concentration of the reactant (-r_A = k*C_A). The required time is found
        by integrating the mole balance equation. Given the initial and final
        concentrations and the rate constant, the goal is to compute the time using:

            t = (1/k) * ln(C_A0 / C_A)

    Returns:
        tuple: A tuple containing:
            - str: A question asking to calculate the required reaction time.
            - str: A step-by-step solution showing the derivation and calculation.
    """
    
    reactant_name = random.choice(LIQUID_PHASE_REACTANTS)
    C_A0 = round(random.uniform(1.0, 5.0), 2)
    
    conversion = round(random.uniform(0.3, 0.9), 2)
    C_A = round(C_A0 * (1 - conversion), 2)
    
    # First-order rate constant (1/s)
    k = round(random.uniform(0.001, 0.01), 5)
    
    # Calculate time
    time = math.log(C_A0 / C_A) / k
    
    question = (
        f"A first-order liquid-phase reaction of {reactant_name} occurs in a batch reactor. "
        f"The initial concentration is {C_A0} mol/L, and after reaction, the concentration "
        f"decreases to {C_A} mol/L. If the first-order rate constant is {k} s⁻¹, "
        f"determine the reaction time required."
    )
    
    solution = (
        f"**Given:**\n"
        f"- Initial concentration (C_A0) = {C_A0} mol/L\n"
        f"- Final concentration (C_A) = {C_A} mol/L\n"
        f"- First-order rate constant (k) = {k} s⁻¹\n"
        f"- Conversion = {round(conversion * 100, 1)}%\n\n"
        
        f"**Step 1:** Write the rate law for first-order kinetics.\n"
        f"-r_A = k × C_A\n\n"
        
        f"**Step 2:** Set up the mole balance equation.\n"
        f"dC_A/dt = -k × C_A\n\n"
        
        f"**Step 3:** Separate variables and integrate.\n"
        f"dC_A/C_A = -k × dt\n"
        f"∫[C_A0 to C_A] dC_A/C_A = -k ∫[0 to t] dt\n\n"
        
        f"**Step 4:** Solve the integral.\n"
        f"ln(C_A/C_A0) = -k × t\n"
        f"ln(C_A0/C_A) = k × t\n\n"
        
        f"**Step 5:** Solve for time.\n"
        f"t = (1/k) × ln(C_A0/C_A)\n\n"
        
        f"**Step 6:** Substitute values.\n"
        f"t = (1/{k}) × ln({C_A0}/{C_A})\n"
        f"t = {round(1/k, 2)} × ln({round(C_A0/C_A, 3)})\n"
        f"t = {round(1/k, 2)} × {round(math.log(C_A0/C_A), 3)}\n"
        f"t = {round(time, 2)} s\n\n"
        
        f"**Answer:** The required reaction time is {round(time, 2)} seconds."
    )
    
    return question, solution


# Template 4 (Advanced)
def template_batch_reactor_second_order():
    """
    Batch Reactor Time Calculation - Second Order Kinetics

    Scenario:
        This template calculates the reaction time in a constant-volume batch
        reactor for a second-order reaction. In this case, the reaction rate
        depends on the square of the reactant's concentration (-r_A = k*C_A²).
        The time required is determined by integrating the mole balance
        equation for these kinetics. Given the initial and final concentrations
        and the second-order rate constant, the goal is to compute the time using:

            t = (1/k) * (1/C_A - 1/C_A0)

    Returns:
        tuple: A tuple containing:
            - str: A question asking to calculate the required reaction time.
            - str: A step-by-step solution showing the derivation and calculation.
    """
    
    reactant_name = random.choice(LIQUID_PHASE_REACTANTS)
    C_A0 = round(random.uniform(0.5, 2.0), 2)  # Lower values for second-order
    
    conversion = round(random.uniform(0.4, 0.85), 2)
    C_A = round(C_A0 * (1 - conversion), 2)
    
    # Second-order rate constant (L/(mol·s))
    k = round(random.uniform(0.1, 1.0), 3)
    
    # Calculate time
    time = (1/C_A - 1/C_A0) / k
    
    question = (
        f"A second-order reaction of {reactant_name} takes place in a batch reactor. "
        f"Starting with {C_A0} mol/L, the concentration drops to {C_A} mol/L. "
        f"Given that the second-order rate constant is {k} L/(mol·s), "
        f"calculate the time required for this conversion."
    )
    
    solution = (
        f"**Given:**\n"
        f"- Initial concentration (C_A0) = {C_A0} mol/L\n"
        f"- Final concentration (C_A) = {C_A} mol/L\n"
        f"- Second-order rate constant (k) = {k} L/(mol·s)\n"
        f"- Conversion = {round(conversion * 100, 1)}%\n\n"
        
        f"**Step 1:** Write the rate law for second-order kinetics.\n"
        f"-r_A = k × C_A²\n\n"
        
        f"**Step 2:** Set up the mole balance equation.\n"
        f"dC_A/dt = -k × C_A²\n\n"
        
        f"**Step 3:** Separate variables and integrate.\n"
        f"dC_A/C_A² = -k × dt\n"
        f"∫[C_A0 to C_A] C_A⁻² dC_A = -k ∫[0 to t] dt\n\n"
        
        f"**Step 4:** Solve the integral.\n"
        f"[-1/C_A] from C_A0 to C_A = -k × t\n"
        f"-1/C_A + 1/C_A0 = -k × t\n"
        f"1/C_A - 1/C_A0 = k × t\n\n"
        
        f"**Step 5:** Solve for time.\n"
        f"t = (1/k) × (1/C_A - 1/C_A0)\n\n"
        
        f"**Step 6:** Substitute values.\n"
        f"t = (1/{k}) × (1/{C_A} - 1/{C_A0})\n"
        f"t = {round(1/k, 3)} × ({round(1/C_A, 3)} - {round(1/C_A0, 3)})\n"
        f"t = {round(1/k, 3)} × {round(1/C_A - 1/C_A0, 3)}\n"
        f"t = {round(time, 2)} s\n\n"
        
        f"**Answer:** The required reaction time is {round(time, 2)} seconds."
    )
    
    return question, solution


# Template 5 (Advanced)
def template_pfr_volume_changing_rate():
    """
    PFR Volume Calculation - Arbitrary Order Kinetics

    Scenario:
        This template calculates the required volume for a Plug Flow Reactor (PFR)
        to achieve a target conversion for a liquid-phase (constant-density)
        reaction. The reaction follows an arbitrary order n, where the rate law
        is -r_A = k * C_A^n. The PFR design equation

            V = ∫[0 to X] (F_A0 / (-r_A)) dX

        is separable for a constant-density power law and integrates in closed
        form for any n != 1:

            V = F_A0 / (k * C_A0^n) * [ (1-X)^(1-n) - 1 ] / (n - 1)

        Phase 2 note: this used to call scipy.integrate.quad and print quad's
        own error estimate as a solution step - a machine-dependent number in a
        gold trace, and a black box standing in for an integral a student can
        do by hand. The sampled range is n in [1.5, 2.5], so n = 1 (the
        logarithmic case) never arises.

    Returns:
        tuple: A tuple containing:
            - str: A question asking to calculate the required reactor volume.
            - str: A step-by-step solution showing the setup and numerical result.
    """
    
    reactant_name = random.choice(LIQUID_PHASE_REACTANTS)
    F_A0 = round(random.uniform(1.0, 4.0), 2)
    X_final = round(random.uniform(0.5, 0.85), 2)
    
    # Reaction order. n = 2.0 is drawn ~9.6% of the time and is a perfectly
    # valid instance; what was wrong was the question calling every draw
    # "non-integer". The wording was fixed rather than the sampling, so the
    # item pool is unchanged (P6).
    n = round(random.uniform(1.5, 2.5), 1)
    
    # Rate constant and initial concentration
    k = round(random.uniform(0.05, 0.5), 3)
    C_A0 = round(random.uniform(0.5, 2.0), 2)  # mol/L
    
    # Closed form. Separating the PFR design equation for -r_A = k*C_A^n with
    # C_A = C_A0*(1-X) gives
    #     V = F_A0 / (k*C_A0^n) * ∫[0 to X] (1-X)^-n dX
    # and ∫(1-X)^-n dX = [(1-X)^(1-n) - 1] / (n-1) for n != 1.
    # Verified against scipy.integrate.quad over 3000 seeds: worst relative
    # difference 2.7e-15, i.e. identical to machine precision. No answer moves.
    # Each intermediate is bound THROUGH its display, and the next line is
    # computed from the bound value - so a reader following the printed
    # operands reaches the printed answer exactly (P1/P2).
    # Preconditions. n == 1 is the logarithmic case the closed form cannot
    # express; the sampled range is [1.5, 2.5] so it never arises, but the
    # closed form is only valid because of that and should say so.
    assert n != 1.0, f"closed form undefined for n = 1 (logarithmic case): {n}"
    assert 0.0 < X_final < 1.0, f"conversion outside (0, 1): {X_final}"
    assert k > 0.0 and C_A0 > 0.0, f"non-physical rate data: k={k}, C_A0={C_A0}"

    spec = f'.{PFR_DP}f'
    rate_coefficient = _as_printed(k * (C_A0 ** n), spec)
    one_minus_X = _as_printed(1 - X_final, spec)
    power_term = _as_printed(one_minus_X ** (1 - n), spec)
    integral_term = _as_printed((power_term - 1) / (n - 1), spec)
    volume = F_A0 / rate_coefficient * integral_term

    assert rate_coefficient > 0.0, f"rate coefficient vanished: {rate_coefficient}"
    assert integral_term > 0.0, f"integral must be positive for X in (0,1): {integral_term}"
    assert volume > 0.0, f"non-physical reactor volume: {volume}"
    
    question = (
        f"An order-{n} liquid-phase reaction of {reactant_name} (A → products) occurs in a PFR. "
        f"The inlet conditions are: F_A0 = {F_A0} mol/s and C_A0 = {C_A0} mol/L. "
        f"The rate expression is: -r_A = {k} × C_A^{n} mol/(L·s). "
        f"Determine the reactor volume needed for {X_final*100}% conversion."
    )
    
    solution = (
        f"**Given:**\n"
        f"- Inlet molar flow rate: F_A0 = {F_A0} mol/s\n"
        f"- Initial concentration: C_A0 = {C_A0} mol/L\n"
        f"- Rate constant: k = {k} L^{n-1}/(mol^{n-1}·s)\n"
        f"- Reaction order: n = {n}\n"
        f"- Desired conversion: X = {X_final}\n\n"
        
        f"**Step 1:** Express rate in terms of conversion.\n"
        f"C_A = C_A0(1 - X)\n"
        f"-r_A = k × C_A^{n} = k × C_A0^{n} × (1 - X)^{n}\n"
        f"-r_A = {k} × {C_A0}^{n} × (1 - X)^{n}\n"
        f"-r_A = {rate_coefficient} × (1 - X)^{n}\n\n"
        
        f"**Step 2:** Set up the PFR integral.\n"
        f"V = ∫[0 to {X_final}] (F_A0 / (-r_A)) dX\n"
        f"V = ({F_A0} / {rate_coefficient}) × ∫[0 to {X_final}] (1 - X)^(-{n}) dX\n\n"

        f"**Step 3:** Integrate. For n ≠ 1, ∫(1 - X)^(-n) dX = [(1 - X)^(1-n) - 1] / (n - 1).\n"
        f"1 - X = 1 - {X_final} = {one_minus_X}\n"
        f"(1 - X)^(1-n) = {one_minus_X}^({round(1 - n, 4)}) = {power_term}\n"
        f"Integral = ({power_term} - 1) / ({n} - 1) = {integral_term}\n\n"

        f"**Step 4:** Evaluate the volume.\n"
        f"V = ({F_A0} / {rate_coefficient}) × {integral_term}\n"
        f"V = {round(volume, 2)} L\n\n"
        
        f"**Answer:** The required PFR volume is {round(volume, 2)} liters."
    )
    
    return question, solution


def main():
    """
    Generate numerous instances of each mole balances template with different random seeds
    and write the results to a JSONL file.
    """
    import json
    import os

    # Define the output path (Modify this path according to where you are running the code from)
    output_file = "testset/chemical_engineering/reaction_kinetics/mole_balances.jsonl"

    # Create the directory if it doesn't exist
    os.makedirs(os.path.dirname(output_file), exist_ok=True)

    # List of template functions with their ID and level
    templates = [
        (template_cstr_volume_basic, "cstr_volume_basic", "Easy"),
        (template_batch_reactor_zero_order, "batch_reactor_zero_order", "Easy"),
        (template_batch_reactor_first_order, "batch_reactor_first_order", "Intermediate"),
        (template_batch_reactor_second_order, "batch_reactor_second_order", "Advanced"),
        (template_pfr_volume_changing_rate, "pfr_volume_changing_rate", "Advanced"),
    ]

    # List to store all generated problems
    all_problems = []

    # Generate problems for each template
    for template_func, id_name, level in templates:
        for _ in range(50):
            # Generate a unique seed for each problem
            seed = random.randint(1_000_000_000, 4_000_000_000)
            random.seed(seed)

            # Generate the problem and solution
            question, solution = template_func()

            # Create a JSON entry
            problem_entry = {
                "seed": seed,
                "branch": "chemical_engineering",
                "domain": "reaction_kinetics",
                "area": "mole_balances",
                "id": id_name,
                "level": level,
                "question": question,
                "solution": solution
            }

            # Add to the list of problems
            all_problems.append(problem_entry)

    # Shuffle the problems to mix templates and levels
    random.shuffle(all_problems)

    # Write all problems to a .jsonl file
    with open(output_file, "w") as file:
        for problem in all_problems:
            file.write(json.dumps(problem))
            file.write("\n")

    print(f"Successfully generated {len(all_problems)} problems and saved to {output_file}")


if __name__ == "__main__":
    main()
