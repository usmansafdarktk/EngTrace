import random
import math
from scipy.integrate import quad
from constants import LIQUID_PHASE_REACTANTS, GENERAL_REACTANTS


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
        reaction. The reaction follows an arbitrary, non-integer order, where the
        rate law is -r_A = k * C_A^n. Because an analytical solution is often
        complex for such cases, the PFR design equation is solved using
        numerical integration:

            V = ∫[0 to X] (F_A0 / (-r_A)) dX

    Returns:
        tuple: A tuple containing:
            - str: A question asking to calculate the required reactor volume.
            - str: A step-by-step solution showing the setup and numerical result.
    """
    
    reactant_name = random.choice(LIQUID_PHASE_REACTANTS)
    F_A0 = round(random.uniform(1.0, 4.0), 2)
    X_final = round(random.uniform(0.5, 0.85), 2)
    
    # Reaction order (non-integer for complexity)
    n = round(random.uniform(1.5, 2.5), 1)
    
    # Rate constant and initial concentration
    k = round(random.uniform(0.05, 0.5), 3)
    C_A0 = round(random.uniform(0.5, 2.0), 2)  # mol/L
    
    def rate_function_advanced(X):
        """Rate as function of conversion: -r_A = k * C_A0^n * (1-X)^n"""
        return k * (C_A0**n) * ((1 - X)**n)
    
    def integrand_advanced(X):
        """Integrand for PFR design equation"""
        return F_A0 / rate_function_advanced(X)
    
    # Numerical integration (analytical solution complex for arbitrary n)
    volume, integration_error = quad(integrand_advanced, 0, X_final)
    
    question = (
        f"An {n}-order liquid-phase reaction of {reactant_name} (A → products) occurs in a PFR. "
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
        f"-r_A = {round(k * (C_A0**n), 4)} × (1 - X)^{n}\n\n"
        
        f"**Step 2:** Set up the PFR integral.\n"
        f"V = ∫[0 to {X_final}] (F_A0 / (-r_A)) dX\n"
        f"V = ∫[0 to {X_final}] ({F_A0} / ({round(k * (C_A0**n), 4)} × (1 - X)^{n})) dX\n\n"
        
        f"**Step 3:** This integral requires numerical methods for n = {n}.\n"
        f"Using numerical integration (scipy.integrate.quad):\n\n"
        
        f"**Step 4:** Numerical result.\n"
        f"V = {round(volume, 2)} L\n"
        f"Integration error: {integration_error:.2e}\n\n"
        
        f"**Note:** For non-integer reaction orders, analytical solutions are complex.\n"
        f"Numerical integration is the standard approach in reactor design.\n\n"
        
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

    # Define the output path
    output_file = "../../../../../testset/reaction_kinetics/mole_balances.jsonl"

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
