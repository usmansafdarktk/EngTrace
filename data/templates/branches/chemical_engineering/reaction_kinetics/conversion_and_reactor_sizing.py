import random
import numpy as np
from data.templates.branches.chemical_engineering.constants import GENERAL_REACTANTS


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

    Returns:
        tuple: A tuple containing:
            - str: A question asking to compute both CSTR and PFR volumes from tabulated data.
            - str: A step-by-step solution showing the calculations, including trapezoidal rule details.
    """
    
    # 1. Generate variable parameters with validation
    reactant_name = random.choice(GENERAL_REACTANTS)
    F_A0 = round(random.uniform(1.0, 5.0), 2)  # mol/s
    
    # Ensure F_A0 is positive
    if F_A0 <= 0:
        F_A0 = 2.0  # Default fallback
    
    # Variable target conversion instead of fixed 0.8
    target_conversion_options = [0.5, 0.6, 0.65, 0.7, 0.75, 0.8, 0.85, 0.9]
    target_conversion = random.choice(target_conversion_options)
    
    # 2. Generate conversion points dynamically based on target conversion
    if target_conversion <= 0.6:
        num_points = 7
    elif target_conversion <= 0.8:
        num_points = 9
    else:
        num_points = 10  # More points for higher conversions
    
    X_values = np.linspace(0.0, target_conversion, num_points)
    X_values = np.round(X_values, 2)  # Clean up floating point artifacts
    
    # 3. Generate realistic 1/(-r_A) values with proper validation
    # Base pattern that increases with conversion (typical behavior)
    base_pattern = 1.5 + 2.0 * X_values + 3.0 * X_values**2
    
    # Add controlled randomness
    random_factor = random.uniform(0.8, 1.4)
    noise = 1 + 0.15 * np.random.uniform(-1, 1, num_points)
    
    # Ensure noise doesn't make values negative or too small
    noise = np.maximum(noise, 0.5)
    
    inv_rate_values = base_pattern * random_factor * noise
    inv_rate_values = np.round(inv_rate_values, 2)
    
    # 4. Data validation and correction
    # Ensure first value is reasonable (at X=0)
    if inv_rate_values[0] < 1.0:
        inv_rate_values[0] = round(random.uniform(1.0, 2.5), 2)
    
    # Ensure all values are positive
    inv_rate_values = np.maximum(inv_rate_values, 0.5)
    
    # Enforce general increasing trend (physical expectation)
    for i in range(1, len(inv_rate_values)):
        if inv_rate_values[i] < inv_rate_values[i-1]:
            inv_rate_values[i] = inv_rate_values[i-1] + round(random.uniform(0.1, 0.3), 2)
    
    # Final validation - ensure no extremely large values
    inv_rate_values = np.minimum(inv_rate_values, 50.0)
    
    # 5. Calculate CSTR volume with error checking
    try:
        inv_rate_at_final = inv_rate_values[-1]  # 1/(-r_A) at final X
        X_final = X_values[-1]
        
        if inv_rate_at_final <= 0:
            inv_rate_at_final = 5.0  # Fallback value
        
        V_CSTR = F_A0 * inv_rate_at_final * X_final
        
        # Ensure positive volume
        if V_CSTR <= 0:
            V_CSTR = F_A0 * 5.0 * X_final  # Fallback calculation
            
    except Exception:
        # Fallback CSTR calculation
        V_CSTR = F_A0 * 8.0 * target_conversion
    
    # 6. Calculate PFR volume with error checking
    try:
        # Validate data before integration
        if len(X_values) != len(inv_rate_values):
            raise ValueError("Data arrays have mismatched lengths")
        
        if np.any(inv_rate_values <= 0):
            raise ValueError("Negative or zero rate values detected")
        
        # Check for monotonic X values
        if not np.all(np.diff(X_values) >= 0):
            raise ValueError("X values not monotonically increasing")
        
        area_under_curve = np.trapezoid(inv_rate_values, X_values)
        
        if area_under_curve <= 0:
            raise ValueError("Negative area under curve")
        
        V_PFR = F_A0 * area_under_curve
        
        # Sanity check on PFR volume
        if V_PFR <= 0 or V_PFR > 10 * V_CSTR:
            raise ValueError("PFR volume unrealistic")
            
    except Exception:
        # Fallback PFR calculation using simple approximation
        avg_inv_rate = np.mean(inv_rate_values) if len(inv_rate_values) > 0 else 5.0
        V_PFR = F_A0 * avg_inv_rate * target_conversion * 0.7  # Approximate PFR advantage
    
    # 7. Final physical validation
    if V_CSTR <= 0:
        V_CSTR = abs(V_CSTR) + 1.0
    if V_PFR <= 0:
        V_PFR = abs(V_PFR) + 1.0
    
    # Ensure PFR is typically more efficient (but allow exceptions)
    if V_PFR > 1.5 * V_CSTR:
        # Unusual case - might indicate data issues, but proceed with warning
        unusual_case = True
    else:
        unusual_case = False
    
    # 8. Create data table for display
    data_table = "X\t1/(-r_A) [L·s/mol]\n"
    data_table += "-" * 25 + "\n"
    for x, inv_r in zip(X_values, inv_rate_values):
        data_table += f"{x:.2f}\t{inv_r:.2f}\n"
    
    # 9. Generate question 
    question = (
        f"A reaction involving {reactant_name} (A → products) is being studied for reactor design. "
        f"Experimental data has been collected and plotted as a Levenspiel plot (1/(-r_A) vs X). "
        f"The inlet molar flow rate is {F_A0} mol/s. Using the data below, calculate:\n\n"
        f"a) The volume of a CSTR required to achieve {target_conversion*100:.0f}% conversion\n"
        f"b) The volume of a PFR required to achieve {target_conversion*100:.0f}% conversion\n\n"
        f"Levenspiel Plot Data:\n"
        f"{data_table}"
    )
    
    # 10. Generate solution 
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
        f"**Step 1:** For a PFR, the design equation is:\n"
        f"V_PFR = F_A0 × ∫[0 to X] (1/(-r_A)) dX\n\n"
        
        f"**Step 2:** The integral represents the area under the Levenspiel plot curve.\n"
        f"Using trapezoidal rule for numerical integration:\n\n"
        
        f"**Step 3:** Apply trapezoidal rule to the data points:\n"
        f"Area = Σ[(X_i+1 - X_i) × (y_i + y_i+1)/2]\n"
        f"where y_i = [1/(-r_A)]_i\n\n"
        
        f"**Step 4:** Calculate individual trapezoid areas:\n"
    )
    
    # Add detailed trapezoidal calculation with error handling
    total_area = 0
    try:
        for i in range(len(X_values)-1):
            dx = X_values[i+1] - X_values[i]
            avg_height = (inv_rate_values[i] + inv_rate_values[i+1]) / 2
            trap_area = dx * avg_height
            
            # Validate each trapezoid calculation
            if dx <= 0 or avg_height <= 0:
                continue  # Skip invalid intervals
            
            total_area += trap_area
            
            solution += (
                f"Interval [{X_values[i]:.2f} to {X_values[i+1]:.2f}]: "
                f"ΔX = {dx:.2f}, Avg height = ({inv_rate_values[i]:.2f} + {inv_rate_values[i+1]:.2f})/2 = {avg_height:.2f}\n"
                f"Area = {dx:.2f} × {avg_height:.2f} = {trap_area:.3f}\n"
            )
    except Exception:
        # Fallback area calculation
        total_area = area_under_curve if 'area_under_curve' in locals() else np.mean(inv_rate_values) * target_conversion
    
    solution += (
        f"\n**Step 5:** Total area under curve:\n"
        f"Total area = {total_area:.3f}\n\n"
        
        f"**Step 6:** Calculate PFR volume:\n"
        f"V_PFR = F_A0 × Area = {F_A0} mol/s × {total_area:.3f}\n"
        f"V_PFR = {round(V_PFR, 2)} L\n\n"
        
        f"**Final Answers:**\n"
        f"a) CSTR Volume = {round(V_CSTR, 2)} L\n"
        f"b) PFR Volume = {round(V_PFR, 2)} L\n\n"
    )
    
    # Add appropriate comparison note based on results
    if unusual_case:
        comparison_note = (
            f"**Note:** In this case, the PFR volume ({round(V_PFR, 2)} L) is unusually large compared to "
            f"the CSTR volume ({round(V_CSTR, 2)} L). This could indicate very flat kinetics or "
            f"experimental measurement uncertainty."
        )
    elif V_PFR < V_CSTR:
        efficiency = ((V_CSTR - V_PFR) / V_CSTR) * 100
        comparison_note = (
            f"**Note:** The PFR requires less volume than the CSTR ({round(V_PFR, 2)} L vs {round(V_CSTR, 2)} L), "
            f"providing a {efficiency:.1f}% volume reduction due to more efficient use of reaction kinetics."
        )
    else:
        comparison_note = (
            f"**Note:** The reactor volumes are similar ({round(V_CSTR, 2)} L vs {round(V_PFR, 2)} L), "
            f"suggesting relatively flat kinetics over this conversion range."
        )
    
    solution += comparison_note
    
    return question, solution


def main():
    """
    Generate numerous instances of the Levenspiel plot template with different random seeds
    and write the results to a JSONL file.
    """

    import json
    import os

    # Define the output path (Modify this path according to where you are running the code from)
    output_file = "../../../../../testset/reaction_kinetics/conversion_and_reactor_sizing.jsonl"

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
