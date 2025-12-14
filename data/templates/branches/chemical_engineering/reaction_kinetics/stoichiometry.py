import random
from data.templates.branches.chemical_engineering.constants import GENERAL_REACTANTS, LIQUID_PHASE_REACTANTS, GAS_PHASE_REACTANTS, PRODUCTS, REACTIONS


# Template 1 (Easy)
def template_batch_moles_vs_conversion():
    """
    Batch System Moles vs. Conversion

    Scenario:
        This template tests the fundamental relationship between conversion ($X_A$)
        and the number of moles of each species in a batch reactor. For a given
        reaction like $aA + bB -> cC + dD$, the goal is to calculate the final
        number of moles of all species based on the initial moles and a specific
        conversion of the limiting reactant using the core stoichiometric relations:

            $N_A = N_{A0}(1 - X_A)$
            $N_B = N_{B0} - (b/a)N_{A0}X_A$
            $N_C = N_{C0} + (c/a)N_{A0}X_A$

    Returns:
        tuple: A tuple containing:
            - str: A question about calculating final moles in a batch reactor.
            - str: A detailed, step-by-step solution.
    """
    # 1. Parameterize Inputs using Valid Reactions
    # Select a real, pre-balanced reaction to ensure chemical plausibility
    reaction_data = random.choice(REACTIONS)
    equation = reaction_data["equation"]

    # Extract reactants and products
    reactant_items = list(reaction_data["reactants"].items())
    product_items = list(reaction_data["products"].items())

    # Map to A, B, C, D variables
    # We assume the reaction has at least 2 reactants.
    reactant_A_name, a = reactant_items[0]
    reactant_B_name, b = reactant_items[1]
    
    product_C_name, c = product_items[0]
    # Handle case where there might be only 1 product
    has_product_D = len(product_items) > 1
    if has_product_D:
        product_D_name, d = product_items[1]
    else:
        product_D_name, d = "None", 0

    # Generate initial moles
    # Ensure A is the limiting reactant (N_A0/a < N_B0/b)
    N_A0 = round(random.uniform(10.0, 25.0), 2)
    min_N_B0 = (N_A0 / a) * b
    N_B0 = round(min_N_B0 * random.uniform(1.2, 2.0), 2)
    
    # Explicitly set initial product moles to zero for clarity
    N_C0 = 0.0
    N_D0 = 0.0

    # Generate a realistic conversion for the limiting reactant A
    X_A = round(random.uniform(0.40, 0.95), 2)

    # 2. Core Calculations
    N_A = N_A0 * (1 - X_A)
    N_B = N_B0 - (b / a) * N_A0 * X_A
    N_C = N_C0 + (c / a) * N_A0 * X_A
    N_D = N_D0 + (d / a) * N_A0 * X_A if has_product_D else 0.0

    # 3. Generate Question and Solution Strings
    question = (
        f"Consider the following reaction carried out in a batch reactor:\n"
        f"**Reaction:** ${equation}$\n\n"
        f"The reactor is initially charged with ${N_A0}$ moles of {reactant_A_name} and "
        f"${N_B0}$ moles of {reactant_B_name}. Assume the initial amount of products is zero.\n"
        f"{reactant_A_name} is the limiting reactant.\n\n"
        f"If the reaction is allowed to proceed until a conversion of ${X_A}$ ($X_A$) is achieved, "
        f"calculate the final number of moles of all species in the reactor."
    )

    solution = (
        f"**Step 1:** Identify Given Information\n"
        f"- Reaction: ${equation}$\n"
        f"- Initial Moles:\n"
        f"  - $N_{{{reactant_A_name},0}} = {N_A0}$ mol\n"
        f"  - $N_{{{reactant_B_name},0}} = {N_B0}$ mol\n"
        f"  - Products = 0 mol\n"
        f"- Conversion of {reactant_A_name}: $X_A = {X_A}$\n\n"
        
        f"**Step 2:** Write Stoichiometric Relations in Terms of Conversion\n"
        f"The number of moles of each species ($N_j$) at any conversion $X_A$ is given by:\n"
        f"- $N_A = N_{{A,0}}(1 - X_A)$\n"
        f"- $N_B = N_{{B,0}} - ({b}/{a}) N_{{A,0}} X_A$\n"
        f"- $N_C = N_{{C,0}} + ({c}/{a}) N_{{A,0}} X_A$\n"
    )
    if has_product_D:
        solution += f"- $N_D = N_{{D,0}} + ({d}/{a}) N_{{A,0}} X_A$\n"
    solution += "\n"
        
    solution += (
        f"**Step 3:** Calculate Final Moles for Each Species\n"
        f"For {reactant_A_name} (A):\n"
        f"$N_A = {N_A0}(1 - {X_A}) = {N_A0}({1-X_A:.2f}) = {round(N_A, 3)}$ mol\n\n"
        f"For {reactant_B_name} (B):\n"
        f"$N_B = {N_B0} - ({b}/{a}) \\times {N_A0} \\times {X_A} = {N_B0} - {round((b/a)*N_A0*X_A, 3)} = {round(N_B, 3)}$ mol\n\n"
        f"For {product_C_name} (C):\n"
        f"$N_C = 0 + ({c}/{a}) \\times {N_A0} \\times {X_A} = {round((c/a)*N_A0*X_A, 3)} = {round(N_C, 3)}$ mol\n\n"
    )
    
    if has_product_D:
        solution += (
            f"For {product_D_name} (D):\n"
            f"$N_D = 0 + ({d}/{a}) \\times {N_A0} \\times {X_A} = {round((d/a)*N_A0*X_A, 3)} = {round(N_D, 3)}$ mol\n\n"
        )
        
    solution += (
        f"**Final Answer**\n"
        f"After reaching a conversion of ${X_A*100} \%$, the final number of moles in the reactor are:\n"
        f"- {reactant_A_name}: ${round(N_A, 3)}$ mol\n"
        f"- {reactant_B_name}: ${round(N_B, 3)}$ mol\n"
        f"- {product_C_name}: ${round(N_C, 3)}$ mol\n"
    )
    if has_product_D:
        solution += f"- {product_D_name}: ${round(N_D, 3)}$ mol"

    return question, solution


# Template 2 (Intermediate)
def template_flow_system_molar_flow_rates():
    """
    Flow System Molar Flow Rates vs. Conversion

    Scenario:
        This template tests the ability to determine the outlet molar flow rates
        of all species in a continuous flow reactor (e.g., CSTR, PFR) given the
        inlet flow rates and the conversion of a limiting reactant. The calculation
        is an application of stoichiometric principles to a flow system, using:

            F_A = F_A0(1 - X_A)
            F_j = F_j0 + nu_j * F_A0 * X_A = F_A0(Theta_j + nu_j * X_A)

        where nu_j is the stoichiometric coefficient and Theta_j = F_j0/F_A0.

    Returns:
        tuple: A tuple containing:
            - str: A question about calculating outlet molar flow rates.
            - str: A detailed, step-by-step solution.
    """
    # 1. Parameterize Inputs using Valid Reactions
    # Select a real, pre-balanced reaction to ensure chemical plausibility
    reaction_data = random.choice(REACTIONS)
    equation = reaction_data["equation"]

    # Extract reactants and products
    # The dictionaries are in the format {"Name": coefficient}
    reactant_items = list(reaction_data["reactants"].items())
    product_items = list(reaction_data["products"].items())

    # Map to A, B, C, D variables
    # We assume the reaction has at least 2 reactants.
    reactant_A_name, a = reactant_items[0]
    reactant_B_name, b = reactant_items[1]
    
    product_C_name, c = product_items[0]
    # Handle case where there might be only 1 product
    has_product_D = len(product_items) > 1
    if has_product_D:
        product_D_name, d = product_items[1]
    else:
        product_D_name, d = "None", 0

    # Generate inlet molar flow rates (mol/min)
    # Ensure A is the limiting reactant by providing excess B
    F_A0 = round(random.uniform(50.0, 150.0), 2)
    min_F_B0 = (F_A0 / a) * b
    F_B0 = round(min_F_B0 * random.uniform(1.2, 2.5), 2)
    
    # Inlet product flow is usually zero or small
    F_C0 = round(random.choice([0.0, random.uniform(5.0, 20.0)]), 2)
    F_D0 = 0.0

    # Generate a realistic conversion
    X_A = round(random.uniform(0.50, 0.90), 2)

    # 2. Core Calculations
    
    # Calculate Theta values
    Theta_B = F_B0 / F_A0
    Theta_C = F_C0 / F_A0

    # Calculate outlet molar flow rates
    F_A = F_A0 * (1 - X_A)
    F_B = F_B0 - (b / a) * F_A0 * X_A
    F_C = F_C0 + (c / a) * F_A0 * X_A
    F_D = F_D0 + (d / a) * F_A0 * X_A if has_product_D else 0.0
    
    # 3. Generate Question and Solution Strings
    question = (
        f"A reaction is carried out in a steady-state Plug Flow Reactor (PFR):\n"
        f"**Reaction:** {equation}\n\n"
        f"The reactor is fed with {reactant_A_name} at a rate of {F_A0} mol/min, "
        f"{reactant_B_name} at {F_B0} mol/min, and {product_C_name} at {F_C0} mol/min. "
        f"{reactant_A_name} is the limiting reactant.\n\n"
        f"The reactor is designed to achieve a conversion of {X_A} (X_A) for the limiting reactant. "
        f"Calculate the molar flow rates of all species exiting the reactor."
    )

    solution = (
        f"**Step 1:** Identify Given Information\n"
        f"- Reaction: {equation}\n"
        f"- Inlet Molar Flow Rates:\n"
        f"  - F_A0 = {F_A0} mol/min\n"
        f"  - F_B0 = {F_B0} mol/min\n"
        f"  - F_C0 = {F_C0} mol/min\n"
    )
    if has_product_D:
        solution += f"  - F_D0 = {F_D0} mol/min\n"
    
    solution += (
        f"- Conversion of {reactant_A_name}: X_A = {X_A}\n\n"
        
        f"**Step 2:** Stoichiometric Relations for a Flow System\n"
        f"The outlet molar flow rate (F_j) of any species is given by:\n"
        f"F_j = F_j0 + nu_j * F_A0 * X_A\n"
        f"Alternatively, in terms of Theta_j = F_j0 / F_A0:\n"
        f"F_j = F_A0(Theta_j + nu_j * X_A)\n"
        f"Here, the normalized stoichiometric coefficients (nu_j) relative to A are:\n"
        f"nu_A = -1, nu_B = -{b}/{a}, nu_C = +{c}/{a}"
    )
    if has_product_D:
        solution += f", nu_D = +{d}/{a}"
    solution += "\n\n"

    solution += (
        f"**Step 3:** Calculate Outlet Molar Flow Rates\n\n"
        f"**For {reactant_A_name} (A):**\n"
        f"F_A = F_A0(1 - X_A) = {F_A0}(1 - {X_A}) = {round(F_A, 2)} mol/min\n\n"
        
        f"**For {reactant_B_name} (B):**\n"
        f"*Using the Direct Method:*\n"
        f"F_B = F_B0 - ({b}/{a}) * F_A0 * X_A = {F_B0} - ({b}/{a}) * {F_A0} * {X_A} = {round(F_B, 2)} mol/min\n"
        f"*Using the Theta Method:*\n"
        f"Theta_B = F_B0 / F_A0 = {F_B0} / {F_A0} = {round(Theta_B, 3)}\n"
        f"F_B = F_A0(Theta_B - ({b}/{a})X_A) = {F_A0}({round(Theta_B, 3)} - ({b}/{a}) * {X_A}) = {round(F_B, 2)} mol/min\n\n"

        f"**For {product_C_name} (C):**\n"
        f"*Using the Direct Method:*\n"
        f"F_C = F_C0 + ({c}/{a}) * F_A0 * X_A = {F_C0} + ({c}/{a}) * {F_A0} * {X_A} = {round(F_C, 2)} mol/min\n"
        f"*Using the Theta Method:*\n"
        f"Theta_C = F_C0 / F_A0 = {F_C0} / {F_A0} = {round(Theta_C, 3)}\n"
        f"F_C = F_A0(Theta_C + ({c}/{a})X_A) = {F_A0}({round(Theta_C, 3)} + ({c}/{a}) * {X_A}) = {round(F_C, 2)} mol/min\n\n"
    )

    if has_product_D:
        solution += (
            f"**For {product_D_name} (D):**\n"
            f"Since F_D0 = 0, the calculation is straightforward:\n"
            f"F_D = F_D0 + ({d}/{a}) * F_A0 * X_A = 0 + ({d}/{a}) * {F_A0} * {X_A} = {round(F_D, 2)} mol/min\n\n"
        )

    solution += (
        f"**Final Answer**\n"
        f"The molar flow rates exiting the reactor are:\n"
        f"- {reactant_A_name} (F_A): {round(F_A, 2)} mol/min\n"
        f"- {reactant_B_name} (F_B): {round(F_B, 2)} mol/min\n"
        f"- {product_C_name} (F_C): {round(F_C, 2)} mol/min\n"
    )
    if has_product_D:
        solution += f"- {product_D_name} (F_D): {round(F_D, 2)} mol/min"

    return question, solution


# Template 3 (Intermediate)
def template_limiting_reactant():
    """
    Finding and Using the Limiting Reactant

    Scenario:
        This template adds a crucial preliminary step to stoichiometric calculations.
        Given a reaction and the initial moles of multiple reactants, the user must
        first identify the limiting reactant by comparing the ratio of initial moles
        to the stoichiometric coefficient for each species:

            Compare: $N_{A0}/a$ vs. $N_{B0}/b$

        The species with the smaller ratio is limiting. Then, using a given
        conversion with respect to that limiting reactant, the user calculates the
        final moles of all species.

    Returns:
        tuple: A tuple containing:
            - str: A question requiring identification of the limiting reactant.
            - str: A detailed solution showing the identification and calculation steps.
    """
    # 1. Parameterize Inputs using Valid Reactions
    # Select a real, pre-balanced reaction to ensure chemical plausibility
    reaction_data = random.choice(REACTIONS)
    
    # Extract reactants and products
    # The dictionaries are in the format {"Name": coefficient}
    reactant_items = list(reaction_data["reactants"].items())
    product_items = list(reaction_data["products"].items())

    # Map to A, B, C, D variables for the template logic
    # We assume the reaction has at least 2 reactants and 1-2 products based on the constants file
    reactant_A_name, a = reactant_items[0]
    reactant_B_name, b = reactant_items[1]
    
    product_C_name, c = product_items[0]
    # Handle case where there might be only 1 product or 2
    if len(product_items) > 1:
        product_D_name, d = product_items[1]
    else:
        # Create a dummy placeholder if only 1 product exists to prevent errors
        product_D_name, d = "Heat/Other", 0

    # Randomly decide which reactant will be limiting to vary the problem type
    is_A_limiting = random.choice([True, False])

    # Generate initial moles based on which reactant is limiting
    if is_A_limiting:
        N_A0 = round(random.uniform(20.0, 50.0), 2)
        # Make B the excess reactant (ratio N_B0/b > N_A0/a)
        min_B = (N_A0 / a) * b
        N_B0 = round(min_B * random.uniform(1.2, 2.0), 2)
    else: # B is limiting
        N_B0 = round(random.uniform(20.0, 50.0), 2)
        # Make A the excess reactant (ratio N_A0/a > N_B0/b)
        min_A = (N_B0 / b) * a
        N_A0 = round(min_A * random.uniform(1.2, 2.0), 2)
    
    # Explicitly define initial product moles as zero
    N_C0 = 0.0
    N_D0 = 0.0

    # Generate a realistic conversion for the limiting reactant
    X = round(random.uniform(0.60, 0.95), 2)

    # 2. Core Calculations & Logic
    
    # Determine the limiting reactant
    ratio_A = N_A0 / a
    ratio_B = N_B0 / b

    if ratio_A <= ratio_B:
        # Case 1: A is the limiting reactant
        limiting_reactant_name = reactant_A_name
        limiting_reactant_sym = 'A'
        
        # Calculate final moles based on X_A
        N_A = N_A0 * (1 - X)
        N_B = N_B0 - (b / a) * N_A0 * X
        N_C = N_C0 + (c / a) * N_A0 * X
        N_D = N_D0 + (d / a) * N_A0 * X
    else:
        # Case 2: B is the limiting reactant
        limiting_reactant_name = reactant_B_name
        limiting_reactant_sym = 'B'

        # Calculate final moles based on X_B
        N_A = N_A0 - (a / b) * N_B0 * X
        N_B = N_B0 * (1 - X)
        N_C = N_C0 + (c / b) * N_B0 * X
        N_D = N_D0 + (d / b) * N_B0 * X

    # 3. Generate Question and Solution Strings
    equation = reaction_data["equation"]

    question = (
        f"The following reaction occurs in a batch reactor:\n"
        f"**Reaction:** {equation}\n\n"
        f"The reactor is initially charged with {N_A0} moles of {reactant_A_name} and "
        f"{N_B0} moles of {reactant_B_name}. Assume the initial amount of products is zero.\n"
        f"The reaction is allowed to proceed until a {X*100}% conversion of the limiting reactant is achieved.\n\n"
        f"First, identify the limiting reactant. Then, calculate the final number of moles for all species."
    )

    solution_step1_identification = (
        f"**Step 1:** Identify the Limiting Reactant\n"
        f"To find the limiting reactant, we compare the ratio of the initial moles to the "
        f"stoichiometric coefficient for each reactant.\n\n"
        f"- For **{reactant_A_name}**: N_A0/a = {N_A0} / {a} = **{round(ratio_A, 2)}**\n"
        f"- For **{reactant_B_name}**: N_B0/b = {N_B0} / {b} = **{round(ratio_B, 2)}**\n\n"
        f"Since {round(min(ratio_A, ratio_B), 2)} < {round(max(ratio_A, ratio_B), 2)}, "
        f"**{limiting_reactant_name} is the limiting reactant**.\n"
    )

    if limiting_reactant_sym == 'A':
        solution_step2_calculation = (
            f"\n**Step 2:** Calculate Final Moles\n"
            f"Using conversion of {reactant_A_name} (X = {X}) as the basis:\n\n"
            f"Moles {reactant_A_name} = {N_A0}(1 - {X}) = **{round(N_A, 2)}** mol\n"
            f"Moles {reactant_B_name} = {N_B0} - ({b}/{a})({N_A0})({X}) = **{round(N_B, 2)}** mol\n"
            f"Moles {product_C_name} = 0 + ({c}/{a})({N_A0})({X}) = **{round(N_C, 2)}** mol\n"
        )
        if d > 0:
            solution_step2_calculation += f"Moles {product_D_name} = 0 + ({d}/{a})({N_A0})({X}) = **{round(N_D, 2)}** mol\n"
    else: 
        solution_step2_calculation = (
            f"\n**Step 2:** Calculate Final Moles\n"
            f"Using conversion of {reactant_B_name} (X = {X}) as the basis:\n\n"
            f"Moles {reactant_A_name} = {N_A0} - ({a}/{b})({N_B0})({X}) = **{round(N_A, 2)}** mol\n"
            f"Moles {reactant_B_name} = {N_B0}(1 - {X}) = **{round(N_B, 2)}** mol\n"
            f"Moles {product_C_name} = 0 + ({c}/{b})({N_B0})({X}) = **{round(N_C, 2)}** mol\n"
        )
        if d > 0:
            solution_step2_calculation += f"Moles {product_D_name} = 0 + ({d}/{b})({N_B0})({X}) = **{round(N_D, 2)}** mol\n"

    solution_final_answer = (
        f"\n**Final Answer**\n"
        f"The limiting reactant is **{limiting_reactant_name}**. The final number of moles are:\n"
        f"- **{reactant_A_name}:** {round(N_A, 2)} mol\n"
        f"- **{reactant_B_name}:** {round(N_B, 2)} mol\n"
        f"- **{product_C_name}:** {round(N_C, 2)} mol\n"
    )
    if d > 0:
        solution_final_answer += f"- **{product_D_name}:** {round(N_D, 2)} mol"

    solution = f"{solution_step1_identification}{solution_step2_calculation}{solution_final_answer}"

    return question, solution


# Template 4 (Advanced)
def template_gas_phase_concentration():
    """
    Gas-Phase Concentration with Volumetric Change

    Scenario:
        For gas-phase reactions with a change in the total number of moles, the
        volumetric flow rate varies with conversion. This template tests the ability
        to calculate outlet concentrations in an isothermal, isobaric flow system
        using the specialized gas-phase concentration equations from Fogler, which
        account for this volumetric change via the parameter epsilon ($\\epsilon$).

        $C_A = C_{A0} \\frac{1 - X_A}{1 + \\epsilon X_A}$
        $C_B = C_{A0} \\frac{\\Theta_B - (b/a)X_A}{1 + \\epsilon X_A}$

    Returns:
        tuple: A tuple containing:
            - str: A question about calculating gas-phase outlet concentrations.
            - str: A detailed solution showing the application of the correct formulas.
    """
    # 1. Randomized Parameters
    
    # Select unique names
    reactant_A_name, reactant_B_name = random.sample(GAS_PHASE_REACTANTS, 2)
    product_C_name = random.choice(PRODUCTS)
    
    # Stoichiometric coefficients
    a = random.choice([1, 2])
    b = random.randint(1, 3)
    c = random.randint(1, 2)

    # Initial concentration (gases have lower concentrations, units are mol/dm^3)
    C_A0 = round(random.uniform(0.05, 0.25), 3)
    
    # Conversion
    X_A = round(random.uniform(0.50, 0.90), 2)
    
    # Theta_B must be in excess of the stoichiometric requirement
    Theta_B = round((b / a) * random.uniform(1.2, 2.5), 2)
    
    # Epsilon (ε) can be positive (expansion) or negative (contraction)
    # Ensure it's not too close to zero to make the problem meaningful
    epsilon = 0
    while abs(epsilon) < 0.1:
        epsilon = round(random.uniform(-0.5, 0.5), 2)

    # 2. Core Calculations
    denominator = 1 + epsilon * X_A
    C_A = C_A0 * (1 - X_A) / denominator
    C_B = C_A0 * (Theta_B - (b / a) * X_A) / denominator

    # 3. Generate Question and Solution Strings
    def format_species(coeff, name):
        return f"{coeff if coeff > 1 else ''}{name}"

    reaction_string = (
        f"{format_species(a, reactant_A_name)} + {format_species(b, reactant_B_name)} -> "
        f"{format_species(c, product_C_name)}"
    )

    question = (
        f"The following elementary gas-phase reaction occurs in a steady-state PFR:\n"
        f"**Reaction:** ${reaction_string}$\n\n"
        f"The reaction is carried out **isothermally** and **isobarically**. The feed enters the reactor with an initial concentration of {reactant_A_name} of $C_{{A0}} = {C_A0}$ mol/dm³.\n\n"
        f"The following parameters are known:\n"
        f"- Molar feed ratio: $\\Theta_B = F_{{B0}}/F_{{A0}} = {Theta_B}$\n"
        f"- Volumetric change parameter: $\\epsilon = {epsilon}$\n\n"
        f"If the reaction achieves a final conversion of $X_A = {X_A}$, what are the outlet concentrations of {reactant_A_name} ($C_A$) and {reactant_B_name} ($C_B$)?"
    )

    solution = (
        f"**Step 1:** Identify Given Information\n"
        f"- Initial Concentration: $C_{{A0}} = {C_A0}$ mol/dm³\n"
        f"- Conversion: $X_A = {X_A}$\n"
        f"- Molar Feed Ratio: $\\Theta_B = {Theta_B}$\n"
        f"- Volumetric Change Parameter: $\\epsilon = {epsilon}$\n"
        f"- Stoichiometric Ratio: $b/a = {b}/{a} = {round(b/a, 2)}$\n\n"
        
        f"**Step 2:** State the Governing Equations\n"
        f"For an isothermal, isobaric gas-phase reaction with a change in the number of moles, the concentration of each species is a function of conversion ($X_A$) and the volumetric change parameter ($\\epsilon$). The term $(1 + \\epsilon X_A)$ in the denominator corrects for the change in volumetric flow rate.\n\n"
        f"For reactant A: \n$C_A = C_{{A0}} \\frac{{1 - X_A}}{{1 + \\epsilon X_A}}$\n\n"
        f"For reactant B: \n$C_B = C_{{A0}} \\frac{{\\Theta_B - (b/a)X_A}}{{1 + \\epsilon X_A}}$\n\n"
        
        f"**Step 3:** Calculate the Denominator Term\n"
        f"First, let's calculate the volumetric correction term, which is common to all species:\n"
        f"$1 + \\epsilon X_A = 1 + ({epsilon})({X_A}) = {round(denominator, 4)}$\n\n"
        
        f"**Step 4:** Calculate Outlet Concentrations\n\n"
        f"**For {reactant_A_name} ($C_A$):**\n"
        f"$C_A = {C_A0} \\frac{{1 - {X_A}}}{{{round(denominator, 4)}}} = {C_A0} \\frac{{{round(1 - X_A, 2)}}}{{{round(denominator, 4)}}} = \\mathbf{{{round(C_A, 4)}}}$ mol/dm³\n\n"
        
        f"**For {reactant_B_name} ($C_B$):**\n"
        f"First, calculate the numerator term for B: \n"
        f"$\\Theta_B - (b/a)X_A = {Theta_B} - ({round(b/a, 2)})({X_A}) = {round(Theta_B - (b/a)*X_A, 4)}$\n"
        f"Now, calculate the concentration:\n"
        f"$C_B = {C_A0} \\frac{{{round(Theta_B - (b/a)*X_A, 4)}}}{{{round(denominator, 4)}}} = \\mathbf{{{round(C_B, 4)}}}$ mol/dm³\n\n"
        
        f"**Final Answer**\n"
        f"The outlet concentrations are:\n"
        f"- **$C_A$:** ${round(C_A, 4)}$ mol/dm³\n"
        f"- **$C_B$:** ${round(C_B, 4)}$ mol/dm³"
    )

    return question, solution


def main():
    """
    Generate numerous instances of each stoichiometry template with different random seeds
    and write the results to a JSONL file.
    """
    import json
    import os

    # Define the output path (Modify this path according to where you are running the code from)
    output_file = "testset/chemical_engineering/reaction_kinetics/stoichiometry.jsonl"

    # Create the directory if it doesn't exist
    os.makedirs(os.path.dirname(output_file), exist_ok=True)

    # List of template functions with their ID and level
    templates = [ 
        (template_batch_moles_vs_conversion, "batch_moles_vs_conversion", "Easy"),
        (template_flow_system_molar_flow_rates, "flow_rates_vs_conversion", "Intermediate"),
        (template_limiting_reactant, "finding_limiting_reactant", "Intermediate"),
        (template_gas_phase_concentration, "gas_phase_concentration", "Advanced"),
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
                "area": "stoichiometry",
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
