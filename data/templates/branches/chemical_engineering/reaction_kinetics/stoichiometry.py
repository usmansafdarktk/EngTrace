import random
import math
from decimal import Decimal, ROUND_HALF_UP
from data.templates.branches.chemical_engineering.constants import GENERAL_REACTANTS, LIQUID_PHASE_REACTANTS, GAS_PHASE_REACTANTS, PRODUCTS, REACTIONS


def _hu(x, places):
    """Round half-up to `places` dp, resolving the tie in DECIMAL.

    `round()` resolves a half-way tie on the binary value and disagrees with a
    reader doing decimal arithmetic (spec P2 as amended, DECISIONS D-012).
    """
    q = Decimal(1).scaleb(-places)
    d = x if isinstance(x, Decimal) else Decimal(repr(x))
    v = d.quantize(q, rounding=ROUND_HALF_UP)
    return int(v) if places == 0 else float(v)


def _as_printed(x, spec):
    """The value a reader recovers from `x` when it is printed with `spec`.

    P2 asks that the stored value and the printed value be the SAME value.
    Rounding alone does not achieve that: a float one ulp away from its own
    printed form puts the template and the reader on opposite sides of a
    display tie (D-016 part 2).
    """
    return float(format(x, spec))


def _is_display_tie(x, places, rel_band=1e-12):
    """Is `x` at, or within a hair of, a half-way tie at `places` dp?

    A tie is the one case where no rounding convention is defensible - a
    decimal reader applying half-up and a binary reader applying `round()`
    disagree, and the printed line closes for only one of them. Such instances
    are resampled rather than resolved (D-016).

    A narrow BAND is quarantined rather than a point, because two independent
    evaluations of the same exact quantity differ by a few ulps and an exact
    rational tie lands on opposite sides of them.
    """
    scaled = abs(x) * 10.0 ** places
    band = max(scaled * rel_band, 1e-9)
    return abs((scaled - math.floor(scaled)) - 0.5) <= band


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

    Trace integrity (Layer 0, 2026-09-23):
        The consumed products (b/a)*N_A0*X_A are the products the solution
        PRINTS: each is bound through its 3-dp display before N_B is formed
        from it, so `N_B0 - <printed product> = <printed N_B>` closes for
        every reader (D-016 part 2); before, N_B was computed from the
        unrounded product and the line missed by 0.001 on 3.4% of draws. The
        products themselves are exact at 4 dp (6 dp for the ammonia
        reaction), so their 3-dp display, like the 3-dp display of every
        answer, still sits on a half-way tie for a fixed fraction of draws;
        removing that requires lengthening the ANSWER display, which is a
        P6 change needing sign-off (D-044) and is recorded, not done, here.

    Screen pass 1 (2026-09-23):
        Two judges flagged `55.00000000000001 %` in the answer sentence
        (X_A*100 interpolated raw; 9% of draws) and N_A printed through
        round() while the other species are bound through their 3-dp
        display. The percentage is now formatted at 1 dp (X_A itself is
        unchanged) and N_A is bound through its 3-dp display like N_B, N_C
        and N_D, so the stored value is the printed one. N_A0*(1 - X_A) is
        exact at 4 dp, so its 3-dp display sits on a half-way tie on 10% of
        draws, exactly as the products' displays do (18.8% of instances
        carry at least one such tie; the census cannot see these lines,
        whose products are implicit or `\\times`). Removing them means
        printing the answers at their exact length: the D-044 decision
        recorded above, still pending. Done 2026-09-24 on the owner's
        decision (the D-044 sign-off, the round-3 policy applied to the
        sibling flow template): the answer display is LENGTHENED to the
        exact precision of the instance, `_p` = 4 dp for a = 1 and 6 dp for
        the ammonia reaction (a = 4, ratios 5/4 and 6/4), and N_A, the three
        products, N_B, N_C and N_D are computed in exact decimal arithmetic
        bound half-up through that display in Step 3 and the Answer block,
        so nothing rounds and no line can tie. The question is unchanged on
        every seed; every gold answer's text lengthens, and its value moves
        on the 75.8% of instances whose 3-dp display had rounded a 4- or
        6-dp exact quantity (500 seeds).

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

    # 2. Core Calculations, in exact decimal arithmetic and bound half-up
    # through the display (D-016 part 2), as the sibling flow template does.
    # The display precision is the exact precision of the instance: N_A0 * X_A
    # and N_A0 * (1 - X_A) are exact at 4 dp, and dividing by a = 4 (ammonia
    # oxidation, ratios 5/4 and 6/4) adds 2 dp; a = 1 for the other three
    # reactions. N_A, the three products and N_B, N_C, N_D are all exact at
    # `_p`, so nothing rounds and no printed line can tie (D-044, the owner's
    # sign-off). The products are PRINTED and then consumed, so N_B is formed
    # from the printed product, which is the exact one.
    _p = {1: 4, 2: 5, 4: 6, 5: 5}[a]
    _dA0, _dB0, _dX = Decimal(repr(N_A0)), Decimal(repr(N_B0)), Decimal(repr(X_A))
    N_A = _hu(_dA0 * (1 - _dX), _p)
    prod_B = _hu(_dA0 * _dX * b / a, _p)
    prod_C = _hu(_dA0 * _dX * c / a, _p)
    prod_D = _hu(_dA0 * _dX * d / a, _p) if has_product_D else 0.0
    N_B = _hu(_dB0 - _dA0 * _dX * b / a, _p)
    N_C = _hu(Decimal(repr(N_C0)) + _dA0 * _dX * c / a, _p)
    N_D = _hu(Decimal(repr(N_D0)) + _dA0 * _dX * d / a, _p) if has_product_D else 0.0

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
        f"$N_A = {N_A0}(1 - {X_A}) = {N_A0}({1-X_A:.2f}) = {N_A:.{_p}f}$ mol\n\n"
        f"For {reactant_B_name} (B):\n"
        f"$N_B = {N_B0} - ({b}/{a}) \\times {N_A0} \\times {X_A} = {N_B0} - {prod_B:.{_p}f} = {N_B:.{_p}f}$ mol\n\n"
        f"For {product_C_name} (C):\n"
        f"$N_C = 0 + ({c}/{a}) \\times {N_A0} \\times {X_A} = {prod_C:.{_p}f} = {N_C:.{_p}f}$ mol\n\n"
    )
    
    if has_product_D:
        solution += (
            f"For {product_D_name} (D):\n"
            f"$N_D = 0 + ({d}/{a}) \\times {N_A0} \\times {X_A} = {prod_D:.{_p}f} = {N_D:.{_p}f}$ mol\n\n"
        )
        
    solution += (
        f"**Answer:**\n"
        f"After reaching a conversion of ${X_A*100:.1f} \\%$, the final number of moles in the reactor are:\n"
        f"- {reactant_A_name}: ${N_A:.{_p}f}$ mol\n"
        f"- {reactant_B_name}: ${N_B:.{_p}f}$ mol\n"
        f"- {product_C_name}: ${N_C:.{_p}f}$ mol\n"
    )
    if has_product_D:
        solution += f"- {product_D_name}: ${N_D:.{_p}f}$ mol"

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

    Trace integrity (Layer 0, 2026-09-23):
        Each outlet flow is a 2-dp inlet flow plus a 2-dp flow times a
        2-dp conversion times a small integer ratio nu_j = b/a: exact at
        4 dp when a = 1 (three of the four reactions) and at 6 dp for the
        ammonia reaction (a = 4, nu_B = 5/4, nu_D = 6/4). Printed at 2 dp,
        F_B, F_C and F_D sat on a half-way tie on 7% of draws
        (748.18 - 5 * 109.26 * 0.75 = 338.455). These are the answers, so
        per the round-3 policy their display is LENGTHENED to the exact
        precision of the instance (`_p` = 4 or 6 dp by the reaction's a)
        and every result is bound half-up in exact decimal arithmetic
        through that display (D-016 part 2, as the sibling batch template
        does), so nothing rounds and no line can tie. Theta_B and Theta_C
        remain 3-dp quotients shown as the cross-check; the direct-method
        line is the one that closes exactly. (The cross-check lines are
        reworked in the Layer 2 fix below.)

    Layer 2 fix (2026-09-25):
        One expert found the Step 3 "Theta Method" lines substituting Theta
        rounded to 3 dp but printing the exact direct-method result, so they
        did not follow from their printed operands: instance 2,
        120.26(1.665 - (5/4)(0.67)) = 99.515, printed 99.572250; instance 3,
        100.93(0.151 + 0.52) = 67.724, printed 67.7736. Confirmed on seeds
        2102 and 2103. Theta_B and Theta_C are now printed at 4 dp and bound
        through that display, and each Theta-method line is computed from
        the printed Theta in exact decimal arithmetic at the 6 dp the product
        has (a 4-dp Theta minus a 2-dp conversion times a small ratio, times
        a 2-dp flow), so the line closes and nothing rounds (D-037); a
        following sentence says the small difference from the direct-method
        value is the rounding of Theta. The direct method remains the answer
        and no gold answer changes. A quotient F_j0/F_A0 whose 4-dp display
        sits on a half-way tie (possible when 100*F_A0 has only the factors
        2 and 5) is redrawn (D-016); no question text changes on a kept draw.

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

    for _attempt in range(200):
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

        # Theta_B and Theta_C are printed at 4 dp and then consumed by the
        # cross-check lines, so they are bound through that display (D-016
        # part 2); a quotient whose display sits on a half-way tie is redrawn.
        if _is_display_tie(F_B0 / F_A0, 4) or _is_display_tie(F_C0 / F_A0, 4):
            continue
        break
    else:
        raise RuntimeError(
            "flow_system_molar_flow_rates: no display-stable sample in 200 draws")

    # 2. Core Calculations

    # Calculate Theta values, bound through their 4-dp display
    Theta_B = _as_printed(F_B0 / F_A0, '.4f')
    Theta_C = _as_printed(F_C0 / F_A0, '.4f')

    # Calculate outlet molar flow rates, in exact decimal arithmetic and
    # bound through their display (D-016 part 2). The display precision is
    # the exact precision of the instance: F_A0 * X_A is exact at 4 dp, and
    # dividing by a = 4 (ammonia oxidation) adds 2 dp; a = 1 otherwise.
    _p = {1: 4, 2: 5, 4: 6, 5: 5}[a]
    _dA0, _dX = Decimal(repr(F_A0)), Decimal(repr(X_A))
    F_A = _hu(_dA0 * (1 - _dX), _p)
    F_B = _hu(Decimal(repr(F_B0)) - _dA0 * _dX * b / a, _p)
    F_C = _hu(Decimal(repr(F_C0)) + _dA0 * _dX * c / a, _p)
    F_D = _hu(Decimal(repr(F_D0)) + _dA0 * _dX * d / a, _p) if has_product_D else 0.0

    # The Theta-method cross-check lines are computed from the PRINTED
    # 4-dp Theta (Layer 2): a 4-dp Theta minus a 2-dp conversion times a
    # ratio with a in {1, 2, 4, 5}, times a 2-dp flow, is exact at 6 dp, so
    # the line is printed at that length and nothing rounds (D-037).
    _dTB, _dTC = Decimal(repr(Theta_B)), Decimal(repr(Theta_C))
    F_B_theta = _hu(_dA0 * (_dTB - _dX * b / a), 6)
    F_C_theta = _hu(_dA0 * (_dTC + _dX * c / a), 6)

    def _theta_note(sym, theta_val, direct_val):
        if theta_val == direct_val:
            return (f"(Theta_{sym} is exact at 4 dp here, so the two methods "
                    f"agree exactly.)")
        return (f"(The small difference from the direct-method value is the "
                f"rounding of Theta_{sym} to 4 dp; the direct-method value is exact.)")
    theta_note_B = _theta_note('B', F_B_theta, F_B)
    theta_note_C = _theta_note('C', F_C_theta, F_C)
    
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
        f"F_A = F_A0(1 - X_A) = {F_A0}(1 - {X_A}) = {F_A:.{_p}f} mol/min\n\n"
        
        f"**For {reactant_B_name} (B):**\n"
        f"*Using the Direct Method:*\n"
        f"F_B = F_B0 - ({b}/{a}) * F_A0 * X_A = {F_B0} - ({b}/{a}) * {F_A0} * {X_A} = {F_B:.{_p}f} mol/min\n"
        f"*Using the Theta Method (cross-check):*\n"
        f"Theta_B = F_B0 / F_A0 = {F_B0} / {F_A0} = {Theta_B:.4f}\n"
        f"F_B = F_A0(Theta_B - ({b}/{a})X_A) = {F_A0}({Theta_B:.4f} - ({b}/{a}) * {X_A}) = {F_B_theta:.6f} mol/min\n"
        f"{theta_note_B}\n\n"

        f"**For {product_C_name} (C):**\n"
        f"*Using the Direct Method:*\n"
        f"F_C = F_C0 + ({c}/{a}) * F_A0 * X_A = {F_C0} + ({c}/{a}) * {F_A0} * {X_A} = {F_C:.{_p}f} mol/min\n"
        f"*Using the Theta Method (cross-check):*\n"
        f"Theta_C = F_C0 / F_A0 = {F_C0} / {F_A0} = {Theta_C:.4f}\n"
        f"F_C = F_A0(Theta_C + ({c}/{a})X_A) = {F_A0}({Theta_C:.4f} + ({c}/{a}) * {X_A}) = {F_C_theta:.6f} mol/min\n"
        f"{theta_note_C}\n\n"
    )

    if has_product_D:
        solution += (
            f"**For {product_D_name} (D):**\n"
            f"Since F_D0 = 0, the calculation is straightforward:\n"
            f"F_D = F_D0 + ({d}/{a}) * F_A0 * X_A = 0 + ({d}/{a}) * {F_A0} * {X_A} = {F_D:.{_p}f} mol/min\n\n"
        )

    solution += (
        f"**Answer:**\n"
        f"The molar flow rates exiting the reactor are:\n"
        f"- {reactant_A_name} (F_A): {F_A:.{_p}f} mol/min\n"
        f"- {reactant_B_name} (F_B): {F_B:.{_p}f} mol/min\n"
        f"- {product_C_name} (F_C): {F_C:.{_p}f} mol/min\n"
    )
    if has_product_D:
        solution += f"- {product_D_name} (F_D): {F_D:.{_p}f} mol/min"

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
        f"\n**Answer:**\n"
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

    Screen pass 1 (2026-09-23):
        All three judges flagged placeholder product names (`Entity 3`,
        `Species VII`) under random coefficients, and one that epsilon was
        drawn independently of the stoichiometry. Both confirmed. The
        reaction is now drawn from 24 real gas-phase reactions with two
        reactants and a change in total moles (19 contractions, 5
        expansions), named and balanced, and epsilon = y_A0*delta with
        delta = (c + d - a - b)/a and y_A0 = 1/(1 + Theta_B) for a feed of A
        and B only, which the question now states; Step 1 shows the check.
        y_A0 is bound through its 4-dp display, y_A0*delta is exact at 4 dp
        (5 dp for a half-integer delta) and printed at that length (D-037),
        epsilon is stated at 2 dp and bound through it, and a draw whose
        y_A0, epsilon, C_A or C_B display sits on a half-way tie is redrawn
        (D-016). Every question and gold answer changes (the reaction list
        and the random stream differ); |epsilon| now spans 0.20-0.91 (500
        seeds) instead of the drawn 0.10-0.50, with 17% of instances
        expanding; the tie screens redraw 2.0% of draws (1.6% on the 2-dp
        epsilon display, scattered over Theta_B).

    Returns:
        tuple: A tuple containing:
            - str: A question about calculating gas-phase outlet concentrations.
            - str: A detailed solution showing the application of the correct formulas.
    """
    # 1. Randomized Parameters

    # Real gas-phase reactions with two reactants and a change in the total
    # number of moles (screen pass 1: the product used to be drawn from a
    # placeholder list, `Acetaldehyde + Ethane -> 2Entity 3`, under random
    # coefficients). Fields: A, B, a, b, products [(name, coeff)], equation.
    reactions = [
        ("Ethylene", "Hydrogen", 1, 1, [("Ethane", 1)], "C2H4(g) + H2(g) → C2H6(g)"),
        ("Propylene", "Hydrogen", 1, 1, [("Propane", 1)], "C3H6(g) + H2(g) → C3H8(g)"),
        ("Acetylene", "Hydrogen", 1, 1, [("Ethylene", 1)], "C2H2(g) + H2(g) → C2H4(g)"),
        ("Ethylene", "Hydrogen Chloride", 1, 1, [("Ethyl Chloride", 1)], "C2H4(g) + HCl(g) → C2H5Cl(g)"),
        ("Acetylene", "Hydrogen Chloride", 1, 1, [("Vinyl Chloride", 1)], "C2H2(g) + HCl(g) → C2H3Cl(g)"),
        ("Ethylene", "Chlorine", 1, 1, [("1,2-Dichloroethane", 1)], "C2H4(g) + Cl2(g) → C2H4Cl2(g)"),
        ("Carbon Monoxide", "Chlorine", 1, 1, [("Phosgene", 1)], "CO(g) + Cl2(g) → COCl2(g)"),
        ("Ethylene", "Water", 1, 1, [("Ethanol", 1)], "C2H4(g) + H2O(g) → C2H5OH(g)"),
        ("Formaldehyde", "Hydrogen", 1, 1, [("Methanol", 1)], "CH2O(g) + H2(g) → CH3OH(g)"),
        ("Acetaldehyde", "Hydrogen", 1, 1, [("Ethanol", 1)], "CH3CHO(g) + H2(g) → C2H5OH(g)"),
        ("Butadiene", "Ethylene", 1, 1, [("Cyclohexene", 1)], "C4H6(g) + C2H4(g) → C6H10(g)"),
        ("Carbon Monoxide", "Hydrogen", 1, 2, [("Methanol", 1)], "CO(g) + 2H2(g) → CH3OH(g)"),
        ("Carbon Monoxide", "Hydrogen", 1, 3, [("Methane", 1), ("Water", 1)], "CO(g) + 3H2(g) → CH4(g) + H2O(g)"),
        ("Carbon Dioxide", "Hydrogen", 1, 3, [("Methanol", 1), ("Water", 1)], "CO2(g) + 3H2(g) → CH3OH(g) + H2O(g)"),
        ("Carbon Dioxide", "Hydrogen", 1, 4, [("Methane", 1), ("Water", 2)], "CO2(g) + 4H2(g) → CH4(g) + 2H2O(g)"),
        ("Nitrogen", "Hydrogen", 1, 3, [("Ammonia", 2)], "N2(g) + 3H2(g) → 2NH3(g)"),
        ("Sulfur Dioxide", "Oxygen", 2, 1, [("Sulfur Trioxide", 2)], "2SO2(g) + O2(g) → 2SO3(g)"),
        ("Nitric Oxide", "Oxygen", 2, 1, [("Nitrogen Dioxide", 2)], "2NO(g) + O2(g) → 2NO2(g)"),
        ("Ethylene", "Oxygen", 2, 1, [("Ethylene Oxide", 2)], "2C2H4(g) + O2(g) → 2C2H4O(g)"),
        ("Methane", "Water", 1, 1, [("Carbon Monoxide", 1), ("Hydrogen", 3)], "CH4(g) + H2O(g) → CO(g) + 3H2(g)"),
        ("Methane", "Carbon Dioxide", 1, 1, [("Carbon Monoxide", 2), ("Hydrogen", 2)], "CH4(g) + CO2(g) → 2CO(g) + 2H2(g)"),
        ("Methane", "Oxygen", 2, 1, [("Carbon Monoxide", 2), ("Hydrogen", 4)], "2CH4(g) + O2(g) → 2CO(g) + 4H2(g)"),
        ("Ethane", "Oxygen", 2, 1, [("Ethylene", 2), ("Water", 2)], "2C2H6(g) + O2(g) → 2C2H4(g) + 2H2O(g)"),
        ("Propane", "Oxygen", 2, 1, [("Propylene", 2), ("Water", 2)], "2C3H8(g) + O2(g) → 2C3H6(g) + 2H2O(g)"),
    ]

    for _attempt in range(200):
        reactant_A_name, reactant_B_name, a, b, products, equation = random.choice(reactions)
        # delta: change in total moles per mole of A reacted, (c + d - a - b)/a.
        # Exact, since a is 1 or 2.
        n_products = sum(coeff for _name, coeff in products)
        delta = (n_products - a - b) / a

        # Initial concentration (gases have lower concentrations, units are mol/dm^3)
        C_A0 = round(random.uniform(0.05, 0.25), 3)

        # Conversion
        X_A = round(random.uniform(0.50, 0.90), 2)

        # Theta_B must be in excess of the stoichiometric requirement
        Theta_B = round((b / a) * random.uniform(1.2, 2.5), 2)

        # Epsilon follows from the stoichiometry and the feed (Fogler): the
        # feed is A and B only, so y_A0 = 1/(1 + Theta_B) and epsilon =
        # y_A0*delta. It used to be drawn independently of both. y_A0 is
        # printed at 4 dp and consumed, so it is bound through that display;
        # y_A0*delta is exact at 4 dp (5 dp when delta is a half-integer) and
        # printed at that length (D-037); epsilon is STATED at 2 dp, bound
        # through it, and consumed downstream as stated. A draw whose y_A0 or
        # epsilon display sits on a half-way tie is redrawn (D-016).
        if _is_display_tie(1 / (1 + Theta_B), 4):
            continue
        y_A0 = _as_printed(1 / (1 + Theta_B), '.4f')
        eps_dp = 4 if delta == int(delta) else 5
        eps_full = _as_printed(y_A0 * delta, f'.{eps_dp}f')
        if _is_display_tie(eps_full, 2):
            continue
        epsilon = _as_printed(eps_full, '.2f')

        # 2. Core Calculations. The denominator is exact at 4 dp and bound
        # through that display; the answers are printed at 4 dp, and a
        # quotient sitting on a half-way tie there is redrawn (D-016).
        denominator = _as_printed(1 + epsilon * X_A, '.4f')
        C_A = C_A0 * (1 - X_A) / denominator
        C_B = C_A0 * (Theta_B - (b / a) * X_A) / denominator
        if _is_display_tie(C_A, 4) or _is_display_tie(C_B, 4):
            continue
        break
    else:
        raise RuntimeError(
            "gas_phase_concentration: no display-stable sample in 200 draws")

    # 3. Generate Question and Solution Strings
    def format_species(coeff, name):
        return f"{coeff} {name}" if coeff > 1 else name

    reaction_string = (
        f"{format_species(a, reactant_A_name)} + {format_species(b, reactant_B_name)} → "
        + " + ".join(format_species(coeff, name) for name, coeff in products)
    )

    question = (
        f"The following gas-phase reaction occurs in a steady-state PFR:\n"
        f"**Reaction:** ${equation}$, i.e. {reaction_string}\n\n"
        f"The reaction is carried out **isothermally** and **isobarically**. The feed contains only {reactant_A_name} (A) and {reactant_B_name} (B), with no inerts, and enters the reactor with an initial concentration of {reactant_A_name} of $C_{{A0}} = {C_A0}$ mol/dm³.\n\n"
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
        f"- Volumetric Change Parameter: $\\epsilon = {epsilon}$ (as stated; it follows from the stoichiometry and the feed: "
        f"$\\delta = (c + d - a - b)/a = ({n_products} - {a} - {b})/{a} = {delta}$ mol of total change per mol of A reacted, "
        f"$y_{{A0}} = 1/(1 + \\Theta_B) = 1/(1 + {Theta_B}) = {y_A0:.4f}$, "
        f"$\\epsilon = y_{{A0}} \\, \\delta = {y_A0:.4f} \\times ({delta}) = {eps_full:.{eps_dp}f} \\approx {epsilon}$)\n"
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
        
        f"**Answer:**\n"
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
