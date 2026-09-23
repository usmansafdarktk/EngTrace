import random
import math
from decimal import Decimal, ROUND_HALF_UP
from data.templates.branches.chemical_engineering.constants import (
    SUBSTANCES_FOR_HEATING, SUBSTANCES_FOR_VAPORIZATION,
    HEATS_OF_FORMATION, REACTIONS, CP_PARAMS, COMBUSTION_REACTIONS,
    CP_PARAMS_COMBUSTION, CP_COMBUSTION_VALID_T_MAX,
)
from data.templates.branches._emission import signed_term, joined_terms


def _as_printed(x, spec):
    """The value a reader recovers from `x` when it is printed with `spec`.

    P2 asks that the stored value and the printed value be the SAME value,
    so every downstream line is computed from the number the reader sees
    (Phase 1, DECISIONS D-016).
    """
    return float(format(x, spec))


def _hu(x, places):
    """Round half-up to `places` dp, resolving the tie in DECIMAL.

    `round()` resolves a half-way tie on the binary value and disagrees with a
    reader doing decimal arithmetic (spec P2 as amended, DECISIONS D-012).
    """
    q = Decimal(1).scaleb(-places)
    d = x if isinstance(x, Decimal) else Decimal(repr(x))
    v = d.quantize(q, rounding=ROUND_HALF_UP)
    return int(v) if places == 0 else float(v)


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


def _exact_spec(x, kind, min_places=3, max_places=12):
    """The shortest `.<p><kind>` spec, p >= min_places, that prints `x` EXACTLY.

    A table coefficient is never coarsened to its display: the display is
    lengthened until a reader who re-parses the printed digits recovers the
    float the template consumed (D-037). Every CP_PARAMS row is reached by
    p = 5; a value not reached by `max_places` is a data error, so this raises
    rather than print a coefficient the solution did not use.
    """
    for places in range(min_places, max_places + 1):
        spec = f'.{places}{kind}'
        if float(format(x, spec)) == x:
            return spec
    raise ValueError(f'{x!r} has no exact {kind!r} display within {max_places} places')


# Template 1 (Easy)
def template_sensible_heat_constant_cp():
    """
    Sensible Heat Calculation (Constant Heat Capacity)

    Scenario:
        This template calculates the sensible heat required to change the
        temperature of a substance without changing its phase. It uses the
        substance's actual specific heat capacity (Cp) for a physically
        accurate problem.

        The governing equation on a mass basis is:
            Q = m * Cp * (T2 - T1)

    Trace integrity (Layer 0, 2026-09-23):
        Q in J (m*Cp*dT, exact at 4-5 dp) is stated at 1 dp and the kJ
        answer is that stated value / 1000, quoted at 2 dp, so the kJ
        display drops one digit of the stated J value and sits exactly on
        a half-way tie whenever the J value ends in 5.0 (625.0 J -> 0.625
        kJ; ~0.8% of draws). The answer's precision is what the item
        quotes, so the tie cannot be removed by lengthening (D-044): as in
        template_sensible_heat_temp_dependent_cp, a draw whose J value or
        whose J/1000 sits on a half-way tie at its display is redrawn, and
        dT, Q in J and Q in kJ are bound through their displays so every
        line is computed from the operands it shows (D-016 part 2).

    Returns:
        tuple: A tuple containing:
            - str: A question asking to compute the sensible heat.
            - str: A step-by-step solution showing the calculation.
    """
    for _attempt in range(200):
        # 1. Parameterize the inputs using the substance list
        substance_data = random.choice(SUBSTANCES_FOR_HEATING)
        substance_name = substance_data["name"]
        substance_state = substance_data["state"]
        Cp = substance_data["Cp"]  # J/(g·K)

        # Retrieve safe temperature limits from the constant data
        # Defaults provided in case keys are missing in legacy data
        min_temp = substance_data.get("min_temp", 20)
        max_temp = substance_data.get("max_temp", 150)

        # Generate a random mass in grams
        m = round(random.uniform(50.0, 1000.0), 1)

        # Generate temperatures within the safe bounds for this specific substance
        # Ensure there's at least a 10 degree window for T2
        safe_max_T1 = max(min_temp, max_temp - 15)
        T1_C = round(random.uniform(min_temp, safe_max_T1), 1)

        # Ensure T2 is higher than T1 but within max limit
        T2_C = round(random.uniform(T1_C + 10, max_temp), 1)

        # 2. Perform the core calculation
        # dT is stated at 1 dp, so the chain consumes the stated value.
        delta_T = _as_printed(T2_C - T1_C, '.1f')
        # Q will be in Joules (g * J/g·K * K)
        Q_J = m * Cp * delta_T
        # Q is stated in J at 1 dp and the kJ answer is that stated value /
        # 1000 at 2 dp. Either display can sit EXACTLY on a half-way tie
        # (Q_J/1000 does whenever the stated J value ends in 5.0: ~0.8% of
        # draws). No rounding convention closes a tie for every reader, and
        # the kJ answer's precision is what the item quotes, so the draw is
        # repeated rather than resolved (D-016, D-044).
        if _is_display_tie(Q_J, 1):
            continue
        Q_J = _as_printed(Q_J, '.1f')
        if _is_display_tie(Q_J / 1000, 2):
            continue
        # Convert to kilojoules for the final answer: the stated J value /
        # 1000, itself bound through its 2-dp display.
        Q_kJ = _as_printed(Q_J / 1000, '.2f')
        break
    else:
        raise RuntimeError(
            "sensible_heat_constant_cp: no display-stable sample in 200 draws")

    # 3. Generate the question and solution strings
    question = (
        f"Calculate the heat in kJ required to raise the temperature of {m} g "
        f"of {substance_state} {substance_name} from {T1_C}°C to {T2_C}°C. "
        f"The specific heat capacity of {substance_name} is {Cp} J/g·K. "
        f"Assume constant pressure and no phase change occurs."
    )

    solution = (
        f"**Step 1:** State the formula.\n"
        f"Q = m * Cp * ΔT\n\n"
        f"**Note:** We assume the specific heat capacity (Cp) is constant over this temperature range.\n"
        f"\n\n"

        f"**Step 2:** List the given values.\n"
        f"- Mass (m) = {m} g\n"
        f"- Specific Heat Capacity (Cp) = {Cp} J/g·K\n"
        f"- Initial Temperature (T1) = {T1_C}°C\n"
        f"- Final Temperature (T2) = {T2_C}°C\n\n"

        f"**Step 3:** Calculate the temperature difference (ΔT).\n"
        f"A change in temperature has the same magnitude in Celsius and Kelvin.\n"
        f"ΔT = T2 - T1 = {T2_C}°C {signed_term(-T1_C)}°C = {round(delta_T, 1)} K\n\n"

        f"**Step 4:** Substitute the values into the formula to find the heat in Joules (J).\n"
        f"Q = {m} g * {Cp} J/g·K * {round(delta_T, 1)} K\n"
        f"Q = {round(Q_J, 1)} J\n\n"

        f"**Step 5:** Convert the heat to kilojoules (kJ) as requested.\n"
        f"Q = {round(Q_J, 1)} J * (1 kJ / 1000 J) = {round(Q_kJ, 2)} kJ\n\n"

        f"**Answer:** The heat required is **{round(Q_kJ, 2)} kJ**."
    )

    return question, solution


# Template 2 (Easy)
def template_latent_heat_vaporization():
    """
    Latent Heat of Vaporization

    Scenario:
        This template calculates the heat required to cause a phase change
        (vaporization) at a constant temperature and pressure. This energy,
        known as latent heat, is used to overcome intermolecular forces
        rather than to increase the substance's kinetic energy (temperature).

        The governing equation is:
            Q = n * ΔH_vap
        Where:
        - Q: Total heat absorbed
        - n: Number of moles
        - ΔH_vap: Molar heat of vaporization

    Returns:
        tuple: A tuple containing:
            - str: A question asking to compute the latent heat.
            - str: A step-by-step solution showing the calculation.
    """
    # 1. Parameterize the inputs using the substance list
    substance_data = random.choice(SUBSTANCES_FOR_VAPORIZATION)
    substance_name = substance_data["name"]
    delta_H_vap = substance_data["delta_H_vap"]  # in kJ/mol

    # Generate a random number of moles
    n = round(random.uniform(0.5, 5.0), 2)  # moles

    # 2. Perform the core calculation
    # The result will be in kJ since ΔH_vap is in kJ/mol
    Q_kJ = n * delta_H_vap

    # 3. Generate the question and solution strings
    question = (
        f"How much heat in kJ is required to completely vaporize {n} moles of "
        f"liquid {substance_name} at its normal boiling point? The molar heat of "
        f"vaporization for {substance_name} is {delta_H_vap} kJ/mol."
    )

    solution = (
        f"**Step 1:** State the formula.\n"
        f"Q = n * ΔH_vap\n\n"

        f"**Step 2:** List the given values.\n"
        f"- Moles (n) = {n} mol\n"
        f"- Molar Heat of Vaporization (ΔH_vap) = {delta_H_vap} kJ/mol\n\n"

        f"**Step 3:** Substitute the values into the formula.\n"
        f"Q = {n} mol * {delta_H_vap} kJ/mol\n\n"

        f"**Step 4:** Calculate the total heat required.\n"
        f"Q = {round(Q_kJ, 2)} kJ\n\n"

        f"**Answer:** The total heat required to vaporize the {substance_name} is **{round(Q_kJ, 2)} kJ**."
    )

    return question, solution


# Template 3 (Easy)
def template_heat_of_reaction_formation():
    """
    Standard Heat of Reaction from Heats of Formation

    Scenario:
        This template applies Hess's Law to find the standard heat of reaction
        at 298.15 K (ΔH_rxn°). It is calculated by subtracting the sum of the
        standard heats of formation (ΔH_f°) of the reactants from the sum of
        the heats of formation of the products, with each being weighted by
        its stoichiometric coefficient (v).

        The governing equation is:
            ΔH_rxn° = Σ(v * ΔH_f°)products - Σ(v * ΔH_f°)reactants

    Returns:
        tuple: A tuple containing:
            - str: A question asking to compute the standard heat of reaction.
            - str: A step-by-step solution showing the calculation.
    """
    # 1. Parameterize the inputs by choosing a random reaction
    reaction_data = random.choice(REACTIONS)
    equation = reaction_data["equation"]
    reactants = reaction_data["reactants"]
    products = reaction_data["products"]

    # 2. Perform the core calculation with full precision
    products_enthalpy = sum(nu * HEATS_OF_FORMATION[species] for species, nu in products.items())
    reactants_enthalpy = sum(nu * HEATS_OF_FORMATION[species] for species, nu in reactants.items())
        
    # Final calculation uses the precise, unrounded intermediate values
    delta_H_rxn = products_enthalpy - reactants_enthalpy

    # 3. Generate the question and solution strings
    required_species = list(reactants.keys()) + list(products.keys())
    hf_list_for_question = "\n".join(
        f"- {s}: {HEATS_OF_FORMATION[s]} kJ/mol" for s in sorted(list(set(required_species)))
    )

    # Question updated for clarity
    question = (
        f"Calculate the standard enthalpy of reaction, ΔH_rxn°, at 298.15 K (in kJ) for the reaction as written:\n"
        f"{equation}\n\n"
        f"Use the standard heats of formation (ΔH_f°) provided below:\n{hf_list_for_question}"
    )

    # Helper functions to format the solution steps
    def format_sum(species_dict):
        return " + ".join([f"({nu} * ΔH_f°[{s}])" for s, nu in species_dict.items()])

    def format_calc(species_dict):
        return " + ".join([f"({nu} * {HEATS_OF_FORMATION[s]})" for s, nu in species_dict.items()])

    # Solution updated to use correct rounding procedures and units
    solution = (
        f"**Step 1:** State the formula.\n"
        f"ΔH_rxn° = Σ(v * ΔH_f°)products - Σ(v * ΔH_f°)reactants\n\n"

        f"**Step 2:** Calculate the total enthalpy of the products.\n"
        f"Σ_products = {format_sum(products)}\n"
        f"Σ_products = {format_calc(products)}\n"
        f"Σ_products = {products_enthalpy:.3f} kJ\n\n"

        f"**Step 3:** Calculate the total enthalpy of the reactants.\n"
        f"Σ_reactants = {format_sum(reactants)}\n"
        f"Σ_reactants = {format_calc(reactants)}\n"
        f"Σ_reactants = {reactants_enthalpy:.3f} kJ\n\n"

        f"**Step 4:** Calculate the standard enthalpy of reaction using the full-precision values.\n"
        f"ΔH_rxn° = ({products_enthalpy:.3f}) - ({reactants_enthalpy:.3f})\n"
        f"ΔH_rxn° = {round(delta_H_rxn, 2)} kJ\n\n"

        f"**Answer:** The standard enthalpy of reaction is **{round(delta_H_rxn, 2)} kJ** for the reaction as written."
    )
    
    return question, solution


# Template 4 (Intermediate)
def template_sensible_heat_temp_dependent_cp():
    """
    Sensible Heat (Temperature-Dependent Heat Capacity)

    Scenario:
        This template provides a more accurate calculation of sensible heat where
        the heat capacity (Cp) is a polynomial function of temperature. The total
        heat required is found by integrating this function over the temperature
        range from T1 to T2.

        The governing equation is:
            Q = n * integral(Cp(T) dT) from T1 to T2

    Trace integrity (Layer 0, 2026-09-23):
        The heat-capacity coefficients are printed at the precision the table
        gives them, never coarser: `.3e` showed B = 39.064e-3 as 3.906e-02 and
        C = -13.301e-6 as -1.330e-05, so a reader recomputing I(T) from the
        printed digits missed the printed value by 2 K (D-037). The
        temperature limits were not the cause: round(., 2) prints exactly, so
        they were already consumed as printed. I(T2) and I(T1) are bound
        through their 4-dp display before they are subtracted, and Q in J
        through its 2-dp display before it is divided by 1000 (the kJ answer
        is likewise bound through its 2-dp display), so every printed line is
        computed from the operands it shows (D-016 part 2).
        A draw whose stated Q in J, or Q/1000, sits exactly on a 2-dp half-way
        tie is redrawn: the kJ answer's precision is what the item quotes, so
        that tie cannot be removed by lengthening (D-044).

    Returns:
        tuple: A tuple containing:
            - str: A question asking to compute the sensible heat via integration.
            - str: A step-by-step solution showing the analytical integration.
    """
    # 1. Parameterize the inputs
    R = 8.314  # J/(mol·K)

    for _attempt in range(200):
        substance_name = random.choice(list(CP_PARAMS.keys()))
        params = CP_PARAMS[substance_name]
        A, B, C, D = params["A"], params["B"], params["C"], params["D"]

        n = round(random.uniform(1.0, 10.0), 2)  # moles

        # Determine valid temperature ranges based on phase to ensure physical plausibility
        if "(l)" in substance_name:
            # Liquids: keep the range lower to avoid boiling (approximate general cap), and
            # START AT THE FIT'S OWN FLOOR. CP_PARAMS declares @domain T_lo=298..1500 K, but
            # uniform(280.0, 300.0) put T1 as low as 281.67 K - the integral evaluated 16 K
            # below the interval the polynomial was fitted over, on 60 of 300 seeds
            # (domain_findings.txt; C3.5 register item 7). The gas branch below already starts
            # at 298.15; this agrees with it.
            # The span stays ~20 K wide rather than clipping to [298.15, 300.0]: collapsing a
            # 20 K draw to 1.85 K would shrink the answer space, which t6_distribution calls a
            # downgrade rather than a fix.
            T1 = round(random.uniform(298.15, 318.0), 2)
            T2 = round(random.uniform(T1 + 20, 350.0), 2)
        else:
            # Gases/Solids: Can handle higher temperatures
            T1 = round(random.uniform(298.15, 500.0), 2)
            T2 = round(random.uniform(T1 + 100, 1200.0), 2)

        # 2. Perform the core calculation via analytical integration
        # integral(A + BT + CT² + DT⁻²)dT = AT + (B/2)T² + (C/3)T³ - D/T
        def integral_mean_cp_over_r(T):
            term_a = A * T
            term_b = (B / 2) * (T**2)
            term_c = (C / 3) * (T**3)
            term_d = -D / T if T != 0 else 0
            return term_a + term_b + term_c + term_d

        # T1 and T2 are stated to 2 dp and round(., 2) IS the stated value, so
        # the limits are consumed as printed. I(T2) and I(T1) are printed at
        # 4 dp and then subtracted, so each is bound through that display
        # first; the difference of two 4-dp values is exact at 4 dp (D-016).
        val_T2 = _as_printed(integral_mean_cp_over_r(T2), '.4f')
        val_T1 = _as_printed(integral_mean_cp_over_r(T1), '.4f')
        integral_K = _as_printed(val_T2 - val_T1, '.4f')
        integral_val = R * integral_K

        Q_J = n * integral_val  # Total heat in Joules
        # Q is stated in J at 2 dp and the kJ answer is that stated value / 1000.
        # Either 2-dp display can sit EXACTLY on a half-way tie (Q_J/1000 does
        # whenever the stated J value ends in 5.00: 0.1% of draws). No rounding
        # convention closes a tie for every reader, and the kJ answer's
        # precision is what the item quotes, so the draw is repeated rather
        # than resolved (D-016, D-044).
        if _is_display_tie(Q_J, 2):
            continue
        Q_J = _as_printed(Q_J, '.2f')
        if _is_display_tie(Q_J / 1000, 2):
            continue
        # The kJ answer is the stated J value / 1000, itself bound through its
        # 2-dp display (no tie can reach it: screened above).
        Q_kJ = _as_printed(Q_J / 1000, '.2f')   # Convert to kilojoules
        break
    else:
        raise RuntimeError(
            "sensible_heat_temp_dependent_cp: no display-stable sample in 200 draws")

    # 3. Generate the question and solution strings
    # Dynamically create the Cp/R equation string for the question
    # Table coefficients are printed at the precision they HAVE, never coarser
    # (D-037): `.3e` showed B = 39.064e-3 as 3.906e-02, and a reader who
    # recomputed I(T) from the printed digits missed the printed value by 2 K.
    # The specs are the shortest that re-parse to the consumed float, and are
    # never shorter than the original `.3e` / `.3g`.
    B_e, C_e, D_e = (_exact_spec(x, 'e') for x in (B, C, D))
    B_g, C_g, D_g = (_exact_spec(x, 'g') for x in (B, C, D))
    cp_eq_parts = [str(A)]
    if B != 0: cp_eq_parts.append(f"{B:{B_e}}*T")
    if C != 0: cp_eq_parts.append(f"{C:{C_e}}*T^2")
    if D != 0: cp_eq_parts.append(f"{D:{D_e}}*T^-2")
    
    cp_eq_str = " + ".join(cp_eq_parts).replace("+ -", "- ")

    question = (
        f"Calculate the heat in kJ required to raise the temperature of {n} moles of "
        f"{substance_name} from {T1} K to {T2} K. The molar heat capacity for "
        f"{substance_name} is given by the equation:\n"
        f"  Cp/R = {cp_eq_str}\n"
        f"Where the constants are: A={A}, B={B:{B_g}}, C={C:{C_g}}, D={D:{D_g}}"
    )

    # Helper to format the polynomial substitution string for the solution
    def format_poly_sub(T_val):
        terms = [f"{A}*({T_val})"]
        if B != 0: terms.append(f"({B:{B_e}}/2)*({T_val})^2")
        if C != 0: terms.append(f"({C:{C_e}}/3)*({T_val})^3")
        if D != 0: terms.append(f"-({D:{D_e}}/{T_val})")
        return " + ".join(terms).replace("+ -", "- ")

    solution = (
        f"**Step 1:** State the integral formula for total heat (Q).\n"
        f"Q = n * R * ∫[from {T1} to {T2}] (Cp/R) dT\n\n"

        f"**Step 2:** Perform the analytical integration.\n"
        f"The indefinite integral of Cp(T)/R is:\n"
        f"I(T) = A*T + (B/2)*T^2 + (C/3)*T^3 - D/T\n\n"

        f"**Step 3:** Evaluate the definite integral.\n"
        f"I({T2}) = {format_poly_sub(T2)} = {val_T2:.4f}\n"
        f"I({T1}) = {format_poly_sub(T1)} = {val_T1:.4f}\n"
        f"Integral Value = I({T2}) - I({T1}) = {val_T2:.4f} - {val_T1:.4f} = {integral_K:.4f} K\n\n"

        f"**Step 4:** Calculate the total heat Q.\n"
        f"Q = n * R * (Integral Value)\n"
        f"Q = {n} mol * 8.314 J/(mol K) * {integral_K:.4f} K\n"
        f"Q = {Q_J:.2f} J\n\n"

        f"**Step 5:** Convert to kilojoules.\n"
        f"Q = {Q_J:.2f} J / 1000 = {Q_kJ:.2f} kJ\n\n"

        f"**Answer:** The heat required is **{Q_kJ:.2f} kJ**."
    )

    return question, solution


# Template 5 (Advanced)
def template_adiabatic_flame_temperature():
    """
    Adiabatic Flame Temperature

    Scenario: Calculates the theoretical maximum temperature of combustion
    products for complete combustion with the theoretical amount of dry air,
    with no dissociation.

    Phase 2 notes:

        STATED ALGORITHM, NOT A SOLVER. This called `scipy.optimize.fsolve` and
        printed "we use a numerical solver" as its Step 4. That is not a
        reasoning trace: it is reproducible only to the solver's tolerance, it
        depends on the SciPy build, and there is no arithmetic in it for a
        reader to follow or a grader to check. It is replaced by the
        mean-heat-capacity iteration from Smith-Van Ness-Abbott - the same
        source CP_PARAMS follows - with a stated initial guess, a fixed
        iteration count, and every intermediate printed:

            T_(k+1) = T0 + (-dH_rxn) / SUM_i n_i <Cp>_i(T0, T_k)

        where the mean heat capacity between T0 and T has the closed form

            <Cp>/R = A + (B/2)(T + T0) + (C/3)(T^2 + T*T0 + T0^2) + D/(T*T0)

        Six iterations from a 2000 K guess bring every one of the eleven
        reactions within 0.035 K of its fixed point, so the answer is quoted to
        the nearest kelvin. Reviewer C measured the tightest rounding-boundary
        margin at 0.116 K (propane) - a 3.3x safety factor - and confirmed the
        map is a local contraction, |g'| about 0.09-0.135, with the same
        answer for any initial guess from 400 K to 10000 K.

        There IS a convergence assert below, and it stays: the iteration count
        being fixed is what makes the trace reproducible, while the assert is
        what stops a silently unconverged answer if the reaction set ever
        changes. Note it is stripped by `python -O`, so it guards development,
        not production (Reviewer C, C-6).

        NO SILENT FALLBACKS. The retry loop caught `Exception` - including the
        solver's own failures - and silently tried a different reaction, so a
        broken calculation became a different question rather than an error.
        Its `while/else` returned a tuple carrying neither steps nor an answer,
        which is not a benchmark item at all. Both are gone: the data
        preconditions are checked explicitly and a violation raises (P5).
        Measured over 2000 seeds, the loop never retried once.

        HIGH-TEMPERATURE HEAT CAPACITIES. This template reaches 2908 K
        (acetylene) while CP_PARAMS is fitted to 1500 K, so it was
        extrapolating ~1400 K past validity and every flame temperature was
        low. It now reads CP_PARAMS_COMBUSTION, fitted over 298-3000 K
        (DECISIONS D-032, D-036), which leaves 92 K of headroom under the
        ceiling. The figure 2844 K appeared here before and was measured with
        the OLD table - correcting the table raised the temperatures, so the
        stale number understated the margin it was describing (Reviewer C,
        C-5).

        PRECISION OF THE ANSWER. Quoted to the nearest kelvin, but the model is
        not good to a kelvin: against a direct NIST Shomate solve the answers
        differ by up to 1.93 K, and on 6 of the 11 reactions the nearest-kelvin
        value differs. The trace is internally exact - it reproduces its own
        arithmetic - so this is a statement about the polynomial, not the
        solver. Step 6 now says so, because a grader marking to +/-1 K would
        otherwise fail a solver who used NIST data directly (Reviewer C, C-2).
    """
    R = 8.314           # J/(mol K)
    T_initial = 298.15  # K
    initial_guess = 2000.0
    n_iterations = 6
    dp = 2              # display precision for the intermediates

    reaction_data = random.choice(COMBUSTION_REACTIONS)
    fuel = reaction_data["fuel"]
    equation = reaction_data["equation"]
    reactants = reaction_data["reactants"]
    products = reaction_data["products"]

    # Preconditions, checked and raised rather than caught and skipped. If a
    # reaction names a species with no thermochemical data, that is a defect in
    # constants.py and the benchmark should stop, not quietly ask a different
    # question (P5).
    missing_hf = [s for s in list(reactants) + list(products)
                  if s not in HEATS_OF_FORMATION]
    assert not missing_hf, (
        f"{fuel}: no heat of formation for {missing_hf}")
    missing_cp = [s for s in products if s not in CP_PARAMS_COMBUSTION]
    assert not missing_cp, (
        f"{fuel}: no high-temperature heat capacity for {missing_cp}")

    # Standard heat of reaction at 298.15 K, bound through its display so the
    # iteration below is computed from the number the reader sees.
    products_enthalpy_298 = sum(nu * HEATS_OF_FORMATION[s]
                                for s, nu in products.items())
    reactants_enthalpy_298 = sum(nu * HEATS_OF_FORMATION[s]
                                 for s, nu in reactants.items())
    delta_H_298_kJ = _as_printed(products_enthalpy_298 - reactants_enthalpy_298,
                                 f'.{dp}f')
    delta_H_298_J = delta_H_298_kJ * 1000.0

    assert delta_H_298_kJ < 0.0, (
        f"{fuel}: combustion must be exothermic, got {delta_H_298_kJ} kJ")

    def mean_cp_over_R(params, T):
        """<Cp>/R for one species, averaged over T_initial to T."""
        return (params["A"]
                + (params["B"] / 2) * (T + T_initial)
                + (params["C"] / 3) * (T * T + T * T_initial + T_initial ** 2)
                + (params["D"] / (T * T_initial) if params["D"] else 0.0))

    # The iteration. Each step's divisor is bound through its display, so the
    # printed line closes exactly (P1/P2).
    T_k = initial_guess
    iterations = []
    for _ in range(n_iterations):
        cp_mixture = _as_printed(
            sum(nu * R * mean_cp_over_R(CP_PARAMS_COMBUSTION[s], T_k)
                for s, nu in products.items()), f'.{dp}f')
        T_next = _as_printed(T_initial + (-delta_H_298_J) / cp_mixture,
                             f'.{dp}f')
        iterations.append((T_k, cp_mixture, T_next))
        T_k = T_next

    adiabatic_temp_kelvin = round(T_k)

    # --- invariants (T7) ---------------------------------------------------
    assert all(cp > 0.0 for _, cp, _ in iterations), (
        f"{fuel}: non-physical mixture heat capacity in the iteration")
    assert abs(iterations[-1][2] - iterations[-1][0]) < 0.5, (
        f"{fuel}: iteration not converged after {n_iterations} passes: "
        f"{iterations[-1][0]} -> {iterations[-1][2]}")
    assert T_initial < adiabatic_temp_kelvin <= CP_COMBUSTION_VALID_T_MAX, (
        f"{fuel}: flame temperature {adiabatic_temp_kelvin} K outside the "
        f"{CP_COMBUSTION_VALID_T_MAX} K validity of CP_PARAMS_COMBUSTION")

    question = (
        f"{fuel} gas enters a furnace at {T_initial} K and is burned completely with "
        f"the theoretical amount of dry air (also at {T_initial} K). Assuming the "
        f"process is adiabatic and there is no shaft work, estimate the "
        f"adiabatic flame temperature in Kelvin."
    )

    solution = (
        f"**Step 1:** Write the balanced chemical equation including inert nitrogen from air.\n"
        f"  {equation}\n\n"

        f"**Step 2:** Calculate the standard heat of reaction at {T_initial} K (ΔH_rxn°).\n"
        f"Using standard heats of formation:\n"
        f"ΔH_rxn° = Σ n·ΔHf°(products) − Σ n·ΔHf°(reactants)\n"
        f"ΔH_rxn° = ({round(products_enthalpy_298, dp)}) − ({round(reactants_enthalpy_298, dp)}) = {delta_H_298_kJ} kJ\n\n"

        f"**Step 3:** Set up the energy balance.\n"
        f"For an adiabatic process the heat released by the reaction is absorbed as sensible heat by the products:\n"
        f"  Σ n_i · <Cp>_i · (T_ad − {T_initial}) + ΔH_rxn° = 0\n"
        f"so   T_ad = {T_initial} + ({-delta_H_298_J}) / Σ n_i·<Cp>_i\n\n"

        f"**Step 4:** State the solution method.\n"
        f"Σ n_i·<Cp>_i depends on T_ad, so the equation is solved by direct iteration:\n"
        f"guess T, evaluate the mean heat capacities over [{T_initial}, T], and recompute T.\n"
        f"The mean heat capacity of a species between {T_initial} and T is\n"
        f"  <Cp>/R = A + (B/2)(T + T0) + (C/3)(T² + T·T0 + T0²) + D/(T·T0)\n"
        f"Initial guess: T = {initial_guess} K. Iterations performed: {n_iterations}.\n\n"

        f"**Step 5:** Iterate.\n"
    )

    for i, (T_in, cp_mixture, T_out) in enumerate(iterations, start=1):
        solution += (
            f"Pass {i}: T = {T_in:.{dp}f} K → Σ n_i·<Cp>_i = {cp_mixture:.{dp}f} J/K\n"
            f"  T = {T_initial} + {-delta_H_298_J} / {cp_mixture:.{dp}f} = {T_out:.{dp}f} K\n"
        )

    solution += (
        f"\n**Step 6:** State the result.\n"
        f"Successive passes agree to within "
        f"{round(abs(iterations[-1][2] - iterations[-1][0]), dp)} K, so the "
        f"iteration has converged at this display precision.\n"
        f"T_ad = {adiabatic_temp_kelvin} K\n"
        f"The iteration is exact for this heat-capacity model, but the model "
        f"itself is good to about ±2 K against NIST reference data, so the "
        f"answer should be read as {adiabatic_temp_kelvin} ± 2 K.\n\n"

        f"**Answer:** The estimated adiabatic flame temperature is **{adiabatic_temp_kelvin} K**."
    )

    return question, solution


def main():
    """
    Generate numerous instances of each heat effects template with different random seeds
    and write the results to a JSONL file.
    """
    import json
    import os

    # Define the output path (Modify this path according to where you are running the code from)
    output_file = "testset/chemical_engineering/thermodynamics/heat_effects.jsonl"

    # Create the directory if it doesn't exist
    os.makedirs(os.path.dirname(output_file), exist_ok=True)

    # List of template functions with their ID and level
    templates = [
        (template_sensible_heat_constant_cp, "sensible_heat_constant_cp", "Easy"),
        (template_latent_heat_vaporization, "latent_heat_vaporization", "Easy"),
        (template_heat_of_reaction_formation, "heat_of_reaction_formation", "Easy"),
        (template_sensible_heat_temp_dependent_cp, "sensible_heat_temp_dependent_cp", "Intermediate"),
        (template_adiabatic_flame_temperature, "adiabatic_flame_temperature", "Advanced"),
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
                "domain": "thermodynamics",
                "area": "heat_effects",
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
