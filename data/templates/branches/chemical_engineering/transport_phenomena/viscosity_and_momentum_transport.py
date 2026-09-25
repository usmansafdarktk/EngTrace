import random
import math
from decimal import Decimal, ROUND_HALF_UP
from data.templates.branches.chemical_engineering.constants import COMMON_LIQUIDS, COMMON_GASES, GAS_MOLECULAR_PARAMS, POWER_LAW_FLUIDS


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


def _decimals(x):
    """Decimal places of the shortest decimal that round-trips `x`.

    0.1 -> 1, 0.012 -> 3, 4.003 -> 3, 32.0 -> 0. Used to print a quantity
    built from table values at the precision it actually has (D-037), never
    coarser.
    """
    exponent = Decimal(repr(float(x))).normalize().as_tuple().exponent
    return max(0, -exponent)


def _is_sci_display_tie(x, mant_dp):
    """Is `x` on a half-way tie when printed with `.{mant_dp}e`?

    The rounding position of a scientific display moves with the decade, so
    the fixed-place test is applied at `mant_dp` places below the leading
    digit (D-016).
    """
    if x == 0:
        return False
    return _is_display_tie(x, mant_dp - math.floor(math.log10(abs(x))))


# Template 1 (Easy)
def template_newtons_law_shear_stress():
    """
    Shear Stress and Force Calculation for Flow Between Parallel Plates

    Scenario:
        This template models Couette flow, a classic fluid dynamics problem where a layer
        of fluid is sheared between two parallel plates. One plate is stationary, and
        the other moves at a constant velocity. Assuming a linear velocity profile,
        we use Newton's Law of Viscosity to find the shear stress in the fluid and
        the force required to move the plate.

        The relevant equations are:
            - Velocity Gradient: dvx/dy = V / Y
            - Shear Stress: tau_yx = mu * (dvx/dy)
            - Force: F = tau_yx * A

    Layer 2 fix (2026-09-25):
        All three experts found the question unsolvable as stated: the
        viscosity the solution uses appeared only in the solution's "Given
        Information", and literature values for the table's fluids differ by
        20% or more (plasma 1.1-1.6 mPa.s), so a reader could not reach the
        gold answer. Confirmed on seeds 2101-2105. The question now states
        mu, with three significant figures as before, or four where the table
        value has them (water 1.002e-3, ethanol 1.074e-3, mercury 1.526e-3)
        and the value used is the displayed one (D-016 part 2, D-037). The
        same experts found tau and F printed at a fixed 3 dp (instance 4:
        tau = 0.0265 Pa printed as 0.027 Pa) and Y printed as a raw float
        ("0.34 cm = 0.0034000000000000002 m"); Y is now bound through its
        exact 4-dp metre display, and tau and F are printed in scientific
        notation at their exact length (never shorter than a 4-dp mantissa):
        with dvx/dy bound at 2 dp they are products of finite decimals, and a
        fixed-width display sat on a half-way tie for half the draws of a
        2-sf viscosity (measured 10.6% of seeds redrawn, 52% for diesel), so
        the display is lengthened rather than the draw resampled (D-037).
        Every displayed-and-consumed operand (Y, dvx/dy, tau) is bound
        through its display, so each printed line is the product of the
        printed operands; a draw whose dvx/dy (a 2-dp over 4-dp quotient)
        sits on a half-way tie at 2 dp is redrawn (D-016; about 1% of
        seeds, not concentrated by fluid). Every gold answer changes in
        format; the question text changes on every seed (mu is added).

    Returns:
        tuple: A tuple containing:
            - str: A question asking to compute shear stress and force.
            - str: A step-by-step solution showing the calculations.
    """
    # 1. Parameterize the inputs with random values
    fluid_name, (density, viscosity) = random.choice(list(COMMON_LIQUIDS.items()))
    # The viscosity is STATED (Layer 2: it was not, and the question could not
    # be solved) with a 3-significant-figure display, or four where the table
    # value has them; the value used is the displayed one (D-016 part 2).
    mu_spec = '.2e' if float(format(viscosity, '.2e')) == viscosity else '.3e'
    viscosity = _as_printed(viscosity, mu_spec)

    for _attempt in range(200):
        V = round(random.uniform(0.1, 2.0), 2)
        Y_cm = round(random.uniform(0.1, 2.5), 2)
        # A 2-dp centimetre value is exactly 4 dp in metres; the float is not
        # (0.34/100 = 0.0034000000000000002), so Y is bound through that display.
        Y_m = _as_printed(Y_cm / 100, '.4f')
        A = round(random.uniform(0.5, 5.0), 2)

        # 2. Perform the core calculation. dvx/dy is printed at 2 dp and then
        # consumed, so it is bound through that display. A quotient of a 2-dp
        # by a 4-dp value can sit exactly on a half-way tie at 2 dp; such a
        # draw is redrawn rather than resolved (D-016).
        if _is_display_tie(V / Y_m, 2):
            continue
        velocity_gradient = _as_printed(V / Y_m, '.2f')
        break
    else:
        raise RuntimeError(
            "newtons_law_shear_stress: no display-stable sample in 200 draws")

    # tau = mu * dvx/dy and F = tau * A are products of finite decimals (a
    # 3-4 sf viscosity, a 2-dp gradient, a 2-dp area), so they are exact at a
    # known length and are printed at that length in scientific notation,
    # never shorter than a 4-dp mantissa: a fixed-width display would sit on
    # a half-way tie for half the draws of a 2-sf viscosity (D-037: lengthen
    # an exact display rather than resample). Both are at most 13 sf, so the
    # float round-trips its own display exactly.
    tau_dec = Decimal(repr(viscosity)) * Decimal(repr(velocity_gradient))
    tau_dp = max(4, len(tau_dec.normalize().as_tuple().digits) - 1)
    shear_stress = float(tau_dec)
    assert _as_printed(shear_stress, f'.{tau_dp}e') == shear_stress
    F_dec = tau_dec * Decimal(repr(A))
    F_dp = max(4, len(F_dec.normalize().as_tuple().digits) - 1)
    force = float(F_dec)
    assert _as_printed(force, f'.{F_dp}e') == force

    # 3. Generate the question and solution strings
    question = (
        f"Two large parallel plates with an area of {A} m^2 each are separated by a "
        f"thin film of {fluid_name} that is {Y_cm} cm thick. The top plate is moved "
        f"at a constant velocity of {V} m/s, while the bottom plate is held stationary. "
        f"The dynamic viscosity of {fluid_name} is {viscosity:{mu_spec}} Pa·s.\n\n"
        f"Assuming the fluid exhibits Newtonian behavior and a linear velocity profile, calculate:\n"
        f"a) The shear stress (tau_yx) exerted on the fluid.\n"
        f"b) The total force (F) required to move the top plate at the given velocity."
    )

    solution = (
        f"**Given Information:**\n"
        f"- Fluid: {fluid_name}\n"
        f"- Dynamic Viscosity of {fluid_name} (mu): {viscosity:{mu_spec}} Pa·s\n"
        f"- Plate Area (A): {A} m^2\n"
        f"- Plate Velocity (V): {V} m/s\n"
        f"- Distance between plates (Y): {Y_cm} cm = {Y_m:.4f} m\n\n"

        f"**Step 1:** Calculate the velocity gradient (dvx/dy).\n"
        f"For a linear velocity profile between a stationary and a moving plate, the gradient is constant:\n"
        f"dvx/dy = V / Y\n"
        f"dvx/dy = {V} m/s / {Y_m:.4f} m = {velocity_gradient:.2f} 1/s\n\n"

        f"**Step 2:** Calculate the shear stress (tau_yx).\n"
        f"Using Newton's Law of Viscosity:\n"
        f"tau_yx = mu * (dvx/dy)\n"
        f"tau_yx = ({viscosity:{mu_spec}} Pa·s) * ({velocity_gradient:.2f} 1/s)\n"
        f"tau_yx = {shear_stress:.{tau_dp}e} Pa\n\n"

        f"**Step 3:** Calculate the force (F).\n"
        f"Force is the shear stress acting over the entire area of the plate:\n"
        f"F = tau_yx * A\n"
        f"F = ({shear_stress:.{tau_dp}e} Pa) * ({A} m^2)\n"
        f"F = {force:.{F_dp}e} N\n\n"

        f"**Answer:**\n"
        f"a) The shear stress in the fluid is **{shear_stress:.{tau_dp}e} Pa**.\n"
        f"b) The force required to move the plate is **{force:.{F_dp}e} N**."
    )

    return question, solution


# Template 2 (Easy)
def template_kinematic_viscosity():
    """
    Kinematic Viscosity Calculation

    Scenario:
        This template tests the definition of kinematic viscosity (ν), which is the
        ratio of the dynamic viscosity (μ) to the density (ρ) of a fluid. It is a
        measure of a fluid's internal resistance to flow under gravitational forces.

        The relevant equations are:
            - Kinematic Viscosity: ν = μ / ρ
            - Unit Conversion: 1 m²/s = 10,000 Stokes (St)

    Returns:
        tuple: A tuple containing:
            - str: A question asking to compute the kinematic viscosity.
            - str: A step-by-step solution showing the calculation.
    """
    # 1. Parameterize the inputs with random values
    
    # Combine liquids and gases into one dictionary for random selection
    all_fluids = {**COMMON_LIQUIDS, **COMMON_GASES}
    fluid_name, (density, viscosity) = random.choice(list(all_fluids.items()))
    
    # Randomly choose the target units for the answer
    target_units = random.choice(["m^2/s", "Stokes (St)"])

    # 2. Perform the core calculation
    
    # Calculate kinematic viscosity in SI units (m²/s)
    nu_si = viscosity / density
    
    # Perform unit conversion if necessary
    final_nu = nu_si
    if "Stokes" in target_units:
        final_nu = nu_si * 10000

    # 3. Generate the question and solution strings
    
    question = (
        f"The dynamic viscosity (μ) of {fluid_name} at standard conditions is approximately "
        f"{viscosity:.2e} Pa·s, and its density (ρ) is {density:.3f} kg/m³.\n\n"
        f"Calculate the kinematic viscosity (ν) of {fluid_name} in units of **{target_units}**."
    )

    solution = (
        f"**Given Information:**\n"
        f"- Fluid: {fluid_name}\n"
        f"- Dynamic Viscosity (μ): {viscosity:.2e} Pa·s\n"
        f"- Density (ρ): {density:.3f} kg/m³\n\n"
        
        f"**Step 1:** State the formula for kinematic viscosity.\n"
        f"Kinematic viscosity (ν) is defined as the ratio of dynamic viscosity to density:\n"
        f"ν = μ / ρ\n\n"
        f"\n\n"
        
        f"**Step 2:** Calculate the kinematic viscosity in SI units (m²/s).\n"
        f"Substitute the given values into the formula:\n"
        f"ν = ({viscosity:.2e} Pa·s) / ({density:.3f} kg/m³)\n"
        f"ν = {nu_si:.3e} m²/s\n\n"
    )
    
    # Add the unit conversion step only if needed
    if "Stokes" in target_units:
        # Fixed formatting here to avoid 0.0000 for small values
        solution += (
            f"**Step 3:** Convert the result to Stokes (St).\n"
            f"The conversion factor is 1 m²/s = 10,000 St.\n"
            f"ν = ({nu_si:.3e} m²/s) * (10,000 St / 1 m²/s)\n"
            f"ν = {final_nu:.4e} St\n\n"
        )
        
    solution += (
        f"**Answer:**\n"
        # Fixed formatting here to avoid 0.0000 for small values
        f"The kinematic viscosity of {fluid_name} is **{final_nu:.4e} {target_units}**."
    )

    return question, solution


# Template 3 (Intermediate)
def template_gas_viscosity_kinetic_theory():
    """
    Gas Viscosity Estimation from Kinetic Theory

    Scenario:
        This template uses the Chapman-Enskog equation, derived from the kinetic
        theory of gases, to estimate the dynamic viscosity (μ) of a low-density gas.
        This model connects a macroscopic property (viscosity) to molecular
        properties (molar mass, size).

        The relevant equation is:
            μ = (2.6693e-6 * sqrt(M*T)) / (σ² * Ω_μ)
        
        Where:
            - μ = viscosity in Pa·s
            - M = Molar Mass in g/mol
            - T = Absolute Temperature in K
            - σ = Lennard-Jones molecular diameter in Angstroms (Å)
            - Ω_μ = Collision integral (dimensionless)

    Trace integrity (Layer 0, 2026-09-23):
        sigma^2 and M*T derive from table values alone, so they are printed at
        the precision they have and never coarser (D-037): sigma^2 at
        2*dec(sigma) dp (13.771521, not 13.772 - the denominator line did not
        close from the 3-dp value on 22% of draws) and M*T at dec(M) dp (M has
        3 dp for helium and hydrogen). The numerator (4-sf) and denominator
        (3 dp) are bound through their displays before the division, so the
        printed viscosity is the quotient of the printed operands (D-016 part
        2). sigma^2 * Omega is exact at 3 dp more than sigma^2, so its 3-dp
        display sits exactly on a half-way tie for four (gas, Omega) pairs -
        xenon and ammonia (2-dp sigma^2) at Omega = 0.950 and 1.050, found by
        exhaustive search of the 17 x 101 input grid; such a draw is redrawn
        rather than resolved (D-016), since the exact display would need 9 dp.
        (Superseded in part by the Layer 2 fix below: Omega is no longer on a
        101-point grid, and the same 3-dp tie screen on sigma^2 * Omega now
        fires only where the 4-dp Omega happens to make the product a tie.)

    Layer 2 fix (2026-09-25):
        All three experts found the collision integral sampled uniformly in
        0.95-1.05, independent of T and of the gas's epsilon/k (loaded from
        the table but never used), although Omega_mu is a function of the
        reduced temperature T* = kT/epsilon. Confirmed on the hand-check
        instance: propane at 338 K has T* = 338/237.1 = 1.4256 and
        Omega_mu = 1.3430, so mu = 9.26e-06 Pa.s, not the 1.303e-05 printed
        from the sampled 0.955; helium at 518 K (T* = 50.7) has Omega_mu =
        0.648, mu = 2.88e-05, not 1.926e-05; argon at 559 K has 0.895, not
        1.042. Omega_mu is now evaluated from T* with the Neufeld-Janzen-Aziz
        correlation for the Lennard-Jones (6-12) potential,
            Omega_mu = 1.16145/T*^0.14874 + 0.52487 exp(-0.7732 T*)
                       + 2.16178 exp(-2.43787 T*),
        valid for 0.3 <= T* <= 100 (Neufeld, Janzen & Aziz, J. Chem. Phys.
        57, 1100 (1972), eq. for Omega(2,2)*; also Bird, Stewart & Lightfoot,
        Transport Phenomena 2nd ed., eq. 1.4-14 and Table E.2), and every
        sampled T* (250-600 K over the table's epsilon/k of 10.22-558.3 K)
        lies in that range. The question keeps its structure: it still GIVES
        Omega_mu, now the computed value at 4 dp, which is the smaller change
        (a reader is still asked to apply the Chapman-Enskog equation, not to
        evaluate the correlation), and the solution's given-information line
        states the provenance (T* from epsilon/k, the correlation). T* is
        bound through its 4-dp display and Omega_mu is computed from the
        bound T*, then bound through its own 4-dp display, so the printed
        Omega follows from the printed T*; a T* whose 4-dp display sits on a
        half-way tie (an integer over a 1-2 dp value) is redrawn (D-016).
        Every gold answer and every question text changes (Omega_mu is a
        different number on every seed).

    Returns:
        tuple: A tuple containing:
            - str: A question asking to estimate the gas viscosity.
            - str: A step-by-step solution showing the calculation.
    """
    # 1. Parameterize the inputs with random values
    for _attempt in range(200):
        gas_name, (molar_mass, sigma, epsilon) = random.choice(list(GAS_MOLECULAR_PARAMS.items()))

        # Temperature in Kelvin
        temperature_K = random.randint(250, 600)

        # Collision integral (dimensionless) from the reduced temperature
        # T* = kT/epsilon with the Neufeld-Janzen-Aziz correlation (Layer 2:
        # it used to be sampled uniformly in 0.95-1.05, independent of T and
        # epsilon/k, so the given Omega contradicted the stated gas and T).
        # T* is printed at 4 dp and consumed, so it is bound through that
        # display; a draw whose display sits on a half-way tie is redrawn.
        if _is_display_tie(temperature_K / epsilon, 4):
            continue
        T_star = _as_printed(temperature_K / epsilon, '.4f')
        assert 0.3 <= T_star <= 100.0, f"T* outside the correlation's range: {T_star}"
        omega_mu = _as_printed(
            1.16145 / T_star ** 0.14874
            + 0.52487 * math.exp(-0.7732 * T_star)
            + 2.16178 * math.exp(-2.43787 * T_star), '.4f')

        # 2. Perform the core calculation

        # M*T and sigma^2 derive from table values alone, so they are printed
        # at the precision they have and never coarser (D-037): M has 3 dp for
        # helium and hydrogen, and sigma^2 has twice sigma's decimals.
        mt_dp = max(2, _decimals(molar_mass))
        MT = _as_printed(molar_mass * temperature_K, f'.{mt_dp}f')
        s2_dp = max(3, 2 * _decimals(sigma))
        sigma2 = _as_printed(sigma ** 2, f'.{s2_dp}f')

        # Numerator of the Chapman-Enskog equation, bound through its 4-sf display
        numerator = _as_printed(2.6693e-6 * math.sqrt(MT), '.4e')

        # Denominator of the Chapman-Enskog equation. sigma^2 * Omega is exact
        # at s2_dp + 3 dp and is displayed at 3 dp; that display sits exactly
        # on a half-way tie for four (gas, Omega) pairs - xenon and ammonia
        # (2-dp sigma^2) at Omega = 0.950 and 1.050 - and such a draw is
        # repeated rather than resolved (D-016).
        if _is_display_tie(sigma2 * omega_mu, 3):
            continue
        denominator = _as_printed(sigma2 * omega_mu, '.3f')

        # Final viscosity calculation, from the printed operands
        viscosity = numerator / denominator
        break
    else:
        raise RuntimeError(
            "gas_viscosity_kinetic_theory: no display-stable sample in 200 draws")

    # 3. Generate the question and solution strings
    
    question = (
        f"Estimate the dynamic viscosity (μ) of {gas_name} gas at a temperature of {temperature_K} K.\n\n"
        f"Use the Chapman-Enskog equation for low-density gases:\n"
        f"μ = (2.6693e-6 * sqrt(M*T)) / (σ^2 * Ω_μ)\n\n"
        f"You are given the following parameters for {gas_name}:\n"
        f"- Molar Mass (M) = {molar_mass} g/mol\n"
        f"- Lennard-Jones diameter (σ) = {sigma} Å\n"
        f"- Collision Integral (Ω_μ) = {omega_mu:.4f}\n\n"
        f"Provide your answer in units of Pascal-seconds (Pa·s)."
    )

    solution = (
        f"**Given Information:**\n"
        f"- Gas: {gas_name}\n"
        f"- Temperature (T): {temperature_K} K\n"
        f"- Molar Mass (M): {molar_mass} g/mol\n"
        f"- Lennard-Jones diameter (σ): {sigma} Å\n"
        f"- Collision Integral (Ω_μ): {omega_mu:.4f} "
        f"(the Lennard-Jones value at the reduced temperature T* = kT/ε = {temperature_K}/{epsilon} = {T_star:.4f}, "
        f"with ε/k = {epsilon} K, from the Neufeld-Janzen-Aziz correlation)\n\n"

        f"**Step 1:** State the Chapman-Enskog equation.\n"
        f"The formula for estimating the viscosity of a low-density gas is:\n"
        f"μ = (2.6693e-6 * sqrt(M*T)) / (σ^2 * Ω_μ)\n\n"

        f"**Step 2:** Substitute the given values into the equation.\n"
        f"μ = (2.6693e-6 * sqrt({molar_mass} * {temperature_K})) / (({sigma})^2 * {omega_mu:.4f})\n\n"

        f"**Step 3:** Calculate the numerator and denominator.\n"
        f"- Numerator = 2.6693e-6 * sqrt({MT:.{mt_dp}f}) = {numerator:.4e}\n"
        f"- Denominator = {sigma2:.{s2_dp}f} * {omega_mu:.4f} = {denominator:.3f}\n\n"
        
        f"**Step 4:** Calculate the final viscosity.\n"
        f"μ = {numerator:.4e} / {denominator:.3f}\n"
        f"μ = {viscosity:.3e} Pa·s\n\n"
        
        f"**Answer:**\n"
        f"The estimated dynamic viscosity of {gas_name} at {temperature_K} K is **{viscosity:.3e} Pa·s**."
    )

    return question, solution


# Template 4 (Intermediate)
def template_reynolds_number_flow_regime():
    """
    Reynolds Number and Flow Regime Calculation

    Scenario:
        This template introduces the Reynolds number (Re), the dimensionless ratio of
        inertial forces to viscous forces in a fluid. Its value is critical for
        predicting whether a flow is smooth (laminar) or chaotic (turbulent).
        The calculation depends on a characteristic length, which varies with geometry.

        The relevant equation is:
            Re = (ρ * v * L) / μ
        
        Where:
            - ρ = density
            - v = average velocity
            - L = characteristic length (e.g., pipe diameter or plate length)
            - μ = dynamic viscosity

    Screen pass 1 (2026-09-23):
        One judge (of three) reported a calculation discrepancy in one
        instance; on the three screened seeds and at 500 seeds Re follows
        from the printed operands to the half-unit, so that claim is
        rejected. Verification found a real defect T1 could not see, because
        the substitution and the result sat on separate lines: the viscosity
        is stated with a 3-significant-figure `.2e` display, but three table
        values (water 1.002e-3, ethanol 1.074e-3, mercury 1.526e-3) carry
        four figures and were consumed at full precision, so the printed Re
        was 0.2-0.4% off the reader's on 4.8% of draws. Such a value is now
        displayed with four figures and the value used is the displayed one
        (D-016 part 2, D-037); Step 2 prints the substitution and the result
        on one line so T1 and the census cover it; a draw whose Re sits on a
        half-way tie at the integer is redrawn (D-016). For a kept draw no
        gold value changes and the question text changes only for those
        three fluids. The tie screen redraws 1.5% of draws, concentrated
        (D-045) on the fluids whose rho/mu terminates: whole milk over a
        plate loses 41% of its draws, kerosene over a plate 27%, propane,
        whole blood and kerosene in a pipe 12-15%.

    Returns:
        tuple: A tuple containing:
            - str: A question asking to compute the Reynolds number and find the flow regime.
            - str: A step-by-step solution.
    """
    # 1. Parameterize the inputs with random values
    all_fluids = {**COMMON_LIQUIDS, **COMMON_GASES}
    fluid_name, (density, viscosity) = random.choice(list(all_fluids.items()))
    # The viscosity is STATED with a 3-significant-figure `.2e` display, but
    # three table values (water 1.002e-3, ethanol 1.074e-3, mercury 1.526e-3)
    # carry four figures and were consumed at full precision, so the printed
    # Re did not follow from the printed operands (0.2-0.4% off on 4.8% of
    # draws). A value the 3-figure display cannot hold is displayed with four
    # figures, and the value used is the displayed one (D-016 part 2, D-037).
    mu_spec = '.2e' if float(format(viscosity, '.2e')) == viscosity else '.3e'
    viscosity = _as_printed(viscosity, mu_spec)

    for _attempt in range(200):
        # Randomly choose a flow geometry
        geometry = random.choice(['pipe', 'flat plate'])

        if geometry == 'pipe':
            # Pipe flow parameters
            characteristic_length = round(random.uniform(0.01, 0.5), 3) # Diameter in m
            velocity = round(random.uniform(0.1, 5.0), 2)
            length_symbol = "D"
            length_name = "diameter"
            scenario_description = f"{fluid_name} is flowing through a smooth circular pipe with an internal {length_name} of {characteristic_length} m."

        else: # flat plate
            # Flat plate flow parameters
            characteristic_length = round(random.uniform(0.1, 2.0), 2) # Length in m
            velocity = round(random.uniform(1.0, 20.0), 1)
            length_symbol = "L"
            length_name = "length"
            scenario_description = f"A flow of {fluid_name} moves over a smooth, thin flat plate with a {length_name} of {characteristic_length} m."

        # 2. Perform the core calculation. Re is printed at the integer; a draw
        # whose value sits exactly on a half-way tie there (possible only for a
        # terminating viscosity such as 8.00e-6) is redrawn (D-016).
        reynolds_number = (density * velocity * characteristic_length) / viscosity
        if _is_display_tie(reynolds_number, 0):
            continue
        break
    else:
        raise RuntimeError(
            "reynolds_number_flow_regime: no display-stable sample in 200 draws")

    # Determine flow regime based on geometry
    if geometry == 'pipe':
        if reynolds_number < 2300:
            regime = "laminar"
            regime_explanation = f"Since Re = {reynolds_number:,.0f} is less than 2300, the flow is in the **laminar** regime."
        elif reynolds_number < 4000:
            regime = "transitional"
            regime_explanation = f"Since Re = {reynolds_number:,.0f} is between 2300 and 4000, the flow is in the **transitional** regime."
        else:
            regime = "turbulent"
            regime_explanation = f"Since Re = {reynolds_number:,.0f} is greater than 4000, the flow is in the **turbulent** regime."
    else: # flat plate
        if reynolds_number < 500000:
            regime = "laminar"
            regime_explanation = f"Since Re = {reynolds_number:,.0f} is less than the critical value of 500,000, the flow is considered **laminar**."
        else:
            regime = "turbulent"
            regime_explanation = f"Since Re = {reynolds_number:,.0f} is greater than the critical value of 500,000, the flow is considered **turbulent**."

    # 3. Generate the question and solution strings
    question = (
        f"{scenario_description} The average velocity of the flow is {velocity} m/s.\n\n"
        f"The properties of {fluid_name} are:\n"
        f"- Density (ρ) = {density} kg/m³\n"
        f"- Dynamic Viscosity (μ) = {viscosity:{mu_spec}} Pa·s\n\n"
        f"Based on this information:\n"
        f"a) Calculate the Reynolds number (Re).\n"
        f"b) Determine the flow regime (laminar, transitional, or turbulent)."
    )

    solution = (
        f"**Given Information:**\n"
        f"- Fluid: {fluid_name}\n"
        f"- Density (ρ): {density} kg/m³\n"
        f"- Dynamic Viscosity (μ): {viscosity:{mu_spec}} Pa·s\n"
        f"- Average Velocity (v): {velocity} m/s\n"
        f"- Characteristic Length ({length_symbol}): {characteristic_length} m ({length_name} of the {geometry})\n\n"
        
        f"**Step 1:** State the formula for the Reynolds number.\n"
        f"The Reynolds number is the ratio of inertial forces to viscous forces.\n"
        f"Re = (ρ * v * {length_symbol}) / μ\n\n"
        
        f"**Step 2:** Substitute the values and calculate Re.\n"
        f"Re = ({density} * {velocity} * {characteristic_length}) / {viscosity:{mu_spec}} = {reynolds_number:,.0f}\n\n"
        
        f"**Step 3:** Determine the flow regime.\n"
        f"For flow in a **{geometry}**, we compare the calculated Re to the standard critical values.\n"
        f"{regime_explanation}\n\n"

        f"**Answer:**\n"
        f"a) The Reynolds number is approximately **{reynolds_number:,.0f}**.\n"
        f"b) The flow regime is **{regime}**."
    )

    return question, solution


# Template 5 (Advanced)
def template_power_law_fluid_shear():
    """
    Shear Stress Calculation for a Power-Law Fluid

    Scenario:
        This template extends Newton's Law of Viscosity to non-Newtonian fluids
        that follow the Ostwald-de Waele (Power Law) model. Unlike Newtonian fluids,
        their viscosity is dependent on the shear rate. This model is crucial for
        describing common substances like paint, ketchup, and suspensions.

        The relevant equations are:
            - Shear Stress: τ_yx = K * (dvx/dy)^n
            - Apparent Viscosity: η = K * |dvx/dy|^(n-1)

    Screen pass 1 (2026-09-23):
        Two judges found every fluid called "non-Newtonian" in the question
        and Step 1 claiming unconditionally that the apparent viscosity
        changes with shear rate, although three fluids in the table (water,
        glycerol, air; 13% of draws) have n = 1. The wording is now
        conditional on n: the question calls the fluid Newtonian when n = 1
        and notes that the power law reduces to Newton's law there, and
        Step 1 says whether the apparent viscosity decreases, increases or
        stays constant with shear rate. No number changes.

    Layer 2 fix (2026-09-25):
        Two experts found the table's blood-plasma entry (K = 0.012 Pa.s^n,
        n = 0.95) about eight times too viscous: on the hand-check instance
        it gives eta = 0.0095 Pa.s at 114.84 1/s, where plasma is essentially
        Newtonian at about 1.2 mPa.s (K looked like a whole-blood value).
        Confirmed. The row is now (0.0012, 1.0), Newtonian at the 1.10-1.30
        mPa.s of normal plasma at 37 degC, sourced in constants.py. Both
        experts also found the answers printed at fixed decimals (tau .3f,
        eta .4f), so the air entry (K = 1.8e-5) always printed
        eta = 0.0000 Pa.s and tau with 0-1 digits (measured: all 40 air
        draws in 1,000); and raw floats printed ("(114.84)^(-0.0500000000
        00000044)", "0.0078000000000000005 m"), on 37% of draws. tau and eta
        are now printed with a 4-decimal-mantissa scientific display, Y is
        bound through its exact 4-dp metre display, and the exponent n - 1
        is formed in decimal at n's own precision. Every displayed-and-
        consumed operand (Y, dvx/dy) is bound through its display, so each
        printed line follows from the printed operands; a draw whose dvx/dy,
        tau or eta sits on a half-way tie at its display is redrawn (D-016).
        Every gold answer changes in format; the question text changes only
        on the plasma draws (about 1 in 22).

    Returns:
        tuple: A tuple containing:
            - str: A question asking for shear stress and apparent viscosity.
            - str: A step-by-step solution.
    """
    # 1. Parameterize the inputs with random values
    fluid_name, (K, n) = random.choice(list(POWER_LAW_FLUIDS.items()))
    # n - 1 is exact at n's own decimals; formed in decimal so that the
    # exponent prints as -0.05, not -0.050000000000000044 (Layer 2).
    n_dp = _decimals(n)
    n_minus_1 = _hu(n - 1, n_dp)

    for _attempt in range(200):
        # Velocity of the moving plate in m/s
        V = round(random.uniform(0.1, 2.0), 2)

        # Distance between plates in cm, converted to meters for calculation.
        # A 2-dp centimetre value is exactly 4 dp in metres; the float is not,
        # so Y is bound through that display (D-016 part 2).
        Y_cm = round(random.uniform(0.5, 5.0), 2)
        Y_m = _as_printed(Y_cm / 100, '.4f')

        # 2. Perform the core calculation. dvx/dy is printed at 2 dp and then
        # consumed, so it is bound through that display; tau and eta through
        # their 4-dp-mantissa displays. A 2-dp over 4-dp quotient, and (for
        # n = 1) a product of K by that quotient, can sit exactly on a
        # half-way tie at their displays; such a draw is redrawn (D-016).
        if _is_display_tie(V / Y_m, 2):
            continue
        velocity_gradient = _as_printed(abs(V / Y_m), '.2f')
        if _is_sci_display_tie(K * (velocity_gradient ** n), 4):
            continue
        shear_stress = _as_printed(K * (velocity_gradient ** n), '.4e')
        if _is_sci_display_tie(K * (velocity_gradient ** (n - 1)), 4):
            continue
        apparent_viscosity = _as_printed(K * (velocity_gradient ** (n - 1)), '.4e')
        break
    else:
        raise RuntimeError(
            "power_law_fluid_shear: no display-stable sample in 200 draws")

    # Determine the fluid behavior for the explanation. The wording is
    # conditional on n (screen pass 1): three fluids in the table are Newtonian
    # (n = 1), and for them the apparent viscosity does NOT change with shear.
    if n < 1:
        behavior = f"shear-thinning (pseudoplastic), because its power-law index n ({n}) is less than 1."
        consequence = "This means its apparent viscosity decreases as the rate of shear increases."
        fluid_class = "non-Newtonian"
        model_note = ""
    elif n > 1:
        behavior = f"shear-thickening (dilatant), because its power-law index n ({n}) is greater than 1."
        consequence = "This means its apparent viscosity increases as the rate of shear increases."
        fluid_class = "non-Newtonian"
        model_note = ""
    else:
        behavior = "Newtonian, because its power-law index n is exactly 1."
        consequence = ("This means its apparent viscosity does not change with the rate of shear: "
                       "it equals the consistency index K, which is then simply the dynamic viscosity.")
        fluid_class = "Newtonian"
        model_note = " (for n = 1 the power law reduces to Newton's law of viscosity)"

    # 3. Generate the question and solution strings
    question = (
        f"A {fluid_class} fluid, {fluid_name.lower()}, is placed between two parallel plates "
        f"separated by {Y_cm} cm. The top plate moves at a constant velocity of {V} m/s, "
        f"creating a linear velocity profile in the fluid.\n\n"
        f"The fluid follows the power-law model{model_note} with the following parameters:\n"
        f"- Consistency Index (K) = {K} Pa·s^n\n"
        f"- Power-Law Index (n) = {n}\n\n"
        f"Calculate:\n"
        f"a) The shear stress (τ_yx) on the fluid.\n"
        f"b) The apparent viscosity (η) at this shear rate."
    )

    solution = (
        f"**Given Information:**\n"
        f"- Fluid: {fluid_name}\n"
        f"- Consistency Index (K): {K} Pa·s^n\n"
        f"- Power-Law Index (n): {n}\n"
        f"- Plate Velocity (V): {V} m/s\n"
        f"- Plate Separation (Y): {Y_cm} cm = {Y_m:.4f} m\n\n"

        f"**Step 1:** Characterize the Fluid Behavior.\n"
        f"The fluid is {behavior} {consequence}\n\n"

        f"**Step 2:** Calculate the Velocity Gradient (Shear Rate).\n"
        f"Assuming a linear velocity profile:\n"
        f"Shear Rate (dvx/dy) = V / Y = {V} m/s / {Y_m:.4f} m = {velocity_gradient:.2f} 1/s\n\n"

        f"**Step 3:** Calculate the Shear Stress (τ_yx).\n"
        f"Using the power-law formula: τ_yx = K * (dvx/dy)^n\n"
        f"τ_yx = {K} * ({velocity_gradient:.2f})^{n}\n"
        f"τ_yx = {shear_stress:.4e} Pa\n\n"

        f"**Step 4:** Calculate the Apparent Viscosity (η).\n"
        f"Apparent viscosity is the effective viscosity at a specific shear rate.\n"
        f"η = K * |dvx/dy|^(n-1)\n"
        f"η = {K} * ({velocity_gradient:.2f})^({n_minus_1:.{n_dp}f})\n"
        f"η = {apparent_viscosity:.4e} Pa·s\n\n"

        f"**Answer:**\n"
        f"a) The shear stress on the fluid is **{shear_stress:.4e} Pa**.\n"
        f"b) The apparent viscosity at this shear rate is **{apparent_viscosity:.4e} Pa·s**."
    )

    return question, solution


def main():
    """
    Generate numerous instances of each viscosity and mechanisms of momentum transport template 
    with different random seeds and write the results to a JSONL file.
    """
    import json
    import os

    # Define the output path (Modify this path according to where you are running the code from)
    output_file = "testset/chemical_engineering/transport_phenomena/viscosity_and_momentum_transport.jsonl"

    # Create the directory if it doesn't exist
    os.makedirs(os.path.dirname(output_file), exist_ok=True)

    # List of template functions with their ID and level
    templates = [
        (template_newtons_law_shear_stress, "newtons_law_shear_stress", "Easy"),
        (template_kinematic_viscosity, "kinematic_viscosity", "Easy"),
        (template_gas_viscosity_kinetic_theory, "gas_viscosity_kinetic_theory", "Intermediate"),
        (template_reynolds_number_flow_regime, "reynolds_number_flow_regime", "Intermediate"),
        (template_power_law_fluid_shear, "power_law_fluid_shear", "Advanced"),
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
                "domain": "transport_phenomena",
                "area": "viscosity_and_momentum_transport",
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
