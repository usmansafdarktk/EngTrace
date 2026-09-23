import random
import math
from decimal import Decimal, ROUND_HALF_UP

from data.templates.branches.mechanical_engineering.constants import SHEAR_MODULUS_VALUES


# Verifier tolerances (D1.3), keyed by template id. Relative.
TEMPLATE_TOLERANCES = {
    # The diameter is quoted to 2 dp in mm over a 19-67 mm range, so half a
    # unit in the last printed place is at most 2.6e-4 relative. 1e-3 sits
    # above that display floor and far below the defect this item carried,
    # where identical printed operands ("d = 2 * 0.0186") produced different
    # printed answers across seeds.
    "template_shaft_design_power": 1e-3,
    # The angle of twist is quoted to five significant figures in radians, so
    # half a unit in the last printed place is at most 5e-5 relative anywhere
    # in the 0.006-0.79 rad range. 2e-3 clears that and the ~5e-7 the stated
    # six-figure J values contribute, with two orders of magnitude to spare.
    "template_composite_shafts_series": 2e-3,
    # Measured over 2,000 seeds: agreement floor 1.40e-3, detection floor 1%.
    # The display floor really is under 2e-5, but display is not what binds
    # here - the trace legitimately consumes a length ratio it states to 3 dp,
    # and the oracle solves the system from the exact lengths instead, so the
    # two differ by that stated rounding. The earlier 1e-4 declaration was 14x
    # below the agreement floor and was argued from display alone.
    "template_statically_indeterminate_shaft": 5e-3,
}


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
def template_shear_stress_torsion():
    """
    Torsion: Shear Stress and Polar Moment of Inertia

    Scenario:
        This template generates a foundational problem testing the ability to calculate
        the polar moment of inertia (J) for a solid or hollow circular shaft.
        It then uses this value to determine the maximum shearing stress (tau_max).

    Core Equations:
        tau_max = (T * c) / J
        J_solid = (pi / 2) * c^4
        J_hollow = (pi / 2) * (c_outer^4 - c_inner^4)

    Trace integrity (Layer 0, 2026-09-23):
        J is displayed to four significant figures and then consumed, so it
        is bound through that display before the stress is computed (D-016
        part 2). The stress is rounded ONCE, at the precision the answer is
        quoted to (3 dp in MPa), and the pascal value is displayed with
        exactly the digits of that answer, so the Pa-to-MPa conversion in
        Step 4 is a decimal shift and not a second rounding: Step 4 used to
        print a 4-figure "5.635e+08 Pa / 1e6" against a 6-figure "563.453
        MPa", which no reader can reproduce (D-016 part 1).

    Returns:
        tuple: A tuple containing:
            - str: A question asking for the maximum shearing stress in a shaft.
            - str: A step-by-step solution to the problem.
    """
    # 1. Parameterize the inputs with random values
    torque = round(random.uniform(100.0, 5000.0), 2)  # Torque in N.m
    shaft_type = random.choice(['solid', 'hollow'])
    
    # Ensure diameters are distinct and practical
    d_outer = random.randint(30, 150)  # Outer diameter in mm
    
    # For hollow shafts, ensure realistic wall thickness
    if shaft_type == 'hollow':
        # Inner diameter is 50-80% of outer diameter (realistic proportions)
        ratio = random.uniform(0.5, 0.8)
        d_inner = round(d_outer * ratio)
    else:
        d_inner = 0

    # Standardize precision for all calculations and outputs
    precision = 3

    # 2. Perform the core calculations for the solution
    
    # Step A: Convert diameters to radii in meters
    c_outer = d_outer / 2000.0  # Convert mm to m and get radius
    if shaft_type == 'hollow':
        c_inner = d_inner / 2000.0
    
    # Step B: Calculate the polar moment of inertia (J). J is displayed to
    # four significant figures and then consumed, so the chain uses the
    # displayed value (D-016 part 2).
    if shaft_type == 'solid':
        polar_moment_J = _as_printed((math.pi / 2) * (c_outer ** 4), '.3e')
        j_calculation_str = f"J = (pi / 2) * c^4 = (pi / 2) * ({c_outer})^4 = {polar_moment_J:.3e} m^4"
    else: # shaft_type == 'hollow'
        polar_moment_J = _as_printed((math.pi / 2) * (c_outer ** 4 - c_inner ** 4), '.3e')
        j_calculation_str = (
            f"J = (pi / 2) * (c_outer^4 - c_inner^4)\n"
            f"   J = (pi / 2) * (({c_outer})^4 - ({c_inner})^4) = {polar_moment_J:.3e} m^4"
        )

    # Step C: Apply the torsion formula to find the stress in Pascals
    # Note: 'c' in the formula refers to the outermost radius, c_outer.
    tau_max_pascals_exact = (torque * c_outer) / polar_moment_J

    # Step D: Convert the final answer to megapascals (MPa). ONE rounding, at
    # the precision the answer is quoted to; the pascal value is then displayed
    # with exactly the digits of that answer (563.453 MPa is 5.63453e+08 Pa),
    # so the decimal conversion in Step 4 is exact (D-016 part 1).
    tau_max_mpa = _as_printed(tau_max_pascals_exact / 1e6, f'.{precision}f')
    tau_pa_spec = '.{}e'.format(
        len(format(tau_max_mpa, f'.{precision}f').replace('.', '').lstrip('0')) - 1)
    tau_max_pascals = _as_printed(tau_max_mpa * 1e6, tau_pa_spec)

    # 3. Generate the question and solution strings
    
    # Construct the part of the question describing the geometry
    if shaft_type == 'solid':
        geometry_desc = f"solid circular shaft with an outer diameter of {d_outer} mm"
    else: # shaft_type == 'hollow'
        geometry_desc = (
            f"hollow circular shaft with an outer diameter of {d_outer} mm "
            f"and an inner diameter of {d_inner} mm"
        )

    question = (
        f"A {geometry_desc} is subjected to a torque of {torque} N.m. "
        f"Determine the maximum shearing stress in the shaft."
    )

    solution = (
        f"**Given:**\n"
        f"Torque (T) = {torque} N.m\n"
        f"Shaft Type = {shaft_type.capitalize()}\n"
        f"Outer Diameter (d_outer) = {d_outer} mm\n"
    )
    if shaft_type == 'hollow':
        solution += f"  - Inner Diameter (d_inner) = {d_inner} mm\n\n"
    else:
        solution += "\n"

    solution += (
        f"**Step 1:** Convert diameters to radii and express in meters.\n"
        f"The maximum stress occurs at the outer surface, so we use the outer radius (c) for the stress calculation.\n"
        f"Outer radius (c) = d_outer / 2 = {d_outer} / 2 = {d_outer/2.0} mm = {c_outer} m\n"
    )
    if shaft_type == 'hollow':
        solution += f"  - Inner radius (c_inner) = d_inner / 2 = {d_inner} / 2 = {d_inner/2.0} mm = {c_inner} m\n"
    solution += "\n"
    
    solution += (
        f"**Step 2:** Calculate the polar moment of inertia (J) for the shaft's cross-section.\n"
        f"For a {shaft_type} shaft:\n"
        f"  {j_calculation_str}\n\n"
        
        f"**Step 3:** Apply the torsion formula to find the maximum shearing stress (tau_max).\n"
        f"The formula is: tau_max = (T * c) / J\n"
        f"\n"
        f"T = {torque} N.m\n"
        f"c = {c_outer} m\n"
        f"J = {polar_moment_J:.3e} m^4\n"
        f"tau_max = ({torque} * {c_outer}) / {polar_moment_J:.3e}\n"
        f"tau_max = {tau_max_pascals:{tau_pa_spec}} Pa\n\n"

        f"**Step 4:** Convert the stress from Pascals (Pa) to Megapascals (MPa).\n"
        f"1 MPa = 1,000,000 Pa\n"
        f"tau_max = {tau_max_pascals:{tau_pa_spec}} Pa / 1e6 = {tau_max_mpa} MPa\n\n"

        f"**Answer:**\n"
        f"The maximum shearing stress in the shaft is {tau_max_mpa} MPa."
    )

    return question, solution


# Template 2 (Easy)
def template_angle_of_twist():
    """
    Torsion: Angle of Twist Calculation

    Scenario:
        This template assesses the ability to calculate the total angle of twist (phi)
        for a uniform solid circular shaft. It requires a correct understanding of the
        relationship between torque, length, material properties (Shear Modulus), and
        the shaft's geometry (Polar Moment of Inertia).

    Core Equations:
        phi = (T * L) / (J * G)
        J_solid = (pi / 2) * c^4

    Trace integrity (Layer 0, 2026-09-23):
        Every quantity that is displayed and then consumed is bound through
        its own display (D-016 part 2): J at four significant figures, G at
        three, and the angle in radians at 4 dp, so the degrees in Step 5 are
        computed from the 4-dp radian value the reader sees. Step 5 used to
        print "0.0579 * (180 / pi) = 3.3181" where the printed operand gives
        3.3174 - the degrees came from an unrounded radian value. A draw whose
        exact radian value sits on a half-way 4-dp tie is redrawn, because no
        rounding closes such a line for every reader (D-016 part 3).

    Returns:
        tuple: A tuple containing:
            - str: A question asking for the angle of twist in a shaft.
            - str: A step-by-step solution to the problem.
    """
    # Standardize precision for final outputs
    precision = 4

    for _attempt in range(200):
        # 1. Parameterize the inputs with random values
        torque = round(random.uniform(500.0, 6000.0), 1)   # Torque in N.m
        length = round(random.uniform(0.5, 4.0), 2)       # Length in m
        diameter = random.randint(25, 100)                # Diameter in mm

        # Randomly select a material and its properties
        material_name, shear_modulus_gpa = random.choice(list(SHEAR_MODULUS_VALUES.items()))

        # 2. Perform the core calculations for the solution

        # Step A: Convert diameter (mm) to radius (m)
        c_radius_m = diameter / 2000.0

        # Step B: Calculate the polar moment of inertia (J). Displayed to four
        # significant figures and then consumed, so bound through that display
        # (D-016 part 2).
        polar_moment_J = _as_printed((math.pi / 2) * (c_radius_m ** 4), '.4e')

        # Step C: Ensure all units are consistent (convert G from GPa to Pa).
        # Displayed to three significant figures, which is exact for every
        # table value; bound so the stored float IS the displayed one.
        shear_modulus_pa = _as_printed(shear_modulus_gpa * 1e9, '.2e')

        # Step D: Apply the angle of twist formula to find the angle in radians.
        # The radian value is displayed at 4 dp and then consumed by Step 5, so
        # the chain continues from the displayed value. A half-way tie at that
        # display has no defensible rounding; the draw is rejected (D-016).
        angle_rad_exact = (torque * length) / (polar_moment_J * shear_modulus_pa)
        if _is_display_tie(angle_rad_exact, precision):
            continue
        angle_rad = _as_printed(angle_rad_exact, f'.{precision}f')

        # Step E: Convert the result from radians to degrees
        angle_deg = math.degrees(angle_rad)
        break
    else:
        raise RuntimeError("angle_of_twist: no closing sample in 200 draws")

    # 3. Generate the question and solution strings
    
    question = (
        f"A solid {material_name.lower()} shaft with a diameter of {diameter} mm and a length of {length} m "
        f"is subjected to a torque of {torque} N.m. "
        f"Given that the shear modulus (G) for {material_name.lower()} is {shear_modulus_gpa} GPa, "
        f"calculate the total angle of twist. Provide the answer in both radians and degrees."
    )

    solution = (
        f"**Given:**\n"
        f"Torque (T) = {torque} N.m\n"
        f"Length (L) = {length} m\n"
        f"Diameter (d) = {diameter} mm\n"
        f"Material = {material_name}\n"
        f"Shear Modulus (G) = {shear_modulus_gpa} GPa\n\n"

        f"**Step 1:** Convert the diameter to a radius in meters.\n"
        f"Radius (c) = d / 2 = {diameter} / 2 = {diameter/2.0} mm = {c_radius_m} m\n\n"

        f"**Step 2:** Calculate the polar moment of inertia (J) for the solid circular shaft.\n"
        f"Formula: J = (pi / 2) * c^4\n"
        f"J = (pi / 2) * ({c_radius_m})^4 = {polar_moment_J:.4e} m^4\n\n"
        
        f"**Step 3:** Ensure consistent units for the angle of twist calculation.\n"
        f"The shear modulus must be in Pascals (Pa) to be consistent with N and m.\n"
        f"G = {shear_modulus_gpa} GPa = {shear_modulus_pa:.2e} Pa\n\n"

        f"**Step 4:** Apply the angle of twist formula to find the angle in radians.\n"
        f"Formula: phi = (T * L) / (J * G)\n"
        f"phi = ({torque} * {length}) / ({polar_moment_J:.4e} * {shear_modulus_pa:.2e})\n"
        f"phi = {angle_rad} radians\n\n"

        f"**Step 5:** Convert the angle from radians to degrees.\n"
        f"Angle in degrees = Angle in radians * (180 / pi)\n"
        f"Angle = {angle_rad} * (180 / pi) = {round(angle_deg, precision)} degrees\n\n"

        f"**Answer:**\n"
        f"The total angle of twist is {angle_rad} radians, which is equivalent to "
        f"{round(angle_deg, precision)} degrees."
    )

    return question, solution


# Template 3 (Intermediate)
def template_shaft_design_power():
    """
    Torsion: Power Transmission and Shaft Design

    Scenario:
        This template generates a classic design problem. It tests the ability to
        first determine the torque on a shaft from its power and rotational speed,
        and then use that torque to calculate the minimum required shaft diameter
        based on an allowable shearing stress.

    Core Equations:
        P = 2 * pi * f * T
        tau_allow = (T * c) / J
        For a solid shaft, this simplifies to: c = ( (2 * T) / (pi * tau) )^(1/3)

    Trace integrity (Phase 1, P2 as amended):
        Every quantity that is displayed and then consumed is bound through its
        own display, so the operand a solver reads IS the operand the chain
        uses. Step 4 previously printed "d = 2 * <c to 4 dp> = <d from the
        UNROUNDED c>", which made 30% of its probed lines non-closing and, more
        damningly, let two seeds print the identical operand "d = 2 * 0.0186"
        and different answers (0.0371 and 0.0373). The radius is now carried at
        5 dp, so doubling it is exact and the metre-to-millimetre conversion is
        exact as well: 5 dp in metres IS 2 dp in millimetres.

    Returns:
        tuple: A tuple containing:
            - str: A question asking for the minimum shaft diameter.
            - str: A step-by-step solution to the design problem.
    """
    precision = 2          # mm, the answer's precision
    # The radius is displayed at 6 dp, not 5. At 5 dp, doubling it put the
    # diameter on a 0.02 mm grid, so only EVEN hundredths of a millimetre were
    # reachable and the item's distinct-answer count fell 29% at 5,000 seeds -
    # a P6 regression that a 1,000-seed measurement cannot see, because the
    # count saturates. At 6 dp the diameter steps by 0.002 mm and every
    # hundredth is reachable.
    C_DP = 6

    for _attempt in range(200):
        power_kw = round(random.uniform(10.0, 500.0), 1)
        allowable_stress_mpa = random.randint(40, 120)

        use_rpm = random.choice([True, False])
        if use_rpm:
            frequency_rpm = random.randint(500, 3000)
            # The frequency is STATED to 2 dp, so 2 dp is what the chain uses.
            frequency_hz = _hu(frequency_rpm / 60.0, 2)
            frequency_str = f"{frequency_rpm} RPM"
        else:
            frequency_hz = _hu(random.randint(10, 100), 2)
            frequency_rpm = frequency_hz * 60
            frequency_str = f"{frequency_hz:.0f} Hz"

        power_w = _hu(power_kw * 1000, 0)
        # An integer number of MPa scaled by 1e6 is exact and ".1e" displays it
        # losslessly, so no display binding is needed here.
        allowable_stress_pa = allowable_stress_mpa * 1e6

        # The torque and c^3 are displayed and then consumed, so their display
        # precision has to be fine enough not to move the answer. At 2 dp and
        # four significant figures respectively they moved the last digit of
        # the 2-dp diameter on 13.6% of draws - c ~ T^(1/3) and c ~ (c^3)^(1/3)
        # so both feed straight through. 4 dp and six significant figures put
        # the combined effect three orders below one display step.
        torque = _hu(power_w / (2 * math.pi * frequency_hz), 4)
        c_cubed = _as_printed((2 * torque) / (math.pi * allowable_stress_pa),
                              ".5e")
        c_exact = c_cubed ** (1 / 3)
        if _is_display_tie(c_exact, C_DP):
            continue
        c_radius_m = _hu(c_exact, C_DP)
        # Doubling a 6-dp value is exact at 6 dp. The millimetre conversion is
        # then a genuine rounding (6 dp in metres is 3 dp in mm), so it carries
        # its own tie guard.
        d_diameter_m = _hu(c_radius_m * 2, C_DP)
        if _is_display_tie(d_diameter_m * 1000, precision):
            continue
        d_diameter_mm = _hu(d_diameter_m * 1000, precision)
        # P1 in one line: the answer reachable from the stated GIVENS must be
        # the answer reachable through the stated INTERMEDIATES. The trace
        # legitimately carries a 2-dp torque and a 4-significant-figure c^3
        # (both printed before use, so P3 holds), but on a value sitting near
        # a 2-dp boundary those roundings can move the last digit: seed 11230
        # gave 9.13 mm through the intermediates and 9.12 mm from the givens
        # (Phase 1 review A, finding F-2). Resample rather than tolerate it -
        # widening the verifier tolerance to cover this would have raised the
        # smallest detectable reasoning error from 0.2% to 0.5%.
        _d_from_givens = 2 * ((2 * (power_w / (2 * math.pi * frequency_hz)))
                              / (math.pi * allowable_stress_pa)) ** (1 / 3) * 1000
        if _hu(_d_from_givens, precision) != d_diameter_mm:
            continue
        break
    else:
        raise RuntimeError("shaft_design_power: no closing sample in 200 draws")

    # --- invariants (T7) ---------------------------------------------------
    assert torque > 0.0, f"non-positive torque: {torque}"
    assert 0.003 <= c_radius_m <= 0.10, f"shaft radius implausible: {c_radius_m}"
    # The designed shaft must actually meet the stress limit it was sized to:
    # tau = 2*T/(pi*c^3) should return the allowable stress. Allow 1% for the
    # display roundings the trace deliberately carries.
    assert abs(2 * torque / (math.pi * c_radius_m ** 3)
               - allowable_stress_pa) <= 0.01 * allowable_stress_pa, (
        f"designed shaft does not meet the allowable stress: c={c_radius_m}")

    question = (
        f"A motor is required to transmit {power_kw} kW of power at a rotational speed of {frequency_str}. "
        f"If the solid circular shaft is to be made from a material with an allowable shearing stress of {allowable_stress_mpa} MPa, "
        f"determine the minimum required diameter for the shaft."
    )

    solution = (
        f"**Given:**\n"
        f"Power (P) = {power_kw} kW\n"
        f"Rotational Speed = {frequency_str}\n"
        f"Allowable Shearing Stress (tau_allow) = {allowable_stress_mpa} MPa\n\n"

        f"**Step 1:** Convert the given values to base SI units (Watts, Pascals, Hertz).\n"
        f"Power (P) = {power_kw} kW = {power_w} W\n"
        f"Allowable Stress (tau_allow) = {allowable_stress_mpa} MPa = {allowable_stress_pa:.1e} Pa\n"
    )

    if use_rpm:
        solution += (
            f"  - Frequency (f) = {frequency_rpm} RPM = {frequency_rpm} / 60 = {frequency_hz:.2f} Hz\n\n"
        )
    else:
        solution += f"  - Frequency (f) = {frequency_str}\n\n"

    solution += (
        f"**Step 2:** Calculate the torque (T) exerted on the shaft.\n"
        f"The relationship between power, torque, and frequency is P = 2 * pi * f * T.\n"
        f"Rearranging for torque: T = P / (2 * pi * f)\n"
        f"T = {power_w} / (2 * pi * {frequency_hz:.2f}) = {torque:.4f} N.m\n\n"

        f"**Step 3:** Determine the required shaft radius (c) using the torsion formula.\n"
        f"The formula for maximum stress in a solid shaft is tau = (T * c) / J, where J = (pi/2) * c^4.\n"
        f"This simplifies to tau = 2 * T / (pi * c^3).\n"
        f"Rearranging to solve for the radius: c^3 = (2 * T) / (pi * tau_allow)\n"
        f"c^3 = (2 * {torque:.4f}) / (pi * {allowable_stress_pa:.1e}) = {c_cubed:.5e} m^3\n"
        f"c = ({c_cubed:.5e})^(1/3) = {c_radius_m:.{C_DP}f} m\n\n"

        f"**Step 4:** Calculate the minimum diameter from the radius.\n"
        f"Diameter (d) = 2 * c\n"
        f"d = 2 * {c_radius_m:.{C_DP}f} = {d_diameter_m:.{C_DP}f} m\n"
        f"In millimeters, d = {d_diameter_m:.{C_DP}f} * 1000 = {d_diameter_mm:.{precision}f} mm\n\n"

        f"**Answer:**\n"
        f"The minimum required diameter for the solid circular shaft is {d_diameter_mm:.{precision}f} mm."
    )

    return question, solution


# Template 4 (Intermediate)
def template_composite_shafts_series():
    """
    Torsion: Composite Shafts in Series

    Scenario:
        This problem involves a shaft made of two different segments joined end-to-end.
        It tests the understanding that the total angle of twist at the free end is
        the algebraic sum of the angles of twist of each individual segment.

    Core Equations:
        phi_total = sum( (T_i * L_i) / (J_i * G_i) )
        For this case: phi_total = phi_1 + phi_2
        J_solid = (pi / 2) * c^4

    Returns:
        tuple: A tuple containing:
            - str: A question about a composite shaft in series.
            - str: A step-by-step solution.
    """
    # 1. Parameterize the inputs with physically realistic values
    
    # Reduce max torque to prevent plastic deformation/failure
    torque = round(random.uniform(200.0, 3000.0), 1)

    # Filter for strong materials (Metals, G > 20 GPa) to ensure elastic behavior
    STRONG_MATERIALS = {k: v for k, v in SHEAR_MODULUS_VALUES.items() if v > 20.0}
    
    assert STRONG_MATERIALS, "SHEAR_MODULUS_VALUES has no material above 20 GPa"

    # Properties for Segment 1 (AB)
    l1 = round(random.uniform(0.5, 2.5), 2)
    d1 = random.randint(40, 120)
    mat1_name, g1_gpa = random.choice(list(STRONG_MATERIALS.items()))

    # Properties for Segment 2 (BC)
    l2 = round(random.uniform(0.5, 2.5), 2)
    d2 = random.randint(30, d1) # Ensure d2 is not larger than d1
    mat2_name, g2_gpa = random.choice(list(STRONG_MATERIALS.items()))
    
    # The segment angles are displayed at PHI_DP and then summed, so the
    # sum is taken over the DISPLAYED values and is exact at PHI_DP. Step 4
    # previously printed the two segments at 5 dp and their total at 4 dp,
    # rounded from the unrounded sum - a double rounding whose printed line
    # disagreed with its own operands verbatim on 89.9% of instances and
    # failed closure on 4.85% (phase0_baseline.md claim 6).
    #  Calculations for Segment 1 (AB)
    c1_m = d1 / 2000.0                       # exact: an integer over 2000
    g1_pa = _as_printed(g1_gpa * 1e9, ".2e")
    # J is displayed and then consumed, so it is bound through its display.
    # At ".3e" that display carried only four significant figures, and since
    # phi ~ 1/J it put ~5e-4 of relative error into every segment angle -
    # enough to move the last printed digit of the smallest totals. ".5e"
    # costs nothing and drops it to ~5e-7.
    j1 = _as_printed((math.pi / 2) * (c1_m ** 4), ".5e")
    phi1_exact = (torque * l1) / (j1 * g1_pa)

    #  Calculations for Segment 2 (BC)
    c2_m = d2 / 2000.0
    g2_pa = _as_printed(g2_gpa * 1e9, ".2e")
    j2 = _as_printed((math.pi / 2) * (c2_m ** 4), ".5e")
    phi2_exact = (torque * l2) / (j2 * g2_pa)

    #  Total Angle of Twist
    #
    # The total twist spans 0.006-0.8 rad, so a fixed 5 dp gave the smallest
    # answers three significant figures and made one instance in a hundred
    # unreachable from the question at that precision. Quote five significant
    # figures instead, and give the segments the SAME precision so their
    # printed sum is exact.
    PHI_DP = max(5, 4 - math.floor(math.log10(phi1_exact + phi2_exact)))
    phi1_rad = _hu(phi1_exact, PHI_DP)
    phi2_rad = _hu(phi2_exact, PHI_DP)
    phi_total_rad = _hu(phi1_rad + phi2_rad, PHI_DP)      # exact at PHI_DP
    DEG_DP = max(4, 4 - math.floor(math.log10(math.degrees(phi_total_rad))))
    phi_total_deg = _hu(math.degrees(phi_total_rad), DEG_DP)

    # --- invariants (T7) ---------------------------------------------------
    # The series total is the algebraic sum of the segment twists, exactly.
    assert phi_total_rad == _hu(phi1_rad + phi2_rad, PHI_DP), (
        f"total twist is not the sum of the segments: "
        f"{phi1_rad} + {phi2_rad} != {phi_total_rad}")
    # The segment that twists more is the one with the larger L/(d^4*G) - but
    # only assert it where the two are actually separable. At seed 14556 the
    # segments differ by 5e-5 relative (0.06868014 vs 0.06868376 rad), so they
    # round to the SAME 5-dp value while the exact proxy still orders them, and
    # a strict correspondence fires on a tie rather than on a defect. Found at
    # 20,000 seeds; invisible at 1,000.
    _k1 = l1 / (d1 ** 4 * g1_gpa)
    _k2 = l2 / (d2 ** 4 * g2_gpa)
    if abs(phi1_exact - phi2_exact) > 1e-6 * max(phi1_exact, phi2_exact):
        assert (phi1_exact > phi2_exact) == (_k1 > _k2), (
            "segment twist ordering contradicts L/(d^4*G)")
    assert 0.0 < phi_total_rad < 2 * math.pi, (
        f"total angle of twist implausible: {phi_total_rad} rad")

    # 3. Generate the question and solution strings
    
    question = (
        f"A composite shaft consists of two segments, AB and BC, rigidly connected at B. "
        f"Segment AB is a solid {mat1_name.lower()} shaft of length {l1} m and diameter {d1} mm. "
        f"Segment BC is a solid {mat2_name.lower()} shaft of length {l2} m and diameter {d2} mm. "
        f"The shaft is fixed at end A, and a torque of {torque} N.m is applied at the free end C. "
        f"Calculate the total angle of twist at end C. "
        f"(Use G = {g1_gpa} GPa for {mat1_name.lower()} and G = {g2_gpa} GPa for {mat2_name.lower()})."
    )

    solution = (
        f"**Given:**\n"
        f"Applied Torque (T) = {torque} N.m\n"
        f"Segment AB: L1={l1} m, d1={d1} mm, Material={mat1_name} (G1={g1_gpa} GPa)\n"
        f"Segment BC: L2={l2} m, d2={d2} mm, Material={mat2_name} (G2={g2_gpa} GPa)\n\n"
        
        f"**Step 1:** Analyze the System\n"
        f"The shaft is fixed at A, so the torque T = {torque} N.m is transmitted through both segments AB and BC. "
        f"The total angle of twist at C is the sum of the twist in segment AB and the twist in segment BC.\n"
        f"phi_total = phi_AB + phi_BC\n"
        f"\n\n"

        f"**Step 2:** Calculate Angle of Twist for Segment AB (phi_AB)\n"
        f"Convert units: c1 = {d1}/2 mm = {c1_m} m; G1 = {g1_gpa} GPa = {g1_pa:.2e} Pa\n"
        f"Polar Moment of Inertia (J1) = (pi/2) * c1^4 = (pi/2) * ({c1_m})^4 = {j1:.5e} m^4\n"
        f"Angle of Twist (phi_AB) = (T * L1) / (J1 * G1)\n"
        f"phi_AB = ({torque} * {l1}) / ({j1:.5e} * {g1_pa:.2e}) = {phi1_rad:.{PHI_DP}f} radians\n\n"
        
        f"**Step 3:** Calculate Angle of Twist for Segment BC (phi_BC)\n"
        f"Convert units: c2 = {d2}/2 mm = {c2_m} m; G2 = {g2_gpa} GPa = {g2_pa:.2e} Pa\n"
        f"Polar Moment of Inertia (J2) = (pi/2) * c2^4 = (pi/2) * ({c2_m})^4 = {j2:.5e} m^4\n"
        f"Angle of Twist (phi_BC) = (T * L2) / (J2 * G2)\n"
        f"phi_BC = ({torque} * {l2}) / ({j2:.5e} * {g2_pa:.2e}) = {phi2_rad:.{PHI_DP}f} radians\n\n"
        
        f"**Step 4:** Calculate Total Angle of Twist at End C\n"
        f"phi_total = phi_AB + phi_BC\n"
        f"phi_total = {phi1_rad:.{PHI_DP}f} + {phi2_rad:.{PHI_DP}f} = {phi_total_rad:.{PHI_DP}f} radians\n"
        f"To convert to degrees: Angle_deg = Angle_rad * (180 / pi)\n"
        f"phi_total = {phi_total_rad:.{PHI_DP}f} * (180 / pi) = {phi_total_deg:.{DEG_DP}f} degrees\n\n"
        
        f"**Answer:**\n"
        f"The total angle of twist at the free end C is {phi_total_rad:.{PHI_DP}f} radians, "
        f"or {phi_total_deg:.{DEG_DP}f} degrees."
    )

    return question, solution


# Template 5 (Advanced)
def template_statically_indeterminate_shaft():
    """
    Torsion: Statically Indeterminate Shaft

    Scenario:
        This problem deals with a shaft that is fixed at both ends and has a torque
        applied at an intermediate point. Since there are two unknown reaction torques
        but only one static equilibrium equation, the problem is statically
        indeterminate.

    Core Equations:
        1. Statics: T_A + T_B = T_applied
        2. Compatibility: T_A * L_AC = T_B * L_BC

    Returns:
        tuple: A tuple containing:
            - str: A question asking for the reaction torques.
            - str: A step-by-step solution using statics and compatibility.
    """
    # 1. Parameterize the inputs with physically realistic values
    
    precision = 2          # N.m, the precision the reactions are stated to
    total_length = round(random.uniform(2.0, 5.0), 2)
    # Reduce max torque to 5000 Nm to ensure elastic behavior for given diameters
    applied_torque = round(random.uniform(500.0, 5000.0), -2) 
    diameter = random.randint(50, 150) # mm
    
    # Ensure the torque is applied at a non-trivial location. The draw is
    # repeated if it puts T_A exactly on a display tie - see below.
    for _attempt in range(200):
        pos_L_AC = round(
            random.uniform(0.2 * total_length, 0.8 * total_length), 2)
        # Bound through its own 2-dp display: a float subtraction of two 2-dp
        # values is not exactly 2 dp (3.5 - 1.71 gives 1.7899999999999998).
        pos_L_BC = _hu(total_length - pos_L_AC, 2)
        if pos_L_BC <= 0:
            continue
        _divisor = _hu(1 + _hu(pos_L_AC / pos_L_BC, 3), 3)
        if not _is_display_tie(applied_torque / _divisor, 2):
            break
    else:
        raise RuntimeError(
            "statically_indeterminate_shaft: no closing position in 200 draws")

    # FILTER: Select only strong materials (Metals, G > 20 GPa) to avoid failure
    # This prevents assigning 5000 Nm of torque to a Nylon shaft.
    strong_materials = {k: v for k, v in SHEAR_MODULUS_VALUES.items() if v > 20.0}
    # P5: no silent fallback - an empty table is a bug, not a default.
    assert strong_materials, "SHEAR_MODULUS_VALUES has no material above 20 GPa"
    material_name, shear_modulus_gpa = random.choice(list(strong_materials.items()))
    
    # 2. Perform the core calculations
    #
    # Step 3 prints "T_A = T / <divisor> = <T_A>" with the divisor displayed at
    # 3 dp while T_A was derived from the UNROUNDED ratio, which left 47% of
    # probed lines non-closing (phase0_baseline.md NEW-2; seed 0 printed
    # "3900 / 1.291 = 3021.85" against an exact 3020.91). The ratio is now
    # bound to the 3 dp it is stated at and the chain consumes that value, so
    # T_A is reachable from the printed operands; T_B then follows from
    # statics on the stated T_A, exactly.
    length_ratio = _hu(pos_L_AC / pos_L_BC, 3)
    divisor = _hu(1 + length_ratio, 3)                  # exact at 3 dp
    reaction_torque_A = _hu(applied_torque / divisor, precision)
    reaction_torque_B = _hu(applied_torque - reaction_torque_A, precision)

    # --- invariants (T7) ---------------------------------------------------
    # Statics: the two reactions must carry the applied torque exactly, at the
    # precision they are stated to.
    assert reaction_torque_A + reaction_torque_B == applied_torque, (
        f"statics violated: {reaction_torque_A} + {reaction_torque_B} "
        f"!= {applied_torque}")
    assert 0.0 < reaction_torque_A < applied_torque, (
        f"reaction at A outside (0, T): {reaction_torque_A}")
    # Compatibility: T_A * L_AC = T_B * L_BC, to within the 3-dp rounding of
    # the length ratio that the trace deliberately carries.
    assert abs(reaction_torque_A * pos_L_AC - reaction_torque_B * pos_L_BC) \
        <= 0.002 * applied_torque * total_length, (
        "compatibility violated: T_A*L_AC != T_B*L_BC")

    # 3. Generate the question and solution strings
    
    question = (
        f"A solid {material_name.lower()} shaft of diameter {diameter} mm and length {total_length} m is "
        f"fixed at both ends, A and B. A torque of {int(applied_torque)} N.m is applied at point C, "
        f"located {pos_L_AC} m from end A. The shear modulus for the material is {shear_modulus_gpa} GPa. "
        f"Determine the reaction torques at the fixed supports, T_A and T_B."
    )

    solution = (
        f"**Given:**\n"
        f"Total Length (L) = {total_length} m\n"
        f"Applied Torque (T_applied) = {int(applied_torque)} N.m\n"
        f"Location of Torque from A (L_AC) = {pos_L_AC} m\n"
        f"Diameter (d) = {diameter} mm, Material = {material_name}\n\n"
        
        f"**Analysis:**\n"
        f"This problem is statically indeterminate because there are two unknown reaction torques (T_A and T_B) "
        f"and only one relevant equation from statics. We need an additional equation from the deformation of the shaft.\n"
        f"\n\n"

        f"**Step 1:** Statics Equilibrium Equation\n"
        f"For the shaft to be in rotational equilibrium, the sum of all torques must be zero.\n"
        f"(1) T_A + T_B = T_applied = {int(applied_torque)} N.m\n\n"

        f"**Step 2:** Compatibility Equation\n"
        f"Since the shaft is fixed at both ends, the total angle of twist from A to B must be zero. "
        f"The twist caused by T_A acting on length AC must equal the twist caused by T_B acting on length BC.\n"
        f"phi_AC = phi_CB\n"
        f"(T_A * L_AC) / (J*G) = (T_B * L_BC) / (J*G)\n"
        f"Since J and G are constant, they cancel out:\n"
        f"(2) T_A * L_AC = T_B * L_BC\n\n"

        f"**Step 3:** Solve the System of Two Equations\n"
        f"From equation (2), express T_B in terms of T_A:\n"
        f"T_B = T_A * ({pos_L_AC} / {pos_L_BC:.2f})\n"
        f"Substitute this into equation (1):\n"
        f"T_A + T_A * ({pos_L_AC} / {pos_L_BC:.2f}) = {int(applied_torque)}\n"
        f"T_A * (1 + {length_ratio:.3f}) = {int(applied_torque)}\n"
        f"T_A = {int(applied_torque)} / {divisor:.3f} = {reaction_torque_A:.{precision}f} N.m\n"
        f"Now find T_B:\n"
        f"T_B = {int(applied_torque)} - {reaction_torque_A:.{precision}f} = {reaction_torque_B:.{precision}f} N.m\n\n"

        f"**Answer:**\n"
        f"The reaction torques at the supports are:\n"
        f"T_A = {reaction_torque_A:.{precision}f} N.m\n"
        f"T_B = {reaction_torque_B:.{precision}f} N.m"
    )

    return question, solution


def main():
    """
    Generate numerous instances of each torsion template 
    with different random seeds and write the results to a JSONL file.
    """
    import json
    import os

    # Define the output path (Modify this path according to where you are running the code from)
    output_file = "testset/mechanical_engineering/mechanics_of_materials/torsion.jsonl"

    # Create the directory if it doesn't exist
    os.makedirs(os.path.dirname(output_file), exist_ok=True)

    # List of template functions with their ID and level
    templates = [
        (template_shear_stress_torsion, "shear_stress_torsion", "Easy"),
        (template_angle_of_twist, "angle_of_twist", "Easy"),
        (template_shaft_design_power, "shaft_design_power", "Intermediate"),
        (template_composite_shafts_series, "composite_shafts_series", "Intermediate"),
        (template_statically_indeterminate_shaft, "statically_indeterminate_shaft", "Advanced"),
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
                "branch": "mechanical_engineering",
                "domain": "mechanics_of_materials",
                "area": "torsion",
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
