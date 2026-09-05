import random
import math
from decimal import Decimal, ROUND_HALF_UP


# Verifier tolerances (D1.3), keyed by template id. Relative.
TEMPLATE_TOLERANCES = {
    # The amplitude answer is quoted to four significant figures, so half a
    # unit in the last printed place is at most 5e-4 relative - the worst case
    # is a mantissa just above 1.000, e.g. 1.0744976 printed as "1.074", and
    # the measured worst over 1,000 seeds is 4.63e-4 at seed 220, which IS
    # that bound rather than any residual chain error. 1e-3 is twice that hard
    # display bound, so nothing fires on presentation, and it is still two
    # orders of magnitude below the error the pre-fix chain injected by
    # re-deriving omega from an unrounded shaft speed (up to 8.4%).
    "template_rotating_unbalance": 1e-3,
    # Measured over 2,000 seeds: worst oracle-vs-trace disagreement on a
    # correct instance ("agreement floor") is 1.81e-3, and the smallest
    # injected answer error caught on >=99% of instances ("detection floor")
    # is 1%. The declared value must sit between the two. 1e-3 was BELOW the
    # agreement floor and would have fired on presentation alone.
    "template_vibration_transmissibility": 3e-3,
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


def _is_display_tie(x, places, rel_band=1e-12):
    """Is `x` at, or within a hair of, a half-way tie at `places` dp?

    A tie is the one case where NO rounding convention is defensible. A reader
    doing decimal arithmetic and applying half-up reads 0.02325 m as 23.3 mm; a
    reader using binary floats and `round()` reads it as 23.2 mm; and the
    printed line closes for exactly one of them whichever the template picks.
    Breaking the tie in decimal (P2 as amended) does not remove the ambiguity,
    it only moves it to the other reader - which is why such instances are
    RESAMPLED rather than resolved (D-016).

    Testing for an EXACT tie is not enough. Two independent evaluations of the
    same exact quantity - this template's, in kN and kN/m^2, and a solver's, in
    N and Pa - differ by a few ulps, so a rational tie such as 29.25 mm lands
    as 29.249999999999996 on one side and 29.250000000000004 on the other.
    Neither is exactly a tie, and the two then round in opposite directions:
    that is how eight non-closing instances per 60,000 seeds survived the first
    version of this guard (Phase 1 review A, finding F-2).

    So a narrow BAND is quarantined rather than a point. The band is a few
    thousand ulps wide where the tie lattice is 10^-places apart, so it removes
    nothing that is not genuinely ambiguous.
    """
    scaled = abs(x) * 10.0 ** places
    band = max(scaled * rel_band, 1e-9)
    return abs((scaled - math.floor(scaled)) - 0.5) <= band


# Template 1 (Easy)
def template_system_properties():
    """
    Vibrations: System Properties (Natural Frequency and Damping Ratio)

    Scenario:
        This template generates a foundational problem for analyzing a single-degree-of-freedom
        spring-mass-damper system. Given the system's physical parameters (mass, stiffness,
        and damping coefficient), the objective is to calculate its key dynamic properties:
        the undamped natural frequency, the critical damping coefficient, and the damping ratio.
        Finally, the system is classified based on its damping level.

    Core Equations:
        omega_n = sqrt(k / m)
        c_cr = 2 * sqrt(k * m)
        zeta = c / c_cr

    Returns:
        tuple: A tuple containing:
            - str: A question asking for the system's dynamic properties and classification.
            - str: A step-by-step solution to the problem.
    """
    # 1. Parameterize the inputs with random values
    
    # Mass in kg, ensuring a practical range
    mass = round(random.uniform(2.0, 250.0), 2)
    
    # Stiffness in N/m
    stiffness = round(random.uniform(1000.0, 150000.0), 1)
    
    # To generate diverse and meaningful scenarios (underdamped, critically damped, overdamped),
    # we first determine the critical damping and then set the actual damping based on it.
    
    # We select a random damping ratio first to define the system type
    target_zeta = round(random.uniform(0.1, 2.5), 3) 
    
    # Handle the edge case of m=0 or k=0, though randomization range prevents it.
    if mass <= 0 or stiffness <= 0:
        # Assign default safe values in the unlikely event of non-positive inputs
        mass = 10.0
        stiffness = 20000.0

    critical_damping = 2 * math.sqrt(stiffness * mass)
    
    # Calculate the actual damping coefficient based on the target zeta
    damping_coeff = round(target_zeta * critical_damping, 2)

    # Standardize precision for all calculations and outputs
    precision = 4

    # 2. Perform the core calculations for the solution
    
    # Step A: Calculate the undamped natural frequency (omega_n)
    omega_n = math.sqrt(stiffness / mass)
    
    # Step B: The critical damping coefficient (c_cr) was already calculated
    # We re-calculate here to show the step clearly in the solution.
    c_critical = 2 * math.sqrt(stiffness * mass)
    
    # Step C: Calculate the damping ratio (zeta)
    damping_ratio = damping_coeff / c_critical
    
    # Step D: Classify the system based on the damping ratio
    if abs(damping_ratio - 1.0) < 1e-9: # Use tolerance for floating point comparison
        system_type = "critically damped"
    elif damping_ratio < 1.0:
        system_type = "underdamped"
    else:
        system_type = "overdamped"

    # 3. Generate the question and solution strings
    
    question = (
        f"A spring-mass-damper system has the following properties:\n"
        f"Mass (m) = {mass} kg\n"
        f"Spring Stiffness (k) = {stiffness} N/m\n"
        f"Damping Coefficient (c) = {damping_coeff} N.s/m\n\n"
        f"Determine the following:\n"
        f"  1. The undamped natural frequency (in rad/s).\n"
        f"  2. The critical damping coefficient.\n"
        f"  3. The damping ratio.\n"
        f"  4. Classify the system as underdamped, critically damped, or overdamped."
    )

    solution = (
        f"**Given:**\n"
        f"Mass (m) = {mass} kg\n"
        f"Stiffness (k) = {stiffness} N/m\n"
        f"Damping Coefficient (c) = {damping_coeff} N.s/m\n\n"
        
        f"**Step 1:** Calculate the undamped natural frequency (omega_n).\n"
        f"The formula is: omega_n = sqrt(k / m)\n"
        f"omega_n = sqrt({stiffness} / {mass})\n"
        f"omega_n = {round(omega_n, precision)} rad/s\n\n"

        f"**Step 2:** Calculate the critical damping coefficient (c_cr).\n"
        f"The formula is: c_cr = 2 * sqrt(k * m)\n"
        f"c_cr = 2 * sqrt({stiffness} * {mass})\n"
        f"c_cr = {round(c_critical, precision)} N.s/m\n\n"

        f"**Step 3:** Calculate the damping ratio (zeta).\n"
        f"The formula is: zeta = c / c_cr\n"
        f"zeta = {damping_coeff} / {round(c_critical, precision)}\n"
        f"zeta = {round(damping_ratio, precision)}\n\n"

        f"**Step 4:** Classify the system based on the damping ratio.\n"
        f"The damping ratio is {round(damping_ratio, precision)}.\n"
        f"Since zeta is {'less than 1' if system_type == 'underdamped' else ('equal to 1' if system_type == 'critically damped' else 'greater than 1')}, "
        f"the system is classified as **{system_type}**.\n\n"
        
        f"**Answer:**\n"
        f"Undamped Natural Frequency: {round(omega_n, precision)} rad/s\n"
        f"Critical Damping Coefficient: {round(c_critical, precision)} N.s/m\n"
        f"Damping Ratio: {round(damping_ratio, precision)}\n"
        f"System Type: {system_type.capitalize()}"
    )

    return question, solution


# Template 2 (Intermediate)
def template_rotating_unbalance():
    """
    Vibrations: Response to Rotating Unbalance

    Scenario:
        This template generates a problem involving a machine with a rotating component
        that is out of balance. This common engineering scenario creates a harmonic
        excitation force whose magnitude is dependent on the operating speed. The goal
        is to determine the resulting steady-state vibration amplitude.

    Core Equations:
        omega_n = sqrt(k / m)
        zeta = c / (2 * sqrt(k * m))
        r = omega / omega_n
        F_0 = m_e * e * omega^2
        X = F_0 / sqrt((k - m * omega^2)^2 + (c * omega)^2)

    Trace integrity (Phase 1, P3):
        The sampler still CHOOSES a frequency ratio and a damping ratio to
        place the instance in a useful regime, but neither leaks into the
        solution. The shaft speed is stated as a whole number of RPM and
        omega is then recomputed FORWARD from that stated speed; the damping
        coefficient is stated to 2 dp and the chain consumes the stated
        value; and zeta is derived from the stated c and the computed c_cr
        rather than printed from the sampled pre-image. Every intermediate
        that is displayed and then consumed is rounded to its display
        precision first, with half-way ties resampled (D-016).

    Returns:
        tuple: A tuple containing:
            - str: A question asking for the steady-state amplitude of vibration.
            - str: A step-by-step solution to the problem.
    """
    precision = 5

    for _attempt in range(200):
        # Total mass of the machine in kg (INCLUDING the unbalanced mass -
        # the question now says so explicitly; it previously did not, and a
        # solver who subtracted m_e was marked wrong).
        m_total = round(random.uniform(50.0, 500.0), 1)

        # Eccentric mass in kg (a small fraction of the total mass)
        m_eccentric = round(random.uniform(0.1, m_total * 0.05), 2)

        # Eccentricity in mm
        eccentricity_mm = random.randint(10, 150)

        # Stiffness in N/m
        stiffness = round(random.uniform(5e4, 2e6), 0)

        # Sampler-side choices: a damping ratio and a frequency ratio, picked
        # to spread the instances below, near and above resonance. NEITHER is
        # printed - both are re-derived below from the stated quantities.
        damping_ratio_zeta = round(random.uniform(0.05, 0.6), 3)
        freq_ratio_r = round(random.uniform(0.3, 3.0), 3)

        omega_n_exact = math.sqrt(stiffness / m_total)
        c_critical = _hu(2 * math.sqrt(stiffness * m_total), precision)
        # c is STATED to 2 dp, so 2 dp is what the chain consumes.
        damping_coeff = _hu(damping_ratio_zeta * c_critical, 2)
        zeta = _hu(damping_coeff / c_critical, precision)
        omega_n = _hu(omega_n_exact, precision)

        # The speed is STATED as a whole number of RPM; omega is recomputed
        # forward from it. The sampled freq_ratio_r reaches no further.
        operating_speed_rpm = round(freq_ratio_r * omega_n_exact * 60
                                    / (2 * math.pi), 0)
        omega_exact = operating_speed_rpm * 2 * math.pi / 60
        omega = _hu(omega_exact, precision)

        eccentricity_m = eccentricity_mm / 1000.0

        force_magnitude_F0 = _hu(m_eccentric * eccentricity_m * omega ** 2,
                                 precision)
        term1 = _hu(stiffness - m_total * omega ** 2, precision)
        term2 = _hu(damping_coeff * omega, precision)
        denominator = _hu(math.sqrt(term1 ** 2 + term2 ** 2), precision)
        if denominator <= 0:
            continue
        amplitude_exact = force_magnitude_F0 / denominator          # m

        # A final answer carried to one significant figure ("0.002 mm") is not
        # gradeable, and this item's amplitudes span four orders of magnitude.
        # Quote the answer to FOUR significant figures wherever that needs more
        # than 3 dp, and keep the metre display exactly 3 dp finer so the mm
        # conversion is exact for every reader.
        mm_dp = max(3, 3 - math.floor(math.log10(abs(amplitude_exact * 1000))))
        m_dp = mm_dp + 3

        if (_is_display_tie(omega_exact, precision)
                or _is_display_tie(amplitude_exact, m_dp)):
            continue                       # no defensible gold answer; redraw
        amplitude_m = _hu(amplitude_exact, m_dp)
        amplitude_mm = _hu(amplitude_m * 1000, mm_dp)
        break
    else:
        raise RuntimeError("rotating_unbalance: no closing sample in 200 draws")

    # --- invariants (T7) ---------------------------------------------------
    assert 0.0 < zeta < 1.0, f"system not underdamped: zeta = {zeta}"
    assert amplitude_mm > 0.0, f"non-positive amplitude: {amplitude_mm}"
    # The rotating-unbalance response is bounded by its resonant peak,
    # X_max = (m_e*e/M) / (2*zeta*sqrt(1 - zeta^2)), for every frequency ratio.
    assert amplitude_m <= 1.001 * (m_eccentric * eccentricity_m / m_total) / (
        2 * zeta * math.sqrt(1 - zeta ** 2)), (
        f"amplitude {amplitude_m} exceeds the resonant bound")

    question = (
        f"A machine with a total mass of {m_total} kg (this total includes the "
        f"rotating unbalanced mass) is supported by a spring and damper system. "
        f"The system has an equivalent stiffness of {stiffness:,.0f} N/m and an equivalent damping "
        f"coefficient of {damping_coeff} N.s/m.\n\n"
        f"The machine contains a rotating component that has an unbalance equivalent to a mass of "
        f"{m_eccentric} kg located at an eccentricity of {eccentricity_mm} mm. "
        f"If the machine operates at a speed of {operating_speed_rpm:,.0f} RPM, "
        f"determine the steady-state amplitude of vibration."
    )

    solution = (
        f"**Given:**\n"
        f"Total Mass (m) = {m_total} kg (includes the unbalanced mass)\n"
        f"Stiffness (k) = {stiffness:,.0f} N/m\n"
        f"Damping Coefficient (c) = {damping_coeff} N.s/m\n"
        f"Eccentric Mass (m_e) = {m_eccentric} kg\n"
        f"Eccentricity (e) = {eccentricity_mm} mm = {eccentricity_m} m\n"
        f"Operating Speed = {operating_speed_rpm:,.0f} RPM\n\n"

        f"**Step 1:** Convert the operating speed from RPM to rad/s.\n"
        f"omega = (Speed in RPM) * (2 * pi / 60)\n"
        f"omega = {operating_speed_rpm:,.0f} * (2 * pi / 60) = {omega} rad/s\n\n"

        f"**Step 2:** Calculate the system's natural frequency (omega_n) and damping ratio (zeta).\n"
        f"omega_n = sqrt(k / m) = sqrt({stiffness:,.0f} / {m_total}) = {omega_n} rad/s\n"
        f"c_cr = 2 * sqrt(k * m) = 2 * sqrt({stiffness:,.0f} * {m_total}) = {c_critical} N.s/m\n"
        f"zeta = c / c_cr = {damping_coeff} / {c_critical} = {zeta}\n\n"

        f"**Step 3:** Calculate the magnitude of the unbalanced force (F_0).\n"
        f"The force from a rotating unbalance is given by F_0 = m_e * e * omega^2.\n"
        f"F_0 = {m_eccentric} kg * {eccentricity_m} m * ({omega} rad/s)^2\n"
        f"F_0 = {force_magnitude_F0} N\n\n"

        f"**Step 4:** Calculate the steady-state amplitude of vibration (X).\n"
        f"The formula for amplitude is: X = F_0 / sqrt((k - m * omega^2)^2 + (c * omega)^2)\n"
        f"Numerator = F_0 = {force_magnitude_F0} N\n"
        f"Denominator Part 1: (k - m * omega^2) = ({stiffness:,.0f} - {m_total} * {omega}^2) = {term1}\n"
        f"Denominator Part 2: (c * omega) = ({damping_coeff} * {omega}) = {term2}\n"
        f"Denominator = sqrt(({term1})^2 + ({term2})^2) = {denominator}\n"
        f"X = {force_magnitude_F0} / {denominator} = {amplitude_m:.{m_dp}f} m\n\n"

        f"**Step 5:** Convert the amplitude to millimeters.\n"
        f"Amplitude in mm = {amplitude_m:.{m_dp}f} m * 1000 mm/m = "
        f"{amplitude_mm:.{mm_dp}f} mm\n\n"

        f"**Answer:**\n"
        f"The steady-state amplitude of vibration is "
        f"**{amplitude_mm:.{mm_dp}f} mm**."
    )

    return question, solution


# Template 3 (Intermediate)
def template_vibration_transmissibility():
    """
    Vibrations: Displacement Transmissibility from Base Excitation

    Scenario:
        This template addresses the problem of a system mounted on a vibrating foundation.
        It is a key concept in vibration isolation, where the goal is to minimize the
        motion transmitted from a vibrating source to a sensitive component. The template
        calculates the ratio of the output amplitude to the input (base) amplitude
        and the absolute amplitude of the system.

    Core Equations:
        omega_n = sqrt(k / m)
        zeta = c / (2 * sqrt(k * m))
        r = omega / omega_n
        TR = X / Y = sqrt( (1 + (2*zeta*r)^2) / ( (1 - r^2)^2 + (2*zeta*r)^2 ) )

    Trace integrity (Phase 1, P3):
        The base excitation frequency is STATED to 2 dp, so omega is
        recomputed forward from the stated frequency and the frequency ratio
        is derived as omega/omega_n from the stated quantities. The sampler
        still chooses a frequency ratio to place the instance below, near or
        above resonance, but that choice no longer reaches the solution: the
        solution previously printed the SAMPLED ratio at 3 dp (r = 3.099
        where the stated givens imply 3.09930) and carried it into
        (1 - r^2)^2, which made 39% of answers not the correct rounding of
        what the question implies.

    Returns:
        tuple: A tuple containing:
            - str: A question asking for the transmissibility and absolute amplitude.
            - str: A step-by-step solution to the problem.
    """
    precision = 4

    for _attempt in range(200):
        mass = round(random.uniform(5.0, 150.0), 1)                 # kg
        stiffness = round(random.uniform(2e3, 5e5), 0)              # N/m
        base_amplitude_Y_mm = round(random.uniform(0.1, 8.0), 2)    # mm

        # Sampler-side choice only: spreads the instances across amplification
        # (r ~ 1) and isolation (r > sqrt(2)). Never printed.
        freq_ratio_r = round(random.uniform(0.2, 5.0), 3)
        # zeta IS a stated given, at 3 dp, so the chain may use it directly.
        damping_ratio_zeta = round(random.uniform(0.05, 0.7), 3)

        omega_n_exact = math.sqrt(stiffness / mass)
        # The base frequency is STATED to 2 dp; omega follows from it.
        base_freq_hz = round(freq_ratio_r * omega_n_exact / (2 * math.pi), 2)

        omega_n = _hu(omega_n_exact, precision)
        omega = _hu(base_freq_hz * 2 * math.pi, precision)
        r = _hu(omega / omega_n, precision)

        base_amplitude_Y_m = base_amplitude_Y_mm / 1000.0

        two_zeta_r_sq = (2 * damping_ratio_zeta * r) ** 2
        tr_num = _hu(1 + two_zeta_r_sq, precision)
        tr_den = _hu((1 - r ** 2) ** 2 + two_zeta_r_sq, precision)
        if tr_den <= 0:
            continue
        tr_exact = math.sqrt(tr_num / tr_den)

        # A transmissibility of 0.056 printed to 3 dp is two significant
        # figures - not gradeable. Quote both answers to at least four
        # significant figures.
        tr_dp = max(3, 3 - math.floor(math.log10(abs(tr_exact))))
        if _is_display_tie(tr_exact, tr_dp):
            continue
        transmissibility_ratio = _hu(tr_exact, tr_dp)

        x_exact = transmissibility_ratio * base_amplitude_Y_mm      # mm
        x_dp = max(3, 3 - math.floor(math.log10(abs(x_exact))))
        if (_is_display_tie(x_exact, x_dp)
                or _is_display_tie(omega / omega_n, precision)):
            continue
        amplitude_X_mm = _hu(x_exact, x_dp)
        break
    else:
        raise RuntimeError(
            "vibration_transmissibility: no closing sample in 200 draws")

    # --- invariants (T7) ---------------------------------------------------
    assert 0.0 < damping_ratio_zeta < 1.0, (
        f"system not underdamped: zeta = {damping_ratio_zeta}")
    assert transmissibility_ratio > 0.0, (
        f"non-positive transmissibility: {transmissibility_ratio}")
    # TR = 1 exactly at r = sqrt(2) for every zeta, and the isolation region
    # r > sqrt(2) is precisely where TR < 1. Allow one display step of slack
    # so an instance sitting on the crossover is not a false failure.
    _step = 10.0 ** (-tr_dp)
    if r > math.sqrt(2) + _step:
        assert transmissibility_ratio < 1.0 + _step, (
            f"r = {r} is in the isolation region but TR = "
            f"{transmissibility_ratio}")
    elif r < math.sqrt(2) - _step:
        assert transmissibility_ratio > 1.0 - _step, (
            f"r = {r} is below the crossover but TR = "
            f"{transmissibility_ratio}")

    question = (
        f"A sensitive instrument of mass {mass} kg is supported by an isolation mount. "
        f"The mount has an effective stiffness of {stiffness:,.0f} N/m and provides a damping "
        f"ratio of {damping_ratio_zeta}.\n\n"
        f"The foundation on which the instrument is placed is vibrating harmonically at a frequency of "
        f"{base_freq_hz} Hz with an amplitude of {base_amplitude_Y_mm} mm.\n\n"
        f"Determine:\n"
        f"1. The displacement transmissibility ratio.\n"
        f"2. The absolute amplitude of vibration of the instrument in millimeters."
    )

    solution = (
        f"**Given:**\n"
        f"Mass (m) = {mass} kg\n"
        f"Stiffness (k) = {stiffness:,.0f} N/m\n"
        f"Damping Ratio (zeta) = {damping_ratio_zeta}\n"
        f"Base Vibration Frequency (f) = {base_freq_hz} Hz\n"
        f"Base Vibration Amplitude (Y) = {base_amplitude_Y_mm} mm = {base_amplitude_Y_m} m\n\n"

        f"**Step 1:** Calculate the system's undamped natural frequency (omega_n).\n"
        f"omega_n = sqrt(k / m) = sqrt({stiffness:,.0f} / {mass}) = {omega_n} rad/s\n\n"

        f"**Step 2:** Convert the base vibration frequency to rad/s and find the frequency ratio (r).\n"
        f"Base frequency (omega) = f * 2 * pi = {base_freq_hz} * 2 * pi = {omega} rad/s\n"
        f"Frequency ratio (r) = omega / omega_n = {omega} / {omega_n} = {r}\n\n"

        f"**Step 3:** Calculate the displacement transmissibility ratio (TR).\n"
        f"The formula is: TR = sqrt( (1 + (2*zeta*r)^2) / ( (1 - r^2)^2 + (2*zeta*r)^2 ) )\n"
        f"Let's calculate the terms:\n"
        f"r = {r}\n"
        f"zeta = {damping_ratio_zeta}\n"
        f"Numerator = 1 + (2 * {damping_ratio_zeta} * {r})^2 = {tr_num}\n"
        f"Denominator = (1 - ({r})^2)^2 + (2 * {damping_ratio_zeta} * {r})^2 = {tr_den}\n"
        f"TR = sqrt({tr_num} / {tr_den}) = {transmissibility_ratio:.{tr_dp}f}\n\n"

        f"**Step 4:** Calculate the absolute amplitude of the instrument (X).\n"
        f"The relationship is X = TR * Y.\n"
        f"X = {transmissibility_ratio:.{tr_dp}f} * {base_amplitude_Y_mm} mm = "
        f"{amplitude_X_mm:.{x_dp}f} mm\n\n"

        f"**Answer:**\n"
        f"The displacement transmissibility ratio is "
        f"**{transmissibility_ratio:.{tr_dp}f}**.\n"
        f"The absolute amplitude of the instrument's vibration is "
        f"**{amplitude_X_mm:.{x_dp}f} mm**."
    )

    return question, solution


# Template 4 (Advanced)
def template_vibration_isolator_design():
    """
    Vibrations: Vibration Isolator Design (Inverse Problem)

    Scenario:
        This template creates an advanced design problem. An engine generating a harmonic
        force needs to be mounted on isolators to limit the force transmitted to its
        foundation. Given a maximum allowable force transmission percentage, the task
        is to determine the necessary stiffness of the isolation system. This is an
        inverse problem, requiring algebraic manipulation to solve for a system parameter.

    Core Equations:
        TR = sqrt( (1 + (2*zeta*r)^2) / ( (1 - r^2)^2 + (2*zeta*r)^2 ) )
        This is solved for r, which is then used to find omega_n and finally k.
        k = m * omega_n^2

    Returns:
        tuple: A tuple containing:
            - str: A question asking for the required stiffness of an isolator.
            - str: A step-by-step solution to the design problem.
    """
    # 1. Parameterize the inputs with random values
    
    # Mass of the engine in kg
    mass = round(random.uniform(100.0, 2000.0), 1)
    
    # Operating speed in RPM
    operating_speed_rpm = random.randint(500, 3000)
    
    # Desired force transmissibility in percent
    transmissibility_percent = random.randint(5, 20)
    
    # Assumed damping ratio for the isolators (typically low)
    damping_ratio_zeta = round(random.uniform(0.05, 0.25), 3)

    # Standardize precision for final outputs
    precision = 4

    # 2. Perform the core calculations for the solution
    
    # Step A: Convert inputs to consistent units
    omega = operating_speed_rpm * (2 * math.pi / 60)
    transmissibility_ratio_TR = transmissibility_percent / 100.0

    # Step B: Solve for the frequency ratio (r)
    # The equation TR^2 = (1 + (2*zeta*r)^2) / ((1-r^2)^2 + (2*zeta*r)^2)
    # rearranges into a quadratic equation in terms of r^2: A*(r^2)^2 + B*(r^2) + C = 0
    TR_sq = transmissibility_ratio_TR**2
    
    A = TR_sq
    B = 4 * (damping_ratio_zeta**2) * (TR_sq - 1) - 2 * TR_sq
    C = TR_sq - 1

    # Calculate the discriminant
    discriminant = B**2 - 4 * A * C
    
    # Ensure a real solution exists (handles complex roots)
    if discriminant < 0:
        # This case is highly unlikely with TR < 1, but it's good practice to handle it.
        return ("Error: No real solution for frequency ratio (complex roots).",
                "The design parameters are not physically achievable.")
        
    # Solve the quadratic equation for r^2
    r_sq_sol1 = (-B + math.sqrt(discriminant)) / (2 * A)
    r_sq_sol2 = (-B - math.sqrt(discriminant)) / (2 * A)
    
    # For effective isolation (transmissibility TR < 1), the frequency ratio 'r'
    # must be greater than sqrt(2). We therefore need the larger, positive root for r^2.
    r_squared = max(r_sq_sol1, r_sq_sol2)

    # Add validation check for non-physical results (negative roots for r^2)
    if r_squared < 0:
        return ("Error: Design parameters result in a non-physical solution (r^2 < 0).",
                "The required transmissibility cannot be achieved with the given damping.")
    
    freq_ratio_r = math.sqrt(r_squared)
    
    # Step C: Calculate the required natural frequency (omega_n)
    omega_n_req = omega / freq_ratio_r
    
    # Step D: Calculate the required stiffness (k)
    stiffness_req = mass * (omega_n_req**2)

    # 3. Generate the question and solution strings
    
    question = (
        f"An engine with a mass of {mass} kg operates at a constant speed of {operating_speed_rpm} RPM. "
        f"It needs to be mounted on a set of vibration isolators. The design specification requires that no more than "
        f"{transmissibility_percent}% of the engine's unbalanced force is transmitted to the foundation.\n\n"
        f"Assuming the isolators have a combined damping ratio of {damping_ratio_zeta}, "
        f"determine the total required stiffness (k) of the isolation system."
    )

    solution = (
        f"**Given:**\n"
        f"Mass (m) = {mass} kg\n"
        f"Operating Speed = {operating_speed_rpm} RPM\n"
        f"Maximum Transmissibility (TR) = {transmissibility_percent}% = {transmissibility_ratio_TR}\n"
        f"Damping Ratio (zeta) = {damping_ratio_zeta}\n\n"

        f"**Step 1:** Convert the operating speed to rad/s.\n"
        f"omega = {operating_speed_rpm} RPM * (2 * pi / 60) = {round(omega, precision)} rad/s\n\n"

        f"**Step 2:** Set up the force transmissibility equation to solve for the frequency ratio (r).\n"
        f"The formula is TR^2 = [1 + (2*zeta*r)^2] / [(1 - r^2)^2 + (2*zeta*r)^2]\n"
        f"Rearranging this gives a quadratic equation in the form A(r^2)^2 + B(r^2) + C = 0.\n"
        f"A = TR^2 = {round(TR_sq, precision)}\n"
        f"B = 4*zeta^2*(TR^2 - 1) - 2*TR^2 = 4*({damping_ratio_zeta}^2)*({round(TR_sq, precision)} - 1) - 2*{round(TR_sq, precision)} = {round(B, precision)}\n"
        f"C = TR^2 - 1 = {round(TR_sq, precision)} - 1 = {round(C, precision)}\n\n"
        
        f"**Step 3:** Solve the quadratic equation for r^2 using the formula r^2 = (-B +/- sqrt(B^2 - 4AC)) / 2A.\n"
        f"Discriminant (D) = B^2 - 4AC = ({round(B, precision)})^2 - 4*({round(A, precision)})*({round(C, precision)}) = {round(discriminant, precision)}\n"
        f"The two solutions for r^2 are: {round(r_sq_sol1, precision)} and {round(r_sq_sol2, precision)}.\n"
        f"For effective vibration isolation, the frequency ratio 'r' must be greater than sqrt(2) (approx 1.414). This requires us to select the larger of the two positive solutions for r^2.\n"
        f"Required r^2 = {round(r_squared, precision)}\n\n"
        
        f"**Step 4:** Calculate the required frequency ratio (r) and validate the isolation condition.\n"
        f"r = sqrt({round(r_squared, precision)}) = {round(freq_ratio_r, precision)}\n"
        f"Check: Is r > sqrt(2)? Yes, {round(freq_ratio_r, precision)} > 1.414. The condition for isolation is met.\n\n"
        
        f"**Step 5:** Determine the required natural frequency (omega_n) of the system.\n"
        f"Since r = omega / omega_n, the required omega_n = omega / r.\n"
        f"omega_n = {round(omega, precision)} / {round(freq_ratio_r, precision)} = {round(omega_n_req, precision)} rad/s\n\n"

        f"**Step 6:** Calculate the total required stiffness (k).\n"
        f"The natural frequency is defined by omega_n = sqrt(k / m). Therefore, k = m * omega_n^2.\n"
        f"k = {mass} kg * ({round(omega_n_req, precision)} rad/s)^2\n"
        f"k = {round(stiffness_req, 0):,.0f} N/m\n\n"
        
        f"**Answer:**\n"
        f"The total required stiffness for the isolation system is **{round(stiffness_req, 0):,.0f} N/m**.\n\n"
        f"*Note: In a practical application, this total stiffness would be distributed among several individual isolator mounts.*"
    )

    return question, solution


def main():
    """
    Generate numerous instances of each free harmonically excited vibrations 
    template with different random seeds and write the results to a JSONL file.
    """
    import json
    import os

    # Define the output path (Modify this path according to where you are running the code from)
    output_file = "testset/mechanical_engineering/vibrations_and_acoustics/harmonically_excited_vibrations.jsonl"

    # Create the directory if it doesn't exist
    os.makedirs(os.path.dirname(output_file), exist_ok=True)

    # List of template functions with their ID and level
    templates = [
        (template_system_properties, "system_properties", "Easy"),
        (template_rotating_unbalance, "rotating_unbalance", "Intermediate"),
        (template_vibration_transmissibility, "vibration_transmissibility", "Intermediate"),
        (template_vibration_isolator_design, "vibration_isolator_design", "Advanced"),
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
                "domain": "vibrations_and_acoustics",
                "area": "harmonically_excited_vibrations",
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
