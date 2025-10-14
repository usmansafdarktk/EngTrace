import random
import math


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

    Returns:
        tuple: A tuple containing:
            - str: A question asking for the steady-state amplitude of vibration.
            - str: A step-by-step solution to the problem.
    """
    # 1. Parameterize the inputs with random values
    
    # Total mass of the machine in kg
    m_total = round(random.uniform(50.0, 500.0), 1)
    
    # Eccentric mass in kg (should be a small fraction of total mass)
    m_eccentric = round(random.uniform(0.1, m_total * 0.05), 2)
    
    # Eccentricity in mm
    eccentricity_mm = random.randint(10, 150)
    
    # Stiffness in N/m
    stiffness = round(random.uniform(5e4, 2e6), 0)
    
    # Define system properties by choosing a damping ratio first for better control
    # Lightly damped systems are common for this problem type.
    damping_ratio_zeta = round(random.uniform(0.05, 0.6), 3)
    
    # Define operating speed by choosing a frequency ratio first
    # This ensures we get a good spread of cases (below, near, and above resonance)
    freq_ratio_r = round(random.uniform(0.3, 3.0), 3)

    # Standardize precision for final outputs
    precision = 5
    
    # 2. Perform the core calculations for the solution
    
    # Handle potential edge cases from inputs
    if m_total <= 0 or stiffness <= 0:
        return "Error: Mass and stiffness must be positive.", "Invalid input parameters."

    # Step A: Calculate fundamental system properties
    omega_n = math.sqrt(stiffness / m_total)
    c_critical = 2 * math.sqrt(stiffness * m_total)
    damping_coeff = damping_ratio_zeta * c_critical
    
    # Step B: Determine the operating speed from the chosen frequency ratio
    omega = freq_ratio_r * omega_n
    operating_speed_rpm = omega * 60 / (2 * math.pi)
    
    # Step C: Convert units for calculation consistency
    eccentricity_m = eccentricity_mm / 1000.0
    
    # Step D: Calculate the magnitude of the unbalanced force
    force_magnitude_F0 = m_eccentric * eccentricity_m * (omega ** 2)
    
    # Step E: Calculate the steady-state amplitude (X) in meters
    # Denominator components
    term1 = stiffness - m_total * (omega ** 2)
    term2 = damping_coeff * omega
    denominator = math.sqrt(term1**2 + term2**2)
    
    amplitude_m = force_magnitude_F0 / denominator if denominator != 0 else float('inf')
    
    # Step F: Convert final amplitude to millimeters for a more intuitive answer
    amplitude_mm = amplitude_m * 1000.0

    # 3. Generate the question and solution strings
    
    question = (
        f"A machine with a total mass of {m_total} kg is supported by a spring and damper system. "
        f"The system has an equivalent stiffness of {stiffness:,.0f} N/m and an equivalent damping "
        f"coefficient of {round(damping_coeff, 2)} N.s/m.\n\n"
        f"The machine contains a rotating component that has an unbalance equivalent to a mass of "
        f"{m_eccentric} kg located at an eccentricity of {eccentricity_mm} mm. "
        f"If the machine operates at a speed of {round(operating_speed_rpm, 0):,.0f} RPM, "
        f"determine the steady-state amplitude of vibration."
    )

    solution = (
        f"**Given:**\n"
        f"Total Mass (m) = {m_total} kg\n"
        f"Stiffness (k) = {stiffness:,.0f} N/m\n"
        f"Damping Coefficient (c) = {round(damping_coeff, 2)} N.s/m\n"
        f"Eccentric Mass (m_e) = {m_eccentric} kg\n"
        f"Eccentricity (e) = {eccentricity_mm} mm = {eccentricity_m} m\n"
        f"Operating Speed = {round(operating_speed_rpm, 0):,.0f} RPM\n\n"

        f"**Step 1:** Convert the operating speed from RPM to rad/s.\n"
        f"omega = (Speed in RPM) * (2 * pi / 60)\n"
        f"omega = {round(operating_speed_rpm, 0):,.0f} * (2 * pi / 60) = {round(omega, precision)} rad/s\n\n"

        f"**Step 2:** Calculate the system's natural frequency (omega_n) and damping ratio (zeta).\n"
        f"omega_n = sqrt(k / m) = sqrt({stiffness:,.0f} / {m_total}) = {round(omega_n, precision)} rad/s\n"
        f"c_cr = 2 * sqrt(k * m) = 2 * sqrt({stiffness:,.0f} * {m_total}) = {round(c_critical, precision)} N.s/m\n"
        f"zeta = c / c_cr = {round(damping_coeff, 2)} / {round(c_critical, precision)} = {round(damping_ratio_zeta, precision)}\n\n"

        f"**Step 3:** Calculate the magnitude of the unbalanced force (F_0).\n"
        f"The force from a rotating unbalance is given by F_0 = m_e * e * omega^2.\n"
        f"F_0 = {m_eccentric} kg * {eccentricity_m} m * ({round(omega, precision)} rad/s)^2\n"
        f"F_0 = {round(force_magnitude_F0, precision)} N\n\n"
        
        f"**Step 4:** Calculate the steady-state amplitude of vibration (X).\n"
        f"The formula for amplitude is: X = F_0 / sqrt((k - m * omega^2)^2 + (c * omega)^2)\n"
        f"Numerator = F_0 = {round(force_magnitude_F0, precision)} N\n"
        f"Denominator Part 1: (k - m * omega^2) = ({stiffness:,.0f} - {m_total} * {round(omega, precision)}^2) = {round(term1, precision)}\n"
        f"Denominator Part 2: (c * omega) = ({round(damping_coeff, 2)} * {round(omega, precision)}) = {round(term2, precision)}\n"
        f"Denominator = sqrt(({round(term1, precision)})^2 + ({round(term2, precision)})^2) = {round(denominator, precision)}\n"
        f"X = {round(force_magnitude_F0, precision)} / {round(denominator, precision)}\n"
        f"X = {round(amplitude_m, precision + 2)} m\n\n"
        
        f"**Step 5:** Convert the amplitude to millimeters.\n"
        f"Amplitude in mm = {round(amplitude_m, precision + 2)} m * 1000 mm/m = {round(amplitude_mm, precision)} mm\n\n"
        
        f"**Answer:**\n"
        f"The steady-state amplitude of vibration is **{round(amplitude_mm, 3)} mm**."
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

    Returns:
        tuple: A tuple containing:
            - str: A question asking for the transmissibility and absolute amplitude.
            - str: A step-by-step solution to the problem.
    """
    # 1. Parameterize the inputs with random values
    
    # Mass of the sensitive instrument in kg
    mass = round(random.uniform(5.0, 150.0), 1)
    
    # Stiffness of the isolation mount in N/m
    stiffness = round(random.uniform(2e3, 5e5), 0)
    
    # Amplitude of the base vibration in mm
    base_amplitude_Y_mm = round(random.uniform(0.1, 8.0), 2)
    
    # Use frequency ratio and damping ratio to control the problem's outcome
    # This ensures we test conditions of amplification (r~1) and isolation (r > sqrt(2))
    freq_ratio_r = round(random.uniform(0.2, 5.0), 3)
    damping_ratio_zeta = round(random.uniform(0.05, 0.7), 3)
    
    # Standardize precision for final outputs
    precision = 4

    # 2. Perform the core calculations for the solution
    
    # Handle potential edge cases
    if mass <= 0 or stiffness <= 0:
        return "Error: Mass and stiffness must be positive.", "Invalid input parameters."

    # Step A: Calculate fundamental system properties
    omega_n = math.sqrt(stiffness / mass)
    c_critical = 2 * math.sqrt(stiffness * mass)
    damping_coeff = damping_ratio_zeta * c_critical

    # Step B: Determine the base excitation frequency from the chosen frequency ratio
    omega = freq_ratio_r * omega_n
    base_freq_hz = omega / (2 * math.pi)

    # Step C: Convert base amplitude to meters for calculation
    base_amplitude_Y_m = base_amplitude_Y_mm / 1000.0

    # Step D: Calculate the displacement transmissibility ratio (TR)
    tr_num = 1 + (2 * damping_ratio_zeta * freq_ratio_r)**2
    tr_den = (1 - freq_ratio_r**2)**2 + (2 * damping_ratio_zeta * freq_ratio_r)**2

    # Avoid division by zero, although highly unlikely with these random ranges
    if tr_den == 0:
        transmissibility_ratio = float('inf')
    else:
        transmissibility_ratio = math.sqrt(tr_num / tr_den)
    
    # Step E: Calculate the absolute amplitude of the instrument's vibration
    amplitude_X_m = transmissibility_ratio * base_amplitude_Y_m
    amplitude_X_mm = amplitude_X_m * 1000.0

    # 3. Generate the question and solution strings
    
    question = (
        f"A sensitive instrument of mass {mass} kg is supported by an isolation mount. "
        f"The mount has an effective stiffness of {stiffness:,.0f} N/m and provides a damping "
        f"ratio of {damping_ratio_zeta}.\n\n"
        f"The foundation on which the instrument is placed is vibrating harmonically at a frequency of "
        f"{round(base_freq_hz, 2)} Hz with an amplitude of {base_amplitude_Y_mm} mm.\n\n"
        f"Determine:\n"
        f"1. The displacement transmissibility ratio.\n"
        f"2. The absolute amplitude of vibration of the instrument in millimeters."
    )

    solution = (
        f"**Given:**\n"
        f"Mass (m) = {mass} kg\n"
        f"Stiffness (k) = {stiffness:,.0f} N/m\n"
        f"Damping Ratio (zeta) = {damping_ratio_zeta}\n"
        f"Base Vibration Frequency (f) = {round(base_freq_hz, 2)} Hz\n"
        f"Base Vibration Amplitude (Y) = {base_amplitude_Y_mm} mm = {base_amplitude_Y_m} m\n\n"

        f"**Step 1:** Calculate the system's undamped natural frequency (omega_n).\n"
        f"omega_n = sqrt(k / m) = sqrt({stiffness:,.0f} / {mass}) = {round(omega_n, precision)} rad/s\n\n"

        f"**Step 2:** Convert the base vibration frequency to rad/s and find the frequency ratio (r).\n"
        f"Base frequency (omega) = f * 2 * pi = {round(base_freq_hz, 2)} * 2 * pi = {round(omega, precision)} rad/s\n"
        f"Frequency ratio (r) = omega / omega_n = {round(omega, precision)} / {round(omega_n, precision)} = {round(freq_ratio_r, precision)}\n\n"

        f"**Step 3:** Calculate the displacement transmissibility ratio (TR).\n"
        f"The formula is: TR = sqrt( (1 + (2*zeta*r)^2) / ( (1 - r^2)^2 + (2*zeta*r)^2 ) )\n"
        f"Let's calculate the terms:\n"
        f"r = {round(freq_ratio_r, precision)}\n"
        f"zeta = {damping_ratio_zeta}\n"
        f"Numerator = 1 + (2 * {damping_ratio_zeta} * {round(freq_ratio_r, precision)})^2 = {round(tr_num, precision)}\n"
        f"Denominator = (1 - ({round(freq_ratio_r, precision)})^2)^2 + (2 * {damping_ratio_zeta} * {round(freq_ratio_r, precision)})^2 = {round(tr_den, precision)}\n"
        f"TR = sqrt({round(tr_num, precision)} / {round(tr_den, precision)}) = {round(transmissibility_ratio, precision)}\n\n"

        f"**Step 4:** Calculate the absolute amplitude of the instrument (X).\n"
        f"The relationship is X = TR * Y.\n"
        f"X = {round(transmissibility_ratio, precision)} * {base_amplitude_Y_mm} mm\n"
        f"X = {round(amplitude_X_mm, precision)} mm\n\n"

        f"**Answer:**\n"
        f"The displacement transmissibility ratio is **{round(transmissibility_ratio, 3)}**.\n"
        f"The absolute amplitude of the instrument's vibration is **{round(amplitude_X_mm, 3)} mm**."
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
