import random
import math


# Template 1 (Easy)
def template_undamped_natural_frequency_translational():
    """
    Free Vibration: Undamped Natural Frequency of a Translational System

    Scenario:
        This template generates a fundamental problem on a simple spring-mass
        system. It tests the ability to calculate the key characteristics of
        its undamped free vibration: the natural frequency in both rad/s (omega_n)
        and Hertz (f_n), and the natural period (tau_n).

    Core Equations:
        omega_n = sqrt(k / m)
        f_n = omega_n / (2 * pi)
        tau_n = 1 / f_n

    Returns:
        tuple: A tuple containing:
            - str: A question asking for the system's natural frequencies and period.
            - str: A step-by-step solution to the problem.
    """
    # 1. Parameterize the inputs with random values for high diversity
    
    # Mass (m): Chosen from a wide range to represent different scales,
    # from small mechanical parts to larger objects.
    mass = round(random.uniform(0.5, 750.0), 2)  # in kg

    # Stiffness (k): Integer values representing a typical range for mechanical springs.
    stiffness = random.randint(500, 250000)      # in N/m
    
    # Standardize precision for all calculations and final outputs
    precision = 3

    # 2. Perform the core calculations for the solution
    
    # Step A: Calculate the undamped natural frequency in radians per second
    omega_n = math.sqrt(stiffness / mass)
    
    # Step B: Convert the natural frequency from rad/s to Hertz (Hz)
    f_n = omega_n / (2 * math.pi)
    
    # Step C: Calculate the natural period of oscillation
    tau_n = 1 / f_n

    # 3. Generate the question and solution strings
    
    question = (
        f"An undamped single-degree-of-freedom system consists of a mass of {mass} kg "
        f"and a spring with a stiffness of {stiffness} N/m. "
        f"Determine the system's undamped natural frequency in both radians per second (rad/s) and Hertz (Hz). "
        f"Also, calculate the natural period of vibration."
    )

    solution = (
        f"**Given:**\n"
        f"Mass (m) = {mass} kg\n"
        f"Spring Stiffness (k) = {stiffness} N/m\n\n"

        f"**Step 1:** Calculate the undamped natural frequency in radians per second (omega_n).\n"
        f"The formula is: omega_n = sqrt(k / m)\n"
        f"omega_n = sqrt({stiffness} / {mass})\n"
        f"omega_n = {round(omega_n, precision)} rad/s\n\n"

        f"**Step 2:** Convert the natural frequency to Hertz (f_n).\n"
        f"The formula is: f_n = omega_n / (2 * pi)\n"
        f"f_n = {round(omega_n, precision)} / (2 * pi)\n"
        f"f_n = {round(f_n, precision)} Hz\n\n"

        f"**Step 3:** Calculate the natural period of vibration (tau_n).\n"
        f"The period is the reciprocal of the frequency in Hz.\n"
        f"The formula is: tau_n = 1 / f_n\n"
        f"tau_n = 1 / {round(f_n, precision)}\n"
        f"tau_n = {round(tau_n, precision)} s\n\n"

        f"**Answer:**\n"
        f"The undamped natural frequency is {round(omega_n, precision)} rad/s or {round(f_n, precision)} Hz.\n"
        f"The natural period of vibration is {round(tau_n, precision)} seconds."
    )

    return question, solution


# Template 2 (Easy)
def template_undamped_natural_frequency_torsional():
    """
    Free Vibration: Undamped Natural Frequency of a Torsional System

    Scenario:
        This template generates a problem for a simple torsional system, which
        typically consists of a disc or flywheel attached to a shaft. It is the
        rotational equivalent of the spring-mass system. The goal is to calculate
        the natural frequency and period of torsional oscillation.

    Core Equations:
        omega_n = sqrt(k_t / J_0)
        f_n = omega_n / (2 * pi)
        tau_n = 1 / f_n

    Returns:
        tuple: A tuple containing:
            - str: A question asking for the system's torsional natural frequencies and period.
            - str: A step-by-step solution to the problem.
    """
    # 1. Parameterize the inputs with random values
    
    # Mass moment of inertia (J_0): Represents the rotational inertia of the disc.
    # The range covers small rotors to medium-sized flywheels.
    mass_moment_of_inertia = round(random.uniform(0.05, 25.0), 3)  # in kg-m^2

    # Torsional stiffness (k_t): Represents the shaft's resistance to twisting.
    torsional_stiffness = random.randint(200, 75000)                # in N-m/rad
    
    # Standardize precision for all calculations and final outputs
    precision = 3

    # 2. Perform the core calculations for the solution
    
    # Step A: Calculate the undamped natural frequency in radians per second
    omega_n = math.sqrt(torsional_stiffness / mass_moment_of_inertia)
    
    # Step B: Convert the natural frequency from rad/s to Hertz (Hz)
    f_n = omega_n / (2 * math.pi)
    
    # Step C: Calculate the natural period of oscillation
    tau_n = 1 / f_n

    # 3. Generate the question and solution strings
    
    question = (
        f"A torsional pendulum consists of a disc with a mass moment of inertia of "
        f"{mass_moment_of_inertia} kg-m^2 attached to a shaft with a torsional stiffness of "
        f"{torsional_stiffness} N-m/rad. The system is undamped. "
        f"Calculate the natural frequency of torsional vibration in both rad/s and Hz, "
        f"and determine the corresponding period."
    )

    solution = (
        f"**Given:**\n"
        f"Mass Moment of Inertia (J_0) = {mass_moment_of_inertia} kg-m^2\n"
        f"Torsional Stiffness (k_t) = {torsional_stiffness} N-m/rad\n\n"

        f"**Step 1:** Calculate the undamped natural frequency in radians per second (omega_n).\n"
        f"The formula for a torsional system is: omega_n = sqrt(k_t / J_0)\n"
        f"omega_n = sqrt({torsional_stiffness} / {mass_moment_of_inertia})\n"
        f"omega_n = {round(omega_n, precision)} rad/s\n\n"

        f"**Step 2:** Convert the natural frequency to Hertz (f_n).\n"
        f"The relationship is: f_n = omega_n / (2 * pi)\n"
        f"f_n = {round(omega_n, precision)} / (2 * pi)\n"
        f"f_n = {round(f_n, precision)} Hz\n\n"

        f"**Step 3:** Calculate the natural period of vibration (tau_n).\n"
        f"The period is the reciprocal of the frequency in Hz.\n"
        f"The formula is: tau_n = 1 / f_n\n"
        f"tau_n = 1 / {round(f_n, precision)}\n"
        f"tau_n = {round(tau_n, precision)} s\n\n"

        f"**Answer:**\n"
        f"The undamped torsional natural frequency is {round(omega_n, precision)} rad/s or {round(f_n, precision)} Hz.\n"
        f"The natural period of vibration is {round(tau_n, precision)} seconds."
    )

    return question, solution


# Template 3 (Intermediate)
def template_undamped_response_initial_conditions():
    """
    Free Vibration: Response of an Undamped System to Initial Conditions

    Scenario:
        This template assesses the ability to determine the specific equation
        of motion, x(t), for an undamped spring-mass system given a set of
        initial conditions (displacement and velocity). This requires first
        finding the natural frequency and then solving for the two constants
        in the general solution.

    Core Equations:
        omega_n = sqrt(k / m)
        General Solution: x(t) = A1 * cos(omega_n * t) + A2 * sin(omega_n * t)
        From initial conditions: A1 = x(0) and A2 = v(0) / omega_n

    Returns:
        tuple: A tuple containing:
            - str: A question asking for the equation of motion.
            - str: A step-by-step solution deriving the equation.
    """
    # 1. Parameterize the inputs with random values
    mass = round(random.uniform(1.0, 500.0), 2)       # in kg
    stiffness = random.randint(1000, 200000)          # in N/m

    # Randomize initial conditions, including cases where one might be zero.
    # Displacement is given in mm for the question, velocity in m/s.
    initial_disp_mm = 0
    initial_vel_ms = 0.0
    
    # Use a chooser to create varied scenarios
    # 1: Only displacement, 2: Only velocity, 3: Both are non-zero
    scenario_choice = random.randint(1, 3)
    if scenario_choice == 1:
        initial_disp_mm = random.randint(-100, 100)
        while initial_disp_mm == 0: initial_disp_mm = random.randint(-100, 100)
    elif scenario_choice == 2:
        initial_vel_ms = round(random.uniform(-5.0, 5.0), 2)
        while initial_vel_ms == 0.0: initial_vel_ms = round(random.uniform(-5.0, 5.0), 2)
    else: # scenario_choice == 3
        initial_disp_mm = random.randint(-100, 100)
        initial_vel_ms = round(random.uniform(-5.0, 5.0), 2)
        while initial_disp_mm == 0: initial_disp_mm = random.randint(-100, 100)
        while initial_vel_ms == 0.0: initial_vel_ms = round(random.uniform(-5.0, 5.0), 2)
    
    # Convert initial displacement from mm to meters for calculations
    initial_disp_m = initial_disp_mm / 1000.0
    
    precision = 4

    # 2. Perform the core calculations for the solution
    
    # Step A: Calculate the natural frequency
    omega_n = math.sqrt(stiffness / mass)
    
    # Step B: Determine the constants A1 and A2 from initial conditions
    A1 = initial_disp_m
    A2 = initial_vel_ms / omega_n

    # 3. Generate the question and solution strings
    
    question = (
        f"An undamped spring-mass system has a mass of {mass} kg and a spring stiffness of {stiffness} N/m. "
        f"The mass is given an initial displacement of {initial_disp_mm} mm and an initial velocity of {initial_vel_ms} m/s. "
        f"Determine the equation of motion, x(t), for the system."
    )

    solution = (
        f"**Given:**\n"
        f"Mass (m) = {mass} kg\n"
        f"Stiffness (k) = {stiffness} N/m\n"
        f"Initial Displacement (x(0)) = {initial_disp_mm} mm = {initial_disp_m} m\n"
        f"Initial Velocity (v(0)) = {initial_vel_ms} m/s\n\n"

        f"**Step 1:** Calculate the undamped natural frequency (omega_n).\n"
        f"omega_n = sqrt(k / m) = sqrt({stiffness} / {mass}) = {round(omega_n, precision)} rad/s\n\n"

        f"**Step 2:** State the general form of the solution for an undamped system.\n"
        f"The general solution is: x(t) = A1 * cos(omega_n * t) + A2 * sin(omega_n * t)\n\n"

        f"**Step 3:** Apply the initial conditions to find the constants A1 and A2.\n"
        f"First, apply the initial displacement at t=0:\n"
        f"x(0) = A1 * cos(0) + A2 * sin(0) = A1\n"
        f"Therefore, A1 = x(0) = {initial_disp_m} m\n\n"
        
        f"Next, find the derivative of x(t) to get the velocity, v(t):\n"
        f"v(t) = dx/dt = -A1 * omega_n * sin(omega_n * t) + A2 * omega_n * cos(omega_n * t)\n"
        f"Now apply the initial velocity at t=0:\n"
        f"v(0) = -A1 * omega_n * sin(0) + A2 * omega_n * cos(0) = A2 * omega_n\n"
        f"Therefore, A2 = v(0) / omega_n = {initial_vel_ms} / {round(omega_n, precision)} = {round(A2, precision)} m\n\n"

        f"**Step 4:** Substitute the constants and omega_n into the general solution.\n"
        f"x(t) = {round(A1, precision)} * cos({round(omega_n, precision)} * t) + ({round(A2, precision)}) * sin({round(omega_n, precision)} * t)\n\n"
    )

    # Clean up the final equation for better readability
    cos_term = ""
    if abs(A1) > 1e-9: # Use a small tolerance to handle floating point inaccuracies
        cos_term = f"{round(A1, precision)}*cos({round(omega_n, precision)}*t)"

    sin_term = ""
    if abs(A2) > 1e-9:
        sign = "-" if A2 < 0 else "+"
        # If the cos_term is empty, don't start with a plus sign
        if not cos_term and A2 > 0:
            sign = ""
        
        sin_term = f" {sign} {abs(round(A2, precision))}*sin({round(omega_n, precision)}*t)"
        # Remove leading space if it's the first term
        if not cos_term:
            sin_term = sin_term.lstrip()
            
    final_equation = cos_term + sin_term
    if not final_equation: final_equation = "0"
    
    solution += (
        f"**Answer:**\n"
        f"The final equation of motion is:\n"
        f"x(t) = {final_equation} (m)"
    )

    return question, solution


# Template 4 (Intermediate)
def template_damping_classification():
    """
    Free Vibration: Damping Parameters and System Classification

    Scenario:
        This template introduces viscous damping and requires the calculation of
        key damping parameters. Given a system's physical properties (mass,
        stiffness, damping coefficient), the goal is to determine the damping
        ratio (zeta), classify the system's behavior (underdamped, critically
        damped, or overdamped), and calculate the damped natural frequency
        (omega_d) if the system is underdamped.

    Core Equations:
        omega_n = sqrt(k / m)
        c_c = 2 * sqrt(k * m)
        zeta = c / c_c
        omega_d = omega_n * sqrt(1 - zeta^2)  (for zeta < 1)

    Returns:
        tuple: A tuple containing:
            - str: A question asking for system classification and key parameters.
            - str: A step-by-step solution.
    """
    # 1. Parameterize the inputs to ensure all cases are generated
    mass = round(random.uniform(2.0, 600.0), 2)       # in kg
    stiffness = random.randint(2000, 300000)          # in N/m

    # To ensure a good distribution of outcomes, we first select a damping ratio
    # and then calculate the corresponding damping coefficient 'c'.
    
    # Choose a scenario: 1 for underdamped, 2 for critically damped, 3 for overdamped
    scenario_choice = random.randint(1, 3)

    if scenario_choice == 1: # Underdamped
        zeta = round(random.uniform(0.15, 0.85), 3)
    elif scenario_choice == 2: # Critically Damped
        zeta = 1.0
    else: # Overdamped
        zeta = round(random.uniform(1.2, 2.5), 3)

    # Standardize precision for final outputs
    precision = 4

    # 2. Perform the core calculations for the solution
    
    # Step A: Calculate the critical damping coefficient (c_c)
    critical_damping_c = 2 * math.sqrt(stiffness * mass)
    
    # Step B: Calculate the actual damping coefficient (c) based on the desired zeta
    damping_coefficient_c = zeta * critical_damping_c
    
    # Step C: Determine the system classification based on zeta
    if zeta < 1:
        classification = "Underdamped"
        # Also calculate natural and damped frequencies for the solution
        omega_n = math.sqrt(stiffness / mass)
        omega_d = omega_n * math.sqrt(1 - zeta**2)
    elif zeta == 1:
        classification = "Critically Damped"
    else:
        classification = "Overdamped"

    # 3. Generate the question and solution strings
    
    question = (
        f"A damped single-degree-of-freedom system has a mass of {mass} kg, "
        f"a spring stiffness of {stiffness} N/m, and a viscous damping coefficient of "
        f"{round(damping_coefficient_c, 2)} N-s/m. "
        f"Calculate the damping ratio (zeta) and classify the system as underdamped, "
        f"critically damped, or overdamped. If the system is underdamped, also "
        f"calculate its damped natural frequency (omega_d)."
    )

    solution = (
        f"**Given:**\n"
        f"Mass (m) = {mass} kg\n"
        f"Stiffness (k) = {stiffness} N/m\n"
        f"Damping Coefficient (c) = {round(damping_coefficient_c, 2)} N-s/m\n\n"

        f"**Step 1:** Calculate the critical damping coefficient (c_c).\n"
        f"The formula is: c_c = 2 * sqrt(k * m)\n"
        f"c_c = 2 * sqrt({stiffness} * {mass})\n"
        f"c_c = {round(critical_damping_c, precision)} N-s/m\n\n"

        f"**Step 2:** Calculate the damping ratio (zeta).\n"
        f"The formula is: zeta = c / c_c\n"
        f"zeta = {round(damping_coefficient_c, 2)} / {round(critical_damping_c, precision)}\n"
        f"zeta = {round(zeta, precision)}\n\n"

        f"**Step 3:** Classify the system based on the value of zeta.\n"
        f"Since zeta {'>' if zeta > 1 else '<' if zeta < 1 else '='} 1, the system is **{classification}**.\n\n"
    )

    # Add the final step only if the system is underdamped
    if classification == "Underdamped":
        solution += (
            f"**Step 4:** Since the system is underdamped, calculate the damped natural frequency (omega_d).\n"
            f"First, find the undamped natural frequency (omega_n):\n"
            f"omega_n = sqrt(k / m) = sqrt({stiffness} / {mass}) = {round(omega_n, precision)} rad/s\n\n"
            f"Now, use the formula: omega_d = omega_n * sqrt(1 - zeta^2)\n"
            f"omega_d = {round(omega_n, precision)} * sqrt(1 - {round(zeta, precision)}^2)\n"
            f"omega_d = {round(omega_d, precision)} rad/s\n\n"
            
            f"**Answer:**\n"
            f"The damping ratio is {round(zeta, precision)}. The system is **{classification}**.\n"
            f"The damped natural frequency is {round(omega_d, precision)} rad/s."
        )
    else:
        solution += (
            f"**Answer:**\n"
            f"The damping ratio is {round(zeta, precision)}. The system is **{classification}**."
        )

    return question, solution


# Template 5 (Intermediate)
def template_logarithmic_decrement():
    """
    Free Vibration: Logarithmic Decrement and Damping Ratio

    Scenario:
        This template models a practical method for determining the damping in
        an underdamped system. By observing the amplitude of vibration at two
        distinct points in time, separated by a known number of cycles, one can
        calculate the logarithmic decrement, which directly relates to the system's
        damping ratio.

    Core Equations:
        delta = (1/n) * ln(x1 / x_{n+1})
        zeta = delta / sqrt((2*pi)^2 + delta^2)

    Returns:
        tuple: A tuple containing:
            - str: A question asking for the logarithmic decrement and damping ratio.
            - str: A step-by-step solution.
    """
    # 1. Parameterize the inputs with random values
    
    # Initialize variables to enter the loop
    final_amplitude_x_n_plus_1 = 0.0
    
    # Use a loop to ensure the final rounded amplitude is never zero.
    # This prevents a division-by-zero error in cases of high decay.
    while final_amplitude_x_n_plus_1 <= 0.0:
        # To ensure physically consistent values, we first generate a realistic
        # damping ratio (zeta) and work backward to find the amplitudes.
        zeta_actual = random.uniform(0.05, 0.25)
        
        # Calculate the corresponding logarithmic decrement
        delta_actual = (2 * math.pi * zeta_actual) / math.sqrt(1 - zeta_actual**2)
        
        # Choose a number of cycles for the measurement
        num_cycles_n = random.randint(2, 8) # Reduced max cycles to lower chance of zeroing out
        
        # Set a random initial amplitude in mm
        initial_amplitude_x1 = round(random.uniform(20.0, 100.0), 1)
        
        # Calculate the final amplitude based on the decrement and number of cycles
        calculated_final_amplitude = initial_amplitude_x1 / math.exp(num_cycles_n * delta_actual)
        
        # Round the final amplitude to make it look like a measured value
        final_amplitude_x_n_plus_1 = round(calculated_final_amplitude, 1)

    # Standardize precision for final outputs in the solution
    precision = 4

    # 2. Perform the core calculations for the solution (as the user would)
    
    # Step A: Calculate the logarithmic decrement from the given amplitudes
    log_decrement_delta = (1 / num_cycles_n) * math.log(initial_amplitude_x1 / final_amplitude_x_n_plus_1)
    
    # Step B: Calculate the damping ratio from the decrement
    damping_ratio_zeta = log_decrement_delta / math.sqrt((2 * math.pi)**2 + log_decrement_delta**2)

    # 3. Generate the question and solution strings
    
    question = (
        f"The amplitude of free vibration of an underdamped system is observed to decay over time. "
        f"The initial amplitude is measured to be {initial_amplitude_x1} mm. After {num_cycles_n} complete cycles, "
        f"the amplitude is measured to be {final_amplitude_x_n_plus_1} mm. "
        f"Based on these observations, calculate the logarithmic decrement (delta) and the damping ratio (zeta) of the system."
    )

    solution = (
        f"**Given:**\n"
        f"Initial Amplitude (x1) = {initial_amplitude_x1} mm\n"
        f"Amplitude after {num_cycles_n} cycles (x{num_cycles_n + 1}) = {final_amplitude_x_n_plus_1} mm\n"
        f"Number of cycles (n) = {num_cycles_n}\n\n"

        f"**Step 1:** Calculate the logarithmic decrement (delta).\n"
        f"The formula is: delta = (1/n) * ln(x1 / x_{num_cycles_n + 1})\n"
        f"delta = (1 / {num_cycles_n}) * ln({initial_amplitude_x1} / {final_amplitude_x_n_plus_1})\n"
        f"delta = (1 / {num_cycles_n}) * ln({round(initial_amplitude_x1 / final_amplitude_x_n_plus_1, precision)})\n"
        f"delta = {round(log_decrement_delta, precision)}\n\n"

        f"**Step 2:** Calculate the damping ratio (zeta) from the logarithmic decrement.\n"
        f"The formula relating zeta and delta is: zeta = delta / sqrt((2*pi)^2 + delta^2)\n"
        f"zeta = {round(log_decrement_delta, precision)} / sqrt((2*pi)^2 + {round(log_decrement_delta, precision)}^2)\n"
        f"zeta = {round(damping_ratio_zeta, precision)}\n\n"

        f"**Answer:**\n"
        f"The logarithmic decrement is {round(log_decrement_delta, precision)}.\n"
        f"The damping ratio of the system is {round(damping_ratio_zeta, precision)}."
    )

    return question, solution


# Template 6 (Advanced)
def template_equivalent_stiffness_frequency():
    """
    Free Vibration: Equivalent Stiffness and Natural Frequency

    Scenario:
        This template addresses systems with multiple springs. It requires first
        simplifying a spring network (either in series or parallel) into a single
        equivalent spring. After finding the equivalent stiffness (k_eq), the
        natural frequency of the overall system is calculated. This adds a
        critical modeling step to the standard frequency calculation.

    Core Equations:
        Parallel Stiffness: k_eq = k1 + k2
        Series Stiffness: k_eq = (k1 * k2) / (k1 + k2)
        Natural Frequency: omega_n = sqrt(k_eq / m)

    Returns:
        tuple: A tuple containing:
            - str: A question asking for the equivalent stiffness and natural frequency.
            - str: A step-by-step solution.
    """
    # 1. Parameterize the inputs with random values
    mass = round(random.uniform(5.0, 100.0), 2)       # in kg
    stiffness_1 = random.randint(1000, 50000)         # in N/m
    stiffness_2 = random.randint(1000, 50000)         # in N/m
    
    # Randomly choose the spring configuration
    configuration = random.choice(['series', 'parallel'])
    
    # Standardize precision for final outputs
    precision = 3

    # 2. Perform the core calculations for the solution
    
    # Step A: Calculate the equivalent stiffness (k_eq) based on the configuration
    if configuration == 'parallel':
        k_eq = stiffness_1 + stiffness_2
        formula_str = "k_eq = k1 + k2"
        calculation_str = f"k_eq = {stiffness_1} + {stiffness_2} = {k_eq} N/m"
    else: # configuration == 'series'
        k_eq = (stiffness_1 * stiffness_2) / (stiffness_1 + stiffness_2)
        formula_str = "k_eq = (k1 * k2) / (k1 + k2)"
        calculation_str = (
            f"k_eq = ({stiffness_1} * {stiffness_2}) / ({stiffness_1} + {stiffness_2}) = "
            f"{round(k_eq, 2)} N/m"
        )
        
    # Step B: Calculate the natural frequency using the equivalent stiffness
    omega_n = math.sqrt(k_eq / mass)

    # 3. Generate the question and solution strings
    
    question = (
        f"A mass of {mass} kg is attached to a system of two springs. The first spring has a "
        f"stiffness (k1) of {stiffness_1} N/m, and the second spring has a stiffness (k2) of {stiffness_2} N/m. "
        f"The two springs are arranged in {configuration}. "
        f"Determine the equivalent stiffness of the spring system and the resulting undamped "
        f"natural frequency (omega_n) in rad/s."
    )

    solution = (
        f"**Given:**\n"
        f"Mass (m) = {mass} kg\n"
        f"Stiffness 1 (k1) = {stiffness_1} N/m\n"
        f"Stiffness 2 (k2) = {stiffness_2} N/m\n"
        f"Configuration = {configuration.capitalize()}\n\n"

        f"**Step 1:** Calculate the equivalent stiffness (k_eq) for the {configuration} configuration.\n"
        f"The formula for springs in {configuration} is: {formula_str}\n"
        f"{calculation_str}\n\n"

        f"**Step 2:** Calculate the undamped natural frequency (omega_n) using the equivalent stiffness.\n"
        f"The formula is: omega_n = sqrt(k_eq / m)\n"
        f"omega_n = sqrt({round(k_eq, 2)} / {mass})\n"
        f"omega_n = {round(omega_n, precision)} rad/s\n\n"

        f"**Answer:**\n"
        f"The equivalent stiffness of the system is {round(k_eq, 2)} N/m.\n"
        f"The undamped natural frequency is {round(omega_n, precision)} rad/s."
    )

    return question, solution


def main():
    """
    Generate numerous instances of each free vibration of single-degree-of-freedom systems 
    template with different random seeds and write the results to a JSONL file.
    """
    import json
    import os

    # Define the output path (Modify this path according to where you are running the code from)
    output_file = "testset/mechanical_engineering/vibrations_and_acoustics/single_degree_of_freedom_systems.jsonl"

    # Create the directory if it doesn't exist
    os.makedirs(os.path.dirname(output_file), exist_ok=True)

    # List of template functions with their ID and level
    templates = [
        (template_undamped_natural_frequency_translational, "undamped_natural_frequency_translational", "Easy"),
        (template_undamped_natural_frequency_torsional, "undamped_natural_frequency_torsional", "Easy"),
        (template_undamped_response_initial_conditions, "undamped_response_initial_conditions", "Intermediate"),
        (template_damping_classification, "damping_classification", "Intermediate"),
        (template_logarithmic_decrement, "logarithmic_decrement", "Intermediate"),
        (template_equivalent_stiffness_frequency, "equivalent_stiffness_frequency", "Advanced"),
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
