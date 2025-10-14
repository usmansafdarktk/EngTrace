import random
import math
from data.templates.branches.electrical_engineering.constants import C0, MEDIA_VELOCITIES 


# Template 1 (Easy)
def template_wave_parameters_basic():
    """
    Wave Parameter Fundamentals Calculation

    Scenario:
        This template explores the fundamental relationships that define a sinusoidal
        traveling wave. These parameters (frequency, wavelength, period, etc.) describe
        a wave's oscillatory behavior in both time and space. The problem provides the
        wave's propagation speed (defined by its medium) and one other key parameter
        (either frequency or wavelength), requiring the calculation of all other
        related properties.

    Core Equations:
        u_p = f * lambda  
        omega = 2 * pi * f  
        T = 1 / f          
        k = 2 * pi / lambda 

    Returns:
        tuple: A tuple containing:
            - str: A question asking to compute various wave parameters.
            - str: A step-by-step solution showing the calculations.
    """
    # 1. Parameterize the inputs with random values
    medium_name, u_p = random.choice(list(MEDIA_VELOCITIES.items()))
    start_with_frequency = random.choice([True, False])

    if start_with_frequency:
        # Scenario 1: Given frequency and medium
        f_mhz = round(random.uniform(50, 500), 1)
        f = f_mhz * 1e6

        # 2. Perform the core calculations
        omega = 2 * math.pi * f
        T = 1 / f
        lambda_ = u_p / f
        k = 2 * math.pi / lambda_

        # 3. Generate the question and solution strings (Plain Text)
        question = (
            f"A sinusoidal wave with a frequency of {f_mhz} MHz is propagating "
            f"through {medium_name}. Given that the phase velocity in this medium is "
            f"{u_p:.2e} m/s, calculate the wave's:\n"
            f"a) Angular frequency (omega)\n"
            f"b) Period (T)\n"
            f"c) Wavelength (lambda)\n"
            f"d) Wavenumber (k)"
        )

        solution = (
            f"**Given:**\n"
            f"- Frequency (f) = {f_mhz} MHz = {f:.2e} Hz\n"
            f"- Medium = {medium_name}\n"
            f"- Phase Velocity (u_p) = {u_p:.2e} m/s\n\n"
            
            f"**Step 1:** Calculate the Angular Frequency (omega)\n"
            f"The angular frequency is related to frequency by omega = 2 * pi * f.\n"
            f"   omega = 2 * {math.pi:.5f} * ({f:.2e} Hz) = {omega:.2e} rad/s\n\n"

            f"**Step 2:** Calculate the Period (T)\n"
            f"The period is the inverse of the frequency, T = 1 / f.\n"
            f"   T = 1 / ({f:.2e} Hz) = {T:.2e} s\n\n"

            f"**Step 3:** Calculate the Wavelength (lambda)\n"
            f"The wavelength is found using the phase velocity, lambda = u_p / f.\n"
            f"   lambda = ({u_p:.2e} m/s) / ({f:.2e} Hz) = {round(lambda_, 3)} m\n\n"

            f"**Step 4:** Calculate the Wavenumber (k)\n"
            f"The wavenumber (or phase constant) is given by k = 2 * pi / lambda.\n"
            f"   k = 2 * {math.pi:.5f} / {round(lambda_, 3)} m = {round(k, 2)} rad/m\n\n"

            f"**Answer:**\n"
            f"- Angular Frequency (omega): {omega:.2e} rad/s\n"
            f"- Period (T): {T:.2e} s\n"
            f"- Wavelength (lambda): {round(lambda_, 3)} m\n"
            f"- Wavenumber (k): {round(k, 2)} rad/m"
        )

    else:
        # Scenario 2: Given wavelength and medium
        lambda_ = round(random.uniform(0.1, 2.0), 2)

        # 2. Perform the core calculations
        f = u_p / lambda_
        omega = 2 * math.pi * f
        T = 1 / f
        k = 2 * math.pi / lambda_

        # 3. Generate the question and solution strings (Plain Text)
        question = (
            f"An electromagnetic wave traveling in {medium_name} is observed to have "
            f"a wavelength of {lambda_} m. Given that the phase velocity in this "
            f"medium is {u_p:.2e} m/s, determine the wave's:\n"
            f"a) Frequency (f)\n"
            f"b) Angular frequency (omega)\n"
            f"c) Period (T)\n"
            f"d) Wavenumber (k)"
        )

        solution = (
            f"**Given:**\n"
            f"- Wavelength (lambda) = {lambda_} m\n"
            f"- Medium = {medium_name}\n"
            f"- Phase Velocity (u_p) = {u_p:.2e} m/s\n\n"
            
            f"**Step 1:** Calculate the Frequency (f)\n"
            f"Frequency is found using the relation u_p = f * lambda, so f = u_p / lambda.\n"
            f"   f = ({u_p:.2e} m/s) / ({lambda_} m) = {f:.2e} Hz = {f/1e6:.2f} MHz\n\n"
            
            f"**Step 2:** Calculate the Angular Frequency (omega)\n"
            f"The angular frequency is omega = 2 * pi * f.\n"
            f"   omega = 2 * {math.pi:.5f} * ({f:.2e} Hz) = {omega:.2e} rad/s\n\n"

            f"**Step 3:** Calculate the Period (T)\n"
            f"The period is the inverse of the frequency, T = 1 / f.\n"
            f"   T = 1 / ({f:.2e} Hz) = {T:.2e} s\n\n"

            f"**Step 4:** Calculate the Wavenumber (k)\n"
            f"The wavenumber is given by k = 2 * pi / lambda.\n"
            f"   k = 2 * {math.pi:.5f} / {lambda_} m = {round(k, 2)} rad/m\n\n"
            
            f"**Answer:**\n"
            f"- Frequency (f): {f/1e6:.2f} MHz\n"
            f"- Angular Frequency (omega): {omega:.2e} rad/s\n"
            f"- Period (T): {T:.2e} s\n"
            f"- Wavenumber (k): {round(k, 2)} rad/m"
        )

    return question, solution


# Template 2 (Easy)
def template_time_to_phasor():
    """
    Time-Domain to Phasor Conversion

    Scenario:
        This template tests the ability to convert a sinusoidal time-domain function
        into its corresponding phasor representation. This is a foundational skill for
        AC circuit analysis and wave analysis. The conversion requires identifying the
        amplitude and phase, and applying a phase shift if the function is a sine wave.

    Core Equations:
        For v(t) = A * cos(omega*t + phi):
        Phasor V = A * exp(j*phi) = A < phi
        
        Identity: sin(theta) = cos(theta - 90 deg)

    Returns:
        tuple: A tuple containing:
            - str: A question asking for the phasor form of a given time-domain signal.
            - str: A step-by-step solution showing the conversion.
    """
    # 1. Parameterize the inputs with random values
    amplitude = round(random.uniform(5.0, 150.0), 2)
    omega = random.randint(100, 1000)
    func_type = random.choice(['cos', 'sin'])
    use_degrees = random.choice([True, False])

    # Standardize precision for all calculations and outputs
    precision = 2

    # Generate phase in either degrees or radians
    if use_degrees:
        phi_deg = round(random.uniform(-180.0, 180.0), 1)
        phi_rad = math.radians(phi_deg)
        phi_str = f"{phi_deg} deg"
    else:
        phi_rad = round(random.uniform(-math.pi, math.pi), precision)
        phi_deg = math.degrees(phi_rad)
        phi_str = f"{phi_rad} rad"

    # Construct the time-domain function string for the question
    time_domain_expr = f"{amplitude} * {func_type}({omega}*t + {phi_str})"

    # 2. Perform the core calculation
    
    # The final amplitude of the phasor is the amplitude of the signal
    phasor_amplitude = amplitude
    
    # The final phase depends on whether the function is cos or sin
    if func_type == 'cos':
        phasor_phi_deg = phi_deg
        phasor_phi_rad = phi_rad
        conversion_step = (
            "**Step 2:** Identify the function type.\n"
            "   The function is a cosine, which is the standard reference for phasors. "
            "No phase adjustment is needed.\n"
            f"   The initial phase is {round(phi_deg, precision)} degrees.\n"
        )
    else: # func_type == 'sin'
        phasor_phi_deg = phi_deg - 90
        phasor_phi_rad = math.radians(phasor_phi_deg)
        conversion_step = (
            "**Step 2:** Convert the sine function to cosine.\n"
            "   The sine function leads the cosine function by 90 degrees. To convert, "
            "we use the identity sin(theta) = cos(theta - 90 deg).\n"
            f"   New phase = (Initial Phase) - 90 deg = {round(phi_deg, precision)} - 90 = {round(phasor_phi_deg, precision)} degrees.\n"
        )
        
    # Normalize the final phase angle to be between -180 and 180 degrees
    phasor_phi_deg = (phasor_phi_deg + 180) % 360 - 180
    phasor_phi_rad = math.radians(phasor_phi_deg)
    
    # Calculate rectangular components
    real_part = phasor_amplitude * math.cos(phasor_phi_rad)
    imag_part = phasor_amplitude * math.sin(phasor_phi_rad)

    # 3. Generate the question and solution strings
    question = (
        f"A signal is described by the time-domain function:\n"
        f"   v(t) = {time_domain_expr}\n\n"
        f"Determine the phasor representation of this signal in both polar (A < phi) "
        f"and rectangular (x + jy) forms. Use degrees for the phase angle."
    )

    solution = (
        f"**Given:**\n"
        f"   The time-domain signal is v(t) = {time_domain_expr}.\n\n"
        
        f"**Step 1:** Extract the amplitude and initial phase.\n"
        f"   By inspection, the amplitude (A) is {amplitude}.\n"
        f"   The initial phase is {phi_str}.\n\n"
        
        f"{conversion_step}\n"
        
        f"**Step 3:** Write the phasor in polar form.\n"
        f"   Phase angles are conventionally expressed in the range [-180°, 180°].\n"
        f"   The phasor has the signal's amplitude and the adjusted phase.\n"
        f"   Phasor V = {round(phasor_amplitude, precision)} < {round(phasor_phi_deg, precision)} degrees.\n\n"
        
        f"**Step 4:** Convert the polar form to rectangular form (x + jy).\n"
        f"   x (real part) = A * cos(phi) = {round(phasor_amplitude, precision)} * cos({round(phasor_phi_deg, precision)} deg) = {round(real_part, precision)}\n"
        f"   y (imaginary part) = A * sin(phi) = {round(phasor_amplitude, precision)} * sin({round(phasor_phi_deg, precision)} deg) = {round(imag_part, precision)}\n"
        f"   Phasor V = {round(real_part, precision)} + j{round(imag_part, precision)}\n\n"

        f"**Answer:**\n"
        f"   The phasor representation is {round(phasor_amplitude, precision)} < {round(phasor_phi_deg, precision)} degrees, "
        f"which is equivalent to {round(real_part, precision)} + j{round(imag_part, precision)}."
    )

    return question, solution


# Template 3 (Intermediate)
def template_wave_equation_interpretation():
    """
    Wave Equation Interpretation

    Scenario:
        This template tests the ability to extract fundamental wave parameters directly
        from the mathematical expression of a traveling wave. It requires matching the
        given equation to the standard form to identify key coefficients (amplitude,
        angular frequency, wavenumber) and then using them to calculate related properties.

    Core Equations:
        f = omega / (2 * pi)
        lambda = 2 * pi / k
        u_p = omega / k

    Returns:
        tuple: A tuple containing:
            - str: A question presenting a wave equation and asking for its properties.
            - str: A step-by-step solution showing the analysis and calculations.
    """
    # 1. Parameterize the inputs with random values
    amplitude = random.randint(10, 200)
    
    # Generate a realistic angular frequency (omega)
    omega_multiple = random.randint(2, 9)
    omega = omega_multiple * math.pi * 1e8
    
    # Generate a realistic phase velocity (u_p) by choosing a refractive index
    # This ensures omega and k are physically consistent.
    refractive_index = round(random.uniform(1.0, 2.5), 2)
    u_p = C0 / refractive_index
    
    # Calculate wavenumber (k) based on omega and u_p
    k = omega / u_p
    
    # Randomly determine the direction of propagation and phase
    direction_sign_str = random.choice(['+', '-'])
    phi_deg = random.randint(-180, 180)
    
    # Determine direction text for the solution
    if direction_sign_str == '-':
        direction_text = "the positive z-direction"
    else:
        direction_text = "the negative z-direction"

    # Construct the wave equation string for the question
    # Format omega for better readability in the question
    omega_str = f"{omega_multiple}pi x 10^8"
    wave_equation = (
        f"E(z, t) = {amplitude} * cos({omega_str}*t {direction_sign_str} {round(k, 2)}*z + {phi_deg} deg) V/m"
    )

    # 2. Perform the core calculations for the solution
    frequency = omega / (2 * math.pi)
    wavelength = (2 * math.pi) / k

    # 3. Generate the question and solution strings
    question = (
        f"An electric field wave is described by the following equation:\n"
        f"   {wave_equation}\n\n"
        f"Based on this expression, determine the following properties of the wave:\n"
        f"a) Amplitude (A)\n"
        f"b) Direction of propagation\n"
        f"c) Frequency (f)\n"
        f"d) Wavelength (lambda)\n"
        f"e) Phase velocity (u_p)"
    )

    solution = (
        f"**Given:**\n"
        f"   The wave equation is E(z, t) = {wave_equation}.\n\n"
        
        f"**Step 1:** Compare the equation to the standard form.\n"
        f"   The standard form for a traveling wave is A * cos(omega*t +/- k*z + phi).\n"
        f"   By matching the terms, we can extract the coefficients:\n"
        f"   - Amplitude (A) = {amplitude} V/m\n"
        f"   - Angular Frequency (omega) = {omega_str} rad/s = {omega:.3e} rad/s\n"
        f"   - Wavenumber (k) = {round(k, 2)} rad/m\n"
        f"   - The sign between the t and z terms is '{direction_sign_str}'.\n\n"
        
        f"**Step 2:** Determine the Amplitude and Direction of Propagation.\n"
        f"   The amplitude (A) is the value multiplying the cosine function, which is {amplitude} V/m.\n"
        f"   The direction is determined by the sign in front of the 'k*z' term. A minus sign (-) indicates propagation in the positive direction, while a plus sign (+) indicates the negative direction.\n"
        f"   Therefore, the wave is traveling in {direction_text}.\n\n"
        
        f"**Step 3:** Calculate the Frequency (f).\n"
        f"   Frequency is related to angular frequency by f = omega / (2 * pi).\n"
        f"   f = ({omega:.3e} rad/s) / (2 * pi) = {frequency:.3e} Hz = {frequency/1e6:.2f} MHz.\n\n"
        
        f"**Step 4:** Calculate the Wavelength (lambda).\n"
        f"   Wavelength is related to the wavenumber by lambda = 2 * pi / k.\n"
        f"   lambda = (2 * pi) / ({round(k, 2)} rad/m) = {round(wavelength, 3)} m.\n\n"
        
        f"**Step 5:** Calculate the Phase Velocity (u_p).\n"
        f"   Phase velocity is the speed of the wave, given by u_p = omega / k.\n"
        f"   u_p = ({omega:.3e} rad/s) / ({round(k, 2)} rad/m) = {u_p:.3e} m/s.\n\n"
        
        f"**Answer:**\n"
        f"   - Amplitude: {amplitude} V/m\n"
        f"   - Direction of Propagation: {direction_text}\n"
        f"   - Frequency: {frequency/1e6:.2f} MHz\n"
        f"   - Wavelength: {round(wavelength, 3)} m\n"
        f"   - Phase Velocity: {u_p:.3e} m/s"
    )

    return question, solution


# Template 4 (Intermediate)
def template_phasor_addition():
    """
    Phasor Addition of Sinusoidal Waves

    Scenario:
        This template assesses the ability to add two sinusoidal waves of the same
        frequency. The standard method is to convert both signals into their phasor
        representations, perform complex number addition in rectangular form, convert
        the result back to polar form, and finally express it as a time-domain signal.

    Core Equations:
        V_phasor = A < phi
        x = A * cos(phi)
        y = A * sin(phi)
        A_total = sqrt(x_total^2 + y_total^2)
        phi_total = atan2(y_total, x_total)
        v_total(t) = A_total * cos(omega*t + phi_total)

    Returns:
        tuple: A tuple containing:
            - str: A question asking for the sum of two time-domain signals.
            - str: A step-by-step solution using the phasor addition method.
    """
    # 1. Parameterize the inputs
    precision = 2
    
    # Amplitudes
    A1 = round(random.uniform(10.0, 50.0), precision)
    A2 = round(random.uniform(10.0, 50.0), precision)
    
    # Phases in degrees
    phi1_deg = round(random.uniform(-180.0, 180.0), 1)
    phi2_deg = round(random.uniform(-180.0, 180.0), 1)
    
    # Shared angular frequency
    omega = random.randint(100, 500)
    
    # Function types (cos or sin)
    func_type1 = random.choice(['cos', 'sin'])
    func_type2 = random.choice(['cos', 'sin'])

    # 2. Generate the question string
    v1_str = f"{A1} * {func_type1}({omega}*t + {phi1_deg} deg)"
    v2_str = f"{A2} * {func_type2}({omega}*t + {phi2_deg} deg)"
    
    question = (
        f"Two signals, v1(t) and v2(t), are defined as:\n"
        f"   v1(t) = {v1_str}\n"
        f"   v2(t) = {v2_str}\n\n"
        f"Find the sum, v_total(t) = v1(t) + v2(t), using the phasor method. "
        f"Express the final answer as a single cosine function."
    )

    # 3. Perform the calculations for the solution
    
    # --- Phasor 1 Conversion ---
    phi1_adj_deg = phi1_deg - 90 if func_type1 == 'sin' else phi1_deg
    phi1_rad = math.radians(phi1_adj_deg)
    x1 = A1 * math.cos(phi1_rad)
    y1 = A1 * math.sin(phi1_rad)

    # --- Phasor 2 Conversion ---
    phi2_adj_deg = phi2_deg - 90 if func_type2 == 'sin' else phi2_deg
    phi2_rad = math.radians(phi2_adj_deg)
    x2 = A2 * math.cos(phi2_rad)
    y2 = A2 * math.sin(phi2_rad)

    # --- Addition ---
    x_total = x1 + x2
    y_total = y1 + y2
    
    # --- Convert Sum back to Polar ---
    A_total = math.hypot(x_total, y_total) # sqrt(x^2 + y^2)
    phi_total_rad = math.atan2(y_total, x_total)
    phi_total_deg = math.degrees(phi_total_rad)

    # --- Final time-domain expression ---
    v_total_str = f"{round(A_total, precision)} * cos({omega}*t + {round(phi_total_deg, precision)} deg)"

    # 4. Generate the solution string
    
    # Helper strings for conversion steps
    conv1_step = f"   Since v1(t) is a {func_type1} function, we adjust the phase: {phi1_deg} - 90 = {round(phi1_adj_deg, 1)} deg." if func_type1 == 'sin' else f"   Since v1(t) is a cosine function, no phase adjustment is needed. Phase = {phi1_deg} deg."
    conv2_step = f"   Since v2(t) is a {func_type2} function, we adjust the phase: {phi2_deg} - 90 = {round(phi2_adj_deg, 1)} deg." if func_type2 == 'sin' else f"   Since v2(t) is a cosine function, no phase adjustment is needed. Phase = {phi2_deg} deg."
    
    solution = (
        f"**Given:**\n"
        f"   v1(t) = {v1_str}\n"
        f"   v2(t) = {v2_str}\n\n"
        
        f"**Step 1:** Convert v1(t) to its phasor form V1.\n"
        f"{conv1_step}\n"
        f"   The polar form is V1 = {A1} < {round(phi1_adj_deg, 1)} deg.\n"
        f"   Convert to rectangular form (x1 + jy1):\n"
        f"   x1 = {A1} * cos({round(phi1_adj_deg, 1)}) = {round(x1, precision)}\n"
        f"   y1 = {A1} * sin({round(phi1_adj_deg, 1)}) = {round(y1, precision)}\n"
        f"   So, V1 = {round(x1, precision)} + j{round(y1, precision)}.\n\n"
        
        f"**Step 2:** Convert v2(t) to its phasor form V2.\n"
        f"{conv2_step}\n"
        f"   The polar form is V2 = {A2} < {round(phi2_adj_deg, 1)} deg.\n"
        f"   Convert to rectangular form (x2 + jy2):\n"
        f"   x2 = {A2} * cos({round(phi2_adj_deg, 1)}) = {round(x2, precision)}\n"
        f"   y2 = {A2} * sin({round(phi2_adj_deg, 1)}) = {round(y2, precision)}\n"
        f"   So, V2 = {round(x2, precision)} + j{round(y2, precision)}.\n\n"
        
        f"**Step 3:** Add the phasors in rectangular form: V_total = V1 + V2.\n"
        f"   V_total = ({round(x1, precision)} + j{round(y1, precision)}) + ({round(x2, precision)} + j{round(y2, precision)})\n"
        f"   V_total = ({round(x1, precision)} + {round(x2, precision)}) + j({round(y1, precision)} + {round(y2, precision)})\n"
        f"   V_total = {round(x_total, precision)} + j{round(y_total, precision)}.\n\n"

        f"**Step 4:** Convert the resultant phasor V_total back to polar form (A < phi).\n"
        f"   A_total = sqrt({round(x_total, precision)}^2 + {round(y_total, precision)}^2) = {round(A_total, precision)}\n"
        f"   phi_total = atan2({round(y_total, precision)}, {round(x_total, precision)}) = {round(phi_total_deg, precision)} deg\n"
        f"   So, V_total = {round(A_total, precision)} < {round(phi_total_deg, precision)} deg.\n\n"

        f"**Step 5:** Convert the final phasor back to the time domain.\n"
        f"   The resulting signal has amplitude A_total, phase phi_total, and the original angular frequency omega.\n"
        f"   v_total(t) = {v_total_str}\n\n"

        f"**Answer:**\n"
        f"   The sum of the two signals is v_total(t) = {v_total_str}."
    )
    
    return question, solution


# Template 5 (Advanced)
def template_standing_wave_formation():
    """
    Standing Wave Formation and Properties

    Scenario:
        This template models the superposition of two identical waves traveling in
        opposite directions, which creates a standing wave. It requires using a
        trigonometric identity to derive the mathematical form of the standing wave
        and then interpreting this form to find the locations of nodes (zero
        amplitude) and antinodes (maximum amplitude).

    Core Equations:
        cos(a) + cos(b) = 2 * cos((a-b)/2) * cos((a+b)/2)
        Standing Wave: y(x,t) = [2*A*cos(k*x)] * cos(omega*t)
        Nodes: k*x = (n + 1/2)*pi
        Antinodes: k*x = n*pi

    Returns:
        tuple: A tuple containing:
            - str: A question asking to derive and analyze a standing wave.
            - str: A step-by-step solution showing the derivation and calculations.
    """
    # 1. Parameterize the inputs
    precision = 2
    A = random.randint(5, 50)
    omega = random.randint(100, 500)
    k = round(random.uniform(1.0, 5.0), precision)

    # 2. Generate the question string
    y1_str = f"{A} * cos({omega}*t - {k}*x)"
    y2_str = f"{A} * cos({omega}*t + {k}*x)"
    
    question = (
        f"Two waves of the same amplitude, frequency, and wavelength travel in opposite directions, given by:\n"
        f"   y1(x,t) = {y1_str}\n"
        f"   y2(x,t) = {y2_str}\n\n"
        f"a) Find the superposition of these waves, y(x,t) = y1(x,t) + y2(x,t), and show that it represents a standing wave.\n"
        f"b) Determine the locations of the first three nodes and antinodes for x >= 0."
    )

    # 3. Perform the calculations for the solution
    
    # Standing wave amplitude
    standing_wave_amplitude = 2 * A
    
    # Wavelength
    wavelength = (2 * math.pi) / k
    
    # Node locations
    nodes = [(n + 0.5) * math.pi / k for n in range(3)]
    
    # Antinode locations
    antinodes = [n * math.pi / k for n in range(3)]

    # 4. Generate the solution string
    
    standing_wave_eq = f"{standing_wave_amplitude} * cos({k}*x) * cos({omega}*t)"
    
    solution = (
        f"**Given:**\n"
        f"   Incident wave y1(x,t) = {y1_str}\n"
        f"   Reflected wave y2(x,t) = {y2_str}\n\n"
        
        f"**Step 1:** Sum the two traveling waves using a trigonometric identity.\n"
        f"   We use the identity: cos(alpha) + cos(beta) = 2 * cos((alpha - beta)/2) * cos((alpha + beta)/2).\n"
        f"   Let alpha = {omega}*t - {k}*x and beta = {omega}*t + {k}*x.\n"
        f"   (alpha - beta)/2 = (({omega}*t - {k}*x) - ({omega}*t + {k}*x))/2 = -{k}*x\n"
        f"   (alpha + beta)/2 = (({omega}*t - {k}*x) + ({omega}*t + {k}*x))/2 = {omega}*t\n"
        f"   Substituting these into the identity and using cos(-z) = cos(z), we get:\n"
        f"   y(x,t) = {A} * [2 * cos(-{k}*x) * cos({omega}*t)] = {standing_wave_eq}\n"
        f"   This is the equation of a standing wave because the spatial dependence (cos({k}*x)) is separate from the time dependence (cos({omega}*t)).\n\n"

        f"**Step 2:** Find the locations of the nodes.\n"
        f"   Nodes are points of zero amplitude, which occur when the spatial term is zero: cos({k}*x) = 0.\n"
        f"   This condition is met when k*x = (n + 1/2)*pi, for n = 0, 1, 2, ...\n"
        f"   Solving for x: x = (n + 1/2) * pi / k.\n"
        f"   - For n=0: x = (0.5 * pi) / {k} = {round(nodes[0], precision)}\n"
        f"   - For n=1: x = (1.5 * pi) / {k} = {round(nodes[1], precision)}\n"
        f"   - For n=2: x = (2.5 * pi) / {k} = {round(nodes[2], precision)}\n\n"

        f"**Step 3:** Find the locations of the antinodes.\n"
        f"   Antinodes are points of maximum amplitude, which occur when |cos({k}*x)| = 1.\n"
        f"   This condition is met when k*x = n*pi, for n = 0, 1, 2, ...\n"
        f"   Solving for x: x = n * pi / k.\n"
        f"   - For n=0: x = (0 * pi) / {k} = {round(antinodes[0], precision)}\n"
        f"   - For n=1: x = (1 * pi) / {k} = {round(antinodes[1], precision)}\n"
        f"   - For n=2: x = (2 * pi) / {k} = {round(antinodes[2], precision)}\n\n"

        f"**Answer:**\n"
        f"   - The standing wave equation is y(x,t) = {standing_wave_eq}.\n"
        f"   - The first three nodes are at x = {round(nodes[0], precision)}, {round(nodes[1], precision)}, and {round(nodes[2], precision)}.\n"
        f"   - The first three antinodes are at x = {round(antinodes[0], precision)}, {round(antinodes[1], precision)}, and {round(antinodes[2], precision)}."
    )
    
    return question, solution


def main():
    """
    Generate numerous instances of each waves and phasors template 
    with different random seeds and write the results to a JSONL file.
    """
    import json
    import os

    # Define the output path (Modify this path according to where you are running the code from)
    output_file = "testset/electrical_engineering/electromagnetics_and_waves/waves_and_phasors.jsonl"

    # Create the directory if it doesn't exist
    os.makedirs(os.path.dirname(output_file), exist_ok=True)

    # List of template functions with their ID and level
    templates = [
        (template_wave_parameters_basic, "wave_parameters_basic", "Easy"),
        (template_time_to_phasor, "time_to_phasor", "Easy"),
        (template_wave_equation_interpretation, "wave_equation_interpretation", "Intermediate"),
        (template_phasor_addition, "phasor_addition", "Intermediate"),
        (template_standing_wave_formation, "standing_wave_formation", "Advanced"),
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
                "branch": "electrical_engineering",
                "domain": "electromagnetics_and_waves",
                "area": "waves_and_phasors",
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
