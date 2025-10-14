import random
import math


# Template 1 (Easy)
def template_signal_energy_power():
    """
    Signal Energy and Power Calculation

    Scenario:
        This template tests the ability to classify a signal as either an energy
        signal (finite energy) or a power signal (finite average power) and to
        calculate the corresponding value. It covers time-limited, time-unlimited
        decaying, and periodic signals.

    Core Equations:
        Signal Energy: Eg = integral from -inf to inf of |g(t)|^2 dt
        Signal Power: Pg = lim(T->inf) (1/T) * integral from -T/2 to T/2 of |g(t)|^2 dt

    Returns:
        tuple: A tuple containing:
            - str: A question asking to classify and calculate the energy/power of a signal.
            - str: A step-by-step solution.
    """
    # 1. Parameterize the inputs with random values
    signal_type = random.choice(['rect_pulse', 'exp_decay', 'periodic'])
    amplitude = random.randint(2, 20)
    precision = 2

    # Common text for all questions to ensure unit clarity
    unit_context = "Assume the signal g(t) represents a voltage in Volts (V) across a 1-ohm resistor."

    if signal_type == 'rect_pulse':
        # --- Energy Signal: Rectangular Pulse ---
        duration = random.randint(2, 10)
        start_time = 0
        end_time = start_time + duration
        signal_expr = f"g(t) = {amplitude} for {start_time} <= t <= {end_time}, and 0 otherwise"
        energy = (amplitude ** 2) * duration

        question = (
            f"Consider the signal defined as:\n"
            f"{signal_expr}\n\n"
            f"{unit_context}\n\n"
            f"Is this an energy signal or a power signal? Calculate its total energy or "
            f"average power accordingly."
        )

        solution = (
            f"**Given:**\n"
            f"The signal is {signal_expr}.\n\n"

            f"**Step 1:** Classify the signal.\n"
            f"The signal is non-zero only for a finite duration ({duration} seconds). "
            f"Time-limited signals have finite energy and zero average power. "
            f"Therefore, this is an **energy signal**.\n\n"

            f"**Step 2:** State the formula for signal energy.\n"
            f"The energy (Eg) of a signal g(t) is:\n"
            f"Eg = integral from -inf to inf of |g(t)|^2 dt\n\n"

            f"**Step 3:** Apply the formula to the given signal.\n"
            f"Since the signal is non-zero only between t = {start_time} and t = {end_time}, the integral becomes:\n"
            f"Eg = integral from {start_time} to {end_time} of ({amplitude})^2 dt = integral from {start_time} to {end_time} of {amplitude**2} dt\n\n"

            f"**Step 4:** Calculate the integral.\n"
            f"Eg = {amplitude**2} * [t] (from t={start_time} to {end_time})\n"
            f"Eg = {amplitude**2} * ({end_time} - {start_time}) = {amplitude**2} * {duration} = {energy}\n\n"

            f"**Answer:**\n"
            f"The signal is an **energy signal** with a total energy of **{round(energy, precision)} Joules (V^2 * s)**."
        )

    elif signal_type == 'exp_decay':
        # --- Energy Signal: Exponential Decay ---
        alpha = random.randint(2, 10)
        signal_expr = f"g(t) = {amplitude} * exp(-{alpha}*t) * u(t)"
        energy = (amplitude ** 2) / (2 * alpha)

        question = (
            f"Consider the signal defined as:\n"
            f"{signal_expr}\n"
            f"where u(t) is the unit step function.\n\n"
            f"{unit_context}\n\n"
            f"Is this an energy signal or a power signal? Calculate its total energy or "
            f"average power accordingly."
        )

        solution = (
            f"**Given:**\n"
            f"The signal is {signal_expr}.\n\n"

            f"**Step 1:** Classify the signal.\n"
            f"The signal is not time-limited, but it decays to zero as t approaches infinity. "
            f"This means its total energy is finite. Therefore, this is an **energy signal**.\n\n"

            f"**Step 2:** State the formula for signal energy.\n"
            f"Eg = integral from -inf to inf of |g(t)|^2 dt\n\n"

            f"**Step 3:** Apply the formula to the given signal.\n"
            f"Due to the unit step function u(t), the integral's limits become 0 to infinity:\n"
            f"Eg = integral from 0 to inf of ({amplitude} * exp(-{alpha}*t))^2 dt\n"
            f"Eg = {amplitude**2} * integral from 0 to inf of exp(-{2*alpha}*t) dt\n\n"

            f"**Step 4:** Calculate the integral.\n"
            f"Eg = {amplitude**2} * [-1/{2*alpha} * exp(-{2*alpha}*t)] (from t=0 to inf)\n"
            f"Eg = {amplitude**2} * ( (0) - (-1/{2*alpha} * exp(0)) ) = {amplitude**2} * (1/{2*alpha})\n"
            f"Eg = {amplitude**2} / {2*alpha} = {energy}\n\n"

            f"**Answer:**\n"
            f"The signal is an **energy signal** with a total energy of **{round(energy, precision)} Joules (V^2 * s)**."
        )

    else: # signal_type == 'periodic'
        # --- Power Signal: Periodic Sinusoid ---
        frequency = random.randint(10, 100)
        signal_expr = f"g(t) = {amplitude} * cos(2 * pi * {frequency} * t)"
        power = (amplitude ** 2) / 2

        question = (
            f"Consider the signal defined as:\n"
            f"{signal_expr}\n\n"
            f"{unit_context}\n\n"
            f"Is this an energy signal or a power signal? Calculate its total energy or "
            f"average power accordingly."
        )

        solution = (
            f"**Given:**\n"
            f"The signal is {signal_expr}.\n\n"

            f"**Step 1:** Classify the signal.\n"
            f"The signal is periodic and continues for all time without decaying. Its total energy would be infinite. "
            f"However, its average power over time is finite and non-zero. "
            f"Therefore, this is a **power signal**.\n\n"

            f"**Step 2:** State the formula for signal power.\n"
            f"For a periodic signal with period T0, the average power (Pg) is calculated over one period:\n"
            f"Pg = (1/T0) * integral over T0 of |g(t)|^2 dt\n\n"

            f"**Step 3:** Apply the formula to the given signal.\n"
            f"For any sinusoidal signal of the form A*cos(w*t + phi), the average power is a standard result:\n"
            f"Pg = A^2 / 2\n"
            f"This result is derived from integrating (A * cos(...))^2 over one full period.\n\n"

            f"**Step 4:** Calculate the power.\n"
            f"In this case, the amplitude A is {amplitude}.\n"
            f"Pg = ({amplitude})^2 / 2\n"
            f"Pg = {amplitude**2} / 2 = {power}\n\n"

            f"**Answer:**\n"
            f"The signal is a **power signal** with an average power of **{round(power, precision)} Watts (V^2)**."
        )
        
    return question, solution


# Template 2 (Easy)
def template_mean_variance():
    """
    Mean and Variance of a Random Variable

    Scenario:
        This template assesses the understanding of basic statistical properties of a
        discrete random variable, specifically its central tendency (mean) and
        dispersion (variance).

    Core Equations:
        Mean (Expected Value): mu_X = E[X] = sum(x_i * P(x_i))
        Variance: sigma_X^2 = E[(X - mu_X)^2] = sum((x_i - mu_X)^2 * P(x_i))

    Returns:
        tuple: A tuple containing:
            - str: A question asking for the mean and variance of a discrete random variable.
            - str: A step-by-step solution.
    """
    # 1. Parameterize the inputs with random values
    num_values = random.randint(3, 5)
    precision = 3

    # Generate a set of unique, sorted integer values for the random variable
    values = sorted(random.sample(range(-10, 21), num_values))

    # Generate clean, rational probabilities that sum to 1
    # Create random integer parts, sum them, and use that sum as the denominator
    parts = [random.randint(1, 10) for _ in range(num_values)]
    total_sum = sum(parts)
    probabilities = [p / total_sum for p in parts]

    # 2. Perform the core calculations
    # Calculate the mean (expected value)
    mean = sum(v * p for v, p in zip(values, probabilities))

    # Calculate the variance
    variance = sum(((v - mean) ** 2) * p for v, p in zip(values, probabilities))

    # 3. Generate the question and solution strings
    
    # Format values and probabilities for display
    values_str = "{" + ", ".join(map(str, values)) + "}"
    probs_str = "{" + ", ".join([f"{p:.{precision}f}" for p in probabilities]) + "}"

    question = (
        f"A discrete random variable X can take the values {values_str} with "
        f"corresponding probabilities {probs_str}, respectively.\n\n"
        f"Calculate the mean (mu_X) and variance (sigma_X^2) of X."
    )

    # --- Build the solution string step-by-step ---

    # Mean calculation steps
    mean_calc_lhs = "mu_X = E[X]"
    mean_calc_rhs = " + ".join([f"({v})({p:.{precision}f})" for v, p in zip(values, probabilities)])
    mean_interm_vals = " + ".join([f"{v * p:.{precision}f}" for v, p in zip(values, probabilities)])
    
    # Variance calculation steps
    var_calc_lhs = "sigma_X^2 = E[(X - mu_X)^2]"
    var_calc_rhs = " + ".join(
        [f"({v} - {mean:.{precision}f})^2({p:.{precision}f})" for v, p in zip(values, probabilities)]
    )
    var_interm_vals = " + ".join(
        [f"({(v - mean):.{precision}f})^2({p:.{precision}f})" for v, p in zip(values, probabilities)]
    )
    var_final_vals = " + ".join(
        [f"{((v - mean)**2 * p):.{precision}f}" for v, p in zip(values, probabilities)]
    )

    solution = (
        f"**Given:**\n"
        f"Values: X = {values_str}\n"
        f"Probabilities: P(X) = {probs_str}\n\n"
        
        f"**Step 1:** Calculate the Mean (mu_X).\n"
        f"The formula for the mean (expected value) is:\n"
        f"mu_X = sum(x_i * P(x_i))\n\n"
        f"Plugging in the given values:\n"
        f"{mean_calc_lhs} = {mean_calc_rhs}\n"
        f"{mean_calc_lhs} = {mean_interm_vals}\n"
        f"{mean_calc_lhs} = {mean:.{precision}f}\n\n"

        f"**Step 2:** Calculate the Variance (sigma_X^2).\n"
        f"The formula for variance is:\n"
        f"sigma_X^2 = sum((x_i - mu_X)^2 * P(x_i))\n\n"
        f"Using the calculated mean (mu_X = {mean:.{precision}f}):\n"
        f"{var_calc_lhs} = {var_calc_rhs}\n"
        f"{var_calc_lhs} = {var_interm_vals}\n"
        f"{var_calc_lhs} = {var_final_vals}\n"
        f"{var_calc_lhs} = {variance:.{precision}f}\n\n"

        f"**Answer:**\n"
        f"The mean of the random variable X is **{mean:.{precision}f}**.\n"
        f"The variance of the random variable X is **{variance:.{precision}f}**."
    )

    return question, solution


# Template 3 (Intermediate)
def template_autocorrelation_rect_pulse():
    """
    Autocorrelation of a Deterministic Signal

    Scenario:
        This template tests the ability to calculate the autocorrelation function of
        a simple deterministic signal, a rectangular pulse. This measures the
        similarity of the signal with a time-shifted version of itself and
        introduces the concept of deriving a new function shape (a triangle) from
        this operation.

    Core Equation:
        Autocorrelation: R_g(tau) = integral from -inf to inf of g(t) * g(t - tau) dt

    Returns:
        tuple: A tuple containing:
            - str: A question asking for the autocorrelation of a rectangular pulse.
            - str: A step-by-step solution.
    """
    # 1. Parameterize the inputs with random values
    amplitude = random.randint(2, 10)
    
    # Use a half-duration to ensure the total duration is an even number
    # and the pulse limits are clean integers.
    half_duration = random.randint(1, 5)
    duration = 2 * half_duration

    # 2. Perform the core calculation
    # The autocorrelation of a rectangular pulse is a triangular function.
    # Peak value R(0) = A^2 * T
    peak_value = (amplitude ** 2) * duration
    
    # The resulting function is R(tau) = A^2 * (T - |tau|) for |tau| <= T
    result_function_str = f"{amplitude**2}*({duration} - |tau|)"

    # 3. Generate the question and solution strings
    signal_expr = f"g(t) = {amplitude} for -{half_duration} <= t <= {half_duration}, and 0 otherwise"

    question = (
        f"Find the autocorrelation function R_g(tau) for the rectangular pulse signal defined as:\n"
        f"    {signal_expr}"
    )

    solution = (
        f"**Given:**\n"
        f"The signal is a rectangular pulse: {signal_expr}.\n"
        f"It has an amplitude A = {amplitude} and a total duration T = {duration}.\n\n"
        
        f"**Step 1:** Understand the Autocorrelation Integral.\n"
        f"The autocorrelation function is defined as R_g(tau) = integral of g(t) * g(t - tau) dt.\n"
        f"Conceptually, this measures the overlapping area between the signal g(t) and a version of itself shifted by 'tau'.\n\n"

        f"**Step 2:** Analyze the Overlap of the Pulses.\n"
        f"The original pulse g(t) exists from t = -{half_duration} to t = {half_duration}.\n"
        f"The shifted pulse g(t - tau) exists from t = tau - {half_duration} to t = tau + {half_duration}.\n"
        f"The product g(t) * g(t - tau) is non-zero only where these two pulses overlap.\n"
        f"In the overlapping region, the value of the product is ({amplitude}) * ({amplitude}) = {amplitude**2}.\n"
        f"The pulses only overlap when the shift 'tau' is between -{duration} and {duration}. Outside this range, the autocorrelation is zero.\n\n"

        f"**Step 3:** Determine the Area of the Overlap.\n"
        f"The integral is simply the area of the rectangular overlapping region.\n"
        f"Area = (Height of Overlap) * (Width of Overlap)\n\n"
        f"The height is constant at {amplitude**2}.\n"
        f"The width (duration) of the overlap depends on the shift 'tau'. When tau = 0, the overlap is the entire pulse duration, {duration}. As |tau| increases, the overlap width decreases linearly. The width is given by the expression: ({duration} - |tau|).\n\n"
        
        f"**Step 4:** Formulate the Autocorrelation Function.\n"
        f"Combining the height and width, the area of overlap is:\n"
        f"R_g(tau) = {amplitude**2} * ({duration} - |tau|)\n\n"
        f"This equation is valid for the range where overlap occurs, which is |tau| <= {duration}. For |tau| > {duration}, R_g(tau) = 0.\n\n"

        f"**Answer:**\n"
        f"The autocorrelation function is a **triangular function** described by:\n"
        f"**R_g(tau) = {result_function_str}**, for **|tau| <= {duration}**, and **0** otherwise.\n"
        f"This triangle has a peak value of **{peak_value}** at tau = 0."
    )

    return question, solution


# Template 4 (Intermediate)
def template_ft_esd_rect_pulse():
    """
    Fourier Transform and Energy Spectral Density (ESD)

    Scenario:
        This template requires finding the Fourier Transform of a simple energy
        signal (a rectangular pulse) and then calculating its Energy Spectral
        Density (ESD). It demonstrates the fundamental relationship between a
        signal's time-domain representation and its frequency-domain energy distribution.

    Core Equations:
        Fourier Transform: G(f) = integral from -inf to inf of g(t)*exp(-j*2*pi*f*t) dt
        Energy Spectral Density (ESD): Psi_g(f) = |G(f)|^2

    Returns:
        tuple: A tuple containing:
            - str: A question asking for the FT and ESD of a rectangular pulse.
            - str: A step-by-step solution.
    """
    # 1. Parameterize the inputs with random values
    amplitude = random.randint(5, 20)
    duration = random.choice([0.5, 1.0, 1.5, 2.0, 4.0])
    half_duration = duration / 2
    
    # 2. Perform the core calculation and format results
    ft_amplitude = amplitude * duration
    esd_amplitude = ft_amplitude ** 2
    
    # Handle sinc argument formatting for cleaner output
    if duration == 1.0:
        sinc_arg_str = "f"
    else:
        sinc_arg_str = f"{duration}*f"

    ft_result_str = f"{ft_amplitude} * sinc({sinc_arg_str})"
    esd_result_str = f"{esd_amplitude} * sinc^2({sinc_arg_str})"

    # 3. Generate the question and solution strings
    signal_expr = f"g(t) with amplitude {amplitude} extending from t = -{half_duration} s to t = {half_duration} s"

    question = (
        f"Given a rectangular pulse {signal_expr}, find its Fourier Transform G(f) and "
        f"its Energy Spectral Density Psi_g(f).\n\n"
        f"Note: Use the definition sinc(x) = sin(pi*x) / (pi*x)."
    )

    solution = (
        f"**Given:**\n"
        f"A rectangular pulse g(t) with:\n"
        f"Amplitude (A) = {amplitude}\n"
        f"Total Duration (T) = {duration} s (from -{half_duration} to {half_duration})\n\n"

        f"**Step 1:** Identify the Fourier Transform Pair for a Rectangular Pulse.\n"
        f"A rectangular pulse centered at the origin, g(t) = A for |t| <= T/2, has a well-known Fourier Transform which is a sinc function:\n"
        f"G(f) = A * T * sinc(f * T)\n"
        f"where sinc(x) is defined as sin(pi*x) / (pi*x).\n\n"

        f"**Step 2:** Calculate the Fourier Transform G(f).\n"
        f"We plug in our given values A = {amplitude} and T = {duration} into the formula:\n"
        f"G(f) = ({amplitude}) * ({duration}) * sinc(f * {duration})\n"
        f"G(f) = {ft_result_str}\n\n"

        f"**Step 3:** State the Formula for Energy Spectral Density (ESD).\n"
        f"The ESD, denoted Psi_g(f), describes how the energy of the signal is distributed over frequency. It is calculated as the squared magnitude of the Fourier Transform:\n"
        f"Psi_g(f) = |G(f)|^2\n\n"

        f"**Step 4:** Calculate the Energy Spectral Density Psi_g(f).\n"
        f"We take the magnitude squared of the G(f) we found in Step 2. Since our G(f) is real-valued, this is equivalent to just squaring the function:\n"
        f"Psi_g(f) = |{ft_result_str}|^2\n"
        f"Psi_g(f) = ({amplitude} * {duration})^2 * sinc^2({sinc_arg_str})\n"
        f"Psi_g(f) = {esd_result_str}\n\n"

        f"**Answer:**\n"
        f"The Fourier Transform is **G(f) = {ft_result_str}**.\n"
        f"The Energy Spectral Density is **Psi_g(f) = {esd_result_str}**."
    )

    return question, solution


# Template 5 (Advanced)
def template_lowpass_equivalent_bandpass():
    """
    Lowpass Equivalent of a Bandpass Signal

    Scenario:
        A fundamental concept in communications is representing a high-frequency bandpass
        signal using a lower-frequency complex equivalent. This template tests the
        ability to derive this lowpass equivalent signal from a given bandpass
        signal by expanding it into its canonical in-phase and quadrature components.

    Core Equations:
        Bandpass Signal: g(t) = g_I(t)*cos(2*pi*fc*t) - g_Q(t)*sin(2*pi*fc*t)
        Lowpass Equivalent Signal: g_tilde(t) = g_I(t) + j*g_Q(t)
        Trigonometric Identity: cos(a + b) = cos(a)*cos(b) - sin(a)*sin(b)

    Returns:
        tuple: A tuple containing:
            - str: A question asking for the components of a bandpass signal.
            - str: A step-by-step solution showing the derivation.
    """
    # 1. Parameterize the inputs with random values
    amplitude = random.randint(10, 100)
    fc_exp = random.randint(6, 8)  # For MHz to hundreds of MHz
    fc_mult = random.randint(1, 9)
    phase_deg = random.choice([30, 45, 60, 90, 120, 135, 150])
    precision = 2
    
    # Create a readable string for the carrier frequency
    fc_str = f"{fc_mult} x 10^{fc_exp}"
    
    # 2. Perform the core calculations
    phase_rad = math.radians(phase_deg)
    cos_phi = math.cos(phase_rad)
    sin_phi = math.sin(phase_rad)

    # From the identity, g(t) = A*cos(phi)*cos(w_c*t) - A*sin(phi)*sin(w_c*t)
    # We can identify the in-phase and quadrature components
    g_I = amplitude * cos_phi
    g_Q = amplitude * sin_phi

    # 3. Generate the question and solution strings
    signal_expr = f"g(t) = {amplitude} * cos(2*pi*({fc_str})*t + {phase_deg} deg)"
    
    question = (
        f"A bandpass signal is described by the following equation:\n"
        f"{signal_expr}\n\n"
        f"Express this signal in its canonical form to determine its in-phase component (g_I(t)), "
        f"quadrature component (g_Q(t)), and its complex lowpass equivalent signal (g_tilde(t))."
    )
    
    solution = (
        f"**Given:**\n"
        f"The bandpass signal is {signal_expr}.\n\n"
        
        f"**Step 1:** State the Goal and the Required Identity.\n"
        f"Our goal is to match the given signal to the canonical bandpass form:\n"
        f"g(t) = g_I(t)*cos(2*pi*fc*t) - g_Q(t)*sin(2*pi*fc*t)\n"
        f"To do this, we use the trigonometric angle addition identity:\n"
        f"cos(a + b) = cos(a)*cos(b) - sin(a)*sin(b)\n\n"

        f"**Step 2:** Expand the Given Signal Expression.\n"
        f"Let a = 2*pi*fc*t and b = {phase_deg} deg. Applying the identity to g(t):\n"
        f"g(t) = {amplitude} * [cos(2*pi*fc*t)*cos({phase_deg} deg) - sin(2*pi*fc*t)*sin({phase_deg} deg)]\n"
        f"g(t) = ({amplitude}*cos({phase_deg} deg))*cos(2*pi*fc*t) - ({amplitude}*sin({phase_deg} deg))*sin(2*pi*fc*t)\n\n"

        f"**Step 3:** Identify the In-Phase (g_I) and Quadrature (g_Q) Components.\n"
        f"By comparing our expanded signal to the canonical form, we can identify g_I(t) and g_Q(t):\n"
        f"g_I(t) is the term multiplying cos(2*pi*fc*t).\n"
        f"g_Q(t) is the term multiplying sin(2*pi*fc*t) (note the minus sign is part of the formula).\n\n"
        f"g_I(t) = {amplitude} * cos({phase_deg} deg) = {amplitude} * ({cos_phi:.{precision+1}f}) = {g_I:.{precision+1}f}\n"
        f"g_Q(t) = {amplitude} * sin({phase_deg} deg) = {amplitude} * ({sin_phi:.{precision+1}f}) = {g_Q:.{precision+1}f}\n"
        f"Note that for this signal, g_I and g_Q are constants because the amplitude and phase are constant.\n\n"
        
        f"**Step 4:** Construct the Complex Lowpass Equivalent Signal (g_tilde(t)).\n"
        f"The complex lowpass equivalent is defined as g_tilde(t) = g_I(t) + j*g_Q(t).\n"
        f"g_tilde(t) = {round(g_I, precision)} + j*{round(g_Q, precision)}\n\n"

        f"**Answer:**\n"
        f"In-Phase Component: **g_I(t) = {round(g_I, precision)}**\n"
        f"Quadrature Component: **g_Q(t) = {round(g_Q, precision)}**\n"
        f"Complex Lowpass Equivalent: **g_tilde(t) = {round(g_I, precision)} + j*{round(g_Q, precision)}**"
    )

    return question, solution


def main():
    """
    Generate numerous instances of each deterministic and random signal analysis template 
    with different random seeds and write the results to a JSONL file.
    """
    import json
    import os

    # Define the output path (Modify this path according to where you are running the code from)
    output_file = "testset/electrical_engineering/digital_communications/deterministic_and_random_signal_analysis.jsonl"

    # Create the directory if it doesn't exist
    os.makedirs(os.path.dirname(output_file), exist_ok=True)

    # List of template functions with their ID and level
    templates = [
        (template_signal_energy_power, "signal_energy_power", "Easy"),
        (template_mean_variance, "mean_variance_discrete_rv", "Easy"),
        (template_autocorrelation_rect_pulse, "autocorrelation_rect_pulse", "Intermediate"),
        (template_ft_esd_rect_pulse, "ft_esd_rect_pulse", "Intermediate"),
        (template_lowpass_equivalent_bandpass, "lowpass_equivalent_bandpass", "Advanced"),
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
                "domain": "digital_communications",
                "area": "deterministic_and_random_signal_analysis",
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
