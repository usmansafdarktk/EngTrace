import random
import math
from fractions import Fraction
from data.templates.branches.electrical_engineering.constants import FREQUENCY_RANGE_HZ, AMPLITUDE_RANGE, PHASE_RANGE_DEG, PHASE_RANGE_RAD, OMEGA_MULTIPLIER_RANGE, SAMPLING_FREQ_RANGE_HZ, F0_RANGE_HZ, GAIN_K_RANGE, DELAY_N0_RANGE, DECIMATION_FACTOR_M_RANGE, OMEGA_DENOMINATOR_RANGE


# Template 1 (Easy)
def template_nyquist_rate_determination():
    """
    Nyquist Rate Determination

    Scenario:
        This template tests the fundamental understanding of the Nyquist-Shannon
        sampling theorem. The user must identify the highest frequency component in a
        continuous-time signal and use it to calculate the minimum sampling rate
        required to prevent aliasing.

    Core Equations:
        F_nyquist = 2 * F_max

    Returns:
        tuple: A tuple containing:
            - str: A question asking for the Nyquist rate of a given signal.
            - str: A step-by-step solution.
    """
    # 1. Parameterize the inputs with random values

    # Determine the number of sinusoidal components in the signal
    num_components = random.choice([2, 3])

    # Generate a set of unique, random frequencies to serve as the basis for the problem.
    # Using random.sample ensures all generated frequencies are distinct.
    frequency_pool = range(FREQUENCY_RANGE_HZ[0], FREQUENCY_RANGE_HZ[1])
    frequencies_hz = sorted(random.sample(frequency_pool, num_components), reverse=True)

    signal_terms = []
    for freq in frequencies_hz:
        # For each frequency, generate a random amplitude, phase, and function type (sin/cos)
        amplitude = round(random.uniform(*AMPLITUDE_RANGE), 1)
        phase_deg = random.randint(*PHASE_RANGE_DEG)
        func_type = random.choice(['cos', 'sin'])

        # Create the mathematical term as a string
        # Example: "15.5 * cos(2*pi*750*t + 45 deg)"
        term_str = f"{amplitude} * {func_type}(2*pi*{freq}*t + {phase_deg} deg)"
        signal_terms.append(term_str)

    # Join the individual terms with a " + " to form the full signal equation
    signal_expression = " + ".join(signal_terms)

    # 2. Perform the core calculation for the solution
    # The highest frequency determines the Nyquist rate.
    f_max = max(frequencies_hz)

    # The Nyquist rate is twice the maximum frequency.
    f_nyquist = 2 * f_max

    # 3. Generate the question and solution strings
    question = (
        f"A continuous-time signal is defined as:\n"
        f"  x_c(t) = {signal_expression}\n\n"
        f"What is the minimum sampling rate (Fs), also known as the Nyquist rate, "
        f"required to sample this signal without loss of information?"
    )

    solution = (
        f"**Given Signal:**\n"
        f"x_c(t) = {signal_expression}\n\n"

        f"**Step 1:** Identify Frequencies\n"
        f"First, we identify all the unique frequency components in the signal's expression. "
        f"By inspecting the '2*pi*F*t' part of each term, we find the frequencies are:\n"
        f"Frequencies = {', '.join(map(str, frequencies_hz))} Hz.\n\n"

        f"**Step 2:** Find Maximum Frequency\n"
        f"The Nyquist rate is determined by the highest frequency component in the signal. "
        f"We find the maximum value from the list of frequencies.\n"
        f"F_max = max({', '.join(map(str, frequencies_hz))}) = {f_max} Hz.\n\n"

        f"**Step 3:** Apply Nyquist Theorem\n"
        f"The Nyquist-Shannon sampling theorem states that to avoid aliasing and be able "
        f"to perfectly reconstruct a signal, the sampling rate (Fs) must be at least "
        f"twice its highest frequency component (F_max).\n"
        f"Fs_min = 2 * F_max.\n\n"

        f"**Step 4:** Calculate the Nyquist Rate\n"
        f"Using the maximum frequency found in Step 2, we calculate the Nyquist rate:\n"
        f"FNyquist = 2 * {f_max} = {f_nyquist} Hz.\n\n"

        f"**Answer:**\n"
        f"The minimum sampling rate required to sample the signal without loss of "
        f"information is {f_nyquist} Hz."
    )

    return question, solution


# Template 2 (Easy)
def template_continuous_to_discrete_conversion():
    """
    Continuous-to-Discrete Signal Conversion

    Scenario:
        This template tests the direct application of the sampling process. The user is
        given a continuous-time sinusoid and a sampling rate and must derive the
        corresponding discrete-time sequence by substituting t=nT. It reinforces the
        relationship between continuous and discrete frequencies (omega = Omega * T).

    Core Equations:
        x[n] = x_c(nT)
        T = 1 / Fs
        omega = Omega * T

    Returns:
        tuple: A tuple containing:
            - str: A question asking for the discrete-time representation of a signal.
            - str: A step-by-step solution.
    """
    # 1. Parameterize the inputs with random values
    amplitude = round(random.uniform(*AMPLITUDE_RANGE), 2)

    # Generate a continuous-time angular frequency as an integer multiple of pi
    omega_multiplier = random.randint(*OMEGA_MULTIPLIER_RANGE)
    omega_continuous = omega_multiplier * math.pi
    omega_continuous_str = f"{omega_multiplier}*pi"

    # Randomly choose to express the phase in degrees or radians for variety
    use_degrees = random.choice([True, False])
    if use_degrees:
        phi_deg = random.randint(*PHASE_RANGE_DEG)
        phi_str = f"{phi_deg} deg"
    else:
        phi_rad = round(random.uniform(*PHASE_RANGE_RAD), 2)
        phi_str = f"{phi_rad} rad"

    # Generate a random sampling frequency
    sampling_freq_hz = random.randint(*SAMPLING_FREQ_RANGE_HZ)

    # Construct the full continuous-time signal expression for the question
    signal_expression = f"{amplitude} * cos({omega_continuous_str}*t + {phi_str})"

    # 2. Perform the core calculation for the solution
    sampling_period = 1 / sampling_freq_hz

    # The discrete-time angular frequency
    omega_discrete = omega_continuous * sampling_period

    # To create a clean, readable solution, simplify the fraction (multiplier / Fs)
    omega_fraction = Fraction(omega_multiplier, sampling_freq_hz).limit_denominator()
    num = omega_fraction.numerator
    den = omega_fraction.denominator

    # Format the omega string for the solution
    if den == 1:
        omega_discrete_str = f"{num}*pi" if num != 1 else "pi"
    elif num == 1:
        omega_discrete_str = f"pi/{den}"
    else:
        omega_discrete_str = f"({num}*pi)/{den}"

    # 3. Generate the question and solution strings
    question = (
        f"A continuous-time signal is given by the expression:\n"
        f"x_c(t) = {signal_expression}\n\n"
        f"This signal is sampled at a frequency of {sampling_freq_hz} Hz. "
        f"Determine the resulting discrete-time signal, x[n], and identify "
        f"its discrete-time angular frequency, omega."
    )

    solution = (
        f"**Given:**\n"
        f"Continuous-time signal: x_c(t) = {signal_expression}\n"
        f"Sampling frequency: Fs = {sampling_freq_hz} Hz\n\n"

        f"**Step 1:** State the Sampling Rule\n"
        f"A discrete-time signal x[n] is obtained by replacing the continuous time variable 't' "
        f"with its discrete counterpart 'nT', where T is the sampling period. The rule is:\n"
        f"x[n] = x_c(nT)\n\n"

        f"**Step 2:** Calculate Sampling Period (T)\n"
        f"The sampling period is the reciprocal of the sampling frequency.\n"
        f"T = 1 / Fs = 1 / {sampling_freq_hz} seconds.\n\n"

        f"**Step 3:** Substitute t = nT into the Equation\n"
        f"We replace every 't' in the original expression with 'nT'.\n"
        f"x[n] = {amplitude} * cos({omega_continuous_str}*(nT) + {phi_str})\n\n"

        f"**Step 4:** Identify Discrete Frequency (omega)\n"
        f"Rearrange the terms to match the standard form x[n] = A * cos(omega*n + phi). "
        f"The discrete-time angular frequency, omega, is the coefficient of 'n'.\n"
        f"omega = Omega * T\n"
        f"omega = ({omega_continuous_str}) * (1 / {sampling_freq_hz})\n\n"

        f"**Step 5:** Calculate and Conclude\n"
        f"Now, we calculate the numerical value for omega and write the final expression.\n"
        f"omega = ({omega_multiplier} / {sampling_freq_hz}) * pi = {omega_discrete_str} rad/sample\n"
        f"The final discrete-time signal is:\n"
        f"x[n] = {amplitude} * cos({omega_discrete_str}*n + {phi_str})\n\n"

        f"**Answer:**\n"
        f"The resulting discrete-time signal is x[n] = {amplitude} * cos({omega_discrete_str}*n + {phi_str}), "
        f"and its discrete-time angular frequency is omega = {omega_discrete_str} rad/sample."
    )

    return question, solution


# Template 3 (Intermediate)
def template_aliased_frequency_identification():
    """
    Identifying Aliased Frequencies

    Scenario:
        This template assesses the ability to predict the outcome of undersampling (aliasing).
        When a signal is sampled below its Nyquist rate, its frequency appears as a
        lower "folded" frequency. This function generates a problem to calculate this
        apparent frequency.

    Core Equations:
        Fa = abs(F0 - k * Fs), where k = round(F0 / Fs)

    Returns:
        tuple: A tuple containing:
            - str: A question asking for the apparent (aliased) frequency.
            - str: A step-by-step solution showing the calculation.
    """
    # 1. Parameterize the inputs with random values

    # Generate a random original frequency for the continuous-time signal.
    f0 = random.randint(*F0_RANGE_HZ)

    # Calculate the corresponding Nyquist rate for this signal.
    f_nyquist = 2 * f0

    # To guarantee aliasing, generate a sampling frequency (Fs) that is strictly
    # less than the Nyquist rate. We set the range to be between 50% and 95% of
    # the Nyquist rate to create a non-trivial undersampling problem.
    fs = random.randint(int(0.5 * f_nyquist), int(0.95 * f_nyquist))

    # 2. Perform the core calculation for the solution

    # Find the integer multiple 'k' which represents the nearest multiple of the
    # sampling frequency to the original frequency.
    k = round(f0 / fs)

    # Calculate the aliased frequency (Fa) using the folding formula.
    fa = abs(f0 - k * fs)

    # Standardize the precision for all final outputs for consistency.
    precision = 2

    # 3. Generate the question and solution strings
    question = (
        f"A continuous-time signal x_c(t) = cos(2*pi*{f0}*t) is sampled at a rate "
        f"of {fs} samples per second.\n\n"
        f"Since this rate is below the Nyquist rate, aliasing occurs. What is the "
        f"apparent frequency, Fa, of the resulting discrete-time signal? "
        f"(The apparent frequency must be in the range 0 <= Fa <= {fs/2} Hz)."
    )

    solution = (
        f"**Given:**\n"
        f"Original Frequency (F0): {f0} Hz\n"
        f"Sampling Frequency (Fs): {fs} Hz\n\n"

        f"**Step 1:** State the Condition\n"
        f"The Nyquist rate required to sample this signal without aliasing is "
        f"2 * F0 = 2 * {f0} = {f_nyquist} Hz. "
        f"Since the sampling frequency Fs = {fs} Hz is less than {f_nyquist} Hz, "
        f"the original frequency will be aliased.\n\n"

        f"**Step 2:** Explain Frequency Folding\n"
        f"When aliasing occurs, the perceived frequency (Fa) is the absolute "
        f"difference between the original frequency and the nearest integer multiple "
        f"of the sampling frequency. The formula is:\n"
        f"Fa = abs(F0 - k * Fs), where 'k' is an integer.\n\n"

        f"**Step 3:** Find the Correct Multiple (k)\n"
        f"We find the integer 'k' that 'folds' F0 into the principal frequency "
        f"range [0, Fs/2]. This is found by rounding the ratio of F0 to Fs.\n"
        f"k = round(F0 / Fs) = round({f0} / {fs}) = round({round(f0/fs, 4)}) = {k}.\n\n"

        f"**Step 4:** Calculate Aliased Frequency (Fa)\n"
        f"Now, we compute Fa using the value of k found in the previous step.\n"
        f"Fa = abs({f0} - {k} * {fs})\n"
        f"Fa = abs({f0} - {k * fs})\n"
        f"Fa = {round(fa, precision)} Hz.\n\n"

        f"**Answer:**\n"
        f"The apparent (aliased) frequency of the resulting signal is {round(fa, precision)} Hz."
    )

    return question, solution


# Template 4 (Intermediate)
def template_cd_dc_system_analysis():
    """
    C/D and D/C System Analysis

    Scenario:
        This template models a full signal processing chain: a continuous-time signal is
        sampled (C/D), processed by a simple discrete-time system (a gain and delay),
        and perfectly reconstructed back to continuous-time (D/C). It tests the
        understanding of how each stage affects the final signal.

    Core Equations:
        x[n] = x_c(nT)
        y_c(t) is the ideal reconstruction of y[n].

    Returns:
        tuple: A tuple containing:
            - str: A question describing the system and asking for the final output.
            - str: A step-by-step solution.
    """
    # 1. Parameterize the inputs with random values

    # --- Signal Parameters ---
    amplitude = round(random.uniform(*AMPLITUDE_RANGE), 2)
    omega_multiplier = random.randint(*OMEGA_MULTIPLIER_RANGE)
    omega_continuous = omega_multiplier * math.pi
    omega_continuous_str = f"{omega_multiplier}*pi"

    # --- Sampling Parameters ---
    # Ensure the Nyquist criterion is satisfied by a wide margin.
    # F0 = Omega0 / (2*pi) = (multiplier*pi) / (2*pi) = multiplier / 2.
    # F_nyquist = 2 * F0 = multiplier.
    f_nyquist = omega_multiplier
    # Sample at 3 to 8 times the Nyquist rate to make it a clear "ideal" case.
    sampling_factor = random.randint(3, 8)
    sampling_freq_hz = f_nyquist * sampling_factor
    sampling_period = 1 / sampling_freq_hz

    # --- System Parameters ---
    gain_k = round(random.uniform(*GAIN_K_RANGE), 2)
    delay_n0 = random.randint(*DELAY_N0_RANGE)

    # 2. Perform the core calculations for the solution

    # The final output amplitude is the original amplitude scaled by the system gain.
    output_amplitude = round(amplitude * gain_k, 2)

    # The total time delay is the discrete sample delay multiplied by the sampling period.
    time_delay_sec = delay_n0 * sampling_period

    # 3. Generate the question and solution strings

    # Define the strings for the mathematical expressions
    x_c_t_expr = f"{amplitude} * cos({omega_continuous_str}*t)"
    system_expr = f"y[n] = {gain_k} * x[n - {delay_n0}]"
    
    question = (
        f"A continuous-time signal is given by:\n"
        f"x_c(t) = {x_c_t_expr}\n\n"
        f"The signal is sampled with a period T = {sampling_period:.2e} s to produce the "
        f"discrete-time sequence x[n]. This sequence is then processed by a system "
        f"defined by the equation:\n"
        f"{system_expr}\n\n"
        f"If the resulting output sequence, y[n], is converted back to a continuous-time "
        f"signal, y_c(t), using an ideal D/C converter, what is the expression for y_c(t)?"
    )

    solution = (
        f"**Given:**\n"
        f"Continuous-time signal: x_c(t) = {x_c_t_expr}\n"
        f"Sampling period: T = {sampling_period:.2e} s\n"
        f"Discrete-time system: {system_expr}\n\n"

        f"**Step 1:** Continuous-to-Discrete (C/D) Conversion\n"
        f"We find the discrete-time signal x[n] by substituting t = nT into x_c(t).\n"
        f"x[n] = x_c(nT) = {amplitude} * cos({omega_continuous_str} * nT)\n\n"

        f"**Step 2:** Discrete-Time Processing\n"
        f"Next, we apply the system's equation to x[n]. This involves replacing "
        f"x[n] with the expression from Step 1 and shifting it by n0 = {delay_n0} samples.\n"
        f"y[n] = {gain_k} * x[n - {delay_n0}]\n"
        f"y[n] = {gain_k} * ({amplitude} * cos({omega_continuous_str} * (n - {delay_n0})T))\n"
        f"y[n] = {output_amplitude} * cos({omega_continuous_str}*T*(n - {delay_n0}))\n\n"

        f"**Step 3: ** Discrete-to-Continuous (D/C) Conversion*\n"
        f"An ideal D/C converter maps the discrete-time signal back to a continuous-time "
        f"sinusoid with the original frequency Omega_0 = {omega_continuous_str}. The discrete "
        f"time 'n' is mapped back to continuous time 't' via the relation t = nT.\n"
        f"y_c(t) = {output_amplitude} * cos({omega_continuous_str} * (t - {delay_n0}T))\n"
        f"The term {delay_n0}*T represents a time delay. Let's calculate its value:\n"
        f"Time Delay = {delay_n0} * T = {delay_n0} * {sampling_period:.2e} = {time_delay_sec:.2e} s\n"
        f"Substituting this back, we get the final expression:\n"
        f"y_c(t) = {output_amplitude} * cos({omega_continuous_str}*(t - {time_delay_sec:.2e}))\n\n"

        f"**Step 4:** Interpret the Result\n"
        f"The output signal y_c(t) is a modified version of the input signal x_c(t). "
        f"Its amplitude has been scaled by a factor of K = {gain_k}, and it has been "
        f"time-delayed by {time_delay_sec:.2e} seconds.\n\n"

        f"**Answer:**\n"
        f"The final continuous-time output signal is y_c(t) = {output_amplitude} * cos({omega_continuous_str}*(t - {time_delay_sec:.2e}))."
    )

    return question, solution


# Template 5 (Advanced)
def template_decimation_aliasing_analysis():
    """
    Decimation (Downsampling) with Aliasing Analysis

    Scenario:
        This template tests the understanding of decimation and its potential to cause
        aliasing if the signal is not properly band-limited. The generated problem may
        or may not result in aliasing, requiring a full analysis.

    Core Equations:
        y[n] = x[Mn]
        Aliasing is avoided if omega_0 < pi / M.
        omega_aliased = abs(omega_prime - 2*pi*k), where omega_prime = omega_0 * M.

    Returns:
        tuple: A tuple containing:
            - str: A multi-part question about the downsampled signal.
            - str: A step-by-step solution with conditional aliasing analysis.
    """
    # 1. Parameterize the inputs with random values
    M = random.randint(*DECIMATION_FACTOR_M_RANGE)
    
    # Randomly decide whether to create a problem that results in aliasing or not.
    will_alias = random.choice([True, False])

    # Generate an initial frequency omega_0 = (num/den)*pi based on the aliasing choice.
    # The condition for NO aliasing is omega_0 < pi/M, which means num/den < 1/M.
    denominator = random.randint(*OMEGA_DENOMINATOR_RANGE)
    
    if not will_alias:
        # We need num/den < 1/M => num < den/M.
        max_numerator = math.floor(denominator / M) - 1
        if max_numerator < 1: max_numerator = 1 # Ensure at least one valid choice
        numerator = random.randint(1, max_numerator)
    else:
        # We need num/den > 1/M => num > den/M.
        min_numerator = math.ceil(denominator / M)
        if min_numerator >= denominator: min_numerator = denominator - 1 # Ensure valid range
        numerator = random.randint(min_numerator, denominator - 1)

    # The original discrete-time frequency
    omega_0 = (numerator / denominator) * math.pi

    # --- In-line formatting for omega_0_str ---
    frac_0 = Fraction(numerator, denominator)
    n0, d0 = frac_0.numerator, frac_0.denominator
    if d0 == 1:
        omega_0_str = f"{n0}*pi" if n0 != 1 else "pi"
    elif n0 == 1:
        omega_0_str = f"pi/{d0}"
    else:
        omega_0_str = f"({n0}*pi)/{d0}"

    # 2. Perform the core calculations for the solution
    aliasing_boundary = math.pi / M
    
    # --- In-line formatting for aliasing_boundary_str ---
    aliasing_boundary_str = f"pi/{M}"

    # The new frequency after substitution before considering the principal range.
    omega_prime = omega_0 * M

    # --- In-line formatting for omega_prime_str ---
    frac_prime = Fraction(numerator * M, denominator)
    np, dp = frac_prime.numerator, frac_prime.denominator
    if dp == 1:
        omega_prime_str = f"{np}*pi" if np != 1 else "pi"
    elif np == 1:
        omega_prime_str = f"pi/{dp}"
    else:
        omega_prime_str = f"({np}*pi)/{dp}"
    
    # Calculate the final perceived frequency.
    if not will_alias:
        final_omega_str = omega_prime_str
        k_val = 0 # No folding needed
    else:
        # The frequency must be folded back into the principal range [-pi, pi].
        k_val = round(omega_prime / (2 * math.pi))
        final_omega = abs(omega_prime - 2 * math.pi * k_val)
        
        # --- In-line formatting for final_omega_str ---
        frac_final = Fraction(final_omega / math.pi).limit_denominator(100)
        nf, df = frac_final.numerator, frac_final.denominator
        if df == 1:
            final_omega_str = f"{nf}*pi" if nf != 1 else "pi"
        elif nf == 0:
            final_omega_str = "0"
        elif nf == 1:
            final_omega_str = f"pi/{df}"
        else:
            final_omega_str = f"({nf}*pi)/{df}"

    # 3. Generate the question and solution strings
    question = (
        f"A discrete-time signal is given by x[n] = cos({omega_0_str}*n).\n"
        f"This signal is downsampled by a factor of M = {M} to produce the signal y[n] = x[Mn].\n\n"
        f"a) Find an expression for the output signal y[n] in the form cos(omega'*n).\n"
        f"b) The condition to avoid aliasing during downsampling is omega_0 < {aliasing_boundary_str}. Based on the given values, did aliasing occur?\n"
        f"c) What is the final perceived angular frequency, omega_a, in the signal y[n]? (The frequency must be in the range 0 <= omega_a <= pi)."
    )
    
    # Build the solution string piece by piece
    solution = (
        f"**Given:**\n"
        f"Signal: x[n] = cos({omega_0_str}*n)\n"
        f"Downsampling Factor: M = {M}\n\n"
        
        f"**Step 1:** Find the Output Signal Expression\n"
        f"We substitute 'n' with 'Mn' in the original expression for x[n]:\n"
        f"  y[n] = x[Mn] = cos({omega_0_str} * (Mn)) = cos(({omega_0_str} * {M}) * n)\n"
        f"The new angular frequency, omega', is:\n"
        f"  omega' = {omega_0_str} * {M} = {omega_prime_str}\n"
        f"So, the expression is y[n] = cos({omega_prime_str}*n).\n\n"

        f"**Step 2:** Check for Aliasing\n"
        f"Aliasing is avoided if the original frequency omega_0 is less than pi/M.\n"
        f"Aliasing Boundary = pi / M = {aliasing_boundary_str} approx {round(aliasing_boundary, 3)}\n"
        f"Original Frequency = {omega_0_str} approx {round(omega_0, 3)}\n"
        f"Comparing the two, we find that {omega_0_str} is {'NOT less' if will_alias else 'less'} than {aliasing_boundary_str}. "
        f"Therefore, aliasing **{'occurred' if will_alias else 'did not occur'}**.\n\n"

        f"**Step 3:** Find the Final Perceived Frequency (omega_a)\n"
    )
    
    if not will_alias:
        solution += (
            f"Since no aliasing occurred, the frequency omega' = {omega_prime_str} is already in the principal range [0, pi].\n"
            f"Therefore, the final perceived frequency is the same.\n"
            f"omega_a = {final_omega_str}"
        )
    else:
        # For clarity, re-calculate final_omega inside the solution string
        final_omega_val = abs(omega_prime - 2 * math.pi * k_val)
        solution += (
            f"Because aliasing occurred, the frequency omega' = {omega_prime_str} is outside the principal range [0, pi] and must be 'folded' back.\n"
            f"We find an integer 'k' such that omega_a = abs(omega' - 2*pi*k) is in the range [0, pi].\n"
            f"k = round(omega' / (2*pi)) = round({omega_prime/math.pi:.2f}*pi / (2*pi)) = round({omega_prime/(2*math.pi):.2f}) = {k_val}\n"
            f"Now we calculate omega_a:\n"
            f"omega_a = abs({omega_prime_str} - 2*pi*{k_val})\n"
            f"omega_a = abs({round(omega_prime, 3)} - {round(2*math.pi*k_val, 3)}) = {round(final_omega_val, 3)}\n"
            f"In terms of pi, this is omega_a = {final_omega_str}.\n\n"
        )
        
    solution += (
        f"\n\n**Answer:**\n"
        f"a) The output signal is y[n] = cos({omega_prime_str}*n).\n"
        f"b) Aliasing **{'did' if will_alias else 'did not'}** occur.\n"
        f"c) The final perceived frequency is omega_a = {final_omega_str}."
    )

    return question, solution


def main():
    """
    Generate numerous instances of each continuous time signals template 
    with different random seeds and write the results to a JSONL file.
    """
    import json
    import os

    # Define the output path (Modify this path according to where you are running the code from)
    output_file = "testset/electrical_engineering/signals_and_systems/continuous_time_signals.jsonl"

    # Create the directory if it doesn't exist
    os.makedirs(os.path.dirname(output_file), exist_ok=True)

    # List of template functions with their ID and level
    templates = [
        (template_nyquist_rate_determination, "nyquist_rate_determination", "Easy"),
        (template_continuous_to_discrete_conversion, "continuous_to_discrete_conversion", "Easy"),
        (template_aliased_frequency_identification, "aliased_frequency_identification", "Intermediate"),
        (template_cd_dc_system_analysis, "cd_dc_system_analysis", "Intermediate"),
        (template_decimation_aliasing_analysis, "decimation_aliasing_analysis", "Advanced"),
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
