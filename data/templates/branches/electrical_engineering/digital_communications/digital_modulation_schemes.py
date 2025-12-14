import random
import math


# Template 1 (Easy)
def template_bpsk_energy_basis():
    """
    BPSK/BASK Signal Energy and Basis Function

    Scenario:
        This template tests the understanding of how to analyze a basic modulated
        signal, which is a fundamental first step in signal space analysis. It
        involves calculating the energy of the signal (Eb) over one bit period
        and deriving its corresponding orthonormal basis function (psi(t)). These
        two components are essential for representing signals as vectors.

    Core Equations:
        Signal: s(t) = A * cos(2*pi*fc*t)
        Energy Integral: Eb = integral from 0 to Tb of [s(t)]^2 dt
        Simplified Energy: Eb = (A^2 * Tb) / 2
        Basis Function: psi_1(t) = s(t) / sqrt(Eb) = sqrt(2 / Tb) * cos(2*pi*fc*t)

    Returns:
        tuple: A tuple containing:
            - str: A question about BPSK signal properties.
            - str: A step-by-step solution.
    """
    # 1. Parameterize the inputs with random values
    
    # Set a standard precision for rounding in all calculations and outputs
    precision = 2

    amplitude = round(random.uniform(1.0, 10.0), precision)
    
    # Generate bit duration in a range from 1 microsecond to 10 milliseconds
    bit_duration_s = random.uniform(1e-6, 1e-2)
    
    # Choose a carrier frequency that is an integer multiple of the bit rate (1/Tb)
    k = random.randint(2, 20)
    carrier_freq_hz = k / bit_duration_s

    # Inline formatting for bit_duration_s 
    prefixes = {6: 'M', 3: 'k', 0: '', -3: 'm', -6: 'u', -9: 'n'}
    exponent_td = int(math.floor(math.log10(abs(bit_duration_s)) / 3.0) * 3)
    prefix_td = prefixes.get(exponent_td, '')
    scaled_td = bit_duration_s / (10**exponent_td)
    bit_duration_str = f"{round(scaled_td, precision)} {prefix_td}s"

    # Inline formatting for carrier_freq_hz 
    exponent_fc = int(math.floor(math.log10(abs(carrier_freq_hz)) / 3.0) * 3)
    prefix_fc = prefixes.get(exponent_fc, '')
    scaled_fc = carrier_freq_hz / (10**exponent_fc)
    carrier_freq_str = f"{round(scaled_fc, precision)} {prefix_fc}Hz"

    # 2. Perform the core calculation
    
    # Calculate energy per bit: Eb = (A^2 * Tb) / 2
    energy_joules = (amplitude**2 * bit_duration_s) / 2
    
    # Calculate the amplitude of the basis function: sqrt(2 / Tb)
    basis_amplitude = math.sqrt(2 / bit_duration_s)

    # Inline formatting for energy_joules 
    exponent_e = int(math.floor(math.log10(abs(energy_joules)) / 3.0) * 3)
    prefix_e = prefixes.get(exponent_e, '')
    scaled_e = energy_joules / (10**exponent_e)
    energy_str = f"{round(scaled_e, precision)} {prefix_e}J"
    
    # 3. Generate the question and solution strings

    question = (
        f"A binary ASK (BASK) signal is defined by the following equation for one bit interval:\n"
        f"s(t) = {amplitude} * cos(2*pi*{carrier_freq_str}*t)\n"
        f"for the time interval 0 <= t <= {bit_duration_str}.\n\n"
        f"Calculate the following:\n"
        f"a) The energy per bit (Eb).\n"
        f"b) The orthonormal basis function psi_1(t)."
    )

    solution = (
        f"**Given Signal:**\n"
        f"s(t) = {amplitude} * cos(2*pi*{carrier_freq_str}*t) for 0 <= t <= {bit_duration_str}\n\n"
        
        f"**Step 1:** Identify Parameters\n"
        f"From the signal definition, we can extract the following parameters:\n"
        f"Amplitude (A): {amplitude} V\n"
        f"Bit Duration (Tb): {bit_duration_str} ({bit_duration_s:.4e} s)\n"
        f"Carrier Frequency (fc): {carrier_freq_str} ({carrier_freq_hz:.4e} Hz)\n\n"
        
        f"**Step 2:** Calculate Energy per Bit (Eb)\n"
        f"The energy of a signal is the integral of its squared magnitude over its duration.\n"
        f"Eb = integral from 0 to Tb of [s(t)]^2 dt\n"
        f"For a sinusoidal signal of the form A*cos(...), the energy over an interval Tb simplifies to:\n"
        f"Eb = (A^2 * Tb) / 2\n"
        f"Substituting the values:\n"
        f"Eb = ({amplitude}^2 * {bit_duration_s:.4e}) / 2\n"
        f"Eb = {energy_joules:.4e} J\n"
        f"Eb = {energy_str}\n\n"

        f"**Step 3:** Determine the Basis Function (psi_1(t))\n"
        f"The orthonormal basis function is found by normalizing the signal by the square root of its energy.\n"
        f"psi_1(t) = s(t) / sqrt(Eb)\n"
        f"We can use the general formula derived from this relationship:\n"
        f"psi_1(t) = sqrt(2 / Tb) * cos(2*pi*fc*t)\n"
        f"First, calculate the amplitude of the basis function:\n"
        f"Amplitude_psi = sqrt(2 / {bit_duration_s:.4e} s) = {round(basis_amplitude, precision)}\n"
        f"Now, construct the full basis function:\n"
        f"psi_1(t) = {round(basis_amplitude, precision)} * cos(2*pi*{carrier_freq_str}*t)\n\n"

        f"**Answer:**\n"
        f"a) The energy per bit (Eb) is {energy_str}.\n"
        f"b) The basis function (psi_1(t)) is {round(basis_amplitude, precision)} * cos(2*pi*{carrier_freq_str}*t)."
    )
    
    return question, solution


# Template 2 (Intermediate)
def template_euclidean_distance_binary():
    """
    Euclidean Distance for Binary Modulation

    Scenario:
        This template tests the ability to represent signals as vectors in a signal
        space and calculate the Euclidean distance between them. This distance is a
        key factor in determining the noise immunity of a modulation scheme.

    Core Equations:
        For BPSK (antipodal): d = 2 * sqrt(Eb)
        For BFSK (orthogonal): d = sqrt(2 * Eb)

    Returns:
        tuple: A tuple containing:
            - str: A question about Euclidean distance.
            - str: A step-by-step solution.
    """
    # 1. Parameterize the inputs with random values
    
    precision = 3
    modulation_type = random.choice(['BPSK', 'Orthogonal BFSK'])
    
    # Generate energy per bit, Eb, in a range from picojoules to nanojoules
    energy_joules = random.uniform(1e-12, 1e-9)
    
    # Inline formatting for energy_joules 
    prefixes = {0: '', -3: 'm', -6: 'u', -9: 'n', -12: 'p'}
    if energy_joules > 0:
        exponent_e = int(math.floor(math.log10(abs(energy_joules)) / 3.0) * 3)
        prefix_e = prefixes.get(exponent_e, '')
        scaled_e = energy_joules / (10**exponent_e)
        energy_str = f"{scaled_e:.{precision}f} {prefix_e}J"
    else:
        energy_str = "0 J"
    
    # 2. Perform the core calculation based on modulation type
    
    sqrt_eb = math.sqrt(energy_joules)

    # Use scientific notation formatting for small numbers instead of round()
    sqrt_eb_fmt = f"{sqrt_eb:.{precision}e}"

    if modulation_type == 'BPSK':
        # BPSK signals are antipodal (180 degrees apart)
        constellation_type = "antipodal"
        dimension = "one-dimensional"
        basis_count = "one basis function, psi_1(t)"
        
        s1_str = f"({sqrt_eb_fmt})"
        s2_str = f"(-{sqrt_eb_fmt})"
        
        distance = 2 * sqrt_eb
        dist_fmt = f"{distance:.{precision}e}"
        
        derivation_steps = (
            f"The signal vectors are s1 = (sqrt(Eb)) and s2 = (-sqrt(Eb)).\n"
            f"s1 = ({sqrt_eb_fmt})\n"
            f"s2 = (-{sqrt_eb_fmt})\n\n"
            

            f"**Step 2:** Calculate Euclidean Distance (d)\n"
            f"The distance d is the magnitude of the difference vector (s1 - s2).\n"
            f"s1 - s2 = (sqrt(Eb) - (-sqrt(Eb))) = (2*sqrt(Eb))\n"
            f"d = ||s1 - s2|| = 2*sqrt(Eb)\n"
            f"d = 2 * {sqrt_eb_fmt}\n"
            f"d = {dist_fmt}"
        )

    else: # Orthogonal BFSK
        # BFSK signals are orthogonal (90 degrees apart)
        constellation_type = "orthogonal"
        dimension = "two-dimensional"
        basis_count = "two basis functions, psi_1(t) and psi_2(t)"

        s1_str = f"({sqrt_eb_fmt}, 0)"
        s2_str = f"(0, {sqrt_eb_fmt})"

        distance = math.sqrt(2 * energy_joules)
        dist_fmt = f"{distance:.{precision}e}"

        derivation_steps = (
            f"The signal vectors are s1 = (sqrt(Eb), 0) and s2 = (0, sqrt(Eb)).\n"
            f"s1 = ({sqrt_eb_fmt}, 0)\n"
            f"s2 = (0, {sqrt_eb_fmt})\n\n"
            

            f"**Step 2: ** Calculate Euclidean Distance (d)\n"
            f"The distance d is the magnitude of the difference vector (s1 - s2).\n"
            f"s1 - s2 = (sqrt(Eb) - 0, 0 - sqrt(Eb)) = (sqrt(Eb), -sqrt(Eb))\n"
            f"d = ||s1 - s2|| = sqrt( (sqrt(Eb))^2 + (-sqrt(Eb))^2 ) = sqrt(2*Eb)\n"
            f"d = sqrt(2 * {energy_joules:.{precision}e})\n"
            f"d = {dist_fmt}"
        )

    # 3. Generate the question and solution strings
    
    question = (
        f"A digital communication system uses {modulation_type} modulation.\n"
        f"The energy per bit is Eb = {energy_str}.\n\n"
        f"Determine the following:\n"
        f"a) The signal constellation points (s1 and s2) in vector form.\n"
        f"b) The Euclidean distance (d) between these two points."
    )
    
    solution = (
        f"**Given Information:**\n"
        f"Modulation Scheme: {modulation_type}\n"
        f"Energy per Bit (Eb): {energy_str} ({energy_joules:.{precision}e} J)\n\n"
        
        f"**Step 1:** Define Signal Vectors\n"
        f"For {modulation_type} modulation, the signals are {constellation_type}. This means we can represent them in a {dimension} signal space using {basis_count}.\n"
        f"First, we calculate sqrt(Eb) = sqrt({energy_joules:.{precision}e}) = {sqrt_eb_fmt}.\n"
        f"{derivation_steps}\n\n"
        
        f"**Answer:**\n"
        f"a) The signal constellation points are:\n"
        f"s1 = {s1_str}\n"
        f"s2 = {s2_str}\n"
        f"b) The Euclidean distance between the points is {dist_fmt}."
    )

    return question, solution


# Template 3 (Intermediate)
def template_average_energy_mqam():
    """
    Average Energy of an M-QAM Constellation

    Scenario:
        This template tests the ability to calculate the average symbol energy by
        analyzing the geometry of the constellation.

    Core Equations:
        General Formula: E_avg = (2/3) * (M - 1) * A^2

    Returns:
        tuple: A tuple containing:
            - str: A question asking for the average energy.
            - str: A step-by-step solution.
    """
    # 1. Parameterize the inputs with random values
    precision = 2
    # Add 256-QAM to the list of possible modulation orders
    M = random.choice([4, 16, 64, 256])
    
    # Randomly decide whether A is an integer or a float
    if random.choice([True, False]):
        A = random.randint(1, 10)
    else:
        A = round(random.uniform(0.2, 10.0), precision)
    
    # 2. Perform the core calculation based on M
    
    # The general formula for average energy in a square M-QAM is (2/3)(M-1)A^2
    avg_energy_coeff = (2/3) * (M - 1)
    avg_energy = avg_energy_coeff * (A**2)

    # Build the solution steps string based on the value of M
    if M == 4:
        coordinate_set = "{+/-A}"
        solution_steps = (
            f"**Step 2:** List Signal Point Energies\n"
            f"For 4-QAM, there is only one type of signal point, located at coordinates (+/-A, +/-A).\n"
            f"There are 4 points, all with energy E = A^2 + A^2 = 2A^2.\n\n"
            
            f"**Step 3:** Calculate Average Energy (E_avg)\n"
            f"We sum the energies of all points and divide by M.\n"
            f"Total Energy = 4 * (2 * A^2) = 8 * A^2\n"
            f"E_avg = (Total Energy) / M = (8 * A^2) / 4 = 2 * A^2\n"
        )
    
    elif M == 16:
        coordinate_set = "{+/-A, +/-3A}"
        solution_steps = (
            f"**Step 2:** List Signal Point Energies\n"
            f"For 16-QAM, the points can be grouped by their distance from the origin (their energy).\n"
            f"4 points at (+/-A, +/-A) have energy E1 = A^2 + A^2 = 2A^2.\n"
            f"8 points at (+/-A, +/-3A) or (+/-3A, +/-A) have energy E2 = A^2 + (3A)^2 = 10A^2.\n"
            f"4 points at (+/-3A, +/-3A) have energy E3 = (3A)^2 + (3A)^2 = 18A^2.\n\n"

            f"**Step 3:** Calculate Average Energy (E_avg)\n"
            f"We sum the energies of all points and divide by M.\n"
            f"Total Energy = 4*(2A^2) + 8*(10A^2) + 4*(18A^2) = 160A^2\n"
            f"E_avg = (Total Energy) / M = (160 * A^2) / 16 = 10 * A^2\n"
        )
        
    elif M == 64:
        coordinate_set = "{+/-A, +/-3A, +/-5A, +/-7A}"
        solution_steps = (
            "**Step 2:** List Signal Point Energies\n"
            "For 64-QAM, there are many groups of points with the same energy. We list a few examples:\n"
            "- The 4 innermost points at (+/-A, +/-A) have energy E = 2A^2.\n"
            "- The 4 outermost points at (+/-7A, +/-7A) have energy E = (7A)^2 + (7A)^2 = 98A^2.\n"
            "This process is continued for all 64 points.\n\n"

            "**Step 3:** Calculate Average Energy (E_avg)\n"
            "Summing the energies for all 64 points and dividing by M gives the average. We use the general formula for square M-QAM: E_avg = (2/3)*(M-1)*A^2.\n"
            "E_avg = (2/3) * (64 - 1) * A^2\n"
            "E_avg = (2/3) * 63 * A^2 = 42 * A^2\n"
        )

    else: # M == 256 (Explicit handling for 256-QAM)
        coordinate_set = "{+/-A, +/-3A, ..., +/-15A}"
        solution_steps = (
            "**Step 2:** List Signal Point Energies\n"
            "For 256-QAM, the grid extends from -15A to +15A. Calculating individual point energies is tedious, so we rely on the general formula derived from the sum of squares.\n"
            "- Innermost points: (+/-A, +/-A)\n"
            "- Outermost points: (+/-15A, +/-15A)\n\n"

            "**Step 3:** Calculate Average Energy (E_avg)\n"
            "Using the general formula for square M-QAM: E_avg = (2/3)*(M-1)*A^2.\n"
            "E_avg = (2/3) * (256 - 1) * A^2\n"
            "E_avg = (2/3) * 255 * A^2 = 170 * A^2\n"
        )

    # 3. Generate the question and solution strings
    
    question = (
        f"For a square {M}-QAM constellation, the signal points are located at (xi, yj) where the coordinates xi and yj are chosen from the set {coordinate_set}.\n\n"
        f"The fundamental distance parameter is A = {A}.\n\n"
        f"Calculate the average energy per symbol (E_avg) for this constellation."
    )
    
    solution = (
        f"**Given Information:**\n"
        f"Modulation Scheme: {M}-QAM\n"
        f"Distance Parameter (A): {A}\n"
        f"\n\n"
        
        f"**Step 1:** Understand the Constellation Structure\n"
        f"The constellation is a square grid where coordinates are odd multiples of A. The energy of any point (x, y) is simply x^2 + y^2.\n\n"
        
        f"{solution_steps}"
        
        f"**Step 4:** Final Calculation\n"
        f"Now, we substitute the value of A = {A}.\n"
        f"E_avg = {round(avg_energy_coeff, precision)} * ({A})^2\n"
        f"E_avg = {round(avg_energy_coeff, precision)} * {round(A**2, precision)}\n"
        f"E_avg = {round(avg_energy, precision)}\n\n"

        f"**Answer:**\n"
        f"The average energy per symbol is {round(avg_energy, precision)}."
    )
    
    return question, solution


# Template 4 (Advanced)
def template_ber_estimation_mary():
    """
    Bit Error Rate (BER) Estimation for M-ary Schemes

    Scenario:
        Predicting the performance of a digital communication system in the presence of
        noise is a critical engineering task. This template tests the ability to
        estimate the bit error rate (BER) from a given signal-to-noise ratio per
        bit (Eb/N0) for common M-ary modulation schemes. It requires using standard
        approximation formulas that are expressed in terms of the complementary
        Gaussian distribution function, Q(x).

    Core Equations:
        Linear SNR: (Eb/N0)_lin = 10^((Eb/N0)_dB / 10)
        Bits per symbol: k = log2(M)
        M-PSK SER: Ps approx 2 * Q(sqrt(2 * Es/N0) * sin(pi/M))
        M-QAM SER: Ps approx 4 * (1 - 1/sqrt(M)) * Q(sqrt(3*k/(M-1) * Eb/N0))
        Gray Code BER: Pb approx Ps / k

    Returns:
        tuple: A tuple containing:
            - str: A question asking for the BER estimation.
            - str: A step-by-step solution showing the calculation.
    """
    # 1. Enhanced Parameter Randomization
    precision = random.randint(1, 2)
    modulation_scheme = random.choice(['M-PSK', 'M-QAM'])

    if modulation_scheme == 'M-PSK':
        M = random.choice([4, 8, 16, 32])
    else:  # M-QAM
        # 4-QAM is identical to 4-PSK (QPSK), so we start at 16 for variety.
        M = random.choice([16, 64, 256])

    eb_n0_db = round(random.uniform(5.0, 25.0), precision)

    # 2. Perform the core calculation
    eb_n0_lin = 10**(eb_n0_db / 10)
    k = math.log2(M)

    # Logic varies based on modulation scheme
    if modulation_scheme == 'M-PSK':
        # For M-PSK, the formula uses Es/N0
        es_n0_lin = k * eb_n0_lin
        q_arg = math.sqrt(2 * es_n0_lin) * math.sin(math.pi / M)
        ser_coeff = 2
        ber_coeff = ser_coeff / k

        ser_approx_str = f"{ser_coeff} * Q({q_arg:.{precision+1}f})"
        ber_approx_str = f"{ber_coeff:.{precision+1}f} * Q({q_arg:.{precision+1}f})"
        
        calculation_steps = (
            f"**Step 2:** Calculate Symbol SNR (Es/N0)\n"
            f"The number of bits per symbol is k = log2({M}) = {k:.0f}.\n"
            f"The energy per symbol (Es) is k times the energy per bit (Eb).\n"
            f"Es/N0 = k * (Eb/N0) = {k:.0f} * {eb_n0_lin:.{precision}f} = {es_n0_lin:.{precision}f}\n\n"
            
            f"**Step 3:** Calculate Symbol Error Rate (SER)\n"
            f"For M-PSK, the SER is approximated by: Ps approx 2*Q(sqrt(2*Es/N0)*sin(pi/M)).\n"
            f"First, calculate the argument of the Q-function:\n"
            f"arg = sqrt(2 * {es_n0_lin:.{precision}f}) * sin(pi / {M})\n"
            f"arg = {q_arg:.{precision+1}f}\n"
            f"So, the SER is: Ps approx {ser_approx_str}\n"
        )

    else:  # M-QAM
        # For square M-QAM, the formula directly uses Eb/N0
        q_arg = math.sqrt((3 * k / (M - 1)) * eb_n0_lin)
        ser_coeff = 4 * (1 - 1/math.sqrt(M))
        ber_coeff = ser_coeff / k
        
        ser_approx_str = f"{ser_coeff:.{precision+1}f} * Q({q_arg:.{precision+1}f})"
        ber_approx_str = f"{ber_coeff:.{precision+1}f} * Q({q_arg:.{precision+1}f})"

        calculation_steps = (
            f"**Step 2:** Calculate Bits per Symbol (k)\n"
            f"The number of bits per symbol is k = log2({M}) = {k:.0f}.\n\n"

            f"**Step 3:** Calculate Symbol Error Rate (SER)\n"
            f"For square M-QAM, the SER is approximated by: Ps approx 4*(1-1/sqrt(M))*Q(sqrt(3*k/(M-1) * Eb/N0)).\n"
            f"First, calculate the argument of the Q-function:\n"
            f"arg = sqrt( (3*{k:.0f} / ({M}-1)) * {eb_n0_lin:.{precision}f} )\n"
            f"arg = {q_arg:.{precision+1}f}\n"
            f"So, the SER is: Ps approx {ser_approx_str}\n"
        )
        
    # 3. Generate the question and solution strings
    question = (
        f"A digital communication system uses {M}-{modulation_scheme.split('-')[1]} modulation over an AWGN channel.\n"
        f"The system operates with a signal-to-noise ratio per bit of Eb/N0 = {eb_n0_db} dB.\n\n"
        "Assuming a Gray code mapping is used, provide a mathematical expression for the estimated Bit Error Rate (BER)."
    )

    solution = (
        f"**Given Information:**\n"
        f"Modulation: {M}-{modulation_scheme.split('-')[1]}\n"
        f"Eb/N0: {eb_n0_db} dB\n\n"
        
        f"**Step 1:** Convert SNR from dB to Linear Scale\n"
        f"The linear value of Eb/N0 is required for the error rate formulas.\n"
        f"Eb/N0 (linear) = 10^( ({eb_n0_db}) / 10 ) = {eb_n0_lin:.{precision}f}\n\n"

        f"{calculation_steps}\n"
        
        f"**Step 4:** Approximate Bit Error Rate (BER)\n"
        f"For Gray-coded constellations, the BER can be approximated from the SER by Pb approx Ps / k.\n"
        f"Pb approx ({ser_approx_str}) / {k:.0f}\n"
        f"Pb approx {ber_approx_str}\n\n"
        
        f"**Answer:**\n"
        f"The estimated Bit Error Rate (BER) is expressed as:\n"
        f"BER approx {ber_approx_str}"
    )

    return question, solution


# Template 5 (Advanced)
def template_null_to_null_bandwidth():
    """
    Null-to-Null Bandwidth Calculation

    Scenario:
        Understanding the bandwidth requirements of a signal is essential for designing
        and analyzing communication systems. The null-to-null bandwidth of the main
        spectral lobe is a fundamental, albeit theoretical, measure of the spectrum
        occupied by a signal shaped with basic rectangular pulses. This template tests
        the ability to calculate this bandwidth by relating the user data rate (Rb)
        to the transmission symbol rate (Rs) for various M-ary modulation schemes.

    Core Equations:
        Bits per Symbol: k = log2(M)
        Symbol Rate: Rs = Rb / k
        Null-to-Null RF Bandwidth: B_null = 2 * Rs

    Returns:
        tuple: A tuple containing:
            - str: A question asking for the null-to-null bandwidth.
            - str: A step-by-step solution showing the calculation.
    """
    # 1. Enhanced Parameter Randomization
    precision = 2
    modulation_scheme = random.choice(['M-PSK', 'M-QAM', 'M-PAM'])

    # Choose M from a scheme-appropriate set
    if modulation_scheme == 'M-PSK':
        M = random.choice([2, 4, 8, 16, 32, 64, 128])  
    elif modulation_scheme == 'M-QAM':
        M = random.choice([4, 16, 32, 64, 128, 256, 512, 1024])  
    else:  # M-PAM
        M = random.choice([2, 4, 8, 16, 32, 64]) 

    # Randomize bit rate value and unit for diversity
    unit_choice = random.choice(['kbps', 'Mbps', 'Gbps'])
    if unit_choice == 'kbps':
        value = random.randint(100, 999)
        multiplier = 1e3
    elif unit_choice == 'Mbps':
        value = random.choice([1, 1.5, 2, 5, 10, 25, 50, 100])
        multiplier = 1e6
    else:  # Gbps
        value = random.choice([1, 2.5, 10])
        multiplier = 1e9

    bit_rate_str = f"{value} {unit_choice}"
    bit_rate_hz = value * multiplier

    # 2. Perform the core calculation
    k = math.log2(M)
    symbol_rate_rs = bit_rate_hz / k
    bandwidth_hz = 2 * symbol_rate_rs

    #  Start: Inline formatting for bandwidth_hz 
    prefixes = {12: 'T', 9: 'G', 6: 'M', 3: 'k', 0: ''}
    if bandwidth_hz > 0:
        exponent_b = int(math.floor(math.log10(abs(bandwidth_hz)) / 3.0) * 3)
        prefix_b = prefixes.get(exponent_b, '')
        scaled_b = bandwidth_hz / (10**exponent_b)
        bandwidth_str = f"{round(scaled_b, precision)} {prefix_b}Hz"
    else:
        bandwidth_str = "0 Hz"
    #  End: Inline formatting 

    #  Start: Inline formatting for symbol_rate_rs 
    if symbol_rate_rs > 0:
        exponent_rs = int(math.floor(math.log10(abs(symbol_rate_rs)) / 3.0) * 3)
        prefix_rs = prefixes.get(exponent_rs, '')
        scaled_rs = symbol_rate_rs / (10**exponent_rs)
        symbol_rate_str = f"{round(scaled_rs, precision)} {prefix_rs}symbols/s"
    else:
        symbol_rate_str = "0 symbols/s"
    #  End: Inline formatting 


    # 3. Generate the question and solution strings
    question = (
        f"A digital communication system transmits data at a rate of {bit_rate_str} using {M}-{modulation_scheme.split('-')[1]} modulation.\n\n"
        "Assuming the signal is shaped with rectangular baseband pulses, what is the null-to-null RF bandwidth of the transmitted signal's main spectral lobe?"
    )

    solution = (
        f"**Given Information:**\n"
        f"Modulation: {M}-{modulation_scheme.split('-')[1]}\n"
        f"Bit Rate (Rb): {bit_rate_str} ({bit_rate_hz:.2e} bps)\n\n"

        f"**Step 1:** Calculate Bits per Symbol (k)\n"
        f"This value determines how many bits are grouped together to form a single symbol.\n"
        f"k = log2(M) = log2({M}) = {k:.0f} bits/symbol\n\n"

        f"**Step 2:** Calculate Symbol Rate (Rs)\n"
        f"The symbol rate (or baud rate) is the rate at which symbols are transmitted. It is found by dividing the bit rate by the number of bits per symbol.\n"
        f"Rs = Rb / k\n"
        f"Rs = {bit_rate_hz:.2e} bps / {k:.0f} bits/symbol\n"
        f"Rs = {symbol_rate_rs:.2e} symbols/s = {symbol_rate_str}\n\n"

        f"**Step 3:** Calculate Null-to-Null Bandwidth (B_null)\n"
        f"For a signal using rectangular pulses of duration Ts, the baseband spectrum is a sinc function with its first null at frequency 1/Ts. The symbol rate Rs is equal to 1/Ts.\n"
        f"For a passband RF signal, the main lobe bandwidth is twice the baseband null frequency.\n"
        f"B_null = 2 * (1/Ts) = 2 * Rs\n"
        f"B_null = 2 * {symbol_rate_rs:.2e} symbols/s\n"
        f"B_null = {bandwidth_hz:.2e} Hz = {bandwidth_str}\n\n"

        f"**Answer:**\n"
        f"The null-to-null RF bandwidth is {bandwidth_str}."
    )

    return question, solution


def main():
    """
    Generate numerous instances of each digital modulation schemes template 
    with different random seeds and write the results to a JSONL file.
    """
    import json
    import os

    # Define the output path (Modify this path according to where you are running the code from)
    output_file = "testset/electrical_engineering/digital_communications/digital_modulation_schemes.jsonl"

    # Create the directory if it doesn't exist
    os.makedirs(os.path.dirname(output_file), exist_ok=True)

    # List of template functions with their ID and level
    templates = [
        (template_bpsk_energy_basis, "bpsk_energy_basis", "Easy"),
        (template_euclidean_distance_binary, "euclidean_distance_binary", "Intermediate"),
        (template_average_energy_mqam, "average_energy_mqam", "Intermediate"),
        (template_ber_estimation_mary, "ber_estimation_mary", "Advanced"),
        (template_null_to_null_bandwidth, "null_to_null_bandwidth", "Advanced"),
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
                "area": "digital_modulation_schemes",
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
