import random
import math
import numpy as np
from data.templates.branches.electrical_engineering.constants import EPSILON_0


# Template 1 (Easy)
def template_coulombs_law():
    """
    Coulomb's Law for Two Point Charges 

    Scenario:
        This template tests the fundamental calculation of the electrostatic force
        between two point charges in a dielectric medium. It requires finding the
        separation vector between the charges and applying the vector form of
        Coulomb's Law.

    Core Equations:
        R_vec = P2 - P1
        F_vec = (1 / (4 * pi * epsilon_0 * epsilon_r)) * (q1 * q2 / |R_vec|^3) * R_vec

    Returns:
        tuple: A tuple containing:
            - str: A question asking for the force vector between two charges.
            - str: A step-by-step solution showing the calculation.
    """
    # 1. Parameterize the inputs with random values
    
    # Generate random integer charges between -15 and 15 uC (excluding 0)
    q1_val = random.choice([i for i in range(-15, 16) if i != 0])
    q2_val = random.choice([i for i in range(-15, 16) if i != 0])
    
    # Convert from microcoulombs to coulombs for calculation
    q1_calc = q1_val * 1e-6
    q2_calc = q2_val * 1e-6

    # Generate two unique 3D points with integer coordinates from -5 to 5
    p1 = np.array([random.randint(-5, 5) for _ in range(3)])
    p2 = np.array([random.randint(-5, 5) for _ in range(3)])
    while np.array_equal(p1, p2):
        p2 = np.array([random.randint(-5, 5) for _ in range(3)])

    # Choose a relative permittivity for the medium
    epsilon_r = round(random.uniform(1.0, 5.0), 2)
    
    # Standardize precision for final outputs
    precision = 3

    # 2. Perform the core calculation
    
    # Calculate the separation vector and its magnitude
    R_vec = p2 - p1
    R_mag = np.linalg.norm(R_vec)

    # Calculate the full permittivity of the medium
    epsilon = epsilon_r * EPSILON_0
    
    # Calculate the scalar part of Coulomb's Law
    scalar_coefficient = (q1_calc * q2_calc) / (4 * math.pi * epsilon * (R_mag**3))
    
    # Calculate the final force vector
    F_vec = scalar_coefficient * R_vec

    # 3. Generate the question and solution strings
    
    question = (
        f"Two point charges are located in a medium with a relative permittivity of {epsilon_r}.\n"
        f"Charge q1 = {q1_val} uC is at position P1 = {p1} m.\n"
        f"Charge q2 = {q2_val} uC is at position P2 = {p2} m.\n\n"
        f"Determine the electrostatic force vector F12 exerted by charge q1 on charge q2."
    )

    solution = (
        f"**Given:**\n"
        f"  - Charge q1 = {q1_val} uC = {q1_val}e-6 C at P1 = {p1} m\n"
        f"  - Charge q2 = {q2_val} uC = {q2_val}e-6 C at P2 = {p2} m\n"
        f"  - Relative permittivity epsilon_r = {epsilon_r}\n\n"
        f"**Constants:**\n"
        f"  - Permittivity of free space, epsilon_0 = {EPSILON_0:.4e} F/m\n\n"
        
        f"**Step 1:** Find the separation vector from q1 to q2.\n"
        f"  The vector R12 points from P1 to P2.\n"
        f"  R12 = P2 - P1 = {p2} - {p1} = {R_vec} m.\n\n"
        
        f"**Step 2:** Calculate the magnitude of the separation vector.\n"
        f"  |R12| = sqrt({R_vec[0]}^2 + {R_vec[1]}^2 + {R_vec[2]}^2)\n"
        f"  |R12| = sqrt({R_vec[0]**2} + {R_vec[1]**2} + {R_vec[2]**2}) = sqrt({R_mag**2:.2f}) = {R_mag:.{precision}f} m.\n\n"
        
        f"**Step 3:** Apply the vector form of Coulomb's Law.\n"
        f"  The formula for the force F12 is:\n"
        f"  F12 = (1 / (4 * pi * epsilon)) * (q1 * q2 / |R12|^3) * R12\n"
        f"  where epsilon = epsilon_r * epsilon_0.\n\n"
        
        f"  First, calculate the total permittivity:\n"
        f"  epsilon = {epsilon_r} * {EPSILON_0:.4e} = {epsilon:.4e} F/m.\n\n"
        
        f"  Now, substitute all values into the formula:\n"
        f"  F12 = (1 / (4 * pi * {epsilon:.4e})) * (({q1_calc:.2e}) * ({q2_calc:.2e}) / ({R_mag:.{precision}f})^3) * {R_vec}\n"
        f"  F12 = ({1/(4 * math.pi * epsilon):.3e}) * ({(q1_calc * q2_calc):.3e} / {R_mag**3:.3e}) * {R_vec}\n"
        f"  F12 = ({scalar_coefficient:.3e}) * {R_vec}\n"
        f"  F12 = <{F_vec[0]:.{precision}e}, {F_vec[1]:.{precision}e}, {F_vec[2]:.{precision}e}> N.\n\n"
        
        f"**Answer:**\n"
        f"  The electrostatic force vector exerted by q1 on q2 is "
        f"<{F_vec[0]:.{precision}e}, {F_vec[1]:.{precision}e}, {F_vec[2]:.{precision}e}> N."
    )

    return question, solution


# Template 2 (Intermediate)
def template_superposition_electric_field():
    """
    Superposition of Electric Fields (Intermediate)

    Scenario:
        This template requires calculating the net electric field at a target
        point due to two or three discrete point charges. It builds on the
        basic E-field calculation by applying the principle of superposition,
        which involves vector addition. The problem is constrained to 2D
        to keep the vector calculations straightforward.

    Core Equations:
        E_i = (1 / (4 * pi * epsilon_0)) * (q_i / |R_i|^3) * R_i
        E_net = E_1 + E_2 + ...

    Returns:
        tuple: A tuple containing:
            - str: A question asking for the net electric field at a point.
            - str: A step-by-step solution showing the calculation for each charge
                   and the final vector sum.
    """
    # 1. Parameterize the inputs with random values
    
    num_charges = random.choice([2, 3])
    charges_val = []
    positions = []
    
    # Use a set to ensure all points (charges and target) are unique
    used_points = set()
    
    for i in range(num_charges):
        # Generate random integer charges between -25 and 25 nC (excluding 0)
        q_val = random.choice([i for i in range(-25, 26) if i != 0])
        charges_val.append(q_val)

        # Generate a unique 2D point for each charge
        while True:
            pos = tuple(random.randint(-10, 10) for _ in range(2))
            if pos not in used_points:
                positions.append(np.array(pos))
                used_points.add(pos)
                break
    
    # Generate a unique 2D target point
    while True:
        target_pos = tuple(random.randint(-10, 10) for _ in range(2))
        if target_pos not in used_points:
            P_target = np.array(target_pos)
            used_points.add(target_pos)
            break
            
    # Convert from nanocoulombs to coulombs for calculation
    charges_calc = [q * 1e-9 for q in charges_val]

    # Standardize precision for final outputs
    precision = 2

    # 2. Perform the core calculation
    
    E_net = np.array([0.0, 0.0])
    solution_steps = ""
    
    for i in range(num_charges):
        q_calc = charges_calc[i]
        q_val = charges_val[i]
        pos = positions[i]
        
        # Calculate separation vector and its magnitude
        R_vec = P_target - pos
        R_mag = np.linalg.norm(R_vec)
        
        # Calculate the E-field vector from this charge
        scalar_part = q_calc / (4 * math.pi * EPSILON_0 * R_mag**3)
        E_vec = scalar_part * R_vec
        
        # Add to the net field
        E_net += E_vec
        
        # Build the solution string for this step
        solution_steps += (
            f"**Step {i+1}:** Calculate the Electric Field from q{i+1} ({q_val} nC).\n"
            f"  The separation vector from q{i+1} to the target point P is:\n"
            f"  R{i+1} = P - P{i+1} = {P_target} - {pos} = {R_vec} m.\n"
            f"  The magnitude is |R{i+1}| = sqrt({R_vec[0]**2} + {R_vec[1]**2}) = {R_mag:.{precision+1}f} m.\n"
            f"  The electric field E{i+1} is given by E = (k * q / |R|^3) * R, where k = 1/(4*pi*eps0).\n"
            f"  E{i+1} = ({1/(4*math.pi*EPSILON_0):.2e}) * ({q_calc:.2e} / {R_mag:.{precision+1}f}^3) * {R_vec}\n"
            f"  E{i+1} = <{E_vec[0]:.{precision}e}, {E_vec[1]:.{precision}e}> N/C.\n\n"
        )

    # 3. Generate the question and solution strings
    
    charge_descriptions = ""
    for i in range(num_charges):
        charge_descriptions += f"  - Charge q{i+1} = {charges_val[i]} nC is at position P{i+1} = {positions[i]} m.\n"

    question = (
        f"Several point charges are located in a vacuum on a 2D plane:\n"
        f"{charge_descriptions}"
        f"An observation point is located at P = {P_target} m.\n\n"
        f"Using the principle of superposition, determine the net electric field vector E_net at point P."
    )

    solution = (
        f"**Given:**\n"
        f"{charge_descriptions}"
        f"  - Observation point P = {P_target} m.\n"
        f"  - The medium is a vacuum (epsilon_r = 1.0).\n\n"
        f"**Principle of Superposition:**\n"
        f"  The total electric field at a point is the vector sum of the electric fields produced by each individual charge.\n"
        f"  E_net = E1 + E2{' + E3' if num_charges == 3 else ''}\n\n"

        f"{solution_steps}"
        f"**Step {num_charges+1}:** Sum the vectors to find the net electric field.\n"
        f"  E_net = E1 + E2{' + E3' if num_charges == 3 else ''}\n"
        f"  E_net = <E1_x + E2_x{' + E3_x' if num_charges == 3 else ''}, E1_y + E2_y{' + E3_y' if num_charges == 3 else ''}>\n"
        f"  E_net = <{E_net[0]:.{precision}e}, {E_net[1]:.{precision}e}> N/C.\n\n"

        f"**Answer:**\n"
        f"  The net electric field vector at point P is <{E_net[0]:.{precision}e}, {E_net[1]:.{precision}e}> N/C."
    )

    return question, solution


# Template 3 (Intermediate)
def template_gauss_law_symmetric():
    """
    Gauss's Law for Symmetric Charge Distributions (Intermediate)

    Scenario:
        This template tests the application of Gauss's Law to find the electric
        field (E) and electric flux density (D) for a highly symmetric charge
        distribution (an infinite line or an infinite sheet). The solution emphasizes
        the conceptual steps: choosing an appropriate Gaussian surface, applying
        Gauss's Law to find D, and then finding E.

    Core Equations:
        Gauss's Law: Integral(D . dS) = Q_enc
        For an infinite line: D = rho_l / (2 * pi * r)
        For an infinite sheet: D = rho_s / 2
        Relation: D = epsilon * E

    Returns:
        tuple: A tuple containing:
            - str: A question about the E and D fields from a charge distribution.
            - str: A detailed, step-by-step solution explaining the derivation.
    """
    # 1. Parameterize the inputs
    dist_type = random.choice(['line', 'sheet'])
    epsilon_r = round(random.uniform(1.0, 6.0), 2)
    epsilon = epsilon_r * EPSILON_0
    precision = 3

    # Initialize variables to be populated in the if/else block
    question = ""
    solution = ""

    # --- Logic for an Infinite Line Charge ---
    if dist_type == 'line':
        rho_l_val = random.choice([i for i in range(-50, 51) if i != 0])
        rho_l_calc = rho_l_val * 1e-9  # Convert nC/m to C/m
        
        r_dist = round(random.uniform(0.1, 2.5), 2)
        point = np.array([r_dist, 0, 0])

        # Core Calculation
        D_mag = rho_l_calc / (2 * math.pi * r_dist)
        D_vec = D_mag * np.array([1, 0, 0]) # Direction is radial (a_r)
        E_vec = D_vec / epsilon

        # Generate Question and Solution Strings
        question = (
            f"An infinite line of charge with a uniform density rho_l = {rho_l_val} nC/m is located "
            f"on the z-axis in a medium with relative permittivity epsilon_r = {epsilon_r}.\n\n"
            f"Using Gauss's Law, find the electric flux density vector (D) and the electric field "
            f"intensity vector (E) at the point P = {point} m."
        )

        solution = (
            f"**Given:**\n"
            f"  - Infinite line charge with rho_l = {rho_l_val} nC/m = {rho_l_calc:.2e} C/m.\n"
            f"  - Relative permittivity epsilon_r = {epsilon_r}.\n"
            f"  - Observation point P = {point} m.\n\n"

            f"**Step 1:** Choose a Gaussian Surface.\n"
            f"  For an infinite line charge, the electric field is purely radial. We choose a closed "
            f"cylindrical surface of radius r = {r_dist} m and length L, coaxial with the z-axis.\n\n"

            f"**Step 2:** Apply Gauss's Law.\n"
            f"  Gauss's Law states: Integral(D . dS) = Q_enc.\n"
            f"  - The flux D is perpendicular to the side surface and parallel to the top/bottom caps. "
            f"Therefore, flux only passes through the curved side surface.\n"
            f"  - Integral(D . dS) = D_r * (Area of side) = D_r * (2 * pi * r * L).\n"
            f"  - The charge enclosed is Q_enc = rho_l * L.\n"
            f"  - Equating them: D_r * (2 * pi * r * L) = rho_l * L.\n"
            f"  - Solving for D_r: D_r = rho_l / (2 * pi * r).\n\n"

            f"**Step 3:** Calculate the Electric Flux Density (D).\n"
            f"  The magnitude of D at r = {r_dist} m is:\n"
            f"  |D| = ({rho_l_calc:.2e}) / (2 * pi * {r_dist}) = {abs(D_mag):.{precision}e} C/m^2.\n"
            f"  At point P = {point} m, the direction is radial, which corresponds to the x-direction (a_x).\n"
            f"  Therefore, D = <{D_vec[0]:.{precision}e}, 0, 0> C/m^2.\n\n"

            f"**Step 4:** Calculate the Electric Field Intensity (E).\n"
            f"  E = D / epsilon, where epsilon = epsilon_r * epsilon_0.\n"
            f"  epsilon = {epsilon_r} * {EPSILON_0:.3e} = {epsilon:.3e} F/m.\n"
            f"  E = <{D_vec[0]:.{precision}e}, 0, 0> / {epsilon:.3e}\n"
            f"  E = <{E_vec[0]:.{precision}e}, 0, 0> V/m.\n\n"
        )

    # --- Logic for an Infinite Sheet of Charge ---
    else: # dist_type == 'sheet'
        rho_s_val = random.choice([i for i in range(-50, 51) if i != 0])
        rho_s_calc = rho_s_val * 1e-9  # Convert nC/m^2 to C/m^2
        
        z_dist = round(random.uniform(0.1, 2.5), 2)
        point = np.array([0, 0, z_dist])
        
        # Core Calculation
        D_mag = rho_s_calc / 2.0
        D_vec = D_mag * np.array([0, 0, 1]) # Direction is normal (a_z)
        E_vec = D_vec / epsilon

        # Generate Question and Solution Strings
        question = (
            f"An infinite sheet of charge with a uniform surface density rho_s = {rho_s_val} nC/m^2 is "
            f"located on the x-y plane (z=0) in a medium with relative permittivity epsilon_r = {epsilon_r}.\n\n"
            f"Using Gauss's Law, find the electric flux density vector (D) and the electric field "
            f"intensity vector (E) at the point P = {point} m."
        )

        solution = (
            f"**Given:**\n"
            f"  - Infinite sheet charge with rho_s = {rho_s_val} nC/m^2 = {rho_s_calc:.2e} C/m^2.\n"
            f"  - Relative permittivity epsilon_r = {epsilon_r}.\n"
            f"  - Observation point P = {point} m.\n\n"

            f"**Step 1:** Choose a Gaussian Surface.\n"
            f"  For an infinite sheet, the electric field is purely normal to the sheet. We choose a "
            f"small cylindrical 'pillbox' or a rectangular box of area A, piercing the sheet and centered at z=0.\n\n"

            f"**Step 2:** Apply Gauss's Law.\n"
            f"  Gauss's Law states: Integral(D . dS) = Q_enc.\n"
            f"  - The flux D is parallel to the sides of the pillbox, so flux only passes through the top and bottom surfaces.\n"
            f"  - Integral(D . dS) = D_z * (Top Area) + D_z * (Bottom Area) = D_z*A + D_z*A = 2*D_z*A.\n"
            f"  - The charge enclosed is Q_enc = rho_s * A.\n"
            f"  - Equating them: 2 * D_z * A = rho_s * A.\n"
            f"  - Solving for D_z: D_z = rho_s / 2.\n\n"

            f"**Step 3:** Calculate the Electric Flux Density (D).\n"
            f"  The magnitude of D is independent of the distance from the sheet:\n"
            f"  |D| = |{rho_s_calc:.2e}| / 2 = {abs(D_mag):.{precision}e} C/m^2.\n"
            f"  At P = {point} m (where z > 0), the direction is normal to the sheet (a_z).\n"
            f"  Therefore, D = <0, 0, {D_vec[2]:.{precision}e}> C/m^2.\n\n"

            f"**Step 4:** Calculate the Electric Field Intensity (E).\n"
            f"  E = D / epsilon, where epsilon = epsilon_r * epsilon_0.\n"
            f"  epsilon = {epsilon_r} * {EPSILON_0:.3e} = {epsilon:.3e} F/m.\n"
            f"  E = <0, 0, {D_vec[2]:.{precision}e}> / {epsilon:.3e}\n"
            f"  E = <0, 0, {E_vec[2]:.{precision}e}> V/m.\n\n"
        )

    # Final Summary for both cases
    final_answer = (
        f"**Answer:**\n"
        f"  The electric flux density is D = <{D_vec[0]:.{precision}e}, {D_vec[1]:.{precision}e}, {D_vec[2]:.{precision}e}> C/m^2.\n"
        f"  The electric field intensity is E = <{E_vec[0]:.{precision}e}, {E_vec[1]:.{precision}e}, {E_vec[2]:.{precision}e}> V/m."
    )
    
    solution += final_answer
    
    return question, solution


# Template 4 (Advanced)
def template_coaxial_capacitance():
    """
    Capacitance of a Coaxial Cable (Advanced)

    Scenario:
        This template tests the ability to derive and calculate the capacitance
        per unit length of a coaxial cable. The solution requires a logical,
        multi-step derivation: first, finding the electric field using the
        premise of Gauss's Law; second, integrating the electric field to find
        the potential difference between the conductors; and finally, using the
        fundamental definition of capacitance (C = Q/V).

    Core Equations:
        1. E = rho_l / (2 * pi * epsilon * r)
        2. V_ab = Integral(E . dl) = (rho_l / (2 * pi * epsilon)) * ln(b/a)
        3. C' = rho_l / V_ab = (2 * pi * epsilon) / ln(b/a)

    Returns:
        tuple: A tuple containing:
            - str: A question asking to derive and calculate the capacitance per unit length.
            - str: A step-by-step solution showing the full derivation and calculation.
    """
    # 1. Parameterize the inputs
    
    # Inner radius in mm, ensuring it's not too small
    a_mm = round(random.uniform(0.5, 2.5), 2)
    
    # Outer radius is 2 to 5 times larger than inner radius
    b_mm = round(a_mm * random.uniform(2, 5), 2)
    
    # Convert radii to meters for calculation
    a_m = a_mm * 1e-3
    b_m = b_mm * 1e-3
    
    # Relative permittivity of a common dielectric like polyethylene or teflon
    epsilon_r = round(random.uniform(2.0, 4.0), 2)
    
    # Standardize precision for final outputs
    precision = 3

    # 2. Perform the core calculation
    
    # Total permittivity of the dielectric
    epsilon = epsilon_r * EPSILON_0
    
    # Capacitance per unit length in F/m
    C_prime = (2 * math.pi * epsilon) / math.log(b_m / a_m)
    
    # Convert to a more readable unit, pF/m
    C_prime_pF = C_prime * 1e12

    # 3. Generate the question and solution strings
    
    question = (
        f"A coaxial cable consists of an inner conductor of radius a = {a_mm} mm and an "
        f"outer conductor of radius b = {b_mm} mm.\n"
        f"The space between the conductors is filled with a dielectric material with a "
        f"relative permittivity of epsilon_r = {epsilon_r}.\n\n"
        f"Derive the general expression for the capacitance per unit length (C') and then "
        f"calculate its numerical value for this cable."
    )

    solution = (
        f"**Given:**\n"
        f"  - Inner radius, a = {a_mm} mm = {a_m:.2e} m\n"
        f"  - Outer radius, b = {b_mm} mm = {b_m:.2e} m\n"
        f"  - Relative permittivity, epsilon_r = {epsilon_r}\n\n"

        f"**Derivation:**\n"
        f"The derivation involves three main steps: finding the electric field E, finding the "
        f"potential difference V, and then finding the capacitance C' = rho_l / V.\n\n"

        f"**Step 1:** Find the Electric Field (E) between the conductors.\n"
        f"  Assume a line charge of +rho_l (C/m) on the inner conductor and -rho_l on the outer. "
        f"By applying Gauss's Law to a cylindrical surface with radius 'r' (where a < r < b), "
        f"the electric field is found to be purely radial:\n"
        f"  E = (rho_l / (2 * pi * epsilon * r)) * a_r\n"
        f"  where epsilon = epsilon_r * epsilon_0.\n\n"

        f"**Step 2:** Find the Potential Difference (V_ab) between the conductors.\n"
        f"  The potential difference is the work done to move a unit charge from the outer "
        f"conductor (at r=b) to the inner conductor (at r=a) against the E field.\n"
        f"  V_ab = -Integral(from b to a, E . dl)\n"
        f"  Since the path is radial, dl = dr * a_r. The integral becomes:\n"
        f"  V_ab = -Integral(from b to a, [rho_l / (2 * pi * epsilon * r)] dr)\n"
        f"  V_ab = -[rho_l / (2 * pi * epsilon)] * [ln(r)](from b to a)\n"
        f"  V_ab = -[rho_l / (2 * pi * epsilon)] * (ln(a) - ln(b))\n"
        f"  Using log properties, ln(a) - ln(b) = -ln(b/a), this simplifies to:\n"
        f"  V_ab = (rho_l / (2 * pi * epsilon)) * ln(b/a)\n\n"

        f"**Step 3:** Derive the Capacitance per Unit Length (C').\n"
        f"  Capacitance per unit length is defined as the charge per unit length (rho_l) "
        f"divided by the potential difference (V_ab).\n"
        f"  C' = rho_l / V_ab\n"
        f"  C' = rho_l / [(rho_l / (2 * pi * epsilon)) * ln(b/a)]\n"
        f"  The rho_l terms cancel, leaving the final expression:\n"
        f"  C' = (2 * pi * epsilon) / ln(b/a)\n\n"

        f"**Step 4:** Calculate the final value.\n"
        f"  First, find the permittivity of the dielectric:\n"
        f"  epsilon = {epsilon_r} * {EPSILON_0:.4e} = {epsilon:.4e} F/m.\n"
        f"  Now, plug the values into the derived formula:\n"
        f"  C' = (2 * pi * {epsilon:.4e}) / ln({b_mm} / {a_mm})\n"
        f"  C' = (2 * pi * {epsilon:.4e}) / ln({b_m/a_m:.{precision}f})\n"
        f"  C' = (2 * pi * {epsilon:.4e}) / {math.log(b_m/a_m):.{precision}f}\n"
        f"  C' = {C_prime:.{precision}e} F/m.\n\n"

        f"**Answer:**\n"
        f"  The derived expression for capacitance per unit length is C' = (2 * pi * epsilon) / ln(b/a).\n"
        f"  The numerical value is {C_prime:.{precision}e} F/m, which is equivalent to "
        f"**{C_prime_pF:.{precision-1}f} pF/m**."
    )

    return question, solution


def main():
    """
    Generate numerous instances of each electrostatics template 
    with different random seeds and write the results to a JSONL file.
    """
    import json
    import os

    # Define the output path (Modify this path according to where you are running the code from)
    output_file = "testset/electrical_engineering/electromagnetics_and_waves/electrostatics.jsonl"

    # Create the directory if it doesn't exist
    os.makedirs(os.path.dirname(output_file), exist_ok=True)

    # List of template functions with their ID and level
    templates = [
        (template_coulombs_law, "coulombs_law_two_charges", "Easy"),
        (template_superposition_electric_field, "superposition_electric_field", "Intermediate"),
        (template_gauss_law_symmetric, "gauss_law_symmetric_charge", "Intermediate"),
        (template_coaxial_capacitance, "coaxial_cable_capacitance", "Advanced"),
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
