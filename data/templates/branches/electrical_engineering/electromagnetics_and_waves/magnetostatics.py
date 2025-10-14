import random
import math


# Template 1 (Easy)
def template_lorentz_force():
    """
    Lorentz Force on a Moving Charge

    Scenario:
        This template tests the fundamental calculation of the magnetic force on a moving
        point charge in a uniform magnetic field. It requires unit conversion and the
        correct application of the vector cross product, which is a core skill in
        electromagnetics.

    Core Equation:
        F_m = q * (u x B)

    Returns:
        tuple: A tuple containing:
            - str: A question asking for the force vector and magnitude.
            - str: A step-by-step solution.
    """
    # 1. Parameterize the inputs with random values
    
    # Generate charge in microcoulombs (uC)
    q_uC = round(random.uniform(1.0, 100.0), 2)
    if random.choice([True, False]):
        q_uC *= -1
    
    # Generate velocity vector components in m/s
    u_vec = [random.randint(-50, 50) for _ in range(3)]
    
    # Generate magnetic field vector components in millitesla (mT)
    B_vec_mT = [round(random.uniform(-200.0, 200.0), 2) for _ in range(3)]

    # Standardize precision for final outputs
    precision = 3

    # 2. Perform the core calculation
    
    # Convert units for calculation
    q_C = q_uC * 1e-6
    B_vec_T = [b * 1e-3 for b in B_vec_mT]
    
    # Calculate the cross product: u x B
    cross_product_x = u_vec[1] * B_vec_T[2] - u_vec[2] * B_vec_T[1]
    cross_product_y = u_vec[2] * B_vec_T[0] - u_vec[0] * B_vec_T[2]
    cross_product_z = u_vec[0] * B_vec_T[1] - u_vec[1] * B_vec_T[0]
    
    # Calculate the force vector: F = q * (u x B)
    force_x = q_C * cross_product_x
    force_y = q_C * cross_product_y
    force_z = q_C * cross_product_z
    
    # Calculate the magnitude of the force
    force_magnitude = math.sqrt(force_x**2 + force_y**2 + force_z**2)
    
    # 3. Generate the question and solution strings
    
    # Helper function to format vectors for display
    def format_vector(v, units=""):
        return f"({v[0]} x_hat + {v[1]} y_hat + {v[2]} z_hat) {units}".strip()

    question = (
        f"A point charge of {q_uC} uC has a velocity of u = {format_vector(u_vec, 'm/s')} "
        f"in a uniform magnetic field described by B = {format_vector(B_vec_mT, 'mT')}.\n\n"
        f"Determine the magnetic force vector F_m acting on the charge and its magnitude."
    )

    solution = (
        f"**Given:**\n"
        f"  - Charge (q): {q_uC} uC\n"
        f"  - Velocity (u): {format_vector(u_vec, 'm/s')}\n"
        f"  - Magnetic Field (B): {format_vector(B_vec_mT, 'mT')}\n\n"
        
        f"**Step 1:** Convert Units to SI\n"
        f"  First, we convert the given values to standard SI units for the calculation.\n"
        f"  - Charge in Coulombs: q = {q_uC} * 1e-6 = {q_C:.2e} C\n"
        f"  - Magnetic Field in Tesla: B = ({B_vec_T[0]:.2e} x_hat + {B_vec_T[1]:.2e} y_hat + {B_vec_T[2]:.2e} z_hat) T\n\n"
        
        f"**Step 2:** Calculate the Cross Product (u x B)\n"
        f"  The force is determined by the formula F_m = q * (u x B). We start by calculating the cross product.\n"
        f"  u x B = [ (u_y * B_z - u_z * B_y) x_hat + (u_z * B_x - u_x * B_z) y_hat + (u_x * B_y - u_y * B_x) z_hat ]\n"
        f"  (u x B)_x = ({u_vec[1]}) * ({B_vec_T[2]:.2e}) - ({u_vec[2]}) * ({B_vec_T[1]:.2e}) = {cross_product_x:.4f}\n"
        f"  (u x B)_y = ({u_vec[2]}) * ({B_vec_T[0]:.2e}) - ({u_vec[0]}) * ({B_vec_T[2]:.2e}) = {cross_product_y:.4f}\n"
        f"  (u x B)_z = ({u_vec[0]}) * ({B_vec_T[1]:.2e}) - ({u_vec[1]}) * ({B_vec_T[0]:.2e}) = {cross_product_z:.4f}\n"
        f"  So, u x B = ({round(cross_product_x, precision)} x_hat + {round(cross_product_y, precision)} y_hat + {round(cross_product_z, precision)} z_hat) T*m/s\n\n"

        f"**Step 3:** Calculate the Force Vector (F_m)\n"
        f"  Now, multiply the cross product by the charge q.\n"
        f"  F_m = ({q_C:.2e} C) * ({round(cross_product_x, precision)} x_hat + {round(cross_product_y, precision)} y_hat + {round(cross_product_z, precision)} z_hat)\n"
        f"  F_m = ({force_x:.{precision}e} x_hat + {force_y:.{precision}e} y_hat + {force_z:.{precision}e} z_hat) N\n\n"
        
        f"**Step 4:** Calculate the Magnitude of the Force\n"
        f"  The magnitude is the square root of the sum of the squares of the components.\n"
        f"  |F_m| = sqrt( ({force_x:.2e})^2 + ({force_y:.2e})^2 + ({force_z:.2e})^2 )\n"
        f"  |F_m| = {force_magnitude:.{precision}e} N\n\n"
        
        f"**Answer:**\n"
        f"  The magnetic force vector is F_m = ({force_x:.{precision}e} x_hat + {force_y:.{precision}e} y_hat + {force_z:.{precision}e} z_hat) N.\n"
        f"  The magnitude of the force is |F_m| = {force_magnitude:.{precision}e} N."
    )

    return question, solution


def main():
    """
    Generate numerous instances of each magnetostatics template 
    with different random seeds and write the results to a JSONL file.
    """
    import json
    import os

    # Define the output path (Modify this path according to where you are running the code from)
    output_file = "testset/electrical_engineering/electromagnetics_and_waves/magnetostatics.jsonl"

    # Create the directory if it doesn't exist
    os.makedirs(os.path.dirname(output_file), exist_ok=True)

    # List of template functions with their ID and level
    templates = [
        (template_lorentz_force, "lorentz_force", "Easy"),
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
                "area": "magnetostatics",
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
