import random
import math
from decimal import Decimal, ROUND_HALF_UP
from data.templates.branches._emission import signed_term, joined_terms


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

    Trace integrity (Layer 0, 2026-09-23):
        B is stated in the question in mT to 2 dp, so in tesla it has at most
        five significant figures. The trace displayed it at `.2e` (three) and
        consumed the full value, so the three cross-product lines did not
        close on their printed operands. Binding B through the `.2e` display
        instead would move |F_m| by up to 7.8% (measured over 20,000 draws;
        the cross product cancels), so the tesla display is lengthened to
        `.4e`, which is exact for every sampled B, and B is bound through it
        (D-016 part 2, D-037). An integer velocity times a 5-dp field is
        exact at 5 dp; displayed at `.4f` it sat on a half-way tie for 22%
        of instances, so the cross product is displayed at `.5f` everywhere
        it appears and bound through that display. q is stated to 2 dp in
        uC (four significant figures) and was displayed at `.2e`; it is
        displayed at `.3e`, exact for every sampled q, and bound through it,
        so Step 3 closes for a reader too. The question text and the gold
        answer are unchanged.

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
    # q is stated to 2 dp in uC (four significant figures); `.3e` is exact for
    # every sampled q, so q is displayed at that precision and bound through
    # it, and Step 3 closes for a reader (D-016 part 2, D-037).
    q_C = _as_printed(q_uC * 1e-6, '.3e')
    # B is stated in the question in mT to 2 dp, so in tesla it has at most
    # five significant figures: display it at that precision (`.4e`) and bind
    # it through the display, so the value a reader recovers from the text is
    # the value the cross product consumes (D-016 part 2, D-037).
    B_vec_T = [_as_printed(b * 1e-3, '.4e') for b in B_vec_mT]
    
    # Calculate the cross product: u x B
    # An integer velocity times a 5-dp field is exact at 5 dp, so the cross
    # product is displayed at `.5f` everywhere it appears and bound through it.
    cross_product_x = _as_printed(u_vec[1] * B_vec_T[2] - u_vec[2] * B_vec_T[1], '.5f')
    cross_product_y = _as_printed(u_vec[2] * B_vec_T[0] - u_vec[0] * B_vec_T[2], '.5f')
    cross_product_z = _as_printed(u_vec[0] * B_vec_T[1] - u_vec[1] * B_vec_T[0], '.5f')
    
    # Calculate the force vector: F = q * (u x B)
    force_x = q_C * cross_product_x
    force_y = q_C * cross_product_y
    force_z = q_C * cross_product_z
    
    # Calculate the magnitude of the force
    force_magnitude = math.sqrt(force_x**2 + force_y**2 + force_z**2)
    
    # 3. Generate the question and solution strings
    
    # Helper function to format vectors for display
    def format_vector(v, units=""):
        return f"({v[0]} x_hat {signed_term(v[1], 'y_hat')} {signed_term(v[2], 'z_hat')}) {units}".strip()

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
        f"  - Charge in Coulombs: q = {q_uC} * 1e-6 = {q_C:.3e} C\n"
        f"  - Magnetic Field in Tesla: B = ({B_vec_T[0]:.4e} x_hat {signed_term(B_vec_T[1], 'y_hat', '{:.4e}'.format)} {signed_term(B_vec_T[2], 'z_hat', '{:.4e}'.format)}) T\n\n"
        
        f"**Step 2:** Calculate the Cross Product (u x B)\n"
        f"  The force is determined by the formula F_m = q * (u x B). We start by calculating the cross product.\n"
        f"  u x B = [ (u_y * B_z - u_z * B_y) x_hat + (u_z * B_x - u_x * B_z) y_hat + (u_x * B_y - u_y * B_x) z_hat ]\n"
        f"  (u x B)_x = ({u_vec[1]}) * ({B_vec_T[2]:.4e}) - ({u_vec[2]}) * ({B_vec_T[1]:.4e}) = {cross_product_x:.5f}\n"
        f"  (u x B)_y = ({u_vec[2]}) * ({B_vec_T[0]:.4e}) - ({u_vec[0]}) * ({B_vec_T[2]:.4e}) = {cross_product_y:.5f}\n"
        f"  (u x B)_z = ({u_vec[0]}) * ({B_vec_T[1]:.4e}) - ({u_vec[1]}) * ({B_vec_T[0]:.4e}) = {cross_product_z:.5f}\n"
        f"  So, u x B = ({cross_product_x:.5f} x_hat {signed_term(cross_product_y, 'y_hat', '{:.5f}'.format)} {signed_term(cross_product_z, 'z_hat', '{:.5f}'.format)}) T*m/s\n\n"

        f"**Step 3:** Calculate the Force Vector (F_m)\n"
        f"  Now, multiply the cross product by the charge q.\n"
        f"  F_m = ({q_C:.3e} C) * ({cross_product_x:.5f} x_hat {signed_term(cross_product_y, 'y_hat', '{:.5f}'.format)} {signed_term(cross_product_z, 'z_hat', '{:.5f}'.format)})\n"
        f"  F_m = ({force_x:.{precision}e} x_hat {signed_term(force_y, 'y_hat', ('{:.' + str(precision) + 'e}').format)} {signed_term(force_z, 'z_hat', ('{:.' + str(precision) + 'e}').format)}) N\n\n"
        
        f"**Step 4:** Calculate the Magnitude of the Force\n"
        f"  The magnitude is the square root of the sum of the squares of the components.\n"
        f"  |F_m| = sqrt( ({force_x:.2e})^2 + ({force_y:.2e})^2 + ({force_z:.2e})^2 )\n"
        f"  |F_m| = {force_magnitude:.{precision}e} N\n\n"
        
        f"**Answer:**\n"
        f"  The magnetic force vector is F_m = ({force_x:.{precision}e} x_hat {signed_term(force_y, 'y_hat', ('{:.' + str(precision) + 'e}').format)} {signed_term(force_z, 'z_hat', ('{:.' + str(precision) + 'e}').format)}) N.\n"
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
