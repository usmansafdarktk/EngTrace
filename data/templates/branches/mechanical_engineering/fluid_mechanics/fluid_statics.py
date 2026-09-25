import random
import math
from decimal import Decimal, ROUND_HALF_UP
from data.templates.branches.mechanical_engineering.constants import GRAVITY, ATMOSPHERIC_PRESSURE_KPA, FLUID_DENSITIES, MATERIAL_DENSITIES, OBJECT_SHAPES, OBJECT_MATERIALS, PIPE_FLUIDS, MANOMETER_FLUIDS


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
def template_hydrostatic_pressure_at_depth():
    """
    Fluid Statics: Hydrostatic Pressure at Depth

    Scenario:
        This template generates a fundamental problem to calculate the absolute
        pressure at a specified depth within a fluid. It considers both the
        pressure exerted by the fluid column and the pressure at the free surface.

    Core Equation:
        p_absolute = p_surface_absolute + (rho * g * h)

    Trace integrity (Layer 0, 2026-09-23):
        The pascal chain (surface pressure, rho*g*h, their sum) is displayed
        at 1 dp and each value is bound through that display before it is
        consumed (D-016 part 2). The answer is quoted at 4 dp in kPa, which is
        exactly the 1-dp pascal value shifted three places: at 3 dp, Step 4
        dropped a digit across the unit conversion and "359181.5 Pa / 1000 =
        359.181 kPa" sat on a half-way tie on every instance whose pascal
        value ended in .5 (D-016 part 1). A depth whose exact rho*g*h sits on
        a 1-dp tie is redrawn (D-016 part 3).

    Returns:
        tuple: A tuple containing:
            - str: A question asking for the absolute pressure at a certain depth.
            - str: A step-by-step solution to the problem.
    """
    # 1. Parameterize the inputs with random values

    # Randomly select a fluid and its properties
    fluid_name, density_rho = random.choice(list(FLUID_DENSITIES.items()))

    # Randomize depth in meters (Restricted to realistic tank depths). The
    # hydrostatic term rho*g*h is exact at 4 dp and is displayed at 1 dp; a
    # depth whose exact term sits on a half-way 1-dp tie has no defensible
    # rounding and is redrawn (D-016 part 3).
    for _attempt in range(200):
        depth_h = round(random.uniform(2.0, 25.0), 2)
        pressure_increase_exact = density_rho * GRAVITY * depth_h
        if not _is_display_tie(pressure_increase_exact, 1):
            break
    else:
        raise RuntimeError("hydrostatic_pressure_at_depth: no closing sample in 200 draws")

    # Randomly decide if the surface pressure is atmospheric or a specified gauge pressure
    is_surface_atmospheric = random.choice([True, False])
    
    # Define Atmospheric Pressure
    P_atm_kpa = 101.325

    if is_surface_atmospheric:
        # Surface is open to atmosphere
        surface_pressure_input_kpa = P_atm_kpa
        surface_pressure_abs_kpa = P_atm_kpa
        pressure_type_desc = "is exposed to the atmosphere"
        surface_calc_note = "Since the surface is open to the atmosphere, the surface absolute pressure is just P_atm."
    else:
        # Surface is pressurized (Gauge Pressure given)
        surface_gauge_kpa = round(random.uniform(10.0, 200.0), 2)
        surface_pressure_input_kpa = surface_gauge_kpa
        
        # FIX: Absolute Surface Pressure = Gauge + Atm
        surface_pressure_abs_kpa = surface_gauge_kpa + P_atm_kpa
        
        pressure_type_desc = f"is maintained at a gauge pressure of {surface_gauge_kpa} kPa"
        surface_calc_note = (
            f"The problem gives gauge pressure. To get absolute pressure at the surface, we add atmospheric pressure:\n"
            f"P_surface_abs = P_gauge + P_atm = {surface_gauge_kpa} + {P_atm_kpa} = {surface_pressure_abs_kpa:.3f} kPa"
        )

    # Standardize precision for final outputs. The pascal chain is carried and
    # displayed at 1 dp, and 1 dp in Pa IS 4 dp in kPa, so the answer is quoted
    # at 4 dp and the conversion in Step 4 is exact; at 3 dp it dropped a digit
    # across the conversion and sat on a tie whenever the Pa value ended in .5
    # (D-016 part 1).
    precision = 4

    # 2. Perform the core calculations for the solution. Every value that is
    # displayed and then consumed is bound through its display (D-016 part 2).

    # Step A: Convert absolute surface pressure to Pascals
    surface_pressure_abs_kpa = _as_printed(surface_pressure_abs_kpa, '.3f')
    surface_pressure_pa = _as_printed(surface_pressure_abs_kpa * 1000, '.1f')

    # Step B: Calculate the pressure increase due to the fluid column (Hydrostatic Pressure)
    # p_hydro = rho * g * h
    pressure_increase_pa = _as_printed(pressure_increase_exact, '.1f')

    # Step C: Calculate the final absolute pressure in Pascals
    absolute_pressure_pa = _as_printed(surface_pressure_pa + pressure_increase_pa, '.1f')

    # Step D: Convert the final answer back to kilopascals (kPa)
    absolute_pressure_kpa = _as_printed(absolute_pressure_pa / 1000, f'.{precision}f')

    # 3. Generate the question and solution strings

    question = (
        f"An object is located {depth_h} m below the surface of a tank containing "
        f"{fluid_name.lower()}. The pressure at the free surface {pressure_type_desc}. "
        f"Assuming the density of {fluid_name.lower()} is {density_rho} kg/m^3 and standard atmospheric pressure is {P_atm_kpa} kPa, "
        f"what is the absolute pressure at this depth?"
    )

    solution = (
        f"**Given:**\n"
        f"Fluid: {fluid_name} (Density rho = {density_rho} kg/m^3)\n"
        f"Depth (h): {depth_h} m\n"
        f"Surface Condition: {pressure_type_desc}\n"
        f"Standard Atmospheric Pressure (P_atm): {P_atm_kpa} kPa\n"
        f"Acceleration due to Gravity (g): {GRAVITY} m/s^2\n\n"

        f"**Step 1:** Determine the Absolute Pressure at the Surface\n"
        f"{surface_calc_note}\n"
        f"Convert to Pascals: P_surface_abs = {surface_pressure_abs_kpa:.3f} kPa * 1000 = {surface_pressure_pa:.1f} Pa\n\n"

        f"**Step 2:** Calculate the Hydrostatic Pressure Increase\n"
        f"The pressure exerted by the fluid column is calculated using: P_hydro = rho * g * h.\n"
        f"P_hydro = {density_rho} kg/m^3 * {GRAVITY} m/s^2 * {depth_h} m\n"
        f"P_hydro = {pressure_increase_pa:.1f} Pa\n"
        f"\n\n"

        f"**Step 3:** Calculate Total Absolute Pressure at Depth\n"
        f"The absolute pressure at depth is the sum of the absolute surface pressure and the hydrostatic pressure.\n"
        f"P_absolute = P_surface_abs + P_hydro\n"
        f"P_absolute = {surface_pressure_pa:.1f} Pa + {pressure_increase_pa:.1f} Pa\n"
        f"P_absolute = {absolute_pressure_pa:.1f} Pa\n\n"

        f"**Step 4:** Convert Final Answer to kPa\n"
        f"P_absolute = {absolute_pressure_pa:.1f} Pa / 1000 = {round(absolute_pressure_kpa, precision)} kPa\n\n"

        f"**Answer:**\n"
        f"The absolute pressure at a depth of {depth_h} m is **{round(absolute_pressure_kpa, precision)} kPa**."
    )

    return question, solution


# Template 2 (Easy)
def template_basic_buoyant_force():
    """
    Fluid Statics: Basic Buoyant Force on a Submerged Object

    Scenario:
        This template applies Archimedes' principle to calculate the buoyant
        force on a fully submerged object with a known volume. It reinforces
        the direct relationship between buoyant force and the weight of the
        displaced fluid.

    Core Equation:
        F_buoyant = rho_fluid * g * V_displaced

    Screen pass 1 (2026-09-23):
        All three judges flagged nonsensical objects from the independent
        draw of a material and a shape ("solid bronze wooden block", a
        "titanium boat hull" fully submerged), the flag the December 2025
        Tribunal raised. Measured over 503 seeds: 31.8% of draws named a
        shape that carries its own material word, 26.6% contradicted the
        drawn material ("uranium concrete piling", "ice boulder", "balsa
        wood metal cylinder"), and 25.0% named a hollow or floating body as
        a solid, fully submerged one. The shape is now drawn from the eight
        geometric solids in OBJECT_SHAPES (see the sampling comment); the
        materials, fluids, volume range, physics and steps are unchanged.
        Item-pool effect, accepted: the shorter choice list re-indexes the
        draw, so the question, and with it the material, the volume and
        the answer, moved on 96.2% of 501 seeds (0-500); the fluid is drawn
        before the shape and is unchanged.

    Layer 2 fix (2026-09-25):
        One expert of the branch rejected the template: material and fluid
        are drawn independently, so a body lighter than the fluid was
        called "solid ... fully submerged" with nothing holding it there -
        aluminum in mercury (seed 2101, which would float carrying its own
        19 kN weight, not 96.5 kN), cork in ethanol (2104), balsa in engine
        oil (2105). Of the two remedies the expert offered, the question
        now states the restraint: the body "is held fully submerged ... by
        a rigid clamp". Chosen over rejecting rho_object < rho_fluid draws
        because (a) Archimedes' principle gives rho_fluid * g * V for any
        body held fully under the surface, so the restraint makes every
        drawn pair physically exact without touching the arithmetic, and
        (b) the rejection would need a density for each of the 30
        descriptive OBJECT_MATERIALS names ("foam", "plastic", "fiberglass"
        and "ceramic" have no MATERIAL_DENSITIES row) and would remove
        every light material in every fluid and every material in mercury,
        a large and material-selective cut of the pool. Wording only: the
        gold answer is unchanged on every seed.

    Returns:
        tuple: A tuple containing:
            - str: A question asking for the buoyant force on an object.
            - str: A step-by-step solution to the problem.
    """
    # 1. Parameterize the inputs with random values

    # Randomly select a fluid and its properties
    fluid_name, density_rho = random.choice(list(FLUID_DENSITIES.items()))

    # Randomly select descriptive properties for the object.
    #
    # Screen pass 1 (2026-09-23): only the geometric solids in OBJECT_SHAPES
    # are drawn, in table order (D-031). The other entries either carry a
    # material word of their own ("wooden block", "metal rod", "concrete
    # piling", "stone") that contradicts the drawn material, or name a hollow
    # or floating body ("boat hull", "buoy", "storage tank") that is neither
    # solid nor fully submerged.
    geometric_solids = ("sphere", "cube", "irregular block", "cylinder",
                        "rectangular prism", "cone", "pyramid", "ball")
    shapes = [s for s in OBJECT_SHAPES if s in geometric_solids]
    assert len(shapes) == len(geometric_solids), "OBJECT_SHAPES lost a geometric solid"
    shape = random.choice(shapes)
    material = random.choice(OBJECT_MATERIALS)

    # Randomize the object's volume in cubic meters
    object_volume = round(random.uniform(0.05, 2.5), 3)

    # Standardize precision for final outputs
    precision = 3

    # 2. Perform the core calculations for the solution

    # Step A: Calculate the buoyant force in Newtons (N)
    # Since the object is fully submerged, the displaced volume equals the object's volume.
    buoyant_force_n = density_rho * GRAVITY * object_volume

    # Step B: Convert to kilonewtons (kN) if the value is large, for better readability
    buoyant_force_kn = buoyant_force_n / 1000

    # 3. Generate the question and solution strings

    question = (
        f"A solid {material} {shape} with a total volume of {object_volume} m^3 is "
        f"held fully submerged in a tank of {fluid_name.lower()} by a rigid clamp. "
        f"Given that the density of {fluid_name.lower()} is {density_rho} kg/m^3, "
        f"calculate the buoyant force acting on the {shape}."
    )

    solution = (
        f"**Given:**\n"
        f"Object Volume (V): {object_volume} m^3\n"
        f"Fluid: {fluid_name}\n"
        f"Density of Fluid (rho): {density_rho} kg/m^3\n"
        f"Acceleration due to Gravity (g): {GRAVITY} m/s^2\n\n"

        f"**Step 1:** Identify the Principle\n"
        f"According to Archimedes' principle, the buoyant force (F_buoyant) on a submerged "
        f"object is equal to the weight of the fluid it displaces.\n"
        f"Since the object is fully submerged, the volume of displaced fluid is equal to the "
        f"volume of the object itself.\n\n"

        f"**Step 2:** Apply the Buoyant Force Formula\n"
        f"The formula for buoyant force is: F_buoyant = rho_fluid * g * V_displaced\n"
        f"F_buoyant = {density_rho} kg/m^3 * {GRAVITY} m/s^2 * {object_volume} m^3\n"
        f"F_buoyant = {round(buoyant_force_n, precision)} N\n\n"

        f"**Step 3:** Convert the Result to Kilonewtons (kN) (Optional)\n"
        f"For larger force values, it is common to express the result in kilonewtons.\n"
        f"F_buoyant = {round(buoyant_force_n, precision)} N / 1000 = {round(buoyant_force_kn, precision)} kN\n\n"

        f"**Answer:**\n"
        f"The buoyant force acting on the {shape} is {round(buoyant_force_n, precision)} N, "
        f"which is equivalent to {round(buoyant_force_kn, precision)} kN."
    )

    return question, solution


# Template 3 (Intermediate)
def template_utube_manometer():
    """
    Fluid Statics: U-Tube Manometer Pressure Measurement

    Scenario:
        This template generates a classic problem involving a U-tube manometer
        used to measure the gauge pressure of a fluid inside a pipe. The solution
        requires applying the principle of hydrostatic equilibrium by balancing
        the pressure contributions from different fluid columns.

    Core Equation:
        P_gauge = (rho_2 * g * h2) - (rho_1 * g * h1)

    Trace integrity (Layer 0, 2026-09-23):
        The two column terms and the gauge pressure in pascals are displayed
        at 4 dp instead of 2: an integer density times g = 9.81 times a 2-dp
        height is exact at 4 dp, so the tie that "1030 * 9.81 * 0.35 =
        3536.51 Pa" sat on is removed by construction (D-037), and Step 5
        subtracts the displayed terms rather than unrounded ones (D-016 part
        2). A gas density carries more decimals; a draw whose exact term sits
        on a 4-dp tie, or whose gauge pressure sits on the 3-dp kPa tie of
        Step 6, is redrawn (D-016 part 3). The kPa answer keeps its 3 dp.

    Layer 2 fix (2026-09-25):
        All three experts of the branch rejected the template because
        MANOMETER_FLUIDS carried zinc (6570 kg/m^3, the density of the
        molten metal, solid below 420 C) and it was the manometer liquid in
        seeds 2102-2104 against engine oil, gasoline and isopropyl alcohol
        at ambient conditions (16.182, 35.099, 35.35 kPa). The table's
        "Molten Metals" group - gallium, tin and zinc - is deleted in
        constants.py (see the note there); nothing in this function
        changes. The hydrostatic balance and arithmetic were confirmed
        correct by all three experts. Deleting rows re-indexes
        random.choice over the table (D-031), so the item-pool effect is
        wholesale; it is measured in the Layer 2 fix report.

    Returns:
        tuple: A tuple containing:
            - str: A question asking for the gauge pressure in a pipe.
            - str: A step-by-step solution to the problem.
    """
    # 1. Parameterize the inputs with a loop to ensure positive gauge pressure
    # This prevents physically confusing scenarios where the description says the open arm 
    # level is "higher" but the math yields suction (negative pressure).
    # Standardize precision for final outputs
    precision = 3

    for _attempt in range(200):
        # Randomly select a fluid for the pipe
        pipe_fluid_name, rho1 = random.choice(list(PIPE_FLUIDS.items()))

        # Randomly select a fluid for the manometer
        manometer_fluid_name, rho2 = random.choice(list(MANOMETER_FLUIDS.items()))

        #  Ensure physical realism: manometer fluid must be denser
        # If the chosen pipe fluid is denser than the manometer fluid, re-select
        # the manometer fluid until it is denser.
        if rho2 <= rho1:
            continue

        # Randomize the vertical heights, in meters
        # h1: distance from pipe centerline down to the fluid interface
        h1 = round(random.uniform(0.1, 0.5), 2)
        # h2: the height difference between the two manometer fluid columns
        h2 = round(random.uniform(0.05, 0.75), 2)

        # 2. Perform the core calculations for verification. An integer
        # density times g = 9.81 times a 2-dp height is exact at 4 dp, so the
        # terms are displayed at 4 dp and nothing is rounded before Step 6
        # (D-037); at 2 dp, "1030 * 9.81 * 0.35 = 3536.51" sat on a half-way
        # tie. A gas density carries more decimals: a draw whose exact term
        # sits on a 4-dp tie is redrawn, as is one whose gauge pressure sits on
        # the 3-dp kPa tie of Step 6 (D-016 part 3).

        # Step A: Calculate the pressure contribution from the pipe fluid column (rho1*g*h1)
        pressure_term1_exact = rho1 * GRAVITY * h1

        # Step B: Calculate the pressure contribution from the manometer fluid column (rho2*g*h2)
        pressure_term2_exact = rho2 * GRAVITY * h2

        if _is_display_tie(pressure_term1_exact, 4) or _is_display_tie(pressure_term2_exact, 4):
            continue
        pressure_term1_pa = _as_printed(pressure_term1_exact, '.4f')
        pressure_term2_pa = _as_printed(pressure_term2_exact, '.4f')

        # Step C: Calculate the gauge pressure in Pascals (Pa) from the
        # displayed terms (D-016 part 2)
        gauge_pressure_pa = _as_printed(pressure_term2_pa - pressure_term1_pa, '.4f')

        # Condition: Pressure must be positive to match the problem description
        # (open arm level is "higher" implies P_pipe > P_atm)
        if gauge_pressure_pa <= 0:
            continue
        if _is_display_tie(gauge_pressure_pa / 1000, precision):
            continue
        break
    else:
        raise RuntimeError("utube_manometer: no closing sample in 200 draws")

    # Step D: Convert the final answer to kilopascals (kPa) for readability
    gauge_pressure_kpa = _as_printed(gauge_pressure_pa / 1000, f'.{precision}f')

    # 3. Generate the question and solution strings

    question = (
        f"A U-tube manometer using {manometer_fluid_name.lower()} (density = {rho2} kg/m^3) is "
        f"connected to a pipe carrying {pipe_fluid_name.lower()} (density = {rho1} kg/m^3). "
        f"The interface between the two fluids is {h1} m below the centerline of the pipe. "
        f"The level of {manometer_fluid_name.lower()} in the arm open to the atmosphere is {h2} m "
        f"higher than the interface. What is the gauge pressure in the pipe?"
    )

    solution = (
        f"**Given:**\n"
        f"Pipe Fluid: {pipe_fluid_name}, Density (rho1) = {rho1} kg/m^3\n"
        f"Manometer Fluid: {manometer_fluid_name}, Density (rho2) = {rho2} kg/m^3\n"
        f"Height from pipe centerline to interface (h1): {h1} m\n"
        f"Height difference in manometer fluid (h2): {h2} m\n"
        f"Acceleration due to Gravity (g): {GRAVITY} m/s^2\n\n"

        f"**Step 1:** State the Principle of Manometry\n"
        f"We can determine the pressure in the pipe by starting at the pipe's centerline, moving "
        f"through the fluid columns to the open end, and balancing the pressures. The pressure at the "
        f"same level within a continuous fluid at rest is equal.\n\n"

        f"**Step 2:** Formulate the Pressure Balance Equation\n"
        f"Let's establish a pressure equation starting from the pipe (P_pipe) and ending at the atmosphere (P_atm). "
        f"The reference level is the interface between the pipe fluid and the manometer fluid.\n"
        f"Pressure from pipe side at interface: P_pipe + (rho1 * g * h1)\n"
        f"Pressure from atmosphere side at interface: P_atm + (rho2 * g * h2)\n"
        f"Equating these gives: P_pipe + (rho1 * g * h1) = P_atm + (rho2 * g * h2)\n\n"

        f"**Step 3:** Solve for the Gauge Pressure\n"
        f"Gauge pressure is P_gauge = P_pipe - P_atm. Rearranging the equation:\n"
        f"P_gauge = (rho2 * g * h2) - (rho1 * g * h1)\n\n"

        f"**Step 4:** Calculate the Individual Pressure Terms\n"
        f"Pressure from manometer fluid column = {rho2} * {GRAVITY} * {h2} = {pressure_term2_pa:.4f} Pa\n"
        f"Pressure from pipe fluid column = {rho1} * {GRAVITY} * {h1} = {pressure_term1_pa:.4f} Pa\n\n"

        f"**Step 5:** Calculate the Final Gauge Pressure\n"
        f"P_gauge = {pressure_term2_pa:.4f} Pa - {pressure_term1_pa:.4f} Pa = {gauge_pressure_pa:.4f} Pa\n\n"

        f"**Step 6:** Convert the Answer to Kilopascals (kPa)\n"
        f"P_gauge = {gauge_pressure_pa:.4f} Pa / 1000 = {round(gauge_pressure_kpa, precision)} kPa\n\n"

        f"**Answer:**\n"
        f"The gauge pressure in the pipe is {round(gauge_pressure_kpa, precision)} kPa."
    )

    return question, solution


# Template 4 (Intermediate)
def template_floating_object_submersion_depth():
    """
    Fluid Statics: Floating Object Submersion Depth

    Scenario:
        This template applies the principle of buoyancy to determine how deep
        an object with a uniform cross-section (like a rectangular block or
        cylinder) will float in a liquid. The solution requires equating the
        object's total weight to the buoyant force acting on its submerged part.

    Core Equations:
        - Weight (W) = rho_object * g * V_total
        - Buoyant Force (F_B) = rho_fluid * g * V_submerged
        - At equilibrium (floating): W = F_B
        - This simplifies to: rho_object * V_total = rho_fluid * V_submerged

    Screen pass 1 (2026-09-23):
        One judge flagged that "stable, upright" floating is assumed and
        not checked, citing a wide, flat aerogel block (2.07 x 1.32 x
        1.09 m) that "would tip over"; another noted that random aspect
        ratios "can still be hydrostatically unstable". Checked with the
        metacentric criterion GM = BM - BG, BM = I/V_sub (see the sampling
        comment): the cited block is very stable, GM = +59.4 m, because a
        nearly weightless block has a huge metacentric radius, so that
        example is rejected; but 40.6% of first draws over 503 seeds,
        blocks and cylinders alike, typically tall sections at a mid-range
        density ratio, had GM <= 0 and would capsize. The stability
        condition is now a sampling constraint: shape and dimensions are
        redrawn until GM > 0 (smallest GM kept over 503 seeds: 1 mm).
        Question wording, physics and step structure are unchanged; the
        question and answer moved on 40.7% of 501 seeds (0-500), exactly
        the redrawn ones.

    Returns:
        tuple: A tuple containing:
            - str: A question asking for the submersion depth of a floating object.
            - str: A step-by-step solution to the problem.
    """
    # 1. Parameterize the inputs to ensure a valid floating scenario

    # Randomly select an object material and a fluid
    obj_material, rho_object = random.choice(list(MATERIAL_DENSITIES.items()))
    fluid_name, rho_fluid = random.choice(list(FLUID_DENSITIES.items()))

    # CRITICAL: Ensure the object will actually float 
    # Re-select until the object's density is less than the fluid's density.
    max_attempts = 20
    attempt = 0
    while rho_object >= rho_fluid and attempt < max_attempts:
        obj_material, rho_object = random.choice(list(MATERIAL_DENSITIES.items()))
        fluid_name, rho_fluid = random.choice(list(FLUID_DENSITIES.items()))
        attempt += 1
    
    # Fallback in the rare case a valid pair isn't found quickly
    if rho_object >= rho_fluid:
        # C3.3: the fallback names two table rows and reads their values FROM
        # the tables, so a corrected density reaches it. It hard-coded copies
        # (500, 998) that a correction would have left behind.
        obj_material = "Pine Wood"
        rho_object = MATERIAL_DENSITIES[obj_material]
        fluid_name = "Fresh Water"
        rho_fluid = FLUID_DENSITIES[fluid_name]

    # Randomly choose the object's shape (uniform cross-section) and its
    # dimensions, keeping only a draw that floats upright in stable
    # equilibrium, as the question asserts.
    #
    # Screen pass 1 (2026-09-23): the metacentric criterion for a floating
    # body, GM = BM - BG > 0 with BM = I / V_sub and, for a homogeneous body
    # of uniform section, BG = (H - h_sub) / 2. For the block the weaker axis
    # is the shorter horizontal side, I = L_long * W_short^3 / 12 over
    # V_sub = L * W * h_sub, so BM = W_short^2 / (12 h_sub); for the cylinder
    # I = pi r^4 / 4 over V_sub = pi r^2 h_sub, so BM = r^2 / (4 h_sub). A
    # draw with GM <= 0 would capsize rather than float upright; it is
    # redrawn (40.6% of first draws over 503 seeds).
    for _attempt in range(200):
        shape = random.choice(["rectangular block", "cylinder"])

        # Randomize dimensions based on shape
        if shape == "rectangular block":
            length = round(random.uniform(0.5, 3.0), 2)
            width = round(random.uniform(0.2, 2.0), 2)
            height = round(random.uniform(0.1, 1.5), 2) # This is the total vertical height
            shape_dims_str = f"dimensions {length} m (length) x {width} m (width) x {height} m (height)"
            shape_dims_given_str = (f"  - Length (L): {length} m\n"
                                    f"  - Width (W): {width} m\n"
                                    f"  - Total Height (H): {height} m")
            bm_times_h_sub = min(length, width) ** 2 / 12.0
        else:  # shape == "cylinder"
            radius = round(random.uniform(0.1, 1.5), 2)
            height = round(random.uniform(0.2, 2.5), 2) # This is the total vertical height
            shape_dims_str = f"a radius of {radius} m and a total height of {height} m"
            shape_dims_given_str = (f"  - Radius (r): {radius} m\n"
                                    f"  - Total Height (H): {height} m")
            bm_times_h_sub = radius ** 2 / 4.0

        h_sub = (rho_object / rho_fluid) * height
        metacentric_height = bm_times_h_sub / h_sub - (height - h_sub) / 2.0
        if metacentric_height > 0.0:
            break
    else:
        raise RuntimeError(
            "floating_object_submersion_depth: no upright-stable draw in 200 attempts")

    # Standardize precision for final outputs
    precision = 4

    # 2. Perform the core calculations for the solution

    # The core insight is that for a uniform cross-section, the area cancels out.
    # W = F_B  =>  rho_obj * g * (Area * H) = rho_fluid * g * (Area * h_sub)
    # Simplifying gives: rho_obj * H = rho_fluid * h_sub
    submersion_depth = (rho_object / rho_fluid) * height

    # 3. Generate the question and solution strings

    question = (
        f"A {shape} made of {obj_material.lower()} (density = {rho_object} kg/m^3) has "
        f"{shape_dims_str}. The object is placed in a large tank of {fluid_name.lower()} "
        f"(density = {rho_fluid} kg/m^3). Assuming the object floats in a stable, upright "
        f"position, calculate its submersion depth (the vertical height of the object "
        f"that is below the fluid surface)."
    )

    solution = (
        f"**Given:**\n"
        f"Object Shape: {shape.capitalize()}\n"
        f"{shape_dims_given_str}\n"
        f"Object Density (rho_obj): {rho_object} kg/m^3\n"
        f"Fluid Density (rho_fluid): {rho_fluid} kg/m^3\n\n"

        f"**Step 1:** State the Principle of Flotation\n"
        f"For an object to float, its total weight (W) must be equal to the buoyant force (F_B) "
        f"exerted by the fluid. The buoyant force is the weight of the displaced fluid.\n"
        f"Equilibrium Condition: W = F_B\n\n"

        f"**Step 2:** Express Weight and Buoyant Force using Densities\n"
        f"Weight (W) = rho_obj * g * V_total\n"
        f"Buoyant Force (F_B) = rho_fluid * g * V_submerged\n\n"
        f"Setting them equal: \n"
        f"rho_obj * g * V_total = rho_fluid * g * V_submerged\n\n"

        f"**Step 3:** Simplify for a Uniform Cross-Section\n"
        f"The acceleration of gravity (g) cancels from both sides. For an object with a uniform "
        f"cross-sectional area (A), the volumes can be expressed as V = A * height.\n"
        f"V_total = A * H_total\n"
        f"V_submerged = A * h_submerged\n\n"
        f"Substituting these into the equation:\n"
        f"rho_obj * (A * H_total) = rho_fluid * (A * h_submerged)\n\n"
        f"The cross-sectional area (A) also cancels out, leaving a simple ratio:\n"
        f"rho_obj * H_total = rho_fluid * h_submerged\n\n"

        f"**Step 4:** Solve for the Submersion Depth (h_submerged)\n"
        f"Rearranging the formula to solve for the unknown depth:\n"
        f"h_submerged = (rho_obj / rho_fluid) * H_total\n\n"
        f"Plugging in the given values:\n"
        f"h_submerged = ({rho_object} / {rho_fluid}) * {height}\n"
        f"h_submerged = {round(submersion_depth, precision)} m\n\n"

        f"**Answer:**\n"
        f"The submersion depth of the {shape} is {round(submersion_depth, precision)} m."
    )

    return question, solution


# Template 5 (Advanced)
def template_hydrostatic_force_on_plane():
    """
    Fluid Statics: Hydrostatic Force and Center of Pressure on a Submerged Plane

    Scenario:
        This template generates a comprehensive problem requiring the calculation of
        the magnitude of the resultant hydrostatic force on a submerged plane
        surface (e.g., a gate, window) and the location where this force acts,
        known as the center of pressure.

    Core Equations:
        - Resultant Force (F_R) = rho * g * h_c * A
        - Center of Pressure (y_R) = y_c + (I_xc / (y_c * A))
          where:
            - rho = fluid density
            - g = acceleration of gravity
            - h_c = vertical depth from free surface to the centroid of the area
            - A = area of the submerged surface
            - y_c = inclined distance from free surface to the centroid
            - I_xc = area moment of inertia about the centroidal axis parallel to the surface

    Returns:
        tuple: A tuple containing:
            - str: A question about hydrostatic force and center of pressure.
            - str: A detailed, step-by-step solution.
    """
    # 1. Parameterize the inputs with random values

    # Select fluid and shape
    fluid_name, rho = random.choice(list(FLUID_DENSITIES.items()))
    shape = random.choice(["rectangle", "circle"])

    # Set geometric and positional parameters
    angle_deg = random.choice([30, 45, 60, 90]) # Angle with the horizontal
    angle_rad = math.radians(angle_deg)
    h_top = round(random.uniform(0.5, 5.0), 2) # Vertical depth to top edge of the plane

    # Randomize dimensions based on shape and calculate geometric properties
    if shape == "rectangle":
        b = round(random.uniform(1.0, 4.0), 2)  # width
        h = round(random.uniform(1.0, 4.0), 2)  # height
        area = b * h
        I_xc = (b * h**3) / 12
        dist_to_centroid_along_plane = h / 2
        shape_desc = f"a rectangular gate with a width of {b} m and a height of {h} m"
        shape_given = f"  - Shape: Rectangle (width b = {b} m, height h = {h} m)"
    else:  # shape == "circle"
        r = round(random.uniform(0.5, 2.0), 2) # radius
        area = math.pi * r**2
        I_xc = (math.pi * r**4) / 4
        dist_to_centroid_along_plane = r
        shape_desc = f"a circular viewport with a radius of {r} m"
        shape_given = f"  - Shape: Circle (radius r = {r} m)"

    precision = 4

    # 2. Perform the core calculations for the solution

    # Step A: Determine the position of the centroid (yc and hc)
    # yc is the inclined distance from the free surface to the centroid
    # hc is the vertical depth from the free surface to the centroid
    if angle_deg == 90:
        y_top = h_top
        angle_desc = "is positioned vertically"
    else:
        y_top = h_top / math.sin(angle_rad)
        angle_desc = f"is inclined at an angle of {angle_deg} degrees to the horizontal"

    y_c = y_top + dist_to_centroid_along_plane
    h_c = y_c * math.sin(angle_rad)

    # Step B: Calculate the resultant force (F_R)
    force_resultant_N = rho * GRAVITY * h_c * area
    force_resultant_kN = force_resultant_N / 1000

    # Step C: Calculate the location of the center of pressure (y_R)
    # y_R is the inclined distance from the free surface to the center of pressure
    y_R = y_c + (I_xc / (y_c * area))

    # 3. Generate the question and solution strings

    question = (
        f"A submerged {shape_desc} {angle_desc} in a tank of {fluid_name.lower()} "
        f"(density = {rho} kg/m^3). The top edge of the gate is {h_top} m vertically below "
        f"the free surface. Calculate:\n"
        f"  a) The magnitude of the resultant hydrostatic force on the gate.\n"
        f"  b) The location of the center of pressure, measured along the incline of the gate "
        f"from the free surface."
    )

    solution = (
        f"**Given:**\n"
        f"{shape_given}\n"
        f"Fluid: {fluid_name} (rho = {rho} kg/m^3)\n"
        f"Vertical depth to top edge (h_top): {h_top} m\n"
        f"Angle of inclination (theta): {angle_deg} degrees\n\n"

        f"**Step 1:** Calculate Geometric Properties of the Gate\n"
        f"Area (A): {area:.{precision}f} m^2\n"
        f"Moment of Inertia about centroid (I_xc): {I_xc:.{precision}f} m^4\n\n"

        f"**Step 2:** Determine the Location of the Centroid (yc and hc)\n"
        f"The centroid is the geometric center of the gate. We need its position relative to the free surface.\n"
        f"  - 'y' distances are measured along the plane's incline.\n"
        f"  - 'h' distances are measured vertically.\n\n"
        f"Inclined distance from surface to top edge (y_top) = h_top / sin(theta)\n"
        f"  y_top = {h_top} / sin({angle_deg}°) = {y_top:.{precision}f} m\n"
        f"Distance from top edge to centroid along plane = {dist_to_centroid_along_plane:.{precision}f} m\n"
        f"Inclined distance from surface to centroid (y_c) = y_top + dist_to_centroid\n"
        f"  y_c = {y_top:.{precision}f} + {dist_to_centroid_along_plane:.{precision}f} = {y_c:.{precision}f} m\n\n"
        f"Vertical depth to centroid (h_c) = y_c * sin(theta)\n"
        f"  h_c = {y_c:.{precision}f} * sin({angle_deg}°) = {h_c:.{precision}f} m\n\n"

        f"**Step 3:** Calculate the Resultant Hydrostatic Force (F_R)\n"
        f"The force is the pressure at the centroid multiplied by the total area.\n"
        f"F_R = rho * g * h_c * A\n"
        f"F_R = {rho} * {GRAVITY} * {h_c:.{precision}f} * {area:.{precision}f}\n"
        f"F_R = {force_resultant_N:.2f} N\n"
        f"F_R = {force_resultant_kN:.{precision}f} kN\n\n"

        f"**Step 4:** Calculate the Center of Pressure (y_R)\n"
        f"The center of pressure is the point where the resultant force acts. It is always below the centroid.\n"
        f"y_R = y_c + (I_xc / (y_c * A))\n"
        f"y_R = {y_c:.{precision}f} + ({I_xc:.{precision}f} / ({y_c:.{precision}f} * {area:.{precision}f}))\n"
        f"y_R = {y_R:.{precision}f} m\n\n"

        f"**Answer:**\n"
        f"a) The magnitude of the resultant hydrostatic force is **{force_resultant_kN:.{precision}f} kN**.\n"
        f"b) The center of pressure is located **{y_R:.{precision}f} m** from the free surface, measured down along the angle of the gate."
    )

    return question, solution


def main():
    """
    Generate numerous instances of each fluid statics template 
    with different random seeds and write the results to a JSONL file.
    """
    import json
    import os

    # Define the output path (Modify this path according to where you are running the code from)
    output_file = "testset/mechanical_engineering/fluid_mechanics/fluid_statics.jsonl"

    # Create the directory if it doesn't exist
    os.makedirs(os.path.dirname(output_file), exist_ok=True)

    # List of template functions with their ID and level
    templates = [
        (template_hydrostatic_pressure_at_depth, "hydrostatic_pressure_at_depth", "Easy"),
        (template_basic_buoyant_force, "basic_buoyant_force", "Easy"),
        (template_utube_manometer, "utube_manometer", "Intermediate"),
        (template_floating_object_submersion_depth, "floating_object_submersion_depth", "Intermediate"),
        (template_hydrostatic_force_on_plane, "hydrostatic_force_on_plane", "Advanced"),
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
                "domain": "fluid_mechanics",
                "area": "fluid_statics",
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
