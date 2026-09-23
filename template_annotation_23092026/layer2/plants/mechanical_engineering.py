"""Planted defects for the mechanical engineering branch (Layer 2). See CONTRACT.md.

plant_mec_1  constant    template_shear_stress_torsion: J computed with the area-moment factor (pi/4) c^4 instead of (pi/2) c^4, so tau_max is doubled.
plant_mec_2  unit        template_hydrostatic_pressure_at_depth: kPa <-> Pa converted by 100 instead of 1000, in Step 1 and Step 4 alike.
plant_mec_3  sign        template_utube_manometer: the pipe-fluid column is subtracted on the pipe side, so P_gauge = rho2*g*h2 + rho1*g*h1.
plant_mec_4  arithmetic  template_undamped_natural_frequency_translational: Step 1 prints omega_n 4% above sqrt(k/m); the rest of the trace keeps the true value.
"""

PLANTS = [
    {
        'plant_id': 'plant_mec_1',
        'base': 'template_shear_stress_torsion',
        'defect_class': 'constant',
        'description': (
            "Polar moment of inertia computed with the area-moment factor: "
            "J = (pi/4) c^4 for the solid shaft and (pi/4)(c_o^4 - c_i^4) for the hollow "
            "one, instead of (pi/2) c^4 = pi d^4/32. Used consistently in the docstring, "
            "the printed formula lines and the computation, so every printed J is half "
            "the true value and tau_max = T c / J is twice the correct stress."
        ),
        'detectable_by': (
            "Step 2 prints J = (pi / 4) * c^4 (pi/4 is the area moment I of a circle); the "
            "polar moment is (pi / 2) * c^4 = pi d^4 / 32. Recomputing J from the printed "
            "radius gives twice the printed value, and the tau_max of Step 3 and the MPa "
            "answer are twice the correct stress."
        ),
        'edits': [
            ("J_solid = (pi / 2) * c^4",
             "J_solid = (pi / 4) * c^4"),
            ("J_hollow = (pi / 2) * (c_outer^4 - c_inner^4)",
             "J_hollow = (pi / 4) * (c_outer^4 - c_inner^4)"),
            ("polar_moment_J = _as_printed((math.pi / 2) * (c_outer ** 4), '.3e')",
             "polar_moment_J = _as_printed((math.pi / 4) * (c_outer ** 4), '.3e')"),
            ('j_calculation_str = f"J = (pi / 2) * c^4 = (pi / 2) * ({c_outer})^4',
             'j_calculation_str = f"J = (pi / 4) * c^4 = (pi / 4) * ({c_outer})^4'),
            ("polar_moment_J = _as_printed((math.pi / 2) * (c_outer ** 4 - c_inner ** 4), '.3e')",
             "polar_moment_J = _as_printed((math.pi / 4) * (c_outer ** 4 - c_inner ** 4), '.3e')"),
            (r'f"J = (pi / 2) * (c_outer^4 - c_inner^4)\n"',
             r'f"J = (pi / 4) * (c_outer^4 - c_inner^4)\n"'),
            ('f"   J = (pi / 2) * (({c_outer})^4 - ({c_inner})^4)',
             'f"   J = (pi / 4) * (({c_outer})^4 - ({c_inner})^4)'),
        ],
        'reskin': [
            ('f"A {geometry_desc} is subjected to a torque of {torque} N.m. "',
             'f"A {geometry_desc} in a pump drive line is subjected to a torque of {torque} N.m. "'),
            ("d_outer = random.randint(30, 150)  # Outer diameter in mm",
             "d_outer = random.randint(40, 140)  # Outer diameter in mm"),
        ],
    },
    {
        'plant_id': 'plant_mec_2',
        'base': 'template_hydrostatic_pressure_at_depth',
        'defect_class': 'unit',
        'description': (
            "Kilopascal-to-pascal conversion done with a factor of 100 instead of 1000, "
            "consistently in both directions: Step 1 converts the surface pressure as "
            "kPa * 100 = Pa and Step 4 converts the pascal total back as Pa / 100 = kPa "
            "(computation and printed lines alike). The hydrostatic contribution to the "
            "quoted kPa answer is therefore ten times too large."
        ),
        'detectable_by': (
            "Step 1 prints 'P_surface_abs = 101.325 kPa * 100 = 10132.5 Pa' (1 kPa is "
            "1000 Pa) and Step 4 divides the pascal total by 100 to return to kPa, so the "
            "answer (e.g. about 1000 kPa at 10 m of water) is ten times what rho*g*h "
            "implies."
        ),
        'edits': [
            ("surface_pressure_pa = _as_printed(surface_pressure_abs_kpa * 1000, '.1f')",
             "surface_pressure_pa = _as_printed(surface_pressure_abs_kpa * 100, '.1f')"),
            ("absolute_pressure_kpa = _as_printed(absolute_pressure_pa / 1000, f'.{precision}f')",
             "absolute_pressure_kpa = _as_printed(absolute_pressure_pa / 100, f'.{precision}f')"),
            (r'f"Convert to Pascals: P_surface_abs = {surface_pressure_abs_kpa:.3f} kPa * 1000 = {surface_pressure_pa:.1f} Pa\n\n"',
             r'f"Convert to Pascals: P_surface_abs = {surface_pressure_abs_kpa:.3f} kPa * 100 = {surface_pressure_pa:.1f} Pa\n\n"'),
            (r'f"P_absolute = {absolute_pressure_pa:.1f} Pa / 1000 = {round(absolute_pressure_kpa, precision)} kPa\n\n"',
             r'f"P_absolute = {absolute_pressure_pa:.1f} Pa / 100 = {round(absolute_pressure_kpa, precision)} kPa\n\n"'),
        ],
        'reskin': [
            ('f"An object is located {depth_h} m below the surface of a tank containing "',
             'f"A pressure transducer is mounted {depth_h} m below the free surface of a storage tank containing "'),
            ("depth_h = round(random.uniform(2.0, 25.0), 2)",
             "depth_h = round(random.uniform(3.0, 20.0), 2)"),
        ],
    },
    {
        'plant_id': 'plant_mec_3',
        'base': 'template_utube_manometer',
        'defect_class': 'sign',
        'description': (
            "Sign of the pipe-fluid column flipped consistently: the interface pressure on "
            "the pipe side is written P_pipe - rho1*g*h1 (Step 2), so the gauge pressure "
            "becomes P_gauge = rho2*g*h2 + rho1*g*h1 in the docstring, Step 3, Step 5 and "
            "the computation. The kPa answer is too high by 2*rho1*g*h1."
        ),
        'detectable_by': (
            "Step 2 writes the interface pressure on the pipe side as P_pipe - (rho1 * g * h1) "
            "although the interface is h1 BELOW the centerline, so the pipe-fluid column must "
            "add to P_pipe; Step 3 and Step 5 then add the two column terms (X Pa + Y Pa) "
            "where the balance requires their difference."
        ),
        'edits': [
            ("P_gauge = (rho_2 * g * h2) - (rho_1 * g * h1)",
             "P_gauge = (rho_2 * g * h2) + (rho_1 * g * h1)"),
            ("gauge_pressure_pa = _as_printed(pressure_term2_pa - pressure_term1_pa, '.4f')",
             "gauge_pressure_pa = _as_printed(pressure_term2_pa + pressure_term1_pa, '.4f')"),
            (r'f"Pressure from pipe side at interface: P_pipe + (rho1 * g * h1)\n"',
             r'f"Pressure from pipe side at interface: P_pipe - (rho1 * g * h1)\n"'),
            (r'f"Equating these gives: P_pipe + (rho1 * g * h1) = P_atm + (rho2 * g * h2)\n\n"',
             r'f"Equating these gives: P_pipe - (rho1 * g * h1) = P_atm + (rho2 * g * h2)\n\n"'),
            (r'f"P_gauge = (rho2 * g * h2) - (rho1 * g * h1)\n\n"',
             r'f"P_gauge = (rho2 * g * h2) + (rho1 * g * h1)\n\n"'),
            (r'f"P_gauge = {pressure_term2_pa:.4f} Pa - {pressure_term1_pa:.4f} Pa = {gauge_pressure_pa:.4f} Pa\n\n"',
             r'f"P_gauge = {pressure_term2_pa:.4f} Pa + {pressure_term1_pa:.4f} Pa = {gauge_pressure_pa:.4f} Pa\n\n"'),
        ],
        'reskin': [
            ('f"connected to a pipe carrying {pipe_fluid_name.lower()} (density = {rho1} kg/m^3). "',
             'f"connected to a horizontal process pipe carrying {pipe_fluid_name.lower()} (density = {rho1} kg/m^3). "'),
            ("h1 = round(random.uniform(0.1, 0.5), 2)",
             "h1 = round(random.uniform(0.12, 0.45), 2)"),
        ],
    },
    {
        'plant_id': 'plant_mec_4',
        'base': 'template_undamped_natural_frequency_translational',
        'defect_class': 'arithmetic',
        'description': (
            "Step 1's printed result altered on that single line: omega_n is printed as "
            "1.04 times the sqrt(k/m) the line's operands give, while f_n, tau_n and the "
            "Answer are carried on the true omega_n."
        ),
        'detectable_by': (
            "Step 1 prints 'omega_n = sqrt(k / m)' with the given k and m and then a result "
            "4% above that square root; Step 2 divides a different omega_n (the true one) "
            "by 2*pi and the Answer repeats it, so the Step 1 line follows neither from its "
            "operands nor from the rest of the trace."
        ),
        'edits': [
            (r'f"omega_n = {round(omega_n, precision)} rad/s\n\n"',
             r'f"omega_n = {round(omega_n * 1.04, precision)} rad/s\n\n"'),
        ],
        'reskin': [
            ('f"An undamped single-degree-of-freedom system consists of a mass of {mass} kg "',
             'f"An instrument package, modelled as an undamped single-degree-of-freedom system, consists of a mass of {mass} kg "'),
            ("mass = round(random.uniform(0.5, 750.0), 2)  # in kg",
             "mass = round(random.uniform(2.0, 400.0), 2)  # in kg"),
        ],
    },
]
