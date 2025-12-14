import random
import math
from data.templates.branches.mechanical_engineering.constants import MATERIAL_PROPERTIES


# Template 1 (Easy)
def template_basic_stress_strain():
    """
    Basic Normal Stress and Strain

    Scenario:
        This template tests the fundamental definitions of normal stress (sigma) and
        normal strain (epsilon). Given a simple prismatic bar with a known geometry
        and an applied axial load, the user must calculate these two basic quantities.

    Core Equations:
        Normal Stress: sigma = P / A
        Normal Strain: epsilon = delta / L

    Returns:
        tuple: A tuple containing:
            - str: A question asking for the normal stress and strain.
            - str: A step-by-step solution showing the calculations.
    """
    # 1. Parameterize inputs with random values, choosing a unit system first.
    use_si_units = random.choice([True, False])
    shape = random.choice(['circular', 'square'])
    precision = 3  # Standardize precision for numerical stability and formatting

    if use_si_units:
        #  SI Unit System 
        load = random.randint(10, 500)  # in kN
        length = round(random.uniform(0.5, 4.0), 2)  # in meters
        elongation = round(random.uniform(0.5, 5.0), 2)  # in mm
        
        if shape == 'circular':
            dimension_val = random.randint(20, 100)  # diameter in mm
            dimension_name = "diameter"
            area = math.pi * (dimension_val / 2)**2  # mm²
            dim_str = f"{dimension_val} mm"
        else:  # square
            dimension_val = random.randint(20, 100)  # side length in mm
            dimension_name = "side length"
            area = dimension_val**2  # mm²
            dim_str = f"{dimension_val} mm"
            
        load_str = f"{load} kN"
        length_str = f"{length} m"
        elongation_str = f"{elongation} mm"
        
        # Core Calculations (using units for clarity: N, mm, MPa)
        # 1 MPa = 1 N/mm^2
        stress = (load * 1000) / area  # Stress in MPa
        # Strain (unitless, but convert units to be consistent)
        # Elongation in mm, Length in m -> convert length to mm
        strain = elongation / (length * 1000)
        
        stress_unit = "MPa"
        strain_unit_explanation = "mm/m or unitless"

    else:
        #  US Customary Unit System 
        load = random.randint(5, 100)  # in kips
        length = round(random.uniform(24.0, 120.0), 1)  # in inches
        elongation = round(random.uniform(0.05, 0.25), 3)  # in inches
        
        if shape == 'circular':
            dimension_val = round(random.uniform(1.0, 5.0), 2)  # diameter in inches
            dimension_name = "diameter"
            area = math.pi * (dimension_val / 2)**2  # in²
            dim_str = f"{dimension_val} in"
        else:  # square
            dimension_val = round(random.uniform(1.0, 5.0), 2)  # side length in inches
            dimension_name = "side length"
            area = dimension_val**2
            dim_str = f"{dimension_val} in"
            
        load_str = f"{load} kips"
        length_str = f"{length} in"
        elongation_str = f"{elongation} in"
        
        # Core Calculations (using units for clarity: kips, in, ksi)
        # 1 ksi = 1 kip/in^2
        stress = load / area  # Stress in ksi
        # Strain (unitless, since both length and elongation are in inches)
        strain = elongation / length
        
        stress_unit = "ksi"
        strain_unit_explanation = "in/in or unitless"

    # 2. Generate the question and solution strings
    question = (
        f"A prismatic bar with a length of {length_str} has a {shape} cross-section "
        f"with a {dimension_name} of {dim_str}. The bar is subjected to an axial "
        f"tensile load of {load_str}, causing it to elongate by {elongation_str}.\n\n"
        f"Determine the following:\n"
        f"a) The normal stress in the bar.\n"
        f"b) The normal strain in the bar."
    )

    solution = (
        f"**Given:**\n"
        f"Load (P): {load_str}\n"
        f"Length (L): {length_str}\n"
        f"Elongation (delta): {elongation_str}\n"
        f"Cross-Section: {shape.capitalize()} with {dimension_name} = {dim_str}\n\n"
        
        f"**Step 1:** Calculate the Cross-Sectional Area (A)\n"
        f"The area of a {shape} cross-section is calculated as follows:\n"
    )
    
    if shape == 'circular':
        solution += (
            f"A = pi * (d/2)^2 = pi * ({dimension_val}/2)^2 = {round(area, precision+1)} "
            f"{'mm^2' if use_si_units else 'in^2'}\n\n"
        )
    else: # square
        solution += (
            f"A = side^2 = ({dimension_val})^2 = {round(area, precision+1)} "
            f"{'mm^2' if use_si_units else 'in^2'}\n\n"
        )
        
    solution += (
        f"**Step 2:** Calculate the Normal Stress (sigma)\n"
        f"Normal stress is defined as the force per unit area: sigma = P / A.\n"
    )
    
    if use_si_units:
        solution += (
            f"To get the result in Megapascals (MPa), we use the load in Newtons (N) and the area in mm^2 (1 MPa = 1 N/mm^2).\n"
            f"P = {load} kN = {load * 1000} N\n"
            f"A = {round(area, precision+1)} mm^2\n"
            f"sigma = ({load * 1000} N) / ({round(area, precision+1)} mm^2) = {round(stress, precision)} MPa\n\n"
        )
    else: # US units
        solution += (
            f"To get the result in kips per square inch (ksi), we use the load in kips and the area in in^2 (1 ksi = 1 kip/in^2).\n"
            f"P = {load} kips\n"
            f"A = {round(area, precision+1)} in^2\n"
            f"sigma = ({load} kips) / ({round(area, precision+1)} in^2) = {round(stress, precision)} ksi\n\n"
        )

    solution += (
        f"**Step 3:** Calculate the Normal Strain (epsilon)\n"
        f"Normal strain is the change in length per unit of original length: epsilon = delta / L.\n"
        f"It is a dimensionless quantity, but units must be consistent for calculation.\n"
    )

    if use_si_units:
        solution += (
            f"We convert the length from meters to millimeters to match the elongation units.\n"
            f"delta = {elongation} mm\n"
            f"L = {length} m = {length * 1000} mm\n"
            f"epsilon = ({elongation} mm) / ({length * 1000} mm) = {strain:.3e}\n\n"
        )
    else: # US units
        solution += (
            f"Both elongation and length are already in inches, so we can divide directly.\n"
            f"delta = {elongation} in\n"
            f"L = {length} in\n"
            f"epsilon = ({elongation} in) / ({length} in) = {strain:.3e}\n\n"
        )

    solution += (
        f"**Answer:**\n"
        f"a) Normal Stress (sigma) = **{round(stress, precision)} {stress_unit}**\n"
        f"b) Normal Strain (epsilon) = **{strain:.3e}** ({strain_unit_explanation})"
    )

    return question, solution


# Template 2 (Easy)
def template_axial_deformation():
    """
    Axial Deformation (Hooke's Law)

    Scenario:
        This template combines stress, strain, and the material's Modulus of
        Elasticity (E). It tests the ability to calculate the total elongation
        or compression of a rod under an axial load. It can also be inverted
        to find the maximum load for a given allowable deformation.

    Core Equation:
        Deformation: delta = (P * L) / (A * E)

    Returns:
        tuple: A tuple containing:
            - str: A question about axial deformation or allowable load.
            - str: A step-by-step solution showing the calculations.
    """
    # 1. Parameterize inputs with random values
    use_si_units = random.choice([True, False])
    shape = random.choice(['circular', 'square'])
    material_name = random.choice(list(MATERIAL_PROPERTIES.keys()))
    solve_for_load = random.choice([True, False])  # The variation
    precision = 3

    if use_si_units:
        #  SI Unit System 
        material_E = MATERIAL_PROPERTIES[material_name]['E_GPa'] # in GPa
        length = round(random.uniform(0.5, 5.0), 2) # in m
        
        if shape == 'circular':
            dimension_val = random.randint(25, 120) # diameter in mm
            dimension_name = "diameter"
            area = math.pi * (dimension_val / 2)**2
            dim_str = f"{dimension_val} mm"
        else: # square
            dimension_val = random.randint(25, 120) # side in mm
            dimension_name = "side length"
            area = dimension_val**2
            dim_str = f"{dimension_val} mm"
        
        E_str = f"{material_E} GPa"
        length_str = f"{length} m"
        
        if solve_for_load:
            # We are solving for P given delta_max
            max_elongation = round(random.uniform(1.0, 6.0), 2) # in mm
            max_elongation_str = f"{max_elongation} mm"
            # Core Calculation: P = (delta * A * E) / L
            # Units: mm * mm^2 * (N/mm^2) / mm = N
            load_N = (max_elongation * area * (material_E * 1000)) / (length * 1000)
            load = load_N / 1000 # convert to kN
            load_str = f"{round(load, precision)} kN"
        else:
            # We are solving for delta given P
            load = random.randint(50, 800) # in kN
            load_str = f"{load} kN"
            # Core Calculation: delta = (P * L) / (A * E)
            # Units: (N * mm) / (mm^2 * N/mm^2) = mm
            elongation = ((load * 1000) * (length * 1000)) / (area * (material_E * 1000))
            elongation_str = f"{round(elongation, precision)} mm"

    else:
        #  US Customary Unit System 
        material_E = MATERIAL_PROPERTIES[material_name]['E_ksi'] # in ksi
        length = round(random.uniform(12.0, 150.0), 1) # in inches
        
        if shape == 'circular':
            dimension_val = round(random.uniform(0.5, 6.0), 2) # diameter in inches
            dimension_name = "diameter"
            area = math.pi * (dimension_val / 2)**2
            dim_str = f"{dimension_val} in"
        else: # square
            dimension_val = round(random.uniform(0.5, 6.0), 2) # side in inches
            dimension_name = "side length"
            area = dimension_val**2
            dim_str = f"{dimension_val} in"
            
        E_str = f"{material_E} ksi"
        length_str = f"{length} in"

        if solve_for_load:
            # Solving for P given delta_max
            max_elongation = round(random.uniform(0.02, 0.2), 4) # in inches
            max_elongation_str = f"{max_elongation} in"
            # Core Calculation: P = (delta * A * E) / L
            # Units: in * in^2 * (kips/in^2) / in = kips
            load = (max_elongation * area * material_E) / length
            load_str = f"{round(load, precision)} kips"
        else:
            # Solving for delta given P
            load = random.randint(10, 200) # in kips
            load_str = f"{load} kips"
            # Core Calculation: delta = (P * L) / (A * E)
            # Units: (kips * in) / (in^2 * kips/in^2) = in
            elongation = (load * length) / (area * material_E)
            elongation_str = f"{round(elongation, precision)} in"
    
    # 2. Generate the question and solution strings
    if solve_for_load:
        question = (
            f"A rod made of {material_name} has a length of {length_str} and a {shape} "
            f"cross-section with a {dimension_name} of {dim_str}. The modulus of "
            f"elasticity for {material_name} is {E_str}.\n\n"
            f"What is the maximum allowable axial load the rod can support if the "
            f"total elongation is not to exceed {max_elongation_str}?"
        )
        
        solution = (
            f"**Given:**\n"
            f"Material: {material_name} (E = {E_str})\n"
            f"Length (L): {length_str}\n"
            f"Cross-Section: {shape.capitalize()}, {dimension_name} = {dim_str}\n"
            f"Maximum Elongation (delta_max): {max_elongation_str}\n\n"
            
            f"**Step 1:** Calculate the Cross-Sectional Area (A)\n"
        )
        if shape == 'circular':
            solution += (f"  A = pi * (d/2)^2 = pi * ({dimension_val}/2)^2 = {round(area, precision+1)} {'mm^2' if use_si_units else 'in^2'}\n\n")
        else: # square
            solution += (f"  A = side^2 = ({dimension_val})^2 = {round(area, precision+1)} {'mm^2' if use_si_units else 'in^2'}\n\n")

        solution += (
            f"**Step 2:** Calculate the Allowable Load (P)\n"
            f"The formula for axial deformation is delta = (P * L) / (A * E).\n"
            f"We can rearrange this to solve for the load P: P = (delta * A * E) / L.\n\n"
            f"**Unit Conversion and Calculation:**\n"
        )
        if use_si_units:
            solution += (
                f"For consistency, we will use units of Newtons (N) and millimeters (mm).\n"
                f"delta = {max_elongation} mm\n"
                f"A = {round(area, precision+1)} mm^2\n"
                f"E = {material_E} GPa = {material_E * 1000} MPa = {material_E * 1000} N/mm^2\n"
                f"L = {length} m = {length * 1000} mm\n\n"
                f"P = ({max_elongation} * {round(area, precision+1)} * {material_E * 1000}) / {length * 1000}\n"
                f"P = {round(load_N, precision)} N = {round(load, precision)} kN\n\n"
            )
        else: # US units
            solution += (
                f"The units are already consistent (kips and inches).\n"
                f"delta = {max_elongation} in\n"
                f"A = {round(area, precision+1)} in^2\n"
                f"E = {material_E} ksi\n"
                f"L = {length} in\n\n"
                f"P = ({max_elongation} * {round(area, precision+1)} * {material_E}) / {length}\n"
                f"P = {round(load, precision)} kips\n\n"
            )
        solution += f"**Answer:**\n  The maximum allowable load is **{load_str}**."
    else: # solve_for_deformation
        question = (
            f"A {material_name} rod with a length of {length_str} is subjected to an axial "
            f"tensile load of {load_str}. The rod has a {shape} cross-section with a "
            f"{dimension_name} of {dim_str}. The modulus of elasticity for "
            f"{material_name} is {E_str}.\n\n"
            f"Calculate the total elongation of the rod."
        )
        solution = (
            f"**Given:**\n"
            f"Material: {material_name} (E = {E_str})\n"
            f"Load (P): {load_str}\n"
            f"Length (L): {length_str}\n"
            f"Cross-Section: {shape.capitalize()}, {dimension_name} = {dim_str}\n\n"

            f"**Step 1:** Calculate the Cross-Sectional Area (A)\n"
        )
        if shape == 'circular':
            solution += (f"  A = pi * (d/2)^2 = pi * ({dimension_val}/2)^2 = {round(area, precision+1)} {'mm^2' if use_si_units else 'in^2'}\n\n")
        else: # square
            solution += (f"  A = side^2 = ({dimension_val})^2 = {round(area, precision+1)} {'mm^2' if use_si_units else 'in^2'}\n\n")
        
        solution += (
            f"**Step 2:** Calculate the Elongation (delta)\n"
            f"The formula for axial deformation is delta = (P * L) / (A * E).\n\n"
            f"**Unit Conversion and Calculation:**\n"
        )
        if use_si_units:
            solution += (
                f"For consistency, we will use units of Newtons (N) and millimeters (mm).\n"
                f"P = {load} kN = {load * 1000} N\n"
                f"L = {length} m = {length * 1000} mm\n"
                f"A = {round(area, precision+1)} mm^2\n"
                f"E = {material_E} GPa = {material_E * 1000} MPa = {material_E * 1000} N/mm^2\n\n"
                f"delta = (({load * 1000}) * ({length * 1000})) / (({round(area, precision+1)}) * ({material_E * 1000}))\n"
                f"delta = {round(elongation, precision)} mm\n\n"
            )
        else: # US units
            solution += (
                f"The units are already consistent (kips and inches).\n"
                f"P = {load} kips\n"
                f"L = {length} in\n"
                f"A = {round(area, precision+1)} in^2\n"
                f"E = {material_E} ksi\n\n"
                f"delta = ({load} * {length}) / ({round(area, precision+1)} * {material_E})\n"
                f"delta = {round(elongation, precision)} in\n\n"
            )
        solution += f"**Answer:**\n  The total elongation of the rod is **{elongation_str}**."

    return question, solution


# Template 3 (Intermediate)
def template_multi_segment_rod():
    """
    Deformation of a Multi-Segment Rod

    Scenario:
        This template involves a composite rod made of segments joined end-to-end and
        fixed at one end. The segments can have different materials, lengths, or
        cross-sectional areas. External loads are applied at the junctions. The user
        must calculate the total deformation by summing the deformations of each segment,
        which requires first finding the internal force in each section.

    Core Equation:
        Total Deformation: delta_total = sum(delta_i) = sum( (P_i * L_i) / (A_i * E_i) )

    Returns:
        tuple: A tuple containing:
            - str: A question about the total deformation of a composite rod.
            - str: A step-by-step solution.
    """
    # 1. Parameterize inputs with validation loop
    precision = 4
    
    # Define a safe list of structural materials (exclude rubber for this high-load problem)
    structural_materials = [k for k in MATERIAL_PROPERTIES.keys() if 'Rubber' not in k]

    while True:
        use_si_units = random.choice([True, False])
        num_segments = random.choice([2, 3])
        
        segments = []
        loads = {} # Loads applied at nodes {node_index: load_value}
        
        # Generate properties for each segment
        for i in range(num_segments):
            material_name = random.choice(structural_materials)
            segment = {'material_name': material_name}
            
            if use_si_units:
                segment['length'] = round(random.uniform(0.2, 1.5), 2) # meters
                diameter = round(random.uniform(20, 75), 1) # mm
                segment['area'] = math.pi * (diameter / 2)**2
                segment['E'] = MATERIAL_PROPERTIES[material_name]['E_GPa']
                segment['dim_str'] = f"{diameter} mm diameter"
                loads[i+1] = random.randint(-250, 250) # kN
            else:
                segment['length'] = round(random.uniform(10.0, 60.0), 1) # inches
                diameter = round(random.uniform(0.75, 3.0), 1) # inches
                segment['area'] = math.pi * (diameter / 2)**2
                segment['E'] = MATERIAL_PROPERTIES[material_name]['E_ksi']
                segment['dim_str'] = f"{diameter} in diameter"
                loads[i+1] = random.randint(-60, 60) # kips

            segments.append(segment)

        # 2. Perform Core Calculations
        total_deformation = 0
        cumulative_load = 0
        max_strain = 0
        
        # Calculate internal forces and deformations
        # Iterate backwards from the last segment to the first
        for i in range(num_segments - 1, -1, -1):
            node_index = i + 1
            cumulative_load += loads[node_index]
            segments[i]['internal_load'] = cumulative_load

        for i in range(num_segments):
            P = segments[i]['internal_load']
            L = segments[i]['length']
            A = segments[i]['area']
            E = segments[i]['E']
            
            if use_si_units:
                # delta = (P_N * L_mm) / (A_mm2 * E_MPa)
                # P in kN -> *1000 -> N
                # L in m -> *1000 -> mm
                # E in GPa -> *1000 -> MPa
                delta = (P * 1000 * L * 1000) / (A * E * 1000) # mm
                current_strain = abs(delta / (L * 1000))
            else:
                delta = (P * L) / (A * E) # inches
                current_strain = abs(delta / L)
            
            segments[i]['deformation'] = delta
            total_deformation += delta
            if current_strain > max_strain:
                max_strain = current_strain

        # Validate: Ensure max strain is < 1% (0.01) for linear elasticity to hold
        if max_strain < 0.01:
            break # Valid problem generated

    # 3. Generate Question and Solution Strings
    question = (
        f"A composite rod, fixed at one end, consists of {num_segments} segments "
        f"joined end-to-end as described below:\n"
    )
    for i, seg in enumerate(segments):
        unit_L = "m" if use_si_units else "in"
        question += (
            f"  - Segment {i+1}: Made of {seg['material_name']}, has a length of {seg['length']} {unit_L}, "
            f"and a cross-section that is {seg['dim_str']}.\n"
        )
    
    question += f"\nThe rod is subjected to the following external axial loads (positive = tension, negative = compression):\n"
    for i in range(num_segments):
        node_index = i + 1
        load_val = loads[node_index]
        unit_P = "kN" if use_si_units else "kips"
        
        if node_index < num_segments:
            location_desc = f"at the junction between Segment {node_index} and {node_index+1}"
        else:
            location_desc = "at the free end"
            
        question += f"  - A load of {load_val} {unit_P} is applied {location_desc}.\n"

    question += (
        f"\nGiven the Moduli of Elasticity:\n"
    )
    unique_materials = {s['material_name'] for s in segments}
    for mat in sorted(list(unique_materials)):
        E_val = MATERIAL_PROPERTIES[mat]['E_GPa' if use_si_units else 'E_ksi']
        unit_E = "GPa" if use_si_units else "ksi"
        question += f"  - E_{mat} = {E_val} {unit_E}\n"

    question += f"\nDetermine the total deformation of the composite rod."

    # Solution
    solution = f"**Given:** The properties and loads as described in the question.\n\n"
    solution += (
        f"**Step 1:** Calculate the Internal Force (P) in Each Segment\n"
        f"To find the internal force in any segment, we make a virtual 'cut' and sum all external forces acting on one side of the cut. We will work from the free end (right) to the fixed wall (left). A positive force indicates tension, and a negative force indicates compression.\n"
    )
    unit_P = "kN" if use_si_units else "kips"
    for i in range(num_segments - 1, -1, -1):
        seg_num = i + 1
        internal_load = segments[i]['internal_load']
        load_desc = " (Tension)" if internal_load > 0 else " (Compression)" if internal_load < 0 else " (No force)"
        
        if seg_num == num_segments: # Last segment
             solution += f"  - **Segment {seg_num}:** P{seg_num} = {loads[seg_num]} {unit_P}{load_desc}\n"
        else:
            prev_load_desc = " (T)" if segments[i+1]['internal_load'] > 0 else " (C)" if segments[i+1]['internal_load'] < 0 else " (0)"
            solution += f"  - **Segment {seg_num}:** P{seg_num} = P{seg_num+1} + {loads[seg_num]} = {segments[i+1]['internal_load']}{prev_load_desc} + {loads[seg_num]} = {internal_load} {unit_P}{load_desc}\n"

    solution += f"\nSummary of Internal Forces:\n"
    for i, seg in enumerate(segments):
        load_desc = " (Tension)" if seg['internal_load'] > 0 else " (Compression)"
        solution += f"  - P{i+1} = {seg['internal_load']} {unit_P}{load_desc}\n"

    solution += (
        f"\n**Step 2:** Calculate the Deformation (δ) of Each Segment\n"
        f"Using the formula δ = (P × L) / (A × E) for each segment.\n"
    )
    unit_L = "m" if use_si_units else "in"
    unit_D = "mm" if use_si_units else "in"
    unit_E = "GPa" if use_si_units else "ksi"

    for i, seg in enumerate(segments):
        solution += f"  - **Segment {i+1} ({seg['material_name']}):**\n"
        
        if use_si_units:
            # Explicitly show unit conversions in the solution string
            P_disp = f"{seg['internal_load']} kN ({seg['internal_load']*1000:.0f} N)"
            L_disp = f"{seg['length']} m ({seg['length']*1000:.0f} mm)"
            E_disp = f"{seg['E']} GPa ({seg['E']*1000:.0f} MPa)"
            
            solution += (
                f"P{i+1} = {P_disp}\n"
                f"L{i+1} = {L_disp}\n"
                f"A{i+1} = {round(seg['area'], precision)} mm²\n"
                f"E{i+1} = {E_disp}\n"
                f"δ{i+1} = ({seg['internal_load']*1000:.0f} × {seg['length']*1000:.0f}) / ({round(seg['area'], precision)} × {seg['E']*1000:.0f}) = {round(seg['deformation'], precision)} mm\n"
            )
        else:
            solution += (
                f"P{i+1} = {seg['internal_load']} kips\n"
                f"L{i+1} = {seg['length']} in\n"
                f"A{i+1} = {round(seg['area'], precision)} in²\n"
                f"E{i+1} = {seg['E']} ksi\n"
                f"δ{i+1} = ({seg['internal_load']} × {seg['length']}) / ({round(seg['area'], precision)} × {seg['E']}) = {round(seg['deformation'], precision)} in\n"
            )

    solution += (
        f"\n**Step 3:** Calculate the Total Deformation\n"
        f"The total deformation is the algebraic sum of the individual segment deformations.\n"
        f"δ_total = "
    )
    delta_sum_str = " + ".join([f"({round(s['deformation'], precision)})" for s in segments])
    solution += delta_sum_str + "\n"
    solution += f"δ_total = {round(total_deformation, precision)} {unit_D}\n\n"
    
    final_desc = "elongation" if total_deformation > 0 else "contraction"
    solution += (
        f"**Answer:**\n"
        f"The total deformation of the rod is **{round(total_deformation, precision)} {unit_D}** "
        f"(a net {final_desc})."
    )

    return question, solution


# Template 4 (Intermediate)
def template_poissons_ratio():
    """
    Poisson's Ratio and Change in Diameter

    Scenario:
        This template tests the understanding of Poisson's Ratio (nu). A rod is
        subjected to an axial load, causing it to elongate and contract laterally.
        The user must calculate the change in the rod's diameter.

    Core Equations:
        Axial Strain: epsilon_axial = sigma_x / E = P / (A * E)
        Lateral Strain: epsilon_lateral = -nu * epsilon_axial
        Change in Diameter: delta_d = d_initial * epsilon_lateral

    Returns:
        tuple: A tuple containing:
            - str: A question asking for the change in diameter.
            - str: A step-by-step solution.
    """
    # 1. Parameterize inputs
    use_si_units = random.choice([True, False])
    material_name = random.choice(list(MATERIAL_PROPERTIES.keys()))
    material = MATERIAL_PROPERTIES[material_name]
    precision = 5 

    # Generate dimensions first
    if use_si_units:
        initial_length_m = round(random.uniform(0.5, 3.0), 2)
        initial_diameter_mm = random.randint(25, 125)
        area_mm2 = math.pi * (initial_diameter_mm / 2)**2
        E_GPa = material['E_GPa']
        nu = material['nu']
        
        # Ensure Elastic Behavior 
        # Instead of random load, pick a safe target strain (0.1% to 0.4%)
        target_strain = random.uniform(0.001, 0.004) * random.choice([-1, 1])
        
        # Calculate Load required to achieve this strain: P = E * A * epsilon
        # E in MPa (GPa * 1000), A in mm2 -> P in Newtons
        required_load_N = (E_GPa * 1000) * area_mm2 * target_strain
        
        # Convert to kN and round to look like a "given" problem value
        load_kN = round(required_load_N / 1000.0)
        # Ensure load is not zero
        if load_kN == 0: load_kN = 1 if target_strain > 0 else -1
            
        # Now recalculate exact stress/strain from this rounded load
        stress_MPa = (load_kN * 1000) / area_mm2
        axial_strain = stress_MPa / (E_GPa * 1000)
        
        # Formatting
        load_str = f"{abs(load_kN)} kN"
        load_type_str = "tension" if load_kN > 0 else "compression"
        dim_str = f"{initial_diameter_mm} mm"
        unit_d = "mm"
        
        # Secondary calculations
        lateral_strain = -nu * axial_strain
        delta_diameter_mm = initial_diameter_mm * lateral_strain
        final_diameter_mm = initial_diameter_mm + delta_diameter_mm
        
    else:
        # US Customary Unit System
        initial_length_in = round(random.uniform(20.0, 100.0), 1)
        initial_diameter_in = round(random.uniform(1.0, 5.0), 2)
        area_in2 = math.pi * (initial_diameter_in / 2)**2
        E_ksi = material['E_ksi']
        nu = material['nu']
        
        # Ensure Elastic Behavior 
        target_strain = random.uniform(0.001, 0.004) * random.choice([-1, 1])
        
        # Calculate Load: P = E * A * epsilon
        required_load_kips = E_ksi * area_in2 * target_strain
        
        # Round to integer kips
        load_kips = round(required_load_kips)
        if load_kips == 0: load_kips = 1 if target_strain > 0 else -1
            
        # Recalculate exact values
        stress_ksi = load_kips / area_in2
        axial_strain = stress_ksi / E_ksi
        
        # Formatting
        load_str = f"{abs(load_kips)} kips"
        load_type_str = "tension" if load_kips > 0 else "compression"
        dim_str = f"{initial_diameter_in} in"
        unit_d = "in"
        
        # Secondary calculations
        lateral_strain = -nu * axial_strain
        delta_diameter_in = initial_diameter_in * lateral_strain
        final_diameter_in = initial_diameter_in + delta_diameter_in

    # 2. Generate Question and Solution Strings
    question = (
        f"A solid circular rod made of {material_name} has an initial diameter of {dim_str}. "
        f"The rod is subjected to an axial load of {load_str} ({load_type_str}).\n\n"
        f"Given the material properties for {material_name}:\n"
        f"Modulus of Elasticity (E) = {material['E_GPa' if use_si_units else 'E_ksi']} {'GPa' if use_si_units else 'ksi'}\n"
        f"Poisson's Ratio (nu) = {material['nu']}\n\n"
        f"Determine the following:\n"
        f"a) The change in the rod's diameter.\n"
        f"b) The final diameter of the rod under the load."
    )

    solution = (
        f"**Given:**\n"
        f"Material: {material_name} (E = {material['E_GPa' if use_si_units else 'E_ksi']} {'GPa' if use_si_units else 'ksi'}, nu = {nu})\n"
        f"Initial Diameter (d_0): {dim_str}\n"
        f"Axial Load (P): {load_str} ({load_type_str})\n\n"

        f"**Step 1:** Calculate Axial Strain (epsilon_axial)\n"
        f"First, we need the cross-sectional area (A) and the axial stress (sigma).\n"
        f"A = pi * (d_0 / 2)^2 = pi * ({initial_diameter_mm if use_si_units else initial_diameter_in} / 2)^2 = {round(area_mm2 if use_si_units else area_in2, 4)} {'mm^2' if use_si_units else 'in^2'}\n"
        f"\n"
    )

    if use_si_units:
        solution += (
            f"sigma = P / A = ({load_kN * 1000} N) / ({round(area_mm2, 4)} mm^2) = {round(stress_MPa, 3)} MPa\n"
            f"Now, calculate axial strain using Hooke's Law: epsilon_axial = sigma / E.\n"
            f"E = {E_GPa} GPa = {E_GPa * 1000} MPa\n"
            f"epsilon_axial = {round(stress_MPa, 3)} MPa / {E_GPa * 1000} MPa = {axial_strain:.4e}\n\n"
        )
    else: # US units
        solution += (
            f"sigma = P / A = {load_kips} kips / {round(area_in2, 4)} in^2 = {round(stress_ksi, 3)} ksi\n"
            f"Now, calculate axial strain using Hooke's Law: epsilon_axial = sigma / E.\n"
            f"epsilon_axial = {round(stress_ksi, 3)} ksi / {E_ksi} ksi = {axial_strain:.4e}\n\n"
        )

    solution += (
        f"**Step 2:** Calculate Lateral Strain (epsilon_lateral)\n"
        f"Lateral strain is related to axial strain by Poisson's ratio: epsilon_lateral = -nu * epsilon_axial.\n"
        f"The negative sign indicates that for positive axial strain (elongation), the lateral strain is negative (contraction), and vice-versa.\n"
        f"epsilon_lateral = -({nu}) * ({axial_strain:.4e}) = {lateral_strain:.4e}\n\n"
        
        f"**Step 3:** Calculate the Change in Diameter (delta_d)\n"
        f"The change in diameter is the lateral strain multiplied by the initial diameter.\n"
        f"delta_d = d_0 * epsilon_lateral\n"
        f"delta_d = ({initial_diameter_mm if use_si_units else initial_diameter_in} {unit_d}) * ({lateral_strain:.4e}) = {round(delta_diameter_mm if use_si_units else delta_diameter_in, precision)} {unit_d}\n\n"

        f"**Step 4:** Calculate the Final Diameter (d_final)\n"
        f"The final diameter is the initial diameter plus the change.\n"
        f"d_final = d_0 + delta_d\n"
        f"d_final = ({initial_diameter_mm if use_si_units else initial_diameter_in} {unit_d}) + ({round(delta_diameter_mm if use_si_units else delta_diameter_in, precision)} {unit_d}) = {round(final_diameter_mm if use_si_units else final_diameter_in, precision)} {unit_d}\n\n"

        f"**Answer:**\n"
        f"a) The change in diameter is **{round(delta_diameter_mm if use_si_units else delta_diameter_in, precision)} {unit_d}**.\n"
        f"b) The final diameter is **{round(final_diameter_mm if use_si_units else final_diameter_in, precision)} {unit_d}**."
    )
    
    return question, solution


# Template 5 (Advanced)
def template_statically_indeterminate():
    """
    Statically Indeterminate Axially Loaded Member

    Scenario:
        A uniform rod is fixed between two rigid walls. An external load is applied
        at an intermediate point. This is statically indeterminate because there are
        two unknown reaction forces (at the walls) but only one static equilibrium
        equation. The problem is solved by using a compatibility equation based on
        the fact that the total deformation of the rod must be zero.

    Core Equations:
        Equilibrium: Sum(F_x) = 0 => R_A + R_C = P
        Compatibility: delta_total = 0 => delta_AB + delta_BC = 0
                       (P_AB * L_AB / AE) + (P_BC * L_BC / AE) = 0

    Returns:
        tuple: A tuple containing:
            - str: A question asking for reaction forces and stresses.
            - str: A step-by-step solution.
    """
    # 1. Parameterize inputs
    use_si_units = random.choice([True, False])
    material_name = random.choice(list(MATERIAL_PROPERTIES.keys()))
    material = MATERIAL_PROPERTIES[material_name]
    shape = random.choice(['circular', 'square'])
    precision = 3

    if use_si_units:
        #  SI Unit System 
        total_length = round(random.uniform(1.0, 3.0), 2) # m
        len_AB = round(random.uniform(0.2 * total_length, 0.8 * total_length), 2) # m
        len_BC = total_length - len_AB
        load_P = random.randint(100, 500) # kN
        
        if shape == 'circular':
            diameter = random.randint(50, 150) # mm
            area = math.pi * (diameter / 2)**2
            dim_str = f"a diameter of {diameter} mm"
        else: # square
            side = random.randint(50, 150) # mm
            area = side**2
            dim_str = f"a side length of {side} mm"
        
        E_val = material['E_GPa']
        unit_L, unit_P, unit_S, unit_A, unit_E = "m", "kN", "MPa", "mm^2", "GPa"

    else:
        #  US Customary Unit System 
        total_length = round(random.uniform(40.0, 120.0), 1) # inches
        len_AB = round(random.uniform(0.2 * total_length, 0.8 * total_length), 1) # inches
        len_BC = total_length - len_AB
        load_P = random.randint(50, 250) # kips
        
        if shape == 'circular':
            diameter = round(random.uniform(2.0, 6.0), 2) # inches
            area = math.pi * (diameter / 2)**2
            dim_str = f"a diameter of {diameter} in"
        else: # square
            side = round(random.uniform(2.0, 6.0), 2) # inches
            area = side**2
            dim_str = f"a side length of {side} in"
            
        E_val = material['E_ksi']
        unit_L, unit_P, unit_S, unit_A, unit_E = "in", "kips", "ksi", "in^2", "ksi"

    # 2. Core Calculations
    # From compatibility: R_A = P * L_BC / L_AC
    # From equilibrium: R_C = P - R_A
    R_A = load_P * (len_BC / total_length)
    R_C = load_P - R_A

    # Internal forces
    P_AB = R_A  # Tension
    P_BC = -R_C # Compression

    # Stresses
    stress_AB = P_AB / area
    stress_BC = P_BC / area

    # Unit conversion for SI stress
    if use_si_units:
        stress_AB_MPa = stress_AB * 1000 # Convert kN/mm^2 to MPa
        stress_BC_MPa = stress_BC * 1000 # Convert kN/mm^2 to MPa

    # 3. Generate Question and Solution Strings
    question = (
        f"A solid {material_name} rod with a uniform {shape} cross-section is fixed at both ends, A and C. "
        f"The rod has a total length of {total_length} {unit_L} and {dim_str}.\n\n"
        f"An external axial load of {load_P} {unit_P} is applied to the right at point B, "
        f"which is located at a distance of {len_AB} {unit_L} from end A.\n\n"
        f"The Modulus of Elasticity for {material_name} is {E_val} {unit_E}.\n\n"
        f"Determine the following:\n"
        f"a) The reaction forces at supports A and C.\n"
        f"b) The normal stress in segments AB and BC of the rod."
    )

    solution = (
        f"**Given:**\n"
        f"Total Length (L_AC): {total_length} {unit_L}\n"
        f"Length of segment AB (L_AB): {len_AB} {unit_L}\n"
        f"Length of segment BC (L_BC): {total_length} - {len_AB} = {round(len_BC, 2)} {unit_L}\n"
        f"Applied Load (P): {load_P} {unit_P}\n"
        f"Area (A): {round(area, 2)} {unit_A}\n"
        f"Modulus of Elasticity (E): {E_val} {unit_E}\n\n"

        f"This is a statically indeterminate problem because there are two unknown support reactions (R_A and R_C) and only one equation of static equilibrium.\n\n"

        f"**Step 1:** Equation of Equilibrium\n"
        f"Summing forces in the x-direction (assuming right is positive):\n"
        f"Sum(F_x) = 0  =>  -R_A + P - R_C = 0\n"
        f"R_A + R_C = {load_P} {unit_P}  -(Equation 1)\n\n"

        f"**Step 2:** Equation of Compatibility\n"
        f"Since the rod is fixed between unyielding supports, the total deformation must be zero.\n"
        f"delta_total = delta_AB + delta_BC = 0\n"
        f"The internal force in segment AB is the reaction R_A (in tension, P_AB = R_A).\n"
        f"The internal force in segment BC is the reaction R_C (in compression, P_BC = -R_C).\n\n"
        f"Using the deformation formula delta = PL/AE:\n"
        f"(R_A * L_AB) / (A * E) - (R_C * L_BC) / (A * E) = 0\n"
        f"Since A and E are constant, they cancel out:\n"
        f"R_A * L_AB = R_C * L_BC  -(Equation 2)\n\n"

        f"**Step 3:** Solve for the Reaction Forces\n"
        f"From Equation 1, we can express R_C as: R_C = {load_P} - R_A.\n"
        f"Substitute this into the simplified Equation 2:\n"
        f"R_A * ({len_AB}) = ({load_P} - R_A) * ({round(len_BC, 2)})\n"
        f"{round(len_AB, 2)}*R_A = {round(load_P * len_BC, 2)} - {round(len_BC, 2)}*R_A\n"
        f"({round(len_AB, 2)} + {round(len_BC, 2)})*R_A = {round(load_P * len_BC, 2)}\n"
        f"({total_length})*R_A = {round(load_P * len_BC, 2)}\n"
        f"R_A = {round(load_P * len_BC, 2)} / {total_length} = {round(R_A, precision)} {unit_P}\n\n"
        f"Now, find R_C using Equation 1:\n"
        f"R_C = {load_P} - R_A = {load_P} - {round(R_A, precision)} = {round(R_C, precision)} {unit_P}\n\n"

        f"**Step 4:** Calculate the Normal Stresses\n"
        f"The stress in each segment is sigma = P_internal / A.\n"
        f"**Segment AB:** The internal force is P_AB = R_A = {round(R_A, precision)} {unit_P} (Tension).\n"
    )
    if use_si_units:
        solution += f"sigma_AB = ({round(R_A, precision)} kN) / ({round(area, 2)} mm^2) * 1000 = {round(stress_AB_MPa, precision)} MPa (Tension)\n"
    else:
        solution += f"sigma_AB = {round(R_A, precision)} kips / {round(area, 2)} in^2 = {round(stress_AB, precision)} ksi (Tension)\n"
    
    solution += f"  - **Segment BC:** The internal force is P_BC = -R_C = -{round(R_C, precision)} {unit_P} (Compression).\n"

    if use_si_units:
        solution += f"sigma_BC = (-{round(R_C, precision)} kN) / ({round(area, 2)} mm^2) * 1000 = {round(stress_BC_MPa, precision)} MPa (Compression)\n\n"
    else:
        solution += f"sigma_BC = -{round(R_C, precision)} kips / {round(area, 2)} in^2 = {round(stress_BC, precision)} ksi (Compression)\n\n"

    solution += (
        f"**Answer:**\n"
        f"a) Reaction Forces: R_A = **{round(R_A, precision)} {unit_P}** and R_C = **{round(R_C, precision)} {unit_P}**.\n"
        f"b) Normal Stresses: sigma_AB = **{round(stress_AB_MPa if use_si_units else stress_AB, precision)} {unit_S} (Tension)** and sigma_BC = **{round(abs(stress_BC_MPa if use_si_units else stress_BC), precision)} {unit_S} (Compression)**."
    )
    return question, solution


def main():
    """
    Generate numerous instances of each stres and strain - axial loading template 
    with different random seeds and write the results to a JSONL file.
    """
    import json
    import os

    # Define the output path (Modify this path according to where you are running the code from)
    output_file = "testset/mechanical_engineering/mechanics_of_materials/stress_and_strain.jsonl"

    # Create the directory if it doesn't exist
    os.makedirs(os.path.dirname(output_file), exist_ok=True)

    # List of template functions with their ID and level
    templates = [
        (template_basic_stress_strain, "basic_stress_strain", "Easy"),
        (template_axial_deformation, "axial_deformation", "Easy"),
        (template_multi_segment_rod, "multi_segment_rod", "Intermediate"),
        (template_poissons_ratio, "poissons_ratio", "Intermediate"),
        (template_statically_indeterminate, "statically_indeterminate", "Advanced"),
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
                "domain": "mechanics_of_materials",
                "area": "stress_and_strain",
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
