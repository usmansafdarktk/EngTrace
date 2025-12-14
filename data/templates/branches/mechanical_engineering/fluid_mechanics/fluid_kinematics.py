import random
import math


# Template 1 (Easy)
def template_fluid_particle_acceleration():
    """
    Fluid Kinematics: Fluid Particle Acceleration in a 2D Field

    Scenario:
        This template generates a problem that tests the ability to calculate the
        acceleration of a fluid particle at a specific point and time. It requires
        applying the material derivative to a given 2D unsteady velocity field.
        The problem distinguishes between the local (time-based) and convective
        (space-based) components of acceleration.

    Core Equations:
        Velocity Field: V(x, y, t) = u(x, y, t)i + v(x, y, t)j
        Acceleration (ax) = du/dt + u*(du/dx) + v*(du/dy)
        Acceleration (ay) = dv/dt + u*(dv/dx) + v*(dv/dy)
        Total Acceleration (a) = sqrt(ax^2 + ay^2)

    Returns:
        tuple: A tuple containing:
            - str: A question asking for the acceleration components and magnitude.
            - str: A step-by-step solution to the problem.
    """
    # 1. Parameterize the inputs with random values for high diversity
    
    # Coefficients for the velocity field polynomials
    # u = A*x*t + B*y**2
    # v = C*x*y + D*t**2
    A = random.randint(1, 10)
    B = random.randint(1, 10)
    C = random.randint(1, 10)
    D = random.randint(1, 10)

    # Specific point (x, y) and time (t) for evaluation
    x_point = round(random.uniform(0.5, 5.0), 1)
    y_point = round(random.uniform(0.5, 5.0), 1)
    t_point = round(random.uniform(0.1, 3.0), 1)

    # Standardize precision for all calculations and outputs
    precision = 3

    # 2. Perform the core calculations for the solution

    # Step A: Evaluate velocity components u and v at the given point and time
    u_val = A * x_point * t_point + B * y_point**2
    v_val = C * x_point * y_point + D * t_point**2

    # Step B: Calculate all necessary partial derivatives
    # For u(x, y, t) = A*x*t + B*y^2
    du_dt = A * x_point
    du_dx = A * t_point
    du_dy = 2 * B * y_point

    # For v(x, y, t) = C*x*y + D*t^2
    dv_dt = 2 * D * t_point
    dv_dx = C * y_point
    dv_dy = C * x_point
    
    # Step C: Calculate the acceleration components (ax and ay)
    # ax = (local acceleration in x) + (convective acceleration in x)
    local_ax = du_dt
    convective_ax = u_val * du_dx + v_val * du_dy
    ax = local_ax + convective_ax

    # ay = (local acceleration in y) + (convective acceleration in y)
    local_ay = dv_dt
    convective_ay = u_val * dv_dx + v_val * dv_dy
    ay = local_ay + convective_ay

    # Step D: Calculate the magnitude of the total acceleration
    a_magnitude = math.sqrt(ax**2 + ay**2)


    # 3. Generate the question and solution strings
    
    question = (
        f"The velocity field of a fluid is given by the equations:\n"
        f"u = {A}xt + {B}y^2\n"
        f"v = {C}xy + {D}t^2\n"
        f"where u and v are in m/s, and x and y are in meters, and t is in seconds.\n\n"
        f"Determine the acceleration components (ax and ay) and the magnitude of the "
        f"acceleration for a fluid particle at the point ({x_point}, {y_point}) m at time t = {t_point} s."
    )

    solution = (
        f"**Given:**\n"
        f"Velocity component u = {A}xt + {B}y^2\n"
        f"Velocity component v = {C}xy + {D}t^2\n"
        f"Point of interest: (x, y) = ({x_point}, {y_point}) m\n"
        f"Time of interest: t = {t_point} s\n\n"

        f"**Step 1:** Calculate the acceleration component in the x-direction (ax).\n"
        f"The formula is: ax = du/dt + u*(du/dx) + v*(du/dy)\n\n"
        
        f"First, find the required partial derivatives of u:\n"
        f"du/dt = d/dt({A}xt + {B}y^2) = {A}x  -> At ({x_point}, {t_point}), du/dt = {A} * {x_point} = {round(du_dt, precision)}\n"
        f"du/dx = d/dx({A}xt + {B}y^2) = {A}t  -> At ({x_point}, {t_point}), du/dx = {A} * {t_point} = {round(du_dx, precision)}\n"
        f"du/dy = d/dy({A}xt + {B}y^2) = {2*B}y  -> At ({y_point}), du/dy = {2*B} * {y_point} = {round(du_dy, precision)}\n\n"

        f"Next, evaluate u and v at the specified point and time:\n"
        f"u = {A}({x_point})({t_point}) + {B}({y_point})^2 = {round(u_val, precision)} m/s\n"
        f"v = {C}({x_point})({y_point}) + {D}({t_point})^2 = {round(v_val, precision)} m/s\n\n"

        f"Now, substitute these values into the acceleration formula for ax:\n"
        f"ax = (du/dt) + u*(du/dx) + v*(du/dy)\n"
        f"ax = {round(du_dt, precision)} + ({round(u_val, precision)})*({round(du_dx, precision)}) + ({round(v_val, precision)})*({round(du_dy, precision)})\n"
        f"ax = {round(local_ax, precision)} + {round(convective_ax, precision)} = {round(ax, precision)} m/s^2\n\n"

        f"**Step 2:** Calculate the acceleration component in the y-direction (ay).\n"
        f"The formula is: ay = dv/dt + u*(dv/dx) + v*(dv/dy)\n\n"

        f"First, find the required partial derivatives of v:\n"
        f"dv/dt = d/dt({C}xy + {D}t^2) = {2*D}t  -> At ({t_point}), dv/dt = {2*D} * {t_point} = {round(dv_dt, precision)}\n"
        f"dv/dx = d/dx({C}xy + {D}t^2) = {C}y  -> At ({y_point}), dv/dx = {C} * {y_point} = {round(dv_dx, precision)}\n"
        f"dv/dy = d/dy({C}xy + {D}t^2) = {C}x  -> At ({x_point}), dv/dy = {C} * {x_point} = {round(dv_dy, precision)}\n\n"

        f"Using the already calculated u and v values:\n"
        f"ay = (dv/dt) + u*(dv/dx) + v*(dv/dy)\n"
        f"ay = {round(dv_dt, precision)} + ({round(u_val, precision)})*({round(dv_dx, precision)}) + ({round(v_val, precision)})*({round(dv_dy, precision)})\n"
        f"ay = {round(local_ay, precision)} + {round(convective_ay, precision)} = {round(ay, precision)} m/s^2\n\n"

        f"**Step 3:** Calculate the magnitude of the total acceleration.\n"
        f"The formula is: a = sqrt(ax^2 + ay^2)\n"
        f"a = sqrt(({round(ax, precision)})^2 + ({round(ay, precision)})^2)\n"
        f"a = {round(a_magnitude, precision)} m/s^2\n\n"
        
        f"**Answer:**\n"
        f"The acceleration components are ax = {round(ax, precision)} m/s^2 and ay = {round(ay, precision)} m/s^2.\n"
        f"The magnitude of the acceleration is {round(a_magnitude, precision)} m/s^2."
    )

    return question, solution


# Template 2 (Intermediate)
def template_volumetric_flow_rate():
    """
    Fluid Kinematics: Volumetric Flow Rate from a Velocity Profile

    Scenario:
        This template assesses the ability to calculate the volumetric flow rate (Q)
        by integrating a given velocity profile over a cross-sectional area. It also
        calculates the average velocity. The problem randomly selects between flow
        in a circular pipe (parabolic profile) and flow in a rectangular channel
        (linear profile).

    Core Equations:
        Volumetric Flow Rate: Q = integral(u dA)
        Average Velocity: V_avg = Q / A
        Pipe Profile: u(r) = U_max * (1 - (r/R)^2), dA = 2*pi*r dr
        Channel Profile: u(y) = U_max * (y/H), dA = W dy

    Returns:
        tuple: A tuple containing:
            - str: A question asking for the flow rate and average velocity.
            - str: A step-by-step solution to the problem.
    """
    # 1. Parameterize the inputs with random values
    geometry = random.choice(['pipe', 'channel'])
    u_max = round(random.uniform(0.5, 10.0), 2)
    precision = 4

    # 2. Perform calculations and generate strings based on the chosen geometry
    if geometry == 'pipe':
        #  Pipe-specific parameters 
        diameter_mm = random.randint(20, 200)
        radius_m = diameter_mm / 2000.0

        #  Core calculations for the pipe 
        area = math.pi * radius_m**2
        # For a parabolic profile in a pipe, the exact integral yields Q = (1/2) * U_max * A
        flow_rate = 0.5 * u_max * area
        avg_velocity = u_max / 2.0

        #  Generate question and solution strings for the pipe 
        question = (
            f"The velocity profile for a fluid flowing through a circular pipe with a diameter of {diameter_mm} mm "
            f"is given by u(r) = U_max * (1 - (r/R)^2), where R is the radius of the pipe and the maximum "
            f"velocity U_max is {u_max} m/s at the centerline (r=0).\n\n"
            f"Determine the volumetric flow rate (Q) and the average velocity (V_avg) of the fluid."
        )

        solution = (
            f"**Given:**\n"
            f"Geometry: Circular Pipe\n"
            f"Diameter (D) = {diameter_mm} mm\n"
            f"Maximum Velocity (U_max) = {u_max} m/s\n"
            f"Velocity Profile: u(r) = U_max * (1 - (r/R)^2)\n\n"

            f"**Step 1:** Determine pipe dimensions in meters.\n"
            f"Radius (R) = D / 2 = {diameter_mm} / 2 = {diameter_mm/2.0} mm = {radius_m} m\n\n"

            f"**Step 2:** Set up the integral for volumetric flow rate (Q).\n"
            f"The formula is Q = integral(u dA) over the cross-sectional area.\n"
            f"For a circular pipe, the differential area is a thin ring: dA = 2*pi*r dr.\n"
            f"Q = integral from 0 to R of [U_max * (1 - (r/R)^2)] * (2*pi*r dr)\n"
            f"Q = 2*pi*U_max * integral from 0 to R of (r - r^3/R^2) dr\n\n"

            f"**Step 3:** Evaluate the integral.\n"
            f"Q = 2*pi*U_max * [r^2/2 - r^4/(4*R^2)] from 0 to R\n"
            f"Q = 2*pi*U_max * [(R^2/2 - R^4/(4*R^2)) - 0]\n"
            f"Q = 2*pi*U_max * [R^2/4] = (pi*R^2*U_max) / 2\n\n"

            f"**Step 4:** Substitute numerical values to find Q.\n"
            f"Q = (pi * ({radius_m})^2 * {u_max}) / 2\n"
            f"Q = {round(flow_rate, precision)} m^3/s\n\n"

            f"**Step 5:** Calculate the average velocity (V_avg).\n"
            f"First, calculate the cross-sectional area (A) = pi*R^2\n"
            f"A = pi * ({radius_m})^2 = {round(area, precision)} m^2\n"
            f"V_avg = Q / A = {round(flow_rate, precision)} / {round(area, precision)}\n"
            f"Alternatively, for this profile, V_avg is known to be U_max / 2.\n"
            f"V_avg = {u_max} / 2 = {round(avg_velocity, precision)} m/s\n\n"

            f"**Answer:**\n"
            f"The volumetric flow rate is {round(flow_rate, precision)} m^3/s, and the average velocity is {round(avg_velocity, precision)} m/s."
        )

    else: # geometry == 'channel'
        #  Channel-specific parameters 
        height_cm = random.randint(5, 50)
        width_cm = random.randint(10, 100)
        height_m = height_cm / 100.0
        width_m = width_cm / 100.0

        #  Core calculations for the channel 
        area = width_m * height_m
        # For a linear profile in a channel, the exact integral yields Q = (1/2) * U_max * A
        flow_rate = 0.5 * u_max * area
        avg_velocity = u_max / 2.0

        #  Generate question and solution strings for the channel 
        question = (
            f"A fluid flows through a rectangular channel that is {width_cm} cm wide and {height_cm} cm high. "
            f"The velocity profile is given by u(y) = U_max * (y/H), where H is the channel height and the "
            f"maximum velocity U_max is {u_max} m/s at the top surface (y=H).\n\n"
            f"Determine the volumetric flow rate (Q) and the average velocity (V_avg) of the fluid."
        )

        solution = (
            f"**Given:**\n"
            f"Geometry: Rectangular Channel\n"
            f"Width (W) = {width_cm} cm\n"
            f"Height (H) = {height_cm} cm\n"
            f"Maximum Velocity (U_max) = {u_max} m/s\n"
            f"Velocity Profile: u(y) = U_max * (y/H)\n\n"

            f"**Step 1:** Determine channel dimensions in meters.\n"
            f"Width (W) = {width_cm} cm = {width_m} m\n"
            f"Height (H) = {height_cm} cm = {height_m} m\n\n"

            f"**Step 2:** Set up the integral for volumetric flow rate (Q).\n"
            f"The formula is Q = integral(u dA) over the cross-sectional area.\n"
            f"For a rectangular channel, the differential area is a thin horizontal strip: dA = W dy.\n"
            f"Q = integral from 0 to H of [U_max * (y/H)] * (W dy)\n"
            f"Q = (W*U_max/H) * integral from 0 to H of (y) dy\n\n"

            f"**Step 3:** Evaluate the integral.\n"
            f"Q = (W*U_max/H) * [y^2/2] from 0 to H\n"
            f"Q = (W*U_max/H) * [(H^2/2) - 0]\n"
            f"Q = (W*U_max*H) / 2\n\n"

            f"**Step 4:** Substitute numerical values to find Q.\n"
            f"Q = ({width_m} * {u_max} * {height_m}) / 2\n"
            f"Q = {round(flow_rate, precision)} m^3/s\n\n"

            f"**Step 5:** Calculate the average velocity (V_avg).\n"
            f"First, calculate the cross-sectional area (A) = W * H\n"
            f"A = {width_m} * {height_m} = {round(area, precision)} m^2\n"
            f"V_avg = Q / A = {round(flow_rate, precision)} / {round(area, precision)}\n"
            f"Alternatively, for this profile, V_avg is known to be U_max / 2.\n"
            f"V_avg = {u_max} / 2 = {round(avg_velocity, precision)} m/s\n\n"

            f"**Answer:**\n"
            f"The volumetric flow rate is {round(flow_rate, precision)} m^3/s, and the average velocity is {round(avg_velocity, precision)} m/s."
        )

    return question, solution


# Template 3 (Intermediate)
def template_vorticity_check():
    """
    Fluid Kinematics: Vorticity and Rotational Flow Check

    Scenario:
        This template introduces the concept of fluid element rotation. It requires
        calculating the curl of a given velocity field to find the vorticity vector.
        This enhanced version uses more complex velocity fields where the vorticity
        itself depends on the spatial coordinates, making the evaluation point critical
        to the final answer. The problem randomly generates either a 2D or 3D field.

    Core Equations:
        Vorticity Vector: zeta = curl(V) = (dw/dy - dv/dz)i + (du/dz - dw/dx)j + (dv/dx - du/dy)k
        Irrotational Flow: A flow is irrotational if the vorticity vector is the zero vector.

    Returns:
        tuple: A tuple containing:
            - str: A question asking for the vorticity vector and a rotationality check.
            - str: A step-by-step solution to the problem.
    """
    # 1. Parameterize the inputs with random values
    flow_type = random.choice(['2D', '3D'])
    precision = 3

    # Random coefficients for the velocity field polynomials.
    # Allowing zeros increases variety, sometimes creating simpler or irrotational cases.
    A = random.randint(-4, 4)
    B = random.randint(-4, 4)
    C = random.randint(-4, 4)

    # Evaluation point
    x_point = round(random.uniform(0.5, 3.0), 1)
    y_point = round(random.uniform(0.5, 3.0), 1)
    
    # 2. Perform calculations and generate strings based on the flow type
    if flow_type == '2D':
        #  Define a more complex 2D velocity field 
        # u = A*x*y, v = B*x^2 + C*y^2
        # This makes derivatives dependent on the point (x, y).

        #  Core calculations for 2D flow 
        du_dy_val = A * x_point
        dv_dx_val = 2 * B * x_point
        
        zeta_z = dv_dx_val - du_dy_val
        is_irrotational = abs(zeta_z) < 1e-9 # Check for zero with floating point tolerance

        #  Generate question and solution strings for 2D flow 
        question = (
            f"A 2D fluid flow is described by the velocity field:\n"
            f"u = {A}xy\n"
            f"v = {B}x^2 + {C}y^2\n"
            f"where u and v are in m/s, and x and y are in meters.\n\n"
            f"Calculate the vorticity at the point ({x_point}, {y_point}) m and "
            f"determine if the flow is rotational or irrotational at that point."
        )

        solution = (
            f"**Given:**\n"
            f"Velocity component u = {A}xy\n"
            f"Velocity component v = {B}x^2 + {C}y^2\n"
            f"Point of interest: (x, y) = ({x_point}, {y_point}) m\n\n"

            f"**Step 1:** State the formula for vorticity in a 2D flow.\n"
            f"For a 2D flow in the x-y plane, the vorticity vector is (0, 0, zeta_z).\n"
            f"The formula for the z-component is: zeta_z = (dv/dx - du/dy)\n\n"

            f"**Step 2:** Calculate the necessary partial derivatives.\n"
            f"  du/dy = d/dy({A}xy) = {A}x\n"
            f"  dv/dx = d/dx({B}x^2 + {C}y^2) = {2*B}x\n\n"

            f"**Step 3:** Evaluate the partial derivatives at the point ({x_point}, {y_point}).\n"
            f"  du/dy at x={x_point} -> {A}({x_point}) = {round(du_dy_val, precision)}\n"
            f"  dv/dx at x={x_point} -> {2*B}({x_point}) = {round(dv_dx_val, precision)}\n\n"

            f"**Step 4:** Substitute the derivative values to find the vorticity component.\n"
            f"  zeta_z = ({round(dv_dx_val, precision)}) - ({round(du_dy_val, precision)}) = {round(zeta_z, precision)}\n"
            f"The vorticity vector at this point is (0, 0, {round(zeta_z, precision)}) rad/s.\n\n"

            f"**Step 5:** Determine if the flow is rotational at this point.\n"
            f"A flow is irrotational at a point if its vorticity is zero there.\n"
            f"Since the vorticity component zeta_z = {round(zeta_z, precision)}, which is {'zero' if is_irrotational else 'not zero'},\n"
            f"the flow is {'irrotational' if is_irrotational else 'rotational'} at this specific point.\n\n"
            
            f"**Answer:**\n"
            f"The vorticity at ({x_point}, {y_point}) is (0, 0, {round(zeta_z, precision)}) rad/s. The flow is {'irrotational' if is_irrotational else 'rotational'} at this point."
        )

    else: # flow_type == '3D'
        D = random.randint(-4, 4)
        E = random.randint(-4, 4)
        F = random.randint(-4, 4)
        z_point = round(random.uniform(0.5, 3.0), 1)

        #  Define a more complex 3D velocity field 
        # u = A*y^2, v = B*z^2 + C*x, w = D*x^2 + E*y
        # This allows all components of vorticity to be non-zero and position-dependent.

        #  Core calculations for 3D flow 
        du_dy_val = 2 * A * y_point
        du_dz_val = 0
        
        dv_dx_val = C
        dv_dz_val = 2 * B * z_point
        
        dw_dx_val = 2 * D * x_point
        dw_dy_val = E
        
        zeta_x = dw_dy_val - dv_dz_val
        zeta_y = du_dz_val - dw_dx_val
        zeta_z = dv_dx_val - du_dy_val
        
        is_irrotational = (abs(zeta_x) < 1e-9 and abs(zeta_y) < 1e-9 and abs(zeta_z) < 1e-9)

        #  Generate question and solution strings for 3D flow 
        question = (
            f"A 3D fluid flow is described by the velocity field:\n"
            f"u = {A}y^2\n"
            f"v = {B}z^2 + {C}x\n"
            f"w = {D}x^2 + {E}y\n"
            f"where u, v, w are in m/s, and x, y, z are in meters.\n\n"
            f"Calculate the vorticity vector at the point ({x_point}, {y_point}, {z_point}) m and "
            f"determine if the flow is rotational or irrotational at that point."
        )

        solution = (
            f"**Given:**\n"
            f"Velocity field: u={A}y^2, v={B}z^2 + {C}x, w={D}x^2 + {E}y\n"
            f"Point of interest: ({x_point}, {y_point}, {z_point}) m\n\n"

            f"**Step 1:** State the formula for the 3D vorticity vector.\n"
            f"zeta = (dw/dy - dv/dz)i + (du/dz - dw/dx)j + (dv/dx - du/dy)k\n\n"

            f"**Step 2:** Calculate all necessary partial derivatives.\n"
            f"From u = {A}y^2: du/dy = {2*A}y, du/dz = 0\n"
            f"From v = {B}z^2 + {C}x: dv/dx = {C}, dv/dz = {2*B}z\n"
            f"From w = {D}x^2 + {E}y: dw/dx = {2*D}x, dw/dy = {E}\n\n"
            
            f"**Step 3:** Evaluate the derivatives at the point ({x_point}, {y_point}, {z_point}).\n"
            f"  du/dy = {2*A}({y_point}) = {round(du_dy_val, precision)}\n"
            f"  dv/dx = {C}\n"
            f"  dv/dz = {2*B}({z_point}) = {round(dv_dz_val, precision)}\n"
            f"  dw/dx = {2*D}({x_point}) = {round(dw_dx_val, precision)}\n"
            f"  dw/dy = {E}\n"
            f"  (du/dz is always 0)\n\n"

            f"**Step 4:** Calculate each component of the vorticity vector.\n"
            f"zeta_x = dw/dy - dv/dz = {E} - ({round(dv_dz_val, precision)}) = {round(zeta_x, precision)}\n"
            f"zeta_y = du/dz - dw/dx = 0 - ({round(dw_dx_val, precision)}) = {round(zeta_y, precision)}\n"
            f"zeta_z = dv/dx - du/dy = {C} - ({round(du_dy_val, precision)}) = {round(zeta_z, precision)}\n\n"
            
            f"**Step 5:** Assemble the vorticity vector and determine if the flow is rotational.\n"
            f"The vorticity vector at this point is ({round(zeta_x, precision)}, {round(zeta_y, precision)}, {round(zeta_z, precision)}) rad/s.\n"
            f"A flow is irrotational at a point if its vorticity is the zero vector (0, 0, 0).\n"
            f"Since the vector is {'the zero vector' if is_irrotational else 'not the zero vector'},\n"
            f"the flow is {'irrotational' if is_irrotational else 'rotational'} at this specific point.\n\n"

            f"**Answer:**\n"
            f"The vorticity vector at ({x_point}, {y_point}, {z_point}) is ({round(zeta_x, precision)}, {round(zeta_y, precision)}, {round(zeta_z, precision)}) rad/s. "
            f"The flow is {'irrotational' if is_irrotational else 'rotational'} at this point."
        )

    return question, solution


# Template 4 (Advanced)
def template_particle_pathline():
    """
    Fluid Kinematics: Finding a Particle's Pathline in a Steady Flow

    Scenario:
        This template distinguishes between the Eulerian description (velocity field)
        and the Lagrangian description (individual particle path). It requires solving
        a system of ordinary differential equations (ODEs).

    Core Equations:
        Pathline Definition: dx/dt = u(x, y) and dy/dt = v(x, y)
        Velocity Field Used: u = A*x, v = -B*y

    Returns:
        tuple: A tuple containing:
            - str: A question asking for the particle's final coordinates.
            - str: A step-by-step solution showing the derivation.
    """
    # 1. Parameterize the inputs with physically realistic values
    # We reduce the magnitude of A and B to ensure the exponential growth doesn't
    # produce physically absurd velocities (e.g. keeping v < 100 m/s).
    # Units for A and B are 1/s.
    A = round(random.uniform(0.1, 0.5), 1)
    B = round(random.uniform(0.1, 0.5), 1)

    x0 = round(random.uniform(0.5, 2.0), 1)
    y0 = round(random.uniform(0.5, 2.0), 1)
    tf = round(random.uniform(1.0, 4.0), 1)

    precision = 3

    # 2. Perform the core calculations for the solution
    # x(t) = x0 * exp(A*t)
    # y(t) = y0 * exp(-B*t)
    xf = x0 * math.exp(A * tf)
    yf = y0 * math.exp(-B * tf)

    # 3. Generate the question and solution strings
    question = (
        f"A steady, 2D velocity field is defined by the equations u = {A}x and v = -{B}y, "
        f"where u and v are in m/s, x and y are in meters, and the constants have units of s^-1.\n\n"
        f"A fluid particle is located at the initial position (x0, y0) = ({x0} m, {y0} m) at time t = 0 s.\n"
        f"Determine the particle's coordinates (x, y) at time t = {tf} s."
    )

    solution = (
        f"**Given:**\n"
        f"Velocity field: u = {A}x, v = -{B}y\n"
        f"Initial position (x0, y0) = ({x0}, {y0})\n"
        f"Final time (tf) = {tf} s\n\n"

        f"**Step 1:** Set up the differential equation for the x-coordinate.\n"
        f"The pathline is defined by the rate of change of the particle's position: dx/dt = u.\n"
        f"  dx/dt = {A}x\n"
        f"\n\n"

        f"**Step 2:** Solve the ODE for x(t) by separating variables.\n"
        f"  (1/x) dx = {A} dt\n"
        f"Integrating both sides gives:\n"
        f"  ln(x) = {A}t + C1\n"
        f"Solving for x by exponentiating:\n"
        f"  x(t) = exp({A}t + C1) = C * exp({A}t)\n\n"

        f"**Step 3:** Apply the initial condition x(0) = {x0} to find the constant C.\n"
        f"  {x0} = C * exp({A}*0) = C * 1\n"
        f"So, C = {x0}. The specific solution for the x-coordinate is:\n"
        f"  x(t) = {x0} * exp({A}t)\n\n"

        f"**Step 4:** Set up and solve the differential equation for the y-coordinate.\n"
        f"The pathline is defined by dy/dt = v.\n"
        f"  dy/dt = -{B}y\n"
        f"Separating variables: (1/y) dy = -{B} dt\n"
        f"Integrating both sides gives:\n"
        f"  ln(y) = -{B}t + D1\n"
        f"  y(t) = D * exp(-{B}t)\n\n"

        f"**Step 5:** Apply the initial condition y(0) = {y0} to find the constant D.\n"
        f"  {y0} = D * exp(-{B}*0) = D * 1\n"
        f"So, D = {y0}. The specific solution for the y-coordinate is:\n"
        f"  y(t) = {y0} * exp(-{B}t)\n\n"

        f"**Step 6:** Calculate the final position at t = {tf} s.\n"
        f"Substitute t = {tf} into the pathline equations:\n"
        f"  x({tf}) = {x0} * exp({A} * {tf}) = {round(xf, precision)} m\n"
        f"  y({tf}) = {y0} * exp(-{B} * {tf}) = {round(yf, precision)} m\n\n"

        f"**Answer:**\n"
        f"At t = {tf} s, the particle's coordinates are ({round(xf, precision)} m, {round(yf, precision)} m)."
    )

    return question, solution


# Template 5 (Advanced)
def template_incompressible_continuity():
    """
    Fluid Kinematics: Deriving a Velocity Component for Incompressible Flow

    Scenario:
        This problem uses the differential form of the conservation of mass
        (continuity equation) for a 2D, incompressible flow. Given one velocity
        component, the user must find the other by performing partial differentiation
        and integration. The problem asks for the simplest possible expression.

    Core Equation:
        2D Incompressible Continuity: du/dx + dv/dy = 0

    Returns:
        tuple: A tuple containing:
            - str: A question asking for the unknown velocity component.
            - str: A step-by-step solution showing the derivation.
    """
    # Helper to format algebraic terms cleanly (e.g. handles + -3x -> - 3x)
    def fmt_term(coeff, var, is_first=False):
        if coeff == 0:
            return ""
        
        # Determine sign prefix
        if coeff < 0:
            sign = "-" if is_first else "- "
            val = -coeff
        else:
            sign = "" if is_first else "+ "
            val = coeff
            
        # Format number (remove .0 if integer)
        if isinstance(val, float) and val.is_integer():
            val = int(val)
        
        # Construct string
        # If coeff is 1 and there is a variable, omit the '1'
        if var:
            if val == 1:
                return f"{sign}{var}"
            else:
                return f"{sign}{val}{var}"
        else:
            return f"{sign}{val}"

    # Helper to join two terms into an expression
    def build_poly(c1, v1, c2, v2):
        t1 = fmt_term(c1, v1, is_first=True)
        t2 = fmt_term(c2, v2, is_first=False)
        return f"{t1} {t2}".strip()

    # 1. Parameterize the inputs with random values
    given_component = random.choice(['u', 'v'])
    
    # Random non-zero coefficients
    A = random.choice([i for i in range(-5, 6) if i != 0])
    B = random.choice([i for i in range(-5, 6) if i != 0])

    # 2. Perform calculations and generate strings based on the chosen component
    if given_component == 'u':
        #  The u-component is given, find v 
        # u = Ax^2 + Bxy
        u_expr = build_poly(A, "x^2", B, "xy")
        
        # du/dx = 2Ax + By
        du_dx_c1, du_dx_c2 = 2*A, B
        du_dx_expr = build_poly(du_dx_c1, "x", du_dx_c2, "y")
        
        # dv/dy = -du/dx = -2Ax - By
        dv_dy_c1, dv_dy_c2 = -2*A, -B
        dv_dy_expr = build_poly(dv_dy_c1, "x", dv_dy_c2, "y")
        
        # v = integral(-2Ax - By) dy = -2Axy - (B/2)y^2
        v_c1, v_c2 = -2*A, -B/2.0
        v_expr = build_poly(v_c1, "xy", v_c2, "y^2")
        
        question = (
            f"A 2D, steady, incompressible flow has a velocity component in the x-direction "
            f"given by u = {u_expr}.\n\n"
            f"Determine the simplest possible expression for the y-component of velocity, v(x, y)."
        )

        solution = (
            f"**Given:**\n"
            f"The flow is 2D, steady, and incompressible.\n"
            f"The u-component of velocity is u = {u_expr}\n\n"
            
            f"**Step 1:** State the 2D incompressible continuity equation.\n"
            f"The equation is: du/dx + dv/dy = 0\n\n"
            
            
            f"**Step 2:** Calculate the partial derivative of the given u-component with respect to x.\n"
            f"  du/dx = d/dx({u_expr})\n"
            f"  du/dx = {du_dx_expr}\n\n"
            
            f"**Step 3:** Rearrange the continuity equation to solve for dv/dy.\n"
            f"  dv/dy = -du/dx\n"
            f"  dv/dy = -({du_dx_expr}) = {dv_dy_expr}\n\n"
            
            f"**Step 4:** Integrate dv/dy with respect to y to find v(x, y).\n"
            f"  v(x, y) = integral({dv_dy_expr}) dy\n"
            f"  v(x, y) = {v_expr} + f(x)\n"
            f"The 'constant' of integration, f(x), is an arbitrary function of x.\n\n"
            
            f"**Step 5:** Provide the simplest form for v(x, y).\n"
            f"To find the simplest expression, we set the integration constant f(x) to zero.\n"
            f"  v(x, y) = {v_expr}\n\n"
            
            f"**Answer:**\n"
            f"The simplest expression for the y-component of velocity is v = {v_expr}."
        )

    else: # given_component == 'v'
        #  The v-component is given, find u 
        # v = Ay^2 + Bxy
        v_expr = build_poly(A, "y^2", B, "xy")

        # dv/dy = 2Ay + Bx
        dv_dy_c1, dv_dy_c2 = 2*A, B
        dv_dy_expr = build_poly(dv_dy_c1, "y", dv_dy_c2, "x")
        
        # du/dx = -dv/dy = -2Ay - Bx
        # Usually written as x term first: -Bx - 2Ay
        du_dx_c1, du_dx_c2 = -B, -2*A
        du_dx_expr = build_poly(du_dx_c1, "x", du_dx_c2, "y")
        
        # u = integral(-Bx - 2Ay) dx = -(B/2)x^2 - 2Axy
        u_c1, u_c2 = -B/2.0, -2*A
        u_expr_sol = build_poly(u_c1, "x^2", u_c2, "xy")

        question = (
            f"A 2D, steady, incompressible flow has a velocity component in the y-direction "
            f"given by v = {v_expr}.\n\n"
            f"Determine the simplest possible expression for the x-component of velocity, u(x, y)."
        )

        solution = (
            f"**Given:**\n"
            f"The flow is 2D, steady, and incompressible.\n"
            f"The v-component of velocity is v = {v_expr}\n\n"
            
            f"**Step 1:** State the 2D incompressible continuity equation.\n"
            f"The equation is: du/dx + dv/dy = 0\n\n"
            
            f"**Step 2:** Calculate the partial derivative of the given v-component with respect to y.\n"
            f"  dv/dy = d/dy({v_expr})\n"
            f"  dv/dy = {dv_dy_expr}\n\n"
            
            f"**Step 3:** Rearrange the continuity equation to solve for du/dx.\n"
            f"  du/dx = -dv/dy\n"
            f"  du/dx = -({dv_dy_expr}) = {du_dx_expr}\n\n"
            
            f"**Step 4:** Integrate du/dx with respect to x to find u(x, y).\n"
            f"  u(x, y) = integral({du_dx_expr}) dx\n"
            f"  u(x, y) = {u_expr_sol} + g(y)\n"
            f"The 'constant' of integration, g(y), is an arbitrary function of y.\n\n"
            
            f"**Step 5:** Provide the simplest form for u(x, y).\n"
            f"To find the simplest expression, we set the integration constant g(y) to zero.\n"
            f"  u(x, y) = {u_expr_sol}\n\n"
            
            f"**Answer:**\n"
            f"The simplest expression for the x-component of velocity is u = {u_expr_sol}."
        )

    return question, solution


def main():
    """
    Generate numerous instances of each fluid kinematics template 
    with different random seeds and write the results to a JSONL file.
    """
    import json
    import os

    # Define the output path (Modify this path according to where you are running the code from)
    output_file = "testset/mechanical_engineering/fluid_mechanics/fluid_kinematics.jsonl"

    # Create the directory if it doesn't exist
    os.makedirs(os.path.dirname(output_file), exist_ok=True)

    # List of template functions with their ID and level
    templates = [
        (template_fluid_particle_acceleration, "fluid_particle_acceleration", "Easy"),
        (template_volumetric_flow_rate, "volumetric_flow_rate", "Intermediate"),
        (template_vorticity_check, "vorticity_check", "Intermediate"),
        (template_particle_pathline, "particle_pathline", "Advanced"),
        (template_incompressible_continuity, "incompressible_continuity", "Advanced"),
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
                "area": "fluid_kinematics",
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
