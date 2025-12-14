import random
import math


# Template 1 (Easy)
def template_signal_operations():
    """
    Signal Operations: Shifting and Reversal

    Scenario:
        This template tests the fundamental understanding of transformations on the
        independent variable (time index 'n'). It requires applying a time-shift
        (x[n-n0]) or a time-reversal (x[-n]) to a given finite-length sequence.
        This is a core skill for understanding more complex operations like convolution.

    Core Equations:
        1. Time-Shift: y[n] = x[n - n0]
        2. Time-Reversal: z[n] = x[-n]

    Returns:
        tuple: A tuple containing:
            - str: A question asking to perform a time-shift or time-reversal on a sequence.
            - str: A step-by-step solution showing the transformation of the sequence.
    """
    # 1. Parameterize the inputs with random values
    
    # Generate a sequence x[n] as a dictionary {index: value}
    seq_length = random.randint(4, 7)
    # Ensure the origin (n=0) is not at the absolute ends of the sequence
    origin_position = random.randint(1, seq_length - 2)
    start_index = -origin_position
    
    x_n = {start_index + i: random.randint(-10, 10) for i in range(seq_length)}
    
    # Choose the operation
    operation = random.choice(['shift', 'reversal'])

    # 2. Perform the core calculation and generate explanatory text
    
    if operation == 'shift':
        # Select a non-zero shift amount
        n0 = random.choice([-3, -2, -1, 1, 2, 3])
        output_var = "y"
        
        # Core calculation: Shift the keys of the dictionary
        result_n = {k + n0: v for k, v in x_n.items()}
        
        # Explanations for the solution
        op_str_symbolic = f"x[n - ({n0})]" if n0 < 0 else f"x[n - {n0}]"
        
        direction = "left" if n0 < 0 else "right"
        explanation = (
            f"The operation is a time shift, {output_var}[n] = {op_str_symbolic}.\n"
            f"This corresponds to a shift of the sequence to the {direction} by {abs(n0)} sample(s).\n"
            f"The value at any original index 'k' is moved to a new index 'k + {n0}'."
        )

    else: # operation == 'reversal'
        n0 = 0 # Not used, but needed for variable scope
        output_var = "z"
        
        # Core calculation: Negate the keys of the dictionary
        result_n = {-k: v for k, v in x_n.items()}
        
        # Explanations for the solution
        op_str_symbolic = "x[-n]"
        explanation = (
            f"The operation is a time reversal, {output_var}[n] = {op_str_symbolic}.\n"
            f"This corresponds to flipping the sequence about the origin (n=0).\n"
            f"The value at any original index 'k' is moved to a new index '-k'."
        )

    # Inlined logic to format the input sequence x_n into a string
    if not x_n:
        x_n_str = "{}"
    else:
        min_idx = min(x_n.keys())
        max_idx = max(x_n.keys())
        parts = []
        for i in range(min_idx, max_idx + 1):
            val = x_n.get(i, 0)
            if i == 0:
                parts.append(f"*{val}*") # Asterisk denotes the origin n=0
            else:
                parts.append(str(val))
        x_n_str = f"{{{', '.join(parts)}}}"

    # Inlined logic to format the result sequence into a string
    if not result_n:
        result_n_str = "{}"
    else:
        min_idx = min(result_n.keys())
        max_idx = max(result_n.keys())
        parts = []
        for i in range(min_idx, max_idx + 1):
            val = result_n.get(i, 0)
            if i == 0:
                parts.append(f"*{val}*")
            else:
                parts.append(str(val))
        result_n_str = f"{{{', '.join(parts)}}}"

    # 3. Generate the question and solution strings
    
    question = (
        f"A discrete-time signal is defined by the sequence:\n"
        f"x[n] = {x_n_str}\n\n"
        f"Determine the resulting sequence, {output_var}[n], after applying the following transformation:\n"
        f"{output_var}[n] = {op_str_symbolic}"
    )
    
    # Inlined logic to build the transformation table for the solution
    table_header = f"""| {'Original Index (k)':<20} | {'Original Value x[k]':<22} | {"New Index (k')":<18} | {"New Value y[k')":<20} |"""
    table_separator = "-" * len(table_header)

    table_rows_list = []
    for k in sorted(x_n.keys()):
        val = x_n[k]
        if operation == 'shift':
            new_k = k + n0
        else: # reversal
            new_k = -k
        table_rows_list.append(f"| {k:<20} | {val:<22} | {new_k:<18} | {val:<20} |")

    # Combine all parts of the table with newlines
    transformation_table = f"{table_header}\n{table_separator}\n" + "\n".join(table_rows_list)

    solution = (
        f"**Given:**\n"
        f"The original sequence is x[n] = {x_n_str}.\n"
        f"The operation is {output_var}[n] = {op_str_symbolic}.\n\n"
        
        f"**Step 1:** Understand the Transformation\n"
        f"{explanation}\n\n"
        
        f"**Step 2:** Apply the Transformation to Each Index\n"
        f"We can create a table to track where each value moves:\n\n"
        f"{transformation_table}\n\n"
        
        f"**Step 3:** Construct the Final Sequence\n"
        f"By collecting the values at their new indices, we get the final sequence.\n"
        f"Remember that any index not explicitly calculated has a value of 0.\n\n"
        
        f"**Answer:**\n"
        f"The resulting sequence is {output_var}[n] = {result_n_str}"
    )

    return question, solution


# Template 2 (Easy)
def template_system_properties_memory_causality():
    """
    System Properties: Memory and Causality 

    Scenario:
        This template tests the ability to analyze a discrete-time system's
        input-output equation to determine two fundamental properties: whether it
        is memoryless and whether it is causal. This requires applying the formal
        definitions of these properties to the given equation.

    Core Definitions:
        1. Memoryless: Output y[n] depends only on the current input x[n].
        2. Causal: Output y[n] depends only on present and past inputs x[k], where k <= n.

    Returns:
        tuple: A tuple containing:
            - str: A question asking to determine if a system is memoryless and causal.
            - str: A step-by-step solution analyzing each property with justifications.
    """
    # 1. Parameterize the inputs by generating a random system type and equation
    
    system_type = random.choice(['memoryless_causal', 'memory_causal', 'memory_noncausal'])
    
    is_memoryless = False
    is_causal = False
    equation_str = ""
    memory_reason = ""
    causality_reason = ""

    if system_type == 'memoryless_causal':
        is_memoryless = True
        is_causal = True
        
        form = random.choice(['power', 'affine', 'scaled_by_n'])
        if form == 'power':
            power = random.randint(2, 3)
            equation_str = f"(x[n])**{power}"
            memory_reason = f"The output y[n] is the input x[n] raised to the power of {power} at the exact same time 'n'. The system does not need to store any other input values."
        elif form == 'affine':
            gain = random.randint(2, 9)
            offset = random.randint(1, 10)
            equation_str = f"{gain} * x[n] + {offset}"
            memory_reason = f"The output y[n] at time 'n' is a scaled and shifted version of the input x[n] at the exact same time 'n'. No past or future values of the input are needed."
        else: # scaled_by_n
            equation_str = f"n * x[n]"
            memory_reason = "The output y[n] depends on the input x[n] at the same time 'n', scaled by the time index 'n' itself. The time index 'n' is not an input signal value that needs to be stored."

        causality_reason = "The output y[n] depends only on the input at the present time 'n'. Since no future inputs (k > n) are required, the system is causal."

    elif system_type == 'memory_causal':
        is_memoryless = False
        is_causal = True
        
        delay = random.randint(1, 4)
        form = random.choice(['simple_delay', 'with_present_term'])
        
        if form == 'simple_delay':
            equation_str = f"x[n - {delay}]"
            memory_reason = f"To compute the output y[n], the system needs to know the input from {delay} time step(s) in the past, x[n - {delay}]. Therefore, the system must have memory."
        else: # with_present_term
            equation_str = f"x[n] + x[n - {delay}]"
            memory_reason = f"The output y[n] depends on both the current input x[n] and a past input x[n - {delay}]. The system must store or 'remember' the value of x[n - {delay}] to calculate the current output."

        causality_reason = f"The output y[n] depends on the present input and/or past inputs (up to time n - {delay}). Since it does not depend on any future inputs, the system is causal."

    else: # system_type == 'memory_noncausal'
        is_memoryless = False
        is_causal = False
        
        advance = random.randint(1, 4)
        form = random.choice(['simple_advance', 'with_past_term'])

        if form == 'simple_advance':
            equation_str = f"x[n + {advance}]"
            causality_reason = f"The output y[n] depends on the input at a future time, x[n + {advance}]. To calculate y[n], the system must know what the input will be {advance} time step(s) in the future. This violates the definition of causality."
        else: # with_past_term
            delay = random.randint(1, 2)
            equation_str = f"x[n - {delay}] + x[n + {advance}]"
            causality_reason = f"The output y[n] depends on a future input, x[n + {advance}]. Even though it also depends on a past input, the dependency on a future value makes the system non-causal."
        
        memory_reason = "The output y[n] depends on input values at times other than 'n'. Therefore, the system is not memoryless."

    # 2. Generate the question and solution strings
    
    question = (
        f"A discrete-time system is described by the input-output equation:\n"
        f"y[n] = {equation_str}\n\n"
        f"Determine if the system is:\n"
        f"a) Memoryless\n"
        f"b) Causal\n\n"
        f"Justify your answers based on the definitions of these properties."
    )

    solution = (
        f"**Given:**\n"
        f"The system's input-output equation is y[n] = {equation_str}.\n\n"

        f"**Step 1:** Analyze the Memoryless Property\n"
        f"**Definition:** A system is **memoryless** if its output at any time 'n' depends *only* on the input at the same time 'n'.\n"
        f"**Analysis:** {memory_reason}\n"
        f"**Conclusion:** Based on the analysis, the system is **{'memoryless' if is_memoryless else 'not memoryless'}**.\n\n"

        f"**Step 2:** Analyze the Causality Property\n"
        f"**Definition:** A system is **causal** if its output at any time 'n' depends *only* on the input at the present time 'n' and past times (i.e., on x[k] for k <= n).\n"
        f"**Analysis:** {causality_reason}\n"
        f"**Conclusion:** Based on the analysis, the system is **{'causal' if is_causal else 'not causal'}**.\n\n"

        f"**Answer:**\n"
        f"a) Memoryless: **{'Yes' if is_memoryless else 'No'}**\n"
        f"b) Causal: **{'Yes' if is_causal else 'No'}**"
    )

    return question, solution


# Template 3 (Intermediate)
def template_finite_convolution():
    """
    Convolution of Finite-Length Sequences 

    Scenario:
        This template tests the direct application of the convolution sum, which is
        the fundamental operation for determining the output of a Linear
        Time-Invariant (LTI) system.

    Core Equations:
        1. Convolution Sum: y[n] = sum(x[k] * h[n-k]) for all k

    Returns:
        tuple: A tuple containing:
            - str: A question asking to find the output of an LTI system.
            - str: A step-by-step solution demonstrating the convolution process.
    """
    # 1. Parameterize the inputs with random values

    # Generate sequence x[n]
    x_len = random.randint(3, 4)
    x_origin_pos = random.randint(0, x_len - 1)
    x_start_idx = -x_origin_pos
    x_n = {x_start_idx + i: random.randint(-3, 3) for i in range(x_len)}
    
    # Ensure the sequence isn't all zeros
    if all(v == 0 for v in x_n.values()):
        x_n[random.choice(list(x_n.keys()))] = random.randint(1, 3)

    # Generate sequence h[n]
    h_len = random.randint(3, 4)
    h_origin_pos = random.randint(0, h_len - 1)
    h_start_idx = -h_origin_pos
    h_n = {h_start_idx + i: random.randint(-2, 2) for i in range(h_len)}
    
    if all(v == 0 for v in h_n.values()):
        h_n[random.choice(list(h_n.keys()))] = random.randint(1, 2)

    # Format x[n] string
    min_idx_x = min(x_n.keys())
    max_idx_x = max(x_n.keys())
    x_parts = []
    for i in range(min_idx_x, max_idx_x + 1):
        val = x_n.get(i, 0)
        x_parts.append(f"*{val}*" if i == 0 else str(val))
    x_n_str = f"{{{', '.join(x_parts)}}}"

    # Format h[n] string
    min_idx_h = min(h_n.keys())
    max_idx_h = max(h_n.keys())
    h_parts = []
    for i in range(min_idx_h, max_idx_h + 1):
        val = h_n.get(i, 0)
        h_parts.append(f"*{val}*" if i == 0 else str(val))
    h_n_str = f"{{{', '.join(h_parts)}}}"


    # 2. Perform the core calculation (Convolution)
    y_n = {}
    y_start_idx = min_idx_x + min_idx_h
    y_end_idx = max_idx_x + max_idx_h

    all_k_indices = sorted(list(x_n.keys()))
    
    for n in range(y_start_idx, y_end_idx + 1):
        current_sum = 0
        for k in all_k_indices:
            x_k = x_n.get(k, 0)
            h_nk = h_n.get(n - k, 0)
            current_sum += x_k * h_nk
        # Store non-zero values (or boundary zeros if needed, but dict sparse is fine)
        if current_sum != 0:
            y_n[n] = current_sum
    
    # 3. Generate the question and solution strings
    
    question = (
        f"An LTI system has an impulse response h[n] = {h_n_str}.\n\n"
        f"Determine the system's output, y[n] = x[n] * h[n], for the input x[n] = {x_n_str}."
    )

    # Build the detailed calculation steps for the solution
    calculation_steps = []
    y_indices_to_show = sorted(y_n.keys()) if y_n else [y_start_idx]
    
    # To avoid too much output, show first, middle, and last
    if len(y_indices_to_show) > 3:
        middle_index = y_indices_to_show[len(y_indices_to_show)//2]
        y_indices_to_show = [y_indices_to_show[0], middle_index, y_indices_to_show[-1]]

    for n in y_indices_to_show:
        step_str = f"For n = {n}:\n"
        step_str += f"y[{n}] = sum( x[k] * h[{n}-k] )\n"
        
        sum_expr_terms = []
        val_expr_terms = []
        
        # Iterate over ALL valid k in x[n] to show full expansion, ensuring consistency
        # Removed the "if x_val != 0" check to show all terms explicitly
        for k in sorted(x_n.keys()):
            x_val = x_n[k]
            h_val = h_n.get(n - k, 0)
            
            # Identify h index for clarity
            h_idx = n - k
            
            sum_expr_terms.append(f"x[{k}]h[{h_idx}]")
            val_expr_terms.append(f"({x_val})({h_val})")
        
        step_str += f"y[{n}] = {' + '.join(sum_expr_terms)}\n"
        step_str += f"y[{n}] = {' + '.join(val_expr_terms)}\n"
        step_str += f"y[{n}] = {y_n.get(n, 0)}\n"
        calculation_steps.append(step_str)

    # Format the final sequence y_n
    if not y_n:
        y_n_str = "{*0*}"
    else:
        min_idx_y = y_start_idx
        max_idx_y = y_end_idx
        y_parts = []
        for i in range(min_idx_y, max_idx_y + 1):
            val = y_n.get(i, 0)
            y_parts.append(f"*{val}*" if i == 0 else str(val))
        y_n_str = f"{{{', '.join(y_parts)}}}"

    y_values_list = [f"    y[{n}] = {y_n.get(n, 0)}" for n in range(y_start_idx, y_end_idx + 1)]
    y_values_str = "\n".join(y_values_list)

    calculation_steps_str = "\n".join(calculation_steps)
        
    solution = (
        f"**Given:**\n"
        f"Input Signal: x[n] = {x_n_str}\n"
        f"Impulse Response: h[n] = {h_n_str}\n\n"

        f"**Step 1:** State the Convolution Formula\n"
        f"The output y[n] of an LTI system is the convolution of the input x[n] with the impulse response h[n]. The convolution sum is defined as:\n"
        f"y[n] = sum over all k of (x[k] * h[n-k])\n\n"

        f"**Step 2:** Apply the Flip-and-Slide Method\n"
        f"We can visualize this process by flipping the impulse response h[k] to get h[-k], and then sliding it by 'n' positions. For each slide 'n', we calculate the sum of the products of the overlapping samples.\n"
        f"\n\n"
        f"Let's calculate a few points explicitly:\n\n"
        f"{calculation_steps_str}\n"
        
        f"**Step 3: Calculate All Output Values**\n"
        f"By continuing this process for all values of 'n' where the sequences overlap (from n={y_start_idx} to n={y_end_idx}), we get the full output sequence:\n"
        f"{y_values_str}\n\n" 

        f"**Answer:**\n"
        f"The complete output sequence is y[n] = {y_n_str}"
    )

    return question, solution


# Template 4 (Intermediate)
def template_system_property_linearity():
    """
    System Properties: Linearity

    Scenario:
        This template tests the ability to formally prove or disprove if a system
        is linear by checking the two defining properties: additivity and
        homogeneity (scaling).

    Core Definitions:
        1. Additivity: T{x1[n] + x2[n]} = T{x1[n]} + T{x2[n]}
        2. Homogeneity: T{a * x[n]} = a * T{x[n]}

    Returns:
        tuple: A tuple containing:
            - str: A question asking to determine if a system is linear.
            - str: A step-by-step solution showing the formal proof.
    """
    # 1. Parameterize the inputs
    
    system_type = random.choice(['linear_gain', 'linear_delay', 'nonlinear_offset', 'nonlinear_power'])
    
    equation_str = ""
    is_linear = False
    
    # We will define these explicitly for each case to avoid "replace" bugs
    y1_str = ""
    y2_str = ""
    
    # Strings for the RHS of the proof steps (excluding "y = " prefix)
    add_sum_str = ""
    add_y3_str = ""
    add_comparison_str = ""
    
    hom_ay1_str = ""
    hom_ya_str = ""
    hom_comparison_str = ""

    if system_type == 'linear_gain':
        is_linear = True
        k = random.randint(1, 10)
        equation_str = f"{k} * x[n]"
        
        y1_str = f"{k} * x1[n]"
        y2_str = f"{k} * x2[n]"
        
        # Additivity
        add_sum_str = f"({k} * x1[n]) + ({k} * x2[n]) = {k} * (x1[n] + x2[n])"
        add_y3_str = f"{k} * (x3[n]) = {k} * (x1[n] + x2[n])"
        add_comparison_str = "Since y3[n] is equal to y1[n] + y2[n], the system satisfies the additivity property."
        
        # Homogeneity
        hom_ay1_str = f"a * ({k} * x1[n])"
        hom_ya_str = f"{k} * (xa[n]) = {k} * (a * x1[n])"
        hom_comparison_str = "Since ya[n] is equal to a * y1[n], the system satisfies the homogeneity property."
        
    elif system_type == 'linear_delay':
        is_linear = True
        d = random.randint(1, 10)
        equation_str = f"x[n - {d}]"
        
        # Explicitly constructing strings handles the 'n-d' index correctly
        y1_str = f"x1[n - {d}]"
        y2_str = f"x2[n - {d}]"

        # Additivity
        add_sum_str = f"x1[n - {d}] + x2[n - {d}]"
        add_y3_str = f"x3[n - {d}] = x1[n - {d}] + x2[n - {d}]"
        add_comparison_str = "Since y3[n] is equal to y1[n] + y2[n], the system satisfies the additivity property."

        # Homogeneity
        hom_ay1_str = f"a * x1[n - {d}]"
        hom_ya_str = f"xa[n - {d}] = a * x1[n - {d}]"
        hom_comparison_str = "Since ya[n] is equal to a * y1[n], the system satisfies the homogeneity property."
        
    elif system_type == 'nonlinear_offset':
        is_linear = False
        C = random.randint(1, 5) * random.choice([-1, 1])
        C_str = f"+ {C}" if C > 0 else f"- {abs(C)}"
        equation_str = f"x[n] {C_str}"
        
        y1_str = f"x1[n] {C_str}"
        y2_str = f"x2[n] {C_str}"
        
        # Additivity
        add_sum_str = f"(x1[n] {C_str}) + (x2[n] {C_str}) = x1[n] + x2[n] + {2*C}"
        add_y3_str = f"(x3[n]) {C_str} = (x1[n] + x2[n]) {C_str}"
        add_comparison_str = f"Since x1[n] + x2[n] + {2*C} is not equal to x1[n] + x2[n] {C_str}, the system fails the additivity test."

        # Homogeneity
        hom_ay1_str = f"a * (x1[n] {C_str}) = a*x1[n] + {C}*a"
        hom_ya_str = f"(xa[n]) {C_str} = (a * x1[n]) {C_str}"
        hom_comparison_str = f"Since a*x1[n] + {C}*a is not equal to a*x1[n] {C_str} (for a != 1), the system fails the homogeneity test."

    elif system_type == 'nonlinear_power':
        is_linear = False
        p = random.randint(2, 3)
        equation_str = f"(x[n])^{p}"
        
        y1_str = f"(x1[n])^{p}"
        y2_str = f"(x2[n])^{p}"
        
        # Additivity
        add_sum_str = f"(x1[n])^{p} + (x2[n])^{p}"
        add_y3_str = f"(x3[n])^{p} = (x1[n] + x2[n])^{p}"
        add_comparison_str = f"In general, (x1[n] + x2[n])^{p} is not equal to (x1[n])^{p} + (x2[n])^{p}. Therefore, the system fails the additivity test."

        # Homogeneity
        hom_ay1_str = f"a * (x1[n])^{p}"
        hom_ya_str = f"(xa[n])^{p} = (a * x1[n])^{p} = (a^{p}) * (x1[n])^{p}"
        hom_comparison_str = f"Since a * (x1[n])^{p} is not equal to (a^{p}) * (x1[n])^{p} (for a != 1), the system fails the homogeneity test."
        
    # 2. Generate the question and solution strings
    
    question = (
        f"A discrete-time system is governed by the equation:\n"
        f"y[n] = {equation_str}\n\n"
        f"Determine if this system is linear by testing the properties of additivity and homogeneity."
    )

    solution = (
        f"**Given:**\n"
        f"The system equation is y[n] = {equation_str}.\n\n"
        
        f"**Step 1:** State the Conditions for Linearity\n"
        f"For a system to be linear, it must satisfy two properties:\n"
        f"1. **Additivity:** T{{x1[n] + x2[n]}} = T{{x1[n]}} + T{{x2[n]}}\n"
        f"2. **Homogeneity (Scaling):** T{{a*x[n]}} = a*T{{x[n]}}\n"
        f"We must test both properties.\n"
        f"\n\n"

        f"**Step 2:** Test for Additivity\n"
        f"Let's define two arbitrary inputs, x1[n] and x2[n]. The corresponding outputs are:\n"
        f"y1[n] = {y1_str}\n"
        f"y2[n] = {y2_str}\n\n"
        f"The sum of these outputs is:\n"
        f"y1[n] + y2[n] = {add_sum_str}\n\n"
        f"Now, let's define a third input x3[n] = x1[n] + x2[n]. The output y3[n] is:\n"
        f"y3[n] = {add_y3_str}\n\n"
        f"**Comparison:** {add_comparison_str}\n\n"

        f"**Step 3:** Test for Homogeneity (Scaling)\n"
        f"Let's define an input x1[n] and a constant 'a'. The output is y1[n] = {y1_str}.\n"
        f"The scaled output is:\n"
        f"a * y1[n] = {hom_ay1_str}\n\n"
        f"Now, let's define a new input xa[n] = a * x1[n]. The output ya[n] is:\n"
        f"ya[n] = {hom_ya_str}\n\n"
        f"**Comparison:** {hom_comparison_str}\n\n"
        
        f"**Answer:**\n"
        f"{'The system satisfies both additivity and homogeneity, therefore, the system is **linear**.' if is_linear else 'The system fails at least one of the tests (additivity or homogeneity), therefore, the system is **not linear**.'}"
    )
    
    return question, solution


# Template 5 (Advanced)
def template_impulse_response_from_lccde():
    """
    Finding the Impulse Response from a Difference Equation (Advanced)

    Scenario:
        This template tests the ability to find the impulse response h[n] for a
        system described by a Linear Constant-Coefficient Difference Equation
        (LCCDE). The process involves setting the input to the unit impulse,
        delta[n], and solving the resulting recurrence relation for h[n].

    Core Definitions:
        1. Impulse Response: h[n] = T{delta[n]}
        2. Causality: h[n] = 0 for n < 0

    Returns:
        tuple: A tuple containing:
            - str: A question asking for the impulse response.
            - str: A step-by-step solution detailing the recursive method.
    """
    # 1. Parameterize by randomly choosing system order and coefficients
    
    order = random.choice(['first', 'second'])
    
    # Generate coefficients with small integer values for clarity
    b0 = random.randint(1, 4) * random.choice([-1, 1])
    a1 = random.randint(1, 5) * random.choice([-1, 1])
    
    # Helper lambda to format coefficients into strings like "+ 3" or "- 2"
    fmt = lambda c, v: f"+ {c}" if c > 0 else f"- {abs(c)}"
    
    if order == 'first':
        # For first order, we can have b0*x[n] and b1*x[n-1] terms
        b1 = random.randint(1, 4) * random.choice([-1, 1])
        
        # Build equation strings
        y_terms = f"y[n] {fmt(a1, 'y[n-1]')}y[n-1]"
        # Randomly omit b1 term for variety
        if random.random() < 0.4:
            b1 = 0
        x_terms = f"{b0}x[n]"
        if b1 != 0:
            x_terms += f" {fmt(b1, 'x[n-1]')}x[n-1]"
        
        equation_str = f"{y_terms} = {x_terms}"
        
        # 2. Perform the core calculation for a first-order system
        h0 = b0
        h1 = b1 - a1 * h0
        
        term1_str = f"{h0}*delta[n]"
        
        h_decay_expr = ""
        if h1 != 0:
            base = -a1
            # Use fmt helper to correctly sign the term
            h_decay_expr = f" {fmt(h1, '')}({base})^(n-1) * u[n-1]"
        
        final_h_n = f"{term1_str}{h_decay_expr}"
        
        # 3. Generate the solution string for a first-order system
        solution_steps = (
            f"**Step 3:** Solve Recursively for Initial Conditions\n"
            f"We assume the system is causal, so **h[n] = 0 for n < 0**.\n\n"
            f"**For n = 0:**\n"
            f"h[0] {fmt(a1, 'h[-1]')}h[-1] = {b0}*delta[0] {'+ 0' if b1==0 else ' ' + fmt(b1, 'delta[-1]')+'*delta[-1]'}\n"
            f"h[0] {fmt(a1, '0')}*(0) = {b0}*(1) + {b1}*(0)\n"
            f"**h[0] = {h0}**\n\n"
            
            f"**For n = 1:**\n"
            f"h[1] {fmt(a1, 'h[0]')}h[0] = {b0}*delta[1] {fmt(b1, 'delta[0]')}*delta[0]\n"
            f"h[1] {fmt(a1, h0)}*({h0}) = {b0}*(0) + {b1}*(1)\n"
            f"h[1] = {b1} - ({a1*h0})\n"
            f"**h[1] = {h1}**\n\n"
            
            f"**Step 4:** Find the Homogeneous Solution\n"
            f"For n >= 2, the input delta terms are zero. The equation becomes homogeneous:\n"
            f"h[n] {fmt(a1, 'h[n-1]')}h[n-1] = 0  =>  h[n] = {-a1}*h[n-1]\n\n"
            f"The solution to this recurrence for n >= 1 is a decaying exponential that starts at n=1 with value h[1].\n"
            f"This part of the response can be written as h[1]*({-a1})^(n-1)*u[n-1].\n\n"
            
            f"**Step 5:** Combine Results for the Final Expression\n"
            f"The total impulse response is the sum of the value at n=0 and the response for n >= 1:\n"
            f"h[n] = h[0]*delta[n] + h[1]*({-a1})^(n-1)*u[n-1]\n"
        )

    else: # order == 'second'
        # To guarantee real, distinct roots: D = a1^2 - 4*a2 > 0
        a2 = random.randint(-5, 5)
        if a2 == 0: a2 = 1 # Avoid trivial case
        while a1**2 - 4*a2 <= 0:
            a1 = random.randint(-6, 6)
            if a1 == 0: a1 = 1
        
        # Randomize RHS
        b1 = random.randint(-3, 3) if random.random() > 0.4 else 0
        # Force b2 = 0 to ensure the coefficient fitting method is valid
        # If b2 != 0, there is an impulse at n=2 that breaks the homogeneous assumption for n>=2
        b2 = 0 

        # Build equation strings
        x_terms = f"{b0}x[n]"
        if b1 != 0: x_terms += f" {fmt(b1, 'x[n-1]')}x[n-1]"
        # b2 term removed
        equation_str = f"y[n] {fmt(a1, 'y[n-1]')}y[n-1] {fmt(a2, 'y[n-2]')}y[n-2] = {x_terms}"
        
        # 2. Perform the core calculation for a second-order system
        h0 = b0
        h1 = b1 - a1 * h0
        
        # Correctly solve r^2 + a1*r + a2 = 0
        discriminant = a1**2 - 4*a2
        r1 = (-a1 + math.sqrt(discriminant)) / 2
        r2 = (-a1 - math.sqrt(discriminant)) / 2
        r1_str, r2_str = f"{r1:.2f}", f"{r2:.2f}"
        
        # Correctly solve for C1 and C2 with general h[0], h[1]
        if abs(r1 - r2) < 1e-9: # Failsafe for very close roots
            r1, r2 = round(r1, 2), round(r2, 2)
        C1 = (h1 - h0 * r2) / (r1 - r2)
        C2 = (h0 * r1 - h1) / (r1 - r2)
        C1_str, C2_str = f"{C1:.2f}", f"{C2:.2f}"
        
        final_h_n = f"({C1_str}({r1_str})^n {fmt(C2, C2_str)}({r2_str})^n) * u[n]"
        
        # 3. Generate the solution string for a second-order system
        h_eq = f"h[n] {fmt(a1, 'h[n-1]')}h[n-1] {fmt(a2, 'h[n-2]')}h[n-2]"
        d_eq = f"{b0}*delta[n]"
        if b1 != 0: d_eq += f" {fmt(b1, 'd[n-1]')}*delta[n-1]"

        solution_steps = (
            f"**Step 3:** Solve Recursively for Initial Conditions\n"
            f"We assume the system is causal, so **h[n] = 0 for n < 0**.\n\n"
            f"**For n = 0:**\n"
            f"h[0] {fmt(a1, 'h[-1]')}h[-1] {fmt(a2, 'h[-2]')}h[-2] = {b0}*delta[0] ... (other delta terms are 0)\n"
            f"h[0] {fmt(a1, '0')}*(0) {fmt(a2, '0')}*(0) = {b0}*(1)\n"
            f"**h[0] = {h0}**\n\n"
            f"**For n = 1:**\n"
            f"h[1] {fmt(a1, 'h[0]')}h[0] {fmt(a2, 'h[-1]')}h[-1] = ... {fmt(b1, 'd[0]')}*delta[0] ...\n"
            f"h[1] {fmt(a1, h0)}*({h0}) {fmt(a2, '0')}*(0) = {b0}*(0) + {b1}*(1)\n"
            f"h[1] = {b1} - {a1*h0}\n"
            f"**h[1] = {h1}**\n\n"

            f"**Step 4:** Find the Homogeneous Solution\n"
            f"For n >= 2, the input is zero (since b2=0), so the equation becomes homogeneous:\n"
            f"h[n] {fmt(a1, 'h[n-1]')}h[n-1] {fmt(a2, 'h[n-2]')}h[n-2] = 0\n\n"
            f"We solve this by finding the roots of the characteristic equation: r^2 {fmt(a1, 'r')}r {fmt(a2, '')} = 0\n"
            f"Using the quadratic formula, the roots are r1 = {r1_str}, r2 = {r2_str}.\n"
            f"The general solution for n >= 0 is h[n] = C1*({r1_str})^n + C2*({r2_str})^n.\n\n"
            
            f"**Step 5:** Use Initial Conditions to Find Coefficients\n"
            f"We use h[0] and h[1] to create a system of two equations:\n"
            f"1) For n=0: h[0] = C1 + C2  =>  {h0} = C1 + C2\n"
            f"2) For n=1: h[1] = C1*({r1_str}) + C2*({r2_str})  =>  {h1} = {r1_str}*C1 + {r2_str}*C2\n\n"
            f"Solving this system yields:\n"
            f"**C1 = {C1_str}** and **C2 = {C2_str}**\n"
        )
        
    question = (
        f"A causal LTI system is described by the difference equation:\n"
        f"{equation_str}\n\n"
        f"Find the impulse response h[n] of the system."
    )
    
    solution = (
        f"**Given:**\n"
        f"The system equation is {equation_str}\n\n"

        f"**Step 1:** Set Input to the Unit Impulse\n"
        f"By definition, the impulse response h[n] is the output y[n] when the input x[n] is the unit impulse, delta[n].\n"
        f"\n\n"

        f"**Step 2:** Substitute h[n] and delta[n] into the Equation\n"
        f"Replacing y[n] with h[n] and x[n] with delta[n], we get:\n"
        f"{equation_str.replace('y', 'h').replace('x', 'delta')}\n\n"
        
        f"{solution_steps}\n"
        
        f"**Answer:**\n"
        f"Substituting the coefficients back into the general form, the impulse response is:\n"
        f"h[n] = {final_h_n}"
    )

    return question, solution


def main():
    """
    Generate numerous instances of each discrete time signals template 
    with different random seeds and write the results to a JSONL file.
    """
    import json
    import os

    # Define the output path (Modify this path according to where you are running the code from)
    output_file = "testset/electrical_engineering/signals_and_systems/discrete_time_signals.jsonl"

    # Create the directory if it doesn't exist
    os.makedirs(os.path.dirname(output_file), exist_ok=True)

    # List of template functions with their ID and level
    templates = [
        (template_signal_operations, "signal_operations", "Easy"),
        (template_system_properties_memory_causality, "system_properties_memory_causality", "Easy"),
        (template_finite_convolution, "finite_convolution", "Intermediate"),
        (template_system_property_linearity, "system_property_linearity", "Intermediate"),
        (template_impulse_response_from_lccde, "impulse_response_from_lccde", "Advanced"),
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
                "domain": "signals_and_systems",
                "area": "discrete_time_signals",
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
