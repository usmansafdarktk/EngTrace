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

    Screen pass 1 (2026-09-23):
        Two judges' findings, both about the printed sequences. The asterisk
        that marks n = 0 was never defined in the question; the question now
        says so. The answer dropped the marker whenever a shift moved the
        support off the origin (D-050 measured 16.4%; 79/500 seeds here), so
        the printed support is widened to include n = 0, the fix D-050 named,
        and Step 3 says the origin is marked. The transformation, the table
        and the sampled values are unchanged; the question text changes on
        every seed (one added clause) and the gold answer on the seeds whose
        support excluded n = 0. The "(D4.5)" the other judge saw is a source
        comment, not emitted text (0/500 instances contain it).

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
            f"The value at any original index 'k' is moved to a new index 'k {signed_term(n0)}'."
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
        # The printed support always includes n = 0, so the origin marker is
        # present even when a shift moves every sample off it (D-050).
        min_idx = min(min(result_n.keys()), 0)
        max_idx = max(max(result_n.keys()), 0)
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
        f"x[n] = {x_n_str}\n"
        f"where the value at the origin n = 0 is enclosed in asterisks and x[n] = 0 outside the listed indices.\n\n"
        f"Determine the resulting sequence, {output_var}[n], after applying the following transformation:\n"
        f"{output_var}[n] = {op_str_symbolic}"
    )
    
    # Inlined logic to build the transformation table for the solution
    # D4.5: the header read `New Value y[k')` -- it opened a bracket and closed
    # a parenthesis, and it named `y` on the reversal branch, whose output
    # variable is `z`.  The sibling columns fix both conventions: an index is
    # parenthesised (`Original Index (k)`), a value is bracketed and carries its
    # own variable (`Original Value x[k]`).
    new_value_col = f"New Value {output_var}[k']"
    table_header = f"""| {'Original Index (k)':<20} | {'Original Value x[k]':<22} | {"New Index (k')":<18} | {new_value_col:<20} |"""
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
        f"Remember that any index not explicitly calculated has a value of 0, and the origin n = 0 is marked with asterisks.\n\n"
        
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
            equation_str = f"{gain} * x[n] {signed_term(offset)}"
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

    Layer 2 fix (2026-09-25):
        A branch expert rejected the template as unsolvable as posed: the
        origins of x[n] and h[n] are drawn at random, but the question
        marked n = 0 only by wrapping one sample in asterisks, a convention
        it never defined (and one that renders as italics in markdown). A
        reader placing both sequences at n = 0 gets the right values at the
        wrong indices (seed 2101: y[0]..y[5] for the gold y[-3]..y[2]). The
        question now defines the convention in the words
        template_signal_operations uses, states the two origin samples
        explicitly (h[0] = ..., x[0] = ...), and asks for y[n] in the same
        notation; Step 3 restates it before the answer. No sampled value,
        formula or gold value changes; the question text changes on every
        seed and the gold answer on none.

    Layer 2 fix (2026-09-26):
        The same expert rejected the round-1 wording (the other two
        approved): the origin was still carried by asterisks, which the
        review app's Markdown renders as italics, and the fallback "so h[0]
        = v and x[0] = v" cannot locate n = 0 when v repeats. Verified on
        the seeds the expert saw: 2201 has x = {2, *-2*, -2}, and a reader
        taking the last -2 as x[0] gets y[0] = 2 for the gold 8; 2203 and
        2204 repeat the origin value too. At HEAD it repeats within its own
        sequence on 59.4% of 500 seeds. Models read the raw asterisks, so
        only a human reader was misled; the convolution values were right.
        Every sequence (question, Given, answer) is now followed by its
        index list, "x[n] = {2, *-2*, -2} for n = -1, 0, 1"; the asterisk
        stays as a redundant marker, called "also marked" so the sentence
        is true rendered or not, and the h[0]/x[0] sentence is gone. Step 3
        names the start index of y[n]. Sampling, arithmetic and gold values
        are unchanged: the question and the answer block change on 100% of
        500 seeds, the answer only by the appended index list, and the gold
        index-to-value map is identical on 500 of 500.

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

    # Each sequence is printed with its index list, which locates n = 0 in
    # words: Markdown renders a starred value as italics and drops the
    # asterisks, so the marker must never carry the origin alone.
    x_idx_str = ", ".join(str(i) for i in range(min_idx_x, max_idx_x + 1))
    h_idx_str = ", ".join(str(i) for i in range(min_idx_h, max_idx_h + 1))


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
    
    # Both sequences contain n = 0 by construction (the origin position is
    # drawn inside each support), so the origin samples are always present.
    question = (
        f"An LTI system has an impulse response h[n] = {h_n_str} for n = {h_idx_str}.\n\n"
        f"Determine the system's output, y[n] = x[n] * h[n] (the convolution of x[n] with h[n]), "
        f"when the input is x[n] = {x_n_str} for n = {x_idx_str}.\n\n"
        f"In each sequence the listed values belong, in order, to the listed indices n, "
        f"and the sequence is 0 at every other n; the value at n = 0 is also marked. "
        f"Give y[n] in the same form, with its indices."
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
        y_list_idx = [0]
    else:
        min_idx_y = y_start_idx
        max_idx_y = y_end_idx
        y_list_idx = list(range(min_idx_y, max_idx_y + 1))
        y_parts = []
        for i in y_list_idx:
            val = y_n.get(i, 0)
            y_parts.append(f"*{val}*" if i == 0 else str(val))
        y_n_str = f"{{{', '.join(y_parts)}}}"
    y_idx_str = ", ".join(str(i) for i in y_list_idx)

    y_values_list = [f"    y[{n}] = {y_n.get(n, 0)}" for n in range(y_start_idx, y_end_idx + 1)]
    y_values_str = "\n".join(y_values_list)

    calculation_steps_str = "\n".join(calculation_steps)
        
    solution = (
        f"**Given:**\n"
        f"Input Signal: x[n] = {x_n_str} for n = {x_idx_str}.\n"
        f"Impulse Response: h[n] = {h_n_str} for n = {h_idx_str}.\n\n"

        f"**Step 1:** State the Convolution Formula\n"
        f"The output y[n] of an LTI system is the convolution of the input x[n] with the impulse response h[n]. The convolution sum is defined as:\n"
        f"y[n] = sum over all k of (x[k] * h[n-k])\n\n"

        f"**Step 2:** Apply the Flip-and-Slide Method\n"
        f"We can visualize this process by flipping the impulse response h[k] to get h[-k], and then sliding it by 'n' positions. For each slide 'n', we calculate the sum of the products of the overlapping samples.\n"
        f"\n\n"
        f"Let's calculate a few points explicitly:\n\n"
        f"{calculation_steps_str}\n"
        
        f"**Step 3:** Calculate All Output Values\n"
        f"By continuing this process for all values of 'n' where the sequences overlap (from n={y_start_idx} to n={y_end_idx}), we get the full output sequence:\n"
        f"{y_values_str}\n"
        f"Collected in order of increasing n, starting at n = {y_list_idx[0]}, with the value at n = 0 also marked as in the question, this gives the answer.\n\n"

        f"**Answer:**\n"
        f"The complete output sequence is y[n] = {y_n_str} for n = {y_idx_str}"
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
        C_str = signed_term(C)   # C is never zero: randint(1,5) * choice([-1,1])
        equation_str = f"x[n] {C_str}"
        
        y1_str = f"x1[n] {C_str}"
        y2_str = f"x2[n] {C_str}"
        
        # Additivity
        add_sum_str = f"(x1[n] {C_str}) + (x2[n] {C_str}) = x1[n] + x2[n] {signed_term(2*C)}"
        add_y3_str = f"(x3[n]) {C_str} = (x1[n] + x2[n]) {C_str}"
        add_comparison_str = f"Since x1[n] + x2[n] {signed_term(2*C)} is not equal to x1[n] + x2[n] {C_str}, the system fails the additivity test."

        # Homogeneity
        hom_ay1_str = f"a * (x1[n] {C_str}) = a*x1[n] {signed_term(C)}*a"
        hom_ya_str = f"(xa[n]) {C_str} = (a * x1[n]) {C_str}"
        hom_comparison_str = f"Since a*x1[n] {signed_term(C)}*a is not equal to a*x1[n] {C_str} (for a != 1), the system fails the homogeneity test."

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

    Screen pass 1 (2026-09-23):
        Two judges reported the answer's formatting and one the stability
        wording. The second-order answer printed C2 unrounded
        ("0.06066017177982138(0.59)^n") because the local sign helper ignored
        its rounded argument (242/246 second-order seeds in 500), and the 2x2
        system was stated with roots at 2 dp but solved with the exact roots.
        The helper is replaced by the shared signed_term; each root is bound
        through its 2-dp display before the solve (D-016 part 2), so solving
        the printed system gives the printed C1, C2, which are rounded
        half-up at 2 dp (D-012); a draw whose C1 or C2 sits on a 2-dp tie is
        redrawn (D-016 part 3; 0/246 second-order seeds in 500). Products are
        written with an explicit "*" ("12*(-5)^(n-1)*u[n-1]"). The
        first-order text called (-a1)^(n-1) "a decaying exponential" on every
        seed although |a1| >= 1 always; the wording now follows |base|. The
        judges' "undefined signed_term" is the module import they were not
        shown. Question text is unchanged on every seed that is not redrawn;
        the gold coefficients move in the last digit where the bound roots
        round C1 or C2 differently (31/246 second-order seeds).

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

    # A signed coefficient renders as an operator and a magnitude ("+ 3", "- 2")
    # through the shared signed_term helper; _2f formats a 2-dp magnitude.
    _2f = lambda v: f"{v:.2f}"

    if order == 'first':
        # For first order, we can have b0*x[n] and b1*x[n-1] terms
        b1 = random.randint(1, 4) * random.choice([-1, 1])

        # Build equation strings
        y_terms = f"y[n] {signed_term(a1)}y[n-1]"
        # Randomly omit b1 term for variety
        if random.random() < 0.4:
            b1 = 0
        x_terms = f"{b0}x[n]"
        if b1 != 0:
            x_terms += f" {signed_term(b1)}x[n-1]"

        equation_str = f"{y_terms} = {x_terms}"

        # 2. Perform the core calculation for a first-order system
        h0 = b0
        h1 = b1 - a1 * h0
        base = -a1

        term1_str = f"{h0}*delta[n]"

        h_decay_expr = ""
        if h1 != 0:
            h_decay_expr = f" {signed_term(h1)}*({base})^(n-1)*u[n-1]"

        final_h_n = f"{term1_str}{h_decay_expr}"

        # (base)^(n-1) decays only for |base| < 1; with |a1| >= 1 it never
        # does, so the description follows the magnitude of the base.
        if abs(base) < 1:
            behaviour = "a decaying exponential"
        elif abs(base) == 1:
            behaviour = f"an exponential of constant magnitude (|{base}| = 1)"
        else:
            behaviour = f"a growing exponential (|{base}| > 1, so the system is unstable)"

        # 3. Generate the solution string for a first-order system
        solution_steps = (
            f"**Step 3:** Solve Recursively for Initial Conditions\n"
            f"We assume the system is causal, so **h[n] = 0 for n < 0**.\n\n"
            f"**For n = 0:**\n"
            f"h[0] {signed_term(a1)}h[-1] = {b0}*delta[0] {'+ 0' if b1==0 else signed_term(b1)+'*delta[-1]'}\n"
            f"h[0] {signed_term(a1)}*(0) = {b0}*(1) {signed_term(b1)}*(0)\n"
            f"**h[0] = {h0}**\n\n"

            f"**For n = 1:**\n"
            f"h[1] {signed_term(a1)}h[0] = {b0}*delta[1] {signed_term(b1)}*delta[0]\n"
            f"h[1] {signed_term(a1)}*({h0}) = {b0}*(0) {signed_term(b1)}*(1)\n"
            f"h[1] = {b1} - ({a1*h0})\n"
            f"**h[1] = {h1}**\n\n"

            f"**Step 4:** Find the Homogeneous Solution\n"
            f"For n >= 2, the input delta terms are zero. The equation becomes homogeneous:\n"
            f"h[n] {signed_term(a1)}h[n-1] = 0  =>  h[n] = {base}*h[n-1]\n\n"
            f"The solution to this recurrence for n >= 1 is {behaviour} that starts at n=1 with value h[1].\n"
            f"This part of the response can be written as h[1]*({base})^(n-1)*u[n-1].\n\n"

            f"**Step 5:** Combine Results for the Final Expression\n"
            f"The total impulse response is the sum of the value at n=0 and the response for n >= 1:\n"
            f"h[n] = h[0]*delta[n] + h[1]*({base})^(n-1)*u[n-1]\n"
        )

    else: # order == 'second'
        # Draw, and redraw if a coefficient sits on a 2-dp display tie (D-016
        # part 3). The first pass consumes the random stream exactly as it did
        # before, so a seed that never ties is unchanged.
        for _attempt in range(200):
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

            # 2. Perform the core calculation for a second-order system
            h0 = b0
            h1 = b1 - a1 * h0

            # Correctly solve r^2 + a1*r + a2 = 0. The roots are stated to 2 dp
            # and then consumed by the 2x2 solve for C1 and C2, so each is bound
            # through that display (D-016 part 2): solving the printed system
            # reproduces the printed coefficients.
            discriminant = a1**2 - 4*a2
            r1 = _as_printed((-a1 + math.sqrt(discriminant)) / 2, '.2f')
            r2 = _as_printed((-a1 - math.sqrt(discriminant)) / 2, '.2f')

            # Correctly solve for C1 and C2 with general h[0], h[1]
            C1 = (h1 - h0 * r2) / (r1 - r2)
            C2 = (h0 * r1 - h1) / (r1 - r2)
            if not (_is_display_tie(C1, 2) or _is_display_tie(C2, 2)):
                break

        r1_str, r2_str = _2f(r1), _2f(r2)
        # Decimal half-up at 2 dp (D-012); a value that rounds to zero prints as 0.00, not -0.00
        C1_d = _hu(C1, 2) or 0.0
        C2_d = _hu(C2, 2) or 0.0
        C1_str, C2_str = _2f(C1_d), _2f(C2_d)

        final_h_n = f"({C1_str}*({r1_str})^n {signed_term(C2_d, fmt=_2f)}*({r2_str})^n) * u[n]"

        # Build equation strings
        x_terms = f"{b0}x[n]"
        if b1 != 0: x_terms += f" {signed_term(b1)}x[n-1]"
        # b2 term removed
        equation_str = f"y[n] {signed_term(a1)}y[n-1] {signed_term(a2)}y[n-2] = {x_terms}"

        # 3. Generate the solution string for a second-order system
        h_eq = f"h[n] {signed_term(a1)}h[n-1] {signed_term(a2)}h[n-2]"
        d_eq = f"{b0}*delta[n]"
        if b1 != 0: d_eq += f" {signed_term(b1)}*delta[n-1]"

        solution_steps = (
            f"**Step 3:** Solve Recursively for Initial Conditions\n"
            f"We assume the system is causal, so **h[n] = 0 for n < 0**.\n\n"
            f"**For n = 0:**\n"
            f"h[0] {signed_term(a1)}h[-1] {signed_term(a2)}h[-2] = {b0}*delta[0] ... (other delta terms are 0)\n"
            f"h[0] {signed_term(a1)}*(0) {signed_term(a2)}*(0) = {b0}*(1)\n"
            f"**h[0] = {h0}**\n\n"
            f"**For n = 1:**\n"
            f"h[1] {signed_term(a1)}h[0] {signed_term(a2)}h[-1] = ... {signed_term(b1)}*delta[0] ...\n"
            f"h[1] {signed_term(a1)}*({h0}) {signed_term(a2)}*(0) = {b0}*(0) {signed_term(b1)}*(1)\n"
            f"h[1] = {b1} {signed_term(-(a1*h0))}\n"
            f"**h[1] = {h1}**\n\n"

            f"**Step 4:** Find the Homogeneous Solution\n"
            f"For n >= 2, the input is zero (since b2=0), so the equation becomes homogeneous:\n"
            f"h[n] {signed_term(a1)}h[n-1] {signed_term(a2)}h[n-2] = 0\n\n"
            f"We solve this by finding the roots of the characteristic equation: r^2 {signed_term(a1)}r {signed_term(a2)} = 0\n"
            f"Using the quadratic formula, the roots are r1 = {r1_str}, r2 = {r2_str}.\n"
            f"The general solution for n >= 0 is h[n] = C1*({r1_str})^n + C2*({r2_str})^n.\n\n"

            f"**Step 5:** Use Initial Conditions to Find Coefficients\n"
            f"We use h[0] and h[1] to create a system of two equations:\n"
            f"1) For n=0: h[0] = C1 + C2  =>  {h0} = C1 + C2\n"
            f"2) For n=1: h[1] = C1*({r1_str}) + C2*({r2_str})  =>  {h1} = {r1_str}*C1 {signed_term(r2, fmt=_2f)}*C2\n\n"
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
