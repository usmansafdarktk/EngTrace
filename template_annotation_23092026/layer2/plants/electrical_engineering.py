"""Planted defects for the electrical engineering branch (plants/CONTRACT.md).

plant_ele_1  constant    template_coulombs_law             epsilon_0 = 8.854e-11 F/m (ten times the true value) bound locally and used everywhere it appears
plant_ele_2  unit        template_wave_parameters_basic    MHz <-> Hz converted with 1e5 instead of 1e6 in both scenarios
plant_ele_3  sign        template_cd_dc_system_analysis    the D/C step writes the delay n0*T as an advance, t + n0 T, in every expression
plant_ele_4  arithmetic  template_mean_variance            the mu_X line prints a sum 0.2 above its printed terms; Step 2 and the answer keep the true mean

Every base is Easy or Intermediate and its solution carries numeric steps an expert can
recompute. Three of the four bases were edited on 2026-09-23 and bind intermediates through
`_as_printed` / `_hu` / `_hu_dec` and format with `signed_term`, `joined_terms`, `paren_neg`
and `rect_str` from `data/templates/branches/_emission.py`. Those helpers, `Decimal`, `np`
and the branch constants are imported at module level in each base file, so they are in the
namespace the mutated source is executed in (build_tasks.mutated_function copies the base
module's `__dict__`); the edits below leave every helper call untouched.

Old substrings that include an f-string's literal `\\n` are written as raw strings so they
match the two-character escape in the base source, not a real newline.
"""

PLANTS = [
    {
        'plant_id': 'plant_ele_1',
        'base': 'template_coulombs_law',
        'defect_class': 'constant',
        'description': ('epsilon_0 = 8.854e-11 F/m, ten times the true 8.854e-12, bound to a local '
                        'eps0 and used in the Constants block, the epsilon line and the F12 line; '
                        'every force component is ten times too small'),
        'detectable_by': ('the Constants block states epsilon_0 = 8.854e-11 F/m; a reader who knows '
                          'epsilon_0 = 8.854e-12 F/m, or checks 1/(4*pi*epsilon) in Step 3, sees the '
                          'tenfold error carried through to F12'),
        'edits': [
            ('epsilon = epsilon_r * EPSILON_0',
             'eps0 = 8.854e-11   # permittivity of free space, F/m\n    epsilon = epsilon_r * eps0'),
            ('epsilon_0 = {EPSILON_0:.4e} F/m',
             'epsilon_0 = {eps0:.4e} F/m'),
            ('{epsilon_r} * {EPSILON_0:.4e} = ',
             '{epsilon_r} * {eps0:.4e} = '),
        ],
        'reskin': [
            ('Two point charges are located in a medium with a relative permittivity of',
             'Two point charges are embedded in a dielectric block with a relative permittivity of'),
            ('epsilon_r = round(random.uniform(1.0, 5.0), 2)',
             'epsilon_r = round(random.uniform(1.5, 6.0), 2)'),
        ],
    },
    {
        'plant_id': 'plant_ele_2',
        'base': 'template_wave_parameters_basic',
        'defect_class': 'unit',
        'description': ('MHz converted to Hz with 1e5 instead of 1e6 in the frequency-first scenario, '
                        'and Hz to MHz with /1e5 in the wavelength-first scenario; omega, T, lambda and '
                        'k follow consistently from the wrong f'),
        'detectable_by': ('the Given line reads e.g. "f = 123.4 MHz = 1.23e+07 Hz", a factor of ten '
                          'low; in the wavelength scenario Step 1 quotes a MHz value ten times the '
                          'Hz value divided by 1e6'),
        'edits': [
            ("f = _as_printed(f_mhz * 1e6, '.2e')",
             "f = _as_printed(f_mhz * 1e5, '.2e')"),
            ("f_mhz_out = _as_printed(f_exact / 1e6, '.2f')",
             "f_mhz_out = _as_printed(f_exact / 1e5, '.2f')"),
        ],
        'reskin': [
            ('A sinusoidal wave with a frequency of {f_mhz} MHz is propagating',
             'A radio-frequency carrier at {f_mhz} MHz is propagating'),
            ('An electromagnetic wave traveling in {medium_name} is observed to have',
             'A plane electromagnetic wave traveling in {medium_name} is measured to have'),
            ('f_mhz = round(random.uniform(50, 500), 1)',
             'f_mhz = round(random.uniform(60, 600), 1)'),
        ],
    },
    {
        'plant_id': 'plant_ele_3',
        'base': 'template_cd_dc_system_analysis',
        'defect_class': 'sign',
        'description': ('the D/C conversion maps (n - n0)T to t + n0 T: the delay is written as an '
                        'advance in the symbolic Step 3 expression, the numeric one and the answer, '
                        'while Step 2 and the prose still describe a delay'),
        'detectable_by': ('Step 2 ends with cos(Omega*T*(n - n0)) and Step 3 turns it into '
                          'cos(Omega*(t + n0 T)) under t = nT; the text calls n0*T a time delay '
                          'while every expression adds it'),
        'edits': [
            ('cos({omega_continuous_str} * (t - {delay_n0}T))',
             'cos({omega_continuous_str} * (t + {delay_n0}T))'),
            (r'f"y_c(t) = {output_amplitude} * cos({omega_continuous_str}*(t - {time_delay_sec:.6f}))\n\n"',
             r'f"y_c(t) = {output_amplitude} * cos({omega_continuous_str}*(t + {time_delay_sec:.6f}))\n\n"'),
            ('output signal is y_c(t) = {output_amplitude} * cos({omega_continuous_str}*(t - {time_delay_sec:.6f}))',
             'output signal is y_c(t) = {output_amplitude} * cos({omega_continuous_str}*(t + {time_delay_sec:.6f}))'),
        ],
        'reskin': [
            ('A continuous-time signal is given by:',
             'A transducer output is modelled by the continuous-time signal:'),
            ('sampling_factor = random.randint(3, 8)',
             'sampling_factor = random.randint(4, 9)'),
        ],
    },
    {
        'plant_id': 'plant_ele_4',
        'base': 'template_mean_variance',
        'defect_class': 'arithmetic',
        'description': ('the printed result of the mu_X line is the true mean plus 0.200; the '
                        'printed terms, the deviations in Step 2 and the answer all carry the true mean'),
        'detectable_by': ('adding the printed terms of the mu_X line gives a sum 0.200 below the '
                          'stated result, and Step 2 then says "Using the calculated mean" with a '
                          'different value from the one Step 1 stated'),
        'edits': [
            (r'= {mean_terms} = {mean_d:.{precision}f}\n\n',
             r"= {mean_terms} = {mean_d + Decimal('0.2'):.{precision}f}\n\n"),
        ],
        'reskin': [
            ('A discrete random variable X can take the values',
             'The symbol level at a receiver decision device is a discrete random variable X that can take the values'),
            ('values = sorted(random.sample(range(-10, 21), num_values))',
             'values = sorted(random.sample(range(-12, 25), num_values))'),
        ],
    },
]
