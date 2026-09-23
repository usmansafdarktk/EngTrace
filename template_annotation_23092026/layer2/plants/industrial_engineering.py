"""Planted defects for the industrial engineering branch (Layer 2 certification).

plant_ind_1  constant    safety_stock_reorder_point: z is the two-sided quantile Phi^-1((1+alpha)/2) for a one-sided cycle-service level, so every level's z is one step too high (95% -> 1.9600 instead of 1.6449).
plant_ind_2  unit        mm1_time_in_system: hours are converted to minutes with a factor of 100 instead of 60, in the computation, the tie screen and the printed Step 4.
plant_ind_3  sign        xbar_r_control_limits: LCL_x is written as x-double-bar + A2*R-bar (centre plus instead of minus), so the printed LCL_x coincides with UCL_x.
plant_ind_4  arithmetic  cp_cpk_from_specs: the printed Cp of Step 2 is 6% above the quotient of its printed operands; sigma-hat, Cpu, Cpl, Cpk and the answer are untouched.

Shape and rules: CONTRACT.md in this directory. Every `old` substring below occurs exactly once
in inspect.getsource() of the base function; the four bases are four different templates from
three sub-areas (production and inventory, stochastic operations, quality and reliability, the
last with one variables-chart and one process-capability template).
"""

PLANTS = [
    {
        'plant_id': 'plant_ind_1',
        'base': 'template_safety_stock_reorder_point',
        'defect_class': 'constant',
        'description': (
            'z taken as the two-sided quantile Phi^-1((1 + alpha)/2) instead of the one-sided '
            'Phi^-1(alpha) for the stated cycle-service level, used throughout: 90% -> 1.6449, '
            '95% -> 1.9600, 98% -> 2.3263, 99% -> 2.5758 (each one level too high)'
        ),
        'detectable_by': (
            'the question and the Given line state a z one level too high for the named '
            'service level (z = 1.9600 for 95%, where a Type-1 level needs Phi^-1(0.95) = '
            '1.6449); SS and R then follow consistently from the wrong z, and the source reads '
            'Z_QUANTILES[round((1 + alpha) / 2, 3)]'
        ),
        'edits': [
            ('z = Z_QUANTILES[alpha]',
             'z = Z_QUANTILES[round((1 + alpha) / 2, 3)]'),
        ],
        'reskin': [
            ('mu = random.randint(50, 2000)',
             'mu = random.randint(100, 1600)'),
            ('A regional warehouse reorders a product whenever its inventory ',
             'An automotive parts distribution centre reorders a fast-moving SKU '
             'whenever its inventory '),
        ],
    },
    {
        'plant_id': 'plant_ind_2',
        'base': 'template_mm1_time_in_system',
        'defect_class': 'unit',
        'description': (
            'hours converted to minutes by a factor of 100 instead of 60, consistently: in '
            'W_min, in the display-tie screen on the minute answer, and in the printed Step 4'
        ),
        'detectable_by': (
            "Step 4 reads 'W = x.xxxx hours * 100 minutes/hour = ... minutes'; the minute "
            'answer is 100 times the hours value (e.g. 0.2500 hours -> 25.0 minutes instead '
            'of 15.0), while rho, L and W in hours in Steps 1-3 are consistent'
        ),
        'edits': [
            ('W_min = round(W_hr * 60, 1)',
             'W_min = round(W_hr * 100, 1)'),
            ('_is_display_tie(W_hr * 60, 1)',
             '_is_display_tie(W_hr * 100, 1)'),
            ('hours * 60 minutes/hour',
             'hours * 100 minutes/hour'),
        ],
        'reskin': [
            # utilization window [0.55, 0.80], inside the base's [0.55, 0.92]
            ('math.floor(0.92 * mu)',
             'math.floor(0.80 * mu)'),
            ('Customers arrive at {setting} according to a Poisson process at ',
             'During the afternoon peak, customers arrive at {setting} according to '
             'a Poisson process at '),
        ],
    },
    {
        'plant_id': 'plant_ind_3',
        'base': 'template_xbar_r_control_limits',
        'defect_class': 'sign',
        'description': (
            'LCL_x written as x-double-bar + A2*R-bar (centre plus instead of minus), in the '
            'docstring formula, the Decimal chain and the printed Step 2, so LCL_x = UCL_x'
        ),
        'detectable_by': (
            'Step 2 prints LCL_x = x-double-bar + A2 * R-bar and its value coincides with '
            'UCL_x; the lower limit must be x-double-bar - A2*R-bar (the R-chart limits and '
            'sigma-hat are unaffected)'
        ),
        'edits': [
            ('LCL_x = xbb - A2 * Rbar',
             'LCL_x = xbb + A2 * Rbar'),
            ('dx - dA2 * dr',
             'dx + dA2 * dr'),
            ('LCL_x = {xbb:.{dp}f} - {A2} *',
             'LCL_x = {xbb:.{dp}f} + {A2} *'),
        ],
        'reskin': [
            ('A quality engineer is setting up X-bar and R charts for ',
             'A process engineer requalifying a production line is setting up X-bar and '
             'R charts for '),
            # subgroup count [24, 30], inside the base's [20, 30]
            ('m = random.randint(20, 30)',
             'm = random.randint(24, 30)'),
        ],
    },
    {
        'plant_id': 'plant_ind_4',
        'base': 'template_cp_cpk_from_specs',
        'defect_class': 'arithmetic',
        'description': (
            'the printed Step 2 result for Cp is cp * 1.06 (one f-string changed); sigma-hat, '
            'Cpu, Cpl, Cpk, the verdict and the answer use the true values'
        ),
        'detectable_by': (
            'recomputing Step 2 from its printed operands, (USL - LSL) / (6 * sigma-hat), '
            'gives a value about 6% below the printed Cp; Step 3 closes on the same operands, '
            'and the printed Cp no longer equals (Cpu + Cpl)/2 as it must for two-sided limits'
        ),
        'edits': [
            ('= {cp:.2f}',
             '= {cp * 1.06:.2f}'),
        ],
        'reskin': [
            ("An X-bar/R chart study of {cfg['phrase']} shows the process ",
             "A six-week X-bar/R chart study of {cfg['phrase']} shows the process "),
            # target one-sided ratios [0.85, 1.95], inside the base's [0.75, 2.05]
            ('cpl_t = random.uniform(0.75, 2.05)',
             'cpl_t = random.uniform(0.85, 1.95)'),
            ('cpu_t = random.uniform(0.75, 2.05)',
             'cpu_t = random.uniform(0.85, 1.95)'),
        ],
    },
]
