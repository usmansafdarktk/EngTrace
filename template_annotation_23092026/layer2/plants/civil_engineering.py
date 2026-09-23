"""Planted defects for the civil engineering branch (Layer 2; see CONTRACT.md).

plant_civ_1  constant    manning_rectangular_discharge: Manning's R^(2/3) replaced by R^(3/4) throughout
plant_civ_2  unit        cantilever_double_integration: metres converted to millimetres by x100 instead of x1000
plant_civ_3  sign        effective_stress_profile: pore pressure added to the total stress (sigma' = sigma + u)
plant_civ_4  arithmetic  beam_support_reactions: Step 2 prints By 4% high; Step 3 and the answer use the true By
"""

PLANTS = [
    {
        'plant_id': 'plant_civ_1',
        'base': 'template_manning_rectangular_discharge',
        'defect_class': 'constant',
        'description': "Manning's equation applied with R^(3/4) instead of R^(2/3), consistently "
                       "in the docstring, the computation, Step 3 and Step 4 (variable renamed R34 "
                       "to match); the discharge is off by the factor R^(1/12).",
        'detectable_by': "Step 3 raises the hydraulic radius to the 3/4 power and Step 4 writes "
                         "Q = (1/n) * A * R^(3/4) * S^(1/2); Manning's exponent on R is 2/3, and a "
                         "hand solution with 2/3 differs from the printed discharge by a few per cent.",
        'edits': [
            ('R^(2/3) * S^(1/2)\n', 'R^(3/4) * S^(1/2)\n'),
            ('R23 = round(R ** (2.0 / 3.0), 4)', 'R34 = round(R ** (3.0 / 4.0), 4)'),
            ('Q = round((1.0 / n) * A * R23 * sqS, 2)', 'Q = round((1.0 / n) * A * R34 * sqS, 2)'),
            ('f"R^(2/3) = ({R:.4f})^(2/3) = {R23:.4f}', 'f"R^(3/4) = ({R:.4f})^(3/4) = {R34:.4f}'),
            ('f"Q = (1/n) * A * R^(2/3) * S^(1/2) "', 'f"Q = (1/n) * A * R^(3/4) * S^(1/2) "'),
            ('* {A:.3f} * {R23:.4f} * {sqS:.5f}', '* {A:.3f} * {R34:.4f} * {sqS:.5f}'),
        ],
        'reskin': [
            ('width b in [2.0, 6.0] m', 'width b in [2.5, 6.0] m'),
            ('b = round(random.uniform(2.0, 6.0), 1)', 'b = round(random.uniform(2.5, 6.0), 1)'),
            ('f"A rectangular open channel is {b:.1f} m wide and has "',
             'f"A rectangular irrigation canal is {b:.1f} m wide and has "'),
        ],
    },
    {
        'plant_id': 'plant_civ_2',
        'base': 'template_cantilever_double_integration',
        'defect_class': 'unit',
        'description': "The tip deflection in metres is converted to millimetres by multiplying by "
                       "100 instead of 1000, consistently (conversion, its exactness assert, the tie "
                       "guard and the deflection band all scaled with it; Step 5 shows metres at 3 dp "
                       "so it still closes); the answer is ten times too small.",
        'detectable_by': "Step 6 reads 'delta = 0.0xx * 100 = x.x mm': a metre-to-millimetre "
                         "conversion by 100. A hand solution gives ten times the printed answer, and "
                         "a few millimetres is implausibly small for the stated load and length.",
        'edits': [
            ('recomputed delta in [3, 45] mm.', 'recomputed delta in [0.3, 4.5] mm.'),
            ('if not _is_display_tie(delta * 1000, 1):', 'if not _is_display_tie(delta * 100, 1):'),
            ('delta_mm = _hu(delta * 1000, 1)', 'delta_mm = _hu(delta * 100, 1)'),
            ('delta = _hu(delta_mm / 1000.0, 4)            # exact at 4 dp; what Step 5 states',
             'delta = _hu(delta_mm / 100.0, 3)            # exact at 3 dp; what Step 5 states'),
            ('{sub_expr} = {delta:.4f} m', '{sub_expr} = {delta:.3f} m'),
            ('assert 3.0 <= delta_mm <= 45.0, f"tip deflection out of bounds: {delta_mm}"',
             'assert 0.3 <= delta_mm <= 4.5, f"tip deflection out of bounds: {delta_mm}"'),
            ('# The x1000 conversion Step 6 prints', '# The x100 conversion Step 6 prints'),
            ('assert Decimal(f"{delta:.4f}") * 1000 == Decimal(f"{delta_mm:.1f}")',
             'assert Decimal(f"{delta:.3f}") * 100 == Decimal(f"{delta_mm:.1f}")'),
            ('f"delta = {delta:.4f} * 1000 = {delta_mm:.1f} mm', 'f"delta = {delta:.3f} * 100 = {delta_mm:.1f} mm'),
        ],
        'reskin': [
            ('cantilever length from L/d in [8, 14]', 'cantilever length from L/d in [9, 14]'),
            ('ld_ratio = random.uniform(8.0, 14.0)', 'ld_ratio = random.uniform(9.0, 14.0)'),
            ('f"A cantilever steel beam is built from a {shape} section (SI "',
             'f"A cantilever steel beam supporting a service platform is built from a {shape} section (SI "'),
        ],
    },
    {
        'plant_id': 'plant_civ_3',
        'base': 'template_effective_stress_profile',
        'defect_class': 'sign',
        'description': "The pore water pressure is added to the total stress instead of subtracted: "
                       "sigma' = sigma + u in the docstring, the computation and Step 5; the answer "
                       "is too high by 2u. The reskin's narrower depths (z1 2-3 m, z2 3-4.5 m) are "
                       "load-bearing: they keep sigma + u inside the template's own [40, 220] kPa "
                       "assert, which the base depths would breach.",
        'detectable_by': "Step 5 states sigma' = sigma + u and adds the pore pressure to the total "
                         "stress; the effective stress principle is sigma' = sigma - u, so the printed "
                         "effective stress exceeds the total stress computed in Step 3.",
        'edits': [
            ("            sigma' = sigma - u\n", "            sigma' = sigma + u\n"),
            ('sigma - u is the exact difference of those two 3-dp values',
             'sigma + u is the exact sum of those two 3-dp values'),
            ('# The difference of two 3-dp values is exact at 3 dp; the answer',
             '# The sum of two 3-dp values is exact at 3 dp; the answer'),
            ('sigma_eff = _hu(sigma - u, 3)', 'sigma_eff = _hu(sigma + u, 3)'),
            ('f"sigma\' = sigma - u = {sigma:.3f} - {u:.3f} = {sigma_eff:.3f} kPa',
             'f"sigma\' = sigma + u = {sigma:.3f} + {u:.3f} = {sigma_eff:.3f} kPa'),
        ],
        'reskin': [
            ('z1 in [2.0, 5.0] m, z2 in [3.0, 8.0] m', 'z1 in [2.0, 3.0] m, z2 in [3.0, 4.5] m'),
            ('z1 = round(random.uniform(2.0, 5.0), 1)', 'z1 = round(random.uniform(2.0, 3.0), 1)'),
            ('z2 = round(random.uniform(3.0, 8.0), 1)', 'z2 = round(random.uniform(3.0, 4.5), 1)'),
            ('f"At a site, a deep deposit of sand has its groundwater table at a "',
             'f"At a proposed warehouse site, a deep deposit of sand has its groundwater table at a "'),
        ],
    },
    {
        'plant_id': 'plant_civ_4',
        'base': 'template_beam_support_reactions',
        'defect_class': 'arithmetic',
        'description': "Step 2 prints the right reaction By 4% higher than (P*a + W*x_bar)/L "
                       "evaluates to; Step 3 and the answer use the true By, so only that one "
                       "line is wrong.",
        'detectable_by': "Recomputing Step 2 from its printed operands gives a By about 4% below the "
                         "printed value, and Step 3 then subtracts a By that differs from the one "
                         "Step 2 printed.",
        'edits': [
            ('= {By:.2f} kN', '= {By * 1.04:.2f} kN'),
        ],
        'reskin': [
            ('point load P in [10, 50] kN', 'point load P in [15, 50] kN'),
            ('P = random.randint(10, 50)', 'P = random.randint(15, 50)'),
            ('f"A simply supported beam AB has a span of L = {L:.1f} m, with a "',
             'f"A simply supported steel floor beam AB has a span of L = {L:.1f} m, with a "'),
        ],
    },
]
