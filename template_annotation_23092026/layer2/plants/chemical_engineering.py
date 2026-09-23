"""Planted defects for the chemical engineering branch (Layer 2 expert certification).

plant_che_1  constant    template_vdw_solve_for_pressure: R = 0.08134 L·bar/(mol·K) used throughout instead of 0.08314.
plant_che_2  unit        template_annulus_flowrate: the pressure drop is converted between kPa and Pa by 100 instead of 1000.
plant_che_3  sign        template_heat_of_reaction_formation: Hess's law reversed (reactants minus products), so ΔH_rxn° has the wrong sign.
plant_che_4  arithmetic  template_batch_reactor_first_order: the last Step 6 line prints t 4% above the product of its printed factors.

Shape and rules: CONTRACT.md in this directory. Each plant is a copy of a real
template of this branch with one defect written into its source; the real
template files are never touched.
"""

PLANTS = [
    {
        'plant_id': 'plant_che_1',
        'base': 'template_vdw_solve_for_pressure',
        'defect_class': 'constant',
        'description': ("gas constant R = 0.08134 L·bar/(mol·K) (digits transposed) used for b, a, "
                        "R*T and P instead of 0.08314"),
        'detectable_by': ("Step 2 lists the gas constant as R = 0.08134 L·bar/(mol·K); the correct value "
                          "is 0.08314, and b (Step 3), a (Step 4), R*T and P (Steps 5-6) all follow "
                          "from the wrong one"),
        'edits': [
            ('R = 0.08314', 'R = 0.08134'),
        ],
        'reskin': [
            ('f"Using the van der Waals equation of state, calculate the pressure in bar "',
             'f"A storage cylinder is charged with {substance_name}. Using the van der Waals '
             'equation of state, calculate the pressure in bar "'),
            ('random.uniform(0.8 * Tc, 2.5 * Tc)', 'random.uniform(0.9 * Tc, 2.3 * Tc)'),
        ],
    },
    {
        'plant_id': 'plant_che_2',
        'base': 'template_annulus_flowrate',
        'defect_class': 'unit',
        'description': ("pressure drop converted between kPa and Pa with a factor of 100 instead of "
                        "1000: the question states Pa/100 as the kPa value and Step 2 multiplies it "
                        "back by 100, so the stated kPa is ten times the Pa value actually used"),
        'detectable_by': ("Step 2 reads 'P0 - PL = X kPa * 100 = Y Pa'; 1 kPa is 1000 Pa, so the Pa "
                          "value carried into the main term and Q is one tenth of the stated pressure "
                          "drop"),
        'edits': [
            ('pressure_drop_kPa = pressure_drop_Pa / 1000.0', 'pressure_drop_kPa = pressure_drop_Pa / 100.0'),
            ('kPa * 1000 = {pressure_drop_Pa} Pa', 'kPa * 100 = {pressure_drop_Pa} Pa'),
        ],
        'reskin': [
            ('f"Laminar flow of {fluid_name} occurs in the annular space between two concentric pipes. "',
             'f"{fluid_name} flows in laminar fashion through the annular gap of a double-pipe heat exchanger. "'),
            ('pipe_length_m = round(random.uniform(10.0, 50.0), 1)',
             'pipe_length_m = round(random.uniform(12.0, 40.0), 1)'),
        ],
    },
    {
        'plant_id': 'plant_che_3',
        'base': 'template_heat_of_reaction_formation',
        'defect_class': 'sign',
        'description': ("Hess's law reversed: ΔH_rxn° stated and computed as Σ(reactants) - Σ(products), "
                        "so every heat of reaction has the opposite sign (exothermic reactions come "
                        "out endothermic)"),
        'detectable_by': ("Step 1 states ΔH_rxn° = Σ(v * ΔH_f°)reactants - Σ(v * ΔH_f°)products and "
                          "Step 4 subtracts in that order; the correct form is products minus "
                          "reactants, so the answer's sign is wrong (a combustion is reported as "
                          "absorbing heat)"),
        'edits': [
            (r'f"ΔH_rxn° = Σ(v * ΔH_f°)products - Σ(v * ΔH_f°)reactants\n\n"',
             r'f"ΔH_rxn° = Σ(v * ΔH_f°)reactants - Σ(v * ΔH_f°)products\n\n"'),
            (r'f"ΔH_rxn° = ({products_enthalpy:.3f}) - ({reactants_enthalpy:.3f})\n"',
             r'f"ΔH_rxn° = ({reactants_enthalpy:.3f}) - ({products_enthalpy:.3f})\n"'),
            ('delta_H_rxn = products_enthalpy - reactants_enthalpy',
             'delta_H_rxn = reactants_enthalpy - products_enthalpy'),
        ],
        'reskin': [
            (r'f"Calculate the standard enthalpy of reaction, ΔH_rxn°, at 298.15 K (in kJ) for the reaction as written:\n"',
             r'f"A process engineer needs the standard enthalpy of reaction, ΔH_rxn°, at 298.15 K (in kJ) for the reaction as written below:\n"'),
            # heats of formation listed in reaction order (reactants, then products) instead of alphabetically
            ('sorted(list(set(required_species)))', 'list(dict.fromkeys(required_species))'),
        ],
    },
    {
        'plant_id': 'plant_che_4',
        'base': 'template_batch_reactor_first_order',
        'defect_class': 'arithmetic',
        'description': ("the last line of Step 6 prints t multiplied by 1.04; the two printed factors "
                        "on the line above and the Answer line carry the true value"),
        'detectable_by': ("in Step 6, 't = (1/k) × ln(C_A0/C_A)' is followed by a result 4% larger than "
                          "the product of the two printed factors, and the Answer line then quotes a "
                          "different (the correct) time"),
        'edits': [
            (r'f"t = {round(time, 2)} s\n\n"', r'f"t = {round(time * 1.04, 2)} s\n\n"'),
        ],
        'reskin': [
            ('f"A first-order liquid-phase reaction of {reactant_name} occurs in a batch reactor. "',
             'f"A first-order liquid-phase decomposition of {reactant_name} is carried out in a '
             'stirred, jacketed batch vessel. "'),
            ('C_A0 = round(random.uniform(1.0, 5.0), 2)', 'C_A0 = round(random.uniform(1.5, 4.5), 2)'),
        ],
    },
]
