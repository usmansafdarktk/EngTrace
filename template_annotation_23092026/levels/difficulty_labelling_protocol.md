# Difficulty Labelling Protocol

We used the following protocol for difficulty rating of the templates.

## 1. Rating scale

Rate each template on the paper's three dimensions, on a 1–3 scale. Rate only what the solver must supply: a formula or method stated in the question does not count.

| Dimension | 1 | 2 | 3 |
|---|---|---|---|
| **Conceptual complexity** | One principle, applied directly | Several principles from one domain, or one modelling choice (regime, configuration, control volume) | Principles from several domains, or assumptions the solver must justify |
| **Mathematical sophistication** | Direct substitution into formulas | Rearrangement, simultaneous or nonlinear equations, probability beyond one formula | Calculus, differential equations, iterative or numerical solution |
| **Procedural depth** | 1–2 dependent intermediate quantities | 3–5 | 6 or more |

The level is set by the sum of the three scores, using this rule, fixed before rating:

| Sum | Level |
|---|---|
| 3–4 | Easy |
| 5–6 | Intermediate |
| 7–9 | Advanced |

Templates with several reasoning paths are rated per path. The template takes the median level of its paths.

## 2. Procedure

1. **Author.** The author records scores on all three dimensions, with a one-line rationale.
2. **Independent raters.** Two experts from the branch who did not write the template rate it independently. They see the template code and one instance per path. They do not see the author's scores or any model results. Raters calibrate first on a shared anchor set.
3. **LLM panel check.** Run the difficulty level validation prompt on every template with a panel of three LLM judges from outside the evaluated model families. Use the LLM-screen panel (Grok 4.6, MiniMax M3, MiMo-V2.5-Pro) at temperature 0. The panel's level is the median of the three judges' levels, which equals the majority vote whenever two judges agree. The panel only flags templates for review; it never sets a label.
4. **Adjudication.** Each dimension's final score is the median of the three expert scores. A template goes back to its three experts for discussion if their scores differ by 2 points on a dimension, or if the panel's level differs from the experts'.
5. **Agreement.** Report Krippendorff's α (ordinal) and Gwet's AC2, per dimension and for the final level. Target α ≥ 0.80, with 0.67 as the minimum. Report the same statistics for the panel, and the panel's agreement with the experts.
6. **Validation and freeze.** Check the procedural depth ratings against the milestone count computed from each gold trace; a large mismatch sends the template back to its experts. Check that hand-solve time recorded during certification rises with level. Hash the label file before any model runs, and never relabel based on model results.

## 3. Difficulty level validation prompt

```text
You are an expert academic and senior engineering educator at a top-tier engineering university (like MIT or Caltech). You have decades of experience designing undergraduate courses and examinations for ABET-accredited programs.

Your task is to rate the difficulty of an engineering problem template for a student who has just completed the core undergraduate course in its domain. Rate only what the student must supply: if the question states a formula, method, or principle, do not count it.

Rate the template on three dimensions, each on a scale of 1 to 3:

Conceptual complexity:
1 = one principle, applied directly
2 = several principles from one domain, or one modelling choice (e.g., flow regime, configuration, control volume)
3 = principles from several domains, or assumptions the student must choose and justify

Mathematical sophistication:
1 = direct substitution into formulas (arithmetic, powers, roots, logarithms)
2 = rearrangement, simultaneous or nonlinear equations, or probability beyond a single formula
3 = calculus, differential equations, or an iterative or numerical solution

Procedural depth:
1 = 1-2 dependent intermediate quantities
2 = 3-5 dependent intermediate quantities
3 = 6 or more dependent intermediate quantities

Sum the three scores and assign a level: 3-4 = Easy, 5-6 = Intermediate, 7-9 = Advanced.

If the instances follow different reasoning paths, rate each path separately.

Please evaluate the following template:
Domain: [DOMAIN_NAME_HERE]
Area: [AREA_NAME_HERE]
Template source code: [TEMPLATE_CODE_HERE]
Instances (one per reasoning path): [INSTANCES_HERE]

Return your answer only in a strict JSON format as given below.

{
  "template": "[TEMPLATE_NAME_HERE]",
  "paths": [
    {
      "path": "<short description of the reasoning path>",
      "conceptual_complexity": <integer 1-3>,
      "mathematical_sophistication": <integer 1-3>,
      "procedural_depth": <integer 1-3>,
      "level": "<Easy | Intermediate | Advanced>",
      "justification": "<one sentence explaining the scores>"
    }
  ]
}

Do not include any preamble, conversational text, explanations, or markdown formatting around the JSON block.
```
