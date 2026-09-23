# Planted defects for Layer 2 - the contract every plants/<branch>.py must meet

A plant is a real template of the same branch with ONE deliberate defect written into a copy
of its source. Experts see plants mixed into their queue, unmarked, and their verdict on each
plant is the measurable part of certification: an expert who approves a planted defect is
measurably not reviewing. Plants never touch the real template files.

Each `plants/<branch>.py` defines `PLANTS`, a list of dicts:

    {
      'plant_id':      'plant_civ_1',                 # unique; <branch-prefix>_<n>
      'base':          'template_beam_internal_moment', # a template of THIS branch
      'defect_class':  'constant' | 'unit' | 'sign' | 'arithmetic',
      'description':   'g = 9.18 used instead of 9.81 throughout',   # keyfile only, never shipped
      'detectable_by': 'the Given line states g = 9.18 m/s^2',        # what a careful expert sees
      'edits':         [('exact old substring', 'new substring'), ...],
      'reskin':        [('exact old substring', 'new substring'), ...],   # optional wording changes
    }

Rules the build enforces (build_tasks.py refuses a plant that breaks one):
1. Every `old` substring in `edits` and `reskin` occurs EXACTLY ONCE in the base function's
   source (inspect.getsource of the function), so the mutation is unambiguous.
2. The mutated function is executed in the base module's namespace (imports, helpers and
   constants available), with `random.seed(s)` for s in 0..49, and must run without raising.
3. Its output must differ from the base's on every seed (the defect is live), and must be
   deterministic (same seed, same text, in two processes).
4. Exactly one defect class per plant. The trace must stay internally consistent EXCEPT for
   the planted defect: a 'constant' plant uses the wrong constant everywhere it appears, a
   'unit' plant applies the wrong conversion consistently, a 'sign' plant flips one sign or
   direction consistently, and an 'arithmetic' plant changes ONE printed step result so that
   line no longer follows from its printed operands (the December Tribunal's kind of catch).
5. The defect must be detectable from what the expert is shown: the rendered instances and
   the source. It must not be hidden in a helper the expert cannot see.
6. `reskin` may change scenario wording, object or material names, or shift a sampling range
   within the template's own validity, so the plant does not read as a duplicate of its base.
   It must not change the physics or the step structure.

Four plants per branch, one of each defect class, from four different templates and, where
the branch allows, different sub-areas. Verify with `python -m template_annotation_23092026.layer2.build_tasks --check-plants`.
