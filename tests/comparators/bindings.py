"""D5.8 / D5.10 -- what each template's answer IS, declared per template.

Phase 4 delivered the comparator contract and bound **4 of 150** templates to
it. Everything else was scored with defaults, and a multipart answer read by a
single-part comparator returned whatever number came first -- which is why
``euclidean_distance_binary`` matched 132 of 132 cross-paired instances of
itself (D-058, N3). **N3 was larger than N1 and N2 together, and it is not a
missing rule: it is a missing declaration.**

Two things live here, and both are *problem data* in |S|3.8's sense -- stated once,
on the gold side, never restated by a candidate (D4.1 |S|2.1):

``DECLARED_UNITS``  D5.10. The unit gold's answer carries, where it carries one
    invariantly. Declared, never inferred: inferring a unit from gold's own
    string makes the comparison agree with whatever gold happens to say, and it
    would reject ``7.65 litres`` against a gold of ``7.65 L`` (D-052).
``BINDINGS``        D5.8. The comparator ``kind`` and its options.

**Both tables are generated, not typed** -- ``python -m tests.comparators.bindings
--regenerate`` rebuilds this file's data from the corpus and the validation run,
and ``--check`` fails if the committed table disagrees with what the corpus now
produces. A binding typed by hand beside the rule it describes is D-034's failure
mode, and the phase that wrote D-034 does not get to commit it.

## What "bound" means, and why it is not "declared"

> A binding counts as **bound** when it produces **zero false accepts over
> N >= 50 instances cross-paired within its own template**, and **zero errors**.

50 instances is 2,450 ordered pairs and resolves ~0.12% per template. Twelve
instances -- Phase 4's figure -- resolves nothing rarer than 25% (D-024/D-026),
and generation is deterministic and cheap, so the bar goes up rather than the
sample down. **A template that raises has not passed** (D5.7b): errors are
counted separately from verdicts, because N4 satisfied "zero false accepts" by
crashing on every pair.

*An unvalidated binding is a false accept waiting to happen; a named unbound
template with a written reason is a result.* ``UNBOUND`` carries the reason for
every template that did not clear the bar, and those reasons are measurements,
not judgements.

## The inventory's ``answer_type`` is NOT the comparator's ``kind``

It is used here only to report coverage. Three examples of why it cannot be used
as a binding:

* ``normal_depth_iteration`` is typed ``array`` and answers ``The normal depth
  is 1.380 m`` -- a scalar.
* ``two_state_steady_state`` is typed ``vector`` and answers a single fraction.
* ``critical_depth_froude_classification`` is typed ``classification`` and
  answers ``The Froude number is 2.361`` -- a number.

Binding from the inventory would have produced three wrong bindings that each
validate cleanly against themselves, which is the exact shape D-034 exists to
stop.
"""

from __future__ import annotations

import re
from typing import Any

from .answer import AnswerSpec, Part, compare_answer
from .extract import NUM_RE, numbers
from .kinds import COMPARATORS, PropertySlot, compare_label_tuple
from .normalize import answer_span, prepare
from .verdict import Verdict

# --------------------------------------------------------------------------
# Part selection for a composite answer.
# --------------------------------------------------------------------------


def nth_quantity(text: str, i: int, span: bool = True) -> str:
    """The slice of ``text`` carrying its ``i``-th asserted quantity.

    A ``multipart`` gold states its parts in order -- ``- CH4(g): 9.83 mol`` and
    so on -- so part ``i`` is the region running from the end of quantity
    ``i-1`` to the end of quantity ``i``, and it contains **exactly one**
    asserted quantity.

    **The slice ends at the number, not 16 characters past it** (Reviewer E,
    E-8). The tail was there to carry the unit along with the value; what it
    actually carried, whenever quantity ``i+1`` fell inside it, was the *next*
    number -- so the part comparator's answer rule selected that instead, two
    parts compared the same quantity, and one asserted quantity was compared by
    nothing at all. **6 of the 31 multipart bindings did this on real gold**,
    and two real gold pairs were accepted for each other because of it:
    ``hydrostatic_pressure_at_depth`` seeds 157/237 compared the pressure twice
    and never the depth; ``rackett_equation_volume`` accepted Carbon dioxide's
    answer for Methane.

    Ending the slice at the number costs nothing now, because parts no longer
    carry a declared unit (E-4).

    Uses the same filtered tokeniser as the extraction rule (``extract.numbers``),
    so a digit inside ``m^2`` or ``C4H10`` cannot be mistaken for a part.
    """
    t = prepare(answer_span(text)[0]) if span else text
    ms = numbers(t)
    if i >= len(ms):
        return ""
    start = ms[i - 1].end() if i else 0
    return t[start:ms[i].end()]


def numeric_parts(n: int, unit: str | None = None) -> AnswerSpec:
    """An ordered, all-required answer of ``n`` numeric parts.

    ``unit`` is accepted and **deliberately ignored** (Reviewer E, E-4). A
    multipart answer's parts carry *different* units -- that is what makes them
    different parts -- so one declared unit applied to all ``n`` of them is
    wrong for at least ``n-1``. It made 13 multipart templates MISMATCH a
    verbatim copy of their own gold: ``volumetric_flow_rate`` declared ``m/s``
    and then rejected its own flow rate in ``m^3/s`` for not being it.

    The parameter is kept rather than removed so the call sites read the same
    and the omission is visible here instead of being silent at each of them.
    Per-part units are a real thing to want and they need a per-part
    declaration, which is Phase 6's (D6.11).
    """
    return AnswerSpec(
        parts=[
            Part(name=f"q{i + 1}", kind="numeric", options={"extract": False},
                 select=(lambda s, i=i: nth_quantity(s, i)))
            for i in range(n)
        ],
        mode="all",
    )


def composite_spec(parts_desc: list[dict[str, Any]], mode: str = "all") -> AnswerSpec:
    """Build an ``AnswerSpec`` from a **declarative** part list (D6.9).

    ``numeric_parts`` builds ``n`` parts that are all ``numeric``.  A mixed
    answer -- a quantity *and* a class -- needs parts of different kinds, which
    is D-047's composition read literally rather than a seventh ``kind``.  No
    comparator changes: each part still dispatches to one of the six.

    Each descriptor is ``{name, kind, options?, select?}`` and every value is a
    plain string, dict or list, because this table is **generated** and is
    written out through ``repr`` -- a lambda would not survive the round trip.
    ``select`` names the slicer rather than being one:

    ``"answer"`` (default)
        the part sees the whole span and its own comparator locates the answer
        -- the numeric answer rule (D5.6), or a categorical label set.
    ``"quantity:i"``
        the part sees the slice carrying the ``i``-th asserted quantity, via
        ``nth_quantity``.  ``extract`` is forced off, exactly as
        ``numeric_parts`` does, because the slice is already extracted.
    """
    parts: list[Part] = []
    for d in parts_desc:
        opts = dict(d.get("options") or {})
        sel = d.get("select") or "answer"
        if sel == "answer":
            select = None
        elif sel.startswith("quantity:"):
            i = int(sel.split(":", 1)[1])
            select = (lambda s, i=i: nth_quantity(s, i))
            opts["extract"] = False
        else:
            raise ValueError(f"unknown part selector {sel!r}")
        parts.append(Part(name=d["name"], kind=d["kind"], options=opts, select=select))
    return AnswerSpec(parts=parts, mode=mode)


# --------------------------------------------------------------------------
# Structural comparators used by the derived bindings.
#
# These are not new `kind`s.  `vector` and `sequence` below dispatch to the six
# via `COMPARATORS`; the two helpers here only *locate* the answer's structure
# inside gold's prose so the right comparator sees the right substring.
# --------------------------------------------------------------------------

#: The component list.  ``<...>`` is what gold writes; ``[...]`` is what models
#: write, and both are the same answer (Reviewer E, R2-F1).
#: A component list.  ``<...>`` is what gold writes and ``[...]`` is what
#: models write -- but a bracket group is only a vector if it is a LIST of
#: numbers, so a comma and numeric-only content are both required.  Without
#: that, ordinary prose brackets read as vectors.
_ANGLE_RE = re.compile(
    r"<[^>]*\d[^>]*>"
    r"|\[[^\]a-zA-Z]*\d[^\]a-zA-Z]*,[^\]a-zA-Z]*\d[^\]a-zA-Z]*\]")

#: A unit-vector term.  **Gold writes ``x_hat``; every archived model answer
#: writes ``x̂``** -- combining circumflex U+0302, or the precomposed ``ŷ``
#: U+0177 and ``ẑ`` U+1E91.  Reading only the ASCII form meant the comparator
#: accepted **0 of 83** real answers across the three vector templates while
#: matching gold against gold perfectly, so the identity gate could not see it
#: *by construction*.
#:
#: This is a missing **surface declaration**, which is the mechanism D4.1 |S|4.2
#: already uses for categorical labels -- not a new rule and not undecidability.
#: A unit-vector term.  The ``_hat`` marker is REQUIRED: without it the
#: pattern also matches ``9.83 mol`` and every scalar answer with a unit
#: becomes a vector.  Measured when it was optional: several templates bound
#: as `vector` and matched everything, 2,450 false accepts in 2,450 pairs.
_HAT_RE = re.compile(r"[-+]?\s*[\d.eE+-]+\s*[a-z]_hat\b")

#: Observed unit-vector notations, folded ONTO the ``_hat`` marker rather than
#: stripped.  Gold writes ``x_hat``; archived model answers write ``x̂`` -- a
#: combining circumflex (U+0302) after NFD, or a precomposed ``ŷ`` / ``ẑ``.
#: Folding the mark *away* would leave a bare letter, and a number followed by a
#: bare letter is just a quantity with a unit.
_HAT_FOLD = re.compile(r"([a-z])\u0302")
_HAT_PRECOMPOSED = {"\u0177": "y_hat", "\u1e91": "z_hat"}


def _fold_hats(text: str) -> str:
    """Fold the observed unit-vector notations onto the ``_hat`` surface."""
    import unicodedata  # noqa: PLC0415
    for src, dst in _HAT_PRECOMPOSED.items():
        text = text.replace(src, dst)
    t = unicodedata.normalize("NFD", text)
    return _HAT_FOLD.sub(r"\1_hat", t)


def vector_components(text: str) -> list[str] | None:
    """The ordered components of a vector answer, or None."""
    t = _fold_hats(prepare(answer_span(text)[0]))
    m = _ANGLE_RE.search(t)
    if m:
        return [x.group(0) for x in NUM_RE.finditer(m.group(0))]
    hats = _HAT_RE.findall(t)
    if len(hats) >= 2:
        # Two or more terms, so this is a component list rather than a stray
        # "2 x" in ordinary algebra.  A one-term match is not a vector.
        return [x.group(0) for h in hats for x in NUM_RE.finditer(h)]
    return None


def compare_vector(gold: str, candidate: str, **_: Any) -> Verdict:
    """Component-wise equality of a vector answer, at gold's display precision."""
    from .kinds import compare_numeric  # noqa: PLC0415
    from .verdict import match, mismatch, unresolved

    g, c = vector_components(gold), vector_components(candidate)
    if g is None:
        return unresolved("vector", "no vector in the gold answer span")
    if c is None:
        return unresolved("vector", "no vector in the candidate answer span")
    if len(g) != len(c):
        return mismatch("vector", f"{len(c)} components against gold's {len(g)}")
    for i, (a, b) in enumerate(zip(g, c)):
        v = compare_numeric(a, b, extract=False)
        if v.outcome == "UNRESOLVED":
            return unresolved("vector", f"component {i}: {v.reason}")
        if not v.is_match:
            return mismatch("vector", f"component {i}: {v.reason}")
    return match("vector", gold_canonical=",".join(g), cand_canonical=",".join(c))


#: Bindings dispatch through these.  `multipart` is not a kind (D-047): it is an
#: `AnswerSpec` built by `numeric_parts`.
EXTRA_COMPARATORS = {"vector": compare_vector}

#: Referenced only by the GENERATED tables below, whose reprs embed it.
_ = PropertySlot


def compare_template(template_id: str, gold: str, candidate: str) -> Verdict:
    """Score a candidate for ``template_id`` under its declared binding.

    Raises ``KeyError`` for an unbound template, deliberately: scoring an
    unbound template with defaults is exactly what produced N3, and a silent
    default here would put it back.
    """
    b = BINDINGS[template_id]
    kind, opts = b["kind"], dict(b.get("options", {}))
    if kind == "multipart":
        return compare_answer(gold, candidate,
                              numeric_parts(opts["n"], opts.get("unit")))
    if kind == "composite":
        # D6.9: parts of different kinds.  Still not a seventh `kind` -- every
        # part below dispatches to one of the six through `compare_answer`.
        return compare_answer(gold, candidate,
                              composite_spec(opts["parts"], opts.get("mode", "all")))
    if kind in EXTRA_COMPARATORS:
        return EXTRA_COMPARATORS[kind](gold, candidate, **opts)
    if kind == "categorical[tuple]":
        return compare_label_tuple(gold, candidate, **opts)
    return COMPARATORS[kind](gold, candidate, **opts)


# --------------------------------------------------------------------------
# GENERATED DATA -- do not edit by hand.
# `python -m tests.comparators.bindings --regenerate` rewrites everything below.
# --------------------------------------------------------------------------

# BEGIN GENERATED
DECLARED_UNITS: dict[str, str] = {
    'template_absorbing_chain_time_to_failure': 'weeks',
    'template_adiabatic_flame_temperature': 'K',
    'template_aliased_frequency_identification': 'Hz',
    'template_angle_of_twist': 'degrees',
    'template_annulus_flowrate': 'm^3/s',
    'template_aoq_ati_rectifying': 'units',
    'template_arl_beta_mean_shift': 'subgroups',
    'template_autocorrelation_rect_pulse': 'otherwise',
    'template_basic_buoyant_force': 'kN',
    'template_basic_eoq': 'units',
    'template_batch_moles_vs_conversion': 'mol',
    'template_batch_reactor_first_order': 'seconds',
    'template_batch_reactor_second_order': 'seconds',
    'template_batch_reactor_zero_order': 'seconds',
    'template_beam_internal_moment': 'kN*m',
    'template_beam_support_reactions': 'kN',
    'template_best_hydraulic_rectangular_section': 'm',
    'template_borrow_pit_fill_volume': 'm^3',
    'template_cantilever_double_integration': 'mm',
    'template_chase_vs_level_aggregate': 'dollars',
    'template_coaxial_capacitance': 'pF/m',
    'template_composite_shafts_series': 'degrees',
    'template_constant_head_permeability': 'cm/s',
    'template_continuous_to_discrete_conversion': 'rad/sample',
    'template_coulombs_law': 'N',
    'template_cstr_volume_basic': 'liters',
    'template_effective_stress_profile': 'kPa',
    'template_epq_finite_production': 'units',
    'template_equivalent_stiffness_frequency': 'rad/s',
    'template_euclidean_distance_binary': 'e-05',
    'template_exponential_mttf_topology': 'hours',
    'template_falling_film_max_velocity': 'm/s',
    'template_floating_object_submersion_depth': 'm',
    'template_flow_system_molar_flow_rates': 'mol/min',
    'template_fluid_particle_acceleration': 'm/s^2',
    'template_force_method_continuous_beam': 'kN',
    'template_gas_phase_concentration': 'mol/dm^3',
    'template_gas_viscosity_kinetic_theory': 'Pa*s',
    'template_gauss_law_symmetric': 'V/m',
    'template_hagen_poiseuille_flowrate': 'm^3/s',
    'template_hydraulic_jump_energy_loss': 'm',
    'template_hydrostatic_force_on_plane': 'kN',
    'template_hydrostatic_pressure_at_depth': 'kPa',
    'template_ideal_gas_volume': 'liters',
    'template_influence_line_max_reaction': 'kN',
    'template_latent_heat_vaporization': 'kJ',
    'template_levenspiel_plot_interpretation': 'L',
    'template_limiting_reactant': 'mol',
    'template_line_balancing_heuristic': 'percent',
    'template_linear_reservoir_routing_step': 'm^3/s',
    'template_lorentz_force': 'N',
    'template_manning_rectangular_discharge': 'm^3/s',
    'template_manning_trapezoidal_velocity': 'm/s',
    'template_max_hump_height_no_choking': 'm',
    'template_mm1_time_in_system': 'minutes',
    'template_mm1k_finite_capacity': 'minutes',
    'template_mmc_waiting_time': 'minutes',
    'template_newsvendor_normal_demand': 'units',
    'template_newtons_law_shear_stress': 'N',
    'template_normal_depth_iteration': 'm',
    'template_nyquist_rate_determination': 'Hz',
    'template_particle_pathline': 'm',
    'template_pfr_volume_changing_rate': 'liters',
    'template_phasor_addition': 'deg',
    'template_power_law_fluid_shear': 'Pa*s',
    'template_primary_consolidation_settlement': 'mm',
    'template_qr_policy_one_iteration': 'units',
    'template_quantity_discount_all_units': 'units',
    'template_rackett_equation_volume': 'cm^3/mol',
    'template_rational_method_peak_flow': 'm^3/s',
    'template_reorder_point_lead_time': 'units',
    'template_rotating_unbalance': 'mm',
    'template_safety_stock_reorder_point': 'units',
    'template_scs_curve_number_runoff': 'mm',
    'template_sensible_heat_constant_cp': 'kJ',
    'template_sensible_heat_temp_dependent_cp': 'kJ',
    'template_server_configuration_selection': 'minutes',
    'template_shaft_design_power': 'mm',
    'template_shear_stress_torsion': 'MPa',
    'template_sigma_reduction_for_cpk': 'percent',
    'template_slope_deflection_end_moment': 'kN*m',
    'template_statically_indeterminate_shaft': 'N.m',
    'template_superposition_electric_field': 'N/C',
    'template_system_properties': 'N.s/m',
    'template_takt_time_line_efficiency': 'percent',
    'template_terzaghi_strip_footing_bearing': 'kPa',
    'template_time_rate_of_consolidation': 'years',
    'template_time_to_phasor': 'degrees',
    'template_truss_method_of_joints': 'kN',
    'template_two_phase_specific_volume': 'm^3/kg',
    'template_undamped_natural_frequency_torsional': 'seconds',
    'template_undamped_natural_frequency_translational': 'seconds',
    'template_utube_manometer': 'kPa',
    'template_vdw_solve_for_pressure': 'bar',
    'template_vdw_solve_for_volume': 'L/mol',
    'template_vibration_isolator_design': 'N/m',
    'template_vibration_transmissibility': 'mm',
    'template_virtual_work_truss_deflection': 'mm',
    'template_volumetric_flow_rate': 'm/s',
    'template_vorticity_check': 'rad/s',
    'template_wave_equation_interpretation': 'm/s',
    'template_wave_parameters_basic': 'rad/m',
    'template_work_isothermal_virial': 'J/mol',
}

BINDINGS: dict[str, dict[str, Any]] = {
    'template_absorbing_chain_time_to_failure': {'kind': 'numeric', 'options': {}, 'note': 'exactly one asserted number', 'answer_type': 'scalar'},
    'template_adiabatic_flame_temperature': {'kind': 'numeric', 'options': {}, 'note': 'exactly one asserted number', 'answer_type': 'scalar'},
    'template_aliased_frequency_identification': {'kind': 'numeric', 'options': {}, 'note': 'exactly one asserted number', 'answer_type': 'scalar'},
    'template_angle_of_twist': {'kind': 'multipart', 'options': {'n': 2}, 'note': '2 asserted numbers in every gold span', 'answer_type': 'scalar'},
    'template_annulus_flowrate': {'kind': 'numeric', 'options': {}, 'note': 'exactly one asserted number', 'answer_type': 'scalar'},
    'template_aoq_ati_rectifying': {'kind': 'numeric', 'options': {}, 'note': 'exactly one asserted number', 'answer_type': 'multipart'},
    'template_arl_beta_mean_shift': {'kind': 'numeric', 'options': {}, 'note': 'exactly one asserted number', 'answer_type': 'multipart'},
    'template_average_energy_mqam': {'kind': 'numeric', 'options': {}, 'note': 'exactly one asserted number', 'answer_type': 'scalar'},
    'template_axial_deformation': {'kind': 'numeric', 'options': {}, 'note': 'exactly one asserted number', 'answer_type': 'scalar'},
    'template_basic_buoyant_force': {'kind': 'multipart', 'options': {'n': 2}, 'note': '2 asserted numbers in every gold span', 'answer_type': 'scalar'},
    'template_basic_eoq': {'kind': 'numeric', 'options': {}, 'note': 'exactly one asserted number', 'answer_type': 'scalar'},
    'template_basic_stress_strain': {'kind': 'multipart', 'options': {'n': 2}, 'note': '2 asserted numbers in every gold span', 'answer_type': 'multipart'},
    'template_batch_reactor_first_order': {'kind': 'numeric', 'options': {}, 'note': 'exactly one asserted number', 'answer_type': 'scalar'},
    'template_batch_reactor_second_order': {'kind': 'numeric', 'options': {}, 'note': 'exactly one asserted number', 'answer_type': 'scalar'},
    'template_batch_reactor_zero_order': {'kind': 'numeric', 'options': {}, 'note': 'exactly one asserted number', 'answer_type': 'scalar'},
    'template_beam_deflection_formula': {'kind': 'numeric', 'options': {}, 'note': 'exactly one asserted number', 'answer_type': 'scalar'},
    'template_beam_internal_moment': {'kind': 'numeric', 'options': {}, 'note': 'exactly one asserted number', 'answer_type': 'scalar'},
    'template_beam_support_reactions': {'kind': 'numeric', 'options': {}, 'note': 'exactly one asserted number', 'answer_type': 'scalar'},
    'template_ber_estimation_mary': {'kind': 'symbolic', 'options': {'symbols': ('t', 'f', 'x', 'n', 'tau')}, 'note': 'a function call in every gold span', 'answer_type': 'symbolic'},
    'template_best_hydraulic_rectangular_section': {'kind': 'numeric', 'options': {}, 'note': 'exactly one asserted number', 'answer_type': 'scalar'},
    'template_borrow_pit_fill_volume': {'kind': 'numeric', 'options': {}, 'note': 'exactly one asserted number', 'answer_type': 'scalar'},
    'template_c_chart_revision': {'kind': 'numeric', 'options': {}, 'note': 'exactly one asserted number', 'answer_type': 'scalar'},
    'template_cantilever_double_integration': {'kind': 'numeric', 'options': {}, 'note': 'exactly one asserted number', 'answer_type': 'scalar'},
    'template_cd_dc_system_analysis': {'kind': 'symbolic', 'options': {'symbols': ('t', 'f', 'x', 'n', 'tau')}, 'note': 'a function call in every gold span', 'answer_type': 'symbolic'},
    'template_chart_pair_selection': {'kind': 'numeric', 'options': {}, 'note': 'exactly one asserted number', 'answer_type': 'scalar'},
    'template_chase_vs_level_aggregate': {'kind': 'numeric', 'options': {}, 'note': 'exactly one asserted number', 'answer_type': 'multipart'},
    'template_composite_shafts_series': {'kind': 'multipart', 'options': {'n': 2}, 'note': '2 asserted numbers in every gold span', 'answer_type': 'scalar'},
    'template_constant_head_permeability': {'kind': 'numeric', 'options': {}, 'note': 'exactly one asserted number', 'answer_type': 'scalar'},
    'template_coulombs_law': {'kind': 'vector', 'options': {}, 'note': 'an angle-bracket or unit-vector form in every gold span', 'answer_type': 'vector'},
    'template_cp_cpk_from_specs': {'kind': 'numeric', 'options': {}, 'note': 'exactly one asserted number', 'answer_type': 'scalar'},
    'template_critical_depth_froude_classification': {'kind': 'numeric', 'options': {}, 'note': 'exactly one asserted number', 'answer_type': 'classification', 'partial': 'checks the quantity only; gold also states a class label, so a right number with a wrong class is credited (D6.9)'},
    'template_cstr_volume_basic': {'kind': 'numeric', 'options': {}, 'note': 'exactly one asserted number', 'answer_type': 'scalar'},
    'template_damping_classification': {'kind': 'composite', 'options': {'parts': [{'name': 'zeta', 'kind': 'numeric', 'select': 'quantity:0'}, {'name': 'regime', 'kind': 'categorical', 'options': {'labels': {'Underdamped': ['underdamped'], 'Critically Damped': ['critically damped'], 'Overdamped': ['overdamped']}}}, {'name': 'omega_d', 'kind': 'numeric', 'select': 'answer'}]}, 'note': "declared: a label set is not derivable from gold's structure", 'answer_type': 'classification'},
    'template_effective_stress_profile': {'kind': 'numeric', 'options': {}, 'note': 'exactly one asserted number', 'answer_type': 'scalar'},
    'template_epq_finite_production': {'kind': 'numeric', 'options': {}, 'note': 'exactly one asserted number', 'answer_type': 'scalar'},
    'template_equivalent_stiffness_frequency': {'kind': 'multipart', 'options': {'n': 2}, 'note': '2 asserted numbers in every gold span', 'answer_type': 'multipart'},
    'template_exponential_mttf_topology': {'kind': 'numeric', 'options': {}, 'note': 'exactly one asserted number', 'answer_type': 'scalar'},
    'template_finite_convolution': {'kind': 'sequence', 'options': {}, 'note': 'a braced sequence in every gold span', 'answer_type': 'array'},
    'template_floating_object_submersion_depth': {'kind': 'numeric', 'options': {}, 'note': 'exactly one asserted number', 'answer_type': 'scalar'},
    'template_fluid_particle_acceleration': {'kind': 'multipart', 'options': {'n': 3}, 'note': '3 asserted numbers in every gold span', 'answer_type': 'vector'},
    'template_force_method_continuous_beam': {'kind': 'numeric', 'options': {}, 'note': 'exactly one asserted number', 'answer_type': 'scalar'},
    'template_ft_esd_rect_pulse': {'kind': 'symbolic', 'options': {'symbols': ('t', 'f', 'x', 'n', 'tau')}, 'note': 'a function call in every gold span', 'answer_type': 'symbolic'},
    'template_gas_phase_concentration': {'kind': 'multipart', 'options': {'n': 2}, 'note': '2 asserted numbers in every gold span', 'answer_type': 'multipart'},
    'template_heat_of_reaction_formation': {'kind': 'numeric', 'options': {}, 'note': 'exactly one asserted number', 'answer_type': 'scalar'},
    'template_hydraulic_jump_energy_loss': {'kind': 'numeric', 'options': {}, 'note': 'exactly one asserted number', 'answer_type': 'scalar'},
    'template_hydrostatic_force_on_plane': {'kind': 'multipart', 'options': {'n': 2}, 'note': '2 asserted numbers in every gold span', 'answer_type': 'multipart'},
    'template_hydrostatic_pressure_at_depth': {'kind': 'multipart', 'options': {'n': 2}, 'note': '2 asserted numbers in every gold span', 'answer_type': 'scalar'},
    'template_ideal_gas_volume': {'kind': 'numeric', 'options': {}, 'note': 'exactly one asserted number', 'answer_type': 'scalar'},
    'template_incompressible_continuity': {'kind': 'symbolic', 'options': {'symbols': ('x', 'y'), 'allow_arbitrary': True}, 'note': "declared: a label set is not derivable from gold's structure", 'answer_type': 'symbolic'},
    'template_infinite_slope_factor_of_safety': {'kind': 'numeric', 'options': {}, 'note': 'exactly one asserted number', 'answer_type': 'scalar'},
    'template_influence_line_max_reaction': {'kind': 'numeric', 'options': {}, 'note': 'exactly one asserted number', 'answer_type': 'scalar'},
    'template_latent_heat_vaporization': {'kind': 'numeric', 'options': {}, 'note': 'exactly one asserted number', 'answer_type': 'scalar'},
    'template_levenspiel_plot_interpretation': {'kind': 'multipart', 'options': {'n': 2}, 'note': '2 asserted numbers in every gold span', 'answer_type': 'array'},
    'template_line_balancing_heuristic': {'kind': 'numeric', 'options': {}, 'note': 'exactly one asserted number', 'answer_type': 'array'},
    'template_linear_reservoir_routing_step': {'kind': 'numeric', 'options': {}, 'note': 'exactly one asserted number', 'answer_type': 'scalar'},
    'template_logarithmic_decrement': {'kind': 'multipart', 'options': {'n': 2}, 'note': '2 asserted numbers in every gold span', 'answer_type': 'multipart'},
    'template_lorentz_force': {'kind': 'vector', 'options': {}, 'note': 'an angle-bracket or unit-vector form in every gold span', 'answer_type': 'vector'},
    'template_lowpass_equivalent_bandpass': {'kind': 'multipart', 'options': {'n': 4}, 'note': '4 asserted numbers in every gold span', 'answer_type': 'multipart'},
    'template_manning_rectangular_discharge': {'kind': 'numeric', 'options': {}, 'note': 'exactly one asserted number', 'answer_type': 'scalar'},
    'template_manning_trapezoidal_velocity': {'kind': 'numeric', 'options': {}, 'note': 'exactly one asserted number', 'answer_type': 'scalar'},
    'template_max_hump_height_no_choking': {'kind': 'numeric', 'options': {}, 'note': 'exactly one asserted number', 'answer_type': 'scalar'},
    'template_mean_variance': {'kind': 'multipart', 'options': {'n': 2}, 'note': '2 asserted numbers in every gold span', 'answer_type': 'multipart'},
    'template_mm1_time_in_system': {'kind': 'numeric', 'options': {}, 'note': 'exactly one asserted number', 'answer_type': 'scalar'},
    'template_mm1k_finite_capacity': {'kind': 'numeric', 'options': {}, 'note': 'exactly one asserted number', 'answer_type': 'scalar'},
    'template_mmc_waiting_time': {'kind': 'numeric', 'options': {}, 'note': 'exactly one asserted number', 'answer_type': 'scalar'},
    'template_multi_segment_rod': {'kind': 'numeric', 'options': {}, 'note': 'exactly one asserted number', 'answer_type': 'scalar'},
    'template_newsvendor_normal_demand': {'kind': 'numeric', 'options': {}, 'note': 'exactly one asserted number', 'answer_type': 'scalar'},
    'template_newtons_law_shear_stress': {'kind': 'multipart', 'options': {'n': 2}, 'note': '2 asserted numbers in every gold span', 'answer_type': 'multipart'},
    'template_normal_depth_iteration': {'kind': 'numeric', 'options': {}, 'note': 'exactly one asserted number', 'answer_type': 'array'},
    'template_nyquist_rate_determination': {'kind': 'numeric', 'options': {}, 'note': 'exactly one asserted number', 'answer_type': 'scalar'},
    'template_p_chart_limits_floor': {'kind': 'numeric', 'options': {}, 'note': 'exactly one asserted number', 'answer_type': 'scalar'},
    'template_particle_pathline': {'kind': 'multipart', 'options': {'n': 3}, 'note': '3 asserted numbers in every gold span', 'answer_type': 'vector'},
    'template_pfr_volume_changing_rate': {'kind': 'numeric', 'options': {}, 'note': 'exactly one asserted number', 'answer_type': 'scalar'},
    'template_phase_relations_degree_of_saturation': {'kind': 'numeric', 'options': {}, 'note': 'exactly one asserted number', 'answer_type': 'scalar'},
    'template_phasor_addition': {'kind': 'symbolic', 'options': {'symbols': ('t', 'f', 'x', 'n', 'tau')}, 'note': 'a function call in every gold span', 'answer_type': 'symbolic'},
    'template_poisson_event_count': {'kind': 'numeric', 'options': {}, 'note': 'exactly one asserted number', 'answer_type': 'scalar'},
    'template_poissons_ratio': {'kind': 'multipart', 'options': {'n': 2}, 'note': '2 asserted numbers in every gold span', 'answer_type': 'multipart'},
    'template_power_law_fluid_shear': {'kind': 'multipart', 'options': {'n': 2}, 'note': '2 asserted numbers in every gold span', 'answer_type': 'multipart'},
    'template_primary_consolidation_settlement': {'kind': 'numeric', 'options': {}, 'note': 'exactly one asserted number', 'answer_type': 'scalar'},
    'template_qr_policy_one_iteration': {'kind': 'numeric', 'options': {}, 'note': 'exactly one asserted number', 'answer_type': 'scalar'},
    'template_quantity_discount_all_units': {'kind': 'numeric', 'options': {}, 'note': 'exactly one asserted number', 'answer_type': 'scalar'},
    'template_rackett_equation_volume': {'kind': 'multipart', 'options': {'n': 2}, 'note': '2 asserted numbers in every gold span', 'answer_type': 'scalar'},
    'template_rational_method_peak_flow': {'kind': 'numeric', 'options': {}, 'note': 'exactly one asserted number', 'answer_type': 'scalar'},
    'template_relative_density_of_sand': {'kind': 'numeric', 'options': {}, 'note': 'exactly one asserted number', 'answer_type': 'scalar'},
    'template_reorder_point_lead_time': {'kind': 'numeric', 'options': {}, 'note': 'exactly one asserted number', 'answer_type': 'scalar'},
    'template_reynolds_number_flow_regime': {'kind': 'composite', 'options': {'parts': [{'name': 'Re', 'kind': 'numeric', 'select': 'answer'}, {'name': 'regime', 'kind': 'categorical', 'options': {'labels': {'laminar': ['laminar'], 'transitional': ['transitional'], 'turbulent': ['turbulent']}}}]}, 'note': "declared: a label set is not derivable from gold's structure", 'answer_type': 'classification'},
    'template_rotating_unbalance': {'kind': 'numeric', 'options': {}, 'note': 'exactly one asserted number', 'answer_type': 'scalar'},
    'template_safety_stock_reorder_point': {'kind': 'numeric', 'options': {}, 'note': 'exactly one asserted number', 'answer_type': 'scalar'},
    'template_scs_curve_number_runoff': {'kind': 'numeric', 'options': {}, 'note': 'exactly one asserted number', 'answer_type': 'scalar'},
    'template_sensible_heat_constant_cp': {'kind': 'numeric', 'options': {}, 'note': 'exactly one asserted number', 'answer_type': 'scalar'},
    'template_sensible_heat_temp_dependent_cp': {'kind': 'numeric', 'options': {}, 'note': 'exactly one asserted number', 'answer_type': 'scalar'},
    'template_server_configuration_selection': {'kind': 'numeric', 'options': {}, 'note': 'exactly one asserted number', 'answer_type': 'multipart'},
    'template_shaft_design_power': {'kind': 'numeric', 'options': {}, 'note': 'exactly one asserted number', 'answer_type': 'scalar'},
    'template_shear_stress_torsion': {'kind': 'numeric', 'options': {}, 'note': 'exactly one asserted number', 'answer_type': 'scalar'},
    'template_sigma_reduction_for_cpk': {'kind': 'numeric', 'options': {}, 'note': 'exactly one asserted number', 'answer_type': 'scalar'},
    'template_signal_operations': {'kind': 'sequence', 'options': {'require_origin': True}, 'note': "declared: a label set is not derivable from gold's structure", 'answer_type': 'array'},
    'template_single_sampling_oc_point': {'kind': 'numeric', 'options': {}, 'note': 'exactly one asserted number', 'answer_type': 'scalar'},
    'template_slope_deflection_end_moment': {'kind': 'numeric', 'options': {}, 'note': 'exactly one asserted number', 'answer_type': 'scalar'},
    'template_statically_indeterminate': {'kind': 'multipart', 'options': {'n': 4}, 'note': '4 asserted numbers in every gold span', 'answer_type': 'multipart'},
    'template_statically_indeterminate_shaft': {'kind': 'multipart', 'options': {'n': 2}, 'note': '2 asserted numbers in every gold span', 'answer_type': 'multipart'},
    'template_superposition_electric_field': {'kind': 'vector', 'options': {}, 'note': 'an angle-bracket or unit-vector form in every gold span', 'answer_type': 'vector'},
    'template_system_properties': {'kind': 'multipart', 'options': {'n': 3}, 'note': '3 asserted numbers in every gold span', 'answer_type': 'multipart'},
    'template_system_properties_memory_causality': {'kind': 'categorical[tuple]', 'options': {'slots': [PropertySlot(name='memoryless', positive=('memoryless', 'no memory', 'without memory', 'memory-less'), negative=('has memory', 'with memory', 'have memory')), PropertySlot(name='causal', positive=('causal',), negative=('noncausal', 'non-causal', 'anticausal', 'anti-causal'))]}, 'note': "declared: a label set is not derivable from gold's structure", 'answer_type': 'classification'},
    'template_system_property_linearity': {'kind': 'categorical', 'options': {'labels': {'linear': ['linear'], 'not linear': ['nonlinear', 'non-linear', 'nonlinearity', 'non-linearity']}}, 'note': "declared: a label set is not derivable from gold's structure", 'answer_type': 'classification'},
    'template_system_reliability_topology': {'kind': 'numeric', 'options': {}, 'note': 'exactly one asserted number', 'answer_type': 'scalar'},
    'template_takt_time_line_efficiency': {'kind': 'numeric', 'options': {}, 'note': 'exactly one asserted number', 'answer_type': 'scalar'},
    'template_terzaghi_strip_footing_bearing': {'kind': 'numeric', 'options': {}, 'note': 'exactly one asserted number', 'answer_type': 'scalar'},
    'template_time_rate_of_consolidation': {'kind': 'numeric', 'options': {}, 'note': 'exactly one asserted number', 'answer_type': 'scalar'},
    'template_time_to_phasor': {'kind': 'multipart', 'options': {'n': 4}, 'note': '4 asserted numbers in every gold span', 'answer_type': 'multipart'},
    'template_truss_method_of_joints': {'kind': 'numeric', 'options': {}, 'note': 'exactly one asserted number', 'answer_type': 'scalar'},
    'template_two_phase_specific_volume': {'kind': 'numeric', 'options': {}, 'note': 'exactly one asserted number', 'answer_type': 'scalar'},
    'template_two_state_steady_state': {'kind': 'numeric', 'options': {}, 'note': 'exactly one asserted number', 'answer_type': 'vector'},
    'template_two_step_transition_probability': {'kind': 'numeric', 'options': {}, 'note': 'exactly one asserted number', 'answer_type': 'scalar'},
    'template_undamped_natural_frequency_torsional': {'kind': 'multipart', 'options': {'n': 3}, 'note': '3 asserted numbers in every gold span', 'answer_type': 'multipart'},
    'template_undamped_natural_frequency_translational': {'kind': 'multipart', 'options': {'n': 3}, 'note': '3 asserted numbers in every gold span', 'answer_type': 'multipart'},
    'template_undamped_response_initial_conditions': {'kind': 'symbolic', 'options': {'symbols': ('t', 'f', 'x', 'n', 'tau')}, 'note': 'a function call in every gold span', 'answer_type': 'symbolic'},
    'template_upward_seepage_quick_condition': {'kind': 'numeric', 'options': {}, 'note': 'exactly one asserted number', 'answer_type': 'scalar'},
    'template_utube_manometer': {'kind': 'numeric', 'options': {}, 'note': 'exactly one asserted number', 'answer_type': 'scalar'},
    'template_vibration_isolator_design': {'kind': 'numeric', 'options': {}, 'note': 'exactly one asserted number', 'answer_type': 'scalar'},
    'template_vibration_transmissibility': {'kind': 'multipart', 'options': {'n': 2}, 'note': '2 asserted numbers in every gold span', 'answer_type': 'multipart'},
    'template_virtual_work_truss_deflection': {'kind': 'numeric', 'options': {}, 'note': 'exactly one asserted number', 'answer_type': 'scalar'},
    'template_volumetric_flow_rate': {'kind': 'multipart', 'options': {'n': 2}, 'note': '2 asserted numbers in every gold span', 'answer_type': 'multipart'},
    'template_wave_equation_interpretation': {'kind': 'multipart', 'options': {'n': 4}, 'note': '4 asserted numbers in every gold span', 'answer_type': 'multipart'},
    'template_wave_parameters_basic': {'kind': 'multipart', 'options': {'n': 4}, 'note': '4 asserted numbers in every gold span', 'answer_type': 'multipart'},
    'template_work_isothermal_virial': {'kind': 'numeric', 'options': {}, 'note': 'exactly one asserted number', 'answer_type': 'scalar'},
    'template_xbar_known_sigma_classification': {'kind': 'numeric', 'options': {}, 'note': 'exactly one asserted number', 'answer_type': 'scalar'},
    'template_xbar_r_control_limits': {'kind': 'numeric', 'options': {}, 'note': 'exactly one asserted number', 'answer_type': 'scalar'},
}

UNBOUND: dict[str, str] = {
    'template_autocorrelation_rect_pulse': 'part(s) [4, 6] read a constant on every instance -- compared by nothing (E-9)',
    'template_batch_moles_vs_conversion': 'the asserted-number count varies across instances: [7, 8]',
    'template_bpsk_energy_basis': 'rejects a verbatim copy of gold on 50/50 seeds (UNRESOLVED: could not parse an expression: SympifyError)',
    'template_coaxial_capacitance': 'part(s) [1] read a constant on every instance -- compared by nothing (E-9)',
    'template_continuous_to_discrete_conversion': 'rejects a verbatim copy of gold on 50/50 seeds (UNRESOLVED: expression did not parse to a single expression (gold tuple, candidate tuple))',
    'template_decimation_aliasing_analysis': '34 false accepts in 2450 pairs',
    'template_euclidean_distance_binary': 'the asserted-number count varies across instances: [5, 7]',
    'template_falling_film_max_velocity': 'the asserted-number count varies across instances: [1, 2]',
    'template_flow_system_molar_flow_rates': 'the asserted-number count varies across instances: [6, 7]',
    'template_gas_viscosity_kinetic_theory': 'the asserted-number count varies across instances: [2, 3]',
    'template_gauss_law_symmetric': '4 false accepts in 2450 pairs',
    'template_hagen_poiseuille_flowrate': 'the asserted-number count varies across instances: [1, 2]',
    'template_impulse_response_from_lccde': 'the asserted-number count varies across instances: [1, 4, 5]',
    'template_kinematic_viscosity': 'the asserted-number count varies across instances: [1, 2]',
    'template_limiting_reactant': 'the asserted-number count varies across instances: [6, 7, 8]',
    'template_null_to_null_bandwidth': '24 false accepts in 2450 pairs',
    'template_pitzer_correlation_z': '1 false accepts in 2450 pairs',
    'template_signal_energy_power': '2 false accepts in 2450 pairs',
    'template_standing_wave_formation': 'rejects a verbatim copy of gold on 50/50 seeds (UNRESOLVED: expression did not parse to a single expression (gold tuple, candidate tuple))',
    'template_truss_method_of_sections': '2 false accepts in 2450 pairs',
    'template_vdw_solve_for_pressure': '2 false accepts in 2450 pairs',
    'template_vdw_solve_for_volume': 'part(s) [1, 3] read a constant on every instance -- compared by nothing (E-9)',
    'template_vorticity_check': 'the asserted-number count varies across instances: [5, 6]',
}

VALIDATION: dict[str, dict[str, int]] = {
    'template_absorbing_chain_time_to_failure': {'pairs': 2450, 'decided': 2450, 'identity_failures': 0, 'constant_parts': 0},
    'template_adiabatic_flame_temperature': {'pairs': 2450, 'decided': 2450, 'identity_failures': 0, 'constant_parts': 0},
    'template_aliased_frequency_identification': {'pairs': 2450, 'decided': 2450, 'identity_failures': 0, 'constant_parts': 0},
    'template_angle_of_twist': {'pairs': 2450, 'decided': 2450, 'identity_failures': 0, 'constant_parts': 0},
    'template_annulus_flowrate': {'pairs': 2450, 'decided': 2450, 'identity_failures': 0, 'constant_parts': 0},
    'template_aoq_ati_rectifying': {'pairs': 2450, 'decided': 2450, 'identity_failures': 0, 'constant_parts': 0},
    'template_arl_beta_mean_shift': {'pairs': 2450, 'decided': 2450, 'identity_failures': 0, 'constant_parts': 0},
    'template_autocorrelation_rect_pulse': {'pairs': 2450, 'decided': 2450, 'identity_failures': 0, 'constant_parts': 2},
    'template_average_energy_mqam': {'pairs': 2450, 'decided': 2450, 'identity_failures': 0, 'constant_parts': 0},
    'template_axial_deformation': {'pairs': 2450, 'decided': 2450, 'identity_failures': 0, 'constant_parts': 0},
    'template_basic_buoyant_force': {'pairs': 2450, 'decided': 2450, 'identity_failures': 0, 'constant_parts': 0},
    'template_basic_eoq': {'pairs': 2450, 'decided': 2450, 'identity_failures': 0, 'constant_parts': 0},
    'template_basic_stress_strain': {'pairs': 2450, 'decided': 2450, 'identity_failures': 0, 'constant_parts': 0},
    'template_batch_reactor_first_order': {'pairs': 2450, 'decided': 2450, 'identity_failures': 0, 'constant_parts': 0},
    'template_batch_reactor_second_order': {'pairs': 2450, 'decided': 2450, 'identity_failures': 0, 'constant_parts': 0},
    'template_batch_reactor_zero_order': {'pairs': 2450, 'decided': 2450, 'identity_failures': 0, 'constant_parts': 0},
    'template_beam_deflection_formula': {'pairs': 2450, 'decided': 2450, 'identity_failures': 0, 'constant_parts': 0},
    'template_beam_internal_moment': {'pairs': 2450, 'decided': 2450, 'identity_failures': 0, 'constant_parts': 0},
    'template_beam_support_reactions': {'pairs': 2450, 'decided': 2450, 'identity_failures': 0, 'constant_parts': 0},
    'template_ber_estimation_mary': {'pairs': 2450, 'decided': 2450, 'identity_failures': 0, 'constant_parts': 0},
    'template_best_hydraulic_rectangular_section': {'pairs': 2450, 'decided': 2450, 'identity_failures': 0, 'constant_parts': 0},
    'template_borrow_pit_fill_volume': {'pairs': 2450, 'decided': 2450, 'identity_failures': 0, 'constant_parts': 0},
    'template_bpsk_energy_basis': {'pairs': 2450, 'unresolved': 2450, 'decided': 0, 'identity_failures': 50, 'constant_parts': 0},
    'template_c_chart_revision': {'pairs': 2450, 'decided': 2450, 'identity_failures': 0, 'constant_parts': 0},
    'template_cantilever_double_integration': {'pairs': 2450, 'decided': 2450, 'identity_failures': 0, 'constant_parts': 0},
    'template_cd_dc_system_analysis': {'pairs': 2450, 'decided': 2450, 'identity_failures': 0, 'constant_parts': 0},
    'template_chart_pair_selection': {'pairs': 2450, 'decided': 2450, 'identity_failures': 0, 'constant_parts': 0},
    'template_chase_vs_level_aggregate': {'pairs': 2450, 'decided': 2450, 'identity_failures': 0, 'constant_parts': 0},
    'template_coaxial_capacitance': {'pairs': 2450, 'decided': 2450, 'identity_failures': 0, 'constant_parts': 1},
    'template_composite_shafts_series': {'pairs': 2450, 'decided': 2450, 'identity_failures': 0, 'constant_parts': 0},
    'template_constant_head_permeability': {'pairs': 2450, 'decided': 2450, 'identity_failures': 0, 'constant_parts': 0},
    'template_continuous_to_discrete_conversion': {'pairs': 2450, 'unresolved': 2450, 'decided': 0, 'identity_failures': 50, 'constant_parts': 0},
    'template_coulombs_law': {'pairs': 2450, 'decided': 2450, 'identity_failures': 0, 'constant_parts': 0},
    'template_cp_cpk_from_specs': {'pairs': 2450, 'decided': 2450, 'identity_failures': 0, 'constant_parts': 0},
    'template_critical_depth_froude_classification': {'pairs': 2450, 'decided': 2450, 'identity_failures': 0, 'constant_parts': 0},
    'template_cstr_volume_basic': {'pairs': 2450, 'decided': 2450, 'identity_failures': 0, 'constant_parts': 0},
    'template_damping_classification': {'pairs': 2450, 'decided': 2450, 'identity_failures': 0, 'constant_parts': 0},
    'template_decimation_aliasing_analysis': {'pairs': 2450, 'false_accepts': 34, 'decided': 2450, 'identity_failures': 0, 'constant_parts': 0},
    'template_effective_stress_profile': {'pairs': 2450, 'decided': 2450, 'identity_failures': 0, 'constant_parts': 0},
    'template_epq_finite_production': {'pairs': 2450, 'decided': 2450, 'identity_failures': 0, 'constant_parts': 0},
    'template_equivalent_stiffness_frequency': {'pairs': 2450, 'decided': 2450, 'identity_failures': 0, 'constant_parts': 0},
    'template_exponential_mttf_topology': {'pairs': 2450, 'decided': 2450, 'identity_failures': 0, 'constant_parts': 0},
    'template_finite_convolution': {'pairs': 2450, 'decided': 2450, 'identity_failures': 0, 'constant_parts': 0},
    'template_floating_object_submersion_depth': {'pairs': 2450, 'decided': 2450, 'identity_failures': 0, 'constant_parts': 0},
    'template_fluid_particle_acceleration': {'pairs': 2450, 'decided': 2450, 'identity_failures': 0, 'constant_parts': 0},
    'template_force_method_continuous_beam': {'pairs': 2450, 'decided': 2450, 'identity_failures': 0, 'constant_parts': 0},
    'template_ft_esd_rect_pulse': {'pairs': 2450, 'decided': 2450, 'identity_failures': 0, 'constant_parts': 0},
    'template_gas_phase_concentration': {'pairs': 2450, 'decided': 2450, 'identity_failures': 0, 'constant_parts': 0},
    'template_gauss_law_symmetric': {'pairs': 2450, 'false_accepts': 4, 'decided': 2450, 'identity_failures': 0, 'constant_parts': 0},
    'template_heat_of_reaction_formation': {'pairs': 2450, 'decided': 2450, 'identity_failures': 0, 'constant_parts': 0},
    'template_hydraulic_jump_energy_loss': {'pairs': 2450, 'decided': 2450, 'identity_failures': 0, 'constant_parts': 0},
    'template_hydrostatic_force_on_plane': {'pairs': 2450, 'decided': 2450, 'identity_failures': 0, 'constant_parts': 0},
    'template_hydrostatic_pressure_at_depth': {'pairs': 2450, 'decided': 2450, 'identity_failures': 0, 'constant_parts': 0},
    'template_ideal_gas_volume': {'pairs': 2450, 'decided': 2450, 'identity_failures': 0, 'constant_parts': 0},
    'template_incompressible_continuity': {'pairs': 2450, 'false_rejects': 4, 'decided': 2450, 'identity_failures': 0, 'constant_parts': 0},
    'template_infinite_slope_factor_of_safety': {'pairs': 2450, 'decided': 2450, 'identity_failures': 0, 'constant_parts': 0},
    'template_influence_line_max_reaction': {'pairs': 2450, 'decided': 2450, 'identity_failures': 0, 'constant_parts': 0},
    'template_latent_heat_vaporization': {'pairs': 2450, 'decided': 2450, 'identity_failures': 0, 'constant_parts': 0},
    'template_levenspiel_plot_interpretation': {'pairs': 2450, 'decided': 2450, 'identity_failures': 0, 'constant_parts': 0},
    'template_line_balancing_heuristic': {'pairs': 2450, 'decided': 2450, 'identity_failures': 0, 'constant_parts': 0},
    'template_linear_reservoir_routing_step': {'pairs': 2450, 'decided': 2450, 'identity_failures': 0, 'constant_parts': 0},
    'template_logarithmic_decrement': {'pairs': 2450, 'decided': 2450, 'identity_failures': 0, 'constant_parts': 0},
    'template_lorentz_force': {'pairs': 2450, 'decided': 2450, 'identity_failures': 0, 'constant_parts': 0},
    'template_lowpass_equivalent_bandpass': {'pairs': 2450, 'decided': 2450, 'identity_failures': 0, 'constant_parts': 0},
    'template_manning_rectangular_discharge': {'pairs': 2450, 'decided': 2450, 'identity_failures': 0, 'constant_parts': 0},
    'template_manning_trapezoidal_velocity': {'pairs': 2450, 'decided': 2450, 'identity_failures': 0, 'constant_parts': 0},
    'template_max_hump_height_no_choking': {'pairs': 2450, 'decided': 2450, 'identity_failures': 0, 'constant_parts': 0},
    'template_mean_variance': {'pairs': 2450, 'decided': 2450, 'identity_failures': 0, 'constant_parts': 0},
    'template_mm1_time_in_system': {'pairs': 2450, 'decided': 2450, 'identity_failures': 0, 'constant_parts': 0},
    'template_mm1k_finite_capacity': {'pairs': 2450, 'decided': 2450, 'identity_failures': 0, 'constant_parts': 0},
    'template_mmc_waiting_time': {'pairs': 2450, 'decided': 2450, 'identity_failures': 0, 'constant_parts': 0},
    'template_multi_segment_rod': {'pairs': 2450, 'decided': 2450, 'identity_failures': 0, 'constant_parts': 0},
    'template_newsvendor_normal_demand': {'pairs': 2450, 'decided': 2450, 'identity_failures': 0, 'constant_parts': 0},
    'template_newtons_law_shear_stress': {'pairs': 2450, 'decided': 2450, 'identity_failures': 0, 'constant_parts': 0},
    'template_normal_depth_iteration': {'pairs': 2450, 'decided': 2450, 'identity_failures': 0, 'constant_parts': 0},
    'template_null_to_null_bandwidth': {'pairs': 2450, 'false_accepts': 24, 'decided': 2450, 'identity_failures': 0, 'constant_parts': 0},
    'template_nyquist_rate_determination': {'pairs': 2450, 'decided': 2450, 'identity_failures': 0, 'constant_parts': 0},
    'template_p_chart_limits_floor': {'pairs': 2450, 'decided': 2450, 'identity_failures': 0, 'constant_parts': 0},
    'template_particle_pathline': {'pairs': 2450, 'decided': 2450, 'identity_failures': 0, 'constant_parts': 0},
    'template_pfr_volume_changing_rate': {'pairs': 2450, 'decided': 2450, 'identity_failures': 0, 'constant_parts': 0},
    'template_phase_relations_degree_of_saturation': {'pairs': 2450, 'decided': 2450, 'identity_failures': 0, 'constant_parts': 0},
    'template_phasor_addition': {'pairs': 2450, 'decided': 2450, 'identity_failures': 0, 'constant_parts': 0},
    'template_pitzer_correlation_z': {'pairs': 2450, 'false_accepts': 1, 'decided': 2450, 'identity_failures': 0, 'constant_parts': 0},
    'template_poisson_event_count': {'pairs': 2450, 'decided': 2450, 'identity_failures': 0, 'constant_parts': 0},
    'template_poissons_ratio': {'pairs': 2450, 'decided': 2450, 'identity_failures': 0, 'constant_parts': 0},
    'template_power_law_fluid_shear': {'pairs': 2450, 'decided': 2450, 'identity_failures': 0, 'constant_parts': 0},
    'template_primary_consolidation_settlement': {'pairs': 2450, 'decided': 2450, 'identity_failures': 0, 'constant_parts': 0},
    'template_qr_policy_one_iteration': {'pairs': 2450, 'decided': 2450, 'identity_failures': 0, 'constant_parts': 0},
    'template_quantity_discount_all_units': {'pairs': 2450, 'decided': 2450, 'identity_failures': 0, 'constant_parts': 0},
    'template_rackett_equation_volume': {'pairs': 2450, 'decided': 2450, 'identity_failures': 0, 'constant_parts': 0},
    'template_rational_method_peak_flow': {'pairs': 2450, 'decided': 2450, 'identity_failures': 0, 'constant_parts': 0},
    'template_relative_density_of_sand': {'pairs': 2450, 'decided': 2450, 'identity_failures': 0, 'constant_parts': 0},
    'template_reorder_point_lead_time': {'pairs': 2450, 'decided': 2450, 'identity_failures': 0, 'constant_parts': 0},
    'template_reynolds_number_flow_regime': {'pairs': 2450, 'decided': 2450, 'identity_failures': 0, 'constant_parts': 0},
    'template_rotating_unbalance': {'pairs': 2450, 'decided': 2450, 'identity_failures': 0, 'constant_parts': 0},
    'template_safety_stock_reorder_point': {'pairs': 2450, 'decided': 2450, 'identity_failures': 0, 'constant_parts': 0},
    'template_scs_curve_number_runoff': {'pairs': 2450, 'decided': 2450, 'identity_failures': 0, 'constant_parts': 0},
    'template_sensible_heat_constant_cp': {'pairs': 2450, 'decided': 2450, 'identity_failures': 0, 'constant_parts': 0},
    'template_sensible_heat_temp_dependent_cp': {'pairs': 2450, 'decided': 2450, 'identity_failures': 0, 'constant_parts': 0},
    'template_server_configuration_selection': {'pairs': 2450, 'decided': 2450, 'identity_failures': 0, 'constant_parts': 0},
    'template_shaft_design_power': {'pairs': 2450, 'decided': 2450, 'identity_failures': 0, 'constant_parts': 0},
    'template_shear_stress_torsion': {'pairs': 2450, 'decided': 2450, 'identity_failures': 0, 'constant_parts': 0},
    'template_sigma_reduction_for_cpk': {'pairs': 2450, 'decided': 2450, 'identity_failures': 0, 'constant_parts': 0},
    'template_signal_energy_power': {'pairs': 2450, 'false_accepts': 2, 'decided': 2450, 'identity_failures': 0, 'constant_parts': 0},
    'template_signal_operations': {'pairs': 2450, 'decided': 2450, 'identity_failures': 0, 'constant_parts': 0},
    'template_single_sampling_oc_point': {'pairs': 2450, 'decided': 2450, 'identity_failures': 0, 'constant_parts': 0},
    'template_slope_deflection_end_moment': {'pairs': 2450, 'decided': 2450, 'identity_failures': 0, 'constant_parts': 0},
    'template_standing_wave_formation': {'pairs': 2450, 'unresolved': 2450, 'decided': 0, 'identity_failures': 50, 'constant_parts': 0},
    'template_statically_indeterminate': {'pairs': 2450, 'decided': 2450, 'identity_failures': 0, 'constant_parts': 0},
    'template_statically_indeterminate_shaft': {'pairs': 2450, 'decided': 2450, 'identity_failures': 0, 'constant_parts': 0},
    'template_superposition_electric_field': {'pairs': 2450, 'decided': 2450, 'identity_failures': 0, 'constant_parts': 0},
    'template_system_properties': {'pairs': 2450, 'decided': 2450, 'identity_failures': 0, 'constant_parts': 0},
    'template_system_properties_memory_causality': {'pairs': 2450, 'decided': 2450, 'identity_failures': 0, 'constant_parts': 0},
    'template_system_property_linearity': {'pairs': 2450, 'decided': 2450, 'identity_failures': 0, 'constant_parts': 0},
    'template_system_reliability_topology': {'pairs': 2450, 'decided': 2450, 'identity_failures': 0, 'constant_parts': 0},
    'template_takt_time_line_efficiency': {'pairs': 2450, 'decided': 2450, 'identity_failures': 0, 'constant_parts': 0},
    'template_terzaghi_strip_footing_bearing': {'pairs': 2450, 'decided': 2450, 'identity_failures': 0, 'constant_parts': 0},
    'template_time_rate_of_consolidation': {'pairs': 2450, 'decided': 2450, 'identity_failures': 0, 'constant_parts': 0},
    'template_time_to_phasor': {'pairs': 2450, 'decided': 2450, 'identity_failures': 0, 'constant_parts': 0},
    'template_truss_method_of_joints': {'pairs': 2450, 'decided': 2450, 'identity_failures': 0, 'constant_parts': 0},
    'template_truss_method_of_sections': {'pairs': 2450, 'false_accepts': 2, 'decided': 2450, 'identity_failures': 0, 'constant_parts': 0},
    'template_two_phase_specific_volume': {'pairs': 2450, 'decided': 2450, 'identity_failures': 0, 'constant_parts': 0},
    'template_two_state_steady_state': {'pairs': 2450, 'decided': 2450, 'identity_failures': 0, 'constant_parts': 0},
    'template_two_step_transition_probability': {'pairs': 2450, 'decided': 2450, 'identity_failures': 0, 'constant_parts': 0},
    'template_undamped_natural_frequency_torsional': {'pairs': 2450, 'decided': 2450, 'identity_failures': 0, 'constant_parts': 0},
    'template_undamped_natural_frequency_translational': {'pairs': 2450, 'decided': 2450, 'identity_failures': 0, 'constant_parts': 0},
    'template_undamped_response_initial_conditions': {'pairs': 2450, 'decided': 2450, 'identity_failures': 0, 'constant_parts': 0},
    'template_upward_seepage_quick_condition': {'pairs': 2450, 'decided': 2450, 'identity_failures': 0, 'constant_parts': 0},
    'template_utube_manometer': {'pairs': 2450, 'decided': 2450, 'identity_failures': 0, 'constant_parts': 0},
    'template_vdw_solve_for_pressure': {'pairs': 2450, 'false_accepts': 2, 'decided': 2450, 'identity_failures': 0, 'constant_parts': 0},
    'template_vdw_solve_for_volume': {'pairs': 2450, 'decided': 2450, 'identity_failures': 0, 'constant_parts': 2},
    'template_vibration_isolator_design': {'pairs': 2450, 'decided': 2450, 'identity_failures': 0, 'constant_parts': 0},
    'template_vibration_transmissibility': {'pairs': 2450, 'decided': 2450, 'identity_failures': 0, 'constant_parts': 0},
    'template_virtual_work_truss_deflection': {'pairs': 2450, 'decided': 2450, 'identity_failures': 0, 'constant_parts': 0},
    'template_volumetric_flow_rate': {'pairs': 2450, 'decided': 2450, 'identity_failures': 0, 'constant_parts': 0},
    'template_wave_equation_interpretation': {'pairs': 2450, 'decided': 2450, 'identity_failures': 0, 'constant_parts': 0},
    'template_wave_parameters_basic': {'pairs': 2450, 'decided': 2450, 'identity_failures': 0, 'constant_parts': 0},
    'template_work_isothermal_virial': {'pairs': 2450, 'decided': 2450, 'identity_failures': 0, 'constant_parts': 0},
    'template_xbar_known_sigma_classification': {'pairs': 2450, 'decided': 2450, 'identity_failures': 0, 'constant_parts': 0},
    'template_xbar_r_control_limits': {'pairs': 2450, 'decided': 2450, 'identity_failures': 0, 'constant_parts': 0},
}
# END GENERATED
