"""C1.3 - census of every numeric constants table, and the given-values measurement.

    python -m tests.constants_integrity.census                    # print the census
    python -m tests.constants_integrity.census --seeds 40         # probe depth
    python -m tests.constants_integrity.census --no-probe         # static part only (fast)
    python -m tests.constants_integrity.census --json out.json    # machine-readable
    python -m tests.constants_integrity.census --selftest         # planted defects

Phase C1.3 classifies every constants table as *needs a citation* or *needs only a
plausibility window*, and the spec calls it "the deliverable that sizes the rest
of the track". A classification is only as good as the measurement under it, and
the series' standing rule is that every count ships with its predicate or a
script that regenerates it. This is that script; the predicates are below, in
the words the counts use.

PREDICATES
----------
P-TABLE    A *numeric table* is an UPPER_CASE top-level assignment in a branch's
           constants.py - `ast.Assign` with exactly one `ast.Name` target
           matching ^[A-Z][A-Z0-9_]*$, or `ast.AnnAssign` - whose value contains
           at least one numeric literal.
P-LITERAL  A *numeric literal* is an `ast.Constant` whose value is an int or a
           float, bool excluded. `-1.5` is one literal (the minus is a UnaryOp).
           Dict KEYS count: `{0.90: 1.2816}` holds two literals.
P-TAGGED   A table is *tagged* if a provenance tag (TAG_RE: an opening bracket
           and a class word) appears on any source line of the assignment, or in
           the contiguous run of comment and blank lines immediately above it.
           This is the Prompt 07 predicate. A tag is a CLAIM, not evidence:
           resolution is test_citations_resolve's job, not this file's.
P-CONSUMER A template consumes a table if its `template_*` function - or any
           module-level function it calls, transitively, inside its own module -
           loads the table's name, or loads any of these ALIASES of it:
             * a local assigned from an expression that loads it
               (`all_fluids = {**COMMON_LIQUIDS, **COMMON_GASES}`);
             * a name bound at module level, in the template's module, by a
               statement that loads it (`_T26R_E24_5PCT = RESISTOR_SERIES_BY_
               TOLERANCE[5]`, or a module-level loop filling a list);
             * a name defined in constants.py - an assignment or a FUNCTION -
               whose body loads it (`chart_factor` reads CONTROL_CHART_FACTORS).
           Aliases chain to a fixpoint.
P-DRAW     A consumer *draws by order* if it calls `random.choice` or
           `random.sample` with an argument that loads the table (or an alias)
           and contains no `sorted(` call. For a dict the order that matters is
           KEY order; for a list or tuple, ELEMENT order (D-031).
P-GIVEN    The given-values measurement, by PERTURBATION rather than by string
           matching. For one FIELD of a table (a column: every leaf reached by
           the same path, with row keys and list indices wildcarded) every leaf
           is nudged - floats x1.0371 (a zero becomes 0.0371), ints +1 - the
           consuming modules are RELOADED so import-time copies see the nudge,
           and each consumer is re-run on the same seed. Per seed:
             RESTATED    the question changed, and it is identical to the
                         original once every number in it is masked: the field's
                         value is printed in the question, possibly rounded or
                         converted.
             HIDDEN      the question is byte-identical and the solution changed:
                         the value moves the gold answer without appearing in the
                         question. THIS is what needs a citation for correctness.
             STRUCTURAL  the question changed in more than its numbers (a
                         rejection loop took another path): indeterminate.
             UNUSED      neither changed.
           A field is HIDDEN for a consumer if ANY seed was HIDDEN. Masking is
           NUM_RE -> `#`. A module that fails to import under the nudge is
           IMPORT-ERROR: the table is validated at import time.

LIMITS, stated so that nobody reads more into a count than it holds
------------------------------------------------------------------
* A value printed in the question at a rounding coarse enough that a 3.71 %
  nudge does not move it reads HIDDEN. That is the conservative direction, and
  it is also the D-025 representability hazard, which C3.6 reports separately.
* Dict keys are never perturbed, so `CONTROL_CHART_FACTORS`' subgroup sizes and
  `Z_QUANTILES`' probabilities are not probed; they are labels a template prints.
* The probe sees one seed range. A field consumed only on a rare branch can read
  UNUSED; `--seeds` raises the depth, and the seed count is printed with every
  result.
* Consumption through `getattr`, `globals()` or a string name is invisible to
  the static part. None was found by grep when this was written.

CLASSIFICATION (C1.3) is computed from the measurement and the table's declared
`@kind`, never from the measurement alone - see `classify()` and the provenance
convention in docs/references/README.md.
"""
from __future__ import annotations

import argparse
import ast
import importlib
import json
import os
import random
import re
import sys

REPO = os.path.abspath(os.path.join(os.path.dirname(__file__), '..', '..'))
if REPO not in sys.path:
    sys.path.insert(0, REPO)

BRANCHES_DIR = os.path.join(REPO, 'data', 'templates', 'branches')
BRANCHES = ('chemical_engineering', 'electrical_engineering',
            'mechanical_engineering', 'civil_engineering',
            'industrial_engineering')

#: Every provenance tag either vocabulary has used - Prompt 07's list plus the
#: C1.2 local-only qualifier. Matching is on the OPENING bracket and the class
#: word only. The first version also required the closing `]` on the same line,
#: and so missed civil's CV_RANGES_M2_YR, whose `[POLICY: sampling-only. ...`
#: tag closes three comment lines later: it reported civil 20 tagged / 495
#: literals against the brief's 21 / 499. A detector narrower than the class it
#: is named after, found by reproducing someone else's count.
TAG_RE = re.compile(
    r'\[(ON-DISK|DERIVED|BY-DEFINITION|KNOWN-DEFECTIVE|VERIFY|POLICY|REALISM|'
    r'DERIVABLE|UNVERIFIED|LOCAL-ONLY)\b')
UPPER_RE = re.compile(r'^[A-Z][A-Z0-9_]*$')
NUM_RE = re.compile(r'[-+]?\d[\d,]*\.?\d*(?:[eE][-+]?\d+)?')

PERTURB_FACTOR = 1.0371


# ==========================================================================
# Static part: tables, literals, tags
# ==========================================================================

def _is_num(node):
    return (isinstance(node, ast.Constant) and not isinstance(node.value, bool)
            and isinstance(node.value, (int, float)))


def static_tables(src):
    """Apply P-TABLE, P-LITERAL and P-TAGGED to one constants.py source."""
    tree = ast.parse(src)
    lines = src.splitlines()
    out = []
    for node in tree.body:
        if isinstance(node, ast.Assign):
            if len(node.targets) != 1 or not isinstance(node.targets[0], ast.Name):
                continue
            name, value = node.targets[0].id, node.value
        elif isinstance(node, ast.AnnAssign) and isinstance(node.target, ast.Name):
            name, value = node.target.id, node.value
        else:
            continue
        if not UPPER_RE.match(name) or value is None:
            continue
        literals = sum(1 for n in ast.walk(value) if _is_num(n))
        if literals == 0:
            continue
        body = lines[node.lineno - 1:node.end_lineno]
        above = []
        i = node.lineno - 2
        while i >= 0 and (not lines[i].strip() or lines[i].lstrip().startswith('#')):
            above.append(lines[i])
            i -= 1
        above.reverse()
        tags = sorted({m.group(1) for ln in body + above for m in TAG_RE.finditer(ln)})
        out.append(dict(name=name, lineno=node.lineno, end_lineno=node.end_lineno,
                        literals=literals, tagged=bool(tags), tags=tags,
                        header='\n'.join(above)))
    return out


# ==========================================================================
# Static part: consumers and draws
# ==========================================================================

def _names_loaded(node):
    return {n.id for n in ast.walk(node) if isinstance(n, ast.Name)}


_SCOPES = (ast.ListComp, ast.SetComp, ast.DictComp, ast.GeneratorExp, ast.Lambda)


def _stores_outside_scopes(node):
    """Names bound by `node` in ITS OWN scope - not inside a comprehension or lambda."""
    out, stack = set(), [node]
    while stack:
        n = stack.pop()
        if isinstance(n, _SCOPES) and n is not node:
            continue
        if isinstance(n, ast.Name) and isinstance(n.ctx, ast.Store):
            out.add(n.id)
        stack.extend(ast.iter_child_nodes(n))
    return out


def _free_loads(fn):
    """Names a function loads that it does not bind itself - its GLOBALS.

    The first version of P-CONSUMER matched module aliases against every name a
    function loaded. A module-level `{k: 2 * v for k, v in LEVELS.items()}` then
    made `k` and `v` "aliases of LEVELS", and any template with its own local `k`
    became a consumer. The self-test's planted consumer set caught it: an
    over-broad detector, the other half of the shape Phase 5 kept meeting.
    """
    bound = {a.arg for a in ast.walk(fn) if isinstance(a, ast.arg)}
    bound |= {n.id for n in ast.walk(fn) if isinstance(n, ast.Name)
              and isinstance(n.ctx, ast.Store)}
    return _names_loaded(fn) - bound


def _module_functions(tree):
    return {n.name: n for n in tree.body if isinstance(n, ast.FunctionDef)}


def _reach(fn, funcs):
    """The function plus every module-level function it calls, transitively."""
    seen, stack = {}, [fn]
    while stack:
        f = stack.pop()
        if f.name in seen:
            continue
        seen[f.name] = f
        for n in ast.walk(f):
            if isinstance(n, ast.Call) and isinstance(n.func, ast.Name) \
                    and n.func.id in funcs:
                stack.append(funcs[n.func.id])
    return list(seen.values())


def _fixpoint(alias):
    """Chain aliases: if a is built from b and b from TABLE, a is built from TABLE."""
    changed = True
    while changed:
        changed = False
        for a, srcs in alias.items():
            extra = set().union(*(alias.get(s, set()) for s in srcs)) - srcs
            if extra:
                srcs |= extra
                changed = True
    return alias


def constants_aliases(src, tables):
    """Names in constants.py - assignments or functions - built from a table."""
    alias = {}
    for node in ast.parse(src).body:
        if isinstance(node, ast.Assign) and len(node.targets) == 1 \
                and isinstance(node.targets[0], ast.Name):
            name, body = node.targets[0].id, node.value
        elif isinstance(node, ast.FunctionDef):
            name, body = node.name, node
        else:
            continue
        loads = _names_loaded(body) - {name}
        if loads:
            alias[name] = set(loads)
    alias = _fixpoint(alias)
    return {a: s & set(tables) for a, s in alias.items() if s & set(tables)}


def _module_aliases(tree, tables, calias):
    """Names bound at module level, in a template module, by a table-loading statement."""
    skip = (ast.FunctionDef, ast.ClassDef, ast.Import, ast.ImportFrom)
    stmts = [st for st in tree.body if not isinstance(st, skip)]
    bound = set().union(*(_stores_outside_scopes(st) for st in stmts)) if stmts else set()
    alias = {}
    for st in stmts:
        loaded = _names_loaded(st)
        for nm in (loaded | _stores_outside_scopes(st)) & bound:
            alias.setdefault(nm, set()).update(loaded - {nm})
    alias = _fixpoint(alias)
    known = set(tables) | set(calias)
    out = {}
    for a, srcs in alias.items():
        hit = (srcs & set(tables)) | set().union(*(calias[c] for c in srcs & set(calias)))
        if hit and a not in known:
            out[a] = hit
    return out


def _local_aliases(fns, names):
    """Locals assigned from an expression that loads one of `names`."""
    alias = {}
    for f in fns:
        for n in ast.walk(f):
            if isinstance(n, (ast.Assign, ast.For, ast.comprehension)):
                value = n.value if isinstance(n, ast.Assign) else n.iter
                targets = n.targets if isinstance(n, ast.Assign) else [n.target]
                src = _names_loaded(value)
                for tgt in targets:
                    for nm in _names_loaded(tgt):
                        alias.setdefault(nm, set()).update(src - {nm})
    alias = _fixpoint(alias)
    return {a: s & set(names) for a, s in alias.items() if s & set(names)}


def consumers_and_draws(branch_dir, module_root, tables, calias):
    """Apply P-CONSUMER and P-DRAW to every template module under `branch_dir`."""
    use = {t: {} for t in tables}
    for dirpath, _d, files in os.walk(branch_dir):
        if '__pycache__' in dirpath:
            continue
        for fn in sorted(files):
            if not fn.endswith('.py') or fn in ('constants.py', '__init__.py'):
                continue
            path = os.path.join(dirpath, fn)
            modname = os.path.relpath(path, module_root)[:-3].replace(os.sep, '.')
            tree = ast.parse(open(path, encoding='utf-8').read())
            funcs = _module_functions(tree)
            # table -> every name that stands for it in this module
            here = {a: set(s) for a, s in calias.items()}
            for a, s in _module_aliases(tree, tables, calias).items():
                here.setdefault(a, set()).update(s)
            for name, f in funcs.items():
                if not name.startswith('template_'):
                    continue
                reach = _reach(f, funcs)
                loaded = set().union(*(_names_loaded(r) for r in reach))
                free = set().union(*(_free_loads(r) for r in reach))
                local = _local_aliases(reach, (set(tables) | set(here)) & free)
                for a, s in list(local.items()):
                    local[a] = (s & set(tables)) | set().union(
                        *(here[x] for x in s & set(here)))
                for t in tables:
                    globals_for_t = {t} | {a for a, s in here.items() if t in s}
                    locals_for_t = {a for a, s in local.items() if t in s}
                    stand_ins = globals_for_t | locals_for_t
                    # a table or a module/constants alias counts only as a
                    # GLOBAL the function reads, never as a local it shadows
                    present = (globals_for_t & free) | (locals_for_t & loaded)
                    via = sorted(('direct' if x == t else f'alias {x}') for x in present)
                    if not via:
                        continue
                    draws = set()
                    for r in reach:
                        for n in ast.walk(r):
                            if not (isinstance(n, ast.Call) and isinstance(n.func, ast.Attribute)
                                    and isinstance(n.func.value, ast.Name)
                                    and n.func.value.id == 'random'
                                    and n.func.attr in ('choice', 'sample') and n.args):
                                continue
                            arg = n.args[0]
                            if not (_names_loaded(arg) & stand_ins):
                                continue
                            is_sorted = any(isinstance(c, ast.Call)
                                            and isinstance(c.func, ast.Name)
                                            and c.func.id == 'sorted' for c in ast.walk(arg))
                            draws.add('sorted' if is_sorted else f'random.{n.func.attr}')
                    use[t][name] = dict(via=via, draws=sorted(draws), module=modname)
    return use


# ==========================================================================
# Dynamic part: the P-GIVEN perturbation probe
# ==========================================================================

def _leaves(obj, path=()):
    if isinstance(obj, bool):
        return
    if isinstance(obj, (int, float)):
        yield path, obj
    elif isinstance(obj, dict):
        for k, v in obj.items():
            yield from _leaves(v, path + (('d', k),))
    elif isinstance(obj, list):
        for i, v in enumerate(obj):
            yield from _leaves(v, path + (('l', i),))
    elif isinstance(obj, tuple):
        for i, v in enumerate(obj):
            yield from _leaves(v, path + (('t', i),))


def field_of(path):
    """Wildcard row keys and list indices; keep tuple indices and string subkeys."""
    parts = []
    for depth, (kind, key) in enumerate(path):
        if kind == 'l' or (kind == 'd' and (depth == 0 or not isinstance(key, str))):
            parts.append('[*]')
        elif kind == 't':
            parts.append(f'[{key}]')
        else:
            parts.append(f'[{key!r}]')
    return ''.join(parts)


def fields(obj):
    out = {}
    for path, _v in _leaves(obj):
        out[field_of(path)] = out.get(field_of(path), 0) + 1
    return out


def _nudge(v):
    # Integers are nudged RELATIVELY too, kept integral. The first version added
    # 1, which moved C0 = 299792458 by 3e-9 - below every display precision - so
    # the probe reported C0 as having NO EFFECT on two templates that consume it.
    if isinstance(v, int):
        return v + max(1, int(round(abs(v) * (PERTURB_FACTOR - 1))))
    return v * PERTURB_FACTOR if v != 0 else 0.0371


def perturbed(obj, field, path=()):
    if isinstance(obj, bool):
        return obj
    if isinstance(obj, (int, float)):
        return _nudge(obj) if field_of(path) == field else obj
    if isinstance(obj, dict):
        return {k: perturbed(v, field, path + (('d', k),)) for k, v in obj.items()}
    if isinstance(obj, list):
        return [perturbed(v, field, path + (('l', i),)) for i, v in enumerate(obj)]
    if isinstance(obj, tuple):
        return tuple(perturbed(v, field, path + (('t', i),)) for i, v in enumerate(obj))
    return obj


class Patch:
    """Swap a table for its nudged copy inside the constants module.

    The constants module holds the table, and may hold a container that holds it
    (`RESISTOR_SERIES_BY_TOLERANCE`); both are swapped, identity-matched, and
    restored on exit. Template modules are NOT patched here - they are reloaded
    by `probe`, which is what reaches a value they copied at import time.
    """

    def __init__(self, cmod, original, replacement):
        self.cmod, self.o, self.r, self.undo = cmod, original, replacement, []

    def __enter__(self):
        for attr, val in list(vars(self.cmod).items()):
            if val is self.o:
                new = self.r
            elif isinstance(val, dict) and any(x is self.o for x in val.values()):
                new = {k: (self.r if x is self.o else x) for k, x in val.items()}
            elif isinstance(val, (list, tuple)) and any(x is self.o for x in val):
                new = type(val)(self.r if x is self.o else x for x in val)
            else:
                continue
            self.undo.append((attr, val))
            setattr(self.cmod, attr, new)
        return self

    def __exit__(self, *exc):
        for attr, val in reversed(self.undo):
            setattr(self.cmod, attr, val)
        return False


def _gen(fn, seed):
    random.seed(seed)
    try:
        import numpy as np
        np.random.seed(seed)
    except ImportError:
        pass
    try:
        q, s = fn()
        return str(q), str(s)
    except Exception as exc:                                    # noqa: BLE001
        return None, f'ERROR {type(exc).__name__}: {exc}'


def _reload_all(modnames):
    mods, err = {}, None
    for m in modnames:
        try:
            mods[m] = importlib.reload(sys.modules[m]) if m in sys.modules \
                else importlib.import_module(m)
        except Exception as exc:                                # noqa: BLE001
            err = f'{m}: {type(exc).__name__}: {exc}'
    return mods, err


def probe(cmod, table_name, field, consumers, seeds, baseline):
    """P-GIVEN for one field against each consumer. Returns tid -> counts."""
    orig = getattr(cmod, table_name)
    modnames = sorted(set(consumers.values()))
    out = {}
    try:
        with Patch(cmod, orig, perturbed(orig, field)):
            mods, err = _reload_all(modnames)
            for tid, m in consumers.items():
                c = dict(RESTATED=0, HIDDEN=0, STRUCTURAL=0, UNUSED=0, ERROR=0,
                         IMPORT_ERROR=0)
                if err or m not in mods:
                    c['IMPORT_ERROR'] = seeds
                    out[tid] = c
                    continue
                fn = getattr(mods[m], tid)
                for s in range(seeds):
                    q0, s0 = baseline[tid][s]
                    q1, s1 = _gen(fn, s)
                    if q0 is None or q1 is None:
                        c['ERROR' if (q0 is None) != (q1 is None) else 'UNUSED'] += 1
                        continue
                    if q0 == q1:
                        c['HIDDEN' if s0 != s1 else 'UNUSED'] += 1
                    elif NUM_RE.sub('#', q0) == NUM_RE.sub('#', q1):
                        c['RESTATED'] += 1
                    else:
                        c['STRUCTURAL'] += 1
                out[tid] = c
    finally:
        _reload_all(modnames)          # against the restored table
    return out


def field_verdict(c):
    for v in ('HIDDEN', 'RESTATED', 'STRUCTURAL', 'IMPORT_ERROR', 'ERROR'):
        if c.get(v):
            return v.replace('_', '-')
    return 'UNUSED'


# ==========================================================================
# The census
# ==========================================================================

HEADER_FIELD_RE = re.compile(r'#\s*@(units|domain|kind)\s*:\s*(.*)')


def header_fields(header):
    out = {}
    for ln in header.splitlines():
        m = HEADER_FIELD_RE.match(ln.strip())
        if m:
            out[m.group(1)] = (out.get(m.group(1), '') + ' ' + m.group(2)).strip()
    return out


def classify(kind, measured):
    """C1.3's rule, from the declared @kind and the P-GIVEN measurement.

    The given-values rule excuses a value from CORRECTNESS-dependence: restated
    in the question, a wrong number cannot make the gold answer disagree with
    the question. It does not excuse a value from being TRUE. A question that
    says mercury has a density of 5000 kg/m^3 is self-consistent and false, and
    Phase C2 shipped 274 items of exactly that (phaseC2_summary.md §6). So the
    rule licenses a plausibility window only for a `range` - a sampling window
    that asserts nothing about a named real entity - and only when every
    consumer is measured to restate it.
    """
    if kind is None:
        return 'UNDECLARED', 'no @kind declared'
    kind = kind.split()[0]
    if kind not in KINDS:
        return 'UNDECLARED', f'unknown @kind {kind!r}'
    if measured == 'UNCONSUMED':
        return 'UNCONSUMED', 'no template consumes it'
    if kind == 'range':
        if measured == 'ALL-RESTATED':
            return 'PLAUSIBILITY', 'a sampling window every consumer restates'
        if measured in ('NO-EFFECT', 'INDETERMINATE'):
            # Consumed as a GUARD: an assert or a screen over values drawn
            # elsewhere. A nudge cannot move an item through a guard, only
            # reject it; the window shapes plausibility, not correctness.
            return 'PLAUSIBILITY', f'{measured}: consumed as a guard or screen'
        return 'CITATION', (f'declared a range but measured {measured}: the '
                            f'given-values rule does not apply')
    if kind in ('property', 'standard', 'measured-constant'):
        return 'CITATION', ('consumed without being stated' if measured == 'SOME-HIDDEN'
                            else f'{measured}: asserts a fact about a named real entity')
    if kind == 'mathematical':
        return 'DERIVATION', 'computable from its definition'
    if kind == 'defined':
        return 'DEFINITION', 'exact by definition or convention'
    return 'DOMAIN', 'a validity bound for another table; checked by its C3.2 suite'


#: The declared kinds, one meaning each. See docs/references/README.md.
KINDS = ('property', 'standard', 'measured-constant', 'defined', 'mathematical',
         'range', 'validity')

#: Probe limits found while measuring the real corpus, recorded where the
#: numbers are produced rather than in a summary that can drift from them:
#:
#:  * GUARD CONSUMPTION. Many civil/industrial templates read a table only in an
#:    assert or a screen, or enumerate candidates over it and then filter. A
#:    nudge then reads ERROR (the assert fires), NO-EFFECT (the filter absorbs
#:    it) or STRUCTURAL. `classify` treats that as a guard for a `range`; for a
#:    `property` or `standard` it changes nothing, because a guard that asserts
#:    a template-side COPY equals the table (hydrology.py's `_SCS_COMBOS`) makes
#:    the table's truth the copy's truth.
#:  * PARALLEL STRINGS. REACTIONS carries each reaction twice: numeric
#:    coefficients, and an `equation` string the question prints. The probe
#:    nudges only numbers, so the coefficients read HIDDEN while the question
#:    does state them. Their warrant is element balance, which is why the table
#:    is `mathematical`, not a measurement this probe could settle.
#:  * CONSTANTS-LEVEL COPIES. MEDIA_VELOCITIES["Vacuum"] is C0 copied at import
#:    of constants.py itself, which the probe does not reload. C0's consumption
#:    through it is seen statically (P-CONSUMER) and not dynamically.


def census(seeds=40, probe_on=True, branches=BRANCHES, root=BRANCHES_DIR,
           module_root=REPO, package='data.templates.branches'):
    report = []
    for branch in branches:
        src = open(os.path.join(root, branch, 'constants.py'), encoding='utf-8').read()
        tables = static_tables(src)
        names = sorted(t['name'] for t in tables)
        calias = constants_aliases(src, names)
        use = consumers_and_draws(os.path.join(root, branch), module_root, names, calias)
        cmod = importlib.import_module(f'{package}.{branch}.constants' if package
                                       else f'{branch}.constants')
        baseline = {}
        for t in tables:
            obj = getattr(cmod, t['name'])
            t.update(branch=branch, type=type(obj).__name__, fields=fields(obj),
                     consumers=use[t['name']], declared=header_fields(t.pop('header')),
                     probe={}, seeds=seeds if probe_on else 0)
            t['draws_by_order'] = sorted(tid for tid, u in t['consumers'].items()
                                         if any(d.startswith('random.') for d in u['draws']))
            if probe_on and t['consumers']:
                cons = {tid: u['module'] for tid, u in t['consumers'].items()}
                for tid, m in cons.items():
                    if tid not in baseline:
                        fn = getattr(importlib.import_module(m), tid)
                        baseline[tid] = [_gen(fn, s) for s in range(seeds)]
                for fld in t['fields']:
                    res = probe(cmod, t['name'], fld, cons, seeds, baseline)
                    t['probe'][fld] = {tid: dict(counts=c, verdict=field_verdict(c))
                                       for tid, c in res.items()}
            verdicts = [v['verdict'] for f in t['probe'].values() for v in f.values()]
            if not t['consumers']:
                t['measured'] = 'UNCONSUMED'
            elif not probe_on:
                t['measured'] = 'NOT-PROBED'
            elif 'HIDDEN' in verdicts:
                t['measured'] = 'SOME-HIDDEN'
            elif 'RESTATED' in verdicts and all(v in ('RESTATED', 'UNUSED') for v in verdicts):
                t['measured'] = 'ALL-RESTATED'
            elif 'STRUCTURAL' in verdicts or 'IMPORT-ERROR' in verdicts:
                t['measured'] = 'INDETERMINATE'
            else:
                t['measured'] = 'NO-EFFECT'
            t['class'], t['class_reason'] = classify(t['declared'].get('kind'), t['measured'])
            report.append(t)
    return report


def print_report(report):
    print('PROVENANCE COVERAGE (P-TABLE, P-LITERAL, P-TAGGED)')
    print(f"{'branch':24s} {'tables':>6s} {'tagged':>6s} {'literals':>8s} {'in tagged':>9s}")
    for b in BRANCHES:
        rows = [t for t in report if t['branch'] == b]
        if rows:
            print(f"{b:24s} {len(rows):6d} {sum(t['tagged'] for t in rows):6d} "
                  f"{sum(t['literals'] for t in rows):8d} "
                  f"{sum(t['literals'] for t in rows if t['tagged']):9d}")
    print()
    print('PER TABLE - measured is P-GIVEN over every consumer and field, '
          f"at {report[0]['seeds'] if report else 0} seeds per template")
    for t in report:
        per = {}
        for f in t['probe'].values():
            for tid, v in f.items():
                per.setdefault(v['verdict'], set()).add(tid)
        print(f"  {t['branch'][:4]} {t['name']:30s} lit={t['literals']:4d} "
              f"tag={','.join(t['tags']) or '-':20s} cons={len(t['consumers']):2d} "
              f"order={len(t['draws_by_order']):2d} {t['measured']:13s} "
              f"units={'Y' if 'units' in t['declared'] else '-'} "
              f"kind={t['declared'].get('kind', '-').split(' ')[0]:17s} {t['class']}")
        for verdict in ('HIDDEN', 'STRUCTURAL', 'IMPORT-ERROR', 'ERROR'):
            if verdict in per:
                print(f"       {verdict.lower()} in: {', '.join(sorted(per[verdict]))}")
    print()
    counts = {}
    for t in report:
        counts[t['class']] = counts.get(t['class'], 0) + 1
    print('CLASSIFICATION', ', '.join(f'{k} {v}' for k, v in sorted(counts.items())))


# ==========================================================================
# Self-test: planted defects, written from the predicate definitions above,
# and run through the SAME functions the census uses - never a copy of them.
# ==========================================================================

_SELFTEST_CONSTANTS = '''
# [ON-DISK] somewhere
DENSITY = {"Alpha": 1000.0, "Beta": 850.0, "Gamma": 13550.0}

PAIRS = {"Alpha": (1000.0, 1.0e-3), "Beta": (850.0, 2.0e-3)}

# [POLICY: sampling-only, and this tag's bracket
# closes on a later line]

WINDOW = (10.0, 20.0)

# [ON-DISK] far away - does not reach ORPHAN, because code intervenes
UNRELATED = 1
ORPHAN = {"x": 5.0}
NAMES = ["a", "b"]
LEVELS = {"low": 1.5, "high": 2.5}
SPREAD = {"x": 4.0, "y": 9.0}


def density_of(name):
    return DENSITY[name]
'''

_SELFTEST_TEMPLATES = '''
import random
from selftest_branch.constants import (DENSITY, PAIRS, WINDOW, ORPHAN, LEVELS,
                                       SPREAD, density_of)

_COPIED = {k: 2 * v for k, v in LEVELS.items()}

def _helper(rows):
    return {k: v for k, v in rows.items()}

def template_restates():
    name, rho = random.choice(list(DENSITY.items()))
    return f"{name} has density {rho} kg/m3. Find 2 rho.", f"**Answer:** {2 * rho}"

def template_hides_other_field_of_restated_row():
    local = _helper(PAIRS)
    name, (rho, mu) = random.choice(list(local.items()))
    return f"{name} at rho = {rho}. Find mu/rho.", f"**Answer:** {mu / rho:.6e}"

def template_window():
    lo, hi = WINDOW
    x = round(random.uniform(lo, hi), 1)
    return f"x = {x}. Find x.", f"**Answer:** {x}"

def template_import_time_copy():
    k = random.choice(sorted(_COPIED))
    return f"Level {k}.", f"**Answer:** {_COPIED[k]}"

def template_constants_accessor():
    name = random.choice(["Alpha", "Beta"])
    return f"A block of {name}.", f"**Answer:** {3 * density_of(name)}"

def template_loop_alias():
    total = 0.0
    for key in SPREAD:
        total += SPREAD[key]
    return "Sum the spread.", f"**Answer:** {total}"
'''


def selftest():
    import tempfile
    failures = []
    tmp = tempfile.mkdtemp(prefix='census_selftest_')
    pkg = os.path.join(tmp, 'selftest_branch')
    os.makedirs(pkg)
    open(os.path.join(pkg, '__init__.py'), 'w').close()
    open(os.path.join(pkg, 'constants.py'), 'w', encoding='utf-8').write(_SELFTEST_CONSTANTS)
    open(os.path.join(pkg, 'tmpl.py'), 'w', encoding='utf-8').write(_SELFTEST_TEMPLATES)
    sys.path.insert(0, tmp)
    try:
        # -- P-TABLE / P-LITERAL / P-TAGGED ----------------------------------
        st = {t['name']: t for t in static_tables(_SELFTEST_CONSTANTS)}
        for n, want in {'DENSITY': True, 'PAIRS': False, 'WINDOW': True,
                        'ORPHAN': False, 'UNRELATED': True}.items():
            if st[n]['tagged'] != want:
                failures.append(f'P-TAGGED {n}: got {st[n]["tagged"]}, planted {want}')
        if 'NAMES' in st:
            failures.append('P-TABLE counted a table with no numeric literal')
        if st['PAIRS']['literals'] != 4:
            failures.append(f'P-LITERAL PAIRS: {st["PAIRS"]["literals"]} != 4')

        # -- P-CONSUMER / P-DRAW, through the census's own code path ----------
        rep = {t['name']: t for t in census(seeds=12, branches=('selftest_branch',),
                                            root=tmp, module_root=tmp, package='')}
        want_cons = {
            'DENSITY': {'template_restates', 'template_constants_accessor'},
            'PAIRS': {'template_hides_other_field_of_restated_row'},
            'WINDOW': {'template_window'},
            'LEVELS': {'template_import_time_copy'},
            'SPREAD': {'template_loop_alias'},
            'ORPHAN': set(),
        }
        for tbl, want in want_cons.items():
            got = set(rep[tbl]['consumers'])
            if got != want:
                failures.append(f'P-CONSUMER {tbl}: got {sorted(got)}, planted {sorted(want)}')
        want_draw = {('DENSITY', 'template_restates'): ['random.choice'],
                     ('PAIRS', 'template_hides_other_field_of_restated_row'): ['random.choice'],
                     ('LEVELS', 'template_import_time_copy'): ['sorted'],
                     ('DENSITY', 'template_constants_accessor'): []}
        for (tbl, tid), want in want_draw.items():
            got = rep[tbl]['consumers'].get(tid, {}).get('draws')
            if got != want:
                failures.append(f'P-DRAW {tbl} in {tid}: got {got}, planted {want}')

        # -- P-GIVEN ----------------------------------------------------------
        want_given = {
            ('DENSITY', '[*]', 'template_restates'): 'RESTATED',
            ('DENSITY', '[*]', 'template_constants_accessor'): 'HIDDEN',
            ('PAIRS', '[*][0]', 'template_hides_other_field_of_restated_row'): 'RESTATED',
            ('PAIRS', '[*][1]', 'template_hides_other_field_of_restated_row'): 'HIDDEN',
            ('WINDOW', '[0]', 'template_window'): 'RESTATED',
            ('LEVELS', '[*]', 'template_import_time_copy'): 'HIDDEN',
            ('SPREAD', '[*]', 'template_loop_alias'): 'HIDDEN',
        }
        for (tbl, fld, tid), want in want_given.items():
            got = rep[tbl]['probe'].get(fld, {}).get(tid, {}).get('verdict')
            if got != want:
                failures.append(f'P-GIVEN {tbl}{fld} in {tid}: got {got}, planted {want}')
        if rep['DENSITY']['measured'] != 'SOME-HIDDEN':
            failures.append(f"table verdict DENSITY: {rep['DENSITY']['measured']}")

        # -- the probe must leave everything as it found it ------------------
        import selftest_branch.constants as sc
        import selftest_branch.tmpl as stt
        if sc.DENSITY['Gamma'] != 13550.0 or stt._COPIED != {'low': 3.0, 'high': 5.0}:
            failures.append('the probe did not restore the tables or the reloaded module')
    finally:
        sys.path.remove(tmp)
        for m in [m for m in sys.modules if m.startswith('selftest_branch')]:
            del sys.modules[m]

    for f in failures:
        print('  - ' + f)
    print(f'selftest: {len(failures)} failure(s)')
    return 1 if failures else 0


def write_markdown(report, path, seeds):
    """The committed C1.3 report. Regenerated, never hand-edited."""
    out = [
        '# Phase C1.3 — constants-table census and classification',
        '',
        '**Generated, not written.** Regenerate with',
        '',
        f'    python -m tests.constants_integrity.census --seeds {seeds} --markdown '
        'docs/re-implementation-sep/phaseC1_census.md',
        '',
        'Every term below - *numeric table*, *literal*, *tagged*, *consumer*, *draws by',
        'order*, *restated*, *hidden* - is a predicate defined in the docstring of',
        '`tests/constants_integrity/census.py`; the classification rule is spec §C1.3.',
        '',
        '## Coverage',
        '',
        '| branch | numeric tables | tagged | numeric literals | in tagged tables |',
        '|---|---:|---:|---:|---:|',
    ]
    for b in BRANCHES:
        rows = [t for t in report if t['branch'] == b]
        out.append(f"| {b} | {len(rows)} | {sum(t['tagged'] for t in rows)} | "
                   f"{sum(t['literals'] for t in rows)} | "
                   f"{sum(t['literals'] for t in rows if t['tagged'])} |")
    out.append(f"| **all** | **{len(report)}** | **{sum(t['tagged'] for t in report)}** | "
               f"**{sum(t['literals'] for t in report)}** | "
               f"**{sum(t['literals'] for t in report if t['tagged'])}** |")
    out += ['', '## Classification', '',
            '| class | tables | literals |', '|---|---:|---:|']
    for cls in ('CITATION', 'PLAUSIBILITY', 'DERIVATION', 'DEFINITION', 'DOMAIN',
                'UNCONSUMED', 'UNDECLARED'):
        rows = [t for t in report if t['class'] == cls]
        if rows:
            out.append(f"| {cls} | {len(rows)} | {sum(t['literals'] for t in rows)} |")
    out += ['', f'## Per table ({seeds} seeds per consuming template)', '',
            '| branch | table | literals | tags | consumers | draw by order | P-GIVEN | '
            '@kind | class | why | hidden in |',
            '|---|---|---:|---|---:|---:|---|---|---|---|---|']
    for t in report:
        hidden = sorted({tid.replace('template_', '') for f in t['probe'].values()
                         for tid, v in f.items() if v['verdict'] == 'HIDDEN'})
        out.append(f"| {t['branch'].split('_')[0]} | `{t['name']}` | {t['literals']} | "
                   f"{', '.join(t['tags']) or '-'} | {len(t['consumers'])} | "
                   f"{len(t['draws_by_order'])} | {t['measured']} | "
                   f"{t['declared'].get('kind', '-').split(' ')[0]} | **{t['class']}** | "
                   f"{t['class_reason']} | {', '.join(hidden) or '-'} |")
    with open(path, 'w', encoding='utf-8', newline='\n') as fh:
        fh.write('\n'.join(out) + '\n')


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument('--seeds', type=int, default=40)
    ap.add_argument('--no-probe', action='store_true')
    ap.add_argument('--json')
    ap.add_argument('--markdown')
    ap.add_argument('--branch', action='append')
    ap.add_argument('--selftest', action='store_true')
    args = ap.parse_args(argv)
    if args.selftest:
        return selftest()
    report = census(args.seeds, not args.no_probe, tuple(args.branch or BRANCHES))
    print_report(report)
    if args.json:
        with open(args.json, 'w', encoding='utf-8') as fh:
            json.dump(report, fh, indent=1, default=str)
    if args.markdown:
        write_markdown(report, args.markdown, args.seeds)
    return 0


if __name__ == '__main__':
    sys.exit(main())
