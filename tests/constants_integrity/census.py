"""C1.3 - census of every numeric constants table, and the given-values measurement.

    python -m tests.constants_integrity.census                    # print the census
    python -m tests.constants_integrity.census --seeds 40         # probe depth
    python -m tests.constants_integrity.census --no-probe         # static part only (fast)
    python -m tests.constants_integrity.census --json out.json    # machine-readable
    python -m tests.constants_integrity.census --markdown out.md  # the committed report
    python -m tests.constants_integrity.census --check            # exit 1 on REVIEW/UNDECLARED
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
           at least one numeric literal. This is the Prompt 07 predicate, and the
           COVERAGE table counts only these.
P-TABLE-LIVE
           Such an assignment whose value holds a number but no numeric literal
           - `(-math.pi, math.pi)`. Not counted in coverage (so the brief's figures
           still reproduce); CLASSIFIED and held to @kind/@units like any other.
           Added after C1 Reviewer G, finding G-4.
P-LITERAL  A *numeric literal* is an `ast.Constant` whose value is an int or a
           float, bool excluded. `-1.5` is one literal (the minus is a UnaryOp).
           Dict KEYS count: `{0.90: 1.2816}` holds two literals.
P-TAGGED   A table is *tagged* if a provenance tag (TAG_RE: an opening bracket
           and a class word) appears on any source line of the assignment, or in
           the contiguous run of comment and blank lines immediately above it.
           A tag is a CLAIM, not evidence: resolution is test_citations_resolve's
           job, not this file's.
P-CONSUMER A template consumes a table if its `template_*` function - or any
           module-level function it calls, transitively, inside its own module -
           reads, as a global, the table's name or a name that STANDS FOR it:
             * a name defined in constants.py - an assignment or a FUNCTION -
               whose body reads it (`chart_factor` reads CONTROL_CHART_FACTORS);
             * a name bound at module level in the template's module whose value,
               at the END of the module's import, derives from it.
           An expression derives from a table through the names it loads AND
           through every call it makes to a module-level function that reads the
           table (transitively). The first version stopped at the called name, so
           `_T27_PLANS = _t27_admissible_plans()` derived from nothing and the
           p-chart template's `random.choice(_T27_PLANS)` was an invisible draw
           (found by the C3.10 guard detector calling two tables guard-only there;
           a corpus-wide scan found that one site).
           Derivation is tracked FLOW-SENSITIVELY, statement by statement: an
           assignment replaces what a name stands for, a container fill
           (`X.append(...)`, `X[k] = ...`) adds to it, a loop is walked to a
           fixpoint and merged with the path that skips it, and an `if` merges
           both branches. The first two versions were flow-insensitive, and a
           helper name like `_lo` bound from a table in one module-level loop and
           rebound from unrelated data in the next "stood for" the table in both
           (C1 Reviewer G §5; `QUEUE_SCENARIOS` was credited with 4 consumers
           where 2 read it).
P-COPY     A template that writes a table's value as a number instead of reading
           it consumes the table invisibly to P-CONSUMER. A copy found is DECLARED
           in the table's header, `# @copied-in: <template_id> <literal>`, and the
           census VERIFIES the literal is a constant in that template's source. A
           declaration that no longer verifies is an error (G-1).
P-DRAW     A consumer *draws by order* if it calls `random.choice` or
           `random.sample` with an argument that derives from the table and
           contains no `sorted(` call. For a dict the order that matters is KEY
           order; for a list or tuple, ELEMENT order (D-031).
P-GIVEN    The given-values measurement, by PERTURBATION rather than by string
           matching. For one FIELD of a table (a column: every leaf reached by
           the same path, with row keys and list indices wildcarded) every leaf
           is nudged - floats x1.0371 (a zero becomes 0.0371), integers by the
           same relative step, kept integral - the consuming modules are RELOADED
           so import-time copies see the nudge, and each consumer is re-run on
           the same seed. Per seed:
             RESTATED    the question changed, and it is identical to the
                         original once every number in it is masked.
             HIDDEN      the question is byte-identical and the solution changed.
             STRUCTURAL  the question changed in more than its numbers.
             ERROR       the consumer raised under the nudge (or its module no
                         longer imports). The table is READ; how, the probe
                         cannot say.
             UNUSED      neither changed.
           Rolled up per table, in this order: SOME-HIDDEN if any HIDDEN;
           ALL-RESTATED if only RESTATED/UNUSED; ERROR if any ERROR; INDETERMINATE
           if any STRUCTURAL; COPIED if the only consumers are declared copies;
           else NO-EFFECT. The first version had no ERROR case, so an all-crash
           table rolled up as NO-EFFECT and was granted PLAUSIBILITY "as a guard"
           (C1 Reviewer G, G-2).

LIMITS, stated so that nobody reads more into a count than it holds
------------------------------------------------------------------
* A value printed at a rounding coarse enough that a 3.71 % nudge does not move
  it reads HIDDEN - the conservative direction, and D-025's hazard (C3.6).
* Dict keys are never perturbed; they are labels a template prints.
* The probe sees one seed range; `--seeds` raises the depth.
* Consumption through `getattr`, `globals()` or a string name is invisible. C1
  Reviewer G's sweep found none.
* An UNDECLARED literal copy is invisible. P-COPY checks the copies declared;
  finding the rest is C3's literal-copy sweep, so the UNCONSUMED count is an
  upper bound.
* GUARD CONSUMPTION reads NO-EFFECT, and only NO-EFFECT - nothing changed and
  nothing raised - is accepted as a guard without a human verdict.
* PARALLEL STRINGS: REACTIONS restates its coefficients through an `equation`
  string the probe does not nudge, so the coefficients read HIDDEN.
* CONSTANTS-LEVEL COPIES (MEDIA_VELOCITIES["Vacuum"] from C0) are seen statically
  and not dynamically: constants.py itself is not reloaded.

CLASSIFICATION (C1.3) is computed from the measurement, the table's declared
`@kind`, and - for a `range` the probe cannot settle - a declared, evidenced
`@given`; never from the measurement alone. See `classify()` and spec §C1.3.
"""
from __future__ import annotations

import argparse
import ast
import copy
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

#: The declared kinds, one meaning each. See spec §C1.2 and §C1.3.
KINDS = ('property', 'standard', 'measured-constant', 'defined', 'mathematical',
         'range', 'validity')


# ==========================================================================
# Static part: tables, literals, tags
# ==========================================================================

def _is_num(node):
    return (isinstance(node, ast.Constant) and not isinstance(node.value, bool)
            and isinstance(node.value, (int, float)))


def _assignments(tree):
    for node in tree.body:
        if isinstance(node, ast.Assign):
            if len(node.targets) != 1 or not isinstance(node.targets[0], ast.Name):
                continue
            name, value = node.targets[0].id, node.value
        elif isinstance(node, ast.AnnAssign) and isinstance(node.target, ast.Name):
            name, value = node.target.id, node.value
        else:
            continue
        if UPPER_RE.match(name) and value is not None:
            yield node, name, value


def _entry(node, name, literals, lines):
    body = lines[node.lineno - 1:node.end_lineno]
    above = []
    i = node.lineno - 2
    while i >= 0 and (not lines[i].strip() or lines[i].lstrip().startswith('#')):
        above.append(lines[i])
        i -= 1
    above.reverse()
    tags = sorted({m.group(1) for ln in body + above for m in TAG_RE.finditer(ln)})
    return dict(name=name, lineno=node.lineno, end_lineno=node.end_lineno,
                literals=literals, tagged=bool(tags), tags=tags,
                header='\n'.join(above), live_only=False)


def static_tables(src):
    """Apply P-TABLE, P-LITERAL and P-TAGGED to one constants.py source."""
    tree = ast.parse(src)
    lines = src.splitlines()
    out = []
    for node, name, value in _assignments(tree):
        literals = sum(1 for n in ast.walk(value) if _is_num(n))
        if literals:
            out.append(_entry(node, name, literals, lines))
    return out


def numeric_tables(src, ns=None):
    """P-TABLE plus P-TABLE-LIVE: every UPPER_CASE table whose VALUE holds a number.

    `ns` is the executed namespace of `src`; it is built here when not given.
    """
    if ns is None:
        ns = {}
        exec(compile(src, '<constants>', 'exec'), ns)            # noqa: S102
    tables = static_tables(src)
    have = {t['name'] for t in tables}
    lines = src.splitlines()
    for node, name, _value in _assignments(ast.parse(src)):
        obj = ns.get(name)
        if name in have or obj is None or callable(obj):
            continue
        if any(True for _ in _leaves(obj)):
            e = _entry(node, name, 0, lines)
            e['live_only'] = True
            tables.append(e)
    tables.sort(key=lambda t: t['lineno'])
    return tables


# ==========================================================================
# Static part: consumers, copies and draws
# ==========================================================================

def _names_loaded(node):
    return {n.id for n in ast.walk(node) if isinstance(n, ast.Name)}


_SCOPES = (ast.ListComp, ast.SetComp, ast.DictComp, ast.GeneratorExp, ast.Lambda)


def _free_loads(fn):
    """Names a function loads that it does not bind itself - its GLOBALS."""
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
    """Names in constants.py - assignments or functions - built from a table.

    Flow-insensitive is sound here: each name is bound once, by one top-level
    statement (constants.py has no module-level loops).
    """
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


_FILLERS = ('append', 'extend', 'add', 'update', 'insert', 'setdefault', 'appendleft')


class _Flow:
    """Flow-sensitive table derivation over a statement list (P-CONSUMER).

    `env` maps a name to the tables its CURRENT value derives from; `ever` is the
    union over the whole walk, used where a call site's exact state is not
    tracked (draw detection). Names absent from `env` fall back to `sources` -
    the tables themselves, constants-level aliases, and module aliases.
    """

    def __init__(self, sources, calls=None):
        self.sources = sources
        self.calls = calls or {}
        self.ever = {}

    def _look(self, env, name):
        return env[name] if name in env else self.sources.get(name, set())

    def tables_in(self, node, env):
        out = set()
        for n in ast.walk(node):
            if isinstance(n, ast.Name):
                out |= self._look(env, n.id)
            # A call to a module-level function yields what that function READS,
            # not just its name. The first version stopped at the name, so
            # `_T27_PLANS = _t27_admissible_plans()` derived from nothing and the
            # p-chart template's random.choice(_T27_PLANS) was an invisible draw.
            if isinstance(n, ast.Call) and isinstance(n.func, ast.Name):
                out |= self.calls.get(n.func.id, set())
        return out

    def _bind(self, env, name, val):
        env[name] = set(val)
        self.ever.setdefault(name, set()).update(val)

    def _fill(self, env, name, val):
        env[name] = self._look(env, name) | set(val)
        self.ever.setdefault(name, set()).update(val)

    def assign(self, env, target, val):
        if isinstance(target, ast.Name):
            self._bind(env, target.id, val)
        elif isinstance(target, (ast.Tuple, ast.List)):
            for elt in target.elts:
                self.assign(env, elt, val)
        elif isinstance(target, ast.Starred):
            self.assign(env, target.value, val)
        elif isinstance(target, (ast.Subscript, ast.Attribute)):
            base = target.value
            while isinstance(base, (ast.Subscript, ast.Attribute)):
                base = base.value
            if isinstance(base, ast.Name):
                extra = self.tables_in(target.slice, env) \
                    if isinstance(target, ast.Subscript) else set()
                self._fill(env, base.id, set(val) | extra)

    @staticmethod
    def _merge(a, b, keys_from):
        out = {}
        for k in set(a) | set(b):
            out[k] = keys_from._look(a, k) | keys_from._look(b, k)
        return out

    def walk(self, stmts, env):
        for st in stmts:
            if isinstance(st, (ast.FunctionDef, ast.AsyncFunctionDef, ast.ClassDef,
                               ast.Import, ast.ImportFrom)):
                continue
            if isinstance(st, ast.Assign):
                val = self.tables_in(st.value, env)
                for t in st.targets:
                    self.assign(env, t, val)
            elif isinstance(st, ast.AnnAssign) and st.value is not None:
                self.assign(env, st.target, self.tables_in(st.value, env))
            elif isinstance(st, ast.AugAssign):
                self.assign(env, st.target,
                            self.tables_in(st.value, env) | self.tables_in(st.target, env))
            elif isinstance(st, (ast.For, ast.AsyncFor, ast.While)):
                before = dict(env)
                for _ in range(2):
                    if not isinstance(st, ast.While):
                        self.assign(env, st.target, self.tables_in(st.iter, env))
                    self.walk(st.body, env)
                merged = self._merge(before, env, self)
                self.walk(st.orelse, merged)
                env.clear()
                env.update(merged)
            elif isinstance(st, ast.If):
                a, b = dict(env), dict(env)
                self.walk(st.body, a)
                self.walk(st.orelse, b)
                merged = self._merge(a, b, self)
                env.clear()
                env.update(merged)
            elif isinstance(st, (ast.With, ast.AsyncWith)):
                for item in st.items:
                    if item.optional_vars is not None:
                        self.assign(env, item.optional_vars,
                                    self.tables_in(item.context_expr, env))
                self.walk(st.body, env)
            elif isinstance(st, ast.Try):
                self.walk(st.body, env)
                for h in st.handlers:
                    self.walk(h.body, env)
                self.walk(st.orelse, env)
                self.walk(st.finalbody, env)
            elif isinstance(st, ast.Expr) and isinstance(st.value, ast.Call):
                f = st.value.func
                if isinstance(f, ast.Attribute) and f.attr in _FILLERS:
                    base = f.value
                    while isinstance(base, (ast.Subscript, ast.Attribute)):
                        base = base.value
                    if isinstance(base, ast.Name):
                        val = set()
                        for arg in list(st.value.args) + [k.value for k in st.value.keywords]:
                            val |= self.tables_in(arg, env)
                        self._fill(env, base.id, val)


def _base_sources(tables, calias):
    sources = {t: {t} for t in tables}
    for a, s in calias.items():
        # A name that is BOTH a table and built from tables stands for itself too.
        # The first flow-sensitive version replaced a table's own entry with the
        # tables it is built from, so RESISTOR_SERIES_BY_TOLERANCE (built from the
        # IEC lists) lost its consumer and MEDIA_VELOCITIES (built from C0) its draw.
        sources[a] = set(s) | ({a} if a in sources else set())
    return sources


def _call_tables(funcs, sources):
    """{module function: tables its reach reads}, through any name in `sources`."""
    out = {}
    for name, f in funcs.items():
        free = set().union(*(_free_loads(r) for r in _reach(f, funcs)))
        got = set().union(*(sources[x] for x in free if x in sources)) if free else set()
        if got:
            out[name] = got
    return out


def module_calls_and_aliases(tree, tables, calias):
    """(call map, module aliases) for one template module, to a fixpoint: a function
    that reads a module alias reaches that alias's tables, and a module alias bound
    from calling a function reaches the tables the function reads."""
    funcs = _module_functions(tree)
    base = _base_sources(tables, calias)
    aliases = {}
    for _ in range(8):
        sources = dict(base)
        for a, s in aliases.items():
            sources[a] = sources.get(a, set()) | s
        calls = _call_tables(funcs, sources)
        flow = _Flow(base, calls)
        env = {}
        flow.walk(tree.body, env)
        new = {n: s for n, s in env.items() if s and n not in base}
        if new == aliases:
            return calls, aliases
        aliases = new
    return calls, aliases


def module_aliases(tree, tables, calias):
    """Module-level names whose value, at the end of import, derives from a table."""
    return module_calls_and_aliases(tree, tables, calias)[1]


def _literals_in(fns):
    out = set()
    for f in fns:
        for n in ast.walk(f):
            if _is_num(n):
                out.add(float(n.value))
    return out


def template_reaches(branch_dir, module_root, tables, calias):
    """Yield (module name, template name, reach, sources, here, calls) for every
    template_* function under `branch_dir`.

    The one walk P-CONSUMER, P-DRAW and the C3.10 guard detector share, so the
    detector reads tables through exactly the aliases the census does rather than
    a copy of them - C1's selftest once re-implemented a detector and agreed with
    itself. `here` maps an alias to the tables it stands for; `sources` adds each
    table standing for itself; `calls` maps a module-level function to the
    tables it reads, which a name bound from calling it inherits.
    """
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
            calls, maliases = module_calls_and_aliases(tree, tables, calias)
            here = {a: set(s) | ({a} if a in tables else set()) for a, s in calias.items()}
            for a, s in maliases.items():
                here.setdefault(a, set()).update(s)
            sources = {t: {t} for t in tables}
            for a, s in here.items():
                sources[a] = sources.get(a, set()) | s
            for name, f in funcs.items():
                if name.startswith('template_'):
                    yield modname, name, _reach(f, funcs), sources, here, calls


def consumers_and_draws(branch_dir, module_root, tables, calias):
    """Apply P-CONSUMER and P-DRAW to every template module under `branch_dir`.

    Also returns, per template, its function source's numeric literals, which
    P-COPY checks a declaration against.
    """
    use = {t: {} for t in tables}
    literals = {}
    for modname, name, reach, sources, here, calls in template_reaches(branch_dir, module_root, tables, calias):
        literals[name] = (modname, _literals_in(reach))
        free = set().union(*(_free_loads(r) for r in reach))
        flow = _Flow(sources, calls)
        for r in reach:
            flow.walk(r.body, {})
        for t in tables:
            stand_ins = {t} | {a for a, s in here.items() if t in s}
            present = stand_ins & free
            if not present:
                continue
            via = sorted('direct' if x == t else f'alias {x}' for x in present)
            draws = set()
            for r in reach:
                for n in ast.walk(r):
                    if not (isinstance(n, ast.Call) and isinstance(n.func, ast.Attribute)
                            and isinstance(n.func.value, ast.Name)
                            and n.func.value.id == 'random'
                            and n.func.attr in ('choice', 'sample') and n.args):
                        continue
                    arg = n.args[0]
                    arg_tables = set()
                    for nm in _names_loaded(arg):
                        arg_tables |= flow.ever.get(nm, set()) | sources.get(nm, set())
                    if t not in arg_tables:
                        continue
                    is_sorted = any(isinstance(c, ast.Call)
                                    and isinstance(c.func, ast.Name)
                                    and c.func.id == 'sorted' for c in ast.walk(arg))
                    draws.add('sorted' if is_sorted else f'random.{n.func.attr}')
            use[t][name] = dict(via=via, draws=sorted(draws), module=modname)
    return use, literals


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
    # Precedence: what the probe CANNOT settle outranks what it can. The first
    # version checked RESTATED before STRUCTURAL and ERROR, so a field restated on
    # some seeds and crashing or rewording the question on others reported
    # RESTATED - and G-2's escape (a crash read as harmless) stayed open one level
    # below the rollup that was fixed for it. Caught by the WINDOW2 plant.
    for v in ('HIDDEN', 'IMPORT_ERROR', 'ERROR', 'STRUCTURAL', 'RESTATED'):
        if c.get(v):
            return v.replace('_', '-')
    return 'UNUSED'


def rollup(verdicts, has_name_consumers, has_copies):
    """Table-level measurement from field verdicts (see P-GIVEN)."""
    if not has_name_consumers:
        return 'COPIED' if has_copies else 'UNCONSUMED'
    if 'HIDDEN' in verdicts:
        return 'SOME-HIDDEN'
    if verdicts and 'RESTATED' in verdicts and all(v in ('RESTATED', 'UNUSED') for v in verdicts):
        return 'ALL-RESTATED'
    if 'ERROR' in verdicts or 'IMPORT-ERROR' in verdicts:
        return 'ERROR'
    if 'STRUCTURAL' in verdicts:
        return 'INDETERMINATE'
    return 'NO-EFFECT'


# ==========================================================================
# The census
# ==========================================================================

HEADER_FIELD_RE = re.compile(r'#\s*@(units|domain|kind|given|copied-in)\s*:\s*(.*)')


def header_fields(header):
    out = {}
    for ln in header.splitlines():
        m = HEADER_FIELD_RE.match(ln.strip())
        if m:
            out[m.group(1)] = (out.get(m.group(1), '') + ' ' + m.group(2)).strip()
    return out


def parse_copies(decl):
    """`@copied-in: template_a 0.2; template_b 57000` -> [(template_a, 0.2), ...]."""
    out = []
    for part in (decl or '').split(';'):
        bits = part.split()
        if not bits:
            continue
        try:
            out.append((bits[0], float(bits[1])))
        except (IndexError, ValueError):
            out.append((bits[0], None))
    return out


def classify(kind, measured, given=None):
    """C1.3's rule, from the declared @kind, the P-GIVEN measurement and @given.

    The given-values rule excuses a value from CORRECTNESS-dependence: restated
    in the question, a wrong number cannot make the gold answer disagree with
    the question. It does not excuse a value from being TRUE - C2 shipped 274
    items that were self-consistent and false (phaseC2_summary.md §6). So only a
    `range` can be a plausibility window, and only when every consumer is
    measured to restate it, or reads it without anything changing or raising
    (a guard). A crash, a structural change or a copy tells the probe the table
    is read and not how, so it needs a human verdict, declared with its
    evidence as @given (C1 Reviewer G, G-2). A declaration never overrides a
    measured HIDDEN.
    """
    if kind is None:
        return 'UNDECLARED', 'no @kind declared'
    kind = kind.split()[0]
    if kind not in KINDS:
        return 'UNDECLARED', f'unknown @kind {kind!r}'
    if measured == 'UNCONSUMED':
        return 'UNCONSUMED', 'no template consumes it (named or declared copy)'
    if kind == 'range':
        if measured == 'SOME-HIDDEN':
            return 'CITATION', ('declared a range but measured SOME-HIDDEN: the '
                                'given-values rule does not apply'
                                + (' - the @given declaration is contradicted' if given else ''))
        if measured == 'ALL-RESTATED':
            return 'PLAUSIBILITY', 'a sampling window every consumer restates'
        if measured == 'NO-EFFECT':
            return 'PLAUSIBILITY', 'NO-EFFECT: read without changing or raising - a guard'
        verdict = (given or '').split()[0] if given else ''
        if verdict in ('stated', 'guard'):
            return 'PLAUSIBILITY', f'{measured}; declared @given: {given}'
        return 'REVIEW', (f'{measured}: the probe cannot tell a guard from a crash or a '
                          f'copy - declare @given: stated|guard (evidence)')
    if kind in ('property', 'standard', 'measured-constant'):
        reasons = {'SOME-HIDDEN': 'consumed without being stated',
                   'COPIED': 'consumed through a copied literal'}
        return 'CITATION', reasons.get(measured,
                                       f'{measured}: asserts a fact about a named real entity')
    if kind == 'mathematical':
        return 'DERIVATION', 'computable from its definition'
    if kind == 'defined':
        return 'DEFINITION', 'exact by definition or convention'
    return 'DOMAIN', 'a validity bound for another table; checked by its C3.2 suite'


def census(seeds=40, probe_on=True, branches=BRANCHES, root=BRANCHES_DIR,
           module_root=REPO, package='data.templates.branches'):
    report = []
    for branch in branches:
        src = open(os.path.join(root, branch, 'constants.py'), encoding='utf-8').read()
        cmod = importlib.import_module(f'{package}.{branch}.constants' if package
                                       else f'{branch}.constants')
        tables = numeric_tables(src, dict(vars(cmod)))
        names = sorted(t['name'] for t in tables)
        calias = constants_aliases(src, names)
        use, lits = consumers_and_draws(os.path.join(root, branch), module_root, names, calias)
        baseline = {}
        for t in tables:
            obj = getattr(cmod, t['name'])
            declared = header_fields(t.pop('header'))
            t.update(branch=branch, type=type(obj).__name__, fields=fields(obj),
                     consumers=dict(use[t['name']]), declared=declared,
                     probe={}, seeds=seeds if probe_on else 0, copies=[], copy_errors=[])
            for tid, lit in parse_copies(declared.get('copied-in')):
                if lit is None:
                    t['copy_errors'].append(f'@copied-in entry {tid!r} names no numeric literal')
                elif tid not in lits:
                    t['copy_errors'].append(f'@copied-in names {tid}, which is not a template here')
                elif lit not in lits[tid][1]:
                    t['copy_errors'].append(f'{tid} contains no literal {lit} - the declared '
                                            f'copy is gone; remove the declaration')
                else:
                    t['copies'].append((tid, lit))
            named = {tid: u['module'] for tid, u in t['consumers'].items()}
            for tid, lit in t['copies']:
                t['consumers'].setdefault(tid, dict(via=[f'copied literal {lit}'], draws=[],
                                                    module=lits[tid][0]))
            # Order is a property of a container; a scalar table has none to depend on.
            t['draws_by_order'] = (sorted(tid for tid, u in t['consumers'].items()
                                          if any(d.startswith('random.') for d in u['draws']))
                                   if isinstance(obj, (dict, list, tuple)) else [])
            if probe_on and named:
                for tid, m in named.items():
                    if tid not in baseline:
                        fn = getattr(importlib.import_module(m), tid)
                        baseline[tid] = [_gen(fn, s) for s in range(seeds)]
                for fld in t['fields']:
                    res = probe(cmod, t['name'], fld, named, seeds, baseline)
                    t['probe'][fld] = {tid: dict(counts=c, verdict=field_verdict(c))
                                       for tid, c in res.items()}
            verdicts = [v['verdict'] for f in t['probe'].values() for v in f.values()]
            if named and not probe_on:
                t['measured'] = 'NOT-PROBED'
            else:
                t['measured'] = rollup(verdicts, bool(named), bool(t['copies']))
            t['class'], t['class_reason'] = classify(declared.get('kind'), t['measured'],
                                                     declared.get('given'))
            if t['copy_errors']:
                t['class'], t['class_reason'] = 'REVIEW', '; '.join(t['copy_errors'])
            report.append(t)
    return report


def print_report(report):
    print('PROVENANCE COVERAGE (P-TABLE, P-LITERAL, P-TAGGED; P-TABLE-LIVE tables excluded)')
    print(f"{'branch':24s} {'tables':>6s} {'tagged':>6s} {'literals':>8s} {'in tagged':>9s}")
    for b in BRANCHES:
        rows = [t for t in report if t['branch'] == b and not t['live_only']]
        if rows:
            print(f"{b:24s} {len(rows):6d} {sum(t['tagged'] for t in rows):6d} "
                  f"{sum(t['literals'] for t in rows):8d} "
                  f"{sum(t['literals'] for t in rows if t['tagged']):9d}")
    live = [t['name'] for t in report if t['live_only']]
    print(f"P-TABLE-LIVE (classified, not counted above): {', '.join(live) or 'none'}")
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
        for tid, lit in t['copies']:
            print(f"       copied literal {lit} in: {tid}")
        if t['declared'].get('given'):
            print(f"       @given: {t['declared']['given'][:100]}")
        for e in t['copy_errors']:
            print(f"       COPY DECLARATION ERROR: {e}")
    print()
    counts = {}
    for t in report:
        counts[t['class']] = counts.get(t['class'], 0) + 1
    print('CLASSIFICATION', ', '.join(f'{k} {v}' for k, v in sorted(counts.items())))


def write_markdown(report, path, seeds):
    """The committed C1.3 report. Regenerated, never hand-edited."""
    lit = [t for t in report if not t['live_only']]
    out = [
        '# Phase C1.3 — constants-table census and classification',
        '',
        '**Generated, not written.** Regenerate with',
        '',
        f'    python -m tests.constants_integrity.census --seeds {seeds} --markdown '
        'docs/re-implementation-sep/track-b/phaseC1_census.md',
        '',
        'Every term below - *numeric table*, *literal*, *tagged*, *consumer*, *copy*,',
        '*draws by order*, *restated*, *hidden* - is a predicate defined in the docstring',
        'of `tests/constants_integrity/census.py`; the classification rule is spec §C1.3.',
        '',
        '## Coverage (P-TABLE)',
        '',
        '| branch | numeric tables | tagged | numeric literals | in tagged tables |',
        '|---|---:|---:|---:|---:|',
    ]
    for b in BRANCHES:
        rows = [t for t in lit if t['branch'] == b]
        out.append(f"| {b} | {len(rows)} | {sum(t['tagged'] for t in rows)} | "
                   f"{sum(t['literals'] for t in rows)} | "
                   f"{sum(t['literals'] for t in rows if t['tagged'])} |")
    out.append(f"| **all** | **{len(lit)}** | **{sum(t['tagged'] for t in lit)}** | "
               f"**{sum(t['literals'] for t in lit)}** | "
               f"**{sum(t['literals'] for t in lit if t['tagged'])}** |")
    live = [f"`{t['branch'].split('_')[0]}.{t['name']}`" for t in report if t['live_only']]
    out += ['', f"Plus **{len(live)}** P-TABLE-LIVE table(s), classified below and not "
                f"counted above: {', '.join(live) or 'none'}.",
            '', f'## Classification ({len(report)} tables)', '',
            '| class | tables | literals |', '|---|---:|---:|']
    for cls in ('CITATION', 'PLAUSIBILITY', 'DERIVATION', 'DEFINITION', 'DOMAIN',
                'UNCONSUMED', 'REVIEW', 'UNDECLARED'):
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
        why = t['class_reason'].replace('|', '/')
        if t['copies']:
            why += ' [copy: ' + ', '.join(f"{tid.replace('template_', '')} {lit}"
                                          for tid, lit in t['copies']) + ']'
        out.append(f"| {t['branch'].split('_')[0]} | `{t['name']}` | {t['literals']} | "
                   f"{', '.join(t['tags']) or '-'} | {len(t['consumers'])} | "
                   f"{len(t['draws_by_order'])} | {t['measured']} | "
                   f"{t['declared'].get('kind', '-').split(' ')[0]} | **{t['class']}** | "
                   f"{why} | {', '.join(hidden) or '-'} |")
    with open(path, 'w', encoding='utf-8', newline='\n') as fh:
        fh.write('\n'.join(out) + '\n')


# ==========================================================================
# Self-test: planted defects, written from the predicate definitions above,
# and run through the SAME functions the census uses - never a copy of them.
# ==========================================================================

_SELFTEST_CONSTANTS = '''
import math

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
TABLE_A = {"x": 1.5, "y": 2.5}

# @kind: range
HIDDEN_KEYS = [0.90, 0.95]

# @kind: range
# @given: stated (planted: the drawn level is printed in the question)
DECLARED_KEYS = [0.90, 0.95]

# @kind: range
WINDOW2 = (4.9, 5.1)

# @kind: standard
# @copied-in: template_copies_ratio 0.37
RATIO = 0.37

# @kind: standard
# @copied-in: template_restates 0.41
MISDECLARED = 0.41

# @kind: range
ANGLE = (-math.pi, math.pi)

BASE = [1.0, 2.0]
SERIES = {5: BASE}
SCALE = 3.0
PLAN_WIN = (20, 23)
ROW_SET = {"a": 1.5, "b": 2.5}


def density_of(name):
    return DENSITY[name]
'''

_SELFTEST_TEMPLATES = '''
import random
from selftest_branch.constants import (DENSITY, PAIRS, WINDOW, ORPHAN, LEVELS,
                                       SPREAD, TABLE_A, HIDDEN_KEYS, DECLARED_KEYS,
                                       WINDOW2, ANGLE, SERIES, SCALE, density_of,
                                       PLAN_WIN, ROW_SET)

_COPIED = {k: 2 * v for k, v in LEVELS.items()}
_Z = {0.90: 1.28, 0.95: 1.64}

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

def template_hidden_lookup():
    a = random.choice(HIDDEN_KEYS)
    return "Find z.", f"**Answer:** {_Z[a]}"

def template_declared_lookup():
    a = random.choice(DECLARED_KEYS)
    return f"Service level {a}. Find z.", f"**Answer:** {_Z[a]}"

def template_structural():
    x = random.uniform(*WINDOW2)
    label = "high" if x > 5 else "low"
    return f"A {label} case with x = {x:.3f}.", f"**Answer:** {x:.3f}"

def template_copies_ratio():
    s = round(random.uniform(1, 2), 2)
    return f"s = {s}", f"**Answer:** {0.37 * s:.4f}"

def template_angle():
    ph = round(random.uniform(*ANGLE), 3)
    return f"phase {ph}", f"**Answer:** {ph}"

_S5 = SERIES[5]

def template_series():
    return "S", f"**Answer:** {random.choice(_S5)}"

def template_scaled():
    xs = [SCALE * 1, SCALE * 2]
    return "X", f"**Answer:** {random.choice(xs)}"

_FILLED_A = []
for _k in TABLE_A:
    _v = TABLE_A[_k]
    _FILLED_A.append(_v)
_FILLED_B = []
for _v in (7.0, 8.0):
    _FILLED_B.append(_v)

def template_reads_a():
    return "A", f"**Answer:** {random.choice(_FILLED_A)}"

def template_reads_b_only():
    return "B", f"**Answer:** {random.choice(_FILLED_B)}"


def _build_plans():
    plans = []
    for m in range(PLAN_WIN[0], PLAN_WIN[1] + 1):
        plans.append(m)
    return plans


_PLANS = _build_plans()            # module level: the table enters only through the call


def template_module_call_draw():
    m = random.choice(_PLANS)
    return f"m = {m}", f"**Answer:** {2 * m}"


def _pick_rows():
    return sorted(ROW_SET.values())


def template_in_template_call_draw():
    rows = _pick_rows()            # inside the template: the result is drawn from
    v = random.choice(rows)
    return f"v = {v}", f"**Answer:** {v}"
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
        # -- P-TABLE / P-TABLE-LIVE / P-LITERAL / P-TAGGED ---------------------
        st = {t['name']: t for t in static_tables(_SELFTEST_CONSTANTS)}
        for n, want in {'DENSITY': True, 'PAIRS': False, 'WINDOW': True,
                        'ORPHAN': False, 'UNRELATED': True}.items():
            if st[n]['tagged'] != want:
                failures.append(f'P-TAGGED {n}: got {st[n]["tagged"]}, planted {want}')
        if 'NAMES' in st:
            failures.append('P-TABLE counted a table with no numeric literal')
        if 'ANGLE' in st:
            failures.append('P-TABLE counted ANGLE, which has no literal')
        if st['PAIRS']['literals'] != 4:
            failures.append(f'P-LITERAL PAIRS: {st["PAIRS"]["literals"]} != 4')
        live = {t['name']: t for t in numeric_tables(_SELFTEST_CONSTANTS)}
        if 'ANGLE' not in live or not live['ANGLE']['live_only'] or 'NAMES' in live:
            failures.append('P-TABLE-LIVE: ANGLE must be a live-only table and NAMES none')

        # -- the census's own code path --------------------------------------
        rep = {t['name']: t for t in census(seeds=12, branches=('selftest_branch',),
                                            root=tmp, module_root=tmp, package='')}

        # P-CONSUMER, including name reuse across two module-level loops
        want_cons = {
            'DENSITY': {'template_restates', 'template_constants_accessor'},
            'PAIRS': {'template_hides_other_field_of_restated_row'},
            'WINDOW': {'template_window'},
            'LEVELS': {'template_import_time_copy'},
            'SPREAD': {'template_loop_alias'},
            'TABLE_A': {'template_reads_a'},
            'ORPHAN': set(),
            'RATIO': {'template_copies_ratio'},
            'ANGLE': {'template_angle'},
            'SERIES': {'template_series'},       # a table read through a module alias
            'BASE': {'template_series'},         # ... and the table it is built from
        }
        for tbl, want in want_cons.items():
            got = set(rep[tbl]['consumers'])
            if got != want:
                failures.append(f'P-CONSUMER {tbl}: got {sorted(got)}, planted {sorted(want)}')

        # P-DRAW
        want_draw = {('DENSITY', 'template_restates'): ['random.choice'],
                     ('PAIRS', 'template_hides_other_field_of_restated_row'): ['random.choice'],
                     ('LEVELS', 'template_import_time_copy'): ['sorted'],
                     ('DENSITY', 'template_constants_accessor'): [],
                     ('TABLE_A', 'template_reads_a'): ['random.choice'],
                     ('SERIES', 'template_series'): ['random.choice'],
                     ('BASE', 'template_series'): ['random.choice'],
                     # a table reaching a draw only through a function call - at
                     # module level, and inside the template (census fix, C3.10)
                     ('PLAN_WIN', 'template_module_call_draw'): ['random.choice'],
                     ('ROW_SET', 'template_in_template_call_draw'): ['random.choice']}
        for (tbl, tid), want in want_draw.items():
            got = rep[tbl]['consumers'].get(tid, {}).get('draws')
            if got != want:
                failures.append(f'P-DRAW {tbl} in {tid}: got {got}, planted {want}')

        if rep['SCALE']['draws_by_order']:
            failures.append(f"P-DRAW: scalar SCALE reported drawn by order: {rep['SCALE']['draws_by_order']}")

        # P-GIVEN field verdicts
        want_given = {
            ('DENSITY', '[*]', 'template_restates'): 'RESTATED',
            ('DENSITY', '[*]', 'template_constants_accessor'): 'HIDDEN',
            ('PAIRS', '[*][0]', 'template_hides_other_field_of_restated_row'): 'RESTATED',
            ('PAIRS', '[*][1]', 'template_hides_other_field_of_restated_row'): 'HIDDEN',
            ('WINDOW', '[0]', 'template_window'): 'RESTATED',
            ('LEVELS', '[*]', 'template_import_time_copy'): 'HIDDEN',
            ('SPREAD', '[*]', 'template_loop_alias'): 'HIDDEN',
            ('HIDDEN_KEYS', '[*]', 'template_hidden_lookup'): 'ERROR',
        }
        for (tbl, fld, tid), want in want_given.items():
            got = rep[tbl]['probe'].get(fld, {}).get(tid, {}).get('verdict')
            if got != want:
                failures.append(f'P-GIVEN {tbl}{fld} in {tid}: got {got}, planted {want}')

        # rollup + classify: a crash is not a guard; a structural change is not
        # a guard; a declaration settles it; a copy is consumption
        want_class = {
            'DENSITY': ('SOME-HIDDEN', None),
            'HIDDEN_KEYS': ('ERROR', 'REVIEW'),
            'DECLARED_KEYS': ('ERROR', 'PLAUSIBILITY'),
            'WINDOW2': ('INDETERMINATE', 'REVIEW'),
            'RATIO': ('COPIED', 'CITATION'),
            'MISDECLARED': ('UNCONSUMED', 'REVIEW'),
            'ANGLE': ('ALL-RESTATED', 'PLAUSIBILITY'),
        }
        for tbl, (m, c) in want_class.items():
            got = (rep[tbl]['measured'], rep[tbl]['class'])
            if got[0] != m or (c is not None and got[1] != c):
                failures.append(f'rollup/classify {tbl}: got {got}, planted ({m}, {c})')
        if not rep['MISDECLARED']['copy_errors']:
            failures.append('P-COPY: a declared copy absent from the template was not reported')

        # tag vs class (C3.1): a POLICY tag on each measured route to a
        # non-PLAUSIBILITY class, on the real rows; controls must stay silent.
        unconsumed = dict(rep['ANGLE'], measured='UNCONSUMED',
                          **dict(zip(('class', 'class_reason'), classify('range', 'UNCONSUMED'))))
        for label, row, want in (
                ('a copied literal (COPIED -> CITATION)', rep['RATIO'], 1),
                ('a crash (ERROR -> REVIEW)', rep['HIDDEN_KEYS'], 1),
                ('control: a restated window', rep['ANGLE'], 0),
                ('control: a declared @given', rep['DECLARED_KEYS'], 0),
                ('control: an unconsumed range', unconsumed, 0)):
            got = len(tag_class_conflicts([dict(row, tags=['POLICY'])]))
            if got != want:
                failures.append(f'tag/class: POLICY on {label}: {got} conflict(s), planted {want}')

        # the probe must leave everything as it found it
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


def tag_class_conflicts(report):
    """A tag that contradicts the measured class (C3.1).

    `[POLICY: sampling-only]` claims a table is a sampling window whose values
    cannot make gold disagree with the question. The census MEASURES that claim:
    only PLAUSIBILITY supports it. A POLICY tag on a CITATION or REVIEW table is
    the record contradicting the evidence - the exact shape of C2's green suite -
    and it stays silent unless something compares the two, because R7 in
    test_citations_resolve.py compares the tag with the DECLARED @kind only.
    UNCONSUMED is exempt: with no consumer there is nothing to contradict.
    """
    return [f"{t['branch']}.{t['name']}: tagged [POLICY: sampling-only] but classified "
            f"{t['class']} - {t['class_reason']}"
            for t in report
            if 'POLICY' in t['tags'] and t['class'] not in ('PLAUSIBILITY', 'UNCONSUMED')]


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument('--seeds', type=int, default=40)
    ap.add_argument('--no-probe', action='store_true')
    ap.add_argument('--json')
    ap.add_argument('--markdown')
    ap.add_argument('--branch', action='append')
    ap.add_argument('--check', action='store_true',
                    help='exit 1 if any table is REVIEW or UNDECLARED')
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
    if args.check:
        bad = [f"{t['branch']}.{t['name']}: {t['class']} - {t['class_reason']}"
               for t in report if t['class'] in ('REVIEW', 'UNDECLARED')]
        bad += tag_class_conflicts(report)
        for b in bad:
            print('  - ' + b)
        print('check: all classified' if not bad else f'check: {len(bad)} FAILURE(S)')
        return 1 if bad else 0
    return 0


if __name__ == '__main__':
    sys.exit(main())
