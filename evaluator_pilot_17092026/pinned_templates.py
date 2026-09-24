"""Run a pilot command against the templates as they stood when the slice was frozen.

    python -m evaluator_pilot_17092026.pinned_templates --check
    python -m evaluator_pilot_17092026.pinned_templates run_evaluator e4

WHY THIS EXISTS. E3 and E4 do not read milestones out of a file. `evaluators/
milestones.py` regenerates each frozen item from the repo's own templates at its
recorded seed, checks the regenerated question and solution BYTE-IDENTICAL to
slice/manifest.jsonl, and only then derives milestones from the template's captured
internals. That check is the guarantee that a milestone is the item's own value and
not something read off the text.

The repo's templates have moved on since the freeze. `3fad887` (2026-09-23, "bind
displayed operands and remove display ties in 45 templates") and `fc1a6dc`
(2026-09-24, "Screen pass 1 fixes") rewrote 30 template files, five of which the
pilot's 60 items come from. The edits are to how a derivation is displayed - the
Reynolds template now writes `Re = (654.8 * 13.4 * 0.73) / 3.13e-04 = 20,464,069`
on one line where the frozen text writes two - so 17 of the 60 items no longer
reproduce, and milestones.py refuses, correctly: the experts annotated the frozen
text, and a milestone set derived from different text would not be the one they
were scored against.

So the five files are read back out of git at PINNED_COMMIT - the repo's state on
the day the pilot scored E3 and E4 - and installed under their own module names
before anything imports them. Nothing in the working tree is touched: `git show`
writes nothing, and the extracted copies live in the (uncommitted) score cache.

This does not relax milestones.py's check - it feeds it the right templates. With
the pin in place all 60 items reproduce byte-identically, which `--check` prints.
"""
from __future__ import annotations

import argparse
import importlib.util
import json
import os
import runpy
import subprocess
import sys

_HERE = os.path.dirname(os.path.abspath(__file__))
_ROOT = os.path.dirname(_HERE)
CACHE = os.path.join(_HERE, 'scores', '_cache', 'pinned_templates')

# The repo's state on 2026-09-19, when the pilot ran E3 and E4 (26f9048, "E1 done").
PINNED_COMMIT = '26f9048'

# The five template files the frozen items come from that later commits rewrote.
# Any other file the templates import is taken from the working tree as it is; the
# byte-identity check on all 60 items is what says that is enough.
FILES = (
    'data/templates/branches/chemical_engineering/reaction_kinetics/stoichiometry.py',
    'data/templates/branches/chemical_engineering/transport_phenomena/viscosity_and_momentum_transport.py',
    'data/templates/branches/electrical_engineering/electromagnetics_and_waves/magnetostatics.py',
    'data/templates/branches/electrical_engineering/electromagnetics_and_waves/waves_and_phasors.py',
    'data/templates/branches/industrial_engineering/stochastic_operations/queueing_systems.py',
)


def _extract(path: str) -> str:
    """The file as of PINNED_COMMIT, cached on disk. Read-only with respect to git."""
    os.makedirs(CACHE, exist_ok=True)
    out = os.path.join(CACHE, '%s__%s' % (PINNED_COMMIT, path.replace('/', '__')))
    if not os.path.exists(out):
        blob = subprocess.run(['git', 'show', '%s:%s' % (PINNED_COMMIT, path)],
                              cwd=_ROOT, capture_output=True, check=True).stdout
        with open(out, 'wb') as fh:
            fh.write(blob)
    return out


def install() -> list:
    """Put the pinned modules in sys.modules, so importlib hands them out instead."""
    if _ROOT not in sys.path:
        sys.path.insert(0, _ROOT)
    done = []
    for path in FILES:
        mod = path[:-3].replace('/', '.')
        spec = importlib.util.spec_from_file_location(mod, _extract(path))
        m = importlib.util.module_from_spec(spec)
        sys.modules[mod] = m
        spec.loader.exec_module(m)
        done.append(mod)
    return done


def check() -> int:
    """How many of the 60 frozen items reproduce, without the pin and with it."""
    items = [json.loads(l) for l in
             open(os.path.join(_HERE, 'slice', 'manifest.jsonl'), encoding='utf-8')]

    def misses():
        from tests.template_integrity.core import discover, generate
        refs = {r.template_id: r for r in discover()}
        out = []
        for it in items:
            inst = generate(refs[it['template_id']], it['seed'], capture=True)
            if inst.question != it['question'] or inst.solution != it['solution']:
                out.append(it['item_id'])
        return out

    if _ROOT not in sys.path:
        sys.path.insert(0, _ROOT)
    before = misses()
    print('working tree as it is : %d of %d items do not reproduce' % (len(before), len(items)))
    for x in before:
        print('    %s' % x)
    print('pinned at %s    : ' % PINNED_COMMIT, end='')
    install()
    after = misses()
    print('%d of %d items do not reproduce' % (len(after), len(items)))
    for x in after:
        print('    %s' % x)
    return 1 if after else 0


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument('--check', action='store_true',
                    help='report reproduction with and without the pin, and stop')
    ap.add_argument('module', nargs='?', help='a pilot module to run, e.g. run_evaluator')
    ap.add_argument('args', nargs=argparse.REMAINDER)
    a = ap.parse_args()
    if a.check or not a.module:
        return check()
    install()
    target = 'evaluator_pilot_17092026.' + a.module
    sys.argv = [target] + a.args
    runpy.run_module(target, run_name='__main__', alter_sys=True)
    return 0


if __name__ == '__main__':
    raise SystemExit(main())
