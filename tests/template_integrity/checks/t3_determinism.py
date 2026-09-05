"""T3 - Determinism.

`random.seed(s)` must fully determine a template's output.

The check generates in **two separate processes** and compares digests. Separate
processes is not incidental: a template that draws from numpy's global
generator, or that depends on per-process state, agrees with itself inside one
process and disagrees across two. `template_levenspiel_plot_interpretation`
calls `np.random.uniform`, which `random.seed()` does not control - the audit's
finding is that the seed recorded beside every generated item cannot regenerate
it.

Performance note
----------------
Only TWO independent process states are needed, not two per template. An
earlier version spawned a subprocess pair per template - 300 interpreter starts,
each re-importing numpy and scipy, for ~500s on the full corpus. Batching every
template into one child and running that child twice does the identical
comparison in roughly 10s. The rigour is unchanged: the same digests are
compared across the same two independent numpy states.
"""
from __future__ import annotations

import dataclasses
import json
import os
import subprocess
import tempfile
import sys
from dataclasses import dataclass, field

from ..core import REPO_ROOT, TemplateRef

# Child program: generate every requested template x seed and print one digest
# per (template, seed). Run twice; compare.
_CHILD = r'''
import hashlib, json, sys
sys.path.insert(0, sys.argv[1])
from tests.template_integrity.core import TemplateRef, generate

# Refs are passed in full rather than rediscovered in the child. Rediscovering
# would silently exclude anything outside data/templates/branches/, so a
# candidate or scratch template could not be checked for determinism before it
# was committed - exactly when the check is most useful (adversary review, F2).
with open(sys.argv[2], encoding="utf-8") as _f:
    specs = json.load(_f)
seeds = [int(x) for x in sys.argv[3].split(",")]
out = {}
for spec in specs:
    ref = TemplateRef(**spec)
    tid = ref.template_id
    per = {}
    for s in seeds:
        # seed_numpy=False is the whole point: the property under test is that
        # random.seed(s) ALONE fixes the output. Seeding numpy here would make
        # an unseeded np.random draw reproduce across both processes, hiding it.
        inst = generate(ref, s, capture=False, seed_numpy=False)
        if inst.error:
            per[s] = "ERROR:" + inst.error
        else:
            per[s] = hashlib.sha256(
                (inst.question + "\x1e" + inst.solution).encode("utf-8")).hexdigest()
    out[tid] = per
sys.stdout.write("\x02" + json.dumps(out) + "\x03")
'''


@dataclass
class DeterminismResult:
    template_id: str
    seeds_checked: int = 0
    mismatched_seeds: list[int] = field(default_factory=list)
    errored_seeds: list[int] = field(default_factory=list)
    note: str = ''

    @property
    def passed(self) -> bool:
        return not self.mismatched_seeds and not self.errored_seeds


def _run_child(refs: list[TemplateRef], seeds: list[int]) -> dict[str, dict[int, str]]:
    # Ref specs go via a temp file, not argv: 150 of them exceed the Windows
    # command-line length limit (WinError 206).
    fd, spec_path = tempfile.mkstemp(suffix='.json', prefix='t3_specs_')
    try:
        with os.fdopen(fd, 'w', encoding='utf-8') as f:
            json.dump([dataclasses.asdict(r) for r in refs], f)
        proc = subprocess.run(
            [sys.executable, '-c', _CHILD, REPO_ROOT,
             spec_path, ','.join(str(s) for s in seeds)],
            capture_output=True, text=True, cwd=REPO_ROOT,
            env={**os.environ, 'PYTHONIOENCODING': 'utf-8',
                 # vary the hash seed so dict/set iteration-order dependence surfaces
                 'PYTHONHASHSEED': 'random'},
            timeout=3600,
        )
    finally:
        try:
            os.unlink(spec_path)
        except OSError:
            pass
    if proc.returncode != 0:
        raise RuntimeError(f'T3 child failed: {proc.stderr[-1500:]}')
    body = proc.stdout
    start, end = body.find('\x02'), body.rfind('\x03')
    if start < 0 or end < 0:
        raise RuntimeError(f'T3 child produced no payload: {proc.stdout[-800:]}')
    raw = json.loads(body[start + 1:end])
    return {tid: {int(k): v for k, v in per.items()} for tid, per in raw.items()}


def run_many(refs: list[TemplateRef],
             seeds: range | list[int] = range(200)) -> dict[str, DeterminismResult]:
    """Check every template in two process runs total."""
    seeds = list(seeds)
    ids = [r.template_id for r in refs]
    a = _run_child(refs, seeds)
    b = _run_child(refs, seeds)
    results: dict[str, DeterminismResult] = {}
    for tid in ids:
        res = DeterminismResult(template_id=tid, seeds_checked=len(seeds))
        pa, pb = a.get(tid, {}), b.get(tid, {})
        for s in seeds:
            va, vb = pa.get(s), pb.get(s)
            if va is None or vb is None or str(va).startswith('ERROR:'):
                res.errored_seeds.append(s)
            elif va != vb:
                res.mismatched_seeds.append(s)
        if res.mismatched_seeds:
            res.note = (f'{len(res.mismatched_seeds)}/{len(seeds)} seeds differ '
                        f'between two processes')
        results[tid] = res
    return results


def run(ref: TemplateRef, seeds: range | list[int] = range(200)) -> DeterminismResult:
    """Single-template convenience wrapper. Prefer run_many for a sweep."""
    return run_many([ref], seeds)[ref.template_id]
