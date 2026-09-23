"""Layer 2 - build every expert's queue: their branch's templates plus planted defects.

    python -m template_annotation_23092026.layer2.build_tasks --check-plants   # verify plants only
    python -m template_annotation_23092026.layer2.build_tasks                  # build tasks/

Produces, under layer2/tasks/ (git-ignored, rebuilt on demand):

    pool.json         what experts see: {code: {branch, area, source, instances[5], gold_numbers}}
    assignment.json   {annotator_id: {branch, codes[...] in the order to work them}}
    keyfile.jsonl     code -> real template id, or plant id + base + defect; NEVER shipped

Design (D-095):
  * each expert reviews all 30 templates of their own branch plus the branch's four
    planted defects (plants/CONTRACT.md), 34 items, shuffled by a seed derived from
    their id so fatigue does not align across the three experts of a branch;
  * instance 0 of every item is the hand-check instance, the same for all three
    experts of a branch, so their hand checks are comparable; instances 1-4 are
    there to cycle through;
  * codes are opaque (BLAKE2b of a build secret and the item id) so a plant cannot
    be told from a real template by its name;
  * instances are precomputed here, so the bundle an expert receives runs no
    template code; it is JSON plus the app.

Plants are verified against the contract on every build: exactly-once edits, 50 seeds
without an exception, output different from the base on every seed, deterministic in
two processes.
"""
from __future__ import annotations

import argparse
import hashlib
import importlib
import inspect
import json
import os
import random
import re
import secrets
import subprocess
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
REPO = HERE.parent.parent
if str(REPO) not in sys.path:
    sys.path.insert(0, str(REPO))

from tests.template_integrity.core import NUM_RE, answer_block, discover, generate, seed_all  # noqa: E402

TASKS = HERE / 'tasks'
PLANTS_DIR = HERE / 'plants'
ROSTER = HERE / 'annotators.json'
PILOT_ROSTER = REPO / 'evaluator_pilot_17092026' / 'annotation' / 'annotators.json'
INSTANCE_SEEDS = (2101, 2102, 2103, 2104, 2105)   # 0 is the hand-check instance
PLANT_SEEDS = (3101, 3102, 3103, 3104, 3105)
PREFIX = {'chemical_engineering': 'che', 'civil_engineering': 'civ', 'electrical_engineering': 'ele',
          'industrial_engineering': 'ind', 'mechanical_engineering': 'mec'}


def roster() -> list[dict]:
    path = ROSTER if ROSTER.exists() else PILOT_ROSTER
    data = json.loads(path.read_text(encoding='utf8'))
    return data['annotators'] if isinstance(data, dict) else data


def gold_numbers(solution: str) -> list[float]:
    block, _ = answer_block(solution)
    out = []
    for tok in NUM_RE.findall(block):
        try:
            out.append(float(tok.replace(',', '')))
        except ValueError:
            pass
    return out


# ------------------------------------------------------------------- plants

def apply_edits(source: str, edits: list, label: str) -> str:
    for old, new in edits:
        n = source.count(old)
        if n != 1:
            raise SystemExit(f'{label}: edit {old!r} occurs {n} times in the base source (must be exactly 1)')
        source = source.replace(old, new)
    return source


def mutated_function(ref, plant: dict):
    """The plant as a callable, executed in the base module's namespace."""
    fn = ref.load()
    src = inspect.getsource(fn)
    src = apply_edits(src, plant.get('edits', []), plant['plant_id'])
    src = apply_edits(src, plant.get('reskin', []), plant['plant_id'] + ' reskin')
    src = re.sub(r'^\s+', '', src, count=1) if src.startswith(' ') else src
    module = sys.modules[ref.module]
    ns = dict(module.__dict__)
    exec(compile(src, f'<plant {plant["plant_id"]}>', 'exec'), ns)
    return ns[fn.__name__], src, fn


def check_plants(verbose: bool = True) -> dict[str, list[dict]]:
    refs = {r.template_id: r for r in discover()}
    out: dict[str, list[dict]] = {}
    for branch in PREFIX:
        path = PLANTS_DIR / f'{branch}.py'
        if not path.exists():
            if verbose:
                print(f'{branch}: no plants file yet')
            continue
        mod = importlib.import_module(f'template_annotation_23092026.layer2.plants.{branch}')
        plants = list(mod.PLANTS)
        classes = [p['defect_class'] for p in plants]
        if sorted(classes) != ['arithmetic', 'constant', 'sign', 'unit']:
            raise SystemExit(f'{branch}: defect classes are {classes}, expected one each of '
                             f'arithmetic, constant, sign, unit')
        bases = [p['base'] for p in plants]
        if len(set(bases)) != len(bases):
            raise SystemExit(f'{branch}: plants must come from four different templates: {bases}')
        for p in plants:
            ref = refs.get(p['base'])
            if ref is None or ref.branch != branch:
                raise SystemExit(f"{p['plant_id']}: base {p['base']} is not a template of {branch}")
            fn_p, src, fn_b = mutated_function(ref, p)
            for s in range(50):
                seed_all(s)
                try:
                    qp, sp = fn_p()
                except Exception as exc:                            # noqa: BLE001
                    raise SystemExit(f"{p['plant_id']}: raised at seed {s}: {type(exc).__name__}: {exc}")
                seed_all(s)
                qb, sb = fn_b()
                if (qp, sp) == (qb, sb):
                    raise SystemExit(f"{p['plant_id']}: identical to its base at seed {s}; the defect is not live")
                seed_all(s)
                if fn_p() != (qp, sp):
                    raise SystemExit(f"{p['plant_id']}: not deterministic at seed {s}")
            p['_source'] = src
            p['_ref'] = ref
            if verbose:
                print(f"  ok {p['plant_id']:14s} {p['defect_class']:10s} from {p['base']}")
        out[branch] = plants
    # determinism across processes: one plant per branch, seed 7, in a child process
    for branch, plants in out.items():
        p = plants[0]
        code = ('import sys, random; sys.path.insert(0, %r); '
                'from template_annotation_23092026.layer2.build_tasks import mutated_function, discover; '
                'import template_annotation_23092026.layer2.plants.%s as m; '
                'ref=[r for r in discover() if r.template_id==%r][0]; fn,_,_=mutated_function(ref, m.PLANTS[0]); '
                'random.seed(7); print(hash(fn()))' % (str(REPO), branch, p['base']))
        child = subprocess.run([sys.executable, '-c', code], capture_output=True, text=True, cwd=REPO)
        fn_p, _, _ = mutated_function(p['_ref'], p)
        seed_all(7)
        here = hash(fn_p())
        # hash() of a str tuple is per-process salted; compare the texts instead
        code2 = code.replace('print(hash(fn()))', 'q,s=fn(); print(len(q), len(s), s[-60:])')
        child = subprocess.run([sys.executable, '-c', code2], capture_output=True, text=True, cwd=REPO)
        seed_all(7)
        q, s = fn_p()
        if child.stdout.strip() != f'{len(q)} {len(s)} {s[-60:]}':
            raise SystemExit(f"{p['plant_id']}: output differs between processes\n{child.stdout}\n{child.stderr}")
    return out


# --------------------------------------------------------------------- build

def build() -> None:
    plants = check_plants(verbose=True)
    missing = [b for b in PREFIX if b not in plants]
    if missing:
        raise SystemExit(f'plants missing for: {missing}')
    secret = secrets.token_hex(8)
    code_of = lambda item_id: 'T-' + hashlib.blake2b((secret + item_id).encode(), digest_size=4).hexdigest()  # noqa: E731

    pool, key = {}, []
    by_branch: dict[str, list[str]] = {b: [] for b in PREFIX}
    for ref in discover():
        inst = [generate(ref, s, capture=False) for s in INSTANCE_SEEDS]
        bad = [i for i in inst if not i.ok]
        if bad:
            raise SystemExit(f'{ref.template_id}: generation error {bad[0].error}')
        code = code_of(ref.template_id)
        pool[code] = {'branch': ref.branch, 'area': ref.area,
                      'source': inspect.getsource(ref.load()),
                      'instances': [{'seed': i.seed, 'question': i.question, 'solution': i.solution} for i in inst],
                      'gold_numbers': gold_numbers(inst[0].solution)}
        key.append({'code': code, 'kind': 'template', 'template_id': ref.template_id, 'branch': ref.branch})
        by_branch[ref.branch].append(code)
    for branch, plist in plants.items():
        for p in plist:
            fn_p, src, _ = mutated_function(p['_ref'], p)
            inst = []
            for s in PLANT_SEEDS:
                seed_all(s)
                q, sol = fn_p()
                inst.append({'seed': s, 'question': q, 'solution': sol})
            code = code_of(p['plant_id'])
            pool[code] = {'branch': branch, 'area': p['_ref'].area, 'source': src, 'instances': inst,
                          'gold_numbers': gold_numbers(inst[0]['solution'])}
            key.append({'code': code, 'kind': 'plant', 'plant_id': p['plant_id'], 'base': p['base'],
                        'defect_class': p['defect_class'], 'description': p['description'],
                        'detectable_by': p['detectable_by'], 'branch': branch})
            by_branch[branch].append(code)

    assignment = {}
    for a in roster():
        codes = list(by_branch[a['branch']])
        rng = random.Random(hashlib.blake2b((secret + a['id']).encode(), digest_size=8).digest())
        rng.shuffle(codes)
        assignment[a['id']] = {'branch': a['branch'], 'codes': codes}

    TASKS.mkdir(exist_ok=True)
    (TASKS / 'pool.json').write_text(json.dumps(pool, ensure_ascii=False, indent=0), encoding='utf8')
    (TASKS / 'assignment.json').write_text(json.dumps(assignment, indent=1), encoding='utf8')
    with (TASKS / 'keyfile.jsonl').open('w', encoding='utf8') as fh:
        for k in key:
            fh.write(json.dumps(k, ensure_ascii=False) + '\n')
    head = subprocess.run(['git', 'rev-parse', 'HEAD'], capture_output=True, text=True, cwd=REPO).stdout.strip()
    (TASKS / 'BUILD.json').write_text(json.dumps({
        'git_head': head, 'templates': sum(1 for k in key if k['kind'] == 'template'),
        'plants': sum(1 for k in key if k['kind'] == 'plant'), 'experts': len(assignment),
        'items_per_expert': {a: len(v['codes']) for a, v in assignment.items()},
        'instance_seeds': INSTANCE_SEEDS, 'plant_seeds': PLANT_SEEDS}, indent=1), encoding='utf8')
    n_t = sum(1 for k in key if k['kind'] == 'template')
    n_p = len(key) - n_t
    print(f'built tasks/: {n_t} templates + {n_p} plants for {len(assignment)} experts, '
          f'{len(next(iter(assignment.values()))["codes"])} items each; keyfile has {len(key)} rows and is never shipped')


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument('--check-plants', action='store_true')
    a = ap.parse_args()
    if a.check_plants:
        found = check_plants()
        print(f'{sum(len(v) for v in found.values())} plants verified in {len(found)} branches')
    else:
        build()


if __name__ == '__main__':
    main()
