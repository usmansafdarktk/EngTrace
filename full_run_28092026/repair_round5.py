"""Round 5: the re-drawn items of the two repaired chemical templates, carried through every store.

    python -m full_run_28092026.repair_round5 --list                  # the item ids the record covers, old and current sha256
    python -m full_run_28092026.repair_round5 --verify-manifest       # manifest.jsonl against its pre-round-5 copy
    python -m full_run_28092026.repair_round5 --archive [--variant V ...] [--dry-run]
    python -m full_run_28092026.repair_round5 --archive-paraphrases [--dry-run]
    python -m full_run_28092026.repair_round5 --paraphrase-kit [--expert che-1]
    python -m full_run_28092026.repair_round5 --paraphrase-score <folder of returned <id>.jsonl>
    python -m full_run_28092026.repair_round5 --rerun                 # BILLS: the re-run, one chain per model in parallel
    python -m full_run_28092026.repair_round5 --stages --variant V    # BILLS: judge then router on a re-scored store
    python -m full_run_28092026.repair_round5 --rescore --variant V [--dry-run]
    python -m full_run_28092026.repair_round5 --verify-rescore --variant V
    python -m full_run_28092026.repair_round5 --check                 # self-test on a temporary copy; touches nothing

Every mode but --rerun is free. --rerun bills: it runs the harness's own command for each model and store
(CHAINS below), after --archive and the harness's dry runs.

WHY. Layer 2 round 5 (template_annotation_23092026/layer2/fixes_round5.md) changed the questions of
work_isothermal_virial and adiabatic_flame_temperature, and round 6 capped the virial's organics; freeze.py
--only-templates drew their items again from the same seed. FREEZE.json's `round5` record lists every item id
either template held before or holds now, with its sha256 before round 5 (None for an id first frozen now) and
now (None for an id the selection no longer takes). A row written for another question is no longer a row of
the benchmark, and the harness resumes on item_id alone (run_traces.existing), so such rows leave the stores
before the current questions are run.

ARCHIVE. For every store (default: main and every arm that holds the subsample's items; --variant adds or
narrows), every row of a covered item that was not written for the item's CURRENT question moves out. A trace
row's question is named by its original_sha256 (an arm with a question of its own, paraphrase or openbook) or
else its item_sha256, and it is current when that hash is the item's sha256 in manifest.jsonl; a row with no hash
is stale, since every row the harness writes carries item_sha256, and every row of an id the selection no
longer takes is stale. Trace rows go to traces/_replaced_round5/<variant>/<model>.jsonl; score rows and stage
rows (e5/, router/, and main's e5_grok-4-6/) go to scores/_replaced/round5_<utc>/<variant>/... A score or stage
row carries no item hash, so it moves only for an item that has no current trace row: once an item has been run
again, its rows are current. Every move is logged per store, model, item and old sha256 in MOVED.jsonl of that
run's archive; kept lines are written back byte for byte through a temporary file, and a file that changed
while it was read is left alone and reported. A second run moves nothing.

PARAPHRASES. --archive-paraphrases moves the covered items' writer attempts that were not written for the
current question out of paraphrase/attempts.jsonl, and their stale verdicts out of paraphrase/accepted.json
(a verdict is current only when it names the question it judged, as --paraphrase-score records), into the same
archive; then it rebuilds paraphrase/manifest.jsonl, pool.jsonl and PARAPHRASE.md with paraphrase.rebuild(),
so the arm's manifest agrees with the re-frozen manifest (run_traces refuses a paraphrase or openbook item
whose original_sha256 is not the manifest's) and the items are unwritten, to be written again from the current
questions (paraphrase.py). An item with no verdict is not kept by the analysis (analyze.accepted_pairs).
--paraphrase-kit builds the expert check of the current passing pairs for one chemical expert with
paraphrase_kit.build, into paraphrase/dist/round5/ (local); --paraphrase-score scores the returned file there
and merges the verdicts into paraphrase/accepted.json, each naming the original and paraphrase it judged.

RE-SCORE (A9). score.py ties a store to the hashes of manifest.jsonl and diversity.json and refuses a store
scored with other inputs; with --replace it archives the whole store, e5/ and router/ included. --rescore
therefore (1) refuses while the variant's traces hold a stale row; (2) moves every stage folder of the store
(e5/, router/, e5_<judge>/) aside into the round-5 archive; (3) runs score.py --variant V --replace; (4) moves
the stage folders back, whatever (3) did; (5) runs --verify-rescore. The stage rows of unchanged items stay
valid because their inputs are unchanged; judge.py and router.py rebuild every row from their reply stores and
bill only the new items. Each re-score is logged in scores/_replaced/round5_rescore.jsonl, so WS-G can repeat
it: the same command, after the evaluator is final.

VERIFY-RESCORE. In every model file of the re-scored store, every row of an item the record does not cover
equals its row in the store score.py archived (logged by --rescore), field by field, and both hold the same
such items; the store's CONFIG.json names the current manifest and diversity.json. Rows carry no provenance
field of their own, so equality is exact; where an evaluator file differs between the two CONFIG.json files,
the fields that moved are listed so the cause is visible.
"""
from __future__ import annotations

import argparse
import collections
import datetime as dt
import hashlib
import json
import os
import shutil
import subprocess
import sys
import tempfile
from pathlib import Path

HERE = Path(__file__).resolve().parent
REPO = HERE.parent
if str(REPO) not in sys.path:
    sys.path.insert(0, str(REPO))

RECORD = 'round5'
DEFAULT_VARIANTS = ('main', 'reasoning-medium', 'flagship', 'flagship-reasoning-medium', 'openbook2', 'openbook',
                    'tool', 'paraphrase', 'repeat1', 'repeat2', 'repeat3')
CHEMICAL = 'chemical_engineering'


class Paths:
    """Where the stores live; the self-test points one at a temporary copy."""

    def __init__(self, root: Path):
        self.root = root
        self.freeze = root / 'FREEZE.json'
        self.manifest = root / 'manifest.jsonl'
        self.traces = root / 'traces'
        self.scores = root / 'scores'
        self.trace_archive = root / 'traces' / '_replaced_round5'
        self.rescore_log = root / 'scores' / '_replaced' / 'round5_rescore.jsonl'
        self.paraphrase = root / 'paraphrase'


P = Paths(HERE)


def now() -> str:
    return dt.datetime.now(dt.timezone.utc).strftime('%Y%m%dT%H%M%SZ')


def sha(b: bytes) -> str:
    return hashlib.sha256(b).hexdigest()


# ------------------------------------------------------------------ the items

def record(p: Paths = P) -> dict:
    doc = json.loads(p.freeze.read_text(encoding='utf-8'))
    if RECORD not in doc:
        raise SystemExit(f'FREEZE.json has no {RECORD} record: run freeze.py --only-templates ... --write --record {RECORD}')
    return doc[RECORD]


def shas(p: Paths = P) -> tuple[dict[str, str | None], dict[str, str | None]]:
    """({item_id: sha256 before round 5, None for an id first frozen since},
        {item_id: sha256 in manifest.jsonl now, None for an id the selection no longer takes})
    over every item id the record covers."""
    items = record(p)['items']
    man = {json.loads(x)['item_id']: json.loads(x)['sha256'] for x in p.manifest.read_text(encoding='utf-8').splitlines()}
    return {x['item_id']: x['old_sha256'] for x in items}, {x['item_id']: man.get(x['item_id']) for x in items}


def list_items() -> int:
    old, cur = shas()
    new = {x['item_id']: x['new_sha256'] for x in record()['items']}
    for i in old:
        flag = '' if cur[i] == new[i] else '   <- the manifest disagrees with the record'
        state = 'left' if cur[i] is None else 'new id' if old[i] is None else 'changed'
        print(f'{i:34s} {state:8s} old {(old[i] or "-")[:12]:12s}  now {(cur[i] or "-")[:12]}{flag}')
    print(f'{len(old)} item ids, {sum(v is not None for v in cur.values())} frozen now; '
          f'templates {", ".join(record()["templates"])}')
    return 0


def verify_manifest() -> int:
    """The manifest against the copy freeze.py archived at the record's first write: every row outside the
    record's templates identical; inside, every row as the record says."""
    rec = record()
    first = (rec.get('refreezes') or [rec])[0]['archive']
    old_rows = [json.loads(x) for x in (HERE / first / 'manifest.jsonl').read_text(encoding='utf-8').splitlines()]
    body = P.manifest.read_text(encoding='utf-8')
    new_rows = [json.loads(x) for x in body.splitlines()]
    old_sha, _cur = shas()
    new_sha = {x['item_id']: x['new_sha256'] for x in rec['items']}
    templates = set(rec['templates'])
    bad = []
    keep_o = [r for r in old_rows if r['template_id'] not in templates]
    keep_n = [r for r in new_rows if r['template_id'] not in templates]
    if keep_o != keep_n:
        bad.append("a row outside the record's templates changed")
    if [r['template_id'] for r in old_rows] != [r['template_id'] for r in new_rows]:
        bad.append('the template blocks moved')
    o_t = {r['item_id']: r for r in old_rows if r['template_id'] in templates}
    n_t = {r['item_id']: r for r in new_rows if r['template_id'] in templates}
    for i in sorted(set(o_t) | set(n_t)):
        if o_t.get(i, {}).get('sha256') != old_sha.get(i) or n_t.get(i, {}).get('sha256') != new_sha.get(i):
            bad.append(f'{i}: sha256 not as recorded')
        if i in o_t and i in n_t and {k for k in set(o_t[i]) | set(n_t[i]) if o_t[i].get(k) != n_t[i].get(k)} - {'sha256'}:
            bad.append(f'{i}: a field other than sha256 differs')
    if json.loads(P.freeze.read_text(encoding='utf-8'))['manifest_sha256'] != sha(body.encode('utf-8')):
        bad.append("FREEZE.json manifest_sha256 is not the manifest's")
    changed = sum(1 for i in set(o_t) & set(n_t) if o_t[i] != n_t[i])
    print(f'manifest: {len(new_rows)} rows; {len(keep_n)} rows outside the record\'s {len(templates)} templates identical '
          f'to the pre-round-5 copy; inside them {changed} rows changed in sha256 only, {len(set(o_t) - set(n_t))} ids '
          f'left and {len(set(n_t) - set(o_t))} came in, each as the record says' + (' - FAILED' if bad else ''))
    for b in bad[:20]:
        print('  FAIL', b)
    return 1 if bad else 0


# ------------------------------------------------------------------ archive

def trace_files(p: Paths, variant: str) -> list[Path]:
    d = p.traces if variant == 'main' else p.traces / variant
    return sorted(d.glob('*.jsonl')) if d.exists() else []


def score_files(p: Paths, variant: str) -> list[tuple[str, Path]]:
    """(kind, path): the model files and every stage folder's model files of one store."""
    d = p.scores / variant
    if not d.exists():
        return []
    out = [('score', f) for f in sorted(d.glob('*.jsonl'))]
    for sub in sorted(x for x in d.iterdir() if x.is_dir()):
        out += [(sub.name, f) for f in sorted(sub.glob('*.jsonl'))]
    return out


def classify(raw: bytes, cur: dict) -> str:
    """'keep' or 'old' for one trace line: a row of a covered item is old unless written for its current question."""
    try:
        r = json.loads(raw)
    except (json.JSONDecodeError, UnicodeDecodeError):
        return 'keep'
    i = r.get('item_id') if isinstance(r, dict) else None
    if i not in cur:
        return 'keep'
    h = r.get('original_sha256') or r.get('item_sha256')
    return 'keep' if cur[i] is not None and h == cur[i] else 'old'


def split_file(path: Path, decide) -> tuple[bytes, list[bytes], list[bytes]]:
    """(original bytes, kept lines, moved lines). `decide(raw) -> 'keep' | 'old'`; blank lines are kept."""
    data = path.read_bytes()
    lines = data.split(b'\n')
    if lines and lines[-1] == b'':
        lines.pop()
    keep, moved = [], []
    for raw in lines:
        if raw.strip() and decide(raw.rstrip(b'\r')) == 'old':
            moved.append(raw)
        else:
            keep.append(raw)
    return data, keep, moved


def rewrite(path: Path, data: bytes, keep: list[bytes]) -> bool:
    """Write the kept lines back through a temporary file, unless the file changed since it was read."""
    body = b'\n'.join(keep) + (b'\n' if data.endswith(b'\n') and keep else b'')
    if path.read_bytes() != data:
        return False
    tmp = path.with_name(path.name + '.round5.tmp')
    tmp.write_bytes(body)
    os.replace(tmp, path)
    return True


def item_of(raw: bytes) -> str | None:
    try:
        r = json.loads(raw)
    except (json.JSONDecodeError, UnicodeDecodeError):
        return None
    return r.get('item_id') if isinstance(r, dict) else None


def archive(variants: list[str], dry: bool, p: Paths = P, quiet: bool = False) -> dict:
    old, cur = shas(p)
    run_dir = p.scores / '_replaced' / f'{RECORD}_{now()}'
    log, totals, problems = [], collections.Counter(), []

    def note(variant, model, kind, path, moved_lines):
        per = collections.defaultdict(list)
        for raw in moved_lines:
            r = json.loads(raw.rstrip(b'\r'))
            per[r.get('item_id')].append(r.get('status'))
        for i, st in sorted(per.items()):
            log.append({'variant': variant, 'model': model, 'kind': kind,
                        'path': path.relative_to(p.root).as_posix(), 'item_id': i, 'old_sha256': old.get(i),
                        'current_sha256': cur.get(i), 'rows': len(st),
                        'statuses': dict(collections.Counter(s or '-' for s in st))})
        totals[f'{variant}:{kind}'] += len(moved_lines)

    for variant in variants:
        current = collections.defaultdict(set)              # model -> covered items with a current trace row
        for f in trace_files(p, variant):
            model = f.stem
            data, keep, moved = split_file(f, lambda raw: classify(raw, cur))
            for raw in keep:
                i = item_of(raw.rstrip(b'\r'))
                if i in cur:
                    current[model].add(i)
            if not moved:
                continue
            note(variant, model, 'trace', f, moved)
            if dry:
                continue
            dest = p.trace_archive / variant / f.name
            dest.parent.mkdir(parents=True, exist_ok=True)
            with open(dest, 'ab') as fh:
                fh.write(b''.join(raw + b'\n' for raw in moved))
            if not rewrite(f, data, keep):
                problems.append(f'{f.relative_to(p.root).as_posix()}: changed while being read; not rewritten '
                                f'(its {len(moved)} stale rows were copied to the archive; run again)')
        for kind, f in score_files(p, variant):
            done = current.get(f.stem, set())

            def decide(raw, done=done):
                i = item_of(raw)
                return 'old' if (i in cur and i not in done) else 'keep'
            data, keep, moved = split_file(f, decide)
            if not moved:
                continue
            note(variant, f.stem, kind, f, moved)
            if dry:
                continue
            dest = run_dir / variant / ('' if kind == 'score' else kind) / f.name
            dest.parent.mkdir(parents=True, exist_ok=True)
            with open(dest, 'ab') as fh:
                fh.write(b''.join(raw + b'\n' for raw in moved))
            if not rewrite(f, data, keep):
                problems.append(f'{f.relative_to(p.root).as_posix()}: changed while being read; not rewritten')
    if log and not dry:
        run_dir.mkdir(parents=True, exist_ok=True)
        with open(run_dir / 'MOVED.jsonl', 'a', encoding='utf-8', newline='\n') as fh:
            for x in log:
                fh.write(json.dumps(x) + '\n')
    if not quiet:
        what = 'would move' if dry else 'moved'
        print(f'{what} {sum(totals.values())} rows' + (f'; archive {run_dir.relative_to(p.root).as_posix()}'
                                                         if log and not dry else ''))
        for k, n in sorted(totals.items()):
            print(f'  {k:42s} {n}')
        for x in problems:
            print('  NOTE', x)
    return {'moved': dict(totals), 'log': log, 'problems': problems, 'run_dir': run_dir if log and not dry else None}


# ------------------------------------------------------------------ paraphrases

def archive_paraphrases(dry: bool, p: Paths = P, rebuild: bool = True) -> dict:
    _old, cur = shas(p)
    att = p.paraphrase / 'attempts.jsonl'
    acc_path = p.paraphrase / 'accepted.json'
    run_dir = p.scores / '_replaced' / f'{RECORD}_{now()}' / 'paraphrase'

    def decide(raw):
        try:
            r = json.loads(raw)
        except (json.JSONDecodeError, UnicodeDecodeError):
            return 'keep'
        i = r.get('item_id')
        if i not in cur:
            return 'keep'
        return 'keep' if cur[i] is not None and r.get('original_sha256') == cur[i] else 'old'
    data, keep, moved = split_file(att, decide) if att.exists() else (b'', [], [])
    acc = json.loads(acc_path.read_text(encoding='utf-8')) if acc_path.exists() else {}
    acc_moved = {i: v for i, v in acc.items()
                 if i in cur and not (cur[i] is not None and v.get('original_sha256') == cur[i])}
    print(f"paraphrase attempts {'to move' if dry else 'moved'}: {len(moved)} rows of "
          f"{len({item_of(r.rstrip(b'\r')) for r in moved})} items; verdicts: {len(acc_moved)}")
    if dry or (not moved and not acc_moved):
        return {'attempts': len(moved), 'verdicts': len(acc_moved)}
    run_dir.mkdir(parents=True, exist_ok=True)
    for name in ('manifest.jsonl', 'pool.jsonl'):
        src = p.paraphrase / name
        if src.exists():
            shutil.copy2(src, run_dir / name)
    if moved:
        with open(run_dir / 'attempts.jsonl', 'ab') as fh:
            fh.write(b''.join(raw + b'\n' for raw in moved))
        if not rewrite(att, data, keep):
            raise SystemExit('paraphrase/attempts.jsonl changed while being read; nothing rewritten, run again')
    if acc_moved:
        (run_dir / 'accepted.json').write_text(json.dumps(acc_moved, indent=1) + '\n', encoding='utf-8')
        rest = {i: v for i, v in acc.items() if i not in acc_moved}
        acc_path.write_text(json.dumps(rest, indent=1) + '\n', encoding='utf-8')
    if rebuild:
        from full_run_28092026 import paraphrase
        paraphrase.rebuild()
    return {'attempts': len(moved), 'verdicts': len(acc_moved), 'archive': run_dir}


def current_paraphrases() -> dict[str, tuple[str, str, str]]:
    """{item_id: (paraphrase text, its sha256, the original's sha256)} for the covered items' current passing pairs."""
    _old, cur = shas()
    man = {r['item_id']: r for r in map(json.loads, (P.paraphrase / 'manifest.jsonl').read_text(encoding='utf-8').splitlines())}
    out = {}
    for r in map(json.loads, (P.paraphrase / 'pool.jsonl').read_text(encoding='utf-8').splitlines()):
        m = man.get(r['item_id'])
        if (r['item_id'] in cur and cur[r['item_id']] is not None and m and m['passed']
                and m['original_sha256'] == cur[r['item_id']] and sha(r['question'].encode('utf-8')) == m['sha256']):
            out[r['item_id']] = (r['question'], m['sha256'], m['original_sha256'])
    return out


def paraphrase_kit(expert: str) -> int:
    from full_run_28092026 import paraphrase_kit as pk, score
    paras = current_paraphrases()
    if not paras:
        raise SystemExit('no current passing paraphrase of a covered item: write them first (paraphrase.py)')
    roster = {e['id']: e for e in pk.roster()}
    if roster.get(expert, {}).get('branch') != CHEMICAL:
        raise SystemExit(f'{expert} is not a chemical expert of the layer-2 roster')
    out = P.paraphrase / 'dist' / RECORD
    info = pk.build(out, {i: v[0] for i, v in paras.items()}, score.pool_items(), [roster[expert]])
    (out / 'judged.json').write_text(json.dumps({i: {'sha256': v[1], 'original_sha256': v[2]} for i, v in paras.items()},
                                                indent=1) + '\n', encoding='utf-8')
    print(f"{info['items']} paraphrases for {expert}; kit in {(out / 'dist').relative_to(REPO).as_posix()}: send "
          f"app.py, guide.md, README.txt and kit_{expert}/; the return goes to paraphrase/returned/{RECORD}/")
    return 0


def paraphrase_score(returned: Path) -> int:
    from full_run_28092026 import paraphrase_kit as pk
    out = P.paraphrase / 'dist' / RECORD
    judged = json.loads((out / 'judged.json').read_text(encoding='utf-8'))
    now_pairs = current_paraphrases()
    stale = [i for i, v in judged.items() if now_pairs.get(i, (None, None, None))[1:] != (v['sha256'], v['original_sha256'])]
    if stale:
        raise SystemExit(f'the kit judged pairs that are no longer current: {stale}; rebuild the kit')
    res = pk.score(returned, out, None)
    new = json.loads((out / 'accepted.json').read_text(encoding='utf-8'))
    acc_path = P.paraphrase / 'accepted.json'
    acc = json.loads(acc_path.read_text(encoding='utf-8'))
    clash = [i for i in new if i in acc]
    if clash:
        raise SystemExit(f'accepted.json already holds {clash}: archive the stale verdicts first (--archive-paraphrases)')
    acc.update({i: {**v, **judged[i]} for i, v in new.items()})
    acc_path.write_text(json.dumps(acc, indent=1) + '\n', encoding='utf-8')
    kept = sum(1 for v in acc.values() if v.get('kept') is not False)
    print(f"round 5: assigned {res['assigned']}, returned {res['returned']}, kept {res['kept']}; accepted.json now "
          f"{len(acc)} verdicts, {kept} pairs kept or outstanding")
    return 0


# ------------------------------------------------------------------ re-run (A8, BILLS)

# One chain per model, the chains in parallel; within a chain the main store first (the critical path), then
# the model's arms. Each model writes only its own trace file per store, so no two steps can collide, and no
# model is called from two stores at once. The tool arm runs the two models it holds (gpt-oss-20b is not
# servable there, D-184); the anchors run only in their arms.
CHAINS = {
    'gpt-oss-20b': ['main', 'openbook2', 'repeat1', 'repeat2', 'repeat3', 'paraphrase'],
    'gemma-4-26b-a4b': ['main', 'repeat1', 'repeat2', 'repeat3', 'paraphrase'],
    'deepseek-v4.1-flash': ['main', 'paraphrase'],
    'qwen3-235b-a22b-2507': ['main', 'repeat1', 'repeat2', 'repeat3', 'paraphrase'],
    'glm-5.3-flash': ['main', 'paraphrase'],
    'glm-5.3': ['main', 'paraphrase'],
    'muse-glimmer-30b': ['main', 'paraphrase'],
    'kimi-k3': ['main', 'paraphrase'],
    'gpt-5.4-mini': ['main', 'reasoning-medium', 'openbook2', 'tool', 'paraphrase'],
    'gemini-3.1-flash-lite': ['main', 'reasoning-medium', 'repeat1', 'repeat2', 'repeat3', 'paraphrase'],
    'claude-sonnet-5': ['main', 'openbook2', 'tool', 'paraphrase'],
    'deepseek-v4-pro': ['flagship'],
    'gpt-5.4': ['flagship', 'flagship-reasoning-medium'],
}
PASSES = 2          # the harness resumes; a second pass calls only service failures and unrun items


def rerun() -> int:
    """BILLS. Every chain's steps through the harness's own command (run_traces --variant V --model M --yes), each
    run PASSES times; logs per chain in traces/_round5_rerun/, and DONE.json when every chain has ended."""
    import concurrent.futures as cf
    log_dir = P.traces / '_round5_rerun'
    log_dir.mkdir(parents=True, exist_ok=True)

    def chain(model: str, steps: list[str]):
        out = []
        with open(log_dir / f'{model}.log', 'a', encoding='utf-8') as fh:
            for variant in steps:
                for k in range(1, PASSES + 1):
                    cmd = [sys.executable, '-m', 'full_run_28092026.run_traces', '--variant', variant, '--model', model, '--yes']
                    fh.write(f'=== {now()} {" ".join(cmd[1:])} (pass {k})\n')
                    fh.flush()
                    rc = subprocess.run(cmd, cwd=REPO, stdout=fh, stderr=subprocess.STDOUT).returncode
                    fh.write(f'=== {now()} exit {rc}\n')
                    fh.flush()
                    out.append({'variant': variant, 'pass': k, 'exit': rc})
        return model, out

    print(f'{now()} re-run: {len(CHAINS)} chains, {sum(len(s) for s in CHAINS.values())} steps', flush=True)
    with cf.ThreadPoolExecutor(len(CHAINS)) as ex:
        results = dict(ex.map(lambda kv: chain(*kv), CHAINS.items()))
    (log_dir / 'DONE.json').write_text(json.dumps({'finished_at_utc': now(), 'chains': results}, indent=1) + '\n',
                                       encoding='utf-8')
    bad = {m: [s for s in st if s['exit'] != 0] for m, st in results.items()}
    print(f'{now()} re-run finished; steps with a non-zero exit: {sum(len(v) for v in bad.values())}', flush=True)
    return 0


def stages(variant: str, e5_budget: float, router_budget: float) -> int:
    """BILLS. The judge (E5) and then the step-check (the router, where the store has one) on one re-scored store,
    each resuming on its reply store, so only the new items' jobs are sent. One after the other: both call the
    same judge model, which refuses more than a few calls at once, and each keeps one reply file. Each cap is the
    reply store's running total plus the budget given, since --max-usd counts the whole store (D-148)."""
    from full_run_28092026.judge_calls import Store
    from full_run_28092026 import judge, router
    runs = [('judge', judge.REPLIES, e5_budget, '4')]
    if (P.scores / variant / 'router').exists():
        runs.append(('router', router.REPLIES, router_budget, '8'))
    for mod, replies, budget, workers in runs:
        cap = round(Store(replies).summary()['billed_usd_all_lines'] + budget, 2)
        cmd = [sys.executable, '-m', f'full_run_28092026.{mod}', '--variant', variant, '--yes',
               '--max-usd', str(cap), '--workers', workers]
        print(f'{now()} {" ".join(cmd[1:])}', flush=True)
        rc = subprocess.run(cmd, cwd=REPO).returncode
        print(f'{now()} {mod} on {variant}: exit {rc}', flush=True)
        if rc:
            return rc
    return 0


# ------------------------------------------------------------------ re-score

def stage_dirs(store: Path) -> list[Path]:
    return sorted(x for x in store.iterdir() if x.is_dir()) if store.exists() else []


def rescore(variant: str, dry: bool) -> int:
    _old, cur = shas()
    left = [f.name for f in trace_files(P, variant)
            if any(classify(raw.rstrip(b'\r'), cur) == 'old' for raw in f.read_bytes().split(b'\n') if raw.strip())]
    if left:
        raise SystemExit(f'{variant}: traces still hold stale rows ({", ".join(left)}); run --archive first')
    store = P.scores / variant
    if not (store / 'CONFIG.json').exists():
        raise SystemExit(f'scores/{variant} has no CONFIG.json: nothing to re-score')
    stages = stage_dirs(store)
    side = P.scores / '_replaced' / f'{RECORD}_{now()}' / f'{variant}_stages'
    cmd = [sys.executable, '-m', 'full_run_28092026.score', '--variant', variant, '--replace']
    print(f'{variant}: stage folders aside {[s.name for s in stages]} -> {side.relative_to(HERE).as_posix()}; '
          f'then {" ".join(cmd[1:])}')
    if dry:
        return 0
    side.mkdir(parents=True, exist_ok=True)
    for s in stages:
        shutil.move(str(s), str(side / s.name))
    started = now()
    proc = None
    try:
        proc = subprocess.run(cmd, cwd=REPO, capture_output=True, text=True, encoding='utf-8')
        sys.stdout.write(proc.stdout)
        sys.stderr.write(proc.stderr)
    finally:
        for s in stage_dirs(side):
            dest = store / s.name
            if dest.exists():
                raise SystemExit(f'{dest} exists after the re-score; the stage rows stay in {s}: move them by hand')
            store.mkdir(parents=True, exist_ok=True)
            shutil.move(str(s), str(dest))
    archived = None
    for line in proc.stdout.splitlines():
        if line.startswith('archived scores') and ' -> ' in line and line.split(' -> ')[0].rstrip().endswith(variant):
            archived = line.split(' -> ', 1)[1].strip().replace('\\', '/')
    entry = {'variant': variant, 'rescored_at_utc': started, 'returncode': proc.returncode,
             'store_archived_to': archived, 'stages_moved': [s.name for s in stages]}
    P.rescore_log.parent.mkdir(parents=True, exist_ok=True)
    with open(P.rescore_log, 'a', encoding='utf-8', newline='\n') as fh:
        fh.write(json.dumps(entry) + '\n')
    if proc.returncode != 0:
        print(f'{variant}: score.py exited {proc.returncode}; the stage folders are back in place')
        return proc.returncode
    return verify_rescore(variant)


def compare_rows(old_rows: dict, new_rows: dict, skip: set) -> tuple[int, int, collections.Counter, list]:
    """(rows compared, rows equal, fields that differ, item ids present in only one side) outside `skip`.

    A field absent on one side and None on the other is the same value: score.py added the tool arm's fields
    (D-184) to every row, None outside that arm, so a store scored before then lacks them. Any other change,
    a value that becomes None included, is a difference."""
    a = {i for i in old_rows if i not in skip}
    b = {i for i in new_rows if i not in skip}
    fields, equal = collections.Counter(), 0
    for i in sorted(a & b):
        o, n = old_rows[i], new_rows[i]
        diff = [k for k in set(o) | set(n) if o.get(k) != n.get(k)]
        if not diff:
            equal += 1
        else:
            fields.update(diff)
    return len(a & b), equal, fields, sorted(a ^ b)


def format_only_keys(old_rows: dict, new_rows: dict) -> list[str]:
    """Keys present in one side's rows only, with None wherever they are present: a format change."""
    ko = {k for r in old_rows.values() for k in r}
    kn = {k for r in new_rows.values() for k in r}
    return sorted(k for k in ko ^ kn
                  if all(r.get(k) is None for r in list(old_rows.values()) + list(new_rows.values())))


def load_rows(path: Path) -> dict:
    if not path.exists():
        return {}
    return {r['item_id']: r for r in (json.loads(x) for x in path.read_text(encoding='utf-8').splitlines() if x.strip())}


def verify_rescore(variant: str) -> int:
    if not P.rescore_log.exists():
        raise SystemExit('no re-score logged: run --rescore first')
    entries = [json.loads(x) for x in P.rescore_log.read_text(encoding='utf-8').splitlines() if x.strip()]
    mine = [e for e in entries if e['variant'] == variant and e.get('store_archived_to')]
    if not mine:
        raise SystemExit(f'no logged re-score of {variant} archived a store')
    prev = HERE / mine[-1]['store_archived_to']
    store = P.scores / variant
    old, _cur = shas()
    bad = []
    cfg_new = json.loads((store / 'CONFIG.json').read_text(encoding='utf-8'))
    cfg_old = json.loads((prev / 'CONFIG.json').read_text(encoding='utf-8')) if (prev / 'CONFIG.json').exists() else {}
    for rel in ('full_run_28092026/manifest.jsonl', 'full_run_28092026/diversity.json'):
        if cfg_new.get('inputs', {}).get(rel) != sha((REPO / rel).read_bytes()):
            bad.append(f'CONFIG.json does not name the current {rel}')
    moved_eval = sorted(k for k in set(cfg_new.get('evaluators', {})) | set(cfg_old.get('evaluators', {}))
                        if cfg_new.get('evaluators', {}).get(k, {}).get('sha256_lf')
                        != cfg_old.get('evaluators', {}).get(k, {}).get('sha256_lf'))
    total = same = 0
    for f in sorted(store.glob('*.jsonl')):
        n_rows, o_rows = load_rows(f), load_rows(prev / f.name)
        if not o_rows:
            bad.append(f'{f.name}: no archived rows to compare with')
            continue
        n, eq, fields, only = compare_rows(o_rows, n_rows, set(old))
        total += n
        same += eq
        covered = sum(1 for i in n_rows if i in old)
        fmt = format_only_keys(o_rows, n_rows)
        print(f'  {variant}/{f.stem}: {eq} of {n} rows of uncovered items equal; covered items in the store: {covered}'
              + (f'; differing fields {dict(fields)}' if fields else '') + (f'; in one side only {only[:5]}' if only else '')
              + (f'; keys added as None (format only): {fmt}' if fmt else ''))
        if eq != n or only:
            bad.append(f'{f.name}: {n - eq} rows differ, {len(only)} items in one side only')
    if moved_eval:
        print(f'  evaluator files that differ between the two CONFIG.json: {moved_eval}')
    print(f'{variant}: {same} of {total} rows of uncovered items equal the archived store '
          f'({prev.relative_to(HERE).as_posix()})' + (' - FAILED' if bad else ''))
    for b in bad[:20]:
        print('  FAIL', b)
    return 1 if bad else 0


# ------------------------------------------------------------------ self-test

def self_test() -> int:
    """Archive, idempotence, the stage rule, byte preservation, an id that left the selection, the paraphrase move
    and the row comparison, on a temporary tree; nothing real is read or written."""
    bad = []
    with tempfile.TemporaryDirectory() as tmp:
        root = Path(tmp)
        p = Paths(root)
        # t_a#0 changed, t_a#1 left the selection, t_a#5 came in
        items = [{'item_id': 't_a#0', 'old_sha256': 'old0', 'new_sha256': 'cur0'},
                 {'item_id': 't_a#5', 'old_sha256': None, 'new_sha256': 'cur5'},
                 {'item_id': 't_a#1', 'old_sha256': 'old1', 'new_sha256': None}]
        p.freeze.write_text(json.dumps({RECORD: {'templates': ['template_t_a'], 'archive': 'x', 'items': items}}),
                            encoding='utf-8')
        p.manifest.write_text('\n'.join(json.dumps({'item_id': i, 'sha256': h}) for i, h in
                                        (('other#3', 'zzz'), ('t_a#0', 'cur0'), ('t_a#5', 'cur5'))) + '\n', encoding='utf-8')
        (p.traces / 'paraphrase').mkdir(parents=True)
        (p.scores / 'main' / 'e5').mkdir(parents=True)
        (p.scores / 'paraphrase').mkdir(parents=True)
        main_lines = [
            json.dumps({'item_id': 'other#3', 'item_sha256': 'zzz', 'status': 'answered'}),
            json.dumps({'item_id': 't_a#0', 'item_sha256': 'old0', 'status': 'answered'}),
            json.dumps({'item_id': 't_a#0', 'item_sha256': 'interim', 'status': 'error'}),
            json.dumps({'item_id': 't_a#1', 'item_sha256': 'old1', 'status': 'answered'}),
            '{not json',
            json.dumps({'item_id': 't_a#0', 'status': 'error'}),
            json.dumps({'item_id': 'other#4', 'item_sha256': 'yyy', 'status': 'empty', 'text': 'café ≤ 1'}),
        ]
        tf = p.traces / 'm1.jsonl'
        tf.write_bytes(('\n'.join(main_lines) + '\n').encode('utf-8'))
        sf = p.scores / 'main' / 'm1.jsonl'
        sf.write_text('\n'.join(json.dumps({'item_id': i}) for i in ['other#3', 't_a#0', 't_a#1']) + '\n', encoding='utf-8')
        ef = p.scores / 'main' / 'e5' / 'm1.jsonl'
        ef.write_text('\n'.join(json.dumps({'item_id': i}) for i in ['other#3', 't_a#0', 't_a#1']) + '\n', encoding='utf-8')
        pf = p.traces / 'paraphrase' / 'm1.jsonl'
        pf.write_text(json.dumps({'item_id': 't_a#0', 'item_sha256': 'para', 'original_sha256': 'old0', 'status': 'answered'})
                      + '\n' + json.dumps({'item_id': 't_a#5', 'item_sha256': 'para5', 'original_sha256': 'cur5',
                                           'status': 'answered'}) + '\n', encoding='utf-8')
        (p.scores / 'paraphrase' / 'm1.jsonl').write_text(
            json.dumps({'item_id': 't_a#0'}) + '\n' + json.dumps({'item_id': 't_a#5'}) + '\n', encoding='utf-8')
        r1 = archive(['main', 'paraphrase'], dry=False, p=p, quiet=True)
        kept = tf.read_bytes().decode('utf-8').splitlines()
        if kept != [main_lines[0], main_lines[4], main_lines[6]]:
            bad.append(f'main trace kept {kept}')
        if r1['moved'].get('main:trace') != 4:
            bad.append(f"main trace moved {r1['moved'].get('main:trace')}, not 4")
        if [json.loads(x)['item_id'] for x in sf.read_text(encoding='utf-8').splitlines()] != ['other#3']:
            bad.append('score rows of stale or departed items stayed')
        if [json.loads(x)['item_id'] for x in ef.read_text(encoding='utf-8').splitlines()] != ['other#3']:
            bad.append('stage rows of stale or departed items stayed')
        if [json.loads(x)['item_id'] for x in pf.read_text(encoding='utf-8').splitlines()] != ['t_a#5']:
            bad.append('paraphrase trace: original_sha256 was not read first')
        if [json.loads(x)['item_id'] for x in (p.scores / 'paraphrase' / 'm1.jsonl').read_text(encoding='utf-8').splitlines()] != ['t_a#5']:
            bad.append('a score row of an item with a current trace row was moved')
        if not (p.trace_archive / 'main' / 'm1.jsonl').exists() or not list((p.scores / '_replaced').glob('round5_*/MOVED.jsonl')):
            bad.append('archive or MOVED.jsonl missing')
        r2 = archive(['main', 'paraphrase'], dry=False, p=p, quiet=True)
        if sum(r2['moved'].values()) != 0:
            bad.append(f"second run moved {r2['moved']}")
        with open(tf, 'ab') as fh:                          # the re-run writes a current row and its stage row
            fh.write((json.dumps({'item_id': 't_a#0', 'item_sha256': 'cur0', 'status': 'answered'}) + '\n').encode())
        with open(ef, 'a', encoding='utf-8') as fh:
            fh.write(json.dumps({'item_id': 't_a#0'}) + '\n')
        r3 = archive(['main'], dry=False, p=p, quiet=True)
        if sum(r3['moved'].values()) != 0 or 't_a#0' not in ef.read_text(encoding='utf-8'):
            bad.append('a current row or its stage row was moved')
        p.paraphrase.mkdir()
        (p.paraphrase / 'attempts.jsonl').write_text(
            json.dumps({'item_id': 'other#3', 'original_sha256': 'zzz'}) + '\n'
            + json.dumps({'item_id': 't_a#0', 'original_sha256': 'old0'}) + '\n'
            + json.dumps({'item_id': 't_a#0', 'status': 'service_failure'}) + '\n'
            + json.dumps({'item_id': 't_a#0', 'original_sha256': 'cur0'}) + '\n', encoding='utf-8')
        (p.paraphrase / 'accepted.json').write_text(json.dumps({'other#3': {'kept': True}, 't_a#0': {'kept': True},
                                                                't_a#5': {'kept': True, 'original_sha256': 'cur5'}}),
                                                     encoding='utf-8')
        rp = archive_paraphrases(False, p=p, rebuild=False)
        left = [json.loads(x).get('original_sha256') for x in (p.paraphrase / 'attempts.jsonl').read_text(encoding='utf-8').splitlines()]
        if rp['attempts'] != 2 or left != ['zzz', 'cur0']:
            bad.append(f'paraphrase attempts left {left}')
        if set(json.loads((p.paraphrase / 'accepted.json').read_text(encoding='utf-8'))) != {'other#3', 't_a#5'}:
            bad.append('stale verdicts stayed in accepted.json, or a current one left')
        n, eq, fields, only = compare_rows({'a': {'x': 1}, 'b': {'x': 1}, 'c': {'x': 1}, 't_a#0': {'x': 1}},
                                           {'a': {'x': 1, 'turns': None}, 'b': {'x': 2}, 'c': {'x': None},
                                            't_a#0': {'x': 9}}, {'t_a#0', 't_a#1', 't_a#5'})
        if (n, eq, dict(fields), only) != (3, 1, {'x': 2}, []):
            bad.append(f'compare_rows gave {(n, eq, dict(fields), only)}')
        if format_only_keys({'a': {'x': 1}}, {'a': {'x': 1, 'turns': None}}) != ['turns']:
            bad.append('format_only_keys missed a key added as None')
    print('SELF-TEST ' + ('OK' if not bad else 'FAILED: ' + '; '.join(bad)))
    return 1 if bad else 0


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__.split('\n')[0])
    g = ap.add_mutually_exclusive_group(required=True)
    g.add_argument('--list', action='store_true')
    g.add_argument('--verify-manifest', action='store_true')
    g.add_argument('--archive', action='store_true')
    g.add_argument('--archive-paraphrases', action='store_true')
    g.add_argument('--paraphrase-kit', action='store_true')
    g.add_argument('--paraphrase-score', metavar='DIR')
    g.add_argument('--rerun', action='store_true', help='BILLS: the chains of CHAINS in parallel through the harness')
    g.add_argument('--stages', action='store_true', help='BILLS: judge, then router, on one re-scored --variant')
    ap.add_argument('--e5-budget', type=float, default=4.0, help='--stages: dollars over the judge store\'s total')
    ap.add_argument('--router-budget', type=float, default=4.0, help='--stages: dollars over the router store\'s total')
    g.add_argument('--rescore', action='store_true')
    g.add_argument('--verify-rescore', action='store_true')
    g.add_argument('--check', action='store_true')
    ap.add_argument('--variant', action='append', help='a store; repeatable (--archive: default every store)')
    ap.add_argument('--expert', default='che-1', help='--paraphrase-kit: the chemical expert who checks the pairs')
    ap.add_argument('--dry-run', action='store_true')
    a = ap.parse_args()
    if a.check:
        return self_test()
    if a.list:
        return list_items()
    if a.verify_manifest:
        return verify_manifest()
    if a.archive:
        r = archive(a.variant or list(DEFAULT_VARIANTS), a.dry_run)
        return 1 if r['problems'] else 0
    if a.archive_paraphrases:
        archive_paraphrases(a.dry_run)
        return 0
    if a.paraphrase_kit:
        return paraphrase_kit(a.expert)
    if a.paraphrase_score:
        return paraphrase_score(Path(a.paraphrase_score))
    if a.rerun:
        return rerun()
    if a.stages:
        if not a.variant or len(a.variant) != 1:
            raise SystemExit('--stages takes exactly one --variant')
        return stages(a.variant[0], a.e5_budget, a.router_budget)
    if not a.variant or len(a.variant) != 1:
        raise SystemExit('--rescore and --verify-rescore take exactly one --variant')
    return rescore(a.variant[0], a.dry_run) if a.rescore else verify_rescore(a.variant[0])


if __name__ == '__main__':
    raise SystemExit(main())
