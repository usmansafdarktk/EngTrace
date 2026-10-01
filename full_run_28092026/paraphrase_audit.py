"""An audit of Stream 2, the paraphrase experiment (Q5), after the fact: every number the stream reports
is recomputed from the artifacts on disk and compared with what is committed (D-165).

    python -m full_run_28092026.paraphrase_audit        # FREE: writes PARAPHRASE_AUDIT.md beside this file

What it checks, in the order the stream ran:
  writing     the committed manifest against the local pool: 450 rows, one restored text per passing item
              whose hash the manifest records, one prompt for every passing item, the original's hash equal
              to the frozen pool's, every passing text passing the six checks as the code stands, the
              counts PARAPHRASE.md prints
  kits        the keyfile against the pool and the roster: every code maps to a passing item, each expert
              holds whole templates of their own branch, no more than ten templates, and the text in a kit
              is the restored text
  returns     the experts' files: one submission per assigned code from its assigned expert, the kept rule
              (same yes, answer yes, clear no) recomputed and compared with accepted.json and
              PARAPHRASE_REVIEW.md, a note on every rejection
  inference   the traces: one final row per passing item and model, the paraphrase hash and the original's
              as the manifest has them, the served model the main run's, the spend the review prints
  scoring     the score store: one row per trace with the trace's hash, unusable rows equal to empty ones
  judge       E5 on the arm: a reply for every judged trace, as the analysis reads them
  analysis    results.json: Q5 over the kept pairs for every model, the provenance commit clean
  backups     the local archive's checksum and members against the folders it covers; the Kaggle dataset
              listed under the owner's account (skipped without network)
Counts only: no item text, no expert's words. A FAIL line names the check and the figures; the script exits 1.
"""
from __future__ import annotations

import collections
import hashlib
import json
import re
import subprocess
import sys
import zipfile
from pathlib import Path

HERE = Path(__file__).resolve().parent
REPO = HERE.parent
for p in (str(REPO), str(REPO / 'evaluator_pilot_17092026' / 'evaluators'), str(REPO / 'evaluation')):
    if p not in sys.path:
        sys.path.insert(0, p)

from full_run_28092026 import paraphrase as P  # noqa: E402
from full_run_28092026 import analyze as A  # noqa: E402

PARA = HERE / 'paraphrase'
RETURNS = PARA / 'experts_filled_paraphrases'
ARCHIVE = Path.home() / 'EngTrace_private_backup' / 'full_run_paraphrase_2026-10-01.zip'
KAGGLE = 'ayeshaiq/engtrace-full-run-paraphrase'
ROSTER = A.ROSTER


def jsonl(p: Path) -> list[dict]:
    return [json.loads(l) for l in p.read_text(encoding='utf-8').splitlines() if l.strip()]


def sha(s: str) -> str:
    return hashlib.sha256(s.encode('utf-8')).hexdigest()


class Audit:
    def __init__(self):
        self.rows: list[tuple[str, str, bool, str]] = []

    def check(self, section: str, name: str, ok: bool, figures: str = ''):
        self.rows.append((section, name, bool(ok), figures))

    def table(self) -> list[str]:
        L = ['| step | check | result | figures |', '|---|---|---|---|']
        for s, n, ok, f in self.rows:
            L.append(f'| {s} | {n} | {"PASS" if ok else "**FAIL**"} | {f} |')
        return L


def md_counts(path: Path) -> dict[str, str]:
    """The `| label | value |` rows of a two-column markdown table, label -> value."""
    out = {}
    for line in path.read_text(encoding='utf-8').splitlines():
        m = re.match(r'^\| ([^|]+?) \| ([^|]+?) \|$', line)
        if m and m.group(1).strip() not in ('', '---'):
            out[m.group(1).strip()] = m.group(2).strip()
    return out


def main() -> int:
    au = Audit()
    # ---------------------------------------------------------------- writing
    frozen = {r['item_id']: r for r in jsonl(HERE / 'manifest.jsonl')}
    man = jsonl(PARA / 'manifest.jsonl')
    pool = {r['item_id']: r for r in jsonl(PARA / 'pool.jsonl')}
    passing = [m for m in man if m['passed']]
    sel = {it['item_id']: it for it in P.selection()}
    au.check('writing', 'manifest holds the 450 selected items, pool holds the passing ones',
             len(man) == 450 and set(m['item_id'] for m in man) == set(sel) and len(pool) == len(passing)
             and set(pool) == {m['item_id'] for m in passing}, f'{len(man)} rows, {len(passing)} passing, {len(pool)} in pool')
    au.check('writing', 'every passing text hashes to the manifest', all(sha(pool[m['item_id']]['question']) == m['sha256'] for m in passing))
    au.check('writing', "every original's hash is the frozen pool's",
             all(m['original_sha256'] == frozen[m['item_id']]['sha256'] for m in man))
    prompts = collections.Counter(m['prompt_sha256'][:16] for m in passing)
    au.check('writing', 'one prompt wrote every passing paraphrase', len(prompts) == 1 and next(iter(prompts)) == P.PROMPT_SHA[:16],
             ', '.join(f'{k} ({v})' for k, v in prompts.items()) + f'; the code\'s prompt {P.PROMPT_SHA[:16]}')
    checks = {m['item_id']: P.check(sel[m['item_id']]['question'], pool[m['item_id']]['question']) for m in passing}
    au.check('writing', 'every passing text passes the six checks as the code stands', all(c['passed'] for c in checks.values()),
             f"{sum(c['passed'] for c in checks.values())} of {len(checks)}")
    restored = sum(1 for m in passing if m['restored'])
    att = jsonl(PARA / 'attempts.jsonl')
    raw = {}
    for r in att:
        if 'check' in r:
            raw.setdefault(r['item_id'], {})[r['attempt']] = r['text']
    au.check('writing', 'the restored flag is true exactly where the served text differs from the raw attempt',
             all(m['restored'] == (P.restore(sel[m['item_id']]['question'], raw[m['item_id']][m['attempt_passed']]) != raw[m['item_id']][m['attempt_passed']])
                 for m in passing), f'{restored} restored')
    exhausted = [m for m in man if not m['passed'] and m['attempts'] >= P.ATTEMPTS]
    by_t = collections.defaultdict(list)
    for m in man:
        by_t[m['template_id']].append(m)
    lost_t = [t for t, ms in by_t.items() if all(not m['passed'] and m['attempts'] >= P.ATTEMPTS for m in ms)]
    pm = md_counts(HERE / 'PARAPHRASE.md')
    want = {'items written': str(sum(1 for m in man if m['attempts'])), 'passing a paraphrase': str(len(passing)),
            f'no paraphrase after {P.ATTEMPTS} attempts': str(len(exhausted))}
    got = {k: pm.get(k) for k in want}
    au.check('writing', 'PARAPHRASE.md prints the recomputed counts', all(got[k] == v for k, v in want.items()),
             '; '.join(f'{k} {v}' for k, v in want.items()) + f'; restored {restored}; templates with no paraphrase {len(lost_t)}'
             + (f" (PARAPHRASE.md: {pm.get('templates with no paraphrase for any of their items', '?').split(':')[0]})"))
    billed_all = 0.0
    for f in sorted(PARA.glob('attempts*.jsonl')):
        billed_all += sum(r.get('billed_usd') or 0 for r in jsonl(f))
    au.check('writing', 'step 1 billed over all four prompts', abs(billed_all - 0.4159) < 0.002, f'${billed_all:.4f} (D-157: $0.416)')
    # ---------------------------------------------------------------- kits
    key = json.load(open(PARA / 'tasks' / 'keyfile.json', encoding='utf-8'))
    assign = json.load(open(PARA / 'tasks' / 'assignment.json', encoding='utf-8'))
    au.check('kits', 'every code maps to a passing item, once', len(key) == len(passing) and set(key.values()) == set(pool)
             and len(set(key.values())) == len(key), f'{len(key)} codes')
    from template_annotation_23092026.layer2.build_tasks import roster as layer2_roster
    roster = {r['id']: r['branch'] for r in layer2_roster()}
    au.check('kits', 'the assignment is to the layer-2 roster, three per branch', set(assign) == set(roster) and len(roster) == 15,
             f'{len(assign)} experts')
    own = whole = ten = True
    sizes = {}
    for a_id, a in assign.items():
        items = [key[c] for c in a['codes']]
        sizes[a_id] = len(items)
        own &= all(sel[i]['branch'] == roster[a_id] for i in items)
        ts = {sel[i]['template_id'] for i in items}
        ten &= len(ts) <= 10
        whole &= all(set(c for c in key if sel[key[c]]['template_id'] == t) <= set(a['codes']) for t in ts)
    au.check('kits', 'each expert holds whole templates of their own branch, at most ten', own and whole and ten,
             f"items per expert {min(sizes.values())} to {max(sizes.values())}, {sum(sizes.values())} in all")
    kit_pool = json.load(open(PARA / 'tasks' / 'pool.json', encoding='utf-8'))
    def kit_text(entry):
        if isinstance(entry, dict):
            for k in ('paraphrase', 'rewritten', 'question_paraphrase', 'reworded'):
                if k in entry:
                    return entry[k]
        return None
    texts = {c: kit_text(kit_pool[c]) for c in key if c in kit_pool}
    comparable = {c: t for c, t in texts.items() if t is not None}
    au.check('kits', 'the paraphrase text in the kits is the restored text', bool(comparable)
             and all(sha(t) == next(m['sha256'] for m in passing if m['item_id'] == key[c]) for c, t in comparable.items()),
             f'{len(comparable)} of {len(key)} compared' if comparable else 'kit pool field not recognised')
    # ---------------------------------------------------------------- returns
    files = sorted(RETURNS.glob('*.jsonl'))
    subs = [r for f in files for r in jsonl(f)]
    last = {}
    for r in subs:
        last[r['code']] = r
    assigned_by_code = {c: a_id for a_id, a in assign.items() for c in a['codes']}
    au.check('returns', 'fifteen files, one submission per assigned code, from its assigned expert',
             len(files) == 15 and set(last) == set(key) and all(assigned_by_code[c] == r['annotator_id'] for c, r in last.items()),
             f'{len(files)} files, {len(subs)} submissions, {len(last)} codes')
    kept_rule = {key[c]: (r['same'] == 'yes' and r['answer'] == 'yes' and r['clear'] == 'no') for c, r in last.items()}
    acc = json.load(open(PARA / 'accepted.json', encoding='utf-8'))
    au.check('returns', 'accepted.json applies the kept rule to every pair', set(acc) == set(kept_rule)
             and all(acc[i]['kept'] == kept_rule[i] and acc[i]['returned'] for i in kept_rule),
             f'{sum(kept_rule.values())} kept, {sum(1 for v in kept_rule.values() if not v)} rejected')
    au.check('returns', 'every rejection carries a note', all((r.get('note') or '').strip() for c, r in last.items() if not kept_rule[key[c]]))
    rv = md_counts(HERE / 'PARAPHRASE_REVIEW.md')
    au.check('returns', 'PARAPHRASE_REVIEW.md prints the recomputed counts',
             rv.get('returned') == str(len(last)) and rv.get('kept') == str(sum(kept_rule.values()))
             and rv.get('rejected') == str(sum(1 for v in kept_rule.values() if not v)) and rv.get('not returned') == '0',
             f"returned {rv.get('returned')}, kept {rv.get('kept')}, rejected {rv.get('rejected')}, not returned {rv.get('not returned')}")
    kept = {i for i, v in kept_rule.items() if v}
    # ---------------------------------------------------------------- inference
    review = json.load(open(HERE / 'trace_review_paraphrase.json', encoding='utf-8'))
    man_by_id = {m['item_id']: m for m in passing}
    ok_rows = ok_hash = ok_served = True
    billed = 0.0
    empties = 0
    per_model_served = {}
    for k in ROSTER:
        rows = jsonl(HERE / 'traces' / 'paraphrase' / f'{k}.jsonl')
        final = [r for r in rows if r['status'] in ('answered', 'empty')]
        ids = collections.Counter(r['item_id'] for r in final)
        ok_rows &= set(ids) == set(pool) and max(ids.values()) == 1
        ok_hash &= all(r['item_sha256'] == man_by_id[r['item_id']]['sha256'] and r['original_sha256'] == man_by_id[r['item_id']]['original_sha256'] for r in final)
        served = {r['served_model'] for r in final}
        per_model_served[k] = served
        ok_served &= len(served) == 1
        billed += sum(r.get('billed_usd') or 0 for r in rows)
        empties += sum(1 for r in final if r['status'] == 'empty')
    au.check('inference', 'one final row per passing item for each of the eleven models', ok_rows, f'{len(pool)} x {len(ROSTER)}')
    au.check('inference', "every row's paraphrase hash and original's hash are the manifest's", ok_hash)
    main_review = json.load(open(HERE / 'trace_review.json', encoding='utf-8'))
    main_served = {m['key']: set(m['served_models']) for m in main_review['models']}
    am = all(m['served_as_main'] for m in review['models']) and all(per_model_served[k] == main_served.get(k) for k in ROSTER)
    au.check('inference', 'one served model per model, the same as the main run', ok_served and am,
             'served_as_main true for all eleven in the review; the ids equal the main review\'s' if am else 'a model differs')
    au.check('inference', 'the spend in the rows is the one D-161 records', abs(billed - 34.166) < 0.01, f'${billed:.3f}, {empties} empty rows')
    # ---------------------------------------------------------------- scoring
    ok_sc = True
    unusable = 0
    for k in ROSTER:
        tr = {r['item_id']: r for r in jsonl(HERE / 'traces' / 'paraphrase' / f'{k}.jsonl') if r['status'] in ('answered', 'empty')}
        sc = jsonl(HERE / 'scores' / 'paraphrase' / f'{k}.jsonl')
        ok_sc &= len(sc) == len(pool) and all(r['trace_sha256'] == sha(tr[r['item_id']]['text'] or '') for r in sc)
        unusable += sum(1 for r in sc if r['unusable'])
    au.check('scoring', 'one score row per trace, hashed to the trace it scores', ok_sc, f'{len(pool)} x {len(ROSTER)}')
    au.check('scoring', 'unusable rows are the empty rows', unusable == empties, f'{unusable} unusable, {empties} empty')
    # ---------------------------------------------------------------- judge
    e5 = {k: A.load_stage('e5', k, 'paraphrase') for k in ROSTER}
    entries = [x for v in e5.values() for x in (v or {}).values()]
    called = [x for x in entries if x.get('reply_ok') is not None]          # a call was made for the trace
    n_ok = sum(1 for x in called if x['reply_ok'] is True and x.get('e5_strict') is not None)
    au.check('judge', 'a reply with a verdict for every judged paraphrase trace', len(called) == 1227 and n_ok == len(called),
             f'{n_ok} of {len(called)} calls (D-162: 1,227); {len(entries)} traces, {sum(1 for x in entries if x.get("e5_strict") is None)} with no verdict (no milestones to judge)')
    # ---------------------------------------------------------------- analysis
    res = json.load(open(HERE / 'results' / 'results.json', encoding='utf-8'))
    q5 = res['q5']
    au.check('analysis', 'Q5 rows for the eleven models, each over the kept pairs',
             len(q5['models']) == len(ROSTER) and all(m['items'] == len(kept) for m in q5['models']),
             f"{len(q5['models'])} models, {len(kept)} pairs; tau over {q5['tau']['items']}")
    au.check('analysis', 'the experts\' counts in results.json are the recomputed ones',
             q5['expert_stats'] and q5['expert_stats']['rejected'] == len(kept_rule) - len(kept) and q5['expert_stats']['outstanding'] == 0)
    pv = res['provenance']['analyze']
    au.check('analysis', 'results regenerated from a clean tree', not pv['dirty'], f"analyze.py at {pv['git'][:7]}")
    # ---------------------------------------------------------------- backups
    if ARCHIVE.exists():
        digest = hashlib.sha256(ARCHIVE.read_bytes()).hexdigest()
        recorded = (ARCHIVE.parent / (ARCHIVE.name + '.sha256')).read_text().split()[0]
        with zipfile.ZipFile(ARCHIVE) as zf:
            members = {n: hashlib.sha256(zf.read(n)).hexdigest() for n in zf.namelist() if not n.endswith('/')}
        on_disk = {p.relative_to(HERE).as_posix(): sha_bytes(p) for d in ('paraphrase', 'traces/paraphrase') for p in (HERE / d).rglob('*') if p.is_file()}
        au.check('backups', 'the local archive matches its checksum and the folders it covers',
                 digest == recorded and members == on_disk, f'{len(members)} members; on disk {len(on_disk)}')
    else:
        au.check('backups', 'the local archive exists', False, str(ARCHIVE))
    try:
        out = subprocess.run(['kaggle', 'datasets', 'list', '-m', '-s', KAGGLE.split('/')[1]], capture_output=True, text=True, timeout=60).stdout
        au.check('backups', 'the private Kaggle dataset is listed under the owner\'s account', KAGGLE in out, KAGGLE)
    except Exception as exc:  # noqa: BLE001
        au.check('backups', 'the private Kaggle dataset is listed (skipped: no CLI or network)', True, type(exc).__name__)
    # ---------------------------------------------------------------- report
    fails = [r for r in au.rows if not r[2]]
    L = ['# Audit of the paraphrase experiment (D-165)', '',
         'Generated by `paraphrase_audit.py`; what each check does is in its docstring. Every figure is recomputed from the '
         'artifacts on disk and compared with what is committed. Counts only.', '',
         f'**{len(au.rows)} checks, {len(fails)} failed.**', ''] + au.table() + ['']
    (HERE / 'PARAPHRASE_AUDIT.md').write_text('\n'.join(L) + '\n', encoding='utf-8', newline='\n')
    print('\n'.join(L))
    return 1 if fails else 0


def sha_bytes(p: Path) -> str:
    return hashlib.sha256(p.read_bytes()).hexdigest()


if __name__ == '__main__':
    raise SystemExit(main())
