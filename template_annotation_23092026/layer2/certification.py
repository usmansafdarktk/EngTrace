"""Layer 2 - the certification status of every template, across all rounds.

    python -m template_annotation_23092026.layer2.certification \\
        --labels <round-1 folder> <round-2 folder> <round-3 folder>

The label folders are given in round order. For each template this finds the latest round
that reviewed it, takes the three experts' verdicts from that round, and regenerates the
instances the experts were shown in that round from the CURRENT code, comparing question and
solution text byte for byte. A template is certified when its latest round is a unanimous
approval AND the current code still produces exactly what was reviewed; a template whose
code or constants moved after its review is not, however it was judged.

Reads, for round 1, tasks/{keyfile.jsonl, pool.json} and for round N, tasks_round<N>/ (all
git-ignored, built by build_tasks.py); plants are skipped. Writes CERTIFICATION.md beside
this file.
"""
from __future__ import annotations

import argparse
import collections
import datetime as dt
import json
import subprocess
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
REPO = HERE.parent.parent
if str(REPO) not in sys.path:
    sys.path.insert(0, str(REPO))

from tests.template_integrity.core import discover, generate  # noqa: E402

from template_annotation_23092026.layer2.score import read_rows, tasks_dir  # noqa: E402


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument('--labels', nargs='+', required=True, help='label folders, round 1 first')
    a = ap.parse_args()

    refs = {r.template_id: r for r in discover()}
    latest: dict[str, dict] = {}                  # template_id -> the latest round that reviewed it
    rounds = []
    for n, folder in enumerate(a.labels, start=1):
        tdir = tasks_dir(n)
        key = {}
        for ln in (tdir / 'keyfile.jsonl').open(encoding='utf8'):
            if ln.strip():
                k = json.loads(ln)
                key[k['code']] = k
        pool = json.loads((tdir / 'pool.json').read_text(encoding='utf8'))
        build = json.loads((tdir / 'BUILD.json').read_text(encoding='utf8'))
        verdicts = collections.defaultdict(dict)
        for r in read_rows(Path(folder)):
            k = key.get(r['code'])
            if k is None:
                raise SystemExit(f'round {n}: code {r["code"]} is not in {tdir.name}/keyfile.jsonl')
            if k['kind'] == 'template':
                verdicts[k['template_id']][r['annotator_id']] = r['decision']
        n_templates = 0
        for code, k in key.items():
            if k['kind'] != 'template':
                continue
            n_templates += 1
            latest[k['template_id']] = {'round': n, 'build': build['git_head'][:7],
                                        'verdicts': dict(sorted(verdicts[k['template_id']].items())),
                                        'instances': pool[code]['instances'], 'branch': k['branch']}
        rounds.append((n, Path(folder).name, build['git_head'][:7], n_templates,
                       sum(len(v) for v in verdicts.values())))

    rows, not_cert = [], []
    for tid in sorted(refs):
        info = latest.get(tid)
        if info is None:
            not_cert.append((tid, 'never reviewed'))
            continue
        same = True
        for inst in info['instances']:
            g = generate(refs[tid], inst['seed'], capture=False)
            if not g.ok or g.question != inst['question'] or g.solution != inst['solution']:
                same = False
                break
        v = list(info['verdicts'].values())
        unanimous = len(v) == 3 and all(d == 'Approve' for d in v)
        ok = unanimous and same
        reason = '' if ok else ('; '.join(x for x in (
            '' if len(v) == 3 else f'{len(v)} verdicts',
            '' if all(d == 'Approve' for d in v) else f"{sum(d == 'Reject' for d in v)} rejection(s)",
            '' if same else 'current output differs from the reviewed instances') if x))
        rows.append((tid, info, v, same, ok))
        if not ok:
            not_cert.append((tid, reason))

    head = subprocess.run(['git', 'rev-parse', 'HEAD'], capture_output=True, text=True, cwd=REPO).stdout.strip()
    by_round = collections.Counter(info['round'] for _, info, _, _, ok in rows if ok)
    out = ['# Layer 2 - certification status of every template\n',
           f"Generated {dt.datetime.now(dt.timezone.utc).isoformat(timespec='seconds')} by `certification.py` at git "
           f"`{head[:10]}` over {len(refs)} templates. A template is certified when the latest round that "
           f"reviewed it is a unanimous approval by its three own-branch experts and the current code still "
           f"produces, byte for byte, the five instances they were shown.\n",
           '| Round | Labels | Built at | Templates reviewed | Verdicts |\n|---:|---|---|---:|---:|']
    for n, name, bh, nt, nv in rounds:
        out.append(f'| {n} | `{name}` | `{bh}` | {nt} | {nv} |')
    out += ['', '| | Templates |\n|---|---:|',
            f'| certified | {sum(1 for r in rows if r[4])} of {len(refs)} |']
    for n in sorted(by_round):
        out.append(f'| of which last reviewed in round {n} | {by_round[n]} |')
    out.append(f'| not certified | {len(not_cert)} |\n')
    if not_cert:
        out.append('| Template | Why not |\n|---|---|')
        out += [f'| `{t}` | {why} |' for t, why in not_cert]
        out.append('')
    out.append('| Template | Branch | Last reviewed | Verdicts | Current output = reviewed | Certified |\n'
               '|---|---|---:|---|---|---|')
    for tid, info, v, same, ok in rows:
        out.append(f"| `{tid}` | {info['branch'].split('_')[0]} | round {info['round']} | "
                   f"{' '.join(d[0] for d in v)} | {'yes' if same else 'NO'} | {'yes' if ok else 'no'} |")
    (HERE / 'CERTIFICATION.md').write_text('\n'.join(out) + '\n', encoding='utf8')
    print(f"certified {sum(1 for r in rows if r[4])} of {len(refs)}; by last round {dict(sorted(by_round.items()))}; "
          f"not certified {len(not_cert)}")
    for t, why in not_cert:
        print('  ', t, '-', why)


if __name__ == '__main__':
    main()
