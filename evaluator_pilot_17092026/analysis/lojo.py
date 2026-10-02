"""Leave-one-judge-out on the published Tribunal, with controls, and the per-judge bias against the experts (X2):
the ablation the July rebuttal promised (docs/EVALUATION_NEXT_STEPS.md A3 and A4; D-174).

    PY=evaluator_pilot_17092026/.venv/Scripts/python
    $PY evaluator_pilot_17092026/analysis/lojo.py            # FREE, offline: RESULTS_LOJO.md and analysis/out/lojo.json

WHAT IT REPLAYS. E0-3J is the published framework with its three judges connected (GPT-5, Claude Opus 4.5, Gemini 3.1
Pro) on the pilot's 300 traces (RESULTS_E0, RESULTS_E1). For every judged trace the framework's own arithmetic is
redone here from what the run left on disk, and nothing is called: the Tier 1 validity matrix V from the scorer cache
(`scores/_cache`, keyed by the step lists and the library versions, as `evaluators/e0_tribunal.py` keys it); the
mismatched predicted steps from the Hungarian assignment on V; each judge's category per step from its stored raw
reply, parsed with the framework's own logic and mapped to its scalar (Alternative Correct 1.0, Calculation Error
0.5, anything else 0.0); the panel's consensus (a strict majority, else the conservative minimum); the gold step each
recovered score is written to (the cross-encoder's single-pair similarity, served from the same cache); the Hungarian
assignment again; precision, recall and recovered F1. Under the full panel the replay must reproduce the stored
recovered_f1 of every judged trace, and the report says on how many it does. Then the same arithmetic runs with one
judge's votes removed: for each trace model the drop of its own family's judge (GPT-5 traces without GPT-5, Claude
Opus 4.7 traces without Claude Opus 4.5, Gemini 3.1 Pro traces without Gemini 3.1 Pro), the two other drops as
placebos, and DeepSeek R1 and Llama 3.1 70B, which share no family with any judge, as controls under every drop.

WHAT IT REPORTS. Per trace model and drop: the mean recovered F1 over its 60 traces (an unjudged trace keeps its
stored score, since no judge touched it), the mean change from the full panel with a template-resampled interval (15
templates, B = 2,000) and a trace-resampled one, how many traces move, and the smallest mean change 15 templates
detect at 80% power. Then the panel's step verdicts against the experts' step labels (version 2, with the
adjudication): on the judged steps that carry a label, whether the panel's consensus calls the step an error
(scalar below 1.0) against whether the experts call it incorrect, as accuracy, precision and recall per model under
the full panel and each drop. X2: per judge (the three of E0-3J and the three of E1, from the E1 store) and per trace
model, the mean of the judge's scalar minus the experts' (1.0 for correct or alternative correct, 0.0 for incorrect;
not-a-claim steps left out) with a trace-resampled interval, and the judge's lenience on expert-incorrect steps
(the share it calls Alternative Correct) and harshness on expert-correct ones (the share it calls an error); the
family cells marked, and the family effect (own-family bias minus other-traces bias) with its interval.

WHAT IT DOES NOT SHOW. This is the published framework as it was run on 15 templates and five models; the full run
does not use it (D-105). A null here says the panel's score on a family's traces does not depend on that family's
judge at the size this slice resolves, which the detectable-difference column states; it does not validate the
judges' accuracy, which RESULTS_X1 measures.
"""
from __future__ import annotations

import collections
import glob
import hashlib
import json
import os
import statistics
import sys
import zlib

import numpy as np
from scipy.optimize import linear_sum_assignment

_HERE = os.path.dirname(os.path.abspath(__file__))
_PILOT = os.path.dirname(_HERE)
_REPO = os.path.dirname(_PILOT)
for p in (os.path.join(_REPO, 'evaluation'), _HERE, os.path.join(_PILOT, 'annotation')):
    if p not in sys.path:
        sys.path.insert(0, p)

from engineering_parser import extract_steps  # noqa: E402
from judge_probe import parse  # noqa: E402
import score_against_labels as S  # noqa: E402

OUT_MD = os.path.join(_PILOT, 'RESULTS_LOJO.md')
OUT_JSON = os.path.join(_HERE, 'out', 'lojo.json')
CACHE_DIR = os.path.join(_PILOT, 'scores', '_cache')
LABELS = os.path.join(_PILOT, 'experts_filled_labels', 'version_2', 'labels')
PROVIDERS = ('openai', 'anthropic', 'google')
FAMILY = {'gpt-5': 'openai', 'claude-opus-4.7': 'anthropic', 'gemini-3.1-pro': 'google',
          'deepseek-r1': None, 'llama-3.1-70b': None}
MODELS = ['gpt-5', 'claude-opus-4.7', 'gemini-3.1-pro', 'deepseek-r1', 'llama-3.1-70b']
B = 2000
Z80, Z975 = 0.8416212335729143, 1.959963984540054


# ------------------------------------------------------------------ inputs

def libraries() -> dict:
    """The scorer library versions the cache keys carry (evaluators/e0_tribunal.py, verbatim list)."""
    from importlib.metadata import version, PackageNotFoundError
    out = {}
    for pkg in ('transformers', 'sentence-transformers', 'bert-score', 'rouge-score', 'tokenizers', 'scipy', 'numpy'):
        try:
            out[pkg] = version(pkg)
        except PackageNotFoundError:
            out[pkg] = None
    return out


LIBS = json.dumps(libraries(), sort_keys=True)


def key(kind, *args) -> str:
    return hashlib.sha256(json.dumps([kind, LIBS, args], sort_keys=True).encode()).hexdigest()


def rows(path):
    return [json.loads(l) for l in open(path, encoding='utf-8') if l.strip()]


def current(d):
    """The store's rows under its latest config (e1_analysis.current)."""
    rs = [r for f in glob.glob(os.path.join(_PILOT, 'scores', d, '*.jsonl')) for r in rows(f) if not r.get('error')]
    latest = max(rs, key=lambda r: r['ts'])['config_sha256']
    return {(r['model_key'], r['item_id']): r for r in rs if r['config_sha256'] == latest}


def cache() -> dict:
    store = {}
    for fn in sorted(os.listdir(CACHE_DIR)):
        if not fn.endswith('.jsonl') or '_judge_replies' in fn:
            continue
        for ln in open(os.path.join(CACHE_DIR, fn), encoding='utf-8'):
            try:
                rec = json.loads(ln)
                store[rec['k']] = rec['v']
            except (json.JSONDecodeError, KeyError):
                continue
    return store


def texts() -> dict:
    out = {}
    for f in os.listdir(os.path.join(_PILOT, 'traces')):
        if f.endswith('.jsonl'):
            for r in rows(os.path.join(_PILOT, 'traces', f)):
                if r.get('ok'):
                    out[(r['model_key'], r['item_id'])] = r['text']
    return out


# ------------------------------------------------------------------ the framework's arithmetic, offline

def scalar(cat) -> float:
    c = str(cat).lower()
    return 1.0 if 'alternative' in c else 0.5 if 'calculation' in c else 0.0


def consensus(votes: list[float]):
    if not votes:
        return None
    counts = collections.Counter(votes)
    for sc, n in counts.most_common():
        if n > len(votes) / 2:
            return sc
    return min(votes)


def f1_of(V: np.ndarray):
    M, N = V.shape
    if M == 0 or N == 0:
        return 0.0, 0.0, 0.0
    r, c = linear_sum_assignment(V, maximize=True)
    R = float(V[r, c].sum())
    p, rc = R / N, R / M
    return (2 * p * rc / (p + rc) if p + rc > 0 else 0.0), p, rc


class Trace:
    """One judged trace: everything the framework's Tier 2 arithmetic needs, read once."""

    def __init__(self, row, item, text, store):
        self.row = row
        gt_txt, gt_val, _ = extract_steps(item['solution'])
        pr_txt, pr_val, _ = extract_steps(text)
        self.gt, self.pr = gt_txt, pr_txt
        self.nonempty = [i for i, s in enumerate(pr_txt) if s.strip()]
        v = store.get(key('tier1', gt_txt, gt_val, pr_txt, pr_val))
        self.ok = v is not None
        if not self.ok:
            return
        self.V = np.array(v['V'], dtype=float).reshape(v['shape'])
        M, N = self.V.shape
        r, c = linear_sum_assignment(self.V, maximize=True)
        matched = set(c[self.V[r, c] == 1.0])
        self.mismatch = [j for j in range(N) if j not in matched]
        self.votes = {j: {} for j in self.mismatch}          # step -> {provider: scalar}
        self.model_of = {}
        for call in row['calls']:
            self.model_of[call['provider']] = call.get('model')
            res = parse(call.get('text') or '') if call.get('ok') else None
            for x in res or []:
                if isinstance(x, dict) and x.get('step_index') in self.votes:
                    self.votes[x['step_index']][call['provider']] = scalar(x.get('category'))
        self.sims = {}
        self.missing_pairs = 0
        for j in self.mismatch:
            if not self.votes[j]:
                continue
            sims = []
            for i in range(M):
                s = store.get(key('pair', gt_txt[i], pr_txt[j]))
                if s is None:
                    self.missing_pairs += 1
                sims.append(s)
            self.sims[j] = sims

    def f1(self, keep: set[str]):
        """Recovered F1 with only the judges in `keep` voting; None when a needed pair similarity is absent."""
        V = self.V.copy()
        for j in self.mismatch:
            sc = consensus([v for p, v in self.votes[j].items() if p in keep])
            if sc is None or sc <= 0:
                continue
            sims = self.sims[j]
            if any(s is None for s in sims):
                return None
            best = int(np.argmax(sims))
            V[best, j] = max(V[best, j], sc)
        return f1_of(V)[0]

    def verdicts(self, keep: set[str]) -> dict:
        """{raw step index: consensus scalar} under the panel `keep`, for the steps with a vote."""
        out = {}
        for j in self.mismatch:
            sc = consensus([v for p, v in self.votes[j].items() if p in keep])
            if sc is not None:
                out[j] = sc
        return out


# ------------------------------------------------------------------ statistics

def seed_of(*parts: str) -> int:
    """A seed that is a function of its labels alone: Python's hash() is salted per process and would make the
    intervals differ between runs."""
    return zlib.crc32('|'.join(parts).encode('utf-8')) % (2 ** 31)


def boot_mean(values_by_group: dict, seed: int, n=B):
    """Mean over all values with a cluster bootstrap over the groups (templates); and the trace bootstrap."""
    groups = [np.array(v, dtype=float) for v in values_by_group.values() if len(v)]
    allv = np.concatenate(groups) if groups else np.array([])
    if not len(allv):
        return float('nan'), [float('nan')] * 2, [float('nan')] * 2
    rng = np.random.default_rng(seed)
    idx = rng.integers(0, len(groups), size=(n, len(groups)))
    sums = np.array([g.sum() for g in groups])
    cnts = np.array([len(g) for g in groups])
    cl = sums[idx].sum(axis=1) / cnts[idx].sum(axis=1)
    tr = allv[rng.integers(0, len(allv), size=(n, len(allv)))].mean(axis=1)
    return float(allv.mean()), [float(np.percentile(cl, 2.5)), float(np.percentile(cl, 97.5))], \
        [float(np.percentile(tr, 2.5)), float(np.percentile(tr, 97.5))]


def detectable(values_by_group: dict) -> float:
    means = [float(np.mean(v)) for v in values_by_group.values() if len(v)]
    if len(means) < 2:
        return float('nan')
    return (Z975 + Z80) * float(np.std(means, ddof=1)) / np.sqrt(len(means))


# ------------------------------------------------------------------ the analysis

def expert_steps():
    truth, _ = S.build_truth(S.read_labels(LABELS), S.read_consensus(LABELS))
    keyfile = {r['code']: (r['model_key'], r['item_id']) for r in rows(os.path.join(S.TASKS, 'keyfile.jsonl'))}
    return {keyfile[c]: t for c, t in truth.items() if c in keyfile}


def main() -> int:
    store = cache()
    items = {r['item_id']: r for r in rows(os.path.join(_PILOT, 'slice', 'manifest.jsonl'))}
    tt = texts()
    truth = expert_steps()
    e03 = current('e0_3j')
    e1 = current('e1')
    traces = {}
    missing_v = missing_pairs = 0
    for k, row in e03.items():
        if not row.get('calls'):
            continue
        t = Trace(row, items[k[1]], tt[k], store)
        if not t.ok:
            missing_v += 1
            continue
        missing_pairs += t.missing_pairs
        traces[k] = t
    judged = len(traces)
    # replication under the full panel
    full = set(PROVIDERS)
    repro = {k: (t.f1(full), float(t.row['scores']['recovered_f1'])) for k, t in traces.items()}
    repl_ok = {k for k, (a, b) in repro.items() if a is not None and abs(a - b) < 1e-9}
    max_dev = max((abs(a - b) for a, b in repro.values() if a is not None), default=0.0)
    unreplayed = [k for k, (a, _) in repro.items() if a is None]
    res = {'judged': judged, 'tier1_missing': missing_v, 'pair_similarities_missing': missing_pairs,
           'unreplayed': len(unreplayed), 'reproduced': len(repl_ok), 'max_abs_deviation': max_dev,
           'libraries': json.loads(LIBS), 'models': {}, 'x2': {}}
    drops = {f'drop_{p}': full - {p} for p in PROVIDERS}
    for m in MODELS:
        ks = sorted(k for k in e03 if k[0] == m)
        base = {k: float(e03[k]['scores']['recovered_f1']) for k in ks}
        f_full = {k: (traces[k].f1(full) if k in traces and k in repl_ok else base[k]) for k in ks}
        entry = {'traces': len(ks), 'judged': sum(k in traces for k in ks), 'replayed': sum(k in repl_ok for k in ks),
                 'family_judge': FAMILY[m], 'f1_full': float(np.mean(list(f_full.values()))),
                 'f1_stored': float(np.mean(list(base.values()))), 'drops': {}}
        for name, keep in drops.items():
            f_var = {k: (traces[k].f1(keep) if k in traces and k in repl_ok else base[k]) for k in ks}
            change_by_t = collections.defaultdict(list)
            for k in ks:
                change_by_t[e03[k]['template_id']].append((f_var[k] if f_var[k] is not None else f_full[k]) - f_full[k])
            mean, ci_t, ci_tr = boot_mean(change_by_t, seed_of('lojo', m, name))
            moved = [k for k in ks if f_var[k] is not None and abs(f_var[k] - f_full[k]) > 1e-9]
            dropped = name.removeprefix('drop_')
            entry['drops'][name] = {
                'role': 'family' if FAMILY[m] == dropped else ('control' if FAMILY[m] is None else 'placebo'),
                'f1': float(np.mean([f_var[k] if f_var[k] is not None else f_full[k] for k in ks])),
                'change': mean, 'ci_templates': ci_t, 'ci_traces': ci_tr, 'moved': len(moved),
                'largest_move': max((abs(f_var[k] - f_full[k]) for k in moved), default=0.0),
                'detectable_templates': detectable(change_by_t)}
        # the panel's step verdicts against the experts, per variant
        entry['steps'] = {}
        aligned = 0
        for name, keep in [('full', full)] + list(drops.items()):
            c = collections.Counter()
            for k in ks:
                t = traces.get(k)
                tr = truth.get(k)
                if t is None or tr is None or k not in repl_ok or len(t.nonempty) != len(tr['steps']):
                    continue
                if name == 'full':
                    aligned += 1
                back = {raw: j for j, raw in enumerate(t.nonempty)}
                for raw, sc in t.verdicts(keep).items():
                    j = back.get(raw)
                    if j is None:
                        continue
                    lab = tr['steps'][j]['label']
                    if lab not in ('correct', 'alternative_correct', 'incorrect'):
                        continue
                    says_err, is_err = sc < 1.0, lab == 'incorrect'
                    c['tp' if says_err and is_err else 'fp' if says_err else 'fn' if is_err else 'tn'] += 1
            n = sum(c.values())
            entry['steps'][name] = {'steps': n, 'accuracy': (c['tp'] + c['tn']) / n if n else None,
                                    'precision': c['tp'] / (c['tp'] + c['fp']) if c['tp'] + c['fp'] else None,
                                    'recall': c['tp'] / (c['tp'] + c['fn']) if c['tp'] + c['fn'] else None,
                                    **dict(c)}
        entry['traces_aligned_with_labels'] = aligned
        res['models'][m] = entry
    # X2: per judge and trace model, the judge's scalar minus the experts', on judged labelled steps
    for store_name, st in (('E0-3J', e03), ('E1', e1)):
        per = collections.defaultdict(lambda: collections.defaultdict(list))      # judge -> model -> [(trace key, diff, lenient, harsh)]
        names = {}
        for k, row in st.items():
            if not row.get('calls'):
                continue
            tr = truth.get(k)
            if tr is None:
                continue
            t = traces.get(k) if store_name == 'E0-3J' else None
            if t is None:                      # E1 (or an E0-3J row the cache cannot replay): parse the votes here
                item, text = items[k[1]], tt[k]
                pr_txt = extract_steps(text)[0]
                nonempty = [i for i, s in enumerate(pr_txt) if s.strip()]
                votes = collections.defaultdict(dict)
                for call in row['calls']:
                    names[call['provider']] = call.get('model')
                    res_ = parse(call.get('text') or '') if call.get('ok') else None
                    for x in res_ or []:
                        if isinstance(x, dict) and isinstance(x.get('step_index'), int):
                            votes[x['step_index']][call['provider']] = scalar(x.get('category'))
            else:
                nonempty, votes = t.nonempty, t.votes
                names.update(t.model_of)
            if len(nonempty) != len(tr['steps']):
                continue
            back = {raw: j for j, raw in enumerate(nonempty)}
            for raw, by_judge in votes.items():
                j = back.get(raw)
                if j is None:
                    continue
                lab = tr['steps'][j]['label']
                if lab not in ('correct', 'alternative_correct', 'incorrect'):
                    continue
                ex = 0.0 if lab == 'incorrect' else 1.0
                for prov, sc in by_judge.items():
                    per[prov][k[0]].append((k, sc - ex, (ex == 0.0 and sc == 1.0), (ex == 1.0 and sc < 1.0), ex))
        out = {}
        for prov, by_model in per.items():
            jm = names.get(prov, prov)
            out[jm] = {'slot': prov, 'models': {}}
            for m in MODELS:
                vals = by_model.get(m, [])
                by_trace = collections.defaultdict(list)
                for k, d, _, _, _ in vals:
                    by_trace[k].append(d)
                mean, _, ci_tr = boot_mean(by_trace, seed_of('x2', store_name, prov, m))
                inc = [x for x in vals if x[4] == 0.0]
                cor = [x for x in vals if x[4] == 1.0]
                out[jm]['models'][m] = {
                    'steps': len(vals), 'traces': len(by_trace), 'bias': mean, 'ci_traces': ci_tr,
                    'family': (store_name == 'E0-3J' and FAMILY[m] == prov),
                    'lenient_on_incorrect': (sum(x[2] for x in inc) / len(inc)) if inc else None, 'incorrect_steps': len(inc),
                    'harsh_on_correct': (sum(x[3] for x in cor) / len(cor)) if cor else None, 'correct_steps': len(cor)}
            if store_name == 'E0-3J':
                fam = [m for m in MODELS if FAMILY[m] == prov]
                own = [x for m in fam for x in by_model.get(m, [])]
                other = [x for m in MODELS if m not in fam for x in by_model.get(m, [])]
                if own and other:
                    rng = np.random.default_rng(3)
                    def bt(vals):
                        bt_ = collections.defaultdict(list)
                        for k, d, *_ in vals:
                            bt_[k].append(d)
                        return [np.array(v) for v in bt_.values()]
                    go, gt_ = bt(own), bt(other)
                    def draw(gs):
                        i = rng.integers(0, len(gs), size=len(gs))
                        return float(np.concatenate([gs[x] for x in i]).mean())
                    diffs = [draw(go) - draw(gt_) for _ in range(B)]
                    out[jm]['family_effect'] = {'own_minus_other': float(np.mean([d for _, d, *_ in own]) - np.mean([d for _, d, *_ in other])),
                                                'ci_traces': [float(np.percentile(diffs, 2.5)), float(np.percentile(diffs, 97.5))],
                                                'own_steps': len(own), 'other_steps': len(other)}
        res['x2'][store_name] = out
    os.makedirs(os.path.dirname(OUT_JSON), exist_ok=True)
    with open(OUT_JSON, 'w', encoding='utf-8') as fh:
        json.dump(res, fh, indent=1)
    text = render(res)
    with open(OUT_MD, 'w', encoding='utf-8', newline='\n') as fh:
        fh.write(text)
    print(text)
    return 0


def f3(v):
    return '-' if v is None or (isinstance(v, float) and np.isnan(v)) else f'{v:.3f}'


def ci(c):
    return '-' if c is None or any(x is None or (isinstance(x, float) and np.isnan(x)) for x in c) else f'{c[0]:+.3f} to {c[1]:+.3f}'


def render(res) -> str:
    L = ['# Leave-one-judge-out on the published Tribunal, with controls; and the per-judge bias (X2)', '',
         'Generated by `analysis/lojo.py`; what it replays and reports is in its docstring. E0-3J, the published framework with '
         'its three judges connected, on the pilot\'s 300 traces; every number is recomputed offline from the stored Tier 1 '
         'matrices, the stored judge replies and the cached similarities, with nothing called (D-174; next steps A3 and A4).', '',
         '## The replay reproduces the run', '', '| | |', '|---|---|',
         f"| judged traces in E0-3J | {res['judged'] + res['tier1_missing']} |",
         f"| Tier 1 matrix found in the cache | {res['judged']} |",
         f"| traces whose recovered F1 the replay reproduces under the full panel (to 1e-9) | {res['reproduced']} of {res['judged']} |",
         f"| largest deviation | {res['max_abs_deviation']:.2e} |",
         f"| traces left out (a similarity the cache does not hold) | {res['unreplayed']} |", '',
         '## Recovered F1 with one judge dropped', '',
         'Per trace model, the mean recovered F1 over its 60 traces under the full panel and with each judge\'s votes removed, '
         'the change with a template-resampled 95% interval (15 templates) and a trace-resampled one, the traces that move, and '
         'the smallest mean change 15 templates detect at 80% power. "family" marks the drop of the judge from the trace model\'s '
         'own family, "placebo" a non-family judge, "control" a model with no judge in its family.', '',
         '| trace model | judged / 60 | F1, full panel | drop | role | F1 | change | 95% CI, templates | 95% CI, traces | traces moved | largest move | detectable |',
         '|---|---:|---:|---|---|---:|---:|---:|---:|---:|---:|---:|']
    for m, e in res['models'].items():
        for name, d in e['drops'].items():
            L.append(f"| `{m}` | {e['judged']} | {e['f1_full']:.3f} | {name.removeprefix('drop_')} | {d['role']} | {d['f1']:.3f} | "
                     f"{d['change']:+.4f} | {ci(d['ci_templates'])} | {ci(d['ci_traces'])} | {d['moved']} | {d['largest_move']:.3f} | "
                     f"{f3(d['detectable_templates'])} |")
    L += ['', "## The panel's step verdicts against the experts, under each panel", '',
          'On the judged steps that carry an expert label (traces whose step split matches the labels\'), whether the panel\'s '
          'consensus calls the step an error (scalar below 1.0) against whether the experts call it incorrect.', '',
          '| trace model | traces aligned | panel | steps | accuracy | precision on "error" | recall on "error" |',
          '|---|---:|---|---:|---:|---:|---:|']
    for m, e in res['models'].items():
        for name, s in e['steps'].items():
            L.append(f"| `{m}` | {e['traces_aligned_with_labels']} | {name.replace('drop_', 'without ')} | {s['steps']} | "
                     f"{f3(s['accuracy'])} | {f3(s['precision'])} | {f3(s['recall'])} |")
    for store_name, judges in res['x2'].items():
        L += ['', f'## X2, {store_name}: each judge\'s step verdict minus the experts\' label, per trace model', '',
              'Bias is the mean of the judge\'s scalar (1.0 Alternative Correct, 0.5 Calculation Error, 0.0 otherwise) minus the '
              'experts\' (1.0 correct or alternative correct, 0.0 incorrect) over the judged labelled steps, with a trace-resampled '
              '95% interval; lenient: the share of expert-incorrect steps the judge calls Alternative Correct; harsh: the share of '
              'expert-correct steps it calls an error. A family cell is the judge on its own family\'s traces.', '',
              '| judge | trace model | steps | bias | 95% CI | lenient on incorrect (n) | harsh on correct (n) | family |',
              '|---|---|---:|---:|---:|---:|---:|---|']
        for jm, j in judges.items():
            for m, c in j['models'].items():
                L.append(f"| `{jm}` | `{m}` | {c['steps']} | {f3(c['bias'])} | {ci(c['ci_traces'])} | "
                         f"{f3(c['lenient_on_incorrect'])} ({c['incorrect_steps']}) | {f3(c['harsh_on_correct'])} ({c['correct_steps']}) | "
                         f"{'yes' if c['family'] else ''} |")
        fe = [(jm, j['family_effect']) for jm, j in judges.items() if j.get('family_effect')]
        if fe:
            L += ['', '| judge | bias on own family\'s traces minus on the others | 95% CI, traces | own / other steps |', '|---|---:|---:|---:|']
            for jm, f in fe:
                L.append(f"| `{jm}` | {f['own_minus_other']:+.3f} | {ci(f['ci_traces'])} | {f['own_steps']} / {f['other_steps']} |")
    L.append('')
    return '\n'.join(L)


if __name__ == '__main__':
    sys.exit(main())
