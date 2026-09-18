"""Do candidate judges DISCRIMINATE? A probe for choosing E1/E5's judges.

E0's judges answered "Alternative Correct" to 93% of the steps sent to them
(RESULTS_E0). A judge that rubber-stamps is useless whatever family it comes from,
so the property that matters most is tested directly, on steps whose status is
already known without asking any model:

  SLIP   a step whose shown arithmetic is WRONG - found by E4's sympy checker and
         confirmed by reading (RESULTS_E4 audit). A discriminating judge must NOT
         call it "Alternative Correct".
  CLEAN  a step from a trace with the right final answer whose every checkable
         claim is arithmetically consistent. A judge should not flag it.

Each prompt is the framework's own Tribunal prompt, built by its own code
(`_tier2_tribunal_batch`), with exactly one step under review, so each verdict
maps to one known label. Replies are parsed with a copy of the framework's own
parsing logic.

    PY=evaluator_pilot_17092026/.venv/Scripts/python
    $PY evaluator_pilot_17092026/analysis/judge_probe.py build    # writes the probe set
    $PY evaluator_pilot_17092026/analysis/judge_probe.py run      # calls every judge (paid)
    $PY evaluator_pilot_17092026/analysis/judge_probe.py report
"""
import os as _os
_ANALYSIS = _os.path.dirname(_os.path.abspath(__file__))
_PILOT = _os.path.dirname(_ANALYSIS)
_REPO = _os.path.dirname(_PILOT)

import concurrent.futures as cf
import glob
import io
import json
import re
import sys
import time
from contextlib import redirect_stdout

OUT = _os.path.join(_PILOT, 'scores', 'judge_probe')
PROBE = _os.path.join(OUT, 'probe_set.json')

# id, route model, family, role in the probe
JUDGES = [
    ('gpt-5',          'openai/gpt-5',                     'OpenAI (E0 judge)',  'reference'),
    ('opus-4.5',       'anthropic/claude-opus-4.5',        'Anthropic (E0 judge)', 'reference'),
    ('glm-5.3',        'z-ai/glm-5.3',                     'Z.ai',       'candidate'),
    ('kimi-k3',        'moonshotai/kimi-k3',               'Moonshot',   'candidate'),
    ('grok-4.6',       'x-ai/grok-4.6',                    'xAI',        'candidate'),
    ('minimax-m3',     'minimax/minimax-m3',               'MiniMax',    'candidate'),
    ('nemotron-ultra', 'nvidia/nemotron-3-ultra-550b-a55b', 'NVIDIA',    'candidate'),
    ('mimo-v2.5-pro',  'xiaomi/mimo-v2.5-pro',             'Xiaomi',     'candidate'),
    ('seed-2.1',       'bytedance-seed/seed-2-1-turbo',    'ByteDance',  'candidate'),
]

# Slips confirmed REAL by reading (RESULTS_E4 audit): (model, item, text of the
# slipped EXPRESSION). The needle must be the erroneous expression itself, not a
# number in it: a first version matched '-64.16', '0.046' and '6.846' where those
# numbers are first computed CORRECTLY, and labelled three correct steps as slips.
SLIPS = [
    ('llama-3.1-70b', 'normal_depth_iteration#1', '(3.8 + 2*1.500*2)*1.500'),
    ('gpt-5', 'normal_depth_iteration#1', '1.853 − (0.2892 ×'),
    ('llama-3.1-70b', 'aoq_ati_rectifying#1', '(1 - 0.059)^13'),
    ('llama-3.1-70b', 'rackett_equation_volume#2', '0.276^(1 + 0.1745)'),
    ('llama-3.1-70b', 'fluid_particle_acceleration#1', '16(2.9)'),
    ('gemma-4-31b', 'aoq_ati_rectifying#0', '0.4525 \\cdot 0.089 \\cdot 80}{100} = 0.32218'),
    # deepseek-r1 normal_depth_iteration#2 is NOT used: the step states 5.071 and then
    # corrects itself in the same step ("2*1.8820 = 3.764 ... = 10.136"), so its label
    # is ambiguous and a judge could not be scored against it fairly.
    ('llama-3.1-70b', 'lorentz_force#0', '(-1924.8 + 3032.56)'),
    ('llama-3.1-70b', 'aoq_ati_rectifying#2', 'AOQ1 = 0.046 / (1 + (1+1) / (50-2+1))'),
    ('llama-3.1-70b', 'normal_depth_iteration#3', '6.846*0.500 / 8.494'),
    ('llama-3.1-70b', 'fluid_particle_acceleration#0', '4(3.4)(2.3)'),
    ('llama-3.1-70b', 'normal_depth_iteration#0', '4.4*(0.559)^(2/3) - 20.56'),
]
N_CLEAN = 12

# Labels corrected after reading, with the evidence. The CLEAN selection trusts E4's
# 1% arithmetic tolerance, and that tolerance passed a real error:
RELABEL = {
    ('llama-3.1-70b', 'damping_classification#3'): (
        'SLIP', 'writes sqrt(124343 * 582.96) = sqrt(72311151.28); the product is '
                '72486995, a 0.24% error inside both E4 (1%) and E3 (0.5%) tolerances. '
                'Every judge flagged it, correctly; the CLEAN label was wrong.'),
}
CATS = ('alternative correct', 'calculation error', 'conceptual error', 'other')


def _jsonl(path):
    return [json.loads(l) for l in open(path, encoding='utf-8')]


def _current(d):
    rows = [json.loads(l) for f in glob.glob(_os.path.join(_PILOT, 'scores', d, '*.jsonl'))
            for l in open(f, encoding='utf-8')]
    rows = [r for r in rows if not r.get('error')]
    latest = max(rows, key=lambda r: r['ts'])['config_sha256']
    return {(r['model_key'], r['item_id']): r for r in rows if r['config_sha256'] == latest}


def build():
    sys.path.insert(0, _os.path.join(_REPO, 'evaluation'))
    os_env = _os.environ
    os_env.setdefault('HF_HUB_OFFLINE', '1')
    with redirect_stdout(io.StringIO()):
        import engtrace_evaluation_framework as fw
    from engineering_parser import extract_steps

    items = {r['item_id']: r for r in _jsonl(_os.path.join(_PILOT, 'slice', 'manifest.jsonl'))}
    traces = {}
    for f in glob.glob(_os.path.join(_PILOT, 'traces', '*.jsonl')):
        k = _os.path.basename(f)[:-6]
        for r in _jsonl(f):
            if r['ok'] and r.get('finish_reason') in ('stop', 'end_turn', 'eos', None):
                traces[(k, r['item_id'])] = r

    def tribunal_prompt(question, gt_steps, pred_steps, idx):
        """The framework's own prompt, captured from its own method."""
        captured = {}

        class Fake:
            client_openai, client_anthropic = True, None

            def _call_single_judge(self, provider, prompt):
                captured['prompt'] = prompt
                return []
        fw.EngTraceFramework._tier2_tribunal_batch(Fake(), question, gt_steps, pred_steps, [idx])
        return captured['prompt']

    def step_of(steps, needle):
        hits = [i for i, s in enumerate(steps) if needle in s.replace(' ', '') or needle in s]
        if not hits:
            squash = needle.replace(' ', '')
            hits = [i for i, s in enumerate(steps) if squash in s.replace(' ', '')]
        return hits[0] if hits else None

    probe = []
    for mk, iid, needle in SLIPS:
        tr = traces.get((mk, iid))
        if not tr:
            print('  skip slip %s %s: no clean-finish trace' % (mk, iid))
            continue
        gt, pred = extract_steps(items[iid]['solution'])[0], extract_steps(tr['text'])[0]
        idx = step_of(pred, needle)
        if idx is None:
            print('  skip slip %s %s: needle %r not found in any step' % (mk, iid, needle))
            continue
        probe.append({'label': 'SLIP', 'model': mk, 'item_id': iid, 'step': idx,
                      'step_text': pred[idx][:400],
                      'prompt': tribunal_prompt(items[iid]['question'], gt, pred, idx)})

    e4, e0 = _current('e4'), _current('e0')
    clean = []
    for k, r in sorted(e4.items()):
        if k not in e0 or e0[k]['scores']['final_answer_acc'] != 1:
            continue
        m = r['meta']
        if m['claims_checked'] < 3 or m['claims_consistent'] != m['claims_checked']:
            continue
        verified = [h for h in m['milestones'] if h['status'] == 'verified']
        if not verified:
            continue
        clean.append((k, verified[0]['value']))
    # spread across models, deterministic
    by_model = {}
    for (k, v) in clean:
        by_model.setdefault(k[0], []).append((k, v))
    pick = []
    while len(pick) < N_CLEAN and any(by_model.values()):
        for mk in sorted(by_model):
            if by_model[mk] and len(pick) < N_CLEAN:
                pick.append(by_model[mk].pop(0))
    sys.path.insert(0, _os.path.join(_PILOT, 'evaluators'))
    import milestones as ms
    for (mk, iid), val in pick:
        tr = traces[(mk, iid)]
        gt, pred = extract_steps(items[iid]['solution'])[0], extract_steps(tr['text'])[0]
        idx = next((i for i, s in enumerate(pred) if ms.stated(val, ms.numbers(s))), None)
        if idx is None:
            continue
        probe.append({'label': 'CLEAN', 'model': mk, 'item_id': iid, 'step': idx,
                      'step_text': pred[idx][:400],
                      'prompt': tribunal_prompt(items[iid]['question'], gt, pred, idx)})

    _os.makedirs(OUT, exist_ok=True)
    json.dump(probe, open(PROBE, 'w', encoding='utf-8'), indent=1, ensure_ascii=False)
    print('probe set: %d SLIP, %d CLEAN -> %s' % (
        sum(p['label'] == 'SLIP' for p in probe), sum(p['label'] == 'CLEAN' for p in probe), PROBE))


def parse(raw):
    """The framework's own extraction logic (_call_single_judge), copied."""
    clean = re.sub(r'^```json\s*', '', raw or '', flags=re.MULTILINE)
    clean = re.sub(r'\s*```$', '', clean, flags=re.MULTILINE)
    try:
        data = json.loads(clean)
    except Exception:
        data = None
        key = clean.find('"results"')
        if key == -1:
            key = clean.find("'results'")
        if key != -1:
            start = clean.rfind('{', 0, key)
            if start != -1:
                snip = clean[start:]
                last = snip.rfind('}')
                while last != -1:
                    try:
                        data = json.loads(snip[:last + 1])
                        break
                    except Exception:
                        try:
                            import ast
                            data = ast.literal_eval(snip[:last + 1])
                            break
                        except Exception:
                            last = snip.rfind('}', 0, last)
    if data is None:
        return None
    res = data.get('results', data.get('steps', [])) if isinstance(data, dict) else data
    if isinstance(res, dict):
        res = [res]
    return res


def call(model, prompt):
    from openai import OpenAI
    from dotenv import load_dotenv
    load_dotenv(_os.path.join(_REPO, '.env'))
    cli = OpenAI(api_key=_os.environ['OPENROUTER_API_KEY'], base_url='https://openrouter.ai/api/v1',
                 timeout=300.0, max_retries=2)
    t = time.time()
    try:
        r = cli.chat.completions.create(model=model, max_tokens=8192,
                                        messages=[{'role': 'user', 'content': prompt}])
        u = r.usage
        return {'ok': True, 'text': r.choices[0].message.content or '', 'served': r.model,
                'finish': r.choices[0].finish_reason, 'in': u.prompt_tokens,
                'out': u.completion_tokens, 'seconds': round(time.time() - t, 1)}
    except Exception as exc:                                       # noqa: BLE001
        return {'ok': False, 'error': str(exc)[:300], 'seconds': round(time.time() - t, 1)}


def run():
    probe = json.load(open(PROBE, encoding='utf-8'))
    path = _os.path.join(OUT, 'replies.jsonl')
    done = set()
    if _os.path.exists(path):
        done = {(r['judge'], r['probe']) for r in _jsonl(path) if r.get('ok')}
    jobs = [(j, i) for j in JUDGES for i in range(len(probe)) if (j[0], i) not in done]
    print('%d calls to make' % len(jobs))

    def one(job):
        (jid, model, fam, role), i = job
        return dict(call(model, probe[i]['prompt']), judge=jid, model=model, probe=i)

    with cf.ThreadPoolExecutor(max_workers=12) as pool, open(path, 'a', encoding='utf-8') as fh:
        for res in pool.map(one, jobs):
            fh.write(json.dumps(res, ensure_ascii=False) + '\n')
            fh.flush()
    print('done')


def report():
    from collections import defaultdict
    probe = json.load(open(PROBE, encoding='utf-8'))
    for p in probe:
        fix = RELABEL.get((p['model'], p['item_id']))
        if fix and p['label'] != fix[0]:
            print('relabel %s %s: %s -> %s (%s)' % (p['model'], p['item_id'], p['label'], fix[0], fix[1][:90]))
            p['label'] = fix[0]
    replies = {}
    for r in _jsonl(_os.path.join(OUT, 'replies.jsonl')):
        if r.get('ok') or (r['judge'], r['probe']) not in replies:
            replies[(r['judge'], r['probe'])] = r
    import urllib.request
    from dotenv import load_dotenv
    load_dotenv(_os.path.join(_REPO, '.env'))
    cat = {m['id']: m['pricing'] for m in json.load(urllib.request.urlopen(urllib.request.Request(
        'https://openrouter.ai/api/v1/models',
        headers={'Authorization': 'Bearer ' + _os.environ['OPENROUTER_API_KEY'], 'User-Agent': 'engtrace'}),
        timeout=60))['data']}

    n_slip = sum(p['label'] == 'SLIP' for p in probe)
    n_clean = len(probe) - n_slip
    print('probe: %d SLIP steps (known arithmetic errors), %d CLEAN steps\n' % (n_slip, n_clean))
    print('%-15s %-22s %6s %7s %7s %8s %9s %8s %7s' % (
        'judge', 'family', 'parsed', 'caught', 'false', 'balanced', 'cost/call', 'med s', 'trunc'))
    rows = []
    for jid, model, fam, role in JUDGES:
        caught = false = parsed = n = trunc = 0
        cost, secs = 0.0, []
        for i, p in enumerate(probe):
            r = replies.get((jid, i))
            if not r:
                continue
            n += 1
            if not r.get('ok'):
                continue
            pin = float(cat.get(model, {}).get('prompt') or 0)
            pout = float(cat.get(model, {}).get('completion') or 0)
            cost += (r['in'] or 0) * pin + (r['out'] or 0) * pout
            secs.append(r['seconds'])
            trunc += r.get('finish') == 'length'
            res = parse(r['text'])
            if not res:
                continue
            verdict = [str(x.get('category', '')).lower() for x in res
                       if isinstance(x, dict) and x.get('step_index') == p['step']]
            if not verdict:
                continue
            parsed += 1
            flagged = 'alternative' not in verdict[0]
            if p['label'] == 'SLIP' and flagged:
                caught += 1
            if p['label'] == 'CLEAN' and flagged:
                false += 1
        tpr = caught / n_slip if n_slip else 0
        fpr = false / n_clean if n_clean else 0
        bal = (tpr + (1 - fpr)) / 2
        rows.append((jid, bal))
        print('%-15s %-22s %3d/%-2d %4d/%-2d %4d/%-2d %8.2f %9.4f %8.1f %7d' % (
            jid, fam[:22], parsed, n, caught, n_slip, false, n_clean, bal,
            cost / max(len(secs), 1), sorted(secs)[len(secs) // 2] if secs else 0, trunc))
    print('\ncaught = SLIP steps NOT called "Alternative Correct"; false = CLEAN steps flagged;')
    print('balanced = mean of catch rate and pass rate (0.5 = coin flip / rubber stamp).')


if __name__ == '__main__':
    {'build': build, 'run': run, 'report': report}[sys.argv[1]]()
