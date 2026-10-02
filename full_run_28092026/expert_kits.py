"""The experts' reading request for the full run: B1 to B4 of docs/EVALUATION_NEXT_STEPS.md, one kit per expert (D-172).

    python -m full_run_28092026.expert_kits --build [--seed 20261002] [--b1 150] [--b3 100] [--b2-per-model 40] [--b2-models K ...]
                                                                   # FREE: expert_request/tasks/ and expert_request/dist/, local
    python -m full_run_28092026.expert_kits --score DIR            # FREE: the returned <id>.jsonl files; EXPERT_REQUEST.md, counts only
    python -m full_run_28092026.expert_kits --selftest             # FREE: a small build in a temp folder, the app driven through it

WHAT IS DRAWN, from the main store as re-scored under D-169 and the traces it scored (score.texts_matching refuses a
trace the store did not score):

  B1  final answers    `answer` items. Half the check's incorrect-or-partial verdicts and half its correct ones, from
                       the top five models, each half taken round-robin over templates (one per template per pass, a
                       seeded shuffle within), so no template dominates and the forms the check reads are spread. The
                       expert sees the problem, the reference answer (the gold's Answer segment), the model's final
                       answer (the trace's Answer segment), and both in full one click away; not the model, not the
                       check's verdict.
  B3  milestones       `milestone` items. Half REACHED and half MISSING verdicts of E5's judge, round-robin over the
                       eleven models and, within a model, over templates; one milestone per trace. The expert sees
                       the problem, the milestone's name and value, the reference, and the full trace; not the judge's
                       verdict.
  B2  wrong answers    `error` items. Per model in --b2-models, --b2-per-model answered, readable traces the check
                       scored incorrect, allocated over Easy / Intermediate / Advanced in proportion to the model's
                       wrong answers (at least one per level where it has one), round-robin over templates within a
                       level. The expert sees the problem, the trace and the reference, and assigns one of the May
                       version's six categories by its stop-at-first-yes hierarchy, or "No error" or "Incomplete".
  B4  templates        `template` items, chemical and civil only. The two Advanced chemical templates of D-171, each
                       with one problem and three top-model answers to it, and six templates whose wrong answers
                       cluster within 5% of the reference, each with one problem and one such answer; two to four
                       questions each on whether the wording pins a single answer at the check's tolerance.

WHO READS WHAT. The experts are the layer-2 roster (the pilot's 15, three per branch). An item goes to experts of
its own branch: `answer` and `milestone` items to two of the three (rotating pairs, so every pair occurs), `error`
and `template` items to all three. Codes are opaque; the keyfile from code to item, model and hidden verdict stays
local.

HOW IT IS DELIVERED. One bundle per expert, expert_request/dist/<id>/, with a sub-folder per task:
B1_final_answers, B3_milestones, B2_wrong_answers and, for chemical and civil, B4_templates. Each sub-folder is
complete on its own: app.py (the same reading_app.py, which takes its title from the one kind it finds), README.txt,
guide.md (the general part of EXPERT_READING_GUIDE.md plus that task's section) and tasks/ with the expert's items
of that kind, shuffled for that expert from the seed. Each produces its own <id>.jsonl, so the tasks can be sent or
dropped separately; --score reads a folder tree.

WHAT IS KEPT. expert_request/ (tasks, dist, build.json, scored.json, returns) holds pool text and the experts'
labels and is gitignored. Committed: this script, reading_app.py, EXPERT_READING_GUIDE.md and EXPERT_REQUEST.md,
which --build writes with the composition (counts) and --score rewrites with the results (counts; no item text, no
expert's words).

WHAT --score REPORTS. B1: agreement between the experts' verdict and the check's, three-way and on the check's
correct and incorrect sides, per branch, and between the two readers (Cohen's kappa). B3: the share of REACHED
milestones the experts say the working obtains, and of MISSING ones they say it does not, with the "not needed"
share, and the readers' agreement. B2: the category distribution per model and level, item majorities, the share of
"No error", and Fleiss' kappa over the three readers. B4: per template and question, the three experts' answers as
counts. The last submission per reader and item counts.
"""
from __future__ import annotations

import argparse
import collections
import datetime as dt
import hashlib
import json
import random
import shutil
import statistics
import sys
import tempfile
from pathlib import Path

HERE = Path(__file__).resolve().parent
REPO = HERE.parent
for p in (str(REPO), str(REPO / 'evaluator_pilot_17092026' / 'evaluators'), str(REPO / 'evaluation')):
    if p not in sys.path:
        sys.path.insert(0, p)

import answer  # noqa: E402

from full_run_28092026 import residual_incorrect as ri  # noqa: E402
from full_run_28092026 import score  # noqa: E402
from full_run_28092026.analyze import ROSTER  # noqa: E402

OUT = HERE / 'expert_request'
REPORT = HERE / 'EXPERT_REQUEST.md'
GUIDE = HERE / 'EXPERT_READING_GUIDE.md'
APP = HERE / 'reading_app.py'
SALT = 'engtrace-expert-request-2026-10'
SEED = 20261002
TOP5 = ['deepseek-v4.1-flash', 'kimi-k3', 'claude-sonnet-5', 'glm-5.3-flash', 'muse-glimmer-30b']
B2_MODELS = ['claude-sonnet-5', 'gpt-5.4-mini', 'gpt-oss-20b', 'gemma-4-26b-a4b']
LEVELS = ['Easy', 'Intermediate', 'Advanced']
READERS = {'template': 3, 'error': 3, 'answer': 2, 'milestone': 2}
ORDER = ['template', 'answer', 'milestone', 'error']
FOLDER = {'answer': 'B1_final_answers', 'error': 'B2_wrong_answers', 'milestone': 'B3_milestones',
          'template': 'B4_templates'}
TASK = {'answer': "B1, final answers: is a model's final answer correct? A minute or two each.",
        'milestone': "B3, milestones: does a model's working obtain a given quantity? Two or three minutes each.",
        'error': 'B2, wrong answers: why did a wrong answer go wrong? Three to five minutes each.',
        'template': "B4, templates: does a problem's wording pin a single answer? Five to ten minutes each."}
B4_MAIN = {'template_work_isothermal_virial': 'virial', 'template_adiabatic_flame_temperature': 'flame'}
B4_NEAR = ['template_vdw_solve_for_volume', 'template_pfr_volume_changing_rate', 'template_annulus_flowrate',
           'template_pitzer_correlation_z', 'template_best_hydraulic_rectangular_section',
           'template_manning_rectangular_discharge']
NEAR = (0.002, 0.05)

ANSWER_OPTIONS = ['correct', 'partially correct', 'incorrect', 'no answer stated']
MILESTONE_OPTIONS = ['yes, the working obtains it', 'no, it never obtains it', 'no, and its route does not need it']
ERROR_OPTIONS = ['1. Hallucination', '2. Setup / Assumption Error', '3. Formula / Principle Error',
                 '4. Unit / Dimensional Error', '5. Sign / Direction Error', '6. Calculation Error',
                 'No error: the answer is correct, or the question admits it',
                 'Incomplete: the working stops before an answer']
EXPERT_TO_CHECK = {'correct': 'correct', 'partially correct': 'partial', 'incorrect': 'incorrect',
                   'no answer stated': 'incorrect'}

README_TASK = """EngTrace reading request, October 2026 - {task}

This folder holds one task: app.py, this README, guide.md, and a tasks folder with your items.

1. Read guide.md first: its last section but one describes this task.

2. Run the app. You need Python 3.11 or newer. In this folder:

       pip install streamlit
       streamlit run app.py

   Your browser opens the app. Pick your id in the sidebar.

3. Work through your items. The app saves after every item. To pause, close the browser tab and
   stop the app; run "streamlit run app.py" again to continue where you left off.

4. When the app says you are done, send the file <your id>.jsonl, which the app has written in
   this folder next to app.py, to the coordinator.

Questions about the task go to the coordinator, not to the other reviewers. Please do not share
the folder: the problems are the benchmark's private test set until the paper is published.
"""

README_EXPERT = """EngTrace reading request, October 2026 - your tasks

Each sub-folder is one task, with its own app, guide and items; each returns its own file:

{folders}

Do them in the order listed if you can. Every sub-folder's README says how to run it; the guides
share their first two sections and differ in the task's own.
"""


def guide_for(kind: str) -> str:
    """The general part of the guide plus the one section for this task."""
    head, *parts = GUIDE.read_text(encoding='utf-8').split('\n## ')
    sections = {p.split('\n', 1)[0].strip(): p for p in parts}
    own = next(t for t in sections if t.startswith(FOLDER[kind][:2] + ':'))
    return head.rstrip() + '\n\n' + '\n\n'.join('## ' + sections[t].rstrip() for t in
                                                ('What this is', 'How to run', own, 'Please')) + '\n'


# ------------------------------------------------------------------ inputs

def experts() -> list[dict]:
    from template_annotation_23092026.layer2.build_tasks import roster as layer2_roster
    return [{'id': r['id'], 'branch': r['branch']} for r in layer2_roster()]


def load(models: list[str]):
    """Items, milestones, store rows, trace texts and E5 rows for the models given."""
    items = score.pool_items()
    ms_all = score.milestone_sets(items)
    rows, texts, e5 = {}, {}, {}
    for key in models:
        rows[key] = {r['item_id']: r for r in map(json.loads, (score.SCORES / 'main' / f'{key}.jsonl')
                                                   .read_text(encoding='utf-8').splitlines())}
        texts[key] = score.texts_matching('main', key, rows[key].values())
        p = score.SCORES / 'main' / 'e5' / f'{key}.jsonl'
        e5[key] = {r['item_id']: r for r in map(json.loads, p.read_text(encoding='utf-8').splitlines())} \
            if p.exists() else {}
    return items, ms_all, rows, texts, e5


def code_of(*parts: str) -> str:
    return 'ER-' + hashlib.sha256('|'.join((SALT,) + parts).encode('utf-8')).hexdigest()[:8]


def answered(r: dict) -> bool:
    return r['status'] == 'answered' and not r['unusable']


# ------------------------------------------------------------------ drawing

def round_robin(groups: dict, n: int, rng: random.Random, tag=None) -> list:
    """Up to n members, one per group per pass, groups in a seeded order and each group's members shuffled;
    `tag(member)` names a key of which at most one member is taken (one milestone per trace)."""
    keys = sorted(groups)
    rng.shuffle(keys)
    pools = {k: rng.sample(list(groups[k]), len(groups[k])) for k in keys}
    out, used = [], set()
    while len(out) < n and any(pools.values()):
        took = False
        for k in keys:
            while pools[k] and len(out) < n:
                x = pools[k].pop()
                t = tag(x) if tag else None
                if t is not None:
                    if t in used:
                        continue
                    used.add(t)
                out.append(x)
                took = True
                break
        if not took:
            break
    return out


def allocate(counts: dict, n: int) -> dict:
    """n split over the keys in proportion to counts, at least one where a key has any, never above its count."""
    total = sum(counts.values())
    if total <= n:
        return dict(counts)
    raw = {k: n * c / total for k, c in counts.items()}
    alloc = {k: min(counts[k], int(raw[k])) for k in counts}
    for k in counts:
        if counts[k] and not alloc[k]:
            alloc[k] = 1
    order = sorted(counts, key=lambda k: raw[k] - int(raw[k]), reverse=True)
    i = 0
    while sum(alloc.values()) < n and any(alloc[k] < counts[k] for k in counts):
        k = order[i % len(order)]
        if alloc[k] < counts[k]:
            alloc[k] += 1
        i += 1
    while sum(alloc.values()) > n:
        k = max(alloc, key=lambda k: alloc[k])
        alloc[k] -= 1
    return alloc


def draw_b1(items, rows, texts, n: int, rng: random.Random) -> list[dict]:
    wrong, right = collections.defaultdict(list), collections.defaultdict(list)
    for key in TOP5:
        for i, r in rows[key].items():
            if answered(r):
                (right if r['answer']['label'] == 'correct' else wrong)[r['template_id']].append((key, i))
    picked = round_robin(wrong, n // 2, rng) + round_robin(right, n - n // 2, rng)
    out = []
    for key, i in picked:
        it, r, text = items[i], rows[key][i], texts[key][i]
        out.append({'pool': {'kind': 'answer', 'branch': it['branch'], 'question': it['question'],
                             'reference': answer.segment(it['solution']).strip(), 'solution': it['solution'],
                             'stated': answer.segment(text).strip(), 'trace': text},
                    'key': {'kind': 'answer', 'item_id': i, 'model': key, 'template_id': it['template_id'],
                            'answer_type': it['answer_type'], 'level': it['level'], 'check': r['answer']['label']},
                    'code': code_of('answer', key, i)})
    return out


def draw_b3(items, ms_all, rows, texts, e5, n: int, rng: random.Random) -> list[dict]:
    out = []
    for verdict in ('REACHED', 'MISSING'):
        per_model = {}
        for key in ROSTER:
            by_t = collections.defaultdict(list)
            for i, e in e5.get(key, {}).items():
                if not e.get('sent') or not answered(rows[key][i]):
                    continue
                for idx, s in enumerate(e['sources']):
                    if s == verdict:
                        by_t[items[i]['template_id']].append((key, i, idx))
            if by_t:
                per_model[key] = round_robin(by_t, sum(map(len, by_t.values())), rng, tag=lambda x: (x[0], x[1]))
        picked = round_robin({k: v for k, v in per_model.items()}, n // 2, rng, tag=lambda x: (x[0], x[1]))
        for key, i, idx in picked:
            it, m = items[i], ms_all[i][idx]
            out.append({'pool': {'kind': 'milestone', 'branch': it['branch'], 'question': it['question'],
                                 'solution': it['solution'], 'milestone_id': m['id'],
                                 'milestone_value': f"{m['value']:.6g}", 'trace': texts[key][i]},
                        'key': {'kind': 'milestone', 'item_id': i, 'model': key, 'template_id': it['template_id'],
                                'index': idx, 'judge': verdict},
                        'code': code_of('milestone', key, i, str(idx))})
    return out


def draw_b2(items, rows, texts, models: list[str], per_model: int, rng: random.Random) -> list[dict]:
    out = []
    for key in models:
        cands = {lv: collections.defaultdict(list) for lv in LEVELS}
        for i, r in rows[key].items():
            if answered(r) and r['answer']['label'] == 'incorrect':
                cands[r['level']][r['template_id']].append(i)
        alloc = allocate({lv: sum(map(len, cands[lv].values())) for lv in LEVELS}, per_model)
        for lv in LEVELS:
            for i in round_robin(cands[lv], alloc.get(lv, 0), rng):
                it = items[i]
                out.append({'pool': {'kind': 'error', 'branch': it['branch'], 'question': it['question'],
                                     'solution': it['solution'], 'trace': texts[key][i], 'level': lv},
                            'key': {'kind': 'error', 'item_id': i, 'model': key, 'template_id': it['template_id'],
                                    'level': lv},
                            'code': code_of('error', key, i)})
    return out


def distances(rows, texts, template_id: str) -> tuple[list[tuple], int, int]:
    """(model, item, closest) for the template's answered incorrect traces, the rows on it, the unusable ones."""
    wrong, total, unusable = [], 0, 0
    for key in ROSTER:
        for i, r in rows[key].items():
            if r['template_id'] != template_id:
                continue
            total += 1
            if r['unusable'] or r['status'] != 'answered':
                unusable += 1
            elif r['answer']['label'] == 'incorrect':
                t = r['answer']['targets']
                d = ri.closest(t['numbers'], answer.values(answer.segment(texts[key][i]))) if t['numbers'] else None
                wrong.append((key, i, d))
    return wrong, total, unusable


def near(d) -> bool:
    return d is not None and NEAR[0] < d <= NEAR[1]


def draw_b4(items, rows, texts) -> list[dict]:
    out = []
    virial = ri.virial_readings(ROSTER, items)
    for tid, role in list(B4_MAIN.items()) + [(t, 'near') for t in B4_NEAR]:
        wrong, total, unusable = distances(rows, texts, tid)
        within = sum(near(d) or (d is not None and d <= NEAR[0]) for _, _, d in wrong)
        short = tid.replace('template_', '')
        if role == 'near':
            cands = sorted(((0 if k in TOP5 else 1, d, i, k) for k, i, d in wrong if near(d)))
            if not cands:
                continue
            _, _, i, k = cands[0]
            chosen, traces = i, [(k, texts[k][i])]
            context = (f'Template {short}. In the run, {len(wrong)} of the {total} answers to its 15 problems were '
                       f'scored wrong ({unusable} more were empty or unreadable); {within} of the {len(wrong)} lie within '
                       f'5% of the reference value. The model answer below is one of those.')
            asks = [{'key': 'unique', 'text': 'Does the question, as worded, pin a single correct answer to within 0.2%?',
                     'options': ['yes', 'no', 'cannot tell']},
                    {'key': 'trace', 'text': 'Is the model answer shown acceptable as a correct answer to the question as worded?',
                     'options': ['yes', 'no']}]
            note = 'Note: if no, what is left open (a method, a data source, a convention, a rounding)?'
        else:
            per_item = collections.Counter(i for k, i, d in wrong if k in TOP5)
            chosen = min(per_item, key=lambda i: (-per_item[i], items[i]['instance_index'])) if per_item \
                else min({i for _, i, _ in wrong}, key=lambda i: items[i]['instance_index'])
            ranked = sorted(((0 if near(d) else 1, 0 if k in TOP5 else 1, k) for k, i, d in wrong if i == chosen))
            traces, seen = [], set()
            for _, _, k in ranked:
                if k not in seen and len(traces) < 3:
                    seen.add(k)
                    traces.append((k, texts[k][chosen]))
            if role == 'virial':
                context = (f'Template {short}. In the run, {len(wrong)} of the {total} answers to its 15 problems were '
                           f'scored wrong ({unusable} more were empty). Of the {len(wrong)} wrong answers, '
                           f'{virial["flow work, pressure-explicit form"]} state the value RT ln(P2/P1) + B (P2 - P1), the '
                           'flow work (integral of V dP) under Z = 1 + BP/RT; the reference computes the closed-system '
                           'work (minus the integral of P dV) under Z = 1 + B/V. The three model answers below are to '
                           'the same problem.')
                asks = [{'key': 'reading', 'text': "Which work does the question's wording ask for?",
                         'options': ['the closed-system work, -∫P dV', 'the flow (shaft) work, ∫V dP',
                                     'either: the wording does not decide', 'neither is clear']},
                        {'key': 'form', 'text': 'Which truncation does "the virial equation of state truncated to two '
                                                'terms" mean here?',
                         'options': ['Z = 1 + B/V (volume-explicit)', 'Z = 1 + BP/RT (pressure-explicit)',
                                     'either: the wording does not decide']},
                        {'key': 'unique', 'text': 'Is the reference answer the only correct answer to within 0.2%?',
                         'options': ['yes', 'no']},
                        {'key': 'traces', 'text': 'Are the model answers shown acceptable as correct solutions of the '
                                                  'question as worded?',
                         'options': ['all of them', 'some of them', 'none of them']}]
            else:
                context = (f'Template {short}. In the run, {len(wrong)} of the {total} answers to its 15 problems were '
                           f'scored wrong ({unusable} more were empty); {within} of the {len(wrong)} lie within 5% of '
                           'the reference value. The question states no heat-capacity data; the reference iterates a '
                           'polynomial heat capacity with coefficients it does not cite. The three model answers below '
                           'are to the same problem.')
                asks = [{'key': 'data', 'text': 'The question gives no heat-capacity data. Can a correct solution be '
                                                'expected to match the reference to within 0.2%?',
                         'options': ['yes: standard data lead to the same value',
                                     'no: standard sources differ by more than that',
                                     'only if the same data source is used']},
                        {'key': 'method', 'text': "Is the reference's method and data a standard textbook choice?",
                         'options': ['yes', 'no']},
                        {'key': 'traces', 'text': 'Are the model answers shown acceptable as correct estimates?',
                         'options': ['all of them', 'some of them', 'none of them']}]
            note = 'Note: what wording, or what data in the question, would remove the ambiguity, if you see one?'
        it = items[chosen]
        out.append({'pool': {'kind': 'template', 'branch': it['branch'], 'template': short, 'context': context,
                             'question': it['question'], 'solution': it['solution'],
                             'traces': [t for _, t in traces], 'asks': asks, 'note_prompt': note},
                    'key': {'kind': 'template', 'template_id': tid, 'item_id': chosen, 'role': role,
                            'models': [k for k, _ in traces], 'wrong': len(wrong), 'within_5pct': within},
                    'code': code_of('template', tid, chosen)})
    return out


# ------------------------------------------------------------------ the build

def assign(tasks: list[dict], roster: list[dict], seed: int) -> dict:
    by_branch = collections.defaultdict(list)
    for e in roster:
        by_branch[e['branch']].append(e['id'])
    queues = {e['id']: {'branch': e['branch'], 'codes': []} for e in roster}
    per_kind = collections.defaultdict(list)
    for t in sorted(tasks, key=lambda t: t['code']):
        per_kind[(t['pool']['branch'], t['pool']['kind'])].append(t['code'])
    for (branch, kind), codes in sorted(per_kind.items()):
        ids = sorted(by_branch[branch])
        if not ids:
            raise SystemExit(f'no expert for branch {branch}')
        k = READERS[kind]
        for j, c in enumerate(codes):
            for r in range(min(k, len(ids))):
                queues[ids[(j + r) % len(ids)]]['codes'].append(c)
    kind_of = {t['code']: t['pool']['kind'] for t in tasks}
    for aid, q in queues.items():
        rng = random.Random(f'{seed}|{aid}')
        ordered = []
        for kind in ORDER:
            block = [c for c in q['codes'] if kind_of[c] == kind]
            rng.shuffle(block)
            ordered += block
        q['codes'] = ordered
    return queues


def build(out: Path, seed: int, b1: int, b3: int, b2_per_model: int, b2_models: list[str],
          report: Path | None) -> dict:
    rng = random.Random(seed)
    models = sorted(set(ROSTER) | set(b2_models))
    items, ms_all, rows, texts, e5 = load(models)
    tasks = draw_b4(items, rows, texts) + draw_b1(items, rows, texts, b1, rng) \
        + draw_b3(items, ms_all, rows, texts, e5, b3, rng) + draw_b2(items, rows, texts, b2_models, b2_per_model, rng)
    codes = [t['code'] for t in tasks]
    if len(set(codes)) != len(codes):
        raise SystemExit('code collision')
    roster = experts()
    queues = assign(tasks, roster, seed)
    pool = {t['code']: t['pool'] for t in tasks}
    keyfile = {t['code']: t['key'] for t in tasks}
    tasks_dir, dist = out / 'tasks', out / 'dist'
    for d in (tasks_dir, dist):
        if d.exists():
            shutil.rmtree(d)
        d.mkdir(parents=True)
    (tasks_dir / 'pool.json').write_text(json.dumps(pool, ensure_ascii=False), encoding='utf-8')
    (tasks_dir / 'assignment.json').write_text(json.dumps(queues, indent=1), encoding='utf-8')
    (tasks_dir / 'keyfile.json').write_text(json.dumps(keyfile, indent=1), encoding='utf-8')
    guides = {k: guide_for(k) for k in FOLDER}
    for aid, q in queues.items():
        listed = []
        for kind in ORDER:
            mine = [c for c in q['codes'] if pool[c]['kind'] == kind]
            if not mine:
                continue
            folder = dist / aid / FOLDER[kind]
            (folder / 'tasks').mkdir(parents=True)
            shutil.copy(APP, folder / 'app.py')
            (folder / 'guide.md').write_text(guides[kind], encoding='utf-8')
            (folder / 'README.txt').write_text(README_TASK.format(task=TASK[kind]), encoding='utf-8')
            (folder / 'tasks' / 'pool.json').write_text(json.dumps({c: pool[c] for c in mine}, ensure_ascii=False),
                                                       encoding='utf-8')
            (folder / 'tasks' / 'assignment.json').write_text(
                json.dumps({aid: {'branch': q['branch'], 'codes': mine}}, indent=1), encoding='utf-8')
            listed.append(f'  {FOLDER[kind]}/   {len(mine)} items.  {TASK[kind]}')
        (dist / aid / 'README.txt').write_text(README_EXPERT.format(folders='\n'.join(listed)), encoding='utf-8')
    comp = composition(tasks, queues, seed, b2_models)
    (out / 'build.json').write_text(json.dumps(comp, indent=1), encoding='utf-8')
    if report:
        report.write_text('\n'.join(report_lines(comp, None)) + '\n', encoding='utf-8', newline='\n')
    return comp


def composition(tasks: list[dict], queues: dict, seed: int, b2_models: list[str]) -> dict:
    c = {'built': dt.date.today().isoformat(), 'seed': seed, 'store_commit': store_commit(), 'items': len(tasks),
         'readings': sum(len(q['codes']) for q in queues.values()), 'kinds': {}, 'experts': {}}
    by_kind = collections.defaultdict(list)
    for t in tasks:
        by_kind[t['pool']['kind']].append(t)
    a = by_kind['answer']
    c['kinds']['answer'] = {
        'items': len(a), 'readers': READERS['answer'],
        'by_check': dict(collections.Counter(t['key']['check'] for t in a)),
        'by_branch': dict(collections.Counter(t['pool']['branch'] for t in a)),
        'by_model': dict(collections.Counter(t['key']['model'] for t in a)),
        'by_answer_type': dict(collections.Counter(t['key']['answer_type'] for t in a)),
        'templates': len({t['key']['template_id'] for t in a})}
    m = by_kind['milestone']
    c['kinds']['milestone'] = {
        'items': len(m), 'readers': READERS['milestone'],
        'by_judge': dict(collections.Counter(t['key']['judge'] for t in m)),
        'by_branch': dict(collections.Counter(t['pool']['branch'] for t in m)),
        'by_model': dict(collections.Counter(t['key']['model'] for t in m)),
        'templates': len({t['key']['template_id'] for t in m})}
    e = by_kind['error']
    c['kinds']['error'] = {
        'items': len(e), 'readers': READERS['error'], 'models': b2_models,
        'by_model_level': {k: dict(collections.Counter(t['key']['level'] for t in e if t['key']['model'] == k))
                           for k in b2_models},
        'by_branch': dict(collections.Counter(t['pool']['branch'] for t in e)),
        'templates': len({t['key']['template_id'] for t in e})}
    t4 = by_kind['template']
    c['kinds']['template'] = {
        'items': len(t4), 'readers': READERS['template'],
        'templates': {t['key']['template_id']: {'role': t['key']['role'], 'branch': t['pool']['branch'],
                                                'answers_shown': len(t['pool']['traces']),
                                                'wrong_in_run': t['key']['wrong'], 'within_5pct': t['key']['within_5pct']}
                      for t in t4}}
    kind_of = {t['code']: t['pool']['kind'] for t in tasks}
    for aid, q in sorted(queues.items()):
        c['experts'][aid] = {'branch': q['branch'], **dict(collections.Counter(kind_of[x] for x in q['codes'])),
                             'total': len(q['codes'])}
    return c


def store_commit() -> str:
    try:
        return json.loads((score.SCORES / 'main' / 'CONFIG.json').read_text(encoding='utf-8')).get('git', '?')[:7]
    except Exception:  # noqa: BLE001
        return '?'


def report_lines(comp: dict, results: dict | None) -> list[str]:
    k = comp['kinds']
    L = ['# The experts\' reading request (B1 to B4)', '',
         f"Generated by `expert_kits.py`: `--build` writes what was sent, `--score` what came back. Counts only; the "
         f"kits, the keyfile and the experts' files stay local. Built {comp['built']} from the main store at "
         f"`{comp['store_commit']}`, seed {comp['seed']}. Delivered as one bundle per expert with a sub-folder per task "
         f"(`B1_final_answers`, `B3_milestones`, `B2_wrong_answers`, `B4_templates`), each with its own app, guide and "
         f"return file.", '',
         '## What was sent', '',
         '| kind | items | readers per item | readings | composition |', '|---|---:|---:|---:|---|']
    a, m, e, t = k['answer'], k['milestone'], k['error'], k['template']
    def fmt(d):
        return ', '.join(f'{x} {n}' for x, n in sorted(d.items(), key=lambda kv: (-kv[1], kv[0])))
    L.append(f"| B1 final answers | {a['items']} | {a['readers']} | {a['items'] * a['readers']} | by the check's verdict: "
             f"{fmt(a['by_check'])}; by branch: {fmt({b.replace('_engineering', ''): n for b, n in a['by_branch'].items()})}; "
             f"by answer type: {fmt(a['by_answer_type'])}; {a['templates']} templates; models: {fmt(a['by_model'])} |")
    L.append(f"| B3 milestones | {m['items']} | {m['readers']} | {m['items'] * m['readers']} | by the judge's verdict: "
             f"{fmt(m['by_judge'])}; by branch: {fmt({b.replace('_engineering', ''): n for b, n in m['by_branch'].items()})}; "
             f"{m['templates']} templates; models: {fmt(m['by_model'])} |")
    L.append(f"| B2 wrong answers | {e['items']} | {e['readers']} | {e['items'] * e['readers']} | per model and level: "
             + '; '.join(f"{mk}: {fmt(lv)}" for mk, lv in e['by_model_level'].items())
             + f"; by branch: {fmt({b.replace('_engineering', ''): n for b, n in e['by_branch'].items()})}; {e['templates']} templates |")
    L.append(f"| B4 templates | {t['items']} | {t['readers']} | {t['items'] * t['readers']} | "
             + '; '.join(f"`{tid.replace('template_', '')}` ({v['role']}, {v['answers_shown']} answer{'s' if v['answers_shown'] > 1 else ''} shown; "
                         f"{v['wrong_in_run']} wrong in the run, {v['within_5pct']} within 5%)" for tid, v in t['templates'].items()) + ' |')
    L += ['', f"Items {comp['items']}, readings {comp['readings']}.", '', '| expert | branch | template | answer | milestone | error | total |',
          '|---|---|---:|---:|---:|---:|---:|']
    for aid, q in comp['experts'].items():
        L.append(f"| `{aid}` | {q['branch'].replace('_engineering', '')} | {q.get('template', 0)} | {q.get('answer', 0)} | "
                 f"{q.get('milestone', 0)} | {q.get('error', 0)} | {q['total']} |")
    L.append('')
    if results:
        L += results_lines(results)
    return L


# ------------------------------------------------------------------ scoring the returns

def cohen(pairs: list[tuple]) -> float | None:
    n = len(pairs)
    if n < 2:
        return None
    po = sum(a == b for a, b in pairs) / n
    cats = {x for p in pairs for x in p}
    pe = sum((sum(a == c for a, _ in pairs) / n) * (sum(b == c for _, b in pairs) / n) for c in cats)
    return (po - pe) / (1 - pe) if pe < 1 else None


def fleiss(rows: list[collections.Counter]) -> float | None:
    rows = [r for r in rows if sum(r.values()) >= 2]
    if len(rows) < 2:
        return None
    k = min(sum(r.values()) for r in rows)
    cats = {c for r in rows for c in r}
    n = len(rows)
    p_i = []
    for r in rows:
        tot = sum(r.values())
        p_i.append((sum(v * (v - 1) for v in r.values())) / (tot * (tot - 1)))
    p_bar = sum(p_i) / n
    tot_all = sum(sum(r.values()) for r in rows)
    p_j = {c: sum(r[c] for r in rows) / tot_all for c in cats}
    pe = sum(v * v for v in p_j.values())
    return (p_bar - pe) / (1 - pe) if pe < 1 and k >= 2 else None


def read_returns(returned: Path, out: Path) -> tuple[dict, dict, dict]:
    keyfile = json.loads((out / 'tasks' / 'keyfile.json').read_text(encoding='utf-8'))
    queues = json.loads((out / 'tasks' / 'assignment.json').read_text(encoding='utf-8'))
    owners = collections.defaultdict(set)
    for aid, q in queues.items():
        for c in q['codes']:
            owners[c].add(aid)
    rows = {}
    for f in sorted(returned.rglob('*.jsonl')):
        for line in f.read_text(encoding='utf-8').splitlines():
            if not line.strip():
                continue
            r = json.loads(line)
            if r['code'] not in keyfile:
                raise SystemExit(f"{f.name}: code {r['code']} is not in this build's keyfile")
            if r['annotator_id'] not in owners[r['code']]:
                raise SystemExit(f"{f.name}: {r['code']} was not assigned to {r['annotator_id']}")
            rows[(r['annotator_id'], r['code'])] = r
    return keyfile, queues, rows


def score_returns(returned: Path, out: Path, report: Path | None) -> dict:
    keyfile, queues, rows = read_returns(returned, out)
    by_code = collections.defaultdict(dict)
    for (aid, c), r in rows.items():
        by_code[c][aid] = r
    res = {'returned_readings': len(rows), 'assigned_readings': sum(len(q['codes']) for q in queues.values()),
           'files': sorted({aid for aid, _ in rows}), 'kinds': {}}
    # B1
    agree = collections.Counter()
    by_check = collections.defaultdict(collections.Counter)
    by_branch = collections.defaultdict(collections.Counter)
    dis_t = collections.Counter()
    pairs, notes = [], []
    for c, readers in by_code.items():
        k = keyfile[c]
        if k['kind'] != 'answer':
            continue
        labs = {}
        for aid, r in readers.items():
            v = r['answers']['verdict']
            labs[aid] = v
            e = EXPERT_TO_CHECK[v]
            by_check[k['check']][e] += 1
            same = e == k['check']
            agree['same' if same else 'different'] += 1
            by_branch[queues[aid]['branch']]['same' if same else 'different'] += 1
            if not same:
                dis_t[k['template_id']] += 1
            if r['note']:
                notes.append({'code': c, 'reader': aid, 'check': k['check'], 'expert': v, 'template': k['template_id'],
                              'model': k['model'], 'note': r['note']})
        if len(labs) == 2:
            pairs.append(tuple(labs[a] for a in sorted(labs)))
    res['kinds']['answer'] = {
        'readings': sum(agree.values()), 'agreement': agree['same'] / sum(agree.values()) if agree else None,
        'by_check': {ck: dict(v) for ck, v in by_check.items()},
        'by_branch': {b: dict(v) for b, v in by_branch.items()},
        'disagreements_by_template': dict(dis_t.most_common()),
        'items_read_twice': len(pairs), 'reader_agreement': (sum(a == b for a, b in pairs) / len(pairs)) if pairs else None,
        'reader_kappa': cohen(pairs), 'notes': notes}
    # B3
    by_judge = collections.defaultdict(collections.Counter)
    pairs, notes = [], []
    for c, readers in by_code.items():
        k = keyfile[c]
        if k['kind'] != 'milestone':
            continue
        labs = {}
        for aid, r in readers.items():
            v = r['answers']['obtained']
            labs[aid] = v
            by_judge[k['judge']][v] += 1
            if r['note']:
                notes.append({'code': c, 'reader': aid, 'judge': k['judge'], 'expert': v, 'template': k['template_id'],
                              'model': k['model'], 'note': r['note']})
        if len(labs) == 2:
            pairs.append(tuple(labs[a] for a in sorted(labs)))
    yes, no1, no2 = MILESTONE_OPTIONS
    rj = by_judge['REACHED']
    mj = by_judge['MISSING']
    res['kinds']['milestone'] = {
        'readings': sum(sum(v.values()) for v in by_judge.values()),
        'by_judge': {j: dict(v) for j, v in by_judge.items()},
        'reached_precision': rj[yes] / sum(rj.values()) if rj else None,
        'missing_precision': (mj[no1] + mj[no2]) / sum(mj.values()) if mj else None,
        'missing_not_needed_share': mj[no2] / sum(mj.values()) if mj else None,
        'items_read_twice': len(pairs), 'reader_agreement': (sum(a == b for a, b in pairs) / len(pairs)) if pairs else None,
        'reader_kappa': cohen(pairs), 'notes': notes}
    # B2
    per_model = collections.defaultdict(collections.Counter)
    per_level = collections.defaultdict(collections.Counter)
    majority = collections.defaultdict(collections.Counter)
    item_rows = collections.defaultdict(list)
    notes = []
    for c, readers in by_code.items():
        k = keyfile[c]
        if k['kind'] != 'error':
            continue
        cnt = collections.Counter(r['answers']['category'] for r in readers.values())
        for cat, n in cnt.items():
            per_model[k['model']][cat] += n
            per_level[k['level']][cat] += n
        top, n = cnt.most_common(1)[0]
        majority[k['model']][top if n >= 2 else 'split'] += 1
        item_rows[k['model']].append(cnt)
        for aid, r in readers.items():
            if r['note'] or r['excerpt']:
                notes.append({'code': c, 'reader': aid, 'model': k['model'], 'template': k['template_id'],
                              'category': r['answers']['category'], 'excerpt': r['excerpt'], 'note': r['note']})
    res['kinds']['error'] = {
        'readings': sum(sum(v.values()) for v in per_model.values()),
        'by_model': {mk: dict(v) for mk, v in per_model.items()},
        'by_level': {lv: dict(v) for lv, v in per_level.items()},
        'majority_by_model': {mk: dict(v) for mk, v in majority.items()},
        'no_error_share_by_model': {mk: v[ERROR_OPTIONS[6]] / sum(v.values()) for mk, v in per_model.items() if sum(v.values())},
        'fleiss_by_model': {mk: fleiss(v) for mk, v in item_rows.items()},
        'fleiss_all': fleiss([r for v in item_rows.values() for r in v]), 'notes': notes}
    # B4
    t4 = {}
    notes = []
    for c, readers in by_code.items():
        k = keyfile[c]
        if k['kind'] != 'template':
            continue
        asks = collections.defaultdict(collections.Counter)
        for aid, r in readers.items():
            for q, v in r['answers'].items():
                asks[q][v] += 1
            if r['note']:
                notes.append({'code': c, 'reader': aid, 'template': k['template_id'], 'note': r['note']})
        t4[k['template_id']] = {'role': k['role'], 'readers': len(readers), 'asks': {q: dict(v) for q, v in asks.items()}}
    res['kinds']['template'] = {'templates': t4, 'notes': notes}
    dwell = []
    for r in rows.values():
        try:
            dwell.append((dt.datetime.fromisoformat(r['submitted_at']) - dt.datetime.fromisoformat(r['opened_at'])).total_seconds())
        except (KeyError, ValueError):
            pass
    res['median_seconds_per_item'] = statistics.median(dwell) if dwell else None
    (out / 'scored.json').write_text(json.dumps(res, indent=1, ensure_ascii=False), encoding='utf-8')
    comp = json.loads((out / 'build.json').read_text(encoding='utf-8'))
    lines = report_lines(comp, res)
    if report:
        report.write_text('\n'.join(lines) + '\n', encoding='utf-8', newline='\n')
    print('\n'.join(results_lines(res)))
    return res


def results_lines(res: dict) -> list[str]:
    def pct(x):
        return '-' if x is None else f'{x:.3f}'
    k = res['kinds']
    L = ['## What came back', '',
         f"Readings returned {res['returned_readings']} of {res['assigned_readings']} assigned, from "
         f"{len(res['files'])} experts; median seconds per item "
         f"{'-' if res['median_seconds_per_item'] is None else round(res['median_seconds_per_item'])}.", '']
    a = k['answer']
    L += ['### B1, the answer check against the experts', '', '| | |', '|---|---|',
          f"| readings | {a['readings']} |", f"| three-way agreement with the check | {pct(a['agreement'])} |"]
    for ck, v in sorted(a['by_check'].items()):
        L.append(f"| the check said {ck}: experts said | " + ', '.join(f'{x} {n}' for x, n in sorted(v.items())) + ' |')
    L += [f"| by branch (same / different) | " + '; '.join(f"{b.replace('_engineering', '')} {v.get('same', 0)} / {v.get('different', 0)}"
                                                       for b, v in sorted(a['by_branch'].items())) + ' |',
          f"| disagreements by template | " + (', '.join(f'`{t}` {n}' for t, n in a['disagreements_by_template'].items()) or 'none') + ' |',
          f"| items read by two experts; their agreement; Cohen's kappa | {a['items_read_twice']}; {pct(a['reader_agreement'])}; {pct(a['reader_kappa'])} |",
          f"| notes left (local) | {len(a['notes'])} |", '']
    m = k['milestone']
    L += ['### B3, the judge against the experts', '', '| | |', '|---|---|', f"| readings | {m['readings']} |"]
    for j, v in sorted(m['by_judge'].items()):
        L.append(f"| the judge said {j}: experts said | " + ', '.join(f'"{x}" {n}' for x, n in sorted(v.items())) + ' |')
    L += [f"| REACHED confirmed (precision) | {pct(m['reached_precision'])} |",
          f"| MISSING confirmed (not obtained, by either answer) | {pct(m['missing_precision'])} |",
          f"| of MISSING, \"route does not need it\" | {pct(m['missing_not_needed_share'])} |",
          f"| items read by two experts; their agreement; Cohen's kappa | {m['items_read_twice']}; {pct(m['reader_agreement'])}; {pct(m['reader_kappa'])} |",
          f"| notes left (local) | {len(m['notes'])} |", '']
    e = k['error']
    L += ['### B2, the error categories', '', f"Readings {e['readings']}; Fleiss' kappa over the three readers, all models: {pct(e['fleiss_all'])}.", '',
          '| model | readings by category | item majorities | "No error" share | Fleiss |', '|---|---|---|---:|---:|']
    for mk, v in e['by_model'].items():
        L.append(f"| `{mk}` | " + ', '.join(f'{c.split(":")[0]} {n}' for c, n in sorted(v.items(), key=lambda kv: -kv[1])) + ' | '
                 + ', '.join(f'{c.split(":")[0]} {n}' for c, n in sorted(e['majority_by_model'].get(mk, {}).items(), key=lambda kv: -kv[1]))
                 + f" | {pct(e['no_error_share_by_model'].get(mk))} | {pct(e['fleiss_by_model'].get(mk))} |")
    L += ['', '| level | readings by category |', '|---|---|']
    for lv, v in e['by_level'].items():
        L.append(f"| {lv} | " + ', '.join(f'{c.split(":")[0]} {n}' for c, n in sorted(v.items(), key=lambda kv: -kv[1])) + ' |')
    L.append('')
    t = k['template']
    L += ['### B4, the templates', '', '| template | role | readers | question | answers |', '|---|---|---:|---|---|']
    for tid, v in t['templates'].items():
        for q, c in v['asks'].items():
            L.append(f"| `{tid.replace('template_', '')}` | {v['role']} | {v['readers']} | {q} | "
                     + ', '.join(f'"{x}" {n}' for x, n in sorted(c.items(), key=lambda kv: -kv[1])) + ' |')
    L += ['', f"Notes left on the templates (local): {len(t['notes'])}.", '']
    return L


# ------------------------------------------------------------------ selftest

def selftest() -> int:
    tmp = Path(tempfile.mkdtemp(prefix='engtrace-expert-kits-'))
    try:
        comp = build(tmp, SEED, b1=40, b3=20, b2_per_model=5, b2_models=B2_MODELS[:2], report=None)
        print(f"built {comp['items']} items, {comp['readings']} readings in {tmp}")
        pool = json.loads((tmp / 'tasks' / 'pool.json').read_text(encoding='utf-8'))
        queues = json.loads((tmp / 'tasks' / 'assignment.json').read_text(encoding='utf-8'))
        kinds_present = {t['kind'] for t in pool.values()}
        assert kinds_present == {'template', 'answer', 'milestone', 'error'}, kinds_present
        for c, t in pool.items():
            for forbidden in ('model', 'check', 'judge', 'item_id'):
                assert forbidden not in t, (c, forbidden)
        from streamlit.testing.v1 import AppTest
        preferred = max(queues, key=lambda a: (len({pool[c]['kind'] for c in queues[a]['codes']}), -len(queues[a]['codes'])))
        seen = []
        for kind in ORDER:
            holders = [a for a in [preferred] + sorted(queues) if (tmp / 'dist' / a / FOLDER[kind]).exists()]
            assert holders, f'no expert received a {kind} folder'
            aid = holders[0]
            folder = tmp / 'dist' / aid / FOLDER[kind]
            assert (folder / 'guide.md').read_text(encoding='utf-8').count('\n## ') == 4, kind
            at = AppTest.from_file(str(folder / 'app.py'), default_timeout=180)
            at.run()
            assert not at.exception, (kind, at.exception)
            at.sidebar.selectbox[0].select(aid).run()
            assert not at.exception, (kind, at.exception)
            first = json.loads((folder / 'tasks' / 'assignment.json').read_text(encoding='utf-8'))[aid]['codes'][0]
            assert any(first in w.value for w in at.subheader), kind
            seen.append((kind, aid))
            if kind == 'answer':   # submit one item and check the return file lands in that folder
                at.radio[0].set_value('correct').run()
                at.button[0].click().run()
                assert not at.exception, at.exception
                returned = folder / f'{aid}.jsonl'
                assert returned.exists(), 'no return file written'
                row = json.loads(returned.read_text(encoding='utf-8').splitlines()[0])
                assert row['code'] == first and row['answers']['verdict'] == 'correct', row
        print('each task folder\'s app renders, and the final-answer one takes a submission: ' + ', '.join(f'{k} ({a})' for k, a in seen))
        res = score_returns(tmp / 'dist', tmp, None)      # the returns are read from the folder tree
        assert res['returned_readings'] == 1, res['returned_readings']
        print('SELFTEST OK')
        return 0
    finally:
        shutil.rmtree(tmp, ignore_errors=True)


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument('--build', action='store_true')
    ap.add_argument('--score', metavar='DIR', help='the folder holding the returned <id>.jsonl files')
    ap.add_argument('--selftest', action='store_true')
    ap.add_argument('--seed', type=int, default=SEED)
    ap.add_argument('--b1', type=int, default=150, help='final-answer items (half wrong-side, half correct)')
    ap.add_argument('--b3', type=int, default=100, help='milestone items (half REACHED, half MISSING)')
    ap.add_argument('--b2-per-model', type=int, default=40)
    ap.add_argument('--b2-models', nargs='*', default=B2_MODELS)
    a = ap.parse_args()
    if a.selftest:
        return selftest()
    if a.build:
        comp = build(OUT, a.seed, a.b1, a.b3, a.b2_per_model, a.b2_models, REPORT)
        print('\n'.join(report_lines(comp, None)))
        print(f'\nkits: {OUT / "dist"}  (one folder per expert, a sub-folder per task; send the expert their folder)')
        return 0
    if a.score:
        score_returns(Path(a.score), OUT, REPORT)
        return 0
    ap.print_help()
    return 1


if __name__ == '__main__':
    sys.exit(main())
