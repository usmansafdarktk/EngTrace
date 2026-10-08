"""The experts' reading request for the full run: B1 to B4 of docs/EVALUATION_NEXT_STEPS.md, one kit per expert (D-172).

    python -m full_run_28092026.expert_kits --build [--seed 20261002] [--b1 150] [--b3 100] [--b2-per-model 40] [--b2-models K ...]
                                                                   # FREE: expert_request/tasks/ and expert_request/dist/, local
    python -m full_run_28092026.expert_kits --score DIR            # FREE: the returned <id>.jsonl files; EXPERT_REQUEST.md, counts only
    python -m full_run_28092026.expert_kits --score DIR --topup-returns DIR2     # ... merged with the B2 top-up's returns
    python -m full_run_28092026.expert_kits --b2-topup [--seed 20261007]          # FREE: expert_request/topup/, local
    python -m full_run_28092026.expert_kits --b2-matched VARIANT [--b2-models K ...]   # FREE: expert_request/matched_<VARIANT>/
    python -m full_run_28092026.expert_kits --score-matched VARIANT DIR           # FREE: its returns, a separate B2 section
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

WHEN THE STORE MOVES UNDER THE READINGS (WS-E, after two templates were repaired and re-run). A reading counts only
while the store still holds what was read: --score drops a B1, B2 or B3 reading whose item's question, or whose
trace, differs from the text the kit showed (the item was re-drawn, or the response re-run), and a B2 reading whose
answer the current store no longer scores incorrect (for instance after a change to the answer check). B1 is scored
against the check's CURRENT verdict, and the count of verdicts that changed since the reading is reported. B4 is
kept as read. Every exclusion is counted by kind and reason, and the denominators are stated.

  --b2-topup   For each B2 model, the readings dropped by that rule are replaced from the current main store's wrong
               answers not yet read, so that the sample is again --b2-per-model per model, allocated over the levels
               by the rule above (each level's shortfall against the model's proportional target first); a model
               with fewer wrong answers than that in all takes all of them. Three readers per item, the same app and
               guide section. Written to expert_request/topup/ (tasks, dist, build.json); --score merges its returns
               (--topup-returns) with the kept readings, and until they come back reports the reduced sample.
  --b2-matched The B2 draw (the same rule, a seed of its own) from another store's rows, for the models given:
               the error analysis of a re-run configuration. Written to expert_request/matched_<variant>/; scored by
               --score-matched into a section of its own, so the default-setting readings stay as they are.
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
RETURNED = OUT / 'expert_requests_filled'      # the returns of the 2026-10-02 build, as the owner filed them
TOPUP = OUT / 'topup'
TOPUP_SEED = 20261007
KEPT = 'kept'

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


def load(models: list[str], variant: str = 'main'):
    """Items, milestones, store rows, trace texts and E5 rows for the models given, from one store."""
    items = score.pool_items()
    ms_all = score.milestone_sets(items)
    rows, texts, e5 = {}, {}, {}
    for key in models:
        rows[key] = {r['item_id']: r for r in map(json.loads, (score.SCORES / variant / f'{key}.jsonl')
                                                   .read_text(encoding='utf-8').splitlines())}
        texts[key] = score.texts_matching(variant, key, rows[key].values())
        p = score.SCORES / variant / 'e5' / f'{key}.jsonl'
        e5[key] = {r['item_id']: r for r in map(json.loads, p.read_text(encoding='utf-8').splitlines())} \
            if p.exists() else {}
    return items, ms_all, rows, texts, e5


def code_of(*parts: str) -> str:
    return 'ER-' + hashlib.sha256('|'.join((SALT,) + parts).encode('utf-8')).hexdigest()[:8]


def answered(r: dict) -> bool:
    return r['status'] == 'answered' and not r['unusable']


def still_wrong(r: dict | None) -> bool:
    """An answered, usable response the store scores incorrect: what B2 reads."""
    return r is not None and answered(r) and r['answer']['label'] == 'incorrect'


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


def error_task(items, texts, key: str, i: str, lv: str, code_kind: str = 'error') -> dict:
    it = items[i]
    return {'pool': {'kind': 'error', 'branch': it['branch'], 'question': it['question'],
                     'solution': it['solution'], 'trace': texts[key][i], 'level': lv},
            'key': {'kind': 'error', 'item_id': i, 'model': key, 'template_id': it['template_id'], 'level': lv},
            'code': code_of(code_kind, key, i)}


def draw_b2(items, rows, texts, models: list[str], per_model: int, rng: random.Random,
            code_kind: str = 'error') -> list[dict]:
    out = []
    for key in models:
        cands = {lv: collections.defaultdict(list) for lv in LEVELS}
        for i, r in rows[key].items():
            if still_wrong(r):
                cands[r['level']][r['template_id']].append(i)
        alloc = allocate({lv: sum(map(len, cands[lv].values())) for lv in LEVELS}, per_model)
        for lv in LEVELS:
            for i in round_robin(cands[lv], alloc.get(lv, 0), rng):
                out.append(error_task(items, texts, key, i, lv, code_kind))
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


# ------------------------------------------------------------------ when the store moves: exclusions, top-up, matched

def shown_pool(out: Path) -> dict:
    """What each code of a build showed: tasks/pool.json, or the experts' task folders where that file is gone."""
    p = out / 'tasks' / 'pool.json'
    if p.exists():
        return json.loads(p.read_text(encoding='utf-8'))
    pool = {}
    for f in sorted((out / 'dist').glob('*/*/tasks/pool.json')):
        pool.update(json.loads(f.read_text(encoding='utf-8')))
    if not pool:
        raise SystemExit(f'{out}: neither tasks/pool.json nor dist/*/*/tasks/pool.json: the texts shown are needed')
    return pool


def reading_status(key: dict, shown: dict | None, items: dict, store: dict) -> str:
    """Does a reading still apply to the store? KEPT, or why not. A B4 reading is about a template and stays."""
    if key['kind'] == 'template':
        return KEPT
    if shown is None:
        raise SystemExit(f"no shown text for an item of {key['model']} on {key['item_id']}")
    it = items.get(key['item_id'])
    if it is None or it['question'] != shown['question']:
        return 'question changed'
    r = store.get(key['model'], {}).get(key['item_id'])
    if r is None or score.sha(shown['trace'].encode('utf-8')) != r.get('trace_sha256'):
        return 'trace replaced'
    if key['kind'] == 'answer' and not answered(r):
        return 'no longer answered'
    if key['kind'] == 'error' and not still_wrong(r):
        return 'no longer incorrect'
    return KEPT


def store_rows(models, variant: str = 'main') -> dict:
    return {m: {r['item_id']: r for r in map(json.loads, (score.SCORES / variant / f'{m}.jsonl')
                                              .read_text(encoding='utf-8').splitlines())} for m in models}


def level_need(target: dict, have: dict, avail: dict, n: int) -> dict:
    """n new items over the levels: each level's shortfall against the target first (shared in proportion when the
    shortfalls exceed n), then any remainder where items are left; never more than a level has."""
    short = {lv: min(max(0, target.get(lv, 0) - have.get(lv, 0)), avail.get(lv, 0)) for lv in LEVELS}
    need = allocate({lv: s for lv, s in short.items() if s}, n) if sum(short.values()) > n else dict(short)
    need = {lv: need.get(lv, 0) for lv in LEVELS}
    rest = n - sum(need.values())
    while rest > 0 and any(avail.get(lv, 0) > need[lv] for lv in LEVELS):
        lv = max(LEVELS, key=lambda x: (avail.get(x, 0) - need[x], -LEVELS.index(x)))
        need[lv] += 1
        rest -= 1
    return need


def topup_folders(base: Path) -> list[Path]:
    """The B2 top-up rounds built so far, in order: topup/, topup2/, ... A round is never rebuilt in place (its kit may
    be out with the experts, and their returns are scored against its keyfile): the next build is the next round."""
    out = []
    while True:
        f = base / ('topup' if not out else f'topup{len(out) + 1}')
        if not (f / 'build.json').exists():
            return out
        out.append(f)


def draw_topup(base: Path, models: list[str], per_model: int, seed: int,
               code_kind: str = 'error-topup') -> tuple[dict, list[dict]]:
    """Per B2 model: which readings still apply (the original B2 and every earlier top-up round), and the replacements
    that bring the sample back to per_model."""
    keyfile = json.loads((base / 'tasks' / 'keyfile.json').read_text(encoding='utf-8'))
    shown = shown_pool(base)
    for f in topup_folders(base):
        keyfile.update(json.loads((f / 'tasks' / 'keyfile.json').read_text(encoding='utf-8')))
        shown.update(shown_pool(f))
    items, _ms, rows, texts, _e5 = load(models)
    rng = random.Random(seed)
    plan, tasks = {}, []
    for key in models:
        read = {c: k for c, k in sorted(keyfile.items()) if k['kind'] == 'error' and k['model'] == key}
        status = {c: reading_status(k, shown.get(c), items, rows) for c, k in read.items()}
        kept = [c for c, s in status.items() if s == KEPT]
        kept_items = {read[c]['item_id'] for c in kept}
        wrong = {i: r for i, r in rows[key].items() if still_wrong(r)}
        cands = {lv: collections.defaultdict(list) for lv in LEVELS}
        for i, r in sorted(wrong.items()):
            if i not in kept_items:
                cands[r['level']][r['template_id']].append(i)
        avail = {lv: sum(map(len, cands[lv].values())) for lv in LEVELS}
        have = collections.Counter(rows[key][read[c]['item_id']]['level'] for c in kept)
        if len(wrong) <= per_model:
            need = dict(avail)                                   # fewer wrong answers than the sample: all of them
        else:
            target = allocate(dict(collections.Counter(r['level'] for r in wrong.values())), per_model)
            need = level_need(target, have, avail, max(0, per_model - len(kept)))
        added = collections.Counter()
        for lv in LEVELS:
            for i in round_robin(cands[lv], need.get(lv, 0), rng):
                tasks.append(error_task(items, texts, key, i, lv, code_kind))
                added[lv] += 1
        dropped = {c: s for c, s in status.items() if s != KEPT}
        plan[key] = {'read': len(read), 'kept': len(kept), 'kept_by_level': {lv: have[lv] for lv in LEVELS},
                     'dropped': dict(collections.Counter(dropped.values())),
                     'dropped_templates': dict(collections.Counter(read[c]['template_id'] for c in dropped)),
                     'wrong_in_store': len(wrong),
                     'wrong_by_level': dict(collections.Counter(r['level'] for r in wrong.values())),
                     'added': sum(added.values()), 'added_by_level': {lv: added[lv] for lv in LEVELS},
                     'final': len(kept) + sum(added.values())}
    return plan, tasks


def write_kits(out: Path, tasks: list[dict], queues: dict) -> None:
    """tasks/ (pool, assignment, keyfile) and dist/<id>/<task folder>/ for a list of tasks and their queues."""
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
        if listed:
            (dist / aid / 'README.txt').write_text(README_EXPERT.format(folders='\n'.join(listed)), encoding='utf-8')


def error_composition(tasks: list[dict], queues: dict) -> dict:
    e = [t for t in tasks if t['key']['kind'] == 'error']
    models = sorted({t['key']['model'] for t in e})
    return {'items': len(e), 'readers': READERS['error'],
            'by_model_level': {k: {lv: sum(t['key']['model'] == k and t['key']['level'] == lv for t in e)
                                   for lv in LEVELS} for k in models},
            'by_branch': dict(collections.Counter(t['pool']['branch'] for t in e)),
            'templates': len({t['key']['template_id'] for t in e}),
            'experts': {aid: len(q['codes']) for aid, q in sorted(queues.items())}}


def build_topup(base: Path, seed: int, per_model: int, models: list[str], report: Path | None) -> dict:
    """The next top-up round (topup/ first, then topup2/, ...); an earlier round is never overwritten."""
    rnd = len(topup_folders(base)) + 1
    top = base / ('topup' if rnd == 1 else f'topup{rnd}')
    if top.exists():
        raise SystemExit(f'{top} exists without a build.json: remove it deliberately before building round {rnd}')
    plan, tasks = draw_topup(base, models, per_model, seed + rnd - 1,
                             'error-topup' if rnd == 1 else f'error-topup{rnd}')
    if len({t['code'] for t in tasks}) != len(tasks):
        raise SystemExit('code collision')
    queues = assign(tasks, experts(), seed)
    write_kits(top, tasks, queues)
    comp = {'built': dt.date.today().isoformat(), 'round': rnd, 'seed': seed + rnd - 1, 'store_commit': store_commit(),
            'per_model': per_model, 'models': plan, **error_composition(tasks, queues),
            'readings': sum(len(q['codes']) for q in queues.values())}
    (top / 'build.json').write_text(json.dumps(comp, indent=1), encoding='utf-8')
    if report:
        write_report(base, report)
    return comp


def build_matched(base: Path, variant: str, models: list[str], per_model: int, seed: int,
                  report: Path | None) -> dict:
    """B2 from another store (a re-run configuration), the same draw rule; kits in matched_<variant>/."""
    items, _ms, rows, texts, _e5 = load(models, variant)
    tasks = draw_b2(items, rows, texts, models, per_model, random.Random(seed), code_kind=f'error-{variant}')
    if len({t['code'] for t in tasks}) != len(tasks):
        raise SystemExit('code collision')
    folder = base / f'matched_{variant}'
    queues = assign(tasks, experts(), seed)
    write_kits(folder, tasks, queues)
    try:
        commit = json.loads((score.SCORES / variant / 'CONFIG.json').read_text(encoding='utf-8')).get('git', '?')[:7]
    except Exception:  # noqa: BLE001
        commit = '?'
    comp = {'built': dt.date.today().isoformat(), 'seed': seed, 'variant': variant, 'store_commit': commit,
            'per_model': per_model, 'wrong_in_store': {k: sum(still_wrong(r) for r in rows[k].values()) for k in models},
            **error_composition(tasks, queues), 'readings': sum(len(q['codes']) for q in queues.values())}
    (folder / 'build.json').write_text(json.dumps(comp, indent=1), encoding='utf-8')
    if report:
        write_report(base, report)
    return comp


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
    write_kits(out, tasks, queues)
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


def report_lines(comp: dict, results: dict | None, topup: dict | None = None, matched: list | None = None,
                 cur: dict | None = None) -> list[str]:
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
    if topup:
        L += topup_lines(topup)
    if cur:
        L += current_lines(cur)
    for mc, mr in matched or ():
        L += matched_lines(mc, mr)
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


def b1_stats(keyfile: dict, by_code: dict, branch_of: dict, store: dict | None = None) -> dict:
    """B1: the experts against the check's verdict at the reading (store None), or against its current verdict, the
    store's, with the count of verdicts that changed since the reading."""
    agree = collections.Counter()
    by_check = collections.defaultdict(collections.Counter)
    by_branch = collections.defaultdict(collections.Counter)
    dis_t = collections.Counter()
    pairs, notes, changed, items = [], [], 0, 0
    for c, readers in by_code.items():
        k = keyfile[c]
        if k['kind'] != 'answer' or not readers:
            continue
        items += 1
        check = k['check'] if store is None else store[k['model']][k['item_id']]['answer']['label']
        changed += check != k['check']
        labs = {}
        for aid, r in readers.items():
            v = r['answers']['verdict']
            labs[aid] = v
            e = EXPERT_TO_CHECK[v]
            by_check[check][e] += 1
            same = e == check
            agree['same' if same else 'different'] += 1
            by_branch[branch_of[aid]]['same' if same else 'different'] += 1
            if not same:
                dis_t[k['template_id']] += 1
            if r['note']:
                notes.append({'code': c, 'reader': aid, 'check': check, 'expert': v, 'template': k['template_id'],
                              'model': k['model'], 'note': r['note']})
        if len(labs) == 2:
            pairs.append(tuple(labs[a] for a in sorted(labs)))
    out = {'readings': sum(agree.values()), 'agreement': agree['same'] / sum(agree.values()) if agree else None,
           'by_check': {ck: dict(v) for ck, v in by_check.items()},
           'by_branch': {b: dict(v) for b, v in by_branch.items()},
           'disagreements_by_template': dict(dis_t.most_common()),
           'items_read_twice': len(pairs), 'reader_agreement': (sum(a == b for a, b in pairs) / len(pairs)) if pairs else None,
           'reader_kappa': cohen(pairs), 'notes': notes}
    if store is not None:
        out.update(items=items, verdicts_changed_since_reading=changed)
    return out


def b3_stats(keyfile: dict, by_code: dict) -> dict:
    by_judge = collections.defaultdict(collections.Counter)
    pairs, notes, items = [], [], 0
    for c, readers in by_code.items():
        k = keyfile[c]
        if k['kind'] != 'milestone' or not readers:
            continue
        items += 1
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
    rj, mj = by_judge['REACHED'], by_judge['MISSING']
    return {'readings': sum(sum(v.values()) for v in by_judge.values()),
            'by_judge': {j: dict(v) for j, v in by_judge.items()},
            'reached_precision': rj[yes] / sum(rj.values()) if rj else None,
            'missing_precision': (mj[no1] + mj[no2]) / sum(mj.values()) if mj else None,
            'missing_not_needed_share': mj[no2] / sum(mj.values()) if mj else None,
            'items_read_twice': len(pairs), 'reader_agreement': (sum(a == b for a, b in pairs) / len(pairs)) if pairs else None,
            'reader_kappa': cohen(pairs), 'notes': notes, 'items': items}


def b2_stats(keyfile: dict, by_code: dict, level_of: dict | None = None) -> dict:
    """The error categories over the readings given; `level_of` (code -> level) overrides the keyfile's level."""
    per_model = collections.defaultdict(collections.Counter)
    per_level = collections.defaultdict(collections.Counter)
    majority = collections.defaultdict(collections.Counter)
    item_rows = collections.defaultdict(list)
    notes = []
    for c, readers in by_code.items():
        k = keyfile[c]
        if k['kind'] != 'error' or not readers:
            continue
        lv = (level_of or {}).get(c, k['level'])
        cnt = collections.Counter(r['answers']['category'] for r in readers.values())
        for cat, n in cnt.items():
            per_model[k['model']][cat] += n
            per_level[lv][cat] += n
        top, n = cnt.most_common(1)[0]
        majority[k['model']][top if n >= 2 else 'split'] += 1
        item_rows[k['model']].append(cnt)
        for aid, r in readers.items():
            if r['note'] or r['excerpt']:
                notes.append({'code': c, 'reader': aid, 'model': k['model'], 'template': k['template_id'],
                              'category': r['answers']['category'], 'excerpt': r['excerpt'], 'note': r['note']})
    return {'readings': sum(sum(v.values()) for v in per_model.values()),
            'by_model': {mk: dict(v) for mk, v in per_model.items()},
            'by_level': {lv: dict(v) for lv, v in per_level.items()},
            'majority_by_model': {mk: dict(v) for mk, v in majority.items()},
            'no_error_share_by_model': {mk: v[ERROR_OPTIONS[6]] / sum(v.values()) for mk, v in per_model.items() if sum(v.values())},
            'fleiss_by_model': {mk: fleiss(v) for mk, v in item_rows.items()},
            'fleiss_all': fleiss([r for v in item_rows.values() for r in v]), 'notes': notes,
            'items_read': sum(len(v) for v in item_rows.values())}


def b4_stats(keyfile: dict, by_code: dict) -> dict:
    t4, notes = {}, []
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
    return {'templates': t4, 'notes': notes}


def median_dwell(rows: dict) -> float | None:
    dwell = []
    for r in rows.values():
        try:
            dwell.append((dt.datetime.fromisoformat(r['submitted_at']) - dt.datetime.fromisoformat(r['opened_at'])).total_seconds())
        except (KeyError, ValueError):
            pass
    return statistics.median(dwell) if dwell else None


def as_returned(keyfile: dict, queues: dict, rows: dict) -> dict:
    """The figures as the readings came back, every reading counted, B1 against the verdict the reader saw judged:
    what --score wrote before the store moved (scored.json), and what the report's "What came back" prints."""
    by_code = collections.defaultdict(dict)
    for (aid, c), r in rows.items():
        by_code[c][aid] = r
    branch_of = {aid: q['branch'] for aid, q in queues.items()}
    return {'returned_readings': len(rows), 'assigned_readings': sum(len(q['codes']) for q in queues.values()),
            'files': sorted({aid for aid, _ in rows}),
            'kinds': {'answer': b1_stats(keyfile, by_code, branch_of), 'milestone': b3_stats(keyfile, by_code),
                      'error': b2_stats(keyfile, by_code), 'template': b4_stats(keyfile, by_code)},
            'median_seconds_per_item': median_dwell(rows)}


def current(keyfile: dict, queues: dict, rows: dict, out: Path, topup_returned: Path | None) -> dict:
    """The figures over the readings the store still holds, merged with the B2 top-up where one was built."""
    shown = shown_pool(out)
    items = score.pool_items()
    models = {k['model'] for k in keyfile.values() if 'model' in k}
    rounds = []                               # each top-up round: its keyfile, queues, returns and shown texts
    for n, top in enumerate(topup_folders(out), 1):
        rkey = json.loads((top / 'tasks' / 'keyfile.json').read_text(encoding='utf-8'))
        rqueues = json.loads((top / 'tasks' / 'assignment.json').read_text(encoding='utf-8'))
        found = [top / 'returned'] + sorted(top.glob('experts_filled*'))        # where the owner files a round's returns
        rdir = topup_returned if (n == 1 and topup_returned is not None) else next((d for d in found if d.exists()), found[0])
        rrows = read_returns(rdir, top)[2] if rdir.exists() else {}
        rounds.append((rkey, rqueues, rrows, shown_pool(top)))
        models |= {k['model'] for k in rkey.values()}
    tkey = {c: k for r in rounds for c, k in r[0].items()}
    tqueues = {}
    for _k, rq, _r, _s in rounds:
        for aid, q in rq.items():
            tqueues.setdefault(aid, {'branch': q['branch'], 'codes': []})['codes'] += q['codes']
    trows = {key: row for r in rounds for key, row in r[2].items()}
    store = store_rows(sorted(models))
    status = {c: reading_status(k, shown.get(c), items, store) for c, k in keyfile.items()}
    topup = None
    if rounds:
        tstatus = {c: reading_status(k, rs.get(c), items, store) for rk, _q, _r, rs in rounds for c, k in rk.items()}
        topup = {'rounds': len(rounds), 'items': len(tkey), 'assigned': sum(len(q['codes']) for q in tqueues.values()),
                 'returned': len(trows), 'files': sorted({aid for aid, _ in trows}),
                 'no_longer_apply': sum(s != KEPT for s in tstatus.values())}
        status.update(tstatus)
    keyfile = {**keyfile, **tkey}
    rows = {**rows, **trows}
    branch_of = {aid: q['branch'] for aid, q in {**queues, **tqueues}.items()}
    excluded_items = collections.defaultdict(collections.Counter)
    excluded_templates = collections.defaultdict(collections.Counter)
    for c, s in status.items():
        if s != KEPT:
            excluded_items[keyfile[c]['kind']][s] += 1
            excluded_templates[keyfile[c]['kind']][keyfile[c]['template_id']] += 1
    excluded_readings = collections.defaultdict(collections.Counter)
    by_code = collections.defaultdict(dict)
    for (aid, c), r in rows.items():
        if status[c] == KEPT:
            by_code[c][aid] = r
        else:
            excluded_readings[keyfile[c]['kind']][status[c]] += 1
    for c, s in status.items():              # a B2 item in the sample whose readings are not back yet
        if s == KEPT and keyfile[c]['kind'] == 'error':
            by_code.setdefault(c, {})
    sample = [c for c in by_code if keyfile[c]['kind'] == 'error']
    level_of = {c: store[keyfile[c]['model']][keyfile[c]['item_id']]['level'] for c in sample}
    err = b2_stats(keyfile, by_code, level_of)
    err['composition'] = {
        'items': len(sample), 'items_read': err['items_read'], 'readers': READERS['error'],
        'by_model_level': {m: {lv: sum(keyfile[c]['model'] == m and level_of[c] == lv for c in sample) for lv in LEVELS}
                           for m in sorted({keyfile[c]['model'] for c in sample})},
        'by_branch': dict(collections.Counter(items[keyfile[c]['item_id']]['branch'] for c in sample)),
        'templates': len({keyfile[c]['template_id'] for c in sample}),
        'from_topup': sum(c in tkey for c in sample)}
    return {'basis': 'the readings the current main store still holds (question and trace as shown; a B2 answer '
                     'still scored incorrect), B1 against the check\'s current verdict, merged with the B2 top-up',
            'returned_readings': len(rows),
            'assigned_readings': sum(len(q['codes']) for q in queues.values()) + sum(len(q['codes']) for q in tqueues.values()),
            'files': sorted({aid for aid, _ in rows}),
            'exclusions': {'items': {k: dict(v) for k, v in excluded_items.items()},
                           'templates': {k: dict(v) for k, v in excluded_templates.items()},
                           'readings': {k: dict(v) for k, v in excluded_readings.items()}},
            'topup': topup,
            'kinds': {'answer': b1_stats(keyfile, by_code, branch_of, store), 'milestone': b3_stats(keyfile, by_code),
                      'error': err, 'template': b4_stats(keyfile, by_code)},
            'median_seconds_per_item': median_dwell(rows)}


SCORED_KEYS = {'answer': ('readings', 'agreement', 'by_check', 'by_branch', 'disagreements_by_template',
                          'items_read_twice', 'reader_agreement', 'reader_kappa', 'notes'),
               'milestone': ('readings', 'by_judge', 'reached_precision', 'missing_precision', 'missing_not_needed_share',
                             'items_read_twice', 'reader_agreement', 'reader_kappa', 'notes'),
               'error': ('readings', 'by_model', 'by_level', 'majority_by_model', 'no_error_share_by_model',
                         'fleiss_by_model', 'fleiss_all', 'notes'),
               'template': ('templates', 'notes')}


def same_as_stored(res: dict, stored: dict) -> list[str]:
    """Where a recomputed as-returned result differs from scored.json as first written (the keys it had)."""
    bad = [k for k in ('returned_readings', 'assigned_readings', 'files', 'median_seconds_per_item')
           if json.loads(json.dumps(res[k])) != stored.get(k)]
    for kind, keys in SCORED_KEYS.items():
        for k in keys:
            if json.loads(json.dumps(res['kinds'][kind][k])) != stored['kinds'][kind].get(k):
                bad.append(f'{kind}.{k}')
    return bad


def score_returns(returned: Path, out: Path, report: Path | None, topup_returned: Path | None = None) -> dict:
    """scored.json (as returned) and scored_current.json (what the store still holds, with the top-up merged)."""
    keyfile, queues, rows = read_returns(returned, out)
    res = as_returned(keyfile, queues, rows)
    scored = out / 'scored.json'
    bad = same_as_stored(res, json.loads(scored.read_text(encoding='utf-8'))) if scored.exists() else ['new']
    if bad:                                   # left as it is when the same returns give the same figures
        if scored.exists():
            print(f'NOTE: the returns give figures other than scored.json for: {", ".join(bad)} (rewritten)')
        scored.write_text(json.dumps(res, indent=1, ensure_ascii=False), encoding='utf-8')
    cur = current(keyfile, queues, rows, out, topup_returned)
    (out / 'scored_current.json').write_text(json.dumps(cur, indent=1, ensure_ascii=False), encoding='utf-8')
    if report:
        write_report(out, report)
    print('\n'.join(results_lines(res) + current_lines(cur)))
    return {'as_returned': res, 'current': cur}


def score_matched(base: Path, variant: str, returned: Path, report: Path | None) -> dict:
    folder = base / f'matched_{variant}'
    keyfile, queues, rows = read_returns(returned, folder)
    by_code = collections.defaultdict(dict)
    for (aid, c), r in rows.items():
        by_code[c][aid] = r
    res = {'variant': variant, 'returned_readings': len(rows),
           'assigned_readings': sum(len(q['codes']) for q in queues.values()),
           'files': sorted({aid for aid, _ in rows}), 'error': b2_stats(keyfile, by_code),
           'median_seconds_per_item': median_dwell(rows)}
    (folder / 'scored.json').write_text(json.dumps(res, indent=1, ensure_ascii=False), encoding='utf-8')
    if report:
        write_report(base, report)
    print('\n'.join(b2_lines(res['error'], f'### B2, matched configuration (`{variant}`)')))
    return res


def write_report(base: Path, report: Path) -> None:
    """EXPERT_REQUEST.md from what is on disk: the build, the returns as scored, the top-up, the current figures and
    any matched B2. The first two sections are what --build and --score always wrote."""
    def opt(p: Path):
        return json.loads(p.read_text(encoding='utf-8')) if p.exists() else None
    comp = json.loads((base / 'build.json').read_text(encoding='utf-8'))
    matched = [(json.loads(b.read_text(encoding='utf-8')), opt(b.parent / 'scored.json'))
               for b in sorted(base.glob('matched_*/build.json'))]
    report.write_text('\n'.join(report_lines(comp, opt(base / 'scored.json'), opt(base / 'topup' / 'build.json'), matched,
                                             opt(base / 'scored_current.json'))) + '\n', encoding='utf-8', newline='\n')


def pct(x) -> str:
    return '-' if x is None else f'{x:.3f}'


def results_lines(res: dict) -> list[str]:
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
    L += ['### B2, the error categories', '', f"Readings {e['readings']}; Fleiss' kappa over the three readers, all models: {pct(e['fleiss_all'])}.", '']
    L += b2_table(e)
    t = k['template']
    L += ['### B4, the templates', '', '| template | role | readers | question | answers |', '|---|---|---:|---|---|']
    for tid, v in t['templates'].items():
        for q, c in v['asks'].items():
            L.append(f"| `{tid.replace('template_', '')}` | {v['role']} | {v['readers']} | {q} | "
                     + ', '.join(f'"{x}" {n}' for x, n in sorted(c.items(), key=lambda kv: -kv[1])) + ' |')
    L += ['', f"Notes left on the templates (local): {len(t['notes'])}.", '']
    return L


def b2_table(e: dict) -> list[str]:
    L = ['| model | readings by category | item majorities | "No error" share | Fleiss |', '|---|---|---|---:|---:|']
    for mk, v in e['by_model'].items():
        L.append(f"| `{mk}` | " + ', '.join(f'{c.split(":")[0]} {n}' for c, n in sorted(v.items(), key=lambda kv: -kv[1])) + ' | '
                 + ', '.join(f'{c.split(":")[0]} {n}' for c, n in sorted(e['majority_by_model'].get(mk, {}).items(), key=lambda kv: -kv[1]))
                 + f" | {pct(e['no_error_share_by_model'].get(mk))} | {pct(e['fleiss_by_model'].get(mk))} |")
    L += ['', '| level | readings by category |', '|---|---|']
    for lv, v in e['by_level'].items():
        L.append(f"| {lv} | " + ', '.join(f'{c.split(":")[0]} {n}' for c, n in sorted(v.items(), key=lambda kv: -kv[1])) + ' |')
    L.append('')
    return L


def b2_lines(e: dict, title: str) -> list[str]:
    return [title, '', f"Readings {e['readings']} on {e['items_read']} items; Fleiss' kappa over the three readers, all "
                       f"models: {pct(e['fleiss_all'])}.", ''] + b2_table(e)


def current_lines(cur: dict) -> list[str]:
    """The current figures. Row labels differ from "What came back" on purpose: a parser of that section must not
    pick these up unawares (docs/appendix_evaluation.py reads its rows by their wording)."""
    k = cur['kinds']
    ex = cur['exclusions']
    names = {'answer': 'B1', 'milestone': 'B3', 'error': 'B2'}
    parts = []
    for kind in ('answer', 'milestone', 'error'):
        its = ex['items'].get(kind, {})
        if its:
            rd = ex['readings'].get(kind, {})
            parts.append(f"{names[kind]} {sum(its.values())} items and {sum(rd.values())} readings ("
                         + ', '.join(f'{why} {n}' for why, n in sorted(its.items())) + '; '
                         + ', '.join(f"`{t.replace('template_', '')}` {n}" for t, n in
                                     sorted(ex['templates'][kind].items(), key=lambda kv: (-kv[1], kv[0]))) + ')')
    L = ['## The current figures: the readings the store still holds', '',
         'Generated by `expert_kits.py --score` into `scored_current.json`. A reading counts while the main store still '
         'holds what was read: the item\'s question and the trace as shown, and, for B2, an answer the check still scores '
         'incorrect. B1 is scored against the check\'s current verdict. The B2 sample is merged with the top-up\'s '
         'readings as they come back. ' + ('Left out: ' + '; '.join(parts) + '.' if parts else 'Nothing is left out.'), '']
    if cur.get('topup'):
        t = cur['topup']
        L += [f"B2 top-up: {t['items']} items sent, {t['assigned']} readings; {t['returned']} returned from "
              f"{len(t['files'])} experts.", '']
    a = k['answer']
    L += ['### B1 now', '', '| | |', '|---|---|', f"| items; readings | {a['items']}; {a['readings']} |",
          f"| three-way agreement with the check's current verdict | {pct(a['agreement'])} |"]
    for ck, v in sorted(a['by_check'].items()):
        L.append(f"| the check's current verdict {ck}: experts said | " + ', '.join(f'{x} {n}' for x, n in sorted(v.items())) + ' |')
    L += [f"| by branch, same / different | " + '; '.join(f"{b.replace('_engineering', '')} {v.get('same', 0)} / {v.get('different', 0)}"
                                                     for b, v in sorted(a['by_branch'].items())) + ' |',
          f"| verdicts changed since the reading | {a['verdicts_changed_since_reading']} |",
          f"| items with two readers; agreement; Cohen's kappa | {a['items_read_twice']}; {pct(a['reader_agreement'])}; {pct(a['reader_kappa'])} |", '']
    m = k['milestone']
    L += ['### B3 now', '', '| | |', '|---|---|', f"| items; readings | {m['items']}; {m['readings']} |"]
    for j, v in sorted(m['by_judge'].items()):
        L.append(f"| the judge ruled {j}: experts said | " + ', '.join(f'"{x}" {n}' for x, n in sorted(v.items())) + ' |')
    L += [f"| REACHED confirmed | {pct(m['reached_precision'])} |",
          f"| MISSING confirmed, by either answer | {pct(m['missing_precision'])} |",
          f"| of MISSING, \"route does not need it\" | {pct(m['missing_not_needed_share'])} |",
          f"| items with two readers; agreement; Cohen's kappa | {m['items_read_twice']}; {pct(m['reader_agreement'])}; {pct(m['reader_kappa'])} |", '']
    e = k['error']
    comp = e['composition']
    L += ['### B2 now', '',
          'The sample: ' + '; '.join(f"`{mk}` {sum(lv.values())} (" + ', '.join(f'{x} {n}' for x, n in lv.items() if n) + ')'
                                     for mk, lv in comp['by_model_level'].items())
          + f"; {comp['items']} items on {comp['templates']} templates"
          + (f" ({comp['from_topup']} from the top-up)" if comp['from_topup'] else '')
          + f"; {comp['items_read']} of them read so far.", '',
          f"Readings {e['readings']} on {e['items_read']} items; Fleiss' kappa over the three readers, all models: "
          f"{pct(e['fleiss_all'])}.", '']
    L += b2_table(e)
    return L


def topup_lines(t: dict) -> list[str]:
    L = ['## The B2 top-up: what was sent', '',
         f"Built {t['built']} from the main store at `{t['store_commit']}`, seed {t['seed']}: per model, the readings whose "
         'question or trace the store no longer holds, or whose answer it no longer scores incorrect, are replaced from '
         f"the wrong answers not yet read, so that each model again has {t['per_model']} (all of them where it has fewer), "
         'allocated over the levels by the original rule. Three readers per item.', '',
         '| model | read | dropped | kept | wrong answers in the store | added | sample now |',
         '|---|---:|---|---:|---:|---|---:|']
    for mk, v in t['models'].items():
        dropped = ', '.join(f'{why} {n}' for why, n in sorted(v['dropped'].items())) or '0'
        added = ', '.join(f'{lv} {n}' for lv, n in v['added_by_level'].items() if n) or '0'
        L.append(f"| `{mk}` | {v['read']} | {dropped} | {v['kept']} | {v['wrong_in_store']} | {added} | {v['final']} |")
    L += ['', f"Items {t['items']}, readings {t['readings']}; per expert: "
              + (', '.join(f'`{aid}` {n}' for aid, n in t['experts'].items()) or 'none') + '.', '']
    return L


def matched_lines(comp: dict, res: dict | None) -> list[str]:
    L = [f"## B2 for the matched configuration (`{comp['variant']}`)", '',
         f"Built {comp['built']} from the `{comp['variant']}` store at `{comp['store_commit']}`, seed {comp['seed']}, by "
         f"the B2 rule: {comp['per_model']} wrong answers per model, or all where fewer. Wrong answers in that store: "
         + ', '.join(f'`{m}` {n}' for m, n in comp['wrong_in_store'].items()) + '. Sent: '
         + '; '.join(f"`{m}` " + ', '.join(f'{lv} {n}' for lv, n in v.items() if n) for m, v in comp['by_model_level'].items())
         + f"; {comp['items']} items, {comp['readings']} readings.", '']
    if res:
        L += [f"Returned {res['returned_readings']} of {res['assigned_readings']} readings from {len(res['files'])} experts.", '']
        L += b2_lines(res['error'], '### B2, the error categories (matched configuration)')
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
        res = score_returns(tmp / 'dist', tmp, None)['as_returned']      # the returns are read from the folder tree
        assert res['returned_readings'] == 1, res['returned_readings']
        selftest_moves(tmp)
        print('SELFTEST OK')
        return 0
    finally:
        shutil.rmtree(tmp, ignore_errors=True)


def simulate_errors(folder: Path, out: Path, seed: int = 1) -> int:
    """Synthetic B2 returns for a kit folder: every assigned reader of every error item, a category at random."""
    keyfile = json.loads((folder / 'tasks' / 'keyfile.json').read_text(encoding='utf-8'))
    queues = json.loads((folder / 'tasks' / 'assignment.json').read_text(encoding='utf-8'))
    rng = random.Random(seed)
    out.mkdir(parents=True, exist_ok=True)
    n = 0
    for aid, q in sorted(queues.items()):
        with (out / f'{aid}.jsonl').open('w', encoding='utf-8') as fh:
            for pos, c in enumerate(q['codes'], 1):
                if keyfile[c]['kind'] != 'error':
                    continue
                cat = rng.choice(ERROR_OPTIONS[:7])
                fh.write(json.dumps({'annotator_id': aid, 'branch': q['branch'], 'kind': 'error', 'code': c,
                                     'position': pos, 'opened_at': '2026-10-08T09:00:00+00:00',
                                     'answers': {'category': cat}, 'excerpt': '' if cat in ERROR_OPTIONS[6:] else 'x',
                                     'note': 'n' if cat == ERROR_OPTIONS[6] else '',
                                     'submitted_at': '2026-10-08T09:01:00+00:00', 'source': 'simulated'}) + '\n')
                n += 1
    return n


def selftest_moves(tmp: Path) -> None:
    """The store moving under the readings, simulated by editing what the kit records it showed: two B2 traces
    replaced and one B1 question re-drawn. The top-up replaces the two, --score leaves the three out and merges."""
    pool_path = tmp / 'tasks' / 'pool.json'
    pool = json.loads(pool_path.read_text(encoding='utf-8'))
    keyfile = json.loads((tmp / 'tasks' / 'keyfile.json').read_text(encoding='utf-8'))
    errs = sorted(c for c, k in keyfile.items() if k['kind'] == 'error')
    victim = keyfile[errs[0]]['model']
    other = next(k['model'] for c, k in keyfile.items() if k['kind'] == 'error' and k['model'] != victim)
    hit = [c for c in errs if keyfile[c]['model'] == victim][:2]
    for c in hit:
        pool[c]['trace'] += ' (a later response)'
    b1 = sorted(c for c, k in keyfile.items() if k['kind'] == 'answer')[0]
    pool[b1]['question'] += ' (re-drawn)'
    pool_path.write_text(json.dumps(pool, ensure_ascii=False), encoding='utf-8')
    n_b2 = simulate_errors(tmp, tmp / 'returns')
    comp = build_topup(tmp, TOPUP_SEED, 5, sorted({victim, other}), None)
    v, o = comp['models'][victim], comp['models'][other]
    assert v['dropped'] == {'trace replaced': 2} and v['kept'] == 3 and o['dropped'] == {} and o['added'] == 0, comp
    assert v['final'] == min(5, v['wrong_in_store']) and v['added'] == v['final'] - 3, v
    added = json.loads((tmp / 'topup' / 'tasks' / 'keyfile.json').read_text(encoding='utf-8'))
    kept_items = {keyfile[c]['item_id'] for c in errs if keyfile[c]['model'] == victim and c not in hit}
    assert not kept_items & {k['item_id'] for k in added.values()}, 'a replacement repeats a kept reading'
    n_top = simulate_errors(tmp / 'topup', tmp / 'topup' / 'returned', seed=2)   # the round's own returns folder
    both = score_returns(tmp / 'returns', tmp, None)
    res = both['current']
    assert both['as_returned']['kinds']['error']['readings'] == n_b2 and 'exclusions' not in both['as_returned']
    assert not same_as_stored(both['as_returned'], json.loads((tmp / 'scored.json').read_text(encoding='utf-8')))
    ex = res['exclusions']
    assert ex['items'] == {'error': {'trace replaced': 2}, 'answer': {'question changed': 1}}, ex
    assert ex['readings']['error'] == {'trace replaced': 6}, ex
    e = res['kinds']['error']
    assert e['readings'] == n_b2 - 6 + n_top and e['composition']['items'] == 10 - 2 + v['added'], e['composition']
    assert e['composition']['from_topup'] == v['added'] and e['fleiss_all'] is not None
    assert res['kinds']['answer']['items'] == 0                     # no B1 reading was returned in this simulation
    assert res['assigned_readings'] == both['as_returned']['assigned_readings'] + comp['readings'], res['assigned_readings']
    again = build_topup(tmp, TOPUP_SEED, 5, sorted({victim, other}), None)      # a second build is round 2, not a rebuild
    assert again['round'] == 2 and (tmp / 'topup' / 'tasks' / 'keyfile.json').read_text(encoding='utf-8') == json.dumps(added, indent=1)
    assert again['models'][victim]['dropped'] == {'trace replaced': 2} and again['models'][victim]['added'] == 0, again['models'][victim]
    print(f"store moves: two B2 readings and one B1 item left out; the top-up adds {v['added']} for `{victim}`, none for "
          f"`{other}`; --score merges {n_top} top-up readings with the {n_b2 - 6} kept; a second build is round 2 and "
          'leaves round 1 as it was')
    variant = 'reasoning-medium'
    if (score.SCORES / variant / 'gpt-5.4-mini.jsonl').exists():
        mc = build_matched(tmp, variant, ['gpt-5.4-mini'], 3, TOPUP_SEED, None)
        n_m = simulate_errors(tmp / f'matched_{variant}', tmp / 'matched_returns', seed=3)
        mr = score_matched(tmp, variant, tmp / 'matched_returns', None)
        assert mc['items'] == 3 and mr['error']['readings'] == n_m == 9, (mc['items'], mr['error']['readings'])
        write_report(tmp, tmp / 'EXPERT_REQUEST.md')
        text = (tmp / 'EXPERT_REQUEST.md').read_text(encoding='utf-8')
        assert '## The B2 top-up: what was sent' in text and f'matched configuration (`{variant}`)' in text
        assert '## The current figures' in text and text.index('## What came back') < text.index('## The current figures')
        cur_part = text[text.index('## The current figures'):]
        assert '| the check said' not in cur_part and '| the judge said' not in cur_part   # not the old rows' wording
        print(f'matched configuration: a B2 kit from `{variant}` is built, scored into its own section, and the report '
              'carries the top-up and the matched sections')


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument('--build', action='store_true')
    ap.add_argument('--score', metavar='DIR', help='the folder holding the returned <id>.jsonl files')
    ap.add_argument('--topup-returns', metavar='DIR', help="with --score: the B2 top-up's returned files")
    ap.add_argument('--b2-topup', action='store_true', help='replace the B2 readings the store no longer holds')
    ap.add_argument('--b2-matched', metavar='VARIANT', help='a B2 kit from another store (a re-run configuration)')
    ap.add_argument('--score-matched', nargs=2, metavar=('VARIANT', 'DIR'))
    ap.add_argument('--selftest', action='store_true')
    ap.add_argument('--seed', type=int, default=None, help=f'default {SEED} for --build, {TOPUP_SEED} for the B2 additions')
    ap.add_argument('--b1', type=int, default=150, help='final-answer items (half wrong-side, half correct)')
    ap.add_argument('--b3', type=int, default=100, help='milestone items (half REACHED, half MISSING)')
    ap.add_argument('--b2-per-model', type=int, default=40)
    ap.add_argument('--b2-models', nargs='*', default=B2_MODELS)
    a = ap.parse_args()
    if a.selftest:
        return selftest()
    if a.b2_topup:
        comp = build_topup(OUT, a.seed or TOPUP_SEED, a.b2_per_model, a.b2_models, REPORT)
        print('\n'.join(topup_lines(comp)))
        print(f'kits: {TOPUP / "dist"}  (one folder per expert; send the expert their folder)')
        return 0
    if a.b2_matched:
        models = [m for m in a.b2_models if (score.SCORES / a.b2_matched / f'{m}.jsonl').exists()]
        comp = build_matched(OUT, a.b2_matched, models, a.b2_per_model, a.seed or TOPUP_SEED, REPORT)
        print('\n'.join(matched_lines(comp, None)))
        print(f'kits: {OUT / ("matched_" + a.b2_matched) / "dist"}')
        return 0
    if a.score_matched:
        score_matched(OUT, a.score_matched[0], Path(a.score_matched[1]), REPORT)
        return 0
    if a.build:
        comp = build(OUT, a.seed or SEED, a.b1, a.b3, a.b2_per_model, a.b2_models, REPORT)
        print('\n'.join(report_lines(comp, None)))
        print(f'\nkits: {OUT / "dist"}  (one folder per expert, a sub-folder per task; send the expert their folder)')
        return 0
    if a.score:
        score_returns(Path(a.score), OUT, REPORT, Path(a.topup_returns) if a.topup_returns else None)
        return 0
    ap.print_help()
    return 1


if __name__ == '__main__':
    sys.exit(main())
