"""The simplified annotator kit: one HTML file per expert, nothing to install.

    python -m template_annotation_23092026.layer2.make_kits            # every expert
    python -m template_annotation_23092026.layer2.make_kits --id civ-2

Writes dist/kit_<id>/ with

    <id>.html                          the whole review in one page: their 34 items embedded,
                                       hand check -> solution -> verdict, progress kept in the
                                       browser, "Download my answers" -> <id>.jsonl
    app/app.py + app/tasks/            the Streamlit app as a second route, same items
    EngTrace-certification-guide.pdf   the guide, and guide.md, the same text
    README.txt                         the two routes in a dozen lines

and zips it to dist/kit_<id>.zip. The page runs no code but its own script, reads and
writes nothing but the browser's local storage, and produces rows in exactly the shape
app.py writes (source: "html"), so score.py reads them unchanged. The keyfile is never
touched; the packager check for plant markers applies here too.
"""
from __future__ import annotations

import argparse
import json
import shutil
import zipfile
from pathlib import Path

HERE = Path(__file__).resolve().parent
TASKS = HERE / 'tasks'
DIST = HERE / 'dist'
GUIDE = HERE / 'EngTrace-certification-guide.pdf'

README = """EngTrace template certification - {id}

1. Read the guide: EngTrace-certification-guide.pdf (or guide.md, the same text).

2. Annotate. Two ways; pick one, both record the same thing.

   A. THE PAGE (nothing to install)
      Open {id}.html in Chrome, Edge or Firefox (double-click it). Work through your
      {n} items. The page saves after every item; to pause, close it and reopen the
      same file in the same browser to continue. When it says you are done, click
      "Download my answers" and send the file {id}.jsonl to the coordinator.

   B. THE APP (if you prefer a desktop app; needs Python 3.11 or newer)
      pip install streamlit
      streamlit run app/app.py
      Pick your id in the sidebar. Your answers are written to app/labels/{id}.jsonl -
      send that file back when you finish.

Do not use both routes for the same item. Questions go to the coordinator.
"""

DEFECTS = ['physics or scenario implausible', 'governing equation or formula', 'constant or table value',
           'unit or conversion', 'sign or direction', 'arithmetic: a step does not follow',
           'final answer wrong or wrong unit', 'question ambiguous or unsolvable', 'wording or formatting']

PAGE = r"""<!DOCTYPE html>
<html lang="en"><head><meta charset="utf-8"><title>EngTrace certification __AID__</title>
<meta name="viewport" content="width=device-width, initial-scale=1">
<style>
:root{--ink:#1b2430;--muted:#5f6b7a;--line:#d9dee5;--bg:#f6f7f9;--card:#fff;--accent:#1b3a5c;--ok:#1f7a3d;--warn:#a5641a;--bad:#a3282f}
*{box-sizing:border-box}body{margin:0;font:15px/1.5 system-ui,Segoe UI,Helvetica,Arial,sans-serif;color:var(--ink);background:var(--bg)}
header{background:var(--accent);color:#fff;padding:12px 20px;display:flex;justify-content:space-between;align-items:center;gap:16px;flex-wrap:wrap}
header h1{font-size:17px;margin:0;font-weight:600}header .prog{font-size:14px;opacity:.95}
main{display:grid;grid-template-columns:240px 1fr;gap:0;min-height:calc(100vh - 52px)}
nav{background:#fff;border-right:1px solid var(--line);padding:12px;overflow:auto;max-height:calc(100vh - 52px);position:sticky;top:0}
nav button{display:block;width:100%;text-align:left;border:0;background:none;padding:6px 8px;border-radius:6px;font:inherit;font-size:13px;color:var(--ink);cursor:pointer}
nav button.done{color:var(--muted)}nav button.done::after{content:" \2713";color:var(--ok)}nav button.cur{background:#e8eef6;font-weight:600}
section{padding:18px 24px;max-width:1200px}
.card{background:var(--card);border:1px solid var(--line);border-radius:10px;padding:16px 18px;margin-bottom:14px}
.two{display:grid;grid-template-columns:1fr 1.4fr;gap:16px}@media(max-width:900px){main{grid-template-columns:1fr}nav{position:static;max-height:none}.two{grid-template-columns:1fr}}
.md{white-space:normal;overflow-wrap:anywhere}.md strong{color:var(--accent)}
pre{background:#f1f3f6;border:1px solid var(--line);border-radius:8px;padding:12px;overflow:auto;font-size:12.5px;line-height:1.45}
textarea,input[type=text]{width:100%;font:inherit;padding:8px;border:1px solid var(--line);border-radius:8px}
button.primary{background:var(--accent);color:#fff;border:0;padding:10px 16px;border-radius:8px;font:inherit;font-weight:600;cursor:pointer}
button.primary:disabled{opacity:.45;cursor:not-allowed}button.plain{background:#fff;border:1px solid var(--line);padding:8px 12px;border-radius:8px;font:inherit;cursor:pointer}
.banner{padding:10px 14px;border-radius:8px;margin:10px 0;font-weight:600}.ok{background:#e6f4ea;color:var(--ok)}.warn{background:#fbf1e0;color:var(--warn)}.info{background:#e8eef6;color:var(--accent)}
.row{display:flex;gap:18px;flex-wrap:wrap;align-items:flex-start}.field{min-width:220px;flex:1}.field label{display:block;font-weight:600;margin-bottom:4px}
.field output{font-weight:700;margin-left:8px}input[type=range]{width:100%}
.defects label{display:block;font-weight:400;margin:2px 0}.hint{color:var(--muted);font-size:13px}
.tabs button{margin-right:6px}.tabs button.on{background:#e8eef6;font-weight:600}
details summary{cursor:pointer;font-weight:600}
</style></head><body>
<header><h1>EngTrace template certification &middot; reviewer <span id="aid"></span></h1><div class="prog" id="prog"></div>
<div><button class="plain" id="dl">Download my answers</button></div></header>
<main><nav id="nav"></nav><section id="view"></section></main>
<script type="application/json" id="data">__DATA__</script>
<script>
const D = JSON.parse(document.getElementById('data').textContent);
const AID = D.aid, CODES = D.codes, POOL = D.pool, DEFECTS = D.defects, BRANCH = D.branch;
const KEY = 'engtrace_labels_' + AID;
const NUM = /[-+]?\d[\d,]*\.?\d*(?:[eE][-+]?\d+)?/g;
let labels = {}; try { labels = JSON.parse(localStorage.getItem(KEY) || '{}'); } catch (e) { labels = {}; }
let S = { code: null, stage: 'hand', inst: 0, opened: null, hand: '', handAt: null };

function now() { return new Date().toISOString().replace(/\.\d{3}Z$/, 'Z'); }
function esc(s) { return String(s).replace(/&/g, '&amp;').replace(/</g, '&lt;').replace(/>/g, '&gt;'); }
function md(s) {
  let t = esc(s).replace(/\\times/g, '&times;').replace(/\\%/g, '%').replace(/\$/g, '');
  t = t.replace(/\*\*([^*\n]+)\*\*/g, '<strong>$1</strong>').replace(/`([^`\n]+)`/g, '<code>$1</code>');
  return t.replace(/\n/g, '<br>');
}
function numbers(text) { const out = []; for (const m of (text || '').match(NUM) || []) { const v = parseFloat(m.replace(/,/g, '')); if (!isNaN(v)) out.push(v); } return out; }
function handMatch(text, gold) {
  const mine = numbers(text); if (!mine.length || !gold.length) return null;
  for (const m of mine) for (const g of gold) { if (g === 0 && Math.abs(m) < 1e-9) return true; if (g !== 0 && Math.abs(m - g) / Math.abs(g) <= 0.01) return true; }
  return false;
}
function save() { localStorage.setItem(KEY, JSON.stringify(labels)); }
function nextCode() { return CODES.find(c => !labels[c]) || null; }
function open(code) { S = { code, stage: 'hand', inst: 0, opened: now(), hand: '', handAt: null }; render(); }

function renderNav() {
  const done = Object.keys(labels).length;
  document.getElementById('aid').textContent = AID;
  document.getElementById('prog').textContent = done + ' of ' + CODES.length + ' submitted';
  document.getElementById('nav').innerHTML = CODES.map((c, i) =>
    '<button data-c="' + c + '" class="' + (labels[c] ? 'done ' : '') + (c === S.code ? 'cur' : '') + '">' + (i + 1) + '. ' + c + '</button>').join('');
  document.querySelectorAll('#nav button').forEach(b => b.onclick = () => open(b.dataset.c));
}

function render() {
  renderNav();
  const v = document.getElementById('view');
  if (!S.code) {
    v.innerHTML = '<div class="card"><h2>All ' + CODES.length + ' items submitted. Thank you.</h2><p>Click <b>Download my answers</b> at the top and send the file <code>' + AID + '.jsonl</code> to the coordinator.</p></div>';
    return;
  }
  const item = POOL[S.code], idx = CODES.indexOf(S.code), inst = item.instances[S.inst], prev = labels[S.code];
  let h = '<details class="card"><summary>How this works (read once)</summary><p>For each item: first solve the question yourself and enter your final answer, then the solution appears and tells you whether your number matches within 1%. Read the solution, open the source code if a number looks odd, cycle through other instances if useful, then score the three dimensions and Approve or Reject. A rejection needs a defect type and a note saying what is wrong and where. Some items in your queue are quality-control items with a known defect; review every item the same way. Your work saves in this browser after each item.</p></details>';
  h += '<div class="card"><h2 style="margin:0 0 4px">Item ' + (idx + 1) + ' of ' + CODES.length + ' &middot; ' + S.code + ' &middot; ' + esc(item.area.replace(/_/g, ' ')) + '</h2>';
  if (prev) h += '<div class="banner info">Already submitted (' + prev.decision + ', ' + prev.submitted_at + '). Submitting again replaces it.</div>';
  if (S.stage === 'hand') {
    h += '<h3>Hand check</h3><p class="hint">Solve this instance yourself and enter the final answer(s) with unit before seeing the solution.</p>';
    h += '<div class="md card" style="background:#fbfbfc">' + md(inst.question) + '</div>';
    h += '<textarea id="hand" rows="3" placeholder="e.g. 24.05 kN (tension); or a classification in words">' + esc(S.hand) + '</textarea>';
    h += '<p><button class="primary" id="show" disabled>Show solution</button></p></div>';
    v.innerHTML = h;
    const ta = document.getElementById('hand'), btn = document.getElementById('show');
    ta.oninput = () => { S.hand = ta.value; btn.disabled = !ta.value.trim(); };
    btn.disabled = !S.hand.trim();
    btn.onclick = () => { S.stage = 'review'; S.handAt = now(); S.match = handMatch(S.hand, item.gold_numbers); render(); window.scrollTo(0, 0); };
    return;
  }
  const m = S.match;
  h += m === true ? '<div class="banner ok">Your answer matches the template\'s answer within 1%.</div>'
     : m === false ? '<div class="banner warn">Your answer differs from the template\'s by more than 1%. Look for which side is right.</div>'
     : '<div class="banner info">No number to compare (a non-numeric answer). Judge the solution directly.</div>';
  h += '<p class="hint">You entered: ' + esc(S.hand) + '</p>';
  h += '<div class="tabs"><button class="plain on" id="tabSol">Question and solution</button><button class="plain" id="tabSrc">Source code</button> <button class="plain" id="another">Another instance (' + (S.inst + 1) + ' of ' + item.instances.length + ')</button></div>';
  h += '<div id="paneSol" class="two" style="margin-top:12px"><div class="md card" style="background:#fbfbfc"><strong>Question</strong><br>' + md(inst.question) + '</div><div class="md card" style="background:#fbfbfc"><strong>Solution</strong><br>' + md(inst.solution) + '</div></div>';
  h += '<div id="paneSrc" style="display:none;margin-top:12px"><pre>' + esc(item.source) + '</pre></div></div>';
  h += '<div class="card"><h3 style="margin-top:0">Your verdict</h3><div class="row">';
  for (const [k, lab] of [['phys', 'Physical plausibility'], ['math', 'Mathematical correctness'], ['ped', 'Pedagogical clarity']]) {
    const val = prev ? prev.scores[{phys:'physical_plausibility',math:'mathematical_correctness',ped:'pedagogical_clarity'}[k]] : 5;
    h += '<div class="field"><label>' + lab + '<output id="o_' + k + '">' + val + '</output></label><input type="range" min="1" max="5" step="1" id="s_' + k + '" value="' + val + '"></div>';
  }
  h += '</div><p><label><input type="radio" name="dec" value="Approve" ' + (!prev || prev.decision === 'Approve' ? 'checked' : '') + '> Approve</label> &nbsp; <label><input type="radio" name="dec" value="Reject" ' + (prev && prev.decision === 'Reject' ? 'checked' : '') + '> Reject</label></p>';
  h += '<div class="defects"><b>Defect type</b> (required on Reject)' + DEFECTS.map(d => '<label><input type="checkbox" value="' + esc(d) + '" ' + (prev && prev.defects.includes(d) ? 'checked' : '') + '> ' + esc(d) + '</label>').join('') + '</div>';
  h += '<p><label><b>Note</b> (required on Reject: what is wrong and where; optional on Approve)</label><textarea id="note" rows="3">' + esc(prev ? prev.feedback : '') + '</textarea></p>';
  h += '<div class="field" style="max-width:320px"><label>Your confidence in this verdict<output id="o_conf">' + (prev ? prev.confidence : 4) + '</output></label><input type="range" min="1" max="5" step="1" id="s_conf" value="' + (prev ? prev.confidence : 4) + '"></div>';
  h += '<p class="hint" id="why"></p><p><button class="primary" id="submit">Submit and continue</button></p></div>';
  v.innerHTML = h;
  for (const k of ['phys', 'math', 'ped', 'conf']) document.getElementById('s_' + k).oninput = e => document.getElementById('o_' + k).textContent = e.target.value;
  document.getElementById('tabSol').onclick = () => { document.getElementById('paneSol').style.display = ''; document.getElementById('paneSrc').style.display = 'none'; };
  document.getElementById('tabSrc').onclick = () => { document.getElementById('paneSol').style.display = 'none'; document.getElementById('paneSrc').style.display = ''; };
  document.getElementById('another').onclick = () => { S.inst = (S.inst + 1) % item.instances.length; render(); };
  document.getElementById('submit').onclick = () => {
    const dec = document.querySelector('input[name=dec]:checked').value;
    const defects = [...document.querySelectorAll('.defects input:checked')].map(x => x.value);
    const note = document.getElementById('note').value.trim();
    if (dec === 'Reject' && (!note || !defects.length)) { document.getElementById('why').textContent = 'A rejection needs a note and at least one defect type.'; return; }
    labels[S.code] = { annotator_id: AID, branch: BRANCH, code: S.code, position: idx + 1, opened_at: S.opened,
      hand_answer: S.hand, hand_numbers: numbers(S.hand), hand_match: S.match, hand_submitted_at: S.handAt,
      hand_instance_seed: item.instances[0].seed, instances_viewed: S.inst + 1,
      scores: { physical_plausibility: +document.getElementById('s_phys').value, mathematical_correctness: +document.getElementById('s_math').value, pedagogical_clarity: +document.getElementById('s_ped').value },
      decision: dec, defects, feedback: note, confidence: +document.getElementById('s_conf').value, submitted_at: now(), source: 'html' };
    save();
    const n = nextCode(); if (n) open(n); else { S.code = null; render(); }
    window.scrollTo(0, 0);
  };
}
document.getElementById('dl').onclick = () => {
  const rows = CODES.filter(c => labels[c]).map(c => JSON.stringify(labels[c])).join('\n') + '\n';
  const a = document.createElement('a'); a.href = URL.createObjectURL(new Blob([rows], { type: 'application/json' })); a.download = AID + '.jsonl'; a.click();
};
open(nextCode());
</script></body></html>
"""


def build_one(aid: str, pool: dict, assignment: dict) -> Path:
    mine = assignment[aid]
    sub = {c: pool[c] for c in mine['codes']}
    payload = {'aid': aid, 'branch': mine['branch'], 'codes': mine['codes'], 'pool': sub, 'defects': DEFECTS}
    blob = json.dumps(payload, ensure_ascii=False)
    if 'plant_' in blob or 'keyfile' in blob:
        raise SystemExit(f'refusing: {aid} payload would reveal a plant')
    blob = blob.replace('<', '\\u003c').replace('\u2028', '\\u2028').replace('\u2029', '\\u2029')
    html = PAGE.replace('__AID__', aid).replace('__DATA__', blob)
    kit = DIST / f'kit_{aid}'
    if kit.exists():
        shutil.rmtree(kit)
    kit.mkdir(parents=True)
    (kit / f'{aid}.html').write_text(html, encoding='utf8')
    shutil.copy(GUIDE, kit / 'EngTrace-certification-guide.pdf')
    shutil.copy(HERE / 'guide.md', kit / 'guide.md')
    (kit / 'README.txt').write_text(README.format(id=aid, n=len(mine['codes'])), encoding='utf8')
    # the Streamlit app as the second route, with only this expert's items
    app = kit / 'app'
    (app / 'tasks').mkdir(parents=True)
    shutil.copy(HERE / 'app.py', app / 'app.py')
    (app / 'tasks' / 'pool.json').write_text(json.dumps(sub, ensure_ascii=False), encoding='utf8')
    (app / 'tasks' / 'assignment.json').write_text(json.dumps({aid: mine}, indent=1), encoding='utf8')
    out = DIST / f'kit_{aid}.zip'
    with zipfile.ZipFile(out, 'w', zipfile.ZIP_DEFLATED) as z:
        for p in sorted(kit.rglob('*')):
            if p.is_file():
                z.write(p, f'kit_{aid}/{p.relative_to(kit).as_posix()}')
        if any('keyfile' in n for n in z.namelist()):
            raise SystemExit('refusing: bundle contains a keyfile')
    return out


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument('--id', default=None)
    a = ap.parse_args()
    if not GUIDE.exists():
        raise SystemExit('typeset the guide first: python -m template_annotation_23092026.layer2.make_guide_pdf')
    pool = json.loads((TASKS / 'pool.json').read_text(encoding='utf8'))
    assignment = json.loads((TASKS / 'assignment.json').read_text(encoding='utf8'))
    for aid in ([a.id] if a.id else sorted(assignment)):
        out = build_one(aid, pool, assignment)
        print(f'{aid}: {out.name} ({out.stat().st_size // 1024} KB)')


if __name__ == '__main__':
    main()
