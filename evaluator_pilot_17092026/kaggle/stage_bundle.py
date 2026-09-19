"""Stage the private Kaggle bundle for the evaluator pilot.

    python evaluator_pilot_17092026/kaggle/stage_bundle.py STAGE_DIR

Copies ONLY what a keyless dry run reads, into a directory outside the repo:

    evaluation/engtrace_evaluation_framework.py   the published framework, byte-identical
    evaluation/engineering_parser.py              its parser
    evaluation/run_inference.py                   run_traces reads the deployed prompt from it
    evaluator_pilot_17092026/{run_evaluator,run_traces,verify_traces}.py, models.json
    evaluator_pilot_17092026/evaluators/e0_tribunal.py
    evaluator_pilot_17092026/slice/               the frozen manifest + FREEZE.json
    evaluator_pilot_17092026/traces/*.jsonl       the five live trace files (not superseded/)

Never copied: .env or any key, testset/, data/, pilot_new_branches/, scores/, the venv.

Before copying, the full freeze rebuild and verify_traces must pass HERE. The
kernel can only check the manifest hash; this is what makes that hash mean the
slice was rebuilt from the templates and matched.
"""
import hashlib
import json
import os
import shutil
import subprocess
import sys

ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
P = 'evaluator_pilot_17092026'
FILES = [
    'evaluation/engtrace_evaluation_framework.py',
    'evaluation/engineering_parser.py',
    'evaluation/run_inference.py',
    P + '/run_evaluator.py', P + '/run_traces.py', P + '/verify_traces.py', P + '/models.json',
    P + '/evaluators/e0_tribunal.py',
    P + '/slice/manifest.jsonl', P + '/slice/FREEZE.json',
] + [P + '/traces/%s.jsonl' % k for k in
     ('gpt-5', 'claude-opus-4.7', 'gemini-3.1-pro', 'deepseek-r1', 'llama-3.1-70b')]


MODES = ('dry', 'capture', 'replay')


def main():
    """stage_bundle.py STAGE_DIR [MODE]. MODE is dry (E0 dry run, default), capture
    or replay (the two Kaggle halves of E1). capture and replay also ship E1's code
    and the scorer caches - 724 KB, and they spare the kernel recomputing Tier 1;
    replay additionally ships the judge replies fetched on the laptop and the price
    snapshot, so the kernel needs no key at any point."""
    if len(sys.argv) not in (2, 3):
        raise SystemExit(__doc__)
    stage = os.path.abspath(sys.argv[1])
    mode = sys.argv[2] if len(sys.argv) == 3 else 'dry'
    if mode not in MODES:
        raise SystemExit('mode must be one of %s' % (MODES,))
    files = list(FILES)
    if mode in ('capture', 'replay'):
        files.append(P + '/evaluators/e1_panel.py')
        cache = os.path.join(ROOT, P, 'scores', '_cache')
        files += [P + '/scores/_cache/' + f for f in sorted(os.listdir(cache))
                  if f.endswith('.jsonl') or f == 'openrouter_prices.json']
        if mode == 'replay' and not os.path.exists(os.path.join(cache, 'e1_judge_replies.jsonl')):
            raise SystemExit('replay needs the reply store: run --fetch-captured first')
    if stage.startswith(ROOT):
        raise SystemExit('stage outside the repo, not in it')

    py = sys.executable
    for mod in ('%s.freeze' % P, '%s.verify_traces' % P):
        # Only the gold five are bundled and scored, so the gold five are what must
        # verify; the robustness cohort is reported elsewhere and never shipped.
        args = [py, '-m', mod] + (['--verify'] if mod.endswith('freeze') else ['--cohort', 'gold'])
        r = subprocess.run(args, cwd=ROOT, capture_output=True, text=True,
                           env={**os.environ, 'PYTHONIOENCODING': 'utf-8'})
        print(r.stdout.strip().splitlines()[-1] if r.stdout.strip() else r.stderr[-500:])
        if r.returncode != 0:
            raise SystemExit('%s failed - not staging' % mod)

    if os.path.isdir(stage):
        shutil.rmtree(stage)
    manifest = {}
    for rel in files:
        src = os.path.join(ROOT, rel)
        dst = os.path.join(stage, rel)
        os.makedirs(os.path.dirname(dst), exist_ok=True)
        shutil.copyfile(src, dst)
        manifest[rel] = hashlib.sha256(open(src, 'rb').read()).hexdigest()

    # Load the real keys before scanning for them. Without this the environment
    # holds none, and the scan passes while looking for nothing.
    from dotenv import dotenv_values
    keys = {k: v for k, v in dotenv_values(os.path.join(ROOT, '.env')).items()
            if v and len(v) >= 16 and ('KEY' in k or 'TOKEN' in k)}
    if not keys:
        raise SystemExit('read no keys from .env - the leak scan would check nothing')
    for bad in ('.env', 'testset', 'data'):
        if os.path.exists(os.path.join(stage, bad)):
            raise SystemExit('refusing: %s ended up in the stage' % bad)
    for dirpath, _, names in os.walk(stage):      # not `files`: that is the bundle list
        for f in names:
            text = open(os.path.join(dirpath, f), 'rb').read()
            for name, v in keys.items():
                if v.encode() in text:
                    raise SystemExit('refusing: the value of %s is inside %s' % (name, f))
    print('leak scan: %d key values from .env, none found in the stage' % len(keys))

    with open(os.path.join(stage, 'ENGTRACE_BUNDLE'), 'w', encoding='utf-8') as fh:
        json.dump(manifest, fh, indent=2)
    with open(os.path.join(stage, 'KAGGLE_MODE'), 'w', encoding='utf-8') as fh:
        fh.write(mode)
    with open(os.path.join(stage, 'dataset-metadata.json'), 'w', encoding='utf-8') as fh:
        json.dump({'title': 'engtrace-evaluator-pilot', 'id': 'ayeshaiq/engtrace-evaluator-pilot',
                   'licenses': [{'name': 'other'}]}, fh, indent=2)
    size = sum(os.path.getsize(os.path.join(d, f)) for d, _, fs in os.walk(stage) for f in fs)
    print('staged %d files, %.1f MB, mode %s -> %s' % (len(files), size / 1e6, mode, stage))


if __name__ == '__main__':
    main()
