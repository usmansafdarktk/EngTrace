"""E1 on Kaggle, in two halves, with the API key never leaving the laptop.

    capture   run E1's framework path with the judges stubbed and write every
              judge prompt it would send (captured.jsonl). No key needed.
    replay    run E1 for real, every judge reply served from the reply store the
              laptop fetched. No key needed; a prompt missing from the store is an
              ERROR on its row, never a silently dropped judge.

Between them, on the laptop:
    run_evaluator e1 --fetch-captured captured.jsonl     (plain HTTP, the key lives here)

The mode is read from the bundle's KAGGLE_MODE, written by stage_bundle.py. As in the
E0 dry run: bundle files are checked against the SHA-256 recorded at staging, every
library in the scorer cache key is pinned to the local venv's version, and a real
CUDA op is proved before anything runs.
"""
import hashlib
import json
import os
import shutil
import subprocess
import sys
import time

OUT = '/kaggle/working/out'
os.makedirs(OUT, exist_ok=True)

root = next(d for d, _, f in os.walk('/kaggle/input') if 'ENGTRACE_BUNDLE' in f)
shutil.copytree(root, '/tmp/repo', dirs_exist_ok=True)
for name in os.listdir('/tmp/repo'):
    if name.endswith('.zip') and not os.path.isdir(os.path.join('/tmp/repo', name[:-4])):
        shutil.unpack_archive(os.path.join('/tmp/repo', name), os.path.join('/tmp/repo', name[:-4]))

expected = json.load(open('/tmp/repo/ENGTRACE_BUNDLE'))
moved = [rel for rel, sha in expected.items()
         if hashlib.sha256(open(os.path.join('/tmp/repo', rel), 'rb').read()).hexdigest() != sha]
print('bundle: %d files, %d mismatched' % (len(expected), len(moved)), flush=True)
if moved:
    raise SystemExit('bundle files differ from what was staged: %s' % moved)
mode = open('/tmp/repo/KAGGLE_MODE').read().strip()
print('mode:', mode, flush=True)
if mode not in ('capture', 'replay'):
    raise SystemExit('this kernel runs capture or replay, not %r' % mode)

PINS = [
    'transformers==4.57.3', 'sentence-transformers==5.1.2', 'bert-score==0.3.13',
    'rouge-score==0.1.2', 'tokenizers==0.22.2', 'scipy==1.18.1', 'numpy==2.5.3',
    'google-generativeai==0.8.6', 'python-dotenv==1.2.3', 'openai', 'anthropic',
]
subprocess.run([sys.executable, '-m', 'pip', 'install', '-q', *PINS], check=True)

probe = subprocess.run([sys.executable, '-c', (
    'import torch,json;ok=torch.cuda.is_available();'
    'x=(torch.ones(4,device="cuda")*2).sum().item() if ok else None;'
    'print(json.dumps({"torch":torch.__version__,"cuda":ok,"device":torch.cuda.get_device_name(0) if ok else None,"op":x}))'
)], capture_output=True, text=True)
info = json.loads(probe.stdout.strip().splitlines()[-1])
json.dump(info, open(os.path.join(OUT, 'gpu.json'), 'w'), indent=2)
print(info, flush=True)
if not info.get('cuda') or info.get('op') != 8.0:
    raise SystemExit('no usable CUDA device - refusing to run on CPU')

env = {**os.environ, 'PYTHONPATH': '/tmp/repo', 'PYTHONIOENCODING': 'utf-8', 'PYTHONUNBUFFERED': '1',
       'HF_HUB_OFFLINE': '0', 'TRANSFORMERS_OFFLINE': '0'}
env.pop('OPENROUTER_API_KEY', None)          # belt and braces: this kernel must run key-less
cmd = [sys.executable, '-m', 'evaluator_pilot_17092026.run_evaluator', 'e1', '--freeze-check', 'hash']
cmd += (['--capture', os.path.join(OUT, 'captured.jsonl')] if mode == 'capture' else ['--no-prefetch'])
t0 = time.time()
rc = subprocess.run(cmd, cwd='/tmp/repo', env=env).returncode

scores = '/tmp/repo/evaluator_pilot_17092026/scores'
if mode == 'replay' and os.path.isdir(os.path.join(scores, 'e1')):
    shutil.copytree(os.path.join(scores, 'e1'), os.path.join(OUT, 'e1'), dirs_exist_ok=True)
json.dump({'mode': mode, 'returncode': rc, 'seconds': round(time.time() - t0)},
          open(os.path.join(OUT, 'run.json'), 'w'), indent=2)
print('%s exited %d after %.0fs' % (mode, rc, time.time() - t0), flush=True)
