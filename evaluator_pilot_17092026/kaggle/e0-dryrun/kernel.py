"""E0 dry run on a Kaggle GPU: Tier 1 and BERTScore for all 300 traces, no judges.

Pushed with:
    python -m kaggle kernels push -p evaluator_pilot_17092026/kaggle/e0-dryrun --accelerator NvidiaTeslaT4

WHY HERE. Tier 1 runs the stsb-roberta-large cross-encoder over every (gold step,
trace step) pair and Longformer BERTScore over every trace. On the laptop CPU that
is 0.5-2s per pair and hours overall; on a GPU it is minutes. The framework picks
CUDA by itself (`device = 'cuda' if torch.cuda.is_available()`), so nothing in the
scoring code changes for this.

WHAT DOES NOT COME HERE. No API key: a dry run never reaches a judge, and the paid
Tribunal runs locally. The bundle holds the evaluation code, the frozen slice and
the traces - no .env, no testset, no template tree.

WHAT HAS TO MATCH. Every library in the scorer cache key is pinned to the local
venv's exact version, so a cache entry written here is keyed identically to one
written on the laptop. torch stays Kaggle's CUDA build; device and torch version
are recorded per row, and `run_evaluator --import-kaggle` refuses the output unless
it agrees with rows already scored on the laptop CPU.
"""
import json
import os
import shutil
import subprocess
import sys
import time

OUT = '/kaggle/working/out'
os.makedirs(OUT, exist_ok=True)


def sh(cmd, **kw):
    print('$ ' + ' '.join(cmd), flush=True)
    return subprocess.run(cmd, check=True, **kw)


# Kaggle may or may not unzip the bundle, so find the marker instead of assuming.
root = next(d for d, _, f in os.walk('/kaggle/input') if 'ENGTRACE_BUNDLE' in f)
shutil.copytree(root, '/tmp/repo', dirs_exist_ok=True)
# --dir-mode zip uploads each subdirectory as a zip; extract any that arrived packed.
for name in os.listdir('/tmp/repo'):
    if name.endswith('.zip') and not os.path.isdir(os.path.join('/tmp/repo', name[:-4])):
        shutil.unpack_archive(os.path.join('/tmp/repo', name), os.path.join('/tmp/repo', name[:-4]))

# The bundle records a SHA-256 per file; check the kernel is reading exactly those bytes.
import hashlib
with open(os.path.join('/tmp/repo', 'ENGTRACE_BUNDLE')) as fh:
    expected = json.load(fh)
moved = [rel for rel, sha in expected.items()
         if hashlib.sha256(open(os.path.join('/tmp/repo', rel), 'rb').read()).hexdigest() != sha]
print('bundle: %d files, %d mismatched' % (len(expected), len(moved)), flush=True)
if moved:
    raise SystemExit('bundle files differ from what was staged: %s' % moved)

PINS = [
    'transformers==4.57.3', 'sentence-transformers==5.1.2', 'bert-score==0.3.13',
    'rouge-score==0.1.2', 'tokenizers==0.22.2', 'scipy==1.18.1', 'numpy==2.5.3',
    'google-generativeai==0.8.6', 'python-dotenv==1.2.3', 'openai', 'anthropic',
]
sh([sys.executable, '-m', 'pip', 'install', '-q', *PINS])

# A GPU that torch cannot actually use would silently send the framework to CPU,
# or fail mid-run. Prove a real CUDA op first, and stop if it does not work.
probe = subprocess.run([sys.executable, '-c', (
    'import torch,json;'
    'ok=torch.cuda.is_available();'
    'x=(torch.ones(4,device="cuda")*2).sum().item() if ok else None;'
    'print(json.dumps({"torch":torch.__version__,"cuda":ok,"device":torch.cuda.get_device_name(0) if ok else None,"op":x}))'
)], capture_output=True, text=True)
print(probe.stdout, probe.stderr[-2000:], flush=True)
info = json.loads(probe.stdout.strip().splitlines()[-1])
with open(os.path.join(OUT, 'gpu.json'), 'w') as fh:
    json.dump(info, fh, indent=2)
if not info.get('cuda') or info.get('op') != 8.0:
    raise SystemExit('no usable CUDA device - refusing to run on CPU')

with open(os.path.join(OUT, 'pip_freeze.txt'), 'w') as fh:
    fh.write(subprocess.run([sys.executable, '-m', 'pip', 'freeze'], capture_output=True, text=True).stdout)

env = {**os.environ, 'PYTHONPATH': '/tmp/repo', 'PYTHONIOENCODING': 'utf-8',
       'HF_HUB_OFFLINE': '0', 'TRANSFORMERS_OFFLINE': '0'}
t0 = time.time()
rc = subprocess.run([sys.executable, '-m', 'evaluator_pilot_17092026.run_evaluator', 'e0',
                     '--dry-run', '--freeze-check', 'hash'], cwd='/tmp/repo', env=env).returncode

scores = '/tmp/repo/evaluator_pilot_17092026/scores'
if os.path.isdir(scores):
    shutil.copytree(scores, os.path.join(OUT, 'scores'), dirs_exist_ok=True)
with open(os.path.join(OUT, 'run.json'), 'w') as fh:
    json.dump({'returncode': rc, 'seconds': round(time.time() - t0)}, fh, indent=2)
print('run_evaluator exited %d after %.0fs' % (rc, time.time() - t0), flush=True)
