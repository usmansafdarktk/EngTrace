"""E2 on a GPU: score every step of every trace with an open process reward model.

    python e2_score.py --prm qwen72|qwen7|versa --inputs e2_inputs.jsonl --out OUT.jsonl
                       [--limit N] [--selfcheck-only]

Runs on HiPerGator (hpg/e2.sbatch). No API key is involved anywhere: the PRM runs
locally on the node, and weights are fetched to node-local scratch, never to /blue.

EACH PRM IS RUN EXACTLY AS ITS MODEL CARD SPECIFIES (the code below follows the
cards' published usage line by line):

  qwen72 / qwen7   Qwen/Qwen2.5-Math-PRM-{72B,7B}. Chat template with the card's system
                   prompt; the assistant turn is the steps joined by "<extra_0>" and
                   ending with it; per-step reward = softmax over the 2 logits at each
                   "<extra_0>" position, index 1.
  versa            UW-Madison-Lee-Lab/VersaPRM, a LoRA adapter on Llama-PRM800K (base pinned too).
                   question + ' \\n\\n' + ' \\n\\n\\n\\n'.join(steps) + ' \\n\\n\\n\\n';
                   per-step score = softmax over logits of tokens [12, 10] at each
                   token 23535, index 1.

SELF-CHECK FIRST. Both Qwen cards print the rewards their example must produce. The run
refuses to score anything unless this machine gives every step the card's verdict at 0.5,
and it records how far each reward is from the card's (see VERDICT_AT below for why the
check is on verdicts, not on a 0.02 tolerance). VersaPRM's card prints no reference numbers; its check is
that the example yields exactly one score per step and that switching the adapter off
changes the scores (the adapter is really loaded).

WHAT IS RECORDED, per trace: the PRM, its repo and pinned revision, n_steps, n_scores,
every step probability, input token count against the model's context limit, and the
item and trace hashes. A trace longer than the context is marked over_length and NOT
scored - a silent truncation would drop its last steps and look like a clean score. A
trace whose score count differs from its step count is marked step_mismatch.
"""
import argparse
import hashlib
import json
import os
import platform
import sys
import time

import torch
import torch.nn.functional as F

REPOS = {
    'qwen72': ('Qwen/Qwen2.5-Math-PRM-72B', '9df429b02adb5f764cd6e30e76a0cca16d501ae1'),
    'qwen7': ('Qwen/Qwen2.5-Math-PRM-7B', '0610740060112df12585d00a1c5f4624d2f59051'),
    'versa': ('UW-Madison-Lee-Lab/VersaPRM', 'cd537f770b8077e31b10c7dc6cfd27f582e5de51'),
}
QWEN_SYSTEM = "Please reason step by step, and put your final answer within \\boxed{}."
CARD_EXAMPLE = {
    'query': "Sue lives in a fun neighborhood.  One weekend, the neighbors decided to play a prank on Sue.  On Friday morning, the neighbors placed 18 pink plastic flamingos out on Sue's front yard.  On Saturday morning, the neighbors took back one third of the flamingos, painted them white, and put these newly painted white flamingos back out on Sue's front yard.  Then, on Sunday morning, they added another 18 pink plastic flamingos to the collection. At noon on Sunday, how many more pink plastic flamingos were out than white plastic flamingos?",
    'response': [
        "To find out how many more pink plastic flamingos were out than white plastic flamingos at noon on Sunday, we can break down the problem into steps. First, on Friday, the neighbors start with 18 pink plastic flamingos.",
        "On Saturday, they take back one third of the flamingos. Since there were 18 flamingos, (1/3 \\times 18 = 6) flamingos are taken back. So, they have (18 - 6 = 12) flamingos left in their possession. Then, they paint these 6 flamingos white and put them back out on Sue's front yard. Now, Sue has the original 12 pink flamingos plus the 6 new white ones. Thus, by the end of Saturday, Sue has (12 + 6 = 18) pink flamingos and 6 white flamingos.",
        "On Sunday, the neighbors add another 18 pink plastic flamingos to Sue's front yard. By the end of Sunday morning, Sue has (18 + 18 = 36) pink flamingos and still 6 white flamingos.",
        "To find the difference, subtract the number of white flamingos from the number of pink flamingos: (36 - 6 = 30). Therefore, at noon on Sunday, there were 30 more pink plastic flamingos out than white plastic flamingos. The answer is (\\boxed{30}).",
    ],
}
CARD_EXPECTED = {'qwen72': [0.9921875, 0.0047607421875, 0.32421875, 0.8203125],
                 'qwen7': [1.0, 0.1904296875, 0.9765625, 1.0]}
# The self-check PASSES when every step of the card example gets the card's verdict at
# the 0.5 threshold PRMs are used with; the continuous gap is recorded beside it
# (max_abs_diff, within_tol). The first design demanded 0.02 on every reward. On this
# hardware (Blackwell, no flash-attn) the 72B misses by up to 0.26 on the card's two
# borderline steps while every verdict agrees, and the 7B, with the same code, lands
# within 0.033 of its card. That gap was traced to kernel numerics, not to the pipeline
# (hpg/README.md, D-090), and the user chose to proceed with it recorded as a deviation.
SELFCHECK_TOL = 0.02
VERDICT_AT = 0.5


class Qwen:
    def __init__(self, repo, rev, attn=None):
        from transformers import AutoModel, AutoTokenizer
        self.tok = AutoTokenizer.from_pretrained(repo, revision=rev, trust_remote_code=True)
        kw = {'attn_implementation': attn} if attn else {}
        self.model = AutoModel.from_pretrained(repo, revision=rev, device_map='auto',
                                               torch_dtype=torch.bfloat16, trust_remote_code=True, **kw).eval()
        self.attn = getattr(self.model.config, '_attn_implementation', None)
        self.sep = self.tok.encode('<extra_0>')[0]
        self.max_len = getattr(self.model.config, 'max_position_embeddings', None)
        # Recorded because it decides fidelity: the card's modeling_qwen2_rm.py builds
        # positions as arange(...).type_as(inv_freq). Under transformers 4.57 this buffer
        # ended up bf16 (positions past 256 not exact) and the self-check failed with an
        # error growing along the sequence; see e2.sbatch for the pinned version.
        self.rotary = sorted({str(b.dtype) for n, b in self.model.named_buffers() if n.endswith('inv_freq')})

    def ids(self, question, steps):
        messages = [{'role': 'system', 'content': QWEN_SYSTEM},
                    {'role': 'user', 'content': question},
                    {'role': 'assistant', 'content': '<extra_0>'.join(steps) + '<extra_0>'}]
        conv = self.tok.apply_chat_template(messages, tokenize=False, add_generation_prompt=False)
        return self.tok.encode(conv, return_tensors='pt')

    @torch.no_grad()
    def score(self, input_ids, card_exact=False):
        input_ids = input_ids.to(self.model.device)
        logits = self.model(input_ids=input_ids)[0]
        # the card's make_step_rewards, for one sample; the card takes the softmax in
        # bf16, we take it in fp32 (card_exact reproduces the card's arithmetic)
        probs = F.softmax(logits if card_exact else logits.float(), dim=-1) * (input_ids == self.sep).unsqueeze(-1)
        sample = probs[0]
        return sample[sample != 0].view(-1, 2)[:, 1].cpu().tolist()


class Versa:
    STEP_TOKEN, CANDIDATES = 23535, [12, 10]
    # The adapter's adapter_config.json names its base by repo id only. Loading the
    # adapter through AutoModelForCausalLM passes the ADAPTER's revision on to the base
    # repo, where that commit does not exist (smoke attempt 3 stopped on exactly that).
    # So the base is loaded explicitly at its own pinned commit, then the adapter on top.
    BASE = ('UW-Madison-Lee-Lab/Llama-PRM800K', '1973a85d64b7a00e50a278b31a6558804bf952eb')

    def __init__(self, repo, rev):
        from peft import PeftConfig, PeftModel
        from transformers import AutoTokenizer, LlamaForCausalLM
        cfg = PeftConfig.from_pretrained(repo, revision=rev)
        if cfg.base_model_name_or_path != self.BASE[0]:
            raise SystemExit('adapter names base %r, expected %r' % (cfg.base_model_name_or_path, self.BASE[0]))
        self.tok = AutoTokenizer.from_pretrained(repo, revision=rev)
        self.tok.pad_token = self.tok.eos_token
        base = LlamaForCausalLM.from_pretrained(self.BASE[0], revision=self.BASE[1],
                                                torch_dtype=torch.bfloat16, device_map='auto')
        self.model = PeftModel.from_pretrained(base, repo, revision=rev).eval()
        self.max_len = getattr(base.config, 'max_position_embeddings', None)
        self.base = {'repo': self.BASE[0], 'revision': self.BASE[1]}

    def ids(self, question, steps):
        text = question + ' \n\n' + ' \n\n\n\n'.join(steps) + ' \n\n\n\n'
        return torch.tensor([self.tok.encode(text)])

    @torch.no_grad()
    def score(self, input_ids):
        input_ids = input_ids.to(self.model.device)
        logits = self.model(input_ids).logits[:, :, self.CANDIDATES]
        scores = logits.float().softmax(dim=-1)[:, :, 1]
        return scores[input_ids == self.STEP_TOKEN].cpu().tolist()


def versions():
    import importlib
    out = {'python': platform.python_version(), 'torch': torch.__version__, 'cuda': torch.version.cuda,
           'gpus': [torch.cuda.get_device_name(i) for i in range(torch.cuda.device_count())]}
    for m in ('transformers', 'peft', 'accelerate', 'huggingface_hub', 'safetensors', 'tokenizers'):
        try:
            out[m] = importlib.import_module(m).__version__
        except Exception:                                          # noqa: BLE001
            out[m] = None
    return out


def selfcheck(prm, runner):
    ids = runner.ids(CARD_EXAMPLE['query'], CARD_EXAMPLE['response'])
    got = runner.score(ids)
    rep = {'prm': prm, 'got': got, 'n_steps': len(CARD_EXAMPLE['response']), 'n_tokens': int(ids.shape[1])}
    if prm in CARD_EXPECTED:
        rep['got_bf16_softmax'] = runner.score(ids, card_exact=True)
        exp = CARD_EXPECTED[prm]
        diff = max(abs(a - b) for a, b in zip(got, exp)) if len(got) == len(exp) else None
        same = len(got) == len(exp) and all((g > VERDICT_AT) == (e > VERDICT_AT) for g, e in zip(got, exp))
        rep.update(expected=exp, max_abs_diff=diff, verdicts_match=same,
                   within_tol=diff is not None and diff <= SELFCHECK_TOL, passed=same)
    else:
        # No reference values on VersaPRM's card. Checked instead: one score per step,
        # each in [0, 1], AND the adapter is really applied - the same input scored with
        # the LoRA adapter switched off must give different scores, or a silently
        # dropped adapter would pass as the bare Llama-PRM800K base.
        with runner.model.disable_adapter():
            base = runner.score(runner.ids(CARD_EXAMPLE['query'], CARD_EXAMPLE['response']))
        shift = max(abs(a - b) for a, b in zip(got, base)) if len(base) == len(got) else None
        rep.update(base_only=base, adapter_shift=shift)
        rep['passed'] = (len(got) == len(CARD_EXAMPLE['response']) and all(0 <= x <= 1 for x in got)
                         and shift is not None and shift > 1e-3)
        rep['note'] = ('no reference values on the card; checked one score per step in [0, 1] '
                       'and that the adapter changes the scores')
    return rep


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--prm', required=True, choices=sorted(REPOS))
    ap.add_argument('--inputs', required=True)
    ap.add_argument('--out', required=True)
    ap.add_argument('--limit', type=int)
    ap.add_argument('--selfcheck-only', action='store_true')
    ap.add_argument('--attn', help='attention implementation (Qwen only); default = the library default')
    a = ap.parse_args()

    repo, rev = REPOS[a.prm]
    t0 = time.time()
    # Fetch every file first with parallel workers: old transformers downloads the 37
    # shards of the 72B one at a time (9 min); this takes a fraction of that.
    from huggingface_hub import snapshot_download
    for r, v in [(repo, rev)] + ([Versa.BASE] if a.prm == 'versa' else []):
        snapshot_download(r, revision=v, max_workers=16)
    runner = Versa(repo, rev) if a.prm == 'versa' else Qwen(repo, rev, a.attn)
    meta = {'prm': a.prm, 'repo': repo, 'revision': rev, 'max_len': runner.max_len,
            'attn_implementation': getattr(runner, 'attn', None) or getattr(runner.model.config, '_attn_implementation', None),
            'base': getattr(runner, 'base', None), 'rotary_inv_freq_dtype': getattr(runner, 'rotary', None),
            'load_seconds': round(time.time() - t0, 1), 'versions': versions(),
            'scorer_md5': hashlib.md5(open(__file__, 'rb').read()).hexdigest(),
            'slurm_job': os.environ.get('SLURM_JOB_ID'), 'host': platform.node()}
    print(json.dumps(meta), flush=True)

    sc = selfcheck(a.prm, runner)
    print('SELFCHECK', json.dumps(sc), flush=True)
    os.makedirs(os.path.dirname(os.path.abspath(a.out)), exist_ok=True)
    with open(a.out + '.meta.json', 'w') as fh:
        json.dump(dict(meta, selfcheck=sc), fh, indent=2)
    if not sc['passed']:
        print('SELFCHECK FAILED - refusing to score: this run does not match the official implementation')
        sys.exit(2)
    if a.selfcheck_only:
        return

    rows = [json.loads(l) for l in open(a.inputs, encoding='utf-8')]
    if a.limit:
        rows = rows[:a.limit]
    n_over = n_mis = 0
    with open(a.out, 'w', encoding='utf-8') as fh:
        for i, r in enumerate(rows, 1):
            ids = runner.ids(r['question'], r['steps'])
            rec = {k: r[k] for k in ('item_id', 'model_key', 'item_sha256', 'trace_sha256')}
            rec.update(prm=a.prm, repo=repo, revision=rev, n_steps=len(r['steps']),
                       n_tokens=int(ids.shape[1]), max_len=runner.max_len)
            if runner.max_len and ids.shape[1] > runner.max_len:
                rec.update(over_length=True, step_probs=None, n_scores=0)
                n_over += 1
            else:
                t1 = time.time()
                probs = runner.score(ids)
                rec.update(over_length=False, step_probs=probs, n_scores=len(probs),
                           step_mismatch=len(probs) != len(r['steps']), seconds=round(time.time() - t1, 2))
                n_mis += rec['step_mismatch']
            fh.write(json.dumps(rec) + '\n')
            if i % 25 == 0 or i == len(rows):
                print('  %d/%d scored, %d over length, %d step mismatches' % (i, len(rows), n_over, n_mis),
                      flush=True)
    print('DONE %s: %d traces, %d over length, %d step mismatches, %.0fs'
          % (a.prm, len(rows), n_over, n_mis, time.time() - t0), flush=True)


if __name__ == '__main__':
    main()
