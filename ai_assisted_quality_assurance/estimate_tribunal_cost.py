"""Estimate what one full AI-Tribunal pass over every template would cost.

Builds the exact prompt `run_ai_tribunal.py` sends (function source + three
generated instances, read from that file's FULL_QA_PROMPT without importing it,
so no API client is constructed and no output directory is created), counts
tokens, and prices one pass per candidate judge.

    python -m ai_assisted_quality_assurance.estimate_tribunal_cost

Prices are OpenRouter list prices in $/M tokens as recorded on 2026-09-19 in
evaluator_pilot_17092026/JUDGE_SELECTION.md and docs/inference_pricing/; edit
PRICES when they move. Token counts use tiktoken's o200k_base when installed
and chars/3.6 otherwise (source code and prose both sit near that ratio), so
treat the totals as +/-20%.
"""
from __future__ import annotations

import ast
import collections
import random
import statistics
import sys
from pathlib import Path

REPO = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(REPO))

from ai_assisted_quality_assurance.template_loader import (  # noqa: E402
    get_areas, get_branches, get_source_code, get_template_files,
    load_template_functions)

OUTPUT_TOKENS_PER_CALL = 200  # the JSON reply with its one-sentence explanation
INSTANCES_PER_TEMPLATE = 3
SEED = 20260922

PRICES = {  # $/M input, $/M output
    'gpt-5 (published judge)': (1.25, 10.00),
    'claude-opus-4.5 (published judge)': (5.00, 25.00),
    'gemini-3.1-pro (successor of the retired gemini-3-pro-preview)': (2.00, 12.00),
    'minimax-m3 (E1 non-suite panel)': (0.30, 1.20),
    'glm-5.3 (E1 non-suite panel)': (1.40, 4.40),
    'kimi-k3 (E1 non-suite panel)': (2.10, 10.95),
}


def prompt_template() -> str:
    src = (REPO / 'ai_assisted_quality_assurance' / 'run_ai_tribunal.py').read_text(encoding='utf8')
    for node in ast.parse(src).body:
        if isinstance(node, ast.Assign) and any(
                isinstance(t, ast.Name) and t.id == 'FULL_QA_PROMPT' for t in node.targets):
            return ast.literal_eval(node.value)
    raise RuntimeError('FULL_QA_PROMPT not found in run_ai_tribunal.py')


def tokenizer():
    try:
        import tiktoken
        enc = tiktoken.get_encoding('o200k_base')
        return (lambda s: len(enc.encode(s))), 'tiktoken o200k_base'
    except Exception:
        return (lambda s: len(s) / 3.6), 'chars/3.6 (tiktoken not installed)'


def main() -> None:
    tok, how = tokenizer()
    prompt = prompt_template()
    random.seed(SEED)
    per_template = []
    for branch in get_branches():
        for area in get_areas(branch):
            for file in get_template_files(branch, area):
                for name, fn in load_template_functions(branch, area, file):
                    inst = [fn() for _ in range(INSTANCES_PER_TEMPLATE)]
                    text = prompt.format(
                        template_code=get_source_code(fn),
                        q1=inst[0][0], s1=inst[0][1], q2=inst[1][0], s2=inst[1][1],
                        q3=inst[2][0], s3=inst[2][1])
                    per_template.append((branch, name, tok(text)))

    n = len(per_template)
    total_in = sum(t for _, _, t in per_template)
    total_out = OUTPUT_TOKENS_PER_CALL * n
    by_branch = collections.Counter()
    for b, _, t in per_template:
        by_branch[b] += t

    print(f'tokenizer: {how}')
    print(f'templates: {n}; instances per template: {INSTANCES_PER_TEMPLATE}')
    print(f'input tokens, one pass, one judge: {total_in:,.0f} '
          f'(median {statistics.median(t for _, _, t in per_template):,.0f} per template, '
          f'max {max(t for _, _, t in per_template):,.0f})')
    for b, t in sorted(by_branch.items()):
        print(f'  {b:24s} {t:>10,.0f}')
    print(f'output tokens, one pass, one judge: {total_out:,} (assumed {OUTPUT_TOKENS_PER_CALL}/call)')
    print('\ncost of one full pass, per judge:')
    for model, (p_in, p_out) in PRICES.items():
        cost = total_in / 1e6 * p_in + total_out / 1e6 * p_out
        print(f'  ${cost:6.2f}  {model}')
    trio = sum(total_in / 1e6 * PRICES[m][0] + total_out / 1e6 * PRICES[m][1]
               for m in list(PRICES)[:3])
    panel = sum(total_in / 1e6 * PRICES[m][0] + total_out / 1e6 * PRICES[m][1]
                for m in list(PRICES)[3:])
    print(f'\npublished trio, one pass: ${trio:.2f}   non-suite panel, one pass: ${panel:.2f}')


if __name__ == '__main__':
    main()
