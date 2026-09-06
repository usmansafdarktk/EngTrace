"""Phase 3 R4 ADOPT-NOW measurements — the reviewer suggestions actioned in-phase.

    python -m tests.trace_schema.phase3_r4_actions

Five measurements, each closing a Reviewer §5 item that was cheap enough to
action rather than defer. Committed rather than left in a scratch script,
because Phase 2 Reviewer A could not check a fallback rate whose evidence lived
outside the repository (`phase2_summary.md` §11) and that gap is not repeated.

  A  Reviewer B §5.1 — the S=0.0016 depletion at 20,000 seeds, not 600.
  B  Reviewer B §5.2.1 — a LEARNED upper bound on predicting the station count
     from question text alone. B scored rules it could think of; a decision tree
     searches the space. If this clears ~85% the lookup verdict flips.
  C  Reviewer B §5.2.3 — is S=0.0016 the only tie-reachable slope on a finer grid?
  D  Reviewer B §5.2.4 — does the secant/true-root disagreement concentrate in
     the low-update instances (a stopping artefact) or spread (a rounding one)?
  E  Reviewer B §5.3 — the `_froude_capped_slope` exactness pathology is a
     property of the HELPER, so every template in that file which divides by
     round(sqrt(S), 5) inherits it. Enumerate them.
"""
from __future__ import annotations

import math
import random
import sys
from collections import Counter
from fractions import Fraction

sys.path.insert(0, __file__.rsplit('tests', 1)[0])

from tests.template_integrity.core import discover, generate  # noqa: E402

SEEDS = int(sys.argv[1]) if len(sys.argv) > 1 else 20000


def _refs():
    return {r.template_id: r for r in discover()}


def measurement_a(refs, n=SEEDS):
    """Slope depletion at S=0.0016, before vs after, at scale."""
    import data.templates.branches.civil_engineering.water_resources.uniform_flow as uf
    kept = Counter()
    rejected_slopes = Counter()
    for s in range(n):
        random.seed(s)
        q, _ = uf.template_normal_depth_iteration()
        val = q.split('slope of S = ')[1].split('.')[0] + '.' + \
            q.split('slope of S = ')[1].split('.')[1].split(' ')[0].rstrip('.')
        kept[float(val)] += 1
    share = 100 * kept[0.0016] / n
    print(f'A  S=0.0016 share of the ACCEPTED pool at {n} seeds: {share:.3f}%')
    print(f'   (Reviewer B measured 3.50% at 600 seeds; master was 4.83%)')
    return share


def measurement_b(refs, n=8000):
    """A learned upper bound on question-text-only prediction of n."""
    from sklearn.tree import DecisionTreeClassifier
    from sklearn.model_selection import cross_val_score
    import numpy as np

    ref = refs['template_line_balancing_heuristic']
    X, y = [], []
    for s in range(n):
        inst = generate(ref, s, capture=True)
        if not inst.ok:
            continue
        v = inst.values
        ts = sorted(int(v[f'ts[{i}]']) for i in range(5))
        CT, total, n_st = int(v['CT']), int(v['total']), int(v['n'])
        nmin = math.ceil(total / CT)
        feats = ts + [CT, total, nmin, total / CT,
                      sum(1 for t in ts if 2 * t > CT),
                      sum(1 for i in range(5) for j in range(i + 1, 5)
                          if ts[i] + ts[j] > CT),
                      max(ts) / CT, min(ts) / CT]
        X.append(feats)
        y.append(n_st)
    X, y = np.array(X), np.array(y)
    floor = 100 * max(Counter(y).values()) / len(y)
    print(f'B  learned upper bound on predicting n from question text, {len(y)} instances')
    print(f'   blind-constant floor                     : {floor:.2f}%')
    for depth in (2, 3, 5, None):
        clf = DecisionTreeClassifier(max_depth=depth, random_state=0)
        acc = 100 * cross_val_score(clf, X, y, cv=5).mean()
        label = f'depth {depth}' if depth else 'unbounded'
        print(f'   decision tree, {label:<10s}                : {acc:.2f}%')
    print(f'   (Reviewer B\'s best hand-built rule: 64.97%; its stated flip threshold ~85%)')


def measurement_c(refs):
    """Is S=0.0016 the only tie-reachable slope? Sweep a 5-dp grid."""
    reachable = []
    for i in range(80, 301):
        S = round(i / 100000 * 10, 5)          # the 4-dp grid the sampler uses
        sq = round(math.sqrt(S), 5)
        if sq * sq == S or Fraction(str(sq)) ** 2 == Fraction(str(S)):
            reachable.append((S, sq))
    print('C  slopes in [0.0008, 0.0030] whose sqrt is EXACT at 5 dp:')
    for S, sq in reachable:
        # a tie needs K = Q*n/sq to terminate, i.e. 1/sq must be a
        # terminating decimal: sq must be 2^a*5^b times a power of ten
        f = Fraction(str(sq))
        d = f.denominator
        num = f.numerator
        term = all(p in (2, 5) for p in _prime_factors(num))
        print(f'     S={S:<8} sqrt={sq:<8} 1/sqrt terminates: {term}')
    print('   (only a terminating 1/sqrt makes K terminating, hence tie-reachable)')


def _prime_factors(n):
    out = set()
    d = 2
    while d * d <= n:
        while n % d == 0:
            out.add(d)
            n //= d
        d += 1
    if n > 1:
        out.add(n)
    return out


def measurement_d(refs, n=4000):
    """Does the secant/true-root disagreement concentrate in low-update runs?"""
    import re
    ref = refs['template_normal_depth_iteration']
    by_updates = Counter()
    diff_by_updates = Counter()
    for s in range(n):
        inst = generate(ref, s, capture=True)
        if not inst.ok:
            continue
        v = inst.values
        b, z, K = v['b'], v['z'], v['K']
        rect = (z == 0.0)

        def sf(y):
            A = b * y if rect else (b + z * y) * y
            P = b + 2 * y if rect else b + 2 * y * math.sqrt(1 + z * z)
            return A * (A / P) ** (2.0 / 3.0)

        lo, hi = 0.2, 6.0
        for _ in range(200):
            mid = (lo + hi) / 2
            if sf(mid) < K:
                lo = mid
            else:
                hi = mid
        direct = round((lo + hi) / 2, 3)
        traced = float(re.search(r'normal depth is ([\d.]+) m', inst.solution).group(1))
        u = len(re.findall(r'Update \d+:', inst.solution))
        by_updates[u] += 1
        if abs(direct - traced) > 1e-9:
            diff_by_updates[u] += 1
    print(f'D  secant-vs-true-root disagreement by update count, {n} seeds:')
    for u in sorted(by_updates):
        tot, bad = by_updates[u], diff_by_updates[u]
        print(f'     {u} update(s): {bad:5d} / {tot:5d}  ({100*bad/tot:5.2f}%)')
    tot, bad = sum(by_updates.values()), sum(diff_by_updates.values())
    print(f'     overall    : {bad:5d} / {tot:5d}  ({100*bad/tot:5.2f}%)')


def measurement_e():
    """Which templates inherit the _froude_capped_slope exactness pathology?"""
    import ast
    import os
    path = os.path.join('data', 'templates', 'branches', 'civil_engineering',
                        'water_resources', 'uniform_flow.py')
    src = open(path, encoding='utf-8').read()
    tree = ast.parse(src)
    print('E  templates in uniform_flow.py that call _froude_capped_slope AND '
          'divide by round(sqrt(S), 5):')
    for node in tree.body:
        if not isinstance(node, ast.FunctionDef) or not node.name.startswith('template_'):
            continue
        body = ast.get_source_segment(src, node) or ''
        calls_helper = '_froude_capped_slope' in body
        rounds_sqrt = 'round(math.sqrt(S), 5)' in body
        divides = '/ sqS' in body or '/ {sqS' in body or 'sqS' in body
        if calls_helper and rounds_sqrt:
            flag = 'DIVIDES BY sqS' if divides else 'multiplies only'
            print(f'     {node.name:<48s} {flag}')
    print('   A template that only MULTIPLIES by sqS cannot form a terminating')
    print('   quotient this way, so the pathology is confined to the dividers.')


if __name__ == '__main__':
    refs = _refs()
    measurement_a(refs)
    print()
    measurement_b(refs)
    print()
    measurement_c(refs)
    print()
    measurement_d(refs)
    print()
    measurement_e()
