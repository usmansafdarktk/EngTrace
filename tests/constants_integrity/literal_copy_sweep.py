"""C3.8 - literal copies of constants-table leaves inside template code.

    python -m tests.constants_integrity.literal_copy_sweep

A template that writes a table's value as a literal instead of reading the table
does not follow a correction to that table (C1 Reviewer G, G-1: SCS_IA_RATIO's 0.2).
This lists CANDIDATES for a human to triage into one of three outcomes (spec C3.8):
a declared `@copied-in`, a template changed to read the table (a P6 event with a
before/after dump), or a coincidence. It gates nothing; the triage record is
docs/re-implementation-sep/phaseC3_literal_copies.md.

THE PREDICATE, and how it got here. A leaf is a candidate against a template literal
when all of these hold:
  * the leaf is a non-integer written with a decimal point and >= 3 significant
    figures as written ('0.0188', '350.0'),
  * the template module imports its branch's constants module,
  * the literal sits inside a template_* function, and
  * its source token is IDENTICAL to the table's written token.
The first design matched every leaf VALUE against every literal and listed 2,147 hits,
dominated by integers (temp_C 25 against any 25). The second added the 3-s.f. and
import conditions and still listed 590, matched by value ('5000.0' against 5000,
helpers and main() included). This third lists 36 (table, template) pairs - a size a
person can read line by line, which is what triage needs.

KNOWN MISS, stated: a copy below 3 significant figures is invisible to this rule. The
case that motivated C3.8, SCS_IA_RATIO's 0.2, is one; it is held by its @copied-in
declaration, which census P-COPY verifies instead.
"""
from __future__ import annotations

import ast
import os
import sys
from collections import defaultdict

REPO = os.path.abspath(os.path.join(os.path.dirname(__file__), '..', '..'))
if REPO not in sys.path:
    sys.path.insert(0, REPO)

from tests.constants_integrity.census import (  # noqa: E402
    BRANCHES, BRANCHES_DIR, _leaves, header_fields, numeric_tables, parse_copies)

MIN_SF = 3


def _sig(tok):
    """Significant figures as WRITTEN: '0.0188' -> 3, '350.0' -> 4, '25' -> 2."""
    t = tok.lower().lstrip('+-').split('e')[0].replace('_', '')
    digits = t.replace('.', '').lstrip('0')
    return len(digits) if '.' in t else len(digits.rstrip('0') or '0')


def _literal_tokens(src):
    """(table, leaf path) -> the token the table writes for that leaf."""
    out = {}
    for node in ast.parse(src).body:
        if not (isinstance(node, ast.Assign) and isinstance(node.targets[0], ast.Name)):
            continue
        name = node.targets[0].id

        def walk(v, path):
            if isinstance(v, ast.Constant) and isinstance(v.value, (int, float)) \
                    and not isinstance(v.value, bool):
                out[(name, path)] = ast.get_source_segment(src, v)
            elif isinstance(v, ast.Dict):
                for k, x in zip(v.keys, v.values):
                    if isinstance(k, ast.Constant):
                        walk(x, path + (('d', k.value),))
            elif isinstance(v, (ast.List, ast.Tuple)):
                for i, x in enumerate(v.elts):
                    walk(x, path + (('l' if isinstance(v, ast.List) else 't', i),))
        walk(node.value, ())
    return out


def _template_literals(branch):
    """[(rel path, line, token, value, owning top-level function)] outside constants.py."""
    out = []
    for dirpath, _dirs, files in os.walk(os.path.join(BRANCHES_DIR, branch)):
        for fn in files:
            if not fn.endswith('.py') or fn in ('constants.py', '__init__.py'):
                continue
            path = os.path.join(dirpath, fn)
            src = open(path, encoding='utf-8').read()
            if f'branches.{branch}.constants' not in src:
                continue
            tree = ast.parse(src)
            rel = os.path.relpath(path, REPO).replace(os.sep, '/')
            for top in tree.body:
                if not (isinstance(top, ast.FunctionDef) and top.name.startswith('template_')):
                    continue
                for node in ast.walk(top):
                    if isinstance(node, ast.Constant) and isinstance(node.value, (int, float)) \
                            and not isinstance(node.value, bool):
                        tok = ast.get_source_segment(src, node) or repr(node.value)
                        out.append((rel, node.lineno, tok, node.value, top.name))
    return out


def sweep():
    """{(branch, table, template, rel): [(declared?, leaf key, token, line)]}."""
    pairs = defaultdict(list)
    for branch in BRANCHES:
        src = open(os.path.join(BRANCHES_DIR, branch, 'constants.py'), encoding='utf-8').read()
        toks = _literal_tokens(src)
        ns = {}
        exec(src, ns)                                            # noqa: S102
        by_token = defaultdict(list)
        for rel, line, tok, _val, fn in _template_literals(branch):
            by_token[tok].append((rel, line, fn))
        for t in numeric_tables(src):
            declared = {(tid, lit) for tid, lit in parse_copies(header_fields(t['header']).get('copied-in'))}
            for path, v in _leaves(ns.get(t['name'])):
                tok = toks.get((t['name'], path))
                if tok is None or isinstance(v, int) or '.' not in tok or _sig(tok) < MIN_SF:
                    continue
                for rel, line, fn in by_token.get(tok, []):
                    key = ''.join(f'[{k!r}]' for _kind, k in path)
                    pairs[(branch, t['name'], fn, rel)].append(((fn, float(v)) in declared, key, tok, line))
    return pairs


def main():
    pairs = sweep()
    for (branch, table, fn, rel), hits in sorted(pairs.items()):
        tag = 'DECLARED' if all(h[0] for h in hits) else 'candidate'
        print(f'{tag:9s} {branch[:5]}.{table} -> {fn} ({rel}): '
              + '; '.join(f'{k} = {tok} at line {ln}' for _d, k, tok, ln in hits))
    print(f'{len(pairs)} (table, template) pairs, {sum(len(h) for h in pairs.values())} literals')
    return 0


if __name__ == '__main__':
    sys.exit(main())
