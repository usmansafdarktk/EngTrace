"""Read E, Ec, G, mu and density out of a MIL-HDBK-5J design-property table page.

    from tests.constants_integrity.milhdbk import design_values
    design_values(page_text) -> {'elastic': {'E': ['9.9'], 'G': ['3.8'], ...},
                                 'density': '0.098', 'layout': 'block', 'columns': 1}

MIL-HDBK-5J (2003, Distribution Statement A) closes each "Design Mechanical and
Physical Properties" table with the elastic constants, in 10^3 ksi, and the density,
in lb/in^3. Its PDF text layer loses the table's geometry, and it does so in two
ways, both measured on the on-disk file rather than assumed:

  INTERLEAVED - each label is followed by its own value, digits often spaced out:
      "E, 103 k s i ....... 6 . 5 Ec, 103 k s i ...... 6 . 5 G, 103 k s i ... 2 . 4 µ ... 0 . 3 5"
  BLOCK - every label first, then every value, one table column after another:
      "E, 103 k s i ... Ec, 103 k s i ... G, 103 k s i ... µ ... 10.5 10.7 4.0 0.33 10.7 10.9 4.0 0.33"

The labels are separated by leader dots (five or more); a missing value is printed
"...", which is kept as None. A page whose segment fits neither layout, whose value
count is not a multiple of its label count, or whose E/Ec/G label does not say
10^3 ksi, raises ValueError: a value is never guessed from a layout it does not fit.

Values are returned as the STRINGS the page prints, so a caller compares digits.
"""
from __future__ import annotations

import re

LEADER = re.compile(r'\.{5,}')
VALUE = re.compile(r'(?<![\w.])(\d+\.\d+|\.\.\.)(?![\w.])')


def normalize(seg):
    s = re.sub(r'\s*\n\s*', ' ', seg)
    # a run of single digits and dots separated by single spaces is one number
    # spaced out by the text layer: "2 9 . 0" -> "29.0", ". . ." -> "..."
    s = re.sub(r'(?<!\S)[\d.](?: [\d.])+(?!\S)', lambda m: m.group(0).replace(' ', ''), s)
    s = re.sub(r'(?<!\S)(\d+\.) (\d+)(?!\S)', r'\1\2', s)       # "0 . 0639" -> "0.0639"
    return re.sub(r' {2,}', ' ', s)


def _key(label, prev_main):
    s = label.replace(' ', '')
    if s.startswith(('µ', 'μ')):
        return 'mu', prev_main
    if ',' in s:
        main = s.split(',', 1)[0]
        if main in ('E', 'Ec', 'G') and '103ksi' not in s:
            raise ValueError(f'label {label!r} does not state 10^3 ksi')
        sub = s.split(':', 1)[1] if ':' in s else ''
        return main + (f':{sub}' if sub else ''), main
    if prev_main and s:
        return f'{prev_main}:{s}', prev_main
    raise ValueError(f'cannot name the label {label!r}')


def design_values(text):
    j = text.find('Physical Properties:')
    i = text.rfind('E, 10', 0, j if j >= 0 else len(text))
    if i < 0 or j < 0:
        raise ValueError('no "E, 10 ... Physical Properties:" segment on this page')
    parts = LEADER.split(normalize(text[i:j]))
    n = len(parts) - 1
    if n < 1:
        raise ValueError('no leader dots in the elastic-constant segment')
    vals = [[m.group(1) for m in VALUE.finditer(p)] for p in parts]
    if vals[0]:
        raise ValueError(f'a value before the first label: {vals[0]}')
    if all(len(v) == 1 for v in vals[1:]):
        layout, cols = 'interleaved', 1
        labels = [parts[0]] + [VALUE.sub('', p) for p in parts[1:-1]]
        table = [[v[0]] for v in vals[1:]]
    elif not any(vals[1:-1]) and vals[-1] and len(vals[-1]) % n == 0:
        layout, cols = 'block', len(vals[-1]) // n
        labels = parts[:-1]
        table = [[vals[-1][c * n + k] for c in range(cols)] for k in range(n)]
    else:
        raise ValueError(f'{n} labels and values {[len(v) for v in vals[1:]]} fit neither layout')
    out, prev = {}, None
    for label, row in zip(labels, table):
        key, prev = _key(label.strip(), prev)
        if key in out:
            raise ValueError(f'label {key!r} occurs twice')
        out[key] = [None if v == '...' else v for v in row]
    rest = normalize(text[j:j + 400])
    k = rest.find('lb/in')
    dens = VALUE.search(rest, k) if k >= 0 else None
    return dict(elastic=out, density=dens.group(1) if dens else None, layout=layout, columns=cols)
