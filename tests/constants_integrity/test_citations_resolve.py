"""C1.6 (was C2.5) - every provenance tag resolves to the artefact it names.

    python -m tests.constants_integrity.test_citations_resolve
    python -m tests.constants_integrity.test_citations_resolve --selftest
    python -m tests.constants_integrity.test_citations_resolve --verbose

The C2.2 plausibility suite asks *"is this value right?"*. This one asks the
prior question: *"does the thing the comment points at actually exist, and is
it what the comment says it is?"* Phase C2 shipped a green 215-check suite that
nonetheless certified three heats of formation with no on-disk backing at all,
because the reference values lived in the test file rather than in the artefact
the citations named (Reviewer G, findings G-1, G-2, G-12; DECISIONS D-034).

Phase C2 wrote four properties for two tables. Reviewer G called generalising
them "blocking for C3" (G-12), and Phase C1.6 does it: every tagged table in
every branch, and every artefact type under docs/references/. C2's four are kept
exactly, for the tables they were written for:

  P1  every value-bearing row of CP_PARAMS / HEATS_OF_FORMATION carries a tag
  P2  the CAS NUMBER WRITTEN IN THE TAG resolves - not merely the species key
  P3  the tag's claim class matches the evidence
  P4  no reference value is hardcoded in a test file (D-034)

and seven are added, one per clause of the C1.2 vocabulary (spec §C1.2):

  R1  every tag parses under C1.2. A DEPRECATED tag - [VERIFY: X], [REALISM],
      [DERIVABLE], [ON-DISK: xlsx], [ON-DISK visual], or an [ON-DISK] with no
      `artefact @ locator` - is LEGACY: counted and listed, not failed. The
      list is C3.7's worklist and its exit gate is zero.
  R2  an [ON-DISK] artefact exists under docs/references/, MANIFEST.json vouches
      for it (status acquired or present), and its bytes match the SHA-256 the
      manifest recorded. A file the manifest does not vouch for is not
      something a clone can re-acquire, whatever is written in it.
  R3  the locator resolves inside the artefact, by type: `quantity=` (CODATA
      text), `T=` + `col=` (NIST fluid TSV - the row must be ON the grid, which
      is the clamp lesson of the acquisition), `cas=` (NIST WebBook JSON),
      `page=` [+ `text=`] (PDF), `member=` (zip), `text=` (HTML / text), and
      `member=` + `wavelength=<x>nm|um` (the refractiveindex.info archive: the
      index EVALUATED at that wavelength, refused outside the dataset's range;
      tests/constants_integrity/refractiveindex.py).
  R4  the RELATION holds. `precision=exact|<n>sf|<n>dp`: the constant equals the
      artefact value at that precision, exactly. `tol=<x>%`: it lies within.
      Neither: the value is not machine-compared, and the citation is counted
      LOCATOR-ONLY rather than passed as if it had been. `via="NAME/x"` says the
      table stores NAME/x rather than x (MEDIA_VELOCITIES stores C0/n), and the
      relation is checked on x = NAME/constant. `scale=<f>` says the table stores
      the artefact's quantity times f (ATMOSPHERIC_PRESSURE_KPA: kPa against
      CODATA's Pa, scale=1e-3), and the relation is checked in the artefact's
      unit, so `precision=` counts the artefact's digits. A value interpolated between two
      tabulated rows must satisfy the relation at BOTH rows too, so a verdict
      never rests on the interpolation.
  R5  an [ON-DISK:LOCAL-ONLY] path is listed in MANIFEST.json's
      local_only_copyrighted. Present on this machine: its page exists. Absent:
      UNRESOLVABLE-FROM-CLONE, counted - never passed silently, never failed.
  R6  [DERIVED], [BY-DEFINITION], [KNOWN-DEFECTIVE] and [UNVERIFIED] state what
      they must: the derivation, the definition, the measured error, the reason.
  R7  the tag agrees with the table's declared @kind: [POLICY: sampling-only]
      only on a `range`, and a `property`, `standard` or `measured-constant` -
      a claim about the world - never [POLICY] or [BY-DEFINITION].

Every count printed at the end is one of these outcomes; nothing is summed
across them.
"""
from __future__ import annotations

import hashlib
import json
import math
import os
import re
import sys
import zipfile

REPO = os.path.abspath(os.path.join(os.path.dirname(__file__), '..', '..'))
if REPO not in sys.path:
    sys.path.insert(0, REPO)

from tests.constants_integrity.census import (  # noqa: E402
    BRANCHES, BRANCHES_DIR, header_fields, numeric_tables)

REFS = os.path.join(REPO, 'docs', 'references')
CONSTANTS = os.path.join(BRANCHES_DIR, 'chemical_engineering', 'constants.py')
REF_PATH = os.path.join(REFS, 'nist_webbook', 'shomate_coefficients.json')

# ==========================================================================
# C2's P1-P3, unchanged in behaviour, for the two tables they were written for
# ==========================================================================

CITED_TABLES = ('CP_PARAMS', 'HEATS_OF_FORMATION')
TAG_RE = re.compile(r'#\s*\[(ON-DISK|DERIVED|KNOWN-DEFECTIVE|BY-DEFINITION)\]\s*(.*)')
CAS_RE = re.compile(r'\b(\d{2,7}-\d{2}-\d)\b')
ROW_RE = re.compile(r'^\s*"([^"]+)"\s*:')


def parse_table(src, name):
    """Yield (species, (tag_class, tag_text)) for each row of a dict literal."""
    lines = src.splitlines()
    start = next(i for i, l in enumerate(lines) if l.startswith(name + ' = {'))
    depth, pending = 0, []
    for line in lines[start:]:
        depth += line.count('{') - line.count('}')
        m = TAG_RE.search(line)
        if m:
            pending.append((m.group(1), m.group(2).strip()))
            continue
        if line.lstrip().startswith('#'):
            if pending:
                pending[-1] = (pending[-1][0],
                               pending[-1][1] + ' ' + line.lstrip('# ').strip())
            continue
        row = ROW_RE.match(line)
        if row:
            yield row.group(1), (pending[-1] if pending else (None, ''))
            pending = []
        if depth == 0 and line.startswith('}'):
            break


def c2_checks(src, ref):
    by_cas = {}
    for sp, e in ref.items():
        if not sp.startswith('_') and e.get('cas'):
            by_cas.setdefault(e['cas'], []).append(sp)
    failures, checks = [], 0
    for table in CITED_TABLES:
        for sp, (cls, text) in parse_table(src, table):
            checks += 1
            if cls is None:
                failures.append(f'P1 {table}[{sp}] carries no provenance tag')
                continue
            if cls in ('DERIVED', 'BY-DEFINITION'):
                if not text:
                    failures.append(f'P3 {table}[{sp}] is [{cls}] but states nothing')
                continue
            if cls == 'KNOWN-DEFECTIVE':
                continue
            entry = ref.get(sp)
            if entry is None:
                failures.append(f'P3 {table}[{sp}] is tagged [ON-DISK] but no entry for '
                                f'it exists in {os.path.basename(REF_PATH)}')
                continue
            cas = CAS_RE.search(text)
            if not cas:
                failures.append(f'P2 {table}[{sp}] is [ON-DISK] but names no CAS number')
                continue
            cas = cas.group(1)
            if cas not in by_cas:
                failures.append(f'P2 {table}[{sp}] cites CAS {cas}, which appears nowhere '
                                f'in {os.path.basename(REF_PATH)}')
            elif entry.get('cas') != cas:
                failures.append(f'P2 {table}[{sp}] cites CAS {cas}, but the artefact entry '
                                f'for {sp} is CAS {entry.get("cas")} '
                                f'(that CAS is {", ".join(by_cas[cas])})')
    return checks, failures


def p4_check():
    embedded = []
    test_dir = os.path.dirname(os.path.abspath(__file__))
    for fn in sorted(os.listdir(test_dir)):
        if not fn.startswith('test_') or not fn.endswith('.py') \
                or fn == os.path.basename(__file__):
            continue
        body = open(os.path.join(test_dir, fn), encoding='utf-8').read()
        for m in re.finditer(r'^([A-Z][A-Z0-9_]*(?:REF|_VALUES|_TABLE))\s*=\s*\{', body, re.M):
            embedded.append(f'{fn}:{m.group(1)}')
    if embedded:
        return ['P4 reference values are hardcoded in a test file, so the suite checks '
                'itself rather than the artefact: ' + ', '.join(embedded)]
    return []


# ==========================================================================
# C1.2 tags: extraction and attachment
# ==========================================================================

OPEN_RE = re.compile(r'\[(ON-DISK|DERIVED|BY-DEFINITION|KNOWN-DEFECTIVE|VERIFY|POLICY|'
                     r'REALISM|DERIVABLE|UNVERIFIED)\b')
KV_RE = re.compile(r'([A-Za-z_]+)=("([^"]*)"|\S+)')


def _comment(line):
    """The comment part of a source line, or None. Strings containing '#' are rare
    in constants.py; a '#' inside a quoted row key is skipped by requiring the
    comment to start outside the last string literal on the line."""
    s = line.lstrip()
    if s.startswith('#'):
        return s[1:].strip()
    i = line.rfind('#')
    if i > 0 and line[:i].count('"') % 2 == 0 and line[:i].count("'") % 2 == 0:
        return line[i + 1:].strip()
    return None


def extract_tags(src):
    """Every tag in a constants.py source, attached to (table, row or None).

    A tag is attached to a ROW when it sits in the comment block directly above
    a `"key":` line inside the table, or trails that line; otherwise to the
    TABLE when it sits in the table's header block or body. A bracket may close
    on a later comment line (civil's CV_RANGES_M2_YR). Tags outside every table
    and header - module docstrings describing the vocabulary - are ignored.
    """
    lines = src.split('\n')
    tables = numeric_tables(src)
    owner = {}                                  # line index -> table name
    for t in tables:
        i = t['lineno'] - 2
        while i >= 0 and (not lines[i].strip() or lines[i].lstrip().startswith('#')):
            owner[i] = t['name']
            i -= 1
        for j in range(t['lineno'] - 1, t['end_lineno']):
            owner[j] = t['name']
    out = []
    pending = []                                # tags awaiting a row line
    for idx, line in enumerate(lines):
        name = owner.get(idx)
        row = ROW_RE.match(line)
        com = _comment(line)
        found = []
        if com is not None and name is not None:
            pos = 0
            while True:
                m = OPEN_RE.search(com, pos)
                if not m:
                    break
                inner, end = com[m.start() + 1:], None
                close = inner.find(']')
                extra, k = '', idx
                while close < 0 and k + 1 < len(lines) and _comment(lines[k + 1]) is not None:
                    k += 1
                    inner += ' ' + _comment(lines[k])
                    close = inner.find(']')
                # A tag STARTS its comment, or follows only other tags on it.
                # Anything else is a mid-line mention - legacy syntax when a
                # table has nothing better, prose when it does (see below).
                lead = re.sub(r'\[[^\]]*\]', '', com[:m.start()]).strip()
                if close < 0:
                    found.append(dict(raw=inner, bracket=inner, body='', line=idx + 1,
                                      unclosed=True, line_start=not lead))
                    break
                bracket = inner[:close]
                rest = inner[close + 1:]
                nxt = OPEN_RE.search(rest)
                body = rest[:nxt.start()] if nxt else rest
                found.append(dict(raw='[' + bracket + ']' + body, bracket=bracket,
                                  body=body.strip(), line=idx + 1, unclosed=False,
                                  line_start=not lead))
                pos = m.start() + 1 + close + 1 if k == idx else len(com)
        # A tag waiting for a row belongs to ITS table: when that table ends
        # first, the tag was about the table. The first version never flushed at
        # a table boundary, so a tag trailing the scalar WATER_DENSITY_KG_M3 was
        # reported against row 'ASTM A992' of STEEL_FY_KSI, twenty lines later.
        if pending and name != pending[0]['table']:
            for f in pending:
                f['row'] = None
                out.append(f)
            pending = []
        body_line = name is not None and idx >= _table(tables, name)['lineno'] - 1
        pure_comment = line.lstrip().startswith('#')
        if row and body_line:
            for f in pending:
                f['row'] = row.group(1)
                out.append(f)
            pending = []
        for f in found:
            f['table'] = name
            if row and body_line:
                f['row'] = row.group(1)            # trails its own row
                out.append(f)
            elif body_line and pure_comment:
                pending.append(f)                  # sits above the next row
            else:
                f['row'] = None                    # header, or trails a non-row line
                out.append(f)
        if body_line and not row and not pure_comment and line.strip() and pending:
            for f in pending:                      # e.g. a comment just above `}`
                f['row'] = None
                out.append(f)
            pending = []
    for f in pending:
        f['row'] = None
        out.append(f)
    # A mid-line bracket is a legacy tag in a table that has nothing better, and
    # PROSE in a table that carries line-start tags ("each row carries an
    # [ON-DISK] citation"). The first version counted CP_PARAMS' header prose as
    # two LEGACY tags; the flag was added and this filter forgotten, once.
    tagged = {f['table'] for f in out if f['line_start']}
    return [f for f in out if f['line_start'] or f['table'] not in tagged]


def _table(tables, name):
    return next(t for t in tables if t['name'] == name)


def classify_tag(tag):
    """(class, detail) under C1.2, or ('LEGACY', why)."""
    b, body = tag['bracket'].strip(), tag['body']
    if tag.get('unclosed'):
        return 'MALFORMED', 'the bracket never closes'
    if b == 'ON-DISK':
        if ' @ ' in f' {body} ' or body.startswith('@'):
            return 'ON-DISK', ''
        if CAS_RE.search(body) and body.startswith('NIST'):
            return 'ON-DISK-C2', ''
        return 'LEGACY', 'an [ON-DISK] with no `artefact @ locator`'
    if b == 'ON-DISK:LOCAL-ONLY':
        return ('ON-DISK:LOCAL-ONLY', '') if ' @ ' in f' {body} ' else \
            ('MALFORMED', 'local-only citation with no `path @ locator`')
    if b.startswith('ON-DISK'):
        return 'LEGACY', f'[{b}]'
    if b.startswith('POLICY'):
        return ('POLICY', '') if 'sampling-only' in b else ('LEGACY', f'[{b}]')
    if b in ('DERIVED', 'BY-DEFINITION', 'KNOWN-DEFECTIVE', 'UNVERIFIED'):
        return b, ''
    if b.startswith(('VERIFY', 'REALISM', 'DERIVABLE')):
        return 'LEGACY', f'[{b}]'
    if b.startswith(('DERIVED', 'BY-DEFINITION', 'KNOWN-DEFECTIVE', 'UNVERIFIED')):
        return b.split(':')[0], ''
    return 'MALFORMED', f'[{b}] is not a C1.2 class'


def parse_citation(body):
    """`artefact @ k=v k="v v" flag ...` -> (artefact, {k: v}, {flags}, note)."""
    head, _, tail = body.partition(' @ ')
    if body.startswith('@'):
        head, tail = '', body[1:]
    note = ''
    if ';' in tail:
        tail, note = tail.split(';', 1)
    kv, flags = {}, set()
    pos = 0
    for m in KV_RE.finditer(tail):
        kv[m.group(1)] = m.group(3) if m.group(3) is not None else m.group(2)
        for word in tail[pos:m.start()].split():
            flags.add(word)
        pos = m.end()
    flags.update(tail[pos:].split())
    return head.strip(), kv, flags, note.strip()


# ==========================================================================
# Artefact readers: one per locator type. Each returns (value_or_None, error).
# ==========================================================================

def _codata(path, quantity):
    for ln in open(path, encoding='utf-8'):
        cols = re.split(r'\s{2,}', ln.strip())
        if cols and cols[0] == quantity:
            txt = cols[1].replace(' ', '').replace('...', '')
            try:
                return float(txt), None
            except ValueError:
                return None, f'value {cols[1]!r} for {quantity!r} does not parse'
    return None, f'quantity {quantity!r} is not in {os.path.basename(path)}'


def _tsv(path, T, col):
    lines = open(path, encoding='utf-8').read().splitlines()
    header = lines[0].split('\t')
    if col not in header:
        return None, f'column {col!r} not in the header of {os.path.basename(path)}'
    j = header.index(col)
    for ln in lines[1:]:
        cells = ln.split('\t')
        try:
            if abs(float(cells[0]) - float(T)) < 0.005:
                return float(cells[j]), None
        except (ValueError, IndexError):
            continue
    return None, (f'no row at T={T} K in {os.path.basename(path)} - NOT on the grid, '
                  f'and the nearest row is not the row cited')


def _pdf_page(path, page, text=None):
    try:
        from pypdf import PdfReader
    except ImportError:
        return None, 'pypdf unavailable: PDF locators cannot be checked on this machine'
    reader = PdfReader(path)
    n = len(reader.pages)
    if not 1 <= int(page) <= n:
        return None, f'page {page} does not exist ({n} pages)'
    if text:
        got = re.sub(r'\s+', ' ', reader.pages[int(page) - 1].extract_text() or '')
        if re.sub(r'\s+', ' ', text) not in got:
            return None, f'text {text!r} is not on page {page}'
    return 'page', None


def _zip_member(path, member):
    with zipfile.ZipFile(path) as z:
        names = z.namelist()
    if not any(n == member or n.endswith('/' + member) for n in names):
        return None, f'member {member!r} is not in {os.path.basename(path)}'
    return 'member', None


def _index_at(path, member, wavelength):
    """(n, brackets or None, error) - a refractive index evaluated at a wavelength."""
    from tests.constants_integrity.refractiveindex import index_at
    m = re.fullmatch(r'(\d+(?:\.\d+)?)(nm|um)', wavelength)
    if not m:
        return None, None, f'wavelength={wavelength!r} is not <number>nm or <number>um'
    lam = float(m.group(1)) * (1e-3 if m.group(2) == 'nm' else 1.0)
    try:
        value, _how, _cond, brackets = index_at(member, lam, archive=path)
    except (KeyError, ValueError) as exc:
        return None, None, f'{member}: {exc}'
    return value, brackets, None


def _solve_via(ns, via, const):
    """The cited quantity x, from a table that stores NAME/x (or x itself)."""
    if via == 'x':
        return const, None
    m = re.fullmatch(r'([A-Za-z_]\w*)/x', via)
    if not m:
        return None, f'via={via!r} is not a form this resolver solves ("x" or "NAME/x")'
    num = ns.get(m.group(1))
    if isinstance(num, bool) or not isinstance(num, (int, float)):
        return None, f'via={via!r} names {m.group(1)!r}, which the module does not define as a number'
    if const == 0:
        return None, f'via={via!r} cannot be solved for a stored value of 0'
    return num / const, None


def _text_snippet(path, text):
    body = open(path, encoding='utf-8', errors='replace').read()
    body = re.sub(r'<[^>]+>', ' ', body)
    if re.sub(r'\s+', ' ', text) not in re.sub(r'\s+', ' ', body):
        return None, f'text {text!r} is not in {os.path.basename(path)}'
    return 'text', None


def _sha256(path):
    h = hashlib.sha256()
    with open(path, 'rb') as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b''):
            h.update(chunk)
    return h.hexdigest()


def _rounded(a, precision):
    if precision == 'exact':
        return a
    m = re.fullmatch(r'(\d+)(sf|dp)', precision)
    if not m:
        raise ValueError(f'precision={precision!r} is not exact, <n>sf or <n>dp')
    n = int(m.group(1))
    if m.group(2) == 'dp':
        return round(a, n)
    return float(f'{a:.{n - 1}e}')


def _same(x, y):
    return x == y or math.isclose(x, y, rel_tol=1e-12, abs_tol=0.0)


def _constant(ns, table, row, kv):
    obj = ns.get(table)
    if row is not None:
        if not isinstance(obj, dict) or row not in obj:
            return None, f'row {row!r} is not in {table}'
        obj = obj[row]
    for key in re.findall(r"\[([^\]]+)\]|([A-Za-z_]\w*)", kv.get('field', '')) \
            if kv.get('field') else []:
        k = key[0] or key[1]
        try:
            obj = obj[int(k)] if k.isdigit() else obj[k.strip("'\"")]
        except (KeyError, IndexError, TypeError):
            return None, f'field {kv["field"]!r} does not exist in {table}'
    if isinstance(obj, bool) or not isinstance(obj, (int, float)):
        return None, (f'{table}{"[" + repr(row) + "]" if row else ""} is not a single '
                      f'number - name the value with field=')
    return obj, None


# ==========================================================================
# The generalised check
# ==========================================================================

def check_source(branch, src, refs=REFS, manifest=None, kinds=None):
    """R1-R7 over one constants.py source. Returns (counts, failures, legacy)."""
    manifest = manifest if manifest is not None else json.load(
        open(os.path.join(refs, 'MANIFEST.json'), encoding='utf-8'))
    by_dest = {e['dest'].replace('\\', '/'): e for e in manifest['sources'].values()
               if e.get('dest')}
    local_only = {f['path'] for f in manifest.get('local_only_copyrighted', {}).get('files', [])}
    ns = {}
    exec(compile(src, f'{branch}/constants.py', 'exec'), ns)        # noqa: S102
    if kinds is None:
        kinds = {t['name']: header_fields(t['header']).get('kind', '').split(' ')[0]
                 for t in numeric_tables(src, ns)}
    counts = dict(tags=0, resolved=0, locator_only=0, unresolvable_from_clone=0,
                  legacy=0, stated=0)
    failures, legacy = [], []
    hashed = {}

    for tag in extract_tags(src):
        counts['tags'] += 1
        where = f"{branch}.{tag['table']}" + (f"[{tag['row']!r}]" if tag.get('row') else '') \
            + f" (line {tag['line']})"
        cls, why = classify_tag(tag)
        kind = kinds.get(tag['table'], '')

        # R7 - the tag agrees with the declared @kind
        if cls == 'POLICY' and kind != 'range':
            failures.append(f'R7 {where}: [POLICY: sampling-only] on a table declared '
                            f'@kind {kind or "(none)"} - only a range can be sampling policy')
        # A TABLE of facts cannot be true by definition; a ROW of one can - an
        # element's heat of formation is 0 because of how the scale is defined,
        # in a table that is otherwise measurements. The first version of R7
        # applied this to rows too and failed HEATS_OF_FORMATION's O2/H2/N2 -
        # the very example spec §C1.2 gives for [BY-DEFINITION].
        if cls == 'BY-DEFINITION' and tag.get('row') is None \
                and kind in ('property', 'standard', 'measured-constant'):
            failures.append(f'R7 {where}: a table-level [BY-DEFINITION] on a {kind} - '
                            f'a table of claims about the world is not true by definition')

        if cls == 'MALFORMED':
            failures.append(f'R1 {where}: {why}: {tag["raw"][:80]!r}')
            continue
        if cls == 'LEGACY':
            counts['legacy'] += 1
            legacy.append(f'{where}: {why}')
            continue
        if cls == 'ON-DISK-C2':
            counts['resolved'] += 1         # resolved by c2_checks, P2
            continue
        if cls in ('DERIVED', 'BY-DEFINITION', 'KNOWN-DEFECTIVE', 'UNVERIFIED'):
            payload = tag['body'] or tag['bracket'].partition(':')[2].strip()
            if not payload:
                failures.append(f'R6 {where}: [{cls}] states no '
                                f'{ {"DERIVED": "derivation", "BY-DEFINITION": "definition", "KNOWN-DEFECTIVE": "measured error", "UNVERIFIED": "reason"}[cls] }')
            else:
                counts['stated'] += 1
            continue
        if cls == 'POLICY':
            counts['stated'] += 1
            continue

        artefact, kv, flags, _note = parse_citation(tag['body'])
        if cls == 'ON-DISK:LOCAL-ONLY':
            if artefact not in local_only:
                failures.append(f'R5 {where}: {artefact!r} is not listed in MANIFEST.json '
                                f'local_only_copyrighted')
                continue
            full = os.path.join(REPO, artefact)
            if not os.path.exists(full):
                counts['unresolvable_from_clone'] += 1
                continue
            if 'page' in kv:
                _v, err = _pdf_page(full, kv['page'], kv.get('text'))
                if err:
                    failures.append(f'R5 {where}: {err}')
                    continue
            counts['resolved'] += 1
            continue

        # ---- [ON-DISK] ----------------------------------------------------
        rel = artefact.replace('\\', '/')
        full = os.path.join(refs, rel)
        entry = by_dest.get(rel)
        if not os.path.exists(full):
            failures.append(f'R2 {where}: cites {rel!r}, which does not exist under '
                            f'docs/references/')
            continue
        if entry is None or entry.get('status') not in ('acquired', 'present'):
            failures.append(f'R2 {where}: {rel!r} exists but MANIFEST.json does not vouch '
                            f'for it, so a clone cannot re-acquire it')
            continue
        if rel not in hashed:
            hashed[rel] = _sha256(full)
        if hashed[rel] != entry.get('sha256'):
            failures.append(f'R2 {where}: {rel!r} does not match the SHA-256 MANIFEST.json '
                            f'recorded - the file changed after it was cited')
            continue

        # R3 - the locator resolves
        ext = rel.rsplit('.', 1)[-1].lower()
        value, err, brackets = None, None, None
        if 'member' in kv and 'wavelength' in kv and ext == 'zip':
            value, brackets, err = _index_at(full, kv['member'], kv['wavelength'])
        elif 'quantity' in kv:
            value, err = _codata(full, kv['quantity'])
        elif 'T' in kv and 'col' in kv:
            value, err = _tsv(full, kv['T'], kv['col'])
        elif 'cas' in kv:
            ref = json.load(open(full, encoding='utf-8'))
            if not any(isinstance(e, dict) and e.get('cas') == kv['cas'] for e in ref.values()):
                err = f'CAS {kv["cas"]} is not in {os.path.basename(full)}'
        elif 'page' in kv and ext == 'pdf':
            value, err = _pdf_page(full, kv['page'], kv.get('text'))
            value = None
        elif 'member' in kv and ext == 'zip':
            _v, err = _zip_member(full, kv['member'])
        elif 'text' in kv:
            _v, err = _text_snippet(full, kv['text'])
        else:
            err = f'no locator this resolver knows for a .{ext} artefact: {tag["body"]!r}'
        if err:
            failures.append(f'R3 {where}: {err}')
            continue

        # R4 - the relation
        if 'precision' not in kv and 'tol' not in kv:
            counts['locator_only'] += 1
            continue
        if value is None:
            failures.append(f'R4 {where}: states a relation, but a {ext} locator yields no '
                            f'machine-readable value to compare')
            continue
        const, err = _constant(ns, tag['table'], tag.get('row'), kv)
        if not err and 'via' in kv:
            const, err = _solve_via(ns, kv['via'], const)
        if not err and 'scale' in kv:
            # the table stores the artefact's quantity in another unit: compare in
            # the artefact's own unit, so precision= keeps counting ITS digits
            try:
                scale = float(kv['scale'])
                if scale <= 0:
                    raise ValueError
                const = const / scale
            except ValueError:
                err = f'scale={kv["scale"]!r} is not a positive number'
        if err:
            failures.append(f'R4 {where}: {err}')
            continue
        # A tabulated value between two rows: the relation must also hold at both
        # rows, or the verdict is the interpolation's rather than the artefact's.
        rows = [('the artefact', value)] + [(f'the tabulated row at {w}', v)
                                            for w, v in (brackets or ())]
        if 'precision' in kv:
            try:
                target = _rounded(value, kv['precision'])
            except ValueError as exc:
                failures.append(f'R4 {where}: {exc}')
                continue
            if not _same(const, target):
                failures.append(f'R4 {where}: the constant is {const!r}, but the artefact '
                                f'gives {value!r}, which at precision={kv["precision"]} is '
                                f'{target!r}')
                continue
            split = [(lbl, v) for lbl, v in rows[1:] if not _same(_rounded(v, kv['precision']), target)]
            if split:
                failures.append(f'R4 {where}: the interpolated {value!r} is {target!r} at '
                                f'precision={kv["precision"]}, but {split[0][0]} gives '
                                f'{split[0][1]!r} - the rounding depends on the interpolation')
                continue
        else:
            tol = float(kv['tol'].rstrip('%'))
            out = [(lbl, v) for lbl, v in rows if abs(const - v) > abs(v) * tol / 100]
            if out:
                lbl, v = out[0]
                failures.append(f'R4 {where}: the constant {const!r} is '
                                f'{100 * (const - v) / v:+.3f}% from {lbl} '
                                f'{v!r}, outside tol={kv["tol"]}')
                continue
        counts['resolved'] += 1
    return counts, failures, legacy


def run(verbose=False):
    failures = []
    src = open(CONSTANTS, encoding='utf-8').read()
    ref = json.load(open(REF_PATH, encoding='utf-8'))
    c2_n, c2_f = c2_checks(src, ref)
    failures += c2_f
    failures += p4_check()
    total = dict(tags=0, resolved=0, locator_only=0, unresolvable_from_clone=0,
                 legacy=0, stated=0)
    legacy_all = []
    for branch in BRANCHES:
        bsrc = open(os.path.join(BRANCHES_DIR, branch, 'constants.py'),
                    encoding='utf-8').read().replace('\r\n', '\n')
        counts, f, legacy = check_source(branch, bsrc)
        failures += f
        legacy_all += legacy
        for k, v in counts.items():
            total[k] += v
        print(f"  {branch:24s} tags {counts['tags']:3d}  resolved {counts['resolved']:3d}  "
              f"locator-only {counts['locator_only']:2d}  stated {counts['stated']:3d}  "
              f"local-only-absent {counts['unresolvable_from_clone']:2d}  "
              f"LEGACY {counts['legacy']:3d}")
    if verbose:
        for x in legacy_all:
            print('    LEGACY ' + x)
    for x in failures:
        print('  - ' + x)
    print(f"{c2_n + 1} C2 checks (P1-P4); {total['tags']} tags across {len(BRANCHES)} "
          f"branches: {total['resolved']} resolved, {total['locator_only']} locator-only, "
          f"{total['stated']} stated, {total['unresolvable_from_clone']} unresolvable from a "
          f"clone, {total['legacy']} LEGACY (C3.7 worklist)")
    print('all pass' if not failures else f'{len(failures)} FAILURES')
    return 1 if failures else 0


# ==========================================================================
# Self-test: planted defects of materially different surface form, written from
# R1-R7's definitions, against the REAL on-disk artefacts. A plant counts only
# if it produces a failure the clean fixture does not.
# ==========================================================================

_CLEAN = '''
# @kind: defined
# @units: m/s
# [ON-DISK] codata_2022/allascii.txt @ quantity="speed of light in vacuum" precision=exact
C0 = 299792458

# @kind: measured-constant
# @units: F/m
# [ON-DISK] codata_2022/allascii.txt @ quantity="vacuum electric permittivity" precision=4sf
EPSILON_0 = 8.854e-12

# @kind: property
# @units: [0]=kg/m^3
FLUIDS = {
    # [ON-DISK] nist_fluid_properties/water_C7732185_isobar_1atm.tsv @ T=293.15 col="Density (kg/m3)" field=[0] precision=5sf
    "Water": (998.21, 1.0e-3),
    # [ON-DISK] nist_fluid_properties/water_C7732185_isobar_1atm.tsv @ T=293.15 col="Density (kg/m3)" field=[0] tol=0.5%
    "Water again": (1000.0, 1.0e-3),
}

# @kind: range
# @units: Hz
# [POLICY: sampling-only]
FREQ = (50, 2000)

# @kind: property
# @units: 1
# [VERIFY: CRC Handbook]
LEGACY_ROW = 1.0

# @kind: property
# @units: 1
# [UNVERIFIED] no source on disk; candidate NASA TR R-132 table 2
GAP = 2.0

# @kind: property
# @units: m/s
MEDIA = {
    # [ON-DISK] refractiveindex_info/refractiveindex.info-database-main.zip @ member="database/data/main/H2O/nk/Daimon-20.0C.yml" wavelength=589nm via="C0/x" precision=4sf
    "Water": C0 / 1.333,
    # [ON-DISK] refractiveindex_info/refractiveindex.info-database-main.zip @ member="database/data/main/H2O/nk/Warren-2008.yml" wavelength=589nm via="C0/x" precision=3sf
    "Ice": C0 / 1.31,
    # [ON-DISK] refractiveindex_info/refractiveindex.info-database-main.zip @ member="database/data/organic/C6H6 - benzene/nk/Chang.yml" wavelength=0.589um via="C0/x" tol=0.1%
    "Benzene": C0 / 1.501,
}

# @kind: defined
# @units: kPa
# [ON-DISK] codata_2022/allascii.txt @ quantity="standard atmosphere" scale=1e-3 precision=exact
ATM_KPA = 101.325
'''

_PLANTS = [
    # A unit scale (C3.1, ATMOSPHERIC_PRESSURE_KPA: kPa against CODATA's Pa). Two
    # forms of a unit error: the factor inverted, and the factor left out.
    ('R4', 'a unit scale in the wrong direction', 'scale=1e-3 precision=exact',
     'scale=1e3 precision=exact'),
    ('R4', 'a unit scale left out', 'scale=1e-3 precision=exact', 'precision=exact'),
    # Wavelength-evaluated indices (C3.1, MEDIA_VELOCITIES). Written from what the
    # locator means: an index is a function of wavelength, valid over a stated
    # range; `tabulated n2` is a different quantity; C0/n is not n; and a value
    # between two tabulated rows is the interpolation's unless both rows agree.
    ('R3', 'wavelength outside the dataset\'s stated range',
     'Daimon-20.0C.yml" wavelength=589nm', 'Daimon-20.0C.yml" wavelength=100nm'),
    ('R3', 'a nonlinear-index (n2) dataset cited as the index',
     '"database/data/main/H2O/nk/Warren-2008.yml"', '"database/data/other/mixed gases/air/n2/Geints.yml"'),
    ('R4', 'stored C0/n disagrees with the index at its precision',
     '"Water": C0 / 1.333', '"Water": C0 / 1.334'),
    ('R4', 'a 3-s.f. rounding that rests on the interpolation (rows 1.48 and 1.47)',
     'main/H2O/nk/Warren-2008.yml" wavelength=589nm via="C0/x" precision=3sf\n    "Ice": C0 / 1.31',
     'organic/C3H8O3 - glycerol/nk/Birkhoff.yml" wavelength=589nm via="C0/x" precision=3sf\n    "Ice": C0 / 1.47'),
    ('R4', 'via names a constant the module does not define',
     'wavelength=0.589um via="C0/x"', 'wavelength=0.589um via="CO/x"'),
    # R2 - the file does not exist (the brief's first mandated plant)
    ('R2', 'cites a file that does not exist',
     'codata_2022/allascii.txt @ quantity="speed', 'codata_2018/allascii.txt @ quantity="speed'),
    # R4 - a real file whose content at the locator disagrees (the second)
    ('R4', 'real file, content disagrees at the locator', 'EPSILON_0 = 8.854e-12',
     'EPSILON_0 = 8.855e-12'),
    # R4 - same class, different surface: a TSV row, and a tolerance breach
    ('R4', 'tolerance breached on a TSV value', '"Water again": (1000.0,',
     '"Water again": (1010.0,'),
    # R3 - the locator names nothing in a real file
    ('R3', 'quantity name not in the file', 'quantity="vacuum electric permittivity"',
     'quantity="electric permittivity of vacuum"'),
    # R3 - off-grid TSV row: near a real row, not on it (the NIST clamp lesson)
    ('R3', 'TSV temperature near a grid row but not on it', '@ T=293.15 col="Density (kg/m3)" field=[0] precision',
     '@ T=293.0 col="Density (kg/m3)" field=[0] precision'),
    # R6 - a class that must state something, stating nothing
    ('R6', 'unverified with no reason', '[UNVERIFIED] no source on disk; candidate NASA TR R-132 table 2', '[UNVERIFIED]'),
    # R7 - sampling policy claimed for a property
    ('R7', 'POLICY on a property table', '# @kind: range', '# @kind: property'),
    # R1 - a bracket that is not a C1.2 class
    ('R1', 'a bracket that never closes', '[POLICY: sampling-only]', '[POLICY: sampling-only'),
    # R4 - a relation demanded of a row holding more than one number, no field=
    ('R4', 'no field= on a multi-number row', 'field=[0] tol=0.5%', 'tol=0.5%'),
]


def selftest():
    bad = []
    _c, clean_f, clean_legacy = check_source('plant', _CLEAN)
    if clean_f:
        print('  the clean fixture itself fails - no plant can be judged:')
        for x in clean_f:
            print('    ' + x)
        return 1
    if len(clean_legacy) != 1:
        bad.append(f'R1 legacy: expected exactly 1 LEGACY tag in the fixture, got {clean_legacy}')
    for code, label, old, new in _PLANTS:
        if _CLEAN.count(old) != 1:
            bad.append(f'{code} {label}: plant anchor occurs {_CLEAN.count(old)} times')
            continue
        _c, f, _l = check_source('plant', _CLEAN.replace(old, new))
        fresh = [x for x in f if x not in clean_f and x.startswith(code)]
        status = 'ok' if fresh else 'FAIL'
        print(f'  [{status}] {code} {label}' + (f' -> {fresh[0][:110]}' if fresh else f' (got {f})'))
        if not fresh:
            bad.append(f'{code} {label}: not detected')
    # R2 - a file that exists but whose bytes are not the ones the manifest recorded:
    # a different surface form from a missing file, built on a temporary copy.
    import shutil
    import tempfile
    tmp = tempfile.mkdtemp(prefix='resolve_selftest_')
    try:
        os.makedirs(os.path.join(tmp, 'codata_2022'))
        dst = os.path.join(tmp, 'codata_2022', 'allascii.txt')
        shutil.copyfile(os.path.join(REFS, 'codata_2022', 'allascii.txt'), dst)
        with open(dst, 'a', encoding='utf-8') as fh:
            fh.write('\n')
        man = json.load(open(os.path.join(REFS, 'MANIFEST.json'), encoding='utf-8'))
        _c, f, _l = check_source('plant', _CLEAN.split('# @kind: property')[0],
                                 refs=tmp, manifest=man)
        fresh = [x for x in f if x.startswith('R2') and 'SHA-256' in x]
        print(f'  [{"ok" if fresh else "FAIL"}] R2 file present but altered after citation')
        if not fresh:
            bad.append(f'R2 altered file: not detected (got {f})')
    finally:
        shutil.rmtree(tmp, ignore_errors=True)

    # -- ATTACHMENT. Where a tag lands decides what it is checked against, and
    #    the first real run attached two wrongly while every failure plant above
    #    passed - there was no plant for attachment at all. Both were found by
    #    reading the LEGACY list the run printed. Plants from the rules:
    got = sorted(((f['table'], f['row'], f['bracket'].split(':')[0])
                  for f in extract_tags(_ATTACH)), key=lambda x: (x[0], x[1] or '', x[2]))
    want = sorted([
        ('SCALAR', None, 'VERIFY'),        # trails a scalar line; must not drift to a row
        ('ROWS', 'a', 'ON-DISK'),          # comment line above a row
        ('ROWS', 'b', 'UNVERIFIED'),       # trails its own row
        ('ROWS', None, 'UNVERIFIED'),      # above the closing brace: the table
        ('MIDLINE', None, 'ON-DISK'),      # mid-line, table has nothing better: kept
        ('MIDLINE', None, 'VERIFY'),
    ], key=lambda x: (x[0], x[1] or '', x[2]))   # ROWS' header prose mention: dropped
    print(f'  [{"ok" if got == want else "FAIL"}] attachment: six tags placed, one prose '
          f'mention dropped')
    if got != want:
        bad.append(f'attachment: got {got}, planted {want}')

    for x in bad:
        print('  - ' + x)
    print(f'selftest: {len(bad)} failure(s)')
    return 1 if bad else 0


_ATTACH = '''
SCALAR = 1.0   # [VERIFY: somewhere]
OTHER = 2.0

# the rows below carry line-start tags, so this [ON-DISK] mention is prose
ROWS = {
    # [ON-DISK] codata_2022/allascii.txt @ quantity="x"
    "a": 1.0,
    "b": 2.0,  # [UNVERIFIED] trailing its own row
    # [UNVERIFIED] sits above the closing brace
}

# see NAVFAC Ch. 3 [ON-DISK]  [VERIFY: Das]
MIDLINE = {"c": 3.0}
'''


if __name__ == '__main__':
    if '--selftest' in sys.argv:
        sys.exit(selftest())
    sys.exit(run(verbose='--verbose' in sys.argv))
