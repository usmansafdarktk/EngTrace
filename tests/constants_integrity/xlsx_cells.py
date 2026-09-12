"""Read named cells out of an .xlsx with the standard library (C3.7, AISC Shapes Database).

    from tests.constants_integrity.xlsx_cells import Workbook
    wb = Workbook('docs/references/civil/aisc_shapes_database_v16.xlsx')
    wb.cell('Database v16.0', label_col='AISC_Manual_Label', label='W8X24', col='Ix', block=1)
      -> '82.7'

An .xlsx is a zip of XML parts. This reads the workbook's sheet map, the shared-string
table and one worksheet, and nothing else - no formulas are evaluated, no styles read:
a cell's value is the string its <v> element stores. openpyxl is installed on the
authoring machine but is not in requirements.txt, and a constants check should not add
a dependency.

The AISC database repeats every header name twice: a US-customary block (columns
B..CF) and an SI block (CG..FJ), measured on the on-disk file. `block=1` is the first
occurrence of a header name reading left to right, `block=2` the second. A header
name or row label that is absent - or, for a label, present in more than one row -
raises LookupError: a cell is never guessed.
"""
from __future__ import annotations

import re
import zipfile


def _col_key(ref):
    letters = re.match(r'[A-Z]+', ref).group(0)
    n = 0
    for ch in letters:
        n = n * 26 + (ord(ch) - 64)
    return n


class Workbook:
    def __init__(self, path):
        self.path = path
        with zipfile.ZipFile(path) as z:
            self._wb = z.read('xl/workbook.xml').decode('utf-8')
            rels = z.read('xl/_rels/workbook.xml.rels').decode('utf-8')
            try:
                ss = z.read('xl/sharedStrings.xml').decode('utf-8')
            except KeyError:
                ss = ''
            self._targets = dict(re.findall(r'<Relationship [^>]*Id="([^"]+)"[^>]*Target="([^"]+)"', rels))
            self._sheets = dict(re.findall(r'<sheet [^>]*name="([^"]+)"[^>]*r:id="([^"]+)"', self._wb))
            self._strings = [''.join(re.findall(r'<t[^>]*>(.*?)</t>', si, re.S))
                             for si in re.findall(r'<si>(.*?)</si>', ss, re.S)]
            self._xml = {}
            for name, rid in self._sheets.items():
                target = self._targets[rid].lstrip('/')
                part = target if target.startswith('xl/') else 'xl/' + target
                self._xml[name] = z.read(part).decode('utf-8')
        self._rows = {}

    def sheet_names(self):
        return list(self._sheets)

    def rows(self, sheet):
        """[(row number, {column letters: value string})] in file order."""
        if sheet not in self._xml:
            raise LookupError(f'no sheet {sheet!r} (sheets: {self.sheet_names()})')
        if sheet not in self._rows:
            out = []
            for rnum, rxml in re.findall(r'<row r="(\d+)"[^>]*>(.*?)</row>', self._xml[sheet], re.S):
                cells = {}
                for ref, attrs, inner in re.findall(r'<c r="([A-Z]+)\d+"([^>]*)>(.*?)</c>', rxml, re.S):
                    v = re.search(r'<v>(.*?)</v>', inner, re.S)
                    if v is None:
                        continue
                    val = v.group(1)
                    if 't="s"' in attrs:
                        val = self._strings[int(val)]
                    cells[ref] = val
                out.append((int(rnum), cells))
            self._rows[sheet] = out
        return self._rows[sheet]

    def column(self, sheet, header, block=1):
        """The column letters of the `block`-th occurrence of `header` in row 1."""
        head = self.rows(sheet)[0][1]
        hits = sorted((c for c, v in head.items() if v == header), key=_col_key)
        if len(hits) < block:
            raise LookupError(f'header {header!r} occurs {len(hits)} time(s) in {sheet!r}, not block {block}')
        return hits[block - 1]

    def cell(self, sheet, label_col, label, col, block=1):
        """The value in column `col` (block-th occurrence) of the ONE row whose
        `label_col` (same block) equals `label`."""
        lc = self.column(sheet, label_col, block)
        vc = self.column(sheet, col, block)
        found = [cells for _r, cells in self.rows(sheet)[1:] if cells.get(lc) == label]
        if len(found) != 1:
            raise LookupError(f'{label!r} in {label_col!r} (block {block}) matches {len(found)} rows of {sheet!r}')
        if vc not in found[0]:
            raise LookupError(f'{label!r} has no value in {col!r} (block {block})')
        return found[0][vc]
