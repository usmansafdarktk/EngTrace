"""Read a refractive index out of the on-disk refractiveindex.info archive.

    from tests.constants_integrity.refractiveindex import index_at
    n, how, conditions = index_at("database/data/main/H2O/nk/Daimon-20.0C.yml", 0.589)

The archive is `docs/references/refractiveindex_info/refractiveindex.info-database-main.zip`
(CC0-1.0, pinned by commit and content fingerprint in MANIFEST.json). It stores most
materials as DISPERSION FORMULAS, not as tabulated indices, so an index at a
wavelength is a computation. The formulas are transcribed here from the database's
OWN definition, `database/doc/Dispersion formulas.pdf` (RefractiveIndex.INFO,
2014-06-29), inside the same archive - not from memory. lambda is in micrometres.

    1  Sellmeier     n^2-1 = C1 + sum C(2i) l^2 / (l^2 - C(2i+1)^2)
    2  Sellmeier-2   n^2-1 = C1 + sum C(2i) l^2 / (l^2 - C(2i+1))
    3  Polynomial    n^2   = C1 + sum C(2i) l^C(2i+1)
    4  RI.INFO       n^2   = C1 + C2 l^C3/(l^2 - C4^C5) + C6 l^C7/(l^2 - C8^C9) + sum_{i>=5} C(2i) l^C(2i+1)
    5  Cauchy        n     = C1 + sum C(2i) l^C(2i+1)
    6  Gases         n-1   = C1 + sum C(2i) / (C(2i+1) - l^-2)
    7  Herzberger    n     = C1 + C2/(l^2-0.028) + C3 (1/(l^2-0.028))^2 + C4 l^2 + C5 l^4 + C6 l^6
    8  Retro         (n^2-1)/(n^2+2) = C1 + C2 l^2/(l^2 - C3) + C4 l^2
    9  Exotic        n^2   = C1 + C2/(l^2 - C3) + C4 (l - C5)/((l - C5)^2 + C6)

A wavelength outside a formula's stated range is REFUSED, never extrapolated. A
tabulated dataset returns the linearly interpolated value AND its two bracketing
rows, so a caller can require that the rounding it cites does not depend on the
interpolation.
"""
from __future__ import annotations

import math
import os
import zipfile

REPO = os.path.abspath(os.path.join(os.path.dirname(__file__), '..', '..'))
ARCHIVE = os.path.join(REPO, 'docs', 'references', 'refractiveindex_info',
                       'refractiveindex.info-database-main.zip')


def _pairs(C, start):
    i = start
    while i + 1 <= len(C):
        yield C[i - 1], C[i]
        i += 2


def formula(kind, C, lam):
    C = list(C) + [0.0] * (17 - len(C))
    lam2 = lam * lam
    if kind == 1:
        return math.sqrt(1 + C[0] + sum(a * lam2 / (lam2 - b * b) for a, b in _pairs(C, 2) if a))
    if kind == 2:
        return math.sqrt(1 + C[0] + sum(a * lam2 / (lam2 - b) for a, b in _pairs(C, 2) if a))
    if kind == 3:
        return math.sqrt(C[0] + sum(a * lam ** b for a, b in _pairs(C, 2) if a))
    if kind == 4:
        s = C[0]
        if C[1]:
            s += C[1] * lam ** C[2] / (lam2 - C[3] ** C[4])
        if C[5]:
            s += C[5] * lam ** C[6] / (lam2 - C[7] ** C[8])
        s += sum(a * lam ** b for a, b in _pairs(C, 10) if a)
        return math.sqrt(s)
    if kind == 5:
        return C[0] + sum(a * lam ** b for a, b in _pairs(C, 2) if a)
    if kind == 6:
        return 1 + C[0] + sum(a / (b - lam ** -2) for a, b in _pairs(C, 2) if a)
    if kind == 7:
        L = 1 / (lam2 - 0.028)
        return C[0] + C[1] * L + C[2] * L * L + C[3] * lam2 + C[4] * lam ** 4 + C[5] * lam ** 6
    if kind == 8:
        r = C[0] + C[1] * lam2 / (lam2 - C[2]) + C[3] * lam2
        return math.sqrt((1 + 2 * r) / (1 - r))
    if kind == 9:
        return math.sqrt(C[0] + C[1] / (lam2 - C[2]) + C[3] * (lam - C[4]) / ((lam - C[4]) ** 2 + C[5]))
    raise ValueError(f'unknown dispersion formula {kind}')


def _load(member, archive=ARCHIVE):
    import yaml
    with zipfile.ZipFile(archive) as z:
        names = z.namelist()
        root = names[0].split('/')[0] + '/'
        path = member if member in names else root + member
        if path not in names:
            raise KeyError(f'{member!r} is not in {os.path.basename(archive)}')
        return yaml.safe_load(z.read(path).decode('utf-8', 'replace'))


def index_at(member, lam_um, archive=ARCHIVE):
    """Return (n, how, conditions, brackets) or raise ValueError/KeyError.

    `brackets` is None for a formula, or ((l0, n0), (l1, n1)) for tabulated data.
    """
    doc = _load(member, archive)
    cond = doc.get('CONDITIONS') or {}
    for d in doc.get('DATA', []):
        kind = str(d.get('type', ''))
        if kind.startswith('formula'):
            lo, hi = (float(x) for x in str(d['wavelength_range']).split())
            if not lo <= lam_um <= hi:
                raise ValueError(f'{kind} is valid {lo}-{hi} um; {lam_um} um is outside it')
            C = [float(x) for x in str(d['coefficients']).split()]
            return formula(int(kind.split()[1]), C, lam_um), f'{kind} ({lo}-{hi} um)', cond, None
        # Exactly `tabulated n` or `tabulated nk`. The first version tested
        # startswith('tabulated n'), which also admits `tabulated n2` - the
        # NONLINEAR index, ~1e-20 m^2/W - and "evaluated" air at n = 0.0000000.
        if kind in ('tabulated n', 'tabulated nk'):
            rows = sorted([float(x) for x in ln.split()]
                          for ln in str(d['data']).strip().splitlines())
            for (w0, n0, *_), (w1, n1, *_) in zip(rows, rows[1:]):
                if w0 <= lam_um <= w1:
                    n = n0 if w1 == w0 else n0 + (n1 - n0) * (lam_um - w0) / (w1 - w0)
                    return n, f'{kind}, interpolated {w0}-{w1} um', cond, ((w0, n0), (w1, n1))
            raise ValueError(f'{kind} has no rows bracketing {lam_um} um')
    raise ValueError('no DATA block with a formula or n data')
