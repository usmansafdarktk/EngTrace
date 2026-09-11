"""C3.2 - every consumer of a table with a declared validity domain stays inside it (D-032).

    python -m tests.constants_integrity.test_consumer_domains
    python -m tests.constants_integrity.test_consumer_domains --selftest

D-032: "A correct value used outside its range is a defect. Declare the domain;
assert the consumers stay inside it." C2's plausibility suite was green because it
never looked where the data was used (CP_PARAMS integrated 90% past its fit).

For every numeric table whose `# @domain:` is not `none (...)`:

  D1  every consumer the census names is covered: by a PROBE - a function of the
      consumer's frame locals at return that yields the domain value the instance
      used, in the domain's unit - or by an IMPLICIT claim that the consumer states
      no such condition, checked against the consumer's own source. A consumer
      with neither, a probe that observes nothing in SEEDS seeds, or an implicit
      claim its source contradicts, is a failure: an unasserted consumer is the
      gap D-032 is about.
  D2  every observed value lies inside the declared domain: `key=lo..hi unit`
      inclusive, `key=v unit` exactly, a word value as a label. A domain key no
      probe observes is reported as UNOBSERVED - counted, never passed as inside.
  D3  an excursion (table, consumer, key) fails unless domain_findings.txt lists it
      with a referral - and a listed excursion that no longer occurs fails too.
      Listing is not acceptance; the register ratchets like the @domain worklist.
"""
from __future__ import annotations

import importlib
import inspect
import os
import re
import sys

REPO = os.path.abspath(os.path.join(os.path.dirname(__file__), '..', '..'))
if REPO not in sys.path:
    sys.path.insert(0, REPO)

from tests.constants_integrity.census import census  # noqa: E402
from tests.constants_integrity.given_evidence import capture  # noqa: E402
from tests.constants_integrity.test_table_metadata import _DOMAIN_PAIR, _split_top  # noqa: E402

SEEDS = 60
FINDINGS = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'domain_findings.txt')


def parse_domain(decl):
    """{key: (kind, lo, hi, unit)}, or None for `none (reason)`."""
    decl = decl.strip()
    if re.match(r'^none \(.+\)$', decl):
        return None
    out = {}
    for p in _split_top(decl):
        m = _DOMAIN_PAIR.match(p)
        if not m:
            raise ValueError(f'@domain entry {p!r} is not key=value [unit]')
        val, unit = m.group('val'), (m.group('unit') or '').strip()
        if re.match(r'^-?\d', val):
            lo, _, hi = val.partition('..')
            out[m.group('key')] = ('range', float(lo), float(hi or lo), unit)
        else:
            out[m.group('key')] = ('label', val, val, unit)
    return out


def inside(spec, value):
    kind, lo, hi, _unit = spec
    return value == lo if kind == 'label' else lo <= value <= hi


# ---------------------------------------------------------------- probes and claims
def _wave_parameters_basic(loc, cmod):
    """The VACUUM wavelength of the instance's wave, in um. refractiveindex.info
    tabulates n against vacuum wavelength; both branches of the template compute f."""
    f = loc.get('f')
    return {'wavelength': cmod.C0 / f * 1e6} if f else {}


PROBES = {
    ('electrical_engineering', 'MEDIA_VELOCITIES', 'template_wave_parameters_basic'): _wave_parameters_basic,
}

_TEMPERATURE = r'temperature|°C|\bdeg(?:rees)? ?C\b|kelvin|\d ?K\b'
IMPLICIT = {
    ('mechanical_engineering', 'FLUID_DENSITIES', tid): (
        _TEMPERATURE, 'states no fluid temperature, so the density is the table\'s 20 °C value; '
                      'P is not asserted (the template computes pressures)')
    for tid in ('template_hydrostatic_pressure_at_depth', 'template_basic_buoyant_force',
                'template_floating_object_submersion_depth', 'template_hydrostatic_force_on_plane')
}


def evaluate(branch, table, domain, fns, cmod, seeds=SEEDS, probes=PROBES, implicit=IMPLICIT):
    """{template_id: result} for one table's consumers. fns: {template_id: function}."""
    results = {}
    for tid, fn in fns.items():
        key3 = (branch, table, tid)
        if key3 in implicit:
            pattern, reason = implicit[key3]
            hit = re.search(pattern, inspect.getsource(fn), re.I)
            results[tid] = (dict(status='IMPLICIT', detail=reason) if not hit else
                            dict(status='IMPLICIT-CLAIM-FALSE',
                                 detail=f'the implicit claim is contradicted: its source says {hit.group(0)!r}'))
            continue
        probe = probes.get(key3)
        if probe is None:
            results[tid] = dict(status='NO-PROBE', detail='no probe and no implicit claim')
            continue
        observed, outside_, used = {}, {}, 0
        for s in range(seeds):
            _q, loc = capture(fn, s)
            vals = probe(loc, cmod)
            if not vals:
                continue
            used += 1
            for k, v in vals.items():
                observed.setdefault(k, []).append(v)
                if k in domain and not inside(domain[k], v):
                    outside_.setdefault(k, []).append((s, v))
        results[tid] = dict(status='PROBED' if used else 'OBSERVED-NOTHING', used=used,
                            observed=observed, outside=outside_,
                            unobserved=sorted(set(domain) - set(observed)))
    return results


def failures_for(branch, table, results):
    fails, excursions = [], set()
    for tid, r in results.items():
        if r['status'] in ('NO-PROBE', 'IMPLICIT-CLAIM-FALSE', 'OBSERVED-NOTHING'):
            fails.append(f'D1 {branch}.{table} {tid}: {r["status"]} - {r.get("detail", "")}')
        for k in r.get('outside', {}):
            excursions.add((f'{branch}.{table}', tid, k))
    return fails, excursions


def read_register(path=FINDINGS):
    out = {}
    if not os.path.exists(path):
        return out
    for ln in open(path, encoding='utf-8'):
        ln = ln.strip()
        if not ln or ln.startswith('#'):
            continue
        head, _, why = ln.partition('|')
        parts = head.split()
        if len(parts) != 3 or not why.strip():
            raise ValueError(f'register line is not "branch.TABLE template key | referral": {ln!r}')
        out[tuple(parts)] = why.strip()
    return out


def reconcile(excursions, registered):
    fails = [f'D3 {t} {tid} {k}: outside its declared domain and not in domain_findings.txt'
             for t, tid, k in sorted(excursions - set(registered))]
    fails += [f'D3 {t} {tid} {k}: listed in domain_findings.txt but no longer occurs - remove it'
              for t, tid, k in sorted(set(registered) - excursions)]
    return fails


def run(seeds=SEEDS):
    report = census(seeds=0, probe_on=False)
    fails, excursions, n_tables = [], set(), 0
    for t in report:
        decl = t['declared'].get('domain')
        dom = parse_domain(decl) if decl else None
        if not dom:
            continue
        n_tables += 1
        cmod = importlib.import_module(f'data.templates.branches.{t["branch"]}.constants')
        fns = {tid: getattr(importlib.import_module(u['module']), tid)
               for tid, u in t['consumers'].items() if 'module' in u}
        res = evaluate(t['branch'], t['name'], dom, fns, cmod, seeds)
        print(f'  {t["branch"][:5]}.{t["name"]}  @domain: {decl}')
        for tid, r in res.items():
            if r['status'] == 'PROBED':
                span = {k: (min(v), max(v)) for k, v in r['observed'].items()}
                print(f'    {tid:44s} PROBED {r["used"]}/{seeds}  observed {span}'
                      + (f'  OUTSIDE {sorted(r["outside"])}' if r['outside'] else '')
                      + (f'  unobserved {r["unobserved"]}' if r['unobserved'] else ''))
            else:
                print(f'    {tid:44s} {r["status"]}: {r.get("detail", "")}')
        f, e = failures_for(t['branch'], t['name'], res)
        fails += f
        excursions |= e
    registered = read_register()
    fails += reconcile(excursions, registered)
    print(f'{n_tables} tables with a declared domain; {len(excursions)} excursion(s), '
          f'{len(registered)} registered')
    for x in fails:
        print('  - ' + x)
    print('all pass' if not fails else f'{len(fails)} FAILURES')
    return 1 if fails else 0


# ---------------------------------------------------------------- self-test
class _FakeConstants:
    C0 = 299792458


def _fx_optical():
    f = 5.093e14                          # a vacuum wavelength of 0.5886 um
    return 'q', f'{f}'


def _fx_radio():
    f = 1.0e8
    return 'q', f'{f}'


def _fx_hot():
    T = 350.0
    return 'q', f'{T}'


def _fx_quiet():
    return 'The block floats in water.', 'a'


def _fx_warm():
    return 'Water at 40 °C.', 'a'


def selftest():
    bad = []
    point = parse_domain('wavelength=0.589 um')
    rng = parse_domain('T=280..300 K')

    def _t_probe(loc, _c):
        return {'T': loc['T']} if 'T' in loc else {}

    def _w_probe(loc, c):
        return {'wavelength': round(c.C0 / loc['f'] * 1e6, 3)} if 'f' in loc else {}

    probes = {('b', 'W', '_fx_optical'): _w_probe, ('b', 'W', '_fx_radio'): _w_probe,
              ('b', 'R', '_fx_hot'): _t_probe}
    implicit = {('b', 'R', '_fx_quiet'): (_TEMPERATURE, 'planted'),
                ('b', 'R', '_fx_warm'): (_TEMPERATURE, 'planted')}
    rw = evaluate('b', 'W', point, {'_fx_optical': _fx_optical, '_fx_radio': _fx_radio}, _FakeConstants,
                  seeds=2, probes=probes, implicit=implicit)
    # '_fx_unprobed' is a real function covered by neither a probe nor a claim
    rr = evaluate('b', 'R', rng, {'_fx_hot': _fx_hot, '_fx_quiet': _fx_quiet, '_fx_warm': _fx_warm,
                                  '_fx_unprobed': _fx_quiet}, _FakeConstants,
                  seeds=2, probes=probes, implicit=implicit)
    fw, ew = failures_for('b', 'W', rw)
    fr, er = failures_for('b', 'R', rr)
    cases = [
        ('control: a consumer inside a POINT domain', not rw['_fx_optical']['outside']),
        ('a POINT domain left by a consumer (radio for an optical index)', ('b.W', '_fx_radio', 'wavelength') in ew),
        ('a RANGE domain left above its top (350 K for 280..300 K)', ('b.R', '_fx_hot', 'T') in er),
        ('control: an implicit claim the source supports', rr['_fx_quiet']['status'] == 'IMPLICIT'),
        ('an implicit claim the source contradicts ("40 °C")', rr['_fx_warm']['status'] == 'IMPLICIT-CLAIM-FALSE'),
        ('a consumer with no probe and no claim', rr['_fx_unprobed']['status'] == 'NO-PROBE'
         and any('_fx_unprobed' in x for x in fr)),
    ]
    reg = {('b.W', '_fx_radio', 'wavelength'): 'planted referral', ('b.R', '_fx_gone', 'T'): 'planted'}
    rec = reconcile(ew | er, reg)
    cases += [
        ('an excursion not in the register', any('_fx_hot' in x and 'not in domain_findings' in x for x in rec)),
        ('a registered excursion that no longer occurs', any('_fx_gone' in x and 'no longer occurs' in x for x in rec)),
        ('control: a registered excursion that occurs', not any('_fx_radio' in x for x in rec)),
    ]
    for label, ok in cases:
        print(f"  [{'ok' if ok else 'FAIL'}] {label}")
        if not ok:
            bad.append(label)
    print(f'selftest: {len(bad)} failure(s)')
    return 1 if bad else 0


if __name__ == '__main__':
    sys.exit(selftest() if '--selftest' in sys.argv else run())
