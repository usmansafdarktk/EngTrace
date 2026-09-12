"""Read a density out of a PubChem PUG-View record fetched with ?heading=Density.

    from tests.constants_integrity.pubchem import density_entries, density_value
    density_value(doc, ref=46) -> (8.94, {'temperature': None, 'source': 'Hazardous ...'})

PubChem (NIH/NLM) is an AGGREGATOR, and that is the whole reason a citation to it has to
name WHICH entry it means. A single record carries several densities that disagree:
copper is 8.94, 8.9 and 8.96 @25 C across five sources; mercury is 13.55 @68 F, 13.534
@25 C, 13.5 and 13.6. Taking the first, or averaging, would state a number no source
states. So a locator names the ReferenceNumber, and this reader returns that entry only.

The schema, mapped from the files rather than assumed (three attempts - the first two read
Record.Reference[], which carries the SourceName but never the citation):

    Record.RecordNumber                                     -> the CID
    Record.Section[].Section[].Section[].TOCHeading         -> "Density"
    Record.Section[]...Information[].ReferenceNumber        -> joins to Record.Reference[]
    Record.Section[]...Information[].Description            -> "PEER REVIEWED", sometimes
    Record.Section[]...Information[].Reference[]            -> the PRIMARY, as a string
    Record.Section[]...Information[].ExtendedReference[].Citation   -> the same, structured
    Record.Section[]...Information[].Value.StringWithMarkup[].String -> "0.7893 g/cu cm at 20 °C"

So the primary hangs off the INFORMATION entry, not off the Reference object. HSDB entries
routinely name CRC, Merck or Kirk-Othmer there; a tag can therefore say whose measurement
it is rather than implying PubChem made one.

Values are returned as the number the string LEADS with, together with the conditions the
string states, because those conditions are the point: "0.7893 g/cu cm at 20 °C" is usable
against a table declaring 20 C, and "6.0947 g/cu cm at 29.8 °C" is gallium ABOVE its
29.76 C melting point - a liquid density, not the solid one a materials table means. This
module reports the temperature; it does not decide whether the caller should want it.

Densities here are in g/cm^3. A table holding kg/m^3 cites scale=1000, so the comparison
happens in the artefact's own unit and precision= counts the digits PubChem prints.
"""
from __future__ import annotations

import re

#: the number a value string leads with
LEAD = re.compile(r'^\s*[<>~]?\s*(\d+(?:\.\d+)?)')
#: a RANGE, which is not a value: "18.7-19.3 @ 20 °C", "1.5-1.8", "1.8-3.51".
#: Tungsten's HSDB entry reads "18.7-19.3 @ 20 °C/4 °C; depends on extent of
#: working." Taking the leading number silently turns that into 18.7 while five
#: other sources in the same record say 19.3 - a wrong value with no signal that
#: a choice was made. A range is reported as a range and refused as a value.
RANGE = re.compile(r'^\s*\d+(?:\.\d+)?\s*[-‐-―]\s*\d+(?:\.\d+)?')
#: a stated temperature: "at 20 °C", "@25 °C", "at 68 °F"
TEMP = re.compile(r'(?:at|@)\s*(-?\d+(?:\.\d+)?)\s*°?\s*([CF])\b', re.I)


def _sections(node, path=()):
    if isinstance(node, dict):
        h = node.get('TOCHeading')
        p = path + (h,) if h else path
        if h:
            yield p, node
        for v in node.values():
            yield from _sections(v, p)
    elif isinstance(node, list):
        for v in node:
            yield from _sections(v, path)


def _primary(info):
    """The primary source this entry cites, or None."""
    for r in info.get('Reference') or []:
        if isinstance(r, str) and r.strip():
            return r.strip()
    for er in info.get('ExtendedReference') or []:
        if er.get('Citation'):
            return er['Citation']
    return None


def record_cid(doc):
    return (doc.get('Record') or {}).get('RecordNumber')


def density_entries(doc, heading='Density'):
    """Every entry under `heading`, in file order.

    Each is {'ref', 'source', 'peer', 'primary', 'string', 'value', 'temperature',
    'unit'} - 'value' is None when the string leads with no number (graphite's
    entries do), and 'temperature' is (number, 'C'|'F') or None.
    """
    rec = doc.get('Record') or {}
    names = {r.get('ReferenceNumber'): r.get('SourceName')
             for r in rec.get('Reference') or []}
    out = []
    for path, sec in _sections(rec):
        if path[-1] != heading:
            continue
        for info in sec.get('Information') or []:
            ref = info.get('ReferenceNumber')
            for swm in (info.get('Value') or {}).get('StringWithMarkup') or []:
                s = (swm.get('String') or '').strip()
                if not s:
                    continue
                m, t = LEAD.match(s), TEMP.search(s)
                is_range = bool(RANGE.match(s))
                out.append({
                    'ref': ref,
                    'source': names.get(ref),
                    'peer': (info.get('Description') or '') == 'PEER REVIEWED',
                    'primary': _primary(info),
                    'string': s,
                    'range': is_range,
                    'value': None if is_range else (float(m.group(1)) if m else None),
                    'temperature': (float(t.group(1)), t.group(2).upper()) if t else None,
                })
    return out


def density_value(doc, ref, heading='Density', cid=None):
    """(value, meta) for ONE cited entry. Raises ValueError rather than choosing.

    `cid`, when given, must equal the record's own RecordNumber: a locator that names
    a CID should fail loudly if the file on disk is a different substance.
    """
    if cid is not None and str(record_cid(doc)) != str(cid):
        raise ValueError(f'record is CID {record_cid(doc)}, not the cited CID {cid}')
    hits = [e for e in density_entries(doc, heading) if str(e['ref']) == str(ref)]
    if not hits:
        refs = sorted({str(e['ref']) for e in density_entries(doc, heading)})
        raise ValueError(f'no {heading} entry with ReferenceNumber {ref} (have: {refs})')
    if len(hits) > 1:
        raise ValueError(f'ReferenceNumber {ref} has {len(hits)} {heading} strings; '
                         f'a citation must name one value, not a set')
    e = hits[0]
    if e.get('range'):
        raise ValueError(f'entry {ref} states a RANGE, not a value: {e["string"]!r} - '
                         f'cite a source that prints one number, not the low end of a span')
    if e['value'] is None:
        raise ValueError(f'entry {ref} states no leading number: {e["string"]!r}')
    return e['value'], e
