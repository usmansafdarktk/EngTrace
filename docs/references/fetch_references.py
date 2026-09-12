"""Acquire every reference source Phases C1 and C3 are allowed to cite.

    python docs/references/fetch_references.py              # fetch anything missing
    python docs/references/fetch_references.py --only nist_fluid
    python docs/references/fetch_references.py --verify     # re-hash what is on disk
    python docs/references/fetch_references.py --list       # print the plan, fetch nothing

**Why this file exists.** C2 established the rule that makes a citation worth
anything: *every `[ON-DISK]` citation resolves to a file a reviewer can open
without network access and without a library* (D-030, D-034). C3 carries about
400 values across four branches, and the sources for most of them were never
acquired. This script is the acquisition step, made reproducible: it records
where each file came from, when, and its SHA-256, so a later reader can prove the
file they are citing is the file that was fetched.

**What it will not do.** It never writes a value into `constants.py`, and it
never extracts, converts or summarises a source. Transcription is C3's work and
must be reviewable against the raw artefact, so the raw artefact is what gets
stored. It also never copies a copyrighted book; those are listed, by path, as
local-only (see `LOCAL_ONLY_COPYRIGHTED`).

**How a download is judged.** An HTTP 200 is not evidence. NIST's saturation
endpoint silently ignores a temperature window, a mirror serves a "403
Forbidden" page with status 200, and a PDF link can return an HTML landing page.
So every file is checked by **content**: a PDF must start `%PDF`, a zip `PK`, a
NIST table must carry its column header and at least three data rows, a WebBook
page's `<title>` must name the species that CAS number was requested for, and a
JANAF table must declare the element in its reference state, and an HTML page's
title must name the page requested. Where the host publishes a hash, it is
checked too. A file that fails is recorded as a failure and is not written.

**Git, and why every copy is pinned.** The large PDFs, zips and `.xlsx` are NOT
committed (repo owner's decision; see `.gitignore`), so a fresh checkout gets
them from this script. Every citation into one was written against the file
MANIFEST.json recorded, so a re-acquired copy is held to that record -- its
SHA-256, or for a GitHub archive its commit id and a content fingerprint -- and a
mismatch is refused rather than recorded over (`_pin_ok`). Files already on disk
are held to it too.

Every source is public domain or openly licensed:

* NIST Chemistry WebBook, NIST-JANAF, CODATA — works of the US Government
* NASA TR R-132 (Svehla 1962) — NASA technical report, US Government work
* MIL-HDBK-5J — Distribution Statement A, approved for public release
* USDA FPL Wood Handbook GTR-282 — US Government work
* USGS, FHWA, NRCS, NAVFAC, MIL-STD-105E, NIST/SEMATECH — US Government works
* AISC Shapes Database v16.0 — distributed free of charge by AISC
* refractiveindex.info database — CC0-1.0
"""
from __future__ import annotations

import argparse
import hashlib
import io
import json
import os
import re
import shutil
import sys
import time
import urllib.error
import urllib.parse
import urllib.request
import zipfile
from datetime import datetime, timezone

HERE = os.path.dirname(os.path.abspath(__file__))
REPO = os.path.abspath(os.path.join(HERE, "..", ".."))
MANIFEST = os.path.join(HERE, "MANIFEST.json")
PILOT_PUBLIC = os.path.join(REPO, "pilot", "references", "public")

UA = "EngTrace-reference-acquisition/1.0 (research benchmark provenance; contact via repository)"
NIST_DELAY_S = 0.6


# ==========================================================================
# The plan
# ==========================================================================

#: NIST Chemistry WebBook fluid properties (SRD 69, thermophysical section).
#: Only fluids the database actually carries -- measured against its fluid list,
#: not assumed. Each gets its full saturation curve (which also yields Tc, Pc and
#: the critical density) and an isobar at 1 atm through 20 C.
NIST_FLUIDS = {
    "C7732185": "Water", "C7727379": "Nitrogen", "C1333740": "Hydrogen",
    "C7782447": "Oxygen", "C630080": "Carbon monoxide", "C124389": "Carbon dioxide",
    "C67561": "Methanol", "C74828": "Methane", "C74840": "Ethane",
    "C74851": "Ethene", "C74986": "Propane", "C115071": "Propene",
    "C106978": "Butane", "C75285": "Isobutane", "C109660": "Pentane",
    "C78784": "2-Methylbutane", "C110543": "Hexane", "C110827": "Cyclohexane",
    "C142825": "Heptane", "C111659": "Octane", "C7440597": "Helium",
    "C7440019": "Neon", "C7440371": "Argon", "C7439909": "Krypton",
    "C7440633": "Xenon", "C7664417": "Ammonia", "C75694": "R11",
    "C75718": "R12", "C75456": "R22", "C306832": "R123", "C811972": "R134a",
    "C71432": "Benzene", "C108883": "Toluene", "C7446095": "Sulfur dioxide",
    "C7783064": "Hydrogen sulfide", "C2551624": "Sulfur hexafluoride",
    # Added after the C3 review: the fuel and oil rows of FLUID_DENSITIES and
    # PIPE_FLUIDS named products with no single composition (kerosene, diesel,
    # jet fuel, SAE grades). These three ARE in the NIST fluid database and are
    # the hydrocarbons those rows can honestly be re-pointed at. CAS numbers read
    # off the WebBook name search: nonane 111-84-2, decane 124-18-5,
    # dodecane 112-40-3.
    "C111842": "Nonane", "C124185": "Decane", "C112403": "Dodecane",
}

#: Saturation states the constants tables STATE at a temperature. An unbounded
#: SatT request comes back on NIST's adaptive grid, which carries none of these
#: rows (phaseC1_summary.md section 10: 100 C water's nearest row is 375.68 K), and
#: reading one off it would be an interpolation. So each is requested on a 1-K
#: grid that STARTS at the stated temperature, and the file is refused unless that
#: row is present and two-phase. Each T is the table's own temperature + 273.15:
#:   REAL_FLUID_DATA temp_C (chemical) - every row whose fluid NIST carries;
#:     R-410A (a blend), ethanol and acetone are not in NIST_FLUIDS.
#:   FLUID_DENSITIES (mechanical) - "Liquid Nitrogen (at -196 C)", "Liquid Oxygen
#:     (at -183 C)", "Liquid Hydrogen (at -253 C)", "R-134a ... (Saturated Liquid
#:     at 25 C)", and "Liquid Propane" (its comment: "At 25 C, under pressure").
NIST_SATURATION_POINTS = {
    "C7732185": [373.15], "C7664417": [298.15], "C124389": [293.15],
    "C7446095": [298.15], "C74828": [112.15], "C74840": [184.15],
    "C74986": [298.15], "C106978": [298.15], "C75285": [298.15],
    "C109660": [298.15], "C78784": [298.15], "C75694": [298.15],
    "C75718": [298.15], "C75456": [298.15], "C811972": [298.15],
    "C306832": [298.15], "C108883": [384.15], "C71432": [353.15],
    "C67561": [338.15], "C110543": [342.15], "C111659": [399.15],
    "C110827": [354.15], "C7727379": [77.15], "C7782447": [90.15],
    "C1333740": [20.15],
}

#: WebBook species pages for substances the fluid database does NOT carry.
#: Mask=4 is phase-change data (Tc, Pc, Tboil, dHvap); Mask=2 is condensed-phase
#: thermochemistry (liquid/solid Cp). Each CAS is paired with the words its page
#: title must contain -- a wrong CAS is the error this check exists to catch.
WEBBOOK_SPECIES = {
    "C64175": ("Ethanol", ["ethanol"]),
    "C67641": ("Acetone", ["acetone"]),
    "C67630": ("Isopropanol", ["isopropyl", "propanol"]),
    "C56235": ("Carbon tetrachloride", ["tetrachloride"]),
    "C7439976": ("Mercury", ["mercury"]),
    "C7782505": ("Chlorine", ["chlorine"]),
    "C106423": ("p-Xylene", ["xylene"]),
    "C56815": ("Glycerol", ["glycer"]),
    "C107211": ("Ethylene glycol", ["ethanediol", "glycol"]),
    "C60297": ("Diethyl ether", ["ether"]),
    "C74862": ("Acetylene", ["acetylene"]),
    "C10043922": ("Radon", ["radon"]),
    # Added at C3 for MANOMETER_FLUIDS (the brief: establish each one's phase at
    # manometer conditions from a source). The CAS numbers were read off the
    # WebBook's own name search, not recalled: "tungsten hexafluoride" -> 7783-82-6.
    # "Tellurium Mercury" is Name Not Found; "mercury telluride" -> 12068-90-5 is
    # fetched as the nearest named species, which is NOT a finding that the row
    # means it.
    "C7783826": ("Tungsten hexafluoride", ["hexafluoride"]),
    "C12068905": ("Mercury telluride", ["telluride"]),
}

#: NIST-JANAF elements for solid heat capacities. File numbers are NOT
#: predictable from the symbol (Fe-001 is wustite, C-001 is niobium carbide), so
#: each is found by scanning IDs and accepting only a header that names the
#: element in its reference state.
JANAF_ELEMENTS = {
    "Fe": "Iron", "Cu": "Copper", "Al": "Aluminum", "Pb": "Lead",
    "W": "Tungsten", "Si": "Silicon", "C": "Carbon", "Hg": "Mercury",
    "Au": "Gold", "Ag": "Silver",
}

#: Single-file sources. `local` names a file in pilot/references/public/ to fall
#: back to when the canonical URL no longer serves the document (the NRCS link to
#: TR-55 is known broken); the fallback is recorded as such, never silently.
FILES = [
    dict(id="codata_2022", group="constants",
         url="https://physics.nist.gov/cuu/Constants/Table/allascii.txt",
         dest="codata_2022/allascii.txt", kind="text",
         must_contain="speed of light in vacuum",
         licence="US Government work (NIST)",
         grounds=["electrical_engineering: C0, EPSILON_0"]),
    # Added at C3: MIL-HDBK-5J states moduli in 10^3 ksi and densities in lb/in^3,
    # and the mechanical tables state GPa and kg/m^3. A conversion factor is a
    # value like any other: read from an artefact, not typed from memory.
    dict(id="nist_sp811_2008", group="constants",
         url="https://nvlpubs.nist.gov/nistpubs/Legacy/SP/nistspecialpublication811e2008.pdf",
         dest="nist_sp811/nistspecialpublication811e2008.pdf", kind="pdf",
         licence="US Government work (NIST)",
         grounds=["mechanical_engineering: unit conversions for MIL-HDBK-5J values "
                  "(ksi to GPa, lb/in^3 to kg/m^3) in the C3.2 consistency checks"]),
    dict(id="nasa_tr_r132", group="transport",
         url="https://ntrs.nasa.gov/api/citations/19630012982/downloads/19630012982.pdf",
         dest="nasa_tr_r132/svehla_1962_nasa_tr_r132.pdf", kind="pdf",
         licence="US Government work (NASA)",
         grounds=["chemical_engineering: GAS_MOLECULAR_PARAMS (Lennard-Jones sigma, eps/k)"]),
    dict(id="mil_hdbk_5j", group="materials",
         url="https://archive.org/download/milhdbk-5-j/MILHDBK5J.pdf",
         dest="mil_hdbk_5j/MIL-HDBK-5J_2003-01-31.pdf", kind="pdf",
         sha1="640dbd16ebf62af555c85d13a812adf8d02a6dd0",
         licence="Distribution Statement A: approved for public release",
         grounds=["mechanical_engineering: MATERIAL_PROPERTIES, SHEAR_MODULUS_VALUES, "
                  "MATERIAL_DENSITIES (metals)"]),
    dict(id="usda_wood_handbook", group="materials",
         url="https://research.fs.usda.gov/treesearch/download/62200.pdf",
         dest="usda_wood_handbook/FPL-GTR-282_2021.pdf", kind="pdf",
         licence="US Government work (USDA Forest Service)",
         grounds=["mechanical_engineering: MATERIAL_DENSITIES (woods)",
                  "chemical_engineering: SUBSTANCES_FOR_HEATING (wood)"]),
    # Pinned to a COMMIT, not to `main`, and verified by content rather than by
    # bytes. GitHub builds archives on request and does not guarantee byte-stable
    # output; measured 2026-09-11, the `main` archive first acquired and this
    # commit's archive differ in SHA-256 yet carry the same commit id in the zip
    # comment and the same 4198 files (path below the root folder, CRC-32, size).
    # The file keeps its original `-main` name so MANIFEST.json and citations
    # into it stay valid.
    dict(id="refractiveindex_info", group="optics",
         url="https://github.com/polyanskiy/refractiveindex.info-database/archive/"
             "c5c2f188e848453def5970e347399d653df2ffc2.zip",
         dest="refractiveindex_info/refractiveindex.info-database-main.zip", kind="zip",
         zip_commit="c5c2f188e848453def5970e347399d653df2ffc2",
         zip_fingerprint="cc8f4dccb6704634e7b9f39c2b0552a10fcfad1251fc99fe18c830568c465669",
         licence="CC0-1.0",
         grounds=["electrical_engineering: MEDIA_VELOCITIES (refractive indices)"]),
    # --- civil / industrial: already cited by those branches, but only ever
    #     present in the gitignored pilot/references/public/ on one machine.
    dict(id="usgs_wsp2339", group="civil",
         url="https://pubs.usgs.gov/wsp/2339/report.pdf",
         local="usgs_wsp2339_mannings_n.pdf",
         dest="civil/usgs_wsp2339_mannings_n.pdf", kind="pdf",
         licence="US Government work (USGS)",
         grounds=["civil_engineering: MANNINGS_N_CHANNELS"]),
    dict(id="fhwa_hds4", group="civil",
         url="https://www.fhwa.dot.gov/engineering/hydraulics/pubs/08090/HDS4_608.pdf",
         local="fhwa_hds4_highway_hydraulics.pdf",
         dest="civil/fhwa_hds4_highway_hydraulics.pdf", kind="pdf",
         licence="US Government work (FHWA)",
         grounds=["civil_engineering: MANNINGS_N_CHANNELS, MANNINGS_N_CONDUITS"]),
    dict(id="fhwa_hec22", group="civil",
         url="https://www.fhwa.dot.gov/engineering/hydraulics/pubs/10009/10009.pdf",
         local="fhwa_hec22_urban_drainage.pdf",
         dest="civil/fhwa_hec22_urban_drainage.pdf", kind="pdf",
         licence="US Government work (FHWA)",
         grounds=["civil_engineering: RATIONAL_C"]),
    # THIRD-PARTY MIRRORS. TR-55, both NAVFAC manuals and MIL-STD-105E were first
    # copied from pilot/references/public/ because no working URL was known. The
    # binaries are not committed (.gitignore), so a URL is what lets a fresh clone
    # have them at all. Each `url` below is a mirror of a public-domain document,
    # NOT the issuing agency, and is accepted for one reason only: its download
    # was measured BYTE-IDENTICAL (SHA-256) to the file MANIFEST.json already
    # recorded, on 2026-09-11. _pin_ok() holds every later download to that hash,
    # so a mirror that changes its copy is refused rather than trusted.
    #
    # The canonical NRCS link is BROKEN: pilot/references/public/MANIFEST.md
    # records it moving, curl returned HTTP 404, and urllib hung on it for over
    # ten minutes during the first full run, blocking every source behind it.
    # It is kept as a fact about the document, not retried. An NRC ADAMS link
    # (ML14219A437) tried as a mirror returned an HTML page, not the PDF.
    dict(id="nrcs_tr55", group="civil",
         url="https://www.oregon.gov/odot/hydraulics/Docs_Hydraulics_Manual/"
             "Hydraulics-07-G-Urban-Hydrology-Small-Watersheds.pdf",
         mirror_of="USDA NRCS TR-55 (1986); hosted by Oregon DOT in its Hydraulics Manual",
         url_verified="2026-09-11: download byte-identical (SHA-256) to the recorded file",
         canonical_url_broken="https://www.nrcs.usda.gov/sites/default/files/2022-10/TR-55%20%281986%29.pdf",
         local="nrcs_tr55_urban_hydrology.pdf",
         dest="civil/nrcs_tr55_urban_hydrology.pdf", kind="pdf",
         licence="US Government work (USDA NRCS)",
         grounds=["civil_engineering: SCS_CURVE_NUMBERS, SCS_IA_RATIO"]),
    dict(id="navfac_dm7_01", group="civil",
         url="https://vulcanhammer.net/wp-content/uploads/2017/01/dm7_01.pdf",
         mirror_of="NAVFAC DM-7.01 Soil Mechanics; hosted by vulcanhammer.net",
         url_verified="2026-09-11: download byte-identical (SHA-256) to the recorded file",
         local="navfac_dm7_01_soil_mechanics.pdf",
         dest="civil/navfac_dm7_01_soil_mechanics.pdf", kind="pdf",
         licence="US Government work (US Navy)",
         grounds=["civil_engineering: SPECIFIC_GRAVITY_RANGES, PERMEABILITY_RANGES_CM_S"]),
    dict(id="navfac_dm7_02", group="civil",
         url="https://vulcanhammer.net/wp-content/uploads/2017/01/dm7_02.pdf",
         mirror_of="NAVFAC DM-7.02 Foundations and Earth Structures; hosted by vulcanhammer.net",
         url_verified="2026-09-11: download byte-identical (SHA-256) to the recorded file",
         local="navfac_dm7_02_foundations.pdf",
         dest="civil/navfac_dm7_02_foundations.pdf", kind="pdf",
         licence="US Government work (US Navy)",
         grounds=["civil_engineering: TERZAGHI_BEARING_FACTORS (cross-check)"]),
    dict(id="aisc_shapes_v16", group="civil",
         url="https://cloud.aisc.org/biggie_bin/aisc-shapes-database-v160-2.xlsx",
         url_verified="2026-09-11: download byte-identical (SHA-256) to the recorded file",
         local="aisc_shapes_database_v16.xlsx",
         dest="civil/aisc_shapes_database_v16.xlsx", kind="zip",
         licence="Distributed free of charge by AISC",
         grounds=["civil_engineering: AISC_W_SHAPES"]),
    dict(id="mil_std_105e", group="industrial",
         url="https://www.expresscorp.com/wp-content/uploads/2023/02/MIL-STD-105E.pdf",
         mirror_of="MIL-STD-105E (1989); hosted by expresscorp.com",
         url_verified="2026-09-11: download byte-identical (SHA-256) to the recorded file",
         local="mil_std_105e_sampling.pdf",
         dest="industrial/mil_std_105e_sampling.pdf", kind="pdf",
         licence="Approved for public release; distribution unlimited",
         grounds=["industrial_engineering: MIL_STD_105E_*"]),
    # The whole-handbook zip is NOT recoverable from a clone. Measured 2026-09-11,
    # its URL answers 302 -> 301 -> NIST ITL's home page (HTML, 0.11 MB), and no
    # other copy is known, so the 174.8 MB file stays gitignored on the one
    # machine that has it. Its only use here was a secondary check on
    # CONTROL_CHART_FACTORS, so the two live pages that DEFINE those factors are
    # fetched instead: small, committed, and enough for a derivation citation.
    # They are the CURRENT pages, not the zip's copies -- the live pmc321.htm
    # (20 955 B) differs from the zip's (19 656 B) even with whitespace
    # normalised -- so cite the pages, not the zip.
    dict(id="nist_sematech_ehandbook", group="industrial",
         url="",
         canonical_url_broken="https://www.itl.nist.gov/div898/handbook/handbook.zip",
         local="nist_sematech_ehandbook.zip",
         dest="industrial/nist_sematech_ehandbook.zip", kind="zip",
         licence="US Government work (NIST)",
         grounds=["industrial_engineering: CONTROL_CHART_FACTORS (secondary check only; "
                  "factors are MathJax, not tables; local-only -- cite sematech_pmc32 and "
                  "sematech_pmc321 instead)"]),
    dict(id="sematech_pmc32", group="industrial",
         url="https://www.itl.nist.gov/div898/handbook/pmc/section3/pmc32.htm",
         dest="industrial/nist_sematech_pmc32.htm", kind="html",
         must_contain="What are Variables Control Charts?",
         licence="US Government work (NIST)",
         grounds=["industrial_engineering: CONTROL_CHART_FACTORS (c4, in closed form)"]),
    dict(id="sematech_pmc321", group="industrial",
         url="https://www.itl.nist.gov/div898/handbook/pmc/section3/pmc321.htm",
         dest="industrial/nist_sematech_pmc321.htm", kind="html",
         must_contain="Shewhart X-bar and R and S Control Charts",
         licence="US Government work (NIST)",
         grounds=["industrial_engineering: CONTROL_CHART_FACTORS (A2, D3, D4 defined from d2 "
                  "and d3; the page does NOT tabulate d2 or d3)"]),
]

#: Copyrighted books that exist on this machine only. Listed so C3 knows what can
#: be cited locally and what that citation cannot survive: a fresh clone. Never
#: copied, never committed, never redistributed (pilot/references/public/MANIFEST.md).
LOCAL_ONLY_COPYRIGHTED = [
    "full_books_civil_engineering",
    "full_books_industrial_engineering",
]


# ==========================================================================
# Machinery
# ==========================================================================

def _now():
    return datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")


def _hash(path, algo="sha256"):
    h = hashlib.new(algo)
    with open(path, "rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


class DeadlineExceeded(OSError):
    """A download that was still arriving when its wall-clock budget ran out."""


def _get(url, timeout=90, deadline_s=None):
    """Fetch ``url`` with a WALL-CLOCK deadline, not only a socket timeout.

    ``urllib``'s ``timeout`` bounds each socket operation, not the download: a
    server that keeps sending a byte every few seconds never trips it. The first
    full acquisition run hung for over ten minutes on the NRCS TR-55 link -- a
    URL already known to be broken -- and blocked every source queued behind it.
    Reading in chunks against a deadline turns a hang into a recorded failure,
    and the local-copy fallback then runs as designed.
    """
    deadline_s = deadline_s if deadline_s is not None else max(timeout * 4, 120)
    start = time.monotonic()
    req = urllib.request.Request(url, headers={"User-Agent": UA})
    with urllib.request.urlopen(req, timeout=timeout) as resp:
        chunks = []
        while True:
            if time.monotonic() - start > deadline_s:
                raise DeadlineExceeded(f"still downloading after {deadline_s:.0f} s")
            chunk = resp.read(1 << 20)
            if not chunk:
                break
            chunks.append(chunk)
        return resp.status, resp.headers.get("Content-Type", ""), b"".join(chunks)


def _content_ok(kind, data, must_contain=None):
    """Judge a download by what it is, never by its HTTP status."""
    head = data[:8]
    if kind == "pdf" and not head.startswith(b"%PDF"):
        return False, f"not a PDF (starts {head!r})"
    if kind == "zip" and not head.startswith(b"PK"):
        return False, f"not a zip/xlsx (starts {head!r})"
    if kind == "text":
        txt = data[:4096].decode("utf-8", "replace").lower()
        if "<html" in txt or "<!doctype" in txt:
            return False, "served an HTML page instead of text"
        if must_contain and must_contain.lower() not in data.decode("utf-8", "replace").lower():
            return False, f"does not contain {must_contain!r}"
    if kind == "json":
        # Without this a JSON source is judged only by its length, and PubChem
        # answers a bad request with a JSON *fault* object that is comfortably
        # over the size floor. So it must parse AND name what was asked for.
        try:
            json.loads(data.decode("utf-8"))
        except (UnicodeDecodeError, ValueError) as exc:
            return False, f"not JSON ({exc})"
        if must_contain and must_contain not in data.decode("utf-8", "replace"):
            return False, f"JSON does not contain {must_contain!r}"
    if kind == "html":
        # An HTML page is exactly what a redirect to a home page also is, so the
        # title must name the page that was asked for.
        title = re.findall(r"<title[^>]*>(.*?)</title>", data.decode("utf-8", "replace"), re.S | re.I)
        if not title or not must_contain or must_contain not in title[0]:
            return False, f"title {title[:1]!r} does not name {must_contain!r}"
    if len(data) < 256:
        return False, f"only {len(data)} bytes"
    return True, "ok"


def _zip_fingerprint(data):
    """Return (zip comment, SHA-256 over each member's path, CRC-32 and size).

    Paths are taken below the archive's root folder, because GitHub names that
    folder after the ref requested (`-main` versus `-<commit>`) while the files
    inside are the same.
    """
    with zipfile.ZipFile(io.BytesIO(data)) as z:
        names = z.namelist()
        root = names[0].split("/")[0] + "/" if names else ""
        rows = sorted((i.filename[len(root):], i.CRC, i.file_size)
                      for i in z.infolist() if not i.is_dir())
        comment = z.comment.decode("ascii", "replace")
    body = "\n".join(f"{n}\t{c:08x}\t{s}" for n, c, s in rows)
    return comment, hashlib.sha256(body.encode()).hexdigest()


def _pin_ok(spec, prior, data):
    """Hold a new copy of a source to the file MANIFEST.json already recorded.

    The large binaries are gitignored, so a fresh clone re-acquires them, and
    every citation into one was written against the RECORDED file. Without this
    check a host that has since changed its copy -- a revised PDF under the same
    URL, a mirror serving something else -- would be accepted and recorded over
    the old hash, and the citations would silently point at a document nobody
    checked them against. So a mismatch fails and nothing is written.
    """
    if spec.get("zip_fingerprint"):
        try:
            comment, fp = _zip_fingerprint(data)
        except zipfile.BadZipFile as exc:
            return False, f"not a readable zip: {exc}"
        if spec.get("zip_commit") and comment != spec["zip_commit"]:
            return False, f"archive is commit {comment!r}, pinned {spec['zip_commit']}"
        if fp != spec["zip_fingerprint"]:
            return False, f"archive content fingerprint {fp} != pinned {spec['zip_fingerprint']}"
        return True, "ok"
    want = prior.get("sha256") or prior.get("expected_sha256")
    got = hashlib.sha256(data).hexdigest()
    if want and got != want:
        return False, (f"sha256 {got} differs from the {want} MANIFEST.json recorded -- the "
                       f"source has changed, so citations into it need re-checking first")
    return True, "ok"


def _write(rel, data):
    path = os.path.join(HERE, rel)
    os.makedirs(os.path.dirname(path), exist_ok=True)
    with open(path, "wb") as fh:
        fh.write(data)
    return path


def _save(manifest):
    """Write the manifest atomically.

    Called after EVERY recorded source, not once at the end. The first version
    saved only on completion, so a run killed partway through -- a hung mirror, a
    closed laptop -- kept the downloaded files but lost their retrieval time,
    HTTP status and provenance route, and the next run could only re-describe
    them as "present". A file whose origin the manifest cannot state is exactly
    the kind of citation this project refuses to accept.
    """
    manifest["generated_utc"] = _now()
    manifest["generator"] = "docs/references/fetch_references.py"
    tmp = MANIFEST + ".tmp"
    # newline="\n": text mode on Windows would write CRLF, so the same manifest
    # would diff line-for-line between a Windows and a Linux run.
    with open(tmp, "w", encoding="utf-8", newline="\n") as fh:
        json.dump(manifest, fh, indent=1, sort_keys=True)
    # On Windows os.replace fails with PermissionError while ANY other process
    # holds MANIFEST.json open - an editor, an indexer, a virus scanner. Phase C1's
    # first run of this script died that way on the 78th of ~150 saves, leaving
    # the .tmp behind. A lock like that is transient, so retry before giving up.
    for attempt in range(20):
        try:
            os.replace(tmp, MANIFEST)
            return
        except PermissionError:
            if attempt == 19:
                raise
            time.sleep(0.25 * (attempt + 1))


def _record(manifest, entry):
    prior = manifest["sources"].get(entry["id"], {})
    if entry.get("status") == "present":
        # Keep what the run that acquired the file knew about it - INCLUDING the
        # URL. A present file was fetched from whichever URL succeeded, and for a
        # NIST isobar that is often the fallback grid, not the first-choice URL
        # the caller rebuilds. Before Phase C1 this list omitted "url", so every
        # re-run re-recorded the water and cyclohexane isobars under the URL NIST
        # had clamped - a manifest naming a request whose response is not the file
        # on disk.
        for key in ("url", "retrieved_utc", "via", "http", "content_type",
                    "published_sha1_verified", "retrieved_utc_note"):
            if prior.get(key) is not None and prior.get("status") != "failed":
                if key == "url" and entry.get("url_from_content"):
                    continue           # the caller proved which URL made the file
                if key == "url" or entry.get(key) is None:
                    entry[key] = prior[key]
        # A file acquired by a run that never wrote its manifest has no recorded
        # retrieval time. The first full run was stopped partway, after it had
        # written ~25 files, so this is not hypothetical. Say so, and give the
        # file's mtime as the best available time -- labelled as exactly that,
        # never presented as a retrieval timestamp. Applies to every source
        # type, not only the single-file ones.
        if not entry.get("retrieved_utc") and not entry.get("retrieved_utc_note"):
            p = os.path.join(HERE, entry["dest"])
            if os.path.exists(p):
                mt = datetime.fromtimestamp(os.path.getmtime(p), timezone.utc)
                entry["retrieved_utc_note"] = (
                    "acquired by a run stopped before its manifest was written; file mtime "
                    + mt.strftime("%Y-%m-%dT%H:%M:%SZ") + " is the best available retrieval time")
    manifest["sources"][entry["id"]] = entry
    _save(manifest)
    status = entry["status"]
    mark = {"acquired": "OK  ", "present": "have", "failed": "FAIL"}.get(status, status)
    extra = entry.get("reason") or entry.get("via", "")
    print(f"  [{mark}] {entry['id']:44s} {extra}", flush=True)


def fetch_file(spec, manifest):
    dest = os.path.join(HERE, spec["dest"])
    base = dict(id=spec["id"], group=spec["group"], dest=spec["dest"], url=spec["url"],
                licence=spec["licence"], grounds=spec["grounds"], kind=spec["kind"])
    for key in ("canonical_url_broken", "mirror_of", "url_verified", "zip_commit", "zip_fingerprint"):
        if spec.get(key):
            base[key] = spec[key]
    prior = manifest["sources"].get(spec["id"], {})
    want = prior.get("sha256") or prior.get("expected_sha256")

    def refuse(reason):
        # Keep the recorded hash in a failed entry, or the next run would have
        # nothing to hold a download to and would accept whatever arrives.
        base.update(status="failed", reason=reason)
        if want:
            base["expected_sha256"] = want
        _record(manifest, base)

    if os.path.exists(dest):
        base.update(status="present", bytes=os.path.getsize(dest), sha256=_hash(dest))
        # A present file is held to the manifest too: a file edited or replaced
        # on disk must not be re-recorded under its new hash as if it were the
        # source. A pinned archive is judged by content, since its bytes may
        # legitimately differ after a re-download.
        if spec.get("zip_fingerprint"):
            with open(dest, "rb") as fh:
                ok, why = _pin_ok(spec, prior, fh.read())
            if not ok:
                refuse(f"file on disk: {why}")
                return
        elif want and base["sha256"] != want:
            refuse(f"file on disk has sha256 {base['sha256']}, MANIFEST.json recorded {want}")
            return
        # A file already on disk is re-checked against the host's PUBLISHED hash
        # every run, not trusted because an earlier run once checked it. The
        # first run was killed before it wrote a manifest, so "an earlier run
        # verified this" was unrecordable anyway -- the check has to be
        # reproducible from the file alone.
        if spec.get("sha1"):
            got = _hash(dest, "sha1")
            if got != spec["sha1"]:
                refuse(f"file on disk has sha1 {got}, host publishes {spec['sha1']}")
                return
            base["published_sha1_verified"] = True
        if not prior.get("retrieved_utc"):
            base["retrieved_utc_note"] = (
                "acquired by a run that was stopped before its manifest was written; "
                "file mtime " + datetime.fromtimestamp(os.path.getmtime(dest), timezone.utc)
                .strftime("%Y-%m-%dT%H:%M:%SZ") + " is the best available retrieval time")
        _record(manifest, base)
        return
    reasons = []
    if spec["url"]:
        try:
            status, ctype, data = _get(spec["url"], timeout=600)
            ok, why = _content_ok(spec["kind"], data, spec.get("must_contain"))
            if ok and spec.get("sha1"):
                got = hashlib.sha1(data).hexdigest()
                if got != spec["sha1"]:
                    ok, why = False, f"sha1 {got} != published {spec['sha1']}"
            if ok:
                ok, why = _pin_ok(spec, prior, data)
            if ok:
                path = _write(spec["dest"], data)
                base.update(status="acquired", via="download", http=status, content_type=ctype,
                            bytes=len(data), sha256=_hash(path), retrieved_utc=_now(),
                            published_sha1_verified=bool(spec.get("sha1")))
                _record(manifest, base)
                return
            reasons.append(f"download rejected: {why}")
        except (urllib.error.URLError, TimeoutError, OSError) as exc:
            reasons.append(f"download failed: {type(exc).__name__}: {exc}")
    local = spec.get("local")
    if local:
        src = os.path.join(PILOT_PUBLIC, local)
        if os.path.exists(src):
            with open(src, "rb") as fh:
                data = fh.read()
            ok, why = _content_ok(spec["kind"], data)
            if ok:
                ok, why = _pin_ok(spec, prior, data)
            if ok:
                os.makedirs(os.path.dirname(dest), exist_ok=True)
                shutil.copyfile(src, dest)
                base.update(status="acquired",
                            via=f"local copy of pilot/references/public/{local} "
                                f"(canonical URL unavailable: {'; '.join(reasons) or 'none recorded'})",
                            bytes=len(data), sha256=_hash(dest), retrieved_utc=_now())
                _record(manifest, base)
                return
            reasons.append(f"local copy rejected: {why}")
        else:
            reasons.append(f"no local copy at pilot/references/public/{local}")
    refuse("; ".join(reasons) or "no URL and no local copy")


def _fluid_url(cas, **params):
    q = dict(Action="Data", Wide="on", ID=cas, Digits="5", RefState="DEF",
             TUnit="K", PUnit="MPa", DUnit="kg/m3", HUnit="kJ/mol",
             WUnit="m/s", VisUnit="Pa*s", STUnit="N/m")
    q.update(params)
    return "https://webbook.nist.gov/cgi/fluid.cgi?" + urllib.parse.urlencode(q)


def _tsv_ok(data, require_T=None):
    """Judge a NIST fluid table by its content -- including the grid it came back on.

    ``require_T``: a temperature (K) that MUST appear as a row. The isobar exists
    to supply 20 C values, and the first run showed why this has to be checked
    rather than assumed: **NIST does not reject a lower bound below the triple
    point, it silently clamps the grid to start there.** Water came back on
    273.16, 283.16, 293.16 ... and cyclohexane on 279.86, 289.86 ..., both with a
    valid header, a full set of rows, HTTP 200 -- and no 293.15 K row at all. The
    fallback grid never ran, because it was written to catch an error that the
    server never raises.
    """
    text = data.decode("utf-8", "replace")
    lines = [ln for ln in text.splitlines() if ln.strip()]
    if not lines or not lines[0].startswith("Temperature (K)"):
        return False, f"no NIST column header (starts {text[:60]!r})"
    if len(lines) < 4:
        return False, f"only {len(lines) - 1} data rows"
    if require_T is not None:
        temps = []
        for ln in lines[1:]:
            try:
                temps.append(float(ln.split("\t", 1)[0]))
            except ValueError:
                continue
        if not any(abs(t - require_T) < 0.005 for t in temps):
            span = f"{min(temps):.2f}..{max(temps):.2f}" if temps else "no numeric rows"
            return False, (f"grid has no {require_T} K row (came back on {span}) -- "
                           f"NIST clamped the requested bounds")
    return True, f"{len(lines) - 1} rows"


def _first_T(data):
    """The temperature of the first data row of a NIST fluid table, or None."""
    for ln in data.decode("utf-8", "replace").splitlines()[1:]:
        try:
            return float(ln.split("\t", 1)[0])
        except ValueError:
            continue
    return None


def _tlow(url):
    """The TLow a NIST isobar request asked for, or NaN if it names none."""
    q = urllib.parse.parse_qs(urllib.parse.urlparse(url).query)
    try:
        return float(q["TLow"][0])
    except (KeyError, IndexError, ValueError):
        return float("nan")


def _isobar_url_ok(entry, data):
    """Does the URL a manifest entry records actually produce the file it names?

    An isobar's first row is the TLow of the request that made it -- unless NIST
    clamped the request, which is exactly the case where the recorded URL and the
    file part company. So the two must agree. This is the check that would have
    caught the manifest recording the clamped water and cyclohexane requests.
    """
    first = _first_T(data)
    tl = _tlow(entry.get("url", ""))
    if first is None:
        return False, "file has no numeric data row"
    if tl != tl:                                   # NaN: the URL names no TLow
        return False, f"recorded URL names no TLow, so nothing ties it to a file starting at {first} K"
    if abs(tl - first) >= 0.005:
        return False, (f"recorded URL asks for TLow={tl} K but the file starts at {first} K -- "
                       f"that request did not produce this file")
    return True, "ok"


def fetch_nist_fluids(manifest):
    for cas, name in NIST_FLUIDS.items():
        slug = re.sub(r"[^a-z0-9]+", "_", name.lower()).strip("_")
        jobs = [
            ("saturation", _fluid_url(cas, Type="SatT"), None),
            # 20 C sits exactly on this grid. Fluids whose triple point is above
            # 253.15 K (water, benzene, cyclohexane) reject the lower bound, so
            # the fallback grid starts at 293.15 K.
            ("isobar_1atm", _fluid_url(cas, Type="IsoBar", P="0.101325",
                                       TLow="253.15", THigh="453.15", TInc="10"),
             _fluid_url(cas, Type="IsoBar", P="0.101325",
                        TLow="293.15", THigh="453.15", TInc="10")),
        ]
        if cas in NIST_ISOBAR_298:
            # 298.15 K is not on the 10 K grid above; this one lands on it exactly.
            # THigh spans several rows on purpose: _content_ok refuses a NIST table
            # with fewer than three data rows, and TLow/THigh 10 K apart gives two.
            jobs.append(("isobar_298K",
                         _fluid_url(cas, Type="IsoBar", P="0.101325",
                                    TLow="298.15", THigh="348.15", TInc="10"),
                         None))
        for what, url, fallback in jobs:
            rel = f"nist_fluid_properties/{slug}_{cas}_{what}.tsv"
            eid = f"nist_fluid:{slug}:{what}"
            # The isobar is only worth having if it carries the 20 C row the
            # tables need. The saturation curve has no fixed grid to check.
            require_T = 293.15 if what == "isobar_1atm" else None
            base = dict(id=eid, group="nist_fluid", dest=rel, url=url,
                        licence="US Government work (NIST SRD 69)", kind="tsv",
                        grounds=["chemical/mechanical fluid tables: density, viscosity, "
                                 "saturation volumes, Tc/Pc/rho_c"], species=name, cas=cas)
            p = os.path.join(HERE, rel)
            if os.path.exists(p):
                # A file on disk is re-judged, not trusted because it exists: the
                # first run wrote water and cyclohexane isobars on a clamped grid
                # with no 293.15 K row, and "present" would have kept them forever.
                with open(p, "rb") as fh:
                    ok_disk, why_disk = _tsv_ok(fh.read(), require_T)
                if ok_disk:
                    base.update(status="present", bytes=os.path.getsize(p), sha256=_hash(p))
                    if require_T is not None:
                        # Record the request that actually produced THIS file,
                        # proven from its content: an isobar's first row is the
                        # TLow of the request that made it. See _isobar_url_ok.
                        with open(p, "rb") as fh:
                            first = _first_T(fh.read())
                        for u in [url] + ([fallback] if fallback else []):
                            if first is not None and abs(_tlow(u) - first) < 0.005:
                                base["url"] = u
                                base["url_from_content"] = (
                                    f"the file's first row is {first} K, the TLow of this request")
                                break
                    _record(manifest, base)
                    continue
                print(f"  [redo] {eid:44s} file on disk rejected: {why_disk}", flush=True)
            tried = []
            for u in [url] + ([fallback] if fallback else []):
                try:
                    _, _, data = _get(u)
                    ok, why = _tsv_ok(data, require_T)
                except (urllib.error.URLError, TimeoutError, OSError) as exc:
                    ok, why = False, f"{type(exc).__name__}: {exc}"
                time.sleep(NIST_DELAY_S)
                if ok:
                    path = _write(rel, data)
                    base.update(status="acquired", via=why, url=u, bytes=len(data),
                                sha256=_hash(path), retrieved_utc=_now())
                    break
                tried.append(why)
            else:
                base.update(status="failed", reason="; ".join(tried))
            _record(manifest, base)


def _sat_row_ok(data, T):
    """A targeted saturation file must carry the row at T, and that row must be a
    two-phase state. The row check is `_tsv_ok`'s (NIST clamps bounds silently).
    The second is new: a request at or above the critical temperature returns the
    critical point, where the liquid and vapour columns are the same state."""
    ok, why = _tsv_ok(data, require_T=T)
    if not ok:
        return ok, why
    lines = [ln for ln in data.decode("utf-8", "replace").splitlines() if ln.strip()]
    head = lines[0].split("\t")
    try:
        jl, jv = head.index("Volume (l, m3/kg)"), head.index("Volume (v, m3/kg)")
    except ValueError:
        return False, "no liquid and vapour volume columns - not a saturation table"
    for ln in lines[1:]:
        cells = ln.split("\t")
        try:
            if abs(float(cells[0]) - T) < 0.005:
                if float(cells[jl]) == float(cells[jv]):
                    return False, (f"at {T} K the liquid and vapour volumes are equal "
                                   f"({cells[jl]}) - the critical point, not two phases")
                return True, why
        except (ValueError, IndexError):
            continue
    return False, f"no parseable row at {T} K"


def fetch_nist_saturation_points(manifest):
    for cas, temps in NIST_SATURATION_POINTS.items():
        name = NIST_FLUIDS[cas]
        slug = re.sub(r"[^a-z0-9]+", "_", name.lower()).strip("_")
        for T in temps:
            tag = f"{T:.2f}K"
            rel = f"nist_fluid_properties/{slug}_{cas}_saturation_{tag}.tsv"
            eid = f"nist_fluid:{slug}:saturation_{tag}"
            # NIST's `SatP` is the saturation table in TEMPERATURE increments and
            # `SatT` the one in pressure increments - read off the two form pages
            # (Action=Page: "Saturation Properties - Temperature Incremnts" carries
            # TLow/THigh/TInc under Type=SatP; the SatT page carries PLow/PHigh/PInc).
            # The first version sent TLow/THigh to SatT; NIST ignored them, returned
            # its adaptive curve for all 25 requests, and _sat_row_ok refused every one.
            url = _fluid_url(cas, Type="SatP", TLow=f"{T:.2f}", THigh=f"{T + 2:.2f}", TInc="1")
            base = dict(id=eid, group="nist_fluid", dest=rel, url=url,
                        licence="US Government work (NIST SRD 69)", kind="tsv",
                        grounds=["chemical_engineering: REAL_FLUID_DATA; mechanical_engineering: "
                                 "FLUID_DENSITIES saturated-liquid rows (a stated temperature)"],
                        species=name, cas=cas, saturation_T_K=T)
            p = os.path.join(HERE, rel)
            if os.path.exists(p):
                with open(p, "rb") as fh:
                    data = fh.read()
                ok_disk, why_disk = _sat_row_ok(data, T)
                ok_url, why_url = _isobar_url_ok(base, data)
                if ok_disk and ok_url:
                    base.update(status="present", bytes=os.path.getsize(p), sha256=_hash(p))
                    _record(manifest, base)
                    continue
                print(f"  [redo] {eid:44s} file on disk rejected: {why_disk if not ok_disk else why_url}",
                      flush=True)
            try:
                _, _, data = _get(url)
                ok, why = _sat_row_ok(data, T)
                if ok:
                    ok, why = _isobar_url_ok(base, data)
            except (urllib.error.URLError, TimeoutError, OSError) as exc:
                ok, why = False, f"{type(exc).__name__}: {exc}"
            time.sleep(NIST_DELAY_S)
            if ok:
                path = _write(rel, data)
                base.update(status="acquired", via=f"row at {T} K present and two-phase",
                            bytes=len(data), sha256=_hash(path), retrieved_utc=_now())
            else:
                base.update(status="failed", reason=why)
            _record(manifest, base)


def fetch_webbook_species(manifest):
    for cas, (name, title_words) in WEBBOOK_SPECIES.items():
        slug = re.sub(r"[^a-z0-9]+", "_", name.lower()).strip("_")
        for mask, what in (("4", "phase_change"), ("2", "condensed_phase")):
            rel = f"nist_webbook_species/{slug}_{cas}_{what}.html"
            url = f"https://webbook.nist.gov/cgi/cbook.cgi?ID={cas}&Units=SI&Mask={mask}"
            base = dict(id=f"webbook:{slug}:{what}", group="webbook", dest=rel, url=url,
                        licence="US Government work (NIST SRD 69)", kind="html",
                        grounds=["CRITICAL_PROPERTIES, SUBSTANCES_FOR_VAPORIZATION, "
                                 "SUBSTANCES_FOR_HEATING (species the fluid database lacks)"],
                        species=name, cas=cas)
            if os.path.exists(os.path.join(HERE, rel)):
                p = os.path.join(HERE, rel)
                base.update(status="present", bytes=os.path.getsize(p), sha256=_hash(p))
                _record(manifest, base)
                continue
            try:
                _, _, data = _get(url)
                time.sleep(NIST_DELAY_S)
                m = re.search(rb"<title>([^<]*)</title>", data, re.I)
                title = m.group(1).decode("utf-8", "replace") if m else ""
                if not any(w in title.lower() for w in title_words):
                    base.update(status="failed",
                                reason=f"page title {title!r} does not name {name} -- wrong CAS?")
                else:
                    path = _write(rel, data)
                    base.update(status="acquired", via=f"title {title!r}", bytes=len(data),
                                sha256=_hash(path), retrieved_utc=_now())
            except (urllib.error.URLError, TimeoutError, OSError) as exc:
                base.update(status="failed", reason=f"{type(exc).__name__}: {exc}")
            _record(manifest, base)


#: Substances whose density no on-disk artefact carries: NIST's fluid database
#: does not list them (checked against its own 74-fluid index) and the WebBook
#: species pages hold no liquid density. PubChem does, with a temperature and a
#: reference. The CID is resolved from the NAME by PubChem itself and recorded,
#: so a wrong CID is a fetch failure rather than a silent wrong substance.
PUBCHEM_DENSITY = {
    "ethanol": "FLUID_DENSITIES, PIPE_FLUIDS: Ethanol",
    "acetone": "FLUID_DENSITIES, PIPE_FLUIDS: Acetone",
    "2-propanol": "FLUID_DENSITIES, PIPE_FLUIDS: Isopropyl alcohol",
    "glycerol": "FLUID_DENSITIES, MANOMETER_FLUIDS: Glycerin",
    "ethylene glycol": "FLUID_DENSITIES, PIPE_FLUIDS: Ethylene glycol",
    "carbon tetrachloride": "FLUID_DENSITIES, MANOMETER_FLUIDS: Carbon tetrachloride",
    "chloroform": "FLUID_DENSITIES, MANOMETER_FLUIDS: Chloroform",
    "bromine": "FLUID_DENSITIES, MANOMETER_FLUIDS: Bromine",
    "mercury": "FLUID_DENSITIES, MANOMETER_FLUIDS, MATERIAL_DENSITIES: Mercury",
    "p-xylene": "FLUID_DENSITIES: Xylene",
    "diiodomethane": "MANOMETER_FLUIDS: Diiodomethane",
    "gallium": "MANOMETER_FLUIDS, MATERIAL_DENSITIES: Gallium",
    "tin": "MANOMETER_FLUIDS, MATERIAL_DENSITIES: Tin",
    "zinc": "MANOMETER_FLUIDS, MATERIAL_DENSITIES: Zinc",
    "copper": "MATERIAL_DENSITIES: Copper",
    "gold": "MATERIAL_DENSITIES: Gold",
    "silver": "MATERIAL_DENSITIES: Silver",
    "lead": "MATERIAL_DENSITIES: Lead",
    "nickel": "MATERIAL_DENSITIES: Nickel",
    "platinum": "MATERIAL_DENSITIES: Platinum",
    "tungsten": "MATERIAL_DENSITIES: Tungsten",
    "titanium": "MATERIAL_DENSITIES: Titanium",
    "aluminum": "MATERIAL_DENSITIES: Aluminum",
    "osmium": "MATERIAL_DENSITIES: Osmium",
    "uranium": "MATERIAL_DENSITIES: Uranium",
    "graphite": "MATERIAL_DENSITIES: Graphite",
}

#: Fluids that also need a grid landing on 298.15 K. The 1-atm isobar steps 10 K
#: from 253.15, so 25 degC is structurally absent from every file on disk - which is
#: why the COMMON_LIQUIDS viscosity diagnosis could not be tested either way (C3.5
#: register item 4). The four organics carry the 4.62-7.16% gap at the declared
#: 293.15 K; WATER is the control, because it is water agreeing at 20 degC that makes
#: the column non-uniform and so makes a 25 degC match meaningful.
NIST_ISOBAR_298 = {
    "C7732185": "Water",
    "C67561": "Methanol",
    "C71432": "Benzene",
    "C108883": "Toluene",
    "C110543": "Hexane",
}

PUBCHEM_NAME_URL = ("https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/name/"
                    "{name}/cids/JSON")
PUBCHEM_VIEW_URL = ("https://pubchem.ncbi.nlm.nih.gov/rest/pug_view/data/compound/"
                    "{cid}/JSON?heading=Density")


def fetch_pubchem_density(manifest):
    """One JSON per substance: PubChem's Density section, with its references.

    PubChem is an AGGREGATOR - its density entries cite CRC, Merck and USCG - so a
    citation to it is weaker than a primary source and stronger than an untagged
    value. The tags say so. The CID is resolved from the name by PubChem and
    recorded in the manifest entry, so the substance a file describes is never a
    guess (the same rule WEBBOOK_SPECIES follows for its CAS numbers).
    """
    for name, grounds in PUBCHEM_DENSITY.items():
        slug = re.sub(r"[^a-z0-9]+", "_", name.lower()).strip("_")
        rel = f"pubchem/{slug}_density.json"
        eid = f"pubchem:{slug}:density"
        base = dict(id=eid, group="pubchem", dest=rel, kind="json",
                    licence="US Government work (NIH/NLM PubChem); entries cite "
                            "third-party primaries (CRC, Merck) named in each record",
                    grounds=[grounds], species=name)
        p = os.path.join(HERE, rel)
        if os.path.exists(p):
            with open(p, "rb") as fh:
                ok_disk, why_disk = _content_ok("json", fh.read())
            if ok_disk:
                base.update(status="present", bytes=os.path.getsize(p), sha256=_hash(p),
                            url=base.get("url", ""))
                _record(manifest, base)
                continue
        try:
            _, _, cid_data = _get(PUBCHEM_NAME_URL.format(
                name=urllib.parse.quote(name)), timeout=60)
            cid = json.loads(cid_data.decode())["IdentifierList"]["CID"][0]
        except Exception as exc:                                   # noqa: BLE001
            base.update(status="failed", reason=f"name -> CID: {type(exc).__name__}: {exc}")
            _record(manifest, base)
            continue
        url = PUBCHEM_VIEW_URL.format(cid=cid)
        try:
            _, _, data = _get(url, timeout=90)
            ok, why = _content_ok("json", data, must_contain=str(cid))
        except (urllib.error.URLError, TimeoutError, OSError) as exc:
            ok, why = False, f"{type(exc).__name__}: {exc}"
        time.sleep(NIST_DELAY_S)
        if ok:
            path = _write(rel, data)
            base.update(status="acquired", via=why, url=url, cid=cid,
                        bytes=len(data), sha256=_hash(path), retrieved_utc=_now())
        else:
            base.update(status="failed", reason=why, url=url, cid=cid)
        _record(manifest, base)


def fetch_janaf(manifest, max_id=90):
    for sym, name in JANAF_ELEMENTS.items():
        rel = f"nist_janaf/{sym}_ref.txt"
        base = dict(id=f"janaf:{sym}", group="janaf", dest=rel, species=name,
                    licence="US Government work (NIST-JANAF, NSRDS-NBS 37)", kind="text",
                    grounds=["chemical_engineering: SUBSTANCES_FOR_HEATING (solid Cp)"])
        if os.path.exists(os.path.join(HERE, rel)):
            p = os.path.join(HERE, rel)
            base.update(status="present", bytes=os.path.getsize(p), sha256=_hash(p))
            _record(manifest, base)
            continue
        want = f"{name} ({sym})".lower()
        found, misses = None, 0
        for i in range(1, max_id + 1):
            url = f"https://janaf.nist.gov/tables/{sym}-{i:03d}.txt"
            try:
                _, _, data = _get(url, timeout=30)
            except (urllib.error.URLError, TimeoutError, OSError):
                data = b""
            time.sleep(NIST_DELAY_S / 2)
            first = data[:200].decode("utf-8", "replace").lower()
            if "<html" in first or not first.strip():
                misses += 1
                if misses >= 3 and i >= 3:
                    break              # the ID sequence for this element has ended
                continue
            misses = 0
            if want in first and "(ref)" in first:
                found = (url, data)
                break
        if found:
            path = _write(rel, found[1])
            base.update(status="acquired", url=found[0], via="header names element in ref state",
                        bytes=len(found[1]), sha256=_hash(path), retrieved_utc=_now())
        else:
            base.update(status="failed", url=f"https://janaf.nist.gov/tables/{sym}-NNN.txt",
                        reason=f"no table whose header names '{name} ({sym})' in its reference "
                               f"state -- NIST-JANAF does not carry it")
        _record(manifest, base)


def list_local_only(manifest):
    out = []
    for d in LOCAL_ONLY_COPYRIGHTED:
        root = os.path.join(PILOT_PUBLIC, d)
        if not os.path.isdir(root):
            continue
        for fn in sorted(os.listdir(root)):
            p = os.path.join(root, fn)
            if os.path.isfile(p) and fn.lower().endswith(".pdf"):
                out.append(dict(path=f"pilot/references/public/{d}/{fn}",
                                bytes=os.path.getsize(p)))
    manifest["local_only_copyrighted"] = dict(
        note=("Copyrighted textbooks present on ONE machine, under the gitignored "
              "pilot/references/. A citation to one of these cannot be resolved from a "
              "fresh clone. Never copy, commit or redistribute."),
        files=out)
    print(f"  [list] {len(out)} copyrighted books present locally, recorded by path only")


def verify(manifest):
    bad = 0
    for sid, e in sorted(manifest["sources"].items()):
        if e.get("status") not in ("acquired", "present"):
            continue
        p = os.path.join(HERE, e["dest"])
        if not os.path.exists(p):
            print(f"  [MISSING] {sid:44s} {e['dest']}")
            bad += 1
        elif _hash(p) != e.get("sha256"):
            print(f"  [CHANGED] {sid:44s} {e['dest']}")
            bad += 1
        elif e.get("group") == "nist_fluid" and sid.endswith(":isobar_1atm"):
            with open(p, "rb") as fh:
                ok, why = _isobar_url_ok(e, fh.read())
            if not ok:
                print(f"  [URL-NOT-FILE] {sid:40s} {why}")
                bad += 1
        elif e.get("group") == "nist_fluid" and "saturation_T_K" in e:
            with open(p, "rb") as fh:
                data = fh.read()
            for label, (ok, why) in (("URL-NOT-FILE", _isobar_url_ok(e, data)),
                                     ("NO-ROW", _sat_row_ok(data, e["saturation_T_K"]))):
                if not ok:
                    print(f"  [{label}] {sid:40s} {why}")
                    bad += 1
    print(f"verify: {bad} problem(s)")
    if bad:
        print("  The large PDFs, zips and .xlsx are gitignored. On a fresh checkout run this "
              "script without --verify first; it re-acquires them and refuses any copy that "
              "differs from what MANIFEST.json recorded.")
    return 1 if bad else 0


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--only", default="", help="group: constants, transport, materials, optics, "
                                              "civil, industrial, nist_fluid, webbook, "
                                              "pubchem, janaf")
    ap.add_argument("--verify", action="store_true")
    ap.add_argument("--list", action="store_true")
    ap.add_argument("--selftest", action="store_true",
                    help="planted defects for the manifest consistency checks")
    args = ap.parse_args(argv)
    if args.selftest:
        return selftest()

    manifest = {"sources": {}}
    if os.path.exists(MANIFEST):
        with open(MANIFEST, encoding="utf-8") as fh:
            manifest = json.load(fh)
        manifest.setdefault("sources", {})

    if args.verify:
        return verify(manifest)
    if args.list:
        for f in FILES:
            print(f"  {f['group']:10s} {f['id']:26s} {f['url'] or '(local fallback only)'}")
        print(f"  nist_fluid {len(NIST_FLUIDS)} fluids x 2 tables, plus "
              f"{sum(len(v) for v in NIST_SATURATION_POINTS.values())} targeted saturation rows")
        print(f"  webbook    {len(WEBBOOK_SPECIES)} species x 2 pages")
        print(f"  pubchem    {len(PUBCHEM_DENSITY)} substances x Density section")
        print(f"  janaf      {len(JANAF_ELEMENTS)} elements")
        return 0

    only = args.only
    for spec in FILES:
        if not only or spec["group"] == only:
            fetch_file(spec, manifest)
    if not only or only == "nist_fluid":
        fetch_nist_fluids(manifest)
        fetch_nist_saturation_points(manifest)
    if not only or only == "webbook":
        fetch_webbook_species(manifest)
    if not only or only == "pubchem":
        fetch_pubchem_density(manifest)
    if not only or only == "janaf":
        fetch_janaf(manifest)
    list_local_only(manifest)
    _save(manifest)                    # the same atomic, LF, lock-tolerant write

    st = [e["status"] for e in manifest["sources"].values()]
    print(f"\n{st.count('acquired')} acquired, {st.count('present')} already present, "
          f"{st.count('failed')} failed -- see MANIFEST.json")
    return 0


def selftest():
    """Planted defects for `_isobar_url_ok`, written from the defect it names.

    The defect: a manifest entry whose recorded URL is not the request that
    produced the file on disk. Two plants of different surface form - a URL that
    asks for the wrong grid, and a URL that names no grid at all - plus a file
    with no data, and a consistent pair that must pass. Per SPEC-CHANGE 17 a
    check that has only ever been shown the case it was built from proves
    nothing about the class.
    """
    header = b"Temperature (K)\tPressure (MPa)\tDensity (kg/m3)\n"
    file_293 = header + b"293.15\t0.10132\t998.21\n303.15\t0.10132\t995.65\n313.15\t0.10132\t992.22\n"
    iso = dict(Type="IsoBar", P="0.101325", THigh="453.15", TInc="10")
    cases = [
        ("recorded the clamped request: TLow 253.15, file starts 293.15",
         dict(url=_fluid_url("C7732185", TLow="253.15", **iso)), file_293, False),
        ("recorded a URL that names no TLow at all",
         dict(url=_fluid_url("C7732185", Type="SatT")), file_293, False),
        ("file carries no numeric data row",
         dict(url=_fluid_url("C7732185", TLow="293.15", **iso)), header + b"n/a\n", False),
        ("consistent request and file (negative control)",
         dict(url=_fluid_url("C7732185", TLow="293.15", **iso)), file_293, True),
    ]
    bad = 0
    for label, entry, data, want in cases:
        got, why = _isobar_url_ok(entry, data)
        ok = got == want
        bad += not ok
        print(f"  [{'ok' if ok else 'FAIL'}] {label}: {'passes' if got else 'flagged'} ({why})")
    # _sat_row_ok: a targeted saturation file (C3). Two defects of different form -
    # a grid without the stated row (a clamp), and a row that is the critical point
    # (a request at or above Tc) - plus a good file.
    sat_head = (b"Temperature (K)\tPressure (MPa)\tDensity (l, kg/m3)\tVolume (l, m3/kg)"
                b"\tDensity (v, kg/m3)\tVolume (v, m3/kg)\n")
    good = sat_head + (b"298.15\t1.0\t1206.7\t0.00082871\t32.4\t0.030857\n"
                       b"299.15\t1.0\t1203.0\t0.00083126\t33.4\t0.029940\n"
                       b"300.15\t1.1\t1199.2\t0.00083389\t34.4\t0.029062\n")
    clamped = sat_head + (b"299.15\t1.0\t1203.0\t0.00083126\t33.4\t0.029940\n"
                          b"300.15\t1.1\t1199.2\t0.00083389\t34.4\t0.029062\n"
                          b"301.15\t1.1\t1195.4\t0.00083655\t35.4\t0.028219\n")
    critical = sat_head + (b"374.21\t4.06\t511.9\t0.0019535\t511.9\t0.0019535\n"
                           b"375.21\t4.06\t511.9\t0.0019535\t511.9\t0.0019535\n"
                           b"376.21\t4.06\t511.9\t0.0019535\t511.9\t0.0019535\n")
    for label, data, T, want in (
            ("saturation grid without the stated row (clamped)", clamped, 298.15, False),
            ("saturation row at the critical point (equal volumes)", critical, 374.21, False),
            ("saturation file with a two-phase row at T (negative control)", good, 298.15, True)):
        got, why = _sat_row_ok(data, T)
        ok = got == want
        bad += not ok
        print(f"  [{'ok' if ok else 'FAIL'}] {label}: {'passes' if got else 'flagged'} ({why})")
    print(f"selftest: {bad} failure(s)")
    return 1 if bad else 0


if __name__ == "__main__":
    sys.exit(main())
