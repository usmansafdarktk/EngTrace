"""Build the appendix's template listings from the current template code, and check the .tex against them.

    python docs/appendix_listings.py           # print the three lstlisting blocks
    python docs/appendix_listings.py --check   # exit 1 unless overleaf_source_04102026/appendices/template_examples.tex
                                               # holds every block exactly (whitespace aside)

Each listing is the template function as it stands in data/templates/branches, shortened the way the May appendix
shortened its listings: the docstring is left out and nothing else. The one change is typographic: a non-ASCII em
dash inside a string becomes "--", because the pythonstyle listing has no mapping for multi-byte characters; the
other non-ASCII characters a listing holds (a middle dot in a unit, superscripts and subscripts in a comment, an
approximately-equal sign) are printed through the listing's own `literate` option, so the code is shown unchanged.

The three templates, one per level from three branches: an Easy template whose question does not state the formula
it needs (the natural frequency of a torsional pendulum), the Intermediate template whose sampled topology changes
the derivation, and an Advanced
template that meets the paper's description of the level (coupled quantities, an implicit solve, and six or more
milestones in every instance's gold trace).
"""
from __future__ import annotations

import ast
import sys
from pathlib import Path

sys.stdout.reconfigure(encoding="utf-8", errors="replace")

ROOT = Path(__file__).resolve().parents[1]
TEX = ROOT / "overleaf_source_04102026/appendices/template_examples.tex"
LISTINGS = ["template_undamped_natural_frequency_torsional",  # Easy, mechanical
            "template_exponential_mttf_topology",            # Intermediate, industrial
            "template_vdw_solve_for_volume"]                 # Advanced, chemical
LITERATE = {"·": r"$\cdot$", "²": r"$^2$", "³": r"$^3$", "₀": r"$_0$", "₁": r"$_1$",
            "₂": r"$_2$", "≈": r"$\approx$"}  # printed by the listing, not changed in the code


def listing(name: str) -> str:
    for f in (ROOT / "data/templates/branches").rglob("*.py"):
        src = f.read_text(encoding="utf-8")
        if f"def {name}(" not in src:
            continue
        lines = src.splitlines()
        node = next(n for n in ast.parse(src).body if isinstance(n, ast.FunctionDef) and n.name == name)
        doc = node.body[0]
        skip = (set(range(doc.lineno, doc.end_lineno + 1))
                if isinstance(doc, ast.Expr) and isinstance(doc.value, ast.Constant) else set())
        body = [lines[i - 1].rstrip() for i in range(node.lineno, node.end_lineno + 1) if i not in skip]
        text = "\n".join(body).replace("\u2014", "--")
        wide = sorted({c for c in text if ord(c) > 127})
        bad = [c for c in wide if c not in LITERATE]
        if bad:
            raise SystemExit(f"{name}: non-ASCII characters left: {bad}")
        options = "style=pythonstyle" + (", literate=" + " ".join(f"{{{c}}}{{{{{LITERATE[c]}}}}}1" for c in wide)
                                         if wide else "")
        return f"\\begin{{lstlisting}}[{options}]\n" + text + "\n\\end{lstlisting}"
    raise SystemExit(f"{name} not found")


blocks = [listing(n) for n in LISTINGS]
# The text says that every evaluation-set instance of the Advanced listing has six or more milestones in its gold
# trace; the milestone file is local (full_run_28092026/scores/ is not committed), so the claim is checked where it
# exists
MILESTONES = ROOT / "full_run_28092026/scores/milestones.json"
if MILESTONES.exists():
    import json
    ms = json.loads(MILESTONES.read_text(encoding="utf-8"))
    adv = LISTINGS[2].replace("template_", "")
    counts = [len(v) for k, v in ms.items() if k.split("#")[0] == adv]
    assert len(counts) == 15 and min(counts) >= 6, (adv, counts)
if "--check" in sys.argv:
    tex = " ".join(TEX.read_text(encoding="utf-8").split())
    missing = [n for n, b in zip(LISTINGS, blocks) if " ".join(b.split()) not in tex]
    for n in missing:
        print("MISSING or changed:", n)
    print(f"{len(blocks) - len(missing)} of {len(blocks)} listings match the current template code")
    sys.exit(1 if missing else 0)
print("\n\n".join(blocks))
