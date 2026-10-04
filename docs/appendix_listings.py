"""Build the appendix's template listings from the current template code, and check the .tex against them.

    python docs/appendix_listings.py           # print the three lstlisting blocks
    python docs/appendix_listings.py --check   # exit 1 unless overleaf_source_04102026/appendices/template_examples.tex
                                               # holds every block exactly (whitespace aside)

Each listing is the template function as it stands in data/templates/branches, shortened the way the May appendix
shortened its listings: the docstring is left out and nothing else. The one change is typographic: a non-ASCII em
dash inside a string becomes "--", because the pythonstyle listing has no mapping for multi-byte characters.
"""
from __future__ import annotations

import ast
import sys
from pathlib import Path

sys.stdout.reconfigure(encoding="utf-8", errors="replace")

ROOT = Path(__file__).resolve().parents[1]
TEX = ROOT / "overleaf_source_04102026/appendices/template_examples.tex"
LISTINGS = ["template_rational_method_peak_flow",            # Easy, civil
            "template_exponential_mttf_topology",            # Intermediate, industrial
            "template_equivalent_stiffness_frequency"]       # Advanced, mechanical


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
        bad = sorted({c for c in text if ord(c) > 127})
        if bad:
            raise SystemExit(f"{name}: non-ASCII characters left: {bad}")
        return "\\begin{lstlisting}[style=pythonstyle]\n" + text + "\n\\end{lstlisting}"
    raise SystemExit(f"{name} not found")


blocks = [listing(n) for n in LISTINGS]
if "--check" in sys.argv:
    tex = " ".join(TEX.read_text(encoding="utf-8").split())
    missing = [n for n, b in zip(LISTINGS, blocks) if " ".join(b.split()) not in tex]
    for n in missing:
        print("MISSING or changed:", n)
    print(f"{len(blocks) - len(missing)} of {len(blocks)} listings match the current template code")
    sys.exit(1 if missing else 0)
print("\n\n".join(blocks))
