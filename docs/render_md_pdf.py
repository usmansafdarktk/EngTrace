"""Render a Markdown document to PDF through headless Microsoft Edge.

    python docs/render_md_pdf.py docs/PAPER_PLAN_OCT2026.md [docs/PAPER_PLAN_OCT2026.pdf]

The Markdown is converted with the `markdown` package (tables, fenced code, table of contents), wrapped in a
print stylesheet, written to the system temp directory, and printed by Edge's headless mode, which needs no
TeX or browser automation library. The script refuses to run without Edge and reports the page count (pypdf)
when the PDF is written. It is for working documents such as the paper plan, not for the paper itself.
"""
from __future__ import annotations

import os
import subprocess
import sys
import tempfile
from pathlib import Path

import markdown

EDGE_CANDIDATES = [
    r"C:\Program Files (x86)\Microsoft\Edge\Application\msedge.exe",
    r"C:\Program Files\Microsoft\Edge\Application\msedge.exe",
]

CSS = """
@page { size: A4; margin: 18mm 17mm 20mm 17mm; }
html { font-size: 10.5pt; }
body { font-family: Georgia, 'Times New Roman', serif; line-height: 1.38; color: #111; max-width: none; margin: 0; }
h1 { font-size: 1.75em; margin: 0 0 0.4em 0; line-height: 1.2; }
h2 { font-size: 1.3em; margin: 1.3em 0 0.4em 0; border-bottom: 1px solid #999; padding-bottom: 2px; page-break-after: avoid; }
h3 { font-size: 1.08em; margin: 1.1em 0 0.3em 0; page-break-after: avoid; }
p { margin: 0.45em 0; text-align: left; }
blockquote { margin: 0.6em 0 0.6em 1.2em; padding-left: 0.8em; border-left: 3px solid #bbb; color: #222; }
code { font-family: Consolas, 'Courier New', monospace; font-size: 0.86em; background: #f3f3f3; padding: 0 2px; }
pre { background: #f3f3f3; padding: 6px 8px; font-size: 0.82em; overflow-x: auto; white-space: pre-wrap; }
table { border-collapse: collapse; width: 100%; margin: 0.6em 0; font-size: 0.86em; page-break-inside: auto; }
th, td { border: 1px solid #999; padding: 3px 5px; vertical-align: top; text-align: left; }
th { background: #e8e8e8; }
tr { page-break-inside: auto; }
ul, ol { margin: 0.3em 0 0.5em 1.4em; padding-left: 0.6em; }
li { margin: 0.18em 0; }
hr { border: 0; border-top: 1px solid #999; margin: 1em 0; }
"""


def find_edge() -> str:
    for candidate in EDGE_CANDIDATES:
        if os.path.exists(candidate):
            return candidate
    sys.exit("Microsoft Edge was not found at its usual locations; install it or edit EDGE_CANDIDATES.")


def render(md_path: Path, pdf_path: Path) -> int:
    text = md_path.read_text(encoding="utf-8")
    body = markdown.markdown(
        text,
        extensions=["tables", "fenced_code", "sane_lists", "toc", "attr_list"],
        output_format="html5",
    )
    title = next((line.lstrip("# ").strip() for line in text.splitlines() if line.startswith("# ")), md_path.stem)
    html = (
        "<!doctype html><html><head><meta charset='utf-8'>"
        f"<title>{title}</title><style>{CSS}</style></head><body>{body}</body></html>"
    )
    html_path = Path(tempfile.gettempdir()) / (pdf_path.stem + ".render.html")
    html_path.write_text(html, encoding="utf-8")
    edge = find_edge()
    cmd = [
        edge,
        "--headless=new",
        "--disable-gpu",
        "--no-pdf-header-footer",
        f"--print-to-pdf={pdf_path.resolve()}",
        html_path.resolve().as_uri(),
    ]
    subprocess.run(cmd, check=False, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL, timeout=120)
    if not pdf_path.exists():
        sys.exit("Edge produced no PDF.")
    try:
        from pypdf import PdfReader

        pages = len(PdfReader(str(pdf_path)).pages)
    except Exception:  # pypdf absent or unreadable output
        pages = -1
    return pages


def main() -> int:
    if len(sys.argv) < 2:
        print(__doc__)
        return 2
    md_path = Path(sys.argv[1])
    pdf_path = Path(sys.argv[2]) if len(sys.argv) > 2 else md_path.with_suffix(".pdf")
    pages = render(md_path, pdf_path)
    print(f"{pdf_path} written" + (f", {pages} pages" if pages > 0 else ""))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
