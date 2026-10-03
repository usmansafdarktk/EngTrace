"""Extract the nine module icon badges from the 3-branch overview figure.

The icons in figures_oct_12/engtrace-overview.pdf are Type3 vector glyphs placed
by draw.io, so they cannot be recovered as Unicode emoji. This script crops each
badge out of the rendered page at 600 dpi and applies a circular alpha mask, so
the five-branch rebuild reuses the *same* artwork for the three original
branches instead of substituting look-alike icons.

Badge boxes were located by masking the module-pill fill (#FFF1CC) on the
rendered page and then finding the white badge disc inside each pill; the
resulting coordinates (PDF points, page 2850x900) are frozen below so the
extraction is reproducible without the detection pass.

Usage:  python figures_oct_12/extract_overview_icons.py
"""

import pathlib

import fitz
from PIL import Image, ImageDraw

HERE = pathlib.Path(__file__).resolve().parent
SRC = HERE / "engtrace-overview.pdf"
OUT = HERE / "assets" / "icons"
DPI = 600

# name -> badge box in PDF points (x0, y0, x1, y1) on the single page
BADGES = {
    "chem_transport": (131.0, 258.0, 174.0, 301.0),
    "chem_kinetics": (531.25, 256.0, 574.0, 298.75),
    "chem_thermo": (902.25, 259.75, 944.75, 302.5),
    "ee_signals": (1678.25, 272.0, 1721.0, 314.75),
    "ee_emwaves": (2031.0, 272.0, 2073.75, 314.75),
    "ee_digicomm": (2370.25, 274.5, 2413.0, 317.25),
    "me_fluid": (796.25, 658.75, 838.75, 701.5),
    "me_materials": (1211.25, 657.0, 1254.0, 699.75),
    "me_vibrations": (1682.5, 657.75, 1725.0, 700.5),
}


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    page = fitz.open(SRC)[0]
    for name, (x0, y0, x1, y1) in BADGES.items():
        pix = page.get_pixmap(dpi=DPI, clip=fitz.Rect(x0, y0, x1, y1))
        img = Image.frombytes("RGB", (pix.width, pix.height), pix.samples)
        size = min(img.size)
        img = img.crop((0, 0, size, size)).convert("RGBA")
        mask = Image.new("L", (size, size), 0)
        # inset by 1px so the pill fill never bleeds in at the rim
        ImageDraw.Draw(mask).ellipse((1, 1, size - 2, size - 2), fill=255)
        img.putalpha(mask)
        img.save(OUT / f"{name}.png")
        print(f"{name}.png  {size}x{size}")


if __name__ == "__main__":
    main()
