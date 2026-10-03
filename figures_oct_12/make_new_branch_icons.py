"""Draw the six module icons for the two branches added after the 3-branch figure.

The nine original icons are recovered verbatim by extract_overview_icons.py. The
Civil and Industrial modules have no artwork in the old figure, so these are
drawn here in the same visual language: a white disc with a #FFF1CC inner ring,
a flat two- or three-colour mark with rounded chunky strokes, filling roughly
two thirds of the disc.

Usage:  python figures_oct_12/make_new_branch_icons.py
"""

import pathlib

OUT = pathlib.Path(__file__).resolve().parent / "assets" / "icons"

# white disc + the pale-yellow inner ring the extracted badges have
BADGE = (
    '<circle cx="50" cy="50" r="50" fill="#FFFFFF"/>'
    '<circle cx="50" cy="50" r="44.5" fill="none" stroke="#FFF1CC" stroke-width="2.5"/>'
)

MARKS = {
    # Civil — Structural Analysis: a Warren truss with pinned joints
    "civil_structural": """
      <g stroke="#5B4EC9" stroke-width="4.6" stroke-linecap="round"
         stroke-linejoin="round" fill="none">
        <path d="M17 64 H83"/>
        <path d="M28 33 H72"/>
        <path d="M17 64 L28 33 L39 64 L50 33 L61 64 L72 33 L83 64"/>
      </g>
      <g fill="#9B8CF2">
        <circle cx="17" cy="64" r="5"/><circle cx="39" cy="64" r="5"/>
        <circle cx="61" cy="64" r="5"/><circle cx="83" cy="64" r="5"/>
        <circle cx="28" cy="33" r="5"/><circle cx="50" cy="33" r="5"/>
        <circle cx="72" cy="33" r="5"/>
      </g>
      <g fill="#5B4EC9">
        <path d="M11 74 l6 -8 l6 8 z"/><path d="M77 74 l6 -8 l6 8 z"/>
      </g>
    """,
    # Civil — Geotechnical: layered strata with a borehole sample
    # soil strata with a cone penetrometer driven into them
    "civil_geotech": """
      <path d="M14 38 h72 v13 h-72 z" fill="#E0B94A"/>
      <path d="M14 51 h72 v13 h-72 z" fill="#B5804F"/>
      <path d="M14 64 h72 v10 a4 4 0 0 1 -4 4 h-64 a4 4 0 0 1 -4 -4 z" fill="#7A6250"/>
      <path d="M14 38 h72 v36 a4 4 0 0 1 -4 4 h-64 a4 4 0 0 1 -4 -4 z"
            fill="none" stroke="#4A3B30" stroke-width="3.6" stroke-linejoin="round"/>
      <g stroke="#4A3B30" stroke-width="3" stroke-linecap="round">
        <path d="M14 51 h72"/><path d="M14 64 h72"/>
      </g>
      <path d="M50 16 v40" stroke="#3C3C46" stroke-width="5" stroke-linecap="round"/>
      <path d="M43 56 h14 l-7 12 z" fill="#D7D9E0" stroke="#3C3C46"
            stroke-width="3.2" stroke-linejoin="round"/>
      <path d="M40 16 h20" stroke="#3C3C46" stroke-width="5" stroke-linecap="round"/>
    """,
    # Civil — Water Resources: open-channel flow in a trapezoidal section
    "civil_water": """
      <path d="M22 44 L30 70 H70 L78 44 Z" fill="#8FD3F4"/>
      <path d="M14 22 L28 72 H72 L86 22" fill="none" stroke="#2F6FB5"
            stroke-width="5" stroke-linecap="round" stroke-linejoin="round"/>
      <path d="M21 44 q7 -7 14 0 t14 0 t14 0 t14 0" fill="none" stroke="#1BA0C4"
            stroke-width="4.6" stroke-linecap="round"/>
      <g stroke="#2F6FB5" stroke-width="4" stroke-linecap="round" stroke-linejoin="round"
         fill="none">
        <path d="M37 58 h16"/><path d="M53 58 l-5 -4.5"/><path d="M53 58 l-5 4.5"/>
      </g>
    """,
    # Industrial — Production & Inventory: stacked cartons on a pallet
    "ie_production": """
      <rect x="20" y="44" width="26" height="24" rx="2.5" fill="#F2A03D"/>
      <rect x="54" y="44" width="26" height="24" rx="2.5" fill="#F2A03D"/>
      <rect x="37" y="20" width="26" height="24" rx="2.5" fill="#FFC773"/>
      <g stroke="#8A5314" stroke-width="3.2" stroke-linejoin="round" fill="none">
        <rect x="20" y="44" width="26" height="24" rx="2.5"/>
        <rect x="54" y="44" width="26" height="24" rx="2.5"/>
        <rect x="37" y="20" width="26" height="24" rx="2.5"/>
      </g>
      <g stroke="#8A5314" stroke-width="3" stroke-linecap="round">
        <path d="M33 44 v24"/><path d="M67 44 v24"/><path d="M50 20 v24"/>
      </g>
      <g stroke="#6E5849" stroke-width="4" stroke-linecap="round">
        <path d="M16 74 h68"/><path d="M24 74 v6"/><path d="M50 74 v6"/><path d="M76 74 v6"/>
      </g>
    """,
    # Industrial — Quality & Reliability: a Shewhart chart with 3-sigma limits
    "ie_quality": """
      <g stroke="#E8505B" stroke-width="3.4" stroke-linecap="round" stroke-dasharray="7 6">
        <path d="M18 30 H82"/><path d="M18 70 H82"/>
      </g>
      <path d="M18 50 H82" stroke="#9AA0A6" stroke-width="3" stroke-linecap="round"/>
      <path d="M22 58 L34 42 L46 56 L58 38 L70 52 L80 34" fill="none"
            stroke="#2E9E6B" stroke-width="4.4" stroke-linecap="round" stroke-linejoin="round"/>
      <g fill="#2E9E6B">
        <circle cx="22" cy="58" r="4"/><circle cx="34" cy="42" r="4"/>
        <circle cx="46" cy="56" r="4"/><circle cx="58" cy="38" r="4"/>
        <circle cx="70" cy="52" r="4"/>
      </g>
      <circle cx="80" cy="34" r="4.6" fill="#E8505B"/>
    """,
    # Industrial — Stochastic Operations: arrivals queueing into a server
    "ie_stochastic": """
      <g fill="#7B5BE0">
        <circle cx="17" cy="50" r="6"/><circle cx="33" cy="50" r="6"/><circle cx="49" cy="50" r="6"/>
      </g>
      <path d="M58 50 H66" stroke="#4A3B9E" stroke-width="3.8" stroke-linecap="round"/>
      <path d="M66 50 l-5 -4.5" stroke="#4A3B9E" stroke-width="3.8" stroke-linecap="round"/>
      <path d="M66 50 l-5 4.5" stroke="#4A3B9E" stroke-width="3.8" stroke-linecap="round"/>
      <circle cx="80" cy="50" r="13" fill="#25C4C9" stroke="#116E78" stroke-width="3.4"/>
      <path d="M80 43 v7 l5 4" fill="none" stroke="#FFFFFF" stroke-width="3.2"
            stroke-linecap="round" stroke-linejoin="round"/>
      <path d="M10 34 q20 -12 40 0" fill="none" stroke="#B9A9F5"
            stroke-width="3.2" stroke-linecap="round" stroke-dasharray="6 5"/>
      <path d="M10 66 q20 12 40 0" fill="none" stroke="#B9A9F5"
            stroke-width="3.2" stroke-linecap="round" stroke-dasharray="6 5"/>
    """,
}


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    for name, mark in MARKS.items():
        svg = (
            '<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 100 100" width="100" height="100">'
            f"{BADGE}{mark.strip()}</svg>"
        )
        (OUT / f"{name}.svg").write_text(svg, encoding="utf-8")
        print(f"{name}.svg")


if __name__ == "__main__":
    main()
