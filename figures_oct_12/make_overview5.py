"""Rebuild the EngTrace overview figure with all five branches.

The original figure (figures_oct_12/engtrace-overview.pdf, drawn in draw.io) shows
three branches. Civil and Industrial Engineering were added later, under
data/templates/branches/, so this script regenerates the same figure from the
taxonomy as data:

  branch  ->  data/templates/branches/<branch>/
  module  ->  a package directory inside it
  topic   ->  a template module (.py) inside that directory

Palette, geometry and typography are taken from the original PDF: page fill
#FFFFFF, panel #FAFAF4 on a #EBDFD0 rule, topic pill #F2ECE3, module pill
#FFF1CC, branch pill #F6DBE2, title pill #EBDFD0, hub ring #8C6E62 / #EBDFD0,
Times New Roman (bold for headings, regular for topics).

The nine original module icons are reused verbatim (extract_overview_icons.py);
the six new ones are drawn to match (make_new_branch_icons.py).

Layout: two rows of three columns, wide rather than tall, so the figure spans
the full text width at the same type size and takes less of the page:

    Chemical      [title / hub]   Electrical
    Civil          Mechanical     Industrial

The arrows are axis-aligned. Chemical and Electrical get straight runs to the
middle of their inner edge; Civil and Industrial fork off those runs and drop
down the gutters; Mechanical gets a straight drop onto its branch pill. Each
head is drawn as its own triangle so the tip lands exactly at the edge it
points to, and every tail starts under the hub ring.

Renders HTML, then prints it to PDF with headless Chrome.

Usage:  python figures_oct_12/make_overview5.py
"""

import math
import pathlib
import subprocess
import sys

from PIL import ImageFont

HERE = pathlib.Path(__file__).resolve().parent
ICONS = HERE / "assets" / "icons"
HTML_OUT = HERE / "engtrace-overview-5branch.html"
PDF_OUT = HERE / "engtrace-overview-5branch.pdf"

CHROME = r"C:\Program Files\Google\Chrome\Application\chrome.exe"
FONT_REG = r"C:\Windows\Fonts\times.ttf"
FONT_BOLD = r"C:\Windows\Fonts\timesbd.ttf"

# ---------------------------------------------------------------- palette ---
C_PAGE = "#FFFFFF"
C_FRAME = "#EBDFD0"
C_PANEL = "#FAFAF4"
C_TOPIC = "#F2ECE3"
C_TOPIC_EDGE = "#E4DACB"
C_MODULE = "#FFF1CC"
C_MODULE_SHADOW = "#DCD4C6"
C_BRANCH = "#F6DBE2"
C_TITLE = "#EBDFD0"
C_HUB_DARK = "#8C6E62"
C_HUB_LIGHT = "#EBDFD0"
C_ARROW = "#808080"
C_INK = "#000000"

# --------------------------------------------------------------- geometry ---
PANEL_W = 1180
MARGIN_X = 80       # page edge to the outer panels
GUTTER = 150        # between columns; the forked arrows drop through it
COL_X = [MARGIN_X + i * (PANEL_W + GUTTER) for i in range(3)]
PAGE_W = 2 * MARGIN_X + 3 * PANEL_W + 2 * GUTTER   # 4000
CENTER_X = PAGE_W // 2

FRAME_X, FRAME_Y, FRAME_R, FRAME_STROKE = 30, 62, 24, 5
FRAME_W = PAGE_W - 2 * FRAME_X
FRAME_PAD_B = 50    # bottom panels to the frame, matching the 45px inner margin at the sides
PAGE_PAD_B = 30     # frame to the page edge, as on the left and right

ROW_A_TOP = 176
ROW_GAP = 86

PANEL_BORDER = 9
PILL_TOP = -34      # branch pill offset from the panel's padding box

PANEL_PAD_T = 82
PANEL_PAD_B = 40
PANEL_PAD_X = 28
COL_GAP = 22
MOD_GAP = 18
PILL_GAP = 12

FS_TITLE = 44       # "EngTrace Taxonomy"
FS_HUB = 60         # "EngTrace"
FS_BRANCH = 38      # branch name
FS_MODULE = 26      # module name
FS_TOPIC = 24       # topic name
FS_BADGE = 26       # the circled branch number

MOD_LINE_H = 34
MOD_PAD_Y = 7
MOD_ICON = 44
MOD_ICON_GAP = 12
TOPIC_LINE_H = 30
TOPIC_PAD_Y = 5
TOPIC_PAD_X = 13

HUB_R = 190
HUB_RING_R = 152

ARROW_W = 8
HEAD_L = 38         # arrow head length along the shaft
HEAD_W = 36         # arrow head width at its base
ARROW_GAP = 4       # tip to the edge it points at
CORNER_R = 28       # radius of each bend
FORK_DX = 100       # a fork's vertical run, measured in from the outer panel's edge

# -------------------------------------------------------------- taxonomy ----
# (icon file, module label, [topic labels]) -- the display names for the
# directories and .py files under data/templates/branches/.
BRANCHES = [
    dict(
        slot="A-left",
        name="Chemical Engineering",
        modules=[
            ("chem_transport.png", "Transport Phenomena",
             ["Shell Momentum Balances", "Viscosity &amp; Momentum Transport"]),
            ("chem_kinetics.png", "Reaction Kinetics",
             ["Mole Balances", "Stoichiometry", "Conversion &amp; Reactor Sizing"]),
            ("chem_thermo.png", "Thermodynamics",
             ["Heat Effects", "Volumetric Properties of Pure Fluids"]),
        ],
    ),
    dict(
        slot="A-right",
        name="Electrical Engineering",
        modules=[
            ("ee_signals.png", "Signals &amp; Systems",
             ["Discrete Time Signals", "Continuous Time Signals"]),
            ("ee_emwaves.png", "Electromagnetics &amp; Waves",
             ["Electrostatics", "Magnetostatics", "Waves &amp; Phasors"]),
            ("ee_digicomm.png", "Digital Communications",
             ["Deterministic &amp; Random Signal Analysis", "Digital Modulation Schemes"]),
        ],
    ),
    dict(
        slot="B-left",
        name="Civil Engineering",
        modules=[
            ("civil_structural.svg", "Structural Analysis",
             ["Determinate Structures", "Influence Lines", "Deflections",
              "Indeterminate Analysis"]),
            ("civil_geotech.svg", "Geotechnical Engineering",
             ["Phase &amp; Index Properties",
              "Permeability, Seepage &amp; Effective Stress",
              "Stress Distribution &amp; Consolidation", "Strength &amp; Stability"]),
            ("civil_water.svg", "Water Resources",
             ["Uniform Flow", "Energy &amp; Rapidly Varied Flow", "Hydrology"]),
        ],
    ),
    dict(
        slot="B-right",
        name="Industrial Engineering",
        modules=[
            ("ie_production.svg", "Production &amp; Inventory",
             ["Deterministic Lot Sizing", "Stochastic Inventory", "Production Planning"]),
            ("ie_quality.svg", "Quality &amp; Reliability Control",
             ["Variables Control Charts", "Attributes Control Charts",
              "Process Capability", "Acceptance Sampling"]),
            ("ie_stochastic.svg", "Stochastic Operations",
             ["Queueing Systems", "Markov Chains", "Poisson Processes",
              "System Reliability"]),
        ],
    ),
    dict(
        slot="C",
        name="Mechanical Engineering",
        modules=[
            ("me_fluid.png", "Fluid Mechanics", ["Fluid Statics", "Fluid Kinematics"]),
            ("me_materials.png", "Mechanics of Materials", ["Torsion", "Stress &amp; Strain"]),
            ("me_vibrations.png", "Vibrations &amp; Acoustics",
             ["Harmonically Excited Vibrations", "Single Degree of Freedom Systems"]),
        ],
    ),
]

COL_W = (PANEL_W - 2 * PANEL_PAD_X - 2 * COL_GAP) // 3


# ------------------------------------------------------------ text metrics --
_fonts = {}


def _font(path, size):
    key = (path, size)
    if key not in _fonts:
        _fonts[key] = ImageFont.truetype(path, size)
    return _fonts[key]


def _plain(s):
    return s.replace("&amp;", "&")


def line_count(s, path, size, avail):
    """Greedy word wrap, matching how the browser breaks these labels.

    The 1.02 factor guards against Chrome laying a line out a hair wider than
    FreeType measures it, which would silently overflow a pill by one line.
    """
    font = _font(path, size)
    lines, cur = 1, ""
    for w in _plain(s).split():
        trial = "{} {}".format(cur, w).strip()
        if cur and font.getlength(trial) * 1.02 > avail:
            lines += 1
            cur = w
        else:
            cur = trial
    return lines


def modpill_height(label):
    avail = COL_W - MOD_ICON - MOD_ICON_GAP - 2 * TOPIC_PAD_X
    h = line_count(label, FONT_BOLD, FS_MODULE, avail) * MOD_LINE_H + 2 * MOD_PAD_Y
    return max(h, MOD_ICON + 2)  # the icon badge sets a floor on the pill height


def module_height(mod):
    _, label, topics = mod
    h = modpill_height(label) + MOD_GAP
    avail = COL_W - 2 * TOPIC_PAD_X - 2
    for i, t in enumerate(topics):
        if i:
            h += PILL_GAP
        h += line_count(t, FONT_REG, FS_TOPIC, avail) * TOPIC_LINE_H + 2 * TOPIC_PAD_Y
    return h


def panel_height(branch):
    return PANEL_PAD_T + max(module_height(m) for m in branch["modules"]) + PANEL_PAD_B


# ----------------------------------------------------------------- render ---
def module_html(mod):
    icon, label, topics = mod
    pills = "".join('<div class="topic">{}</div>'.format(t) for t in topics)
    return (
        '<div class="module">'
        '<div class="modpill" style="height:{h}px">'
        '<img class="modicon" src="assets/icons/{icon}" alt="">'
        '<span class="modlabel">{label}</span>'
        "</div>{pills}</div>"
    ).format(h=modpill_height(label), icon=icon, label=label, pills=pills)


def panel_html(number, branch, x, y, h):
    mods = "".join(module_html(m) for m in branch["modules"])
    return (
        '<div class="panel" style="left:{x}px;top:{y}px;width:{w}px;height:{h}px">'
        '<div class="branchpill"><span class="badge">{n}</span>{name}</div>'
        '<div class="cols">{mods}</div>'
        "</div>"
    ).format(x=x, y=y, w=PANEL_W, h=h, n=number, name=branch["name"], mods=mods)


def hub_html(cx, cy):
    circ = 2 * math.pi * HUB_RING_R
    period = circ / 3          # three dark arcs and three light ones, interleaved
    dash = period * 0.40       # the rest of each period is the white gap
    gap = period - dash
    d = HUB_R * 2
    return """
<div class="hub" style="left:{lx}px;top:{ly}px;width:{d}px;height:{d}px">
  <svg viewBox="0 0 {d} {d}" width="{d}" height="{d}">
    <circle cx="{r}" cy="{r}" r="{ro}" fill="#FFFFFF" stroke="{ink}" stroke-width="11"/>
    <circle cx="{r}" cy="{r}" r="{rr}" fill="none" stroke="{dark}"
            stroke-width="10" stroke-linecap="round"
            stroke-dasharray="{dash:.1f} {gap:.1f}" transform="rotate(-72 {r} {r})"/>
    <circle cx="{r}" cy="{r}" r="{rr}" fill="none" stroke="{light}"
            stroke-width="10" stroke-linecap="round"
            stroke-dasharray="{dash:.1f} {gap:.1f}" stroke-dashoffset="{off:.1f}"
            transform="rotate(-72 {r} {r})"/>
  </svg>
  <div class="hublabel">EngTrace</div>
</div>""".format(
        lx=cx - HUB_R, ly=cy - HUB_R, d=d, r=HUB_R, ro=HUB_R - 6, rr=HUB_RING_R,
        ink=C_INK, dark=C_HUB_DARK, light=C_HUB_LIGHT,
        dash=dash, gap=gap, off=-period / 2,
    )


def arrow_svg(points):
    """An axis-aligned arrow through `points`, with rounded bends.

    The last point is the edge the arrow points at: the head's tip stops
    ARROW_GAP short of it, and the shaft runs a little way under the head so
    no seam shows between them.
    """
    def unit(a, b):
        n = math.hypot(b[0] - a[0], b[1] - a[1])
        return (b[0] - a[0]) / n, (b[1] - a[1]) / n

    tx, ty = points[-1]
    ux, uy = unit(points[-2], points[-1])
    tipx, tipy = tx - ux * ARROW_GAP, ty - uy * ARROW_GAP
    basex, basey = tipx - ux * HEAD_L, tipy - uy * HEAD_L
    pts = points[:-1] + [(basex + ux * 2, basey + uy * 2)]

    d = ["M {:.1f} {:.1f}".format(*pts[0])]
    for a, (bx, by), c in zip(pts, pts[1:], pts[2:]):
        (ix, iy), (ox, oy) = unit(a, (bx, by)), unit((bx, by), c)
        sweep = 1 if ix * oy - iy * ox > 0 else 0   # clockwise on screen
        d.append("L {:.1f} {:.1f}".format(bx - ix * CORNER_R, by - iy * CORNER_R))
        d.append("A {r} {r} 0 0 {s} {:.1f} {:.1f}".format(
            bx + ox * CORNER_R, by + oy * CORNER_R, r=CORNER_R, s=sweep))
    d.append("L {:.1f} {:.1f}".format(*pts[-1]))

    nx, ny = -uy * HEAD_W / 2, ux * HEAD_W / 2
    head = "{:.1f},{:.1f} {:.1f},{:.1f} {:.1f},{:.1f}".format(
        tipx, tipy, basex + nx, basey + ny, basex - nx, basey - ny)
    return (
        '<path d="{d}" fill="none" stroke="{c}" stroke-width="{w}"/>'
        '<polygon points="{h}" fill="{c}"/>'
    ).format(d=" ".join(d), c=C_ARROW, w=ARROW_W, h=head)


def arrows_svg(routes, page_h):
    return """
<svg class="arrows" viewBox="0 0 {w} {h}" width="{w}" height="{h}">
  {p}
</svg>""".format(w=PAGE_W, h=page_h, p="".join(arrow_svg(r) for r in routes))


def layout():
    """Panel boxes, hub centre and arrow routes, shared with make_overview5_drawio.py."""
    by_slot = {b["slot"]: b for b in BRANCHES}
    h_a = max(panel_height(by_slot["A-left"]), panel_height(by_slot["A-right"]))
    h_b = max(panel_height(by_slot["B-left"]), panel_height(by_slot["B-right"]))
    h_c = panel_height(by_slot["C"])

    y_a = ROW_A_TOP
    y_b = y_a + h_a + ROW_GAP
    # Mechanical is shorter than its neighbours, so it sits on the row's
    # baseline; the space above it leaves room for the hub's downward arrow
    y_c = y_b + max(0, h_b - h_c)
    bottom = max(y_b + h_b, y_c + h_c)
    hub_cy = y_a + h_a / 2

    frame_h = bottom + FRAME_PAD_B - FRAME_Y
    page_h = FRAME_Y + frame_h + PAGE_PAD_B

    x_l, x_m, x_r = COL_X
    # (number, branch, x, y, height), numbered as in the original figure
    panels = [
        (1, by_slot["A-left"], x_l, y_a, h_a),
        (2, by_slot["A-right"], x_r, y_a, h_a),
        (3, by_slot["B-left"], x_l, y_b, h_b),
        (4, by_slot["B-right"], x_r, y_b, h_b),
        (5, by_slot["C"], x_m, y_c, h_c),
    ]

    # tails start under the hub ring, so they emerge cleanly from its edge
    tail = HUB_R - 8
    hub_l, hub_r = (CENTER_X - tail, hub_cy), (CENTER_X + tail, hub_cy)
    fork_l, fork_r = x_l + PANEL_W + FORK_DX, x_r - FORK_DX
    y_mid_b = y_b + h_b / 2
    routes = [
        [hub_l, (x_l + PANEL_W, hub_cy)],                                       # -> Chemical
        [hub_r, (x_r, hub_cy)],                                                 # -> Electrical
        [hub_l, (fork_l, hub_cy), (fork_l, y_mid_b), (x_l + PANEL_W, y_mid_b)],  # -> Civil
        [hub_r, (fork_r, hub_cy), (fork_r, y_mid_b), (x_r, y_mid_b)],            # -> Industrial
        [(CENTER_X, hub_cy + tail), (CENTER_X, y_c + PANEL_BORDER + PILL_TOP)],  # -> Mechanical
    ]
    return dict(panels=panels, routes=routes, hub_cy=hub_cy, frame_h=frame_h, page_h=page_h)


def build():
    lay = layout()
    page_h, frame_h, hub_cy = lay["page_h"], lay["frame_h"], lay["hub_cy"]
    panels = [panel_html(*p) for p in lay["panels"]]

    css = """
@page {{ size: {pw}px {ph}px; margin: 0; }}
* {{ box-sizing: border-box; }}
html, body {{ margin: 0; padding: 0; background: {page}; }}
body {{ font-family: "Times New Roman", Times, serif; color: {ink};
        -webkit-font-smoothing: antialiased; }}
.page {{ position: relative; width: {pw}px; height: {ph}px; overflow: hidden; }}
.frame {{ position: absolute; left: {fx}px; top: {fy}px; width: {fw}px; height: {fh}px;
          border: {fs}px solid {frame}; border-radius: {fr}px; }}
.titlepill {{ position: absolute; left: 50%; top: {fy}px; transform: translate(-50%, -50%);
              background: {title}; border-radius: 14px; padding: 9px 46px;
              font-weight: bold; font-size: {fst}px; line-height: 1.18; white-space: nowrap; }}
.panel {{ position: absolute; background: {panel}; border: {pb}px solid {frame};
          border-radius: 26px; }}
.branchpill {{ position: absolute; left: 50%; top: {pt}px; transform: translateX(-50%);
               background: {branch}; border-radius: 14px; padding: 6px 30px 8px 30px;
               font-weight: bold; font-size: {fsb}px; line-height: 1.22; white-space: nowrap; }}
.badge {{ position: absolute; left: -26px; top: -30px; width: 48px; height: 48px;
          border-radius: 50%; background: #FFFFFF; border: 2.5px solid {ink};
          font-size: {fsn}px; font-weight: bold; line-height: 43px; text-align: center; }}
/* the hatched drop shadow behind each circled branch number */
.badge::after {{ content: ""; position: absolute; left: 3px; top: 3px; z-index: -1;
                 width: 48px; height: 48px; border-radius: 50%;
                 background: repeating-linear-gradient(45deg, {ink} 0 1.2px,
                                                       transparent 1.2px 4.5px); }}
.cols {{ position: absolute; left: {ppx}px; top: {ppt}px; right: {ppx}px;
         display: flex; gap: {cgap}px; align-items: flex-start; }}
.module {{ width: {colw}px; }}
/* pills shrink to fit their label, as in the three-branch original, and wrap
   only once they would exceed the column */
.modpill {{ position: relative; background: {module}; box-shadow: 5px 5px 0 {modshadow};
            display: flex; align-items: center; gap: {igap}px;
            width: fit-content; max-width: 100%;
            padding: {mpy}px {tpx}px; margin-bottom: {mgap}px; }}
.modicon {{ width: {ic}px; height: {ic}px; flex: 0 0 auto; }}
.modlabel {{ font-weight: bold; font-size: {fsm}px; line-height: {mlh}px; }}
.topic {{ background: {topic}; border: 1px solid {topicedge};
          width: fit-content; max-width: 100%;
          padding: {tpy}px {tpx}px; margin-bottom: {pgap}px;
          font-size: {fstop}px; line-height: {tlh}px; }}
.topic:last-child {{ margin-bottom: 0; }}
.hub {{ position: absolute; }}
.hublabel {{ position: absolute; left: 0; top: 50%; width: 100%; transform: translateY(-50%);
             text-align: center; font-weight: bold; font-size: {fsh}px; line-height: 1; }}
.arrows {{ position: absolute; left: 0; top: 0; }}
""".format(
        pw=PAGE_W, ph=page_h, page=C_PAGE, ink=C_INK,
        fx=FRAME_X, fy=FRAME_Y, fw=FRAME_W, fh=frame_h, fs=FRAME_STROKE, fr=FRAME_R,
        frame=C_FRAME, title=C_TITLE, panel=C_PANEL, branch=C_BRANCH,
        module=C_MODULE, modshadow=C_MODULE_SHADOW, topic=C_TOPIC, topicedge=C_TOPIC_EDGE,
        fst=FS_TITLE, fsb=FS_BRANCH, fsn=FS_BADGE, fsm=FS_MODULE, fstop=FS_TOPIC, fsh=FS_HUB,
        ppx=PANEL_PAD_X, ppt=PANEL_PAD_T, cgap=COL_GAP, colw=COL_W, igap=MOD_ICON_GAP,
        mpy=MOD_PAD_Y, tpx=TOPIC_PAD_X, mgap=MOD_GAP, ic=MOD_ICON, mlh=MOD_LINE_H,
        tpy=TOPIC_PAD_Y, pgap=PILL_GAP, tlh=TOPIC_LINE_H, pb=PANEL_BORDER, pt=PILL_TOP,
    )

    return """<!DOCTYPE html>
<html><head><meta charset="utf-8"><title>EngTrace Taxonomy</title>
<style>{css}</style></head>
<body><div class="page">
  <div class="frame"></div>
  {arrows}
  {panels}
  {hub}
  <div class="titlepill">EngTrace Taxonomy</div>
</div></body></html>""".format(
        css=css,
        arrows=arrows_svg(lay["routes"], page_h),
        panels="".join(panels),
        hub=hub_html(CENTER_X, hub_cy),
    )


def main():
    HTML_OUT.write_text(build(), encoding="utf-8")
    print("wrote " + HTML_OUT.name)
    if not pathlib.Path(CHROME).exists():
        sys.exit("Chrome not found at " + CHROME)
    subprocess.run(
        [CHROME, "--headless", "--disable-gpu", "--no-pdf-header-footer",
         "--print-to-pdf=" + str(PDF_OUT), HTML_OUT.as_uri()],
        check=True, capture_output=True,
    )
    print("wrote " + PDF_OUT.name)


if __name__ == "__main__":
    main()
