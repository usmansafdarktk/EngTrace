"""Write the five-branch overview as an editable draw.io file.

Every piece is its own shape on one flat layer: the frame, the five panels,
branch pills, number badges and their hatched shadows, module pills, icons,
topic pills, the title, the hub ring, its six arcs and its label, and the five
arrows. Nothing is grouped, nested in a container or locked, so each shape can
be clicked and moved on its own. The frame ignores clicks inside it
(pointerEvents=0), so dragging on empty space starts a selection box instead
of moving the background, and the diagram has no page boundary.

Box positions and line breaks are measured from the HTML that make_overview5.py
renders (headless Chrome), so the draw.io file matches the PDF. The arrows are
connected to the hub ring and to their panels, so they follow when those move.

Usage:  python figures_oct_12/make_overview5_drawio.py
"""

import base64
import html
import json
import pathlib
import re
import subprocess
import sys
import xml.etree.ElementTree as ET

import make_overview5 as m

OUT = m.HERE / "engtrace-overview-5branch.drawio"
MEASURE_HTML = m.HERE / "_measure_tmp.html"

UNLOCKED = "movable=1;resizable=1;rotatable=1;deletable=1;editable=1;locked=0;connectable=1;"
FONT = "fontFamily=Times New Roman;fontColor={};".format(m.C_INK)

# Collects every box the figure draws, with the line breaks Chrome chose, as JSON.
MEASURE_JS = r"""
<script>
window.addEventListener('load', () => {
  const box = el => { const r = el.getBoundingClientRect(); return [r.left, r.top, r.width, r.height]; };
  function lines(el, ownTextOnly) {
    const out = []; let cur = [], top = null;
    const walk = document.createTreeWalker(el, NodeFilter.SHOW_TEXT); let n;
    while ((n = walk.nextNode())) {
      if (ownTextOnly && n.parentElement !== el) continue;
      const re = /\S+/g; let w;
      while ((w = re.exec(n.data))) {
        const r = document.createRange(); r.setStart(n, w.index); r.setEnd(n, w.index + w[0].length);
        const t = r.getClientRects()[0].top;
        if (top !== null && Math.abs(t - top) > 3) { out.push(cur.join(' ')); cur = []; }
        top = t; cur.push(w[0]);
      }
    }
    if (cur.length) out.push(cur.join(' '));
    return out;
  }
  const data = {
    frame: box(document.querySelector('.frame')),
    title: {box: box(document.querySelector('.titlepill')), text: document.querySelector('.titlepill').textContent},
    hub: box(document.querySelector('.hub')),
    hublabel: {box: box(document.querySelector('.hublabel')), text: document.querySelector('.hublabel').textContent},
    panels: [...document.querySelectorAll('.panel')].map(p => {
      const pill = p.querySelector('.branchpill'), badge = p.querySelector('.badge');
      return {
        box: box(p),
        pill: {box: box(pill), lines: lines(pill, true)},
        badge: {box: box(badge), text: badge.textContent},
        modules: [...p.querySelectorAll('.module')].map(mod => ({
          pill: {box: box(mod.querySelector('.modpill')), lines: lines(mod.querySelector('.modlabel'))},
          icon: {box: box(mod.querySelector('.modicon')), src: mod.querySelector('.modicon').getAttribute('src')},
          topics: [...mod.querySelectorAll('.topic')].map(t => ({box: box(t), lines: lines(t)})),
        })),
      };
    }),
  };
  const pre = document.createElement('pre'); pre.id = 'MEASURE';
  pre.textContent = JSON.stringify(data); document.body.appendChild(pre);
});
</script>"""

# The hatched drop shadow behind each number badge, as in the CSS
# repeating-linear-gradient(45deg, ink 0 1.2px, transparent 1.2px 4.5px).
HATCH_SVG = (
    '<svg xmlns="http://www.w3.org/2000/svg" width="48" height="48" viewBox="0 0 48 48">'
    '<defs><pattern id="h" width="4.5" height="4.5" patternUnits="userSpaceOnUse" '
    'patternTransform="rotate(-45)"><rect width="1.2" height="4.5" fill="{}"/></pattern></defs>'
    '<circle cx="24" cy="24" r="24" fill="url(#h)"/></svg>'
).format(m.C_INK)


def measure():
    page = m.build().replace("</body>", MEASURE_JS + "</body>")
    MEASURE_HTML.write_text(page, encoding="utf-8")
    try:
        dom = subprocess.run(
            [m.CHROME, "--headless", "--disable-gpu", "--window-size=4200,1400",
             "--virtual-time-budget=5000", "--dump-dom", MEASURE_HTML.as_uri()],
            capture_output=True, text=True, encoding="utf-8", check=True,
        ).stdout
    finally:
        MEASURE_HTML.unlink()
    found = re.search(r'<pre id="MEASURE">(.*?)</pre>', dom, re.S)
    if not found:
        sys.exit("Chrome returned no measurements")
    return json.loads(html.unescape(found.group(1)))


def data_uri(path):
    # draw.io styles are ';'-separated, so embedded images drop the ';base64' marker
    kind = "svg+xml" if path.suffix == ".svg" else "png"
    return "data:image/{},{}".format(kind, base64.b64encode(path.read_bytes()).decode("ascii"))


def label(lines, line_h):
    return '<div style="line-height:{}">{}</div>'.format(
        line_h, "<br>".join(html.escape(s) for s in lines))


class Diagram:
    def __init__(self):
        self.root = ET.Element("root")
        ET.SubElement(self.root, "mxCell", id="0")
        ET.SubElement(self.root, "mxCell", id="1", parent="0")

    def vertex(self, cid, style, box, inset=0.0, value=""):
        """A shape whose outer edge matches `box`; `inset` is half its stroke width."""
        x, y, w, h = box
        cell = ET.SubElement(self.root, "mxCell", id=cid, value=value,
                             style=style + UNLOCKED, vertex="1", parent="1")
        ET.SubElement(cell, "mxGeometry", {
            "x": "{:.2f}".format(x + inset), "y": "{:.2f}".format(y + inset),
            "width": "{:.2f}".format(w - 2 * inset), "height": "{:.2f}".format(h - 2 * inset),
            "as": "geometry"})

    def edge(self, cid, source, target, style, waypoints):
        cell = ET.SubElement(self.root, "mxCell", id=cid, value="", style=style,
                             edge="1", parent="1", source=source, target=target)
        geo = ET.SubElement(cell, "mxGeometry", relative="1", **{"as": "geometry"})
        if waypoints:
            pts = ET.SubElement(geo, "Array", **{"as": "points"})
            for x, y in waypoints:
                ET.SubElement(pts, "mxPoint", x="{:.2f}".format(x), y="{:.2f}".format(y))

    def write(self, path, page_h):
        model = ET.Element("mxGraphModel", {
            "grid": "1", "gridSize": "10", "guides": "1", "tooltips": "1", "connect": "1",
            "arrows": "1", "fold": "1", "page": "0", "pageScale": "1",
            "pageWidth": str(m.PAGE_W), "pageHeight": str(int(page_h)), "math": "0", "shadow": "0"})
        model.append(self.root)
        mxfile = ET.Element("mxfile", host="app.diagrams.net")
        diagram = ET.SubElement(mxfile, "diagram", name="EngChain Taxonomy", id="engchain-taxonomy")
        diagram.append(model)
        ET.indent(mxfile)
        ET.ElementTree(mxfile).write(path, encoding="utf-8", xml_declaration=False)


def hub_arcs():
    """The six ring arcs as draw.io arc angles (fractions of a turn, clockwise from 12 o'clock).

    Mirrors hub_html(): three dark and three light dashes, each 40% of a
    120-degree period, the dark set rotated -72 degrees and the light set
    offset by half a period.
    """
    period, dash = 120.0, 48.0
    arcs = []
    for color, offset in ((m.C_HUB_DARK, 0.0), (m.C_HUB_LIGHT, period / 2)):
        for k in range(3):
            start = -72 + offset + k * period + 90   # from 3 o'clock to 12 o'clock
            arcs.append((color, (start % 360) / 360, ((start + dash) % 360) / 360))
    return arcs


def main():
    lay = m.layout()
    got = measure()
    d = Diagram()

    d.vertex("frame", "rounded=1;absoluteArcSize=1;arcSize={};whiteSpace=wrap;html=1;"
             "fillColor=none;strokeColor={};strokeWidth={};pointerEvents=0;".format(
                 2 * (m.FRAME_R - m.FRAME_STROKE / 2), m.C_FRAME, m.FRAME_STROKE),
             got["frame"], inset=m.FRAME_STROKE / 2)

    # Arrows come before the hub in z-order, so their tails hide under the ring.
    gap = m.PANEL_BORDER / 2 + m.ARROW_GAP
    arrow = ("edgeStyle=none;html=1;rounded=1;arcSize={};endArrow=block;endFill=1;endSize=18;"
             "startArrow=none;strokeColor={};strokeWidth={};"
             "exitX={};exitY={};exitDx=0;exitDy=0;exitPerimeter=0;"
             "entryX={};entryY={};entryDx={};entryDy={};entryPerimeter=0;")
    ends = [  # (target, exit x/y on the ring, entry x/y/dx/dy on the target)
        ("panel-1", (0, 0.5), (1, 0.5, gap, 0)),
        ("panel-2", (1, 0.5), (0, 0.5, -gap, 0)),
        ("panel-3", (0, 0.5), (1, 0.5, gap, 0)),
        ("panel-4", (1, 0.5), (0, 0.5, -gap, 0)),
        ("pill-5", (0.5, 1), (0.5, 0, 0, -m.ARROW_GAP)),
    ]
    for (target, ex, en), route in zip(ends, lay["routes"]):
        d.edge("arrow-" + target.split("-")[1], "hub-ring", target,
               arrow.format(2 * m.CORNER_R, m.C_ARROW, m.ARROW_W, *ex, *en), route[1:-1])

    hatch = "data:image/svg+xml," + base64.b64encode(HATCH_SVG.encode()).decode("ascii")
    for n, p in enumerate(got["panels"], 1):
        d.vertex("panel-{}".format(n), "rounded=1;absoluteArcSize=1;arcSize={};whiteSpace=wrap;html=1;"
                 "fillColor={};strokeColor={};strokeWidth={};".format(
                     2 * (26 - m.PANEL_BORDER / 2), m.C_PANEL, m.C_FRAME, m.PANEL_BORDER),
                 p["box"], inset=m.PANEL_BORDER / 2)
        d.vertex("pill-{}".format(n), "rounded=1;absoluteArcSize=1;arcSize=28;whiteSpace=wrap;html=1;"
                 "fillColor={};strokeColor=none;fontStyle=1;fontSize={};align=center;verticalAlign=middle;"
                 "spacing=0;spacingLeft=28;spacingRight=28;spacingTop=6;spacingBottom=8;".format(
                     m.C_BRANCH, m.FS_BRANCH) + FONT,
                 p["pill"]["box"], value=label(p["pill"]["lines"], 1.22))
        bx, by, bw, bh = p["badge"]["box"]
        border = 2.5
        d.vertex("badge-shadow-{}".format(n), "shape=image;html=1;imageAspect=0;aspect=fixed;"
                 "image={};".format(hatch), (bx + border + 3, by + border + 3, 48, 48))
        d.vertex("badge-{}".format(n), "ellipse;whiteSpace=wrap;html=1;aspect=fixed;fillColor=#FFFFFF;"
                 "strokeColor={};strokeWidth={};fontStyle=1;fontSize={};align=center;"
                 "verticalAlign=middle;spacing=0;".format(m.C_INK, border, m.FS_BADGE) + FONT,
                 p["badge"]["box"], inset=border / 2, value=html.escape(p["badge"]["text"]))

        for j, mod in enumerate(p["modules"], 1):
            tag = "{}-{}".format(n, j)
            # spacingRight is 2px under the CSS padding, so draw.io never re-wraps a label
            d.vertex("module-" + tag, "rounded=0;whiteSpace=wrap;html=1;fillColor={};strokeColor=none;"
                     "shadow=1;shadowColor={};shadowOffsetX=5;shadowOffsetY=5;shadowBlur=0;shadowOpacity=100;"
                     "fontStyle=1;fontSize={};align=left;verticalAlign=middle;spacing=0;"
                     "spacingLeft={};spacingRight={};".format(
                         m.C_MODULE, m.C_MODULE_SHADOW, m.FS_MODULE,
                         m.TOPIC_PAD_X + m.MOD_ICON + m.MOD_ICON_GAP, m.TOPIC_PAD_X - 2) + FONT,
                     mod["pill"]["box"], value=label(mod["pill"]["lines"], "{}px".format(m.MOD_LINE_H)))
            d.vertex("icon-" + tag, "shape=image;html=1;imageAspect=0;aspect=fixed;image={};".format(
                         data_uri(m.HERE / mod["icon"]["src"])), mod["icon"]["box"])
            for k, t in enumerate(mod["topics"], 1):
                d.vertex("topic-{}-{}".format(tag, k),
                         "rounded=0;whiteSpace=wrap;html=1;fillColor={};strokeColor={};strokeWidth=1;"
                         "fontSize={};align=left;verticalAlign=middle;spacing=0;spacingLeft={};"
                         "spacingRight={};spacingTop={};spacingBottom={};".format(
                             m.C_TOPIC, m.C_TOPIC_EDGE, m.FS_TOPIC, m.TOPIC_PAD_X, m.TOPIC_PAD_X - 2,
                             m.TOPIC_PAD_Y, m.TOPIC_PAD_Y) + FONT,
                         t["box"], inset=0.5, value=label(t["lines"], "{}px".format(m.TOPIC_LINE_H)))

    hx, hy, hw, hh = got["hub"]
    cx, cy = hx + hw / 2, hy + hh / 2
    ring = m.HUB_R - 6
    d.vertex("hub-ring", "ellipse;whiteSpace=wrap;html=1;aspect=fixed;fillColor=#FFFFFF;"
             "strokeColor={};strokeWidth=11;".format(m.C_INK),
             (cx - ring, cy - ring, 2 * ring, 2 * ring))
    r = m.HUB_RING_R
    for i, (color, start, end) in enumerate(hub_arcs(), 1):
        d.vertex("hub-arc-{}".format(i), "shape=mxgraph.basic.arc;html=1;aspect=fixed;fillColor=none;"
                 "strokeColor={};strokeWidth=10;linecap=round;pointerEvents=0;"
                 "startAngle={:.4f};endAngle={:.4f};".format(color, start, end),
                 (cx - r, cy - r, 2 * r, 2 * r))
    d.vertex("hub-label", "text;html=1;align=center;verticalAlign=middle;fontStyle=1;fontSize={};"
             "spacing=0;".format(m.FS_HUB) + FONT,
             got["hublabel"]["box"], value=label([got["hublabel"]["text"]], 1))

    d.vertex("title", "rounded=1;absoluteArcSize=1;arcSize=28;whiteSpace=wrap;html=1;fillColor={};"
             "strokeColor=none;fontStyle=1;fontSize={};align=center;verticalAlign=middle;"
             "spacing=0;spacingLeft=44;spacingRight=44;spacingTop=9;spacingBottom=9;".format(
                 m.C_TITLE, m.FS_TITLE) + FONT,
             got["title"]["box"], value=label([got["title"]["text"]], 1.18))

    d.write(OUT, lay["page_h"])
    print("wrote {} ({} shapes)".format(OUT.name, len(d.root) - 2))


if __name__ == "__main__":
    main()
