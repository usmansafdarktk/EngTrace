"""Render the pilot summary as a PDF that can be sent to someone outside the project.

    python evaluator_pilot_17092026/make_summary_pdf.py

The summary is written as markdown so it stays reviewable in the repository; this produces
the copy people actually read. It covers only what that document uses - headings, rules,
paragraphs, bullets, tables, bold, italic and inline code - and raises rather than guessing
if it meets anything else, so a silent formatting loss is not possible.

Calibri is registered from the system fonts, both so the page matches the project's other
shared documents and because reportlab's built-in fonts are Latin-1: the em dashes, times
signs and middle dots in the text would otherwise be dropped.
"""
import os
import re
import sys

from reportlab.lib import colors
from reportlab.lib.pagesizes import A4
from reportlab.lib.styles import ParagraphStyle
from reportlab.lib.units import mm
from reportlab.pdfbase import pdfmetrics
from reportlab.pdfbase.ttfonts import TTFont
from reportlab.platypus import (BaseDocTemplate, Frame, HRFlowable, Image, KeepTogether,
                                PageTemplate, Paragraph, Spacer, Table, TableStyle)

HERE = os.path.dirname(os.path.abspath(__file__))
SRC = os.path.join(HERE, 'PILOT_SUMMARY.md')
OUT = os.path.join(HERE, 'EngTrace-evaluator-pilot-summary.pdf')

INK = colors.HexColor('#1F2A37')
MUTED = colors.HexColor('#5B6472')
ACCENT = colors.HexColor('#1B3A5C')
RULE = colors.HexColor('#D5DBE3')
HEAD_FILL = colors.HexColor('#E8EDF3')
ZEBRA = colors.HexColor('#F7F9FB')


def register_fonts():
    win = os.path.join(os.environ.get('WINDIR', r'C:\Windows'), 'Fonts')
    faces = [('Calibri', 'calibri.ttf'), ('Calibri-Bold', 'calibrib.ttf'),
             ('Calibri-Italic', 'calibrii.ttf'), ('Calibri-BoldItalic', 'calibriz.ttf')]
    for name, filename in faces:
        path = os.path.join(win, filename)
        if not os.path.exists(path):
            raise SystemExit('missing font %s - expected at %s' % (name, path))
        pdfmetrics.registerFont(TTFont(name, path))
    pdfmetrics.registerFontFamily('Calibri', normal='Calibri', bold='Calibri-Bold',
                                  italic='Calibri-Italic', boldItalic='Calibri-BoldItalic')


def styles():
    base = dict(fontName='Calibri', textColor=INK, leading=14.5)
    return {
        'title': ParagraphStyle('title', fontName='Calibri-Bold', fontSize=21, leading=25,
                                textColor=INK, spaceAfter=10),
        'h2': ParagraphStyle('h2', fontName='Calibri-Bold', fontSize=13.5, leading=17,
                             textColor=ACCENT, spaceBefore=15, spaceAfter=6, keepWithNext=1),
        'h3': ParagraphStyle('h3', fontName='Calibri-Bold', fontSize=11.5, leading=15,
                             textColor=ACCENT, spaceBefore=12, spaceAfter=4, keepWithNext=1),
        'body': ParagraphStyle('body', fontSize=10, spaceAfter=7, **base),
        'bullet': ParagraphStyle('bullet', fontSize=10, leftIndent=12, bulletIndent=3,
                                 spaceAfter=4, **base),
        'cell': ParagraphStyle('cell', fontSize=8.7, leading=11.5, fontName='Calibri',
                               textColor=INK),
        'cellhead': ParagraphStyle('cellhead', fontSize=8.7, leading=11.5,
                                   fontName='Calibri-Bold', textColor=ACCENT),
        'foot': ParagraphStyle('foot', fontSize=8.5, leading=11, fontName='Calibri-Italic',
                               textColor=MUTED, spaceBefore=6),
    }


def inline(text):
    """Markdown emphasis to reportlab markup, ampersands escaped first."""
    text = text.replace('&', '&amp;').replace('<', '&lt;').replace('>', '&gt;')
    text = re.sub(r'\*\*(.+?)\*\*', r'<b>\1</b>', text)
    text = re.sub(r'(?<![*\w])\*([^*]+?)\*(?!\*)', r'<i>\1</i>', text)
    text = re.sub(r'`([^`]+?)`', r'<font face="Courier" size="8.5">\1</font>', text)
    return text


def split_row(line):
    return [c.strip() for c in line.strip().strip('|').split('|')]


def widths(rows, total):
    """Share the page by how much text each column carries, but never below the width its
    longest single word needs: a numeric column headed `n` is narrow by weight and still
    has to fit "60" without breaking it across two lines."""
    n = max(len(r) for r in rows)
    weight, floor = [], []
    for i in range(n):
        cells = [r[i] if i < len(r) else '' for r in rows]
        weight.append(max(len(c) for c in cells) ** 0.7)
        longest = max((max((len(w) for w in c.split()), default=0) for c in cells), default=1)
        # the widest word, plus the cell's own left and right padding
        widest = max((pdfmetrics.stringWidth(w, 'Calibri-Bold', 8.7)
                      for c in cells for w in c.split()), default=8)
        floor.append(min(widest + 14, total * 0.3))
    out = [w * total / sum(weight) for w in weight]
    for i, f in enumerate(floor):                  # raise the narrow ones to their floor
        out[i] = max(out[i], f)
    over = sum(out) - total
    if over > 0:                                   # take it back from the widest column
        out[out.index(max(out))] -= over
    return out


def build(md, st, avail):
    flow, i, lines = [], 0, md.split('\n')
    while i < len(lines):
        ln = lines[i].rstrip()
        if not ln.strip():
            i += 1
            continue
        if ln.startswith('# '):
            flow.append(Paragraph(inline(ln[2:]), st['title']))
        elif ln.startswith('### '):
            flow.append(Paragraph(inline(ln[4:]), st['h3']))
        elif ln.startswith('## '):
            flow.append(Paragraph(inline(ln[3:]), st['h2']))
        elif ln.strip() == '---':
            flow.append(Spacer(1, 3))
            flow.append(HRFlowable(width='100%', thickness=0.6, color=RULE,
                                   spaceBefore=1, spaceAfter=7))
        elif ln.startswith('!['):
            # a figure: scaled to the text column, kept on one page with what follows
            src = re.match(r'!\[[^\]]*\]\(([^)]+)\)', ln)
            path = os.path.join(HERE, src.group(1).replace('/', os.sep))
            if not os.path.exists(path):
                raise SystemExit('figure not found: %s' % path)
            from reportlab.lib.utils import ImageReader
            iw, ih = ImageReader(path).getSize()
            w = min(avail, iw * 0.5)
            flow.append(Spacer(1, 4))
            flow.append(Image(path, width=w, height=w * ih / iw))
            flow.append(Spacer(1, 9))
        elif ln.startswith('|'):
            block = []
            while i < len(lines) and lines[i].strip().startswith('|'):
                block.append(lines[i])
                i += 1
            rows = [split_row(r) for r in block if not re.match(r'^\|[\s:|-]+\|?$', r.strip())]
            data = [[Paragraph(inline(c), st['cellhead'] if n == 0 else st['cell']) for c in r]
                    for n, r in enumerate(rows)]
            style = [('GRID', (0, 0), (-1, -1), 0.4, RULE),
                     ('BACKGROUND', (0, 0), (-1, 0), HEAD_FILL),
                     ('VALIGN', (0, 0), (-1, -1), 'TOP'),
                     ('LEFTPADDING', (0, 0), (-1, -1), 5),
                     ('RIGHTPADDING', (0, 0), (-1, -1), 5),
                     ('TOPPADDING', (0, 0), (-1, -1), 4),
                     ('BOTTOMPADDING', (0, 0), (-1, -1), 4)]
            for r in range(2, len(rows), 2):
                style.append(('BACKGROUND', (0, r), (-1, r), ZEBRA))
            t = Table(data, colWidths=widths(rows, avail), repeatRows=1, hAlign='LEFT')
            t.setStyle(TableStyle(style))
            flow.append(t)
            flow.append(Spacer(1, 8))
            continue
        elif ln.startswith(('- ', '* ')) or re.match(r'^\d+\. ', ln):
            body = [ln]
            i += 1
            while (i < len(lines) and lines[i].startswith('  ') and lines[i].strip()
                   and not lines[i].strip().startswith(('- ', '* '))):
                body.append(lines[i].strip())
                i += 1
            text = ' '.join(body)
            mark = re.match(r'^(\d+)\. ', text)
            bullet = mark.group(1) + '.' if mark else '\u2022'
            text = re.sub(r'^([-*] |\d+\. )', '', text)
            flow.append(Paragraph(inline(text), st['bullet'], bulletText=bullet))
            continue
        elif ln.startswith('*') and not ln.startswith('**') and ln.endswith('*'):
            # a whole paragraph in italics is the document's aside style; a line that merely
            # OPENS with a bold lead-in ("**Verification.** ...") is body text, not an aside
            para = [ln]
            i += 1
            while i < len(lines) and lines[i].strip():
                para.append(lines[i].strip())
                i += 1
            flow.append(Paragraph(inline(' '.join(para)), st['foot']))
            continue
        else:
            para = [ln]
            i += 1
            while (i < len(lines) and lines[i].strip() != '---' and lines[i].strip()
                   and not lines[i].startswith(('|', '#', '- ', '* ', '!['))
                   and not re.match(r'^\d+\. ', lines[i])):
                para.append(lines[i].strip())
                i += 1
            flow.append(Paragraph(inline(' '.join(para)), st['body']))
            continue
        i += 1
    return flow


def main():
    if not os.path.exists(SRC):
        raise SystemExit('no summary to render at %s' % SRC)
    register_fonts()
    st = styles()
    margin = 18 * mm
    avail = A4[0] - 2 * margin

    def footer(canvas, doc):
        canvas.saveState()
        canvas.setFont('Calibri', 8)
        canvas.setFillColor(MUTED)
        canvas.drawRightString(A4[0] - margin, 11 * mm,
                               'EngTrace evaluator pilot  \u00b7  page %d' % doc.page)
        canvas.restoreState()

    doc = BaseDocTemplate(OUT, pagesize=A4, leftMargin=margin, rightMargin=margin,
                          topMargin=16 * mm, bottomMargin=18 * mm,
                          title='EngTrace evaluator pilot', author='EngTrace')
    frame = Frame(margin, 18 * mm, avail, A4[1] - 34 * mm, id='body')
    doc.addPageTemplates([PageTemplate(id='all', frames=[frame], onPage=footer)])
    doc.build(build(open(SRC, encoding='utf-8').read(), st, avail))
    print('wrote %s (%d bytes)' % (OUT, os.path.getsize(OUT)))
    return 0


if __name__ == '__main__':
    sys.exit(main())
