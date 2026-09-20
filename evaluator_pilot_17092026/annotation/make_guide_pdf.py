"""Render guide.md to the PDF the annotators are sent.

    python evaluator_pilot_17092026/annotation/make_guide_pdf.py

guide.md is the source of truth; this only typesets it. The markdown subset used
there is headings (#, ##, ###), paragraphs, "-" bullets, ```code blocks```,
**bold**, *italic* and `inline code`.
"""
from __future__ import annotations

import os
import re

from reportlab.lib.enums import TA_JUSTIFY
from reportlab.lib.pagesizes import A4
from reportlab.lib.styles import ParagraphStyle
from reportlab.lib.units import mm
from reportlab.platypus import (BaseDocTemplate, Frame, ListFlowable, ListItem,
                                PageTemplate, Paragraph, Spacer)

HERE = os.path.dirname(os.path.abspath(__file__))
SRC = os.path.join(HERE, 'guide.md')
OUT = os.path.join(HERE, 'EngTrace-annotation-guide.pdf')

BODY = ParagraphStyle('body', fontName='Helvetica', fontSize=10.2, leading=15,
                      alignment=TA_JUSTIFY, spaceAfter=7)
H1 = ParagraphStyle('h1', fontName='Helvetica-Bold', fontSize=19, leading=23, spaceAfter=14)
H2 = ParagraphStyle('h2', fontName='Helvetica-Bold', fontSize=13, leading=17,
                    spaceBefore=15, spaceAfter=7, textColor='#1b3a5c')
H3 = ParagraphStyle('h3', fontName='Helvetica-Bold', fontSize=11, leading=15,
                    spaceBefore=10, spaceAfter=5)
BULLET = ParagraphStyle('bullet', parent=BODY, spaceAfter=4)
CODE = ParagraphStyle('code', fontName='Courier', fontSize=8.4, leading=11.4,
                      leftIndent=8, spaceBefore=4, spaceAfter=9,
                      backColor='#f4f5f7', borderPadding=6, borderColor='#dfe3e8',
                      borderWidth=0.5)


def inline(text: str) -> str:
    """The markdown subset -> reportlab markup, on already-escaped text."""
    text = text.replace('&', '&amp;').replace('<', '&lt;').replace('>', '&gt;')
    text = re.sub(r'`([^`]+)`', r'<font face="Courier">\1</font>', text)
    text = re.sub(r'\*\*([^*]+)\*\*', r'<b>\1</b>', text)
    text = re.sub(r'(?<![*\w])\*([^*]+)\*(?!\w)', r'<i>\1</i>', text)
    return text


def blocks(md: str):
    """Yield (kind, payload) for each block: h1/h2/h3, para, bullets."""
    lines = md.splitlines()
    i, n = 0, len(lines)
    while i < n:
        line = lines[i]
        if not line.strip():
            i += 1
        elif line.startswith('```'):
            i += 1
            code = []
            while i < n and not lines[i].startswith('```'):
                code.append(lines[i]); i += 1
            i += 1
            yield 'code', '\n'.join(code)
        elif line.startswith('### '):
            yield 'h3', line[4:]; i += 1
        elif line.startswith('## '):
            yield 'h2', line[3:]; i += 1
        elif line.startswith('# '):
            yield 'h1', line[2:]; i += 1
        elif line.lstrip().startswith('- '):
            items, indent = [], []
            while i < n and (lines[i].lstrip().startswith('- ') or
                             (lines[i].startswith('  ') and lines[i].strip() and items)):
                if lines[i].lstrip().startswith('- '):
                    items.append(lines[i].lstrip()[2:].strip())
                    indent.append(len(lines[i]) - len(lines[i].lstrip()))
                else:
                    items[-1] += ' ' + lines[i].strip()
                i += 1
            yield 'bullets', list(zip(items, indent))
        else:
            para = []
            while i < n and lines[i].strip() and not lines[i].startswith(('#', '- ', '```')):
                para.append(lines[i].strip()); i += 1
            yield 'para', ' '.join(para)


def story(md: str):
    out = []
    for kind, payload in blocks(md):
        if kind == 'h1':
            out.append(Paragraph(inline(payload), H1))
        elif kind == 'h2':
            out.append(Paragraph(inline(payload), H2))
        elif kind == 'h3':
            out.append(Paragraph(inline(payload), H3))
        elif kind == 'code':
            body = (payload.replace('&', '&amp;').replace('<', '&lt;').replace('>', '&gt;')
                    .replace(' ', '&nbsp;').replace('\n', '<br/>'))
            out.append(Paragraph(body, CODE))
        elif kind == 'para':
            out.append(Paragraph(inline(payload), BODY))
        else:
            out.append(ListFlowable(
                [ListItem(Paragraph(inline(t), BULLET), leftIndent=10 + ind * 1.2)
                 for t, ind in payload],
                bulletType='bullet', bulletFontSize=7, start='square',
                leftIndent=12, spaceAfter=7))
    return out


def footer(canvas, doc):
    canvas.saveState()
    canvas.setFont('Helvetica', 8)
    canvas.setFillColor('#5b6472')
    canvas.drawString(20 * mm, 12 * mm, 'EngTrace annotation guide')
    canvas.drawRightString(A4[0] - 20 * mm, 12 * mm, 'Page %d' % doc.page)
    canvas.restoreState()


def main():
    md = open(SRC, encoding='utf-8').read()
    doc = BaseDocTemplate(OUT, pagesize=A4, title='EngTrace annotation guide',
                          author='EngTrace', leftMargin=20 * mm, rightMargin=20 * mm,
                          topMargin=18 * mm, bottomMargin=20 * mm)
    frame = Frame(doc.leftMargin, doc.bottomMargin, doc.width, doc.height, id='body')
    doc.addPageTemplates([PageTemplate(id='page', frames=[frame], onPage=footer)])
    doc.build(story(md) + [Spacer(1, 4)])
    print('%s (%.0f KB)' % (OUT, os.path.getsize(OUT) / 1024))


if __name__ == '__main__':
    main()
