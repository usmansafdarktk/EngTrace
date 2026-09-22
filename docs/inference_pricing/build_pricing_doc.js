// Build EngTrace-inference-pricing.docx: model prices and what 2,250 items would cost.
//
//   npm install docx            (once, anywhere on NODE_PATH)
//   node build_pricing_doc.js   -> EngTrace-inference-pricing.docx beside this script
//
// Prices are per million tokens, read on OpenRouter 19-21 September 2026. The cost column
// assumes every model writes as much as GPT-5 did on the evaluator pilot - 270 input and
// 5,232 output tokens per problem, reasoning included (evaluator_pilot_17092026/analysis/
// inference_cost.py) - over 150 templates x 15 instances = 2,250 runs per model, no retries.

const fs = require('fs');
const path = require('path');
const {
  Document, Packer, Paragraph, TextRun, Table, TableRow, TableCell, WidthType,
  ShadingType, AlignmentType, HeadingLevel, BorderStyle, LevelFormat, Footer, PageNumber,
} = require('docx');

const ITEMS = 2250;
const REASON = { in: 270, out: 5232 };

const GROUPS = [
  ['Open-weight models, through OpenRouter', 'Cheapest upstream serving fp8 or better.', [
    ['openai/gpt-oss-20b', 0.018, 0.090], ['google/gemma-4-26b-a4b-it', 0.090, 0.300],
    ['google/gemma-4-31b-it', 0.130, 0.380], ['deepseek/deepseek-v4-flash', 0.036, 0.071],
    ['deepseek/deepseek-v4.1-flash', 0.140, 0.420], ['qwen/qwen3-235b-a22b', 0.087, 0.350],
    ['qwen/qwen3.8-27b', 0.150, 1.875], ['z-ai/glm-5.3-flash', 0.075, 0.250],
    ['z-ai/glm-5.3', 1.400, 4.400], ['meta/muse-glimmer-30b', 0.300, 1.100],
    ['moonshotai/kimi-k3', 2.100, 10.950]]],
  ['OpenAI', '', [
    ['gpt-5-nano †', 0.05, 0.40], ['gpt-5.6-luna', 0.20, 1.20], ['gpt-5.4-nano', 0.20, 1.25],
    ['gpt-5-mini †', 0.25, 2.00], ['gpt-5.4-mini', 0.75, 4.50], ['gpt-5 †', 1.25, 10.00],
    ['gpt-5.6-terra', 2.00, 12.00], ['gpt-5.6-sol', 4.00, 20.00], ['gpt-5.5', 5.00, 30.00],
    ['gpt-6-astra', 10.00, 50.00]]],
  ['Google', '', [
    ['gemini-2.5-flash-lite †', 0.10, 0.40], ['gemini-3.1-flash-lite', 0.25, 1.50],
    ['gemini-2.5-flash †', 0.30, 2.50], ['gemini-3.6 / 3.7 / 3.8-flash', 0.75, 3.75],
    ['gemini-2.5-pro †', 1.25, 10.00], ['gemini-3.5-flash', 1.50, 9.00],
    ['gemini-3.1-pro-preview', 2.00, 12.00]]],
  ['Anthropic', '', [
    ['claude-haiku-4.5', 1.00, 5.00], ['claude-sonnet-5', 2.00, 10.00],
    ['claude-opus-5', 5.00, 25.00], ['claude-fable-5-1', 10.00, 50.00]]],
];

const cost = (p, t) => ITEMS * (t.in * p[1] + t.out * p[2]) / 1e6;
const money = (v) => '$' + (v < 10 ? v.toFixed(2) : v.toFixed(0));
// Prices exactly as quoted: three decimals for the open-weight list, two elsewhere.
const price = (v, dp) => v.toFixed(dp);

const FONT = 'Calibri';
const INK = '1F2A37';
const MUTED = '5B6472';
const ACCENT = '1B3A5C';
const HEAD_FILL = 'E8EDF3';
const ZEBRA = 'F7F9FB';
const W = [4560, 1500, 1500, 1800];                 // sums to 9360 = US Letter minus 1" margins
const border = { style: BorderStyle.SINGLE, size: 4, color: 'D5DBE3' };
const borders = { top: border, bottom: border, left: border, right: border };

const run = (text, o = {}) => new TextRun({ text, font: FONT, size: o.size || 20, bold: o.bold,
  italics: o.italics, color: o.color || INK });
const para = (children, o = {}) => new Paragraph({ children, spacing: { after: o.after ?? 120, before: o.before ?? 0 },
  alignment: o.align, keepNext: o.keepNext });

function cell(text, i, o = {}) {
  // keepNext on every row but the last keeps a table on one page (Word honours it per paragraph)
  return new TableCell({
    borders, width: { size: W[i], type: WidthType.DXA },
    shading: o.fill ? { fill: o.fill, type: ShadingType.CLEAR, color: 'auto' } : undefined,
    margins: { top: 48, bottom: 48, left: 110, right: 110 },
    children: [new Paragraph({ alignment: i === 0 ? AlignmentType.LEFT : AlignmentType.RIGHT, keepNext: o.keepNext,
      children: [run(text, { bold: o.bold, size: o.size || 19, color: o.color })] })],
  });
}

function table(rows, dp) {
  const head = new TableRow({ tableHeader: true, cantSplit: true, children: [
    'Model', 'Input $/M', 'Output $/M', 'Est. cost',
  ].map((h, i) => cell(h, i, { bold: true, fill: HEAD_FILL, size: 18, color: ACCENT, keepNext: true })) });
  const body = rows.map((r, k) => {
    const o = { fill: k % 2 ? ZEBRA : undefined, keepNext: k < rows.length - 1 };
    return new TableRow({ cantSplit: true, children: [
      cell(r[0], 0, o), cell(price(r[1], dp), 1, o), cell(price(r[2], dp), 2, o),
      cell(money(cost(r, REASON)), 3, { ...o, bold: true }),
    ] });
  });
  return new Table({ width: { size: 9360, type: WidthType.DXA }, columnWidths: W, rows: [head, ...body] });
}

const bullet = (children) => new Paragraph({ numbering: { reference: 'dots', level: 0 }, children,
  spacing: { after: 90 } });

const children = [
  new Paragraph({ heading: HeadingLevel.TITLE, spacing: { after: 60 },
    children: [new TextRun({ text: 'EngTrace inference pricing', font: FONT, size: 44, bold: true, color: INK })] }),
  para([run('Inference is each model solving every EngTrace problem once: 150 templates × 15 instances = 2,250 runs per model.', { color: MUTED, size: 21 })], { after: 200 }),
];

for (const [title, note, rows] of GROUPS) {
  const dp = title.startsWith('Open-weight') ? 3 : 2;
  children.push(new Paragraph({ heading: HeadingLevel.HEADING_2, keepNext: true, spacing: { before: 300, after: 60 },
    children: [new TextRun({ text: title, font: FONT, size: 26, bold: true, color: ACCENT })] }));
  if (note) children.push(para([run(note, { color: MUTED, size: 18, italics: true })], { after: 100, keepNext: true }));
  children.push(table(rows, dp));
}

children.push(
  para([run('† Retires before publication: gpt-5, gpt-5-mini and gpt-5-nano on 11 December 2026; the Gemini 2.5 models from 20 October 2026.', { size: 18, color: MUTED })], { before: 200 }),
  new Paragraph({ heading: HeadingLevel.HEADING_2, spacing: { before: 240, after: 100 },
    children: [new TextRun({ text: 'How the estimates are made', font: FONT, size: 26, bold: true, color: ACCENT })] }),
  bullet([run('Est. cost = 2,250 × (input tokens × input price + output tokens × output price) ÷ 1,000,000, the full cost of all 2,250 runs.')]),
  bullet([run('It assumes each model writes as much as GPT-5 did in EngTrace’s evaluator pilot: about 270 input and 5,232 output tokens per problem, reasoning included. A model that answers without reasoning writes about 800 output tokens and would cost roughly a sixth of the figure shown.')]),
  bullet([run('Not included: re-runs of truncated or failed answers (0–15% extra for most models in the pilot) and half-price batch endpoints.')]),
  bullet([run('Prices as read on OpenRouter, 19–21 September 2026; they move between upstreams from day to day.')]),
);

const doc = new Document({
  creator: 'EngTrace', title: 'EngTrace inference pricing',
  styles: { default: { document: { run: { font: FONT, size: 20 } } } },
  numbering: { config: [{ reference: 'dots', levels: [{ level: 0, format: LevelFormat.BULLET, text: '•',
    alignment: AlignmentType.LEFT, style: { paragraph: { indent: { left: 460, hanging: 260 } } } }] }] },
  sections: [{
    properties: { page: { size: { width: 12240, height: 15840 }, margin: { top: 1300, bottom: 1300, left: 1440, right: 1440 } } },
    footers: { default: new Footer({ children: [new Paragraph({ alignment: AlignmentType.RIGHT, children: [
      new TextRun({ children: ['EngTrace inference pricing · page ', PageNumber.CURRENT], font: FONT, size: 16, color: MUTED })] })] }) },
    children,
  }],
});

const out = path.join(__dirname, 'EngTrace-inference-pricing.docx');
Packer.toBuffer(doc).then((buf) => { fs.writeFileSync(out, buf); console.log('wrote', out, buf.length, 'bytes'); });
