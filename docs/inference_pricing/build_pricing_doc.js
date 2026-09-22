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
  ShadingType, AlignmentType, HeadingLevel, BorderStyle, Footer, PageNumber,
} = require('docx');

const ITEMS = 2250;
const REASON = { in: 270, out: 5232 };

// The roster, chosen 2026-09-22. Excluded: models the evaluator pilot generated with
// (gpt-5, claude-opus-4.7, gemini-3.1-pro, deepseek-r1, llama-3.1-70b, and the robustness
// pair gemma-4-31b-it and qwen3.8-27b) and models used as judges (gpt-5, claude-opus-4.5,
// gemini-3.1-pro, grok, minimax, mimo). deepseek-v4-flash dropped for v4.1-flash; the
// frontier list trimmed to one or two models per provider.
const GROUPS = [
  ['Open-source models', 'Through OpenRouter, cheapest upstream serving fp8 or better.', 3, [
    ['openai/gpt-oss-20b', 0.018, 0.090], ['google/gemma-4-26b-a4b-it', 0.090, 0.300],
    ['deepseek/deepseek-v4.1-flash', 0.140, 0.420], ['qwen/qwen3-235b-a22b', 0.087, 0.350],
    ['z-ai/glm-5.3-flash', 0.075, 0.250], ['z-ai/glm-5.3', 1.400, 4.400],
    ['meta/muse-glimmer-30b', 0.300, 1.100], ['moonshotai/kimi-k3', 2.100, 10.950]]],
  ['Frontier models', '', 2, [
    ['openai/gpt-5.4-nano', 0.20, 1.25], ['openai/gpt-5.4-mini', 0.75, 4.50],
    ['google/gemini-3.1-flash-lite', 0.25, 1.50],
    ['anthropic/claude-haiku-4.5', 1.00, 5.00], ['anthropic/claude-sonnet-5', 2.00, 10.00]]],
];

const cost = (p, t) => ITEMS * (t.in * p[1] + t.out * p[2]) / 1e6;
// To the cent throughout, so each table's rows add up to its total exactly.
const money = (v) => '$' + v.toFixed(2);
// Prices exactly as quoted: three decimals for the open-source list, two for frontier.
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
  return new TableCell({
    borders, width: { size: W[i], type: WidthType.DXA },
    shading: o.fill ? { fill: o.fill, type: ShadingType.CLEAR, color: 'auto' } : undefined,
    margins: { top: 48, bottom: 48, left: 110, right: 110 },
    children: [new Paragraph({ alignment: i === 0 ? AlignmentType.LEFT : AlignmentType.RIGHT, keepNext: o.keepNext,
      children: [run(text, { bold: o.bold, size: o.size || 19, color: o.color })] })],
  });
}

// keepNext on the header and every model row keeps each table, total included, on one page.
function table(rows, dp) {
  const head = new TableRow({ tableHeader: true, cantSplit: true, children: [
    'Model', 'Input $/M', 'Output $/M', 'Est. cost',
  ].map((h, i) => cell(h, i, { bold: true, fill: HEAD_FILL, size: 18, color: ACCENT, keepNext: true })) });
  const body = rows.map((r, k) => {
    const o = { fill: k % 2 ? ZEBRA : undefined, keepNext: true };
    return new TableRow({ cantSplit: true, children: [
      cell(r[0], 0, o), cell(price(r[1], dp), 1, o), cell(price(r[2], dp), 2, o),
      cell(money(cost(r, REASON)), 3, { ...o, bold: true }),
    ] });
  });
  const total = rows.reduce((t, r) => t + cost(r, REASON), 0);
  const foot = new TableRow({ cantSplit: true, children: [
    cell('Total', 0, { bold: true, fill: HEAD_FILL, color: ACCENT }), cell('', 1, { fill: HEAD_FILL }),
    cell('', 2, { fill: HEAD_FILL }), cell(money(total), 3, { bold: true, fill: HEAD_FILL, color: ACCENT }),
  ] });
  return new Table({ width: { size: 9360, type: WidthType.DXA }, columnWidths: W, rows: [head, ...body, foot] });
}

const children = [
  new Paragraph({ heading: HeadingLevel.TITLE, spacing: { after: 60 },
    children: [new TextRun({ text: 'EngTrace inference pricing', font: FONT, size: 44, bold: true, color: INK })] }),
  para([run('Inference is each model solving every EngTrace problem once; evaluating the answers is not included. '
    + 'EngTrace covers five engineering branches (chemical, civil, electrical, industrial and mechanical), each with '
    + '30 problem templates. Every template is instantiated 15 times with different sampled values, giving '
    + '5 × 30 × 15 = 2,250 problems (870 easy, 870 intermediate, 510 advanced), so each model makes 2,250 runs. '
    + 'The estimated cost is the total over those runs, assuming each model writes as much as GPT-5 did in '
    + 'EngTrace’s evaluator pilot: about 270 input and 5,232 output tokens per problem, reasoning included.',
    { color: MUTED, size: 21 })], { after: 200 }),
];

for (const [title, note, dp, rows] of GROUPS) {
  children.push(new Paragraph({ heading: HeadingLevel.HEADING_2, keepNext: true, spacing: { before: 300, after: 60 },
    children: [new TextRun({ text: title, font: FONT, size: 26, bold: true, color: ACCENT })] }));
  if (note) children.push(para([run(note, { color: MUTED, size: 18, italics: true })], { after: 100, keepNext: true }));
  children.push(table(rows, dp));
}

const doc = new Document({
  creator: 'EngTrace', title: 'EngTrace inference pricing',
  styles: { default: { document: { run: { font: FONT, size: 20 } } } },
  sections: [{
    properties: { page: { size: { width: 12240, height: 15840 }, margin: { top: 1300, bottom: 1300, left: 1440, right: 1440 } } },
    footers: { default: new Footer({ children: [new Paragraph({ alignment: AlignmentType.RIGHT, children: [
      new TextRun({ children: ['EngTrace inference pricing · page ', PageNumber.CURRENT], font: FONT, size: 16, color: MUTED })] })] }) },
    children,
  }],
});

const out = path.join(__dirname, 'EngTrace-inference-pricing.docx');
Packer.toBuffer(doc).then((buf) => { fs.writeFileSync(out, buf); console.log('wrote', out, buf.length, 'bytes'); });
