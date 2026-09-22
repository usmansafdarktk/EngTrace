// Build EngTrace-inference-pricing.docx: model prices and what 2,250 items would cost.
//
//   npm install docx            (once, anywhere on NODE_PATH)
//   node build_pricing_doc.js   -> EngTrace-inference-pricing.docx beside this script
//
// Prices are per million tokens, read on OpenRouter 19-21 September 2026. The two cost
// columns use EngTrace's own token counts, measured on the evaluator pilot's traces
// (evaluator_pilot_17092026/analysis/inference_cost.py):
//   direct answer        318 input / 813 output tokens per item (pilot's non-reasoning models)
//   reasons like GPT-5   270 input / 5,232 output tokens per item (GPT-5, reasoning included)
// for the full benchmark of 150 templates x 15 instances = 2,250 items, no retries.

const fs = require('fs');
const path = require('path');
const {
  Document, Packer, Paragraph, TextRun, Table, TableRow, TableCell, WidthType,
  ShadingType, AlignmentType, HeadingLevel, BorderStyle, LevelFormat, Footer, PageNumber,
} = require('docx');

const ITEMS = 2250;
const DIRECT = { in: 318, out: 813 };
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
const W = [3060, 1200, 1200, 1950, 1950];          // sums to 9360 = US Letter minus 1" margins
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
    'Model', 'Input $/M', 'Output $/M', 'Direct answer', 'Reasons like GPT-5',
  ].map((h, i) => cell(h, i, { bold: true, fill: HEAD_FILL, size: 18, color: ACCENT, keepNext: true })) });
  const body = rows.map((r, k) => {
    const o = { fill: k % 2 ? ZEBRA : undefined, keepNext: k < rows.length - 1 };
    return new TableRow({ cantSplit: true, children: [
      cell(r[0], 0, o), cell(price(r[1], dp), 1, o), cell(price(r[2], dp), 2, o),
      cell(money(cost(r, DIRECT)), 3, o), cell(money(cost(r, REASON)), 4, { ...o, bold: true }),
    ] });
  });
  return new Table({ width: { size: 9360, type: WidthType.DXA }, columnWidths: W, rows: [head, ...body] });
}

const bullet = (children) => new Paragraph({ numbering: { reference: 'dots', level: 0 }, children,
  spacing: { after: 90 } });

const children = [
  new Paragraph({ heading: HeadingLevel.TITLE, spacing: { after: 60 },
    children: [new TextRun({ text: 'EngTrace inference pricing', font: FONT, size: 44, bold: true, color: INK })] }),
  para([run('Model prices per million tokens and the estimated cost of one full benchmark run', { color: MUTED, size: 22 })], { after: 60 }),
  para([run('Prices read on OpenRouter, 19–21 September 2026', { color: MUTED, size: 18, italics: true })], { after: 280 }),

  para([run('The tables list each model’s input and output price. The last two columns estimate what '),
    run('one full EngTrace run', { bold: true }),
    run(' costs for that model: 150 templates × 15 instances = 2,250 items, inference only, no evaluation. '
      + 'They use token counts measured on EngTrace’s own evaluator pilot, because a model’s cost '
      + 'depends on how much it writes, and that is not known in advance for a new model:')]),
  bullet([run('Direct answer: ', { bold: true }), run('about 318 input and 813 output tokens per item, the average of the pilot’s non-reasoning models.')]),
  bullet([run('Reasons like GPT-5: ', { bold: true }), run('about 270 input and 5,232 output tokens per item, GPT-5’s measured average including its reasoning.')]),
  para([run('A model’s real cost will usually fall between the two. The heaviest reasoner in the pilot, '
    + 'DeepSeek-R1, wrote about 10,000 output tokens per item, roughly twice the second column.')], { before: 60 }),

  new Paragraph({ heading: HeadingLevel.HEADING_2, spacing: { before: 240, after: 100 },
    children: [new TextRun({ text: 'Two things to watch', font: FONT, size: 26, bold: true, color: ACCENT })] }),
  bullet([run('Output is where the money goes for reasoning models. ', { bold: true }),
    run('Output prices are four to eight times input, and engineering problems produce long reasoning, so the cost tracks how much each model writes.')]),
  bullet([run('OpenRouter prices move between upstreams from day to day. ', { bold: true }),
    run('These were read around 19 to 21 September 2026 and should be re-read before the run.')]),
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
  bullet([run('Cost = 2,250 × (input tokens × input price + output tokens × output price) ÷ 1,000,000.')]),
  bullet([run('Token counts come from the pilot’s answers to 60 items covering all five branches and all three difficulty levels, under the benchmark’s own prompt.')]),
  bullet([run('Re-runs after truncated or failed answers are not included; the pilot needed 0–15% extra for most models and up to 45% for a model whose reasoning ran away.')]),
  bullet([run('Half-price batch endpoints, where a provider offers them, would roughly halve the frontier rows.')]),
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
