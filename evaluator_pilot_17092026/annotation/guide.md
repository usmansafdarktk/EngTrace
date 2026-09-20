# EngTrace annotation guide

## What you are doing, and why it matters

You will read solutions to engineering problems written by large language models, and
judge them **step by step**. Each solution states a problem, works through it, and
gives a final answer. Your job is to say, for every step, whether the reasoning in it
is right.

The reason this is worth your time: several automated systems claim to judge whether
model reasoning is sound, and at present there is no way to tell which of them is
correct, because each one has only ever been compared against the others. Your labels
become the reference they are measured against. A final answer can be right for the
wrong reasons, and wrong after almost entirely correct work, so only a step-level
judgement settles it.

You will annotate **68 solutions**: the 60 in your own discipline, and 8 shared ones
from other disciplines that every annotator labels so we can check we are all applying
the same standard. Expect roughly **8 to 12 minutes each**, so about 10 hours in total.
Work in as many sittings as you like.

## What you will not be told, on purpose

You will not see which model wrote a solution, and you will not see what any automated
judge thought of it. Both would influence your judgement. Each solution appears under a
code such as `T-4f2a9c`.

## Start with the calibration set

Ten of your solutions are marked as the calibration set, and they come first. Every
annotator labels these same ten. When they are done, we compare and discuss any
systematic disagreement before the bulk of the work. If something in this guide turns
out to be ambiguous, that is when we fix it.

## What you are given for each solution

**The problem**, exactly as the model received it.

**The reference solution.** One worked solution, not the only valid one. Consult it
whenever it helps, but do not treat departures from it as errors: there is a separate
label for a different valid route.

**The steps**, as the solution was split. Each step gets a label.

**The milestones.** Quantities a correct derivation reaches, with their values. For
each, say whether this solution obtains it.

**The solution as a whole.** Final answer, overall soundness, and your confidence.

## The step labels

**correct** — What the step states and computes holds. The method is sound and the
numbers are right.

**alternative_correct** — The step is right, but it is not how the reference solution
proceeds: a different valid correlation, a different but legitimate assumption, a
different order of work. Use this rather than marking it wrong. It matters that we can
distinguish "wrong" from "not what we expected".

**incorrect** — Something in the step is wrong. You then choose the kind:

- **calculation** — the method is right, the arithmetic is not. The step sets up the
  correct expression and then evaluates it wrongly.
- **conceptual** — the wrong equation, principle, or assumption. Using the wrong
  correlation, ignoring a term that matters, mixing incompatible units.
- **unsupported** — a value or claim asserted with no basis: a coefficient that appears
  from nowhere, a table value that was never looked up, a conclusion that does not
  follow.

**not_a_claim** — The step restates the question, announces a plan, or lists what is
given, without asserting anything that can be right or wrong. Do not label such a step
correct just because nothing in it is wrong.

## The rule that matters most: where an error starts

When a step makes an error and later steps carry the bad number forward, **label only
the step where the error originates as incorrect**. A later step that uses the wrong
value correctly is *correct*: its own reasoning and arithmetic are sound, given what it
was handed.

The reason is that we need to identify *where* reasoning breaks down. If every step
after an error were marked wrong, every judgement would collapse into "the answer was
wrong", which we already know.

The exception: a later step that introduces a *new* error of its own is also incorrect.

## Judging numbers

- **Rounding is not an error.** Carrying 3 significant figures instead of 5, or
  rounding at a different point, is fine as long as the result is what the stated
  method gives. Judge as you would a student's working: is the number the one this
  method produces, to the precision shown?
- **A wrong digit is a calculation error**, even if the final answer survives it.
- **Units matter.** A value correct in the wrong unit, or a conversion that was never
  done, is an error — conceptual if the method mishandles units, calculation if it is a
  slipped factor.
- **A step that corrects itself.** If a step states a wrong value, visibly fixes it
  within the same step, and ends with the right value, label it correct and add a note.
  These are rare, and the note lets us find them.

## Milestones

For each milestone, one of:

- **reached** — the solution arrives at this quantity. It counts whichever form, unit or
  rounding it appears in, and whether or not it is labelled the same way. What matters
  is that the solution computed it.
- **not_reached** — the solution never arrives at this quantity.
- **not_needed** — a valid alternative method that legitimately never needs this
  quantity. Only for methods that are genuinely valid, not for solutions that simply
  skipped a step.

## The solution as a whole

**final_answer** — one of `correct`, `incorrect`, `partial` (a problem asking for
several quantities where some are right), `not_stated`.

**reasoning_sound** — `yes` or `no`, your holistic judgement. A solution can reach the
right answer through flawed reasoning, and that is a `no`.

**confidence** — `high`, `medium` or `low`. Low confidence is useful information, not an
admission of anything. Use it whenever the problem sits at the edge of your expertise or
the solution is hard to follow.

**comment** and **flagged** — use the comment for anything a reader would need to know.
Flag a solution when it raises a question this guide does not answer: an ambiguous
problem statement, a reference solution you believe is wrong, a step you genuinely
cannot resolve. Flagged solutions are reviewed together.

## Things to avoid

- **Do not reward or penalise style.** Terse and verbose solutions are judged the same
  way. A step that is brief but right is correct.
- **Do not treat the reference solution as the only correct method.**
- **Do not let the final answer decide the steps.** A right answer does not make the
  work right, and a wrong answer does not make every step wrong. This is the single most
  common way step-level annotation goes wrong.
- **Do not guess at which model wrote a solution**, and do not let a hunch about it
  affect the labels.
- **Judge independently.** If you discuss a specific solution with another annotator
  while labelling, say so in the comment, since our agreement statistics assume
  independent judgements.

## Recording your answers

Two ways, and they produce the same thing. Use whichever you prefer.

**A. The annotation app.** It shows one solution at a time with the labels as buttons,
saves as you go, and will not let you submit a solution with anything unanswered. Start
it with the command you were sent; work is resumable, and a submitted solution can be
reopened and changed.

**B. Your workbook file.** You were sent `<your-id>.json`, which contains all 68
solutions with empty fields for you to fill in, in any text or JSON editor. Fill in the
quoted values and leave everything else untouched:

```json
{
  "code": "T-4f2a9c",
  "steps": [
    {"index": 0, "text": "...", "label": "correct", "error_type": "", "note": ""},
    {"index": 1, "text": "...", "label": "incorrect", "error_type": "calculation",
     "note": "16*2.9 evaluated as 44.4"}
  ],
  "milestones": [{"id": "Tr", "value": 0.6042, "status": "reached"}],
  "final_answer": "correct",
  "reasoning_sound": "yes",
  "confidence": "high",
  "comment": "",
  "flagged": false
}
```

The permitted values are exactly those named above: `label` is one of `correct`,
`alternative_correct`, `incorrect`, `not_a_claim`; `error_type` is `calculation`,
`conceptual` or `unsupported` and is required when, and only when, the label is
`incorrect`; `status` is `reached`, `not_reached` or `not_needed`; `final_answer` is
`correct`, `incorrect`, `partial` or `not_stated`; `reasoning_sound` is `yes` or `no`;
`confidence` is `high`, `medium` or `low`.

Send the file back when you are done, or in parts as you go. It is checked on arrival
and you will be told about anything missing or misspelled. You can switch between the
app and the file: your workbook can be produced from what you have already submitted,
and a filled workbook can be loaded back into the app.
