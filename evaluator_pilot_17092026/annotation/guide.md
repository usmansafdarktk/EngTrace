# EngTrace annotation guide

## What you are doing, and why it matters

You will read solutions to engineering problems written by AI systems, and judge them
**step by step**. Each solution states a problem, works through it, and gives a final
answer. Your job is to say, for every step, whether the reasoning in it is right.

The reason this is worth your time: there are several automated systems that claim to
judge whether AI reasoning is sound, and at present there is no way to tell which of
them is correct, because each one has only ever been compared against the others. Your
labels become the reference they are measured against. A final answer can be right for
the wrong reasons, and wrong after almost entirely correct work, so only a step-level
judgement settles it.

You will annotate **68 solutions**: the 60 in your own discipline, and 8 shared ones
from other disciplines that every annotator labels so we can check we are all applying
the same standard. Expect roughly **8 to 12 minutes each**, so about 10 hours in total.
Work in as many sittings as you like; progress is saved after every solution.

## What you will not be told, on purpose

You will not see which AI system wrote a solution, and you will not see what any
automated judge thought of it. Both would influence your judgement. Each solution
appears under a code such as `T-4f2a9c`.

## Start with the calibration set

The first 10 solutions in your list are the calibration set, marked as such. Everyone
labels these. When they are done, we compare and discuss any systematic disagreement
before the bulk of the work. If something in this guide turns out to be ambiguous, that
is when we fix it.

## The screen

**Problem.** The question, exactly as the AI system received it.

**Reference solution.** Collapsed by default. It is one worked solution, not the only
valid one. Open it whenever it helps, but do not treat departures from it as errors:
there is a separate label for a different valid route.

**Steps.** The solution split into its steps, each with its own label.

**Milestones.** Quantities a correct derivation reaches. For each, say whether this
solution obtains it.

**The solution as a whole.** Final answer, overall soundness, and your confidence.

## The step labels

**Correct.** What the step states and computes holds. The method is sound and the
numbers are right.

**Correct, different route.** The step is right, but it is not how the reference
solution proceeds: a different valid correlation, a different but legitimate
assumption, a different order of work. Use this rather than marking it wrong. It is
important that we can distinguish "wrong" from "not what we expected".

**Incorrect.** Something in the step is wrong. You then choose the kind:

- *Calculation error* — the method is right, the arithmetic is not. The step sets up
  the correct expression and then evaluates it wrongly.
- *Conceptual error* — the wrong equation, principle, or assumption. Using the wrong
  correlation, ignoring a term that matters, mixing incompatible units.
- *Unsupported* — a value or claim asserted with no basis: a coefficient that appears
  from nowhere, a table value that was never looked up, a conclusion that does not
  follow.

**No claim to check.** The step restates the question, announces a plan, or lists what
is given, without asserting anything that can be right or wrong. Do not label such a
step correct just because nothing in it is wrong.

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
  done, is an error — conceptual if the method mishandles units, calculation if it is
  a slipped factor.
- **A step that corrects itself.** If a step states a wrong value and then visibly
  fixes it within the same step, and ends with the right value, label it correct and
  add a note. These are rare, and the note lets us find them.

## Milestones

Each problem has a handful of quantities that a correct derivation passes through,
shown with their values. For each one:

- **Obtained** — the solution arrives at this quantity. It counts whichever form, unit
  or rounding it appears in, and whether or not it is labelled the same way. What
  matters is that the solution computed it.
- **Not obtained** — the solution never arrives at this quantity.
- **Not needed on this route** — a valid alternative method that legitimately never
  needs this quantity. Only for methods that are genuinely valid, not for solutions
  that simply skipped a step.

## The solution as a whole

**Final answer.** Correct, incorrect, partly correct (when a problem asks for several
quantities and some are right), or no final answer stated.

**Is the reasoning sound overall?** Your holistic judgement. A solution can reach the
right answer through flawed reasoning, and that is a "no".

**Confidence.** Low confidence is useful information, not an admission of anything. Use
it whenever the problem sits at the edge of your expertise or the solution is hard to
follow.

**Comments and the discussion flag.** Use the comment box for anything a reader would
need to know. Tick "flag for discussion" when a solution raises a question the guide
does not answer — an ambiguous problem statement, a reference solution you believe is
wrong, a step you genuinely cannot resolve. Flagged solutions are reviewed together.

## Things to avoid

- **Do not reward or penalise style.** Terse and verbose solutions are judged the same
  way. A step that is brief but right is correct.
- **Do not treat the reference solution as the only correct method.**
- **Do not let the final answer decide the steps.** A right answer does not make the
  work right, and a wrong answer does not make every step wrong. This is the single
  most common way step-level annotation goes wrong.
- **Do not guess at which system wrote a solution**, and do not let a hunch about it
  affect the labels.

## Practical notes

- Progress saves when you submit each solution. You can close the app and return.
- You can reopen a solution you have already submitted and change it; the latest
  version is the one we use.
- Solutions are in a different order for each annotator, so do not expect your list to
  match a colleague's.
- If two of you discuss a specific solution while labelling, say so in the comment,
  since our agreement statistics assume independent judgements.
