# E5 — deterministic first, a judge only where E3 falls short

Run 2026-09-19. **300 of 300 traces, 142 sent to the judge, 0 judge failures, $0.47**
(plus $0.43 validating the judge first). Evaluator `evaluators/e5_hybrid.py`; every
figure is printed by `analysis/e5_analysis.py report`.

## How it works

1. **E3 runs first**, unchanged: which of the item's milestones does the trace state?
2. **For each milestone E3 did not find**, one call to **MiMo-V2.5-Pro** (D-088) per
   trace asks: `REACHED` (the trace obtains it, perhaps in another form, unit or
   rounding), `NOT_NEEDED` (a valid route that does not pass through it), or `MISSING`.
   The judge sees the question, the trace and each milestone's name and value — not
   the gold solution.
3. **E5-strict** credits E3's matches plus `REACHED`. E5-lenient also credits
   `NOT_NEEDED`, and is **not a valid score** — see below.

## The judge was validated before any score was trusted

On answers known without asking it (D-087): milestones E3 *did* find, with their true
values, and the same milestones with every value ×1.37 — quantities the trace never
states. 40 traces, 88 milestones each way.

| | REACHED | MISSING | NOT_NEEDED |
|---|---|---|---|
| **True values** — should be REACHED | **67 (76%)** | 20 | 1 |
| **Values ×1.37** — should be MISSING | **0** | 66 | **22 (25%)** |

- **REACHED is trustworthy.** Never given to a fake value — 0 of 88, a 95% upper bound
  of about 3.4% — and given to 76% of true ones. E5-strict is therefore valid and, if
  anything, **conservative**: it under-credits by roughly a quarter of the milestones it
  sends to the judge.
- **NOT_NEEDED is not.** MiMo used it to excuse a quarter of values the trace never
  computed. A lenient score would credit made-up numbers. Reported below only as a
  diagnostic of how much the judge would excuse.

## Results

| Model | E3 | E4 | **E5 strict** | *E5 lenient (invalid)* | E0 reasoning F1 |
|---|---|---|---|---|---|
| gpt-5 | 0.867 | 0.863 | **0.962** | *0.981* | 0.474 |
| claude-opus-4.7 | 0.877 | 0.877 | **0.946** | *0.983* | 0.449 |
| gemini-3.1-pro | 0.860 | 0.860 | **0.929** | *0.939* | 0.454 |
| deepseek-r1 | 0.855 | 0.855 | **0.923** | *0.961* | 0.402 |
| llama-3.1-70b | 0.285 | 0.281 | **0.321** | *0.341* | 0.115 |

**The judged fraction — the statistic the Suggested Actions asked for: 23.6%.** Of 1,245
milestones across the 300 traces, 951 (76.4%) were settled deterministically by E3 and
294 needed the judge. Of those 294 the judge found 73 reached (25%), excused 36, and
ruled 185 missing (63%).

Two things the table shows:

1. **E5 recovers what E3's number-matching cannot see.** 73 milestones the trace did
   obtain — in a form, unit or rounding E3 misses — are credited, and the validation
   says those credits are sound.
2. **It compresses the frontier.** The four frontier models land between 0.92 and 0.96,
   close to the ceiling, where E3 had them between 0.855 and 0.877. That is what a
   milestone-coverage score should do when the models mostly do reach the right
   quantities; it also means E5 separates the frontier models less than E0 does. The
   gap to Llama (0.321) is untouched.

## A caveat that applies to every evaluator so far

Correlation with final-answer correctness:

| | corr |
|---|---|
| E0 reasoning F1 | 0.741 |
| E3 | 0.471 |
| E5 lenient | 0.431 |
| E5 strict | 0.384 |

**This does not rank the evaluators, and must not be read as E0 being best.** "Final-answer
correctness" here is E0's own check, which reads the wrong number on `rackett` (E0-F1)
and has no real answer on 8 of 15 templates (E0-F2). And E0's reasoning score is
**mechanically coupled** to it: E0 decides which traces reach its judges, and so which
steps are recovered, *from* that check. E0 agrees with its own yardstick by construction.

E5 plausibly correlates less precisely because it credits correct intermediate reasoning
on traces that the flawed check marks wrong. The only fair referee is the expert
annotation (X1). Until it exists, these correlations are a sanity check — each evaluator
separates Llama from the rest — and nothing more.

## Cost and run

$0.47 for 142 judge calls, fetched concurrently (32 at a time) before scoring and
replayed from the store, on the laptop's system Python: E5 needs no scorer models. With
validation, E5 cost about $0.90.
