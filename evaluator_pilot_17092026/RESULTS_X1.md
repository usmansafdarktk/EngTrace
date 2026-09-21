# X1 — the evaluators against the expert labels

Scored 2026-09-22. **300 traces, each labelled by three experts from the trace's own
branch** (15 experts, 900 trace-annotations for the ground truth, plus 120 shared
calibration annotations). The labels are in `experts_filled_labels/`. The headline
comparison is `annotation/score_against_labels.py --labels experts_filled_labels`; the
intervals, baselines and step-level analysis are `analysis/x1_analysis.py`.

This is the referee the pilot was built for. Until now every comparison was evaluators
against each other, or against E0's own final-answer check.

## The ground truth

| | steps | milestones | final answer |
|---|---|---|---|
| Fleiss κ, pooled over branches | **0.70** | 0.81 | 0.87 |
| chemical / civil / electrical / industrial / mechanical | 0.31 / 0.65 / 0.58 / 0.87 / 0.81 | 0.39 / 0.99 / 0.78 / 0.87 / 0.66 | 0.42 / 0.97 / 0.96 / 0.98 / 0.90 |
| calibration set, all 15 experts across branches | 0.72 | — | 0.76 |

Majority vote leaves no step, milestone or verdict disputed. Chemical engineering
agrees least: one of its three experts marked 98% of steps correct, against 77–85%
for the other two. Every result below is also reported **without chemical
engineering**, and no conclusion changes. The app recorded widely different entry
times by branch (a median of 0.3 s to 128 s per solution), so the app's timing does not
document how long the judgement took.

The experts judged **236 traces sound and 64 unsound**, and found 279 incorrect steps
among 2,091.

## Finding 1 — at the trace level, the final answer decides almost everything

**The experts' own final-answer verdict predicts their "reasoning sound" verdict with
AUROC 0.958**, better than every evaluator. Of the 228 traces with a correct final
answer, the experts judged only **3** unsound.

| Separating sound from unsound traces | AUROC | 95% CI | minus E0 |
|---|---|---|---|
| E0 — the published framework | 0.851 | 0.805–0.893 | — |
| E0-3J | 0.792 | 0.737–0.844 | **−0.059** (significant) |
| E1 | 0.849 | 0.802–0.891 | −0.002 |
| E2, 72B, fraction of steps ok | 0.839 | 0.782–0.889 | −0.011 |
| E2, 72B, lowest step reward | 0.872 | 0.825–0.913 | +0.022 |
| E3 | 0.819 | 0.754–0.876 | −0.032 |
| E4 | 0.819 | 0.754–0.876 | −0.032 |
| E5 | 0.869 | 0.812–0.921 | +0.018 |
| *baseline: E0's own answer check* | 0.798 | 0.751–0.841 | **−0.053** (significant) |
| *baseline: the experts' answer verdict* | **0.958** | 0.925–0.984 | **+0.107** (significant) |

Three consequences:

- **No evaluator beats E0 at the trace level.** Every interval for the difference
  includes zero except E0-3J, which is significantly *worse*. Without chemical
  engineering, E3 and E4 are also significantly worse (−0.09).
- **The trace-level test cannot rank reasoning evaluators on this slice.** "Right
  answer, flawed reasoning" is the case a reasoning evaluator exists for, and the
  slice has 3 such traces. Ranking the evaluators on it needs a set built to contain
  that case: harder items, weaker models.
- **The single largest available improvement is a correct final-answer check.** E0's
  check agrees with the experts on **only 76%** of traces, which puts a number on
  E0-F1 (a restated input read as the answer) and E0-F2 (no number on the Answer
  line). As a predictor of soundness it scores 0.798, where an accurate answer check
  scores 0.958.

## Finding 2 — milestone coverage is accurate, and E5's judge earns its place

| Milestones, against the experts' "obtained" | precision | recall | F1 |
|---|---|---|---|
| E3 — deterministic matching | 0.923 | 0.915 | 0.919 |
| E4 — E3 plus stated-arithmetic checks | 0.923 | 0.915 | 0.919 |
| **E5 — E3, then a judge on what E3 misses** | **0.928** | **0.990** | **0.958** |

E5's judge recovers nearly every milestone E3's number-matching misses (recall 0.915
→ 0.990) at no cost in precision. That confirms, against experts, what the synthetic
validation in RESULTS_E5 suggested: its REACHED credits are sound. E4 adds nothing
over E3 at this level, as RESULTS_E4 found.

## Finding 3 — the 72B process reward model finds wrong steps, but over-flags badly where it matters most

| Steps, against the experts' "incorrect" | AUROC | 95% CI | precision | recall | F1 |
|---|---|---|---|---|---|
| Qwen2.5-Math-PRM-72B | **0.860** | 0.833–0.887 | 0.469 | 0.629 | 0.538 |
| Qwen2.5-Math-PRM-7B | 0.798 | 0.764–0.830 | 0.387 | 0.513 | 0.441 |
| VersaPRM | 0.636 | 0.594–0.678 | 0.354 | 0.082 | 0.134 |
| *inside correct-answer traces only (81 incorrect steps):* | | | | | |
| Qwen2.5-Math-PRM-72B | 0.759 | 0.707–0.805 | **0.117** | 0.290 | 0.167 |

- **The 72B ranks steps well overall** (AUROC 0.86), but flags a third more steps
  than the experts do (358 against 267), and under half its flags are real errors.
- **Inside correct-answer traces, where a step-level evaluator would add most, 88% of
  its flags are false alarms** and it finds 29% of the real errors. This is the
  over-flagging RESULTS_E2 could only suspect.
- **VersaPRM's failure is confirmed:** it finds 8% of the incorrect steps.

## Finding 4 — the experts contradict E2's reordering of the frontier

| | experts: traces sound | E0 | E1 | E5 | E2 (72B) |
|---|---|---|---|---|---|
| gpt-5 | **0.983** | 0.474 | 0.470 | 0.962 | 0.821 |
| claude-opus-4.7 | 0.917 | 0.449 | 0.448 | 0.946 | 0.864 |
| gemini-3.1-pro | 0.933 | 0.454 | 0.452 | 0.929 | 0.876 |
| deepseek-r1 | 0.933 | 0.402 | 0.401 | 0.923 | 0.921 |
| llama-3.1-70b | 0.167 | 0.115 | 0.117 | 0.321 | 0.560 |

The experts put **GPT-5 first**, and so do E0, E1 and E5. E2 put it **fourth** —
RESULTS_E2 flagged that as possibly a penalty on terse steps, and the experts say it
was. Every evaluator separates Llama from the frontier. Among the four frontier
models the experts' figures sit within about 7 points of each other, too close for
this slice to rank the evaluators on them.

## What this means

1. **Fix the final-answer check first.** It is the dominant trace-level signal, and
   E0's parser is wrong on a quarter of traces.
2. **For reasoning progress, milestone coverage with a residual judge (E5) is the
   strongest candidate:** F1 0.958 against experts, and cheap ($0.47 for 300 traces).
3. **Step-level error detection is not solved by an off-the-shelf process reward
   model.** The 72B is the best available but over-flags badly on correct-answer
   traces. Adapting it to engineering, or recalibrating its threshold against these
   labels, is the obvious next experiment.
4. **The next labelled set should be built for the hard case**: correct answers with
   flawed reasoning. This slice has 3 such traces, too few to rank reasoning
   evaluators where it matters most.

## Follow-ups this turned up

- The judge probe labelled Claude's `aoq_ati_rectifying#3` step as CLEAN. The experts
  unanimously found a real slip in it: 0.93^49 = 0.0285538, not 0.02857, so P(X=1) is
  0.0999, not 0.1000. It should be relabelled SLIP in `judge_probe.RELABEL`, and E2's
  probe validation re-run.
- The guide does not say how to treat a notational slip whose computed result is
  right (GPT-5's `normal_depth_iteration#1`: a dropped minus sign in a displayed
  expression, correct value). The experts called it correct; the guide should say so
  explicitly.
