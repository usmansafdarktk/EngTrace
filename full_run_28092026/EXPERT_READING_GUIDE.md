# EngTrace reading request, October 2026: the guide

## What this is

The benchmark has been run. Eleven language models answered 2,250 engineering problems, one answer each, and every
answer was scored by automatic checks that were validated against the labels you gave in September. Before the
results go into the paper, four things need an engineer's reading rather than a check. Your items are from your own
branch only, and nothing has to be solved from scratch: every item shows the problem and its reference solution.

Each task is its own folder with its own app, and the folder's `guide.md` describes that task only:

- `B1_final_answers`: is a model's final answer correct? About 20 items, a minute or two each.
- `B3_milestones`: does a model's working obtain a given quantity? About 13 items, two or three minutes each.
- `B2_wrong_answers`: why did a wrong answer go wrong? 20 to 41 items, three to five minutes each.
- `B4_templates` (chemical and civil engineers only): does a problem's wording pin a single answer? Two to six
  items, five to ten minutes each.

If you cannot do all of them, the order that helps most is B4, then B1, then B3, then B2.

## How to run

Each task folder holds `app.py`, `README.txt`, `guide.md` and a `tasks` folder. As for the paraphrase check:

1. You need Python 3.11 or newer. In the task's folder: `pip install streamlit` (once on your machine), then
   `streamlit run app.py`.
2. Your browser opens the app. Pick your id in the sidebar.
3. Work through the items. The app saves after every item; to pause, close the tab and stop the app, and run it
   again to continue where you left off. You can jump back to any item and submit it again.
4. When the app says you are done, send the file `<your id>.jsonl`, which the app writes next to that folder's
   `app.py`, to the coordinator. Each task folder produces its own file.

Everything is shown as plain text, exactly as the model wrote it or read it, never rendered.

## B1: a final answer

You see the problem, the reference answer, and the final answer a model stated. The model's full working and the
full reference solution are one click away. Say whether the model's final answer is **correct**, **partially
correct**, **incorrect**, or **no answer stated**.

- Judge the model's final answer against the reference answer, part by part.
- A correct rounding at the precision the model shows is correct. The same quantity in another unit, as a fraction,
  in scientific notation, or with π left symbolic is correct.
- Where the question prescribes a rounding ("to four decimals", "to the nearest whole unit, round half up"), hold the
  answer to it: a last digit that does not match the prescribed rounding is incorrect.
- Partially correct: an answer with several parts, some right and some wrong or missing.
- No answer stated: the model gives no final value for what was asked.
- If you think the reference itself is wrong, or that the question allows the model's answer as well as the
  reference's, still give the verdict the reference implies, and say so in the note. The note is what we act on.

## B3: a quantity in the working

You see the problem, one quantity from the reference solution (the name the solution uses for it, and its value),
and a model's full working. Say whether the model's working obtains that quantity:

- **yes, the working obtains it**: the model states that value, allowing another unit, a different rounding, or
  another form, including an expression it then evaluates;
- **no, it never obtains it**;
- **no, and its route does not need it**: the model solves the problem by a valid route that never passes through
  this quantity. Say which route in the note.

Use the reference solution to see what the quantity is; it appears there under the same name or value. An
intermediate the model computes but writes wrongly counts as not obtained.

## B2: why a wrong answer went wrong

You see the problem and a model's full working whose final answer the check scored as wrong; the reference solution
is one click away. Reading the working from the top, choose the **first** of these six questions you would answer
yes to. An error at an earlier stage invalidates everything after it, so a model that chose the wrong governing
equation is not labelled for the arithmetic that followed.

1. **Hallucination**: does the model use a value, constant or equation that does not exist in any engineering
   reference?
2. **Setup / Assumption Error**: is the problem framed wrongly before any equation is applied (a wrong boundary
   condition, an ignored constraint, a misidentified system configuration)?
3. **Formula / Principle Error**: is the governing equation or principle wrong, given that the setup is right?
4. **Unit / Dimensional Error**: does dimensional analysis fail (a missing conversion, mixed unit systems, a value
   substituted in the wrong unit)?
5. **Sign / Direction Error**: is the sign or physical direction of a quantity wrong, given that the formula and
   units are right?
6. **Calculation Error**: does the arithmetic or algebra produce the wrong result, given that every earlier step is
   right?

Two more options cover what the six do not. **No error: the answer is correct, or the question admits it** is for a
working you find sound whose answer differs from the reference because the question leaves a method, a data source
or a convention open, or because the reference is wrong; say which in the note. **Incomplete: the working stops
before an answer** is for a trace cut off before it concludes.

Then paste, into the excerpt box, the shortest piece of the working that shows the error: the line with the wrong
equation, the wrong substitution or the wrong arithmetic. The excerpt is required for the six error categories.

## B4: questions about a template

You see a problem as the benchmark posed it, its reference solution, and one or three model answers to it, with a
short note on what the run found for that template. Each item asks two to four specific questions, mostly about
whether the question's wording determines a single correct answer at the precision the check uses, 0.2%, and whether
the model answers shown are acceptable. There is no expected answer; your reading is the evidence. The note box asks
what wording, if any, would remove an ambiguity you see.

## Please

Work alone, and do not discuss items with the other reviewers until everyone has submitted. Answer every item if you
can; an item you cannot judge can be submitted with a note saying so. Questions about the task go to the
coordinator, not to the other reviewers. Please do not share the folders: the problems are the benchmark's private
test set until the paper is published.
