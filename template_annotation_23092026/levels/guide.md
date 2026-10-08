# EngTrace: rating the difficulty of templates

## What this is

You will see the 30 templates of your branch, one at a time. Each is shown as one problem drawn from it, with
its reference solution. For each template you make two choices:

1. **Its level**: Easy, Intermediate or Advanced, judged on the three dimensions below.
2. **Whether the problem states the governing formula or names the method**: yes or no.

The other two experts of your branch rate the same 30 templates independently. There is no answer to match:
we measure how far independent ratings agree. Please work on your own and do not discuss the items with
colleagues until you have sent your file. The whole task takes 30 to 45 minutes.

## How to run

You need Python 3.11 or newer. In the folder that holds `app.py`:

    pip install streamlit
    streamlit run app.py

Your browser opens the app. Pick your id in the sidebar. The app saves after every item; to pause, close the
tab and stop the app, and run `streamlit run app.py` again to continue where you left off. When the app says
you are done, send the file `<your id>.jsonl`, which it has written next to `app.py`, to the coordinator.

## The rubric

Judge the template, not the particular numbers of the problem shown: another draw from the same template
changes the values, not the problem.

| dimension | toward Easy | toward Advanced |
|---|---|---|
| Conceptual complexity | one isolated principle | the synthesis of several principles or domains |
| Mathematical sophistication | direct algebraic substitution | differential equations, or an iterative or implicit solution |
| Procedural depth | a short chain of steps | a long chain of interdependent steps, each needing the ones before |

- **Easy**: low on all three dimensions: one principle, applied by direct substitution, in a short chain.
- **Intermediate**: in between.
- **Advanced**: a synthesis of principles or domains, a differential equation or an iterative or implicit
  solution, or a long chain of interdependent steps.

The level is your overall judgement over the three dimensions, as you would set it for a homework or
examination problem in your field.

**The formula question.** Answer **yes** if the problem's text gives the equation to use (for example "using
P = ρgh") or names the method, model or correlation to apply (for example "using the Rackett equation", "by
the phasor method"). Answer **no** if the solver has to choose the principle or the method. Data the problem
gives (constants, property values, a table) are not a formula, and neither is a definition of notation (such
as "sinc(x) = sin(πx)/(πx)").

## Please

- Rate each template on its own; you can revisit an item from the sidebar, and submitting again replaces your
  earlier answer.
- Use the note for anything that made the choice hard, for example a problem whose level depends on what the
  solver is expected to know.
- Questions about the task go to the coordinator, not to the other reviewers. Please do not share the folder:
  the problems are the benchmark's private test set until the paper is published.
