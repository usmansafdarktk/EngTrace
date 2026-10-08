# EngTrace: grading symbolic final answers

## What this is

Some problems in the benchmark ask for an expression rather than a single number: a signal y(t), an impulse
response h[n], a velocity field v(x, y), a bit error rate written with the Q-function. The automatic scorer reads
the numbers in a model's final answer and compares them with the reference, so an answer written in another but
equal form (y²/2 for 0.5y², a phase shifted by a full turn, exact roots where the reference rounds them) can be
scored wrong. You will see model final answers to such problems, each beside the reference answer, and say whether
the two are equivalent. Some of the answers were scored right and some wrong; you are not told which, nor which
model wrote them. Your grades decide where an automatic equivalence check may be used, so please grade what is on
the page, not what the model may have meant.

The items come in groups: all the answers to one kind of problem together. The whole task takes about an hour.

## How to run

You need Python 3.11 or newer. In the folder that holds `app.py`:

    pip install streamlit
    streamlit run app.py

Your browser opens the app. Pick your id in the sidebar. The app saves after every item; to pause, close the
tab and stop the app, and run `streamlit run app.py` again to continue where you left off. When the app says
you are done, send the file `<your id>.jsonl`, which it has written next to `app.py`, to the coordinator.

## The grades

Each item shows the problem, the reference answer, and the model's final answer as the scorer reads it (from the
response's final-answer heading, up to 700 characters). The full response and the full reference solution are one
click away; use them to see what the final answer says if the box cuts it off, not to credit working that the
final answer does not state.

**Equivalent**: the final answer states the same quantity as the reference, to the precision the reference
shows. For example:

- algebraically equal: `-6xy - y²/2` for `-6xy - 0.5y²`;
- equal through a standard identity: `Q(x) = ½·erfc(x/√2)`; `sin θ = cos(θ − 90°)`; a phase that differs by a
  whole number of turns (2π rad, 360°);
- exact where the reference rounds: `(1 + √2)ⁿ` where the reference writes `2.41ⁿ`; a value with more digits
  that rounds to the reference's;
- in another unit, where the unit is stated: a displacement in mm where the reference is in m;
- a quantity the problem gives left as a symbol (`Eb/N0` where the problem states its value), if substituting
  the problem's value gives the reference;
- where the reference is one value written as an expression (`BER ≈ 0.750·Q(8.78)`), the same value, written as
  an expression or as a number. For these items the app also prints the reference's value. Because the Q-function
  falls steeply, a value computed from a slightly different rounding of its argument can differ in its second
  digit; judge whether the difference comes from rounding.

Every part the problem asks for must be there and equivalent (for example both G(f) and Ψ(f), or both x[n] and
ω).

**Not equivalent**: a different quantity, a missing or extra term or factor, a missing part, a different
approximation that gives a different value, or numbers that differ from the reference by more than its rounding.

**Unreadable**: no single final answer can be read: the response stops before one, states two different final
answers, or the box and the full response hold none.

The note is optional. Please use it when the reference itself looks wrong to you, or when an answer uses a
different standard approximation that you would accept as correct although it is not equivalent.

## Please

- Grade each item on its own; you can revisit an item from the sidebar, and submitting again replaces your
  earlier grade.
- Questions about the task go to the coordinator, not to the other reviewers. Please do not share the folder:
  the problems are the benchmark's private test set until the paper is published.
