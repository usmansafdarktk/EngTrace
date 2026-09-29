# EngTrace paraphrase check: the guide

## What this is

The benchmark tests whether a model's answer survives a rewording of the problem. Each problem in
your kit was rewritten by a language model. It was told to keep every number, unit, symbol and
technical term, and the order of the parts, and to change only the wording. A script has already
checked the numbers, the units, the symbols and the order of the parts.

You check what a script cannot: that the rewritten version is still the same problem, with the same
answer. Your kit holds up to 30 items from your own branch, whole templates, about one or two
minutes each.

## What you see

For each item:
- the original problem and the rewritten version, side by side;
- the original's answer;
- the original's full solution, one click away.

All of it is shown as plain text, exactly as a model receives it.

## The three questions

**1. Is the rewritten version the same problem?** Answer yes when it gives the same data, asks for the
same quantities and states the same conditions and assumptions. Answer no when anything given, asked
for or assumed has changed, gone missing or been added, even a single word. Words like these matter:
"gauge", "absolute", "steady", "per unit length", "initially", "at least", "approximately",
"neglect friction".

**2. Does the original's answer still answer it, exactly?** Answer yes when a correct solution of the
rewritten version gives exactly the original's answer: the same values, units and classification.
Answer no when any part of the answer would differ.

**3. Does the rewritten version add an ambiguity, an error or a hint?** Answer yes in three cases:
- it can reasonably be read in two ways where the original cannot;
- it contains a mistake the original does not, such as a wrong term or grammar that changes the meaning;
- it hints at the method or the answer, for example by naming a formula the original leaves to the
  solver.

**Keep the pair when the answers are yes, yes, no.** Any other answer needs a short note: what
changed, and where.

## What does not count as a change

Different wording, a different order of sentences, and a question turned into an instruction ("Find
..." for "What is ...?") are fine. Judge the rewritten version against the original, not against
what the original should have said. If the original itself looks wrong, say so in the note, but still
answer the questions about the rewrite.

## Please

Work alone, and do not discuss items with the other reviewers until everyone has submitted. The app
saves after every item, so you can stop and continue later. When it says you are done, send back the
file `<your id>.jsonl` from the folder you ran the app in.
