# Checking an arithmetic checker: flagged calculations

Thank you for helping. The workbook `flag_review_reader.xlsx` has one row per calculation. Each row is
one calculation taken from an AI model's written solution to an engineering problem. An automatic checker decided that
the result the model wrote is **wrong at the precision it shows**. We need a person to say whether the
checker is right.

No specialist knowledge of any branch is needed: each row is arithmetic and rounding. It should take a
minute or two per row.

## What to do on each row

1. **Recompute the left side** of `claim` from the numbers shown; a calculator is fine.
   `checker's left value` is what the checker got, so this also shows whether it read the expression
   correctly.
2. **Compare** your result with the number on the right side of `claim`, at the precision in
   `displayed precision`. For example, 0.001 means the number is shown to three decimals, so it should
   equal your result rounded to three decimals.
3. **Choose a verdict** from the list in the `verdict` column:
   - `slip`: the model's arithmetic is wrong at the digit it shows.
   - `checker`: the arithmetic is right and the checker misread it. For example: a unit it did not
     recognise; a value carried over, rounded, from an earlier step; a sentence split in the wrong
     place; a number taken from the wrong side.
   - `unsure`: you cannot tell.
4. **For every `checker` verdict, write briefly why in `note`.** The notes are what the checker will be
   fixed from.

## The columns

- **code**: the row's identifier. Please leave it unchanged: it is how your verdicts are matched back.
- **problem**: the type of problem the calculation comes from.
- **claim**: the calculation as the checker split it, left side = right side, with brackets showing how
  it grouped the terms.
- **checker's left value**: what the checker computed for the left side.
- **right value**: the number the model wrote on the right side.
- **displayed precision**: the place value of the last digit the model shows (0.001 means three decimals).
- **units (left | right)**: the units the checker read on each side, where it found any; `-` means none
  on that side.
- **step**: the full step the claim came from, exactly as the model wrote it.
- **verdict** and **note**: yours to fill in.

The `step` text is the model's raw output, with its maths in LaTeX: `\frac{a}{b}` is a/b, `\sqrt{x}` is
the square root of x, `\times` is multiplication, `\approx` means "approximately" and `x^{2}` is x
squared.

## Returning it

Please leave the other columns as they are and send the workbook back as it is (.xlsx). A row left
blank simply counts as not read.

## Confidential

The material comes from an unreleased test set. Please keep both files to yourself, and delete your
copies once you have sent the workbook back.
