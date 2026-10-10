# SWITCH W1B: coverage variants at matched settings

Script: `full_run_28092026/coverage_variants.py` (only file edited). Outputs: `results/coverage_variants.json`, `results/sections/coverage_variants.md` (new "Matched settings" subsection).

## The `matched` key
`matched = {templates: {reading: n}, models: {model: {reading: summary, store, route_adjusted_left_out, sources}}, pairs: {reading: [55 pairs]}, separation_rule: {claude_vs_deepseek: {...}}, verbosity: {...}, reasoning_tokens: {model: {...} | null}, stores: {model: store}}`. Inner shapes match the default family's keys. The three re-run models read `reasoning-medium-full`; the other eight read `main`, and their summaries are asserted equal to their default values.

## Results at matched settings (148 / 148 / 148 / 127 templates)
| reading | gemma-4-26b-a4b | gpt-5.4-mini | gemini-3.1-flash-lite | pairs separated (of 55) | Claude − DeepSeek, Holm p |
|---|---|---|---|---|---|
| as scored | 0.884 | 0.892 | 0.881 | 15 (default 28) | +0.023, p = 0.069 |
| matching alone | 0.850 | 0.847 | 0.855 | 20 (26) | +0.042, p = 0.0005 |
| route-adjusted | 0.959 | 0.972 | 0.953 | 17 (29) | +0.003, p = 1.0 |
| intermediate-only | 0.813 | 0.817 | 0.815 | 9 (18) | +0.037, p = 0.032 |

- **The separation does not hold under every reading.** It also no longer holds as scored: the Holm p is 0.049 at default and 0.069 at matched settings. The difference is the same; the family changed. It holds only under matching alone and intermediate-only.
- **What the judge adds over matching alone:** 0.019 (glm-5.3) to 0.045 (gpt-5.4-mini); at default the range is 0.019 to 0.038.
- **Verbosity slope at matched settings:** +0.0020 MC per 10 numeric values [+0.0005, +0.0041], which excludes zero; at default it is +0.0015 [−0.0001, +0.0034]. Per 1,000 visible tokens: −0.0010 [−0.0114, +0.0076]. Within model and template: +0.0025 [+0.0005, +0.0049]. The fit covers 23,862 responses.
- **Reasoning-token quartiles:** identical to the default key, because the default already reads each re-run model at its reasoning store.

## Default family unchanged
I diffed the JSON against the pre-run copy. Every existing key is byte-identical: `quick`, `draws`, `rule`, `keys`, `templates`, `models`, `pairs`, `separation_rule`, `verbosity` and `reasoning_tokens`. Only the provenance fields `git`, `script_sha256_lf` and `written_at_utc` changed. The default part of the markdown section is identical too, and `--render` reproduces the whole section from the JSON.

**Runtime:** 2 min 43 s for the full pass, both families. **Self-test:** passes, including the new matched-family block on hand-made data.

## Not computed
Nothing was left out.
