# Phase C3 — Reviewer H (domain: thermochemistry / materials / fluids), ELECTRICAL branch

Frozen ref `2fe50c72153892d0b1fa67e7ed486ebb36abca0f`. Working tree clean at that commit.
Scope: whether the NUMBERS in `data/templates/branches/electrical_engineering/constants.py`
are right, and whether each row's provenance CLASS is the honest one. Tag syntax, `@domain`
presence and register completeness are Reviewer G's; the other four branches are the other
H reviewers'; the `MEDIA_VELOCITIES` optical-vs-radio excursion is already registered in
`tests/constants_integrity/domain_findings.txt` and is not re-derived here.

Method: I did **not** use `tests/constants_integrity/refractiveindex.py`. I wrote my own
evaluator of the RefractiveIndex.INFO dispersion formulas from their published definitions
and read the YAML members straight out of the archive. The implementation is anchored
externally: it reproduces SCHOTT's own published catalogue d-line indices to six decimals
(N-BK7 nd = 1.516800 vs catalogue 1.51680; N-SF2 nd = 1.647690 vs catalogue 1.64769), so a
transcription error in my formulas would have shown up there before it reached any row.

## 1. Verdict

**PASS WITH FINDINGS** — every value I re-derived is correct to every digit the tag prints
and all four claim classes are honest, but the `tol=` class is fitted to its own residual
and the two gas rows are judged by a refractivity comparison that is condition-corrected
for neither (5 findings, none blocking).

## 2. Independent re-derivation

29 datasets evaluated at 0.589 um; the 8-row mandate is exceeded. "Committed" is the number
written in the tag or group comment.

| Row / claim | Dataset consulted | My number | Committed | Agree? |
|---|---|---|---|---|
| Helium (defect) | `main/He/nk/Ermolov.yml` — **formula 2** | 1.0000348831 | 1.0000349 | yes |
| Helium (defect) | `main/He/nk/Mansfield.yml` — **formula 6** | 1.0000349109 | 1.0000349 | yes |
| Helium Δ on n | vs row 1.000036 | +1.117e-06 | +1.1e-06 | yes |
| Helium Δ on n−1 | vs Ermolov / Mansfield | +3.202% / +3.120% | +3.2% | yes (Ermolov) |
| GaP (defect) | `main/GaP/nk/Adachi.yml` (tab, interp 0.58591–0.59424) | 3.40685 | 3.4068 | yes |
| GaP | `Aspnes.yml` (tab, interp 0.5636–0.5904) | 3.37740 | 3.3774 | yes |
| GaP | `Bond.yml` (tab n, interp 0.5–0.6) | 3.36160 | 3.3616 | yes |
| GaP | `Jellison.yml` (tab, interp 0.585–0.590) | 3.37680 | 3.3768 | yes |
| GaP | `Khmelevskaia.yml` (tab, interp 0.588–0.589) | 3.39155 | 3.3916 | yes |
| GaP spread vs 3.5 | — | +2.74% … +4.12% | +2.7% to +4.1% | yes |
| Air | `other/mixed gases/air/nk/Borzsonyi.yml` — formula 2 | 1.0002882101 | 1.0002882 | yes |
| Air alt | `Ciddor.yml` — formula 6, 288.15 K | 1.0002771520 | 1.0002772 | yes |
| Carbon Dioxide | `main/CO2/nk/Bideau-Mehu.yml` — formula 6 | 1.0004488812 | 1.0004489 | yes |
| Water | `main/H2O/nk/Daimon-20.0C.yml` — formula 2 | 1.3333583993 | 1.3333584 | yes |
| Ethanol | `organic/C2H6O - ethanol/nk/Kedenburg.yml` | 1.3615136671 | 1.3615137 | yes |
| Glycerine | `organic/C3H8O3 - glycerol/nk/Gupta.yml` — **formula 5** | 1.4712766417 | 1.4712766 | yes |
| Glycerine alt | `Birkhoff.yml` (tab, **interp** 0.41328–0.61992) | 1.4714963 | 1.4714963 | yes |
| Benzene | `organic/C6H6 - benzene/nk/Chang.yml` — formula 1 | 1.4995124717 | 1.4995125 | yes |
| Benzene alt | `Moutzouris.yml` | 1.4956337 | 1.4956337 | yes |
| Carbon Disulfide | `main/CS2/nk/Chemnitz.yml` | 1.6281390184 | 1.6281390 | yes |
| Ice | `main/H2O/nk/Warren-2008.yml` (tab, **interp** 0.58–0.59) | 1.30973 | 1.30973 | yes |
| Ice brackets | same | 1.31 @0.58, 1.3097 @0.59 | same | yes |
| Fused Silica | `main/SiO2/nk/Malitson.yml` — formula 1 | 1.4584132065 | 1.4584132 | yes |
| Crown Glass | `specs/schott/optical/N-BK7.yml` — formula 2 | 1.5167401205 | 1.5167401 | yes |
| Polystyrene | `organic/(C8H8)n - polystyrene/nk/Sultanova.yml` | 1.5914840706 | 1.5914841 | yes |
| Flint Glass | `specs/schott/optical/N-SF2.yml` — formula 2 | 1.6475512458 | 1.6475512 | yes |
| Sapphire (o) | `main/Al2O3/nk/Malitson-o.yml` — formula 1 | 1.7680928242 | 1.7680928 | yes |
| Sapphire (e) | `main/Al2O3/nk/Malitson-e.yml` | 1.7600181802 | 1.7600182 | yes |
| Diamond | `main/C/nk/Peter.yml` — formula 1 | 2.4172981917 | 2.4172982 | yes |
| Diamond alt | `Phillip.yml` (tab, interp) | 2.4166020 | 2.4166020 | yes |
| `C0` | `codata_2022/allascii.txt` | 299 792 458 (exact) | 299792458 | yes |
| `EPSILON_0` | same, 8.854 187 8188e-12 → 4 s.f. | 8.854e-12, −2.121e-05 rel | −2.12e-5 | yes |

**No disagreement at any printed digit, in any row.** Every "also rounds to" remark in the
group comments is true and every "does not" is true (checked: CS2 Chang 1.623394 → 1.623,
correctly excluded; Rheims glycerol formula starts at 0.5893 um and is genuinely refused at
0.589).

**Class judgments on the four rows the brief singles out.**

- **Helium `[KNOWN-DEFECTIVE]` — correct, and it survives the strongest attack I could
  mount.** 1.000036 is not a rounding of either dataset at its own 1e-06 precision (both
  round to 1.000035). More decisively, the discrepancy cannot be explained away by
  conditions: to make 1.000036 right you would need 104 569 Pa at 273.15 K, or 264.68 K at
  101 325 Pa, and both datasets are already at exactly the group comment's 273.15 K /
  101 325 Pa. Invisible in n (1.1e-06, below every display in the corpus) but 3.2% in n−1 —
  and n−1 is the physically meaningful quantity for a gas, being the density-proportional
  one (Gladstone–Dale). Calling that "defective" rather than absorbing it in a tolerance is
  the honest call for a domain reader, even though no emitted item can see it.
- **GaP `[KNOWN-DEFECTIVE]` — correct.** All five datasets covering 589 nm land in
  3.3616–3.4068 and **every one rounds to 3.4 at the row's own 1-dp precision**, so the row
  fails even at its own stated coarseness. The sixth GaP member, `Parsons.yml`, is far-IR
  (54–366 um) and its exclusion from "all five … that cover 589 nm" is exactly right.
- **The `[UNVERIFIED]` rows — all four honest.** PTFE: no PTFE/teflon/C2F4 dataset exists;
  the nearest fluoropolymer, `(C2ClF3)n`, indeed has no row bracketing 589 nm (the other
  fluorine members are small liquids — perfluorohexane 1.2516, trifluoroacetic anhydride
  1.2694 — not PTFE). Polyethylene: both `(C2H4)n` datasets are genuinely out of band
  (David 1.99996–12.26 um, spin-coated LLDPE with sub-2 um data removed by the authors;
  Smith 40–200 um HDPE). "Glass (amorphous semiconductor)": exactly **18** amorphous
  datasets cover 589 nm, spanning **1.7400** (`other/amorphous/HAC/nk/Smith-Td250C.yml`) to
  **3.8830** (`main/Ge2Sb2Te5/nk/Frantz-amorphous.yml`) — the claim is exact to both
  endpoints and both member names. Tissue: 7.14² = 50.98 and √51 = 7.1414, so the row is
  demonstrably √ε_r and not an index; the archive's human-body datasets are all optical,
  none is muscle, and they span 1.3368–1.5825, which rounds to the claimed 1.34–1.58.
- **Family rows and Sapphire — honest, because the comment names the chooser.** N-BK7 is
  the canonical crown and 1.5167401 → 1.52 at the row's 3 s.f.; N-SF2 is a genuine dense
  flint (SF = Schwerflint) and 1.6475512 → 1.65. Both comments say in terms that the member
  was "chosen by the implementer, not by the row" and refer the question onward, which is
  the disclosure that keeps a family row from overstating its evidence. Sapphire likewise
  states the ray, gives the extraordinary value (1.7600182 → 1.76), and picks the ordinary
  ray — which is the ray conventional tables quote for an unqualified sapphire index. I
  would not call any of these dishonest. See F-5 for the residual domain caveat.

## 3. Findings

### F-1 CONFIRMED — the Air row's "+1.7% on the refractivity" is an uncorrected cross-condition comparison; the row is ~4x better than its own tag says

`Borzsonyi.yml` is at 273 K / **100 000 Pa**, which the tag itself notes is not the group
comment's 1 atm. Comparing n−1 across a 1.3% pressure difference without scaling is not a
comparison a domain reader should accept. Scaling by Gladstone–Dale (n−1 ∝ ρ ∝ P/T) to
273.15 K / 101 325 Pa gives n = 1.000291869, against the row's 1.000293 — a residual of
1.13e-06 on n and **+0.387%** on n−1, not +1.7%. The three independent air datasets
(Peck, Birch, Ciddor, all 1.000277 at 288.15 K) scale to ≈1.000292 at 0 °C and agree.

```
python -c "import zipfile,yaml,math;Z=zipfile.ZipFile('docs/references/refractiveindex_info/refractiveindex.info-database-main.zip');R='refractiveindex.info-database-c5c2f188e848453def5970e347399d653df2ffc2/';d=yaml.safe_load(Z.read(R+'database/data/other/mixed gases/air/nk/Borzsonyi.yml'));C=[float(x) for x in str(d['DATA'][0]['coefficients']).split()];l2=0.589**2;n=math.sqrt(1+C[0]+sum(C[i]*l2/(l2-C[i+1]) for i in range(1,len(C)-1,2)));r=(n-1)*(101325/100000)*(273/273.15);print(n, 1+r, 1.000293-(1+r), (2.93e-4/r-1)*100)"
```

**Impact:** documentation-honesty, not a wrong value — the row is right and the tag
understates how right it is. It matters only because the *same* refractivity framing, used
one row later, is what condemns Helium. A reader comparing the two tags sees "+1.7%
tolerated, +3.2% defective" and will infer a threshold that does not exist; the real
distinction is that He's gap survives condition-correction and air's does not.

### F-2 CONFIRMED — every `tol=` in the branch is the observed residual rounded up, so the tol class cannot fail

| Row | residual I measured | `tol=` written |
|---|---|---|
| Air | 0.000479% | 0.0005% |
| Glycerine | 0.11705% | 0.12% |
| Benzene | 0.09918% | 0.1% |

Three for three, the tolerance is the smallest round number above the discrepancy it must
admit. A tolerance chosen after seeing the residual is an accommodation, not an uncertainty
budget: it is unfalsifiable by construction, and it silently converts "how far off is this?"
into "it passed". Reproduce with the per-row numbers in §2 (e.g. `1.473/1.4712766417 - 1`).

**Impact:** moderate and structural rather than numeric. It is the one mechanism in the
branch by which a future wrong value could be absorbed rather than surfaced — precisely the
failure mode `[KNOWN-DEFECTIVE]` exists to prevent. The C3 rule quoted in the Helium tag
("tol= is only for an artefact whose conditions differ from the table's") is the right
instinct but constrains *when* a tolerance may be used, never *how large* it may be.

### F-3 PLAUSIBLE — the Helium tag names 2 of the 3 datasets that cover 589 nm

Five He members exist; `Cuthbertson.yml` (0.2753–0.5462 um) and `Smith.yml` (0.168–0.288
um) are out of band and `Shaw.yml` is an n2 (nonlinear) file, but **`Borzsonyi.yml` covers
589 nm** and gives n = 1.000034387 at 273 K / 100 000 Pa. The tag does not mention it.

```
python -c "import zipfile,yaml,math;Z=zipfile.ZipFile('docs/references/refractiveindex_info/refractiveindex.info-database-main.zip');R='refractiveindex.info-database-c5c2f188e848453def5970e347399d653df2ffc2/';d=yaml.safe_load(Z.read(R+'database/data/main/He/nk/Borzsonyi.yml'));C=[float(x) for x in str(d['DATA'][0]['coefficients']).split()];l2=0.589**2;print(math.sqrt(1+C[0]+sum(C[i]*l2/(l2-C[i+1]) for i in range(1,len(C)-1,2))), d.get('CONDITIONS'))"
```

Scaled to 101 325 Pa it is ≈1.0000348, so it **does not rescue the row** and the
`[KNOWN-DEFECTIVE]` class is untouched. The finding is one of asymmetric thoroughness: the
GaP tag one screen below says "all five on-disk datasets that cover 589 nm" and is exactly
right, while the He tag names two without claiming completeness and leaves a third
unmentioned. **Impact:** low — cosmetic, but it is the difference between a defect claim a
reader can audit and one they must re-derive.

### F-4 PLAUSIBLE — Glycerine's `tol=` is applied to a dataset that states no conditions at all

`Gupta.yml` has `CONDITIONS: {}` — no temperature, no composition. The C3 rule the Helium
tag cites licenses `tol=` "only for an artefact whose conditions differ from the table's";
unknown conditions are not differing conditions, and the tag's own parenthetical ("no
temperature stated in the file") concedes the point. Domain-wise this is the row where it
matters most: glycerol is strongly hygroscopic and its index is a steep function of water
content (anhydrous ≈1.4746 at 20 °C, falling ~0.001 per ~0.5 wt% water), so a 0.12%
tolerance is absorbing a *composition* ambiguity — which no tolerance on a pure-substance
citation should do — rather than a temperature offset. The same "no temperature stated"
concession appears on Kedenburg (CS2), Arosa (SiO2), Zhang (polystyrene) and Phillip
(diamond), but those rows use `precision=` and are unaffected. **Impact:** low-moderate;
the value is right, the warrant is weaker than the tag implies.

### F-5 PLAUSIBLE (domain) — the tissue row is not merely the wrong band, it is the wrong microwave quantity

The `[UNVERIFIED]` tag's claim that 7.14 is √ε_r and not an optical index is correct and
well argued. But √ε_r is the **lossless** phase index, and muscle is a conspicuously lossy
medium. Using Gabriel-type tissue parameters, the true phase index
n = √(ε′/2 · (√(1+tan²δ)+1)):

| f | ε′ | σ (S/m) | tan δ | lossy phase index | row's √ε_r |
|---|---|---|---|---|---|
| 3 GHz (the row's own label) | 52.0 | 2.14 | 0.247 | **7.27** | 7.14 |
| 500 MHz (consumer's top) | 56.9 | 1.00 | 0.632 | **7.88** | 7.14 |
| 50 MHz (consumer's bottom) | 72.2 | 0.69 | 3.436 | **12.86** | 7.14 |

So even granting the row its microwave reading and its own 3 GHz label, it is ~1.8% high on
velocity; at the band the single consumer actually draws it is high by up to 80%, and at
50 MHz the medium is conduction-dominated (tan δ > 3) where a phase velocity from √ε_r is
not meaningful at all. **Impact:** low for the pool (the row is already `[UNVERIFIED]` and
already inside the registered D-032 excursion), but it should be recorded, because the
obvious "fix" for that excursion — re-sourcing the rows as radio-frequency permittivities —
would reproduce this error unless it carries the loss term. A re-sourcing that writes
√ε_r(50 MHz) is still wrong by 34%.

## 4. Falsification attempts that failed

- **Cherry-picking sweep.** For all 11 cited material directories I evaluated *every*
  sibling dataset at 589 nm (~90 datasets) and checked whether a member closer to the row's
  value was passed over, or whether a chosen member is an outlier. Nothing. The picks are
  the defensible ones throughout, and several are actively discriminating: `main/SiO2`
  contains crystalline quartz (Ghosh-o/e, 1.544/1.553) alongside fused silica and the
  "Fused Silica" row correctly takes Malitson; `main/C` contains 28 carbon datasets from
  graphite and nanotubes (1.33) to diamond (2.77) and the "Diamond" row correctly takes
  Peter, whose own COMMENTS read "Cubic carbon (diamond)".
- **The Helium defect.** Tried to dissolve it three ways — a third on-disk dataset (F-3, it
  does not), a rounding at the row's own precision (both round to 1.000035), and a
  conditions story (would need 104.6 kPa or 264.7 K against a file that already states
  273.15 K / 101 325 Pa). The class holds.
- **The GaP defect.** Tried the sixth member (`Parsons`, far-IR, correctly excluded) and
  tried rounding at the row's own 1 dp (all five give 3.4, not 3.5). Holds.
- **The `[UNVERIFIED]` rows as concealment.** For each I searched the archive independently
  for a dataset the implementer might have missed: PTFE by five name patterns, polyethylene
  by band, amorphous by member scan, tissue by human-body scan. In all four cases the
  archive really does lack what the row needs, and in two cases the tag's stated numbers
  (18 datasets; 1.74–3.88; 1.34–1.58) reproduce **exactly**.
- **Interpolation-dependence.** The three rows resting on tabulated data were checked at
  both bracketing rows, not just the interpolant: Ice (1.31 @0.58 and 1.3097 @0.59, both
  1.31 at 3 s.f. — the tag's claim is true), Phillip diamond (2.4209/2.4114, both 2.42),
  and Birkhoff glycerol, which interpolates across 0.41328–0.61992 between 1.48 and 1.47 so
  its 3 s.f. value *would* depend on the interpolation — and the tag says exactly that and
  declines to cite it. That is the most careful single judgment in the branch.
- **The evaluator as a claim under test.** My from-scratch implementation agrees with every
  committed digit, and independently reproduces SCHOTT's published catalogue nd for both
  glasses to 6 dp and the literature n_o(589) = 1.7682 for sapphire. I found no formula
  transcription error. (I did not exercise formula 4/7/8/9 — no electrical row uses them.)
- **Suite corroboration.** `python -m tests.constants_integrity.test_citations_resolve`
  passes: electrical 34 tags, 16 resolved, 0 locator-only, 0 LEGACY, 0 unresolvable.

## 5. Further probing and improvements (separately time-boxed, 15 min)

1. **Cap `tol=`, don't just gate it (F-2).** The C3 rule says when a tolerance is allowed;
   it should also say what may set its size — a stated experimental uncertainty, a spread
   across datasets, or a declared condition offset — and the resolver should reject a
   `tol=` that exceeds, say, 2x the residual it admits. As written, F-2's three tolerances
   are indistinguishable from curve-fitting.
2. **Condition-normalise gas comparisons (F-1).** Where a gas dataset's P/T differs from
   the table's, the tag should print the Gladstone–Dale-scaled residual alongside the raw
   one. This costs two multiplications, and it is what makes the Helium verdict legible:
   He's gap is the one that *survives* scaling.
3. **Adopt the GaP pattern everywhere (F-3).** Defect tags should enumerate every on-disk
   dataset covering the wavelength, as GaP does, rather than a sufficient subset.
4. **Carry the loss term into any RF re-sourcing (F-5).** If the registered D-032 excursion
   is resolved by re-sourcing rows as radio-frequency properties, the stored quantity must
   be the complex/lossy phase index, not √ε_r, or the lossiest rows will be wrong by tens of
   percent while looking better-sourced than they are.
5. **Consider naming the ray and the grade in the key, not only the comment.** "Sapphire
   (ordinary ray)", "Crown Glass (N-BK7)", "Flint Glass (N-SF2)" would make the
   implementer's pick visible at the point of *use* rather than only at the point of
   citation. Domain caveat that motivates this: N-SF2's nd = 1.64769 sits at the very
   bottom of the dense-flint range (the SF family runs to ~1.95), so "Flint Glass (dense)"
   would resolve just as defensibly against SF11 at ~1.785 — a >8% swing in the row's value
   that is invisible to anyone reading only the key.
6. **Not attempted, worth a later pass:** whether `MEDIA_VELOCITIES`' *ordering* (the
   group comments' "Gases / Liquids / Solids / Special Cases") is load-bearing for any
   template that iterates the dict, and whether the two defective rows are drawn at
   different rates than the rest — both are item-pool questions owned by
   `phaseC3_item_pool_impact.md` and deliberately out of scope here.
