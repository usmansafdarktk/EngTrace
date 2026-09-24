"""X5 - E2's flagging threshold, fitted on one half of the traces and reported on the other.

    python evaluator_pilot_17092026/analysis/prm_threshold.py [--labels DIR]

E2 calls a step flagged when its process reward falls below 0.5. That 0.5 is the number
the PRM model cards print in their own examples; this pilot never chose it and never
checked it against the expert labels. RESULTS_X1 closes by asking for exactly that:
recalibrate it "on a split so the threshold is not fitted and reported on the same
data". A reviewer will ask the same question, and the only answer worth giving is a
number measured on traces the threshold never saw.

So the 300 labelled traces are split in two by a hash of the trace code - the same
deterministic split `analysis/answer_check.py` uses for its tolerance, with no random
seed. The threshold that maximises step F1 against the experts' "incorrect" labels on
one half is reported on the other, both ways round, for each of the three PRMs, over
all traces and over the hard case: steps inside traces whose final answer the experts
call correct (RESULTS_X1 Finding 5), the one place where reasoning is judged apart from
the answer.

Printed:
  1. HELD OUT    precision, recall and F1 on the half the threshold never saw, beside
                 what the stock 0.5 gives on that same half. The gap is the whole value
                 of calibrating, and it carries a 95% bootstrap interval over resampled
                 held-out traces. If there is no gap, that is the answer and it closes
                 the question.
  2. STABILITY   the value chosen on each half and the width of the near-maximum band
                 around it. A threshold at the top of a flat plateau is a claim about
                 the PRM; a threshold on a sharp peak is a claim about 150 traces.
  3. OPERATING   the threshold that buys precision ~0.80 on held-out steps and what
     POINT       recall it costs. Over-flagging is E2's documented weakness (RESULTS_E2:
                 only 33-57% of frontier traces pass with nothing flagged), so this is
                 the operating point a user would actually ask for.
  4. IN SAMPLE   the same fit made and reported on all 300 traces at once, so the
                 optimism the split removes is visible as a number rather than assumed.

Steps are aligned with rewards exactly as `x1_analysis.steps()` aligns them: a trace
whose PRM step count does not match the labels is skipped (2 traces for each Qwen PRM,
both over the context limit), and `not_a_claim` steps are not scored. Nothing here
calls a model or touches a GPU - every reward is already in `scores/e2`.
"""
import os as _os
_ANALYSIS = _os.path.dirname(_os.path.abspath(__file__))
_PILOT = _os.path.dirname(_ANALYSIS)

import argparse
import bisect
import hashlib
import json
import random
import sys

sys.path[:0] = [_os.path.join(_PILOT, 'annotation')]
import score_against_labels as S  # noqa: E402

STOCK = 0.5
GRID = [i / 1000.0 for i in range(1, 1000)]     # every threshold to 0.001
SCORED = ('correct', 'alternative_correct', 'incorrect')
BAND = 0.01          # "as good as the peak" for the plateau width
OFFSETS = (-0.10, -0.05, 0.05, 0.10)
TARGETS = (0.80, 0.60, 0.50)
MIN_FLAG = 10        # a threshold must flag this many steps to claim a precision
B = 2000             # bootstrap resamples, as x1_analysis uses


def half(code):
    """The deterministic split `analysis/answer_check.py` uses: a hash of the trace
    code, no random seed, so the two halves are the same on every run and every
    machine."""
    return int(hashlib.blake2b(code.encode(), digest_size=4).hexdigest(), 16) % 2


def counts(rows, t):
    """(TP, FP, FN) at threshold `t` for one trace's steps."""
    tp = fp = fn = 0
    for v, y in rows:
        if v < t:
            tp, fp = tp + y, fp + (not y)
        else:
            fn += y
    return tp, fp, fn


class Curve:
    """One PRM's whole precision/recall curve over one set of steps.

    A step is flagged when its reward is below the threshold, so the rewards themselves
    are the only decision boundaries. Sorting once gives the exact contingency table at
    any threshold from a prefix count, which is what sweeping 999 thresholds over
    twelve curves needs. The steps stay grouped by trace as well, because that is the
    unit the bootstrap resamples.
    """

    def __init__(self, groups):
        self.groups = groups
        rows = sorted((r for g in groups for r in g), key=lambda r: r[0])
        self.v = [v for v, _ in rows]
        self.cum = [0]
        for _, y in rows:
            self.cum.append(self.cum[-1] + (1 if y else 0))
        self.n, self.bad, self.traces = len(rows), self.cum[-1], len(groups)
        self.base = self.bad / self.n if self.n else float('nan')

    def at(self, t):
        """precision, recall, F1 and how many steps are flagged at threshold `t`."""
        i = bisect.bisect_left(self.v, t)       # steps [0, i) have reward < t
        tp = self.cum[i]
        p, r, f = S.prf(tp, i - tp, self.bad - tp)
        return p, r, f, i

    def best(self):
        """The threshold maximising F1, taken as the midpoint of the widest run of grid
        points that attain the maximum - the middle of the winning plateau rather than
        whichever reward happens to sort first. Returns the threshold, its F1 and the
        grid indices of that run."""
        f1 = [self.at(t)[2] for t in GRID]
        top = max([x for x in f1 if x == x] or [float('nan')])
        if top != top:
            return float('nan'), float('nan'), 0, 0
        lo = hi = None
        i = 0
        while i < len(GRID):
            if f1[i] == top:
                j = i
                while j + 1 < len(GRID) and f1[j + 1] == top:
                    j += 1
                if lo is None or j - i > hi - lo:
                    lo, hi = i, j
                i = j + 1
            else:
                i += 1
        return (GRID[lo] + GRID[hi]) / 2, top, lo, hi

    def band(self, top, lo, hi):
        """The unbroken stretch of thresholds around the peak whose F1 stays within
        BAND of it. Wide means the choice of threshold hardly matters; narrow means the
        chosen value is carrying the result."""
        a, b = lo, hi
        while a > 0 and self.at(GRID[a - 1])[2] >= top - BAND:
            a -= 1
        while b + 1 < len(GRID) and self.at(GRID[b + 1])[2] >= top - BAND:
            b += 1
        return GRID[a], GRID[b]

    def at_precision(self, target):
        """The largest threshold whose precision reaches `target` - largest, because a
        higher threshold flags more steps and so recalls more. At least MIN_FLAG steps
        must be flagged, so one lucky flag cannot claim precision 1.0."""
        for t in reversed(GRID):
            p, _, _, n = self.at(t)
            if n >= MIN_FLAG and p == p and p >= target:
                return t
        return None

    def boot_gain(self, t, seed=1):
        """95% interval for F1(t) - F1(0.5) on this set of steps, resampling TRACES.
        Steps inside a trace are not independent, and whether calibrating is worth
        anything is the one claim this script rests on, so it gets an interval."""
        pre = [(counts(g, t), counts(g, STOCK)) for g in self.groups]
        rnd, out, n = random.Random(seed), [], len(pre)
        for _ in range(B):
            a, s = [0, 0, 0], [0, 0, 0]
            for _ in range(n):
                ca, cs = pre[rnd.randrange(n)]
                for k in (0, 1, 2):
                    a[k] += ca[k]
                    s[k] += cs[k]
            out.append(S.prf(*a)[2] - S.prf(*s)[2])
        out = sorted(x for x in out if x == x)
        if not out:
            return float('nan'), float('nan')
        return out[int(.025 * len(out))], out[int(.975 * len(out)) - 1]


def scored_steps(truth, per, codes, prm):
    """Each trace's (reward, the experts call this step incorrect) pairs, kept grouped."""
    groups = []
    for c in codes:
        pr = (per.get(c) or {}).get(prm)
        if pr is None or len(pr) != len(truth[c]['steps']):
            continue
        groups.append([(v, s['label'] == 'incorrect')
                       for s, v in zip(truth[c]['steps'], pr) if s['label'] in SCORED])
    return Curve(groups)


def prepare(truth, per, codes):
    """Per PRM: the whole-subset curve, and per fit half the curve the threshold is
    chosen on, the curve it is reported on and the threshold itself."""
    out = {}
    for p in S.PRMS:
        out[(p, 'all')] = scored_steps(truth, per, codes, p)
        for fit in (0, 1):
            cf = scored_steps(truth, per, [c for c in codes if half(c) == fit], p)
            ch = scored_steps(truth, per, [c for c in codes if half(c) != fit], p)
            out[(p, fit)] = (cf, ch) + cf.best()
    return out


def held_out(fit, title):
    print('\n  %s' % title)
    print('  %-8s %3s %4s %6s %5s %4s %5s %6s  %-20s  %-20s %7s  %s'
          % ('PRM', 'fit', 'held', 'traces', 'steps', 'bad', 'base', 't*',
             'held out at t*', 'held out at 0.500', 'dF1', '95% CI'))
    print('  %-8s %3s %4s %6s %5s %4s %5s %6s  %6s %6s %6s  %6s %6s %6s'
          % ('', '', '', '', '', '', '', '', 'prec', 'rec', 'F1', 'prec', 'rec', 'F1'))
    for p in S.PRMS:
        for f in (0, 1):
            _, ch, t, _, _, _ = fit[(p, f)]
            pa, ra, fa, _ = ch.at(t)
            p0, r0, f0, _ = ch.at(STOCK)
            lo, hi = ch.boot_gain(t)
            print('  %-8s %3d %4d %6d %5d %4d %5.3f %6.3f  %6.3f %6.3f %6.3f  %6.3f %6.3f %6.3f'
                  ' %+7.3f  (%+.3f, %+.3f)%s'
                  % (p, f, 1 - f, ch.traces, ch.n, ch.bad, ch.base, t,
                     pa, ra, fa, p0, r0, f0, fa - f0, lo, hi, '' if lo <= 0 <= hi else ' *'))


def stability(fit, title):
    print('\n  %s' % title)
    print('  %-8s %4s %6s %7s  %-18s %6s  %s'
          % ('PRM', 'fit', 't*', 'F1 t*', 'F1 within %.2f' % BAND, 'width',
             ' '.join('%8s' % ('F1 %+.2f' % o) for o in OFFSETS)))
    for p in S.PRMS:
        for f in (0, 1):
            cf, _, t, top, lo, hi = fit[(p, f)]
            a, b = cf.band(top, lo, hi)
            off = [cf.at(min(max(t + o, GRID[0]), GRID[-1]))[2] for o in OFFSETS]
            print('  %-8s %4d %6.3f %7.3f  %6.3f - %-9.3f %6.3f  %s'
                  % (p, f, t, top, a, b, b - a, ' '.join('%8.3f' % x for x in off)))


def operating(fit, title):
    print('\n  %s' % title)
    print('  %-8s %4s %6s %6s  %-15s  %s'
          % ('PRM', 'fit', 'target', 't', 'on the fit half', 'on the held-out half'))
    print('  %-8s %4s %6s %6s  %7s %7s  %7s %7s %7s %7s'
          % ('', '', 'prec', '', 'prec', 'rec', 'prec', 'rec', 'F1', 'flagged'))
    for p in S.PRMS:
        for f in (0, 1):
            cf, ch, _, _, _, _ = fit[(p, f)]
            for target in TARGETS:
                t = cf.at_precision(target)
                if t is None:
                    reach = []
                    for u in GRID:
                        pu, _, _, nu = cf.at(u)
                        if nu >= MIN_FLAG and pu == pu:
                            reach.append((pu, u))
                    bp, bt = max(reach) if reach else (float('nan'), float('nan'))
                    print('  %-8s %4d %6.2f %6s  never reached (best on the fit half %.3f, at t=%.3f)'
                          % (p, f, target, '-', bp, bt))
                    continue
                fp_, fr_, _, _ = cf.at(t)
                hp, hr, hf, hn = ch.at(t)
                print('  %-8s %4d %6.2f %6.3f  %7.3f %7.3f  %7.3f %7.3f %7.3f %7d'
                      % (p, f, target, t, fp_, fr_, hp, hr, hf, hn))


def in_sample(fit, title):
    print('\n  %s' % title)
    print('  %-8s %5s %5s %6s %7s %7s %7s  %7s %9s %9s'
          % ('PRM', 'steps', 'bad', 't*', 'prec', 'rec', 'F1', 'F1 @0.5', 'held out', 'optimism'))
    for p in S.PRMS:
        c = fit[(p, 'all')]
        t, _, _, _ = c.best()
        pa, ra, fa, _ = c.at(t)
        held = [fit[(p, f)][1].at(fit[(p, f)][2])[2] for f in (0, 1)]
        mean_held = sum(held) / 2
        print('  %-8s %5d %5d %6.3f %7.3f %7.3f %7.3f  %7.3f %9.3f %+9.3f'
              % (p, c.n, c.bad, t, pa, ra, fa, c.at(STOCK)[2], mean_held, fa - mean_held))


def verdict(fits):
    print('\nVERDICT - does calibrating the threshold buy anything, held out?')
    print('  %-44s %8s %8s %9s %14s'
          % ('subset / PRM', 'F1 t*', 'F1 0.500', 'gain', 't* (half 0, 1)'))
    for title, fit in fits:
        for p in S.PRMS:
            cal = [fit[(p, f)][1].at(fit[(p, f)][2])[2] for f in (0, 1)]
            stock = [fit[(p, f)][1].at(STOCK)[2] for f in (0, 1)]
            ts = [fit[(p, f)][2] for f in (0, 1)]
            print('  %-44s %8.3f %8.3f %+9.3f %7.3f %6.3f'
                  % ('%s / %s' % (title.split(' (')[0].lower(), p),
                     sum(cal) / 2, sum(stock) / 2, (sum(cal) - sum(stock)) / 2, ts[0], ts[1]))
    print('\n  gain is the mean over both directions of the split; section 1 carries the interval')
    print('  on each direction separately, and marks it * when it excludes zero. A gain near zero')
    print('  means the stock 0.5 already sits on the flat top of the F1 curve and E2 needs no')
    print('  recalibration - which says nothing about whether the PRM is any good, only that this')
    print('  knob is not what is wrong with it.')


def load(labels_dir):
    truth, _ = S.build_truth(S.read_labels(labels_dir), S.read_consensus(labels_dir))
    keyfile = {r['code']: (r['model_key'], r['item_id'])
               for r in map(json.loads, open(_os.path.join(S.TASKS, 'keyfile.jsonl'), encoding='utf-8'))}
    e2 = S.current('e2')
    if not e2:
        raise SystemExit('no E2 rows in scores/e2 - run the evaluator first')
    per = {c: {p: e2[keyfile[c]]['meta']['prm'][p]['step_probs'] for p in S.PRMS}
           for c in truth if keyfile[c] in e2}
    return truth, per


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--labels', default=_os.path.join(_PILOT, 'experts_filled_labels', 'version_2', 'labels'))
    a = ap.parse_args()
    truth, per = load(a.labels)
    right = sorted(c for c in truth if truth[c]['final_answer'] == 'correct')
    n0 = sum(1 for c in truth if half(c) == 0)
    print('expert ground truth: %d traces, %d of them with a correct final answer'
          % (len(truth), len(right)))
    print('split on blake2b(code) mod 2: half 0 holds %d traces, half 1 holds %d - no random seed'
          % (n0, len(truth) - n0))
    print('a step is flagged when its reward is below the threshold; E2 ships %.1f' % STOCK)
    print('base = the share of scorable steps the experts call incorrect, the precision floor')

    sets = [('ALL TRACES', sorted(truth)),
            ('INSIDE CORRECT-ANSWER TRACES (the hard case, RESULTS_X1 Finding 5)', right)]
    fits = [(title, prepare(truth, per, codes)) for title, codes in sets]

    print('\n1. HELD OUT - the threshold maximises F1 on the fit half and is reported on the other')
    for title, fit in fits:
        held_out(fit, title)
    print('\n  * = the 95%% interval on dF1 excludes zero (%d resamples of the held-out traces)' % B)

    print('\n2. STABILITY - the chosen value, and how flat F1 is around it on the half it was fitted on')
    for title, fit in fits:
        stability(fit, title)

    print('\n3. THE OPERATING POINT - precision bought on the fit half, paid for on the held-out half')
    print('   (a threshold must flag at least %d steps on the fit half to claim a precision)' % MIN_FLAG)
    for title, fit in fits:
        operating(fit, title)

    print('\n4. IN SAMPLE, for reference - fitted and reported on all the labelled traces at once')
    print('   optimism = this in-sample F1 minus the mean held-out F1 from section 1')
    for title, fit in fits:
        in_sample(fit, title)

    verdict(fits)


if __name__ == '__main__':
    main()
