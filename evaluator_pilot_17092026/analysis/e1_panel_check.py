"""Is the proposed E1 panel available on OpenRouter today, and what would E1 cost?

Availability, per judge: listed in the catalogue, not scheduled for expiry, at
least one serving endpoint, structured-output support, and a LIVE call with a
Tribunal-shaped JSON request that parses.

Cost: from what each judge actually used in the discrimination probe
(scores/judge_probe/replies.jsonl - real Tribunal prompts of the same size as
E0's), priced at today's catalogue rates, times the number of traces E0's trigger
sends to judges on the shared seed (the E1 run uses the same trigger).
"""
import os as _os
_ANALYSIS = _os.path.dirname(_os.path.abspath(__file__))
_PILOT = _os.path.dirname(_ANALYSIS)
_REPO = _os.path.dirname(_PILOT)

import glob
import json
import os
import statistics as st
import time
import urllib.request

from dotenv import load_dotenv

load_dotenv(_os.path.join(_REPO, '.env'))
PANEL = [('minimax-m3', 'minimax/minimax-m3'),
         ('mimo-v2.5-pro', 'xiaomi/mimo-v2.5-pro'),
         ('grok-4.6', 'x-ai/grok-4.6')]
H = {'Authorization': 'Bearer ' + os.environ['OPENROUTER_API_KEY'], 'User-Agent': 'engtrace'}


def get(url):
    return json.load(urllib.request.urlopen(urllib.request.Request(url, headers=H), timeout=60))


cat = {m['id']: m for m in get('https://openrouter.ai/api/v1/models')['data']}
ASK = ('Return a raw JSON object only, nothing else: '
       '{"results": [{"step_index": 2, "category": "Calculation Error"}]}')

print('AVAILABILITY (checked %s UTC)\n' % time.strftime('%Y-%m-%d %H:%M', time.gmtime()))
prices = {}
for key, mid in PANEL:
    m = cat.get(mid)
    print('== %s  (%s)' % (key, mid))
    if not m:
        print('   NOT LISTED')
        continue
    pin, pout = float(m['pricing']['prompt']) * 1e6, float(m['pricing']['completion']) * 1e6
    prices[key] = (pin, pout)
    sp = set(m.get('supported_parameters') or [])
    print('   listed      yes; expiration: %s' % (m.get('expiration_date') or 'none announced'))
    print('   price       $%.2f in / $%.2f out per M tokens' % (pin, pout))
    print('   weights     %s' % ('open (%s)' % m['hugging_face_id'] if m.get('hugging_face_id') else 'closed'))
    print('   JSON mode   %s' % ('supported' if {'response_format', 'structured_outputs'} & sp else 'not advertised'))
    eps = get('https://openrouter.ai/api/v1/models/%s/endpoints' % mid)['data']['endpoints']
    for e in eps:
        up = e.get('uptime_last_30m')
        print('   endpoint    %-22s status=%-3s uptime30m=%-7s max_out=%s' % (
            e.get('provider_name'), e.get('status'),
            '%.1f%%' % up if isinstance(up, (int, float)) else up, e.get('max_completion_tokens')))
    body = json.dumps({'model': mid, 'max_tokens': 4096,
                       'messages': [{'role': 'user', 'content': ASK}]}).encode()
    t = time.time()
    try:
        r = json.load(urllib.request.urlopen(urllib.request.Request(
            'https://openrouter.ai/api/v1/chat/completions', data=body,
            headers={**H, 'Content-Type': 'application/json'}), timeout=300))
        txt = (r['choices'][0]['message'].get('content') or '').strip()
        try:
            parsed = json.loads(txt.strip('`').replace('json\n', '', 1))
            ok = parsed.get('results', [{}])[0].get('category') == 'Calculation Error'
        except Exception:
            ok = False
        print('   LIVE CALL   %s  served=%s via %s  %.1fs  JSON %s' % (
            'OK' if txt else 'EMPTY', r.get('model'), r.get('provider'), time.time() - t,
            'parses correctly' if ok else 'DID NOT PARSE: %r' % txt[:80]))
    except Exception as exc:
        print('   LIVE CALL   FAIL %s' % str(exc)[:160])
    print()

# ---- cost ------------------------------------------------------------------
rows = [json.loads(l) for l in open(_os.path.join(_PILOT, 'scores', 'judge_probe', 'replies.jsonl'), encoding='utf-8')]
e0 = [json.loads(l) for f in glob.glob(_os.path.join(_PILOT, 'scores', 'e0', '*.jsonl')) for l in open(f, encoding='utf-8')]
latest = max((r for r in e0 if not r.get('error')), key=lambda r: r['ts'])['config_sha256']
judged = sum(1 for r in e0 if r['config_sha256'] == latest and not r.get('error') and r['calls'])
e0_real = [c for r in e0 if r['config_sha256'] == latest and not r.get('error') for c in r['calls'] if c.get('ok')]
e0_in = st.mean(c['prompt_tokens'] for c in e0_real if c.get('prompt_tokens'))

print('COST OF E1 ON THE PILOT')
print('traces E0 sends to judges on the shared seed: %d (E1 uses the same trigger)' % judged)
print('real E0 Tribunal prompts average %.0f input tokens\n' % e0_in)
print('%-15s %9s %9s %11s %10s %11s' % ('judge', 'in tok', 'out tok', '$ / call', 'E1 total', 'if 2x out'))
tot = tot_hi = 0.0
for key, _ in PANEL:
    mine = [r for r in rows if r['judge'] == key and r.get('ok')]
    tin = max(st.mean(r['in'] for r in mine), e0_in)          # never below real prompt size
    tout = st.mean(r['out'] for r in mine)
    pin, pout = prices[key]
    per = (tin * pin + tout * pout) / 1e6
    hi = (tin * pin + 2 * tout * pout) / 1e6
    tot += per * judged
    tot_hi += hi * judged
    print('%-15s %9.0f %9.0f %11.4f %10.2f %11.2f' % (key, tin, tout, per, per * judged, hi * judged))
print('%-15s %9s %9s %11s %10.2f %11.2f' % ('PANEL', '', '', '', tot, tot_hi))
print('\n"if 2x out": a real E0 prompt asks about several steps at once where the probe asked')
print('about one, so the reply can be longer. Doubling the output is the planning ceiling.')
