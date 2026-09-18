"""Are the pilot's models listed, un-deprecated, served, and answering - today?"""
import os as _os
_ANALYSIS = _os.path.dirname(_os.path.abspath(__file__))
_PILOT_DIR = _os.path.dirname(_ANALYSIS)
_REPO = _os.path.dirname(_PILOT_DIR)
import datetime as dt
import json
import os
import time
import urllib.request

from dotenv import load_dotenv

load_dotenv(_os.path.join(_REPO, '.env'))
OR = os.environ['OPENROUTER_API_KEY']
GK = os.environ['GEMINI_API_KEY']


def get(url, headers=None):
    req = urllib.request.Request(url, headers={'User-Agent': 'engtrace', **(headers or {})})
    with urllib.request.urlopen(req, timeout=60) as fh:
        return json.load(fh)


H = {'Authorization': 'Bearer ' + OR}
cat = {m['id']: m for m in get('https://openrouter.ai/api/v1/models', H)['data']}


def ts(x):
    return dt.datetime.fromtimestamp(x, dt.timezone.utc).strftime('%Y-%m-%d') if x else None


def openrouter(mid):
    m = cat.get(mid)
    print('\n=== %s (OpenRouter)' % mid)
    if not m:
        print('  NOT LISTED')
        return
    print('  listed since   %s' % ts(m.get('created')))
    print('  expiration     %s' % (m.get('expiration_date') or 'none announced'))
    print('  price in/out   $%.2f / $%.2f per M' % (float(m['pricing']['prompt']) * 1e6,
                                                   float(m['pricing']['completion']) * 1e6))
    eps = get('https://openrouter.ai/api/v1/models/%s/endpoints' % mid, H)['data']['endpoints']
    for e in eps:
        up = e.get('uptime_last_30m')
        print('  endpoint %-16s status=%-4s uptime30m=%-6s max_out=%s quant=%s'
              % (e.get('provider_name'), e.get('status'),
                 '%.1f%%' % up if isinstance(up, (int, float)) else up,
                 e.get('max_completion_tokens'), e.get('quantization')))
    fam = mid.split('/')[0] + '/'
    base = mid.split('/')[1].split('-')[0]
    newer = sorted((k, ts(v.get('created'))) for k, v in cat.items()
                   if k.startswith(fam) and base in k and ':' not in k
                   and (v.get('created') or 0) > (m.get('created') or 0))
    print('  newer same-line ids listed: %s' % (', '.join('%s (%s)' % kv for kv in newer[:8]) or 'none'))
    probe(mid)


def probe(mid):
    body = json.dumps({'model': mid, 'max_tokens': 2048,
                       'messages': [{'role': 'user', 'content': 'Reply with the single word OK.'}]}).encode()
    req = urllib.request.Request('https://openrouter.ai/api/v1/chat/completions', data=body,
                                 headers={**H, 'Content-Type': 'application/json', 'User-Agent': 'engtrace'})
    t = time.time()
    try:
        with urllib.request.urlopen(req, timeout=300) as fh:
            r = json.load(fh)
        print('  LIVE CALL      OK  served=%s via %s in %.1fs, reply %r'
              % (r.get('model'), r.get('provider'), time.time() - t,
                 (r['choices'][0]['message'].get('content') or '').strip()[:20]))
    except Exception as exc:
        print('  LIVE CALL      FAIL %s' % str(exc)[:200])


def google(mid):
    print('\n=== %s (Google, direct)' % mid)
    try:
        m = get('https://generativelanguage.googleapis.com/v1beta/models/%s?key=%s' % (mid, GK))
        print('  listed         %s  version=%s' % (m.get('displayName'), m.get('version')))
        print('  description    %s' % (m.get('description') or '')[:160])
    except Exception as exc:
        print('  NOT LISTED     %s' % str(exc)[:160])
    allm = get('https://generativelanguage.googleapis.com/v1beta/models?pageSize=1000&key=%s' % GK)['models']
    pro = sorted(x['name'].split('/')[-1] for x in allm if 'gemini-3' in x['name'] and 'pro' in x['name'])
    print('  gemini-3* pro ids Google lists now: %s' % ', '.join(pro))
    body = json.dumps({'contents': [{'parts': [{'text': 'Reply with the single word OK.'}]}]}).encode()
    req = urllib.request.Request(
        'https://generativelanguage.googleapis.com/v1beta/models/%s:generateContent?key=%s' % (mid, GK),
        data=body, headers={'Content-Type': 'application/json'})
    t = time.time()
    try:
        with urllib.request.urlopen(req, timeout=300) as fh:
            r = json.load(fh)
        txt = r['candidates'][0]['content']['parts'][0].get('text', '')
        print('  LIVE CALL      OK  modelVersion=%s in %.1fs, reply %r'
              % (r.get('modelVersion'), time.time() - t, txt.strip()[:20]))
    except Exception as exc:
        print('  LIVE CALL      FAIL %s' % str(exc)[:200])


print('checked %s UTC' % dt.datetime.now(dt.timezone.utc).strftime('%Y-%m-%d %H:%M'))
print('\n##### TRACE SOURCES (stage 1)')
for mid in ('openai/gpt-5', 'anthropic/claude-opus-4.7', 'deepseek/deepseek-r1-0528',
            'meta-llama/llama-3.1-70b-instruct'):
    openrouter(mid)
google('gemini-3.1-pro-preview')
print('\n##### E0 JUDGES (stage 2) - gpt-5 and gemini-3.1-pro-preview are above')
openrouter('anthropic/claude-opus-4.5')
google('gemini-3-pro-preview')
